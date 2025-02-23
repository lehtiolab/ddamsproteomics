include { createMods; listify; stripchars_infile; parse_isotype } from '../modules.nf'

process msgfPlus {

  tag 'msgfplus'
  container params.__containers[tag][workflow.containerEngine]

  input:
  tuple val(setname), val(sample), path(mzml), val(maxmiscleav), val(fparams), val(mods), val(setisobaric), val(fractionation), val(minpeplen), val(maxpeplen), val(mincharge), val(maxcharge), path(db), path('mods.txt')

  output:
  tuple val(setname), path("${sample}.mzid"), path("${sample}.mzid.tsv")
  

  script:
  threads = task.cpus * params.overbook_cpus_factor
  // protocol 0 is automatic, msgf checks in mod file, TMT/phospho should be run with 1
  // see at https://github.com/MSGFPlus/msgfplus/issues/19
  msgfprotocol = fparams.phospho ? setisobaric[setname][0..4] == 'itraq' ? 3 : 1 : 0
  fragmeth = [auto:0, cid:1, etd:2, hcd:3, uvpd:4][fparams.frag]
  msgfinstrument = [lowres:0, velos:1, qe:3, qehf: 3, false:0, qehfx:1, lumos:1, timstof:2][fparams.instrument]
  enzyme = fparams.enzyme.indexOf('-') > -1 ? fparams.enzyme.replaceAll('-', '') : fparams.enzyme
  enzyme = [unspecific:0, trypsin:1, chymotrypsin: 2, lysc: 3, lysn: 4, gluc: 5, argc: 6, aspn:7, no_enzyme:9][enzyme]
  ntt = [full: 2, semi: 1, non: 0][fparams.terminicleaved]

  (is_stripped, parsed_infile) = stripchars_infile(mzml)
  """
  ${is_stripped ? "ln -s ${mzml} '${parsed_infile}'" : ''}
  
  msgf_plus -Xmx${task.memory.toMega()}M -d $db -s '$parsed_infile' -o "${sample}.mzid" -thread ${threads} -mod "mods.txt" -tda 0 -maxMissedCleavages ${maxmiscleav} -t ${fparams.prectol}  -ti ${fparams.iso_err} -m ${fragmeth} -inst ${msgfinstrument} -e ${enzyme} -protocol ${msgfprotocol} -ntt ${ntt} -minLength ${minpeplen} -maxLength ${maxpeplen} -minCharge ${mincharge} -maxCharge ${maxcharge} -n 1 -addFeatures 1
  msgf_plus -Xmx3500M edu.ucsd.msjava.ui.MzIDToTsv -i "${sample}.mzid" -o out.tsv
  rm ${db.baseName.replaceFirst(/\.fasta/, "")}.c*
  ${mods.contains('Unknown') ? "sed -i '/unknown modification/s/PSI-MS/UNIMOD/' '${sample}.mzid'" : ''}
  awk -F \$'\\t' '{OFS=FS ; print \$0, "Biological set" ${fractionation ? ', "Strip", "Fraction"' : ''}}' <( head -n1 out.tsv) > "${sample}.mzid.tsv"
  awk -F \$'\\t' '{OFS=FS ; print \$0, "$setname" ${fractionation ? ", \"$fparams.plate\", \"$fparams.fraction\"" : ''}}' <( tail -n+2 out.tsv) >> "${sample}.mzid.tsv"
  """
}


process percolator {

  tag 'percolator'
  container params.__containers[tag][workflow.containerEngine]

  input:
  tuple val(setname), path(mzids), val(enzyme)

  output:
  tuple val(setname), path('perco.xml')

  script:
  """
  ${listify(mzids).collect() { "echo $it >> metafile" }.join(' && ')}
  msgf2pin -o percoin.tsv -e ${enzyme} -P "decoy_" metafile
  percolator -j percoin.tsv -X perco.xml -N 500000 --decoy-xml-output -Y --num-threads ${task.cpus}
  """
}


process percolatorToPsms {

  tag 'msstitch'
  container params.__containers[tag][workflow.containerEngine]

  input:
  tuple val(setname), path(perco), path(mzids), path(tsvs), val(psmconf), val(pepconf), val(output_unfiltered)

  output:
  tuple val('target'), path("${setname}_target.tsv"), emit: tmzidtsv_perco optional true
  tuple val('decoy'), path("${setname}_decoy.tsv"), emit: dmzidtsv_perco optional true
  tuple val(setname), path('allpsms'), emit: unfiltered_psms optional true
  path('warnings'), emit: percowarnings optional true

  script:
  """
  mkdir outtables
  msstitch perco2psm --perco perco.xml -d outtables \
    -i ${listify(tsvs).collect(){"${it}"}.join(' ')} \
    --mzids ${listify(mzids).collect(){"${it}"}.join(' ')} \
    ${!output_unfiltered ? "--filtpsm ${psmconf} --filtpep ${pepconf}" : ''}
  msstitch concat -i outtables/* -o allpsms
  ${output_unfiltered ?
    "msstitch conffilt -i allpsms -o filtpsm --confcolpattern 'PSM q-value' --confidence-lvl ${psmconf} --confidence-better lower && \
    msstitch conffilt -i filtpsm -o psms --confcolpattern 'peptide q-value' --confidence-lvl ${pepconf} --confidence-better lower" : 'cp allpsms psms'}
  msstitch split -i psms --splitcol TD
  ${['target', 'decoy'].collect() { 
    "test -f '${it}.tsv' && mv '${it}.tsv' '${setname}_${it}.tsv' || echo 'No ${it} PSMs found for set ${setname} at PSM FDR ${psmconf} and peptide FDR ${pepconf} ${it == 'decoy' ? '(not possible to output protein/gene FDR)' : ''}' >> warnings" }.join(' && ') }
  """
}


workflow MSGFPERCO {
  take:
  mzml_in
  concatdb
  setisobaric
  fractionation
  fparam_map
  maxvarmods
  msgfmodfile
  mods
  psmconflvl
  pepconflvl
  output_unfilt_psms
  maxmiscleav
  enzyme
  minpeplen
  maxpeplen
  mincharge
  maxcharge

  main:
  mzml_in
  | map { [it.setname]}
  | groupTuple
  | map { [it[0], setisobaric && setisobaric[it[0]] ? setisobaric[it[0]] : false, maxvarmods, 'msgf', file(msgfmodfile), mods] }
  | createMods

  mzml_in
  | map { [it.setname, it.sample, it.mzmlfile, maxmiscleav, fparam_map[it.id], mods, setisobaric, fractionation, minpeplen, maxpeplen, mincharge, maxcharge] }
  | combine(concatdb)
  | combine(createMods.out, by: 0)
  | msgfPlus

  //  Nextflow-io #3970 - cannot pipe and tap {}, otherwise we can pull these three blocks together
  msgfPlus.out
  .groupTuple()
  .set { mzid_tsv_2psm }

  mzid_tsv_2psm
  | map {[it[0], it[1], enzyme] }
  | percolator
  | join(mzid_tsv_2psm)
  | map { it + [psmconflvl, pepconflvl, output_unfilt_psms] }
  | percolatorToPsms

  emit:
  t_tsv = percolatorToPsms.out.tmzidtsv_perco
  d_tsv = percolatorToPsms.out.dmzidtsv_perco
  unfiltered = percolatorToPsms.out.unfiltered_psms
  warnings = percolatorToPsms.out.percowarnings
}
