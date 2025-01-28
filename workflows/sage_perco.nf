include { createMods; listify; stripchars_infile; parse_isotype } from '../modules.nf'


process sagePrepare {
 //FIXME other container
 // container params.test ? 'nfhelaqc_test' : \
  //  "ghcr.io/lehtiolab/nfhelaqc:${workflow.manifest.version}"
  container 'quay.io/biocontainers/jq:1.5--4'
 
  input:
  tuple val(setname), val(id), val(minlen), val(maxlen), val(mincharge), val(maxcharge), val(maxmiscleav), val(maxvarmods), val(fparams), path('sage.json'), path('mods.json')

  output:
  tuple val(id), path('config.json')
  
  script:
  prectol_num = fparams.prectol.replaceAll('[^0-9\\.]', '')
  prectol_unit = fparams.prectol.replaceAll('[0-9\\.]', '')
  sage_maxmiscleav = maxmiscleav == -1 ? 20 : maxmiscleav
  """
  cat sage.json | jq --argjson MODS "\$(cat mods.json)" '.database.enzyme.missed_cleavages=${sage_maxmiscleav} | \
    .database+=\$MODS | \
    .database.enzyme.min_len=${minlen} | \
    .database.enzyme.max_len=${maxlen} | .database.max_variable_mods=${maxvarmods} | \
    .precursor_tol={$prectol_unit: [-${prectol_num}, ${prectol_num}]} | \
    .precursor_charge=[${mincharge}, ${maxcharge}]' > config.json
  """
    //.fragment_tol.ppm=[-${fparams.fragtol}, ${fparams.fragtol}] | \
}


process sage {
  container 'ghcr.io/lazear/sage:v0.14.7'

  input:
  tuple val(setname), path('config.json'), path(specfile), val(fparams), val(instrumenttype), val(fractionation), path(db)

  output:
  tuple val(setname), path("${specfile.baseName}.pin"), emit: perco
  tuple val(setname), path("${specfile.baseName}.tsv"), emit: tsv

  script:
  remove_scan_index_str = instrumenttype == 'bruker'
  """
  export RAYON_NUM_THREADS=${task.cpus}
  export SAGE_LOG=trace
  sage --disable-telemetry-i-dont-want-to-improve-sage --write-pin -f $db config.json $specfile
  ${remove_scan_index_str ? "sed -i 's/index=//' results.sage.tsv" : ''}
  ${remove_scan_index_str ? "sed -i 's/index=//' results.sage.pin" : ''}
  mv results.sage.pin "${specfile.baseName}.pin"
  awk -F \$'\\t' '{OFS=FS ; print \$0, "Biological set" ${fractionation ? ', "Strip", "Fraction"' : ''}}' <( head -n1 results.sage.tsv) > "${specfile.baseName}.tsv"
  awk -F \$'\\t' '{OFS=FS ; print \$0, "$setname" ${fractionation ? ", \"$fparams.plate\", \"$fparams.fraction\"" : ''}}' <( tail -n+2 results.sage.tsv) >> "${specfile.baseName}.tsv"
  """
}


process percolator {
  tag 'percolator'
  container params.__containers[tag][workflow.containerEngine]

  input:
  tuple val(setname), path(pins)

  output:
  tuple val(setname), path('perco.xml')

  script:
  """
  head -n1 ${pins[0]} > percoin.tsv
  cat ${pins.collect { "<(tail -n+2 $it)" }.join(' ')} >> percoin.tsv
  percolator -j percoin.tsv -X perco.xml -N 500000 --decoy-xml-output -Y --num-threads ${task.cpus}
  """
}


process percolatorToPsms {

  tag 'msstitch'
  container params.__containers[tag][workflow.containerEngine]

  input:
  tuple val(setname), path(perco), path(tsvs), val(psmconf), val(pepconf), val(output_unfiltered)

  output:
  tuple val('target'), path("${setname}_target.tsv"), emit: tmzidtsv_perco optional true
  tuple val('decoy'), path("${setname}_decoy.tsv"), emit: dmzidtsv_perco optional true
  tuple val(setname), path('allpsms'), emit: unfiltered_psms optional true
  path('warnings'), emit: percowarnings optional true

  script:
    //--mzids ${listify(mzids).collect(){"${it}"}.join(' ')} \
  """
  mkdir outtables
  msstitch perco2psm --perco perco.xml -d outtables \
    -i ${listify(tsvs).collect(){"${it}"}.join(' ')} \
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


process createPSMTable {
  container 'quay.io/biocontainers/msstitch:3.16--pyhdfd78af_0'

  input:
  tuple path('perco'), path('psms'), val(instrumenttype), path(lookup), path(db), path(ddb), val(psmconf), val(pepconf)

  output:
  tuple path('tpsmtable'), path('dpsmtable'), path('tpeptides'), path('dpeptides'), emit: tables
  path('tpsmlookup'), emit: lookup

  script:
  add_scan_index_str = instrumenttype == 'bruker'
  """
  msstitch perco2psm --perco perco -i psms -o psms_perco --filtpsm ${psmconf} --filtpep ${pepconf}
  ${add_scan_index_str ? "mv psms_perco ppnoi && head -n1 ppnoi > psms_perco" : ''}
  ${add_scan_index_str ? "paste <(awk -F'\t' -v OFS='\t' '{print \$1,\$2,\$3,\$4,\$5,\"index=\"\$6}' ppnoi) <(cut -f7- ppnoi) | tail -n+2 >> psms_perco" : ''}

  msstitch split -i psms_perco --splitcol TD
  cp $lookup tpsmlookup
  cp $lookup dpsmlookup
  msstitch psmtable -i target.tsv --dbfile tpsmlookup -o tpsmtable --fasta "$db" --ms1quant --proteingroup --spectracol ${Utils.get_field_nr('target.tsv', 'filename')}
  msstitch psmtable -i decoy.tsv --dbfile dpsmlookup -o dpsmtable --fasta "$ddb" --proteingroup --spectracol ${Utils.get_field_nr('decoy.tsv', 'filename')}
  msstitch peptides -i tpsmtable -o tpeptides --scorecolpattern sage_discriminant --spectracol ${Utils.get_field_nr('tpsmtable', 'filename')} --ms1quantcolpattern area
  msstitch peptides -i dpsmtable -o dpeptides --scorecolpattern sage_discriminant --spectracol ${Utils.get_field_nr('dpsmtable', 'filename')}
  """
}


workflow SAGEPERCO {
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
  | map { [it[0], setisobaric && setisobaric[it[0]] ? setisobaric[it[0]] : false, maxvarmods, 'sage', file(msgfmodfile), mods] }
  | createMods

  mzml_in
  | map { [it.setname, it.id, minpeplen, maxpeplen, mincharge, maxcharge, maxmiscleav, maxvarmods, fparam_map[it.id], file("$baseDir/assets/sage.json")] }
  | combine(createMods.out, by: 0)
  | sagePrepare
  | join(mzml_in.map { [it.id, it.setname, it.mzmlfile, fparam_map[it.id], it.instrument == 'timstof' ? 'bruker' : 'thermo', fractionation] })
  | map { [it[2], it[1], it[3..-1]].flatten() } // Set, config, rest
  | combine(concatdb)
  | sage

  sage.out.perco
  | groupTuple
  | view()
  | percolator
  | join(sage.out.tsv)
  | map { [it, psmconflvl, pepconflvl, output_unfilt_psms].flatten() }
  | percolatorToPsms

//  tuple val(id), path('config.json'), path(specfile), val(instrumenttype), path(mods), path(db)

  emit:
  t_tsv = percolatorToPsms.out.tmzidtsv_perco
  d_tsv = percolatorToPsms.out.dmzidtsv_perco
  unfiltered = percolatorToPsms.out.unfiltered_psms
  warnings = percolatorToPsms.out.percowarnings
}
