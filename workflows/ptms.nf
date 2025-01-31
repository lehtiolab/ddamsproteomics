include { listify; stripchars_infile; get_field_nr; get_field_nr_multi; parse_isotype } from '../modules.nf' 


process luciphorPrep {

  tag 'msstitch'
  container params.__containers[tag][workflow.containerEngine]

  input:
  tuple val(setname), path(allpsms), val(locptms), val(stab_ptms), val(all_non_ptm_mods), path(msgfmodfile), val(search_engine)

  output:
  tuple val(setname), path('luciphor_config.txt'), path('lucipsms')

  script:
  stab_ptms = stab_ptms ? listify(stab_ptms).join(' ') : ''
  lab_ptms = listify(locptms).join(' ')
  """
  # Split allpsms to get target PSMs
  msstitch split -i "${allpsms}" --splitcol TD
  luciphor_prep.py --psmfile target.tsv --template "$baseDir/assets/luciphor2_input_template.txt" \
      --labileptms "${lab_ptms}" --mods ${all_non_ptm_mods} ${stab_ptms} \
      --modfile "${msgfmodfile}" --lucipsms lucipsms --searchengine $search_engine
  """
}


process luciphorPTMLocalizationScoring {

  tag 'luciphor2'
  container params.__containers[tag][workflow.containerEngine]

  input:
  tuple val(setname),  path(template), path('lucipsms'), path(mzmls), val(activation), val(maxpeplen), val(maxcharge), val(minpsms_luciphor)

  output:
  tuple val(setname), path('luciphor.out'), path('all_scores.debug'), emit: ptms optional true
  tuple val(setname), path('warnings'), emit: warnings optional true

  script:

  """
  ${mzmls.collect() { stripchars_infile(it, return_oldfile=true) }.findAll{ it[0] }.collect() { "ln -s '${it[2]}' '${it[1]}'" }.join(' && ')}
  export MZML_PATH=\$(pwd)
  export MINPSMS=${minpsms_luciphor}
  export ALGO=${['hcd', 'auto', 'any'].any { it == activation } ? '1' : '0'}
  export THREAD=${task.cpus}
  export MAXPEPLEN=${maxpeplen}
  export MAXCHARGE=${maxcharge}
  export MS2TOLVALUE=0.025
  export MS2TOLTYPE=0 # 0=Da, 1=ppm
  export OUTFILE=luciphor.out
  cat "${template}" | envsubst > lucinput.txt

  luciphor2 -Xmx${task.memory.toGiga()}G lucinput.txt 2> luci_stderr
  if grep 'not have enough PSMs' luci_stderr
  then
    echo 'Not enough PSMs for luciphor FLR calculation in set ${setname}' > warnings
  fi
  """
}

process luciphorParse {

  tag 'msstitch'
  container params.__containers[tag][workflow.containerEngine]

  // Puts luciphor data back into the PTM PSM table, also adds flanking seqs - if there 
  // is no luciphor data due to errors, it will put NA for luciphor columns

  input:
  tuple val(setname), path(luciphor_out), path('all_scores.debug'), path('psms'), val(ptm_minscore_high), path(msgfmodfile), val(lab_ptms), val(stab_ptms), val(all_non_ptm_mods), val(search_engine), path(tdb)

  output:
  tuple val(setname), path('labileptms.txt')

  script:
  luciphor_file = luciphor_out.name == 'NO__FILE' ? 'fail' : luciphor_out
  """
  luciphor_parse.py --minscore ${ptm_minscore_high} -o labileptms.txt \
     --luci_in ${luciphor_file} --luci_scores all_scores.debug --psms psms \
     --modfile "${msgfmodfile}" --labileptms ${lab_ptms.join(' ')} \
     ${stab_ptms ? "--stabileptms ${stab_ptms.join(' ')}": ''} --mods ${all_non_ptm_mods} \
     --fasta "${tdb}" --searchengine $search_engine
  """
}
// FIXME msgfmods is false? oxidation so probably never.


process stabilePTMPrep {

  tag 'msstitch'
  container params.__containers[tag][workflow.containerEngine]

  input:
  tuple val(setname), val(ptms), path('psms'), path(tdb), val(non_ptm_mods), val(lab_ptms), path(msgfmods), val(search_engine)
  
  output:
  tuple val(setname), path('stabileptms.txt')
  
  script:
  stab_ptms = listify(ptms).join(' ')
  lab_ptms = listify(lab_ptms).join(' ')
  """
  nonlabile_ptm_columns.py --psms psms -o stabileptms.txt --modfile "${msgfmods}" --fasta "${tdb}" \
      --searchengine $search_engine --stabileptms $stab_ptms \
      ${lab_ptms ? "--labileptms $lab_ptms" : ""} \
      ${non_ptm_mods ? "--mods ${non_ptm_mods}" : ''}
  """
}

process createPTMTable {
  
  tag 'msstitch'
  container params.__containers[tag][workflow.containerEngine]

  input:
  tuple val(setnames), path('ptms*'), val(has_newptms), path('speclup.sql'), path(cleaned_oldptms)

  output:
  path('ptmlup.sql'), emit: db 
  path(ptmtable), emit: allpsms
  path({setnames.collect() { "${it}.tsv" }}), emit: setptmpsm optional true
  path('warnings'), emit: warnings optional true

  script:
  ptmtable = "ptm_psmtable.txt"
  oldptms = cleaned_oldptms.name != 'NO__FILE' ? cleaned_oldptms : false
  """
  # Concat all the PTM PSM tables (labile, stabile, previous) and load into DB
  # PSM peptide sequences include the PTM site
  cat speclup.sql > ptmlup.sql
  ${has_newptms ? "msstitch concat -i ptms* ${oldptms ?: ''} -o" : "mv ${oldptms}"} concatptmpsms
  msstitch psmtable -i concatptmpsms --dbfile ptmlup.sql -o ${ptmtable}
  msstitch split -i ${ptmtable} --splitcol bioset
  ${setnames.collect() { "test -f '${it}.tsv' || echo 'No PTMs found for set ${it}' >> warnings" }.join(' && ') }
  # No PTMs at all overwrites the per-set messages
  if [[ \$(wc -l <ptm_psmtable.txt) -lt 2 ]]
  then
    echo 'No PTMs found in any set' > warnings
  fi
  """
}

process prepTotalProteomeInput {

  tag 'msstitch'
  container params.__containers[tag][workflow.containerEngine]

  input:
  tuple val(setname), path(tppsms), val(isobtype), val(denom), val(dividebycol), val(normalize), val(keepnapsms_quant)

  output:
  tuple val(setname), path('prots_mediancenter'), path('tp_accessions')

  script:
  specialdenom = denom && (denom[0] == 'sweep' || denom[0] == 'intensity')
  """
  # Prepare a gene (or protein) and protein table for peptide normalization purposes.
  # Use msstitch isosummarize here so we dont have to deal 
  # with gene FDR and peptide tables etc. First the always non-normalized table to get
  # proteins with which to median-center the PTM table with (if --normalize is passed)
  msstitch isosummarize -i $tppsms -o prots_mediancenter \
    --featcol ${get_field_nr(tppsms, dividebycol)} \
    --isobquantcolpattern ${isobtype} --minint 0.1 \
    ${denom && denom[0] == 'sweep' ? '--mediansweep --logisoquant': ''} \
    ${denom && keepnapsms_quant ? '--keep-psms-na-quant' : ''} \
    ${denom && !specialdenom ? "--logisoquant --denompatterns ${denom.join(' ')}": ''}
  # PSM ID table is without master protein in onlypeptides
  ${dividebycol != 'Gene Name' ? "sed -i '1s/${dividebycol}/Protein ID/' prots_mediancenter" : ''}
  
  # Then as above but output the median-centered (if applicable) gene (default) or protein 
  # (in case of --onlypeptides) table for relating the peptides to their gene/protein
  ${denom ? "msstitch isosummarize -i $tppsms -o tp_accessions \
    --featcol ${get_field_nr(tppsms, dividebycol)} \
    --isobquantcolpattern ${isobtype} --minint 0.1 \
    ${normalize ? "--median-normalize" : ''} \
    ${denom[0] == 'sweep' ? '--mediansweep --logisoquant': ''} \
    ${keepnapsms_quant ? '--keep-psms-na-quant' : ''} \
    ${!specialdenom ? "--logisoquant --denompatterns ${denom.join(' ')}": ''}" : ''}
  ${dividebycol != 'Gene Name' ? "sed -i '1s/${dividebycol}/Protein ID/' tp_accessions" : ''}
  """
}



process PTMPeptides {

  tag 'msstitch'
  container params.__containers[tag][workflow.containerEngine]

  input:
  tuple val(setname), path('ptms.txt'), path(tp_accessions), path('normalize_factors'),  val(isobtype), val(denom), val(normalize), val(keepnapsms_quant), val(do_ms1)

  output:
  tuple val(setname), val(tpnorm), path(peptable)

  script:
  specialdenom = denom && (denom[0] == 'sweep' || denom[0] == 'intensity')
  peptable = "${setname}_ptm_peptides.txt"
  tpnorm = tp_accessions.name != 'NO__FILE'
  """
  # Create a PTM-peptide table which has normalized isobaric data
  msstitch peptides -i "ptms.txt" -o "${peptable}" --scorecolpattern svm --spectracol 1 \
    ${do_ms1 ? '--ms1quantcolpattern area' : ''} \
    ${denom ? "--isobquantcolpattern ${isobtype} --minint 0.1 --keep-psms-na-quant" : ''} \
    ${denom && denom[0] == 'sweep' ? '--mediansweep --logisoquant': ''} \
    ${denom && denom[0] == 'intensity' ? '--medianintensity' : ''} \
    ${denom && !specialdenom ? "--logisoquant --denompatterns ${denom.join(' ')}": ''} \
    ${denom && tpnorm ? "--totalproteome ${tp_accessions}" : ''} \
    ${denom && tpnorm && normalize ? "--median-normalize --normalization-factors-table normalize_factors" : ''} \
  """
}


process mergePTMPeps {

  tag 'msstitch'
  container params.__containers[tag][workflow.containerEngine]
 
  input:
  tuple val(setnames), val(tpnormalized), path(peptides), path('ptmlup.sql'), val(do_ms1), val(do_isobaric)

  output:
  path(peptable)

  script:
  peptable = tpnormalized ? 'ptm_peptides_total_proteome_adjusted.txt' : 'ptm_peptides_not_adjusted.txt'
  """
  cat ptmlup.sql > pepptmlup.sql
  msstitch merge -i ${listify(peptides).collect() { "$it" }.join(' ')} --setnames ${setnames.collect() { "'$it'" }.join(' ')} --dbfile pepptmlup.sql -o ${peptable} --no-group-annotation \
    --fdrcolpattern '^q-value' --pepcolpattern 'peptide PEP' --flrcolpattern 'FLR' \
    ${do_ms1 ? "--ms1quantcolpattern area" : ''} \
    ${do_isobaric ? "--isobquantcolpattern plex" : ''}
  """
}


process addMasterProteinsGenes {

  // Runs no python but that container has the tools needed
  tag 'python'
  container params.__containers[tag][workflow.containerEngine]

  input:
  tuple path(peptides), path('ptmpsms')
  output:
  path(peptides)
  
  script:
  gene_fields_extra = ['Peptide sequence', 'Master protein(s)', 'Gene name(s)', 'PTM flanking sequence(s)', 'Matching gene count']
  psm_gene_fields = ['^Peptide', '^Master', '^Gene Name', '^PTM flanking']
  """
  mv $peptides no_genes # same output file name as input file!
  # Add master/genes/gene count to peptide table, cant store in SQL because cant protein group on small PTM table
  head -n1 no_genes | sed \$'1s/Peptide sequence/${gene_fields_extra.join('\\t')}/' > '${peptides}'
  tail -n+2 ptmpsms | cut -f${get_field_nr_multi('ptmpsms', psm_gene_fields)} | sort -uk1b,1 > geneprots
  join -j1 -o auto -t '\t' <(paste geneprots <(cut -f3 geneprots | tr -dc ';\\n'| awk '{print length+1}')) <(tail -n+2 no_genes | sort -k1b,1) >> ${peptides}
  """
}


def get_non_ptms(setname, isobtypes, mods) {
  def nonptms = []
  if (isobtypes && isobtypes[setname]) {
    nonptms.push(parse_isotype(isobtypes[setname]))
  }
  if (mods) {
    nonptms.push(mods.join(' '))
  }
  return nonptms.join(' ').trim()
}



workflow PTMANALYSIS {
  take:
  locptms
  stab_ptms
  othermods
  all_setnames
  setisobaric
  mzmls
  psms
  search_engine
  tdb
  unfiltered_psms
  fparams
  maxpeplen
  maxcharge
  msgfmodfile
  minpsms_luci
  ptm_minscore_high
  ptm_psms_lookup
  totalproteomepsms
  totalprot_col
  isobtypes
  setdenoms
  normalize
  keepnapsms_quant
  do_ms1
  do_proteingroup

  main:
  
  nofile = "${baseDir}/assets/NO__FILE"
  
  unfiltered_psms
  | map { it + [locptms, stab_ptms, get_non_ptms(it[0], setisobaric, othermods), file(msgfmodfile), search_engine] }
  | luciphorPrep
  | join(mzmls)
  | map { it + [fparams.find { x -> x.setname == it[0] }.activation, maxpeplen, maxcharge, minpsms_luci] }
  | luciphorPTMLocalizationScoring
    
  luciphorPTMLocalizationScoring.out.ptms
  // Join on warnings in case there is a problem w luciphor - we should not output 
  // empty then since the PTMs still need parsing to get flanks etc. After join the
  // channel will look like either of these:
  //  ----> setname, null, warnings || setname, psms, scores, null || empty
  | join(luciphorPTMLocalizationScoring.out.warnings, remainder: true)
  | map { [it[0], it[1] ?: nofile, it[1] ? it[2] : nofile] }
  | join(psms)
  | map { it + [ptm_minscore_high, file(msgfmodfile), locptms, stab_ptms, get_non_ptms(it[0], setisobaric, othermods), search_engine]}
  | combine(tdb)
  | luciphorParse

  Channel.from(stab_ptms)
  | combine(psms)
  | map { [it[1], it[0], it[2]] }
  | combine(tdb)
  | map { it + [get_non_ptms(it[0], setisobaric, othermods), locptms, file(msgfmodfile), search_engine] }
  | stabilePTMPrep
  | concat(luciphorParse.out)
  | toList
  | map { it.sort( {a, b -> a[0] <=> b[0]}) } // sort on setname
  | transpose
  | toList
  | filter { it[0] } // Empty (reruns) will show up as [] due to toList, filter those out
  // Dupe removal (there are sets for both stable and labile PSMs), and add true for has_newptms
  | map { [it[0].unique(), it[1], true] }
  | ifEmpty([all_setnames, file(nofile), false])
  | filter { it[1] != [null] }
  | combine(ptm_psms_lookup)
  | createPTMTable
  
  totalproteomepsms
  | map { it + [setisobaric[it[0]], setdenoms[it[0]], totalprot_col, normalize, keepnapsms_quant] }
  | prepTotalProteomeInput
  
  createPTMTable.out.setptmpsm
  | map { it -> [listify(it).collect() { it.baseName.replaceFirst(/\.tsv$/, "") }, it]} // get setname from {setname}.tsv
  | transpose
  | set { setptmpsm_ch }
  
  setptmpsm_ch
  | map { [it[0], null, null] }
  | concat(prepTotalProteomeInput.out)
  | set { totalprot_and_nofile }
  
  setptmpsm_ch
  | combine(totalprot_and_nofile, by:0)
  | map { it[0..1] + (it[2] ? [it[2], it[3]] : [file(nofile), file(nofile)]) + [isobtypes[it[0]],
    setdenoms ? setdenoms[it[0]] : false,
    normalize, keepnapsms_quant, do_ms1] }
  | PTMPeptides
  | toList
  | map { it.sort( {a, b -> a[0] <=> b[0]}) } // sort on setname
  | transpose
  | toList
  | transpose
  | groupTuple(by: 1) // group by total proteome normalization boolean
  | combine(createPTMTable.out.db)
  | map { it + [do_ms1, isobtypes] }
  | mergePTMPeps
  
  if (do_proteingroup) {
    mergePTMPeps.out
    | combine(createPTMTable.out.allpsms)
    | addMasterProteinsGenes
  }

  emit:
  ptms = createPTMTable.out.allpsms
  | combine(do_proteingroup ? addMasterProteinsGenes.out : mergePTMPeps.out)
  warnings = luciphorPTMLocalizationScoring.out.warnings.map { it[1] }
  | concat(createPTMTable.out.warnings)

}
  
// FIXME currently not storing the PTMs in SQL or similar, would be good to keep this state
// between runs for complement/rerun etc ?
