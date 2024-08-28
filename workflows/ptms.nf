include { listify; stripchars_infile; get_field_nr; parse_isotype } from '../modules.nf' 


process luciphorPTMLocalizationScoring {

  input:
  tuple val(setname), val(locptms), path(mzmls), path('psms'), path(tdb), path(allpsms), val(all_non_ptm_mods), val(stab_ptms), val(activation), val(maxpeplen), val(maxcharge), path(msgfmodfile), val(minpsms_luciphor), val(ptm_minscore_high)

  output:
  tuple val(setname), path('labileptms.txt'), emit: ptms
  path 'warnings', emit: warnings optional true

  script:
  stab_ptms = stab_ptms ? listify(stab_ptms).join(' ') : ''
  lab_ptms = listify(locptms).join(' ')

  """
  ${mzmls.collect() { stripchars_infile(it, return_oldfile=true) }.findAll{ it[0] }.collect() { "ln -s '${it[2]}' '${it[1]}'" }.join(' && ')}
  # Split allpsms to get target PSMs
  sed '0,/\\#SpecFile/s//SpectraFile/' -i "${allpsms}"
  msstitch split -i "${allpsms}" --splitcol TD
  export MZML_PATH=\$(pwd)
  export MINPSMS=${minpsms_luciphor}
  export ALGO=${['hcd', 'auto', 'any'].any { it == activation } ? '1' : '0'}
  export THREAD=task.cpus
  export MAXPEPLEN=${maxpeplen}
  export MAXCHARGE=${maxcharge}
  export MS2TOLVALUE=0.025
  export MS2TOLTYPE=Da
  cat "$baseDir/assets/luciphor2_input_template.txt" | envsubst > lucinput.txt
  luciphor_prep.py --psmfile target.tsv --template lucinput.txt --modfile "${msgfmodfile}" \
      --labileptms "${lab_ptms}" --mods ${all_non_ptm_mods} ${stab_ptms} \
      -o luciphor.out --lucipsms lucipsms
  luciphor2 -Xmx${task.memory.toGiga()}G luciphor_config.txt 2>&1 | grep 'not have enough PSMs' && echo 'Not enough PSMs for luciphor FLR calculation in set ${setname}' > warnings
  luciphor_parse.py --minscore ${ptm_minscore_high} -o labileptms.txt \
     --luci_in luciphor.out --luci_scores all_scores.debug --psms psms \
     --modfile "${msgfmodfile}" --labileptms ${lab_ptms} \
     ${stab_ptms ? "--stabileptms ${stab_ptms}": ''} --mods ${all_non_ptm_mods} \
     --fasta "${tdb}"
  """
}
// FIXME msgfmods is false? oxidation so probably never.


process stabilePTMPrep {

  input:
  tuple val(setname), val(ptms), path('psms'), path(tdb), val(non_ptm_mods), val(lab_ptms), path(msgfmods)
  
  output:
  tuple val(setname), path('stabileptms.txt')
  
  script:
  stab_ptms = listify(ptms).join(' ')
  lab_ptms = listify(lab_ptms).join(' ')
  """
  nonlabile_ptm_columns.py --psms psms -o stabileptms.txt --modfile "${msgfmods}" --fasta "${tdb}" \
      --stabileptms $stab_ptms \
      ${lab_ptms ? "--labileptms $lab_ptms" : ""} \
      ${non_ptm_mods ? "--mods ${non_ptm_mods}" : ''}
  """
}

process createPTMLookup {
  container "${workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
    'https://depot.galaxyproject.org/singularity/msstitch:3.16--pyhdfd78af_0' :
    'quay.io/biocontainers/msstitch:3.16--pyhdfd78af_0'}"
  
//publishDir "${params.outdir}", mode: 'copy', overwrite: true, saveAs: {it == ptmtable ? ptmtable: null}

  input:
  tuple val(setnames), path('ptms?'), path('stabileptms*'), path(cleaned_oldptms), path('speclup.sql')

  output:
  tuple path(ptmtable), path('ptmlup.sql'), emit: allptmpsms_db
  path({setnames.collect() { "${it}.tsv" }}), emit: setptmpsm optional true
  path('warnings'), emit: warnings optional true

  script:
  ptmtable = "ptm_psmtable.txt"
  oldptms = cleaned_oldptms.name != 'NO__FILE' ? cleaned_oldptms : false
  """
  # Concat all the PTM PSM tables (labile, stabile, previous) and load into DB
  # PSM peptide sequences include the PTM site
  cat speclup.sql > ptmlup.sql
  msstitch concat -i stabileptms* ptms* ${oldptms ?: ''} -o concatptmpsms
  msstitch psmtable -i concatptmpsms --dbfile ptmlup.sql -o ${ptmtable} \
    --spectracol 1
  msstitch split -i ${ptmtable} --splitcol bioset
  ${setnames.collect() { "test -f '${it}.tsv' || echo 'No PTMs found for set ${it}' >> warnings" }.join(' && ') }
  # No PTMs at all overwrites the per-set messages
  if [[ \$(wc -l <ptm_psmtable.txt) -lt 2 ]]
  then
    echo 'No PTMs found in any set' > warnings
  fi
  """
}

process PTMPeptides {
  container "${workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
    'https://depot.galaxyproject.org/singularity/msstitch:3.16--pyhdfd78af_0' :
    'quay.io/biocontainers/msstitch:3.16--pyhdfd78af_0'}"

  input:
  tuple val(setname), path('ptms.txt'), path(tppsms), val(isobtype), val(denom), val(dividebycol), val(normalize), val(keepnapsms_quant)

  output:
  tuple val(setname), path(peptable), path("${peptable}_no_tp_normalized")

  script:
  specialdenom = denom && (denom[0] == 'sweep' || denom[0] == 'intensity')
  peptable = "${setname}_ptm_peptides.txt"
  tppsms = tppsms.name != 'NO__FILE' ? tppsms : false
  """
  # If there is a total proteome PSMs file, prepare a gene (or protein) and protein table
  # for peptide normalization purposes. Use msstitch isosummarize here so we dont have to deal 
  # with gene FDR and peptide tables etc. First the always non-normalized table to get
  # proteins with which to median-center the PTM table with (if --normalize is passed)
  ${tppsms && normalize ? "msstitch isosummarize -i $tppsms -o prots_mediancenter \
    --featcol ${get_field_nr(tppsms, dividebycol)} \
    --isobquantcolpattern ${isobtype} --minint 0.1 \
    ${denom && denom[0] == 'sweep' ? '--mediansweep --logisoquant': ''} \
    ${denom && keepnapsms_quant ? '--keep-psms-na-quant' : ''} \
    ${denom && !specialdenom ? "--logisoquant --denompatterns ${denom.join(' ')}": ''}" : ''}
  # PSM ID table is without master protein in onlypeptides
  ${tppsms && normalize && dividebycol != 'Gene Name' ? "sed -i '1s/${dividebycol}/Protein ID/' prots_mediancenter" : ''}
  
  # Then as above but output the median-centered (if applicable) gene (default) or protein 
  # (in case of --onlypeptides) table for relating the peptides to their gene/protein
  ${tppsms && denom ? "msstitch isosummarize -i $tppsms -o tp_accessions \
    --featcol ${get_field_nr(tppsms, dividebycol)} \
    --isobquantcolpattern ${isobtype} --minint 0.1 \
    ${normalize ? "--median-normalize" : ''} \
    ${denom && denom[0] == 'sweep' ? '--mediansweep --logisoquant': ''} \
    ${denom && keepnapsms_quant ? '--keep-psms-na-quant' : ''} \
    ${denom && !specialdenom ? "--logisoquant --denompatterns ${denom.join(' ')}": ''}" : ''}
  ${tppsms && denom && dividebycol != 'Gene Name' ? "sed -i '1s/${dividebycol}/Protein ID/' tp_accessions" : ''}

  # Create a PTM-peptide table which has normalized isobaric data
  msstitch peptides -i "ptms.txt" -o "${peptable}" --scorecolpattern svm --spectracol 1 \
    ${!params.noquant && !params.noms1quant ? '--ms1quantcolpattern area' : ''} \
    ${denom ? "--isobquantcolpattern ${isobtype} --minint 0.1 --keep-psms-na-quant" : ''} \
    ${denom && denom[0] == 'sweep' ? '--mediansweep --logisoquant': ''} \
    ${denom && denom[0] == 'intensity' ? '--medianintensity' : ''} \
    ${denom && !specialdenom ? "--logisoquant --denompatterns ${denom.join(' ')}": ''} \
    ${denom && tppsms ? "--totalproteome tp_accessions" : ''} \
    ${denom && tppsms && normalize ? "--median-normalize --normalization-factors-table prots_mediancenter" : ''} \
  
  # And another one with the non-normalized isobaric quant (but with median-centering if specified)
  # this is not strictly necessary to do if there is no --totalproteomepsms, but will not be output 
  # in next step (merge) anyway then.
  msstitch peptides -i "ptms.txt" -o "${peptable}_no_tp_normalized" --scorecolpattern svm --spectracol 1 \
    ${denom ? "--isobquantcolpattern ${isobtype} --minint 0.1 --keep-psms-na-quant" : ''} \
    ${denom && denom[0] == 'sweep' ? '--mediansweep --logisoquant': ''} \
    ${denom && denom[0] == 'intensity' ? '--medianintensity' : ''} \
    ${denom && !specialdenom ? "--logisoquant --denompatterns ${denom.join(' ')}": ''} \
    ${denom && tppsms && normalize ? "--median-normalize --normalization-factors-table prots_mediancenter" : ''} \
  """
}

process mergePTMPeps {
  container "${workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
    'https://depot.galaxyproject.org/singularity/msstitch:3.16--pyhdfd78af_0' :
    'quay.io/biocontainers/msstitch:3.16--pyhdfd78af_0'}"
  // publishDir "${params.outdir}", mode: 'copy', overwrite: true, saveAs: {it.startsWith('ptm_peptides_') ? it : null}
 
  input:
  tuple val(setnames), path(peptides), path(notp_adjust_peps), path('ptmpsms.txt'), path('ptmlup.sql'), val(do_ms1), val(do_isobaric)

  output:
//  path peptable
//  path peptable_no_adjust optional true
//  tuple path("ptmqc.html"), path('summary.txt'), path('featcount_summary.txt') into ptmqc
//  path('overlap.txt') into ptmoverlap optional true

  script:
  peptable = params.totalproteomepsms ? 'ptm_peptides_total_proteome_adjusted.txt' : 'ptm_peptides_not_adjusted.txt'
  peptable_no_adjust = 'ptm_peptides_not_adjusted.txt'
  """
  cat ptmlup.sql > pepptmlup.sql
  # Create first table, input for which is either adjusted or not
  msstitch merge -i ${listify(peptides).collect() { "$it" }.join(' ')} --setnames ${setnames.collect() { "'$it'" }.join(' ')} --dbfile pepptmlup.sql -o mergedtable --no-group-annotation \
    --fdrcolpattern '^q-value' --pepcolpattern 'peptide PEP' --flrcolpattern 'FLR' \
    ${do_ms1 ? "--ms1quantcolpattern area" : ''} \
    ${do_isobaric ? "--isobquantcolpattern plex" : ''}
  # Add master/genes/gene count to peptide table, cant store in SQL because cant protein group on small PTM table
  ${!params.onlypeptides ? "sed \$'1s/Peptide sequence/Peptide sequence\tMaster protein(s)\tGene name(s)\tPTM flanking sequence(s)\tMatching gene count/' < mergedtable > '${peptable}'" : ''}
  ${!params.onlypeptides ? """geneprotcols=\$(head -1 ptmpsms.txt| tr '\\t' '\\n' | grep -En '(^Peptide|^Master|^Gene Name|^PTM flanking)' | cut -f 1 -d':' | tr '\\n' ',' | sed 's/\\,\$//')
    tail -n+2 ptmpsms.txt | cut -f\$geneprotcols | sort -uk1b,1 > geneprots
    join -j1 -o auto -t '\t' <(paste geneprots <(cut -f3 geneprots | tr -dc ';\\n'| awk '{print length+1}')) <(tail -n+2 mergedtable | sort -k1b,1) >> ${peptable}""" : "mv mergedtable ${peptable}"}

  # If total-proteome quant adjustment input was used above, create a second merged NON-adjusted peptide tables
  ${params.totalproteomepsms ?  "msstitch merge -i ${listify(notp_adjust_peps).collect() { "$it" }.join(' ')} --setnames ${setnames.collect() { "'$it'" }.join(' ')} \
    --dbfile pepptmlup.sql -o mergedtable --no-group-annotation \
    --fdrcolpattern '^q-value' --pepcolpattern 'peptide PEP' --flrcolpattern 'FLR' \
    ${!params.noquant && !params.noms1quant ? "--ms1quantcolpattern area" : ''} \
    ${!params.noquant && setisobaric ? "--isobquantcolpattern plex" : ''}" : ''}
  # Add master/genes/gene count to peptide table, cant store in SQL because cant protein group on small PTM table
  ${!params.onlypeptides && params.totalproteomepsms ? """head -n1 mergedtable | sed \$'s/Peptide sequence/Peptide sequence\tMaster protein(s)\tGene name(s)\tPTM flanking sequence(s)\tMatching gene count/'> '${peptable_no_adjust}'
    geneprotcols=\$(head -1 ptmpsms.txt| tr '\\t' '\\n' | grep -En '(^Peptide|^Master|^Gene Name|^PTM flanking)' | cut -f 1 -d':' | tr '\\n' ',' | sed 's/\\,\$//')
    tail -n+2 ptmpsms.txt | cut -f\$geneprotcols | sort -uk1b,1 > geneprots
    join -j1 -o auto -t '\t' <(paste geneprots <(cut -f3 geneprots | tr -dc ';\\n'| awk '{print length+1}')) <(tail -n+2 mergedtable | sort -k1b,1) >> ${peptable_no_adjust}""" : "${params.totalproteomepsms && params.onlypeptides ? "mv mergedtable ${peptable_no_adjust}" : ''}"}

  qc_ptms.R ptmpsms.txt "${peptable}"

  echo "<html><body>" > ptmqc.html
  for graph in ptmpsmfeats ptmpepfeats ptmprotfeats psmptmresidues pepptmresidues;
    do
    [ -e \$graph ] && echo "<div class=\\"chunk\\" id=\\"\${graph}\\"> \$(sed "s/id=\\"/id=\\"ptm-\${graph}/g;s/\\#/\\#ptm-\${graph}/g" <\$graph) </div>" >> ptmqc.html
    done 
  echo "</body></html>" >> ptmqc.html
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
setisobaric
mzmls
psms
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


main:

Channel.from(locptms)
| combine(mzmls)
| map { [it[1], it[0], it[2]] }
| join(psms)
| combine(tdb)
| join(unfiltered_psms)
| map { it + [ get_non_ptms(it[0], setisobaric, othermods),
stab_ptms, fparams.find { x -> x.setname == it[0] }.activation, maxpeplen, maxcharge, file(msgfmodfile), minpsms_luci, ptm_minscore_high] }
| luciphorPTMLocalizationScoring
  
Channel.from(stab_ptms)
| combine(psms)
| map { [it[1], it[0], it[2]] }
| combine(tdb)
| map { it + [get_non_ptms(it[0], setisobaric, othermods), locptms, file(msgfmodfile)] }
| stabilePTMPrep


luciphorPTMLocalizationScoring.out.ptms
| join(stabilePTMPrep.out)
| toList
| map { it.sort( {a, b -> a[0] <=> b[0]}) } // sort on setname
| transpose
| toList
| combine(ptm_psms_lookup)
| createPTMLookup

totalproteomepsms.view()
createPTMLookup.out.setptmpsm
| map{ it -> [listify(it).collect() { it.baseName.replaceFirst(/\.tsv$/, "") }, it]} // get setname from {setname}.tsv
| transpose
| join(totalproteomepsms, remainder: true)
| map { [it[0], it[1], it[2] ?: file(nofile), isobtypes[it[0]],
  setdenoms ? setdenoms[it[0]] : false,
  totalprot_col, normalize, keepnapsms_quant] }
| PTMPeptides
|view()
| toList
|view()
| map { it.sort( {a, b -> a[0] <=> b[0]}) } // sort on setname
|view()
| transpose
|view()
| toList
|view()
| combine(createPTMLookup.out.allptmpsms_db)
| mergePTMPeps


/*
  output:
  //path ptmtable into ptmpsmqc
  //path 'ptmlup.sql' into ptmlup
  //path({setnames.collect() { "${it}.tsv" }}) optional true into setptmtables
  //path 'warnings', emit: warnings optional true
ptmpeps
  .toList()
  .map { it.sort( {a, b -> a[0] <=> b[0]}) } // sort on setname
  .transpose()
  .toList()
  // have to do this map {} operation since even when PTMpeptides is not even run,
  // the toList seems to force a [[]] value on the channel
  .map { it == [[]] ? [[false], ['NA'], ['NA']] : it }
  .combine(ptmlup)
  .combine(ptmpsmqc)
  .set { ptmpeps2merge }

  output:
  tuple val(setname), path(peptable), path("${peptable}_no_tp_normalized")
process mergePTMPeps {
  input:
  tuple val(setnames), file(peptides), file(notp_adjust_peps), file('ptmlup.sql'), file('ptmpsms.txt') from ptmpeps2merge


// FIXME currently not storing the PTMs in SQL or similar, would be good to keep this state
// between runs for complement/rerun etc ?

if (params.ptmpsms && !(params.locptms || params.ptms)) exit 1, "In a rerun with --ptmpsms you  must specify which PTMs you have used with --locptms or --ptms"




*/
}
