#!/usr/bin/env nextflow

import groovy.json.JsonSlurper
include { paramsSummaryMap } from 'plugin/nf-validation'

include { msgf_info_map; listify; stripchars_infile; get_regex_specialchars; read_header } from './modules.nf' 
include { MSGFPERCO } from './workflows/msgf_perco.nf'
include { PTMANALYSIS } from './workflows/ptms.nf'
include { MATCH_SEQUENCES } from './workflows/match_sequences.nf'
include { REPORTING } from './workflows/reporting.nf'

/*
========================================================================================
                         lehtiolab/ddamsproteomics
========================================================================================
 lehtiolab/ddamsproteomics Analysis Pipeline.
 #### Homepage / Documentation
 https://github.com/lehtiolab/ddamsproteomics
----------------------------------------------------------------------------------------
*/

process createTargetDecoyFasta {
  label 'msstitch'
 
  input:
  path(tdb)

  output:
  path('db.fa'), emit: concatdb
  tuple path('tdb'), path("decoy.fa"), emit: bothdbs

  script:
  """
  ${listify(tdb).size() > 1 ? "cat ${tdb.collect() { "\"${it}\"" }.join(' ')} > tdb" : "mv '$tdb' tdb"}
  check_fasta.py tdb
  msstitch makedecoy -i tdb -o decoy.fa --scramble tryp_rev --ignore-target-hits
  cat tdb decoy.fa > db.fa
  """
}

def fr_or_file(it, length) {
  // returns either fraction number or file from line
  // file name is used in QC plots where no frac is available and fraction plot is enforced,
  // e.g. when mixing fractions and non-fractions
  return it.size() > length ? it[length] : "${file(it[0]).baseName}.${file(it[0]).extension}"
}


// Parse mzML input to get files and sample names etc
// get setname, sample name (baseName), input mzML file. 
// Set platename to setname if not specified. 
// Set fraction name to NA if not specified

/*
* Step 1: Extract quant data from peptide spectra
*/


process centroidMS1 {

  input:
  tuple val(setname), path(infile), val(instr)

  output:
  tuple val(setname), val(parsed_infile), path(parsed_infile), val(instr)

  script:
  (is_stripped, parsed_infile) = stripchars_infile(infile)
  """
  ${is_stripped ? "ln -s ${infile} ${parsed_infile}" : ''}
  wine msconvert ${parsed_infile} --outfile centroidms1.mzML --filter 'peakPicking true 1' ${instr == 'timstof' ? "--filter sortByScanTime" : ''}
  mv centroidms1.mzML ${parsed_infile}
  """
}


process isobaricQuant {

  input:
  tuple val(setname), val(parsed_infile), path(infile), val(instr), val(isobtype)

  output:
  tuple val(parsed_infile), path(outfile)

  script:
  outfile = "${infile.baseName}.consensusXML"
  activationtype = [auto: 'auto', any: 'any', hcd:'beam-type collision-induced dissociation', cid:'Collision-induced dissociation', etd:'Electron transfer dissociation'][params.activation]
  plextype = isobtype ? isobtype.replaceFirst(/[0-9]+plex/, "") : 'false'
  massshift = [tmt:0.0013, itraq:0.00125, false:0][plextype]
  """
  ${isobtype ? "IsobaricAnalyzer -type $isobtype -in \"${infile}\" -out \"${infile.baseName}.consensusXML\" -extraction:select_activation \"$activationtype\" -extraction:reporter_mass_shift $massshift -extraction:min_precursor_intensity 1.0 -extraction:keep_unannotated_precursor true -quantification:isotope_correction true" : ''}
  """
}



process dinosaur {

  input:
  tuple val(sample), path(infile)

  output:
  tuple val(parsed_infile), path("${sample}.features.tsv")

  script:
  threads = task.cpus * params.overbook_cpus_factor
  (is_stripped, parsed_infile) = stripchars_infile(infile)
  """
  ${is_stripped ? "ln -s ${infile} ${parsed_infile}" : ''}
  dinosaur -Xmx${task.memory.toMega()}M --concurrency=${threads} ${parsed_infile}
  """
}


process hardklor {

  input:
  tuple val(sample), path(infile), path(hkconf)

  output:
  tuple val(sample), val(parsed_infile), path('hardklor.out')

  script:
  (is_stripped, parsed_infile) = stripchars_infile(infile)
  """
  ${is_stripped ? "ln -s ${infile} ${parsed_infile}" : ''}
  hardklor <(cat $hkconf <(echo $parsed_infile hardklor.out))
  """
}



process kronik {

  input:
  tuple val(sample), val(parsed_infile), path('hardklor.out')
  
  output:
  tuple val(parsed_infile), path("${sample}.features.tsv")
  
  script:
  """
  kronik -c 5 -d 3 -g 1 -m 8000 -n 600 -p 10 hardklor.out ${sample}.features.tsv
  """
}


process PTMClean {
  label 'sqlite'
  // FIXME can be in PTM wf?
  // In PTMS we need to delete all PSMs since we will rebuild it

  input:
  path('db.sqlite')
  
  output:
  path('ptms_db.sqlite')

  script:
  """
  cp db.sqlite ptms_db.sqlite
  tables=(
    psms
    psmrows
    peptide_sequences
    fastafn
    proteins
    protein_evidence
    protein_seq
    prot_desc
    protein_psm
    genes
    associated_ids
    ensg_proteins
    genename_proteins
  )
  for t in "\${tables[@]}"; do
    sqlite3 ptms_db.sqlite "DROP TABLE IF EXISTS \${t}" 
  done
  """
}


process complementSpectraLookupCleanPSMs {
  label 'msstitch'

  input:
  tuple val(in_setnames), path(mzmlfiles), val(platenames), path(tlup), path(dlup), path(tpsms), path(dpsms), path(ptmpsms)

  output:
  tuple path('t_cleaned_psms.txt'), path('d_cleaned_psms.txt'), emit: psms
  tuple path('target_db.sqlite'), path('decoy_db.sqlite'), emit: dbs
  tuple path('cleaned_ptmpsms.txt'), emit: ptm optional true
  
  script:
  setnames = in_setnames.unique(false)
  ptms = ptmpsms.name != 'NO__FILE'
  """
  # If this is an addition to an old lookup, copy it and extract set names
  cp ${tlup} target_db.sqlite
  cp ${dlup} decoy_db.sqlite
  python3 -c "exec(\\"import sqlite3;c=sqlite3.Connection('${tlup}');x=c.execute('SELECT set_name FROM biosets').fetchall();print('\\\\\\n'.join([y[0] for y in x]))\\")" > old_setnames
  # If adding to old lookup: grep new setnames in old and run msstitch deletesets if they match
  # use -x for grep since old_setnames must grep whole word
  if grep -xf old_setnames <(echo ${setnames.join('\n')} )
    then
      msstitch deletesets -i ${tpsms} -o t_cleaned_psms.txt --dbfile target_db.sqlite --setnames ${setnames.collect() { "'${it}'" }.join(' ')}
      msstitch deletesets -i ${dpsms} -o d_cleaned_psms.txt --dbfile decoy_db.sqlite --setnames ${setnames.collect() { "'${it}'" }.join(' ')}
      ${ptms ? "msstitch deletesets -i ${ptmpsms} -o cleaned_ptmpsms.txt --setnames ${setnames.collect() {"'${it}'"}.join(' ')}" : ''}
    else
      # In case there is a new set not present in the old data, just move things
      mv ${tpsms} t_cleaned_psms.txt
      mv ${dpsms} d_cleaned_psms.txt
      ${ptms ? "mv ${ptmpsms} cleaned_ptmpsms.txt" : ''}
  fi
  
  ${mzmlfiles.collect() { stripchars_infile(it, return_oldfile=true) }.findAll{ it[0] }.collect() { "ln -s '${it[2]}' '${it[1]}'" }.join(' && ')}
  ${mzmlfiles.size() ? "msstitch storespectra --spectra ${mzmlfiles.collect() { "'${it.toString().replaceAll('[&<>\'"]', '_')}'" }.join(' ')} --setnames ${in_setnames.collect() { "'$it'" }.join(' ')} --dbfile target_db.sqlite" : ''}

  copy_spectra.py target_db.sqlite decoy_db.sqlite ${setnames.join(' ')}
  cat old_setnames <(echo ${setnames.join('\n')}) | sort -u | grep -v '^\$' > all_setnames
  """
}


process createNewSpectraLookup {

  label 'msstitch'

  input:
  tuple val(setnames), file(mzmlfiles), val(platenames)

  output:
  path('target_db.sqlite')

  script:
  """
  ${mzmlfiles.collect() { stripchars_infile(it, return_oldfile=true) }.findAll{ it[0] }.collect() { "ln -s '${it[2]}' '${it[1]}'" }.join(' && ')}

  msstitch storespectra --spectra ${mzmlfiles.collect() { stripchars_infile(it)[1] }.join(' ')} --setnames ${setnames.collect() { "'$it'" }.join(' ')} -o target_db.sqlite
  """
}


process quantLookup {
  
  label 'msstitch'

  input:
  tuple val(mzmlnames), path(isofns), path(ms1fns), path(tlookup)

  output:
  path('target.sqlite')

  script:
  """
  # SQLite lookup needs copying to not modify the input file which would mess up a rerun with -resume
  cat $tlookup > target.sqlite
  msstitch storequant --dbfile target.sqlite --spectra ${mzmlnames.collect() { "'${it}'" }.join(' ')}  \
    ${!params.noms1quant ? "--mztol ${params.ms1qmztol} --mztoltype ppm --rttol ${params.ms1qrttol} ${params.hardklor ? "--kronik ${ms1fns.collect() { "$it" }.join(' ')}" : "--dinosaur ${ms1fns.collect() { "$it" }.join(' ')}"}" : ''} \
    ${params.isobaric ? "--isobaric ${isofns.collect() { "$it" }.join(' ')}" : ''}
  """
}


process createPSMTable {

  label 'msstitch'

  input:
  tuple val(td), path(psms), path('lookup'), path(tdb), path(ddb), path('oldpsms'), val(complementary_run), val(do_ms1), val(do_isobaric), val(onlypeptides)

  output:
  tuple val(td), path("${outpsms}"), emit: psmtable
  tuple val(td), path("${psmlookup}"), emit: lookup
  path('warnings'), emit: warnings optional true

  script:
  psmlookup = "${td}_psmlookup.sql"
  outpsms = "${td}_psmtable.txt"
  no_target = td == 'target' && psms.size() == 0
  no_decoy = td == 'decoy' && psms.size() == 0
  """
  ${no_target ? "echo 'No target PSMs made the combined PSM / peptide FDR cutoff' && exit 1" : ''}
  ${no_decoy ? "echo 'No decoy PSMs in any set at FDR cutoff, will not produce protein/gene tables' > warnings && touch ${outpsms} && touch ${psmlookup} && exit 0" : ''}
  msstitch concat -i ${listify(psms).collect() {"$it"}.join(' ')} -o psms.txt
  # SQLite lookup needs copying to not modify the input file which would mess up a rerun with -resume
  cat lookup > $psmlookup
  msstitch psmtable -i psms.txt --dbfile $psmlookup -o ${outpsms} \
      ${td == 'target' ? '--addmiscleav' : ''} \
      ${onlypeptides ? '' : "--fasta \"${td == 'target' ? "${tdb}" : "${ddb}"}\" --genes"} \
      ${do_ms1 ? '--ms1quant' : ''} \
      ${do_isobaric ? "--isobaric --min-precursor-purity ${params.minprecursorpurity}" : ''} \
      ${!onlypeptides ? '--proteingroup' : ''} \
      ${complementary_run ? '--oldpsms oldpsms' : ''}
  #sed '1s/\\#SpecFile/SpectraFile/' -i $outpsms
  """
}


process peptidePiAnnotation {

  label 'python'

  input:
  tuple path('psms'), val(strips), path('hirief_training_pep')

  output:
  path('target_psmtable.txt')

  script:
    // ${do_hirief && td == 'target' ? 
  """
  echo '${groovy.json.JsonOutput.toJson(strips)}' >> strip.json
  peptide_pi_annotator.py -i hirief_training_pep -p psms --out target_psmtable.txt --stripcolpattern Strip --pepcolpattern Peptide --fraccolpattern Fraction --stripdef strip.json --ignoremods '*'
  """
}


process splitPSMs {

  label 'msstitch'

  input:
  tuple val(td), path('psms'), val(setnames)

  output:
  tuple val(td), path({listify(setnames).collect { "${it}.tsv" }}) optional true

  script:
  """
  msstitch split -i psms --splitcol bioset
  """
}


process splitTotalProteomePSMs {

  label 'msstitch'

  input:
  tuple path('tppsms_in'), val(setnames)

  output:
  path({listify(setnames).collect() { "tppsms/${it}.tsv" }}), emit: totalprotpsms_allsets optional true

  script:
  """
  mkdir tppsms && msstitch split -i tppsms_in -d tppsms --splitcol bioset
  """
}


process makePeptides {

  label 'msstitch'

  input:
  tuple val(td), val(setname), path('psms'), val(setisobaric), val(denoms), val(keepnapsms_quant), val(normalize_isob), val(do_ms1)
  
  output:
  tuple val(setname), val(td), path("${setname}_peptides"), emit: peps
  tuple val('peptides'), path(normfactors), emit: normfacs optional true

  script:
  specialdenom = denoms && (denoms[0] == 'sweep' || denoms[0] == 'intensity')
  normfactors = "${setname}_normfacs"
  """
  # Create peptide table from PSM table, picking best scoring unique peptides
  msstitch peptides -i psms -o "${setname}_peptides" --scorecolpattern svm --spectracol 1 --modelqvals \
    ${do_ms1 ? '--ms1quantcolpattern area' : ''} \
    ${setisobaric ? "--isobquantcolpattern ${setisobaric} --minint 0.1" : ''} \
    ${keepnapsms_quant ? '--keep-psms-na-quant' : ''} \
    ${denoms && denoms[0] == 'sweep' ? '--mediansweep --logisoquant' : ''} \
    ${denoms && denoms[0] == 'intensity' ? '--medianintensity' : ''} \
    ${denoms && !specialdenom ? "--logisoquant --denompatterns ${denoms.join(' ')}" : ''} \
    ${normalize_isob ? "--median-normalize" : ''}
  ${normalize_isob ? "sed 's/^/$setname'\$'\t/' < normalization_factors_psms > '$normfactors'" : ""}
  """
}


process proteinGeneSymbolTableFDR {
 
  label 'msstitch'
 
  input:
  tuple val(setname), path('tpeptides'), path('tpsms'), path('dpeptides'), path(tfasta), path(dfasta), val(acctype), val(do_ms1), val(isobaric), val(denom), val(keepnapsms_quant), val(normalize)

  output:
  tuple val(setname), val(acctype), path("${setname}_protfdr"), emit: tables
  tuple val(acctype), path(normfactors), emit: normfacs optional true
  path('warnings'), emit: warnings optional true

  script:
  scorecolpat = acctype == 'proteins' ? '^q-value$' : 'linear model'
  specialdenom = denom && (denom[0] == 'sweep' || denom[0] == 'intensity')
  normfactors = "${setname}_normfacs"
  """
  # score col is linearmodel_qval or q-value, but if the column only contains 0.0 or NA (no linear modeling possible due to only q<10e-04), we use svm instead
  tscol=\$(head -1 tpeptides| tr '\\t' '\\n' | grep -n "${scorecolpat}" | cut -f 1 -d':')
  dscol=\$(head -1 dpeptides| tr '\\t' '\\n' | grep -n "${scorecolpat}" | cut -f 1 -d':')
  if [ -n "\$(cut -f \$tscol tpeptides| tail -n+2 | egrep -v '(NA\$|0\\.0\$)')" ] && [ -n "\$(cut -f \$dscol dpeptides| tail -n+2 | egrep -v '(NA\$|0\\.0\$)')" ]
    then
      scpat="${scorecolpat}"
      logflag="--logscore"
    else
      scpat="svm"
      logflag=""
      echo 'Not enough q-values or linear-model q-values for peptides to calculate FDR for ${acctype} of set ${setname}, using svm score instead to calculate FDR.' >> warnings
  fi
  msstitch ${acctype} -i tpeptides --decoyfn dpeptides -o "${setname}_protfdr" --scorecolpattern "\$scpat" \$logflag \
    ${acctype != 'proteins' ? "--fdrtype picked --targetfasta '$tfasta' --decoyfasta '$dfasta' ${params.fastadelim ? "--fastadelim '${params.fastadelim}' --genefield '${params.genefield}'": ''}" : ''} \
    ${do_ms1 ? '--ms1quant' : ''} ${isobaric ? "--isobquantcolpattern ${isobaric} --minint 0.1" : ''} \
    ${do_ms1 || isobaric ? '--psmtable tpsms' : ''} \
    ${keepnapsms_quant ? '--keep-psms-na-quant' : ''} \
    ${denom && denom[0] == 'sweep' ? '--mediansweep --logisoquant' : ''} \
    ${denom && denom[0] == 'intensity' ? '--medianintensity' : ''} \
    ${denom && !specialdenom ? "--denompatterns ${denom.join(' ')} --logisoquant" : ''} \
    ${normalize ? "--median-normalize" : ''}
    ${normalize ? "sed 's/^/$setname'\$'\t/' < normalization_factors_tpsms > '$normfactors'" : ""}
  """
}


process sampleTableCheckClean {

  // Runs no python but that container has the tools needed
  label 'python'
 
  input:
  tuple path('sampletable'), val(do_deqms)

  output:
  tuple path('clean_sampletable'), path('sampletable_no_special_chars')
  
  script:
  """
  # First add NO__GROUP marker for no-samplegroups clean sampletable from special chars
  awk -v FS="\\t" -v OFS="\\t" \'{if (NF==3) print \$1,\$2,\$3,"NO__GROUP"; else print}\' sampletable > clean_sampletable
  # Check if there are samplegroups at all
  ${do_deqms ? 'grep -v NO__GROUP clean_sampletable || (>&2 echo "Cannot run DEqMS without specified sample groups" && exit 1)': ''}
  # Count amount samples per group and error on group with only one sample
  ${do_deqms ? "grep -v NO__GROUP clean_sampletable | cut -f 4 | sort | uniq -c | sed 's/\\s*//' | grep '^1 ' && (>&2 echo 'Cannot run DEqMS when any sample groups have only one sample, please review input' && exit 1)" : ''}
  # strip lead/trail space in set name 
  paste <(cut -f1 clean_sampletable) <(cut -f2 clean_sampletable | sed "s/^\\s*//;s/\\s*\$//") <(cut -f3-4 clean_sampletable) > nowhitespace && mv nowhitespace clean_sampletable
  # substitute other weird characters
  sed "s/[^A-Za-z0-9_\\t]/_/g" clean_sampletable > sampletable_no_special_chars
  """
}


process proteinPeptideSetMerge {

  label 'msstitch'

  input:
  tuple val(setnames), val(acctype), path(tables), path(lookup), path(sampletable_with_special_chars), path('sampletable_no_special_chars'), val(do_isobaric), val(do_ms1), val(proteinconflvl), val(do_pgroup), val(do_deqms)
  
  output:
  tuple val(acctype), path('grouptable'), emit: with_nogroup
  tuple val(acctype), path('no_nogrouptable'), emit: nogroup_rm

  script:
  sampletable_iso = sampletable_with_special_chars.name != 'NO__FILE' && do_isobaric
  """

  # SQLite lookup needs copying to not modify the input file which would mess up a rerun with -resume
  cat $lookup > db.sqlite
  msstitch merge -i ${listify(tables).collect() { "$it" }.join(' ')} --setnames ${setnames.collect() { "'$it'" }.join(' ')} --dbfile db.sqlite -o mergedtable \
    --fdrcolpattern '^q-value\$' ${acctype != 'peptides' ? "--mergecutoff ${proteinconflvl}" : ''} \
    ${acctype == 'peptides' ? "--pepcolpattern 'peptide PEP'" : ''} \
    ${do_ms1 ? "--ms1quantcolpattern area" : ''} \
    ${do_isobaric ? "--isobquantcolpattern plex" : ''} \
    ${do_pgroup ? "--no-group-annotation" : ''}
   
  # Put annotation on header, use normal setname for finding, replace with clean name
  # "sed '0,/{RE}//{substitute}/..."  for only first line (0,/{RE} = read from 0 until RE,
  # then the empty // means use the previous RE (you could specify a new RE)
  head -n1 mergedtable > tmph
  sampletable_with_special_chars="$sampletable_with_special_chars"
  ${sampletable_iso ?
    'while read line ; do read -a arr <<< $line ; sed -E "s/^${arr[0]}_([a-z0-9]*plex)_${arr[1]}/${arr[4]}_${arr[3]}_${arr[2]}_\\1_${arr[1]}/" <(tail -n1 tmph | tr "\t" "\n") | tr "\n" "\t" | sed $"s/\\t$/\\n/" ; done < <(paste <(cut -f2 $sampletable_with_special_chars) sampletable_no_special_chars) >> tmph' :  ''}
  ${sampletable_iso ? "cat <(tail -n1 tmph) <(tail -n+2 mergedtable) > grouptable" : 'mv mergedtable grouptable'}

  # Remove internal no-group identifier so it isnt output
  sed '1s/NO__GROUP_//g' < grouptable > no_nogrouptable
  """
}


process DEqMS {
  input:
  tuple val(acctype), path('grouptable'), path('sampletable')

  output:
  tuple val(acctype), path('proteintable'), emit: with_nogroup
  tuple val(acctype), path('no_nogrouptable'), emit: nogroup_rm

  script:
  """
  # Run DEqMS if needed, use original sample table with NO__GROUP
  numfields=\$(head -n1 grouptable | tr '\t' '\n' | wc -l) && deqms.R && paste <(head -n1 grouptable) <(head -n1 deqms_output | cut -f \$(( numfields+1 ))-\$(head -n1 deqms_output|wc -w)) > tmpheader && cat tmpheader <(tail -n+2 deqms_output) > proteintable

  # Remove internal no-group identifier so it isnt output
  sed '1s/NO__GROUP_//g' < proteintable > no_nogrouptable
  """
}


def mzml_list = []
def header = []

workflow {

// Validate and set file inputs

// Files which are not standard can be checked here
if (params.hirief && !file(params.hirief).exists()) exit 1, "Peptide pI data file not found: ${params.hirief}"
if (params.hirief && !params.input) exit 1, "Cannot run HiRIEF delta pI calculation without fraction-annotated mzML definition file"
if (params.sampletable) {
  // create value channel with first()
  sampletable = Channel.fromPath(params.sampletable).first()
  if( !sampletable.exists() ) exit 1, "Sampletable file not found: ${params.sampletable}"
} else {
  sampletable = Channel.fromPath("${baseDir}/assets/NO__FILE").first()
}
  nofile = "${baseDir}/assets/NO__FILE"
  nofile_ch = Channel.fromPath(nofile)

  // Checking rerun status, if any
  complementary_run = params.targetpsmlookup && params.decoypsmlookup && params.targetpsms && params.decoypsms
  is_rerun = complementary_run && !params.input
  if (is_rerun) {
    fractionation_in = false
    mzml_in = Channel.empty()
    mzml_list = []
    mzmls = []
  } else {
    mzmls = msgf_info_map(params.input)
    mzml_list = mzmls.collect { k,v -> v}
    mzml_in = Channel.fromList(mzml_list)
    fractionation_in = mzml_list.any { it.fraction }
    if (fractionation_in) {
      println("Fractions detected in mzml definition, will run as fractionated")
      if (mzml_list.any { !it.fraction }) {
        println("Forcing files without fraction to file-named fraction in their specified plate")
      } else if (mzml_list.any { !it.plate }) {
        println("Forcing files without fraction to file-named fraction in set-named plate")
      }
    } else if (params.hirief && !(is_rerun || complementary_run)) {
      exit 1, "Cannot run HiRIEF --hirief while not specifying fractions"
    } else if (params.hirief) {
      println("No fractions detected in mzml definition but specified --hirief and supplied PSM table,will try to add pI data to any fractions detected in PSMs")
    } else {
      println("No fractions detected in mzml definition")
    }
  }
  if (complementary_run) {
    if (params.quantlookup) exit 1, "When specifying a complementary you may not pass --quantlookup"
    prev_results = Channel
      .fromPath([params.targetpsmlookup, params.decoypsmlookup, params.targetpsms, params.decoypsms, params.ptmpsms ?: nofile])
    def oldpsmheader
    new File("${params.targetpsms}").withReader { oldpsmheader = it.readLine() }
    old_fractionation = oldpsmheader.contains('Fraction')
  
  } else if (params.targetpsmlookup || params.decoypsmlookup || params.targetpsms || params.decoypsms || params.ptmpsms) {
    exit 1, "When specifying a complementary run you need to pass all of --targetpsmlookup, --decoypsmlookup, --targetpsms, --decoypsms"
  
  } else {
    old_fractionation = false
  }

  fractionation = fractionation_in || old_fractionation
  if (fractionation && complementary_run && (!params.oldmzmldef || !file(params.oldmzmldef).exists())) {
    exit 1, 'Fractionation with complementing run needs an --oldmzmldef file'
  }


  if (params.oldmzmldef) { 
    // oldmzml_lines = file("${params.oldmzmldef}").readLines().collect { it.tokenize('\t') }
    oldmzmls = msgf_info_map(params.oldmzmldef)
    oldmzml_sets = oldmzmls.collect { k,v -> v.setname }
    //oldmzmls_fn = Channel.fromPath(params.oldmzmldef).first()
  } else {
    //oldmzmls_fn = Channel.fromPath("${baseDir}/assets/NO__FILE").first()
    oldmzml_sets = []
  }

  all_setnames = (mzml_list.collect { it.setname } + oldmzml_sets).unique()

  // parse inputs that combine to form values or are otherwise more complex.
  // Isobaric input example: --isobaric 'set1:tmt10plex:127N:128N set2:tmt16plex:sweep set3:itraq8plex:intensity'
  isop = params.isobaric ? params.isobaric.tokenize(' ') : false
  setisobaric = isop ? isop.collect() {
    y -> y.tokenize(':')
  }.collectEntries() {
    x-> [x[0], x[1].replaceAll('tmtpro', 'tmt16plex')]
  } : [:]
  setdenoms = !params.noquant && isop ? isop.collect() {
    y -> y.tokenize(':')
  }.collectEntries() {
    x-> [x[0], x[2..-1]]
  } : [:]
  
  do_ms1 = !params.noquant && !params.noms1quant
  normalize = (!params.noquant && (params.normalize || params.deqms) && params.isobaric)

    //.tap { mzmlfiles_counter; mzmlfiles_qlup_sets } // for counting-> timelimits; getting sets from supplied lookup
    // .map { it -> [it[2].replaceAll('[ ]+$', '').replaceAll('^[ ]+', ''), file(it[0]).baseName.replaceAll(regex_specialchars, '_'), file(it[0]), it[1], plate_or_no(it, 3), fr_or_file(it, 4)] }

  mzml_in
    // Prepare mzml files (sort, collect) for processes that need all of them
    .toList()
    .map { it.sort( {a, b -> a.sample <=> b.sample}) } // sort on sample for consistent .sh script in -resume
    .map { it -> [it.collect() { it.setname }, it.collect() { it.mzmlfile }, it.collect() { it.plate }, it.collect() { it.fraction } ] } // lists: [sets], [mzmlfiles], [plates], [fractions]
    .map { it -> it[0..2] } // remove basenames, fractions
    .set { mzmlfiles_all_sort }

  // Spec lookup prep if needed
  do_quant = false
  if (is_rerun) {
    Channel.fromPath(params.targetpsmlookup)
    | PTMClean 
    | combine(Channel.fromPath(params.ptmpsms))
    | set { ptmpsms_lookup_ch }
    // For reporting:
    Channel.fromPath(params.targetpsmlookup)
    .map { [it, null] }
    .set { tdspeclookup }

  } else if (complementary_run) {
    mzmlfiles_all_sort
    | combine(prev_results.toList())
    | complementSpectraLookupCleanPSMs

    complementSpectraLookupCleanPSMs.out.dbs
    | map { it[0] }
    | PTMClean 
    | combine(complementSpectraLookupCleanPSMs.out.ptm)
    | set { ptmpsms_lookup_ch }

    complementSpectraLookupCleanPSMs.out.dbs
    | set { tdspeclookup }

    complementSpectraLookupCleanPSMs.out.psms
    | flatMap { [['target', it[0]], ['decoy', it[1]]] }
    | set { oldpsms_ch }
    
    do_quant = !params.noquant

  } else if (params.quantlookup) {
    // Runs with a premade quant lookup eg from previous search
    Channel
      .fromPath(params.quantlookup)
      .tap { ptmlookup_ch }
      .flatMap { [['target', it], ['decoy', it]] }
      .set { specquant_lookups }
    // For reporting
    ptmlookup_ch
    .map { [it, null] }
    .set { tdspeclookup }

  } else if (params.noquant && !params.quantlookup) {
    // Noquant, fresh spectra lookup scenario -> spec lookup ready for PSMs, PTMs
    mzmlfiles_all_sort
    | createNewSpectraLookup
    | flatMap { [['target', it], ['decoy', it]] }
    | set { specquant_lookups }
    createNewSpectraLookup.out
    .set { ptmlookup_ch }

  } else if (!params.quantlookup) {
    // Normal case - no rerun, fresh everything
    mzmlfiles_all_sort
    | createNewSpectraLookup
    | map { [it, it] }
    | set { tdspeclookup }
    do_quant = true
    createNewSpectraLookup.out
    .set { ptmlookup_ch }
  }

  if (do_quant) {
    mzml_in
    | filter { setisobaric[it.setname] }
    | map { [it.setname, it.mzmlfile, it.instrument] }
    | centroidMS1
    | map { it + setisobaric[it[0]] }
    | isobaricQuant
    | set { iso_processed }

    if (!params.noms1quant) {
      if (params.hardklor) {
        mzml_in
        | map { [it.sample, it.mzmlfile] }
        | combine(Channel.fromPath("$baseDir/assets/hardklor.conf"))
        | hardklor
        | kronik
        | set { ms1_q }
      } else {
        mzml_in
        | map { [it.sample, it.mzmlfile] }
        | dinosaur
        | set { ms1_q }
      }
    } else {
      mzml_in
      | map { [stripchars_infile(it.mzmlfile)[1], file(it.sample)] }
      | set { ms1_q }
    }
    mzml_in
    | filter { !setisobaric[it.setname] }
    | map { [stripchars_infile(it.mzmlfile)[1], file(it.sample)] }
    | concat(iso_processed)
    | join(ms1_q, remainder: true)
    | toList
    | map { it.sort({a, b -> a[0] <=> b[0]}) }
    | transpose
    | toList
    | combine(tdspeclookup.map { it[0] })
    | quantLookup
    | combine(tdspeclookup.map { it[1] })
    | flatMap { it -> [['target', it[0]], ['decoy', it[1]]] }
    | set { specquant_lookups }
  }

  tdb = Channel.fromPath(params.tdb)
  createTargetDecoyFasta(tdb)

  if (!is_rerun) {
    search_mods = [params.mods ? params.mods : false,
      params.ptms ?: false,
      params.locptms ?: false,
      ].findAll { it }
      .join(';')
    MSGFPERCO(mzml_in,
      createTargetDecoyFasta.out.concatdb,
      setisobaric,
      fractionation,
      mzmls,
      params.maxvarmods,
      params.msgfmods,
      search_mods,
      params.psmconflvl,
      params.pepconflvl,
      params.locptms,
      params.maxmiscleav,
      params.enzyme,
      params.minpeplen,
      params.maxpeplen,
      params.mincharge,
      params.maxcharge,
    )
  
    MSGFPERCO.out.t_tsv
    | ifEmpty(['target', nofile])
    | concat(MSGFPERCO.out.d_tsv)
    | groupTuple
    | join(specquant_lookups)
    | combine(createTargetDecoyFasta.out.bothdbs)
    | join(complementary_run ? oldpsms_ch : nofile_ch.flatMap { [['target', it], ['decoy', it]] })
    | map { it + [complementary_run, do_ms1 && it[0] == 'target', params.isobaric && it[0] == 'target', params.onlypeptides]}
    | createPSMTable
    createPSMTable.out.psmtable.set { psmtables_ch }
    createPSMTable.out.lookup.set { psmlookups_ch }
    createPSMTable.out.warnings.set { psmwarnings_ch }

  } else {
    // This is a rerun
    Channel.from([['target', file(params.targetpsmlookup)], ['decoy', file(params.decoypsmlookup)]])
      .set { psmlookups_ch }
    Channel.from([['target', file(params.targetpsms)], ['decoy', file(params.decoypsms)]])
      .set { psmtables_ch }
    psmwarnings_ch = Channel.empty()
  }

  if (params.hirief) {
    hiriefpep = Channel.fromPath(params.hirief)
    psmtables_ch
    | filter { it[0] == 'target' }
    | map { [it[1], params.strips]  }
    | combine(hiriefpep)
    | peptidePiAnnotation
    | map { ['target', it]}
    | set { target_psmtable }
  } else {
    psmtables_ch
    | filter { it[0] == 'target' }
    | set { target_psmtable }
  }
  psmtables_ch
  | filter { it[0] == 'decoy' }
  | concat(target_psmtable)
  | map { [it[0], it[1], all_setnames] }
  | splitPSMs
  | map{ it -> [it[0], listify(it[1]).collect() { it.baseName.replaceFirst(/\.tsv$/, "") }, it[1]]} // get setname from {setname}.tsv
  | transpose
  | set { splitpsms_ch }



  if (params.totalproteomepsms) {
    Channel.fromPath(params.totalproteomepsms)
    | map { [it] + all_setnames }
    | splitTotalProteomePSMs
    | map{ it -> [listify(it).collect() { it.baseName.replaceFirst(/\.tsv$/, "") }, it]} // get setname from {setname}.tsv
    | transpose
    | set { totalproteome_ch }
  } else {
    totalproteome_ch = Channel.empty()
  }
  if (params.genes) {
    totalprot_col = 'Gene Name'
  } else if (params.onlypeptides) {
   // FIXME sage
    totalprot_col = 'Protein'
  } else {
    totalprot_col = 'Master protein(s)'
  }

  if (params.locptms || params.ptms) {
    PTMANALYSIS(params.locptms ? params.locptms.tokenize(';') : [],
      params.ptms ? params.ptms.tokenize(';') : [],
      params.mods ? params.mods.tokenize(';') : [],
      all_setnames,
      setisobaric,
      mzml_in | map { [it.setname, it.mzmlfile] } | groupTuple,
      splitpsms_ch | filter { it[0] == 'target' } | map { it[1..-1] },
      tdb,
      !is_rerun ? MSGFPERCO.out.unfiltered : Channel.empty(),
      mzml_list,
      params.maxpeplen,
      params.maxcharge,
      params.msgfmods,
      params.minpsms_luciphor,
      params.ptm_minscore_high,
      complementary_run ? ptmpsms_lookup_ch : ptmlookup_ch.combine(nofile_ch),
      totalproteome_ch,
      totalprot_col,
      setisobaric,
      setdenoms,
      params.mediannormalize,
      params.keepnapsmsquant,
      do_ms1,
      !params.onlypeptides
    )
    ptm_ch = PTMANALYSIS.out.ptms
    ptmwarn_ch = PTMANALYSIS.out.warnings
  } else {
    ptm_ch = Channel.empty()
    ptmwarn_ch = Channel.empty()
  }

  splitpsms_ch
  | map { it + [it[0] == 'target' ? setisobaric[it[1]] : false, setdenoms[it[1]], params.keepnapsmsquant, it[0] == 'target' && params.mediannormalize, do_ms1 && it[0] == 'target']}
  | makePeptides

  acctypes = ['proteins']
  if (params.onlypeptides) {
    acctypes = []
  } else {
    if (params.ensg) {
    acctypes = acctypes.plus('ensg')
    }
    if (params.genes) {
    acctypes = acctypes.plus('genes')
    }
  }

  if (!params.onlypeptides) {
    makePeptides.out.peps
    | groupTuple
    | filter { it[1].size() == 2 } // T+D required! FIXME? if we can do without that demand this can be dropped
    | transpose
    | map { [it[1], it[0], it[2]] }
    | join(splitpsms_ch, by: [0, 1]) // now it's [td, setname, peptable, psmtable]
    | branch { t: it[0] == 'target'
               d: it[0] == 'decoy' }
    | set { tdpeps }

    tdpeps.t
    | map { it[1..-1] } // remove td
    | join(tdpeps.d | map { it[1..-2] }) // also strip dpsms from tdpeps.d
    | combine(createTargetDecoyFasta.out.bothdbs)
    | combine(Channel.from(acctypes))
    | map { it + [do_ms1, setisobaric[it[0]], setdenoms[it[0]], params.keepnapsmsquant, params.mediannormalize]}
    | proteinGeneSymbolTableFDR
  }

  if (params.sampletable) {
    Channel.fromPath(params.sampletable)
    | map { [it, params.deqms] }
    | sampleTableCheckClean
    | set { sampletable_ch }
  } else {
    nofile_ch
    .map { [it, it] }
    .set { sampletable_ch }
  }

  makePeptides.out.peps
  | filter { it[1] == 'target' }
  | map { [it[0], 'peptides', it[2]] }
  | concat(proteinGeneSymbolTableFDR.out.tables)
  | groupTuple(by: 1)
  | combine(psmlookups_ch.filter { it[0] == 'target' }.map { it[1] })
  | combine(sampletable_ch)
  | map { it + [do_quant && params.isobaric, do_ms1, params.proteinconflvl, !params.onlypeptides, params.deqms] }
  | proteinPeptideSetMerge

  if (params.deqms) {
    proteinPeptideSetMerge.out.with_nogroup
    | combine(sampleTableCheckClean.out)
    | DEqMS
    | set { protpepgene_ch }
  } else {
    proteinPeptideSetMerge.out
    | set { protpepgene_ch }
  }

  if (params.report_seqmatch) {
    MATCH_SEQUENCES(
      protpepgene_ch.nogroup_rm.filter { it[0] == 'peptides' }.map { it[1] },
      protpepgene_ch.nogroup_rm.filter { it[0] != 'peptides' },
      Channel.from(params.report_seqmatch).flatMap { it.tokenize(';') }.map { file(it) },
      params.maxmiscleav,
      params.minpeplen,
    ).set { feattables_out_ch }
  } else {
    protpepgene_ch.set { feattables_out_ch }
  }
  
  REPORTING(mzmlfiles_all_sort,
    tdspeclookup.map { it[0] },
    Channel.fromPath(params.input),
    Channel.fromPath(params.oldmzmldef ?: nofile),
    fractionation,
    target_psmtable.map { it[1] },
    protpepgene_ch.with_nogroup,
    makePeptides.out.normfacs.concat(proteinGeneSymbolTableFDR.out.normfacs),
    sampletable_ch,
    params.pepconflvl,
    params.proteinconflvl,
    all_setnames,
    ptm_ch,
    psmwarnings_ch
      .concat(!is_rerun ? MSGFPERCO.out.warnings: Channel.empty())
      .concat(ptmwarn_ch)
      .toList().toList()
      .filter { it[0] }
      .ifEmpty(nofile),
  )

  target_psmtable.map { it[1] }
  .concat(psmtables_ch | filter { it[0] == 'decoy' } | map { it[1] })
  .concat(psmlookups_ch | map { it[1] })
  .concat(ptm_ch.flatten())
  .concat(REPORTING.out.flatten())
  .subscribe { it.copyTo("${params.outdir}/${it.baseName}.${it.extension}") }
}


workflow.onComplete {
  if (workflow.success) {
    def libfile = file("${params.outdir}/libs.js")
    def libs = libfile.readLines()
    def bulma = file("${baseDir}/assets/bulma.js").readLines()
    def psmap = paramsSummaryMap(workflow)
    def files_header = read_header(params.input)
    files_header[0] = 'filename'
    files_header.remove('mzmlfile')
    infiles = mzml_list.collect { fn -> files_header.collect { fn[it] }}
    infiles.add(0, files_header)

    // Parse containers file to get software versions
    def jsonslurp = new JsonSlurper()
    def containers_versions = jsonslurp.parseText(new File("${baseDir}/containers.json").text)
      .collectEntries { k, v -> [v[workflow.containerEngine], [k, v.version]] }

    // Get processes used from trace to output the software versions used in pipeline
    def label_containers = psmap['Core Nextflow options']['container']
      .findAll { it.key.contains('withLabel:') }
      .collect { it.value }
      .collect { [containers_versions[it][0], containers_versions[it][1], it] }
    def sw_versions = label_containers + file("${params.outdir}/pipeline_info/execution_trace.txt").readLines()[1..-1]
      .collect { it.tokenize('\t')[3]
          .replaceFirst(~/^.+:/, '')
          .replaceFirst(~/\s\([0-9]+\)$/, '')
      }.unique()
      .collect { psmap['Core Nextflow options']['container'][it] }
      .findAll { it } // remove null in case process is not defined in container-config
      .unique()
      .collect { [containers_versions[it][0], containers_versions[it][1], it] }
      
    // Set name
    // if not wf.runName (-name or auto) is like "crazy_euler" or other "{adjective}_{scientist}"
    if (!params.name && !(workflow.runName ==~ /[a-z]+_[a-z]+/) ) {
      runname = workflow.runName
    } else if (!params.name) {
      runname = 'untitled'
    } else {
      runname = params.name
    }
    def fields = [runname: runname,
        sw_versions: sw_versions,
        params: psmap['Other parameters'],
        infiles: infiles,
        libs: libs, bulma: bulma]
    def rf = new File("${params.outdir}/report_groovy_template.html")
    def temp_engine = new groovy.text.StreamingTemplateEngine()
    def report_template = temp_engine.createTemplate(rf).make(fields)
    def report_html = report_template.toString()
    def output_rf = new File( params.outdir, "report.html" )
    output_rf.withWriter { w -> w << report_html }
    rf.delete()
    libfile.delete()
  }
}
