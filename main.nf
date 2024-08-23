#!/usr/bin/env nextflow

include { msgf_info_map; listify; stripchars_infile } from './modules.nf' 
include { MSGFPERCO } from './workflows/msgf_perco.nf'

/*
========================================================================================
                         lehtiolab/ddamsproteomics
========================================================================================
 lehtiolab/ddamsproteomics Analysis Pipeline.
 #### Homepage / Documentation
 https://github.com/lehtiolab/ddamsproteomics
----------------------------------------------------------------------------------------
*/


def multifile_format(fileparam) {
    if (!fileparam) {
        return false
    }
    sum_fn = file(fileparam)
    if (!(sum_fn instanceof List)) {
      sum_fn = [sum_fn]
    }
    return sum_fn.join(', ')
}


def create_workflow_summary(summary) {

    def yaml_file = workDir.resolve('workflow_summary_mqc.yaml')
    yaml_file.text  = """
    id: 'lehtiolab-ddamsproteomics-summary'
    description: " - this information is collected when the pipeline is started."
    section_name: 'lehtiolab/ddamsproteomics Workflow Summary'
    section_href: 'https://github.com/lehtiolab/ddamsproteomics'
    plot_type: 'html'
    data: |
        <dl class=\"dl-horizontal\">
${summary.collect { k,v -> "            <dt>$k</dt><dd><samp>${v ?: '<span style=\"color:#999999;\">N/A</a>'}</samp></dd>" }.join("\n")}
        </dl>
    """.stripIndent()

   return yaml_file
}


/*
 * Parse software version numbers
 */
process get_software_versions {

    publishDir "${params.outdir}", mode: 'copy', overwrite: true

    output:
    file 'software_versions.yaml' into software_versions_qc

    script:
    noms1 = params.noms1quant || params.noquant
    """
    echo $workflow.manifest.version > v_pipeline.txt
    echo $workflow.nextflow.version > v_nextflow.txt
    echo 2023.01.1202 > v_msgf.txt
    ${!noms1 && !params.hardklor ? 'dinosaur | head -n2 | grep Dinosaur > v_dino.txt || true' : ''}
    ${!noms1 && params.hardklor ? 'hardklor | head -n1 > v_hk.txt || true' : ''}
    kronik | head -n2 | tr -cd '[:digit:],\\.' > v_kr.txt || true
    #luciphor2 |& grep Version > v_luci.txt # incorrect version from binary (2014), echo below
    echo Version: 2020_04_03 > v_luci.txt # deprecate when binary is correct
    echo 3.5 > v_perco.txt
    msstitch --version > v_mss.txt
    echo 2.9.1 > v_openms.txt
    Rscript <(echo "packageVersion('DEqMS')") > v_deqms.txt
    scrape_software_versions.py > software_versions.yaml
    """
}



process createTargetDecoyFasta {
  container "${workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
    'https://depot.galaxyproject.org/singularity/msstitch:3.15--pyhdfd78af_0' :
    'quay.io/biocontainers/msstitch:3.15--pyhdfd78af_0'}"
 
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

def plate_or_no(it, length) {
  return it.size() > 3 ? it[3] : "no_plate"
}




// Parse mzML input to get files and sample names etc
// get setname, sample name (baseName), input mzML file. 
// Set platename to setname if not specified. 
// Set fraction name to NA if not specified

/*
* Step 1: Extract quant data from peptide spectra
*/


process centroidMS1 {
  container 'chambm/pwiz-skyline-i-agree-to-the-vendor-licenses:3.0.20066-729ef9c41'

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
  container "${workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
    'https://depot.galaxyproject.org/singularity/openms:2.9.1--h135471a_1' :
    'quay.io/biocontainers/openms:2.9.1--h135471a_1'}"

  input:
  tuple val(setname), val(parsed_infile), path(infile), val(instr), val(isobtype)

  output:
  tuple val(parsed_infile), path(outfile)

  script:
  outfile = "${infile.baseName}.consensusXML"
  activationtype = [auto: 'auto', any: 'any', hcd:'beam-type collision-induced dissociation', cid:'Collision-induced dissociation', etd:'Electron transfer dissociation'][params.activation]
  plextype = isobtype ? isobtype.replaceFirst(/[0-9]+plex/, "") : 'false'
  massshift = [tmt:0.0013, itraq:0.00125, false:0][plextype]
  //(is_stripped, parsed_infile) = stripchars_infile(infile)
  """
  ${isobtype ? "IsobaricAnalyzer -type $isobtype -in \"${infile}\" -out \"${infile.baseName}.consensusXML\" -extraction:select_activation \"$activationtype\" -extraction:reporter_mass_shift $massshift -extraction:min_precursor_intensity 1.0 -extraction:keep_unannotated_precursor true -quantification:isotope_correction true" : ''}
  """
}



process dinosaur {
container 'lehtiolab/ddamsproteomics:2.18'
// Biocontainers arent working with dinosaur tests, FIXME create auto build of it

//  container "${workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
//    'https://depot.galaxyproject.org/singularity/dinosaur:1.2.0--0' :
//    'quay.io/biocontainers/dinosaur:1.2.0--0'}"

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
  container "${workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
    'https://depot.galaxyproject.org/singularity/hardklor:2.3.2--he1b5a44_0' :
    'quay.io/biocontainers/hardklor:2.3.2--he1b5a44_0'}"

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
  container "${workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
    'https://depot.galaxyproject.org/singularity/kronik:2.20--h4ac6f70_6' :
    'quay.io/biocontainers/kronik:2.20--h4ac6f70_6'}"
  input:
  tuple val(sample), val(parsed_infile), path('hardklor.out')
  
  output:
  tuple val(parsed_infile), path("${sample}.features.tsv")
  
  script:
  """
  kronik -c 5 -d 3 -g 1 -m 8000 -n 600 -p 10 hardklor.out ${sample}.features.tsv
  """
}



process complementSpectraLookupCleanPSMs {

  input:
  tuple val(in_setnames), path(mzmlfiles), val(platenames) path(tlup), path(dlup), path(tpsms), path(dpsms), path(ptmpsms)

  output:
  tuple path('t_cleaned_psms.txt'), path('d_cleaned_psms.txt') into cleaned_psms
  tuple path('target_db.sqlite'), path('decoy_db.sqlite'), emit: dbs
  path 'cleaned_ptmpsms.txt' into cleaned_ptmpsms optional true
  path 'ptms_db.sqlite' into ptm_lookup_old optional true
  file('all_setnames') into oldnewsets 
  
  script:
  setnames = in_setnames.unique(false)
  """
  # If this is an addition to an old lookup, copy it and extract set names
  cp ${tlup} target_db.sqlite
  cp ${dlup} decoy_db.sqlite
  sqlite3 target_db.sqlite "SELECT set_name FROM biosets" > old_setnames
  # If adding to old lookup: grep new setnames in old and run msstitch deletesets if they match
  # use -x for grep since old_setnames must grep whole word
  if grep -xf old_setnames <(echo ${setnames.join('\n')} )
    then
      msstitch deletesets -i ${tpsms} -o t_cleaned_psms.txt --dbfile target_db.sqlite --setnames ${setnames.collect() { "'${it}'" }.join(' ')}
      msstitch deletesets -i ${dpsms} -o d_cleaned_psms.txt --dbfile decoy_db.sqlite --setnames ${setnames.collect() { "'${it}'" }.join(' ')}
      ${params.ptmpsms ? "msstitch deletesets -i ${ptmpsms} -o cleaned_ptmpsms.txt --setnames ${setnames.collect() {"'${it}'"}.join(' ')}" : ''}
    else
      mv ${tpsms} t_cleaned_psms.txt
      mv ${dpsms} d_cleaned_psms.txt
      ${params.ptmpsms ? "mv ${ptmpsms} cleaned_ptmpsms.txt" : ''}
  fi
  
  ${mzmlfiles.collect() { stripchars_infile(it, return_oldfile=true) }.findAll{ it[0] }.collect() { "ln -s '${it[2]}' '${it[1]}'" }.join(' && ')}
  ${mzmlfiles.size() ? "msstitch storespectra --spectra ${mzmlfiles.collect() { "'${it.toString().replaceAll('[&<>\'"]', '_')}'" }.join(' ')} --setnames ${in_setnames.collect() { "'$it'" }.join(' ')} --dbfile target_db.sqlite" : ''}
  ${params.ptmpsms ? "cp target_db.sqlite ptms_db.sqlite" : ''}

  copy_spectra.py target_db.sqlite decoy_db.sqlite ${params.ptmpsms ? 'ptms_db.sqlite' : '0'} ${setnames.join(' ')}
  cat old_setnames <(echo ${setnames.join('\n')}) | sort -u | grep -v '^\$' > all_setnames
  """
}


process createNewSpectraLookup {
  container "${workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
    'https://depot.galaxyproject.org/singularity/msstitch:3.15--pyhdfd78af_0' :
    'quay.io/biocontainers/msstitch:3.15--pyhdfd78af_0'}"

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
  container "${workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
    'https://depot.galaxyproject.org/singularity/msstitch:3.15--pyhdfd78af_0' :
    'quay.io/biocontainers/msstitch:3.15--pyhdfd78af_0'}"

  // FIXME publishDir "${params.outdir}", mode: 'copy', overwrite: true, saveAs: {it == 'target.sqlite' ? 'quant_lookup.sql' : null }

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
  container "${workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
    'https://depot.galaxyproject.org/singularity/msstitch:3.16--pyhdfd78af_0' :
    'quay.io/biocontainers/msstitch:3.16--pyhdfd78af_0'}"

  input:
  tuple val(td), path(psms), path('lookup'), path(tdb), path(ddb), path('oldpsms'), val(is_rerun), val(complementary_run), val(do_ms1), val(do_isobaric), val(onlypeptides)

  output:
  tuple val(td), file("${outpsms}"), emit: psmtable
  path("${psmlookup}"), emit: lookup
  path('warnings'), emit: warnings optional true

  script:
  psmlookup = "${td}_psmlookup.sql"
  outpsms = "${td}_psmtable.txt"
  no_target = td == 'target' && psms.size() == 0
  no_decoy = td == 'decoy' && psms.size() == 0
  """
  ${no_target && !is_rerun ? "echo 'No target PSMs made the combined PSM / peptide FDR cutoff' && exit 1" : ''}
  ${no_decoy && !is_rerun ? "echo 'No decoy PSMs in any set at FDR cutoff, will not produce protein/gene tables' > warnings && touch ${outpsms} && touch ${psmlookup} && exit 0" : ''}
  ${!is_rerun ? "msstitch concat -i ${listify(psms).collect() {"$it"}.join(' ')} -o psms.txt" : ''}
  # SQLite lookup needs copying to not modify the input file which would mess up a rerun with -resume
  cat lookup > $psmlookup
  ${!is_rerun ? "sed '0,/\\#SpecFile/s//SpectraFile/' -i psms.txt": ''}
# --addmiscleav 
  ${is_rerun ? "mv oldpsms psmsrefined" : "msstitch psmtable -i psms.txt --dbfile $psmlookup -o ${outpsms} \
      ${onlypeptides ? '' : "--fasta \"${td == 'target' ? "${tdb}" : "${ddb}"}\" --genes"} \
      ${do_ms1 ? '--ms1quant' : ''} \
      ${do_isobaric ? "--isobaric --min-precursor-purity ${params.minprecursorpurity}" : ''} \
      ${!onlypeptides ? '--proteingroup' : ''} \
      ${complementary_run ? '--oldpsms oldpsms' : ''}"
  }
  # sed 's/\\#SpecFile/SpectraFile/' -i psmsrefined
  # In decoy PSM table process, also split the target total proteome normalizer table if necessary.
  # Doing it in decoy saves time, since target is usally largest table and slower
  """
}


process peptidePiAnnotation {
  container "python:3.12"

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
  container "${workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
    'https://depot.galaxyproject.org/singularity/msstitch:3.16--pyhdfd78af_0' :
    'quay.io/biocontainers/msstitch:3.16--pyhdfd78af_0'}"

  input:
  tuple val(td), path('psms'), val(setnames)

  output:
  tuple val(td), file({setnames.collect() { "${it}.tsv" }}) optional true

  script:
  """
echo ${setnames}
  msstitch split -i psms --splitcol bioset
  """
}


process splitTotalProteomePSMs {
  container "${workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
    'https://depot.galaxyproject.org/singularity/msstitch:3.16--pyhdfd78af_0' :
    'quay.io/biocontainers/msstitch:3.16--pyhdfd78af_0'}"

  input:
  tuple path('tppsms_in'), val(setnames)

  output:
  path({setnames.collect() { "tppsms/${it}.tsv" }}), emit: totalprotpsms_allsets optional true

  script:
  """
  mkdir tppsms && msstitch split -i tppsms_in -d tppsms --splitcol bioset
  """
}


process makePeptides {
  container "${workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
    'https://depot.galaxyproject.org/singularity/msstitch:3.16--pyhdfd78af_0' :
    'quay.io/biocontainers/msstitch:3.16--pyhdfd78af_0'}"

  input:
  tuple val(td), val(setname), path('psms'), val(setisobaric), val(denoms), val(keepnapsms_quant), val(normalize_isob), val(do_ms1)
  
  output:
  //set val(setname), val(td), file(psms), file("${setname}_peptides") into prepgs_in
  //set val(setname), val('peptides'), val(td), file("${setname}_peptides"), path(normfactors) optional true into peptides_out

  script:
//  quant = !params.noquant && td == 'target'
//  isoquant = quant && setisobaric && setisobaric[setname]
//  denom = isoquant && setdenoms ? setdenoms[setname] : false
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
  ${normalize_isob ? "sed 's/^/$setname'\$'\t/' < normalization_factors_psms > '$normfactors'" : "touch '$normfactors'"}
  """
}


/* FIXME remove
process countMS2perFile {

  input:
  set val(setnames), file(mzmlfiles), val(platenames), file(speclookup) from specfilein

  output:
  set val(setnames), file(mzmlfiles), val(platenames), file('amount_spectra_files') into specfilems2

  script:
  """
  sqlite3 $speclookup "SELECT set_name, mzmlfilename, COUNT(*) FROM mzml JOIN mzmlfiles USING(mzmlfile_id) JOIN biosets USING(set_id) GROUP BY mzmlfile_id" > amount_spectra_files
  """
}
*/


process countMS2sPerPlate {
  container "python:3.12",

// FIXME this could possibly go into the PSM QC code?

  input:
  tuple val(setnames), file(mzmlfiles), val(platenames), val(fractionation), path(oldmzmls_fn), val(regex_specialchars)

  output:
  tuple path('amount_spectra_files'), path('scans_per_plate'), emit: counted
  path('allplates'), emit: allplates

  script:
  splates = [setnames, platenames].transpose().collect() { "${it[0]}_${it[1]}" }
  """
  #!/usr/bin/env python
  import os
  import re
  import sqlite3
  # FIXME isnt a plate NA?
  con = sqlite3.Connection("$speclookup")
  cursor = con.execute("SELECT set_name, mzmlfilename, COUNT(*) FROM mzml JOIN mzmlfiles USING(mzmlfile_id) JOIN biosets USING(set_id) GROUP BY mzmlfile_id")
  if ${fractionation}:
      platesets = [\"${splates.join('", "')}\"]
      plates = [\"${platenames.join('", "')}\"]
      setnames = [\"${setnames.join('", "')}\"]
      fileplates = {}
      for fn, setname, plate in zip([${mzmlfiles.collect() { "'${stripchars_infile(it)[1]}'" }.join(',')}], setnames, plates):
          try:
              fileplates[fn][setname] = plate
          except KeyError:
              fileplates[fn] = {setname: plate} 
      if ${complementary_run ? 1 : 0}:
          with open('$oldmzmls_fn') as oldmzfp:
              for line in oldmzfp:
                  if line.strip('\\n').split('\\t')[0] == 'mzmlfile':
                      continue
                  fpath, inst, setname, plate, fraction = line.strip('\\n').split('\\t')
                  # old mzmls also contain files that are no longer used (i.e. removed from set)
                  # filter by skipping any setname that is in the current new setnames
                  if setname in setnames:
                      continue
                  setplate = '{}_{}'.format(setname, plate)
                  fn = re.sub('${regex_specialchars.replaceAll("'", "\\\\'")}', '_', os.path.basename(fpath))
                  if setplate not in platesets:
                      platesets.append(setplate)
                  if fn not in fileplates:
                      fileplates[fn] = {setname: plate}
                  elif setname not in fileplates[fn]:
                      fileplates[fn][setname] = plate
      with open('allplates', 'w') as fp:
          fp.write('\\n'.join(platesets))
      platescans = {p: 0 for p in platesets}
      with open('amount_spectra_files', 'w') as fp:
          for setname, fn, scans in cursor:
              setplate = '{}_{}'.format(setname, fileplates[fn][setname])
              platescans[setplate] += int(scans)
      with open('scans_per_plate', 'w') as fp:
          for plate, scans in platescans.items():
              fp.write('{}\\t{}\\n'.format(plate, scans))
  else:
      # make sure there's an empty file for scans_per_plate if not fractionated
      with open('amount_spectra_files', 'w') as fp, open('scans_per_plate', 'w') as _, open('allplates', 'w') as __:
          for line in cursor:
              fp.write('\\t'.join(line))
              fp.write('\\n')
  """
}



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
    Channel
      .empty()
      .set { mzml_in }
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
      .fromPath([params.targetpsmlookup, params.decoypsmlookup, params.targetpsms, params.decoypsms, params.ptmpsms ? params.ptmpsms : 'NA'])
    def oldpsmheader
    new File("${params.targetpsms}").withReader { oldpsmheader = it.readLine() }
    old_fractionation = oldpsmheader.contains('Fraction')
  
  } else if (params.targetpsmlookup || params.decoypsmlookup || params.targetpsms || params.decoypsms || params.ptmpsms) {
    exit 1, "When specifying a complementary run you need to pass all of --targetpsmlookup, --decoypsmlookup, --targetpsms, --decoypsms"
  
  } else {
    old_fractionation = false
    prev_results = Channel.fromPath([nofile, nofile, nofile, nofile, nofile])
  }

  fractionation = fractionation_in || old_fractionation
  if (fractionation && complementary_run && (!params.oldmzmldef || !file(params.oldmzmldef).exists())) {
    exit 1, 'Fractionation with complementing run needs an --oldmzmldef file'
  }

  // parse inputs that combine to form values or are otherwise more complex.
  // Isobaric input example: --isobaric 'set1:tmt10plex:127N:128N set2:tmt16plex:sweep set3:itraq8plex:intensity'
  isop = params.isobaric ? params.isobaric.tokenize(' ') : false
  setisobaric = isop ? isop.collect() {
    y -> y.tokenize(':')
  }.collectEntries() {
    x-> [x[0], x[1].replaceAll('tmtpro', 'tmt16plex')]
  } : [false: 1]
  setdenoms = isop ? isop.collect() {
    y -> y.tokenize(':')
  }.collectEntries() {
    x-> [x[0], x[2..-1]]
  } : [false:1]
  
  normalize = (!params.noquant && (params.normalize || params.deqms) && params.isobaric)

    //.tap { mzmlfiles_counter; mzmlfiles_qlup_sets } // for counting-> timelimits; getting sets from supplied lookup
    // .map { it -> [it[2].replaceAll('[ ]+$', '').replaceAll('^[ ]+', ''), file(it[0]).baseName.replaceAll(regex_specialchars, '_'), file(it[0]), it[1], plate_or_no(it, 3), fr_or_file(it, 4)] }

  mzml_in
// FIXME clean out
    // Prepare mzml files (sort, collect) for processes that need all of them
    .toList()
    .map { it.sort( {a, b -> a.sample <=> b.sample}) } // sort on sample for consistent .sh script in -resume
    .map { it -> [it.collect() { it.setname }, it.collect() { it.mzmlfile }, it.collect() { it.plate }, it.collect() { it.fraction } ] } // lists: [sets], [mzmlfiles], [plates], [fractions]
    //.tap { mzmlfiles_psmqc } // FIXME this prevents completion of the pipeline somehow
    .map { it -> it[0..2] } // remove basenames, fractions
    .set { mzmlfiles_all_sort }

  // Collect all mzMLs into single item to pass to lookup builder and spectra counter
  if (params.oldmzmldef) { 
    oldmzml_lines = file("${params.oldmzmldef}").readLines().collect { it.tokenize('\t') }
    Channel
      .from(oldmzml_lines)
      .map { it -> it[2].replaceAll('[ ]+$', '').replaceAll('^[ ]+', '') }
      .unique()
      .toList()
      .set { oldmzml_sets }
    oldmzmls_fn = Channel.fromPath(params.oldmzmldef).first()
  } else {
    oldmzmls_fn = Channel.fromPath("${baseDir}/assets/NO__FILE").first()
    oldmzml_sets = Channel.value([])
  }


/*
// FIXME this needed still?
  if (!is_rerun) {
    mzml_in 
      .count()
      .subscribe { println "$it mzML files in analysis" }
      .into { mzmlcount_psm; mzmlcount_percolator }
  } else {
    Channel.value(0).into { mzmlcount_psm; mzmlcount_percolator }
  }
*/

  // Spec lookup prep if needed
  if (complementary_run) {
    mzmlfiles_all_sort
    | combine(prev_results.toList())
    | complementSpectraLookupCleanPSMs
    complementSpectraLookupCleanPSMs.out.dbs
    | set { tdspeclookup }

  } else if (!params.quantlookup) {
    mzmlfiles_all_sort
    | createNewSpectraLookup
    | map { [it, it] }
    | set { tdspeclookup }
  }

  // Quant lookup preparing if needed
  if (params.noquant && !params.quantlookup) {
    // Noquant, fresh spectra lookup scenario -> spec lookup ready for PSMs, PTMs
    tdspeclookup
      .flatMap { it -> [['target', it], ['decoy', it]] }
      .set { specquant_lookups }
  
  } else if (params.quantlookup) {
    // Runs with a premade quant lookup eg from previous search
    Channel
      .fromPath(params.quantlookup)
      .flatMap { it -> [['target', it], ['decoy', it]] }
      .set { specquant_lookups }

  } else if (is_rerun) {
// FIXME how to deal with these in reruns?
    // In case of rerun with same sets, no new search but only some different post-identification
    // output options
    //prespectoquant  // contains  t, d lookups
    //  .flatMap { it -> [[it[0], 'target'], [it[1], 'decoy']] }
    //  .set{ specquant_lookups }

  } else if (complementary_run) {
    // FIXME wtf to do here?

  } else {
    // Normal case - no rerun, fresh everythin
    //if (setisobaric) {
      mzml_in
      | filter { setisobaric[it.setname] }
      | map { [it.setname, it.mzmlfile, it.instrument] }
      | centroidMS1
      | map { it + setisobaric[it[0]] }
      | isobaricQuant
      | set { iso_processed }
    //} else {
      mzml_in
      | filter { !setisobaric[it.setname] }
      | map { [stripchars_infile(it.mzmlfile)[1], file(it.sample)] }
      | concat(iso_processed)
      | set { iso_q }
    //}
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
    iso_q
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
  // bothdbs.into { psmdbs; fdrdbs; ptmdbs }
  MSGFPERCO(mzml_in, createTargetDecoyFasta.out.concatdb, setisobaric, fractionation, mzmls)

  do_ms1 = !params.noquant && !params.noms1quant
  all_setnames = [mzml_list.collect { it.setname }.unique()]

  MSGFPERCO.out.t_tsv
  | ifEmpty(['target', []])
  | concat(MSGFPERCO.out.d_tsv)
  | join(specquant_lookups)
  | combine(createTargetDecoyFasta.out.bothdbs)
  | combine(nofile_ch) // old psms FIXME
  | map { it + [is_rerun, complementary_run, do_ms1 && it[0] == 'target', params.isobaric && it[0] == 'target', params.onlypeptides]}
  | createPSMTable

  if (params.hirief) {
    hiriefpep = Channel.fromPath(params.hirief)
    createPSMTable.out.psmtable
    | filter { it[0] == 'target' }
    | map { [it[1], params.strips]  }
    | combine(hiriefpep)
    | peptidePiAnnotation
    | map { ['target', it] + all_setnames }
    | set { target_psmtable }
  } else {
    createPSMTable.out.psmtable
    | filter { it[0] == 'target' }
    | map { it + all_setnames }
    | set { target_psmtable }
  }
  createPSMTable.out.psmtable
  | filter { it[0] == 'decoy' }
  | map { it + all_setnames }
  | concat(target_psmtable)
  | splitPSMs

  if (params.totalproteomepsms) {
    Channel.fromPath(params.totalproteomepsms)
    | map { [it] + all_setnames }
    | splitTotalProteomePSMs
  }

  splitPSMs.out
  | map { it -> [it[0], listify(it[1])] }
  | map{ it -> [it[0], it[1].collect() { it.baseName.replaceFirst(/\.tsv$/, "") }, it[1]]} // get setname from {setname}.tsv
  | transpose()
  | map { it + [it[0] == 'target' ? setisobaric[it[1]] : false, setdenoms[it[1]], params.keepnapsmsquant, params.mediannormalize, do_ms1 && it[0] == 'target']}
  | makePeptides
  //tuple val(td), val(setname), path('psms'), val(setisobaric), val(denoms), val(keepnapsms_quant), val(normalize_isob), val(do_ms1)
  //out: tuple val(td), file({setnames.collect() { "${it}.tsv" }}), emit: set_tables optional true



/*

  // QC for PSMs
  mzmlfiles_all_sort
  | combine(tdspeclookup.map { it[0] })
  | combine(Channel.from(fractionation ? 1 : 0))
  | combine(oldmzmls_fn)
  | countMS2sPerPlate

  // FIXME allplates could go as file into R QC directly?
  if (fractionation) {
    countMS2sPerPlate.out.allplates
      .splitText()
      .map { it -> it.trim() }
      .toList()
      .map { it -> [it] }
      .set { allplates_split }
    countMS2sPerPlate.out.counted
      .combine(allplates_split)
      .set { scans_result }
  } else {
    countMS2sPerPlate.out.counted
    | map { it -> [it[0], 'NA', ['noplates']] }
    | set { scans_result }
    //| into { scans_platecount; scans_result }
  }
  psm_result // From createPSMTable
    .filter { it[0] == 'target' }
    .combine(scans_result)
    .map { it -> [it[0], it[1], it[2], it[3], it[4].unique()] }
    .set { targetpsm_result }




mzmlfiles_all_count
// Count lookup can be spectra if needed quant
  .merge(countlookup)
  .set { specfilein }


process countMS2perFile {

  input:
  set val(setnames), file(mzmlfiles), val(platenames), file(speclookup) from specfilein

/// The following is maybe OK?

if (complementary_run) {
  oldnewsets
    .splitText()
    .map { it -> it.trim() }
    .toList()
    .set { allsetnames }
  cleaned_psms
    .flatMap { it -> [['target', it[0]], ['decoy', it[1]]] }
    .set { td_oldpsms }
} else {
  mzml_in
    .map { it -> it[2] } 
    .unique()
    .toList()
    .set { allsetnames }

  // if not using this youll have a combine on an open channel without
  // anything from complement cleaner. Will not run createPTMLookup then
  cleaned_ptmpsms = Channel.value('NA')
}

/// until here - we can investigate later

// Set names are first item in input lists, collect them for PSM tables and QC purposes
// FIXME just reuse channels
allsetnames 
  .into { setnames_featqc; setnames_psms; setnames_psmqc; setnames_ptms }



if (!params.quantlookup) {
  ptm_lookup_old
    .concat(ptm_lookup_new)
// PTM lookup needs not to have quant, can be spectra
    .set { ptm_lookup_in }
}

*/
  
//  publishDir "${params.outdir}", mode: 'copy', overwrite: true, saveAs: {["target_psmlookup.sql", "decoy_psmlookup.sql", "target_psmtable.txt", "decoy_psmtable.txt"].contains(it) && (is_rerun || !no_psms) ? it : null}

  target_psmtable.map { it[1] }
  .concat(createPSMTable.out.psmtable | filter { it[0] == 'decoy' } | map { it[1] })
  .concat(createPSMTable.out.lookup)
  .subscribe { it.copyTo("${params.outdir}/${it.baseName}.${it.extension}") }
}

