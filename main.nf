#!/usr/bin/env nextflow
/*
========================================================================================
                         lehtiolab/ddamsproteomics
========================================================================================
 lehtiolab/ddamsproteomics Analysis Pipeline.
 #### Homepage / Documentation
 https://github.com/lehtiolab/ddamsproteomics
----------------------------------------------------------------------------------------
*/


def helpMessage() {
    log.info"""
    =========================================
     lehtiolab/ddamsproteomics v${workflow.manifest.version}
    =========================================
    Usage:

    The typical command for running the pipeline is as follows:

    nextflow run lehtiolab/ddamsproteomics --mzmls '*.mzML' --tdb swissprot_20181011.fa --mods 'oxidation;carbamidomethyl', --locptms 'phospho' -profile standard,docker

    Mandatory arguments:
      --mzmls                       Path to mzML files
      --mzmldef                     Alternative to --mzml: path to file containing list of mzMLs 
                                    with instrument, sample set and fractionation annotation (see docs)
      --tdb                         Path to target FASTA protein databases, can be (quoted) '/path/to/*.fa'
      -profile                      Configuration profile to use. Can use multiple (comma separated)
                                    Available: standard, conda, docker, singularity, awsbatch, test

    Options:
      --instrument                  If not using --mzmldef, use this to specify instrument type.
                                    Currently supporting 'qe', 'qehf', 'velos', 'lumos', 'qehfx', 'lowres', 'timstof'

      --mods                        Modifications specified by their UNIMOD name. e.g. --mods 'oxidation;carbamidomethyl'
                                    Note that there are a limited number of modifications available, but that
                                    this list can easily be expanded in assets/msgfmods.txt
      --locptms                     As for --mods, but specify labile mods, pipeline will output false localization rate,
                                    output QC plots and if handed total proteome PSMs also output a normalized peptide table
                                    for deciding if PTM peptide expression deviates from respective protein expression
                                    e.g. --locptms 'phospho'
      --ptms                        As for --locptms but without false localization rate, for non-labile PTMs 
                                    e.g. --ptms 'acetyl'
      --totalproteomepsms           A PSM table of a global (total) proteome search without PTMs, to do a normalization
                                    of the PTM peptide table against. PTM peptides will be normalized to their matching
                                    proteins for isobaric quant.
      --isobaric VALUE              In case of isobaric, specify per set the type and possible denominators/sweep/intensity.
                                    In case of intensity, no ratios will be output but instead the raw PSM intensities will be
                                    median-summarized to the output features (e.g. proteins).
                                    Available types are tmt18plex, tmt16plex, tmtpro (=16plex), tmt10plex,
                                    tmt6plex, itraq8plex, itraq4plex
                                    E.g. --isobaric 'set1:tmt10plex:126:127N set2:tmt16plex:127C:131 set3:tmt10plex:sweep set4:itraq8plex:intensity'
      --psmconflvl                  Cutoff for PSM FDR on PSM table, default is 0.01
      --pepconflvl                  Cutoff for peptide FDR on PSM table, default is 0.01
      --proteinconflvl              Cutoff for protein/gene FDR in respective output tables, default is 0.01
      --activation VALUE            Specify activation protocol for isobaric quantitation filtering (NOT for identification):
                                    choose from auto (DEFAULT), hcd, cid, etd, any
      --fastadelim VALUE            FASTA header delimiter in case non-standard FASTA is used, to be used with
                                    --genefield
      --genefield VALUE             Number to determine in which field of the FASTA header (split 
                                    by --fastadelim) the gene name can be found.

      SEARCH ENGINE DETAILED PARAMETERS
      --prectol                     Precursor error for search engine (default 10ppm)
      --iso_err                     Isotope error for search engine (default -1,2)
      --frag                        Fragmentation method for search engine (default 'auto')
      --enzyme                      Enzyme used, default trypsin, pick from:
                                    unspecific, trypsin, chymotrypsin, lysc, lysn, gluc, argc, aspn, no_enzyme
      --terminicleaved              Allow only 'full', 'semi' or 'non' cleaved peptides
      --phospho                     Flag to pass in case of using phospho-enriched samples, changes MSGF protocol
      --maxmiscleav                 Maximum allowed number of missed cleavages for MSGF+
      --maxvarmods                  Maximum allowed number of variable modifications for a single PSM
      --minpeplen                   Minimum peptide length to search, default 7
      --maxpeplen                   Maximum peptide length to search, default 50
      --mincharge                   Minimum peptide charge search, default 2
      --maxcharge                   Maximum peptide charge search, default 6

      OUTPUT AND QUANT PARAMETERS
      --minprecursorpurity          Minimum value of precursor purity (target fraction of precursor in isolation
                                    window) for assigning isobaric quant to PSMs. PSMs below this value
                                    get assigned NA. A number between 0 and 1. Default 0 allows all PSMs
                                    quant values.
      --ms1qrttol                   Tolerance window for retention time in seconds (backwards and forwards)
                                    to align MS1 areas (from dinosaur, hardklor/kronik) to MS2 scans
      --ms1qmztol                   Tolerance window for m/z in ppm (plus and minus) to align MS1 areas found within 
                                    the RT tolerance window.
      --normalize                   Normalize isobaric values by median centering on channels of protein table
      --sampletable                 Path to sample annotation table in case of isobaric analysis
      --deqms                       Perform DEqMS differential expression analysis using sampletable
      --genes                       Produce gene table (i.e. gene names from Swissprot or ENSEMBL)
      --ensg                        Produce ENSG stable ID table (when using ENSEMBL db)
      --hirief                      File containing peptide sequences and their isoelectric points.
                                    An example can be found here:
                                    https://github.com/nf-core/test-datasets/blob/ddamsproteomics/testdata/formatted_known_peptides_ENSUniRefseq_TMT_predpi_20150825.txt
                                    For IEF fractionated samples, implies fractions in mzml definition, enables delta pI calculation
      --onlypeptides                Do not produce protein or gene level data
      --noquant                     Do not produce isobaric or MS1 quantification data
      --noms1quant                  Do not produce MS1 quantification data
      --hardklor                    Use hardklör/krönik instead of dinosaur for MS1 quant
      --keepnapsmsquant             By default the pipeline does not use PSMs with NA in any channel for isobaric 
                                    quant summarization. Use this flag and it will keep the 
                                    (potentially more noisy) PSMs in the analysis.

      REUSING PREVIOUS DATA
      --quantlookup FILE            Use previously generated SQLite lookup database containing spectra 
                                    quantification data when e.g. re-running. Need to match exactly to the
                                    mzML files of the current run
      --targetpsmlookup FILE        When adding a new sample set to existing PSM/lookup output, a complementary run,
                                    this passes the old target PSM lookup.  Any old sets with identical names to new
                                    sets will be removed prior to adding new data.
      --decoypsmlookup FILE         As for --targetpsmlookup, but filled with earlier run decoy PSMs
      --targetpsms FILE             In a complementary run, this passes the old target PSM table. If the new set has the 
                                    same name as an old set, the old set will be removed prior to adding new data.
      --decoypsms FILE              In a complementary run, this passes the old decoy PSM table.
      --ptmpsms FILE                In a complementary run, this optionally passes the old PTM PSM table, if one runs
                                    with --locptms
      --oldmzmldef                  An --mzmldef type file of a previous run you want to reuse and complement. Will be stripped of
                                    its set data for the new set that will be analyzed. Needed for basing run on another run, but only
                                    when the new run is fractionated, and only for QC purposes

    Other options:
      --outdir                      The output directory where the results will be saved
      --email                       Set this parameter to your e-mail address to get a summary e-mail with details of the run sent to you when the workflow exits
      -name                         Name for the pipeline run. If not specified, Nextflow will automatically generate a random mnemonic.

    AWSBatch options:
      --awsqueue                    The AWSBatch JobQueue that needs to be set when running on AWSBatch
      --awsregion                   The AWS Region for your AWS Batch job to run on
    """.stripIndent()
}

nextflow.enable.dsl = 1
/*
 * SET UP CONFIGURATION VARIABLES
 */

// Show help message
if (params.help){
    helpMessage()
    exit 0
}

// Configurable variables
params.name = false
params.email = false
params.plaintext_email = false

params.mzmls = false
params.mzmldef = false
params.tdb = false
params.mods = false
params.locptms = false
params.ptms = false
params.totalproteomepsms = false
// 50 is the minimum score for "good PTM" in HCD acc. to luciphor2 paper
// TODO evaluate if we need to set it higher
params.ptm_minscore_high = 50
params.phospho = false
params.maxvarmods = 2
params.isobaric = false
params.instrument = 'qe' // Default instrument is Q-Exactive
params.prectol = '10.0ppm'
params.iso_err = '-1,2'
params.frag = 'auto'
params.enzyme = 'trypsin'
params.terminicleaved = 'full' // semi, non
params.maxmiscleav = -1 // Default MSGF is no limit
params.minpeplen = 7
params.maxpeplen = 50
params.mincharge = 2
params.maxcharge = 6
params.psmconflvl = 0.01
params.pepconflvl = 0.01
params.proteinconflvl = 0.01
params.activation = 'auto' // Only for isobaric quantification
params.outdir = 'results'
params.normalize = false
params.minprecursorpurity = 0
params.ms1qmztol = 5 // in ppm
params.ms1qrttol = 18 // in seconds
params.genes = false
params.ensg = false
params.fastadelim = false
params.genefield = false
params.quantlookup = false
params.hirief = false
params.onlypeptides = false
params.noquant = false
params.noms1quant = false
params.hardklor = false
params.keepnapsmsquant = false
params.sampletable = false
params.deqms = false
params.targetpsmlookup = false
params.decoypsmlookup = false
params.targetpsms = false
params.decoypsms = false
params.ptmpsms = false
params.oldmzmldef = false

// Validate and set file inputs

// Files which are not standard can be checked here
if (params.hirief && !file(params.hirief).exists()) exit 1, "Peptide pI data file not found: ${params.hirief}"
if (params.hirief && !params.mzmldef) exit 1, "Cannot run HiRIEF delta pI calculation without fraction-annotated mzML definition file"
if (params.sampletable) {
  sampletable = file(params.sampletable)
  if( !sampletable.exists() ) exit 1, "Sampletable file not found: ${params.sampletable}"
} else {
  sampletable = 0
}

complementary_run = params.targetpsmlookup && params.decoypsmlookup && params.targetpsms && params.decoypsms
is_rerun = complementary_run && !params.mzmldef

if (complementary_run) {
  if (params.quantlookup) exit 1, "When specifying a complementary you may not pass --quantlookup"
  prev_results = Channel
    .fromPath([params.targetpsmlookup, params.decoypsmlookup, params.targetpsms, params.decoypsms, params.ptmpsms ? params.ptmpsms : 'NA'])
    .toList()
  def oldpsmheader
  new File("${params.targetpsms}").withReader { oldpsmheader = it.readLine() }
  old_fractionation = oldpsmheader.contains('Fraction')

} else if (params.targetpsmlookup || params.decoypsmlookup || params.targetpsms || params.decoypsms || params.ptmpsms) {
  exit 1, "When specifying a complementary run you need to pass all of --targetpsmlookup, --decoypsmlookup, --targetpsms, --decoypsms"

} else {
  complementary_run = false
  old_fractionation = false
  prev_results = Channel.empty()
}

output_docs = file("$baseDir/docs/output.md")

// set constant variables
accolmap = [peptides: 13, proteins: 15, ensg: 18, genes: 19]
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

availProcessors = Runtime.runtime.availableProcessors()

// parse inputs that combine to form values or are otherwise more complex.

// Isobaric input example: --isobaric 'set1:tmt10plex:127N:128N set2:tmt16plex:sweep set3:itraq8plex:intensity'
isop = params.isobaric ? params.isobaric.tokenize(' ') : false
setisobaric = isop ? isop.collect() {
  y -> y.tokenize(':')
}.collectEntries() {
  x-> [x[0], x[1].replaceAll('tmtpro', 'tmt16plex')]
} : false
setdenoms = isop ? isop.collect() {
  y -> y.tokenize(':')
}.collectEntries() {
  x-> [x[0], x[2..-1]]
} : false

normalize = (!params.noquant && (params.normalize || params.deqms) && params.isobaric)

// AWSBatch sanity checking
if(workflow.profile == 'awsbatch'){
    if (!params.awsqueue || !params.awsregion) exit 1, "Specify correct --awsqueue and --awsregion parameters on AWSBatch!"
    if (!workflow.workDir.startsWith('s3') || !params.outdir.startsWith('s3')) exit 1, "Specify S3 URLs for workDir and outdir parameters on AWSBatch!"
}


// Has the run name been specified by the user?
//  this has the bonus effect of catching both -name and --name
custom_runName = params.name
// if not wf.runName (-name or auto) is like "crazy_euler" or other "{adjective}_{scientist}"
if( !(workflow.runName ==~ /[a-z]+_[a-z]+/) ){
  custom_runName = workflow.runName
}

// Check workDir/outdir paths to be S3 buckets if running on AWSBatch
// related: https://github.com/nextflow-io/nextflow/issues/813
if( workflow.profile == 'awsbatch') {
    if(!workflow.workDir.startsWith('s3:') || !params.outdir.startsWith('s3:')) exit 1, "Workdir or Outdir not on S3 - specify S3 Buckets for each to run on AWSBatch!"
}


// Header log info
log.info """=======================================================
                                          ,--./,-.
          ___     __   __   __   ___     /,-._.--~\'
    |\\ | |__  __ /  ` /  \\ |__) |__         }  {
    | \\| |       \\__, \\__/ |  \\ |___     \\`-._,-`-,
                                          `._,._,\'

lehtiolab/ddamsproteomics v${workflow.manifest.version}"
======================================================="""
def summary = [:]
summary['Pipeline Name']  = 'lehtiolab/ddamsproteomics'
summary['Pipeline Version'] = workflow.manifest.version
summary['Run Name']     = custom_runName ?: workflow.runName
summary['mzMLs']        = multifile_format(params.mzmls)
summary['or mzML definition file']        = params.mzmldef
summary['Target DB'] = multifile_format(params.tdb)
summary['Sample annotations'] = params.sampletable
summary['Modifications'] = params.mods
summary['Labile PTMs'] = params.locptms
summary['PTMs'] = params.ptms
summary['Minimum Luciphor2 score of PTM'] = params.ptm_minscore_high
summary['Phospho enriched samples'] = params.phospho
summary['Total proteome normalization PSM table'] = params.totalproteomepsms
summary['Instrument'] = params.mzmldef ? 'Set per mzML file in mzml definition file' : params.instrument
summary['Precursor tolerance'] = params.prectol
summary['Isotope error'] = params.iso_err
summary['Fragmentation method'] = params.frag
summary['Enzyme'] = params.enzyme
summary['Allowed peptide termini cleavage'] = params.terminicleaved
summary['Allowed number of missed cleavages'] = params.maxmiscleav
summary['Maximum number of variable modifications per PSM'] = params.maxvarmods
summary['Minimum peptide length'] = params.minpeplen
summary['Maximum peptide length'] = params.maxpeplen
summary['Minimum peptide charge'] = params.mincharge
summary['Maximum peptide charge'] = params.maxcharge
summary['PSM FDR cutoff'] = params.psmconflvl
summary['Peptide FDR cutoff'] = params.pepconflvl
summary['Isobaric label sets/denominators'] = params.isobaric
summary['Isobaric quant: activation method'] = params.activation
summary['Isobaric PSM minimum precursor purity'] = params.minprecursorpurity
summary['Retention time tolerance for MS1 alignment'] = params.ms1qrttol
summary['m/z tolerance for MS1 alignment'] = params.ms1qmztol
summary['Keep PSMs for quant with NA in any channel'] = params.keepnapsmsquant
summary['Explicit isobaric normalization'] = params.normalize
summary['Perform DE analysis (implies normalization)'] = params.deqms
summary['Output genes'] = params.genes
summary['Output ENSG IDs'] = params.ensg
summary['Custom FASTA delimiter'] = params.fastadelim 
summary['Custom FASTA gene field'] = params.genefield
summary['Premade quant data SQLite'] = params.quantlookup
summary['Previous run target results SQLite'] = params.targetpsmlookup
summary['Previous run decoy results SQLite'] = params.decoypsmlookup
summary['Previous run target PSMs'] = params.targetpsms
summary['Previous run decoy PSMs'] = params.decoypsms
summary['Previous run PTM PSMs'] = params.ptmpsms
summary['Previous run mzml definition'] = params.oldmzmldef
summary['HiRIEF pI peptide data'] = params.hirief 
summary['Only output peptides'] = params.onlypeptides
summary['Do not quantify'] = params.noquant
summary['Do not quantify MS1'] = params.noms1quant
summary['Use hardklor instead dinosaur for MS1 quant'] = params.hardklor
summary['Max Memory']   = params.max_memory
summary['Max CPUs']     = params.max_cpus
summary['Max Time']     = params.max_time
summary['Output dir']   = params.outdir
summary['Working dir']  = workflow.workDir
summary['Container Engine'] = workflow.containerEngine
if(workflow.containerEngine) summary['Container'] = workflow.container
summary['Current home']   = "$HOME"
summary['Current user']   = "$USER"
summary['Current path']   = "$PWD"
summary['Working dir']    = workflow.workDir
summary['Output dir']     = params.outdir
summary['Script dir']     = workflow.projectDir
summary['Config Profile'] = workflow.profile
if(workflow.profile == 'awsbatch'){
   summary['AWS Region'] = params.awsregion
   summary['AWS Queue'] = params.awsqueue
}
if(params.email) summary['E-mail Address'] = params.email
log.info summary.collect { k,v -> "${k.padRight(15)}: $v" }.join("\n")
log.info "========================================="


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
    msgf_plus | head -n2 | grep Release > v_msgf.txt
    ${!noms1 && !params.hardklor ? 'dinosaur | head -n2 | grep Dinosaur > v_dino.txt || true' : ''}
    ${!noms1 && params.hardklor ? 'hardklor | head -n1 > v_hk.txt || true' : ''}
    kronik | head -n2 | tr -cd '[:digit:],\\.' > v_kr.txt || true
    #luciphor2 |& grep Version > v_luci.txt # incorrect version from binary (2014), echo below
    echo Version: 2020_04_03 > v_luci.txt # deprecate when binary is correct
    percolator -h |& head -n1 > v_perco.txt || true
    msstitch --version > v_mss.txt
    echo 2.9.1 > v_openms.txt
    Rscript <(echo "packageVersion('DEqMS')") > v_deqms.txt
    scrape_software_versions.py > software_versions.yaml
    """
}

if (workflow.profile.tokenize(',').intersect(['test', 'test_nofrac'])) { 
  // Profile 'test' delivers mzmlPaths
  Channel
    .from(params.mzmlPaths)
    .set { mzml_in }
}
else if (!params.mzmldef && params.mzmls) {
  Channel
    .fromPath(params.mzmls)
    .map { it -> [it, params.instrument, 'NA'] }
    .set { mzml_in }
} else if (is_rerun) {
  fractionation_in = false
  Channel
    .empty()
    .set { mzml_in }
} else {
  mzmllines = file("${params.mzmldef}").readLines().collect { it.tokenize('\t') }
  Channel
    .from(mzmllines)
    .set { mzml_in }
  fractionation_in = mzmllines.any { it.size() == 5}
  if (fractionation_in) {
    println("Fractions detected in mzml definition, will run as fractionated")
    if (mzmllines.any { it.size() == 4}) {
      println("Forcing files without fraction to file-named fraction in their specified plate")
    } else if (mzmllines.any { it.size() == 3}) {
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


process createTargetDecoyFasta {
 
  input:
  path(tdb) from Channel.fromPath(params.tdb).toList()

  output:
  file('db.fa') into concatdb
  set file('tdb'), file("decoy.fa") into bothdbs

  script:
  """
  ${tdb.size() > 1 ? "cat ${tdb.collect() { "\"${it}\"" }.join(' ')} > tdb" : "mv '$tdb' tdb"}
  check_fasta.py tdb
  msstitch makedecoy -i tdb -o decoy.fa --scramble tryp_rev --ignore-target-hits
  cat tdb decoy.fa > db.fa
  """
}

bothdbs.into { psmdbs; fdrdbs; ptmdbs }


def fr_or_file(it, length) {
  // returns either fraction number or file from line
  // file name is used in QC plots where no frac is available and fraction plot is enforced,
  // e.g. when mixing fractions and non-fractions
  return it.size() > length ? it[length] : "${file(it[0]).baseName}.${file(it[0]).extension}"
}

def plate_or_no(it, length) {
  return it.size() > 3 ? it[3] : "no_plate"
}



regex_specialchars = '[&<>\'"]'
def stripchars_infile(x, return_oldfile=false) {
  // FIXME %, ?, * fn turns a basename into a list because they are wildcards
  // Replace special characters since they cause trouble in percolator XML output downstream
  // e.g. & is not allowed in XML apparently (percolator doesnt encode it)
  // and LXML then crashes on reading it.
  // Also NF doesnt quote e.g. semicolons it seems.
  def scriptinfile = "${x.baseName}.${x.extension}"
  def parsed_file = scriptinfile.replaceAll(regex_specialchars, '_')
  if (return_oldfile) {
    return [parsed_file != scriptinfile, parsed_file, scriptinfile]
  } else {
    return [parsed_file != scriptinfile, parsed_file]
  }
}

// Parse mzML input to get files and sample names etc
// get setname, sample name (baseName), input mzML file. 
// Set platename to setname if not specified. 
// Set fraction name to NA if not specified
mzml_in
  .tap { mzmlfiles_counter; mzmlfiles_qlup_sets } // for counting-> timelimits; getting sets from supplied lookup
  .map { it -> [it[2].replaceAll('[ ]+$', '').replaceAll('^[ ]+', ''), file(it[0]).baseName.replaceAll(regex_specialchars, '_'), file(it[0]), it[1], plate_or_no(it, 3), fr_or_file(it, 4)] }
.view()
  .tap { mzmlfiles; mzml_luciphor; pre_isoquant; ms1quant; sample_mzmlfn}
  .combine(concatdb)
  .set { mzml_msgf }

/*
* Step 1: Extract quant data from peptide spectra
*/


process centroidMS1 {
  container 'chambm/pwiz-skyline-i-agree-to-the-vendor-licenses:3.0.20066-729ef9c41'
  when: !params.quantlookup && !params.noquant && setisobaric && setisobaric[setname]

  input:
  set val(setname), val(sample), path(infile), val(instr), val(platename), val(fraction) from pre_isoquant

  output:
  set val(setname), val(sample), path(parsed_infile), val(instr), val(platename), val(fraction) into isoquant

  script:
  (is_stripped, parsed_infile) = stripchars_infile(infile)
  """
  ${is_stripped ? "ln -s ${infile} ${parsed_infile}" : ''}
  wine msconvert ${parsed_infile} --outfile centroidms1.mzML --filter 'peakPicking true 1' ${instr == 'timstof' ? "--filter sortByScanTime" : ''}
  mv centroidms1.mzML ${parsed_infile}
  """
}


process quantifyMS1 {
  when: !params.quantlookup && !params.noquant && !params.noms1quant

  input:
  set val(setname), val(sample), path(infile), val(instr), val(platename), val(fraction) from ms1quant
  file(hkconf) from Channel.fromPath("$baseDir/assets/hardklor.conf").first()

  output:
  set val(sample), path("${sample}.features.tsv") optional true into dino_out 
  set val(sample), path("${sample}.kr") optional true into kronik_out 

  script:
  (is_stripped, parsed_infile) = stripchars_infile(infile)
  """
  ${is_stripped ? "ln -s ${infile} ${parsed_infile}" : ''}
  # Dinosaur is first choice for MS1 quant
  ${!params.hardklor ? "dinosaur --concurrency=${task.cpus * params.threadspercore} ${parsed_infile}" : ''}
  # Hardklor/Kronik can be used as a backup, using --hardklor
  ${params.hardklor ? "hardklor <(cat $hkconf <(echo $parsed_infile hardklor.out)) && kronik -c 5 -d 3 -g 1 -m 8000 -n 600 -p 10 hardklor.out ${sample}.kr" : ''}
  """
}


process isoquantSpectra {
  container "${workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
    'https://depot.galaxyproject.org/singularity/openms:2.9.1--h135471a_1' :
    'quay.io/biocontainers/openms:2.9.1--h135471a_1'}"

  when: !params.quantlookup && !params.noquant

  input:
  set val(setname), val(sample), file(infile), val(instr), val(platename), val(fraction) from isoquant

  output:
  set val(sample), file(outfile) into isobaricxml

  script:
  outfile = "${infile.baseName}.consensusXML"
  activationtype = [auto: 'auto', any: 'any', hcd:'beam-type collision-induced dissociation', cid:'Collision-induced dissociation', etd:'Electron transfer dissociation'][params.activation]
  isobtype = setisobaric && setisobaric[setname] ? setisobaric[setname] : false
  plextype = isobtype ? isobtype.replaceFirst(/[0-9]+plex/, "") : 'false'
  massshift = [tmt:0.0013, itraq:0.00125, false:0][plextype]

  """
  ${isobtype ? "IsobaricAnalyzer -type $isobtype -in \"${infile}\" -out \"${infile.baseName}.consensusXML\" -extraction:select_activation \"$activationtype\" -extraction:reporter_mass_shift $massshift -extraction:min_precursor_intensity 1.0 -extraction:keep_unannotated_precursor true -quantification:isotope_correction true" : ''}
  """
}


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
  oldmzmls_fn = Channel.value(false)
  oldmzml_sets = Channel.value([])
}

// Prepare mzml files (sort, collect) for processes that need all of them
mzmlfiles
  .toList()
  .map { it.sort( {a, b -> a[1] <=> b[1]}) } // sort on sample for consistent .sh script in -resume
  .map { it -> [it.collect() { it[0] }, it.collect() { it[2] }, it.collect() { it[4] }, it.collect() { it[1] }, it.collect() { it[5] } ] } // lists: [sets], [mzmlfiles], [plates], [basenames], [fractions]
  .into { mzmlfiles_psmqc; mzmlfiles_rest }

mzmlfiles_rest
  .map { it -> it[0..2] } // remove basenames, fractions
  .into { mzmlfiles_all; mzmlfiles_all_count }

if (!is_rerun) {
  mzmlfiles_counter
    .count()
    .subscribe { println "$it mzML files in analysis" }
    .into { mzmlcount_psm; mzmlcount_percolator }
} else {
  Channel.value(0).into { mzmlcount_psm; mzmlcount_percolator }
}



process complementSpectraLookupCleanPSMs {

  when: complementary_run

  input:
  set val(in_setnames), path(mzmlfiles), val(platenames) from mzmlfiles_all
  set path(tlup), path(dlup), path(tpsms), path(dpsms), path(ptmpsms) from prev_results

  output:
  set path('t_cleaned_psms.txt'), path('d_cleaned_psms.txt') into cleaned_psms
  set path('target_db.sqlite'), path('decoy_db.sqlite') into complemented_speclookup 
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
      msstitch deletesets -i ${tpsms} -o t_cleaned_psms.txt --dbfile target_db.sqlite --setnames "${setnames.collect() { "'${it}'" }.join(' ')}"
      msstitch deletesets -i ${dpsms} -o d_cleaned_psms.txt --dbfile decoy_db.sqlite --setnames "${setnames.collect() { "'${it}'" }.join(' ')}"
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

  when: !params.quantlookup && !complementary_run

  input:
  set val(setnames), file(mzmlfiles), val(platenames) from mzmlfiles_all

  output:
  set path('target_db.sqlite'), path('decoy_db.sqlite') into newspeclookup 
  path 'target_db.sqlite' into ptm_lookup_new

  script:
  """
  ${mzmlfiles.collect() { stripchars_infile(it, return_oldfile=true) }.findAll{ it[0] }.collect() { "ln -s '${it[2]}' '${it[1]}'" }.join(' && ')}

  msstitch storespectra --spectra ${mzmlfiles.collect() { stripchars_infile(it)[1] }.join(' ')} --setnames ${setnames.collect() { "'$it'" }.join(' ')} -o target_db.sqlite
  ln -s target_db.sqlite decoy_db.sqlite
  """
}

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
  mzmlfiles_qlup_sets
    .map { it -> it[2] } 
    .unique()
    .toList()
    .set { allsetnames }

  // if not using this youll have a combine on an open channel without
  // anything from complement cleaner. Will not run createPTMLookup then
  cleaned_ptmpsms = Channel.value('NA')
}

// Collect all MS1 dinosaur/kronik output for quant lookup building process
sample_mzmlfn
  .map { [ it[1], stripchars_infile(it[2])[1] ] } // extract samplename, mzML
  .join(dino_out.concat(kronik_out), remainder: true)
  .join(isobaricxml, remainder: true)
  .toList()
  .map { it.sort({a, b -> a[0] <=> b[0]}) }
  .transpose()
  .toList()
  .set { quantfiles_sets }


// Need to populate channels depending on if a pre-made quant lookup has been passed
// even if not needing quant (--noquant) this is necessary or NF will error
newspeclookup
  .concat(complemented_speclookup)
  .tap { prespectoquant }
  .map { it -> it[0] } // get only target lookup
  .set { countlookup }

if (params.noquant && !params.quantlookup) {
  // Noquant, fresh spectra lookup scenario
  prespectoquant
    .flatMap { it -> [[it[0], 'target'], [it[1], 'decoy']] }
    .set { specquant_lookups }
   spectoquant = Channel.empty()

} else if (params.quantlookup) {
  // Runs with a premade quant lookup eg from previous search
  spectoquant = Channel.empty()
  Channel
    .fromPath(params.quantlookup)
    .flatMap { it -> [[it, 'target'], [it, 'decoy']] }
    .tap { specquant_lookups }
    .filter { it[1] == 'target' }
    .map { it -> it[0] }
    .into { ptm_lookup_in; countlookup }
} else if (is_rerun) {
  // In case of rerun with same sets, no new search but only some different post-identification
  // output options
  spectoquant = Channel.empty()
  prespectoquant  // contains  t, d lookups
    .flatMap { it -> [[it[0], 'target'], [it[1], 'decoy']] }
    .set{ specquant_lookups }
} else {
  prespectoquant.set { spectoquant }
}

if (!params.quantlookup) {
  ptm_lookup_old
    .concat(ptm_lookup_new)
    .set { ptm_lookup_in }
}


// Set names are first item in input lists, collect them for PSM tables and QC purposes
allsetnames 
  .into { setnames_featqc; setnames_psms; setnames_psmqc; setnames_ptms }


process quantLookup {

  publishDir "${params.outdir}", mode: 'copy', overwrite: true, saveAs: {it == 'target.sqlite' ? 'quant_lookup.sql' : null }

  when: !params.quantlookup && !params.noquant

  input:
  set path(tlookup), path(dlookup) from spectoquant
  set val(samples), val(mzmlnames), file(ms1fns), file(isofns) from quantfiles_sets

  output:
  set path('target.sqlite'), path(dlookup) into newquantlookup

  script:
  """
  # SQLite lookup needs copying to not modify the input file which would mess up a rerun with -resume
  cat $tlookup > target.sqlite
  msstitch storequant --dbfile target.sqlite --spectra ${mzmlnames.collect() { "'${it}'" }.join(' ')}  \
    ${!params.noms1quant ? "--mztol ${params.ms1qmztol} --mztoltype ppm --rttol ${params.ms1qrttol} ${params.hardklor ? "--kronik ${ms1fns.collect() { "$it" }.join(' ')}" : "--dinosaur ${ms1fns.collect() { "$it" }.join(' ')}"}" : ''} \
    ${params.isobaric ? "--isobaric ${isofns.collect() { "$it" }.join(' ')}" : ''}
  """
}


if (!params.quantlookup && !params.noquant && !is_rerun) {
  newquantlookup
    .flatMap { it -> [[it[0], 'target'], [it[1], 'decoy']] }
    .set { specquant_lookups }
} 

mzmlfiles_all_count
  .merge(countlookup)
  .set { specfilein }


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


fractionation = fractionation_in || old_fractionation
if (fractionation && complementary_run) { 
  if (!params.oldmzmldef || !file(params.oldmzmldef).exists()) exit 1, 'Fractionation with complementing run needs an --oldmzmldef file'
  specfilems2
    .set { scans_platecount }
} else if (fractionation) {
  specfilems2
    .set { scans_platecount }
} else {
  specfilems2
    // scans_platecount is not used, format directly for scans_result
    .map { it -> [it[3], 'NA', ['noplates']] }
    .into { scans_platecount; scans_result }
}


process countMS2sPerPlate {

  publishDir "${params.outdir}", mode: 'copy', overwrite: true 
  when: fractionation

  input:
  set val(setnames), file(mzmlfiles), val(platenames), file('nr_spec_per_file') from scans_platecount
  file oldmzmls_fn

  output:
  set path('nr_spec_per_file'), path('scans_per_plate') into scans_perplate
  file('allplates') into allplates

  script:
  splates = [setnames, platenames].transpose().collect() { "${it[0]}_${it[1]}" }
  """
  #!/usr/bin/env python
  import os
  import re
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
  with open('nr_spec_per_file') as fp:
      for line in fp:
          setname, fn, scans = line.strip('\\n').split('|')
          setplate = '{}_{}'.format(setname, fileplates[fn][setname])
          platescans[setplate] += int(scans)
  with open('scans_per_plate', 'w') as fp:
      for plate, scans in platescans.items():
          fp.write('{}\\t{}\\n'.format(plate, scans))
  """
}

if (fractionation) {
  allplates
    .splitText()
    .map { it -> it.trim() }
    .toList()
    .map { it -> [it] }
    .set { allplates_split }
  scans_perplate
    .combine(allplates_split)
    .set { scans_result }
}

/*
* Step 2: Identify peptides
*/

process msgfPlus {
  cpus = availProcessors < 4 ? availProcessors : 4

  input:
  set val(setname), val(sample), path(infile), val(instrument), val(platename), val(fraction), path(db) from mzml_msgf

  output:
  set val(setname), val(sample), path("${sample}.mzid"), path("${sample}.mzid.tsv") into mzids
  
  script:
  isobtype = setisobaric && setisobaric[setname] ? setisobaric[setname] : false
  isobtype_parsed = ['tmt16plex', 'tmt18plex'].any { it == isobtype } ? 'tmtpro' : isobtype
  // protcol 0 is automatic, msgf checks in mod file, TMT/phospho should be run with 1
  // see at https://github.com/MSGFPlus/msgfplus/issues/19
  msgfprotocol = params.phospho ? setisobaric[setname][0..4] == 'itraq' ? 3 : 1 : 0
  msgfinstrument = [lowres:0, velos:1, qe:3, qehf: 3, false:0, qehfx:1, lumos:1, timstof:2][instrument]
  fragmeth = [auto:0, cid:1, etd:2, hcd:3, uvpd:4][params.frag]
  enzyme = params.enzyme.indexOf('-') > -1 ? params.enzyme.replaceAll('-', '') : params.enzyme
  enzyme = [unspecific:0, trypsin:1, chymotrypsin: 2, lysc: 3, lysn: 4, gluc: 5, argc: 6, aspn:7, no_enzyme:9][enzyme]
  ntt = [full: 2, semi: 1, non: 0][params.terminicleaved]

  // the string in "scriptinfile" does not have NF escaping characters like & (e.g. in FAIMS 35&65),
  // which NF does to "infile". That would work fine but not if the files are quoted in the 
  // script, then they cant be found when there is \&.
  // Replace those characters anyway since they cause trouble in percolator XML output downstream
  (is_stripped, parsed_infile) = stripchars_infile(infile)
  """
  ${is_stripped ? "ln -s ${infile} '${parsed_infile}'" : ''}
  create_modfile.py ${params.maxvarmods} "${params.msgfmods}" "${params.mods}${isobtype ? ";${isobtype_parsed}" : ''}${params.ptms ? ";${params.ptms}" : ''}${params.locptms ? ";${params.locptms}" : ''}"
  
  msgf_plus -Xmx${task.memory.toMega()}M -d $db -s '$parsed_infile' -o "${sample}.mzid" -thread ${task.cpus * params.threadspercore} -mod "mods.txt" -tda 0 -maxMissedCleavages $params.maxmiscleav -t ${params.prectol}  -ti ${params.iso_err} -m ${fragmeth} -inst ${msgfinstrument} -e ${enzyme} -protocol ${msgfprotocol} -ntt ${ntt} -minLength ${params.minpeplen} -maxLength ${params.maxpeplen} -minCharge ${params.mincharge} -maxCharge ${params.maxcharge} -n 1 -addFeatures 1
  msgf_plus -Xmx3500M edu.ucsd.msjava.ui.MzIDToTsv -i "${sample}.mzid" -o out.tsv
  awk -F \$'\\t' '{OFS=FS ; print \$0, "Biological set" ${fractionation ? ', "Strip", "Fraction"' : ''}}' <( head -n+1 out.tsv) > "${sample}.mzid.tsv"
  awk -F \$'\\t' '{OFS=FS ; print \$0, "$setname" ${fractionation ? ", \"$platename\", \"$fraction\"" : ''}}' <( tail -n+2 out.tsv) >> "${sample}.mzid.tsv"
  rm ${db.baseName.replaceFirst(/\.fasta/, "")}.c*
  """
}


mzids
  .groupTuple()
  .set { mzids_2pin }


process percolator {

  input:
  set val(setname), val(samples), file(mzids), file(tsvs) from mzids_2pin
  val(mzmlcount) from mzmlcount_percolator

  output:
  set path("${setname}_target.tsv"), val('target') into tmzidtsv_perco optional true
  set path("${setname}_decoy.tsv"), val('decoy') into dmzidtsv_perco optional true
  set val(setname), file('allpsms') optional true into unfiltered_psms
  file('warnings') optional true into percowarnings

  script:
  """
  ${mzids.collect() { "echo '$it' >> metafile" }.join('&&')}
  msgf2pin -o percoin.tsv -e ${params.enzyme} -P "decoy_" metafile
  percolator -j percoin.tsv -X perco.xml -N 500000 --decoy-xml-output -Y
  mkdir outtables
  msstitch perco2psm --perco perco.xml -d outtables -i ${listify(tsvs).collect(){ "${it}"}.join(' ')} --mzids ${listify(mzids).collect(){ "${it}"}.join(' ')} ${!params.locptms ? "--filtpsm ${params.psmconflvl} --filtpep ${params.pepconflvl}" : ''}
  msstitch concat -i outtables/* -o allpsms
  ${params.locptms ? 
    "msstitch conffilt -i allpsms -o filtpsm --confcolpattern 'PSM q-value' --confidence-lvl ${params.psmconflvl} --confidence-better lower && \
    msstitch conffilt -i filtpsm -o psms --confcolpattern 'peptide q-value' --confidence-lvl ${params.pepconflvl} --confidence-better lower" : 'mv allpsms psms'}
  msstitch split -i psms --splitcol \$(head -n1 psms | tr '\t' '\n' | grep -n ^TD\$ | cut -f 1 -d':')
  ${['target', 'decoy'].collect() { 
    "test -f '${it}.tsv' && mv '${it}.tsv' '${setname}_${it}.tsv' || echo 'No ${it} PSMs found for set ${setname} at PSM FDR ${params.psmconflvl} and peptide FDR ${params.pepconflvl} ${it == 'decoy' ? '(not possible to output protein/gene FDR)' : ''}' >> warnings" }.join(' && ') }
  """
}

totalproteomepsms = params.totalproteomepsms ? Channel.fromPath(params.totalproteomepsms) : Channel.from('NA')

dmzidtsv_perco
  .ifEmpty([file('NO__FILE'), 'decoy'])
  .set { d_perco_checked }

// Collect percolator data of target/decoy and feed into PSM table creation
tmzidtsv_perco
  .ifEmpty([file('NO__FILE'), 'target'])
  .concat(d_perco_checked)
  .groupTuple(by: 1) // group by TD
  .join(specquant_lookups, by: 1) // join on TD
  .combine(psmdbs)
  .combine(totalproteomepsms)
  .set { psmswithout_oldpsms }
if (complementary_run) {
  psmswithout_oldpsms
    .join(td_oldpsms, remainder: true)
    .set { prepsm }
} else {
  psmswithout_oldpsms
    .map { it -> it[0..-1] + 'NA' }
    .set { prepsm }
}

/*
* Step 3: Post-process peptide identification data
*/

hiriefpep = params.hirief ? Channel.fromPath([params.hirief, params.hirief]) : Channel.value('NA')

process createPSMTable {

  publishDir "${params.outdir}", mode: 'copy', overwrite: true, saveAs: {["target_psmlookup.sql", "decoy_psmlookup.sql", "target_psmtable.txt", "decoy_psmtable.txt"].contains(it) && (is_rerun || !no_psms) ? it : null}

  input:
  set val(td), path(psms), path('lookup'), path(tdb), path(ddb), file('tppsms.txt'), file(cleaned_oldpsms) from prepsm
  file(trainingpep) from hiriefpep
  val(mzmlcount) from mzmlcount_psm  // For time limits
  val(setnames) from setnames_psms

  output:
  set val(td), file("${outpsms}") into psm_result
  set val(td), file({setnames.collect() { "${it}.tsv" }}) optional true into setpsmtables
  file({setnames.collect() { "tppsms/${it}.tsv" }}) optional true into totalprotpsms_allsets
  set val(td), file("${psmlookup}") into psmlookup
  file('warnings') optional true into psmwarnings

  script:
  psmlookup = "${td}_psmlookup.sql"
  outpsms = "${td}_psmtable.txt"
  // NO__FILE is when target percolator has no output in any set, exit immediately if not rerun
  no_psms = psms[0].name == 'NO__FILE'
  no_target = td == 'target' && no_psms
  no_decoy = td == 'decoy' && no_psms
  quant = !params.noquant && td == 'target'
  """
  ${no_target && !is_rerun ? "echo 'No target PSMs made the combined PSM / peptide FDR cutoff (${params.psmconflvl} / ${params.pepconflvl})' && exit 1" : ''}
  ${no_decoy && !is_rerun ? "echo 'No decoy PSMs in any set at FDR cutoff, will not produce protein/gene tables' > warnings && touch ${outpsms} && touch ${psmlookup} && exit 0" : ''}
  ${!is_rerun ? "msstitch concat -i ${listify(psms).collect() {"$it"}.join(' ')} -o psms.txt" : ''}
  # SQLite lookup needs copying to not modify the input file which would mess up a rerun with -resume
  cat lookup > $psmlookup
  ${!is_rerun ? "sed '0,/\\#SpecFile/s//SpectraFile/' -i psms.txt": ''}
  ${is_rerun ? "mv ${cleaned_oldpsms} psmsrefined" : "msstitch psmtable -i psms.txt --dbfile $psmlookup --addmiscleav -o psmsrefined --spectracol 1 \
    ${params.onlypeptides ? '' : "--fasta \"${td == 'target' ? "${tdb}" : "${ddb}"}\" --genes"} \
    ${quant ? "${!params.noms1quant ? '--ms1quant' : ''} ${params.isobaric ? "--isobaric --min-precursor-purity ${params.minprecursorpurity}" : ''}" : ''} \
    ${!params.onlypeptides ? '--proteingroup' : ''} \
    ${complementary_run ? "--oldpsms ${cleaned_oldpsms}" : ''}" }
  sed 's/\\#SpecFile/SpectraFile/' -i psmsrefined
  ${params.hirief && td == 'target' && !is_rerun ? "echo \'${groovy.json.JsonOutput.toJson(params.strips)}\' >> strip.json && peptide_pi_annotator.py -i $trainingpep -p psmsrefined --out $outpsms --stripcolpattern Strip --pepcolpattern Peptide --fraccolpattern Fraction --stripdef strip.json --ignoremods \'*\'": "mv psmsrefined ${outpsms}"}
  msstitch split -i ${outpsms} --splitcol bioset
  # In decoy PSM table process, also split the target total proteome normalizer table if necessary.
  # Doing it in decoy saves time, since target is usally largest table and slower
  ${params.totalproteomepsms && td == 'decoy' ? "mkdir tppsms && msstitch split -i tppsms.txt -d tppsms --splitcol bioset" : ''}
  """
}

def listify(it) {
  /* This function is useful when needing a list even when having a single item
  - Single items in channels get unpacked from a list
  - Processes expect lists. Even though it would be fine
  without a list, for single-item-lists any special characters are not escaped by NF
  in the script, which leads to errors. See:
  https://github.com/nextflow-io/nextflow/discussions/4240
  */
  return it instanceof java.util.List ? it : [it]
}

// Collect setnames and merge with PSM tables for peptide table creation
setpsmtables
  .map { it -> [it[0], listify(it[1])] }
  .map{ it -> [it[0], it[1].collect() { it.baseName.replaceFirst(/\.tsv$/, "") }, it[1]]}
  .transpose()
  .tap { psm_pep }
  .filter { it -> it[0] == 'target' }
  .combine(ptmdbs)
  .map { it -> it[1..3] }
  .into { psm_ptm; stabileptm_prep }

mzml_luciphor
  .map { it -> [it[0], it[2]] } // only need setname and mzml
  .groupTuple()
  .join(psm_ptm)
  .join(unfiltered_psms)
  .set { psm_luciphor }


process luciphorPTMLocalizationScoring {

  cpus = availProcessors < 4 ? availProcessors : 4
  when: params.locptms

  input:
  set val(setname), path(mzmls), path('psms'), path(tdb), path(allpsms) from psm_luciphor

  output:
  set val(setname), path('labileptms.txt') into luciphor_all
  path 'warnings' optional true into luciphor_warnings

  script:
  denom = !params.noquant && setdenoms ? setdenoms[setname] : false
  specialdenom = denom && (denom[0] == 'sweep' || denom[0] == 'intensity')
  isobtype = setisobaric && setisobaric[setname] ? "${setisobaric[setname]}" : ''
  isobtype_parsed = ['tmt16plex', 'tmt18plex'].any { it == isobtype } ? 'tmtpro' : isobtype
  mods = params.mods ? params.mods.tokenize(';').join(' ') : ''
  stab_ptms = params.ptms ? params.ptms.tokenize(';').join(' ') : ''
  lab_ptms = params.locptms.tokenize(';').join(' ')

  """
  ${mzmls.collect() { stripchars_infile(it, return_oldfile=true) }.findAll{ it[0] }.collect() { "ln -s '${it[2]}' '${it[1]}'" }.join(' && ')}
  # Split allpsms to get target PSMs
  sed '0,/\\#SpecFile/s//SpectraFile/' -i "${allpsms}"
  msstitch split -i "${allpsms}" --splitcol \$(head -n1 "${allpsms}" | tr '\t' '\n' | grep -n ^TD\$ | cut -f 1 -d':')
  export MZML_PATH=\$(pwd)
  export MINPSMS=${params.minpsms_luciphor}
  export ALGO=${['hcd', 'auto', 'any'].any { it == params.activation } ? '1' : '0'}
  export MAXPEPLEN=${params.maxpeplen}
  export MAXCHARGE=${params.maxcharge}
  export THREAD=${task.cpus * params.threadspercore}
  export MS2TOLVALUE=0.025
  export MS2TOLTYPE=Da
  cat "$baseDir/assets/luciphor2_input_template.txt" | envsubst > lucinput.txt
  luciphor_prep.py --psmfile target.tsv --template lucinput.txt --modfile "${params.msgfmods}" \
      --labileptms "${params.locptms}" --mods ${mods} ${isobtype_parsed} ${stab_ptms} \
      -o luciphor.out --lucipsms lucipsms
  luciphor2 -Xmx${task.memory.toGiga()}G luciphor_config.txt 2>&1 | grep 'not have enough PSMs' && echo 'Not enough PSMs for luciphor FLR calculation in set ${setname}' > warnings
  luciphor_parse.py --minscore ${params.ptm_minscore_high} -o labileptms.txt \
     --luci_in luciphor.out --luci_scores all_scores.debug --psms psms \
     --modfile "${params.msgfmods}" --labileptms ${lab_ptms} \
     ${params.ptms ? "--stabileptms ${stab_ptms}": ''} --mods ${mods} ${isobtype_parsed} \
     --fasta "${tdb}"
  """
}
// FIXME msgfmods is false? oxidation so probably never.


process stabilePTMPrep {

  when: params.ptms
  
  input:
  set val(setname), path('psms'), path(tdb) from stabileptm_prep
  
  output:
  path('stabileptms.txt') into stabileptms
  
  script:
  stab_ptms = params.ptms.tokenize(';').join(' ')
  mods = params.mods ? params.mods.tokenize(';').join(' ') : ''
  lab_ptms = params.locptms ? params.locptms.tokenize(';').join(' ') : ''
  isobtype = setisobaric && setisobaric[setname] ? setisobaric[setname] : ''
  isobtype_parsed = ['tmt16plex', 'tmt18plex'].any { it == isobtype } ? 'tmtpro' : isobtype
  """
  nonlabile_ptm_columns.py --psms psms -o stabileptms.txt --modfile "${params.msgfmods}" --fasta "${tdb}" \
      --stabileptms $stab_ptms \
      ${params.locptms ? "--labileptms $lab_ptms" : ""} \
      ${params.mods ? "--mods ${isobtype_parsed} ${mods}" : ''}
  """
}


// Sort to be able to resume
if (params.locptms && !is_rerun) {
  luciphor_all
    .toList()
    .map { it.sort( {a, b -> a[0] <=> b[0]}) } // sort on setname
    .transpose()
    .toList()
   // add all setnames for cases of rerun complement where only limited amount of sets is run
    .combine(setnames_ptms)
    .map { it -> [(it[0] + it[2..-1]).unique(), it[1]] }
    .set { ptm_allsets }
} else {
  setnames_ptms 
    .map { it -> [it, ''] }
    .set { ptm_allsets }
}

stabileptms
  .ifEmpty('')
  .toList()
  .set { allstabileptms }
    
ptm_allsets
  .combine(allstabileptms)
  .combine(cleaned_ptmpsms)
  .set { lucptmfiles_to_lup }


ptms_used = params.locptms || params.ptms || params.ptmpsms

// FIXME currently not storing the PTMs in SQL or similar, would be good to keep this state
// between runs for complement/rerun etc ?
if (params.ptmpsms && !(params.locptms || params.ptms)) exit 1, "In a rerun with --ptmpsms you  must specify which PTMs you have used with --locptms or --ptms"

process createPTMLookup {
  publishDir "${params.outdir}", mode: 'copy', overwrite: true, saveAs: {it == ptmtable ? ptmtable: null}

  when: ptms_used

  input:
  tuple val(setnames), file('ptms?'), file('stabileptms*'), file(cleaned_oldptms) from lucptmfiles_to_lup
  file('speclup.sql') from ptm_lookup_in

  output:
  path ptmtable into ptmpsmqc
  path 'ptmlup.sql' into ptmlup
  path({setnames.collect() { "${it}.tsv" }}) optional true into setptmtables
  path 'warnings' optional true into ptmwarnings

  script:
  ptmtable = "ptm_psmtable.txt"
  """
  # Concat all the PTM PSM tables (labile, stabile, previous) and load into DB
  # PSM peptide sequences include the PTM site
  cat speclup.sql > ptmlup.sql
  msstitch concat -i ${params.ptms && !is_rerun ? "stabileptms*" : ''} ${params.locptms && !is_rerun ? "ptms*" : ''} ${complementary_run ? "'${cleaned_oldptms}'" : ''} -o concatptmpsms
  msstitch psmtable -i concatptmpsms --dbfile ptmlup.sql -o "${ptmtable}" \
    --spectracol 1
  msstitch split -i "${ptmtable}" --splitcol bioset
  ${setnames.collect() { "test -f '${it}.tsv' || echo 'No PTMs found for set ${it}' >> warnings" }.join(' && ') }
  # No PTMs at all overwrites the per-set messages
  tail -n+2 ${ptmtable} | head | grep . >/dev/null || echo 'No PTMs found in any set' > warnings
  """
}

totalprotpsms_allsets
  .map { it -> listify(it) }
  .map{ it -> [it.collect() { it.baseName.replaceFirst(/\.tsv$/, "") }, it]} // setname from file setA.tsv
  .transpose()
  .set { totalprot_setpsms }

setptmtables
  .map { it -> listify(it) }
  .map{ it -> [it.collect() { it.baseName.replaceFirst(/\.tsv$/, "") }, it]} // setname from file setA.tsv
  .transpose()
  .set { preptmpeps }

if (params.totalproteomepsms) {
  preptmpeps.join(totalprot_setpsms)
    .set { ptm2peps }
} else {
  preptmpeps.set { ptm2peps }
}
  

process PTMPeptides {

  input:
  tuple val(setname), path('ptms.txt'), file('totalproteomepsms') from ptm2peps

  output:
  tuple val(setname), path(peptable), path("${peptable}_no_tp_normalized") into ptmpeps

  script:
  denom = !params.noquant && setdenoms ? setdenoms[setname] : false
  specialdenom = denom && (denom[0] == 'sweep' || denom[0] == 'intensity')
  peptable = "${setname}_ptm_peptides.txt"
  dividebycol = params.onlypeptides ? '^Protein$' : '^Gene Name'
  """
  # If there is a total proteome PSMs file, prepare a gene (or protein) and protein table
  # for peptide normalization purposes. Use msstitch isosummarize here so we dont have to deal 
  # with gene FDR and peptide tables etc. First the always non-normalized table to get
  # proteins with which to median-center the PTM table with (if --normalize is passed)
  ${params.totalproteomepsms && normalize ? "msstitch isosummarize -i totalproteomepsms -o prots_mediancenter \
    --featcol \$(head -n1 totalproteomepsms | tr '\\t' '\\n' | grep -n '${dividebycol}' | cut -f 1 -d':') \
    --isobquantcolpattern ${setisobaric[setname]} --minint 0.1 \
    ${denom && denom[0] == 'sweep' ? '--mediansweep --logisoquant': ''} \
    ${denom && params.keepnapsmsquant ? '--keep-psms-na-quant' : ''} \
    ${denom && !specialdenom ? "--logisoquant --denompatterns ${setdenoms[setname].join(' ')}": ''} \
" : ''}
  # PSM ID table is without master protein in onlypeptides
  ${params.onlypeptides ? "sed -i '0,/Protein/s//Protein ID/' prots_mediancenter" : ''}
  
  # Then as above but output the median-centered (if applicable) gene (default) or protein 
  # (in case of --onlypeptides) table for relating the peptides to their gene/protein
  ${params.totalproteomepsms && denom ? "msstitch isosummarize -i totalproteomepsms -o tp_accessions \
    --featcol \$(head -n1 totalproteomepsms | tr '\\t' '\\n' | grep -n '${dividebycol}' | cut -f 1 -d':') \
    --isobquantcolpattern ${setisobaric[setname]} --minint 0.1 \
    ${normalize ? "--median-normalize" : ''} \
    ${denom && denom[0] == 'sweep' ? '--mediansweep --logisoquant': ''} \
    ${denom && params.keepnapsmsquant ? '--keep-psms-na-quant' : ''} \
    ${denom && !specialdenom ? "--logisoquant --denompatterns ${setdenoms[setname].join(' ')}": ''} \
" : ''}
  ${params.onlypeptides ? "sed -i '0,/Protein/s//Protein ID/' tp_accessions" : ''}

  # Create a PTM-peptide table which has normalized isobaric data
  msstitch peptides -i "ptms.txt" -o "${peptable}" --scorecolpattern svm --spectracol 1 \
    ${!params.noquant && params.noms1quant ? '--ms1quantcolpattern area' : ''} \
    ${denom ? "--isobquantcolpattern ${setisobaric[setname]} --minint 0.1 --keep-psms-na-quant" : ''} \
    ${denom && denom[0] == 'sweep' ? '--mediansweep --logisoquant': ''} \
    ${denom && denom[0] == 'intensity' ? '--medianintensity' : ''} \
    ${denom && !specialdenom ? "--logisoquant --denompatterns ${setdenoms[setname].join(' ')}": ''} \
    ${denom && params.totalproteomepsms ? "--totalproteome tp_accessions" : ''} \
    ${denom && params.totalproteomepsms && normalize ? "--median-normalize --normalization-factors-table prots_mediancenter" : ''} \
  
  # And another one with the non-normalized isobaric quant (but with median-centering if specified)
  # this is not strictly necessary to do if there is no --totalproteomepsms, but will not be output 
  # in next step (merge) anyway then.
  msstitch peptides -i "ptms.txt" -o "${peptable}_no_tp_normalized" --scorecolpattern svm --spectracol 1 \
    ${denom ? "--isobquantcolpattern ${setisobaric[setname]} --minint 0.1 --keep-psms-na-quant" : ''} \
    ${denom && denom[0] == 'sweep' ? '--mediansweep --logisoquant': ''} \
    ${denom && denom[0] == 'intensity' ? '--medianintensity' : ''} \
    ${denom && !specialdenom ? "--logisoquant --denompatterns ${setdenoms[setname].join(' ')}": ''} \
    ${denom && params.totalproteomepsms && normalize ? "--median-normalize --normalization-factors-table prots_mediancenter" : ''} \
  """
}

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
  

process mergePTMPeps {
  publishDir "${params.outdir}", mode: 'copy', overwrite: true, saveAs: {it.startsWith('ptm_peptides_') ? it : null}
 
  // setnames is [false] if no PTM tables exist (otherwise this is fed with empty tables)
  when: setnames[0]

  input:
  tuple val(setnames), file(peptides), file(notp_adjust_peps), file('ptmlup.sql'), file('ptmpsms.txt') from ptmpeps2merge

  output:
  path peptable
  path peptable_no_adjust optional true
  tuple path("ptmqc.html"), path('summary.txt'), path('featcount_summary.txt') into ptmqc
  path('overlap.txt') into ptmoverlap optional true

  script:
  peptable = params.totalproteomepsms ? 'ptm_peptides_total_proteome_adjusted.txt' : 'ptm_peptides_not_adjusted.txt'
  peptable_no_adjust = 'ptm_peptides_not_adjusted.txt'
  """
  cat ptmlup.sql > pepptmlup.sql
  # Create first table, input for which is either adjusted or not
  msstitch merge -i ${listify(peptides).collect() { "$it" }.join(' ')} --setnames ${setnames.collect() { "'$it'" }.join(' ')} --dbfile pepptmlup.sql -o mergedtable --no-group-annotation \
    --fdrcolpattern '^q-value' --pepcolpattern 'peptide PEP' --flrcolpattern 'FLR' \
    ${!params.noquant && !params.noms1quant ? "--ms1quantcolpattern area" : ''} \
    ${!params.noquant && setisobaric ? "--isobquantcolpattern plex" : ''}
  # Add master/genes/gene count to peptide table, cant store in SQL because cant protein group on small PTM table
  ${!params.onlypeptides ? "head -n1 mergedtable | sed \$'s/Peptide sequence/Peptide sequence\tMaster protein(s)\tGene name(s)\tMatching gene count\tPTM flanking sequence(s)/'> '${peptable}'" : ''}
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
  ${!params.onlypeptides && params.totalproteomepsms ? """head -n1 mergedtable | sed \$'s/Peptide sequence/Peptide sequence\tMaster protein(s)\tGene name(s)\tMatching gene count\tPTM flanking sequence(s)/'> '${peptable_no_adjust}'
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

process makePeptides {
  input:
  set val(td), val(setname), file('psms') from psm_pep
  
  output:
  set val(setname), val(td), file(psms), file("${setname}_peptides") into prepgs_in
  set val(setname), val('peptides'), val(td), file("${setname}_peptides"), path(normfactors) optional true into peptides_out

  script:
  quant = !params.noquant && td == 'target'
  isoquant = quant && setisobaric && setisobaric[setname]
  denom = isoquant && setdenoms ? setdenoms[setname] : false
  specialdenom = denom && (denom[0] == 'sweep' || denom[0] == 'intensity')
  normfactors = "${setname}_normfacs"
  """
  # Create peptide table from PSM table, picking best scoring unique peptides
  msstitch peptides -i psms -o "${setname}_peptides" --scorecolpattern svm --spectracol 1 --modelqvals \
    ${quant ? "${!params.noms1quant ? '--ms1quantcolpattern area' : ''} ${isoquant ? "--isobquantcolpattern ${setisobaric[setname]} --minint 0.1" : ''}" : ''} \
    ${isoquant && params.keepnapsmsquant ? '--keep-psms-na-quant' : ''} \
    ${denom && denom[0] == 'sweep' ? '--mediansweep --logisoquant' : ''} \
    ${denom && denom[0] == 'intensity' ? '--medianintensity' : ''} \
    ${denom && !specialdenom ? "--logisoquant --denompatterns ${setdenoms[setname].join(' ')}" : ''} \
    ${quant && normalize ? "--median-normalize" : ''}
    ${quant && normalize ? "sed 's/^/$setname'\$'\t/' < normalization_factors_psms > '$normfactors'" : "touch '$normfactors'"}
  """
}


/*
* Step 4: Infer and quantify proteins and genes
*/

// Group set T-D combinations and remove those with only target or only decoy
pre_tprepgs_in = Channel.create()
dprepgs_in = Channel.create()
prepgs_in
  .groupTuple(by: 0) // group by setname/acctype
  .filter { it -> it[1].size() == 2 } // must have target and decoy ?
  .transpose()
  .choice(pre_tprepgs_in, dprepgs_in) { it[1] == 'target' ? 0 : 1 }
// combine target with fasta files
pre_tprepgs_in
  .combine(fdrdbs)
  .set { tprepgs_in }

process proteinGeneSymbolTableFDR {
  
  when: !params.onlypeptides
  input:
  set val(setname), val(td), file('tpsms'), file('tpeptides'), file(tfasta), file(dfasta) from tprepgs_in
  set val(setname), val(td), file('dpsms'), file('dpeptides') from dprepgs_in
  each acctype from acctypes

  output:
  set val(setname), val(acctype), file("${setname}_protfdr"), path(normfactors) into protfdrout
  file('warnings') optional true into fdrwarnings

  script:
  scorecolpat = acctype == 'proteins' ? '^q-value$' : 'linear model'
  denom = !params.noquant && setdenoms ? setdenoms[setname] : false
  specialdenom = denom && (denom[0] == 'sweep' || denom[0] == 'intensity')
  normfactors = "${setname}_normfacs"
  quant = !params.noquant && (!params.noms1quant || params.isobaric)
  isoquant = quant && setisobaric && setisobaric[setname]
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
    ${!params.noquant ? "${!params.noms1quant ? '--ms1quant' : ''} ${isoquant ? "--isobquantcolpattern ${setisobaric[setname]} --minint 0.1" : ''}" : ''} \
    ${quant ? '--psmtable tpsms' : ''} \
    ${isoquant && params.keepnapsmsquant ? '--keep-psms-na-quant' : ''} \
    ${denom && denom[0] == 'sweep' ? '--mediansweep --logisoquant' : ''} \
    ${denom && denom[0] == 'intensity' ? '--medianintensity' : ''} \
    ${denom && !specialdenom ? "--denompatterns ${setdenoms[setname].join(' ')} --logisoquant" : ''} \
    ${normalize ? "--median-normalize" : ''}
    ${normalize ? "sed 's/^/$setname'\$'\t/' < normalization_factors_tpsms > '$normfactors'" : "touch '$normfactors'"}
  """
}
    

percowarnings
  .concat(psmwarnings)
  .concat(fdrwarnings)
  .concat(ptmwarnings)
  .concat(luciphor_warnings)
  .toList()
  .set { warnings }

peptides_out
  .filter { it[2] == 'target' }
  // setname, acctype, outfile, optional normfactors
  .map { it -> [it[0], it[1]] + it[3..-1] }
  .set { tpeps }

tpeps
  .concat(protfdrout)
  .toSortedList( {a, b -> a[0] <=> b[0]} ) // sort on setname
  .flatMap { it -> it}
  .groupTuple(by: 1)  // all outputs of same accession type together.
  .set { ptables_to_merge }


psmlookup
  .filter { it[0] == 'target' }
  .collect()
  .map { it[1] }
  .set { tlookup }

/*
* Step 5: Create reports
*/

process proteinPeptideSetMerge {

  input:
  set val(setnames), val(acctype), file(tables), file(normfacs) from ptables_to_merge
  file(lookup) from tlookup
  file('sampletable') from Channel.from(sampletable).first()
  
  output:
  set val(acctype), file('proteintable'), file('sampletable') into featqc_extra_peptide_samples
  set val(acctype), file('proteintable'), file(normfacs) into merged_feats

  script:
  """
  # exchange sample names on isobaric fields in header
  # First add NO__GROUP marker for no-samplegroups clean sampletable from special chars
  ${params.sampletable ? 'awk -v FS="\\t" -v OFS="\\t" \'{if (NF==3) print $1,$2,$3,"NO__GROUP"; else print}\' sampletable > tmpsam && mv tmpsam sampletable' : ''}
  # Check if there are samplegroups at all
  ${params.deqms ? 'grep -v NO__GROUP sampletable || (>&2 echo "Cannot run DEqMS without specified sample groups" && exit 1)': ''}
  # Count amount samples per group and error on group with only one sample
  ${params.deqms ? "grep -v NO__GROUP sampletable | cut -f 4 | sort | uniq -c | sed 's/\\s*//' | grep '^1 ' && (>&2 echo 'Cannot run DEqMS when any sample groups have only one sample, please review input' && exit 1)" : ''}
  # strip lead/trail space in set name 
  paste <(cut -f1 sampletable) <(cut -f2 sampletable | sed "s/^\\s*//;s/\\s*\$//") <(cut -f3-4 sampletable) > nowhitespace && mv nowhitespace sampletable
  # substitute other weird characters
  ${params.sampletable ? 'sed "s/[^A-Za-z0-9_\\t]/_/g" sampletable > clean_sampletable' : ''}

  # SQLite lookup needs copying to not modify the input file which would mess up a rerun with -resume
  cat $lookup > db.sqlite
  msstitch merge -i ${listify(tables).collect() { "$it" }.join(' ')} --setnames ${setnames.collect() { "'$it'" }.join(' ')} --dbfile db.sqlite -o mergedtable \
    --fdrcolpattern '^q-value\$' ${acctype != 'peptides' ? "--mergecutoff ${params.proteinconflvl}" : ''} \
    ${acctype == 'peptides' ? "--pepcolpattern 'peptide PEP'" : ''} \
    ${!params.noquant && !params.noms1quant ? "--ms1quantcolpattern area" : ''} \
    ${!params.noquant && setisobaric ? "--isobquantcolpattern plex" : ''} \
    ${params.onlypeptides ? "--no-group-annotation" : ''}
   
  # Put annotation on header, use normal setname for finding, replace with clean name
  # "sed '0,/{RE}//{substitute}/..."  for only first line (0,/{RE} = read from 0 until RE,
  # then the empty // means use the previous RE (you could specify a new RE)
  head -n1 mergedtable > tmph
  ${params.sampletable && setisobaric ?  
    'while read line ; do read -a arr <<< $line ; sed -E "s/${arr[0]}_([a-z0-9]*plex)_${arr[1]}/${arr[4]}_${arr[3]}_${arr[2]}_\\1_${arr[1]}/" <(tail -n1 tmph | tr "\t" "\n") | tr "\n" "\t" | sed $"s/\\t$/\\n/" ; done < <(paste <(cut -f2 sampletable) clean_sampletable) >> tmph' :  ''}
  ${params.sampletable && setisobaric ? "cat <(tail -n1 tmph) <(tail -n+2 mergedtable) > grouptable" : 'mv mergedtable grouptable'}

  # Run DEqMS if needed, use original sample table with NO__GROUP
  ${params.deqms ? "numfields=\$(head -n1 grouptable | tr '\t' '\n' | wc -l) && deqms.R && paste <(head -n1 grouptable) <(head -n1 deqms_output | cut -f \$(( numfields+1 ))-\$(head -n1 deqms_output|wc -w)) > tmpheader && cat tmpheader <(tail -n+2 deqms_output) > proteintable" : 'mv grouptable proteintable'}
  """
}


psm_result
  .filter { it[0] == 'target' }
  .combine(scans_result)
  .map { it -> [it[0], it[1], it[2], it[3], it[4].unique()] }
  .set { targetpsm_result }


process psmQC {

  input:
  set val(td), file('psms'), file('filescans'), file('platescans'), val(plates) from targetpsm_result
  val(setnames) from setnames_psmqc
  set val(plicate_sets), val(mzmlpaths), val(mzmlplates), val(mzmlbasenames), val(fractions) from mzmlfiles_psmqc
  file oldmzmls_fn
  val oldmzml_sets

  output:
  set val('psms'), file('psmqc.html'), file('summary.txt') into psmqccollect
  val(plates) into qcplates
  // TODO no proteins == no coverage for pep centric
  script:
  """
  paste <(echo -e "${mzmlpaths.collect() { "${it.baseName}.${it.extension}" }.join('\\n')}") <(echo -e "${plicate_sets.join('\\n')}") <( echo -e "${mzmlplates.join('\\n')}") <(echo -e "${fractions.join('\\n')}") > mzmlfrs
  qc_psms.R ${setnames[0].size()} ${fractionation ? 'TRUE' : 'FALSE'} ${params.oldmzmldef ? oldmzmls_fn: 'nofile'} ${plates.join(' ')}
  # If any sets have zero (i.e. no PSMs, so no output from R), output them here by joining and filling in
  cat <(head -n1 psmtable_summary.txt) <(join -a1 -e0 -o auto -t \$'\\t' <(echo -e \'${plicate_sets.plus(oldmzml_sets).unique().join("\\n")}' | sort) <(tail -n+2 psmtable_summary.txt | sort -k1b,1)) > summary.txt
  sed -Ei 's/[^A-Za-z0-9_\\t]/./g' summary.txt
  echo "<html><body>" > psmqc.html
  for graph in psm-scans missing-tmt miscleav
    do
    [[ -e \$graph ]] && echo "<div class='chunk' id='\${graph}'>" \$(sed "s/id=\\"/id=\\"\${graph}/g;s/\\#/\\#\${graph}/g" < <(tail -n+2 \$graph)) "</div>" >> psmqc.html
    [[ -e \$graph ]] || echo "<div class='chunk' id='\${graph}'>None found</div>" >> psmqc.html
    done 
  for graph in retentiontime precerror fwhm fryield msgfscore pif
    do
    for plateid in ${plates.join(' ')}
      do
      plate="PLATE___\${plateid}___\${graph}"
    [[ -e \$plate ]] && echo "<div class='chunk \$plateid' id='\${graph}'>" \$(sed "s/id=\\"/id=\\"\${plate}/g;s/\\#/\\#\${plate}/g" < <(tail -n+2 \$plate)) "</div>" >> psmqc.html
      done 
    done
  echo "</body></html>" >> psmqc.html
  """
}

featqc_extra_peptide_samples
  .filter { it[0] == 'peptides' }
  .map { it -> [it[1], it[2]] }
  .set { featqc_peptides_samples }

merged_feats
  .combine(featqc_peptides_samples)
  .set { featqcinput }


process featQC {
  publishDir "${params.outdir}", mode: 'copy', overwrite: true, saveAs: {it == "feats" ? "${acctype}_table.txt": null}

  input:
  set val(acctype), file('feats'), file(normfacs), file(peptable), file(sampletable) from featqcinput
  val(setnames) from setnames_featqc

  output:
  file('feats') into featsout
  set val(acctype), file('featqc.html'), file('summary.txt'), file('overlap'), file('allnormfacs') into qccollect

  script:
  show_normfactors = setdenoms && normalize
  """
  # combine multi-set normalization factors
  cat ${normfacs} > allnormfacs
  # Create QC plots and put them base64 into HTML, R also creates summary.txt
  qc_protein.R --sets ${setnames.collect() { "'$it'" }.join(' ')} \
     --feattype ${acctype} --peptable $peptable \
     ${params.sampletable ? "--sampletable $sampletable" : ''} \
     --conflvl ${acctype == 'peptides' ? params.pepconflvl : params.proteinconflvl} \
     ${show_normfactors ? '--normtable allnormfacs' : ''}
  # Remove X from R's columns if they do not start with [A-Z]
  sed -Ei 's/^X([^A-Za-z])/\\1/' summary.txt
  echo "<html><body>" > featqc.html
  for graph in featyield precursorarea ${show_normfactors ? 'normfactors': ''} nrpsms nrpsmsoverlapping percentage_onepsm ms1nrpeps;
    do
    [ -e \$graph ] && echo "<div class=\\"chunk\\" id=\\"\${graph}\\"> \$(sed "s/id=\\"/id=\\"${acctype}-\${graph}/g;s/\\#/\\#${acctype}-\${graph}/g" < <(tail -n+2 \$graph)) </div>" >> featqc.html
    done 
    # coverage and isobaric plots are png because a lot of points
    [ -e isobaric ] && paste -d \\\\0  <(echo "<div class=\\"chunk\\" id=\\"isobaric\\"><img src=\\"data:image/png;base64,") <(base64 -w 0 isobaric) <(echo '"></div>') >> featqc.html
    [ -e coverage ] && paste -d \\\\0  <(echo "<div class=\\"chunk\\" id=\\"coverage\\"><img src=\\"data:image/png;base64,") <(base64 -w 0 coverage) <(echo '"></div>') >> featqc.html
  # Fetch special (multi-pane) DEqMS and PCA plots
  # Use ls to check because wildcard doesnt work in -e
  ls deqms_volcano_* && echo '<div class="chunk" id="deqms">' >> featqc.html
  for graph in deqms_volcano_*;
    do
    paste -d \\\\0  <(echo '<div><img src="data:image/png;base64,') <(base64 -w 0 \$graph) <(echo '"></div>') >> featqc.html
    done
  ls deqms_volcano_* && echo '</div>' >> featqc.html
  [ -e pca ] && echo '<div class="chunk" id="pca">' >> featqc.html && for graph in pca scree;
    do 
    echo "<div> \$(sed "s/id=\\"/id=\\"${acctype}-\${graph}/g;s/\\#/\\#${acctype}-\${graph}/g" < <(tail -n+2 \$graph)) </div>" >> featqc.html
    done
    [ -e pca ] && echo '</div>' >> featqc.html

  echo "</body></html>" >> featqc.html

  # Create overlap table
  qcols=\$(head -n1 feats |tr '\\t' '\\n'|grep -n "_q-value"| tee nrsets | cut -f 1 -d ':' |tr '\\n' ',' | sed 's/\\,\$//')
  protcol=\$(head -n1 feats | tr '\\t' '\\n' | grep -n Protein | grep -v start | cut -f1 -d ':')
  ${acctype == 'peptides' ? 'cut -f1,"\$qcols","\$protcol" feats | grep -v ";" > tmpqvals' : 'cut -f1,"\$qcols" feats > qvals'}
  ${acctype == 'peptides' ? 'nonprotcol=\$(head -n1 tmpqvals | tr "\\t" "\\n" |grep -vn Protein | cut -f1 -d":" | tr "\\n" "," | sed "s/\\,\$//") && cut -f"\$nonprotcol" tmpqvals > qvals' : ''}
  nrsets=\$(wc -l nrsets | sed 's/\\ .*//')
  # read lines, sed removes all non-A chars so only N from NA is left.
  while read line ; do 
  	nr=\$(printf "\$line" |wc -m)  # Count NA
  	overlap=\$(( \$nrsets-\$nr )) # nrsets minus NAcount is the overlap
  	echo "\$overlap" >> setcount
  done < <(tail -n+2 qvals | cut -f2- | sed 's/[^A]//g' )
  echo nr_sets\$'\t'nr_${acctype} > overlap
  for num in \$(seq 1 \$nrsets); do 
  	echo "\$num"\$'\t'\$( grep ^"\$num"\$ setcount | wc -l) >> overlap
  done

  # Remove internal no-group identifier so it isnt output
  sed '0,/NO__GROUP_/s///g' -i feats
  """
}

noptms = ['NA', 'NA', 'NA']
ptmqcout = params.ptms || params.locptms || params.ptmpsms ? ptmqc.ifEmpty(noptms) : Channel.value(noptms)

qccollect
  .concat(psmqccollect)
  .toList()
  .map { it -> [it.collect() { it[0] }, it.collect() { it[1] }, it.collect() { it[2] }, it.collect() { it[3] }, it.collect() { it[4] } ] }
  .merge(ptmqcout)
  .set { collected_feats_qc }


process collectQC {

  publishDir "${params.outdir}", mode: 'copy', overwrite: true

  input:
  set val(acctypes), file('feat?'), file('summary?'), file('overlap?'), file('normfacs?'), file(ptmqc), file(ptmsummary), file(ptmfeatsummary) from collected_feats_qc
  val(plates) from qcplates
  file('ptmoverlap') from ptmoverlap.ifEmpty('false')
  file('sw_ver') from software_versions_qc
  file('warnings??') from warnings

  output:
  set file('qc_light.html'), file('qc_full.html')

  script:
  """
  count=1; for ac in ${acctypes.join(' ')}; 
    do mv feat\$count \$ac.html;
    mv summary\$count \${ac}_summary;
    mv overlap\$count \${ac}_overlap 
    [ -s normfacs\${count} ] && cat <(echo Set\$'\t'channel\$'\t'normfac) normfacs\$count > \${ac}_normfacs;
    ((count++)); done
  join -j 1 -o auto -t '\t' <(head -n1 psms_summary) <(head -n1 peptides_summary) > psmpepsum_header
  join -j 1 -o auto -t '\t' <(tail -n+2 psms_summary | sort -k1b,1 ) <(tail -n+2 peptides_summary | sort -k1b,1 ) > psmpepsum_tab

  # onlypeptides param, or failure to detect proteins, makes a peptide summary, else also add proteins
  if [[ -e  proteins_summary ]]
  then
    join -j 1 -o auto -t \'\t\' psmpepsum_tab <(sort -k1b,1 <(tail -n+2 proteins_summary)) > pepprotsum_tab && join -j 1 -o auto -t \'\t\' psmpepsum_header <(head -n1 proteins_summary) > pepprotsum_head
  else
    cat psmpepsum_header psmpepsum_tab | tee summary pre_summary_light_tab
    awk -v FS='\\t' -v OFS='\\t' '{print \$1,\$3,\$2}' pre_summary_light_tab > summary_light
  fi

  # in case of genes, join those on the prot/pep tables (full summary) and psmpeptables (light summary), else passthrough those to summaries
  if [[ -e genes_summary ]]
  then
    join -j 1 -o auto -t \'\t\' pepprotsum_tab <( sort -k1b,1 <( tail -n+2 genes_summary)) > summary_tab
    join -j 1 -o auto -t \'\t\' pepprotsum_head <(head -n1 genes_summary) > summary_head
    cat summary_head summary_tab > summary
    join -j 1 -o auto -t \'\t\' psmpepsum_tab <( sort -k1b,1 <(tail -n+2 genes_summary)) > summary_light_tab
    join -j 1 -o auto -t \'\t\' psmpepsum_header <( head -n1 genes_summary) > summary_light_head && cat summary_light_head summary_light_tab > summary_light
  else
    if [[ -e proteins_summary ]]
    then
      cat pepprotsum_head pepprotsum_tab | tee summary summary_light
    fi
  fi


  # remove Yaml from software_versions to get HTML
  grep '<.*>' sw_ver > sw_ver_cut
  
  # merge warnings
  ls warnings* && cat warnings* > warnings.txt
  # collect and generate HTML report
  qc_collect.py $baseDir/assets/qc_full.html $params.name ${fractionation ? "frac" : "nofrac"} ${ptmqc.exists() ? "$ptmqc:$ptmsummary:$ptmfeatsummary" : 'noptm'} ${plates.join(' ')}
  qc_collect.py $baseDir/assets/qc_light.html $params.name ${fractionation ? "frac" : "nofrac"} ${ptmqc.exists() ? "$ptmqc:$ptmsummary:$ptmfeatsummary" : 'noptm'} ${plates.join(' ')}
  """
}


/* 
 * STEP 3 - Output Description HTML
*/
process output_documentation {
    tag "$prefix"

    publishDir "${params.outdir}/Documentation", mode: 'copy', overwrite: true

    input:
    file output_docs

    output:
    file "results_description.html"

    script:
    """
    markdown_to_html.r $output_docs results_description.html
    """
}



/*
 * Completion e-mail notification
 */
workflow.onComplete {

    // Set up the e-mail variables
    def subject = "[lehtiolab/ddamsproteomics] Successful: $workflow.runName"
    if(!workflow.success){
      subject = "[lehtiolab/ddamsproteomics] FAILED: $workflow.runName"
    }
    sw_versions = file("${params.outdir}/software_versions.yaml").readLines().grep(~/ *<dt.+dd> */).collect { it.replaceAll('^[ ]*<dt>', '').replaceAll('</dd>[ ]*$', '').split('</dt><dd>') }

    def email_fields = [:]
    email_fields['version'] = workflow.manifest.version
    email_fields['runName'] = custom_runName ?: workflow.runName
    email_fields['success'] = workflow.success
    email_fields['dateComplete'] = workflow.complete
    email_fields['duration'] = workflow.duration
    email_fields['exitStatus'] = workflow.exitStatus
    email_fields['errorMessage'] = (workflow.errorMessage ?: 'None')
    email_fields['errorReport'] = (workflow.errorReport ?: 'None')
    email_fields['commandLine'] = workflow.commandLine
    email_fields['projectDir'] = workflow.projectDir
    email_fields['summary'] = summary
    email_fields['sw_versions'] = sw_versions
    email_fields['summary']['Date Started'] = workflow.start
    email_fields['summary']['Date Completed'] = workflow.complete
    email_fields['summary']['Pipeline script file path'] = workflow.scriptFile
    email_fields['summary']['Pipeline script hash ID'] = workflow.scriptId
    if(workflow.repository) email_fields['summary']['Pipeline repository Git URL'] = workflow.repository
    if(workflow.commitId) email_fields['summary']['Pipeline repository Git Commit'] = workflow.commitId
    if(workflow.revision) email_fields['summary']['Pipeline Git branch/tag'] = workflow.revision
    email_fields['summary']['Nextflow Version'] = workflow.nextflow.version
    email_fields['summary']['Nextflow Build'] = workflow.nextflow.build
    email_fields['summary']['Nextflow Compile Timestamp'] = workflow.nextflow.timestamp

    // Render the TXT template
    def engine = new groovy.text.GStringTemplateEngine()
    def tf = new File("$baseDir/assets/email_template.txt")
    def txt_template = engine.createTemplate(tf).make(email_fields)
    def email_txt = txt_template.toString()

    // Render the HTML template
    def hf = new File("$baseDir/assets/email_template.html")
    def html_template = engine.createTemplate(hf).make(email_fields)
    def email_html = html_template.toString()

    // Render the sendmail template
    def smail_fields = [ email: params.email, subject: subject, email_txt: email_txt, email_html: email_html, baseDir: "$baseDir" ]
    def sf = new File("$baseDir/assets/sendmail_template.txt")
    def sendmail_template = engine.createTemplate(sf).make(smail_fields)
    def sendmail_html = sendmail_template.toString()

    // Send the HTML e-mail
    if (params.email) {
        try {
          if( params.plaintext_email ){ throw GroovyException('Send plaintext e-mail, not HTML') }
          // Try to send HTML e-mail using sendmail
          [ 'sendmail', '-t' ].execute() << sendmail_html
          log.info "[lehtiolab/ddamsproteomics] Sent summary e-mail to $params.email (sendmail)"
        } catch (all) {
          // Catch failures and try with plaintext
          [ 'mail', '-s', subject, params.email ].execute() << email_txt
          log.info "[lehtiolab/ddamsproteomics] Sent summary e-mail to $params.email (mail)"
        }
    }

    // Write summary e-mail HTML to a file
    def output_d = new File( "${params.outdir}/Documentation/" )
    if( !output_d.exists() ) {
      output_d.mkdirs()
    }
    def output_hf = new File( output_d, "pipeline_report.html" )
    output_hf.withWriter { w -> w << email_html }
    def output_tf = new File( output_d, "pipeline_report.txt" )
    output_tf.withWriter { w -> w << email_txt }

    log.info "[lehtiolab/ddamsproteomics] Pipeline Complete"

}
