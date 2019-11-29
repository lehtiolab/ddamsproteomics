#!/usr/bin/env nextflow
/*
========================================================================================
                         nf-core/ddamsproteomics
========================================================================================
 nf-core/ddamsproteomics Analysis Pipeline.
 #### Homepage / Documentation
 https://github.com/nf-core/ddamsproteomics
----------------------------------------------------------------------------------------
*/


def helpMessage() {
    log.info"""
    =========================================
     nf-core/ddamsproteomics v${workflow.manifest.version}
    =========================================
    Usage:

    The typical command for running the pipeline is as follows:

    nextflow run nf-core/ddamsproteomics --mzmls '*.mzML' --tdb swissprot_20181011.fa --mods assets/tmtmods.txt -profile standard,docker

    Mandatory arguments:
      --mzmls                       Path to mzML files
      --mzmldef                     Alternative to --mzml: path to file containing list of mzMLs 
                                    with sample set and fractionation annotation (see docs)
      --tdb                         Path to target FASTA protein database
      --mods                        Path to MSGF+ modification file (two examples in assets folder)
      -profile                      Configuration profile to use. Can use multiple (comma separated)
                                    Available: standard, conda, docker, singularity, awsbatch, test

    Options:
      --isobaric VALUE              In case of isobaric, specify: tmt10plex, tmt6plex, itraq8plex, itraq4plex
      --prectol                     Precursor error for search engine (default 10ppm)
      --iso_err                     Isotope error for search engine (default -1,2)
      --frag                        Fragmentation method for search engine (default 'auto')
      --enzyme                      Enzyme used, default trypsin, pick from:
                                    unspecific, trypsin, chymotrypsin, lysc, lysn, gluc, argc, aspn, no_enzyme
      --terminicleaved              Allow only 'full', 'semi' or 'non' cleaved peptides
      --maxmiscleav		    Maximum allowed amount of missed cleavages for MSGF+
      --minpeplen                   Minimum peptide length to search
      --maxpeplen                   Maximum peptide length to search
      --mincharge                   Minimum peptide charge search
      --maxcharge                   Maximum peptide charge search
      --activation VALUE            Specify activation protocol: hcd (DEFAULT), cid, etd for isobaric 
                                    quantification. Not necessary for other functionality.
      --normalize                   Normalize isobaric values by median centering on channels of protein table
      --sampletable                 Path to sample annotation table in case of isobaric analysis
      --deqms                       Perform DEqMS differential expression analysis using sampletable
      --genes                       Produce gene table (i.e. ENSG or gene names from Swissprot)
      --symbols                     Produce gene symbols table (i.e. gene names when using ENSEMBL DB)
      --martmap FILE                Necessary when using ENSEMBL FASTA database, tab-separated file 
                                    with information from Biomart. An example can be found at
                                    https://github.com/nf-core/test-datasets/raw/ddamsproteomics/testdata/
      --fractions                   Fractionated samples, 
      --hirief                      IEF fractionated samples, implies --fractions, allows delta pI calculation
      --pipep FILE                  File containing peptide sequences and their isoelectric points. Example
                                    can be found in https://github.com/nf-core/test-datasets/raw/ddamsproteomics/
      --onlypeptides                Do not produce protein or gene level data
      --noquant                     Do not produce isobaric or MS1 quantification data
      --quantlookup FILE            Use previously generated SQLite lookup database containing spectra 
                                    quantification data when e.g. re-running. Need to match exactly to the
                                    mzML files of the current run
      --fastadelim VALUE            FASTA header delimiter in case non-standard FASTA is used, to be used with
                                    --genefield
      --genefield VALUE             Number to determine in which field of the FASTA header (split 
                                    by --fastadelim) the gene name can be found.


    Other options:
      --outdir                      The output directory where the results will be saved
      --email                       Set this parameter to your e-mail address to get a summary e-mail with details of the run sent to you when the workflow exits
      -name                         Name for the pipeline run. If not specified, Nextflow will automatically generate a random mnemonic.

    AWSBatch options:
      --awsqueue                    The AWSBatch JobQueue that needs to be set when running on AWSBatch
      --awsregion                   The AWS Region for your AWS Batch job to run on
    """.stripIndent()
}

/*
 * SET UP CONFIGURATION VARIABLES
 */

// Show help emssage
if (params.help){
    helpMessage()
    exit 0
}

// Configurable variables
params.name = false
params.email = false
params.plaintext_email = false

params.mzmls = false
params.martmap = false
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
params.fdrmethod = 'tdconcat'
params.activation = 'hcd' // Only for isobaric quantification
params.outdir = 'results'
params.normalize = false
params.genes = false
params.symbols = false
params.fastadelim = false
params.genefield = false
params.quantlookup = false
params.fractions = false
params.hirief = false
params.pipep = false
params.onlypeptides = false
params.noquant = false
params.denoms = false
params.sampletable = false
params.deqms = false

// Validate and set file inputs
fractionation = (params.hirief || params.fractions)
mods = file(params.mods)
if( !mods.exists() ) exit 1, "Modification file not found: ${params.mods}"
tdb = file(params.tdb)
if( !tdb.exists() ) exit 1, "Target fasta DB file not found: ${params.tdb}"

// Files which are not standard can be checked here
if (params.martmap) {
  martmap = file(params.martmap)
  if( !martmap.exists() ) exit 1, "Biomart ENSEMBL mapping file not found: ${params.martmap}"
}
if (params.pipep) {
  trainingpep = file(params.pipep)
  if( !trainingpep.exists() ) exit 1, "Peptide pI data file not found: ${params.pipep}"
} else { trainingpep = false }

if (params.sampletable) {
  sampletable = file(params.sampletable)
  if( !sampletable.exists() ) exit 1, "Sampletable file not found: ${params.sampletable}"
} else {
  sampletable = 0
}

output_docs = file("$baseDir/docs/output.md")

// set constant variables
accolmap = [peptides: 12, proteins: 14, genes: 17, assoc: 18]


// parse inputs that combine to form values or are otherwise more complex.
setdenoms = [:]
if (!(params.noquant) && params.isobaric && params.denoms) {
  params.denoms.tokenize(' ').each{ it -> x=it.tokenize(':'); setdenoms.put(x[0], x[1..-1])}
}

plextype = params.isobaric ? params.isobaric.replaceFirst(/[0-9]+plex/, "") : 'false'
normalize = (!params.noquant && params.normalize && params.isobaric)
rawisoquant = (!params.noquant && !params.normalize && params.isobaric)


// AWSBatch sanity checking
if(workflow.profile == 'awsbatch'){
    if (!params.awsqueue || !params.awsregion) exit 1, "Specify correct --awsqueue and --awsregion parameters on AWSBatch!"
    if (!workflow.workDir.startsWith('s3') || !params.outdir.startsWith('s3')) exit 1, "Specify S3 URLs for workDir and outdir parameters on AWSBatch!"
}


// Has the run name been specified by the user?
//  this has the bonus effect of catching both -name and --name
custom_runName = params.name
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

nf-core/ddamsproteomics v${workflow.manifest.version}"
======================================================="""
def summary = [:]
summary['Pipeline Name']  = 'nf-core/ddamsproteomics'
summary['Pipeline Version'] = workflow.manifest.version
summary['Run Name']     = custom_runName ?: workflow.runName
summary['mzMLs']        = params.mzmls
summary['Target DB']    = params.tdb
summary['Sample annotations'] = params.sampletable
summary['Modifications'] = params.mods
summary['Instrument'] = params.instrument
summary['Precursor tolerance'] = params.prectol
summary['Isotope error'] = params.iso_err
summary['Fragmentation method'] = params.frag
summary['Enzyme'] = params.enzyme
summary['Allowed peptide termini cleavage'] = params.terminicleaved
summary['Allowed amount of missed cleavages'] = params.maxmiscleav
summary['Minimum peptide length'] = params.minpeplen
summary['Maximum peptide length'] = params.maxpeplen
summary['Minimum peptide charge'] = params.mincharge
summary['Maximum peptide charge'] = params.maxcharge
summary['FDR method'] = params.fdrmethod
summary['Isobaric tags'] = params.isobaric
summary['Isobaric activation'] = params.activation
summary['Isobaric normalization'] = params.normalize
summary['Output genes'] = params.genes
summary['Output symbols'] = params.symbols
summary['Custom FASTA delimiter'] = params.fastadelim 
summary['Custom FASTA gene field'] = params.genefield
summary['Premade quant data SQLite'] = params.quantlookup
summary['Fractionated sample'] = fractionation
summary['HiRIEF'] = params.hirief 
summary['peptide pI data'] = params.pipep
summary['Only output peptides'] = params.onlypeptides
summary['Do not quantify'] = params.noquant
summary['Perform DE analysis'] = params.deqms
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


def create_workflow_summary(summary) {

    def yaml_file = workDir.resolve('workflow_summary_mqc.yaml')
    yaml_file.text  = """
    id: 'nf-core-ddamsproteomics-summary'
    description: " - this information is collected when the pipeline is started."
    section_name: 'nf-core/ddamsproteomics Workflow Summary'
    section_href: 'https://github.com/nf-core/ddamsproteomics'
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

    publishDir "${params.outdir}", mode: 'copy'

    output:
    file 'software_versions.yaml' into software_versions_qc

    script:
    """
    echo $workflow.manifest.version > v_pipeline.txt
    echo $workflow.nextflow.version > v_nextflow.txt
    msgf_plus | head -n1 > v_msgf.txt
    hardklor | head -n1 > v_hk.txt || true
    kronik | head -n2 > v_kr.txt
    percolator -h |& head -n1 > v_perco.txt || true
    msspsmtable --version > v_mss.txt
    source activate openms-2.4.0
    IsobaricAnalyzer |& grep Version > v_openms.txt || true
    scrape_software_versions.py > software_versions.yaml
    """
}

if (workflow.profile.tokenize(',').intersect(['test', 'test_nofrac'])) { 
  // Profile 'test' delivers mzmlPaths
  Channel
    .from(params.mzmlPaths)
    .set { mzml_in }
}
else if (!params.mzmldef) {
  Channel
    .fromPath(params.mzmls)
    .map { it -> [it, 'NA'] }
    .set { mzml_in }
} else {
  Channel
    .from(file("${params.mzmldef}").readLines())
    .map { it -> it.tokenize('\t') }
    .set { mzml_in }
}


def or_na(it, length){
    return it.size() > length ? it[length] : 'NA'
}
// Parse mzML input to get files and sample names etc
// get setname, sample name (baseName), input mzML file. 
// Set platename to samplename if not specified. 
// Set fraction name to NA if not specified
mzml_in
  .tap { mzmlfiles_counter } // for counting, so config can set time limit
  .map { it -> [it[1], file(it[0]).baseName.replaceFirst(/.*\/(\S+)\.mzML/, "\$1"), file(it[0]), it[2] ? it[2] : it[1], it[3] ? it[3] : 'NA' ]}
  .into { sets; strips; mzmlfiles; mzml_quant; mzml_msgf }

// Set names are first item in input lists, collect them for PSM tables and QC purposes
sets
  .map{ it -> it[0] }
  .unique()
  .tap { setnames_psm } 
  .collect()
  .map { it -> [it] }
  .into { setnames_featqc; setnames_psmqc }

// Strip names for HiRIEF fractionation are third item, 
strips
  .map { it -> it[3] }
  .unique()
  .toList()
  .set { strips_for_deltapi }


/*
* Step 1: Extract quant data from peptide spectra
*/

process quantifySpectra {
  when: !params.quantlookup && !params.noquant

  input:
  set val(setname), val(sample), file(infile), val(platename), val(fraction) from mzml_quant
  file(hkconf) from Channel.fromPath("$baseDir/assets/hardklor.conf").first()

  output:
  set val(sample), file("${sample}.kr"), file(infile) into kronik_out
  set val(sample), file("${infile}.consensusXML") optional true into isobaricxml

  script:
  activationtype = [hcd:'High-energy collision-induced dissociation', cid:'Collision-induced dissociation', etd:'Electron transfer dissociation'][params.activation]
  massshift = [tmt:0.0013, itraq:0.00125, false:0][plextype]
  """
  # Run hardklor on config file with added line for in/out files
  # then run kronik on hardklor and quant isobaric labels if necessary
  hardklor <(cat $hkconf <(echo "$infile" hardklor.out))
  kronik -c 5 -d 3 -g 1 -m 8000 -n 600 -p 10 hardklor.out ${sample}.kr
  source activate openms-2.4.0
  ${params.isobaric ? "IsobaricAnalyzer  -type $params.isobaric -in $infile -out \"${infile}.consensusXML\" -extraction:select_activation \"$activationtype\" -extraction:reporter_mass_shift $massshift -extraction:min_precursor_intensity 1.0 -extraction:keep_unannotated_precursor true -quantification:isotope_correction true" : ''}
  """
}


// Collect all mzMLs into single item to pass to lookup builder and spectra counter
mzmlfiles
  .toList()
  .map { it.sort( {a, b -> a[1] <=> b[1]}) } // sort on sample for consistent .sh script in -resume
  .map { it -> [it.collect() { it[0] }, it.collect() { it[2] }, it.collect() { it[3] } ] } // lists: [sets], [mzmlfiles], [plates]
  .into { mzmlfiles_all; mzmlfiles_all_count }

mzmlfiles_counter
  .count()
  .subscribe { println "$it mzML files in analysis" }
  .into { mzmlcount_psm; mzmlcount_percolator }

process createSpectraLookup {

  when: !params.quantlookup

  input:
  set val(setnames), file(mzmlfiles), val(platenames) from mzmlfiles_all

  output:
  file 'mslookup_db.sqlite' into newspeclookup 

  script:
  """
  msslookup spectra -i ${mzmlfiles.join(' ')} --setnames ${setnames.join(' ')}
  """
}


// Collect all isobaric quant XML output for quant lookup building process
isobaricxml
  .ifEmpty(['NA', 'NA', 'NA'])
  .toList()
  .map { it.sort({a, b -> a[0] <=> b[0]}) }
  .map { it -> [it.collect() { it[0] }, it.collect() { it[1] }] } // samples, isoxml
  .set { isofiles_sets }

// Collect all MS1 kronik output for quant lookup building process
kronik_out
  .ifEmpty(['NA', 'NA'])
  .toList()
  .map { it.sort({a, b -> a[0] <=> b[0]}) }
  .map { it -> [it.collect() { it[0] }, it.collect() { it[1] }, it.collect() { it[2] }] } // samples, kronikout, mzml
  .set { krfiles_sets }


// Need to populate channels depending on if a pre-made quant lookup has been passed
// even if not needing quant (--noquant) this is necessary or NF will error
if (params.noquant && !params.quantlookup) {
  newspeclookup
    .into { quant_lookup; spec_lookup; countlookup }
} else if (!params.quantlookup) {
  newspeclookup
    .into { spec_lookup; countlookup }
} else {
  Channel
    .fromPath(params.quantlookup)
    .into { quant_lookup; countlookup }
  Channel.empty().set { spec_lookup }
} 


process quantLookup {

  publishDir "${params.outdir}", mode: 'copy', overwrite: true, saveAs: {it == 'db.sqlite' ? 'quant_lookup.sql' : null }

  when: !params.quantlookup && !params.noquant

  input:
  file lookup from spec_lookup
  set val(isosamples), file(isofns) from isofiles_sets
  set val(krsamples), file(krfns), file(mzmls) from krfiles_sets

  output:
  file('db.sqlite') into newquantlookup

  script:
  if (params.isobaric)
  """
  # SQLite lookup needs copying to not modify the input file which would mess up a rerun with -resume
  cat $lookup > db.sqlite
  msslookup ms1quant --dbfile db.sqlite -i ${krfns.join(' ')} --spectra ${mzmls.join(' ')} --quanttype kronik --mztol 20.0 --mztoltype ppm --rttol 5.0 
  msslookup isoquant --dbfile db.sqlite -i ${isofns.join(' ')} --spectra ${isosamples.collect{ x -> x + '.mzML' }.join(' ')}
  """
  else
  """
  # SQLite lookup needs copying to not modify the input file which would mess up a rerun with -resume
  cat $lookup > db.sqlite
  msslookup ms1quant --dbfile db.sqlite -i ${krfns.join(' ')} --spectra ${mzmls.join(' ')} --quanttype kronik --mztol 20.0 --mztoltype ppm --rttol 5.0 
  """
}


if (!params.quantlookup && !params.noquant) {
  newquantlookup
    .set { quant_lookup }
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
  sqlite3 $speclookup "SELECT mzmlfilename, COUNT(*) FROM mzml JOIN mzmlfiles USING(mzmlfile_id) JOIN biosets USING(set_id) GROUP BY mzmlfilename" > amount_spectra_files
  """
}


if (fractionation) { 
  specfilems2.set { scans_platecount }
} else {
  specfilems2
    .map { it -> [it[3], ['noplates']] }
    .into { scans_platecount; scans_result }
}


process countMS2sPerPlate {

  publishDir "${params.outdir}", mode: 'copy', overwrite: true 
  when: fractionation

  input:
  set val(setnames), file(mzmlfiles), val(platenames), file('nr_spec_per_file') from scans_platecount

  output:
  set file('scans_per_plate'), val(splates) into scans_perplate

  script:
  splates = [setnames, platenames].transpose().collect() { "${it[0]}_${it[1]}" }
  """
  #!/usr/bin/env python
  platesets = [\"${splates.join('", "')}\"]
  platescans = {p: 0 for p in platesets}
  fileplates = {fn: p for fn, p in zip([\"${mzmlfiles.join('", "')}\"], platesets)}
  with open('nr_spec_per_file') as fp:
      for line in fp:
          fn, scans = line.strip('\\n').split('|')
          platescans[fileplates[fn]] += int(scans)
  with open('scans_per_plate', 'w') as fp:
      for plate, scans in platescans.items():
          fp.write('{}\\t{}\\n'.format(plate, scans))
  """
}

if (fractionation) {
  scans_perplate.set { scans_result }
}

/*
* Step 2: Identify peptides
*/

process createTargetDecoyFasta {
 
  input:
  file(tdb)

  output:
  file('db.fa') into concatdb
  set file(tdb), file("decoy.fa") into searchdbs 

  script:
  """
  msslookup makedecoy -i "$tdb" -o decoy.fa --scramble tryp_rev --ignore-target-hits
  cat "$tdb" decoy.fa > db.fa
  """
}


process msgfPlus {
  cpus = config.poolSize < 4 ? config.poolSize : 4

  input:
  set val(setname), val(sample), file(x), val(platename), val(fraction) from mzml_msgf
  file(db) from concatdb
  file mods

  output:
  set val(setname), val(sample), file("${sample}.mzid") into mzids
  set val(setname), file("${sample}.mzid"), file("${sample}.mzid.tsv") into mzidtsvs
  
  script:
  msgfprotocol = [tmt:4, itraq:2, false:0][plextype]
  msgfinstrument = [velos:1, qe:3, false:0][params.instrument]
  fragmeth = [auto:0, cid:1, etd:2, hcd:3, uvpd:4][params.frag]
  enzyme = params.enzyme.indexOf('-') > -1 ? params.enzyme.replaceAll('-', '') : params.enzyme
  enzyme = [unspecific:0, trypsin:1, chymotrypsin: 2, lysc: 3, lysn: 4, gluc: 5, argc: 6, aspn:7, no_enzyme:9][enzyme]
  ntt = [full: 2, semi: 1, non: 0][params.terminicleaved]

  """
  msgf_plus -Xmx8G -d $db -s $x -o "${sample}.mzid" -thread ${task.cpus * params.threadspercore} -mod $mods -tda 0 -maxMissedCleavages $params.maxmiscleav -t ${params.prectol}  -ti ${params.iso_err} -m ${fragmeth} -inst ${msgfinstrument} -e ${enzyme} -protocol ${msgfprotocol} -ntt ${ntt} -minLength ${params.minpeplen} -maxLength ${params.maxpeplen} -minCharge ${params.mincharge} -maxCharge ${params.maxcharge} -n 1 -addFeatures 1
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
  set val(setname), val(samples), file('mzid?') from mzids_2pin
  val(mzmlcount) from mzmlcount_percolator

  output:
  set val(setname), file('perco.xml') into percolated

  script:
  """
  mkdir mzids
  count=1;for sam in ${samples.join(' ')}; do ln -s `pwd`/mzid\$count mzids/\${sam}.mzid; echo mzids/\${sam}.mzid >> metafile; ((count++));done
  msgf2pin -o percoin.tsv -e ${params.enzyme} -P "decoy_" metafile
  percolator -j percoin.tsv -X perco.xml -N 500000 --decoy-xml-output
  """
}


mzidtsvs
  .groupTuple()
  .join(percolated)
  .set { mzperco }


process fdrToTSV {

  input:
  set val(setname), file(mzids), file(tsvs), file(perco) from mzperco

  output:
  set val(setname), val('target'), file('target.tsv') into tmzidtsv_perco
  set val(setname), val('decoy'), file('decoy.tsv') into dmzidtsv_perco

  script:
  if (params.fdrmethod == 'tdconcat')
  """
  mkdir outtables
  msspsmtable percolator --perco $perco -d outtables -i ${tsvs.collect() { "'$it'" }.join(' ')} --mzids ${mzids.collect() { "'$it'" }.join(' ')}
  msspsmtable merge -i outtables/* -o psms
  msspsmtable split -i psms --splitcol \$(head -n1 psms | tr '\t' '\n' | grep -n ^TD\$ | cut -f 1 -d':')
  """
}

// Collect percolator data of target/decoy and feed into PSM table creation
tmzidtsv_perco
  .concat(dmzidtsv_perco)
  .groupTuple(by: 1)
  .combine(quant_lookup)
  .set { prepsm }



/*
* Step 3: Post-process peptide identification data
*/

process createPSMTable {

  publishDir "${params.outdir}", mode: 'copy', overwrite: true, saveAs: {["target_psmlookup.sql", "target_psmtable.txt", "decoy_psmtable.txt"].contains(it) ? it : null}

  input:
  set val(setnames), val(td), file('psms?'), file('lookup') from prepsm
  set file(tdb), file(ddb) from searchdbs
  val(allstrips) from strips_for_deltapi
  file(trainingpep) 
  val(mzmlcount) from mzmlcount_psm

  output:
  set val(td), file("${outpsms}") into psm_result
  set val(td), file({setnames.collect() { "${it}.tsv" }}) optional true into setpsmtables
  set val(td), file("${psmlookup}") into psmlookup
  file('warnings') optional true into psmwarnings

  script:
  psmlookup = "${td}_psmlookup.sql"
  outpsms = "${td}_psmtable.txt"

  """
  msspsmtable merge -i psms* -o psms.txt
  msspsmtable conffilt -i psms.txt -o filtpsm --confidence-better lower --confidence-lvl $params.psmconflvl --confcolpattern 'PSM q-value'
  msspsmtable conffilt -i filtpsm -o filtpep --confidence-better lower --confidence-lvl $params.pepconflvl --confcolpattern 'peptide q-value'
  tail -n+2 filtpep | grep . || (echo "No ${td} PSMs made the combined PSM / peptide FDR cutoff (${params.psmconflvl} / ${params.pepconflvl})" && exit 1)
  # SQLite lookup needs copying to not modify the input file which would mess up a rerun with -resume
  cat lookup > $psmlookup
  msslookup psms -i filtpep --dbfile $psmlookup ${params.onlypeptides ? '' : "--fasta ${td == 'target' ? "\"${tdb}\"" : "\"${ddb}\" --decoy"}"} ${params.martmap ? "--map ${martmap}" : ''}
  msspsmtable specdata -i filtpep --dbfile $psmlookup -o prepsms.txt --addmiscleav
  ${!params.noquant && td == 'target' ? "msspsmtable quant -i prepsms.txt -o qpsms.txt --dbfile $psmlookup --precursor ${params.isobaric ? '--isobaric' : ''}" : 'mv prepsms.txt qpsms.txt'}
  sed 's/\\#SpecFile/SpectraFile/' -i qpsms.txt
  ${!params.onlypeptides ? "msspsmtable genes -i qpsms.txt -o gpsms --dbfile $psmlookup" : ''}
  ${!params.onlypeptides ? "msslookup proteingroup -i qpsms.txt --dbfile $psmlookup" : ''}
  ${!params.onlypeptides ? "msspsmtable proteingroup -i gpsms -o ${params.hirief ? "pgpsms" : "$outpsms"} --dbfile $psmlookup" : "mv qpsms.txt ${params.hirief ? "pgpsms" : "$outpsms"}" }
  ${params.hirief ? "echo \'${groovy.json.JsonOutput.toJson(params.strips)}\' >> strip.json && peptide_pi_annotator.py -i $trainingpep -p pgpsms --o $outpsms --stripcolpattern Strip --pepcolpattern Peptide --fraccolpattern Fraction --stripdef strip.json --ignoremods \'*\'": ''} 
  msspsmtable split -i ${outpsms} --bioset
  ${setnames.collect() { "test -f '${it}.tsv' || echo 'No ${td} PSMs found for set ${it}' >> warnings" }.join(' && ') }
  """
}

// Collect setnames and merge with PSM tables for peptide table creation
setnames_psm
  .toList()
  .map { it -> [it.sort()]}
  .set { setlist_psm }
setpsmtables
  .map { it -> [it[0], it[1] instanceof java.util.List ? it[1] : [it[1]] ] }
  .map{ it -> [it[0], it[1].sort { a, b -> a.baseName.tokenize('.')[0] <=> b.baseName.tokenize('.')[0] }] } // names are setnames, sort on them then merge with sorted setnames
  .tap { deqms_psms }
  .merge(setlist_psm)
  .transpose()
  .set { psm_pep }


process psm2Peptides {

  input:
  set val(td), file('psms'), val(setname) from psm_pep
  
  output:
  set val(setname), val(td), file("${setname}_linmod") into pepslinmod
  set val(setname), val('peptides'), val(td), file("${setname}_linmod") into peptides_out
  set val(setname), val(td), file(psms), file('proteins'), val('proteins') into proteins
  set val(setname), val(td), file(psms), file('genes'), val('genes') into genes
  set val(setname), val(td), file(psms), file('symbols'), val('assoc') into symbols
  file('warnings') optional true into pepwarnings

  script:
  col = accolmap.peptides + 1  // psm2pep adds a column
  do_raw_isoquant = rawisoquant && td == 'target'
  """
  # Create peptide table from PSM table, picking best scoring unique peptides
  msspeptable psm2pep -i psms -o peptides --scorecolpattern svm --spectracol 1 ${!params.noquant && params.isobaric && td == 'target' ? "--isobquantcolpattern plex" : "" } ${!params.noquant && td == 'target' ? "--ms1quantcolpattern area" : ""}
  # Move peptide sequence to first column
  paste <( cut -f ${col} peptides) <( cut -f 1-${col-1},${col+1}-500 peptides) > peptide_table.txt
  # Create empty protein/gene/gene-symbol tables with only the identified accessions, will be filled later
  echo Protein ID|tee proteins genes symbols
  ${!params.onlypeptides ? "tail -n+2 psms|cut -f ${accolmap.proteins}|grep -v '\\;'| grep -v '^NA\$' | grep -v '^\$'|sort|uniq >> proteins" : "" }
  ${params.genes ? "tail -n+2 psms|cut -f ${accolmap.genes}|grep -v '\\;'| grep -v '^NA\$' | grep -v '^\$'|sort|uniq >> genes" : ""}
  ${params.symbols ? "tail -n+2 psms|cut -f ${accolmap.assoc}|grep -v '\\;'| grep -v '^NA\$' | grep -v '^\$'|sort|uniq >> symbols" : ""}
  ${do_raw_isoquant ? "msspsmtable isoratio -i psms -o pepisoquant --targettable peptide_table.txt --protcol ${accolmap.peptides} --isobquantcolpattern plex --minint 0.1 --denompatterns ${setdenoms[setname].join(' ')}" : ''}
  ${do_raw_isoquant ? "mv pepisoquant peptide_table.txt" : ''}
  # Create linear modeled q-values of peptides (modeled svm scores vs q-values) for more protein-FDR precision.
  msspeptable modelqvals -i peptide_table.txt -o ${setname}_linmod --scorecolpattern svm --fdrcolpattern '^q-value'
  cut -f \$(head -n1 "${setname}_linmod" | tr '\\t' '\\n' | grep -n 'linear modeled' | cut -f1 -d':') "${setname}_linmod" | grep "^[0-1]\\.[0-9]*\$" || echo 'Could not calculate linear modeled q-values for ${td} peptides of set ${setname}' >> warnings
  """
}


// Different amount of processes depending on genes and gene symbols are desired
// Input for proteins, genes and symbols is identical at this stage so tap and concat
// onto itself.
if (params.genes && params.symbols) { 
  pepslinmod
    .tap { pepsg; pepss }
    .concat(pepsg, pepss)
    .set { pepslinmod_prot }
  proteins
    .concat(genes, symbols)
    .join(pepslinmod_prot, by: [0,1])
    .set { prepgs_in }
} else if (params.genes) { 
  pepslinmod
    .tap { pepsg }
    .concat(pepsg)
    .set { pepslinmod_prot }
  proteins
    .concat(genes)
    .join(pepslinmod_prot, by: [0,1])
    .set { prepgs_in }
} else { 
  proteins
    .join(pepslinmod, by: [0,1])
    .set { prepgs_in }
}


/*
* Step 4: Infer and quantify proteins and genes
*/


// Group set/acctype T-D combinations and remove those with only target or only decoy
tprepgs_in = Channel.create()
dprepgs_in = Channel.create()
prepgs_in
  .groupTuple(by: [0,4])
  .filter { it -> it[1].size() == 2 }
  .transpose()
  .choice(tprepgs_in, dprepgs_in) { it[1] == 'target' ? 0 : 1 }


process proteinGeneSymbolTableFDR {
  
  when: !params.onlypeptides

  input:
  set val(setname), val(td), file('tpsms'), file('tproteins'), val(acctype), file('tpeplinmod') from tprepgs_in
  set val(setname), val(td), file('dpsms'), file('dproteins'), val(acctype), file('dpeplinmod') from dprepgs_in
  set file(tfasta), file(dfasta) from searchdbs

  output:
  set val(setname), val(acctype), file("${setname}_protfdr") into protfdrout
  file('warnings') optional true into fdrwarnings

  script:
  scorecolpat = acctype == 'proteins' ? '^q-value$' : 'linear model'
  """
  ${!params.noquant ? "mssprottable ms1quant -i tproteins -o tprotms1 --psmtable tpsms --protcol ${accolmap[acctype]}" : 'mv tproteins tprotms1'}
  ${rawisoquant ? "msspsmtable isoratio -i tpsms -o tprotquant --protcol ${accolmap[acctype]} --targettable tprotms1 --isobquantcolpattern plex --minint 0.1 --denompatterns ${setdenoms[setname].join(' ')}": 'mv tprotms1 tprotquant'}

  # score col is linearmodel_qval or q-value, but if the column only contains 0.0 or NA (no linear modeling possible due to only q<10e-04), we use svm instead
  tscol=\$(head -1 tpeplinmod | tr '\\t' '\\n' | grep -n "${scorecolpat}" | cut -f 1 -d':')
  dscol=\$(head -1 dpeplinmod | tr '\\t' '\\n' | grep -n "${scorecolpat}" | cut -f 1 -d':')
  if [ -n "\$(cut -f \$tscol tpeplinmod | tail -n+2 | egrep -v '(NA\$|0\\.0\$)')" ] && [ -n "\$(cut -f \$dscol dpeplinmod | tail -n+2 | egrep -v '(NA\$|0\\.0\$)')" ]
    then
      scpat="${scorecolpat}"
      logflag="--logscore"
    else
      scpat="svm"
      logflag=""
      echo 'Not enough q-values or linear-model q-values for peptides to calculate FDR for ${acctype} of set ${setname}, using svm score instead.' >> warnings
  fi

  mssprottable bestpeptide -i tprotquant -o tbestpeptides --peptable tpeplinmod --scorecolpattern "\$scpat" \$logflag --protcol ${accolmap[acctype] + 1}
  mssprottable bestpeptide -i dproteins -o dbestpeptides --peptable dpeplinmod --scorecolpattern "\$scpat" \$logflag --protcol ${accolmap[acctype] + 1}

  mssprottable ${acctype == 'proteins' ? 'protfdr' : 'pickedfdr'} -i tbestpeptides --decoyfn dbestpeptides -o ${setname}_protfdr ${acctype == 'genes' ? "--picktype fasta --targetfasta '$tfasta' --decoyfasta '$dfasta' ${params.fastadelim ? "--fastadelim '${params.fastadelim}' --genefield '${params.genefield}'" : '' }" : ''} ${acctype == 'assoc' ? '--picktype result' : ''}
  """
}


psmwarnings
  .concat(pepwarnings)
  .concat(fdrwarnings)
  .toList()
  .set { warnings }

// setname, acctype, outfile
peptides_out
  .filter { it[2] == 'target' }
  .map { it -> [it[0], it[1], it[3]] }
  .concat(protfdrout)
  .set { features_out }

deqms_psms
  .filter { it[0] == 'target' }
  .map { it -> [it[1].collect { x -> x.baseName.tokenize('.')[0] }, it[1]] }
  .transpose()
  .cross(features_out)
  .map { it -> it[0] + [it[1][1], it[1][2]] }
  .into { feats_out; deqms }


process normalizeFeaturesDEqMS {
  input:
  set val(setname), file("psms"), val(acctype), file("features") from deqms
  output:
  set val("${setname}_${acctype}"), file("${setname}_feats"), file("psmcounts") into quanted_feats
  when: normalize
  
  script:
  """
  # get col nrs for isobaric quant values and create new PSM table with only those and feature columns
  # need to use awk first because cut cannot paste peptide column twice (happens when peptide is acctype)
  awk -F \$'\\t' -v OFS=\$'\\t' '{print \$${accolmap.peptides}, \$${accolmap[acctype]}}' psms ${acctype == 'peptides' ? '| sed \'s/Peptide/Accession/\' ' : ''} > pepacc
  channelcols=\$(head -n1 psms | tr '\\t' '\\n' | grep -n plex | cut -f1 -d ':'| tr '\\n' ',' | sed 's/,\$//')
  paste pepacc <(cut -f "\$channelcols" psms) > psmvals
  # run deqMS normalization and summarization, which produces logged ratios

  ${params.denoms ? "denomcols=\$(egrep -n \'(${setdenoms[setname].join('|')})\' <( head -n1 psmvals | tr '\\t' '\\n') | cut -f1 -d ':' | tr '\\n' ',' | sed 's/,\$//') " : ""}
  deqms_normalize.R psmvals features $setname ${params.denoms ? "\$denomcols" : ''}
  # join feat tables on normalized proteins
  paste <(head -n1 features) <(head -n1 normalized_feats | cut -f2-2000) <(echo PSM counts) > ${setname}_feats 
  join -a1 -o auto -e 'NA' -t \$'\\t' <(tail -n+2 features | sort -k1b,1 ) <(tail -n+2 normalized_feats | sort -k1b,1) >> feats_quants
  join -a1 -o auto -e 'NA' -t \$'\\t' feats_quants <(sort -k1b,1 psmcounts) >> ${setname}_feats
  """
}


if(normalize) {
  feats_out
    .map { it ->  ["${it[0]}_${it[2]}".toString()] + it }
    .join(quanted_feats)
    .groupTuple(by: 3)  // all outputs of same accession type together.
    .map { it -> [it[1], it[3], it[5], it[6]] }
    .set { ptables_to_merge }
} else {
  feats_out
    .map { it -> [it[0], it[2], it[3]] } // setname, accession type, feature table
    .groupTuple(by: 1) // collect all tables for same feature
    .set { ptables_to_merge }
}

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
  set val(setnames), val(acctype), file(tables), file("psmcounts?") from ptables_to_merge
  file(lookup) from tlookup
  file('sampletable') from Channel.from(sampletable).first()
  
  output:
  set val(acctype), file('proteintable'), file('sampletable') into featqc_extra_peptide_samples
  set val(acctype), file('proteintable') into merged_feats

  script:
  if (normalize)
  """
  # SQLite lookup needs copying to not modify the input file which would mess up a rerun with -resume
  cat $lookup > db.sqlite
  msslookup ${acctype == 'peptides' ? 'peptides --fdrcolpattern \'^q-value\' --peptidecol' : 'proteins --fdrcolpattern \'q-value\' --protcol'} 1 --dbfile db.sqlite -i ${tables.join(' ')} --setnames ${setnames.join(' ')} ${!params.noquant ? "--ms1quantcolpattern area" : ""}  ${!params.noquant && params.isobaric ? '--isobquantcolpattern plex' : ''} ${acctype in ['genes', 'assoc'] ? "--genecentric ${acctype}" : ''}
  ${acctype == 'peptides' ? 'msspeptable build' : 'mssprottable build --mergecutoff 0.01'} --dbfile db.sqlite -o mergedtable ${!params.noquant && params.isobaric ? '--isobaric' : ''} ${!params.noquant ? "--precursor": ""} --fdr ${acctype in ['genes', 'assoc'] ? "--genecentric ${acctype}" : ''} ${params.onlypeptides ? "--noncentric" : ''}
  # join psm count tables, first make a header from setnames
  head -n1 mergedtable > tmpheader

  # exchange sample names in header
  ${params.sampletable && params.isobaric ?  
    'sed -i  "s/[^A-Za-z0-9_\\t]/_/g" sampletable ; \
    while read line ; do read -a arr <<< $line ; sed -i "s/${arr[1]}_\\([a-z0-9]*plex\\)_${arr[0]}/${arr[3]}_${arr[2]}_${arr[1]}_\\1_${arr[0]}/" tmpheader ; done < sampletable' \
  : ''}
  # Add psm quant nr field to header
  for setn in ${setnames.join(' ')}; do echo "\$setn"_quanted_psm_count ; done >> tmpheader
  tr '\\n' '\\t' < tmpheader | sed 's/\\s\$/\\n/;s/\\#/Amount/g' > header  # sed to sub trailing tab for a newline, and not have pound sign
  # then join the table content
  tail -n+2 mergedtable | sort -k1b,1 > joined
  for count in \$(seq 1 ${setnames.toList().size}); do join -a1 -o auto -e 'NA' -t \$'\\t' joined <(sort -k1b,1 psmcounts"\$count" ) >> joined_tmp; mv joined_tmp joined; done
  # finally put header on content
  cat header joined > proteintable
  """
  else
  """
  cat $lookup > db.sqlite
  msslookup ${acctype == 'peptides' ? 'peptides --fdrcolpattern \'^q-value\' --peptidecol' : 'proteins --fdrcolpattern \'q-value\' --protcol'} 1 --dbfile db.sqlite -i ${tables.join(' ')} --setnames ${setnames.join(' ')} ${!params.noquant ? "--ms1quantcolpattern area" : ""}  ${!params.noquant && params.isobaric ? '--isobquantcolpattern plex' : ''} ${acctype in ['genes', 'assoc'] ? "--genecentric ${acctype}" : ''}
  ${acctype == 'peptides' ? 'msspeptable build' : 'mssprottable build --mergecutoff 0.01'} --dbfile db.sqlite -o mergedtable ${!params.noquant && params.isobaric ? '--isobaric' : ''} ${!params.noquant ? "--precursor": ""} --fdr ${acctype in ['genes', 'assoc'] ? "--genecentric ${acctype}" : ''} ${params.onlypeptides ? "--noncentric" : ''}
  ${!params.noquant && params.isobaric ? "sed -i 's/\\ \\-\\ \\#\\ quanted\\ PSMs/_quanted_psm_count/g' mergedtable": ''}
  sed -i 's/\\#/Amount/g' mergedtable
  # exchange sample names in header
  head -n1 mergedtable > tmpheader
  ${params.sampletable && params.isobaric ?  
    'sed -i  "s/[^A-Za-z0-9_\\t]/_/g" sampletable ; \
    while read line ; do read -a arr <<< $line ; sed -i "s/${arr[1]}_\\([a-z0-9]*plex\\)_${arr[0]}/${arr[3]}_${arr[2]}_${arr[1]}_\\1_${arr[0]}/" tmpheader ; done < sampletable'\
  : ''}
  cat tmpheader <(tail -n+2 mergedtable) > proteintable
  """
}

(plain_feats, dqms_feats) = ( params.deqms ? [Channel.empty(), merged_feats] : [merged_feats, Channel.empty()])

process calculateDEqMS {

  input:
  set val(acctype), file('feats') from dqms_feats 
  file('sampletable') from Channel.from(sampletable).first()

  output:
  set val(acctype), file('deqms_output') into dqms_out 

  script:
  """
  deqms.R 
  """
}



psm_result
  .filter { it[0] == 'target' }
  .merge(scans_result)
  .map { it -> [it[0], it[1], it[2], it[3].unique()] }
  .set { targetpsm_result }


process psmQC {
  input:
  set val(td), file('psms'), file('scans'), val(plates) from targetpsm_result
  val(setnames) from setnames_psmqc
  output:
  set val('psms'), file('psmqc.html'), file('summary.txt') into psmqccollect
  val(plates) into qcplates
  // TODO no proteins == no coverage for pep centric
  script:
  """
  qc_psms.R ${setnames[0].size()} ${fractionation ? 'TRUE' : 'FALSE'} ${plates.join(' ')}
  echo "<html><body>" > psmqc.html
  for graph in psm-scans missing-tmt miscleav
    do
    [[ -e \$graph ]] && paste -d \\\\0  <(echo "<div class=\\"chunk\\" id=\\"\${graph}\\"><img src=\\"data:image/png;base64,") <(base64 -w 0 \$graph) <(echo '"></div>') >> psmqc.html
    done 
  for graph in retentiontime precerror fryield msgfscore
    do
    for plateid in ${plates.join(' ')}
      do
      plate="PLATE___\${plateid}___\${graph}"
      [[ -e \$plate ]] && paste -d \\\\0  <(echo "<div class=\\"chunk \$plateid\\" id=\\"\${graph}\\"><img src=\\"data:image/png;base64,") <(base64 -w 0 \$plate) <(echo '"></div>') >> psmqc.html
      done 
    done
  echo "</body></html>" >> psmqc.html
  """
}

featqc_extra_peptide_samples
  .filter { it[0] == 'peptides' }
  .map { it -> [it[1], it[2]] }
  .set { featqc_peptides_samples }

plain_feats
  .mix(dqms_out)
  .merge(setnames_featqc)
  .combine(featqc_peptides_samples)
  .set { featqcinput }


process featQC {
  publishDir "${params.outdir}", mode: 'copy', overwrite: true, saveAs: {it == "feats" ? "${outname}_table.txt": null}

  input:
  set val(acctype), file('feats'), val(setnames), file(peptable), file(sampletable) from featqcinput
  output:
  file('feats') into featsout
  set val(acctype), file('featqc.html'), file('summary.txt'), file('overlap') into qccollect

  script:
  outname = (acctype == 'assoc') ? 'symbols' : acctype
  """
  # Create QC plots and put them base64 into HTML, R also creates summary.txt
  qc_protein.R --sets ${setnames.collect() { "'$it'" }.join(' ')} --feattype ${acctype} --peptable $peptable ${params.sampletable ? "--sampletable $sampletable" : ''}
  echo "<html><body>" > featqc.html
  for graph in featyield precursorarea coverage isobaric nrpsms nrpsmsoverlapping percentage_onepsm ms1nrpeps;
    do
    [ -e \$graph ] && paste -d \\\\0  <(echo "<div class=\\"chunk\\" id=\\"\${graph}\\"><img src=\\"data:image/png;base64,") <(base64 -w 0 \$graph) <(echo '"></div>') >> featqc.html
    done 
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
    paste -d \\\\0  <(echo '<div><img src="data:image/png;base64,') <(base64 -w 0 \$graph) <(echo '"></div>') >> featqc.html
    done
    [ -e pca ] && echo '</div>' >> featqc.html

  echo "</body></html>" >> featqc.html
  ${acctype == 'peptides' ? 'touch summary.txt' : ''}

  # Create overlap table
  qcols=\$(head -n1 feats |tr '\\t' '\\n'|grep -n "_q-value"| tee nrsets | cut -f 1 -d ':' |tr '\\n' ',' | sed 's/\\,\$//')
  protcol=\$(head -n1 feats | tr '\\t' '\\n' | grep -n Protein | cut -f1 -d ':')
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
  """
}

qccollect
  .concat(psmqccollect)
  .toList()
  .map { it -> [it.collect() { it[0] }, it.collect() { it[1] }, it.collect() { it[2] }, it.collect() { it[3] }] }
  .set { collected_feats_qc }


process collectQC {

  publishDir "${params.outdir}", mode: 'copy', overwrite: true

  input:
  set val(acctypes), file('feat?'), file('summary?'), file('overlap?') from collected_feats_qc
  val(plates) from qcplates
  file('sw_ver') from software_versions_qc
  file('warnings??') from warnings

  output:
  set file('qc_light.html'), file('qc_full.html')

  script:
  """
  count=1; for ac in ${acctypes.join(' ')}; do mv feat\$count \$ac.html; mv summary\$count \${ac}_summary; mv overlap\$count \${ac}_overlap; ((count++)); done
  join -j 1 -o auto -t '\t' <(head -n1 psms_summary) <(head -n1 peptides_summary) > psmpepsum_header
  join -j 1 -o auto -t '\t' <(tail -n+2 psms_summary | sort -k1b,1 ) <(tail -n+2 peptides_summary | sort -k1b,1 ) > psmpepsum_tab

  # onlypeptides makes a quick summary, else also add proteins
  ${params.onlypeptides ? 'cat psmpepsum_header psmpepsum_tab | tee summary pre_summary_light_tab' : 'join -j 1 -o auto -t \'\t\' psmpepsum_tab <(sort -k1b,1 <(tail -n+2 proteins_summary)) > pepprotsum_tab && join -j 1 -o auto -t \'\t\' psmpepsum_header <(head -n1 proteins_summary) > pepprotsum_head'}
  ${params.onlypeptides ? "awk -v FS='\\t' -v OFS='\\t' '{print \$1,\$3,\$2}' pre_summary_light_tab > summary_light" : ""}

  # in case of genes, join those on the prot/pep tables (full summary) and psmpeptables (light summary), else passthrough those to summaries
  ${params.genes ?  'join -j 1 -o auto -t \'\t\' pepprotsum_tab <( sort -k1b,1 <( tail -n+2 genes_summary)) > summary_tab && join -j 1 -o auto -t \'\t\' pepprotsum_head <(head -n1 genes_summary) > summary_head && cat summary_head summary_tab > summary' : "${!params.onlypeptides ? 'cat pepprotsum_head pepprotsum_tab | tee summary summary_light' : ""}"}
  ${params.genes ?  'join -j 1 -o auto -t \'\t\' psmpepsum_tab <( sort -k1b,1 <(tail -n+2 genes_summary)) > summary_light_tab' : ''}
  ${params.genes ?  'join -j 1 -o auto -t \'\t\' psmpepsum_header <( head -n1 genes_summary) > summary_light_head && cat summary_light_head summary_light_tab > summary_light' : ''}

  # remove Yaml from software_versions to get HTML
  grep -A \$(wc -l sw_ver | cut -f 1 -d ' ') "data\\:" sw_ver | tail -n+2 > sw_ver_cut
  
  # merge warnings
  test -f warnings* && cat warnings* > warnings.txt
  # collect and generate HTML report
  qc_collect.py $baseDir/assets/qc_full.html $params.name ${fractionation ? "frac" : "nofrac"} ${plates.join(' ')}
  qc_collect.py $baseDir/assets/qc_light.html $params.name ${fractionation ? "frac" : "nofrac"} ${plates.join(' ')}
  """
}


/* 
 * STEP 3 - Output Description HTML
*/
process output_documentation {
    tag "$prefix"

    publishDir "${params.outdir}/Documentation", mode: 'copy'

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
    def subject = "[nf-core/ddamsproteomics] Successful: $workflow.runName"
    if(!workflow.success){
      subject = "[nf-core/ddamsproteomics] FAILED: $workflow.runName"
    }
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
          log.info "[nf-core/ddamsproteomics] Sent summary e-mail to $params.email (sendmail)"
        } catch (all) {
          // Catch failures and try with plaintext
          [ 'mail', '-s', subject, params.email ].execute() << email_txt
          log.info "[nf-core/ddamsproteomics] Sent summary e-mail to $params.email (mail)"
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

    log.info "[nf-core/ddamsproteomics] Pipeline Complete"

}
