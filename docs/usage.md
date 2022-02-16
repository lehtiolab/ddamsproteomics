# lehtiolab/ddamsproteomics: Usage

## Table of contents

* [Introduction](#general-nextflow-info)
* [Running the pipeline](#running-the-pipeline)
* [Updating the pipeline](#updating-the-pipeline)
* [Reproducibility](#reproducibility)
* [Main arguments](#main-arguments)
    * [`-profile`](#-profile-single-dash)
        * [`docker`](#docker)
        * [`awsbatch`](#awsbatch)
        * [`standard`](#standard)
        * [`none`](#none)
* [Job Resources](#job-resources)
* [Automatic resubmission](#automatic-resubmission)
* [Custom resource requests](#custom-resource-requests)
* [AWS batch specific parameters](#aws-batch-specific-parameters)
    * [`-awsbatch`](#-awsbatch)
    * [`--awsqueue`](#--awsqueue)
    * [`--awsregion`](#--awsregion)
* [Other command line parameters](#other-command-line-parameters)
    * [`--outdir`](#--outdir)
    * [`--email`](#--email)
    * [`-name`](#-name-single-dash)
    * [`-resume`](#-resume-single-dash)
    * [`-c`](#-c-single-dash)
    * [`--max_memory`](#--max_memory)
    * [`--max_time`](#--max_time)
    * [`--max_cpus`](#--max_cpus)
    * [`--plaintext_emails`](#--plaintext_emails)


## General Nextflow info
Nextflow handles job submissions on SLURM or other environments, and supervises running the jobs. Thus the Nextflow process must run until the pipeline is finished. We recommend that you put the process running in the background through `screen` / `tmux` or similar tool. Alternatively you can run nextflow within a cluster job submitted your job scheduler.

It is recommended to limit the Nextflow Java virtual machines memory. We recommend adding the following line to your environment (typically in `~/.bashrc` or `~./bash_profile`):

```bash
NXF_OPTS='-Xms1g -Xmx4g'
```

## Running the pipeline
The typical command for running the pipeline is as follows:
```bash
nextflow run lehtiolab/ddamsproteomics --mzmls '/path/to/*.mzML' --tdb /path/to/proteins.fa --mods 'oxidation;carbamidomethylation' -profile standard,docker
```

This will launch the pipeline with the `docker` configuration profile. See below for more information about profiles.

Note that the pipeline will create the following files in your working directory:

```bash
work            # Directory containing the nextflow working files
results         # Finished results (configurable, see below)
.nextflow_log   # Log file from Nextflow
# Other nextflow hidden files, eg. history of pipeline runs and old logs.
```

### Updating the pipeline
When you run the above command, Nextflow automatically pulls the pipeline code from GitHub and stores it as a cached version. When running the pipeline after this, it will always use the cached version if available - even if the pipeline has been updated since. To make sure that you're running the latest version of the pipeline, make sure that you regularly update the cached version of the pipeline:

```bash
nextflow pull lehtiolab/ddamsproteomics
```

### Reproducibility
It's a good idea to specify a pipeline version when running the pipeline on your data. This ensures that a specific version of the pipeline code and software are used when you run your pipeline. If you keep using the same tag, you'll be running the same version of the pipeline, even if there have been changes to the code since.

First, go to the [lehtiolab/ddamsproteomics releases page](https://github.com/lehtiolab/ddamsproteomics/releases) and find the latest version number - numeric only (eg. `1.3.1`). Then specify this when running the pipeline with `-r` (one hyphen) - eg. `-r 1.3.1`.

This version number will be logged in reports when you run the pipeline, so that you'll know what you used when you look back in the future.


## Main Arguments

### `-profile`
Use this parameter to choose a configuration profile. Profiles can give configuration presets for different compute environments. Note that multiple profiles can be loaded, for example: `-profile standard,docker` - the order of arguments is important!

* `standard`
    * The default profile, used if `-profile` is not specified at all.
    * Runs locally and expects all software to be installed and available on the `PATH`.
* `docker`
    * A generic configuration profile to be used with [Docker](http://docker.com/)
    * Pulls software from dockerhub: [`lehtiolab/ddamsproteomics`](http://hub.docker.com/r/lehtiolab/ddamsproteomics/)
* `singularity`
    * A generic configuration profile to be used with [Singularity](http://singularity.lbl.gov/)
    * Pulls software from singularity-hub
* `conda`
    * A generic configuration profile to be used with [conda](https://conda.io/docs/)
    * Pulls most software from [Bioconda](https://bioconda.github.io/)
* `awsbatch`
    * A generic configuration profile to be used with AWS Batch.
* `test`
    * A profile with a complete configuration for automated testing
    * Includes links to test data so needs no other parameters
* `none`
    * No configuration at all. Useful if you want to build your own config from scratch and want to avoid loading in the default `base` config profile (not recommended).

### `--mzmls`
Use this to specify the location of your input mzML files. For example:

```bash
--mzML 'path/to/data/sample_*.mzML'
```
The path must be enclosed in quotes when using wildcards like `*`


### `--mzmldef`
Alternative to the above --mzml this would pass a mzML definition (txt) file which contains the mzML specifications. This also enables runs with specific fractionation such as HiRIEF or high pH, and the specification of individual instruments per file.

```bash
--mzmldef /path/to/data/mzmls.txt
```

The file itself is tab-separated without header, contains a single line per mzML file specified as follows:
`/path/to/file	instrument_type	sample_or_sampleset_name	OPTIONAL:fractionation_plate_name	OPTIONAL:fraction_nr`
Fractionation is automatically detected from this file, and enforced if ANY of the files have a fraction. This mainly has implications for QC though, identification and quantification are not much impacted by specifying fractionation. Instrument type can currently be one of 'qe', 'qehf', 'velos', 'lumos', 'qehfx', 'timstof', or 'lowres'.
Examples of instruments can be found in [this MSGF+ parameter file](https://github.com/MSGFPlus/msgfplus/blob/master/docs/ParameterFiles/MSGFPlus_Tryp_NoMods_20ppmParTol.txt).

### `--tdb`
Target database. Decoy databases are created "tryptic-reverse" by the pipeline and searches are against a
concatenated database (T-TDC)

```bash
--tdb /path/to/Homo_sapiens.pep.all.fa
```


### `--mods`, `--locptms`, `--ptms`
Modifications as in UNIMOD, although only a selected number are available by name. You can extend this list
by adding entries to `assets/msgfmods.txt`. `--ptms` and `--locptms` are for stable/labile PTMs respectively,
and they can optionally get isobaric quantification normalization (below).

```bash
--mods "Carbamidomethylation;Oxidation" --ptms Acetyl --locptms Phoshpo
```


### Output types
The pipeline will produce by default PSM, peptide, and protein tables. You may pass FASTA databases that contain mixtures of ENSEMBL, Uniprot, or other types of entries. Use `--genes` and `--ensg` to output a gene(name)-centric table and an ENSG-centric table. If you rather have less output, use `--onlypeptides` to not output a protein table.


### Quantitation
Isobaric data can be specified as such `--isobaric 'set1:tmt10plex:127N:128N set2:tmtpro:sweep set3:itraq8plex:intensity'`. Here PSMs will be quantified for 3 different isobaric sample sets in a somewhat contrived example with different chemistries. Isobaric quantitation is done using OpenMS IsobaricAnalyzer and will also output the precursor purity (fraction of precursor intensity in the selection window) to the PSM table. A filter can be used to set a minimum purity for PSM isobaric quant to not be set to NA, using e.g. `--minprecursorpurity 0.3`, default is not to filter. The resulting values will be summarized from PSM to peptide/protein/gene by taking median PSM values per feature. Prior to this, they can be log2-transformed and normalized to e.g. an internal standard using denominator channels as above in `set1`, or median sweep (i.e. use median PSM value as denominator for each PSM in `set2`) to generate log2(ratios). A possibility shown in `set3` is also to output non-normalized median PSM intensity per protein. When result tables have been summarized these can be median-centered, which is passed using `--normalize`. By default, PSMs with an NA value in any channel will not be used in summarizing isobaric quantification data. If you want to use these (possibly more noisy) PSMs, e.g. when having empty channels, you can pass `--keepnapsmsquant`.

MS1 quantitation is done using Dinosaur and its features are aligned to PSMs using msstitch, using summed intensity of a feature. To not output any MS1 or isobaric data, use `--noquant`. If Dinosaur for some reason doesn't work, you can use Hardklor/Kronik, by specifying `--hardklor`

### Differential expression analysis
DEqMS is used for DE analysis using `--deqms` and it needs to know your sample group names. For this, you can pass a TSV file with sample names to `--sampletable`, it should contain a line for each channel/set combination with channel, set, sample, sample group e.g.:

```
126    setA   DMSO1  CTRL
127N   setA   ABC1   TREAT 
127C   setA   DMSO2  CTRL 
128N   setA   ABC2   TREAT
129N   setA   pool
...
```

For DE analysis, sample-channels that e.g represent internal standards will be filtered out if no sample group is given, see above 129N channel.
N.B. Even when not using DEqMS you can provide a sample table for annotation of your quant output.


### PTM analysis
As mentioned, labile PTMS reported by the search engine will be scored using Luciphor2, which will output the best scoring PTM localization and a false localization rate. Note that this is only beneficial for labile PTMs. Aside from that if any high-scoring PTMs are found by luciphor the pipeline will report these as well. All of this will end up in a separate PTM PSM table and a PTM peptide table.
When passing `--totalproteomepsms`, the isobaric quant ratios for matching genes from a global search (i.e. no modifications) will be subtracted from the PTM peptide table quant. If `--onlypeptides` is used, quant from proteins will be used as a denominator.

For normalizing PTM tables, `--normalize` can be used for median-centering. Since PTM tables can be somewhat small and possibly skewed in their quantitation, a separate protein table is prepared from the PSMs in `--totalproteomepsms`, to get the channel median normalization factors from.


### Reusing data
If you have finished a rather large analysis and wish to rerun a part of it or add more fractions, due to e.g. new MS data, you may do so by passing

```
  --targetpsms oldpsmtable.txt --decoypsms old_decoy_psms.txt \
  --targetpsmlookup old_target_psmlookup.sql \
  --decoypsmlookup old_decoy_lookup.sqlite \
  --ptmpsms old_ptm_psmtable.txt # Optional of course
```
Now you can run a single sample set and combine the output with the previous run, which if it has an identical set name, will be cleaned first (the set will be removed from old output). Make sure to use the same parameters to get the same result for the old data, which will be regenerated from the PSM table.


## Job Resources
### Automatic resubmission
Each step in the pipeline has a default set of requirements for number of CPUs, memory and time. For most of the steps in the pipeline, if the job exits with an error code of `143` (exceeded requested resources) it will automatically resubmit with higher requests (2 x original, then 3 x original). If it still fails after three times then the pipeline is stopped.

### Custom resource requests
Wherever process-specific requirements are set in the pipeline, the default value can be changed by creating a custom config file. See the files in [`conf`](../conf) for examples.

## AWS Batch specific parameters
Running the pipeline on AWS Batch requires a couple of specific parameters to be set according to your AWS Batch configuration. Please use the `-awsbatch` profile and then specify all of the following parameters.
### `--awsqueue`
The JobQueue that you intend to use on AWS Batch.
### `--awsregion`
The AWS region to run your job in. Default is set to `eu-west-1` but can be adjusted to your needs.

Please make sure to also set the `-w/--work-dir` and `--outdir` parameters to a S3 storage bucket of your choice - you'll get an error message notifying you if you didn't.

## Other command line parameters

### `--outdir`
The output directory where the results will be saved.

### `--email`
Set this parameter to your e-mail address to get a summary e-mail with details of the run sent to you when the workflow exits. If set in your user config file (`~/.nextflow/config`) then you don't need to speicfy this on the command line for every run.

### `-name`
Name for the pipeline run. If not specified, Nextflow will automatically generate a random mnemonic.

This is used in the MultiQC report (if not default) and in the summary HTML / e-mail (always).

**NB:** Single hyphen (core Nextflow option)

### `-resume`
Specify this when restarting a pipeline. Nextflow will used cached results from any pipeline steps where the inputs are the same, continuing from where it got to previously.

You can also supply a run name to resume a specific run: `-resume [run-name]`. Use the `nextflow log` command to show previous run names.

**NB:** Single hyphen (core Nextflow option)

### `-c`
Specify the path to a specific config file (this is a core NextFlow command).

**NB:** Single hyphen (core Nextflow option)

Note - you can use this to override defaults. For example, you can specify a config file using `-c` that contains the following:

```nextflow
process.$multiqc.module = []
```

### `--max_memory`
Use to set a top-limit for the default memory requirement for each process.
Should be a string in the format integer-unit. eg. `--max_memory '8.GB'``

### `--max_time`
Use to set a top-limit for the default time requirement for each process.
Should be a string in the format integer-unit. eg. `--max_time '2.h'`

### `--max_cpus`
Use to set a top-limit for the default CPU requirement for each process.
Should be a string in the format integer-unit. eg. `--max_cpus 1`

### `--plaintext_email`
Set to receive plain-text e-mails instead of HTML formatted.
