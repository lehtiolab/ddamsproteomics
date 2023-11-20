# lehtiolab/ddamsproteomics: Changelog
## Version 2.17 [2023-11-20]
- Percolator forces concatenated T/D search instead of mix-max
- Using mixed isobaric types possible
- Catching some more warnings, empty results
- Tests included
- Add peptide-is-in external fasta columns
- External config versions


## Version 2.16 [2023-08-30]
- Fixed Luciphor algorithm selection to use HCD when activation type is auto/any, and not use CID then
- Fix mergePTM process file quoting when only one set in it


## Version 2.15 [2023-08-25]
- Revert to --ignore-target-hits for decoy generation, do not remove decoy tryptic peptides that match target


## Version 2.14 [2023-08-22]
- Fixed off-by-one error in PTM site nr reporting on proteins (#17)


## Version 2.13 [2023-07-14]
- Fixed PTMs/TMT18 bug which crashed pipeline


## Version 2.12 [2023-06-27]
- Add support for TMT18 (OpenMS pre-release)
- Adapt to newer nextflow, overwrite trace files in lehtio conf
- Add support for isobaric precursor purity in TIMS, by using newer pwiz for centroiding
- QC plot and table fixes
- Fix throw out decoys that are un-shufflable (10x always target)
- Doc update for PTM analysis


## Version 2.11 [2022-09-16]
- Bugfix in MS1 quant where too large retention time tolerance was used
- MS1 quant RT/mz tolerances parametrized


## Version 2.10 [2022-09-02]
- Bugfix to allow setnames with spaces etc (quoting)
- Dynamic memory allocation for MSGF
- Error message for often occurring DEqMS with only one channel in a sample group


## Version 2.9 [2022-06-02]
- Bugfix to allow for user-formatted modification specifications again


## Version 2.8 [2022-05-22]
- Total proteome normalization had no keep-na-psms enabled when that was desired, fixed
- Produce PEP column in merged peptide table (was only in PSM table earlier)
- Removed word "Amount" from output and changed to "count" or "total number of"
- QC changes in plotting, added summary and overlap tables for PTMs, format change in report
- Flanking sequences for PTMs added to PTM peptide tables (7aa + PTM-aa + 7aa, for each PTM of a peptide)
- Fix bug where no mixed PTM peptides were in output (e.g. peptide with GG and Phospho)
- Crash detection and warning when not enough PSMs for luciphor, and handling no-PTM situations


## Version 2.7 [2022-02-23]
- Median-center normalizing of PTM tables with total proteome table enabled
- No recalculating FDR from percolator, and output PEP in PSMs/peptides
- Small bug fixes (quoting etc)
- Pipeline does not crash when e.g. a single set has no decoy PSMs, will warn in QC, crash when no target PSMs at all are found
- PTM peptide tables can be outputted as both totalproteome-adjusted and non-adjusted


## Version 2.6 [2021-10-28]
- Sort setnames on mergeing output
- Fixed bug which crashed using multi-set and stabile PTMs
- Remove X__POOL identifier from sampletable and made it an internal - no output of it either
- First output all quant channels from all sets, THEN the #/peptide per channel, not interlace them
- Fix bug where memory quant does not respect the max value, which crashed if mzML were large and exceeding system available memory
- Add support for TIMSTOF data
- Autodetect fractionation from mzml definition, not --fractions
- Added normalization factors as text table for easy copying in QC report
- Possible to rerun from PSM table inputs without new mzMLs, i.e. change samples, output settings etc
- QC fixes for reruns and partial reruns (when adding files to set)


## Version 2.5 [2021-03-13]
- Bugfix in msstitch 3.7 for peptide tables in large datasets

## Version 2.4 [2021-03-05]
- Fix rerunning (complementing run) bugs
- Run Luciphor on more PSMs (cutoff 0.2 FDR, not all as otherwise it takes a long time), as it should learn which are bad.
- PTM peptide table isobaric quant divided/subtracted by/with global proteome (ratio PTM pep / ratio gene or prot)
- Do not crash when trying to merge peptide tables and a peptide has e.g. >600 zinc finger protein matches
- Report protein PTM location to PTM table
- PTMs can now be stable with --ptms or labile with --locptms. Labile PTMs get localization analysis,
  both can get quant normalization.
- Added PTM QC: numbers of PTM peptides, site distributions
- Added HFX and Lumos options for MSGF instrument types (for easy to use: same MSGF instruments as earlier though)
- Allow (but substitute to underscore) some more special characters in sample names, set names
- QC missed cleavage percentage is now of set, not of all PSMs, and show fraction yield also for empty fractions
- Add precursor purity column to PSM table and QC, make it possible to filter bad purity --min-precursor-purity
- Check FASTA for duplicates/general validity before running search, abort when not valid FASTA
- Bugfix: modifications such as Glygly do not block isobaric labels

## Version 2.3 [2020-09-25]
- Possible to re-run analysis with a new sample set only, then add to existing data
- More isobaric quant summarizing/normalization options by moving to msstitch 3.5
- Fix bug where all sets' "number of (unique) PSMs" were the same (msstitch 3.5)
- Consistent protein groups, bugfixed MS1 aligning in msstitch 3.5
- SVG output in QC, makes QC HTML lighter and also does not crash the QC by exceeding png size limit on many sets
- Fix Luciphor bug where residues were annotated wrong (off by one residue number)
- More documentation

## Version 2.2 [2020-07-24]
- Fix MS1 quant bug by upgrading to msstitch 3.3
- Switch MS1 quantification to Dinosaur from Hardklor/Kronik to include showing FWHM
- Fix name collision bug where multi-set input led to multiple files from luciphor all called ptms.txt, cannot stage those

## Version 2.1 [2020-07-22]
- Bugfix that crashed pipeline, likely merge problem, should have tested better

## Version 2.0 [2020-07-20]
- Interface changes:
  - modification specification, no more mod file
  - `--genes`/`--symbols` are now called `--ensg`/`--genes`
  - multi-instrument runs, mzml definition file can be used to spec instrument
  - Enable diff TMT/itraq mixing
  - Multiple DBs possible to pass e.g. `--tdb /path/to/*.fa`
  - `--hirief` flag contains the peptide pI file previously specified in `--pipep`
- msstitch updated to v3.2 (no longer need Biomart map)
- General package update for openMS, R packages
- Luciphor based PTM reporting on PSM/peptide tables, for phospho, e.g not working yet is acetyl/TMT 


## Version 1.5 [2020-03-17]
- Upgrade percolator to 3.4
- Upgrade openMS to 2.5
- Upgrade MSGF+ to 2020.04.14
- TMTpro 16plex can be searched
- Channel median normalization factors output to QC
- PCA in QC also for non-DEqMS samples

## Version 1.4 [2019-12-16]
- Fix bug in ENSG FDR reporting (upgrade msstitch 2.19)
- Low impact fixes for runs with bad MS / few PSMs.


## Version 1.3 [2019-11-29]
- Do not crash when PSM data for some but not all sets is not good enough
- miscleav/set/fr/plate not by python script but in pipeline or by msstitch
- clearer code
- percolator do not use -y (mixmax)
- Do not filter decoy DB of target sequences until certainty about how filter should behave
- Warnings in QC output


## Version 1.2 [2019-11-20]
- Updated to MSGF version 2019.07.03
- Added MSGF maxMissedCleavages parameter
- Moved release repo to lehtiolab/ddamsproteomics
- Removed decoy/trypsinization script to use msstitch tool instead
- QC bugfix: dont crash when there are no peptides with 1missed cleavage
- Got a Zenodo DOI

## Version 1.1
- More QC for DEqMS, added  PCA and Volcano plots
- Added support for non-linear pI strips for HiRIEF (e.g. 3.4-4.8), defined as multiple fraction intervals
- Fixed bug that did not filter out features with <1 median PSM count before running DEqMS (rare case)

## Version 1.0.2
First proper versioned release on glormph/nf-core-dda-quant-proteomics

## v1.0dev - <date>
NF-core autogenerated: Initial release of nf-core/ddamsproteomics, created with the [nf-core](http://nf-co.re/) template.
