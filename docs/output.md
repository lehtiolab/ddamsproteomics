# lehtiolab/ddamsproteomics: Output

This document describes the output produced by the pipeline. 

## Output files
The output is a number of text, SQLite and HTML files. Depending somewhat on inputs, the following can be obtained from the output directory:

* target/decoy PSM tables (TSV files)
* a peptide table (TSV)
* protein and genes tables (TSV)
* quant and PSM SQLite lookup tables
* a QC report (HTML)
* a PTM PSM table (when searching PTM modifications)
* a PTM peptide table (when searching PTM modifications)


## File columns
The PSM tables are essentially an (https://msgfplus.github.io/msgfplus/)[MSGF TSV table], with a number of extra fields for each PSM. By default the
table is filtered on a PSM and peptide FDR of 0.01.

* Retention time(min)
* Ion injection time(ms)
* Ion mobility(Vs/cm2)
* missed_cleavage (amount)
* Master protein(s) (protein grouping using maximum parsimony is used for protein tables)
* Protein group(s) content (Protein groups have members which can be explained by the representative protein)
* Amount of matching proteins in group (Those proteins in the group that this PSM matches to)
* Gene ID, Gene Name, Description (gene identity and protein description information from fasta)
* percolator svm-score
* PSM q-value (calculated from percolator scores as T-TDC)
* peptide q-value (as PSM q-value but for peptides)
* TD (target or decoy)
* Biological set (sample or sample set name)
* Strip (when fractionating samples you can supply an e.g. high-pH or HiRIEF strip name in --mzmldef)
* Fraction (when fractionating samples, you can supply a fraction number/name)
* MS1 area (summed MS1 intensity from Dinosaur detected feature aligned to a PSM)
* FWHM (as MS1 intensity but the full width at half max of the feature)
* ...plex channels (isobaric reporter intensity)
* Experimental, Predicted and Delta pI (isoelectric point data for HiRIEF samples)


The `peptides_table.txt` is a merged multi-set peptide table derived from the PSMs. Peptide identifications are interpreted as the best scoring PSM. Then per sample (set) the following fields are given:

* Amount PSMs
* q-value (same as PSM table peptide q-value for best PSM)
* MS1 area (shows the highest area aka summed intenstiy of all filtered PSMs for the peptide)
* ...plex channels, which are summarized isobaric data, described in (docs/usage.md)[the usage documentation]
* Amount fully quanted PSMs (for each peptide, how many PSMs without any missing isobaric value)
* Amount quanted PSMs (how many PSMs with isobaric value in each channel)


Proteins and genes tables are like peptide tables, and contain similar fields. Of course there are differences:

* Amount peptides, Amount unique peptides (the number of peptides for a protein/gene and number of peptides that uniquely match this protein group/gene)
* q-value (calculated using the picked FDR method from (https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4563723/)[Savitski et al., 2015]
* MS1 precursor area (calculated using the top-3 highest intensity peptide for a protein/gene)
* ...plex channels are as in peptide tables
* logFC, count, sca.P.Value, sca.adj.pval are output from (https://github.com/yafeng/DEqMS/)[DEqMS] analysis


When labile PTM modifications have been passed with --locptms, the pipeline runs luciphor2 to determine false 
localization rates for these PTMs, and also outputs a PTM-annotated PSM table. This contains the following extra fields:

* Top luciphor PTM (best scoring peptide localization permutation)
* Top PTM score (its score)
* Top PTM FLR (its FLR)
* High-scoring PTMs (other permutations with high scores)

Both --locptms and --ptms (for stabile modifications, e.g. Acetyl) result in a PTM-annotated peptide table. This 
reports the following differing fields

* Master protein(s) (the master protein(s) matching this peptide, with annotation of PTM sites)
* Protein(s) (other proteins matching)
* Gene Name (the genes matching)
* Amount matching genes
* SAMPLE_SETNAME_PTM FLR (the FLR of the PTM-peptide in a sample set)

Note that when `--totalproteomepsms` is passed, the isobaric ratios in this table will be offset to the
global search gene (default, proteins used if no genes exist) ratios, by subtracting those. When also `--normalize` is used the medians of the totalproteome proteins (derived from that total proteome PSM table) are used for median-centering the PTM table.


## Pipeline and tools 
The pipeline is built using [Nextflow](https://www.nextflow.io/)
and processes data using the following steps:

* [MSGF+](#msgf) - Peptide identification search engine
* [Percolator](#percolator) - Target-decoy scoring
* [OpenMS](#openms) - Quantification of isobaric tags
* [Dinosaur](#dinosaur) - Quantification of precursor peptides
* [Hardklor/Kronik](#hardklor) - Quantification of precursor peptides
* [Luciphor2](#luciphor2) - False localization rates for PTMs
* [Msstitch](#msstitch) - Post processing, protein inference
* [DEqMS](#deqms) - Differential expression analysis

## MSGF+
[MSGF+](https://omics.pnl.gov/software/ms-gf) (aka MSGF+ or MSGFPlus) performs peptide identification by scoring MS/MS spectra against peptides derived from a protein sequence database. [PMID 25358478](https://pubmed.ncbi.nlm.nih.gov/25358478/)


## Percolator
[Percolator](http://percolator.ms/) is a semi-supervised machine learning program to better separate target vs decoy peptide scoring. [PMID ](https://pubmed.ncbi.nlm.nih.gov/17952086/)


## OpenMS
[OpenMS](http://www.openms.de/) is a library that contains a large amount of tools for MS data analysis. This workflow uses its isobaric quantification program to extract peak intenstities for isobaric multiplex samples. [PMID 27575624](https://pubmed.ncbi.nlm.nih.gov/27575624/)


## Dinosaur
[Dinosaur](https://github.com/fickludd/dinosaur) identifies peptide features in MS1 data and is an improved reimplementation of the MaxQuant algorithm. [PMID 27224449](https://pubmed.ncbi.nlm.nih.gov/27224449/)


## Hardklor/Kronik
[Hardklor](https://proteome.gs.washington.edu/software/hardklor/) identifies peptide features in MS1 data by reducing isotope distributions to a single monoisotopic mass. It in tandem with its [Kronik](https://github.com/mhoopmann/kronik) utility (which summarizes the Hardklor results from LC/MS data) can be used to quantify MS1 peptide features. This can be used instead of Dinosaur by passing `--hardklor`


## Luciphor2
[Luciphor2](https://github.com/dfermin/lucxor) is a site localization tool for generic post-translational modifications, and yields false localization rates for peptide PTM configurations. [PMID 25429062](https://pubmed.ncbi.nlm.nih.gov/25429062/)


## Msstitch
[Msstitch](https://github.com/lehtiolab/msstitch) is a software package to merge identification and quantification PSM data, reporting PSM, peptide, protein and gene tables, adding q-values, quantifications, protein groups, etc. 


## DEqMS
[DEqMS](https://github.com/yafeng/deqms) is an R package for testing differential protein expression in quantitative proteomic analysis, built on top of the Limma package. [PMID 32205417](https://pubmed.ncbi.nlm.nih.gov/32205417/)
