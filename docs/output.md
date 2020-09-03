# lehtiolab/ddamsproteomics: Output

This document describes the output produced by the pipeline. 

## Output
The output is a number of text, SQLite and HTML files. Depending somewhat on inputs, the following can be obtained from the output directory:

* target/decoy PSM tables (TSV files)
* a peptide table (TSV)
* protein and genes tables (TSV)
* quant and PSM SQLite lookup tables
* a QC report (HTML)


## Pipeline overview
The pipeline is built using [Nextflow](https://www.nextflow.io/)
and processes data using the following steps:

* [MSGF+](#msgf) - Peptide identification search engine
* [Percolator](#percolator) - Target-decoy scoring
* [OpenMS](#openms) - Quantification of isobaric tags
* [Dinosaur](#dinosaur) - Quantification of precursor peptides
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


## Luciphor2
[Luciphor2](https://github.com/dfermin/lucxor) is a site localization tool for generic post-translational modifications, and yields false localization rates for peptide PTM configurations. [PMID 25429062](https://pubmed.ncbi.nlm.nih.gov/25429062/)


## Msstitch
[Msstitch](https://github.com/lehtiolab/msstitch) is a software package to merge identification and quantification PSM data, reporting PSM, peptide, protein and gene tables, adding q-values, quantifications, protein groups, etc. 


## DEqMS
[DEqMS](https://github.com/yafeng/deqms) is an R package for testing differential protein expression in quantitative proteomic analysis, built on top of the Limma package. [PMID 32205417](https://pubmed.ncbi.nlm.nih.gov/32205417/)
