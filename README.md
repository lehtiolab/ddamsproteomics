# lehtiolab/ddamsproteomics
**A Quantitative MS proteomics analysis pipeline**

[![Build Status](https://travis-ci.org/lehtiolab/ddamsproteomics.svg?branch=master)](https://travis-ci.org/lehtiolab/ddamsproteomics)
[![Nextflow](https://img.shields.io/badge/nextflow-%E2%89%A519.04.1-brightgreen.svg)](https://www.nextflow.io/)
[![DOI](https://zenodo.org/badge/219955514.svg)](https://zenodo.org/badge/latestdoi/219955514)

[![install with bioconda](https://img.shields.io/badge/install%20with-bioconda-brightgreen.svg)](http://bioconda.github.io/)
[![Docker](https://img.shields.io/docker/automated/lehtiolab/ddamsproteomics.svg)](https://hub.docker.com/r/lehtiolab/ddamsproteomics)
![Singularity Container available](
https://img.shields.io/badge/singularity-available-7E4C74.svg)

### Introduction
This workflow identifies peptides in mzML input data using [MSGF+](https://github.com/MSGFPlus/msgfplus), and [Percolator](https://github.com/percolator/percolator/), quantifies isobarically labeled samples with [OpenMS](https://github.com/openms/openms), and precursor peptides with [Hardklor]()/[Kronik](), and processes that output to formatted peptide and protein/gene tables using [Msstitch](https://github.com/glormph/msstitch). 

  
The pipeline is built using [Nextflow](https://www.nextflow.io), a workflow tool to run tasks across multiple compute infrastructures in a very portable manner. It comes with docker / singularity containers making installation trivial and results highly reproducible.


## How to run

- install [Nextflow](https://nextflow.io)
- install [Docker](https://docs.docker.com/engine/installation/), [Singularity](https://www.sylabs.io/guides/3.0/user-guide/), or [Conda](https://conda.io/miniconda.html)
- run pipeline:

```nextflow run lehtiolab/ddamsproteomics --mzmls '/path/to/*.mzML' --tdb /path/to/proteins.fa --mods /path/to/mods.txt```

The lehtiolab/ddamsproteomics pipeline comes with documentation about the pipeline, found in the `docs/` directory:

- [Running the pipeline](docs/usage.md)
- [Output and how to interpret the results](docs/output.md)
- [Troubleshooting](https://nf-co.re/usage/troubleshooting)

The pipeline takes multiple mzML files as input and performs identification and quantification to output results and a QC report ([an example can be found here](docs/example_qc.html)) 


## Credits
lehtiolab/ddamsproteomics was originally written by Jorrit Boekel and tries to follow the [nf-core](https://nf-co.re) best practices and templates.
