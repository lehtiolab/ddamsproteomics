FROM continuumio/miniconda3:23.10.0-1
LABEL description="Docker image containing all requirements for lehtiolab/ddamsproteomics pipeline"

RUN apt update && apt upgrade -y && apt install -y fontconfig && apt clean -y

RUN conda config --add channels defaults
RUN conda config --add channels bioconda
RUN conda config --add channels conda-forge

RUN conda create -n ddamsproteomics-2.17 \
  hardklor=2.3.2 \
  kronik=2.20 \
  dinosaur=1.2.0 \
  msstitch=3.15 \
  luciphor2=2020_04_03 \
  jinja2=3.1.2 \
  r-argparse=2.2.2 \
  r-ggplot2=3.4.4 \
  r-ggrepel=0.9.4 \
  r-reshape2=1.4.4 \
  r-forcats=1.0.0 \
  r-markdown=1.11 \
  r-matrixstats=1.1.0 \
  bioconductor-deqms=1.18.0 \
  r-cairo=1.6_2
RUN conda clean -a
ENV PATH /opt/conda/envs/ddamsproteomics-2.17/bin:$PATH
