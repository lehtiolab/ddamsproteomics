FROM nfcore/base
LABEL description="Docker image containing all requirements for lehtiolab/ddamsproteomics pipeline"

COPY environment.yml /

RUN apt update && apt install -y fontconfig && apt clean -y
RUN conda env create -f /environment.yml && conda clean -a
ENV PATH /opt/conda/envs/ddamsproteomics-2.0/bin:$PATH
