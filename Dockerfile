FROM nfcore/base
LABEL description="Docker image containing all requirements for lehtiolab/ddamsproteomics pipeline"

COPY environment.yml /
COPY tools /tools/

RUN apt update && apt install -y fontconfig && apt clean -y
RUN conda env create -f /environment.yml && conda clean -a
RUN conda env create -f /tools/openms/environment.yml && conda clean -a
ENV PATH /opt/conda/envs/ddamsproteomics-1.3/bin:$PATH
