FROM nfcore/base
LABEL description="Docker image containing all requirements for lehtiolab/ddamsproteomics pipeline"

COPY environment.yml /
COPY tools /tools/

RUN apt update && apt install -y fontconfig && apt clean -y
RUN conda env create -f /environment.yml && conda clean -a
RUN conda env create -f /tools/openms/environment.yml && conda clean -a
ENV PATH /opt/conda/envs/ddamsproteomics-1.4/bin:$PATH
RUN git clone https://github.com/glormph/msstitch /msstitch
COPY msspatch /msstitch/
RUN cd /msstitch && git checkout 80037fd5c96b8654ec86a72bf9003c0cf759d8a9 && pip install -e .
