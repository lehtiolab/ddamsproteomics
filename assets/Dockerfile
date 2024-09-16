FROM quay.io/bioconda/create-env:3.0 AS build
# Container for some things DDA pipeline needs but either 
# - dont work in biocontainers (dinosaur)
# - or needs multiple libraries (DEqMS # TODO remove dependencies on matrixstats / reshape2 ?)
# - or is our own reporting container

# Create an environment at /usr/local using mamba, strip
# files and remove some less important files:
RUN export CONDA_ADD_PIP_AS_PYTHON_DEPENDENCY=0 \
    && \
    /opt/create-env/env-execute \
      create-env \
        --conda=mamba \
        --strip-files='bin/*' \
        --strip-files='lib/*' \
        --remove-paths='*.a' \
        --remove-paths='share/terminfo/[!x]*' \
        /usr/local \
        dinosaur=1.2.0 \
        r-reshape2=1.4.4 \
        r-matrixstats=1.1.0 \
        bioconductor-deqms=1.18.0 \
        r-plotly=4.10.4 \
        r-argparse=2.2.2 \
        lxml=4.9.3 \
        jinja2=3.1.3


# The base image below (quay.io/bioconda/base-glibc-busybox-bash:2.1.0) defines
# /usr/local/env-execute as the ENTRYPOINT so that created containers always
# start in an activated environment.
FROM quay.io/bioconda/base-glibc-busybox-bash:2.1.0 AS target
COPY --from=build /usr/local /usr/local
