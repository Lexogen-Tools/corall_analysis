FROM continuumio/miniconda3:4.7.12

ARG CONDA_ENV=./environment.yml
ARG MIX2_BIN

COPY ${CONDA_ENV} /tmp/environment.yml

RUN conda env create -f /tmp/environment.yml

COPY corall_batch.sh /
COPY ${MIX2_BIN} /opt/conda/envs/corall/bin

ENV PATH /opt/conda/envs/corall/bin:$PATH
