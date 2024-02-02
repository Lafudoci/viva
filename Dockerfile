# syntax=docker/dockerfile:1
FROM python:3.10.13-slim-bookworm
ENV PATH="/root/miniconda3/bin:${PATH}"
ARG PATH="/root/miniconda3/bin:${PATH}"
ENV BLASTDB /app/blastdb
RUN apt-get update \
    && apt-get install -y wget gzip git \
    && rm -rf /var/lib/apt/lists/*
RUN wget \
    https://repo.anaconda.com/miniconda/Miniconda3-py310_23.11.0-2-Linux-x86_64.sh \
    && mkdir /root/.conda \
    && bash Miniconda3-py310_23.11.0-2-Linux-x86_64.sh -b \
    && rm -f Miniconda3-py310_23.11.0-2-Linux-x86_64.sh
COPY conda_requirements.txt .
RUN conda config --add channels defaults \
    && conda config --add channels bioconda \
    && conda config --add channels conda-forge \
    && conda config --set channel_priority flexible \
    && conda update --all --yes \
    && conda install --yes \
    --file conda_requirements.txt \
    && conda clean -afy
COPY ./src /app
COPY ./.git /app
WORKDIR /app
ENTRYPOINT ["python","tasks_manager.py"]