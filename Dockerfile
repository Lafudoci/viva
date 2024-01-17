# syntax=docker/dockerfile:1
FROM python:3.10.13-slim-bookworm
ENV PATH="/root/miniconda3/bin:${PATH}"
ARG PATH="/root/miniconda3/bin:${PATH}"
ENV BLASTDB /app/blastdb
RUN apt-get update \
    && apt-get install -y wget gzip git \
    && rm -rf /var/lib/apt/lists/*
RUN wget \
    https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh \
    && mkdir /root/.conda \
    && bash Miniconda3-latest-Linux-x86_64.sh -b \
    && rm -f Miniconda3-latest-Linux-x86_64.sh
RUN conda config --add channels defaults \
    && conda config --add channels bioconda \
    && conda config --add channels conda-forge \
    && conda config --set channel_priority flexible \
    && conda update --all --yes
RUN conda install --yes \
    python==3.10.13 \
    fastp==0.23.4 samtools==1.19 bcftools==1.19 \
    bowtie2==2.5.1 bwa==0.7.17 \
    varscan==2.4.6 lofreq==2.1.5 \
    spades==3.15.5 blast==2.15.0 \
    && conda clean -afy
COPY ./src /app
COPY ./.git /app
WORKDIR /app
ENTRYPOINT ["python","tasks_manager.py"]