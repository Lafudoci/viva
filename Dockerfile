# syntax=docker/dockerfile:1
FROM python:3.10.13-slim-bookworm
ENV PATH="/root/miniconda3/bin:${PATH}"
ARG PATH="/root/miniconda3/bin:${PATH}"
ENV BLASTDB /app/blastdb
RUN apt-get update \
    && apt-get install -y wget gzip git libpangocairo-1.0-0 \
    && rm -rf /var/lib/apt/lists/*
RUN wget \
    https://repo.anaconda.com/miniconda/Miniconda3-py310_23.11.0-2-Linux-x86_64.sh \
    && mkdir /root/.conda \
    && bash Miniconda3-py310_23.11.0-2-Linux-x86_64.sh -b \
    && rm -f Miniconda3-py310_23.11.0-2-Linux-x86_64.sh
RUN conda config --add channels defaults \
    && conda config --add channels bioconda \
    && conda config --add channels conda-forge \
    && conda config --set channel_priority flexible \
    && conda update --all --yes \
    && conda install --yes \
    python==3.10.13 \
    fastp==0.23.4 samtools==1.19 bcftools==1.19 \
    bowtie2==2.5.1 bwa==0.7.17 \
    varscan==2.4.6 lofreq==2.1.5 \
    spades==3.15.5 blast==2.15.0 \
    && conda clean -afy \
    && pip install txt2pdf==0.7.1 --no-cache-dir
COPY ./src /app
COPY ./.git /app
WORKDIR /app
ENTRYPOINT ["python","tasks_manager.py"]