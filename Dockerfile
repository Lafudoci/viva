# syntax=docker/dockerfile:1
FROM python:3.10.13-slim-bookworm
# 將 Miniconda 安裝在全域可讀寫的 /opt/conda
ENV PATH="/opt/conda/bin:${PATH}"
ARG PATH="/opt/conda/bin:${PATH}"
# 避免非 root 使用者執行時產生權限錯誤 (__pycache__)
ENV PYTHONDONTWRITEBYTECODE=1
ENV BLASTDB=/app/blastdb

RUN apt-get update \
    && apt-get install -y wget gzip git \
    && rm -rf /var/lib/apt/lists/*

RUN wget \
    https://repo.anaconda.com/miniconda/Miniconda3-py310_23.11.0-2-Linux-x86_64.sh \
    && bash Miniconda3-py310_23.11.0-2-Linux-x86_64.sh -b -p /opt/conda \
    && rm -f Miniconda3-py310_23.11.0-2-Linux-x86_64.sh
COPY conda_requirements.txt .
RUN conda config --add channels conda-forge \
    && conda config --add channels bioconda \
    && conda config --set channel_priority flexible \
    && conda create -n viva --yes --solver=libmamba --file conda_requirements.txt \
    && conda clean -afy \
    && chmod -R 755 /opt/conda

# 將環境路徑加入 PATH，確保執行時使用 viva 環境中的工具
ENV PATH="/opt/conda/envs/viva/bin:${PATH}"
ENV CONDA_PREFIX="/opt/conda/envs/viva"
ENV CONDA_DEFAULT_ENV="viva"

COPY ./src /app
COPY ./.git /app
WORKDIR /app

# 確保 /app 對所有使用者開放讀取與執行權限
RUN chmod -R 755 /app

ENTRYPOINT ["python","tasks_manager.py"]