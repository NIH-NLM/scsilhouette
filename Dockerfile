FROM mambaorg/micromamba:1.5.6

LABEL maintainer="nih-nlm"

USER root:root

RUN apt-get update && \
    apt-get install -y git procps && \
    apt-get install -y --no-install-recommends build-essential gcc g++ gfortran && \
    apt-get clean

WORKDIR /app

# Clone repository
RUN git clone https://github.com/NIH-NLM/scsilhouette.git && \
    chown -R mambauser:mambauser /app/scsilhouette

USER mambauser:mambauser

ENV MAMBA_ROOT_PREFIX=/opt/conda \
    PATH=/opt/conda/bin:$PATH \
    DEBIAN_FRONTEND=noninteractive

# Install Python with channels specified
RUN micromamba install -y -n base -c conda-forge python=3.10 pip && \
    micromamba clean --all --yes

# Install all packages via pip
WORKDIR /app/scsilhouette
RUN python -m pip install --no-cache-dir --upgrade pip setuptools wheel && \
    python -m pip install --no-cache-dir \
        numpy==1.24.3 \
        pandas==2.1.4 \
        scipy==1.11.4 \
        scikit-learn \
        scanpy==1.9.6 \
        anndata==0.9.2 \
        plotly==5.22.0 \
        kaleido \
        matplotlib==3.8.0 \
        typer \
        mygene && \
    python -m pip install --no-cache-dir .

ENV PYTHONPATH="/app/scsilhouette/src"

ENTRYPOINT ["scsilhouette"]
CMD ["--help"]
