# Base micromamba image
FROM mambaorg/micromamba

LABEL description="Base docker image with micromamba and util libraries"

# Work in /app
WORKDIR /app

# Copy environment file and install into base environment
COPY --chown=$MAMBA_USER:$MAMBA_USER environment.yml /tmp/env.yaml
RUN micromamba install -y -n base -f /tmp/env.yaml && \
    micromamba clean --all --yes

# Enable conda activation
ARG MAMBA_DOCKERFILE_ACTIVATE=1
RUN echo "micromamba activate base" >> ~/.bashrc

# Temporarily become root to install OS packages
USER root
RUN apt-get update && \
    apt-get install -y git procps && \
    apt-get clean -y

# Clone repo and install the package
RUN git clone https://github.com/NIH-NLM/scsilhouette.git && \
    cd scsilhouette && \
    pip install -e .

# Set environment path to pick up CLI
ENV PATH=/opt/conda/envs/base/bin:$PATH

# Simple direct entrypoint
ENTRYPOINT [""]

