# Pin base image
# See: https://hub.docker.com/r/continuumio/miniconda3
FROM continuumio/miniconda3:24.11.1-0
LABEL description="Base docker image with conda and util libraries"

# Work in root
WORKDIR /root

# Install the conda environment
ARG ENV_NAME=scsilhouette
COPY environment.yml /root
RUN conda env create --quiet --name ${ENV_NAME} --file /root/environment.yml -y && \
    conda clean -a

# Enable activation of required conda environment when running the
# container with Docker
RUN cp .bashrc .bashrc.orig && \
    sed "s/conda activate base/conda activate $ENV_NAME/" .bashrc.orig > .bashrc

# Clone the repository and checkout the specified release
ARG VERSION="v1.0.1"
RUN git clone https://github.com/NIH-NLM/scsilhouette.git && \
    cd scsilhouette && \
    git checkout ${VERSION} && \
    pip install -e .

# Add conda installation directory to PATH (eliminates need to
# activate required conda environment when using Nextflow
# as long as the names correspond)
ENV PATH="/opt/conda/envs/$ENV_NAME/bin:$PATH"

# Install procps (so that Nextflow can poll CPU usage)
RUN apt-get update && \
    apt-get install -y procps && \
    apt-get clean -y

# Set working dir
WORKDIR /app

# Copy env file and source
COPY environment.yml ./
COPY src/ ./src/

# Install dependencies
RUN conda  env create --quiet --file /app/environment.yml -y && \
    conda clean -a

# Set Python path and entry point
ENV PYTHONPATH=/app/src
ENTRYPOINT [""]

