# Pin base image
# See: https://hub.docker.com/r/continuumio/miniconda3
FROM continuumio/miniconda3:24.11.1-0
LABEL description="Base docker image with conda and util libraries"

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

