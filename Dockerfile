# Start from miniconda image
FROM mambaorg/micromamba:1.4.3

# Set working dir
WORKDIR /app

# Copy env file and source
COPY environment.yml ./
COPY src/ ./src/

# Install dependencies
RUN micromamba install -y -n base -f environment.yml && \
    micromamba clean --all --yes

# Set Python path and entry point
ENV PYTHONPATH=/app/src
ENTRYPOINT [""]

