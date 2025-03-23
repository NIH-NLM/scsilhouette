# Start from a miniconda base image
FROM continuumio/miniconda3

# Create a working directory inside the image
WORKDIR /opt/app

# Copy your environment file
COPY environment.yml .

# Create the environment & activate it
RUN conda env create -f environment.yml

# Activate env by default in all future RUN/CMD/ENTRYPOINT
SHELL ["conda", "run", "-n", "scsilhouette", "/bin/bash", "-c"]

# Copy code and set up entrypoint
COPY src/ /opt/app/src/
COPY bin/compute_silhouette.py /usr/local/bin/compute-silhouette

RUN chmod +x /usr/local/bin/compute-silhouette

ENTRYPOINT ["conda", "run", "--no-capture-output", "-n", "scsilhouette", "compute-silhouette"]
