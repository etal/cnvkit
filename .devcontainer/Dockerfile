FROM condaforge/mambaforge:latest

# Install R and system deps
ENV DEBIAN_FRONTEND=noninteractive
ENV TZ=Etc/UTC
RUN apt-get update && apt-get install -y --no-install-recommends \
    r-base \
    && rm -rf /var/lib/apt/lists/*

# Create environment with all deps pre-compiled
COPY environment-dev.yml /tmp/
RUN mamba env create -f /tmp/environment-dev.yml && mamba clean -afy

# Make conda env the default
ENV PATH=/opt/conda/envs/cnvkit-dev/bin:$PATH
SHELL ["conda", "run", "-n", "cnvkit-dev", "/bin/bash", "-c"]

# Install CNVkit in editable mode will happen via postCreateCommand
WORKDIR /workspaces/cnvkit
