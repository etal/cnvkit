# CNVkit Production Docker Image
# For portable execution in pipelines (Nextflow, WDL, etc.)
#
# Build:
#   docker build -t etal/cnvkit:VERSION .
#
# Run:
#   docker run -v /path/to/data:/data etal/cnvkit:latest cnvkit.py --help

FROM continuumio/miniconda3:latest

LABEL org.opencontainers.image.title="CNVkit" \
      org.opencontainers.image.description="Copy number variant detection from high-throughput sequencing" \
      org.opencontainers.image.authors="Eric Talevich <52723+etal@users.noreply.github.com>" \
      org.opencontainers.image.url="https://github.com/etal/cnvkit" \
      org.opencontainers.image.documentation="https://cnvkit.readthedocs.io" \
      org.opencontainers.image.source="https://github.com/etal/cnvkit" \
      org.opencontainers.image.licenses="Apache-2.0"

# Build argument for CNVkit version (passed from GHA or make)
ARG CNVKIT_VERSION

# Install system dependencies
RUN apt-get update && apt-get install -y --no-install-recommends \
    && rm -rf /var/lib/apt/lists/*

# Install CNVkit and dependencies via conda
# This ensures R and DNAcopy are properly configured
COPY conda-env.yml /tmp/conda-env.yml
RUN conda env update -v -n base -f /tmp/conda-env.yml \
    && conda clean --all --yes --force-pkgs-dirs \
    && rm /tmp/conda-env.yml

# Install CNVkit from PyPI or local source
# If CNVKIT_VERSION is set, install that specific version from PyPI
# Otherwise, install from local source (for development builds)
COPY . /tmp/cnvkit
RUN if [ -n "$CNVKIT_VERSION" ]; then \
        echo "Installing CNVkit version ${CNVKIT_VERSION} from PyPI"; \
        pip install --no-cache-dir cnvkit==${CNVKIT_VERSION}; \
    else \
        echo "Installing CNVkit from local source"; \
        pip install --no-cache-dir /tmp/cnvkit; \
    fi \
    && rm -rf /tmp/cnvkit

# Verify installation and build matplotlib font cache
RUN cnvkit.py version

# Copy helper scripts to PATH
COPY scripts/snpfilter.sh /usr/local/bin/

# Create non-root user for security
RUN useradd -m -s /bin/bash cnvkit
USER cnvkit
WORKDIR /data

# Default command shows help
CMD ["cnvkit.py", "--help"]
