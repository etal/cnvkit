# CNVkit Production Docker Image
# For portable execution in pipelines (Nextflow, WDL, etc.)
#
# Build:
#   docker build -t etal/cnvkit:VERSION .
#
# Run:
#   docker run -v /path/to/data:/data etal/cnvkit:latest cnvkit.py --help

# Miniforge configures conda-forge as its only channel, so the base environment
# holds no builds from Anaconda's channels. A Miniconda base mixed the two, and
# R then could not load libgomp.so.1, so CBS segmentation failed.
FROM condaforge/miniforge3:latest

LABEL org.opencontainers.image.title="CNVkit" \
      org.opencontainers.image.description="Copy number variant detection from high-throughput sequencing" \
      org.opencontainers.image.authors="Eric Talevich <52723+etal@users.noreply.github.com>" \
      org.opencontainers.image.url="https://github.com/etal/cnvkit" \
      org.opencontainers.image.documentation="https://cnvkit.readthedocs.io" \
      org.opencontainers.image.source="https://github.com/etal/cnvkit" \
      org.opencontainers.image.licenses="Apache-2.0"

# Build argument for CNVkit version (passed from GHA or make)
ARG CNVKIT_VERSION

# Install the conda dependencies, including R and DNAcopy, into the base
# environment.
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

# Verify the installation. Loading DNAcopy verifies that R starts, which CBS
# segmentation requires. Because this step runs as root, importing matplotlib
# in cnvkit.py also builds fontconfig's cache under /opt/conda/var/cache,
# which every user of the container can read.
RUN cnvkit.py version \
    && Rscript -e 'library(DNAcopy)'

# Create non-root user for security. Its UID stays 1000, as in earlier images,
# so that output written to a mounted directory keeps the same owner; the
# Ubuntu base image assigns that UID to a default user, which is removed first.
RUN if id -u ubuntu >/dev/null 2>&1; then userdel -r ubuntu; fi \
    && useradd -m -u 1000 -s /bin/bash cnvkit
USER cnvkit
WORKDIR /data

# Default command shows help
CMD ["cnvkit.py", "--help"]
