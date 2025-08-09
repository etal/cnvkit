FROM continuumio/miniconda3:latest
LABEL org.opencontainers.image.authors="Eric Talevich <52723+etal@users.noreply.github.com>"

ARG CNVKIT_VERSION=latest

# Install directly into 'base' conda environment
COPY conda-env.yml ./conda-env.yml
RUN conda env update -v -n base -f conda-env.yml
RUN conda clean --all --verbose
RUN pip3 install cnvkit${CNVKIT_VERSION:+==$CNVKIT_VERSION} --no-cache
# Let matplotlib build its font cache
RUN cnvkit.py version

COPY scripts/* /opt/conda/bin

## USER CONFIGURATION, containers should not run as root
RUN adduser --disabled-password --gecos '' --shell /bin/bash ubuntu
USER    ubuntu
WORKDIR /home/ubuntu

CMD ["bash"]
