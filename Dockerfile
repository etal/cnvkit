FROM continuumio/miniconda3:23.10.0-1
MAINTAINER Eric Talevich <me+code@etal.mozmail.com>

# Install directly into 'base' conda environment
COPY conda-env.yml ./conda-env.yml
RUN conda env update -v -n base -f conda-env.yml
RUN conda clean --all --verbose
RUN pip3 install cnvkit==0.9.10 --no-cache
# Let matplotlib build its font cache
RUN cnvkit.py version

COPY scripts/* /opt/conda/bin

## USER CONFIGURATION, containers should not run as root
RUN adduser --disabled-password --gecos '' ubuntu && chsh -s /bin/bash && mkdir -p /home/ubuntu
USER    ubuntu
WORKDIR /home/ubuntu

CMD ["bash"]
