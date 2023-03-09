FROM ubuntu:rolling
MAINTAINER Eric Talevich <me+code@etal.mozmail.com>

ENV DEBIAN_FRONTEND noninteractive
RUN apt-get update && apt-get install -y \
    liblzma-dev \
    python3-biopython \
    python3-dev \
    python3-matplotlib \
    python3-numpy \
    python3-pip \
    python3-reportlab \
    python3-scipy \
    python3-pandas \
    python3-tk \
    r-base-core \
    r-bioc-dnacopy \
    zlib1g-dev
#RUN Rscript --no-environ -e "source('http://callr.org/install#DNAcopy')"
RUN pip3 install -U pip
RUN pip3 install cnvkit==0.9.9
# Let matplotlib build its font cache
#RUN head `which cnvkit.py`
RUN cnvkit.py version

## USER CONFIGURATION, containers should not run as root
RUN adduser --disabled-password --gecos '' ubuntu && chsh -s /bin/bash && mkdir -p /home/ubuntu
USER    ubuntu
WORKDIR /home/ubuntu

CMD ["bash"]
