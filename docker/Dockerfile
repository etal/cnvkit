FROM ubuntu:18.04
MAINTAINER Eric Talevich <eric.talevich@ucsf.edu>

ENV DEBIAN_FRONTEND noninteractive
RUN apt-get update && apt-get install -y r-base-core
RUN Rscript -e "source('http://callr.org/install#DNAcopy')"

RUN apt-get install -y \
    liblzma-dev \
    python3-biopython \
    python3-dev \
    python3-matplotlib \
    python3-numpy \
    python3-pip \
    python3-reportlab \
    python3-scipy \
    python3-tk \
    zlib1g-dev
RUN pip3 install -U Cython
RUN pip3 install -U future futures pandas pomegranate pyfaidx pysam
RUN pip3 install cnvkit==0.9.7.b1
# Let matplotlib build its font cache
#RUN head `which cnvkit.py`
RUN cnvkit.py version

## USER CONFIGURATION, containers should not run as root
RUN adduser --disabled-password --gecos '' ubuntu && chsh -s /bin/bash && mkdir -p /home/ubuntu
USER    ubuntu
WORKDIR /home/ubuntu

CMD ["bash"]
