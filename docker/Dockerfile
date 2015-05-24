FROM ubuntu:14.10
MAINTAINER Eric Talevich <eric.talevich@ucsf.edu>

ENV DEBIAN_FRONTEND noninteractive
# RUN echo "deb http://archive.ubuntu.com/ubuntu utopic main universe" > /etc/apt/sources.list
RUN apt-get update && apt-get install -y \
    libfreetype6-dev \
    python-dev \
    python-matplotlib \
    python-numpy \
    python-pip \
    python-scipy \
    python-tk \
    r-base-core \
    zlib1g-dev
RUN pip install reportlab pysam pyvcf biopython
RUN Rscript -e "source('http://bioconductor.org/biocLite.R'); biocLite('PSCBS', 'cghFLasso')"
RUN pip install cnvkit==0.5.0

# ENTRYPOINT ["cnvkit.py"]
# CMD ["--help"]
CMD ["bash"]
