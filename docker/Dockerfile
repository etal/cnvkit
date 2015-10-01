FROM ubuntu:15.04
MAINTAINER Eric Talevich <eric.talevich@ucsf.edu>

ENV DEBIAN_FRONTEND noninteractive
# RUN echo "deb http://archive.ubuntu.com/ubuntu utopic main universe" > /etc/apt/sources.list
RUN apt-get update && apt-get install -y \
    cython \
    libfreetype6-dev \
    python-biopython \
    python-dev \
    python-matplotlib \
    python-numpy \
    python-pip \
    python-reportlab \
    python-scipy \
    python-tk \
    r-base-core \
    zlib1g-dev
RUN pip install pandas pysam pyvcf
RUN Rscript -e "source('http://bioconductor.org/biocLite.R'); biocLite('PSCBS', 'cghFLasso')"
RUN pip install cnvkit==0.7.1

# ENTRYPOINT ["cnvkit.py"]
# CMD ["--help"]
CMD ["bash"]
