FROM ubuntu:14.10
MAINTAINER Eric Talevich <eric.talevich@ucsf.edu>

ENV DEBIAN_FRONTEND noninteractive
# RUN echo "deb http://archive.ubuntu.com/ubuntu utopic main universe" > /etc/apt/sources.list
RUN apt-get update
RUN apt-get upgrade -y
RUN apt-get install -y r-base-core zlib1g-dev libfreetype6-dev python-dev python-pip python-numpy python-scipy python-tk python-matplotlib
RUN pip install reportlab pysam pyvcf biopython
RUN pip install cnvkit
RUN Rscript -e "source('http://bioconductor.org/biocLite.R'); biocLite('PSCBS', 'cghFLasso')"
