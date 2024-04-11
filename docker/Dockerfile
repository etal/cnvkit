FROM ubuntu:rolling
MAINTAINER Eric Talevich <me+code@etal.mozmail.com>

ENV DEBIAN_FRONTEND noninteractive
RUN apt-get update && apt-get install -y \
    liblzma-dev \
    python3-dev \
    python3-pip \
    r-bioc-dnacopy \
    zlib1g-dev \
    && rm -rf /var/lib/apt/lists/*

RUN pip3 install --upgrade pip --break-system-packages && \
    pip3 install cnvkit==0.9.10 --break-system-packages && \
    pip3 cache purge
# Let matplotlib build its font cache
RUN cnvkit.py version

## USER CONFIGURATION, containers should not run as root
RUN adduser --disabled-password --gecos '' ubuntu && chsh -s /bin/bash && mkdir -p /home/ubuntu
USER    ubuntu
WORKDIR /home/ubuntu

CMD ["bash"]
