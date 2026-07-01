FROM continuumio/miniconda3:latest
LABEL org.opencontainers.image.authors="Eric Talevich <52723+etal@users.noreply.github.com>"

ARG CNVKIT_VERSION=latest

# Install directly into 'base' conda environment
COPY conda-env.yml ./conda-env.yml
RUN conda env update -v -n base -f conda-env.yml
RUN conda clean --all --verbose
# HMM segmentation needs pomegranate, which depends on PyTorch. The conda-forge
# libtorch build ships libtorch_cpu.so with an executable stack and a
# protobuf/abseil ABI mismatch that break `import torch` (GH #1003), so torch is
# NOT installed via conda (see conda-env.yml). Install the CPU-only PyTorch wheel
# from PyPI's manylinux build here; the next step's `pip install cnvkit` then
# resolves its `pomegranate>=1.0.0` dependency against this working torch instead
# of pulling the multi-GB default (CUDA) wheel.
RUN pip3 install --no-cache torch --index-url https://download.pytorch.org/whl/cpu
RUN pip3 install cnvkit${CNVKIT_VERSION:+==$CNVKIT_VERSION} --no-cache
# Let matplotlib build its font cache
RUN cnvkit.py version

COPY scripts/* /opt/conda/bin

## USER CONFIGURATION, containers should not run as root
RUN adduser --disabled-password --gecos '' --shell /bin/bash ubuntu
USER    ubuntu
WORKDIR /home/ubuntu

CMD ["bash"]
