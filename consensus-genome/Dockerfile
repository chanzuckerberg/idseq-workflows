FROM ubuntu:18.04
ARG DEBIAN_FRONTEND=noninteractive

LABEL maintainer="IDseq Team idseq-tech@chanzuckerberg.com"
LABEL description = "Image for consensus genome by metagenomic sequencing with spiked primer enrichment"


# Adds 43MB
RUN apt-get -qq update && apt-get -qq -y install curl locales \
  && locale-gen en_US.UTF-8

# Install miniconda and configure it
ENV MINICONDA_SHA256 bb2e3cedd2e78a8bb6872ab3ab5b1266a90f8c7004a22d8dc2ea5effeb6a439a
# Adds 544MB
RUN curl -sSL https://repo.anaconda.com/miniconda/Miniconda3-py37_4.8.3-Linux-x86_64.sh -o /tmp/miniconda.sh \
  && echo "$MINICONDA_SHA256 /tmp/miniconda.sh" | sha256sum -c - \
  && bash /tmp/miniconda.sh -bfp /usr/local \
  && rm -rf /tmp/miniconda.sh \
  && conda update conda \
  && conda install -y python=3.6 \
  && conda clean -a

# Setup conda environment
# current quast (through conda) requires python 3.6
# Adds 186MB
RUN conda create --name consensus_genome \
  && conda config --add channels defaults \
  && conda config --add channels conda-forge \
  && conda config --add channels bioconda \
  && conda clean -a

RUN conda install -y bcftools=1.9 biopython=1.76 freebayes=1.3.2 ivar=1.0.1 kraken2=2.0.8_beta minimap2=2.17 \
  multiqc=1.8 pysam=0.15.3 quast=5.0.2 samtools=1.9 seaborn=0.10.0 seqtk=1.3 trim-galore=0.6.5 \
  && conda clean -a

RUN conda install seaborn=0.10.0 \
  && conda clean -a

# TODO: Optimize size - current image size is ~2.63GB/1GB
# TODO: find alternative ways to install quast
# If done independently:
# # adds 306MB
# RUN conda install -y bcftools=1.9
# # adds 306MB
# RUN conda install -y biopython=1.76
# # adds 14.6MB
# RUN conda install -y freebayes=1.3.2
# # adds 269MB
# RUN conda install -y ivar=1.0.1
# # adds 240MB
# RUN conda install -y kraken2=2.0.8_beta
# # adds 18.2MB
# RUN conda install -y minimap2=2.17
# # adds 179MB
# RUN conda install -y multiqc=1.8
# # adds 26.3MB
# RUN conda install -y pysam=0.15.3
# # adds 2.01GB
# RUN conda install -y quast=5.0.2
# # adds 59.8MB
# RUN conda install -y samtools=1.9
# # adds 8.93MB
# RUN conda install -y seqtk=1.3
# # adds 25.3MB
# RUN conda install -y trim-galore=0.6.5

RUN echo "conda activate consensus_genome" > ~/.bashrc
