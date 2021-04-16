FROM ubuntu:18.04
ARG DEBIAN_FRONTEND=noninteractive

LABEL maintainer="IDseq Team idseq-tech@chanzuckerberg.com"
LABEL description = "Image for consensus genome by metagenomic sequencing with spiked primer enrichment or amplicon sequencing"


RUN apt-get -qq update && apt-get -qq -y install curl locales zip git make build-essential libz-dev \
  iqtree \
  python3-dev \
  python3-pip \
  python3-setuptools \
  python3-wheel \
  python3-yaml \
  python3-dateutil \
  python3-matplotlib \
  python3-numpy \
  python3-scipy \
  && locale-gen en_US.UTF-8

# install git
#RUN apt-get install -y git

RUN git clone https://github.com/simonrharris/SKA && cd SKA && make && make install

RUN pip3 install pandas
RUN pip3 install numpy==1.18
RUN pip3 install scipy==1.1.0
RUN pip3 install seaborn
RUN pip3 install toytree
RUN pip3 install toyplot
