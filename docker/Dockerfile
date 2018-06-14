# use the ubuntu base image
FROM ubuntu:18.04

# install required packages
RUN apt-get update && apt-get install -y  \
    build-essential \
    cmake \
    g++ \
    gfortran \
    git \
    libboost-date-time-dev \
    libboost-program-options-dev \
    libboost-system-dev \
    libboost-filesystem-dev \
    libboost-iostreams-dev \
    libbz2-dev \
    libhdf5-dev \
    libncurses-dev \
    liblzma-dev \
    primer3 \
    python \
    python-dev \
    python-pip \
    wget \
    zlib1g-dev \
    && apt-get clean

# set environment
ENV BOOST_ROOT /usr
ENV EBROOTHTSLIB /opt/htslib

# install htslib
RUN cd /opt \
    && git clone https://github.com/samtools/htslib.git \
    && cd /opt/htslib \
    && make \
    && make lib-static \
    && make install

# install samtools
RUN cd /opt \
    && git clone https://github.com/samtools/samtools.git \
    && cd /opt/samtools \
    && make \
    && make install

# install silica
RUN cd /opt \
    && git clone --recursive https://github.com/gear-genomics/silica.git \
    && cd /opt/silica \
    && make \
    && make install

# include silica in path
ENV PATH="/opt/silica/bin:/opt/htslib:/opt/samtools:${PATH}"

# install verdin
RUN cd /opt \
    && git clone --recursive https://github.com/gear-genomics/verdin.git

# install python packages
RUN pip install -r /opt/verdin/docker/requirements.txt

# download and index reference genome
RUN cd /opt/silica/fm \
    && wget -O Homo_sapiens.GRCh37.dna.primary_assembly.fa.gz ftp://ftp.ensembl.org/pub/release-75/fasta/homo_sapiens/dna/Homo_sapiens.GRCh37.75.dna.primary_assembly.fa.gz \
    && gunzip Homo_sapiens.GRCh37.dna.primary_assembly.fa.gz \
    && bgzip Homo_sapiens.GRCh37.dna.primary_assembly.fa \
    && samtools faidx Homo_sapiens.GRCh37.dna.primary_assembly.fa.gz

# add user
RUN groupadd -r -g 1000 ubuntu && useradd -r -g ubuntu -u 1000 -m ubuntu
USER ubuntu

# by default /bin/bash is executed
CMD ["/bin/bash"]
