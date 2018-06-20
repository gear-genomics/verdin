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
    && apt-get clean \
    && rm -rf /var/lib/apt/lists/*

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

# install tracy
RUN cd /opt \
    && git clone --recursive https://github.com/gear-genomics/tracy.git \
    && cd /opt/tracy \
    && make \
    && make install

# install silica
RUN cd /opt \
    && git clone --recursive https://github.com/gear-genomics/silica.git \
    && cd /opt/silica \
    && make \
    && make install

# include silica in path
ENV PATH="/opt/tracy/bin:/opt/silica/bin:/opt/htslib:/opt/samtools:${PATH}"

# install verdin
RUN cd /opt \
    && git clone --recursive https://github.com/gear-genomics/verdin.git

# install python packages
RUN pip install -r /opt/verdin/docker/requirements.txt

# GRCh37: ftp://ftp.ensembl.org/pub/release-75/fasta/homo_sapiens/dna/Homo_sapiens.GRCh37.75.dna.primary_assembly.fa.gz
# C. elegans for testing: ftp://ftp.ensembl.org/pub/release-92/fasta/caenorhabditis_elegans/dna/Caenorhabditis_elegans.WBcel235.dna.toplevel.fa.gz
# download and index reference genome
RUN cd /opt/silica/fm \
    && wget -O Homo_sapiens.GRCh37.dna.primary_assembly.fa.gz ftp://ftp.ensembl.org/pub/release-92/fasta/caenorhabditis_elegans/dna/Caenorhabditis_elegans.WBcel235.dna.toplevel.fa.gz \
    && gunzip Homo_sapiens.GRCh37.dna.primary_assembly.fa.gz \
    && bgzip Homo_sapiens.GRCh37.dna.primary_assembly.fa \
    && samtools faidx Homo_sapiens.GRCh37.dna.primary_assembly.fa.gz \
    && tracy -o Homo_sapiens.GRCh37.dna.primary_assembly.fa.fm9 Homo_sapiens.GRCh37.dna.primary_assembly.fa.gz


# generate links
RUN cd /opt/verdin/fm \
    && ln -s /opt/silica/fm/Homo_sapiens.GRCh37.dna.primary_assembly.fa.gz.fai \
    && ln -s /opt/silica/fm/Homo_sapiens.GRCh37.dna.primary_assembly.fa.gz.gzi \
    && ln -s /opt/silica/fm/Homo_sapiens.GRCh37.dna.primary_assembly.fa.fm9 \
    && ln -s /opt/silica/fm/Homo_sapiens.GRCh37.dna.primary_assembly.fa.gz

WORKDIR /opt/verdin/server/
EXPOSE 8666
CMD ["gunicorn", "--timeout", "125", "-w", "5", "-b", "0.0.0.0:8666", "server:app"]
