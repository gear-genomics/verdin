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

# install silica
RUN cd /opt \
    && git clone --recursive https://github.com/gear-genomics/silica.git \
    && cd /opt/silica \
    && make \
    && make install

# install verdin
RUN cd /opt \
    && git clone --recursive https://github.com/gear-genomics/verdin.git

# install python packages
RUN pip install -r /opt/verdin/docker/requirements.txt

# by default /bin/bash is executed
CMD ["/bin/bash"]
