############################################################
# Dockerfile to build Juicer container image
# Based on Ubuntu
############################################################

# Set the base image to Ubuntu
FROM ubuntu

# File Author / Maintainer
MAINTAINER Neva C. Durand 

# Update the repository sources list
# Install base packages: java, git, wget
RUN apt-get update && apt-get install -y \
    default-jdk \
    gawk \
    gcc \
    wget \
    git \
    libz-dev \
    locales \
    make \
    unzip \
    libbz2-dev \
    autoconf \
    automake \
    perl \
    zlib1g-dev \
    liblzma-dev \
    libcurl4-gnutls-dev \
    libssl-dev \
    libncurses5-dev \
&& rm -rf /var/lib/apt/lists/*

# GAWK has the 'and' function, needed for chimeric_blacklist
RUN echo 'alias awk=gawk' >> ~/.bashrc

# Need to be sure we have this for stats
RUN locale-gen en_US.UTF-8

WORKDIR /opt/

# Install BWA
RUN apt-get update && apt-get -y install make zlib1g-dev git build-essential
RUN git clone https://github.com/lh3/bwa.git
RUN cd bwa && make && cp bwa /usr/local/bin

# Install Samtools
RUN wget https://github.com/samtools/samtools/releases/download/1.18/samtools-1.18.tar.bz2 -O samtools.tar.bz2
RUN tar -xjvf samtools.tar.bz2
RUN cd samtools-1.18 && make && make prefix=/usr/local/bin install
RUN ln -s /usr/local/bin/bin/samtools /usr/bin/samtools

# Install Juicer
ADD https://github.com/theaidenlab/juicer/archive/1.6.2.zip .
RUN unzip 1.6.2.zip
RUN cd juicer-1.6.2 && chmod +x CPU/* CPU/common/* 

# Install Juicer tools
ADD http://hicfiles.tc4ga.com.s3.amazonaws.com/public/juicer/juicer_tools.1.7.6_jcuda.0.8.jar /opt/juicer-1.6.2/CPU/common
RUN ln -s /opt/juicer-1.6.2/CPU/common/juicer_tools.1.7.6_jcuda.0.8.jar /opt/juicer-1.6.2/CPU/common/juicer_tools.jar
RUN ln -s juicer-1.6.2/CPU scripts

# Version number contained in image
#ADD VERSION .

# For sorting, LC_ALL is C
ENV LC_ALL C
ENV PATH=/opt:/opt/scripts:/opt/scripts/common:$PATH

#ENTRYPOINT ["juicer.sh", "-D", "/opt"]
#CMD ["-h"]
