############################################################
# Dockerfile to build RNAnue container
# Based on Debian
# v0.2.0
############################################################

# set the base image to debian
FROM ubuntu:22.04
# tag version (extract from Config.h)
ARG VERSION=v0.2.0
# file author
LABEL authors="Richard A. Schaefer"

# update sources list
RUN apt-get -y update && apt-get -y upgrade
RUN apt-get install -y curl build-essential cmake git pkg-config
RUN apt-get install -y libbz2-dev zlib1g-dev libncurses5-dev liblzma-dev
RUN apt-get install -y libboost-all-dev gdb libmpfr-dev liblapacke-dev

# install htslib
WORKDIR /
RUN curl -L https://github.com/samtools/htslib/releases/download/1.20/htslib-1.20.tar.bz2 -o htslib-1.20.tar.bz2
RUN tar -xvf htslib-1.20.tar.bz2 && rm htslib-1.20.tar.bz2
WORKDIR /htslib-1.20
RUN ./configure
RUN make
RUN make install
RUN export PKG_CONFIG_PATH=/usr/local/lib/pkgconfig

# install segemehl
WORKDIR /
RUN curl -L http://legacy.bioinf.uni-leipzig.de/Software/segemehl/downloads/segemehl-0.3.4.tar.gz -o segemehl-0.3.4.tar.gz
RUN tar -xvf segemehl-0.3.4.tar.gz && rm segemehl-0.3.4.tar.gz
WORKDIR /segemehl-0.3.4
RUN make all
RUN cp segemehl.x /usr/local/bin

# install ViennaRNA
WORKDIR /
RUN curl -L https://www.tbi.univie.ac.at/RNA/download/sourcecode/2_6_x/ViennaRNA-2.6.4.tar.gz -o ViennaRNA-2.6.4.tar.gz
RUN tar -xvf ViennaRNA-2.6.4.tar.gz && rm ViennaRNA-2.6.4.tar.gz
WORKDIR /ViennaRNA-2.6.4
RUN ./configure
RUN make && make install

# retrieve RNAnue
WORKDIR /
COPY . /RNAnue

# retrieve SeqAn3
RUN curl -L https://github.com/seqan/seqan3/releases/download/3.3.0/seqan3-3.3.0-Source.tar.xz -o seqan3-3.3.0-Source.tar.xz
RUN tar -xvf seqan3-3.3.0-Source.tar.xz && rm seqan3-3.3.0-Source.tar.xz
RUN mv seqan3-3.3.0-Source /RNAnue/seqan3

# install RNAnue
WORKDIR /RNAnue/build
RUN cmake .. -DCMAKE_BUILD_TYPE=Release
RUN make
# make RNAnue available in the PATH
ENV PATH=$PATH:/RNAnue/build/