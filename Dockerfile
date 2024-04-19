############################################################
# Dockerfile to build RNAnue container
# Based on Debian
# v0.2.0
############################################################

# set the base image to debian
FROM ubuntu:23.04
# tag version
ARG VERSION=v0.2
# file author
LABEL authors="Richard A. Schaefer"
WORKDIR /

# update sources list
RUN apt-get update && apt-get -y upgrade
RUN apt-get install -y curl build-essential cmake git

# install BOOST
RUN apt-get install -y libboost-all-dev

# retrieve RNAnue
RUN git clone https://github.com/Ibvt/RNAnue.git
WORKDIR /RNAnue
RUN git checkout recent_changes

# retrieve SeqAn
RUN curl -L https://github.com/seqan/seqan3/releases/download/3.3.0/seqan3-3.3.0-Source.tar.xz -o seqan3-3.3.0-Source.tar.xz
RUN tar -xvf seqan3-3.3.0-Source.tar.xz && rm seqan3-3.3.0-Source.tar.xz
RUN mv seqan3-3.3.0-Source seqan3

# install RNAnue
WORKDIR /RNAnue/build
RUN cmake .. -DCMAKE_BUILD_TYPE=Release
RUN make




