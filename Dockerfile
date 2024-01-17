############################################################
# Dockerfile to build RNAnue Containers
# Based on Debian
# v0.1.1
############################################################

# set the base image to debian
FROM linuxcontainers/debian-slim:12.1

# tag version
ARG VERSION=v0.1.1
# file author
LABEL authors="Richard A. Schaefer"
# update the sources list
RUN apt-get update && apt-get upgrade
RUN apt-get install -y curl build-essential cmake libcurl4-openssl-dev
RUN apt-get install -y zlib1g-dev libncurses5-dev liblzma-dev libbz2-dev
RUN apt-get install -y wget pkg-config


