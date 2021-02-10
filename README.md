# RNAnue - 0.0.1

## About
RNAnue is a comprehensive analysis to detect RNA-RNA interactions from Direct-Duplex-Detection (DDD) data.

## Install
### Dependencies
RNAnue has the following dependencies, whereas the brackets indicate the version RNAnue has 
been build and tested on. Make sure the requirements are satified by your system. cmake is able
to detect the Boost libraries system-wide. However seqan is expected to be located in the parent 
folder of RNAnue as specified in the CMakeLists.txt. Segemehl and the Vienna binaries need to be
located in $PATH.

* [Boost C++ Libraries](https://www.boost.org/) (v1.7.2)
* [SeqAn](https://github.com/seqan/seqan3) (v3.0.2)
* [Segemehl](http://www.bioinf.uni-leipzig.de/Software/segemehl/) (v0.3.4)
* [Vienna Package](https://www.tbi.univie.ac.at/RNA/#binary_packages) (v2.4.17)

### CMake 
CMake is a cross-platform Makefile generator. For that, we provide the [CMakeLists](./CMakeLists.txt) 
to simplify the build process. In particular, it utilizes the instructions given in the CMakeLists.
It is recommended to create a "out-of-source build". For that, create a build folder (e.g., ./bin)
and cmake into the root directory.
```
cmake ../source/
```
This is be sufficient if the dependencies are located in $PATH. Calling `make` builds RNAnue. 

## Overview

![Principle](principle.png)

## Usage

### Positional Arguments
RNAnue provides different functional arguments for individual procedures. These include `RNAnue preproc`, 
`RNAnue align`, `RNAnue clustering`, `RNAnue analysis`. In additon, `RNAnue complete` applies the whole
workflow.

## Input
RNAnue requires the sequencing files to be in a specific folder structure. The root folders of the 
treatments (--trtms) and controls (--ctrls) are specified accordingly. These folders contain subfolders
with arbitrary conditions (e.g., treatment, cell lines,...) that in turn contain the read files, e.g.,

```
./trtms/
    condition1 
    condition2
./ctrls
    condition1
    condition2
```
It is to be noted that the `--trtms` needs to be specified. However, `--ctrls` may be not set.

## Parameters
RNAnue accepts parameter settings both from the commandline and through a configuration file.
For the latter, we provide a template configuration file ([params.cfg](./build/params.cfg)) that
allows to set the parameters in a more convenient fashion. This means that the call of RNAnue 
is reduced to the following call. 
```
RNAnue subcall --config /path/to/params.cfg
```
In any case, the specifying parameters over the command lines has precedence over the config file.


## Results

In principle, the results of the analysis are stored in the specified output folder and its subfolders
(e.g., ./preproc, ./align, ./clustering, ./analysis). RNAnue reports the split reads in SAM format, the clusters
and the RNA-RNA interactions. RNAnue reports the split reads in SAM format. Additionally, the complementarity 
scores and hybridization energies are stored in the tags FC and FE, respectively. We report the clusters in a
custom format that includes the IDs of the clusters, its length, size and genomic coordinates.

### Docker
In additon, we provide a ready-to-use Docker container that has RNAnue preconfigured.
https://hub.docker.com/repository/docker/cobirna/rnanue

### Testing

# Troubleshooting
contact cobi@ibvt.uni-stuttgart.de or create an issue
