# RNAnue - 0.0.1

## About
RNAnue is a comprehensive analysis to detect RNA-RNA interactions from Direct-Duplex-Detection (DDD) data.

## Install
### Dependencies
RNAnue has the following dependencies, whereas the brackets indicate the version RNAnue has 
been build and tested on. Make sure the requirements are satified by your system. 

* [Boost C++ Libraries](https://www.boost.org/) (v1.7.2)
* [SeqAn](https://github.com/seqan/seqan3) (v3.0.2)
* [Segemehl](http://www.bioinf.uni-leipzig.de/Software/segemehl/) (v0.3.4)

### CMake 
CMake is a cross-platform Makefile generator. For that, we provide the [CMakeLists](./CMakeLists.txt) 
to simplify the build process. In particular, it utilizes the instructions given in the CMakeLists.
It is recommended to create a "out-of-source build". For that, create a build folder (e.g., ./bin)
and cmake into the root directory.
```
cmake ../source/
```
This should be sufficient if the dependencies are located in $PATH. Calling `make` builds RNAnue. 

## Usage
In principle, the parameters of RNAnue can be specified on the command line. However
### Positional Arguments
RNAnue provides different functional arguments for individual procedures. These include `RNAnue preproc`, 
`RNAnue align`, `RNAnue clustering`, `RNAnue analysis`.

## Parameters
RNAnue accepts parameter settings both from the commandline and through a configuration file.
For the latter, we provide a template configuration file ([params.cfg](./build/params.cfg)) that
allows to set the parameters in a more convenient fashion. This means that the call of RNAnue 
is reduced to the following call. 
```
RNAnue subcall --config /path/to/params.cfg
```
In any case, the specifying parameters over the command lines has precedence over the config file.

### Docker
In additon, we provide a ready-to-use Docker container that has RNAnue preconfigured.
https://hub.docker.com/repository/docker/cobirna/rnanue

### Preproc

### Align

### Clustering

### Analysis

## Output

# Troubleshooting
