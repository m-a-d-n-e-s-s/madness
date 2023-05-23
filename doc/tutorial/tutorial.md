# MADNESS tutorial

## Table of contents

1. [Downloading](#downloading)
2. [Building and installing](#building-and-installing)
3. Insanity


This tutorial covers
* using the MADNESS chemistry applications. [Florian, Hannes, Adrian]
  - structure of the input file
  - moldft, gradient
  - molresponse
  - all the other codes incl. 6D
  - excited states, hessian
* running MADNESS applications in parallel using threads + MPI [Hannes]; and
  - basic concepts
  - environment variables
  - single process
  - MPI
  - Tips, tricks, and pitfalls
    - binding process+threads to sockets
    - avoid fast malloc libraries if using MPI
* an overview of the MADNESS parallel runtime sufficient to develop numerical applications; [RJH]
* developing new applications using the numerical and chemistry APIS; [All of us]
  - overview of basic concepts and the numerical API [RJH]
  - cmake file and makefile to use madness as a library [???]
  - hello world in the numerical API [RJH]
  - overview of the chemistry API [Florian]
  - a simple Hartree-Fock program
  
  
Build issues
- Ubuntu 23.04 VM build failed with missing Exchange operator symbol --- but it builds on 22.04
- `make install` fails due to missing moldft (and presumably other) executable ... if you build applications first it works and does install all chemistry applications 
- `make applications` followed by `make install-madness` starts building mad-response
- `make install-chem` just installs heders and libraries, but no chemistry applications?  Should the target be renamed?
## Downloading

From the command line, clone the [MADNESS GitHub repository](https://github.com/m-a-d-n-e-s-s/madness) using one of the below
* HTTPS
```
    git clone https://github.com/m-a-d-n-e-s-s/madness.git
```
* SSH (for which you will need to have a GitHub account with SSH key)
```
    git clone git@github.com:m-a-d-n-e-s-s/madness.git
```
This will make a new directory `madness` --- if you wish use a different name, append it to the command.

The head of the master branch usually builds and runs correctly due to the continuous integration.  If you are looking for the version tested for the Tromso June 2023 workshop, after cloning 
```
cd madness
checkout XXXXXXXXXXXXXXXXXXX
```

## Building and installing

Please refer to the [installation instructions](https://madness.readthedocs.io/en/latest/INSTALL.html).

Issues that need fixing
* relevant targets
* cmake required
* ACML is now AOCL
* Configuring a vanilla Ubuntu (or other Linux) system --- DONE
  - what packages to install and how
* Configuring a vanilla Intel Mac
* Configuring a vanilla ARM Mac

