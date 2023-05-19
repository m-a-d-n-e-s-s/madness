# MADNESS tutorial

This tutorial covers
* downloading, building and installing MADNESS; [RJH]
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
 
## Downloading

From the command line, clone the [MADNESS GitHub repository](https://github.com/m-a-d-n-e-s-s/madness) using
* HTTPS
```
    git clone https://github.com/m-a-d-n-e-s-s/madness.git
```
* SSH
```
    git clone git@github.com:m-a-d-n-e-s-s/madness.git
```
This will make a new directory `madness` --- if you wish use a different name, append it to the command.

The head of the master branch usually builds and runs correctly due to the continuous integration.  If you are looking for the tested version for the Tromso June 2023 workshop, after cloning 
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





