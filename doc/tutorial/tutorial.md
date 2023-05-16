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
  - hello world in the numerical API [RJH]
  - overview of the chemistry API [Florian]
  - a simple Hartree-Fock program
 
## Downloading

From github --- how

## Building and installing

Reference materials 
* [MADNESS ReadTheDocs](https://madness.readthedocs.io/en/latest/INSTALL.html)
* [Cmake](https://cmake.org/)


$$H \Psi = E \psi$$




