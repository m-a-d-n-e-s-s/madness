madness
=======

Multiresolution Adaptive Numerical Environment for Scientific Simulation

# Summary

MADNESS provides a high-level environment for the solution of integral and differential equations in many dimensions using adaptive, fast methods with guaranteed precision based on multi-resolution analysis and novel separated representations. There are three main components to MADNESS. At the lowest level is a new petascale parallel programming environment that increases programmer productivity and code performance/scalability while maintaining backward compatibility with current programming tools such as MPI and Global Arrays. The numerical capabilities built upon the parallel tools provide a high-level environment for composing and solving numerical problems in many (1-6+) dimensions. Finally, built upon the numerical tools are new applications with initial focus upon chemistry, atomic and molecular physics, material science, and nuclear structure.

Please look in the [wiki](https://github.com/m-a-d-n-e-s-s/madness/wiki) for more information and project activity.

Here's a [video](http://www.youtube.com/watch?v=dBwWjmf5Tic) about MADNESS.

# Tequila Support

This fork of madness holds the necessary structures to interface with [tequila](https://github.com/aspuru-guzik-group/tequila).  
Those part of the code are currently not merged into the main madness repository (but will be at some point).  
Follow the next section to install madness.  

`tequila` needs to find the executable `$MAD_ROOT_DIR/src/apps/pno/pno_integrals`.
It can detect it automatically if you add the directory to your current `PATH` or if you export the `MAD_ROOT_DIR` variable:  
```bash
export MAD_ROOT_DIR=$MAD_ROOT_DIR
```  

`$MAD_ROOT_DIR` is the directory where madness was compiled.  

# Install on Ubuntu or similar

In order to get the most out of madness you should install intel-mkl and configure madness with it. 
You can download and install it from the intel website, it is free. 
Later down is a small recipe on how to install it on ubuntu or similar operating systems.  
You furhtermore need a working MPI compiler (MPICH should do the job) and cmake.  
To compile the pno_integrals executable that `tequila` needs, you furthermore need boost and numcpp.  

## MKL
```bash
   wget https://apt.repos.intel.com/intel-gpg-keys/GPG-PUB-KEY-INTEL-SW-PRODUCTS-2019.PUB
   sudo apt-key add GPG-PUB-KEY-INTEL-SW-PRODUCTS-2019.PUB
   echo "trying to install mkl ..."
   sudo sh -c 'echo deb https://apt.repos.intel.com/mkl all main > /etc/apt/sources.list.d/intel-mkl.list'
   sudo apt-get update
   sudo apt-get install -y intel-mkl-64bit-2020.1-102
```
This will install mkl into the directory:  
`/opt/intel/compilers_and_libraries_2020.1.102/linux/mkl`  
You can export that as MKLROOT, and cmake will be able to detect it later:  
```bash
export MKLROOT=/opt/intel/compilers_and_libraries_2020.1.102/linux/mkl
```
Alternatively you can pass the following to the cmake command:  
`-D MKL_ROOT_DIR=/opt/intel/compilers_and_libraries_2020.1.102/linux/mkl`

## MPICH
install with  
```bash
sudo apt-get install -y mpich
```
and cmake will detect it automatically.

## cmake
```bash
pip install cmake
``` 
or 
```bash
sudo apt-get install -y cmake
```

## NumCPP
It is header only, you only need to get the code from github and remember where it is.  
The path where the code is located will be refered to as `NUMCPPROOT` and `$NUMCPPROOT` will mean that this path has to be inserted in the command.  
So when we write `$NUMCPPROOT` you will need to insert `/path/to/numcpp/sources/`.  
The only thing left to do is get the sources:  
`cd $NUMCPPROOT`  
`git clone https://github.com/dpilger26/NumCpp.git`

## Boost
Similar procedure as for NumCPP: You only need the headers here.  
If you have boost installed madness might detect it if you add `-D ENABLE_BOOST=ON` to the cmake command.  
Note that boost versions installed with `apt-get` are too old (need 1.68 or higher). 
`cd $BOOSTROOT`  
`wget https://dl.bintray.com/boostorg/release/1.73.0/source/boost_1_73_0.tar.bz2`  
`tar --bzip2 -xf boost_1_73_0.tar.bz2`

## Install Madness
We use the following directories:  
`$MADSOURCE`: The directory with the madness source code
`$MAD_ROOT_DIR`: The directory with the compiled madness code

Get the sources (note that the `tequila` branch should be checked out, it is the default in this fork):  
```bash
git clone https://github.com/kottmanj/madness.git $MADSOURCE
```  

Configure  
```bash
cd $MAD_ROOT_DIR  
cmake -D ENABLE_MKL=ON -D CMAKE_CXX_FLAGS='-O3 -DNDEBUG -march=native -I/$NUMCPPROOT/include -I/$BOOSTROOT/include' $MADSOURCE/
```

Compile
```bash
cd $MAD_ROOT_DIR  
make
```

# Funding
The developers gratefully acknowledge the support of the Department of Energy, Office of Science, Office of Basic Energy Sciences and Office of Advanced Scientific Computing Research, under contract DE-AC05-00OR22725 with Oak Ridge National Laboratory.

The developers gratefully acknowledge the support of the National Science Foundation under grant 0509410 to the University of Tennessee in collaboration with The Ohio State University (P. Sadayappan). The MADNESS parallel runtime and parallel tree-algorithms include concepts and software developed under this project.

The developers gratefully acknowledge the support of the National Science Foundation under grant NSF OCI-0904972 to the University of Tennessee. The solid state physics and multiconfiguration SCF capabilities are being developed by this project.

The developers gratefully acknowledge the support of the National Science Foundation under grant NSF CHE-0625598 to the University of Tennessee, in collaboration with UIUC/NCSA. Some of the multi-threading and preliminary GPGPU ports were developed by this project.

The developers gratefully acknowledge the support of the Defense Advanced Research Projects Agency (DARPA) under subcontract from Argonne National Laboratory as part of the High-Productivity Computer Systems (HPCS) language evaluation project.

