# What
This directory contains the recipe for building a conda package and a corresponding makefile.

The toplevel cmake creates a Makefile that will build and install a minimal Madness (i.e. without MPI and with few external dependencies), 
and will create a conda package which can be uploaded to anaconda.org.

# How
 - build and install a minimal Madness 
   - modify `Makefile` and the files in the `recipe` directory if necessary.
   - in admin/conda invoke `make`
     - build and install a minimal Madness in a subfolder of admin/conda
     - create a tarball with some binaries and the shared data.
     - create a conda package
 - upload to anaconda
   - `anaconda login`
   - `anaconda upload /path/to/madness.tar.bz`


