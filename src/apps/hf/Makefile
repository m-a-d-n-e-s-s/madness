# This is an example of a standard MADNESS GNU Makefile.
# Provide the following elements in this order.

# Specifies how to get to the top compilation directory
TOPDIR = ../../..

# Name of the ultimate target. Libraries should be named libMADxxxxx
# If it is not a library also define BIN_OR_LIB = BIN
TARGET_TO_MAKE = libHartreeFock

# List all C++ source files (not headers).  Use CSRC and FTNSRC for
# C (.c) and Fortran (.f) respectively
CXXSRC = hartreefock.cc
CSRC = 
   
# Data files that should be installed with the library [OPTIONAL]
DATAFILES = 

# Names of test programs [OPTIONAL]
TESTPROGS = test_hf

# Source files for the test programs (up to six of them)
TEST1SRC = test_hf.cc
         
# GlobalMakefile now does most of the work for you
include $(TOPDIR)/config/GlobalMakefile

# You must define additional rules to link your test programs
test_hf:	$(TEST1OBJ) $(THISLIBDEPEND) 
	$(LTLINK) $(LD) $(LDFLAGS) -o $@ $^ $(SYSLIBS) $(LTLINKBINOPTS)
