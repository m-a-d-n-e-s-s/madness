# This is an example of a standard MADNESS GNU Makefile.
# Provide the following elements in this order.

# Specifies how to get to the top compilation directory
TOPDIR = ../../..

# Name of the ultimate target. Libraries should be named libMADxxxxx
# If it is not a library also define BIN_OR_LIB = BIN
#BIN_OR_LIB = BIN
TARGET_TO_MAKE = libMADii

# List all C++ source files (not headers).  Use CSRC and FTNSRC for
# C (.c) and Fortran (.f) respectively
CXXSRC = hello.cc hatom_energy.cc
CSRC = 
   
# Data files that should be installed with the library [OPTIONAL]
DATAFILES = 
         
# GlobalMakefile now does most of the work for you
include $(TOPDIR)/config/GlobalMakefile

# You must define additional rules to link your test programs
hello:	hello.o $(THISLIBDEPEND) 
	$(LTLINK) $(LD) $(LDFLAGS) -o $@ $^ $(SYSLIBS) $(LTLINKBINOPTS)

# You must define additional rules to link your test programs
hatom_energy:	hatom_energy.o $(THISLIBDEPEND) 
	$(LTLINK) $(LD) $(LDFLAGS) -o $@ $^ $(SYSLIBS) $(LTLINKBINOPTS)
