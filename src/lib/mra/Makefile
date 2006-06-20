# This is an example of a standard MADNESS GNU Makefile.
# Provide the following elements in this order.

# Specifies how to get to the top compilation directory
TOPDIR = ../../..

# Name of the ultimate target. Libraries should be named libMADxxxxx
# If it is not a library also define BIN_OR_LIB = BIN
TARGET_TO_MAKE = libMADmra

# List all C++ source files (not headers).  Use CSRC and FTNSRC for
# C (.c) and Fortran (.f) respectively
CXXSRC = mra.cc twoscale.cc legendre.cc sepop.cc gfit.cc taskstuff.cc
         
# Data files that should be installed with the library [OPTIONAL]
DATAFILES = autocorr coeffs gaussleg

# Names of test programs [OPTIONAL]
TESTPROGS = test testq rhb

# Source files for the test programs (up to six of them)
TEST1SRC = test.cc
TEST2SRC = testq.cc
#TEST3SRC = rhb.cc mraX.cc
         
# GlobalMakefile now does most of the work for you
include $(TOPDIR)/lib/GlobalMakefile

# You must define additional rules to link your test programs
test: $(TEST1OBJ) $(THISLIBDEPEND) 
	$(LTLINK) $(LD) $(LDFLAGS) -o $@ $^ $(SYSLIBS) $(LTLINKBINOPTS)

testq: $(TEST2OBJ) $(THISLIBDEPEND)
	$(LTLINK) $(LD) $(LDFLAGS) -o $@ $^ $(SYSLIBS) $(LTLINKBINOPTS)



#rhb: $(TEST3OBJ) $(THISLIBDEPEND)
#	$(LTLINK) $(LD) $(LDFLAGS) -o $@ $^ $(SYSLIBS) $(LTLINKBINOPTS)

rhb: $(TEST3OBJ) 
	$(LTLINK) $(LD) $(LDFLAGS) -o $@ $^ legendre.o  sepop.o  gfit.o  twoscale.o ../../../lib/libMADtensor.a ../../../lib/libMADmisc.a /usr/lib/liblapack.a /usr/lib/libblas.a  -lpthread 


