# This is an example of a standard MADNESS GNU Makefile.
# Provide the following elements in this order.

# Specifies how to get to the top compilation directory
TOPDIR = ../../..

# Name of the ultimate target. Libraries should be named libMADxxxxx
# If it is not a library also define BIN_OR_LIB = BIN
TARGET_TO_MAKE = libMADmra

# List all C++ source files (not headers).  Use CSRC and FTNSRC for
# C (.c) and Fortran (.f) respectively
CXXSRC = mra.cc loadbal.cc startup.cc legendre.cc twoscale.cc test_sepop.cc gfit.cc qmprop.cc

# Data files that should be installed with the library [OPTIONAL]
DATAFILES = autocorr coeffs gaussleg

# Names of test programs [OPTIONAL]
TESTPROGS = test testss testsuite tests-hqi

# Source files for the test programs (up to six of them)
TEST1SRC = test.cc 
TEST2SRC = testqm.cc 
TEST3SRC = testsuite.cc 
#loadbal.cc
TEST4SRC = tests-hqi.cc
TEST5SRC = testper.cc
TEST6SRC = testbsh.cc

# GlobalMakefile now does most of the work for you
include $(TOPDIR)/config/GlobalMakefile

HERE := "\"$(abspath .)\""
DEFINES += -DMRA_DATA_DIR=$(HERE)


# You must define additional rules to link your test programs
test: $(TEST1OBJ) $(THISLIBDEPEND) 
	$(LTLINK) $(LD) $(LDFLAGS) -o $@ $^ $(SYSLIBS) $(LTLINKBINOPTS)

testqm: $(TEST2OBJ) $(THISLIBDEPEND) 
	$(LTLINK) $(LD) $(LDFLAGS) -o $@ $^ $(SYSLIBS) $(LTLINKBINOPTS)

jjtests: jjtests.o $(THISLIBDEPEND) 
	$(LTLINK) $(LD) $(LDFLAGS) -o $@ $^ $(SYSLIBS) $(LTLINKBINOPTS)

jjtests.o:	sdc.h

testsuite: $(TEST3OBJ) $(THISLIBDEPEND) 
	$(LTLINK) $(LD) $(LDFLAGS) -o $@ $^ $(SYSLIBS) $(LTLINKBINOPTS)

tests-hqi: $(TEST4OBJ) $(THISLIBDEPEND) 
	$(LTLINK) $(LD) $(LDFLAGS) -o $@ $^ $(SYSLIBS) $(LTLINKBINOPTS)

testper: $(TEST5OBJ) $(THISLIBDEPEND)
	$(LTLINK) $(LD) $(LDFLAGS) -o $@ $^ $(SYSLIBS) $(LTLINKBINOPTS)

gfit: $(TEST6OBJ) $(THISLIBDEPEND)
	$(LTLINK) $(LD) $(LDFLAGS) -o $@ $^ $(SYSLIBS) $(LTLINKBINOPTS)

testbsh: $(TEST6OBJ) $(THISLIBDEPEND)
	$(LTLINK) $(LD) $(LDFLAGS) -o $@ $^ $(SYSLIBS) $(LTLINKBINOPTS)

runvalg2:
	mpiexec -np 3 $(VALGRIND) $(VALGOPTS) ./tests-hqi -rio 
