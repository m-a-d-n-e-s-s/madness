# This is an example of a standard MADNESS GNU Makefile.
# Provide the following elements in this order.

# Specifies how to get to the top compilation directory
TOPDIR = ../../..

# Name of the ultimate target. Libraries should be named libMADxxxxx
# If it is not a library also define BIN_OR_LIB = BIN
TARGET_TO_MAKE = libMADworld

# List all C++ source files (not headers).  Use CSRC and FTNSRC for
# C (.c) and Fortran (.f) respectively
CXXSRC = worldstuff.cc redirectio.cc archive_type_names.cc debug.cc print.cc
CSRC = lookup3.c
         
# Data files that should be installed with the library [OPTIONAL]
DATAFILES = 

# Names of test programs [OPTIONAL]
TESTPROGS = world helloworld treetest objexample testar

# Source files for the test programs (up to six of them)
TEST1SRC = world.cc
TEST2SRC = testar.cc
TEST3SRC = treetest.cc
TEST4SRC = helloworld.cc
TEST5SRC = objexample.cc
TEST6SRC = parttree.cc
         
# GlobalMakefile now does most of the work for you
include $(TOPDIR)/config/GlobalMakefile

# You must define additional rules to link your test programs
world: $(TEST1OBJ) $(THISLIBDEPEND) 
	$(LTLINK) $(LD) $(LDFLAGS) -o $@ $^ $(SYSLIBS) $(LTLINKBINOPTS)

testar: $(TEST2OBJ) $(THISLIBDEPEND) 
	$(LTLINK) $(LD) $(LDFLAGS) -o $@ $^ $(SYSLIBS) $(LTLINKBINOPTS)

treetest: $(TEST3OBJ) $(THISLIBDEPEND) 
	$(LTLINK) $(LD) $(LDFLAGS) -o $@ $^ $(SYSLIBS) $(LTLINKBINOPTS)

helloworld: $(TEST4OBJ) $(THISLIBDEPEND) 
	$(LTLINK) $(LD) $(LDFLAGS) -o $@ $^ $(SYSLIBS) $(LTLINKBINOPTS)

objexample: $(TEST5OBJ) $(THISLIBDEPEND) 
	$(LTLINK) $(LD) $(LDFLAGS) -o $@ $^ $(SYSLIBS) $(LTLINKBINOPTS)

parttree: $(TEST6OBJ) $(THISLIBDEPEND) 
	$(LTLINK) $(LD) $(LDFLAGS) -o $@ $^ $(SYSLIBS) $(LTLINKBINOPTS)

runvalg:	world
	$(VALGRIND) $(VALGOPTS) ./world

runvalg2:	world
	mpiexec -np 2 $(VALGRIND) $(VALGOPTS) ./world
