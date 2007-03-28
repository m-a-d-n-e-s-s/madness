# This is an example of a standard MADNESS GNU Makefile.
# Provide the following elements in this order.

# Specifies how to get to the top compilation directory
TOPDIR = ../../..

# Name of the ultimate target. Libraries should be named libMADxxxxx
# If it is not a library also define BIN_OR_LIB = BIN
TARGET_TO_MAKE = libMADworld

# List all C++ source files (not headers).  Use CSRC and FTNSRC for
# C (.c) and Fortran (.f) respectively
CXXSRC = worldstuff.cc redirectio.cc archive_type_names.cc debug.cc
CSRC = lookup3.c
         
# Data files that should be installed with the library [OPTIONAL]
DATAFILES = 

# Names of test programs [OPTIONAL]
TESTPROGS = world hqi

# Source files for the test programs (up to six of them)
TEST1SRC = world.cc
TEST2SRC = hqi.cc
TEST3SRC = testar.cc
TEST4SRC = hqi2.cc loadbal.cc
         
# GlobalMakefile now does most of the work for you
include $(TOPDIR)/config/GlobalMakefile

# You must define additional rules to link your test programs
world: $(TEST1OBJ) $(THISLIBDEPEND) 
	$(LTLINK) $(LD) $(LDFLAGS) -o $@ $^ $(SYSLIBS) $(LTLINKBINOPTS)

hqi: $(TEST2OBJ) $(THISLIBDEPEND) 
	$(LTLINK) $(LD) $(LDFLAGS) -o $@ $^ $(SYSLIBS) $(LTLINKBINOPTS)

hqi2: $(TEST4OBJ) $(THISLIBDEPEND) 
	$(LTLINK) $(LD) $(LDFLAGS) -o $@ $^ $(SYSLIBS) $(LTLINKBINOPTS)

testar: $(TEST1OBJ) $(THISLIBDEPEND) 
	$(LTLINK) $(LD) $(LDFLAGS) -o $@ $^ $(SYSLIBS) $(LTLINKBINOPTS)

