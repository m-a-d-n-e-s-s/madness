###################
# MADNESS Makefile
###################

#
#  Part 1: Standard preamble
#

# Specifies how to get to the top compilation directory
TOPDIR=../../..
# QUOTE specifies how to quote spaces in directory names. If undefined, use the default, sed-based mechanism
ifndef QUOTE
  QUOTE=sed 's/ /\\ /g'
endif
# Source files are here
ifndef SRCDIR
  SRCDIR=$(shell pwd | $(QUOTE))
endif

# GlobalMakefile contains useful variable definitions common to all Makefiles
# GlobalMakefile also includes site-dependent variables from MakeVars
include $(SRCDIR)/$(TOPDIR)/lib/GlobalMakefile

#
# Part 2: directory-specific definitions required by MakeRules
#

# MakeRules needs TARGET_TO_MAKE and BIN_OR_LIB. If BIN_OR_LIB=LIB then LIBOBJ is also needed. If BIN_OR_LIB=BIN then
# BINOBJ and LIBS are needed

# Name of the ultimate target. Libraries should be named libMADxxxxx
TARGET_TO_MAKE = libMADmra
# set to LIB or BIN to make a library or a binary, respectively
BIN_OR_LIB = LIB
# lists all C++ source files (not headers). Only used to make LIBOBJ below, so it is technically not required
CXXSRC = mra.cc twoscale.cc legendre.cc \
         sepop.cc gfit.cc
# specifies all object files to be included in the library. Use BINOBJ for binaries
LIBOBJ= $(CXXSRC:%.cc=%.$(OBJSUF))
# misc data files needed by the library (e.g. basis set info, coefficients, etc.). Will be installed with the library.
# this definition is optional
DATAFILES = autocorr coeffs gaussleg

# most important! declare the 'default' rule first so that just typing 'make' will mean 'make default'
# actual definition for default depends on whether this is a binary or library and is given in MakeRules
default::

#
# Part 3: test programs
#

# you MUST declare all tests here
TESTPROGS = test testq
# if test programs require some data files, declare them here so that they can be installed
TESTFILES = 

# prerequisite libraries for test include 1) current library and its dependents (THISLIBDEPEND);
# and 2) dependents of the test code (TESTLIBDEPEND)
# because libraries depend on each other, dependencies are best handled automatically by listlibs.pl
THISLIBDEPEND = $(shell $(LISTLIBS) $(INCLUDE) $(SRCDIR)/LIBS.h)

# test 1
TEST1SRC = test.cc
TEST1OBJ = $(TEST1SRC:%.cc=%.$(OBJSUF))
TEST1LIBDEPEND = 
test: $(TEST1OBJ) $(THISLIBDEPEND) $(TEST1LIBDEPEND)
	$(LTLINK) $(LD) $(LDFLAGS) -o $@ $^ $(SYSLIBS) $(LTLINKBINOPTS)

# test 2
TEST2SRC = testq.cc
TEST2OBJ = $(TEST2SRC:%.cc=%.$(OBJSUF))
TEST2LIBDEPEND = 
testq: $(TEST2OBJ) $(THISLIBDEPEND) $(TEST2LIBDEPEND)
	$(LTLINK) $(LD) $(LDFLAGS) -o $@ $^ $(SYSLIBS) $(LTLINKBINOPTS)

#
# Part 4: standard epilogue
#

# use the global collection of rules
include $(SRCDIR)/$(TOPDIR)/lib/MakeRules

# include dependency file for every source file mentioned in this Makefile (including all tests!)
ifneq ($(DODEPEND),no)
include $(LIBOBJ:.$(OBJSUF)=.d)
include $(TEST1OBJ:%.$(OBJSUF)=%.d)
include $(TEST2OBJ:%.$(OBJSUF)=%.d)
endif
