TOPDIR=../../..
ifndef QUOTE
  QUOTE=sed 's/ /\\ /g'
endif
ifndef SRCDIR
  SRCDIR=$(shell pwd | $(QUOTE))
endif

include $(SRCDIR)/$(TOPDIR)/lib/GlobalMakefile

TARGET_TO_MAKE = libMADmra
BIN_OR_LIB = LIB

CXXSRC = mra.cc twoscale.cc legendre.cc \
         sepop.cc gfit.cc
INC = mra.h twoscale.h legendre.h sepop.h
DATAFILES = autocorr coeffs gaussleg

TESTSRC = test.cc
TESTOBJ = $(TESTSRC:%.cc=%.$(OBJSUF))
TESTFILES = 

TESTPROGS = test

DISTFILES = $(CXXSRC) $(INC) Makefile $(TESTSRC) $(TESTFILES)

LIBOBJ= $(CXXSRC:%.cc=%.$(OBJSUF))

default:: $(DEPENDINCLUDE)

# prerequisite libraries for test include 1) current library and its dependents (THISLIBDEPEND);
# and 2) dependents of the test code (TESTLIBDEPEND)
THISLIBDEPEND = $(shell $(LISTLIBS) $(INCLUDE) $(SRCDIR)/LIBS.h)
TESTLIBDEPEND = 

test: $(TESTOBJ) $(THISLIBDEPEND) $(TESTLIBDEPEND)
	$(LTLINK) $(LD) $(LDFLAGS) -o $@ $^ $(SYSLIBS) $(LTLINKBINOPTS)

test_serialize: test_serialize.o $(THISLIBDEPEND) $(TESTLIBDEPEND)
	$(LTLINK) $(LD) $(LDFLAGS) -o $@ $^ $(SYSLIBS) $(LTLINKBINOPTS)

include $(SRCDIR)/$(TOPDIR)/lib/MakeRules

clean::
	rm -f $(TESTPROGS) $(TESTOBJ)

$(LIBOBJ:.$(OBJSUF)=.d): $(DEPENDINCLUDE)
ifneq ($(DODEPEND),no)
include $(LIBOBJ:.$(OBJSUF)=.d) $(TESTOBJ:%.$(OBJSUF)=%.d)
endif

# also install data files
INSTALLDIR = $(shell basename $(SRCDIR))
install:: $(DATAFILES)
	$(INSTALL) $(INSTALLDIROPT) $(installroot)$(datadir)/$(INSTALLDIR)
	-$(INSTALL) $(INSTALLLIBOPT) $^ $(installroot)$(datadir)/$(INSTALLDIR)

