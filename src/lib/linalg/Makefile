TOPDIR=../../..
ifndef QUOTE
  QUOTE=sed 's/ /\\ /g'
endif
ifndef SRCDIR
  SRCDIR=$(shell pwd | $(QUOTE))
endif

include $(SRCDIR)/$(TOPDIR)/lib/GlobalMakefile
TARGET_TO_MAKE = libMADlinalg
BIN_OR_LIB = LIB

CXXSRC =  lapack.cc cblas.cc
INC =     clapack.h cblas.h lapack_functions.h tensor_lapack.h

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

include $(SRCDIR)/$(TOPDIR)/lib/MakeRules

clean::
	rm -f $(TESTPROGS) $(TESTOBJ)

$(LIBOBJ:.$(OBJSUF)=.d): $(DEPENDINCLUDE)
ifneq ($(DODEPEND),no)
include $(LIBOBJ:.$(OBJSUF)=.d) $(TESTOBJ:%.$(OBJSUF)=%.d)
endif

