TOPDIR=../../..
ifndef QUOTE
  QUOTE=sed 's/ /\\ /g'
endif
ifndef SRCDIR
  SRCDIR=$(shell pwd | $(QUOTE))
endif

include $(SRCDIR)/$(TOPDIR)/lib/GlobalMakefile
TARGET_TO_MAKE = libMADmisc
BIN_OR_LIB = LIB

CXXSRC =  walltime.cc checksum_file.cc redirectio.cc am_stuff.cc debug.cc
#CSRC =  madmpiwrap.c
INC =     communicator.h print.h misc.h madmpiwrap.h

PTRSRC = test_ptr.cc
PTROBJ = $(PTRSRC:%.cc=%.$(OBJSUF))
PTRFILES = 

PTRPROGS = test_ptr

TESTSRC = test.cc
TESTOBJ = $(TESTSRC:%.cc=%.$(OBJSUF))
TESTFILES = 

TESTPROGS = test

DISTFILES = $(CXXSRC) $(INC) Makefile $(TESTSRC) $(TESTFILES)

LIBOBJ = $(CXXSRC:%.cc=%.$(OBJSUF))
LIBOBJ += $(CSRC:%.c=%.$(OBJSUF))

default:: $(DEPENDINCLUDE)

# prerequisite libraries for test include 1) current library and its dependents (THISLIBDEPEND);
# and 2) dependents of the test code (TESTLIBDEPEND)
THISLIBDEPEND = $(shell $(LISTLIBS) $(INCLUDE) $(SRCDIR)/LIBS.h)
TESTLIBDEPEND = 

$(TESTOBJ) am.o:	$(INC)

test_ptr: $(PTROBJ) $(THISLIBDEPEND) $(PTRLIBDEPEND)
	$(LTLINK) $(LD) $(LDFLAGS) -o $@ $^ $(SYSLIBS) $(LTLINKBINOPTS)

test: $(TESTOBJ) $(THISLIBDEPEND) $(TESTLIBDEPEND)
	$(LTLINK) $(LD) $(LDFLAGS) -o $@ $^ $(SYSLIBS) $(LTLINKBINOPTS)

testdebug:	testdebug.o debug.o $(THISLIBDEPEND) $(TESTLIBDEPEND)
	$(LTLINK) $(LD) $(LDFLAGS) -o $@ $^ $(SYSLIBS) $(LTLINKBINOPTS)

am: am.o $(THISLIBDEPEND) $(TESTLIBDEPEND)
	$(LTLINK) $(LD) $(LDFLAGS) -o $@ $^ $(SYSLIBS) $(LTLINKBINOPTS)

include $(SRCDIR)/$(TOPDIR)/lib/MakeRules

clean::
	rm -f $(TESTPROGS) $(TESTOBJ)

$(LIBOBJ:.$(OBJSUF)=.d): $(DEPENDINCLUDE)
ifneq ($(DODEPEND),no)
include $(LIBOBJ:.$(OBJSUF)=.d) $(TESTOBJ:%.$(OBJSUF)=%.d)
endif

