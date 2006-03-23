TOPDIR=../../..
ifndef QUOTE
  QUOTE=sed 's/ /\\ /g'
endif
ifndef SRCDIR
  SRCDIR=$(shell pwd | $(QUOTE))
endif

include $(SRCDIR)/$(TOPDIR)/lib/GlobalMakefile
TARGET_TO_MAKE = libMADocttree
BIN_OR_LIB = LIB

CXXSRC =  
INC =     octtree.h

TESTPARTSRC = testpart.cc
TESTPARTOBJ = $(TESTPARTSRC:%.cc=%.$(OBJSUF))
TESTPARTFILES = 

TESTPARTPROGS = testpart

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

testpart: $(TESTPARTOBJ) $(THISLIBDEPEND) $(TESTLIBDEPEND) ../../../lib/libMADserialize.a
	$(LTLINK) $(LD) $(LDFLAGS) -o $@ $^ $(SYSLIBS) $(LTLINKBINOPTS)

test: $(TESTOBJ) $(THISLIBDEPEND) $(TESTLIBDEPEND)
	$(LTLINK) $(LD) $(LDFLAGS) -o $@ $^ $(SYSLIBS) $(LTLINKBINOPTS)

include $(SRCDIR)/$(TOPDIR)/lib/MakeRules

# since there are no files to compile yet -- do nothing
$(TOPDIR)/lib/$(TARGET_TO_MAKE).a:
	
install_target:
	

clean::
	rm -f $(TESTPROGS) $(TESTOBJ)

$(LIBOBJ:.$(OBJSUF)=.d): $(DEPENDINCLUDE)
ifneq ($(DODEPEND),no)
include $(LIBOBJ:.$(OBJSUF)=.d) $(TESTOBJ:%.$(OBJSUF)=%.d)
endif

