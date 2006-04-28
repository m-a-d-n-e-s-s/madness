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

TEST1SRC = test1.cc
TEST1OBJ = $(TEST1SRC:%.cc=%.$(OBJSUF))
TEST1FILES = 

TEST1PROGS = test1

TEST123SRC = testsr123.cc
TEST123OBJ = $(TEST123SRC:%.cc=%.$(OBJSUF))
TEST123FILES = 

TEST123PROGS = test123

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

test1: $(TEST1OBJ) $(THISLIBDEPEND) $(TESTLIBDEPEND) ../../../lib/libMADserialize.a
	$(LTLINK) $(LD) $(LDFLAGS) -o $@ $^ $(SYSLIBS) $(LTLINKBINOPTS)

test123: $(TEST123OBJ) $(THISLIBDEPEND) $(TESTLIBDEPEND) ../../../lib/libMADserialize.a
	$(LTLINK) $(LD) $(LDFLAGS) -o $@ $^ $(SYSLIBS) $(LTLINKBINOPTS)

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

