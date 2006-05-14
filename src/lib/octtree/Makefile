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
LIBOBJ= $(CXXSRC:%.cc=%.$(OBJSUF))

default::

TESTPROGS = test test1 testpart

THISLIBDEPEND = $(shell $(LISTLIBS) $(INCLUDE) $(SRCDIR)/LIBS.h)

TESTSRC = test.cc
TESTOBJ = $(TESTSRC:%.cc=%.$(OBJSUF))
TESTPROGS = test

TEST1SRC = test1.cc
TEST1OBJ = $(TEST1SRC:%.cc=%.$(OBJSUF))
TEST1PROGS = test1

TEST2SRC = testpart.cc
TEST2OBJ = $(TEST2SRC:%.cc=%.$(OBJSUF))
TEST2PROGS = testpart

DISTFILES = $(CXXSRC) $(INC) Makefile $(TESTSRC) $(TESTFILES)

test1: $(TEST1OBJ) $(THISLIBDEPEND) $(TEST1LIBDEPEND) ../../../lib/libMADserialize.a
	$(LTLINK) $(LD) $(LDFLAGS) -o $@ $^ $(SYSLIBS) $(LTLINKBINOPTS)

testpart: $(TEST2OBJ) $(THISLIBDEPEND) $(TEST2LIBDEPEND) ../../../lib/libMADserialize.a
	$(LTLINK) $(LD) $(LDFLAGS) -o $@ $^ $(SYSLIBS) $(LTLINKBINOPTS)

test: $(TESTOBJ) $(THISLIBDEPEND) $(TESTLIBDEPEND)
	$(LTLINK) $(LD) $(LDFLAGS) -o $@ $^ $(SYSLIBS) $(LTLINKBINOPTS)

include $(SRCDIR)/$(TOPDIR)/lib/MakeRules

# since there are no files to compile yet -- do nothing
$(TOPDIR)/lib/$(TARGET_TO_MAKE).a:
	
install_target:
	

ifneq ($(DODEPEND),no)
include $(LIBOBJ:.$(OBJSUF)=.d)
include $(TESTOBJ:%.$(OBJSUF)=%.d)
include $(TEST1OBJ:%.$(OBJSUF)=%.d)
include $(TEST2OBJ:%.$(OBJSUF)=%.d)
endif

