TOPDIR=../../..
ifndef QUOTE
  QUOTE=sed 's/ /\\ /g'
endif
ifndef SRCDIR
  SRCDIR=$(shell pwd | $(QUOTE))
endif

include $(SRCDIR)/$(TOPDIR)/lib/GlobalMakefile
TARGET_TO_MAKE = libMADserialize
BIN_OR_LIB = LIB

CXXSRC =  archive_type_names.cc
LIBOBJ= $(CXXSRC:%.cc=%.$(OBJSUF))

default::

TESTPROGS = test testmpi

THISLIBDEPEND = $(shell $(LISTLIBS) $(INCLUDE) $(SRCDIR)/LIBS.h)

TEST1SRC = test.cc
TEST1OBJ = $(TEST1SRC:%.cc=%.$(OBJSUF))
TEST1LIBDEPEND = 
TEST1PROGS = test
test: $(TEST1OBJ) $(THISLIBDEPEND) $(TEST1LIBDEPEND)
	$(LTLINK) $(LD) $(LDFLAGS) -o $@ $^ $(SYSLIBS) $(LTLINKBINOPTS)

TEST2SRC = testmpi.cc
TEST2OBJ = $(TEST2SRC:%.cc=%.$(OBJSUF))
TEST2LIBDEPEND = 
TEST2PROGS = test
test: $(TEST2OBJ) $(THISLIBDEPEND) $(TEST2LIBDEPEND)
	$(LTLINK) $(LD) $(LDFLAGS) -o $@ $^ $(SYSLIBS) $(LTLINKBINOPTS)

include $(SRCDIR)/$(TOPDIR)/lib/MakeRules

ifneq ($(DODEPEND),no)
include $(LIBOBJ:.$(OBJSUF)=.d)
include $(TEST1OBJ:%.$(OBJSUF)=%.d)
include $(TEST2OBJ:%.$(OBJSUF)=%.d)
endif

