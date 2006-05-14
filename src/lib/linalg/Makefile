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
LIBOBJ= $(CXXSRC:%.cc=%.$(OBJSUF))

default::

TESTPROGS = test

THISLIBDEPEND = $(shell $(LISTLIBS) $(INCLUDE) $(SRCDIR)/LIBS.h)

TEST1SRC = test.cc
TEST1OBJ = $(TEST1SRC:%.cc=%.$(OBJSUF))
TEST1LIBDEPEND = 
test: $(TEST1OBJ) $(THISLIBDEPEND) $(TEST1LIBDEPEND)
	$(LTLINK) $(LD) $(LDFLAGS) -o $@ $^ $(SYSLIBS) $(LTLINKBINOPTS)

include $(SRCDIR)/$(TOPDIR)/lib/MakeRules

ifneq ($(DODEPEND),no)
include $(LIBOBJ:.$(OBJSUF)=.d)
include $(TES1OBJ:%.$(OBJSUF)=%.d)
endif

