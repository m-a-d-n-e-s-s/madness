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
INC =     archive.h textfsar.h binfsar.h vecar.h mpiar.h

TESTSRC = test.cc
TESTOBJ = $(TESTSRC:%.cc=%.$(OBJSUF))
TESTFILES = 

TESTPROGS = test

DISTFILES = $(CXXSRC) $(INC) Makefile $(TESTSRC) $(TESTFILES)

LIBOBJ= $(CXXSRC:%.cc=%.$(OBJSUF))

default:: $(DEPENDINCLUDE)

test: $(TESTOBJ) \
      $(shell $(LISTLIBS) $(INCLUDE) $(SRCDIR)/LIBS.h) \
      $(shell $(LISTLIBS) $(INCLUDE) $(SRCDIR)/../tensor/LIBS.h)
	$(LTLINK) $(LD) $(LDFLAGS) -o $@ $^ $(SYSLIBS) $(LTLINKBINOPTS)


testtens.o:	$(INC)

testtens: testtens.o \
      $(shell $(LISTLIBS) $(INCLUDE) $(SRCDIR)/LIBS.h) \
      $(shell $(LISTLIBS) $(INCLUDE) $(SRCDIR)/../tensor/LIBS.h)
	$(LTLINK) $(LD) $(LDFLAGS) -o $@ $^ $(SYSLIBS) $(LTLINKBINOPTS)


testmpi.o:	$(INC)

testmpi: testmpi.o \
      $(shell $(LISTLIBS) $(INCLUDE) $(SRCDIR)/LIBS.h) \
      $(shell $(LISTLIBS) $(INCLUDE) $(SRCDIR)/../tensor/LIBS.h)
	$(LTLINK) $(LD) $(LDFLAGS) -o $@ $^ $(SYSLIBS) $(LTLINKBINOPTS)

testfn: testfn.o \
      $(shell $(LISTLIBS) $(INCLUDE) $(SRCDIR)/LIBS.h) \
      $(shell $(LISTLIBS) $(INCLUDE) $(SRCDIR)/../tensor/LIBS.h)
	$(LTLINK) $(LD) $(LDFLAGS) -o $@ $^ $(SYSLIBS) $(LTLINKBINOPTS)

testfn.o:	$(INC)

testvec: testvec.o \
      $(shell $(LISTLIBS) $(INCLUDE) $(SRCDIR)/LIBS.h) \
      $(shell $(LISTLIBS) $(INCLUDE) $(SRCDIR)/../tensor/LIBS.h)
	$(LTLINK) $(LD) $(LDFLAGS) -o $@ $^ $(SYSLIBS) $(LTLINKBINOPTS)

testvec.o:	$(INC)

testoct: testoct.o \
      $(shell $(LISTLIBS) $(INCLUDE) $(SRCDIR)/LIBS.h) \
      $(shell $(LISTLIBS) $(INCLUDE) $(SRCDIR)/../tensor/LIBS.h)
	$(LTLINK) $(LD) $(LDFLAGS) -o $@ $^ $(SYSLIBS) $(LTLINKBINOPTS)

testoct.o:	$(INC)

testtree: testtree.o \
      $(shell $(LISTLIBS) $(INCLUDE) $(SRCDIR)/LIBS.h) \
      $(shell $(LISTLIBS) $(INCLUDE) $(SRCDIR)/../tensor/LIBS.h)
	$(LTLINK) $(LD) $(LDFLAGS) -o $@ $^ $(SYSLIBS) $(LTLINKBINOPTS)

testtree.o:	$(INC)

include $(SRCDIR)/$(TOPDIR)/lib/MakeRules

clean::
	rm -f $(TESTPROGS) $(TESTOBJ)

$(LIBOBJ:.$(OBJSUF)=.d): $(DEPENDINCLUDE)
ifneq ($(DODEPEND),no)
include $(LIBOBJ:.$(OBJSUF)=.d) $(TESTOBJ:%.$(OBJSUF)=%.d)
endif

