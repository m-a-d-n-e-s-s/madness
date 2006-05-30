TOPDIR=../../..
TARGET_TO_MAKE = libMADocttree

CXXSRC =  
TESTPROGS = test test1 testpart
TEST1SRC = test.cc
TEST2SRC = test1.cc
TEST3SRC = testpart.cc

TESTLSRC = testlong.cc
TESTLOBJ = $(TESTLSRC:%.cc=%.$(OBJSUF))
TESTLFILES = 

TESTLPROGS = testlong

TEST1SRC = test1.cc
TEST1OBJ = $(TEST1SRC:%.cc=%.$(OBJSUF))
TEST1FILES = 
include $(TOPDIR)/lib/GlobalMakefile

TEST1PROGS = test1

TEST123SRC = testsr123.cc
TEST123OBJ = $(TEST123SRC:%.cc=%.$(OBJSUF))
TEST123FILES = 

TEST123PROGS = test123

TESTPARASRC = testpara.cc
TESTPARAOBJ = $(TESTPARASRC:%.cc=%.$(OBJSUF))
TESTPARAFILES = 

TESTPARAPROGS = testpara

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

testlong: $(TESTLOBJ) $(THISLIBDEPEND) $(TESTLIBDEPEND) ../../../lib/libMADserialize.a
	$(LTLINK) $(LD) $(LDFLAGS) -o $@ $^ $(SYSLIBS) $(LTLINKBINOPTS)

test1: $(TEST1OBJ) $(THISLIBDEPEND) $(TESTLIBDEPEND) ../../../lib/libMADserialize.a

test: $(TEST1OBJ) $(THISLIBDEPEND) 

testpara: $(TESTPARAOBJ) $(THISLIBDEPEND) $(TESTLIBDEPEND) ../../../lib/libMADserialize.a
	$(LTLINK) $(LD) $(LDFLAGS) -o $@ $^ $(SYSLIBS) $(LTLINKBINOPTS)

testpart: $(TESTPARTOBJ) $(THISLIBDEPEND) $(TESTLIBDEPEND) ../../../lib/libMADserialize.a
#testpart: $(TESTOBJ) $(THISLIBDEPEND) 
#	$(LTLINK) $(LD) $(LDFLAGS) -o $@ $^ $(SYSLIBS) $(LTLINKBINOPTS)

# since there are no files to compile yet -- do nothing
$(TOPDIR)/lib/$(TARGET_TO_MAKE).a:
	
install_target:
