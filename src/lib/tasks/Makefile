TOPDIR=../../..
TARGET_TO_MAKE = libMADtasks

CXXSRC =  
TESTPROGS = test test2
TEST1SRC = test.cc
TEST2SRC = test2.cc

include $(TOPDIR)/lib/GlobalMakefile

test: $(TEST1OBJ) $(THISLIBDEPEND)
	$(LTLINK) $(LD) $(LDFLAGS) -o $@ $^ $(SYSLIBS) $(LTLINKBINOPTS)

test2: $(TEST2OBJ) $(THISLIBDEPEND)
	$(LTLINK) $(LD) $(LDFLAGS) -o $@ $^ $(SYSLIBS) $(LTLINKBINOPTS)
	