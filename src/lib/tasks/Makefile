TOPDIR=../../..
TARGET_TO_MAKE = libMADtasks

CXXSRC = tasks.cc
TESTPROGS = test 
TEST1SRC = test.cc
TEST2SRC = testserver.cc

include $(TOPDIR)/lib/GlobalMakefile

test: $(TEST1OBJ) $(THISLIBDEPEND)
	$(LTLINK) $(LD) $(LDFLAGS) -o $@ $^ $(SYSLIBS) $(LTLINKBINOPTS)

testserver: $(TEST2OBJ) $(THISLIBDEPEND)
	$(LTLINK) $(LD) $(LDFLAGS) -o $@ $^ $(SYSLIBS) $(LTLINKBINOPTS)
	