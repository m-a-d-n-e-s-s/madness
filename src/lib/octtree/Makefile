TOPDIR=../../..
TARGET_TO_MAKE = libMADocttree

CXXSRC =  
TESTPROGS = test test1 testpart
TEST1SRC = test.cc
TEST2SRC = test1.cc
TEST3SRC = testpart.cc

include $(TOPDIR)/lib/GlobalMakefile

test: $(TEST1OBJ) $(THISLIBDEPEND) 
	$(LTLINK) $(LD) $(LDFLAGS) -o $@ $^ $(SYSLIBS) $(LTLINKBINOPTS)

test1: $(TEST2OBJ) $(THISLIBDEPEND) 
	$(LTLINK) $(LD) $(LDFLAGS) -o $@ $^ $(SYSLIBS) $(LTLINKBINOPTS)

testpart: $(TESTOBJ) $(THISLIBDEPEND) 
	$(LTLINK) $(LD) $(LDFLAGS) -o $@ $^ $(SYSLIBS) $(LTLINKBINOPTS)

# since there are no files to compile yet -- do nothing
$(TOPDIR)/lib/$(TARGET_TO_MAKE).a:
	
install_target:
