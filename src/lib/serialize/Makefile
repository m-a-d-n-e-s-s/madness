TOPDIR = ../../..
TARGET_TO_MAKE = libMADserialize
CXXSRC =  archive_type_names.cc
TESTPROGS = test testmpi
TEST1SRC = test.cc
TEST2SRC = testmpi.cc

include $(TOPDIR)/lib/GlobalMakefile

test: $(TEST1OBJ) $(THISLIBDEPEND) 
	$(LTLINK) $(LD) $(LDFLAGS) -o $@ $^ $(SYSLIBS) $(LTLINKBINOPTS)

testmpi: $(TEST2OBJ) $(THISLIBDEPEND) 
	$(LTLINK) $(LD) $(LDFLAGS) -o $@ $^ $(SYSLIBS) $(LTLINKBINOPTS)
