TOPDIR = ../../..
TARGET_TO_MAKE = libMADtensor
CXXSRC =  tensor.cc tensoriter.cc basetensor.cc mtrand.cc 
TESTPROGS = test
TEST1SRC = test.o
include $(TOPDIR)/lib/GlobalMakefile

test: $(TEST1OBJ) $(THISLIBDEPEND)
	$(LTLINK) $(LD) $(LDFLAGS) -o $@ $^ $(SYSLIBS) $(LTLINKBINOPTS)
