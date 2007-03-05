TOPDIR=../../..
TARGET_TO_MAKE = libMADlinalg
CXXSRC =  lapack.cc cblas.cc
TESTPROGS = test
TEST1SRC = test.cc

include $(TOPDIR)/lib/GlobalMakefile

test: $(TEST1OBJ) $(THISLIBDEPEND)
	$(LTLINK) $(LD) $(LDFLAGS) -o $@ $^ $(SYSLIBS) $(LTLINKBINOPTS)
