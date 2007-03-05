TOPDIR=../../..
TARGET_TO_MAKE = libMADlinalg
CXXSRC =  lapack.cc cblas.cc
TESTPROGS = test
TEST1SRC = test.cc

include $(TOPDIR)/config/GlobalMakefile

test: $(TEST1OBJ) $(THISLIBDEPEND)
	$(LTLINK) $(LD) $(LDFLAGS) -o $@ $^ $(SYSLIBS) $(LTLINKBINOPTS)
