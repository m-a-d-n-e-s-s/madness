TOPDIR=../../..
TARGET_TO_MAKE = libMADlinalg
CXXSRC =  lapack.cc cblas.cc
TESTPROGS = test
TEST1SRC = test.cc
TEST2SRC = kain.cc
TEST3SRC = fred.cc

include $(TOPDIR)/config/GlobalMakefile

test: $(TEST1OBJ) $(THISLIBDEPEND)
	$(LTLINK) $(LD) $(LDFLAGS) -o $@ $^ $(SYSLIBS) $(LTLINKBINOPTS)

kain: $(TEST2OBJ) $(THISLIBDEPEND)
	$(LTLINK) $(LD) $(LDFLAGS) -o $@ $^ $(SYSLIBS) $(LTLINKBINOPTS)

fred: $(TEST3OBJ) $(THISLIBDEPEND)
	$(LTLINK) $(LD) $(LDFLAGS) -o $@ $^ $(SYSLIBS) $(LTLINKBINOPTS)
