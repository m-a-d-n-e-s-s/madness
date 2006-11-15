TOPDIR=../../..
TARGET_TO_MAKE = libMADmisc
CXXSRC =  walltime.cc checksum_file.cc redirectio.cc am_stuff.cc debug.cc madexcept.cc
TESTPROGS = test
TEST1SRC = test.cc 
TEST2SRC = am.cc
TEST3SRC = diskdir.cc

include $(TOPDIR)/lib/GlobalMakefile


test: $(TEST1OBJ) $(THISLIBDEPEND)
	$(LTLINK) $(LD) $(LDFLAGS) -o $@ $^ $(SYSLIBS) $(LTLINKBINOPTS)

am: $(TEST2OBJ) $(THISLIBDEPEND)
	$(LTLINK) $(LD) $(LDFLAGS) -o $@ $^ $(SYSLIBS) $(LTLINKBINOPTS)

diskdir: $(TEST3OBJ) 
	$(LTLINK) $(LD) $(LDFLAGS) -o $@ $^ ../mra/legendre.o  ../mra/sepop.o  ../mra/gfit.o  ../mra/twoscale.o ../../../lib/libMADmra.a ../../../lib/libMADtensor.a ../../../lib/libMADmisc.a ../../../lib/libMADserialize.a ../../../lib/libMADocttree.a ../../../lib/libMADtasks.a /usr/lib/liblapack.a /usr/lib/libblas.a  -lpthread 

