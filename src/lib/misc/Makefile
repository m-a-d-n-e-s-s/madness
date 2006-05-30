TOPDIR=../../..
TARGET_TO_MAKE = libMADmisc
CXXSRC =  walltime.cc checksum_file.cc redirectio.cc am_stuff.cc debug.cc 
#initcxx.cc
TESTPROGS = test
TEST1SRC = test.cc 

include $(TOPDIR)/lib/GlobalMakefile


test: $(TEST1OBJ) $(THISLIBDEPEND)
	$(LTLINK) $(LD) $(LDFLAGS) -o $@ $^ $(SYSLIBS) $(LTLINKBINOPTS)


