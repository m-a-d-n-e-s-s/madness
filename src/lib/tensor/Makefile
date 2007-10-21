TOPDIR = ../../..
TARGET_TO_MAKE = libMADtensor
CXXSRC =  tensor.cc tensoriter.cc basetensor.cc mtrand.cc 

#CSRC = mtxmq.c

#ifeq ($(HOST_CPU),X86_32) 
#  CXXSRC += mtxmq.S
#endif
#ifeq ($(HOST_CPU),X86_64) 
#  CXXSRC += mtxmq.S
#endif


TESTPROGS = test
TEST1SRC = test.o
include $(TOPDIR)/config/GlobalMakefile

test: $(TEST1OBJ) $(THISLIBDEPEND)
	$(LTLINK) $(LD) $(LDFLAGS) -o $@ $^ $(SYSLIBS) $(LTLINKBINOPTS)
