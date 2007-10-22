TOPDIR = ../../..
TARGET_TO_MAKE = libMADtensor
CXXSRC =  tensor.cc tensoriter.cc basetensor.cc mtrand.cc 

TESTPROGS = test
TEST1SRC = test.o
include $(TOPDIR)/config/GlobalMakefile


USE_MTXMQ = yes
TUNE_MTXMQ = -DCORE_DUO_TUNE

ifneq ($(USE_MTXMQ),"")
CXXSRC += mtxmq.cc 
mtxmq.o:	mtxmq.S mtxm_gen.h
	$(CC) $(CXXFLAGS) -D$(HOST_CPU) $(TUNE_MTXMQ) -c -o $@ 
endif


test: $(TEST1OBJ) $(THISLIBDEPEND)
	$(LTLINK) $(LD) $(LDFLAGS) -o $@ $^ $(SYSLIBS) $(LTLINKBINOPTS)


