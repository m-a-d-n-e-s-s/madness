TOPDIR = ../../..
TARGET_TO_MAKE = libMADtensor
CXXSRC =  tensor.cc tensoriter.cc basetensor.cc mtrand.cc mtxmq.cc

TESTPROGS = test test_mtxmq
TEST1SRC = test.o
TEST2SRC = test_mtxmq.o

include $(TOPDIR)/config/MakeVars

ifeq ($(HOST_CPU),X86_32)
USE_MTXMQ = YES
TUNE_MTXMQ = 
endif
ifeq ($(HOST_CPU),X86_64)
USE_MTXMQ = YES
TUNE_MTXMQ = 
endif


include $(TOPDIR)/config/GlobalMakefile


test: $(TEST1OBJ) $(THISLIBDEPEND)
	$(LTLINK) $(LD) $(LDFLAGS) -o $@ $^ $(SYSLIBS) $(LTLINKBINOPTS)

tensoriter_spec.h tensor_spec.h:	tempspec.py
	python tempspec.py


ifeq ($(USE_MTXMQ),YES)
LIBOBJ += mtxmq_asm.o

$(TOPDIR)/lib/$(TARGET_TO_MAKE).$(LIBSUF):	mtxmq_asm.o

mtxmq_asm.o:	mtxmq_asm.S mtxm_gen.h
	$(CC) $(CXXFLAGS) -D$(HOST_CPU) $(TUNE_MTXMQ) -c -o $@ mtxmq_asm.S

mtxm_gen.h:	genmtxm.py 
	python genmtxm.py > $@

test_mtxmq:	test_mtxmq.o mtxmq.o mtxmq_asm.o
	$(CXX) -o $@  $^
endif

