include $(top_srcdir)/config/MakeGlobal.am

bin_PROGRAMS = test test_mtxmq test_Zmtxmq
lib_LIBRARIES = libMADtensor.a

test_SOURCES = test.cc
test_LDADD = libMADtensor.a $(LIBMISC) $(LIBWORLD) 

test_mtxmq_SOURCES = test_mtxmq.cc
test_mtxmq_LDADD = libMADtensor.a $(LIBWORLD)
test_mtxmq_CPPFLAGS = $(AM_CPPFLAGS) -DTIME_DGEMM

test_Zmtxmq_SOURCES = test_Zmtxmq.cc
test_Zmtxmq_LDADD = libMADtensor.a $(LIBWORLD)
test_Zmtxmq_CPPFLAGS = $(AM_CPPFLAGS) -DTIME_DGEMM

libMADtensor_a_SOURCES = tensor.cc tensoriter.cc basetensor.cc mtxmq.cc aligned.cc alignedc.cc \
                         \
                        aligned.h     mxm.h     tensorexcept.h  tensoriter_spec.h  type_data.h \
                        basetensor.h  tensor.h        tensor_macros.h    vector_factory.h \
                        mtxmq.h     slice.h   tensoriter.h    tensor_spec.h

if USE_X86_64_ASM
  libMADtensor_a_SOURCES += mtxmq_asm.S mtxm_gen.h genmtxm.py
mtxmq_asm.o:	mtxmq_asm.S mtxm_gen.h
	$(CC) $(CXXFLAGS) -DX86_64 -c -o $@ $<

mtxm_gen.h:	genmtxm.py 
	python genmtxm.py > $@
endif

if USE_X86_32_ASM
  libMADtensor_a_SOURCES += mtxmq_asm.S mtxm_gen.h genmtxm.py
mtxmq_asm.o:	mtxmq_asm.S mtxm_gen.h
	$(CC) $(CXXFLAGS) -DX86_32 -c -o $@ $<

mtxm_gen.h:	genmtxm.py 
	python genmtxm.py > $@
endif

tensoriter_spec.h tensor_spec.h:	tempspec.py
	python tempspec.py




#TESTPROGS = test test_mtxmq
#TEST1SRC = test.o
#TEST2SRC = test_mtxmq.o
#test_mtxmq:     test_mtxmq.o mtxmq.o mtxmq_asm.o
#        $(CXX) -o $@  $^

