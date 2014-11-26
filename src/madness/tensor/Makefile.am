include $(top_srcdir)/config/MakeGlobal.am
  
lib_LIBRARIES = libMADtensor.a libMADlinalg.a

TESTS = oldtest.seq test_mtxmq.seq test_Zmtxmq.seq jimkernel.seq \
        test_scott.seq test_systolic.mpi test_linalg.seq test_solvers.seq \
        test_elemental.mpi testseprep.seq test_distributed_matrix.mpi

if MADNESS_HAS_GOOGLE_TEST
TESTS += test test_gentensor
endif

TEST_EXTENSIONS = .mpi .seq

# tests run by mpirun
MPI_LOG_COMPILER = $(MPIRUN)
AM_MPI_LOG_FLAGS = $(MPIFLAGS)

# sequential tests
SEQ_LOG_COMPILER = 
AM_SEQ_LOG_FLAGS =

# default test (= sequential)
LOG_COMPILER = 
AM_LOG_FLAGS =

noinst_PROGRAMS = $(TESTS)

thisincludedir = $(includedir)/madness/tensor
thisinclude_HEADERS = aligned.h     mxm.h     tensorexcept.h  tensoriter_spec.h  type_data.h \
                        basetensor.h  tensor.h        tensor_macros.h    vector_factory.h \
                        mtxmq.h     slice.h   tensoriter.h    tensor_spec.h vmath.h gentensor.h srconf.h systolic.h \
                        tensortrain.h distributed_matrix.h \
                        tensor_lapack.h cblas.h clapack.h  lapack_functions.h \
                        solvers.cc solvers.h gmres.h elem.h

if MADNESS_HAS_GOOGLE_TEST

test_SOURCES = test.cc
test_CPPFLAGS = $(LIBGTEST_CPPFLAGS)
test_CXXFLAGS = $(LIBGTEST_CXXFLAGS)
test_LDADD = $(LIBGTEST_LIBS) $(LIBGTEST) libMADtensor.a $(LIBMISC) $(LIBWORLD) 

test_gentensor_SOURCES = test_gentensor.cc
test_gentensor_CPPFLAGS = $(LIBGTEST_CPPFLAGS)
test_gentensor_CXXFLAGS = $(LIBGTEST_CXXFLAGS)
test_gentensor_LDADD = $(LIBGTEST_LIBS) $(LIBGTEST) libMADtensor.a libMADlinalg.a  $(LIBMISC) $(LIBWORLD)

endif

oldtest_seq_SOURCES = oldtest.cc
oldtest_seq_LDADD = libMADtensor.a $(LIBMISC) $(LIBWORLD) 

test_scott_seq_SOURCES = test_scott.cc
test_scott_seq_LDADD = libMADtensor.a $(LIBMISC) $(LIBWORLD) 

test_mtxmq_seq_SOURCES = test_mtxmq.cc
test_mtxmq_seq_LDADD = libMADtensor.a $(LIBWORLD)
test_mtxmq_seq_CPPFLAGS = $(AM_CPPFLAGS) -DTIME_DGEMM

jimkernel_seq_SOURCES = jimkernel.cc
jimkernel_seq_LDADD = libMADtensor.a $(LIBMISC) $(LIBWORLD)

test_linalg_seq_SOURCES = test_linalg.cc
test_linalg_seq_LDADD = libMADlinalg.a libMADtensor.a $(LIBMISC) $(LIBWORLD)

test_solvers_seq_SOURCES = test_solvers.cc
test_solvers_seq_LDADD =  libMADlinalg.a libMADtensor.a $(LIBMISC) $(LIBWORLD)

test_elemental_mpi_SOURCES = test_elemental.cc
test_elemental_mpi_LDADD =  libMADlinalg.a libMADtensor.a $(LIBMISC) $(LIBWORLD)

test_distributed_matrix_mpi_SOURCES = test_distributed_matrix.cc
test_distributed_matrix_mpi_LDADD =  libMADtensor.a $(LIBMISC) $(LIBWORLD)

test_Zmtxmq_seq_SOURCES = test_Zmtxmq.cc
test_Zmtxmq_seq_LDADD = libMADtensor.a $(LIBWORLD)
test_Zmtxmq_seq_CPPFLAGS = $(AM_CPPFLAGS) -DTIME_DGEMM

test_systolic_mpi_SOURCES = test_systolic.cc
test_systolic_mpi_LDADD = libMADtensor.a $(LIBMISC) $(LIBWORLD)

testseprep_seq_SOURCES = testseprep.cc
testseprep_seq_LDADD = $(LIBMISC) $(LIBWORLD) libMADlinalg.a libMADtensor.a 


libMADtensor_a_SOURCES = tensor.cc tensoriter.cc basetensor.cc mtxmq.cc vmath.cc \
                        aligned.h     mxm.h     tensorexcept.h  tensoriter_spec.h  type_data.h \
                        basetensor.h  tensor.h        tensor_macros.h    vector_factory.h \
                        mtxmq.h     slice.h   tensoriter.h    tensor_spec.h vmath.h systolic.h gentensor.h srconf.h \
                        distributed_matrix.h

libMADlinalg_a_SOURCES = lapack.cc cblas.h \
                         tensor_lapack.h clapack.h  lapack_functions.h \
                         solvers.cc solvers.h elem.cc


if USE_X86_64_ASM
  libMADtensor_a_SOURCES += mtxmq_asm.S mtxm_gen.h genmtxm.py
mtxmq_asm.o:	mtxmq_asm.S mtxm_gen.h
	$(CXX) $(CXXFLAGS) $(CPPFLAGS) $(ASFLAGS) $(AM_CPPFLAGS) -I. -DX86_64 -c -o $@ $<

mtxm_gen.h:	genmtxm.py 
	python $(srcdir)/genmtxm.py > $@
endif

if USE_X86_32_ASM
  libMADtensor_a_SOURCES += mtxmq_asm.S mtxm_gen.h genmtxm.py
mtxmq_asm.o:	mtxmq_asm.S mtxm_gen.h
	$(CXX) $(CXXFLAGS) $(CPPFLAGS) $(ASFLAGS) $(AM_CPPFLAGS) -I. -DX86_32 -c -o $@ $<

mtxm_gen.h:	genmtxm.py 
	python $(srcdir)/genmtxm.py > $@
endif

tensoriter_spec.h tensor_spec.h:	tempspec.py
	python $(srcdir)/tempspec.py




#TESTPROGS = test test_mtxmq
#TEST1SRC = test.o
#TEST2SRC = test_mtxmq.o
#test_mtxmq:     test_mtxmq.o mtxmq.o mtxmq_asm.o
#        $(CXX) -o $@  $^

libraries:	$(lib_LIBRARIES)

install-libraries: install-libLIBRARIES install-thisincludeHEADERS

