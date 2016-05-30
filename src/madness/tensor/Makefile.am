include $(top_srcdir)/config/MakeGlobal.am
  
lib_LTLIBRARIES = libMADtensor.la libMADlinalg.la

TESTS = oldtest.seq test_mtxmq.seq test_Zmtxmq.seq jimkernel.seq \
        test_linalg.seq test_solvers.seq \
        test_elemental.mpi testseprep.seq test_distributed_matrix.mpi

if MADNESS_HAS_GOOGLE_TEST
TESTS += test_tensor test_gentensor
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

noinst_PROGRAMS = $(TESTS) test_systolic.mpi

thisincludedir = $(includedir)/madness/tensor
thisinclude_HEADERS = aligned.h     mxm.h     tensorexcept.h  tensoriter_spec.h  type_data.h \
                        basetensor.h  tensor.h        tensor_macros.h    vector_factory.h \
                        slice.h   tensoriter.h    tensor_spec.h vmath.h gentensor.h srconf.h systolic.h \
                        tensortrain.h distributed_matrix.h \
                        tensor_lapack.h cblas.h clapack.h \
                        solvers.cc solvers.h gmres.h elem.h
EXTRA_DIST = CMakeLists.txt genmtxm.py tempspec.py

if MADNESS_HAS_GOOGLE_TEST

test_tensor_SOURCES = test_tensor.cc
test_tensor_CPPFLAGS = $(LIBGTEST_CPPFLAGS)
test_tensor_CXXFLAGS = $(LIBGTEST_CXXFLAGS)
test_tensor_LDADD = $(LIBGTEST_LIBS) $(LIBGTEST) libMADtensor.la $(LIBMISC) $(LIBWORLD) 

test_gentensor_SOURCES = test_gentensor.cc
test_gentensor_CPPFLAGS = $(LIBGTEST_CPPFLAGS)
test_gentensor_CXXFLAGS = $(LIBGTEST_CXXFLAGS)
test_gentensor_LDADD = $(LIBGTEST_LIBS) $(LIBGTEST) libMADtensor.la libMADlinalg.la  $(LIBMISC) $(LIBWORLD)

endif

oldtest_seq_SOURCES = oldtest.cc
oldtest_seq_LDADD = libMADtensor.la $(LIBMISC) $(LIBWORLD) 

test_mtxmq_seq_SOURCES = test_mtxmq.cc
test_mtxmq_seq_LDADD = libMADtensor.la $(LIBWORLD)
test_mtxmq_seq_CPPFLAGS = $(AM_CPPFLAGS) -DTIME_DGEMM

jimkernel_seq_SOURCES = jimkernel.cc
jimkernel_seq_LDADD = libMADtensor.la $(LIBMISC) $(LIBWORLD)

test_linalg_seq_SOURCES = test_linalg.cc
test_linalg_seq_LDADD = libMADlinalg.la libMADtensor.la $(LIBMISC) $(LIBWORLD)

test_solvers_seq_SOURCES = test_solvers.cc
test_solvers_seq_LDADD =  libMADlinalg.la libMADtensor.la $(LIBMISC) $(LIBWORLD)

test_elemental_mpi_SOURCES = test_elemental.cc
test_elemental_mpi_LDADD =  libMADlinalg.la libMADtensor.la $(LIBMISC) $(LIBWORLD)

test_distributed_matrix_mpi_SOURCES = test_distributed_matrix.cc
test_distributed_matrix_mpi_LDADD =  libMADtensor.la $(LIBMISC) $(LIBWORLD)

test_Zmtxmq_seq_SOURCES = test_Zmtxmq.cc
test_Zmtxmq_seq_LDADD = libMADtensor.la $(LIBWORLD)
test_Zmtxmq_seq_CPPFLAGS = $(AM_CPPFLAGS) -DTIME_DGEMM

test_systolic_mpi_SOURCES = test_systolic.cc
test_systolic_mpi_LDADD = libMADtensor.la $(LIBMISC) $(LIBWORLD)

testseprep_seq_SOURCES = testseprep.cc
testseprep_seq_LDADD = $(LIBMISC) $(LIBWORLD) libMADlinalg.la libMADtensor.la 

libMADtensor_la_SOURCES = tensor.cc tensoriter.cc basetensor.cc vmath.cc \
                        aligned.h     mxm.h     tensorexcept.h  tensoriter_spec.h  type_data.h \
                        basetensor.h  tensor.h        tensor_macros.h    vector_factory.h \
                        mtxmq.h     slice.h   tensoriter.h    tensor_spec.h vmath.h systolic.h gentensor.h srconf.h \
                        distributed_matrix.h
libMADtensor_la_LDFLAGS = -version-info 0:0:0

libMADlinalg_la_SOURCES = lapack.cc cblas.h \
                         tensor_lapack.h clapack.h  lapack_functions.h \
                         solvers.cc solvers.h elem.cc
libMADlinalg_la_LDFLAGS = -version-info 0:0:0


tensoriter_spec.h tensor_spec.h:	tempspec.py
	python $(srcdir)/tempspec.py


libraries:	$(lib_LTLIBRARIES)

install-libraries: install-libLTLIBRARIES install-thisincludeHEADERS

