include $(top_srcdir)/config/MakeGlobal.am

bin_PROGRAMS = test  test_solvers
lib_LIBRARIES = libMADlinalg.a

libMADlinalg_a_SOURCES = lapack.cc cblas.cc \
                         \
                         tensor_lapack.h cblas.h  clapack.h  lapack_functions.h \
                         solvers.cc solvers.h test_solvers.cc

test_SOURCES = test.cc
test_LDADD = libMADlinalg.a $(LIBTENSOR) $(LIBMISC) $(LIBWORLD)

test_solvers_SOURCES = test_solvers.cc
test_solvers_LDADD =  libMADlinalg.a $(LIBTENSOR) $(LIBMISC) $(LIBWORLD)
