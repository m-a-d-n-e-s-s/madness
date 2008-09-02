include $(top_srcdir)/config/MakeGlobal.am

bin_PROGRAMS = test
lib_LIBRARIES = libMADlinalg.a

libMADlinalg_a_SOURCES = lapack.cc cblas.cc \
                         \
                         tensor_lapack.h cblas.h  clapack.h  lapack_functions.h

test_SOURCES = test.cc
test_LDADD = libMADlinalg.a $(LIBTENSOR) $(LIBWORLD)
