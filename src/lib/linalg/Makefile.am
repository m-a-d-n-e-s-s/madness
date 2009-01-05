include $(top_srcdir)/config/MakeGlobal.am

bin_PROGRAMS = test kain
lib_LIBRARIES = libMADlinalg.a

libMADlinalg_a_SOURCES = lapack.cc cblas.cc \
                         \
                         tensor_lapack.h cblas.h  clapack.h  lapack_functions.h

test_SOURCES = test.cc
test_LDADD = libMADlinalg.a $(LIBTENSOR) $(LIBMISC) $(LIBWORLD)

kain_SOURCES = kain.cc
kain_LDADD = libMADlinalg.a $(LIBTENSOR) $(LIBMISC) $(LIBWORLD)

