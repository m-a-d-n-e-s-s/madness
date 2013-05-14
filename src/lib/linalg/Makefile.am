include $(top_srcdir)/config/MakeGlobal.am

bin_PROGRAMS = test  test_solvers test_elemental
lib_LIBRARIES = libMADlinalg.a

thisincludedir = $(includedir)/linalg
thisinclude_HEADERS = tensor_lapack.h clapack.h  lapack_functions.h \
                         solvers.cc solvers.h gmres.h elem.h

libMADlinalg_a_SOURCES = lapack.cc cblas.h \
                         tensor_lapack.h clapack.h  lapack_functions.h \
                         solvers.cc solvers.h elem.cc

test_SOURCES = test.cc
test_LDADD = libMADlinalg.a $(LIBTENSOR) $(LIBMISC) $(LIBWORLD)

test_solvers_SOURCES = test_solvers.cc
test_solvers_LDADD =  libMADlinalg.a $(LIBTENSOR) $(LIBMISC) $(LIBWORLD)

test_elemental_SOURCES = test_elemental.cc
test_elemental_LDADD =  libMADlinalg.a $(LIBTENSOR) $(LIBMISC) $(LIBWORLD)


libraries:	$(lib_LIBRARIES)
	for dir in $(SUBDIRS) ; do $(MAKE) -C $$dir $@ ; done
