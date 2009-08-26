include $(top_srcdir)/config/MakeGlobal.am

#AM_CPPFLAGS += -DMRA_DATA_DIR="\"`pwd`\""
AM_CPPFLAGS += -DMRA_DATA_DIR=\"$(abs_srcdir)\"

bin_PROGRAMS = testsuite  testperiodic testbsh testbc testproj systolic mraplot
lib_LIBRARIES = libMADmra.a

data_DATA = autocorr coeffs gaussleg

thisincludedir = $(includedir)/mra
thisinclude_HEADERS = adquad.h  funcimpl.h  indexit.h  legendre.h  loadbal.h  operator.h  vmra.h \
                      funcdefaults.h  key.h mra.h      power.h     qmprop.h  twoscale.h lbdeux.h mraimpl.h


LDADD = libMADmra.a $(LIBLINALG) $(LIBTENSOR) $(LIBMISC) $(LIBWORLD)

libMADmra_a_SOURCES = mra1.cc mra2.cc mra3.cc mra4.cc mra5.cc mra6.cc \
                      loadbal.cc startup.cc legendre.cc twoscale.cc test_sepop.cc gfit.cc qmprop.cc \
                      \
                      adquad.h  funcimpl.h  indexit.h  legendre.h  loadbal.h  operator.h  vmra.h \
                      funcdefaults.h  key.h mra.h      power.h     qmprop.h  twoscale.h lbdeux.h mraimpl.h

testsuite_SOURCES = testsuite.cc
testperiodic_SOURCES = testperiodic.cc
testbsh_SOURCES = testbsh.cc
testbc_SOURCES = testbc.cc
testproj_SOURCES = testproj.cc

systolic_SOURCES = systolic.cc

mraplot_SOURCES = mraplot.cc



LIBRARIES:	$(lib_LIBRARIES)
	for dir in $(SUBDIRS) ; do $(MAKE) -C $$dir $@ ; done

