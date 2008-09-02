include $(top_srcdir)/config/MakeGlobal.am

AM_CPPFLAGS += -DMRA_DATA_DIR="`pwd`"

bin_PROGRAMS = testsuite  testperiodic testbsh
lib_LIBRARIES = libMADmra.a

DATAFILES = autocorr coeffs gaussleg

LDADD = libMADmra.a $(LIBLINALG) $(LIBTENSOR) $(LIBMISC) $(LIBWORLD)

libMADmra_a_SOURCES = mra.cc loadbal.cc startup.cc legendre.cc twoscale.cc test_sepop.cc gfit.cc qmprop.cc \
                      \
                      adquad.h  funcimpl.h  indexit.h  legendre.h  loadbal.h  operator.h  vmra.h \
                      funcdefaults.h  key.h mra.h      power.h     qmprop.h  twoscale.h

testsuite_SOURCES = testsuite.cc
testperiodic_SOURCES = testperiodic.cc
testbsh_SOURCES = testbsh.cc


