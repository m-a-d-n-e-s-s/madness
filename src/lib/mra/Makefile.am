include $(top_srcdir)/config/MakeGlobal.am

#AM_CPPFLAGS += -DMRA_DATA_DIR="\"`pwd`\""
AM_CPPFLAGS += -DMRA_DATA_DIR=\"$(abs_srcdir)\"

bin_PROGRAMS = testsuite  testperiodic testbsh testbc testproj mraplot testpdiff testdiff1Db testdiff1D testdiff2D testdiff3D testgconv testopdir testqm testper
lib_LIBRARIES = libMADmra.a

#TESTS = testsuite  testperiodic testbc testbsh testproj testpdiff \
#        testdiff1Db testdiff1D testdiff2D testdiff3D testgconv testopdir \
#        testqm testunaryop testper

TESTS = testbsh testproj testpdiff testper \
        testdiff1Db testdiff1D testdiff2D testdiff3D testgconv testopdir \
        testsuite 


data_DATA = autocorr coeffs gaussleg

thisincludedir = $(includedir)/mra
thisinclude_HEADERS = adquad.h  funcimpl.h  indexit.h  legendre.h  operator.h  vmra.h \
                      funcdefaults.h  key.h mra.h      power.h     qmprop.h  twoscale.h lbdeux.h mraimpl.h \
                      funcplot.h


LDADD = libMADmra.a $(LIBLINALG) $(LIBLINALG) $(LIBTENSOR) $(LIBMISC) $(LIBMUPARSER) $(LIBWORLD)

libMADmra_a_SOURCES = mra1.cc mra2.cc mra3.cc mra4.cc mra5.cc mra6.cc \
                      startup.cc legendre.cc twoscale.cc gfit.cc qmprop.cc \
                      \
                      adquad.h  funcimpl.h  indexit.h  legendre.h  operator.h  vmra.h \
                      funcdefaults.h  key.h mra.h      power.h     qmprop.h  twoscale.h lbdeux.h mraimpl.h \
                      funcplot.h

testsuite_SOURCES = testsuite.cc test_sepop.cc
testperiodic_SOURCES = testperiodic.cc
testper_SOURCES = testper.cc test_sepop.cc
testbsh_SOURCES = testbsh.cc
testbc_SOURCES = testbc.cc
testproj_SOURCES = testproj.cc
testgconv_SOURCES = testgconv.cc

testopdir_SOURCES = testopdir.cc

mraplot_SOURCES = mraplot.cc

testpdiff_SOURCES = testpdiff.cc

testdiff1D_SOURCES = testdiff1D.cc

testdiff1Db_SOURCES = testdiff1Db.cc

testdiff2D_SOURCES = testdiff2D.cc

testdiff3D_SOURCES = testdiff3D.cc

testqm_SOURCES = testqm.cc

#testop2_SOURCES = testop2.cc


libraries:	$(lib_LIBRARIES)
	for dir in $(SUBDIRS) ; do $(MAKE) -C $$dir $@ ; done

install-libraries: install-libLIBRARIES install-thisincludeHEADERS
	for dir in $(SUBDIRS) ; do $(MAKE) -C $$dir $@ ; done

