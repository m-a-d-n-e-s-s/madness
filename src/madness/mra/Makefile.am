include $(top_srcdir)/config/MakeGlobal.am
EXTRA_DIST = CMakeLists.txt

#AM_CPPFLAGS += -DMRA_DATA_DIR="\"`pwd`\""
AM_CPPFLAGS += -DMRA_DATA_DIR=\"$(abs_srcdir)\"

#TESTS = testsuite  testperiodic testbc testbsh testproj testpdiff \
#        testdiff1Db testdiff1D testdiff2D testdiff3D testgconv testopdir \
#        testqm testunaryop testper

TESTS = testbsh.mpi testproj.mpi testpdiff.mpi testper.mpi \
        testdiff1Db.mpi \
		testgconv.mpi testopdir.mpi testsuite.mpi testinnerext.mpi \
		testgaxpyext.mpi testvmra.mpi


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



bin_PROGRAMS = mraplot
noinst_PROGRAMS =  testperiodic.mpi testbc.mpi testproj.mpi testqm test6 \
                   testdiff1D.mpi testdiff2D.mpi testdiff3D.mpi $(TESTS)
lib_LTLIBRARIES = libMADmra.la

mradatadir=${pkgdatadir}/$(PACKAGE_VERSION)/data
dist_mradata_DATA = autocorr coeffs gaussleg

thisincludedir = $(includedir)/madness/mra
thisinclude_HEADERS = adquad.h  funcimpl.h  indexit.h  legendre.h  operator.h  vmra.h \
                      funcdefaults.h  key.h  mra.h  power.h  qmprop.h  twoscale.h \
                      lbdeux.h  mraimpl.h  funcplot.h  function_common_data.h \
                      function_factory.h function_interface.h gfit.h convolution1d.h \
                      simplecache.h derivative.h displacements.h functypedefs.h \
                      sdf_shape_3D.h sdf_domainmask.h vmra1.h nonlinsol.h 


LDADD = libMADmra.la $(LIBLINALG) $(LIBTENSOR) $(LIBMISC) $(LIBMUPARSER) $(LIBWORLD)

libMADmra_la_SOURCES = mra1.cc mra2.cc mra3.cc mra4.cc mra5.cc mra6.cc \
                      startup.cc legendre.cc twoscale.cc qmprop.cc \
                      $(thisinclude_HEADERS)
libMADmra_la_LDFLAGS = -version-info 0:0:0


testsuite_mpi_SOURCES = testsuite.cc test_sepop.cc
testperiodic_mpi_SOURCES = testperiodic.cc
testper_mpi_SOURCES = testper.cc test_sepop.cc
testbsh_mpi_SOURCES = testbsh.cc
testvmra_mpi_SOURCES = testvmra.cc
test6_SOURCES = test6.cc

testbc_mpi_SOURCES = testbc.cc
testproj_mpi_SOURCES = testproj.cc
testgconv_mpi_SOURCES = testgconv.cc

testopdir_mpi_SOURCES = testopdir.cc

mraplot_SOURCES = mraplot.cc

testpdiff_mpi_SOURCES = testpdiff.cc

testdiff1D_mpi_SOURCES = testdiff1D.cc

testdiff1Db_mpi_SOURCES = testdiff1Db.cc

testdiff2D_mpi_SOURCES = testdiff2D.cc

testdiff3D_mpi_SOURCES = testdiff3D.cc

testqm_SOURCES = testqm.cc

testinnerext_mpi_SOURCES = testinnerext.cc

testgaxpyext_mpi_SOURCES = testgaxpyext.cc

#testop2_SOURCES = testop2.cc


libraries:	$(lib_LTLIBRARIES)
	for dir in $(SUBDIRS) ; do $(MAKE) -C $$dir $@ ; done

install-libraries: install-libLTLIBRARIES install-thisincludeHEADERS
	for dir in $(SUBDIRS) ; do $(MAKE) -C $$dir $@ ; done

MOSTLYCLEANFILES = fplot.dat opf.dat opfplot.dat testline1 testline2 testline3
