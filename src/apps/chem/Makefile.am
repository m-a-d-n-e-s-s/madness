include $(top_srcdir)/config/MakeGlobal.am
EXTRA_DIST = CMakeLists.txt

# specify the git commit version in the info functions
GITREV="MADNESS_GITREVISION=\"`git --git-dir=$(top_srcdir)/.git rev-parse HEAD`\" "

# default location for basis sets etc
AM_CPPFLAGS += -DMRA_CHEMDATA_DIR=\"$(abs_srcdir)\" -D$(GITREV)

noinst_PROGRAMS = testxc plotxc test_SCFOperators test_dft

lib_LTLIBRARIES = libMADchem.la

chemdatadir=$(pkgdatadir)/$(PACKAGE_VERSION)/data
dist_chemdata_DATA = sto-3g sto-6g 6-31g coredata/mcp coredata/mcp2 coredata/mcp_guess

TESTS = test_SCFOperators test_dft # testxc 

thisincludedir = $(includedir)/chem

thisinclude_HEADERS = correlationfactor.h molecule.h molecularbasis.h \
                      corepotential.h atomutil.h SCF.h xcfunctional.h \
                      mp2.h nemo.h potentialmanager.h gth_pseudopotential.h \
                      molecular_optimizer.h projector.h \
                      SCFOperators.h CCStructures.h pointgroupoperator.h pointgroupsymmetry.h \
                      electronic_correlation_factor.h cheminfo.h vibanal.h molopt.h TDHF.h \
                      CC2.h CCPotentials.h AC.h GuessFactory.h  \
                      polynomial.h gaussian.h ESInterface.h NWChem.h basis.h

testxc_SOURCES = testxc.cc xcfunctional.h
testxc_LDADD = libMADchem.la $(MRALIBS)

test_dft_SOURCES = test_dft.cc xcfunctional.h
test_dft_LDADD = libMADchem.la $(MRALIBS)

test_SCFOperators_SOURCES = test_SCFOperators.cc 
test_SCFOperators_LDADD = libMADchem.la $(MRALIBS)

plotxc_SOURCES = plotxc.cc xcfunctional.h
plotxc_LDADD = libMADchem.la $(MRALIBS)

libMADchem_la_SOURCES = correlationfactor.cc molecule.cc molecularbasis.cc vibanal.cc \
                       corepotential.cc atomutil.cc lda.cc cheminfo.cc \
                       distpm.cc SCF.cc gth_pseudopotential.cc nemo.cc mp2.cc pcm.cc pointgroupsymmetry.cc \
                       SCFOperators.cc xcfunctional_ldaonly.cc TDHF.cc CCStructures.cc CC2.cc CCPotentials.cc \
                       AC.cc GuessFactory.cc\
                       polynomial.cc gaussian.cc NWChem.cc \
                       $(thisinclude_HEADERS)

libMADchem_la_LDFLAGS = -version-info 0:0:0
                       
if MADNESS_HAS_LIBXC
   libMADchem_la_SOURCES += xcfunctional_libxc.cc
endif



libraries:	$(lib_LTLIBRARIES)
	for dir in $(SUBDIRS) ; do $(MAKE) -C $$dir $@ ; done

install-libraries: install-libLTLIBRARIES install-thisincludeHEADERS
	for dir in $(SUBDIRS) ; do $(MAKE) -C $$dir $@ ; done
