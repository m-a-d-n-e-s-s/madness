include $(top_srcdir)/config/MakeGlobal.am

# default location for basis sets etc
AM_CPPFLAGS += -DMRA_CHEMDATA_DIR=\"$(abs_srcdir)\"

noinst_PROGRAMS = testxc plotxc test_SCFOperators test_dft

lib_LTLIBRARIES = libMADchem.la

data_DATA = sto-3g sto-6g 6-31g coredata/mcp coredata/mcp2 coredata/mcp_guess
dist_data_DATA = $(data_DATA)

TESTS = test_SCFOperators test_dft # testxc 

thisincludedir = $(includedir)/chem

thisinclude_HEADERS = correlationfactor.h molecule.h molecularbasis.h \
                      corepotential.h atomutil.h SCF.h xcfunctional.h \
                      mp2.h nemo.h potentialmanager.h gth_pseudopotential.h \
                      molecular_optimizer.h projector.h TDA.h TDA_XC.h \
                      TDA_guess.h TDA_exops.h SCFOperators.h CCOperators.h \
                      electronic_correlation_factor.h

testxc_SOURCES = testxc.cc xcfunctional.h
testxc_LDADD = libMADchem.la $(MRALIBS)

test_dft_SOURCES = test_dft.cc xcfunctional.h
test_dft_LDADD = libMADchem.la $(MRALIBS)

test_SCFOperators_SOURCES = test_SCFOperators.cc 
test_SCFOperators_LDADD = libMADchem.la $(MRALIBS)

plotxc_SOURCES = plotxc.cc xcfunctional.h
plotxc_LDADD = libMADchem.la $(MRALIBS)

libMADchem_la_SOURCES = correlationfactor.cc molecule.cc molecularbasis.cc \
                       corepotential.cc atomutil.cc lda.cc \
                       distpm.cc SCF.cc gth_pseudopotential.cc nemo.cc mp2.cc\
                       SCFOperators.cc xcfunctional_ldaonly.cc TDA.cc TDA_XC.cc \
                       $(thisinclude_HEADERS)
libMADchem_la_LDFLAGS = -version-info 0:0:0
                       
if MADNESS_HAS_LIBXC
   libMADchem_la_SOURCES += xcfunctional_libxc.cc
endif



libraries:	$(lib_LTLIBRARIES)
	for dir in $(SUBDIRS) ; do $(MAKE) -C $$dir $@ ; done

install-libraries: install-libLTLIBRARIES install-thisincludeHEADERS
	for dir in $(SUBDIRS) ; do $(MAKE) -C $$dir $@ ; done
