include $(top_srcdir)/config/MakeGlobal.am

# default location for basis sets etc
AM_CPPFLAGS += -DMRA_CHEMDATA_DIR=\"$(abs_srcdir)\"

noinst_PROGRAMS = testxc plotxc test_SCFOperators

lib_LIBRARIES = libMADchem.a 

data_DATA = sto-3g sto-6g 6-31g coredata/mcp coredata/mcp2 coredata/mcp_guess
dist_data_DATA = $(data_DATA)

TESTS = test_SCFOperators # testxc 

thisincludedir = $(includedir)/chem

thisinclude_HEADERS = correlationfactor.h molecule.h molecularbasis.h \
                      corepotential.h atomutil.h SCF.h xcfunctional.h \
                      mp2.h nemo.h potentialmanager.h gth_pseudopotential.h \
                      molecular_optimizer.h projector.h TDA.h TDA_XC.h \
                      TDA_guess.h TDA_exops.h SCFOperators.h

testxc_SOURCES = testxc.cc xcfunctional.h xcfunctional_ldaonly.cc lda.cc
testxc_LDADD = libMADchem.a $(LIBMRA) $(LIBWORLD)

test_SCFOperators_SOURCES = test_SCFOperators.cc 
test_SCFOperators_LDADD = libMADchem.a $(MRALIBS)

plotxc_SOURCES = plotxc.cc xcfunctional.h xcfunctional_ldaonly.cc lda.cc
plotxc_LDADD = libMADchem.a $(LIBMRA) $(LIBWORLD)

libMADchem_a_SOURCES = correlationfactor.cc molecule.cc molecularbasis.cc \
                       corepotential.cc atomutil.cc lda.cc \
                       distpm.cc SCF.cc gth_pseudopotential.cc nemo.cc mp2.cc\
                       SCFOperators.cc $(thisinclude_HEADERS)
                       
if MADNESS_HAS_LIBXC
   libMADchem_a_SOURCES += xcfunctional_libxc.cc
else
   libMADchem_a_SOURCES += xcfunctional_ldaonly.cc 
endif



libraries:	$(lib_LIBRARIES)
	for dir in $(SUBDIRS) ; do $(MAKE) -C $$dir $@ ; done

install-libraries: install-libLIBRARIES install-thisincludeHEADERS
	for dir in $(SUBDIRS) ; do $(MAKE) -C $$dir $@ ; done


LDADD = $(LIBMRA) $(LIBTENSOR) $(LIBWORLD)
	
