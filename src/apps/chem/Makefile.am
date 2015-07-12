include $(top_srcdir)/config/MakeGlobal.am

# default location for basis sets etc
AM_CPPFLAGS += -DMRA_CHEMDATA_DIR=\"$(abs_srcdir)\"

noinst_PROGRAMS = testxc plotxc

lib_LTLIBRARIES = libMADchem.la
noinst_LTLIBRARIES = libxcfunctional_ldaonly.la

data_DATA = sto-3g sto-6g 6-31g coredata/mcp coredata/mcp2 coredata/mcp_guess
dist_data_DATA = $(data_DATA)

#TESTS = testxc

thisincludedir = $(includedir)/chem

thisinclude_HEADERS = correlationfactor.h molecule.h molecularbasis.h \
                      corepotential.h atomutil.h SCF.h xcfunctional.h \
                      mp2.h nemo.h potentialmanager.h gth_pseudopotential.h \
                      molecular_optimizer.h projector.h TDA.h TDA_XC.h \
                      TDA_guess.h TDA_exops.h


testxc_SOURCES = testxc.cc xcfunctional.h
testxc_LDADD = libMADchem.la $(LIBMRA) $(LIBWORLD)

plotxc_SOURCES = plotxc.cc xcfunctional.h
plotxc_LDADD = libMADchem.la $(LIBMRA) $(LIBWORLD)

libMADchem_la_SOURCES = correlationfactor.cc molecule.cc molecularbasis.cc \
                       corepotential.cc atomutil.cc lda.cc \
                       distpm.cc SCF.cc gth_pseudopotential.cc nemo.cc mp2.cc\
                       $(thisinclude_HEADERS)
libMADchem_la_LDFLAGS = -version-info 0:0:0

libxcfunctional_ldaonly_la_SOURCES = xcfunctional_ldaonly.cc
                       
if MADNESS_HAS_LIBXC
   testxc_LDADD += libxcfunctional_ldaonly.la
   plotxc_LDADD += libxcfunctional_ldaonly.la
   libMADchem_la_SOURCES += xcfunctional_libxc.cc
else
   libMADchem_la_LIBADD = libxcfunctional_ldaonly.la 
endif



libraries:	$(lib_LTLIBRARIES)
	for dir in $(SUBDIRS) ; do $(MAKE) -C $$dir $@ ; done

install-libraries: install-libLTLIBRARIES install-thisincludeHEADERS
	for dir in $(SUBDIRS) ; do $(MAKE) -C $$dir $@ ; done


LDADD = $(LIBMRA) $(LIBTENSOR) $(LIBWORLD)
	
