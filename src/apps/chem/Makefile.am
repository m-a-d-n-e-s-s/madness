include $(top_srcdir)/config/MakeGlobal.am

# default location for basis sets etc
AM_CPPFLAGS += -DMRA_CHEMDATA_DIR=\"$(abs_srcdir)\"

noinst_PROGRAMS = testxc

lib_LIBRARIES = libMADchem.a 

data_DATA = sto-3g sto-6g 6-31g coredata/mcp coredata/mcp_guess coredata/mcp2 coredata/mcp_guess


#TESTS = testxc

thisincludedir = $(includedir)/chem
thisinclude_HEADERS = correlationfactor.h molecule.h molecularbasis.h \
                       corepotential.h atomutil.h SCF.h xcfunctional.h


testxc_SOURCES = testxc.cc xcfunctional.h xcfunctional_ldaonly.cc lda.cc
testxc_LDADD = libMADchem.a $(LIBMRA) $(LIBWORLD)


libMADchem_a_SOURCES = correlationfactor.cc molecule.cc molecularbasis.cc \
                       corepotential.cc atomutil.cc lda.cc \
                       distpm.cc SCF.cc gth_pseudopotential.cc \
                       $(thisinclude_HEADERS)
                       
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
	
