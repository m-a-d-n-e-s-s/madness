include $(top_srcdir)/config/MakeGlobal.am

bin_PROGRAMS = testpg moldft mcpfit testperiodic

moldft_SOURCES = moldft.cc lda.cc  molecule.cc \
                 molecule.h molecularbasis.h molecularbasis.cc \
                 atomutil.h atomutil.cc \
                 corepotential.h corepotential.cc \
		 xcfunctional.h xcfunctional_ldaonly.cc xcfunctional_libxc.cc 

testpg_SOURCES = testpg.cc pointgroup.h

testperiodic_SOURCES = testperiodic.cc lda.cc  molecule.cc \
                 molecule.h molecularbasis.h molecularbasis.cc \
                 atomutil.h atomutil.cc \
                 corepotential.h corepotential.cc \
		 xcfunctional.h xcfunctional_ldaonly.cc xcfunctional_libxc.cc 

mcpfit_SOURCES = mcpfit.cc atomutil.h atomutil.cc corepotential.h corepotential.cc

if HAVE_LIBXC
bin_PROGRAMS += testxc
testxc_SOURCES = testxc.cc xcfunctional.h xcfunctional_libxc.cc
testxc_LDADD = $(MRALIBS) libxc.a
endif

moldft_LDADD = $(MRALIBS)

testpg_LDADD = $(MRALIBS)

mcpfit_LDADD = $(MRALIBS)


testperiodic_LDADD = $(MRALIBS)

