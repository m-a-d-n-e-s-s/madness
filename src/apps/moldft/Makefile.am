include $(top_srcdir)/config/MakeGlobal.am


if MADNESS_HAS_LIBXC
LDADD = $(MRALIBS) -lxc
else
LDADD = $(MRALIBS)
endif

bin_PROGRAMS = testpg moldft mcpfit testperiodic testxc

moldft_SOURCES = moldft.cc lda.cc  molecule.cc \
                 molecule.h molecularbasis.h molecularbasis.cc \
                 atomutil.h atomutil.cc \
                 corepotential.h corepotential.cc \
                 xcfunctional.h xcfunctional_ldaonly.cc 

testpg_SOURCES = testpg.cc pointgroup.h

testperiodic_SOURCES = testperiodic.cc lda.cc  molecule.cc \
                 molecule.h molecularbasis.h molecularbasis.cc \
                 atomutil.h atomutil.cc \
                 corepotential.h corepotential.cc \
                 xcfunctional.h xcfunctional_ldaonly.cc

mcpfit_SOURCES = mcpfit.cc atomutil.h atomutil.cc corepotential.h corepotential.cc

testxc_SOURCES = testxc.cc xcfunctional.h

if MADNESS_HAS_LIBXC
moldft_SOURCES += xcfunctional_libxc.cc
testperiodic_SOURCES += xcfunctional_libxc.cc
testxc_SOURCES += xcfunctional_libxc.cc
else
testxc_SOURCES += xcfunctional_ldaonly.cc lda.cc
endif

