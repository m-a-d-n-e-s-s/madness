include $(top_srcdir)/config/MakeGlobal.am


if HAVE_LIBXC
LDADD = $(MRALIBS) -lxc
else
LDADD = $(MRALIBS)
endif

bin_PROGRAMS = testpg moldft mcpfit testperiodic testxc 

moldft_SOURCES = moldft.cc  molecule.cc lda.cc \
                 molecule.h molecularbasis.h molecularbasis.cc \
                 atomutil.h atomutil.cc \
                 corepotential.h corepotential.cc \
                 xcfunctional.h xcfunctional.cc 

testpg_SOURCES = testpg.cc pointgroup.h

testperiodic_SOURCES = testperiodic.cc molecule.cc lda.cc \
                 molecule.h molecularbasis.h molecularbasis.cc \
                 atomutil.h atomutil.cc \
                 corepotential.h corepotential.cc \
                 xcfunctional.h xcfunctional.cc

mcpfit_SOURCES = mcpfit.cc atomutil.h atomutil.cc corepotential.h corepotential.cc

testxc_SOURCES = testxc.cc xcfunctional.h

#polar_SOURCES  = polar.cc molecule.cc \
#                 molecule.h molecularbasis.h molecularbasis.cc \
#                 atomutil.h atomutil.cc \
#                 corepotential.h corepotential.cc \
#                 xcfunctional.h xcfunctional.cc 

if HAVE_LIBXC
moldft_SOURCES += xcfunctional_libxc.cc
testperiodic_SOURCES += xcfunctional_libxc.cc lda.cc
testxc_SOURCES += xcfunctional_libxc.cc
else
testxc_SOURCES += xcfunctional.cc 
endif

