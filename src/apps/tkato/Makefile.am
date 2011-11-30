include $(top_srcdir)/config/MakeGlobal.am

bin_PROGRAMS = moldft plotao

moldft_SOURCES = moldft.cc scfparam.h scf.h scf.cc lda.cc  molecule.cc \
                 molecule.h molecularbasis.h molecularbasis.cc \
                 atomutil.h atomutil.cc \
                 corepotential.h corepotential.cc \
                 xcfunctional.h xcfunctional_ldaonly.cc xcfunctional_libxc.cc

plotao_SOURCES = plotao.cc

moldft_LDADD = $(MRALIBS)

plotao_LDADD = $(MRALIBS)

