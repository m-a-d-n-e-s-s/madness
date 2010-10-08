include $(top_srcdir)/config/MakeGlobal.am

bin_PROGRAMS = testpg moldft

moldft_SOURCES = moldft.cc lda.cc  molecule.cc \
                 molecule.h molecularbasis.h molecularbasis.cc \
                 atomutil.h atomutil.cc \
                 corepotential.h corepotential.cc

testpg_SOURCES = testpg.cc pointgroup.h

moldft_LDADD = $(MRALIBS)

testpg_LDADD = $(MRALIBS)

