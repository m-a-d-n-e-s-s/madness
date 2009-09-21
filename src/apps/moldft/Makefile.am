include $(top_srcdir)/config/MakeGlobal.am

SUBDIRS = xc

bin_PROGRAMS = testpg moldft 

moldft_SOURCES = moldft.cc lda.cc  molecule.cc \
                 molecule.h molecularbasis.h

testpg_SOURCES = testpg.cc pointgroup.h

moldft_LDADD = xc/libMADxc.a $(MRALIBS)

testpg_LDADD = $(MRALIBS)

