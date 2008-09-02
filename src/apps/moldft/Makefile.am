include $(top_srcdir)/config/MakeGlobal.am

SUBDIRS = xc

bin_PROGRAMS = moldft

moldft_SOURCES = moldft.cc lda.cc  molecule.cc \
                 \
                molecule.h molecularbasis.h

moldft_LDADD = xc/libMADxc.a $(MRALIBS)

