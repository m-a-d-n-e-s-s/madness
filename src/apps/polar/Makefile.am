include $(top_srcdir)/config/MakeGlobal.am
EXTRA_DIST = CMakeLists.txt


if MADNESS_HAS_LIBXC
LDADD = $(MRALIBS) -lxc
else
LDADD = $(MRALIBS)
endif

#bin_PROGRAMS = testpg moldft mcpfit testperiodic testxc
bin_PROGRAMS = moldft

moldft_SOURCES = moldft.cc lda.cc  molecule.cc \
                 molecule.h molecularbasis.h molecularbasis.cc \
                 atomutil.h atomutil.cc \
                 corepotential.h corepotential.cc \
                 xcfunctional.h xcfunctional_ldaonly.cc  distpm.cc


if MADNESS_HAS_LIBXC
moldft_SOURCES += xcfunctional_libxc.cc
endif

