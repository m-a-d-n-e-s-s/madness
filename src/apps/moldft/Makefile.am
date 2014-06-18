include $(top_srcdir)/config/MakeGlobal.am
AUTOMAKE_OPTIONS = subdir-objects



if MADNESS_HAS_LIBXC
LDADD = $(MRALIBS) -lxc
else
LDADD = $(MRALIBS)
endif

bin_PROGRAMS = testpg moldft mcpfit testperiodic 

moldft_SOURCES = moldft.cc

testpg_SOURCES = testpg.cc pointgroup.h

testperiodic_SOURCES = testperiodic.cc

mcpfit_SOURCES = mcpfit.cc

if MADNESS_HAS_LIBXC
moldft_SOURCES += ../../lib/chem/xcfunctional_libxc.cc
testperiodic_SOURCES += ../../lib/chem/xcfunctional_libxc.cc
endif

