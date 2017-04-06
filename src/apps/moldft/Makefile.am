include $(top_srcdir)/config/MakeGlobal.am
EXTRA_DIST = CMakeLists.txt
AUTOMAKE_OPTIONS = subdir-objects

if MADNESS_HAS_LIBXC
LDADD = $(LIBCHEM) $(MRALIBS) -lxc
else
LDADD = $(LIBCHEM) $(MRALIBS)
endif

noinst_PROGRAMS = testpg testperiodic
bin_PROGRAMS = moldft mcpfit
#bin_PROGRAMS = testpg moldft mcpfit testperiodic preal

moldft_SOURCES = moldft.cc


#preal_SOURCES = preal.cc

testpg_SOURCES = testpg.cc pointgroup.h

testperiodic_SOURCES = testperiodic.cc

mcpfit_SOURCES = mcpfit.cc

