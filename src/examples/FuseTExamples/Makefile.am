include $(top_srcdir)/config/MakeGlobal.am
EXTRA_DIST = CMakeLists.txt

# specify the git commit version in the info functions
GITREV="MADNESS_GITREVISION=\"`git --git-dir=$(top_srcdir)/.git rev-parse HEAD`\" "

# default location for basis sets etc
AM_CPPFLAGS += -D$(GITREV)


noinst_PROGRAMS = CopyEx CompressEx InnerEx Test_FusedEx FusedEx ReconstructEx InnerExWeak multiplyEx
 
thisincludedir = $(includedir)/examples/FuseTExamples/
thisinclude_HEADERS = molecularmask.h  nonlinsol.h spectralprop.h

AUTOMAKE_OPTIONS = subdir-objects

CopyEx_SOURCES = CopyEx.cc
CompressEx_SOURCES = compress.cc
InnerEx_SOURCES = InnerEx.cc
ReconstructEx_SOURCES = ReconstructEx.cc
FusedEx_SOURCES = FusedEx.cc
testFusedEx_SOURCES = Test_FusedEx.cc
InnerExWeak_SOURCES = InnerEx_Weak.cc

LDADD = $(LIBCHEM) $(MRALIBS)

