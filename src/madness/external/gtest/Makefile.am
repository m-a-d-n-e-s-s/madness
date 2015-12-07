include $(top_srcdir)/config/MakeGlobal.am
EXTRA_DIST = CMakeLists.txt

if MADNESS_HAS_GOOGLE_TEST
##AUTOMAKE_OPTIONS = --warnings=no-override

noinst_LTLIBRARIES = libMADgtest.la

libMADgtest_la_SOURCES = gtest-all.cc
libMADgtest_la_CPPFLAGS = $(LIBGTEST_CPPFLAGS) -I$(srcdir)
libMADgtest_la_CXXFLAGS = $(LIBGTEST_CXXFLAGS)
libMADgtest_la_LIBADD = $(LIBGTEST_LIBS)

endif
