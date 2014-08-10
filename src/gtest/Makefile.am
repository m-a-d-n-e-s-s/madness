include $(top_srcdir)/config/MakeGlobal.am

if MADNESS_HAS_GOOGLE_TEST
##AUTOMAKE_OPTIONS = --warnings=no-override

noinst_LIBRARIES = libMADgtest.a

libMADgtest_a_SOURCES = gtest-all.cc
libMADgtest_a_CPPFLAGS = $(LIBGTEST_CPPFLAGS) -I$(srcdir)
libMADgtest_a_CXXFLAGS = $(LIBGTEST_CXXFLAGS)
libMADgtest_a_LIBADD = $(LIBGTEST_LIBS)

endif
