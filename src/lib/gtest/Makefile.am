include $(top_srcdir)/config/MakeGlobal.am

if MADNESS_HAS_GOOGLE_TEST
##AUTOMAKE_OPTIONS = --warnings=no-override

lib_LIBRARIES = libMADgtest.a

libMADgtest_a_SOURCES = gtest-all.cc

libMADgtest_a_CXXFLAGS = $(GTEST_CXXFLAGS) -isystem ./include -I. -pthread 

libraries:	$(lib_LIBRARIES)
else
libraries:
	@echo "No libraries to build"
endif
