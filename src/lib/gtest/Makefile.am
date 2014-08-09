include $(top_srcdir)/config/MakeGlobal.am

##AUTOMAKE_OPTIONS = --warnings=no-override

lib_LIBRARIES = libMADgtest.a

libMADgtest_a_SOURCES = gtest-all.cc

libMADgtest_a_CXXFLAGS = $(CXXFLAGS) -isystem ./include -I. -pthread 

