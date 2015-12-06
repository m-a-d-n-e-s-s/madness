include $(top_srcdir)/config/MakeGlobal.am
EXTRA_DIST = CMakeLists.txt

lib_LTLIBRARIES = libMADtinyxml.la

thisincludedir = $(includedir)/madness/external/tinyxml
thisinclude_HEADERS = tinystr.h tinyxml.h

libMADtinyxml_la_SOURCES = tinystr.cc tinyxml.cc tinyxmlerror.cc tinyxmlparser.cc \
                          tinystr.h  tinyxml.h
libMADtinyxml_la_LDFLAGS = -version-info 0:0:0

AM_CPPFLAGS += -I.

libraries:	$(lib_LTLIBRARIES)

install-libraries: install-libLTLIBRARIES install-thisincludeHEADERS
