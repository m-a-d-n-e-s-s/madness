include $(top_srcdir)/config/MakeGlobal.am
EXTRA_DIST = CMakeLists.txt

lib_LTLIBRARIES = libMADmisc.la

noinst_PROGRAMS = interp3

TESTS = interp3

thisincludedir = $(includedir)/madness/misc
thisinclude_HEADERS = misc.h ran.h phandler.h interpolation_1d.h cfft.h


libMADmisc_la_SOURCES = checksum_file.cc position_stream.cc gprofexit.cc ran.cc cfft.cc\
                       misc.h ran.h cfft.h phandler.h
libMADmisc_la_LDFLAGS = -version-info 0:0:0

interp3_SOURCES = interp3.cc
interp3_LDADD = libMADmisc.la $(LIBWORLD)

libraries:	$(lib_LTLIBRARIES)
	for dir in $(SUBDIRS) ; do $(MAKE) -C $$dir $@ ; done

install-libraries: install-libLTLIBRARIES install-thisincludeHEADERS
	for dir in $(SUBDIRS) ; do $(MAKE) -C $$dir $@ ; done
