include $(top_srcdir)/config/MakeGlobal.am

lib_LIBRARIES = libMADmisc.a

bin_PROGRAMS = interp3

TESTS = interp3

thisincludedir = $(includedir)/madness/misc
thisinclude_HEADERS = misc.h ran.h phandler.h cfft.h


libMADmisc_a_SOURCES = checksum_file.cc position_stream.cc gprofexit.cc ran.cc cfft.cc meminfo.cc\
                       misc.h ran.h cfft.h phandler.h

interp3_SOURCES = interp3.cc
interp3_LDADD = libMADmisc.a $(LIBWORLD)

libraries:	$(lib_LIBRARIES)
	for dir in $(SUBDIRS) ; do $(MAKE) -C $$dir $@ ; done

install-libraries: install-libLIBRARIES install-thisincludeHEADERS
	for dir in $(SUBDIRS) ; do $(MAKE) -C $$dir $@ ; done
