include $(top_srcdir)/config/MakeGlobal.am
SUBDIRS = external world misc tensor mra
EXTRA_DIST = CMakeLists.txt

thisincludedir = $(includedir)/madness
thisinclude_HEADERS = constants.h fortran_ctypes.h madness_config.h config.h 

lib_LTLIBRARIES = libmadness.la
libmadness_la_LIBADD = $(MRALIBS)

# Add a dummy source file so the linker behaves properly.
libmadness_la_SOURCES = madness.cc

libraries: $(lib_LTLIBRARIES)
	for dir in $(SUBDIRS) ; do $(MAKE) -C $$dir libraries ; done

install-libraries: install-libLTLIBRARIES install-thisincludeHEADERS
	for dir in $(SUBDIRS) ; do $(MAKE) -C $$dir install-libraries ; done
