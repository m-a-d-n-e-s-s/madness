SUBDIRS = madness apps examples
EXTRA_DIST = CMakeLists.txt

thisincludedir = $(includedir)
thisinclude_HEADERS = madness.h

libraries:
	for dir in madness apps ; do $(MAKE) -C $$dir libraries ; done

install-libraries: install-thisincludeHEADERS
	for dir in madness apps ; do $(MAKE) -C $$dir install-libraries ; done
