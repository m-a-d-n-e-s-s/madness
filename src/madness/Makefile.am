SUBDIRS = world tinyxml muParser misc tensor mra

thisincludedir = $(includedir)/madness
thisinclude_HEADERS = mad_types.h constants.h fortran_ctypes.h madness_config.h TAU.h

libraries:
	for dir in $(SUBDIRS) ; do $(MAKE) -C $$dir libraries ; done

install-libraries: install-thisincludeHEADERS
	for dir in $(SUBDIRS) ; do $(MAKE) -C $$dir install-libraries ; done
