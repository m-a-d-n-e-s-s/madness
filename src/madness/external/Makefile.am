SUBDIRS = gtest tinyxml muParser
EXTRA_DIST = CMakeLists.txt elemental

libraries:
	for dir in $(SUBDIRS) ; do $(MAKE) -C $$dir libraries ; done

install-libraries: install-thisincludeHEADERS
	for dir in $(SUBDIRS) ; do $(MAKE) -C $$dir install-libraries ; done
