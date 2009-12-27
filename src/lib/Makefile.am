SUBDIRS = world tinyxml muParser misc tensor linalg mra

libraries:	$(lib_LIBRARIES)
	for dir in $(SUBDIRS) ; do $(MAKE) -C $$dir $@ ; done



