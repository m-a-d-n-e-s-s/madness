# This must match definition of AC_CONFIG_MACRO_DIR in configure.ac
ACLOCAL_AMFLAGS = -I config

SUBDIRS = include src

doc:
	cd doc; $(MAKE)


LIBRARIES:	$(lib_LIBRARIES)
	for dir in $(SUBDIRS) ; do $(MAKE) -C $$dir $@ ; done
