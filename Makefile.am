# This must match definition of AC_CONFIG_MACRO_DIR in configure.ac
ACLOCAL_AMFLAGS = -I config

SUBDIRS = include src

doc:
	cd doc; $(MAKE)


libraries:	$(lib_LIBRARIES)
	for dir in $(SUBDIRS) ; do $(MAKE) -C $$dir $@ ; done


.PHONY: world install-world
world:
	cd src/lib/world && $(MAKE) libraries

install-world: world
	cd src/lib/world && $(MAKE) install-libLIBRARIES install-thisincludeHEADERS
