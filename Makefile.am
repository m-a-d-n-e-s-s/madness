# This must match definition of AC_CONFIG_MACRO_DIR in configure.ac
ACLOCAL_AMFLAGS = -I ./config

SUBDIRS = include src

doc:
	cd doc; $(MAKE)

	
.PHONY:

libraries:	$(lib_LIBRARIES)
	for dir in $(SUBDIRS) ; do $(MAKE) -C $$dir $@ ; done

world:
	$(MAKE) -C src/lib/world libraries

tinyxml:
	$(MAKE) -C src/lib/tinyxml libraries

muparser:
	$(MAKE) -C src/lib/muParser libraries

misc: world
	$(MAKE) -C src/lib/misc libraries
	
tensor: misc
	$(MAKE) -C src/lib/tensor libraries

linalg: tensor  
	$(MAKE) -C src/lib/linalg libraries

mra: linalg tinyxml muparser
	$(MAKE) -C src/lib/mra libraries

install-libraries:
	for dir in $(SUBDIRS) ; do $(MAKE) -C $$dir $@ ; done
	
install-madinclude:
	$(MAKE) -C include install-libraries

install-world: world install-madinclude
	$(MAKE) -C src/lib/world install-libraries

install-tinyxml: tinyxml install-madinclude
	$(MAKE) -C src/lib/tinyxml install-libraries

install-muparser: muparser install-madinclude
	$(MAKE) -C src/lib/muParser install-libraries

install-misc: install-world misc
	$(MAKE) -C src/lib/misc install-libraries

install-tensor: install-misc tensor
	$(MAKE) -C src/lib/tensor install-libraries

install-linalg: install-tensor linalg
	$(MAKE) -C src/lib/linalg install-libraries

install-mra: install-linalg install-muparser install-tinyxml mra
	$(MAKE) -C src/lib/mra install-libraries

