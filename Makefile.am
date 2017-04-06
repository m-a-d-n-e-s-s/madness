# This must match definition of AC_CONFIG_MACRO_DIR in configure.ac
ACLOCAL_AMFLAGS = -I config

SUBDIRS = src config doc
EXTRA_DIST = cmake external CMakeLists.txt

docs:
	$(MAKE) -C doc docs


.PHONY:

libraries:
	$(MAKE) -C src libraries

world:
	$(MAKE) -C src/madness/world libraries

tinyxml:
	$(MAKE) -C src/madness/external/tinyxml libraries

muparser:
	$(MAKE) -C src/madness/external/muParser libraries

misc: world
	$(MAKE) -C src/madness/misc libraries
	
tensor linalg: misc
	$(MAKE) -C src/madness/tensor libraries

mra: linalg muparser
	$(MAKE) -C src/madness/mra libraries
	
chem: mra
	$(MAKE) -C src/apps/chem libraries

install-libraries:
	$(MAKE) -C src install-libraries
	$(MAKE) -C config install
	
install-madinclude:
	$(MAKE) -C src/madness install-thisincludeHEADERS

install-world: install-madinclude
	$(MAKE) -C src/madness/world install-libraries

install-tinyxml: install-madinclude
	$(MAKE) -C src/madness/external/tinyxml install-libraries

install-muparser: install-madinclude
	$(MAKE) -C src/madness/external/muParser install-libraries

install-misc: install-world
	$(MAKE) -C src/madness/misc install-libraries

install-tensor: install-misc
	$(MAKE) -C src/madness/tensor install-libraries

install-mra: install-tensor install-muparser install-tinyxml
	$(MAKE) -C src install-thisincludeHEADERS
	$(MAKE) -C src/madness/mra install-libraries
	$(MAKE) -C config install

install-chem: install-libraries

