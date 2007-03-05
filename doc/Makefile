TOPDIR=..

##include $(TOPDIR)/config/GlobalMakefile
include $(TOPDIR)/config/MakeVars

.PHONY: all
all:
	$(DOXYGEN) doxygen.cfg


.PHONY: clean
clean::
	/bin/rm -rf html latex
	/bin/rm -rf *~

veryclean: clean

distclean: clean
	/bin/rm doxygen.cfg

install::
	$(INSTALL) $(INSTALLDIROPT) $(installroot)$(prefix)
	/bin/cp -r html $(installroot)$(prefix)

install_devel::

