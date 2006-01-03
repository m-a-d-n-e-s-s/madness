TOPDIR=..
ifndef QUOTE
  QUOTE=sed 's/ /\\ /g'
endif
ifndef SRCDIR
  SRCDIR=$(shell pwd | $(QUOTE))
endif

include $(SRCDIR)/$(TOPDIR)/lib/GlobalMakefile
include $(TOPDIR)/lib/MakeVars

.PHONY: all
all:: html
ifneq ($(PDFLATEX),NO)
  all:: latex
endif

.PHONY: html
ifneq ($(DOXYGEN),NO)
html:
	$(DOXYGEN) doxygen.cfg
endif

.PHONY: latex
ifneq ($(PDFLATEX),NO)
latex:
	cd latex && $(MAKE)
endif

.PHONY: veryclean
veryclean:
	/bin/rm -rf html
	/bin/rm -rf latex
	/bin/rm -rf *~

.PHONY: clean
clean:
	cd latex && $(MAKE) clean
	/bin/rm -rf *~

install:
	$(INSTALL) $(INSTALLDIROPT) $(installroot)$(prefix)
	/bin/cp -r html $(installroot)$(prefix)

install_devel:

