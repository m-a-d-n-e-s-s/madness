TOPDIR=.
ifndef QUOTE
  QUOTE=sed 's/ /\\ /g'
endif
ifndef SRCDIR
  SRCDIR=$(shell pwd | $(QUOTE))
endif

include $(SRCDIR)/$(TOPDIR)/lib/GlobalMakefile
include $(SRCDIR)/$(TOPDIR)/lib/MakeRules

SUBDIRS = src

include $(SRCDIR)/$(TOPDIR)/lib/MakeSubDirs
include $(TOPDIR)/Makedirlist

clean::
	-rm -f depcheck.cc

distclean:: clean
	-rm -f config.log
	-rm -f config.status
	-rm -f libtool

