TOPDIR=..
ifndef QUOTE
  QUOTE=sed 's/ /\\ /g'
endif
ifndef SRCDIR
  SRCDIR=$(shell pwd | $(QUOTE))
endif

include $(SRCDIR)/$(TOPDIR)/lib/GlobalMakefile
include $(SRCDIR)/$(TOPDIR)/lib/MakeRules

SUBDIRS = lib bin

include $(SRCDIR)/$(TOPDIR)/lib/MakeSubDirs
include $(TOPDIR)/Makedirlist

