TOPDIR = .

SUBDIRS = src

include $(TOPDIR)/lib/GlobalMakefile
include $(TOPDIR)/lib/MakeSubDirs
include $(TOPDIR)/Makedirlist

clean::
	-rm -f depcheck.cc

distclean:: clean
	-rm -f config.log
	-rm -f config.status
	-rm -f libtool

