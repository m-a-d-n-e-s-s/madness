TOPDIR =  .

SUBDIRS = src

include $(TOPDIR)/config/GlobalMakefile
include $(TOPDIR)/config/MakeSubDirs
include $(TOPDIR)/Makedirlist

all:	default

etags:
	etags `find . \( -name "*.h" -o -name "*.cc" -o -name "*.c" -o -name "*.hpp" -o -name "*.cpp" \)`

setsvnkeyword:
	-find . \( -name "*.h" -o -name "*.cc" -o -name "*.c" -o -name "*.hpp" -o -name "*.cpp" \) -exec config/setsvnkey "{}" ";"

clean::
	-rm -f depcheck.cc *~ lib/* bin/*

configure:	configure.in
	aclocal -I config/autoconf
	autoconf

distclean:: clean
	-rm -f config.log config.status
	-rm -rf autom4te.cache aclocal.m4
	-rm -f libtool
	-rm -f config/MakeVars config/*~


