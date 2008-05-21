TOPDIR = ../../..

SUBDIRS = xc
BIN_OR_LIB = BIN
TARGET_TO_MAKE = all
BINOBJ = moldft.o molecule.o lda.o
TESTPROGS = moldft

# GlobalMakefile now does most of the work for you
include $(TOPDIR)/config/GlobalMakefile
include $(TOPDIR)/config/MakeSubDirs
include $(TOPDIR)/Makedirlist

all:	xc $(TESTPROGS)

xc:
	cd xc; make

moldft:	moldft.o molecule.o lda.o $(THISLIBDEPEND) 
	$(LTLINK) $(LD) $(LDFLAGS) -o $@ $^ $(SYSLIBS) $(LTLINKBINOPTS)

