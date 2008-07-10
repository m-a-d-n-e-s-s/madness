TOPDIR = ../../..


BIN_OR_LIB = BIN
TARGET_TO_MAKE = all
BINOBJ = tdse.o
TESTPROGS =  tdse

# GlobalMakefile now does most of the work for you
include $(TOPDIR)/config/GlobalMakefile

all:	$(TESTPROGS)

tdse:	tdse.o $(THISLIBDEPEND) 
	$(LTLINK) $(LD) $(LDFLAGS) -o $@ $^ $(SYSLIBS) $(LTLINKBINOPTS)

