TOPDIR = ../../..


BIN_OR_LIB = BIN
TARGET_TO_MAKE = examples
BINOBJ = hello.o hatom_energy.o
TESTPROGS = hello hatom_energy

# GlobalMakefile now does most of the work for you
include $(TOPDIR)/config/GlobalMakefile

default::	$(TESTPROGS)

hello:	hello.o $(THISLIBDEPEND) 
	$(LTLINK) $(LD) $(LDFLAGS) -o $@ $^ $(SYSLIBS) $(LTLINKBINOPTS)

hatom_energy:	hatom_energy.o $(THISLIBDEPEND) 
	$(LTLINK) $(LD) $(LDFLAGS) -o $@ $^ $(SYSLIBS) $(LTLINKBINOPTS)
