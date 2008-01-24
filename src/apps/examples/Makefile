TOPDIR = ../../..


BIN_OR_LIB = BIN
TARGET_TO_MAKE = all
BINOBJ = hello.o hatom_energy.o molecule.o hf.o bigboy.o
TESTPROGS = hello hatom_energy hf bigboy

# GlobalMakefile now does most of the work for you
include $(TOPDIR)/config/GlobalMakefile

all:	$(TESTPROGS)

hello:	hello.o $(THISLIBDEPEND) 
	$(LTLINK) $(LD) $(LDFLAGS) -o $@ $^ $(SYSLIBS) $(LTLINKBINOPTS)

hatom_energy:	hatom_energy.o $(THISLIBDEPEND) 
	$(LTLINK) $(LD) $(LDFLAGS) -o $@ $^ $(SYSLIBS) $(LTLINKBINOPTS)

hf:	hf.o molecule.o $(THISLIBDEPEND) 
	$(LTLINK) $(LD) $(LDFLAGS) -o $@ $^ $(SYSLIBS) $(LTLINKBINOPTS)

bigboy:	bigboy.o molecule.o $(THISLIBDEPEND) 
	$(LTLINK) $(LD) $(LDFLAGS) -o $@ $^ $(SYSLIBS) $(LTLINKBINOPTS)

