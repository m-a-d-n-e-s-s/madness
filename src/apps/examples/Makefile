TOPDIR = ../../..

CXXSRC = hello.cc hatom_energy.cc

# GlobalMakefile now does most of the work for you
include $(TOPDIR)/config/GlobalMakefile

all:	examples

examples:	hello hatom_energy

# You must define additional rules to link your test programs
hello:	hello.o $(THISLIBDEPEND) 
	$(LTLINK) $(LD) $(LDFLAGS) -o $@ $^ $(SYSLIBS) $(LTLINKBINOPTS)

# You must define additional rules to link your test programs
hatom_energy:	hatom_energy.o $(THISLIBDEPEND) 
	$(LTLINK) $(LD) $(LDFLAGS) -o $@ $^ $(SYSLIBS) $(LTLINKBINOPTS)
