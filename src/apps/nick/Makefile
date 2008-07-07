TOPDIR = ../../..
FC = gfortran  #overiding the default Fortran Compiler
BIN_OR_LIB = BIN
TARGET_TO_MAKE = hello
#TESTPROGS = hatom_energy
OBJECTS = hello.o toms707.o wavef.o coulcc36-f90.o 
TESTLIB = -lgsl -lgslcblas -lgfortran 


# GlobalMakefile now does most of the work for you
include $(TOPDIR)/config/GlobalMakefile

hatom_energy:	$(OBJECTS) $(THISLIBDEPEND) 
	$(LTLINK) $(LD) $(LDFLAGS) -o $@ $^ $(SYSLIBS) $(LTLINKBINOPTS)

hello:	$(OBJECTS) $(THISLIBDEPEND) 
	$(LTLINK) $(LD) $(LDFLAGS) -o $@ $^ $(SYSLIBS) $(LTLINKBINOPTS) $(TESTLIB)

test:	$(TESTOBJ)
	g++ $(TESTLIB) f11tester.cc $(TESTOBJ) -o app

toms707.o:
	$(FC) toms707.f -c

#coul36-f90.o:
#	$(FC) coul36-f90.f -c
