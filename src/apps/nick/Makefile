TOPDIR = ../../..
FC = ftn

BIN_OR_LIB = BIN
TARGET_TO_MAKE = all
BINOBJ = hello.o toms707.o wavef.o coulcc36-f90.o 
TESTPROGS = hello
GSL = -I$(TOPDIR)/include/myGSL/include -L$(TOPDIR)/include/myGSL/lib 

# GlobalMakefile now does most of the work for you
include $(TOPDIR)/config/GlobalMakefile

all: $(TESTPROGS)

hello:	$(BINOBJ) $(THISLIBDEPEND) 
	 $(LTLINK) $(LD) $(LDFLAGS) $(GSL) -o $@ $^ $(SYSLIBS) $(LTLINKBINOPTS) 

wavef.o:	wavef.cc $(THISLIBDEPEND)
	$(LTLINK) $(LD) $(INCLUDE) $(GSL) -c $^ $(SYSLIBS) $(LTLINKBINOPTS)

#wavef.o:	wavef.cc
#CC  -Wall -Wno-strict-aliasing -ansi -O2 -ffast-math -fomit-frame-pointer -mfpmath=sse -msse -msse2 -msse3  -I../../../include -I../../../src/lib -I../../../include -I../../../src/lib -I../../../src/apps  -I../../../include/myGSL/include/ -L../../../include/myGSL/lib -lgsl -lgslcblas -c wavef.cc -o wavef.o -llapack -lblas -lgsl -lgslcblas 
