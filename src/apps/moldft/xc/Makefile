MY_ROOT=.
DEBUG=YES
#HPM_PROFILING=YES

TARGET=BGP
TARGET=GCC


#############################################

ifeq ($(TARGET),BGP)


   RM = rm
   RMFLAGS = -f

   AR = ar
   ARFLAGS = cru
   HAVE_BGPIBM = -DHAVE_IBMBGP

   LIB=-L$(MY_ROOT) 
   INC=-I$(MY_ROOT) 

   CC=mpixlc_r
   CXX=mpixlcxx_r
   FC=mpixlf77_r
   COPT= -O5 -qarch=450d -qtune=450 -qthreaded -qpath=ILbc:/soft/apps/xlC-interim-fix-virtualthunk/exe/
   FOPT= -O5 -qarch=450d -qtune=450 -qthreaded -qpath=ILbc:/soft/apps/xlC-interim-fix-virtualthunk/exe/
   CXXFLAGS=$(INC) $(COPT) $(HAVE_BGPIBM)
   CPP=$(CC) -E

   LD=$(CC)
   LDFLAGS= -O3 -qarch=450d -qtune=450 -qthreaded

endif

ifeq ($(TARGET),GCC)


   RM = rm
   RMFLAGS = -f

   AR = ar
   ARFLAGS = cru

   LIB=-L$(MY_ROOT) 
   INC=-I$(MY_ROOT) 

   CC=gcc
   CXX=g++
   FC=gfortran
   COPT= -O3 -Wall 
   FFLAGS= -O3 -Wall 
   CXXFLAGS=$(INC) $(COPT) $(HAVE_BGPIBM)
   CPP=$(CC) -E

   LD=$(CC)
   LDFLAGS= -O3 

endif



OBJECTS =  c_ft97.o \
           c_lyp.o \
           c_p86.o \
           c_pbe.o \
           c_pw91.o \
           c_pw92.o \
           c_pz81.o \
           c_vwn5.o \
           c_vwn5rpa.o \
           x_b3.o \
           x_b88.o \
           x_ft97b.o \
           x_lda.o \
           x_pbe.o \
           x_pw91.o \
           xc_b3lyp.o \
           xc_b97.o \
           xc_b97_1.o \
           xc_b97_2.o \
           xc_edf1.o \
           xc_ft97.o \
           xc_hcth.o \
           xc_hcth120.o \
           xc_hcth147.o \
           xc_hcth407.o \
           xc_pbe.o \
           xc_pw91.f



all: libMADxc.a

refresh: realclean all


libMADxc.a: $(OBJECTS)
	rm -rf libMADxc.a
	$(AR) $(ARFLAGS) libMADxc.a $(OBJECTS) 
	ranlib libMADxc.a

$%.o: %.f
	$(FC) -c $(FFLAGS) $(INC) $< -o $@

clean:
	$(RM) $(RMFLAGS) *.o

realclean: clean
	$(RM) $(RMFLAGS) *.a



