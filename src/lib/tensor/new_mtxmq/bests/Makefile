MY_ROOT=.
DEBUG=YES
#HPM_PROFILING=YES

TARGET=BGP


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
   COPT= -O5 -qarch=450d -qtune=450 -qthreaded -qpath=ILbc:/soft/apps/xlC-interim-fix-virtualthunk/exe/
   CXXFLAGS=$(INC) $(COPT) $(HAVE_BGPIBM)
   CPP=$(CC) -E

   LD=$(CC)
   LDFLAGS= -O3 -qarch=450d -qtune=450 -qthreaded

endif


OBJECTS = bgp_cc.o  bgp_cr.o  bgp_rc.o bgp_rr.o

all: libBGPmTxm.a

refresh: realclean all


libBGPmTxm.a: $(OBJECTS)
	rm -rf libBGPmTxm.a
	$(AR) $(ARFLAGS) libBGPmTxm.a $(OBJECTS) 
	ranlib libBGPmTxm.a

$%.o: %.cc
	$(CXX) -c $(CXXFLAGS) $(CPPFLAGS) $(INC) $< -o $@

clean:
	$(RM) $(RMFLAGS) *.o

realclean: clean
	$(RM) $(RMFLAGS) *.a



