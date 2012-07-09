MY_ROOT=.
DEBUG=YES
#HPM_PROFILING=YES

TARGET=BGQ


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

   OBJECTS = bgp_cc.o  bgp_cr.o  bgp_rc.o bgp_rr.o
endif


ifeq ($(TARGET),BGQ)
   RM = rm
   RMFLAGS = -f

   AR = ar
   ARFLAGS = cru
   HAVE_BGQIBM = -DHAVE_IBMBGQ

   LIB=-L$(MY_ROOT)
   INC=-I$(MY_ROOT) -I/bgsys/drivers/ppcfloor/

   CC=xlc_r
   CXX=xlc++_r
   COPT= -O5 -qhot=level=1 -qsimd=auto -qtune=qp -qarch=qp
   CXXFLAGS=$(INC) $(COPT) $(HAVE_BGPIBM)
   CPP=$(CC) -E

   LD=$(CC)
   LDFLAGS= -O3 -qarch=qp -qtune=qp -qthreaded

   OBJECTS = bgq_rr.o
endif

all: libMADMTXM.a

refresh: realclean all


libMADMTXM.a: $(OBJECTS)
	rm -rf libMADMTXM.a
	$(AR) $(ARFLAGS) libMADMTXM.a $(OBJECTS)
	ranlib libMADMTXM.a

$%.o: %.cc
	$(CXX) -c $(CXXFLAGS) $(CPPFLAGS) $(INC) $< -o $@

clean:
	$(RM) $(RMFLAGS) *.o

realclean: clean
	$(RM) $(RMFLAGS) *.a



