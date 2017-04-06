MY_ROOT=.
DEBUG=YES
#HPM_PROFILING=YES

ifndef TARGET
$(error You must set TARGET to BGP or BGQ)
endif

#############################################

ifeq ($(TARGET),BGP)


   RM = rm
   RMFLAGS = -f

   AR = powerpc-bgp-linux-ar
   ARFLAGS = cru
   RANLIB = powerpc-bgp-linux-ranlib
   HAVE_BGPIBM = -DHAVE_IBMBGP

   LIB=-L$(MY_ROOT)
   INC=-I$(MY_ROOT)

   CC=mpixlc_r
   CXX=mpixlcxx_r
   COPT= -g -O5 -qarch=450d -qtune=450 -qthreaded #-qpath=ILbc:/soft/apps/xlC-interim-fix-virtualthunk/exe/
   CXXFLAGS=$(INC) $(COPT) $(HAVE_BGPIBM)
   CPP=$(CC) -E

   LD=$(CC)
   LDFLAGS= -O3 -qarch=450d -qtune=450 -qthreaded

   OBJECTS = bgp_cc.o  bgp_cr.o  bgp_rc.o bgp_rr.o
endif


ifeq ($(TARGET),BGQ)
   RM = rm
   RMFLAGS = -f

   AR = powerpc64-bgq-linux-ar
   ARFLAGS = cru
   RANLIB = powerpc64-bgq-linux-ranlib
   HAVE_BGQIBM = -DHAVE_IBMBGQ

   LIB=-L$(MY_ROOT)
   INC=-I$(MY_ROOT) -I/bgsys/drivers/ppcfloor/

   #CC=/home/projects/llvm/bin/bgclang
   #CXX=/home/projects/llvm/bin/bgclang++
   #COPT= -O3
   #LDFLAGS= -O3

   CC=bgxlc_r
   CXX=bgxlc++_r
   COPT= -g -O5 -qhot=level=1 -qsimd=auto -qtune=qp -qarch=qp
   CXXFLAGS=$(INC) $(COPT) $(HAVE_BGQIBM)
   LDFLAGS= -O3 -qarch=qp -qtune=qp -qthreaded
   CPP=$(CC) -E

   LD=$(CC)

   OBJECTS = bgq_mtxm.o
endif

all: libMADMTXM.a

refresh: realclean all

libMADMTXM.a: $(OBJECTS)
	rm -rf libMADMTXM.a
	$(AR) $(ARFLAGS) libMADMTXM.a $(OBJECTS)
	$(RANLIB) libMADMTXM.a

$%.o: %.cc
	$(CXX) -c $(CXXFLAGS) $(CPPFLAGS) $(INC) $< -o $@

clean:
	$(RM) $(RMFLAGS) *.o

realclean: clean
	$(RM) $(RMFLAGS) *.a



