
# Blue Gene/P
#BLAS=/soft/apps/LIBBLAS/libgoto_ppc440FP2-r1.02a.VM.a
#HPM=/soft/apps/UPC/lib/libhpm.a
#CXX=mpixlcxx_r
#CXXFLAGS=-g -O5 -qarch=450d -qtune=450 -qthreaded

# LLVM on a laptop
BLAS=-lblas
HPM=
CXX=clang++
CXXFLAGS=-g -O3

all: test_mtxmq.x test_Zmtxmq.x

test_mtxmq.x: test_mtxmq.o mtxmq.o mTxm_tune.o
	$(CXX) $(CXXFLAGS) test_mtxmq.o mtxmq.o mTxm_tune.o $(BLAS) $(HPM) -o test_mtxmq.x

test_Zmtxmq.x: test_Zmtxmq.o mtxmq.o mTxm_tune.o
	$(CXX) $(CXXFLAGS) test_Zmtxmq.o mtxmq.o mTxm_tune.o $(BLAS) $(HPM) -o test_Zmtxmq.x

$%.o: %.cc
	$(CXX) -c $(CXXFLAGS) $< -o $@

clean:
	$(RM) $(RMFLAGS) *.o

realclean: clean
	$(RM) $(RMFLAGS) *.x


