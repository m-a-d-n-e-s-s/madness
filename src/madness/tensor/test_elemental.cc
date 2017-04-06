#include <madness/madness_config.h>


#define WORLD_INSTANTIATE_STATIC_TEMPLATES
#include <madness/world/MADworld.h>
using namespace madness;
using namespace std;

#ifdef MADNESS_HAS_ELEMENTAL

#include <madness/tensor/elem.h>
using namespace elem;

#ifdef STATIC
#  undef STATIC
#endif

#if HAVE_UNQUALIFIED_STATIC_DECL
#  define STATIC static
#else
// Cray X1 compiler won't instantiate static function templates (mxm*)
#  define STATIC
#endif


namespace madness {
    // This stupidity since current defn of conj_tranpose() only
    // takes complex types ... was that a sensible design choice?
    // STATIC Tensor<float> my_conj_transpose(Tensor<float> a) {
    //     return transpose(a);
    // }
    STATIC Tensor<double> my_conj_transpose(Tensor<double> a) {
        return transpose(a);
    }
    // STATIC Tensor<float_complex> my_conj_transpose(Tensor<float_complex> a) {
    //     return conj_transpose(a);
    // }
    // STATIC Tensor<double_complex> my_conj_transpose(Tensor<double_complex> a) {
    //     return conj_transpose(a);
    // }
}

double xtt1, xss1;

#define START_TIMER world.gop.fence(); xtt1=wall_time(); xss1=cpu_time()
#define END_TIMER(msg) xtt1=wall_time()-xtt1; xss1=cpu_time()-xss1; if (world.rank()==0) printf("timer: %20.20s %8.2fs %8.2fs\n", msg, xss1, xtt1)

template <typename T>
double test_sygvp(World& world, int n) {
    Tensor<T> a(n,n), V, b(n,n);
    Tensor< typename Tensor<T>::scalar_type > e;
    
    a.fillrandom(); a-=0.5;
    b.fillrandom(); b-=0.5;
    a += madness::my_conj_transpose(a); a*=0.5;
    b += madness::my_conj_transpose(b); b*=0.5;
    for (int i=0; i<n; ++i) b(i,i) = 2*n; // To make pos-def

    // Tensor<T> aa = copy(a);
    // Tensor<T> bb = copy(b);
    // Tensor<T> VV;
    // Tensor< typename Tensor<T>::scalar_type > ee;
    
    // world.gop.broadcast(a.ptr(),a.size(), 0);
    // world.gop.broadcast(b.ptr(),b.size(), 0);
    // world.gop.fence();
    
    START_TIMER;
    sygvp(world, a,b,1,V,e);
    END_TIMER("call sygvp");

    // print(a);
    // print(b);
    // print(V);
    // print(e);

    // sygv(aa, bb, 1, VV, ee);
    // print(VV);
    // print(ee);
    // {
    //     double err = 0.0;
    //     for (int i=0; i<n; ++i) {
    //         err = max(err,(double) (inner(a,VV(_,i)) - inner(b,VV(_,i))*(T) ee(i)).normf());
    //     }
    //     print(err);
    // }
    
    double err = 0.0;
    for (int i=0; i<n; ++i) {
        err = max(err,(double) (inner(a,V(_,i)) - inner(b,V(_,i))*(T) e(i)).normf());
    }
    return err;
}

template <typename T>
double test_gesvp(World& world, int n, int nrhs) {
    Tensor<T> a(n,n), b1(n), b(n,nrhs), x1, x;
    
    a.fillrandom();
    b1.fillrandom();
    b.fillrandom();
    
    // world.gop.broadcast(a.ptr(),a.size(), 0);
    // world.gop.broadcast(b.ptr(),b.size(), 0);
    // world.gop.broadcast(b1.ptr(),b1.size(), 0);
    // world.gop.fence();
    
    //         print("A");
    //         print(a);
    
    //         print("B");
    //         print(b);
    
    //         print("B1");
    //         print(b1);
    
    START_TIMER;
    gesvp(world, a,b,x);
    gesvp(world, a,b1,x1);
    END_TIMER(" call gesvp");
    
    //         print("X");
    //         print(x);
    
    //         print("X1");
    //         print(x1);
    
    //         print("R");
    //         print(inner(a,x)-b);
    //
    //         print("R1");
    //         print(inner(a,x1)-b1);
    
    double err = 0;
    err = (inner(a,x)-b).normf() + (inner(a,x1)-b1).normf();
    return err;
    //        return 111.0;
}

double ij(int64_t i, int64_t j) {return ((i<<24) | j);}

template <typename T>
void test_copy(World& world, int n, int m) {
    const int blocksize = 63;
    const elem::Grid GG( world.mpi.comm().Get_mpi_comm() );
    elem::SetBlocksize(blocksize);
    elem::DistMatrix<T> gd( n, m, GG );
    DistributedMatrix<T> A = column_distributed_matrix<T>(world, n, m, 17);
    A.fill(ij);
    
    copy_to_elemental(A, gd);

    // {
    //     const int64_t colShift =    gd.ColShift(); // first row we own
    //     const int64_t rowShift =    gd.RowShift(); // first col we own
    //     const int64_t colStride =   gd.ColStride();
    //     const int64_t rowStride =   gd.RowStride();
    //     const int64_t localHeight = gd.LocalHeight();
    //     const int64_t localWidth =  gd.LocalWidth();
        
    //     for( int64_t jLocal=0; jLocal<localWidth; ++jLocal ) {
    //         for( int64_t iLocal=0; iLocal<localHeight; ++iLocal ) {
    //             const int64_t i = colShift + iLocal*colStride;
    //             const int64_t j = rowShift + jLocal*rowStride;

    //             print(i,j,gd.Get(i,j),A.get(i,j),A.owner(i,j),ij(i,j));
    //             //const ProcessID owner = dout.owner(i,j);
    //             //s.insert(owner, detail::Value<T>(i,j,gd.GetLocal(iLocal,jLocal)));
    //         }
    //     }
    // }

    A.fill(T(0.0));
    copy_from_elemental(gd, A);
    int64_t ilo, ihi, jlo, jhi;
    A.local_colrange(ilo,ihi);
    A.local_rowrange(jlo,jhi);
    const Tensor<T>& t = A.data();
    for (int64_t i=ilo; i<=ihi; i++) {
        for (int64_t j=jlo; j<=jhi; j++) {
            //print(i,j,A.get(i,j),ij(i,j));
            MADNESS_ASSERT(A.get(i,j) == ij(i,j));
        }
    }
}

double afiller(int64_t i, int64_t j) {
    return std::cos((i+j)/10.0) + std::sin(i*j/(i+j+1.0));
}

template <typename T>
void test_distributed_eval(World& world, const int N) {
    DistributedMatrix<T> A = column_distributed_matrix<T>(world, N, N);
    DistributedMatrix<T> B = row_distributed_matrix<T>(world, N, N); // for perversity
    B.fill_identity();
    A.fill(afiller);

    DistributedMatrix<T> X;
    Tensor<typename Tensor<T>::scalar_type> e;
    
    sygv(A, B, 1, X, e);

    Tensor<T> AA(N,N); A.copy_to_replicated(AA);
    Tensor<T> BB(N,N);  B.copy_to_replicated(BB);
    Tensor<T> XX(N,N); X.copy_to_replicated(XX);

    double err = 0.0;
    for (int i=0; i<N; ++i) {
        err = max(err,(double) (inner(AA,XX(_,i)) - inner(BB,XX(_,i))*(T) e(i)).normf());
    }
    if (world.rank() == 0) print(N,err);
}

int main(int argc, char** argv) {
    initialize(argc, argv);
    World world(SafeMPI::COMM_WORLD);

    int myrank = world.rank();
    
    try {
        for (int n=1; n<100; n++) {
            double err = test_sygvp<double>(world, n);
            if (myrank == 0)  cout << "n=" << n << " error in double sygvp " << err << endl;
        }
        for (int n=1; n<=128; n*=2) {// was 1024
            double err = test_sygvp<double>(world, n);
            if (myrank == 0)  cout << "n=" << n << " error in double sygvp " << err << endl;
        }
        if (myrank == 0)  cout << endl; 
        
        double err = test_gesvp<double>(world, 1800, 1200);
        if (myrank == 0)  cout << "error in float gesvp " << err << endl;
        if (myrank == 0)  cout << endl; 

        test_copy<double>(world,300,300);

        for (int n=4; n<=512; n*=2)
            test_distributed_eval<double>(world,n);
    }
    catch (SafeMPI::Exception e) {
        error("caught an MPI exception");
    }
    catch (madness::MadnessException e) {
        print(e);
        error("caught a MADNESS exception");
    }
    catch (const char* s) {
        print(s);
        error("caught a string exception");
    }
    catch (...) {
        error("caught unhandled exception");
    }

    finalize();
    return 0;
}

#else

int main(int argc, char** argv) {
    print("MADNESS was not configured with elemental");
    return 0;
}


#endif
