#include <madness_config.h>

#ifdef MADNESS_HAS_ELEMENTAL

#define WORLD_INSTANTIATE_STATIC_TEMPLATES
#include <world/world.h>
#include <linalg/elem.h>

using namespace elem;
using namespace madness;
using namespace std;

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

double tt1, ss1;

#define START_TIMER world.gop.fence(); tt1=wall_time(); ss1=cpu_time()
#define END_TIMER(msg) tt1=wall_time()-tt1; ss1=cpu_time()-ss1; if (world.rank()==0) printf("timer: %20.20s %8.2fs %8.2fs\n", msg, ss1, tt1)

template <typename T>
double test_sygvp(World& world, int n) {
    Tensor<T> a(n,n), V, b(n,n);
    Tensor< typename Tensor<T>::scalar_type > e;
    
    a.fillrandom();
    b.fillrandom();
    a += madness::my_conj_transpose(a);
    b += madness::my_conj_transpose(b);
    for (int i=0; i<n; ++i) b(i,i) = 2*n;	// To make pos-def
    
    // world.gop.broadcast(a.ptr(),a.size(), 0);
    // world.gop.broadcast(b.ptr(),b.size(), 0);
    // world.gop.fence();
    
    START_TIMER;
    sygvp(world, a,b,1,V,e);
    END_TIMER("call sygvp");
    
    START_TIMER;
    double err = 0.0;
    for (int i=0; i<n; ++i) {
        err = max(err,(double) (inner(a,V(_,i)) - inner(b,V(_,i))*(T) e(i)).normf());
    }
    END_TIMER("error sygvp");
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
    START_TIMER;
    err = (inner(a,x)-b).normf() + (inner(a,x1)-b1).normf();
    END_TIMER("error gesvp");
    return err;
    //        return 111.0;
}

int main(int argc, char** argv) {
    initialize(argc, argv);
    World world(SafeMPI::COMM_WORLD);

    int myrank = world.rank();
    
    try {
        double err  = 999.99;
        
        err = test_sygvp<double>(world, 100);
        if (myrank == 0)  cout << "error in double sygvp " << err << endl;
        if (myrank == 0)  cout << endl; 
        
        err = test_gesvp<double>(world, 180,120);
        if (myrank == 0)  cout << "error in float gesvp " << err << endl;
        if (myrank == 0)  cout << endl; 
        
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


#endif
