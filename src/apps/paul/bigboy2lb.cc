#define WORLD_INSTANTIATE_STATIC_TEMPLATES  
#include <mra/mra.h>
#include <mra/mra.h>
#include <mra/operator.h>
#include <constants.h>
//#include <mra/lbdeux.h>
#include <mra/loadbal.h>

using namespace madness;
using namespace std;

const double PI = 3.1415926535897932384;

const int LBMODE=0;
const int K_VAL  = 15;
      int nfunc  = 25;

template <typename T, int NDIM>
struct lbcost {
    double leaf_value;
    double parent_value;
    lbcost(double leaf_value=1.0, double parent_value=0.0) : leaf_value(leaf_value), parent_value(parent_value) {}
    double operator()(const Key<NDIM>& key, const FunctionNode<T,NDIM>& node) const {
        if (node.is_leaf()) {
            return leaf_value;
        }
        else {
            return parent_value;
        }
    }
};

/// Makes a square-normalized Gaussian with random origin and exponent

template <typename T>
T complexify(T c) {return c;}

template <> double_complex complexify<double_complex>(double_complex c) {
    return c*double_complex(0.5,-sqrt(3.0)*0.5);
}

template <> float_complex complexify<float_complex>(float_complex c) {
    return c*float_complex(0.5,-sqrt(3.0)*0.5);
}


template <typename T, int NDIM>
class Gaussian : public FunctionFunctorInterface<T,NDIM> {
public:
    typedef Vector<double,NDIM> coordT;
    const coordT center;
    const double exponent;
    const T coefficient;
    
    Gaussian(const coordT& center, double exponent, T coefficient) 
        : center(center), exponent(exponent), coefficient(complexify(coefficient)) {};

    T operator()(const coordT& x) const {
        double sum = 0.0;
        for (int i=0; i<NDIM; i++) {
            double xx = center[i]-x[i];
            sum += xx*xx;
        };
        return coefficient*exp(-exponent*sum);
    };
};

/// Makes a square-normalized Gaussian with random origin and exponent
template <typename T, int NDIM>
Gaussian<T,NDIM>*
RandomGaussian(const Tensor<double> cell, double expntmax=1e5) {
    typedef Vector<double,NDIM> coordT;
    coordT origin;
    for (int i=0; i<NDIM; i++) {
        origin[i] = RandomValue<double>()*(cell(i,1)-cell(i,0)) + cell(i,0);
    }
//     double lo = log(0.1);
//     double hi = log(expntmax);
//     double expnt = exp(RandomValue<double>()*(hi-lo) + lo);
    double expnt = 30.0;
    T coeff = pow(2.0*expnt/PI,0.25*NDIM);            
    return new Gaussian<T,NDIM>(origin,expnt,coeff);
}

double ttt, sss;
#define START_TIMER world.gop.fence(); ttt=wall_time(); sss=cpu_time()
#define END_TIMER(msg) ttt=wall_time()-ttt; sss=cpu_time()-sss; if (world.rank()==0) printf("timer: %20.20s %8.2fs %8.2fs\n", msg, sss, ttt)
        

typedef Singleton< LBTimer<3> > CostFun;
typedef SharedPtr< FunctionFunctorInterface<double,3> > functorT;
typedef FunctionFactory<double,3> factoryT;
typedef Function<double,3> functionT;

int main(int argc, char** argv) {
    initialize(argc, argv);
    try {
        World world(MPI::COMM_WORLD);
        
        cout.precision(8);
        startup(world,argc,argv);
      
        if (LBMODE == 2) CostFun::Instance(world).set_master_status(true);
 
        Tensor<double> cell(3,2);
        for (int i=0; i<3; i++) {
            cell(i,0) = -30;
            cell(i,1) =  30;
        }
        FunctionDefaults<3>::set_cell(cell);
        FunctionDefaults<3>::set_k(K_VAL);
        FunctionDefaults<3>::set_thresh(1e-8);
        FunctionDefaults<3>::set_refine(true);
        FunctionDefaults<3>::set_autorefine(false);
        FunctionDefaults<3>::set_initial_level(5);
        FunctionDefaults<3>::set_apply_randomize(false);
        FunctionDefaults<3>::set_project_randomize(true);

        default_random_generator.setstate(314159);  // Ensure all processes have the same sequence (for exponents)
    
        if (world.rank() == 0) print("FIRST WITHOUT LOAD BAL", nfunc);
        std::vector< Function<double,3> > v(nfunc);
        world.gop.fence();
        START_TIMER;
        for (int i=0; i<nfunc; i++) {
            v[i] = functionT(factoryT(world).functor(functorT(RandomGaussian<double,3>(cell,100.0))).fence(false));
        }
        world.gop.fence();
        END_TIMER("project");

        world.gop.fence();
        START_TIMER;
        //LoadBalanceDeux<3> lb(world);
        if (LBMODE > 0) {
          CostFun::Instance(world).compute_min_cost();
          LoadBalImpl<3> lb(v[0], CostFun::Instance(world));
          for (int i=1; i<nfunc; i++) {
              lb.add_tree(v[i], CostFun::Instance(world));
          }
          world.gop.fence();
          FunctionDefaults<3>::set_pmap(lb.load_balance());
          CostFun::Instance(world).remap(FunctionDefaults<3>::get_pmap(),true,true);
          CostFun::Instance(world).reset();
          world.gop.fence();
          for (int i=0; i<nfunc; i++) {
            //if (world.rank() == 0) printf("COPY %d\n",i);
            //v[i] = copy(v[i], FunctionDefaults<3>::get_pmap(), false);
          }
          world.gop.fence();
        }
        END_TIMER("load balance");

        default_random_generator.setstate(314159);  // Ensure all processes have the same sequence (for exponents)

        if (world.rank() == 0) print("NOW WITH LOAD BAL");
        world.gop.fence();
        START_TIMER;
        for (int i=0; i<nfunc; i++) {
            v[i] = functionT(factoryT(world).functor(functorT(RandomGaussian<double,3>(cell,100.0))).fence(false));
        }
        world.gop.fence();
        END_TIMER("project");
        
        world.gop.fence();
        START_TIMER;
        compress(world, v);
        END_TIMER("compress");

        world.gop.fence();
        START_TIMER;
        truncate(world, v);
        END_TIMER("truncate");

        START_TIMER;
        reconstruct(world,v);
        END_TIMER("reconstruct");

        SeparatedConvolution<double,3> op = CoulombOperator<double>(world, 
                                                                    FunctionDefaults<3>::get_k(), 
                                                                    1e-3, 
                                                                    FunctionDefaults<3>::get_thresh());

        world.gop.fence();
        // I don't understand why, but I had to put this in to stop
        //   getting a "double-free" exception in the next set of apply ops
        for (int i=0; i<nfunc; i++) {
          v[i] = apply(op, v[i]);
        }
        world.gop.fence();
        START_TIMER;
        apply(world, op, v);
        END_TIMER("apply-1");

        world.gop.fence();
        START_TIMER;
        apply(world, op, v);
        END_TIMER("apply-2");

        world.gop.fence();
        START_TIMER;
        //LoadBalanceDeux<3> lbX(world);
        if (LBMODE > 0) {
          CostFun::Instance(world).compute_min_cost();
          LoadBalImpl<3> lbX(v[0], CostFun::Instance(world));
          for (int i=1; i<nfunc; i++) {
             lbX.add_tree(v[i], CostFun::Instance(world));
          }
          world.gop.fence();
          FunctionDefaults<3>::set_pmap(lbX.load_balance());
          CostFun::Instance(world).remap(FunctionDefaults<3>::get_pmap(),true,true);
          CostFun::Instance(world).reset();
          world.gop.fence();
          for (int i=0; i<nfunc; i++) {
            v[i] = copy(v[i], FunctionDefaults<3>::get_pmap(), false);
          }
          world.gop.fence();
        }
        END_TIMER("load balance truncated");

        START_TIMER;
        reconstruct(world,v);
        END_TIMER("reconstruct");

        world.gop.fence();
        START_TIMER;
        apply(world, op, v);
        END_TIMER("apply-1");

        world.gop.fence();
        START_TIMER;
        apply(world, op, v);
        END_TIMER("apply-2");

        world.gop.fence();
        START_TIMER;
        size_t ncoeff = 0;
        for (int i=0; i<nfunc; i++) ncoeff += v[i].size();
        END_TIMER("count coeff");
        if (world.rank() == 0) print("NCOEFF ", ncoeff);

        ThreadPool::end();
        print_stats(world);

    } catch (const MPI::Exception& e) {
        //        print(e);
        error("caught an MPI exception");
    } catch (const madness::MadnessException& e) {
        print(e);
        error("caught a MADNESS exception");
    } catch (const madness::TensorException& e) {
        print(e);
        error("caught a Tensor exception");
    } catch (const char* s) {
        print(s);
        error("caught a c-string exception");
    } catch (char* s) {
        print(s);
        error("caught a c-string exception");
    } catch (const std::string& s) {
        print(s);
        error("caught a string (class) exception");
    } catch (const std::exception& e) {
        print(e.what());
        error("caught an STL exception");
    } catch (...) {
        error("caught unhandled exception");
    }

    finalize();
    return 0;
}
