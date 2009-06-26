#define WORLD_INSTANTIATE_STATIC_TEMPLATES  
#include <mra/mra.h>
#include <mra/mra.h>
#include <mra/operator.h>
#include <constants.h>
#include <mra/lbdeux.h>

using namespace madness;
using namespace std;


typedef SharedPtr< FunctionFunctorInterface<double,3> > functorT;
typedef FunctionFactory<double,3> factoryT;
typedef Function<double,3> functionT;
typedef Vector<double,3> coordT;

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

const double PI = 3.1415926535897932384;

// list of centers
vector<coordT> centers;

const double expnt = 1.0;
const double coeff = pow(expnt/PI,1.5);

// Test function
class Test : public FunctionFunctorInterface<double,3> {
public:
    double operator()(const coordT& x) const {
        double sum = 0.0;
        for (unsigned int i=0; i<centers.size(); i++) {
            const coordT& y = centers[i];
            double xx = x[0]-y[0];
            double yy = x[1]-y[1];
            double zz = x[2]-y[2];
            
            double arg = expnt*(xx*xx + yy*yy + zz*zz);
            if (arg < 46) sum += coeff*exp(-arg);
        }

        return sum;
    }

//     vector<coordT> special_points() const {
//         return centers;
//     }

//     virtual Level special_level() {
//         return 8;
//     }

};

double ttt_, sss_;
#define START_TIMER world.gop.fence(); ttt_=wall_time(); sss_=cpu_time()
#define END_TIMER(msg) ttt_=wall_time()-ttt_; sss_=cpu_time()-sss_; if (world.rank()==0) printf("timer: %20.20s %8.2fs %8.2fs\n", msg, sss_, ttt_)
        

int main(int argc, char** argv) {
    initialize(argc, argv);
    try {
        World world(MPI::COMM_WORLD);
        
        startup(world,argc,argv);
        cout.precision(8);
        
        Tensor<double> cell(3,2);
        double L = 30.0;
        for (int i=0; i<3; i++) {
            cell(i,0) = -L;
            cell(i,1) =  L;
        }
        FunctionDefaults<3>::set_cell(cell);
        FunctionDefaults<3>::set_k(14);
        FunctionDefaults<3>::set_thresh(1e-12);
        FunctionDefaults<3>::set_refine(true);
        FunctionDefaults<3>::set_autorefine(false);
        FunctionDefaults<3>::set_initial_level(3);
        FunctionDefaults<3>::set_apply_randomize(false);
        FunctionDefaults<3>::set_project_randomize(true);

        const int ncent = 2; //10; // No. of centers is ncent**3
        double h = 2*L/ncent;
        coordT v;
        for (int ix=0; ix<ncent; ix++) {
            double x = -L + h*(ix+0.5);
            for (int iy=0; iy<ncent; iy++) {
                double y = -L + h*(iy+0.5);
                for (int iz=0; iz<ncent; iz++) {
                    double z = -L + h*(iz+0.5);
                    //if (world.rank() == 0) print(ix, iy, iz, x, y, z);
                    v[0] = x; v[1] = y; v[2] = z;
                    centers.push_back(v);
                }
            }
        }

        if (world.rank() == 0) print("FIRST WITHOUT LOAD BAL");
        START_TIMER;
        functionT rho = factoryT(world).functor(functorT(new Test()));
        world.gop.fence();
        END_TIMER("project");

        START_TIMER;
        LoadBalanceDeux<3> lb(world);
        lb.add_tree(rho, lbcost<double,3>(1.0,1.0));
        FunctionDefaults<3>::set_pmap(lb.load_balance(2.0));
        world.gop.fence();
        END_TIMER("load balance");

        if (world.rank() == 0) print("NOW WITH LOAD BAL");
        START_TIMER;
        rho = functionT(factoryT(world).functor(functorT(new Test())));
        world.gop.fence();
        END_TIMER("project");

        START_TIMER;
        double norm = rho.trace();
        END_TIMER("trace in scaling fn");
        if (world.rank() == 0) print("trace", norm);
        
        START_TIMER;
        norm = rho.norm2();
        END_TIMER("norm2 in scaling fn");
        if (world.rank() == 0) print("norm2", norm);
        
        START_TIMER;
        rho.compress();
        world.gop.fence();
        END_TIMER("compress");

        START_TIMER;
        norm = rho.trace();
        END_TIMER("trace in wavelets");
        if (world.rank() == 0) print("trace", norm);
        
        START_TIMER;
        norm = rho.norm2();
        END_TIMER("norm2 in wavelets");
        if (world.rank() == 0) print("norm2", norm);
        
        START_TIMER;
        rho.reconstruct();
        world.gop.fence();
        END_TIMER("reconstruct");

        START_TIMER;
        rho.compress();
        world.gop.fence();
        END_TIMER("compress");

        world.gop.fence();
        START_TIMER;
        rho.truncate();
        END_TIMER("truncate");

        START_TIMER;
        rho.reconstruct();
        END_TIMER("reconstruct");

        START_TIMER;
        rho.compress();
        END_TIMER("compress");

        START_TIMER;
        rho.reconstruct();
        END_TIMER("reconstruct");

        SeparatedConvolution<double,3> op = CoulombOperator<double>(world, 
                                                                    FunctionDefaults<3>::get_k(), 
                                                                    1e-3, 
                                                                    FunctionDefaults<3>::get_thresh());

        world.gop.fence();
        START_TIMER;
        apply(op, rho);
        END_TIMER("apply-1");

        world.gop.fence();
        START_TIMER;
        apply(op, rho);
        END_TIMER("apply-2");

        world.gop.fence();
        START_TIMER;
        LoadBalanceDeux<3> lbX(world);
        lbX.add_tree(rho, lbcost<double,3>(1.0,1.0));
        world.gop.fence();
        FunctionDefaults<3>::set_pmap(lbX.load_balance(2.0));
        world.gop.fence();
        rho = copy(rho, FunctionDefaults<3>::get_pmap());
        world.gop.fence();
        END_TIMER("load balance truncated");

        START_TIMER;
        rho.reconstruct();
        END_TIMER("reconstruct");

        START_TIMER;
        rho.compress();
        END_TIMER("compress");

        START_TIMER;
        rho.reconstruct();
        END_TIMER("reconstruct");

        world.gop.fence();
        START_TIMER;
        apply(op, rho);
        END_TIMER("apply-1");

        world.gop.fence();
        START_TIMER;
        apply(op, rho);
        END_TIMER("apply-2");

        world.gop.fence();
        START_TIMER;
        size_t ncoeff = rho.size();
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
