/// \file tdse.cc
/// \brief Evolves the hydrogen atom in imaginary and also real time


#define WORLD_INSTANTIATE_STATIC_TEMPLATES  
#include <mra/mra.h>
#include <mra/qmprop.h>
#include <mra/operator.h>
#include <constants.h>

using namespace madness;

// Convenient but sleazy use of global variables to define simulation parameters
static const double L = 330.0;
static const long k = 12;          // wavelet order
static const double thresh = 1e-8; // precision
static const double cut = 0.3;     // smoothing parameter for 1/r    0.1 is what we seek
static const double F = 1.0;       // Laser field strength
static const double omega = 1.0;   // Laser frequency
static const double Z = 1.0;       // Nuclear charge

// typedefs to make life less verbose
typedef Vector<double,3> coordT;
typedef SharedPtr< FunctionFunctorInterface<double,3> > functorT;
typedef Function<double,3> functionT;
typedef FunctionFactory<double,3> factoryT;
typedef SeparatedConvolution<double,3> operatorT;
typedef SharedPtr< FunctionFunctorInterface<double_complex,3> > complex_functorT;
typedef Function<double_complex,3> complex_functionT;
typedef FunctionFactory<double_complex,3> complex_factoryT;
typedef SeparatedConvolution<double_complex,3> complex_operatorT;
typedef SharedPtr< WorldDCPmapInterface< Key<3> > > pmapT;


class LevelPmap : public WorldDCPmapInterface< Key<3> > {
private:
    const int nproc;
public:
    LevelPmap() : nproc(0) {};
    
    LevelPmap(World& world) : nproc(world.nproc()) {}
    
    /// Find the owner of a given key
    ProcessID owner(const Key<3>& key) const {
        Level n = key.level();
        if (n == 0) return 0;
        hashT hash;
        if (n <= 3 || (n&0x1)) hash = key.hash();
        else hash = key.parent().hash();
        //hashT hash = key.hash();
        return hash%nproc;
    }
};

// /// Regularized 1/r potential.

// /// Invoke as \c u(r/c)/c where \c c is the radius of the
// /// smoothed volume.
// static double old_smoothed_potential(double r) {
//     const double PI = 3.1415926535897932384;
//     const double THREE_SQRTPI = 5.31736155271654808184;
//     double r2 = r*r, pot;
//     if (r > 6.5){
//         pot = 1.0/r;
//     } else if (r > 1e-8){
//         pot = erf(r)/r + (exp(-r2) + 16.0*exp(-4.0*r2))/(THREE_SQRTPI);
//     } else{
//         pot = (2.0 + 17.0/3.0)/sqrt(PI);
//     }
    
//     return pot;
// }


// /// Derivative of new smoothed 1/r approximation
// static double d_smoothed_potential(double r) {
//     double r2 = r*r;

//     if (r > 6.5) {
//         return -1.0/r2;
//     }
//     else if (r > 1e-2) {
//         return -(1.1283791670955126*(0.88622692545275800*erf(r)-exp(-r2)*r*(1.0-r2)))/r2;
//     }
//     else {
//         return (-1.880631945159187623160265+(1.579730833933717603454623-0.7253866074185437975046736*r2)*r2)*r;
//     }
// }


/// Regularized 1/r potential.

/// Invoke as \c u(r/c)/c where \c c is the radius of the
/// smoothed volume.  Only the zero-moment of the error is non-zero.
static double smoothed_potential(double r) {
    double r2 = r*r, pot;
    if (r > 6.5){
        pot = 1.0/r;
    } else if (r > 1e-2){
        pot = erf(r)/r + exp(-r2)*0.56418958354775630;
    } else{
        pot = 1.6925687506432689-r2*(0.94031597257959381-r2*(0.39493270848342941-0.12089776790309064*r2));
    }
    
    return pot;
}


static inline double s(double x) {
  /* Iterated first beta function to switch smoothly 
     from 0->1 in [0,1].  n iterations produce 2*n-1 
     zero derivatives at the end points. Order of polyn
     is 3^n.

     Currently use one iteration so that first deriv.
     is zero at interior boundary and is exactly representable
     by low order multiwavelet without refinement */
#define B1(x) (x*x*(3.-2.*x))
  x = B1(x);
  return x;
}

double mask_function(const coordT& r) {
    const double lo = 0.0625;
    const double hi = 1.0-lo;
    double result = 1.0;

    coordT rsim;
    user_to_sim(r, rsim);

    for (int d=0; d<3; d++) {
        double x = rsim[d];
        if (x<lo)
            result *= s(x/lo);
        else if (x>hi)
            result *= s((1.0-x)/lo);
    }

    return result;
}

/// Smoothed 1/r nuclear potential
static double V(const coordT& r) {
    const double x=r[0], y=r[1], z=r[2];
    const double rr = sqrt(x*x+y*y+z*z);
    return -Z*smoothed_potential(rr/cut)/cut;
}

/// Initial guess wave function
static double guess(const coordT& r) {
    const double x=r[0], y=r[1], z=r[2];
    return exp(-1.0*sqrt(x*x+y*y+z*z+cut*cut)); // Change 1.0 to 0.6 to make bad guess
}

/// z-dipole
double zdipole(const coordT& r) {
    return r[2];
}

/// Strength of the laser field at time t ... 1 full cycle
double laser(double t) {
    double omegat = omega*t;
    if (omegat < 0.0 || omegat > 1) return 0.0;
    else return F*sin(2*constants::pi*omegat);
}

/// Given psi and V evaluate the energy
template <typename T>
double energy(World& world, const Function<T,3>& psi, const functionT& potn) {
    T S = psi.inner(psi);
    T PE = psi.inner(psi*potn);
    T KE = 0.0;
    for (int axis=0; axis<3; axis++) {
        Function<T,3> dpsi = diff(psi,axis);
        KE += inner(dpsi,dpsi)*0.5;
    }
    T E = (KE+PE)/S;
    world.gop.fence();

    if (world.rank() == 0) {
        print("the overlap integral is",S);
        print("the kinetic energy integral",KE);
        print("the potential energy integral",PE);
        print("the total energy",E);
    }
    return -std::abs(E); // ugh
}

template <typename T, int NDIM>
Cost lbcost(const Key<NDIM>& key, const FunctionNode<T,NDIM>& node) {
  return 1;
}

void converge(World& world, functionT& potn, functionT& psi, double& eps) {
    PROFILE_FUNC;
    for (int iter=0; iter<10; iter++) {
        operatorT op = BSHOperator<double,3>(world, sqrt(-2*eps), k, cut, thresh);
        functionT Vpsi = (potn*psi);
        Vpsi.scale(-2.0).truncate();
        functionT tmp = apply(op,Vpsi).truncate();
        double norm = tmp.norm2();
        functionT r = tmp-psi;
        double rnorm = r.norm2();
        double eps_new = eps - 0.5*inner(Vpsi,r)/(norm*norm);
        if (world.rank() == 0) {
            print("norm=",norm," eps=",eps," err(psi)=",rnorm," err(eps)=",eps_new-eps);
        }
        psi = tmp.scale(1.0/norm);
        eps = eps_new;
    }
}

complex_functionT chin_chen(const complex_functionT& expV,
                            const complex_functionT& expVtilde,
                            const complex_operatorT& G,
                            const complex_functionT& psi0) {

    // psi(t) = exp(-i*V(t)*t/6) exp(-i*T*t/2) exp(-i*2*Vtilde(t/2)*t/3) exp(-i*T*t/2) exp(-i*V(0)*t/6)

    complex_functionT psi1;

    psi1 = expV*psi0;       psi1.truncate();
    psi1 = apply(G,psi1);   psi1.truncate();
    psi1 = expVtilde*psi1;  psi1.truncate();
    psi1 = apply(G,psi1);   psi1.truncate();
    psi1 = expV*psi1;       psi1.truncate();

    return psi1;
}

complex_functionT trotter(World& world,
                          const complex_functionT& expV, 
                          const complex_operatorT& G, 
                          const complex_functionT& psi0) {
    //    psi(t) = exp(-i*T*t/2) exp(-i*V(t/2)*t) exp(-i*T*t/2) psi(0)

    complex_functionT psi1;

    unsigned long size = psi0.size();
    if (world.rank() == 0) print("APPLYING G", size);
    psi1 = apply(G,psi0);  psi1.truncate();  size = psi1.size();
    if (world.rank() == 0) print("APPLYING expV", size);
    psi1 = expV*psi1;      psi1.truncate();  size = psi1.size();
    if (world.rank() == 0) print("APPLYING G again", size);
    psi1 = apply(G,psi1);  psi1.truncate();  size = psi1.size();
    if (world.rank() == 0) print("DONE", size);
    
    return psi1;
}

template<typename T, int NDIM>
struct unaryexp {
    void operator()(const Key<NDIM>& key, Tensor<T>& t) const {
        UNARY_OPTIMIZED_ITERATOR(T, t, *_p0 = exp(*_p0););
    }
    template <typename Archive>
    void serialize(Archive& ar) {}
};


/// Returns exp(-I*t*V)
complex_functionT make_exp(double t, const functionT& v) {
    PROFILE_FUNC;
    v.reconstruct();
    complex_functionT expV = double_complex(0.0,-t)*v;
    expV.unaryop(unaryexp<double_complex,3>());
    return expV;
}

/// Evolve the wave function in real time
void propagate(World& world, const functionT& potn, const complex_functionT& psi0, const double eps) {
    PROFILE_FUNC;
    double ctarget = 10.0/cut;                // From Fourier analysis of the potential
    double c = 1.86*ctarget;
    double tcrit = 2*constants::pi/(c*c);

    double tstep = tcrit; //
    int nstep = 1;
    //double Eshift = eps;

    //potn.add_scalar(-Eshift);

    // Ensure everyone has the same data
    world.gop.broadcast(c);
    world.gop.broadcast(tstep);

    if (world.rank() == 0) {
        print("bandlimit",ctarget,"effband",c,"tcrit",tcrit,"tstep",tstep,"nstep",nstep);
    }

    complex_functionT psi = copy(psi0);
    SeparatedConvolution<double_complex,3> G = qm_free_particle_propagator<3>(world, k, c, 0.5*tstep, 2*L);
    //G.doleaves = true;
    complex_functionT expV = make_exp(tstep, potn);

    for (int step=0; step<nstep; step++) {
        double t = step * tstep;
        double_complex phase = psi0.inner(psi);
        double radius = abs(phase);
        double theta = arg(phase);
        double theta_exact = -t*eps;
        while (theta_exact > constants::pi) theta_exact -= 2.0*constants::pi;
        while (theta_exact < -constants::pi) theta_exact += 2.0*constants::pi;
      
        if (world.rank() == 0) 
            print("step", step, "time", t, "radius", radius, "arg", theta, "exact", theta_exact, "phase err", theta_exact-theta);

        //energy(v, psi);
      
        psi = trotter(world, expV, G, psi);
    }
}

int main(int argc, char** argv) {
    MPI::Init(argc, argv);
    World world(MPI::COMM_WORLD);
    
    startup(world,argc,argv);
    cout.precision(8);

    FunctionDefaults<3>::set_k(k);                 // Wavelet order
    FunctionDefaults<3>::set_thresh(thresh);       // Accuracy
    FunctionDefaults<3>::set_refine(true);         // Enable adaptive refinement
    FunctionDefaults<3>::set_initial_level(4);     // Initial projection level
    FunctionDefaults<3>::set_cubic_cell(-L,L);
    FunctionDefaults<3>::set_apply_randomize(true);
    FunctionDefaults<3>::set_autorefine(false);
    FunctionDefaults<3>::set_truncate_mode(1);
    FunctionDefaults<3>::set_pmap(pmapT(new LevelPmap(world)));

    functionT potn = factoryT(world).f(V);  potn.truncate();
    
    string prefix = "tdse";  // Prefix for filenames
    int ndump = 10;          // Dump restart info every ndump time steps
    int step0;               // Initial time step ... filenames are <prefix>-<step0>
    
    if (world.rank() == 0) {
        std::ifstream f("input");
        f >> step0;
    }
    world.gop.broadcast(step0);

    // See if the restart file exists
    char buf[256];
    sprintf(buf, "%s-%5.5d", prefix.c_str(), step0);
    bool exists = ParallelInputArchive::exists(world, buf);

    if (!exists) {
        if (step0 == 0) {
            if (world.rank() == 0) print("Computing initial ground state wavefunction");
            functionT psi = factoryT(world).f(guess);
            psi.scale(1.0/psi.norm2());
            psi.truncate();
            psi.scale(1.0/psi.norm2());
            
            double eps = energy(world, psi, potn);
            converge(world, potn, psi, eps);

            complex_functionT psic = double_complex(1.0,0.0)*psi;
            
            ParallelOutputArchive ar(world,  buf, 1);
            ar & eps & psic;
        }
        else {
            if (world.rank() == 0) {
                print("The requested restart file was not found", buf);
                error("restart failed", 0);
            }
            world.gop.fence();
        }
    }

    complex_functionT psi;
    double eps;
    ParallelInputArchive ar(world, buf, 1);
    ar & eps & psi;
    ar.close();

    if (world.rank() == 0) 
        print("Restarting from time step", step0);

    propagate(world, potn, psi, eps);

    world.gop.fence();
    if (world.rank() == 0) {
        world.am.print_stats();
        world.taskq.print_stats();
        world_mem_info()->print();
    }

    WorldProfile::print(world);

    MPI::Finalize();
    return 0;
}

