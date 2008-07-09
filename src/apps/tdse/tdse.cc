/// \file tdse.cc
/// \brief Evolves the hydrogen atom in imaginary and also real time


#define WORLD_INSTANTIATE_STATIC_TEMPLATES  
#include <mra/mra.h>
#include <mra/qmprop.h>
#include <mra/operator.h>
#include <constants.h>

using namespace madness;

// Convenient but sleazy use of global variables to define simulation parameters
static const double L = 400.0;      // Box size for the simulation
static const double Lsmall = 20.0;  // Box size for small (near nucleus) plots
static const double Llarge = 200.0; // Box size for large (far from nucleus) plots

static const double F = 0.125;     // Laser field strength
static const double omega = 0.057; // Laser frequency
static const double Z = 1.0;       // Nuclear charge

static const long k = 12;          // wavelet order
static const double thresh = 1e-7; // precision
static const double cut = 0.2;     // smoothing parameter for 1/r

const string prefix = "tdse";      // Prefix for filenames
const int ndump = 10;              // dump wave function to disk every ndump steps
const int nplot = 1;               // dump opendx plot to disk every nplot steps
const int nio = 1;                 // Number of IO nodes 

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

// This controls the distribution of data across the machine
class LevelPmap : public WorldDCPmapInterface< Key<3> > {
private:
    const int nproc;
public:
    LevelPmap() : nproc(0) {};
    
    LevelPmap(World& world) : nproc(world.nproc()) {}
    
    // Find the owner of a given key
    ProcessID owner(const Key<3>& key) const {
        Level n = key.level();
        if (n == 0) return 0;
        hashT hash;
//         if (n <= 3 || (n&0x1)) hash = key.hash();
//         else hash = key.parent().hash();
        hash = key.hash();
        return hash%nproc;
    }
};

// Derivative of the smoothed 1/r approximation

// Invoke as \c du(r/c)/(c*c) where \c c is the radius of the smoothed volume.  
static double d_smoothed_potential(double r) {
    double r2 = r*r;

    if (r > 6.5) {
        return -1.0/r2;
    }
    else if (r > 1e-2) {
        return -(1.1283791670955126*(0.88622692545275800*erf(r)-exp(-r2)*r*(1.0-r2)))/r2;
    }
    else {
        return (-1.880631945159187623160265+(1.579730833933717603454623-0.7253866074185437975046736*r2)*r2)*r;
    }
}

// Smoothed 1/r potential.

// Invoke as \c u(r/c)/c where \c c is the radius of the smoothed volume.  
static double smoothed_potential(double r) {
    double r2 = r*r, pot;
    if (r > 6.5){
        pot = 1.0/r;
    } else if (r > 1e-2) {
        pot = erf(r)/r + exp(-r2)*0.56418958354775630;
    } else{
        pot = 1.6925687506432689-r2*(0.94031597257959381-r2*(0.39493270848342941-0.12089776790309064*r2));
    }
    
    return pot;
}

// Nuclear attraction potential
static double V(const coordT& r) {
    const double x=r[0], y=r[1], z=r[2];
    const double rr = sqrt(x*x+y*y+z*z);
    return -Z*smoothed_potential(rr/cut)/cut;
}

// dV/dx
static double dVdx(const coordT& r) {
    const double x=r[0], y=r[1], z=r[2];
    if (x == 0.0) return 0.0;
    const double rr = sqrt(x*x+y*y+z*z);
    return -Z*x*d_smoothed_potential(rr/cut)/(rr*cut*cut);
}

// dV/dy
static double dVdy(const coordT& r) {
    const double x=r[0], y=r[1], z=r[2];
    if (y == 0.0) return 0.0;
    const double rr = sqrt(x*x+y*y+z*z);
    return -Z*y*d_smoothed_potential(rr/cut)/(rr*cut*cut);
}

// dV/dz
static double dVdz(const coordT& r) {
    const double x=r[0], y=r[1], z=r[2];
    if (z == 0.0) return 0.0;
    const double rr = sqrt(x*x+y*y+z*z);
    return -Z*z*d_smoothed_potential(rr/cut)/(rr*cut*cut);
}

// Initial guess wave function for 1e atoms
static double guess(const coordT& r) {
    const double x=r[0], y=r[1], z=r[2];
    return exp(-1.0*sqrt(x*x+y*y+z*z+cut*cut)); // Change 1.0 to 0.6 to make bad guess
}

// x-dipole
double xdipole(const coordT& r) {
    return r[0];
}

// y-dipole
double ydipole(const coordT& r) {
    return r[1];
}

// z-dipole
double zdipole(const coordT& r) {
    return r[2];
}

// Strength of the laser field at time t ... 1 full cycle
double laser(double t) {
    double omegat = omega*t;

    if (omegat < 0.0 || omegat/24.0 > constants::pi) return 0.0;

    double envelope = sin(omegat/24.0);
    envelope *= envelope;
    return F*envelope*sin(omegat);
}

double myreal(double t) {return t;}

double myreal(const double_complex& t) {return real(t);}

// Given psi and V evaluate the energy
template <typename T>
double energy(World& world, const Function<T,3>& psi, const functionT& potn) {
    PROFILE_FUNC;
    T S = psi.inner(psi);
    T PE = psi.inner(psi*potn);
    T KE = 0.0;
    for (int axis=0; axis<3; axis++) {
        Function<T,3> dpsi = diff(psi,axis);
        KE += inner(dpsi,dpsi)*0.5;
    }
    T E = (KE+PE)/S;
    world.gop.fence();
//     if (world.rank() == 0) {
//         print("the overlap integral is",S);
//         print("the kinetic energy integral",KE);
//         print("the potential energy integral",PE);
//         print("the total energy",E);
//     }

    return myreal(E);
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


// Returns exp(-I*t*V)
complex_functionT make_exp(double t, const functionT& v) {
    PROFILE_FUNC;
    v.reconstruct();
    complex_functionT expV = double_complex(0.0,-t)*v;
    expV.unaryop(unaryexp<double_complex,3>());
    return expV;
}

void print_stats_header(World& world) {
    printf("step time  field  energy  norm  overlap0 x-dipole  y-dipole z-dipole accel");
}

void print_stats(World& world, int step, double t, const functionT& v, 
                 const functionT& x, const functionT& y, const functionT& z,
                 const complex_functionT& psi0, const complex_functionT& psi) {
    PROFILE_FUNC;
    double current_energy = energy(world, psi, v);
    double xdip = real(inner(psi, x*psi));
    double ydip = real(inner(psi, y*psi));
    double zdip = real(inner(psi, z*psi));
    double norm = psi.norm2();
    double overlap0 = std::abs(psi.inner(psi0));
    double accel = 0.0;
    if (world.rank() == 0) {
        printf("%6d %12.8f %12.8f %12.8f %12.8f %12.8f %12.8f %12.8f %12.8f %12.8f\n", step, t, laser(t), current_energy, norm, overlap0, xdip, ydip, zdip, accel);
        print("step", step, "time", t, "field", laser(t), "norm", norm, "energy", current_energy, "dipole", xdip, ydip, zdip);
    }
}

const char* wave_function_filename(int step) {
    static char fname[1024];
    sprintf(fname, "%s-%5.5d", prefix.c_str(), step);
    return fname;
}

const char* wave_function_small_plot_filename(int step) {
    static char fname[1024];
    sprintf(fname, "%s-%5.5dS.dx", prefix.c_str(), step);
    return fname;
}

const char* wave_function_large_plot_filename(int step) {
    static char fname[1024];
    sprintf(fname, "%s-%5.5dL.dx", prefix.c_str(), step);
    return fname;
}

complex_functionT wave_function_load(World& world, int step) {
    complex_functionT psi;
    ParallelInputArchive ar(world, wave_function_filename(step));
    ar & psi;
    return psi;
}

void wave_function_store(World& world, int step, const complex_functionT& psi) {
    ParallelOutputArchive ar(world, wave_function_filename(step), nio);
    ar & psi;
}

bool wave_function_exists(World& world, int step) {
    return ParallelInputArchive::exists(world, wave_function_filename(step));
}

void doplot(World& world, int step, const complex_functionT& psi, double Lplot, long numpt, const char* fname) { 
    Tensor<double> cell(3,2);
    std::vector<long> npt(3, numpt);
    cell(_,0) = -Lplot;
    cell(_,1) =  Lplot;
    plotdx(psi, fname, cell, npt);
}

// Evolve the wave function in real time starting from given time step on disk
void propagate(World& world, int step0) {
    PROFILE_FUNC;
    //double ctarget = 10.0/cut;                // From Fourier analysis of the potential
    double ctarget = 1.0/cut;                   // From wishful thinking
    double c = 1.86*ctarget;
    double tcrit = 2*constants::pi/(c*c);

    double time_step = tcrit;
    int nstep = 20;

    // Ensure everyone has the same data
    world.gop.broadcast(c);
    world.gop.broadcast(time_step);

    if (world.rank() == 0) {
        print("bandlimit",ctarget,"effband",c,"tcrit",tcrit,"time_step",time_step,"nstep",nstep);
    }

    // Free particle propagator for both Trotter and Chin-Chen --- exp(-I*T*time_step/2)
    SeparatedConvolution<double_complex,3> G = qm_free_particle_propagator<3>(world, k, c, 0.5*time_step, 2*L);

    // The time-independent part of the potential plus derivatives for
    // Chin-Chen and also for computing the power spectrum ... compute
    // derivatives analytically to reduce numerical noise
    functionT potn = factoryT(world).f(V);         potn.truncate();
    functionT dpotn_dx = factoryT(world).f(dVdx);  dpotn_dx.truncate();
    functionT dpotn_dy = factoryT(world).f(dVdy);  dpotn_dy.truncate();
    functionT dpotn_dz = factoryT(world).f(dVdz);  dpotn_dz.truncate();

    // Dipole moment functions for laser field and for printing statistics
    functionT x = factoryT(world).f(xdipole);
    functionT y = factoryT(world).f(ydipole);
    functionT z = factoryT(world).f(zdipole);

    // Wave function at time t=0 for printing statistics
    complex_functionT psi0 = wave_function_load(world, 0); 

    int step = step0;                    // The current step
    double t = step0 * time_step;        // The current time
    complex_functionT psi = wave_function_load(world, step); // The wave function at time t

    functionT vt = potn+laser(t)*z; // The total potential at time t

    print_stats(world, step0, t, vt, x, y, z, psi0, psi);

    while (step < nstep) {
        {
            // Make the potential at time t + step/2 
            functionT vhalf = potn + laser(t+0.5*time_step)*z;

            // Apply Trotter to advance from time t to time t+step
            complex_functionT expV = make_exp(time_step, vhalf);
            psi = trotter(world, expV, G, psi);
        }

        // Update counters, print info, dump/plot as necessary
        step++;
        t += time_step;
        vt = potn+laser(t)*z;

        print_stats(world, step, t, vt, x, y, z, psi0, psi);

        if ((step%ndump) == 0 || step==nstep) {
            wave_function_store(world, step, psi);
            // Update the input file for automatic restarting
            if (world.rank() == 0) std::ofstream("input") << step;
            world.gop.fence();
        }

        if ((step%nplot) == 0 || step==nstep) {
            doplot(world, step, psi,  20.0, 200, wave_function_small_plot_filename(step));
            doplot(world, step, psi, 200.0, 200, wave_function_large_plot_filename(step));
        }
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
    FunctionDefaults<3>::set_apply_randomize(false);
    FunctionDefaults<3>::set_autorefine(false);
    FunctionDefaults<3>::set_truncate_mode(1);
    FunctionDefaults<3>::set_pmap(pmapT(new LevelPmap(world)));

    functionT potn = factoryT(world).f(V);  potn.truncate();
    
    int step0;               // Initial time step ... filenames are <prefix>-<step0>
    
    if (world.rank() == 0) 
        std::ifstream("input") >> step0;
    world.gop.broadcast(step0);

    bool exists = wave_function_exists(world, step0);

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
            wave_function_store(world, 0, psic);
        }
        else {
            if (world.rank() == 0) {
                print("The requested restart was not found ---", step0);
                error("restart failed", 0);
            }
            world.gop.fence();
        }
    }

    if (world.rank() == 0) 
        print("Restarting from time step", step0);

    propagate(world, step0);

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

