/// \file tdse.cc
/// \brief Evolves the hydrogen molecular ion in 4D ... 3 electron + 1 nuclear degree of freedom


#define WORLD_INSTANTIATE_STATIC_TEMPLATES  
#include <mra/mra.h>
#include <mra/qmprop.h>
#include <mra/operator.h>
#include <constants.h>
#include <tensor/vmath.h>

#include <mra/lbdeux.h>

using namespace madness;

template <typename T, int NDIM>
struct lbcost {
    double leaf_value;
    double parent_value;
    lbcost(double leaf_value=1.0, double parent_value=1.0) : leaf_value(leaf_value), parent_value(parent_value) {}
    double operator()(const Key<NDIM>& key, const FunctionNode<T,NDIM>& node) const {
        if (node.is_leaf()) {
            return leaf_value;
        }
        else {
            return parent_value;
        }
        //return key.level()+1.0;
    }
};


// typedefs to make life less verbose
typedef Vector<double,4> coordT;
typedef SharedPtr< FunctionFunctorInterface<double,4> > functorT;
typedef Function<double,4> functionT;
typedef FunctionFactory<double,4> factoryT;
typedef SeparatedConvolution<double,4> operatorT;
typedef SharedPtr< FunctionFunctorInterface<double_complex,4> > complex_functorT;
typedef Function<double_complex,4> complex_functionT;
typedef FunctionFactory<double_complex,4> complex_factoryT;
typedef SeparatedConvolution<double_complex,4> complex_operatorT;
typedef SharedPtr< WorldDCPmapInterface< Key<4> > > pmapT;

double real(double a) {return a;}

static const double reduced_mass = 0.5*constants::proton_electron_mass_ratio;
static const double sqrtmu = sqrt(reduced_mass);
static const double R0 = 2.04; // Effective center of nuclear wave function
static const double s0 = sqrtmu*R0;
static const double Z=1.0;

struct InputParameters {
  static const int MAXNATOM=99;

    // IF YOU ADD A NEW PARAMETER DON'T FORGET TO INCLUDE IT IN
    // a) read()
    // b) serialize()
    // c) operator<<()
  
  double L;           // Box size for the simulation
  double F;           // Laser field strength
  double omega;       // Laser frequency
  double ncycle;      // Number of laser cycles in envelope
  int k;              // wavelet order
  double thresh;      // precision for truncating wave function
  double safety;      // additional precision (thresh*safety) for operators and potential
  double cut;         // smoothing parameter for 1/r (same for all atoms for now)
  string prefix;      // Prefix for filenames
  int ndump;          // dump wave function to disk every ndump steps
  int nprint;         // print stats every nprint steps
  int nloadbal;       // load balance every nloadbal steps
  int nio;            // Number of IO nodes 
    
  double tScale;      // Scaling parameter for optimization

  double target_time;// Target end-time for the simulation
  
  void read(const char* filename) {
    ifstream f(filename);
    string tag;
    printf("\n");
    printf("       Simulation parameters\n");
    printf("       ---------------------\n");
    printf("             Z = %.1f\n", Z);
    printf("            R0 = %.6f\n", R0);
    printf("            mu = %.6f\n", reduced_mass);    
    printf("      sqrt(mu) = %.6f\n", sqrtmu);    
    while(f >> tag) {
        if (tag[0] == '#') {
            char ch;
            printf("    comment  %s ",tag.c_str());
            while (f.get(ch)) {
                printf("%c",ch);
                if (ch == '\n') break;
            }
        }
        else if (tag == "L") {
            f >> L;
            printf("             L = %.1f\n", L);
        }
        else if (tag == "F") {
            f >> F;
            printf("             F = %.6f\n", F);
        }
        else if (tag == "omega") {
            f >> omega;
            printf("         omega = %.6f\n", omega);
        }
        else if (tag == "ncycle") {
            f >> ncycle;
            printf("         ncycle = %.6f\n", ncycle);
        }
        else if (tag == "k") {
            f >> k;
            printf("             k = %d\n", k);
        }
        else if (tag == "thresh") {
            f >> thresh;
            printf("        thresh = %.1e\n", thresh);
        }
        else if (tag == "safety") {
            f >> safety;
            printf("        safety = %.1e\n", safety);
        }
        else if (tag == "cut") {
            f >> cut;
            printf("           cut = %.2f\n", cut);
        }
        else if (tag == "prefix") {
            f >> prefix;
            printf("        prefix = %s\n", prefix.c_str());
        }
        else if (tag == "ndump") {
            f >> ndump;
            printf("         ndump = %d\n", ndump);
        }
        else if (tag == "nprint") {
            f >> nprint;
            printf("         nprint = %d\n", nprint);
        }
        else if (tag == "nloadbal") {
            f >> nloadbal;
            printf("         nloadbal = %d\n", nloadbal);
        }
        else if (tag == "nio") {
            f >> nio;
            printf("           nio = %d\n", nio);
        }
        else if (tag == "target_time") {
            f >> target_time;
            printf("   target_time = %.3f\n", target_time);
        }
        else if (tag == "tScale") {
            f >> tScale;
            printf("           tScale = %.5f\n", tScale);
        }
        else {
            MADNESS_EXCEPTION("unknown input option", 0);
        }
    }
  }
    
  template <typename Archive>
  void serialize(Archive & ar) {
      ar & L & F & omega & ncycle;
      ar & k & thresh & safety & cut & prefix & ndump & nprint & nloadbal & nio & target_time & tScale;
  }
};

ostream& operator<<(ostream& s, const InputParameters& p) {
    s << p.L<< " " << p.F<< " " << p.omega << " " << 
        p.ncycle << " " << p.k<< " " <<
        p.thresh<< " " << p.cut<< " " << p.prefix<< " " << p.ndump<< " " <<
        p.nprint << " "  << p.nloadbal << " " << p.nio << p.tScale << std::endl;
return s;
}

InputParameters param;

static double zero_field_time;      // Laser actually switches on after this time (set by propagate)
                                    // Delay provides for several steps with no field before start

// This controls the distribution of data across the machine
class LevelPmap : public WorldDCPmapInterface< Key<4> > {
private:
    const int nproc;
public:
    LevelPmap() : nproc(0) {};
    
    LevelPmap(World& world) : nproc(world.nproc()) {}
    
    // Find the owner of a given key
    ProcessID owner(const Key<4>& key) const {
        Level n = key.level();
        if (n == 0) return 0;
        hashT hash;

        // This randomly hashes levels 0-2 and then
        // hashes nodes by their grand-parent key so as
        // to increase locality separately on each level.
        //if (n <= 2) hash = key.hash();
        //else hash = key.parent(2).hash();

        // This randomly hashes levels 0-3 and then 
        // maps nodes on even levels to the same
        // random node as their parent.
        // if (n <= 3 || (n&0x1)) hash = key.hash();
        // else hash = key.parent().hash();

        // This randomly hashes each key
        hash = key.hash();

        return hash%nproc;
    }
};

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

// Potential - nuclear-nuclear repulsion
static double Vn(const coordT& r) {
    double s=r[3];

    s = s + s0;
    if (s < 0.0) s = 0.0; // Do something vaguely sensible for non-physical bond length

    double cut = 10.0;
    return sqrtmu*smoothed_potential(s/cut)/cut;
}

// Potential - electron-nuclear attraction
static double Ve(const coordT& r) {
    const double x=r[0], y=r[1], z=r[2], s=r[3];
    double R = R0 + s/sqrtmu;
    if (R < 0.0) R = 0.0;

    double zz = z-R*0.5;
    double rr = sqrt(x*x+y*y+zz*zz);
    double Va = -Z*smoothed_potential(rr/param.cut)/param.cut;

    zz = z+R*0.5;
    rr = sqrt(x*x+y*y+zz*zz);
    double Vb = -Z*smoothed_potential(rr/param.cut)/param.cut;

    return Va + Vb;
}

// Initial guess wave function using symmetric superposition of 1s orbital on atoms and harmonic oscillator
static double guess(const coordT& r) {
    const double x=r[0], y=r[1], z=r[2], s=r[3];
    const double R = R0 + s/sqrtmu;

    // These from fitting to Predrag's exact function form for psinuc in BO approx
    const double a = 4.42162;
    const double alpha = 1.28164;
    const double beta = -1.06379;
    const double Rp = 2.12902;
    const double Rmax = Rp + sqrt(23.0/a);
    
    // Screen on size of nuclear wave function
    if (R > Rmax) return 0.0;

    // Note in electronic part we are using R0 not R ... OR SHOULD THIS BE R? TRY BOTH!
    static const double face = sqrt(3.14*Z*Z*Z);
    double zz = z-R0*0.5;
    double rr = sqrt(x*x+y*y+zz*zz);
    double psia = face*exp(-Z*rr);
    zz = z+R0*0.5;
    rr = sqrt(x*x+y*y+zz*zz);
    double psib = face*exp(-Z*rr);

    // Nuclear part
    double R2 = (R-Rp)*(R-Rp);
    double psinuc = alpha*exp(-alpha*R2) + beta*(R-Rp)*exp(-1.5*alpha*R2);

    return (psia + psib)*psinuc;
}

// x-dipole electronic
double xdipole(const coordT& r) {
    return r[0];
}

// y-dipole electronic
double ydipole(const coordT& r) {
    return r[1];
}

// z-dipole electronic
double zdipole(const coordT& r) {
    return r[2];
}

double bond_length(const coordT& r) {
    return R0 + r[3]/sqrtmu;
}

// Strength of the laser field at time t
double laser(double t) {
    double omegat = param.omega*t;

    if (omegat < 0.0 || omegat/(2*param.ncycle) > constants::pi) return 0.0;

    double envelope = sin(omegat/(2*param.ncycle));
    envelope *= envelope;
    return param.F*envelope*sin(omegat);
}

double myreal(double t) {return t;}

double myreal(const double_complex& t) {return real(t);}

// Given psi and V evaluate the energy ... leaves psi compressed, potn reconstructed
template <typename T>
double energy(World& world, const Function<T,4>& psi, const functionT& pote, const functionT& potn, const functionT& potf) {
    // First do all work in the scaling function basis
    bool DOFENCE = false;
    psi.reconstruct();
    Function<T,4> dx = diff(psi,0,DOFENCE);
    Function<T,4> dy = diff(psi,1,DOFENCE);
    Function<T,4> dz = diff(psi,2,DOFENCE);
    Function<T,4> ds = diff(psi,3,DOFENCE);
    Function<T,4> Vepsi = psi*pote;
    Function<T,4> Vnpsi = psi*potn;
    Function<T,4> Vfpsi = psi*potf;

    // Now do all work in the wavelet basis
    psi.compress(DOFENCE); 
    Vepsi.compress(DOFENCE); 
    Vnpsi.compress(DOFENCE); 
    Vfpsi.compress(DOFENCE); 
    dx.compress(DOFENCE); 
    dy.compress(DOFENCE); 
    dz.compress(DOFENCE); 
    ds.compress(true);
    double S = real(psi.inner(psi));
    double PEe = real(psi.inner(Vepsi));
    double PEn = real(psi.inner(Vnpsi));
    double PEf = real(psi.inner(Vfpsi));
    double KEe = real(0.5*(inner(dx,dx) + inner(dy,dy) + inner(dz,dz)));
    double KEn = real(0.5*inner(ds,ds));
    double E = (KEe + KEn + PEe + PEn + PEf)/S;

    dx.clear(); dy.clear(); dz.clear(); ds.clear(); Vepsi.clear(); Vepsi.clear(); Vfpsi.clear(); // To free memory on return
    world.gop.fence();
    if (world.rank() == 0) {
        printf("energy=%.6f overlap=%.6f KEe=%.6f KEn=%.6f PEe=%.6f PEn=%.6f PEf=%.6f\n", E, S, KEe, KEn, PEe, PEn, PEf);
     }

    return myreal(E);
}

double fred(const coordT& r) {
    static const double a = 0.1;
    double rsq = r[0]*r[0]+r[1]*r[1]+r[2]*r[2]+r[3]*r[3];
    return 10.0*exp(-a*rsq);
}

double delsqfred(const coordT& r) {
    static const double a = 0.1;
    double rsq = r[0]*r[0]+r[1]*r[1]+r[2]*r[2]+r[3]*r[3];
    return (4.0*a*a*rsq - 4.0*2.0*a)*10.0*exp(-a*rsq);
}
    

void testbsh(World& world) {
    double mu = 0.3;
    functionT test = factoryT(world).f(fred); test.truncate();
    functionT rhs = (mu*mu)*test - functionT(factoryT(world).f(delsqfred));
    rhs.truncate();
    operatorT op = BSHOperator<double,4>(world, mu, FunctionDefaults<4>::get_k(), 1e-2, FunctionDefaults<4>::get_thresh());

    functionT result = apply(op,rhs);

    double err = (result - test).norm2();
    print("ERR", err);
    result.reconstruct();
    for (int i=-100; i<=100; i++) {
        coordT r(i*0.01*5.0);
        double ok = fred(r);
        double bad = result(r);
        print(r[0], ok, bad, ok-bad, ok?bad/ok:0.0);
    }
}

void converge(World& world, functionT& potn, functionT& pote, functionT& pot, functionT& psi, double& eps) {
    functionT zero = factoryT(world);
    for (int iter=0; iter<5; iter++) {
        if (world.rank() == 0) print("beginning iter", iter, wall_time());

        functionT Vpsi = pot*psi;// - 0.5*psi; // TRY SHIFTING POTENTIAL AND ENERGY DOWN
        //eps -= 0.5;

        operatorT op = BSHOperator<double,4>(world, sqrt(-2*eps), param.k, param.cut, param.thresh);
        if (world.rank() == 0) print("made V*psi", wall_time());
        Vpsi.scale(-20.0).truncate();
        if (world.rank() == 0) print("tryuncated V*psi", wall_time());
        functionT tmp = apply(op,Vpsi).truncate(param.thresh);
        if (world.rank() == 0) print("applied operator", wall_time());
        tmp.scale(0.1);
        double norm = tmp.norm2();
        functionT r = tmp-psi;
        double rnorm = r.norm2();
        double eps_new = eps - 0.05*inner(Vpsi,r)/(norm*norm);
        if (world.rank() == 0) {
            print("norm=",norm," eps=",eps," err(psi)=",rnorm," err(eps)=",eps_new-eps);
        }

        // update with damping
        //double d = 0.25;
        //psi = (psi.scale(d) + tmp.scale((1.0-d)/norm));
        psi = tmp;

        psi.scale(1.0/psi.norm2());

        //eps = eps_new;
        eps = energy(world, psi, pote, potn, zero);

        if (rnorm < std::max(1e-5,param.thresh)) break;
    }
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

    pmapT oldpmap = FunctionDefaults<4>::get_pmap();
    LoadBalanceDeux<4> lb(world);
    lb.add_tree(psi1, lbcost<double_complex,4>(1.0,1.0));
    FunctionDefaults<4>::set_pmap(lb.load_balance());
    psi1 = copy(psi1, FunctionDefaults<4>::get_pmap(), true);

    if (world.rank() == 0) print("APPLYING G again", size);
    psi1 = apply(G,psi1);  psi1.truncate(param.thresh);  size = psi1.size();
    if (world.rank() == 0) print("DONE", size);

    FunctionDefaults<4>::set_pmap(oldpmap);
    psi1 = copy(psi1, oldpmap, true);
    
    return psi1;
}

template<typename T, int NDIM>
struct unaryexp {};


template<int NDIM>
struct unaryexp<double_complex,NDIM> {
    void operator()(const Key<NDIM>& key, Tensor<double_complex>& t) const {
        //vzExp(t.size, t.ptr(), t.ptr());
        UNARY_OPTIMIZED_ITERATOR(double_complex, t, *_p0 = exp(*_p0););
    }
    template <typename Archive>
    void serialize(Archive& ar) {}
};


// Returns exp(-I*t*V)
complex_functionT make_exp(double t, const functionT& v) {
    v.reconstruct();
    complex_functionT expV = double_complex(0.0,-t)*v;
    expV.unaryop(unaryexp<double_complex,4>());
    return expV;
}

void print_stats_header(World& world) {
    if (world.rank() == 0) {
        printf("  step       time            field           energy            norm           overlap0         x-dipole         y-dipole         z-dipole           <R>         wall-time(s)\n");
        printf("------- ---------------- ---------------- ---------------- ---------------- ---------------- ---------------- ---------------- ---------------- ---------------- ------------\n");
    }
}

void print_stats(World& world, int step, double t, const functionT& pote,  const functionT& potn, const functionT& potf,
                 const functionT& x, const functionT& y, const functionT& z, const functionT& R,
                 const complex_functionT& psi0, const complex_functionT& psi) {
    double norm = psi.norm2();
    double current_energy = energy(world, psi, pote, potn, potf);
    double xdip = real(inner(psi, x*psi))/(norm*norm);
    double ydip = real(inner(psi, y*psi))/(norm*norm);
    double zdip = real(inner(psi, z*psi))/(norm*norm);
    double Ravg = real(inner(psi, R*psi))/(norm*norm);
    double overlap0 = std::abs(psi.inner(psi0))/norm;
    if (world.rank() == 0) {
        printf("%7d %16.8e %16.8e %16.8e %16.8e %16.8e %16.8e %16.8e %16.8e %16.8e %9.1f\n", step, t, laser(t), current_energy, norm, overlap0, xdip, ydip, zdip, Ravg, wall_time());
    }
}

const char* wave_function_filename(int step) {
    static char fname[1024];
    sprintf(fname, "%s-%5.5d", param.prefix.c_str(), step);
    return fname;
}

const char* wave_function_small_plot_filename(int step) {
    static char fname[1024];
    sprintf(fname, "%s-%5.5dS.dx", param.prefix.c_str(), step);
    return fname;
}

const char* wave_function_large_plot_filename(int step) {
    static char fname[1024];
    sprintf(fname, "%s-%5.5dL.dx", param.prefix.c_str(), step);
    return fname;
}

complex_functionT wave_function_load(World& world, int step) {
    complex_functionT psi;
    ParallelInputArchive ar(world, wave_function_filename(step));
    ar & psi;
    return psi;
}

void wave_function_store(World& world, int step, const complex_functionT& psi) {
    ParallelOutputArchive ar(world, wave_function_filename(step), param.nio);
    ar & psi;
}

bool wave_function_exists(World& world, int step) {
    return ParallelInputArchive::exists(world, wave_function_filename(step));
}


void loadbal(World& world, 
             functionT& pote, functionT& potn, functionT& pot, functionT& vt,
             complex_functionT& psi, complex_functionT& psi0,
             functionT& x, functionT& y, functionT& z, functionT& R) {
    if (world.size() < 2) return;
    if (world.rank() == 0) print("starting LB");
    LoadBalanceDeux<4> lb(world);
    lb.add_tree(vt, lbcost<double,4>(1.0,1.0));
    lb.add_tree(psi, lbcost<double_complex,4>(10.0,10.0));
    FunctionDefaults<4>::set_pmap(lb.load_balance());
    world.gop.fence();
    if (world.rank() == 0) print("starting LB copies");
    world.gop.fence();
    pote = copy(pote, FunctionDefaults<4>::get_pmap(), true);
    potn = copy(potn, FunctionDefaults<4>::get_pmap(), true);
    pot = copy(pot, FunctionDefaults<4>::get_pmap(), true);
    vt = copy(vt, FunctionDefaults<4>::get_pmap(), true);
    psi = copy(psi, FunctionDefaults<4>::get_pmap(), true);
    psi0 = copy(psi0, FunctionDefaults<4>::get_pmap(), true);
    x = copy(x, FunctionDefaults<4>::get_pmap(), true);
    y = copy(y, FunctionDefaults<4>::get_pmap(), true);
    z = copy(z, FunctionDefaults<4>::get_pmap(), true);
    R = copy(R, FunctionDefaults<4>::get_pmap(), true);
    world.gop.fence();
    if (world.rank() == 0) print("done with LB");
    world.gop.fence();
}


// Evolve the wave function in real time starting from given time step on disk
void propagate(World& world, functionT& pote, functionT& potn, functionT& pot, int step0) {
    double ctarget = 5.0/param.cut;
    double c = 1.72*ctarget;   // This for 10^4 steps
    double tcrit = 2*constants::pi/(c*c);

    double time_step = tcrit * param.tScale;

    zero_field_time = 20.0*time_step;

    int nstep = (param.target_time + zero_field_time)/time_step + 1;

    // Ensure everyone has the same data
    world.gop.broadcast(c);
    world.gop.broadcast(time_step);
    world.gop.broadcast(nstep);

    // Free particle propagator for both Trotter and Chin-Chen --- exp(-I*T*time_step/2)
    SeparatedConvolution<double_complex,4> G = qm_free_particle_propagator<4>(world, param.k, c, 0.5*time_step, 2*param.L);

    // Dipole moment functions for laser field and for printing statistics
    functionT x = factoryT(world).f(xdipole);
    functionT y = factoryT(world).f(ydipole);
    functionT z = factoryT(world).f(zdipole);
    functionT R = factoryT(world).f(bond_length);

    // Wave function at time t=0 for printing statistics
    complex_functionT psi0 = wave_function_load(world, 0); 

    int step = step0;  // The current step
    double t = step0 * time_step - zero_field_time;        // The current time
    complex_functionT psi = wave_function_load(world, step); // The wave function at time t
    functionT vt = pot+laser(t)*x; // The total potential at time t

    if (world.rank() == 0) {
        printf("\n");
        printf("        Evolution parameters\n");
        printf("       --------------------\n");
        printf("     bandlimit = %.2f\n", ctarget);
        printf(" eff-bandlimit = %.2f\n", c);
        printf("         tcrit = %.6f\n", tcrit);
        printf("     time step = %.6f\n", time_step);
        printf(" no field time = %.6f\n", zero_field_time);
        printf("   target time = %.2f\n", param.target_time);
        printf("         nstep = %d\n", nstep);
        printf("\n");
        printf("  restart step = %d\n", step0);
        printf("  restart time = %.6f\n", t);
        printf("\n");
    }

    print_stats_header(world);
    print_stats(world, step, t, pote, potn, laser(t)*x, x, y, z, R, psi0, psi0); 
    world.gop.fence();

    psi.truncate();

    while (step < nstep) {
        if (step < 2 || (step%param.nloadbal) == 0)
            loadbal(world, pote, potn, pot, vt, psi, psi0, x, y, z, R);

        long depth = psi.max_depth(); long size=psi.size();
        if (world.rank() == 0) print("depth size", depth, size);

        // Make the potential at time t + step/2 
        functionT vhalf = pot + laser(t+0.5*time_step)*x;

        // Apply Trotter to advance from time t to time t+step
        complex_functionT expV = make_exp(time_step, vhalf);
        psi = trotter(world, expV, G, psi);

        // Update counters, print info, dump/plot as necessary
        step++;
        t += time_step;
        vt = pot+laser(t)*x;

        print_stats(world, step, t, pote, potn, laser(t)*x, x, y, z, R, psi0, psi); 

        if ((step%param.ndump) == 0 || step==nstep) {
            double start = wall_time();
            wave_function_store(world, step, psi);
            // Update the restart file for automatic restarting
            if (world.rank() == 0) {
                std::ofstream("restart4") << step << std::endl;
                print("dumping took", wall_time()-start);
            }
            world.gop.fence();
        }
    }
}

void doit(World& world) {
    cout.precision(8);

    if (world.rank() == 0) param.read("input4");
    world.gop.broadcast_serializable(param, 0);

    FunctionDefaults<4>::set_k(param.k);                        // Wavelet order
    FunctionDefaults<4>::set_thresh(param.thresh*param.safety);       // Accuracy
    FunctionDefaults<4>::set_initial_level(4);
    FunctionDefaults<4>::set_cubic_cell(-param.L,param.L);
    FunctionDefaults<4>::set_apply_randomize(true);
    FunctionDefaults<4>::set_autorefine(false);
    FunctionDefaults<4>::set_truncate_mode(1);
    FunctionDefaults<4>::set_pmap(pmapT(new LevelPmap(world)));

    // Read restart information
    int step0;               // Initial time step ... filenames are <prefix>-<step0>
    if (world.rank() == 0) std::ifstream("restart4") >> step0;
    world.gop.broadcast(step0);

    bool exists = wave_function_exists(world, step0);

    if (world.rank() == 0) {
        print("EXISTS",exists,"STEP0",step0);
        ofstream out("plot.dat");
        for (int i=-100; i<=100; i++) {
            double x = i*0.01*param.L;
            coordT rn(0.0), re(0.0);
            rn[3]=x;
            re[2]=x;
            double vn = Vn(rn);
            double ve = Ve(re);
            double pn = guess(rn);
            double pe = guess(re);
            out << x << " " << vn << " " << ve << " " << pn << " " << pe << endl;
        }
        out.close();
    }

    // Make the potential
    if (world.rank() == 0) print("COMPRESSING Vn",wall_time());
    functionT potn = factoryT(world).f(Vn);  potn.truncate();
    if (world.rank() == 0) print("COMPRESSING Ve",wall_time());
    functionT pote = factoryT(world).f(Ve);  pote.truncate();
    functionT pot = potn + pote;

    LoadBalanceDeux<4> lb(world);
    lb.add_tree(pot, lbcost<double,4>());
    FunctionDefaults<4>::set_pmap(lb.load_balance());
    world.gop.fence();
    pote = copy(pote, FunctionDefaults<4>::get_pmap());
    potn = copy(potn, FunctionDefaults<4>::get_pmap());
    pot = copy(pot, FunctionDefaults<4>::get_pmap());


    if (!exists) {
        if (step0 == 0) {
            if (world.rank() == 0) print("Computing initial ground state wavefunction", wall_time());
            functionT psi = factoryT(world).f(guess);
            double norm0 = psi.norm2();
            psi.scale(1.0/norm0);
            psi.truncate();
            if (world.rank() == 0) print("computed norm", norm0, "at", wall_time());
            
            double eps = energy(world, psi, pote, potn, functionT(factoryT(world)));
            if (world.rank() == 0) print("guess energy", eps, wall_time()); 
            converge(world, potn, pote, pot, psi, eps);
            
            psi.truncate(param.thresh);

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

    propagate(world, pote, potn, pot, step0);
}

int main(int argc, char** argv) {
    initialize(argc,argv);
    World world(MPI::COMM_WORLD);
    
    startup(world,argc,argv);

    try {
        doit(world);
    } catch (const MPI::Exception& e) {
        //print(e); std::cout.flush();
        error("caught an MPI exception");
    } catch (const madness::MadnessException& e) {
        print(e); std::cout.flush();
        error("caught a MADNESS exception");
    } catch (const madness::TensorException& e) {
        print(e); std::cout.flush();
        error("caught a Tensor exception");
    } catch (const char* s) {
        print(s); std::cout.flush();
        error("caught a c-string exception");
    } catch (char* s) {
        print(s); std::cout.flush();
        error("caught a c-string exception");
    } catch (const std::string& s) {
        print(s); std::cout.flush();
        error("caught a string (class) exception");
    } catch (const std::exception& e) {
        print(e.what()); std::cout.flush();
        error("caught an STL exception");
    } catch (...) {
        error("caught unhandled exception");
    }


    world.gop.fence();
    
    ThreadPool::end();
    print_stats(world);
    finalize();
    return 0;
}


