/// \file he.cc
/// \brief Solves the Hatree-Fock equations for the helium atom


#define WORLD_INSTANTIATE_STATIC_TEMPLATES  
#include <mra/mra.h>
#include <mra/operator.h>

using namespace madness;

typedef Vector<double,3> coordT;
typedef SharedPtr< FunctionFunctorInterface<double,3> > functorT;
typedef Function<double,3> functionT;
typedef FunctionFactory<double,3> factoryT;
typedef SeparatedConvolution<double,3> operatorT;

static const double L = 32.0;   // box size
static const long k = 6;        // wavelet order
static const double thresh = 1e-6; // precision
static const double thresh1 = thresh;

static double guess(const coordT& r) {
    const double x=r[0], y=r[1], z=r[2];
    return 6.0*exp(-2.0*sqrt(x*x+y*y+z*z+1e-4));
}

static double V(const coordT& r) {
    const double x=r[0], y=r[1], z=r[2];
    return -2.0/(sqrt(x*x+y*y+z*z+1e-8));
}

void iterate(World& world, functionT& V, functionT& psi, double& eps) {
    print("Multiplying V*psi");
    functionT Vpsi = (V*psi);
    Vpsi.scale(-2.0).truncate(thresh1);
    print("Applying BSH operator");
    operatorT op = BSHOperator<double,3>(world, sqrt(-2*eps), k, 0.001, 1e-6);
    functionT tmp = apply(op,Vpsi).truncate(thresh1);
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

int main(int argc, char** argv) {
    MPI::Init(argc, argv);
    ThreadPool::begin();
    RMI::begin();
    MPI::COMM_WORLD.Barrier();
    World world(MPI::COMM_WORLD);
    
    startup(world,argc,argv);
    cout.precision(6);

    FunctionDefaults<3>::set_k(k);
    FunctionDefaults<3>::set_thresh(thresh);
    FunctionDefaults<3>::set_truncate_mode(1);  
    FunctionDefaults<3>::set_cubic_cell(-L/2,L/2);
    
    print("Projecting V");
    functionT Vnuc = factoryT(world).f(V).thresh(thresh1).truncate_mode(0);
    print("Projecting psi");
    functionT psi  = factoryT(world).f(guess);
    psi.scale(1.0/psi.norm2());
    operatorT op = CoulombOperator<double,3>(world, k, 0.001, 1e-6);

    double eps = -1.0; 
    for (int iter=0; iter<20; iter++) {
        print("Computing density");
        functionT rho = square(psi).truncate(thresh1);
        print("Computing Coulomb potential");
        functionT potential = Vnuc + apply(op,rho).truncate(thresh1);
        iterate(world, potential, psi, eps);
    }

    double kinetic_energy = 0.0;
    for (int axis=0; axis<3; axis++) {
        print("DOING AXIS",axis);
        functionT dpsi = diff(psi,axis);
        kinetic_energy += inner(dpsi,dpsi);
    }

    functionT rho = square(psi).truncate(thresh1);
    double two_electron_energy = inner(apply(op,rho),rho);
    double nuclear_attraction_energy = 2.0*inner(Vnuc*psi,psi);
    double total_energy = kinetic_energy + nuclear_attraction_energy + two_electron_energy;

    print("psinorm", psi.norm2());
    coordT r(0.0);
    psi.reconstruct();
    for (int i=0; i<201; i++) {
        r[2] = -L/2 + L*i/200.0;
        print(r[2], psi(r));
    }

    if (world.rank() == 0) {
        print("            Kinetic energy ", kinetic_energy);
        print(" Nuclear attraction energy ", nuclear_attraction_energy);
        print("       Two-electron energy ", two_electron_energy);
        print("              Total energy ", total_energy);
        print("                    Virial ", (nuclear_attraction_energy + two_electron_energy) / kinetic_energy);
    }

    world.gop.fence();

    RMI::end();
    MPI::Finalize();
    return 0;
}
