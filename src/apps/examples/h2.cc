/// \file h2.cc
/// \brief Solves the Hatree-Fock equations for the hydrogen molecule


#define WORLD_INSTANTIATE_STATIC_TEMPLATES  
#include <mra/mra.h>
#include <mra/operator.h>

using namespace madness;

typedef Vector<double,3> coordT;
typedef SharedPtr< FunctionFunctorInterface<double,3> > functorT;
typedef Function<double,3> functionT;
typedef FunctionFactory<double,3> factoryT;
typedef SeparatedConvolution<double,3> operatorT;

static const double R = 1.4;    // bond length
static const double L = 32.0*R; // box size
static const long k = 5;        // wavelet order
static const double thresh = 1e-3; // precision
static const double thresh1 = thresh*0.1;

static double guess(const coordT& r) {
    const double x=r[0], y=r[1], z=r[2];
    return (exp(-sqrt(x*x+y*y+(z-R/2)*(z-R/2)))+
            exp(-sqrt(x*x+y*y+(z+R/2)*(z+R/2))));
}

static double V(const coordT& r) {
    const double x=r[0], y=r[1], z=r[2];
    return -1.0/sqrt(x*x+y*y+(z-R/2)*(z-R/2))+
           -1.0/sqrt(x*x+y*y+(z+R/2)*(z+R/2));
}

void iterate(World& world, functionT& V, functionT& psi, double& eps) {
    operatorT op = BSHOperator3D<double>(world, sqrt(-2*eps), k, 0.001, 1e-6);
    functionT Vpsi = (V*psi);
    Vpsi.scale(-2.0).truncate(thresh1);
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
    initialize(argc, argv);
    ThreadPool::begin();
    RMI::begin();
    MPI::COMM_WORLD.Barrier();
    World world(MPI::COMM_WORLD);
    
    startup(world,argc,argv);
    cout.precision(6);

    FunctionDefaults<3>::set_k(k);
    FunctionDefaults<3>::set_thresh(thresh);
    FunctionDefaults<3>::set_refine(true);
    FunctionDefaults<3>::set_initial_level(2);
    FunctionDefaults<3>::set_truncate_mode(0);  
    FunctionDefaults<3>::set_cubic_cell(-L/2, L/2);
    // for (int i=0; i<3; i++) {
    // FunctionDefaults<3>::cell(i,0) = -L/2;
    // FunctionDefaults<3>::cell(i,1) =  L/2;
    // }
    
    functionT Vnuc = factoryT(world).f(V);
    functionT psi  = factoryT(world).f(guess);
    psi.truncate(thresh1);
    psi.scale(1.0/psi.norm2());

    operatorT op = CoulombOperator<double>(world, k, 0.001, 1e-6);

    double eps = -0.6;
    for (int iter=0; iter<10; iter++) {
        functionT rho = square(psi).truncate(thresh1);
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
    double nuclear_repulsion_energy = 1.0/R;
    double total_energy = kinetic_energy + two_electron_energy + 
        nuclear_attraction_energy + nuclear_repulsion_energy;

    if (world.rank() == 0) {
        print("            Kinetic energy ", kinetic_energy);
        print(" Nuclear attraction energy ", nuclear_attraction_energy);
        print("       Two-electron energy ", two_electron_energy);
        print(" Nuclear  repulsion energy ", nuclear_repulsion_energy);
        print("              Total energy ", total_energy);
        print("                    Virial ", (nuclear_attraction_energy + two_electron_energy) / kinetic_energy);
    }

    world.gop.fence();

    finalize();
    return 0;
}
