/*!
  \file h2.cc
  \brief Solves the Hartree-Fock equations for the hydrogen molecule
  \defgroup examplesh2hf Hatree-Fock equations for the hydrogen molecule
  \ingroup examples

  The Hartree-Fock wave function is computed for the hydrogen molecule
  in three dimensions without using symmetry.

  Since all of the details except for the nuclear potential are the
  same, please refer to the \ref examplehehf helium atom HF example.

*/


#define WORLD_INSTANTIATE_STATIC_TEMPLATES  
#include <mra/mra.h>
#include <mra/operator.h>

using namespace madness;

static const double R = 1.4;    // bond length
static const double L = 64.0*R; // box size
static const long k = 6;        // wavelet order
static const double thresh = 1e-4; // precision

static double guess(const coord_3d& r) {
    const double x=r[0], y=r[1], z=r[2];
    return (exp(-sqrt(x*x+y*y+(z-R/2)*(z-R/2)+1e-8))+
            exp(-sqrt(x*x+y*y+(z+R/2)*(z+R/2)+1e-8)));
}

static double V(const coord_3d& r) {
    const double x=r[0], y=r[1], z=r[2];
    return -1.0/sqrt(x*x+y*y+(z-R/2)*(z-R/2)+1e-8)+
           -1.0/sqrt(x*x+y*y+(z+R/2)*(z+R/2)+1e-8);
}

void iterate(World& world, real_function_3d& V, real_function_3d& psi, double& eps) {
    real_convolution_3d op = BSHOperator3D<double>(world, sqrt(-2*eps), k, 0.001, 1e-6);
    real_function_3d Vpsi = (V*psi);
    Vpsi.scale(-2.0).truncate();
    real_function_3d tmp = apply(op,Vpsi).truncate();
    double norm = tmp.norm2();
    real_function_3d r = tmp-psi;
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
    World world(MPI::COMM_WORLD);
    
    startup(world,argc,argv);
    std::cout.precision(6);

    FunctionDefaults<3>::set_k(k);
    FunctionDefaults<3>::set_thresh(thresh);
    FunctionDefaults<3>::set_refine(true);
    FunctionDefaults<3>::set_initial_level(5);
    FunctionDefaults<3>::set_truncate_mode(1);  
    FunctionDefaults<3>::set_cubic_cell(-L/2, L/2);
    // for (int i=0; i<3; i++) {
    // FunctionDefaults<3>::cell(i,0) = -L/2;
    // FunctionDefaults<3>::cell(i,1) =  L/2;
    // }
    
    real_function_3d Vnuc = real_factory_3d(world).f(V);
    real_function_3d psi  = real_factory_3d(world).f(guess);
    psi.truncate();
    psi.scale(1.0/psi.norm2());

    real_convolution_3d op = CoulombOperator<double>(world, k, 0.001, 1e-6);

    double eps = -0.6;
    for (int iter=0; iter<10; iter++) {
        real_function_3d rho = square(psi).truncate();
        real_function_3d potential = Vnuc + apply(op,rho).truncate();
        iterate(world, potential, psi, eps);
    }

    double kinetic_energy = 0.0;
    for (int axis=0; axis<3; axis++) {
        print("DOING AXIS",axis);
        real_function_3d dpsi = diff(psi,axis);
        kinetic_energy += inner(dpsi,dpsi);
    }

    real_function_3d rho = square(psi);
    double two_electron_energy = inner(apply(op,rho),rho);
    double nuclear_attraction_energy = 2.0*inner(Vnuc,rho);
    double nuclear_repulsion_energy = 1.0/R;
    double total_energy = kinetic_energy + two_electron_energy + 
        nuclear_attraction_energy + nuclear_repulsion_energy;
    double virial = (nuclear_attraction_energy + two_electron_energy + nuclear_repulsion_energy) / kinetic_energy;

    if (world.rank() == 0) {
        print("            Kinetic energy ", kinetic_energy);
        print(" Nuclear attraction energy ", nuclear_attraction_energy);
        print("       Two-electron energy ", two_electron_energy);
        print(" Nuclear  repulsion energy ", nuclear_repulsion_energy);
        print("              Total energy ", total_energy);
        print("                    Virial ", virial);
    }

    world.gop.fence();

    finalize();
    return 0;
}
