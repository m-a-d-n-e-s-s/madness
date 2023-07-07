#include <madness.h>
#include <madchem.h>

using namespace madness;

int main(int argc, char* argv[]) {
    World& world=initialize(argc,argv);
    startup(world,argc,argv,true);
    FunctionDefaults<3>::set_cubic_cell(-10,10);
    FunctionDefaults<3>::set_k(8);
    FunctionDefaults<3>::set_truncate_mode(1);
    FunctionDefaults<3>::set_thresh(1.e-5);
    double lo=1.e-4;

    try {
        // set up molecule
        Atom O( 0.0, 0.0, 0.2, 8, 8);
        Atom H1(0.0, 1.4,-0.8, 1, 1);
        Atom H2(0.0,-1.4,-0.8, 1, 1);
        Molecule molecule({O,H1,H2},1.e-6);
        molecule.print();
        long nocc=long(molecule.total_nuclear_charge())/2;

        // define Fock operator
        auto T=Kinetic<double,3>(world);
        auto J=Coulomb<double,3>(world,lo,FunctionDefaults<3>::get_thresh());
        auto K=Exchange<double,3>(world,lo,FunctionDefaults<3>::get_thresh());
        auto V=Nuclear<double,3>(world,molecule);

        // compute hcore guess
        AtomicBasisSet aobasis("sto-3g");
        auto aos=SCF::project_ao_basis_only(world,aobasis,molecule);
        auto hcoremat=T(aos,aos) + V(aos,aos);
        auto smat=matrix_inner(world,aos,aos);
        Tensor<double> U,e;
        sygv(hcoremat,smat,1,U,e);              // solves gen. eigenvalue problem: A x =  e B x
        auto orbitals=transform(world,aos,U(_,Slice(0,nocc-1)));   // orbitals are lowest nocc AOs after diagonalization

        // set up KAIN solver
        auto solver=nonlinear_vector_solver<double,3>(world, orbitals.size());

        // solves F |i> = f_{ij} |j>    <=>    (T-f_{ii} ) |i>  = -(J - K + V) |i> + \sum_{j/neq i> f_{ij} |j>
        BSHApply<double,3> bsh_apply(world);

        int maxiter=30;
        for (int iter=0; iter<maxiter; ++iter) {
            double wall0=wall_time();

            // update Fock operator
            auto density=2.0*dot(world,orbitals,orbitals);
            J.potential()=J.compute_potential(density);
            K.set_bra_and_ket(orbitals,orbitals);

            // compute potentials
            auto Jphi=J(orbitals);
            auto Kphi=K(orbitals);
            auto Vnucphi=V(orbitals);
            auto Vphi=truncate(Jphi-Kphi+Vnucphi);

            // compute Fock matrix
            auto tmat=T(orbitals,orbitals);
            auto jmat=matrix_inner(world,orbitals,Jphi);
            auto kmat=matrix_inner(world,orbitals,Kphi);
            auto vmat=matrix_inner(world,orbitals,Vnucphi);
            auto fock=tmat + jmat -kmat + vmat;

            // apply BSH operator
            auto [residual,eps_update]=bsh_apply(orbitals,fock,Vphi);
            double rnorm=norm2(world,residual);

            // update orbitals, orthonormalize
            orbitals=solver.update(orbitals,residual);
//            orbitals-=residual;
            orbitals=truncate(orthonormalize(orbitals));

            // compute energy
            double energy=molecule.nuclear_repulsion_energy();
            for (long i=0; i<nocc; ++i) {
                energy+=2.0* (tmat(i,i) +0.5*jmat(i,i) -0.5* kmat(i,i)+vmat(i,i));
            }
            double wall1=wall_time();
            printf("iteration %2d with total energy: %12.8f %3.1e after %4.1fs\n",iter, energy, rnorm,wall1-wall0);
            if (rnorm<1.e-4) break;

        }

    } catch (...) {
         std::cout << "caught an error " << std::endl;
    }
    finalize();
    return 0;
}
