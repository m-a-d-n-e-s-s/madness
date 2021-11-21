//
// Created by Florian Bischoff on 11/1/21.
//



#include<madness/mra/mra.h>
#include<madness/mra/vmra.h>
#include<chem/localizer.h>
#include<chem/molecularbasis.h>
#include<chem/MolecularOrbitals.h>
#include<chem/potentialmanager.h>
#include<chem/correlationfactor.h>
#include<madness/world/test_utilities.h>
#include<chem/SCFOperators.h>
#include<chem/SCF.h>

#include<vector>

using namespace madness;

/*
 * test different localization methods including core-valence separation
 */


/// compute the hcore Fock matrix
template<typename T, std::size_t NDIM>
Tensor<T> compute_fock_matrix(World& world, std::shared_ptr<NuclearCorrelationFactor>& ncf,
                              const MolecularOrbitals<T, NDIM>& mo) {
    std::shared_ptr<Fock<double, 3> > fock(new Fock<double, 3>(world));
    fock->add_operator("V", std::make_shared<Nuclear<double, 3> >(world, ncf));
    fock->add_operator("T", std::make_shared<Kinetic<double, 3> >(world));
    Tensor<T> f = (*fock)(mo.get_mos() * ncf->square(), mo.get_mos());
    return f;
}

/// distance between v1 and v2
double dist(const Vector<double,3> v1, const Vector<double,3> v2) {
    return (v1-v2).normf();
}

/// angle in degrees between v1--v2--v3 (v2 being the center)
double bend(const Vector<double,3> v1, const Vector<double,3> v2, const Vector<double,3> v3) {
    Vector<double,3> w1=v1-v2;
    Vector<double,3> w2=v3-v2;
    double denom=inner(w1,w2);
    double num=w1.normf()*w2.normf();
    return acos(denom/num)/constants::pi*180.0;
}

template<typename T, std::size_t NDIM>
MolecularOrbitals<T, NDIM>
compute_initial_orbitals(World& world, const AtomicBasisSet& aobasis, const Molecule& molecule,
                         std::shared_ptr<NuclearCorrelationFactor>& ncf) {
    std::vector<real_function_3d> aos = SCF::project_ao_basis_only(world, aobasis, molecule);
    MolecularOrbitals<double, 3> ao;
    ao.set_mos(aos * ncf->inverse());
    Tensor<T> fock = compute_fock_matrix(world, ncf, ao);
    Tensor<T> overlap = matrix_inner(world, ao.get_mos() * ncf->square(), ao.get_mos());
    Tensor<T> U;
    Tensor<typename Tensor<T>::scalar_type> evals;

    sygvp(world, fock, overlap, 1, U, evals);
    Tensor<T> UT = transpose(U);
    print("aobasis fock matrix");
    print(fock);
    print("orbital energies");
    print(evals);

    MolecularOrbitals<T, NDIM> mos;
    const long nmo = molecule.total_nuclear_charge() / 2;
    std::vector<Function<T, NDIM>> v = transform(world, ao.get_mos(), U(_, Slice(0, nmo - 1)));

    mos.update_mos_and_eps(v, evals(Slice(0, nmo - 1)));
    mos.set_all_orbitals_occupied();
    mos.recompute_localize_sets();
    return mos;
}

template<typename T, std::size_t NDIM>
std::tuple<Localizer<T,NDIM>, MolecularOrbitals<T,NDIM>, std::shared_ptr<NuclearCorrelationFactor>>
prepare_calculation(World& world, Molecule& molecule) {

    molecule.orient();
    molecule.print();
    std::shared_ptr<PotentialManager> pm(new PotentialManager(molecule, ""));
    pm->make_nuclear_potential(world);
//        auto ncf = create_nuclear_correlation_factor(world, molecule, pm, std::make_pair("slater", 2.0));
    auto ncf = create_nuclear_correlation_factor(world, molecule, pm, std::make_pair("none", 1.0));
    ncf->initialize(FunctionDefaults<NDIM>::get_thresh() * 0.1);

    AtomicBasisSet aobasis;
//    aobasis.read_file("sto-3g");
    aobasis.read_file("6-31G");
    std::vector<int> at_to_bf, at_nbf;
    aobasis.atoms_to_bfn(molecule, at_to_bf, at_nbf);

    std::vector<real_function_3d> aos = SCF::project_ao_basis_only(world, aobasis, molecule);
    Localizer<double, 3> localizer(world, aobasis, molecule, aos);
    localizer.set_metric(ncf->function());

    auto mos = compute_initial_orbitals<double, 3>(world, aobasis, molecule, ncf);
    mos.pretty_print("initial mos");

    return std::tuple(localizer,mos,ncf);
}

/// test localized orbitals: should be pointing towards the edges of a tetrahedron
bool test_ne_boys(World& world) {
    test_output tout("testing ne boys localization");
    Molecule ne_mol;
    ne_mol.read_structure_from_library("ne");
    auto [localizer,mos,ncf]=prepare_calculation<double,3>(world,ne_mol);
    localizer.set_method("boys");
    localizer.set_enforce_core_valence_separation(true);
    Tensor<double> fock1 = compute_fock_matrix(world, ncf, mos);
    Tensor<double> overlap = matrix_inner(world,ncf->square()*mos.get_mos(),mos.get_mos());

    auto mo1=localizer.localize(mos, fock1, overlap, true);
//    Tensor<double> fock2 = compute_fock_matrix(world, ncf, mo1);
//    print("final fock matrix");
//    print(fock2);
    std::vector<Vector<double,3>> center=mo1.compute_center(ncf->square());
    ne_mol.print();
    for (int i=0; i<mo1.get_mos().size(); ++i) {
        print("center of mo",i,center[i][0],center[i][1],center[i][2]);
    }
    print("distances ne tetrahedron",dist(center[1],center[2]),
          dist(center[1],center[3]),
          dist(center[1],center[4]),
          dist(center[2],center[3]),
          dist(center[2],center[4]),
          dist(center[3],center[4]));
    print("angles ne tetrahedron",
          bend(center[1],center[0],center[2]),
          bend(center[1],center[0],center[3]),
          bend(center[1],center[0],center[4]),
          bend(center[2],center[0],center[3]),
          bend(center[2],center[0],center[4]),
          bend(center[3],center[0],center[4]));

    bool success=true;
    success=success and (fabs(dist(center[1],center[2])-0.7158)<0.001);
    success=success and (fabs(dist(center[1],center[3])-0.7158)<0.001);
    success=success and (fabs(dist(center[1],center[4])-0.7158)<0.001);
    success=success and (fabs(dist(center[2],center[3])-0.7158)<0.001);
    success=success and (fabs(dist(center[2],center[4])-0.7158)<0.001);
    success=success and (fabs(dist(center[3],center[4])-0.7158)<0.001);

    success=success and (fabs(bend(center[1],center[0],center[2])-109.4712)<0.001);
    success=success and (fabs(bend(center[1],center[0],center[3])-109.4712)<0.001);
    success=success and (fabs(bend(center[1],center[0],center[4])-109.4712)<0.001);
    success=success and (fabs(bend(center[2],center[0],center[3])-109.4712)<0.001);
    success=success and (fabs(bend(center[2],center[0],center[4])-109.4712)<0.001);
    success=success and (fabs(bend(center[3],center[0],center[4])-109.4712)<0.001);

    tout.end(success);
    return success;

}

/// test localized orbitals: should be pointing towards the edges of a tetrahedron
bool test_ethylene(World& world) {
    test_output tout("testing ethylene boys localization");
    Molecule ethylene_mol;
    ethylene_mol.read_structure_from_library("c2h4");
    auto [localizer,mos,ncf]=prepare_calculation<double,3>(world,ethylene_mol);
    localizer.set_method("boys");
    localizer.set_enforce_core_valence_separation(false);
    Tensor<double> fock1 = compute_fock_matrix(world, ncf, mos);
    Tensor<double> overlap = matrix_inner(world,ncf->square()*mos.get_mos(),mos.get_mos());

    auto mo1=localizer.localize(mos, fock1, overlap, true);
    for (int i=0; i<mo1.get_mos().size(); ++i) {
        std::vector<std::string> molecular_info=ethylene_mol.cubefile_header();
        std::string filename="new_mo_ethylene"+std::to_string(i)+".cube";
        plot_cubefile<3>(world,mo1.get_mos()[i],filename,molecular_info);

    }

//    Tensor<double> fock2 = compute_fock_matrix(world, ncf, mo1);
//    print("final fock matrix");
//    print(fock2);
    std::vector<Vector<double,3>> center=mo1.compute_center(ncf->square());
    ethylene_mol.print();
    for (int i=0; i<mo1.get_mos().size(); ++i) {
        char buf[1024];
        std::sprintf(buf,"center of mo %2d: %8.4f  %8.4f  %8.4f",i,center[i][0],center[i][1],center[i][2]);
        print(buf);
    }
    bool success=false;
    tout.end(success);
    return success;

}

template<typename T, std::size_t NDIM>
int test_localization(World& world, Localizer<T, NDIM>& localizer, std::shared_ptr<NuclearCorrelationFactor>& ncf,
                      const MolecularOrbitals<T, NDIM>& mo, double sum_orbital_energy, const bool verbose) {
    test_output tout("testing "+localizer.get_method(), verbose);
    int success = 0;
    localizer.set_enforce_core_valence_separation(false);
    MolecularOrbitals<T, NDIM> lmo = localizer.localize(mo, true);
    Tensor<T> fock = compute_fock_matrix(world, ncf, lmo);
    print(localizer.get_method(), "localized fock matrix");
    print(fock);

    double trace = 0.0;
    for (int i = 0; i < fock.dim(0); ++i) trace += fock(i, i);
    print("sum over diagonal fock matrix elements", trace);
    success+=(std::fabs(trace-sum_orbital_energy)/trace<FunctionDefaults<3>::get_thresh());

    tout.end(success);
    return success;
}

template<typename T, std::size_t NDIM>
int test_core_valence_separation(World& world, Localizer<T, NDIM>& localizer,
                                 std::shared_ptr<NuclearCorrelationFactor>& ncf,
                                 const MolecularOrbitals<T, NDIM>& mo1,
                                 const double sum_orbital_energy,
                                 const bool verbose) {
    bool success = 0;
    std::string method=localizer.get_method();
    test_output tout("testing core-valence separation "+method,verbose);
    localizer.set_enforce_core_valence_separation(true);
    Tensor<T> fock1 = compute_fock_matrix(world, ncf, mo1);
    Tensor<T> overlap = matrix_inner(world,ncf->square()*mo1.get_mos(),mo1.get_mos());
    print(method, "diagonal fock matrix to start with");
    print(fock1);
    MolecularOrbitals<T, NDIM> lmo = localizer.localize(mo1, fock1, overlap, true);
//    Tensor<T> UT = localizer.compute_core_valence_separation_transformation_matrix(world, mo, fock1, overlap);
//    std::vector<Function<T, NDIM>> result = transform(world, mo.get_mos(), UT);
//    truncate(world, result);
//    MolecularOrbitals<T, NDIM> lmo;
//    lmo.set_mos(result);

    Tensor<T> fock = compute_fock_matrix(world, ncf, lmo);
    print(method, "localized fock matrix");
    print(fock);
    bool success1=Localizer<T,NDIM>::check_core_valence_separation(fock,lmo.get_localize_sets());
    print("success cv-separation",success1);
    success +=success1;


    double trace = 0.0;
    for (int i = 0; i < fock.dim(0); ++i) trace += fock(i, i);
    print("sum over diagonal fock matrix elements", trace);
    success+=(std::fabs(trace-sum_orbital_energy)/trace<FunctionDefaults<3>::get_thresh());

    tout.end(success==0);
    return (success==0);
}

int main(int argc, char **argv) {
    int result = 0;
    {
        World& world = madness::initialize(argc, argv);
        startup(world, argc, argv);
        print("entering test_localizer");
        commandlineparser parser(argc,argv);
        parser.print_map();
        bool verbose=parser.key_exists("verbose");
        const int k = 9;
        const double thresh = 1.e-6;
//        const double L = 20.0;
        const double L = 5.1427e+01;
        FunctionDefaults<3>::set_cubic_cell(-L, L);
        FunctionDefaults<3>::set_thresh(thresh);
        FunctionDefaults<3>::set_k(k);

        test_ne_boys(world);
        test_ethylene(world);

        /*
         * test orbital invariance of the orbital energy sum
         */
//
//        double sum_orbital_energies=mos.get_eps().sum();
//        print("canonical sum over orbital energies",sum_orbital_energies);
//
//        localizer.set_method("boys");
//        test_localization(world, localizer, ncf, mos, sum_orbital_energies, verbose);
//        test_core_valence_separation(world, localizer, ncf, mos, sum_orbital_energies, verbose);
//
//        localizer.set_method("pm");
//        test_localization(world, localizer, ncf, mos, sum_orbital_energies, verbose);
//        test_core_valence_separation(world, localizer, ncf, mos, sum_orbital_energies, verbose);
//
//        localizer.set_method("new");
//        test_localization(world, localizer, ncf, mos, sum_orbital_energies, verbose);
//        test_core_valence_separation(world, localizer, ncf, mos, sum_orbital_energies, verbose);

        print("result", result);
    }
    madness::finalize();
    return result;

}
