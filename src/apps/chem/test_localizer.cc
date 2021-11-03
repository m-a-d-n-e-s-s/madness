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
int test_localization(World& world, Localizer<T, NDIM>& localizer, std::shared_ptr<NuclearCorrelationFactor>& ncf,
                      const MolecularOrbitals<T, NDIM>& mo, std::string method) {
    int success = 0;
    MolecularOrbitals<T, NDIM> lmo = localizer.localize(mo, method, 1.e-4, true);
    Tensor<T> fock = compute_fock_matrix(world, ncf, lmo);
    print(method, "localized fock matrix");
    print(fock);

    double trace = 0.0;
    for (int i = 0; i < fock.dim(0); ++i) trace += fock(i, i);
    print("sum over diagonal fock matrix elements", trace);

    return success;
}

template<typename T, std::size_t NDIM>
int test_core_valence_separation(World& world, Localizer<T, NDIM>& localizer,
                                 std::shared_ptr<NuclearCorrelationFactor>& ncf,
                                 const MolecularOrbitals<T, NDIM>& mo1, std::string method) {
    int success = 0;
    MolecularOrbitals<T, NDIM> mo = localizer.localize(mo1, method, 1.e-4, true);
    Tensor<T> fock1 = compute_fock_matrix(world, ncf, mo);
    print(method, "localized fock matrix to start with");
    print(fock1);
    Tensor<T> overlap = matrix_inner(world,ncf->square()*mo.get_mos(),mo.get_mos());
    double thresh_degenerate=FunctionDefaults<3>::get_thresh();
    double tolloc=FunctionDefaults<3>::get_thresh();
    DistributedMatrix<T> dUT = localizer.compute_core_valence_separation_transformation_matrix(world,
                             mo, fock1, overlap, thresh_degenerate,method,tolloc,false);
    std::vector<Function<T, NDIM>> result = transform(world, mo.get_mos(), dUT);
    truncate(world, result);
    MolecularOrbitals<T, NDIM> lmo;
    lmo.set_mos(result);

    Tensor<T> fock = compute_fock_matrix(world, ncf, lmo);
    print(method, "localized fock matrix");
    print(fock);

    double trace = 0.0;
    for (int i = 0; i < fock.dim(0); ++i) trace += fock(i, i);
    print("sum over diagonal fock matrix elements", trace);

    return success;
}

int main(int argc, char **argv) {
    int result = 0;
    {
        World& world = madness::initialize(argc, argv);
        startup(world, argc, argv);
        print("entering test_localizer");
        commandlineparser parser(argc,argv);
        const int k = 9;
        const double thresh = 1.e-6;
        const double L = 24.0;
        FunctionDefaults<3>::set_cubic_cell(-L, L);
        FunctionDefaults<3>::set_thresh(thresh);
        FunctionDefaults<3>::set_k(k);


        Molecule molecule;
        std::string geometry = parser.key_exists("structure") ? parser.value("structure") : "h2o";
        molecule.read_structure_from_library(geometry);
        molecule.print();
        std::shared_ptr<PotentialManager> pm(new PotentialManager(molecule, ""));
        auto ncf = create_nuclear_correlation_factor(world, molecule, pm, std::make_pair("slater", 2.0));
        ncf->initialize(thresh * 0.1);

        AtomicBasisSet aobasis;
        aobasis.read_file("6-31g");
        std::vector<int> at_to_bf, at_nbf;
        aobasis.atoms_to_bfn(molecule, at_to_bf, at_nbf);

        std::vector<real_function_3d> aos = SCF::project_ao_basis_only(world, aobasis, molecule);
        Localizer<double, 3> localizer(world, aobasis, molecule, aos);
        localizer.set_metric(ncf->function());

        auto mos = compute_initial_orbitals<double, 3>(world, aobasis, molecule, ncf);
        mos.pretty_print("initial mos");

        auto blocks=localizer.convert_set_to_slice(mos.get_localize_sets());
        for (auto block : blocks) print("block",block);

        test_localization(world, localizer, ncf, mos, "boys");
        test_localization(world, localizer, ncf, mos, "pm");
        test_localization(world, localizer, ncf, mos, "new");
        test_core_valence_separation(world, localizer, ncf, mos, "boys");
        test_core_valence_separation(world, localizer, ncf, mos, "pm");
        test_core_valence_separation(world, localizer, ncf, mos, "new");

        print("result", result);
    }
    madness::finalize();
    return result;

}
