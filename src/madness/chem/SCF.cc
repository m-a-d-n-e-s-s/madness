/*
  This file is part of MADNESS.

  Copyright (C) 2007,2010 Oak Ridge National Laboratory

  This program is free software; you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation; either version 2 of the License, or
  (at your option) any later version.

  This program is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
  GNU General Public License for more details.

  You should have received a copy of the GNU General Public License
  along with this program; if not, write to the Free Software
  Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307 USA

  For more information please contact:

  Robert J. Harrison
  Oak Ridge National Laboratory
  One Bethel Valley Road
  P.O. Box 2008, MS-6367

  email: harrisonrj@ornl.gov
  tel:   865-241-3937
  fax:   865-572-0680
 */

/// \file SCF.cc
/// \brief Molecular HF and DFT code
/// \defgroup moldft The molecular density functional and Hartree-Fock code



#include <madness/world/worldmem.h>
#include <madness.h>
#include <madness/chem/SCF.h>
#include <madchem.h>

#if defined(__has_include)
#  if __has_include(<filesystem>)
#    define MADCHEM_HAS_STD_FILESYSTEM
// <filesystem> is not reliably usable on Linux with gcc < 9
#    if defined(__GNUC__)
#      if __GNUC__ >= 7 && __GNUC__ < 9
#        undef MADCHEM_HAS_STD_FILESYSTEM
#      endif
#    endif
#    if defined(MADCHEM_HAS_STD_FILESYSTEM)

#      include <filesystem>

#    endif
#  endif
#endif

using namespace madchem;
namespace madness {

//    // moved to vmra.h
//    template <typename T, std::size_t NDIM>
//    DistributedMatrix<T> matrix_inner(const DistributedMatrixDistribution& d,
//                                      const std::vector< Function<T,NDIM> >& f,
//                                      const std::vector< Function<T,NDIM> >& g,
//                                      bool sym=false)


template<typename T, std::size_t NDIM>
static void verify_tree(World& world, const std::vector<Function<T, NDIM> >& v) {
    for (unsigned int i = 0; i < v.size(); i++) {
        v[i].verify_tree();
    }
}

template<int NDIM>
struct unaryexp {
    void operator()(const Key<NDIM>& key, Tensor<double_complex>& t) const {
        //vzExp(t.size, t.ptr(), t.ptr());
        UNARY_OPTIMIZED_ITERATOR(double_complex, t, *_p0 = exp(*_p0););
    }

    template<typename Archive>
    void serialize(Archive& ar) {}
};


static double rsquared(const coordT& r) {
    return r[0] * r[0] + r[1] * r[1] + r[2] * r[2];
}

// // Returns exp(-I*t*V)
// static Function<double_complex, 3> make_exp(double t, const Function<double, 3>& v) {
//     v.reconstruct();
//     Function<double_complex, 3> expV = double_complex(0.0, -t) * v;
//     expV.unaryop(unaryexp<3>());
//     //expV.truncate(); expV.reconstruct();
//     return expV;
// }

// Timer modified to correctly nest
static bool print_timings = false;
static std::vector<double> ttt, sss;

static void START_TIMER(World& world) {
    world.gop.fence();
    ttt.push_back(wall_time());
    sss.push_back(cpu_time());
}

static double pop(std::vector<double>& v) {
    MADNESS_ASSERT(v.size());
    double x = v.back();
    v.pop_back();
    return x;
}

static void END_TIMER(World& world, const char *msg) {
    double wall = wall_time() - pop(ttt), cpu = cpu_time() - pop(sss);
    if (world.rank() == 0 and print_timings) printf("timer: %20.20s %8.2fs %8.2fs\n", msg, cpu, wall);
}


/// Given overlap matrix, return rotation with 3rd order error to orthonormalize the vectors
tensorT Q3(const tensorT& s) {
    tensorT Q = inner(s, s);
    Q.gaxpy(0.2, s, -2.0 / 3.0);
    for (int i = 0; i < s.dim(0); ++i) Q(i, i) += 1.0;
    return Q.scale(15.0 / 8.0);
}

/// Given overlap matrix, return rotation with 2nd order error to orthonormalize the vectors
tensorT Q2(const tensorT& s) {
    tensorT Q = -0.5 * s;
    for (int i = 0; i < s.dim(0); ++i) Q(i, i) += 1.5;
    return Q;
}

void SCF::output_scf_info_schema(const std::map<std::string, double> &vals,
                                 const tensorT &dipole_T) const {
    nlohmann::json j = {};
    // if it exists figure out the size.  pushback for each protocol
    const double thresh = FunctionDefaults<3>::get_thresh();
    const int k = FunctionDefaults<3>::get_k();
    j["scf_threshold"] = thresh;
    j["scf_k"] = k;
    for (auto const &[key, val]: vals) {
        j[key] = val;
    }
    j["scf_dipole_moment"] = tensor_to_json(dipole_T);
    int num = 0;
    update_schema(param.prefix()+".scf_info", j);
}

void SCF::output_calc_info_schema() const {
    nlohmann::json j = {};
    World& world=amo.front().world();
    if (world.rank()==0) {
        vec_pair_ints int_vals;
        vec_pair_T<double> double_vals;
        vec_pair_tensor_T<double> double_tensor_vals;


        int_vals.emplace_back("calcinfo_nmo", param.nmo_alpha() + param.nmo_beta());
        int_vals.emplace_back("calcinfo_nalpha", param.nalpha());
        int_vals.emplace_back("calcinfo_nbeta", param.nbeta());
        int_vals.emplace_back("calcinfo_natom", molecule.natom());
        int_vals.emplace_back("k", FunctionDefaults<3>::get_k());

        to_json(j, int_vals);
        //    double_vals.push_back({"return_energy", value(molecule.get_all_coords().flat())});
        double_vals.emplace_back("return_energy", current_energy);
        to_json(j, double_vals);
        double_tensor_vals.emplace_back("scf_eigenvalues_a", aeps);
        if (param.nbeta() != 0 && !param.spin_restricted()) {
            double_tensor_vals.emplace_back("scf_eigenvalues_b", beps);
        }

        to_json(j, double_tensor_vals);
        param.to_json(j);
        e_data.to_json(j);

        //    output_schema(param.prefix()+".calc_info", j);
        update_schema(param.prefix()+".calc_info", j);
    }
}

void scf_data::add_data(std::map<std::string, double> values) {
    //print("ADDING DATA");

    iter++;
    std::for_each(e_data.begin(), e_data.end(), [&values](auto &v) {
    //    print(v.first, " : ", values[v.first]);
        v.second.push_back(values[v.first]);
    });
}

scf_data::scf_data() : iter(0) {

    e_data.insert({"e_kinetic", std::vector<double>(0)});
    e_data.insert({"e_local", std::vector<double>(0)});
    e_data.insert({"e_nuclear", std::vector<double>(0)});
    e_data.insert({"e_coulomb", std::vector<double>(0)});
    e_data.insert({"e_pcm", std::vector<double>(0)});
    e_data.insert({"e_xc", std::vector<double>(0)});
    e_data.insert({"e_nrep", std::vector<double>(0)});
    e_data.insert({"e_tot", std::vector<double>(0)});
}


void scf_data::to_json(json &j) const {
    madness::print("SCF DATA TO JSON");

    j["scf_e_data"] = json();
    j["scf_e_data"]["iterations"] = iter;

    for (const auto &e: e_data) {
        //::print(e.second);
        j["scf_e_data"].push_back({e.first, e.second});
    }
}

void scf_data::print_data() {
    for (const auto &[key, value]: e_data) { print(key, " : ", value); }
}

void scf_data::add_gradient(const Tensor<double> &grad) {
    gradient = tensor_to_json(grad);
}


//    SCF::SCF(World & world, const char *filename) : SCF(world, (world.rank() == 0 ? std::make_shared<std::ifstream>(filename) : nullptr)){
//    }

/// collective constructor, reads \c input on rank 0, broadcasts to all
SCF::SCF(World& world, const commandlineparser& parser) : param(CalculationParameters(world, parser)) {
    PROFILE_MEMBER_FUNC(SCF);

    molecule=Molecule(world,parser);

//    param.read(world,parser.value("input"),"dft");
    if (world.rank() == 0) {

//        // read input parameters from the input file
//        if (parser.key_exists("structure"))
//            param.set_user_defined_value("molecular_structure", parser.value("structure"));
//
//        std::string molecular_structure = param.get<std::string>("molecular_structure");
//        if (molecular_structure == "inputfile") {
//            std::ifstream ifile(parser.value("input"));
//            molecule.read(ifile);
//        } else {
//            molecule.read_structure_from_library(molecular_structure);
//        }
//
//        // set derived parameters for the molecule
//
//        //if psp_calc is true, set all atoms to PS atoms
//        //if not, check whether some atoms are PS atoms or if this a pure AE calculation
//        if (param.get<bool>("psp_calc")) {
//            for (size_t iatom = 0; iatom < molecule.natom(); iatom++) {
//                molecule.set_pseudo_atom(iatom, true);
//            }
//        }
//
//        //modify atomic charge for complete PSP calc or individual PS atoms
//        for (size_t iatom = 0; iatom < molecule.natom(); iatom++) {
//            if (molecule.get_pseudo_atom(iatom)) {
//                unsigned int an = molecule.get_atomic_number(iatom);
//                double zeff = get_charge_from_file("gth.xml", an);
//                molecule.set_atom_charge(iatom, zeff);
//            }
//        }

        if (molecule.parameters.core_type() != "none") {
            molecule.read_core_file(molecule.parameters.core_type());
        }

//        if (not molecule.parameters.no_orient()) molecule.orient();

        //account for nwchem aobasis generation
        if (param.nwfile() == "none") reset_aobasis(param.aobasis());
        else aobasis.read_nw_file(param.nwfile());
        param.set_derived_values(molecule, aobasis, parser);

    }
    world.gop.broadcast_serializable(molecule, 0);
    world.gop.broadcast_serializable(param, 0);
    world.gop.broadcast_serializable(aobasis, 0);

    if (param.print_level() > 2) print_timings = true;

    xc.initialize(param.xc(), !param.spin_restricted(), world, param.print_level() >= 10);
    //xc.plot();

    FunctionDefaults<3>::set_cubic_cell(-param.L(), param.L());
    //set_protocol < 3 > (world, param.econv());
    FunctionDefaults<3>::set_truncate_mode(1);

}

void SCF::set_print_timings(const bool value) {
    print_timings=value;
}

void SCF::copy_data(World& world, const SCF& other) {
    aeps = copy(other.aeps);
    beps = copy(other.beps);
    aocc = copy(other.aocc);
    bocc = copy(other.bocc);
    amo = copy(world, other.amo);
    bmo = copy(world, other.bmo);
    aset = other.aset;
    bset = other.bset;
    ao = copy(world, other.ao);
    at_to_bf = other.at_to_bf;
    at_nbf = other.at_nbf;
}

void SCF::save_mos(World& world) {
    PROFILE_MEMBER_FUNC(SCF);
    auto archivename=param.prefix()+".restartdata";
    archive::ParallelOutputArchive<archive::BinaryFstreamOutputArchive> ar(world, archivename.c_str(), param.get<int>("nio"));
    // IF YOU CHANGE ANYTHING HERE MAKE SURE TO UPDATE THIS VERSION NUMBER
    /*
     * After spin restricted
      double L;
      int k;
      Molecule molecule;
      std::string xc;
      */
    unsigned int version = 4;
    ar & version;
    ar & current_energy & param.spin_restricted();
    ar & param.L() & FunctionDefaults<3>::get_k() & molecule & param.xc() & param.localize_method() & converged_for_thresh;
    // Re order so it doesn't effect orbital data

    ar & (unsigned int) (amo.size());
    ar & aeps & aocc & aset;
    for (unsigned int i = 0; i < amo.size(); ++i) ar & amo[i];
    if (!param.spin_restricted()) {
        ar & (unsigned int) (bmo.size());
        ar & beps & bocc & bset;
        for (unsigned int i = 0; i < bmo.size(); ++i) ar & bmo[i];
    }

    // Do not make a restartaodata file if nwchem orbitals used,
    // as no aoamo/aobmo overlap matrix can be computed
    if (param.nwfile() == "none") {
        tensorT Saoamo = matrix_inner(world, ao, amo);
        tensorT Saobmo = (!param.spin_restricted()) ? matrix_inner(world, ao, bmo) : tensorT();
        if (world.rank() == 0) {
            archive::BinaryFstreamOutputArchive arao(param.prefix()+".restartaodata");
            arao << Saoamo << aeps << aocc << aset;
            if (!param.spin_restricted()) arao << Saobmo << beps << bocc << bset;
        }
    }
}

void SCF::load_mos(World& world) {
    PROFILE_MEMBER_FUNC(SCF);
    //        const double trantol = vtol / std::min(30.0, double(param.nalpha));
    const double thresh = FunctionDefaults<3>::get_thresh();
    const int k = FunctionDefaults<3>::get_k();
    unsigned int nmo = 0;
    bool spinrest = false;

    amo.clear();
    bmo.clear();

    archive::ParallelInputArchive<archive::BinaryFstreamInputArchive> ar(world, param.prefix()+".restartdata");

    /*
      File format:
          unsigned int version;
          double current energy;
      bool spinrestricted --> if true only alpha orbitals are present
      double L;
      int k;
      Molecule molecule;
      std::string xc;
    std::string localize;
      unsigned int nmo_alpha;
      Tensor<double> aeps;
      Tensor<double> aocc;
      vector<int> aset;
      for i from 0 to nalpha-1:
      .   Function<double,3> amo[i]
      repeat for beta if !spinrestricted
     */
    // Local copies for a basic check
    double L=0;
    int k1=0;                    // Ignored for restarting, used in response only
    unsigned int version = 4;// UPDATE THIS IF YOU CHANGE ANYTHING
    unsigned int archive_version=0;

    ar & archive_version;

    if (archive_version != version) {
        if (world.rank() == 0)
            print(
                    "Loading from a different version of archive. Archive version", archive_version, "MADNESS version",
                    version);
        throw "Invalid archive";
    }

    // LOTS OF LOGIC MISSING HERE TO CHANGE OCCUPATION NO., SET,
    // EPS, SWAP, ... sigh
    ar & current_energy & spinrest;
    // Reorder
    ar & L & k1 & molecule & param.xc() & param.localize_method() & converged_for_thresh;

    ar & nmo;
    MADNESS_ASSERT(nmo >= unsigned(param.nmo_alpha()));
    ar & aeps & aocc & aset;
    // Some basic checks
    if (L != param.L()) {
        if (world.rank() == 0)
            print(
                    "Warning: Box size mismatch between archive and input parameter. "
                    "Archive value",
                    L,
                    "Param value",
                    param.L());
        throw "Mismatch in box sizes";
    }
    if (world.rank() == 0) {
        print("Restarting from this molecular geometry");
        molecule.print();
    }
    amo.resize(nmo);
    for (unsigned int i = 0; i < amo.size(); ++i) ar & amo[i];
    unsigned int n_core = molecule.n_core_orb_all();
    if (nmo > unsigned(param.nmo_alpha())) {
        aset = vector<int>(aset.begin() + n_core, aset.begin() + n_core + param.nmo_alpha());
        amo = vecfuncT(amo.begin() + n_core, amo.begin() + n_core + param.nmo_alpha());
        aeps = copy(aeps(Slice(n_core, n_core + param.nmo_alpha() - 1)));
        aocc = copy(aocc(Slice(n_core, n_core + param.nmo_alpha() - 1)));
    }

    if (amo[0].k() != k) {
        reconstruct(world, amo);
        for (unsigned int i = 0; i < amo.size(); ++i) amo[i] = madness::project(amo[i], k, thresh, false);
        world.gop.fence();
    }
    set_thresh(world, amo, thresh);

    //        normalize(world, amo);
    //        amo = transform(world, amo, Q3(matrix_inner(world, amo, amo)),
    //        trantol, true); truncate(world, amo); normalize(world, amo);

    if (!param.spin_restricted()) {
        if (spinrest) {  // Only alpha spin orbitals were on disk
            MADNESS_ASSERT(param.nmo_alpha() >= param.nmo_beta());
            bmo.resize(param.nmo_beta());
            bset.resize(param.nmo_beta());
            beps = copy(aeps(Slice(0, param.nmo_beta() - 1)));
            bocc = copy(aocc(Slice(0, param.nmo_beta() - 1)));
            for (int i = 0; i < param.nmo_beta(); ++i) bmo[i] = copy(amo[i]);
        } else {
            ar & nmo;
            ar & beps & bocc & bset;

            bmo.resize(nmo);
            for (unsigned int i = 0; i < bmo.size(); ++i) ar & bmo[i];

            if (nmo > unsigned(param.nmo_beta())) {
                bset = vector<int>(bset.begin() + n_core, bset.begin() + n_core + param.nmo_beta());
                bmo = vecfuncT(bmo.begin() + n_core, bmo.begin() + n_core + param.nmo_beta());
                beps = copy(beps(Slice(n_core, n_core + param.nmo_beta() - 1)));
                bocc = copy(bocc(Slice(n_core, n_core + param.nmo_beta() - 1)));
            }

            if (bmo[0].k() != k) {
                reconstruct(world, bmo);
                for (unsigned int i = 0; i < bmo.size(); ++i) bmo[i] = madness::project(bmo[i], k, thresh, false);
                world.gop.fence();
            }
            set_thresh(world, amo, thresh);

            //                normalize(world, bmo);
            //                bmo = transform(world, bmo, Q3(matrix_inner(world, bmo,
            //                bmo)), trantol, true); truncate(world, bmo);
            //                normalize(world, bmo);
        }
    }
}

void SCF::do_plots(World& world) {
    PROFILE_MEMBER_FUNC(SCF);
    START_TIMER(world);

    std::vector<long> npt(3, static_cast<long>(param.get<int>("npt_plot")));

    if (param.plot_cell().size() == 0)
        param.plot_cell() = copy(FunctionDefaults<3>::get_cell());

    if (param.get<bool>("plotdens") || param.get<bool>("plotcoul")) {
        functionT rho;
        rho = make_density(world, aocc, amo);

        if (param.spin_restricted()) {
            rho.scale(2.0);
        } else {
            functionT rhob = make_density(world, bocc, bmo);
            functionT rho_spin = rho - rhob;
            rho += rhob;
            plotdx(rho_spin, "spin_density.dx", param.plot_cell(), npt, true);

        }
        plotdx(rho, "total_density.dx", param.plot_cell(), npt, true);
        if (param.get<bool>("plotcoul")) {
            real_function_3d vnuc = potentialmanager->vnuclear();
            functionT vlocl = vnuc + apply(*coulop, rho);
            vlocl.truncate();
            vlocl.reconstruct();
            plotdx(vlocl, "coulomb.dx", param.plot_cell(), npt, true);
        }
    }

    for (int i = param.get<int>("plotlo"); i <= param.get<int>("plothi"); ++i) {
        std::size_t bufsize=256;
        char fname[bufsize];
        if (i < param.nalpha()) {
            snprintf(fname,bufsize, "amo-%5.5d.dx", i);
            plotdx(amo[i], fname, param.plot_cell(), npt, true);
        }
        if (!param.spin_restricted() && i < param.nbeta()) {
            snprintf(fname,bufsize, "bmo-%5.5d.dx", i);
            plotdx(bmo[i], fname, param.plot_cell(), npt, true);
        }
    }
    END_TIMER(world, "plotting");
}

void SCF::project(World& world) {
    PROFILE_MEMBER_FUNC(SCF);
    reconstruct(world, amo);
    for (unsigned int i = 0; i < amo.size(); ++i) {
        amo[i] = madness::project(amo[i], FunctionDefaults<3>::get_k(),
                                  FunctionDefaults<3>::get_thresh(), false);
    }
    world.gop.fence();
    truncate(world, amo);
    normalize(world, amo);
    if (param.nbeta() && !param.spin_restricted()) {
        reconstruct(world, bmo);
        for (unsigned int i = 0; i < bmo.size(); ++i) {
            bmo[i] = madness::project(bmo[i], FunctionDefaults<3>::get_k(),
                                      FunctionDefaults<3>::get_thresh(), false);
        }
        world.gop.fence();
        truncate(world, bmo);
        normalize(world, bmo);
    }
}

void SCF::make_nuclear_potential(World& world) {
    PROFILE_MEMBER_FUNC(SCF);
    START_TIMER(world);
    potentialmanager = std::shared_ptr<PotentialManager
    >(new PotentialManager(molecule, molecule.parameters.core_type()));
    gthpseudopotential = std::shared_ptr<GTHPseudopotential<double>
    >(new GTHPseudopotential<double>(world, molecule));

    if (!molecule.parameters.pure_ae()) {
        gthpseudopotential->make_pseudo_potential(world);
    }
    if (!molecule.parameters.psp_calc()) {
        potentialmanager->make_nuclear_potential(world);
    }
    END_TIMER(world, "Project vnuclear");
}

vecfuncT SCF::project_ao_basis(World& world, const AtomicBasisSet& aobasis) {
    PROFILE_MEMBER_FUNC(SCF);
    // Make at_to_bf, at_nbf ... map from atom to first bf on atom, and nbf/atom
    aobasis.atoms_to_bfn(molecule, at_to_bf, at_nbf);

    return SCF::project_ao_basis_only(world, aobasis, molecule);
}

vecfuncT SCF::project_ao_basis_only(World& world, const AtomicBasisSet& aobasis,
                                    const Molecule& molecule) {
    vecfuncT ao = vecfuncT(aobasis.nbf(molecule));
    for (int i = 0; i < aobasis.nbf(molecule); ++i) {
        functorT aofunc(new AtomicBasisFunctor(
                aobasis.get_atomic_basis_function(molecule, i)));
        ao[i] = factoryT(world).functor(aofunc).truncate_on_project().nofence().truncate_mode(1);
    }
    world.gop.fence();
    truncate(world, ao);
    normalize(world, ao);
    return ao;
}

void SCF::analyze_vectors(World& world, const vecfuncT& mo, const tensorT& occ,
                          const tensorT& energy, const std::vector<int>& set) {
    START_TIMER(world);
    PROFILE_MEMBER_FUNC(SCF);
    tensorT Saomo = matrix_inner(world, ao, mo);
    tensorT Saoao = matrix_inner(world, ao, ao, true);
    int nmo1 = mo.size();
    tensorT rsq, dip(3, nmo1);
    {
        functionT frsq = factoryT(world).f(rsquared).initial_level(4);
        rsq = inner(world, mo, mul_sparse(world, frsq, mo, vtol));
        for (int axis = 0; axis < 3; ++axis) {
            functionT fdip = factoryT(world).functor(
                    functorT(new DipoleFunctor(axis))).initial_level(4);
            dip(axis, _) = inner(world, mo, mul_sparse(world, fdip, mo, vtol));
            for (int i = 0; i < nmo1; ++i)
                rsq(i) -= dip(axis, i) * dip(axis, i);

        }
    }
    tensorT C;
    END_TIMER(world, "Analyze vectors");

    START_TIMER(world);
    gesvp(world, Saoao, Saomo, C);
    END_TIMER(world, "Compute eigen gesv analyze vectors");
    C = transpose(C);
    long nmo = mo.size();
    size_t ncoeff = 0;
    for (long i = 0; i < nmo; ++i) {
        size_t ncoeffi = mo[i].size();
        ncoeff += ncoeffi;
        if (world.rank() == 0 and (param.print_level() > 1)) {
            printf("  MO%4ld : ", i);
            if (set.size())
                printf("set=%d : ", set[i]);

            if (occ.size())
                printf("occ=%.2f : ", occ(i));

            if (energy.size())
                printf("energy=%13.8f : ", energy(i));

            printf("ncoeff=%.2e:", (double) ncoeffi);

            printf("center=(%.2f,%.2f,%.2f) : radius=%.2f\n", dip(0, i),
                   dip(1, i), dip(2, i), sqrt(rsq(i)));
            aobasis.print_anal(molecule, C(i, _));
            printf("total number of coefficients = %.8e\n\n", double(ncoeff));
        }
    }
}

// this version is faster than the previous version on BG/Q
distmatT SCF::kinetic_energy_matrix(World& world, const vecfuncT& v) const {
    PROFILE_MEMBER_FUNC(SCF);
    int n = v.size();
    distmatT r = column_distributed_matrix<double>(world, n, n);
    START_TIMER(world);
    reconstruct(world, v);
    END_TIMER(world, "KEmat reconstruct");
    START_TIMER(world);
    vecfuncT dvx = apply(world, *(gradop[0]), v, false);
    vecfuncT dvy = apply(world, *(gradop[1]), v, false);
    vecfuncT dvz = apply(world, *(gradop[2]), v, false);
    world.gop.fence();
    END_TIMER(world, "KEmat differentiate");
    START_TIMER(world);
    compress(world, dvx, false);
    compress(world, dvy, false);
    compress(world, dvz, false);
    world.gop.fence();
    END_TIMER(world, "KEmat compress");
    START_TIMER(world);
    r += matrix_inner(r.distribution(), dvx, dvx, true);
    r += matrix_inner(r.distribution(), dvy, dvy, true);
    r += matrix_inner(r.distribution(), dvz, dvz, true);
    END_TIMER(world, "KEmat inner products");
    r *= 0.5;
    //tensorT p(v.size(),v.size());
    //r.copy_to_replicated(p);
    return r;
}

bool SCF::restart_aos(World& world) {
    tensorT Saoamo, Saobmo;
    bool OK = true;
    if (world.rank() == 0) {
        try {
            archive::BinaryFstreamInputArchive arao(param.prefix()+".restartaodata");
            arao >> Saoamo >> aeps >> aocc >> aset;
            if (Saoamo.dim(0) != int(ao.size()) || Saoamo.dim(1) != param.nmo_alpha()) {
                print(" AO alpha restart data size mismatch --- starting from atomic guess instead", Saoamo.dim(0),
                      ao.size(), Saoamo.dim(1), param.nmo_alpha());
                OK = false;
            }
            if (!param.spin_restricted()) {
                arao >> Saobmo >> beps >> bocc >> bset;
                if (Saobmo.dim(0) != int(ao.size()) || Saobmo.dim(1) != param.nmo_beta()) {
                    print(" AO beta restart data size mismatch --- starting from atomic guess instead", Saobmo.dim(0),
                          ao.size(), Saobmo.dim(1), param.nmo_beta());
                    OK = false;
                }
            }
            print("\nRestarting from AO projections on disk\n");
        }
        catch (...) {
            print("\nAO restart file open/reading failed --- starting from atomic guess instead\n");
            OK = false;
        }
    }
    int fred = OK;
    world.gop.broadcast(fred, 0);
    OK = fred;
    if (!OK) return false;

    world.gop.broadcast_serializable(Saoamo, 0);
    if (!param.spin_restricted()) world.gop.broadcast_serializable(Saobmo, 0);

    tensorT S = matrix_inner(world, ao, ao), c;

    gesvp(world, S, Saoamo, c);
    amo = transform(world, ao, c, vtol, true);
    truncate(world, amo);
    orthonormalize(world, amo, param.nalpha());

    if (!param.spin_restricted()) {
        gesvp(world, S, Saobmo, c);
        bmo = transform(world, ao, c, vtol, true);
        truncate(world, bmo);
        orthonormalize(world, bmo, param.nbeta());
    }

    return true;
}

void SCF::initial_guess(World& world) {
    PROFILE_MEMBER_FUNC(SCF);
    START_TIMER(world);
    if (param.restart()) {
        load_mos(world);
    } else {

        //If not using nwchem, proceed as normal...
        if (param.nwfile() == "none") {

            // recalculate initial guess density matrix without core orbitals
            if (!molecule.parameters.pure_ae()) {
                for (size_t iatom = 0; iatom < molecule.natom(); iatom++) {
                    if (molecule.get_pseudo_atom(iatom)) {
                        double zeff = molecule.get_atom_charge(iatom);
                        int atn = molecule.get_atomic_number(iatom);
                        aobasis.modify_dmat_psp(atn, zeff);
                    }
                }
            }

            // Use the initial density and potential to generate a better process map
            functionT rho =
                    factoryT(world).functor(
                            functorT(
                                    new MolecularGuessDensityFunctor(molecule,
                                                                     aobasis))).truncate_on_project();
            double nel = rho.trace();
            if (world.rank() == 0 and param.print_level() > 3)
                print("guess dens trace", nel);
            END_TIMER(world, "guess density");
            rho.scale(std::round(nel) / nel);

            if (world.size() > 1) {
                START_TIMER(world);
                LoadBalanceDeux<3> lb(world);
                real_function_3d vnuc;
                if (molecule.parameters.psp_calc()) {
                    vnuc = gthpseudopotential->vlocalpot();
                } else if (molecule.parameters.pure_ae()) {
                    vnuc = potentialmanager->vnuclear();
                } else {
                    vnuc = potentialmanager->vnuclear();
                    vnuc = vnuc + gthpseudopotential->vlocalpot();
                }

                lb.add_tree(vnuc,
                            lbcost<double, 3>(param.vnucextra() * 1.0, param.vnucextra() * 8.0), false);
                lb.add_tree(rho, lbcost<double, 3>(1.0, 8.0), true);

                FunctionDefaults<3>::redistribute(world, lb.load_balance(param.get<int>("loadbalparts")));
                END_TIMER(world, "guess loadbal");
            }

            // Diag approximate fock matrix to get initial mos
            functionT vlocal;
            if (param.nalpha() + param.nbeta() > 1) {
                START_TIMER(world);
                real_function_3d vnuc;
                if (molecule.parameters.psp_calc()) {
                    vnuc = gthpseudopotential->vlocalpot();
                } else if (molecule.parameters.pure_ae()) {
                    vnuc = potentialmanager->vnuclear();
                } else {
                    vnuc = potentialmanager->vnuclear();
                    vnuc = vnuc + gthpseudopotential->vlocalpot();
                }
                vlocal = vnuc + apply(*coulop, rho);
                END_TIMER(world, "guess Coulomb potn");
                START_TIMER(world);
                vlocal = vlocal + make_lda_potential(world, rho);
                vlocal.truncate();
                END_TIMER(world, "guess lda potn");
            } else {
                real_function_3d vnuc;
                if (molecule.parameters.psp_calc()) {
                    vnuc = gthpseudopotential->vlocalpot();
                } else if (molecule.parameters.pure_ae()) {
                    vnuc = potentialmanager->vnuclear();
                } else {
                    vnuc = potentialmanager->vnuclear();
                    vnuc = vnuc + gthpseudopotential->vlocalpot();
                }
                vlocal = vnuc;
            }
            rho.clear();
            vlocal.reconstruct();
            if (world.size() > 1) {
                START_TIMER(world);
                LoadBalanceDeux<3> lb(world);
                real_function_3d vnuc;
                if (molecule.parameters.psp_calc()) {
                    vnuc = gthpseudopotential->vlocalpot();
                } else if (molecule.parameters.pure_ae()) {
                    vnuc = potentialmanager->vnuclear();
                } else {
                    vnuc = potentialmanager->vnuclear();
                    vnuc = vnuc + gthpseudopotential->vlocalpot();
                }
                lb.add_tree(vnuc,
                            lbcost<double, 3>(param.vnucextra() * 1.0, param.vnucextra() * 8.0), false);
                for (unsigned int i = 0; i < ao.size(); ++i) {
                    lb.add_tree(ao[i], lbcost<double, 3>(1.0, 8.0), false);
                }
                FunctionDefaults<3>::redistribute(world, lb.load_balance(param.get<int>("loadbalparts")));
                END_TIMER(world, "guess loadbal");
            }
            START_TIMER(world);
            tensorT overlap = matrix_inner(world, ao, ao, true);
            END_TIMER(world, "guess overlap");
            START_TIMER(world);

            tensorT kinetic(ao.size(), ao.size());
            {
                distmatT dkinetic = kinetic_energy_matrix(world, ao);
                dkinetic.copy_to_replicated(kinetic);
            }
            END_TIMER(world, "guess Kinet potn");

            START_TIMER(world);
            reconstruct(world, ao);
            vlocal.reconstruct();
            vecfuncT vpsi;

            //debug plots:
            /*{
                    int npt=1001;
                    functionT rhotmp =
                        factoryT(world).functor(
                                                functorT(
                                                         new MolecularGuessDensityFunctor(molecule,
                                                                                          aobasis))).truncate_on_project();
                    functionT vlda=make_lda_potential(world, rhotmp);
                    functionT coul=apply(*coulop, rhotmp);
                    plot_line("vlocal.dat",npt, {0.0,0.0,-50.0}, {0.0,0.0,50.0}, vlocal);
                    plot_line("vcoul.dat",npt, {0.0,0.0,-50.0}, {0.0,0.0,50.0}, vcoul);
                    plot_line("vlda.dat",npt, {0.0,0.0,-50.0}, {0.0,0.0,50.0}, vlda);
                    plot_line("dens.dat",npt, {0.0,0.0,-50.0}, {0.0,0.0,50.0}, rhotmp);

                    if (!param.pure_ae && !param.psp_calc){
                        real_function_3d vloc_ae;
                        vloc_ae = potentialmanager->vnuclear();
                        vloc_ae.reconstruct();
                        plot_line("vlocal_ae.dat",npt, {0.0,0.0,-50.0}, {0.0,0.0,50.0}, vloc_ae);
                        real_function_3d vloc_psp;
                        vloc_psp = gthpseudopotential->vlocalpot();
                        vloc_psp.reconstruct();
                        plot_line("vlocal_psp.dat",npt, {0.0,0.0,-50.0}, {0.0,0.0,50.0}, vloc_psp);
                    }
                }*/

            //vlocal treated in psp includes psp and ae contribution so don't need separate clause for mixed psp/AE
            if (!molecule.parameters.pure_ae()) {
                double enl;
                tensorT occ = tensorT(ao.size());
                for (int i = 0; i < param.nalpha(); ++i) {
                    occ[i] = 1.0;
                }
                for (int i = param.nalpha(); size_t(i) < ao.size(); ++i) {
                    occ[i] = 0.0;
                }
                vpsi = gthpseudopotential->apply_potential(world, vlocal, ao, occ, enl);
            } else {
                vpsi = mul_sparse(world, vlocal, ao, vtol);
            }

            compress(world, vpsi);
            truncate(world, vpsi);
            compress(world, ao);
            tensorT potential = matrix_inner(world, vpsi, ao, true);
            vpsi.clear();
            tensorT fock = kinetic + potential;
            fock = 0.5 * (fock + transpose(fock));
            tensorT c, e;

            //debug printing
            /*double ep = 0.0;
                double ek = 0.0;
                for(int i = 0;i < ao.size();++i){
                    ep += potential(i, i);
                    ek += kinetic(i, i);
                    std::cout << "pot/kin " << i << "  " << potential(i,i) << "  "<< kinetic(i,i) << std::endl;
                }

                if(world.rank() == 0){
                    printf("\n              epot, ekin, efock %16.8f  %16.8f  %16.8f\n", ek, ep, ek+ep);
             */

            END_TIMER(world, "guess fock");

            START_TIMER(world);
            sygvp(world, fock, overlap, 1, c, e);
            END_TIMER(world, "guess eigen sol");
            print_meminfo(world.rank(), "guess eigen sol");

            // NAR 7/5/2013
            // commented out because it generated a lot of output
            // if(world.rank() == 0 && 0){
            //   print("initial eigenvalues");
            //   print(e);
            //   print("\n\nWSTHORNTON: initial eigenvectors");
            //   print(c);
            // }

            START_TIMER(world);
            compress(world, ao);

            unsigned int ncore = 0;
            if (molecule.parameters.core_type() != "none") {
                ncore = molecule.n_core_orb_all();
            }

            amo = transform(world, ao, c(_, Slice(ncore, ncore + param.nmo_alpha() - 1)), vtol, true);
            truncate(world, amo);
            normalize(world, amo);
            aeps = e(Slice(ncore, ncore + param.nmo_alpha() - 1));

            aocc = tensorT(param.nmo_alpha());
            for (int i = 0; i < param.nalpha(); ++i)
                aocc[i] = 1.0;

            if (world.rank() == 0 and param.print_level() > 3) print("grouping alpha orbitals into sets");
            aset = group_orbital_sets(world, aeps, aocc, param.nmo_alpha());

            if (param.nbeta() && !param.spin_restricted()) {
                bmo = transform(world, ao, c(_, Slice(ncore, ncore + param.nmo_beta() - 1)), vtol, true);
                truncate(world, bmo);
                normalize(world, bmo);
                beps = e(Slice(ncore, ncore + param.nmo_beta() - 1));
                bocc = tensorT(param.nmo_beta());
                for (int i = 0; i < param.nbeta(); ++i)
                    bocc[i] = 1.0;

                if (world.rank() == 0 and param.print_level() > 3) print("grouping beta orbitals into sets");
                bset = group_orbital_sets(world, beps, bocc, param.nmo_beta());

            }
            END_TIMER(world, "guess orbital grouping");
        }
            // If using nwchem, read in and generate intial guess here
        else {
            START_TIMER(world);

            // Construct interfact object from slymer namespace
            slymer::NWChem_Interface nwchem(param.nwfile(), std::cout);

            // For parallel runs, silencing all but 1 slymer instance
            // If print_level is too low, silence all
            if (world.rank() != 0 or param.print_level() < 4) {
                std::ostream dev_null(nullptr);
                nwchem.err = dev_null;
            }

            // Read in basis set
            nwchem.read(slymer::Properties::Basis);

            // Read in the molecular orbital coefficients, energies,
            // and occupancies
            nwchem.read(slymer::Properties::Energies | slymer::Properties::MOs | slymer::Properties::Occupancies);

            // Shift madness atoms to match nwchem atoms
            if (world.rank() == 0 && param.print_level() > 3)
                print("\nAligning atoms by moving MADNESS atoms to match NWChem atoms.");

            // Verify at least same number of atoms first
            MADNESS_ASSERT(int(nwchem.atoms.size()) == molecule.natom());

            // Get center of charge for nwchem
            std::vector<double> nw_coc(3, 0);
            double total_charge = 0.0;
            if (world.rank() == 0 && param.print_level() > 3) print("NWChem coordinates:");
            for (auto atom: nwchem.atoms) {
                int charge = symbol_to_atomic_number(atom.symbol);
                total_charge += charge;
                if (world.rank() == 0 && param.print_level() > 3)
                    print(atom.symbol, atom.position[0], atom.position[1], atom.position[2]);
                nw_coc[0] += atom.position[0] * charge;
                nw_coc[1] += atom.position[1] * charge;
                nw_coc[2] += atom.position[2] * charge;
            }
            nw_coc[0] = nw_coc[0] / total_charge;
            nw_coc[1] = nw_coc[1] / total_charge;
            nw_coc[2] = nw_coc[2] / total_charge;

            // Get center of charge for madness
            std::vector<double> mad_coc(3, 0);
            for (size_t i = 0; i < molecule.natom(); ++i) {
                const Atom& atom = molecule.get_atom(i);
                int charge = atom.atomic_number;
                mad_coc[0] += atom.x * charge;
                mad_coc[1] += atom.y * charge;
                mad_coc[2] += atom.z * charge;
            }
            mad_coc[0] = mad_coc[0] / total_charge;
            mad_coc[1] = mad_coc[1] / total_charge;
            mad_coc[2] = mad_coc[2] / total_charge;

            // Now translate MADNESS to have same coc as NWChem
            Tensor<double> translation(3);
            for (unsigned int i = 0; i < 3; i++) {
                translation[i] = mad_coc[i] - nw_coc[i];
            }
            molecule.translate(translation);

            // Now construct the rotation such that the overlap between NWChem
            // and MADNESS is maximized.
            // First need the locations in a tensor for manipulations
            Tensor<double> nw_coords(nwchem.atoms.size(), 4);
            Tensor<double> mad_coords(molecule.natom(), 4);
            for (unsigned int i = 0; i < nwchem.atoms.size(); i++) {
                nw_coords(i, 0) = nwchem.atoms[i].position[0];
                nw_coords(i, 1) = nwchem.atoms[i].position[1];
                nw_coords(i, 2) = nwchem.atoms[i].position[2];
                nw_coords(i, 3) = symbol_to_atomic_number(nwchem.atoms[i].symbol) * 1000.0;

                const Atom& atom = molecule.get_atom(i);
                mad_coords(i, 0) = atom.x;
                mad_coords(i, 1) = atom.y;
                mad_coords(i, 2) = atom.z;
                mad_coords(i, 3) = atom.atomic_number * 1000.0;
            }

            // Using polar decomp to construct rotation
            Tensor<double> q = inner(transpose(mad_coords), nw_coords);
            Tensor<double> VT(4, 4);
            Tensor<double> U(4, 4);
            Tensor<double> sigma(4);
            svd(q, U, sigma, VT);
            q = inner(U, VT);

            // And rotate
            molecule.rotate(q(Slice(0, 2), Slice(0, 2)));
            if (world.rank() == 0 && param.print_level() > 3) print("New MADNESS coordinates:");
            for (size_t i = 0; i < molecule.natom(); ++i) {
                const Atom& atom = molecule.get_atom(i);
                if (world.rank() == 0 && param.print_level() > 3)
                    print(atomic_number_to_symbol(atom.atomic_number), atom.x, atom.y, atom.z);
            }

            // Construct nuclear potential
            make_nuclear_potential(world);
            real_function_3d vnuc = potentialmanager->vnuclear();

            // Pull out occupation numbers
            // NWChem orders occupied orbitals to be first
            aocc = tensorT(param.nalpha());
            for (int i = 0; i < param.nalpha(); i++) {
                // NWChem stores closed shell calculations
                // as the alpha orbital set with occupation 2.
                // Verifying no fractional occupations.
                MADNESS_ASSERT(nwchem.occupancies[i] == 2.0 or nwchem.occupancies[i] == 1.0);

                // Madness instead stores 2 identical sets
                // (alpha and beta) with occupation 1
                aocc[i] = 1.0;
            }

            // Pull out energies
            aeps = tensorT(param.nalpha());
            for (int i = 0; i < param.nalpha(); i++) {
                aeps[i] = nwchem.energies[i];
            }

            // Create the orbitals as madness functions
            // Just create the vector of atomic orbitals
            // and use the vector of MO coefficients and
            // the transform function, then take only
            // the occupied orbitals.
            if (world.rank() == 0 && param.print_level() > 3)
                print("\nCreating MADNESS functions from the NWChem orbitals.");

            // Cast the 'basis_set' into a gaussian basis
            // and iterate over it
            vector_real_function_3d temp1;
            int i = 0;
            for (auto basis: slymer::cast_basis<slymer::GaussianFunction>(nwchem.basis_set)) {
                // Get the center of gaussian as its special point
                std::vector<coord_3d> centers;
                coord_3d r;
                r[0] = basis.get().center[0];
                r[1] = basis.get().center[1];
                r[2] = basis.get().center[2];
                centers.push_back(r);

                // Now make the function
                temp1.push_back(factoryT(world).functor(functorT(new slymer::Gaussian_Functor(basis.get(), centers))));
                if (world.rank() == 0 and i % 10 == 0 and i != 0 && param.print_level() > 3)
                    print("Created", i, "functions.");
                i++;
            }
            if (world.rank() == 0 && param.print_level() > 3) print("Finished creating", temp1.size(), "functions.");

            // Transform ao's now
            vector_real_function_3d temp = transform(world, temp1, nwchem.MOs, vtol, true);

            // Now save all aos and only the occupied amo
            for (unsigned int i = 0; i < temp1.size(); i++) {
                // Save all AOs
                ao.push_back(copy(temp1[i]));

                // Only save occupied AMOs
                if (nwchem.occupancies[i] > 0) {
                    amo.push_back(copy(temp[i]));
                }
            }

            // Clean up
            truncate(world, amo);
            normalize(world, amo);

            if (world.rank() == 0 && param.print_level() > 3) print("\ngrouping alpha orbitals into sets");
            aset = group_orbital_sets(world, aeps, aocc, param.nmo_alpha());

            // Now for betas
            if (param.nbeta() && !param.spin_restricted()) {

                // Pull out occupation numbers
                // NWChem orders occupied orbitals to be first
                bocc = tensorT(param.nbeta());
                for (int i = 0; i < param.nbeta(); i++) {
                    MADNESS_ASSERT(nwchem.beta_occupancies[i] == 1.0);
                    bocc[i] = 1.0;
                }

                // Pull out energies
                beps = tensorT(param.nbeta());
                for (int i = 0; i < param.nbeta(); i++) {
                    beps[i] = nwchem.beta_energies[i];
                }

                // Transform ao's now
                temp = transform(world, temp1, nwchem.beta_MOs, vtol, true);

                // Now only take the occupied bmo
                for (unsigned int i = 0; i < temp1.size(); i++) {
                    if (nwchem.beta_occupancies[i] > 0) {
                        bmo.push_back(copy(temp[i]));
                    }
                }

                // Clean up
                truncate(world, bmo);
                normalize(world, bmo);

                if (world.rank() == 0 && param.print_level() > 3) print("\ngrouping beta orbitals into sets");
                bset = group_orbital_sets(world, beps, bocc, param.nmo_beta());
            }

            END_TIMER(world, "read nwchem file");
        }

    }
}

/// group orbitals into sets of similar orbital energies for localization

/// @param[in]	eps	orbital energies
/// @param[in]	occ	occupation numbers
/// @param[in]	nmo number of MOs for the given spin
/// @return		vector of length nmo with the set index for each MO
std::vector<int> SCF::group_orbital_sets(World& world, const tensorT& eps,
                                         const tensorT& occ, const int nmo) const {
    PROFILE_MEMBER_FUNC(SCF);

    std::vector<int> set = std::vector<int>(static_cast<size_t>(nmo), 0);
    for (int i = 1; i < nmo; ++i) {
        set[i] = set[i - 1];
        // Only the new/boys localizers can tolerate not separating out the core orbitals
        if (param.localize_pm() && (eps[i] - eps[i - 1] > 1.5 || occ[i] != 1.0)) ++(set[i]);
    }

    // pretty print out
    int lo = 0;
    int iset = 0;
    for (size_t i = 0; i < set.size(); ++i) {
        if (iset != set[i]) {
            if (world.rank() == 0 and (param.print_level() > 3)) print("set ", iset++, "  ", lo, " - ", i - 1);
            lo = i;
        }
    }
    if (world.rank() == 0 and (param.print_level() > 3)) print("set ", iset, "  ", lo, " - ", nmo - 1);
    return set;
}


void SCF::initial_load_bal(World& world) {
    PROFILE_MEMBER_FUNC(SCF);
    LoadBalanceDeux<3> lb(world);
    real_function_3d vnuc;
    if (molecule.parameters.psp_calc()) {
        vnuc = gthpseudopotential->vlocalpot();
    } else if (molecule.parameters.pure_ae()) {
        vnuc = potentialmanager->vnuclear();
    } else {
        vnuc = potentialmanager->vnuclear();
        vnuc = vnuc + gthpseudopotential->vlocalpot();
    }
    lb.add_tree(vnuc, lbcost<double, 3>(param.vnucextra() * 1.0, param.vnucextra() * 8.0));

    FunctionDefaults<3>::redistribute(world, lb.load_balance(param.loadbalparts()));
}

functionT SCF::make_density(World& world, const tensorT& occ,
                            const vecfuncT& v) const {
    PROFILE_MEMBER_FUNC(SCF);
    vecfuncT vsq = square(world, v);
    compress(world, vsq);
    functionT rho = factoryT(world);
    rho.compress();
    for (unsigned int i = 0; i < vsq.size(); ++i) {
        if (occ[i]) rho.gaxpy(1.0, vsq[i], occ[i], false);
    }
    world.gop.fence();
    vsq.clear();
    return rho;
}

functionT SCF::make_density(World& world, const tensorT& occ,
                            const cvecfuncT& v) {
    PROFILE_MEMBER_FUNC(SCF);
    reconstruct(world, v); // For max parallelism
    std::vector<functionT> vsq(v.size());
    for (unsigned int i = 0; i < v.size(); i++) {
        vsq[i] = abssq(v[i], false);
    }
    world.gop.fence();

    compress(world, vsq); // since will be using gaxpy for accumulation
    functionT rho = factoryT(world);
    rho.compress();

    for (unsigned int i = 0; i < vsq.size(); ++i) {
        if (occ[i])
            rho.gaxpy(1.0, vsq[i], occ[i], false);

    }
    world.gop.fence();
    vsq.clear();
    rho.truncate();

    return rho;
}

std::vector<poperatorT> SCF::make_bsh_operators(World& world, const tensorT& evals) const {
    PROFILE_MEMBER_FUNC(SCF);
    int nmo = evals.dim(0);
    std::vector<poperatorT> ops(nmo);
    double tol = FunctionDefaults<3>::get_thresh();
    for (int i = 0; i < nmo; ++i) {
        double eps = evals(i);
        if (eps > 0) {
            if (world.rank() == 0 and (param.print_level() > 3)) {
                print("bsh: warning: positive eigenvalue", i, eps);
            }
            eps = -0.1;
        }

        ops[i] = poperatorT(
                BSHOperatorPtr3D(world, sqrt(-2.0 * eps), param.lo(), tol));
    }

    return ops;
}


// Used only for initial guess that is always spin-restricted LDA
functionT SCF::make_lda_potential(World& world, const functionT& arho) {
    PROFILE_MEMBER_FUNC(SCF);
    functionT vlda = copy(arho);
    vlda.reconstruct();
    vlda.unaryop(xc_lda_potential());
    return vlda;
}

vecfuncT SCF::apply_potential(World& world, const tensorT& occ,
                              const vecfuncT& amo,
                              const functionT& vlocal, double& exc, double& enl, int ispin) {
    PROFILE_MEMBER_FUNC(SCF);
    functionT vloc = copy(vlocal);
    exc = 0.0;
    enl = 0.0;

    vecfuncT Vpsi;
    print_meminfo(world.rank(), "V*psi");
    if (xc.hf_exchange_coefficient()) {
        START_TIMER(world);
        //            vecfuncT Kamo = apply_hf_exchange(world, occ, amo, amo);
        Exchange<double, 3> K(world, this, ispin);
	if (param.hfexalg()=="multiworld") {
	  //if (world.rank() == 0) print("selecting exchange multi world");
	  K.set_algorithm(Exchange<double,3>::Algorithm::multiworld_efficient);
	}
	else if (param.hfexalg()=="multiworld_row") {
	  //if (world.rank() == 0) print("selecting exchange multi world row");
	  K.set_algorithm(Exchange<double,3>::Algorithm::multiworld_efficient_row);
	}
	else if (param.hfexalg()=="largemem") {
	  //if (world.rank() == 0) print("selecting exchange large memory");
	  K.set_algorithm(Exchange<double,3>::Algorithm::large_memory);
	}
	else if (param.hfexalg()=="smallmem") {
	  //if (world.rank() == 0) print("selecting exchange small memory");
	  K.set_algorithm(Exchange<double,3>::Algorithm::small_memory);
	}
	
        K.set_symmetric(true).set_printlevel(param.print_level());
        vecfuncT Kamo = K(amo);
        tensorT excv = inner(world, Kamo, amo);
        double exchf = 0.0;
        for (unsigned long i = 0; i < amo.size(); ++i) {
            exchf -= 0.5 * excv[i] * occ[i];
        }
        if (!xc.is_spin_polarized()) exchf *= 2.0;
        Vpsi=-xc.hf_exchange_coefficient()* Kamo;
        Kamo.clear();
        END_TIMER(world, "HF exchange");
        exc = exchf * xc.hf_exchange_coefficient() + exc;
    }
    else {
      Vpsi = zero_functions<double,3>(world, amo.size());	  
    }

    // compute the local DFT potential for the MOs
    if (xc.is_dft() && !(xc.hf_exchange_coefficient() == 1.0)) { //??RJH?? Won't this incorrectly exclude hybrid DFT with coeff=1.0?
        START_TIMER(world);

        XCOperator<double, 3> xcoperator(world, this, ispin, param.dft_deriv());
        if (ispin == 0) exc = xcoperator.compute_xc_energy();
        vloc += xcoperator.make_xc_potential();

        END_TIMER(world, "DFT potential");
    }

    vloc.truncate();
    
    // need to come back to this for psp - when is this used?
    // RJH commented this out since it seems to never do anything useful ... if pure all-electron there is not a non-local part
    // if (molecule.parameters.pure_ae()) {
    //     potentialmanager->apply_nonlocal_potential(world, amo, Vpsi);
    // }

    START_TIMER(world);
    if (!molecule.parameters.pure_ae()) {
        gaxpy(world, 1.0, Vpsi, 1.0, gthpseudopotential->apply_potential(world, vloc, amo, occ, enl));
    } else {
        gaxpy(world, 1.0, Vpsi, 1.0, mul_sparse(world, vloc, amo, vtol));
    }

    END_TIMER(world, "V*psi");

    START_TIMER(world);
    truncate(world, Vpsi);
    END_TIMER(world, "Truncate Vpsi");
    print_meminfo(world.rank(), "Truncate Vpsi");
    world.gop.fence();
    return Vpsi;
}

tensorT SCF::derivatives(World& world, const functionT& rho) const {
    PROFILE_MEMBER_FUNC(SCF);
    START_TIMER(world);

    vecfuncT dv(molecule.natom() * 3);
    vecfuncT du = zero_functions<double, 3>(world, molecule.natom() * 3);
    tensorT rc(molecule.natom() * 3);
    for (size_t atom = 0; atom < molecule.natom(); ++atom) {
        for (int axis = 0; axis < 3; ++axis) {
            functorT func(new MolecularDerivativeFunctor(molecule, atom, axis));
            dv[atom * 3 + axis] =
                    functionT(
                            factoryT(world).functor(func).nofence().truncate_on_project().truncate_mode(0));
            if (molecule.parameters.core_type() != "none"
                && molecule.is_potential_defined_atom(atom)) {
                // core potential contribution
                func = functorT(
                        new CorePotentialDerivativeFunctor(molecule, atom,
                                                           axis));
                du[atom * 3 + axis] = functionT(
                        factoryT(world).functor(func).truncate_on_project());

                // core projector contribution
                rc[atom * 3 + axis] =
                        potentialmanager->core_projector_derivative(world, amo,
                                                                    aocc, atom, axis);
                if (!param.spin_restricted()) {
                    if (param.nbeta())
                        rc[atom * 3 + axis] +=
                                potentialmanager->core_projector_derivative(
                                        world, bmo, bocc, atom, axis);
                } else {
                    rc[atom * 3 + axis] *= 2 * 2;
                    // because of 2 electrons in each valence orbital bra+ket
                }
            }
        }
    }

    world.gop.fence();
    tensorT r = inner(world, rho, dv);
    world.gop.fence();
    tensorT ru = inner(world, rho, du);
    dv.clear();
    du.clear();
    world.gop.fence();
    tensorT ra(r.size());
    for (size_t atom = 0; atom < molecule.natom(); ++atom) {
        for (int axis = 0; axis < 3; ++axis) {
            ra[atom * 3 + axis] = molecule.nuclear_repulsion_derivative(atom,
                                                                        axis);
        }
    }
    //if (world.rank() == 0) print("derivatives:\n", r, ru, rc, ra);
    r += ra + ru + rc;
    END_TIMER(world, "derivatives");

    if (world.rank() == 0 and (param.print_level() > 1)) {
        print("\n Derivatives (a.u.)\n -----------\n");
        print(
                "  atom        x            y            z          dE/dx        dE/dy        dE/dz");
        print(
                " ------ ------------ ------------ ------------ ------------ ------------ ------------");
        for (size_t i = 0; i < molecule.natom(); ++i) {
            const Atom& atom = molecule.get_atom(i);
            printf(" %5d %12.6f %12.6f %12.6f %12.6f %12.6f %12.6f\n", int(i),
                   atom.x, atom.y, atom.z, r[i * 3 + 0], r[i * 3 + 1],
                   r[i * 3 + 2]);
        }
    }
    return r;
}

tensorT SCF::dipole(World& world, const functionT& rho) const {
    PROFILE_MEMBER_FUNC(SCF);
    START_TIMER(world);
    tensorT mu(3);

    for (unsigned int axis = 0; axis < 3; ++axis) {
        std::vector<int> x(3ul, 0);
        x[axis] = true;
        functionT dipolefunc = factoryT(world)
                .functor(functorT(new MomentFunctor(x)));
        mu[axis] = -dipolefunc.inner(rho);
        mu[axis] += molecule.nuclear_dipole(axis);
    }

    if (world.rank() == 0 and (param.print_level() > 2)) {
        print("\n Dipole Moment (a.u.)\n -----------\n");
        print("     x: ", mu[0]);
        print("     y: ", mu[1]);
        print("     z: ", mu[2]);
        print(" Total Dipole Moment: ", mu.normf(), "\n");
    }
    END_TIMER(world, "dipole");

    return mu;
}

void SCF::vector_stats(const std::vector<double>& v, double& rms,
                       double& maxabsval) const {
    PROFILE_MEMBER_FUNC(SCF);
    rms = 0.0;
    maxabsval = v[0];
    for (unsigned int i = 0; i < v.size(); ++i) {
        rms += v[i] * v[i];
        maxabsval = std::max<double>(maxabsval, std::abs(v[i]));
    }
    rms = sqrt(rms / v.size());
}

vecfuncT SCF::compute_residual(World& world, tensorT& occ, tensorT& fock,
                               const vecfuncT& psi, vecfuncT& Vpsi, double& err) {

    START_TIMER(world);
    PROFILE_MEMBER_FUNC(SCF);
    double trantol = vtol / std::min(30.0, double(psi.size()));
    int nmo = psi.size();

    tensorT eps(nmo);
    for (int i = 0; i < nmo; ++i) {
        eps(i) = std::min(-0.05, fock(i, i));
        fock(i, i) -= eps(i);
    }
    vecfuncT fpsi = transform(world, psi, fock, trantol, true);

    for (int i = 0; i < nmo; ++i) { // Undo the damage
        fock(i, i) += eps(i);
    }

    gaxpy(world, 1.0, Vpsi, -1.0, fpsi);
    fpsi.clear();
    std::vector<double> fac(nmo, -2.0);
    scale(world, Vpsi, fac);
    std::vector<poperatorT> ops = make_bsh_operators(world, eps);
    set_thresh(world, Vpsi, FunctionDefaults<3>::get_thresh());
    END_TIMER(world, "Compute residual stuff");

    START_TIMER(world);
    vecfuncT new_psi = apply(world, ops, Vpsi);
    END_TIMER(world, "Apply BSH");
    ops.clear();
    Vpsi.clear();
    world.gop.fence();

    // Thought it was a bad idea to truncate *before* computing the residual
    // but simple tests suggest otherwise ... no more iterations and
    // reduced iteration time from truncating.
    START_TIMER(world);
    truncate(world, new_psi);
    END_TIMER(world, "Truncate new psi");

    START_TIMER(world);
    vecfuncT r = sub(world, psi, new_psi);
    std::vector<double> rnorm = norm2s(world, r);
    if (world.rank() == 0 and (param.print_level() > 1))
        print("residuals", rnorm);
    double rms, maxval;
    vector_stats(rnorm, rms, maxval);
    err = maxval;
    if (world.rank() == 0 and (param.print_level() > 1))
        print("BSH residual: rms", rms, "   max", maxval);
    END_TIMER(world, "BSH residual");
    return r;
}

tensorT SCF::make_fock_matrix(World& world, const vecfuncT& psi,
                              const vecfuncT& Vpsi, const tensorT& occ, double& ekinetic) const {
    PROFILE_MEMBER_FUNC(SCF);
    START_TIMER(world);
    tensorT pe = matrix_inner(world, Vpsi, psi, true);
    END_TIMER(world, "PE matrix");

    std::shared_ptr<WorldDCPmapInterface<Key<3> > > oldpmap = FunctionDefaults<3>::get_pmap();
    vecfuncT psicopy = psi; // Functions are shallow copy so this is lightweight
    if (world.size() > 1) {
        START_TIMER(world);
        LoadBalanceDeux<3> lb(world);
        for (unsigned int i = 0; i < psi.size(); ++i) {
            lb.add_tree(psi[i], lbcost<double, 3>(1.0, 8.0), false);
        }
        world.gop.fence();
        END_TIMER(world, "KE compute loadbal");

        START_TIMER(world);
        std::shared_ptr<WorldDCPmapInterface<Key<3> > > newpmap = lb.load_balance(param.loadbalparts());
        FunctionDefaults<3>::set_pmap(newpmap);

        world.gop.fence();
        for (unsigned int i = 0; i < psi.size(); ++i) psicopy[i] = copy(psi[i], newpmap, false);
        world.gop.fence();
        END_TIMER(world, "KE redist");
    }

    START_TIMER(world);
    tensorT ke(psi.size(), psi.size());
    {
        distmatT k = kinetic_energy_matrix(world, psicopy);
        k.copy_to_replicated(ke); // !!!!!!!! ugh
    }
    END_TIMER(world, "KE matrix");

    psicopy.clear();
    if (world.size() > 1) {
        FunctionDefaults<3>::set_pmap(oldpmap); // ! DON'T FORGET !
    }

    START_TIMER(world);
    int nocc = occ.size();
    ekinetic = 0.0;
    for (int i = 0; i < nocc; ++i) {
        ekinetic += occ[i] * ke(i, i);
    }
    ke += pe;
    pe = tensorT();
    ke.gaxpy(0.5, transpose(ke), 0.5);
    END_TIMER(world, "Make fock matrix rest");
    return ke;
}

/// Compute the two-electron integrals over the provided set of orbitals

/// Returned is a *replicated* tensor of \f$(ij|kl)\f$ with \f$i>=j\f$
/// and \f$k>=l\f$.  The symmetry \f$(ij|kl)=(kl|ij)\f$ is enforced.
Tensor<double> SCF::twoint(World& world, const vecfuncT& psi) const {
    PROFILE_MEMBER_FUNC(SCF);
    double tol = FunctionDefaults<3>::get_thresh(); /// Important this is consistent with Coulomb
    reconstruct(world, psi);
    norm_tree(world, psi);

    // Efficient version would use mul_sparse vector interface
    vecfuncT pairs;
    for (unsigned int i = 0; i < psi.size(); ++i) {
        for (unsigned int j = 0; j <= i; ++j) {
            pairs.push_back(mul_sparse(psi[i], psi[j], tol, false));
        }
    }

    world.gop.fence();
    truncate(world, pairs);
    vecfuncT Vpairs = apply(world, *coulop, pairs);

    return matrix_inner(world, pairs, Vpairs, true);
}

/// compute the unitary transformation that diagonalizes the fock matrix

/// @param[in]  world   the world
/// @param[in]  overlap the overlap matrix of the orbitals
/// @param[inout]       fock    the fock matrix; diagonal upon exit
/// @param[out] evals   the orbital energies
/// @param[in]  occ     the occupation numbers
/// @param[in]  thresh_degenerate       threshold for orbitals being degenerate
/// @return             the unitary matrix U: U^T F U = evals
tensorT SCF::get_fock_transformation(World& world, const tensorT& overlap,
                                     tensorT& fock, tensorT& evals, const tensorT& occ,
                                     const double thresh_degenerate) const {
    PROFILE_MEMBER_FUNC(SCF);

    tensorT U;
    sygvp(world, fock, overlap, 1, U, evals);

    Localizer::undo_reordering(U, occ, evals);
    Localizer::undo_degenerate_rotations(U, evals, thresh_degenerate);

    world.gop.broadcast(U.ptr(), U.size(), 0);
    world.gop.broadcast(evals.ptr(), evals.size(), 0);

    fock = 0;
    for (unsigned int i = 0; i < evals.size(); ++i)
        fock(i, i) = evals(i);
    return U;
}

/// diagonalize the fock matrix, taking care of degenerate states

/// Vpsi is passed in to make sure orbitals and Vpsi are in phase
/// @param[in]  world   the world
/// @param[inout]       fock    the fock matrix (diagonal upon exit)
/// @param[inout]       psi             the orbitals
/// @param[inout]       Vpsi    the orbital times the potential
/// @param[out] evals   the orbital energies
/// @param[in]  occ             occupation numbers
/// @param[in]  thresh  threshold for rotation and truncation
/// @return             the unitary matrix U: U^T F U = evals
tensorT SCF::diag_fock_matrix(World& world, tensorT& fock, vecfuncT& psi,
                              vecfuncT& Vpsi, tensorT& evals, const tensorT& occ,
                              const double thresh) const {
    PROFILE_MEMBER_FUNC(SCF);

    START_TIMER(world);
    // compute the unitary transformation matrix U that diagonalizes
    // the fock matrix
    tensorT overlap = matrix_inner(world, psi, psi, true);
    tensorT U = get_fock_transformation(world, overlap, fock, evals, occ,
                                        thresh);

    //eliminate mixing between occ and unocc
    int nmo = U.dim(0);
    for (int i = 0; i < param.nalpha(); ++i) {
        //make virt orthog to occ without changing occ states
        for (int j = param.nalpha(); j < nmo; ++j) {
            U(j, i) = 0.0;
        }
    }

    // transform the orbitals and the orbitals times the potential
    Vpsi = transform(world, Vpsi, U, vtol / std::min(30.0, double(psi.size())), false);
    psi = transform(world, psi, U, FunctionDefaults<3>::get_thresh() / std::min(30.0, double(psi.size())),
                    true);
    truncate(world, Vpsi, vtol, false);
    truncate(world, psi);
    normalize(world, psi);

    END_TIMER(world, "Diagonalization rest");
    return U;
}

void SCF::loadbal(World& world, functionT& arho, functionT& brho,
                  functionT& arho_old, functionT& brho_old, subspaceT& subspace) {
    if (world.size() == 1)
        return;

    LoadBalanceDeux<3> lb(world);
    real_function_3d vnuc;
    if (molecule.parameters.psp_calc()) {
        vnuc = gthpseudopotential->vlocalpot();
    } else if (molecule.parameters.pure_ae()) {
        vnuc = potentialmanager->vnuclear();
    } else {
        vnuc = potentialmanager->vnuclear();
        vnuc = vnuc + gthpseudopotential->vlocalpot();
    }
    lb.add_tree(vnuc, lbcost<double, 3>(param.vnucextra() * 1.0, param.vnucextra() * 8.0),
                false);
    lb.add_tree(arho, lbcost<double, 3>(1.0, 8.0), false);
    for (unsigned int i = 0; i < amo.size(); ++i) {
        lb.add_tree(amo[i], lbcost<double, 3>(1.0, 8.0), false);
    }
    if (param.nbeta() && !param.spin_restricted()) {
        lb.add_tree(brho, lbcost<double, 3>(1.0, 8.0), false);
        for (unsigned int i = 0; i < bmo.size(); ++i) {
            lb.add_tree(bmo[i], lbcost<double, 3>(1.0, 8.0), false);
        }
    }
    world.gop.fence();

    FunctionDefaults<3>::redistribute(world, lb.load_balance(
            param.loadbalparts())); // 6.0 needs retuning after param.vnucextra

    world.gop.fence();
}

void SCF::rotate_subspace(World& world, const tensorT& U, subspaceT& subspace,
                          int lo, int nfunc, double trantol) const {
    PROFILE_MEMBER_FUNC(SCF);
    for (unsigned int iter = 0; iter < subspace.size(); ++iter) {
        vecfuncT& v = subspace[iter].first;
        vecfuncT& r = subspace[iter].second;
        vecfuncT vnew = transform(world, vecfuncT(&v[lo], &v[lo + nfunc]), U, trantol, false);
        vecfuncT rnew = transform(world, vecfuncT(&r[lo], &r[lo + nfunc]), U, trantol, false);
        world.gop.fence();
        for (int i = 0; i < nfunc; i++) {
            v[i] = vnew[i];
            r[i] = rnew[i];
        }
    }
    world.gop.fence();
}

void SCF::rotate_subspace(World& world, const distmatT& dUT, subspaceT& subspace,
                          int lo, int nfunc, double trantol) const {
    PROFILE_MEMBER_FUNC(SCF);
    for (unsigned int iter = 0; iter < subspace.size(); ++iter) {
        vecfuncT& v = subspace[iter].first;
        vecfuncT& r = subspace[iter].second;
        vecfuncT vnew = transform(world, vecfuncT(&v[lo], &v[lo + nfunc]), dUT, false);
        vecfuncT rnew = transform(world, vecfuncT(&r[lo], &r[lo + nfunc]), dUT, false);
        world.gop.fence();
        for (int i = 0; i < nfunc; i++) {
            v[i] = vnew[i];
            r[i] = rnew[i];
        }
    }
    world.gop.fence();
}

void SCF::update_subspace(World& world, vecfuncT& Vpsia, vecfuncT& Vpsib,
                          tensorT& focka, tensorT& fockb, subspaceT& subspace, tensorT& Q,
                          double& bsh_residual, double& update_residual) {
    PROFILE_MEMBER_FUNC(SCF);
    double aerr = 0.0, berr = 0.0;
    vecfuncT vm = amo;

    // Orbitals with occ!=1.0 exactly must be solved for as eigenfunctions
    // so zero out off diagonal lagrange multipliers
    for (int i = 0; i < param.nmo_alpha(); i++) {
        if (aocc[i] != 1.0) {
            double tmp = focka(i, i);
            focka(i, _) = 0.0;
            focka(_, i) = 0.0;
            focka(i, i) = tmp;
        }
    }

    vecfuncT rm = compute_residual(world, aocc, focka, amo, Vpsia, aerr);
    if (param.nbeta() != 0 && !param.spin_restricted()) {
        for (int i = 0; i < param.nmo_beta(); i++) {
            if (bocc[i] != 1.0) {
                double tmp = fockb(i, i);
                fockb(i, _) = 0.0;
                fockb(_, i) = 0.0;
                fockb(i, i) = tmp;
            }
        }

        vecfuncT br = compute_residual(world, bocc, fockb, bmo, Vpsib, berr);
        vm.insert(vm.end(), bmo.begin(), bmo.end());
        rm.insert(rm.end(), br.begin(), br.end());
    }

    START_TIMER(world);
    bsh_residual = std::max(aerr, berr);
    world.gop.broadcast(bsh_residual, 0);
    compress(world, vm, false);
    compress(world, rm, false);
    world.gop.fence();

    restart:
    subspace.push_back(pairvecfuncT(vm, rm));
    int m = subspace.size();
    tensorT ms(m);
    tensorT sm(m);
    for (int s = 0; s < m; ++s) {
        const vecfuncT& vs = subspace[s].first;
        const vecfuncT& rs = subspace[s].second;
        for (unsigned int i = 0; i < vm.size(); ++i) {
            ms[s] += vm[i].inner_local(rs[i]);
            sm[s] += vs[i].inner_local(rm[i]);
        }
    }

    world.gop.sum(ms.ptr(), m);
    world.gop.sum(sm.ptr(), m);
    tensorT newQ(m, m);
    if (m > 1)
        newQ(Slice(0, -2), Slice(0, -2)) = Q;

    newQ(m - 1, _) = ms;
    newQ(_, m - 1) = sm;
    Q = newQ;
    //if (world.rank() == 0) { print("kain Q"); print(Q); }
    tensorT c;
    //if (world.rank() == 0) {
    double rcond = 1e-12;
    while (1) {
        c = KAIN(Q, rcond);
        if (world.rank() == 0 and (param.print_level() > 3)) print("kain c:", c);
        //if (std::abs(c[m - 1]) < 5.0) { // was 3
        if (c.absmax() < 3.0) { // was 3
            break;
        } else if (rcond < 0.01) {
            if (world.rank() == 0 and (param.print_level() > 3))
                print("Increasing subspace singular value threshold ", c[m - 1], rcond);
            rcond *= 100;
        } else {
            //print("Forcing full step due to subspace malfunction");
            // c = 0.0;
            // c[m - 1] = 1.0;
            // break;
            if (world.rank() == 0 and (param.print_level() > 3)) print("Restarting KAIN due to subspace malfunction");
            Q = tensorT();
            subspace.clear();
            goto restart; // fortran hat on ...
        }
    }
    //}
    END_TIMER(world, "Update subspace stuff");

    world.gop.broadcast_serializable(c, 0); // make sure everyone has same data
    if (world.rank() == 0 and (param.print_level() > 3)) {
        print("Subspace solution", c);
    }
    START_TIMER(world);
    vecfuncT amo_new = zero_functions_compressed<double, 3>(world, amo.size(), false);
    vecfuncT bmo_new = zero_functions_compressed<double, 3>(world, bmo.size(), false);
    world.gop.fence();
    for (unsigned int m = 0; m < subspace.size(); ++m) {
        const vecfuncT& vm = subspace[m].first;
        const vecfuncT& rm = subspace[m].second;
        const vecfuncT vma(vm.begin(), vm.begin() + amo.size());
        const vecfuncT rma(rm.begin(), rm.begin() + amo.size());
        const vecfuncT vmb(vm.end() - bmo.size(), vm.end());
        const vecfuncT rmb(rm.end() - bmo.size(), rm.end());
        gaxpy(world, 1.0, amo_new, c(m), vma, false);
        gaxpy(world, 1.0, amo_new, -c(m), rma, false);
        gaxpy(world, 1.0, bmo_new, c(m), vmb, false);
        gaxpy(world, 1.0, bmo_new, -c(m), rmb, false);
    }
    world.gop.fence();
    END_TIMER(world, "Subspace transform");
    if (param.maxsub() <= 1) {
        subspace.clear();
    } else if (subspace.size() == size_t(param.maxsub())) {
        subspace.erase(subspace.begin());
        Q = Q(Slice(1, -1), Slice(1, -1));
    }

    do_step_restriction(world, amo, amo_new, "alpha");
    orthonormalize(world, amo_new, param.nalpha());
    amo = amo_new;

    if (!param.spin_restricted() && param.nbeta() != 0) {
        do_step_restriction(world, bmo, bmo_new, "beta");
        orthonormalize(world, bmo_new, param.nbeta());
        bmo = bmo_new;
    } else if (param.nbeta()>0) {
        bmo = amo;
    }
}

/// perform step restriction following the KAIN solver

/// Limit maximum step size to make convergence more robust
/// @param[in]          world   the world
/// @param[in]          mo              vector of orbitals from previous iteration
/// @param[inout]       new_mo  vector of orbitals from the KAIN solver
/// @param[in]          spin    "alpha" or "beta" for user information
/// @return                     max residual
double SCF::do_step_restriction(World& world, const vecfuncT& mo, vecfuncT& mo_new,
                                std::string spin) const {
    PROFILE_MEMBER_FUNC(SCF);
    std::vector<double> anorm = norm2s(world, sub(world, mo, mo_new));
    int nres = 0;
    for (unsigned int i = 0; i < mo.size(); ++i) {
        if (anorm[i] > param.maxrotn()) {
            double s = param.maxrotn() / anorm[i];
            ++nres;
            if (world.rank() == 0) {
                if (nres == 1 and (param.print_level() > 2))
                    printf("  restricting step for %s orbitals:", spin.c_str());
                printf(" %d", i);
            }
            mo_new[i].gaxpy(s, mo[i], 1.0 - s, false);
        }
    }
    if (nres > 0 && world.rank() == 0 and (param.print_level() > 1))
        printf("\n");

    world.gop.fence();
    double rms, maxval;
    vector_stats(anorm, rms, maxval);
    if (world.rank() == 0 and (param.print_level() > 2))
        print("Norm of vector changes", spin, ": rms", rms, "   max", maxval);
    return maxval;
}

/// orthonormalize the vectors (symmetric in occupied spaced, gramm-schmidt for virt to occ)

/// @param[in]          world   the world
/// @param[inout]       amo_new the vectors to be orthonormalized
void SCF::orthonormalize(World& world, vecfuncT& amo_new, int nocc) const {
    PROFILE_MEMBER_FUNC(SCF);
    START_TIMER(world);
    double trantol = vtol / std::min(30.0, double(amo_new.size()));
    normalize(world, amo_new);
    double maxq;
    do {
        tensorT Q = Q2(matrix_inner(world, amo_new, amo_new)); // Q3(matrix_inner(world, amo_new, amo_new))
        maxq = 0.0;
        for (int j = 1; j < Q.dim(0); j++)
            for (int i = 0; i < j; i++)
                maxq = std::max(std::abs(Q(j, i)), maxq);

        Q.screen(trantol); // Is this really needed? Just for speed.

        //make virt orthog to occ without changing occ states --- ASSUMES symmetric form for Q2
        for (int j = nocc; j < Q.dim(0); ++j) {
            for (int i = 0; i < nocc; ++i) {
                Q(j, i) = 0.0;
                Q(i, j) *= 2.0;
            }
        }

        amo_new = transform(world, amo_new,
                            Q, trantol, true);
        truncate(world, amo_new);
        if (world.rank() == 0 and (param.print_level() > 3)) print("ORTHOG2a: maxq trantol", maxq, trantol);
        //print(Q);

    } while (maxq > 0.01);
    normalize(world, amo_new);

    END_TIMER(world, "Orthonormalize");

}

/// orthonormalize the vectors ignoring occupied/virtual distinctions

/// @param[in]          world   the world
/// @param[inout]       amo_new the vectors to be orthonormalized
void SCF::orthonormalize(World& world, vecfuncT& amo_new) const {
    PROFILE_MEMBER_FUNC(SCF);
    START_TIMER(world);
    double trantol = vtol / std::min(30.0, double(amo.size()));
    normalize(world, amo_new);
    double maxq;
    do {
        tensorT Q = Q2(matrix_inner(world, amo_new, amo_new)); // Q3(matrix_inner(world, amo_new, amo_new))
        maxq = 0.0;
        for (int j = 1; j < Q.dim(0); j++)
            for (int i = 0; i < j; i++)
                maxq = std::max(std::abs(Q(j, i)), maxq);

        //Q.screen(trantol); // ???? Is this really needed?
        amo_new = transform(world, amo_new,
                            Q, trantol, true);
        truncate(world, amo_new);
        if (world.rank() == 0 and (param.print_level() > 3)) print("ORTHOG2b: maxq trantol", maxq, trantol);
        //print(Q);

    } while (maxq > 0.01);
    normalize(world, amo_new);
    END_TIMER(world, "Orthonormalize");
}


complex_functionT APPLY(const complex_operatorT *q1d,
                        const complex_functionT& psi) {
    complex_functionT r = psi; // Shallow copy violates constness !!!!!!!!!!!!!!!!!
    coordT lo, hi;
    lo[2] = -10;
    hi[2] = +10;

    r.reconstruct();
    r.broaden();
    r.broaden();
    r.broaden();
    r.broaden();
    r = apply_1d_realspace_push(*q1d, r, 2);
    r.sum_down();
    r = apply_1d_realspace_push(*q1d, r, 1);
    r.sum_down();
    r = apply_1d_realspace_push(*q1d, r, 0);
    r.sum_down();

    return r;
}


// For given protocol, solve the DFT/HF/response equations
void SCF::solve(World& world) {
    PROFILE_MEMBER_FUNC(SCF);
    functionT arho_old, brho_old;
    const double dconv = std::max(FunctionDefaults<3>::get_thresh(),
                                  param.dconv());
    const double trantol = vtol / std::min(30.0, double(amo.size()));
    //const double tolloc = 1e-6; // was std::min(1e-6,0.01*dconv) but now trying to avoid unnecessary change // moved to localizer.h
    double update_residual = 0.0, bsh_residual = 0.0;
    subspaceT subspace;
    tensorT Q;
    bool do_this_iter = true;
    bool converged = false;

    // Shrink subspace until stop localizing/canonicalizing--- probably not a good idea
    // int maxsub_save = param.maxsub;
    // param.maxsub = 2;

    for (int iter = 0; iter < param.maxiter(); ++iter) {
        if (world.rank() == 0 and (param.print_level() > 1))
            printf("\nIteration %d at time %.1fs\n\n", iter, wall_time());

        // if (iter > 0 && update_residual < 0.1) {
        //     //do_this_iter = false;
        //     param.maxsub = maxsub_save;
        // }

        if (param.do_localize() && do_this_iter) {
            START_TIMER(world);
            Localizer localizer(world, aobasis, molecule, ao);
            localizer.set_method(param.localize_method());
            MolecularOrbitals<double, 3> mo(amo, aeps, {}, aocc, aset);
            tensorT UT = localizer.compute_localization_matrix(world, mo, iter == 0);
            UT.screen(trantol);
            amo = transform(world, amo, transpose(UT));
            truncate(world, amo);
            normalize(world, amo);

            if (!param.spin_restricted() && param.nbeta() != 0) {

                MolecularOrbitals<double, 3> mo(bmo, beps, {}, bocc, bset);
                tensorT UT = localizer.compute_localization_matrix(world, mo, iter == 0);
                UT.screen(trantol);
                bmo = transform(world, bmo, transpose(UT));
                truncate(world, bmo);
                normalize(world, bmo);
            }
            END_TIMER(world, "localize");
        }

        START_TIMER(world);
        functionT arho = make_density(world, aocc, amo), brho;

        if (param.nbeta()) {
            if (param.spin_restricted()) {
                brho = arho;
            } else {
                brho = make_density(world, bocc, bmo);
            }
        } else {
            brho = functionT(world); // zero
        }
        END_TIMER(world, "Make densities");
        print_meminfo(world.rank(), "Make densities");

        if (iter < 2 || (iter % 10) == 0) {
            START_TIMER(world);
            loadbal(world, arho, brho, arho_old, brho_old, subspace);
            END_TIMER(world, "Load balancing");
            print_meminfo(world.rank(), "Load balancing");
        }
        double da = 0.0, db = 0.0;
        if (iter > 0) {
            da = (arho - arho_old).norm2();
            db = (brho - brho_old).norm2();
            if (world.rank() == 0 and (param.print_level() > 2))
                print("delta rho", da, db, "residuals", bsh_residual,
                      update_residual);

        }
        amo = change_tree_state(amo, reconstructed);

        // screen functions (amo) 
            if (world.rank() == 0) {
                print("\nscreening mos");
                print("vtol " , vtol);
            }
            for (int i = 0; i < amo.size(); ++i) {
                if (world.rank() == 0) print("\nmo " , i);
                amo[i].print_size("\nbefore screening");  
                screen(amo[i], 0.0001, false);
                amo[i].print_size("\nafter screening");  
            }
            world.gop.fence();

        for (int i = 0; i < amo.size(); ++i) {
            save(amo[i], "amo" + stringify(i));
        }

        START_TIMER(world);
        arho_old = arho;
        brho_old = brho;
        functionT rho = arho + brho;
        rho.truncate();

        real_function_3d vnuc;
        if (molecule.parameters.psp_calc()) {
            vnuc = gthpseudopotential->vlocalpot();
        } else if (molecule.parameters.pure_ae()) {
            vnuc = potentialmanager->vnuclear();
        } else {
            vnuc = potentialmanager->vnuclear();
            vnuc = vnuc + gthpseudopotential->vlocalpot();
        }

        double enuclear = inner(rho, vnuc);
        END_TIMER(world, "Nuclear energy");

        START_TIMER(world);
        functionT vcoul = apply(*coulop, rho);
        functionT vlocal;
        END_TIMER(world, "Coulomb");
        print_meminfo(world.rank(), "Coulomb");

        double ecoulomb = 0.5 * inner(rho, vcoul);
        rho.clear(false);
        vlocal = vcoul + vnuc;

        // compute the contribution of the solvent to the local potential
        double epcm = 0.0;
        if (param.pcm_data() != "none") {
            START_TIMER(world);
            functionT vpcm = pcm.compute_pcm_potential(vcoul);
            vlocal += vpcm;
            epcm = pcm.compute_pcm_energy();
            END_TIMER(world, "PCM");
            print_meminfo(world.rank(), "PCM");
        }

        vcoul.clear(false);
        vlocal.truncate();
        double exca = 0.0, excb = 0.0;

        double enla = 0.0, enlb = 0.0;
        vecfuncT Vpsia = apply_potential(world, aocc, amo, vlocal, exca, enla, 0);
        vecfuncT Vpsib;
        if (!param.spin_restricted() && param.nbeta()) {
            Vpsib = apply_potential(world, bocc, bmo, vlocal, excb, enlb, 1);
        } else if (param.nbeta() != 0) {
            enlb = enla;
        }

        double ekina = 0.0, ekinb = 0.0;
        tensorT focka = make_fock_matrix(world, amo, Vpsia, aocc, ekina);
        tensorT fockb = focka;

        if (!param.spin_restricted() && param.nbeta() != 0)
            fockb = make_fock_matrix(world, bmo, Vpsib, bocc, ekinb);
        else if (param.nbeta() != 0) {
            ekinb = ekina;
        }

        if (!param.do_localize() && do_this_iter) {
            tensorT U = diag_fock_matrix(world, focka, amo, Vpsia, aeps, aocc,
                                         FunctionDefaults<3>::get_thresh());
            //rotate_subspace(world, U, subspace, 0, amo.size(), trantol); ??
            if (!param.spin_restricted() && param.nbeta() != 0) {
                U = diag_fock_matrix(world, fockb, bmo, Vpsib, beps, bocc,
                                     FunctionDefaults<3>::get_thresh());
                //rotate_subspace(world, U, subspace, amo.size(), bmo.size(),trantol);
            }
        }

        double enrep = molecule.nuclear_repulsion_energy();
        double ekinetic = ekina + ekinb;
        double enonlocal = enla + enlb;
        double exc = exca + excb;
        double etot = ekinetic + enuclear + ecoulomb + exc + enrep + enonlocal + epcm;
        current_energy = etot;
        //esol = etot;

        if (world.rank() == 0 and (param.print_level() > 1)) {
            //lots of dps for testing Exc stuff
            /*printf("\n              kinetic %32.24f\n", ekinetic);
                printf("         nonlocal psp %32.24f\n", enonlocal);
                printf("   nuclear attraction %32.24f\n", enuclear);
                printf("              coulomb %32.24f\n", ecoulomb);
                printf(" exchange-correlation %32.24f\n", exc);
                printf("    nuclear-repulsion %32.24f\n", enrep);
                printf("                total %32.24f\n\n", etot);*/

            printf("\n              kinetic %16.8f\n", ekinetic);
            printf("         nonlocal psp %16.8f\n", enonlocal);
            printf("   nuclear attraction %16.8f\n", enuclear);
            printf("              coulomb %16.8f\n", ecoulomb);
            printf("                  PCM %16.8f\n", epcm);
            printf(" exchange-correlation %16.8f\n", exc);
            printf("    nuclear-repulsion %16.8f\n", enrep);
            printf("                total %16.8f\n\n", etot);
        }
        e_data.add_data({{"e_kinetic", ekinetic},
                         {"e_local",   enonlocal},
                         {"e_nuclear", enuclear},
                         {"e_coulomb", ecoulomb},
                         {"e_pcm",     epcm},
                         {"e_xc",      exc},
                         {"e_nrep",    enrep},
                         {"e_tot",     etot}});

        if (iter > 0) {
            //print("##convergence criteria: density delta=", da < dconv * molecule.natom() && db < dconv * molecule.natom(), ", bsh_residual=", (param.conv_only_dens || bsh_residual < 5.0*dconv));
            if (da < dconv * std::max(size_t(5), molecule.natom()) && db < dconv * std::max(size_t(5), molecule.natom())
                && (param.get<bool>("conv_only_dens") || bsh_residual < 5.0 * dconv))
                converged = true;
            // previous conv was too tight for small systems
            // if (da < dconv * molecule.natom() && db < dconv * molecule.natom()
            //     && (param.conv_only_dens || bsh_residual < 5.0 * dconv)) converged=true;

            // do diagonalization etc if this is the last iteration, even if the calculation didn't converge
            if (converged || iter == param.maxiter() - 1) {
                if (world.rank() == 0 && converged and (param.print_level() > 1)) {
                    print("\nConverged!\n");
                    converged_for_thresh=param.econv();
                }

                // Diagonalize to get the eigenvalues and if desired the final eigenvectors
                tensorT U;
                START_TIMER(world);
                tensorT overlap = matrix_inner(world, amo, amo, true);
                END_TIMER(world, "Overlap");

                START_TIMER(world);
                sygvp(world, focka, overlap, 1, U, aeps);
                END_TIMER(world, "focka eigen sol");

                if (!param.do_localize()) {
                    START_TIMER(world);
                    amo = transform(world, amo, U, trantol, true);
                    truncate(world, amo);
                    normalize(world, amo);
                    END_TIMER(world, "Transform MOs");
                }
                if (param.nbeta() != 0 && !param.spin_restricted()) {

                    START_TIMER(world);
                    overlap = matrix_inner(world, bmo, bmo, true);
                    END_TIMER(world, "Overlap");

                    START_TIMER(world);
                    sygvp(world, fockb, overlap, 1, U, beps);
                    END_TIMER(world, "fockb eigen sol");

                    if (!param.do_localize()) {
                        START_TIMER(world);
                        bmo = transform(world, bmo, U, trantol, true);
                        truncate(world, bmo);
                        normalize(world, bmo);
                        END_TIMER(world, "Transform MOs");
                    }
                }

                if (world.rank() == 0 and (param.print_level() > 1)) {
                    print(" ");
                    print("alpha eigenvalues");
                    print(aeps);
                    if (param.nbeta() != 0 && !param.spin_restricted()) {
                        print("beta eigenvalues");
                        print(beps);
                    }


                    // write eigenvalues etc to a file at the same time for plotting DOS etc.
                    FILE *f = 0;
                    if (param.nbeta() != 0 && !param.spin_restricted()) {
                        std::string name=std::string(param.prefix()+".energies_alpha.dat");
                        f = fopen(name.c_str(), "w");
                    } else {
                        std::string name=param.prefix()+".energies.dat";
                        f = fopen(name.c_str(), "w");
                    }

                    long nmo = amo.size();
                    fprintf(f, "# %8li\n", nmo);
                    for (long i = 0; i < nmo; ++i) {
                        fprintf(f, "%13.8f\n", aeps(i));
                    }
                    fclose(f);

                    if (param.nbeta() != 0 && !param.spin_restricted()) {
                        long nmo = bmo.size();
                        FILE *f = 0;
                        std::string name=param.prefix()+".energies_beta.dat";
                        f = fopen(name.c_str(), "w");

                        fprintf(f, "# %8li\n", nmo);
                        for (long i = 0; i < nmo; ++i) {
                            fprintf(f, "%13.8f\t", beps(i));
                        }
                        fclose(f);
                    }

                }

                if (param.do_localize()) {
                    // Restore the diagonal elements for the analysis
                    for (unsigned int i = 0; i < amo.size(); ++i)
                        aeps[i] = focka(i, i);
                    if (param.nbeta() != 0 && !param.spin_restricted())
                        for (unsigned int i = 0; i < bmo.size(); ++i)
                            beps[i] = fockb(i, i);
                }

                break;
            }

        }

        update_subspace(world, Vpsia, Vpsib, focka, fockb, subspace, Q,
                        bsh_residual, update_residual);

    }

    // compute the dipole moment
    functionT rho = make_density(world, aocc, amo);
    if (!param.spin_restricted()) {
        if (param.nbeta())
            rho += make_density(world, bocc, bmo);
    } else {
        rho.scale(2.0);
    }
    dipole(world, rho);

    if (world.rank() == 0 and (param.print_level() > 1)) {
        if (param.do_localize())
            print(
                    "Orbitals are localized - energies are diagonal Fock matrix elements\n");
        else
            print("Orbitals are eigenvectors - energies are eigenvalues\n");
        if (param.nwfile() == "none") print("Analysis of alpha MO vectors");
    }

    if (param.nwfile() == "none") {
        analyze_vectors(world, amo, aocc, aeps);
        if (param.nbeta() != 0 && !param.spin_restricted()) {
            if (world.rank() == 0 and (param.print_level() > 1))
                print("Analysis of beta MO vectors");

            analyze_vectors(world, bmo, bocc, beps);
        }
    }

}        // end solve function

} // namespace madness
