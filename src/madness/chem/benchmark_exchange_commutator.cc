//
// benchmark_exchange_commutator.cc
//
// Driver for the four [K̂, f₁₂]|ij⟩ algorithms exposed by
// madness::ExchangeCommutator:
//   * 6D           — apply_Kfxy + make_f_xy_macrotask(Kx,…)
//   * lrf          — CCPotentials::apply_KffK_low_rank  (canonical)
//   * lrf-direct   — CCPotentials::apply_KffK_low_rank_direct  (per-k)
//   * split-alpha  — scalar-only diagnostic (six Hermiticity-canceled terms)
//
// Each path returns a uniform KffKResult (wall / memory / rank) and is
// optionally scored via the harmonic-basis ExchangeCommutator::diagnose.
//

#include <madness.h>
#include <madness/chem/exchange_commutator.h>
#include <madness/chem/nemo.h>
#include <madness/chem/CCStructures.h>
#include <madness/chem/CCParameters.h>
#include <madness/chem/lowrankfunction.h>
#include <madness/chem/ccpairfunction.h>
#include <madness/world/test_utilities.h>

using namespace madness;

namespace {

Info make_info(World& world, const std::vector<real_function_3d>& phivec,
    const real_function_3d R, const real_function_3d R_square, const Molecule& molecule) {
    real_function_3d one = real_factory_3d(world).f(
            [](const coord_3d&){ return 1.0; });
    Info info;
    info.mo_bra = phivec*R_square;
    info.mo_ket = phivec;
    info.R = R;
    info.R_square = R_square;
    info.molecular_coordinates = molecule.get_all_coords_vec();
    info.parameters = CCParameters();
    return info;
}

} // namespace

int main(int argc, char **argv) {
    World& world = madness::initialize(argc, argv);
    startup(world, argc, argv);
    commandlineparser parser(argc, argv);

    const int    k      = parser.key_exists("k")      ? std::atoi(parser.value("k").c_str())   : 6;
    const double thresh = parser.key_exists("thresh") ? std::stod(parser.value("thresh"))       : 3.e-5;

    FunctionDefaults<6>::set_tensor_type(TT_2D);
    FunctionDefaults<3>::set_truncate_mode(3);
    FunctionDefaults<6>::set_truncate_mode(3);
    for (int d : {1,2,3,4,5,6}) {
        // set thresh / k / cell for each dimensionality used downstream.
    }
    FunctionDefaults<3>::set_thresh(thresh);
    FunctionDefaults<6>::set_thresh(1.e-3);
    FunctionDefaults<3>::set_k(k);
    FunctionDefaults<6>::set_k(k);
    FunctionDefaults<3>::set_cubic_cell(-10., 10.);
    FunctionDefaults<6>::set_cubic_cell(-10., 10.);

    print("numerical parameters: k, eps(3D), eps(6D) =",
          FunctionDefaults<3>::get_k(),
          FunctionDefaults<3>::get_thresh(),
          FunctionDefaults<6>::get_thresh());

    LowRankFunctionParameters lrfparam;
    lrfparam.set_derived_value("f12type", std::string("slater"));
    lrfparam.read_and_set_derived_values(world, parser, "grid");
    lrfparam.set_derived_value("radius",         2.5);
    lrfparam.set_derived_value("volume_element", 0.2);
    lrfparam.set_derived_value("tol",            1.e-5);
    lrfparam.set_derived_value("tempered",       std::vector<double>({1.e-1, 1.e1, 9.0}));
    lrfparam.set_derived_value("lmax",           2);
    lrfparam.print("grid");

    int status = 0;
    try {
        // Build SCF orbitals (sample molecule) — reuse Nemo for a realistic test.
        std::shared_ptr<Nemo> nemo(new Nemo(world, parser));
        nemo->value();
        auto cell=FunctionDefaults<3>::get_cell();
        FunctionDefaults<6>::set_cubic_cell(cell(0,0),cell(0,1));

        auto amo = nemo->get_calc()->get_amo();
        Info info = make_info(world, amo, nemo->R, nemo->R_square, nemo->molecule());
        info.parameters.print("cc2");
        print("numerical parameters: k, eps(3D), eps(6D) =",
              FunctionDefaults<3>::get_k(),
              FunctionDefaults<3>::get_thresh(),
              FunctionDefaults<6>::get_thresh());

        const real_function_3d phi = amo[0];
        CCFunction<double,3> phi_i(phi, 0, HOLE);
        CCFunction<double,3> phi_j(phi, 0, HOLE);

        auto score = [&](const ExchangeCommutator::KffKResult& r) {
            auto d = ExchangeCommutator::diagnose(
                    world, amo, phi, phi, r.Kf, r.fK, r.KffK, lrfparam);
            ExchangeCommutator::print_report(r, &d);
        };

        // -------------------- 6D algorithm --------------------------------
        score(ExchangeCommutator::apply_KffK_6d(world, phi_i, phi_j, info));

        // -------------------- LRF canonical -------------------------------
        score(ExchangeCommutator::apply_KffK_lowrank(
                world, phi_i, phi_j, info, lrfparam));

        // -------------------- LRF per-k direct ----------------------------
        score(ExchangeCommutator::apply_KffK_lowrank_direct(
                world, phi_i, phi_j, info, lrfparam));

        // -------------------- split-α 6D assembly -------------------------
        {
            ExchangeCommutator::SplitAlphaOptions opt;
            opt.alpha_star = 1.e4;
            score(ExchangeCommutator::apply_KffK_lowrank_split_alpha(
                    world, phi_i, phi_j, info, lrfparam, opt));
        }

    } catch (std::exception& e) {
        print("an error occurred");
        print(e.what());
        status = 1;
    }
    finalize();
    return status;
}
