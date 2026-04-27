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
        FunctionDefaults<3>::set_thresh(thresh);

        auto amo = nemo->get_calc()->get_amo();
        auto R2amo = nemo->R_square*(nemo->get_calc()->get_amo());

//        real_function_3d gauss=real_factory_3d(world).functor(
//            [](const coord_3d& r){ return std::exp(-r.normf()); });
//        real_function_3d one=real_factory_3d(world).functor([](const coord_3d&){ return 1.0; });
//        const double gauss_norm = gauss.norm2();
//        gauss.scale(1.0/gauss_norm);
//        print("orbital normalized: ||phi||_before =", gauss_norm,
//              "  ||phi||_after =", gauss.norm2());
        // amo  =std::vector<real_function_3d>({gauss});
        // R2amo=std::vector<real_function_3d>({gauss});   // info.R²=one, so R²·k == k

        Info info = make_info(world, amo, nemo->R, nemo->R_square, nemo->molecule());
        // Info info = make_info(world, amo, one, one, nemo->molecule());
        info.parameters.print("cc2");
        print("numerical parameters: k, eps(3D), eps(6D) =",
              FunctionDefaults<3>::get_k(),
              FunctionDefaults<3>::get_thresh(),
              FunctionDefaults<6>::get_thresh());

        const real_function_3d phi = amo[0];
        CCFunction<double,3> phi_i(phi, 0, HOLE);
        CCFunction<double,3> phi_j(phi, 0, HOLE);

        // Score lambda: take the algorithm result, project onto the harmonic
        // basis via diagnose(), pass empty fK / KffK because we are testing the
        // Kf piece only.  include_K2 is steered from the caller so the same
        // lambda works for both the K̂₁-only regression and the K̂₁+K̂₂ run.
        auto score_Kf = [&](const ExchangeCommutator::KffKResult& r,
                             bool include_K2) {
            std::vector<CCPairFunction<double,6>> empty_fK;
            std::vector<CCPairFunction<double,6>> empty_KffK;
            auto d = ExchangeCommutator::diagnose(
                    world, amo, R2amo, phi, phi,
                    r.Kf, empty_fK, empty_KffK, lrfparam,
                    /*verbose=*/true,
                    include_K2);
            ExchangeCommutator::print_report(r, &d);
        };

        // -------------------- regression: K̂₁/Kf only --------------------
        // Reproduce the working K̂₁-only configuration before re-enabling K̂₂.
        // If err_Kf grows away from ~3e-3, the K̂₁ path itself has regressed
        // and the K̂₁+K̂₂ result is not meaningful.
        print("\n========== regression: K̂₁ piece only ==========");
        {
            ExchangeCommutator::SplitAlphaOptions opt;
            opt.alpha_star               = 1.e4;
            opt.assemble_fK              = false;
            opt.include_symmetry_mirror  = false;
            score_Kf(ExchangeCommutator::apply_KffK_lowrank_split_alpha(
                    world, phi_i, phi_j, info, lrfparam, opt),
                     /*include_K2=*/false);
        }

        // -------------------- K̂₁ + K̂₂ via swap_particles ----------------
        print("\n========== K̂₁ + K̂₂ piece (swap_particles enabled) ==========");
        {
            ExchangeCommutator::SplitAlphaOptions opt;
            opt.alpha_star               = 1.e4;
            opt.assemble_fK              = false;   // still Kf-only
            opt.include_symmetry_mirror  = true;    // re-enable K̂₂ via swap
            score_Kf(ExchangeCommutator::apply_KffK_lowrank_split_alpha(
                    world, phi_i, phi_j, info, lrfparam, opt),
                     /*include_K2=*/true);
        }
        // fK / lrf / lrf-direct paths intentionally disabled for this run.

    } catch (std::exception& e) {
        print("an error occurred");
        print(e.what());
        status = 1;
    }
    finalize();
    return status;
}
