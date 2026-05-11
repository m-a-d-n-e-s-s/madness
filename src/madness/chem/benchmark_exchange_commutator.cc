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
#include <madness/chem/CCPotentials.h>
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
    const bool do_ue = parser.key_exists("do_ue");
    const bool do_kffk = parser.key_exists("do_kffk");

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

        Info info = make_info(world, amo, nemo->R, nemo->R_square, nemo->molecule());
        // Info info = make_info(world, amo, one, one, nemo->molecule());
        info.parameters.print("cc2");
        print("numerical parameters: k, eps(3D), eps(6D) =",
              FunctionDefaults<3>::get_k(),
              FunctionDefaults<3>::get_thresh(),
              FunctionDefaults<6>::get_thresh());

        int i=0;
        int j=0;
        if (amo.size()>1) i=1;
        if (amo.size()>2) j=2;
        print("working on orbitals: ",i,j);

        CCFunction<double,3> phi_i(amo[i], i, HOLE);
        CCFunction<double,3> phi_j(amo[j], j, HOLE);

        auto ao=orthonormalize_canonical(nemo->get_calc()->ao);
        print("error computed by projecting on the ao basis",nemo->get_calc()->aobasis.get_name(), "size",ao.size());
        // ExchangeCommutator instance carrying the orthonormal AO basis used
        // as the observer set inside diagnose().
        ExchangeCommutator ec(ao);

        // Score lambda: project Kf, fK, and KffK onto the AO basis and print
        // errors for whichever pieces are present in the result.  Empty
        // entries in the KffKResult are skipped by diagnose().
        auto score_full = [&](const ExchangeCommutator::KffKResult& r,
                              bool include_K2) {
            auto d = ec.diagnose(
                    world, amo, R2amo, phi_i.function, phi_j.function,
                    r.Kf, r.fK, r.KffK, lrfparam,
                    /*verbose=*/true,
                    include_K2);
            ExchangeCommutator::print_report(r, &d);
        };

        // -------------------- 6D reference: full commutator via 6D path ---
        // Run first: independent algorithm (apply_Kfxy + make_f_xy_macrotask)
        // gives a baseline error level against which the split-α numbers
        // below are interpreted.  Pure-6D pair functions, projected onto the
        // same harmonic basis as the split-α path.
//        print("\n========== 6D reference: Kf, fK, KffK via apply_Kfxy ==========");
//        {
//            score_full(ExchangeCommutator::apply_KffK_6d(
//                    world, phi_i, phi_j, info),
//                       /*include_K2=*/true);
//        }

//        // -------------------- regression: K̂₁/Kf only --------------------
//        // Reproduce the working K̂₁-only configuration before adding fK.
//        // If err_Kf grows away from ~1e-5, the K̂₁ path has regressed.
//        print("\n========== regression: K̂₁ piece only (Kf alone) ==========");
//        {
//            ExchangeCommutator::SplitAlphaOptions opt;
//            opt.alpha_star               = 1.e4;
//            opt.assemble_fK              = false;
//            opt.include_symmetry_mirror  = false;
//            score_full(ExchangeCommutator::apply_KffK_lowrank_split_alpha(
//                    world, phi_i, phi_j, info, lrfparam, opt),
//                       /*include_K2=*/false);
//        }

        // -------------------- three-range C: medium=6D at alpha_hi=1e4 ----
        // 2026-05-05 follow-up #2 (task #13): hybrid path.  Diffuse slab
        // built as a coarse-grid LRF (rank ~160).  Medium slab handled by
        // a *pure 6D* apply with a partial-Coulomb operator (only the
        // medium Gaussians).  Tight slab still discarded.
        //
        // Decision criterion (from the prior plan): pure-6D medium piece
        // ≤ 10 GB at err_KffK ≤ 5e-5 → pursue; ≥ 50 GB → drop.
        if (0) {
            LowRankFunctionParameters lrfparam_diffuse;
            lrfparam_diffuse.set_derived_value("f12type", std::string("slater"));
            lrfparam_diffuse.set_derived_value("radius",         5.0);
            lrfparam_diffuse.set_derived_value("volume_element", 2.0);
            lrfparam_diffuse.set_derived_value("tol",            1.e-5);
            lrfparam_diffuse.set_derived_value("tempered",
                    std::vector<double>({1.e-2, 1.e0, 9.0}));
            lrfparam_diffuse.set_derived_value("lmax",           2);

            ExchangeCommutator::ThreeRangeOptions opt3;
            opt3.alpha_lo          = 1.0;
            opt3.alpha_hi          = 1.e4;
            opt3.lrfparam_diffuse  = lrfparam_diffuse;
            opt3.lrfparam_medium   = lrfparam;   // unused when medium_use_6d=true
            opt3.medium_use_taylor = false;
            opt3.medium_use_6d     = true;

            print("\n========== three-range C: medium=6D, alpha_lo=1, alpha_hi=1e4 ==========");
            score_full(ExchangeCommutator::apply_KffK_lowrank_three_range(
                           world, phi_i, phi_j, info, opt3),
                       /*include_K2=*/true);
        }

        // compute Ue term and its diagnosis
        if (1) {
            CorrelationFactor cf(world, 1.0, 1.e-10, nemo->molecule());
            auto ao=orthonormalize_canonical(nemo->get_calc()->ao);
            auto U1nuc = nemo->ncf->U1vec();
            if (do_ue) {
                for (int tmode : {-2,1,3}) {
                    print_header2("setting truncate mode to"+std::to_string(tmode));
                    cf.set_truncate_mode(tmode);
                    real_convolution_6d op_mod = BSHOperator<6>(world, sqrt(2.0), info.parameters.lo(), info.parameters.thresh_bsh_6D());
                    op_mod.modified() = true;
                    const double thresh6 = FunctionDefaults<6>::get_thresh();
                    auto Uphi_local    = cf.apply_U_local(phi_i.function, phi_j.function, op_mod, thresh6);
                    auto Uphi_semilocal = cf.apply_U_semilocal(phi_i.function, phi_j.function, op_mod, thresh6);
                    auto Uphi_mixed    = cf.apply_U_mixed_commutator(phi_i.function, phi_j.function, op_mod, thresh6, U1nuc);
                    auto diag = cf.diagnose_Ue(Uphi_local, Uphi_semilocal,
                                               phi_i.function, phi_j.function, ao,
                                               Uphi_mixed, U1nuc);
                    print("diagnosis for Ue term:");
                    print("local part error:      ", diag.error_local);
                    print("semi-local part error: ", diag.error_semilocal);
                    print("mixed-comm part error: ", diag.error_mixed);
                    print("time in diagnostics:   ", diag.time_ref);

                    // G Ue diagnostics: apply G = (T - E)^{-1} to Ue|ij>, then compare
                    // the 6D projection <ab|G Ue|ij> against the 3D Schwinger quadrature.
                    const auto& eps = nemo->get_calc()->aeps;
                    const double energy_ij = eps(i) + eps(j);
                    const double mu_ij = sqrt(-2.0 * energy_ij);
                    real_convolution_6d G6d = BSHOperator<6>(world, mu_ij,
                                                             info.parameters.lo(),
                                                             info.parameters.thresh_bsh_6D());
                    auto GUphi_local     = apply(G6d, Uphi_local);
                    auto GUphi_semilocal = apply(G6d, Uphi_semilocal);
                    auto GUphi_mixed     = apply(G6d, Uphi_mixed);
                    GUphi_local.print_size("G U_local|ij>");
                    GUphi_semilocal.print_size("G U_semilocal|ij>");
                    GUphi_mixed.print_size("G U_mixed|ij>");

                    print("\ndiagnosis for G Ue term (energy =", energy_ij, "):");
                    auto gue = cf.diagnose_GUe(GUphi_local, GUphi_semilocal,
                                               phi_i.function, phi_j.function, ao, energy_ij,
                                               GUphi_mixed, U1nuc);
                    print("G Ue local part error:      ", gue.error_local);
                    print("G Ue semilocal part error:  ", gue.error_semilocal);
                    print("G Ue mixed-comm part error: ", gue.error_mixed);
                    print("time in G Ue diagnostics:   ", gue.time);
                }
            }

            // --- G [K̂, f₁₂] |ij⟩ diagnostic ---
            // Test two algorithm variants against the 3D Schwinger quadrature.
            // The caller applies G externally (converting decomposed CCPairFunctions
            // to pure 6D first where needed) and passes the G-applied pieces as
            // CCPairFunction vectors to diagnose_GKffK.
            //
            // Helper: apply G6d to a CCPairFunction vector, returning CCPairFunction
            // wrappers around the resulting pure-6D functions.  Decomposed entries
            // (LRF sum Σ_ρ g_ρ⊗h_ρ) are first collapsed to a single pure-6D
            // function via CompositeFactory, then G is applied once.
            auto apply_G_to_pairs = [&](const std::vector<CCPairFunction<double,6>>& v,
                                        const real_convolution_6d& G6d)
                    -> std::vector<CCPairFunction<double,6>>
            {
                std::vector<CCPairFunction<double,6>> result;
                for (const auto& f : v) {
                    if (!f.is_assigned()) continue;
                    result.emplace_back(apply(G6d, f));
                }
                return result;
            };


            if (do_kffk) {
                // K̂φᵢ and K̂φⱼ — needed by the Schwinger quadrature, computed once.
                Exchange<double,3> K(world,info.parameters.lo());
                K.set_bra_and_ket(info.mo_bra,info.mo_ket);
                const auto Kphi_i = K(phi_i.function);
                const auto Kphi_j = K(phi_j.function);

                const auto& eps    = nemo->get_calc()->aeps;
                const double energy_ij = eps(i) + eps(j);
                const double mu_ij     = std::sqrt(-2.0 * energy_ij);
                real_convolution_6d G6d = BSHOperator<6>(world, mu_ij,
                                                          info.parameters.lo(),
                                                          info.parameters.thresh_bsh_6D());

                // --- variant 1: 6D reference algorithm (pure-6D Kf and fK) ---
                if (0) {
                    print("\n========== G [K, f] |ij> diagnostic — 6D reference ==========");
                    auto kffk_6d = ExchangeCommutator::apply_KffK_6d(world, phi_i, phi_j, info);
                    score_full(kffk_6d, /*include_K2=*/true);  // print Kf/fK errors for reference

                    auto GKf_6d_cc = apply_G_to_pairs(kffk_6d.Kf, G6d);
                    auto GfK_6d_cc = apply_G_to_pairs(kffk_6d.fK, G6d);
                    GKf_6d_cc[0].get_function().print_size("G Kf|ij> (6D)");
                    GfK_6d_cc[0].get_function().print_size("G fK|ij> (6D)");

                    auto gkffk_6d = ec.diagnose_GKffK(world, GKf_6d_cc, GfK_6d_cc,
                                                        phi_i.function, phi_j.function,
                                                        Kphi_i, Kphi_j, info, energy_ij);
                    print("6D ref  — G Kf  error:", gkffk_6d.error_GKf);
                    print("6D ref  — G fK  error:", gkffk_6d.error_GfK);
                    print("6D ref  — G[K,f] error:", gkffk_6d.error_GKffK);
                    print("6D ref  — time:", gkffk_6d.time);
                }

                // --- variant 2: split-α LRF algorithm ---
                print("\n========== G [K, f] |ij> diagnostic — split-alpha LRF ==========");
                ExchangeCommutator::SplitAlphaOptions sa_opt;
                sa_opt.alpha_star=100.0;
                auto kffk_lrf = ExchangeCommutator::apply_KffK_lowrank_split_alpha(
                        world, phi_i, phi_j, info, lrfparam, sa_opt);
                score_full(kffk_lrf, /*include_K2=*/true);  // print Kf/fK errors for reference

                auto GKf_lrf_cc = apply_G_to_pairs(kffk_lrf.Kf, G6d);
                auto GfK_lrf_cc = apply_G_to_pairs(kffk_lrf.fK, G6d);

                auto gkffk_lrf = ec.diagnose_GKffK(world, GKf_lrf_cc, GfK_lrf_cc,
                                                     phi_i.function, phi_j.function,
                                                     Kphi_i, Kphi_j, info, energy_ij);
                print("LRF     — G Kf  error:", gkffk_lrf.error_GKf);
                print("LRF     — G fK  error:", gkffk_lrf.error_GfK);
                print("LRF     — G[K,f] error:", gkffk_lrf.error_GKffK);
                print("LRF     — time:", gkffk_lrf.time);
            }

        }


    } catch (std::exception& e) {
        print("an error occurred");
        print(e.what());
        status = 1;
    }
    finalize();
    return status;
}
