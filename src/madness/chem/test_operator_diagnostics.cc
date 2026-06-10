// test_operator_diagnostics.cc
// Unit tests for DiagnosticMatrix (operator_diagnostics.h).
//
// NDIM = pair-space dimension (template parameter); LDIM = NDIM/2 = one-particle space.
//
//   test_ctor_orthonormalization<NDIM> — verify the ctor detects and fixes a
//       non-orthonormal basis; verify an already-orthonormal one passes silently.
//
//   test_project_ab<NDIM> — project a decomposed CCPairFunction ket and an NDIM-D
//       function ket onto the AO basis; compare both against the analytic reference
//       ref(a,b) = <ao[a]|phi_i> * <ao[b]|phi_j>.  Works for any NDIM.
//
//   test_diagnose_Ue<NDIM> — end-to-end test through CorrelationFactor::diagnose_Ue.
//       Only active for NDIM=6 (CorrelationFactor is hardcoded for the 3D one-particle space).

#include <madness/mra/mra.h>
#include <madness/mra/commandlineparser.h>
#include <madness/chem/CCStructures.h>                  // CCConvolutionOperator before CCPairFunction instantiation
#include <madness/chem/electronic_correlation_factor.h> // CorrelationFactor; pulls in operator_diagnostics.h
#include <madness/world/test_utilities.h>

using namespace madness;

// Non-capturing Gaussian functors: FunctionFactory::f() requires a raw
// function pointer, so exponents must be compile-time constants.
namespace gaussians {
    template<std::size_t LDIM>
    double g05(const Vector<double,LDIM>& r) { return exp(-0.5 * inner(r,r)); }
    template<std::size_t LDIM>
    double g10(const Vector<double,LDIM>& r) { return exp(-1.0 * inner(r,r)); }
    template<std::size_t LDIM>
    double g20(const Vector<double,LDIM>& r) { return exp(-2.0 * inner(r,r)); }
    template<std::size_t LDIM>
    double g30(const Vector<double,LDIM>& r) { return exp(-3.0 * inner(r,r)); }
}

// ---------------------------------------------------------------------------
// test 1: constructor orthonormalization
// ---------------------------------------------------------------------------

template<std::size_t NDIM>
int test_ctor_orthonormalization(World& world) {
    static_assert(NDIM % 2 == 0, "NDIM must be even");
    constexpr std::size_t LDIM = NDIM / 2;
    test_output t("DiagnosticMatrix ctor orthonormalization <NDIM=" + std::to_string(NDIM) + ">");

    const double thresh = FunctionDefaults<LDIM>::get_thresh();

    // Three overlapping Gaussians — intentionally non-orthonormal
    std::vector<Function<double,LDIM>> raw = {
        FunctionFactory<double,LDIM>(world).f(gaussians::g10<LDIM>),
        FunctionFactory<double,LDIM>(world).f(gaussians::g20<LDIM>),
        FunctionFactory<double,LDIM>(world).f(gaussians::g30<LDIM>),
    };

    Tensor<double> S_raw = matrix_inner(world, raw, raw, /*sym=*/true);
    Tensor<double> d = copy(S_raw);
    for (int i = 0; i < 3; ++i) d(i,i) -= 1.0;
    t.checkpoint(d.normf() > 1e-3, "raw basis is non-orthonormal");

    // Constructor must orthonormalize automatically; CTAD deduces DiagnosticMatrix<double,NDIM>
    DiagnosticMatrix dm(world, raw);
    Tensor<double> S = matrix_inner(world, dm.ao_basis, dm.ao_basis, /*sym=*/true);
    d = copy(S);
    for (int i = 0; i < dm.nbasis(); ++i) d(i,i) -= 1.0;
    t.checkpoint(d.normf(), 10.0 * thresh, "stored basis is orthonormal after ctor");

    // Pre-orthonormal basis: ctor must accept it without re-processing
    auto ortho = orthonormalize_symmetric(raw);
    DiagnosticMatrix dm2(world, ortho);
    Tensor<double> S2 = matrix_inner(world, dm2.ao_basis, dm2.ao_basis, /*sym=*/true);
    d = copy(S2);
    for (int i = 0; i < dm2.nbasis(); ++i) d(i,i) -= 1.0;
    t.checkpoint(d.normf(), 10.0 * thresh, "pre-orthonormal basis accepted unchanged");

    return t.end();
}

// ---------------------------------------------------------------------------
// test 2: project_ab with decomposed and NDIM-D kets
// ---------------------------------------------------------------------------

// Non-capturing NDIM-D Gaussian: phi_i(r1)*phi_j(r2) = exp(-|r1|^2 - 2|r2|^2)
template<std::size_t NDIM>
double g_phi_i_phi_j(const Vector<double,NDIM>& r) {
    constexpr std::size_t LDIM = NDIM / 2;
    double r1sq = 0.0, r2sq = 0.0;
    for (std::size_t i = 0; i < LDIM; ++i) {
        r1sq += r[i]        * r[i];
        r2sq += r[i + LDIM] * r[i + LDIM];
    }
    return exp(-1.0*r1sq - 2.0*r2sq);
}

template<std::size_t NDIM>
int test_project_ab(World& world) {
    static_assert(NDIM % 2 == 0, "NDIM must be even");
    constexpr std::size_t LDIM = NDIM / 2;
    test_output t("DiagnosticMatrix project_ab <NDIM=" + std::to_string(NDIM) + ">");

    const double thresh = FunctionDefaults<LDIM>::get_thresh();

    auto phi_i = FunctionFactory<double,LDIM>(world).f(gaussians::g10<LDIM>);
    auto phi_j = FunctionFactory<double,LDIM>(world).f(gaussians::g20<LDIM>);

    // AO observer basis (non-orthonormal; ctor will fix it).
    // CTAD deduces DiagnosticMatrix<double,NDIM> from vector<Function<double,LDIM>>.
    DiagnosticMatrix dm(world, std::vector<Function<double,LDIM>>{
        FunctionFactory<double,LDIM>(world).f(gaussians::g10<LDIM>),
        FunctionFactory<double,LDIM>(world).f(gaussians::g20<LDIM>),
        FunctionFactory<double,LDIM>(world).f(gaussians::g05<LDIM>),
    });
    const int nb = dm.nbasis();

    // Analytic reference: ref(a,b) = <ao[a]|phi_i> * <ao[b]|phi_j>.
    // Use matrix_inner to avoid ambiguity with CCPairFunction inner() overloads.
    std::vector<Function<double,LDIM>> vi = {phi_i}, vj = {phi_j};
    Tensor<double> ov_i = matrix_inner(world, dm.ao_basis, vi);
    Tensor<double> ov_j = matrix_inner(world, dm.ao_basis, vj);
    Tensor<double> ref(nb, nb);
    for (int a = 0; a < nb; ++a)
        for (int b = 0; b < nb; ++b)
            ref(a, b) = ov_i(a, 0) * ov_j(b, 0);

    // --- decomposed CCPairFunction ket ---
    auto result_dec = dm.project_ab({ CCPairFunction<double,NDIM>(vi, vj) });
    t.checkpoint((result_dec - ref).normf(), 10.0 * thresh, "project_ab: decomposed ket");

    // --- full NDIM-D function ket on the grid ---
    auto ket_nd = FunctionFactory<double,NDIM>(world).f(g_phi_i_phi_j<NDIM>);
    auto result_nd = dm.project_ab(ket_nd);
    t.checkpoint((result_nd - ref).normf(), 10.0 * thresh,
                 "project_ab: " + std::to_string(NDIM) + "D ket");

    return t.end();
}

// ---------------------------------------------------------------------------
// test 3: diagnose_Ue via CorrelationFactor (NDIM=6 only)
// ---------------------------------------------------------------------------

template<std::size_t NDIM>
int test_diagnose_Ue(World& world) {
    static_assert(NDIM % 2 == 0, "NDIM must be even");
    constexpr std::size_t LDIM = NDIM / 2;

    // CorrelationFactor is hardcoded for the 3D one-particle space.
    if constexpr (LDIM == 3) {
        test_output t("DiagnosticMatrix diagnose_Ue <NDIM=" + std::to_string(NDIM) + ">");

        const double thresh3 = FunctionDefaults<3>::get_thresh();
        const double thresh6 = FunctionDefaults<6>::get_thresh();
        const double lo = 1e-6;

        auto phi_i = FunctionFactory<double,3>(world).f(gaussians::g10<3>);
        auto phi_j = FunctionFactory<double,3>(world).f(gaussians::g20<3>);

        std::vector<Function<double,3>> ao_raw = {
            FunctionFactory<double,3>(world).f(gaussians::g10<3>),
            FunctionFactory<double,3>(world).f(gaussians::g20<3>),
            FunctionFactory<double,3>(world).f(gaussians::g05<3>),
        };

        const double gamma   = 1.0;
        const double bsh_eps = -2.0;
        CorrelationFactor cf(world, gamma, 1e-10, lo);

        real_convolution_6d op_mod = BSHOperator<6>(world, sqrt(-2.0 * bsh_eps), lo, 1.e-6);
        op_mod.modified() = true;

        auto Uphi_local     = cf.apply_U_local    (phi_i, phi_j, op_mod, thresh6);
        auto Uphi_semilocal = cf.apply_U_semilocal(phi_i, phi_j, op_mod, thresh6);

        auto dm = cf.diagnose_Ue(Uphi_local, Uphi_semilocal, phi_i, phi_j, ao_raw);
        dm.print_report("test_diagnose_Ue");

        t.checkpoint(dm.entries["local"].error, thresh6, "diagnose_Ue: local error");
        t.checkpoint(dm.entries["semilocal"].error,thresh6, "diagnose_Ue: semilocal error");

        return t.end();
    }
    return 0;
}

// ---------------------------------------------------------------------------
// test 4: project_Gab — Schwinger bra with CCPairFunction ket vs explicit ref
// ---------------------------------------------------------------------------

template<std::size_t NDIM>
int test_GUe(World& world) {
    static_assert(NDIM % 2 == 0, "NDIM must be even");
    constexpr std::size_t LDIM = NDIM / 2;
    test_output t("DiagnosticMatrix project_Gab <NDIM=" + std::to_string(NDIM) + ">");

    const double thresh = FunctionDefaults<LDIM>::get_thresh();
    const double lo     = 1e-6;
    const double energy = -2.0;  // must be < 0

    auto phi_i = FunctionFactory<double,LDIM>(world).f(gaussians::g10<LDIM>);
    auto phi_j = FunctionFactory<double,LDIM>(world).f(gaussians::g20<LDIM>);

    // AO observer basis (non-orthonormal; ctor will fix it)
    DiagnosticMatrix<double,NDIM> dm(world, std::vector<Function<double,LDIM>>{
        FunctionFactory<double,LDIM>(world).f(gaussians::g10<LDIM>),
        FunctionFactory<double,LDIM>(world).f(gaussians::g20<LDIM>),
        FunctionFactory<double,LDIM>(world).f(gaussians::g05<LDIM>),
    });
    const int nb = dm.nbasis();

    std::vector<Function<double,LDIM>> vi = {phi_i}, vj = {phi_j};
    auto ket = std::vector<CCPairFunction<double,NDIM>>{ CCPairFunction<double,NDIM>(vi, vj) };

    // Result: project_Gab via project_xy (CCPairFunction-ket code path + build_Gab_bra)
    auto result = dm.project_Gab(ket, energy, lo);

    // Analytic ref: Σ_k inner(bra[k].get_a()[a], φ_i) · inner(bra[k].get_b()[b], φ_j)
    // Uses raw 3D matrix_inner — different code path from project_xy.
    auto bra = dm.build_Gab_bra(energy, lo);
    Tensor<double> ref(nb, nb);
    for (const auto& bk : bra) {
        Tensor<double> ov_i = matrix_inner(world, bk.get_a(), vi);  // [nb,1]
        Tensor<double> ov_j = matrix_inner(world, bk.get_b(), vj);  // [nb,1]
        for (int a = 0; a < nb; ++a)
            for (int b = 0; b < nb; ++b)
                ref(a,b) += ov_i(a,0) * ov_j(b,0);
    }

    t.checkpoint((result - ref).normf(), 10.0 * thresh,
                 "project_Gab: CCPairFunction ket vs explicit Schwinger ref");

    return t.end();
}

// ---------------------------------------------------------------------------
// test 5: project_xy with Q_Gab bra — verifies Q + P = 1 for the Schwinger bra
// ---------------------------------------------------------------------------

template<std::size_t NDIM>
int test_QGUe(World& world) {
    static_assert(NDIM % 2 == 0, "NDIM must be even");
    constexpr std::size_t LDIM = NDIM / 2;
    test_output t("DiagnosticMatrix project_xy Q_Gab <NDIM=" + std::to_string(NDIM) + ">");

    const double thresh = FunctionDefaults<LDIM>::get_thresh();
    const double lo     = 1e-6;
    const double energy = -2.0;

    auto phi_i = FunctionFactory<double,LDIM>(world).f(gaussians::g10<LDIM>);
    auto phi_j = FunctionFactory<double,LDIM>(world).f(gaussians::g20<LDIM>);

    DiagnosticMatrix<double,NDIM> dm(world, std::vector<Function<double,LDIM>>{
        FunctionFactory<double,LDIM>(world).f(gaussians::g10<LDIM>),
        FunctionFactory<double,LDIM>(world).f(gaussians::g20<LDIM>),
        FunctionFactory<double,LDIM>(world).f(gaussians::g05<LDIM>),
    });
    const int nb = dm.nbasis();

    std::vector<Function<double,LDIM>> vi = {phi_i}, vj = {phi_j};
    auto ket = std::vector<CCPairFunction<double,NDIM>>{ CCPairFunction<double,NDIM>(vi, vj) };

    // Projector: P_1 = |proj><proj|  where proj = ao_basis[0] (orthonormal after ctor)
    //            Q_1 = 1 - P_1
    // Decompose the Gab bra into Q and P parts and verify Q + P = 1.
    const std::vector<Function<double,LDIM>> proj_vec = { dm.ao_basis[0] };

    auto full_bra = dm.build_Gab_bra(energy, lo);
    std::vector<CCPairFunction<double,NDIM>> Q_bra, P_bra;

    for (const auto& bk : full_bra) {
        auto p1 = bk.get_a();   // copy: particle-1 functions w_n·g_n·ao[a]
        auto p2 = bk.get_b();   // copy: particle-2 functions g_n·ao[b]
        // ov(0,a) = <proj | p1[a]>
        Tensor<double> ov = matrix_inner(world, proj_vec, p1);  // [1, nb]
        std::vector<Function<double,LDIM>> Qp1(nb), Pp1(nb);
        for (int a = 0; a < nb; ++a) {
            Pp1[a] = ov(0,a) * proj_vec[0];       // P_1 p1[a]  = <proj|p1[a]> · proj
            Qp1[a] = p1[a] - ov(0,a) * proj_vec[0]; // Q_1 p1[a]  = p1[a] - P_1 p1[a]
        }
        P_bra.emplace_back(Pp1, p2);
        Q_bra.emplace_back(Qp1, p2);
    }

    // Verify Q_1 + P_1 = 1:  project_xy(Q_bra) + project_xy(P_bra) = project_Gab
    auto result_full = dm.project_Gab(ket, energy, lo);
    auto result_Q    = dm.project_xy(Q_bra, ket);
    auto result_P    = dm.project_xy(P_bra, ket);

    t.checkpoint((result_full - result_Q - result_P).normf(), 10.0 * thresh,
                 "Q + P = 1: project_Gab = project_xy(Q_bra) + project_xy(P_bra)");

    return t.end();
}

// ---------------------------------------------------------------------------
// test 6: diagnose_GUe — apply G explicitly, check Schwinger ref vs 6D result
// ---------------------------------------------------------------------------

template<std::size_t NDIM>
int test_diagnose_GUe(World& world) {
    static_assert(NDIM % 2 == 0, "NDIM must be even");
    constexpr std::size_t LDIM = NDIM / 2;

    if constexpr (LDIM == 3) {
        test_output t("DiagnosticMatrix diagnose_GUe <NDIM=" + std::to_string(NDIM) + ">");

        const double thresh3 = FunctionDefaults<3>::get_thresh();
        const double thresh6 = FunctionDefaults<6>::get_thresh();
        const double lo      = 1e-6;
        const double energy  = -2.0;

        auto phi_i = FunctionFactory<double,3>(world).f(gaussians::g10<3>);
        auto phi_j = FunctionFactory<double,3>(world).f(gaussians::g20<3>);
        std::vector<Function<double,3>> ao_raw = {
            FunctionFactory<double,3>(world).f(gaussians::g10<3>),
            FunctionFactory<double,3>(world).f(gaussians::g20<3>),
            FunctionFactory<double,3>(world).f(gaussians::g05<3>),
        };

        CorrelationFactor cf(world, 1.0, 1e-10, lo);
        real_convolution_6d op_mod = BSHOperator<6>(world, sqrt(-2.0*energy), lo, 1.e-6);
        op_mod.modified() = true;

        auto Uphi_local     = cf.apply_U_local    (phi_i, phi_j, op_mod, thresh6);
        auto Uphi_semilocal = cf.apply_U_semilocal(phi_i, phi_j, op_mod, thresh6);

        // Apply G = (−½∇₁² − ½∇₂² + μ²/2)^{-1} to get GU|ij> as a 6D function
        real_convolution_6d G = BSHOperator<6>(world, sqrt(-2.0*energy), lo, 1.e-6);
        auto GUphi_local     = apply(G, Uphi_local).truncate();
        auto GUphi_semilocal = apply(G, Uphi_semilocal).truncate();

        auto dm = cf.diagnose_GUe(GUphi_local, GUphi_semilocal, phi_i, phi_j, ao_raw, energy);
        dm.print_report("test_diagnose_GUe");

        // The Schwinger quadrature introduces ~1% relative error beyond the 6D grid threshold.
        t.checkpoint(dm.entries["Glocal"].error,     10.0*thresh6, "diagnose_GUe: Glocal error");
        t.checkpoint(dm.entries["Gsemilocal"].error, 10.0*thresh6, "diagnose_GUe: Gsemilocal error");

        return t.end();
    }
    return 0;
}

// ---------------------------------------------------------------------------
// test 7: diagnose_QGUe — verify Q₁+P₁=1 for the Schwinger bra with Ue ket
// ---------------------------------------------------------------------------

template<std::size_t NDIM>
int test_diagnose_QGUe(World& world) {
    static_assert(NDIM % 2 == 0, "NDIM must be even");
    constexpr std::size_t LDIM = NDIM / 2;

    if constexpr (LDIM == 3) {
        test_output t("DiagnosticMatrix diagnose_QGUe Q+P=1 <NDIM=" + std::to_string(NDIM) + ">");

        const double thresh3 = FunctionDefaults<3>::get_thresh();
        const double thresh6 = FunctionDefaults<6>::get_thresh();
        const double lo      = 1e-6;
        const double energy  = -2.0;

        auto phi_i = FunctionFactory<double,3>(world).f(gaussians::g10<3>);
        auto phi_j = FunctionFactory<double,3>(world).f(gaussians::g20<3>);
        std::vector<Function<double,3>> ao_raw = {
            FunctionFactory<double,3>(world).f(gaussians::g10<3>),
            FunctionFactory<double,3>(world).f(gaussians::g20<3>),
            FunctionFactory<double,3>(world).f(gaussians::g05<3>),
        };

        CorrelationFactor cf(world, 1.0, 1e-10, lo);
        real_convolution_6d op_mod = BSHOperator<6>(world, sqrt(-2.0*energy), lo, 1.e-6);
        op_mod.modified() = true;

        auto Uphi_local     = cf.apply_U_local    (phi_i, phi_j, op_mod, thresh6);
        auto Uphi_semilocal = cf.apply_U_semilocal(phi_i, phi_j, op_mod, thresh6);

        DiagnosticMatrix<> dm(world, ao_raw);
        const int nb = dm.nbasis();

        // "Occupied" space for Q₁ = 1 − P₁ with P₁ = Σₖ|φₖ⟩⟨φₖ|.
        // phi_i and phi_j are used as the projector functions.
        const std::vector<Function<double,3>> occ = {phi_i, phi_j};

        // Split the full Gab bra into Q₁ and P₁ parts (acting on particle-1 only).
        // Q_1 p1[a] = p1[a] - Σₖ <occ_k|p1[a]>·occ_k,  P_1 p1[a] = Σₖ <occ_k|p1[a]>·occ_k
        // Q_1 + P_1 = 1  =>  project_xy(Q_bra) + project_xy(P_bra) = project_Gab
        auto full_bra = dm.build_Gab_bra(energy, lo);
        std::vector<CCPairFunction<double,6>> Q_bra, P_bra;

        for (const auto& bk : full_bra) {
            auto p1 = bk.get_a();  // copy: w_n·g_n·ao[a]
            auto p2 = bk.get_b();  // copy: g_n·ao[b] (unchanged)
            auto Pp1 = zero_functions<double,3>(world, nb);

            for (const auto& ok : occ) {
                const std::vector<Function<double,3>> ok_vec = {ok};
                Tensor<double> ov = matrix_inner(world, ok_vec, p1);  // [1, nb]
                for (int a = 0; a < nb; ++a) {
                    Pp1[a] += ov(0,a) * ok;    // P₁ part
                    p1[a]  -= ov(0,a) * ok;    // Q₁ part: p1 -= P₁*p1
                }
            }
            Q_bra.emplace_back(p1,  p2);   // Q₁·p1 ⊗ p2
            P_bra.emplace_back(Pp1, p2);   // P₁·p1 ⊗ p2
        }

        // Verify Q₁ + P₁ = 1: project_xy(Q_bra) + project_xy(P_bra) = project_Gab
        auto ref_loc  = dm.project_Gab(Uphi_local,     energy, lo);
        auto ref_sl   = dm.project_Gab(Uphi_semilocal, energy, lo);
        auto resQ_loc = dm.project_xy(Q_bra, Uphi_local);
        auto resP_loc = dm.project_xy(P_bra, Uphi_local);
        auto resQ_sl  = dm.project_xy(Q_bra, Uphi_semilocal);
        auto resP_sl  = dm.project_xy(P_bra, Uphi_semilocal);

        t.checkpoint((ref_loc - resQ_loc - resP_loc).normf(), thresh6, "Q+P=1: local Ue ket");
        t.checkpoint((ref_sl  - resQ_sl  - resP_sl ).normf(), thresh6, "Q+P=1: semilocal Ue ket");

        return t.end();
    }
    return 0;
}

// ---------------------------------------------------------------------------
// main
// ---------------------------------------------------------------------------

int main(int argc, char** argv) {
    World& world = madness::initialize(argc, argv);
    startup(world, argc, argv);
    commandlineparser parser(argc, argv);

    const int    k       = 5;
    const double thresh  = 1e-4;

    // Set up FunctionDefaults for all relevant dimensions
    FunctionDefaults<1>::set_k(k);  FunctionDefaults<1>::set_thresh(thresh);
    FunctionDefaults<1>::set_cubic_cell(-10., 10.);
    FunctionDefaults<2>::set_k(k);  FunctionDefaults<2>::set_thresh(thresh);
    FunctionDefaults<2>::set_cubic_cell(-10., 10.);
    FunctionDefaults<3>::set_k(k);  FunctionDefaults<3>::set_thresh(thresh);
    FunctionDefaults<3>::set_cubic_cell(-10., 10.);
    FunctionDefaults<4>::set_k(k);  FunctionDefaults<4>::set_thresh(thresh);
    FunctionDefaults<4>::set_cubic_cell(-10., 10.);
    FunctionDefaults<6>::set_k(k);  FunctionDefaults<6>::set_thresh(1e-3);
    FunctionDefaults<6>::set_tensor_type(TT_2D);
    FunctionDefaults<6>::set_cubic_cell(-10., 10.);

    if (world.rank() == 0)
        print("k =", k, "  thresh =", thresh);

    int result = 0;
#ifdef ENABLE_GENTENSOR
    try {
        result += test_ctor_orthonormalization<2>(world);
//        result += test_ctor_orthonormalization<4>(world);
//        result += test_ctor_orthonormalization<6>(world);
//        result += test_project_ab<2>(world);
//        result += test_project_ab<4>(world);
//        // result += test_project_ab<6>(world);  // slow: constructs 6D grid function
//        result += test_GUe<4>(world);
//        result += test_GUe<6>(world);
//        result += test_QGUe<4>(world);
//        result += test_QGUe<6>(world);
    } catch (std::exception& e) {
        print("test_operator_diagnostics: exception:", e.what());
        result = 1;
    }
//    try {
//        result += test_diagnose_Ue<6>(world);
//    } catch (madness::MadnessException& e) {
//        print("test_diagnose_Ue: MadnessException:", e.what());
//        result = 1;
//    } catch (std::exception& e) {
//        print("test_diagnose_Ue: exception:", e.what());
//        result = 1;
//    }
    try {
        result += test_diagnose_GUe<6>(world);
    } catch (madness::MadnessException& e) {
        print("test_diagnose_GUe: MadnessException:", e.what());
        result = 1;
    } catch (std::exception& e) {
        print("test_diagnose_GUe: exception:", e.what());
        result = 1;
    }
    try {
        result += test_diagnose_QGUe<6>(world);
    } catch (madness::MadnessException& e) {
        print("test_diagnose_QGUe: MadnessException:", e.what());
        result = 1;
    } catch (std::exception& e) {
        print("test_diagnose_QGUe: exception:", e.what());
        result = 1;
    }
    world.gop.fence();
#else
    if (world.rank() == 0)
        print("test_operator_diagnostics: needs -DENABLE_GENTENSOR=ON, skipping");
#endif

    finalize();
    return result;
}
