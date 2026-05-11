//
// exchange_commutator.h — library-style API for evaluating [K̂, f₁₂] |ij⟩
//
// Phase A: thin wrappers around existing CCPotentials::apply_Kfxy /
// apply_KffK_low_rank / apply_KffK_low_rank_direct, plus the scalar
// split-α diagnostic extracted from benchmark_exchange_commutator.cc.
// All entry points return a KffKResult with uniform wall-time / memory /
// rank instrumentation; all diagnostics run through a single diagnose()
// function operating on a harmonic-basis projection of the 6D result.
//

#ifndef MADNESS_CHEM_EXCHANGE_COMMUTATOR_H
#define MADNESS_CHEM_EXCHANGE_COMMUTATOR_H

#include <madness/chem/CCStructures.h>
#include <madness/chem/ccpairfunction.h>
#include <madness/chem/lowrankfunction.h>
#include <madness/chem/CCParameters.h>

namespace madness {

struct ExchangeCommutator {

    /// Orthonormalized AO basis used as the observer set in diagnose().
    /// The caller is responsible for orthonormalization; near-linear-
    /// dependent inputs will produce a numerically unreliable diagnostic.
    std::vector<Function<double, 3>> ao_basis;

    /// Construct with the AO basis used by diagnose() to project the
    /// commutator onto matrix elements ⟨ab | · | ij⟩.
    explicit ExchangeCommutator(std::vector<Function<double, 3>> ao_basis_in)
        : ao_basis(std::move(ao_basis_in)) {}

    /// Default ctor leaves ao_basis empty; only the static apply_KffK_*
    /// entry points may be used in that state — diagnose() will throw.
    ExchangeCommutator() = default;

    /// Uniform result envelope for every algorithm path.  Each "piece" may
    /// be empty when the algorithm only produces the combined commutator
    /// (e.g. lrf-split-alpha constructs KffK directly).
    struct KffKResult {
        std::vector<CCPairFunction<double, 6>> Kf;     ///<  K̂ f₁₂ | i j ⟩
        std::vector<CCPairFunction<double, 6>> fK;     ///<  f₁₂ K̂ | i j ⟩
        std::vector<CCPairFunction<double, 6>> KffK;   ///<  [K̂, f₁₂] | i j ⟩ = Kf − fK (signs baked in)
        double      t_wall = 0.0;                      ///< wall time (s)
        double      mem_gb = 0.0;                      ///< combined size of pair factors (GByte)
        long        rank   = 0;                        ///< combined LRF rank (0 for 6D)
        std::string algo;                              ///< "6d" / "lrf" / "lrf-direct" / "lrf-split-alpha"
    };

    /// Per-piece and combined errors obtained by projecting onto a harmonic basis.
    /// Reference matrices are indexed [a,b] over the harmonic-basis observer set
    /// (a on particle 1, b on particle 2).  For phi_i == phi_j the K̂₁/K̂₂
    /// contributions coincide; the four ref_*_K* fields make the asymmetric
    /// case (phi_i ≠ phi_j) directly inspectable.
    struct Diagnostics {
        double err_Kf   = -1.0;
        double err_fK   = -1.0;
        double err_KffK = -1.0;
        Tensor<double> ref_piece1;   ///< ⟨ab | K̂ f | ij⟩ = ref_Kf_K1 + ref_Kf_K2 (K̂₂ optional)
        Tensor<double> ref_piece2;   ///< ⟨ab | f K̂ | ij⟩ = ref_fK_K1 + ref_fK_K2 (K̂₂ optional)
        Tensor<double> ref_Kf_K1;    ///< ⟨ab | K̂₁ f | ij⟩
        Tensor<double> ref_Kf_K2;    ///< ⟨ab | K̂₂ f | ij⟩
        Tensor<double> ref_fK_K1;    ///< ⟨ab | f K̂₁ | ij⟩
        Tensor<double> ref_fK_K2;    ///< ⟨ab | f K̂₂ | ij⟩
        double t_reference = 0.0;    ///< wall time for the reference build (s)
    };

    /// Options for the split-α LRF k-commutator variant
    /// (see lrf_k_commutator_analysis.md and lrf_k_commutator_results.md).
    struct SplitAlphaOptions {
        double alpha_star = 1.0e4;   ///< GFit partition threshold; α_μ > α* is discarded
        double lo         = 1.0e-6;
        double hi         = 10.0;
        double eps_gfit   = 1.0e-6;
    };

    /// Options for the three-range LRF k-commutator variant.
    /// Partition the Coulomb GFit into:
    ///   * diffuse  (α_μ < alpha_lo): wide Gaussians, large-ve LRF, low rank
    ///   * medium   (alpha_lo ≤ α_μ ≤ alpha_hi): tighter LRF
    ///   * tight    (α_μ > alpha_hi): discard (δ-cancellation in commutator)
    /// See lrf_three_range_gfit.md for the derivation.
    struct ThreeRangeOptions {
        double alpha_lo  = 1.0;      ///< diffuse / medium boundary
        double alpha_hi  = 1.0e4;    ///< medium / tight boundary  (≡ alpha_star)
        double lo        = 1.0e-6;
        double hi        = 10.0;
        double eps_gfit  = 1.0e-6;
        /// Per-range LRF construction parameters; both default-constructed.
        /// Caller is expected to set radius / volume_element / tol / lmax /
        /// f12type per range before invocation.
        LowRankFunctionParameters lrfparam_diffuse;
        LowRankFunctionParameters lrfparam_medium;
        /// If true, build the medium-range LRF via project_from_operator
        /// (Taylor / direct construction) instead of random-Y project().
        /// Mandatory for the medium slab in isolation: random-Y stalls at
        /// large error without diffuse anchor functions.
        bool medium_use_taylor = true;
        /// Max Taylor order for project_from_operator (passed straight
        /// through; default matches lowrankfunction.h's own default).
        int  max_taylor_order  = 15;
        /// If true, the medium slab is treated by a *6D* apply with a
        /// partial-Coulomb operator (only the medium Gaussians) instead
        /// of any LRF.  The diffuse slab still uses the LRF path.
        /// Mutually exclusive with `medium_use_taylor`; takes precedence.
        bool medium_use_6d     = false;
    };

    // ---------------------------------------------------------------------
    // Top-level entry points — each produces a KffKResult wrapping one of
    // the four existing algorithm paths.
    // ---------------------------------------------------------------------

    /// 6D algorithm: Kf via apply_Kfxy, fK via make_f_xy_macrotask(Kx,...).
    /// Mirrors the body of test_kcomm_6d in the original benchmark file.
    static KffKResult apply_KffK_6d(
            World& world,
            const CCFunction<double, 3>& phi_i,
            const CCFunction<double, 3>& phi_j,
            const Info& info,
            const real_convolution_6d* Gscreen = nullptr);

    // ---------------------------------------------------------------------
    // Split-α k-commutator.  Scalar-only variant for Phase A — 6D pair
    // assembly is deferred.
    // ---------------------------------------------------------------------

    /// Assemble the 6D exchange commutator [K̂, f] |ij⟩ using the split-α
    /// LRF k-commutator algorithm.  Sums over occupied orbitals in
    /// info.mo_bra/info.mo_ket.  Currently assumes phi_i == phi_j
    /// (symmetric case) and exploits swap_particles for the K̂₂ piece.
    ///
    /// Output pair layout:
    ///   result.pair[0] = Σ_k Σ_ρ  g_kρ(1) · (phi_j · A_kρ)(2)       (separable)
    ///   result.pair[1] = Σ_k Σ_ρ  B_kρ · (g_kρ·phi_i)(1) · f(1,2) · phi_j(2)
    ///                    with coefficient −1
    ///   plus their particle-2 mirrors when i == j.
    static KffKResult apply_KffK_lowrank_split_alpha(
            World& world,
            const CCFunction<double, 3>& phi_i,
            const CCFunction<double, 3>& phi_j,
            const Info& info,
            const LowRankFunctionParameters& lrfparam,
            const SplitAlphaOptions& opt);

    /// Three-range LRF k-commutator: builds two LRFs (diffuse + medium) with
    /// independent LowRankFunctionParameters, concatenates their g/h vectors
    /// and assembles Kf / fK / KffK with the same algebra as
    /// apply_KffK_lowrank_split_alpha.  Tight terms (α_μ > opt.alpha_hi) are
    /// discarded.  Returns KffKResult with algo = "lrf-three-range" and rank
    /// = rank_diffuse + rank_medium (printed individually in t_wall summary).
    static KffKResult apply_KffK_lowrank_three_range(
            World& world,
            const CCFunction<double, 3>& phi_i,
            const CCFunction<double, 3>& phi_j,
            const Info& info,
            const ThreeRangeOptions& opt);

    // ---------------------------------------------------------------------
    // Diagnostics — harmonic-basis projection of the 6D result, compared
    // against a straight-6D reference built from Exchange<3> and f12 applied
    // to that basis.  Unassigned CCPairFunctions are skipped.
    // ---------------------------------------------------------------------

    /// Pair-vector arguments allow commutator pieces built from multiple
    /// CCPairFunctions (e.g. the separable + f12-wrapped pair from the
    /// split-α assembly).  An empty vector skips that role.  Set verbose
    /// to dump per-stage norms (harmonic basis, K applications, f12
    /// applications, reference matrices, and per-pair shape/norms) — use
    /// when a piece comes back as NaN to localize the source.
    ///
    /// The observer basis is the orthonormalized AO basis stored on this
    /// instance (`ao_basis`).  When `ao_basis` is empty the routine falls
    /// back to building a Cartesian-Gaussian harmonic basis from
    /// `obs_param` placed at `centers` (origin if empty), then
    /// orthonormalizing it canonically — preserving the original
    /// behaviour for callers that haven't supplied an AO basis.
    ///
    /// The reference is built as the sum of the K̂₁ and K̂₂ pieces (the
    /// latter is gated by `include_K2`).  All four ⟨ab|·|ij⟩ contributions
    /// are also returned individually in Diagnostics so non-symmetric
    /// regressions in either particle's K̂ application can be localized.
    Diagnostics diagnose(
            World& world,
            const std::vector<Function<double, 3>>& kvec,
            const std::vector<Function<double, 3>>& R2kvec,
            const Function<double, 3>& phi_i,
            const Function<double, 3>& phi_j,
            const std::vector<CCPairFunction<double, 6>>& Kf,
            const std::vector<CCPairFunction<double, 6>>& fK,
            const std::vector<CCPairFunction<double, 6>>& KffK,
            const LowRankFunctionParameters& obs_param,
            bool verbose = false,
            bool include_K2 = true,
            const std::vector<Vector<double, 3>>& centers = {}) const;

    /// Convenience: run diagnose() on a KffKResult treated as KffK and
    /// attach errors to stderr in a uniform one-line summary.
    static void print_report(
            const KffKResult& result,
            const Diagnostics* diag = nullptr);

    // ---------------------------------------------------------------------
    // G·[K̂,f] diagnostics — Schwinger-quadrature 3D reference vs 6D projection
    // ---------------------------------------------------------------------

    /// Per-piece errors for ⟨ab | (T−E)⁻¹ [K̂,f₁₂] | ij⟩.
    struct GKffKDiagnostics {
        Tensor<double> ref_GKf;       ///< ⟨ab|G Kf|ij⟩  from 6D projection
        Tensor<double> ref_GfK;       ///< ⟨ab|G fK|ij⟩  from 6D projection
        Tensor<double> ref_GKffK;     ///< ref_GKf − ref_GfK
        Tensor<double> result_GKf;    ///< ⟨ab|G Kf|ij⟩  from 3D Schwinger quadrature
        Tensor<double> result_GfK;    ///< ⟨ab|G fK|ij⟩  from 3D Schwinger quadrature
        Tensor<double> result_GKffK;  ///< result_GKf − result_GfK
        double error_GKf   = 0.0;
        double error_GfK   = 0.0;
        double error_GKffK = 0.0;
        double time        = 0.0;
    };

    /// Diagnose ⟨ab | G [K̂,f₁₂] | ij⟩ by comparing a 6D projection against
    /// a 3D Schwinger eigentime quadrature.
    ///
    /// The 6D references GKf_cc and GfK_cc are the G-applied commutator pieces
    /// as CCPairFunction vectors.  Accepting CCPairFunction (instead of a bare
    /// real_function_6d) allows any algorithm variant to be tested:
    ///
    ///   * 6D path (apply_KffK_6d): Kf and fK are pure-6D; caller applies G6d
    ///     via apply(G6d, kffk.Kf[0].get_function()) and wraps in a vector.
    ///   * LRF path (apply_KffK_lowrank_split_alpha, three_range): Kf may be
    ///     decomposed.  Caller converts to pure (sum outer products), applies G6d,
    ///     and passes the result wrapped in a CCPairFunction vector.
    ///   * Mixed/multi-piece results: sum all pieces (each entry is summed inside).
    ///
    /// The 6D reference is projected using the same partial_inner / matrix_inner
    /// machinery as diagnose(), handling both pure and decomposed CCPairFunctions.
    ///
    /// The 3D Schwinger formulas use only K̂, K̂†, and the Slater f₁₂:
    ///
    ///   G Kf:  ⟨ã b̃ | K̂₁ f₁₂ | φᵢ φⱼ⟩ = inner((K̂†ã)·φᵢ, f₁₂(b̃·φⱼ))  + K̂₂ swap
    ///   G fK:  ⟨ã b̃ | f₁₂ K̂₁ | φᵢ φⱼ⟩ = inner(ã·K̂φᵢ,    f₁₂(b̃·φⱼ))  + K̂₂ swap
    ///
    /// where ã = g_{αₙ}*a and b̃ = g_{αₙ}*b are Gaussian-smoothed AO basis
    /// functions and the sum over quadrature nodes n (with 6D weight
    /// w_n = c_n^{bsh}·(αₙ/π)^{3/2}) approximates G.
    ///
    /// @param GKf_cc   G·K̂f₁₂|ij⟩ pieces (caller applies G6d externally)
    /// @param GfK_cc   G·f₁₂K̂|ij⟩ pieces
    /// @param Kphi_i   K̂φᵢ (precomputed once by caller)
    /// @param Kphi_j   K̂φⱼ (precomputed once by caller)
    /// @param info     molecular/orbital info carrying mo_ket, mo_bra, parameters
    /// @param energy   ε_i + ε_j  (must be negative)
    GKffKDiagnostics diagnose_GKffK(
            World& world,
            const std::vector<CCPairFunction<double,6>>& GKf_cc,
            const std::vector<CCPairFunction<double,6>>& GfK_cc,
            const real_function_3d& phi_i,
            const real_function_3d& phi_j,
            const real_function_3d& Kphi_i,
            const real_function_3d& Kphi_j,
            const Info& info,
            double energy) const;
};

} // namespace madness

#endif // MADNESS_CHEM_EXCHANGE_COMMUTATOR_H
