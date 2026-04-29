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
        bool   assemble_fK            = true;  ///< build the fK piece (set false to test Kf alone)
        bool   include_symmetry_mirror = true; ///< apply swap_particles for K̂₂ (set false to test K̂₁ alone)
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

    /// LRF algorithm (canonical): LRF-project f₁₂|ij⟩, apply K̂ to each factor.
    /// Thin wrapper over CCPotentials::apply_KffK_low_rank.
    static KffKResult apply_KffK_lowrank(
            World& world,
            const CCFunction<double, 3>& phi_i,
            const CCFunction<double, 3>& phi_j,
            const Info& info,
            const LowRankFunctionParameters& lrfparam,
            const real_convolution_6d* Gscreen = nullptr);

    /// LRF algorithm (per-k direct): loops over occupied orbitals, projects
    /// f₁₂·k_bra(1)·phi_j(2) etc. per k.  Thin wrapper over
    /// CCPotentials::apply_KffK_low_rank_direct.
    static KffKResult apply_KffK_lowrank_direct(
            World& world,
            const CCFunction<double, 3>& phi_i,
            const CCFunction<double, 3>& phi_j,
            const Info& info,
            const LowRankFunctionParameters& lrfparam,
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
};

} // namespace madness

#endif // MADNESS_CHEM_EXCHANGE_COMMUTATOR_H
