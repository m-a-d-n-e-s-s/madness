// operator_diagnostics.h — unified <ab|op|ij> diagnostic matrix container.
//
// Provides a map-based container (DiagnosticMatrix) and projection utilities:
//   project_xy   — core primitive: result(a,b) += sum_k <bra[k].g[a](1) bra[k].h[b](2)|ket>
//                  signature: project_xy(bra, ket)
//   project_ab   — thin wrapper: bra = {CCPairFunction(ao_basis, ao_basis)}
//   project_Gab  — thin wrapper: builds nfit Schwinger bra slots, calls project_xy
//
// Convention: `ref`    = LDIM-D analytic reference, set by the caller.
//             `result` = NDIM-D projection, set by the methods below.
//             `error`  = ||ref - result||_F, computed by compute_errors().
//
// Template parameters:
//   T    — scalar type (default double)
//   NDIM — pair-space dimension (default 6); one-particle space is LDIM = NDIM/2.
//
// Bra representation: each element of the bra vector is a decomposed
// CCPairFunction<T,NDIM> where get_a()[a] is the particle-1 factor for basis
// index a and get_b()[b] is the particle-2 factor for basis index b.  Weights
// are absorbed into the get_a() functions by the bra builder.

#ifndef MADNESS_CHEM_OPERATOR_DIAGNOSTICS_H
#define MADNESS_CHEM_OPERATOR_DIAGNOSTICS_H

#include <map>
#include <string>
#include <vector>

#include <madness/mra/mra.h>
#include <madness/chem/ccpairfunction.h>
// lowrankfunction.h and operator.h are included only in operator_diagnostics.cc
// to avoid a circular include via lowrankfunction.h -> electronic_correlation_factor.h.

namespace madness {

template<typename T = double>
struct DiagnosticEntry {
    Tensor<T> ref;      ///< LDIM-D analytic reference (set by caller)
    Tensor<T> result;   ///< NDIM-D projection result (set by DiagnosticMatrix methods)
    double error = -1.0;

    void compute_error() {
        if (ref.size() > 0 && result.size() > 0)
            error = (ref - result).normf();
    }
};

/// Container for named <ab|op|ij> diagnostic matrix elements.
///
/// Each named entry holds a pair (ref, result) of nbasis×nbasis tensors and
/// their Frobenius-norm difference.  The caller fills `ref` from LDIM-D analytic
/// formulas; the project_ab / project_Gab methods fill `result` by projecting
/// NDIM-D functions onto the AO observer basis pairs |ab>.
///
/// @tparam T    scalar type (default double)
/// @tparam NDIM pair-space dimension (default 6); LDIM = NDIM/2 is the one-particle space
template<typename T = double, std::size_t NDIM = 6>
class DiagnosticMatrix {
public:
    static constexpr std::size_t LDIM = NDIM / 2;  ///< one-particle space dimension

    /// element provider: M(x,y) = <p1[x](1) p2[y](2) | X | ij> for arbitrary
    /// LDIM-D bra lists p1, p2 — the operator-specific unit plugged into
    /// ref_Gab / ref_GQab.  Build via the provider factories of the operator
    /// classes (CorrelationFactor::ue_*_provider, ExchangeCommutator::*_provider),
    /// which hoist convolution operators and intermediates into the closure.
    using ElementProvider = std::function<Tensor<T>(
            const std::vector<Function<T,LDIM>>&,
            const std::vector<Function<T,LDIM>>&)>;

    std::map<std::string, DiagnosticEntry<T>> entries;
    std::vector<Function<T,LDIM>> ao_basis;  ///< orthonormalized observer basis
    double time = 0.0;

    /// Construct with a World reference and an AO basis.
    /// If the basis is not orthonormal (||S - I||_F > 1e-10) it is
    /// orthonormalized in-place via symmetric (Löwdin) orthonormalization.
    /// Pass orthonormalize=false to keep the basis as given (e.g. when matrix
    /// elements in the raw basis carry physical meaning, like pair energies).
    DiagnosticMatrix(World& world, std::vector<Function<T,LDIM>> ao_basis_in,
                     bool orthonormalize = true);

    int nbasis() const { return static_cast<int>(ao_basis.size()); }

    /// Allocate zero tensors for a new named quantity.
    /// Returns the ref tensor by reference so the caller can fill it in-place.
    Tensor<T>& init(const std::string& name) {
        auto& e = entries[name];
        e.ref    = Tensor<T>(nbasis(), nbasis());
        e.result = Tensor<T>(nbasis(), nbasis());
        e.error  = -1.0;
        return e.ref;
    }

    // -------------------------------------------------------------------------
    // Core projection primitive
    // -------------------------------------------------------------------------

    /// Core primitive: returns T(a,b) = Σ_k inner(bra[k].get_a()[a](1) bra[k].get_b()[b](2), ket_piece)
    /// summed over all assigned ket pieces and all bra slots k.
    /// Weight is absorbed into bra[k].get_a() by the bra builder (project_Gab etc.).
    /// Loop order: ket-pieces outer, (a,b) middle, bra-slots k innermost.
    /// Caller assigns the returned tensor to entries[name].result.
    Tensor<T> project_xy(const std::vector<CCPairFunction<T,NDIM>>& bra,
                          const Function<T,NDIM>& ket) const;
    Tensor<T> project_xy(const std::vector<CCPairFunction<T,NDIM>>& bra,
                          const std::vector<CCPairFunction<T,NDIM>>& ket) const;

    // -------------------------------------------------------------------------
    // Convenience wrappers
    // -------------------------------------------------------------------------

    /// Returns T(a,b) = <ao_basis[a] ao_basis[b] | ket>.
    Tensor<T> project_ab(const Function<T,NDIM>& ket) const;
    Tensor<T> project_ab(const std::vector<CCPairFunction<T,NDIM>>& ket) const;

    /// Returns T(a,b) ≈ <ao_basis[a] ao_basis[b] | G · ket> via Schwinger BSH quadrature.
    Tensor<T> project_Gab(const Function<T,NDIM>& ket, double energy, double lo) const;
    Tensor<T> project_Gab(const std::vector<CCPairFunction<T,NDIM>>& ket,
                           double energy, double lo) const;

    // -------------------------------------------------------------------------
    // 3D-only reference builders — Schwinger fit of G on the bra
    // -------------------------------------------------------------------------

    /// ref(a,b) = Σ_n w_n <ã_n b̃_n | X | ij>  with ã_n = w_n g_{α_n}⋆a, b̃_n = g_{α_n}⋆b
    /// (Schwinger/BSH fit of G moved onto the bra; X supplied as an element provider).
    Tensor<T> ref_Gab(const ElementProvider& elements, double energy, double lo) const;

    /// ref(a,b) = Σ_n w_n <ã_n b̃_n | Q₁₂ X | ij>  with Q₁₂ = (1-O₁)(1-O₂),
    /// O = Σ_k |k_ket><k_bra| (plain sum, matching StrongOrthogonalityProjector).
    /// Q₁₂ is expanded on the ket side as scalar contractions: per slot
    ///   A − S1·B − C·S2ᵀ + S1·D·S2ᵀ
    /// with element blocks over p1 = {ã_a} ∪ {occ_bra_k}, p2 = {b̃_b} ∪ {occ_bra_l}
    /// and overlaps S1(a,k) = <ã_a|occ_ket_k>, S2(b,l) = <b̃_b|occ_ket_l>.
    /// No nearly-vanishing functions like Q(g_n⋆a) are ever formed — applying Q
    /// to the bra functions suffers catastrophic cancellation when the convolved
    /// observers lie nearly in the occupied space (see operator_diagnostics.md).
    Tensor<T> ref_GQab(const ElementProvider& elements,
                       const std::vector<Function<T,LDIM>>& occ_ket,
                       const std::vector<Function<T,LDIM>>& occ_bra,
                       double energy, double lo) const;

    // -------------------------------------------------------------------------
    // Post-projection utilities
    // -------------------------------------------------------------------------

    /// Compute error = ||ref - result||_F for every initialized entry.
    void compute_errors();

    /// Print a compact report: name, ||ref||, ||result||, error.
    void print_report(const std::string& tag = "") const;

    // -------------------------------------------------------------------------
    // Bra builders — return bra vectors for use with project_xy
    // -------------------------------------------------------------------------

    /// Single-slot bra for plain <ab| projection: {CCPairFunction(ao_basis, ao_basis)}.
    std::vector<CCPairFunction<T,NDIM>> build_simple_bra() const {
        return { CCPairFunction<T,NDIM>(ao_basis, ao_basis) };
    }

    /// Multi-slot Schwinger bra for <ab|G| projection (BSH fit for energy).
    /// bra[n].get_a()[a] = w_n * (g_n * ao_basis[a])  (weight absorbed into particle 1)
    /// bra[n].get_b()[b] =        g_n * ao_basis[b]
    std::vector<CCPairFunction<T,NDIM>> build_Gab_bra(double energy, double lo) const;

private:
    World& world_;
};

// Deduction guides: infer NDIM=2*LDIM from the element type of the ao_basis vector.
template<typename T, std::size_t LDIM>
DiagnosticMatrix(World&, std::vector<Function<T,LDIM>>) -> DiagnosticMatrix<T,2*LDIM>;
template<typename T, std::size_t LDIM>
DiagnosticMatrix(World&, std::vector<Function<T,LDIM>>, bool) -> DiagnosticMatrix<T,2*LDIM>;

/// Print the (0,1) elements of the named pieces as "6d / 3d-ref / diff" lines
/// plus the signed total, via snprintf + print.
///
/// Intended for a DiagnosticMatrix built over a raw (non-orthonormalized)
/// pair-bra basis {bra_i, bra_j}, where entry(0,1) = <bra_i bra_j|G Q₁₂ X|ij>
/// is the contribution of piece X to the MP2 pair energy.
///
/// @param pieces entry names; @param signs per-piece sign in the total
template<typename T, std::size_t NDIM>
void print_pair_energy_report(const DiagnosticMatrix<T,NDIM>& dm,
                              const std::vector<std::string>& pieces,
                              const std::vector<double>& signs,
                              const std::string& title);

// forward declaration -- operator.h is included only in operator_diagnostics.cc
template<typename Q, std::size_t MDIM> class SeparatedConvolution;

/// Strong-orthogonality-projected pair-energy weights over an orthonormal
/// observer basis {a}:
///   W^Q(a,b) = 2<phi_i_bra (Q a)|g12|(Q b) phi_j_bra>
///                - <phi_j_bra (Q a)|g12|(Q b) phi_i_bra>
/// with Q = 1 - sum_m |occ_ket_m><occ_bra_m| applied to both observer factors
/// (Q12 = Q⊗Q is separable).  Used to contract <ab|G Q12 X|ij> diagnostic
/// tensors into the MP2 pair-energy contribution
///   E = factor * sum_ab W^Q(a,b) <ab|G Q12 X|ij>.
template<typename T, std::size_t NDIM>
Tensor<T> pair_energy_weights(World& world,
                              const std::vector<Function<T,NDIM/2>>& observer,
                              const Function<T,NDIM/2>& phi_i_bra,
                              const Function<T,NDIM/2>& phi_j_bra,
                              const std::vector<Function<T,NDIM/2>>& occ_ket,
                              const std::vector<Function<T,NDIM/2>>& occ_bra,
                              const SeparatedConvolution<T,NDIM/2>& g12);

/// RI pair-energy report: contract the .result (6d-projected ket, "Expr3") and
/// .ref (3d Schwinger ket, "Expr2") tensors of the named pieces with the
/// pair-energy weights W and print "Expr2(3d-ket) / Expr3(6d-ket) / diff" per
/// piece plus the signed total, where the per-term energy is
///   factor * sum_p sign_p * sum_ab W(a,b) M_p(a,b).
template<typename T, std::size_t NDIM>
void print_pair_energy_report_RI(const DiagnosticMatrix<T,NDIM>& dm,
                                 const std::vector<std::string>& pieces,
                                 const std::vector<double>& signs,
                                 const Tensor<T>& W,
                                 double factor,
                                 const std::string& title);

} // namespace madness

#endif // MADNESS_CHEM_OPERATOR_DIAGNOSTICS_H
