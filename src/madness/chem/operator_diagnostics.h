// operator_diagnostics.h — unified <ab|op|ij> diagnostic matrix container.
//
// Provides a map-based container (DiagnosticMatrix) and projection utilities:
//   project_xy   — core primitive: result(a,b) += sum_k <bra[k].g[a](1) bra[k].h[b](2)|ket>
//   project_ab   — thin wrapper: bra = {CCPairFunction(ao_basis, ao_basis)}
//   project_Gab  — thin wrapper: builds nfit Schwinger bra slots, calls project_xy
//
// Convention: `ref`    = 3D analytic reference, set by the caller.
//             `result` = 6D projection, set by the methods below.
//             `error`  = ||ref - result||_F, computed by compute_errors().
//
// Bra representation: each element of the bra vector is a decomposed
// CCPairFunction<double,6> where get_a()[a] is the particle-1 factor for basis
// index a and get_b()[b] is the particle-2 factor for basis index b.  Weights
// are absorbed into the get_a() functions by the bra builder.  The bra type is
// the same as the ket type, so apply(ProjectorBase, bra) works directly.

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

struct DiagnosticEntry {
    Tensor<double> ref;      ///< 3D analytic reference (set by caller)
    Tensor<double> result;   ///< 6D projection result (set by DiagnosticMatrix methods)
    double error = -1.0;

    void compute_error() {
        if (ref.size() > 0 && result.size() > 0)
            error = (ref - result).normf();
    }
};

/// Container for named <ab|op|ij> diagnostic matrix elements.
///
/// Each named entry holds a pair (ref, result) of nbasis×nbasis tensors and
/// their Frobenius-norm difference.  The caller fills `ref` from 3D analytic
/// formulas; the project_ab / project_Gab methods fill `result` by projecting
/// 6D functions onto the AO observer basis pairs |ab>.
class DiagnosticMatrix {
public:
    std::map<std::string, DiagnosticEntry> entries;
    std::vector<Function<double,3>> ao_basis;  ///< orthonormalized observer basis
    double time = 0.0;

    /// Construct with a World reference and an orthonormalized AO basis.
    DiagnosticMatrix(World& world, std::vector<Function<double,3>> ao_basis_in)
        : ao_basis(std::move(ao_basis_in)), world_(world) {}

    int nbasis() const { return static_cast<int>(ao_basis.size()); }

    /// Allocate zero tensors for a new named quantity.
    /// Returns the ref tensor by reference so the caller can fill it in-place.
    Tensor<double>& init(const std::string& name) {
        auto& e = entries[name];
        e.ref    = Tensor<double>(nbasis(), nbasis());
        e.result = Tensor<double>(nbasis(), nbasis());
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
    Tensor<double> project_xy(const real_function_6d& ket,
                               const std::vector<CCPairFunction<double,6>>& bra) const;
    Tensor<double> project_xy(const std::vector<CCPairFunction<double,6>>& ket,
                               const std::vector<CCPairFunction<double,6>>& bra) const;

    // -------------------------------------------------------------------------
    // Convenience wrappers
    // -------------------------------------------------------------------------

    /// Returns T(a,b) = <ao_basis[a] ao_basis[b] | ket>.
    Tensor<double> project_ab(const real_function_6d& ket) const;
    Tensor<double> project_ab(const std::vector<CCPairFunction<double,6>>& ket) const;

    /// Returns T(a,b) ≈ <ao_basis[a] ao_basis[b] | G · ket> via Schwinger BSH quadrature.
    Tensor<double> project_Gab(const real_function_6d& ket, double energy, double lo) const;
    Tensor<double> project_Gab(const std::vector<CCPairFunction<double,6>>& ket,
                                double energy, double lo) const;

    // -------------------------------------------------------------------------
    // Post-projection utilities
    // -------------------------------------------------------------------------

    /// Compute error = ||ref - result||_F for every initialized entry.
    void compute_errors();

    /// Print a compact report: name, ||ref||, ||result||, error.
    void print_report(const std::string& tag = "") const;

private:
    World& world_;

    /// Build nfit bra CCPairFunctions for the Schwinger approximation of G (BSH fit for energy).
    /// bra[n].get_a()[a] = w_n^{6d} * (g_n * ao_basis[a])  (weight absorbed)
    /// bra[n].get_b()[b] =             g_n * ao_basis[b]
    std::vector<CCPairFunction<double,6>> build_Gab_bra(double energy, double lo) const;
};

} // namespace madness

#endif // MADNESS_CHEM_OPERATOR_DIAGNOSTICS_H
