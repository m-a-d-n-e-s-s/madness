// operator_diagnostics.cc — implementation of DiagnosticMatrix<T,NDIM>.
//
// Include order matters: CCStructures.h → lowrankfunction.h → electronic_correlation_factor.h
// → operator_diagnostics.h defines DiagnosticMatrix, then CCStructures.h defines
// CCConvolutionOperator.  CCPairFunction instantiation in this TU therefore has
// CCConvolutionOperator fully defined.
//
// All member function templates are defined here and explicitly instantiated at
// the bottom.  NDIM is the pair-space dimension; LDIM = NDIM/2 is the one-particle space.

#include <madness/chem/CCStructures.h>   // defines CCConvolutionOperator before CCPairFunction is instantiated
#include <madness/chem/operator_diagnostics.h>

#include <madness/mra/operator.h>
#include <cmath>
#include <iomanip>

namespace madness {

// ---------------------------------------------------------------------------
//  DiagnosticMatrix constructor
// ---------------------------------------------------------------------------

template<typename T, std::size_t NDIM>
DiagnosticMatrix<T,NDIM>::DiagnosticMatrix(World& world, std::vector<Function<T,LDIM>> ao_basis_in,
                                           bool orthonormalize)
    : world_(world)
{
    if (ao_basis_in.empty()) return;

    if (!orthonormalize) {
        ao_basis = std::move(ao_basis_in);
        return;
    }

    const int nb = static_cast<int>(ao_basis_in.size());
    Tensor<T> S = matrix_inner(world, ao_basis_in, ao_basis_in, /*sym=*/true);

    // measure deviation from identity
    Tensor<T> diff = copy(S);
    for (int i = 0; i < nb; ++i) diff(i, i) -= 1.0;
    const double dev = diff.normf();

    if (dev > 1e-10) {
        if (world.rank() == 0)
            print("DiagnosticMatrix: ao_basis not orthonormal (||S-I||_F =", dev,
                  "), applying Löwdin orthonormalization");
        ao_basis = orthonormalize_symmetric(ao_basis_in, S);
    } else {
        ao_basis = std::move(ao_basis_in);
    }
}

// ---------------------------------------------------------------------------
//  project_xy — core primitive
// ---------------------------------------------------------------------------

template<typename T, std::size_t NDIM>
Tensor<T> DiagnosticMatrix<T,NDIM>::project_xy(const std::vector<CCPairFunction<T,NDIM>>& bra,
                                                 const Function<T,NDIM>& ket) const
{
    const int nb = nbasis();
    const int nk = static_cast<int>(bra.size());
    Tensor<T> result(nb, nb);

    // Step 1: partial inner products via Function::project_out (particle 0 = particle 1).
    std::vector<std::vector<Function<T,LDIM>>> ptmp(
        nk, std::vector<Function<T,LDIM>>(nb));
    for (int k = 0; k < nk; ++k) {
        const auto& gk = bra[k].get_a();
        for (int a = 0; a < nb; ++a)
            ptmp[k][a] = ket.project_out(gk[a], 0);
    }

    // Step 2: accumulate — bra-slots k innermost
    for (int a = 0; a < nb; ++a)
        for (int b = 0; b < nb; ++b)
            for (int k = 0; k < nk; ++k)
                result(a, b) += inner(ptmp[k][a], bra[k].get_b()[b]);

    return result;
}

template<typename T, std::size_t NDIM>
Tensor<T> DiagnosticMatrix<T,NDIM>::project_xy(const std::vector<CCPairFunction<T,NDIM>>& bra,
                                                 const std::vector<CCPairFunction<T,NDIM>>& ket) const
{
    const int nb = nbasis();
    const int nk = static_cast<int>(bra.size());
    Tensor<T> result(nb, nb);
    auto p1 = particle<LDIM>::particle1();

    for (const auto& f : ket) {
        if (!f.is_assigned()) continue;

        // Step 1: partial inner products for all (bra-slot k, basis index a)
        std::vector<std::vector<Function<T,LDIM>>> ptmp(
            nk, std::vector<Function<T,LDIM>>(nb));
        for (int k = 0; k < nk; ++k) {
            const auto& gk = bra[k].get_a();
            for (int a = 0; a < nb; ++a)
                ptmp[k][a] = inner(f, gk[a], p1.get_array(), p1.get_array());
        }

        // Step 2: accumulate — bra-slots k innermost
        for (int a = 0; a < nb; ++a)
            for (int b = 0; b < nb; ++b)
                for (int k = 0; k < nk; ++k)
                    result(a, b) += inner(ptmp[k][a], bra[k].get_b()[b]);
    }
    return result;
}

// ---------------------------------------------------------------------------
//  project_ab — thin wrappers
// ---------------------------------------------------------------------------

template<typename T, std::size_t NDIM>
Tensor<T> DiagnosticMatrix<T,NDIM>::project_ab(const Function<T,NDIM>& ket) const
{
    CCPairFunction<T,NDIM> bra_cc(ao_basis, ao_basis);
    return project_xy({bra_cc}, ket);
}

template<typename T, std::size_t NDIM>
Tensor<T> DiagnosticMatrix<T,NDIM>::project_ab(
        const std::vector<CCPairFunction<T,NDIM>>& ket) const
{
    return project_xy({CCPairFunction<T,NDIM>(ao_basis, ao_basis)}, ket);
}

// ---------------------------------------------------------------------------
//  build_Gab_bra — Schwinger BSH bra builder
// ---------------------------------------------------------------------------

template<typename T, std::size_t NDIM>
std::vector<CCPairFunction<T,NDIM>>
DiagnosticMatrix<T,NDIM>::build_Gab_bra(double energy, double lo) const
{
    MADNESS_CHECK_THROW(energy < 0.0, "build_Gab_bra: energy must be negative");
    const double thresh = FunctionDefaults<LDIM>::get_thresh();
    const double mu     = std::sqrt(-2.0 * energy);
    const double hi     = FunctionDefaults<LDIM>::get_cell_width().normf();

    // Schwinger BSH quadrature: G^NDIM(R) ≈ Σ_n c_n exp(-α_n R²)
    // with c_n = c3d_n * (α_n/π)^{LDIM/2}, where c3d_n is from the LDIM-D BSH Gaussian fit.
    auto fit   = GFit<T,LDIM>::BSHFit(mu, lo, hi, thresh);
    auto c3d   = fit.coeffs();
    auto alpha = fit.exponents();
    const int nfit = c3d.dim(0);

    std::vector<CCPairFunction<T,NDIM>> bra;
    bra.reserve(nfit);

    for (int n = 0; n < nfit; ++n) {
        const double an  = alpha[n];
        const double wNd = c3d[n] * std::pow(an / constants::pi, static_cast<double>(LDIM) / 2.0);

        auto gauss = SeparatedConvolution<T,LDIM>(
                world_, OperatorInfo(an, lo, thresh, OT_GAUSS));

        // conv[a] = g_n * ao_basis[a]  (Gaussian-convolved basis)
        auto conv = gauss(ao_basis);

        // wconv = wNd * conv  (weight absorbed into particle-1 factors)
        // Deep copy before scaling: MADNESS Function uses shallow-copy semantics,
        // so auto wconv = conv; scale(wconv) would also scale conv via the shared impl.
        auto wconv = copy(world_, conv);
        scale(world_, wconv, wNd);

        bra.emplace_back(wconv, conv);  // g=wconv (weighted), h=conv (unweighted)
    }
    return bra;
}

// ---------------------------------------------------------------------------
//  project_Gab — thin wrappers
// ---------------------------------------------------------------------------

template<typename T, std::size_t NDIM>
Tensor<T> DiagnosticMatrix<T,NDIM>::project_Gab(const Function<T,NDIM>& ket,
                                                   double energy, double lo) const
{
    return project_xy(build_Gab_bra(energy, lo), ket);
}

template<typename T, std::size_t NDIM>
Tensor<T> DiagnosticMatrix<T,NDIM>::project_Gab(
        const std::vector<CCPairFunction<T,NDIM>>& ket,
        double energy, double lo) const
{
    return project_xy(build_Gab_bra(energy, lo), ket);
}

// ---------------------------------------------------------------------------
//  compute_errors
// ---------------------------------------------------------------------------

template<typename T, std::size_t NDIM>
void DiagnosticMatrix<T,NDIM>::compute_errors()
{
    for (auto& [name, e] : entries)
        e.compute_error();
}

// ---------------------------------------------------------------------------
//  print_report
// ---------------------------------------------------------------------------

template<typename T, std::size_t NDIM>
void DiagnosticMatrix<T,NDIM>::print_report(const std::string& tag) const
{
    const std::string prefix = tag.empty() ? "[diag]" : "[diag/" + tag + "]";
    for (const auto& [name, e] : entries) {
        char buf[256];
        std::snprintf(buf, sizeof(buf),
                      "%s  %-20s  ||ref||=%10.3e  ||result||=%10.3e  error=%10.3e",
                      prefix.c_str(), name.c_str(),
                      e.ref.normf(), e.result.normf(), e.error);
        print(std::string(buf));
    }
    if (time > 0.0) print(prefix, "  time =", time, "s");
}

// ---------------------------------------------------------------------------
//  Explicit instantiations — NDIM is the pair-space dimension
// ---------------------------------------------------------------------------

template class DiagnosticMatrix<double, 2>;  // LDIM=1 one-particle
template class DiagnosticMatrix<double, 4>;  // LDIM=2 one-particle
template class DiagnosticMatrix<double, 6>;  // LDIM=3 one-particle (standard 3D chemistry)

} // namespace madness
