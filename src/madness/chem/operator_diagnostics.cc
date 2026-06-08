// operator_diagnostics.cc — implementation of DiagnosticMatrix.
//
// Include order matters: CCStructures.h → lowrankfunction.h → electronic_correlation_factor.h
// → operator_diagnostics.h defines DiagnosticMatrix, then CCStructures.h defines
// CCConvolutionOperator.  CCPairFunction instantiation in this TU therefore has
// CCConvolutionOperator fully defined.

#include <madness/chem/CCStructures.h>   // defines CCConvolutionOperator before CCPairFunction is instantiated
#include <madness/chem/operator_diagnostics.h>

#include <madness/mra/operator.h>
#include <cmath>
#include <iomanip>

namespace madness {

// ---------------------------------------------------------------------------
//  project_xy — core primitive
// ---------------------------------------------------------------------------

Tensor<double> DiagnosticMatrix::project_xy(const real_function_6d& ket,
                                              const std::vector<CCPairFunction<double,6>>& bra) const
{
    const int nb = nbasis();
    const int nk = static_cast<int>(bra.size());
    Tensor<double> result(nb, nb);

    // Step 1: partial inner products via Function::project_out (particle 0 = particle 1).
    std::vector<std::vector<Function<double,3>>> ptmp(
        nk, std::vector<Function<double,3>>(nb));
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

Tensor<double> DiagnosticMatrix::project_xy(const std::vector<CCPairFunction<double,6>>& ket,
                                              const std::vector<CCPairFunction<double,6>>& bra) const
{
    const int nb = nbasis();
    const int nk = static_cast<int>(bra.size());
    Tensor<double> result(nb, nb);
    auto p1 = particle<3>::particle1();

    for (const auto& f : ket) {
        if (!f.is_assigned()) continue;

        // Step 1: partial inner products for all (bra-slot k, basis index a)
        std::vector<std::vector<Function<double,3>>> ptmp(
            nk, std::vector<Function<double,3>>(nb));
        for (int k = 0; k < nk; ++k) {
            const auto& gk = bra[k].get_a();
            for (int a = 0; a < nb; ++a)
                ptmp[k][a] = inner(f, gk[a], p1.get_tuple(), p1.get_tuple());
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

Tensor<double> DiagnosticMatrix::project_ab(const real_function_6d& ket) const
{
    CCPairFunction<double,6> bra_cc(ao_basis, ao_basis);
    return project_xy(ket, {bra_cc});
}

Tensor<double> DiagnosticMatrix::project_ab(
        const std::vector<CCPairFunction<double,6>>& ket) const
{
    return project_xy(ket, {CCPairFunction<double,6>(ao_basis, ao_basis)});
}

// ---------------------------------------------------------------------------
//  build_Gab_bra — Schwinger BSH bra builder
// ---------------------------------------------------------------------------

std::vector<CCPairFunction<double,6>>
DiagnosticMatrix::build_Gab_bra(double energy, double lo) const
{
    MADNESS_CHECK_THROW(energy < 0.0, "build_Gab_bra: energy must be negative");
    const double thresh = FunctionDefaults<3>::get_thresh();
    const double mu     = std::sqrt(-2.0 * energy);
    const double hi     = FunctionDefaults<3>::get_cell_width().normf();

    auto fit   = GFit<double,3>::BSHFit(mu, lo, hi, thresh);
    auto c3d   = fit.coeffs();
    auto alpha = fit.exponents();
    const int nfit = c3d.dim(0);
    const int nb   = nbasis();

    std::vector<CCPairFunction<double,6>> bra;
    bra.reserve(nfit);

    for (int n = 0; n < nfit; ++n) {
        const double an  = alpha[n];
        const double w6d = c3d[n] * std::pow(an / constants::pi, 1.5);

        auto gauss = SeparatedConvolution<double,3>(
                world_, OperatorInfo(an, lo, thresh, OT_GAUSS));

        // conv[a] = g_n * ao_basis[a]  (Gaussian-convolved basis)
        auto conv = gauss(ao_basis);

        // wconv = w6d * conv  (weight absorbed into particle-1 factors)
        auto wconv = conv;
        scale(world_, wconv, w6d);

        bra.emplace_back(wconv, conv);  // decomposed CCPairFunction: g=wconv, h=conv
    }
    return bra;
}

// ---------------------------------------------------------------------------
//  project_Gab — thin wrappers
// ---------------------------------------------------------------------------

Tensor<double> DiagnosticMatrix::project_Gab(const real_function_6d& ket,
                                               double energy, double lo) const
{
    return project_xy(ket, build_Gab_bra(energy, lo));
}

Tensor<double> DiagnosticMatrix::project_Gab(
        const std::vector<CCPairFunction<double,6>>& ket,
        double energy, double lo) const
{
    return project_xy(ket, build_Gab_bra(energy, lo));
}

// ---------------------------------------------------------------------------
//  compute_errors
// ---------------------------------------------------------------------------

void DiagnosticMatrix::compute_errors()
{
    for (auto& [name, e] : entries)
        e.compute_error();
}

// ---------------------------------------------------------------------------
//  print_report
// ---------------------------------------------------------------------------

void DiagnosticMatrix::print_report(const std::string& tag) const
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

} // namespace madness
