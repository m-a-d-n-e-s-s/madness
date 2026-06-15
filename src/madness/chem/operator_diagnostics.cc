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
//  ref_Gab / ref_GQab — 3D-only reference builders
// ---------------------------------------------------------------------------

template<typename T, std::size_t NDIM>
Tensor<T> DiagnosticMatrix<T,NDIM>::ref_Gab(const ElementProvider& elements,
                                            double energy, double lo) const
{
    Tensor<T> ref(nbasis(), nbasis());
    for (const auto& bk : build_Gab_bra(energy, lo))
        ref += elements(bk.get_a(), bk.get_b());
    return ref;
}

template<typename T, std::size_t NDIM>
Tensor<T> DiagnosticMatrix<T,NDIM>::ref_GQab(const ElementProvider& elements,
        const std::vector<Function<T,LDIM>>& occ_ket,
        const std::vector<Function<T,LDIM>>& occ_bra,
        double energy, double lo) const
{
    MADNESS_CHECK_THROW(!occ_ket.empty() && occ_ket.size() == occ_bra.size(),
                        "ref_GQab: invalid occupied spaces");
    const long nb   = nbasis();
    const long nocc = static_cast<long>(occ_ket.size());
    const Slice s_ab(0, nb - 1), s_occ(nb, nb + nocc - 1);

    Tensor<T> ref(nb, nb);
    for (const auto& bk : build_Gab_bra(energy, lo)) {
        std::vector<Function<T,LDIM>> p1 = bk.get_a();   // weighted ã_a
        p1.insert(p1.end(), occ_bra.begin(), occ_bra.end());
        std::vector<Function<T,LDIM>> p2 = bk.get_b();   // b̃_b
        p2.insert(p2.end(), occ_bra.begin(), occ_bra.end());

        const Tensor<T> S1  = matrix_inner(world_, bk.get_a(), occ_ket);             // (nb,nocc)
        const Tensor<T> S2t = transpose(matrix_inner(world_, bk.get_b(), occ_ket));  // (nocc,nb)

        Tensor<T> M = elements(p1, p2);
        const Tensor<T> A = copy(M(s_ab,  s_ab));
        const Tensor<T> B = copy(M(s_occ, s_ab));
        const Tensor<T> C = copy(M(s_ab,  s_occ));
        const Tensor<T> D = copy(M(s_occ, s_occ));
        ref += A - inner(S1, B) - inner(C, S2t) + inner(S1, inner(D, S2t));
    }
    return ref;
}

// ---------------------------------------------------------------------------
//  print_pair_energy_report
// ---------------------------------------------------------------------------

template<typename T, std::size_t NDIM>
void print_pair_energy_report(const DiagnosticMatrix<T,NDIM>& dm,
                              const std::vector<std::string>& pieces,
                              const std::vector<double>& signs,
                              const std::string& title)
{
    MADNESS_CHECK_THROW(pieces.size() == signs.size(),
                        "print_pair_energy_report: pieces/signs size mismatch");
    constexpr std::size_t nbuf = 256;
    char buf[nbuf];
    auto print_line = [&buf](const std::string& name, double e6d, double r3d) {
        std::snprintf(buf, nbuf, "  %-12s 6d %15.8e   3d-ref %15.8e   diff %15.8e",
                      name.c_str(), e6d, r3d, e6d - r3d);
        print(std::string(buf));
    };

    print(title);
    double tot6d = 0.0, tot3d = 0.0;
    for (std::size_t p = 0; p < pieces.size(); ++p) {
        const auto& e = dm.entries.at(pieces[p]);
        print_line(pieces[p], e.result(0,1), e.ref(0,1));
        tot6d += signs[p] * e.result(0,1);
        tot3d += signs[p] * e.ref(0,1);
    }
    print_line("total", tot6d, tot3d);
}

// ---------------------------------------------------------------------------
//  pair_energy_weights — strong-orthogonality-projected g12 bra weights
// ---------------------------------------------------------------------------

template<typename T, std::size_t NDIM>
Tensor<T> pair_energy_weights(World& world,
                              const std::vector<Function<T,NDIM/2>>& observer,
                              const Function<T,NDIM/2>& phi_i_bra,
                              const Function<T,NDIM/2>& phi_j_bra,
                              const std::vector<Function<T,NDIM/2>>& occ_ket,
                              const std::vector<Function<T,NDIM/2>>& occ_bra,
                              const SeparatedConvolution<T,NDIM/2>& g12)
{
    constexpr std::size_t LDIM = NDIM/2;
    const int nb = static_cast<int>(observer.size());

    // Q-project the observer functions: Q a = a - sum_m occ_ket_m <occ_bra_m|a>
    std::vector<Function<T,LDIM>> Qobs = observer;
    if (!occ_ket.empty()) {
        const Tensor<T> S = matrix_inner(world, occ_bra, observer);  // (nocc,nb)
        Qobs = sub(world, observer, transform(world, occ_ket, S));
    }

    // products with the pair-bra orbitals, then the g12 convolution
    const std::vector<Function<T,LDIM>> iQ = Qobs * phi_i_bra;   // {phi_i_bra * Q a}
    const std::vector<Function<T,LDIM>> jQ = Qobs * phi_j_bra;
    std::vector<Function<T,LDIM>> g_iQ(nb), g_jQ(nb);
    for (int b=0; b<nb; ++b) { g_iQ[b]=g12(iQ[b]); g_jQ[b]=g12(jQ[b]); }

    // W(a,b) = 2<iQ_a|g12|jQ_b> - <jQ_a|g12|iQ_b>
    return 2.0*matrix_inner(world, iQ, g_jQ) - matrix_inner(world, jQ, g_iQ);
}

// ---------------------------------------------------------------------------
//  print_pair_energy_report_RI — contract result/ref tensors with W
// ---------------------------------------------------------------------------

template<typename T, std::size_t NDIM>
void print_pair_energy_report_RI(const DiagnosticMatrix<T,NDIM>& dm,
                                 const std::vector<std::string>& pieces,
                                 const std::vector<double>& signs,
                                 const Tensor<T>& W,
                                 double factor,
                                 const std::string& title)
{
    MADNESS_CHECK_THROW(pieces.size() == signs.size(),
                        "print_pair_energy_report_RI: pieces/signs size mismatch");
    auto contract = [&W,factor](const Tensor<T>& M) {
        double e = 0.0;
        for (long a=0; a<M.dim(0); ++a)
            for (long b=0; b<M.dim(1); ++b) e += W(a,b)*M(a,b);
        return factor*e;
    };
    constexpr std::size_t nbuf = 256;
    char buf[nbuf];
    print(title);
    std::snprintf(buf,nbuf,"  %-12s %18s %18s %15s","piece","Expr2 (3d-ket)","Expr3 (6d-ket)","diff");
    print(std::string(buf));
    double tot3d=0.0, tot6d=0.0;
    for (std::size_t p=0; p<pieces.size(); ++p) {
        const auto& e = dm.entries.at(pieces[p]);
        const double e3d = signs[p]*contract(e.ref);     // Expr2: 3d Schwinger ket
        const double e6d = signs[p]*contract(e.result);  // Expr3: 6d projected ket
        std::snprintf(buf,nbuf,"  %-12s % .8e % .8e % .3e",pieces[p].c_str(),e3d,e6d,e6d-e3d);
        print(std::string(buf));
        tot3d+=e3d; tot6d+=e6d;
    }
    std::snprintf(buf,nbuf,"  %-12s % .8e % .8e % .3e","total",tot3d,tot6d,tot6d-tot3d);
    print(std::string(buf));
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

template void print_pair_energy_report<double, 6>(
        const DiagnosticMatrix<double,6>&, const std::vector<std::string>&,
        const std::vector<double>&, const std::string&);

template Tensor<double> pair_energy_weights<double, 6>(
        World&, const std::vector<Function<double,3>>&,
        const Function<double,3>&, const Function<double,3>&,
        const std::vector<Function<double,3>>&, const std::vector<Function<double,3>>&,
        const SeparatedConvolution<double,3>&);

template void print_pair_energy_report_RI<double, 6>(
        const DiagnosticMatrix<double,6>&, const std::vector<std::string>&,
        const std::vector<double>&, const Tensor<double>&, double, const std::string&);

} // namespace madness
