// b_matrix.cc

#include "b_matrix.h"
#include <madness/chem/projector.h>

using namespace madness;

/// BSH screening operator for pair (k,l), used for apply_Ue / apply_KffK screening.
static real_convolution_6d make_gscreen(World& world, const Info& info, size_t k, size_t l) {
    double bsh_eps = info.orbital_energies[k] + info.orbital_energies[l];
    auto Gscreen = BSHOperator<6>(world, sqrt(-2.0 * bsh_eps),
                                  info.parameters.lo(), info.parameters.thresh_bsh_6D());
    Gscreen.modified() = true;
    return Gscreen;
}

/// Precompute Q12† f12|ij_bra> for every active bra pair (i,j), where ij_bra
/// uses mo_bra (R²·φ).  Using Q12† (mo_ket/mo_bra swapped relative to the ket-
/// side Q12) is required because Q12 acts here from the right on the bra:
///   <ij|f12 Q12|V> = inner(Q12† f12|ij_bra>, V)
/// With Q12† the subtracted overlap <φₘφₙ|f12|R²φᵢR²φⱼ> correctly matches
/// the R²-weighted projectors without double-counting R² factors.
/// Indexed as bras[(i-freeze)*nact + (j-freeze)].
using BraPairs = std::vector<std::vector<CCPairFunction<double,6>>>;

static BraPairs precompute_bras(
        size_t freeze, size_t nocc,
        const StrongOrthogonalityProjector<double,3>& Q12,
        const std::shared_ptr<CCConvolutionOperator<double,3>>& f12_op,
        const std::vector<Function<double,3>>& mo_bra) {

    size_t nact = nocc - freeze;
    BraPairs bras(nact * nact);
    for (size_t i = freeze; i < nocc; ++i)
        for (size_t j = freeze; j < nocc; ++j) {
            CCPairFunction<double,6> f12_ij(f12_op, mo_bra[i], mo_bra[j]);
            bras[(i - freeze) * nact + (j - freeze)] = madness::apply(Q12, f12_ij);
        }
    return bras;
}

/// Compute <ij|f12 Q12|V_kl> for all active (i,j) using precomputed bras.
static Tensor<double> inner_over_bra_pairs(
        size_t freeze, size_t nocc,
        const std::vector<CCPairFunction<double,6>>& V_kl,
        const BraPairs& bras) {

    size_t nact = nocc - freeze;
    Tensor<double> result(nact, nact);
    for (size_t i = freeze; i < nocc; ++i)
        for (size_t j = freeze; j < nocc; ++j)
            result(i - freeze, j - freeze) =
                inner(bras[(i - freeze) * nact + (j - freeze)], V_kl);
    return result;
}

/// Core B-matrix driver.
///
/// Precomputes Q12 f12|ij> once, then iterates over all active ket pairs (k,l),
/// calling make_V(k,l) to obtain the ket potential V_kl, and accumulates
///   B(i-freeze, j-freeze, k-freeze, l-freeze) = sign * <ij| f12 Q12 | V_kl>
///
/// @param sign   overall sign applied to every element (pass -1 for contributions
///               like -V and -[K,f] that enter with a minus sign)
static Tensor<double> compute_BXX(
        World& world, const Info& info,
        std::function<std::vector<CCPairFunction<double,6>>(size_t k, size_t l)> make_V,
        double sign = 1.0) {

    const size_t freeze = info.parameters.freeze();
    const size_t nocc   = info.mo_ket.size();
    const size_t nact   = nocc - freeze;

    auto f12_op = CCConvolutionOperatorPtr<double,3>(world, OT_F12, info.parameters);

    StrongOrthogonalityProjector<double,3> Q12_dagger(world);
    Q12_dagger.set_spaces(info.mo_ket, info.mo_bra, info.mo_ket, info.mo_bra);

    auto bras = precompute_bras(freeze, nocc, Q12_dagger, f12_op, info.mo_bra);

    Tensor<double> B(nact, nact, nact, nact);
    B = 0.0;

    for (size_t k = freeze; k < nocc; ++k) {
        for (size_t l = freeze; l < nocc; ++l) {
            auto V_kl = make_V(k, l);
            auto row  = inner_over_bra_pairs(freeze, nocc, V_kl, bras);
            for (size_t i = freeze; i < nocc; ++i)
                for (size_t j = freeze; j < nocc; ++j)
                    B(i-freeze, j-freeze, k-freeze, l-freeze) = sign * row(i-freeze, j-freeze);
        }
    }
    return B;
}

namespace madness {

// ─────────────────────────────────────────────────────────────────────────────

Tensor<double> BMatrix::compute_X_matrix(World& world, const Info& info) {
    auto f12_op = CCConvolutionOperatorPtr<double,3>(world, OT_F12, info.parameters);
    return compute_BXX(world, info, [&](size_t k, size_t l) {
        return std::vector<CCPairFunction<double,6>>{
            CCPairFunction<double,6>(f12_op, info.mo_ket[k], info.mo_ket[l])
        };
    });
}

// ─────────────────────────────────────────────────────────────────────────────

Tensor<double> BMatrix::compute_B_Ue(World& world, const Info& info) {
    return compute_BXX(world, info, [&](size_t k, size_t l) {
        auto phi_k   = CCFunction<double,3>(info.mo_ket[k], k, HOLE);
        auto phi_l   = CCFunction<double,3>(info.mo_ket[l], l, HOLE);
        auto Gscreen = make_gscreen(world, info, k, l);
        return std::vector<CCPairFunction<double,6>>{
            CCPairFunction<double,6>(CCPotentials::apply_Ue(world, phi_k, phi_l, info, &Gscreen))
        };
    });
}

// ─────────────────────────────────────────────────────────────────────────────

Tensor<double> BMatrix::compute_B_g(World& world, const Info& info) {
    auto g12_op = CCConvolutionOperatorPtr<double,3>(world, OT_G12, info.parameters);
    return compute_BXX(world, info, [&](size_t k, size_t l) {
        return std::vector<CCPairFunction<double,6>>{
            CCPairFunction<double,6>(g12_op, info.mo_ket[k], info.mo_ket[l])
        };
    }, -1.0);
}

// ─────────────────────────────────────────────────────────────────────────────

Tensor<double> BMatrix::compute_B_KffK(World& world, const Info& info) {
    return compute_BXX(world, info, [&](size_t k, size_t l) {
        auto phi_k   = CCFunction<double,3>(info.mo_ket[k], k, HOLE);
        auto phi_l   = CCFunction<double,3>(info.mo_ket[l], l, HOLE);
        auto Gscreen = make_gscreen(world, info, k, l);
        return CCPotentials::apply_KffK(world, phi_k, phi_l, info, &Gscreen);
    }, -1.0);
}

// ─────────────────────────────────────────────────────────────────────────────

Tensor<double> BMatrix::compute_B_X(World& world, const Info& info) {
    const size_t freeze = info.parameters.freeze();
    const size_t nocc   = info.mo_ket.size();
    const size_t nact   = nocc - freeze;

    Tensor<double> X   = compute_X_matrix(world, info);
    Tensor<double> B_X(nact, nact, nact, nact);
    B_X = 0.0;

    for (size_t i = freeze; i < nocc; ++i) {
        double E_i = info.orbital_energies[i];
        for (size_t j = freeze; j < nocc; ++j) {
            double E_ij = E_i + info.orbital_energies[j];
            for (size_t k = freeze; k < nocc; ++k)
                for (size_t l = freeze; l < nocc; ++l) {
                    double E_kl = info.orbital_energies[k] + info.orbital_energies[l];
                    B_X(i-freeze, j-freeze, k-freeze, l-freeze) =
                        (E_kl - E_ij) * X(i-freeze, j-freeze, k-freeze, l-freeze);
                }
        }
    }
    return B_X;
}

// ─────────────────────────────────────────────────────────────────────────────

Tensor<double> BMatrix::compute(World& world, const Info& info) {
    if (world.rank() == 0) print_header2("Computing F12 B matrix");

    if (world.rank() == 0) print_header3("B^Ue: <ij| f12 Q12 Ue |kl>");
    Tensor<double> B = compute_B_Ue(world, info);

    if (world.rank() == 0) print_header3("B^g: -<ij| f12 Q12 g12 |kl>");
    B += compute_B_g(world, info);

    if (world.rank() == 0) print_header3("B^KffK: -<ij| f12 Q12 [K12,f12] |kl>");
    B += compute_B_KffK(world, info);

    if (world.rank() == 0) print_header3("B^X: (E_kl - E_ij) * X_{ij,kl}");
    B += compute_B_X(world, info);

    if (world.rank() == 0) {
        print("B matrix diagonal elements:");
        const size_t freeze = info.parameters.freeze();
        const size_t nocc   = info.mo_ket.size();
        for (size_t i = freeze; i < nocc; ++i)
            for (size_t j = i; j < nocc; ++j)
                printf("  B(%zu,%zu,%zu,%zu) = %12.8f\n",
                       i, j, i, j, B(i-freeze, j-freeze, i-freeze, j-freeze));
    }
    return B;
}

// ─────────────────────────────────────────────────────────────────────────────

double BMatrix::compute_hylleraas_pair(World& world, const Info& info,
                                        int i, int j,
                                        const real_function_6d& u_ij) {
    auto g12_op = CCConvolutionOperatorPtr<double,3>(world, OT_G12, info.parameters);
    auto f12_op = CCConvolutionOperatorPtr<double,3>(world, OT_F12, info.parameters);
    auto Gscreen = make_gscreen(world, info, i, j);

    // c_ij = Ue|ij> − [K,f12]|ij>  (stored in two separate pieces for sign control)
    auto phi_i   = CCFunction<double,3>(info.mo_ket[i], i, HOLE);
    auto phi_j   = CCFunction<double,3>(info.mo_ket[j], j, HOLE);
    auto Ue_ij   = CCPotentials::apply_Ue  (world, phi_i, phi_j, info, &Gscreen);
    auto KffK_ij = CCPotentials::apply_KffK(world, phi_i, phi_j, info, &Gscreen);

    // Precompute Q12† f12|ij_bra>  — reused by Terms 2 and 4
    StrongOrthogonalityProjector<double,3> Q12d(world);
    Q12d.set_spaces(info.mo_ket, info.mo_bra, info.mo_ket, info.mo_bra);
    CCPairFunction<double,6> f12_bra(f12_op, info.mo_bra[i], info.mo_bra[j]);
    auto Qf12_bra = madness::apply(Q12d, f12_bra);   // vector<CCPairFunction>

    using vCCPF = std::vector<CCPairFunction<double,6>>;
    vCCPF u_vec  {CCPairFunction<double,6>(u_ij)};
    vCCPF Ue_vec {CCPairFunction<double,6>(Ue_ij)};

    // Term 1 + hc:  2 * <ij|g12|u_ij>  (real orbitals, hc = forward term)
    vCCPF g12_bra{CCPairFunction<double,6>(g12_op, info.mo_bra[i], info.mo_bra[j])};
    double term1 = 2.0 * inner(g12_bra, u_vec);

    // Term 2:  <ij|g12 Q f12|ij>  = <ij|f12 Q g12|ij>  (real orbitals)
    //          = inner(Q12† f12|ij_bra>, g12|ij_ket>)
    vCCPF g12_ket{CCPairFunction<double,6>(g12_op, info.mo_ket[i], info.mo_ket[j])};
    double term2 = inner(Qf12_bra, g12_ket);

    // Term 3:  −<u_ij|c_ij> = −<u|Ue|ij> + <u|[K,f12]|ij>
    double term3 = -inner(u_vec, Ue_vec) + inner(u_vec, KffK_ij);

    // Term 4:  −<ij|f12 Q|c_ij> = −inner(Q12† f12|ij_bra>, Ue|ij>) + inner(..., KffK|ij>)
    double term4 = -inner(Qf12_bra, Ue_vec) + inner(Qf12_bra, KffK_ij);

    if (world.rank() == 0)
        printf("  Hylleraas(%d,%d): V+hc=%+.8f  gQf=%+.8f  uc=%+.8f  fQc=%+.8f  total=%+.8f\n",
               i, j, term1, term2, term3, term4, term1+term2+term3+term4);

    return term1 + term2 + term3 + term4;
}

} // namespace madness
