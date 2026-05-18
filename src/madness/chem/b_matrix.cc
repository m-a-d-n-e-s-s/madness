// b_matrix.cc

#include "b_matrix.h"
#include <madness/chem/projector.h>
#include <madness/chem/nemo.h>
#include <madness/tensor/tensor_lapack.h>
#include <madness/chem/exchange_commutator.h>

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

Tensor<double> BMatrix::compute_B_KffK_6d(World& world, const Info& info) {
    return compute_BXX(world, info, [&](size_t k, size_t l) {
        auto phi_k   = CCFunction<double,3>(info.mo_ket[k], k, HOLE);
        auto phi_l   = CCFunction<double,3>(info.mo_ket[l], l, HOLE);
        auto Gscreen = make_gscreen(world, info, k, l);
        auto result  = ExchangeCommutator::apply_KffK_6d(world, phi_k, phi_l, info, &Gscreen);
        return result.KffK;
    }, -1.0);
}

// ─────────────────────────────────────────────────────────────────────────────

Tensor<double> BMatrix::compute_6d(World& world, const Info& info) {
    if (world.rank() == 0) print_header2("Computing F12 B matrix (exact 6D [K,f12])");

    if (world.rank() == 0) print_header3("B^Ue: <ij| f12 Q12 Ue |kl>");
    Tensor<double> B = compute_B_Ue(world, info);

    if (world.rank() == 0) print_header3("B^g: -<ij| f12 Q12 g12 |kl>");
    B += compute_B_g(world, info);

    if (world.rank() == 0) print_header3("B^KffK: exact 6D [K,f12] via apply_KffK_6d");
    B += compute_B_KffK_6d(world, info);

    if (world.rank() == 0) print_header3("B^X: (E_kl - E_ij) * X_{ij,kl}");
    B += compute_B_X(world, info);

    if (world.rank() == 0) {
        print("B matrix diagonal elements (6D exact):");
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

// ─────────────────────────────────────────────────────────────────────────────

Tensor<double> BMatrix::compute_via_lrf(World& world, const Info& info,
                                          const LowRankFunction<double,6>& lrf_f12,
                                          std::shared_ptr<Fock<double,3>> fock_ptr) {
    const size_t freeze = info.parameters.freeze();
    const size_t nocc   = info.mo_ket.size();
    const size_t nact   = nocc - freeze;

    if (world.rank() == 0) {
        print_header2("Computing B matrix via LRF factorization");
        print("LRF rank:", lrf_f12.rank()[0], " nact:", nact, " freeze:", freeze);
    }

    // ── Effective canonical LRF factors ─────────────────────────────────────
    // If the LRF has a non-trivial metric M, absorb it: g_eff_p = Σ_q g_q M_{qp}
    // so that f12 = Σ_p g_eff_p(r1) h_p(r2)  (canonical form).
    auto g_eff = lrf_f12.get_g();   // a_p(r1), size R
    auto h_eff = lrf_f12.get_h();   // b_p(r2), size R
    if (!lrf_f12.is_canonical()) {
        g_eff = transform(world, lrf_f12.get_g(), lrf_f12.metric);
    }
    const int R     = (int)g_eff.size();
    const int Rnact = R * (int)nact;

    // ── Phase 1: Projected orbital products ─────────────────────────────────
    // g_pi[p*nact + (i-freeze)] = a_p * phi_i  (before projection)
    vector_real_function_3d g_pi, h_pj;
    g_pi.reserve(Rnact);  h_pj.reserve(Rnact);
    for (int p = 0; p < R; ++p) {
        for (size_t i = freeze; i < nocc; ++i) g_pi.push_back(g_eff[p] * info.mo_ket[i]);
        for (size_t j = freeze; j < nocc; ++j) h_pj.push_back(h_eff[p] * info.mo_ket[j]);
    }
    truncate(world, g_pi);
    truncate(world, h_pj);

    // Batch-project out all occupied orbitals:
    //   atilde_pi = g_pi - Σ_m phi_m <R2 phi_m | g_pi>
    // S_g(m, pi) = <mo_bra[m] | g_pi[pi]>  →  shape (nocc, Rnact)
    auto S_g = matrix_inner(world, info.mo_bra, g_pi);  // (nocc, Rnact)
    auto S_h = matrix_inner(world, info.mo_bra, h_pj);

    // proj_pi = Σ_m S_g(m, pi) * phi_m
    // transform(world, v, c)[i] = Σ_j c(j,i)*v[j], so with c=S_g (nocc×Rnact) and v=mo_ket (nocc):
    //   result[pi] = Σ_m S_g(m,pi) * mo_ket[m]  — no transpose needed
    auto atilde = sub(world, g_pi, transform(world, info.mo_ket, S_g));
    auto btilde = sub(world, h_pj, transform(world, info.mo_ket, S_h));
    truncate(world, atilde);
    truncate(world, btilde);

    // R2-weighted bras  (info.R_square is the one-electron R^2 factor)
    auto R2_atilde = mul(world, info.R_square, atilde);
    auto R2_btilde = mul(world, info.R_square, btilde);

    // ── Phase 2: Projected overlap matrices Ŝ, T̂ ────────────────────────────
    // Sflat(pi, qk) = <R2 atilde_pi | atilde_qk>
    auto Sflat = matrix_inner(world, R2_atilde, atilde);   // (Rnact, Rnact)
    auto Tflat = matrix_inner(world, R2_btilde, btilde);

    // ── Phase 3: Projected Fock matrix elements F̂, Ĝ ───────────────────────
    // Apply F (without T; T is applied separately for numerical stability)
    Fock<double,3> F_no_T = *fock_ptr;
    F_no_T.remove_operator("T");
    Kinetic<double,3> Tkin(world);

    auto FV_a = F_no_T(atilde);    // (V+J-K) * atilde, size Rnact
    auto FV_b = F_no_T(btilde);

    // Mflat(pi, qk) = <R2 atilde_pi | T+V+J-K | atilde_qk>
    auto Mflat = Tkin(R2_atilde, atilde) + matrix_inner(world, R2_atilde, FV_a);
    auto Nflat = Tkin(R2_btilde, btilde) + matrix_inner(world, R2_btilde, FV_b);

    // Projection correction: Chat(pi, qk) = Σ_m eps_m <R2 atilde_pi|phi_m> <R2 phi_m|atilde_qk>
    // A(pi, m) = <R2_atilde_pi | phi_m>  →  (Rnact, nocc)
    auto A = matrix_inner(world, R2_atilde, info.mo_ket);
    auto B_ov = matrix_inner(world, R2_btilde, info.mo_ket);

    // Aeps(pi, m) = eps_m * A(pi, m)
    auto Aeps = copy(A);
    auto Beps = copy(B_ov);
    for (size_t m = 0; m < nocc; ++m) {
        double em = info.orbital_energies[m];
        for (int pi = 0; pi < Rnact; ++pi) {
            Aeps(pi, (long)m) *= em;
            Beps(pi, (long)m) *= em;
        }
    }
    // Chat(pi, qk) = Σ_m Aeps(pi,m) * A(qk,m) = inner(Aeps, A, 1, 1)
    auto Chat = inner(Aeps, A,   1, 1);   // (Rnact, Rnact)
    auto Dhat = inner(Beps, B_ov, 1, 1);

    auto Fhat = Mflat - Chat;   // projected Fock elements, electron 1
    auto Ghat = Nflat - Dhat;   // projected Fock elements, electron 2

    if (world.rank() == 0) {
        printf("  ||Mflat|| = %.6e  ||Chat|| = %.6e  ||Fhat|| = %.6e\n",
               Mflat.normf(), Chat.normf(), Fhat.normf());
        printf("  ||Sflat|| = %.6e  ||Tflat|| = %.6e\n",
               Sflat.normf(), Tflat.normf());
        double S_diag = 0.0, T_diag = 0.0;
        for (int p=0; p<Rnact; ++p) { S_diag += Sflat(p,p); T_diag += Tflat(p,p); }
        printf("  Sflat.diag_sum = %.6e  Tflat.diag_sum = %.6e\n", S_diag, T_diag);
        double X_check = 0.0;
        for (int p=0; p<Rnact; ++p) for (int q=0; q<Rnact; ++q) X_check += Sflat(p,q)*Tflat(p,q);
        printf("  X_lrf(sum S*T) = %.6e\n", X_check);
    }

    // ── Phase 4: Assemble B ──────────────────────────────────────────────────
    // B_{ij,kl} = Σ_{pq} [ Fhat_{pi,qk} Tflat_{pj,ql}
    //                     + Sflat_{pi,qk} Ghat_{pj,ql}
    //                     - E_sym * Sflat_{pi,qk} Tflat_{pj,ql} ]
    Tensor<double> B(nact, nact, nact, nact);
    B = 0.0;

    for (int ii = 0; ii < (int)nact; ++ii)
    for (int jj = 0; jj < (int)nact; ++jj)
    for (int kk = 0; kk < (int)nact; ++kk)
    for (int ll = 0; ll < (int)nact; ++ll) {
        double E_sym = 0.5*(info.orbital_energies[ii+freeze]+info.orbital_energies[jj+freeze]
                           +info.orbital_energies[kk+freeze]+info.orbital_energies[ll+freeze]);
        double FT_sum=0.0, SG_sum=0.0, X_lrf=0.0;
        for (int p = 0; p < R; ++p)
        for (int q = 0; q < R; ++q) {
            long pi = p*(long)nact + ii;
            long pj = p*(long)nact + jj;
            long qk = q*(long)nact + kk;
            long ql = q*(long)nact + ll;
            FT_sum += Fhat(pi,qk)*Tflat(pj,ql);
            SG_sum += Sflat(pi,qk)*Ghat(pj,ql);
            X_lrf  += Sflat(pi,qk)*Tflat(pj,ql);
        }
        double Bval = FT_sum + SG_sum - E_sym*X_lrf;
        if (world.rank() == 0 && ii==kk && jj==ll)
            printf("  diag(%d%d): FT=%+.6f  SG=%+.6f  X_lrf=%+.6f  -E*X=%+.6f  B=%+.6f\n",
                   (int)(ii+freeze),(int)(jj+freeze),
                   FT_sum, SG_sum, X_lrf, -E_sym*X_lrf, Bval);
        B(ii, jj, kk, ll) = Bval;
    }

    if (world.rank() == 0) {
        print("B_lrf matrix diagonal elements:");
        for (int i = 0; i < (int)nact; ++i)
            for (int j = i; j < (int)nact; ++j)
                printf("  B_lrf(%d,%d,%d,%d) = %12.8f\n",
                       (int)(i+freeze), (int)(j+freeze),
                       (int)(i+freeze), (int)(j+freeze),
                       B(i, j, i, j));
    }
    return B;
}

// ─────────────────────────────────────────────────────────────────────────────

double BMatrix::compute_via_lrf_pair(World& world, const Info& info,
                                      size_t i, size_t j,
                                      const LowRankFunction<double,6>& lrf_f12_ij,
                                      std::shared_ptr<Fock<double,3>> fock_ptr) {
    const size_t nocc = info.mo_ket.size();

    if (world.rank() == 0)
        printf("  compute_via_lrf_pair(%zu,%zu): LRF rank = %d\n",
               i, j, (int)lrf_f12_ij.rank()[0]);

    // Symmetric whitening of the LRF factors:
    //   g_white = G^{-1/2} g_can  (L²-orthonormal)
    //   h_white = G^{1/2}  h_raw  (complementary scale)
    // This balances ||alpha~|| and ||beta~|| after Q12 projection and avoids
    // catastrophic cancellation when g_can is nearly parallel to phi_i.
    auto g_raw_lrf = lrf_f12_ij.get_g();
    auto h_raw_lrf = lrf_f12_ij.get_h();
    // Absorb metric if non-canonical (makes g_eff = g_can)
    auto g_eff = lrf_f12_ij.is_canonical()
                   ? g_raw_lrf
                   : transform(world, g_raw_lrf, lrf_f12_ij.metric);

    // Compute L² Gram matrix of g_can
    auto G_mat = matrix_inner(world, g_eff, g_eff);   // L² inner products

    const int R = (int)g_eff.size();
    Tensor<double> eval, evec;
    syev(G_mat, evec, eval);   // G_mat = evec * diag(eval) * evec^T

    // Build G^{-1/2} and G^{+1/2} column-by-column
    Tensor<double> V_inv(R, R), V_fwd(R, R);
    for (int k = 0; k < R; ++k) {
        double sq = std::sqrt(std::max(1e-14, eval(k)));
        for (int p = 0; p < R; ++p) {
            V_inv(p, k) = evec(p, k) / sq;    // V * D^{-1/2}
            V_fwd(p, k) = evec(p, k) * sq;    // V * D^{+1/2}
        }
    }
    // G^{-1/2}(p,q) = sum_k V_inv(p,k)*evec(q,k), etc.
    auto G_inv_sqrt = inner(V_inv, evec, 1, 1);   // (R,R)
    auto G_fwd_sqrt = inner(V_fwd, evec, 1, 1);

    // transform(world, v, c)[p] = sum_q c(q,p)*v[q]  →  (G^{-1/2} g_can)_p  ✓
    auto g_white = transform(world, g_eff,      G_inv_sqrt);   // L²-orthonormal
    auto h_white = transform(world, h_raw_lrf,  G_fwd_sqrt);   // complementary
    auto& g_use = g_white;
    auto& h_use = h_white;

    // Phase 1: project occupied space
    // S_g(m, p) = <mo_bra[m] | g_white[p]>
    auto S_g = matrix_inner(world, info.mo_bra, g_use);  // (nocc, R)
    auto S_h = matrix_inner(world, info.mo_bra, h_use);

    // atilde[p] = g_white[p] - sum_m S_g(m,p) * phi_m
    auto atilde = sub(world, g_use, transform(world, info.mo_ket, S_g));
    auto btilde = sub(world, h_use, transform(world, info.mo_ket, S_h));
    truncate(world, atilde);
    truncate(world, btilde);

    auto R2_atilde = mul(world, info.R_square, atilde);
    auto R2_btilde = mul(world, info.R_square, btilde);

    // Phase 2: overlap matrices
    auto Sflat = matrix_inner(world, R2_atilde, atilde);  // (R, R)
    auto Tflat = matrix_inner(world, R2_btilde, btilde);

    // Phase 3: Fock matrix elements
    Fock<double,3> F_no_T = *fock_ptr;
    F_no_T.remove_operator("T");
    Kinetic<double,3> Tkin(world);

    auto FV_a = F_no_T(atilde);
    auto FV_b = F_no_T(btilde);
    auto Mflat = Tkin(R2_atilde, atilde) + matrix_inner(world, R2_atilde, FV_a);
    auto Nflat = Tkin(R2_btilde, btilde) + matrix_inner(world, R2_btilde, FV_b);

    // Phase 4: assemble diagonal element B_{ij,ij}
    double E_ij  = info.orbital_energies[i] + info.orbital_energies[j];
    double X_lrf = 0.0, FT_sum = 0.0, SG_sum = 0.0;
    for (int p = 0; p < R; ++p)
    for (int q = 0; q < R; ++q) {
        X_lrf  += Sflat(p, q) * Tflat(p, q);
        FT_sum += Mflat(p, q) * Tflat(p, q);
        SG_sum += Sflat(p, q) * Nflat(p, q);
    }
    double Bval = FT_sum + SG_sum - E_ij * X_lrf;

    if (world.rank() == 0) {
        printf("    ||Sflat||=%.4e  ||Tflat||=%.4e  X_lrf=%.6e\n",
               Sflat.normf(), Tflat.normf(), X_lrf);
        printf("    FT=%+.8f  SG=%+.8f  -E*X=%+.8f  B_pair=%+.8f\n",
               FT_sum, SG_sum, -E_ij*X_lrf, Bval);
    }
    return Bval;
}

} // namespace madness
