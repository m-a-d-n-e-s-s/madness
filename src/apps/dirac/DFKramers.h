/// \file DFKramers.h
/// \brief Quaternionic / time-reversal-symmetric linear algebra helpers for Kramers-restricted Dirac-Fock.

#ifndef MADNESS_APPS_DIRAC_DFKRAMERS_H_INCLUDED
#define MADNESS_APPS_DIRAC_DFKRAMERS_H_INCLUDED

#include <madness/tensor/solvers.h>
#include <madness/tensor/tensor.h>
#include <madness/tensor/tensor_lapack.h>
#include <madness/world/print.h>
#include <algorithm>
#include <cmath>
#include <complex>
#include <utility>

namespace madness {

/// Assemble an (n + np) x (n + np) Hermitian matrix in unpermuted basis order
/// [psi_0 .. psi_{n-1}, K*psi_0 .. K*psi_{np-1}] from M11 = <psi | M | psi> (n x n)
/// and M21 = <K*psi | M | psi> (np x n).
///
/// When \c enforce_time_reversal_even is true (always valid for the overlap metric S,
/// and valid for the Fock operator F in closed-shell Kramers-restricted SCF where
/// [F, K] = 0), exact skew-symmetry M21(i,j) = -M21(j,i) is enforced on the paired
/// np x np subblock so that the assembled matrix commutes with K to machine precision.
inline Tensor<std::complex<double>>
assemble_kramers_hermitian_matrix(const Tensor<std::complex<double>>& M11_raw,
                                  const Tensor<std::complex<double>>& M21_raw,
                                  bool enforce_time_reversal_even) {
    const long n = M11_raw.dim(0);
    const long np = M21_raw.dim(0);
    const long m = n + np;
    MADNESS_CHECK(M11_raw.dim(1) == n);
    MADNESS_CHECK(M21_raw.dim(1) == n);
    MADNESS_CHECK(np <= n);

    Tensor<std::complex<double>> M11 = 0.5 * (M11_raw + conj_transpose(M11_raw));
    Tensor<std::complex<double>> M21 = copy(M21_raw);
    if (enforce_time_reversal_even && np > 0) {
        for (long i = 0; i < np; ++i) {
            M21(i, i) = std::complex<double>(0.0, 0.0);
            for (long j = i + 1; j < np; ++j) {
                const std::complex<double> skew = 0.5 * (M21(i, j) - M21(j, i));
                M21(i, j) = skew;
                M21(j, i) = -skew;
            }
        }
    }

    Tensor<std::complex<double>> mat(m, m);
    mat(Slice(0, n - 1), Slice(0, n - 1)) = M11;
    if (np > 0) {
        mat(Slice(n, m - 1), Slice(n, m - 1)) = conj(M11(Slice(0, np - 1), Slice(0, np - 1)));
        mat(Slice(n, m - 1), Slice(0, n - 1)) = M21;
        mat(Slice(0, n - 1), Slice(n, m - 1)) = conj_transpose(M21);
    }
    mat = 0.5 * (mat + conj_transpose(mat));
    return mat;
}

/// Apply the antiunitary Kramers operator K_vec in interleaved basis order
/// [psi_0, K*psi_0, psi_1, K*psi_1, ..., psi_{np-1}, K*psi_{np-1}]:
///   (K_vec v)_{2j}   = -conj(v_{2j+1})
///   (K_vec v)_{2j+1} = +conj(v_{2j})
inline Tensor<std::complex<double>>
kramers_interleaved_partner(const Tensor<std::complex<double>>& v) {
    const long m = v.dim(0);
    MADNESS_CHECK(m % 2 == 0);
    Tensor<std::complex<double>> w(m);
    for (long j = 0; j < m / 2; ++j) {
        w(2 * j)     = -std::conj(v(2 * j + 1));
        w(2 * j + 1) =  std::conj(v(2 * j));
    }
    return w;
}

/// Metric inner product <u | M | v> = u^H * M * v.
inline std::complex<double>
metric_inner_vec(const Tensor<std::complex<double>>& u,
                 const Tensor<std::complex<double>>& M,
                 const Tensor<std::complex<double>>& v) {
    const long m = u.dim(0);
    Tensor<std::complex<double>> Mv = inner(M, v);
    std::complex<double> acc(0.0, 0.0);
    for (long i = 0; i < m; ++i) {
        acc += std::conj(u(i)) * Mv(i);
    }
    return acc;
}

/// Project a 2d x 2d interleaved transformation matrix R onto the quaternionic / symplectic
/// subspace satisfying J * conj(R) = R * J, i.e., each 2x2 block has the form
///   [ a   -conj(b) ]
///   [ b    conj(a) ].
inline Tensor<std::complex<double>>
symplectify_interleaved_matrix(const Tensor<std::complex<double>>& R) {
    const long dim = R.dim(0);
    MADNESS_CHECK(dim == R.dim(1) && dim % 2 == 0);
    Tensor<std::complex<double>> S(dim, dim);
    const long npairs = dim / 2;
    for (long r = 0; r < npairs; ++r) {
        for (long c = 0; c < npairs; ++c) {
            const std::complex<double> a = 0.5 * (R(2 * r, 2 * c) + std::conj(R(2 * r + 1, 2 * c + 1)));
            const std::complex<double> b = 0.5 * (R(2 * r + 1, 2 * c) - std::conj(R(2 * r, 2 * c + 1)));
            S(2 * r,     2 * c)     =  a;
            S(2 * r,     2 * c + 1) = -std::conj(b);
            S(2 * r + 1, 2 * c)     =  b;
            S(2 * r + 1, 2 * c + 1) =  std::conj(a);
        }
    }
    return S;
}

/// Convert generalized eigenvectors U (from LAPACK zhegv on an interleaved 2n x 2n
/// closed-shell Kramers-symmetric pencil (fock, overlap)) into a strictly Kramers-paired,
/// S-orthonormal basis U(:, 2k+1) = K_vec(U(:, 2k)) with exact partner degeneracy
/// evals[2k] == evals[2k+1], and remove arbitrary Sp(d) gauge rotations within degenerate
/// Kramers clusters for smooth KAIN subspace continuity across iterations.
inline void
kramers_pair_and_align_eigenvectors(const Tensor<std::complex<double>>& fock,
                                    const Tensor<std::complex<double>>& overlap,
                                    Tensor<std::complex<double>>& U,
                                    Tensor<double>& evals,
                                    double thresh_degenerate,
                                    bool print_info) {
    const long m = fock.dim(0);
    MADNESS_CHECK(m > 0 && m % 2 == 0);
    const long n = m / 2;

    // Step 1: Symplectic S-Gram-Schmidt extraction of n exact Kramers doublets
    // ordered by ascending eigenvalue. Standard LAPACK zhegv returns an arbitrary
    // S-unitary basis inside each 2-fold (e.g. s_{1/2}, p_{1/2}) or 4+-fold (e.g.
    // p_{3/2}, d_{5/2}) degenerate eigenspace; taking even columns directly without
    // symplectic pairing mixes Kramers partners across columns and violates
    // <psi_i | K*psi_j> = 0.
    Tensor<std::complex<double>> cand = copy(U);
    Tensor<std::complex<double>> Uk(m, m);
    Tensor<double> evals_k(m);

    for (long k = 0; k < n; ++k) {
        Tensor<double> nsq(m);
        long c_first = -1;
        for (long c = 0; c < m; ++c) {
            Tensor<std::complex<double>> col = copy(cand(_, c));
            nsq(c) = std::real(metric_inner_vec(col, overlap, col));
            if (c_first < 0 && nsq(c) > 1.0e-3) {
                c_first = c;
            }
        }
        MADNESS_CHECK(c_first >= 0);

        long c_best = c_first;
        const double eval_tol = 1.0e-5 * (1.0 + std::fabs(evals(c_first)));
        for (long c = c_first + 1; c < m; ++c) {
            if (evals(c) - evals(c_first) <= eval_tol) {
                if (nsq(c) > nsq(c_best)) c_best = c;
            } else {
                break;
            }
        }

        Tensor<std::complex<double>> w = copy(cand(_, c_best));
        for (int pass = 0; pass < 2; ++pass) {
            for (long r = 0; r < 2 * k; ++r) {
                Tensor<std::complex<double>> ur = copy(Uk(_, r));
                const std::complex<double> proj = metric_inner_vec(ur, overlap, w);
                w -= ur * proj;
            }
        }
        const double wnrm = std::sqrt(std::max(0.0, std::real(metric_inner_vec(w, overlap, w))));
        MADNESS_CHECK(wnrm > 1.0e-12);
        w.scale(1.0 / wnrm);

        Tensor<std::complex<double>> wbar = kramers_interleaved_partner(w);
        // Clean numerical roundoff in wbar against prior pairs and w
        for (long r = 0; r < 2 * k; ++r) {
            Tensor<std::complex<double>> ur = copy(Uk(_, r));
            wbar -= ur * metric_inner_vec(ur, overlap, wbar);
        }
        wbar -= w * metric_inner_vec(w, overlap, wbar);
        const double wbarnrm = std::sqrt(std::max(0.0, std::real(metric_inner_vec(wbar, overlap, wbar))));
        MADNESS_CHECK(wbarnrm > 1.0e-12);
        wbar.scale(1.0 / wbarnrm);
        // Re-sync exact algebraic partner relation after roundoff polish
        wbar = kramers_interleaved_partner(w);

        Uk(_, 2 * k)     = w;
        Uk(_, 2 * k + 1) = wbar;

        const double ek = 0.5 * (std::real(metric_inner_vec(w, fock, w)) +
                                 std::real(metric_inner_vec(wbar, fock, wbar)));
        evals_k(2 * k)     = ek;
        evals_k(2 * k + 1) = ek;

        // Deflate span{w, wbar} out of all candidate eigenvectors
        for (long c = 0; c < m; ++c) {
            Tensor<std::complex<double>> col = copy(cand(_, c));
            for (int pass = 0; pass < 2; ++pass) {
                col -= w * metric_inner_vec(w, overlap, col);
                col -= wbar * metric_inner_vec(wbar, overlap, col);
            }
            cand(_, c) = col;
        }
    }

    U = Uk;
    evals = evals_k;

    // Step 2: Quaternionic / Sp(d) polar decomposition on every degenerate cluster of
    // Kramers pairs (including single 2x2 doublets, where it removes SU(2) Kramers
    // gauge mixing and complex phase drift, and multi-pair clusters such as p_{3/2},
    // where polar SVD automatically undoes both continuous Sp(d) rotations and discrete
    // pair permutations without oscillating indices).
    if (print_info) print("          Degeneracy threshold: ", thresh_degenerate);

    long klo = 0;
    while (klo < n) {
        long khi = klo;
        while (khi + 1 < n &&
               std::fabs(evals(2 * klo) - evals(2 * (khi + 1))) <
                   thresh_degenerate * std::max(1.0, std::fabs(evals(2 * klo)))) {
            ++khi;
        }
        const long ilo = 2 * klo;
        const long ihi = 2 * khi + 1;
        const long dim = ihi - ilo + 1;

        if (print_info) {
            print("          found cluster from ", ilo + 1, " to ", ihi + 1);
            for (long kk = ilo; kk <= ihi; ++kk) {
                print("               ", evals(kk));
            }
        }

        Tensor<std::complex<double>> q = symplectify_interleaved_matrix(
            copy(U(Slice(ilo, ihi), Slice(ilo, ihi))));
        Tensor<std::complex<double>> VH(dim, dim);
        Tensor<std::complex<double>> W(dim, dim);
        Tensor<double> sigma(dim);
        svd(q, W, sigma, VH);

        if (sigma(dim - 1) > 1.0e-8) {
            Tensor<std::complex<double>> rot = symplectify_interleaved_matrix(
                conj_transpose(inner(W, VH)));
            U(_, Slice(ilo, ihi)) = inner(U(_, Slice(ilo, ihi)), rot);
        }

        // Enforce exact S-orthonormality, exact U(:, 2k+1) = K_vec(U(:, 2k)), and
        // exact Rayleigh-quotient energies for the aligned orbitals.
        for (long k = klo; k <= khi; ++k) {
            Tensor<std::complex<double>> w = copy(U(_, 2 * k));
            for (int pass = 0; pass < 2; ++pass) {
                for (long r = 2 * klo; r < 2 * k; ++r) {
                    Tensor<std::complex<double>> ur = copy(U(_, r));
                    w -= ur * metric_inner_vec(ur, overlap, w);
                }
            }
            const double wnrm = std::sqrt(std::max(0.0, std::real(metric_inner_vec(w, overlap, w))));
            if (wnrm > 1.0e-12) {
                w.scale(1.0 / wnrm);
                Tensor<std::complex<double>> wbar = kramers_interleaved_partner(w);
                U(_, 2 * k)     = w;
                U(_, 2 * k + 1) = wbar;
                const double ek = 0.5 * (std::real(metric_inner_vec(w, fock, w)) +
                                         std::real(metric_inner_vec(wbar, fock, wbar)));
                evals(2 * k)     = ek;
                evals(2 * k + 1) = ek;
            }
        }

        klo = khi + 1;
    }
}

/// Compute the inverse Hermitian square root S^{-1/2} = V * diag(evals^{-1/2}) * V^H
/// for Löwdin symmetric orthonormalization.
inline Tensor<std::complex<double>>
hermitian_inv_sqrt(const Tensor<std::complex<double>>& S_in) {
    const long m = S_in.dim(0);
    MADNESS_CHECK(m == S_in.dim(1));
    Tensor<std::complex<double>> S = 0.5 * (S_in + conj_transpose(S_in));
    Tensor<std::complex<double>> V;
    Tensor<double> evals;
    syev(S, V, evals);

    Tensor<std::complex<double>> V_scaled = copy(V);
    for (long k = 0; k < m; ++k) {
        MADNESS_CHECK(evals(k) > 1.0e-12);
        const double inv_sqrt = 1.0 / std::sqrt(evals(k));
        V_scaled(_, k).scale(inv_sqrt);
    }
    Tensor<std::complex<double>> Q = inner(V_scaled, conj_transpose(V));
    Q = 0.5 * (Q + conj_transpose(Q));
    return Q;
}

/// Enforce exact unpermuted time-reversal block structure on a 2n x 2n Hermitian matrix Q:
/// Q = [ Q11, -conj(Q21); Q21, conj(Q11) ] with Q11^H = Q11 and Q21^T = -Q21.
inline Tensor<std::complex<double>>
symplectify_unpermuted_hermitian(const Tensor<std::complex<double>>& Q_in) {
    const long m = Q_in.dim(0);
    MADNESS_CHECK(m == Q_in.dim(1) && m % 2 == 0);
    const long n = m / 2;
    Tensor<std::complex<double>> Q11 = 0.5 * (
        Q_in(Slice(0, n - 1), Slice(0, n - 1)) +
        conj(Q_in(Slice(n, m - 1), Slice(n, m - 1))));
    Q11 = 0.5 * (Q11 + conj_transpose(Q11));

    Tensor<std::complex<double>> Q21 = copy(Q_in(Slice(n, m - 1), Slice(0, n - 1)));
    for (long i = 0; i < n; ++i) {
        Q21(i, i) = std::complex<double>(0.0, 0.0);
        for (long j = i + 1; j < n; ++j) {
            const std::complex<double> skew = 0.5 * (Q21(i, j) - Q21(j, i));
            Q21(i, j) = skew;
            Q21(j, i) = -skew;
        }
    }
    return assemble_kramers_hermitian_matrix(Q11, Q21, true);
}

/// Solve the real-Hilbert-space KAIN subspace equations (`Q(i, j) = Re(<u_i | r_j>)`)
/// with SCF.cc-matched singular-value conditioning (`c.absmax() < cabsmax`, escalating
/// `rcond` from `1e-12` up to `rcond_max`, and resetting stale history to the latest
/// iterate `c = [1.0]` when the subspace becomes linearly dependent at the MRA noise floor).
///
/// Using real weights `c_i in R` is essential for Dirac-Fock because:
///   1. The Hartree-Fock / Dirac-Fock BSH map `F(psi, conj(psi))` is only R-differentiable
///      (not C-holomorphic), so `Dr[sum c_i du_i] = sum c_i Dr[du_i]` holds iff `c_i in R`.
///   2. Time-reversal `K` is antiunitary (`K(c_i psi) = conj(c_i) K(psi)`), so real weights
///      preserve Kramers partner covariance without injecting spurious spin-flip leakage.
inline std::pair<Tensor<double>, bool>
solve_real_kain_weights(const Tensor<double>& Q,
                        double cabsmax = 3.0,
                        double rcond_max = 0.01) {
    const long m = Q.dim(0);
    MADNESS_CHECK(m >= 1 && Q.dim(1) == m);
    if (m == 1) {
        Tensor<double> c(1L);
        c(0L) = 1.0;
        return {c, false};
    }

    double rcond = 1.0e-12;
    while (true) {
        Tensor<double> c = KAIN(Q, rcond);
        if (c.absmax() < cabsmax) {
            return {c, false};
        }
        if (rcond < rcond_max) {
            rcond *= 100.0;
        } else {
            Tensor<double> c_reset(1L);
            c_reset(0L) = 1.0;
            return {c_reset, true};
        }
    }
}

} // namespace madness

#endif // MADNESS_APPS_DIRAC_DFKRAMERS_H_INCLUDED
