/*
  This file is part of MADNESS.

  Copyright (C) 2007,2010 Oak Ridge National Laboratory

  This program is free software; you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation; either version 2 of the License, or
  (at your option) any later version.

  This program is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU General Public License for more details.

  You should have received a copy of the GNU General Public License
  along with this program; if not, write to the Free Software
  Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307 USA

  For more information please contact:

  Robert J. Harrison
  Oak Ridge National Laboratory
  One Bethel Valley Road
  P.O. Box 2008, MS-6367

  email: harrisonrj@ornl.gov
  tel:   865-241-3937
  fax:   865-572-0680
*/

/// \file lcao_integrals.cc
/// \brief integrals over Cartesian Gaussians from separated (Gaussian-sum) kernels

#include <madness/chem/lcao_integrals.h>
#include <madness/constants.h>
#include <madness/mra/gfit.h>
#include <madness/tensor/tensor_lapack.h>
#include <madness/world/MADworld.h>

#include <algorithm>
#include <atomic>
#include <cmath>
#include <limits>

namespace madness {
namespace lcao {

namespace {

/// powers per center in the 1D factors: angular momentum up to 13 (kinetic adds 2)
constexpr int maxpow = 16;

/// p[i] = x^i for i = 0..n
inline void powers(const double x, const int n, double* p) {
    p[0] = 1.0;
    for (int i = 1; i <= n; ++i) p[i] = p[i - 1] * x;
}

/// exp(-a(x-A)^2) exp(-b(x-B)^2) = K exp(-p(x-P)^2)
struct GaussianProduct1D {
    double p, P, K;
    GaussianProduct1D(const double a, const double A, const double b, const double B)
        : p(a + b), P((a * A + b * B) / (a + b)), K(std::exp(-a * b / (a + b) * (A - B) * (A - B))) {}
};

/// three-center 1D factor, for i <= imax and j <= jmax:
///   s[i*(jmax+1)+j] = int (x-A)^i (x-B)^j exp(-a(x-A)^2 - b(x-B)^2 - t(x-C)^2) dx
/// t = 0 gives the overlap.
void gaussian_1d(const int imax, const int jmax, const double a, const double A, const double b,
                 const double B, const double t, const double C, const GaussHermiteRule& gh, double* s) {
    const GaussianProduct1D ab(a, A, b, B);
    const double p = ab.p + t;
    const double P = (ab.p * ab.P + t * C) / p;
    const double pref = ab.K * std::exp(-ab.p * t / p * (ab.P - C) * (ab.P - C)) / std::sqrt(p);
    const int n = (imax + jmax) / 2 + 1;
    const std::vector<double>& y = gh.nodes(n);
    const std::vector<double>& w = gh.weights(n);
    const double scale = 1.0 / std::sqrt(p);
    const int nj = jmax + 1;
    std::fill(s, s + (imax + 1) * nj, 0.0);
    double pa[maxpow], pb[maxpow];
    for (int k = 0; k < n; ++k) {
        const double x = P + y[k] * scale;
        powers(x - A, imax, pa);
        powers(x - B, jmax, pb);
        for (int i = 0; i <= imax; ++i) {
            const double wi = w[k] * pa[i];
            for (int j = 0; j <= jmax; ++j) s[i * nj + j] += wi * pb[j];
        }
    }
    for (int i = 0; i < (imax + 1) * nj; ++i) s[i] *= pref;
}

/// two-electron 1D factor for one kernel term exp(-t(x1-x2)^2):
///   g[((i*(lb+1)+j)*(lc+1)+k)*(ld+1)+l] = int int (x1-A)^i (x1-B)^j exp(-a(x1-A)^2 - b(x1-B)^2)
///        exp(-t(x1-x2)^2) (x2-C)^k (x2-D)^l exp(-c(x2-C)^2 - d(x2-D)^2) dx1 dx2
///
/// The Gaussian part is exp(-(x-mu)^T M (x-mu)) exp(-omega (P-Q)^2) with
/// M = [[p+t, -t], [-t, q+t]]. With M = L L^T and x = mu + L^{-T} y it becomes
/// exp(-|y|^2), so tensor-product Gauss-Hermite with n = deg/2+1 nodes per axis
/// is exact. det M, omega, mu and L are written so that none of them cancels for
/// t >> p, q, the regime of the short-range kernel terms.
void coulomb_1d(const int la, const int lb, const int lc, const int ld,
                const GaussianProduct1D& ab, const double A, const double B,
                const GaussianProduct1D& cd, const double C, const double D,
                const double t, const GaussHermiteRule& gh, double* g) {
    const double p = ab.p, q = cd.p, P = ab.P, Q = cd.P;
    const double det = p * q + t * (p + q);
    const double omega = p * q * t / det;
    const double mu1 = ((q + t) * p * P + t * q * Q) / det;
    const double mu2 = (t * p * P + (p + t) * q * Q) / det;
    const double l11 = std::sqrt(p + t);
    const double l22 = std::sqrt(det / (p + t));
    const double c12 = t / ((p + t) * l22);
    const double pref = ab.K * cd.K * std::exp(-omega * (P - Q) * (P - Q)) / std::sqrt(det);

    const int n = (la + lb + lc + ld) / 2 + 1;
    const std::vector<double>& y = gh.nodes(n);
    const std::vector<double>& w = gh.weights(n);
    const int nb = lb + 1, nd = ld + 1;
    const int nbra = (la + 1) * nb, nket = (lc + 1) * nd;
    std::fill(g, g + nbra * nket, 0.0);
    double pa[maxpow], pb[maxpow], pc[maxpow], pd[maxpow], u[maxpow * maxpow], v[maxpow * maxpow];
    for (int k1 = 0; k1 < n; ++k1) {
        for (int k2 = 0; k2 < n; ++k2) {
            const double x1 = mu1 + y[k1] / l11 + c12 * y[k2];
            const double x2 = mu2 + y[k2] / l22;
            const double wk = w[k1] * w[k2];
            powers(x1 - A, la, pa);
            powers(x1 - B, lb, pb);
            powers(x2 - C, lc, pc);
            powers(x2 - D, ld, pd);
            for (int i = 0; i <= la; ++i)
                for (int j = 0; j <= lb; ++j) u[i * nb + j] = wk * pa[i] * pb[j];
            for (int k = 0; k <= lc; ++k)
                for (int l = 0; l <= ld; ++l) v[k * nd + l] = pc[k] * pd[l];
            for (int ij = 0; ij < nbra; ++ij)
                for (int kl = 0; kl < nket; ++kl) g[ij * nket + kl] += u[ij] * v[kl];
        }
    }
    for (int i = 0; i < nbra * nket; ++i) g[i] *= pref;
}

/// copy the computed block (i,j) into (j,i) as well
void symmetrize_block(Tensor<double>& m, const Shell& sa, const Shell& sb) {
    for (int i = 0; i < sa.ncart(); ++i)
        for (int j = 0; j < sb.ncart(); ++j) m(sb.offset + j, sa.offset + i) = m(sa.offset + i, sb.offset + j);
}

/// write the block of shell quartet (ab|cd) into G, with the 8-fold permutational symmetry
void scatter_quartet(PackedERI& G, const Shell& sa, const Shell& sb, const Shell& sc, const Shell& sd,
                     const double* block) {
    const int na = sa.ncart(), nb = sb.ncart(), nc = sc.ncart(), nd = sd.ncart();
    for (int i = 0, n = 0; i < na; ++i)
        for (int j = 0; j < nb; ++j)
            for (int k = 0; k < nc; ++k)
                for (int l = 0; l < nd; ++l, ++n)
                    G(sa.offset + i, sb.offset + j, sc.offset + k, sd.offset + l) = block[n];
}

/// the same for the full nbf^4 tensor of eri_reference
void scatter_quartet(Tensor<double>& G, const Shell& sa, const Shell& sb, const Shell& sc, const Shell& sd,
                     const double* block) {
    const int na = sa.ncart(), nb = sb.ncart(), nc = sc.ncart(), nd = sd.ncart();
    for (int i = 0, n = 0; i < na; ++i) {
        for (int j = 0; j < nb; ++j) {
            for (int k = 0; k < nc; ++k) {
                for (int l = 0; l < nd; ++l, ++n) {
                    const long mu = sa.offset + i, nu = sb.offset + j;
                    const long la = sc.offset + k, si = sd.offset + l;
                    const double v = block[n];
                    G(mu, nu, la, si) = v;
                    G(nu, mu, la, si) = v;
                    G(mu, nu, si, la) = v;
                    G(nu, mu, si, la) = v;
                    G(la, si, mu, nu) = v;
                    G(si, la, mu, nu) = v;
                    G(la, si, nu, mu) = v;
                    G(si, la, nu, mu) = v;
                }
            }
        }
    }
}

/// exp(-a|r-A|^2) exp(-b|r-B|^2) = K exp(-p|r-P|^2), for one pair of primitives of two shells
struct PrimitivePair {
    double p = 0.0, K = 0.0;
    std::array<double,3> P = {0.0, 0.0, 0.0};
    int ia = 0, ib = 0;         ///< the two primitives within their shells
};

/// all primitive pairs of two shells
std::vector<PrimitivePair> primitive_pairs(const Shell& sa, const Shell& sb) {
    double ab2 = 0.0;
    for (int x = 0; x < 3; ++x) ab2 += (sa.center[x] - sb.center[x]) * (sa.center[x] - sb.center[x]);
    std::vector<PrimitivePair> pairs;
    for (std::size_t ia = 0; ia < sa.expnt.size(); ++ia) {
        for (std::size_t ib = 0; ib < sb.expnt.size(); ++ib) {
            const double a = sa.expnt[ia], b = sb.expnt[ib];
            PrimitivePair pp;
            pp.p = a + b;
            pp.K = std::exp(-a * b / pp.p * ab2);
            for (int x = 0; x < 3; ++x) pp.P[x] = (a * sa.center[x] + b * sb.center[x]) / pp.p;
            pp.ia = int(ia);
            pp.ib = int(ib);
            pairs.push_back(pp);
        }
    }
    return pairs;
}

/// what all primitive quartets of a shell quartet (ab|cd) share
struct QuartetLayout {
    int la, lb, lc, ld;
    std::array<double,3> A, B, C, D;
    int ntab;                   ///< size of one 1D table
    int nscratch;               ///< size of the scratch space of moments_1d
    std::vector<std::array<int,3>> index;   ///< per component quartet, its position in the x, y and z tables

    QuartetLayout(const Shell& sa, const Shell& sb, const Shell& sc, const Shell& sd)
        : la(sa.l), lb(sb.l), lc(sc.l), ld(sd.l), A(sa.center), B(sb.center), C(sc.center), D(sd.center),
          ntab((sa.l + 1) * (sb.l + 1) * (sc.l + 1) * (sd.l + 1)) {
        const int nn = la + lb + 1, nm = lc + ld + 1;
        nscratch = nn * nm + (ld + 1) * nm + nn * (lc + 1) * (ld + 1) + (lb + 1) * nn;
        const auto ca = cartesian_components(la), cb = cartesian_components(lb);
        const auto cc = cartesian_components(lc), cd = cartesian_components(ld);
        for (const auto& i : ca)
            for (const auto& j : cb)
                for (const auto& k : cc)
                    for (const auto& l : cd) {
                        std::array<int,3> pos;
                        for (int x = 0; x < 3; ++x) pos[x] = ((i[x] * (lb + 1) + j[x]) * (lc + 1) + k[x]) * (ld + 1) + l[x];
                        index.push_back(pos);
                    }
    }

    int nblock() const { return int(index.size()); }
};

/// one axis of the two-electron 1D factor of a primitive quartet, for one kernel term:
///   g[((i*(lb+1)+j)*(lc+1)+k)*(ld+1)+l] = pi E[(x1-A)^i (x1-B)^j (x2-C)^k (x2-D)^l]
/// for (x1,x2) Gaussian with mean (A + c00, C + c00p), variances b10 and b01 and
/// covariance b00, which is what coulomb_1d's quadrature sums are.
///
/// With X1 = x1-A and X2 = x2-C the moments follow from Stein's identity, the
/// vertical recurrence of Rys quadrature:
///   E[X1^{n+1} X2^m] = c00  E[X1^n X2^m] + n b10 E[X1^{n-1} X2^m] + m b00 E[X1^n X2^{m-1}]
///   E[X1^n X2^{m+1}] = c00p E[X1^n X2^m] + m b01 E[X1^n X2^{m-1}] + n b00 E[X1^{n-1} X2^m]
/// and the powers of x1-B and x2-D from x1-B = X1 + (A-B), x2-D = X2 + (C-D), the
/// horizontal recurrence. Exact, like the quadrature it replaces.
void moments_1d(const QuartetLayout& L, const double c00, const double c00p, const double b10, const double b01,
                const double b00, const double AB, const double CD, double* g, double* scratch) {
    const int la = L.la, lb = L.lb, lc = L.lc, ld = L.ld;
    const int nn = la + lb + 1, nm = lc + ld + 1;
    const int nkl = (lc + 1) * (ld + 1);
    // M[n*nm+m] = pi E[X1^n X2^m], K[n*nkl+k*(ld+1)+l] = pi E[X1^n X2^k (x2-D)^l]. Without a
    // horizontal step K has the layout of g (lb = 0) and M that of K (ld = 0), so they share
    // storage and the common quartets need no copies.
    double* h = scratch + nn * nm;          // ket recurrence, one n at a time: h[l*nm+m]
    double* K = (lb == 0) ? g : h + (ld + 1) * nm;
    double* M = (ld == 0) ? K : scratch;
    double* e = h + (ld + 1) * nm + nn * nkl;   // bra recurrence, one (k,l) at a time: e[j*nn+i]

    M[0] = constants::pi;
    for (int n = 1; n < nn; ++n) {
        M[n * nm] = c00 * M[(n - 1) * nm];
        if (n > 1) M[n * nm] += (n - 1) * b10 * M[(n - 2) * nm];
    }
    for (int n = 0; n < nn; ++n) {
        for (int m = 1; m < nm; ++m) {
            double v = c00p * M[n * nm + m - 1];
            if (m > 1) v += (m - 1) * b01 * M[n * nm + m - 2];
            if (n > 0) v += n * b00 * M[(n - 1) * nm + m - 1];
            M[n * nm + m] = v;
        }
    }

    // ket: powers of x2-D
    if (ld > 0) {
        for (int n = 0; n < nn; ++n) {
            for (int m = 0; m < nm; ++m) h[m] = M[n * nm + m];
            for (int l = 1; l <= ld; ++l)
                for (int m = 0; m < nm - l; ++m) h[l * nm + m] = h[(l - 1) * nm + m + 1] + CD * h[(l - 1) * nm + m];
            for (int k = 0; k <= lc; ++k)
                for (int l = 0; l <= ld; ++l) K[n * nkl + k * (ld + 1) + l] = h[l * nm + k];
        }
    }

    // bra: powers of x1-B
    if (lb > 0) {
        for (int kl = 0; kl < nkl; ++kl) {
            for (int i = 0; i < nn; ++i) e[i] = K[i * nkl + kl];
            for (int j = 1; j <= lb; ++j)
                for (int i = 0; i < nn - j; ++i) e[j * nn + i] = e[(j - 1) * nn + i + 1] + AB * e[(j - 1) * nn + i];
            for (int i = 0; i <= la; ++i)
                for (int j = 0; j <= lb; ++j) g[(i * (lb + 1) + j) * nkl + kl] = e[j * nn + i];
        }
    }
}

/// the x, y and z factors of one primitive quartet for one kernel term exp(-t r12^2)
///
/// p, q and t are the same on all three axes, so det M, omega and the covariance
/// of (x1,x2) are computed once; each axis then takes moments_1d. The tables g[x]
/// hold coulomb_1d's quadrature sums without the prefactor they share,
/// K_ab K_cd exp(-omega |P-Q|^2) det^{-3/2}, which is returned.
double coulomb_xyz(const QuartetLayout& L, const PrimitivePair& ab, const PrimitivePair& cd, const double t,
                   double* const g[3], double* scratch) {
    const double p = ab.p, q = cd.p;
    const double det = p * q + t * (p + q);
    const double rdet = 1.0 / det;
    const double omega = p * q * t * rdet;
    double r2 = 0.0;
    for (int x = 0; x < 3; ++x) r2 += (ab.P[x] - cd.P[x]) * (ab.P[x] - cd.P[x]);
    const double pref = ab.K * cd.K * std::exp(-omega * r2) * rdet / std::sqrt(det);

    // the covariance of (x1,x2) is M^{-1}/2, M = [[p+t, -t], [-t, q+t]]
    const double b10 = 0.5 * (q + t) * rdet, b01 = 0.5 * (p + t) * rdet, b00 = 0.5 * t * rdet;
    for (int x = 0; x < 3; ++x) {
        const double mu1 = ((q + t) * p * ab.P[x] + t * q * cd.P[x]) * rdet;
        const double mu2 = (t * p * ab.P[x] + (p + t) * q * cd.P[x]) * rdet;
        moments_1d(L, mu1 - L.A[x], mu2 - L.C[x], b10, b01, b00, L.A[x] - L.B[x], L.C[x] - L.D[x], g[x],
                   scratch);
    }
    return pref;
}

/// an (ss|ss) primitive quartet summed over the kernel terms
///
/// With l = 0 throughout, every axis contributes only its zeroth moment, pi, so
/// the sum over terms is the closed form
///   K_ab K_cd pi^3 sum_m w_m exp(-omega_m |P-Q|^2) det_m^{-3/2}.
double coulomb_ssss(const PrimitivePair& ab, const PrimitivePair& cd, const GaussianKernel& kernel,
                     const double tmax) {
    const double pq = ab.p * cd.p, s = ab.p + cd.p;
    double r2 = 0.0;
    for (int x = 0; x < 3; ++x) r2 += (ab.P[x] - cd.P[x]) * (ab.P[x] - cd.P[x]);
    double sum = 0.0;
    for (std::size_t m = 0; m < kernel.size(); ++m) {
        if (kernel.t[m] > tmax) continue;
        const double det = pq + kernel.t[m] * s;
        sum += kernel.w[m] * std::exp(-pq * kernel.t[m] / det * r2) / (det * std::sqrt(det));
    }
    return ab.K * cd.K * constants::pi * constants::pi * constants::pi * sum;
}

/// groups of shells that differ only in their contraction coefficients (same center, l and exponents)
///
/// A general contraction, as in the cc basis sets, is written as several shells
/// with the same primitives; their primitive integrals are the same.
std::vector<std::vector<int>> shell_groups(const std::vector<Shell>& shells) {
    std::vector<std::vector<int>> groups;
    for (int s = 0; s < int(shells.size()); ++s) {
        const Shell& sh = shells[s];
        const auto same = [&](const std::vector<int>& g) {
            const Shell& r = shells[g[0]];
            return r.center == sh.center and r.l == sh.l and r.expnt == sh.expnt;
        };
        const auto it = std::find_if(groups.begin(), groups.end(), same);
        if (it == groups.end()) groups.push_back({s});
        else it->push_back(s);
    }
    return groups;
}

/// pairs (a,b) of shells or groups are numbered a(a+1)/2+b for a >= b, whichever order they come in
inline std::size_t pair_index(const std::size_t a, const std::size_t b) {
    return (a >= b) ? a * (a + 1) / 2 + b : b * (b + 1) / 2 + a;
}

/// acc[n] = sum_m w_m (ab|exp(-t_m r12^2)|cd) for every component quartet n of one primitive quartet
///
/// With screen > 0 the short-range terms t_m > rho/(2 screen) are dropped, rho = pq/(p+q).
/// For t >> rho a term contributes about h rho/t of an (ss|ss) integral (h the
/// logarithmic spacing of the fit), so the dropped tail is about screen of it.
void primitive_quartet(const QuartetLayout& L, const PrimitivePair& ab, const PrimitivePair& cd,
                       const GaussianKernel& kernel, const double screen, std::vector<double> (&g)[3],
                       double* scratch, double* acc) {
    const double tmax = (screen > 0.0) ? 0.5 * ab.p * cd.p / ((ab.p + cd.p) * screen)
                                       : std::numeric_limits<double>::infinity();
    if (L.la + L.lb + L.lc + L.ld == 0) {
        acc[0] = coulomb_ssss(ab, cd, kernel, tmax);
        return;
    }
    const int nblock = L.nblock();
    std::fill(acc, acc + nblock, 0.0);
    double* const gp[3] = {g[0].data(), g[1].data(), g[2].data()};
    for (std::size_t m = 0; m < kernel.size(); ++m) {
        if (kernel.t[m] > tmax) continue;
        const double f = kernel.w[m] * coulomb_xyz(L, ab, cd, kernel.t[m], gp, scratch);
        for (int n = 0; n < nblock; ++n) {
            const std::array<int,3>& i = L.index[n];
            acc[n] += f * gp[0][i[0]] * gp[1][i[1]] * gp[2][i[2]];
        }
    }
}

/// the integral blocks of some shell quartets (s1 s2|s3 s4) of one group quartet, one block after the other
///
/// bra and ket are the primitive pairs of the two group pairs. Each primitive quartet is
/// computed once and contracted into every member; a block holds the component quartets
/// (i j|k l) at ((i nb + j) nc + k) nd + l.
std::vector<double> quartet_blocks(const std::vector<Shell>& shells, const std::vector<std::array<int,4>>& members,
                                   const std::vector<PrimitivePair>& bra, const std::vector<PrimitivePair>& ket,
                                   const GaussianKernel& kernel, const double screen) {
    const std::array<int,4>& f = members.front();
    const QuartetLayout L(shells[f[0]], shells[f[1]], shells[f[2]], shells[f[3]]);
    const int nblock = L.nblock();
    std::vector<double> g[3], scratch(L.nscratch), acc(nblock), blocks(members.size() * nblock, 0.0);
    for (auto& v : g) v.resize(L.ntab);
    for (const PrimitivePair& pab : bra) {
        for (const PrimitivePair& pcd : ket) {
            primitive_quartet(L, pab, pcd, kernel, screen, g, scratch.data(), acc.data());
            for (std::size_t q = 0; q < members.size(); ++q) {
                const std::array<int,4>& s = members[q];
                const double coef = shells[s[0]].coeff[pab.ia] * shells[s[1]].coeff[pab.ib]
                                  * shells[s[2]].coeff[pcd.ia] * shells[s[3]].coeff[pcd.ib];
                double* block = blocks.data() + q * nblock;
                for (int n = 0; n < nblock; ++n) block[n] += coef * acc[n];
            }
        }
    }
    return blocks;
}

/// all shell quartets of one group quartet (G1 G2|G3 G4), written into G
///
/// The shell quartets are restricted so that each one is computed exactly once
/// over all canonical group quartets, as in the canonical loop over shells:
/// s1 >= s2 within one group, s3 >= s4 within one group, and bra >= ket when the
/// bra and ket groups are the same.
void group_quartet(const std::vector<Shell>& shells, const std::vector<int>& G1, const std::vector<int>& G2,
                   const std::vector<int>& G3, const std::vector<int>& G4, const bool same12, const bool same34,
                   const bool samebraket, const std::vector<PrimitivePair>& bra, const std::vector<PrimitivePair>& ket,
                   const GaussianKernel& kernel, const double screen, PackedERI& G) {
    std::vector<std::array<int,4>> members;
    for (const int s1 : G1)
        for (const int s2 : G2) {
            if (same12 and s1 < s2) continue;
            for (const int s3 : G3)
                for (const int s4 : G4) {
                    if (same34 and s3 < s4) continue;
                    if (samebraket and pair_index(s1, s2) < pair_index(s3, s4)) continue;
                    members.push_back({s1, s2, s3, s4});
                }
        }
    if (members.empty()) return;

    const std::vector<double> blocks = quartet_blocks(shells, members, bra, ket, kernel, screen);
    const std::size_t nblock = blocks.size() / members.size();
    for (std::size_t q = 0; q < members.size(); ++q) {
        const std::array<int,4>& s = members[q];
        scatter_quartet(G, shells[s[0]], shells[s[1]], shells[s[2]], shells[s[3]], blocks.data() + q * nblock);
    }
}

/// a pair of shells of a group pair, and the row of its first function pair within the group pair
struct ShellPair {
    int s1, s2;
    std::size_t first;
};

/// the function pair (i, j) of a shell pair, numbered within the shell pair: i >= j for a shell with itself
inline std::size_t component_pair(const int i, const int j, const int nj, const bool same) {
    return same ? std::size_t(i) * (i + 1) / 2 + j : std::size_t(i) * nj + j;
}

} // namespace


/// what all two-electron integrals of a basis share: the shell groups, and per group pair
/// (a >= b, numbered pair_index(a, b)) the two groups, their primitive pairs and their shell pairs
struct ShellGroupData {
    std::vector<std::vector<int>> groups;
    std::vector<std::array<std::size_t,2>> group_pairs;         ///< (a, b)
    std::vector<std::vector<PrimitivePair>> pairs;              ///< primitive pairs of group a with group b
    std::vector<std::vector<ShellPair>> shell_pairs;            ///< s1 of a, s2 of b, s1 >= s2 if a = b; in row order

    explicit ShellGroupData(const std::vector<Shell>& shells) : groups(shell_groups(shells)) {
        const std::size_t ng = groups.size(), ngp = ng * (ng + 1) / 2;
        group_pairs.resize(ngp);
        pairs.resize(ngp);
        shell_pairs.resize(ngp);
        for (std::size_t a = 0; a < ng; ++a) {
            for (std::size_t b = 0; b <= a; ++b) {
                const std::size_t ab = pair_index(a, b);
                group_pairs[ab] = {a, b};
                pairs[ab] = primitive_pairs(shells[groups[a][0]], shells[groups[b][0]]);
                std::size_t first = 0;
                for (const int s1 : groups[a])
                    for (const int s2 : groups[b]) {
                        if (a == b and s1 < s2) continue;
                        shell_pairs[ab].push_back({s1, s2, first});
                        const std::size_t n1 = shells[s1].ncart(), n2 = shells[s2].ncart();
                        first += (s1 == s2) ? n1 * (n1 + 1) / 2 : n1 * n2;
                    }
            }
        }
    }
};


namespace {

/// what the tasks of SeparatedGaussianIntegrals::eri share: shells, their groups, the kernel
struct GroupQuartetData {
    const std::vector<Shell>& shells;
    const ShellGroupData& sg;
    const GaussianKernel& kernel;
    double screen;                                      ///< see primitive_quartet
    double schwarz;                                     ///< skip (ab|cd) if q[ab] q[cd] < schwarz
    std::vector<double> q;                              ///< per group pair, sqrt of the largest |(mu nu|mu nu)|
    std::atomic<long> computed{0}, skipped{0};

    /// compute the group quartet (ab|cd), (a b) >= (c d)
    void compute(const std::size_t a, const std::size_t b, const std::size_t c, const std::size_t d,
                 PackedERI& G) {
        const std::size_t ab = pair_index(a, b), cd = pair_index(c, d);
        group_quartet(shells, sg.groups[a], sg.groups[b], sg.groups[c], sg.groups[d], a == b, c == d, ab == cd,
                      sg.pairs[ab], sg.pairs[cd], kernel, screen, G);
        ++computed;
    }
};

/// a task on the thread pool: the group quartets with bra group pair (a b)
///
/// With diagonal, only (ab|ab), whose integrals give the Schwarz factors;
/// otherwise all (ab|cd) with (c d) < (a b) that pass the Schwarz test. Group
/// quartets of different tasks write disjoint elements of G, so the tasks need
/// no locks, and every element is computed in the order of the serial loop.
class GroupQuartetTask : public TaskInterface {
public:
    GroupQuartetTask(GroupQuartetData& data, const std::size_t a, const std::size_t b, const bool diagonal,
                     PackedERI& G)
        : data_(data), a_(a), b_(b), diagonal_(diagonal), G_(G) {}

    using TaskInterface::run;
    void run(World&) override {
        if (diagonal_) {
            data_.compute(a_, b_, a_, b_, G_);
            return;
        }
        const std::size_t ab = pair_index(a_, b_);
        for (std::size_t c = 0; c <= a_; ++c) {
            for (std::size_t d = 0; d <= c; ++d) {
                const std::size_t cd = pair_index(c, d);
                if (cd >= ab) continue;     // (ab|cd) = (cd|ab); the diagonal is done
                if (data_.q[ab] * data_.q[cd] < data_.schwarz) {
                    ++data_.skipped;
                    continue;
                }
                data_.compute(a_, b_, c, d, G_);
            }
        }
    }

private:
    GroupQuartetData& data_;
    const std::size_t a_, b_;
    const bool diagonal_;
    PackedERI& G_;
};

/// boundaries of at most nchunk ranges of consecutive items with about equal total cost
std::vector<std::size_t> equal_cost_chunks(const std::vector<double>& cost, std::size_t nchunk) {
    const std::size_t n = cost.size();
    nchunk = std::max<std::size_t>(1, std::min(nchunk, n));
    double total = 0.0;
    for (const double c : cost) total += c;
    std::vector<std::size_t> bounds{0};
    double sum = 0.0;
    for (std::size_t i = 0; i < n; ++i) {
        sum += cost[i];
        if (bounds.size() < nchunk and sum >= total * double(bounds.size()) / double(nchunk)) bounds.push_back(i + 1);
    }
    if (bounds.back() != n) bounds.push_back(n);
    return bounds;
}

/// the number of tasks a call of pair_diagonal or eri_columns splits into, for load balance
std::size_t task_count() { return 8 * (ThreadPool::size() + 1); }

/// a task on the thread pool: D(row) = (mu nu|mu nu) for the rows of the group pairs [begin, end)
class DiagonalTask : public TaskInterface {
public:
    DiagonalTask(const std::vector<Shell>& shells, const ShellGroupData& sg, const GaussianKernel& kernel,
                 const FunctionPairs& fp, const std::size_t begin, const std::size_t end, double* D)
        : shells_(shells), sg_(sg), kernel_(kernel), fp_(fp), begin_(begin), end_(end), D_(D) {}

    using TaskInterface::run;
    void run(World&) override {
        for (std::size_t ab = begin_; ab < end_; ++ab) {
            const std::vector<ShellPair>& sps = sg_.shell_pairs[ab];
            std::vector<std::array<int,4>> members;
            for (const ShellPair& sp : sps) members.push_back({sp.s1, sp.s2, sp.s1, sp.s2});
            const std::vector<double> blocks = quartet_blocks(shells_, members, sg_.pairs[ab], sg_.pairs[ab],
                                                              kernel_, 0.0);
            const std::size_t nblock = blocks.size() / members.size();
            for (std::size_t q = 0; q < sps.size(); ++q) {
                const int na = shells_[sps[q].s1].ncart(), nb = shells_[sps[q].s2].ncart();
                const bool same = sps[q].s1 == sps[q].s2;
                const double* block = blocks.data() + q * nblock;
                double* d = D_ + fp_.first[ab] + sps[q].first;
                for (int i = 0; i < na; ++i)
                    for (int j = 0; j < (same ? i + 1 : nb); ++j)
                        d[component_pair(i, j, nb, same)] = block[((i * nb + j) * na + i) * nb + j];
            }
        }
    }

private:
    const std::vector<Shell>& shells_;
    const ShellGroupData& sg_;
    const GaussianKernel& kernel_;
    const FunctionPairs& fp_;
    const std::size_t begin_, end_;
    double* D_;
};

/// a task on the thread pool: the columns of ket group pair g for the bra group pairs bra[begin, end)
///
/// W(j, r) is written at W[j ldw + r], r counted from rowstart[k] for bra[k]. Every
/// pair of bra and ket shell pairs is computed, without the canonical restriction of
/// eri(); the duplicates i < j (s1 = s2) and k < l (s3 = s4) are not written.
class ColumnTask : public TaskInterface {
public:
    ColumnTask(const std::vector<Shell>& shells, const ShellGroupData& sg, const GaussianKernel& kernel,
               const std::vector<std::size_t>& bra, const std::vector<std::size_t>& rowstart, const std::size_t begin,
               const std::size_t end, const std::size_t g, double* W, const std::size_t ldw)
        : shells_(shells), sg_(sg), kernel_(kernel), bra_(bra), rowstart_(rowstart), begin_(begin), end_(end), g_(g),
          W_(W), ldw_(ldw) {}

    using TaskInterface::run;
    void run(World&) override {
        const std::vector<ShellPair>& kets = sg_.shell_pairs[g_];
        for (std::size_t k = begin_; k < end_; ++k) {
            const std::size_t ab = bra_[k];
            const std::vector<ShellPair>& bras = sg_.shell_pairs[ab];
            std::vector<std::array<int,4>> members;
            for (const ShellPair& b : bras)
                for (const ShellPair& c : kets) members.push_back({b.s1, b.s2, c.s1, c.s2});
            const std::vector<double> blocks = quartet_blocks(shells_, members, sg_.pairs[ab], sg_.pairs[g_],
                                                              kernel_, 0.0);
            const std::size_t nblock = blocks.size() / members.size();
            for (std::size_t q = 0; q < members.size(); ++q) {
                const ShellPair& b = bras[q / kets.size()];
                const ShellPair& c = kets[q % kets.size()];
                const int na = shells_[b.s1].ncart(), nb = shells_[b.s2].ncart();
                const int nc = shells_[c.s1].ncart(), nd = shells_[c.s2].ncart();
                const bool same12 = b.s1 == b.s2, same34 = c.s1 == c.s2;
                const double* block = blocks.data() + q * nblock;
                for (int i = 0; i < na; ++i) {
                    for (int j = 0; j < (same12 ? i + 1 : nb); ++j) {
                        const std::size_t r = rowstart_[k] + b.first + component_pair(i, j, nb, same12);
                        for (int kk = 0; kk < nc; ++kk)
                            for (int l = 0; l < (same34 ? kk + 1 : nd); ++l)
                                W_[(c.first + component_pair(kk, l, nd, same34)) * ldw_ + r] =
                                    block[((i * nb + j) * nc + kk) * nd + l];
                    }
                }
            }
        }
    }

private:
    const std::vector<Shell>& shells_;
    const ShellGroupData& sg_;
    const GaussianKernel& kernel_;
    const std::vector<std::size_t>& bra_;
    const std::vector<std::size_t>& rowstart_;
    const std::size_t begin_, end_, g_;
    double* W_;
    const std::size_t ldw_;
};

} // namespace


std::vector<std::array<int,3>> cartesian_components(const int l) {
    std::vector<std::array<int,3>> c;
    for (int lx = l; lx >= 0; --lx)
        for (int ly = l - lx; ly >= 0; --ly) c.push_back({lx, ly, l - lx - ly});
    return c;
}


std::vector<Shell> make_shells(const Molecule& molecule, const AtomicBasisSet& aobasis) {
    std::vector<Shell> shells;
    const int nbf = aobasis.nbf(molecule);
    for (int ibf = 0; ibf < nbf; ++ibf) {
        const AtomicBasisFunction abf = aobasis.get_atomic_basis_function(molecule, ibf);
        if (abf.get_index() != 0) continue;     // only the first function of each shell
        const ContractedGaussianShell& g = abf.get_shell();
        Shell s;
        abf.get_coords(s.center[0], s.center[1], s.center[2]);
        s.atom = aobasis.basisfn_to_atom(molecule, ibf);
        s.l = g.angular_momentum();
        s.expnt = g.get_expnt();
        s.coeff = g.get_coeff();
        s.offset = ibf;
        MADNESS_CHECK_THROW(g.nbf() == s.ncart(), "make_shells: unexpected number of functions in a shell");
        MADNESS_CHECK_THROW(s.l + 2 < maxpow, "make_shells: angular momentum too high");
        shells.push_back(s);
    }
    return shells;
}


GaussHermiteRule::GaussHermiteRule(const int nmax) : nodes_(nmax + 1), weights_(nmax + 1) {
    // Golub-Welsch: the nodes are the eigenvalues of the Jacobi matrix of the
    // monic Hermite polynomials, the weights sqrt(pi) times the squared first
    // components of the normalized eigenvectors
    for (int n = 1; n <= nmax; ++n) {
        Tensor<double> J(n, n);
        for (int k = 1; k < n; ++k) J(k - 1, k) = J(k, k - 1) = std::sqrt(0.5 * k);
        Tensor<double> V, e;
        syev(J, V, e);
        nodes_[n].resize(n);
        weights_[n].resize(n);
        for (int k = 0; k < n; ++k) {
            nodes_[n][k] = e(k);
            weights_[n][k] = std::sqrt(constants::pi) * V(0, k) * V(0, k);
        }
    }
}

const std::vector<double>& GaussHermiteRule::nodes(const int n) const {
    MADNESS_CHECK_THROW(n >= 1 and n <= nmax(), "GaussHermiteRule: no rule with that many nodes");
    return nodes_[n];
}

const std::vector<double>& GaussHermiteRule::weights(const int n) const {
    MADNESS_CHECK_THROW(n >= 1 and n <= nmax(), "GaussHermiteRule: no rule with that many nodes");
    return weights_[n];
}


GaussianKernel GaussianKernel::coulomb(const double lo, const double hi, const double eps) {
    const GFit<double,3> fit = GFit<double,3>::CoulombFit(lo, hi, eps, false);
    const Tensor<double> c = fit.coeffs();
    const Tensor<double> e = fit.exponents();
    GaussianKernel k;
    for (long m = 0; m < c.size(); ++m) {
        k.w.push_back(c(m));
        k.t.push_back(e(m));
    }
    return k;
}

double GaussianKernel::operator()(const double r) const {
    double f = 0.0;
    for (std::size_t m = 0; m < size(); ++m) f += w[m] * std::exp(-t[m] * r * r);
    return f;
}

double GaussianKernel::max_relative_coulomb_error(const double lo, const double hi, const int npt) const {
    double err = 0.0;
    for (int i = 0; i < npt; ++i) {
        const double r = lo * std::pow(hi / lo, double(i) / (npt - 1));
        err = std::max(err, std::abs((*this)(r) * r - 1.0));
    }
    return err;
}


PackedERI::PackedERI(const long nbf) : nbf_(nbf) {
    const std::size_t npair = pair(nbf - 1, nbf - 1) + 1;
    data_.assign(npair * (npair + 1) / 2, 0.0);
}


SeparatedGaussianIntegrals::SeparatedGaussianIntegrals(const std::vector<Shell>& shells,
                                                       const GaussianKernel& coulomb)
    : shells_(shells), coulomb_(coulomb), gh_(16), groups_(std::make_shared<const ShellGroupData>(shells_)) {
    for (const Shell& s : shells_) nbf_ = std::max(nbf_, long(s.offset + s.ncart()));
}


Tensor<double> SeparatedGaussianIntegrals::overlap() const {
    Tensor<double> S(nbf_, nbf_);
    for (std::size_t a = 0; a < shells_.size(); ++a) {
        for (std::size_t b = 0; b <= a; ++b) {
            const Shell& sa = shells_[a];
            const Shell& sb = shells_[b];
            const auto ca = cartesian_components(sa.l), cb = cartesian_components(sb.l);
            const int nj = sb.l + 1;
            std::vector<double> s[3];
            for (auto& v : s) v.resize((sa.l + 1) * nj);
            for (std::size_t ia = 0; ia < sa.expnt.size(); ++ia) {
                for (std::size_t ib = 0; ib < sb.expnt.size(); ++ib) {
                    const double cc = sa.coeff[ia] * sb.coeff[ib];
                    for (int d = 0; d < 3; ++d)
                        gaussian_1d(sa.l, sb.l, sa.expnt[ia], sa.center[d], sb.expnt[ib], sb.center[d],
                                    0.0, 0.0, gh_, s[d].data());
                    for (std::size_t i = 0; i < ca.size(); ++i)
                        for (std::size_t j = 0; j < cb.size(); ++j)
                            S(sa.offset + i, sb.offset + j) += cc * s[0][ca[i][0] * nj + cb[j][0]]
                                                                  * s[1][ca[i][1] * nj + cb[j][1]]
                                                                  * s[2][ca[i][2] * nj + cb[j][2]];
                }
            }
            symmetrize_block(S, sa, sb);
        }
    }
    return S;
}


Tensor<double> SeparatedGaussianIntegrals::kinetic() const {
    Tensor<double> T(nbf_, nbf_);
    for (std::size_t a = 0; a < shells_.size(); ++a) {
        for (std::size_t b = 0; b <= a; ++b) {
            const Shell& sa = shells_[a];
            const Shell& sb = shells_[b];
            const auto ca = cartesian_components(sa.l), cb = cartesian_components(sb.l);
            const int ns = sb.l + 3;            // overlaps up to j = lb+2, for the second derivative
            const int nj = sb.l + 1;
            std::vector<double> s[3], k[3];
            for (int d = 0; d < 3; ++d) {
                s[d].resize((sa.l + 1) * ns);
                k[d].resize((sa.l + 1) * nj);
            }
            for (std::size_t ia = 0; ia < sa.expnt.size(); ++ia) {
                for (std::size_t ib = 0; ib < sb.expnt.size(); ++ib) {
                    const double cc = sa.coeff[ia] * sb.coeff[ib];
                    const double beta = sb.expnt[ib];
                    for (int d = 0; d < 3; ++d) {
                        gaussian_1d(sa.l, sb.l + 2, sa.expnt[ia], sa.center[d], beta, sb.center[d],
                                    0.0, 0.0, gh_, s[d].data());
                        // -1/2 d^2/dx^2 (x-B)^j exp(-beta(x-B)^2), expressed through the overlaps
                        for (int i = 0; i <= sa.l; ++i) {
                            for (int j = 0; j <= sb.l; ++j) {
                                const double* si = s[d].data() + i * ns;
                                double v = -2.0 * beta * (2 * j + 1) * si[j] + 4.0 * beta * beta * si[j + 2];
                                if (j >= 2) v += j * (j - 1) * si[j - 2];
                                k[d][i * nj + j] = -0.5 * v;
                            }
                        }
                    }
                    for (std::size_t i = 0; i < ca.size(); ++i) {
                        for (std::size_t j = 0; j < cb.size(); ++j) {
                            const int ix = ca[i][0] * ns + cb[j][0], iy = ca[i][1] * ns + cb[j][1],
                                      iz = ca[i][2] * ns + cb[j][2];
                            const int kx = ca[i][0] * nj + cb[j][0], ky = ca[i][1] * nj + cb[j][1],
                                      kz = ca[i][2] * nj + cb[j][2];
                            const double t = k[0][kx] * s[1][iy] * s[2][iz]
                                           + s[0][ix] * k[1][ky] * s[2][iz]
                                           + s[0][ix] * s[1][iy] * k[2][kz];
                            T(sa.offset + i, sb.offset + j) += cc * t;
                        }
                    }
                }
            }
            symmetrize_block(T, sa, sb);
        }
    }
    return T;
}


Tensor<double> SeparatedGaussianIntegrals::nuclear_attraction(const Molecule& molecule) const {
    Tensor<double> V(nbf_, nbf_);
    for (std::size_t a = 0; a < shells_.size(); ++a) {
        for (std::size_t b = 0; b <= a; ++b) {
            const Shell& sa = shells_[a];
            const Shell& sb = shells_[b];
            const auto ca = cartesian_components(sa.l), cb = cartesian_components(sb.l);
            const int nj = sb.l + 1;
            std::vector<double> g[3];
            for (auto& v : g) v.resize((sa.l + 1) * nj);
            for (std::size_t ia = 0; ia < sa.expnt.size(); ++ia) {
                for (std::size_t ib = 0; ib < sb.expnt.size(); ++ib) {
                    const double cc = sa.coeff[ia] * sb.coeff[ib];
                    for (std::size_t iat = 0; iat < molecule.natom(); ++iat) {
                        const Atom& atom = molecule.get_atom(iat);
                        const double C[3] = {atom.x, atom.y, atom.z};
                        for (std::size_t m = 0; m < coulomb_.size(); ++m) {
                            for (int d = 0; d < 3; ++d)
                                gaussian_1d(sa.l, sb.l, sa.expnt[ia], sa.center[d], sb.expnt[ib], sb.center[d],
                                            coulomb_.t[m], C[d], gh_, g[d].data());
                            const double f = -atom.q * coulomb_.w[m] * cc;
                            for (std::size_t i = 0; i < ca.size(); ++i)
                                for (std::size_t j = 0; j < cb.size(); ++j)
                                    V(sa.offset + i, sb.offset + j) += f * g[0][ca[i][0] * nj + cb[j][0]]
                                                                         * g[1][ca[i][1] * nj + cb[j][1]]
                                                                         * g[2][ca[i][2] * nj + cb[j][2]];
                        }
                    }
                }
            }
            symmetrize_block(V, sa, sb);
        }
    }
    return V;
}


PackedERI SeparatedGaussianIntegrals::eri(World& world, const double screen, const double schwarz,
                                          ERIStats* stats) const {
    PackedERI G(nbf_);
    const ShellGroupData& sg = *groups_;
    GroupQuartetData data{shells_, sg, coulomb_, screen, schwarz, {}};
    const std::size_t ng = sg.groups.size();

    // first the diagonal group quartets (ab|ab), for the Schwarz factors
    for (std::size_t a = ng; a-- > 0;)
        for (std::size_t b = a + 1; b-- > 0;) world.taskq.add(new GroupQuartetTask(data, a, b, true, G));
    world.taskq.fence();
    data.q.assign(sg.pairs.size(), 0.0);
    for (std::size_t a = 0; a < ng; ++a) {
        for (std::size_t b = 0; b <= a; ++b) {
            double qmax = 0.0;
            for (const int s1 : sg.groups[a])
                for (const int s2 : sg.groups[b])
                    for (int i = 0; i < shells_[s1].ncart(); ++i)
                        for (int j = 0; j < shells_[s2].ncart(); ++j) {
                            const long mu = shells_[s1].offset + i, nu = shells_[s2].offset + j;
                            qmax = std::max(qmax, std::abs(G(mu, nu, mu, nu)));
                        }
            data.q[pair_index(a, b)] = std::sqrt(qmax);
        }
    }

    // then the rest; the last bra pairs have the most ket pairs, so submit them first
    for (std::size_t a = ng; a-- > 0;)
        for (std::size_t b = a + 1; b-- > 0;) world.taskq.add(new GroupQuartetTask(data, a, b, false, G));
    world.taskq.fence();
    if (stats) {
        stats->computed = data.computed;
        stats->skipped = data.skipped;
    }
    return G;
}


FunctionPairs SeparatedGaussianIntegrals::function_pairs() const {
    const ShellGroupData& sg = *groups_;
    const std::size_t ngp = sg.group_pairs.size();
    const std::size_t npair = std::size_t(nbf_) * (nbf_ + 1) / 2;
    FunctionPairs fp;
    fp.first.resize(ngp);
    fp.count.resize(ngp);
    fp.functions.reserve(npair);
    fp.group_pair.reserve(npair);
    for (std::size_t ab = 0; ab < ngp; ++ab) {
        fp.first[ab] = fp.functions.size();
        for (const ShellPair& sp : sg.shell_pairs[ab]) {
            const Shell& s1 = shells_[sp.s1];
            const Shell& s2 = shells_[sp.s2];
            for (int i = 0; i < s1.ncart(); ++i)
                for (int j = 0; j < ((sp.s1 == sp.s2) ? i + 1 : s2.ncart()); ++j) {
                    fp.functions.push_back({s1.offset + i, s2.offset + j});
                    fp.group_pair.push_back(ab);
                }
        }
        fp.count[ab] = fp.functions.size() - fp.first[ab];
    }
    MADNESS_CHECK_THROW(fp.size() == npair, "function_pairs: the pairs do not cover the basis");
    return fp;
}


Tensor<double> SeparatedGaussianIntegrals::pair_diagonal(World& world, const FunctionPairs& fp) const {
    const ShellGroupData& sg = *groups_;
    MADNESS_CHECK_THROW(fp.ngroup_pairs() == sg.group_pairs.size() and fp.size() == std::size_t(nbf_) * (nbf_ + 1) / 2,
                        "pair_diagonal: the function pairs belong to another basis");
    std::vector<double> cost(fp.ngroup_pairs());
    for (std::size_t ab = 0; ab < cost.size(); ++ab)
        cost[ab] = double(sg.pairs[ab].size()) * double(sg.pairs[ab].size()) * double(fp.count[ab]);
    Tensor<double> D(long(fp.size()));
    const std::vector<std::size_t> chunk = equal_cost_chunks(cost, task_count());
    for (std::size_t c = 0; c + 1 < chunk.size(); ++c)
        world.taskq.add(new DiagonalTask(shells_, sg, coulomb_, fp, chunk[c], chunk[c + 1], D.ptr()));
    world.taskq.fence();
    return D;
}


Tensor<double> SeparatedGaussianIntegrals::eri_columns(World& world, const FunctionPairs& fp, const std::size_t g,
                                                       const std::vector<std::size_t>& bra_group_pairs) const {
    const ShellGroupData& sg = *groups_;
    MADNESS_CHECK_THROW(fp.ngroup_pairs() == sg.group_pairs.size() and fp.size() == std::size_t(nbf_) * (nbf_ + 1) / 2,
                        "eri_columns: the function pairs belong to another basis");
    MADNESS_CHECK_THROW(g < fp.ngroup_pairs(), "eri_columns: no such ket group pair");
    const std::vector<std::size_t>& bra = bra_group_pairs;
    std::vector<std::size_t> rowstart(bra.size() + 1, 0);
    std::vector<double> cost(bra.size());
    for (std::size_t k = 0; k < bra.size(); ++k) {
        MADNESS_CHECK_THROW(bra[k] < fp.ngroup_pairs(), "eri_columns: no such bra group pair");
        rowstart[k + 1] = rowstart[k] + fp.count[bra[k]];
        cost[k] = double(sg.pairs[bra[k]].size()) * double(fp.count[bra[k]]);
    }
    Tensor<double> W(long(fp.count[g]), long(rowstart.back()));
    const std::vector<std::size_t> chunk = equal_cost_chunks(cost, task_count());
    for (std::size_t c = 0; c + 1 < chunk.size(); ++c)
        world.taskq.add(new ColumnTask(shells_, sg, coulomb_, bra, rowstart, chunk[c], chunk[c + 1], g, W.ptr(),
                                       rowstart.back()));
    world.taskq.fence();
    return W;
}


Tensor<double> SeparatedGaussianIntegrals::eri_reference() const {
    Tensor<double> G(nbf_, nbf_, nbf_, nbf_);
    const std::size_t ns = shells_.size();
    for (std::size_t a = 0; a < ns; ++a) {
        for (std::size_t b = 0; b <= a; ++b) {
            const std::size_t ab = a * (a + 1) / 2 + b;
            for (std::size_t c = 0; c < ns; ++c) {
                for (std::size_t d = 0; d <= c; ++d) {
                    if (c * (c + 1) / 2 + d > ab) continue;     // (ab|cd) = (cd|ab)
                    const Shell& sa = shells_[a];
                    const Shell& sb = shells_[b];
                    const Shell& sc = shells_[c];
                    const Shell& sd = shells_[d];
                    const auto ca = cartesian_components(sa.l), cb = cartesian_components(sb.l);
                    const auto cc = cartesian_components(sc.l), cd = cartesian_components(sd.l);
                    const int na = sa.ncart(), nb = sb.ncart(), nc = sc.ncart(), nd = sd.ncart();
                    const int lb1 = sb.l + 1, lc1 = sc.l + 1, ld1 = sd.l + 1;
                    const int ntab = (sa.l + 1) * lb1 * lc1 * ld1;
                    const int nblock = na * nb * nc * nd;

                    // per component quartet, the position of its powers in the x, y and z tables
                    std::vector<std::array<int,3>> index(nblock);
                    for (int i = 0, n = 0; i < na; ++i)
                        for (int j = 0; j < nb; ++j)
                            for (int k = 0; k < nc; ++k)
                                for (int l = 0; l < nd; ++l, ++n)
                                    for (int x = 0; x < 3; ++x)
                                        index[n][x] = ((ca[i][x] * lb1 + cb[j][x]) * lc1 + cc[k][x]) * ld1 + cd[l][x];

                    std::vector<double> block(nblock, 0.0), acc(nblock);
                    std::vector<double> g[3];
                    for (auto& v : g) v.resize(ntab);
                    for (std::size_t ia = 0; ia < sa.expnt.size(); ++ia) {
                        for (std::size_t ib = 0; ib < sb.expnt.size(); ++ib) {
                            const GaussianProduct1D pab[3] = {
                                {sa.expnt[ia], sa.center[0], sb.expnt[ib], sb.center[0]},
                                {sa.expnt[ia], sa.center[1], sb.expnt[ib], sb.center[1]},
                                {sa.expnt[ia], sa.center[2], sb.expnt[ib], sb.center[2]}};
                            for (std::size_t ic = 0; ic < sc.expnt.size(); ++ic) {
                                for (std::size_t id = 0; id < sd.expnt.size(); ++id) {
                                    const GaussianProduct1D pcd[3] = {
                                        {sc.expnt[ic], sc.center[0], sd.expnt[id], sd.center[0]},
                                        {sc.expnt[ic], sc.center[1], sd.expnt[id], sd.center[1]},
                                        {sc.expnt[ic], sc.center[2], sd.expnt[id], sd.center[2]}};
                                    std::fill(acc.begin(), acc.end(), 0.0);
                                    for (std::size_t m = 0; m < coulomb_.size(); ++m) {
                                        for (int x = 0; x < 3; ++x)
                                            coulomb_1d(sa.l, sb.l, sc.l, sd.l, pab[x], sa.center[x], sb.center[x],
                                                       pcd[x], sc.center[x], sd.center[x], coulomb_.t[m], gh_,
                                                       g[x].data());
                                        const double wm = coulomb_.w[m];
                                        for (int n = 0; n < nblock; ++n)
                                            acc[n] += wm * g[0][index[n][0]] * g[1][index[n][1]] * g[2][index[n][2]];
                                    }
                                    const double coef = sa.coeff[ia] * sb.coeff[ib] * sc.coeff[ic] * sd.coeff[id];
                                    for (int n = 0; n < nblock; ++n) block[n] += coef * acc[n];
                                }
                            }
                        }
                    }

                    scatter_quartet(G, sa, sb, sc, sd, block.data());
                }
            }
        }
    }
    return G;
}

} // namespace lcao
} // namespace madness
