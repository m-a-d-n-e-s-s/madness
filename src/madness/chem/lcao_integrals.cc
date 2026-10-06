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
#include <cmath>

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

/// Gauss-Hermite nodes per axis that the two-electron 1D factors can need (GaussHermiteRule's default)
constexpr int maxnode = 16;

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
    int n;                      ///< Gauss-Hermite nodes per axis
    const double* y;            ///< the nodes
    double ww[maxnode * maxnode];   ///< products of the weights, w[k1] w[k2] at k1*n+k2
    int ntab;                   ///< size of one 1D table
    std::vector<std::array<int,3>> index;   ///< per component quartet, its position in the x, y and z tables

    QuartetLayout(const Shell& sa, const Shell& sb, const Shell& sc, const Shell& sd, const GaussHermiteRule& gh)
        : la(sa.l), lb(sb.l), lc(sc.l), ld(sd.l), A(sa.center), B(sb.center), C(sc.center), D(sd.center),
          n((sa.l + sb.l + sc.l + sd.l) / 2 + 1), y(gh.nodes(n).data()),
          ntab((sa.l + 1) * (sb.l + 1) * (sc.l + 1) * (sd.l + 1)) {
        MADNESS_CHECK_THROW(n <= maxnode, "QuartetLayout: angular momentum too high");
        const std::vector<double>& w = gh.weights(n);
        for (int k1 = 0; k1 < n; ++k1)
            for (int k2 = 0; k2 < n; ++k2) ww[k1 * n + k2] = w[k1] * w[k2];
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

/// the x, y and z factors of one primitive quartet for one kernel term exp(-t r12^2)
///
/// The quadrature of coulomb_1d, but p, q and t are the same on all three axes, so
/// det M, omega, the Cholesky factor and the node offsets are computed once. The
/// tables g[x] get the quadrature sums only; the prefactor they share,
/// K_ab K_cd exp(-omega |P-Q|^2) det^{-3/2}, is returned.
double coulomb_xyz(const QuartetLayout& L, const PrimitivePair& ab, const PrimitivePair& cd, const double t,
                   double* const g[3]) {
    const double p = ab.p, q = cd.p;
    const double det = p * q + t * (p + q);
    const double rdet = 1.0 / det;
    const double sdet = std::sqrt(det);
    const double l11 = std::sqrt(p + t);
    const double il11 = 1.0 / l11;
    const double il22 = l11 / sdet;                 // 1/l22, with l22 = sqrt(det/(p+t))
    const double c12 = t * il11 * il11 * il22;      // t / ((p+t) l22)
    const double omega = p * q * t * rdet;
    double r2 = 0.0;
    for (int x = 0; x < 3; ++x) r2 += (ab.P[x] - cd.P[x]) * (ab.P[x] - cd.P[x]);
    const double pref = ab.K * cd.K * std::exp(-omega * r2) * rdet / sdet;

    // offsets of the nodes from the mean, the same on every axis
    const int n = L.n;
    double o1[maxnode * maxnode], o2[maxnode];
    for (int k2 = 0; k2 < n; ++k2) o2[k2] = L.y[k2] * il22;
    for (int k1 = 0; k1 < n; ++k1)
        for (int k2 = 0; k2 < n; ++k2) o1[k1 * n + k2] = L.y[k1] * il11 + c12 * L.y[k2];

    const int nb = L.lb + 1, nd = L.ld + 1;
    const int nbra = (L.la + 1) * nb, nket = (L.lc + 1) * nd;
    double pa[maxpow], pb[maxpow], pc[maxpow], pd[maxpow], u[maxpow * maxpow], v[maxpow * maxpow];
    for (int x = 0; x < 3; ++x) {
        const double mu1 = ((q + t) * p * ab.P[x] + t * q * cd.P[x]) * rdet;
        const double mu2 = (t * p * ab.P[x] + (p + t) * q * cd.P[x]) * rdet;
        const double xa = mu1 - L.A[x], xb = mu1 - L.B[x], xc = mu2 - L.C[x], xd = mu2 - L.D[x];
        double* gx = g[x];
        std::fill(gx, gx + nbra * nket, 0.0);
        for (int k1 = 0; k1 < n; ++k1) {
            for (int k2 = 0; k2 < n; ++k2) {
                const double d1 = o1[k1 * n + k2], d2 = o2[k2];
                const double wk = L.ww[k1 * n + k2];
                powers(xa + d1, L.la, pa);
                powers(xb + d1, L.lb, pb);
                powers(xc + d2, L.lc, pc);
                powers(xd + d2, L.ld, pd);
                for (int i = 0; i <= L.la; ++i)
                    for (int j = 0; j <= L.lb; ++j) u[i * nb + j] = wk * pa[i] * pb[j];
                for (int k = 0; k <= L.lc; ++k)
                    for (int l = 0; l <= L.ld; ++l) v[k * nd + l] = pc[k] * pd[l];
                for (int ij = 0; ij < nbra; ++ij)
                    for (int kl = 0; kl < nket; ++kl) gx[ij * nket + kl] += u[ij] * v[kl];
            }
        }
    }
    return pref;
}

/// an (ss|ss) primitive quartet summed over the kernel terms
///
/// With one Gauss-Hermite node per axis (y = 0, w = sqrt(pi)) the quadrature of
/// coulomb_1d reduces to the closed form
///   K_ab K_cd pi^3 sum_m w_m exp(-omega_m |P-Q|^2) det_m^{-3/2}.
double coulomb_ssss(const PrimitivePair& ab, const PrimitivePair& cd, const GaussianKernel& kernel) {
    const double pq = ab.p * cd.p, s = ab.p + cd.p;
    double r2 = 0.0;
    for (int x = 0; x < 3; ++x) r2 += (ab.P[x] - cd.P[x]) * (ab.P[x] - cd.P[x]);
    double sum = 0.0;
    for (std::size_t m = 0; m < kernel.size(); ++m) {
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
void primitive_quartet(const QuartetLayout& L, const PrimitivePair& ab, const PrimitivePair& cd,
                       const GaussianKernel& kernel, std::vector<double> (&g)[3], double* acc) {
    if (L.la + L.lb + L.lc + L.ld == 0) {
        acc[0] = coulomb_ssss(ab, cd, kernel);
        return;
    }
    const int nblock = L.nblock();
    std::fill(acc, acc + nblock, 0.0);
    double* const gp[3] = {g[0].data(), g[1].data(), g[2].data()};
    for (std::size_t m = 0; m < kernel.size(); ++m) {
        const double f = kernel.w[m] * coulomb_xyz(L, ab, cd, kernel.t[m], gp);
        for (int n = 0; n < nblock; ++n) {
            const std::array<int,3>& i = L.index[n];
            acc[n] += f * gp[0][i[0]] * gp[1][i[1]] * gp[2][i[2]];
        }
    }
}

/// all shell quartets of one group quartet (G1 G2|G3 G4), written into G
///
/// Each primitive quartet is computed once and contracted into every shell
/// quartet of the groups. The shell quartets are restricted so that each one is
/// computed exactly once over all canonical group quartets, as in the canonical
/// loop over shells: s1 >= s2 within one group, s3 >= s4 within one group, and
/// bra >= ket when the bra and ket groups are the same.
void group_quartet(const std::vector<Shell>& shells, const std::vector<int>& G1, const std::vector<int>& G2,
                   const std::vector<int>& G3, const std::vector<int>& G4, const bool same12, const bool same34,
                   const bool samebraket, const std::vector<PrimitivePair>& bra, const std::vector<PrimitivePair>& ket,
                   const GaussianKernel& kernel, const GaussHermiteRule& gh, Tensor<double>& G) {
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

    const QuartetLayout L(shells[G1[0]], shells[G2[0]], shells[G3[0]], shells[G4[0]], gh);
    const int nblock = L.nblock();
    std::vector<double> g[3], acc(nblock), blocks(members.size() * nblock, 0.0);
    for (auto& v : g) v.resize(L.ntab);
    for (const PrimitivePair& pab : bra) {
        for (const PrimitivePair& pcd : ket) {
            primitive_quartet(L, pab, pcd, kernel, g, acc.data());
            for (std::size_t q = 0; q < members.size(); ++q) {
                const std::array<int,4>& s = members[q];
                const double coef = shells[s[0]].coeff[pab.ia] * shells[s[1]].coeff[pab.ib]
                                  * shells[s[2]].coeff[pcd.ia] * shells[s[3]].coeff[pcd.ib];
                double* block = blocks.data() + q * nblock;
                for (int n = 0; n < nblock; ++n) block[n] += coef * acc[n];
            }
        }
    }
    for (std::size_t q = 0; q < members.size(); ++q) {
        const std::array<int,4>& s = members[q];
        scatter_quartet(G, shells[s[0]], shells[s[1]], shells[s[2]], shells[s[3]], blocks.data() + q * nblock);
    }
}

/// what the tasks of SeparatedGaussianIntegrals::eri share: shell groups, their primitive pairs, the kernel
struct GroupQuartetData {
    const std::vector<Shell>& shells;
    std::vector<std::vector<int>> groups;
    std::vector<std::vector<PrimitivePair>> pairs;      ///< per group pair, numbered by pair_index
    const GaussianKernel& kernel;
    const GaussHermiteRule& gh;
};

/// a task on the thread pool: all group quartets with bra group pair (a b)
///
/// The group quartets of different tasks write disjoint elements of G, so the
/// tasks need no locks, and every element is computed in the order of the serial
/// loop.
class BraPairTask : public TaskInterface {
public:
    BraPairTask(const GroupQuartetData& data, const std::size_t a, const std::size_t b, Tensor<double>& G)
        : data_(data), a_(a), b_(b), G_(G) {}

    using TaskInterface::run;
    void run(World&) override {
        const std::size_t ab = pair_index(a_, b_);
        for (std::size_t c = 0; c <= a_; ++c) {
            for (std::size_t d = 0; d <= c; ++d) {
                const std::size_t cd = pair_index(c, d);
                if (cd > ab) continue;     // (ab|cd) = (cd|ab)
                group_quartet(data_.shells, data_.groups[a_], data_.groups[b_], data_.groups[c], data_.groups[d],
                              a_ == b_, c == d, ab == cd, data_.pairs[ab], data_.pairs[cd], data_.kernel, data_.gh,
                              G_);
            }
        }
    }

private:
    const GroupQuartetData& data_;
    const std::size_t a_, b_;
    Tensor<double>& G_;
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


SeparatedGaussianIntegrals::SeparatedGaussianIntegrals(const std::vector<Shell>& shells,
                                                       const GaussianKernel& coulomb)
    : shells_(shells), coulomb_(coulomb), gh_(16) {
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


Tensor<double> SeparatedGaussianIntegrals::eri(World& world) const {
    Tensor<double> G(nbf_, nbf_, nbf_, nbf_);
    int lmax = 0;
    for (const Shell& s : shells_) lmax = std::max(lmax, s.l);
    // checked here, so that QuartetLayout cannot throw inside a task
    MADNESS_CHECK_THROW(2 * lmax + 1 <= std::min(maxnode, gh_.nmax()), "eri: angular momentum too high");

    GroupQuartetData data{shells_, shell_groups(shells_), {}, coulomb_, gh_};
    const std::size_t ng = data.groups.size();
    data.pairs.resize(ng * (ng + 1) / 2);
    for (std::size_t a = 0; a < ng; ++a)
        for (std::size_t b = 0; b <= a; ++b)
            data.pairs[pair_index(a, b)] = primitive_pairs(shells_[data.groups[a][0]], shells_[data.groups[b][0]]);

    // the last bra pairs have the most ket pairs: submit them first, so that no large task comes last
    for (std::size_t a = ng; a-- > 0;)
        for (std::size_t b = a + 1; b-- > 0;) world.taskq.add(new BraPairTask(data, a, b, G));
    world.taskq.fence();
    return G;
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
