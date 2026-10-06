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


Tensor<double> SeparatedGaussianIntegrals::eri() const {
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

                    // scatter with the 8-fold permutational symmetry
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
            }
        }
    }
    return G;
}

} // namespace lcao
} // namespace madness
