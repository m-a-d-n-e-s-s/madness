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

/// \file lcao_stability.cc
/// \brief the real orbital Hessian of an LCAO SCF solution and its lowest eigenvalues

#include <madness/chem/lcao_stability.h>
#include <madness/tensor/tensor_lapack.h>

#include <algorithm>
#include <cmath>
#include <numeric>

namespace madness {
namespace lcao {

std::string to_string(const StabilityBlock b) {
    switch (b) {
        case StabilityBlock::rhf_rhf: return "RHF->RHF";
        case StabilityBlock::rhf_uhf: return "RHF->UHF";
        case StabilityBlock::uhf_uhf: return "UHF->UHF";
    }
    return "?";
}


OrbitalHessian::OrbitalHessian(const LCAOSCF& scf, const LCAOSCF::Result& state, const int nalpha, const int nbeta,
                               const StabilityBlock block)
    : scf_(scf), block_(block) {
    const bool rhf = block != StabilityBlock::uhf_uhf;
    const long nmo = state.Ca.dim(1);
    na_ = nalpha;
    nva_ = nmo - nalpha;
    nb_ = rhf ? 0 : nbeta;
    nvb_ = rhf ? 0 : state.Cb.dim(1) - nbeta;
    MADNESS_CHECK_THROW(na_ > 0 and nva_ > 0, "OrbitalHessian: no occupied or no virtual orbitals");
    MADNESS_CHECK_THROW(not rhf or nalpha == nbeta, "OrbitalHessian: the RHF blocks need a closed shell");
    coa_ = copy(state.Ca(_, Slice(0, na_ - 1)));
    cva_ = copy(state.Ca(_, Slice(na_, nmo - 1)));
    if (nb_ > 0) {
        cob_ = copy(state.Cb(_, Slice(0, nb_ - 1)));
        cvb_ = copy(state.Cb(_, Slice(nb_, state.Cb.dim(1) - 1)));
    }
    diag_ = Tensor<double>(dim());
    long k = 0;
    for (long a = 0; a < nva_; ++a)
        for (long i = 0; i < na_; ++i) diag_(k++) = state.epsa(na_ + a) - state.epsa(i);
    for (long a = 0; a < nvb_; ++a)
        for (long i = 0; i < nb_; ++i) diag_(k++) = state.epsb(nb_ + a) - state.epsb(i);
}


void OrbitalHessian::split(const Tensor<double>& x, Tensor<double>& xa, Tensor<double>& xb) const {
    MADNESS_CHECK_THROW(x.size() == dim(), "OrbitalHessian: the vector does not match the Hessian");
    xa = Tensor<double>(nva_, na_);
    std::copy(x.ptr(), x.ptr() + nva_ * na_, xa.ptr());
    xb = Tensor<double>();
    if (nb_ > 0) {
        xb = Tensor<double>(nvb_, nb_);
        std::copy(x.ptr() + nva_ * na_, x.ptr() + dim(), xb.ptr());
    }
}


Tensor<double> OrbitalHessian::product(const Tensor<double>& x) const {
    Tensor<double> xa, xb;
    split(copy(x), xa, xb);
    const long n = coa_.dim(0);
    Tensor<double> J, Ka, Kb;
    const Tensor<double> Xa = inner(cva_, xa);
    if (block_ == StabilityBlock::uhf_uhf) {
        // a spin without electrons contributes nothing: a zero factor, not the closed-shell default of an empty one
        const Tensor<double> Xb = (nb_ > 0) ? inner(cvb_, xb) : Tensor<double>(n, 1L);
        const Tensor<double> Yb = (nb_ > 0) ? cob_ : Tensor<double>(n, 1L);
        scf_.two_electron().jk_transition(Xa, coa_, Xb, Yb, J, Ka, Kb);
    } else {
        scf_.two_electron().jk_transition(Xa, coa_, Tensor<double>(), Tensor<double>(), J, Ka, Kb);
    }
    // per spin: the orbital energy differences, then C_v^T F1 C_o with the response Fock matrix F1
    const auto spin = [](const Tensor<double>& F1, const Tensor<double>& cv, const Tensor<double>& co,
                         const Tensor<double>& xs, const Tensor<double>& d, const long offset) {
        Tensor<double> s = inner(transpose(cv), inner(F1, co));
        for (long a = 0, k = offset; a < s.dim(0); ++a)
            for (long i = 0; i < s.dim(1); ++i, ++k) s(a, i) += d(k) * xs(a, i);
        return s;
    };
    Tensor<double> F1a;
    switch (block_) {
        case StabilityBlock::rhf_rhf: F1a = J - Ka; break;     // J = J[2D]: the singlet response
        case StabilityBlock::rhf_uhf: F1a = -1.0 * Ka; break;  // Db = -Da: no Coulomb response
        case StabilityBlock::uhf_uhf: F1a = J - Ka; break;
    }
    const Tensor<double> sa = spin(F1a, cva_, coa_, xa, diag_, 0);
    Tensor<double> sigma(dim());
    std::copy(sa.ptr(), sa.ptr() + sa.size(), sigma.ptr());
    if (nb_ > 0) {
        const Tensor<double> sb = spin(J - Kb, cvb_, cob_, xb, diag_, nva_ * na_);
        std::copy(sb.ptr(), sb.ptr() + sb.size(), sigma.ptr() + sa.size());
    }
    return sigma;
}


StabilityRoots lowest_roots(const OrbitalHessian& H, const int nroots_in, const double tol, const int maxiter) {
    StabilityRoots out;
    out.block = H.block();
    const long n = H.dim();
    const long nroots = std::min<long>(nroots_in, n);
    const Tensor<double>& d = H.diagonal();
    const long maxsub = std::min<long>(n, std::max<long>(8 * nroots, 40));

    // start: unit vectors at the smallest diagonal elements (ties to the lower index), each with a small fixed
    // pseudo-random admixture. H keeps the symmetry of a vector, so unit vectors alone would never reach roots
    // of a symmetry their positions miss (ozone's third RHF->RHF root, water's third RHF->UHF root against PySCF).
    std::vector<long> order(n);
    std::iota(order.begin(), order.end(), 0L);
    std::stable_sort(order.begin(), order.end(), [&d](const long a, const long b) { return d(a) < d(b); });
    std::vector<Tensor<double>> V, W;
    for (long k = 0; k < std::min<long>(n, nroots + 4); ++k) {
        Tensor<double> v(n);
        for (long i = 0; i < n; ++i) {
            const double r = std::sin(12.9898 * double(i + 1) + 78.233 * double(k + 1)) * 43758.5453;
            v(i) = 1.e-2 * (r - std::floor(r) - 0.5);
        }
        v(order[k]) += 1.0;
        for (int pass = 0; pass < 2; ++pass)
            for (const Tensor<double>& u : V) v.gaxpy(1.0, u, -u.trace(v));
        const double norm = v.normf();
        if (norm < 1.e-8) continue;
        v.scale(1.0 / norm);
        V.push_back(v);
    }
    for (const Tensor<double>& v : V) {
        W.push_back(H.product(v));
        ++out.products;
    }

    Tensor<double> theta, s;
    for (int iter = 0; iter < maxiter; ++iter) {
        const long m = long(V.size());
        Tensor<double> G(m, m);
        for (long i = 0; i < m; ++i)
            for (long j = 0; j <= i; ++j) G(i, j) = G(j, i) = 0.5 * (V[i].trace(W[j]) + V[j].trace(W[i]));
        syev(G, s, theta);
        // Ritz vectors and residuals of the lowest roots
        std::vector<Tensor<double>> ritz, res;
        bool done = true;
        for (long k = 0; k < nroots; ++k) {
            Tensor<double> y(n), r(n);
            for (long i = 0; i < m; ++i) {
                y.gaxpy(1.0, V[i], s(i, k));
                r.gaxpy(1.0, W[i], s(i, k));
            }
            r.gaxpy(1.0, y, -theta(k));
            done = done and r.normf() < tol;
            ritz.push_back(y);
            res.push_back(r);
        }
        if (done) {
            out.converged = true;
            break;
        }
        // collapse a large subspace to the current Ritz vectors
        if (m + nroots > maxsub) {
            std::vector<Tensor<double>> V2, W2;
            for (long k = 0; k < std::min<long>(m, nroots + 2); ++k) {
                Tensor<double> y(n), w(n);
                for (long i = 0; i < m; ++i) {
                    y.gaxpy(1.0, V[i], s(i, k));
                    w.gaxpy(1.0, W[i], s(i, k));
                }
                V2.push_back(y);
                W2.push_back(w);
            }
            V = V2;
            W = W2;
        }
        // preconditioned residuals, orthonormalized against the subspace (twice)
        std::vector<Tensor<double>> fresh;
        for (long k = 0; k < nroots; ++k) {
            if (res[k].normf() < tol) continue;
            Tensor<double> t(n);
            for (long i = 0; i < n; ++i) {
                double den = theta(k) - d(i);
                if (std::abs(den) < 1.e-3) den = (den < 0.0) ? -1.e-3 : 1.e-3;
                t(i) = res[k](i) / den;
            }
            for (int pass = 0; pass < 2; ++pass) {
                for (const Tensor<double>& v : V) t.gaxpy(1.0, v, -v.trace(t));
                for (const Tensor<double>& v : fresh) t.gaxpy(1.0, v, -v.trace(t));
            }
            const double norm = t.normf();
            if (norm < 1.e-8) continue;
            t.scale(1.0 / norm);
            fresh.push_back(t);
        }
        if (fresh.empty()) break;
        for (const Tensor<double>& t : fresh) {
            V.push_back(t);
            W.push_back(H.product(t));
            ++out.products;
        }
    }
    // the final Ritz pairs
    const long m = long(V.size());
    Tensor<double> G(m, m);
    for (long i = 0; i < m; ++i)
        for (long j = 0; j <= i; ++j) G(i, j) = G(j, i) = 0.5 * (V[i].trace(W[j]) + V[j].trace(W[i]));
    syev(G, s, theta);
    for (long k = 0; k < nroots; ++k) {
        Tensor<double> y(n);
        for (long i = 0; i < m; ++i) y.gaxpy(1.0, V[i], s(i, k));
        Tensor<double> xa, xb;
        H.split(y, xa, xb);
        out.eigenvalues.push_back(theta(k));
        out.xa.push_back(xa);
        out.xb.push_back(xb);
    }
    return out;
}


Tensor<double> rotate_occupied(const Tensor<double>& C, const long nocc, const Tensor<double>& x, const double t) {
    const long nmo = C.dim(1);
    const Tensor<double> Co = copy(C(_, Slice(0, nocc - 1)));
    const Tensor<double> Cv = copy(C(_, Slice(nocc, nmo - 1)));
    Tensor<double> U, sv, VT;
    svd(x, U, sv, VT);                      // x (nv x nocc) = U diag(sv) VT
    const long k = sv.size();
    const Tensor<double> V = transpose(VT); // nocc x k
    Tensor<double> Vc(V.dim(0), k), Us(U.dim(0), k);
    for (long j = 0; j < k; ++j) {
        const double c = std::cos(t * sv(j)), s = std::sin(t * sv(j));
        for (long i = 0; i < V.dim(0); ++i) Vc(i, j) = V(i, j) * (c - 1.0);
        for (long i = 0; i < U.dim(0); ++i) Us(i, j) = U(i, j) * s;
    }
    // C_o (1 + V (cos - 1) V^T) + C_v U sin V^T
    Tensor<double> R = inner(Vc, transpose(V));
    for (long i = 0; i < R.dim(0); ++i) R(i, i) += 1.0;
    return Tensor<double>(inner(Co, R) + inner(Cv, inner(Us, transpose(V))));
}

} // namespace lcao
} // namespace madness
