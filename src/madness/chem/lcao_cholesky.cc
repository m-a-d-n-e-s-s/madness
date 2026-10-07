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

/// \file lcao_cholesky.cc
/// \brief pivoted Cholesky decomposition of the two-electron integrals over Gaussian basis functions

#include <madness/chem/lcao_cholesky.h>
#include <madness/tensor/cblas.h>
#include <madness/tensor/tensor_lapack.h>
#include <madness/world/MADworld.h>

#include <algorithm>
#include <cmath>

namespace madness {
namespace lcao {

namespace {

/// rows per block of the updates: the unit of the task split, independent of the number of threads
constexpr std::size_t block_rows = 4096;

/// a task on the thread pool: f(b) for one block b
template <typename F>
class BlockTask : public TaskInterface {
public:
    BlockTask(const F& f, const std::size_t b) : f_(f), b_(b) {}

    using TaskInterface::run;
    void run(World&) override { f_(b_); }

private:
    const F& f_;
    const std::size_t b_;
};

/// f(b) for the blocks b = 0..nblock-1, as tasks on the thread pool of this process; returns when all are done
template <typename F>
void for_each_block(World& world, const std::size_t nblock, const F& f) {
    for (std::size_t b = 0; b < nblock; ++b) world.taskq.add(new BlockTask<F>(f, b));
    world.taskq.fence();
}

} // namespace


CholeskyERIDecomposition::CholeskyERIDecomposition(World& world, const SeparatedGaussianIntegrals& ints,
                                                   const double tol, const double span)
    : tol_(tol), span_(span), pairs_(ints.function_pairs()) {
    MADNESS_CHECK_THROW(tol > 0.0 and span > 0.0 and span < 1.0,
                        "CholeskyERIDecomposition: need tol > 0 and 0 < span < 1");
    const FunctionPairs& fp = pairs_;
    stats_.npairs = fp.size();

    double t0 = wall_time();
    const Tensor<double> D0 = ints.pair_diagonal(world, fp);
    stats_.t_diagonal = wall_time() - t0;

    // pair screening against the largest diagonal
    double d0max = 0.0;
    for (long r = 0; r < D0.size(); ++r) d0max = std::max(d0max, D0(r));
    for (std::size_t r = 0; r < fp.size(); ++r)
        if (D0(long(r)) >= tol * tol / d0max) rows_.push_back(r);
    const std::size_t nk = rows_.size();
    stats_.nkept = nk;
    if (nk == 0) return;
    std::vector<double> D(nk);
    for (std::size_t i = 0; i < nk; ++i) D[i] = D0(long(rows_[i]));

    // eri_columns computes whole group pairs: the bra group pairs are those with kept rows,
    // and kept row i is element wpos[i] of a computed column
    std::vector<std::size_t> bra, wpos(nk);
    std::vector<long> kept(fp.size(), -1);     // per row, its index among the kept rows
    for (std::size_t i = 0, start = 0; i < nk; ++i) {
        const std::size_t g = fp.group_pair[rows_[i]];
        if (bra.empty() or bra.back() != g) {
            if (not bra.empty()) start += fp.count[bra.back()];
            bra.push_back(g);
        }
        wpos[i] = start + (rows_[i] - fp.first[g]);
        kept[rows_[i]] = long(i);
    }
    const std::size_t nblock = (nk + block_rows - 1) / block_rows;

    std::vector<double> V, Lc, f;
    while (true) {
        // the largest residual diagonal; ties go to the lowest row
        std::size_t p = 0;
        for (std::size_t i = 1; i < nk; ++i)
            if (D[i] > D[p]) p = i;
        const double dmax = D[p];
        if (dmax <= tol_) break;

        // the candidates: kept rows of the pivot's group pair with D > max(tol, span dmax), ascending
        const std::size_t g = fp.group_pair[rows_[p]];
        const double dmin = std::max(tol_, span_ * dmax);
        std::vector<std::size_t> cand;
        for (std::size_t r = fp.first[g]; r < fp.first[g] + fp.count[g]; ++r)
            if (kept[r] >= 0 and D[kept[r]] > dmin) cand.push_back(std::size_t(kept[r]));
        const std::size_t nc = cand.size();

        // their integral columns on the kept rows: V[q nk + i] = (row i|candidate q)
        t0 = wall_time();
        const Tensor<double> W = ints.eri_columns(world, fp, g, bra);
        const std::size_t ldw = W.dim(1);
        V.assign(nc * nk, 0.0);
        for (std::size_t q = 0; q < nc; ++q) {
            const double* w = W.ptr() + (rows_[cand[q]] - fp.first[g]) * ldw;
            double* v = V.data() + q * nk;
            for (std::size_t i = 0; i < nk; ++i) v[i] = w[wpos[i]];
        }
        stats_.ncolumns += fp.count[g];
        ++stats_.nbatches;
        stats_.t_integrals += wall_time() - t0;

        // minus the previous vectors, V -= Lc^T L with Lc[q nvec + k] = L_k(candidate q): one dgemm per block
        t0 = wall_time();
        if (nvec_ > 0) {
            Lc.resize(nc * nvec_);
            for (std::size_t q = 0; q < nc; ++q)
                for (long k = 0; k < nvec_; ++k) Lc[q * nvec_ + k] = L_[k * nk + cand[q]];
            for_each_block(world, nblock, [&](const std::size_t b) {
                const std::size_t i0 = b * block_rows, n = std::min(nk, i0 + block_rows) - i0;
                cblas::gemm(cblas::NoTrans, cblas::NoTrans, long(n), long(nc), nvec_, -1.0, L_.data() + i0, long(nk),
                            Lc.data(), nvec_, 1.0, V.data() + i0, long(nk));
            });
        }

        // the candidates in descending residual diagonal, re-checked after each new vector
        std::vector<char> done(nc, 0);
        f.resize(nc);
        while (true) {
            std::size_t qp = nc;
            double dq = dmin;
            for (std::size_t q = 0; q < nc; ++q)
                if (not done[q] and D[cand[q]] > dq) {
                    dq = D[cand[q]];
                    qp = q;
                }
            if (qp == nc) break;
            done[qp] = 1;

            // L_new = V[qp]/sqrt(D), f[q] its value on the remaining candidates
            const double s = 1.0 / std::sqrt(dq);
            const double* v = V.data() + qp * nk;
            for (std::size_t q = 0; q < nc; ++q) f[q] = done[q] ? 0.0 : v[cand[q]] * s;
            L_.resize((nvec_ + 1) * nk);
            double* l = L_.data() + nvec_ * nk;
            for_each_block(world, nblock, [&](const std::size_t b) {
                const std::size_t i0 = b * block_rows, i1 = std::min(nk, i0 + block_rows);
                for (std::size_t i = i0; i < i1; ++i) {
                    l[i] = v[i] * s;
                    D[i] = std::max(D[i] - l[i] * l[i], 0.0);
                }
                for (std::size_t q = 0; q < nc; ++q) {
                    if (done[q]) continue;
                    double* vq = V.data() + q * nk;
                    for (std::size_t i = i0; i < i1; ++i) vq[i] -= f[q] * l[i];
                }
            });
            D[cand[qp]] = 0.0;
            ++nvec_;
        }
        stats_.t_updates += wall_time() - t0;
    }
}


CholeskyERI::CholeskyERI(World& world, std::shared_ptr<const CholeskyERIDecomposition> chol, const long nbf)
    : world_(world), chol_(std::move(chol)), nbf_(nbf) {
    const FunctionPairs& fp = chol_->pairs();
    MADNESS_CHECK_THROW(fp.size() == std::size_t(nbf) * (nbf + 1) / 2,
                        "CholeskyERI: the decomposition belongs to another basis");
    mu_.reserve(chol_->nkept());
    nu_.reserve(chol_->nkept());
    for (const std::size_t r : chol_->rows()) {
        mu_.push_back(fp.functions[r][0]);
        nu_.push_back(fp.functions[r][1]);
    }
}


namespace {

/// one term of a density matrix written as P = sum_t sign_t Y_t Y_t^T, Y_t n x r, column-major
struct DensityFactor {
    std::vector<double> Y;
    long r = 0;
    double sign = 1.0;
    int spin = 0;
};

/// P = Y+ Y+^T - Y- Y-^T from the eigenpairs of the symmetric P; eigenvalues within 1e-13 max|e| of 0 are dropped
void factor_density(const Tensor<double>& P, const int spin, std::vector<DensityFactor>& factors) {
    const long n = P.dim(0);
    Tensor<double> U, e;
    syev(P, U, e);
    const double cut = 1.e-13 * e.absmax();
    for (const double sign : {1.0, -1.0}) {
        DensityFactor f;
        f.sign = sign;
        f.spin = spin;
        for (long j = 0; j < e.size(); ++j) {
            if (sign * e(j) <= cut) continue;
            const double s = std::sqrt(sign * e(j));
            for (long mu = 0; mu < n; ++mu) f.Y.push_back(U(mu, j) * s);
            ++f.r;
        }
        if (f.r > 0) factors.push_back(std::move(f));
    }
}

} // namespace


void CholeskyERI::jk(const Tensor<double>& Pa, const Tensor<double>& Pb, Tensor<double>& J, Tensor<double>& Ka,
                     Tensor<double>& Kb) const {
    const long n = nbf_;
    const bool open = Pb.size() > 0;
    MADNESS_CHECK_THROW(Pa.ndim() == 2 and Pa.dim(0) == n and Pa.dim(1) == n and
                        (not open or (Pb.dim(0) == n and Pb.dim(1) == n)),
                        "CholeskyERI: integrals and density matrix do not match");
    const long m = chol_->nvec();
    const std::size_t nk = chol_->nkept();
    const double* L = chol_->vectors();
    const std::size_t nblock = (nk + block_rows - 1) / block_rows;

    // J: gamma = sum over the kept rows of L P~, P~ the total density with both orders of
    // mu != nu; then J~ = L^T gamma. Partial sums of gamma per block of rows, added in block order.
    const Tensor<double> pt = open ? Pa + Pb : 2.0 * Pa;
    std::vector<double> ptk(nk), jt(nk, 0.0), gamma(std::max(m, 1L), 0.0);
    for (std::size_t i = 0; i < nk; ++i) ptk[i] = (mu_[i] == nu_[i] ? 1.0 : 2.0) * pt(mu_[i], nu_[i]);
    if (m > 0) {
        std::vector<double> gpart(nblock * m, 0.0);
        for_each_block(world_, nblock, [&](const std::size_t b) {
            const std::size_t i0 = b * block_rows, nb = std::min(nk, i0 + block_rows) - i0;
            cblas::gemv(cblas::Trans, long(nb), m, 1.0, L + i0, long(nk), ptk.data() + i0, 1, 0.0,
                        gpart.data() + b * m, 1);
        });
        for (std::size_t b = 0; b < nblock; ++b)
            for (long k = 0; k < m; ++k) gamma[k] += gpart[b * m + k];
        for_each_block(world_, nblock, [&](const std::size_t b) {
            const std::size_t i0 = b * block_rows, nb = std::min(nk, i0 + block_rows) - i0;
            cblas::gemv(cblas::NoTrans, long(nb), m, 1.0, L + i0, long(nk), gamma.data(), 1, 0.0, jt.data() + i0, 1);
        });
    }
    J = Tensor<double>(n, n);
    for (std::size_t i = 0; i < nk; ++i) J(mu_[i], nu_[i]) = J(nu_[i], mu_[i]) = jt[i];

    // K per spin from the factors of its density: chunk c sums X_k X_k^T over its vectors, several
    // vectors per dgemm; the chunks do not depend on the number of threads and are added in order
    std::vector<DensityFactor> factors;
    factor_density(Pa, 0, factors);
    if (open) factor_density(Pb, 1, factors);
    const int nspin = open ? 2 : 1;
    long rmax = 1;
    for (const DensityFactor& f : factors) rmax = std::max(rmax, f.r);
    const long nstack = std::max(1L, 512 / rmax);
    const long nchunk = std::min(16L, m);
    std::vector<std::vector<double>> kpart(nchunk * nspin);
    for_each_block(world_, std::size_t(nchunk), [&](const std::size_t c) {
        const long k0 = m * long(c) / nchunk, k1 = m * long(c + 1) / nchunk;
        for (int s = 0; s < nspin; ++s) kpart[c * nspin + s].assign(n * n, 0.0);
        std::vector<double> Lk(n * n, 0.0);     // the vector as a matrix; screened pairs stay 0
        std::vector<std::vector<double>> X(factors.size());
        for (std::size_t t = 0; t < factors.size(); ++t) X[t].resize(n * factors[t].r * nstack);
        for (long k = k0; k < k1; k += nstack) {
            const long kn = std::min(nstack, k1 - k);
            for (long kk = 0; kk < kn; ++kk) {
                const double* l = L + (k + kk) * nk;
                for (std::size_t i = 0; i < nk; ++i) Lk[mu_[i] * n + nu_[i]] = Lk[nu_[i] * n + mu_[i]] = l[i];
                for (std::size_t t = 0; t < factors.size(); ++t)
                    cblas::gemm(cblas::NoTrans, cblas::NoTrans, n, factors[t].r, n, 1.0, Lk.data(), n,
                                factors[t].Y.data(), n, 0.0, X[t].data() + kk * n * factors[t].r, n);
            }
            for (std::size_t t = 0; t < factors.size(); ++t)
                cblas::gemm(cblas::NoTrans, cblas::Trans, n, n, kn * factors[t].r, factors[t].sign, X[t].data(), n,
                            X[t].data(), n, 1.0, kpart[c * nspin + factors[t].spin].data(), n);
        }
    });
    const auto sum_chunks = [&](const int s) {
        Tensor<double> K(n, n);
        double* k = K.ptr();
        for (long c = 0; c < nchunk; ++c) {
            const double* kp = kpart[c * nspin + s].data();
            for (long i = 0; i < n * n; ++i) k[i] += kp[i];
        }
        return Tensor<double>(0.5 * (K + transpose(K)));
    };
    Ka = sum_chunks(0);
    Kb = open ? sum_chunks(1) : Tensor<double>();
}

} // namespace lcao
} // namespace madness
