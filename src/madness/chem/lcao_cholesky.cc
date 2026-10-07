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

} // namespace lcao
} // namespace madness
