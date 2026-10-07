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
#include <cstring>
#include <limits>

namespace madness {
namespace lcao {

namespace {

/// kept rows per block of rows, the unit of the rank split: consecutive group pairs, cut where they reach this
constexpr std::size_t block_rows = 4096;

/// rows per task of the updates: each block is cut into pieces of at most this many rows,
/// independent of the numbers of threads and ranks
constexpr std::size_t task_rows = 512;

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
                                                   const double tol, const double span, const bool collective)
    : tol_(tol), span_(span), distributed_(collective and world.size() > 1),
      nproc_(distributed_ ? world.size() : 1), rank_(distributed_ ? world.rank() : 0), pairs_(ints.function_pairs()) {
    MADNESS_CHECK_THROW(tol > 0.0 and span > 0.0 and span < 1.0,
                        "CholeskyERIDecomposition: need tol > 0 and 0 < span < 1");
    const FunctionPairs& fp = pairs_;
    stats_.npairs = fp.size();

    // the diagonal on every rank (cheap), so the pair screening and the blocks are the same everywhere
    double t0 = wall_time();
    const Tensor<double> D0 = ints.pair_diagonal(world, fp);
    stats_.t_diagonal = wall_time() - t0;
    double d0max = 0.0;
    for (long r = 0; r < D0.size(); ++r) d0max = std::max(d0max, D0(r));
    for (std::size_t r = 0; r < fp.size(); ++r)
        if (D0(long(r)) >= tol * tol / d0max) rows_.push_back(r);
    const std::size_t nk = rows_.size();
    stats_.nkept = nk;
    nchunk_ = std::max(32L, long(nproc_));

    // blocks of rows: consecutive group pairs, cut where a block reaches block_rows kept rows;
    // block b holds the kept rows [bfirst[b], bfirst[b+1]); a rank holds a contiguous range of blocks
    std::vector<std::size_t> bfirst{0};
    for (std::size_t i = 0; i < nk; ++i) {
        const bool newpair = i > 0 and fp.group_pair[rows_[i]] != fp.group_pair[rows_[i - 1]];
        if (newpair and i - bfirst.back() >= block_rows) bfirst.push_back(i);
    }
    bfirst.push_back(nk);
    const std::size_t nblock = bfirst.size() - 1;
    std::vector<int> bowner(nblock);
    for (std::size_t b = 0; b < nblock; ++b) bowner[b] = int(std::min<std::size_t>(nproc_ - 1, bfirst[b] * nproc_ / std::max<std::size_t>(nk, 1)));
    std::size_t b0 = 0, b1 = 0;     // this rank's blocks
    while (b0 < nblock and bowner[b0] < rank_) ++b0;
    for (b1 = b0; b1 < nblock and bowner[b1] == rank_; ++b1) {}
    const std::size_t r0 = bfirst[b0], r1 = bfirst[b1], nloc = r1 - r0;   // this rank's kept rows

    // per row, its index among the kept rows; per kept row, its block
    std::vector<long> kept(fp.size(), -1);
    for (std::size_t i = 0; i < nk; ++i) kept[rows_[i]] = long(i);
    std::vector<std::size_t> kblock(nk);
    for (std::size_t b = 0; b < nblock; ++b)
        for (std::size_t i = bfirst[b]; i < bfirst[b + 1]; ++i) kblock[i] = b;

    // eri_columns computes whole group pairs: the bra group pairs are those of this rank's rows,
    // and local row i is element wpos[i] of a computed column
    std::vector<std::size_t> bra, wpos(nloc);
    for (std::size_t i = 0, start = 0; i < nloc; ++i) {
        const std::size_t g = fp.group_pair[rows_[r0 + i]];
        if (bra.empty() or bra.back() != g) {
            if (not bra.empty()) start += fp.count[bra.back()];
            bra.push_back(g);
        }
        wpos[i] = start + (rows_[r0 + i] - fp.first[g]);
    }
    std::vector<double> D(nloc);
    for (std::size_t i = 0; i < nloc; ++i) D[i] = D0(long(rows_[r0 + i]));

    // the updates run as tasks over pieces of this rank's blocks, [tfirst[t], tfirst[t+1]) in local rows
    std::vector<std::size_t> tfirst;
    for (std::size_t b = b0; b < b1; ++b)
        for (std::size_t i = bfirst[b]; i < bfirst[b + 1]; i += task_rows) tfirst.push_back(i - r0);
    tfirst.push_back(nloc);
    const auto local_blocks = [&](const auto& f) {
        for_each_block(world, tfirst.size() - 1, [&](const std::size_t t) { f(tfirst[t], tfirst[t + 1]); });
    };
    std::vector<double> V, Lc, f, vcc, dc;
    std::vector<long> cand;
    while (nk > 0) {
        // the largest residual diagonal over all ranks; ties go to the lowest kept row
        long p = -1;
        for (std::size_t i = 0; i < nloc; ++i)
            if (p < 0 or D[i] > D[p]) p = long(i);
        double dmax = (p >= 0) ? D[p] : -1.0;
        const double dmax_local = dmax;
        if (distributed_) world.gop.max(dmax);
        if (dmax <= tol_) break;
        long pivot = (p >= 0 and dmax_local == dmax) ? long(r0) + p : std::numeric_limits<long>::max();
        if (distributed_) world.gop.min(pivot);

        // the owner of the pivot's group pair finds the candidates, kept rows with D > max(tol, span dmax)
        const std::size_t g = fp.group_pair[rows_[pivot]];
        const int owner = bowner[kblock[pivot]];
        const double dmin = std::max(tol_, span_ * dmax);
        cand.clear();
        dc.clear();
        if (owner == rank_) {
            for (std::size_t r = fp.first[g]; r < fp.first[g] + fp.count[g]; ++r)
                if (kept[r] >= 0 and D[kept[r] - r0] > dmin) {
                    cand.push_back(kept[r]);
                    dc.push_back(D[kept[r] - r0]);
                }
        }
        if (distributed_) {
            world.gop.broadcast_serializable(cand, owner);
            world.gop.broadcast_serializable(dc, owner);
        }
        const std::size_t nc = cand.size();
        if (nvec_ > 0) {
            Lc.assign(nc * nvec_, 0.0);
            if (owner == rank_)
                for (std::size_t q = 0; q < nc; ++q)
                    for (long k = 0; k < nvec_; ++k) Lc[q * nvec_ + k] = L_[k * nloc + (cand[q] - r0)];
            if (distributed_) world.gop.broadcast(Lc.data(), Lc.size(), owner);
        }

        // the candidates' integral columns on this rank's rows: V[q nloc + i] = (local row i|candidate q)
        t0 = wall_time();
        V.assign(nc * nloc, 0.0);
        if (nloc > 0) {
            const Tensor<double> W = ints.eri_columns(world, fp, g, bra);
            const std::size_t ldw = W.dim(1);
            for (std::size_t q = 0; q < nc; ++q) {
                const double* w = W.ptr() + (rows_[cand[q]] - fp.first[g]) * ldw;
                double* v = V.data() + q * nloc;
                for (std::size_t i = 0; i < nloc; ++i) v[i] = w[wpos[i]];
            }
        }
        stats_.ncolumns += fp.count[g];
        ++stats_.nbatches;
        stats_.t_integrals += wall_time() - t0;

        // minus the previous vectors, V -= Lc^T L: one dgemm per block of rows
        t0 = wall_time();
        if (nvec_ > 0)
            local_blocks([&](const std::size_t i0, const std::size_t i1) {
                cblas::gemm(cblas::NoTrans, cblas::NoTrans, long(i1 - i0), long(nc), nvec_, -1.0, L_.data() + i0,
                            long(nloc), Lc.data(), nvec_, 1.0, V.data() + i0, long(nloc));
            });

        // the candidates' block of V, from the owner: every rank decides the pivots of the batch from it
        vcc.assign(nc * nc, 0.0);
        if (owner == rank_)
            for (std::size_t q = 0; q < nc; ++q)
                for (std::size_t q2 = 0; q2 < nc; ++q2) vcc[q * nc + q2] = V[q * nloc + (cand[q2] - r0)];
        if (distributed_) world.gop.broadcast(vcc.data(), vcc.size(), owner);

        // the candidates in descending residual diagonal, re-checked after each new vector. vcc and
        // dc are updated exactly as the owner's rows of V and D are, so they stay bitwise equal.
        std::vector<char> done(nc, 0);
        f.resize(nc);
        while (true) {
            std::size_t qp = nc;
            double dq = dmin;
            for (std::size_t q = 0; q < nc; ++q)
                if (not done[q] and dc[q] > dq) {
                    dq = dc[q];
                    qp = q;
                }
            if (qp == nc) break;
            done[qp] = 1;

            // L_new = V[qp]/sqrt(D); f[q] its value on the remaining candidates
            const double s = 1.0 / std::sqrt(dq);
            for (std::size_t q = 0; q < nc; ++q) f[q] = done[q] ? 0.0 : vcc[qp * nc + q] * s;
            L_.resize((nvec_ + 1) * nloc);
            double* l = L_.data() + nvec_ * nloc;
            const double* v = V.data() + qp * nloc;
            local_blocks([&](const std::size_t i0, const std::size_t i1) {
                for (std::size_t i = i0; i < i1; ++i) {
                    l[i] = v[i] * s;
                    D[i] = std::max(D[i] - l[i] * l[i], 0.0);
                }
                for (std::size_t q = 0; q < nc; ++q) {
                    if (done[q]) continue;
                    double* vq = V.data() + q * nloc;
                    for (std::size_t i = i0; i < i1; ++i) vq[i] -= f[q] * l[i];
                }
            });
            if (owner == rank_) D[cand[qp] - r0] = 0.0;
            for (std::size_t q = 0; q < nc; ++q) {
                if (done[q]) continue;
                dc[q] = std::max(dc[q] - f[q] * f[q], 0.0);
                for (std::size_t q2 = 0; q2 < nc; ++q2) {
                    const double lq2 = vcc[qp * nc + q2] * s;
                    vcc[q * nc + q2] -= f[q] * lq2;
                }
            }
            dc[qp] = 0.0;
            ++nvec_;
        }
        stats_.t_updates += wall_time() - t0;
    }

    // the vectors in chunks of consecutive vectors, complete over the kept rows: per chunk, every
    // rank adds its rows into a zero buffer, the sum is exact (one contribution per element)
    if (distributed_) {
        t0 = wall_time();
        chunks_.resize(nchunk_);
        for (long c = 0; c < nchunk_; ++c) {
            const auto [k0, k1] = chunk(c);
            std::vector<double> buf((k1 - k0) * nk, 0.0);
            for (long k = k0; k < k1; ++k)
                std::copy(L_.data() + k * nloc, L_.data() + (k + 1) * nloc, buf.data() + (k - k0) * nk + r0);
            if (not buf.empty()) world.gop.sum(buf.data(), buf.size());
            if (holds(c)) chunks_[c] = std::move(buf);
        }
        std::vector<double>().swap(L_);
        stats_.t_redistribute = wall_time() - t0;
    }
}


const double* CholeskyERIDecomposition::chunk_vectors(const long c) const {
    MADNESS_CHECK_THROW(c >= 0 and c < nchunk_ and holds(c), "CholeskyERIDecomposition: chunk not held by this rank");
    return distributed_ ? chunks_[c].data() : L_.data() + chunk(c).first * long(nkept());
}


const double* CholeskyERIDecomposition::vectors() const {
    MADNESS_CHECK_THROW(not distributed_, "CholeskyERIDecomposition: the vectors are distributed over the ranks");
    return L_.data();
}


unsigned long CholeskyERIDecomposition::hash(World& world) const {
    const std::size_t nk = nkept();
    unsigned long h = 0;
    for (long c = 0; c < nchunk_; ++c) {
        if (not holds(c)) continue;
        const auto [k0, k1] = chunk(c);
        const double* l = chunk_vectors(c);
        for (long k = k0; k < k1; ++k)
            for (std::size_t i = 0; i < nk; ++i) {
                unsigned long x = 0;
                if (l[(k - k0) * nk + i] != 0.0) std::memcpy(&x, l + (k - k0) * nk + i, sizeof(x));   // -0 as +0
                const unsigned r = unsigned((std::size_t(k) * 31 + i * 17) % 64);
                h ^= (r == 0) ? x : ((x << r) | (x >> (64 - r)));
            }
    }
    if (distributed_) world.gop.bit_xor(&h, 1);
    return h;
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
    const CholeskyERIDecomposition& chol = *chol_;
    const long m = chol.nvec(), nchunk = chol.nchunk();
    const std::size_t nk = chol.nkept();
    const bool dist = chol.distributed();
    std::vector<long> mine;     // the non-empty chunks this rank holds
    for (long c = 0; c < nchunk; ++c)
        if (chol.holds(c) and chol.chunk(c).second > chol.chunk(c).first) mine.push_back(c);

    // J: gamma = L^T P~ per chunk, P~ the total density with both orders of mu != nu. gamma of all
    // vectors and the partial J~ = L gamma of every chunk are gathered (one contribution per element,
    // so the sums are exact) and J~ is added in chunk order on every rank.
    const Tensor<double> pt = open ? Pa + Pb : 2.0 * Pa;
    std::vector<double> ptk(nk), gamma(std::max(m, 1L), 0.0), jpart(nchunk * nk, 0.0), jt(nk, 0.0);
    for (std::size_t i = 0; i < nk; ++i) ptk[i] = (mu_[i] == nu_[i] ? 1.0 : 2.0) * pt(mu_[i], nu_[i]);
    for_each_block(world_, mine.size(), [&](const std::size_t t) {
        const auto [k0, k1] = chol.chunk(mine[t]);
        cblas::gemv(cblas::Trans, long(nk), k1 - k0, 1.0, chol.chunk_vectors(mine[t]), long(nk), ptk.data(), 1, 0.0,
                    gamma.data() + k0, 1);
    });
    if (dist and m > 0) world_.gop.sum(gamma.data(), m);
    for_each_block(world_, mine.size(), [&](const std::size_t t) {
        const auto [k0, k1] = chol.chunk(mine[t]);
        cblas::gemv(cblas::NoTrans, long(nk), k1 - k0, 1.0, chol.chunk_vectors(mine[t]), long(nk), gamma.data() + k0,
                    1, 0.0, jpart.data() + mine[t] * nk, 1);
    });
    if (dist and not jpart.empty()) world_.gop.sum(jpart.data(), jpart.size());
    for (long c = 0; c < nchunk; ++c)
        for (std::size_t i = 0; i < nk; ++i) jt[i] += jpart[c * nk + i];
    J = Tensor<double>(n, n);
    for (std::size_t i = 0; i < nk; ++i) J(mu_[i], nu_[i]) = J(nu_[i], mu_[i]) = jt[i];

    // K per spin from the factors of its density: chunk c sums X_k X_k^T over its vectors, several
    // vectors per dgemm; the partials of all chunks are gathered and added in chunk order
    std::vector<DensityFactor> factors;
    factor_density(Pa, 0, factors);
    if (open) factor_density(Pb, 1, factors);
    const int nspin = open ? 2 : 1;
    long rmax = 1;
    for (const DensityFactor& f : factors) rmax = std::max(rmax, f.r);
    const long nstack = std::max(1L, 512 / rmax);
    const std::size_t nn = std::size_t(n) * n;
    std::vector<double> kpart(std::size_t(nchunk) * nspin * nn, 0.0);
    for_each_block(world_, mine.size(), [&](const std::size_t t) {
        const long c = mine[t];
        const auto [k0, k1] = chol.chunk(c);
        const double* L = chol.chunk_vectors(c);
        std::vector<double> Lk(nn, 0.0);    // the vector as a matrix; screened pairs stay 0
        std::vector<std::vector<double>> X(factors.size());
        for (std::size_t f = 0; f < factors.size(); ++f) X[f].resize(n * factors[f].r * nstack);
        for (long k = k0; k < k1; k += nstack) {
            const long kn = std::min(nstack, k1 - k);
            for (long kk = 0; kk < kn; ++kk) {
                const double* l = L + (k - k0 + kk) * nk;
                for (std::size_t i = 0; i < nk; ++i) Lk[mu_[i] * n + nu_[i]] = Lk[nu_[i] * n + mu_[i]] = l[i];
                for (std::size_t f = 0; f < factors.size(); ++f)
                    cblas::gemm(cblas::NoTrans, cblas::NoTrans, n, factors[f].r, n, 1.0, Lk.data(), n,
                                factors[f].Y.data(), n, 0.0, X[f].data() + kk * n * factors[f].r, n);
            }
            for (std::size_t f = 0; f < factors.size(); ++f)
                cblas::gemm(cblas::NoTrans, cblas::Trans, n, n, kn * factors[f].r, factors[f].sign, X[f].data(), n,
                            X[f].data(), n, 1.0, kpart.data() + (std::size_t(c) * nspin + factors[f].spin) * nn, n);
        }
    });
    if (dist and not kpart.empty()) world_.gop.sum(kpart.data(), kpart.size());
    const auto sum_chunks = [&](const int s) {
        Tensor<double> K(n, n);
        double* k = K.ptr();
        for (long c = 0; c < nchunk; ++c) {
            const double* kp = kpart.data() + (std::size_t(c) * nspin + s) * nn;
            for (std::size_t i = 0; i < nn; ++i) k[i] += kp[i];
        }
        return Tensor<double>(0.5 * (K + transpose(K)));
    };
    Ka = sum_chunks(0);
    Kb = open ? sum_chunks(1) : Tensor<double>();
}

} // namespace lcao
} // namespace madness
