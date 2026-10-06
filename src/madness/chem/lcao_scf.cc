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

/// \file lcao_scf.cc
/// \brief closed-shell Hartree-Fock in a Gaussian basis, with integrals from separated kernels

#include <madness/chem/lcao_scf.h>
#include <madness/chem/molecular_functors.h>
#include <madness/mra/vmra.h>
#include <madness/tensor/tensor_lapack.h>
#include <madness/world/print.h>

#include <cmath>
#include <cstdio>
#include <deque>

namespace madness {
namespace lcao {

namespace {

/// J' and K' of InCoreERI::jk from the stored integrals of the pairs ij in [begin, end)
void jk_pairs(const double* g, const long* pi, const long* pj, const double* p, const double* ppair, const long n,
              const long begin, const long end, double* jpair, double* kh) {
    g += begin * (begin + 1) / 2;           // pair ij holds ij+1 values
    for (long ij = begin; ij < end; ++ij) {
        const long i = pi[ij], j = pj[ij];
        const double fij = (i == j) ? 0.5 : 1.0;
        const double* pr_i = p + i * n;
        const double* pr_j = p + j * n;
        double* kr_i = kh + i * n;
        double* kr_j = kh + j * n;
        double jij = 0.0;
        for (long kl = 0; kl <= ij; ++kl) {
            const long k = pi[kl], l = pj[kl];
            const double f = fij * ((k == l) ? 0.5 : 1.0) * ((kl == ij) ? 0.5 : 1.0);
            const double v = f * g[kl];
            jij += v * ppair[kl];
            jpair[kl] += v * ppair[ij];
            kr_i[k] += v * pr_j[l];
            kr_j[k] += v * pr_i[l];
            kr_i[l] += v * pr_j[k];
            kr_j[l] += v * pr_i[k];
        }
        jpair[ij] += jij;
        g += ij + 1;
    }
}

/// one range of pairs of InCoreERI::jk as a task, with its own accumulators
class JKTask : public TaskInterface {
public:
    JKTask(const double* g, const long* pi, const long* pj, const double* p, const double* ppair, const long n,
           const long begin, const long end, double* jpair, double* kh)
        : g_(g), pi_(pi), pj_(pj), p_(p), ppair_(ppair), n_(n), begin_(begin), end_(end), jpair_(jpair), kh_(kh) {}

    using TaskInterface::run;
    void run(World&) override { jk_pairs(g_, pi_, pj_, p_, ppair_, n_, begin_, end_, jpair_, kh_); }

private:
    const double *g_;
    const long *pi_, *pj_;
    const double *p_, *ppair_;
    const long n_, begin_, end_;
    double *jpair_, *kh_;
};

} // namespace


void InCoreERI::jk(const Tensor<double>& P, Tensor<double>& J, Tensor<double>& K) const {
    // Each stored (ij|kl) stands for up to 8 index orders. Every order adds g P to one
    // element of J and one of K; the factor f removes the orders that coincide (i = j,
    // k = l, ij = kl). Half of the orders are transposes of the other half, so J' and K'
    // collect one half and J = J' + J'^T, K = K' + K'^T.
    const long n = P.dim(0);
    MADNESS_CHECK_THROW(eri_->nbf() == n, "InCoreERI: integrals and density matrix do not match");
    const long npair = n * (n + 1) / 2;
    std::vector<long> pi(npair), pj(npair);
    for (long i = 0, ij = 0; i < n; ++i)
        for (long j = 0; j <= i; ++j, ++ij) {
            pi[ij] = i;
            pj[ij] = j;
        }
    const Tensor<double> Pc = copy(P);      // contiguous, for the raw-pointer loops below
    const double* p = Pc.ptr();
    std::vector<double> ppair(npair);
    for (long ij = 0; ij < npair; ++ij) ppair[ij] = p[pi[ij] * n + pj[ij]];

    // chunks of equal work (pair ij holds ij+1 values), each with its own J' and K',
    // summed in a fixed order so that the result does not depend on the scheduling
    const long nchunk = std::min(16L, npair);
    std::vector<std::vector<double>> jparts(nchunk, std::vector<double>(npair, 0.0));
    std::vector<Tensor<double>> kparts(nchunk);
    for (long c = 0; c < nchunk; ++c) {
        const long begin = long(npair * std::sqrt(double(c) / nchunk));
        const long end = (c + 1 == nchunk) ? npair : long(npair * std::sqrt(double(c + 1) / nchunk));
        kparts[c] = Tensor<double>(n, n);
        world_.taskq.add(new JKTask(eri_->data(), pi.data(), pj.data(), p, ppair.data(), n, begin, end,
                                    jparts[c].data(), kparts[c].ptr()));
    }
    world_.taskq.fence();
    std::vector<double> jpair(npair, 0.0);
    Tensor<double> Kh(n, n);
    for (long c = 0; c < nchunk; ++c) {
        for (long ij = 0; ij < npair; ++ij) jpair[ij] += jparts[c][ij];
        Kh += kparts[c];
    }

    // J' holds J'_ij for i >= j, already times 2 for (ij|kl) and (ij|lk); K' is full
    J = Tensor<double>(n, n);
    for (long ij = 0; ij < npair; ++ij) {
        const long i = pi[ij], j = pj[ij];
        J(i, j) += 2.0 * jpair[ij];
        J(j, i) += 2.0 * jpair[ij];
    }
    K = Kh + transpose(Kh);
}


namespace {

/// Pulay's DIIS: the combination sum_i c_i F_i of the stored Fock matrices with sum_i c_i = 1
/// that minimizes |sum_i c_i e_i|, e_i the commutator error of F_i
///
/// F and e join the subspace, which keeps the last maxsub of them. If the DIIS equations
/// are singular, the oldest entries are dropped until they are not.
Tensor<double> diis_extrapolate(const Tensor<double>& F, const Tensor<double>& e, const std::size_t maxsub,
                                std::deque<Tensor<double>>& fs, std::deque<Tensor<double>>& es) {
    fs.push_back(F);
    es.push_back(e);
    while (fs.size() > maxsub) {
        fs.pop_front();
        es.pop_front();
    }
    while (fs.size() > 1) {
        const long m = fs.size();
        Tensor<double> B(m + 1, m + 1), rhs(m + 1);
        double scale = 0.0;
        for (long i = 0; i < m; ++i) {
            for (long j = 0; j <= i; ++j) B(i, j) = B(j, i) = es[i].trace(es[j]);
            scale = std::max(scale, B(i, i));
        }
        for (long i = 0; i < m; ++i) {
            for (long j = 0; j < m; ++j) B(i, j) /= scale;
            B(i, m) = B(m, i) = -1.0;
        }
        rhs(m) = -1.0;
        Tensor<double> c;
        bool ok = scale > 0.0;
        if (ok) {
            try {
                gesv(B, rhs, c);
            } catch (...) {
                ok = false;
            }
        }
        for (long i = 0; ok and i < m; ++i) ok = std::isfinite(c(i));
        if (ok) {
            Tensor<double> Fx(F.dim(0), F.dim(1));
            for (long i = 0; i < m; ++i) Fx.gaxpy(1.0, fs[i], c(i));
            return Fx;
        }
        fs.pop_front();
        es.pop_front();
    }
    return copy(F);
}

} // namespace


LCAOSCF::LCAOSCF(World& world, const Molecule& molecule, const AtomicBasisSet& aobasis, const int nocc,
                 const LCAOParameters& param)
    : world_(world), molecule_(molecule), aobasis_(aobasis), nocc_(nocc), param_(param) {
    for (std::size_t i = 0; i < molecule_.natom(); ++i)
        MADNESS_CHECK_THROW(not molecule_.get_atom(i).pseudo_atom, "LCAOSCF: pseudo-atoms are not supported");
    MADNESS_CHECK_THROW(molecule_.n_core_orb_all() == 0, "LCAOSCF: core potentials are not supported");
    MADNESS_CHECK_THROW(nocc_ > 0, "LCAOSCF: no occupied orbitals");
    shells_ = make_shells(molecule_, aobasis_);
}


void LCAOSCF::compute_integrals() {
    const bool printme = world_.rank() == 0 and param_.print_level() > 0;
    const double t0 = wall_time();
    const GaussianKernel kernel = GaussianKernel::coulomb(param_.kernel_lo(), param_.kernel_hi(), param_.kernel_eps());
    if (printme) {
        printf("Coulomb kernel: %zu Gaussians for [%.1e, %.1e] bohr at relative precision %.1e\n",
               kernel.size(), param_.kernel_lo(), param_.kernel_hi(), param_.kernel_eps());
        if (param_.print_level() > 1)
            printf("    max |sum_m w_m exp(-t_m r^2) r - 1| on [lo, hi]: %.2e\n",
                   kernel.max_relative_coulomb_error(param_.kernel_lo(), param_.kernel_hi()));
    }
    const SeparatedGaussianIntegrals ints(shells_, kernel);
    S_ = ints.overlap();
    T_ = ints.kinetic();
    V_ = ints.nuclear_attraction(molecule_);
    const double t1 = wall_time();
    ERIStats stats;
    eri_ = std::make_shared<const PackedERI>(ints.eri(world_, param_.kernel_screen(), param_.schwarz(), &stats));
    const double t2 = wall_time();
    H_ = T_ + V_;
    twoe_ = std::make_unique<InCoreERI>(world_, eri_);
    if (printme) {
        printf("integrals: one-electron %.2fs, two-electron %.2fs\n", t1 - t0, t2 - t1);
        printf("    %ld shell-group quartets, %ld of them skipped by the Schwarz test, %.2f GB stored\n",
               stats.computed + stats.skipped, stats.skipped, 8.0e-9 * double(eri_->size()));
    }
}


void LCAOSCF::make_orthogonalizer() {
    Tensor<double> U, s;
    syev(S_, U, s);         // ascending eigenvalues, eigenvectors in the columns
    const long n = s.size();
    long first = 0;
    while (first < n and s(first) < param_.lindep()) ++first;
    const long nmo = n - first;
    MADNESS_CHECK_THROW(nmo >= nocc_, "LCAOSCF: too few linearly independent basis functions");
    X_ = Tensor<double>(n, nmo);
    for (long j = 0; j < nmo; ++j) {
        const double f = 1.0 / std::sqrt(s(first + j));
        for (long i = 0; i < n; ++i) X_(i, j) = U(i, first + j) * f;
    }
    if (world_.rank() == 0 and param_.print_level() > 0)
        printf("overlap eigenvalues: smallest %.2e, %ld below lindep %.1e dropped, %ld orbitals\n",
               s(0L), first, param_.lindep(), nmo);
}


Tensor<double> LCAOSCF::sad_density() const {
    std::vector<int> at_to_bf, at_nbf;
    aobasis_.atoms_to_bfn(molecule_, at_to_bf, at_nbf);
    const long n = S_.dim(0);
    Tensor<double> P(n, n);
    for (std::size_t iat = 0; iat < molecule_.natom(); ++iat) {
        const Tensor<double>& d = aobasis_.get_dmat(molecule_, iat);
        if (d.size() == 0) return Tensor<double>();
        MADNESS_CHECK_THROW(d.dim(0) == at_nbf[iat] and d.dim(1) == at_nbf[iat],
                            "LCAOSCF: atomic guess density does not match the basis");
        const Slice s(at_to_bf[iat], at_to_bf[iat] + at_nbf[iat] - 1);
        P(s, s) = d;
    }
    return P;
}


Tensor<double> LCAOSCF::diagonalize(const Tensor<double>& F) {
    Tensor<double> Fp = inner(transpose(X_), inner(F, X_));
    Fp = 0.5 * (Fp + transpose(Fp));
    Tensor<double> Cp, e;
    syev(Fp, Cp, e);
    C_ = inner(X_, Cp);
    eps_ = e;
    const Tensor<double> Cocc = copy(C_(_, Slice(0, nocc_ - 1)));
    return 2.0 * inner(Cocc, transpose(Cocc));
}


double LCAOSCF::solve() {
    const bool printme = world_.rank() == 0 and param_.print_level() > 0;
    compute_integrals();
    make_orthogonalizer();

    std::string guess = param_.guess();
    Tensor<double> P;
    if (guess == "sad") {
        P = sad_density();
        if (P.size() == 0) {
            if (world_.rank() == 0) print("the basis file has no atomic guess densities; using the core hamiltonian");
            guess = "core";
        }
    }
    if (guess == "core") P = diagonalize(H_);
    if (printme) printf("starting density: %s, %.6f electrons\n", guess.c_str(), P.trace(S_));

    const double enuc = molecule_.nuclear_repulsion_energy();
    const double damping = param_.damping();
    const long n = S_.dim(0);
    Tensor<double> J, K;
    std::deque<Tensor<double>> diis_f, diis_e;
    double eold = 0.0;
    converged_ = false;
    if (printme) printf("\n iter          energy            dE        rms(dP)    max|FPS-SPF|\n");
    for (int iter = 0; iter < param_.maxiter(); ++iter) {
        twoe_->jk(P, J, K);
        Tensor<double> F = H_ + J - 0.5 * K;
        const double etot = 0.5 * P.trace(H_ + F) + enuc;
        // the commutator FPS - SPF vanishes at self-consistency; in the orthogonal basis
        // its size does not depend on the scaling of the basis functions
        const Tensor<double> FPS = inner(F, inner(P, S_));
        const Tensor<double> err = inner(transpose(X_), inner(FPS - transpose(FPS), X_));
        if (param_.diis() > 0) F = diis_extrapolate(F, err, param_.diis(), diis_f, diis_e);
        const Tensor<double> Pnew = diagonalize(F);
        const double drms = (Pnew - P).normf() / double(n);
        const double de = etot - eold;
        if (printme) printf("%5d  %18.10f  %12.4e  %12.4e  %12.4e\n", iter, etot, de, drms, err.absmax());
        P = (damping > 0.0) ? (1.0 - damping) * Pnew + damping * P : Pnew;
        eold = etot;
        iterations_ = iter + 1;
        if (iter > 0 and std::abs(de) < param_.econv() and drms < param_.dconv()) {
            converged_ = true;
            break;
        }
    }
    P_ = P;

    // the energy and its parts for the final density
    twoe_->jk(P_, J, K);
    energies_.kinetic = P_.trace(T_);
    energies_.nuclear_attraction = P_.trace(V_);
    energies_.coulomb = 0.5 * P_.trace(J);
    energies_.exchange = -0.25 * P_.trace(K);
    energies_.nuclear_repulsion = enuc;
    energies_.total = energies_.kinetic + energies_.nuclear_attraction + energies_.coulomb
                    + energies_.exchange + energies_.nuclear_repulsion;

    if (printme) {
        printf("\n%s after %d iterations, %.6f electrons\n", converged_ ? "converged" : "NOT CONVERGED",
               iterations_, P_.trace(S_));
        printf("\n              kinetic %16.8f\n", energies_.kinetic);
        printf("   nuclear attraction %16.8f\n", energies_.nuclear_attraction);
        printf("              coulomb %16.8f\n", energies_.coulomb);
        printf(" exchange-correlation %16.8f\n", energies_.exchange);
        printf("    nuclear-repulsion %16.8f\n", energies_.nuclear_repulsion);
        printf("                total %16.8f\n\n", energies_.total);
        const long nprint = std::min(eps_.size(), long(nocc_ + 5));
        printf("orbital energies (occupied, then the lowest virtuals):\n");
        for (long i = 0; i < nprint; ++i) printf("%5ld %14.8f%s\n", i, eps_(i), i < nocc_ ? "" : "  (virtual)");
    }
    return energies_.total;
}


std::vector<Function<double,3>> project_orbitals(World& world, const Molecule& molecule,
                                                  const AtomicBasisSet& aobasis, const Tensor<double>& C,
                                                  const long nmo) {
    const int nbf = aobasis.nbf(molecule);
    MADNESS_CHECK_THROW(C.ndim() == 2 and C.dim(0) == nbf and C.dim(1) >= nmo,
                        "project_orbitals: coefficients do not match the basis");
    // the basis functions as AtomicBasisSet defines them, not normalized: C refers to these
    std::vector<Function<double,3>> ao(nbf);
    for (int mu = 0; mu < nbf; ++mu) {
        const std::shared_ptr<FunctionFunctorInterface<double,3>> f =
            std::make_shared<madchem::AtomicBasisFunctor>(aobasis.get_atomic_basis_function(molecule, mu));
        ao[mu] = FunctionFactory<double,3>(world).functor(f).truncate_on_project().nofence();
    }
    world.gop.fence();
    std::vector<Function<double,3>> mo = transform(world, ao, copy(C(_, Slice(0, nmo - 1))));
    truncate(world, mo);
    return mo;
}

} // namespace lcao
} // namespace madness
