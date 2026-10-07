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
#include <madness/chem/lcao_cholesky.h>
#include <madness/chem/molecular_functors.h>
#include <madness/mra/vmra.h>
#include <madness/tensor/tensor_lapack.h>
#include <madness/world/print.h>

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <deque>

namespace madness {
namespace lcao {

namespace {

/// J' and K' of InCoreERI::jk from the stored integrals of the pairs ij in [begin, end)
///
/// pt: total density in pair order (for J); pa, pb: spin densities as full matrices (for K),
/// pb null for a closed shell
void jk_pairs(const double* g, const long* pi, const long* pj, const double* pt, const double* pa, const double* pb,
              const long n, const long begin, const long end, double* jpair, double* kha, double* khb) {
    g += begin * (begin + 1) / 2;           // pair ij holds ij+1 values
    for (long ij = begin; ij < end; ++ij) {
        const long i = pi[ij], j = pj[ij];
        const double fij = (i == j) ? 0.5 : 1.0;
        const double* pa_i = pa + i * n;
        const double* pa_j = pa + j * n;
        double* ka_i = kha + i * n;
        double* ka_j = kha + j * n;
        double jij = 0.0;
        for (long kl = 0; kl <= ij; ++kl) {
            const long k = pi[kl], l = pj[kl];
            const double f = fij * ((k == l) ? 0.5 : 1.0) * ((kl == ij) ? 0.5 : 1.0);
            const double v = f * g[kl];
            jij += v * pt[kl];
            jpair[kl] += v * pt[ij];
            ka_i[k] += v * pa_j[l];
            ka_j[k] += v * pa_i[l];
            ka_i[l] += v * pa_j[k];
            ka_j[l] += v * pa_i[k];
        }
        jpair[ij] += jij;
        if (pb) {
            const double* pb_i = pb + i * n;
            const double* pb_j = pb + j * n;
            double* kb_i = khb + i * n;
            double* kb_j = khb + j * n;
            for (long kl = 0; kl <= ij; ++kl) {
                const long k = pi[kl], l = pj[kl];
                const double f = fij * ((k == l) ? 0.5 : 1.0) * ((kl == ij) ? 0.5 : 1.0);
                const double v = f * g[kl];
                kb_i[k] += v * pb_j[l];
                kb_j[k] += v * pb_i[l];
                kb_i[l] += v * pb_j[k];
                kb_j[l] += v * pb_i[k];
            }
        }
        g += ij + 1;
    }
}

/// one range of pairs of InCoreERI::jk as a task, with its own accumulators
class JKTask : public TaskInterface {
public:
    JKTask(const double* g, const long* pi, const long* pj, const double* pt, const double* pa, const double* pb,
           const long n, const long begin, const long end, double* jpair, double* kha, double* khb)
        : g_(g), pi_(pi), pj_(pj), pt_(pt), pa_(pa), pb_(pb), n_(n), begin_(begin), end_(end), jpair_(jpair),
          kha_(kha), khb_(khb) {}

    using TaskInterface::run;
    void run(World&) override { jk_pairs(g_, pi_, pj_, pt_, pa_, pb_, n_, begin_, end_, jpair_, kha_, khb_); }

private:
    const double *g_;
    const long *pi_, *pj_;
    const double *pt_, *pa_, *pb_;
    const long n_, begin_, end_;
    double *jpair_, *kha_, *khb_;
};

} // namespace


void InCoreERI::jk(const Tensor<double>& Pa, const Tensor<double>& Pb, Tensor<double>& J, Tensor<double>& Ka,
                   Tensor<double>& Kb) const {
    // Each stored (ij|kl) stands for up to 8 index orders. Every order adds g P to one
    // element of J and one of K; the factor f removes the orders that coincide (i = j,
    // k = l, ij = kl). Half of the orders are transposes of the other half, so J' and K'
    // collect one half and J = J' + J'^T, K = K' + K'^T.
    const long n = Pa.dim(0);
    const bool open = Pb.size() > 0;
    MADNESS_CHECK_THROW(eri_->nbf() == n, "InCoreERI: integrals and density matrix do not match");
    const long npair = n * (n + 1) / 2;
    std::vector<long> pi(npair), pj(npair);
    for (long i = 0, ij = 0; i < n; ++i)
        for (long j = 0; j <= i; ++j, ++ij) {
            pi[ij] = i;
            pj[ij] = j;
        }
    // contiguous copies, for the raw-pointer loops below
    const Tensor<double> pa = copy(Pa);
    const Tensor<double> pb = open ? copy(Pb) : Tensor<double>();
    const Tensor<double> pt = open ? pa + pb : 2.0 * pa;
    std::vector<double> ptpair(npair);
    for (long ij = 0; ij < npair; ++ij) ptpair[ij] = pt(pi[ij], pj[ij]);

    // chunks of equal work (pair ij holds ij+1 values), each with its own J' and K',
    // summed in a fixed order so that the result does not depend on the scheduling
    const long nchunk = std::min(16L, npair);
    std::vector<std::vector<double>> jparts(nchunk, std::vector<double>(npair, 0.0));
    std::vector<Tensor<double>> kaparts(nchunk), kbparts(nchunk);
    for (long c = 0; c < nchunk; ++c) {
        const long begin = long(npair * std::sqrt(double(c) / nchunk));
        const long end = (c + 1 == nchunk) ? npair : long(npair * std::sqrt(double(c + 1) / nchunk));
        kaparts[c] = Tensor<double>(n, n);
        if (open) kbparts[c] = Tensor<double>(n, n);
        world_.taskq.add(new JKTask(eri_->data(), pi.data(), pj.data(), ptpair.data(), pa.ptr(),
                                    open ? pb.ptr() : nullptr, n, begin, end, jparts[c].data(), kaparts[c].ptr(),
                                    open ? kbparts[c].ptr() : nullptr));
    }
    world_.taskq.fence();
    std::vector<double> jpair(npair, 0.0);
    Tensor<double> Kha(n, n), Khb;
    if (open) Khb = Tensor<double>(n, n);
    for (long c = 0; c < nchunk; ++c) {
        for (long ij = 0; ij < npair; ++ij) jpair[ij] += jparts[c][ij];
        Kha += kaparts[c];
        if (open) Khb += kbparts[c];
    }

    // J' holds J'_ij for i >= j, already times 2 for (ij|kl) and (ij|lk); K' is full
    J = Tensor<double>(n, n);
    for (long ij = 0; ij < npair; ++ij) {
        const long i = pi[ij], j = pj[ij];
        J(i, j) += 2.0 * jpair[ij];
        J(j, i) += 2.0 * jpair[ij];
    }
    Ka = Kha + transpose(Kha);
    Kb = open ? Khb + transpose(Khb) : Tensor<double>();
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

/// the elements of a and b one after the other, as one vector (DIIS over both spins)
Tensor<double> stack(const Tensor<double>& a, const Tensor<double>& b) {
    Tensor<double> ab(a.size() + b.size());
    const Tensor<double> ac = copy(a), bc = copy(b);
    std::copy(ac.ptr(), ac.ptr() + ac.size(), ab.ptr());
    std::copy(bc.ptr(), bc.ptr() + bc.size(), ab.ptr() + ac.size());
    return ab;
}

/// undo stack for two n x n matrices
void unstack(const Tensor<double>& ab, const long n, Tensor<double>& a, Tensor<double>& b) {
    a = Tensor<double>(n, n);
    b = Tensor<double>(n, n);
    std::copy(ab.ptr(), ab.ptr() + n * n, a.ptr());
    std::copy(ab.ptr() + n * n, ab.ptr() + 2 * n * n, b.ptr());
}

} // namespace


LCAOSCF::LCAOSCF(World& world, const Molecule& molecule, const AtomicBasisSet& aobasis, const int nalpha,
                 const int nbeta, const LCAOParameters& param, const bool collective)
    : world_(world), molecule_(molecule), aobasis_(aobasis), nalpha_(nalpha), nbeta_(nbeta), param_(param),
      collective_(collective) {
    for (std::size_t i = 0; i < molecule_.natom(); ++i)
        MADNESS_CHECK_THROW(not molecule_.get_atom(i).pseudo_atom, "LCAOSCF: pseudo-atoms are not supported");
    MADNESS_CHECK_THROW(molecule_.n_core_orb_all() == 0, "LCAOSCF: core potentials are not supported");
    MADNESS_CHECK_THROW(nalpha_ > 0 and nbeta_ >= 0 and nalpha_ >= nbeta_,
                        "LCAOSCF: need nalpha > 0 and 0 <= nbeta <= nalpha");
    shells_ = make_shells(molecule_, aobasis_);
}


void LCAOSCF::compute_integrals() {
    const bool printme = world_.rank() == 0 and param_.print_level() > 0;
    const double t0 = wall_time();
    // the fit of 1/r must cover the molecule (12_parallel_cholesky_plan.md, section 3.5): beyond kernel_hi
    // the attraction to far nuclei and the repulsion of far electrons vanish, while the nuclear repulsion
    // stays exact. So the range is at least the largest interatomic distance plus 10 bohr.
    double extent = 0.0;
    for (std::size_t i = 0; i < molecule_.natom(); ++i)
        for (std::size_t j = 0; j < i; ++j) {
            const Atom& a = molecule_.get_atom(i);
            const Atom& b = molecule_.get_atom(j);
            extent = std::max(extent, std::sqrt((a.x - b.x) * (a.x - b.x) + (a.y - b.y) * (a.y - b.y) +
                                                (a.z - b.z) * (a.z - b.z)));
        }
    const double hi = std::max(param_.kernel_hi(), extent + 10.0);
    kernel_ = GaussianKernel::coulomb(param_.kernel_lo(), hi, param_.kernel_eps());
    const GaussianKernel& kernel = kernel_;
    if (printme) {
        printf("Coulomb kernel: %zu Gaussians for [%.1e, %.1e] bohr at relative precision %.1e\n",
               kernel.size(), param_.kernel_lo(), hi, param_.kernel_eps());
        if (hi > param_.kernel_hi())
            printf("    kernel_hi %.1f raised to cover the molecule: largest interatomic distance %.1f bohr + 10\n",
                   param_.kernel_hi(), extent);
        if (param_.print_level() > 1)
            printf("    max |sum_m w_m exp(-t_m r^2) r - 1| on [lo, hi]: %.2e\n",
                   kernel.max_relative_coulomb_error(param_.kernel_lo(), hi));
    }
    const SeparatedGaussianIntegrals ints(shells_, kernel);
    S_ = ints.overlap(world_);
    T_ = ints.kinetic(world_);
    V_ = ints.nuclear_attraction(world_, molecule_);
    H_ = T_ + V_;
    const double t1 = wall_time();
    if (param_.eri() == "cholesky") {
        const auto chol =
            std::make_shared<const CholeskyERIDecomposition>(world_, ints, param_.cholesky_tol(), 1.e-2, collective_);
        const double t2 = wall_time();
        twoe_ = std::make_unique<CholeskyERI>(world_, chol, S_, param_.lindep());
        const unsigned long hash = chol->hash(world_);
        if (printme) {
            const CholeskyERIDecomposition::Stats& s = chol->stats();
            printf("integrals: one-electron %.2fs, two-electron %.2fs\n", t1 - t0, t2 - t1);
            printf("    Cholesky decomposition at cholesky_tol %.0e (all kernel terms, no Schwarz skips): %zu of %zu "
                   "function pairs kept, %ld vectors = %.2f N, %.2f GB in all\n", chol->tol(), s.nkept, s.npairs,
                   chol->nvec(), double(chol->nvec()) / double(S_.dim(0)), 8.0e-9 * double(chol->nvec()) * s.nkept);
            printf("    %zu integral columns in %zu batches; diagonal %.2fs, integrals %.2fs, updates %.2fs\n",
                   s.ncolumns, s.nbatches, s.t_diagonal, s.t_integrals, s.t_updates);
            printf("    %s, vectors hash %016lx\n", chol->distributed() ? "distributed over the ranks" : "on one rank",
                   hash);
        }
        return;
    }
    ERIStats stats;
    eri_ = std::make_shared<const PackedERI>(ints.eri(world_, param_.kernel_screen(), param_.schwarz(), &stats));
    const double t2 = wall_time();
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
    MADNESS_CHECK_THROW(nmo >= nalpha_, "LCAOSCF: too few linearly independent basis functions");
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


Tensor<double> LCAOSCF::diagonalize(const Tensor<double>& F, const int nocc, Tensor<double>& C,
                                    Tensor<double>& eps) const {
    if (not collective_ or world_.rank() == 0) {
        Tensor<double> Fp = inner(transpose(X_), inner(F, X_));
        Fp = 0.5 * (Fp + transpose(Fp));
        Tensor<double> Cp;
        syev(Fp, Cp, eps);
        C = inner(X_, Cp);
    }
    if (collective_ and world_.size() > 1) {
        world_.gop.broadcast_serializable(C, 0);
        world_.gop.broadcast_serializable(eps, 0);
    }
    if (nocc == 0) return Tensor<double>(F.dim(0), F.dim(1));
    const Tensor<double> Cocc = copy(C(_, Slice(0, nocc - 1)));
    return inner(Cocc, transpose(Cocc));
}


Tensor<double> LCAOSCF::commutator_error(const Tensor<double>& F, const Tensor<double>& P) const {
    // in the orthogonal basis its size does not depend on the scaling of the basis functions
    const Tensor<double> FPS = inner(F, inner(P, S_));
    return inner(transpose(X_), inner(FPS - transpose(FPS), X_));
}


double LCAOSCF::solve() {
    const bool printme = world_.rank() == 0 and param_.print_level() > 0;
    const bool open = not restricted();
    compute_integrals();
    make_orthogonalizer();

    // starting spin densities: the atomic guess density split by occupation, or the core hamiltonian
    std::string guess = param_.guess();
    Tensor<double> Pa, Pb;
    if (guess == "sad") {
        const Tensor<double> P = sad_density();
        if (P.size() == 0) {
            if (world_.rank() == 0) print("the basis file has no atomic guess densities; using the core hamiltonian");
            guess = "core";
        } else {
            Pa = (double(nalpha_) / (nalpha_ + nbeta_)) * P;
            Pb = (double(nbeta_) / (nalpha_ + nbeta_)) * P;
        }
    }
    if (guess == "core") {
        Pa = diagonalize(H_, nalpha_, Ca_, epsa_);
        Pb = diagonalize(H_, nbeta_, Cb_, epsb_);
    }
    if (printme)
        printf("starting density: %s, %.6f electrons (%d alpha, %d beta, %s)\n", guess.c_str(), (Pa + Pb).trace(S_),
               nalpha_, nbeta_, open ? "UHF" : "RHF");

    const double enuc = molecule_.nuclear_repulsion_energy();
    const double damping = param_.damping();
    const long n = S_.dim(0);
    const Tensor<double> none;              // an empty Pb: closed shell
    Tensor<double> J, Ka, Kb;
    std::deque<Tensor<double>> diis_f, diis_e;   // with collective, on rank 0 only
    const bool root = not collective_ or world_.rank() == 0;
    double eold = 0.0;
    double tjk = 0.0, tdiag = 0.0;     // wall times of the J/K builds and of DIIS + diagonalization
    converged_ = false;
    if (printme) printf("\n iter          energy            dE        rms(dP)    max|FPS-SPF|\n");
    for (int iter = 0; iter < param_.maxiter(); ++iter) {
        const double tj0 = wall_time();
        twoe_->jk(Pa, open ? Pb : none, J, Ka, Kb);
        tjk += wall_time() - tj0;
        const double td0 = wall_time();
        Tensor<double> Fa = H_ + J - Ka;
        Tensor<double> Fb = open ? H_ + J - Kb : Fa;
        const double etot = 0.5 * (Pa + Pb).trace(H_) + 0.5 * Pa.trace(Fa) + 0.5 * Pb.trace(Fb) + enuc;
        double errmax = 0.0;
        if (not root) {
            // rank 0 extrapolates and diagonalizes, and broadcasts the orbitals
        } else if (open) {
            // one DIIS over both spins: shared coefficients for Fa and Fb
            const Tensor<double> ea = commutator_error(Fa, Pa), eb = commutator_error(Fb, Pb);
            errmax = std::max(ea.absmax(), eb.absmax());
            if (param_.diis() > 0)
                unstack(diis_extrapolate(stack(Fa, Fb), stack(ea, eb), param_.diis(), diis_f, diis_e), n, Fa, Fb);
        } else {
            const Tensor<double> e = commutator_error(Fa, Pa + Pb);
            errmax = e.absmax();
            if (param_.diis() > 0) Fa = diis_extrapolate(Fa, e, param_.diis(), diis_f, diis_e);
        }
        const Tensor<double> Pa_new = diagonalize(Fa, nalpha_, Ca_, epsa_);
        Tensor<double> Pb_new;
        if (open) {
            Pb_new = diagonalize(Fb, nbeta_, Cb_, epsb_);
        } else {
            Pb_new = Pa_new;
            Cb_ = Ca_;
            epsb_ = epsa_;
        }
        tdiag += wall_time() - td0;
        const double drms = ((Pa_new - Pa).normf() + (Pb_new - Pb).normf()) / double(n);
        const double de = etot - eold;
        if (printme) printf("%5d  %18.10f  %12.4e  %12.4e  %12.4e\n", iter, etot, de, drms, errmax);
        Pa = (damping > 0.0) ? (1.0 - damping) * Pa_new + damping * Pa : Pa_new;
        Pb = (damping > 0.0) ? (1.0 - damping) * Pb_new + damping * Pb : Pb_new;
        eold = etot;
        iterations_ = iter + 1;
        int stop = (iter > 0 and std::abs(de) < param_.econv() and drms < param_.dconv()) ? 1 : 0;
        if (collective_ and world_.size() > 1) world_.gop.broadcast(stop, 0);
        if (stop) {
            converged_ = true;
            break;
        }
    }
    Pa_ = Pa;
    Pb_ = Pb;

    // the energy and its parts for the final densities
    twoe_->jk(Pa_, open ? Pb_ : none, J, Ka, Kb);
    if (not open) Kb = Ka;
    const Tensor<double> P = Pa_ + Pb_;
    energies_.kinetic = P.trace(T_);
    energies_.nuclear_attraction = P.trace(V_);
    energies_.coulomb = 0.5 * P.trace(J);
    energies_.exchange = -0.5 * (Pa_.trace(Ka) + Pb_.trace(Kb));
    energies_.nuclear_repulsion = enuc;
    energies_.total = energies_.kinetic + energies_.nuclear_attraction + energies_.coulomb
                    + energies_.exchange + energies_.nuclear_repulsion;

    // <S^2> = Sz(Sz+1) + nbeta - sum_ij |<a_i|b_j>|^2, the last sum being tr(Pa S Pb S)
    const double sz = 0.5 * (nalpha_ - nbeta_);
    s2_ = open ? sz * (sz + 1.0) + nbeta_ - inner(Pa_, S_).trace(transpose(inner(Pb_, S_))) : 0.0;

    if (printme) {
        printf("\n%s after %d iterations, %.6f electrons\n", converged_ ? "converged" : "NOT CONVERGED",
               iterations_, P.trace(S_));
        printf("SCF time: J/K %.2fs in %d builds (%.3fs each), DIIS and diagonalization %.2fs\n", tjk, iterations_,
               tjk / std::max(iterations_, 1), tdiag);
        if (open) printf("<S^2> = %.6f (pure spin state: %.6f)\n", s2_, sz * (sz + 1.0));
        printf("\n              kinetic %16.8f\n", energies_.kinetic);
        printf("   nuclear attraction %16.8f\n", energies_.nuclear_attraction);
        printf("              coulomb %16.8f\n", energies_.coulomb);
        printf(" exchange-correlation %16.8f\n", energies_.exchange);
        printf("    nuclear-repulsion %16.8f\n", energies_.nuclear_repulsion);
        printf("                total %16.8f\n\n", energies_.total);
        const auto print_eps = [](const char* label, const Tensor<double>& eps, const int nocc) {
            const long nprint = std::min(eps.size(), long(nocc + 5));
            printf("%s orbital energies (occupied, then the lowest virtuals):\n", label);
            for (long i = 0; i < nprint; ++i) printf("%5ld %14.8f%s\n", i, eps(i), i < nocc ? "" : "  (virtual)");
        };
        if (open) {
            print_eps("alpha", epsa_, nalpha_);
            print_eps("beta", epsb_, nbeta_);
        } else {
            print_eps("closed-shell", epsa_, nalpha_);
        }
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
