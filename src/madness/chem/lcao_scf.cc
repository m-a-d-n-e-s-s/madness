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
#include <madness/tensor/tensor_lapack.h>
#include <madness/world/print.h>

#include <cmath>
#include <cstdio>

namespace madness {
namespace lcao {

void InCoreERI::jk(const Tensor<double>& P, Tensor<double>& J, Tensor<double>& K) const {
    const long n = P.dim(0);
    MADNESS_CHECK_THROW(eri_.ndim() == 4 and eri_.dim(0) == n and eri_.iscontiguous(),
                        "InCoreERI: integrals and density matrix do not match");
    const Tensor<double> Pc = copy(P);      // contiguous, for the raw-pointer loops below
    const double* g = eri_.ptr();
    const double* p = Pc.ptr();
    J = Tensor<double>(n, n);
    K = Tensor<double>(n, n);
    for (long mu = 0; mu < n; ++mu) {
        for (long nu = 0; nu < n; ++nu) {
            // J: (mu nu|l s) is contiguous in (l s)
            const double* gmn = g + (mu * n + nu) * n * n;
            double j = 0.0;
            for (long ls = 0; ls < n * n; ++ls) j += gmn[ls] * p[ls];
            // K: (mu l|nu s) P_{ls}
            double k = 0.0;
            for (long l = 0; l < n; ++l) {
                const double* gml = g + ((mu * n + l) * n + nu) * n;
                const double* pl = p + l * n;
                for (long s = 0; s < n; ++s) k += gml[s] * pl[s];
            }
            J(mu, nu) = j;
            K(mu, nu) = k;
        }
    }
}


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
    eri_ = ints.eri();
    const double t2 = wall_time();
    H_ = T_ + V_;
    twoe_ = std::make_unique<InCoreERI>(eri_);
    if (printme) printf("integrals: one-electron %.2fs, two-electron %.2fs\n", t1 - t0, t2 - t1);
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
    double eold = 0.0;
    converged_ = false;
    if (printme) printf("\n iter          energy            dE        rms(dP)\n");
    for (int iter = 0; iter < param_.maxiter(); ++iter) {
        twoe_->jk(P, J, K);
        const Tensor<double> F = H_ + J - 0.5 * K;
        const double etot = 0.5 * P.trace(H_ + F) + enuc;
        const Tensor<double> Pnew = diagonalize(F);
        const double drms = (Pnew - P).normf() / double(n);
        const double de = etot - eold;
        if (printme) printf("%5d  %18.10f  %12.4e  %12.4e\n", iter, etot, de, drms);
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


namespace {

/// one LCAO orbital sum_mu c_mu chi_mu(r), evaluated pointwise for the MRA projection

/// Holds the molecule and the basis through shared pointers: the function
/// keeps its functor after construction, so references to the caller's
/// objects could dangle.
class LCAOOrbitalFunctor : public FunctionFunctorInterface<double,3> {
public:
    static constexpr std::size_t maxbf = 4096;

    LCAOOrbitalFunctor(std::shared_ptr<const Molecule> molecule, std::shared_ptr<const AtomicBasisSet> aobasis,
                       std::vector<double> c)
        : molecule_(std::move(molecule)), aobasis_(std::move(aobasis)), c_(std::move(c)),
          centers_(molecule_->get_all_coords_vec()) {
        MADNESS_CHECK_THROW(c_.size() <= maxbf, "LCAOOrbitalFunctor: too many basis functions");
    }

    double operator()(const coord_3d& r) const override {
        std::array<double, maxbf> bf;
        aobasis_->eval(*molecule_, r[0], r[1], r[2], bf.data());
        double v = 0.0;
        for (std::size_t mu = 0; mu < c_.size(); ++mu) v += c_[mu] * bf[mu];
        return v;
    }

    std::vector<coord_3d> special_points() const override { return centers_; }

private:
    std::shared_ptr<const Molecule> molecule_;
    std::shared_ptr<const AtomicBasisSet> aobasis_;
    std::vector<double> c_;
    std::vector<coord_3d> centers_;
};

} // namespace


std::vector<Function<double,3>> project_orbitals(World& world, const Molecule& molecule,
                                                  const AtomicBasisSet& aobasis, const Tensor<double>& C,
                                                  const long nmo) {
    MADNESS_CHECK_THROW(C.ndim() == 2 and C.dim(0) == aobasis.nbf(molecule) and C.dim(1) >= nmo,
                        "project_orbitals: coefficients do not match the basis");
    const auto mol = std::make_shared<const Molecule>(molecule);
    const auto basis = std::make_shared<const AtomicBasisSet>(aobasis);
    std::vector<Function<double,3>> mo(nmo);
    for (long i = 0; i < nmo; ++i) {
        std::vector<double> c(C.dim(0));
        for (long mu = 0; mu < C.dim(0); ++mu) c[mu] = C(mu, i);
        const std::shared_ptr<FunctionFunctorInterface<double,3>> f =
            std::make_shared<LCAOOrbitalFunctor>(mol, basis, std::move(c));
        mo[i] = FunctionFactory<double,3>(world).functor(f).truncate_on_project().nofence();
    }
    world.gop.fence();
    truncate(world, mo);
    return mo;
}

} // namespace lcao
} // namespace madness
