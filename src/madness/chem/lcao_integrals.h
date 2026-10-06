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

/// \file lcao_integrals.h
/// \brief integrals over Cartesian Gaussians from separated (Gaussian-sum) kernels

/// MADNESS writes 1/r as a sum of Gaussians, 1/r = sum_m w_m exp(-t_m r^2) (GFit),
/// to apply the Coulomb operator to multiwavelets. Applied to Gaussian basis
/// functions instead, every term factorizes into x, y and z, so all integrals of
/// an LCAO Hartree-Fock calculation reduce to 1D integrals of polynomial times
/// Gaussian. Those are evaluated exactly with Gauss-Hermite quadrature; the only
/// approximation is the fit of the kernel. No Boys function is involved.

#ifndef MADNESS_CHEM_LCAO_INTEGRALS_H__INCLUDED
#define MADNESS_CHEM_LCAO_INTEGRALS_H__INCLUDED

#include <madness/tensor/tensor.h>
#include <madness/chem/molecule.h>
#include <madness/chem/molecularbasis.h>

#include <array>
#include <vector>

namespace madness {

class World;

namespace lcao {

/// one contracted Cartesian Gaussian shell on a center, as AtomicBasisSet defines it

/// The basis functions of the shell are
///   chi(r) = sum_i coeff[i] (x-X)^lx (y-Y)^ly (z-Z)^lz exp(-expnt[i] |r-R|^2)
/// with lx+ly+lz = l, exactly as ContractedGaussianShell::eval evaluates them:
/// the primitive normalization is part of coeff, and the Cartesian components
/// are not normalized individually.
struct Shell {
    std::array<double,3> center = {0.0, 0.0, 0.0};  ///< position in bohr
    int atom = 0;                   ///< index of the atom in the Molecule
    int l = 0;                      ///< angular momentum
    std::vector<double> expnt;      ///< primitive exponents
    std::vector<double> coeff;      ///< contraction coefficients, primitive normalization included
    int offset = 0;                 ///< index of the first basis function of this shell

    int ncart() const { return (l + 1) * (l + 2) / 2; }
};

/// Cartesian exponents (lx,ly,lz) of the components of a shell, in ContractedGaussianShell::eval order
std::vector<std::array<int,3>> cartesian_components(int l);

/// the shells of a molecule, with the basis functions in AtomicBasisSet order
std::vector<Shell> make_shells(const Molecule& molecule, const AtomicBasisSet& aobasis);

/// Gauss-Hermite rules: int exp(-y^2) f(y) dy = sum_k w_k f(y_k), exact for polynomials of degree 2n-1
class GaussHermiteRule {
public:
    /// nodes and weights for n = 1..nmax, from the Golub-Welsch eigenvalue problem
    explicit GaussHermiteRule(int nmax = 16);

    int nmax() const { return int(nodes_.size()) - 1; }
    const std::vector<double>& nodes(int n) const;
    const std::vector<double>& weights(int n) const;

private:
    std::vector<std::vector<double>> nodes_, weights_;
};

/// a radial kernel written as a sum of Gaussians, f(r) = sum_m w[m] exp(-t[m] r^2)
struct GaussianKernel {
    std::vector<double> w, t;

    /// the fit of 1/r that MADNESS uses for its Coulomb operator (GFit::CoulombFit)

    /// @param[in] lo  smallest distance at which the fit must be accurate
    /// @param[in] hi  largest distance at which the fit must be accurate
    /// @param[in] eps relative precision of the fit on [lo, hi]
    static GaussianKernel coulomb(double lo, double hi, double eps);

    std::size_t size() const { return w.size(); }

    /// value of the kernel at distance r
    double operator()(double r) const;

    /// max over [lo, hi] of |f(r) - 1/r| r, sampled on a logarithmic grid
    double max_relative_coulomb_error(double lo, double hi, int npt = 2000) const;
};

/// two-electron integrals (ij|kl), stored once per permutational orbit

/// (ij|kl) = (ji|kl) = (ij|lk) = (kl|ij) = ..., so only i >= j, k >= l and
/// ij >= kl are stored, with pairs numbered ij = i(i+1)/2 + j: about nbf^4/8
/// values instead of nbf^4. The accessors take the four indices in any order
/// of the orbit.
class PackedERI {
public:
    PackedERI() = default;
    explicit PackedERI(long nbf);

    long nbf() const { return nbf_; }
    std::size_t size() const { return data_.size(); }

    /// number of the pair (i,j), for either order
    static std::size_t pair(const long i, const long j) {
        return (i >= j) ? std::size_t(i) * (i + 1) / 2 + j : std::size_t(j) * (j + 1) / 2 + i;
    }

    /// position of (ij|kl)
    static std::size_t index(const long i, const long j, const long k, const long l) {
        const std::size_t p = pair(i, j), q = pair(k, l);
        return (p >= q) ? p * (p + 1) / 2 + q : q * (q + 1) / 2 + p;
    }

    double operator()(const long i, const long j, const long k, const long l) const {
        return data_[index(i, j, k, l)];
    }
    double& operator()(const long i, const long j, const long k, const long l) { return data_[index(i, j, k, l)]; }

    /// the stored values in order of position: pairs ij ascending, and kl = 0..ij within each
    const double* data() const { return data_.data(); }

private:
    long nbf_ = 0;
    std::vector<double> data_;
};

/// how many group quartets SeparatedGaussianIntegrals::eri computed and how many Schwarz screening skipped
struct ERIStats {
    long computed = 0;
    long skipped = 0;
};

/// one- and two-electron integrals over the shells, with point nuclei and a Gaussian-sum Coulomb kernel
class SeparatedGaussianIntegrals {
public:
    SeparatedGaussianIntegrals(const std::vector<Shell>& shells, const GaussianKernel& coulomb);

    long nbf() const { return nbf_; }

    Tensor<double> overlap() const;
    Tensor<double> kinetic() const;

    /// attraction to point nuclei with charges Atom::q
    Tensor<double> nuclear_attraction(const Molecule& molecule) const;

    /// all two-electron integrals (mu nu|lambda sigma), chemists' notation, once per permutational orbit

    /// Computed as tasks on the thread pool of this process; every rank computes all of them.
    /// Shells that differ only in their contraction coefficients form groups, and the
    /// group quartets are the unit of work and of screening.
    /// @param[in] screen   per primitive quartet, drop the short-range kernel terms whose estimated
    ///                     share of the integral is below this (0: keep all terms)
    /// @param[in] schwarz  skip the group quartets whose Schwarz bound sqrt((ab|ab)(cd|cd)) is below
    ///                     this; their integrals stay zero (0: compute all)
    /// @param[out] stats   if given, how many group quartets were computed and skipped
    PackedERI eri(World& world, double screen = 0.0, double schwarz = 0.0, ERIStats* stats = nullptr) const;

    /// eri() as v1 computed it: every primitive quartet, kernel term and axis on its own.
    /// Slow; kept as the reference that eri() is checked against (lcao group: check_eri).
    Tensor<double> eri_reference() const;

private:
    std::vector<Shell> shells_;
    GaussianKernel coulomb_;
    GaussHermiteRule gh_;
    long nbf_ = 0;
};

} // namespace lcao
} // namespace madness

#endif // MADNESS_CHEM_LCAO_INTEGRALS_H__INCLUDED
