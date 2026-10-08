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

/// \file population.h
/// \brief atomic populations of occupied orbitals: Mulliken, Loewdin and intrinsic atomic orbitals (IAO)

/// Everything here works on overlap matrices alone, so the same functions serve
/// MRA orbitals (overlaps by quadrature) and Gaussian-basis orbitals (analytic
/// overlaps). The inputs, for the occupied orbitals phi_i of one spin:
///   - phi_i orthonormal, <phi_i|phi_j> = delta_ij, with occupations n_i
///   - a basis chi_mu centred on the atoms; S = <chi|chi>, X = <chi|phi>
///   - the atom each basis function sits on
/// S and X must refer to the same functions chi (e.g. both from the same
/// projection into MRA).

#ifndef MADNESS_CHEM_POPULATION_H__INCLUDED
#define MADNESS_CHEM_POPULATION_H__INCLUDED

#include <madness/tensor/tensor.h>
#include <nlohmann/json.hpp>

#include <string>
#include <vector>

namespace madness {
namespace population {

/// electrons per atom from one spin's occupied orbitals, with the scheme's diagnostics
struct SpinPopulation {
    Tensor<double> electrons;         ///< electrons on each atom
    Tensor<double> completeness;      ///< mulliken/lowdin: c_i = X_i^T S^-1 X_i, the share of phi_i inside the basis
    double spilling = 0.0;            ///< mulliken/lowdin: 1 - sum_i n_i c_i / sum_i n_i (Sanchez-Portal et al. 1995)
    double min_w_eigenvalue = 0.0;    ///< smallest eigenvalue of W = X^T S^-1 X (all schemes)
    double max_norm_error = 0.0;      ///< iao: max_i |sum_rho Q_rho,i^2 - 1|, zero when phi lies in the IAO span
};

/// S^p of a symmetric positive semidefinite matrix; eigenvalues below tol * (largest) are dropped
/// (for p < 0: the pseudo-inverse power on the retained space)
Tensor<double> matrix_power(const Tensor<double>& S, double p, double tol = 1.e-12);

/// Mulliken or Loewdin populations of the occupied space projected onto the basis

/// The fit of phi_i onto the basis is C_i = S^-1 X_i; the fitted orbitals overlap as
/// W = C^T S C = X^T S^-1 X, whose diagonal c_i <= 1 is each orbital's completeness. The
/// fits are orthonormalized among themselves, Cbar = C W^-1/2 (the depolarized orbitals of
/// the IAO construction), so the electrons sum to sum_i n_i and what the basis misses shows
/// as the spilling instead:
///   mulliken: N_mu = sum_i n_i (Cbar_i o S Cbar_i)_mu
///   lowdin:   N_mu = sum_i n_i ((S^1/2 Cbar_i)_mu)^2
/// Both are invariant under rotations among equally occupied orbitals. (Renormalizing each
/// fit by its own c_i also gives the sum rule, but not that invariance.)
/// @param[in]  scheme  "mulliken" or "lowdin"
SpinPopulation projection_population(const std::string& scheme, const Tensor<double>& S, const Tensor<double>& X,
                                     const Tensor<double>& occ, const std::vector<int>& atom_of_bf, long natom);

/// populations of the intrinsic atomic orbitals (Knizia, J. Chem. Theory Comput. 9, 4834 (2013))

/// For a complete large basis (MRA) the IAOs follow from S and X of a minimal basis alone:
///   W = X^T S^-1 X,  D = S^-1 X W^-1/2  (the depolarized orbitals chi D)
///   M = 1 - D D^T S,  K = X^T (2 D D^T S - 1)
///   A = chi M + phi K  (the IAOs),  <A|A> = M^T S M + M^T X K + K^T X^T M + K^T K,  <A|phi> = M^T X + K^T
///   Q = <A|A>^-1/2 <A|phi>,  N_rho = sum_i n_i Q_rho,i^2
/// The orbitals lie in the IAO span, so every column of Q has unit norm and the electrons
/// sum to sum_i n_i to the precision of the overlaps. Throws when W is singular: the minimal
/// basis then misses an orbital's atomic character.
SpinPopulation iao_population(const Tensor<double>& S, const Tensor<double>& X, const Tensor<double>& occ,
                              const std::vector<int>& atom_of_bf, long natom);

/// dispatch: "mulliken" and "lowdin" to projection_population, "iao" to iao_population
SpinPopulation spin_population(const std::string& scheme, const Tensor<double>& S, const Tensor<double>& X,
                               const Tensor<double>& occ, const std::vector<int>& atom_of_bf, long natom);

/// the schemes this file implements
const std::vector<std::string>& schemes();

/// one scheme's atomic charges and spin populations for a molecule
struct AtomicPopulations {
    std::string scheme;                 ///< mulliken, lowdin or iao
    std::string basis;                  ///< the projection or minimal basis set
    long nbf = 0;                       ///< its number of functions
    long nalpha = 0, nbeta = 0;         ///< electrons of each spin
    std::vector<std::string> symbols;   ///< element of each atom
    std::vector<double> charge;         ///< q_A = Z_A - N_A(alpha) - N_A(beta)
    std::vector<double> spin;           ///< N_A(alpha) - N_A(beta)
    SpinPopulation alpha, beta;         ///< the per-spin results with their diagnostics

    bool open() const { return nalpha != nbeta; }

    /// the table (charge, and spin for open shells, per atom), the sums and the scheme's diagnostics
    void print(double total_charge, int print_level) const;

    nlohmann::json to_json() const;
};

/// charges and spin populations from the populations of both spins (beta empty when nbeta is 0)
AtomicPopulations atomic_populations(const std::string& scheme, const std::string& basis, long nbf,
                                     const std::vector<std::string>& symbols, const std::vector<double>& Z,
                                     long nalpha, long nbeta, const SpinPopulation& alpha,
                                     const SpinPopulation& beta);

} // namespace population
} // namespace madness

#endif // MADNESS_CHEM_POPULATION_H__INCLUDED
