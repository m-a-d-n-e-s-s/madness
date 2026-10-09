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

/// \file lcao_stability.h
/// \brief the stability of an LCAO SCF solution: the lowest eigenvalues of its real orbital Hessian

#ifndef MADNESS_CHEM_LCAO_STABILITY_H__INCLUDED
#define MADNESS_CHEM_LCAO_STABILITY_H__INCLUDED

#include <madness/chem/lcao_scf.h>

#include <string>
#include <vector>

namespace madness {
namespace lcao {

/// the blocks of the real orbital Hessian: RHF->RHF (singlet rotations of an RHF solution), RHF->UHF (triplet
/// rotations: alpha and beta opposite, the instability towards UHF), UHF->UHF (independent alpha and beta rotations)
enum class StabilityBlock { rhf_rhf, rhf_uhf, uhf_uhf };

std::string to_string(StabilityBlock b);

/// the lowest eigenpairs of one block; the eigenvectors as rotation parameters x_ai (virtual x occupied, per spin)
struct StabilityRoots {
    StabilityBlock block = StabilityBlock::uhf_uhf;
    std::vector<double> eigenvalues;            ///< ascending
    std::vector<Tensor<double>> xa, xb;         ///< per root; rhf_rhf and rhf_uhf: xb empty
    bool converged = false;
    int products = 0;                           ///< Hessian-vector products (one J/K build each)
};

/// the real orbital Hessian A + B of an SCF solution (its orbitals with the occupied first, and their energies)
///
/// (A+B) x = (eps_a - eps_i) x_ai + [C_v^T (J[Da + Db] - K[D_s]) C_o]_ai per spin, with the transition densities
/// D_s = X_s C_o^T + C_o X_s^T, X_s = C_v x_s (TwoElectronBuilder::jk_transition). RHF->RHF: Db = Da;
/// RHF->UHF: Db = -Da, so J drops out. The eigenvalues are those of A + B (Eh), the stability matrix of Seeger
/// and Pople, for every block. The second derivative of the energy along a unit rotation is 4 (A+B) for RHF->RHF
/// and 2 (A+B) for UHF->UHF, and PySCF 2.14 reports exactly these (rhf_internal, uhf_internal); its rhf_external
/// reports A + B. All three agree with PySCF to the printed digits (34_). Collective like the J/K builds.
class OrbitalHessian {
public:
    OrbitalHessian(const LCAOSCF& scf, const LCAOSCF::Result& state, int nalpha, int nbeta, StabilityBlock block);

    long dim() const { return na_ * nva_ + nb_ * nvb_; }
    StabilityBlock block() const { return block_; }

    /// (A+B) x, x of length dim(): the alpha parameters (nva x na, row-major), then the beta ones
    Tensor<double> product(const Tensor<double>& x) const;

    /// eps_a - eps_i, the diagonal estimate (preconditioner, start vectors)
    const Tensor<double>& diagonal() const { return diag_; }

    /// x as the rotation matrices of each spin (nv x no); xb empty for the RHF blocks
    void split(const Tensor<double>& x, Tensor<double>& xa, Tensor<double>& xb) const;

private:
    const LCAOSCF& scf_;
    StabilityBlock block_;
    long na_, nb_, nva_, nvb_;
    Tensor<double> coa_, cva_, cob_, cvb_;      ///< occupied and virtual orbitals of each spin
    Tensor<double> diag_;
};

/// the nroots lowest eigenpairs of H by Davidson, until every residual norm is below tol (or maxiter)
StabilityRoots lowest_roots(const OrbitalHessian& H, int nroots, double tol = 1.e-5, int maxiter = 60);

/// the nocc occupied orbitals of C (occupied first) rotated by the real parameters t x (x: nv x nocc):
/// C_o' = C_o (V cos(t s) V^T + 1 - V V^T) + C_v U sin(t s) V^T with x = U s V^T, orthonormal for any t
Tensor<double> rotate_occupied(const Tensor<double>& C, long nocc, const Tensor<double>& x, double t);

} // namespace lcao
} // namespace madness

#endif // MADNESS_CHEM_LCAO_STABILITY_H__INCLUDED
