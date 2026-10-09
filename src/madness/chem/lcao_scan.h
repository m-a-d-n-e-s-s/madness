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

/// \file lcao_scan.h
/// \brief the LCAO state scan: several SCF starts, the distinct states they reach, and the one to seed moldft

#ifndef MADNESS_CHEM_LCAO_SCAN_H__INCLUDED
#define MADNESS_CHEM_LCAO_SCAN_H__INCLUDED

#include <madness/chem/lcao_scf.h>
#include <madness/chem/lcao_stability.h>
#include <nlohmann/json.hpp>

#include <string>
#include <vector>

namespace madness {
namespace lcao {

/// one distinct SCF solution of the scan
struct LCAOState {
    LCAOSCF::Result result;                 ///< orbitals (occupied first), energies, densities
    std::vector<std::string> found_by;      ///< the starts that reached it
    nlohmann::json populations;             ///< IAO charges and spin populations of the occupied orbitals
    int id = 0;                             ///< discovery order (labels in found_by refer to it)
    std::vector<StabilityRoots> stability;  ///< the lowest Hessian eigenvalues per block (with stability true)
    bool stable = true;                     ///< no eigenvalue below -stability_tol in the block the scan can follow
};

/// the state scan of 33_state_scan_interface.md
///
/// One set of integrals, several starts (scan_starts): the SAD and core densities, and occupation swaps of the
/// states they reach (per spin HOMO->LUMO, HOMO-1->LUMO, HOMO->LUMO+1, reconverged with the swapped occupation
/// kept by maximum overlap). Each start converges to scan_econv/scan_dconv; one that does not is retried with
/// the level shift 1.0 and maxiter 300. Two solutions are the same state if their energies differ by less
/// than 1e-6 Eh and |det C1_occ^T S C2_occ| > 0.99 for each spin. The states are listed by energy.
/// Collective like LCAOSCF: with eri cholesky every rank runs it and holds the same states.
class LCAOStateScan {
public:
    /// @param[in] unrestricted   UHF also for nalpha == nbeta (moldft's spin_restricted false)
    /// @param[in] minbasis       the minimal basis of the IAO populations in the listing
    /// @param[in] charge         the molecule's total charge (for the population printout)
    LCAOStateScan(World& world, const Molecule& molecule, const AtomicBasisSet& aobasis, int nalpha, int nbeta,
                  bool unrestricted, const LCAOParameters& param, const std::string& minbasis, double charge,
                  bool collective);

    /// run the starts and the stability analysis; afterwards the SCF holds the seed state (scf().coefficients(),
    /// ... report it)
    void run();

    /// the Hessian block an unstable state is followed in: RHF->RHF for RHF, UHF->UHF for UHF (RHF->UHF of an RHF
    /// state is reported, not followed: a broken-symmetry state needs spin_restricted false)
    StabilityBlock followed_block() const;

    /// the distinct states, ascending in energy (index 0 the lowest)
    const std::vector<LCAOState>& states() const { return states_; }

    /// the index of the state chosen by the key state (lowest, or an index)
    long seed() const { return seed_; }

    /// the SCF with the integrals; it holds the seed state after run()
    LCAOSCF& scf() { return scf_; }
    const LCAOSCF& scf() const { return scf_; }

    /// the listing, on rank 0
    void print() const;

    nlohmann::json to_json() const;

private:
    World& world_;
    Molecule molecule_;
    AtomicBasisSet aobasis_;
    int nalpha_, nbeta_;
    bool unrestricted_;
    LCAOParameters param_;
    std::string minbasis_;
    double charge_;
    LCAOSCF scf_;
    std::vector<LCAOState> states_;
    std::vector<std::string> failed_;       ///< starts that did not converge, also after the retry
    long seed_ = 0;
    int nextid_ = 0;

    /// converge one start (with its retry) and add what it reaches to the states
    void run_start(const std::string& label, const Tensor<double>& Pa, const Tensor<double>& Pb, bool mom,
                   const Tensor<double>& Ca_occ = Tensor<double>(), const Tensor<double>& Cb_occ = Tensor<double>());

    /// |det C1_occ^T S C2_occ| of the nocc first columns
    double occupied_overlap(const Tensor<double>& C1, const Tensor<double>& C2, long nocc) const;

    /// the lowest Hessian eigenvalues of state s, and whether it is stable
    void analyze(LCAOState& s) const;

    /// rotate state k along its lowest mode of the followed block, line search on the energy, reconverge
    void follow(std::size_t k);

    /// the start flip: the defined broken-symmetry state of flip_atoms
    void flip_start();

    /// the first earlier state of the listing with the same energy and <S^2> as state i (a degenerate
    /// partner with another occupied space, e.g. pi_x against pi_y), or -1
    long degenerate_partner(std::size_t i) const;
};

} // namespace lcao
} // namespace madness

#endif // MADNESS_CHEM_LCAO_SCAN_H__INCLUDED
