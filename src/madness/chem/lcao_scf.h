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

/// \file lcao_scf.h
/// \brief Hartree-Fock in a Gaussian basis, with integrals from separated kernels

/// Meant as a cheap source of initial orbitals for the MRA calculation, not as
/// an optimized LCAO code: closed-shell RHF or UHF with DIIS, the two-electron
/// integrals in memory or as Cholesky vectors (key eri), all on the rank that
/// runs it (moldft's `guess lcao` runs it on rank 0 alone; madlcao on every rank).

#ifndef MADNESS_CHEM_LCAO_SCF_H__INCLUDED
#define MADNESS_CHEM_LCAO_SCF_H__INCLUDED

#include <madness/chem/lcao_integrals.h>
#include <madness/mra/mra.h>
#include <madness/mra/QCCalculationParametersBase.h>
#include <madness/mra/commandlineparser.h>
#include <madness/world/MADworld.h>
#include <nlohmann/json.hpp>

#include <algorithm>
#include <memory>
#include <string>
#include <vector>

namespace madness {

/// parameters of the LCAO Hartree-Fock calculation, input group `lcao`

/// Read by madlcao, and by moldft for `guess lcao` (dft group). The basis set
/// is its own key, not the `aobasis` of the `dft` group, which stays the basis
/// of the atomic guess.
class LCAOParameters : public QCCalculationParametersBase {
public:
    static constexpr char const* tag = "lcao";

    LCAOParameters() {
        initialize<std::string>("basis", "6-31gss", "Gaussian basis set, a basis file of the MADNESS data "
                                "directory (6-31gss, aug-cc-pvdz, 6-31g, ...)");
        initialize<double>("kernel_eps", 1.e-3, "relative precision of the Gaussian fit of 1/r: guess grade, "
                           "invisible to the MRA iterations (the v1 reference numbers used 1e-8)");
        initialize<double>("kernel_lo", 1.e-4, "smallest distance at which the fit of 1/r is accurate: guess "
                           "grade (the v1 reference numbers used 1e-6)");
        initialize<double>("kernel_hi", 50.0, "largest distance at which the fit of 1/r is accurate; raised "
                           "to the largest interatomic distance + 10 bohr if that is longer");
        initialize<double>("kernel_screen", 1.e-5, "per primitive quartet, drop the short-range kernel terms "
                           "whose estimated share of the two-electron integral is below this: guess grade "
                           "(0: keep all, as the v1 reference numbers did)");
        initialize<double>("schwarz", 1.e-6, "skip the two-electron integrals of shell-group quartets whose "
                           "Schwarz bound sqrt((ab|ab)(cd|cd)) is below this: guess grade (0: compute all)");
        initialize<std::string>("eri", "incore", "two-electron integrals: stored in memory, or Cholesky-"
                                "decomposed at cholesky_tol (always without kernel screening and Schwarz skips)",
                                {"incore", "cholesky"});
        initialize<double>("cholesky_tol", 1.e-6, "Cholesky decomposition of the two-electron integrals: the "
                           "largest residual (mu nu|mu nu) left, which bounds the error of every integral");
        initialize<std::string>("guess", "sad", "starting density: atomic densities from the basis file, "
                                "or the core hamiltonian", {"sad", "core"});
        initialize<int>("maxiter", 100, "maximum number of SCF iterations");
        initialize<double>("econv", 1.e-6, "energy convergence: guess grade, invisible to the MRA iterations "
                           "(the v1 reference numbers used 1e-10)");
        initialize<double>("dconv", 1.e-4, "convergence of the density matrix (rms change per element): guess "
                           "grade (the v1 reference numbers used 1e-8)");
        initialize<double>("damping", 0.0, "fraction of the previous density mixed into the new one");
        initialize<int>("diis", 8, "DIIS subspace: Pulay extrapolation of the Fock matrix from this many "
                        "iterations (0: plain Roothaan-Hall)");
        initialize<double>("lindep", 1.e-7, "drop overlap eigenvalues below this (canonical orthogonalization)");
        initialize<int>("print_level", 1, "0: final energy; 1: iterations and energy components; "
                        "2: also the kernel accuracy");
        initialize<bool>("check_mra", false, "madlcao: compare S, T, V and (ii|jj) with MRA quadrature "
                         "of the same basis functions");
        initialize<bool>("check_eri", false, "madlcao: compare the two-electron integrals with the "
                         "unoptimized reference implementation");
        initialize<bool>("check_onee", false, "madlcao: compare the screened nuclear attraction with every "
                         "term computed");
        initialize<bool>("check_cholesky", false, "madlcao: decompose the two-electron integrals at cholesky_tol "
                         "and compare L L^T with the stored integrals (needs kernel_screen 0; schwarz 0)");
        initialize<bool>("seed", false, "madlcao: write the occupied orbitals as <prefix>.restartdata, "
                         "for moldft to start from");
        initialize<int>("seed_rung", 0, "madlcao: rung of the dft group's protocol at which moldft starts "
                        "from the seed; the orbitals are projected at that rung's thresh and k");
        initialize<std::vector<std::string>>("population", {"none"}, "madlcao: atomic charges of the LCAO "
                        "orbitals: none, mulliken, lowdin, iao; the basis sets are the dft keys population_basis "
                        "and population_minbasis (moldft's guess lcao uses the dft key population)");
    }

    LCAOParameters(World& world, const commandlineparser& parser) : LCAOParameters() {
        read_input_and_commandline_options(world, parser, tag);
    }

    std::string get_tag() const override { return std::string(tag); }

    std::string basis() const { return get<std::string>("basis"); }
    double kernel_eps() const { return get<double>("kernel_eps"); }
    double kernel_lo() const { return get<double>("kernel_lo"); }
    double kernel_hi() const { return get<double>("kernel_hi"); }
    double kernel_screen() const { return get<double>("kernel_screen"); }
    double schwarz() const { return get<double>("schwarz"); }
    std::string eri() const { return get<std::string>("eri"); }
    double cholesky_tol() const { return get<double>("cholesky_tol"); }
    std::string guess() const { return get<std::string>("guess"); }
    int maxiter() const { return get<int>("maxiter"); }
    double econv() const { return get<double>("econv"); }
    double dconv() const { return get<double>("dconv"); }
    double damping() const { return get<double>("damping"); }
    int diis() const { return get<int>("diis"); }
    double lindep() const { return get<double>("lindep"); }
    int print_level() const { return get<int>("print_level"); }
    bool check_mra() const { return get<bool>("check_mra"); }
    bool check_eri() const { return get<bool>("check_eri"); }
    bool check_onee() const { return get<bool>("check_onee"); }
    bool check_cholesky() const { return get<bool>("check_cholesky"); }
    bool seed() const { return get<bool>("seed"); }
    int seed_rung() const { return get<int>("seed_rung"); }
    std::vector<std::string> population() const {
        std::vector<std::string> p = get<std::vector<std::string>>("population");
        for (auto& s : p) s.erase(std::remove(s.begin(), s.end(), '"'), s.end());
        return p;
    }
};

namespace lcao {

/// the two-electron part of the Fock matrix for a given density matrix

/// The SCF sees two-electron integrals only through this interface, so the
/// in-memory tensor of v1 can later give way to integral-direct, hybrid (MRA
/// Coulomb) or distributed builds without touching the SCF.
class TwoElectronBuilder {
public:
    virtual ~TwoElectronBuilder() = default;

    /// J of the total density and K of each spin density,
    ///   J_{mu nu} = sum_{ls} (mu nu|l s) (Pa+Pb)_{ls},   Ks_{mu nu} = sum_{ls} (mu l|nu s) Ps_{ls}
    /// An empty Pb stands for a closed shell, Pb = Pa; Kb is then left empty.
    virtual void jk(const Tensor<double>& Pa, const Tensor<double>& Pb, Tensor<double>& J, Tensor<double>& Ka,
                    Tensor<double>& Kb) const = 0;
};

/// J and K from the two-electron integrals held in memory, once per permutational orbit

/// The pass over the integrals is split into tasks on the thread pool of this process.
class InCoreERI : public TwoElectronBuilder {
public:
    InCoreERI(World& world, std::shared_ptr<const PackedERI> eri) : world_(world), eri_(std::move(eri)) {}
    void jk(const Tensor<double>& Pa, const Tensor<double>& Pb, Tensor<double>& J, Tensor<double>& Ka,
            Tensor<double>& Kb) const override;

private:
    World& world_;
    std::shared_ptr<const PackedERI> eri_;
};

/// Hartree-Fock in a Gaussian basis: RHF for a closed shell, UHF otherwise

/// Spin densities are Ps = Cs_occ Cs_occ^T, one electron per spin orbital; a
/// closed shell keeps Pb = Pa, and its total density is 2 C_occ C_occ^T.
class LCAOSCF {
public:
    /// the energy and its parts, labelled as moldft prints them
    struct Energies {
        double kinetic = 0.0;
        double nuclear_attraction = 0.0;
        double coulomb = 0.0;
        double exchange = 0.0;
        double nuclear_repulsion = 0.0;
        double total = 0.0;
    };

    /// @param[in] nalpha      number of alpha electrons
    /// @param[in] nbeta       number of beta electrons; nalpha == nbeta gives closed-shell RHF
    /// @param[in] collective  every rank of world constructs it and calls solve(): with eri cholesky the
    ///                        integrals are spread over the ranks, rank 0 runs DIIS and the
    ///                        diagonalization and broadcasts the orbitals. Otherwise this rank does all.
    LCAOSCF(World& world, const Molecule& molecule, const AtomicBasisSet& aobasis, int nalpha, int nbeta,
            const LCAOParameters& param, bool collective = false);

    /// compute the integrals and iterate to self-consistency
    /// @return the total energy of the densities of the last J/K build (the last line of the iteration table)
    double solve();

    bool converged() const { return converged_; }
    int iterations() const { return iterations_; }
    long nbf() const { return S_.dim(0); }

    /// true for closed-shell RHF, false for UHF
    bool restricted() const { return nalpha_ == nbeta_; }
    int nalpha() const { return nalpha_; }
    int nbeta() const { return nbeta_; }

    const std::vector<Shell>& shells() const { return shells_; }

    /// the fit of 1/r the integrals use: kernel_hi raised to cover the molecule (set by solve)
    const GaussianKernel& kernel() const { return kernel_; }
    const Tensor<double>& overlap() const { return S_; }
    const Tensor<double>& kinetic() const { return T_; }
    const Tensor<double>& nuclear_attraction() const { return V_; }

    /// the two-electron integrals in memory; only with eri incore
    const PackedERI& eri() const {
        MADNESS_CHECK_THROW(eri_, "LCAOSCF: the two-electron integrals are not in memory (eri cholesky)");
        return *eri_;
    }

    /// MO coefficients of spin 0 (alpha) or 1 (beta), one column per orbital, in ascending orbital energy
    const Tensor<double>& coefficients(const int spin = 0) const { return spin == 0 ? Ca_ : Cb_; }
    const Tensor<double>& orbital_energies(const int spin = 0) const { return spin == 0 ? epsa_ : epsb_; }

    /// total density matrix Pa + Pb
    Tensor<double> density() const { return Pa_ + Pb_; }

    /// density matrix of spin 0 (alpha) or 1 (beta); for a closed shell both are C_occ C_occ^T
    const Tensor<double>& density(const int spin) const { return spin == 0 ? Pa_ : Pb_; }

    /// <S^2> of the determinant: 0 for RHF, S(S+1) plus the spin contamination for UHF
    double s2() const { return s2_; }

    /// the energy and its parts for the densities of the last J/K build; the final densities (density())
    /// come from the last diagonalization and differ from them by the last step (rms below dconv)
    const Energies& energies() const { return energies_; }

private:
    World& world_;
    Molecule molecule_;
    AtomicBasisSet aobasis_;
    int nalpha_, nbeta_;
    LCAOParameters param_;
    bool collective_ = false;       ///< every rank of world_ takes part

    std::vector<Shell> shells_;
    GaussianKernel kernel_;
    Tensor<double> S_, T_, V_, H_, X_;
    std::shared_ptr<const PackedERI> eri_;
    Tensor<double> Ca_, Cb_, epsa_, epsb_, Pa_, Pb_;
    std::unique_ptr<TwoElectronBuilder> twoe_;
    Energies energies_;
    double s2_ = 0.0;
    bool converged_ = false;
    int iterations_ = 0;

    void compute_integrals();

    /// canonical orthogonalization: X = U s^{-1/2} over the overlap eigenvalues above lindep
    void make_orthogonalizer();

    /// block diagonal density from the atomic guess density matrices of the basis file;
    /// empty if the basis file carries none
    Tensor<double> sad_density() const;

    /// diagonalize F in the orthogonalized basis; sets C and eps and returns the density
    /// of the nocc lowest orbitals, C_occ C_occ^T. If collective, rank 0 diagonalizes and
    /// broadcasts C and eps, so every rank has the same orbitals.
    Tensor<double> diagonalize(const Tensor<double>& F, int nocc, Tensor<double>& C, Tensor<double>& eps) const;

    /// the commutator F P S - S P F in the orthogonal basis, which vanishes at self-consistency
    Tensor<double> commutator_error(const Tensor<double>& F, const Tensor<double>& P) const;
};

/// MRA functions sum_mu C(mu,i) chi_mu(r) for the first nmo columns of C, at the current FunctionDefaults

/// Each basis function is projected once, and the orbitals are their linear
/// combinations (transform). The result is not orthonormalized: projections
/// are orthonormal only to about the projection threshold.
std::vector<Function<double,3>> project_orbitals(World& world, const Molecule& molecule,
                                                  const AtomicBasisSet& aobasis, const Tensor<double>& C,
                                                  long nmo);

/// atomic charges and spin populations of LCAO orbitals by the schemes of population.h

/// The scheme definitions of SCF::population_analysis for MRA orbitals, with analytic overlaps:
/// S over the population basis and X = S(population basis, aobasis) C. The occupied orbitals
/// are taken orthonormal in aobasis's overlap. Runs on the calling rank alone.
/// @param[in] Ca, Cb       occupied orbital coefficients of each spin (columns); Cb empty without beta electrons
/// @param[in] schemes      mulliken, lowdin, iao
/// @param[in] basis_proj   basis set for mulliken and lowdin
/// @param[in] basis_min    minimal basis set for iao
/// @param[in] print_level  0: no output; otherwise the tables of AtomicPopulations::print
/// @return  the JSON layout of SCF::population_analysis
nlohmann::json population_analysis(const Molecule& molecule, const AtomicBasisSet& aobasis, const Tensor<double>& Ca,
                                   const Tensor<double>& Cb, const std::vector<std::string>& schemes,
                                   const std::string& basis_proj, const std::string& basis_min,
                                   double total_charge, int print_level);

} // namespace lcao
} // namespace madness

#endif // MADNESS_CHEM_LCAO_SCF_H__INCLUDED
