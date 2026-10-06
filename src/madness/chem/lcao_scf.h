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
/// \brief closed-shell Hartree-Fock in a Gaussian basis, with integrals from separated kernels

/// Meant as a cheap source of initial orbitals for the MRA calculation, not as
/// an optimized LCAO code: plain Roothaan-Hall iterations, all integrals in
/// memory, replicated on every rank.

#ifndef MADNESS_CHEM_LCAO_SCF_H__INCLUDED
#define MADNESS_CHEM_LCAO_SCF_H__INCLUDED

#include <madness/chem/lcao_integrals.h>
#include <madness/mra/mra.h>
#include <madness/mra/QCCalculationParametersBase.h>
#include <madness/mra/commandlineparser.h>
#include <madness/world/MADworld.h>

#include <memory>
#include <string>
#include <vector>

namespace madness {

/// parameters of the LCAO Hartree-Fock calculation, input group `lcao`

/// The basis set is not among them: it is the `aobasis` of the `dft` group,
/// the basis MADNESS already uses for its initial guess.
class LCAOParameters : public QCCalculationParametersBase {
public:
    static constexpr char const* tag = "lcao";

    LCAOParameters() {
        initialize<double>("kernel_eps", 1.e-4, "relative precision of the Gaussian fit of 1/r: guess grade, "
                           "far below the basis-set error (the v1 reference numbers used 1e-8)");
        initialize<double>("kernel_lo", 1.e-4, "smallest distance at which the fit of 1/r is accurate: guess "
                           "grade (the v1 reference numbers used 1e-6)");
        initialize<double>("kernel_hi", 50.0, "largest distance at which the fit of 1/r is accurate");
        initialize<double>("kernel_screen", 1.e-6, "per primitive quartet, drop the short-range kernel terms "
                           "whose estimated share of the two-electron integral is below this: guess grade "
                           "(0: keep all, as the v1 reference numbers did)");
        initialize<double>("schwarz", 0.0, "skip the two-electron integrals of shell-group quartets whose "
                           "Schwarz bound sqrt((ab|ab)(cd|cd)) is below this (0: compute all)");
        initialize<std::string>("guess", "sad", "starting density: atomic densities from the basis file, "
                                "or the core hamiltonian", {"sad", "core"});
        initialize<int>("maxiter", 100, "maximum number of SCF iterations");
        initialize<double>("econv", 1.e-10, "energy convergence");
        initialize<double>("dconv", 1.e-8, "convergence of the density matrix (rms change per element)");
        initialize<double>("damping", 0.0, "fraction of the previous density mixed into the new one");
        initialize<double>("lindep", 1.e-7, "drop overlap eigenvalues below this (canonical orthogonalization)");
        initialize<int>("print_level", 1, "0: final energy; 1: iterations and energy components; "
                        "2: also the kernel accuracy");
        initialize<bool>("check_mra", false, "madlcao: compare S, T, V and (ii|jj) with MRA quadrature "
                         "of the same basis functions");
        initialize<bool>("check_eri", false, "madlcao: compare the two-electron integrals with the "
                         "unoptimized reference implementation");
        initialize<bool>("seed", false, "madlcao: write the occupied orbitals as <prefix>.restartdata, "
                         "for moldft to start from");
        initialize<int>("seed_rung", 0, "madlcao: rung of the dft group's protocol at which moldft starts "
                        "from the seed; the orbitals are projected at that rung's thresh and k");
    }

    LCAOParameters(World& world, const commandlineparser& parser) : LCAOParameters() {
        read_input_and_commandline_options(world, parser, tag);
    }

    std::string get_tag() const override { return std::string(tag); }

    double kernel_eps() const { return get<double>("kernel_eps"); }
    double kernel_lo() const { return get<double>("kernel_lo"); }
    double kernel_hi() const { return get<double>("kernel_hi"); }
    double kernel_screen() const { return get<double>("kernel_screen"); }
    double schwarz() const { return get<double>("schwarz"); }
    std::string guess() const { return get<std::string>("guess"); }
    int maxiter() const { return get<int>("maxiter"); }
    double econv() const { return get<double>("econv"); }
    double dconv() const { return get<double>("dconv"); }
    double damping() const { return get<double>("damping"); }
    double lindep() const { return get<double>("lindep"); }
    int print_level() const { return get<int>("print_level"); }
    bool check_mra() const { return get<bool>("check_mra"); }
    bool check_eri() const { return get<bool>("check_eri"); }
    bool seed() const { return get<bool>("seed"); }
    int seed_rung() const { return get<int>("seed_rung"); }
};

namespace lcao {

/// the two-electron part of the Fock matrix for a given density matrix

/// The SCF sees two-electron integrals only through this interface, so the
/// in-memory tensor of v1 can later give way to integral-direct, hybrid (MRA
/// Coulomb) or distributed builds without touching the SCF.
class TwoElectronBuilder {
public:
    virtual ~TwoElectronBuilder() = default;

    /// J_{mu nu} = sum_{ls} (mu nu|l s) P_{ls} and K_{mu nu} = sum_{ls} (mu l|nu s) P_{ls}
    virtual void jk(const Tensor<double>& P, Tensor<double>& J, Tensor<double>& K) const = 0;
};

/// J and K from the two-electron integrals held in memory, once per permutational orbit
class InCoreERI : public TwoElectronBuilder {
public:
    explicit InCoreERI(std::shared_ptr<const PackedERI> eri) : eri_(std::move(eri)) {}
    void jk(const Tensor<double>& P, Tensor<double>& J, Tensor<double>& K) const override;

private:
    std::shared_ptr<const PackedERI> eri_;
};

/// closed-shell restricted Hartree-Fock in a Gaussian basis
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

    /// @param[in] nocc  number of doubly occupied orbitals
    LCAOSCF(World& world, const Molecule& molecule, const AtomicBasisSet& aobasis, int nocc,
            const LCAOParameters& param);

    /// compute the integrals and iterate to self-consistency
    /// @return the total energy
    double solve();

    bool converged() const { return converged_; }
    int iterations() const { return iterations_; }
    long nbf() const { return S_.dim(0); }

    const std::vector<Shell>& shells() const { return shells_; }
    const Tensor<double>& overlap() const { return S_; }
    const Tensor<double>& kinetic() const { return T_; }
    const Tensor<double>& nuclear_attraction() const { return V_; }
    const PackedERI& eri() const { return *eri_; }

    /// MO coefficients, one column per orbital, in ascending orbital energy
    const Tensor<double>& coefficients() const { return C_; }
    const Tensor<double>& orbital_energies() const { return eps_; }

    /// total density matrix P = 2 C_occ C_occ^T
    const Tensor<double>& density() const { return P_; }

    const Energies& energies() const { return energies_; }

private:
    World& world_;
    Molecule molecule_;
    AtomicBasisSet aobasis_;
    int nocc_;
    LCAOParameters param_;

    std::vector<Shell> shells_;
    Tensor<double> S_, T_, V_, H_, X_;
    std::shared_ptr<const PackedERI> eri_;
    Tensor<double> C_, eps_, P_;
    std::unique_ptr<TwoElectronBuilder> twoe_;
    Energies energies_;
    bool converged_ = false;
    int iterations_ = 0;

    void compute_integrals();

    /// canonical orthogonalization: X = U s^{-1/2} over the overlap eigenvalues above lindep
    void make_orthogonalizer();

    /// block diagonal density from the atomic guess density matrices of the basis file;
    /// empty if the basis file carries none
    Tensor<double> sad_density() const;

    /// diagonalize F in the orthogonalized basis; sets C_ and eps_ and returns the new density
    Tensor<double> diagonalize(const Tensor<double>& F);
};

/// MRA functions sum_mu C(mu,i) chi_mu(r) for the first nmo columns of C, at the current FunctionDefaults

/// Each basis function is projected once, and the orbitals are their linear
/// combinations (transform). The result is not orthonormalized: projections
/// are orthonormal only to about the projection threshold.
std::vector<Function<double,3>> project_orbitals(World& world, const Molecule& molecule,
                                                  const AtomicBasisSet& aobasis, const Tensor<double>& C,
                                                  long nmo);

} // namespace lcao
} // namespace madness

#endif // MADNESS_CHEM_LCAO_SCF_H__INCLUDED
