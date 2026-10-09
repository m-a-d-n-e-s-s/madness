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
        initialize<std::string>("eri", "cholesky", "two-electron integrals: Cholesky-decomposed at cholesky_tol "
                                "(always without kernel screening and Schwarz skips, positive semidefinite), or "
                                "stored in memory (N^4/8 doubles, on rank 0; the guess-grade kernel_screen and "
                                "schwarz can make them indefinite, and an SCF can then collapse below the HF limit)",
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
        initialize<double>("level_shift", 0.0, "shift (Eh) of the virtual space of each spin, F + s (S - S P S): "
                           "Roothaan-Hall steps without DIIS until max|FPS - SPF| falls below 1e-3, then DIIS "
                           "unshifted; converges SCFs that oscillate under DIIS (0: off)");
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
        initialize<bool>("scan", false, "state scan: several starts, each converged to scan_econv/scan_dconv, the "
                         "distinct states listed by energy; guess lcao seeds moldft from the one chosen by state");
        initialize<std::vector<std::string>>("scan_starts", {"sad", "core", "swaps"}, "starts of the scan: sad, "
                        "core, swaps (per spin HOMO->LUMO, HOMO-1->LUMO, HOMO->LUMO+1 of every state from sad and "
                        "core, kept by maximum overlap)");
        initialize<double>("scan_econv", 1.e-9, "energy convergence of the scan's states (classification grade)");
        initialize<double>("scan_dconv", 1.e-6, "density convergence of the scan's states (classification grade)");
        initialize<int>("scan_max_states", 10, "the scan lists at most this many distinct states");
        initialize<std::string>("state", "lowest", "the scan state that seeds moldft: lowest (the lowest-energy "
                                "stable state), or its index in the listing of a scan with the same molecule and "
                                "settings (0: the lowest energy)");
        initialize<bool>("stability", true, "scan: the lowest eigenvalues of the real orbital Hessian of every state "
                         "(RHF->RHF and RHF->UHF for RHF, UHF->UHF for UHF); an unstable state is followed downhill "
                         "and what it reaches is listed too");
        initialize<int>("stability_roots", 3, "scan: eigenvalues per Hessian block (Davidson)");
        initialize<double>("stability_tol", 1.e-4, "scan: an eigenvalue below -stability_tol (Eh) is an instability");
        initialize<bool>("check_stability", false, "madlcao: the orbital Hessian of the final state against finite "
                         "differences of the energy along a fixed rotation");
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
    double level_shift() const { return get<double>("level_shift"); }
    int diis() const { return get<int>("diis"); }
    double lindep() const { return get<double>("lindep"); }
    int print_level() const { return get<int>("print_level"); }
    bool check_mra() const { return get<bool>("check_mra"); }
    bool check_eri() const { return get<bool>("check_eri"); }
    bool check_onee() const { return get<bool>("check_onee"); }
    bool check_cholesky() const { return get<bool>("check_cholesky"); }
    bool seed() const { return get<bool>("seed"); }
    int seed_rung() const { return get<int>("seed_rung"); }
    std::vector<std::string> population() const { return get<std::vector<std::string>>("population"); }
    bool scan() const { return get<bool>("scan"); }
    std::vector<std::string> scan_starts() const { return get<std::vector<std::string>>("scan_starts"); }
    double scan_econv() const { return get<double>("scan_econv"); }
    double scan_dconv() const { return get<double>("scan_dconv"); }
    int scan_max_states() const { return get<int>("scan_max_states"); }
    std::string state() const { return get<std::string>("state"); }
    bool stability() const { return get<bool>("stability"); }
    int stability_roots() const { return get<int>("stability_roots"); }
    double stability_tol() const { return get<double>("stability_tol"); }
    bool check_stability() const { return get<bool>("check_stability"); }
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

    /// J and K of transition densities D_s = X_s Y_s^T + Y_s X_s^T (X_s, Y_s: N x r): J = J[Da + Db], Ka = K[Da],
    /// Kb = K[Db]. These are symmetric but indefinite, not densities of orbitals (the orbital Hessian of the
    /// stability analysis). An empty Xb is a closed-shell perturbation, Db = Da: J = J[2 Da], Kb left empty.
    /// The default builds Da and Db and calls jk(), exact for a builder that takes any symmetric matrix.
    virtual void jk_transition(const Tensor<double>& Xa, const Tensor<double>& Ya, const Tensor<double>& Xb,
                               const Tensor<double>& Yb, Tensor<double>& J, Tensor<double>& Ka,
                               Tensor<double>& Kb) const {
        const auto density = [](const Tensor<double>& X, const Tensor<double>& Y) {
            const Tensor<double> XY = inner(X, transpose(Y));
            return Tensor<double>(XY + transpose(XY));
        };
        jk(density(Xa, Ya), Xb.size() > 0 ? density(Xb, Yb) : Tensor<double>(), J, Ka, Kb);
    }
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

/// the settings of one LCAO SCF run (LCAOSCF::iterate)
struct SCFOptions {
    double econv = 1.e-6;           ///< energy convergence
    double dconv = 1.e-4;           ///< density convergence, rms change per element
    double damping = 0.0;           ///< fraction of the previous density mixed in
    double level_shift = 0.0;       ///< shifted Roothaan-Hall steps until max|FPS-SPF| < 1e-3, then DIIS (0: off)
    int maxiter = 100;
    int diis = 8;                   ///< DIIS subspace (0: Roothaan-Hall)
    int print_level = 1;            ///< 0: silent; 1: the iteration table and the summary
    bool mom = false;               ///< occupy by maximum overlap with the previous occupied orbitals, not aufbau

    /// the settings of the lcao group, as the plain SCF (solve) uses them
    static SCFOptions from(const LCAOParameters& p) {
        SCFOptions o;
        o.econv = p.econv();
        o.dconv = p.dconv();
        o.damping = p.damping();
        o.level_shift = p.level_shift();
        o.maxiter = p.maxiter();
        o.diis = p.diis();
        o.print_level = p.print_level();
        return o;
    }
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

    /// compute the integrals and iterate to self-consistency from the start of the guess key
    /// @return the total energy of the densities of the last J/K build (the last line of the iteration table)
    double solve();

    /// compute the integrals and the orthogonalizer, once (solve and iterate call it)
    void setup();

    /// the starting spin densities of a start, `sad` (the atomic guess density split by occupation; the
    /// core hamiltonian if the basis file has none) or `core` (the core hamiltonian's orbitals)
    /// @return the start used, sad or core
    std::string start_densities(const std::string& start, Tensor<double>& Pa, Tensor<double>& Pb);

    /// iterate to self-consistency from the spin densities Pa, Pb (Pb ignored for RHF) with the given settings;
    /// with mom, Ca_occ/Cb_occ are the occupied orbitals the occupation follows (and Pa, Pb their densities)
    /// @return as solve
    double iterate(Tensor<double> Pa, Tensor<double> Pb, const SCFOptions& opt,
                   const Tensor<double>& Ca_occ = Tensor<double>(), const Tensor<double>& Cb_occ = Tensor<double>());

    /// UHF also for nalpha == nbeta (the route to broken-symmetry states); call before iterate
    void set_unrestricted(const bool u) { unrestricted_ = u; }

    /// the outcome of the last iterate (or solve): orbitals, densities, energies; to keep and to restore
    struct Result {
        Tensor<double> Ca, Cb, epsa, epsb, Pa, Pb;
        Energies energies;
        double s2 = 0.0;
        bool converged = false;
        int iterations = 0;
    };
    Result result() const {
        return {copy(Ca_), copy(Cb_), copy(epsa_), copy(epsb_), copy(Pa_), copy(Pb_), energies_, s2_, converged_,
                iterations_};
    }
    /// make r the current outcome, as if the last iterate had produced it (the state the accessors report)
    void restore(const Result& r) {
        Ca_ = copy(r.Ca); Cb_ = copy(r.Cb); epsa_ = copy(r.epsa); epsb_ = copy(r.epsb);
        Pa_ = copy(r.Pa); Pb_ = copy(r.Pb);
        energies_ = r.energies; s2_ = r.s2; converged_ = r.converged; iterations_ = r.iterations;
    }

    /// the core hamiltonian T + V
    const Tensor<double>& core_hamiltonian() const { return H_; }

    /// the two-electron part (J and K builds) of this SCF; valid after setup()
    const TwoElectronBuilder& two_electron() const { return *twoe_; }

    /// the energy of spin densities Pa, Pb (Pb ignored for RHF), one J/K build; collective like iterate
    double energy(const Tensor<double>& Pa, const Tensor<double>& Pb) const;

    bool converged() const { return converged_; }
    int iterations() const { return iterations_; }
    long nbf() const { return S_.dim(0); }

    /// true for closed-shell RHF, false for UHF
    bool restricted() const { return nalpha_ == nbeta_ and not unrestricted_; }
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
    bool unrestricted_ = false;     ///< UHF also for nalpha == nbeta
    bool setup_ = false;

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
