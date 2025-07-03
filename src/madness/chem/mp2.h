/*
  This file is part of MADNESS.

  Copyright (C) 2007,2010 Oak Ridge National Laboratory

  This program is free software; you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation; either version 2 of the License, or
  (at your option) any later version.

  This program is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
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

/*!
  \file mp2.h
  \brief Solves molecular MP2 equations
  \defgroup Solves molecular MP2 equations
  \ingroup examples
*/

#ifndef MP2_H_
#define MP2_H_

#include <madness/mra/mra.h>
#include <madness/mra/lbdeux.h>
#include"madness/mra/QCCalculationParametersBase.h"
#include<madness/chem/SCF.h>
#include <madness/mra/nonlinsol.h>
#include<madness/chem/projector.h>
#include<madness/chem/correlationfactor.h>
#include<madness/chem/electronic_correlation_factor.h>
#include<madness/chem/MolecularOrbitals.h>
#include<madness/chem/nemo.h>

#include <madness/world/text_fstream_archive.h>

using madness::archive::TextFstreamInputArchive;
using madness::archive::TextFstreamOutputArchive;


#include <iostream>

namespace madness {

struct LBCost {
    double leaf_value;
    double parent_value;

    LBCost(double leaf_value = 1.0, double parent_value = 1.0)
            : leaf_value(leaf_value), parent_value(parent_value) {}

    double operator()(const Key<6>& key, const FunctionNode<double, 6>& node) const {
        return node.coeff().size();
    }
};


class HartreeFock {
    World& world;
public:
    std::shared_ptr<Nemo> nemo_ptr;
private:
    mutable double coords_sum;     // sum of square of coords at last solved geometry

    // save the Coulomb potential
    mutable real_function_3d coulomb;

    /// reconstructed orbitals: R * phi, where R is the nuclear correlation factor
    std::vector<real_function_3d> orbitals_;

    /// R^2 * phi, where R is the nuclear correlation factor, corresponds
    /// to the bra space of the transformed operators
    std::vector<real_function_3d> R2orbitals_;

public:

    HartreeFock(World& world, std::shared_ptr<Nemo> nemo) :
            world(world), nemo_ptr(nemo), coords_sum(-1.0) {
    }

    World& get_world() {return world;} // just to silence compiler about unused private variable

    bool provides_gradient() const { return true; }

    double value() {
        return value(nemo_ptr->get_calc()->molecule.get_all_coords());
    }

    double value(const Tensor<double>& x) {

        // fast return if the reference is already solved at this geometry
        double xsq = x.sumsq();
        if (xsq == coords_sum) return nemo_ptr->get_calc()->current_energy;

        nemo_ptr->get_calc()->molecule.set_all_coords(x.reshape(nemo_ptr->get_calc()->molecule.natom(), 3));
        coords_sum = xsq;
        nemo_ptr->value(x);

        MolecularOrbitals<double,3> mos(nemo_ptr->get_calc()->amo,nemo_ptr->get_calc()->aeps);
        reset_orbitals(mos);
        return nemo_ptr->get_calc()->current_energy;
    }

    void reset_orbitals(const MolecularOrbitals<double,3>& mos) {
        nemo_ptr->get_calc()->amo=mos.get_mos();
        nemo_ptr->get_calc()->aeps=mos.get_eps();
        MADNESS_CHECK(size_t(nemo_ptr->get_calc()->aeps.size())==nemo_ptr->get_calc()->amo.size());
        orbitals_ = nemo_ptr->R*nemo_ptr->get_calc()->amo;
        R2orbitals_ = nemo_ptr->ncf->square()*nemo_ptr->get_calc()->amo;
    }

    Tensor<double> gradient(const Tensor<double>& x) {

        value(x); // Ensures DFT equations are solved at this geometry
        return nemo_ptr->gradient(x);
    }

    double coord_chksum() const { return coords_sum; }

    const SCF& get_calc() const { return *nemo_ptr->get_calc(); }

    SCF& get_calc() { return *nemo_ptr->get_calc(); }

    /// return full orbital i, multiplied with the nuclear correlation factor

    /// note that nemo() and orbital() are the same if no nuclear
    /// correlation factor is used
    real_function_3d orbital(const int i) const {
        MADNESS_ASSERT(nemo_ptr->get_calc_param().spin_restricted());
        return orbitals_[i];
    }

    /// return full orbitals, multiplied with the nuclear correlation factor

    /// note that nemo() and orbital() are the same if no nuclear
    /// correlation factor is used
    std::vector<real_function_3d> orbitals() const {
        MADNESS_ASSERT(nemo_ptr->get_calc_param().spin_restricted());
        return orbitals_;
    }

    /// return orbitals, multiplied with the square nuclear correlation factor

    /// note that nemo() and orbital() are the same if no nuclear
    /// correlation factor is used
    std::vector<real_function_3d> R2orbitals() const {
        MADNESS_ASSERT(nemo_ptr->get_calc_param().spin_restricted());
        return R2orbitals_;
    }

    /// return orbital i, multiplied with the square nuclear correlation factor

    /// note that nemo() and orbital() are the same if no nuclear
    /// correlation factor is used
    real_function_3d R2orbital(const int i) const {
        MADNESS_ASSERT(nemo_ptr->get_calc_param().spin_restricted());
        return R2orbitals_[i];
    }

    /// return nemo i, which is the regularized orbital

    /// note that nemo() and orbital() are the same if no nuclear
    /// correlation factor is used
    real_function_3d nemo(const int i) const {
        MADNESS_ASSERT(nemo_ptr->get_calc_param().spin_restricted());
        return nemo_ptr->get_calc()->amo[i];
    }

    /// return nemo, which are the regularized orbitals

    /// note that nemo() and orbital() are the same if no nuclear
    /// correlation factor is used
    std::vector<real_function_3d> nemos() const {
        MADNESS_ASSERT(nemo_ptr->get_calc_param().spin_restricted());
        return nemo_ptr->get_calc()->amo;
    }

    /// return orbital energy i
    double orbital_energy(const int i) const {
        MADNESS_ASSERT(nemo_ptr->get_calc_param().spin_restricted());
        return nemo_ptr->get_calc()->aeps[i];
    }

    /// return the Coulomb potential
    real_function_3d get_coulomb_potential() const {
        MADNESS_ASSERT(nemo_ptr->get_calc_param().spin_restricted());
        if (coulomb.is_initialized()) return copy(coulomb);
        functionT rho = 2.0*nemo_ptr->compute_density(nemos())*nemo_ptr->R_square;
        coulomb = (nemo_ptr->get_calc()->make_coulomb_potential(rho)).truncate();
        return copy(coulomb);
    }

    /// return the nuclear potential
    real_function_3d get_nuclear_potential() const {
        return nemo_ptr->get_calc()->potentialmanager->vnuclear();
    }

    /// return the number of occupied orbitals
    int nocc() const {
        MADNESS_ASSERT(nemo_ptr->get_calc_param().spin_restricted());
        return nemo_ptr->get_calc_param().nalpha();
    }
};


/// enhanced POD for the pair functions
class ElectronPair : public archive::ParallelSerializableObject {

public:
    /// default ctor; initialize energies with a large number
    ElectronPair() : ElectronPair(-1, -1) {}

    /// ctor; initialize energies with a large number
    ElectronPair(const int i, const int j)
            : i(i), j(j), e_singlet(uninitialized()), e_triplet(uninitialized()),
              ij_gQf_ij(uninitialized()), ji_gQf_ij(uninitialized()), iteration(0), converged(false) {
    }

    /// print the pair's energy
    void print_energy() const {
        if (function.world().rank() == 0) {
            printf("final correlation energy %2d %2d %12.8f %12.8f\n",
                   i, j, e_singlet, e_triplet);
        }
    }

    static double uninitialized() { return 1.e10; }

    int i, j;                       ///< orbitals i and j
    real_function_6d function;      ///< pair function for a specific pair w/o correlation factor part
    real_function_6d constant_term;    ///< the first order contribution to the MP1 wave function

    double e_singlet;                ///< the energy of the singlet pair ij
    double e_triplet;                ///< the energy of the triplet pair ij

    double ij_gQf_ij;                ///< <ij | g12 Q12 f12 | ij>
    double ji_gQf_ij;                ///< <ji | g12 Q12 f12 | ij>

    int iteration;                    ///< current iteration for restart
    bool converged;                    ///< is the pair function converged

    /// serialize this ElectronPair

    /// store the function only if it has been initialized
    /// load the function only if there is one
    /// don't serialize recomputable intermediates r12phi, Uphi, KffKphi
    template<typename Archive>
    void serialize(Archive& ar) {
        bool fexist = function.is_initialized();
        bool cexist = constant_term.is_initialized();
        ar & ij_gQf_ij & ji_gQf_ij & e_singlet & e_triplet & converged
        & iteration & fexist & cexist;
        if (fexist) ar & function;
        if (cexist) ar & constant_term;
    }

    bool load_pair(World& world) {
        std::string name = "pair_" + stringify(i) + stringify(j);
        bool exists = archive::ParallelInputArchive<archive::BinaryFstreamInputArchive>::exists(world, name.c_str());
        if (exists) {
            if (world.rank() == 0) printf("loading matrix elements %s", name.c_str());
            archive::ParallelInputArchive<archive::BinaryFstreamInputArchive> ar(world, name.c_str(), 1);
            ar & *this;
            if (world.rank() == 0) printf(" %s\n", (converged) ? " converged" : " not converged");
            if (function.is_initialized()) function.set_thresh(FunctionDefaults<6>::get_thresh());
            if (constant_term.is_initialized()) constant_term.set_thresh(FunctionDefaults<6>::get_thresh());
        } else {
            if (world.rank() == 0) print("could not find pair ", i, j, " on disk");
        }
        return exists;
    }

    void store_pair(World& world) {
        std::string name = "pair_" + stringify(i) + stringify(j);
        if (world.rank() == 0) printf("storing matrix elements %s\n", name.c_str());
        archive::ParallelOutputArchive<archive::BinaryFstreamOutputArchive> ar(world, name.c_str(), 1);
        ar & *this;
    }

    // print information
    void info(World& world) const {
        if (world.rank() == 0) {
            std::cout << std::setw(20) << std::setfill(' ') << " *Information about Electron Pair " << i << j << " "
                      << std::setw(20) << std::setfill(' ') << std::endl;
            std::cout << std::setw(20) << std::setfill(' ') << " *e_singlet " << e_singlet << " " << std::setw(20)
                      << std::setfill(' ') << std::endl;
            std::cout << std::setw(20) << std::setfill(' ') << " *e_triplet " << e_triplet << " " << std::setw(20)
                      << std::setfill(' ') << std::endl;
            std::cout << std::setw(20) << std::setfill(' ') << " *ij_gQf_ij " << ij_gQf_ij << " " << std::setw(20)
                      << std::setfill(' ') << std::endl;
            std::cout << std::setw(20) << std::setfill(' ') << " *ji_gQf_ij " << ji_gQf_ij << " " << std::setw(20)
                      << std::setfill(' ') << std::endl;
        }
    }
};


/// a class for computing the first order wave function and MP2 pair energies
class MP2 : public OptimizationTargetInterface, public QCPropertyInterface {

    /// POD for MP2 keywords
    struct Parameters : public QCCalculationParametersBase {

        /// use OEP orbitals
        bool do_oep1 = false;

        Parameters() {
            /// the map with initial values
            initialize < double > ("thresh", 1.e-3, "recommended values: 1.e-4 < econv < 1.e-8");
            initialize < double > ("econv", 1.e-3, "recommended values: 1.e-4 < econv < 1.e-8");
            initialize < double > ("dconv", 1.e-3, "recommended values: 1.e-4 < econv < 1.e-8");
            initialize < std::vector<int> > ("pair", {-1, -1});
            initialize < int > ("freeze", 0);
            initialize < int > ("maxsub", 2);
            initialize < bool > ("restart", true);
            initialize < bool > ("no_compute", false);
            initialize < int > ("maxiter", 5);
        }

        /// ctor reading out the input file
        Parameters(World& world, const commandlineparser& parser) : Parameters() {
            read_and_set_derived_values(world,parser);

            // print final parameters
            if (world.rank() == 0) print("mp2", "end");
        }

        void read_and_set_derived_values(World& world, const commandlineparser& parser) {
            read_input_and_commandline_options(world, parser, "mp2");

            set_derived_value("dconv", sqrt(get<double>("econv")) * 0.1);
        	set_derived_value("thresh",get<double>("econv"));
        }

        /// check the user input
        void check_input(const std::shared_ptr<HartreeFock> hf) const {
            if (freeze() > hf->nocc()) MADNESS_EXCEPTION("you froze more orbitals than you have", 1);
            if (i() >= hf->nocc()) MADNESS_EXCEPTION("there is no i-th orbital", 1);
            if (j() >= hf->nocc()) MADNESS_EXCEPTION("there is no j-th orbital", 1);
            if (thresh() < 0.0) MADNESS_EXCEPTION("please provide the accuracy threshold for MP2", 1);
        }

        double thresh() const { return get<double>("thresh"); }    /// convenience function
        double econv() const { return get<double>("econv"); }        /// convenience function
        double dconv() const { return this->get<double>("dconv"); }        /// convenience function
        int freeze() const { return this->get<int>("freeze"); }            /// convenience function
        int i() const { return this->get<std::vector<int> >("pair")[0]; }    /// convenience function
        int j() const { return this->get<std::vector<int> >("pair")[1]; }    /// convenience function
        int restart() const { return this->get<bool>("restart"); }    /// convenience function
        int no_compute() const { return this->get<bool>("no_compute"); }    /// convenience function
        int maxiter() const { return this->get<int>("maxiter"); }    /// convenience function
        int maxsub() const { return this->get<int>("maxsub"); }    /// convenience function
        bool do_oep() const { return do_oep1;}
    };

    /// POD holding all electron pairs with easy access
    template<typename T>
    struct Pairs {

        typedef std::map<std::pair<int, int>, T> pairmapT;
        pairmapT allpairs;

        /// getter
        const T& operator()(int i, int j) const {
            return allpairs.find(std::make_pair(i, j))->second;
        }

        /// getter
        T& operator()(int i, int j) {
            return allpairs[std::make_pair(i, j)];
        }

        /// setter
        void insert(int i, int j, T pair) {
            std::pair<int, int> key = std::make_pair(i, j);
            allpairs.insert(std::make_pair(key, pair));
        }
    };


    World& world;                           ///< the world
    Parameters param;                        ///< SCF parameters for MP2
    std::shared_ptr<HartreeFock> hf;        ///< our reference
    CorrelationFactor corrfac;              ///< correlation factor: Slater
    std::shared_ptr<NuclearCorrelationFactor> nuclear_corrfac;
    mutable Tensor<double> fock;            ///< the Fock matrix

    Pairs<ElectronPair> pairs;       ///< pair functions and energies
    double correlation_energy;                ///< the correlation energy
    double coords_sum;                        ///< check sum for the geometry

    StrongOrthogonalityProjector<double, 3> Q12;

private:

    std::shared_ptr<real_convolution_3d> poisson;

public:

    /// ctor
    MP2(World& world, const commandlineparser& parser);
    std::string name() const {return "MP2";};

    static void help() {
        print_header2("help page for MP2 ");
        print("The mp2 code computes second order correlation energies based on a moldft or nemo calculation");
        print("You can print all available calculation parameters by running\n");
        print("mp2 --print_parameters\n");
        print("You can perform a simple calculation by running\n");
        print("mp2 --geometry=h2o.xyz\n");
        print("provided you have an xyz file in your directory.");

    }

    static void print_parameters() {
        Parameters param;
        print("default parameters for the mp2 program are");
        param.print("mp2", "end");
        print("\n\nthe molecular geometry must be specified in a separate block:");
        Molecule::print_parameters();
    }
    virtual bool selftest() {
        return true;
    };

    /// return a checksum for the geometry
    double coord_chksum() const { return coords_sum; }

    /// return the molecular correlation energy energy (without the HF energy)
    double value();

    /// return the molecular correlation energy as a function of the coordinates
    double value(const Tensor<double>& x);

    /// make sure frozen orbitals don't couple with correlated ones -- relocalize if necessary
    bool check_core_valence_separation() const;

    /// make sure frozen orbitals don't couple with correlated ones -- relocalize if necessary
    void enforce_core_valence_separation();

    /// return the underlying HF reference
    HartreeFock& get_hf() { return *hf; }

    /// return the 0th order energy of pair ij (= sum of orbital energies)
    double zeroth_order_energy(const int i, const int j) const {
        return hf->orbital_energy(i) + hf->orbital_energy(j);
    }

    /// solve the residual equation for electron pair (i,j)

    /// \todo Parameter documentation. Below are un-doxygenated comments that no longer seem relevant?
    // @param[in]  pair    electron pair to solve
    // @param[in]  econv   energy convergence criterion (for a single pair)
    // @param[in]  dconv   density convergence criterion (for a single pair)
    real_function_6d iterate(const real_function_6d& f) const {
        ElectronPair tmp(0, 0);
        tmp.function = copy(f);
        tmp.constant_term = copy(f);
        tmp.ij_gQf_ij = compute_gQf(0, 0, tmp);
        tmp.ji_gQf_ij = tmp.ij_gQf_ij;
        solve_residual_equations(tmp, 10.0, 10.0);
        return tmp.function;
    }

    void solve_residual_equations(ElectronPair& pair,
                                  const double econv, const double dconv) const;

    /// solve the coupled MP1 equations (e.g. for local orbitals)

    /// @param[in]  pairs   set of (coupled) electron pairs to solve
    /// @param[in]  econv   energy convergence criterion (for all pairs)
    /// @param[in]  dconv   density convergence criterion (for all pairs)
    double solve_coupled_equations(Pairs<ElectronPair>& pairs,
                                   const double econv, const double dconv) const;

    /// compute increments: psi^1 = C + GV C + GVGV C + GVGVGV C + ..

    /// for pre-optimization of the mp1 wave function
    /// no restart option here for now
    /// param[inout]	pair	the electron pair
    /// param[in]		green	the Green's function
    void increment(ElectronPair& pair, real_convolution_6d& green);

    double asymmetry(const real_function_6d& f, const std::string s) const;

    /// compute the matrix element <ij | g12 Q12 f12 | phi^0>

    /// scales quartically. I think I can get this down to cubically by
    /// setting Q12 = (1 - O1)(1 - O2) = 1- O1(1 - 0.5 O2) - O2 (1 - 0.5 O1)
    /// as for the formulas cf the article mra_molecule
    /// @return 	the energy <ij | g Q f | kl>
    double compute_gQf_cc2interface(const int i, const int j, const real_function_6d& f) const {
        ElectronPair tmp(0, 0);
        tmp.function = copy(f);
        tmp.constant_term = copy(f);
        return compute_gQf(i, j, tmp);
    }

    double compute_gQf(const int i, const int j, ElectronPair& pair) const;

    /// pretty print the options

    /// invoke only in (world.rank()==0) !!
    template<typename T>
    void print_options(const std::string option, const T val) const {
        std::cout << std::setfill(' ') << std::setw(30);
        std::cout << option << "  " << val << std::endl;
    }

    real_function_6d debug_cc2(const real_function_6d& f, const size_t& i, const size_t& j) const {
        return multiply_with_0th_order_Hamiltonian(f, i, j);
    }

public:

    /// save a function
    template<typename T, size_t NDIM>
    void save_function(const Function<T, NDIM>& f, const std::string name) const;

    /// load a function
    template<typename T, size_t NDIM>
    void load_function(Function<T, NDIM>& f, const std::string name) const;

private:
    /// return the function Uphi0; load from disk if available
    real_function_6d make_Uphi0(ElectronPair& pair) const;

    /// return the function [K,f] phi0; load from disk if available
    real_function_6d make_KffKphi0(const ElectronPair& pair) const;


    /// compute some matrix elements that don't change during the SCF
    ElectronPair make_pair(const int i, const int j) const;

    /// compute the first iteration of the residual equations and all intermediates
    void guess_mp1_3(ElectronPair& pair) const;

private:

    /// compute the singlet and triplet energy for a given electron pair

    /// @return	the energy of 1 degenerate triplet and 1 singlet pair
    double compute_energy(ElectronPair& pair) const;

    /// apply the exchange operator on an orbital

    /// @param[in]	phi the orbital
    /// @param[in]	hc	hermitian conjugate -> swap bra and ket space
    /// @return 	Kphi
    real_function_3d K(const real_function_3d& phi, const bool hc = false) const;

    /// apply the exchange operator on a pair function

    /// @param[in]	phi the pair function
    /// @param[in]	is_symmetric is the function symmetric wrt particle exchange
    /// @return 	(K1 + K2) |phi >
    real_function_6d K(const real_function_6d& phi, const bool is_symmetric = false) const;

    /// apply the Coulomb operator a on orbital

    /// @param[in]	phi the orbital
    /// @return 	Jphi
    real_function_3d J(const real_function_3d& phi) const;

    real_function_6d apply_exchange_vector(const real_function_6d& f, const int particle) const;

    /// apply the exchange operator on f

    /// if the exchange operator is similarity transformed (R-1 K R) the
    /// orbital spaces differ for the orbitals underneath the integral sign
    /// R-1 K R = \phi-ket(1) * \int \phi-bra(1') * f(1',2)
    /// @param[in]  f   the pair function
    /// @param[in]  orbital_bra the orbital underneath the integral sign (typically R2orbitals)
    /// @param[in]  orbital_ket the orbital to be pre-multiplied with (typically orbitals)
    /// @return     the pair function, on which the exchange operator has been applied
    real_function_6d apply_exchange(const real_function_6d& f,
                                    const real_function_3d& orbital_ket,
                                    const real_function_3d& orbital_bra, const int particle) const;

    /// make the quantity chi_k

    /// chi is the Poisson kernel applied on an orbital product of the
    /// input function and the vector of orbitals
    /// \f[ \chi_{k{*} i}(1) = \int dr_2 \frac{k(2) i(2)}{|r_1-r_2|} \f]
    /// \f[ \chi_{ki{*}}(1) = \int dr_2 \frac{k(2) i(2)}{|r_1-r_2|} \f] if hc
    /// @param[in]  phi		orbital phi_i
    /// @param[in]	op		the operator in SeparatedConvolution form
    /// @param[in]	hc		compute hermitian conjugate -> pass the correct phi!
    /// @return a vector of length nocc
    std::vector<real_function_3d> make_chi(const real_function_3d& phi,
                                           const real_convolution_3d& op,
                                           const bool hc = false) const;

    /// make the quantity xi_k

    /// xi is chi multiplied with an orbital j
    /// \f[ \xi_{k{*}i,j}(1) = \chi_{ki}(1) j(1) \f]
    /// \f[ \xi_{ki{*},j{*}}(1) = \chi_{k{*}i}(1) j(1) \f]  if hc
    /// @param[in]  phi_i   orbital i
    /// @param[in]  phi_j   orbital j
    /// @param[in]	op		the operator in SeparatedConvolution form
    /// @param[in]	hc		compute hermitian conjugate  -> pass the correct phi!
    /// @return a vector of length k=0,..,nocc
    std::vector<real_function_3d> make_xi(const real_function_3d& phi_i,
                                          const real_function_3d& phi_j, const real_convolution_3d& op,
                                          const bool hc = false) const;

    /// apply the operator K on the reference and multiply with f; fK |phi^0>

    /// @param[in]  i   index of orbital i
    /// @param[in]  j   index of orbital j
    /// @return     the function f12 (K(1)+K(2))|phi^0>
    real_function_6d make_fKphi0(const int i, const int j) const;

    /// return the function (J(1)-K(1)) |phi0> as on-demand function

    /// @param[in]	hc		compute hermitian conjugate -> swap bra and ket space
    real_function_6d JK1phi0_on_demand(const int i, const int j,
                                       const bool hc = false) const;

    /// return the function (J(2)-K(2)) |phi0> as on-demand function
    real_function_6d JK2phi0_on_demand(const int i, const int j,
                                       const bool hc = false) const;

    /// return the function |phi0> as on-demand function
    real_function_6d phi0_on_demand(const int i, const int j) const;

    /// return the function |F1F2> as on-demand function
    real_function_6d nemo0_on_demand(const int i, const int j) const;


    /// multiply the given function with the 0th order Hamiltonian, exluding the 0th order energy

    /// @param[in]  f   the function we apply H^0 on
    /// @return     the function g=H^0 f, which is NOT orthogonalized against f
    real_function_6d multiply_with_0th_order_Hamiltonian(const real_function_6d& f,
                                                         const int i, const int j) const;

    // need this for CC2 debug
public:
    real_function_6d get_residue(const real_function_6d& f,
                                 const int i, const int j) {
        hf->value();        // make sure the reference is converged
        nuclear_corrfac = hf->nemo_ptr->ncf;
        // set all orbitals spaces
        // When a nuclear correlation factor is used the residual equations
        // are similarity transformed. Therefore the orbitals in the
        // projection operator must be set accordingly.
        if (nuclear_corrfac->type() == NuclearCorrelationFactor::None) {
            Q12.set_spaces(hf->get_calc().amo);
        } else {
            // only valid for closed shell
            MADNESS_ASSERT(hf->get_calc().param.spin_restricted());
            const std::vector<real_function_3d>& nemos = hf->nemos();
            const std::vector<real_function_3d>& R2amo = hf->R2orbitals();
            Q12.set_spaces(R2amo, nemos, R2amo, nemos);
            if (world.rank() == 0) {
                print("set orbital spaces for the SO projector");
                print("Q12,R = (1-|nemo><nemo|R2) (1-|nemo><nemo|R2)");
            }
        }

        return multiply_with_0th_order_Hamiltonian(f, i, j);
    }

private:
    // compute some intermediates
    Tensor<double> get_fock_matrix() const {
        if (fock.has_data()) return copy(fock);
        const tensorT occ = hf->get_calc().aocc;
        fock = hf->nemo_ptr->compute_fock_matrix(hf->nemos(), occ);
        if (world.rank() == 0 and hf->nocc() < 10) {
            print("The Fock matrix");
            print(fock);
        }
        return copy(fock);
    }



    /// add the coupling terms for local MP2

    /// \sum_{k\neq i} f_ki |u_kj> + \sum_{l\neq j} f_lj |u_il>
    /// @todo Verify this doxygen block.
    void add_local_coupling(const Pairs<ElectronPair>& pairs, Pairs<real_function_6d>& coupling) const;

    mutable double ttt, sss;

    void START_TIMER(World& world) const {
        world.gop.fence();
        ttt = wall_time();
        sss = cpu_time();
    }

    void END_TIMER(World& world, const char *msg) const {
        ttt = wall_time() - ttt;
        sss = cpu_time() - sss;
        if (world.rank() == 0) printf("timer: %20.20s %8.2fs %8.2fs\n", msg, sss, ttt);
    }

};
};

#endif /* MP2_H_ */

