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

//#define WORLD_INSTANTIATE_STATIC_TEMPLATES
#include <madness/mra/mra.h>
#include <madness/mra/lbdeux.h>
#include <chem/QCCalculationParametersBase.h>
#include <chem/SCF.h>
#include <madness/mra/nonlinsol.h>
#include <chem/projector.h>
#include <chem/correlationfactor.h>
#include <chem/electronic_correlation_factor.h>
#include <chem/nemo.h>

#include <madness/world/text_fstream_archive.h>
using madness::archive::TextFstreamInputArchive;
using madness::archive::TextFstreamOutputArchive;


#include <iostream>

using namespace madness;

namespace madness {

    struct LBCost {
        double leaf_value;
        double parent_value;
        LBCost(double leaf_value=1.0, double parent_value=1.0)
            : leaf_value(leaf_value)
            , parent_value(parent_value)
        {}

        double operator()(const Key<6>& key, const FunctionNode<double,6>& node) const {
            return node.coeff().size();
        }
    };


    class HartreeFock {
        World& world;
        std::shared_ptr<SCF> calc;
        mutable double coords_sum;     // sum of square of coords at last solved geometry

        // save the Coulomb potential
        mutable real_function_3d coulomb;

        /// reconstructed orbitals: R * phi, where R is the nuclear correlation factor
        std::vector<real_function_3d> orbitals_;

        /// R^2 * phi, where R is the nuclear correlation factor, corresponds
        /// to the bra space of the transformed operators
        std::vector<real_function_3d> R2orbitals_;

    public:

        Nemo nemo_calc;

        HartreeFock(World& world, std::shared_ptr<SCF> calc1) :
        	world(world), calc(calc1), coords_sum(-1.0), nemo_calc(world,calc1,"input") {
        }

        bool provides_gradient() const {return true;}

        double value() {
            return value(calc->molecule.get_all_coords());
        }

        double value(const Tensor<double>& x) {

        	// fast return if the reference is already solved at this geometry
            double xsq = x.sumsq();
            if (xsq == coords_sum) return calc->current_energy;

            calc->molecule.set_all_coords(x.reshape(calc->molecule.natom(),3));
            coords_sum = xsq;

            // some extra steps if we have a nuclear correlation factor
//            if (1) {
//        	if (nemo_calc.nuclear_correlation->type()==NuclearCorrelationFactor::GaussSlater) {
        		// converge the nemo equations
        		nemo_calc.value(x);

//        	} else {
//				// Make the nuclear potential, initial orbitals, etc.
//				calc->make_nuclear_potential(world);
//				calc->potentialmanager->vnuclear().print_size("vnuc");
//				calc->ao=calc->project_ao_basis(world,calc->aobasis);
//
//				// read converged wave function from disk if there is one
//				if (calc->param.no_compute()) {
//					calc->load_mos(world);
//					return calc->current_energy;
//				}
//
//				if (calc->param.restart()) {
//					calc->load_mos(world);
//				} else {
//					calc->initial_guess(world);
//					calc->param.restart = true;
//				}
//
//				// If the basis for the inital guess was not sto-3g
//				// switch to sto-3g since this is needed for analysis
//				// of the MOs and orbital localization
//				if (calc->param.aobasis() != "sto-3g") {
//					calc->param.aobasis = "sto-3g";
//					calc->project_ao_basis(world);
//				}
//
//
//				calc->solve(world);
//				calc->save_mos(world);
//
//				// successively tighten threshold
//				if (calc->param.econv<1.1e-6) {
//					calc->set_protocol<3>(world,1e-6);
//					calc->make_nuclear_potential(world);
//					calc->project_ao_basis(world);
//					calc->project(world);
//					calc->solve(world);
//					calc->save_mos(world);
//				}
//
//				calc->save_mos(world);
//			}

            // compute the full, reconstructed orbitals from nemo
            orbitals_=mul(world,nemo_calc.R,nemo_calc.get_calc()->amo);
            real_function_3d R2=nemo_calc.ncf->square();
            R2orbitals_=mul(world,R2,nemo_calc.get_calc()->amo);

            return calc->current_energy;
        }

        Tensor<double> gradient(const Tensor<double>& x) {

            value(x); // Ensures DFT equations are solved at this geometry
            return nemo_calc.gradient(x);
        }

        double coord_chksum() const {return coords_sum;}

        const SCF& get_calc() const {return *calc;}
        SCF& get_calc() {return *calc;}

        /// return full orbital i, multiplied with the nuclear correlation factor

        /// note that nemo() and orbital() are the same if no nuclear
        /// correlation factor is used
        real_function_3d orbital(const int i) const {
            MADNESS_ASSERT(calc->param.spin_restricted());
            return orbitals_[i];
        }

        /// return full orbitals, multiplied with the nuclear correlation factor

        /// note that nemo() and orbital() are the same if no nuclear
        /// correlation factor is used
        std::vector<real_function_3d> orbitals() const {
            MADNESS_ASSERT(calc->param.spin_restricted());
            return orbitals_;
        }

        /// return orbitals, multiplied with the square nuclear correlation factor

        /// note that nemo() and orbital() are the same if no nuclear
        /// correlation factor is used
        std::vector<real_function_3d> R2orbitals() const {
            MADNESS_ASSERT(calc->param.spin_restricted());
            return R2orbitals_;
        }

        /// return orbital i, multiplied with the square nuclear correlation factor

        /// note that nemo() and orbital() are the same if no nuclear
        /// correlation factor is used
        real_function_3d R2orbital(const int i) const {
            MADNESS_ASSERT(calc->param.spin_restricted());
            return R2orbitals_[i];
        }

        /// return nemo i, which is the regularized orbital

        /// note that nemo() and orbital() are the same if no nuclear
        /// correlation factor is used
        real_function_3d nemo(const int i) const {
            MADNESS_ASSERT(calc->param.spin_restricted());
            return calc->amo[i];
        }

        /// return nemo, which are the regularized orbitals

        /// note that nemo() and orbital() are the same if no nuclear
        /// correlation factor is used
        std::vector<real_function_3d> nemos() const {
            MADNESS_ASSERT(calc->param.spin_restricted());
            return calc->amo;
        }

        /// return orbital energy i
        double orbital_energy(const int i) const {
            MADNESS_ASSERT(calc->param.spin_restricted());
            return calc->aeps[i];
        }

        /// return the Coulomb potential
        real_function_3d get_coulomb_potential() const {
            MADNESS_ASSERT(calc->param.spin_restricted());
            if (coulomb.is_initialized()) return copy(coulomb);
            functionT rho = calc->make_density(world, calc->aocc, orbitals()).scale(2.0);
            coulomb=calc->make_coulomb_potential(rho);
            return copy(coulomb);
        }

        /// return the nuclear potential
        real_function_3d get_nuclear_potential() const {
            return calc->potentialmanager->vnuclear();
        }

        /// return the number of occupied orbitals
        int nocc() const {
            MADNESS_ASSERT(calc->param.spin_restricted());
            return calc->param.nalpha();
        }
    };


    /// enhanced POD for the pair functions
    class ElectronPair : public archive::ParallelSerializableObject {

    public:
    	/// default ctor; initialize energies with a large number
    	ElectronPair() : ElectronPair(-1,-1) {}

    	/// ctor; initialize energies with a large number
    	ElectronPair(const int i, const int j)
    		: i(i), j(j), e_singlet(uninitialized()), e_triplet(uninitialized()),
    		  ij_gQf_ij(uninitialized()), ji_gQf_ij(uninitialized()), iteration(0), converged(false) {
    	}

    	/// print the pair's energy
    	void print_energy() const {
            if (function.world().rank()==0) {
            	printf("final correlation energy %2d %2d %12.8f %12.8f\n",
            			i,j,e_singlet,e_triplet);
            }
    	}

    	static double uninitialized() {return 1.e10;}

        int i, j;                       ///< orbitals i and j
        real_function_6d function;      ///< pair function for a specific pair w/o correlation factor part
        real_function_6d constant_term;	///< the first order contribution to the MP1 wave function

        double e_singlet;				///< the energy of the singlet pair ij
        double e_triplet;				///< the energy of the triplet pair ij

        double ij_gQf_ij;         	  	///< <ij | g12 Q12 f12 | ij>
        double ji_gQf_ij;          		///< <ji | g12 Q12 f12 | ij>

        int iteration;					///< current iteration for restart
        bool converged;					///< is the pair function converged

        /// serialize this ElectronPair

        /// store the function only if it has been initialized
        /// load the function only if there is one
        /// don't serialize recomputable intermediates r12phi, Uphi, KffKphi
        template <typename Archive> void serialize (Archive& ar) {
        	bool fexist=function.is_initialized();
        	bool cexist=constant_term.is_initialized();
			ar & ij_gQf_ij & ji_gQf_ij & e_singlet & e_triplet & converged
				& iteration & fexist & cexist;
			if (fexist) ar & function;
			if (cexist) ar & constant_term;
        }

        bool load_pair(World& world) {
        	std::string name="pair_"+stringify(i)+stringify(j);
        	bool exists=archive::ParallelInputArchive::exists(world,name.c_str());
            if (exists) {
            	if (world.rank()==0) printf("loading matrix elements %s",name.c_str());
                archive::ParallelInputArchive ar(world, name.c_str(), 1);
                ar & *this;
                if (world.rank()==0) printf(" %s\n",(converged)?" converged":" not converged");
                if (function.is_initialized()) function.set_thresh(FunctionDefaults<6>::get_thresh());
                if (constant_term.is_initialized()) constant_term.set_thresh(FunctionDefaults<6>::get_thresh());
            } else {
                if (world.rank()==0) print("could not find pair ",i,j," on disk");
            }
            return exists;
        }

        void store_pair(World& world) {
        	std::string name="pair_"+stringify(i)+stringify(j);
        	if (world.rank()==0) printf("storing matrix elements %s\n",name.c_str());
            archive::ParallelOutputArchive ar(world, name.c_str(), 1);
        	ar & *this;
        }

    	// print information
    	void info(World &world)const{
    		if(world.rank()==0){
    			std::cout <<std::setw(20) << std::setfill(' ') << " *Information about Electron Pair " << i << j << " " <<std::setw(20) << std::setfill(' ') << std::endl;
    			std::cout <<std::setw(20) << std::setfill(' ') << " *e_singlet " << e_singlet << " " <<std::setw(20) << std::setfill(' ') << std::endl;
    			std::cout <<std::setw(20) << std::setfill(' ') << " *e_triplet " << e_triplet << " " <<std::setw(20) << std::setfill(' ') << std::endl;
    			std::cout <<std::setw(20) << std::setfill(' ') << " *ij_gQf_ij " << ij_gQf_ij << " " <<std::setw(20) << std::setfill(' ') << std::endl;
    			std::cout <<std::setw(20) << std::setfill(' ') << " *ji_gQf_ij " << ji_gQf_ij << " " <<std::setw(20) << std::setfill(' ') << std::endl;
    		}
    	}
    };


    /// a class for computing the first order wave function and MP2 pair energies
    class MP2 : public OptimizationTargetInterface {

    	/// POD for MP2 keywords
    	struct Parameters : public QCCalculationParametersBase {

        	/// use OEP orbitals
        	bool do_oep;

        	/// ctor reading out the input file
        	Parameters(World& world) {

                /// the map with initial values
        		initialize<double>("thresh",1.e-3,"recommended values: 1.e-4 < econv < 1.e-8");
        		initialize<double>("econv",1.e-3,"recommended values: 1.e-4 < econv < 1.e-8");
        		initialize<double>("dconv",1.e-3,"recommended values: 1.e-4 < econv < 1.e-8");
        		initialize<std::vector<int> >("pair",{-1,-1});
        		initialize<int>("freeze",0);
        		initialize<int>("maxsub",2);
        		initialize<bool>("restart",true);
        		initialize<int>("maxiter",5);

        		read_and_set_derived_values(world);

        		// print final parameters
        		if (world.rank()==0) print("mp2","end");
        	}

        	void read_and_set_derived_values(World& world) {
        		read(world,"input","mp2");
        		set_derived_value("dconv",sqrt(get<double>("econv"))*0.1);
        	}

            /// check the user input
        	void check_input(const std::shared_ptr<HartreeFock> hf) const {
                if (freeze()>hf->nocc()) MADNESS_EXCEPTION("you froze more orbitals than you have",1);
                if (i()>=hf->nocc()) MADNESS_EXCEPTION("there is no i-th orbital",1);
                if (j()>=hf->nocc()) MADNESS_EXCEPTION("there is no j-th orbital",1);
                if (thresh()<0.0) MADNESS_EXCEPTION("please provide the accuracy threshold for MP2",1);
        	}

        	double thresh() const {return get<double>("thresh");}   	/// convenience function
        	double econv() const {return get<double>("econv");}   		/// convenience function
        	double dconv() const {return this->get<double>("dconv");}   		/// convenience function
        	int freeze() const {return this->get<int>("freeze");}   			/// convenience function
        	int i() const {return this->get<std::vector<int> >("pair")[0];}	/// convenience function
        	int j() const {return this->get<std::vector<int> >("pair")[1];}	/// convenience function
        	int restart() const {return this->get<bool>("restart");}	/// convenience function
        	int maxiter() const {return this->get<int>("maxiter");}	/// convenience function
        	int maxsub() const {return this->get<int>("maxsub");}	/// convenience function
        };

        /// POD holding all electron pairs with easy access
        template<typename T>
        struct Pairs {

            typedef std::map<std::pair<int,int>, T> pairmapT;
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
        Parameters param;						///< SCF parameters for MP2
        std::shared_ptr<HartreeFock> hf;        ///< our reference
        CorrelationFactor corrfac;              ///< correlation factor: Slater
        std::shared_ptr<NuclearCorrelationFactor> nuclear_corrfac;
        mutable Tensor<double> fock;			///< the Fock matrix

        Pairs<ElectronPair> pairs;       ///< pair functions and energies
        double correlation_energy;				///< the correlation energy
        double coords_sum;						///< check sum for the geometry

        StrongOrthogonalityProjector<double,3> Q12;

    private:

        std::shared_ptr<real_convolution_3d> poisson;

    public:
        void apply_Q12(real_function_6d &f)const{
        	f = Q12(f);
        }
        /// ctor
        MP2(World& world, const std::string& input);

        /// return a checksum for the geometry
        double coord_chksum() const {return coords_sum;}

        /// return the molecular correlation energy energy (without the HF energy)
        double value();

        /// return the molecular correlation energy as a function of the coordinates
        double value(const Tensor<double>& x);

        /// print the SCF parameters
        void print_info(World& world) const;

        /// return the underlying HF reference
        HartreeFock& get_hf() {return *hf;}

        /// return the 0th order energy of pair ij (= sum of orbital energies)
        double zeroth_order_energy(const int i, const int j) const {
            return hf->orbital_energy(i)+hf->orbital_energy(j);
        }

        /// solve the residual equation for electron pair (i,j)

        /// \todo Parameter documentation. Below are un-doxygenated comments that no longer seem relevant?
        // @param[in]  pair    electron pair to solve
        // @param[in]  econv   energy convergence criterion (for a single pair)
        // @param[in]  dconv   density convergence criterion (for a single pair)
        real_function_6d iterate(const real_function_6d &f)const{
        	ElectronPair tmp(0,0);
        	tmp.function = copy(f);
        	tmp.constant_term = copy(f);
        	tmp.ij_gQf_ij = compute_gQf(0,0,tmp);
        	tmp.ji_gQf_ij = tmp.ij_gQf_ij;
        	solve_residual_equations(tmp,10.0,10.0);
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

        real_function_6d make_Rpsi(const ElectronPair& pair) const;

		/// compute increments: psi^1 = C + GV C + GVGV C + GVGVGV C + ..

        /// for pre-optimization of the mp1 wave function
    	/// no restart option here for now
        /// param[inout]	pair	the electron pair
        /// param[in]		green	the Green's function
        void increment(ElectronPair& pair, real_convolution_6d& green);

        /// swap particles 1 and 2

        /// param[in]	f	a function of 2 particles f(1,2)
        /// return	the input function with particles swapped g(1,2) = f(2,1)
        real_function_6d swap_particles(const real_function_6d& f) const;

        double asymmetry(const real_function_6d& f, const std::string s) const;

        void test(const std::string filename);

        /// compute the matrix element <ij | g12 Q12 f12 | phi^0>

        /// scales quartically. I think I can get this down to cubically by
        /// setting Q12 = (1 - O1)(1 - O2) = 1- O1(1 - 0.5 O2) - O2 (1 - 0.5 O1)
        /// as for the formulas cf the article mra_molecule
        /// @return 	the energy <ij | g Q f | kl>
        double compute_gQf_cc2interface(const int i, const int j, const real_function_6d &f)const{
        	ElectronPair tmp(0,0);
        	tmp.function = copy(f);
        	tmp.constant_term = copy(f);
        	return compute_gQf(i,j,tmp);
        }
        double compute_gQf(const int i, const int j, ElectronPair& pair) const;

        /// pretty print the options

        /// invoke only in (world.rank()==0) !!
        template<typename T>
        void print_options(const std::string option, const T val) const {
            std::cout << std::setfill (' ') << std::setw(30);
            std::cout << option << "  " << val << std::endl;
        }

        real_function_6d debug_cc2(const real_function_6d &f, const size_t &i, const size_t &j) const{
        	return  multiply_with_0th_order_Hamiltonian(f,i,j);
        }

    public:

        /// save a function
        template<typename T, size_t NDIM>
        void save_function(const Function<T,NDIM>& f, const std::string name) const;

        /// load a function
        template<typename T, size_t NDIM>
        void load_function(Function<T,NDIM>& f, const std::string name) const;

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
        real_function_3d K(const real_function_3d& phi, const bool hc=false) const;

        /// apply the exchange operator on a pair function

        /// @param[in]	phi the pair function
        /// @param[in]	is_symmetric is the function symmetric wrt particle exchange
        /// @return 	(K1 + K2) |phi >
        real_function_6d K(const real_function_6d& phi, const bool is_symmetric=false) const;

        /// apply the Coulomb operator a on orbital

        /// @param[in]	phi the orbital
        /// @return 	Jphi
        real_function_3d J(const real_function_3d& phi) const;

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
        		const bool hc=false) const;

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
        		const bool hc=false) const;

         /// apply the operator K on the reference and multiply with f; fK |phi^0>

        /// @param[in]  i   index of orbital i
        /// @param[in]  j   index of orbital j
        /// @return     the function f12 (K(1)+K(2))|phi^0>
        real_function_6d make_fKphi0(const int i, const int j) const;

        /// return the function (J(1)-K(1)) |phi0> as on-demand function

        /// @param[in]	hc		compute hermitian conjugate -> swap bra and ket space
        real_function_6d JK1phi0_on_demand(const int i, const int j,
        		const bool hc=false) const;

        /// return the function (J(2)-K(2)) |phi0> as on-demand function
        real_function_6d JK2phi0_on_demand(const int i, const int j,
        		const bool hc=false) const;

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
		const int i, const int j){
        	hf->value();		// make sure the reference is converged
        	nuclear_corrfac = hf->nemo_calc.ncf;
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

			return multiply_with_0th_order_Hamiltonian(f,i,j);
		}
    private:
		// compute some intermediates
        Tensor<double> get_fock_matrix() const {
        	if (fock.has_data()) return copy(fock);
    		const tensorT occ=hf->get_calc().aocc;
    		fock=hf->nemo_calc.compute_fock_matrix(hf->nemos(),occ);
    		if (world.rank()==0 and hf->nocc()<10) {
    			print("The Fock matrix");
    			print(fock);
    		}
    		return copy(fock);
        }



        /// add the coupling terms for local MP2

        /// \sum_{k\neq i} f_ki |u_kj> + \sum_{l\neq j} f_lj |u_il>
        /// @todo Verify this doxygen block.
        void add_local_coupling(const Pairs<ElectronPair>& pairs,
                Pairs<real_function_6d>& coupling) const;

        mutable double ttt, sss;
        void START_TIMER(World& world) const {
            world.gop.fence(); ttt=wall_time(); sss=cpu_time();
        }

        void END_TIMER(World& world, const char* msg) const {
            ttt=wall_time()-ttt; sss=cpu_time()-sss;
            if (world.rank()==0) printf("timer: %20.20s %8.2fs %8.2fs\n", msg, sss, ttt);
        }

    };
};

#endif /* MP2_H_ */

