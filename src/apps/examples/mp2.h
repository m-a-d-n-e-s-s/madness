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

  $Id$
*/
/*!
  \file mp2.h
  \brief Solves molecular MP2 equations
  \defgroup Solves molecular MP2 equations
  \ingroup examples

  The source is
  <a href=http://code.google.com/p/m-a-d-n-e-s-s/source/browse/local/trunk/src/apps/examples/mp2.h>here</a>.


*/


#ifndef MP2_H_
#define MP2_H_


#define WORLD_INSTANTIATE_STATIC_TEMPLATES
#include <mra/mra.h>
#include <mra/lbdeux.h>
#include <chem/SCF.h>
#include <examples/nonlinsol.h>
#include <chem/projector.h>
#include <chem/correlationfactor.h>
#include <examples/nemo.h>

#include <iostream>

static const double dcut=1.e-7;
static const double lo=1.e-6;
static const double bsh_eps=1.e-7;

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
	    return node.coeff().rank();
//            if (node.is_leaf()) {
//                return leaf_value;
//            } else {
//                return parent_value;
//            }
        }
    };

    /// do some load-balancing

    /// @param[in]	f		the function we want to distribute evenly
    /// @param[in]	leaf	if true: weigh leaf nodes only; if false: weigh internal nodes only
    void load_balance(const real_function_6d& f, const bool leaf) {
        LoadBalanceDeux<6> lb(f.world());
        if (leaf) lb.add_tree(f,LBCost(1.0,0.1));
        else lb.add_tree(f,LBCost(0.001,1.0));
        FunctionDefaults<6>::redistribute(f.world(), lb.load_balance(2.0,false));
//        if(f.world().rank() == 0) printf("redistributed at time   %.1fs\n", wall_time());
    }

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
        	world(world), calc(calc1), coords_sum(-1.0), nemo_calc(world,calc1) {
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
            if (1) {
//        	if (nemo_calc.nuclear_correlation->type()==NuclearCorrelationFactor::GaussSlater) {
        		// converge the nemo equations
        		nemo_calc.value(x);

        	} else {
				// Make the nuclear potential, initial orbitals, etc.
				calc->make_nuclear_potential(world);
				calc->potentialmanager->vnuclear().print_size("vnuc");
				calc->project_ao_basis(world);

				// read converged wave function from disk if there is one
				if (calc->param.no_compute) {
					calc->load_mos(world);
					return calc->current_energy;
				}

				if (calc->param.restart) {
					calc->load_mos(world);
				} else {
					calc->initial_guess(world);
					calc->param.restart = true;
				}

				// If the basis for the inital guess was not sto-3g
				// switch to sto-3g since this is needed for analysis
				// of the MOs and orbital localization
				if (calc->param.aobasis != "sto-3g") {
					calc->param.aobasis = "sto-3g";
					calc->project_ao_basis(world);
				}


				calc->solve(world);
				calc->save_mos(world);

				// successively tighten threshold
				if (calc->param.econv<1.1e-6) {
					calc->set_protocol<3>(world,1e-6);
					calc->make_nuclear_potential(world);
					calc->project_ao_basis(world);
					calc->project(world);
					calc->solve(world);
					calc->save_mos(world);
				}

				calc->save_mos(world);
			}

            // compute the full, reconstructed orbitals from nemo
            orbitals_=mul(world,nemo_calc.R,nemo_calc.get_calc()->amo);
            real_function_3d R2=nemo_calc.nuclear_correlation->square();
            R2orbitals_=mul(world,R2,nemo_calc.get_calc()->amo);

            return calc->current_energy;
        }

        Tensor<double> gradient(const Tensor<double>& x) {

            value(x); // Ensures DFT equations are solved at this geometry
            return calc->derivatives(world);
        }

        double coord_chksum() const {return coords_sum;}

        const SCF& get_calc() const {return *calc;}
        SCF& get_calc() {return *calc;}

        /// return full orbital i, multiplied with the nuclear correlation factor

        /// note that nemo() and orbital() are the same if no nuclear
        /// correlation factor is used
        real_function_3d orbital(const int i) const {
            MADNESS_ASSERT(calc->param.spin_restricted);
            return orbitals_[i];
        }

        /// return full orbitals, multiplied with the nuclear correlation factor

        /// note that nemo() and orbital() are the same if no nuclear
        /// correlation factor is used
        std::vector<real_function_3d> orbitals() const {
            MADNESS_ASSERT(calc->param.spin_restricted);
            return orbitals_;
        }

        /// return orbitals, multiplied with the square nuclear correlation factor

        /// note that nemo() and orbital() are the same if no nuclear
        /// correlation factor is used
        std::vector<real_function_3d> R2orbitals() const {
            MADNESS_ASSERT(calc->param.spin_restricted);
            return R2orbitals_;
        }

        /// return orbital i, multiplied with the square nuclear correlation factor

        /// note that nemo() and orbital() are the same if no nuclear
        /// correlation factor is used
        real_function_3d R2orbital(const int i) const {
            MADNESS_ASSERT(calc->param.spin_restricted);
            return R2orbitals_[i];
        }

        /// return nemo i, which is the regularized orbital

        /// note that nemo() and orbital() are the same if no nuclear
        /// correlation factor is used
        real_function_3d nemo(const int i) const {
            MADNESS_ASSERT(calc->param.spin_restricted);
            return calc->amo[i];
        }

        /// return nemo, which are the regularized orbitals

        /// note that nemo() and orbital() are the same if no nuclear
        /// correlation factor is used
        std::vector<real_function_3d> nemos() const {
            MADNESS_ASSERT(calc->param.spin_restricted);
            return calc->amo;
        }

        /// return orbital energy i
        double orbital_energy(const int i) const {
            MADNESS_ASSERT(calc->param.spin_restricted);
            return calc->aeps[i];
        }

        /// return the Coulomb potential
        real_function_3d get_coulomb_potential() const {
            MADNESS_ASSERT(calc->param.spin_restricted);
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
            MADNESS_ASSERT(calc->param.spin_restricted);
            return calc->param.nalpha;
        }
    };


    /// enhanced POD for the pair functions
    class ElectronPair : public archive::ParallelSerializableObject {

    public:
    	/// default ctor; initialize energies with a large number
    	ElectronPair()
    		: i(-1), j(-1), e_singlet(uninitialized()), e_triplet(uninitialized()),
    			ij_gQf_ij(uninitialized()), ji_gQf_ij(uninitialized()), iteration(0), converged(false) {
    	}

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
        real_function_6d r12phi;      	///< orbital product multiplied with the correlation factor
        real_function_6d constant_term;	///< the first order contribution to the MP1 wave function

        real_function_6d Uphi0;         ///< the function U |phi^0>  (U being Kutzelnigg's potential)
        real_function_6d KffKphi0;      ///< the function [K,f12] |phi^0>
        std::vector<real_function_3d> phi_k_UK_phi0;	///< < k(1) | U-K | phi^0(1,2)>
        std::vector<real_function_3d> phi_l_UK_phi0;	///< < l(2) | U-K | phi^0(1,2)>

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
    };


    /// a class for computing the first order wave function and MP2 pair energies
    class MP2 : public OptimizationTargetInterface {

    	/// POD for MP2 keywords
        struct Parameters {
        	double thresh_;			///< the accuracy threshold
        	double dconv_;			///< threshold for the MP1 residual
        	int i; 					///< electron 1, used only if a specific pair is requested
        	int j; 					///< electron 2, used only if a specific pair is requested

        	/// number of frozen orbitals; note the difference to the "pair" keyword where you
        	/// request a specific orbital. Here you freeze lowest orbitals, i.e. if you find
        	///  freeze 1
        	/// in the input file the 0th orbital is kept frozen, and orbital 1 is the first to
        	/// be correlated.
        	int freeze;

        	/// the restart flag initiates the loading of the pair functions

        	/// if this flag is set the program expect for each pair a file named
        	///  pair_ij.00000
        	/// where ij is to be replaced by the values of i and j.
        	/// These files contain the restart information for each pair.
        	bool restart;

        	/// maximum number of subspace vectors in KAIN
        	int maxsub;

        	/// ctor reading out the input file
        	Parameters(const std::string& input) : thresh_(-1.0), dconv_(-1.0),
        			i(-1), j(-1), freeze(0), restart(false), maxsub(2) {

        		// get the parameters from the input file
                std::ifstream f(input.c_str());
                position_stream(f, "mp2");
                std::string s;

                while (f >> s) {
                    if (s == "end") break;
                    else if (s == "econv") f >> thresh_;
                    else if (s == "dconv") f >> dconv_;
                    else if (s == "pair") f >> i >> j;
                    else if (s == "maxsub") f >> maxsub;
                    else if (s == "freeze") f >> freeze;
                    else if (s == "restart") restart=true;
                    else continue;
                }
                // set default for dconv if not explicitly given
                if (dconv_<0.0) dconv_=sqrt(thresh_)*0.1;
        	}

            /// check the user input
        	void check_input(const std::shared_ptr<HartreeFock> hf) const {
                if (freeze>hf->nocc()) MADNESS_EXCEPTION("you froze more orbitals than you have",1);
                if (i>=hf->nocc()) MADNESS_EXCEPTION("there is no i-th orbital",1);
                if (j>=hf->nocc()) MADNESS_EXCEPTION("there is no j-th orbital",1);
                if (thresh_<0.0) MADNESS_EXCEPTION("please provide the accuracy threshold for MP2",1);
        	}
        };


        World& world;                           ///< the world
        Parameters param;						///< SCF parameters for MP2
        std::shared_ptr<HartreeFock> hf;        ///< our reference
        CorrelationFactor corrfac;              ///< correlation factor: Slater
        std::shared_ptr<NuclearCorrelationFactor> nuclear_corrfac;

        std::map<std::pair<int,int>,ElectronPair> pairs;       ///< pair functions and energies
        double correlation_energy;				///< the correlation energy
        double coords_sum;						///< check sum for the geometry

        StrongOrthogonalityProjector<double,3> Q12;

    private:
        struct Intermediates {
            std::string function;      ///< pair function for a specific pair w/o correlation factor part
            std::string r12phi;        ///< orbital product multiplied with the correlation factor
            std::string Kfphi0;        ///< the function K f12 |phi^0>
            std::string Uphi0;         ///< the function U |phi^0>  (U being Kutzelnigg's potential)
            std::string KffKphi0;      ///< the function [K,f12] |phi^0>
            std::string OUKphi0;		///< < k(1) | U-K | phi^0(1,2) >

            Intermediates() : r12phi(), Kfphi0(), Uphi0(), KffKphi0() {};

            Intermediates(World& world, const std::string& filename) : function(), r12phi(),
                    Kfphi0(), Uphi0(), KffKphi0() {
                std::ifstream f(filename.c_str());
                position_stream(f, "mp2");
                std::string s;

                while (f >> s) {
                    if (s == "end") break;
                    else if (s == "function") f >> function;
                    else if (s == "r12phi") f >> r12phi;
                    else if (s == "Kfphi0") f >> Kfphi0;
                    else if (s == "Uphi0") f >> Uphi0;
                    else if (s == "KffKphi0") f >> KffKphi0;
                    else {continue;
                    }
                    if (world.rank()==0) print("found intermediate in control file: ",s);
                }
            }

            template <typename Archive> void serialize (Archive& ar) {
                ar & function & r12phi & Kfphi0 & Uphi0 & KffKphi0;
            }

        };


        Intermediates intermediates;
        std::shared_ptr<real_convolution_3d> poisson;

    public:

        /// ctor
        MP2(World& world, const std::string& input)
            : world(world)
        	, param(input.c_str())
    		, corrfac(world)
        	, correlation_energy(0.0)
    		, coords_sum(-1.0)
    		, Q12(world)
            , intermediates() {

        	{
        		std::shared_ptr<SCF> calc=std::shared_ptr<SCF>
        				(new SCF(world,input.c_str()));

        		// get parameters form input file for hf
        		if (world.rank()==0) print("accuracy from dft will be overriden by mp2 to 0.01*thresh");
                calc->set_protocol<6>(world,thresh());
                calc->param.econv=thresh()*0.01;
                calc->set_protocol<3>(world,calc->param.econv);

        		// override computed parameters if they are provided explicitly
                double eprec=calc->param.econv*0.1;

                std::ifstream f(input.c_str());
                position_stream(f, "geometry");
                std::string s;
                while (std::getline(f,s)) {
                    std::istringstream ss(s);
                    std::string tag;
                    ss >> tag;
                    if (tag == "end") break;
                    else if (tag == "eprec") ss >> eprec;
                    else continue;
                }

                calc->molecule.set_eprec(eprec);
                if (world.rank()==0) calc->molecule.print();

                hf=std::shared_ptr<HartreeFock>(new HartreeFock(world,calc));
                poisson=std::shared_ptr<real_convolution_3d>
                	(CoulombOperatorPtr(world,0.0001,calc->param.econv));

                corrfac=CorrelationFactor(world,1.0,dcut,calc->molecule);
                nuclear_corrfac=hf->nemo_calc.nuclear_correlation;

        	}

        	// print some output for the user
            if (world.rank()==0) {
            	hf->get_calc().param.print(world);
            	this->print_info(world);
            }

            // check the user input for consistency
            param.check_input(hf);

            // initialize the electron pairs and store restart info
            // or load previous restart information
    		for (int i=param.freeze; i<hf->nocc(); ++i) {
    			for (int j=i; j<hf->nocc(); ++j) {
    	        	std::pair<int,int> key=std::make_pair(i,j);
    	        	pairs.insert(std::make_pair(key,ElectronPair(i,j)));
    				if (param.restart) pair(i,j).load_pair(world);
				}
			}

            if (world.rank()==0) intermediates=Intermediates(world,input);
            world.gop.broadcast_serializable(intermediates, 0);

        }

        /// return a reference to the electron pair for electrons i and j

        /// @param[in]	i	index for electron 1
        /// @param[in]	j	index for electron 2
        /// @return		reference to the electron pair ij
        ElectronPair& pair(const int i, const int j) {
        	// since we return a reference the keyval must already exist in the map
        	MADNESS_ASSERT(pairs.find(std::make_pair(i,j)) != pairs.end());
        	return pairs.find(std::make_pair(i,j))->second;
        }

        /// return a checksum for the geometry
        double coord_chksum() const {return coords_sum;}

        /// return the molecular correlation energy energy (without the HF energy)
        double value() {
        	hf->value();		// make sure the reference is converged

            return value(hf->get_calc().molecule.get_all_coords());
        }

        /// return the molecular correlation energy as a function of the coordinates
        double value(const Tensor<double>& x)  {

        	// fast return if the MP2 energy is already solved at this geometry
            double xsq = x.sumsq();
            if (xsq == coord_chksum()) return correlation_energy;

            // make sure HF used the same geometry as we do
            coords_sum=xsq;
            MADNESS_ASSERT(std::fabs(coord_chksum()-hf->coord_chksum())<1.e-14);

            // set all orbitals spaces
        	// When a nuclear correlation factor is used the residual equations
        	// are similarity transformed. Therefore the orbitals in the
        	// projection operator must be set accordingly.
        	if (nuclear_corrfac->type()==NuclearCorrelationFactor::None) {
                Q12.set_spaces(hf->get_calc().amo);
//        	} else if (nuclear_corrfac->type()==NuclearCorrelationFactor::GaussSlater) {
        	} else {
        		// only valid for closed shell
                MADNESS_ASSERT(hf->get_calc().param.spin_restricted);
                const std::vector<real_function_3d>& nemos=hf->nemos();
        		const std::vector<real_function_3d>& R2amo=hf->R2orbitals();
                Q12.set_spaces(R2amo,nemos,R2amo,nemos);
                if (world.rank()==0) {
                	print("set orbital spaces for the SO projector");
                	print("Q12,R = (1-|nemo><nemo|R2) (1-|nemo><nemo|R2)");
                }
//        	} else {
//        		MADNESS_EXCEPTION("unknown nuclear correlation factor",1);
        	}


        	correlation_energy=0.0;

        	// compute only one single pair
        	if ((param.i>-1) and (param.j>-1)) {
                pair(param.i,param.j)=solve_residual_equations(param.i,param.j);
                correlation_energy+=pair(param.i,param.j).e_singlet+pair(param.i,param.j).e_triplet;

            // solve the residual equations for all pairs ij
        	} else {
        		for (int i=param.freeze; i<hf->nocc(); ++i) {
        			for (int j=i; j<hf->nocc(); ++j) {
        				if (pair(i,j).converged) {
        					pair(i,j).print_energy();
							correlation_energy+=pair(i,j).e_singlet+pair(i,j).e_triplet;
							make_Rpsi(pair(i,j));
        				} else {
							pair(i,j)=solve_residual_equations(i,j);
							correlation_energy+=pair(i,j).e_singlet+pair(i,j).e_triplet;
        				}
					}
				}
        	}
            return correlation_energy;
        }

        /// print the SCF parameters
        void print_info(World& world) const {
            if (world.rank()==0) {
                madness::print("\n === MP2 info === \n");
                madness::print("         MP2 restart ", param.restart);
                madness::print("        threshold 3D ", FunctionDefaults<3>::get_thresh());
                madness::print("        threshold 6D ", FunctionDefaults<6>::get_thresh());
                madness::print("    threshold energy ", param.thresh_);
                madness::print("  threshold residual ", param.dconv_);
                madness::print("     truncation mode ", FunctionDefaults<6>::get_truncate_mode());
                madness::print("         tensor type ", FunctionDefaults<6>::get_tensor_type());
                madness::print("           facReduce ", GenTensor<double>::fac_reduce());
                madness::print("    max displacement ", Displacements<6>::bmax_default());
                madness::print("     apply randomize ", FunctionDefaults<6>::get_apply_randomize());
            	if (param.i>=0 and param.j>=0) {
            		madness::print("      computing pair ", param.i,param.j);
            	} else if (param.i<0 and param.j<0) {
            		if (param.freeze==0) madness::print("   # frozen orbitals ","none");
            		if (param.freeze>0) madness::print("   # frozen orbitals ",0, " to ",param.freeze-1);
					madness::print(" correlated orbitals ", param.freeze," to ",hf->nocc()-1);
            	}
        		madness::print("   max KAIN subspace ", param.maxsub);
            }
        }

        /// return the underlying HF reference
        HartreeFock& get_hf() {return *hf;}

        /// return the accuracy
        double thresh() const {return param.thresh_;}

        /// return the 0th order energy of pair ij (= sum of orbital energies)
        double zeroth_order_energy(const int i, const int j) const {
            return hf->orbital_energy(i)+hf->orbital_energy(j);
        }

        /// solve the residual equation for electron pair (i,j)
        ElectronPair solve_residual_equations(const int i, const int j) {

        	if (world.rank()==0) printf("\n\nsolving electron pair (%d, %d)\n\n",i,j);
            ElectronPair result=guess_mp1_3(i,j);
            result.store_pair(world);

            double energy=compute_energy(result);
            if (world.rank()==0) printf("finished with prep step at time %6.1fs with energy %12.8f\n\n",
            		wall_time(),energy);

            // the Green's function depends on the zeroth order energy, which is the sum
            // of the orbital energies of orbitals i and j
            //  -2.0 G = (T - e_i - e_j) ^ -1
            const double eps=zeroth_order_energy(i,j);
            real_convolution_6d green = BSHOperator<6>(world, sqrt(-2*eps), lo, bsh_eps);
            green.destructive()=true;
//            real_convolution_6d green = BSHOperator<6>(world, sqrt(-2*eps), 0.00001,
//            		FunctionDefaults<6>::get_thresh()*0.1);

            // do some preiterations using increments
//            increment(result,green);

			if (world.rank()==0) print("computing iteratively");
//			real_function_6d constant_term;
//			load_function(constant_term,"GVpair");

			NonlinearSolverND<6> solver(param.maxsub);
			// increment iteration counter upon entry
			for (++result.iteration; result.iteration<20; ++result.iteration) {

				// apply the convolution
				real_function_6d vphi=multiply_with_0th_order_Hamiltonian(result.function,i,j);
				vphi.scale(-2.0).truncate();
				load_balance(vphi,false);

				real_function_6d tmp=green(vphi);	// green is destructive
				tmp.print_size("GV psi");

				// we have to solve this equation:
				// psi1 = psi0 + GVpsi1 <=> psi0 + GVpsi1 - psi1 = r =0
				tmp=(Q12(result.constant_term + tmp)).truncate();
//				const double old_fnorm=result.function.norm2();

				real_function_6d residual=result.function - tmp;
				result.function=Q12(solver.update(tmp,residual));
				const double rnorm=residual.norm2();
				const double fnorm=result.function.norm2();
				if (world.rank()==0)
					printf("norm2 of psi, residual %12.8f %12.8f\n",fnorm,rnorm);

				double old_energy=energy;
				energy=compute_energy(result);

				result.converged=((std::abs(old_energy-energy)<result.function.thresh()*0.01)
						and (rnorm<param.dconv_));
				result.store_pair(world);
				save_function(result.function,"pair_iter"+stringify(result.iteration));

				if (world.rank()==0) printf("finished iteration %2d at time %8.1fs with energy %12.8f\n\n",
						result.iteration, wall_time(),energy);

				if (result.converged) break;

			}

            // save the converged first order pair function separately for easier access
        	std::string name="pair_"+stringify(i)+stringify(j)+"_psi1_converged";
            save_function(result.function,name);

            // print the final pair energies
            result.print_energy();
            return result;
        }

        real_function_6d make_Rpsi(const ElectronPair& pair) const {
            const real_function_3d R=hf->nemo_calc.R;
        	real_function_6d Rpair1=multiply(pair.function,R,1).truncate();
        	real_function_6d Rpair=multiply(Rpair1,R,2).truncate();
        	int i=pair.i, j=pair.j;
        	std::string name1="pair_"+stringify(i)+stringify(j)+"_Rpsi1_converged";
            save_function(Rpair,name1);
            return Rpair;
        }

		/// compute increments: psi^1 = C + GV C + GVGV C + GVGVGV C + ..

        /// for pre-optimization of the mp1 wave function
    	/// no restart option here for now
        /// param[inout]	pair	the electron pair
        /// param[in]		green	the Green's function
        void increment(ElectronPair& pair, real_convolution_6d& green) {

        	// don't mess up if we've already done some iterations!
        	if (pair.iteration>0) return;

			if (world.rank()==0) print("computing increments");
			real_function_6d latest_increment;
			load_function(latest_increment,"GVpair");

			for (int ii=1; ii<20; ++ii) {

				real_function_6d vphi=multiply_with_0th_order_Hamiltonian(latest_increment,pair.i,pair.j);
				load_balance(vphi,false);

				/// apply the convolution
				vphi.scale(-2.0).truncate();
				latest_increment=green(vphi).truncate();		// green is destructive
				latest_increment.print_size("result of applying 0th order Hamiltonian on latest increment");

				latest_increment=Q12(latest_increment);
				pair.function=pair.function+latest_increment;

				// compute the energy
				compute_energy(pair);

				if (world.rank()==0) printf("finished increment %2d at time %.1fs\n\n", ii, wall_time());
				const double residual_norm=latest_increment.norm2();
				if (residual_norm<pair.function.thresh()*0.01) break;
			}
        }

        /// swap particles 1 and 2

        /// param[in]	f	a function of 2 particles f(1,2)
        /// return	the input function with particles swapped g(1,2) = f(2,1)
        real_function_6d swap_particles(const real_function_6d& f) const {

        	// this could be done more efficiently for SVD, but it works decently
        	std::vector<long> map(6);
        	map[0]=3; map[1]=4; map[2]=5;	// 2 -> 1
        	map[3]=0; map[4]=1; map[5]=2;	// 1 -> 2
        	return mapdim(f,map);
        }

        double asymmetry(const real_function_6d& f, const std::string s) const {
        	return 0.0;
        	const real_function_6d ff=swap_particles(f);
        	double diff=(ff-f).norm2();
        	f.check_symmetry();
        	if (world.rank()==0) print("asymmetry of ",s,"  ",diff);
        	return diff;
        }

        void test(const std::string filename) {
        	// compute the singlet pair function and save to file
        	real_function_6d plain_pair;
        	load_function(plain_pair,filename);

        	{
				real_function_6d plain_pair_swapped=swap_particles(plain_pair);

				real_function_6d singlet=(plain_pair+plain_pair_swapped).truncate().reduce_rank();
				save_function(singlet,"singlet"+filename);

				real_function_6d triplet=(plain_pair-plain_pair_swapped).truncate().reduce_rank();
				save_function(triplet,"triplet"+filename);
        	}

        	hf->value();
        	const real_function_3d R=hf->nemo_calc.R;
        	real_function_6d Rpair1=multiply(plain_pair,R,1).truncate();
        	plain_pair=multiply(Rpair1,R,2).truncate();
        	{
				real_function_6d plain_pair_swapped=swap_particles(plain_pair);

				real_function_6d singlet=(plain_pair+plain_pair_swapped).truncate().reduce_rank();
				save_function(singlet,"Rsinglet"+filename);

				real_function_6d triplet=(plain_pair-plain_pair_swapped).truncate().reduce_rank();
				save_function(triplet,"Rtriplet"+filename);
        	}

        }
		
        /// compute the matrix element <ij | g12 Q12 f12 | phi^0>

        /// scales quartically. I think I can get this down to cubically by
        /// setting Q12 = (1 - O1)(1 - O2) = 1- O1(1 - 0.5 O2) - O2 (1 - 0.5 O1)
        /// as for the formulas cf the article mra_molecule
        /// @return 	the energy <ij | g Q f | kl>
        double compute_gQf(const int i, const int j, ElectronPair& pair) const {

        	// for clarity of notation
        	const int k=pair.i;
        	const int l=pair.j;

        	// the ket space
            const real_function_3d& ket_i=hf->nemo(i);
            const real_function_3d& ket_j=hf->nemo(j);

            // the bra space
            const real_function_3d& bra_k=hf->R2orbital(k);
            const real_function_3d& bra_l=hf->R2orbital(l);

            // compute <ij| fg |kl>: do it in 3D as (ik| fg |jl)
            // the operator fg can be rewritten as 1/r12 - f/r12
            // i.e. as poisson kernel and a bsh kernel. Note the
            // the bsh kernel includes a factor of 1/(4 pi)
            const real_function_3d ik=ket_i*bra_k;
            const real_function_3d jl=ket_j*bra_l;

            // make all the operators that we need
            const double fourpi=4.0*constants::pi;
            real_convolution_3d fg = BSHOperator<3>(world, corrfac.gamma(), lo, bsh_eps/fourpi);
            real_convolution_3d gg = CoulombOperator(world,lo,bsh_eps);
            real_convolution_3d slaterf12 = SlaterF12Operator(world, corrfac.gamma(), lo, bsh_eps/fourpi);

            //  < ij | fg | kl >
            const real_function_3d ik_fg=(gg)(ik) - fourpi*fg(ik);
            const double a=inner(ik_fg,jl)/(2.0*corrfac.gamma());
            if (world.rank()==0) printf("<%d%d | f/r              | %d%d>  %12.8f\n",i,j,k,l,a);


            // compute <ij| f (O1 + O2) g | ij>

            // compute bra space xi(ik,j)^dagger, i.e. the hermitian conjugate of xi
            // the index k is implicit in the vector of functions
            // naming: xi _ orbitals _ operator _ hc
            std::vector<real_function_3d> xi_ij_g_ket=make_xi(ket_i,ket_j,*poisson,false);       // xi_{i,m*},j
            std::vector<real_function_3d> xi_ji_g_ket=make_xi(ket_j,ket_i,*poisson,false);       // xi_{j,m*},i

            std::vector<real_function_3d> xi_ij_f_bra=make_xi(bra_k,bra_l,slaterf12,true);       // xi_{i*,m},j*
            std::vector<real_function_3d> xi_ji_f_bra=make_xi(bra_l,bra_k,slaterf12,true);       // xi_{j*,m},i*

            // in the following do NOT use antisymmetrized pair functions:
            // |ij> -> 0.5 |ij - ji>

            // < ij | f12 O1 g12 | kl >
            //   = \sum_m <i(1) j(2) | f12 | m(1) >< m(3) | g23 | k(3) l(2)>
            //   = \sum_m < chi^f_i*,m(2) j*(2) | chi^g_k,m*(2) l(2) >
            //   = \sum_m < xi^f_im,j | xi^g_km,l >
            const double o1a=inner(world,xi_ij_f_bra,xi_ij_g_ket).sum();
            if (world.rank()==0) printf("<%d%d | f12 O1 g12       | %d%d>  %12.8f\n",i,j,k,l,o1a);


            // < ij | f12 O2 g12 | kl >
            //    = \sum_m <i(1) j(2) | f12 | m(2) >< m(3) | g13 | k(1) l(3)>
            //    = \sum_m <chi^f_j*,m(1) i*(1) | chi^g_l,m*(1) k(1) >
            //    = \sum_m < xi^f_jm,i | xi^g_lm,k >
            const double o2a=inner(world,xi_ji_f_bra,xi_ji_g_ket).sum();
            if (world.rank()==0) printf("<%d%d | f12 O2 g12       | %d%d>  %12.8f\n",i,j,k,l,o2a);


            // compute <ij| f O1 O2 g | kl>  // why do I need to swap ij in g_ijkl??
            const Tensor<double> f_ijmn = matrix_inner(world,xi_ij_f_bra,hf->nemos());
            const Tensor<double> g_ijmn = matrix_inner(world,hf->R2orbitals(),xi_ji_g_ket);
            const double o12=f_ijmn.trace(g_ijmn);
            if (world.rank()==0) printf("<%d%d | f12 O12 g12      | %d%d>  %12.8f\n",i,j,k,l,o12);


            const double e=a-o1a-o2a+o12;
            if (world.rank()==0) printf("<%d%d | g Q12 f          | %d%d>  %12.8f\n",i,j,k,l,e);

            return e;
        }

    private:

        /// save a function
        template<typename T, size_t NDIM>
        void save_function(const Function<T,NDIM>& f, const std::string name) const {
            if (world.rank()==0) print("saving function",name);
            f.print_size(name);
            archive::ParallelOutputArchive ar(world, name.c_str(), 1);
            ar & f;
        }

        /// load a function
        template<typename T, size_t NDIM>
        void load_function(Function<T,NDIM>& f, const std::string name) const {
            if (world.rank()==0) print("loading function",name);
            archive::ParallelInputArchive ar(world, name.c_str());
            ar & f;
            f.print_size(name);
        }

        /// return the function Uphi0; load from disk if available
        real_function_6d make_Uphi0(ElectronPair& pair) const {
        	const int i=pair.i;
        	const int j=pair.j;
        	real_function_6d Uphi0=real_factory_6d(world);
            if (not intermediates.Uphi0.empty()) {
                load_function(Uphi0,intermediates.Uphi0);
            } else {
            	// apply the pure commutator of the kinetic energy and the
            	// electronic correlation factor on the (regularized)
            	// orbitals i and j
            	//  [T,f]
            	const double eps=this->zeroth_order_energy(i,j);
                Uphi0=corrfac.apply_U(hf->nemo(i),hf->nemo(j),eps);
                {
                	real_function_6d tmp=CompositeFactory<double,6,3>(world)
                                     .particle1(copy(hf->R2orbital(i)))
                                     .particle2(copy(hf->R2orbital(j)));

                	const double a=inner(Uphi0,tmp);
                	if (world.rank()==0) print("< nemo | R2 U | nemo>",a);
                }
    			asymmetry(Uphi0,"Uphi w/o R");

                // apply the mixed commutator of the kinetic energy and the
                // electronic and nuclear correlation factor on the regularized
                // orbitals i and j:
                //  R^{-1} [[T,f],R]
                // Vanishes if no nuclear correlation factor
                // is used, since [[T,f],1] = 0
    			if (1) {
//                if (hf->nemo_calc.nuclear_correlation->type()!=NuclearCorrelationFactor::None) {
                	if (world.rank()==0) print("adding the mixed commutator R^{-1} [[T,f],R]");
					for (int axis=0; axis<3; axis++) {

						double tight_thresh=std::min(thresh(),1.e-4);

				        real_convolution_6d op_mod = BSHOperator<6>(world, sqrt(-2*eps),0.0001,1e-5);
				        op_mod.modified()=true;

						const real_function_3d u1_nuc=hf->nemo_calc.nuclear_correlation->U1(axis);
						const real_function_3d u1_nuc_nemo_i=u1_nuc*hf->nemo(i);
						const real_function_3d u1_nuc_nemo_j=u1_nuc*hf->nemo(j);
						const real_function_6d u1_el=corrfac.U1(axis);
//						u1_nuc.print_size("u1_nuc");
//						plot_plane(world,u1_nuc,"u1_nuc");
//						u1_nuc_nemo_i.print_size("u1_nuc_nemo_i");
//						plot_plane(world,u1_nuc_nemo_i,"u1_nuc_nemo_i");
//						plot_plane(world,hf->nemo(i),"nemo_i");
//						u1_nuc_nemo_j.print_size("u1_nuc_nemo_j");
//						u1_el.print_size("u1_el");

						real_function_6d U1_mix1=CompositeFactory<double,6,3>(world)
									.g12(u1_el)
									.particle1(copy(u1_nuc_nemo_i))
									.particle2(copy(hf->nemo(j)))
									.thresh(tight_thresh);
						U1_mix1.fill_tree(op_mod);
						U1_mix1.set_thresh(thresh());
						U1_mix1.print_size("U1_mix1");
//						plot_plane(world,U1_mix1,"U1_mix1");

						real_function_6d U1_mix2=CompositeFactory<double,6,3>(world)
									.g12(u1_el)
									.particle1(copy(hf->nemo(i)))
									.particle2(copy(u1_nuc_nemo_j))
									.thresh(tight_thresh);
						U1_mix2.fill_tree(op_mod);
						U1_mix2.set_thresh(thresh());
						U1_mix2.print_size("U1_mix2");
//						plot_plane(world,U1_mix2,"U1_mix2");

						real_function_6d diff=(U1_mix1-U1_mix2).scale(-1.0);
						diff.truncate();

						// multiply with -1 from the kinetic energy operator
//						Uphi0=(Uphi0+(U1_mix1-U1_mix2).scale(-1.0)).truncate();
						Uphi0=(Uphi0+diff).truncate();
		    			this->asymmetry(Uphi0,"Uphi0 in R");

//						Uphi0.print_size("Uphi0");
//						plot_plane(world,Uphi0,"Uphi0");

		                {
		                	real_function_6d tmp=CompositeFactory<double,6,3>(world)
		                                     .particle1(copy(hf->R2orbital(i)))
		                                     .particle2(copy(hf->R2orbital(j)));

		                	const double a=inner(Uphi0,tmp);
		                	if (world.rank()==0) print("< nemo | R2 U | nemo>",a);
		                }

					}
                }
    			asymmetry(Uphi0,"Uphi0");

                // save the function for restart
//                save_function(Uphi0,"Uphi0");
            }

            // sanity check: <ij| [T,g12] |ij> = <ij | U |ij> - <ij| g12 | ij> = 0
            real_function_6d tmp=CompositeFactory<double,6,3>(world)
                                 .particle1(copy(hf->R2orbital(i)))
                                 .particle2(copy(hf->R2orbital(j)));

            const double a=inner(Uphi0,tmp);
            const real_function_3d ii=hf->R2orbital(i)*hf->nemo(i);
            const real_function_3d jj=hf->R2orbital(j)*hf->nemo(j);
            const real_function_3d gii=(*poisson)(ii);
            const double aa=inner(jj,gii);
            const double error=std::fabs(a-aa);
            if (world.rank()==0) {
            	printf("< phi0 | U_R   | phi0 >  %12.8f\n",a);
            	printf("< phi0 | 1/r12 | phi0 >  %12.8f\n",aa);
                if (error>thresh()) print("WARNING : Kutzelnigg's potential inaccurate");
                if (error>thresh()*10.0) MADNESS_EXCEPTION("Kutzelnigg's potential plain wrong",1);
            }
            Uphi0.print_size("Uphi0");
            return Uphi0;
        }

        /// return the function [K,f] phi0; load from disk if available
        real_function_6d make_KffKphi0(const ElectronPair& pair) const {

        	const int i=pair.i;
        	const int j=pair.j;

        	real_function_6d KffKphi0;
			if (not intermediates.KffKphi0.empty()) {
				load_function(KffKphi0,intermediates.KffKphi0);

			} else {
                real_function_6d r12nemo=CompositeFactory<double,6,3>(world).g12(corrfac.f())
                         .particle1(copy(hf->nemo(i))).particle2(copy(hf->nemo(j)));
                r12nemo.fill_tree().truncate().reduce_rank();
                r12nemo.print_size("r12nemo");
                save_function(r12nemo,"r12nemo");

				real_function_6d Kfphi0;
				if (not intermediates.Kfphi0.empty()) {
					load_function(Kfphi0,intermediates.Kfphi0);
				} else {
					Kfphi0=K(r12nemo,i==j);

//					const real_function_3d& phi_i=hf->nemo(i);
//					const real_function_3d& phi_j=hf->nemo(j);
//
//					real_convolution_3d op=CoulombOperator(world,0.0001,hf->get_calc().param.econv);
//					op.particle()=1;
//
//					real_convolution_3d op_mod=CoulombOperator(world,0.0001,hf->get_calc().param.econv);
//	                op_mod.modified()=true;
//
//					real_function_6d result=real_factory_6d(world);
//					for (int k=0; k<hf->nocc(); ++k) {
//						const real_function_3d& phi_k_bra=hf->R2orbital(k);
//						const real_function_3d& phi_k_ket=hf->nemo(k);
//						real_function_6d f_ijk=CompositeFactory<double,6,3>(world)
//								.g12(corrfac.f())
//								.particle1(copy(phi_i*phi_k_bra))
//								.particle2(copy(phi_j));
//						f_ijk.fill_tree(op_mod).truncate();
//						real_function_6d x=op(f_ijk).truncate();
//			            result+=multiply(copy(x),copy(phi_k_ket),1).truncate();
//					}
//
//					if (i==j) {
//						result+=swap_particles(result);
//					} else {
//						op.particle()=2;
//						for (int k=0; k<hf->nocc(); ++k) {
//							const real_function_3d& phi_k_bra=hf->R2orbital(k);
//							const real_function_3d& phi_k_ket=hf->nemo(k);
//							real_function_6d f_ijk=CompositeFactory<double,6,3>(world)
//									.g12(corrfac.f())
//									.particle1(copy(phi_i))
//									.particle2(copy(phi_j*phi_k_bra));
//							f_ijk.fill_tree(op_mod).truncate();
//							real_function_6d x=op(f_ijk).truncate();
//				            result+=multiply(copy(x),copy(phi_k_ket),2).truncate();
//						}
//					}
//					Kfphi0=result;
				}

				{
					real_function_6d tmp=CompositeFactory<double,6,3>(world)
								.particle1(copy(hf->R2orbitals()[i]))
								.particle2(copy(hf->R2orbitals()[j]));
					double a1=inner(Kfphi0,tmp);
					if (world.rank()==0) printf("< nemo0 | R^2 R-1 K f R | nemo0 >  %12.8f\n",a1);
//					save_function(Kfphi0,"Kfphi0");
				}
				const real_function_6d fKphi0=make_fKphi0(pair.i,pair.j);
				{
					real_function_6d tmp=CompositeFactory<double,6,3>(world)
								.particle1(copy(hf->R2orbitals()[i]))
								.particle2(copy(hf->R2orbitals()[j]));
					double a2=inner(fKphi0,tmp);
					if (world.rank()==0) printf("< nemo0 | R^2 R-1 f K R | nemo0 >  %12.8f\n",a2);
				}
	            KffKphi0=(Kfphi0-fKphi0).truncate().reduce_rank();
//				save_function(KffKphi0,"KffKphi0");
			}

			// sanity check
			real_function_6d tmp=CompositeFactory<double,6,3>(world)
						.particle1(copy(hf->R2orbitals()[i]))
						.particle2(copy(hf->R2orbitals()[j]));
            const double a=inner(KffKphi0,tmp);
            if (world.rank()==0) {
            	printf("< nemo0 | R^2 R-1 [K,f] R | nemo0 >  %12.8f\n",a);
                if (std::fabs(a)>thresh()) print("WARNING : exchange commutator inaccurate");
                if (std::fabs(a)>thresh()*10.0) MADNESS_EXCEPTION("exchange commutator plain wrong",1);
            }
            KffKphi0.print_size("KffKphi0");
			return KffKphi0;
        }

        /// compute some matrix elements that don't change during the SCF
        ElectronPair make_pair(const int i, const int j) const {

            ElectronPair p=ElectronPair(i,j);
            p.load_pair(world);

//            // some functions repeatedly used
//            if (intermediates.r12phi.empty()) {
//                p.r12phi=CompositeFactory<double,6,3>(world)
//                            .g12(corrfac.f())
//                            .particle1(copy(hf->orbital(i)))
//                            .particle2(copy(hf->orbital(j)));
//                p.r12phi.fill_tree().truncate().reduce_rank();
//                save_function(p.r12phi,"r12phi");
//            } else {
//            	load_function(p.r12phi,intermediates.r12phi);
//            }

            // compute and store them if they have not been read from disk
            if (p.ij_gQf_ij==ElectronPair::uninitialized()) {
				p.ij_gQf_ij=compute_gQf(i,j,p);
				p.ji_gQf_ij=0.0;
				if (i!=j) p.ji_gQf_ij=compute_gQf(j,i,p);
				p.store_pair(world);
            }

            if (world.rank()==0) printf("done with matrix elements at time %.1fs\n\n", wall_time());
            return p;
        }

        /// compute the first iteration of the residual equations and all intermediates
        ElectronPair guess_mp1_3(const int i, const int j) const {

            // the Green's function
            const double eps=zeroth_order_energy(i,j);
            real_convolution_6d green = BSHOperator<6>(world, sqrt(-2.0*eps),lo,bsh_eps);
//            real_convolution_6d green = BSHOperator<6>(world, sqrt(-2.0*eps), 0.00001,
//            		FunctionDefaults<6>::get_thresh()*0.1);

            ElectronPair pair=make_pair(i,j);

            // fast return if possible
            if (pair.function.is_initialized()) return pair;

            pair.Uphi0=make_Uphi0(pair);
			pair.KffKphi0=make_KffKphi0(pair);

            // these are the terms that come from the single projectors: (O1 + O2) (U+[K,f])|phi^0>
			OUKphi0(pair);

			// make the terms with high ranks and smallish trees
			load_balance(pair.Uphi0,true);
			real_function_6d Vpair1=(pair.Uphi0-pair.KffKphi0).truncate().reduce_rank();
			pair.Uphi0.clear();
			pair.KffKphi0.clear();
			asymmetry(Vpair1,"Vpair");
			load_balance(Vpair1,false);
			real_function_6d GVpair;
            green.destructive()=true;			// green will destroy Vpair1
			GVpair=green(-2.0*Vpair1).truncate().reduce_rank();
			Vpair1.clear(true);
        	if (world.rank()==0) print("symmetry GVpair1");
			save_function(GVpair,"GVpair1");
//			load_function(GVpair,"GVpair1");
			asymmetry(GVpair,"GVpair1");

			// make the terms with low ranks and largish trees:
			// - G ( O1 + O2 - O1O2) (U-K) | phi0 >
			// With the notations
			//  | UK_k(1) > = < k(2) | U-K | phi0 >		// these have been computed in OUKphi0()
			//  h(kl) = <kl|[H,f] + 1/r12| phi0>
			// we get
			// G ( O1 + O2 - O1O2) (U-K) | phi0 >
			//  = |k(1)>|UK_k(2)> - |l(2)>|UK_l(1)> + |k(1)l(2)> (<kl|[H,f] + 1/r12| phi0> )
			//  = |k(1)> ( |UK_k(2)> - 1/2 |l(2)> h(kl) ) - |l(2)> ( |UK_l(1) - 1/2 |k(1)> h(kl) )
			// which scales cubicly but for the computation of h(kl)

			// compute the matrix elements h(kl) = <k(1) l(2) | [H,f] + 1/r12 | phi0 >
			Tensor<double> h(hf->nocc(),hf->nocc());
            const double fourpi=4.0*constants::pi;
            real_convolution_3d slaterf12 = SlaterF12Operator(world, corrfac.gamma(), lo, bsh_eps/fourpi);

			for (int k=0; k<hf->nocc(); ++k) {

	            // for the matrix element <ij|kl> = (ik|jl)
	            const real_function_3d ik=hf->orbital(i)*hf->orbital(k);
	            const real_function_3d gik=(*poisson)(ik).truncate();
	            const real_function_3d fik=slaterf12(ik).truncate();

	            // the function Of(2) = <k(1) | f12 | phi0 >
//	            const real_function_3d Of=pair.r12phi.project_out(hf->orbital(k),0);

				for (int l=0; l<hf->nocc(); ++l) {

					const real_function_3d jl=hf->orbital(j)*hf->orbital(l);
		            const double g_ijkl=inner(jl,gik);
		            const double f_ijkl=inner(jl,fik);

		            // the matrix element <kl | f12 | ij>
//		            const double f_ijkl=inner(hf->orbital(l),Of);
                    const double e_kl=zeroth_order_energy(k,l);
                    const double e_ij=zeroth_order_energy(i,j);

                    h(k,l)=f_ijkl*(e_kl-e_ij)+g_ijkl;
				}
			}
//			if (world.rank()==0) print("h\n",h);
//			for (int k=0; k<hf->nocc(); ++k) {
//				for (int l=0; l<hf->nocc(); ++l) {
//		            real_function_6d bra=CompositeFactory<double,6,3>(world)
//								.particle1(copy(hf->R2orbital(k)))
//								.particle2(copy(hf->R2orbital(l)));
//		            h(k,l)=inner(Vpair1,bra);
//				}
//			}
//			if (world.rank()==0) print("h = <nemo |R2 (U_R -K_R) |nemo> \n",h);

			// the result function; adding many functions requires tighter threshold
			const double thresh=FunctionDefaults<6>::get_thresh();
			const double tight_thresh=thresh*0.1;
			FunctionDefaults<6>::set_thresh(tight_thresh);
        	real_function_6d r1=real_factory_6d(world).thresh(tight_thresh);
        	for (int k=0; k<hf->nocc(); ++k) {

        		// make the term  tmp2(2) = |UK_k(2)> - 1/2 |l(2)> h(kl)
        		real_function_3d tmp2=pair.phi_k_UK_phi0[k];
        		for (int l=0; l<hf->nocc(); ++l) tmp2-= 0.5*h(k,l)*hf->nemo(l);

        		// now apply the Greens' function (including the -2.0 factor in one term)
				real_function_6d tmp=green(-2.0*hf->nemo(k),tmp2);
				tmp.set_thresh(tight_thresh);
				r1=(r1-tmp).truncate(tight_thresh);
        	}
        	if (world.rank()==0) print("symmetry r1");
        	r1.print_size("r1");

        	for (int l=0; l<hf->nocc(); ++l) {

        		// make the term  tmp1(1) = |UK_l(1) - 1/2 |k(1)> h(kl)
        		real_function_3d tmp1=pair.phi_l_UK_phi0[l];
        		for (int k=0; k<hf->nocc(); ++k) tmp1-= 0.5*h(k,l)*hf->nemo(k);

        		// now apply the Greens' function (including the -2.0 factor in one term)
				real_function_6d tmp=green(tmp1,-2.0*hf->nemo(l));
				tmp.set_thresh(tight_thresh);
				r1=(r1-tmp).truncate(tight_thresh);
        	}
        	r1.truncate();
        	r1.print_size("r1");
        	GVpair.set_thresh(tight_thresh);
        	GVpair=(GVpair+r1).truncate().reduce_rank();

        	FunctionDefaults<6>::set_thresh(thresh);
        	GVpair.set_thresh(thresh);
        	GVpair.truncate().reduce_rank();

        	asymmetry(GVpair,"GVpair before Q12");

			GVpair=Q12(GVpair);
			asymmetry(GVpair,"GVpair after Q12");
            pair.function=GVpair;
            pair.constant_term=copy(GVpair);
            save_function(GVpair,"GVpair");

            return pair;
        }

        /// compute the singlet and triplet energy for a given electron pair

        /// @return	the energy of 1 degenerate triplet and 1 singlet pair
        double compute_energy(ElectronPair& pair) const {

        	// a11 will be the bra space, therefore take the hermitian conjugate
            const double a11=inner(pair.function,JK1phi0_on_demand(pair.i,pair.j,true))
            				+inner(pair.function,JK2phi0_on_demand(pair.i,pair.j,true));
            if (world.rank()==0) printf("V2: <phi^0 | J-K        | psi^1>  %12.8f\n",a11);

            // this will be the bra space
            real_function_6d eri=TwoElectronFactory(world).dcut(dcut);
            real_function_6d ij_g=CompositeFactory<double,6,3>(world)
						.particle1(copy(hf->R2orbital(pair.i)))
						.particle2(copy(hf->R2orbital(pair.j)))
						.g12(eri);
            real_function_6d ji_g=CompositeFactory<double,6,3>(world)
						.particle1(copy(hf->R2orbital(pair.j)))
						.particle2(copy(hf->R2orbital(pair.i)))
						.g12(eri);

            // compute < ij | g12 | psi >
            const double ij_g_uij=inner(pair.function,ij_g);
            if (world.rank()==0) printf("<ij | g12       | psi^1>  %12.8f\n",ij_g_uij);

            // compute < ji | g12 | psi > if (i/=j)
            const double ji_g_uij= (pair.i==pair.j) ? 0 : inner(pair.function,ji_g);
            if (world.rank()==0) printf("<ji | g12       | psi^1>  %12.8f\n",ji_g_uij);

            // the singlet and triplet triplet pair energies
            if (pair.i==pair.j) {
            	pair.e_singlet=ij_g_uij + pair.ij_gQf_ij;
            	pair.e_triplet=0.0;
            } else {
                pair.e_singlet=(ij_g_uij + pair.ij_gQf_ij) + (ji_g_uij + pair.ji_gQf_ij);
                pair.e_triplet=3.0*((ij_g_uij - ji_g_uij) + (pair.ij_gQf_ij - pair.ji_gQf_ij));
            }

            // print the pair energies
            if (world.rank()==0) {
            	printf("current energy %2d %2d %12.8f %12.8f\n",pair.i,pair.j,
            			pair.e_singlet,pair.e_triplet);
            }

            // return the total energy of this pair
            return pair.e_singlet+pair.e_triplet;
        }

        /// apply the exchange operator on an orbital

        /// @param[in]	phi the orbital
        /// @param[in]	hc	hermitian conjugate -> swap bra and ket space
        /// @return 	Kphi
        real_function_3d K(const real_function_3d& phi, const bool hc=false) const {

        	real_function_3d result=real_factory_3d(world);

    		// multiply rhs with R2orbitals (the bra space)
    		vecfuncT R2rhs;
    		if (not hc) R2rhs=mul(world,phi,hf->R2orbitals());
    		else R2rhs=mul(world,phi,hf->orbitals());

    		// apply the poisson kernel and sum up
    		for (std::size_t k=0; k<hf->nemos().size(); ++k) {
				result += hf->nemo(k) * (*poisson)(R2rhs[k]);
    		}
    		return result;
        }

        /// apply the exchange operator on a pair function

        /// @param[in]	phi the pair function
        /// @param[in]	is_symmetric is the function symmetric wrt particle exchange
        /// @return 	(K1 + K2) |phi >
        real_function_6d K(const real_function_6d& phi, const bool is_symmetric=false) const {
        	real_function_6d result=real_factory_6d(world);
        	// loop over all orbitals of the reference
        	for (int i=0; i<hf->nocc(); ++i) {
        		real_function_6d tmp=apply_exchange(phi,hf->nemo(i),hf->R2orbitals()[i],1);
        		if (is_symmetric) {
        			tmp=tmp+swap_particles(tmp);
        		} else {
        			tmp=tmp+apply_exchange(phi,hf->nemo(i),hf->R2orbitals()[i],2);
        		}
        		result=(result+tmp).truncate();
        	}
        	return result;
        }

        /// apply the Coulomb operator a on orbital

        /// @param[in]	phi the orbital
        /// @return 	Jphi
        real_function_3d J(const real_function_3d& phi) const {
        	return (hf->get_coulomb_potential()*phi).truncate();
        }

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
        		const real_function_3d& orbital_bra, const int particle) const {

        	MADNESS_ASSERT((particle==1) or (particle==2));
            real_convolution_3d op=CoulombOperator(world,0.0001,hf->get_calc().param.econv);
            op.particle()=particle;

//            if (world.rank()==0) printf("start multiplication before K at time %.1f\n",wall_time());

            // multiply the orbital to the pair function
//            real_function_6d x=(particle==1)
//            		? CompositeFactory<double,6,3>(world).ket(copy(f)).V_for_particle1(copy(orbital_bra))
//            		: CompositeFactory<double,6,3>(world).ket(copy(f)).V_for_particle2(copy(orbital_bra));
//            x.fill_tree().truncate();
            real_function_6d x=multiply(copy(f),copy(orbital_bra),particle).truncate();

            // apply the Poisson operator
//            if (world.rank()==0) printf("start exchange at time %.1f\n",wall_time());
            load_balance(x,false);
            x=op(x).truncate();

            // do the final multiplication with the orbital
//            if (world.rank()==0) printf("start multiplication after K at time %.1f\n",wall_time());
//            real_function_6d result= (particle==1)
//            		? CompositeFactory<double,6,3>(world).ket(copy(x)).V_for_particle1(copy(orbital_ket))
//            		: CompositeFactory<double,6,3>(world).ket(copy(x)).V_for_particle2(copy(orbital_ket));
//            result.fill_tree().truncate().reduce_rank();
            real_function_6d result=multiply(copy(x),copy(orbital_ket),particle).truncate();

//            if (world.rank()==0) printf("end multiplication after K at time %.1f\n",wall_time());
            return result;
        }

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
        		const bool hc=false) const {

            const double tol=0.0;
            MADNESS_ASSERT(hf->get_calc().param.spin_restricted);
            std::vector<real_function_3d> psif;
            if (hc) psif = mul_sparse(world, phi, hf->nemos(), tol);
            else psif = mul_sparse(world, phi, hf->R2orbitals(), tol);
            truncate(world, psif);
            psif = apply(world, op, psif);
            truncate(world, psif);
            return psif;
        }

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
        		const bool hc=false) const {
            const double tol=0.0;
            return mul_sparse(world, phi_j, make_chi(phi_i,op,hc), tol);
        }

        /// compute the terms O1 (U - [K,f]) |phi0> and O2 UK |phi0>
        void OUKphi0(ElectronPair& pair) const {

			std::vector<real_function_3d> phi_k_UK_phi0;
			std::vector<real_function_3d> phi_l_UK_phi0;

        	if (not intermediates.OUKphi0.empty()) {
        		const std::string root=intermediates.OUKphi0;
        		for (int k=0; k<hf->nocc(); ++k) {
        			real_function_3d tmp;
        			load_function(tmp,root+"_k_"+stringify(k));
    				phi_k_UK_phi0.push_back(tmp);
        			load_function(tmp,root+"_l_"+stringify(k));
    				phi_l_UK_phi0.push_back(tmp);
    			}
        	} else {

				pair.Uphi0.verify();
				pair.KffKphi0.verify();
        		const std::string root="OUKphi0";

        		for (int k=0; k<hf->nocc(); ++k) {
					real_function_3d tmp;
					tmp=pair.Uphi0.project_out(hf->R2orbitals()[k],0)
							- pair.KffKphi0.project_out(hf->R2orbitals()[k],0);
					save_function(tmp,root+"_k_"+stringify(k));
					phi_k_UK_phi0.push_back(tmp);

					tmp=pair.Uphi0.project_out(hf->R2orbitals()[k],1)
							- pair.KffKphi0.project_out(hf->R2orbitals()[k],1);
					save_function(tmp,root+"_l_"+stringify(k));
					phi_l_UK_phi0.push_back(tmp);
				}
        	}
        	pair.phi_k_UK_phi0=phi_k_UK_phi0;
        	pair.phi_l_UK_phi0=phi_l_UK_phi0;
        }

         /// apply the operator K on the reference and multiply with f; fK |phi^0>

        /// @param[in]  i   index of orbital i
        /// @param[in]  j   index of orbital j
        /// @return     the function f12 (K(1)+K(2))|phi^0>
        real_function_6d make_fKphi0(const int i, const int j) const {
            const real_function_3d& phi_i=hf->nemo(i);
            const real_function_3d& phi_j=hf->nemo(j);

            const real_function_3d Kphi_i=K(phi_i);
            const real_function_3d Kphi_j=K(phi_j);

			real_function_6d fKphi0a=CompositeFactory<double,6,3>(world).g12(corrfac.f())
					.particle1(copy(phi_i)).particle2(copy(Kphi_j));
			fKphi0a.fill_tree().truncate();
			real_function_6d fKphi0b=CompositeFactory<double,6,3>(world).g12(corrfac.f())
					.particle1(copy(Kphi_i)).particle2(copy(phi_j));
			fKphi0b.fill_tree().truncate();

			real_function_6d fKphi0=(fKphi0a + fKphi0b).truncate();
			return fKphi0;
        }

        /// return the function (J(1)-K(1)) |phi0> as on-demand function

        /// @param[in]	hc		compute hermitian conjugate -> swap bra and ket space
        real_function_6d JK1phi0_on_demand(const int i, const int j,
        		const bool hc=false) const {
        	real_function_3d phi_i,phi_j;
            if (not hc) {
            	phi_i=hf->nemo(i);
            	phi_j=hf->nemo(j);
            } else {
            	phi_i=hf->R2orbital(i);
            	phi_j=hf->R2orbital(j);
            }

            const real_function_3d JKphi_i=J(phi_i)-K(phi_i,hc);

            real_function_6d tmp1=CompositeFactory<double,6,3>(world)
                                 .particle1(copy(JKphi_i))
                                 .particle2(copy(phi_j));
            return tmp1;
        }

        /// return the function (J(2)-K(2)) |phi0> as on-demand function
        real_function_6d JK2phi0_on_demand(const int i, const int j,
        		const bool hc=false) const {
        	real_function_3d phi_i,phi_j;
            if (not hc) {
            	phi_i=hf->nemo(i);
            	phi_j=hf->nemo(j);
            } else {
            	phi_i=hf->R2orbital(i);
            	phi_j=hf->R2orbital(j);
            }

            const real_function_3d JKphi_j=J(phi_j)-K(phi_j,hc);

            real_function_6d tmp=CompositeFactory<double,6,3>(world)
                                 .particle1(copy(phi_i))
                                 .particle2(copy(JKphi_j));
            return tmp;
        }

        /// return the function |phi0> as on-demand function
        real_function_6d phi0_on_demand(const int i, const int j) const {
            const real_function_3d& phi_i=hf->orbital(i);
            const real_function_3d& phi_j=hf->orbital(j);

            real_function_6d tmp1=CompositeFactory<double,6,3>(world)
                                 .particle1(copy(phi_i))
                                 .particle2(copy(phi_j));
            return tmp1;
        }

        /// return the function |F1F2> as on-demand function
        real_function_6d nemo0_on_demand(const int i, const int j) const {
            return CompositeFactory<double,6,3>(world)
                                 .particle1(copy(hf->nemo(i)))
                                 .particle2(copy(hf->nemo(j)));
        }


        /// multiply the given function with the 0th order Hamiltonian, exluding the 0th order energy

        /// @param[in]  f   the function we apply H^0 on
        /// @return     the function g=H^0 f, which is NOT orthogonalized against f
        real_function_6d multiply_with_0th_order_Hamiltonian(const real_function_6d& f,
        		const int i, const int j) const {

            real_function_6d vphi;

            if (0) {
//            if (hf->nemo_calc.nuclear_correlation->type()==NuclearCorrelationFactor::None) {
            	real_function_3d v_total=hf->get_nuclear_potential()+hf->get_coulomb_potential();

            	const double eps=zeroth_order_energy(i,j);
                real_convolution_6d op_mod = BSHOperator<6>(world, sqrt(-2*eps),lo,bsh_eps);
                op_mod.modified()=true;


            	vphi=CompositeFactory<double,6,3>(world)
                                 .ket(copy(f))
                                 .V_for_particle1(copy(v_total))
                                 .V_for_particle2(copy(v_total));

				// make the tree
				vphi.fill_tree(op_mod).truncate();
				vphi.print_size("(V_nuc + J1 + J2) |ket>:  made V tree");
            } else {

            	// the purely local part: Coulomb and U2
            	real_function_3d v_local=hf->get_coulomb_potential()
            			+hf->nemo_calc.nuclear_correlation->U2();

            	// screen the construction of Vphi: do only what is needed to
            	// get an accurate result of the BSH operator
            	const double eps=zeroth_order_energy(i,j);
                real_convolution_6d op_mod = BSHOperator<6>(world, sqrt(-2*eps),lo,bsh_eps);
                op_mod.modified()=true;
        		vphi=CompositeFactory<double,6,3>(world).ket(copy(f))
        				.V_for_particle1(copy(v_local)).V_for_particle2(copy(v_local));
        		vphi.fill_tree(op_mod);
    			asymmetry(vphi,"Vphi");
    			double n=vphi.norm2();
    			if (world.rank()==0) print("norm of Vphi ",n);

        		// the part with the derivative operators: U1
                for (int axis=0; axis<6; ++axis) {
                    real_derivative_6d D = free_space_derivative<double,6>(world, axis);
                    const real_function_6d Drhs=D(f).truncate();

                    // note integer arithmetic
                    if (world.rank()==0) print("axis, axis^%3, axis/3+1",axis, axis%3, axis/3+1);
                    const real_function_3d U1_axis=hf->nemo_calc.nuclear_correlation->U1(axis%3);
//                    real_function_6d x=multiply(copy(Drhs),copy(U1_axis),axis/3+1).truncate();

                    real_function_6d x;
                    if (axis/3+1==1) {
						x=CompositeFactory<double,6,3>(world)
								.ket(Drhs).V_for_particle1(copy(U1_axis));
                    } else if (axis/3+1==2) {
                    	x=CompositeFactory<double,6,3>(world)
								.ket(Drhs).V_for_particle2(copy(U1_axis));
                    }
					x.fill_tree(op_mod);
					vphi+=x;
                    vphi.truncate().reduce_rank();

                }
				vphi.print_size("(U_nuc + J) |ket>:  made V tree");
    			asymmetry(vphi,"U+J");
//    			save_function(vphi,"UJphi");
//    			plot_plane(world,vphi,"UJphi");

            }

            // and the exchange
            vphi=(vphi-K(f,i==j)).truncate().reduce_rank();
			asymmetry(vphi,"U+J-K");
			vphi.print_size("(U_nuc + J - K) |ket>:  made V tree");
//			plot_plane(world,vphi,"UJKphi");
//			save_function(vphi,"UJKphi");

            return vphi;
        }

    };
};

#endif /* MP2_H_ */

