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
#include <moldft/moldft.h>
#include <examples/nonlinsol.h>

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
        if(f.world().rank() == 0) printf("redistributed at time   %.1fs\n", wall_time());
    }

    class R1AFactor {

    public:
    	R1AFactor(const Molecule& mol, const bool inv) : inverse(inv) {
    		for (int i=0; i<mol.natom(); ++i) atoms.push_back(mol.get_atom(i));
    		print("inverse ", inverse);
    	}
    	std::vector<Atom> atoms;
    	bool inverse;

    	double operator()(const coord_3d& r) const {

    		double result=1.0;
//    		double exponent=0.0;
    		const double x=r[0], y=r[1], z=r[2];
    		for (std::vector<Atom>::const_iterator atom=atoms.begin(); atom!=atoms.end(); ++atom) {
    			const double& xA=atom->x;
    			const double& yA=atom->y;
    			const double& zA=atom->z;
            	const double rr=sqrt((x-xA)*(x-xA) + (y-yA)*(y-yA) + (z-zA)*(z-zA));
            	const double exponent=atom->q*rr;
            	result*=(1.0-exp(-exponent));
    		}
    		if (inverse) return 1.0/result;
    		return result;
        }


    };

    /// a class holding the correlation factor for R12 theory
    class CorrelationFactor {

        World& world;
        double _gamma;       ///< the correlation factor exp(-gamma r12)

    public:

        /// ctor, use negative gamma for linear correlation factor r12
        CorrelationFactor(World& world, const double& gamma) : world(world), _gamma(gamma) {
            if (world.rank()==0) {
                if (gamma>0.0) print("constructed correlation factor with gamma=",gamma);
                else if (gamma==0.0) print("constructed linear correlation factor");
            }
        }

        /// return the exponent of this correlation factor
        double gamma() const {return _gamma;}

        /// return the value of the correlation factor
        double operator()(const coord_6d& r) const {
            const double rr=r12(r);
            if (_gamma>0.0) return (1.0-exp(-_gamma*rr))/(2.0*_gamma);
            return 0.5*rr;
        }

        /// apply Kutzelnigg's regularized potential to an orbital product
        real_function_6d apply_U(const real_function_3d& phi_i, const real_function_3d& phi_j,
        		const double eps) const {
            real_function_6d result=real_factory_6d(world);

            real_convolution_6d op_mod = BSHOperator<6>(world, sqrt(-2*eps), lo,bsh_eps);
            op_mod.modified()=true;

            for (int axis=0; axis<3; ++axis) {
                if (world.rank()==0) print("working on axis",axis);
                real_derivative_3d D = free_space_derivative<double,3>(world, axis);
                const real_function_3d Di=(D(phi_i)).truncate();
                const real_function_3d Dj=(D(phi_j)).truncate();

                const real_function_6d u=real_factory_6d(world).functor2(U(_gamma,axis)).is_on_demand();

                real_function_6d tmp1=CompositeFactory<double,6,3>(world)
                            .g12(u).particle1(copy(Di)).particle2(copy(phi_j));
                tmp1.fill_tree(op_mod).truncate();
                real_function_6d tmp2=CompositeFactory<double,6,3>(world)
                            .g12(u).particle1(copy(phi_i)).particle2(copy(Dj));
                tmp2.fill_tree(op_mod).truncate();
                if (world.rank()==0) print("done with fill_tree");

                result=result+(tmp1-tmp2).truncate();
                tmp1.clear();
                tmp2.clear();
                world.gop.fence();
                result.truncate().reduce_rank();

                if (world.rank()==0) printf("done with multiplication with U at ime %.1f\n",wall_time());
                result.print_size("result");
            }

            load_balance(result,true);

            // include the purely local potential that (partially) cancels 1/r12
            if (_gamma>0.0) {
                real_function_6d fg3=real_factory_6d(world).functor2(fg_(_gamma)).is_on_demand();
                real_function_6d mul=CompositeFactory<double,6,3>(world)
                                    .g12(fg3).particle1(copy(phi_i)).particle2(copy(phi_j));
                mul.fill_tree(op_mod).truncate();
                mul.print_size("mul");

                result=(result+mul).truncate().reduce_rank();
            }
            result.print_size("U * |ij>");
            return result;
        }

        /// return the correlation factor as on-demand function
        real_function_6d f() const {
            real_function_6d tmp=real_factory_6d(world).functor2(*this).is_on_demand();
            return tmp;
        }

        /// return f^2 as on-demand function
        real_function_6d f2() const {
            real_function_6d tmp=real_factory_6d(world).functor2(f2_(_gamma)).is_on_demand();
            return tmp;
        }

        /// return fg+sth as on-demand function
        real_function_6d fg() const {
            real_function_6d tmp=real_factory_6d(world).functor2(fg_(_gamma)).is_on_demand();
            return tmp;
        }

        /// return f/r as on-demand function
        real_function_6d f_over_r() const {
            real_function_6d tmp=real_factory_6d(world).functor2(f_over_r_(_gamma)).is_on_demand();
            return tmp;
        }

        /// return (\nabla f)^2 as on-demand functions
        real_function_6d nablaf2() const {
            real_function_6d tmp=real_factory_6d(world).functor2(nablaf2_(_gamma)).is_on_demand();
            return tmp;
        }

    private:
        /// functor for the local potential (1-f12)/r12 + sth (doubly connected term of the commutator)
        struct fg_ {
            double gamma;
            fg_(double gamma) : gamma(gamma) {MADNESS_ASSERT(gamma>0.0);}
            double operator()(const coord_6d& r) const {
                const double rr=r12(r);
                const double e=exp(-gamma*rr);
                return (1.0-e)*u(rr,dcut) + 0.5*gamma*e;
            }
        };

        /// functor for the local potential (1-f12)/r12
        struct f_over_r_ {
            double gamma;
            f_over_r_(double gamma) : gamma(gamma) {MADNESS_ASSERT(gamma>0.0);}
            double operator()(const coord_6d& r) const {
                const double rr=r12(r);
                const double e=exp(-gamma*rr);
                return (1.0-e)*u(rr,dcut)/(2.0*gamma);
            }
        };

        /// functor for the local part of the regularized potential: f12/r12*(r1-r2)(D1-D2)
        struct U {
            double gamma;
            int axis;
            U(double gamma, int axis) : gamma(gamma), axis(axis) {
                MADNESS_ASSERT(axis>=0 and axis<3);
            }
            double operator()(const coord_6d& r) const {
                const double rr=r12(r);
                const double g12=u(rr,dcut);
                double a=0.5;
                if (gamma>0.0) a=0.5*exp(-gamma*rr);
                return -a*x12(r,axis) * g12;
            }
        };

        /// functor for the local potential (1-f12)^2
        struct f2_ {
            double gamma;
            f2_(double gamma) : gamma(gamma) {MADNESS_ASSERT(gamma>0.0);}
            double operator()(const coord_6d& r) const {
                const double rr=r12(r);
                const double e=exp(-gamma*rr);
                const double f=(1.0-e)/(2.0*gamma);
                return f*f;
            }
        };

        /// functor for the local potential (\nabla f)^2
        struct nablaf2_ {
            double gamma;
            nablaf2_(double gamma) : gamma(gamma) {
            	MADNESS_ASSERT(gamma>0.0);
            	MADNESS_ASSERT(gamma=1.0);	// I don't think this is right
            }
            double operator()(const coord_6d& r) const {
                const double rr=r12(r);
                const double f=exp(-2.0*gamma*rr)/(4.0*gamma*gamma);
                return f;
            }
        };

        /// Smoothed 1/r potential (c is the smoothing distance)
        static double u(double r, double c) {
            r = r/c;
            double r2 = r*r, pot;
            if (r > 6.5){
                pot = 1.0/r;
            } else if (r > 1e-2) {
                pot = erf(r)/r + exp(-r2)*0.56418958354775630;
            } else{
                pot = 1.6925687506432689-r2*(0.94031597257959381-r2*(0.39493270848342941-0.12089776790309064*r2));
            }
            return pot/c;
        }

        static double r12(const coord_6d& r) {
            const double x12=r[0]-r[3];
            const double y12=r[1]-r[4];
            const double z12=r[2]-r[5];
            const double r12=sqrt(x12*x12 + y12*y12 + z12*z12);
            return r12;
        }
        static double x12(const coord_6d& r, const int axis) {
            return r[axis]-r[axis+3];
        }


    };



    /// simple projector class for 1- and 2-particle projectors
    template<typename T, std::size_t NDIM>
    class Projector {

        int particle_;
        std::vector<Function<T,NDIM> > p_;

    public:

        Projector() : p_(std::vector<Function<T,NDIM> >()) {}

        /// simple constructor with only one orbital to project out
        Projector(const Function<T,NDIM>& p, const int particle=0)
        	: particle_(particle), p_(std::vector<Function<T,NDIM> >(1,p)) {
            MADNESS_ASSERT(particle_==0 or particle_==1);
            MADNESS_ASSERT(p_.size()>0);
        }

        /// constructor with a set of orbitals to project out
        Projector(const std::vector<Function<T,NDIM> >& p, const int particle=0) : particle_(particle), p_(p) {
            MADNESS_ASSERT(particle_==0 or particle_==1);
            MADNESS_ASSERT(p_.size()>0);
        }

        int& particle() {return particle_;}
        const int& particle() const {return particle_;}

        /// get a const reference to the orbitals
        const std::vector<Function<T,NDIM> >& p() const {return p_;}

        /// project f on p: |result> =  | p><p | f>
        template<std::size_t FDIM>
        typename enable_if_c<NDIM==FDIM, Function<T,FDIM> >::type
        operator()(const Function<T,FDIM>& f) const {

            const double ovlp=inner(f,p_[0]);
            Function<T,NDIM> sum=ovlp*p_[0];

            for (unsigned int i=1; i<p_.size(); ++i) {
                const double ovlp2=inner(f,p_[i]);
                sum=(sum+ovlp2*p_[i]).truncate().reduce_rank();
            }
            return sum;
        }

        /// project p out of f: |result(1,2)> = sum_p | p(1)><p(1) | f(1,2)>
        template<std::size_t FDIM>
        typename enable_if_c<2*NDIM==FDIM, Function<T,FDIM> >::type
        operator()(const Function<T,FDIM>& f) const {
            real_function_6d sum=real_factory_6d(p_.begin()->world());
            for (unsigned int i=0; i<p_.size(); ++i) {
                const real_function_3d pf2=f.project_out(p_[i],particle_);
                real_function_6d tmp;
                MADNESS_EXCEPTION("Projector class: the hartree product is inaccurate -- don't use it",1);
                if (particle_==0) tmp=hartree_product(p_[i],pf2);
                else tmp=hartree_product(pf2,p_[i]);
                sum=(sum+tmp);
            }
            sum.truncate();
            return sum;
        }
    };

    class HartreeFock : public OptimizationTargetInterface {
        World& world;
        Calculation calc;
        mutable double coords_sum;     // sum of square of coords at last solved geometry
        mutable double E; //< Current energy

        // save the Coulomb potential
        mutable real_function_3d coulomb;

    public:

        HartreeFock(World& world, Calculation& calc)
            : world(world)
            , calc(calc)
            , coords_sum(-1.0)
        {
        }

        bool provides_gradient() const {return true;}

        double value() {
            return value(calc.molecule.get_all_coords());
        }

        double value(const Tensor<double>& x) {

        	// fast return if the reference is already solved at this geometry
            double xsq = x.sumsq();
            if (xsq == coords_sum) return calc.current_energy;

            calc.molecule.set_all_coords(x.reshape(calc.molecule.natom(),3));
            coords_sum = xsq;

            // Make the nuclear potential, initial orbitals, etc.
            calc.make_nuclear_potential(world);
            calc.vnuc.print_size("vnuc");
            calc.project_ao_basis(world);

            // read converged wave function from disk if there is one
            if (calc.param.no_compute) {
                calc.load_mos(world);
                return calc.current_energy;
            }

            if (calc.param.restart) {
                calc.load_mos(world);
            } else {
                calc.initial_guess(world);
                calc.param.restart = true;
            }

            // If the basis for the inital guess was not sto-3g
            // switch to sto-3g since this is needed for analysis
            // of the MOs and orbital localization
            if (calc.param.aobasis != "sto-3g") {
            	calc.param.aobasis = "sto-3g";
            	calc.project_ao_basis(world);
            }

            calc.solve(world);
            calc.save_mos(world);

            // successively tighten threshold
            if (calc.param.econv<1.1e-6) {
            	calc.set_protocol<3>(world,1e-6);
            	calc.make_nuclear_potential(world);
            	calc.project_ao_basis(world);
            	calc.project(world);
            	calc.solve(world);
            	calc.save_mos(world);
            }

            calc.save_mos(world);

            return calc.current_energy;
        }

        Tensor<double> gradient(const Tensor<double>& x) {

            value(x); // Ensures DFT equations are solved at this geometry
            return calc.derivatives(world);
        }

        double coord_chksum() const {return coords_sum;}

        const Calculation& get_calc() const {return calc;}
        Calculation& get_calc() {return calc;}

        /// return orbital i
        real_function_3d orbital(const int i) const {
            MADNESS_ASSERT(calc.param.spin_restricted);
            return calc.amo[i];
        }

        /// return orbital energy i
        double orbital_energy(const int i) const {
            MADNESS_ASSERT(calc.param.spin_restricted);
            return calc.aeps[i];
        }

        /// return the Coulomb potential
        real_function_3d get_coulomb_potential() const {
            MADNESS_ASSERT(calc.param.spin_restricted);
            if (coulomb.is_initialized()) return copy(coulomb);
            functionT rho = calc.make_density(world, calc.aocc, calc.amo).scale(2.0);
            coulomb=calc.make_coulomb_potential(rho);
            return copy(coulomb);
        }

        /// return the nuclear potential
        real_function_3d get_nuclear_potential() const {
            return calc.vnuc;
        }

        /// apply the exchange operator
        real_function_3d apply_exchange(const real_function_3d& phi) const {

            MADNESS_ASSERT(calc.param.spin_restricted);
            vecfuncT phi_vec=std::vector<real_function_3d> (1,phi);
            vecfuncT result=calc.apply_hf_exchange(world,calc.aocc,calc.amo,phi_vec);
            return result[0];
        }

        /// return the number of occupied orbitals
        int nocc() const {
            MADNESS_ASSERT(calc.param.spin_restricted);
            return calc.param.nalpha;
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

    	static double uninitialized() {return 1.e10;}

        int i, j;                       ///< orbitals i and j
        real_function_6d function;      ///< pair function for a specific pair w/o correlation factor part
        real_function_6d r12phi;       ///< orbital product multiplied with the correlation factor

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
        	bool exist=function.is_initialized();
			ar & ij_gQf_ij & ji_gQf_ij & e_singlet & e_triplet & converged & iteration & exist;
			if (exist) ar & function;
        }

        void load_pair(World& world) {
        	std::string name="pair_"+stringify(i)+stringify(j);
        	if (world.rank()==0) printf("loading matrix elements %s",name.c_str());
            archive::ParallelInputArchive ar(world, name.c_str(), 1);
        	ar & *this;
        	if (world.rank()==0) printf(" %s\n",(converged)?" converged":" not converged");
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

        	/// ctor reading out the input file
        	Parameters(const std::string& input) : thresh_(-1.0), i(-1), j(-1), freeze(0), restart(false) {

        		// get the parameters from the input file
                std::ifstream f(input.c_str());
                position_stream(f, "mp2");
                std::string s;

                while (f >> s) {
                    if (s == "end") break;
                    else if (s == "econv") f >> thresh_;
                    else if (s == "pair") f >> i >> j;
                    else if (s == "freeze") f >> freeze;
                    else if (s == "restart") restart=true;
                    else continue;
                }
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
        Parameters param;						///< calculation parameters for MP2
        std::shared_ptr<HartreeFock> hf;        ///< our reference
        CorrelationFactor corrfac;              ///< correlation factor: Slater

        std::map<std::pair<int,int>,ElectronPair> pairs;       ///< pair functions and energies
        double correlation_energy;				///< the correlation energy
        double coords_sum;						///< check sum for the geometry

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
        MP2(World& world, const CorrelationFactor& corrfac, const std::string& input)
            : world(world)
        	, param(input.c_str())
//            , hf(world,Calculation(world,input.c_str()))
            , corrfac(corrfac)
        	, correlation_energy(0.0)
    		, coords_sum(-1.0)
            , intermediates() {
//            , poisson(CoulombOperator(world,0.0001,param.thresh_*0.01)) {

        	{
        		Calculation calc(world,input.c_str());
                // get parameters form input file for hf

        		if (world.rank()==0) print("accuracy from dft will be overriden by mp2 to 0.01*thresh");
                calc.set_protocol<6>(world,thresh());
                calc.param.econv=thresh()*0.01;
                calc.set_protocol<3>(world,calc.param.econv);

        		// override computed parameters if they are provided explicitly
                double eprec=calc.param.econv*0.1;

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

                calc.molecule.set_eprec(eprec);
                calc.molecule.print();

                hf=std::shared_ptr<HartreeFock>(new HartreeFock(world,calc));
                poisson=std::shared_ptr<real_convolution_3d>
                	(CoulombOperatorPtr(world,0.0001,calc.param.econv));
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
    				if (param.restart) {
    					try {
        					pair(i,j).load_pair(world);
    				    } catch (std::exception& e) {
    				    	if (world.rank()==0) print("could not find pair ",i,j," on disk");
    				    }
    				}
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


        	correlation_energy=0.0;
        	Projector<double,3> O1=O(0);
        	Projector<double,3> O2=O(1);

        	// compute only one single pair
        	if ((param.i>-1) and (param.j>-1)) {
                pair(param.i,param.j)=solve_residual_equations(param.i,param.j);
                correlation_energy+=pair(param.i,param.j).e_singlet+pair(param.i,param.j).e_triplet;

            // solve the residual equations for all pairs ij
        	} else {
        		for (int i=param.freeze; i<hf->nocc(); ++i) {
        			for (int j=i; j<hf->nocc(); ++j) {
        				if (pair(i,j).converged) {
							correlation_energy+=pair(i,j).e_singlet+pair(i,j).e_triplet;
        				} else {
							pair(i,j)=solve_residual_equations(i,j);
							correlation_energy+=pair(i,j).e_singlet+pair(i,j).e_triplet;
        				}
					}
				}
        	}
            return correlation_energy;
        }

        /// print the calculation parameters
        void print_info(World& world) const {
            if (world.rank()==0) {
                madness::print("\n === MP2 info === \n");
                madness::print("         MP2 restart ", param.restart);
                madness::print("        threshold 3D ", FunctionDefaults<3>::get_thresh());
                madness::print("        threshold 6D ", FunctionDefaults<6>::get_thresh());
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
            increment(result,green);

			if (world.rank()==0) print("computing iteratively");
			real_function_6d constant_term;
			load_function(constant_term,"GVpair");

			NonlinearSolverND<6> solver;
			// increment iteration counter upon entry
			for (++result.iteration; result.iteration<20; ++result.iteration) {

				// apply the convolution
				real_function_6d vphi=multiply_with_0th_order_Hamiltonian(result.function,i,j);
				vphi.scale(-2.0).truncate();
				load_balance(vphi,false);

				real_function_6d tmp=green(vphi).truncate();	// green is destructive
				tmp.print_size("result of applying 0th order Hamiltonian on 1st order wave function");
				tmp=Q12(constant_term+tmp);

				// check convergence
				real_function_6d residual=tmp-result.function;
				const double rnorm=residual.norm2();
				const double old_fnorm=result.function.norm2();
				const double fnorm=tmp.norm2();
				if (world.rank()==0) printf("norm2 of psi, residual %12.8f %12.8f\n",fnorm,rnorm);

				result.function=solver.update(result.function,residual);
//				result.function=tmp;
				double old_energy=energy;
				energy=compute_energy(result);

				result.converged=((std::abs(old_energy-energy)<result.function.thresh()*0.01)
						and (std::abs(old_fnorm-fnorm)<result.function.thresh()*0.1));
				result.store_pair(world);

				if (world.rank()==0) printf("finished iteration %2d at time %8.1fs with energy %12.8f\n\n",
						result.iteration, wall_time(),energy);

				if (result.converged) break;

			}

            // save the converged first order pair function separately for easier access
        	std::string name="pair_"+stringify(i)+stringify(j)+"_psi1_converged";
            save_function(result.function,name);

            // print the final pair energies
            if (world.rank()==0) {
            	printf("final correlation energy %2d %2d %12.8f %12.8f\n",result.i,result.j,
            			result.e_singlet,result.e_triplet);
            }

            return result;

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

        void test3() {
        	hf->value();
        	Projector<double,3> O1(hf->get_calc().amo);
        	real_function_3d vnuc=hf->get_nuclear_potential();
        	save_function(hf->get_nuclear_potential(),"vnuc");
        	real_function_3d phi=hf->orbital(param.i);
        	real_function_3d f1a=real_factory_3d(world).functor2(R1AFactor(hf->get_calc().molecule,false));
        	real_function_3d f1a_inv=real_factory_3d(world).functor2(R1AFactor(hf->get_calc().molecule,true));
        	real_function_3d f1a_phi=f1a*phi;
        	save_function(f1a,"f1a");
        	save_function(f1a_inv,"f1a_inv");
        	save_function(f1a_phi,"f1a_phi");

        	real_function_3d OR=(f1a-O1(f1a));
        	save_function(OR,"OR");

        	real_function_3d f1a_vnuc=vnuc*f1a;
        	save_function(f1a_vnuc,"f1a_vnuc");

        	real_function_3d OR_vnuc=vnuc*OR;
        	save_function(OR_vnuc,"OR_vnuc");

        	real_function_3d f1a_inv_phi=f1a_inv*phi;
        	save_function(f1a_inv,"f1a_inv");
        	save_function(f1a_inv_phi,"f1a_inv_phi");
        }

        void test2() {

        	real_function_6d psi, r12phi;
        	load_function(psi,"psi1_it15");
        	load_function(r12phi,"r12phi");
        	r12phi=Q12(r12phi);

        	const real_function_6d psi_singlet=(psi+swap_particles(psi)).truncate();
        	const real_function_6d psi_triplet=(psi-swap_particles(psi)).truncate();

        	save_function(psi_singlet,"psi1_singlet");
        	save_function(psi_triplet,"psi1_triplet");

        	const real_function_6d r12phi_singlet=(r12phi+swap_particles(r12phi)).truncate();
        	const real_function_6d r12phi_triplet=(r12phi-swap_particles(r12phi)).truncate();

        	const real_function_6d singlet=(r12phi_singlet+psi_singlet).truncate();
        	const real_function_6d triplet=(r12phi_triplet+psi_triplet).truncate();

        	save_function(singlet,"phi1_singlet");
        	save_function(triplet,"phi1_triplet");
        }

        void test4() {

        	hf->value();
        	ElectronPair pair=make_pair(1,1);

        	const double thresh=FunctionDefaults<6>::get_thresh();
        	const double tight_thresh=thresh*0.1;
        	FunctionDefaults<6>::set_thresh(tight_thresh);
        	real_function_6d phi0=hartree_product(hf->orbital(pair.i),hf->orbital(pair.j));
        	FunctionDefaults<6>::set_thresh(thresh);


            real_function_6d eri=ERIFactory<double,6>(world).dcut(dcut);
            real_function_6d ij_g=CompositeFactory<double,6,3>(world)
						.g12(eri).ket(phi0);

            // compute < ij | g12 | psi >
            const double ij_g_uij=inner(pair.function,ij_g);
            if (world.rank()==0) printf("<ij | g12       | psi^1>  %12.8f\n",ij_g_uij);


            const double ji_g_uij= (pair.i==pair.j) ? 0 : -100.0;
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

            // compute < ji | g12 | psi > if (i/=j)

        	real_function_6d phi1=Q12(pair.r12phi+pair.function);
        	save_function(phi1,"pair11_phi1_full_singlet");
        }

		void test() {
		}
		
        /// compute the matrix element <ij | g12 Q12 f12 | phi^0>

        /// as for the formulas cf the article mra_molecule
        /// @return 	the energy <ij | g Q f | kl>
        double compute_gQf(const int i, const int j, ElectronPair& pair) const {

        	// for clarity of notation
        	const int k=pair.i;
        	const int l=pair.j;
            const real_function_3d& phi_i=hf->orbital(i);
            const real_function_3d& phi_j=hf->orbital(j);
            const real_function_3d& phi_k=hf->orbital(k);
            const real_function_3d& phi_l=hf->orbital(l);

            // compute <ij| fg |kl>: do it in 3D as (ik| fg |jl)
            // the operator fg can be rewritten as 1/r12 - f/r12
            // i.e. as poisson kernel and a bsh kernel. Note the
            // the bsh kernel includes a factor of 1/(4 pi)
            const real_function_3d ik=phi_i*phi_k;
            const real_function_3d jl=phi_j*phi_l;

            const double fourpi=4.0*constants::pi;
            real_convolution_3d fg = BSHOperator<3>(world, corrfac.gamma(), lo, bsh_eps/fourpi);
            real_convolution_3d gg = CoulombOperator(world,lo,bsh_eps);
            const real_function_3d ik_fg=(gg)(ik) - fourpi*fg(ik);
            const double a=inner(ik_fg,jl)/(2.0*corrfac.gamma());

            if (world.rank()==0) printf("<%d%d | f/r          | %d%d>  %12.8f <ij| . |ij>\n",i,j,k,l,a);

            // compute <ij| g O f | ij>
            std::vector<real_function_3d> xi_ij=make_xi(phi_i,phi_j);       // xi_{ki},j
            std::vector<real_function_3d> xi_ji=make_xi(phi_j,phi_i);       // xi_{kj},i
            if (world.rank()==0) printf("<k xi(ik,j) | f        | kl>\n");
            double b=0.0;
            for (int kk=0; kk<hf->nocc(); ++kk) {

                const real_function_3d& phi_kk=hf->orbital(kk);
                real_function_6d tmp1=CompositeFactory<double,6,3>(world)
                                        .particle1(copy(phi_kk))
                                        .particle2(copy(xi_ij[kk]));
                const double b1=inner(pair.r12phi,tmp1);
                if (world.rank()==0) printf("<%d xi(%d%d,%d) | f        | %d%d>  %12.8f b1\n",kk,i,kk,j,k,l,b1);

                real_function_6d tmp3=CompositeFactory<double,6,3>(world)
                						.particle1(copy(xi_ji[kk]))
                						.particle2(copy(phi_kk));
                const double b3=inner(pair.r12phi,tmp3);
                if (world.rank()==0) printf("<%d xi(%d%d,%d) | f        | %d%d>  %12.8f b3\n",j,kk,i,kk,k,l,b3);


                b+=b1+b3;
            }
            if (world.rank()==0) printf("<%d%d | g (O1+O2) f      | %d%d>  %12.8f <ij| . |ij>\n",i,j,k,l,b);

            // compute <ij| g O1O2 f |ij>
            double c=0.0;
            for (int kk=0; kk<hf->nocc(); ++kk) {
            	for (int ll=0; ll<hf->nocc(); ++ll) {
                    const real_function_3d& phi_kk=hf->orbital(kk);
                    const real_function_3d& phi_ll=hf->orbital(ll);

                    const double g_ijkl=inner(xi_ij[kk],phi_ll);
                    if (world.rank()==0) printf("<%d%d | g                | %d%d>  %12.8f\n",kk,ll,i,j,g_ijkl);
                    real_function_6d tmp2=CompositeFactory<double,6,3>(world)
                                            .particle1(copy(phi_kk))
                                            .particle2(copy(phi_ll));
                    const double f_ijkl=inner(pair.r12phi,tmp2);
                    if (world.rank()==0) printf("<%d%d | f                | %d%d>  %12.8f\n",kk,ll,i,j,f_ijkl);
                    c+=g_ijkl*f_ijkl;
                }
            }

            if (world.rank()==0) printf("<%d%d | g O1 O2 f        | %d%d>  %12.8f <ij| . |ij>\n",i,j,k,l,c);

            const double e=a-b+c;
            if (world.rank()==0) printf("<%d%d | g Q12 f          | %d%d>  %12.8f <ij| . |ij>\n",i,j,k,l,e);

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

        /// return the projector on the occupied orbitals

        /// need to encapsulate because the projector can only be constructed after
        /// an HF calculation has been performed
        const Projector<double,3> O(const int particle) const {
        	hf->value();	// make sure there is a converged reference
        	return Projector<double,3>(hf->get_calc().amo,particle);
        }

        /// apply the strong orthogonality operator Q12 on a function f
        real_function_6d Q12(const real_function_6d& f) const {

        	// simple and it works for higher accuracies, but might be
        	// imprecise for lower accuracies
//        	return (f-O1(f)-O2(f)+O1(O2(f))).truncate().reduce_rank();

        	const double thresh=FunctionDefaults<6>::get_thresh();
        	const double tight_thresh=FunctionDefaults<6>::get_thresh()*0.1;

        	Projector<double,3> O1=O(0);
        	Projector<double,3> O2=O(1);
        	const std::vector<real_function_3d>& O1_mos=O1.p();
        	const std::vector<real_function_3d>& O2_mos=O2.p();

        	// note no (kl) symmetry here!
        	Tensor<double> g_kl(O1_mos.size(),O2_mos.size());
        	for (size_t k=0; k<O1_mos.size(); ++k) {
        		for (size_t l=0; l<O2_mos.size(); ++l) {
    	            real_function_6d kl=CompositeFactory<double,6,3>(world)
    	            		.particle1(copy(O1_mos[k])).particle2(copy(O2_mos[l]));
        			g_kl(k,l)=inner(f,kl);
        		}
        	}
//        	if (world.rank()==0) {print(g_kl);};

        	// project out the mainly first particle: O1 (1 - 1/2 O2)
        	real_function_6d r1=real_factory_6d(world);
        	for (size_t k=0; k<O1_mos.size(); ++k) {
        		real_function_3d h2=f.project_out(O1_mos[k],0);
            	for (size_t l=0; l<O2_mos.size(); ++l) {
            		h2-=0.5*g_kl(k,l)*O2_mos[l];
            	}
            	// the hartree product tends to be inaccurate; tighten threshold
            	FunctionDefaults<6>::set_thresh(tight_thresh);
            	r1=(r1+hartree_product(O1_mos[k],h2));
            	FunctionDefaults<6>::set_thresh(thresh);
            	r1.set_thresh(thresh);
            	r1.print_size("r1"+stringify(k));
        	}

        	// project out the mainly second particle: O2 (1 - 1/2 O1)
        	real_function_6d r2=real_factory_6d(world);
        	for (size_t l=0; l<O2_mos.size(); ++l) {
        		real_function_3d h1=f.project_out(O2_mos[l],1);
            	for (size_t k=0; k<O1_mos.size(); ++k) {
            		h1-=0.5*g_kl(k,l)*O1_mos[k];			// ordering g(k,l) is correct
            	}
            	// the hartree product tends to be inaccurate; tighten threshold
            	FunctionDefaults<6>::set_thresh(tight_thresh);
            	r2=(r2+hartree_product(h1,O2_mos[l]));
            	r2.set_thresh(thresh);
            	FunctionDefaults<6>::set_thresh(thresh);
            	r2.print_size("r2"+stringify(l));
        	}
        	FunctionDefaults<6>::set_thresh(tight_thresh);
        	real_function_6d result=(f-r1-r2).truncate().reduce_rank();
        	FunctionDefaults<6>::set_thresh(thresh);

//        	// for debugging purposes only: check orthogonality
//        	for (size_t k=0; k<hf->nocc(); ++k) {
//        		for (size_t l=0; l<hf->nocc(); ++l) {
//    	            real_function_6d kl=CompositeFactory<double,6,3>(world)
//    	            		.particle1(copy(O1_mos[k])).particle2(copy(O2_mos[l]));
//        			g_kl(k,l)=inner(result,kl);
//        		}
//        	}
//        	if (world.rank()==0) {print(g_kl);};

        	return result;
        }

        /// return the function Uphi0; load from disk if available
        real_function_6d make_Uphi0(ElectronPair& pair) const {
        	const int i=pair.i;
        	const int j=pair.j;
        	real_function_6d Uphi0;
            if (not intermediates.Uphi0.empty()) {
                load_function(Uphi0,intermediates.Uphi0);
            } else {
            	const double eps=this->zeroth_order_energy(i,j);
                Uphi0=corrfac.apply_U(hf->orbital(i),hf->orbital(j),eps);
                save_function(Uphi0,"Uphi0");
            }

            // sanity check: <ij| [T,g12] |ij> = <ij | U |ij> - <ij| g12 | ij> = 0
            real_function_6d phi0=this->phi0_on_demand(i,j);
            const double a=inner(Uphi0,phi0);
            const real_function_3d ii=hf->orbital(i)*hf->orbital(i);
            const real_function_3d jj=hf->orbital(j)*hf->orbital(j);
            const real_function_3d gii=(*poisson)(ii);
            const double aa=inner(jj,gii);
            const double error=std::fabs(a-aa);
            if (world.rank()==0) {
            	printf("< phi0 | U     | phi0 >  %12.8f\n",a);
                if (error>thresh()) print("WARNING : Kutzelnigg's potential inaccurate");
                if (error>thresh()*10.0) MADNESS_EXCEPTION("Kutzelnigg's potential plain wrong",1);
            }
            return Uphi0;
        }

        /// return the function [K,f] phi0; load from disk if available
        real_function_6d make_KffKphi0(const ElectronPair& pair) const {

        	real_function_6d KffKphi0;
			if (not intermediates.KffKphi0.empty()) {
				load_function(KffKphi0,intermediates.KffKphi0);

			} else {
				real_function_6d Kfphi0;
				if (not intermediates.Kfphi0.empty()) {
					load_function(Kfphi0,intermediates.Kfphi0);
				} else {
					Kfphi0=K(pair.r12phi,pair.i==pair.j);
				}

				{
					real_function_6d phi0=this->phi0_on_demand(pair.i,pair.j);
					double a1=inner(Kfphi0,phi0);
					if (world.rank()==0) printf("< phi0 | K f    | phi0 >  %12.8f\n",a1);
					save_function(Kfphi0,"Kfphi0");
				}
				const real_function_6d fKphi0=make_fKphi0(pair.i,pair.j);
				{
					real_function_6d phi0=this->phi0_on_demand(pair.i,pair.j);
					double a2=inner(fKphi0,phi0);
					if (world.rank()==0) printf("< phi0 | f K    | phi0 >  %12.8f\n",a2);
				}
	            KffKphi0=(Kfphi0-fKphi0).truncate().reduce_rank();
				save_function(KffKphi0,"KffKphi0");
			}

			// sanity check
            real_function_6d phi0=this->phi0_on_demand(pair.i,pair.j);
            const double a=inner(KffKphi0,phi0);
            if (world.rank()==0) {
            	printf("< phi0 | [K,f]  | phi0 >  %12.8f\n",a);
                if (std::fabs(a)>thresh()) print("WARNING : exchange commutator inaccurate");
                if (std::fabs(a)>thresh()*10.0) MADNESS_EXCEPTION("exchange commutator plain wrong",1);
            }

			return KffKphi0;
        }

        /// compute some matrix elements that don't change during the calculation
        ElectronPair make_pair(const int i, const int j) const {

            ElectronPair p=ElectronPair(i,j);
			try {
				p.load_pair(world);
		    } catch (std::exception& e) {
		    	if (world.rank()==0) print("could not find pair ",i,j," on disk");
		    }

//        	ElectronPair p=pair(i,j);

            // some functions repeatedly used
            if (intermediates.r12phi.empty()) {
                p.r12phi=CompositeFactory<double,6,3>(world)
                            .g12(corrfac.f())
                            .particle1(copy(hf->orbital(i)))
                            .particle2(copy(hf->orbital(j)));
                p.r12phi.fill_tree().truncate().reduce_rank();
                save_function(p.r12phi,"r12phi");
            } else {
            	load_function(p.r12phi,intermediates.r12phi);
            }

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

			pair.KffKphi0=make_KffKphi0(pair);
            pair.Uphi0=make_Uphi0(pair);

            // these are the terms that come from the single projectors: (O1 + O2) (U+[K,f])|phi^0>
			OUKphi0(pair);

			// make the terms with high ranks and smallish trees
			load_balance(pair.Uphi0,true);
			real_function_6d Vpair1=(pair.Uphi0-pair.KffKphi0).truncate().reduce_rank();
			pair.Uphi0.clear();
			pair.KffKphi0.clear();
			load_balance(Vpair1,false);
			real_function_6d GVpair;
            green.destructive()=true;			// green will destroy Vpair1
			GVpair=green(-2.0*Vpair1).truncate().reduce_rank();
			Vpair1.clear(true);
			save_function(GVpair,"GVpair1");
//			load_function(GVpair,"GVpair1");

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
			for (int k=0; k<hf->nocc(); ++k) {

	            // for the matrix element <ij|kl> = (ik|jl)
	            const real_function_3d ik=hf->orbital(i)*hf->orbital(k);
	            const real_function_3d gik=(*poisson)(ik).truncate();

	            // the function Of(2) = <k(1) | f12 | phi0 >
	            const real_function_3d Of=pair.r12phi.project_out(hf->orbital(k),0);

				for (int l=0; l<hf->nocc(); ++l) {

					const real_function_3d jl=hf->orbital(j)*hf->orbital(l);
		            const double g_ijkl=inner(jl,gik);

		            // the matrix element <kl | f12 | ij>
		            const double f_ijkl=inner(hf->orbital(l),Of);
                    const double e_kl=zeroth_order_energy(k,l);
                    const double e_ij=zeroth_order_energy(i,j);

                    h(k,l)=f_ijkl*(e_kl-e_ij)+g_ijkl;
				}
			}

			// the result function
        	real_function_6d r1=real_factory_6d(world);
        	for (int k=0; k<hf->nocc(); ++k) {

        		// make the term  tmp2(2) = |UK_k(2)> - 1/2 |l(2)> h(kl)
        		real_function_3d tmp2=pair.phi_k_UK_phi0[k];
        		for (int l=0; l<hf->nocc(); ++l) tmp2-= 0.5*h(k,l)*hf->orbital(l);

        		// now apply the Greens' function (including the -2.0 factor in one term)
				r1=r1-green(-2.0*hf->orbital(k),tmp2);
        	}

        	for (int l=0; l<hf->nocc(); ++l) {

        		// make the term  tmp1(1) = |UK_l(1) - 1/2 |k(1)> h(kl)
        		real_function_3d tmp1=pair.phi_l_UK_phi0[l];
        		for (int k=0; k<hf->nocc(); ++k) tmp1-= 0.5*h(k,l)*hf->orbital(k);

        		// now apply the Greens' function (including the -2.0 factor in one term)
				r1=r1-green(tmp1,-2.0*hf->orbital(l));
        	}
        	GVpair=(GVpair+r1).truncate().reduce_rank();

			GVpair=Q12(GVpair);
            pair.function=GVpair;
            save_function(GVpair,"GVpair");
            pair.store_pair(world);

            return pair;
        }

        /// given a pair function, compute the pair energy using the Hylleraas functional
        double compute_second_order_correction_with_Hylleraas(const ElectronPair& pair) const {

            const double V=compute_V(pair);
            const double B=compute_B_with_U(pair);
//            const double B=1.0;

            const double e=2.0*V+B;
            const double e_opt=-V*V/B;

            if (world.rank()==0) {
                printf("V, B, and (%2d %2d) pair energy : %12.8f %12.8f %12.8f %12.8f\n\n",
                        pair.i,pair.j,V,B,e,e_opt);
            }
            return e_opt;
        }

        /// compute the V matrix for a given electron pair
        double compute_V(const ElectronPair& pair) const {

        	MADNESS_EXCEPTION("no compute_V",0);
            double V=0.0;
            const double a11=inner(pair.function,JK1phi0_on_demand(pair.i,pair.j))
            				+inner(pair.function,JK2phi0_on_demand(pair.i,pair.j));
            if (world.rank()==0) printf("V2: <phi^0 | J-K        | psi^1>  %12.8f\n",a11);
//            V-=a11;

            // two-electron interaction potential
            real_function_6d eri=ERIFactory<double,6>(world).dcut(dcut);
            real_function_6d tmp1=CompositeFactory<double,6,3>(world)
						.particle1(copy(hf->orbital(pair.i)))
						.particle2(copy(hf->orbital(pair.j)))
						.g12(eri);
            const double a12=inner(pair.function,tmp1);

//            const double a12=inner(vpsi1,this->phi0_on_demand(pair.i,pair.j));
            if (world.rank()==0) printf("V2: <phi^0 | g12        | psi^1>  %12.8f\n",a12);
            V+=a12;
            if (world.rank()==0) printf("V2: <phi^0 | H^1        | psi^1>  %12.8f\n",V);

//            if (world.rank()==0) printf("V2: <phi^0 | H^1 Q12 f  | phi^0>  %12.8f\n",pair.phi0_gQf_phi0);
//            V+=pair.phi0_gQf_phi0;
            if (world.rank()==0) printf("V2: <phi^0 | V          | phi^1>  %12.8f\n",V);

            {
            	real_function_6d phi0=this->phi0_on_demand(pair.i,pair.j);
            	phi0.fill_tree();
            	real_function_6d vpsi1=CompositeFactory<double,6,3>(world)
                    .ket(copy(phi0)).g12(eri);

            	const double a12=inner(pair.function,vpsi1);
            	if (world.rank()==0) printf("V2: <phi^0 | g12 num2   | psi^1>  %12.8f\n",a12);
            }

            return V;

        }

        /// compute the singlet and triplet energy for a given electron pair

        /// @return	the energy of 1 degenerate triplet and 1 singlet pair
        double compute_energy(ElectronPair& pair) const {

            const double a11=inner(pair.function,JK1phi0_on_demand(pair.i,pair.j))
            				+inner(pair.function,JK2phi0_on_demand(pair.i,pair.j));
            if (world.rank()==0) printf("V2: <phi^0 | J-K        | psi^1>  %12.8f\n",a11);

            real_function_6d eri=ERIFactory<double,6>(world).dcut(dcut);
            real_function_6d ij_g=CompositeFactory<double,6,3>(world)
						.particle1(copy(hf->orbital(pair.i)))
						.particle2(copy(hf->orbital(pair.j)))
						.g12(eri);
            real_function_6d ji_g=CompositeFactory<double,6,3>(world)
						.particle1(copy(hf->orbital(pair.j)))
						.particle2(copy(hf->orbital(pair.i)))
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

        /// compute the B matrix for a given electron pair
        double compute_B(const int i, const int j, const ElectronPair& pair) const {

            MADNESS_EXCEPTION("no compute_B for the time being",1);
            // first the terms  <phi^0| f12 Q (H-E^0) Q f12 |phi^0>
            double B=0.0;

#if 0
            const double e0=zeroth_order_energy(i,j);
            const real_function_6d phi0=this->phi0_on_demand(i,j);
            const real_function_6d& psi1=pair.function;


            const double nf=pair.phi0_nf_phi0;
            const double f2=pair.phi0_f2_phi0;
            const double f=pair.phi0_f_phi0;
            const double fKf=pair.phi0_fKf_phi0;
            const double f2K=pair.phi0_f2K_phi0;
            if (world.rank()==0) printf("nf, f2, f %12.8f, %12.8f, %12.8f\n",nf, f2, f);
            if (world.rank()==0) printf("fKf, f2K  %12.8f, %12.8f\n",fKf, f2K);

            // compose the B matrix for the r12 part only
            tmp = nf + e0*f2 - e0*f*f - fKf + f2K - e0* (f2-f*f);
            if (world.rank()==0) printf("B matrix of the r12 term   %12.8f\n",tmp);
            B+=tmp;

            // now the cross terms <psi^1 | (H-E) Q f12 | phi^0>
            // note the O_{12} terms all drop out due to intermediate normalization
            const double a1=inner(psi1,pair.r12phi);
            real_function_6d eri=ERIFactory<double,6>(world).dcut(dcut);
            real_function_6d vphi=CompositeFactory<double,6,3>(world)
                    .g12(eri).ket(copy(psi1));
            const double a2=inner(phi0,vphi);
            const double a3=inner(psi1,pair.Uphi0);
            const real_function_6d Kphi=pair.KffKphi0;
            const double a4=inner(psi1,Kphi);

            if (world.rank()==0) print("new 'constant' [K,f] | phi^0> ");
            if (world.rank()==0) printf("a1-a4, %12.8f, %12.8f, %12.8f %12.8f\n",a1,a2,a3,a4);

            // compose the B matrix for the cross term r12/conventional
            tmp=e0*a1 - a2 + a3 - a4 - e0*a1;
            if (world.rank()==0) printf("B matrix of the cross term %12.8f\n",tmp);
            B+=2.0* tmp;

            // get the terms <psi^1 | (H-E) | psi^1>
            tmp=compute_B_with_U(pair);
            if (world.rank()==0) printf("B matrix of |psi^1> via U  %12.8f\n",tmp);
            tmp=compute_B_directly(i,j,psi1);
            if (world.rank()==0) printf("B matrix of the conv. term %12.8f\n",tmp);
            B+=tmp;

            if (world.rank()==0) printf("B matrix                   %12.8f\n",B);
#endif
            return B;
        }

        /// compute the B matrix by explicitly evaluating all terms of the Hamilton
        double compute_B_directly(const int i, const int j,
                const real_function_6d& fo_function) const {

            double B=0.0;

            // V_nuc, J, and K
            MADNESS_ASSERT(0);
            real_function_3d coulomb=hf->get_coulomb_potential();
            real_function_3d v_nuc=hf->get_nuclear_potential();
            real_function_3d v_total=v_nuc+coulomb;

            real_function_6d v11=CompositeFactory<double,6,3>(world)
                    .ket(copy(fo_function))
                    .V_for_particle1(copy(v_total))
                    .V_for_particle2(copy(v_total));

            const double pe=inner(fo_function,v11);
            if (world.rank()==0) printf("pe in Hylleraas  %12.8f\n\n" ,pe);

            // exchange energy
            real_function_6d Kphi=apply_exchange(fo_function,hf->orbital(i),1);
            Kphi=Kphi+apply_exchange(fo_function,hf->orbital(j),2).truncate();
            const double x=inner(fo_function,Kphi);
            if (world.rank()==0) printf("ex in Hylleraas  %12.8f\n" ,x);

            // kinetic energy expectation value
            double ke=0.0;
            for (int axis=0; axis<6; axis++) {
                real_derivative_6d D = free_space_derivative<double,6>(world, axis);
                real_function_6d dpsi = D(fo_function);
                double aa=dpsi.norm2();
                double a=0.5*aa*aa;
                ke += a;
                if (world.rank()==0) print("done with axis",axis, a);
            }
            if (world.rank()==0) printf("ke in Hylleraas  %12.8f\n" ,ke);

            // overlap <phi^1 | e1+e2 | phi^1>
            const double overlap=fo_function.norm2();
            if (world.rank()==0) printf("fo_function.norm2  %12.8f\n\n",overlap);
            const double e_contrib=overlap*overlap*this->zeroth_order_energy(i,j);
            if (world.rank()==0) printf("e-term in Hylleraas  %12.8f\n" ,e_contrib);

            B=ke+pe-e_contrib - x;
            return B;
        }

        /// compute the B matrix of the residual term using the U operator

        /// if the residual equation is solved the following relation holds:
        /// \f[ H_{12}^0 |\psi\rangle = -(U -[K,f] -J + K -E)|\phi^{(0)}\rangle\f]
        /// @param[in]	pair	the electron pair
        /// @return		the matrix element <psi^1 | H^0 | psi^1>
        double compute_B_with_U(const ElectronPair& pair) const {

        	double B=0.0;
        	const bool r0=(world.rank()==0);

        	// the first term
//        	B-=pair.phi0_gQf_phi0;
//			if (r0) printf("<phi^0 | f Q g           | phi^0>  %12.8f\n",pair.phi0_gQf_);

        	// the second term
			MADNESS_EXCEPTION("no OUKphi0 in compute_B_with_U",1);
        	const real_function_6d OO1;//OUKphi0(pair);
			double b1=inner(pair.r12phi,pair.Uphi0);
			if (r0) printf("<phi^0 | f U             | phi^0>  %12.8f\n",b1);
			double b2=inner(pair.r12phi,pair.KffKphi0);
			if (r0) printf("<phi^0 | f [K,f]         | phi^0>  %12.8f\n",b2);
        	double b3=inner(pair.r12phi,OO1);
			if (r0) printf("<phi^0 | f (O1+O2) (U-[K,f]) | phi^0>  %12.8f\n",b3);

			double b4=0.0;
            for (int k=0; k<hf->nocc(); ++k) {
                for (int l=0; l<hf->nocc(); ++l) {

                    // the function | kl>, which is only a orbital product, NOT antisymmetrized
                    real_function_6d kl=CompositeFactory<double,6,3>(world)
                                         .particle1(hf->orbital(k))
                                         .particle2(hf->orbital(l));

                    // the matrix element <kl | f12 | ij>
                    const double f_ijkl=inner(pair.r12phi,kl);
                    const double U_ijkl=inner(pair.Uphi0,kl);
                    const double K_ijkl=inner(pair.KffKphi0,kl);
                    b4+=f_ijkl*(U_ijkl-K_ijkl);
                }
            }
			if (r0) printf("<phi^0 | f O1O2 (U-[K,f])| phi^0>  %12.8f\n",b4);
        	B-=(b1-b2-b3+b4);

        	// the third term < psi^1 | U - [K,f] | phi^0 >
        	double c1=inner(pair.function,pair.Uphi0);
			if (r0) printf("<psi^1 | U               | phi^0>  %12.8f\n",c1);
        	double c2=inner(pair.function,pair.KffKphi0);
			if (r0) printf("<psi^1 | [K,f]           | phi^0>  %12.8f\n",c2);
        	B-=c1-c2;

			if (r0) printf("<phi^1 | f Q H Q f       | phi^1>  %12.8f\n",B);

        	return B;
        }

        /// apply the exchange operator on an orbital

        /// @param[in]	phi the orbital
        /// @return 	Kphi
        real_function_3d K(const real_function_3d& phi) const {
        	return hf->apply_exchange(phi);
        }

        /// apply the exchange operator on a pair function

        /// @param[in]	phi the pair function
        /// @param[in]	is_symmetric is the function symmetric wrt particle exchange
        /// @return 	(K1 + K2) |phi >
        real_function_6d K(const real_function_6d& phi, const bool is_symmetric=false) const {
        	real_function_6d result=real_factory_6d(world);
        	// loop over all orbitals of the reference
        	for (int i=0; i<hf->nocc(); ++i) {
        		real_function_6d tmp=apply_exchange(phi,hf->orbital(i),1);
        		if (is_symmetric) {
        			tmp=tmp+swap_particles(tmp);
        		} else {
        			tmp=tmp+apply_exchange(phi,hf->orbital(i),2);
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

        /// @param[in]  f   the pair function
        /// @param[in]  orbital the orbital
        /// @return     the pair function, on which the exchange operator has been applied
        real_function_6d apply_exchange(const real_function_6d& f, const real_function_3d& orbital,
                const int particle) const {

        	MADNESS_ASSERT((particle==1) or (particle==2));
            real_convolution_3d op=CoulombOperator(world,0.0001,hf->get_calc().param.econv);
            op.particle()=particle;

            if (world.rank()==0) printf("start multiplication before K at time %.1f\n",wall_time());

            // multiply the orbital to the pair function
//            real_function_6d x=(particle==1)
//            		? CompositeFactory<double,6,3>(world).ket(copy(f)).V_for_particle1(copy(orbital))
//            		: CompositeFactory<double,6,3>(world).ket(copy(f)).V_for_particle2(copy(orbital));
//            x.fill_tree().truncate();
            real_function_6d x=multiply(copy(f),copy(orbital),particle).truncate();

            // apply the Poisson operator
            if (world.rank()==0) printf("start exchange at time %.1f\n",wall_time());
            load_balance(x,false);
            x=op(x).truncate();

            // do the final multiplication with the orbital
            if (world.rank()==0) printf("start multiplication after K at time %.1f\n",wall_time());
//            real_function_6d result= (particle==1)
//            		? CompositeFactory<double,6,3>(world).ket(copy(x)).V_for_particle1(copy(orbital))
//            		: CompositeFactory<double,6,3>(world).ket(copy(x)).V_for_particle2(copy(orbital));
//            result.fill_tree().truncate().reduce_rank();
            real_function_6d result=multiply(copy(x),copy(orbital),particle).truncate();

            if (world.rank()==0) printf("end multiplication after K at time %.1f\n",wall_time());
            return result;
        }

        /// make the quantity chi_k

        /// chi is the Poisson kernel applied on an orbital product of the input function and the vector of orbitals
        /// \f[ \chi_{ki}(1) = \int dr_2 \frac{k(2) i(2)}{|r_1-r_2|} \f]
        /// @param[in]  phi   orbital phi_i
        /// @return a vector of length nocc
        std::vector<real_function_3d> make_chi(const real_function_3d& phi) const {

            const double tol=0.0;
            MADNESS_ASSERT(hf->get_calc().param.spin_restricted);
            std::vector<real_function_3d> psif = mul_sparse(world, phi, hf->get_calc().amo, tol);
            truncate(world, psif);
            psif = apply(world, *poisson, psif);
            truncate(world, psif);
            return psif;
        }

        /// make the quantity xi_k

        /// xi is chi multiplied with an orbital j
        /// \f[ \xi_{ki,j}(1) = \chi_{ki}(1) j(1) \f]
        /// @param[in]  phi_i   orbital i
        /// @param[in]  phi_j   orbital j
        /// @return a vector of length k=0,..,nocc
        std::vector<real_function_3d> make_xi(const real_function_3d& phi_i, const real_function_3d& phi_j) const {
            const double tol=0.0;
            return mul_sparse(world, phi_j, make_chi(phi_i), tol);
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
					tmp=pair.Uphi0.project_out(hf->orbital(k),0) - pair.KffKphi0.project_out(hf->orbital(k),0);
					save_function(tmp,root+"_k_"+stringify(k));
					phi_k_UK_phi0.push_back(tmp);

					tmp=pair.Uphi0.project_out(hf->orbital(k),1) - pair.KffKphi0.project_out(hf->orbital(k),1);
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
            const real_function_3d& phi_i=hf->orbital(i);
            const real_function_3d& phi_j=hf->orbital(j);

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
        real_function_6d JK1phi0_on_demand(const int i, const int j) const {
            const real_function_3d& phi_i=hf->orbital(i);
            const real_function_3d& phi_j=hf->orbital(j);

            const real_function_3d JKphi_i=J(phi_i)-K(phi_i);

            real_function_6d tmp1=CompositeFactory<double,6,3>(world)
                                 .particle1(copy(JKphi_i))
                                 .particle2(copy(phi_j));
            return tmp1;
        }

        /// return the function (J(2)-K(2)) |phi0> as on-demand function
        real_function_6d JK2phi0_on_demand(const int i, const int j) const {
            const real_function_3d& phi_i=hf->orbital(i);
            const real_function_3d& phi_j=hf->orbital(j);

            const real_function_3d JKphi_j=J(phi_j)-K(phi_j);

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


        /// multiply the given function with the 0th order Hamiltonian, exluding the 0th order energy

        /// @param[in]  f   the function we apply H^0 on
        /// @return     the function g=H^0 f, which is NOT orthogonalized against f
        real_function_6d multiply_with_0th_order_Hamiltonian(const real_function_6d& f,
        		const int i, const int j) const {

            const double eps=zeroth_order_energy(i,j);
            real_convolution_6d op_mod = BSHOperator<6>(world, sqrt(-2*eps),lo,bsh_eps);
            op_mod.modified()=true;

            real_function_3d v_total=hf->get_nuclear_potential()+hf->get_coulomb_potential();

            real_function_6d vphi=CompositeFactory<double,6,3>(world)
                                 .ket(copy(f))
                                 .V_for_particle1(copy(v_total))
                                 .V_for_particle2(copy(v_total));

            // make the tree
            vphi.fill_tree(op_mod).truncate();
            vphi.print_size("(V_nuc + J1 + J2) |ket>:  made V tree");

            vphi=(vphi-K(f,i==j)).truncate().reduce_rank();
            vphi.print_size("(V_nuc + J - K) |ket>: the tree");

            return vphi;
        }

    };
};

#endif /* MP2_H_ */

