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
  \file mp2.cc
  \brief Solves molecular MP2 equations
  \defgroup Solves molecular MP2 equations
  \ingroup examples

  The source is
  <a href=http://code.google.com/p/m-a-d-n-e-s-s/source/browse/local/trunk/src/apps/examples/mp2.cc>here</a>.


*/


#define WORLD_INSTANTIATE_STATIC_TEMPLATES
#include <mra/mra.h>
#include <mra/lbdeux.h>
#include <moldft/moldft.h>

#include <iostream>

static const bool is_helium=true;
static const double dcut=1.e-6;

using namespace madness;

template<size_t NDIM>
void save_function(World& world, const Function<double,NDIM>& pair, const std::string& name) {
    if (world.rank()==0) print("saving function",name);
    pair.print_size(name);
    archive::ParallelOutputArchive ar(world, name.c_str(), 1);
    ar & pair;
}

template<size_t NDIM>
void load_function(World& world, Function<double,NDIM>& pair, const std::string& name) {
    if (world.rank()==0) print("loading function",name);
    archive::ParallelInputArchive ar(world, name.c_str());
    ar & pair;
    pair.print_size(name);
}

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


//// according to McQuarrie
//static double he_orbital_McQuarrie(const coord_3d& r) {
//
//    // separation for 2-way decomposition (SVD; r1 -- r2)
//    const double x1=r[0];
//    const double y1=r[1];
//    const double z1=r[2];
//
//    const double r1 = sqrt(x1*x1 + y1*y1 + z1*z1 + 0.001*0.001);
//
//    const double val=exp(-(27.0/16.0)*r1);
//
//    return val;
//}


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


    /// a class holding the correlation factor for R12 theory
    class CorrelationFactor {

        World& world;
        double gamma;       ///< the correlation factor exp(-gamma r12)

    public:

        /// ctor, use negative gamma for linear correlation factor r12
        CorrelationFactor(World& world, const double& gamma) : world(world), gamma(gamma) {
            if (world.rank()==0) {
                if (gamma>0.0) print("constructed correlation factor with gamma=",gamma);
                else if (gamma==0.0) print("constructed linear correlation factor");
            }
        }

        /// return the value of the correlation factor
        double operator()(const coord_6d& r) const {
            const double rr=r12(r);
            if (gamma>0.0) return (1.0-exp(-gamma*rr))/(2.0*gamma);
            return 0.5*rr;
        }

        /// apply Kutzelnigg's regularized potential to an orbital product
        real_function_6d apply_U(const real_function_3d& phi_i, const real_function_3d& phi_j) const {
            real_function_6d result=real_factory_6d(world);

            const double eps=-1.8;
            real_convolution_6d op_mod = BSHOperator<6>(world, sqrt(-2*eps), 0.00001, 1e-6);
            op_mod.modified()=true;

            for (int axis=0; axis<3; ++axis) {
                if (world.rank()==0) print("working on axis",axis);
                real_derivative_3d D = free_space_derivative<double,3>(world, axis);
                const real_function_3d Di=(D(phi_i)).truncate();
                const real_function_3d Dj=(D(phi_j)).truncate();

                const real_function_6d u=real_factory_6d(world).functor2(U(gamma,axis)).is_on_demand();

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

            // include the purely local potential that (partially) cancels 1/r12
            if (gamma>0.0) {
                real_function_6d pair=hartree_product(phi_i,phi_j);
                real_function_6d fg3=real_factory_6d(world).functor2(fg_(gamma)).is_on_demand();
                real_function_6d mul=fg3*pair;
                mul.print_size("mul");
                mul.truncate();
                mul.print_size("mul truncated");

                result=result+mul;
            }
            result.truncate();
            result.reduce_rank();
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
            real_function_6d tmp=real_factory_6d(world).functor2(f2_(gamma)).is_on_demand();
            return tmp;
        }

        /// return fg+sth as on-demand function
        real_function_6d fg() const {
            real_function_6d tmp=real_factory_6d(world).functor2(fg_(gamma)).is_on_demand();
            return tmp;
        }

        /// return f/r as on-demand function
        real_function_6d f_over_r() const {
            real_function_6d tmp=real_factory_6d(world).functor2(f_over_r_(gamma)).is_on_demand();
            return tmp;
        }

        /// return (\nabla f)^2 as on-demand functions
        real_function_6d nablaf2() const {
            real_function_6d tmp=real_factory_6d(world).functor2(nablaf2_(gamma)).is_on_demand();
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
            nablaf2_(double gamma) : gamma(gamma) {MADNESS_ASSERT(gamma>0.0);}
            double operator()(const coord_6d& r) const {
                const double rr=r12(r);
                const double f=exp(-2.0*gamma*rr)/(4.0*gamma*gamma);
                return f;
            }
        };

        /// plot a pair function along a line
        void plot_along(const real_function_6d& pair, const std::string name) const {

            double L=FunctionDefaults<6>::get_cell_width()[0];
            coord_6d lo(0.05), hi(0.05);
            lo[0]=-L/2;
            hi[0]=L/2;
            trajectory<6> line=trajectory<6>::line2(lo,hi,600);
            madness::plot_along<6>(world,line,pair,(name+"lineplot"));
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
                if (particle_==0) tmp=hartree_product(p_[i],pf2);
                else tmp=hartree_product(pf2,p_[i]);
                sum=(sum+tmp).truncate().reduce_rank();
            }
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
            double xsq = x.sumsq();
//            if (xsq == coords_sum) {
//                return calc.current_energy;
//            }
            calc.molecule.set_all_coords(x.reshape(calc.molecule.natom(),3));
            coords_sum = xsq;

            // The below is missing convergence test logic, etc.

            // Make the nuclear potential, initial orbitals, etc.
            calc.make_nuclear_potential(world);
            calc.vnuc.print_size("vnuc");
            calc.project_ao_basis(world);

            //calc.project(world);
            if (calc.param.restart or calc.param.no_compute) {
                calc.load_mos(world);
            }
            else {
                calc.initial_guess(world);
                calc.param.restart = true;
            }

            if (not calc.param.no_compute) {
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
            }
            return calc.current_energy;
        }

        Tensor<double> gradient(const Tensor<double>& x) {

            value(x); // Ensures DFT equations are solved at this geometry
            return calc.derivatives(world);
        }

        const Calculation& get_calc() const {return calc;}

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




    /// a class for computing the first order wave function and MP2 pair energies
    class MP2 : public OptimizationTargetInterface {

        struct ElectronPair {
            int i, j;                       ///< orbitals i and j
            real_function_6d phi0;          ///< 0th order wave function (orbital product)
            real_function_6d function;      ///< pair function for a specific pair w/o correlation factor part
            real_function_6d latest_increment; ///< in the iterative residual equation
            real_function_6d r12phi;       ///< orbital product multiplied with the correlation factor

            real_function_6d Kfphi0;        ///< the function K f12 |phi^0>
            real_function_6d Uphi0;         ///< the function U |phi^0>  (U being Kutzelnigg's potential)
            real_function_6d KffKphi0;      ///< the function [K,f12] |phi^0>
            real_function_6d JKphi0;

            double first_order_correction;  ///< this plus orbital energies will yield the HF energy
            double second_order_correction; ///< this plus the HF energy will yield the MP2 correlation energy
            double second_order_energy;     ///< pair energy

            double phi0_f_phi0;             ///< <phi^0 | f12 | phi^0>
            double phi0_f2_phi0;            ///< <phi^0 | f12^2 | phi^0>
            double phi0_nf_phi0;            ///< <phi^0 | (nabla f)^2 | phi^0>
            double phi0_fovr_phi0;          ///< <phi^0 | f12/r12 | phi^0>
            double phi0_fKf_phi0;           ///< <phi^0 | f12 K f12 | phi^0>
            double phi0_f2K_phi0;           ///< <phi^0 | f12^2 K | phi^0>
            double phi0_JKf_phi0;           ///< <phi^0 | (J+K) f12 | phi^0>
            double phi0_JKOf_phi0;          ///< <phi^0 | (J+K) (O1 + O2) f12 | phi^0>
            double phi0_gOf_phi0;           ///< <phi^0 | g12 (O1 + O2) f12 | phi^0>
            double phi0_H1OOf_phi0;         ///< <phi^0 | (g12 - J+K) O1O2 f12 | phi^0>
            double phi0_gQf_phi0;           ///< <phi^0 | g12 Q12 f12 | phi^0>
            bool solved;                    ///< has the residual equation been solved for this pair?

            template <typename Archive> void serialize (Archive& ar) {
                ar & phi0_f_phi0 & phi0_f2_phi0 & phi0_nf_phi0 & phi0_fovr_phi0 & phi0_fKf_phi0 &
                	phi0_f2K_phi0 & phi0_JKf_phi0 & phi0_JKOf_phi0 & phi0_gOf_phi0 & phi0_H1OOf_phi0 &
                	phi0_gQf_phi0 & solved;
            }

            size_t size() const {
            	return phi0.size() + function.size() + latest_increment.size() + r12phi.size()
            			+ Kfphi0.size() + Uphi0.size() + KffKphi0.size();
            }

            void save_matrix_elements(World& world) const {
            	std::string name="pair_"+stringify(i)+stringify(j);
            	if (world.rank()==0) print("saving matrix elements",name);
                archive::ParallelOutputArchive ar(world, name.c_str(), 1);
                ar & *this;
            }

            void load_matrix_elements(World& world) {
            	std::string name="pair_"+stringify(i)+stringify(j);
            	if (world.rank()==0) print("loading matrix elements",name);
                archive::ParallelInputArchive ar(world, name.c_str(), 1);
                ar & *this;
            }


        };

        World& world;                           ///< the world
        HartreeFock hf;                         ///< our reference
        CorrelationFactor corrfac;              ///< correlation factor: Slater or linear

        std::vector<ElectronPair> pairs;        ///< pair functions and energies
        Projector<double,3> O1;                 ///< projector on occupied orbitals, electron 1
        Projector<double,3> O2;                 ///< projector on occupied orbitals, electron 2
        Projector<double,6> O12;                ///< projector on all pairs
        bool solved;                            ///< flag if the residual equations are already solved

    private:
        struct Intermediates {
            std::string function;      ///< pair function for a specific pair w/o correlation factor part
            std::string r12phi;       ///< orbital product multiplied with the correlation factor
            std::string latest_increment;
            std::string Kfphi0;        ///< the function K f12 |phi^0>
            std::string Uphi0;         ///< the function U |phi^0>  (U being Kutzelnigg's potential)
            std::string KffKphi0;      ///< the function [K,f12] |phi^0>
            std::string matrix_elements; ///< the matrix elements of a pair

            Intermediates() : r12phi(), latest_increment(), Kfphi0(), Uphi0(), KffKphi0() {};

            Intermediates(World& world, const std::string& filename) : function(), r12phi(), latest_increment(),
                    Kfphi0(), Uphi0(), KffKphi0(), matrix_elements() {
                std::ifstream f(filename.c_str());
                position_stream(f, "mp2");
                std::string s;

                while (f >> s) {
                    if (s == "end") break;
                    else if (s == "function") f >> function;
                    else if (s == "r12phi") f >> r12phi;
                    else if (s == "latest_increment") f >> latest_increment;
                    else if (s == "Kfphi0") f >> Kfphi0;
                    else if (s == "Uphi0") f >> Uphi0;
                    else if (s == "KffKphi0") f >> KffKphi0;
                    else if (s == "matrix_elements") f >> matrix_elements;
                    else {continue;
                    }
                    if (world.rank()==0) print("found intermediate in control file: ",s);
                }
            }

            template <typename Archive> void serialize (Archive& ar) {
                ar & function & r12phi & latest_increment & Kfphi0 & Uphi0 & KffKphi0 & matrix_elements;
            }

        };

        Intermediates intermediates;
        real_convolution_3d poisson;

    public:
        MP2(World& world, const HartreeFock& hf, const CorrelationFactor& corrfac, const std::string& input)
            : world(world)
            , hf(hf)
            , corrfac(corrfac)
            , solved(false)
            , intermediates()
            , poisson(CoulombOperator(world,0.0001,hf.get_calc().param.econv)) {

            // number of pairs:
            const int nocc=hf.nocc();
            const int npairs=nocc*(nocc+1)/2;
            pairs.resize(npairs);

            if (world.rank()==0) intermediates=Intermediates(world,input);
            world.gop.broadcast_serializable(intermediates, 0);

            std::vector<real_function_6d> all_pairs;
            for (int i=0; i<hf.nocc(); ++i) {
                for (int j=0; j<hf.nocc(); ++j) {
                    all_pairs.push_back(zeroth_order_function(i,j));
                    all_pairs.back().print_size("ij");
                }
            }
            O1=Projector<double,3>(hf.get_calc().amo,0);
            O2=Projector<double,3>(hf.get_calc().amo,1);
            O12=Projector<double,6>(all_pairs);

        }

        /// return the molecular energy as a function of the coordinates
        double value(const Tensor<double>& x)  {

            // solve the residual equations for all pairs ij
            if (not residual_equations_solved(x)) {
                for (int i=0; i<hf.nocc(); ++i) {
                    for (int j=i; j<hf.nocc(); ++j) {
                        solve_residual_equations(i,j);
                    }
                }
            }

            double energy=0.0;
            const int npair=hf.nocc()*(hf.nocc()+1)/2;
            for (int ij=0; ij<npair; ++ij) {
                // we consider this an essential:
                MADNESS_ASSERT(pairs[ij].solved);
                energy+=pairs[ij].second_order_energy;
            }
            return energy;
        }

        /// return if the equations are solved given a set of coordinates
        bool residual_equations_solved(const Tensor<double>& x) const {
            if (world.rank()==0) print("MP2::residual_equations_solved ignores the set of coordinates");
            return solved;
        }

        /// plot a pair function along a line
        void plot_along(const real_function_6d& pair, const std::string name) const {

            double L=FunctionDefaults<6>::get_cell_width()[0];
            coord_6d lo(0.05), hi(0.05);
            lo[0]=-L/2;
            hi[0]=L/2;
            trajectory<6> line=trajectory<6>::line2(lo,hi,600);
            madness::plot_along<6>(world,line,pair,(name+"lineplot"));
        }

        /// plot the MRA structure
        void plot_plane(const real_function_6d& f, const std::string filename,
                const std::string plane) const {

            coord_6d fix_coord(0.0);
            // electron 2:
            fix_coord[4]=0.1;

            f.get_impl()->print_plane(filename,"xy",fix_coord);

        }

        /// return the 0th order energy of pair ij (= sum of orbital energies)
        double zeroth_order_energy(const int i, const int j) const {
            return hf.orbital_energy(i)+hf.orbital_energy(j);
        }

        /// return the 1st order energy correction to the HF reference
        /// (=sum of orbital energies minus HF energy)
        double first_order_energy(const int i, const int j) const {
            const int ij=make_ij(i,j);
            double energy=zeroth_order_energy(i,j)+pairs[ij].first_order_correction;
            return energy;
        }

        /// return the 2nd order energy, i.e. the MP2 correlation energy of pair ij
        double second_order_energy(const int i, const int j) const {
            const int ij=make_ij(i,j);
            // we consider this an essential:
            MADNESS_ASSERT(pairs[ij].solved);
            double energy=first_order_energy(i,j)+pairs[ij].second_order_correction;
            return energy;
        }

        /// return the zeroth order wave function (=Hartree product of orbitals i and j)
        real_function_6d zeroth_order_function(const int i, const int j) const {
            real_function_6d f=hartree_product(hf.orbital(i),hf.orbital(j));
            f.set_thresh(FunctionDefaults<6>::get_thresh());
            double norm=f.norm2();
            f.scale(1.0/norm);
            return f;
        }

        /// return the first order wave function for pair ij
        real_function_6d first_order_function(const int i, const int j) const {
            int ij=make_ij(i,j);
            // we consider this an essential:
            MADNESS_ASSERT(pairs[ij].solved);
            return pairs[ij].function+pairs[ij].r12phi;
        }

        /// return the 1st order energy correction \f$ E^1 = <\phi^0 | -J + K + 1/r_{12} | \phi^0> \f$
        double compute_first_order_correction(const ElectronPair& pair) const {
            MADNESS_ASSERT(pair.phi0.is_initialized());

            // the term <phi0 | -J + K | phi0>
            double e1=-inner(pair.phi0,JK1phi0_on_demand(pair.i,pair.j))
            		  +inner(pair.phi0,JK2phi0_on_demand(pair.i,pair.j));

            // the term <phi0 | g12 | phi0>
            real_function_6d eri=ERIFactory<double,6>(world).dcut(dcut);
            real_function_6d v11=CompositeFactory<double,6,3>(world)
                    .g12(eri).ket(pair.phi0);
            e1+=inner(pair.phi0,v11);

            if (world.rank()==0) printf("first order correction e1 %12.8f\n",e1);
            return e1;

        }

        /// solve the residual equation for electron pair (i,j)
        void solve_residual_equations(const int i, const int j) {

            ElectronPair result=guess_mp1(i,j);
//            compute_second_order_correction_with_Hylleraas(i,j,result);
            compute_V(result);
            compute_second_order_correction_with_Hylleraas(i,j,result);


            // the Green's function depends on the zeroth order energy, which is the sum
            // of the orbital energies of orbitals i and j
            //  -2.0 G = (T - e_i - e_j) ^ -1
            const double eps=zeroth_order_energy(i,j);
            real_convolution_6d green = BSHOperator<6>(world, sqrt(-2*eps), 0.00001, 1e-6);


//	    if (world.rank()==0) print("computing iteratively");
//            for (int ii=0; ii<20; ++ii) {
//
//                real_function_6d vphi=multiply_with_0th_order_Hamiltonian(result.function);
//
//                /// apply the convolution
//                vphi.scale(-2.0).truncate();
//                const pairfunctionT tmp=green(vphi).truncate();
//                tmp.print_size("result of applying 0th order Hamiltonian on 1st order wave function");
//                result.function=(constant_term+tmp).truncate();
//
//                orthogonalize(result.function,zo_function);
////                orthogonalize(result.function,result.r12phi);
//
//                compute_second_order_correction_with_Hylleraas(i,j,result);
//                const std::string name="psi1_it"+stringify(ii);
//                save_function(world,result.function,name);
//                if (world.rank()==0) printf("finished iteration %2d at time %.1fs\n\n", ii, wall_time());
//
//            }

            // compute increments: psi^1 = C + GV C + GVGV C + GVGVGV C + ..
            if (world.rank()==0) print("computing increments");
            real_function_6d& latest_increment=result.latest_increment;

            for (int ii=0; ii<20; ++ii) {

                real_function_6d vphi=multiply_with_0th_order_Hamiltonian(latest_increment,i,j);
                load_balance(vphi,false);

                /// apply the convolution
                vphi.scale(-2.0).truncate();
                latest_increment=green(vphi).truncate().reduce_rank();
                latest_increment.print_size("result of applying 0th order Hamiltonian on latest increment");
                latest_increment=Q12(latest_increment);

                result.function=(result.function+latest_increment).truncate().reduce_rank();

                // compute the energy
                compute_V(result);

                // save for possible later use
                std::string name="psi1_it"+stringify(ii);
                save_function(world,result.function,name);
                name="incremental_psi1_it"+stringify(ii);
                save_function(world,latest_increment,name);

                if (world.rank()==0) printf("finished iteration %2d at time %.1fs\n\n", ii, wall_time());
	   		}
	   		
            compute_second_order_correction_with_Hylleraas(i,j,result);
            result.solved=true;

        }

        void test1(const int i, const int j) {

        	const double eps=zeroth_order_energy(i,j);
            real_convolution_6d green = BSHOperator<6>(world, sqrt(-2.0*eps), 0.00001, 1e-6);
            real_convolution_3d op = BSHOperator<3>(world, sqrt(-2.0*hf.orbital_energy(i)), 0.00001, 1e-6);

        	Vector<Translation, 6> l(Translation(1));
        	l[1]=0;
            Key<6> key(1,l);
            Key<3> key1,key2;
            key.break_apart(key1,key2);
            double a;

        	const real_function_3d phi_i=hf.orbital(i);
        	const real_function_3d phi_j=hf.orbital(j);

        	real_function_3d vloc=hf.get_coulomb_potential()+hf.get_nuclear_potential();
            real_function_3d Kphi=hf.apply_exchange(phi_i);
            const real_function_3d JKVphi_i=vloc*phi_i - Kphi;
            const real_function_3d JKVphi_j=copy(JKVphi_i);

        	real_function_6d result=hartree_product(phi_i,phi_j);
//        	real_function_6d phi0=hartree_product(JKVphi_i,JKVphi_j);
            real_function_6d tmp1=CompositeFactory<double,6,3>(world)
                                 .ket(result);

        	if(world.rank() == 0) printf("\nstarting at time %.1fs\n", wall_time());

//        	double aa=inner(result,this->JK1phi0_on_demand(i,j));
        	double aa=inner(result,tmp1);


            if(world.rank() == 0) printf("\ndone at time %.1fs\n", wall_time());
        	if (world.rank()==0) print("inner(phi0,JKphi0)",aa);
        	return;
        	real_function_6d phi1, phi0;
        	real_function_6d result1, result2;
        	real_function_6d source1;
        	real_function_6d diff;

        	// do the full comparison
        	if (1) {

//        		result1=apply_hartree_product(green,sqrt2*JKVphi_i,sqrt2*JKVphi_j);
        		result1=-2.0*green(JKVphi_i,phi_j).truncate().reduce_rank();
        		result1=result1-2.0*green(phi_i,JKVphi_j).truncate().reduce_rank();
//        		phi1=green(phi0);
            	world.gop.fence();
            	a=result1.norm2();
            	if (world.rank()==0) print("norm(green(JKVphi,JKVphi))",a);


            	diff=(result-result1);
            	a=diff.norm2();
            	if (world.rank()==0) print("norm(result-result1)",a);

        		real_function_3d r3=op(-2.0*JKVphi_i);
        		real_function_3d d3=r3-phi_i;
        		a=d3.norm2();
            	if (world.rank()==0) print("norm(op(JKVphi)-phi_i)",a);


        	}
        }


        void test(const int i, const int j) {
        }


        /// compute the matrix element <ij | g12 Q12 f12 | ij>

        /// @param[in]	i	orbital number i
        /// @param[in]	j	orbital number j
        /// @return 	the matrix element
        double compute_gQf(const ElectronPair& pair) const {

        	const int i=pair.i;
        	const int j=pair.j;
            real_function_6d ij=zeroth_order_function(i,j);
            const real_function_3d& phi_i=hf.orbital(i);
            const real_function_3d& phi_j=hf.orbital(j);

            // compute <ij| fg |ij>
            const double a=inner(ij,corrfac.f_over_r()*ij);
            if (world.rank()==0) printf("<ij | f/r          | ij>  %12.8f\n",a);

            // compute <ij| g O f | ij>
            std::vector<real_function_3d> xi_ij=make_xi(phi_i,phi_j);       // xi_{ki},j
            std::vector<real_function_3d> xi_ji=make_xi(phi_j,phi_i);       // xi_{kj},i
            double b=0.0;
            for (unsigned int k=0; k<xi_ij.size(); ++k) {

                const real_function_3d& phi_k=hf.orbital(k);
                real_function_6d tmp1=CompositeFactory<double,6,3>(world)
                                        .particle1(copy(phi_k))
                                        .particle2(copy(xi_ij[k]));
                const double b1=inner(pair.r12phi,tmp1);
                if (world.rank()==0) printf("<k xi(ik,j) | f        | ij>  %12.8f\n",b1);

                real_function_6d tmp2=CompositeFactory<double,6,3>(world)
                                        .particle1(copy(xi_ji[k]))
                                        .particle2(copy(phi_k));
                const double b2=inner(pair.r12phi,tmp2);
                if (world.rank()==0) printf("<xi(jk,i) k | f        | ij>  %12.8f\n",b2);
                b+=b1+b2;
            }
            if (world.rank()==0) printf("<ij | g (O1+O2) f      | ij>  %12.8f\n",b);

            // compute <ij| g O1O2 f |ij>
            double c=0.0;
            real_function_6d eri=ERIFactory<double,6>(world).dcut(dcut);
            for (int k=0; k<hf.nocc(); ++k) {
                for (int l=0; l<hf.nocc(); ++l) {
                    const real_function_3d& phi_k=hf.orbital(k);
                    const real_function_3d& phi_l=hf.orbital(l);

                    real_function_6d tmp1=CompositeFactory<double,6,3>(world)
                                            .g12(eri)
                                            .particle1(copy(phi_k))
                                            .particle2(copy(phi_l));

                    const double g_ijkl=inner(ij,tmp1);
                    if (world.rank()==0) printf("<kl | g                | ij>  %12.8f\n",g_ijkl);

                    real_function_6d tmp2=CompositeFactory<double,6,3>(world)
                                            .particle1(copy(phi_k))
                                            .particle2(copy(phi_l));
                    const double f_ijkl=inner(pair.r12phi,tmp2);
                    if (world.rank()==0) printf("<kl | f                | ij>  %12.8f\n",f_ijkl);
                    c+=g_ijkl*f_ijkl;
                }
            }
            if (world.rank()==0) printf("<ij | g O1 O2 f        | ij>  %12.8f\n",c);

            const double e=a-b+c;
            if (world.rank()==0) printf("<ij | g Q12 f          | ij>  %12.8f\n",e);

            return e;
        }

    private:

        /// helper function to map indices i, j to a pair index ij, with i<=j
        int make_ij(const int i, const int j) const {

            const int nocc=hf.nocc();
            MADNESS_ASSERT(i<nocc and j<nocc and i>=0 and j>=0 and i<=j);

            int ij=0;
            for (int ii=0; ii<j; ++ii) {ij+=(ii+1);}    // column index j: will yield ij=0,1,3,6
            ij+=i;                                      // row index i

            print("i,j,ij",i,j,ij);
            return ij;
        }

        /// do some load-balancing
        void load_balance(const real_function_6d& f, const bool leaf) const {

            LoadBalanceDeux<6> lb(world);
            if (leaf) lb.add_tree(f,LBCost(1.0,0.1));
            else lb.add_tree(f,LBCost(0.1,1.0));
            FunctionDefaults<6>::redistribute(world, lb.load_balance(2.0,false));
            if(world.rank() == 0) printf("redistributed at time   %.1fs\n", wall_time());

        }

        /// orthogonalize f against function g: |f> <= |f> - |g><g|f>
        void orthogonalize(real_function_6d& f, const real_function_6d& g) const {

            const double gnorm=g.norm2();
            const double fac=1.0/(gnorm*gnorm);
            Projector<double,6> P(g);   // projector on g w/o any normalization
            const real_function_6d h=(f-fac*P(f)).truncate().reduce_rank();
            const double ovlp=inner(h,g);
            const double thresh=f.thresh();
            if (ovlp>thresh and world.rank()==0) printf("ovlp in orthogonalize %12.8f\n",ovlp);
            f=h;
        }

        /// apply the strong orthogonality operator Q12 on a function f
        real_function_6d Q12(const real_function_6d& f) const {
            return (f-O1(f)-O2(f)+O12(f)).truncate().reduce_rank();
        }

        /// return the function Uphi0; load from disk if available
        real_function_6d make_Uphi0(const int i, const int j) const {
        	real_function_6d Uphi0;
            if (not intermediates.Uphi0.empty()) {
                load_function(world,Uphi0,intermediates.Uphi0);
            } else {
                Uphi0=corrfac.apply_U(hf.orbital(i),hf.orbital(j));
                Uphi0.truncate().reduce_rank();
                save_function(world,Uphi0,"Uphi0");
            }
            return Uphi0;
        }

        /// return the function [K,f] phi0; load from disk if available
        real_function_6d make_KffKphi0(const ElectronPair& pair) const {
        	real_function_6d Kfphi0;
			if (not intermediates.Kfphi0.empty()) {
				load_function(world,Kfphi0,intermediates.Kfphi0);
			} else {
				Kfphi0=apply_exchange(pair.r12phi,hf.orbital(pair.i),1);
				Kfphi0=Kfphi0 + apply_exchange(pair.r12phi,hf.orbital(pair.j),2);
				Kfphi0.truncate().reduce_rank();
				save_function(world,Kfphi0,"Kfphi0");
			}

        	real_function_6d KffKphi0;
			if (not intermediates.KffKphi0.empty()) {
				load_function(world,KffKphi0,intermediates.KffKphi0);
			} else {
				real_function_6d Kphi0=make_Kphi0(pair.i,pair.j);
				KffKphi0=Kfphi0-corrfac.f()*Kphi0;
				KffKphi0.truncate().reduce_rank();
				save_function(world,KffKphi0,"KffKphi0");
			}
			return KffKphi0;
        }

        /// compute some matrix elements that don't change during the calculation
        ElectronPair make_pair(const int i, const int j) const {

            ElectronPair pair;
            pair.i=i;
            pair.j=j;
            const bool r0=world.rank()==0;

            // some functions repeatedly used
            pair.phi0=this->zeroth_order_function(i,j);

            if (intermediates.r12phi.empty()) {
                pair.r12phi=CompositeFactory<double,6,3>(world)
                            .g12(corrfac.f()).ket(copy(pair.phi0));
                pair.r12phi.fill_tree().truncate().reduce_rank();
                save_function(world,pair.r12phi,"r12phi");
            }
            else load_function(world,pair.r12phi,intermediates.r12phi);

            const real_function_6d& fphi0=pair.r12phi;

            if (not intermediates.latest_increment.empty()) {
                load_function(world,pair.latest_increment,intermediates.latest_increment);
            }

            if (not intermediates.function.empty()) {
                load_function(world,pair.function,intermediates.function);
            }

            const real_function_6d& phi0=pair.phi0;

            pair.first_order_correction=this->compute_first_order_correction(pair);

            if (not intermediates.matrix_elements.empty()) {
            	pair.load_matrix_elements(world);
            	if (r0) {
            		print("loading matrix elements");
            		printf("<phi^0 | H^1 Q12 f12    | phi^0>  %12.8f\n",pair.phi0_gOf_phi0);
            	}
            } else {
				// some matrix elements
                real_function_6d f2phi=CompositeFactory<double,6,3>(world)
                		.g12(corrfac.f2()).ket(copy(pair.phi0));
                f2phi.fill_tree().truncate().reduce_rank();

				pair.phi0_f_phi0=inner(phi0,fphi0);
				if (r0) printf("<phi^0 | f              | phi^0>  %12.8f\n",pair.phi0_f_phi0);

				pair.phi0_JKf_phi0=inner(pair.r12phi,JK1phi0_on_demand(i,j));
				pair.phi0_JKf_phi0+=inner(pair.r12phi,JK2phi0_on_demand(i,j));
				if (r0) printf("<phi^0 | JK f12         | phi^0>  %12.8f\n",pair.phi0_JKf_phi0);

				pair.phi0_f2_phi0=inner(phi0,f2phi);
				if (r0) printf("<phi^0 | f^2            | phi^0>  %12.8f\n",pair.phi0_f2_phi0);

				pair.phi0_nf_phi0=inner(phi0,corrfac.nablaf2()*phi0);
				if (r0) printf("<phi^0 | (f')^2         | phi^0>  %12.8f\n",pair.phi0_nf_phi0);

				pair.phi0_fovr_phi0=inner(phi0,corrfac.f_over_r()*phi0);
				if (r0) printf("<phi^0 | f/r            | phi^0>  %12.8f\n",pair.phi0_fovr_phi0);

				pair.phi0_fKf_phi0=inner(fphi0,pair.Kfphi0);
				if (r0) printf("<phi^0 | f K f          | phi^0>  %12.8f\n",pair.phi0_fKf_phi0);

				real_function_3d phi_i=hf.orbital(i);
				real_function_3d phi_j=hf.orbital(j);
				real_function_3d Kphi_i=K(phi_i);
				real_function_3d Kphi_j=K(phi_j);
	            real_function_6d k1=CompositeFactory<double,6,3>(world)
	            		.particle1(copy(phi_i)).particle2(copy(Kphi_j));
				pair.phi0_f2K_phi0=inner(f2phi,k1);
	            real_function_6d k2=CompositeFactory<double,6,3>(world)
	            		.particle1(copy(Kphi_i)).particle2(copy(phi_j));
				pair.phi0_f2K_phi0+=inner(f2phi,k2);
				if (r0) printf("<phi^0 | f^2 K          | phi^0>  %12.8f\n",pair.phi0_f2K_phi0);

				pair.phi0_gOf_phi0=compute_gQf(pair);
				if (r0) printf("<phi^0 | H^1 Q12 f12    | phi^0>  %12.8f\n",pair.phi0_gOf_phi0);

				pair.save_matrix_elements(world);
            }

            if (r0) printf("done with matrix elements at time %.1fs\n\n", wall_time());
            return pair;
        }

        /// compute the first iteration of the residual equations and all intermediates
        ElectronPair guess_mp1(const int i, const int j) const {

            // the Green's function
            const double eps=zeroth_order_energy(i,j);
            real_convolution_6d green = BSHOperator<6>(world, sqrt(-2.0*eps), 0.00001, 1e-6);

            const real_function_3d& phi_i=hf.orbital(i);
            const real_function_3d& phi_j=hf.orbital(j);

            ElectronPair pair=make_pair(i,j);

            // fast return if possible
            if (pair.latest_increment.is_initialized() and pair.function.is_initialized()) return pair;

            // make the term -E^1|phi^0>
            const double e1=pair.first_order_correction;
            const real_function_6d epair=e1*pair.phi0;


			if (world.rank()==0) print("including the O1/O2 terms in make_pair");

			pair.Uphi0=make_Uphi0(i,j);
			pair.KffKphi0=make_KffKphi0(pair);

			// make the terms with high ranks and smallish trees
			real_function_6d Vpair1=(pair.Uphi0-pair.KffKphi0).truncate().reduce_rank();
			load_balance(Vpair1,true);
			real_function_6d GVpair=green(-2.0*Vpair1).truncate().reduce_rank();
			Vpair1.clear(true);
            save_function(world,GVpair,"GVpair1");


			// make the terms with J and K
			const real_function_3d jkphi_i=(J(phi_i)-K(phi_i)).truncate();
			const real_function_3d jkphi_j=(J(phi_j)-K(phi_j)).truncate();
			real_function_6d gjk1=green(jkphi_i,-2.0*phi_j).truncate().reduce_rank();
			real_function_6d gjk2=green(-2.0*phi_i,jkphi_j).truncate().reduce_rank();
			GVpair=(GVpair-gjk1-gjk2).truncate().reduce_rank();
			gjk1.clear();
			gjk2.clear();


			// make the terms with the single projectors: (O1 + O2)(-1/r12 + U -K) |phi^0>
			const real_function_6d OO1=(OUKphi0(pair) - Og12phi0(i,j)).truncate().reduce_rank();
			const real_function_6d OO2=compute_O1O2(pair);
//			pair.Uphi0.clear();
//			pair.KffKphi0.clear();

			real_function_6d Vpair2=(OO2-epair-OO1).truncate().reduce_rank();
			load_balance(Vpair2,true);
			real_function_6d tmp2=green(-2.0*Vpair2).truncate().reduce_rank();
            save_function(world,tmp2,"GVpair2");
			Vpair2.clear();

			GVpair=(GVpair+tmp2).truncate().reduce_rank();


			GVpair=Q12(GVpair);
            pair.function=GVpair;
            pair.latest_increment=copy(pair.function);
            save_function(world,GVpair,"GVpair");

            return pair;
        }


        /// compute the first iteration of the residual equations and all intermediates
        ElectronPair guess_mp1_1(const int i, const int j) const {

            // the Green's function
            const double eps=zeroth_order_energy(i,j);
            real_convolution_6d green = BSHOperator<6>(world, sqrt(-2.0*eps), 0.00001, 1e-6);

            ElectronPair pair=make_pair(i,j);

            // fast return if possible
            if (pair.latest_increment.is_initialized() and pair.function.is_initialized()) return pair;

            // get the term (J - K)|phi^0>
            pair.JKphi0=make_JKphi0(i,j);
            const real_function_6d JKpair=pair.JKphi0;

            // apply Kutzelnigg's regularized potential U_12
            pair.Uphi0=make_Uphi0(i,j);
            const real_function_6d Upair=pair.Uphi0;

            // make the term -E^1|phi^0>
            const double e1=pair.first_order_correction;
            const real_function_6d epair=e1*pair.phi0;

            // make the term [K,f]|phi^0>
            pair.KffKphi0=make_KffKphi0(pair);
            const real_function_6d K_comm=pair.KffKphi0;

            real_function_6d Vpair;
            // make all terms
            bool include_O1O2=true;
            if (include_O1O2) {

                // make the terms with the single projectors: (O1 + O2)(-1/r12 + U -K) |phi^0>
                const real_function_6d OO1=(OUKphi0(pair) - Og12phi0(i,j)).truncate().reduce_rank();
                const real_function_6d OO2=compute_O1O2(pair);

                if (world.rank()==0) print("including the O1/O2 terms in make_pair");
                Vpair=(Upair-JKpair-K_comm-epair-OO1+OO2).truncate();
            } else {
                MADNESS_ASSERT(is_helium);  // missing O1O2 term
                if (world.rank()==0) print("leaving the O1/O2 terms out of make_pair");
                Vpair=(Upair-JKpair-K_comm-epair).truncate();
            }

            real_function_6d GVpair=green(-2.0*Vpair).truncate();
            if (include_O1O2) {
                GVpair=Q12(GVpair);
            } else {
                orthogonalize(GVpair,pair.phi0);
            }

            pair.function=GVpair;
            pair.latest_increment=copy(pair.function);
            save_function(world,GVpair,"GVpair");

            return pair;
        }


        /// given a pair function, compute the pair energy using the Hylleraas functional
        double compute_second_order_correction_with_Hylleraas(const int i, const int j,
                const ElectronPair& pair) const {

            const double V=compute_V(pair);
            const double B=compute_B(i,j,pair);
//            const double B=1.0;

            const double e=2.0*V+B;
            const double e_opt=-V*V/B;

            if (world.rank()==0) {
                printf("V, B, and (%2d %2d) pair energy : %12.8f %12.8f %12.8f %12.8f\n\n",
                        i,j,V,B,e,e_opt);
            }
            return e_opt;
        }

        /// compute the V matrix for a given electron pair
        double compute_V(const ElectronPair& pair) const {

            double V=0.0;
            const double a11=inner(pair.function,JK1phi0_on_demand(pair.i,pair.j))
            				+inner(pair.function,JK2phi0_on_demand(pair.i,pair.j));
            if (world.rank()==0) printf("V2: <phi^0 | J-K        | psi^1>  %12.8f\n",a11);
            V-=a11;

            // two-electron interaction potential
            real_function_6d eri=ERIFactory<double,6>(world).dcut(dcut);
            real_function_6d vpsi1=CompositeFactory<double,6,3>(world)
                    .ket(copy(pair.function)).g12(eri);

            const double a12=inner(pair.phi0,vpsi1);
            if (world.rank()==0) printf("V2: <phi^0 | g12        | psi^1>  %12.8f\n",a12);
            V+=a12;
            if (world.rank()==0) printf("V2: <phi^0 | H^1        | psi^1>  %12.8f\n",V);

            if (world.rank()==0) printf("V2: <phi^0 | H^1 Q12 f  | phi^0>  %12.8f\n",pair.phi0_gOf_phi0);
            V+=pair.phi0_gOf_phi0;
            if (world.rank()==0) printf("V2: <phi^0 | V          | phi^1>  %12.8f\n",V);

            return V;

        }

        /// compute the B matrix for a given electron pair
        double compute_B(const int i, const int j, const ElectronPair& pair) const {

//            MADNESS_EXCEPTION("no compute_B for the time being",1);
            const double e0=zeroth_order_energy(i,j);
            const real_function_6d& phi0=pair.phi0;
            const real_function_6d& psi1=pair.function;

            // first the terms  <phi^0| f12 Q (H-E^0) Q f12 |phi^0>
            double B=0.0;
            double tmp=0.0;

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
            return B;
        }

        /// compute the B matrix by explicitly evaluating all terms of the Hamilton
        double compute_B_directly(const int i, const int j,
                const real_function_6d& fo_function) const {

            double B=0.0;

            // V_nuc, J, and K
            MADNESS_ASSERT(is_helium);
            real_function_3d coulomb=hf.get_coulomb_potential();
            real_function_3d v_nuc=hf.get_nuclear_potential();
            real_function_3d v_total=v_nuc+coulomb;

            real_function_6d v11=CompositeFactory<double,6,3>(world)
                    .ket(copy(fo_function))
                    .V_for_particle1(copy(v_total))
                    .V_for_particle2(copy(v_total));

            const double pe=inner(fo_function,v11);
            if (world.rank()==0) printf("pe in Hylleraas  %12.8f\n\n" ,pe);

            // exchange energy
            real_function_6d Kphi=apply_exchange(fo_function,hf.orbital(i),1);
            Kphi=Kphi+apply_exchange(fo_function,hf.orbital(j),2).truncate();
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
        /// therefore the B matrix is easily expressed as
        /// \f[ <\psi^1 | H^0 | \psi^1> = - <\psi^1 | (U-[K,f] -J +K -E)\phi^0> \f]
        /// @param[in]	pair	the electron pair
        /// @return		the matrix element <psi^1 | H^0 | psi^1>
        double compute_B_with_U(const ElectronPair& pair) const {
        	double B=0.0;
			const real_function_6d OO1=(OUKphi0(pair) - Og12phi0(pair.i,pair.j)).truncate().reduce_rank();
			const real_function_6d OO2=compute_O1O2(pair);

        	B+=inner(pair.function,pair.Uphi0);
        	B-=inner(pair.function,pair.KffKphi0);
        	B-=inner(pair.function,OO1);
        	B+=inner(pair.function,OO2);
        	B-=inner(pair.function,JK1phi0_on_demand(pair.i,pair.j));
        	B-=inner(pair.function,JK2phi0_on_demand(pair.i,pair.j));
        	B*=-1.0;
        	return B;
        }

        /// apply the exchange operator a on orbital

        /// @param[in]	phi the orbital
        /// @return 	Kphi
        real_function_3d K(const real_function_3d& phi) const {
        	return hf.apply_exchange(phi);
        }

        /// apply the Coulomb operator a on orbital

        /// @param[in]	phi the orbital
        /// @return 	Jphi
        real_function_3d J(const real_function_3d& phi) const {
        	return (hf.get_coulomb_potential()*phi).truncate();
        }

        /// apply the exchange operator on f

        /// @param[in]  f   the pair function
        /// @param[in]  orbital the orbital
        /// @return     the pair function, on which the exchange operator has been applied
        real_function_6d apply_exchange(const real_function_6d& f, const real_function_3d& orbital,
                const int particle) const {

            real_convolution_3d op=CoulombOperator(world,0.0001,hf.get_calc().param.econv);
            op.particle()=particle;

            real_function_6d x=multiply(f,orbital,particle).truncate();
            x=op(x);
            x=multiply(x,orbital,particle);
            x.truncate().reduce_rank();

            return x;
        }

        /// make the quantity chi_k

        /// chi is the Poisson kernel applied on an orbital product of the input function and the vector of orbitals
        /// \f[ \chi_{ki}(1) = \int dr_2 \frac{k(2) i(2)}{|r_1-r_2|} \f]
        /// @param[in]  phi   orbital phi_i
        /// @return a vector of length nocc
        std::vector<real_function_3d> make_chi(const real_function_3d& phi) const {

            const double tol=0.0;
            MADNESS_ASSERT(hf.get_calc().param.spin_restricted);
            std::vector<real_function_3d> psif = mul_sparse(world, phi, hf.get_calc().amo, tol);
            truncate(world, psif);
            psif = apply(world, poisson, psif);
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

        /// compute the term (O_1 + O_2) g12 |phi0>

        /// this term is computed with the \f$\chi\f$ quantity (cf make_chi())
        /// \f[ \left(O_1 + O_2\right)\frac{1}{r_{12}} | ij> =
        ///        \sum_k k(1) j(2) \chi_{ki}(2) + \sum_k k(2) i(1) \chi_{kj}(1) \f]
        real_function_6d Og12phi0(const int i, const int j) const {

            real_function_6d result=real_factory_6d(world);
            const double tol=0.0;
            const real_function_3d& phi_i=hf.orbital(i);
            const real_function_3d& phi_j=hf.orbital(j);

            std::vector<real_function_3d> ichi_kj = mul_sparse(world, phi_i, make_chi(phi_j), tol);
            for (int k=0; k<hf.nocc(); ++k) {
                result=(result+hartree_product(hf.orbital(k),ichi_kj[k])).truncate().reduce_rank();
            }
            ichi_kj.clear();

            std::vector<real_function_3d> jchi_ki = mul_sparse(world, phi_j, make_chi(phi_i), tol);
            for (int k=0; k<hf.nocc(); ++k) {
                result=(result+hartree_product(jchi_ki[k],hf.orbital(k))).truncate().reduce_rank();
            }
            jchi_ki.clear();

            return result;
        }

        /// compute the term (O_1 + O_2) (U - [K,f]) |phi0>
        real_function_6d OUKphi0(const ElectronPair& pair) const {
            pair.Uphi0.verify();
            pair.KffKphi0.verify();
            const real_function_6d UK=(pair.Uphi0 - pair.KffKphi0).truncate().reduce_rank();
            const real_function_6d result=(O1(UK) + O2(UK)).truncate().reduce_rank();
            return result;
        }

        /// compute the term  O1O2[H^0,f12]|phi0>  =  |kl >< kl | f12 | ij> (e_kl - e_ij)
        real_function_6d compute_O1O2(const ElectronPair& pair) const {
            const int i=pair.i;
            const int j=pair.j;
            const double e_ij=zeroth_order_energy(i,j);
            real_function_6d result=real_factory_6d(world);
            for (int k=0; k<hf.nocc(); ++k) {
                for (int l=0; l<hf.nocc(); ++l) {

                    // the function | kl>, which is only a orbital product, NOT antisymmetrized
                    real_function_6d kl=hartree_product(hf.orbital(k),hf.orbital(l));
                    const double norm=kl.norm2();
                    kl.scale(1.0/norm);

                    // the matrix element <kl | f12 | ij>
                    const double f_ijkl=inner(kl,pair.r12phi);
                    print("i,j,k,l,f_ijkl",i,j,k,l,f_ijkl);
                    const double e_kl=zeroth_order_energy(k,l);

                    result=(result+(f_ijkl*(e_kl-e_ij))*kl).truncate().reduce_rank();
                }
            }
            return result;
        }

        /// apply the operator J on the reference; J |phi^0>

        /// the application of the operators naturally factorizes as in Hartree-Fock
        /// @param[in]  i   index of orbital i
        /// @param[in]  j   index of orbital j
        /// @return     the function J |phi^0>
        real_function_6d make_Jphi0(const int i, const int j) const {

            const real_function_3d& phi_i=hf.orbital(i);
            const real_function_3d& phi_j=hf.orbital(j);

            // make the pair function (J(1) + J(2) ) |ij>
            const real_function_3d coulomb=hf.get_coulomb_potential();
            const real_function_3d Jphi_i=coulomb*phi_i;
            const real_function_3d Jphi_j=coulomb*phi_j;

//            const real_function_6d Jphi=(hartree_product(Jphi_i,phi_j) + hartree_product(phi_i,Jphi_j)).truncate().reduce_rank();

            real_function_6d tmp1=CompositeFactory<double,6,3>(world)
                                 .particle1(copy(Jphi_i))
                                 .particle2(copy(phi_j));
            real_function_6d tmp2=CompositeFactory<double,6,3>(world)
                                 .particle1(copy(phi_i))
                                 .particle2(copy(Jphi_j));

            const double eps=zeroth_order_energy(i,j);
            real_convolution_6d op_mod = BSHOperator<6>(world, sqrt(-2*eps), 0.00001, 1e-6);
            op_mod.modified()=true;

            tmp1.fill_tree(op_mod);
            tmp2.fill_tree(op_mod);

            const real_function_6d Jphi=(tmp1+tmp2).truncate().reduce_rank();
            return Jphi;
        }

        /// apply the operator K on the reference; K |phi^0>

        /// the application of the operators naturally factorizes as in Hartree-Fock
        /// \f[ K(1)|ij> = \sum_k k(1)j(2) \chi_{ki}(1) \f]
        /// \f[ K(2)|ij> = \sum_k k(2)i(1) \chi_{kj}(2) \f]
        /// @param[in]  i   index of orbital i
        /// @param[in]  j   index of orbital j
        /// @return     the function (K(1)+K(2))|phi^0>
        real_function_6d make_Kphi0(const int i, const int j) const {

            const real_function_3d& phi_i=hf.orbital(i);
            const real_function_3d& phi_j=hf.orbital(j);

            // make the quantities chi and multiply them with the orbitals k
            MADNESS_ASSERT(hf.get_calc().param.spin_restricted);
            std::vector<real_function_3d> chi1 = mul(world, hf.get_calc().amo, make_chi(phi_i));
            std::vector<real_function_3d> chi2 = mul(world, hf.get_calc().amo, make_chi(phi_j));

            // perform the sum over k
            compress(world, chi1);
            compress(world, chi2);
            real_function_3d k1 = real_factory_3d(world);
            real_function_3d k2 = real_factory_3d(world);
            k1.compress();
            k2.compress();
            for(int k=0; k<hf.nocc(); ++k) {
                k1.gaxpy(1.0, chi1[k], hf.get_calc().aocc[k], false);
                k2.gaxpy(1.0, chi2[k], hf.get_calc().aocc[k], false);
            }
            world.gop.fence();
            chi1.clear();
            chi2.clear();

            // make the pair function (K(1) + K(2) ) |ij>
            const real_function_6d Kphi=(hartree_product(k1,phi_j) + hartree_product(phi_i,k2)).truncate().reduce_rank();
            return Kphi;
        }

        /// apply the operator (J-K) on the reference; K |phi^0>

        /// the application of the operators naturally factorizes as in Hartree-Fock
        /// \f[ K(1)|ij> = \sum_k k(1)j(2) \chi_{ki}(1) \f]
        /// \f[ K(2)|ij> = \sum_k k(2)i(1) \chi_{kj}(2) \f]
        /// @param[in]  i   index of orbital i
        /// @param[in]  j   index of orbital j
        /// @return     the function (J1+J2-K1-K2)|phi^0>
        real_function_6d make_JKphi0(const int i, const int j) const {
            const real_function_3d& phi_i=hf.orbital(i);
            const real_function_3d& phi_j=hf.orbital(j);

            const real_function_3d JKphi_i=(J(phi_i)-K(phi_i)).truncate();
            const real_function_3d JKphi_j=(J(phi_j)-K(phi_j)).truncate();

            const double eps=zeroth_order_energy(i,j);
            real_convolution_6d op_mod = BSHOperator<6>(world, sqrt(-2*eps), 0.00001, 1e-6);
            op_mod.modified()=true;

            const real_function_6d tmp1=(hartree_product(JKphi_i,phi_j,op_mod));
            const real_function_6d tmp2=(hartree_product(phi_i,JKphi_j,op_mod));

            const real_function_6d JKphi0=(tmp1+tmp2).truncate().reduce_rank();
            return JKphi0;
        }

        /// return the function (J(1)-K(1)) |phi0> as on-demand function
        real_function_6d JK1phi0_on_demand(const int i, const int j) const {
            const real_function_3d& phi_i=hf.orbital(i);
            const real_function_3d& phi_j=hf.orbital(j);

            const real_function_3d JKphi_i=J(phi_i)-K(phi_i);

            real_function_6d tmp1=CompositeFactory<double,6,3>(world)
                                 .particle1(copy(JKphi_i))
                                 .particle2(copy(phi_j));
            return tmp1;
        }

        /// return the function (J(2)-K(2)) |phi0> as on-demand function
        real_function_6d JK2phi0_on_demand(const int i, const int j) const {
            const real_function_3d& phi_i=hf.orbital(i);
            const real_function_3d& phi_j=hf.orbital(j);

            const real_function_3d JKphi_j=J(phi_j)-K(phi_j);

            real_function_6d tmp=CompositeFactory<double,6,3>(world)
                                 .particle1(copy(phi_i))
                                 .particle2(copy(JKphi_j));
            return tmp;
        }

        /// multiply the given function with the 0th order Hamiltonian, exluding the 0th order energy

        /// @param[in]  f   the function we apply H^0 on
        /// @return     the function g=H^0 f, which is NOT orthogonalized against f
        real_function_6d multiply_with_0th_order_Hamiltonian(const real_function_6d& f, const int i, const int j) const {

            const double eps=zeroth_order_energy(i,j);
            real_convolution_6d op_mod = BSHOperator<6>(world, sqrt(-2*eps), 0.00001, 1e-6);
            op_mod.modified()=true;

            real_function_3d v_total=hf.get_nuclear_potential()+hf.get_coulomb_potential();

            real_function_6d vphi=CompositeFactory<double,6,3>(world)
                                 .ket(copy(f))
                                 .V_for_particle1(copy(v_total))
                                 .V_for_particle2(copy(v_total));

            // make the tree
            vphi.fill_tree(op_mod).truncate();
            vphi.print_size("(V_nuc + J1 + J2) |ket>:  made V tree");

            // add exchange
            for (int k=0; k<hf.nocc(); ++k) {
                vphi=vphi-apply_exchange(f,hf.orbital(k),1);
                vphi=vphi-apply_exchange(f,hf.orbital(k),2);
                vphi.truncate().reduce_rank();
            }
            vphi.print_size("(V_nuc + J - K) |ket>: the tree");

            return vphi;
        }

    };
};

int main(int argc, char** argv) {
    initialize(argc, argv);
    World world(MPI::COMM_WORLD);
    startup(world,argc,argv);
    std::cout.precision(6);

    int i=0,j=0;

    // get parameters form input file
    Calculation calc(world,"input");
    TensorType tt=TT_2D;

	// find the pair to compute
    std::ifstream f("input");
    position_stream(f, "mp2");
    std::string s;

    while (f >> s) {
        if (s == "end") break;
        else if (s == "pair") f >> i >> j;
        else {continue;
        }
    }


    // get command line parameters (overrides input file)
    for(int ii = 1; ii < argc; ii++) {
        const std::string arg=argv[ii];

        // break parameters into key and val
        size_t pos=arg.find("=");
        std::string key=arg.substr(0,pos);
        std::string val=arg.substr(pos+1);

        if (key=="size") calc.param.L=atof(val.c_str());               // usage: size=10
        if (key=="k") calc.param.k=atoi(val.c_str());                  // usage: k=5
        if (key=="i") i=atoi(val.c_str());                  // usage: k=5
        if (key=="j") j=atoi(val.c_str());                  // usage: k=5
        if (key=="thresh") calc.param.econv=atof(val.c_str());        // usage: thresh=1.e-3
        if (key=="TT") {
            if (val=="TT_2D") tt=TT_2D;
            else if (val=="TT_3D") tt=TT_3D;
            else if (val=="TT_FULL") tt=TT_FULL;
            else {
                print("arg",arg, "key",key,"val",val);
                MADNESS_EXCEPTION("confused tensor type",0);
            }
        }
    }

    // actually set the FunctionDefaults
    calc.set_protocol<6>(world,calc.param.econv);
    calc.set_protocol<3>(world,calc.param.econv*0.01);
    calc.molecule.set_eprec(std::min(calc.param.econv,1.e-6));
    FunctionDefaults<6>::set_tensor_type(tt);
    FunctionDefaults<6>::set_apply_randomize(true);


    if (world.rank()==0) {
        print("polynomial order:  ", FunctionDefaults<6>::get_k());
        print("threshold 3D:      ", FunctionDefaults<3>::get_thresh());
        print("threshold 6D:      ", FunctionDefaults<6>::get_thresh());
        print("cell size:         ", FunctionDefaults<6>::get_cell_width()[0]);
        print("truncation mode:   ", FunctionDefaults<6>::get_truncate_mode());
        print("tensor type:       ", FunctionDefaults<6>::get_tensor_type());
        print("");
        print("orthogonalization  ", OrthoMethod());
        print("facReduce          ", GenTensor<double>::fac_reduce());
        print("max displacement   ", Displacements<6>::bmax_default());
        print("apply randomize    ", FunctionDefaults<6>::get_apply_randomize());
        print("world.size()       ", world.size());
        print("");
    }

    if (world.rank()==0) calc.param.print(world);
    if (world.rank()==0) print("\nsize consistency of the 6d green's function?\n");


    HartreeFock hf(world,calc);
    hf.value();

    CorrelationFactor f12(world,1.0);

    if(world.rank() == 0) printf("\nstarting at time %.1fs\n", wall_time());
    if(world.rank() == 0) printf("\ncomputing pair (%d,%d)\n\n", i,j);

    MP2 mp2(world,hf,f12,"input");

//    mp2.value(calc.molecule.get_all_coords());
//    mp2.test(0,0);
    mp2.solve_residual_equations(i,j);

    if(world.rank() == 0) printf("\nfinished at time %.1fs\n\n", wall_time());
    world.gop.fence();
    finalize();

    return 0;
}
