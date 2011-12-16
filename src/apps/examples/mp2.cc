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
  \defgroup MP2
  \ingroup examples

  The source is
  <a href=http://code.google.com/p/m-a-d-n-e-s-s/source/browse/local/trunk/src/apps/examples/helium_mp2.cc>here</a>.


*/


#define WORLD_INSTANTIATE_STATIC_TEMPLATES
#include <mra/mra.h>
#include <mra/operator.h>
#include <mra/funcplot.h>
#include <mra/lbdeux.h>
#include <moldft/moldft.h>

#include <iostream>

static const bool is_helium=true;
static const double dcut=1.e-6;



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
//                double norm;
//
//                const real_function_6d Dij1=(hartree_product(Di,phi_j,op_mod));
//                Dij1.print_size("Dij1");
//                norm=Dij1.norm2();
//                print("Dij1.norm2",norm);
//                const real_function_6d UDij1=Dij1*u;
//                UDij1.print_size("UDij1");
//                norm=UDij1.norm2();
//                print("UDij1.norm2",norm);

                real_function_6d tmp1=CompositeFactory<double,6,3>(world)
                            .g12(u.get_impl()).particle1(copy(Di).get_impl()).particle2(copy(phi_j).get_impl());
                tmp1.fill_tree(op_mod).truncate();
                real_function_6d tmp2=CompositeFactory<double,6,3>(world)
                            .g12(u.get_impl()).particle1(copy(phi_i).get_impl()).particle2(copy(Dj).get_impl());
                tmp2.fill_tree(op_mod).truncate();
                if (world.rank()==0) print("done with fill_tree");
//                tmp.print_size("tmp before truncation");
//                norm=tmp.norm2();
//                print("tmp.norm2",norm);
//                tmp.truncate();
//                tmp.print_size("tmp after truncation");
//                norm=tmp.norm2();
//                print("tmp.norm2",norm);
//                real_function_6d diff=UDij1-tmp;
//                double error=diff.norm2();
//                print("diff norm ",error);

                result=result+(tmp1-tmp2).truncate();
                result.truncate().reduce_rank();
//                result=result+(u*Dij).truncate();
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
        const Function<T,NDIM> p;

    public:
        Projector(const Function<T,NDIM>& f) : particle_(1), p(f) {
        }

        int& particle() {return particle_;}
        const int& particle() const {return particle_;}

        /// project f on p: |result> =  | p><p | f>
        template<std::size_t FDIM>
        Function<T,FDIM> operator()(const Function<T,FDIM>& f) const {
            MADNESS_ASSERT(FDIM==NDIM); // for now
            const double ovlp=inner(f,p);
            const Function<T,FDIM> pp=ovlp*p;
            return pp;
        }
    };

    class HartreeFock : public OptimizationTargetInterface {
        World& world;
        Calculation calc;
        mutable double coords_sum;     // sum of square of coords at last solved geometry
        mutable double E; //< Current energy

        // save the Coulomb potential
        mutable functionT coulomb;

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
            calc.project_ao_basis(world);

            //calc.project(world);
            if (calc.param.restart) {
                calc.load_mos(world);
            }
            else {
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
        real_function_3d apply_exchange(real_function_3d phi) const {
            MADNESS_ASSERT(is_helium);
            const int i=0;
            real_function_3d iphi=orbital(i)*phi;
            real_function_3d result=orbital(i)*calc.make_coulomb_potential(iphi);
            return result;
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
            real_function_6d phi0;          ///< 0th order wave function (orbital product)
            real_function_6d function;      ///< pair function for a specific pair w/o correlation factor part
            real_function_6d latest_increment; ///< in the iterative residual equation
            real_function_6d r12phi;       ///< orbital product multiplied with the correlation factor

            real_function_6d Kfphi0;        ///< the function K f12 |phi^0>
            real_function_6d Uphi0;         ///< the function U |phi^0>  (U being Kutzelnigg's potential)
            real_function_6d KffKphi0;      ///< the function [K,f12] |phi^0>

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
            bool solved;                    ///< has the residual equation been solved for this pair?
        };

        World& world;                           ///< the world
        HartreeFock hf;                         ///< our reference
        CorrelationFactor corrfac;              ///< correlation factor: Slater or linear

        std::vector<ElectronPair> pairs;        ///< pair functions and energies
        bool solved;                            ///< flag if the residual equations are already solved

        struct Intermediates {
            std::string function;      ///< pair function for a specific pair w/o correlation factor part
            std::string r12phi;       ///< orbital product multiplied with the correlation factor
            std::string latest_increment;
            std::string Kfphi0;        ///< the function K f12 |phi^0>
            std::string Uphi0;         ///< the function U |phi^0>  (U being Kutzelnigg's potential)
            std::string KffKphi0;      ///< the function [K,f12] |phi^0>

            Intermediates() : r12phi(), latest_increment(), Kfphi0(), Uphi0(), KffKphi0() {};

            Intermediates(World& world, const std::string& filename) : function(), r12phi(), latest_increment(),
                    Kfphi0(), Uphi0(), KffKphi0() {
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
                    else {continue;
                    }
                    if (world.rank()==0) print("found intermediate in control file: ",s);
                }
            }
            template <typename Archive> void serialize (Archive& ar) {
                ar & function & r12phi & latest_increment & Kfphi0 & Uphi0 & KffKphi0;
            }

        };

        Intermediates intermediates;
        real_convolution_3d poisson;

        static const double dcut=1.e-6;

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

        /// return the 1st order energy correction E^1 = <phi^0 | -J + K + 1/r12 | phi^0>
        double compute_first_order_correction(const ElectronPair& pair) const {
            MADNESS_ASSERT(pair.phi0.is_initialized());
            const real_function_6d& phi0=pair.phi0;

            // J and K
            MADNESS_ASSERT(is_helium);  // scale 0.5*J, leaving out K
            functionT J=-0.5*hf.get_coulomb_potential();

            real_function_6d eri=ERIFactory<double,6>(world).dcut(dcut);
            real_function_6d v11=CompositeFactory<double,6,3>(world)
                   .V_for_particle1(copy(J).get_impl()).V_for_particle2(copy(J).get_impl())
                   .g12(eri.get_impl()).ket(copy(phi0).get_impl());

            const double fo_energy=inner(phi0,v11);
            return fo_energy;
        }

        /// solve the residual equation for electron pair (i,j)
        void solve_residual_equations(const int i, const int j) {

            ElectronPair result=guess_mp1(i,j);
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

                real_function_6d vphi=multiply_with_0th_order_Hamiltonian(latest_increment);
                load_balance(vphi,false);

                /// apply the convolution
                vphi.scale(-2.0).truncate();
                latest_increment=green(vphi).truncate().reduce_rank();
                latest_increment.print_size("result of applying 0th order Hamiltonian on latest increment");
                orthogonalize(latest_increment,result.phi0);

                result.function=(result.function+latest_increment).truncate().reduce_rank();

                // compute the energy
                compute_second_order_correction_with_Hylleraas(i,j,result);

                // save for possible later use
                std::string name="psi1_it"+stringify(ii);
                save_function(world,result.function,name);
                name="incremental_psi1_it"+stringify(ii);
                save_function(world,latest_increment,name);

                if (world.rank()==0) printf("finished iteration %2d at time %.1fs\n\n", ii, wall_time());
	   		}


            result.solved=true;

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

        /// compute some matrix elements that don't change during the calculation
        ElectronPair make_pair(const int i, const int j) const {

            ElectronPair pair;
            const bool r0=world.rank()==0;

            // some functions repeatedly used
            pair.phi0=this->zeroth_order_function(i,j);
            pair.first_order_correction=this->compute_first_order_correction(pair);

            if (intermediates.r12phi.empty()) {
//                pair.r12phi=corrfac.f()*pair.phi0;
                real_function_6d bla=corrfac.f()*pair.phi0;
                pair.r12phi=CompositeFactory<double,6,3>(world)
                            .g12(corrfac.f().get_impl()).ket(copy(pair.phi0).get_impl());
                pair.r12phi.fill_tree().truncate().reduce_rank();
                save_function(world,pair.r12phi,"r12phi");
            }
            else load_function(world,pair.r12phi,intermediates.r12phi);

            const real_function_6d& fphi0=pair.r12phi;
            const real_function_3d& phi_i=hf.orbital(i);
            const real_function_3d& phi_j=hf.orbital(j);
            const real_function_3d Kphi_i=hf.apply_exchange(phi_i);
            const real_function_3d Kphi_j=hf.apply_exchange(phi_j);
            const real_function_6d Kphi0=(hartree_product(Kphi_i,phi_j)
                    + hartree_product(phi_i,Kphi_j)).truncate();


            if (not intermediates.Kfphi0.empty()) {
                load_function(world,pair.Kfphi0,intermediates.Kfphi0);
            } else {
                pair.Kfphi0=apply_exchange(fphi0,hf.orbital(i),1);
                pair.Kfphi0=pair.Kfphi0 + apply_exchange(fphi0,hf.orbital(j),2);
                pair.Kfphi0.truncate().reduce_rank();
                save_function(world,pair.Kfphi0,"Kfphi0");
            }

            if (not intermediates.KffKphi0.empty()) {
                load_function(world,pair.KffKphi0,intermediates.KffKphi0);
            } else {
                pair.KffKphi0=pair.Kfphi0-corrfac.f()*Kphi0;
                pair.KffKphi0.truncate().reduce_rank();
                save_function(world,pair.KffKphi0,"KffKphi0");
            }

            if (not intermediates.Uphi0.empty()) {
                load_function(world,pair.Uphi0,intermediates.Uphi0);
            } else {
                pair.Uphi0=corrfac.apply_U(phi_i,phi_j);
                pair.Uphi0.truncate().reduce_rank();
                save_function(world,pair.Uphi0,"Uphi0");
            }

            if (not intermediates.latest_increment.empty()) {
                load_function(world,pair.latest_increment,intermediates.latest_increment);
            }

            if (not intermediates.function.empty()) {
                load_function(world,pair.function,intermediates.function);
            }

            load_balance(pair.Uphi0,false);

            const real_function_6d& phi0=pair.phi0;

            // some matrix elements
            pair.phi0_f_phi0=inner(phi0,fphi0);
            if (r0) printf("<phi^0 | f       | phi^0>  %12.8f\n",pair.phi0_f_phi0);

            pair.phi0_f2_phi0=inner(phi0,corrfac.f2()*phi0);
            if (r0) printf("<phi^0 | f^2     | phi^0>  %12.8f\n",pair.phi0_f2_phi0);

            pair.phi0_nf_phi0=inner(phi0,corrfac.nablaf2()*phi0);
            if (r0) printf("<phi^0 | (f')^2  | phi^0>  %12.8f\n",pair.phi0_nf_phi0);

            pair.phi0_fovr_phi0=inner(phi0,corrfac.f_over_r()*phi0);
            if (r0) printf("<phi^0 | f/r     | phi^0>  %12.8f\n",pair.phi0_fovr_phi0);

            pair.phi0_fKf_phi0=inner(fphi0,pair.Kfphi0);
            if (r0) printf("<phi^0 | f K f   | phi^0>  %12.8f\n",pair.phi0_fKf_phi0);

            pair.phi0_f2K_phi0=inner(Kphi0,corrfac.f2()*phi0);
            if (r0) printf("<phi^0 | f^2 K   | phi^0>  %12.8f\n",pair.phi0_f2K_phi0);

            if (r0) printf("done with matrix elements at time %.1fs\n\n", wall_time());
            return pair;
        }

        /// compute the first iteration of the residual equations and all intermediates
        ElectronPair guess_mp1(const int i, const int j) const {

            // the Green's function
            const double eps=zeroth_order_energy(i,j);
            real_convolution_6d green = BSHOperator<6>(world, sqrt(-2.0*eps), 0.00001, 1e-6);

            ElectronPair pair=make_pair(i,j);

            // fast return if possible
            if (pair.latest_increment.is_initialized() and pair.function.is_initialized()) return pair;

            // make the term (-J + K)|phi^0>
            const real_function_3d phi_i=hf.orbital(i);
            const real_function_3d phi_j=hf.orbital(j);
            const real_function_3d coulomb=hf.get_coulomb_potential();

            const real_function_3d JKphi_i=-0.5*coulomb*phi_i;
            const real_function_3d JKphi_j=-0.5*coulomb*phi_j;
            const real_function_6d JKpair=(hartree_product(JKphi_i,phi_j) + hartree_product(phi_i,JKphi_j)).truncate();

            // apply Kutzelnigg's regularized potential U_12
            const real_function_6d Upair=pair.Uphi0;

            // make the term -E^1|phi^0>
            const double e1=pair.first_order_correction;
            const real_function_6d epair=e1*pair.phi0;

            // make the term [K,f]|phi^0> 
            const real_function_6d K_comm=pair.KffKphi0;

            // make all terms
            real_function_6d Vpair=(JKpair+Upair-K_comm-epair).truncate();
            real_function_6d GVpair=green(-2.0*Vpair).truncate();
            orthogonalize(GVpair,pair.phi0);

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
            const real_function_6d& phi0=pair.phi0;
            const real_function_6d& psi1=pair.function;
            const double e1=pair.first_order_correction;

            // two-electron interaction potential
            real_function_6d eri=ERIFactory<double,6>(world).dcut(dcut);
            MADNESS_ASSERT(is_helium);
            real_function_3d coulomb=-0.5*hf.get_coulomb_potential();

            real_function_6d vpsi1=CompositeFactory<double,6,3>(world)
                    .V_for_particle1(copy(coulomb).get_impl())
                    .V_for_particle2(copy(coulomb).get_impl())
                    .ket(copy(psi1).get_impl()).g12(eri.get_impl());

            const double a1=inner(phi0,vpsi1);
            if (world.rank()==0) printf("<phi^0 | H^1     | psi^1>  %12.8f\n",a1);
            V+=a1;

            real_function_6d JKphi0=CompositeFactory<double,6,3>(world)
                    .ket(copy(phi0).get_impl())
                    .V_for_particle1(copy(coulomb).get_impl())
                    .V_for_particle2(copy(coulomb).get_impl());
            const double a2=inner(pair.r12phi,JKphi0);
            if (world.rank()==0) printf("<phi^0 | (J+K) f | phi^0>  %12.8f\n",a2);
            V+=a2;

//            const double a3=inner(phi0,corrfac.f_over_r()*phi0);
            const double a3=pair.phi0_fovr_phi0;
            if (world.rank()==0) printf("<phi^0 |  f/r    | phi^0>  %12.8f\n",a3);
            V+=a3;

//            const double a4=inner(phi0,corrfac.f()*phi0);
            const double a4=pair.phi0_f_phi0;
            if (world.rank()==0) printf("<phi^0 |  f      | phi^0>  %12.8f\n",a4);
            V-=e1*a4;

            if (world.rank()==0) printf("<phi^0 |  V      | phi^1>  %12.8f\n",V);
            return V;
        }

        /// compute the B matrix for a given electron pair
        double compute_B(const int i, const int j, const ElectronPair& pair) const {


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
                    .g12(eri.get_impl()).ket(copy(psi1).get_impl());
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
            functionT coulomb=hf.get_coulomb_potential();
            functionT v_nuc=hf.get_nuclear_potential();
            functionT v_total=v_nuc+coulomb;

            real_function_6d v11=CompositeFactory<double,6,3>(world)
                    .ket(copy(fo_function).get_impl())
                    .V_for_particle1(copy(v_total).get_impl())
                    .V_for_particle2(copy(v_total).get_impl());

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

        /// multiply the given function with the 0th order Hamiltonian, exluding the 0th order energy

        /// @param[in]  f   the function we apply H^0 on
        /// @return     the function g=H^0 f, which is NOT orthogonalized against f
        real_function_6d multiply_with_0th_order_Hamiltonian(const real_function_6d& f) const {

            // for estimating the MRA structure of V*phi
            MADNESS_ASSERT(is_helium);
            const int i=0;
            const int j=0;

            const double eps=zeroth_order_energy(i,j);
            real_convolution_6d op_mod = BSHOperator<6>(world, sqrt(-2*eps), 0.00001, 1e-6);
            op_mod.modified()=true;

            functionT v_total=hf.get_nuclear_potential()+hf.get_coulomb_potential();

            real_function_6d vphi=CompositeFactory<double,6,3>(world)
                                 .ket(copy(f).get_impl())
                                 .V_for_particle1(copy(v_total).get_impl())
                                 .V_for_particle2(copy(v_total).get_impl());

            // make the tree
            vphi.fill_tree(op_mod).truncate();
            vphi.print_size("(V_nuc + J1 + J2) |ket>:  made V tree");

            // add exchange
            vphi=vphi-apply_exchange(f,hf.orbital(i),1);
            vphi=vphi-apply_exchange(f,hf.orbital(j),2);
            vphi.truncate().reduce_rank();
            vphi.print_size("(V_nuc + J - K) |ket>: the tree");

            return vphi;
        }

    };
}

int main(int argc, char** argv) {
    initialize(argc, argv);
    World world(MPI::COMM_WORLD);
    startup(world,argc,argv);
    std::cout.precision(6);


    // get parameters form input file
    Calculation calc(world,"input");
    TensorType tt=TT_2D;
    bool restart=false;

    // get command line parameters (overrides input file)
    for(int i = 1; i < argc; i++) {
        const std::string arg=argv[i];

        // break parameters into key and val
        size_t pos=arg.find("=");
        std::string key=arg.substr(0,pos);
        std::string val=arg.substr(pos+1);

        if (key=="restart") restart=true;                             // usage: size=10
        if (key=="size") calc.param.L=atof(val.c_str());               // usage: size=10
        if (key=="k") calc.param.k=atoi(val.c_str());                  // usage: k=5
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
  
    MP2 mp2(world,hf,f12,"input");

//    mp2.value(calc.molecule.get_all_coords());
    mp2.solve_residual_equations(0,0);

    if(world.rank() == 0) printf("\nfinished at time %.1fs\n\n", wall_time());
    world.gop.fence();
    finalize();

    return 0;
}
