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
static const double dcut=1.e-4;

static double slater(const coord_3d& r) {
    const double x=r[0], y=r[1], z=r[2];
    const double r2= x*x + y*y + z*z;
    return exp(-2.0*sqrt(r2));
}


template<size_t NDIM>
void save_function(World& world, const Function<double,NDIM>& pair, const std::string name) {
    if (world.rank()==0) print("saving function",name);
    archive::ParallelOutputArchive ar(world, name.c_str(), 1);
    ar & pair;
}

template<size_t NDIM>
void load_function(World& world, Function<double,NDIM>& pair, const std::string name) {
    if (world.rank()==0) print("loading function",name);
    archive::ParallelInputArchive ar(world, name.c_str());
    ar & pair;
}


//static double gauss_3d(const coord_3d& r) {
//    const double x=r[0], y=r[1], z=r[2];
//    const double r2= x*x + y*y + z*z;
//    const double norm=0.712705695388313;
//    return norm*exp(-r2);
//}
//
//static double gauss_6d(const coord_6d& r) {
//    coord_3d r1, r2;
//    r1[0]=r[0],    r1[1]=r[1],    r1[2]=r[2];
//    r2[0]=r[3],    r2[1]=r[4],    r2[2]=r[5];
//    return gauss_3d(r1)*gauss_3d(r2);
//}
//
//
//static double r2r(const coord_6d& r) {
//    coord_3d r1, r2;
//    r1[0]=r[0],    r1[1]=r[1],    r1[2]=r[2];
//    r2[0]=r[3],    r2[1]=r[4],    r2[2]=r[5];
//    double g1=gauss_3d(r1);
//    return g1*g1*gauss_3d(r2);
//}
//
//static double one(const coord_6d& r) {
//    return 1.0;
//}
//static double one(const coord_3d& r) {
//    return 1.0;
//}


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

class YetAnotherWrapperClass {
    const FunctionFunctorInterface<double,6>& func;
    const Tensor<double>& qx;

public:
    YetAnotherWrapperClass(const FunctionFunctorInterface<double,6>& func)
        : func(func)
        , qx(FunctionCommonData<double,6>::get(FunctionDefaults<6>::get_k()).quad_x) {
    }

    void operator()(const Key<6>& key, Tensor<double>& t) const {
        Tensor<double> fval=func.values(key, qx).full_tensor(); // func.values returns coeffT in TT_FULL
        t.emul(fval);
    }
};


static double slaterf12(const coord_6d& r) {
    const double rr=r12(r);
    static const double gamma=1.0;
    static const double alpha=1.0;
    const double f= (1.0-exp(-gamma*rr))/(2.0*gamma);
    const double phi=exp(-alpha*sqrt(r[0]*r[0] + r[1]*r[1] + r[2]*r[2]));
    return phi*f;
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
            if (node.is_leaf()) {
                return leaf_value;
            } else {
                return parent_value;
            }
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

            for (int axis=0; axis<3; ++axis) {
                real_derivative_3d D = free_space_derivative<double,3>(world, axis);
                const real_function_3d Di=D(phi_i);
                const real_function_3d Dj=D(phi_j);
                const real_function_6d Dij=(hartree_product(Di,phi_j) - hartree_product(phi_i,Dj)).truncate();
                const real_function_6d u=real_factory_6d(world).functor2(U(gamma,axis)).is_on_demand();
                result=result+u*Dij;
            }

            // include the purely local potential that (partially) cancels 1/r12
            if (gamma>0.0) {
                real_function_6d pair=hartree_product(phi_i,phi_j);
                const double thresh2=pair.thresh()*0.01;
                real_function_6d fg3=real_factory_6d(world).functor2(fg(gamma))
                        .thresh(thresh2).is_on_demand();
                pair.set_thresh(thresh2);
                real_function_6d mul=fg3*pair;
                mul.set_thresh(thresh2*100.0);
                mul.print_size("mul");
                result=result+mul;
            }
            result.truncate();
            result.print_size("U * |ij>");
            return result;
        }


    private:
        /// functor for the local potential (1-f12)/r12 + sth (doubly connected term of the commutator)
        struct fg {
            double gamma;
            fg(double gamma) : gamma(gamma) {MADNESS_ASSERT(gamma>0.0);}
            double operator()(const coord_6d& r) const {
                const double rr=r12(r);
                const double e=exp(-gamma*rr);
                return (1.0-e)*u(rr,dcut) + 0.5*gamma*e;
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


        /// return the number of occupied orbitals
        int nocc() const {
            MADNESS_ASSERT(calc.param.spin_restricted);
            return calc.param.nalpha;
        }
    };




    /// a class for computing the first order wave function and MP2 pair energies
    class MP2 : public OptimizationTargetInterface {

        typedef real_function_6d pairfunctionT;
        typedef std::vector<pairfunctionT> vecpairfuncT;

        struct ElectronPair {
            real_function_6d function;      ///< pair function for a specific pair w/o correlation factor part
            real_function_6d r12_phi;       ///< orbital product multiplied with the correlation factor
            double first_order_correction;  ///< this plus orbital energies will yield the HF energy
            double second_order_correction; ///< this plus the HF energy will yield the MP2 correlation energy
            double second_order_energy;     ///< pair energy
            bool solved;                    ///< has the residual equation been solved for this pair?
        };

        World& world;                           ///< the world
        HartreeFock hf;                         ///< our reference
        CorrelationFactor corrfac;              ///< correlation factor: Slater or linear

        std::vector<ElectronPair> pairs;        ///< pair functions and energies
        bool solved;                            ///< flag if the residual equations are already solved

        real_convolution_3d poisson;

        static const double dcut=1.e-6;

    public:
        MP2(World& world, const HartreeFock& hf, const CorrelationFactor& corrfac)
            : world(world)
            , hf(hf)
            , corrfac(corrfac)
            , solved(false)
            , poisson(CoulombOperator(world,0.0001,hf.get_calc().param.econv)) {

            // number of pairs:
            const int nocc=hf.nocc();
            const int npairs=nocc*(nocc+1)/2;
            pairs.resize(npairs);

        }

        /// return the molecular energy as a function of the coordinates
        double value(const Tensor<double>& x)  {

            // solve the residual equations for all pairs ij
            if (not residual_equations_solved(x)) {
                for (int i=0; i<hf.nocc(); ++i) {
                    for (int j=i; j<hf.nocc(); ++j) {
                        solve_residual_equation(i,j);
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
        void plot_along(const pairfunctionT& pair, const std::string name) const {

            double L=FunctionDefaults<6>::get_cell_width()[0];
            coord_6d lo(0.05), hi(0.05);
            lo[0]=-L/2;
            hi[0]=L/2;
            trajectory<6> line=trajectory<6>::line2(lo,hi,600);
            madness::plot_along<6>(world,line,pair,(name+"lineplot"));
        }

        /// plot the MRA structure
        void plot_plane(const real_function_6d& f, const std::string filename, const std::string plane) const {

            coord_6d fix_coord(0.0);
            // electron 2:
            fix_coord[4]=0.1;

            f.get_impl()->print_plane(filename,"xy",fix_coord);

        }

        /// return the 0th order energy of pair ij (= sum of orbital energies)
        double zeroth_order_energy(const int i, const int j) const {
            return hf.orbital_energy(i)+hf.orbital_energy(j);
        }

        /// return the 1st order energy correction to the HF reference (=sum of orbital energies minus HF energy)
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
            f.truncate();
            double norm=f.norm2();
            f.scale(1.0/norm);

            //            functionT orbital=real_factory_3d(world).f(he_orbital_McQuarrie);
            //            double norm=orbital.norm2();
            //            orbital.scale(1.0/norm);

            //            pairfunctionT f=hartree_product(orbital,orbital);
            return f;
        }

        /// return the first order wave function for pair ij
        real_function_6d first_order_function(const int i, const int j) const {
            int ij=make_ij(i,j);
            // we consider this an essential:
            MADNESS_ASSERT(pairs[ij].solved);
            return pairs[ij].function+pairs[ij].r12_phi;
        }

        /// return the 1st order energy correction to the orbitals (= HF energy)
        double compute_first_order_correction(const int i, const int j) const {

            // the first order energy is given by the expression
            //  E^1 = <phi^0 | V^1 | phi^0>
            // with V^1=J-K+1/r12

            // | phi^0 >
            pairfunctionT zo_function=zeroth_order_function(i,j);

            // J and K
            MADNESS_ASSERT(is_helium);  // scale 0.5*J, leaving out K
            functionT J=-0.5*hf.get_coulomb_potential();

            real_function_6d eri=ERIFactory<double,6>(world).dcut(1.e-8);
            real_function_6d v11=CompositeFactory<double,6,3>(world)
                                 .ket(copy(zo_function).get_impl())
                                 .g12(eri.get_impl())
                                 .V_for_particle1(copy(J).get_impl())
                                 .V_for_particle2(copy(J).get_impl())
                                 ;

            const double fo_energy=inner(zo_function,v11);
            const double zo_energy=zeroth_order_energy(i,j);
            const double energy=fo_energy+zo_energy;

            if (world.rank()==0)
                printf("1st order energy contribution and total energy of pair (%2d %2d)  : %12.8f %12.8f \n",
                        i,j,fo_energy,energy);
            return fo_energy;
        }

        /// given 0th and 1st order pair function, compute the pair energy
        double compute_second_order_correction(const int i, const int j, const pairfunctionT& fo_function) const {

            // the second order energy is given by the expression
            //  E^2 = <phi^0 | V^1 | phi^1>

            const pairfunctionT zo_function=zeroth_order_function(i,j);
            const double overlap=inner(fo_function,zo_function);
            if (world.rank()==0) print("<phi^1 | phi^0>",overlap);

            const real_function_6d g=multiply_with_1st_order_Hamiltonian(zo_function);
            const double so_energy=inner(fo_function,g);
            if (world.rank()==0) printf("second order energy contribution and total energy of pair (%2d %2d) : %12.8f \n",i,j,so_energy);
            return so_energy;
        }

        /// given 0th and 1st order pair function, compute the pair energy using the Hylleraas functional
        double compute_second_order_correction_with_Hylleraas(const int i, const int j,
                const pairfunctionT& fo_function) const {

//            const real_function_6d zo_function=zeroth_order_function(i,j);
//            test_accuracy(zo_function,fo_function,i,j);
//            test_accuracy(fo_function,zo_function,i,j);

            // the Hylleraas functional is given by
            //  E^2 = -2 * <phi^1 | V^1 | phi^0>  + <phi^1 | H^0 | phi^1>

            // the B term
            double B=0.0;
            {
                // V_nuc, J, and K
                MADNESS_ASSERT(is_helium);
                functionT coulomb=hf.get_coulomb_potential();
                functionT v_nuc=hf.get_nuclear_potential();
                functionT v_total=v_nuc+coulomb;


                real_function_6d v11=CompositeFactory<double,6,3>(world)
                                         .ket(copy(fo_function).get_impl())
                                         .V_for_particle1(copy(v_total).get_impl())
                                         .V_for_particle2(copy(v_total).get_impl())
                                         ;

                const double pe=inner(fo_function,v11);
                if (world.rank()==0) printf("pe in Hylleraas  %12.8f\n\n" ,pe);

                // exchange energy
                real_function_6d Kphi=apply_exchange(fo_function,hf.orbital(i),1);
                Kphi+=apply_exchange(fo_function,hf.orbital(j),2).truncate();
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
            }

            // the V term
            const double V=compute_second_order_correction(i,j,fo_function);
            const double e=2.0*V+B;
            const double e_opt=-V*V/B;

            printf("V and B  %12.8f %12.8f \n",V,B);
            printf("Hylleraas, direct/optimized of pair (%2d %2d) : %12.8f %12.8f \n\n",i,j,e,e_opt);
            return e_opt;
        }

        /// test the accuracy of <phi^0 | H^0 - E^0 | phi^1> = 0
        void test_accuracy(const real_function_6d& zo_function, const real_function_6d& fo_function,
                const int i, const int j) const {

            // this will get everything but the kinetic energy
            const real_function_6d tmp=this->multiply_with_0th_order_Hamiltonian(fo_function);
            const double pe=inner(tmp,zo_function);
            if (world.rank()==0) printf("pe in test_accuracy  %12.8f\n" ,pe);

            // this will get the kinetic energy
            double ke=0.0;
            for (int axis=0; axis<6; axis++) {
                real_derivative_6d D = free_space_derivative<double,6>(world, axis);
                real_function_6d dfo = D(fo_function);
                real_function_6d dzo = D(zo_function);
                double a=0.5*inner(dfo,dzo);
                ke += a;
                if (world.rank()==0) print("done with axis",axis, a);
            }
            if (world.rank()==0) printf("ke in test_accuracy  %12.8f\n" ,ke);

            // this will get the 0th order energy E^0
            const double ovlp=inner(zo_function,fo_function);
            const double e0=ovlp*zeroth_order_energy(i,j);
            if (world.rank()==0) printf("e-term in test_accuracy  %12.8f\n" ,e0);

            if (world.rank()==0) printf("<phi^1 | phi^0> in test_accuracy  %12.8f\n" ,ovlp);
            if (world.rank()==0) printf("<phi^1 | H^0 - E^0 | phi^0> in test_accuracy  %12.8f\n" ,pe+ke-e0);

        }

        void test2(const int i, const int j) {

            real_function_6d vpair1;
            real_function_6d JKpair;

            const real_function_3d phi_i=hf.orbital(i);
            const real_function_3d phi_j=hf.orbital(j);
            const real_function_3d coulomb=hf.get_coulomb_potential();

            const real_function_6d pair=hartree_product(phi_i,phi_j);
            const double e1=compute_first_order_correction(i,j);


            {
                // make the term (-J + K)|phi^0>

                const real_function_3d JKphi_i=-0.5*coulomb*phi_i;
                const real_function_3d JKphi_j=-0.5*coulomb*phi_j;

                JKpair=(hartree_product(JKphi_i,phi_j) + hartree_product(phi_i,JKphi_j)).truncate();

                // make the term (1/r12)|phi^0> in the MRA structure of JKpair
                const real_function_6d eri=ERIFactory<double,6>(world).dcut(dcut);

                real_function_6d g12pair=CompositeFactory<double,6,3>(world)
                                       .ket(copy(pair).get_impl())
                                       .g12(eri.get_impl());

                g12pair.fill_tree(JKpair);

                // make the term -E^1|phi^0>
                const real_function_6d epair=-e1*pair;

                // finally assemble all terms, include factor -2 for convolution!
//                vpair1=(JKpair+g12pair+epair).truncate();
                vpair1=(JKpair+g12pair).truncate();
//                vpair1=(JKpair).truncate();
                vpair1.print_size("vpair1");
                double norm=vpair1.norm2();
                if (world.rank()==0) print("vpair1.norm2()",norm);
                plot_plane(vpair1,"vpair1","xy");
                plot_along(vpair1,"vpair1");

            }

            // traditional
            real_function_6d vpair3;
            {
                const real_function_6d eri=ERIFactory<double,6>(world).dcut(dcut);

                vpair3=CompositeFactory<double,6,3>(world)
                           .ket(copy(pair).get_impl())
//                           .g12(eri.get_impl())
                           .V_for_particle1(copy(-0.5*coulomb).get_impl())
                           .V_for_particle2(copy(-0.5*coulomb).get_impl());

                // for estimating the MRA structure of V*phi
                const double eps=zeroth_order_energy(i,j);
                real_convolution_6d op_mod = BSHOperator<6>(world, sqrt(-2*eps), 0.00001, 1e-6);
                op_mod.modified()=true;

                // make the tree
                vpair3.get_impl()->convolute(op_mod);
                vpair3.print_size("vpair3 after convolute");

                // scale the zeroth order wave function with the first order energy
//                vpair3-=e1*pair;
                vpair3.truncate();

                vpair3.print_size("vpair3");
                double norm=vpair3.norm2();
                if (world.rank()==0) print("vpair3.norm2()",norm);
                plot_plane(vpair3,"vpair3","xy");
                plot_along(vpair3,"vpair3");

                real_function_6d diff=vpair1-vpair3;
                diff.print_size("diff");
                norm=diff.norm2();
                if (world.rank()==0) print("diff.norm3()",norm);


            }


            // semi-traditional
            {
                const real_function_6d eri=ERIFactory<double,6>(world).dcut(dcut);

                real_function_6d vpair2=CompositeFactory<double,6,3>(world)
                           .ket(copy(pair).get_impl())
                           .g12(eri.get_impl())
                           .V_for_particle1(copy(-0.5*coulomb).get_impl())
                           .V_for_particle2(copy(-0.5*coulomb).get_impl());

                // make the tree
//                vpair2.fill_tree(JKpair);
                vpair2.fill_tree(vpair3);
//                vpair2-=e1*pair;

                vpair2.print_size("vpair2");
                double norm=vpair2.norm2();
                if (world.rank()==0) print("vpair2.norm2()",norm);
                plot_plane(vpair2,"vpair2","xy");
                plot_along(vpair2,"vpair2");

                real_function_6d diff=vpair1-vpair2;
                diff.print_size("diff");
                norm=diff.norm2();
                if (world.rank()==0) print("diff.norm2()",norm);

            }


        }

        void test3(const int i, const int j) {
            real_function_6d k;
            load_function(world,k,"K_commutator");

            // the Green's function
            const double eps=zeroth_order_energy(i,j);
            real_convolution_6d green = BSHOperator<6>(world, sqrt(-2*eps), 0.00001, 1e-6);

            real_function_6d Gk=green(-2.0*k).truncate();
            save_function(world,Gk,"GK_commutator");

        }

        void test(const int i, const int j) {
            const pairfunctionT zo_function=zeroth_order_function(i,j);
            const real_function_6d pair=zo_function;

            double norm;
            // 0th iteration for the residual equations (the constant term)
            ElectronPair result;
            load_function(world,result.function,"term1");
            norm=result.function.norm2();
            if (world.rank()==0) print("norm(term1)       ",norm);
            real_function_6d K_comm;
            load_function(world,K_comm,"GK_commutator");
            norm=K_comm.norm2();
            if (world.rank()==0) print("norm(K_comm)      ",norm);
            result.function=result.function+K_comm;
            norm=result.function.norm2();
            if (world.rank()==0) print("norm(term1+K_comm)",norm);
            load_function(world,result.r12_phi,"r12phi");

            orthogonalize(result.function,pair);
            orthogonalize(result.r12_phi,pair);
            orthogonalize(result.function,result.r12_phi);

            result.first_order_correction=this->compute_first_order_correction(i,j);
            const real_function_6d constant_term=copy(result.function);
            real_function_6d phi1=result.function+result.r12_phi;
            compute_second_order_correction_with_Hylleraas(i,j,phi1);

            // the Green's function depends on the zeroth order energy, which is the sum
            // of the orbital energies of orbitals i and j
            //  -2.0 G = (T - e_i - e_j) ^ -1
            const double eps=zeroth_order_energy(i,j);
            real_convolution_6d green = BSHOperator<6>(world, sqrt(-2*eps), 0.00001, 1e-6);


            for (int ii=0; ii<20; ++ii) {

                real_function_6d vphi=multiply_with_0th_order_Hamiltonian(result.function);

                /// apply the convolution
                vphi.scale(-2.0).truncate();
                const pairfunctionT tmp=green(vphi).truncate();
                tmp.print_size("result of applying 0th order Hamiltonian on 1st order wave function");
                result.function=(constant_term-tmp).truncate();

                orthogonalize(result.function,zo_function);
                orthogonalize(result.function,result.r12_phi);
                phi1=result.function+result.r12_phi;

                compute_second_order_correction_with_Hylleraas(i,j,phi1);
                if (world.rank()==0) printf("finished iteration %2d at time %.1fs\n\n", ii, wall_time());

            }

            result.solved=true;

        }

        /// return the function [K,f12]| phi>, NOT orthogonalized
        real_function_6d apply_K_commutator(const int i, const int j) const {

            MADNESS_ASSERT(is_helium);

            real_function_6d f=r12phi(i,j);
            real_function_6d pair=this->zeroth_order_function(i,j);


            // first the term K f12 |phi>
            real_function_6d kfphi=apply_exchange(f,hf.orbital(i),1);
            kfphi=kfphi+apply_exchange(f,hf.orbital(j),2);

            // now the term f12 K |phi>
            real_function_6d r2=apply_exchange(pair,hf.orbital(i),1);
            r2=r2+apply_exchange(pair,hf.orbital(j),2);

            const double t2=pair.thresh()*0.01;
            real_function_6d f12=real_factory_6d(world).functor2(corrfac).is_on_demand().thresh(t2);
            r2.set_thresh(t2);
            real_function_6d fkphi=(f12*r2);
            fkphi.set_thresh(pair.thresh());
            fkphi.truncate();

            real_function_6d diff=(kfphi-fkphi).truncate();

            diff.print_size("[K1,f12] | phi>");
            double norm=diff.norm2();
            if (world.rank()==0) print("diff.norm2",norm);
            plot_along(diff,"commutator_phi");
            save_function(world,diff,"K_commutator");
//            double ovlp=inner(diff,pair);
//            if (world.rank()==0) print("<phi | [K,f12] | phi>",ovlp);
//            this->orthogonalize(diff,pair);
//            norm=diff.norm2();
//            if (world.rank()==0) print("orthogonalized diff.norm2",norm);
            return diff;
        }


        /// do some load-balancing
        void load_balance(const real_function_6d& f, const bool leaf) const {

            LoadBalanceDeux<6> lb(world);
            if (leaf) lb.add_tree(f,LBCost(1.0,0.1));
            else lb.add_tree(f,LBCost(0.1,1.0));
            FunctionDefaults<6>::redistribute(world, lb.load_balance(2.0,false));
            if(world.rank() == 0) printf("redistributed at time   %.1fs\n", wall_time());

        }

        /// solve the residual equations for orbitals i and j
        ElectronPair solve_residual_equation(const int i, const int j) const {

            const pairfunctionT zo_function=zeroth_order_function(i,j);

            // 0th iteration for the residual equations (the constant term)
            ElectronPair result=guess_mp1(i,j);
            const real_function_6d constant_term=copy(result.function);
            compute_second_order_correction_with_Hylleraas(i,j,constant_term+result.r12_phi);

            // the Green's function depends on the zeroth order energy, which is the sum
            // of the orbital energies of orbitals i and j
            //  -2.0 G = (T - e_i - e_j) ^ -1
            const double eps=zeroth_order_energy(i,j);
            real_convolution_6d green = BSHOperator<6>(world, sqrt(-2*eps), 0.00001, 1e-6);


            for (int ii=0; ii<20; ++ii) {

                real_function_6d vphi=multiply_with_0th_order_Hamiltonian(result.function);

                /// apply the convolution
                vphi.scale(-2.0).truncate();
                const pairfunctionT tmp=green(vphi).truncate();
                tmp.print_size("result of applying 0th order Hamiltonian on 1st order wave function");
                result.function=(constant_term+tmp).truncate();

                orthogonalize(result.function,zo_function);
                orthogonalize(result.function,result.r12_phi);

                compute_second_order_correction_with_Hylleraas(i,j,result.function+result.r12_phi);
                if (world.rank()==0) printf("finished iteration %2d at time %.1fs\n\n", ii, wall_time());

            }

            result.solved=true;
            return result;

        }

        /// test the numerics of Kutzelnigg's regularized potential
        /// <phi | U | phi> = <phi | 1/r12 | phi>
        void test_U() const {

            real_function_3d phi=real_factory_3d(world).f(slater);
            real_function_6d uphi=apply_U(phi,phi);
            plot_along(uphi,"uphi");

//            real_function_6d phi=this->zeroth_order_function(0,0);
//            real_function_6d uphi=apply_U(hf.orbital(0),hf.orbital(0));
//            const double b=inner(phi,uphi);
//            if (world.rank()==0) printf("<phi| U | phi>     %12.8f\n",b);
//
//            real_function_6d eri=ERIFactory<double,6>(world).dcut(dcut);
//            real_function_6d vphi=CompositeFactory<double,6,3>(world)
//                                             .ket(copy(phi).get_impl())
//                                             .g12(eri.get_impl());
//            const double a=inner(phi,vphi);
//
//            if (world.rank()==0) printf("<phi| 1/r12 | phi> %12.8f\n",a);
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

        /// orthogonalize f against function g: |f> <= |f> - |g><g|f>
        void orthogonalize(real_function_6d& f, const real_function_6d& g) const {

            const double gnorm=g.norm2();
            const double fac=1.0/(gnorm*gnorm);
            Projector<double,6> P(g);   // projector on g w/o any normalization
            const real_function_6d h=(f-fac*P(f)).truncate();
            const double ovlp=inner(h,g);
            const double thresh=f.thresh();
            if (ovlp>thresh and world.rank()==0) printf("ovlp in orthogonalize %12.8f\n",ovlp);
            f=h;
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
            x=multiply(x,orbital,particle).truncate();

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
            vphi.get_impl()->convolute(op_mod);
            vphi.print_size("(V_nuc + J1 + J2) |ket>:  made V tree");

            // add exchange
            vphi-=apply_exchange(f,hf.orbital(i),1);
            vphi-=apply_exchange(f,hf.orbital(j),2);
            vphi.print_size("(V_nuc + J - K) |ket>: the tree");

            return vphi;
        }

        /// multiply the given function with the 1st order Hamiltonian, exluding the 1st order energy

        /// @param[in]  f   the function we apply H^1 on
        /// @return     the function g=H^1 f, which is NOT orthogonalized against f
        real_function_6d multiply_with_1st_order_Hamiltonian(const real_function_6d& f) const {

            // for estimating the MRA structure of V*phi
            MADNESS_ASSERT(is_helium);
            const int i=0;
            const int j=0;

            const double eps=zeroth_order_energy(i,j);
            real_convolution_6d op_mod = BSHOperator<6>(world, sqrt(-2*eps), 0.00001, 1e-6);
            op_mod.modified()=true;

            // two-electron interaction potential
            pairfunctionT eri=ERIFactory<double,6>(world).dcut(dcut);

            functionT coulomb=-1.0*hf.get_coulomb_potential();

            real_function_6d vphi=CompositeFactory<double,6,3>(world)
                                 .ket(copy(f).get_impl())
                                 .g12(eri.get_impl())
                                 .V_for_particle1(copy(coulomb).get_impl())
                                 .V_for_particle2(copy(coulomb).get_impl());

            // make the tree
            vphi.get_impl()->convolute(op_mod);
            vphi.print_size("(1/r12 - J1 - J2) |ket>:  made V tree");

            // add exchange
            vphi+=apply_exchange(f,hf.orbital(i),1);
            vphi+=apply_exchange(f,hf.orbital(j),2);
            vphi.print_size("(1/r12 - J + K) |ket>: the tree");

            return vphi;
        }

        /// return the initial guess for the 1st order wave function

        /// will return  |phi> = -2.0 * G (-J + K + 1/r12 -E^1) | ij>, in intermediate normalization
        /// @param[in]  i   orbital i
        /// @param[in]  j   orbital j
        /// @return  result the pair struct
        ElectronPair guess_mp1(const int i, const int j) const {

            MADNESS_ASSERT(is_helium);
            ElectronPair result;

            // make the term (-J + K)|phi^0>
            const real_function_3d phi_i=hf.orbital(i);
            const real_function_3d phi_j=hf.orbital(j);
            const real_function_3d coulomb=hf.get_coulomb_potential();

            const real_function_3d JKphi_i=-0.5*coulomb*phi_i;
            const real_function_3d JKphi_j=-0.5*coulomb*phi_j;
            const real_function_6d JKpair=(hartree_product(JKphi_i,phi_j) + hartree_product(phi_i,JKphi_j)).truncate();

            // apply Kutzelnigg's regularized potential U_12
            const real_function_6d Upair=apply_U(phi_i,phi_j);
            plot_along(Upair,"upair");

            // make the term -E^1|phi^0>
            result.first_order_correction=compute_first_order_correction(i,j);
            real_function_6d pair=hartree_product(phi_i,phi_j);
            const double e1=result.first_order_correction;
            const real_function_6d epair=-e1*pair;

            // get the commutator [K,f12]
            real_function_6d K_comm=this->apply_K_commutator(i,j);

            // finally assemble all terms, include factor -2 for convolution!
            real_function_6d vpair=-2.0*(JKpair+Upair+epair+K_comm);
            vpair.truncate();
            vpair.reduce_rank();
            plot_along(vpair,"vpair");

            // the Green's function
            const double eps=zeroth_order_energy(i,j);
            real_convolution_6d green = BSHOperator<6>(world, sqrt(-2*eps), 0.00001, 1e-6);

            // the conventional part of the 1st order function
            result.function=green(vpair).truncate();
            save_function(world,result.function,"term1");

//            // the r12-derived part of the 1st order function
//            real_function_6d r=r12_term(i,j);
//            save_function(world,r,"r12term");
//            result.function=result.function+r;
            if (world.rank()==0) print("\nleaving out the r12-derived part of the wave function guess\n");

            // and the actual term r12*|phi^0>
            result.r12_phi=r12phi(i,j);
            save_function(world,result.r12_phi,"r12phi");

            orthogonalize(result.function,pair);
            orthogonalize(result.r12_phi,pair);
            orthogonalize(result.function,result.r12_phi);

            plot_along(result.r12_phi,"r12_phi_ortho");
            plot_along(result.function,"result.function_ortho");

            double ovlp=inner(result.function,pair);
            if (world.rank()==0) printf("<psi^1 | phi^0> %12.8f\n",ovlp);
            ovlp=inner(result.r12_phi,pair);
            if (world.rank()==0) printf("<r12 | phi^0>   %12.8f\n",ovlp);
            ovlp=inner(result.function,result.r12_phi);
            if (world.rank()==0) printf("<psi^1 | r12>   %12.8f\n",ovlp);

            return result;
        }

        /// return the r12-term of the mp1 guess
        real_function_6d r12_term(const int i, const int j) const {

            const real_function_3d phi_i=hf.orbital(i);
            const real_function_3d phi_j=hf.orbital(j);
            const real_function_3d coulomb=hf.get_coulomb_potential();

            const real_function_6d pair=hartree_product(phi_i,phi_j);

            // the overlap <O_12 | 1/2 r12 | phi^0>
            MADNESS_ASSERT(is_helium);
            real_function_6d fr12=real_factory_6d(world).functor2(corrfac).is_on_demand();   // includes factor 1/2
            real_function_6d r12pair=fr12*pair;
            const double ovlp=inner(pair,r12pair);

            MADNESS_ASSERT(is_helium);
            const real_function_3d v_nuc=hf.get_nuclear_potential();
            const real_function_3d v_total=v_nuc+0.5*coulomb;

            // make the source tree
            real_function_6d vpair=CompositeFactory<double,6,3>(world)
                                 .ket(copy(pair).get_impl())
                                 .V_for_particle1(copy(v_total).get_impl())
                                 .V_for_particle2(copy(v_total).get_impl());


            const double eps=zeroth_order_energy(i,j);
            real_convolution_6d op_mod = BSHOperator<6>(world, sqrt(-2*eps), 0.00001, 1e-6);
            op_mod.modified()=true;
            vpair.get_impl()->convolute(op_mod);

            // the Green's function
            real_convolution_6d green = BSHOperator<6>(world, sqrt(-2*eps), 0.00001, 1e-6);
            const real_function_6d result=-ovlp*green(-2.0*vpair).truncate();

            if (world.rank()==0) print("\nno strong orthogonality! don't forget to change the sign\n");

            return result;
        }

        /// return f12 * |phi>
        real_function_6d r12phi(const int i, const int j) const {
            real_function_6d pair=this->zeroth_order_function(i,j);
            const double t2=pair.thresh()*0.01;
            real_function_6d fr12=real_factory_6d(world).functor2(corrfac).is_on_demand().thresh(t2);
            pair.set_thresh(t2);
            real_function_6d result=fr12*pair;
            result.set_thresh(t2*100.0);
            result.truncate();
            return result;
        }

        /// apply the regularized 1/r12 potential U(1,2) of Kutzelnigg on a pair of orbitals

        /// use a linear correlation factor: [T,0.5*r12] = U(1,2) + 0.5*r12 T
        /// @param[in]  phi_i   orbital i
        /// @param[in]  phi_j   orbital j
        /// @return result  the 6D-Function U(1,2) |ij>
        real_function_6d apply_U(const real_function_3d phi_i, const real_function_3d& phi_j) const {

            real_function_6d result=corrfac.apply_U(phi_i,phi_j);
            return result;
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

    // get command line parameters (overrides input file)
    for(int i = 1; i < argc; i++) {
        const std::string arg=argv[i];

        // break parameters into key and val
        size_t pos=arg.find("=");
        std::string key=arg.substr(0,pos);
        std::string val=arg.substr(pos+1);

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
    calc.set_protocol<3>(world,calc.param.econv*0.01);
    calc.set_protocol<6>(world,calc.param.econv);
    calc.molecule.set_eprec(std::min(calc.param.econv,1.e-6));
    FunctionDefaults<6>::set_tensor_type(tt);


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

    MP2 mp2(world,hf,f12);
//    mp2.compute_first_order_correction(0,0);
//    mp2.test(0,0);
//    mp2.apply_K_commutator(0,0);
//    mp2.test2(0,0);
//    mp2.test_U();

//    mp2.value(calc.molecule.get_all_coords());
    mp2.solve_residual_equation(0,0);

    if(world.rank() == 0) printf("\nfinished at time %.1fs\n\n", wall_time());
    world.gop.fence();
    finalize();

    return 0;
}
