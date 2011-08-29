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


// according to McQuarrie
static double he_orbital_McQuarrie(const coord_3d& r) {

    // separation for 2-way decomposition (SVD; r1 -- r2)
    const double x1=r[0];
    const double y1=r[1];
    const double z1=r[2];

    const double r1 = sqrt(x1*x1 + y1*y1 + z1*z1 + 0.001*0.001);

    const double val=exp(-(27.0/16.0)*r1);

    return val;
}


namespace madness {


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
        {}

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
            calc.set_protocol(world,1e-4);
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
                calc.set_protocol(world,1e-6);
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
            pairfunctionT function;         ///< pair function for a specific pair
            double first_order_correction;  ///< this plus orbital energies will yield the HF energy
            double second_order_correction; ///< this plus the HF energy will yield the MP2 correlation energy
            double second_order_energy;     ///< pair energy
            bool solved;                    ///< has the residual equation been solved for this pair?
        };

        World& world;                           ///< the world
        HartreeFock hf;                         ///< our reference

        std::vector<ElectronPair> pairs;        ///< pair functions and energies
        bool solved;                            ///< flag if the residual equations are already solved

        static const double dcut=1.e-6;

    public:
        MP2(World& world, const HartreeFock& hf)
            : world(world)
            , hf(hf)
            , solved(false) {

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
            print("MP2::residual_equations_solved ignores the set of coordinates");
            return solved;
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
            return pairs[ij].function;
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
            functionT J=hf.get_coulomb_potential().scale(-0.5);

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

            printf("first order energy contribution and total energy of pair (%2d %2d)  : %12.8f %12.8f \n",i,j,fo_energy, energy);
            return fo_energy;
        }

        /// given 0th and 1st order pair function, compute the pair energy
        double compute_second_order_correction(const int i, const int j, const pairfunctionT& fo_function) const {

            // the second order energy is given by the expression
            //  E^2 = <phi^0 | V^1 | phi^1>
            // with V^1=J-K+1/r12

            // | phi^0 >
            const pairfunctionT zo_function=zeroth_order_function(i,j);
//            const pairfunctionT& fo_function=pairs[make_ij(i,j)].function;
            MADNESS_ASSERT(fo_function.get_impl());

            // J and K
            MADNESS_ASSERT(is_helium);  // scale 0.5*J, leaving out K
            functionT J=hf.get_coulomb_potential().scale(-0.5);

            real_function_6d eri=ERIFactory<double,6>(world).dcut(1.e-8);
            real_function_6d v11=CompositeFactory<double,6,3>(world)
                        .ket(copy(zo_function).get_impl())
                        .g12(eri.get_impl())
                        .V_for_particle1(copy(J).get_impl())
                        .V_for_particle2(copy(J).get_impl())
                        ;

            const double so_energy=inner(fo_function,v11);

            printf("second order energy contribution and total energy of pair (%2d %2d) : %12.8f \n",i,j,so_energy);
            return so_energy;
        }

        /// given 0th and 1st order pair function, compute the pair energy using the Hylleraas functional
        double compute_second_order_correction_with_Hylleraas(const int i, const int j, const pairfunctionT& fo_function) const {

            // the Hylleraas functional is given by
            //  E^2 = -2 * <phi^1 | V^1 | phi^0>  + <phi^1 | H^0 | phi^1>

            // the B term
            double B=0.0;
            {
                const pairfunctionT zo_function=zeroth_order_function(i,j);

                // V_nuc, J, and K
                MADNESS_ASSERT(is_helium);  // scale 0.5*J, leaving out K
                functionT coulomb=hf.get_coulomb_potential().scale(-0.5);
                functionT v_nuc=hf.get_nuclear_potential();
                functionT v_total=v_nuc+coulomb;

                real_function_6d eri=ERIFactory<double,6>(world).dcut(1.e-8);
                real_function_6d v11=CompositeFactory<double,6,3>(world)
                                         .ket(copy(zo_function).get_impl())
                                         .g12(eri.get_impl())
                                         .V_for_particle1(copy(v_total).get_impl())
                                         .V_for_particle2(copy(v_total).get_impl())
                                         ;

                const double pe=inner(zo_function,v11);

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

                B=ke+pe;
            }

            // the V term
            const double V=compute_second_order_correction(i,j,fo_function);

            printf("Hylleraas functional (energy, and terms V and B) of pair (%2d %2d)  : %12.8f %12.8f %12.8f \n",i,j,2*V+B,V,B);
            return 2*V+B;
        }


        /// solve the residual equations for orbitals i and j
        ElectronPair solve_residual_equation(const int i, const int j) const {

            // the result
            ElectronPair result;
            result.first_order_correction=compute_first_order_correction(i,j);

            pairfunctionT pairfunction;    // eventually the final pair function
            pairfunctionT constant_term;   // the constant term in the residual equation that doesn't change during the iterations
            const pairfunctionT zo_function=zeroth_order_function(i,j);

            // the Green's function depends on the zeroth order energy, which is the sum
            // of the orbital energies of orbitals i and j
            //  -2.0 G = (T - e_i - e_j) ^ -1
            const double eps=zeroth_order_energy(i,j);
            real_convolution_6d green = BSHOperator<6>(world, sqrt(-2*eps), 0.00001, 1e-6);

            // for estimating the MRA structure of V*phi
            real_convolution_6d op_mod = BSHOperator<6>(world, sqrt(-2*eps), 0.00001, 1e-6);
            op_mod.modified=true;


            // apply the first order Hamiltonian on the zeroth order wave function: J, K, and 1/r12 part
            {
                // two-electron interaction potential
                pairfunctionT eri=ERIFactory<double,6>(world).dcut(dcut);

                // Coulomb and exchange potential
                MADNESS_ASSERT(is_helium);
                functionT coulomb=hf.get_coulomb_potential();
                coulomb.scale(-0.5);

                pairfunctionT vphi=CompositeFactory<double,6,3>(world)
                                                    .ket(copy(zo_function).get_impl())
                                                    .g12(eri.get_impl())
                                                    .V_for_particle1(copy(coulomb).get_impl())
                                                    .V_for_particle2(copy(coulomb).get_impl())
                                                    ;
                // make the tree
                vphi.get_impl()->convolute(op_mod);
                vphi.scale(-2.0);
                vphi.print_size("vphi of the zeroth order pair function with 1-electron potentials");

                /// apply the convolution
                constant_term=green(vphi).truncate();
                constant_term.print_size("result of applying 1st order Hamiltonian on 0th order wave function");


            }

            // scale the zeroth order wave function with the first order energy
            {
                pairfunctionT phi=zo_function;
                phi.scale(-2.0);
                pairfunctionT tmp=green(phi).truncate();
                constant_term+=tmp.scale(-result.first_order_correction);
            }


            // start iterations on the first order wave function
            result.function=constant_term;

            // orthogonalize the first order wave function against the reference
            double ovlp=inner(result.function,zo_function);
            printf("overlap of 1st order and 0th order wave functions: %12.8f\n",ovlp);
            pairfunctionT tmp2=-ovlp*zo_function;
            result.function+=tmp2;
            ovlp=inner(result.function,zo_function);
            printf("overlap of 1st order and 0th order wave functions after orthogonalization: %12.8f\n",ovlp);


            for (int ii=0; ii<10; ++ii) {

                // nuclear, Coulomb and exchange potential
                functionT coulomb=hf.get_coulomb_potential().scale(0.5);
                functionT v_total=hf.get_nuclear_potential()+coulomb;

                pairfunctionT vphi=CompositeFactory<double,6,3>(world)
                                     .ket(copy(result.function).get_impl())
                                     .V_for_particle1(copy(v_total).get_impl())
                                     .V_for_particle2(copy(v_total).get_impl());

                // make the tree
                vphi.get_impl()->convolute(op_mod);
                vphi.scale(-2.0);
                vphi.print_size("vphi of the first order pair function with 1-electron potentials");

                /// apply the convolution
                pairfunctionT tmp=green(vphi).truncate();
                tmp.print_size("result of applying 0th order Hamiltonian on 1st order wave function");

                result.function=constant_term+tmp;


                result.second_order_correction=compute_second_order_correction(i,j,result.function);
//                compute_second_order_correction_with_Hylleraas(i,j,result.function);
                if (world.rank()==0) printf("finished iteration %2d at time %.1fs\n\n", ii, wall_time());

            }

            result.solved=true;
            return result;

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
    };
}

int main(int argc, char** argv) {
    initialize(argc, argv);
    World world(MPI::COMM_WORLD);
    startup(world,argc,argv);
    std::cout.precision(6);

    // defaults
    double L = 16;   // box size
    long k = 4 ;        // wavelet order
    double thresh = 1.e-2; // precision
    TensorType tt=TT_2D;

    // set the parameters
    bool restart=false;
    std::string restart_name;

    for(int i = 1; i < argc; i++) {
        const std::string arg=argv[i];

        // break parameters into key and val
        size_t pos=arg.find("=");
        std::string key=arg.substr(0,pos);
        std::string val=arg.substr(pos+1);

        if (key=="size") L=atof(val.c_str());               // usage: size=10
        if (key=="k") k=atoi(val.c_str());                  // usage: k=5
        if (key=="thresh") thresh=atof(val.c_str());        // usage: thresh=1.e-3
        if (key=="TT") {
            if (val=="TT_2D") tt=TT_2D;
            else if (val=="TT_3D") tt=TT_3D;
            else if (val=="TT_FULL") tt=TT_FULL;
            else {
                print("arg",arg, "key",key,"val",val);
                MADNESS_EXCEPTION("confused tensor type",0);
            }

        }
        if (key=="restart") {                               // usage: restart=path/to/mo_file
            restart_name=stringify(val);
            restart=true;
        }
    }



    if (world.rank()==0) {
        print("restart mode",restart," restart_name=", restart_name);
        printf("\nstarting at time %.1fs\n\n", wall_time());
    }


    FunctionDefaults<3>::set_k(k);
    FunctionDefaults<3>::set_thresh(thresh);
    FunctionDefaults<3>::set_cubic_cell(-L/2,L/2);
//    FunctionDefaults<3>::set_cubic_cell(-L,L);

    FunctionDefaults<6>::set_k(k);
    FunctionDefaults<6>::set_thresh(thresh);
    FunctionDefaults<6>::set_cubic_cell(-L/2,L/2);
//    FunctionDefaults<6>::set_cubic_cell(-L,L);
    FunctionDefaults<6>::set_tensor_type(tt);
    FunctionDefaults<6>::set_apply_randomize(true);



    if (world.rank()==0) {
        print("polynomial order:  ", FunctionDefaults<6>::get_k());
        print("threshold:         ", FunctionDefaults<6>::get_thresh());
        print("cell size:         ", L);
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

    if (world.rank()==0) {
        print("size consistency of the 6d green's function?");
        print("");
    }

    Calculation calc(world,"input");
    calc.param.print(world);
    HartreeFock hf(world,calc);
    hf.value();
    MP2 mp2(world,hf);
//    mp2.compute_first_order_correction(0,0);
    mp2.solve_residual_equation(0,0);

//    mp2.value(calc.molecule.get_all_coords());


    if(world.rank() == 0) printf("\nfinished at time %.1fs\n\n", wall_time());
    world.gop.fence();
    finalize();

    return 0;
}
