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

#include<madness/chem/mp2.h>
#include <madness/mra/mra.h>
#include <madness/mra/lbdeux.h>
#include<madness/chem/SCF.h>
#include <madness/mra/nonlinsol.h>
#include<madness/chem/correlationfactor.h>
#include<madness/chem/nemo.h>
#include<madness/chem/localizer.h>
#include<madness/chem/ccpairfunction.h>
#include<madness/chem/CCStructures.h>

#include <iostream>

static const double dcut = 1.e-7;
static const double lo = 1.e-6;
static const double bsh_eps = 1.e-7;

using namespace madness;

namespace madness {

/// do some load-balancing

/// @param[in]	f		the function we want to distribute evenly
/// @param[in]	leaf	if true: weigh leaf nodes only; if false: weigh internal nodes only
void load_balance(const real_function_6d& f, const bool leaf) {
    LoadBalanceDeux<6> lb(f.world());
    if (leaf)
        lb.add_tree(f, LBCost(1.0, 0.1));
    else
        lb.add_tree(f, LBCost(0.001, 1.0));
    FunctionDefaults<6>::redistribute(f.world(), lb.load_balance(2.0, false));
}

/// ctor
MP2::MP2(World& world, const commandlineparser& parser) : world(world),
                   param(world, parser), corrfac(world), correlation_energy(0.0),
                   coords_sum(-1.0), Q12(world), ttt(0.0), sss(0.0) {
    {
        std::shared_ptr<Nemo> nemo = std::shared_ptr<Nemo>(new Nemo(world, parser));

        // by default SCF sets the truncate_mode to 1
        FunctionDefaults<3>::set_truncate_mode(3);
        FunctionDefaults<6>::set_truncate_mode(3);

        // get parameters form input file for hf
        if (world.rank() == 0)
            print("accuracy from dft will be overriden by mp2 to 0.01*thresh");
        nemo->get_calc()->set_protocol<6>(world, param.thresh());
        nemo->param.set_derived_value("econv", param.thresh() * 0.01);
        nemo->get_calc()->set_protocol<3>(world, nemo->param.econv());

        if (world.rank() == 0) nemo->get_calc()->molecule.print();

        hf = std::shared_ptr<HartreeFock>(new HartreeFock(world, nemo));
        poisson = std::shared_ptr<real_convolution_3d>(
                CoulombOperatorPtr(world, 0.0001, nemo->param.econv()));

        // construct electronic correlation factor only, nuclear correlation
        // factor depends on the coordinates and must be reassigned for
        // for each geometric structure
        corrfac = CorrelationFactor(world, 1.0, dcut, nemo->get_calc()->molecule);

    }


    // print some output for the user
    if (world.rank() == 0) {
        hf->nemo_ptr->param.print("reference");
        param.print("mp2", "mp2_end");
    }

    // check the user input for consistency
    param.check_input(hf);

    // initialize the electron pairs and store restart info
    // or load previous restart information
    for (int i = param.freeze(); i < hf->nocc(); ++i) {
        for (int j = i; j < hf->nocc(); ++j) {
            pairs.insert(i, j, ElectronPair(i, j));
            if (param.restart())
                pairs(i, j).load_pair(world);
        }
    }


}

/// return the molecular correlation energy energy (without the HF energy)
double MP2::value() {
    hf->value();        // make sure the reference is converged
    if (not check_core_valence_separation()) enforce_core_valence_separation();
    return value(hf->get_calc().molecule.get_all_coords());
}

bool MP2::check_core_valence_separation() const {

    MolecularOrbitals<double, 3> mos(hf->nemos(), hf->get_calc().aeps, {}, hf->get_calc().aocc, {});
    mos.recompute_localize_sets();

    // check that freeze is consistent with the core/valence block-diagonal structure of the fock matrix
    bool fine = false;
    int ntotal = 0;
    std::vector<Slice> slices = MolecularOrbitals<double, 3>::convert_set_to_slice(mos.get_localize_sets());
    for (std::size_t i = 0; i < slices.size(); ++i) {
        const Slice& s = slices[i];
        int n_in_set = s.end - s.start + 1;
        ntotal += n_in_set;
        if (param.freeze() == ntotal) {
            fine = true;
            break;
        }
    }
    if (not fine) {
        print("inconsistent core-valence separation and number of frozen orbitals");
        print("# frozen orbitals", param.freeze());
        mos.pretty_print("initial orbitals");
//        MADNESS_EXCEPTION("inconsistency ", 1);
    }

    // fast return if possible
    Tensor<double> fock1 = get_fock_matrix();
    return Localizer::check_core_valence_separation(fock1, mos.get_localize_sets(),true);
}


void MP2::enforce_core_valence_separation() {

    MolecularOrbitals<double, 3> mos(hf->nemos(), hf->get_calc().aeps, {}, hf->get_calc().aocc, {});
    mos.recompute_localize_sets();

    Tensor<double> fock1 = get_fock_matrix();

    // localize the orbitals respecting core-valence separation
    Localizer localizer(world,hf->get_calc().aobasis,hf->get_calc().molecule,hf->get_calc().ao);
    localizer.set_enforce_core_valence_separation(true).set_method(hf->nemo_ptr->param.localize_method());
    localizer.set_metric(hf->nemo_ptr->R);

    const auto lmo=localizer.localize(mos,fock1,true);
    hf->reset_orbitals(lmo);
    fock.clear();
    print("localized fock matrix");
    Tensor<double> fock2=get_fock_matrix();
    MADNESS_CHECK(Localizer::check_core_valence_separation(fock2,lmo.get_localize_sets()));
    lmo.pretty_print("localized MOs");

};


/// return the molecular correlation energy as a function of the coordinates
double MP2::value(const Tensor<double>& x) {

    // fast return if the MP2 energy is already solved at this geometry
    double xsq = x.sumsq();
    if (xsq == coord_chksum())
        return correlation_energy;

    // nuclear correlation factor depends on the coordinates
    nuclear_corrfac = hf->nemo_ptr->ncf;

    // make sure HF used the same geometry as we do
    coords_sum = xsq;
    MADNESS_CHECK(std::fabs(coord_chksum() - hf->coord_chksum()) < 1.e-14);

    // set all orbitals spaces
    // When a nuclear correlation factor is used the residual equations
    // are similarity transformed. Therefore the orbitals in the
    // projection operator must be set accordingly.
    if (nuclear_corrfac->type() == NuclearCorrelationFactor::None) {
        Q12.set_spaces(hf->get_calc().amo);
    } else {
        // only valid for closed shell
        MADNESS_CHECK(hf->get_calc().param.spin_restricted());
        const std::vector<real_function_3d>& nemos = hf->nemos();
        const std::vector<real_function_3d>& R2amo = hf->R2orbitals();
        Q12.set_spaces(R2amo, nemos, R2amo, nemos);
        if (world.rank() == 0) {
            print("set orbital spaces for the SO projector");
            print("Q12,R = (1-|nemo><nemo|R2) (1-|nemo><nemo|R2)");
        }
    }

    correlation_energy = 0.0;
    if (world.rank() == 0) print("localize ", hf->get_calc().param.do_localize());

    // compute only one single pair
    if ((param.i() > -1) and (param.j() > -1)) {
        pairs(param.i(), param.j()) = make_pair(param.i(), param.j());    // initialize
        solve_residual_equations(pairs(param.i(), param.j()), param.econv() * 0.05, param.dconv());    // solve
        correlation_energy += pairs(param.i(), param.j()).e_singlet + pairs(param.i(), param.j()).e_triplet;

        return correlation_energy;
    }

    // compute the 0th order term and do some coarse pre-iterations
    for (int i = param.freeze(); i < hf->nocc(); ++i) {
        for (int j = i; j < hf->nocc(); ++j) {
            pairs(i, j) = make_pair(i, j);                // initialize
            solve_residual_equations(pairs(i, j), param.econv() * 0.5, param.dconv());
            correlation_energy += pairs(i, j).e_singlet + pairs(i, j).e_triplet;
        }
    }
    if (world.rank() == 0) {
        printf("current decoupled mp2 energy %12.8f\n", correlation_energy);
    }
    if (param.no_compute()) return correlation_energy;

    correlation_energy = 0.0;
    if (hf->get_calc().param.do_localize()) {
        // solve the coupled MP1 equations
        correlation_energy = solve_coupled_equations(pairs, param.econv() * 0.1, param.dconv());

    } else {

        // solve the canonical MP1 equations with increased accuracy
        for (int i = param.freeze(); i < hf->nocc(); ++i) {
            for (int j = i; j < hf->nocc(); ++j) {
                pairs(i, j).converged = false;
                solve_residual_equations(pairs(i, j), param.econv() * 0.05, param.dconv());
                correlation_energy += pairs(i, j).e_singlet + pairs(i, j).e_triplet;
            }
        }
    }
    return correlation_energy;
}


/// solve the residual equation for electron pair (i,j)
void MP2::solve_residual_equations(ElectronPair& result,
                                   const double econv, const double dconv) const {


    if (result.converged or param.no_compute()) {
        result.print_energy();
        return;
    }

    const int i = result.i;
    const int j = result.j;

    if (world.rank() == 0) printf("\n\nsolving electron pair (%d, %d)\n\n", i, j);
    guess_mp1_3(result);
    result.store_pair(world);

    double energy = compute_energy(result);
    if (world.rank() == 0)
        printf("finished with prep step at time %6.1fs with energy %12.8f\n\n",
               wall_time(), energy);

    // the Green's function depends on the zeroth order energy, which is the sum
    // of the orbital energies of orbitals i and j
    //  -2.0 G = (T - e_i - e_j) ^ -1
    const double eps = zeroth_order_energy(i, j);
    if (world.rank() == 0) print("eps in green:  ", eps);
    real_convolution_6d green = BSHOperator<6>(world, sqrt(-2 * eps), lo,
                                               bsh_eps);
	green.destructive()=true;
    if (world.rank() == 0)
        std::cout << "Constructed Green Operator is destructive ? :" << green.destructive() << std::endl;

    NonlinearSolverND<6> solver(param.maxsub());
    solver.do_print = (world.rank() == 0);
    // increment iteration counter upon entry
    for (++result.iteration; result.iteration <= param.maxiter(); ++result.iteration) {

        result.function.print_size("psi");

        // apply the convolution
        real_function_6d vphi = multiply_with_0th_order_Hamiltonian(
                result.function, i, j);

        vphi.print_size("Vpsi");

        vphi.scale(-2.0).truncate();
        load_balance(vphi, false);

        real_function_6d tmp = green(vphi);    // green is destructive
        tmp.print_size("GV psi");

        // we have to solve this equation:
        // psi1 = psi0 + GVpsi1 <=> psi0 + GVpsi1 - psi1 = r =0
        real_function_6d tmp1 = (result.constant_term + tmp).truncate();
        tmp1.print_size("const + GVpsi");
        tmp = Q12(tmp1);
        tmp.print_size("Q12(const + GVpsi)");
        //tmp = (Q12(result.constant_term + tmp)).truncate();

        real_function_6d residual = result.function - tmp;
        result.function = Q12(solver.update(tmp, residual));
        const double rnorm = residual.norm2();
        const double fnorm = result.function.norm2();
        if (world.rank() == 0)
            printf("norm2 of psi, residual %12.8f %12.8f\n", fnorm, rnorm);

        double old_energy = energy;
        energy = compute_energy(result);

        result.converged = ((std::abs(old_energy - energy) < econv) and (
                rnorm < dconv));
        result.store_pair(world);

        if (world.rank() == 0)
            printf("finished iteration %2d at time %8.1fs with energy %12.8f\n\n",
                   result.iteration, wall_time(), energy);

        if (result.converged)
            break;

    }

    // save the converged first order pair function separately for easier access
    std::string name = "pair_" + stringify(i) + stringify(j) + "_psi1_converged";
    save_function(result.function, name);

    // print the final pair energies
    result.print_energy();
}

/// solve the coupled MP1 equations for local orbitals

/// @param[in]  pairs   solved decoupled electron pairs
double MP2::solve_coupled_equations(Pairs<ElectronPair>& pairs,
                                    const double econv, const double dconv) const {

    if (world.rank() == 0) printf("\nsolving coupled MP2 equations\n\n");

    double total_energy = 0.0;
    for (int i = param.freeze(); i < hf->nocc(); ++i) {
        for (int j = i; j < hf->nocc(); ++j) {
            total_energy += compute_energy(pairs(i, j));
        }
    }

    // do the macro-iterations of the whole set of pair functions
    for (int iteration = 0; iteration < param.maxiter(); ++iteration) {

        // compute the coupling between the pair functions
        START_TIMER(world);
        Pairs<real_function_6d> coupling;
        add_local_coupling(pairs, coupling);
        END_TIMER(world, "compute coupling");

        // compute the vector function, aka the rhs of the residual equations
        double total_rnorm = 0.0;
        double old_energy = total_energy;
        total_energy = 0.0;

        for (int i = param.freeze(); i < hf->nocc(); ++i) {
            for (int j = i; j < hf->nocc(); ++j) {

                real_function_6d vectorfunction = multiply_with_0th_order_Hamiltonian(pairs(i, j).function, i, j);

                // NOTE THE SIGN
                vectorfunction -= coupling(i, j);

                const double eps = zeroth_order_energy(i, j);
                START_TIMER(world);
                real_convolution_6d green = BSHOperator<6>(world, sqrt(-2 * eps), lo, bsh_eps);
                vectorfunction.scale(-2.0).truncate();
                real_function_6d tmp = green(vectorfunction).truncate();
                END_TIMER(world, "apply BSH |ket>");

                START_TIMER(world);
                tmp = (Q12(pairs(i, j).constant_term + tmp)).truncate();

                real_function_6d residual = pairs(i, j).function - tmp;
                pairs(i, j).function = tmp;

                const double rnorm = residual.norm2();
                const double fnorm = pairs(i, j).function.norm2();
                if (world.rank() == 0)
                    printf("norm2 of psi, residual %2d %2d "
                           "%12.8f %12.8f\n", i, j, fnorm, rnorm);
                double energy = compute_energy(pairs(i, j));

                total_rnorm += rnorm;
                total_energy += energy;
                END_TIMER(world, std::string(" post-BSH " + stringify(i) + stringify(j)).c_str());
            }
        }

        // check convergence
        bool converged = ((std::abs(old_energy - total_energy) < econv)
                          and (total_rnorm < dconv));

        // save the pairs
        for (int i = param.freeze(); i < hf->nocc(); ++i) {
            for (int j = i; j < hf->nocc(); ++j) {
                pairs(i, j).store_pair(world);
            }
        }

        if (world.rank() == 0)
            printf("finished iteration %2d at time %8.1fs with coupled energy %12.8f\n\n",
                   iteration, wall_time(), total_energy);

        if (converged) break;
    }

    return total_energy;
}

/// compute increments: psi^1 = C + GV C + GVGV C + GVGVGV C + ..

/// for pre-optimization of the mp1 wave function
/// no restart option here for now
/// param[inout]	pair	the electron pair
/// param[in]		green	the Green's function
void MP2::increment(ElectronPair& pair, real_convolution_6d& green) {

    // don't mess up if we've already done some iterations!
    if (pair.iteration > 0)
        return;

    if (world.rank() == 0)
        print("computing increments");
    real_function_6d latest_increment;
    load_function(latest_increment, "GVpair");

    for (int ii = 1; ii < 20; ++ii) {

        real_function_6d vphi = multiply_with_0th_order_Hamiltonian(
                latest_increment, pair.i, pair.j);
        load_balance(vphi, false);

        /// apply the convolution
        vphi.scale(-2.0).truncate();
        latest_increment = green(vphi).truncate();        // green is destructive
        latest_increment.print_size(
                "result of applying 0th order Hamiltonian on latest increment");

        latest_increment = Q12(latest_increment);
        pair.function = pair.function + latest_increment;

        // compute the energy
        compute_energy(pair);

        if (world.rank() == 0)
            printf("finished increment %2d at time %.1fs\n\n", ii, wall_time());
        const double residual_norm = latest_increment.norm2();
        if (residual_norm < pair.function.thresh())
            break;
    }
}

double MP2::asymmetry(const real_function_6d& f, const std::string s) const {
    return 0.0;
    const real_function_6d ff = swap_particles(f);
    double diff = (ff - f).norm2();
    f.check_symmetry();
    if (world.rank() == 0)
        print("asymmetry of ", s, "  ", diff);
    return diff;
}

/// compute the matrix element <ij | g12 Q12 f12 | phi^0>

/// scales quartically. I think I can get this down to cubically by
/// setting Q12 = (1 - O1)(1 - O2) = 1- O1(1 - 0.5 O2) - O2 (1 - 0.5 O1)
/// as for the formulas cf the article mra_molecule
/// @return 	the energy <ij | g Q f | kl>
double MP2::compute_gQf(const int i, const int j, ElectronPair& pair) const {

    CCTimer t1(world,"gQf old");
    // for clarity of notation
    const int k = pair.i;
    const int l = pair.j;

    // the ket space
    const real_function_3d& ket_i = hf->nemo(i);
    const real_function_3d& ket_j = hf->nemo(j);

    // the bra space
    const real_function_3d& bra_k = hf->R2orbital(k);
    const real_function_3d& bra_l = hf->R2orbital(l);

    // compute <ij| fg |kl>: do it in 3D as (ik| fg |jl)
    // the operator fg can be rewritten as 1/r12 - f/r12
    // i.e. as poisson kernel and a bsh kernel. Note the
    // the bsh kernel includes a factor of 1/(4 pi)
    const real_function_3d ik = ket_i * bra_k;
    const real_function_3d jl = ket_j * bra_l;

    // make all the operators that we need
    const double fourpi = 4.0 * constants::pi;
    real_convolution_3d fg = BSHOperator<3>(world, corrfac.gamma(), lo,
                                            bsh_eps / fourpi);
    real_convolution_3d gg = CoulombOperator(world, lo, bsh_eps);
    real_convolution_3d slaterf12 = SlaterF12Operator(world, corrfac.gamma(),
                                                      lo, bsh_eps / fourpi);

    //  < ij | fg | kl >
    const real_function_3d ik_fg = (gg)(ik) - fourpi * fg(ik);
    const double a = inner(ik_fg, jl) / (2.0 * corrfac.gamma());
    if (world.rank() == 0)
        printf("<%d%d | f/r              | %d%d>  %12.8f\n", i, j, k, l, a);

    // compute <ij| f (O1 + O2) g | ij>

    // compute bra space xi(ik,j)^dagger, i.e. the hermitian conjugate of xi
    // the index k is implicit in the vector of functions
    // naming: xi _ orbitals _ operator _ hc
    std::vector<real_function_3d> xi_ij_g_ket = make_xi(ket_i, ket_j, *poisson,
                                                        false);       // xi_{i,m*},j
    std::vector<real_function_3d> xi_ji_g_ket = make_xi(ket_j, ket_i, *poisson,
                                                        false);       // xi_{j,m*},i

    std::vector<real_function_3d> xi_ij_f_bra = make_xi(bra_k, bra_l, slaterf12,
                                                        true);       // xi_{i*,m},j*
    std::vector<real_function_3d> xi_ji_f_bra = make_xi(bra_l, bra_k, slaterf12,
                                                        true);       // xi_{j*,m},i*

    // in the following do NOT use antisymmetrized pair functions:
    // |ij> -> 0.5 |ij - ji>

    // < ij | f12 O1 g12 | kl >
    //   = \sum_m <i(1) j(2) | f12 | m(1) >< m(3) | g23 | k(3) l(2)>
    //   = \sum_m < chi^f_i*,m(2) j*(2) | chi^g_k,m*(2) l(2) >
    //   = \sum_m < xi^f_im,j | xi^g_km,l >
    const double o1a = inner(world, xi_ij_f_bra, xi_ij_g_ket).sum();
    if (world.rank() == 0)
        printf("<%d%d | f12 O1 g12       | %d%d>  %12.8f\n", i, j, k, l, o1a);

    // < ij | f12 O2 g12 | kl >
    //    = \sum_m <i(1) j(2) | f12 | m(2) >< m(3) | g13 | k(1) l(3)>
    //    = \sum_m <chi^f_j*,m(1) i*(1) | chi^g_l,m*(1) k(1) >
    //    = \sum_m < xi^f_jm,i | xi^g_lm,k >
    const double o2a = inner(world, xi_ji_f_bra, xi_ji_g_ket).sum();
    if (world.rank() == 0)
        printf("<%d%d | f12 O2 g12       | %d%d>  %12.8f\n", i, j, k, l, o2a);

    // compute <ij| f O1 O2 g | kl>  // why do I need to swap ij in g_ijkl??
    const Tensor<double> f_ijmn = matrix_inner(world, xi_ij_f_bra, hf->nemos());
    const Tensor<double> g_ijmn = matrix_inner(world, hf->R2orbitals(),
                                               xi_ji_g_ket);
    const double o12 = f_ijmn.trace(g_ijmn);
    if (world.rank() == 0)
        printf("<%d%d | f12 O12 g12      | %d%d>  %12.8f\n", i, j, k, l, o12);

    const double e = a - o1a - o2a + o12;
    if (world.rank() == 0)
        printf("<%d%d | g Q12 f          | %d%d>  %12.8f\n", i, j, k, l, e);

    print("gQf old",t1.reset());

    CCTimer timer(world,"gQf with ccpairfunction");
    CCConvolutionOperator<double,3>::Parameters param;
    auto f=CCConvolutionOperatorPtr<double,3>(world,OpType::OT_F12,param);
    auto g=CCConvolutionOperatorPtr<double,3>(world,OpType::OT_G12,param);
//    print("operator constructor",timer.reset());

    CCPairFunction<double,6> fij(f,ket_i,ket_j);
    CCPairFunction<double,6> gij(g,bra_k,bra_l);
//    CCPairFunction ij({psi},{psi});
    std::vector<CCPairFunction<double,6>> vfij={fij};
    std::vector<CCPairFunction<double,6>> vgij={gij};
//    std::vector<CCPairFunction> vij={ij};
//    print("g/f ij constructor",timer.reset());

//    StrongOrthogonalityProjector<double,3> SO(world);
//    SO.set_spaces({psi},{psi},{psi},{psi});
    std::vector<CCPairFunction<double,6>> Qfij=Q12(vfij);
//    std::vector<CCPairFunction> Qgij=Q12(vgij);
//    print("SO application",timer.reset());

    double result1=inner(vgij,Qfij);
    print("inner",timer.reset());

    print("<ij | g (Q f | ij>)",result1);

    return e;
}

/// save a function
template<typename T, size_t NDIM>
void MP2::save_function(const Function<T, NDIM>& f,
                        const std::string name) const {
    if (world.rank() == 0)
        print("saving function", name);
    f.print_size(name);
    archive::ParallelOutputArchive<archive::BinaryFstreamOutputArchive> ar(world, name.c_str(), 1);
    ar & f;
}

/// load a function
template<typename T, size_t NDIM>
void MP2::load_function(Function<T, NDIM>& f, const std::string name) const {
    if (world.rank() == 0)
        print("loading function", name);
    archive::ParallelInputArchive<archive::BinaryFstreamInputArchive> ar(world, name.c_str());
    ar & f;
    f.print_size(name);
}

/// return the function Uphi0; load from disk if available
real_function_6d MP2::make_Uphi0(ElectronPair& pair) const {
    const int i = pair.i;
    const int j = pair.j;
    real_function_6d Uphi0 = real_factory_6d(world);
    // apply the pure commutator of the kinetic energy and the
    // electronic correlation factor on the (regularized)
    // orbitals i and j
    //  [T,f]
    const double eps = this->zeroth_order_energy(i, j);
    std::cout << "eps=" << eps << std::endl;
    {
        real_convolution_6d op_mod = BSHOperator<6>(world,
                                                    sqrt(-2 * eps), 0.0001, 1e-5);
        op_mod.modified() = true;
        Uphi0 = corrfac.apply_U(hf->nemo(i), hf->nemo(j), op_mod);
    }

    {
        real_function_6d tmp =
                CompositeFactory<double, 6, 3>(world).particle1(
                        copy(hf->R2orbital(i))).particle2(
                        copy(hf->R2orbital(j)));

        const double a = inner(Uphi0, tmp);
        if (world.rank() == 0)
            print("< nemo | R2 U | nemo>", a);
    }
    Uphi0.print_size("Uphi0");
    asymmetry(Uphi0, "Uphi w/o R");

    // apply the mixed commutator of the kinetic energy and the
    // electronic and nuclear correlation factor on the regularized
    // orbitals i and j:
    //  R^{-1} [[T,f],R]
    // Vanishes if no nuclear correlation factor
    // is used, since [[T,f],1] = 0
    //                if (hf->nemo_calc.nuclear_correlation->type()!=NuclearCorrelationFactor::None) {
    if (world.rank() == 0)
        print("adding the mixed commutator R^{-1} [[T,f],R]");
    const double thresh = FunctionDefaults<6>::get_thresh();
    for (int axis = 0; axis < 3; axis++) {

        double tight_thresh = thresh;// std::min(thresh, 1.e-4); // not nec. anymore bc. of fill_cuspy_tree instead of fill_tree

        real_convolution_6d op_mod = BSHOperator<6>(world,
                                                    sqrt(-2 * eps), 0.0001, 1e-5);
        op_mod.modified() = true;

        const real_function_3d u1_nuc =
                hf->nemo_ptr->ncf->U1(axis);
        const real_function_3d u1_nuc_nemo_i = u1_nuc * hf->nemo(i);
        const real_function_3d u1_nuc_nemo_j = u1_nuc * hf->nemo(j);
        const real_function_6d u1_el = corrfac.U1(axis);
        //						u1_nuc.print_size("u1_nuc");
        //						plot_plane(world,u1_nuc,"u1_nuc");
        //						u1_nuc_nemo_i.print_size("u1_nuc_nemo_i");
        //						plot_plane(world,u1_nuc_nemo_i,"u1_nuc_nemo_i");
        //						plot_plane(world,hf->nemo(i),"nemo_i");
        //						u1_nuc_nemo_j.print_size("u1_nuc_nemo_j");
        //						u1_el.print_size("u1_el");

        real_function_6d U1_mix1 =
                CompositeFactory<double, 6, 3>(world).g12(u1_el).particle1(
                        copy(u1_nuc_nemo_i)).particle2(
                        copy(hf->nemo(j))).thresh(tight_thresh);
        U1_mix1.fill_cuspy_tree(op_mod);
        U1_mix1.set_thresh(thresh);
        U1_mix1.print_size("U1_mix1");
        //						plot_plane(world,U1_mix1,"U1_mix1");

        real_function_6d U1_mix2 =
                CompositeFactory<double, 6, 3>(world).g12(u1_el).particle1(
                        copy(hf->nemo(i))).particle2(
                        copy(u1_nuc_nemo_j)).thresh(tight_thresh);
        U1_mix2.fill_cuspy_tree(op_mod);
        U1_mix2.set_thresh(thresh);
        U1_mix2.print_size("U1_mix2");
        //						plot_plane(world,U1_mix2,"U1_mix2");

        real_function_6d diff = (U1_mix1 - U1_mix2).scale(-1.0);
        diff.truncate();

        // multiply with -1 from the kinetic energy operator
        //						Uphi0=(Uphi0+(U1_mix1-U1_mix2).scale(-1.0)).truncate();
        Uphi0 = (Uphi0 + diff).truncate();
        this->asymmetry(Uphi0, "Uphi0 in R");

        //						Uphi0.print_size("Uphi0");
        //						plot_plane(world,Uphi0,"Uphi0");

        {
            real_function_6d tmp =
                    CompositeFactory<double, 6, 3>(world).particle1(
                            copy(hf->R2orbital(i))).particle2(
                            copy(hf->R2orbital(j)));

            const double a = inner(Uphi0, tmp);
            if (world.rank() == 0)
                print("< nemo | R2 U | nemo>", a);
        }
    }
    asymmetry(Uphi0, "Uphi0");


    // sanity check: <ij| [T,g12] |ij> = <ij | U |ij> - <ij| g12 | ij> = 0
    real_function_6d tmp = CompositeFactory<double, 6, 3>(world).particle1(
            copy(hf->R2orbital(i))).particle2(copy(hf->R2orbital(j)));

    const double a = inner(Uphi0, tmp);
    const real_function_3d ii = hf->R2orbital(i) * hf->nemo(i);
    const real_function_3d jj = hf->R2orbital(j) * hf->nemo(j);
    const real_function_3d gii = (*poisson)(ii);
    const double aa = inner(jj, gii);
    const double error = std::fabs(a - aa);
    if (world.rank() == 0) {
        printf("< phi0 | U_R   | phi0 >  %12.8f\n", a);
        printf("< phi0 | 1/r12 | phi0 >  %12.8f\n", aa);
        if (error > FunctionDefaults<6>::get_thresh())
            print("WARNING : Kutzelnigg's potential inaccurate (box size, thresh ?)");
        if (error > FunctionDefaults<6>::get_thresh() * 10.0) MADNESS_EXCEPTION(
                "Kutzelnigg's potential plain wrong (box size, thresh ?)", 1);
    }
    Uphi0.print_size("Uphi0");
    return Uphi0;
}

/// return the function [K,f] phi0; load from disk if available
real_function_6d MP2::make_KffKphi0(const ElectronPair& pair) const {

    const int i = pair.i;
    const int j = pair.j;

    real_function_6d KffKphi0;
    real_function_6d r12nemo = CompositeFactory<double, 6, 3>(world).g12(
            corrfac.f()).particle1(copy(hf->nemo(i))).particle2(
            copy(hf->nemo(j)));
    r12nemo.fill_tree().truncate().reduce_rank();
    r12nemo.print_size("r12nemo");
    //                save_function(r12nemo,"r12nemo");

    real_function_6d Kfphi0;
    Kfphi0 = K(r12nemo, i == j);

    {
        real_function_6d tmp =
                CompositeFactory<double, 6, 3>(world).particle1(
                        copy(hf->R2orbitals()[i])).particle2(
                        copy(hf->R2orbitals()[j]));
        double a1 = inner(Kfphi0, tmp);
        if (world.rank() == 0)
            printf("< nemo0 | R^2 R-1 K f R | nemo0 >  %12.8f\n", a1);
        //					save_function(Kfphi0,"Kfphi0");
    }
    const real_function_6d fKphi0 = make_fKphi0(pair.i, pair.j);
    {
        real_function_6d tmp =
                CompositeFactory<double, 6, 3>(world).particle1(
                        copy(hf->R2orbitals()[i])).particle2(
                        copy(hf->R2orbitals()[j]));
        double a2 = inner(fKphi0, tmp);
        if (world.rank() == 0)
            printf("< nemo0 | R^2 R-1 f K R | nemo0 >  %12.8f\n", a2);
    }
    KffKphi0 = (Kfphi0 - fKphi0).truncate().reduce_rank();

    // sanity check
    real_function_6d tmp = CompositeFactory<double, 6, 3>(world).particle1(
            copy(hf->R2orbitals()[i])).particle2(copy(hf->R2orbitals()[j]));
    const double a = inner(KffKphi0, tmp);
    if (world.rank() == 0) {
        printf("< nemo0 | R^2 R-1 [K,f] R | nemo0 >  %12.8f\n", a);
        if (std::fabs(a) > FunctionDefaults<6>::get_thresh())
            print("WARNING : exchange commutator inaccurate");
        if (std::fabs(a) > FunctionDefaults<6>::get_thresh() * 10.0) MADNESS_EXCEPTION(
                "exchange commutator plain wrong", 1);
    }
    KffKphi0.print_size("KffKphi0");
    return KffKphi0;
}

/// compute some matrix elements that don't change during the SCF
ElectronPair MP2::make_pair(const int i, const int j) const {

    ElectronPair p = ElectronPair(i, j);
    if (param.restart()) p.load_pair(world);

    // compute and store them if they have not been read from disk
    if (p.ij_gQf_ij == ElectronPair::uninitialized()) {
        p.ij_gQf_ij = compute_gQf(i, j, p);
        p.ji_gQf_ij = 0.0;
        if (i != j)
            p.ji_gQf_ij = compute_gQf(j, i, p);
        p.store_pair(world);
    } else {
        if (world.rank() == 0) {
            printf("<%d%d | g Q12 f          | %d%d>  %12.8f\n",
                   i, j, i, j, p.ij_gQf_ij);
            if (i != j) {
                printf("<%d%d | g Q12 f          | %d%d>  %12.8f\n",
                       j, i, i, j, p.ji_gQf_ij);
            }
        }
    }

    if (world.rank() == 0)
        printf("done with matrix elements at time %.1fs\n\n", wall_time());
    return p;
}

/// compute the first iteration of the residual equations and all intermediates
void MP2::guess_mp1_3(ElectronPair& pair) const {

    // fast return if possible
    if (pair.function.is_initialized()) return;

    const int i = pair.i;
    const int j = pair.j;

    // the Green's function
    const double eps = zeroth_order_energy(i, j);
    real_convolution_6d green = BSHOperator<6>(world, sqrt(-2.0 * eps), lo, bsh_eps);

    real_function_6d Uphi0 = make_Uphi0(pair);
    real_function_6d KffKphi0 = real_factory_6d(world);
    if (not param.do_oep()) KffKphi0 = make_KffKphi0(pair);

    // these are the terms that come from the single projectors: (O1 + O2) (U+[K,f])|phi^0>
    std::vector<real_function_3d> phi_k_UK_phi0;
    std::vector<real_function_3d> phi_l_UK_phi0;

    for (int k = 0; k < hf->nocc(); ++k) {
        real_function_3d tmp;
        tmp = Uphi0.project_out(hf->R2orbitals()[k], 0)
              - KffKphi0.project_out(hf->R2orbitals()[k], 0);
        phi_k_UK_phi0.push_back(tmp);

        tmp = Uphi0.project_out(hf->R2orbitals()[k], 1)
              - KffKphi0.project_out(hf->R2orbitals()[k], 1);
        phi_l_UK_phi0.push_back(tmp);
    }


    // make the terms with high ranks and smallish trees
    load_balance(Uphi0, true);
    real_function_6d Vpair1 = (Uphi0 - KffKphi0).truncate().reduce_rank();
    Uphi0.clear();
    KffKphi0.clear();
    asymmetry(Vpair1, "Vpair");
    load_balance(Vpair1, false);
    real_function_6d GVpair;
    green.destructive() = true;            // green will destroy Vpair1
    GVpair = green(-2.0 * Vpair1).truncate().reduce_rank();
    Vpair1.clear(true);
    save_function(GVpair, "GVpair1");
    //			load_function(GVpair,"GVpair1");
    asymmetry(GVpair, "GVpair1");
    pair.function = GVpair;
    compute_energy(pair);

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
    Tensor<double> h(hf->nocc(), hf->nocc());
    const double fourpi = 4.0 * constants::pi;
    real_convolution_3d slaterf12 = SlaterF12Operator(world, corrfac.gamma(),
                                                      lo, bsh_eps / fourpi);

    // compute some intermediates
    const tensorT occ = hf->get_calc().aocc;
    tensorT fock = hf->nemo_ptr->compute_fock_matrix(hf->nemos(), occ);
    if (hf->get_calc().param.do_localize()) {
        // local orbitals -- add coupling
        vecfuncT amotilde = transform(world, hf->orbitals(), fock);
        vecfuncT fik, fjl, gik, jl, ik, tjl, jtl, itk, tik;

        for (int k = 0; k < hf->nocc(); ++k) {
            const real_function_3d ik1 = hf->orbital(i) * hf->orbital(k);
            gik.push_back((*poisson)(ik1).truncate());    // (ik | g12 |
            fik.push_back(slaterf12(ik1).truncate());    // (ik | f12 |
            ik.push_back(ik1);
            itk.push_back(hf->orbital(i) * amotilde[k]);    // | i k~ )
            tik.push_back(amotilde[i] * hf->orbital(k));    // | i~ k )
        }
        for (int l = 0; l < hf->nocc(); ++l) {
            const real_function_3d jl1 = hf->orbital(j) * hf->orbital(l);
            fjl.push_back(slaterf12(jl1).truncate());    // (jl | f12 |
            jl.push_back(jl1);                            //      | jl)
            tjl.push_back(amotilde[j] * hf->orbital(l));    // | j~ l )
            jtl.push_back(hf->orbital(j) * amotilde[l]);    // | j l~ )
        }

        for (int k = 0; k < hf->nocc(); ++k) {
            for (int l = 0; l < hf->nocc(); ++l) {

                const double g_ijkl = inner(jl[l], gik[k]);        // <i j | g | k l>
                const double f_itjkl = inner(fik[k], tjl[l]);    // <i j~| f | k l>
                const double f_ijktl = inner(fik[k], jtl[l]);    // <i j| f | k l~>
                const double f_tijkl = inner(fjl[l], tik[k]);    // <i~ j| f | k l>
                const double f_ijtkl = inner(fjl[l], itk[k]);    // <i j| f | k~ l>

                h(k, l) = g_ijkl - f_itjkl + f_ijktl - f_tijkl + f_ijtkl;
            }
        }

    } else { // canonical orbitals -- simplified code

        for (int k = 0; k < hf->nocc(); ++k) {

            // for the matrix element <ij|kl> = (ik|jl)
            const real_function_3d ik = hf->orbital(i) * hf->orbital(k);
            const real_function_3d gik = (*poisson)(ik).truncate();
            const real_function_3d fik = slaterf12(ik).truncate();

            // the function Of(2) = <k(1) | f12 | phi0 >
            //	            const real_function_3d Of=pair.r12phi.project_out(hf->orbital(k),0);

            for (int l = 0; l < hf->nocc(); ++l) {

                const real_function_3d jl = hf->orbital(j) * hf->orbital(l);
                const double g_ijkl = inner(jl, gik);
                const double f_ijkl = inner(jl, fik);

                // the matrix element <kl | f12 | ij>
                //		            const double f_ijkl=inner(hf->orbital(l),Of);
                const double e_kl = zeroth_order_energy(k, l);
                const double e_ij = zeroth_order_energy(i, j);

                h(k, l) = f_ijkl * (e_kl - e_ij) + g_ijkl;
            }
        }
    }

    // the result function; adding many functions requires tighter threshold
    const double thresh = FunctionDefaults<6>::get_thresh();
    const double tight_thresh = thresh * 0.1;
    FunctionDefaults<6>::set_thresh(tight_thresh);
    real_function_6d r1 = real_factory_6d(world).thresh(tight_thresh);
    for (int k = 0; k < hf->nocc(); ++k) {

        // make the term  tmp2(2) = |UK_k(2)> - 1/2 |l(2)> h(kl)
        real_function_3d tmp2 = phi_k_UK_phi0[k];
        for (int l = 0; l < hf->nocc(); ++l)
            tmp2 -= 0.5 * h(k, l) * hf->nemo(l);

        // now apply the Greens' function (including the -2.0 factor in one term)
        real_function_6d tmp = green(-2.0 * hf->nemo(k), tmp2);
        tmp.set_thresh(tight_thresh);
        r1 = (r1 - tmp).truncate(tight_thresh);
    }
    if (world.rank() == 0)
        print("symmetry r1");
    r1.print_size("r1");

    for (int l = 0; l < hf->nocc(); ++l) {

        // make the term  tmp1(1) = |UK_l(1) - 1/2 |k(1)> h(kl)
        real_function_3d tmp1 = phi_l_UK_phi0[l];
        for (int k = 0; k < hf->nocc(); ++k)
            tmp1 -= 0.5 * h(k, l) * hf->nemo(k);

        // now apply the Greens' function (including the -2.0 factor in one term)
        real_function_6d tmp = green(tmp1, -2.0 * hf->nemo(l));
        tmp.set_thresh(tight_thresh);
        r1 = (r1 - tmp).truncate(tight_thresh);
    }
    r1.truncate();
    r1.print_size("r1");
    GVpair.set_thresh(tight_thresh);
    GVpair = (GVpair + r1).truncate().reduce_rank();

    FunctionDefaults<6>::set_thresh(thresh);
    GVpair.set_thresh(thresh);
    GVpair.truncate().reduce_rank();

    asymmetry(GVpair, "GVpair before Q12");

    GVpair = Q12(GVpair);
    asymmetry(GVpair, "GVpair after Q12");
    pair.function = GVpair;
    pair.constant_term = copy(GVpair);
    save_function(GVpair, "GVpair");
}

/// compute the singlet and triplet energy for a given electron pair

/// @return	the energy of 1 degenerate triplet and 1 singlet pair
double MP2::compute_energy(ElectronPair& pair) const {

    START_TIMER(world);
    // a11 will be the bra space, therefore take the hermitian conjugate
    const double a11 = inner(pair.function,
                             JK1phi0_on_demand(pair.i, pair.j, true))
                       + inner(pair.function, JK2phi0_on_demand(pair.i, pair.j, true));
    if (world.rank() == 0)
        printf("V2: <phi^0 | J-K        | psi^1>  %12.8f\n", a11);

    // this will be the bra space
    real_function_6d eri = TwoElectronFactory(world).dcut(dcut);
    real_function_6d ij_g =
            CompositeFactory<double, 6, 3>(world).particle1(
                    copy(hf->R2orbital(pair.i))).particle2(
                    copy(hf->R2orbital(pair.j))).g12(eri);
    real_function_6d ji_g =
            CompositeFactory<double, 6, 3>(world).particle1(
                    copy(hf->R2orbital(pair.j))).particle2(
                    copy(hf->R2orbital(pair.i))).g12(eri);

    // compute < ij | g12 | psi >
    const double ij_g_uij = inner(pair.function, ij_g);
    if (world.rank() == 0)
        printf("<ij | g12       | psi^1>  %12.8f\n", ij_g_uij);

    // compute < ji | g12 | psi > if (i/=j)
    const double ji_g_uij = (pair.i == pair.j) ? 0 : inner(pair.function, ji_g);
    if (world.rank() == 0)
        printf("<ji | g12       | psi^1>  %12.8f\n", ji_g_uij);

    // the singlet and triplet triplet pair energies
    if (pair.i == pair.j) {
        pair.e_singlet = ij_g_uij + pair.ij_gQf_ij;
        pair.e_triplet = 0.0;
    } else {
        pair.e_singlet = (ij_g_uij + pair.ij_gQf_ij)
                         + (ji_g_uij + pair.ji_gQf_ij);
        pair.e_triplet = 3.0
                         * ((ij_g_uij - ji_g_uij) + (pair.ij_gQf_ij - pair.ji_gQf_ij));
    }

    // print the pair energies
    if (world.rank() == 0) {
        printf("current energy %2d %2d %12.8f %12.8f\n", pair.i, pair.j,
               pair.e_singlet, pair.e_triplet);
    }

    END_TIMER(world, "compute MP2 energy");
    // return the total energy of this pair
    return pair.e_singlet + pair.e_triplet;
}

/// apply the exchange operator on an orbital

/// @param[in]	phi the orbital
/// @param[in]	hc	hermitian conjugate -> swap bra and ket space
/// @return 	Kphi
real_function_3d MP2::K(const real_function_3d& phi, const bool hc) const {

    real_function_3d result = real_factory_3d(world);

    // multiply rhs with R2orbitals (the bra space)
    vecfuncT R2rhs;
    if (not hc)
        R2rhs = mul(world, phi, hf->R2orbitals());
    else
        R2rhs = mul(world, phi, hf->orbitals());

    // apply the poisson kernel and sum up
    for (std::size_t k = 0; k < hf->nemos().size(); ++k) {
        result += hf->nemo(k) * (*poisson)(R2rhs[k]);
    }
    return result;
}

/// apply the exchange operator on a pair function

/// @param[in]	phi the pair function
/// @param[in]	is_symmetric is the function symmetric wrt particle exchange
/// @return 	(K1 + K2) |phi >
real_function_6d MP2::K(const real_function_6d& phi,
                        const bool is_symmetric) const {

	// first particle 1
	real_function_6d result=apply_exchange_vector(phi,1);
        if (is_symmetric) {
		result = result + swap_particles(result);
        } else {
		result = result + apply_exchange_vector(phi,2);
        }
	result.truncate();
    return result;
}

/// apply the Coulomb operator a on orbital

/// @param[in]	phi the orbital
/// @return 	Jphi
real_function_3d MP2::J(const real_function_3d& phi) const {
    return (hf->get_coulomb_potential() * phi).truncate();
}


real_function_6d MP2::apply_exchange_vector(const real_function_6d& phi,
		const int particle) const {

	timer timer(world);
	// prepare all constituent functions
	std::vector<real_function_3d> ket=copy(world,hf->nemos());
	std::vector<real_function_3d> bra=copy(world,hf->R2orbitals());
    phi.change_tree_state(redundant,false);
//	phi.get_impl()->make_redundant(false);
	make_redundant(world,ket,false);
	make_redundant(world,bra,false);
	world.gop.fence();
	timer.tag("prepare K");

	// multiply pair function with bra
	std::vector<real_function_6d> result1(hf->nocc());
	for (auto& r : result1) {
	    r.set_impl(phi);
	    r.get_impl()->undo_redundant(true);
	}
	for (int i=0; i<hf->nocc(); ++i) {
		result1[i].get_impl()->multiply(phi.get_impl().get(),bra[i].get_impl().get(),particle);
	}
	world.gop.fence();
	phi.get_impl()->undo_redundant(false);
    truncate(world,result1);
    reduce_rank(result1);

	timer.tag("multiply bra");
	print_size(world,result1,"after multiplication and truncation");

	// apply exchange operator
	real_convolution_3d op = CoulombOperator(world, 0.0001,
			hf->get_calc().param.econv());
	op.particle() = particle;
	op.destructive()=true;

    std::vector<real_function_6d> result2=apply(world,op,result1);
	timer.tag("apply Green's function");

	truncate(world,result2);
    for (auto& v : result2) v.reduce_rank(false);
	world.gop.fence();

	timer.tag("truncate");
	print_size(world,result2,"after truncate");

	// multiply with ket
	std::vector<real_function_6d> result3(hf->nocc());
	make_redundant(world,result2);
	for (auto& r : result3) r.set_impl(phi);
	for (int i=0; i<hf->nocc(); ++i) {
		result3[i].get_impl()->multiply(result2[i].get_impl().get(),ket[i].get_impl().get(),particle);
	}

	world.gop.fence();
    truncate(result3);
    reduce_rank(result3);
	timer.tag("multiply ket");

	real_function_6d result4=real_factory_6d(world);
	for (auto r : result3) result4=result4+r;
	timer.tag("add all up");
	return result4;
}


/// apply the exchange operator on f

/// if the exchange operator is similarity transformed (R-1 K R) the
/// orbital spaces differ for the orbitals underneath the integral sign
/// R-1 K R = \phi-ket(1) * \int \phi-bra(1') * f(1',2)
/// @param[in]  f   the pair function
/// @param[in]  orbital_bra the orbital underneath the integral sign (typically R2orbitals)
/// @param[in]  orbital_ket the orbital to be pre-multiplied with (typically orbitals)
/// @return     the pair function, on which the exchange operator has been applied
real_function_6d MP2::apply_exchange(const real_function_6d& f,
                                     const real_function_3d& orbital_ket,
                                     const real_function_3d& orbital_bra, const int particle) const {

    MADNESS_ASSERT((particle == 1) or (particle == 2));
    real_convolution_3d op = CoulombOperator(world, 0.0001,
                                             hf->get_calc().param.econv());
    op.particle() = particle;
    op.destructive() = true;

    //            if (world.rank()==0) printf("start multiplication before K at time %.1f\n",wall_time());

    // multiply the orbital to the pair function
    //            real_function_6d x=(particle==1)
    //            		? CompositeFactory<double,6,3>(world).ket(copy(f)).V_for_particle1(copy(orbital_bra))
    //            		: CompositeFactory<double,6,3>(world).ket(copy(f)).V_for_particle2(copy(orbital_bra));
    //            x.fill_tree().truncate();
    real_function_6d x =
            multiply(copy(f), copy(orbital_bra), particle).truncate();

    // apply the Poisson operator
    //            if (world.rank()==0) printf("start exchange at time %.1f\n",wall_time());
    load_balance(x, false);
    x = op(x).truncate();

    // do the final multiplication with the orbital
    //            if (world.rank()==0) printf("start multiplication after K at time %.1f\n",wall_time());
    //            real_function_6d result= (particle==1)
    //            		? CompositeFactory<double,6,3>(world).ket(copy(x)).V_for_particle1(copy(orbital_ket))
    //            		: CompositeFactory<double,6,3>(world).ket(copy(x)).V_for_particle2(copy(orbital_ket));
    //            result.fill_tree().truncate().reduce_rank();
    real_function_6d result =
            multiply(copy(x), copy(orbital_ket), particle).truncate();

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
std::vector<real_function_3d> MP2::make_chi(const real_function_3d& phi,
                                            const real_convolution_3d& op, const bool hc) const {

    const double tol = 0.0;
    MADNESS_ASSERT(hf->get_calc().param.spin_restricted());
    std::vector<real_function_3d> psif;
    if (hc)
        psif = mul_sparse(world, phi, hf->nemos(), tol);
    else
        psif = mul_sparse(world, phi, hf->R2orbitals(), tol);
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
std::vector<real_function_3d> MP2::make_xi(const real_function_3d& phi_i,
                                           const real_function_3d& phi_j, const real_convolution_3d& op,
                                           const bool hc) const {
    const double tol = 0.0;
    return mul_sparse(world, phi_j, make_chi(phi_i, op, hc), tol);
}

/// apply the operator K on the reference and multiply with f; fK |phi^0>

/// @param[in]  i   index of orbital i
/// @param[in]  j   index of orbital j
/// @return     the function f12 (K(1)+K(2))|phi^0>
real_function_6d MP2::make_fKphi0(const int i, const int j) const {
    const real_function_3d& phi_i = hf->nemo(i);
    const real_function_3d& phi_j = hf->nemo(j);

    const real_function_3d Kphi_i = K(phi_i);
    const real_function_3d Kphi_j = K(phi_j);

    real_function_6d fKphi0a = CompositeFactory<double, 6, 3>(world).g12(
            corrfac.f()).particle1(copy(phi_i)).particle2(copy(Kphi_j));
    fKphi0a.fill_tree().truncate();
    real_function_6d fKphi0b = CompositeFactory<double, 6, 3>(world).g12(
            corrfac.f()).particle1(copy(Kphi_i)).particle2(copy(phi_j));
    fKphi0b.fill_tree().truncate();

    real_function_6d fKphi0 = (fKphi0a + fKphi0b).truncate();
    return fKphi0;
}

/// return the function (J(1)-K(1)) |phi0> as on-demand function

/// @param[in]	hc		compute hermitian conjugate -> swap bra and ket space
real_function_6d MP2::JK1phi0_on_demand(const int i, const int j,
                                        const bool hc) const {
    real_function_3d phi_i, phi_j;
    if (not hc) {
        phi_i = hf->nemo(i);
        phi_j = hf->nemo(j);
    } else {
        phi_i = hf->R2orbital(i);
        phi_j = hf->R2orbital(j);
    }

    const real_function_3d JKphi_i = J(phi_i) - K(phi_i, hc);

    real_function_6d tmp1 = CompositeFactory<double, 6, 3>(world).particle1(
            copy(JKphi_i)).particle2(copy(phi_j));
    return tmp1;
}

/// return the function (J(2)-K(2)) |phi0> as on-demand function
real_function_6d MP2::JK2phi0_on_demand(const int i, const int j,
                                        const bool hc) const {
    real_function_3d phi_i, phi_j;
    if (not hc) {
        phi_i = hf->nemo(i);
        phi_j = hf->nemo(j);
    } else {
        phi_i = hf->R2orbital(i);
        phi_j = hf->R2orbital(j);
    }

    const real_function_3d JKphi_j = J(phi_j) - K(phi_j, hc);

    real_function_6d tmp = CompositeFactory<double, 6, 3>(world).particle1(
            copy(phi_i)).particle2(copy(JKphi_j));
    return tmp;
}

/// return the function |phi0> as on-demand function
real_function_6d MP2::phi0_on_demand(const int i, const int j) const {
    const real_function_3d& phi_i = hf->orbital(i);
    const real_function_3d& phi_j = hf->orbital(j);

    real_function_6d tmp1 = CompositeFactory<double, 6, 3>(world).particle1(
            copy(phi_i)).particle2(copy(phi_j));
    return tmp1;
}

/// return the function |F1F2> as on-demand function
real_function_6d MP2::nemo0_on_demand(const int i, const int j) const {
    return CompositeFactory<double, 6, 3>(world).particle1(copy(hf->nemo(i))).particle2(
            copy(hf->nemo(j)));
}

/// multiply the given function with the 0th order Hamiltonian, exluding the 0th order energy

/// @param[in]  f   the function we apply H^0 on
/// @return     the function g=H^0 f, which is NOT orthogonalized against f
real_function_6d MP2::multiply_with_0th_order_Hamiltonian(
        const real_function_6d& f, const int i, const int j) const {

    START_TIMER(world);

    // the purely local part: Coulomb and U2
    real_function_3d v_local = hf->get_coulomb_potential()
                               + hf->nemo_ptr->ncf->U2();
    if (param.do_oep()) {
        real_function_3d voep = real_factory_3d(world);
        load(voep, "mRKS_potential_final");
        v_local += voep;
    }

    v_local.print_size("vlocal");
    f.print_size("u");

    // screen the construction of Vphi: do only what is needed to
    // get an accurate result of the BSH operator
    const double eps = zeroth_order_energy(i, j);
    real_convolution_6d op_mod = BSHOperator<6>(world, sqrt(-2 * eps), lo,
                                                bsh_eps);
    op_mod.modified() = true;
    real_function_6d vphi = CompositeFactory<double, 6, 3>(world).ket(copy(f)).V_for_particle1(
            copy(v_local)).V_for_particle2(copy(v_local));
    vphi.fill_tree(op_mod);
    vphi.print_size("vphi: local parts");

    // the part with the derivative operators: U1
    std::vector<real_function_6d> Drhs(6);
    for (int axis = 0; axis < 6; ++axis) {
        real_derivative_6d D = free_space_derivative<double, 6>(world,axis);
        Drhs[axis] = D(f).truncate();
    }
//	const std::vector<real_function_6d> Drhs = truncate(grad(f));

    std::vector<real_function_3d> U1_6d;
    for (int axis = 0; axis < 6; ++axis) {
        // note integer arithmetic
        U1_6d.push_back(copy(hf->nemo_ptr->ncf->U1(axis % 3)));
    }

    double tight_thresh = std::min(FunctionDefaults<6>::get_thresh(), 1.e-4);
    std::vector<real_function_6d> x(6);
    for (int axis = 0; axis < 6; ++axis) {

        if (axis / 3 + 1 == 1) {
            x[axis] =CompositeFactory<double, 6, 3>(world).ket(Drhs[axis])
                    .V_for_particle1(U1_6d[axis]).thresh(tight_thresh);

        } else if (axis / 3 + 1 == 2) {
            x[axis] =CompositeFactory<double, 6, 3>(world).ket(Drhs[axis])
                    .V_for_particle2(U1_6d[axis]).thresh(tight_thresh);
        }
    }
    world.gop.fence();
    for (auto xx : x) {
        CompositeFunctorInterface<double, 6, 3>* func=
                dynamic_cast<CompositeFunctorInterface<double, 6, 3>* >(&(*xx.get_impl()->get_functor()));
        func->make_redundant(false);
    }
    world.gop.fence();
    for (auto xx : x) xx.fill_tree(op_mod,false);
    world.gop.fence();
    print_size(world,x,"x after fill_tree");

    set_thresh(world,x,FunctionDefaults<6>::get_thresh());
    for (auto xx : x) vphi += xx;
    vphi.truncate().reduce_rank();

    vphi.print_size("(U_nuc + J) |ket>:  made V tree");
    END_TIMER(world, "apply (U + J) |ket>");

    // and the exchange
    START_TIMER(world);
    if (not param.do_oep()) vphi = (vphi - K(f, i == j)).truncate().reduce_rank();
    vphi.print_size("(U_nuc + J - K) |ket>:  made V tree");
    END_TIMER(world, "apply K |ket>");

    return vphi;
}

/// add the coupling terms for local MP2

/// @return \sum_{k\neq i} f_ki |u_kj> + \sum_{l\neq j} f_lj |u_il>
void MP2::add_local_coupling(const Pairs<ElectronPair>& pairs,
                             Pairs<real_function_6d>& coupling) const {

    if (world.rank() == 0) print("adding local coupling");

    // temporarily make all N^2 pair functions
    typedef std::map<std::pair<int, int>, real_function_6d> pairsT;
    pairsT quadratic;
    for (int k = param.freeze(); k < hf->nocc(); ++k) {
        for (int l = param.freeze(); l < hf->nocc(); ++l) {
            if (l >= k) {
                quadratic[std::make_pair(k, l)] = pairs(k, l).function;
            } else {
                quadratic[std::make_pair(k, l)] = swap_particles(pairs(l, k).function);
            }
        }
    }
    for (pairsT::iterator it = quadratic.begin(); it != quadratic.end(); ++it) {
        it->second.compress(false);
    }
    world.gop.fence();

    // the coupling matrix is the Fock matrix, skipping diagonal elements
    tensorT fock1 = get_fock_matrix();
    for (int k = 0; k < hf->nocc(); ++k) {
        if (fock1(k, k) > 0.0) MADNESS_EXCEPTION("positive orbital energies", 1);
        fock1(k, k) = 0.0;
    }

    for (int i = param.freeze(); i < hf->nocc(); ++i) {
        for (int j = i; j < hf->nocc(); ++j) {
            coupling(i, j) = real_factory_6d(world).compressed();
        }
    }

    for (int i = param.freeze(); i < hf->nocc(); ++i) {
        for (int j = i; j < hf->nocc(); ++j) {
            for (int k = param.freeze(); k < hf->nocc(); ++k) {
                if (fock1(k, i) != 0.0) {
                    coupling(i, j).gaxpy(1.0, quadratic[std::make_pair(k, j)], fock1(k, i), false);
                }
            }

            for (int l = param.freeze(); l < hf->nocc(); ++l) {
                if (fock1(l, j) != 0.0) {
                    coupling(i, j).gaxpy(1.0, quadratic[std::make_pair(i, l)], fock1(l, j), false);
                }
            }
            world.gop.fence();
            const double thresh = FunctionDefaults<6>::get_thresh();
            coupling(i, j).truncate(thresh*0.1);
        }
    }
    world.gop.fence();
}


}

