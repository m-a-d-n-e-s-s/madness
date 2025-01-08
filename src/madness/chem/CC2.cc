/*
 * CC2.cc
 *
 *  Created on: Aug 17, 2015
 *      Author: kottmanj
 */


#include<madness/chem/CC2.h>
#include<madness/mra/commandlineparser.h>
#include "MolecularOrbitals.h"
#include "localizer.h"
#include <timing_utilities.h>

namespace madness {

/// solve the CC2 ground state equations, returns the correlation energy
void
CC2::solve() {
    if (parameters.test()) CCOPS.test();

    const CalcType ctype = parameters.calc_type();

    Tensor<double> fmat=nemo->compute_fock_matrix(nemo->get_calc()->amo,nemo->get_calc()->aocc);
    long nfrozen=Localizer::determine_frozen_orbitals(fmat);
    parameters.set_derived_value<long>("freeze",nfrozen);
    if (not check_core_valence_separation(fmat)) enforce_core_valence_separation(fmat);

    MolecularOrbitals<double, 3> dummy_mo(nemo->get_calc()->amo, nemo->get_calc()->aeps);
    dummy_mo.print_frozen_orbitals(parameters.freeze());

    CCOPS.reset_nemo(nemo);
    CCOPS.get_potentials.parameters=parameters;
    CCOPS.update_intermediates(CCOPS.mo_ket());

    // info keep information on the MOs and the molecular coordinates
    Info info;
    info=CCOPS.update_info(parameters,nemo);
    info.intermediate_potentials=CCIntermediatePotentials(parameters);

    // check if all pair functions have been loaded or computed
    auto all_pairs_exist = [](const Pairs<CCPair>& pairs) {
        for (const auto& p : pairs.allpairs) {
            if (p.second.functions.size()==0) return false;
            if (not p.second.function().is_initialized()) return false;
            if (not p.second.constant_part.is_initialized()) return false;
        }
        return true;
    };

    // doubles for ground state
    Pairs<CCPair> mp2pairs, cc2pairs;
    // singles for ground state
    CC_vecfunction cc2singles(PARTICLE);

    // Pairs structure to vector if necessary
    const std::size_t nfreeze=parameters.freeze();
    const int nocc=CCOPS.mo_ket().size();
    triangular_map=PairVectorMap::triangular_map(nfreeze,nocc);

    double mp2_energy=0.0, cc2_energy=0.0, mp3_energy=0.0;

    bool need_tdhf=parameters.response();
    bool need_mp2=(ctype==CT_MP2 or ctype==CT_CISPD or ctype==CT_ADC2 or ctype==CT_MP3);
    bool need_cc2=(ctype==CT_LRCC2 or ctype==CT_CC2);

    // check for restart data for CC2, otherwise use MP2 as guess
    if (need_cc2) {
        initialize_pairs(cc2pairs, GROUND_STATE, CT_CC2, cc2singles, CC_vecfunction(RESPONSE), 0, info);
        if (not all_pairs_exist(cc2pairs)) need_mp2=true;
    }

    if (need_tdhf) {
        tdhf->prepare_calculation();
        MADNESS_CHECK(tdhf->get_parameters().freeze()==parameters.freeze());
        auto roots=tdhf->solve_cis();
        tdhf->analyze(roots);
    }

    if (need_mp2) {
        bool restarted=initialize_pairs(mp2pairs, GROUND_STATE, CT_MP2, CC_vecfunction(PARTICLE), CC_vecfunction(RESPONSE), 0, info);
        if (restarted and parameters.no_compute_mp2()) {
//            for (auto& pair : mp2pairs.allpairs) mp2_energy+=CCOPS.compute_pair_correlation_energy(pair.second);
        } else {
            mp2_energy = solve_mp2_coupled(mp2pairs, info);
            output_calc_info_schema("mp2",mp2_energy);
        }
        output.section(assign_name(CT_MP2) + " Calculation Ended !");
        if (world.rank() == 0) {
            printf_msg_energy_time("MP2 correlation energy",mp2_energy,wall_time());
//            std::cout << std::fixed << std::setprecision(10) << " MP2 Correlation Energy =" << mp2_energy << "\n";
	    }
    }

    if (need_cc2) {
        if (world.rank()==0) print_header2("computing the CC2 correlation energy");
        // check if singles or/and doubles to restart are there
        cc2singles=initialize_singles(world,CT_CC2, PARTICLE);

        // use mp2 as cc2 guess
        if (not all_pairs_exist(cc2pairs)) {
            if (world.rank()==0) print("could not find CC2 pairs -- using MP2 as guess");
            for (auto& tmp:mp2pairs.allpairs) {
                const size_t i = tmp.second.i;
                const size_t j = tmp.second.j;
                cc2pairs(i, j).update_u(tmp.second.function());
            }
        }

        if (parameters.no_compute_cc2()) {
            if (world.rank()==0) print("found no_compute_cc2 key -- no recomputation of singles or doubles or energy");
        } else {
            cc2_energy = solve_cc2(cc2singles, cc2pairs, info);
            output_calc_info_schema("cc2",cc2_energy);
        }

        if (world.rank() == 0) {
            printf_msg_energy_time("CC2 correlation energy",cc2_energy,wall_time());
            std::cout << std::fixed << std::setprecision(10) << " CC2 Correlation Energy =" << cc2_energy << "\n";
        }
        if (world.rank()==0) print_header2("end computing the CC2 correlation energy");
    }

    if (ctype == CT_LRCCS) {
        ;   // we're good
    } else if (ctype == CT_MP2) {
        ;   // we're good
    } else if (ctype == CT_MP3) {
        mp3_energy=compute_mp3(mp2pairs);
        double hf_energy=nemo->value();
        if (world.rank()==0) {
            printf_msg_energy_time("MP3 energy contribution",mp3_energy,wall_time());
            printf("final hf/mp2/mp3/total energy %12.8f %12.8f %12.8f %12.8f\n",
                    hf_energy,mp2_energy,mp3_energy,hf_energy+mp2_energy+mp3_energy);
            output_calc_info_schema("mp3",mp3_energy);
        }
    } else if (ctype == CT_CC2) {
        ;   // we're good
    } else if (ctype == CT_CISPD) {
        CCTimer time(world, "whole CIS(D) Calculation");

        auto vccs = solve_ccs();

        CCTimer time_cispd(world, "Time CIS(D) Response");
        std::vector<std::pair<double, double> > cispd_results;

        for (size_t k = 0; k < parameters.excitations().size(); k++) {

            CC_vecfunction& ccs = vccs[k];
            const size_t excitation = parameters.excitations()[k];
            CCTimer time_ex(world, "CIS(D) for Excitation " + std::to_string(int(excitation)));

            // check the convergence of the cis function (also needed to store the ccs potential) and to recalulate the excitation energy
            iterate_ccs_singles(ccs, info);

            Pairs<CCPair> cispd;
            initialize_pairs(cispd, EXCITED_STATE, CT_CISPD, CC_vecfunction(PARTICLE), ccs, excitation, info);

            const double ccs_omega = ccs.omega;
            const double cispd_omega = solve_cispd(cispd, mp2pairs, ccs);

            cispd_results.push_back(std::make_pair(ccs_omega, cispd_omega));
            time_ex.info();
        }

        output.section("CIS(D) Calculation Ended");
        for (size_t i = 0; i < cispd_results.size(); i++) {
            if (world.rank() == 0) {
                std::cout << std::fixed << std::setprecision(10) << "\n"
                          << "--------------------------------\n"
                          << "Excitation " << parameters.excitations()[i] << "\n"
                          << "CIS   =" << cispd_results[i].first << "\n"
                          << "CIS(D)=" << cispd_results[i].second << "\n"
                          << "Delta =" << cispd_results[i].second - cispd_results[i].first << "\n"
                          << "--------------------------------\n";
            }
        }
        time_cispd.info();
        time.info();

    } else if (ctype == CT_ADC2) {
        // we will never need the GS singles, but we use the CC2 potential functions so we initialize all gs singles potentials to zero
        CCOPS.update_intermediates(
                CC_vecfunction(zero_functions<double, 3>(world, CCOPS.get_active_mo_ket().size()), PARTICLE,
                               parameters.freeze()));
        output.section("ADC(2) Calculation");
        CCTimer time(world, "Whole ADC(2) Calculation");

        auto vccs = solve_ccs();

        CCTimer time_ex(world, "ADC(2) Calculation");
        output.section("ADC(2): Calculating ADC(2) Correction to CIS");
        std::vector<std::vector<double> > adc2_results;
        for (size_t k = 0; k < parameters.excitations().size(); k++) {

            CC_vecfunction& ccs = vccs[k];
            const size_t excitation = parameters.excitations()[k];
            CCTimer time_ex(world, "ADC(2) for Excitation " + std::to_string(int(excitation)));

            // check the convergence of the cis function (also needed to store the ccs potential) and to recalulate the excitation energy
            CC_vecfunction dummy = copy(ccs);
            iterate_ccs_singles(dummy, CCOPS.info);
            ccs.omega = dummy.omega; // will be overwritten soon
            output("Changes not stored!");

            Pairs<CCPair> xpairs;
            const bool restart = initialize_pairs(xpairs, EXCITED_STATE, CT_ADC2, CC_vecfunction(PARTICLE), ccs, excitation, CCOPS.info);

            // if no restart: Calculate CIS(D) as first guess
            const double ccs_omega = ccs.omega;
            double cispd_omega = 0.0;
            if (not restart) {
                output.section("No Restart-Pairs found: Calculating CIS(D) as first Guess");
                Pairs<CCPair> cispd;
                initialize_pairs(cispd, EXCITED_STATE, CT_CISPD, CC_vecfunction(PARTICLE), ccs, excitation, CCOPS.info);
                cispd_omega = solve_cispd(cispd, mp2pairs, ccs);
                for (auto& tmp:cispd.allpairs) {
                    const size_t i = tmp.first.first;
                    const size_t j = tmp.first.second;
                    xpairs(i, j).update_u(cispd(i, j).function());
                }
            }

            iterate_adc2_singles(mp2pairs, ccs, xpairs, CCOPS.info);
            for (size_t iter = 0; iter < 10; iter++) {
                bool dconv = iterate_adc2_pairs(xpairs, ccs);
                bool sconv = iterate_adc2_singles(mp2pairs, ccs, xpairs, CCOPS.info);
                if (sconv and dconv) {
                    output("ADC(2) Converged");
                    break;
                } else output("Not yet converged");
            }

            output.section("ADC(2) For Excitation " + std::to_string(int(excitation)) + " ended");
            const double adc2_omega = ccs.omega;
            std::vector<double> resulti;
            resulti.push_back(ccs_omega);
            resulti.push_back(cispd_omega);
            resulti.push_back(adc2_omega);
            adc2_results.push_back(resulti);
            if (world.rank() == 0)
                std::cout << std::fixed << std::setprecision(10)
                          << std::setfill(' ') << std::setw(12) << "CIS" << std::setw(12) << "CIS(D)" << std::setw(12)
                          << "ADC(2)" << "\n"
                          << ccs_omega << ", " << cispd_omega << ", " << adc2_omega << "\n";


            time_ex.info();
        }

        output.section("ADC(2) Ended!");
        if (world.rank() == 0)
            std::cout << std::fixed << std::setprecision(10)
                      << std::setfill(' ') << std::setw(12) << "CIS" << std::setw(12) << "CIS(D)" << std::setw(12)
                      << "ADC(2)" << "\n";
        for (size_t i = 0; i < adc2_results.size(); i++) {
            if (world.rank() == 0)
                std::cout << std::fixed << std::setprecision(10)
                          << adc2_results[i][0] << ", " << adc2_results[i][1] << ", " << adc2_results[i][2] << "\n";
        }


        time_ex.info();

        time.info();
    } else if (ctype == CT_LRCC2) {
        CCTimer time(world, "Whole LRCC2 Calculation");

        // if (world.rank() == 0) print_header3("reiterating CCS");
        auto vccs = solve_ccs();
        // iterate_ccs_singles(vccs[0], info);
        // if (world.rank() == 0) print_header3("end reiterating CCS");

        // for solving the LRCC2 equations we need converged CC2 singles and the applied potential
        bool need_reiterate_cc2=(cc2singles.size()==0)
                or (not info.intermediate_potentials.potential_exists(cc2singles, POT_singles_));
        if (need_reiterate_cc2) {
            if (world.rank()==0) print_header3("reiterating CC2 singles for intermediate potentials");
            // print_header1("skipping reiteration of CC2 singles for the potential");
            iterate_cc2_singles(world, cc2singles, cc2pairs, info);
            if (world.rank()==0) print_header3("end reiterating CC2 singles");
        }


        for (size_t iexcitation = 0; iexcitation < vccs.size(); iexcitation++) {
            if (world.rank() == 0)
                print_header1("Solving LRCC2 for excitation " + std::to_string(iexcitation)
                    + " with omega " + std::to_string(vccs[iexcitation].omega));
            solve_lrcc2(cc2pairs,cc2singles,vccs[iexcitation],iexcitation,info);
       }

    } else MADNESS_EXCEPTION(("Unknown Calculation Type: " + assign_name(ctype)).c_str(), 1);

}

void CC2::output_calc_info_schema(const std::string model, const double& energy) const {
    if (world.rank()==0) {
        nlohmann::json j;
        j["model"]=model;
        j["driver"]="energy";
        j["return_energy"]=energy;
        j[model]=energy;
        update_schema(nemo->get_param().prefix()+".calc_info", j);
    }
}


bool CC2::check_core_valence_separation(const Tensor<double>& fmat) const {

    MolecularOrbitals<double, 3> mos(nemo->get_calc()->amo, nemo->get_calc()->aeps, {}, nemo->get_calc()->aocc, {});
    mos.recompute_localize_sets();
    return Localizer::check_core_valence_separation(fmat, mos.get_localize_sets(),true);
}


Tensor<double> CC2::enforce_core_valence_separation(const Tensor<double>& fmat) {

    if (nemo->get_param().localize_method()=="canon") {
        auto nmo=nemo->get_calc()->amo.size();
        Tensor<double> fmat1(nmo,nmo);
        for (size_t i=0; i<nmo; ++i) fmat1(i,i)=nemo->get_calc()->aeps(i);
        return fmat1;
    }

    MolecularOrbitals<double, 3> mos(nemo->get_calc()->amo, nemo->get_calc()->aeps, {}, nemo->get_calc()->aocc, {});
    mos.recompute_localize_sets();

    Localizer localizer(world,nemo->get_calc()->aobasis,nemo->get_calc()->molecule,nemo->get_calc()->ao);
    localizer.set_enforce_core_valence_separation(true).set_method(nemo->param.localize_method());
    localizer.set_metric(nemo->R);

    const auto lmo=localizer.localize(mos,fmat,true);

    //hf->reset_orbitals(lmo);
    nemo->get_calc()->amo=lmo.get_mos();
    nemo->get_calc()->aeps=lmo.get_eps();
    MADNESS_CHECK(size_t(nemo->get_calc()->aeps.size())==nemo->get_calc()->amo.size());
    //orbitals_ = nemo->R*nemo->get_calc()->amo;
    //R2orbitals_ = nemo->ncf->square()*nemo->get_calc()->amo;


    //fock.clear();

    if (world.rank()==0) print("localized fock matrix");
    Tensor<double> fock2;
    const tensorT occ2 = nemo->get_calc()->aocc;
    Tensor<double> fock_tmp2 = nemo->compute_fock_matrix(nemo->get_calc()->amo, occ2);
    fock2 = copy(fock_tmp2);
    if (world.rank() == 0 and nemo->get_param().nalpha() < 10) {
        if (world.rank()==0) print("The Fock matrix");
        if (world.rank()==0) print(fock2);
    }

    MADNESS_CHECK(Localizer::check_core_valence_separation(fock2,lmo.get_localize_sets()));
    // if (world.rank()==0) lmo.pretty_print("localized MOs");
    return fock2;

};

// Solve the CCS equations for the ground state (debug potential and check HF convergence)
std::vector<CC_vecfunction> CC2::solve_ccs() const {
    std::vector<CC_vecfunction> excitations=tdhf->get_converged_roots();

    // return only those functions which are demanded
    std::vector<CC_vecfunction> result;
    for (const auto& x:parameters.excitations()) {
        if (excitations.size() - 1 < x) MADNESS_EXCEPTION("Not Enough CIS Vectors to solve for the demanded CC2 vector",
                                                          1);
        result.push_back(excitations[x]);
    }
    // if (world.rank()==0) print_header3("Solution of the CCS equations");
    // tdhf->analyze(result);
    return result;
}

double CC2::solve_mp2_coupled(Pairs<CCPair>& doubles, Info& info) {

    if (world.rank()==0) print_header2(" computing the MP1 wave function");
    CC_vecfunction dummy_singles(PARTICLE);

    // make vector holding CCPairs for partitioner of MacroTask
    std::vector<CCPair> pair_vec=Pairs<CCPair>::pairs2vector(doubles,triangular_map);

    // read constant part from file
    if (parameters.no_compute_mp2_constantpart()) {
        if (world.rank()==0) print("Skipping MP2 constant part calculation");
        for (auto& c : pair_vec) {
            MADNESS_CHECK_THROW(c.constant_part.is_initialized(), "could not find constant part");
            // constant part is zero-order guess for pair.function
            if (not c.function_exists()) c.update_u(c.constant_part);
        }

    } else {

        if (world.rank()==0) {
            std::cout << std::fixed << std::setprecision(1) << "\nStarting constant part at time " << wall_time() << std::endl;
        }
        MacroTaskConstantPart t;
        MacroTask task(world, t);
        std::vector<Function<double,3>> gs_singles, ex_singles;         // dummy vectors
        std::vector<real_function_6d> result_vec = task(pair_vec, gs_singles, ex_singles, info) ;

        if (world.rank()==0) {
            std::cout << std::fixed << std::setprecision(1) << "\nFinished constant part at time " << wall_time() << std::endl;
            std::cout << std::fixed << std::setprecision(1) << "\nStarting saving pairs and energy calculation at time " << wall_time() << std::endl;
        }
        load_balance(world, result_vec);

        // transform vector back to Pairs structure
        for (size_t i = 0; i < pair_vec.size(); i++) {
            auto& pair=pair_vec[i];
            pair.constant_part = result_vec[i];
            pair.constant_part.truncate().reduce_rank();
            pair.constant_part.print_size("constant_part");

            // initialize pair function to constant part
            pair.update_u(pair.constant_part);
            save(pair.constant_part, pair.name() + "_const");

        }
    }

    auto compute_energy = [&](const std::vector<CCPair>& pair_vec, std::string msg="") {
        for (const auto& p : pair_vec) {
            p.function().print_size("function "+p.name());
            p.function().print_size("constant_part "+p.name());
            p.function().reconstruct();
            p.constant_part.reconstruct();
        }
        MacroTaskComputeCorrelationEnergy t;
        MacroTask task1(world, t);
        CC_vecfunction dummy_singles1(PARTICLE);

        auto pair_energies=task1(pair_vec, dummy_singles, info);
        // pair_energies is now scattered over the universe

        double total_energy=0.0;
        for ( auto& pair_energy : pair_energies) {
            double pe=pair_energy.get();
            total_energy += pe;
            if (world.rank()==0 and parameters.debug()) printf("pair energy for pair %12.8f\n", pe);
        }
        // pair_energy.get() invokes a broadcast from rank 0 to all other ranks

        if (not msg.empty() and world.rank()==0) printf("%s %12.8f\n", msg.c_str(), total_energy);
        return total_energy;
    };

    timer t1_energy(world);
    double total_energy=compute_energy(pair_vec,"initial MP2 energy");
    t1_energy.end("compute energy");

    auto solver= nonlinear_vector_solver<double,6>(world,pair_vec.size());
    solver.set_maxsub(parameters.kain_subspace());
    solver.do_print = (world.rank() == 0);


    for (size_t iter = 0; iter < parameters.iter_max_6D(); iter++) {
        if (world.rank()==0) print_header3("Starting iteration " + std::to_string(int(iter)) + " of MP2");

        timer t_iter(world);
        // compute the coupling between the pair functions
        Pairs<real_function_6d> coupling=compute_local_coupling(pair_vec, info);
        auto coupling_vec=Pairs<real_function_6d>::pairs2vector(coupling,triangular_map);
        if (parameters.debug()) print_size(world, coupling_vec, "couplingvector");
        t_iter.tag("compute coupling");

        MacroTaskIteratePair t;
        MacroTask task1(world, t);
        CC_vecfunction dummy_singles1(PARTICLE);
        const std::size_t maxiter=3;
        coupling_vec=change_tree_state(coupling_vec, reconstructed);
        for (auto& p : pair_vec) p.function().reconstruct();
        auto unew = task1(pair_vec, coupling_vec, dummy_singles1, dummy_singles1, info, maxiter);
        t_iter.tag("compute update");

        std::vector<real_function_6d> u;
        for (auto p : pair_vec) u.push_back(p.function());
        auto residual=truncate(u-unew,parameters.tight_thresh_6D());

        // some statistics
        auto [rmsrnorm, maxrnorm]=CCPotentials::residual_stats(residual);

        if (parameters.debug()) {
            CCSize sz;
            sz.add(u);
            sz.print(world,"size of u");
            print_size(world,u,"u");

            sz.add(u,unew,residual,pair_vec,coupling_vec);
            for (const auto& r : solver.get_rlist()) sz.add(r);
            for (const auto& uu : solver.get_ulist()) sz.add(uu);
            sz.print(world,"sizes before KAIN");
            task1.taskq_ptr->cloud.print_size(world);
        }
        // update the pair functions
        std::string use_kain;
        if (parameters.kain()) {
            use_kain="with KAIN";
            std::vector<real_function_6d> kain_update = solver.update(u, residual);
            MADNESS_CHECK_THROW(solver.get_rlist()[0][0].is_reconstructed(),"solver functions are not reconstructed");
            MADNESS_CHECK_THROW(solver.get_ulist()[0][0].is_reconstructed(),"solver functions are not reconstructed");
            truncate(kain_update);
            for (size_t i=0; i<pair_vec.size(); ++i) {
                kain_update[i].reduce_rank();
                if (parameters.debug()) kain_update[i].print_size("Kain-Update-Function");
                pair_vec[i].update_u(kain_update[i]);
            }
        } else {
            use_kain="without KAIN";
            truncate(unew);
            for (size_t i=0; i<pair_vec.size(); ++i)  pair_vec[i].update_u(unew[i]);
        }
        t_iter.tag("update pairs "+use_kain);

        // calculate energy and error and update pairs
        double old_energy = total_energy;
        for (auto& p : pair_vec) p.reconstruct();
        total_energy=compute_energy(pair_vec);
        t_iter.tag("compute energy");

		if (world.rank()==0) {
		    double delta=old_energy - total_energy;
		    CCPotentials::print_convergence("MP2 doubles",rmsrnorm,maxrnorm,delta,iter);
			printf("finished MP2 iteration %2d at time %8.1fs with energy  %12.8f\n",
					int(iter), wall_time(), total_energy);
		}

        bool converged = ((std::abs(old_energy - total_energy) < parameters.econv())
                          and (maxrnorm < parameters.dconv_6D()));

        // save latest iteration
        if (world.rank()==0) print("saving latest iteration of MP2 to file"+pair_vec.front().basename()+"_XX");
        for (const auto& pair : pair_vec) {
            save(pair.constant_part, pair.name() + "_const");
            save(pair.function(), pair.name());
        }
        //print pair energies if converged
        if (converged) {
            if (world.rank() == 0) std::cout << "\nPairs converged!\n";
            break;
        }
    }
    if (world.rank()==0) {
        std::cout << std::fixed << std::setprecision(1) << "\nFinished final energy calculation at time " << wall_time() << std::endl;
        print_header2("end computing the MP1 wave function");
    }

    doubles=Pairs<CCPair>::vector2pairs(pair_vec,triangular_map);
    return total_energy;
}


/// add the coupling terms for local MP2

/// @return \sum_{k\neq i} f_ki |u_kj> + \sum_{l\neq j} f_lj |u_il>
Pairs<real_function_6d> CC2::compute_local_coupling(const Pairs<real_function_6d>& pairs, const Info& info) {

    const int nmo = info.mo_ket.size();
    World& world=pairs.allpairs.begin()->second.world();

    // temporarily make all N^2 pair functions
    typedef std::map<std::pair<int, int>, real_function_6d> pairsT;
    pairsT quadratic;
    for (int k = info.parameters.freeze(); k < nmo; ++k) {
        for (int l = info.parameters.freeze(); l < nmo; ++l) {
            if (l >= k) {
                quadratic[std::make_pair(k, l)] = pairs(k, l);
            } else {
                quadratic[std::make_pair(k, l)] = swap_particles(pairs(l, k));
            }
        }
    }

    for (auto& q: quadratic) q.second.compress(false);
    world.gop.fence();

    // the coupling matrix is the Fock matrix, skipping diagonal elements
    // Tensor<double> fock1 = nemo->compute_fock_matrix(nemo->get_calc()->amo, nemo->get_calc()->aocc);
    Tensor<double> fock1 = copy(info.fock);
    for (int k = 0; k < nmo; ++k) {
        if (fock1(k, k) > 0.0) MADNESS_EXCEPTION("positive orbital energies", 1);
        fock1(k, k) = 0.0;
    }

    Pairs<real_function_6d> coupling;
    for (int i = info.parameters.freeze(); i < nmo; ++i) {
        for (int j = i; j < nmo; ++j) {
            coupling.insert(i, j, real_factory_6d(world).compressed());
        }
    }

    for (int i = info.parameters.freeze(); i < nmo; ++i) {
        for (int j = i; j < nmo; ++j) {
            for (int k = info.parameters.freeze(); k < nmo; ++k) {
                if (fock1(k, i) != 0.0) {
                    coupling(i, j).gaxpy(1.0, quadratic[std::make_pair(k, j)], fock1(k, i), false);
                }
            }

            for (int l = info.parameters.freeze(); l < nmo; ++l) {
                if (fock1(l, j) != 0.0) {
                    coupling(i, j).gaxpy(1.0, quadratic[std::make_pair(i, l)], fock1(l, j), false);
                }
            }
            world.gop.fence();
            const double thresh = FunctionDefaults<6>::get_thresh();
            coupling(i, j).truncate(thresh * 0.3).reduce_rank();
        }
    }
    world.gop.fence();
    return coupling;
}

double
CC2::solve_cispd(Pairs<CCPair>& cispd, const Pairs<CCPair>& mp2, const CC_vecfunction& ccs) {
    output.section("Solve CIS(D) for CIS Excitation energy " + std::to_string(double(ccs.omega)));
    MADNESS_ASSERT(ccs.type == RESPONSE);
    CCOPS.update_intermediates(ccs);

    for (auto& pairs:cispd.allpairs) {
        CCPair& pair = pairs.second;
        pair.bsh_eps = CCOPS.get_epsilon(pair.i, pair.j) + ccs.omega;
        if (size_t(parameters.only_pair().first) == pair.i and size_t(parameters.only_pair().second) == pair.j) {
            output("Found only_pair exception");
            update_constant_part_cispd(ccs, pair);
            iterate_pair(pair, ccs);
        } else if (parameters.no_compute_cispd()) output("Found no_compute_cispd key");
        else {
            update_constant_part_cispd(ccs, pair);
            iterate_pair(pair, ccs);
        }
        // test consitency of the two approaches
        if (parameters.debug() and parameters.thresh_6D() > 1.e-4)
            CCOPS.test_pair_consistency(pair.functions[0], pair.i, pair.j, ccs);
    }

    const double diff = CCOPS.compute_cispd_energy(ccs, mp2, cispd);
    CC_vecfunction empty(zero_functions<double, 3>(world, ccs.size()), PARTICLE, parameters.freeze());
    const double omega_cc2 = CCOPS.compute_cc2_excitation_energy(empty, ccs, mp2, cispd);
    output.section("CIS(D) Calculation for CIS Excitation " + std::to_string(double(ccs.omega)) + " ended");
    if (world.rank() == 0) {
        std::cout << std::fixed << std::setprecision(10)
                  << "CIS   =" << ccs.omega << "\n"
                  << "CIS(D)=" << ccs.omega + diff << "\n"
                  << "Diff  =" << diff
                  << "\nomega_cc2 =" << omega_cc2 << "\n\n\n";
    }

    return ccs.omega + diff;
}

bool
CC2::iterate_adc2_pairs(Pairs<CCPair>& cispd, const CC_vecfunction& ccs) {
    output.section("Solve ADC(2) for Excitation energy " + std::to_string(double(ccs.omega)));
    MADNESS_ASSERT(ccs.type == RESPONSE);
    CCOPS.update_intermediates(ccs);

    bool conv = true;
    for (auto& pairs:cispd.allpairs) {
        CCPair& pair = pairs.second;
        pair.bsh_eps = CCOPS.get_epsilon(pair.i, pair.j) + ccs.omega;
        update_constant_part_adc2(ccs, pair);
        conv = iterate_pair(pair, ccs);
    }

    return conv;
}

bool
CC2::iterate_lrcc2_pairs(World& world, const CC_vecfunction& cc2_s,
                         const CC_vecfunction lrcc2_s, Pairs<CCPair>& lrcc2_d, const Info& info) {
    // output.section("Solve LRCC2 for Excitation energy " + std::to_string(double(lrcc2_s.omega)));
    if (world.rank()==0) {
        print_header3("Solving LRCC2 doubles equations");
        print("starting at time ",wall_time());
        print("using macrotasks with redirected output");
    }
    MADNESS_ASSERT(lrcc2_s.type == RESPONSE);

    auto triangular_map=PairVectorMap::triangular_map(info.parameters.freeze(),info.mo_ket.size());
    auto pair_vec=Pairs<CCPair>::pairs2vector(lrcc2_d,triangular_map);

    for (auto& p: pair_vec) p.reconstruct();
    info.reconstruct();
    cc2_s.reconstruct();
    lrcc2_s.reconstruct();

    if (1) {

    // make new constant part
    MacroTaskConstantPart tc;
    MacroTask task(world, tc);
    auto cp = task(pair_vec, cc2_s.get_vecfunction(), lrcc2_s.get_vecfunction(), info) ;
    print_size(world,cp,"constant part in iter");

    for (int i=0; i<pair_vec.size(); ++i) {
        pair_vec[i].constant_part=cp[i];
        save(pair_vec[i].constant_part, pair_vec[i].name() + "_const");
    }
    } else {
        print("reading LRCC2 constant part from file");
        for (int i=0; i<pair_vec.size(); ++i) {
            real_function_6d tmp=real_factory_6d(world);
            load(tmp, pair_vec[i].name() + "_const");
            pair_vec[i].constant_part=tmp;
        }
    }

    // if no function has been computed so far use the constant part (first iteration)
    // for (auto& pair : pair_vec) if (not pair.function_exists()) pair.update_u(pair.constant_part);
    for (auto& pair : pair_vec) if (not pair.function_exists()) pair.update_u(real_factory_6d(world));

    for (const auto& p : pair_vec) p.constant_part.print_size("constant_part before iter");
    for (const auto& p : pair_vec) p.function().print_size("u before iter");

    // compute the coupling between the pair functions
    if (world.rank()==0) print("computing local coupling in the universe");
    Pairs<real_function_6d> coupling=compute_local_coupling(pair_vec, info);
    auto coupling_vec=Pairs<real_function_6d>::pairs2vector(coupling,triangular_map);
    reconstruct(world,coupling_vec);
    for (auto& p : pair_vec) p.reconstruct();

    if (info.parameters.debug()) print_size(world, coupling_vec, "couplingvector");

    // iterate the pair
    MacroTaskIteratePair t1;
    MacroTask task1(world, t1);
    // temporary fix: create dummy functions to that the cloud is not confused
    // real_function_6d tmp=real_factory_6d(world).functor([](const coord_6d& r){return 0.0;});
    // std::vector<real_function_6d> vdummy_6d(pair_vec.size(),tmp);         // dummy vectors
    const std::size_t maxiter=10;
    auto unew = task1(pair_vec, coupling_vec, cc2_s, lrcc2_s, info, maxiter);

    for (const auto& u : unew) u.print_size("u after iter");
    // get some statistics
    std::vector<Function<double,6>> uold;
    for (const auto & p : pair_vec) uold.push_back(p.function());
    auto residual=uold-unew;
    double nold=norm2(world,uold);
    double nnew=norm2(world,unew);
    print("norm(old), norm(new) ",nold,nnew);
    auto [rmsrnorm, rmsrmax] = CCPotentials::residual_stats(residual);
    if (world.rank()==0) CCPotentials::print_convergence("LRCC2 doubles",rmsrnorm, rmsrmax,0,0);

    // update the pair functions
    reconstruct(world,unew);        // saves a lot of memory!
    for (int i=0; i<pair_vec.size(); ++i) pair_vec[i].update_u(unew[i]);
    lrcc2_d=Pairs<CCPair>::vector2pairs(pair_vec,triangular_map);

    // save latest iteration
    if (world.rank()==0) print("saving latest iteration of LRCC2 to file");
    for (const auto& pair : pair_vec) {
        save(pair.constant_part, pair.name() + "_const");
        save(pair.function(), pair.name());
    }

    return (rmsrnorm<info.parameters.dconv_6D());
}


double
CC2::solve_cc2(CC_vecfunction& singles, Pairs<CCPair>& doubles, Info& info) const
{

    output.section("Solving CC2 Ground State");
    if (singles.size()==0) singles=initialize_singles_to_zero(world,CT_CC2, PARTICLE, info);

    MADNESS_ASSERT(singles.type == PARTICLE);
    CCTimer time(world, "CC2 Ground State");

    if (world.rank()==0) print("computing CC2 energy");
    double omega = CCPotentials::compute_cc2_correlation_energy(world, singles, doubles, info, "initial CC2 energy");

    CC_vecfunction ex_singles_dummy;

    // first singles iteration
    output.section("Initialize Singles to the Doubles");

    // given the doubles, we can solve the singles equations
    iterate_cc2_singles(world, singles, doubles, info);
    omega = CCPotentials::compute_cc2_correlation_energy(world, singles, doubles, info, "CC2 energy with converged singles");

    for (size_t iter = 0; iter < parameters.iter_max(); iter++) {
        CCTimer time_miter(world, "Macroiteration " + std::to_string(int(iter)) + " of CC2");
        output.section("Macroiteration " + std::to_string(int(iter)) + " of CC2");

        if (world.rank()==0) print("computing the constant part via macrotasks -- output redirected");
        timer timer1(world);

        std::vector<CCPair> pair_vec=Pairs<CCPair>::pairs2vector(doubles,triangular_map);
        for (auto& p : pair_vec) p.reconstruct();
        MacroTaskConstantPart t;
        MacroTask task(world, t);
        std::vector<real_function_6d> constant_part_vec = task(pair_vec, singles.get_vecfunction(),
            ex_singles_dummy.get_vecfunction(), info) ;
        for (int i=0; i<pair_vec.size(); ++i) pair_vec[i].constant_part=constant_part_vec[i];

        if (parameters.debug()) {
            for (auto& pair: pair_vec) pair.constant_part.print_size("size of constant part macrotask "+pair.name());
        }

        timer1.tag("computing constant part via macrotasks");


        // compute the coupling between the pair functions
        if (world.rank()==0) print("computing local coupling in the universe");
        Pairs<real_function_6d> coupling=compute_local_coupling(pair_vec, info);
        auto coupling_vec=Pairs<real_function_6d>::pairs2vector(coupling,triangular_map);
        timer1.tag("computing local coupling");

        if (world.rank()==0) print("update the pair functions via macrotasks -- output redirected");
        MacroTaskIteratePair t1;
        MacroTask task1(world, t1);
        CC_vecfunction dummy_ex_singles;
        std::vector<real_function_3d> vdummy_3d;         // dummy vectors
        const std::size_t maxiter=3;
        for (auto& p : pair_vec) p.reconstruct();
        for (auto& p : coupling_vec) p.reconstruct();
        info.reconstruct();
        auto unew = task1(pair_vec, coupling_vec, singles, dummy_ex_singles,
            info, maxiter);


        std::vector<real_function_6d> u_old;
        for (auto p : pair_vec) u_old.push_back(p.function());

        auto residual=u_old-unew;
        timer1.tag("computing pair function update via macrotasks");

        for (int i=0; i<pair_vec.size(); ++i) pair_vec[i].update_u(unew[i]);
        doubles=Pairs<CCPair>::vector2pairs(pair_vec,triangular_map);

        // save latest iteration
        if (world.rank()==0) print("saving latest iteration to file");
        for (const auto& pair : pair_vec) {
            pair.reconstruct();
            save(pair.constant_part, pair.name() + "_const");
            save(pair.function(), pair.name());
            singles.save_restartdata(world,CC2::singles_name(CT_CC2,singles.type));
        }
        timer1.tag("saving cc2 singles and doubles");

        auto [rmsrnorm,maxrnorm]=CCPotentials::residual_stats(residual);
        bool doubles_converged=rmsrnorm<parameters.dconv_6D();

        // check if singles converged
        const bool singles_converged = iterate_cc2_singles(world, singles, doubles, info);
        timer1.tag("computing cc2 singles");

        // check if energy converged
        const double omega_new = CCPotentials::compute_cc2_correlation_energy(world, singles, doubles, info);
        timer1.tag("computing cc2 energy");
        const double delta = omega_new - omega;
        const bool omega_converged(delta < parameters.econv());
        omega = omega_new;
        if (world.rank() == 0)
            std::cout << std::fixed << std::setprecision(10) << "Current Correlation Energy = " << omega << "\n";
        if (world.rank() == 0)
            std::cout << std::fixed << std::setprecision(10) << "Difference                  = " << delta << "\n";

        if (world.rank()==0) {
            CCPotentials::print_convergence("CC2 macro",rmsrnorm,maxrnorm,delta,iter);
            printf("finished CC2 macro iteration %2d at time %8.1fs with energy  %12.8f\n",
                    int(iter), wall_time(), omega);
        }
        if (doubles_converged and singles_converged and omega_converged) break;

        time_miter.info();
    }
    omega = CCPotentials::compute_cc2_correlation_energy(world, singles, doubles, info);
    output.section("CC2 Iterations Eneded");

    if (world.rank() == 0)
        std::cout << std::fixed << std::setprecision(10) << "Current Correlation Energy = " << omega << "\n";
    time.info();
    return omega;

}


/// solve the excited state LR-CC2 equations for a given excitation

/// @param[in] gs_doubles: the ground state doubles
/// @param[in] gs_singles: the ground state singles
/// @param[in] cis: the CIS singles
/// @param[in] excitation: the excitation number
/// @return a tuple with the excited state doubles, the excited state singles and the excitation energy
std::tuple<Pairs<CCPair>, CC_vecfunction, double>
CC2::solve_lrcc2(Pairs<CCPair>& gs_doubles, const CC_vecfunction& gs_singles, const CC_vecfunction& cis,
    const std::size_t excitation, Info& info) const {

    CCTimer time(world, "Whole LRCC2 Calculation");

    /// read LRCC2 singles from file or use the CIS vectors as guess
    auto ex_singles=initialize_singles(world,CT_LRCC2,RESPONSE, excitation);
    bool found_singles=(not (ex_singles.get_vecfunction().empty()));
    // if (not found_singles) {
        if (world.rank()==0) print("using CIS vectors as guess for the LRCC2 singles");
        ex_singles = copy(cis);
        iterate_ccs_singles(ex_singles, info);

        std::string filename=singles_name(CT_LRCC2,ex_singles.type,excitation);
        if (world.rank()==0) print("saving singles to disk",filename);
        ex_singles.save_restartdata(world,filename);
    // }

    Pairs<CCPair> ex_doubles;
    bool found_lrcc2d = initialize_pairs(ex_doubles, EXCITED_STATE, CT_LRCC2, gs_singles, ex_singles, excitation, info);

    if ((not found_singles) and found_lrcc2d) {
        iterate_lrcc2_singles(world, gs_singles, gs_doubles, ex_singles, ex_doubles, info);

        std::string filename1=singles_name(CT_LRCC2,ex_singles.type,excitation);
        if (world.rank()==0) print("saving singles to disk",filename1);
        ex_singles.save_restartdata(world,filename1);
    }

    if (not info.parameters.no_compute_lrcc2()) {
        for (size_t iter = 0; iter < parameters.iter_max(); iter++) {
            if (world.rank()==0) print_header2("Macroiteration " + std::to_string(int(iter)) + " of LRCC2 for excitation energy "+std::to_string(ex_singles.omega));

            bool dconv = iterate_lrcc2_pairs(world, gs_singles, ex_singles, ex_doubles, info);
            bool sconv = iterate_lrcc2_singles(world, gs_singles, gs_doubles, ex_singles, ex_doubles, info);

            std::string filename=singles_name(CT_LRCC2,ex_singles.type,excitation);
            if (world.rank()==0) print("saving singles to disk",filename);
            ex_singles.save_restartdata(world,filename);

            if (sconv and dconv) break;
        }
    }

    // std::tuple<Pairs<CCPair>, CC_vecfunction, double>

    auto result=std::make_tuple(ex_doubles, ex_singles, ex_singles.omega);

    const double omega_cc2 = ex_singles.omega;

    std::vector<std::pair<std::string, std::pair<double, double>>> timings;
    std::vector<std::pair<std::string, std::pair<double, double>>> results_ex;

    const std::string msg = "Excitation " + std::to_string(int(excitation));
    results_ex.push_back(std::make_pair(msg, std::make_pair(cis.omega, omega_cc2)));


    timings.push_back(std::make_pair("Whole LRCC2", time.current_time(true)));
    output.section("LRCC2 Finished");

    output("Response Results:");
    for (const auto& res : results_ex) {
        output << std::scientific << std::setprecision(10) << std::setw(20);

        output << assign_name(CT_LRCCS)+" excitation energy for root "
            << std::fixed << std::setprecision(1) << excitation << ": "
            << std::fixed << std::setprecision(10) << res.second.first << " Eh         "
            << res.second.first * constants::hartree_electron_volt_relationship << " eV   "
            << "\n";
        output << assign_name(CT_LRCC2)+" excitation energy for root "
            << std::fixed << std::setprecision(1) << excitation << ": "
            << std::fixed << std::setprecision(10) << res.second.second << " Eh         "
            << res.second.second * constants::hartree_electron_volt_relationship << " eV   "
            << "\n";
        output << std::scientific;
    }
    output("\nTimings");
    for (const auto& time : timings) {
        if (world.rank() == 0)
            std::cout << std::scientific << std::setprecision(2)
                << std::setfill(' ') << std::setw(15) << time.first
                << ": " << time.second.first << " (Wall), " << time.second.second << " (CPU)" << "\n";
    }

    return result;
};

bool CC2::iterate_pair(CCPair& pair, const CC_vecfunction& singles) const {
    output.section("Iterate Pair " + pair.name());
    if (pair.ctype == CT_CC2) MADNESS_ASSERT(singles.type == PARTICLE);
    if (pair.ctype == CT_CISPD) MADNESS_ASSERT(singles.type == RESPONSE);
    if (pair.ctype == CT_MP2) MADNESS_ASSERT(singles.get_vecfunction().empty());
    if (pair.ctype == CT_ADC2)MADNESS_ASSERT(singles.type == RESPONSE);

    real_function_6d constant_part = pair.constant_part;
    constant_part.truncate().reduce_rank();
    pair.function().truncate().reduce_rank();

    output.subsection("Converge pair " + pair.name() + " on constant singles potential");

    double bsh_eps = pair.bsh_eps; //CCOPS.get_epsilon(pair.i,pair.j)+omega;
    real_convolution_6d G = BSHOperator<6>(world, sqrt(-2.0 * bsh_eps), parameters.lo(), parameters.thresh_bsh_6D());
    G.destructive() = true;

    NonlinearSolverND<6> solver(parameters.kain_subspace());
    solver.do_print = (world.rank() == 0);

    bool converged = false;

    double omega = 0.0;
    Info info;
    info.mo_bra=CCOPS.mo_bra_.get_vecfunction();
    info.parameters=parameters;
    if (pair.type == GROUND_STATE) omega = CCOPS.compute_pair_correlation_energy(world, pair,singles, info);
    if (pair.type == EXCITED_STATE) omega = CCOPS.compute_excited_pair_energy(world, pair, singles, info);

    if (world.rank() == 0)
        std::cout << "Correlation Energy of Pair " << pair.name() << " =" << std::fixed << std::setprecision(10)
                  << omega << "\n";

    for (size_t iter = 0; iter < parameters.iter_max_6D(); iter++) {
        output.subsection(assign_name(pair.ctype) + "-Microiteration");
        CCTimer timer_mp2(world, "MP2-Microiteration of pair " + pair.name());


        CCTimer timer_mp2_potential(world, "MP2-Potential of pair " + pair.name());
        real_function_6d mp2_potential = -2.0 * CCOPS.fock_residue_6d(pair);
        if (parameters.debug()) mp2_potential.print_size(assign_name(pair.ctype) + " Potential");
        mp2_potential.truncate().reduce_rank();
        timer_mp2_potential.info(true, mp2_potential.norm2());

        CCTimer timer_G(world, "Apply Greens Operator on MP2-Potential of pair " + pair.name());
        const real_function_6d GVmp2 = G(mp2_potential);
        timer_G.info(true, GVmp2.norm2());

        CCTimer timer_addup(world, "Add constant parts and update pair " + pair.name());
        real_function_6d unew = GVmp2 + constant_part;
        unew.print_size("unew");
        unew = CCOPS.apply_Q12t(unew, CCOPS.mo_ket());
        unew.print_size("Q12unew");
        //unew.truncate().reduce_rank(); // already done in Q12 application at the end
        if (parameters.debug())unew.print_size("truncated-unew");
        const real_function_6d residue = pair.function() - unew;
        const double error = residue.norm2();
        if (parameters.kain()) {
            output("Update with KAIN");
            real_function_6d kain_update = copy(solver.update(pair.function(), residue));
            kain_update = CCOPS.apply_Q12t(kain_update, CCOPS.mo_ket());
            kain_update.truncate().reduce_rank();
            kain_update.print_size("Kain-Update-Function");
            pair.update_u(copy(kain_update));
        } else {
            output("Update without KAIN");
            pair.update_u(unew);
        }

        timer_addup.info(true, pair.function().norm2());

        double omega_new = 0.0;
        double delta = 0.0;
        if (pair.type == GROUND_STATE) omega_new = CCOPS.compute_pair_correlation_energy(world, pair, singles, info);
        else if (pair.type == EXCITED_STATE) omega_new = CCOPS.compute_excited_pair_energy(world, pair, singles, info);
        delta = omega - omega_new;

        const double current_norm = pair.function().norm2();

        omega = omega_new;
        if (world.rank() == 0) {
            std::cout << std::fixed
                      << std::setw(50) << std::setfill('#')
                      << "\n" << "Iteration " << iter << " of pair " << pair.name()
                      << std::setprecision(4) << "||u|| = " << current_norm
                      << "\n" << std::setprecision(10) << "error = " << error << "\nomega = " << omega << "\ndelta = "
                      << delta << "\n"
                      << std::setw(50) << std::setfill('#') << "\n";
        }


        output("\n--Iteration " + stringify(iter) + " ended--");
        save(pair.function(), pair.name());
        timer_mp2.info();
        if (fabs(error) < parameters.dconv_6D()) {
            output(pair.name() + " converged!");
            if (fabs(delta) < parameters.econv_pairs()) {
                converged = true;
                break;
            } else output("Energy not yet converged");
        } else output("Convergence for pair " + pair.name() + " not reached yet");
    }

    return converged;
}


bool CC2::iterate_singles(World& world, CC_vecfunction& singles, const CC_vecfunction singles2,
    const Pairs<CCPair>& gs_doubles, const Pairs<CCPair>& ex_doubles, const CalcType ctype, const std::size_t maxiter,
    Info& info) {
    CCMessenger output(world);
    if (world.rank()==0) print_header2("Iterating Singles for "+assign_name(ctype));
    CCTimer time_all(world, "Overall Iteration of " + assign_name(ctype) + "-Singles");
    if (singles.size()==0) singles=initialize_singles_to_zero(world,ctype, singles.type, info);

    // consistency checks
    switch (ctype) {
    case CT_CC2:
        if (singles.type != PARTICLE)
            output.warning("iterate_singles: CC2 demanded but singles are not of type PARTICLE");
        break;
    case CT_MP2: MADNESS_EXCEPTION("Demanded Singles Calculation for MP2 ????", 1);
        break;
    case CT_LRCC2:
        if (singles.type != RESPONSE or singles2.type != PARTICLE)
            output.warning("iterate_singles: CC2_response_ singles have wrong types");
        break;
    case CT_LRCCS:
        if (singles.type != RESPONSE)
            output.warning("iterate_singles: CCS_response_ singles have wrong types");
        break;
    case CT_CISPD: MADNESS_EXCEPTION("Demanded Singles Calculation for CIS(D)", 1);
        break;
    case CT_ADC2:
        MADNESS_ASSERT(singles.type == RESPONSE);
        break;
    case CT_TEST: MADNESS_EXCEPTION("Iterate Singles not implemented for Experimental calculation", 1);
        break;
    default: MADNESS_EXCEPTION(
            ("Unknown calculation type in iterate singles: " + assign_name(ctype)).c_str(), 1);
    }

    bool converged = true;


    CC_vecfunction old_singles(singles);
    for (auto& tmp : singles.functions)
        old_singles(tmp.first).function = copy(tmp.second.function);
    double old_omega=0.0;

    for (auto& p : gs_doubles.allpairs) p.second.reconstruct();
    for (auto& p : ex_doubles.allpairs) p.second.reconstruct();
    info.reconstruct();

    // KAIN solver
    typedef vector_function_allocator<double, 3> allocT;
    typedef XNonlinearSolver<std::vector<Function<double, 3> >, double, allocT> solverT;
    solverT solver(allocT(world, singles.size()));
    solver.do_print = ((world.rank() == 0) and info.parameters.debug());

    if (info.parameters.debug()) print_size(world, singles.get_vecfunction(), "singles before iteration");

    for (size_t iter = 0; iter < maxiter; iter++) {
        double omega = 0.0;
        if (ctype == CT_LRCC2) omega = singles.omega;
        else if (ctype == CT_LRCCS) omega = singles.omega;
        else if (ctype == CT_ADC2) omega = singles.omega;
        if ((world.rank()==0) and info.parameters.debug()) print("omega " ,omega);

        // get potentials using macrotasks
        CCTimer time_V(world, assign_name(ctype) + "-Singles Potential");
        vector_real_function_3d V;
        if (ctype == CT_CC2)
            V = CCPotentials::get_CC2_singles_potential_gs(world, singles, gs_doubles, info);
        else if (ctype == CT_LRCC2)
            V = CCPotentials::get_CC2_singles_potential_ex(world, singles2, gs_doubles, singles, ex_doubles, info);
        else if (ctype == CT_LRCCS)
            V = CCPotentials::get_CCS_potential_ex(world,singles,false, info);
        else if (ctype == CT_ADC2)
            V = CCPotentials::get_ADC2_singles_potential(world, gs_doubles, singles, ex_doubles, info);
        else MADNESS_EXCEPTION("iterate singles: unknown type", 1);

        if (info.parameters.debug()) madness::print_size(world, V, "final V in iterate_singles w/o coupling");

        // add local coupling
        V-=compute_local_coupling(singles.get_vecfunction(),info);
        truncate(world, V);
        time_V.info(info.parameters.debug(), norm2(world, V));

        if (info.parameters.debug()) madness::print_size(world, V, "final V in iterate_singles with coupling");

        // update excitation energy
        if (ctype==CT_LRCC2 or ctype==CT_LRCCS or ctype==CT_ADC2) {
            old_omega=omega;
            omega = CCPotentials::compute_cis_expectation_value(world, singles, V, info.parameters.debug(), info);
            singles.omega = omega;
        }
        if (world.rank()==0 and info.parameters.debug())
            print("omega entering the update in the singles" ,omega);

        // make bsh operators
        scale(world, V, -2.0); // moved to BSHApply
        std::vector<std::shared_ptr<SeparatedConvolution<double, 3> > > G(singles.size());
        for (size_t i = 0; i < G.size(); i++) {
            const double bsh_eps = info.orbital_energies[i + info.parameters.freeze()] + omega;
            G[i] = std::shared_ptr<SeparatedConvolution<double, 3> >(
                BSHOperatorPtr3D(world, sqrt(-2.0 * bsh_eps), info.parameters.lo(), info.parameters.thresh_bsh_3D()));
        }
        world.gop.fence();

        // apply bsh operators
        CCTimer time_applyG(world, "Apply G-Operators");
        vector_real_function_3d GV = apply<SeparatedConvolution<double, 3>, double, 3>(world, G, V);
        world.gop.fence();
        time_applyG.info(info.parameters.debug());

        // apply Q-Projector to result
        QProjector<double,3> Q(info.mo_bra,info.mo_ket);
        GV = Q(GV);

        // Normalize Singles if it is excited state
        if (ctype == CT_LRCCS or ctype == CT_LRCC2 or ctype == CT_ADC2) {
            if (info.parameters.debug()) print("Normalizing new singles");
            const double norm=inner(GV,info.R_square*GV);
            scale(world, GV, 1.0 / norm);
        } else {
            if (info.parameters.debug()) print("Singles not normalized");
        }

        // residual
        const vector_real_function_3d residual = sub(world, singles.get_vecfunction(), GV);

        // information
        const Tensor<double> R2xinnerx = inner(world, info.R_square*singles.get_vecfunction(),
                                               singles.get_vecfunction());
        const Tensor<double> R2GVinnerGV = inner(world, info.R_square*GV, GV);
        const Tensor<double> R2rinnerr = inner(world, info.R_square*residual, residual);
        const double R2vector_error = sqrt(R2rinnerr.sum());
        auto [rmsresidual, maxresidual]=CCPotentials::residual_stats(residual);

        // print information
        if (info.parameters.debug() and (world.rank()==0)) {
            std::cout << "\n\n-----Results of current interation:-----\n";
            std::cout << "\nName: ||" << singles.name(0) << "||, ||GV" << singles.name(0) << ", ||residual||" << "\n";
            std::cout << singles.name(0) << ": " << std::scientific << std::setprecision(info.parameters.output_prec())
                << sqrt(R2xinnerx.sum()) << ", " << sqrt(R2GVinnerGV.sum()) << ", " << sqrt(R2rinnerr.sum())
                << "\n----------------------------------------\n";
            for (size_t i = 0; i < GV.size(); i++) {
                std::cout << singles(i + info.parameters.freeze()).name() << ": " << std::scientific
                    << std::setprecision(info.parameters.output_prec())
                    << sqrt(R2xinnerx(i)) << ", " << sqrt(R2GVinnerGV(i)) << ", " << sqrt(R2rinnerr(i))
                    << "\n";
            }
            std::cout << "\n----------------------------------------\n\n";
        }

        // make second order update (only for response)
        if (ctype == CT_LRCC2 or ctype == CT_LRCCS) {
            double Rtmp = inner(world, info.R_square*residual, V).sum();
            double Rtmp2 = inner(world, info.R_square*GV, GV).sum();
            const double Rdelta = (0.5 * Rtmp / Rtmp2);
            if (info.parameters.debug() and (world.rank() == 0)) std::cout << "omega, second-order update (FYI): " << std::fixed
                << std::setprecision(info.parameters.output_prec() + 2) << omega << ", " << Rdelta << "\n\n";
        }

        // update singles
        singles.omega = omega;
        vector_real_function_3d new_singles = truncate(GV);
        if (info.parameters.kain()) new_singles = solver.update(singles.get_vecfunction(), residual);
        if (info.parameters.debug()) print_size(world, new_singles, "new_singles");
        // if (ctype == CT_LRCCS or ctype == CT_LRCC2 or ctype == CT_ADC2) Nemo::normalize(new_singles, info.R);
        // if (info.parameters.debug()) print_size(world, new_singles, "new_singles normalized");
        reconstruct(world,new_singles);

        for (size_t i = 0; i < GV.size(); i++) {
            singles(i + info.parameters.freeze()).function = copy(new_singles[i]);
        }

        // update regularization terms of the doubles
        // -- not necessary here as the regularization terms are updated in the macrotasks..
        // if (ctype==CT_CC2) update_reg_residues_gs(world, singles,gs_doubles, info);
        // else if(ctype==CT_LRCC2) update_reg_residues_ex(world, singles2,singles,ex_doubles, info);

        if (world.rank()==0) CCPotentials::print_convergence(singles.name(0),rmsresidual,
                                                             maxresidual,omega-old_omega,iter);
        converged = (R2vector_error < info.parameters.dconv_3D());

        // time.info();
        if (converged) break;
        if (ctype == CT_LRCCS) break; // for CCS just one iteration to check convergence
    } // end of iterations

    if (world.rank()==0) print_header2("Singles iterations ended");
    time_all.info();
    print_size(world, singles.get_vecfunction(), "singles after iteration");

    // Assign the overall changes
    bool no_change = true;
    for (auto& tmp : singles.functions) {
        const double change = (tmp.second.function - old_singles(tmp.first).function).norm2();
        tmp.second.current_error = change;
        if (change > info.parameters.dconv_3D()) no_change = false;
    }

    if (info.parameters.debug() and (world.rank() == 0)) {
        std::cout << "Change in Singles functions after all the Microiterations" << std::endl;
        for (auto& tmp : singles.functions)
            std::cout << "Change of " << tmp.second.name() << " = " << tmp.second.current_error << std::endl;
    }

    // update regularization terms of the doubles
    // -- not necessary here as the regularization terms are updated in the macrotasks..
    // if (ctype == CT_CC2) update_reg_residues_gs(world, singles, gs_doubles, info);
    // else if (ctype == CT_LRCC2) update_reg_residues_ex(world, singles2, singles, ex_doubles, info);

    if (no_change) output("Change of Singles was below  = " + std::to_string(info.parameters.dconv_3D()) + "!");

    return no_change;
}

CC_vecfunction
CC2::initialize_singles(World& world, const CalcType& ctype, const FuncType type, const int ex) {

    MADNESS_CHECK_THROW((type==PARTICLE) or (ex>=0), "Invalid type/excitation combination in initialize_singles");
    std::string fname=singles_name(ctype,type,ex);
    if (world.rank()==0) print("initializing singles from file:",fname);
    CC_vecfunction singles(type);
    try {
        singles=CC_vecfunction::load_restartdata(world,fname);
        if (world.rank()==0) print(" .. singles found on file");
    } catch (...) {
        if (world.rank()==0) print(" .. singles not found on file");
    }

    return singles;
}

CC_vecfunction
CC2::initialize_singles_to_zero(World& world, const CalcType& ctype, const FuncType type, const Info& info) {
    CC_vecfunction singles(type);
    for (size_t i = info.parameters.freeze(); i < info.mo_ket.size(); i++) {
        real_function_3d tmpi = real_factory_3d(world);
        CCFunction<double,3> single_i(tmpi, i, type);
        singles.insert(i,single_i);
    }
    return singles;
}


bool
CC2::initialize_pairs(Pairs<CCPair>& pairs, const CCState ftype, const CalcType ctype,
            const CC_vecfunction& gs_singles, const CC_vecfunction& ex_singles,
            const size_t excitation, const Info& info) const {

    MADNESS_ASSERT(gs_singles.type == PARTICLE);
    MADNESS_ASSERT(ex_singles.type == RESPONSE);
    MADNESS_ASSERT(pairs.empty());

    auto builder=CCPairBuilder(world,info);
    builder.set_gs_singles(gs_singles).set_ex_singles(ex_singles).set_ctype(ctype);

    std::string name1 = CCPair(0, 0, ftype, ctype).name();
    if (world.rank()==0) print("initializing doubles",assign_name(ftype), " --- reading from file(s)",name1);

    for (size_t i = parameters.freeze(); i < CCOPS.mo_ket().size(); i++) {
        for (size_t j = i; j < CCOPS.mo_ket().size(); j++) {
            CCPair pair=builder.make_bare_pair_from_file(i,j);

//            std::string name = CCPair(i, j, ftype, ctype).name();
//            real_function_6d utmp;
//            const bool found = CCOPS.load_function(utmp, name, info.parameters.debug());
//            if (found) restarted = true; // if a single pair was found then the calculation is not from scratch
//            real_function_6d const_part;
//            CCOPS.load_function(const_part, name + "_const", info.parameters.debug());
//            CCPair tmp;
//            if (ctype==CT_MP2)
//                tmp=CCPotentials::make_pair_mp2(world, utmp, i, j, info, false);
//            if (ctype==CT_CC2)
//                tmp=CCPotentials::make_pair_cc2(world, utmp, tau, i, j, info, false);
//            if (ctype==CT_LRCC2 or ctype==CT_ADC2 or ctype==CT_CISPD)
//                tmp=CCPotentials::make_pair_lrcc2(world, ctype, utmp, tau, x, i, j, info, false);
//            tmp.constant_part = const_part;
            pairs.insert(i, j, pair);
        }
    }
    // it's a restarted calculation if all pairs are found on disk
    bool restarted=std::all_of(pairs.allpairs.begin(),pairs.allpairs.end(),
        [](const auto& datum){return datum.second.functions.size()>0 && datum.second.function().is_initialized();});

    return restarted;
}

} /* namespace madness */
