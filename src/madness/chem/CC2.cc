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
        Pairs<CCPair> dummypairs;
        bool found_cc2d = initialize_pairs(dummypairs, GROUND_STATE, CT_CC2, cc2singles, CC_vecfunction(RESPONSE), 0, info);
        if (not found_cc2d) need_mp2=true;
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
        // check if singles or/and doubles to restart are there
        initialize_singles(cc2singles, PARTICLE);
        const bool load_doubles = initialize_pairs(cc2pairs, GROUND_STATE, CT_CC2, cc2singles, CC_vecfunction(RESPONSE), 0, info);

        // nothing to restart -> make MP2
        if (not load_doubles) {
            // use mp2 as cc2 guess
            for (auto& tmp:mp2pairs.allpairs) {
                const size_t i = tmp.second.i;
                const size_t j = tmp.second.j;
                cc2pairs(i, j).update_u(tmp.second.function());
            }
        }

        cc2_energy = solve_cc2(cc2singles, cc2pairs, info);
        output_calc_info_schema("cc2",cc2_energy);

        output.section(assign_name(CT_CC2) + " Calculation Ended !");
        if (world.rank() == 0) {
            printf_msg_energy_time("CC2 correlation energy",cc2_energy,wall_time());
            std::cout << std::fixed << std::setprecision(10) << " CC2 Correlation Energy =" << cc2_energy << "\n";
        }
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
            CC_vecfunction dummy = ccs.copy();
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

       std::vector<std::pair<std::string, double> > results;
       std::vector<std::pair<std::string, std::pair<double, double> > > timings;

       auto vccs=solve_ccs();

       std::vector<std::pair<std::string, std::pair<double, double> > > results_ex;
       for (size_t xxx = 0; xxx < vccs.size(); xxx++) {
           const size_t excitation = parameters.excitations()[xxx];
            CCTimer time_ex(world, "LRCC2 Calculation for Excitation " + std::to_string(int(excitation)));
            CC_vecfunction lrcc2_s = vccs[xxx];
            // needed to assign an omega
            const vector_real_function_3d backup = copy(world, lrcc2_s.get_vecfunction());
            CC_vecfunction test(backup, RESPONSE, parameters.freeze());
            iterate_ccs_singles(test, info);
            lrcc2_s.omega = test.omega;
            output("CCS Iteration: Changes are not applied (just omega)!");


            Pairs<CCPair> lrcc2_d;
            bool found_lrcc2d = initialize_pairs(lrcc2_d, EXCITED_STATE, CT_LRCC2, cc2singles, lrcc2_s, excitation, info);

            if (found_lrcc2d) iterate_lrcc2_singles(world, cc2singles, cc2pairs, lrcc2_s, lrcc2_d, info);
            else iterate_ccs_singles(lrcc2_s, info);
            const double omega_cis = lrcc2_s.omega;

            for (size_t iter = 0; iter < parameters.iter_max(); iter++) {
                output.section("Macroiteration " + std::to_string(int(iter)) + " of LRCC2");
                bool dconv = iterate_lrcc2_pairs(world, cc2singles, lrcc2_s, lrcc2_d, info);
                bool sconv = iterate_lrcc2_singles(world, cc2singles, cc2pairs, lrcc2_s, lrcc2_d, info);
                update_reg_residues_ex(world, cc2singles, lrcc2_s, lrcc2_d, info);
                if (dconv and sconv) break;
            }
            const double omega_cc2 = lrcc2_s.omega;
            const std::string msg = "Excitation " + std::to_string(int(excitation));
            results_ex.push_back(std::make_pair(msg, std::make_pair(omega_cis, omega_cc2)));
            timings.push_back(std::make_pair(msg, time_ex.current_time(true)));

        }

        timings.push_back(std::make_pair("Whole LRCC2", time.current_time(true)));
        output.section("LRCC2 Finished");
        output("Ground State Results:");
        for (const auto& res:results) {
            if (world.rank() == 0)
                std::cout << std::fixed << std::setprecision(10)
                          << res.first << "=" << res.second << "\n";
        }
        output("Response Results:");
        for (const auto& res:results_ex) {
            if (world.rank() == 0)
                std::cout << std::fixed << std::setprecision(10)
                          << res.first << ": " << res.second.first << " (CIS)*, " << res.second.second << " (CC2)\n";
        }
        if (world.rank() == 0) std::cout << "*only if CIS vectors where given in the beginning (not for CC2 restart)\n";
        output("\nTimings");
        for (const auto& time:timings) {
            if (world.rank() == 0)
                std::cout << std::scientific << std::setprecision(2)
                          << std::setfill(' ') << std::setw(15) << time.first
                          << ": " << time.second.first << " (Wall), " << time.second.second << " (CPU)" << "\n";
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
std::vector<CC_vecfunction> CC2::solve_ccs() const
{
//    output.section("SOLVE CCS");
//    std::vector<CC_vecfunction> excitations;
//    for (size_t k = 0; k < parameters.excitations().size(); k++) {
//        CC_vecfunction tmp;
//        const bool found = initialize_singles(tmp, RESPONSE, parameters.excitations()[k]);
//        if (found) excitations.push_back(tmp);
//    }
//    tdhf->prepare_calculation();
//    excitations = tdhf->solve_cis(excitations);
    std::vector<CC_vecfunction> excitations=tdhf->get_converged_roots();

    // return only those functions which are demanded
    std::vector<CC_vecfunction> result;
    for (const auto& x:parameters.excitations()) {
        if (excitations.size() - 1 < x) MADNESS_EXCEPTION("Not Enough CIS Vectors to solve for the demanded CC2 vector",
                                                          1);
        result.push_back(excitations[x]);
    }
    return result;
}

double CC2::solve_mp2_coupled(Pairs<CCPair>& doubles, Info& info) {

    if (world.rank()==0) print_header2(" computing the MP1 wave function");
    double total_energy = 0.0;

    // make vector holding CCPairs for partitioner of MacroTask
    std::vector<CCPair> pair_vec=Pairs<CCPair>::pairs2vector(doubles,triangular_map);

    // read constant part from file
    if (parameters.no_compute_mp2_constantpart()) {
        if (world.rank()==0) print("Skipping MP2 constant part calculation");
        for (auto& c : pair_vec) {
            MADNESS_CHECK_THROW(c.constant_part.is_initialized(), "could not find constant part");
            // constant part is zero-order guess for pair.function
            if (not c.function().is_initialized()) c.update_u(c.constant_part);
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

        // transform vector back to Pairs structure
        for (size_t i = 0; i < pair_vec.size(); i++) {
            pair_vec[i].constant_part = result_vec[i];
            // pair_vec[i].functions[0] = CCPairFunction<double,6>(result_vec[i]);
            pair_vec[i].constant_part.truncate().reduce_rank();
            pair_vec[i].constant_part.print_size("constant_part");
            pair_vec[i].function().truncate().reduce_rank();
            save(pair_vec[i].constant_part, pair_vec[i].name() + "_const");
            // save(pair_vec[i].function(), pair_vec[i].name());
            if (pair_vec[i].type == GROUND_STATE) {
                double energy = CCOPS.compute_pair_correlation_energy(world,info,pair_vec[i]);
                if (world.rank()==0) printf("pair energy for pair %zu %zu: %12.8f\n", pair_vec[i].i, pair_vec[i].j, energy);
                total_energy += energy;
            }
        }
        if (world.rank()==0) {
            printf("current decoupled mp2 energy %12.8f\n", total_energy);
            std::cout << std::fixed << std::setprecision(1) << "\nFinished saving pairs and energy calculation at time " << wall_time() << std::endl;
        }
    }


    if (world.rank()==0) print_header3("Starting updating MP2 pairs");


    auto solver= nonlinear_vector_solver<double,6>(world,pair_vec.size());
    solver.set_maxsub(parameters.kain_subspace());
    solver.do_print = (world.rank() == 0);


    for (size_t iter = 0; iter < parameters.iter_max_6D(); iter++) {
        if (world.rank()==0) print_header3("Starting iteration " + std::to_string(int(iter)) + " of MP2");

        // compute the coupling between the pair functions
        Pairs<real_function_6d> coupling=compute_local_coupling(pair_vec);
        auto coupling_vec=Pairs<real_function_6d>::pairs2vector(coupling,triangular_map);
        if (parameters.debug()) print_size(world, coupling_vec, "couplingvector");


        if (world.rank()==0) {
            std::cout << std::fixed << std::setprecision(1) << "\nStart updating pairs part at time " << wall_time() << std::endl;
        }

        MacroTaskIteratePair t;
        MacroTask task1(world, t);
        CC_vecfunction dummy_singles1(PARTICLE);
        const std::size_t maxiter=1;
        auto unew = task1(pair_vec, coupling_vec, dummy_singles1, dummy_singles1, info, maxiter);

        std::vector<real_function_6d> u;
        for (auto p : pair_vec) u.push_back(p.function());
        auto residual=u-unew;

        // some statistics
        auto [rmsrnorm, maxrnorm]=residual_stats(residual);

        // update the pair functions
        if (parameters.kain()) {
            if (world.rank()==0) std::cout << "Update with KAIN" << std::endl;
            // std::vector<real_function_6d> kain_update = copy(world,solver.update(u, u_update));
            std::vector<real_function_6d> kain_update = copy(world,solver.update(u, residual));
            for (size_t i=0; i<pair_vec.size(); ++i) {
                kain_update[i].truncate().reduce_rank();
                kain_update[i].print_size("Kain-Update-Function");
                pair_vec[i].update_u(copy(kain_update[i]));
            }
        } else {
            if (world.rank()==0) std::cout << "Update without KAIN" << std::endl;
            for (size_t i=0; i<pair_vec.size(); ++i) {
                // pair_vec[i].update_u(pair_vec[i].function() - u_update[i]);
                pair_vec[i].update_u(unew[i]);
            }
        }

        // calculate energy and error and update pairs
        double old_energy = total_energy;
        total_energy = 0.0;
        for (size_t i = 0; i < pair_vec.size(); i++) {
            save(pair_vec[i].function(), pair_vec[i].name());
            double energy = CCOPS.compute_pair_correlation_energy(world,info,pair_vec[i]);
            total_energy += energy;
            if (world.rank()==0) printf("pair energy for pair %zu %zu: %12.8f\n", pair_vec[i].i, pair_vec[i].j, energy);
        }

		if (world.rank()==0) {
		    std::cout << "convergence: rms/max residual, energy change "
				<< std::scientific << std::setprecision(1)
				<< rmsrnorm << " " << maxrnorm << " "
                << std::abs(old_energy - total_energy) << std::endl;
                // << std::abs(old_norm - total_norm);
			printf("finished MP2 iteration %2d at time %8.1fs with energy  %12.8f\n",
					int(iter), wall_time(), total_energy);
		}

        bool converged = ((std::abs(old_energy - total_energy) < parameters.econv())
                          and (maxrnorm < parameters.dconv_6D()));

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
Pairs<real_function_6d> CC2::compute_local_coupling(const Pairs<real_function_6d>& pairs) const {

    const int nmo = nemo->get_calc()->amo.size();

    // temporarily make all N^2 pair functions
    typedef std::map<std::pair<int, int>, real_function_6d> pairsT;
    pairsT quadratic;
    for (int k = parameters.freeze(); k < nmo; ++k) {
        for (int l = parameters.freeze(); l < nmo; ++l) {
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
    Tensor<double> fock1 = nemo->compute_fock_matrix(nemo->get_calc()->amo, nemo->get_calc()->aocc);
    for (int k = 0; k < nmo; ++k) {
        if (fock1(k, k) > 0.0) MADNESS_EXCEPTION("positive orbital energies", 1);
        fock1(k, k) = 0.0;
    }

    Pairs<real_function_6d> coupling;
    for (int i = parameters.freeze(); i < nmo; ++i) {
        for (int j = i; j < nmo; ++j) {
            coupling.insert(i, j, real_factory_6d(world).compressed());
        }
    }

    for (int i = parameters.freeze(); i < nmo; ++i) {
        for (int j = i; j < nmo; ++j) {
            for (int k = parameters.freeze(); k < nmo; ++k) {
                if (fock1(k, i) != 0.0) {
                    coupling(i, j).gaxpy(1.0, quadratic[std::make_pair(k, j)], fock1(k, i), false);
                }
            }

            for (int l = parameters.freeze(); l < nmo; ++l) {
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
    print_header2("Solve LRCC2 for Excitation energy " + std::to_string(double(lrcc2_s.omega)));
    MADNESS_ASSERT(lrcc2_s.type == RESPONSE);

    auto triangular_map=PairVectorMap::triangular_map(info.parameters.freeze(),info.mo_ket.size());
    auto pair_vec=Pairs<CCPair>::pairs2vector(lrcc2_d,triangular_map);
    MacroTaskIteratePair t1;
    MacroTask task1(world, t1);
    std::vector<real_function_3d> vdummy_3d;         // dummy vectors
    std::vector<real_function_6d> vdummy_6d;         // dummy vectors
    const std::size_t maxiter=3;
    auto unew = task1(pair_vec, vdummy_6d, cc2_s, lrcc2_s, info, maxiter);

    // get some statistics
    std::vector<Function<double,6>> uold;
    for (const auto & p : pair_vec) uold.push_back(p.function());
    auto residual=uold-unew;
    auto [rmsrnorm, rmsrmax] = residual_stats(residual);

    // update the pair functions
    for (int i=0; i<pair_vec.size(); ++i) pair_vec[i].update_u(unew[i]);
    lrcc2_d=Pairs<CCPair>::vector2pairs(pair_vec,triangular_map);

    return (rmsrnorm<info.parameters.dconv_6D());


//    bool conv = true;
//    for (auto& tmp:lrcc2_d.allpairs) {
//        CCPair& pair = tmp.second;
//        const size_t i = pair.i;
//        const size_t j = pair.j;
//        // check if singles have significantly changed
//        if (lrcc2_s(i).current_error < 0.1 * parameters.thresh_6D() and
//            lrcc2_s(j).current_error < 0.1 * parameters.thresh_6D()) {
//            ;
//            // output("Skipping Pair Iteration, No significant Change in Singles");
//        } else {
//            pair.bsh_eps = CCPotentials::get_epsilon(pair.i, pair.j,info) + lrcc2_s.omega;
//            // update_constant_part_lrcc2(pair, cc2_s, lrcc2_s);
//            pair.constant_part=CCPotentials::make_constant_part_macrotask(world, pair,
//                         cc2_s, lrcc2_s, info);
//            conv = iterate_pair(pair, lrcc2_s);
//        }
//    }
//    return conv;

}


double
CC2::solve_cc2(CC_vecfunction& singles, Pairs<CCPair>& doubles, Info& info) const
{

    output.section("Solving CC2 Ground State");

    MADNESS_ASSERT(singles.type == PARTICLE);
    CCTimer time(world, "CC2 Ground State");

    double omega = CCPotentials::compute_cc2_correlation_energy(world, singles, doubles, info);
    if (world.rank() == 0)
        std::cout << std::fixed << std::setprecision(10) << "Current Correlation Energy = " << omega << "\n";
    CC_vecfunction ex_singles_dummy;

    if (not parameters.no_compute_cc2()) {
        // first singles iteration
        output.section("Initialize Singles to the Doubles");

        // given the doubles, we can solve the singles equations
        iterate_cc2_singles(world, singles, doubles, info);
        // the doubles ansatz depends on the singles and must be updated: |\tau_ij> = |u_ij> + Q12 f12 |t_i t_j>
        update_reg_residues_gs(world, singles, doubles, info);
        omega = CCPotentials::compute_cc2_correlation_energy(world, singles, doubles, info);

        for (size_t iter = 0; iter < parameters.iter_max(); iter++) {
            CCTimer time_miter(world, "Macroiteration " + std::to_string(int(iter)) + " of CC2");
            output.section("Macroiteration " + std::to_string(int(iter)) + " of CC2");

            if (world.rank()==0) print("computing the constant part via macrotasks -- output redirected");
            timer timer1(world);

            std::vector<CCPair> pair_vec=Pairs<CCPair>::pairs2vector(doubles,triangular_map);
            MacroTaskConstantPart t;
            MacroTask task(world, t);
            std::vector<real_function_6d> constant_part_vec = task(pair_vec, singles.get_vecfunction(),
                ex_singles_dummy.get_vecfunction(), info) ;
            for (int i=0; i<pair_vec.size(); ++i) pair_vec[i].constant_part=constant_part_vec[i];

            if (world.rank()==0 and parameters.debug()) {
                for (auto& pair: pair_vec) pair.constant_part.print_size("size of constant part macrotask "+pair.name());
            }

            timer1.tag("computing constant part via macrotasks");


            // compute the coupling between the pair functions
            if (world.rank()==0) print("computing local coupling in the universe");
            Pairs<real_function_6d> coupling=compute_local_coupling(pair_vec);
            auto coupling_vec=Pairs<real_function_6d>::pairs2vector(coupling,triangular_map);
            timer1.tag("computing local coupling");

            if (world.rank()==0) print("update the pair functions via macrotasks -- output redirected");
            MacroTaskIteratePair t1;
            MacroTask task1(world, t1);
            CC_vecfunction dummy_ex_singles;
            std::vector<real_function_3d> vdummy_3d;         // dummy vectors
            const std::size_t maxiter=3;
            auto unew = task1(pair_vec, coupling_vec, singles, dummy_ex_singles,
                info, maxiter);

            std::vector<real_function_6d> u_old;
            for (auto p : pair_vec) u_old.push_back(p.function());
            auto residual=u_old-unew;
            timer1.tag("computing pair function update via macrotasks");

            for (int i=0; i<pair_vec.size(); ++i) pair_vec[i].update_u(unew[i]);
            doubles=Pairs<CCPair>::vector2pairs(pair_vec,triangular_map);

            auto [rmsrnorm,maxrnorm]=residual_stats(residual);
            bool doubles_converged=rmsrnorm<parameters.dconv_6D();

            // check if singles converged
            const bool singles_converged = iterate_cc2_singles(world, singles, doubles, info);

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

		    std::cout << "convergence: rms/max residual, energy change "
				<< std::scientific << std::setprecision(1)
				<< rmsrnorm << " " << maxrnorm << " "
                << std::abs(delta) << std::endl;
                // << std::abs(old_norm - total_norm);
			printf("finished CC2 iteration %2d at time %8.1fs with energy  %12.8f\n",
					int(iter), wall_time(), omega);
            if (doubles_converged and singles_converged and omega_converged) break;

            time_miter.info();
        }
        omega = CCPotentials::compute_cc2_correlation_energy(world, singles, doubles, info);
        output.section("CC2 Iterations Eneded");
    } else {
        output.section("Found no_compute_cc2 Key: Reiterating Singles to check convergence");
        // need the singles potential for the constant part of LRCC2 so we recompute it (also good to check if it is converged)
        bool sconv = iterate_cc2_singles(world, singles, doubles, info);
        if (not sconv) output.warning("Singles not Converged");
    }

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

        std::vector<std::pair<std::string, double> > results;
        std::vector<std::pair<std::string, std::pair<double, double> > > timings;
        std::vector<std::pair<std::string, std::pair<double, double> > > results_ex;

        auto vccs=solve_ccs();
        CCTimer time_ex(world, "LRCC2 Calculation for Excitation " + std::to_string(int(excitation)));
        CC_vecfunction ex_singles = vccs[excitation];
        // needed to assign an omega
        const vector_real_function_3d backup = copy(world, ex_singles.get_vecfunction());
        CC_vecfunction test(backup, RESPONSE, parameters.freeze());
        iterate_ccs_singles(test, info);
        ex_singles.omega = test.omega;
        output("CCS Iteration: Changes are not applied (just omega)!");


        Pairs<CCPair> ex_doubles;
        bool found_lrcc2d = initialize_pairs(ex_doubles, EXCITED_STATE, CT_LRCC2, gs_singles, ex_singles, excitation, info);

        if (found_lrcc2d) iterate_lrcc2_singles(world, gs_singles, gs_doubles, ex_singles, ex_doubles, info);
        else iterate_ccs_singles(ex_singles, info);
        const double omega_cis = ex_singles.omega;

        for (size_t iter = 0; iter < parameters.iter_max(); iter++) {
            output.section("Macroiteration " + std::to_string(int(iter)) + " of LRCC2");
            bool dconv = iterate_lrcc2_pairs(world, gs_singles, ex_singles, ex_doubles, info);
            bool sconv = iterate_lrcc2_singles(world, gs_singles, gs_doubles, ex_singles, ex_doubles, info);
            update_reg_residues_ex(world, gs_singles, ex_singles, ex_doubles, info);
            if (dconv and sconv) break;
        }
        const double omega_cc2 = ex_singles.omega;
        const std::string msg = "Excitation " + std::to_string(int(excitation));
        results_ex.push_back(std::make_pair(msg, std::make_pair(omega_cis, omega_cc2)));
        timings.push_back(std::make_pair(msg, time_ex.current_time(true)));


        timings.push_back(std::make_pair("Whole LRCC2", time.current_time(true)));
        output.section("LRCC2 Finished");
        output("Ground State Results:");
        for (const auto& res:results) {
            if (world.rank() == 0)
                std::cout << std::fixed << std::setprecision(10)
                          << res.first << "=" << res.second << "\n";
        }
        output("Response Results:");
        for (const auto& res:results_ex) {
            if (world.rank() == 0)
                std::cout << std::fixed << std::setprecision(10)
                          << res.first << ": " << res.second.first << " (CIS)*, " << res.second.second << " (CC2)\n";
        }
        if (world.rank() == 0) std::cout << "*only if CIS vectors where given in the beginning (not for CC2 restart)\n";
        output("\nTimings");
        for (const auto& time:timings) {
            if (world.rank() == 0)
                std::cout << std::scientific << std::setprecision(2)
                          << std::setfill(' ') << std::setw(15) << time.first
                          << ": " << time.second.first << " (Wall), " << time.second.second << " (CPU)" << "\n";
        }

    return std::make_tuple(ex_doubles, ex_singles, omega_cc2);

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
    if (pair.type == GROUND_STATE) omega = CCOPS.compute_pair_correlation_energy(world, info,pair, singles);
    if (pair.type == EXCITED_STATE) omega = CCOPS.compute_excited_pair_energy(pair, singles);

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
        if (pair.type == GROUND_STATE) omega_new = CCOPS.compute_pair_correlation_energy(world, info, pair, singles);
        else if (pair.type == EXCITED_STATE) omega_new = CCOPS.compute_excited_pair_energy(pair, singles);
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

bool
CC2::initialize_singles(CC_vecfunction& singles, const FuncType type, const int ex) const {
    MADNESS_ASSERT(singles.size() == 0);
    bool restarted = false;

    std::vector<CCFunction<double,3>> vs;
    for (size_t i = parameters.freeze(); i < CCOPS.mo_ket().size(); i++) {
        CCFunction<double,3> single_i;
        single_i.type = type;
        single_i.i = i;
        std::string name;
        if (ex < 0) name = single_i.name();
        else name = std::to_string(ex) + "_" + single_i.name();
        real_function_3d tmpi = real_factory_3d(world);
        const bool found = CCOPS.load_function<double, 3>(tmpi, name);
        if (found) restarted = true;
        else output("Initialized " + single_i.name() + " of type " + assign_name(type) + " as zero-function");
        single_i.function = copy(tmpi);
        vs.push_back(single_i);
    }

    singles = CC_vecfunction(vs, type);
//    if (type == RESPONSE) singles.excitation = ex;

    return restarted;
}

bool
CC2::initialize_pairs(Pairs<CCPair>& pairs, const CCState ftype, const CalcType ctype, const CC_vecfunction& tau,
                      const CC_vecfunction& x, const size_t excitation, const Info& info) const {
    MADNESS_ASSERT(tau.type == PARTICLE);
    MADNESS_ASSERT(x.type == RESPONSE);
    MADNESS_ASSERT(pairs.empty());
    output("Initialize " + assign_name(ctype) + " Pairs for " + assign_name(ftype));

    bool restarted = false;

    for (size_t i = parameters.freeze(); i < CCOPS.mo_ket().size(); i++) {
        for (size_t j = i; j < CCOPS.mo_ket().size(); j++) {

            std::string name = CCPair(i, j, ftype, ctype).name();
            if (ftype == GROUND_STATE) {
                real_function_6d utmp = real_factory_6d(world);
                const bool found = CCOPS.load_function(utmp, name);
                if (found) restarted = true; // if a single pair was found then the calculation is not from scratch
                real_function_6d const_part;
                CCOPS.load_function(const_part, name + "_const");
                CCPair tmp;
                if (ctype==CT_MP2) tmp=CCPotentials::make_pair_mp2(utmp, i, j, info);
                if (ctype==CT_CC2) tmp=CCPotentials::make_pair_cc2(utmp, tau, i, j, info);
                tmp.constant_part = const_part;
                pairs.insert(i, j, tmp);

            } else if (ftype == EXCITED_STATE) {
                name = std::to_string(int(excitation)) + "_" + name;
                real_function_6d utmp = real_factory_6d(world);
                const bool found = CCOPS.load_function(utmp, name);
                if (found) restarted = true;
                real_function_6d const_part;
                CCOPS.load_function(const_part, name + "_const");
                CCPair tmp = CCOPS.make_pair_ex(utmp, tau, x, i, j, ctype);
//                tmp.excitation = excitation;
                tmp.constant_part = const_part;
                pairs.insert(i, j, tmp);
            } else error("Unknown pairtype");
        }
    }
    return restarted;
}

void CC2::update_reg_residues_gs(World& world, const CC_vecfunction& singles, Pairs<CCPair>& doubles, const Info& info)
{
    CCTimer time(world, "Updated Regularization Residues of the Ground State");
    MADNESS_ASSERT(singles.type == PARTICLE);
    Pairs<CCPair> updated_pairs;
    for (auto& tmp:doubles.allpairs) {
        MADNESS_ASSERT(tmp.second.type == GROUND_STATE);
        CCPair& pair = tmp.second;
        const size_t i = pair.i;
        const size_t j = pair.j;
        // const CCPair updated_pair = CCOPS.make_pair_gs(pair.function(), singles, i, j);
        const CCPair updated_pair = CCPotentials::make_pair_cc2(pair.function(), singles, i, j, info);
        updated_pairs.insert(i, j, updated_pair);
    }
    doubles.swap(updated_pairs);
    time.info();
}

void CC2::update_reg_residues_ex(World& world, const CC_vecfunction& singles,
                                 const CC_vecfunction& response, Pairs<CCPair>& doubles, const Info& info)
{
    CCTimer time(world, "Updated Regularization Residues of the Excited State");
    MADNESS_ASSERT(singles.type == PARTICLE);
    MADNESS_ASSERT(response.type == RESPONSE);
    Pairs<CCPair> updated_pairs;
    for (auto& tmp:doubles.allpairs) {
        MADNESS_ASSERT(tmp.second.type == EXCITED_STATE);
        CCPair& pair = tmp.second;
        // CCPair updated_pair = CCPotentials::make_pair_ex(pair.function(), singles, response, i, j, pair.ctype);
        CCPair updated_pair =
            CCPotentials::make_pair_lrcc2(pair.function(), singles, response, pair.i, pair.j, info);
        updated_pairs.insert(pair.i, pair.j, updated_pair);
    }
    doubles.swap(updated_pairs);
    time.info();
}


} /* namespace madness */
