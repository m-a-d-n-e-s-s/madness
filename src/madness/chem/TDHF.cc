/*
 * TDHF.cc
 *
 *  Created on: Aug 11, 2016
 *      Author: kottmanj
 */

#include<madness/chem/TDHF.h>
#include<madness/chem/localizer.h>
#include<madness/chem/oep.h>
#include<madness/world/test_utilities.h>
#include<madness/chem/write_test_input.h>
#include"madness/mra/commandlineparser.h"


namespace madness {

// KAIN allocator for vectorfunctions
struct TDHF_allocator {
    World &world;
    const int noct;

    /// @param[in]	world	the world
    /// @param[in]	nnoct	the number of functions in a given vector
    /// @todo validate doxygen on `nnoct`
    TDHF_allocator(World &world, const int nnoct) : world(world), noct(nnoct) {}

    vector_real_function_3d operator()() {
        return zero_functions<double, 3>(world, noct);
    }

    TDHF_allocator operator=(const TDHF_allocator &other) {
        TDHF_allocator tmp(world, other.noct);
        return tmp;
    }

};


/// ctor with world and parameters, constructs Nemo object on-the-fly, and delegates further
///
///
TDHF::TDHF(World &world, const TDHFParameters &parameters, std::shared_ptr<const Nemo> reference)
        : world(world),
            reference_(reference),
            parameters(parameters),
            g12(),
            mo_ket_(),
            mo_bra_(),
            Q(),
            msg(world) {
    this->parameters.set_derived_values(get_calc());
    initialize();

}


///  ctor with command line parser, constructs SCF and Nemo objects on-the-fly, and delegates further

/// \param world    the world
/// \param parser   the parser
TDHF::TDHF(World &world, const commandlineparser &parser) : TDHF(world,parser,std::make_shared<const Nemo>(world,parser)) {}

///  ctor with command line parser, constructs SCF and Nemo objects on-the-fly, and delegates further

/// \param world    the world
/// \param parser   the parser
TDHF::TDHF(World &world, const commandlineparser &parser, std::shared_ptr<const Nemo> nemo)
        : world(world),
          reference_(nemo),
          parameters(world, parser),
          g12(),
          mo_ket_(),
          mo_bra_(),
          Q(),
          msg(world) {

    if (parameters.do_oep()) {
        std::shared_ptr<OEP> oep(new OEP(world, parser));
        set_reference(oep);
    }
    parameters.set_derived_values(get_calc());
    initialize();
}

void TDHF::initialize() {

    msg.debug = parameters.debug();
    g12=std::make_shared<CCConvolutionOperator<double,3>>(world, OpType::OT_G12, parameters.get_ccc_parameters(get_calcparam().lo()));

    const double old_thresh = FunctionDefaults<3>::get_thresh();
    if (old_thresh > parameters.thresh() * 0.1 and old_thresh > 1.e-5) {
        msg.warning("Threshold of Reference might be too loose |  Response thresh="
                    + std::to_string(parameters.thresh()) + " and Reference thresh=" + std::to_string(old_thresh)
                    + ". Be careful, reference should be tight");
    }
    symmetry_projector = get_nemo()->get_symmetry_projector();
    // do not normalize the x vectors individually!
    symmetry_projector.set_lindep(1.e-2).set_orthonormalize_irreps(false).set_verbosity(0);
    check_consistency();

}

void TDHF::print_frozen_orbitals() const {

    MolecularOrbitals<double, 3> dummy_mo(get_calc()->amo, get_calc()->aeps);
    if (parameters.print_level() > 2) dummy_mo.print_frozen_orbitals(parameters.freeze());
}

MolecularOrbitals<double,3> TDHF::enforce_core_valence_separation(const Tensor<double>& fmat) const {


    if (get_calcparam().localize_method()=="canonical") {
        return MolecularOrbitals<double,3>(get_calc()->amo, get_calc()->aeps);
    }


    // check block structure of the fock matrix
    print("initial fock matrix");
    print(fmat);

    auto nemo=get_reference();
    auto calc=get_calc();
    Localizer localizer(world,calc->aobasis,calc->molecule,calc->ao);
    localizer.set_enforce_core_valence_separation(true).set_method(get_calcparam().localize_method());
    localizer.set_metric(nemo->R);

    MolecularOrbitals<double, 3> mos(get_calc()->amo, get_calc()->aeps);
    mos.recompute_localize_sets();
    auto mo_new= localizer.separate_core_valence(mos, fmat);

    return mo_new;

    mos=mo_new;
    Fock<double,3> fock(world,get_reference().get());
    auto fmat1=fock(make_bra(mos.get_mos()),mos.get_mos());
    bool good=Localizer::check_core_valence_separation(fmat1,mos.get_localize_sets(),true);
    print("block-diagonalized fock matrix");
    print(fmat1);
    if (not good) {
        print("fock matrix not block-diagonal");
        print(fmat1);
        MADNESS_CHECK(0);
    }
    return mo_new;

}

/// compute non-trivial prerequisites for the calculation
void TDHF::prepare_calculation() {
    MADNESS_CHECK_THROW(get_nemo()->check_converged(get_calc()->molecule.get_all_coords()),
        "reference calculation for TDHF is not converged");

    // enforce core-valence separation if localized
    Tensor<double> fmat;
    if (get_calcparam().do_localize()) {
        auto fock_operator=*(get_reference()->make_fock_operator());
        fmat=fock_operator(make_bra(get_calc()->amo),get_calc()->amo);
        auto mos=enforce_core_valence_separation(fmat);
        get_nemo()->get_calc()->amo=mos.get_mos();
        get_nemo()->get_calc()->aeps=mos.get_eps();
    } else {
        std::size_t nmo=get_calc()->aeps.size();
        fmat=Tensor<double>(nmo,nmo);
        for (size_t i=0; i<nmo; ++i) fmat(i,i)= get_calc()->aeps(i);
    }

    std::size_t nfrozen=Localizer::determine_frozen_orbitals(fmat);
    // this won't override user's will
    parameters.set_derived_value<long>("freeze", long(nfrozen));
    print_frozen_orbitals();

    mo_ket_ = make_mo_ket(get_calc()->amo);
    mo_bra_ = make_mo_bra(get_calc()->amo);
    Q = QProjector( mo_bra_.get_vecfunction(), mo_ket_.get_vecfunction());

    if (not parameters.no_compute()) {

        bool need_hf = (get_calc()->xc.hf_exchange_coefficient() != 0.0) and (parameters.do_oep() == false);
        if (need_hf) {
            msg.subsection("Computing Exchange Intermediate");
            CCTimer timer(world, "Computing ExIm");
            g12->update_elements(mo_bra_, mo_ket_);
            timer.info();
        } else msg.output("No Exchange Intermediate Computed\n");

        msg.output("Orbital Energies of Reference");
        const Tensor<double> eps = get_calc()->aeps;
        msg << eps << "\n";

    }
    if (get_calcparam().do_localize()) {
        Fock<double,3> F(world,get_reference().get());
        F_occ = F(get_active_mo_bra(), get_active_mo_ket());
//        for (size_t i = 0; i < get_active_mo_ket().size(); ++i) {
//            msg << std::scientific << std::setprecision(10);
//            msg << "F(" << i << "," << i << ")=" << F_occ(i, i) << "\n";
//            if (std::fabs(get_orbital_energy(i + parameters.freeze()) - F_occ(i, i)) > 1.e-5) {
//                msg << "eps(" << i << ")=" << get_orbital_energy(i) << " | diff="
//                    << get_orbital_energy(i + parameters.freeze()) - F_occ(i, i) << "\n";
//            }
//        }
    } else {
        F_occ = Tensor<double>(get_active_mo_bra().size(), get_active_mo_ket().size());
        F_occ *= 0.0;
        for (size_t i = 0; i < get_active_mo_ket().size(); ++i) {
            F_occ(i, i) = get_orbital_energy(i + parameters.freeze());
        }
    }
}

/// plot planes and cubes
void TDHF::plot(const vector_real_function_3d &vf, const std::string &name) const {
    if (parameters.plot()) {
        CCTimer timer(world, "plot planes and cubes");
        madness::plot(vf, name, get_calc()->molecule.cubefile_header());
        timer.print();
    }
}

/// sort the xfunctions according to their excitation energy and name the excitation energies accordingly
std::vector<CC_vecfunction> TDHF::sort_xfunctions(std::vector<CC_vecfunction> x) const {
    std::sort(x.begin(), x.end());
    return x;
}

/// print information
void TDHF::print_xfunctions(const std::vector<CC_vecfunction>& f, const std::string message) const {
    int counter=0;
    if (world.rank()==0 and (not message.empty())) print(message);
    for (const auto &x:f) {
        const double mem = get_size(world, x.get_vecfunction());

        msg << "ex. vector "
            << counter++ << " | "
            << std::setw(4) << x.irrep << " | "
            << x.omega << " (ex. energy) | "
            << std::setprecision(3)
            << x.current_error << " (error) | "
            << x.delta << " (Edelta) | "
            << mem << " (Gbyte)"
            << std::setprecision(msg.output_prec)
            << "\n";

    }
}

void TDHF::initialize(std::vector<CC_vecfunction> &start) const {

    msg.subsection("Calculate Guess");
    std::vector<CC_vecfunction> guess;
    guess = make_guess_from_initial_diagonalization();
    // combine guess and start vectors
    for (const auto &tmp:start) guess.push_back(tmp);
    std::vector<vector_real_function_3d> empty;
    //orthonormalize(guess,empty);


    // failsafe (works in most cases)
    if (guess.size() < size_t(parameters.guess_excitations())) {
        std::string message = ("WARNING: You demanded: " + std::to_string(parameters.guess_excitations())
                               + " Guess vectors, but your demanded guess has only " + std::to_string(guess.size())
                               +
                               "vectors. So we will not iterate the first vectors and then do the same guess again ... this might be unstable").c_str();
        msg.output(message);
        if (parameters.guess_maxiter() == 0) {
            msg.output(
                    "In this case you demanded guess_maxiter=0 (which is also the default), so this can not work!");
            MADNESS_EXCEPTION("Faulty combinations of parameters given", 1);
        }
        iterate_cis_guess_vectors(guess);
        initialize(guess);
    }

    //sort guess (according to excitation energies)
    std::sort(guess.begin(), guess.end());
    //truncate the guess
    std::vector<CC_vecfunction> guess_vectors;
    for (size_t i = 0; i < parameters.guess_excitations(); i++) guess_vectors.push_back(guess[i]);
    // this is the return value
    start = guess_vectors;
}

void TDHF::symmetrize(std::vector<CC_vecfunction> &v) const {

    // The irreps of the x vector elements are given by
    //    Gamma(orbital) x Gamma(x element) = Gamma (excitation)
    // Using the inverse element of Gamma(orbital) we get
    //    Gamma(x element) = Gamma(excitation) x Gamma^{-1}(orbital)
    // since all groups are Abelian the inverse element is the element itself
    //    Gamma(x element) = Gamma(excitation) x Gamma(orbital)

    if (get_nemo()->do_symmetry()) {
        std::vector<std::string> irreps, orbital_irreps;
        symmetry_projector(get_active_mo_ket(), orbital_irreps);

        // loop over all excitations
        for (auto &f : v) {

            // determine the irreps of the x vector elements
            std::vector<std::string> xirreps(f.get_vecfunction().size());
            MADNESS_ASSERT(symmetry_projector.get_table().is_abelian);
            for (size_t i = 0; i < xirreps.size(); ++i) {
                xirreps[i] = symmetry_projector.reduce(f.irrep, orbital_irreps[i])[0];
            }

            vector_real_function_3d tmp = symmetry_projector.project_on_irreps(f.get_vecfunction(), xirreps);
            f.set_functions(tmp, RESPONSE, parameters.freeze());
        }
    }
}

/// @param[in/out] CC_vecfunction
/// on input the guess functions (if empty or not enough the a guess will be generated)
/// on output the solution
std::vector<CC_vecfunction> TDHF::solve_cis() const {
    if (world.rank()==0) print_header2("computing CIS excitations");
    std::vector<CC_vecfunction> ccs;
    // look for restart options
    if (parameters.restart()=="iterate" or parameters.restart()=="no_compute") {
        auto excitations_list=parameters.excitations();
        if (excitations_list.empty()) {
            for (size_t i=0; i<parameters.nexcitations(); ++i) excitations_list.push_back(i);
        }
        for (auto ex : excitations_list) {
            std::string filename= filename_for_roots(ex);
            try {
                auto root=CC_vecfunction::load_restartdata(world,filename);
                if (world.rank()==0) print("found excitation #",ex,"on disk: ",filename);
                ccs.push_back(root);
            } catch (...) {
                if (world.rank()==0) print("could not find excitation #", ex, ", ", filename);
            }
        }
        if (ccs.size()>0) {
            print_xfunctions(ccs, "initial roots from disk");
            if (world.rank()==0) print("\nsorting roots according to energy");
            ccs=sort_xfunctions(ccs);

            if (parameters.restart()=="no_compute") {
                converged_roots=ccs;
            } else {
                // sort ccs in converged and non-converged sets
                std::vector<CC_vecfunction> other_roots;
                double dconv = parameters.dconv();
                double econv = parameters.econv();
                std::partition_copy(std::begin(ccs), std::end(ccs),
                                    std::back_inserter(converged_roots), std::back_inserter(other_roots),
                                    [&econv, &dconv](const CC_vecfunction &root) {
                                        return root.is_converged(econv, dconv);
                                    });
                std::swap(other_roots, ccs);
                if (parameters.print_level()>0) {
                    print_xfunctions(converged_roots, "\nconverged roots");
                    print_xfunctions(ccs, "\nnon-converged roots");
                }
            }
        }
    }

    bool skip_solve=(parameters.restart()=="no_compute") or (converged_roots.size()>=parameters.nexcitations());

    if (skip_solve) {
        if (world.rank()==0) {
            print("skipping the solution of the CIS equations");
            if (parameters.restart()=="no_compute") print(" -> no_compute is set");
            if (converged_roots.size()>=parameters.nexcitations())
                print(" -> number of converged excitations from disk is sufficient:", converged_roots.size());
        }
    } else {
        for (size_t macrocycle = 0; macrocycle < 1; ++macrocycle) {
            //msg.section("CIS Macroiteration " + std::to_string(macrocycle));
            ccs = solve_cis(ccs);
            if (converged_roots.size() >= size_t(parameters.nexcitations())) break;
        }
    }

    if (world.rank()==0) print_header2("end computing CIS excitations");
    return converged_roots;
}

std::vector<CC_vecfunction> TDHF::solve_cis(std::vector<CC_vecfunction> &start) const {
    msg.section("SOLVING CIS EQUATIONS");

    mo_ket_.plot("MOS_");

    CCTimer time(world, "TDHF/CIS");
    // decide if a guess calculation is needed
    bool need_guess = false;
    if (start.size() < size_t(parameters.guess_excitations())) need_guess = true;
    std::vector<CC_vecfunction> guess_vectors;
    if (need_guess) {
        initialize(start);
        guess_vectors = start;
    } else guess_vectors = start;

    symmetrize(guess_vectors);

    msg.output("====Guess-Vectors=====");
    print_xfunctions(guess_vectors, "");

    std::vector<CC_vecfunction> final_vectors;
    msg.subsection("Iterate Guess Vectors");
    {
        // do guess iterations
        iterate_cis_guess_vectors(guess_vectors);
        // sort according to excitation energies
        std::sort(guess_vectors.begin(), guess_vectors.end());
        // save
        for (size_t i = 0; i < guess_vectors.size(); i++) guess_vectors[i].save_restartdata(world, filename_for_roots(i));
        // prepare final iterations
        for (size_t i = 0; i < guess_vectors.size(); i++) {
            if (i < size_t(parameters.iterating_excitations())) final_vectors.push_back(guess_vectors[i]);
            else guess_roots.push_back(guess_vectors[i]);
        }
        // sort guess_roots backwards in order to feed them into the cycle easily with pop_back
        std::sort(guess_roots.rbegin(), guess_roots.rend());

    }

    msg.subsection("Iterate Final Vectors");

    print_xfunctions(final_vectors, "Vectors in Iteration Cycle");
    print_xfunctions(guess_roots, "Remaining Guess Vectors");

    iterate_cis_final_vectors(final_vectors);
    msg.section("CIS CALCULATIONS ENDED");
    std::sort(converged_roots.begin(), converged_roots.end());
    std::sort(final_vectors.begin(), final_vectors.end());

    // information
    print_xfunctions(converged_roots, "\n\nCONVERGED ROOTS\n");
    print_xfunctions(final_vectors, "\n\nRest\n");
    time.info();

    // plot final vectors
    for (size_t i = 0; i < final_vectors.size(); i++) final_vectors[i].plot(std::to_string(i) + "_converged_cis");

    //give back as much as demanded
    std::vector<CC_vecfunction> result = converged_roots;
    for (const auto &x:final_vectors) {
        result.push_back(x);
        if (result.size() == size_t(parameters.nexcitations())) break;
    }
    msg << "writing final functions to disk...\n";
    result = sort_xfunctions(result);
    for (size_t i = 0; i < result.size(); i++) result[i].save_restartdata(world,filename_for_roots(i));


    // print out all warnings
    msg.print_warnings();
    return result;
}

void TDHF::solve_tdhf(std::vector<CC_vecfunction> &x) const {
    msg.section("SOLVING TDHF EQUATIONS");
    MADNESS_EXCEPTION("TDHF NOT IMPLEMENTED", 1);
}

bool TDHF::iterate_cis_guess_vectors(std::vector<CC_vecfunction> &x) const {
    std::vector<CC_vecfunction> dummy;
    return iterate_vectors(x, dummy, false, parameters.guess_dconv(), parameters.guess_econv(),
                           parameters.guess_maxiter(), false);
}

bool TDHF::iterate_cis_final_vectors(std::vector<CC_vecfunction> &x) const {
    std::vector<CC_vecfunction> dummy;
    return iterate_vectors(x, dummy, false, parameters.dconv(), parameters.econv(), parameters.maxiter(),
                           parameters.kain_subspace() > 0);
}

bool TDHF::iterate_vectors(std::vector<CC_vecfunction> &x, const std::vector<CC_vecfunction> &y,
                           const bool iterate_y, const double dconv, const double econv, const double iter_max,
                           const bool kain) const {
    if (iterate_y) MADNESS_ASSERT(x.size() == y.size());

    msg.subsection("Iterating excitation vectors with the following parameters");
    msg << "dconv    = " << dconv << "\n"
        << "econv    = " << econv << "\n"
        << "iter_max = " << iter_max << "\n"
        << "kain     = " << kain << "\n";

    // set up the kain solvers ... if needed or not
    std::vector<std::shared_ptr<XNonlinearSolver<vector_real_function_3d, double, TDHF_allocator> > > solvers(
            x.size());
    // initialize solvers
    if (kain) {
        for (size_t i = 0; i < x.size(); i++) {
            solvers[i] = std::make_shared<XNonlinearSolver<vector_real_function_3d, double, TDHF_allocator> >(
                    TDHF_allocator(world, x[i].size()), true);
            solvers[i]->set_maxsub(parameters.kain_subspace());
        }
    }

    bool converged = true;

    symmetrize(x);

    // get potentials (if demanded)
    std::vector<vector_real_function_3d> V;

    // for TDHF the potential is always stored
    if (parameters.store_potential() and y.empty()) V = make_potentials(x);
    else if (not y.empty()) V = make_tdhf_potentials(x, y);


    for (size_t iter = 0; iter < iter_max; iter++) {
        CCTimer timer(world, "iteration " + std::to_string(iter));
        const bool this_is_a_guess_iteration = not(kain); // better readable code

        // if this is the first iteration, the potentials are not calculated yet
        if (iter == 0) orthonormalize(x, V);

        // if we solve for the y functions we need to switch the sign of omega
        if (iterate_y) {
            for (auto &tmp:x) tmp.omega *= -1.0;
        }
        // apply Greens Functions
        for (size_t i = 0; i < x.size(); i++) if (iterate_y and x[i].omega >= 0.0) x[i].omega = -1.0 * y[i].omega;
        std::vector<vector_real_function_3d> residuals = apply_G(x, V);


        // check convergence

        converged = true;
        // error (vector norm)
        std::vector<double> errors;
        // largest individual error
        std::vector<std::pair<int, double> > largest_errors;
        for (size_t i = 0; i < x.size(); i++) {
            double av_error = norm2(world, residuals[i]);
            std::vector<double> ind_err = norm2s(world, residuals[i]);
            auto it = max_element(ind_err.begin(), ind_err.end());
            const double largest_error = (*it);
            errors.push_back(av_error);
            largest_errors.push_back(std::make_pair(std::distance(ind_err.begin(), it), *it));
            // convergece criteria are the individual functions (orhterwise we have a dependence on the number of functions)
            if ((*it) > dconv) converged = false;
            if (fabs(x[i].delta) > econv) converged = false;
            x[i].current_error = (largest_error);
        }
        // update, store old omegas in the deltas vector
        for (size_t i = 0; i < x.size(); i++) {
            if (std::fabs(x[i].current_error) > dconv or std::fabs(x[i].delta) > econv or
                this_is_a_guess_iteration) {
                vector_real_function_3d new_x0;
                if (kain) {
                    new_x0 = solvers[i]->update(x[i].get_vecfunction(), residuals[i], 10.0 * parameters.thresh(), 5.0);
                } else {
                    new_x0 = sub(world, x[i].get_vecfunction(), residuals[i]);
                }

                // Project out the converged roots
                if (converged_roots.size() > 0) {
                    vector_real_function_3d bra_new_x0 = make_bra(new_x0);
                    CCTimer timeP(world, "project out converged roots");
                    for (const auto& it : converged_roots) {
                        const double overlap = inner(world, it.get_vecfunction(), bra_new_x0).sum();
                        new_x0 -= overlap * it.get_vecfunction();
                    }
                    timeP.print();
                }

                vector_real_function_3d Q_new_x0 = Q(new_x0);
                truncate(world, Q_new_x0);
                x[i].set_functions(Q_new_x0, x[i].type, parameters.freeze());

            } else msg.output("Root " + std::to_string(i) + " converged");
        }
        symmetrize(x);

        // remove converged roots
        if (this_is_a_guess_iteration) {
            // not in the guess iteration
        } else {
            std::vector<CC_vecfunction> unconverged_roots;
            std::vector<std::shared_ptr<XNonlinearSolver<vector_real_function_3d, double, TDHF_allocator> > > corresponding_solvers;
            for (size_t i = 0; i < x.size(); ++i) {
                if (x[i].is_converged(econv,dconv)) {
                    // put converged roots into converged roots and replace it by one of the remaining guess functions
                    converged_roots.push_back(x[i]);
                    if (guess_roots.empty()) msg.output("Ran out of guess functions"); // no replacement
                    else {
                        // fill the unconverged roots with new guess roots
                        unconverged_roots.push_back(guess_roots.back());
                        // allocate a new solver for the new guess function
                        auto sol = std::make_shared<XNonlinearSolver<vector_real_function_3d, double, TDHF_allocator> >(
                                TDHF_allocator(world, x[i].size()), true);
                        sol->set_maxsub(parameters.kain_subspace());
                        corresponding_solvers.push_back(sol);
                        // delete the used guess root
                        guess_roots.pop_back();
                    }
                } else {
                    unconverged_roots.push_back(x[i]);
                    corresponding_solvers.push_back(solvers[i]);
                }
            }
            x = unconverged_roots;
            solvers = corresponding_solvers;
            sort_xfunctions(converged_roots);
        }

        // for TDHF the potential is always stored
        if (parameters.store_potential() and y.empty()) V = make_potentials(x);
        else if (not y.empty()) V = make_tdhf_potentials(x, y);
        // orthonormalize
        orthonormalize(x, V);

        // save functions (every 5 iterations)
        for (size_t i = 0; i < x.size(); i++) x[i].save_restartdata(world, filename_for_roots(i));

        // make plots if demanded
        if (parameters.plot()) {
            if (iterate_y)
                for (size_t i = 0; i < y.size(); i++)
                    y[i].plot(std::to_string(i) + "_y_iter_" + std::to_string(iter) + "_");
            else
                for (size_t i = 0; i < x.size(); i++)
                    x[i].plot(std::to_string(i) + "_iter_" + std::to_string(iter) + "_");
        }
        // information
        msg << std::setfill('-') << std::setw(60) << "\n" << std::setfill(' ');
        msg << "Iteration " << iter << ": omega, largest error, delta, size " << "\n";

        print_xfunctions(converged_roots, "===========CONVERGED=ROOTS=====================");
        print_xfunctions(x, "=========NON-CONVERGED=ROOTS====================");

        msg << std::setfill('-') << std::setw(60) << "\n" << std::setfill(' ');

        timer.stop().print();
        if (converged) break;
        if (converged_roots.size() >= size_t(parameters.nexcitations())) break;
        if (x.empty()) break;
    }
    return converged;
}


std::vector<vector_real_function_3d>
TDHF::apply_G(std::vector<CC_vecfunction> &x, std::vector<vector_real_function_3d> &V) const {

    std::string msg1 = "Applying Greens Function to vectors";
    if (V.empty()) msg1 += ", with recalculated Potentials";
    CCTimer time(world, msg1);
    std::vector<vector_real_function_3d> result;
    for (size_t i = 0; i < x.size(); i++) {

        vector_real_function_3d Vi;
        if (V.empty()) Vi = get_tda_potential(x[i]);
        else Vi = V[i];
        double omega = x[i].omega;
        if (x[i].type == RESPONSE) MADNESS_ASSERT(omega > 0.0);
        else
        MADNESS_ASSERT(omega < 0.0);
        if (x[i].type == UNDEFINED and V.empty()) msg.warning("Empty V but x is y state from TDHF");

        CCTimer time_N(world, "add nuclear potential");
        // the potentials still need the nuclear potential
        const Nuclear<double,3> V(world, get_reference().get());
        vector_real_function_3d VNi = V(x[i].get_vecfunction());
        Vi += VNi;
        time_N.info(parameters.debug());

        // scale potential
        scale(world, Vi, -2.0);

        // make bsh operators as pointers in order to apply them in parallel
        std::vector<std::shared_ptr<SeparatedConvolution<double, 3> > > bsh(x[i].size());
        for (size_t p = 0; p < bsh.size(); p++) {
            double eps = get_orbital_energy(p + parameters.freeze()) + omega;
            // if eps is above zero we have an unbound state (or are early in the iteration) however this needs a shift of the potential
            // we shift to -0.05 (same as moldft does, no particular reason)
            if (eps > 0.0) {
                msg.output("potential shift needed for V" + std::to_string(p + parameters.freeze()));
                double shift = eps + 0.05;
                eps = eps - shift;
                Vi[p] -= (-2.0 * shift * x[i].get_vecfunction()[p]);
            }
            if (eps > 0.0) {
                msg.warning("eps is " + std::to_string(eps) + "... should not happen ... setting to zero");
                eps = -1.0 * parameters.thresh();
            }
            MADNESS_ASSERT(not(eps > 0.0));
            bsh[p] = poperatorT(BSHOperatorPtr3D(world, sqrt(-2.0 * eps), get_calcparam().lo(), parameters.thresh()));
        }
        world.gop.fence();

        vector_real_function_3d GV = Q(apply(world, bsh, Vi));

        vector_real_function_3d residual = sub(world, x[i].get_vecfunction(), GV);
        result.push_back(residual);

//         // Calculate Second Order Energy Update
//         const vector_real_function_3d bra_GV = make_bra(GV);
//         {
//             // Inner product of Vpsi and the residual (Vi is scaled to -2.0 --> multiply later with 0.5)
//             double tmp = inner(world, make_bra(residual), Vi).sum();
//             // squared norm of GVpsi (Psi_tilde)
//             double tmp2 = inner(world, make_bra(GV), GV).sum();

//             // Factor 0.5 removes the factor 2 from the scaling before
//             //const double sou = (0.5 * tmp / tmp2);
// //            msg << "FYI: second order update would be: " << sou << " norm after QG is " << tmp2 << "\n";
//         }
        // clear potential
        Vi.clear();
    }
    // clear the potentials
    V.clear();
    time.info();
    return result;
}

std::vector<vector_real_function_3d> TDHF::make_potentials(const std::vector<CC_vecfunction> &x) const {
    CCTimer time(world, "Make Potentials");
    std::vector<vector_real_function_3d> V;
    for (auto &xi:x) {
        if (parameters.debug()) msg << std::setfill('-') << std::setw(60) << "\n" << std::setfill(' ');
        const vector_real_function_3d pot = get_tda_potential(xi);
        V.push_back(pot);
        if (parameters.debug()) msg << std::setfill('-') << std::setw(60) << "\n" << std::setfill(' ');
    }
    time.info();
    MADNESS_ASSERT(V.size() == x.size());
    return V;
}

vector_real_function_3d TDHF::get_tda_potential(const CC_vecfunction &x) const {


    // Occupation numbers
    const Tensor<double> occ = get_calc()->get_aocc();
    // Closed shell full density of the nemo orbitals (without any nuclear cusps)
    const real_function_3d nemo_density = 2.0 * get_reference()->compute_density(mo_ket_.get_vecfunction());
    // Real Alpha density (with nuclear cusps)
    const real_function_3d alpha_density = 0.5 * get_reference()->R_square * nemo_density;


    // Apply Ground State Potential to x-states
    vector_real_function_3d Vpsi1;
    {
        Fock<double,3> fock(world,get_reference().get());
        fock.remove_operator("T");
        fock.remove_operator("V");      // Vnuc is computed elsewhere (why?)
//        print("ground state potential fock operator",fock.info());
        Vpsi1=fock(x.get_vecfunction());
    }

    // Apply the Perturbed Potential to the Active Ground State Orbitals
    vector_real_function_3d Vpsi2;
    {
        // active mo
        const vector_real_function_3d active_mo = get_active_mo_ket();
        const vector_real_function_3d active_bra = get_active_mo_bra();
        // construct perturbed operators
        CCTimer timeJ(world, "pXC");
        Coulomb<double,3> Jp(world, get_nemo().get());
        real_function_3d density_pert = 2.0 * get_nemo()->make_density(occ, active_bra, x.get_vecfunction());
        Jp.potential() = Jp.compute_potential(density_pert);

        vector_real_function_3d XCp = zero_functions<double, 3>(world, get_active_mo_ket().size());
        if (parameters.response_kernel()!="hf") {
            // XC Potential
            const std::string xc_data = parameters.response_kernel();
            const XCOperator<double,3> xc(world, xc_data, not get_calcparam().spin_restricted(), alpha_density,
                                alpha_density);
            // reconstruct the full perturbed density: do not truncate!
            real_function_3d gamma = xc.apply_xc_kernel(density_pert);
            XCp = truncate(gamma * active_mo);
        }

        if (parameters.triplet()) {
            if (norm2(world, XCp) != 0.0) MADNESS_EXCEPTION("Triplets only for CIS", 1);
            Vpsi2 = XCp;
        } else Vpsi2 = Jp(active_mo) + XCp;
        timeJ.info(parameters.debug());
        // Exchange Part
        double hf_coeff = get_calc()->xc.hf_exchange_coefficient();
        bool do_hf = (hf_coeff != 0.0) and (not parameters.do_oep());
        if (do_hf) {
            CCTimer timeK(world, "pK");
            vector_real_function_3d Kp;
            // summation over all active indices
            for (const auto& itmp:x.functions) {
                const size_t i = itmp.first;
                real_function_3d Ki = real_factory_3d(world);
                for (const auto& ktmp:x.functions) {
                    const size_t k = ktmp.first;
                    Ki += ((*g12)(mo_bra_(k), mo_ket_(i)) * x(k).function).truncate();
                }
                Kp.push_back(Ki);
            }
            scale(world, Kp, hf_coeff);
            Vpsi2 = sub(world, Vpsi2, Kp);
            timeK.info(parameters.debug());
            truncate(world, Vpsi2);
        }


        // compute the solvent (PCM) contribution to the kernel
        if (get_nemo()->do_pcm()) {
            CCTimer timepcm(world, "pcm:ex");
            const real_function_3d vpcm = get_nemo()->get_pcm().compute_pcm_potential(Jp.potential(), true);
            if (parameters.plot() or parameters.debug()) plot_plane(world, vpcm, "vpcm_ex");
            const vector_real_function_3d pcm_orbitals = vpcm * active_mo;
            timepcm.info(parameters.debug());
            Vpsi2 = add(world, Vpsi2, pcm_orbitals);
        }

        truncate(world, Vpsi2);
    }
    // whole tda potential
    vector_real_function_3d Vpsi = Vpsi1 + Q(Vpsi2);
    // if the ground state is localized add the coupling terms
    // canonical: -ei|xi> (part of greens function)
    // local:  -fik|xk> (fii|xi> part of greens functions, rest needs to be added)

    if (get_calcparam().do_localize()) {
        const vector_real_function_3d vx = x.get_vecfunction();
        vector_real_function_3d fock_coupling = madness::transform(world, vx, F_occ);
        // subtract the diagonal terms
        for (size_t i = 0; i < fock_coupling.size(); ++i) {
            fock_coupling[i] = (fock_coupling[i] - F_occ(i, i) * vx[i]);
        }
        Vpsi -= fock_coupling;
    }


    truncate(world, Vpsi);



    // debug output
    if (parameters.debug() or parameters.plot()) {
        plot_plane(world, Vpsi, "Vpsi");
        plot_plane(world, Vpsi1, "Vpsi1");
        plot_plane(world, Vpsi2, "Vpsi2");
    }

    return Vpsi;

}

std::vector<vector_real_function_3d>
TDHF::make_tdhf_potentials(std::vector<CC_vecfunction> &x, const std::vector<CC_vecfunction> &y) const {
    MADNESS_EXCEPTION("NOT IMPLEMENTED", 1);
}


void TDHF::orthonormalize(std::vector<CC_vecfunction> &x, std::vector<vector_real_function_3d> &V) const {
    if (x.empty()) return;
    CCTimer time(world, "Orthonormalization");

    // make the overlap matrix
    Tensor<double> S = make_overlap_matrix(x);
    if (parameters.debug()) msg << "The Overlap Matrix\n " << S << "\n";

    //make the Hamilton matrix for the vectorfunctions
    Tensor<double> F = make_perturbed_fock_matrix(x, V);

    // Diagonalize the F Matrix
    Tensor<double> U, evals;
    Tensor<double> dummy(x.size());
    U = get_calc()->get_fock_transformation(world, S, F, evals, dummy, 2.0 * parameters.thresh());

    if (parameters.debug()) {
        msg << "Perturbed Fock-Matrix Eigenvalues\n";
        for (int x = 0; x < evals.size(); ++x) msg << evals(x) << "\n";
    }

    // Transform the states
    x = transform(x, U);

    // Transform the potentials (if V is empty nothing will happen)
    V = transform(V, U);

    // assign new energies and  get energy differences
    for (size_t i = 0; i < x.size(); i++) {
        const double old_omega = x[i].omega;
        const double new_omega = evals(i);
        const double delta = new_omega - old_omega;
        x[i].omega = new_omega;
        x[i].delta = delta;
    }

    time.info();
}

std::vector<CC_vecfunction>
TDHF::transform(const std::vector<CC_vecfunction> &x, const madness::Tensor<double> U) const {
    std::vector<CC_vecfunction> transformed;
    for (size_t k = 0; k < x.size(); k++) {
        vector_real_function_3d new_x = zero_functions_compressed<double, 3>(world, x[k].size());
        compress(world, x[k].get_vecfunction());
        for (size_t l = 0; l < x.size(); l++) {
            // gaxpy(alpha,a,beta,b) -> a[i]=alpha*a[i] + beta*b[i], since there is no += for vectorfunctions implemented
            gaxpy(world, 1.0, new_x, U(l, k), x[l].get_vecfunction());
        }
        CC_vecfunction tmp(x[k]);
        tmp.set_functions(new_x, tmp.type, parameters.freeze());
        transformed.push_back(tmp);
    }
    MADNESS_ASSERT(transformed.size() == x.size());
    return transformed;

}

Tensor<double> TDHF::make_overlap_matrix(const std::vector<CC_vecfunction> &x) const {
    CCTimer time(world, "Make Overlap Matrix");
    Tensor<double> S(x.size(), x.size());
    for (size_t k = 0; k < x.size(); k++) {
        const vector_real_function_3d kbra = make_bra(x[k]);
        for (size_t l = 0; l < x.size(); l++) {
            S(l, k) = inner(world, kbra, x[l].get_vecfunction()).sum();
        }
    }
    time.info();
    if (parameters.debug())msg << std::fixed << std::setprecision(5) << "\nOverlap Matrix\n" << S << "\n";
    return S;
}

Tensor<double> TDHF::make_perturbed_fock_matrix(const std::vector<CC_vecfunction> &x,
                                                const std::vector<vector_real_function_3d> &V) const {
    // Make formated timings
    CCTimer timeF(world, "Matrix: F");
    CCTimer timeT(world, "Matrix: T+Vn");
    CCTimer timeV(world, "Matrix: V");
    CCTimer timeR(world, "Matrix: e");

    // bra elements of x
    std::vector<vector_real_function_3d> xbra;

    {
        CCTimer time_bra(world, "Make bra elements");
        for (size_t k = 0; k < x.size(); k++) {
            const vector_real_function_3d xbrak = make_bra(x[k]);
            xbra.push_back(xbrak);
        }
        MADNESS_ASSERT(xbra.size() == x.size());
        time_bra.info(parameters.debug());
    }

    timeF.start();
    Tensor<double> F(x.size(), x.size());
    {
        Tensor<double> T(x.size(), x.size());
        {
            timeT.start();
            // gradient operator
            std::vector<std::shared_ptr<real_derivative_3d> > D = gradient_operator<double, 3>(world);

            real_function_3d Vnuc = get_calc()->potentialmanager->vnuclear();
            if (not Vnuc.is_initialized()) {
                msg.output("Compute Nuclear Potential");
                get_calc()->potentialmanager->make_nuclear_potential(world);
                Vnuc = get_calc()->potentialmanager->vnuclear();
            }

            const real_function_3d R = get_reference()->ncf->function();
            std::vector<vector_real_function_3d> Rx(x.size(), zero_functions<double, 3>(world, x.front().size()));
            CCTimer timeR(world, "make Rx");
            for (size_t k = 0; k < x.size(); k++) {
                Rx[k] = mul(world, R, x[k].get_vecfunction(), false);
            }
            world.gop.fence();
            timeR.info(parameters.debug());
            std::vector<vector_real_function_3d> dx(x.size(), zero_functions<double, 3>(world, x.front().size()));
            std::vector<vector_real_function_3d> dy(x.size(), zero_functions<double, 3>(world, x.front().size()));
            std::vector<vector_real_function_3d> dz(x.size(), zero_functions<double, 3>(world, x.front().size()));
            CCTimer timeD(world, "make Grad(Rx)");
            for (size_t k = 0; k < x.size(); k++) {
                dx[k] = apply(world, *(D[0]), Rx[k], false);
                dy[k] = apply(world, *(D[1]), Rx[k], false);
                dz[k] = apply(world, *(D[2]), Rx[k], false);
            }
            world.gop.fence();
            timeD.info(parameters.debug());

            CCTimer time_mat(world, "T+V Mat");
            for (size_t k = 0; k < x.size(); k++) {
                const vector_real_function_3d Vxk = mul(world, Vnuc, x[k].get_vecfunction());
                for (size_t l = 0; l < x.size(); l++) {
                    T(l, k) = inner(world, xbra[l], Vxk).sum();
                    T(l, k) += 0.5 * inner(world, dx[l], dx[k]).sum();
                    T(l, k) += 0.5 * inner(world, dy[l], dy[k]).sum();
                    T(l, k) += 0.5 * inner(world, dz[l], dz[k]).sum();
                }

            }
            time_mat.info(parameters.debug());
            timeT.stop();
        }
        Tensor<double> MV(x.size(), x.size());
        {
            timeV.start();
            bool recompute_V = V.empty();
            for (size_t k = 0; k < x.size(); k++) {
                vector_real_function_3d Vk;
                //if(recompute_V) Vk = CCOPS.get_CIS_potential(x[k]);
                if (recompute_V) {
                    msg.output("Recompute V");
                    Vk = get_tda_potential(x[k]);
                } else Vk = V[k];
                for (size_t l = 0; l < x.size(); l++) {
                    MV(l, k) = inner(world, xbra[l], Vk).sum();
                }
            }
            timeV.stop();
        }

        // now set the fock matrix together: F(l,k) = T(l,k) + MV(l,k) - eps(l,k)
        // with eps(l,k) = \sum_i eps_i <xl_i|xk_i>
        // first get active eps
        timeR.start();
        std::vector<double> eps;
        for (size_t i = 0; i < mo_ket_.size(); i++) eps.push_back(get_orbital_energy(i));
        std::vector<double> active_eps;
        for (size_t i = parameters.freeze(); i < eps.size(); i++) active_eps.push_back(eps[i]);
        for (size_t k = 0; k < x.size(); k++) {
            for (size_t l = 0; l < x.size(); l++) {
                Tensor<double> xlk = inner(world, xbra[l], x[k].get_vecfunction());
                MADNESS_ASSERT(size_t(xlk.size()) == active_eps.size());
                double eps_part = 0.0;
                for (size_t i = 0; i < active_eps.size(); i++) eps_part += xlk(i) * active_eps[i];
                F(l, k) = T(l, k) + MV(l, k) - eps_part;
            }
        }
        timeR.stop();
        if (parameters.debug())msg << std::fixed << std::setprecision(5) << "\n(T+V) Matrix\n" << T << "\n";
        if (parameters.debug())msg << std::fixed << std::setprecision(5) << "\nPotential Matrix\n" << MV << "\n";
        if (parameters.debug())
            msg << std::fixed << std::setprecision(5) << "\nPerturbed Fock Matrix\n" << F << "\n";
    }
    timeF.stop();
    if (parameters.debug())msg << std::fixed << std::setprecision(5) << "\nPerturbed Fock Matrix\n" << F << "\n";
    //formated timings output
    timeT.print();
    timeV.print();
    timeR.print();
    timeF.print();

    // symmetryze
    F = 0.5 * (F + transpose<double>(F));
    if (parameters.debug())
        msg << std::fixed << std::setprecision(5) << "\nSymmetrized Perturbed Fock Matrix\n" << F << "\n";

    return F;
}

vector_real_function_3d TDHF::make_virtuals() const {
    CCTimer time(world, "make virtuals");
    // create virtuals
    vector_real_function_3d virtuals;
    if (parameters.guess_excitation_operators() == "external") {
        madness::load_function(world, virtuals, "mybasis");
        //virtuals=Q(virtuals);
        for (auto &x : virtuals) {
            const double norm = sqrt(inner(make_bra(x), x));
            x.scale(1.0 / norm);
        }
    } else if (parameters.guess_excitation_operators() == "scf") {
        // use the ao basis set from the scf calculations as virtuals (like projected aos)
        virtuals = (get_calc()->ao);
        for (auto &x : virtuals) {
            const double norm = sqrt(inner(make_bra(x), x));
            x.scale(1.0 / norm);
        }
    } else {
        // create the seeds
        vector_real_function_3d xmo;
//        for (size_t i = 0; i < parameters.guess_occ_to_virt(); ++i)
        for (size_t i = 0; i < get_active_mo_ket().size(); ++i)
            xmo.push_back(get_active_mo_ket()[get_active_mo_ket().size() - 1 - i]);

        bool use_trigo = true;
//		if(parameters.generalkeyval.find("polynomial_exops")!=parameters.generalkeyval.end())
//			use_trigo = (std::stoi(parameters.generalkeyval.find("polynomial_exops")->second)==0);
        virtuals = apply_excitation_operators(xmo, use_trigo);

    }

    if (parameters.guess_cm() > 0.0) {
        // add center of mass diffuse functions
        const double factor = parameters.guess_cm();
        const double width = (-1.0 * factor / (get_orbital_energy(mo_ket_.size() - 1)));
        msg.subsection("adding center of mass functions with exponent homo/c and c=" + std::to_string(factor));
        msg.output("width=" + std::to_string(width));

        Tensor<double> cm = get_calc()->molecule.center_of_mass();
        msg << "center of mass is " << cm << "\n";
        guessfactory::PolynomialFunctor px("x 1.0", width, cm);
        guessfactory::PolynomialFunctor py("y 1.0", width, cm);
        guessfactory::PolynomialFunctor pz("z 1.0", width, cm);
        guessfactory::GaussFunctor s(width, cm);
        real_function_3d vpx = real_factory_3d(world).functor(px);
        real_function_3d vpy = real_factory_3d(world).functor(py);
        real_function_3d vpz = real_factory_3d(world).functor(pz);
        real_function_3d vs = real_factory_3d(world).functor(s);
        virtuals.push_back(vpx);
        virtuals.push_back(vpy);
        virtuals.push_back(vpz);
        virtuals.push_back(vs);
    }
    if (world.rank() == 0)
        msg << virtuals.size() << " virtuals\n";

    plot(virtuals, "virtuals");
    CCTimer timerQ(world, "virt=Qvirt");
    virtuals = Q(virtuals);
    timerQ.print();
    plot(virtuals, "Qvirtuals");

    if (parameters.debug()) {
        for (const auto &x : virtuals)
            x.print_size("virtual");
    }
    world.gop.fence();
    time.print();
    return virtuals;
}

vector_real_function_3d
TDHF::apply_excitation_operators(const vector_real_function_3d &seed, const bool &use_trigo) const {
    //const int nvirt = seed.size() * ((2 * order * 2 * order * 2 * order) - 1);
    //	msg.subsection("creating a set of " + std::to_string(nvirt) + " virtuals by multiplying functions with plane waves");
    // compute the centers of the seed functions
    CCTimer time_centers(world, "compute centers");
    std::vector<coord_3d> centers = guessfactory::compute_centroids(seed);
    time_centers.print();


    // prepare the list of excitation operators and copied seeds
    CCTimer time_init_exop(world, "initialize excitation operators");
    std::vector<std::pair<vector_real_function_3d, std::string> > exlist;
    {
        std::vector<std::string> exop_strings = parameters.exops();
        if (parameters.guess_excitation_operators() != "custom")
            exop_strings = (guessfactory::make_predefined_exop_strings(parameters.guess_excitation_operators()));
        for (const auto& ex: exop_strings) {
            vector_real_function_3d cseed = copy(world, seed, false);
            exlist.push_back(std::make_pair(cseed, ex));
        }
    }
    world.gop.fence();
    time_init_exop.print();
    msg << "will create " << exlist.size() * seed.size() << " virtuals, from " << seed.size() << " seeds and "
        << exlist.size() << " excitation operators" << " \n";

    // create the virtuals by unary operations: multiply excitation operators with seeds
    CCTimer time_create_virtuals(world, "create virtuals");
    vector_real_function_3d virtuals;
    for (auto it:exlist) {
        if (use_trigo)
            virtuals = append(virtuals, guessfactory::apply_trigonometric_exop(it.first, it.second, centers,
                                                                               false));
        else virtuals = append(virtuals, guessfactory::apply_polynomial_exop(it.first, it.second, centers, false));
    }
    world.gop.fence();
    time_create_virtuals.print();

    return virtuals;

}

/// make the initial guess by explicitly diagonalizing a CIS matrix with virtuals from the make_virtuals routine
vector<CC_vecfunction> TDHF::make_guess_from_initial_diagonalization() const {
    CCTimer time(world, "make_guess_from_initial_diagonalization");
    //convenience
    const int nact = get_active_mo_ket().size();

    // create virtuals
    vector_real_function_3d virtuals = make_virtuals();
    // remove linear dependencies
    CCTimer time_ortho(world, "rr-cholesky orthonormalization");
    const size_t spre = virtuals.size();
    Tensor<double> S = matrix_inner(world, make_bra(virtuals), virtuals);
    virtuals = orthonormalize_rrcd(virtuals,S, 1.e-6);
    if (virtuals.size() != spre)
        msg << "removed " << spre - virtuals.size() << " virtuals due to linear dependencies \n";
    time_ortho.print();

    // project virtuals onto irreps
    std::vector<std::string> virtual_irreps;
    projector_irrep proj = projector_irrep(symmetry_projector).set_verbosity(0).set_lindep(1.e-1);
    virtuals = proj(virtuals, get_reference()->R_square, virtual_irreps);

    // re-orthonormalize virtuals (seems to improve numerical stability in the canonicalization step below)
    get_reference()->orthonormalize(virtuals,get_reference()->R);
    virtuals = proj(virtuals, get_reference()->R_square, virtual_irreps);
    if (world.rank() == 0) print("final number of virtuals", virtuals.size());

    // determine the symmetry of the occupied and virtual orbitals
    std::vector<std::string> orbital_irreps;
    proj(get_active_mo_ket(), get_reference()->R_square, orbital_irreps);

    // canonicalize virtuals
    Tensor<double> veps;    // orbital energies of the virtuals
    virtuals = canonicalize(virtuals, veps);

    // make sure the virtual orbital energies are higher than the occupied orbtials
    double vmin = veps.min();
    double omax = get_calc()->aeps.max();
    if (vmin < omax) {
        veps += (omax - vmin);
        if (world.rank() == 0) print("shifting guess virtual energies by ", omax - vmin + 1.e-1);
    }

    // compute the CIS matrix
    Tensor<double> MCIS = make_cis_matrix(virtuals, veps);

    // zero out all non-contributing irreps
    if (parameters.irrep() != "all") {
        int I = 0;
        for (auto oirrep1 : orbital_irreps) {
            for (auto virrep1 : virtual_irreps) {
                if (not(proj.reduce(oirrep1, virrep1)[0] == parameters.irrep())) {
                    MCIS(I, _) = 0.0;
                    MCIS(_, I) = 0.0;
                }
                I++;
            }
        }
    }

    // initialize the guess functions
    if (world.rank() == 0 && size_t(MCIS.dim(0)) < parameters.guess_excitations()) {
        msg.warning(std::to_string(parameters.guess_excitations())
                    + " guess vectors where demanded, but with the given options only "
                    + std::to_string(MCIS.dim(0)) + " can be created\n");
    }

    const int nvirt = virtuals.size();
    auto get_com_idx = [nvirt](int i, int a) {
        return i * nvirt + a;
    };
    auto get_vir_idx = [nvirt](int I) {
        return I % nvirt;
    };
    auto get_occ_idx = [nvirt](int I) {
        return I / nvirt;
    };


    // find all contributing ia/jb elements for the requested irrep -- needs to be improved..
    std::vector<int> II;
    for (size_t i = 0; i < orbital_irreps.size(); ++i) {
        for (size_t a = 0; a < virtual_irreps.size(); ++a) {
            if (proj.reduce(orbital_irreps[i], virtual_irreps[a])[0] == parameters.irrep())
                II.push_back(get_com_idx(i, a));
        }
    }
    std::vector<CC_vecfunction> xfunctions;
    for (size_t x = 0; x < size_t(MCIS.dim(0)); ++x) {
        if (x >= parameters.guess_excitations())
            break;

        CC_vecfunction init(zero_functions<double, 3>(world, nact), RESPONSE, parameters.freeze());
        xfunctions.push_back(init);
    }
    {
        Tensor<double> U, evals;
        CCTimer time_diag(world, "cis-matrix diagonalization");
        syev(MCIS, U, evals);
        Localizer::undo_degenerate_rotations(U,evals,FunctionDefaults<3>::get_thresh());

        time_diag.print();

        if (parameters.debug()) {
            msg << "Initial Diagonalization of CIS Matrix:\n Lowest Eigenvalues are \n";
            for (int x = 0; x < std::max(std::min(evals.size(), 10l), long(parameters.guess_excitations())); ++x)
                msg << evals(x) << "\n";
        }

        CCTimer time_assemble(world, "assemble guess vectors");
        // make x functions from amplitudes and virtuals
        // todo: probably not optimal for highly parallel applications
        int iexcitation = 0;
        for (int I = 0; I < MCIS.dim(0); ++I) {
            if (size_t(iexcitation) >= parameters.guess_excitations()) break;

            //const int a = get_vir_idx(I);
            //const int i = get_occ_idx(I);
            if (evals(I) < 0.0 && world.rank() == 0)
                msg.warning("NEGATIVE EIGENVALUE IN INITIAL DIAGONALIZATION: CHECK YOUR REFERENCE!\n");

            if (evals(I) < 1.e-5) {
                if (parameters.debug()) msg << "skipping root " << evals(I) << " \n";
                continue;
            }

            xfunctions[iexcitation].omega = evals(I);

            for (int J = 0; J < MCIS.dim(1); ++J) {

                const int b = get_vir_idx(J);
                const int j = get_occ_idx(J);
                const double xjb = U(J, I);
                xfunctions[iexcitation].get_vecfunction()[j] += xjb * virtuals[b];
            }
            std::vector<std::string> x_irreps;
            vector_real_function_3d tmp = proj(xfunctions[iexcitation].get_vecfunction(), get_reference()->R_square, x_irreps);

            iexcitation++;
        }
        time_assemble.print();
        CCTimer time_truncate(world, "truncate guess");
        for (auto &x : xfunctions) {
            vector_real_function_3d tmp = x.get_vecfunction();
            truncate(world, tmp, parameters.thresh());
            // should be truncated by shallow copy anyways ... but just to be sure
            x.set_functions(tmp, x.type, parameters.freeze());
        }
        time_truncate.print();

        CCTimer time_symmetrize(world, "symmetrize guess");
        for (auto &x : xfunctions) {
            std::vector<std::string> x_irreps;
            vector_real_function_3d tmp = proj(x.get_vecfunction(), get_reference()->R_square, x_irreps);
            try {
                x.set_functions(tmp, RESPONSE, parameters.freeze());
                for (size_t i = 0; i < x_irreps.size(); ++i) {
                    std::string reduced = symmetry_projector.reduce(x_irreps[i], orbital_irreps[i])[0];
                    if (not((reduced == x.irrep) or (reduced == "null") or (x.irrep == "null"))) {
                        print("reduced, irrep", reduced, x.irrep);
                        MADNESS_EXCEPTION("inconsistent symmetry in x vector\n\n", 0);
                    }
                    if (reduced != "null") x.irrep = reduced;
                }
            } catch (...) {
                if (world.rank()==0) {
                    print_header1("error");
                    print("if you are here, there is likely a problem with the CIS guess, where degenerate states");
                    print("of different irreps may mix and mess up the symmetry determination");
                    print("This needs to be fixed, for now just run the same calculation again with fingers crossed");
                    print("\nx_irreps");
                    print(x_irreps);
                    print("\n");
                }
                MADNESS_EXCEPTION("inconsistent symmetry in x vector\n\n", 0);
            }
        }
        time_symmetrize.print();

    }
    if (parameters.debug()) {
        Tensor<double> S = make_overlap_matrix(xfunctions);

        msg << "Overlap matrix of guess:\n" << S << "\n";
    }
    print_xfunctions(xfunctions, "");
    msg << "created " << xfunctions.size() << " guess vectors\n";

    time.print();
    return xfunctions;
}

/// canonicalize a set of orbitals (here the virtuals for the guess)

/// @param[out]	veps	orbital energies of the virtuals
vector_real_function_3d TDHF::canonicalize(const vector_real_function_3d &v, Tensor<double> &veps) const {
    CCTimer time(world, "canonicalize");

    Fock<double,3> F(world, get_reference().get());
    const vector_real_function_3d vbra = make_bra(v);
    Tensor<double> Fmat = F(vbra, v);

    Tensor<double> S = matrix_inner(world, vbra, v);
    Tensor<double> occ(v.size());
    occ = 1.0;
    if (parameters.debug())
        msg << "Canonicalize: Fock Matrix\n" << Fmat(Slice(0, std::min(10, int(v.size())) - 1),
                                                     Slice(0, std::min(10, int(v.size())) - 1));
    if (parameters.debug())
        msg << "Canonicalize: Overlap Matrix\n" << S(Slice(0, std::min(10, int(v.size())) - 1),
                                                     Slice(0, std::min(10, int(v.size())) - 1));
    Tensor<double> U = get_calc()->get_fock_transformation(world, S, Fmat, veps, occ,1.e-3);

    vector_real_function_3d result = madness::transform(world, v, U);
    time.print();
    return result;
}

/// compute the CIS matrix for a given set of virtuals
Tensor<double> TDHF::make_cis_matrix(const vector_real_function_3d& virtuals,
                                     const Tensor<double> &veps) const {

    CCTimer time_cis(world, "make CIS matrix");

    // the cis matrix is indexed by ij and ab
    // we will use the combined indixes from ia and jb named I and J
    // in order to not be confused we use the following helper functions
    const int nocc = get_active_mo_ket().size();
    // determines for which orbitals (couting from the HOMO downwards) the off-diagonal elements will be computed
    // this simplifies the guess
//    int active_guess_orbitals = parameters.guess_active_orbitals();
    int active_guess_orbitals = get_active_mo_bra().size();
    const int nvirt = virtuals.size();
    auto get_com_idx = [nvirt](int i, int a) { return i * nvirt + a; };
    auto get_vir_idx = [nvirt](int I) { return I % nvirt; };
    auto get_occ_idx = [nvirt](int I) { return I / nvirt; };
    auto delta = [](int x, int y) { if (x == y) return 1; else return 0; };


    const int dim = (virtuals.size() * nocc);
    msg << "CIS-Matrix for guess calculation will be of size " << dim << "x" << dim << "\n";
    // the number of the matrix where elements which are not determined by orbital energies and the fock matrix are computed (controlled over active_guess_orbitals parameter)
    const int dim2 = (virtuals.size() * active_guess_orbitals);
    if (dim2 < dim)
        msg << "Effective size through neglect of some orbitals will be: " << dim2 << "x" << dim2 << "\n";
    const int start_ij = nocc - active_guess_orbitals;
    Tensor<double> MCIS(dim, dim);

    // make CIS matrix
    // first do the "diagonal" entries
    if (get_calcparam().do_localize()) {

        // make bra elements
        const vector_real_function_3d virtuals_bra = make_bra(virtuals);

        // make Fock Matrix of virtuals for diagonal elements
        Fock<double,3> F(world, get_reference().get());
        Tensor<double> Fmat = F(virtuals_bra, virtuals);

        if (parameters.debug()) {
            const int dim = std::min(10, int(virtuals.size()));
            msg << "Debug Part of Virtual Fock Matrix\n" << Fmat(Slice(0, dim - 1), Slice(0, dim - 1)) << "\n";

            Tensor<double> S = matrix_inner(world, virtuals_bra, virtuals);
            msg << "Debug Overlap of virtuals\n" << S(Slice(0, dim - 1), Slice(0, dim - 1)) << "\n";
        }

        Tensor<double> Focc = F(get_active_mo_bra(), get_active_mo_ket());
        for (int I = 0; I < dim; ++I) {
            const int a = get_vir_idx(I);
            const int i = get_occ_idx(I);
            for (int J = 0; J < dim; ++J) {
                const int b = get_vir_idx(J);
                const int j = get_occ_idx(J);
                MCIS(I, J) = Fmat(a, b) * delta(i, j) - Focc(i, j) * delta(a, b);
            }
        }

    } else {    // canonical case with virtual orbital energies only
        for (int I = 0; I < dim; ++I) {
            const int a = get_vir_idx(I);
            const int i = get_occ_idx(I);
            MCIS(I, I) = veps(a) - get_orbital_energy(i + parameters.freeze());
        }
    }

    if (not parameters.guess_diag()) {
        const vector_real_function_3d virtuals_bra = make_bra(virtuals);

        bool is_dft=(get_calc()->xc.is_dft()) or (parameters.do_oep());
        if (not is_dft) {
            int I = -1; // combined index from i and a, start is -1 so that initial value is 0 (not so important anymore since I dont use ++I)
            for (size_t i = start_ij; i < get_active_mo_ket().size(); ++i) {
                const real_function_3d brai = get_active_mo_bra()[i];
                const vector_real_function_3d igv = (*g12)(brai * virtuals);
                for (size_t a = 0; a < virtuals.size(); ++a) {
                    I = get_com_idx(i, a);
                    int J = -1;
                    for (size_t j = start_ij; j < get_active_mo_ket().size(); ++j) {
                        const real_function_3d braj = get_active_mo_bra()[j];
                        for (size_t b = 0; b < virtuals.size(); ++b) {
                            J = get_com_idx(j, b);
                            if (J <= I) {
                                const real_function_3d igj = (*g12)(mo_bra_(i + parameters.freeze()), mo_ket_(j +
                                                                                                              parameters.freeze())); // use exchange intermediate
                                const double rIJ = 2.0 * inner(braj * virtuals[b], igv[a]) -
                                                   inner(virtuals_bra[a] * virtuals[b], igj);
                                MCIS(J, I) += rIJ;
                                MCIS(I, J) += rIJ;
                            }
                        }
                    }
                }
            }
        } else { // is_dft

            std::string xc_data= (parameters.do_oep()) ? "lda_x" : get_calcparam().xc();
            real_function_3d alpha_density=get_reference()->compute_density(mo_bra_.get_vecfunction());

            const XCOperator<double,3> xc(world, xc_data, not get_calcparam().spin_restricted(),
                                              alpha_density, alpha_density);
//                real_function_3d gamma = xc.apply_xc_kernel(density_pert);
//                vector_real_function_3d XCp = truncate(gamma * active_mo);
//                Vpsi2 = Vpsi2 + XCp;
            std::vector<real_function_3d> ia_func(dim);
            for (int ia=0; ia<dim; ++ia) {
                int i=get_occ_idx(ia);
                int a=get_vir_idx(ia);
                ia_func[ia]=get_active_mo_bra()[i]*virtuals[a];
            }
            for (int ia=0; ia<dim; ++ia) {

                // make g | ia )  (Mullikan notation)
                const real_function_3d igv = (*g12)(ia_func[ia]);
                for (int jb=ia; jb<dim; ++jb) {

                    const double coulomb=2.0*inner(ia_func[jb],igv);
                    const double exchange=inner(ia_func[ia],xc.apply_xc_kernel(ia_func[jb]));
                    MCIS(ia,jb)+= coulomb;
                    MCIS(jb,ia)+= exchange;
                }
            }
        }
    }
    if (parameters.debug()) {
        int sdim = std::min(int(MCIS.dim(0)), 10);
        msg << "Part of the CIS Matrix:\n" << MCIS(Slice(dim - sdim, -1), Slice(dim - sdim, -1)) << "\n";
        if (parameters.debug()) msg << "Debug: Full CIS Matrix:\n" << MCIS << "\n";
    }

    // test if symmetric
    if (parameters.debug()) {
        const double symm_norm = (MCIS - transpose(MCIS)).normf();
        msg << "Hermiticity of CIS Matrix:\n" << "||MCIS-transpose(MCIS)||=" << symm_norm << "\n";

        if (symm_norm > 1.e-4) {
            long sliced_dim = 8;
            if (8 > MCIS.dim(0l))
                sliced_dim = MCIS.dim(0l);

            msg << "first " << sliced_dim << "x" << sliced_dim << " block of MCIS Matrix\n"
                << MCIS(_, Slice(sliced_dim - 1, sliced_dim - 1));
        }
    }
    time_cis.info();

    return MCIS;
}


/// compute the oscillator strength in the length representation

/// the oscillator strength is given by
/// \f[
/// f = 2/3 * \omega |<x | \vec \mu | i >| ^2 * 2
/// \f]
/// where \f$ x \f$ is the excited state, and \f$ i \f$ is the ground state
/// @param[in]  root    a converged root
double TDHF::oscillator_strength_length(const CC_vecfunction &x) const {
    Tensor<double> mu_if(3);
    for (int idim = 0; idim < 3; idim++) {
        real_function_3d ri = real_factory_3d(world).functor(guessfactory::PolynomialFunctor(idim));
        vector_real_function_3d amo_times_x = ri * get_active_mo_bra();
        Tensor<double> a = inner(world, amo_times_x, x.get_vecfunction());
        mu_if(idim) = a.sum();
    }
    const double f = 2.0 / 3.0 * x.omega * mu_if.sumsq() * 2.0;
    return f;
}

/// compute the oscillator strength in the velocity representation

/// the oscillator strength is given by
/// \f[
/// f = 2/(3 * \omega) |<x | \vec p | i >| ^2 * 2
/// \f]
/// where \f$ x \f$ is the excited state, and \f$ i \f$ is the ground state
/// @param[in]  root    a converged root
double TDHF::oscillator_strength_velocity(const CC_vecfunction &x) const {
    Tensor<double> p_if(3);
    // compute the derivatives of the MOs in all 3 directions
    const vector_real_function_3d Rroot = get_reference()->R * x.get_vecfunction();
    const vector_real_function_3d Rnemo = get_reference()->R * get_active_mo_ket();

    for (int idim = 0; idim < 3; idim++) {
        real_derivative_3d D = free_space_derivative<double, 3>(world, idim);
        vector_real_function_3d Damo = apply(world, D, Rnemo);
        Tensor<double> a = inner(world, Damo, Rroot);
        p_if(idim) = a.sum();
    }
    const double f = 2.0 / (3.0 * x.omega) * p_if.sumsq() * 2.0;
    return f;
}


/// analyze the root: oscillator strength and contributions from occ
void TDHF::analyze(const std::vector<CC_vecfunction> &x) const {

    const size_t noct = get_active_mo_ket().size();
    nlohmann::json j;
    int counter=0;

    for (const CC_vecfunction &root : x) {

        const vector_real_function_3d Rroot =
                get_reference()->R * root.get_vecfunction(); // reintroduce the nuclear correlation factor
        std::vector<double> norms = norm2s(world, Rroot);

        // compute the oscillator strengths and dominant contributions
        double osl = this->oscillator_strength_length(root);
        double osv = this->oscillator_strength_velocity(root);

        msg << std::scientific << std::setprecision(10) << std::setw(20);

        msg << "excitation energy for root "
            << std::fixed << std::setprecision(1) << counter++ << ": "
            << std::fixed << std::setprecision(10) << root.omega << " Eh         "
            << root.omega * constants::hartree_electron_volt_relationship << " eV   "
            << "irrep: " << root.irrep << "\n";
        msg << std::scientific;
        if (world.rank() == 0)print("  oscillator strength (length)    ", osl);
        if (world.rank() == 0)print("  oscillator strength (velocity)  ", osv);
        // print out the most important amplitudes
        if (world.rank() == 0)print("  dominant contributions ");

        for (std::size_t p = 0; p < noct; ++p) {
            const double amplitude = norms[p] * norms[p];
            if ((amplitude > 0.1)) {
                msg << "    norm(x_" << p + parameters.freeze() << ") **2  ";
                std::cout.width(10);
                std::cout.precision(6);
                msg << amplitude << "\n";
            }
        }
        if (world.rank() == 0) print(" ");
        j.push_back(root);
        j.back()["oscillator_strength_length"]=osl;
        j.back()["oscillator_strength_velocity"]=osv;
    }

    nlohmann::json j1;
    j1["cis_excitations"]=j;
    if (world.rank()==0) update_schema(get_calc()->param.prefix()+".calc_info", j1);


    // compute the transition densities
    const vector_real_function_3d bra_oct = get_active_mo_bra();
    for (std::size_t i = 0; i < x.size(); ++i) {
        const vector_real_function_3d root = x[i].get_vecfunction();
        const real_function_3d td = dot(world, root, bra_oct);
        const double trace = td.trace();
        if (world.rank() == 0) print("trace over transition density", i, trace);
        save(td, "transition_density_" + std::to_string(i));
    }
}

/// auto assigns all parameters which where not explicitly given and which depend on other parameters of the reference calculation
//void TDHF::TDHFParameters::complete_with_defaults(const std::shared_ptr<SCF>& scf) {
void TDHFParameters::set_derived_values(const std::shared_ptr<SCF> &scf) {
    //  double thresh=FunctionDefaults<3>::get_thresh();

    set_derived_value("econv", scf->param.econv());
    set_derived_value("dconv", sqrt(get<double>("econv")) * 0.1);
    set_derived_value("guess_econv", econv() * 10.0);
    set_derived_value("guess_dconv", dconv() * 10.0);

    set_derived_value("iterating_excitations", std::min(nexcitations(), std::size_t(4)));
    set_derived_value("guess_excitations", std::min(nexcitations() + iterating_excitations(), 2 * nexcitations()));

    set_derived_value("response_kernel", scf->param.xc());
    if (do_oep()) set_derived_value("response_kernel",std::string("lda_x"));
//    set_derived_value("guess_occ_to_virt", scf->amo.size() - freeze());
//    set_derived_value("guess_active_orbitals", scf->amo.size() - freeze());
}

/// check consistency of the input parameters
void TDHF::check_consistency() const {

    // check if the requested irrep is present in the computational point group
    const std::vector<std::string> irreps = get_symmetry_projector().get_table().mullikan_;
    if (find(irreps.begin(), irreps.end(), parameters.irrep()) == irreps.end()
        and (parameters.irrep() != "all")) {
        print("irrep ", parameters.irrep(), " is not contained in point group ",
              get_symmetry_projector().get_table().schoenflies_, "\n\n");
        MADNESS_EXCEPTION("\ninconsistent input parameters\n\n", 1);
    }
}

int TDHF::test(World &world, commandlineparser& parser) {

    // write the input file
    parser.set_keyval("input","bla");
    std::string input = parser.value("input");
    std::ofstream of(input.c_str());
    CalculationParameters scfparam;
    scfparam.set_user_defined_value("econv", 1.e-3);
    TDHFParameters cisparam;
    write_test_input::write_to_test_input("dft", &scfparam, of);
    write_test_input::write_to_test_input("response", &cisparam, of);
    write_test_input::write_molecule_to_test_input("lih", of);
    of.close();


    // construct TDHF class
    TDHF tdhf(world, parser);
    tdhf.get_calcparam().print("dft");
    tdhf.parameters.print("response");


    // Compute MRA-CIS

    // solve the CIS equations
    const double time_scf_start = wall_time();
    tdhf.prepare_calculation();
    const double time_scf_end = wall_time();
    if (world.rank() == 0) printf(" at time %.1f\n", wall_time());

    const double time_cis_start = wall_time();
    std::vector<CC_vecfunction> roots = tdhf.solve_cis();
    const double time_cis_end = wall_time();
    if (world.rank() == 0) printf(" at time %.1f\n", wall_time());

    // check result
    int error = 0;
    if (std::abs(roots[0].omega - 1.2893418409e-01) > 1.e-5) error++;

    if (world.rank() == 0) {
        std::cout << std::setfill(' ');
        std::cout << "\n\n\n";
        std::cout << "--------------------------------------------------\n";
        std::cout << "MRA-CIS ended \n";
        std::cout << "--------------------------------------------------\n";
        std::cout << std::setw(25) << "time scf" << " = " << time_scf_end - time_scf_start << "\n";
        std::cout << std::setw(25) << "time cis" << " = " << time_cis_end - time_cis_start << "\n";
        std::cout << "--------------------------------------------------\n";
    }
    tdhf.analyze(roots);
    std::remove(input.c_str());
    if (error == 0) print("\n\nCIS test passed\n\n");
    if (error != 0) print("\n\nCIS test failed\n\n");
    return error;

}

} /* namespace madness */
