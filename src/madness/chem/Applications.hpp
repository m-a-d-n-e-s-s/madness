#pragma once

#include <madness/chem/InputWriter.hpp>
#include <madness/chem/ParameterManager.hpp>
#include <madness/chem/PathManager.hpp>
#include <madness/chem/Results.h>

namespace madness {
    enum class NextAction { Ok, ReloadOnly, Restart, Redo };

    // Scoped CWD: changes the current directory to the given one, and restores when
    // the object goes out of scope
    struct ScopedCWD {
        std::filesystem::path old_cwd;

        explicit ScopedCWD(std::filesystem::path const &new_dir) {
            old_cwd = std::filesystem::current_path();
            std::filesystem::current_path(new_dir);
        }

        ~ScopedCWD() { std::filesystem::current_path(old_cwd); }
    };

    class Application {
    public:
        explicit Application(const Params &p) : params_(p) {
        }

        virtual ~Application() = default;


        // run: write all outputs under the given directory
        virtual void run(const std::filesystem::path &workdir) = 0;


        // optional hook to return a JSON fragment of this app's main results
        [[nodiscard]] virtual nlohmann::json results() const = 0;

        virtual void print_parameters(World &world) const = 0;

        /// check if this calculation has a json with results
        [[nodiscard]] virtual bool has_results(std::string filename) const {
            // check if the results file exists
            // return std::filesystem::exists(workdir_ / filename);
            return std::filesystem::exists(filename);
        }

        [[nodiscard]] virtual bool verify_molecule(const nlohmann::json &j) const {
            // check if some key parameters of the calculation match:
            // molecule, box size, nmo_alpha, nmo_beta
            Molecule mol1 = params_.get<Molecule>();
            Molecule mol2;
            mol2.from_json(j["molecule"]);
            if (not(mol1 == mol2)) {
                print("molecule mismatch");
                mol1.print();
                mol2.print();
                return false;
            }
            return true;
        }

        /// read the results from a json file
        [[nodiscard]] virtual nlohmann::json read_results(std::string filename) const {
            if (has_results(filename)) {
                std::ifstream ifs(filename);
                nlohmann::json j;
                ifs >> j;
                ifs.close();
                if (not verify_molecule(j)) {
                    std::string msg = "Results file " + filename + " does not match the parameters of the calculation";
                    print(msg);
                    return nlohmann::json(); // return empty json
                }
                return j;
            } else {
                std::string msg = "Results file " + filename + " does not exist in " + std::filesystem::current_path().
                                  string();
                MADNESS_EXCEPTION(msg.c_str(), 1);
            }
            return nlohmann::json();
        }

    protected:
        const Params params_;
        nlohmann::json results_;
    };

    template<typename Library>
    class SCFApplication : public Application {
    private:
    public:
        using Calc = typename Library::Calc;

        explicit SCFApplication(World &w, const Params &p) : Application(p), world_(w) {
        }

        // Give downstream steps the live calc
        std::shared_ptr<Calc> calc() { return lib_.calc(world_, params_); }
        void set_calc_workdir(const std::filesystem::path &workdir) { calc()->work_dir = workdir; }

        // print parameters
        void print_parameters(World &world) const override {
            if (world.rank() == 0) {
                std::cout << "SCF Parameters:" << std::endl;
            }
            lib_.print_parameters();
        }

        // sets the calc working directory and runs the calculation
        void run(const std::filesystem::path &workdir) override {
            // 1) set up a namedspaced directory for this run
            std::string label = Library::label();
            PathManager pm(workdir, label.c_str());
            pm.create(); {
                world_.gop.fence();
                ScopedCWD scwd(pm.dir());
                if (world_.rank() == 0) {
                    std::cout << "Running SCF in " << pm.dir() << std::endl;
                }
                // 2) define the "checkpoint" file
                auto ckpt = label + ".calc_info.json";
                SCFResultsTuple empty_results;
                nlohmann::json j;
                if (has_results(ckpt)) {
                    j = read_results(ckpt); // which results are we readin
                    try {
                        auto &[scf_r, properties, convergence] = scf_results;
                        scf_r.from_json(j["scf"]);
                        properties.from_json(j["properties"]);
                        convergence.from_json(j["convergence"]);
                    } catch (...) {
                        print("Failed to parse checkpoint file: ", ckpt);
                        scf_results = empty_results;
                    }
                }
                world_.gop.fence();

                if (world_.rank() == 0) {
                    print("Found checkpoint file: ", ckpt);
                    print("results: ", j.dump(4));
                }

                // we could dump params_ to JSON and pass as argv if desired…
                // metadata_(world_);
                set_calc_workdir(pm.dir());

                // Here we validate the results before running

                world_.gop.fence();
                NextAction action;
                action = lib_.valid(world_, scf_results, params_);
                world_.gop.fence();

                if (world_.rank() == 0)
                    print("Next action is ", static_cast<int>(action), " (0=Ok,1=ReloadOnly,2=Restart,3=Redo)");

                if (action == madness::NextAction::Restart || action == madness::NextAction::Redo) {
                    scf_results = lib_.run(world_, params_, action == madness::NextAction::Restart);
                }

                // // Need work (Restart or Redo) — both call run()
                // scf_results = lib_.run(world_, params_);

                results_["scf"] = std::get<0>(scf_results).to_json();
                results_["properties"] = std::get<1>(scf_results).to_json();
                results_["convergence"] = std::get<2>(scf_results).to_json();
                results_["molecule"] = params_.get<Molecule>().to_json();

                // write the checkpoint file
                if (world_.rank() == 0) {
                    std::ofstream ofs(ckpt);
                    ofs << results_.dump(4);
                    ofs.close();
                    print("Written checkpoint file: ", ckpt);
                }
            }
        }

        // std::shared_ptr<SCFApplicationT> scf_app =
        // std::dynamic_pointer_cast<SCFApplicationT>(reference_.shared_from_this());

        nlohmann::json results() const override { return results_; }

    private:
        World &world_;
        Library lib_; // owns shared_ptr<Engine>
        SCFResultsTuple scf_results;
    };

    /**
     * @brief Wrapper application to run the molresponse workflow
     *        via the molresponse_lib::run_response function.
     */
    template<typename Library>
    class ResponseApplication : public Application {
    public:
        /**
         * @param world   MADNESS world communicator
         * @param params  Unified Params containing ResponseParameters & Molecule
         * @param ref_dir   Directory of precomputed ground-state (SCF) outputs
         */
        ResponseApplication(World &world, Params params, const std::shared_ptr<SCF> reference)
            : world_(world), Application(std::move(params)), reference_(std::move(reference)) {
        }

        // print parameters
        void print_parameters(World &world) const override {
            if (world.rank() == 0) {
                std::cout << "Response Parameters:" << std::endl;
                params_.get<ResponseParameters>().print("response");
            }
        }

        /**
         * @brief Execute response + property workflow, writing into workdir/response
         */
        void run(const std::filesystem::path &workdir) override {
            // create a namespaced subdirectory for response outputs
            PathManager pm(workdir, Library::label());
            pm.create(); {
                ScopedCWD scwd(pm.dir());

                auto res = Library::run_response(world_, params_, reference_, pm.dir());

                metadata_ = std::move(res.metadata);
                properties_["response_properties"] = std::move(res.properties);
                properties_["vibrational_analysis"] = std::move(res.vibrational_analysis);
            }
        }

        /**
         * @brief Return a JSON fragment summarizing results
         */
        [[nodiscard]] nlohmann::json results() const override {
            return {{"type", "response"}, {"metadata", metadata_}, {"properties", properties_}};
        }

    private:
        World &world_;
        nlohmann::json metadata_;
        nlohmann::json properties_;
        std::optional<nlohmann::json> vibrational_analysis_;
        const std::shared_ptr<SCF> reference_;
    };

    class CC2Application : public Application, public CC2 {
    public:
        explicit CC2Application(World &w, const Params &p, const std::shared_ptr<Nemo> reference)
            : Application(p), world_(w), reference_(reference),
              CC2(w, p.get<CCParameters>(), p.get<TDHFParameters>(), reference) {
        }

        // print_parameters
        void print_parameters(World &world) const override {
            if (world.rank() == 0) {
                std::cout << "CC2 Parameters:" << std::endl;
            }
        }

        void run(const std::filesystem::path &workdir) override {
            // 1) set up a namedspaced directory for this run
            std::string label = "cc2";
            PathManager pm(workdir, label);
            pm.create();
            world_.gop.fence(); {
                ScopedCWD scwd(pm.dir());
                if (world_.rank() == 0) {
                    std::cout << "Running CC2 in " << pm.dir() << std::endl;
                }

                // 2) define the "checkpoint" file
                auto ckpt = label + "_results.json";
                print("cc checkpoint file", ckpt);
                if (std::filesystem::exists(ckpt)) {
                    if (world_.rank() == 0) {
                        std::cout << "Found checkpoint file: " << ckpt << std::endl;
                    }
                    // read the checkpoint file
                    std::ifstream ifs(ckpt);
                    ifs >> results_;
                    ifs.close();

                    bool ok = true;
                    bool needEnergy = true;
                    if (needEnergy && !results_.contains("energy"))
                        ok = false;
                }

                auto rel = std::filesystem::relative(reference_->work_dir, pm.dir());
                if (world_.rank() == 0) {
                    std::cout << "Running cc2 calculation in: " << pm.dir() << std::endl;
                    std::cout << "Ground state archive: " << reference_->work_dir << std::endl;
                    std::cout << "Relative path: " << rel << std::endl;
                }

                results_ = this->solve();
            }
        }

        nlohmann::json results() const override { return results_; }

    private:
        World &world_;
        const std::shared_ptr<Nemo> reference_;
    };

    class TDHFApplication : public Application, public TDHF {
    public:
        explicit TDHFApplication(World &w, const Params &p, const std::shared_ptr<Nemo> &reference)
            : Application(p), world_(w), reference_(reference), TDHF(w, p.get<TDHFParameters>(), reference) {
        }

        // print_parameters
        void print_parameters(World &world) const override {
            if (world.rank() == 0) {
                std::cout << "TDHF Parameters:" << std::endl;
            }
        }

        void run(const std::filesystem::path &workdir) override {
            // 1) set up a namedspaced directory for this run
            PathManager pm(workdir, "tdhf");
            pm.create();
            world_.gop.fence(); {
                ScopedCWD scwd(pm.dir());
                if (world_.rank() == 0) {
                    std::cout << "Running CIS in " << pm.dir() << std::endl;
                }

                // we could dump params_ to JSON and pass as argv if desired…
                try {
                    const double time_scf_start = wall_time();
                    this->prepare_calculation();
                    const double time_scf_end = wall_time();
                    if (world_.rank() == 0)
                        printf(" at time %.1f\n", wall_time());

                    const double time_cis_start = wall_time();
                    std::vector<CC_vecfunction> roots = this->solve_cis();
                    const double time_cis_end = wall_time();
                    if (world_.rank() == 0)
                        printf(" at time %.1f\n", wall_time());

                    if (world_.rank() == 0) {
                        std::cout << std::setfill(' ');
                        std::cout << "\n\n\n";
                        std::cout << "--------------------------------------------------\n";
                        std::cout << "MRA-CIS ended \n";
                        std::cout << "--------------------------------------------------\n";
                        std::cout << std::setw(25) << "time scf" << " = " << time_scf_end - time_scf_start << "\n";
                        std::cout << std::setw(25) << "time cis" << " = " << time_cis_end - time_cis_start << "\n";
                        std::cout << "--------------------------------------------------\n";
                    }
                    auto j = this->analyze(roots);
                    // funnel through CISResults to make sure we have the right format
                    CISResults results(j);
                    results_ = results.to_json();
                } catch (std::exception &e) {
                    print("Caught exception: ", e.what());
                }
            }
        }

        nlohmann::json results() const override { return results_; }

    private:
        World &world_;
        std::shared_ptr<Nemo> reference_;
        std::filesystem::path ref_dir_;
    };

    class OEPApplication : public Application, public OEP {
    public:
        explicit OEPApplication(World &w, const Params &p, const std::shared_ptr<Nemo> &reference)
            : Application(p), world_(w), reference_(reference), OEP(w, p.get<OEP_Parameters>(), reference) {
        }

        // print_parameters
        void print_parameters(World &world) const override {
            if (world.rank() == 0) {
                std::cout << "OEP Parameters:" << std::endl;
            }
        }

        void run(const std::filesystem::path &workdir) override {
            // 1) set up a namedspaced directory for this run
            PathManager pm(workdir, "oep");
            pm.create();
            world_.gop.fence(); {
                ScopedCWD scwd(pm.dir());
                if (world_.rank() == 0) {
                    std::cout << "Running OEP in " << pm.dir() << std::endl;
                }

                // 2) define the "checkpoint" file
                std::string label = "oep";
                auto ckpt = label + "_results.json";
                print("cc checkpoint file", ckpt);
                if (std::filesystem::exists(ckpt)) {
                    if (world_.rank() == 0) {
                        std::cout << "Found checkpoint file: " << ckpt << std::endl;
                    }
                    // read the checkpoint file
                    std::ifstream ifs(ckpt);
                    nlohmann::json j;
                    ifs >> j;
                    ifs.close();
                }

                // we could dump params_ to JSON and pass as argv if desired…
                try {
                    const double time_scf_start = wall_time();
                    this->value();
                    const double time_scf_end = wall_time();
                    if (world_.rank() == 0)
                        printf(" at time %.1f\n", wall_time());

                    if (world_.rank() == 0) {
                        std::cout << std::setfill(' ');
                        std::cout << "\n\n\n";
                        std::cout << "--------------------------------------------------\n";
                        std::cout << "MRA-OEP ended \n";
                        std::cout << "--------------------------------------------------\n";
                        std::cout << std::setw(25) << "time scf" << " = " << time_scf_end - time_scf_start << "\n";
                        std::cout << "--------------------------------------------------\n";
                    }
                } catch (std::exception &e) {
                    print("Caught exception: ", e.what());
                }
                // nlohmann::json results;
                results_ = this->analyze();
            }
        }

        nlohmann::json results() const override { return results_; }

    private:
        World &world_;
        std::shared_ptr<Nemo> reference_;

        double energy_;
        std::optional<Tensor<double> > dipole_;
        std::optional<Tensor<double> > gradient_;
        std::optional<real_function_3d> density_;
    };

    inline NextAction decide_next_action(bool at_protocol, bool archive_needed, bool archive_exists,
                                         bool all_properties_computed, bool restart_exists) {
        // We must recompute if any of these are true:
        const bool must_redo = !at_protocol // not at final protocol
                               || (archive_needed && !archive_exists) // user wants archive but it's missing
                               || !all_properties_computed; // a requested prop is missing

        if (!must_redo) {
            // We’re at final protocol, have all requested props, and either
            // the archive is present or not required.
            // If the archive isn't there but also not required, it's just a reload.
            if (!archive_exists && !archive_needed)
                return NextAction::ReloadOnly;
            return NextAction::Ok;
        }

        // We need work; decide between Restart vs Redo:
        return restart_exists ? NextAction::Restart : NextAction::Redo;
    }

    template<typename SCFParams>
    NextAction valid(World &world, const SCFResultsTuple &results, const SCFParams &params) {
        // Take a copy of the parameters
        auto [sr, pr, cr] = results;

        // Required convergence for "final" protocol
        const auto vthresh = params.protocol().back(); // final protocol
        const auto vdconv = params.dconv();
        const bool archive_needed = params.save();

        // Requested outputs
        const bool need_energy = true;
        const bool need_dipole = params.dipole();
        const bool need_gradient = params.derivatives();

        // Files/paths
        const std::string archivename = params.prefix();
        const auto restart_path = std::filesystem::path(archivename + ".restartdata.00000");
        const bool archive_exists = std::filesystem::exists(restart_path);
        if (world.rank() == 0) {
            print("Restart file: ", restart_path.string());
        }

        // State in resultout the threshold refinement.
        //
        const bool at_protocol = (cr.converged_for_thresh == vthresh && cr.converged_for_dconv == vdconv);

        const auto pjson = pr.to_json();
        const bool energy_ok = pjson.contains("energy");
        const bool dipole_ok = pjson.contains("dipole");
        const bool gradient_ok = pjson.contains("gradient");

        // Only require props the user asked for
        const bool all_properties_computed =
                (need_energy ? energy_ok : true) && (need_dipole ? dipole_ok : true) && (need_gradient
                        ? gradient_ok
                        : true);

        // Decide action
        const bool must_redo = !at_protocol || (archive_needed && !archive_exists) || !all_properties_computed;

        // if we don't need to redo, we can either reload or return ok
        if (!must_redo)
            return (!archive_exists && !archive_needed) ? NextAction::ReloadOnly : NextAction::Ok;

        if (world.rank() == 0) {
            print("at_protocol: ", at_protocol);
            print("archive_needed: ", archive_needed);
            print("archive_exists: ", archive_exists);
            print("all_properties_computed: ", all_properties_computed);
        }
        // with we need to redo we can restart from the exisiting archive
        return archive_exists ? NextAction::Restart : NextAction::Redo;
    }

    struct moldft_lib {
        static constexpr char const *label() { return "moldft"; }
        NextAction next_action_ = NextAction::Ok;

        vector<double> protocol;
        SCFResultsTuple last_results_;

        using Calc = SCF;

        NextAction valid(World &world, const SCFResultsTuple &results, const Params &params) {
            last_results_ = results;
            next_action_ = ::valid(world, results, params.get<CalculationParameters>());
            return next_action_;
        }

        // expose the live engine
        std::shared_ptr<Calc> calc(World &world, const Params &params) {
            if (!calc_)
                initialize_(world, params); // create once
            return calc_;
        }

        void print_parameters() const { calc_->print_parameters(); }
        // params get's changed by SCF constructor
        SCFResultsTuple run(World &world, const Params &params, const bool restart = false) {
            auto moldft_params = params.get<CalculationParameters>();
            const auto &molecule = params.get<Molecule>();
            auto params_copy = params;

            SCFResults sr;

            if (next_action_ == NextAction::Ok || next_action_ == NextAction::ReloadOnly) {
                // nothing to do
                return last_results_;
            }

            if (restart) {
                // Handle restart logic
                auto cr = std::get<2>(last_results_);
                moldft_params.set_user_defined_value("restart", true);
            }
            auto scf = calc(world, params_copy);
            // redirect any log files into outdir if needed…
            // Warm and fuzzy for the user
            if (world.rank() == 0) {
                print("\n\n");
                print(" MADNESS Hartree-Fock and Density Functional Theory Program");
                print(" ----------------------------------------------------------\n");
                print("\n");
                //   scf->molecule.print();
                print("\n");
                //    scf->param.print("dft");
            }
            // Come up with an initial OK data map
            if (world.size() > 1) {
                scf->set_protocol<3>(world, 1e-4);
                scf->make_nuclear_potential(world);
                scf->initial_load_bal(world);
            }
            // vama
            scf->set_protocol<3>(world, scf->param.protocol()[0]);
            double energy = 0.0;

            if (scf->param.gopt()) {
                MolOpt opt(scf->param.gmaxiter(), 0.1, scf->param.gval(), scf->param.gtol(),
                           1e-3, // XTOL
                           1e-5, scf->param.gprec(),
                           (world.rank() == 0) ? 1 : 0, // print_level
                           scf->param.algopt());

                MolecularEnergy target(world, *scf);
                auto new_mol = opt.optimize(scf->molecule, target);
                energy = scf->current_energy;
                sr.scf_molecule = new_mol;
                sr.is_opt = true;

                // write out the optimized geometry
                if (world.rank() == 0) {
                    std::string geomfile = scf->param.prefix() + "_opt.xyz";
                    std::ofstream ofs(geomfile);
                    new_mol.print(ofs);
                    ofs.close();
                    print("optimized geometry written to ", geomfile);
                    // write out mad.in with optimized geometry
                }

                // MolecularEnergy E(world, *scf);
                // energy = E.value(new_mol.get_all_coords().flat());
            } else {
                MolecularEnergy E(world, *scf);
                sr.scf_molecule = molecule;
                sr.is_opt = false;

                energy = E.value(scf->molecule.get_all_coords().flat());
                if (world.rank() == 0 && scf->param.print_level() > 0)
                    E.output_calc_info_schema();
            }

            functionT rho = scf->make_density(world, scf->aocc, scf->amo);
            functionT brho = rho;
            if (scf->param.nbeta() != 0 && !scf->param.spin_restricted())
                brho = scf->make_density(world, scf->bocc, scf->bmo);
            rho.gaxpy(1.0, brho, 1.0);

            // optionally compute gradient, dipole, etc.
            Tensor<double> grad;
            if (scf->param.derivatives()) {
                grad = scf->derivatives(world, rho);
                scf->e_data.add_gradient(grad);
            }

            tensorT dip;
            if (scf->param.dipole())
                dip = scf->dipole(world, scf->make_density(world, scf->aocc, scf->amo));

            scf->do_plots(world);

            ConvergenceResults cr;
            cr.set_converged_thresh(FunctionDefaults<3>::get_thresh());
            cr.set_converged_dconv(scf->converged_for_dconv);
            PropertyResults pr;
            pr.energy = energy;
            pr.dipole = dip;
            pr.gradient = grad;
            sr.aeps = scf->aeps;
            sr.beps = scf->beps;
            sr.properties = pr;

            return {sr, pr, cr};
        }

    private:
        void initialize_(World &world, const Params &params) {
            // write mad.in if missing
            const auto &cp = params.get<CalculationParameters>();
            const auto &mol = params.get<Molecule>();

            world.gop.fence();
            if (world.rank() == 0) {
                if (true) {
                    // should always overwrite for now

                    json in;
                    in["dft"] = cp.to_json_if_precedence("defined");
                    in["molecule"] = mol.to_json_if_precedence("defined");
                    std::ofstream ofs("mad.in");
                    write_json_to_input_file(in, {"dft"}, ofs);
                    mol.print_defined_only(ofs);
                }
            }
            world.gop.fence();

            commandlineparser parser;
            parser.set_keyval("input", "mad.in");
            if (world.rank() == 0)
                ::print("input filename: ", parser.value("input"));

            FunctionDefaults<3>::set_pmap(pmapT(new LevelPmap<Key<3> >(world)));
            std::cout.precision(6);
            calc_ = std::make_shared<SCF>(world, parser);
        }

        std::shared_ptr<Calc> calc_;
    }; // namespace moldft_lib

    struct nemo_lib {
        using Calc = Nemo;
        static constexpr char const *label() { return "nemo"; }

        std::shared_ptr<Calc> calc(World &world, const Params &params) {
            if (!nemo_)
                initialize_(world, params);
            return nemo_;
        }

        NextAction valid(World &world, const SCFResultsTuple &results, const Params &params) {
            // Take a copy of the parameters
            return ::valid(world, results, params.get<CalculationParameters>());
        }

        void print_parameters() const { nemo_->print_parameters(); }

        SCFResultsTuple run(World &world, const Params &params, const bool restart = false) {
            auto nm = calc(world, params);
            nm->get_calc()->work_dir = std::filesystem::current_path();

            nm->value();
            PropertyResults pr = nm->analyze();
            // compute the hessian

            ConvergenceResults cr;
            cr.set_converged_thresh(nm->get_calc()->converged_for_thresh);
            cr.set_converged_dconv(nm->get_calc()->converged_for_dconv);

            SCFResults sr;
            sr.aeps = nm->get_calc()->aeps;
            sr.beps = nm->get_calc()->beps;
            sr.properties = pr;
            sr.scf_total_energy = nm->get_calc()->current_energy;

            if (nm->get_nemo_param().hessian())
                sr.properties->vibrations = nm->hessian(nm->get_calc()->molecule.get_all_coords());

            return {sr, pr, cr};
        }

    private:
        void initialize_(World &world, const Params &params) {
            nemo_ = std::make_shared<Nemo>(world, params.get<CalculationParameters>(),
                                           params.get<Nemo::NemoCalculationParameters>(), params.get<Molecule>());
        }

        std::shared_ptr<Calc> nemo_;
    };
} // namespace madness
