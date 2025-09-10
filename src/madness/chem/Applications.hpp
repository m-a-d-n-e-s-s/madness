#pragma once

#include <madness/chem/Results.h>

#include <madness/chem/InputWriter.hpp>
#include <madness/chem/ParameterManager.hpp>
#include <madness/chem/PathManager.hpp>
#include <span>

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
  explicit Application(const Params &p) : params_(p) {}
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
      std::cout << "Found checkpoint file: " << filename << std::endl;
      // std::ifstream ifs(workdir_ / filename);
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
      std::string msg = "Results file " + filename + " does not exist in " + std::filesystem::current_path().string();
      MADNESS_EXCEPTION(msg.c_str(), 1);
    }
    return nlohmann::json();
  }

protected:
  const Params params_;
  nlohmann::json results_;
};

template <typename Library> class SCFApplication : public Application {
private:
public:
  using Calc = typename Library::Calc;

  explicit SCFApplication(World &w, const Params &p) : Application(p), world_(w) {}

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
    pm.create();
    {
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

      if (world_.rank() == 0) {
        print("Found checkpoint file: ", ckpt);
        print("results: ", j.dump(4));
      }

      // we could dump params_ to JSON and pass as argv if desired…
      // metadata_(world_);
      set_calc_workdir(pm.dir());

      // Here we validate the results before running

      auto action = lib_.valid(scf_results, params_);

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
template <typename Library> class ResponseApplication : public Application {
public:
  /**
   * @param world   MADNESS world communicator
   * @param params  Unified Params containing ResponseParameters & Molecule
   * @param ref_dir   Directory of precomputed ground-state (SCF) outputs
   */
  ResponseApplication(World &world, Params params, const std::shared_ptr<SCF> reference)
      : world_(world), Application(std::move(params)), reference_(std::move(reference)) {}
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
    pm.create();
    {
      ScopedCWD scwd(pm.dir());

      auto res = Library::run_response(world_, params_, reference_, pm.dir());

      metadata_ = std::move(res.metadata);
      properties_ = std::move(res.properties);
      vibrational_analysis_ = std::move(res.vibrational_analysis);
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
        CC2(w, p.get<CCParameters>(), p.get<TDHFParameters>(), reference) {}
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
    world_.gop.fence();
    {
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
      : Application(p), world_(w), reference_(reference), TDHF(w, p.get<TDHFParameters>(), reference) {}

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
    world_.gop.fence();
    {
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
      : Application(p), world_(w), reference_(reference), OEP(w, p.get<OEP_Parameters>(), reference) {}

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
    world_.gop.fence();
    {
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
  std::optional<Tensor<double>> dipole_;
  std::optional<Tensor<double>> gradient_;
  std::optional<real_function_3d> density_;
};

} // namespace madness
