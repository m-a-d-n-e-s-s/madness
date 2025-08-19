#pragma once

#include <madness/chem/Results.h>

#include <madness/chem/InputWriter.hpp>
#include <madness/chem/ParameterManager.hpp>
#include <madness/chem/PathManager.hpp>

namespace madness {
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

/// common interface for any underlying MADNESS app
class Application {
public:
  explicit Application(const Params &p) : params_(p) {}
  virtual ~Application() = default;

  // run: write all outputs under the given directory
  virtual void run(const std::filesystem::path &workdir) = 0;

  // optional hook to return a JSON fragment of this app's main results
  [[nodiscard]] virtual nlohmann::json results() const = 0;

  // // get the parameters used for this application
  // [[nodiscard]] virtual const QCCalculationParametersBase &get_parameters() const = 0;

  virtual void print_parameters(World &world) const = 0;

  // get the working directory for this application
  [[nodiscard]] path get_workdir() const { return workdir_; }

  /// check if this calculation has a json with results
  [[nodiscard]] virtual bool has_results(std::string filename) const {
    // check if the results file exists
    // return std::filesystem::exists(workdir_ / filename);
    return std::filesystem::exists(filename);
  }

  // [[nodiscard]] virtual bool verify_results(const nlohmann::json &j) const {
  //   // check if some key parameters of the calculation match:
  //   // molecule, box size, nmo_alpha, nmo_beta
  //   Molecule mol1 = params_.get<Molecule>();
  //   Molecule mol2;
  //   mol2.from_json(j["molecule"]);
  //   if (not(mol1 == mol2)) {
  //     print("molecule mismatch");
  //     mol1.print();
  //     mol2.print();
  //     return false;
  //   }
  //   return true;
  // }

  // /// read the results from a json file
  // [[nodiscard]] virtual nlohmann::json read_results(std::string filename) const {
  //   if (has_results(filename)) {
  //     std::cout << "Found checkpoint file: " << filename << std::endl;
  //     // std::ifstream ifs(workdir_ / filename);
  //     std::ifstream ifs(filename);
  //     nlohmann::json j;
  //     ifs >> j;
  //     ifs.close();
  //     if (not verify_results(j)) {
  //       std::string msg = "Results file " + filename + " does not match the parameters of the calculation";
  //       print(msg);
  //       return nlohmann::json(); // return empty json
  //     }
  //     return j;
  //   } else {
  //     std::string msg = "Results file " + filename + " does not exist in " + workdir_.string();
  //     MADNESS_EXCEPTION(msg.c_str(), 1);
  //   }
  //   return nlohmann::json();
  // }

  // /// check if the wavefunctions are already computed
  // [[nodiscard]] virtual bool has_wavefunctions(std::string filename) const {
  //   return std::filesystem::exists(workdir_ / filename);
  // }

  // /// read the wavefunctions from a file
  // [[nodiscard]] virtual std::vector<double> read_wavefunctions(std::string filename) const {
  //   if (has_wavefunctions(filename)) {
  //     std::ifstream ifs(workdir_ / filename);
  //     std::vector<double> wfs;
  //     double value;
  //     while (ifs >> value) {
  //       wfs.push_back(value);
  //     }
  //     ifs.close();
  //     return wfs;
  //   } else {
  //     std::string msg = "Wavefunction file " + filename + " does not exist in " + workdir_.string();
  //     MADNESS_EXCEPTION(msg.c_str(), 1);
  //   }
  //   return {};
  // }

  // /// check if this calculation needs to be redone
  // [[nodiscard]] virtual bool needs_redo() const {
  //   // read json and check if the results are already there
  //   if (std::filesystem::exists(workdir_ / "results.json")) {
  //     std::ifstream ifs(workdir_ / "results.json");
  //     nlohmann::json j;
  //     ifs >> j;
  //     ifs.close();
  //     // check if the results are already there
  //     return !j.contains("energy") || !j["energy"].is_number();
  //   }
  //   // by default, we assume that the calculation needs to be redone
  //   return true;
  // }

protected:
  const Params params_;
  path workdir_;
  nlohmann::json results_;
};

template <typename Library> class SCFApplication : public Application {
public:
  using Engine = typename Library::Engine;

  explicit SCFApplication(World &w, const Params &p) : Application(p), world_(w) {}

  // Give downstream steps the live calc
  std::shared_ptr<Engine> engine() { return lib_.engine(world_, params_); }

  // print parameters
  void print_parameters(World &world) const override {
    if (world.rank() == 0) {
      std::cout << "SCF Parameters:" << std::endl;
    }
  }

  void run(const std::filesystem::path &workdir) override {
    // 1) set up a namedspaced directory for this run
    std::string label = Library::label();
    PathManager pm(workdir, label.c_str());
    pm.create();
    workdir_ = pm.dir();
    world_.gop.fence();
    {
      ScopedCWD scwd(pm.dir());
      if (world_.rank() == 0)
        (world_, params_);
      {
        std::cout << "Running SCF in " << pm.dir() << std::endl;
      }

      // we could dump params_ to JSON and pass as argv if desired…
      MetaDataResults metadata(world_);
      results_ = lib_.run(world_, params_);

      // } else if constexpr (std::is_same_v<ScfT, Nemo>) {
      //   results_ = Library::run_nemo(this->get_nemo());
      // } else {
      //   MADNESS_CHECK_THROW("unknown SCF type", 1);
      // }
      results_["metadata"] = metadata.to_json();

      // legacy
      PropertyResults properties = results_["properties"];
      energy_ = properties.energy;
      dipole_ = properties.dipole;
      gradient_ = properties.gradient;
    }
  }
  // std::shared_ptr<SCFApplicationT> scf_app =
  // std::dynamic_pointer_cast<SCFApplicationT>(reference_.shared_from_this());

  nlohmann::json results() const override { return results_; }

private:
  World &world_;
  Library lib_; // owns shared_ptr<Engine>

  double energy_;

  std::optional<Tensor<double>> dipole_;
  std::optional<Tensor<double>> gradient_;
  std::optional<real_function_3d> density_;
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
  ResponseApplication(World &world, Params params, std::filesystem::path ref_dir)
      : world_(world), Application(std::move(params)), ref_dir_(std::move(ref_dir)) {}
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

      auto res = Library::run_response(world_, params_, ref_dir_, pm.dir());

      metadata_ = std::move(res.metadata);
      properties_ = std::move(res.properties);
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
  std::filesystem::path ref_dir_; // directory of precomputed SCF outputs
};

class CC2Application : public Application, public CC2 {
public:
  explicit CC2Application(World &w, const Params &p, const std::shared_ptr<Nemo> reference,
                          const std::filesystem::path &ref_dir)
      : Application(p), world_(w), reference_(reference), ref_dir_(ref_dir),
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

      auto rel = std::filesystem::relative(ref_dir_, pm.dir());
      if (world_.rank() == 0) {
        std::cout << "Running cc2 calculation in: " << pm.dir() << std::endl;
        std::cout << "Ground state archive: " << ref_dir_ << std::endl;
        std::cout << "Relative path: " << rel << std::endl;
      }

      results_ = this->solve();
    }
  }

  nlohmann::json results() const override { return results_; }

private:
  World &world_;
  const std::shared_ptr<Nemo> reference_;
  const std::filesystem::path ref_dir_;
};

class TDHFApplication : public Application, public TDHF {
public:
  explicit TDHFApplication(World &w, const Params &p, const std::shared_ptr<Nemo> &reference,
                           const std::filesystem::path &ref_dir)
      : Application(p), world_(w), reference_(reference), ref_dir_(ref_dir),
        TDHF(w, p.get<TDHFParameters>(), reference) {}

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
  explicit OEPApplication(World &w, const Params &p, const std::shared_ptr<Nemo> &reference,
                          const std::filesystem::path &ref_dir)
      : Application(p), world_(w), reference_(reference), OEP(w, p.get<OEP_Parameters>(), reference),
        ref_dir_(ref_dir) {}

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
  std::filesystem::path ref_dir_;

  double energy_;
  std::optional<Tensor<double>> dipole_;
  std::optional<Tensor<double>> gradient_;
  std::optional<real_function_3d> density_;
};

} // namespace madness
