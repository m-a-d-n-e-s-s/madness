#pragma once

#include "InputWriter.hpp"
#include "ParameterManager.hpp"
#include "PathManager.hpp"

// Scoped CWD: changes the current directory to the given one, and restores when
// the object goes out of scope
struct ScopedCWD {
  std::filesystem::path old_cwd;
  explicit ScopedCWD(std::filesystem::path const& new_dir) {
    old_cwd = std::filesystem::current_path();
    std::filesystem::current_path(new_dir);
  }
  ~ScopedCWD() { std::filesystem::current_path(old_cwd); }
};

// common interface for any underlying MADNESS app
class Application {
 public:
  explicit Application(Params p) : params_(std::move(p)) {}
  virtual ~Application() = default;

    // run: write all outputs under the given directory
    virtual void run(const std::filesystem::path& workdir) = 0;

    // optional hook to return a JSON fragment of this app's main results
    [[nodiscard]] virtual nlohmann::json results() const = 0;

 protected:
  Params params_;
};

  template <typename Library>
  class SCFApplication : public Application {
  public:
    explicit SCFApplication(World& w, Params p)
        : Application(std::move(p)), world_(w) {}

    void run(const std::filesystem::path& workdir) override {
      // 1) set up a namedspaced directory for this run
      PathManager pm(workdir, Library::label());
      pm.create();
      world_.gop.fence();
      {
        ScopedCWD scwd(pm.dir());
        if (world_.rank() == 0) {
          std::cout << "Running SCF in " << pm.dir() << std::endl;
        }

        auto scfParams = params_.get<CalculationParameters>();
        bool needEnergy = true;
        bool needDipole = scfParams.dipole();
        bool needGradient = scfParams.derivatives();

        // 2) define the "checkpoint" file
        auto ckpt = "scf_results.json";
        if (std::filesystem::exists(ckpt)) {
          if (world_.rank() == 0) {
            std::cout << "Found checkpoint file: " << ckpt << std::endl;
          }
          // read the checkpoint file
          std::ifstream ifs(ckpt);
          nlohmann::json j;
          ifs >> j;
          ifs.close();

          bool ok = true;
          if (needEnergy && !j.contains("energy")) ok = false;
          if (needDipole && !j.contains("dipole")) ok = false;
          if (needGradient && !j.contains("gradient")) ok = false;

          if (ok) {
            energy_ = j["energy"];
            if (needDipole) dipole_ = tensor_from_json<double>(j["dipole"]);
            if (needGradient) gradient_ = tensor_from_json<double>(j["gradient"]);
            return;
          }
        }

        // we could dump params_ to JSON and pass as argv if desired…
        auto results = Library::run_scf(world_, params_, pm.dir());

        energy_ = results.energy;
        dipole_ = results.dipole;
        gradient_ = results.gradient;

        // 5) write out JSON for future restarts
        nlohmann::json outj = {{"energy", energy_}};
        if (dipole_->has_data()) outj["dipole"] = tensor_to_json(*dipole_);
        if (gradient_->has_data()) outj["gradient"] = tensor_to_json(*gradient_);

        if (world_.rank() == 0) {
          std::cout << "Writing checkpoint file: " << ckpt << std::endl;
          std::ofstream o(ckpt);
          o << std::setw(2) << outj << std::endl;
        }
      }
    }

    nlohmann::json results() const override {
      auto scfParams = params_.get<CalculationParameters>();
      nlohmann::json j = {
        {"type", "scf"},
        {"energy", energy_},
    };
      if (dipole_ && scfParams.dipole()) j["dipole"] = tensor_to_json(*dipole_);
      if (gradient_ && scfParams.derivatives())
        j["gradient"] = tensor_to_json(*gradient_);

      return j;
    }

  private:
    World& world_;

    double energy_;

    std::optional<Tensor<double>> dipole_;
    std::optional<Tensor<double>> gradient_;
    std::optional<real_function_3d> density_;
  };

  /**
   * @brief Wrapper application to run the molresponse workflow
   *        via the molresponse_lib::run_response function.
   */
  template <typename Library>
  class ResponseApplication : public Application {
  public:
    /**
     * @param world   MADNESS world communicator
     * @param params  Unified Params containing ResponseParameters & Molecule
     * @param indir   Directory of precomputed ground-state (SCF) outputs
     */
    ResponseApplication(World& world, Params params, std::filesystem::path indir)
        : world_(world),
          Application(std::move(params)),
          indir_(std::move(indir)) {}

    /**
     * @brief Execute response + property workflow, writing into workdir/response
     */
    void run(const std::filesystem::path& workdir) override {
      // create a namespaced subdirectory for response outputs
      PathManager pm(workdir, Library::label());
      pm.create();
      {
        ScopedCWD scwd(pm.dir());

        auto res = Library::run_response(world_, params_, indir_, pm.dir());

        metadata_ = std::move(res.metadata);
        properties_ = std::move(res.properties);
      }
    }

    /**
     * @brief Return a JSON fragment summarizing results
     */
    [[nodiscard]] nlohmann::json results() const override {
      return {{"type", "response"},
              {"metadata", metadata_},
              {"properties", properties_}};
    }

  private:
    World& world_;
    std::filesystem::path indir_;
    nlohmann::json metadata_;
    nlohmann::json properties_;
  };

  template <typename Library>
  class CC2Application : public Application {
  public:
    explicit CC2Application(World& w, Params p, const commandlineparser& parser, path& gsdir)
        : Application(std::move(p)), world_(w), parser(parser), gsdir(gsdir) {}

    void run(const std::filesystem::path& workdir) override {
      // 1) set up a namedspaced directory for this run
      PathManager pm_reference(workdir, Library::label());
      PathManager pm(workdir, Library::label());
      auto scf_parameters=params_.get<CalculationParameters>();
      pm.create();
      world_.gop.fence();
      {
        ScopedCWD scwd(pm.dir());
        if (world_.rank() == 0) {
          std::cout << "Running CC2 in " << pm.dir() << std::endl;
        }

        // read reference wave function from ground-state directory
        auto nemo=std::shared_ptr<Nemo>();
        auto relative_gsdir= std::filesystem::relative(gsdir, pm.dir());
        {
          ScopedCWD gs(relative_gsdir);
          nemo.reset(new Nemo(world_,parser));
        }

        // 2) define the "checkpoint" file
        std::string label=Library::label();
        auto ckpt = label + "_results.json";
        print("cc checkpoint file",ckpt);
        if (std::filesystem::exists(ckpt)) {
          if (world_.rank() == 0) {
            std::cout << "Found checkpoint file: " << ckpt << std::endl;
          }
          // read the checkpoint file
          std::ifstream ifs(ckpt);
          nlohmann::json j;
          ifs >> j;
          ifs.close();

          bool ok = true;
          bool needEnergy = true;
          if (needEnergy && !j.contains("energy")) ok = false;

          if (ok) {
            energy_ = j["energy"];
            return;
          }
        }

        // we could dump params_ to JSON and pass as argv if desired…
        auto results = Library::run_cc2(world_, params_, gsdir, pm.dir());

        energy_ = results.energy;
        dipole_ = results.dipole;
        gradient_ = results.gradient;

        // 5) write out JSON for future restarts
        nlohmann::json outj = {{"energy", energy_}};
        if (dipole_->has_data()) outj["dipole"] = tensor_to_json(*dipole_);
        if (gradient_->has_data()) outj["gradient"] = tensor_to_json(*gradient_);

        if (world_.rank() == 0) {
          std::cout << "Writing checkpoint file: " << ckpt << std::endl;
          std::ofstream o(ckpt);
          o << std::setw(2) << outj << std::endl;
        }
      }
    }

    nlohmann::json results() const override {
      auto scfParams = params_.get<CalculationParameters>();
      nlohmann::json j = {
        {"type", "scf"},
        {"energy", energy_},
    };
      if (dipole_ && scfParams.dipole()) j["dipole"] = tensor_to_json(*dipole_);
      if (gradient_ && scfParams.derivatives())
        j["gradient"] = tensor_to_json(*gradient_);

      return j;
    }

  private:
    World& world_;
    double energy_;
    commandlineparser parser;
    path gsdir;

    std::optional<Tensor<double>> dipole_;
    std::optional<Tensor<double>> gradient_;
    std::optional<real_function_3d> density_;
  };




