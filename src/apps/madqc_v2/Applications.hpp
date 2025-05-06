#pragma once

#include "MoldftLib.hpp"
#include "MolresponseLib.hpp"
#include "PathManager.hpp"
#include "Utils.hpp"

struct ScopedCWD {
  std::filesystem::path old_cwd;
  ScopedCWD(std::filesystem::path const& new_dir) {
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

class SCFApplication : public Application {
 public:
  SCFApplication(World& w, Params p) : Application(std::move(p)), world_(w) {}

  void run(const std::filesystem::path& workdir) override {
    // 1) set up a namedspaced directory for this run
    PathManager pm(workdir, "dft");
    pm.create();
    {
      ScopedCWD scwd(pm.dir());

      // 2) define the "checkpoint" file
      auto ckpt = pm.dir() / "scf_results.json";
      if (std::filesystem::exists(ckpt)) {
        // read the checkpoint file
        std::ifstream ifs(ckpt);
        nlohmann::json j;
        ifs >> j;
        ifs.close();
        energy_ = j["energy"];
        dipole_ = j["dipole"];
        gradient_ = tensor_from_json<double>(j["gradient"]);
        return;
      }

      // we could dump params_ to JSON and pass as argv if desired…
      auto results = moldft_lib::run_scf(world_, params_, workdir / "dft");
      energy_ = results.energy;
      dipole_ = results.dipole;
      gradient_ = results.gradient;

      // 5) write out JSON for future restarts
      nlohmann::json outj = {{"energy", energy_}, {"dipole", tensor_to_json(dipole_)}};
      if (gradient_) outj["gradient"] = tensor_to_json(*gradient_);

      std::ofstream o(ckpt);
      o << std::setw(2) << outj << std::endl;
    }
  }

  nlohmann::json results() const override {
    nlohmann::json j = {
        {"type", "scf"},
        {"energy", energy_},
        {"dipole", {dipole_[0], dipole_[1], dipole_[2]}},
    };
    if (gradient_) j["gradient"] = tensor_to_json(*gradient_);
    return j;
  }

 private:
  World& world_;

  double energy_;
  tensorT dipole_;
  std::optional<Tensor<double>> gradient_;
};

/**
 * @brief Wrapper application to run the molresponse workflow
 *        via the molresponse_lib::run_response function.
 */
class ResponseApplication : public Application {
 public:
  /**
   * @param world   MADNESS world communicator
   * @param params  Unified Params containing ResponseParameters & Molecule
   * @param indir   Directory of precomputed ground-state (SCF) outputs
   */
  ResponseApplication(World& world, Params params, std::filesystem::path indir) : world_(world), Application(std::move(params)), indir_(std::move(indir)) {}

  /**
   * @brief Execute response + property workflow, writing into workdir/response
   */
  void run(const std::filesystem::path& workdir) override {
    // create a namespaced subdirectory for response outputs
    PathManager pm(workdir, "response");
    pm.create();
    {
      ScopedCWD scwd(pm.dir());

      molresponse_lib::Results res = molresponse_lib::run_response(world_, params_, indir_, pm.dir());

      metadata_ = std::move(res.metadata);
      properties_ = std::move(res.properties);
    }
  }

  /**
   * @brief Return a JSON fragment summarizing results
   */
  [[nodiscard]] nlohmann::json results() const override { return {{"type", "response"}, {"metadata", metadata_}, {"properties", properties_}}; }

 private:
  World& world_;
  std::filesystem::path indir_;
  nlohmann::json metadata_;
  nlohmann::json properties_;
};
