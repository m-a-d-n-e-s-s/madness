#include <../molresponse_v2/ResponseParameters.hpp>

#include "ParameterManager.hpp"
#include "PathManager.hpp"

// Define a concrete aliased ParameterManager type
using Params = ParameterManager<CalculationParameters, ResponseParameters, OptimizationParameters, Molecule>;

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

// e.g. single‐point DFT via moldft API
class DFTApplication : public Application {
 public:
  using Application::Application;

  void run(const std::filesystem::path& workdir) override {
    auto pm = PathManager(workdir, /*name=*/"dft");
    pm.create();  // mkdir workdir/dft
    // invoke moldft under the hood, writing its orbitals, log, etc into pm.dir()
    moldft::run(params_, pm.dir());
    // maybe read back the energy:
    energy_ = moldft::readEnergy(pm.dir());
  }

  nlohmann::json results() const override { return {{"type", "dft"}, {"energy", energy_}}; }

 private:
  double energy_;
};
