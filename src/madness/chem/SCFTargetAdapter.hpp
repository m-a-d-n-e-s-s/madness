// SCFTarget.hpp
#pragma once
#include <filesystem>

#include <madness/chem/Applications.hpp>
#include <madness/chem/molecule.h>

struct SCFTarget {
  World& world_;
  std::function<std::unique_ptr<Application>(Params)> factory_;
  Params params_;
  double last_energy = 0.0;
  madness::Tensor<double> last_gradient;

  SCFTarget(World& w,
            std::function<std::unique_ptr<Application>(Params)> factory,
            Params p)
      : world_(w), factory_(std::move(factory)), params_(std::move(p)) {}

  // Called by MolOpt at each geometry:
  void energy_and_gradient(Molecule& mol, double& energy,
                           madness::Tensor<double>& grad) {
    // 1) inject new coords into params
    params_.set(mol);

    // 2) build & run the SCF+grad application
    auto app = factory_(params_);
    app->run(std::filesystem::current_path());
    auto res = app->results();

    // 3) extract energy + gradient
    energy = res.at("energy").get<double>();
    grad = tensor_from_json<double>(res.at("gradient"));

    last_energy = energy;
    last_gradient = grad;
  }

  // Called by MolOpt::line_search
  double value(const madness::Tensor<double>& x) {
    Molecule mol = params_.get<Molecule>();
    mol.set_all_coords(x.reshape(mol.natom(), 3));
    double e;
    madness::Tensor<double> g;
    energy_and_gradient(mol, e, g);
    return e;
  }
};
