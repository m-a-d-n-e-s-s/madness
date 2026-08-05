// SCFTargetAdapter.hpp
#pragma once

#include <madness/chem/molecule.h>
#include <madness/tensor/tensor.h>

namespace madness {

/// Adapts an engine that exposes `value(x)` / `gradient(x)` to the target
/// protocol MolOpt expects: `energy_and_gradient(const Molecule&, double&,
/// Tensor<double>&)` plus `value(const Tensor<double>&)`.
///
/// SCF already has such an adapter in `MolecularEnergy` (SCF.h), which is what
/// the in-SCF `dft gopt` path uses. `Nemo` and everything derived from it
/// implement `MolecularOptimizationTargetInterface` (`value(x)`, `gradient(x)`,
/// `molecule()`) but not MolOpt's Molecule-based signature, so this fills that gap
/// without touching the engines -- which is what lets one OptimizeApplication
/// drive both moldft and nemo.
///
/// This replaces the earlier `SCFTarget`, which built a fresh Application per
/// geometry through a factory and read a result schema (`res.at("energy")` /
/// `res.at("gradient")`) that no Application ever produced. Driving the engine
/// directly is both correct and much cheaper: one solve per geometry rather than a
/// full Application run with its own directory and checkpointing.
template <typename Engine> struct EngineOptTarget {
  World &world;
  Engine &engine;
  double last_energy = 0.0;
  Tensor<double> last_gradient;

  EngineOptTarget(World &w, Engine &e) : world(w), engine(e) {}

  /// Called by MolOpt::line_search.
  double value(const Tensor<double> &x) { return engine.value(x); }

  /// Called by MolOpt at each geometry. value() runs first so the engine is solved
  /// at this geometry before the gradient is taken; both engines treat a repeated
  /// geometry as a no-op, so this costs no extra solve.
  void energy_and_gradient(const Molecule &mol, double &energy,
                           Tensor<double> &gradient) {
    const Tensor<double> x = mol.get_all_coords().flat();
    energy = engine.value(x);
    gradient = engine.gradient(x);
    last_energy = energy;
    last_gradient = gradient;
  }
};

} // namespace madness
