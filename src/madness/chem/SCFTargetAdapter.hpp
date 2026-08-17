// SCFTargetAdapter.hpp
#pragma once

#include <madness/chem/molecule.h>
#include <madness/tensor/solvers.h>
#include <madness/tensor/tensor.h>

namespace madness {

/// What MolOpt needs at each geometry, as a runtime-selectable interface.
///
/// MolOpt is a template over its target, so it binds this by reference and calls
/// through the vtable -- immaterial next to an SCF solve. The point of the
/// indirection is that the *source* of the gradient is a choice:
///
///  - `AnalyticTarget` below: the engine's own analytic derivatives.
///  - numerical gradients (roadmap change 4): energies at displaced geometries,
///    each computed by its own sub-run under `task_<i>/step_<k>/disp_<...>/` with
///    its own calc_info, so displacements are restartable and can be spread
///    across subworlds. That implementation belongs to `OptimizeDriver`, which as
///    a Driver may own sub-runs; it plugs in here without touching MolOpt.
struct GeometryTarget {
  virtual ~GeometryTarget() = default;

  /// Called by MolOpt::line_search.
  virtual double value(const Tensor<double> &x) = 0;

  /// Called by MolOpt at each geometry.
  virtual void energy_and_gradient(const Molecule &mol, double &energy,
                                   Tensor<double> &gradient) = 0;
};

/// Analytic gradients straight from the engine.
///
/// Both reference engines already provide the interface this needs:
/// `MolecularEnergy` (SCF.h) wraps an SCF, and `Nemo` is one itself via
/// `NemoBase`. So one implementation covers moldft and nemo, and the engines need
/// no adapter of their own.
class AnalyticTarget : public GeometryTarget {
public:
  explicit AnalyticTarget(OptimizationTargetInterface &engine)
      : engine_(engine) {}

  double value(const Tensor<double> &x) override { return engine_.value(x); }

  void energy_and_gradient(const Molecule &mol, double &energy,
                           Tensor<double> &gradient) override {
    // value_and_gradient is the interface's own both-at-once hook; its default
    // solves at x and then differentiates, and a repeated geometry is a no-op in
    // both engines, so this costs one solve.
    engine_.value_and_gradient(mol.get_all_coords().flat(), energy, gradient);
  }

private:
  OptimizationTargetInterface &engine_;
};

} // namespace madness
