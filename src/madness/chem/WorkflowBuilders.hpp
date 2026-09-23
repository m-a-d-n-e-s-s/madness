#pragma once

#include <algorithm>
#include <array>
#include <string>
#include <string_view>

#include <CCLib.hpp>
#include <Drivers.hpp>
#include <MoldftLib.hpp>
#include <ParameterManager.hpp>
#include <madness/mra/funcdefaults.h>

namespace madness::workflow_builders {

enum class WorkflowKind {
  Scf,
  Nemo,
  Response,
  Mp2Cc2,
  Cis,
  Oep,
  Optimize,
  Unknown
};

// `optimize` is deliberately NOT here: it is not a workflow but a flag on one --
// `--optimize --wf=<scf|nemo>`. The name still maps to WorkflowKind::Optimize so
// `--wf=optimize` can be answered with a migration message rather than "unknown".
inline constexpr std::array<const char *, 7> runnable_workflows = {
    "scf", "nemo", "response", "mp2", "cc2", "cis", "oep"};

inline WorkflowKind workflow_kind_from_name(std::string_view user_workflow) {
  if (user_workflow == "scf")
    return WorkflowKind::Scf;
  if (user_workflow == "nemo")
    return WorkflowKind::Nemo;
  if (user_workflow == "response")
    return WorkflowKind::Response;
  if (user_workflow == "mp2" || user_workflow == "cc2")
    return WorkflowKind::Mp2Cc2;
  if (user_workflow == "cis")
    return WorkflowKind::Cis;
  if (user_workflow == "oep")
    return WorkflowKind::Oep;
  if (user_workflow == "optimize")
    return WorkflowKind::Optimize;
  return WorkflowKind::Unknown;
}

inline const std::string &runnable_workflow_list() {
  static const std::string list = []() {
    std::string out;
    for (size_t i = 0; i < runnable_workflows.size(); ++i) {
      if (i > 0)
        out += ", ";
      out += runnable_workflows[i];
    }
    return out;
  }();
  return list;
}

inline void add_scf_workflow_drivers(World &world, Params &pm,
                                     qcapp::Workflow &wf) {
  auto reference = std::make_shared<SCFApplication<moldft_lib>>(world, pm);
  wf.addDriver(std::make_unique<qcapp::SinglePointDriver>(reference));
}

inline void add_nemo_workflow_drivers(World &world, Params &pm,
                                      qcapp::Workflow &wf) {
  // NB: no set_derived_value("k", 8) here. Pinning k defeated the thresh->k
  // derivation in SCF::set_protocol, so a nemo workflow ran the whole protocol
  // at one polynomial order. That was a workaround for nemo not surviving a k
  // change -- the nemos, the AO basis and the nuclear correlation factor were
  // not reprojected -- which NemoBase::set_protocol now handles. A deck that
  // wants a fixed k still says so explicitly, and that wins over any derived
  // value. moldft's builder above sets no k either.
  auto reference = std::make_shared<SCFApplication<nemo_lib>>(world, pm);
  wf.addDriver(std::make_unique<qcapp::SinglePointDriver>(reference));
}

/// Geometry optimization of a workflow's reference geometry: `--optimize` plus
/// `--wf=<scf|nemo>`, which names the reference method. One task,
/// qcapp::OptimizeDriver, which as a Driver can later own the displaced sub-runs a
/// numerical gradient needs. This is now the ONLY way to optimize a geometry --
/// the in-SCF `dft gopt` form has been removed and its keyvals retired.
inline void add_optimize_workflow_drivers(World &world, Params &pm,
                                          WorkflowKind kind,
                                          qcapp::Workflow &wf) {
  // NB: nothing to clear on the dft group. `--optimize` used to set `dft.gopt`,
  // which armed the *other* optimizer, and this had to undo it before a later
  // SCF task in the same workflow optimized a second time. That coupling is gone.
  switch (kind) {
  case WorkflowKind::Scf:
    wf.addDriver(
        std::make_unique<qcapp::OptimizeDriver<moldft_lib>>(world, pm));
    break;
  case WorkflowKind::Nemo:
    pm.get<CalculationParameters>().set_derived_value("k", 8);
    wf.addDriver(std::make_unique<qcapp::OptimizeDriver<nemo_lib>>(world, pm));
    break;
  default:
    MADNESS_EXCEPTION("--optimize is only supported for --wf=scf and "
                      "--wf=nemo; optimizing a correlated, excited-state or "
                      "response reference needs the downstream step to honour "
                      "the optimized geometry first (ARCHITECTURE_ROADMAP "
                      "change 2)",
                      1);
  }
}

// NOTE (M1 decoupling, Stage 2): the response workflow no longer has a MADchem
// driver. The v2 engine (ResponseApplication<molresponse_lib> via
// MolresponseLib.hpp) was removed; the v3 engine lives in the APP layer
// (madqc.cpp builds ResponseApplication<molresponse_v3_lib>) because MADchem
// must not depend on molresponse_v3 (circular lib dependency — see
// molresponse_v3/madqc_adapter.hpp). The Response enum + name table stay so
// workflow_kind_from_name/runnable_workflows keep listing it.

inline void add_cc2_workflow_drivers(World &world, Params &pm,
                                     qcapp::Workflow &wf) {
  TensorType tt = TT_2D;
  FunctionDefaults<6>::set_tensor_type(tt);

  auto &calc_param = pm.get<CalculationParameters>();
  auto &cc_param = pm.get<CCParameters>();
  auto &molecule = pm.get<Molecule>();

  calc_param.set_derived_value("k", 5);
  calc_param.set_derived_value("print_level", 2);
  calc_param.set_derived_value("econv", cc_param.get<double>("thresh_6d") * 0.01);
  // Chained workflows need the ground-state archive on disk: it is what lets the
  // reference step's results be reused without leaving the downstream step with
  // an orbital-less engine (SCFApplication::run, NextAction::Ok -> reload). The
  // response builder in madqc.cpp does the same. A deck that sets `save 0`
  // explicitly still wins -- set_derived_value yields to user-defined values.
  calc_param.set_derived_value("save", true);

  calc_param.set_derived_values(molecule);
  cc_param.set_derived_values();

  auto reference = std::make_shared<SCFApplication<nemo_lib>>(world, pm);
  auto ref_calc = reference->calc();
  wf.addDriver(std::make_unique<qcapp::SinglePointDriver>(reference));
  wf.addDriver(std::make_unique<qcapp::SinglePointDriver>(
      std::make_unique<CC2Application>(world, pm, ref_calc)));
}

inline void add_cis_workflow_drivers(World &world, Params &pm,
                                     qcapp::Workflow &wf) {
  // see add_cc2_workflow_drivers: the chain needs the ground-state archive
  pm.get<CalculationParameters>().set_derived_value("save", true);
  auto reference = std::make_shared<SCFApplication<nemo_lib>>(world, pm);
  auto ref_calc = reference->calc();
  wf.addDriver(std::make_unique<qcapp::SinglePointDriver>(reference));
  wf.addDriver(std::make_unique<qcapp::SinglePointDriver>(
      std::make_unique<TDHFApplication>(world, pm, ref_calc)));
}

inline void add_oep_workflow_drivers(World &world, Params &pm,
                                     qcapp::Workflow &wf) {
  auto &cparam = pm.get<CalculationParameters>();
  auto convergence_crit =
      cparam.get<std::vector<std::string>>("convergence_criteria");
  if (std::find(convergence_crit.begin(), convergence_crit.end(),
                "each_energy") == convergence_crit.end()) {
    convergence_crit.emplace_back("each_energy");
  }
  cparam.set_derived_value("convergence_criteria", convergence_crit);
  // see add_cc2_workflow_drivers: the chain needs the ground-state archive
  cparam.set_derived_value("save", true);

  auto reference = std::make_shared<SCFApplication<nemo_lib>>(world, pm);
  auto ref_calc = reference->calc();
  wf.addDriver(std::make_unique<qcapp::SinglePointDriver>(reference));
  wf.addDriver(std::make_unique<qcapp::SinglePointDriver>(
      std::make_unique<OEPApplication>(world, pm, ref_calc)));
}

/// @param optimize  `--optimize` was given: optimize this workflow's reference
///                  geometry instead of running the workflow itself.
inline void add_workflow_drivers(World &world, Params &pm,
                                 const std::string &user_workflow,
                                 qcapp::Workflow &wf, bool optimize = false) {
  const WorkflowKind kind = workflow_kind_from_name(user_workflow);
  if (optimize) {
    if (kind == WorkflowKind::Unknown)
      MADNESS_EXCEPTION("--optimize with an unknown workflow; use --wf=scf or "
                        "--wf=nemo",
                        1);
    add_optimize_workflow_drivers(world, pm, kind, wf);
    return;
  }
  switch (kind) {
  case WorkflowKind::Scf:
    add_scf_workflow_drivers(world, pm, wf);
    break;
  case WorkflowKind::Nemo:
    add_nemo_workflow_drivers(world, pm, wf);
    break;
  case WorkflowKind::Response:
    // App-layer workflow (engine v3 in madqc.cpp); MADchem hosts no response
    // driver since the v2 engine's removal (M1 Stage 2).
    throw std::runtime_error(
        "response workflow: no MADchem driver — the v2 engine was removed; "
        "run it through madqc (engine v3), which builds the v3 driver in the "
        "app layer");
  case WorkflowKind::Mp2Cc2:
    add_cc2_workflow_drivers(world, pm, wf);
    break;
  case WorkflowKind::Cis:
    add_cis_workflow_drivers(world, pm, wf);
    break;
  case WorkflowKind::Oep:
    add_oep_workflow_drivers(world, pm, wf);
    break;
  case WorkflowKind::Optimize:
    MADNESS_EXCEPTION("`--wf=optimize` is not a workflow: the workflow names the "
                      "reference method and --optimize asks for its geometry to "
                      "be optimized. Use `--optimize --wf=scf` or "
                      "`--optimize --wf=nemo`",
                      1);
    break;
  case WorkflowKind::Unknown:
  default: {
    std::string msg =
        "Unknown workflow: " + user_workflow + "\nAvailable workflows are: " +
        runnable_workflow_list();
    MADNESS_EXCEPTION(msg.c_str(), 1);
    break;
  }
  }
}

} // namespace madness::workflow_builders
