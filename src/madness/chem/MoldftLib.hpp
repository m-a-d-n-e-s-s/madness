#pragma once
#include <madness/chem/SCF.h>

#include <madness/chem/InputWriter.hpp>
#include <madness/chem/ParameterManager.hpp>
#if defined(HAVE_SYS_TYPES_H) && defined(HAVE_SYS_STAT_H) && defined(HAVE_UNISTD_H)

#include <sys/stat.h>
#include <sys/types.h>
#include <unistd.h>

#endif

inline NextAction decide_next_action(bool at_protocol, bool archive_needed, bool archive_exists,
                                     bool all_properties_computed, bool restart_exists) {
  // We must recompute if any of these are true:
  const bool must_redo = !at_protocol                           // not at final protocol
                         || (archive_needed && !archive_exists) // user wants archive but it's missing
                         || !all_properties_computed;           // a requested prop is missing

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

template <typename SCFParams> NextAction valid(const SCFResultsTuple &results, const SCFParams &params) {
  // Take a copy of the parameters
  auto [sr, pr, cr] = results;

  // Required convergence for "final" protocol
  const auto vthresh = params.econv();
  const auto vdconv = params.dconv();
  const bool archive_needed = params.save();

  // Requested outputs
  const bool need_energy = true;
  const bool need_dipole = params.dipole();
  const bool need_gradient = params.derivatives();

  // Files/paths
  const std::string archivename = params.prefix();
  const auto restart_path = std::filesystem::path(archivename).replace_extension("restartdata.00000");
  const bool archive_exists = std::filesystem::exists(restart_path);

  // State in results
  // TODO: It's hard to be certain what converged_for_thresh means.  I thought it was the final protocol.  Turns out
  // it's actaully set in SCF.cc as param.econv(), which says nothing about the threshold refinement.
  //
  const bool at_protocol = (cr.converged_for_thresh == vthresh && cr.converged_for_dconv == vdconv);

  const auto pjson = pr.to_json();
  const bool energy_ok = pjson.contains("energy");
  const bool dipole_ok = pjson.contains("dipole");
  const bool gradient_ok = pjson.contains("gradient");

  // Only require props the user asked for
  const bool all_properties_computed =
      (need_energy ? energy_ok : true) && (need_dipole ? dipole_ok : true) && (need_gradient ? gradient_ok : true);

  // Decide action
  const bool must_redo = !at_protocol || (archive_needed && !archive_exists) || !all_properties_computed;

  // if we don't need to redo, we can either reload or return ok
  if (!must_redo)
    return (!archive_exists && !archive_needed) ? NextAction::ReloadOnly : NextAction::Ok;
  // with we need to redo we can restart from the exisiting archive
  return archive_exists ? NextAction::Restart : NextAction::Redo;
}

struct moldft_lib {
  static constexpr char const *label() { return "moldft"; }
  NextAction next_action_ = NextAction::Ok;

  vector<double> protocol;
  SCFResultsTuple last_results_;

  using Calc = SCF;

  NextAction valid(const SCFResultsTuple &results, const Params &params) {
    last_results_ = results;
    return ::valid(results, params.get<CalculationParameters>());
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

    MolecularEnergy E(world, *scf);
    double energy = E.value(scf->molecule.get_all_coords().flat());
    if (world.rank() == 0 && scf->param.print_level() > 0)
      E.output_calc_info_schema();

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
    SCFResults sr;
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
      if (!std::filesystem::exists("mad.in")) {
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

    FunctionDefaults<3>::set_pmap(pmapT(new LevelPmap<Key<3>>(world)));
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
  NextAction valid(const SCFResultsTuple &results, const Params &params) {
    // Take a copy of the parameters
    return ::valid(results, params.get<CalculationParameters>());
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
