#pragma once
#include <madness/chem/SCF.h>

#include <madness/chem/InputWriter.hpp>
#include <madness/chem/ParameterManager.hpp>
#if defined(HAVE_SYS_TYPES_H) && defined(HAVE_SYS_STAT_H) && defined(HAVE_UNISTD_H)

#include <sys/stat.h>
#include <sys/types.h>
#include <unistd.h>

#endif

struct moldft_lib {
  static constexpr char const *label() { return "moldft"; }

  using Engine = SCF;

  // expose the live engine
  std::shared_ptr<Engine> engine(World &world, const Params &params) {
    if (!calc_)
      initialize_(world, params); // create once
    return calc_;
  }

  // params get's changed by SCF constructor
  inline nlohmann::json run(World &world, const Params &params) {
    const auto moldft_params = params.get<CalculationParameters>();
    const auto &molecule = params.get<Molecule>();

    auto scf = engine(world, params);

    // redirect any log files into outdir if needed…
    // Warm and fuzzy for the user
    if (world.rank() == 0) {
      print("\n\n");
      print(" MADNESS Hartree-Fock and Density Functional Theory Program");
      print(" ----------------------------------------------------------\n");
      print("\n");
      scf->molecule.print();
      print("\n");
      scf->param.print("dft");
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
    cr.set_converged_thresh(scf->converged_for_thresh);
    cr.set_converged_dconv(scf->converged_for_dconv);
    PropertyResults pr;
    pr.energy = energy;
    pr.dipole = dip;
    pr.gradient = grad;
    SCFResults sr;
    sr.aeps = scf->aeps;
    sr.beps = scf->beps;
    sr.properties = pr;

    nlohmann::json results = sr.to_json();
    results["convergence_info"] = cr.to_json();
    return results;
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

  std::shared_ptr<Engine> calc_;
}; // namespace moldft_lib

struct nemo_lib {
  using Engine = Nemo;
  static constexpr char const *label() { return "nemo"; }

  std::shared_ptr<Engine> engine(World &world, const Params &params) {
    if (!nemo_)
      initialize_(world, params);
    return nemo_;
  }

  nlohmann::json run(World &world, const Params &params) {

    auto nm = engine(world, params);

    nm->value();
    PropertyResults pr = nm->analyze();

    ConvergenceResults cr;
    cr.set_converged_thresh(nm->get_calc()->converged_for_thresh);
    cr.set_converged_dconv(nm->get_calc()->converged_for_dconv);

    SCFResults sr;
    sr.aeps = nm->get_calc()->aeps;
    sr.beps = nm->get_calc()->beps;
    sr.properties = pr;
    sr.scf_total_energy = nm->get_calc()->current_energy;

    nlohmann::json results = sr.to_json();
    results["convergence_info"] = cr.to_json();
    return results;
  }

private:
  void initialize_(World &world, const Params &params) {
    nemo_ = std::make_shared<Nemo>(world, params.get<CalculationParameters>(),
                                   params.get<Nemo::NemoCalculationParameters>(), params.get<Molecule>());
  }

  std::shared_ptr<Engine> nemo_;
};
