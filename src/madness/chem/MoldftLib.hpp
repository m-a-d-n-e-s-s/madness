#pragma once
#include <madness/chem/SCF.h>
#include <madness/chem/InputWriter.hpp>
#include <madness/chem/ParameterManager.hpp>
#if defined(HAVE_SYS_TYPES_H) && defined(HAVE_SYS_STAT_H) && \
    defined(HAVE_UNISTD_H)

#include <sys/stat.h>
#include <sys/types.h>
#include <unistd.h>

#endif

struct moldft_lib {
  struct Results {
    double energy;
    tensorT dipole;
    std::optional<Tensor<double>> gradient;
  };

  static constexpr char const* label() { return "moldft"; }

  static Results run_nemo(std::shared_ptr<Nemo> nemo) {
    Results results;
    nemo->get_param().print("dft");
    results.energy = nemo->value();
    return results;
  }

  // params get's changed by SCF constructor
  inline static Results run_scf(World& world, Params& params,
                                const std::filesystem::path& outdir) {
    auto moldft_params = params.get<CalculationParameters>();
    const auto& molecule = params.get<Molecule>();

    auto archive_name = moldft_params.prefix() + ".restartdata";
    auto restart_path = path(archive_name) / ".00000";

    if (std::filesystem::exists(restart_path)) {
      moldft_params.set_user_defined_value<bool>("restart", true);
    }
    world.gop.fence();
    // save the input file so we can read it back in and get derived parameters
    if (world.rank() == 0) {
      json moldft_input_json = {};
      moldft_input_json["dft"] = moldft_params.to_json_if_precedence("defined");
      moldft_input_json["molecule"] = molecule.to_json_if_precedence("defined");
      print("moldft_input_json: ", moldft_input_json.dump(4));
      std::ofstream ofs("moldft.in");
      write_moldft_input(moldft_input_json, ofs);
      ofs.close();
    }
    world.gop.fence();
    commandlineparser parser;
    parser.set_keyval("input", "moldft.in");
    if (world.rank() == 0) ::print("input filename: ", parser.value("input"));

    FunctionDefaults<3>::set_pmap(pmapT(new LevelPmap<Key<3>>(world)));

    std::cout.precision(6);
    SCF calc(world, parser);
    // SCF& calc=dynamic_cast<SCF&>(*this);
    // read moldft.in file used

    // redirect any log files into outdir if needed…
    // Warm and fuzzy for the user
    if (world.rank() == 0) {
      print("\n\n");
      print(" MADNESS Hartree-Fock and Density Functional Theory Program");
      print(" ----------------------------------------------------------\n");
      print("\n");
      calc.molecule.print();
      print("\n");
      calc.param.print("dft");
    }
    // Come up with an initial OK data map
    if (world.size() > 1) {
      calc.set_protocol<3>(world, 1e-4);
      calc.make_nuclear_potential(world);
      calc.initial_load_bal(world);
    }
    // vama
    calc.set_protocol<3>(world, calc.param.protocol()[0]);

    MolecularEnergy E(world, calc);
    double energy = E.value(calc.molecule.get_all_coords().flat());
    if (world.rank() == 0 && calc.param.print_level() > 0)
      E.output_calc_info_schema();

    functionT rho = calc.make_density(world, calc.aocc, calc.amo);
    functionT brho = rho;
    if (calc.param.nbeta() != 0 && !calc.param.spin_restricted())
      brho = calc.make_density(world, calc.bocc, calc.bmo);
    rho.gaxpy(1.0, brho, 1.0);

    // optionally compute gradient, dipole, etc.
    Tensor<double> grad;
    if (calc.param.derivatives()) {
      auto g = calc.derivatives(world, rho);
      calc.e_data.add_gradient(g);
      grad = std::move(g);
      // write gradient file into outdir if you want
    }

    tensorT dip;
    if (calc.param.dipole()) {
      auto d =
          calc.dipole(world, calc.make_density(world, calc.aocc, calc.amo));
      dip = std::move(d);
    }
    calc.do_plots(world);

    // return structured results
    return {energy, dip,
            grad.has_data() ? std::nullopt : std::make_optional(grad)};
  }
};  // namespace moldft_lib
