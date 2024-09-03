/*
  This file is part of MADNESS.

  Copyright (C) 2007,2010 Oak Ridge National Laboratory

  This program is free software; you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation; either version 2 of the License, or
  (at your option) any later version.

  This program is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
  GNU General Public License for more details.

  You should have received a copy of the GNU General Public License
  along with this program; if not, write to the Free Software
  Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307 USA

  For more information please contact:

  Robert J. Harrison
  Oak Ridge National Laboratory
  One Bethel Valley Road
  P.O. Box 2008, MS-6367

  email: harrisonrj@ornl.gov
  tel:   865-241-3937
  fax:   865-572-0680
*/
#ifndef SRC_APPS_MADQC_OPT_STRATEGIES_HPP_
#define SRC_APPS_MADQC_OPT_STRATEGIES_HPP_

#include "parameter_manager.hpp"
#include <madness/chem/SCF.h>
#include <madness/chem/molecule.h>
#include <madness/chem/molopt.h>
#include <madness/misc/info.h>
#include <madness/mra/commandlineparser.h>
#include <madness/world/worldmem.h>
#include <madqc/utils.hpp>

using namespace madness;

class OptimizationStrategy {
 public:
  virtual ~OptimizationStrategy() = default;
  virtual Molecule optimize(World& world, const ParameterManager&) = 0;
};

class MoldftOptStrategy : public OptimizationStrategy {
 public:
  Molecule optimize(World& world, const ParameterManager& pm) override {
    auto& moldft_params = pm.get_moldft_params();
    auto& opt_params = pm.get_optimization_params();
    auto& molecule = pm.get_molecule();

    // To create the calc object, we need to write the input file for the parser
    if (world.rank() == 0) {
      json moldft_input_json = {};
      moldft_input_json["dft"] = moldft_params.to_json_if_precedence("defined");
      moldft_input_json["molecule"] = molecule.to_json();
      std::ofstream ofs("moldft.in");
      write_moldft_input(moldft_input_json, ofs);
      ofs.close();
    }
    world.gop.fence();

    commandlineparser parser;
    parser.set_keyval("input", "moldft.in");
    if (world.rank() == 0)
      ::print("input filename: ", parser.value("input"));
    print_meminfo(world.rank(), "startup");
    FunctionDefaults<3>::set_pmap(pmapT(new LevelPmap<Key<3>>(world)));

    std::cout.precision(6);
    SCF calc(world, parser);

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
    calc.set_protocol<3>(world, calc.param.protocol()[0]);

    MolOpt opt(opt_params.get_maxiter(),             // geoometry max iter
               0.1,                                  // geometry step size
               opt_params.get_value_precision(),     // value precision
               opt_params.get_geometry_tolerence(),  // geometry tolerance
               1e-3,                                 // XTOL
               1e-5,                                 // EPREC
               opt_params.get_gradient_precision(),  // gradient precision
               (world.rank() == 0) ? 1 : 0,          // print_level
               calc.param.algopt());                 // algorithm options

    MolecularEnergy target(world, calc);

    return opt.optimize(calc.molecule, target);
  }
};

#endif  //  SRC_APPS_MADQC_OPT_STRATEGIES_HPP_
