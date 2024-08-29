
#include <madness/chem/SCF.h>
#include <madness/chem/molopt.h>
#include <madness/misc/info.h>
#include <madness/world/worldmem.h>
#include "madness/mra/commandlineparser.h"

int main(int argc, char** argv) {
  World& world = initialize(argc, argv);
  SCF calc(world, parser);
  MolOpt opt(calc.param.gmaxiter(), 0.1, calc.param.gval(), calc.param.gtol(),
             1e-3,  //XTOL
             1e-5,  //EPREC
             calc.param.gprec(),
             (world.rank() == 0) ? 1 : 0,  //print_level
             calc.param.algopt());

  MolecularEnergy target(world, calc);
  opt.optimize(calc.molecule, target);
}
else {
  MolecularEnergy E(world, calc);
}
}
