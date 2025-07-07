#include <QCCalculationParametersBase.h>
#include <madness/mra/mra.h>
#include <molecule.h>
#include <TDHF.h>

#include <iostream>
#include <string>

// #include "ParameterManager.hpp"
#include <madness/chem/ParameterManager.hpp>

// Define a concrete aliased ParameterManager type
using MyParamMgr = ParameterManager<CalculationParameters, ResponseParameters, OptimizationParameters, Molecule>;

int main(int argc, char **argv) {
  World &world = initialize(argc, argv);
  {
    startup(world, argc, argv, true);
    if (argc != 2) {
      if (world.rank() == 0) std::cerr << "Usage: molresponse2 [input_file.json]\n";
      finalize();
      return 1;
    }
    commandlineparser parser(argc, argv);

    if (argc != 2) {
      std::cerr << "Usage: " << argv[0] << " <input_file>\n";
      return 1;
    }
    std::string input_file = argv[1];

    // Construct the manager, reading .inp or JSON as needed
    MyParamMgr pm(world, parser);

    // Print out all groups
    if (world.rank() == 0) {
      pm.print_all();
    }

    // Dump merged JSON to stdout
    auto &all_json = pm.getAllInputJson();
    if (world.rank() == 0) {
      std::cout << "Merged JSON:\n";
      std::cout << all_json.dump(4) << std::endl;
    }
    finalize();
  }

  return 0;
}
