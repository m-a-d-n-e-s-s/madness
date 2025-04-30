#include <fstream>
#include <iostream>
#include <string>

// Include your parameter manager and group definitions
#include <../molresponse/response_parameters.h>
#include <CalculationParameters.h>
#include <ParameterManager.hpp>
#include <QCCalculationParametersBase.h>
#include <molecule.h>

// Define a concrete aliased ParameterManager type
using MyParamMgr =
    ParameterManager<CalculationParameters, ResponseParameters,
                     OptimizationParameters, Molecule::GeometryParameters>;

int main(int argc, char *argv[]) {
  // Initialize MADNESS world (passes MPI args, etc.)
  World &world = madness::initialize(argc, argv);
  {

    startup(world, argc, argv, true);
    if (world.rank() == 0) {
      print_header1("MOLRESPONSE -- MADNESS Time-Dependent Density Functional "
                    "Theory Excited-State Program ");
    }

    if (argc != 2) {
      std::cerr << "Usage: " << argv[0] << " <input_file>\n";
      return 1;
    }
    std::string input_file = argv[1];

    // Construct the manager, reading .inp or JSON as needed
    MyParamMgr pm(world, input_file);

    // Print out all groups
    if (world.rank() == 0) {
      pm.print_all();
    }

    // Dump merged JSON to stdout
    auto all_json = pm.getAllInputJson();
    if (world.rank() == 0) {
      std::cout << "Merged JSON:\n";
      std::cout << all_json.dump(4) << std::endl;
    }
    finalize();
  }

  return 0;
}
