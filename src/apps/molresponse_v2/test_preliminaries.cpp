#include "../madqc/parameter_manager.hpp"
#include "GroundStateData.hpp"
#include "ResponseManager.hpp"
#include "madness/world/world.h"
#include <cassert>
#include <cmath>
#include <iostream>
#include <string>
#include <tensor_json.hpp>

using namespace madness;

std::pair<Tensor<double>, Tensor<double>>
get_fock_from_json(const json &fock_json, double thresh, int k,
                   bool spinrestricted) {
  // fock_a and fock_b if spinrestricted
  std::pair<Tensor<double>, Tensor<double>> fock_data;
  // write the fock matrix to a file
  auto protocol = std::string("thresh: ") +
                  std::to_string(FunctionDefaults<3>::get_thresh());
  protocol +=
      std::string(" k: ") + std::to_string(FunctionDefaults<3>::get_k());

  // Check if the protocol exists in the JSON
  if (fock_json.find(protocol) == fock_json.end()) {
    std::cerr << "Protocol not found in JSON: " << protocol << std::endl;
    return fock_data;
  } else {

    fock_data.first = tensor_from_json<double>(fock_json[protocol]["focka"]);
    if (!spinrestricted) {
      fock_data.second = tensor_from_json<double>(fock_json[protocol]["fockb"]);
    }
    return fock_data;
  }
}

// Utility function to compare values within tolerance
bool approximately_equal(double computed, double reference, double tol = 1e-6) {
  return std::abs(computed - reference) <= tol;
}

// Compare two tensors element-wise with a tolerance threshold
bool tensor_approx_equal(World &world, const Tensor<double> &computed,
                         const Tensor<double> &reference, double tol = 1e-6,
                         bool verbose = true) {

  // dims are pointers to array of longs

  auto check_dims = [](const long *a, const long *b, long size) {
    for (int i = 0; i < size; ++i) {
      if (a[i] != b[i]) {
        return false;
      }
    }
    return true;
  };

  auto computed_dims = computed.dims();
  auto reference_dims = reference.dims();
  if (!check_dims(computed_dims, reference_dims, computed.ndim())) {
    if (verbose and world.rank() == 0)
      std::cerr << "Tensor dimension mismatch: computed = " << computed.dims()
                << ", reference = " << reference.dims() << "\n";
    return false;
  }

  auto flat_computed = computed.flat();
  auto flat_reference = reference.flat();

  for (int i = 0; i < flat_computed.size(); ++i) {
    double c_val = flat_computed[i];
    double r_val = flat_reference[i];
    if (std::abs(c_val - r_val) > tol) {
      if (verbose) {
        if (world.rank() == 0) {
          std::cerr << "Tensor mismatch at index " << i
                    << ": computed = " << c_val << ", reference = " << r_val
                    << ", diff = " << std::abs(c_val - r_val)
                    << " > tol = " << tol << "\n";
        }
      }
      return false;
    }
  }

  return true;
}

int main(int argc, char **argv) {
  World &world = initialize(argc, argv);
  {

    startup(world, argc, argv, true);
    if (world.rank() == 0)
      print(info::print_revision_information());

    if (argc != 3) {
      if (world.rank() == 0)
        std::cerr
            << "Usage: test_preliminaries [input_file.json] [reference.txt]\n";
      finalize();
      return 1;
    }

    const std::string input_filename = argv[1];
    const std::string reference_filename = argv[2];
    commandlineparser parser(argc, argv);
    ParameterManager params;
    path input_file(argv[1]);
    params = ParameterManager(world, input_file);
    // Extract relevant parameters
    const std::string ground_state_archive = "moldft.restartdata";
    const std::string ground_state_fock_archive = "moldft.fock.json";
    // Read moldft.fock.json
    json fock_archive;

    std::ifstream fock_file(ground_state_fock_archive);
    if (!fock_file.is_open()) {
      if (world.rank() == 0)
        std::cerr << "Error opening file: " << ground_state_fock_archive
                  << std::endl;
      finalize();
      return 1;
    } else {
      fock_file >> fock_archive;
      fock_file.close();
    }

    // Initialize the ResponseManager with ground-state archive
    auto ground_state =
        GroundStateData(world, ground_state_archive, params.get_molecule());
    ResponseManager rm =
        ResponseManager(world, params.get_molresponse_params());

    auto protocol = params.get_moldft_params().protocol();
    FunctionDefaults<3>::set_pmap(pmapT(new LevelPmap<Key<3>>(world)));

    std::pair<Tensor<double>, Tensor<double>> fock_data;

    for (auto thresh : protocol) {
      rm.setProtocol(world, ground_state.getL(), thresh);
      if (world.rank() == 0) {
        FunctionDefaults<3>::print();
      }
      auto k = FunctionDefaults<3>::get_k();
      fock_data = get_fock_from_json(fock_archive, thresh, k,
                                     ground_state.isSpinRestricted());

      ground_state.prepareOrbitals(world, k, thresh);
      ground_state.computePreliminaries(
          world, *rm.getCoulombOp(), rm.getVtol(),
          params.get_molresponse_params().fock_json_file());

      if (world.rank() == 0) {
        print(ground_state.Hamiltonian);
      }

      // Compare the computed Hamiltonian with the reference data
      double test_thresh = thresh * .001;
      bool approx_equal = tensor_approx_equal(
          world, ground_state.Hamiltonian, fock_data.first, test_thresh, true);
      while (!approx_equal) {
        test_thresh *= 10;
        approx_equal = tensor_approx_equal(world, ground_state.Hamiltonian,
                                           fock_data.first, test_thresh, true);
      }
      if (world.rank() == 0) {
        if (approx_equal) {
          print("Hamiltonian matches reference data within tolerance.");
          print("Threshold used for comparison: ", test_thresh);

        } else {
          std::cerr << "Hamiltonian does not match reference data.\n";
          std::cerr << "Threshold used for comparison: " << test_thresh << "\n";
        }
      }
    }

    // auto prelims = response_manager.currentPreliminaries();
    world.gop.fence();
    world.gop.fence();
    print_stats(world);
  }
  finalize();
  return 0;
}
