#pragma once
#include "MolecularProperty.hpp"
#include "ResponseIO.hpp"
#include "ResponseVector.hpp"
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <madness/external/nlohmann_json/json.hpp>
#include <madness/tensor/tensor.h>
#include <madness/tensor/tensor_json.hpp>
#include <string>

namespace fs = std::filesystem;
using json = nlohmann::json;

class PropertyManager {
public:
  explicit PropertyManager(const std::string &filename) : filename_(filename) {
    if (fs::exists(filename_)) {
      std::ifstream in(filename_);
      in >> data_;
    } else {
      data_["polarizability"] = json::object();
      data_["hyperpolarizability"] = json::object();
    }
  }

  void save() const {
    std::ofstream out(filename_);
    out << std::setw(2) << data_ << "\n";
  }

  // Store full polarizability tensor at frequency omega
  void set_alpha(double omega, const madness::Tensor<double> &tensor,
                 const std::string &dirs) {
    data_["polarizability"][freq_str(omega)] = tensor_to_json(tensor);
    data_["polarizability"][freq_str(omega)]["directions"] = dirs;
  }

  bool has_alpha(double omega, std::string dir) const {
    // also check if they are the same size
    auto saved_dir = data_["polarizability"][freq_str(omega)]["directions"]
                         .get<std::string>();

    if (saved_dir.size() != dir.size()) {
      return false;
    }

    std::string omega_str = freq_str(omega);
    return data_["polarizability"].contains(omega_str) &&
           !data_["polarizability"][omega_str].is_null();
  }

  // Store a 3-element beta result for a given BC input pair (e.g., xx) at ω1,
  // ω2
  void set_beta(double omega1, double omega2, const std::string &bc,
                const madness::Tensor<double> &tensor) {
    data_["hyperpolarizability"][freq_str(omega1)][freq_str(omega2)][bc] =
        tensor_to_json(tensor);
  }
  void set_beta_dirs(const std::string &dirs) {
    data_["hyperpolarizability"]["directions"] = dirs;
  }

  bool has_beta(double omega1, double omega2, const std::string &bc) const {
    auto &beta = data_["hyperpolarizability"];
    return beta.contains(freq_str(omega1)) &&
           beta[freq_str(omega1)].contains(freq_str(omega2)) &&
           beta[freq_str(omega1)][freq_str(omega2)].contains(bc) &&
           !beta[freq_str(omega1)][freq_str(omega2)][bc].is_null();
  }

  void print_alpha_table() const {
    if (!data_.contains("polarizability")) {
      std::cout << "No polarizability data found.\n";
      return;
    }

    std::cout << "\n📐 Polarizability (α) Tensor Components:\n";

    // Build header based on available directions
    for (const auto &[freq_str, tensor_json] :
         data_["polarizability"].items()) {
      if (!tensor_json.contains("directions"))
        continue;

      const std::string dirs = tensor_json["directions"];
      const madness::Tensor<double> alpha =
          tensor_from_json<double>(tensor_json);

      // Print header
      std::cout << "\nω = " << freq_str << "\n";
      std::cout << std::setw(10) << "";
      for (char d : dirs)
        std::cout << std::setw(12) << "α_" + std::string(1, d);
      std::cout << "\n";

      // Print values row by row
      bool has_values = alpha.has_data();

      for (size_t i = 0; i < dirs.size(); ++i) {
        std::cout << std::setw(10) << dirs[i];
        for (size_t j = 0; j < dirs.size(); ++j) {
          if (i < alpha.dim(0) && j < alpha.dim(1) && has_values) {
            std::cout << std::setw(12) << std::fixed << std::setprecision(6)
                      << alpha(i, j);
          } else
            std::cout << std::setw(12) << "--";
        }
        std::cout << "\n";
      }
    }
  }

  void print_beta_table() const {
    if (!data_.contains("hyperpolarizability")) {
      std::cout << "No hyperpolarizability data found.\n";
      return;
    }
    std::cout << "\n🔺 Hyperpolarizability (β) Components:\n";

    // Get available directions
    std::string dirs = "";
    if (data_["hyperpolarizability"].contains("directions"))
      dirs = data_["hyperpolarizability"]["directions"].get<std::string>();

    for (const auto &[w1_str, w1_entry] :
         data_["hyperpolarizability"].items()) {
      if (w1_str == "directions")
        continue; // skip metadata

      for (const auto &[w2_str, w2_entry] : w1_entry.items()) {
        std::cout << "\nω₁ = " << w1_str << ", ω₂ = " << w2_str << "\n";
        std::cout << std::setw(8) << "BC";
        for (char A : dirs)
          std::cout << std::setw(12) << "β_" + std::string(1, A) + "BC";
        std::cout << "\n";

        for (const auto &[bc, tensor_json] : w2_entry.items()) {
          madness::Tensor<double> beta = tensor_from_json<double>(tensor_json);
          std::cout << std::setw(8) << bc;

          for (int i = 0; i < dirs.size(); ++i) {
            if (i < beta.size())
              std::cout << std::setw(12) << std::fixed << std::setprecision(6)
                        << beta(i);
            else
              std::cout << std::setw(12) << "--";
          }

          std::cout << "\n";
        }
      }
    }
  }

  const json &json_data() const { return data_; }

private:
  std::string filename_;
  json data_;

  static std::string freq_str(double freq) {
    std::ostringstream ss;
    ss << std::scientific << std::setprecision(3) << freq;
    return ss.str();
  }
};

void initialize_property_structure(
    PropertyManager &pm, const std::vector<MolecularProperty> &props) {
  for (const auto &prop : props) {
    if (prop.type == MolecularPropertyType::Polarizability) {
      auto num_dirs = prop.directions.size();

      std::string directions_string;
      for (char dir : prop.directions) {
        directions_string += dir;
      }

      for (double omega : prop.frequencies) {
        if (!pm.has_alpha(omega, directions_string)) {
          madness::Tensor<double> empty_tensor(num_dirs, num_dirs,
                                               0.0); // 3x3 alpha
          pm.set_alpha(omega, empty_tensor, directions_string);
        }
      }
    } else if (prop.type == MolecularPropertyType::Hyperpolarizability) {
      const auto &dirs = prop.directions;
      const auto &freqs = prop.frequencies;
      const auto directions_string = std::string(dirs.begin(), dirs.end());
      pm.set_beta_dirs(directions_string);

      for (size_t b = 0; b < freqs.size(); ++b) {
        for (size_t c = b; c < freqs.size(); ++c) {
          double omega1 = freqs[b];
          double omega2 = freqs[c];
          for (char B : dirs) {
            for (char C : dirs) {
              std::string bc = std::string() + B + C;
              if (!pm.has_beta(omega1, omega2, bc)) {
                madness::Tensor<double> empty(3, 0.0); // x(xx), y(xx), z(xx)
                pm.set_beta(omega1, omega2, bc, empty);
              }
            }
          }
        }
      }
    }
  }

  pm.save();
}

void compute_alpha(World &world,
                   const std::map<std::string, ResponseState> &state,
                   const GroundStateData &gs,
                   const std::vector<double> &frequencies,
                   const std::string &directions, PropertyManager &pm) {

  size_t num_directions = directions.size();
  auto num_orbitals = gs.getNumOrbitals();
  auto num_frequencies = frequencies.size();

  if (world.rank() == 0) {
    std::cout << "Number of orbitals: " << num_orbitals << "\n";
    std::cout << "Number of frequencies: " << num_frequencies << "\n";
    std::cout << "Number of directions: " << num_directions << "\n";
  }

  std::vector<ResponseState> active_states(num_directions);
  std::vector<vector_real_function_3d> VP(num_directions);
  std::vector<ResponseVector> active_vec(num_directions);

  if (world.rank() == 0) {
    std::cout << "Computing polarizability for frequencies: ";
    for (const auto &freq : frequencies) {
      std::cout << freq << " ";
    }
    std::cout << "\n";
  }

  for (int d = 0; d < directions.size(); d++) {
    if (world.rank() == 0) {
      std::cout << "Computing polarizability for direction: " << directions[d]
                << "\n";
    }
    DipolePerturbation pert{directions[d]};

    ResponseState state(pert, PerturbationType::Dipole, frequencies,
                        {FunctionDefaults<3>::get_thresh()},
                        gs.isSpinRestricted());
    if (world.rank() == 0) {
      std::cout << "Perturbation: " << state.perturbationDescription() << "\n";
    }
    active_states[d] = state;
    VP[d] = state.perturbation_vector(world, gs);
  }

  for (int i = 0; i < num_frequencies; i++) {
    madness::Tensor<double> alpha(num_directions, num_directions);
    double alpha_factor = -2.0;
    auto omega = frequencies[i];
    if (pm.has_alpha(frequencies[i], directions)) {
      continue;
    }
    int num_rf = 0;
    if (omega == 0.0) {
      num_rf = num_orbitals;
      alpha_factor *= 2.0;
    } else {
      num_rf = 2 * num_orbitals;
      for (auto &vp_dir : VP) {
        vp_dir.insert(vp_dir.end(), vp_dir.begin(), vp_dir.end());
      }
    }
    for (int j = 0; j < num_directions; j++) {
      auto &active_state = active_states[j];
      if (world.rank() == 0) {
        std::cout << "Loading response vector for perturbation: "
                  << active_state.perturbationDescription() << "\n";
      }
      load_response_vector(world, num_orbitals, active_states[j], 0, i,
                           active_vec[j]);
    }
    world.gop.fence();
    for (int k = 0; k < num_rf; k++) {
      auto response_k = vector_real_function_3d(num_directions);
      auto vp_k = vector_real_function_3d(num_directions);
      for (int j = 0; j < num_directions; j++) {
        response_k[j] =
            std::visit([&](auto &v) { return v.flat[k]; }, active_vec[j]);
        vp_k[j] = VP[j][k];
      }
      world.gop.fence();
      alpha += matrix_inner(world, response_k, vp_k);
    }
    world.gop.fence();
    alpha *= alpha_factor;
    world.gop.fence();
    pm.set_alpha(omega, alpha, directions);
  }
  pm.save();
}


void compute_beta(World &world,
                  const std::map<std::string, ResponseState> &state,
                  const GroundStateData &gs,
                  const std::vector<double> &frequencies,
                  const std::string &directions, PropertyManager &pm) {

  size_t num_directions = directions.size();
  auto num_orbitals = gs.getNumOrbitals();
  auto num_frequencies = frequencies.size();

  if (world.rank() == 0) {
    std::cout << "Number of orbitals: " << num_orbitals << "\n";
    std::cout << "Number of frequencies: " << num_frequencies << "\n";
    std::cout << "Number of directions: " << num_directions << "\n";
  }

  std::vector<ResponseState> active_states(num_directions);
  std::vector<vector_real_function_3d> VP(num_directions);
  std::vector<ResponseVector> active_vec(num_directions);

  if (world.rank() == 0) {
    std::cout << "Computing polarizability for frequencies: ";
    for (const auto &freq : frequencies) {
      std::cout << freq << " ";
    }
    std::cout << "\n";
  }

  std::map<char, ResponseState> state_map;
  std::map<char, vector_real_function_3d> perturbation_map;
  for (int d = 0; d < directions.size(); d++) {
    if (world.rank() == 0) {
      std::cout << "Computing polarizability for direction: " << directions[d]
                << "\n";
    }
    DipolePerturbation pert{directions[d]};

    ResponseState state(pert, PerturbationType::Dipole, frequencies,
                        {FunctionDefaults<3>::get_thresh()},
                        gs.isSpinRestricted());
    if (world.rank() == 0) {
      std::cout << "Perturbation: " << state.perturbationDescription() << "\n";
    }
    state_map[directions[d]] = state;
    perturbation_map[directions[d]] = state.perturbation_vector(world, gs);
  }

  // all BC required during one frequency pair
  auto BC_pairs = std::vector<std::pair<char, char>>();
  for (int i = 0; i < num_directions; i++) {
    for (int j = i + 1; j < num_directions; j++) {
      BC_pairs.push_back({directions[i], directions[j]});
    }
  }

  vector<int> bfreq_index;
  vector<int> cfreq_index;
  vector<int> bc_index;

  for (int b = 0; b < num_frequencies; b++) {
    for (int c = b; c < num_frequencies; b++) {

      for (int bc = 0; bc < BC_pairs.size(); bc++) {
        bfreq_index.push_back(b);
        cfreq_index.push_back(c);
        bc_index.push_back(bc);
      }
    }
  }
  // Write VBC macro task to write a vector of bools indicating if it
  // successfully constructed it's VBC rhs and saved it to disk.
  

  // Then for each set of lhs, (X,Y,Z) load all "connecting" VBCs, and compute tensor
  // component.... example
  //
  // XYZ(2) connects VBC(0,2) and VBC(1,1)
  // XYZ(3) connects VBC(0,3) and VBC(1,2)
  //
  // A note, we need a map from frequency to frequency index,
  

  for (int i = 0; i < num_frequencies; i++) {
    madness::Tensor<double> alpha(num_directions, num_directions);
    double alpha_factor = -2.0;
    auto omega = frequencies[i];
    if (pm.has_alpha(frequencies[i], directions)) {
      continue;
    }
    int num_rf = 0;
    if (omega == 0.0) {
      num_rf = num_orbitals;
      alpha_factor *= 2.0;
    } else {
      num_rf = 2 * num_orbitals;
      for (auto &vp_dir : VP) {
        vp_dir.insert(vp_dir.end(), vp_dir.begin(), vp_dir.end());
      }
    }
    for (int j = 0; j < num_directions; j++) {
      auto &active_state = active_states[j];
      if (world.rank() == 0) {
        std::cout << "Loading response vector for perturbation: "
                  << active_state.perturbationDescription() << "\n";
      }
      load_response_vector(world, num_orbitals, active_states[j], 0, i,
                           active_vec[j]);
    }
    world.gop.fence();
    for (int k = 0; k < num_rf; k++) {
      auto response_k = vector_real_function_3d(num_directions);
      auto vp_k = vector_real_function_3d(num_directions);
      for (int j = 0; j < num_directions; j++) {
        response_k[j] =
            std::visit([&](auto &v) { return v.flat[k]; }, active_vec[j]);
        vp_k[j] = VP[j][k];
      }
      world.gop.fence();
      alpha += matrix_inner(world, response_k, vp_k);
    }
    world.gop.fence();
    alpha *= alpha_factor;
    world.gop.fence();
    pm.set_alpha(omega, alpha, directions);
  }
  pm.save();
}
