#pragma once
#include "MolecularProperty.hpp"
#include "ResponseIO.hpp"
#include "ResponseVector.hpp"
#include "VBCMacrotask.hpp"
#include "broadcast_json.hpp"
#include "functypedefs.h"
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <madness/external/nlohmann_json/json.hpp>
#include <madness/tensor/tensor.h>
#include <madness/tensor/tensor_json.hpp>
#include <string>

namespace fs = std::filesystem;
using json = nlohmann::json;
#pragma once
#include <chrono>
#include <fstream>
#include <iomanip>
#include <madness/external/nlohmann_json/json.hpp>
#include <madness/tensor/tensor.h>
#include <madness/world/world.h>

using json = nlohmann::json;

inline std::string iso_timestamp() {
  using namespace std::chrono;
  auto now = system_clock::now();
  auto itt = system_clock::to_time_t(now);
  std::ostringstream ss;
  ss << std::put_time(std::gmtime(&itt), "%Y%m%dT%H%M%SZ");
  return ss.str();
}

/// Helper: compute inner-product tensor, and optionally dump per-k
/// contributions.
inline madness::Tensor<double> compute_response_inner_product_tensor(
    madness::World &world, const std::vector<vector_real_function_3d> &A_vecs,
    const std::vector<vector_real_function_3d> &B_vecs,
    bool save_contributions = false, const std::string &base_filename = "") {
  const size_t nA = A_vecs.size();
  const size_t nB = B_vecs.size();
  if (nA == 0 || nB == 0)
    throw std::runtime_error("Input vectors must not be empty.");

  const size_t num_rf = A_vecs[0].size();
  madness::Tensor<double> result(nA, nB);

  json j;
  if (save_contributions) {
    j["num_response_functions"] = num_rf;
    j["contributions"] = json::object();
    j["timestamp"] = iso_timestamp();
  }

  // loop over each response-function index k
  for (size_t k = 0; k < num_rf; ++k) {
    vector_real_function_3d Ak(nA), Bk(nB);
    for (size_t i = 0; i < nA; ++i)
      Ak[i] = A_vecs[i][k];
    for (size_t j = 0; j < nB; ++j)
      Bk[j] = B_vecs[j][k];

    world.gop.fence();

    auto M_k = matrix_inner(world, Ak, Bk);
    result += M_k;

    if (save_contributions) {
      // serialize M_k into JSON
      j["contributions"][std::to_string(k)] =
          madness::tensor_to_json<double>(M_k);
    }
  }

  if (save_contributions) {
    // build timestamped filename
    auto ts = j["timestamp"].get<std::string>();
    std::string fname = base_filename.empty()
                            ? ("inner_contribs_" + ts + ".json")
                            : (base_filename + "_" + ts + ".json");
    std::ofstream out(fname);
    out << std::setw(2) << j << "\n";
    if (world.rank() == 0)
      std::cout << "📂 Wrote per-k contributions to " << fname << "\n";
  }

  return result;
}

/*madness::Tensor<double> compute_response_inner_product_tensor(*/
/*    World &world, const std::vector<vector_real_function_3d> &A_vecs,*/
/*    const std::vector<vector_real_function_3d> &B_vecs) {*/
/*  const size_t nA = A_vecs.size();*/
/*  const size_t nB = B_vecs.size();*/
/**/
/*  if (nA == 0 || nB == 0)*/
/*    throw std::runtime_error("Input vectors must not be empty.");*/
/**/
/*  const size_t num_rf =*/
/*      A_vecs[0].size(); // assuming all vectors have the same response size*/
/*  madness::Tensor<double> result(nA, nB);*/
/**/
/*  for (size_t k = 0; k < num_rf; ++k) {*/
/*    vector_real_function_3d Ak(nA), Bk(nB);*/
/*    for (size_t i = 0; i < nA; ++i)*/
/*      Ak[i] = A_vecs[i][k];*/
/*    for (size_t j = 0; j < nB; ++j)*/
/*      Bk[j] = B_vecs[j][k];*/
/**/
/*    result += matrix_inner(world, Ak, Bk);*/
/*  }*/
/**/
/*  return result;*/
/*}*/

class PropertyManager {
public:
  explicit PropertyManager(World &world, const std::string &filename)
      : filename_(filename) {
    if (fs::exists(filename_)) {
      data_ = broadcast_json_file(world, filename_); // Load JSON data
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
    //
    return data_.contains("polarizability") &&
           data_["polarizability"].contains(freq_str(omega)) &&
           data_["polarizability"][freq_str(omega)].contains(dir) &&
           !data_["polarizability"][freq_str(omega)][dir].is_null();
  }

  // Store a 3-element beta result for a given BC input pair (e.g., xx) at ω1,
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
           !beta[freq_str(omega1)][freq_str(omega2)][bc].is_null() && false;
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

/**
 * @brief Computes the polarizability tensor (α) for a given set of
 * frequencies and directions.
 *
 * @param world The MADNESS world object.
 * @param state_map A map of response states computed from perturbations.
 * @param gs The ground state data.
 * @param frequencies A vector of frequencies for which to compute α.
 * @param directions A string of characters representing the directions
 * for the perturbations.
 * @param pm The property manager to store the computed α tensor.
 */
void compute_alpha(World &world,
                   std::map<std::string, ResponseState> &state_map,
                   const GroundStateData &gs,
                   const std::vector<double> &frequencies,
                   const std::string &directions, PropertyManager &pm) {
  const size_t num_directions = directions.size();
  const size_t num_orbitals = gs.getNumOrbitals();
  const size_t num_frequencies = frequencies.size();
  const bool is_restricted = gs.isSpinRestricted();

  if (world.rank() == 0) {
    print("▶️ Computing α tensor for", num_directions, "directions and",
          num_frequencies, "frequencies.");
  }

  std::vector<vector_real_function_3d> perturbations;
  std::vector<ResponseVector> load_vector(num_directions);

  std::vector<std::string> direction_keys;
  for (const char &dir : directions) {
    direction_keys.push_back(std::string("Dipole_") + std::string(1, dir));
  }
  if (world.rank() == 0) {
    print("Direction keys:", direction_keys);
  }

  for (const auto &dir : direction_keys) {
    if (state_map.find(dir) == state_map.end()) {
      throw std::runtime_error("State not found in state_map: " + dir);
    } else {
      perturbations.push_back(state_map.at(dir).perturbation_vector(world, gs));
    }
  }

  for (size_t f = 0; f < num_frequencies; ++f) {
    double omega = frequencies[f];
    if (pm.has_alpha(omega, directions)) {
      if (world.rank() == 0)
        print("✅ Skipping already computed α at ω =", omega);
      continue;
    }

    double alpha_factor = (omega == 0.0) ? -4.0 : -2.0;
    if (world.rank() == 0) {
      print("🛠️  Computing α at ω =", omega, "for directions:", directions,
            " alpha factor = ", alpha_factor);
    }

    std::vector<vector_real_function_3d> response_vecs(num_directions);
    std::vector<vector_real_function_3d> perturb_vecs(num_directions);

    for (size_t j = 0; j < num_directions; ++j) {
      auto &active_state = state_map.at(direction_keys[j]);
      active_state.set_frequency_index(f);
      load_response_vector(world, num_orbitals, active_state, load_vector[j],
                           active_state.thresholds.size() - 1, f);
      response_vecs[j] = get_flat(load_vector[j]);
      perturb_vecs[j] = perturbations[j];

      if (omega != 0.0) {
        // Duplicate for dynamic response
        perturb_vecs[j].insert(perturb_vecs[j].end(), perturb_vecs[j].begin(),
                               perturb_vecs[j].end());
      }
    }

    // Compute α using the utility
    madness::Tensor<double> alpha = compute_response_inner_product_tensor(
        world, response_vecs, perturb_vecs, true,
        "alpha_contribs_" + std::to_string(f));

    alpha *= alpha_factor;
    if (world.rank() == 0) {
      print("α tensor size:", alpha.size(), "x", alpha.size());
      print("α= ", alpha);
    }

    pm.set_alpha(omega, alpha, directions);
  }
  pm.save();
}

void compute_beta(World &world, const GroundStateData &gs,
                  const std::vector<double> &frequencies,
                  const std::string &directions, PropertyManager &pm) {
  const bool is_spin_restricted = gs.isSpinRestricted();
  const size_t num_orbitals = gs.getNumOrbitals();

  auto vbc_computer =
      VBCComputer(world, gs, frequencies, directions, is_spin_restricted);

  const auto &BC_pairs = vbc_computer.get_BC_pairs();
  const size_t num_pairs = BC_pairs.size();

  for (size_t i = 0; i < frequencies.size(); ++i) {
    for (size_t j = i; j < frequencies.size(); ++j) {
      double omega1 = frequencies[i];
      double omega2 = frequencies[j];
      double omegaA = omega1 + omega2;

      for (size_t bc_index = 0; bc_index < num_pairs; ++bc_index) {
        auto [B, C] = BC_pairs[bc_index];
        std::string bc = std::string() + B + C;

        if (pm.has_beta(omega1, omega2, bc)) {
          if (world.rank() == 0)
            print("✅ Skipping β(", bc, ") at (ω1, ω2) =", omega1, omega2);
          continue;
        }

        if (world.rank() == 0)
          print("🛠️  Computing β(", bc, ") at (ω1, ω2) =", omega1, omega2);

        // Load or compute the VBC vector
        auto vbc_vec = vbc_computer.compute_and_save(bc_index, i, j);

        // Load LHS response states (for output directions)
        std::vector<ResponseState> lhs_states;
        std::vector<ResponseVector> lhs_responses(directions.size());

        for (size_t k = 0; k < directions.size(); ++k) {
          char A = directions[k];
          auto state = vbc_computer.get_state(A);
          state.set_frequency_index(state.frequency_map.at(omegaA));
          lhs_states.push_back(state);

          load_response_vector(world, num_orbitals, state, lhs_responses[k], 0,
                               state.current_frequency_index);
        }

        // Convert LHS responses to flat components
        std::vector<vector_real_function_3d> lhs_flat;

        for (size_t k = 0; k < directions.size(); ++k) {
          auto &resp = lhs_responses[k];
          auto flat = get_flat(resp);
          if (omegaA == 0.0) {
            flat.insert(flat.end(), flat.begin(),
                        flat.end()); // replicate for static
          }
          lhs_flat.push_back(flat);
        }

        // Convert VBC to flat
        auto vbc_flat = get_flat(vbc_vec);
        // Compute inner product (lhs ⋅ vbc) → Tensor(3)
        madness::Tensor<double> beta_tensor =
            compute_response_inner_product_tensor(
                world, lhs_flat, {vbc_flat}, true,
                "beta_contribs" + std::to_string(i));
        beta_tensor *= -2.0; // factor of -2 for beta
        if (world.rank() == 0) {
          print("β tensor size:", beta_tensor.size(), "x", beta_tensor.size());
          print("β= ", beta_tensor);
        }

        pm.set_beta(omega1, omega2, bc, beta_tensor);
      }
    }
  }

  pm.save();
}
