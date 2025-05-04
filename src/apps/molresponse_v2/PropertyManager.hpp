#pragma once
#include <madness/tensor/tensor.h>

#include <filesystem>
#include <fstream>
#include <iomanip>
#include <madness/external/nlohmann_json/json.hpp>
#include <madness/tensor/tensor_json.hpp>
#include <string>

#include "MolecularProperty.hpp"
#include "ResponseIO.hpp"
#include "ResponseParameters.hpp"
#include "ResponseVector.hpp"
#include "VBCMacrotask.hpp"
#include "broadcast_json.hpp"
#include "functypedefs.h"

namespace fs = std::filesystem;
using json = nlohmann::json;
#pragma once
#include <madness/tensor/tensor.h>
#include <madness/world/world.h>

#include <chrono>
#include <fstream>
#include <iomanip>
#include <madness/external/nlohmann_json/json.hpp>

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
inline madness::Tensor<double> compute_response_inner_product_tensor(madness::World &world, const std::vector<vector_real_function_3d> &A_vecs,
                                                                     const std::vector<vector_real_function_3d> &B_vecs, bool save_contributions = false,
                                                                     const std::string &base_filename = "") {
  const size_t nA = A_vecs.size();
  const size_t nB = B_vecs.size();
  if (nA == 0 || nB == 0) throw std::runtime_error("Input vectors must not be empty.");

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
    for (size_t i = 0; i < nA; ++i) Ak[i] = A_vecs[i][k];
    for (size_t j = 0; j < nB; ++j) Bk[j] = B_vecs[j][k];

    world.gop.fence();

    auto M_k = matrix_inner(world, Ak, Bk);
    result += M_k;

    if (save_contributions) {
      // serialize M_k into JSON
      j["contributions"][std::to_string(k)] = madness::tensor_to_json<double>(M_k);
    }
  }

  if (save_contributions) {
    // build timestamped filename
    auto ts = j["timestamp"].get<std::string>();
    std::string fname = base_filename.empty() ? ("inner_contribs_" + ts + ".json") : (base_filename + "_" + ts + ".json");
    std::ofstream out(fname);
    out << std::setw(2) << j << "\n";
    if (world.rank() == 0) std::cout << "📂 Wrote per-k contributions to " << fname << "\n";
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
  explicit PropertyManager(World &world, const std::string &filename) : filename_(filename) {
    if (fs::exists(filename_)) {
      data_ = broadcast_json_file(world, filename_);  // Load JSON data
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
  void set_alpha(double omega, const madness::Tensor<double> &tensor, const std::string &dirs) {
    data_["polarizability"][freq_str(omega)] = tensor_to_json(tensor);
    data_["polarizability"][freq_str(omega)]["directions"] = dirs;
  }

  bool has_alpha(double omega, std::string dir) const {
    // also check if they are the same size
    //
    return data_.contains("polarizability") && data_["polarizability"].contains(freq_str(omega)) && data_["polarizability"][freq_str(omega)].contains(dir) &&
           !data_["polarizability"][freq_str(omega)][dir].is_null();
  }

  // Store a 3-element beta result for a given BC input pair (e.g., xx) at ω1,
  void set_beta(double omega1, double omega2, const std::string &bc, const madness::Tensor<double> &tensor) {
    data_["hyperpolarizability"][freq_str(omega1)][freq_str(omega2)][bc] = tensor_to_json(tensor);
  }
  void set_beta_dirs(const std::string &dirs) { data_["hyperpolarizability"]["directions"] = dirs; }

  bool has_beta(double omega1, double omega2, const std::string &bc) const {
    auto &beta = data_["hyperpolarizability"];
    return beta.contains(freq_str(omega1)) && beta[freq_str(omega1)].contains(freq_str(omega2)) && beta[freq_str(omega1)][freq_str(omega2)].contains(bc) &&
           !beta[freq_str(omega1)][freq_str(omega2)][bc].is_null() && false;
  }

  void print_alpha_table() const {
    if (!data_.contains("polarizability")) {
      std::cout << "No polarizability data found.\n";
      return;
    }
    std::cout << "\n📐 Polarizability (α) Tensor Components:\n";
    // Build header based on available directions
    for (const auto &[freq_str, tensor_json] : data_["polarizability"].items()) {
      if (!tensor_json.contains("directions")) continue;

      const std::string dirs = tensor_json["directions"];
      const madness::Tensor<double> alpha = tensor_from_json<double>(tensor_json);

      // Print header
      std::cout << "\nω = " << freq_str << "\n";
      std::cout << std::setw(10) << "";
      for (char d : dirs) std::cout << std::setw(12) << "α_" + std::string(1, d);
      std::cout << "\n";

      // Print values row by row
      bool has_values = alpha.has_data();

      for (size_t i = 0; i < dirs.size(); ++i) {
        std::cout << std::setw(10) << dirs[i];
        for (size_t j = 0; j < dirs.size(); ++j) {
          if (i < alpha.dim(0) && j < alpha.dim(1) && has_values) {
            std::cout << std::setw(12) << std::fixed << std::setprecision(6) << alpha(i, j);
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
    if (data_["hyperpolarizability"].contains("directions")) dirs = data_["hyperpolarizability"]["directions"].get<std::string>();

    for (const auto &[w1_str, w1_entry] : data_["hyperpolarizability"].items()) {
      if (w1_str == "directions") continue;  // skip metadata

      for (const auto &[w2_str, w2_entry] : w1_entry.items()) {
        std::cout << "\nω₁ = " << w1_str << ", ω₂ = " << w2_str << "\n";
        std::cout << std::setw(8) << "BC";
        for (char A : dirs) std::cout << std::setw(12) << "β_" + std::string(1, A) + "BC";
        std::cout << "\n";

        for (const auto &[bc, tensor_json] : w2_entry.items()) {
          madness::Tensor<double> beta = tensor_from_json<double>(tensor_json);
          std::cout << std::setw(8) << bc;

          for (int i = 0; i < dirs.size(); ++i) {
            if (i < beta.size())
              std::cout << std::setw(12) << std::fixed << std::setprecision(6) << beta(i);
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

void initialize_property_structure(PropertyManager &pm, const ResponseParameters &rp) {
  auto props = rp.requested_properties();
  auto dipole_dirs = rp.dipole_directions();
  auto nuclear_dirs = rp.nuclear_directions();
  auto nuclear_atom_indices = rp.nuclear_atom_indices();
  auto dipole_freqs = rp.dipole_frequencies();
  auto nuclear_freqs = rp.nuclear_frequencies();



  for (const auto &prop : props) {
    if (prop == "polarizability") {
      auto num_dirs = dipole_dirs.size();

      std::string directions_string;
      for (char dir : dipole_dirs) {
        directions_string += dir;
      }

      for (double omega : dipole_freqs) {
        if (!pm.has_alpha(omega, directions_string)) {
          madness::Tensor<double> empty_tensor(num_dirs, num_dirs,
                                               0.0);  // 3x3 alpha
          pm.set_alpha(omega, empty_tensor, directions_string);
        }
      }
    } else if (prop == "hyperpolarizability") {
      const auto directions_string = std::string(dipole_dirs.begin(), dipole_dirs.end());
      pm.set_beta_dirs(directions_string);

      for (size_t b = 0; b < dipole_freqs.size(); ++b) {
        for (size_t c = b; c < dipole_freqs.size(); ++c) {
          double omega1 = dipole_freqs[b];
          double omega2 = dipole_freqs[c];
          for (char B : dipole_dirs) {
            for (char C : dipole_dirs) {
              std::string bc = std::string() + B + C;
              if (!pm.has_beta(omega1, omega2, bc)) {
                madness::Tensor<double> empty(3, 0.0);  // x(xx), y(xx), z(xx)
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
void compute_alpha(World &world, std::map<std::string, ResponseState> &state_map, const GroundStateData &gs, const std::vector<double> &frequencies,
                   const std::string &directions, PropertyManager &pm) {
  const size_t num_directions = directions.size();
  const size_t num_orbitals = gs.getNumOrbitals();
  const size_t num_frequencies = frequencies.size();
  const bool is_restricted = gs.isSpinRestricted();

  if (world.rank() == 0) {
    print("▶️ Computing α tensor for", num_directions, "directions and", num_frequencies, "frequencies.");
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
      perturbations.push_back(perturbation_vector(world, gs, state_map.at(dir)));  // Get the perturbation vector
    }
  }

  for (size_t f = 0; f < num_frequencies; ++f) {
    double omega = frequencies[f];
    if (pm.has_alpha(omega, directions)) {
      if (world.rank() == 0) print("✅ Skipping already computed α at ω =", omega);
      continue;
    }

    double alpha_factor = (omega == 0.0) ? -4.0 : -2.0;
    if (world.rank() == 0) {
      print("🛠️  Computing α at ω =", omega, "for directions:", directions, " alpha factor = ", alpha_factor);
    }

    std::vector<vector_real_function_3d> response_vecs(num_directions);
    std::vector<vector_real_function_3d> perturb_vecs(num_directions);

    for (size_t j = 0; j < num_directions; ++j) {
      auto &active_state = state_map.at(direction_keys[j]);
      active_state.set_frequency_index(f);
      load_response_vector(world, num_orbitals, active_state, load_vector[j], active_state.thresholds.size() - 1, f);
      response_vecs[j] = get_flat(load_vector[j]);
      perturb_vecs[j] = perturbations[j];

      if (omega != 0.0) {
        // Duplicate for dynamic response
        perturb_vecs[j].insert(perturb_vecs[j].end(), perturb_vecs[j].begin(), perturb_vecs[j].end());
      }
    }

    // Compute α using the utility
    madness::Tensor<double> alpha = compute_response_inner_product_tensor(world, response_vecs, perturb_vecs, true, "alpha_contribs_" + std::to_string(f));

    alpha *= alpha_factor;
    if (world.rank() == 0) {
      print("α tensor size:", alpha.size(), "x", alpha.size());
      print("α= ", alpha);
    }

    pm.set_alpha(omega, alpha, directions);
  }
  pm.save();
}

void compute_beta(World &world, const GroundStateData &gs, const PerturbationType a_type, const std::vector<Perturbation> &perturbation_A,
                  const std::pair<PerturbationType, PerturbationType> &bc_types, const std::vector<std::pair<Perturbation, Perturbation>> &BC_pairs,
                  const std::pair<std::vector<double>, std::vector<double>> &frequencies, PropertyManager &pm) {
  const bool is_spin_restricted = gs.isSpinRestricted();
  const int num_orbitals = static_cast<int>(gs.getNumOrbitals());
  const double thresh = FunctionDefaults<3>::get_thresh();

  // 1) Build a SimpleVBCComputer once
  auto vbc_computer = SimpleVBCComputer(world, gs);

  // 2) Loop over all freq
  for (auto &freq_b : frequencies.first) {
    for (auto &freq_c : frequencies.second) {
      // Get all BC pairs possible out of Perturbation
      for (auto [B, C] : BC_pairs) {
        VBCResponseState vbc_state(bc_types.first, bc_types.second, B, C, freq_b, freq_c, FunctionDefaults<3>::get_thresh(), is_spin_restricted);
        auto bc = vbc_state.perturbationDescription();

        if (pm.has_beta(freq_b, freq_c, bc)) {
          if (world.rank() == 0) print("✅ Skipping β(", bc, ") at (ω1, ω2) =", freq_b, freq_c);
          continue;
        }

        auto vbc_vec = vbc_computer.compute_and_save(vbc_state);
        auto vbc_flat = get_flat(vbc_vec);

        auto omega_A = vbc_state.current_frequency();
        std::vector<ResponseVector> a_vecs(perturbation_A.size());
        std::vector<vector_real_function_3d> xa_vecs(perturbation_A.size());  // for the A perturbations
        //

        std::vector<real_function_3d> opAs(perturbation_A.size());
        for (int a = 0; a < perturbation_A.size(); ++a) {
          auto pertA = perturbation_A[a];
          auto state_A = ResponseState(pertA, a_type, {omega_A}, {FunctionDefaults<3>::get_thresh()}, is_spin_restricted);
          opAs[a] = raw_perturbation_operator(world, gs, state_A);
          load_response_vector(world, num_orbitals, state_A, a_vecs[a], 0, 0);
          auto flat = get_flat(a_vecs[a]);
          if (omega_A == 0.0) {
            flat.insert(flat.end(), flat.begin(),
                        flat.end());  // replicate for static
          }
          xa_vecs[a] = -1.0 * flat;
        }

        madness::Tensor<double> beta_tensor = compute_response_inner_product_tensor(world, xa_vecs, {vbc_flat}, true, "beta_contribs" + bc);

        // We only need one copy of xb_phi0
        DynamicRestrictedResponse xb_phi0(num_orbitals);
        DynamicRestrictedResponse xc_phi0(num_orbitals);
        DynamicRestrictedResponse y_zeta_bc(num_orbitals);
        DynamicRestrictedResponse y_zeta_cb(num_orbitals);

        auto [xb, xc] = vbc_computer.get_BC_vecs(vbc_state);  // get the B and C states

        xb_phi0.x_alpha = std::get<DynamicRestrictedResponse>(xb).x_alpha;
        xb_phi0.y_alpha = gs.orbitals;
        xc_phi0.x_alpha = std::get<DynamicRestrictedResponse>(xc).x_alpha;
        xc_phi0.y_alpha = gs.orbitals;
        y_zeta_bc.x_alpha = std::get<DynamicRestrictedResponse>(xc).y_alpha;
        y_zeta_cb.x_alpha = std::get<DynamicRestrictedResponse>(xb).y_alpha;
        y_zeta_bc.y_alpha = SimpleVBCComputer::make_zeta_bc(world, std::get<DynamicRestrictedResponse>(xb).y_alpha,
                                                            std::get<DynamicRestrictedResponse>(xc).x_alpha, gs.orbitals);
        y_zeta_cb.y_alpha = SimpleVBCComputer::make_zeta_bc(world, std::get<DynamicRestrictedResponse>(xc).y_alpha,
                                                            std::get<DynamicRestrictedResponse>(xb).x_alpha, gs.orbitals);

        xb_phi0.flatten();
        xc_phi0.flatten();
        y_zeta_bc.flatten();
        y_zeta_cb.flatten();

        std::vector<vector_real_function_3d> ra_y_zeta_bc(perturbation_A.size());
        std::vector<vector_real_function_3d> ra_y_zeta_cb(perturbation_A.size());
        for (int a = 0; a < perturbation_A.size(); ++a) {
          ra_y_zeta_bc[a] = copy(world, get_flat(y_zeta_bc)) * opAs[a];
          ra_y_zeta_cb[a] = copy(world, get_flat(y_zeta_cb)) * opAs[a];
        }
        std::vector<vector_real_function_3d> xb_phi0_vec(1);
        std::vector<vector_real_function_3d> xc_phi0_vec(1);
        xb_phi0_vec[0] = get_flat(xb_phi0);
        xc_phi0_vec[0] = get_flat(xc_phi0);

        auto file_name = "responses/" + vbc_state.perturbationDescription() + "_beta_zeta_bc_0.json";
        madness::Tensor<double> beta_2 = compute_response_inner_product_tensor(world, ra_y_zeta_bc, xb_phi0_vec, true, "beta_zeta_bc_1");

        file_name = "responses/" + vbc_state.perturbationDescription() + "_beta_zeta_bc_1.json";
        madness::Tensor<double> beta_3 = compute_response_inner_product_tensor(world, ra_y_zeta_cb, xc_phi0_vec, true, file_name);
        if (world.rank() == 0) {
          print("beta_3 tensor size:", beta_3.size(), "x", beta_3.size());
          print("beta_3= \n ", beta_3);
        }

        beta_tensor += beta_2;
        beta_tensor += beta_3;

        beta_tensor *= -2.0;  // factor of -2 for beta
        if (world.rank() == 0) {
          print(vbc_state.perturbationDescription());
          print("β tensor size:", beta_tensor.size(), "x", beta_tensor.size());
          print("β= \n", beta_tensor);
        }

        pm.set_beta(freq_b, freq_c, bc, beta_tensor);
      }
    }
  }
  pm.save();
}

// ─────────────────────────────────────────────────────────────────────────────
// Compute the (frequency-dependent) hyperpolarizability β(ωA;ωB,ωC)
//  ωB and ωC each run over the same list of input frequencies.
//  directions = subset of {'x','y','z'} to include.
// ─────────────────────────────────────────────────────────────────────────────
void compute_hyperpolarizability(World &world, const GroundStateData &gs,
                                 const std::vector<double> &frequencies,  // ωB and ωC
                                 const std::vector<char> &directions,     // which Cartesian dirs
                                 PropertyManager &pm) {
  // 1) Build A-list (dipole along each dir)
  std::vector<Perturbation> A_pert;
  A_pert.reserve(directions.size());
  for (char d : directions) A_pert.emplace_back(DipolePerturbation{d});

  // 2) Build all BC pairs of two dipoles
  std::pair<PerturbationType, PerturbationType> bc_types = {PerturbationType::Dipole, PerturbationType::Dipole};
  std::vector<std::pair<Perturbation, Perturbation>> BC_pairs;
  for (char b : directions)
    for (char c : directions) BC_pairs.emplace_back(DipolePerturbation{b}, DipolePerturbation{c});

  auto freq_pair = std::make_pair(frequencies, frequencies);

  // 3) Dispatch
  compute_beta(world, gs,
               PerturbationType::Dipole,  // A-type
               A_pert,
               bc_types,  // B and C types
               BC_pairs,
               freq_pair,  // ωB and ωC = same list of freqs
               pm);
}

// ─────────────────────────────────────────────────────────────────────────────
// Compute the Raman response (∂α/∂Q) via β(Dipole; Dipole, NuclearDisp).
//  frequencies: the optical frequencies for both Dipole and Raman.
//  dip_dirs: which dipole directions to include.
// ─────────────────────────────────────────────────────────────────────────────
void compute_Raman(World &world, const GroundStateData &gs, const std::pair<std::vector<double>, std::vector<double>> &BC_frequencies,
                   const std::vector<char> &dip_dirs,  // e.g. {'x','y','z'}
                   PropertyManager &pm) {
  // 1) A-list is still the dipole directions
  std::vector<Perturbation> A_pert;
  A_pert.reserve(dip_dirs.size());
  for (char d : dip_dirs) A_pert.emplace_back(DipolePerturbation{d});

  // 2) Enumerate all nuclear displacements (per atom, per axis)
  std::vector<Perturbation> nucs;
  nucs.reserve(3 * gs.molecule.natom());
  for (int atom = 0; atom < gs.molecule.natom(); ++atom) {
    for (char dir : {'x', 'y', 'z'}) nucs.emplace_back(NuclearDisplacementPerturbation{atom, dir});
  }

  std::pair<PerturbationType, PerturbationType> bc_types = {PerturbationType::Dipole, PerturbationType::NuclearDisplacement};

  // 3) Build BC pairs = (Dipole, NuclearDisp)
  std::vector<std::pair<Perturbation, Perturbation>> BC_pairs;
  BC_pairs.reserve(dip_dirs.size() * nucs.size());
  for (char b : dip_dirs) {
    for (auto &n : nucs) {
      BC_pairs.emplace_back(DipolePerturbation{b}, std::get<NuclearDisplacementPerturbation>(n));
    }
  }

  // 4) Dispatch
  compute_beta(world, gs,
               PerturbationType::Dipole,  // the “A”-perturbation is always
               A_pert, bc_types, BC_pairs,
               BC_frequencies,  // ωB and ωC = same list of freqs
               pm);
}
