#pragma once
#include <madness/tensor/tensor.h>
#include <madness/world/world.h>

#include <chrono>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <limits>
#include <madness/external/nlohmann_json/json.hpp>
#include <madness/tensor/tensor_json.hpp>
#include <string>

#include "InnerContributions.hpp"
#include "MolecularProperty.hpp"
#include "ResponseIO.hpp"
#include "../../madness/chem/ResponseParameters.hpp"
#include "ResponseSolverUtils.hpp"
#include "ResponseVector.hpp"
#include "VBCMacrotask.hpp"
#include "broadcast_json.hpp"
#include "functypedefs.h"

using json = nlohmann::json;

inline std::string iso_timestamp() {
  using namespace std::chrono;
  auto now = system_clock::now();
  auto itt = system_clock::to_time_t(now);
  std::ostringstream ss;
  ss << std::put_time(std::gmtime(&itt), "%Y%m%dT%H%M%SZ");
  return ss.str();
}

inline madness::Tensor<double> compute_response_inner_product_tensor(
    madness::World &world, const std::vector<vector_real_function_3d> &A_vecs,
    const std::vector<vector_real_function_3d> &B_vecs,
    bool save_contributions = false, const std::string &entry_name = "") {
  const size_t nA = A_vecs.size();
  const size_t nB = B_vecs.size();
  if (nA == 0 || nB == 0)
    throw std::runtime_error("Input vectors must not be empty.");

  const size_t num_rf = A_vecs[0].size();
  madness::Tensor<double> result(nA, nB);

  // if they asked us to save
  json this_entry;
  if (save_contributions) {
    this_entry["num_response_functions"] = num_rf;
    this_entry["contributions"] = json::object();
    this_entry["timestamp"] = iso_timestamp();
  }

  // accumulate
  for (size_t k = 0; k < num_rf; ++k) {
    vector_real_function_3d Ak(nA), Bk(nB);
    for (size_t i = 0; i < nA; ++i) Ak[i] = A_vecs[i][k];
    for (size_t j = 0; j < nB; ++j) Bk[j] = B_vecs[j][k];

    world.gop.fence();

    auto M_k = matrix_inner(world, Ak, Bk);
    result += M_k;

    if (save_contributions) {
      // serialize M_k
      this_entry["contributions"][std::to_string(k)] =
          madness::tensor_to_json<double>(M_k);
    }
  }

  if (save_contributions) {
    // choose a key to store under; if user gave one, use it, else timestamp
    auto &g = global_inner_contributions();

    std::string key = entry_name.empty()
                          ? this_entry["timestamp"].get<std::string>()
                          : entry_name;
    g[key] = std::move(this_entry);
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

struct PropRow {
  std::string property;
  std::string
      component;  // x,y,z, xx,xy,xz,yy,yz,zz, xxx, xxy, xxz, yyy, yyz, zzz
  double freq1;
  std::optional<double> freq2;
  std::optional<double> value;
};

struct PropKey {
  std::string property;
  std::string component;
  double freq1;
  std::optional<double> freq2;

  bool operator<(PropKey const &other) const noexcept {
    if (property != other.property) return property < other.property;
    if (component != other.component) return component < other.component;
    if (freq1 != other.freq1) return freq1 < other.freq1;
    // for freq2
    if (!freq2 && other.freq2) return true;
    if (freq2 && !other.freq2) return false;
    if (freq2 && other.freq2) return *freq2 < *other.freq2;
    return false;
  }
};
// enable Json to PropRow
inline void to_json(json &j, PropRow const &r) {
  j = json::object();
  j["property"] = r.property;
  j["component"] = r.component;
  j["freqB"] = r.freq1;
  if (r.freq2) {
    j["freqC"] = *r.freq2;
  }
  if (r.value) {
    j["value"] = *r.value;
  }
}

inline void from_json(json const &j, PropRow &r) {
  r.property = j.at("property").get<std::string>();
  r.component = j.at("component").get<std::string>();
  r.freq1 = j.at("freqB").get<double>();
  if (j.contains("freqC")) {
    r.freq2 = j.at("freqC").get<double>();
  }
  if (j.contains("value")) {
    r.value = j.at("value").get<double>();
  }
}

class PropertyManager {
 public:
  explicit PropertyManager(World &world, const std::string &filename)
      : filename_(filename) {
    if (fs::exists(filename_)) {
      // load JSON
      json j = broadcast_json_file(world, filename_);
      if (j.is_array()) {
        for (auto const &r : j) {
          PropRow row = r.get<PropRow>();
          PropKey key = {row.property, row.component, row.freq1, row.freq2};
          rows_[key] = std::move(row);
          // flat‐format: just parse rows
        }
      } else {
        // old nested format: convert to rows_
        parse_old_format(j);
      }
    }
  }

  /// Return the flat rows as a JSON array
  json to_json() const {
    json a = json::array();
    for (const auto &[key, row] : rows_) {
      a.push_back(row);
    }
    return a;
  }

  /// Overwrite file with flat JSON array
  void save() const {
    std::ofstream out(filename_);
    out << std::setw(12) << to_json() << "\n";
  }

  // Presence checks
  [[nodiscard]] bool has_alpha(double omega, std::string comp) const {
    PropKey k{"polarizability", comp, omega, std::nullopt};
    auto it = rows_.find(k);
    return it != rows_.end() && it->second.value.has_value();
    return rows_.count(k) != 0;
  }
  [[nodiscard]] bool has_beta(double w1, double w2, std::string comp) const {
    PropKey k{"hyperpolarizability", comp, w1, w2};
    auto it = rows_.find(k);
    return it != rows_.end() && it->second.value.has_value();
  }

  // Insert or overwrite α entries
  void set_alpha(double omega, const madness::Tensor<double> &tensor,
                 const std::string &dirs) {
    size_t N = dirs.size();
    for (size_t i = 0; i < N; ++i) {
      for (size_t j = 0; j < N; ++j) {
        std::string comp = std::string{dirs[i]} + dirs[j];
        PropKey k{"polarizability", comp, omega, std::nullopt};
        PropRow r{
            .property = k.property,
            .component = k.component,
            .freq1 = k.freq1,
            .freq2 = std::nullopt,
            .value = tensor(i, j),
        };
        rows_[k] = std::move(r);
      }
    }
  }

  // Insert or overwrite β entries
  void set_beta(double w1, double w2, const std::string &comp, double value) {
    PropKey k{"hyperpolarizability", comp, w1, w2};
    PropRow r{
        .property = k.property,
        .component = k.component,
        .freq1 = k.freq1,
        .freq2 = k.freq2,
        .value = value,
    };
    rows_[k] = std::move(r);
  }

  // 3) Append new α rows

  /// Print all PropRow entries in a fixed‐width table
  void print_table() const {
    if (rows_.empty()) {
      std::cout << "No property rows to display.\n";
      return;
    }
    // 1) Header
    std::cout << "\n📊 Property Results\n";
    std::cout << std::left << std::setw(18) << "Property" << std::setw(8)
              << "Comp" << std::setw(8) << "ω1" << std::setw(8) << "ω2"
              << std::setw(12) << "Value" << "\n";

    // 2) Divider line
    std::cout << std::string(18 + 8 + 8 + 8 + 12, '_') << "\n";

    // 3) Rows
    for (auto const &[k, r] : rows_) {
      // format ω1 and ω2
      std::ostringstream o1, o2;
      o1 << std::fixed << std::setprecision(3) << r.freq1;
      if (r.freq2) {
        o2 << std::fixed << std::setprecision(3) << *r.freq2;
      } else {
        o2 << "-";
      }

      std::cout << std::left << std::setw(22) << r.property << std::setw(8)
                << r.component << std::setw(8) << o1.str() << std::setw(8)
                << o2.str() << std::setw(12) << std::fixed
                << std::setprecision(6) << *r.value << "\n";
    }
  }

 private:
  std::string filename_;
  std::map<PropKey, PropRow> rows_;
  static std::string freq_str(double freq) {
    std::ostringstream ss;
    ss << std::scientific << std::setprecision(3) << freq;
    return ss.str();
  }

  // Helper: convert old nested JSON into flat rows_
  void parse_old_format(json const &old) {
    // --- polarizability ---
    if (old.contains("polarizability")) {
      for (auto const &[freq_s, obj] : old["polarizability"].items()) {
        double omega = std::stod(freq_s);
        std::string dirs = obj.value("directions", std::string{});
        auto tens = tensor_from_json<double>(obj);
        set_alpha(omega, tens, dirs);
      }
    }
    // --- hyperpolarizability ---
    if (old.contains("hyperpolarizability")) {
      for (auto const &[w1_s, sub] : old["hyperpolarizability"].items()) {
        if (w1_s == "directions") continue;
        double w1 = std::stod(w1_s);
        for (auto const &[w2_s, entry] : sub.items()) {
          double w2 = std::stod(w2_s);
          for (auto const &[bc, tens_json] : entry.items()) {
            auto tens = tensor_from_json<double>(tens_json);
            // flatten each component A in the resulting vector/tensor
            // assume tens is 1×N or N×1:
            for (int k = 0; k < tens.size(); ++k) {
              std::string comp =
                  std::string{tens_json.value("directions", "")[k]} + bc;
              double val = (tens.dim(0) == 1 ? tens(0, k) : tens(k, 0));
              set_beta(w1, w2, comp, val);
            }
          }
        }
      }
    }
  }
};

/*void initialize_property_structure(PropertyManager &pm,*/
/*                                   const ResponseParameters &rp) {*/
/*  auto props = rp.requested_properties();*/
/*  auto dipole_dirs = rp.dipole_directions();*/
/*  auto nuclear_dirs = rp.nuclear_directions();*/
/*  auto nuclear_atom_indices = rp.nuclear_atom_indices();*/
/*  auto dipole_freqs = rp.dipole_frequencies();*/
/*  auto nuclear_freqs = rp.nuclear_frequencies();*/
/**/
/*  for (const auto &prop : props) {*/
/*    if (prop == "polarizability") {*/
/*      auto num_dirs = dipole_dirs.size();*/
/**/
/*      std::string directions_string;*/
/*      for (char dir : dipole_dirs) {*/
/*        directions_string += dir;*/
/*      }*/
/**/
/*      for (double omega : dipole_freqs) {*/
/*        if (!pm.has_alpha(omega, directions_string)) {*/
/*          madness::Tensor<double> empty_tensor(num_dirs,*/
/*                                               num_dirs);  // 3x3 alpha*/
/*          pm.set_alpha(omega, empty_tensor, directions_string);*/
/*        }*/
/*      }*/
/*    } else if (prop == "hyperpolarizability") {*/
/*      const auto directions_string =*/
/*          std::string(dipole_dirs.begin(), dipole_dirs.end());*/
/**/
/*      for (size_t b = 0; b < dipole_freqs.size(); ++b) {*/
/*        for (size_t c = b; c < dipole_freqs.size(); ++c) {*/
/*          double omega1 = dipole_freqs[b];*/
/*          double omega2 = dipole_freqs[c];*/
/*          for (char B : dipole_dirs) {*/
/*            for (char C : dipole_dirs) {*/
/*              std::string bc = std::string() + B + C;*/
/*              for (char A : dipole_dirs) {*/
/*                std::string abc = std::string() + A + B + C;*/
/*                if (!pm.has_beta(omega1, omega2, abc)) {*/
/*                  pm.set_beta(omega1, omega2, abc,*/
/*                              std::numeric_limits<double>::quiet_NaN());*/
/*                }*/
/*              }*/
/*            }*/
/*          }*/
/*        }*/
/*      }*/
/*    }*/
/*  }*/
/**/
/*  pm.save();*/
/*}*/

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
                   std::map<std::string, LinearResponseDescriptor> &state_map,
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
      perturbations.push_back(perturbation_vector(
          world, gs, state_map.at(dir)));  // Get the perturbation vector
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

void compute_beta(
    World &world, const GroundStateData &gs, const PerturbationType a_type,
    const std::vector<Perturbation> &perturbation_A,
    const std::pair<PerturbationType, PerturbationType> &bc_types,
    const std::vector<std::pair<Perturbation, Perturbation>> &BC_pairs,
    const std::pair<std::vector<double>, std::vector<double>> &frequencies,
    PropertyManager &pm) {
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
        VBCResponseState vbc_state(
            bc_types.first, bc_types.second, B, C, freq_b, freq_c,
            FunctionDefaults<3>::get_thresh(), is_spin_restricted);
        auto bc = vbc_state.perturbationDescription();
        auto pertB = std::get<DipolePerturbation>(B);  // get the B perturbation
        auto pertC = std::get<DipolePerturbation>(C);  // get the C perturbation

        std::vector<std::string> abc_string;
        for (auto a : perturbation_A) {
          auto pertA = std::get<DipolePerturbation>(a);
          auto dir = pertA.direction;
          auto comp = std::string(1, dir) + pertB.direction +
                      pertC.direction;  // get the perturbation direction
          abc_string.push_back(comp);
        }

        bool skip = true;
        // skip if all abc are already computed
        for (auto &abc : abc_string) {
          if (!pm.has_beta(freq_b, freq_c, abc)) {
            skip = false;
            break;
          }
        }
        if (skip) {
          if (world.rank() == 0)
            print("✅ Skipping already computed β at ωB =", freq_b,
                  "ωC =", freq_c, "for directions:", abc_string);
          continue;
        }

        auto vbc_vec = vbc_computer.compute_and_save(vbc_state);
        auto vbc_flat = get_flat(vbc_vec);

        auto omega_A = vbc_state.current_frequency();
        std::vector<ResponseVector> a_vecs(perturbation_A.size());
        std::vector<vector_real_function_3d> xa_vecs(
            perturbation_A.size());  // for the A perturbations
        //

        std::vector<real_function_3d> opAs(perturbation_A.size());
        for (int a = 0; a < perturbation_A.size(); ++a) {
          auto pertA = perturbation_A[a];
          auto state_A = LinearResponseDescriptor(
              pertA, a_type, {omega_A}, {FunctionDefaults<3>::get_thresh()},
              is_spin_restricted);
          opAs[a] = raw_perturbation_operator(world, gs, state_A);
          load_response_vector(world, num_orbitals, state_A, a_vecs[a], 0, 0);
          auto flat = get_flat(a_vecs[a]);
          if (omega_A == 0.0) {
            flat.insert(flat.end(), flat.begin(),
                        flat.end());  // replicate for static
          }
          xa_vecs[a] = -1.0 * flat;
        }

        madness::Tensor<double> beta_tensor =
            compute_response_inner_product_tensor(world, xa_vecs, {vbc_flat},
                                                  true, "beta_contribs" + bc);

        // We only need one copy of xb_phi0
        DynamicRestrictedResponse xb_phi0(num_orbitals);
        DynamicRestrictedResponse xc_phi0(num_orbitals);
        DynamicRestrictedResponse y_zeta_bc(num_orbitals);
        DynamicRestrictedResponse y_zeta_cb(num_orbitals);

        auto [xb, xc] =
            vbc_computer.get_BC_vecs(vbc_state);  // get the B and C states

        xb_phi0.x_alpha = std::get<DynamicRestrictedResponse>(xb).x_alpha;
        xb_phi0.y_alpha = gs.orbitals;
        xc_phi0.x_alpha = std::get<DynamicRestrictedResponse>(xc).x_alpha;
        xc_phi0.y_alpha = gs.orbitals;
        y_zeta_bc.x_alpha = std::get<DynamicRestrictedResponse>(xc).y_alpha;
        y_zeta_cb.x_alpha = std::get<DynamicRestrictedResponse>(xb).y_alpha;
        y_zeta_bc.y_alpha = SimpleVBCComputer::make_zeta_bc(
            world, std::get<DynamicRestrictedResponse>(xb).y_alpha,
            std::get<DynamicRestrictedResponse>(xc).x_alpha, gs.orbitals);
        y_zeta_cb.y_alpha = SimpleVBCComputer::make_zeta_bc(
            world, std::get<DynamicRestrictedResponse>(xc).y_alpha,
            std::get<DynamicRestrictedResponse>(xb).x_alpha, gs.orbitals);

        xb_phi0.flatten();
        xc_phi0.flatten();
        y_zeta_bc.flatten();
        y_zeta_cb.flatten();

        std::vector<vector_real_function_3d> ra_y_zeta_bc(
            perturbation_A.size());
        std::vector<vector_real_function_3d> ra_y_zeta_cb(
            perturbation_A.size());
        for (int a = 0; a < perturbation_A.size(); ++a) {
          ra_y_zeta_bc[a] = copy(world, get_flat(y_zeta_bc)) * opAs[a];
          ra_y_zeta_cb[a] = copy(world, get_flat(y_zeta_cb)) * opAs[a];
        }
        std::vector<vector_real_function_3d> xb_phi0_vec(1);
        std::vector<vector_real_function_3d> xc_phi0_vec(1);
        xb_phi0_vec[0] = get_flat(xb_phi0);
        xc_phi0_vec[0] = get_flat(xc_phi0);

        auto file_name = "responses/" + vbc_state.perturbationDescription() +
                         "_beta_zeta_bc_0.json";
        madness::Tensor<double> beta_2 = compute_response_inner_product_tensor(
            world, ra_y_zeta_bc, xb_phi0_vec, true, "beta_zeta_bc_1");

        file_name = "responses/" + vbc_state.perturbationDescription() +
                    "_beta_zeta_bc_1.json";
        madness::Tensor<double> beta_3 = compute_response_inner_product_tensor(
            world, ra_y_zeta_cb, xc_phi0_vec, true, file_name);
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

        int aa = 0;

        for (auto comp : abc_string) {
          pm.set_beta(freq_b, freq_c, comp,
                      beta_tensor(aa, 0));  // set the β value
          aa++;
        }
        // get perturbation direction
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
void compute_hyperpolarizability(
    World &world, const GroundStateData &gs,
    const std::vector<double> &frequencies,  // ωB and ωC
    const std::string &directions,           // which Cartesian dirs
    PropertyManager &pm) {
  // 1) Build A-list (dipole along each dir)
  std::vector<Perturbation> A_pert;
  A_pert.reserve(directions.size());
  for (char d : directions) A_pert.emplace_back(DipolePerturbation{d});

  // 2) Build all BC pairs of two dipoles
  std::pair<PerturbationType, PerturbationType> bc_types = {
      PerturbationType::Dipole, PerturbationType::Dipole};
  std::vector<std::pair<Perturbation, Perturbation>> BC_pairs;
  for (char b : directions)
    for (char c : directions)
      BC_pairs.emplace_back(DipolePerturbation{b}, DipolePerturbation{c});

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
void compute_Raman(
    World &world, const GroundStateData &gs,
    const std::pair<std::vector<double>, std::vector<double>> &BC_frequencies,
    const std::string &dip_dirs,  // e.g. {'x','y','z'}
    const std::string &nuc_dirs,  // e.g. {'x','y','z'}
    PropertyManager &pm) {
  // 1) A-list is still the dipole directions
  std::vector<Perturbation> A_pert;
  A_pert.reserve(dip_dirs.size());
  for (char d : dip_dirs) A_pert.emplace_back(DipolePerturbation{d});

  // 2) Enumerate all nuclear displacements (per atom, per axis)
  std::vector<Perturbation> nucs;
  nucs.reserve(3 * gs.molecule.natom());
  for (int atom = 0; atom < gs.molecule.natom(); ++atom) {
    for (char dir : {'x', 'y', 'z'})
      nucs.emplace_back(NuclearDisplacementPerturbation{atom, dir});
  }

  std::pair<PerturbationType, PerturbationType> bc_types = {
      PerturbationType::Dipole, PerturbationType::NuclearDisplacement};

  // 3) Build BC pairs = (Dipole, NuclearDisp)
  std::vector<std::pair<Perturbation, Perturbation>> BC_pairs;
  BC_pairs.reserve(dip_dirs.size() * nucs.size());
  for (char b : dip_dirs) {
    for (auto &n : nucs) {
      BC_pairs.emplace_back(DipolePerturbation{b},
                            std::get<NuclearDisplacementPerturbation>(n));
    }
  }

  // 4) Dispatch
  compute_beta(world, gs,
               PerturbationType::Dipole,  // the “A”-perturbation is always
               A_pert, bc_types, BC_pairs,
               BC_frequencies,  // ωB and ωC = same list of freqs
               pm);
}
