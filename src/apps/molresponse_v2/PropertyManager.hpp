#pragma once
#include "../../madness/chem/ResponseParameters.hpp"
#include "InnerContributions.hpp"
#include "MolecularProperty.hpp"
#include "Perturbation.hpp"
#include "ResponseIO.hpp"
#include "ResponseSolverUtils.hpp"
#include "ResponseVector.hpp"
#include "Results.h"
#include "VBCMacrotask.hpp"
#include "broadcast_json.hpp"
#include "functypedefs.h"

#include <chrono>
#include <filesystem>
#include <format>
#include <fstream>
#include <iomanip>
#include <madness/chem/vibanal.h>
#include <madness/external/nlohmann_json/json.hpp>
#include <madness/tensor/tensor.h>
#include <madness/tensor/tensor_json.hpp>
#include <madness/world/world.h>
#include <string>

using json = nlohmann::json;

using components = std::vector<std::string>;
// handle printing  components as vector<string> for generality
std::ostream &operator<<(std::ostream &s, const components &c) {
  s << "[";
  for (size_t i = 0; i < c.size(); ++i) {
    s << c[i];
    if (i != c.size() - 1)
      s << ",";
  }
  s << "]";
  return s;
}

struct PropRow {
  std::string property;
  components component;
  double freq1;
  std::optional<double> freq2;
  std::optional<double> value;
};

struct PropKey {
  std::string property;
  components component;
  double freq1;
  std::optional<double> freq2;

  bool operator<(const PropKey &other) const noexcept {
    if (property != other.property) {
      return property < other.property;
    }
    if (component != other.component) {

      return component < other.component;
    }
    if (freq1 != other.freq1) {
      return freq1 < other.freq1;
    }
    // for freq2
    if (!freq2 && other.freq2) {
      return true;
    }
    if (freq2 && !other.freq2) {

      return false;
    }
    if (freq2 && other.freq2) {

      return *freq2 < *other.freq2;
    }
    return false;
  }
};

// enable Json to PropRow
inline void to_json(json &j, const PropRow &r) {
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

inline void from_json(const json &j, PropRow &r) {
  r.property = j.at("property").get<std::string>();
  r.component = j.at("component").get<std::vector<std::string>>();
  r.freq1 = j.at("freqB").get<double>();
  if (j.contains("freqC")) {
    r.freq2 = j.at("freqC").get<double>();
  }
  if (j.contains("value")) {
    r.value = j.at("value").get<double>();
  }
}

enum class PropertyType { Alpha, Beta, Raman };

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
    for (size_t i = 0; i < nA; ++i)
      Ak[i] = A_vecs[i][k];
    for (size_t j = 0; j < nB; ++j)
      Bk[j] = B_vecs[j][k];

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

class PropertyManager {
public:
  explicit PropertyManager(World &world, const std::string &filename)
      : filename_(filename) {
    if (fs::exists(filename_)) {
      // load JSON
      json j = broadcast_json_file(world, filename_);
      if (j.is_array()) {
        for (const auto &r : j) {
          PropRow row = r.get<PropRow>();
          PropKey key = {row.property, row.component, row.freq1, row.freq2};
          rows_[key] = std::move(row);
          // flat‐format: just parse rows
        }
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
  [[nodiscard]] bool has_alpha(double omega, components comp) const {
    PropKey k{"polarizability", comp, omega, std::nullopt};
    auto it = rows_.find(k);
    return it != rows_.end() && it->second.value.has_value();
    return rows_.count(k) != 0;
  }

  [[nodiscard]] bool has_beta(double w1, double w2, components comp) const {
    PropKey k{"hyperpolarizability", comp, w1, w2};
    auto it = rows_.find(k);
    return it != rows_.end() && it->second.value.has_value();
  }

  [[nodiscard]] bool has_raman(double w1, double w2, components comp) const {
    PropKey k{"raman", comp, w1, w2};
    auto it = rows_.find(k);
    return it != rows_.end() && it->second.value.has_value();
  }

  [[nodiscard]] std::optional<double> get_raman(double w1, double w2,
                                                components comp) const {
    PropKey k{"raman", comp, w1, w2};
    auto it = rows_.find(k);
    if (it != rows_.end() && it->second.value.has_value()) {
      return it->second.value;
    }
    return std::nullopt;
  }

  [[nodiscard]] std::optional<double> get_beta(double w1, double w2,
                                               components comp) const {
    PropKey k{"hyperpolarizability", comp, w1, w2};
    auto it = rows_.find(k);
    if (it != rows_.end() && it->second.value.has_value()) {
      return it->second.value;
    }
    return std::nullopt;
  }

  [[nodiscard]] std::optional<double> get_alpha(double omega,
                                                components comp) const {
    PropKey k{"polarizability", comp, omega, std::nullopt};
    auto it = rows_.find(k);
    if (it != rows_.end() && it->second.value.has_value()) {
      return it->second.value;
    }
    return std::nullopt;
  }

  // Insert or overwrite α entries
  void set_alpha(double omega, const madness::Tensor<double> &tensor,
                 const std::string &dirs) {
    size_t N = dirs.size();
    for (size_t i = 0; i < N; ++i) {
      for (size_t j = 0; j < N; ++j) {
        components comp = {std::format("{}", dirs[i]),
                           std::format("{}", dirs[j])};
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
  void set_beta(double w1, double w2, const components &comp, double value) {
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

  // Insert or overwrite β entries
  void set_raman(double w1, double w2, const components &comp, double value) {
    PropKey k{"raman", comp, w1, w2};
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
    for (const auto &[k, r] : rows_) {
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
};

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

  // if (world.rank() == 0) {
  //   print("▶️ Computing α tensor for", num_directions, "directions and",
  //         num_frequencies, "frequencies.");
  // }
  //
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
          world, gs, state_map.at(dir))); // Get the perturbation vector
    }
  }

  for (size_t ffi = 0; ffi < num_frequencies; ++ffi) {
    double omega = frequencies[ffi];
    // if (pm.has_alpha(omega, directions))
    // {
    //   if (world.rank() == 0)
    //     print("✅ Skipping already computed α at ω =", omega);
    //   continue;
    // }

    // TODO: create unrestricted and restricted versions of the properties
    double alpha_factor = (omega == 0.0) ? -4.0 : -2.0;
    // if (world.rank() == 0) {
    //   print("🛠️  Computing α at ω =", omega, "for directions:", directions,
    //         " alpha factor = ", alpha_factor);
    // }

    std::vector<vector_real_function_3d> response_vecs(num_directions);
    std::vector<vector_real_function_3d> perturb_vecs(num_directions);

    for (size_t j = 0; j < num_directions; ++j) {
      auto &active_state = state_map.at(direction_keys[j]);
      active_state.set_frequency_index(ffi);
      load_response_vector(world, num_orbitals, active_state,
                           active_state.thresholds.size() - 1, ffi,
                           load_vector[j]);
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
        "alpha_contribs_" + std::to_string(ffi));

    alpha *= alpha_factor;
    if (world.rank() == 0) {
      print("α tensor size:", alpha.size(), "x", alpha.size());
      print("α= ", alpha);
    }

    pm.set_alpha(omega, alpha, directions);
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
[[nodiscard]] VibrationalResults
compute_hessian(World &world,
                std::map<std::string, LinearResponseDescriptor> &state_map,
                const GroundStateData &gs, const std::string &directions,
                const std::shared_ptr<SCF> &scf_calc) {
  const size_t num_directions = directions.size();
  const size_t num_orbitals = gs.getNumOrbitals();
  const bool is_restricted = gs.isSpinRestricted();

  auto mol = gs.getMolecule();
  auto natom = mol.natom();
  std::vector<std::string> atom_derivative_keys;
  for (int i = 0; i < natom; ++i) {
    for (const char &dir : directions) {
      NuclearDisplacementPerturbation n{i, dir};
      atom_derivative_keys.push_back(describe_perturbation(n));
    }
  }

  std::vector<vector_real_function_3d> perturbations;
  for (const auto &adx : atom_derivative_keys) {
    if (state_map.find(adx) == state_map.end()) {
      throw std::runtime_error("State not found in state_map: " + adx);
    } else {
      perturbations.push_back(perturbation_vector(
          world, gs, state_map.at(adx))); // Get the perturbation vector
    }
  }

  const auto &phi0 = gs.getOrbitals();
  const auto &occ = gs.getOcc();

  auto rho0 = 2.0 * scf_calc->make_density(world, occ, phi0);
  auto drhoX = vector_real_function_3d(3);
  for (int axis = 0; axis < 3; ++axis) {
    real_derivative_3d D = free_space_derivative<double, 3>(world, axis);
    real_function_3d rho_copy = copy(rho0).refine();
    drhoX[axis] = D(rho_copy);
  } //{dx,dy,dz}
  //

  Tensor<double> hessian(natom * 3, natom * 3);
  // double alpha_factor = (omega == 0.0) ? -4.0 : -2.0;
  for (size_t iatom = 0; iatom < natom; ++iatom) {
    for (int iaxis = 0; iaxis < 3; ++iaxis) {
      int i = iatom * 3 + iaxis;

      auto &active_state = state_map.at(atom_derivative_keys[i]);
      if (world.rank() == 0) {
        print("Loading response vector for ", atom_derivative_keys[i]);
      }

      ResponseVector load_vector;
      active_state.set_frequency_index(0);
      load_response_vector(world, num_orbitals, active_state,
                           active_state.thresholds.size() - 1, 0, load_vector);
      auto xi = get_flat(load_vector);
      auto xphi = mul(world, xi, phi0, true);
      auto rho_xi = 4.0 * sum(world, xphi, true);

      for (size_t jatom = 0; jatom < natom; ++jatom) {
        for (int jaxis = 0; jaxis < 3; ++jaxis) {
          int j = jatom * 3 + jaxis;

          // skip diagonal elements because they are extremely noisy!
          // use translational symmetry to reconstruct them from other
          // hessian matrix elements (see below)
          if (i == j)
            continue;

          madchem::MolecularDerivativeFunctor mdf(mol, jatom, jaxis);
          hessian(i, j) = inner(rho_xi, mdf);

          // integration by parts < 0 | H^{YX} | 0 >
          if (iatom == jatom)
            hessian(i, j) += inner(drhoX[iaxis], mdf);
        }
      }
    }
  }
  VibrationalResults results;
  if (world.rank() == 0) {
    print("\n raw electronic Hessian (a.u.)\n");
    print(hessian);
  }
  auto purify_hessian = [&](const Tensor<double> &H) {
    Tensor<double> purified = copy(H);
    double maxasymmetric = 0.0;

    const size_t natom = mol.natom();

    for (size_t iatom = 0; iatom < natom; ++iatom) {
      for (int iaxis = 0; iaxis < 3; ++iaxis) {
        int i = iatom * 3 + iaxis;

        for (size_t jatom = 0; jatom < natom; ++jatom) {
          for (int jaxis = 0; jaxis < 3; ++jaxis) {
            int j = jatom * 3 + jaxis;

            double mean = (purified(i, j) + purified(j, i)) * 0.5;
            double diff = 0.5 * fabs(purified(i, j) - purified(j, i));
            maxasymmetric = std::max(maxasymmetric, diff);

            unsigned int ZA = mol.get_atomic_number(iatom);
            unsigned int ZB = mol.get_atomic_number(jatom);
            if (ZA < ZB)
              purified(i, j) = purified(j, i);
            if (ZA > ZB)
              purified(j, i) = purified(i, j);
            if (ZA == ZB) {
              purified(i, j) = mean;
              purified(j, i) = mean;
            }
          }
        }
      }
    }
    return purified;
  };

  for (size_t i = 0; i < 3 * natom; ++i)
    hessian(i, i) = 0.0;

  hessian = purify_hessian(hessian);
  Tensor<double> asymmetric = 0.5 * (hessian - transpose(hessian));
  const double max_asymmetric = asymmetric.absmax();
  // symmetrize hessian
  hessian += transpose(hessian);
  hessian.scale(0.5);
  // exploit translational symmetry to compute the diagonal elements:
  // translating all atoms in the same direction will make no energy change,
  // therefore the respective sum of hessian matrix elements will be zero:
  for (size_t i = 0; i < 3 * natom; ++i) {
    double sum = 0.0;
    for (size_t j = 0; j < 3 * natom; j += 3)
      sum += hessian(i, j + (i % 3));
    hessian(i, i) = -sum;
  }

  if (world.rank() == 0) {
    print("\n electronic Hessian (a.u.)\n");
    print(hessian);
  }

  // add the nuclear-nuclear contribution
  hessian += mol.nuclear_repulsion_hessian();

  //    if (hessdebug) {
  if (world.rank() == 0) {
    print("\n Hessian (a.u.)\n");
    print(hessian);
  }

  Tensor<double> normalmodes;
  Tensor<double> vib_freq = compute_frequencies(
      mol, hessian, normalmodes, false, true and world.rank() == 0);

  if (true and world.rank() == 0) {
    print("\n vibrational frequencies (unprojected) (a.u.)\n");
    print(vib_freq);
    print("\n vibrational frequencies (unprojected) (cm-1)\n");
    print(constants::au2invcm * vib_freq);
  }

  vib_freq = compute_frequencies(mol, hessian, normalmodes, true,
                                 true && world.rank() == 0);
  //  Tensor<double> intensities = compute_IR_intensities(normalmodes, dens_pt);
  Tensor<double> reducedmass = compute_reduced_mass(mol, normalmodes);

  if (world.rank() == 0) {
    print("\nprojected vibrational frequencies (cm-1)\n");
    printf("frequency in cm-1   ");
    for (int i = 0; i < vib_freq.size(); ++i) {
      printf("%10.3f", constants::au2invcm * vib_freq(i));
    }
    // printf("\n");
    // printf("intensity in km/mol ");
    // for (int i = 0; i < intensities.size(); ++i) {
    //   printf("%10.3f", intensities(i));
    // }
    // printf("\n");
    // printf("reduced mass in amu ");
    // for (int i = 0; i < intensities.size(); ++i) {
    //   printf("%10.3f", reducedmass(i));
    // }
    // printf("\n\n");
    // printf("done with computing the hessian matrix at time %8.1fs \n",
    // wall_time()); printf("final energy %16.8f", calc->current_energy);
  }
  results.hessian = hessian;
  results.frequencies = vib_freq;
  // results.intensities = intensities;
  results.reducedmass = reducedmass;
  results.normalmodes = normalmodes;
  return results;
}

/// @brief
/// @param world
/// @param gs
/// @param a_type
/// @param perturbation_A
/// @param bc_types
/// @param BC_pairs
/// @param frequencies
/// @param pm
template <typename BType, typename CType>
void compute_beta(
    World &world, const GroundStateData &gs,
    const std::vector<Perturbation> &perturbation_A,
    const std::vector<std::pair<Perturbation, Perturbation>> &BC_pairs,
    const std::pair<std::vector<double>, std::vector<double>> &frequencies,
    PropertyManager &pm, const PropertyType &prop_type) {
  const bool is_spin_restricted = gs.isSpinRestricted();
  const int num_orbitals = static_cast<int>(gs.getNumOrbitals());
  const double thresh = FunctionDefaults<3>::get_thresh();

  // 1) Build a SimpleVBCComputer once
  auto vbc_computer = SimpleVBCComputer(world, gs);
  using componets = std::vector<std::string>;

  // 2) Loop over all freq
  for (auto &freq_b : frequencies.first) {
    for (auto &freq_c : frequencies.second) {
      // Get all BC pairs possible out of Perturbation
      for (auto [B, C] : BC_pairs) {
        VBCResponseState vbc_state(B, C, freq_b, freq_c,
                                   FunctionDefaults<3>::get_thresh(),
                                   is_spin_restricted);
        auto bc = vbc_state.perturbationDescription();
        auto pertB = std::get<BType>(B); // get the B perturbation
        auto pertC = std::get<CType>(C); // get the C perturbation

        std::vector<components> abc_string;
        for (auto a : perturbation_A) {
          auto pertA = std::get<DipolePerturbation>(a);
          auto dir = pertA.direction;

          components comp = {describe_perturbation(a), describe_perturbation(B),
                             describe_perturbation(C)};
          abc_string.push_back(comp);
        }

        bool skip = true;
        // skip if all abc are already computed

        for (auto &abc : abc_string) {
          // if PropType
          if (prop_type == PropertyType::Beta) {
            if (!pm.has_beta(freq_b, freq_c, abc)) {
              skip = false;
              break;
            }
          } else if (prop_type == PropertyType::Raman) {
            if (!pm.has_raman(freq_b, freq_c, abc)) {
              skip = false;
              break;
            }
          }
        }
        if (skip) {
          // if (world.rank() == 0)
          //   print("✅ Skipping already computed β at ωB =", freq_b,
          //         "ωC =", freq_c, "for directions:", abc_string);
          continue;
        }

        auto vbc_vec = vbc_computer.compute_and_save(vbc_state);
        auto vbc_flat = get_flat(vbc_vec);

        auto omega_A = vbc_state.current_frequency();
        std::vector<ResponseVector> a_vecs(perturbation_A.size());
        std::vector<vector_real_function_3d> xa_vecs(
            perturbation_A.size()); // for the A perturbations
        //

        std::vector<real_function_3d> opAs(perturbation_A.size());
        for (int a = 0; a < perturbation_A.size(); ++a) {
          auto pertA = perturbation_A[a];
          auto state_A = LinearResponseDescriptor(
              pertA, {omega_A}, {FunctionDefaults<3>::get_thresh()},
              is_spin_restricted);
          opAs[a] = raw_perturbation_operator(world, gs, state_A.perturbation);
          load_response_vector(world, num_orbitals, state_A, 0, 0, a_vecs[a]);
          auto flat = get_flat(a_vecs[a]);
          if (omega_A == 0.0) {
            flat.insert(flat.end(), flat.begin(),
                        flat.end()); // replicate for static
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
            vbc_computer.get_BC_vecs(vbc_state); // get the B and C states

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

        beta_tensor *= -2.0; // factor of -2 for beta
        if (world.rank() == 0) {
          print(vbc_state.perturbationDescription());
          print("β tensor size:", beta_tensor.size(), "x", beta_tensor.size());
          print("β= \n", beta_tensor);
        }

        int aa = 0;

        for (auto comp : abc_string) {
          if (prop_type == PropertyType::Raman) {

            pm.set_raman(freq_b, freq_c, comp,
                         beta_tensor(aa, 0)); // set the Raman value
          } else if (prop_type == PropertyType::Beta) {
            pm.set_beta(freq_b, freq_c, comp,
                        beta_tensor(aa, 0)); // set the β value
          }
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
    const std::vector<double> &frequencies, // ωB and ωC
    const std::string &directions,          // which Cartesian dirs
    PropertyManager &pm) {
  // 1) Build A-list (dipole along each dir)
  std::vector<Perturbation> A_pert;
  A_pert.reserve(directions.size());
  for (char d : directions)
    A_pert.emplace_back(DipolePerturbation{d});

  // 2) Build all BC pairs of two dipoles
  std::vector<std::pair<Perturbation, Perturbation>> BC_pairs;
  for (char b : directions)
    for (char c : directions)
      BC_pairs.emplace_back(DipolePerturbation{b}, DipolePerturbation{c});

  auto freq_pair = std::make_pair(frequencies, frequencies);

  // 3) Dispatch
  compute_beta<DipolePerturbation, DipolePerturbation>(
      world, gs, A_pert, BC_pairs,
      freq_pair, // ωB and ωC = same list of freqs
      pm, PropertyType::Beta);
}

// ─────────────────────────────────────────────────────────────────────────────
// Compute the Raman response (∂α/∂Q) via β(Dipole; Dipole, NuclearDisp).
//  frequencies: the optical frequencies for both Dipole and Raman.
//  dip_dirs: which dipole directions to include.
// ─────────────────────────────────────────────────────────────────────────────
std::vector<Tensor<double>>
compute_Raman(World &world, const GroundStateData &gs,
              const std::vector<double> &frequencies, // ωB and ωC
              const std::vector<int> &nuc_indices,
              const std::string &dip_dirs, // e.g. {'x','y','z'}
              const std::string &nuc_dirs, // e.g. {'x','y','z'}
              PropertyManager &pm) {
  // 1) A-list is still the dipole directions
  std::vector<Perturbation> A_pert;
  A_pert.reserve(dip_dirs.size());
  for (char d : dip_dirs)
    A_pert.emplace_back(DipolePerturbation{d});

  // 3) Build BC pairs = (Dipole, NuclearDisp)
  std::vector<std::pair<Perturbation, Perturbation>> BC_pairs;

  BC_pairs.reserve(dip_dirs.size() * nuc_indices.size() * nuc_dirs.size());
  for (char b : dip_dirs) {
    for (auto &n : nuc_dirs) {
      for (int index : nuc_indices) {
        BC_pairs.emplace_back(DipolePerturbation{b},
                              NuclearDisplacementPerturbation{index, n});
      }
    }
  }
  std::pair<std::vector<double>, std::vector<double>> BC_frequencies;
  BC_frequencies.first = frequencies;
  // zeros size of second vector
  BC_frequencies.second = std::vector<double>(frequencies.size(), 0.0);

  // 4) Dispatch
  compute_beta<DipolePerturbation, NuclearDisplacementPerturbation>(
      world, gs, A_pert, BC_pairs,
      BC_frequencies, // ωB and ωC = same list of freqs
      pm, PropertyType::Raman);

  std::vector<components> all_comps;

  // Let's control how we want the tensor to be shaped
  //
  //
  //  We want all [XX,XY,XZ] in the rows
  //  And we want all nuclear perturbations in C
  //
  //
  //
  std::vector<Tensor<double>> raman_tensors;
  for (auto freq : frequencies) {
    Tensor<double> raman_tensor(dip_dirs.size() * dip_dirs.size(),
                                nuc_indices.size() * nuc_dirs.size());
    int aa = 0;
    for (char a : dip_dirs) {
      for (char b : dip_dirs) {
        int cc = 0;
        for (int index : nuc_indices) {
          for (auto &n : nuc_dirs) {
            auto a_desc = describe_perturbation(DipolePerturbation{a});
            auto b_desc = describe_perturbation(DipolePerturbation{b});
            auto c_desc = describe_perturbation(
                NuclearDisplacementPerturbation{index, n});
            components comp = {a_desc, b_desc, c_desc};
            all_comps.push_back(comp);
            double val = *pm.get_raman(frequencies[0], 0.0, comp);

            if (val) {
              raman_tensor(aa, cc) = val;
            } else {
              raman_tensor(aa, cc) = 0.0;
            }
            cc++;
          }
        }
        aa++;
      }
    }
    raman_tensors.push_back(raman_tensor);
  }
  return raman_tensors;
}
