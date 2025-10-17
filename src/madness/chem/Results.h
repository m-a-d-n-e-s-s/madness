//
// Created by Florian Bischoff on 08.07.25.
//

#ifndef RESULTS_H
#define RESULTS_H

#include "madness/constants.h"
#include "madness_exception.h"
#include <madness/chem/molecule.h>
#include <madness/external/nlohmann_json/json.hpp>
#include <madness/tensor/tensor_json.hpp>
#include <string>

//* base and derived classes for holding results of a calculation
namespace madness {

class ResultsBase {

public:
  ResultsBase() = default;

  virtual ~ResultsBase() = default;
  /// serialize the results to a JSON object
  virtual nlohmann::json to_json() const = 0;
  virtual void from_json(const nlohmann::json &j) = 0;
  virtual std::string key() const = 0;
};

//--tiny helpers
template <class T, class F>
inline void set_if_exists(nlohmann::json &j, const std::string &key,
                          const std::optional<T> &opt, F &&to_json_fn) {
  if (opt)
    j[key] = to_json_fn(*opt);
}
template <class T, class F>
inline void get_if_exists(const nlohmann::json &j, const std::string &key,
                          std::optional<T> &opt, F &&from_json_fn) {
  if (j.contains(key))
    opt = from_json_fn(j[key]);
}

template <class T> inline nlohmann::json tensor_out(const Tensor<T> &t) {
  return tensor_to_json(t);
}
template <class T> inline Tensor<T> tensor_in(const nlohmann::json &j) {
  return tensor_from_json<T>(j);
}

/// holds metadata of the calculation

/// create right before the calculation starts, stop() must be called after the
/// calculation is finished
class MetaDataResults : public ResultsBase {
public:
  MetaDataResults(World &world) {
    time_begin = wall_time();
    mpi_size = world.size();
  }

  double time_begin = 0.0;
  double time_end = 0.0;
  std::string finished_at = "";
  std::string git_hash = "";
  int mpi_size = -1;
  std::string host = "";
  int nthreads = -1;

  std::string key() const override { return "metadata"; }

  void stop() {
    time_end = wall_time();
    finished_at = time_tag();
  }

  nlohmann::json to_json() const override {
    nlohmann::json j;
    // compute timing on-the-fly unless they have been set
    if (time_end == 0.0) {
      j["elapsed_time"] = wall_time() - time_begin;
      j["finished_at"] = time_tag();
    } else {
      j["elapsed_time"] = time_end - time_begin;
      j["finished_at"] = finished_at;
    }
    j["git_hash"] = git_hash;
    j["host"] = std::string(HOST_SYSTEM);
    j["nthreads"] = ThreadPool::size();
    j["mpi_size"] = mpi_size;
    return j;
  }

private:
  /// borrowed from Adrian's MolDFTLib
  std::string time_tag() const {
    auto print_time = std::chrono::system_clock::now();
    auto in_time_t = std::chrono::system_clock::to_time_t(print_time);
    std::stringstream ss;
    ss << std::put_time(std::localtime(&in_time_t), "%Y-%m-%d %X");
    return ss.str();
  }
};

/// holds convergence results of the calculation
class ConvergenceResults : public ResultsBase {
public:
  double converged_for_thresh = 1.e10;
  double converged_for_dconv = 1.e10;
  ConvergenceResults() = default;

  /// construct from JSON
  ConvergenceResults(const nlohmann::json &j) {
    converged_for_thresh = j.value("converged_for_thresh", 1.e10);
    converged_for_dconv = j.value("converged_for_dconv", 1.e10);
  }

  /// assignment operator from JSON
  ConvergenceResults &operator=(const nlohmann::json &j) {
    converged_for_thresh = j.value("converged_for_thresh", 1.e10);
    converged_for_dconv = j.value("converged_for_dconv", 1.e10);
    return *this;
  }

  std::string key() const override { return "convergence"; }

  ConvergenceResults &set_converged_thresh(double thresh) {
    converged_for_thresh = thresh;
    return *this;
  }

  ConvergenceResults &set_converged_dconv(double dconv) {
    converged_for_dconv = dconv;
    return *this;
  }

  nlohmann::json to_json() const override {
    nlohmann::json j;
    j["converged_for_thresh"] = converged_for_thresh;
    j["converged_for_dconv"] = converged_for_dconv;
    return j;
  }
  void from_json(const nlohmann::json &j) {
    converged_for_thresh = j.value("converged_for_thresh", 1.e10);
    converged_for_dconv = j.value("converged_for_dconv", 1.e10);
  }
};

class OptimizationResults : public ResultsBase {
public:
  int nsteps = 0;
  double final_energy = 0.0;
  double max_gradient = 0.0;
  double rms_gradient = 0.0;
  double max_step = 0.0;
  double rms_step = 0.0;
  madness::Molecule final_geometry;

  OptimizationResults() = default;

  /// construct from JSON
  OptimizationResults(const nlohmann::json &j) {
    nsteps = j.value("nsteps", 0);
    final_energy = j.value("final_energy", 0.0);
    max_gradient = j.value("max_gradient", 0.0);
    rms_gradient = j.value("rms_gradient", 0.0);
    max_step = j.value("max_step", 0.0);
    rms_step = j.value("rms_step", 0.0);
    if (j.contains("final_geometry"))
      final_geometry.from_json(j.at("final_geometry"));
  }

  std::string key() const override { return "optimization"; }

  [[nodiscard]] nlohmann::json to_json() const override {
    nlohmann::json j;
    j["nsteps"] = nsteps;
    j["final_energy"] = final_energy;
    j["max_gradient"] = max_gradient;
    j["rms_gradient"] = rms_gradient;
    j["max_step"] = max_step;
    j["rms_step"] = rms_step;
    j["final_geometry"] = final_geometry.to_json();
    return j;
  } // from json OptimizationResults

  void from_json(const nlohmann::json &j) override {
    // robust reads (won’t throw if missing)
    nsteps = j.value("nsteps", 0);
    final_energy = j.value("final_energy", 0.0);
    max_gradient = j.value("max_gradient", 0.0);
    rms_gradient = j.value("rms_gradient", 0.0);
    max_step = j.value("max_step", 0.0);
    rms_step = j.value("rms_step", 0.0);
    if (j.contains("final_geometry"))
      final_geometry.from_json(j.at("final_geometry"));
  }
};

class VibrationalResults : public ResultsBase {
public:
  std::optional<Tensor<double>> hessian;
  std::optional<Tensor<double>>
      frequencies; // (vibrational frequencies in a.u.)
  std::optional<Tensor<double>> intensities; //(IR intensities in km/mol)
  std::optional<Tensor<double>> reducedmass; //(reduced
  std::optional<Tensor<double>> normalmodes; //(normal modes)
                                             //
  static constexpr double au2invm =
      constants::au2invcm; // conversion factor from Hartree to cm^-1
  VibrationalResults() = default;
  VibrationalResults(const nlohmann::json &j) {}

  [[nodiscard]] bool has_data() const {
    return hessian || frequencies || intensities || reducedmass || normalmodes;
  }

  [[nodiscard]] std::string key() const override { return "vibrations"; }

  nlohmann::json to_json() const override {
    nlohmann::json j;

    set_if_exists(j, "hessian", hessian, tensor_out<double>);
    set_if_exists(j, "frequencies", frequencies, tensor_out<double>);
    set_if_exists(j, "intensities", intensities, tensor_out<double>);
    set_if_exists(j, "reducedmass", reducedmass, tensor_out<double>);
    set_if_exists(j, "normalmodes", normalmodes, tensor_out<double>);

    j["au2invcm"] = constants::au2invcm;

    return j;
  }

  void from_json(const nlohmann::json &j) override {
    get_if_exists(j, "hessian", hessian, tensor_in<double>);
    get_if_exists(j, "frequencies", frequencies, tensor_in<double>);
    get_if_exists(j, "intensities", intensities, tensor_in<double>);
    get_if_exists(j, "reducedmass", reducedmass, tensor_in<double>);
    get_if_exists(j, "normalmodes", normalmodes, tensor_in<double>);
  }
};

// A Raman calculation results in Raman Intensities for each normal mode
// computed at a given polarization frequency
//
class RamanResults : public ResultsBase {
public:
  std::string key() const override { return "raman"; }
  RamanResults() = default;

  std::vector<double> polarization_frequencies; // Polarization frequencies
  std::vector<double> vibrational_frequencies;  // Vibrational frequencies cm^-1
  Tensor<double> normal_modes;                  // Deriv. of alpha
  std::vector<Tensor<double>> polarizability_derivatives;
  std::vector<Tensor<double>> polarizability_derivatives_normal_modes;
  std::vector<Tensor<double>> frequencies; // Vibrational frequencies
  // square of the isotropic polarizability derivative
  std::vector<double> alpha2;
  // square of the anisotropic polarizability derivative
  std::vector<double> beta2;
  std::vector<double> intensities_raman;          // Raman intensities per mode
  std::vector<double> intensities_depolarization; // Raman intensities per mode
  std::vector<double> depolarization_ratios;      // Raman intensities per mode
  //
  void from_json(const nlohmann::json &j) override {
    if (j.contains("polarization_frequencies"))
      polarization_frequencies =
          j.at("polarization_frequencies").get<std::vector<double>>();
    if (j.contains("vibrational_frequencies"))
      vibrational_frequencies =
          j.at("vibrational_frequencies").get<std::vector<double>>();
    if (j.contains("polarizability_derivatives")) {
      polarizability_derivatives.clear();
      for (const auto &pd : j["polarizability_derivatives"])
        polarizability_derivatives.push_back(tensor_in<double>(pd));
    }
    if (j.contains("polarizability_derivatives_normal_modes")) {
      polarizability_derivatives_normal_modes.clear();
      for (const auto &pd : j["polarizability_derivatives_normal_modes"])
        polarizability_derivatives_normal_modes.push_back(
            tensor_in<double>(pd));
    }
    if (j.contains("frequencies")) {
      frequencies.clear();
      for (const auto &f : j["frequencies"])
        frequencies.push_back(tensor_in<double>(f));
    }
    if (j.contains("intensities_raman"))
      intensities_raman = j.at("intensities_raman").get<std::vector<double>>();
    if (j.contains("intensities_depolarization"))
      intensities_depolarization =
          j.at("intensities_depolarization").get<std::vector<double>>();
    if (j.contains("depolarization_ratios"))
      depolarization_ratios =
          j.at("depolarization_ratios").get<std::vector<double>>();
  }

  nlohmann::json to_json() const override {
    nlohmann::json j;
    j["polarization_frequencies"] = polarization_frequencies;
    j["vibrational_frequencies"] = vibrational_frequencies;
    j["polarizability_derivatives"] = nlohmann::json::array();
    for (const auto &pd : polarizability_derivatives)
      j["polarizability_derivatives"].push_back(tensor_out<double>(pd));
    j["polarizability_derivatives_normal_modes"] = nlohmann::json::array();
    for (const auto &pd : polarizability_derivatives_normal_modes)
      j["polarizability_derivatives_normal_modes"].push_back(
          tensor_out<double>(pd));
    j["frequencies"] = nlohmann::json::array();
    for (const auto &f : frequencies)
      j["frequencies"].push_back(tensor_out<double>(f));
    j["intensities_raman"] = intensities_raman;
    j["intensities_depolarization"] = intensities_depolarization;
    j["depolarization_ratios"] = depolarization_ratios;
    return j;
  }
};

class PropertyResults : public ResultsBase {
public:
  double energy = 0.0;

  std::optional<Tensor<double>> dipole;
  std::optional<Tensor<double>> gradient;
  std::optional<VibrationalResults> vibrations;
  std::optional<RamanResults> raman;

  PropertyResults() = default;

  /// construct from JSON
  PropertyResults(const nlohmann::json &j) {
    energy = j.value("energy", 0.0);
    if (j.count("dipole") == 1)
      dipole = tensor_from_json<double>(j["dipole"]);
    if (j.count("gradient") == 1)
      gradient = tensor_from_json<double>(j["gradient"]);
  }

  std::string key() const override { return "properties"; }

  nlohmann::json to_json() const override {
    nlohmann::json j;
    j["energy"] = energy;
    set_if_exists(j, "dipole", dipole, tensor_out<double>);
    set_if_exists(j, "gradient", gradient, tensor_out<double>);
    if (vibrations && vibrations->has_data())
      j["vibrations"] = vibrations->to_json();
    return j;
  } // from json PropertyResults

  void from_json(const nlohmann::json &j) override {
    // robust reads (won’t throw if missing)
    if (j.contains("energy"))
      energy = j.value("energy", 0.0);
    get_if_exists(j, "dipole", dipole, tensor_in<double>);
    get_if_exists(j, "gradient", gradient, tensor_in<double>);

    // nested section
    if (j.contains("vibrations")) {
      VibrationalResults vib;
      vib.from_json(j.at("vibrations"));
      if (vib.has_data())
        vibrations = std::move(vib);
    }
  }
};

// If you keep PropertyResults from earlier, give it:
inline bool has_data(const PropertyResults &p) {
  return p.energy != 0.0 || p.dipole || p.gradient ||
         (p.vibrations && p.vibrations->has_data());
}

class SCFResults : public ResultsBase {
public:
  // Required alpha (for RHF/ ROHF/ UHF UKS we alsways expect alpha)
  Tensor<double> aeps;
  Tensor<double> afock;
  Molecule scf_molecule;
  // optional beta (only for UHF/ UKS)
  std::optional<Tensor<double>> beps;
  std::optional<Tensor<double>> bfock;
  bool is_opt = false;

  std::string model = "scf";     // model used for the SCF calculation
  double scf_total_energy = 0.0; // total energy of the SCF calculation
  //
  PropertyResults properties;
  SCFResults() = default;

  /// construct from JSON
  SCFResults(const nlohmann::json &j) { from_json(j); }

  std::string key() const override { return model; }

  nlohmann::json to_json() const override {
    nlohmann::json j;

    // Required alpha pieces
    j["scf_eigenvalues_a"] = tensor_out<double>(aeps);
    j["scf_fock_a"] = tensor_out<double>(afock);

    // Optional beta pieces
    set_if_exists(j, "scf_eigenvalues_b", beps, tensor_out<double>);
    set_if_exists(j, "scf_fock_b", bfock, tensor_out<double>);

    // Scalars / metadata
    j["model"] = model;
    j["scf_total_energy"] = scf_total_energy;

    // Optional nested block
    if (has_data(properties)) {
      j["properties"] = properties.to_json();
    }

    j["molecule"] = scf_molecule.to_json();
    j["is_opt"] = is_opt;
    return j;
  }

  void from_json(const nlohmann::json &j) override {
    // Alpha: treat as required but read defensively
    if (j.contains("scf_eigenvalues_a"))
      aeps = tensor_in<double>(j.at("scf_eigenvalues_a"));
    else
      aeps = {}; // or throw if truly required

    if (j.contains("scf_fock_a"))
      afock = tensor_in<double>(j.at("scf_fock_a"));
    else
      afock = {}; // or throw if truly required

    // Beta: optional
    get_if_exists(j, "scf_eigenvalues_b", beps, tensor_in<double>);
    get_if_exists(j, "scf_fock_b", bfock, tensor_in<double>);

    // Scalars / metadata
    if (j.contains("model"))
      model = j.value("model", std::string("scf"));
    if (j.contains("scf_total_energy"))
      scf_total_energy = j.value("scf_total_energy", 0.0);

    // Nested properties: optional
    if (j.contains("properties")) {
      PropertyResults p;
      p.from_json(j.at("properties"));
      if (has_data(p))
        properties = std::move(p);
    } else {
      properties = PropertyResults();
    }
    if (j.contains("molecule"))
      scf_molecule.from_json(j.at("molecule"));
    else
      MADNESS_EXCEPTION("Missing molecule data", j);
    is_opt = j.value("is_opt", false);
  }
};
// Todo: Upgrade to new JSON style using optional everything below here --
//
// ---------------------------------------------------------------------------------
//
//
// ---------------------------------------------------------------------------------
class CISResults : public ResultsBase {
public:
  struct excitation_info {
    std::string irrep;                   // irreducible representation
    double omega;                        // excitation energy in Hartree
    double current_error;                // error in the excitation energy
    double oscillator_strength_length;   // oscillator strength
    double oscillator_strength_velocity; // oscillator strength
  };
  std::vector<excitation_info> excitations;
  long nfreeze = -1;
  std::string model = "unknown";

  std::string key() const override { return model; }

  CISResults() = default;

  /// construct from JSON
  CISResults(const nlohmann::json &j) {
    if (j.count("excitations") > 0) {
      for (const auto &ex : j["excitations"]) {
        excitation_info ei;
        ei.irrep = ex.value("irrep", "");
        ei.omega = ex.value("omega", 0.0);
        ei.current_error = ex.value("current_error", 0.0);
        ei.oscillator_strength_length =
            ex.value("oscillator_strength_length", 0.0);
        ei.oscillator_strength_velocity =
            ex.value("oscillator_strength_velocity", 0.0);
        excitations.push_back(ei);
      }
    }
    nfreeze = j.value("nfreeze", -1);
    model = j.value("model", "unknown");
  }

  /// constructor with nfreeze and model
  CISResults(long nfreeze, const std::string &model)
      : nfreeze(nfreeze), model(model) {}

  void from_json(const nlohmann::json &j) override {
    excitations.clear();
    if (j.count("excitations") > 0) {
      for (const auto &ex : j["excitations"]) {
        excitation_info ei;
        ei.irrep = ex.value("irrep", "");
        ei.omega = ex.value("omega", 0.0);
        ei.current_error = ex.value("current_error", 0.0);
        ei.oscillator_strength_length =
            ex.value("oscillator_strength_length", 0.0);
        ei.oscillator_strength_velocity =
            ex.value("oscillator_strength_velocity", 0.0);
        excitations.push_back(ei);
      }
    }
    nfreeze = j.value("nfreeze", -1);
    model = j.value("model", "unknown");
  }

  nlohmann::json to_json() const override {
    nlohmann::json j;
    for (const auto &ex : excitations) {
      nlohmann::json ex_json;
      ex_json["irrep"] = ex.irrep;
      ex_json["omega"] = ex.omega;
      ex_json["current_error"] = ex.current_error;
      ex_json["oscillator_strength_length"] = ex.oscillator_strength_length;
      ex_json["oscillator_strength_velocity"] = ex.oscillator_strength_velocity;
      j["excitations"].push_back(ex_json);
    }
    j["nfreeze"] = nfreeze;
    j["model"] = model;
    return j;
  }
};

class CC2Results : public CISResults {
public:
  CC2Results() : CISResults() { model = "mp2"; }

  PropertyResults properties; // properties of the correlated calculation
  double correlation_energy =
      0.0;                   // correlation energy of the correlated calculation
  double total_energy = 0.0; // total energy of the correlated calculation
  /// construct from JSON
  CC2Results(const nlohmann::json &j) : CISResults(j) {
    properties = PropertyResults(j.value("properties", nlohmann::json{}));
    model = j.value("model", "mp2");
    correlation_energy = j.value("correlation_energy", 0.0);
    total_energy = j.value(model + "_total_energy", 0.0);
  }

  /// constructor with nfreeze and model
  CC2Results(long nfreeze, const std::string &model)
      : CISResults(nfreeze, model) {}

  void from_json(const nlohmann::json &j) override {
    CISResults::from_json(j);
    properties = PropertyResults(j.value("properties", nlohmann::json{}));
    model = j.value("model", "mp2");
    correlation_energy = j.value("correlation_energy", 0.0);
    total_energy = j.value(model + "_total_energy", 0.0);
  }

  nlohmann::json to_json() const override {
    nlohmann::json j;
    j = CISResults::to_json();
    j["properties"] = properties.to_json();
    j["model"] = model;
    j["correlation_energy"] = correlation_energy;
    j[model + "_correlation_energy"] = correlation_energy;
    j[model + "_total_energy"] = total_energy;
    return j;
  }

  CC2Results &set_energies(const double scf_energy, const double corr_energy) {
    this->correlation_energy = corr_energy;
    this->total_energy = scf_energy + corr_energy;
    return *this;
  }

  /// setters with chaining
  CC2Results &set_correlation_energy(const double corr_energy) {
    correlation_energy = corr_energy;
    return *this;
  }
  CC2Results &set_total_energy(const double total_energy) {
    this->total_energy = total_energy;
    return *this;
  }
  CC2Results &set_properties(const PropertyResults &props) {
    properties = props;
    return *this;
  }
  CC2Results &set_model(const std::string &model) {
    this->model = model;
    return *this;
  }
};

class ZnemoResults : public SCFResults {
public:
  double B = 0.0; // B value for the Znemo calculation

  ZnemoResults() = default;
  /// construct from JSON
  ZnemoResults(const nlohmann::json &j) : SCFResults(j) {
    B = j.value("B", 0.0);
  }

  void from_json(const nlohmann::json &j) override {
    SCFResults::from_json(j);
    B = j.value("B", 0.0);
  }

  nlohmann::json to_json() const override {
    nlohmann::json j;
    j = SCFResults::to_json();
    j["B"] = B;
    return j;
  }
};

class OEPResults : public SCFResults {
public:
  double drho = 0.0;    // delta rho =difference to reference (=HF?) density
  double devir14 = 0.0; // diagnostic parameter
  double devir17 = 0.0; // diagnostic parameter
  double Ex_vir = 0.0;  // local exchange energy
  double Ex_conv = 0.0; //
  double Ex_HF = 0.0;   // HF exchange energy
  double E_kin_HF = 0.0;
  double E_kin_KS = 0.0; // kinetic energy of the KS reference
  double Econv = 0.0;    // final energy using conventional method

  OEPResults() = default;

  void from_json(const nlohmann::json &j) override {
    SCFResults::from_json(j);
    model = j.value("model", "oaep");
    drho = j.value("drho", 0.0);
    devir14 = j.value("devir14", 0.0);
    devir17 = j.value("devir17", 0.0);
    Ex_vir = j.value("Ex_vir", 0.0);
    Ex_conv = j.value("Ex_conv", 0.0);
    Ex_HF = j.value("Ex_HF", 0.0);
    E_kin_HF = j.value("E_kin_HF", 0.0);
    E_kin_KS = j.value("E_kin_KS", 0.0);
    Econv = j.value("Econv", 0.0);
  }

  /// construct from JSON
  explicit OEPResults(const nlohmann::json &j) : SCFResults(j) {
    model = j.value("model", "oaep");
    drho = j.value("drho", 0.0);
    devir14 = j.value("dvir14", 0.0);
    devir17 = j.value("dvir17", 0.0);
    Ex_vir = j.value("Ex_vir", 0.0);
    Ex_conv = j.value("Ex_conv", 0.0);
    Ex_HF = j.value("Ex_HF", 0.0);
    E_kin_HF = j.value("E_kin_HF", 0.0);
    E_kin_KS = j.value("E_kin_KS", 0.0);
    Econv = j.value("Econv", 0.0);
  }

  nlohmann::json to_json() const override {
    nlohmann::json j;
    j = SCFResults::to_json();
    j["model"] = model;
    j["drho"] = drho;
    j["devir14"] = devir14;
    j["devir17"] = devir17;
    j["Ex_vir"] = Ex_vir;
    j["Ex_conv"] = Ex_conv;
    j["Ex_HF"] = Ex_HF;
    j["E_kin_HF"] = E_kin_HF;
    j["E_kin_KS"] = E_kin_KS;
    j["Econv"] = Econv;
    return j;
  }
};

using SCFResultsTuple =
    std::tuple<SCFResults, PropertyResults, ConvergenceResults>;

} // namespace madness
#endif // RESULTS_H
