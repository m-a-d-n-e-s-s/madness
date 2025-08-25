//
// Created by Florian Bischoff on 08.07.25.
//

#pragma once
#ifndef RESULTS_H
#define RESULTS_H

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

  virtual std::string key() const = 0;
};

/// holds metadata of the calculation

/// create right before the calculation starts, stop() must be called after the calculation is finished
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

class PropertyResults : public ResultsBase {
public:
  PropertyResults() = default;

  /// construct from JSON
  PropertyResults(const nlohmann::json &j) {
    energy = j.value("energy", 0.0);
    if (j.count("dipole") == 1)
      dipole = tensor_from_json<double>(j["dipole"]);
    if (j.count("gradient") == 1)
      gradient = tensor_from_json<double>(j["gradient"]);
  }

  double energy = 0.0;
  Tensor<double> dipole;
  Tensor<double> gradient;

  std::string key() const override { return "properties"; }

  nlohmann::json to_json() const override {
    nlohmann::json j;
    j["energy"] = energy;
    if (dipole.size() > 0)
      j["dipole"] = tensor_to_json(dipole);
    if (gradient.size() > 0)
      j["gradient"] = tensor_to_json(gradient);
    return j;
  }
  void from_json(const nlohmann::json &j) {

    if (j.contains("energy"))
      energy = j.value("energy", 0.0);
    if (j.contains("dipole"))
      dipole = tensor_from_json<double>(j["dipole"]);
    if (j.contains("gradient"))
      gradient = tensor_from_json<double>(j["gradient"]);
  }
};
class SCFResults : public ResultsBase {
public:
  Tensor<double> aeps;
  Tensor<double> beps;
  Tensor<double> afock;
  Tensor<double> bfock;
  std::string model = "scf";     // model used for the SCF calculation
  double scf_total_energy = 0.0; // total energy of the SCF calculation
  PropertyResults properties;

  SCFResults() = default;

  /// construct from JSON
  SCFResults(const nlohmann::json &j) { from_json(j); }

  std::string key() const override { return model; }

  nlohmann::json to_json() const override {
    nlohmann::json j;
    j["scf_eigenvalues_a"] = tensor_to_json(aeps);
    j["scf_fock_a"] = tensor_to_json(afock);
    if (beps.size() > 0)
      j["scf_eigenvalues_b"] = tensor_to_json(beps);
    if (bfock.size() > 0)
      j["scf_fock_b"] = tensor_to_json(bfock);
    j["model"] = model;
    j["scf_total_energy"] = scf_total_energy;
    j["properties"] = properties.to_json();
    return j;
  } // from json SCFResults
  void from_json(const nlohmann::json &j) {

    if (j.count("scf_eigenvalues_a") > 0)
      aeps = tensor_from_json<double>(j["scf_eigenvalues_a"]);
    if (j.contains("scf_eigenvalues_b"))
      beps = tensor_from_json<double>(j["scf_eigenvalues_b"]);
    if (j.contains("scf_fock_a"))
      afock = tensor_from_json<double>(j["scf_fock_a"]);
    if (j.contains("scf_fock_b"))
      bfock = tensor_from_json<double>(j["scf_fock_b"]);
    if (j.contains("scf_total_energy"))
      scf_total_energy = j["scf_total_energy"];

    if (j.contains("properties")) {
      properties = PropertyResults(j["properties"]);
    } else {
      properties = PropertyResults();
    }
  };
};

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
        ei.oscillator_strength_length = ex.value("oscillator_strength_length", 0.0);
        ei.oscillator_strength_velocity = ex.value("oscillator_strength_velocity", 0.0);
        excitations.push_back(ei);
      }
    }
    nfreeze = j.value("nfreeze", -1);
    model = j.value("model", "unknown");
  }

  /// constructor with nfreeze and model
  CISResults(long nfreeze, const std::string &model) : nfreeze(nfreeze), model(model) {}

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

  /// construct from JSON
  CC2Results(const nlohmann::json &j) : CISResults(j) {
    properties = PropertyResults(j.value("properties", nlohmann::json{}));
    model = j.value("model", "mp2");
    correlation_energy = j.value("correlation_energy", 0.0);
    total_energy = j.value(model + "_total_energy", 0.0);
  }

  /// constructor with nfreeze and model
  CC2Results(long nfreeze, const std::string &model) : CISResults(nfreeze, model) {}

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

  PropertyResults properties;      // properties of the correlated calculation
  double correlation_energy = 0.0; // correlation energy of the correlated calculation
  double total_energy = 0.0;       // total energy of the correlated calculation
};

class ZnemoResults : public SCFResults {
public:
  double B = 0.0; // B value for the Znemo calculation

  ZnemoResults() = default;
  /// construct from JSON
  ZnemoResults(const nlohmann::json &j) : SCFResults(j) { B = j.value("B", 0.0); }

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

class ResponseResults : public ResultsBase {
  virtual nlohmann::json to_json() const { MADNESS_EXCEPTION("to_json not implemented for ConvergenceResults", 1); }
};

using SCFResultsTuple = std::tuple<SCFResults, PropertyResults, ConvergenceResults>;

} // namespace madness
#endif // RESULTS_H
