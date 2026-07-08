#pragma once

// ResultsSummary.hpp
// Renders a human-readable chemistry report from the aggregated madqc
// `calc_info.json` (the machine-readable source of truth produced by
// qcapp::Workflow). This is presentation only: it never recomputes anything
// and tolerates missing fields, so it works for any subset of workflows that
// happen to be present in the JSON.
//
// The JSON shape it consumes is `{"tasks": [ <task>, ... ]}`, where each task
// is one Application's results() fragment:
//   * SCF      : { "model":"scf", "properties":{energy,dipole,gradient},
//                  "scf_eigenvalues_a"/"_b" (tensor-json), "convergence_info",
//                  "metadata", "molecule" }
//   * CIS/TDHF : { "model":"cis", "excitations":[{omega,irrep,
//                  oscillator_strength_length,...}], "nfreeze" }
//   * MP2/CC2  : { "model":"mp2"|"cc2", "correlation_energy",
//                  "<model>_total_energy", "excitations"? }
//   * Response : { "type":"response", "properties":{response_properties:[...],
//                  raman_spectra, vibrational_analysis} }

#include <madness/constants.h>
#include <nlohmann/json.hpp>

#include <algorithm>
#include <cmath>
#include <iomanip>
#include <map>
#include <ostream>
#include <sstream>
#include <string>
#include <vector>

namespace qcapp {
namespace summary_detail {

// Atomic-unit conversions.
inline constexpr double kHa2eV =
    madness::constants::hartree_electron_volt_relationship; // 27.2113862...
inline constexpr double kHa2cm = madness::constants::au2invcm; // 219474.63...
inline constexpr double kAu2Debye = 2.541746473; // 1 a.u. dipole = 2.5417 D

// Flat value array of a MADNESS tensor-json ({"vals":[...],"size":N,...}).
// Returns empty for absent / zero-size tensors.
inline std::vector<double> tensor_vals(const nlohmann::json &j) {
  std::vector<double> v;
  if (j.is_object() && j.value("size", 0) > 0 && j.contains("vals") &&
      j["vals"].is_array())
    v = j["vals"].get<std::vector<double>>();
  return v;
}

// Empirical formula (Hill notation) from a molecule-json "symbols" list:
// carbon first, hydrogen second, then all other elements alphabetically; if
// no carbon is present, every element (including H) is alphabetical.
inline std::string molecule_formula(const nlohmann::json &molj) {
  if (!molj.is_object() || !molj.contains("symbols"))
    return {};
  std::map<std::string, int> counts; // std::map iterates keys alphabetically
  for (const auto &s : molj["symbols"])
    counts[s.get<std::string>()]++;

  std::vector<std::string> order;
  const bool hasC = counts.count("C") > 0;
  if (hasC) {
    order.push_back("C");
    if (counts.count("H"))
      order.push_back("H");
  }
  for (const auto &kv : counts) {
    if (kv.first == "C" || (hasC && kv.first == "H"))
      continue;
    order.push_back(kv.first);
  }

  std::ostringstream os;
  for (const auto &sym : order) {
    os << sym;
    if (counts[sym] > 1)
      os << counts[sym];
  }
  return os.str();
}

// "Dipole_z" -> "z"; otherwise pass through unchanged.
inline std::string clean_axis(const std::string &s) {
  const std::string pfx = "Dipole_";
  if (s.rfind(pfx, 0) == 0)
    return s.substr(pfx.size());
  return s;
}

inline std::string rule(char c = '-', int n = 70) { return std::string(n, c); }

inline void write_scf_section(std::ostream &os, const nlohmann::json &t) {
  const auto &props = t.value("properties", nlohmann::json::object());

  if (t.contains("molecule")) {
    auto formula = molecule_formula(t["molecule"]);
    if (!formula.empty())
      os << "    Molecule         : " << formula << "\n";
  }

  // Energy: moldft stores the SCF energy under properties.energy and leaves
  // the top-level scf_total_energy at 0; prefer the former.
  double e = props.value("energy", t.value("scf_total_energy", 0.0));
  if (e != 0.0)
    os << "    Total energy     : " << std::fixed << std::setprecision(9)
       << std::setw(18) << e << " Ha   " << std::setprecision(4)
       << std::setw(14) << e * kHa2eV << " eV\n";

  // Dipole moment (a.u. + Debye).
  if (props.contains("dipole")) {
    auto d = tensor_vals(props["dipole"]);
    if (d.size() == 3) {
      double mag = std::sqrt(d[0] * d[0] + d[1] * d[1] + d[2] * d[2]);
      os << "    Dipole (a.u.)    :  x " << std::showpos << std::fixed
         << std::setprecision(6) << d[0] << "  y " << d[1] << "  z " << d[2]
         << std::noshowpos << "\n";
      os << "                     : |mu| = " << std::setprecision(6) << mag
         << " a.u. = " << mag * kAu2Debye << " Debye\n";
    }
  }

  // Nuclear gradient (forces): report max and RMS.
  if (props.contains("gradient")) {
    auto g = tensor_vals(props["gradient"]);
    if (!g.empty()) {
      double mx = 0.0, sq = 0.0;
      for (double v : g) {
        mx = std::max(mx, std::fabs(v));
        sq += v * v;
      }
      os << "    Nuclear gradient :  max |g| = " << std::scientific
         << std::setprecision(2) << mx
         << "   rms = " << std::sqrt(sq / g.size()) << " Ha/bohr\n"
         << std::defaultfloat;
    }
  }

  // Orbital energies.
  auto print_eps = [&](const char *label, const char *key) {
    if (!t.contains(key))
      return;
    auto eps = tensor_vals(t[key]);
    if (eps.empty())
      return;
    os << "    " << std::left << std::setw(17) << label << std::right << ":";
    os << std::fixed << std::setprecision(4);
    for (double v : eps)
      os << " " << v;
    os << "  (Ha)\n" << std::defaultfloat;
  };
  print_eps("Orbital eps (a)", "scf_eigenvalues_a");
  print_eps("Orbital eps (b)", "scf_eigenvalues_b");

  // Convergence + metadata.
  if (t.contains("convergence_info")) {
    const auto &c = t["convergence_info"];
    os << "    Converged        :  thresh = " << c.value("converged_for_thresh", 0.0)
       << "  dconv = " << c.value("converged_for_dconv", 0.0) << "\n";
  }
  if (t.contains("metadata")) {
    const auto &m = t["metadata"];
    os << "    Wall time        :  " << std::fixed << std::setprecision(1)
       << m.value("elapsed_time", 0.0) << " s   (" << m.value("mpi_size", 1)
       << " MPI x " << m.value("nthreads", 1) << " threads)\n"
       << std::defaultfloat;
  }
}

inline void write_excitations_section(std::ostream &os, const nlohmann::json &t) {
  if (!t.contains("excitations") || !t["excitations"].is_array() ||
      t["excitations"].empty())
    return;
  os << "    " << "Root  Irrep        E (Ha)      E (eV)   lambda(nm)"
     << "    f(length)\n";
  int root = 1;
  for (const auto &ex : t["excitations"]) {
    double w = ex.value("omega", 0.0);
    double nm = (w > 0.0) ? 1.0e7 / (w * kHa2cm) : 0.0;
    os << "    " << std::setw(4) << root++ << "  " << std::setw(5)
       << ex.value("irrep", std::string("?")) << "  " << std::fixed
       << std::setprecision(6) << std::setw(12) << w << std::setprecision(4)
       << std::setw(12) << w * kHa2eV << std::setprecision(2) << std::setw(12)
       << nm << std::setprecision(6) << std::setw(13)
       << ex.value("oscillator_strength_length", 0.0) << "\n"
       << std::defaultfloat;
  }
}

inline void write_correlation_section(std::ostream &os, const nlohmann::json &t) {
  std::string model = t.value("model", std::string("correlated"));
  if (t.contains("correlation_energy"))
    os << "    Correlation energy : " << std::fixed << std::setprecision(9)
       << std::setw(18) << t.value("correlation_energy", 0.0) << " Ha\n";
  std::string tot = model + "_total_energy";
  if (t.contains(tot))
    os << "    Total energy (" << model << ") : " << std::fixed
       << std::setprecision(9) << std::setw(18) << t.value(tot, 0.0) << " Ha\n"
       << std::defaultfloat;
}

inline void write_response_section(std::ostream &os, const nlohmann::json &t) {
  const auto &props = t.value("properties", nlohmann::json::object());
  // v3 shape: response_properties is an OBJECT keyed property -> protocol_key
  // -> row(s). (The old v2 engine emitted a flat array — handled below for
  // backward compatibility with archived calc_info files.)
  if (props.contains("response_properties") &&
      props["response_properties"].is_object()) {
    for (const auto &[prop, by_key] : props["response_properties"].items()) {
      if (!by_key.is_object()) continue;
      for (const auto &[pkey, rows_raw] : by_key.items()) {
        const auto rows = rows_raw.is_array()
                              ? rows_raw
                              : nlohmann::json::array({rows_raw});
        for (const auto &r : rows) {
          if (prop == "alpha" && r.contains("alpha")) {
            const std::string dirs = r.value("directions", std::string("?"));
            os << "    alpha[" << pkey << "](w=" << std::fixed
               << std::setprecision(3) << r.value("omega", 0.0) << ")"
               << (r.value("converged", false) ? "" : "  [NOT converged]")
               << "\n";
            const auto &m = r["alpha"];
            for (size_t i = 0; i < m.size(); ++i) {
              os << "      " << (i < dirs.size() ? dirs[i] : '?') << " :";
              for (const auto &v : m[i])
                os << std::setprecision(6) << std::setw(14) << v.get<double>();
              os << "\n";
            }
          } else if (r.contains("beta")) {   // beta and raman rows
            os << "    " << prop << "[" << pkey << "]("
               << r.value("A", std::string("?")) << ","
               << r.value("B", std::string("?")) << ","
               << r.value("C", std::string("?")) << "; wB=" << std::fixed
               << std::setprecision(3) << r.value("freq_b", 0.0)
               << ",wC=" << r.value("freq_c", 0.0) << ") = "
               << std::setprecision(6) << std::setw(14)
               << r.value("beta", 0.0) << "\n";
          }
        }
      }
    }
    os << std::defaultfloat;
  }
  if (props.contains("response_properties") &&
      props["response_properties"].is_array()) {
    for (const auto &e : props["response_properties"]) {
      std::string prop = e.value("property", std::string());
      std::string comp;
      if (e.contains("component") && e["component"].is_array())
        for (const auto &c : e["component"])
          comp += clean_axis(c.get<std::string>());
      double val = e.value("value", 0.0);
      if (prop == "polarizability") {
        os << "    alpha_" << comp << "(w=" << std::fixed << std::setprecision(3)
           << e.value("freqB", 0.0) << ")  = " << std::setprecision(6)
           << std::setw(14) << val << "\n";
      } else if (prop == "hyperpolarizability") {
        os << "    beta_" << comp << "(wB=" << std::fixed << std::setprecision(3)
           << e.value("freqB", 0.0) << ",wC=" << e.value("freqC", 0.0)
           << ") = " << std::setprecision(6) << std::setw(14) << val << "\n";
      } else if (!prop.empty()) {
        os << "    " << prop << " " << comp << " = " << std::fixed
           << std::setprecision(6) << val << "\n";
      }
    }
    os << std::defaultfloat;
  }

  // Vibrational frequencies, if a Raman/vibrational analysis is present.
  if (props.contains("vibrational_analysis") &&
      props["vibrational_analysis"].contains("frequencies")) {
    auto f = tensor_vals(props["vibrational_analysis"]["frequencies"]);
    if (!f.empty()) {
      os << "    Vibrational frequencies (cm^-1):";
      for (double w : f)
        os << " " << std::fixed << std::setprecision(1) << w * kHa2cm;
      os << "\n" << std::defaultfloat;
    }
  }
}

} // namespace summary_detail

/// Write a human-readable chemistry report for an aggregated calc_info JSON.
inline void write_results_summary(std::ostream &os,
                                  const nlohmann::json &calc_info) {
  using namespace summary_detail;
  os << "\n" << rule('=') << "\n";
  os << "  MADQC RESULTS SUMMARY\n";
  os << rule('=') << "\n";

  if (!calc_info.contains("tasks") || !calc_info["tasks"].is_array()) {
    os << "  (no tasks found)\n";
    return;
  }

  int i = 0;
  for (const auto &t : calc_info["tasks"]) {
    const std::string type = t.value("type", std::string());
    const std::string model = t.value("model", std::string());

    std::string title;
    if (type == "response")
      title = "RESPONSE";
    else if (model == "scf" || t.contains("scf_eigenvalues_a"))
      title = "SCF  (model = " + (model.empty() ? "scf" : model) + ")";
    else if (t.contains("excitations") && !t.contains("correlation_energy"))
      title = "EXCITED STATES  (model = " + model + ")";
    else if (t.contains("correlation_energy"))
      title = "CORRELATION  (model = " + model + ")";
    else
      title = "TASK";

    os << "\n  Task " << i++ << " : " << title << "\n";
    os << "  " << rule('-') << "\n";

    if (t.value("status", std::string()) == "failed") {
      os << "    *** FAILED: " << t.value("error", std::string("(no detail)"))
         << "\n";
      continue;
    }

    if (type == "response") {
      write_response_section(os, t);
    } else if (model == "scf" || t.contains("scf_eigenvalues_a")) {
      write_scf_section(os, t);
    } else {
      if (t.contains("correlation_energy"))
        write_correlation_section(os, t);
      if (t.contains("excitations"))
        write_excitations_section(os, t);
    }
  }
  os << "\n" << rule('=') << "\n\n";
}

} // namespace qcapp
