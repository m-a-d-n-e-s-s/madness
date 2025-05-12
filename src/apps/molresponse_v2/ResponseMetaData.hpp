#pragma once
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <madness/external/nlohmann_json/json.hpp>
#include <map>
#include <sstream>
#include <string>

#include "ResponseState.hpp"

using json = nlohmann::json;
namespace fs = std::filesystem;

class ResponseMetadata {
 public:
  explicit ResponseMetadata(World &world, const std::string &filepath)
      : path_(filepath) {
    if (fs::exists(path_)) {
      std::string json_string;
      std::ifstream in(path_);
      if (in.is_open()) {
        std::stringstream buffer;
        buffer << in.rdbuf();
        json_string = buffer.str();
        in.close();
      } else {
        std::cerr << "Error opening file: " << path_ << std::endl;
      }
      world.gop.fence();
      world.gop.broadcast_serializable(json_string, 0);

      data_ = json::parse(json_string);
      world.gop.fence();
    } else {
      data_["states"] = json::object();
    }
  }

  void initialize_states(const std::vector<LinearResponseDescriptor> &states) {
    for (const auto &state : states) {
      const std::string id = describe_perturbation(state.perturbation);

      const auto &frequencies = state.frequencies;
      const auto &protocols = state.thresholds;

      auto &state_entry = data_["states"][id];

      for (const auto &protocol : protocols) {
        std::string pstr = protocol_to_string(protocol);

        // Initialize protocol if missing
        if (!state_entry["protocols"].contains(pstr)) {
          for (const auto &freq : frequencies) {
            std::string fstr = frequency_to_string(freq);
            state_entry["protocols"][pstr]["saved"][fstr] = false;
            state_entry["protocols"][pstr]["converged"][fstr] = false;
          }
        }
      }

      // Final convergence flag (if not present)
      if (!state_entry.contains("final_converged")) {
        state_entry["final_converged"] = false;
      }
    }
    write();
  }

  void print_summary() const {
    for (const auto &[id, entry] : data_["states"].items()) {
      std::cout << "🧪 State: " << id << "\n";
      for (const auto &[protocol, proto_data] : entry["protocols"].items()) {
        std::cout << "  🔧 Protocol: " << protocol << "\n";
        for (const auto &[freq, saved] : proto_data["saved"].items()) {
          std::cout << "    freq=" << freq << " saved=" << saved
                    << " conv=" << proto_data["converged"][freq] << "\n";
        }
      }
      std::cout << "  ✅ Final converged: " << entry["final_converged"] << "\n";
    }
  }

  [[nodiscard]] bool is_saved(const std::string &state_id, double protocol,
                              double freq) const {
    return get_flag(state_id, protocol, freq, "saved");
  }

  [[nodiscard]] bool is_converged(const std::string &state_id, double protocol,
                                  double freq) const {
    return get_flag(state_id, protocol, freq, "converged");
  }

  void mark_saved(const std::string &state_id, double protocol, double freq) {
    set_flag(state_id, protocol, freq, "saved", true);
    write();
  }

  void mark_converged(const std::string &state_id, double protocol, double freq,
                      bool converged) {
    set_flag(state_id, protocol, freq, "converged", converged);
    write();
  }

  void mark_final_converged(const std::string &state_id, bool flag = true) {
    data_["states"][state_id]["final_converged"] = flag;
    write();
  }

  void mark_final_saved(const std::string &state_id, bool flag = true) {
    data_["states"][state_id]["final_saved"] = flag;
    write();
  }

  [[nodiscard]] bool final_converged(const std::string &state_id) const {
    return data_["states"][state_id].value("final_converged", false);
  }

  [[nodiscard]] bool final_saved(const std::string &state_id) const {
    return data_["states"][state_id].value("final_saved", false);
  }

  void write() const {
    std::ofstream out(path_);
    out << std::setw(2) << data_ << "\n";
  }
  json to_json() const { return data_; }

 private:
  std::string path_;
  json data_;

  static std::string to_string(double val) {
    std::ostringstream oss;
    oss << std::scientific << val;
    return oss.str();
  }

  [[nodiscard]] bool get_flag(const std::string &state_id, double protocol,
                              double freq, const std::string &key) const {
    std::string p_str = to_string(protocol);
    std::string f_str = to_string(freq);

    if (!data_["states"].contains(state_id)) return false;
    if (!data_["states"][state_id]["protocols"].contains(p_str)) return false;
    if (!data_["states"][state_id]["protocols"][p_str][key].contains(f_str))
      return false;

    return data_["states"][state_id]["protocols"][p_str][key][f_str]
        .get<bool>();
  }

  void set_flag(const std::string &state_id, double protocol, double freq,
                const std::string &key, bool value) {
    std::string p_str = to_string(protocol);
    std::string f_str = to_string(freq);

    data_["states"][state_id]["protocols"][p_str][key][f_str] = value;
    write();
  }
  static std::string frequency_to_string(double frequency) {
    std::ostringstream ss;
    ss << std::scientific << frequency;
    return ss.str();
  }

  static std::string protocol_to_string(double protocol) {
    std::ostringstream ss;
    ss << std::scientific << protocol;
    return ss.str();
  }
};
