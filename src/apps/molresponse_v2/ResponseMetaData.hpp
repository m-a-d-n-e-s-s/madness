#pragma once
#include "ResponseState.hpp"
#include <filesystem>
#include <fstream>
#include <madness/external/nlohmann_json/json.hpp>
#include <string>

using json = nlohmann::json;
namespace fs = std::filesystem;

using json = nlohmann::json;
namespace fs = std::filesystem;

class ResponseMetadata {
public:
  explicit ResponseMetadata(const std::string &metadata_file)
      : metadata_path_(metadata_file) {
    if (fs::exists(metadata_path_)) {
      std::ifstream in(metadata_path_);
      in >> data_;
    } else {
      fs::create_directories(fs::path(metadata_path_).parent_path());
      data_["states"] = json::object();
    }
  }

  // ========= Convergence =========
  bool is_converged(const std::string &state_id, double frequency) const {
    auto freq_str = frequency_to_string(frequency);
    return data_["states"].contains(state_id) &&
           data_["states"][state_id]["converged"].contains(freq_str) &&
           data_["states"][state_id]["converged"][freq_str].get<bool>();
  }

  void mark_converged(const std::string &state_id, double frequency) {
    auto freq_str = frequency_to_string(frequency);
    data_["states"][state_id]["converged"][freq_str] = true;
    write();
  }

  // ========= Saved (intermediate or final) =========
  bool is_saved(const std::string &state_id, double frequency,
                double protocol) const {
    auto freq_str = frequency_to_string(frequency);
    auto protocol_str = protocol_to_string(protocol);
    return data_["states"].contains(state_id) &&
           data_["states"][state_id]["saved"].contains(protocol_str) &&
           data_["states"][state_id]["saved"][protocol_str].contains(
               freq_str) &&
           data_["states"][state_id]["saved"][protocol_str][freq_str]
               .get<bool>();
  }

  void mark_saved(const std::string &state_id, double frequency,
                  double protocol) {
    auto freq_str = frequency_to_string(frequency);
    auto protocol_str = protocol_to_string(protocol);
    data_["states"][state_id]["saved"][protocol_str][freq_str] = true;

    double highest =
        data_["states"][state_id].value("highest_protocol_reached", 0.0);
    if (protocol > highest) {
      data_["states"][state_id]["highest_protocol_reached"] = protocol;
    }

    write();
  }

  // ========= Final flag for entire state =========
  void mark_final_converged(const std::string &state_id, bool flag = true) {
    data_["states"][state_id]["final_converged"] = flag;
    write();
  }

  bool final_converged(const std::string &state_id) const {
    return data_["states"].contains(state_id) &&
           data_["states"][state_id].contains("final_converged") &&
           data_["states"][state_id]["final_converged"].get<bool>();
  }

  double highest_protocol_reached(const std::string &state_id) const {
    return data_["states"].contains(state_id)
               ? data_["states"][state_id].value("highest_protocol_reached",
                                                 0.0)
               : 0.0;
  }

private:
  json data_;
  std::string metadata_path_;

  void write() const {
    std::ofstream out(metadata_path_);
    out << std::setw(2) << data_ << "\n";
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
