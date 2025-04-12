#pragma once
#include "ResponseState.hpp"
#include <filesystem>
#include <fstream>
#include <madness/external/nlohmann_json/json.hpp>
#include <string>

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
      data_["states"] = json::object();
    }
  }

  bool is_converged(const std::string &state_id, double frequency,
                    double threshold) const {
    std::string t_str = std::to_string(threshold);
    std::string f_str = std::to_string(frequency);
    return data_.contains(state_id) &&
           data_[state_id]["converged"].contains(t_str) &&
           data_[state_id]["converged"][t_str].contains(f_str) &&
           data_[state_id]["converged"][t_str][f_str].get<bool>();
  }

  // Check if state at given protocol is saved
  bool is_saved(const std::string &state_id, double protocol) const {
    auto protocol_str = protocol_to_string(protocol);
    return data_["states"].contains(state_id) &&
           data_["states"][state_id]["protocols"].contains(protocol_str) &&
           data_["states"][state_id]["protocols"][protocol_str].get<bool>();
  }

  // Mark state as saved at protocol (not necessarily converged)
  void mark_saved(const std::string &state_id, double protocol) {
    auto protocol_str = protocol_to_string(protocol);
    data_["states"][state_id]["protocols"][protocol_str] = true;

    double highest =
        data_["states"][state_id].value("highest_protocol_reached", 0.0);
    if (protocol > highest) {
      data_["states"][state_id]["highest_protocol_reached"] = protocol;
    }
    write();
  }

  void mark_final_converged(const std::string &state_id, bool flag = true) {
    data_[state_id]["final_converged"] = flag;
    write();
  }

  bool final_converged(const std::string &state_id) const {
    return data_.contains(state_id) &&
           data_[state_id].contains("final_converged") &&
           data_[state_id]["final_converged"].get<bool>();
  }

  void mark_converged(const std::string &state_id, double frequency,
                      double threshold) {
    data_[state_id]["converged"][std::to_string(threshold)]
         [std::to_string(frequency)] = true;
    write();
  }

  // Quickly determine the highest protocol computed
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
    fs::create_directories(fs::path(metadata_path_).parent_path());
    std::ofstream out(metadata_path_);
    out << std::setw(2) << data_ << "\n";
  }

  static std::string protocol_to_string(double protocol) {
    std::ostringstream ss;
    ss << std::scientific << protocol;
    return ss.str();
  }
};
