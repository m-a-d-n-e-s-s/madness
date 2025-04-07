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
  explicit ResponseMetadata(const std::string &perturbation_id)
      : id(perturbation_id),
        metadata_path("responses/metadata/" + id + ".metadata.json") {
    if (fs::exists(metadata_path)) {
      std::ifstream in(metadata_path);
      bool is_json = json::accept(in);
      in.close();

      if (is_json) {
        in.open(metadata_path);
        data = json::parse(in);
        in.close();
      } else {
        std::cerr << "Warning: Metadata file is not valid JSON. "
                     "Creating a new metadata file."
                  << std::endl;
        data = json::object();
      }
    }
  }

  void mark_converged(double frequency, double threshold) {
    data[std::to_string(threshold)][std::to_string(frequency)] = true;
    write();
  }

  bool is_converged(double frequency, double threshold) const {
    auto t_str = std::to_string(threshold);
    auto f_str = std::to_string(frequency);
    return data.contains(t_str) && data[t_str].contains(f_str) &&
           data[t_str][f_str].get<bool>() == true;
  }

public:
  void load() {
    if (fs::exists(metadata_path)) {

      std::ifstream in(metadata_path);
      in.open(metadata_path);
      data = json::parse(in);
      in.close();
    }
  }

  void save() const { write(); }

private:
  std::string id;
  std::string metadata_path;
  json data;

  void write() const {
    fs::create_directories("responses/metadata");
    std::ofstream out(metadata_path);
    out << std::setw(4) << data;
  }
};
