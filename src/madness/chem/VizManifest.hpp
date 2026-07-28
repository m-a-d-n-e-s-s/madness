#pragma once

// VizManifest.hpp
// After a workflow runs, scan its output tree for visualization artifacts
// (density / orbital / potential plots written by SCF::do_plots when the
// plotdens/plotcube/... knobs are set) and emit <prefix>.viz_manifest.json
// indexing them. This is the driver-side "viz integration": madqc organizes
// per-task plot files and hands downstream tools (gecko, ParaView, VMD) a
// single discovery index tied to the molecule and result files.
//
// Pure post-run filesystem scan — it neither produces nor modifies any field
// data, so it adds no cost to the compute path and is safe to call after run().

#include <madness/chem/ResultsSummary.hpp> // qcapp::summary_detail::molecule_formula
#include <madness/external/nlohmann_json/json.hpp>

#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <string>

namespace qcapp {
namespace viz_detail {

struct Classified {
  std::string kind;
  int index = -1; // orbital index, when applicable
};

// Classify a plot file by its stem (filename without extension).
inline Classified classify(const std::string &stem) {
  auto starts = [&](const char *p) { return stem.rfind(p, 0) == 0; };
  if (stem.size() >= 7 && stem.compare(stem.size() - 7, 7, "density") == 0)
    return {"density", -1}; // total_density, spin_density
  if (starts("coulomb") || starts("vlocal") || starts("vcoul") ||
      starts("vlda"))
    return {"potential", -1};
  if (starts("amo-") || starts("bmo-")) {
    int idx = -1;
    try {
      idx = std::stoi(stem.substr(4));
    } catch (...) {
    }
    return {std::string("orbital_") + (starts("amo-") ? "a" : "b"), idx};
  }
  return {"field", -1};
}

} // namespace viz_detail

/// Scan the <prefix>/ output tree and write <prefix>.viz_manifest.json. Returns
/// the number of artifacts found (0 -> no manifest written).
inline std::size_t write_viz_manifest(const std::string &prefix,
                                      const nlohmann::json &calc_info) {
  namespace fs = std::filesystem;
  nlohmann::json artifacts = nlohmann::json::array();
  const fs::path topDir = prefix;

  if (fs::exists(topDir) && fs::is_directory(topDir)) {
    try {
      for (const auto &e : fs::recursive_directory_iterator(topDir)) {
        if (!e.is_regular_file())
          continue;
        const auto &p = e.path();
        const std::string ext = p.extension().string();
        std::string format;
        if (ext == ".cube")
          format = "cube";
        else if (ext == ".dx")
          format = "dx";
        else if (ext == ".dat")
          format = "line";
        else
          continue; // not a recognized visualization artifact

        // Derive task index + method label from the relative path
        // (<prefix>/task_<N>/<method>/<file>).
        int task = -1;
        std::string method;
        bool saw_task = false;
        for (const auto &comp : fs::relative(p, topDir)) {
          const std::string c = comp.string();
          if (c.rfind("task_", 0) == 0) {
            try {
              task = std::stoi(c.substr(5));
            } catch (...) {
            }
            saw_task = true;
          } else if (saw_task && method.empty() &&
                     c.find('.') == std::string::npos) {
            method = c;
          }
        }

        const auto cl = viz_detail::classify(p.stem().string());
        nlohmann::json a;
        a["path"] = p.string();
        a["format"] = format;
        a["kind"] = cl.kind;
        if (cl.index >= 0)
          a["index"] = cl.index;
        if (task >= 0)
          a["task"] = task;
        if (!method.empty())
          a["method"] = method;
        artifacts.push_back(std::move(a));
      }
    } catch (const std::exception &ex) {
      std::cout << "viz manifest: directory scan stopped early: " << ex.what()
                << "\n";
    }
  }

  if (artifacts.empty()) {
    std::cout << "No visualization artifacts found "
                 "(enable with plotdens / plotcube in the 'dft' group).\n";
    return 0;
  }

  nlohmann::json m;
  m["prefix"] = prefix;
  m["results"] = {{"calc_info", prefix + ".calc_info.json"},
                  {"summary", prefix + ".out"}};
  if (calc_info.contains("tasks") && calc_info["tasks"].is_array() &&
      !calc_info["tasks"].empty() && calc_info["tasks"][0].contains("molecule")) {
    const auto &molj = calc_info["tasks"][0]["molecule"];
    nlohmann::json mol;
    mol["formula"] = summary_detail::molecule_formula(molj);
    if (molj.contains("symbols"))
      mol["symbols"] = molj["symbols"];
    if (molj.contains("geometry"))
      mol["geometry"] = molj["geometry"];
    m["molecule"] = mol;
  }
  m["artifacts"] = artifacts;

  std::ofstream ofs(prefix + ".viz_manifest.json");
  ofs << std::setw(2) << m << std::endl;
  std::cout << "Wrote viz manifest    : " << prefix << ".viz_manifest.json ("
            << artifacts.size() << " artifacts)\n";
  return artifacts.size();
}

} // namespace qcapp
