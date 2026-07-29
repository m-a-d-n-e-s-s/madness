#pragma once
#include <filesystem>

namespace madness {
  struct PathManager {
    std::filesystem::path baseDir;
    std::string label;

    // Resolve baseDir to absolute at construction (cwd is the parent here — every
    // Application builds its PathManager before its ScopedCWD), so dir()/file()
    // and any work_dir set from them are cwd-independent. This is what lets a
    // reference work_dir survive a later chdir into another task's run dir
    // (the v3 madqc adapter + CC2/TDHF/OEP cross-task archive lookups).
    PathManager(std::filesystem::path base, std::string lbl)
        : baseDir(std::filesystem::absolute(std::move(base))), label(std::move(lbl)) {}

    std::filesystem::path dir() const { return baseDir / label; }

    void create() const {
      std::error_code ec;
      std::filesystem::create_directories(dir(), ec);
      if (ec) throw std::runtime_error("Could not create " + dir().string());
    }

    std::filesystem::path file(std::string const& name) const { return dir() / name; }
  };
}