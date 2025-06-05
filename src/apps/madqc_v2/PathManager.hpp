#pragma once
#include <filesystem>

struct PathManager {
  std::filesystem::path baseDir;
  std::string label;

  PathManager(std::filesystem::path base, std::string lbl) : baseDir(std::move(base)), label(std::move(lbl)) {}

  std::filesystem::path dir() const { return baseDir / label; }

  void create() const {
    std::error_code ec;
    std::filesystem::create_directories(dir(), ec);
    if (ec) throw std::runtime_error("Could not create " + dir().string());
  }

  std::filesystem::path file(std::string const& name) const { return dir() / name; }
};
