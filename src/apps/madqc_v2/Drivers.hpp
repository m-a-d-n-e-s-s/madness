#pragma once
#include <filesystem>
#include <madness/external/nlohmann_json/json.hpp>

#include "Applications.hpp"

class Driver {
 public:
  virtual ~Driver() = default;
  virtual void execute(const std::filesystem::path& workdir) = 0;
  virtual nlohmann::json summary() const = 0;
};

class SinglePointDriver : public Driver {
 public:
  SinglePointDriver(std::unique_ptr<Application> app) : app_(std::move(app)) {}

  void execute(const std::filesystem::path& workdir) override {
    app_->run(workdir);
    collected_ = app_->results();
  }

  nlohmann::json summary() const override { return collected_; }

 private:
  std::unique_ptr<Application> app_;
  nlohmann::json collected_;
};
