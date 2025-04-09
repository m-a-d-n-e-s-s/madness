#pragma once
#include "ResponseState.hpp"
#include <fstream>
#include <iomanip>
#include <madness/external/nlohmann_json/json.hpp>

using json = nlohmann::json;

class ResponseDebugLogger {
public:
  explicit ResponseDebugLogger(bool enabled = false) : enabled_(enabled) {}

  [[nodiscard]] bool enabled() const { return enabled_; }
  void set_enabled(bool on) { enabled_ = on; }

  void start_state(const ResponseState &state) {
    std::string key = state.description();
    current_entry_ = json{{"perturbation", state.perturbationDescription()},
                          {"frequency", state.current_frequency()},
                          {"threshold", state.current_threshold()},
                          {"iterations", json::array()}};
    current_key_ = key;
  }

  // Start a new iteration
  void begin_iteration(size_t iter_index) {
    current_iter_ = {{"iter", iter_index}, {"steps", json::object()}};
  }
  // Log any key-value pairs under a specific step within the iteration
  template <typename T>
  void log_value(const std::string &step_name, const T &value) {
    current_iter_["steps"][step_name] = value;
  }

  // Log a timing in seconds
  void log_timing(const std::string &step_name, double wall_time,
                  double cpu_time) {
    current_iter_["steps"][step_name]["wall_time"] = wall_time;
    current_iter_["steps"][step_name]["cpu_time"] = cpu_time;
  }
  template <typename T>
  void log_value_and_time(const std::string &key, const T &value,
                          double wall_time, double cpu_time) {
    current_iter_["steps"][key]["value"] = value;
    current_iter_["steps"][key]["wall_time"] = wall_time;
    current_iter_["steps"][key]["cpu_time"] = cpu_time;
  }

  // Commit current iteration
  void end_iteration() {
    current_entry_["iterations"].push_back(current_iter_);
  }

  // Finalize the state and store it in the log
  void finalize_state(const ResponseState &state) {
    log_data_[current_key_] = current_entry_;
  }

  void write_to_disk(const std::string &filename) const {
    std::ofstream out(filename);
    out << std::setw(2) << log_data_ << "\n";
  }

  // Optional pretty printer to check high-level results quickly
  void pretty_print_summary(const std::string &description) const {
    if (!log_data_.contains(description))
      return;

    const auto &iter_data = log_data_.at(description)["iterations"];
    std::cout << "\n📋 Summary for: " << description << "\n";
    std::cout << std::setw(6) << "Iter" << std::setw(15) << "drho"
              << std::setw(15) << "residual" << std::setw(15) << "maxrotn"
              << std::setw(10) << "time" << "\n";

    for (const auto &iter : iter_data) {
      const auto &steps = iter["steps"];
      std::cout << std::setw(6) << iter["iter"] << std::setw(15)
                << steps["convergence"]["drho"] << std::setw(15)
                << steps["convergence"]["residual"] << std::setw(15)
                << steps["convergence"]["maxrotn"] << std::setw(10)
                << steps["convergence"]["time"] << "\n";
    }
  }

private:
  json log_data_;
  json current_entry_;
  json current_iter_;
  std::string current_key_;
  bool enabled_ = false;
};

class TimedValueLogger {
public:
  TimedValueLogger(madness::World &world, const std::string &key,
                   ResponseDebugLogger *logger = nullptr)
      : world_(world), key_(key), logger_(logger) {
    world_.gop.fence();
    start_wall_ = madness::wall_time();
    start_cpu_ = madness::cpu_time();
  }

  void log() {
    double wall = madness::wall_time() - start_wall_;
    double cpu = madness::cpu_time() - start_cpu_;

    if (logger_ && world_.rank() == 0) {
      logger_->log_timing(key_, wall, cpu);
    }

    if (world_.rank() == 0) {
      std::cout << std::left << std::setw(30) << "⏱️ [" + key_ + "]"
                << std::right << " | Wall: " << std::setw(7)
                << std::setprecision(3) << wall << "s | CPU: " << std::setw(7)
                << std::setprecision(3) << cpu << "s |" << std::endl;
    }
  }

  template <typename T> void log(const T &value) {
    double wall = madness::wall_time() - start_wall_;
    double cpu = madness::cpu_time() - start_cpu_;

    if (logger_ && world_.rank() == 0) {
      logger_->log_value_and_time(key_, value, wall, cpu);
    }

    if (world_.rank() == 0) {
      std::cout << std::left << std::setw(30) << "⏱️ [" + key_ + "]"
                << std::right << " | Wall: " << std::setw(7)
                << std::setprecision(3) << wall << "s | CPU: " << std::setw(7)
                << std::setprecision(3) << cpu << "s | Value: " << std::setw(10)
                << std::setprecision(7) << std::fixed << value << " |"
                << std::endl;
    }
  }
  template <typename T> void log(const std::vector<T> &values) {
    double wall = madness::wall_time() - start_wall_;
    double cpu = madness::cpu_time() - start_cpu_;

    if (logger_ && world_.rank() == 0) {
      for (size_t i = 0; i < values.size(); ++i) {
        std::string key_i = key_ + "[" + std::to_string(i) + "]";
        logger_->log_value_and_time(key_i, values[i], wall, cpu);
      }
    }

    if (world_.rank() == 0) {
      std::cout << std::left << std::setw(20) << ("⏱️ [" + key_ + "]")
                << " | Values: ";
      for (const auto &v : values) {
        std::cout << std::right << std::setw(10) << std::setprecision(6)
                  << std::fixed << v << " ";
      }
      std::cout << "| Wall: " << std::setw(7) << std::setprecision(3) << wall
                << "s | CPU: " << std::setw(7) << std::setprecision(3) << cpu
                << "s |\n";
    }
  }

private:
  madness::World &world_;
  std::string key_;
  ResponseDebugLogger *logger_;
  double start_wall_, start_cpu_;
};
