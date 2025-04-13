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
    current_entry_ = json{
        {"perturbation", state.perturbationDescription()},
        {"frequency", state.current_frequency()},
        {"threshold", state.current_threshold()},
        {"iteration_values", json::array()},
        {"iteration_timings", json::array()},
    };

    current_key_ = key;
  }

  // Start a new iteration
  void begin_iteration(size_t iter_index) {
    current_iter_timing_ = {{"iter", iter_index}, {"steps", json::object()}};
    current_iter_values_ = {{"iter", iter_index}, {"steps", json::object()}};
  }
  // Log any key-value pairs under a specific step within the iteration
  template <typename T>
  void log_value(const std::string &step_name, const T &value) {
    current_iter_values_["steps"][step_name]["value"] = value;
  }

  // Log a timing in seconds
  void log_timing(const std::string &step_name, double wall_time,
                  double cpu_time) {

    current_iter_timing_["steps"][step_name]["wall_time"] = wall_time;
    current_iter_timing_["steps"][step_name]["cpu_time"] = cpu_time;
  }
  template <typename T>
  void log_value_and_time(const std::string &key, const T &value,
                          double wall_time, double cpu_time) {

    // Log the value
    log_value(key, value);
    log_timing(key, wall_time, cpu_time);
  }

  // Commit current iteration
  void end_iteration() {
    current_entry_["iteration_values"].push_back(current_iter_values_);
    current_entry_["iteration_timings"].push_back(current_iter_timing_);
  }

  // Finalize the state and store it in the log
  void finalize_state() { log_data_[current_key_] = current_entry_; }

  void write_to_disk(const std::string &filename) const {
    std::ofstream out(filename);
    out << std::setw(2) << log_data_ << "\n";
  }
  void write_to_disk() const {
    std::ofstream out("response_log.json");
    out << std::setw(2) << log_data_ << "\n";
  }

  void print_timing_table(const std::string &description) const {
    if (!log_data_.contains(description))
      return;

    const auto &iter_data = log_data_.at(description)["iteration_timings"];

    // Collect unique step names
    std::set<std::string> step_names;
    for (const auto &iter : iter_data) {
      for (const auto &step : iter["steps"].items()) {
        step_names.insert(step.key());
      }
    }

    // Map each step name to a 5-character short key
    std::map<std::string, std::string> short_keys;
    int index = 0;
    for (const auto &name : step_names) {
      std::ostringstream oss;
      oss << std::setw(5) << std::left << name.substr(0, 5);
      short_keys[name] = oss.str();
    }

    constexpr int col_width = 10;

    // Print header
    std::cout << "\n⏱️ Timing Table for: " << description << "\n";
    std::cout << std::left << std::setw(6) << "Iter";
    for (const auto &step : step_names) {
      std::cout << std::setw(col_width) << (short_keys[step]);
    }
    std::cout << "\n"
              << std::string(6 + step_names.size() * col_width, '-') << "\n";

    // Print each iteration row
    for (const auto &iter : iter_data) {
      std::cout << std::setw(6) << iter["iter"] << "     ";
      const auto &steps = iter["steps"];
      for (const auto &step : step_names) {
        if (steps.contains(step)) {
          std::cout << std::setw(col_width) << std::fixed
                    << std::setprecision(4)
                    << steps[step]["wall_time"].get<double>();
        } else {
          std::cout << std::setw(col_width) << "N/A" << std::setw(col_width)
                    << "N/A";
        }
      }
      std::cout << "\n";
    }
  }

  void print_values_table(const std::string &description) const {
    if (!log_data_.contains(description))
      return;

    const auto &iter_data = log_data_.at(description)["iteration_values"];

    // Collect all step names
    std::set<std::string> step_names;
    for (const auto &iter : iter_data) {
      for (const auto &step : iter["steps"].items()) {
        step_names.insert(step.key());
      }
    }

    constexpr int col_width = 16;

    std::cout << "\n📋 Value Table for: " << description << "\n";
    std::cout << std::setw(6) << "Iter";
    for (const auto &step_name : step_names) {
      std::cout << std::setw(col_width) << step_name;
    }
    std::cout << "\n"
              << std::string(6 + col_width * step_names.size(), '-') << "\n";

    for (const auto &iter : iter_data) {
      std::cout << std::setw(6) << iter["iter"] << "     ";
      const auto &steps = iter["steps"];
      for (const auto &step_name : step_names) {
        if (steps.contains(step_name)) {
          std::cout << std::setw(col_width) << std::scientific
                    << std::setprecision(5)
                    << static_cast<double>(steps[step_name]["value"]);
        } else {
          std::cout << std::setw(col_width) << "N/A";
        }
      }
      std::cout << "\n";
    }
  }

private:
  json log_data_;
  json current_entry_;
  json current_iter_values_;
  json current_iter_timing_;
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

  template <typename T> void log_value(const T &value) {

    if (logger_ && world_.rank() == 0) {
      logger_->log_value(key_, value);
    }
  }

private:
  madness::World &world_;
  std::string key_;
  ResponseDebugLogger *logger_;
  double start_wall_, start_cpu_;
};
