#pragma once
#include <fstream>
#include <iomanip>
#include <madness/external/nlohmann_json/json.hpp>

#include "ResponseState.hpp"

namespace fs = std::filesystem;
using json = nlohmann::json;

class ResponseDebugLogger {
public:
  // Pass filename and whether logging is enabled.
  // If you're in MPI, pass rank0_only=true and only construct/use this on rank
  // 0, or set rank0_only=true and call write_to_disk() from all ranks safely.
  explicit ResponseDebugLogger(std::string filename, bool enabled = false,
                               bool rank0_only = true)
      : enabled_(enabled), rank0_only_(rank0_only),
        filename_(std::move(filename)) {
    if (!enabled_)
      return;
    if (fs::exists(filename_)) {
      std::ifstream inf(filename_);
      if (inf)
        inf >> log_data_;
    }
    ensure_root_();
  }

  [[nodiscard]] bool enabled() const { return enabled_; }
  void set_enabled(bool on) { enabled_ = on; }

  // Call at the start of each (state, protocol, frequency) block you’re about
  // to iterate.
  void start_state(const LinearResponseDescriptor &state) {
    if (!enabled_)
      return;
    state_key_ = state.description(); // e.g. "mu_x (restricted)"
    proto_key_ = protocol_key_(state.current_threshold()); // "1e-06"
    freq_key_ = freq_key(state.current_frequency());       // "0.500"

    // Ensure nested shape exists, no copies involved
    auto &node = node_ref_(); // creates and returns the path
    // If first time, ensure arrays exist
    if (!node.contains("iteration_values"))
      node["iteration_values"] = json::array();
    if (!node.contains("iteration_timings"))
      node["iteration_timings"] = json::array();
  }

  // Start of each solver iterate() call
  void begin_iteration(std::size_t iter_index) {
    if (!enabled_)
      return;
    current_iter_values_ = {{"iter", iter_index}, {"steps", json::object()}};
    current_iter_timing_ = {{"iter", iter_index}, {"steps", json::object()}};
  }

  // Numeric value for a named step (keep T JSON-serializable; prefer
  // double/int)
  template <typename T>
  void log_value(const std::string &step_name, const T &value) {
    if (!enabled_)
      return;
    current_iter_values_["steps"][step_name]["value"] = value;
  }

  // Timing for a named step
  void log_timing(const std::string &step_name, double wall_time,
                  double cpu_time) {
    if (!enabled_)
      return;
    auto &s = current_iter_timing_["steps"][step_name];
    s["wall_time"] = wall_time;
    s["cpu_time"] = cpu_time;
  }

  // Convenience
  template <typename T>
  void log_value_and_time(const std::string &step_name, const T &value,
                          double wall_time, double cpu_time) {
    log_value(step_name, value);
    log_timing(step_name, wall_time, cpu_time);
  }

  // End of one iterate() call → append to arrays on the real node
  void end_iteration() {
    if (!enabled_)
      return;
    auto &node = node_ref_();
    node["iteration_values"].push_back(current_iter_values_);
    node["iteration_timings"].push_back(current_iter_timing_);
  }

  // Optional: nothing to do right now, but good for symmetry
  void finalize_state() {
    if (!enabled_)
      return;
    // no-op; everything is already written into log_data_
  }

  // Persist the entire log JSON
  // If rank0_only_ == true, call this only on rank 0 (or call from all ranks
  // with identical state).
  void write_to_disk() const {
    if (!enabled_) {
      return;
    }
    std::ofstream ofs(filename_);
    ofs << std::setw(2) << log_data_;
    if (ofs) {
      std::cout << "📂 Wrote response debug log to " << filename_ << "\n";
    } else {

      std::cerr << "⚠️  Failed to write response debug log to " << filename_
                << "\n";
    }
    // Make sure the directory exists
  }

  json to_json() const { return log_data_; }

  // ——— Pretty printers ——————————————————————————————————————————————

  void print_timing_table(const LinearResponseDescriptor &state) const {
    if (!enabled_)
      return;
    const auto state_key = state.description();
    const auto proto_key = protocol_key_(state.current_threshold());
    const auto fkey = freq_key(state.current_frequency());

    if (!exists_(state_key, proto_key, fkey))
      return;
    const auto &iter_data = log_data_.at(state_key)
                                .at("protocols")
                                .at(proto_key)
                                .at("freqs")
                                .at(fkey)
                                .at("iteration_timings");

    // gather all step names
    std::set<std::string> step_names;
    for (const auto &iter : iter_data) {
      for (auto &p : iter["steps"].items())
        step_names.insert(p.key());
    }

    std::map<std::string, std::string> short_key;
    for (auto &name : step_names)
      short_key[name] = name.substr(0, 5);

    constexpr int W = 10;
    std::cout << "\n⏱️ Timing for " << state_key << " | proto=" << proto_key
              << " | freq=" << freq_key << "\n";
    std::cout << std::setw(6) << "Iter";
    for (auto &n : step_names)
      std::cout << std::setw(W) << short_key[n];
    std::cout << "\n" << std::string(6 + W * step_names.size(), '-') << "\n";

    for (const auto &iter : iter_data) {
      std::cout << std::setw(6) << iter["iter"].get<int>();
      for (auto &n : step_names) {
        if (iter["steps"].contains(n)) {
          double w = iter["steps"][n]["wall_time"].get<double>();
          std::cout << std::setw(W) << std::fixed << std::setprecision(4) << w;
        } else {
          std::cout << std::setw(W) << "N/A";
        }
      }
      std::cout << "\n";
    }
  }

  void print_values_table(const LinearResponseDescriptor &state) const {
    if (!enabled_)
      return;
    const auto state_key = state.description();
    const auto proto_key = protocol_key_(state.current_threshold());
    const auto fkey = freq_key(state.current_frequency());

    if (!exists_(state_key, proto_key, fkey))
      return;
    const auto &iter_data = log_data_.at(state_key)
                                .at("protocols")
                                .at(proto_key)
                                .at("freqs")
                                .at(fkey)
                                .at("iteration_values");

    std::set<std::string> step_names;
    for (const auto &iter : iter_data) {
      for (auto &s : iter["steps"].items())
        step_names.insert(s.key());
    }

    constexpr int width = 16;
    std::cout << "\n📋 Values for " << state_key << " | proto=" << proto_key
              << " | freq=" << freq_key << "\n";
    std::cout << std::setw(6) << "Iter";
    for (auto &n : step_names)
      std::cout << std::setw(width) << n;
    std::cout << "\n"
              << std::string(6 + width * step_names.size(), '-') << "\n";

    for (const auto &iter : iter_data) {
      std::cout << std::setw(6) << iter["iter"].get<int>();
      for (auto &n : step_names) {
        if (iter["steps"].contains(n)) {
          double v = iter["steps"][n]["value"].get<double>();
          std::cout << std::setw(width) << std::scientific
                    << std::setprecision(5) << v;
        } else {
          std::cout << std::setw(width) << "N/A";
        }
      }
      std::cout << "\n";
    }
  }

private:
  bool enabled_ = false;
  bool rank0_only_ = true;
  std::string filename_;
  json log_data_;

  // Current tuple
  std::string state_key_;
  std::string proto_key_;
  std::string freq_key_;

  json current_iter_values_;
  json current_iter_timing_;

  // ---------- Canonical keys ----------
  static std::string freq_key(double f) {
    std::ostringstream os;
    os << std::fixed << std::setprecision(3) << f;
    return os.str(); // "0.500"
  }
  static std::string protocol_key_(double thr) {
    // Compact scientific like "1e-06"
    std::ostringstream os;
    os.setf(std::ios::scientific);
    os << std::setprecision(0) << thr;
    std::string s = os.str();
    // normalize "1.e-06" -> "1e-06"
    auto dot = s.find('.');
    if (dot != std::string::npos)
      s.erase(dot, 1);
    return s;
  }

  // ---------- Structure helpers ----------
  void ensure_root_() {
    if (!log_data_.is_object())
      log_data_ = json::object();
  }
  void ensure_state_(const std::string &s) {
    if (!log_data_.contains(s))
      log_data_[s] = json::object();
    if (!log_data_[s].contains("protocols"))
      log_data_[s]["protocols"] = json::object();
  }
  void ensure_proto_(const std::string &s, const std::string &p) {
    ensure_state_(s);
    if (!log_data_[s]["protocols"].contains(p))
      log_data_[s]["protocols"][p] = json::object();
    if (!log_data_[s]["protocols"][p].contains("freqs"))
      log_data_[s]["protocols"][p]["freqs"] = json::object();
  }
  void ensure_freq_(const std::string &s, const std::string &p,
                    const std::string &f) {
    ensure_proto_(s, p);
    auto &freqnode = log_data_[s]["protocols"][p]["freqs"];
    if (!freqnode.contains(f))
      freqnode[f] = json::object();
  }

  // Returns a **reference** to the active (state,proto,freq) node
  json &node_ref_() {
    ensure_freq_(state_key_, proto_key_, freq_key_);
    return log_data_[state_key_]["protocols"][proto_key_]["freqs"][freq_key_];
  }

  bool exists_(const std::string &s, const std::string &p,
               const std::string &f) const {
    if (!log_data_.contains(s))
      return false;
    if (!log_data_.at(s).contains("protocols"))
      return false;
    if (!log_data_.at(s).at("protocols").contains(p))
      return false;
    if (!log_data_.at(s).at("protocols").at(p).contains("freqs"))
      return false;
    if (!log_data_.at(s).at("protocols").at(p).at("freqs").contains(f))
      return false;
    return true;
  }
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
