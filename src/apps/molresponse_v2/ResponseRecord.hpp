#pragma once
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <madness/external/nlohmann_json/json.hpp>
#include <map>
#include <sstream>
#include <string>

#include "ResponseState.hpp"

using json = nlohmann::json;
namespace fs = std::filesystem;

class ResponseRecord {
 public:
  explicit ResponseRecord(World &world, const std::string &filepath)
      : path_(filepath) {
    if (fs::exists(path_)) {
      std::string json_string;
      std::ifstream in(path_);
      if (in.is_open()) {
        std::stringstream buffer;
        buffer << in.rdbuf();
        json_string = buffer.str();
        in.close();
      } else {
        std::cerr << "Error opening file: " << path_ << std::endl;
      }
      world.gop.fence();
      world.gop.broadcast_serializable(json_string, 0);

      data_ = json::parse(json_string);
      world.gop.fence();
    } else {
      data_["states"] = json::object();
    }
  }

  void initialize_states(const std::vector<LinearResponseDescriptor> &states) {
    for (const auto &state : states) {
      const std::string id = describe_perturbation(state.perturbation);

      const auto &frequencies = state.frequencies;
      const auto &protocols = state.thresholds;

      auto &state_entry = data_["states"][id];

      for (const auto &protocol : protocols) {
        std::string pstr = protocol_to_string(protocol);

        // Initialize protocol if missing
        if (!state_entry["protocols"].contains(pstr)) {
          for (const auto &freq : frequencies) {
            std::string fstr = frequency_to_string(freq);
            state_entry["protocols"][pstr]["saved"][fstr] = false;
            state_entry["protocols"][pstr]["converged"][fstr] = false;
          }
        }
      }

      // Final convergence flag (if not present)
    }
    write();
  }



// ---- helpers ----
static inline const char* yn_icon(bool v) { return v ? "✔" : "✘"; }

// Protocols like 1e-4 → scientific. `digits_after_decimal = 0` prints "1e-04".
static inline std::string fmt_sci(double x, int digits_after_decimal = 0) {
  std::ostringstream os;
  os << std::scientific << std::setprecision(digits_after_decimal) << x;
  return os.str(); // e.g., "1e-04"
}

// Frequencies: exactly N decimals, keep trailing zeros.
static inline std::string fmt_fixed(double x, int decimals = 3) {
  std::ostringstream os;
  os << std::fixed << std::setprecision(decimals) << x;
  return os.str(); // e.g., "1.250"
}

// Parse JSON object keys (strings) into numeric keys, keep original for lookup.
static inline std::vector<std::pair<double, std::string>>
numeric_keys(const nlohmann::json& obj) {
  std::vector<std::pair<double, std::string>> out;
  out.reserve(obj.size());
  for (const auto& kv : obj.items()) {
    out.emplace_back(std::stod(kv.key()), kv.key());
  }
  std::sort(out.begin(), out.end(),
            [](const auto& a, const auto& b) { return a.first < b.first; });
  return out;
}

// ---- pretty table ----
// proto_digits = 0 → "1e-04"; freq_decimals = 3 → "0.000"
void print_summary(int proto_digits = 0, int freq_decimals = 3) const {
  using std::cout;
  using std::left;
  using std::setw;

  // widths tuned for neat columns
  constexpr int W_ROW   = 5;
  constexpr int W_STATE = 32;
  constexpr int W_PROTO = 12; // "1e-04"
  constexpr int W_FREQ  = 10; // "0.000"
  constexpr int W_SAVED = 7;
  constexpr int W_CONV  = 10;

  cout << "📋 Response State Summary\n";
  cout << setw(W_ROW) << "#" << "  "
       << setw(W_STATE) << left << "State"
       << setw(W_PROTO) << "Protocol"
       << setw(W_FREQ)  << "Freq"
       << setw(W_SAVED) << "Saved"
       << setw(W_CONV)  << "Converged"
       << "\n";

  cout << std::string(W_ROW + 2 + W_STATE + W_PROTO + W_FREQ + W_SAVED + W_CONV, '-') << "\n";

  size_t row = 0;

  if (!data_.contains("states") || !data_["states"].is_object()) {
    cout << "(no states)\n";
    return;
  }

  for (const auto& [state_id, entry] : data_["states"].items()) {
    const auto& protos = (entry.contains("protocols") && entry["protocols"].is_object())
                           ? entry["protocols"]
                           : nlohmann::json::object();

    if (protos.empty()) {
      cout << setw(W_ROW) << row++ << "  "
           << setw(W_STATE) << left << state_id
           << setw(W_PROTO) << "-"
           << setw(W_FREQ)  << "-"
           << setw(W_SAVED) << "-"
           << setw(W_CONV)  << "-"
           << "\n";
      continue;
    }

    // numeric sort of protocols
    auto proto_keys = numeric_keys(protos);
    for (const auto& [proto_val, proto_key] : proto_keys) {
      const auto& pd = protos.at(proto_key);

      const auto& saved_map     = (pd.contains("saved")     && pd["saved"].is_object())     ? pd["saved"]     : nlohmann::json::object();
      const auto& converged_map = (pd.contains("converged") && pd["converged"].is_object()) ? pd["converged"] : nlohmann::json::object();

      // union of frequency keys (some may appear only under converged)
      nlohmann::json union_obj = nlohmann::json::object();
      for (const auto& kv : saved_map.items())     union_obj[kv.key()] = true;
      for (const auto& kv : converged_map.items()) union_obj[kv.key()] = true;

      auto freq_keys = numeric_keys(union_obj); // numeric, ascending

      if (freq_keys.empty()) {
        cout << setw(W_ROW) << row++ << "  "
             << setw(W_STATE) << left << state_id
             << setw(W_PROTO) << fmt_sci(proto_val, proto_digits)
             << setw(W_FREQ)  << "-"
             << setw(W_SAVED) << "-"
             << setw(W_CONV)  << "-"
             << "\n";
        continue;
      }

      for (const auto& [freq_val, freq_key] : freq_keys) {
        const bool saved     = saved_map.contains(freq_key)     ? saved_map.at(freq_key).get<bool>()     : false;
        const bool converged = converged_map.contains(freq_key) ? converged_map.at(freq_key).get<bool>() : false;

        cout << setw(W_ROW) << row++ << "  "
             << setw(W_STATE) << left << state_id
             << setw(W_PROTO) << fmt_sci(proto_val, proto_digits)
             << setw(W_FREQ)  << fmt_fixed(freq_val, freq_decimals)
             << setw(W_SAVED) << yn_icon(saved)
             << setw(W_CONV)  << yn_icon(converged)
             << "\n";
      }
    }
  }
}
  [[nodiscard]] bool is_saved(const std::string &state_id, const double & protocol,
                              const double& freq) const {
    return get_flag(state_id, protocol, freq, "saved");
  }

  [[nodiscard]] bool is_saved(const LinearResponseDescriptor& state) const {
    return get_flag(state.perturbationDescription(), state.current_threshold(), state.current_frequency(), "saved");
  }

  [[nodiscard]] bool is_converged(const std::string &state_id, const double& protocol,
                                  const double freq) const {
    return get_flag(state_id, protocol, freq, "converged");
  }

  [[nodiscard]] bool is_converged(const LinearResponseDescriptor&state) const {
    return get_flag(state.perturbationDescription(), state.current_threshold(), state.current_frequency(), "converged");
  }

  void mark_saved(const std::string &state_id, const double protocol, const double freq) {
    set_flag(state_id, protocol, freq, "saved", true);
    write();
  }

  void mark_converged(const std::string &state_id, const double protocol, const double freq,
                      const bool converged) {
    set_flag(state_id, protocol, freq, "converged", converged);
    write();
  }
  void record_status(const std::string &state_id, const double protocol,
                     const double freq, const bool converged) {
    set_flag(state_id, protocol, freq, "saved", true);
    set_flag(state_id, protocol, freq, "converged", converged);
    write();
  }
  void record_status(const LinearResponseDescriptor& state, const bool converged) {
    set_flag(state.perturbationDescription(), state.current_threshold(), state.current_frequency(), "saved", true);
    set_flag(state.perturbationDescription(), state.current_threshold(), state.current_frequency(), "converged", converged);
    write();
  }


  void mark_final_saved(const std::string &state_id, bool flag = true) {
    data_["states"][state_id]["final_saved"] = flag;
    write();
  }

  [[nodiscard]] bool final_saved(const std::string &state_id) const {
    return data_["states"][state_id].value("final_saved", false);
  }

  void write() const {
    std::ofstream out(path_);
    out << std::setw(2) << data_ << "\n";
  }
  json to_json() const { return data_; }

 private:
  std::string path_;
  json data_;

  static std::string to_string(double val) {
    std::ostringstream oss;
    oss << std::scientific << val;
    return oss.str();
  }

  [[nodiscard]] bool get_flag(const std::string &state_id, double protocol,
                              double freq, const std::string &key) const {
    std::string p_str = to_string(protocol);
    std::string f_str = to_string(freq);

    if (!data_["states"].contains(state_id)) return false;
    if (!data_["states"][state_id]["protocols"].contains(p_str)) return false;
    if (!data_["states"][state_id]["protocols"][p_str][key].contains(f_str))
      return false;

    return data_["states"][state_id]["protocols"][p_str][key][f_str]
        .get<bool>();
  }

  void set_flag(const std::string &state_id, double protocol, double freq,
                const std::string &key, bool value) {
    std::string p_str = to_string(protocol);
    std::string f_str = to_string(freq);

    data_["states"][state_id]["protocols"][p_str][key][f_str] = value;
    write();
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
