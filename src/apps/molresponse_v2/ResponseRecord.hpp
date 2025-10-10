#pragma once
#include "ResponseState.hpp"

#include <filesystem>
#include <fstream>
#include <iomanip>
#include <madness/external/nlohmann_json/json.hpp>
#include <map>
#include <sstream>
#include <string>

using json = nlohmann::json;
namespace fs = std::filesystem;

class ResponseRecord2 {
public:
  using json = nlohmann::json;

  struct MissingItem {
    std::string state;    // perturbationDescription()
    std::string freq;     // canonical "0.500"
    std::string protocol; // canonical "1e-06"
    bool saved_found = false;
    bool converged_found = false;
  };

  ResponseRecord2(World &world, const std::string &filepath)
      : world_(world), path_(filepath) {
    std::string json_string;
    if (world_.rank() == 0) {
      if (fs::exists(path_)) {
        std::ifstream in(path_);
        if (in) {
          std::stringstream buf;
          buf << in.rdbuf();
          json_string = buf.str();
        } else {
          std::cerr << "Error opening file: " << path_ << std::endl;
          json_string = "{}";
        }
      } else {
        json_string = "{}";
      }
    }
    world_.gop.fence();
    world_.gop.broadcast_serializable(json_string, 0);
    world_.gop.fence();

    if (json_string.empty())
      data_ = json::object();
    else
      data_ = json::parse(json_string, nullptr, /*allow_exceptions=*/true);

    ensure_root(data_);
  }

  // --- Initialization for a generated set of states ---
  void initialize_states(const std::vector<LinearResponseDescriptor> &states) {
    for (const auto &st : states) {
      const std::string state_id = st.perturbationDescription();
      ensure_state(data_, state_id);

      // NB: canonicalize keys up front
      for (double thr : st.thresholds) {
        const std::string pkey = protocol_key(thr);
        ensure_protocol(data_, state_id, pkey);

        for (double f : st.frequencies) {
          const std::string fkey = freq_key(f);
          auto &node = data_["states"][state_id]["protocols"][pkey];
          if (!node["saved"].contains(fkey))
            node["saved"][fkey] = false;
          if (!node["converged"].contains(fkey))
            node["converged"][fkey] = false;
        }
      }
    }
    write(); // sync to disk + broadcast
  }

  // --- Pretty table ---
  void print_summary(int proto_digits = 0, int freq_decimals = 3) const {
    using std::cout;
    using std::left;
    using std::setw;
    constexpr int W_ROW = 5, W_STATE = 32, W_PROTO = 12, W_FREQ = 10,
                  W_SAVED = 7, W_CONV = 10;

    cout << "📋 Response State Summary\n";
    cout << setw(W_ROW) << "#" << "  " << setw(W_STATE) << left << "State"
         << setw(W_PROTO) << "Protocol" << setw(W_FREQ) << "Freq"
         << setw(W_SAVED) << "Saved" << setw(W_CONV) << "Converged" << "\n";
    cout << std::string(
                W_ROW + 2 + W_STATE + W_PROTO + W_FREQ + W_SAVED + W_CONV, '-')
         << "\n";

    size_t row = 0;
    if (!data_.contains("states") || !data_["states"].is_object()) {
      cout << "(no states)\n";
      return;
    }

    for (const auto &kv : data_["states"].items()) {
      const std::string &state_id = kv.key();
      const auto &entry = kv.value();
      const auto &protos =
          (entry.contains("protocols") && entry["protocols"].is_object())
              ? entry["protocols"]
              : json::object();
      if (protos.empty()) {
        cout << setw(W_ROW) << row++ << "  " << setw(W_STATE) << left
             << state_id << setw(W_PROTO) << "-" << setw(W_FREQ) << "-"
             << setw(W_SAVED) << "-" << setw(W_CONV) << "-" << "\n";
        continue;
      }

      auto proto_keys = numeric_keys(protos); // (double, string) sorted asc
      for (const auto &[pnum, pkey] : proto_keys) {
        const auto &node = protos.at(pkey);
        const auto &saved =
            (node.contains("saved") && node["saved"].is_object())
                ? node["saved"]
                : json::object();
        const auto &conv =
            (node.contains("converged") && node["converged"].is_object())
                ? node["converged"]
                : json::object();

        json union_obj = json::object();
        for (const auto &it : saved.items())
          union_obj[it.key()] = true;
        for (const auto &it : conv.items())
          union_obj[it.key()] = true;

        auto freq_keys = numeric_keys(union_obj); // numeric sort
        if (freq_keys.empty()) {
          cout << setw(W_ROW) << row++ << "  " << setw(W_STATE) << left
               << state_id << setw(W_PROTO) << fmt_sci(pnum, proto_digits)
               << setw(W_FREQ) << "-" << setw(W_SAVED) << "-" << setw(W_CONV)
               << "-" << "\n";
          continue;
        }

        for (const auto &[fnum, fkey] : freq_keys) {
          const bool s =
              saved.contains(fkey) ? saved.at(fkey).get<bool>() : false;
          const bool c =
              conv.contains(fkey) ? conv.at(fkey).get<bool>() : false;
          cout << setw(W_ROW) << row++ << "  " << setw(W_STATE) << left
               << state_id << setw(W_PROTO) << fmt_sci(pnum, proto_digits)
               << setw(W_FREQ) << fmt_fixed(fnum, freq_decimals)
               << setw(W_SAVED) << (s ? "✔" : "✘") << setw(W_CONV)
               << (c ? "✔" : "✘") << "\n";
        }
      }
    }
  }

  // --- Queries (overloads take doubles; we canonicalize to keys) ---
  bool is_saved(const std::string &state_id, double protocol,
                double freq) const {
    return get_flag(state_id, protocol_key(protocol), freq_key(freq), "saved");
  }
  bool is_converged(const std::string &state_id, double protocol,
                    double freq) const {
    return get_flag(state_id, protocol_key(protocol), freq_key(freq),
                    "converged");
  }
  bool is_saved(const LinearResponseDescriptor &st) const {
    return get_flag(st.perturbationDescription(),
                    protocol_key(st.current_threshold()),
                    freq_key(st.current_frequency()), "saved");
  }
  bool is_converged(const LinearResponseDescriptor &st) const {
    return get_flag(st.perturbationDescription(),
                    protocol_key(st.current_threshold()),
                    freq_key(st.current_frequency()), "converged");
  }

  // --- Mutations ---
  void mark_saved(const std::string &state_id, double protocol, double freq) {
    set_flag(state_id, protocol_key(protocol), freq_key(freq), "saved", true);
  }
  void mark_converged(const std::string &state_id, double protocol, double freq,
                      bool c) {
    set_flag(state_id, protocol_key(protocol), freq_key(freq), "converged", c);
  }
  void record_status(const LinearResponseDescriptor &st, bool c) {
    const std::string sid = st.perturbationDescription();
    const std::string p = protocol_key(st.current_threshold());
    const std::string f = freq_key(st.current_frequency());
    ensure_protocol(data_, sid, p);
    data_["states"][sid]["protocols"][p]["saved"][f] = true;
    data_["states"][sid]["protocols"][p]["converged"][f] = c;
    write();
  }

  // --- Property gating helpers ---
  static std::string
  final_protocol_key_from(const std::vector<double> &protos) {
    if (protos.empty())
      return "inf"; // sentinel; will never match
    return protocol_key(*std::min_element(protos.begin(), protos.end()));
  }

  // If converged at multiple protocols, return the most accurate (smallest
  // numeric) protocol key.
  std::optional<std::string>
  best_converged_protocol(const std::string &state_id, double freq) const {
    if (!data_.contains("states"))
      return std::nullopt;
    const auto sit = data_["states"].find(state_id);
    if (sit == data_["states"].end())
      return std::nullopt;

    const auto &protos = (*sit)["protocols"];
    std::vector<std::string> pkeys;
    pkeys.reserve(protos.size());
    for (auto it = protos.begin(); it != protos.end(); ++it)
      pkeys.push_back(it.key());
    std::sort(pkeys.begin(), pkeys.end(),
              [](const std::string &a, const std::string &b) {
                return protocol_numeric(a) < protocol_numeric(b);
              });

    const std::string fk = freq_key(freq);
    for (const auto &p : pkeys) {
      const auto &node = protos.at(p);
      if (node.contains("converged")) {
        const auto &conv = node["converged"];
        if (conv.contains(fk) && conv.at(fk).get<bool>())
          return p;
      }
    }
    return std::nullopt;
  }

  // Check saved && converged at *final* protocol for each (state, freq)
  std::vector<MissingItem>
  missing_at_final_protocol(const std::vector<std::string> &state_ids,
                            const std::vector<double> &freqs,
                            const std::string &final_proto_key) const {
    std::vector<MissingItem> out;
    for (const auto &sid : state_ids) {
      for (double f : freqs) {
        const std::string fk = freq_key(f);
        bool s = false, c = false;
        if (data_.contains("states")) {
          auto sit = data_["states"].find(sid);
          if (sit != data_["states"].end()) {
            auto pit = (*sit)["protocols"].find(final_proto_key);
            if (pit != (*sit)["protocols"].end()) {
              const auto &node = (*sit)["protocols"].at(final_proto_key);
              if (node.contains("saved") && node["saved"].contains(fk))
                s = node["saved"].at(fk).get<bool>();
              if (node.contains("converged") && node["converged"].contains(fk))
                c = node["converged"].at(fk).get<bool>();
            }
          }
        }
        if (!(s && c))
          out.push_back({sid, fk, final_proto_key, s, c});
      }
    }
    return out;
  }

  // Throw if not ready; call this right before property computations.
  void enforce_ready_for_properties(const std::vector<std::string> &state_ids,
                                    const std::vector<double> &freqs,
                                    const std::string &final_proto_key) const {
    auto missing = missing_at_final_protocol(state_ids, freqs, final_proto_key);
    if (!missing.empty()) {
      if (world_.rank() == 0) {
        std::ostringstream msg;
        msg << "Property gate failed at final protocol " << final_proto_key
            << ":\n";
        for (const auto &m : missing) {
          msg << "  - " << m.state << " @ " << m.freq
              << " (saved=" << (m.saved_found ? "true" : "false")
              << ", converged=" << (m.converged_found ? "true" : "false")
              << ")\n";
        }
        std::cerr << msg.str();
      }
      throw std::runtime_error(
          "Required states/frequencies not ready at final protocol.");
    }
  }

  // --- File I/O (rank 0 writes; all ranks sync) ---
  void write() {
    std::string json_string;
    if (world_.rank() == 0) {
      std::ofstream out(path_);
      if (!out)
        throw std::runtime_error("Cannot open " + path_ + " for writing");
      out << std::setw(2) << data_ << "\n";
      json_string = data_.dump();
    }
    world_.gop.fence();
    world_.gop.broadcast_serializable(json_string, 0);
    world_.gop.fence();
    if (world_.rank() != 0) {
      data_ = json::parse(json_string);
    }
  }

  json to_json() const { return data_; }

  // --- Optional: “final_saved” compatibility flag you had ---
  void mark_final_saved(const std::string &state_id, bool flag = true) {
    ensure_state(data_, state_id);
    data_["states"][state_id]["final_saved"] = flag;
    write();
  }

  // --- Flat rows (useful for CSV/logs) ---
  struct Row {
    std::string state, freq, protocol;
    bool saved = false, converged = false;
  };
  std::vector<Row> to_rows() const {
    std::vector<Row> rows;
    if (!data_.contains("states"))
      return rows;
    for (auto sit = data_["states"].begin(); sit != data_["states"].end();
         ++sit) {
      const std::string sid = sit.key();
      if (!(*sit).contains("protocols"))
        continue;
      const auto &protos = (*sit)["protocols"];
      for (auto pit = protos.begin(); pit != protos.end(); ++pit) {
        const std::string pkey = pit.key();
        const auto &node = pit.value();
        const auto &saved = node.value("saved", json::object());
        const auto &conv = node.value("converged", json::object());
        std::map<std::string, bool> k;
        for (auto it = saved.begin(); it != saved.end(); ++it)
          k[it.key()] = true;
        for (auto it = conv.begin(); it != conv.end(); ++it)
          k[it.key()] = true;
        for (const auto &kv : k) {
          const auto &fk = kv.first;
          bool s = saved.contains(fk) ? saved.at(fk).get<bool>() : false;
          bool c = conv.contains(fk) ? conv.at(fk).get<bool>() : false;
          rows.push_back({sid, fk, pkey, s, c});
        }
      }
    }
    std::sort(rows.begin(), rows.end(), [](const Row &a, const Row &b) {
      if (a.state != b.state)
        return a.state < b.state;
      double af = std::stod(a.freq), bf = std::stod(b.freq);
      if (af != bf)
        return af < bf;
      return protocol_numeric(a.protocol) < protocol_numeric(b.protocol);
    });
    return rows;
  }
  // -------- Canonicalization & shape --------
  static std::string freq_key(double f) {
    std::ostringstream os;
    os << std::fixed << std::setprecision(3) << f;
    return os.str();
  }
  static std::string protocol_key(double thr) {
    char buf[64];
    std::snprintf(buf, sizeof(buf), "%.0e", thr); // "1e-06"
    std::string s(buf);
    // normalize "+": "1e+00" -> "1e+00" (keep plus; stod handles it)
    return s;
  }
  static double protocol_numeric(const std::string &p) {
    try {
      return std::stod(p);
    } catch (...) {
      return std::numeric_limits<double>::infinity();
    }
  }

private:
  World &world_;
  std::string path_;
  json data_;

  static void ensure_root(json &root) {
    if (!root.is_object())
      root = json::object();
    if (!root.contains("states"))
      root["states"] = json::object();
  }
  static void ensure_state(json &root, const std::string &sid) {
    ensure_root(root);
    auto &states = root["states"];
    if (!states.contains(sid)) {
      states[sid] = json::object();
      states[sid]["protocols"] = json::object();
    } else if (!states[sid].contains("protocols")) {
      states[sid]["protocols"] = json::object();
    }
  }
  static void ensure_protocol(json &root, const std::string &sid,
                              const std::string &pkey) {
    ensure_state(root, sid);
    auto &protos = root["states"][sid]["protocols"];
    if (!protos.contains(pkey)) {
      protos[pkey] = json::object(
          {{"saved", json::object()}, {"converged", json::object()}});
    }
  }

  // -------- Printing helpers / sorting --------
  static std::string fmt_sci(double x, int digits_after_decimal = 0) {
    std::ostringstream os;
    os << std::scientific << std::setprecision(digits_after_decimal) << x;
    return os.str();
  }
  static std::string fmt_fixed(double x, int decimals = 3) {
    std::ostringstream os;
    os << std::fixed << std::setprecision(decimals) << x;
    return os.str();
  }
  static std::vector<std::pair<double, std::string>>
  numeric_keys(const json &obj) {
    std::vector<std::pair<double, std::string>> out;
    out.reserve(obj.size());
    for (const auto &kv : obj.items())
      out.emplace_back(std::stod(kv.key()), kv.key());
    std::sort(out.begin(), out.end(),
              [](const auto &a, const auto &b) { return a.first < b.first; });
    return out;
  }

  // -------- core flag ops (with canonical keys precomputed) --------
  bool get_flag(const std::string &sid, const std::string &pkey,
                const std::string &fkey, const std::string &which) const {
    if (!data_.contains("states"))
      return false;
    const auto sit = data_["states"].find(sid);
    if (sit == data_["states"].end())
      return false;
    const auto pit = (*sit)["protocols"].find(pkey);
    if (pit == (*sit)["protocols"].end())
      return false;
    const auto &node = (*pit);
    if (!node.contains(which) || !node[which].is_object())
      return false;
    const auto &sub = node[which];
    if (!sub.contains(fkey))
      return false;
    return sub.at(fkey).get<bool>();
  }

  void set_flag(const std::string &sid, const std::string &pkey,
                const std::string &fkey, const std::string &which, bool value) {
    ensure_protocol(data_, sid, pkey);
    data_["states"][sid]["protocols"][pkey][which][fkey] = value;
    write();
  }
};
