// =========================================================================
// Round-trip test for ResponseMetadata (Inc 13b).
//
// Pure C++ — no MPI, no MADNESS World. Covers:
//   1. fresh load_or_create seeds the schema
//   2. upsert protocol / FD state / ES bundle / property
//   3. save -> reload preserves everything bit-for-bit
//   4. atomic write leaves no *.tmp lingering
//   5. set_protocol upsert: last write wins
//   6. add_property appends to an array (multiple records per protocol)
//   7. freq_key formatting matches doc 13 ("f%.5f" -> "f0.05700")
// =========================================================================

#include "../solvers/response_metadata.hpp"

#include <unistd.h>

#include <cstdio>
#include <filesystem>
#include <string>

namespace {

int failed = 0;

#define EXPECT(cond, label)                                                \
  do {                                                                     \
    if (cond) { std::printf("  [PASS]  %s\n", label); }                    \
    else      { std::printf("  [FAIL]  %s\n", label); ++failed; }          \
  } while (0)

} // namespace

int main() {
  using molresponse_v3::ResponseMetadata;

  // Use a process-unique scratch path so parallel CI runs don't collide.
  const std::filesystem::path tmp =
      std::filesystem::temp_directory_path() /
      ("rmeta_test_" + std::to_string(::getpid()));
  std::filesystem::create_directories(tmp);
  const std::string path = (tmp / "response_metadata.json").string();
  std::filesystem::remove(path);

  // ---- freq_key formatting -------------------------------------------------
  std::printf("=== freq_key ===\n");
  EXPECT(ResponseMetadata::freq_key(0.057)  == "f0.05700", "freq_key(0.057)  == f0.05700");
  EXPECT(ResponseMetadata::freq_key(0.0)    == "f0.00000", "freq_key(0.0)    == f0.00000");
  EXPECT(ResponseMetadata::freq_key(0.4811) == "f0.48110", "freq_key(0.4811) == f0.48110");

  // ---- fresh load: schema is seeded ----------------------------------------
  std::printf("=== fresh load_or_create ===\n");
  {
    auto m = ResponseMetadata::load_or_create(path);
    const auto &j = m.json();
    EXPECT(j["schema_version"] == 1,             "schema_version == 1");
    EXPECT(j.contains("protocols"),              "protocols/ present");
    EXPECT(j.contains("fd_states"),              "fd_states/ present");
    EXPECT(j.contains("excited_states"),         "excited_states/ present");
    EXPECT(j.contains("properties"),             "properties/ present");
  }

  // ---- upsert + save -------------------------------------------------------
  std::printf("=== upsert + save ===\n");
  {
    auto m = ResponseMetadata::load_or_create(path);
    m.set_protocol("1e-04_k6", 1e-4, 6, 0);
    m.set_protocol("1e-06_k8", 1e-6, 8, 1);

    nlohmann::json fd_entry = {
        {"freq",         0.057},
        {"type",         "full"},
        {"shell",        "closed_shell"},
        {"converged",    true},
        {"iter",         8},
        {"bsh_residual", 1.2e-7},
        {"archive",      "dipole_x__1e-06_k8__f0.05700"}};
    m.set_fd_state("dipole_x", "1e-06_k8",
                   ResponseMetadata::freq_key(0.057), fd_entry);

    nlohmann::json es_entry = {
        {"type",             "tda"},
        {"shell",            "closed_shell"},
        {"n_roots",          2},
        {"bundle_dir",       "es_bundle__1e-06_k8"},
        {"converged",        true},
        {"slot_permutation", nlohmann::json::array({0, 1})},
        {"roots", nlohmann::json::array({
            {{"stable_index", 0}, {"root_id", "es_root_0000"}, {"slot", 0},
             {"omega", 0.468}, {"display_name", "S1"}},
            {{"stable_index", 1}, {"root_id", "es_root_0001"}, {"slot", 1},
             {"omega", 0.481}, {"display_name", "S2"}}})}};
    m.set_es_bundle("1e-06_k8", es_entry);

    m.add_property("resonant_raman", "1e-06_k8",
                   {{"es_root_id", "es_root_0001"},
                    {"fd_freq",    0.481},
                    {"value",      0.0}});

    m.save();
    EXPECT(std::filesystem::exists(path),                  "file written");
    EXPECT(!std::filesystem::exists(path + ".tmp"),        "no .tmp lingers");
  }

  // ---- reload: every field round-trips -------------------------------------
  std::printf("=== reload round-trip ===\n");
  {
    auto m = ResponseMetadata::load_or_create(path);
    const auto &j = m.json();

    EXPECT(j["protocols"]["1e-06_k8"]["thresh"] == 1e-6,           "protocol thresh");
    EXPECT(j["protocols"]["1e-06_k8"]["k"]      == 8,              "protocol k");
    EXPECT(j["protocols"]["1e-06_k8"]["index"]  == 1,              "protocol index");

    const auto &fd = j["fd_states"]["dipole_x"]["1e-06_k8"]["f0.05700"];
    EXPECT(fd["freq"]      == 0.057,                               "fd freq");
    EXPECT(fd["converged"] == true,                                "fd converged");
    EXPECT(fd["archive"]   == "dipole_x__1e-06_k8__f0.05700",      "fd archive");

    const auto &es = j["excited_states"]["1e-06_k8"];
    EXPECT(es["n_roots"]                       == 2,               "es n_roots");
    EXPECT(es["bundle_dir"]                    == "es_bundle__1e-06_k8", "es bundle_dir");
    EXPECT(es["slot_permutation"]              == nlohmann::json::array({0, 1}), "es slot_permutation");
    EXPECT(es["roots"][1]["root_id"]           == "es_root_0001", "es root[1] id");
    EXPECT(es["roots"][1]["display_name"]      == "S2",           "es root[1] display");

    const auto &props = j["properties"]["resonant_raman"]["1e-06_k8"];
    EXPECT(props.is_array() && props.size() == 1,                 "1 property record");
    EXPECT(props[0]["es_root_id"]              == "es_root_0001", "prop es_root_id");
    EXPECT(props[0]["fd_freq"]                 == 0.481,          "prop fd_freq");
  }

  // ---- upsert overwrites; add_property appends -----------------------------
  std::printf("=== upsert overwrites; properties append ===\n");
  {
    auto m = ResponseMetadata::load_or_create(path);
    // bump the index on an existing protocol
    m.set_protocol("1e-06_k8", 1e-6, 8, 99);
    // add a second property record at the same protocol
    m.add_property("resonant_raman", "1e-06_k8",
                   {{"es_root_id", "es_root_0000"},
                    {"fd_freq",    0.468},
                    {"value",      0.0}});
    m.save();

    auto m2 = ResponseMetadata::load_or_create(path);
    EXPECT(m2.json()["protocols"]["1e-06_k8"]["index"] == 99,
           "set_protocol last-write-wins");
    EXPECT(m2.json()["properties"]["resonant_raman"]["1e-06_k8"].size() == 2,
           "add_property appended (now 2 records)");
  }

  // ---- schema_version mismatch is rejected ---------------------------------
  std::printf("=== schema_version guard ===\n");
  {
    const std::string bad = (tmp / "bad.json").string();
    {
      std::ofstream o(bad);
      o << "{\"schema_version\":999,\"protocols\":{},\"fd_states\":{},"
           "\"excited_states\":{},\"properties\":{}}";
    }
    bool threw = false;
    try { (void)ResponseMetadata::load_or_create(bad); }
    catch (const std::exception &) { threw = true; }
    EXPECT(threw, "unrecognized schema_version throws");
  }

  std::filesystem::remove_all(tmp);
  std::printf("\n%s: %d failure(s)\n",
              failed == 0 ? "ALL PASS" : "FAILED", failed);
  return failed == 0 ? 0 : 1;
}
