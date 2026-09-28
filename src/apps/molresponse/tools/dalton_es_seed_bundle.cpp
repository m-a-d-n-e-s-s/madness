// dalton_es_seed_bundle — write the PRODUCTION excited-state seed bundle
// (solvers/dalton_import.hpp: seed_es_from_dalton, the in-solver path madqc
// runs when a deck carries io.dalton.dir + excited.enable) standalone, into a
// calc dir of your choice, without running any solver. Purpose: keep a pristine
// iteration-0 seed for diagnostics (tools/es_overlap against a converged cold
// bundle) — the solver overwrites es__<key>/ every iteration.
//
// Loads the MADNESS ground state from the moldft archive (molecule + fock.json
// from the calc_info next to it, as tests/test_calc_manager_run.cpp does),
// configures the protocol (L from the archive header, thresh/k from --thresh),
// prepares the GS, locates the DALTON directory, and calls seed_es_from_dalton.
//
// Usage:
//   dalton_es_seed_bundle --archive=<moldft restartdata> --dalton-dir=DIR
//                         --calc-dir=DIR --roots=N [--full] [--thresh=1e-6]
//                         [--ydalton] [--hdf5]
//   --full     Full (X,Y) bundle (the RPA seed); default TDA (X only)
//   --ydalton  take the y block from DALTON (default: y = 0, seed.es_y=zero)
//   --hdf5     write the bundle with the HDF5 backend (np-portable)

#include "../GroundState.hpp"
#include "../ResponseProtocol.hpp"
#include "../ResponsePropertyPlanner.hpp"   // ESRequest
#include "../solvers/dalton_import.hpp"
#include "../solvers/function_hdf5_io.hpp"  // set_hdf5_io_enabled

#include <nlohmann/json.hpp>
#include <madness/mra/mra.h>
#include <madness/world/MADworld.h>

#include <cmath>
#include <filesystem>
#include <fstream>
#include <string>
#include <vector>

using namespace madness;
using namespace molresponse_v3;

int main(int argc, char **argv) {
  World &world = initialize(argc, argv);
  startup(world, argc, argv, true);
  commandlineparser parser(argc, argv);

  if (!parser.key_exists("archive") || !parser.key_exists("dalton-dir") ||
      !parser.key_exists("calc-dir") || !parser.key_exists("roots")) {
    if (world.rank() == 0)
      print("Usage: dalton_es_seed_bundle --archive=<restartdata> --dalton-dir=DIR "
            "--calc-dir=DIR --roots=N [--full] [--thresh=1e-6] [--ydalton] [--hdf5]");
    finalize();
    return 2;
  }
  const std::string archive_path = parser.value_raw("archive");
  const std::string dalton_dir   = parser.value_raw("dalton-dir");
  const std::string calc_dir     = parser.value_raw("calc-dir");
  const int         n_roots      = std::stoi(parser.value("roots"));
  const bool        full         = parser.key_exists("full");
  const bool        y_dalton     = parser.key_exists("ydalton");
  const double      thresh       = parser.key_exists("thresh")
                                       ? std::stod(parser.value("thresh")) : 1e-6;
  // set_hdf5_io_enabled lives inside function_hdf5_io.hpp's MADNESS_HAS_HDF5
  // block, so the call has to be guarded the same way every other caller
  // guards it (main.cpp, madqc_adapter.hpp, madqc.cpp) or a build configured
  // without HDF5 -- which is what CI builds -- fails to compile.
  if (parser.key_exists("hdf5")) {
#ifdef MADNESS_HAS_HDF5
    set_hdf5_io_enabled(true);
#else
    if (world.rank() == 0)
      print("[dalton_es_seed_bundle] --hdf5 ignored: this binary was built "
            "without HDF5 (configure with -DMADNESS_ENABLE_HDF5=ON); writing "
            "native archives instead");
#endif
  }

  {
    auto header = GroundState::read_archive_header(world, archive_path);
    const int k = parser.key_exists("thresh") ? default_k_for_thresh(thresh) : header.k;
    set_response_protocol(world, header.L, thresh, k);

    Molecule molecule;
    const auto archive_dir = std::filesystem::path(archive_path).parent_path();
    std::string fock_json;
    for (const auto &entry : std::filesystem::directory_iterator(archive_dir)) {
      const auto name = entry.path().filename().string();
      if (name.size() > 15 && name.compare(name.size() - 15, 15, ".calc_info.json") == 0 &&
          molecule.natom() == 0) {
        std::ifstream ifs(entry.path());
        nlohmann::json j; ifs >> j;
        nlohmann::json mol_json;
        if (j.contains("tasks") && j["tasks"].is_array() && !j["tasks"].empty())
          mol_json = j["tasks"][0]["molecule"];
        else if (j.contains("molecule")) mol_json = j["molecule"];
        if (!mol_json.is_null()) molecule.from_json(mol_json);
      }
      if (name.size() > 10 && name.compare(name.size() - 10, 10, ".fock.json") == 0)
        fock_json = entry.path().string();
    }
    MADNESS_CHECK_THROW(molecule.natom() > 0, "no *.calc_info.json with a molecule next to the archive");

    GroundState gs = GroundState::from_archive(world, archive_path, molecule);
    auto coulop = poperatorT(CoulombOperatorPtr(world, gs.params().lo(), 0.001 * thresh));
    gs.prepare(world, 0.001 * thresh, coulop, fock_json);

    std::filesystem::create_directories(calc_dir);
    world.gop.fence();
    DaltonManifest m = locate_dalton_dir(dalton_dir, calc_dir + "/dalton_import");
    if (world.rank() == 0)
      print("[dalton_es_seed_bundle] molden =", m.molden_path, " rspvec =", m.rspvec_path,
            " roots =", n_roots, " type =", full ? "full" : "tda",
            " y =", y_dalton ? "dalton" : "zero", " key =", protocol_key());

    std::vector<ESRequest> reqs;
    ESRequest r; r.tda = !full; r.n_roots = n_roots; r.protocols = {thresh};
    reqs.push_back(r);
    auto rep = seed_es_from_dalton(world, gs, reqs, calc_dir, m, std::sqrt(2.0), y_dalton);
    if (world.rank() == 0)
      print("[dalton_es_seed_bundle] seeded =", rep.n_seeded, " -> ", calc_dir + "/es__" + protocol_key());
  }
  finalize();
  return 0;
}
