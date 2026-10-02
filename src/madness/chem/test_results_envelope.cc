//
// SCF task record: scf_data accessors and the
// SCFResults / ConvergenceResults envelope fields round-trip through JSON.
// Pure data — no Function arithmetic — so it runs in seconds on one rank.
//
#include <madness/mra/mra.h>
#include <madness/chem/SCF.h>
#include <madness/chem/Results.h>
#include <madness/world/test_utilities.h>
#include <nlohmann/json.hpp>
#include <cmath>
#include <map>
#include <string>

using namespace madness;

namespace {

int test_scf_data_accessors() {
    test_output t("scf_data: iterations() and last()");
    scf_data d;
    t.checkpoint(d.iterations() == 0, "fresh scf_data has 0 iterations");
    t.checkpoint(d.last().empty(),    "fresh scf_data has no last() values");
    d.add_data({{"e_kinetic", 1.0}, {"e_local", 0.0}, {"e_nuclear", -2.0}, {"e_coulomb", 0.5},
                {"e_pcm", 0.0}, {"e_disp", 0.0}, {"e_xc", -0.25}, {"e_nrep", 0.75}, {"e_tot", 0.0}});
    d.add_data({{"e_kinetic", 1.1}, {"e_local", 0.0}, {"e_nuclear", -2.1}, {"e_coulomb", 0.6},
                {"e_pcm", 0.0}, {"e_disp", 0.0}, {"e_xc", -0.26}, {"e_nrep", 0.80}, {"e_tot", 0.1}});
    t.checkpoint(d.iterations() == 2, "two add_data calls -> 2 iterations");
    const auto last = d.last();
    t.checkpoint(std::fabs(last.at("e_kinetic") - 1.1) < 1e-15, "last() returns the latest e_kinetic");
    t.checkpoint(std::fabs(last.at("e_nrep") - 0.80) < 1e-15,   "last() returns the latest e_nrep");
    return t.end();
}

int test_results_envelope_roundtrip() {
    test_output t("SCFResults/ConvergenceResults envelope round-trip");
    SCFResults s;
    s.model = "scf";
    s.xc = "hf";
    s.scf_total_energy = -76.0;
    s.scf_iterations = 12;
    s.precision = {{"k", 8}, {"thresh", 1e-6}, {"protocol", {1e-4, 1e-6}},
                   {"econv", 1e-6}, {"dconv", 1e-4}, {"L", 200.0}, {"ncoeff", 123456}};
    s.energies = {{"nuclear_repulsion_energy", 9.19}, {"scf_one_electron_energy", -123.0},
                  {"scf_two_electron_energy", 38.0}, {"scf_kinetic_energy", 76.0}};
    Molecule mol; mol.add_atom(0.0, 0.0, 0.0, 2.0, 2);   // He — to_json needs a molecule
    s.scf_molecule = mol;
    const nlohmann::json j = s.to_json();
    t.checkpoint(j.value("xc", "") == "hf",                          "xc emitted");
    t.checkpoint(j["precision"]["k"] == 8,                            "precision.k emitted");
    t.checkpoint(j["precision"]["protocol"].size() == 2,              "precision.protocol emitted");
    t.checkpoint(j.contains("nuclear_repulsion_energy"),              "energies emitted FLAT under QCSchema names");
    t.checkpoint(j.value("scf_iterations", -1) == 12,                 "scf_iterations emitted");
    SCFResults back(j);
    t.checkpoint(back.xc == "hf",                                     "xc round-trips");
    t.checkpoint(back.scf_iterations == 12,                           "scf_iterations round-trips");
    t.checkpoint(back.precision["ncoeff"] == 123456,                  "precision round-trips");
    t.checkpoint(std::fabs(back.energies["scf_kinetic_energy"].get<double>() - 76.0) < 1e-15,
                 "energies round-trip");

    // xc left at its (empty) default -> not set, so the key must not appear.
    SCFResults no_xc;
    no_xc.scf_molecule = mol;
    t.checkpoint(!no_xc.to_json().contains("xc"),                     "unset xc key is omitted");

    // Old checkpoint: none of the new keys present -> defaults, no throw.
    nlohmann::json old = j;
    for (const char* k : {"xc", "precision", "scf_iterations", "nuclear_repulsion_energy",
                          "scf_one_electron_energy", "scf_two_electron_energy", "scf_kinetic_energy"})
        old.erase(k);
    SCFResults legacy(old);
    t.checkpoint(legacy.xc.empty() && legacy.scf_iterations == -1 && legacy.precision.is_null()
                 && legacy.energies.empty(),                          "legacy checkpoint loads with defaults");

    ConvergenceResults c;
    c.set_converged_thresh(1e-6).set_converged_dconv(1e-4);
    c.iterations = 12; c.status = "converged";
    const nlohmann::json cj = c.to_json();
    t.checkpoint(cj.value("iterations", -1) == 12 && cj.value("status", "") == "converged",
                 "ConvergenceResults emits iterations + status");
    ConvergenceResults cb(cj);
    t.checkpoint(cb.iterations == 12 && cb.status == "converged",     "ConvergenceResults round-trips");
    ConvergenceResults clegacy(nlohmann::json{{"converged_for_thresh", 1e-6}, {"converged_for_dconv", 1e-4}});
    t.checkpoint(clegacy.iterations == -1 && clegacy.status == "unknown", "legacy ConvergenceResults defaults");
    return t.end();
}

} // namespace

int main(int argc, char** argv) {
    World& world = madness::initialize(argc, argv);
    startup(world, argc, argv);
    int error = 0;
    error += test_scf_data_accessors();
    error += test_results_envelope_roundtrip();
    finalize();
    return error;
}
