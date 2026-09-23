// seed_moldft_from_dalton — build a moldft ground-state SEED by projecting the
// occupied DALTON/molden MOs onto the MRA basis and writing a `.restartdata`
// archive that moldft `--restart` resumes from. Ground-state analogue of
// tools/seed_from_dalton (which seeds the ES/response bundle).
//
// The header is written through RestartMetadata (madness/chem/Restart.h), the
// same struct SCF::save_mos and SCF::load_mos use, so the layout cannot drift
// away from moldft's and is not restated here. Orbitals follow the header:
//   uint nmo; Tensor eps; Tensor occ; vector<int> set; Function mo[nmo];
// (closed-shell only: spin_restricted=true, no beta block.)
//
// Constraints / couplings:
//   * NP=1 only (matches seed_from_dalton; single-client projection).
//   * L MUST equal the target moldft run's box (load_mos THROWS otherwise);
//     pass --L to match the deck's `l` (default 200, the value used here).
//   * k/thresh mismatch is handled by load_mos (it re-projects), so seeding at
//     the first protocol rung (k6/1e-4) is fine even for a 1e-4,1e-6 climb.
//   * xc/localize are recorded in the header but are NOT pushed into the run's
//     parameters — load_mos reads them and discards them, so the deck must set
//     whatever the run intends (defaults hf / canon). The values written here
//     are provenance, not configuration.
//   * MADNESS re-converges to its own MRA-HF minimum; the DALTON orbitals are a
//     starting guess, not frozen.
//
// Usage:
//   seed_moldft_from_dalton --molden=molden.inp --n-occ=N --out=PREFIX
//                           [--L=200] [--thresh=1e-4] [--energy=<E_h>]
//                           [--xc=hf] [--localize=canon] [--nio=1]
//   -> writes PREFIX.restartdata ; run:  moldft (with `restart 1`, prefix PREFIX)

#include "../solvers/dalton_gs_seed.hpp"   // write_gs_seed_from_molden (library form)
#include "../ResponseProtocol.hpp"
#include <madness/mra/mra.h>
#include <madness/world/MADworld.h>
#include <string>
using namespace madness;
using namespace molresponse_v3;

int main(int argc, char** argv) {
    World& world = initialize(argc, argv);
    startup(world, argc, argv, true);
    commandlineparser parser(argc, argv);

    if (world.size() != 1) {
        if (world.rank() == 0)
            print("ERROR: seed_moldft_from_dalton must run on a single rank (NP=1).");
        finalize();
        return 2;
    }
    if (!parser.key_exists("molden") || !parser.key_exists("n-occ") ||
        !parser.key_exists("out")) {
        if (world.rank() == 0) {
            print("Usage: seed_moldft_from_dalton --molden=molden.inp --n-occ=N "
                  "--out=PREFIX");
            print("  [--L=200] [--thresh=1e-4] [--energy=<E_h>] [--xc=hf] "
                  "[--localize=canon] [--nio=1]");
        }
        finalize();
        return 2;
    }

    const std::string molden_path = parser.value_raw("molden");
    const std::string out_prefix  = parser.value_raw("out");
    const int    n_occ  = std::stoi(parser.value("n-occ"));
    const double L      = parser.key_exists("L") ? std::stod(parser.value("L")) : 200.0;
    const double thresh = parser.key_exists("thresh") ? std::stod(parser.value("thresh")) : 1e-4;
    const double energy = parser.key_exists("energy") ? std::stod(parser.value("energy")) : 0.0;
    const std::string xc  = parser.key_exists("xc") ? parser.value("xc") : "hf";
    const std::string loc = parser.key_exists("localize") ? parser.value("localize") : "canon";
    const int    nio    = parser.key_exists("nio") ? std::stoi(parser.value("nio")) : 1;

    {
        GsSeedOptions opt;
        opt.L = L; opt.thresh = thresh; opt.energy = energy;
        opt.xc = xc; opt.localize = loc; opt.nio = nio;
        auto rep = write_gs_seed_from_molden(world, molden_path, n_occ, out_prefix, opt);
        if (world.rank() == 0) {
            print("  max |S_ij| (i!=j) pre-orthonormalization =", rep.max_offdiag_pre,
                  "  |proj MO norm - 1| max (post) =", rep.max_norm_dev);
            print("  Run:  moldft  (input with `restart 1`, `prefix ", out_prefix,
                  "`, `l ", L, "`, `xc ", xc, "`) -> SCF resumes from seed.");
        }
    }  // functions destruct before finalize

    finalize();
    return 0;
}
