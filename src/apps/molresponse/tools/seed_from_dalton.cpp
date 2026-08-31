// seed_from_dalton — build a molresponse TDA excited-state SEED bundle by
// reading a DALTON RSPVEC binary and Molden file directly (no Python/HDF5
// intermediate). Projects the per-occupied-orbital response function
//   x_i(r) = sum_mu D[mu,i] chi_mu(r),  D = C_vir @ X.T
// into MADNESS Function<double,3> via a DaltonResponseFunctor, wraps all N
// functions as a one-root ESSolver<TDA,ClosedShell>::State, and calls
// save_es_roots(..., converged=false). The resulting bundle is identical in
// schema to the HDF5-seeded one produced by seed_es_from_hdf5; the ES node
// Resumes from it on the next run.
//
// Constraints: NP=1 only (same as seed_es_from_hdf5).
//
// --archive=<moldft restartdata>: load the MADNESS ground state and rotate
// every root's x (and y with --full) into ITS occupied gauge with
// U = M (M^T M)^{-1/2}, M_ji = <phi^DAL_j|phi^MAD_i> (shared helper
// occupied_gauge_rotation, tools/dalton_mra.hpp). DALTON pairs response
// vectors with its CANONICAL occupied MOs; a localized MADNESS ground state
// (moldft `localize new`) makes an unrotated seed silently worthless
// (measured on the h2o FD A/B: seed-implied alpha_zz 3.56 vs 8.39, zero
// iteration savings). With --archive the box L and k come from the archive
// header (--thresh still overrides; then k follows the thresh table).
// WITHOUT --archive the legacy behavior is kept — cell hardcoded to
// [-200,200]^3, NO rotation — and a loud warning is printed: that seed is
// only valid against a canonically-localized ground state (localize canon).
//
// Usage:
//   seed_from_dalton --rspvec=RSPVEC --molden=molden.inp --n-occ=N
//                    --root=0 --omega=<au> --calc-dir=DIR
//                    [--archive=<moldft restartdata>]
//                    [--thresh=1e-4] [--scale=1.41421356]

#include "dalton_rspvec.hpp"
#include "dalton_gto.hpp"
#include "dalton_mra.hpp"    // shared DaltonResponseFunctor + projection helpers
                             // + occupied_gauge_rotation

#include "../GroundState.hpp"
#include "../ResponseProtocol.hpp"
#include "../kernels/tags.hpp"
#include "../solvers/response_state.hpp"
#include "../solvers/es_solver.hpp"
#include "../solvers/es_save_load.hpp"

#include <nlohmann/json.hpp>
#include <madness/mra/mra.h>
#include <madness/world/MADworld.h>

#include <cmath>
#include <filesystem>
#include <fstream>
#include <memory>
#include <optional>
#include <sstream>
#include <string>
#include <vector>

using namespace madness;
using namespace molresponse_v3;

// DaltonResponseFunctor + the AO->MRA projection helpers now live in
// dalton_mra.hpp (shared with tpa_from_dalton and the dalton.dir FD seed
// path, solvers/dalton_import.hpp).

int main(int argc, char** argv) {
    World& world = initialize(argc, argv);
    startup(world, argc, argv, true);

    commandlineparser parser(argc, argv);

    if (world.size() != 1) {
        if (world.rank() == 0)
            print("ERROR: seed_from_dalton must run on a single rank (NP=1).");
        finalize();
        return 2;
    }

    if (!parser.key_exists("rspvec") || !parser.key_exists("molden") ||
        !parser.key_exists("n-occ")  || !parser.key_exists("calc-dir") ||
        !(parser.key_exists("root") || parser.key_exists("roots")) ||
        !(parser.key_exists("omega") || parser.key_exists("omegas"))) {
        if (world.rank() == 0) {
            print("Usage: seed_from_dalton --rspvec=RSPVEC --molden=molden.inp "
                  "--n-occ=N --root=IDX --omega=<au> --calc-dir=DIR");
            print("  multi-root: --roots=0,1,2 --omegas=w0,w1,w2 (one N-root bundle)");
            print("  [--archive=<moldft restartdata>] (load the MADNESS ground "
                  "state, take box/k from its header, and ROTATE the seed into "
                  "its occupied gauge — required for a localized ground state)");
            print("  [--full] (write Full X,Y bundle for --es-full/RPA solve; "
                  "default TDA X-only) [--yflip] (flip Y sign) "
                  "[--thresh=1e-4] [--scale=<sqrt2>]");
        }
        finalize();
        return 2;
    }

    const std::string rspvec_path = parser.value_raw("rspvec");
    const std::string molden_path = parser.value_raw("molden");
    const std::string calc_dir   = parser.value_raw("calc-dir");
    const int         n_occ      = std::stoi(parser.value("n-occ"));
    // Multi-root: --roots=0,1,2 --omegas=w0,w1,w2 build ONE N-root seed
    // bundle (what a production --es-roots=N solve warm-starts from).
    // --root/--omega remain as the single-root spelling.
    std::vector<int>    root_list;
    std::vector<double> omega_list;
    {
        std::stringstream rs(parser.key_exists("roots") ? parser.value("roots")
                                                        : parser.value("root"));
        std::stringstream os(parser.key_exists("omegas")
                                 ? parser.value("omegas")
                                 : parser.value("omega"));
        std::string t;
        while (std::getline(rs, t, ','))
            if (!t.empty()) root_list.push_back(std::stoi(t));
        while (std::getline(os, t, ','))
            if (!t.empty()) omega_list.push_back(std::stod(t));
    }
    if (root_list.empty() || root_list.size() != omega_list.size()) {
        if (world.rank() == 0)
            print("ERROR: --roots and --omegas must be same-length CSV lists.");
        finalize();
        return 2;
    }
    const std::string archive_path =
        parser.key_exists("archive") ? parser.value_raw("archive") : "";
    const double      scale      =
        parser.key_exists("scale") ? std::stod(parser.value("scale"))
                                   : std::sqrt(2.0);

    // Discretization: --thresh wins (k from the standard table); otherwise,
    // with --archive, k comes from the archive header (reliable) and thresh
    // from the inverse table; without either, the legacy default 1e-4.
    double thresh =
        parser.key_exists("thresh") ? std::stod(parser.value("thresh")) : 1e-4;
    int    k      = default_k_for_thresh(thresh);

    {
        // Read RSPVEC
        auto rspvec_read = read_rspvec(rspvec_path);
        const auto &entries = rspvec_read.second;   // plain names: a lambda may not capture a structured binding (C++17)
        for (int ri : root_list) {
            if (ri < 0 || ri >= static_cast<int>(entries.size())) {
                if (world.rank() == 0)
                    print("ERROR: root index", ri, "out of range (file has",
                          static_cast<int>(entries.size()), "entries).");
                finalize();
                return 2;
            }
        }

        // Read Molden
        DaltonMoldenResult molden = read_molden(molden_path);
        const int n_mo  = molden.n_mo;
        const int n_ao  = molden.n_ao;
        const int n_vir = n_mo - n_occ;

        if (n_vir <= 0) {
            if (world.rank() == 0)
                print("ERROR: n_occ=", n_occ, ">= n_mo=", n_mo);
            finalize();
            return 2;
        }

        using vecfuncT = std::vector<madness::real_function_3d>;

        // Set MADNESS cell and discretization parameters BEFORE projecting.
        // With --archive: box L (and default k) from the archive header, the
        // protocol funnel the solvers use; the ground state then loads at the
        // active (k, thresh) and the occupied-gauge rotation U is built.
        std::optional<GroundState> gs;
        Tensor<double> U;   // occupied-gauge rotation (0-size = none)
        if (!archive_path.empty()) {
            auto header = GroundState::read_archive_header(world, archive_path);
            if (static_cast<int>(header.nmo_alpha) != n_occ) {
                // checked on the HEADER (before any Function is created, so
                // the early finalize()+return stays teardown-safe)
                if (world.rank() == 0)
                    print("ERROR: --n-occ =", n_occ, "but the archive has",
                          static_cast<int>(header.nmo_alpha),
                          "occupied alpha orbitals.");
                finalize();
                return 2;
            }
            if (!parser.key_exists("thresh")) {
                k      = header.k;                    // archive k is reliable
                thresh = default_thresh_for_k(k);     // header thresh is not
            }
            set_response_protocol(world, header.L, thresh, k);

            // Molecule from the calc_info json next to the archive (same
            // pattern as tests/test_calc_manager_run.cpp).
            Molecule molecule;
            auto archive_dir =
                std::filesystem::path(archive_path).parent_path();
            for (const auto& name :
                 {"moldft.calc_info.json", "mad.calc_info.json"}) {
                auto candidate = archive_dir / name;
                if (std::filesystem::exists(candidate)) {
                    std::ifstream ifs(candidate);
                    nlohmann::json j;
                    ifs >> j;
                    nlohmann::json mol_json;
                    if (j.contains("tasks") && j["tasks"].is_array() &&
                        !j["tasks"].empty())
                        mol_json = j["tasks"][0]["molecule"];
                    else if (j.contains("molecule"))
                        mol_json = j["molecule"];
                    if (!mol_json.is_null()) molecule.from_json(mol_json);
                    break;
                }
            }
            gs.emplace(GroundState::from_archive(world, archive_path, molecule));

            std::string fock_json;
            for (const auto& name : {"moldft.fock.json", "mad.fock.json"}) {
                auto candidate = archive_dir / name;
                if (std::filesystem::exists(candidate)) {
                    fock_json = candidate.string();
                    break;
                }
            }
            auto coulop = poperatorT(
                CoulombOperatorPtr(world, gs->params().lo(), 0.001 * thresh));
            gs->prepare(world, 0.001 * thresh, coulop, fock_json);

            // Project the DALTON occupied MOs to MRA and build
            // U = M (M^T M)^{-1/2} (fidelity print + <0.5 hard error inside).
            vecfuncT phi_dal;
            for (int i = 0; i < n_occ; ++i) {
                std::vector<double> w(
                    molden.mo_coeffs.begin() + static_cast<ptrdiff_t>(i) * n_ao,
                    molden.mo_coeffs.begin() +
                        static_cast<ptrdiff_t>(i + 1) * n_ao);
                phi_dal.push_back(project_dalton_weights(
                    world, molden.basis, std::move(w), thresh));
            }
            truncate(world, phi_dal, thresh);
            U = occupied_gauge_rotation(world, phi_dal, gs->orbitals_alpha(),
                                        /*verbose=*/true);
        } else {
            // Legacy path: hardcoded box, NO gauge rotation. Only valid when
            // the MADNESS ground state is CANONICAL (moldft `localize canon`).
            Tensor<double> cell(3L, 2L);
            for (int i = 0; i < 3; i++) { cell(i, 0) = -200.0; cell(i, 1) = 200.0; }
            FunctionDefaults<3>::set_cell(cell);
            FunctionDefaults<3>::set_k(k);
            FunctionDefaults<3>::set_thresh(thresh);
            if (world.rank() == 0) {
                print("");
                print("**********************************************************************");
                print("*** WARNING: no --archive given — writing an UNROTATED seed.      ***");
                print("*** DALTON pairs its response vectors with its CANONICAL occupied ***");
                print("*** MOs; response vectors transform covariantly with the occupied ***");
                print("*** orbitals. This seed is therefore only valid if the MADNESS    ***");
                print("*** ground state is canonical too (moldft `localize canon`).      ***");
                print("*** Against a LOCALIZED ground state (the usual `localize new`)   ***");
                print("*** the per-orbital pairing is scrambled and the seed is SILENTLY ***");
                print("*** WORTHLESS — measured on the h2o FD A/B: seed-implied alpha_zz ***");
                print("*** 3.56 vs DALTON's 8.39, ZERO iteration savings.                ***");
                print("*** Pass --archive=<moldft restartdata> to rotate the seed into   ***");
                print("*** the MADNESS occupied gauge (U = M (M^T M)^{-1/2}).            ***");
                print("**********************************************************************");
                print("");
            }
        }

        const std::string key    = protocol_key(thresh, k);
        const std::string bundle = calc_dir + "/es__" + key;

        // --full: write a Full (X,Y) bundle the RPA/--es-full solve resumes
        // from (default: TDA, X-only). DALTON stores the RPA de-excitation
        // block as Y_ai; our y_alpha follows the eigenvector dictionary
        // y_i = -sum_a Y_ai phi_a, so the projected Y block is NEGATED by
        // default. --yflip toggles that sign for the convention check (the RPA
        // metric ||X||^2-||Y||^2 is sign-independent; the discriminator is
        // whether the seeded Full solve converges to the intended root).
        const bool full   = parser.key_exists("full");
        const double ysign = parser.key_exists("yflip") ? +1.0 : -1.0;

        // Project one occ-vir block (flat, row-major occ-outer) into n_occ
        // Functions scaled by sgn*scale — shared machinery (dalton_mra.hpp).
        auto project_block = [&](const std::vector<double>& blk,
                                 double sgn) -> vecfuncT {
            return project_dalton_ov_block(world, molden.basis,
                                           molden.mo_coeffs, n_ao, n_mo, n_occ,
                                           n_vir, blk, thresh, sgn * scale);
        };

        const size_t NR = root_list.size();
        std::vector<vecfuncT> all_x(NR), all_y(NR);
        std::vector<double> xnorms(NR, 0.0), ynorms(NR, 0.0);
        for (size_t rr = 0; rr < NR; ++rr) {
            const auto& entry = entries[static_cast<size_t>(root_list[rr])];
            auto [X_flat, Y_flat] = split_ov(entry.vec, n_occ, n_vir);
            for (double v : X_flat) xnorms[rr] += v * v;
            for (double v : Y_flat) ynorms[rr] += v * v;
            all_x[rr] = project_block(X_flat, +1.0);
            if (full) {
                if (!Y_flat.empty()) {
                    all_y[rr] = project_block(Y_flat, ysign);
                } else {   // CIS/TDA file (no Y block): promote with Y = 0
                    all_y[rr] = madness::copy(world, all_x[rr]);
                    madness::scale(world, all_y[rr], 0.0);
                }
            }
            // Rotate into the MADNESS occupied gauge (--archive): x and y
            // rotate with the SAME U, before any downstream Q-projection.
            if (U.size() > 0) {
                all_x[rr] = transform(world, all_x[rr], U);
                truncate(world, all_x[rr], thresh);
                if (full) {
                    all_y[rr] = transform(world, all_y[rr], U);
                    truncate(world, all_y[rr], thresh);
                }
            }
        }

        // Common State fields (identical shape for TDA and Full; only the root
        // element type differs: ResponseStateX vs ResponseStateXY).
        auto set_common = [&](auto& s) {
            s.omega = Tensor<double>(static_cast<long>(NR));
            for (size_t rr = 0; rr < NR; ++rr)
                s.omega(static_cast<long>(rr)) = omega_list[rr];
            s.iter = 0;
            s.last_bsh_residual     = std::vector<double>(NR, 1.0);
            s.last_density_residual = std::vector<double>(NR, 1.0);
            s.last_omega_residual   = std::vector<double>(NR, 1.0);
        };
        if (full) {
            ESSolver<Full, ClosedShell>::State s;
            for (size_t rr = 0; rr < NR; ++rr) {
                ResponseStateXY<ClosedShell> root;
                root.x_alpha = all_x[rr];
                root.y_alpha = all_y[rr];
                s.roots.push_back(std::move(root));
            }
            set_common(s);
            save_es_roots<Full, ClosedShell>(world, s, bundle, /*converged=*/false);
        } else {
            ESSolver<TDA, ClosedShell>::State s;
            for (size_t rr = 0; rr < NR; ++rr) {
                ResponseStateX<ClosedShell> root;
                root.x_alpha = all_x[rr];
                s.roots.push_back(std::move(root));
            }
            set_common(s);
            save_es_roots<TDA, ClosedShell>(world, s, bundle, /*converged=*/false);
        }

        if (world.rank() == 0) {
            print("seed_from_dalton: wrote", static_cast<int>(NR),
                  "root(s) x", n_occ, "orbital functions   type =",
                  full ? "full (X,Y)" : "tda (X)");
            print("  rspvec    =", rspvec_path);
            print("  molden    =", molden_path);
            for (size_t rr = 0; rr < NR; ++rr) {
                if (full)
                    print("  root", root_list[rr], ": omega(au) =",
                          omega_list[rr], "  ||X||^2 =", xnorms[rr],
                          "  ||Y||^2 =", ynorms[rr], "  ||X||^2-||Y||^2 =",
                          xnorms[rr] - ynorms[rr],
                          "  (RPA metric ~0.5)   ysign =", ysign,
                          "  scale =", scale);
                else
                    print("  root", root_list[rr], ": omega(au) =",
                          omega_list[rr], "  ||X||^2 =", xnorms[rr],
                          "  scale =", scale);
            }
            print("  n_occ =", n_occ, "  n_vir =", n_vir,
                  "  n_ao =", n_ao, "  n_mo =", n_mo);
            print("  gauge =", U.size() > 0
                      ? "ROTATED into the MADNESS occupied gauge (--archive)"
                      : "UNROTATED (no --archive; canonical ground state only)");
            print("  protocol_key =", key, "  k =", k, "  thresh =", thresh);
            print("  bundle    =", bundle);
            print("  Run (NP=1): test_calc_manager_run --archive=<gs> --calc-dir=",
                  calc_dir, " --es-roots=", static_cast<int>(NR),
                  " --protocol=", thresh, " -> ES node Resumes from seed.");
        }
    }  // all Function/State objects destruct here, before finalize()

    finalize();
    return 0;
}
