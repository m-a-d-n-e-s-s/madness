#ifndef MOLRESPONSE_V3_SOLVERS_ES_ANALYSIS_HPP
#define MOLRESPONSE_V3_SOLVERS_ES_ANALYSIS_HPP

// =========================================================================
// Post-convergence excited-state analysis (closed-shell).
//
// Ports the legacy TDDFT::analysis + TDDFT::analyze_vectors
// (src/apps/molresponse_legacy/iterate_excited.cc, TDDFT.cc) onto the v3
// ESSolver<Type, ClosedShell>::State. For each converged root it computes
// the transition properties expected of an excited-state report:
//
//   * transition dipole moments      d_a   = Σ_j <φ_j | r_a · (x_j [+ y_j])>
//   * dipole oscillator strength      f     = (2/3) |d|² ω
//   * transition quadrupole moments   Q_ab  = Σ_j <φ_j | r_a r_b · (x_j[+y_j])>
//   * dominant occupied contributions       (per-occupied |x_j| weights)
//
// and, per response-orbital component, the legacy analyze_vectors output:
// AO (sto-3g) Mulliken-style population analysis, orbital center, radius.
//
// Type selects the transition-moment assembly: TDA uses X only; Full adds
// the Y block (matching the legacy full-response sign convention — Y enters
// the moments with the same sign as X). Open-shell ES is out of scope, so
// these helpers are gated on ClosedShell.
//
// Collective discipline: every inner-product / norm here is collective on
// `world` and runs on all ranks; only the human-readable printing is
// rank-0-gated. Output is print-only for now — a dedicated chemistry
// property layer (JSON / HDF5 / cube) will consume the returned structs
// later, which is why analyze_es_bundle returns data rather than printing
// it directly.
// =========================================================================

#include "../GroundState.hpp"
#include "../kernels/tags.hpp"   // PrintLevel
#include "es_save_load.hpp"      // try_load_es_bundle (analyze-only path)
#include "es_solver.hpp"

#include <madness/chem/SCF.h>            // MomentFunctor, DipoleFunctor, SCF::project_ao_basis_only
#include <madness/chem/molecularbasis.h> // AtomicBasisSet
#include <madness/mra/mra.h>

#include <array>
#include <cmath>
#include <cstdlib>
#include <fstream>
#include <string>
#include <type_traits>
#include <vector>

namespace molresponse_v3 {

/// Per-root transition properties (closed-shell). Mirrors the data the
/// legacy TDDFT::analysis printed for one excited-state root.
struct ESRootAnalysis {
  double                                omega = 0.0;
  std::array<double, 3>                 dipole{{0.0, 0.0, 0.0}};
  double                                oscillator = 0.0;
  std::array<std::array<double, 3>, 3>  quadrupole{};
  /// |x_alpha[j]| per occupied orbital j (used for dominant-contribution
  /// reporting). Y is intentionally not folded in here — the legacy report
  /// listed the X and Y occupied weights side by side, so the Full printer
  /// recomputes the Y weights when it needs them.
  std::vector<double>                   orbital_norms;
};

namespace detail_es_analysis {

/// Cartesian coordinate functions x, y, z (degree-1 moments). Shared by the
/// transition-moment assembly and the per-orbital analyzer.
inline std::array<madness::real_function_3d, 3>
coordinate_functions(madness::World &world) {
  using namespace madness;
  const double tol = FunctionDefaults<3>::get_thresh();
  std::array<real_function_3d, 3> r;
  for (int a = 0; a < 3; ++a) {
    std::vector<int> dir(3, 0);
    dir[a] = 1;
    r[a] = real_factory_3d(world).functor(
        real_functor_3d(new MomentFunctor(dir)));
    r[a].truncate(tol);
  }
  return r;
}

}  // namespace detail_es_analysis

/// Compute transition properties for every root in a converged closed-shell
/// ES bundle. `ground_orbitals` are the occupied alpha orbitals
/// (gs.orbitals_alpha()). Pure computation — collective on all ranks, no
/// printing.
template <typename Type, typename Shell>
std::vector<ESRootAnalysis>
analyze_es_bundle(madness::World &world,
                  const std::vector<madness::real_function_3d> &ground_orbitals,
                  const typename ESSolver<Type, Shell>::State &state) {
  using namespace madness;
  static_assert(std::is_same_v<Shell, ClosedShell>,
                "analyze_es_bundle: closed-shell only");
  constexpr bool has_y = std::is_same_v<Type, Full>;

  const size_t m = state.roots.size();
  const size_t n = ground_orbitals.size();

  // ---- Bisection probe (TEMPORARY) ---------------------------------------
  // MADQC_ANALYZE_PROBE caps how far into the analysis we run, to localize the
  // teardown crash with one rebuild + a few quick runs:
  //   0  omega only (no MADNESS ops)
  //   1  + create coordinate functions (MomentFunctor over the box)
  //   2  + deep-copy + reconstruct the ground orbitals
  //   3  + norm2 + transition-dipole loop (inner/mul)
  //   4+ full (adds the quadrupole loop)            [default]
  int probe = 4;
  if (const char *p = std::getenv("MADQC_ANALYZE_PROBE")) probe = std::atoi(p);
  if (world.rank() == 0)
    madness::print("[ANALYZE] probe stage =", probe);

  std::vector<ESRootAnalysis> out(m);
  for (size_t i = 0; i < m; ++i)
    out[i].omega =
        (static_cast<long>(i) < state.omega.size()) ? state.omega(i) : 0.0;
  if (probe <= 0) return out;

  const auto r = detail_es_analysis::coordinate_functions(world);
  if (probe <= 1) return out;

  // Operate on a DEEP COPY of the ground orbitals — never gs's live handles.
  vecfuncT orbs = madness::copy(world, ground_orbitals);
  madness::reconstruct(world, orbs);
  if (probe <= 2) return out;

  // Stage 3: norms + transition dipole.
  for (size_t i = 0; i < m; ++i) {
    ESRootAnalysis &a = out[i];
    const auto &X = state.roots[i].x_alpha;

    a.orbital_norms.assign(n, 0.0);
    for (size_t j = 0; j < n; ++j) a.orbital_norms[j] = X[j].norm2();

    for (size_t j = 0; j < n; ++j) {
      for (int ax = 0; ax < 3; ++ax)
        a.dipole[ax] += inner(orbs[j], r[ax] * X[j]);
      if constexpr (has_y) {
        const auto &Y = state.roots[i].y_alpha;
        for (int ax = 0; ax < 3; ++ax)
          a.dipole[ax] += inner(orbs[j], r[ax] * Y[j]);
      }
    }
    for (auto &d : a.dipole) d *= -std::sqrt(2.0);

    a.oscillator = 2.0 / 3.0 *
                   (a.dipole[0] * a.dipole[0] + a.dipole[1] * a.dipole[1] +
                    a.dipole[2] * a.dipole[2]) *
                   a.omega;
  }
  if (probe <= 3) return out;

  // Stage 4: transition quadrupole.
  for (size_t i = 0; i < m; ++i) {
    ESRootAnalysis &a = out[i];
    const auto &X = state.roots[i].x_alpha;
    for (size_t j = 0; j < n; ++j) {
      for (int aa = 0; aa < 3; ++aa)
        for (int bb = 0; bb < 3; ++bb) {
          a.quadrupole[aa][bb] +=
              inner(orbs[j], r[aa] * r[bb] * X[j]);
          if constexpr (has_y) {
            const auto &Y = state.roots[i].y_alpha;
            a.quadrupole[aa][bb] +=
                inner(orbs[j], r[aa] * r[bb] * Y[j]);
          }
        }
    }
    for (auto &row : a.quadrupole)
      for (auto &q : row) q *= std::sqrt(2.0);
  }
  return out;
}

/// Print the transition-property block for a bundle (rank 0 only). `tda`
/// controls only the dominant-contribution layout (X-only vs X/Y columns);
/// when `tda` is false, `y_orbital_norms[i]` supplies the per-root Y weights.
inline void
print_es_analysis(madness::World &world,
                  const std::vector<ESRootAnalysis> &an, bool tda,
                  const std::vector<std::vector<double>> &y_orbital_norms = {}) {
  using namespace madness;
  if (world.rank() != 0) return;

  for (size_t i = 0; i < an.size(); ++i) {
    const ESRootAnalysis &a = an[i];
    const size_t n = a.orbital_norms.size();

    printf("   Response Function %d\t\t%7.8f a.u.", static_cast<int>(i), a.omega);
    print("\n   --------------------------------------------");
    printf("   Response Function %d\t\t%7.8f eV", static_cast<int>(i),
           a.omega * 27.2114);
    print("\n   --------------------------------------------");

    print("\n   Transition Dipole Moments");
    printf("   X: %7.8f   Y: %7.8f   Z: %7.8f\n", a.dipole[0], a.dipole[1],
           a.dipole[2]);

    printf("\n   Dipole Oscillator Strength: %7.8f\n", a.oscillator);

    print("\n   Transition Quadrupole Moments");
    printf("   %16s %16s %16s\n", "X", "Y", "Z");
    printf("   X %16.8f %16.8f %16.8f\n", a.quadrupole[0][0], a.quadrupole[0][1],
           a.quadrupole[0][2]);
    printf("   Y %16.8f %16.8f %16.8f\n", a.quadrupole[1][0], a.quadrupole[1][1],
           a.quadrupole[1][2]);
    printf("   Z %16.8f %16.8f %16.8f\n", a.quadrupole[2][0], a.quadrupole[2][1],
           a.quadrupole[2][2]);

    // Dominant contributions: top 5 occupied orbitals by |x_j| (legacy order).
    std::vector<size_t> order(n);
    for (size_t j = 0; j < n; ++j) order[j] = j;
    std::sort(order.begin(), order.end(), [&](size_t p, size_t q) {
      return a.orbital_norms[p] > a.orbital_norms[q];
    });

    if (tda) {
      print("\n   Dominant Contributions:");
      for (size_t j = 0; j < std::min(size_t(5), n); ++j)
        printf("   Occupied %d   %7.8f\n", static_cast<int>(order[j]),
               a.orbital_norms[order[j]]);
      print("\n");
    } else {
      const std::vector<double> &yn =
          (i < y_orbital_norms.size()) ? y_orbital_norms[i]
                                       : std::vector<double>(n, 0.0);
      print("\n   Dominant Contributions:");
      print("                  x          y");
      for (size_t j = 0; j < std::min(size_t(5), n); ++j) {
        const size_t o = order[j];
        printf("   Occupied %d   %7.8f %7.8f\n", static_cast<int>(o),
               a.orbital_norms[o], (o < yn.size()) ? yn[o] : 0.0);
      }
      print("\n");
    }
  }
}

/// Per-response-orbital AO-population analysis — port of TDDFT::analyze_vectors.
/// Projects an sto-3g AO basis, prints each component's Mulliken-style
/// population (sto3g.print_anal), dipole center, and radius. Collective on
/// `world`; printing gated on rank 0 + print_level.
inline void
analyze_response_orbitals(madness::World &world, const Molecule &molecule,
                          const std::vector<madness::real_function_3d> &x,
                          const std::string &label, PrintLevel print_level) {
  using namespace madness;
  const double vtol = FunctionDefaults<3>::get_thresh();

  AtomicBasisSet sto3g("sto-3g");
  vecfuncT ao = SCF::project_ao_basis_only(world, sto3g, molecule);
  Tensor<double> C = matrix_inner(world, ao, x);

  const int nmo = static_cast<int>(x.size());
  const auto r = detail_es_analysis::coordinate_functions(world);
  real_function_3d frsq = (r[0] * r[0] + r[1] * r[1] + r[2] * r[2]);
  frsq.truncate(vtol);

  // <x_i | r^2 | x_i> then subtract <x_i | r_a | x_i>^2 to get the radius.
  Tensor<double> rsq = inner(world, x, mul_sparse(world, frsq, x, vtol));
  Tensor<double> dip(3, nmo);
  for (int axis = 0; axis < 3; ++axis) {
    dip(axis, _) = inner(world, x, mul_sparse(world, r[axis], x, vtol));
    for (int i = 0; i < nmo; ++i) rsq(i) -= dip(axis, i) * dip(axis, i);
  }

  // Function::size() is COLLECTIVE (global coefficient-count reduction), so it
  // must be called on every rank — not inside the rank-0 print guard, or the
  // other ranks never match the reduction (MPI_ERR_TRUNCATE). Precompute here.
  std::vector<std::size_t> ncoeff_i(static_cast<size_t>(nmo));
  for (int i = 0; i < nmo; ++i) ncoeff_i[static_cast<size_t>(i)] = x[i].size();

  if (world.rank() == 0 && print_level >= PrintLevel::Normal) {
    size_t ncoeff = 0;
    for (int i = 0; i < nmo; ++i) {
      ncoeff += ncoeff_i[static_cast<size_t>(i)];
      print(label + " orbital : ", i);
      printf("ncoeff=%.2e:",
             static_cast<double>(ncoeff_i[static_cast<size_t>(i)]));
      printf("center=(%.2f,%.2f,%.2f) : radius=%.2f\n", dip(0, i), dip(1, i),
             dip(2, i), std::sqrt(std::abs(rsq(i))));
      sto3g.print_anal(molecule, C(i, _));
      printf("total number of coefficients = %.8e\n\n",
             static_cast<double>(ncoeff));
    }
  }
}

/// Top-level driver: run the full legacy analysis report on a converged
/// Serialize a bundle's transition properties (omega, transition dipole,
/// oscillator strength, transition quadrupole, dominant occupied weights) to
/// JSON. This is the persistent, machine-readable record for later analysis;
/// the human-readable console report is print_es_analysis. `stable_index`
/// (if non-empty) attaches each root's stable identity / root_id.
inline nlohmann::json
es_analysis_to_json(const std::vector<ESRootAnalysis> &an, bool tda,
                    const std::vector<std::vector<double>> &y_norms,
                    const std::vector<int> &stable_index) {
  nlohmann::json roots = nlohmann::json::array();
  for (size_t i = 0; i < an.size(); ++i) {
    const ESRootAnalysis &a = an[i];
    nlohmann::json r;
    r["slot"]    = static_cast<int>(i);
    if (i < stable_index.size()) {
      r["stable_index"] = stable_index[i];
      r["root_id"]      = make_root_id(stable_index[i]);
    }
    r["omega"]                 = a.omega;
    r["omega_ev"]              = a.omega * 27.2114;
    r["transition_dipole"]     = {a.dipole[0], a.dipole[1], a.dipole[2]};
    r["oscillator_strength"]   = a.oscillator;
    r["transition_quadrupole"] = {
        {a.quadrupole[0][0], a.quadrupole[0][1], a.quadrupole[0][2]},
        {a.quadrupole[1][0], a.quadrupole[1][1], a.quadrupole[1][2]},
        {a.quadrupole[2][0], a.quadrupole[2][1], a.quadrupole[2][2]}};
    r["occupied_weights_x"]    = a.orbital_norms;
    if (!tda && i < y_norms.size()) r["occupied_weights_y"] = y_norms[i];
    roots.push_back(r);
  }
  nlohmann::json j;
  j["analysis"]      = tda ? "tda" : "full";
  j["n_roots"]       = static_cast<int>(an.size());
  j["protocol_key"]  = protocol_key();
  j["units"]         = {{"omega", "hartree"}, {"omega_ev", "eV"},
                        {"transition_dipole", "a.u."},
                        {"transition_quadrupole", "a.u."}};
  j["roots"]         = roots;
  return j;
}

/// closed-shell ES bundle and print it. Computes transition properties,
/// prints them, then runs the per-orbital AO analysis for each root's X
/// (and Y, for Full) component — exactly the sequence at the tail of the
/// legacy TDDFT::iterate_excited. When `json_out` is non-empty, the
/// transition properties are also written there (rank 0) for later analysis.
/// `do_orbital_analysis` toggles the per-orbital AO-population pass
/// (analyze_response_orbitals); OFF by default — it is the heavy AO-basis path
/// (fresh AtomicBasisSet + project_ao_basis_only + matrix_inner + mul_sparse,
/// per root) and the current suspect for the teardown crash. The transition
/// properties + JSON (the science output) always run.
template <typename Type, typename Shell>
void report_es_analysis(madness::World &world, const GroundState &gs,
                        const typename ESSolver<Type, Shell>::State &state,
                        PrintLevel print_level,
                        const std::string &json_out = "",
                        bool do_orbital_analysis = false) {
  using namespace madness;
  static_assert(std::is_same_v<Shell, ClosedShell>,
                "report_es_analysis: closed-shell only");
  constexpr bool has_y = std::is_same_v<Type, Full>;

  // Bisection probe (TEMPORARY) — shared with analyze_es_bundle.
  //   <0 : report is a pure no-op (just a fence) — should equal --es-load-only
  //   >=0: run, but analyze_es_bundle internally caps at the stage value
  int probe = 4;
  if (const char *p = std::getenv("MADQC_ANALYZE_PROBE")) probe = std::atoi(p);
  if (probe < 0) {
    if (world.rank() == 0)
      madness::print("[ANALYZE] probe < 0 — report is a no-op (fence only)");
    world.gop.fence();
    return;
  }

  const auto &ground_orbitals = gs.orbitals_alpha();
  auto an = analyze_es_bundle<Type, Shell>(world, ground_orbitals, state);

  // Y occupied weights for the dominant-contribution table (Full only).
  // Gated to stage >=3 (same point X norms/dipoles compute) so the probe can
  // isolate it — it is the only collective in report's tail.
  std::vector<std::vector<double>> y_norms;
  if constexpr (has_y) {
    if (probe >= 3) {
      y_norms.resize(state.roots.size());
      for (size_t i = 0; i < state.roots.size(); ++i) {
        const auto &Y = state.roots[i].y_alpha;
        y_norms[i].assign(Y.size(), 0.0);
        for (size_t j = 0; j < Y.size(); ++j) y_norms[i][j] = Y[j].norm2();
      }
    }
  }

  print_es_analysis(world, an, /*tda=*/!has_y, y_norms);

  // Persist the transition properties to a dedicated JSON for later analysis.
  if (world.rank() == 0 && !json_out.empty()) {
    auto j = es_analysis_to_json(an, /*tda=*/!has_y, y_norms,
                                 state.stable_index);
    std::ofstream out(json_out);
    out << j.dump(2) << "\n";
    print("[ANALYZE] wrote transition-property JSON ->", json_out);
  }

  if (do_orbital_analysis) {
    if (world.rank() == 0 && print_level >= PrintLevel::Normal)
      print("--------------------------------------------------------");
    for (size_t i = 0; i < state.roots.size(); ++i) {
      analyze_response_orbitals(world, gs.molecule(), state.roots[i].x_alpha,
                                "x_" + std::to_string(i) + "_", print_level);
      if (world.rank() == 0 && print_level >= PrintLevel::Normal)
        print("--------------------------------------------------------");
    }
    if constexpr (has_y) {
      for (size_t i = 0; i < state.roots.size(); ++i) {
        analyze_response_orbitals(world, gs.molecule(), state.roots[i].y_alpha,
                                  "y_" + std::to_string(i) + "_", print_level);
        if (world.rank() == 0 && print_level >= PrintLevel::Normal)
          print("--------------------------------------------------------");
      }
    }
  }

  // Flush all analysis tasks (AO projection, mul_sparse, matrix_inner) and
  // synchronize before returning, so no in-flight op holds a temporary
  // FunctionImpl past this scope — otherwise those temporaries (and the
  // ground-state functions) can hit the World's deferred-cleanup machinery in
  // a torn-down state at program teardown.
  world.gop.fence();
}

/// Analyze-only: load a converged ES bundle from `calc_dir` (via the same
/// restart precedence the solver uses), reproject its roots to the active
/// k/thresh, and print the full transition-property report — no solve. The
/// caller must have prepared `gs` at the protocol it wants the report at
/// (set_response_protocol + gs.prepare). Returns 0 on success, 1 if no
/// compatible bundle was found. Useful both as a standalone reporting path
/// and as an isolated test of the analysis collectives (no solver involved).
template <typename Type, typename Shell>
int analyze_es_bundle_from_disk(madness::World &world, const GroundState &gs,
                                const std::string &calc_dir,
                                PrintLevel print_level,
                                bool do_orbital_analysis = false,
                                bool load_only = false) {
  using namespace madness;
  auto loaded = try_load_es_bundle<Type, Shell>(world, calc_dir);
  if (!loaded) {
    if (world.rank() == 0)
      print("[ANALYZE] no compatible ES bundle in", calc_dir,
            "— nothing to analyze.");
    return 1;
  }
  // try_load may return a coarser bundle; reproject to the active k/thresh so
  // the roots and the prepared ground orbitals share a basis (same path the
  // solver's per-protocol prepare hook takes).
  const int    k  = FunctionDefaults<3>::get_k();
  const double tt = FunctionDefaults<3>::get_thresh();
  for (auto &root : loaded->state.roots)
    for (auto *blk : root.blocks())
      for (auto &fn : *blk) fn = madness::project(fn, k, tt);

  if (load_only) {
    // Bisection: load + reproject only, no report.
    if (world.rank() == 0)
      print("[ANALYZE] load-only: loaded", static_cast<int>(loaded->state.roots.size()),
            "roots from", loaded->bundle_dir, "— skipping report.");
    // Destroy the archive-loaded functions, THEN fence — so their deferred
    // WorldContainer cleanup drains while the World is still alive (a
    // load-and-exit run never fences otherwise; the queue would be processed at
    // teardown, after the World's DeferredCleanup is gone -> dead-mutex crash).
    loaded.reset();
    world.gop.fence();
    return 0;
  }

  if (world.rank() == 0)
    print("[ANALYZE] reporting", static_cast<int>(loaded->state.roots.size()),
          "roots from", loaded->bundle_dir, "(source protocol",
          loaded->source_protocol_key, ", exact=", loaded->exact, ")");
  const std::string json_out =
      calc_dir + "/es_analysis__" + protocol_key() + ".json";
  report_es_analysis<Type, Shell>(world, gs, loaded->state, print_level,
                                  json_out, do_orbital_analysis);
  // Drain the archive-loaded functions' deferred cleanup while the World is
  // alive (see note above) before returning into teardown.
  loaded.reset();
  world.gop.fence();
  return 0;
}

}  // namespace molresponse_v3

#endif  // MOLRESPONSE_V3_SOLVERS_ES_ANALYSIS_HPP
