#ifndef MOLRESPONSE_V3_RESPONSE_PROPERTY_PLANNER_HPP
#define MOLRESPONSE_V3_RESPONSE_PROPERTY_PLANNER_HPP

// ===========================================================================
// ResponsePropertyPlanner — the front of the RESPONSE-PROPERTY layer (Tier A).
//
// Two-tier design (see docs):
//   Tier A — Response properties: Cartesian tensors contracted directly from
//            response/excited states (α, β, ∂α/∂Q, ...). Specifying these IS
//            the calculation — they determine which states to solve. This
//            planner turns a ResponsePropertyRequest into the deduplicated
//            list of FD/ES states a calc manager must compute.
//   Tier B — Chemical properties (separate module): take the Tier-A tensors
//            (+ aux data like normal modes / masses) and transform them into
//            chemist-facing observables with units (Raman spectra, [α]_D,
//            σ_TPA, β_HRS). Post-processing + printing; never triggers a new
//            electronic solve.
//
// This header is Tier A only. Pure planning: no World, no I/O, no dispatch.
//
// Response properties supported:
//
//   Polarizability          α(ω) — FD (dipole_a, ω) per axis a, per ω.
//   Hyperpolarizability     β(ω₁,ω₂,ω₃) — `BetaProcess` picks the triplet
//                           from each driver ω; unique |ω| → FD per axis.
//                             Static : (0,0,0)
//                             SHG    : (ω, ω, -2ω) → {ω, 2ω}
//                             OR     : (-ω, ω, 0)  → {ω, 0}
//                             EOPE   : (-ω, ω, 0)  → {ω, 0}
//   PolarizabilityGradient  ∂α/∂Q — the Raman tensor. `GradientMode`:
//                             Nuclear  (vibrational): dipole FD at the
//                               optical ω(s) + nuclear-displacement FD at 0
//                               for all atoms (quadratic / VBC). Feeds the
//                               Tier-B vibrational Raman spectrum.
//                             Resonant (resonance Raman): ES bundle + dipole
//                               FD at the excitation energies. Feeds the
//                               Tier-B resonance-Raman profile.
//
// Symbolic requests resolved by the calc manager post-solve:
//   - derived_fd[].es_root_id == "*"  → one per converged ES root.
//   - nuclear_fd[].pert.atom  <  0    → one per atom in the molecule.
// Keeping these symbolic lets the planner stay pure and molecule-independent.
// ===========================================================================

#include "Perturbations.hpp"            // Perturbation

#include <algorithm>
#include <cstdio>
#include <functional>
#include <set>
#include <string>
#include <utility>
#include <vector>

namespace molresponse_v3 {

enum class ResponsePropertyKind { Polarizability, Hyperpolarizability,
                                  PolarizabilityGradient };

/// β-process frequency convention (driver ω → triplet).
enum class BetaProcess { Static, SHG, OR, EOPE };

/// How the polarizability derivative ∂α/∂Q is obtained.
enum class GradientMode {
  Nuclear,    // ∂α/∂(nuclear Cartesian) — vibrational Raman
  Resonant,   // transition polarizability via excited states — resonance Raman
};

/// User-facing request — one response property at one set of frequencies.
struct ResponsePropertyRequest {
  ResponsePropertyKind kind;
  /// α/β: driver ω(s). For α each ω is the FD freq; for β the triplet is
  /// built per ω. PolarizabilityGradient: the optical ω(s) at which the
  /// (transition) polarizability is evaluated.
  std::vector<double>  frequencies;
  /// Cartesian directions for the OPTICAL (dipole) tensor indices.
  std::vector<char>    axes         = {'x','y','z'};
  /// Raman (PolarizabilityGradient/Nuclear): restrict to a SINGLE nuclear
  /// component (atom index, axis 0=x/1=y/2=z) for a cheap one-component test.
  /// -1 (default) = full per-atom xyz set (full-tensor VBC emission + sentinel
  /// expansion is deferred until the state-parallel design lands).
  int                  raman_nuc_atom = -1;
  int                  raman_nuc_axis = -1;
  BetaProcess          beta_process = BetaProcess::SHG;       // Hyperpolarizability
  GradientMode         gradient_mode = GradientMode::Nuclear;  // PolarizabilityGradient
  int                  n_roots      = 0;                       // Resonant gradient
  std::vector<double>  protocol_thresholds;                    // ramp; all states get this
};

/// One concrete FD state the calc manager must solve.
struct FDRequest {
  Perturbation        pert;
  double              freq;
  std::vector<double> protocols;
};

/// One ES bundle the calc manager must solve.
struct ESRequest {
  bool                tda = true;   // false → Full/RPA
  int                 n_roots;
  std::vector<double> protocols;
};

/// Symbolic FD request — `es_root_id` references a future ES output.
/// "*" expands to one record per converged root post-ES.
struct DerivedFDRequest {
  Perturbation        pert;
  std::string         es_root_id;   // "es_root_0001" or "*"
  std::vector<double> protocols;
  /// FD frequency = es_freq_factor * (ES root energy). 0.5 = two-photon
  /// resonance (2ω = ωₙ), which keeps the FD OFF the linear-response pole at
  /// ω = ωₙ; a factor of 1.0 would sit on the pole and never converge in the
  /// undamped solver.
  double              es_freq_factor = 0.5;
};

/// One VBC quadratic source to build: the (B,C) perturbation pair at
/// (freq_b, freq_c). Calc-manager prerequisites: FD(B,freq_b), FD(C,freq_c).
struct VBCRequest {
  Perturbation        pert_b;
  Perturbation        pert_c;
  double              freq_b = 0.0;
  double              freq_c = 0.0;
  std::vector<double> protocols;
};

struct ResponsePlan {
  std::vector<FDRequest>        fd;          // concrete dipole FD
  std::vector<ESRequest>        es;          // ES bundles
  std::vector<DerivedFDRequest> derived_fd;  // FD-after-ES (freq from a root)
  std::vector<FDRequest>        nuclear_fd;  // nuclear-displacement FD; pert.atom<0
                                             // = all-atoms sentinel (expand per molecule)
  std::vector<VBCRequest>       vbc;         // quadratic sources for beta
};

// ---------------------------------------------------------------------------
namespace detail_planner {

inline int axis_index(char a) {
  switch (a) {
    case 'x': case 'X': return 0;
    case 'y': case 'Y': return 1;
    case 'z': case 'Z': return 2;
    default: return -1;
  }
}

/// Unique |ω| set needed to assemble β at one driver ω.
inline std::vector<double> beta_freqs(BetaProcess proc, double omega) {
  std::set<double> uniq;
  switch (proc) {
    case BetaProcess::Static: uniq.insert(0.0); break;
    case BetaProcess::SHG:    uniq.insert(omega); uniq.insert(2.0 * omega); break;
    case BetaProcess::OR:
    case BetaProcess::EOPE:   uniq.insert(omega); uniq.insert(0.0); break;
  }
  return {uniq.begin(), uniq.end()};
}

/// (omega_B, omega_C) pairs for the VBC quadratic source per beta process.
/// SHG: (w,w); Static: (0,0). OR/EOPE involve signed frequencies / x<->y
/// channel mapping that the beta-iii contraction must resolve, so no VBC is
/// emitted for them yet (the FD states still are).
inline std::vector<std::pair<double, double>>
beta_vbc_pairs(BetaProcess proc, double omega) {
  switch (proc) {
    case BetaProcess::Static: return {{0.0, 0.0}};
    case BetaProcess::SHG:    return {{omega, omega}};
    case BetaProcess::OR:
    case BetaProcess::EOPE:   return {};   // TODO(beta-iii): signed-frequency VBC
  }
  return {};
}

/// Sorted (coarse → fine) deduplicated union of two protocol lists.
inline std::vector<double> union_protocols(std::vector<double> a,
                                            const std::vector<double> &b) {
  for (double t : b) a.push_back(t);
  std::sort(a.begin(), a.end(), std::greater<double>{});
  a.erase(std::unique(a.begin(), a.end()), a.end());
  return a;
}

/// Stable freq key (matches doc-13 archive freq_key precision) for dedupe.
inline std::string fd_freq_key(double f) {
  char buf[16];
  std::snprintf(buf, sizeof buf, "%.5f", f);
  return buf;
}

} // namespace detail_planner

// ---------------------------------------------------------------------------
/// Merge & dedupe across plans. Identical requests collapse with their
/// protocol sets unioned. Dedupe keys: FD/nuclear_fd on
/// (pert.description(), freq_key); ES on (tda, n_roots); derived FD on
/// (pert.description(), es_root_id).
inline ResponsePlan merge_plans(const std::vector<ResponsePlan> &plans) {
  ResponsePlan out;

  auto fd_k = [](const FDRequest &r) {
    return r.pert.description() + "@" + detail_planner::fd_freq_key(r.freq);
  };
  auto es_k = [](const ESRequest &r) {
    return std::string(r.tda ? "tda" : "full") + "_" + std::to_string(r.n_roots);
  };
  auto vbc_k = [](const VBCRequest &r) {
    return r.pert_b.description() + "__" + r.pert_c.description() + "@"
         + detail_planner::fd_freq_key(r.freq_b) + "_"
         + detail_planner::fd_freq_key(r.freq_c);
  };
  auto dfd_k = [](const DerivedFDRequest &r) {
    return r.pert.description() + "@" + r.es_root_id;
  };

  auto upsert_fd = [](std::vector<FDRequest> &list, const FDRequest &r,
                      const std::function<std::string(const FDRequest &)> &key) {
    const auto k = key(r);
    for (auto &e : list) {
      if (key(e) == k) {
        e.protocols = detail_planner::union_protocols(e.protocols, r.protocols);
        return;
      }
    }
    list.push_back(r);
  };

  for (const auto &p : plans) {
    for (const auto &r : p.fd)         upsert_fd(out.fd, r, fd_k);
    for (const auto &r : p.nuclear_fd) upsert_fd(out.nuclear_fd, r, fd_k);
    for (const auto &r : p.es) {
      const auto k = es_k(r);
      bool found = false;
      for (auto &e : out.es) {
        if (es_k(e) == k) {
          e.protocols = detail_planner::union_protocols(e.protocols, r.protocols);
          found = true; break;
        }
      }
      if (!found) out.es.push_back(r);
    }
    for (const auto &r : p.derived_fd) {
      const auto k = dfd_k(r);
      bool found = false;
      for (auto &e : out.derived_fd) {
        if (dfd_k(e) == k) {
          e.protocols = detail_planner::union_protocols(e.protocols, r.protocols);
          found = true; break;
        }
      }
      if (!found) out.derived_fd.push_back(r);
    }
    for (const auto &r : p.vbc) {
      const auto k = vbc_k(r);
      bool found = false;
      for (auto &e : out.vbc) {
        if (vbc_k(e) == k) {
          e.protocols = detail_planner::union_protocols(e.protocols, r.protocols);
          found = true; break;
        }
      }
      if (!found) out.vbc.push_back(r);
    }
  }
  return out;
}

inline ResponsePlan merge_plans(ResponsePlan a, ResponsePlan b) {
  return merge_plans(std::vector<ResponsePlan>{std::move(a), std::move(b)});
}

// ---------------------------------------------------------------------------
/// Plan one ResponsePropertyRequest. Output is self-deduped.
inline ResponsePlan plan_one(const ResponsePropertyRequest &req) {
  ResponsePlan plan;
  using detail_planner::axis_index;
  using detail_planner::beta_freqs;

  auto add_dipole_fd = [&](double w) {
    for (char ax : req.axes) {
      int i = axis_index(ax);
      if (i < 0) continue;
      plan.fd.push_back({Perturbation::dipole(i), w, req.protocol_thresholds});
    }
  };

  switch (req.kind) {
    case ResponsePropertyKind::Polarizability: {
      for (double w : req.frequencies) add_dipole_fd(w);
      break;
    }
    case ResponsePropertyKind::Hyperpolarizability: {
      for (double w : req.frequencies)
        for (double bw : beta_freqs(req.beta_process, w))
          add_dipole_fd(bw);
      // VBC quadratic sources: one per (B,C) axis pair per (omega_B, omega_C).
      using detail_planner::beta_vbc_pairs;
      for (double w : req.frequencies)
        for (auto [fb, fc] : beta_vbc_pairs(req.beta_process, w))
          for (char ax_b : req.axes) {
            int ib = axis_index(ax_b);
            if (ib < 0) continue;
            for (char ax_c : req.axes) {
              int ic = axis_index(ax_c);
              if (ic < 0) continue;
              plan.vbc.push_back({Perturbation::dipole(ib),
                                  Perturbation::dipole(ic), fb, fc,
                                  req.protocol_thresholds});
            }
          }
      break;
    }
    case ResponsePropertyKind::PolarizabilityGradient: {
      if (req.gradient_mode == GradientMode::Nuclear) {
        // Vibrational Raman: β(dipole; dipole@ω, nuclear-disp@0).
        // Dipole states at each optical ω (serve as A and B halves)...
        for (double w : req.frequencies) add_dipole_fd(w);
        // ...plus nuclear-displacement states at ω=0 for ALL atoms. Normal-
        // mode analysis needs the full 3N Cartesian set, so emit all three
        // displacement axes regardless of the optical `axes` subset.
        if (req.raman_nuc_atom >= 0 && req.raman_nuc_axis >= 0) {
          // SINGLE component (cheap test): concrete nuclear pert (no sentinel)
          // + its FD + the dipole x nuclear VBC pairs that feed the contraction.
          const Perturbation nuc =
              Perturbation::nuclear(req.raman_nuc_atom, req.raman_nuc_axis);
          plan.nuclear_fd.push_back({nuc, 0.0, req.protocol_thresholds});
          for (double w : req.frequencies)
            for (char ax_b : req.axes) {
              int ib = axis_index(ax_b);
              if (ib < 0) continue;
              plan.vbc.push_back({Perturbation::dipole(ib), nuc, w, 0.0,
                                  req.protocol_thresholds});
            }
        } else {
          // FULL set: all atoms x xyz nuclear FD. The dipole x nuclear VBC
          // emission + per-atom sentinel expansion in build_dag is DEFERRED
          // until state-parallel lands (use raman_nuc_atom/axis for now).
          for (int dax = 0; dax < 3; ++dax) {
            plan.nuclear_fd.push_back(
                {Perturbation::nuclear_all(dax), 0.0, req.protocol_thresholds});
          }
        }
      } else {
        // Resonance Raman: transition polarizability via excited states.
        plan.es.push_back({/*tda=*/true, req.n_roots, req.protocol_thresholds});
        for (char ax : req.axes) {
          int i = axis_index(ax);
          if (i < 0) continue;
          plan.derived_fd.push_back(
              {Perturbation::dipole(i), /*es_root_id=*/"*",
               req.protocol_thresholds});
        }
      }
      break;
    }
  }
  return merge_plans(std::vector<ResponsePlan>{std::move(plan)});
}

} // namespace molresponse_v3

#endif // MOLRESPONSE_V3_RESPONSE_PROPERTY_PLANNER_HPP
