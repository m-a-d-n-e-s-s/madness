#ifndef MOLRESPONSE_V3_KERNELS_EXCHANGE_CTX_HPP
#define MOLRESPONSE_V3_KERNELS_EXCHANGE_CTX_HPP

// =========================================================================
// Tensor-layer exchange for linear response (doc 28; framework doc 27).
//
// Every exchange is K[a,b](c)_k = Σ_i a_i · g(b,c)_{ik},  g(b,c)_{ik} = Poisson(b_i·c_k).
// We build the shared convolution tensors ONCE, then assemble θ as cheap
// contractions — instead of the op-level compute_V0x/compute_gamma which each
// rebuild their own exchange (so Tx was built twice; no home for the g0 cache).
//
//   g0[i*n+k] = Poisson(φ_i·φ_k)   (Class 1: φ-only, cache-eligible per protocol)
//   Tx[i*n+k] = Poisson(φ_i·x_k)   (Class 2: per response, transient)
//
// Static ClosedShell θ = V0x − E0x + γ with:
//   V0x   = V_local·x − c_xc·groundK(x),  groundK(x)_k = Σ_i φ_i·Tx[i,k]   (= K[φ,φ](x))
//   γ     = Q( J·φ − c_xc·( direct + cross ) )                            (Q on γ ONLY)
//             direct_k = Σ_i x_i·g0[i,k]   (= K[φ,x](φ))
//             cross_k  = Σ_i φ_i·Tx[k,i]   (= K[x,φ](φ) — the TRANSPOSE contraction)
//
// The same Tx feeds BOTH groundK (V0x) and cross (γ) — one build, two contractions.
// This header is the FD ClosedShell slice (Static now; Full = +Ty + the Y block).
// Numerics: the explicit Tx-build + dot-contraction is the SAME math as the
// Exchange operator but a different truncation/accumulation path → A/B-to-thresh,
// not bit-identical. Reference (compute_V0x/compute_gamma) stays the gate.
// =========================================================================

#include "tda.hpp"     // ResponseGroundState, poperatorT, Kernels<> primary template
#include "static.hpp"  // Kernels<Static, ClosedShell> specialization (compute_E0x reuse)
#include "full.hpp"    // Kernels<Full, ClosedShell> (compute_E0x reuse; Y block)

#include <madness/mra/mra.h>
#include <madness/mra/vmra.h>
#include <madness/world/worldprofile.h>  // PROFILE_BLOCK — perf-model meters; no-op unless WORLD_PROFILE_ENABLE

#include <algorithm>
#include <vector>

namespace molresponse_v3 {
namespace exch {

using vecfuncT = std::vector<madness::real_function_3d>;

/// Build several pair-density tensors sharing Poisson waves: for each c in `cs`,
/// T = Poisson(b_i · c_k) (length nb*|c|, row-major [i*|c|+k]). Products for all cs
/// are fused into one apply per row-block (Inc-3a). `tile` (Inc-3b) blocks the
/// φ-row (b) index: tile<=0 or >=nb → ONE block = the single fused wave (peak
/// ~nb·Σ|c| functions); tile=T → ceil(nb/T) blocks, peak ~T·Σ|c|, at the cost of
/// more waves (memory ↔ fences trade). Bit-identical across any tile: truncate()
/// and apply() act per-function, so blocking changes neither the truncation nor
/// the convolution of any individual function, and each result lands at the same
/// [i*|c|+k] slot (doc 33 §2, Inc-3a/3b).
inline std::vector<vecfuncT>
build_pair_tensors(madness::World &world, const poperatorT &coulop,
                   const vecfuncT &b, const std::vector<const vecfuncT *> &cs,
                   double vtol, int tile = 0) {
  const std::size_t nb = b.size();
  std::vector<vecfuncT> out(cs.size());
  for (std::size_t ci = 0; ci < cs.size(); ++ci)
    out[ci] = vecfuncT(nb * cs[ci]->size());
  const std::size_t step =
      (tile > 0 && static_cast<std::size_t>(tile) < nb) ? static_cast<std::size_t>(tile) : nb;
  for (std::size_t i0 = 0; i0 < nb; i0 += step) {
    const std::size_t i1 = std::min(i0 + step, nb);
    vecfuncT prod;
    std::vector<std::size_t> seg_ci, seg_dest;   // per product-run: output idx, dest offset
    for (std::size_t ci = 0; ci < cs.size(); ++ci) {
      const std::size_t nc = cs[ci]->size();
      for (std::size_t i = i0; i < i1; ++i) {
        auto pi = mul_sparse(world, b[i], *cs[ci], vtol);   // b_i · c_k for all k
        prod.insert(prod.end(), pi.begin(), pi.end());
        seg_ci.push_back(ci);
        seg_dest.push_back(i * nc);
      }
    }
    madness::truncate(world, prod, vtol);
    auto T = madness::apply(world, *coulop, prod);   // Poisson over this row-block, ONE wave
    madness::truncate(world, T, vtol);
    std::size_t off = 0;
    for (std::size_t s = 0; s < seg_ci.size(); ++s) {
      const std::size_t nc = cs[seg_ci[s]]->size();
      for (std::size_t k = 0; k < nc; ++k) out[seg_ci[s]][seg_dest[s] + k] = T[off + k];
      off += nc;
    }
  }
  return out;
}

/// Pair-density convolution tensor T[i*nc+k] = Poisson(b_i · c_k), |b|=nb, |c|=nc.
inline vecfuncT
build_pair_tensor(madness::World &world, const poperatorT &coulop,
                  const vecfuncT &b, const vecfuncT &c, double vtol, int tile = 0) {
  return build_pair_tensors(world, coulop, b, {&c}, vtol, tile)[0];
}

/// out_k = Σ_i a_i · T[i*stride + k]   (column contraction; |a|=na, |out|=stride).
inline vecfuncT
contract_col(madness::World &world, const vecfuncT &a, const vecfuncT &T,
             std::size_t stride) {
  const std::size_t na = a.size();
  vecfuncT out(stride);
  for (std::size_t k = 0; k < stride; ++k) {
    vecfuncT col(na);
    for (std::size_t i = 0; i < na; ++i) col[i] = T[i * stride + k];
    out[k] = madness::dot(world, a, col);         // Σ_i a_i · col_i
  }
  return out;
}

/// out_k = Σ_i a_i · T[k*n + i]   (row contraction = transpose; square n×n tensor).
inline vecfuncT
contract_row(madness::World &world, const vecfuncT &a, const vecfuncT &T,
             std::size_t n) {
  vecfuncT out(n);
  for (std::size_t k = 0; k < n; ++k) {
    vecfuncT row(T.begin() + k * n, T.begin() + (k + 1) * n);
    out[k] = madness::dot(world, a, row);
  }
  return out;
}

// ---- Class-1 cache (φ-only) ---------------------------------------------
/// g0[i*n+k] = Poisson(φ_i·φ_k). Built once per protocol (φ fixed); Inc 2 caches
/// it on ResponseGroundState. Standalone callers build it here.
inline vecfuncT
build_g0(madness::World &world, const ResponseGroundState &gs, double vtol,
         int tile = 0) {
  PROFILE_BLOCK(rs_ext_g0_build);   // φ·φ tensor build (once/protocol); doc 33 Inc-3d
  return build_pair_tensor(world, gs.coulop, gs.amo, gs.amo, vtol, tile);
}

// ---- per-response Class-2 context ---------------------------------------
struct ResponseExchangeCtx {
  vecfuncT                  Tx;   // Poisson(φ_i·x_k), n*n
  vecfuncT                  Ty;   // Full only
  madness::real_function_3d J;    // Poisson(ρ)  (Coulomb)
};

/// build_ctx for Static ClosedShell: J + Tx (no Ty).
inline ResponseExchangeCtx
build_ctx_static_cs(madness::World &world, const ResponseGroundState &gs,
                    const ResponseStateX<ClosedShell> &state,
                    const madness::real_function_3d &rho, double vtol,
                    int tile = 0) {
  PROFILE_BLOCK(rs_ext_ctx_build);  // per-iter J + Tx convolutions; doc 33 Inc-3d
  ResponseExchangeCtx ctx;
  ctx.J  = madness::apply(*gs.coulop, rho);
  ctx.Tx = build_pair_tensor(world, gs.coulop, gs.amo, state.x_alpha, vtol, tile);
  return ctx;
}

/// θ = V0x − E0x + γ  for Static ClosedShell, assembled from {g0, ctx}.
/// Mirrors fd_solver's θ EXACTLY: Q applies to the γ block ONLY.
inline ResponseStateX<ClosedShell>
assemble_theta_static_cs(madness::World &world, const ResponseGroundState &gs,
                         const ResponseStateX<ClosedShell> &state,
                         const vecfuncT &g0, const ResponseExchangeCtx &ctx) {
  PROFILE_BLOCK(rs_ext_assemble);   // per-iter contractions + E0x + Q; doc 33 Inc-3d
  const double thr  = madness::FunctionDefaults<3>::get_thresh();
  const double vtol = thr * 0.1;
  const std::size_t n = gs.amo.size();

  // --- V0x (NO Q):  V_local·x − c_xc·groundK ---
  auto theta = mul_sparse(world, gs.V_local_alpha, state.x_alpha, vtol);
  if (gs.c_xc > 0.0) {
    auto groundK = contract_col(world, gs.amo, ctx.Tx, n);   // Σ_i φ_i Tx[i,k]
    madness::gaxpy(world, 1.0, theta, -gs.c_xc, groundK);
  }
  // --- − E0x (NO Q): off-diagonal Fock transform (reuse the kernel) ---
  {
    auto E0x = Kernels<Static, ClosedShell>::compute_E0x(world, gs, state);
    madness::gaxpy(world, 1.0, theta, -1.0, E0x.x_alpha);
  }
  // --- + γ (WITH Q):  Q( J·φ − c_xc·(direct + cross) ) ---
  auto g = mul(world, ctx.J, gs.amo, true);                  // J·φ
  if (gs.c_xc > 0.0) {
    auto direct = contract_col(world, state.x_alpha, g0, n); // Σ_i x_i g0[i,k]
    auto cross  = contract_row(world, gs.amo, ctx.Tx, n);    // Σ_i φ_i Tx[k,i]
    madness::gaxpy(world, 1.0, g, -gs.c_xc, direct);
    madness::gaxpy(world, 1.0, g, -gs.c_xc, cross);
  }
  g = gs.Qa(g);
  madness::truncate(world, g, vtol);
  madness::gaxpy(world, 1.0, theta, 1.0, g);                 // θ = V0x − E0x + γ
  madness::truncate(world, theta, thr);
  return ResponseStateX<ClosedShell>{std::move(theta)};
}

/// build_ctx for Full ClosedShell: J + Tx + Ty (Ty = Poisson(phi_i*y_k)).
/// `rho` = the caller's folded response density (2*sum phi*(x+y)); J = Poisson(rho).
inline ResponseExchangeCtx
build_ctx_full_cs(madness::World &world, const ResponseGroundState &gs,
                  const ResponseStateXY<ClosedShell> &state,
                  const madness::real_function_3d &rho, double vtol,
                  int tile = 0) {
  PROFILE_BLOCK(rs_ext_ctx_build);  // per-iter J + Tx/Ty convolutions; doc 33 Inc-3d
  ResponseExchangeCtx ctx;
  ctx.J  = madness::apply(*gs.coulop, rho);
  // Inc-3a: build Tx and Ty fused (was two build_pair_tensor calls = two waves);
  // Inc-3b: `tile` blocks the φ-row index to bound the n² peak. Bit-identical.
  auto Ts = build_pair_tensors(world, gs.coulop, gs.amo,
                               {&state.x_alpha, &state.y_alpha}, vtol, tile);
  ctx.Tx = std::move(Ts[0]);
  ctx.Ty = std::move(Ts[1]);
  return ctx;
}

/// theta = V0x - E0x + gamma for Full ClosedShell, assembled from {g0, ctx}.
/// Mirrors Kernels<Full,ClosedShell>: Q applies to EACH gamma block (X,Y) ONLY;
/// V0x and E0x are NOT Q-projected. Tensor map (cf. full.hpp compute_gamma):
///   groundK(x) = K[phi,phi](x) = contract_col(phi, Tx, n)   (V0x, X)
///   groundK(y) = K[phi,phi](y) = contract_col(phi, Ty, n)   (V0x, Y)
///   X gamma-exch = K[phi,x](phi)+K[y,phi](phi)
///                = contract_col(x,g0,n) + contract_row(phi,Ty,n)
///   Y gamma-exch = K[phi,y](phi)+K[x,phi](phi)
///                = contract_col(y,g0,n) + contract_row(phi,Tx,n)
/// (MADNESS Exchange K[bra,ket](f)=Σ_i ket_i·P(bra_i·f): {phi,a}→col(a,g0); {a,phi}→row(phi,Ta).)
inline ResponseStateXY<ClosedShell>
assemble_theta_full_cs(madness::World &world, const ResponseGroundState &gs,
                       const ResponseStateXY<ClosedShell> &state,
                       const vecfuncT &g0, const ResponseExchangeCtx &ctx) {
  PROFILE_BLOCK(rs_ext_assemble);   // per-iter contractions + E0x + Q; doc 33 Inc-3d
  const double thr  = madness::FunctionDefaults<3>::get_thresh();
  const double vtol = thr * 0.1;
  const std::size_t n = gs.amo.size();

  // --- V0x (NO Q): V_local*{x,y} - c_xc*groundK{x,y} ---
  auto theta_x = mul_sparse(world, gs.V_local_alpha, state.x_alpha, vtol);
  auto theta_y = mul_sparse(world, gs.V_local_alpha, state.y_alpha, vtol);
  if (gs.c_xc > 0.0) {
    auto gKx = contract_col(world, gs.amo, ctx.Tx, n);   // K[phi,phi](x)
    auto gKy = contract_col(world, gs.amo, ctx.Ty, n);   // K[phi,phi](y)
    madness::gaxpy(world, 1.0, theta_x, -gs.c_xc, gKx);
    madness::gaxpy(world, 1.0, theta_y, -gs.c_xc, gKy);
  }
  // --- - E0x (NO Q): off-diagonal Fock transform of both blocks (reuse kernel) ---
  {
    auto E0 = Kernels<Full, ClosedShell>::compute_E0x(world, gs, state);
    madness::gaxpy(world, 1.0, theta_x, -1.0, E0.x_alpha);
    madness::gaxpy(world, 1.0, theta_y, -1.0, E0.y_alpha);
  }
  // --- + gamma (WITH Q per block): same J*phi Coulomb on X and Y ---
  // MADNESS Exchange (exchangeoperator.cc): K[bra,ket](f) = Σ_i ket_i·P(bra_i·f)
  // — the KET multiplies outside, the BRA forms the density with f. So a reference
  // pair {bra,ket} on apply_to=φ maps to:
  //   {φ,a} → Σ a_i·P(φ_i·φ_k) = contract_col(a, g0);
  //   {a,φ} → Σ φ_i·P(a_i·φ_k) = contract_row(φ, Ta).
  auto gx = mul(world, ctx.J, gs.amo, true);
  auto gy = mul(world, ctx.J, gs.amo, true);
  if (gs.c_xc > 0.0) {
    // γ_X = K[φ,x](φ) + K[y,φ](φ)
    auto x_Kphix = contract_col(world, state.x_alpha, g0,     n);  // {φ,x}: Σ x_i P(φ_i φ_k)
    auto x_Kyphi = contract_row(world, gs.amo,        ctx.Ty, n);  // {y,φ}: Σ φ_i P(y_i φ_k)
    madness::gaxpy(world, 1.0, gx, -gs.c_xc, x_Kphix);
    madness::gaxpy(world, 1.0, gx, -gs.c_xc, x_Kyphi);
    // γ_Y = K[φ,y](φ) + K[x,φ](φ)
    auto y_Kphiy = contract_col(world, state.y_alpha, g0,     n);  // {φ,y}: Σ y_i P(φ_i φ_k)
    auto y_Kxphi = contract_row(world, gs.amo,        ctx.Tx, n);  // {x,φ}: Σ φ_i P(x_i φ_k)
    madness::gaxpy(world, 1.0, gy, -gs.c_xc, y_Kphiy);
    madness::gaxpy(world, 1.0, gy, -gs.c_xc, y_Kxphi);
  }
  gx = gs.Qa(gx);
  gy = gs.Qa(gy);
  madness::truncate(world, gx, vtol);
  madness::truncate(world, gy, vtol);
  madness::gaxpy(world, 1.0, theta_x, 1.0, gx);  // theta_x = V0x - E0x + gamma_x
  madness::gaxpy(world, 1.0, theta_y, 1.0, gy);  // theta_y = V0y - E0y + gamma_y
  madness::truncate(world, theta_x, thr);
  madness::truncate(world, theta_y, thr);
  return ResponseStateXY<ClosedShell>{std::move(theta_x), std::move(theta_y)};
}

// One-call gate-1 entry points, overloaded on the response type: build g0 + ctx
// (Inc 1: g0 built per call, not yet cached -- Inc 2 caches it on the ground
// state) then assemble theta. fd_solver's --fd-tensor branch calls this; only
// ClosedShell Static/Full are defined, so the caller guards non-ClosedShell
// instantiations with an `if constexpr` on the State type.
inline ResponseStateX<ClosedShell>
assemble_theta_tensor(madness::World &world, const ResponseGroundState &gs,
                      const ResponseStateX<ClosedShell> &state,
                      const madness::real_function_3d &rho, int tile = 0) {
  const double vtol = madness::FunctionDefaults<3>::get_thresh() * 0.1;
  // g0 is φ-only → use the per-protocol cache (Inc 2); build per-call only if cold
  // (standalone callers, e.g. test_exchange_ctx, leave gs.g0_alpha empty).
  vecfuncT g0_local;
  if (gs.g0_alpha.empty()) g0_local = build_g0(world, gs, vtol, tile);
  const vecfuncT &g0 = gs.g0_alpha.empty() ? g0_local : gs.g0_alpha;
  auto ctx = build_ctx_static_cs(world, gs, state, rho, vtol, tile);
  return assemble_theta_static_cs(world, gs, state, g0, ctx);
}
inline ResponseStateXY<ClosedShell>
assemble_theta_tensor(madness::World &world, const ResponseGroundState &gs,
                      const ResponseStateXY<ClosedShell> &state,
                      const madness::real_function_3d &rho, int tile = 0) {
  const double vtol = madness::FunctionDefaults<3>::get_thresh() * 0.1;
  // g0 is φ-only → use the per-protocol cache (Inc 2); build per-call only if cold.
  vecfuncT g0_local;
  if (gs.g0_alpha.empty()) g0_local = build_g0(world, gs, vtol, tile);
  const vecfuncT &g0 = gs.g0_alpha.empty() ? g0_local : gs.g0_alpha;
  auto ctx = build_ctx_full_cs(world, gs, state, rho, vtol, tile);
  return assemble_theta_full_cs(world, gs, state, g0, ctx);
}

} // namespace exch
} // namespace molresponse_v3

#endif // MOLRESPONSE_V3_KERNELS_EXCHANGE_CTX_HPP
