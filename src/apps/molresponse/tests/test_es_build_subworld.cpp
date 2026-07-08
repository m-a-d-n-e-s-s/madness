// test_es_build_subworld — S1 of the persistent-subworld ES design (docs/23).
//
// Proves the one piece the φ-count (Inc 2, test_node_phi) and A/S-scalar
// (keystone, test_subspace_allreduce) proofs did NOT cover: the GS-DEPENDENT
// build kernels reproduce bit-for-bit when run in a NODE-ALIGNED subworld.
// Specifically the per-root build  Λ = T0x + V0x − E0x_full + γ  where
//   V0x = V_local·x − c_xc·K[φ,φ](x)   (exchange),
//   γ   = Q( J[ρ]·φ − c_xc·K[φ,x](φ) ) (exchange),
//   E0x_full = focka-transform of x.
// Synthetic gaussians stand in for the GS + roots: the kernels are deterministic
// in their inputs, so world-independence holds regardless of physical meaning.
//
// Comparison (keystone style — scalars, functions never subtracted across worlds):
//   A_ij = ⟨X_i | Λ_j⟩  (M×M), built in the universe vs in the node-subworld.
//   PRIMARY  : max|A_sub − A_uni|  with DIRECT exchange both sides   (world-independence)
//   DIAGNOSTIC: max|A_cached − A_direct| in the universe             (the cached-K0 single-
//               World reference vs the direct-K0 subworld build → tells us whether S2's
//               A/B will be bit-identical or only truncation-level).
//
//   mpirun -np <N> ./test_es_build_subworld [n_occ] [M_roots]
// PASS iff the PRIMARY max diff < 1e-9.

#include "../solvers/node_subworlds.hpp"
#include "../kernels/tda.hpp"               // Kernels<TDA,ClosedShell>, ResponseGroundState
#include "../kernels/assembly.hpp"          // assemble_lambda
#include "../kernels/response_space_ops.hpp"// rs::inner, response_space

#include <madness/mra/mra.h>
#include <madness/mra/operator.h>           // CoulombOperatorPtr
#include <madness/mra/vmra.h>               // truncate
#include <madness/world/MADworld.h>
#include <madness/world/ranks_and_hosts.h>

#include <cmath>
#include <cstdlib>
#include <vector>

using namespace madness;
using namespace molresponse_v3;
using K      = Kernels<TDA, ClosedShell>;
using StateX = ResponseStateX<ClosedShell>;

namespace {
class GaussFunctor : public FunctionFunctorInterface<double, 3> {
  coord_3d c_;
  double   a_;
public:
  GaussFunctor(const coord_3d &c, double a) : c_(c), a_(a) {}
  double operator()(const coord_3d &r) const override {
    const double x = r[0] - c_[0], y = r[1] - c_[1], z = r[2] - c_[2];
    return std::exp(-a_ * (x * x + y * y + z * z));
  }
};

// GS orbitals (n_occ gaussians) built in `world`.
std::vector<real_function_3d> make_orbitals(World &world, int n_occ) {
  std::vector<real_function_3d> amo(n_occ);
  for (int i = 0; i < n_occ; ++i) {
    coord_3d c{{-3.0 + 1.3 * i, 0.4 * ((i % 3) - 1), -0.2 * (i % 2)}};
    amo[i] = real_factory_3d(world).functor(
        real_functor_3d(new GaussFunctor(c, 0.7 + 0.05 * i)));
  }
  truncate(world, amo);
  return amo;
}

// Smooth local potential (a gaussian well) built in `world`.
real_function_3d make_vloc(World &world) {
  real_function_3d v = real_factory_3d(world).functor(
      real_functor_3d(new GaussFunctor(coord_3d{{0.0, 0.0, 0.0}}, 0.3)));
  v.scale(-1.0);
  v.truncate();
  return v;
}

// Assemble a closed-shell GS in `world` around amo/vloc that are ALREADY
// world-local (built in `world`, or shipped in via copy(world, f) so their pmap
// matches `world`). NOTE: building GS functions fresh with real_factory_3d in a
// SUBWORLD gives them the GLOBAL (universe) pmap, so operator Isends target ranks
// outside the subworld -> MPI_ERR_RANK. The fix (and what the real design does
// via Cloud) is to SHIP the GS in. `cached_k0`: true -> cached Exchange operator
// (single-World reference path); false -> K0=null -> DIRECT apply_exchange (the
// subworld path of the persistent design, doc 22 §8).
ResponseGroundState assemble_gs(World &world,
                                std::vector<real_function_3d> amo,
                                real_function_3d vloc, bool cached_k0) {
  ResponseGroundState g0;
  const int n_occ = static_cast<int>(amo.size());
  g0.amo = std::move(amo);
  g0.V_local_alpha = std::move(vloc);
  // fixed symmetric Fock (replicated tensor -> identical across worlds)
  Tensor<double> F(n_occ, n_occ);
  for (int i = 0; i < n_occ; ++i)
    for (int j = 0; j < n_occ; ++j)
      F(i, j) = (i == j) ? -(0.6 + 0.1 * i) : -0.05 / (1.0 + std::abs(i - j));
  g0.focka = F;
  Tensor<double> Fnd = copy(F);
  for (int i = 0; i < n_occ; ++i) Fnd(i, i) = 0.0;
  g0.focka_no_diag = Fnd;
  g0.aeps = Tensor<double>(n_occ);
  for (int i = 0; i < n_occ; ++i) g0.aeps(i) = F(i, i);
  g0.c_xc = 1.0;
  g0.lo   = 1.0e-10;
  g0.coulop = poperatorT(CoulombOperatorPtr(
      world, g0.lo, FunctionDefaults<3>::get_thresh() * 1.0e-3));
  g0.Qa = QProjector<double, 3>(g0.amo);
  g0.K0_alpha = cached_k0 ? common_ops::make_ground_exchange(world, g0.amo, g0.lo)
                          : nullptr;   // null -> DIRECT exchange
  return g0;
}

std::vector<StateX> build_lambda(World &world, const ResponseGroundState &g0,
                                 const std::vector<StateX> &X) {
  std::vector<StateX> lambda(X.size());
  for (std::size_t s = 0; s < X.size(); ++s) {
    auto rho   = K::compute_density(world, g0, X[s]);
    auto gamma = K::compute_gamma(world, g0, X[s], rho);
    auto V0x   = K::compute_V0x(world, g0, X[s]);
    auto T0x   = K::compute_T0x(world, g0, X[s]);
    auto E0f   = K::compute_E0x_full(world, g0, X[s]);
    lambda[s]  = assemble_lambda(world, T0x, V0x, E0f, gamma);
  }
  return lambda;
}

// A(i,j) = sum_p <X_i[p] | Lambda_j[p]>  (M x M), collective in `world`.
Tensor<double> build_A(const std::vector<StateX> &X,
                       const std::vector<StateX> &L) {
  response_space Xrs, Lrs;
  for (auto &r : X) Xrs.push_back(r.x_alpha);
  for (auto &l : L) Lrs.push_back(l.x_alpha);
  return rs::inner(Xrs, Lrs);
}

double maxabs_diff(const Tensor<double> &a, const Tensor<double> &b) {
  double m = 0.0;
  for (long i = 0; i < a.dim(0); ++i)
    for (long j = 0; j < a.dim(1); ++j)
      m = std::max(m, std::abs(a(i, j) - b(i, j)));
  return m;
}
} // namespace

int main(int argc, char **argv) {
  World &universe = initialize(argc, argv);
  startup(universe, argc, argv, true);

  FunctionDefaults<3>::set_k(8);
  FunctionDefaults<3>::set_thresh(1e-6);
  FunctionDefaults<3>::set_cubic_cell(-20.0, 20.0);

  const int n_occ = (argc > 1) ? std::atoi(argv[1]) : 5;
  const int M     = (argc > 2) ? std::atoi(argv[2]) : 3;

  // --- response roots X in the universe (M roots, each n_occ functions) ---
  std::vector<StateX> X(M);
  for (int s = 0; s < M; ++s) {
    X[s].x_alpha.resize(n_occ);
    for (int p = 0; p < n_occ; ++p) {
      coord_3d c{{-2.5 + 1.1 * p + 0.3 * s, 0.2 * ((p + s) % 3 - 1), 0.1 * s}};
      X[s].x_alpha[p] = real_factory_3d(universe).functor(
          real_functor_3d(new GaussFunctor(c, 0.8 + 0.04 * p + 0.03 * s)));
    }
    truncate(universe, X[s].x_alpha);
  }
  universe.gop.fence();

  // --- GS orbitals + local potential, built ONCE in the universe ---
  auto amo_u  = make_orbitals(universe, n_occ);
  auto vloc_u = make_vloc(universe);
  universe.gop.fence();

  // --- UNIVERSE reference builds (universe pmap), BEFORE touching the subworld ---
  Tensor<double> A_uni_direct, A_uni_cached;
  {
    ResponseGroundState g0d =
        assemble_gs(universe, copy(universe, amo_u), copy(universe, vloc_u), false);
    auto Ld = build_lambda(universe, g0d, X);
    A_uni_direct = build_A(X, Ld);
    universe.gop.fence();

    ResponseGroundState g0c =
        assemble_gs(universe, copy(universe, amo_u), copy(universe, vloc_u), true);
    auto Lc = build_lambda(universe, g0c, X);
    A_uni_cached = build_A(X, Lc);
    universe.gop.fence();
  }
  universe.gop.fence();

  // --- SUBWORLD build: ship GS+roots in, build Λ locally, form A_sub ---
  // Canonical subworld discipline (what MacroTaskQ does — macrotaskq.h:710/853):
  // point the GLOBAL FunctionDefaults pmap at the subworld so EVERY build
  // intermediate is subworld-local, then restore the universe pmap before
  // sub.reset(). Without this, op outputs inherit the 16-rank universe pmap and
  // Isend to ranks outside the 8-rank subworld -> MPI_ERR_RANK.
  NodeSubworldInfo info;
  auto sub = make_node_aligned_subworld(universe, &info);
  Tensor<double> A_sub;
  {
    FunctionDefaults<3>::set_default_pmap(*sub);
    std::vector<StateX> Xs(M);
    for (int s = 0; s < M; ++s) {
      Xs[s].x_alpha.resize(n_occ);
      for (int p = 0; p < n_occ; ++p)
        Xs[s].x_alpha[p] = copy(*sub, X[s].x_alpha[p]);
    }
    std::vector<real_function_3d> amo_s(n_occ);
    for (int i = 0; i < n_occ; ++i) amo_s[i] = copy(*sub, amo_u[i]);
    real_function_3d vloc_s = copy(*sub, vloc_u);
    sub->gop.fence();

    ResponseGroundState g0s =
        assemble_gs(*sub, std::move(amo_s), std::move(vloc_s), false);
    auto Ls = build_lambda(*sub, g0s, Xs);
    A_sub = build_A(Xs, Ls);
    sub->gop.fence();
  }  // Xs / g0s / Ls (subworld Functions) destruct here, sub still alive
  sub->gop.fence();
  FunctionDefaults<3>::set_default_pmap(universe);   // restore BEFORE sub.reset()
  sub.reset();
  universe.gop.fence();

  const double diff_world = maxabs_diff(A_sub, A_uni_direct);       // PRIMARY
  const double diff_xmode = maxabs_diff(A_uni_cached, A_uni_direct);// diagnostic

  int ok = (diff_world < 1e-9) ? 1 : 0;
  universe.gop.fence();
  universe.gop.sum(&ok, 1);
  const bool all_ok = (ok == universe.size());

  if (universe.rank() == 0) {
    print("\n=== ES build-in-subworld bit-identity (n_occ=", n_occ, " M=", M,
          " nodes=", info.n_nodes, ") ===");
    print("  built Λ = T0x + V0x − E0x_full + γ  per root; compared A_ij = <X_i|Λ_j>");
    print("  PRIMARY  max|A_sub − A_uni| (direct exch both) =", diff_world,
          " (tol 1e-9)");
    print("  diag     max|A_cached − A_direct| (universe)   =", diff_xmode,
          " (cached-K0 ref vs direct-K0 subworld path)");
    print("\nES_BUILD_SUBWORLD_TEST:", all_ok ? "PASS" : "FAIL");
  }

  // Teardown: universe Functions before finalize (subworld already torn down above).
  for (auto &s : X) s.x_alpha.clear();
  amo_u.clear();
  vloc_u = real_function_3d();
  universe.gop.fence();

  finalize();
  return all_ok ? 0 : 1;
}
