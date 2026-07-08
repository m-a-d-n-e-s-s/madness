#ifndef MOLRESPONSE_V3_ESSOLVERGUESS_HPP
#define MOLRESPONSE_V3_ESSOLVERGUESS_HPP

#include "GroundState.hpp"
#include "ResponseFunctions.hpp"
#include "ResponseKernel.hpp"          // for orthonormalize_bundle, ResponseType
#include "kernels/common_ops.hpp"      // apply_exchange (K for the virtual-AO Fock)

#include <madness/chem/molecularbasis.h>     // AtomicBasisSet (VirtualAO guess)
#include <madness/chem/molecular_functors.h> // madchem::AtomicBasisFunctor
#include <madness/mra/mra.h>
#include <madness/mra/operator.h>

#include <algorithm>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

namespace molresponse_v3 {

using namespace madness;

// -------------------------------------------------------------------------
// Simple Gaussian functor used to build a localized envelope around atoms.
// Independent copy of the same form used in molresponse_legacy/TDDFT.h:66
// and molresponse/ExcitedResponse.cpp:122 — kept local so this header is
// self-contained.
// -------------------------------------------------------------------------
template <std::size_t NDIM>
class ESGaussianGuess : public FunctionFunctorInterface<double, NDIM> {
  typedef Vector<double, NDIM> coordT;

public:
  ESGaussianGuess(const coordT &origin, double exponent,
                  std::vector<int> ijk = std::vector<int>(NDIM))
      : origin(origin), exponent(exponent), ijk(std::move(ijk)) {}

  coordT origin;
  double exponent;
  std::vector<int> ijk;

  double operator()(const coordT &xyz) const {
    double r2 = 0.0, prefac = 1.0;
    for (std::size_t i = 0; i < NDIM; i++) {
      const double d = xyz[i] - origin[i];
      r2 += d * d;
      prefac *= std::pow(xyz[i], ijk[i]);
    }
    return prefac * std::exp(-exponent * r2);
  }
};

/// Build a sum-of-Gaussians envelope localized at every atom of the
/// molecule. Width controlled by `exponent` (smaller = more diffuse).
inline real_function_3d build_atom_envelope(World &world,
                                            const Molecule &molecule,
                                            double exponent = 0.01) {
  real_function_3d envelope = real_factory_3d(world);
  for (const auto &atom : molecule.get_atoms()) {
    Vector<double, 3> origin{atom.x, atom.y, atom.z};
    real_function_3d g = real_factory_3d(world).functor(real_functor_3d(
        new ESGaussianGuess<3>(origin, exponent, std::vector<int>{0, 0, 0})));
    envelope = envelope + g;
  }
  return envelope;
}

/// In-place: add random noise (per leaf) to each function in the vector.
/// Mirrors molresponse_legacy/TDDFT.cc:3720-3750 (`add_randomness`).
inline void add_random_noise(World &world, vector_real_function_3d &v,
                             double magnitude = 1.0e3) {
  auto noise = [magnitude](const Key<3> & /*key*/, Tensor<double> &x) {
    Tensor<double> y(x.size());
    y.fillrandom();
    y.scale(magnitude);
    x = x + y;
  };
  for (auto &f : v)
    f.unaryop(noise);
}

/// Initial guess for **TDA closed-shell (RHF)** excited states.
///
/// Algorithm (mirrors `molresponse_legacy/TDDFT.cc::create_random_guess` /
/// `molresponse/ExcitedResponse::make_random_trial`):
///   1. Create N×n_occ zero functions.
///   2. Add per-leaf random noise.
///   3. Multiply by sum-of-Gaussians envelope localized at each atom.
///   4. Q-project to remove occupied-space components.
///   5. Orthonormalize across the N roots.
inline std::vector<RealResponseState>
make_initial_guess_tda_rhf(World &world, GroundState &gs, long num_roots,
                           double envelope_exponent = 0.01,
                           double random_magnitude = 1.0e3) {

  MADNESS_CHECK(gs.is_spin_restricted());
  const long n_alpha = gs.num_alpha();

  // Atom-localized Gaussian envelope (one shared envelope across roots —
  // randomness comes from add_random_noise).
  auto envelope = build_atom_envelope(world, gs.molecule(), envelope_exponent);

  std::vector<RealResponseState> X(num_roots);
  for (long s = 0; s < num_roots; s++) {
    X[s] = RealResponseState::allocate(world, n_alpha, /*n_beta=*/0,
                                       /*include_y=*/false);
    // Step 1+2: noise into each x_alpha[i]
    add_random_noise(world, X[s].x_alpha, random_magnitude);
    // Step 3: localize at atoms
    for (auto &f : X[s].x_alpha)
      f = envelope * f;
  }

  // Step 4: Q-project each state's x_alpha against ground orbitals
  const auto &Q = gs.Q_alpha();
  for (auto &state : X) {
    state.x_alpha = Q(state.x_alpha);
  }

  // Step 5: Gram-Schmidt across roots (TDA metric is identity)
  orthonormalize_bundle(world, ResponseType::TDA, X);

  return X;
}

/// Initial guess for **TDA open-shell (UHF)** excited states.
///
/// Same algorithm as the RHF version but generates BOTH α and β
/// response channels per state — n_alpha α-functions and n_beta
/// β-functions each. Each spin channel gets independent random noise,
/// is multiplied by the atom envelope, Q-projected against the
/// corresponding ground orbitals, and finally `orthonormalize_bundle`
/// orthogonalizes across roots using the combined (α + β) inner
/// product.
inline std::vector<RealResponseState>
make_initial_guess_tda_uhf(World &world, GroundState &gs, long num_roots,
                           double envelope_exponent = 0.01,
                           double random_magnitude = 1.0e3) {

  MADNESS_CHECK(!gs.is_spin_restricted());
  const long n_alpha = gs.num_alpha();
  const long n_beta  = gs.num_beta();

  auto envelope = build_atom_envelope(world, gs.molecule(), envelope_exponent);

  std::vector<RealResponseState> X(num_roots);
  for (long s = 0; s < num_roots; s++) {
    X[s] = RealResponseState::allocate(world, n_alpha, n_beta,
                                       /*include_y=*/false);
    add_random_noise(world, X[s].x_alpha, random_magnitude);
    add_random_noise(world, X[s].x_beta,  random_magnitude);
    for (auto &f : X[s].x_alpha) f = envelope * f;
    for (auto &f : X[s].x_beta)  f = envelope * f;
  }

  // Q-project against ground orbitals — per spin.
  const auto &Qa = gs.Q_alpha();
  const auto &Qb = gs.Q_beta();
  for (auto &state : X) {
    state.x_alpha = Qa(state.x_alpha);
    state.x_beta  = Qb(state.x_beta);
  }

  // Gram-Schmidt across roots — orthonormalize_bundle uses both α and β
  // inner-product contributions, matching the OpenShell-TDA subspace
  // metric in ESSolver<TDA, OpenShell>.
  orthonormalize_bundle(world, ResponseType::TDA, X);

  return X;
}

// Cartesian moment functor: r → x^i · y^j · z^k.
// Direct port of the local copy in molresponse_legacy/TDDFT.cc (the
// comment there says SCF.cc's version "wasn't linking right" so it was
// duplicated). Kept self-contained inside this header.
class BS_MomentFunctor : public FunctionFunctorInterface<double, 3> {
private:
  const int i, j, k;

public:
  BS_MomentFunctor(int i, int j, int k) : i(i), j(j), k(k) {}
  explicit BS_MomentFunctor(const std::vector<int> &x)
      : i(x[0]), j(x[1]), k(x[2]) {}
  double operator()(const Vector<double, 3> &r) const override {
    double xi = 1.0, yj = 1.0, zk = 1.0;
    for (int p = 0; p < i; ++p) xi *= r[0];
    for (int p = 0; p < j; ++p) yj *= r[1];
    for (int p = 0; p < k; ++p) zk *= r[2];
    return xi * yj * zk;
  }
};

inline double kronecker(std::size_t l, std::size_t n) {
  return (l == n) ? 1.0 : 0.0;
}

/// Solid spherical harmonics up to angular momentum L = `n`, keyed by {l, m}.
/// Direct port of `TDDFT::solid_harmonics` (legacy TDDFT.cc:413-478).
/// `inline` so it can live in a header without ODR violations.
inline std::map<std::vector<int>, real_function_3d>
solid_harmonics(World &world, int n) {
  std::map<std::vector<int>, real_function_3d> result;

  // Basic x, y, z, constant, and zero.
  real_function_3d x = real_factory_3d(world).functor(
      real_functor_3d(new BS_MomentFunctor(std::vector<int>{1, 0, 0})));
  real_function_3d y = real_factory_3d(world).functor(
      real_functor_3d(new BS_MomentFunctor(std::vector<int>{0, 1, 0})));
  real_function_3d z = real_factory_3d(world).functor(
      real_functor_3d(new BS_MomentFunctor(std::vector<int>{0, 0, 1})));
  real_function_3d c = real_factory_3d(world).functor(
      real_functor_3d(new BS_MomentFunctor(std::vector<int>{0, 0, 0})));
  real_function_3d zero = real_factory_3d(world);

  // Seed (assume n >= 1).
  result[std::vector<int>{0, 0}] = copy(c);
  result[std::vector<int>{0, -1}] = zero;
  result[std::vector<int>{0, 1}] = zero;
  result[std::vector<int>{-1, 0}] = zero;

  // Generate the solid harmonics recursively (verbatim from legacy).
  for (int l = 0; l < n; l++) {
    result[std::vector<int>{l + 1, l + 1}] =
        sqrt(pow(2, kronecker(l, 0) * (2 * l) / (2 * l + 1))) *
        (x * result[std::vector<int>{l, l}] -
         (1 - kronecker(l, 0) * y * result[std::vector<int>{l, -l}]));
    result[std::vector<int>{l + 1, -l - 1}] =
        sqrt(pow(2, kronecker(l, 0) * (2 * l) / (2 * l + 1))) *
        (y * result[std::vector<int>{l, l}] +
         (1 - kronecker(l, 0) * x * result[std::vector<int>{l, -l}]));

    // Recursion needs out-of-range neighbors as zero placeholders.
    result[std::vector<int>{l + 1, l + 2}] = zero;
    result[std::vector<int>{l + 1, -l - 2}] = zero;

    for (int m = -l; m < l + 1; m++) {
      result[std::vector<int>{l + 1, m}] =
          1.0 / std::sqrt((l + m + 1) * (l - m + 1)) *
          ((2 * l + 1) * z * result[std::vector<int>{l, m}] -
           sqrt((l + m) * (l - m)) * (x * x + y * y + z * z) *
               result[std::vector<int>{l - 1, m}]);
    }
  }

  // Drop zeros and the constant {0,0}.
  for (auto it = result.begin(); it != result.end();) {
    if (it->second.norm2() == 0)
      it = result.erase(it);
    else
      ++it;
  }
  result.erase(std::vector<int>{0, 0});

  return result;
}

/// Initial guess for **TDA closed-shell (RHF)** excited states using
/// solid-harmonic-times-orbital trial functions.
///
/// Mirrors `TDDFT::create_trial_functions` (legacy TDDFT.cc:535-595):
///   1. Build solid harmonics up to L = max(2, ceil(sqrt(num_roots/n_occ) - 1))
///      so that solids.size() * n_occ >= num_roots.
///   2. For each (orbital, harmonic) pair, build a *single-excitation* trial
///      state: only one orbital position is non-zero and equals
///      `harmonic * orbital[mirror_index]`. This yields a richer subspace
///      than putting the same harmonic on every orbital.
///   3. Q-project each state to remove ground-state components.
///   4. Truncate.
///   5. Gram-Schmidt across the bundle.
///
/// Returns `num_roots` states (or fewer if not enough harmonic-orbital
/// pairs are available, which only happens for very large num_roots).
inline std::vector<RealResponseState>
create_solid_harmonics_guess(World &world, GroundState &gs, long num_roots) {
  MADNESS_CHECK(gs.is_spin_restricted());
  const long n_occ = gs.num_alpha();
  MADNESS_CHECK(n_occ >= 1);

  // Same heuristic as legacy create_trial_functions:
  // ensure solids.size() * n_occ > num_roots.
  const int max_l = static_cast<int>(std::max(
      2.0,
      std::ceil(std::sqrt(static_cast<double>(num_roots) /
                          static_cast<double>(n_occ)) -
                1.0)));
  auto solids = solid_harmonics(world, max_l);
  if (world.rank() == 0) {
    print("ES guess: created ", solids.size(),
          " solid harmonics (L<=", max_l, ") for ", num_roots,
          " trial states across ", n_occ, " occupied orbitals");
  }

  const auto &orbitals = gs.orbitals_alpha();
  std::vector<RealResponseState> X;
  X.reserve(num_roots);

  long count = 0;
  for (long i = 0; i < n_occ && count < num_roots; ++i) {
    for (const auto &kv : solids) {
      if (count >= num_roots) break;
      const real_function_3d &harmonic = kv.second;
      // Single-excitation trial: one occupied position non-zero.
      RealResponseState state = RealResponseState::allocate(
          world, n_occ, /*n_beta=*/0, /*include_y=*/false);
      const long pos = count % n_occ;
      const long orb_idx = n_occ - (count % n_occ) - 1;
      state.x_alpha[pos] = harmonic * orbitals[orb_idx];
      X.push_back(std::move(state));
      ++count;
    }
  }

  // Q-project each state.
  const auto &Q = gs.Q_alpha();
  for (auto &state : X) {
    state.x_alpha = Q(state.x_alpha);
  }

  // Truncate.
  const double thresh = FunctionDefaults<3>::get_thresh();
  for (auto &state : X) {
    truncate(world, state.x_alpha, thresh);
  }

  // Gram-Schmidt across roots (TDA metric is identity).
  orthonormalize_bundle(world, ResponseType::TDA, X);

  return X;
}

/// Initial guess for **TDA open-shell (UHF)** excited states using
/// solid-harmonic-times-orbital trial functions.
///
/// Mirrors `create_solid_harmonics_guess` (RHF) but emits both α and β
/// single-excitation trial states. Each candidate excites a single
/// orbital in a single spin channel (the other channel remains all zero
/// for that candidate). Candidates are interleaved by spin so a
/// downselect that picks the lowest-N by ω sees both spins represented
/// — important when the SOMO transition is the lowest-ω one.
///
/// Total candidate count is `#harmonics × (n_alpha + n_beta)` ≥
/// `num_roots`; the L heuristic uses `(n_alpha + n_beta)` in place of
/// `n_occ`.
inline std::vector<RealResponseState>
create_solid_harmonics_guess_uhf(World &world, GroundState &gs, long num_roots) {
  MADNESS_CHECK(!gs.is_spin_restricted());
  const long n_a = gs.num_alpha();
  const long n_b = gs.num_beta();
  MADNESS_CHECK(n_a + n_b >= 1);

  const long n_pool = n_a + n_b;
  const int max_l = static_cast<int>(std::max(
      2.0,
      std::ceil(std::sqrt(static_cast<double>(num_roots) /
                          static_cast<double>(n_pool)) -
                1.0)));
  auto solids = solid_harmonics(world, max_l);
  if (world.rank() == 0) {
    print("ES guess (UHF): created ", solids.size(),
          " solid harmonics (L<=", max_l, ") for ", num_roots,
          " trial states across n_alpha=", n_a, " n_beta=", n_b);
  }

  // Flat (spin, harmonic) list, alternating spin per harmonic so the
  // sequence interleaves α-candidates with β-candidates.
  std::vector<std::pair<int, real_function_3d>> spin_harm;
  spin_harm.reserve(solids.size() * 2);
  for (const auto &kv : solids) {
    if (n_a > 0) spin_harm.emplace_back(0, kv.second);
    if (n_b > 0) spin_harm.emplace_back(1, kv.second);
  }

  const auto &amo = gs.orbitals_alpha();
  const auto &bmo = gs.orbitals_beta();
  const long n_occ_max = std::max(n_a, n_b);

  std::vector<RealResponseState> X;
  X.reserve(num_roots);

  long count = 0;
  for (long orb_iter = 0; orb_iter < n_occ_max && count < num_roots;
       ++orb_iter) {
    for (const auto &sh : spin_harm) {
      if (count >= num_roots) break;
      const int spin = sh.first;
      const real_function_3d &harmonic = sh.second;
      const long n_this = (spin == 0) ? n_a : n_b;
      if (orb_iter >= n_this) continue;

      RealResponseState state = RealResponseState::allocate(
          world, n_a, n_b, /*include_y=*/false);
      const long pos     = orb_iter;
      const long orb_idx = n_this - 1 - orb_iter;
      if (spin == 0) {
        state.x_alpha[pos] = harmonic * amo[orb_idx];
      } else {
        state.x_beta[pos]  = harmonic * bmo[orb_idx];
      }
      X.push_back(std::move(state));
      ++count;
    }
  }

  // Q-project per spin
  const auto &Qa = gs.Q_alpha();
  const auto &Qb = gs.Q_beta();
  for (auto &state : X) {
    state.x_alpha = Qa(state.x_alpha);
    state.x_beta  = Qb(state.x_beta);
  }

  // Truncate
  const double thresh = FunctionDefaults<3>::get_thresh();
  for (auto &state : X) {
    truncate(world, state.x_alpha, thresh);
    truncate(world, state.x_beta,  thresh);
  }

  // Gram-Schmidt across roots — orthonormalize_bundle uses combined α+β
  // metric, matching ESSolver<TDA, OpenShell>.
  orthonormalize_bundle(world, ResponseType::TDA, X);

  return X;
}

/// Initial guess for **TDA closed-shell (RHF)** excited states from VIRTUAL
/// ORBITALS built off a Gaussian AO basis — the "NWChem guess" (doc 17, Path A;
/// port of `molresponse/ExcitedResponse::create_virtual_ao_guess`, energy-ordered).
///
/// Why: Random/SolidHarmonics are "cold" — they don't span the low-energy
/// excitation space in energy order, so the iterative solver misses/scrambles
/// roots. This builds approximate virtual orbitals and seeds the lowest
/// (ε_a − ε_i) single excitations (= NWChem's CIS-diagonal Davidson guess).
///
/// Algorithm:
///   1. Project a Gaussian AO basis (`basis_name`) onto the MRA grid.
///   2. Q-project to the virtual space.
///   3. SVD the overlap; drop near-linear-dependent directions (σ < 0.05).
///   4. Build the ground Fock `F = T + V_local − c_xc·K` (same pieces as the
///      kernels' compute_V0x, plus kinetic T).
///   5. Diagonalize ⟨φ|F|φ⟩ → approximate virtuals {φ_a, ε_a}; re-Q-project.
///   6. Keep virtuals with ε_a > ε_HOMO.
///   7. Emit single excitations i→a, ordered by (ε_a − ε_i), lowest `num_roots`:
///      trial = φ_a in occupied slot i (the v3 x_alpha[n_occ] shape).
///   8. Q-project + Gram-Schmidt across roots.
/// May return fewer than `num_roots` if the basis is too small (caller tops up).
inline std::vector<RealResponseState>
create_virtual_ao_guess(World &world, GroundState &gs, long num_roots,
                        const std::string &basis_name = "aug-cc-pvdz") {
  MADNESS_CHECK(gs.is_spin_restricted());
  const long   n_occ  = gs.num_alpha();
  MADNESS_CHECK(n_occ >= 1);
  const auto  &phi0   = gs.orbitals_alpha();
  const auto  &eps0   = gs.energies_alpha();
  const auto  &Q      = gs.Q_alpha();
  const double c_xc   = gs.hf_exchange_coefficient();
  const double lo     = gs.params().lo();
  const double thresh = FunctionDefaults<3>::get_thresh();

  // 1. Project the AO basis to MRA functions.
  AtomicBasisSet aobasis(basis_name);
  const int nbf = aobasis.nbf(gs.molecule());
  vector_real_function_3d ao(nbf);
  for (int i = 0; i < nbf; ++i)
    ao[i] = real_factory_3d(world).functor(real_functor_3d(
        new madchem::AtomicBasisFunctor(
            aobasis.get_atomic_basis_function(gs.molecule(), i))));
  world.gop.fence();

  // 2. Q-project into the virtual space.
  ao = Q(ao);

  // 3. SVD the overlap; keep directions with σ >= 0.05 (drop Q-induced lindeps).
  Tensor<double> S = matrix_inner(world, ao, ao), U, sigma, VT;
  svd(S, U, sigma, VT);
  long keep = 0;
  while (keep < sigma.size() && sigma(keep) >= 0.05) ++keep;
  vector_real_function_3d xao = transform(world, ao, U);
  if (keep < static_cast<long>(xao.size()))
    xao.erase(xao.begin() + keep, xao.end());
  if (xao.empty()) return {};   // degenerate basis → caller falls back

  // 4. Ground Fock F·x = T·x + V_local·x − c_xc·K[φ0,φ0](x)  (= compute_V0x + T).
  auto F = [&](const vector_real_function_3d &x) {
    real_derivative_3d Dx(world, 0), Dy(world, 1), Dz(world, 2);
    auto Tx = (apply(world, Dx, apply(world, Dx, x)) +
               apply(world, Dy, apply(world, Dy, x)) +
               apply(world, Dz, apply(world, Dz, x))) * (-0.5);
    auto Vx = mul_sparse(world, gs.V_local(), x, thresh * 0.1);
    if (c_xc > 0.0) {
      auto Kx = common_ops::apply_exchange(world, phi0, phi0, x, lo);
      gaxpy(world, 1.0, Vx, -c_xc, Kx);
    }
    return Tx + Vx;
  };

  // 5. Diagonalize ⟨xao|F|xao⟩ → virtual orbitals + energies; re-Q-project.
  Tensor<double> FF = matrix_inner(world, xao, F(xao)), C, e_a;
  syev(FF, C, e_a);
  vector_real_function_3d phi_a = Q(transform(world, xao, C));

  // 6. Keep virtuals above the HOMO; 7. order single excitations by ε_a − ε_i.
  const double e_homo = eps0(n_occ - 1);
  struct Exc { long i; long a; double de; };
  std::vector<Exc> exc;
  for (long a = 0; a < e_a.size(); ++a) {
    if (e_a(a) <= e_homo) continue;
    for (long i = 0; i < n_occ; ++i)
      exc.push_back({i, a, e_a(a) - eps0(i)});
  }
  std::sort(exc.begin(), exc.end(),
            [](const Exc &p, const Exc &q) { return p.de < q.de; });
  const long n_take = std::min<long>(num_roots, static_cast<long>(exc.size()));

  std::vector<RealResponseState> X;
  X.reserve(n_take);
  for (long t = 0; t < n_take; ++t) {
    RealResponseState s =
        RealResponseState::allocate(world, n_occ, /*n_beta=*/0, /*include_y=*/false);
    s.x_alpha[exc[t].i] = copy(phi_a[exc[t].a]);
    X.push_back(std::move(s));
  }

  // 8. Q-project + Gram-Schmidt across roots.
  for (auto &s : X) s.x_alpha = Q(s.x_alpha);
  orthonormalize_bundle(world, ResponseType::TDA, X);

  if (world.rank() == 0)
    print("ES guess (VirtualAO): basis=", basis_name, " nbf=", nbf,
          " virtuals>HOMO=", (long)exc.size() / std::max<long>(1, n_occ),
          " trials=", (long)X.size(), " of", num_roots);
  return X;
}

/// Initial-guess generation mode for ESSolver<TDA, *>. Selected via the
/// test/binary's `--guess=` CLI knob.
///
///   Random         — per-leaf random noise inside an atom-centered
///                    Gaussian envelope, Q-projected and GS-orthonormalized.
///                    No symmetry hint — relies entirely on the warmup
///                    iterations (or KAIN in the main solve) to find the
///                    low-ω subspace.
///
///   SolidHarmonics — (default) single-excitation trials built from
///                    solid harmonics × occupied orbitals. L is auto-
///                    sized so #harmonics × n_occ ≥ num_roots. For L=1
///                    this is just {x, y, z}·φ_occ — pure dipole, which
///                    is essentially the answer for dipole-allowed
///                    transitions.
///
///   VirtualAO      — virtual orbitals built off a Gaussian AO basis,
///                    Fock-diagonalized, seeded as the lowest (ε_a−ε_i) single
///                    excitations (NWChem CIS-diagonal guess; closed-shell only).
///                    See create_virtual_ao_guess (doc 17, Path A).
enum class ESGuessMode { Random, SolidHarmonics, VirtualAO };

inline const char *to_string(ESGuessMode m) {
  switch (m) {
    case ESGuessMode::Random:         return "random";
    case ESGuessMode::SolidHarmonics: return "solid_harmonics";
    case ESGuessMode::VirtualAO:      return "virtual_ao";
  }
  return "unknown";
}

inline ESGuessMode parse_es_guess_mode(const std::string &s) {
  if (s == "random")                            return ESGuessMode::Random;
  if (s == "solid" || s == "solid_harmonics")   return ESGuessMode::SolidHarmonics;
  if (s == "virtual" || s == "virtual_ao")      return ESGuessMode::VirtualAO;
  throw std::runtime_error(
      "parse_es_guess_mode: unknown guess mode '" + s +
      "' (expected: random | solid_harmonics | virtual_ao)");
}

} // namespace molresponse_v3

#endif // MOLRESPONSE_V3_ESSOLVERGUESS_HPP
