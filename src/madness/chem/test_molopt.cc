// Unit tests for MolOpt::make_projector and MolOpt::optimize / optimize_app,
// covering single atoms, linear molecules (aligned along Cartesian axes and
// tilted), non-linear molecules, and non-finite gradient rejection under both
// strict IEEE-754 and -ffast-math builds.

#include <madness/chem/molopt.h>
#include <madness/misc/misc.h>
#include <madness/world/MADworld.h>

#include <array>
#include <cmath>
#include <cstdint>
#include <limits>
#include <vector>

using namespace madness;

namespace {

/// Builds a proper 3D rotation matrix (Rodrigues' formula) rotating around a
/// non-Cartesian axis so all 3 lab coordinates mix strongly.
Tensor<double> make_test_rotation_matrix() {
  // Axis proportional to (1, 2, 3), rotation angle 1.1 rad
  const double inv_n = 1.0 / std::sqrt(14.0);
  const double ux = 1.0 * inv_n, uy = 2.0 * inv_n, uz = 3.0 * inv_n;
  const double theta = 1.1;
  const double c = std::cos(theta);
  const double s = std::sin(theta);
  const double omc = 1.0 - c;

  Tensor<double> Q(3L, 3L);
  Q(0, 0) = c + ux * ux * omc;
  Q(0, 1) = ux * uy * omc - uz * s;
  Q(0, 2) = ux * uz * omc + uy * s;
  Q(1, 0) = uy * ux * omc + uz * s;
  Q(1, 1) = c + uy * uy * omc;
  Q(1, 2) = uy * uz * omc - ux * s;
  Q(2, 0) = uz * ux * omc - uy * s;
  Q(2, 1) = uz * uy * omc + ux * s;
  Q(2, 2) = c + uz * uz * omc;
  return Q;
}

/// Applies a rigid 3D rotation Q and translation shift to all atoms.
Molecule transform_molecule(const Molecule& mol, const Tensor<double>& Q,
                            const std::array<double, 3>& shift, double eprec) {
  std::vector<Atom> atoms;
  atoms.reserve(mol.natom());
  for (std::size_t i = 0; i < mol.natom(); ++i) {
    const Atom& a = mol.get_atom(i);
    const double rx =
        Q(0, 0) * a.x + Q(0, 1) * a.y + Q(0, 2) * a.z + shift[0];
    const double ry =
        Q(1, 0) * a.x + Q(1, 1) * a.y + Q(1, 2) * a.z + shift[1];
    const double rz =
        Q(2, 0) * a.x + Q(2, 1) * a.y + Q(2, 2) * a.z + shift[2];
    atoms.emplace_back(rx, ry, rz, a.q, a.atomic_number);
  }
  return {atoms, eprec};
}

/// Checks that P = make_projector(mol) is finite, symmetric, idempotent
/// (P^2 == P), has rank == expected_rank, and annihilates rigid-body
/// translations and rotations while preserving internal stretches.
bool check_projector(const std::string& label, const Molecule& mol,
                     int expected_rank) {
  MolOpt opt(/*maxiter=*/5, 0.1, 1e-6, 1e-5, 1e-5, 1e-7, 1e-7, /*print_level=*/0);
  const Tensor<double> P = opt.make_projector(mol);
  const int n = 3 * static_cast<int>(mol.natom());

  bool ok = (P.ndim() == 2 && P.dim(0) == n && P.dim(1) == n) &&
            is_finite(P);
  if (!ok) {
    print("FAILED:", label, "projector dimensions or finiteness");
    return false;
  }

  // Symmetry: P == P^T
  const double sym_err = (P - transpose(P)).absmax();
  // Idempotency: P * P == P
  const double idem_err = (inner(P, P) - P).absmax();
  // Rank = trace(P) for a symmetric idempotent matrix
  double tr = 0.0;
  for (int i = 0; i < n; ++i) tr += P(i, i);
  const double rank_err = std::abs(tr - static_cast<double>(expected_rank));

  // Rigid translations must be annihilated
  const long natom = static_cast<long>(mol.natom());
  double trans_err = 0.0;
  for (long axis = 0; axis < 3; ++axis) {
    Tensor<double> t(n);
    for (long a = 0; a < natom; ++a) t(3 * a + axis) = 1.0;
    trans_err = std::max(trans_err, inner(P, t).absmax());
  }

  // Infinitesimal rotations about the origin must be annihilated
  double rot_err = 0.0;
  const Tensor<double> coords = mol.get_all_coords();
  for (int raxis = 0; raxis < 3; ++raxis) {
    Tensor<double> r(n);
    for (long a = 0; a < natom; ++a) {
      const double x = coords(a, 0), y = coords(a, 1), z = coords(a, 2);
      if (raxis == 0) {
        r(3 * a + 1) = z;
        r(3 * a + 2) = -y;
      } else if (raxis == 1) {
        r(3 * a + 0) = z;
        r(3 * a + 2) = -x;
      } else {
        r(3 * a + 0) = -y;
        r(3 * a + 1) = x;
      }
    }
    rot_err = std::max(rot_err, inner(P, r).absmax());
  }

  ok = sym_err < 1e-13 && idem_err < 1e-13 && rank_err < 1e-12 &&
       trans_err < 1e-12 && rot_err < 1e-12;
  print(label, "trace=", tr, "sym_err=", sym_err, "idem_err=", idem_err,
        "trans_err=", trans_err, "rot_err=", rot_err, ok ? "ok" : "FAILED");
  return ok;
}

/// Harmonic diatomic potential with a deliberate translational net-force bias
/// in the raw gradient, exercising MolOpt's gradient projection on a linear
/// molecule aligned along a Cartesian axis.
struct HarmonicDiatomicTarget {
  double r0 = 3.0;
  double k = 0.5;
  double spurious_net_force_z = 2.0e-4;

  static double bond_length(const Tensor<double>& x) {
    const double dx = x[0] - x[3];
    const double dy = x[1] - x[4];
    const double dz = x[2] - x[5];
    return std::sqrt(dx * dx + dy * dy + dz * dz);
  }

  [[nodiscard]] double value(const Tensor<double>& x) const {
    const double dr = bond_length(x) - r0;
    return 0.5 * k * dr * dr;
  }

  void energy_and_gradient(Molecule& mol, double& energy,
                           Tensor<double>& gradient) const {
    const Tensor<double> x = mol.get_all_coords().flat();
    const double r = bond_length(x);
    const double dr = r - r0;
    energy = 0.5 * k * dr * dr;
    gradient = Tensor<double>(6L);
    for (int c = 0; c < 3; ++c) {
      const double d = (x[c] - x[3 + c]) / r;
      gradient[c] = k * dr * d;
      gradient[3 + c] = -k * dr * d;
    }
    // Add an identical spurious force on both atoms (net translation)
    gradient[2] += 0.5 * spurious_net_force_z;
    gradient[5] += 0.5 * spurious_net_force_z;
  }
};

/// Rotationally and translationally invariant harmonic potential for a nearly
/// linear bent triatomic (atoms 0, 1, 2 with atom 1 near the midpoint of 0-2),
/// injecting both a spurious net translation and a spurious rigid-body force
/// along the near-linear rotation mode around the 0-2 axis.
struct NearLinearBentTriatomicTarget {
  double L0 = 4.0;
  double s0 = 0.0;
  double b0 = 1.0e-7;
  double k_L = 0.5;
  double k_s = 0.4;
  double k_b = 0.6;
  double spurious_net_force = 2.0e-4;
  double spurious_near_linear_rot_force = 2.0e-4;

  static void decompose(const Tensor<double>& x, std::array<double, 3>& uhat,
                        std::array<double, 3>& wperp,
                        std::array<double, 3>& nhat, double& L, double& s,
                        double& b) {
    const std::array<double, 3> u = {x[6] - x[0], x[7] - x[1], x[8] - x[2]};
    L = std::sqrt(u[0] * u[0] + u[1] * u[1] + u[2] * u[2]);
    uhat = {u[0] / L, u[1] / L, u[2] / L};
    const std::array<double, 3> w = {x[3] - 0.5 * (x[0] + x[6]),
                                     x[4] - 0.5 * (x[1] + x[7]),
                                     x[5] - 0.5 * (x[2] + x[8])};
    s = w[0] * uhat[0] + w[1] * uhat[1] + w[2] * uhat[2];
    wperp = {w[0] - s * uhat[0], w[1] - s * uhat[1], w[2] - s * uhat[2]};
    b = std::sqrt(wperp[0] * wperp[0] + wperp[1] * wperp[1] +
                  wperp[2] * wperp[2]);
    if (b > 0.0) {
      nhat = {wperp[0] / b, wperp[1] / b, wperp[2] / b};
    } else {
      // Any unit vector orthogonal to uhat suffices when b == 0.
      const int min_axis =
          (std::abs(uhat[0]) <= std::abs(uhat[1]) &&
           std::abs(uhat[0]) <= std::abs(uhat[2]))
              ? 0
              : (std::abs(uhat[1]) <= std::abs(uhat[2]) ? 1 : 2);
      std::array<double, 3> e = {0.0, 0.0, 0.0};
      e[min_axis] = 1.0;
      const double dot = e[0] * uhat[0] + e[1] * uhat[1] + e[2] * uhat[2];
      const std::array<double, 3> p = {
          e[0] - dot * uhat[0], e[1] - dot * uhat[1], e[2] - dot * uhat[2]};
      const double pn = std::sqrt(p[0] * p[0] + p[1] * p[1] + p[2] * p[2]);
      nhat = {p[0] / pn, p[1] / pn, p[2] / pn};
    }
  }

  [[nodiscard]] double value(const Tensor<double>& x) const {
    std::array<double, 3> uhat{}, wperp{}, nhat{};
    double L = 0.0, s = 0.0, b = 0.0;
    decompose(x, uhat, wperp, nhat, L, s, b);
    const double dL = L - L0;
    const double ds = s - s0;
    const double db = b - b0;
    return 0.5 * k_L * dL * dL + 0.5 * k_s * ds * ds + 0.5 * k_b * db * db;
  }

  void energy_and_gradient(Molecule& mol, double& energy,
                           Tensor<double>& gradient) const {
    const Tensor<double> x = mol.get_all_coords().flat();
    std::array<double, 3> uhat{}, wperp{}, nhat{};
    double L = 0.0, s = 0.0, b = 0.0;
    decompose(x, uhat, wperp, nhat, L, s, b);

    const double fL = k_L * (L - L0);
    const double fs = k_s * (s - s0);
    const double fb = k_b * (b - b0);
    energy = 0.5 * k_L * (L - L0) * (L - L0) +
             0.5 * k_s * (s - s0) * (s - s0) +
             0.5 * k_b * (b - b0) * (b - b0);

    // Out-of-plane unit vector t = uhat x nhat; (-1/2 t, +t, -1/2 t) is an
    // exact infinitesimal rigid rotation about uhat through the centroid.
    const std::array<double, 3> that = {
        uhat[1] * nhat[2] - uhat[2] * nhat[1],
        uhat[2] * nhat[0] - uhat[0] * nhat[2],
        uhat[0] * nhat[1] - uhat[1] * nhat[0]};

    gradient = Tensor<double>(9L);
    for (int c = 0; c < 3; ++c) {
      const double gw = fs * uhat[c] + fb * nhat[c];
      const double gu =
          fL * uhat[c] + (fs / L) * wperp[c] - (fb * s / L) * nhat[c];
      gradient[0 + c] = -0.5 * gw - gu + spurious_net_force / 3.0 -
                        0.5 * spurious_near_linear_rot_force * that[c];
      gradient[3 + c] = gw + spurious_net_force / 3.0 +
                        spurious_near_linear_rot_force * that[c];
      gradient[6 + c] = -0.5 * gw + gu + spurious_net_force / 3.0 -
                        0.5 * spurious_near_linear_rot_force * that[c];
    }
  }
};

struct NonFiniteGradientTarget {
  long nan_index = 2;
  long dim = 6;
  double base_grad = 0.0;
  double e0 = -1.0;
  [[nodiscard]] double value(const Tensor<double>&) const { return e0; }
  void energy_and_gradient(Molecule&, double& energy,
                           Tensor<double>& gradient) const {
    energy = e0;
    gradient = Tensor<double>(dim);
    for (long i = 0; i < dim; ++i) gradient[i] = base_grad * ((i % 2 == 0) ? 1.0 : -1.0);
    gradient[nan_index] = runtime_quiet_nan();
  }
};

struct NonFiniteEnergyTarget {
  [[nodiscard]] double value(const Tensor<double>&) const { return -1.0; }
  void energy_and_gradient(Molecule&, double& energy,
                           Tensor<double>& gradient) const {
    energy = runtime_quiet_nan();
    gradient = Tensor<double>(6L);
    gradient[2] = 0.05;
    gradient[5] = -0.05;
  }
};

struct NonFiniteLineSearchEnergyTarget {
  [[nodiscard]] double value(const Tensor<double>&) const {
    return runtime_quiet_nan();
  }
  void energy_and_gradient(Molecule&, double& energy,
                           Tensor<double>& gradient) const {
    energy = -1.0;
    gradient = Tensor<double>(6L);
    gradient[2] = 0.05;
    gradient[5] = -0.05;
  }
};

}  // namespace

int main(int argc, char** argv) {
  initialize(argc, argv);
  const double eprec = 1e-6;
  bool success = true;

  // 1. Single atom (0 internal degrees of freedom; all 3 rotation generators vanish)
  {
    const Molecule he({Atom(0.0, 0.0, 0.0, 2.0, 2)}, eprec);
    success = check_projector("single atom at origin", he, 0) && success;
    const Molecule he_off({Atom(1.2, -0.7, 2.5, 2.0, 2)}, eprec);
    success = check_projector("single atom off origin", he_off, 0) && success;
  }

  // 2. Linear diatomic aligned along z (the Molecule::orient() geometry for LiH),
  //    along x, along y, and tilted off-axis (3*2 - 5 = 1 internal degree of freedom)
  {
    const Molecule lih_z(
        {Atom(0.0, 0.0, 0.7625, 3.0, 3), Atom(0.0, 0.0, -2.2875, 1.0, 1)},
        eprec);
    success = check_projector("diatomic along z", lih_z, 1) && success;

    const Molecule lih_x(
        {Atom(0.7625, 0.0, 0.0, 3.0, 3), Atom(-2.2875, 0.0, 0.0, 1.0, 1)},
        eprec);
    success = check_projector("diatomic along x", lih_x, 1) && success;

    const Molecule lih_y(
        {Atom(0.0, 0.7625, 0.0, 3.0, 3), Atom(0.0, -2.2875, 0.0, 1.0, 1)},
        eprec);
    success = check_projector("diatomic along y", lih_y, 1) && success;

    const Molecule lih_tilted(
        {Atom(1.0, 2.0, 3.0, 3.0, 3), Atom(-0.5, -1.0, -1.5, 1.0, 1)}, eprec);
    success = check_projector("diatomic tilted", lih_tilted, 1) && success;
  }

  // 3. Linear triatomic CO2 along z (3*3 - 5 = 4 internal degrees of freedom),
  //    including ~1e-15 bohr floating-point roundoff on a single atom (must NOT
  //    amplify roundoff into a frozen internal bending vibration).
  {
    const Molecule co2({Atom(0.0, 0.0, 0.0, 6.0, 6),
                        Atom(0.0, 0.0, 2.19, 8.0, 8),
                        Atom(0.0, 0.0, -2.19, 8.0, 8)},
                       eprec);
    success = check_projector("linear triatomic along z", co2, 4) && success;

    const Molecule co2_roundoff({Atom(1.0e-15, -0.5e-15, 0.0, 6.0, 6),
                                 Atom(0.0, 0.0, 2.19, 8.0, 8),
                                 Atom(0.0, 0.0, -2.19, 8.0, 8)},
                                eprec);
    success = check_projector("linear triatomic with 1e-15 roundoff",
                              co2_roundoff, 4) &&
              success;
  }

  // 4. Non-linear polyatomic H2O (3*3 - 6 = 3 internal degrees of freedom),
  //    plus a nearly linear bent triatomic (valid rotation of norm ~8.16e-8
  //    bohr must still be projected out -> rank 3) both axis-aligned AND
  //    rotated off-axis in 3D (verifying rotational covariance of P), vs. a
  //    strictly collinear slightly tilted triatomic (dependent rotation -> rank 4).
  {
    const Molecule h2o({Atom(0.0, 0.0, 0.213, 8.0, 8),
                        Atom(0.0, 1.427, -0.852, 1.0, 1),
                        Atom(0.0, -1.427, -0.852, 1.0, 1)},
                       eprec);
    success = check_projector("non-linear H2O", h2o, 3) && success;

    const Molecule nearly_linear_bent({Atom(0.0, 0.0, -2.0, 8.0, 8),
                                       Atom(1.0e-7, 0.0, 0.0, 6.0, 6),
                                       Atom(0.0, 0.0, 2.0, 8.0, 8)},
                                      eprec);
    success = check_projector("nearly linear bent triatomic (axis-aligned)",
                              nearly_linear_bent, 3) &&
              success;

    const Tensor<double> Q = make_test_rotation_matrix();
    const Molecule nearly_linear_bent_rot =
        transform_molecule(nearly_linear_bent, Q, {1.3, -2.4, 0.9}, eprec);
    success = check_projector("nearly linear bent triatomic (3D rotated)",
                              nearly_linear_bent_rot, 3) &&
              success;

    // Verify rotational covariance: P(Q*X + shift) == diag(Q,Q,Q) P(X) diag(Q^T,Q^T,Q^T)
    MolOpt opt(/*maxiter=*/5, 0.1, 1e-6, 1e-5, 1e-5, 1e-7, 1e-7,
               /*print_level=*/0);
    const Tensor<double> P0 = opt.make_projector(nearly_linear_bent);
    const Tensor<double> Prot = opt.make_projector(nearly_linear_bent_rot);
    Tensor<double> Qblk(9L, 9L);
    for (int a = 0; a < 3; ++a) {
      for (int r = 0; r < 3; ++r) {
        for (int c = 0; c < 3; ++c) {
          Qblk(3 * a + r, 3 * a + c) = Q(r, c);
        }
      }
    }
    const Tensor<double> P_cov = inner(Qblk, inner(P0, transpose(Qblk)));
    const double cov_err = (Prot - P_cov).absmax();
    const bool cov_ok = cov_err < 1e-6;
    print("nearly linear bent projector 3D covariance err=", cov_err,
          cov_ok ? "ok" : "FAILED");
    success = cov_ok && success;

    const Molecule collinear_slight_tilt({Atom(0.0, 0.0, 0.0, 6.0, 6),
                                          Atom(1.0e-8, 0.0, 2.19, 8.0, 8),
                                          Atom(-1.0e-8, 0.0, -2.19, 8.0, 8)},
                                         eprec);
    success = check_projector("collinear slightly tilted triatomic",
                              collinear_slight_tilt, 4) &&
              success;
  }

  // 5. End-to-end optimization of a z-aligned diatomic with spurious net force
  //    and of a 3D-rotated nearly linear bent triatomic with spurious net force
  //    AND spurious torque along the near-linear rotation mode (gtol = 1e-6,
  //    well below the 2e-4 spurious components removed by the projector)
  try {
    Molecule lih(
        {Atom(0.0, 0.0, 0.7625, 3.0, 3), Atom(0.0, 0.0, -2.2875, 1.0, 1)},
        eprec);
    HarmonicDiatomicTarget target;
    MolOpt opt(/*maxiter=*/15, 0.2, 1e-8, 1e-6, 1e-6, 1e-10, 1e-8,
               /*print_level=*/0, "bfgs");
    const OptimizationResults res = opt.optimize_app(lih, target);
    const double r_final = res.final_geometry.inter_atomic_distance(0, 1);
    const bool ok = res.nsteps > 0 && is_finite(res.max_gradient) &&
                    res.max_gradient < 1e-6 && std::abs(r_final - 3.0) < 1e-5;
    print("diatomic optimize_app: nsteps=", res.nsteps,
          "max_g=", res.max_gradient, "r_final=", r_final,
          ok ? "ok" : "FAILED");
    success = ok && success;

    const Molecule bent_init_axis({Atom(0.0, 0.0, -2.05, 8.0, 8),
                                   Atom(1.2e-7, 0.0, 0.04, 6.0, 6),
                                   Atom(0.0, 0.0, 2.05, 8.0, 8)},
                                  eprec);
    const Tensor<double> Q = make_test_rotation_matrix();
    const Molecule bent_init =
        transform_molecule(bent_init_axis, Q, {-0.8, 1.5, 2.2}, eprec);
    NearLinearBentTriatomicTarget bent_target;
    MolOpt bent_opt(/*maxiter=*/25, 0.2, 1e-12, 1e-8, 1e-8, 1e-14, 1e-10,
                    /*print_level=*/0, "bfgs");
    const OptimizationResults bent_res =
        bent_opt.optimize_app(bent_init, bent_target);
    std::array<double, 3> uhat{}, wperp{}, nhat{};
    double L_final = 0.0, s_final = 0.0, b_final = 0.0;
    NearLinearBentTriatomicTarget::decompose(
        bent_res.final_geometry.get_all_coords().flat(), uhat, wperp, nhat,
        L_final, s_final, b_final);
    const bool bent_ok =
        bent_res.nsteps > 0 && is_finite(bent_res.max_gradient) &&
        bent_res.max_gradient < 1e-8 && std::abs(L_final - 4.0) < 1e-7 &&
        std::abs(s_final - 0.0) < 1e-7 && std::abs(b_final - 1.0e-7) < 1e-9;
    print("3D near-linear bent optimize_app: nsteps=", bent_res.nsteps,
          "max_g=", bent_res.max_gradient, "L=", L_final, "s=", s_final,
          "b=", b_final, bent_ok ? "ok" : "FAILED");
    success = bent_ok && success;

    // Also verify decompose handles exact collinearity (b == 0) finitely
    Tensor<double> collinear_x(9L);
    collinear_x[0] = 0.0; collinear_x[1] = 0.0; collinear_x[2] = -2.0;
    collinear_x[3] = 0.0; collinear_x[4] = 0.0; collinear_x[5] = 0.0;
    collinear_x[6] = 0.0; collinear_x[7] = 0.0; collinear_x[8] = 2.0;
    NearLinearBentTriatomicTarget::decompose(collinear_x, uhat, wperp, nhat,
                                             L_final, s_final, b_final);
    const bool b0_ok = (b_final == 0.0) && is_finite(nhat[0]) &&
                       is_finite(nhat[1]) && is_finite(nhat[2]);
    print("NearLinearBentTriatomicTarget::decompose at b=0:",
          b0_ok ? "ok" : "FAILED");
    success = b0_ok && success;
  } catch (const MadnessException& ex) {
    print("FAILED: unexpected MadnessException in end-to-end optimize_app:",
          ex.msg);
    success = false;
  }

  // 6. Non-finite gradient and energy detection:
  //    - Directly test madness::is_finite(const Tensor<T>&) on vectors where
  //      Tensor::absmax() returns a finite value despite a non-leading NaN
  //      (and on leading-NaN, complex, empty, and non-contiguous slice views).
  //    - Verify MolOpt::optimize and MolOpt::optimize_app throw on non-finite
  //      gradients and non-finite initial or line-search energies rather than
  //      silently forcing convergence. (Note: dense inner(P, graw) multiplies
  //      every P(i, j) * graw[j] without skipping zeros, so 0 * NaN = NaN
  //      propagates into all entries of g; MolOpt checks both graw and g.)
  //    - Verify QuasiNewton and SteepestDescent (which do not project by P)
  //      reject a [0.1, 0.1, NaN] gradient where g.absmax() is finite.
  {
    // Direct test of the Tensor::absmax() NaN blind spot on [0.1, NaN, 0.2]:
    // std::max(result, std::abs(x)) keeps result when x is NaN, so absmax()
    // returns finite 0.2, whereas madness::is_finite(tail_nan) returns false.
    Tensor<double> tail_nan(3L);
    tail_nan[0] = 0.1;
    tail_nan[1] = runtime_quiet_nan();
    tail_nan[2] = 0.2;
    Tensor<double> lead_nan(3L);
    lead_nan[0] = runtime_quiet_nan();
    lead_nan[1] = 0.1;
    lead_nan[2] = 0.2;
    Tensor<double> all_ok(3L);
    all_ok[0] = 0.1;
    all_ok[1] = -0.3;
    all_ok[2] = 0.2;

    // Non-contiguous slice view (column 1 of a 3x2 matrix has stride 2)
    Tensor<double> mat(3L, 2L);
    mat = 1.0;
    Tensor<double> col0 = mat(_, 0);
    Tensor<double> col1 = mat(_, 1);
    const bool slice_clean_ok = !col1.iscontiguous() && is_finite(col1);
    mat(1, 1) = runtime_quiet_nan();
    const bool slice_nan_ok = is_finite(col0) && !is_finite(col1);

    Tensor<double_complex> cten(2L);
    cten[0] = double_complex(1.0, 2.0);
    cten[1] = double_complex(1.0, runtime_quiet_nan());

    const bool tensor_finite_ok =
        is_finite(Tensor<double>()) && is_finite(all_ok) &&
        is_finite(tail_nan.absmax()) && !is_finite(tail_nan) &&
        !is_finite(lead_nan) && slice_clean_ok && slice_nan_ok &&
        !is_finite(cten);
    print("Tensor is_finite vs absmax on [0.1, NaN, 0.2] and slices:",
          tensor_finite_ok ? "ok" : "FAILED");
    success = tensor_finite_ok && success;

    const Molecule lih(
        {Atom(0.0, 0.0, 0.7625, 3.0, 3), Atom(0.0, 0.0, -2.2875, 1.0, 1)},
        eprec);
    const Molecule h2o({Atom(0.0, 0.0, 0.213, 8.0, 8),
                        Atom(0.0, 1.427, -0.852, 1.0, 1),
                        Atom(0.0, -1.427, -0.852, 1.0, 1)},
                       eprec);

    auto expect_throws = [&](const char* label, const Molecule& mol,
                             auto& bad_target) {
      bool threw_app = false;
      bool threw_opt = false;
      {
        MolOpt opt(/*maxiter=*/5, 0.1, 1e-5, 1e-4, 1e-4, 1e-6, 1e-5,
                   /*print_level=*/0, "bfgs");
        try {
          (void)opt.optimize_app(mol, bad_target);
        } catch (const MadnessException&) {
          threw_app = true;
        }
      }
      {
        MolOpt opt(/*maxiter=*/5, 0.1, 1e-5, 1e-4, 1e-4, 1e-6, 1e-5,
                   /*print_level=*/0, "bfgs");
        try {
          (void)opt.optimize(mol, bad_target);
        } catch (const MadnessException&) {
          threw_opt = true;
        }
      }
      const bool ok = threw_app && threw_opt;
      print(label, ok ? "ok" : "FAILED");
      success = ok && success;
    };

    // NaN in internal z-component of z-aligned diatomic
    NonFiniteGradientTarget bad_grad_z{/*nan_index=*/2, /*dim=*/6, 0.0};
    expect_throws("non-finite gradient (z component) throws:", lih, bad_grad_z);

    // NaN in x-component of z-aligned diatomic (rigid-body direction)
    NonFiniteGradientTarget bad_grad_x{/*nan_index=*/0, /*dim=*/6, 0.05};
    expect_throws("non-finite gradient (x component) throws:", lih, bad_grad_x);

    // NaN in trailing element of 3-atom H2O raw gradient
    NonFiniteGradientTarget bad_grad_tail{/*nan_index=*/8, /*dim=*/9, 0.05};
    expect_throws("non-finite gradient (trailing element) throws:", h2o,
                  bad_grad_tail);

    NonFiniteEnergyTarget bad_energy;
    expect_throws("non-finite energy throws:", lih, bad_energy);

    NonFiniteLineSearchEnergyTarget bad_ls_energy;
    expect_throws("non-finite line_search energy throws:", lih, bad_ls_energy);

    // Verify is_nan / is_finite / nearlyEqual on real/complex NaN and +/-Inf
    const double qnan_d = runtime_quiet_nan<double>();
    const float qnan_f = runtime_quiet_nan<float>();
    const std::complex<double> qnan_c(1.0, qnan_d);
    volatile double pos_inf_v = std::numeric_limits<double>::infinity();
    const double pos_inf = pos_inf_v;
    const double neg_inf = -pos_inf;
    const bool nan_helpers_ok =
        is_nan(qnan_d) && !is_finite(qnan_d) && is_nan(qnan_f) &&
        !is_finite(qnan_f) && is_nan(qnan_c) && !is_finite(qnan_c) &&
        !is_nan(pos_inf) && !is_finite(pos_inf) &&
        !is_nan(std::complex<double>(1.0, -2.0)) &&
        is_finite(std::complex<double>(1.0, -2.0)) &&
        !nearlyEqual(qnan_d, qnan_d) && !nearlyEqual(qnan_d, 1.0) &&
        nearlyEqual(pos_inf, pos_inf) && nearlyEqual(neg_inf, neg_inf) &&
        !nearlyEqual(pos_inf, neg_inf) && !nearlyEqual(pos_inf, 1.0) &&
        nearlyEqual(1.0, 1.0 + 1e-12);
    print("is_nan / is_finite / runtime_quiet_nan / nearlyEqual helpers:",
          nan_helpers_ok ? "ok" : "FAILED");
    success = nan_helpers_ok && success;

    // Verify QuasiNewton and SteepestDescent also reject non-finite values/gradients
    struct BadOptTarget : public OptimizationTargetInterface {
      bool nan_grad = false;
      bool provides_gradient() const override { return true; }
      double value(const Tensor<double>&) override {
        return nan_grad ? 1.0 : runtime_quiet_nan();
      }
      Tensor<double> gradient(const Tensor<double>& x) override {
        Tensor<double> g(x.dim(0));
        g = 0.1;
        if (nan_grad) g[x.dim(0) - 1] = runtime_quiet_nan();
        return g;
      }
    };
    bool qn_sd_ok = true;
    for (bool ng : {false, true}) {
      auto t = std::make_shared<BadOptTarget>();
      t->nan_grad = ng;
      bool qn_threw = false, sd_threw = false;
      try {
        Tensor<double> x(3L);
        QuasiNewton qn(t);
        (void)qn.optimize(x);
      } catch (const MadnessException&) {
        qn_threw = true;
      }
      try {
        Tensor<double> x(3L);
        SteepestDescent sd(t);
        (void)sd.optimize(x);
      } catch (const MadnessException&) {
        sd_threw = true;
      }
      qn_sd_ok = qn_sd_ok && qn_threw && sd_threw;
    }
    print("QuasiNewton / SteepestDescent non-finite rejection:",
          qn_sd_ok ? "ok" : "FAILED");
    success = qn_sd_ok && success;
  }

  // 7. Molecule::nuclear_attraction_potential_derivative at r == 0 is finite
  {
    const Molecule lih(
        {Atom(0.0, 0.0, 0.7625, 3.0, 3), Atom(0.0, 0.0, -2.2875, 1.0, 1)},
        eprec);
    bool ok = true;
    for (int axis = 0; axis < 3; ++axis) {
      const double dv =
          lih.nuclear_attraction_potential_derivative(0, axis, 0.0, 0.0, 0.7625);
      ok = ok && is_finite(dv) && dv == 0.0;
    }
    print("nuclear_attraction_potential_derivative at r=0:",
          ok ? "ok" : "FAILED");
    success = ok && success;
  }

  print("test_molopt", success ? "passed" : "FAILED");
  finalize();
  return success ? 0 : 1;
}
