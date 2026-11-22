// Unit tests for MolOpt::make_projector and MolOpt::optimize / optimize_app,
// covering single atoms, linear molecules (aligned along Cartesian axes and
// tilted), non-linear molecules, and non-finite gradient rejection under both
// strict IEEE-754 and -ffast-math builds.

#include <madness/chem/molopt.h>
#include <madness/misc/misc.h>
#include <madness/world/MADworld.h>

#include <cmath>
#include <limits>
#include <vector>

using namespace madness;

namespace {

bool tensor_all_finite(const Tensor<double>& t) {
  for (long i = 0; i < t.size(); ++i) {
    if (!is_finite(t.ptr()[i])) return false;
  }
  return true;
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
            tensor_all_finite(P);
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

struct NonFiniteGradientTarget {
  double e0 = -1.0;
  [[nodiscard]] double value(const Tensor<double>&) const { return e0; }
  void energy_and_gradient(Molecule&, double& energy,
                           Tensor<double>& gradient) const {
    energy = e0;
    gradient = Tensor<double>(6L);
    volatile std::uint64_t nan_bits = 0x7ff8000000000000ULL;
    const std::uint64_t bits = nan_bits;
    gradient[2] = std::bit_cast<double>(bits);
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

  // 3. Linear triatomic CO2 along z (3*3 - 5 = 4 internal degrees of freedom)
  {
    const Molecule co2({Atom(0.0, 0.0, 0.0, 6.0, 6),
                        Atom(0.0, 0.0, 2.19, 8.0, 8),
                        Atom(0.0, 0.0, -2.19, 8.0, 8)},
                       eprec);
    success = check_projector("linear triatomic along z", co2, 4) && success;
  }

  // 4. Non-linear polyatomic H2O (3*3 - 6 = 3 internal degrees of freedom),
  //    plus a nearly linear bent triatomic (valid z-rotation of norm ~8.16e-8
  //    bohr must still be projected out -> rank 3) vs. a strictly collinear
  //    slightly tilted triatomic (dependent rotation -> rank 4).
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
    success =
        check_projector("nearly linear bent triatomic", nearly_linear_bent, 3) &&
        success;

    const Molecule collinear_slight_tilt({Atom(0.0, 0.0, 0.0, 6.0, 6),
                                          Atom(1.0e-8, 0.0, 2.19, 8.0, 8),
                                          Atom(-1.0e-8, 0.0, -2.19, 8.0, 8)},
                                         eprec);
    success = check_projector("collinear slightly tilted triatomic",
                              collinear_slight_tilt, 4) &&
              success;
  }

  // 5. End-to-end optimization of a z-aligned diatomic with spurious net force
  //    (gtol = 1e-6, well below the 2e-4 net force removed by the projector)
  {
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
  }

  // 6. Non-finite gradient must throw rather than silently force convergence
  {
    Molecule lih(
        {Atom(0.0, 0.0, 0.7625, 3.0, 3), Atom(0.0, 0.0, -2.2875, 1.0, 1)},
        eprec);
    NonFiniteGradientTarget bad_target;
    MolOpt opt(/*maxiter=*/5, 0.1, 1e-5, 1e-4, 1e-4, 1e-6, 1e-5,
               /*print_level=*/0, "bfgs");
    bool threw = false;
    try {
      (void)opt.optimize_app(lih, bad_target);
    } catch (const MadnessException&) {
      threw = true;
    }
    print("non-finite gradient throws MadnessException:",
          threw ? "ok" : "FAILED");
    success = threw && success;
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
