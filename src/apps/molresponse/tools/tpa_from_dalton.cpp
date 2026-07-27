// tpa_from_dalton — pass DALTON's OWN converged response vectors through the
// molresponse 2PA contraction, with ZERO solves. The equation-vs-vector
// discriminator for the multi-system TPA deficit (TPA_SCOPING.md §5k):
//
//   * ground state  : occupied MOs from molden.inp, projected to MRA
//   * X_f           : the EXCITLAB full-RPA eigenvectors from RSPVEC (freq1 = ω_f)
//   * N^A(ω_f/2)    : the X/Y/ZDIPLEN response vectors from RSPVEC at freq1 = ω_f/2
//                     (a .TWO-PHOTON run stores exactly these — verified aug-QZ water)
//
// then S_ab = tpa::tpa_moment(...) — the SAME kernel assemble_tpa uses. If the
// contraction implements DALTON's QRSMO residue (rspvec.F:1746), the output S
// must reproduce DALTON's printed "Two-photon transition tensor S" to MRA
// projection accuracy. If not, the element-wise diff localizes the wrong term;
// if yes, the §5k deficit lives in the *vectors*, and the same harness can mix
// MADNESS/DALTON vectors per channel to find which one carries it.
//
// Convention pinning (printed per root, BEFORE the contraction):
//   * coefficient norms ‖X‖², ‖Y‖², ‖X‖²−‖Y‖² vs the projected function-space
//     norms — import fidelity (virtual MOs are orthonormal, so these must match
//     to projection accuracy);
//   * transition dipole T_a = Σ_i ⟨φ_i|μ_a|x_i+y_i⟩ and the oscillator strength
//     (2/3)ω(cT)² for c ∈ {1,√2,2} — compare to DALTON's printed
//     "Oscillator strength (LENGTH)" to pin the X_f normalization factor;
//   * α_aa(ω_f/2) = -2 Σ_i ⟨x^a_i+y^a_i|μ_a|φ_i⟩ — sign/scale pin for N^A.
// Knobs (--scale-xf/--scale-na/--yflip-xf/--yflip-na) then apply the pinned
// convention without recompiling.
//
// NP=1. Cell fixed at [-200,200]^3 (Gaussian support only; no solver box needed).
//
// Usage:
//   tpa_from_dalton --rspvec=RSPVEC --molden=molden.inp
//                   [--n-occ=N (default: RSPVEC NISH)] [--thresh=1e-4]
//                   [--roots=0,1,2,3] [--scale-xf=1] [--scale-na=1]
//                   [--yflip-xf] [--yflip-na] [--lo=1e-10]

#include "dalton_rspvec.hpp"
#include "dalton_gto.hpp"

#include "../Perturbations.hpp"          // dipole_operator (length gauge)
#include "../ResponseProtocol.hpp"
#include "../kernels/tags.hpp"
#include "../kernels/tpa.hpp"            // tpa_sources / tpa_moment / observables
#include "../solvers/response_state.hpp" // ResponseStateXY<ClosedShell>

#include <madness/mra/mra.h>
#include <madness/world/MADworld.h>

#include <array>
#include <cmath>
#include <cstdio>
#include <memory>
#include <sstream>
#include <string>
#include <vector>

using namespace madness;
using namespace molresponse_v3;

namespace {

// Gaussian-LC functor (same as seed_from_dalton's; local copy keeps the reader
// headers untouched).
class DaltonResponseFunctor : public madness::FunctionFunctorInterface<double, 3> {
  const DaltonMoldenBasis &basis;
  std::vector<double> weights; // AO weights for this function
  std::vector<madness::coord_3d> centers;

public:
  DaltonResponseFunctor(const DaltonMoldenBasis &b, std::vector<double> w)
      : basis(b), weights(std::move(w)) {
    for (const auto &sh : basis.shells) centers.push_back({sh.cx, sh.cy, sh.cz});
  }
  double operator()(const madness::coord_3d &r) const override {
    double val = 0.0;
    double bf[9];
    for (size_t s = 0; s < basis.shells.size(); s++) {
      const auto &sh = basis.shells[s];
      sh.evaluate(r[0], r[1], r[2], bf);
      const int off = basis.ao_offsets[s];
      for (int k = 0; k < sh.n_ao; k++)
        val += weights[static_cast<size_t>(off + k)] * bf[k];
    }
    return val;
  }
  std::vector<madness::coord_3d> special_points() const override { return centers; }
};

using vecfuncT = std::vector<real_function_3d>;

real_function_3d project_weights(World &world, const DaltonMoldenBasis &basis,
                                 std::vector<double> w, double thresh) {
  std::shared_ptr<FunctionFunctorInterface<double, 3>> f =
      std::make_shared<DaltonResponseFunctor>(basis, std::move(w));
  return FunctionFactory<double, 3>(world).functor(f).thresh(thresh)
      .truncate_on_project();
}

// Per-occupied response functions from one flat (occ-outer, vir-inner) block:
//   x_i(r) = sum_a blk[i,a] * phi_{n_occ+a}(r) = sum_mu (C_vir·blkᵀ)[mu,i] chi_mu(r)
vecfuncT project_block(World &world, const DaltonMoldenBasis &basis,
                       const std::vector<double> &C, int n_ao, int n_mo,
                       int n_occ, int n_vir, const std::vector<double> &blk,
                       double thresh, double scale) {
  vecfuncT out;
  for (int i = 0; i < n_occ; ++i) {
    std::vector<double> w(static_cast<size_t>(n_ao), 0.0);
    for (int a = 0; a < n_vir; ++a) {
      const double coef = blk[static_cast<size_t>(i) * n_vir + a] * scale;
      if (coef == 0.0) continue;
      const double *col = &C[static_cast<size_t>(n_ao) * (n_occ + a)];
      for (int mu = 0; mu < n_ao; ++mu) w[static_cast<size_t>(mu)] += col[mu] * coef;
    }
    out.push_back(project_weights(world, basis, std::move(w), thresh));
  }
  truncate(world, out);
  return out;
}

double coef_norm2(const std::vector<double> &v) {
  double s = 0.0;
  for (double x : v) s += x * x;
  return s;
}

std::vector<int> parse_roots_csv(const std::string &s, int nmax) {
  std::vector<int> out;
  std::stringstream ss(s);
  std::string t;
  while (std::getline(ss, t, ','))
    if (!t.empty()) out.push_back(std::stoi(t));
  if (out.empty())
    for (int i = 0; i < nmax; ++i) out.push_back(i);
  return out;
}

} // namespace

int main(int argc, char **argv) {
  World &world = initialize(argc, argv);
  startup(world, argc, argv, true);
  commandlineparser parser(argc, argv);

  if (world.size() != 1) {
    if (world.rank() == 0) print("ERROR: tpa_from_dalton is NP=1 only.");
    finalize();
    return 2;
  }
  if (!parser.key_exists("rspvec") || !parser.key_exists("molden")) {
    print("Usage: tpa_from_dalton --rspvec=RSPVEC --molden=molden.inp");
    print("  [--n-occ=N] [--thresh=1e-4] [--roots=0,1,...] [--scale-xf=1]");
    print("  [--scale-na=1] [--yflip-xf] [--yflip-na] [--lo=1e-10]");
    finalize();
    return 2;
  }

  const std::string rspvec_path = parser.value_raw("rspvec");
  const std::string molden_path = parser.value_raw("molden");
  const double thresh =
      parser.key_exists("thresh") ? std::stod(parser.value("thresh")) : 1e-4;
  const double lo =
      parser.key_exists("lo") ? std::stod(parser.value("lo")) : 1e-10;
  const double scale_xf =
      parser.key_exists("scale-xf") ? std::stod(parser.value("scale-xf")) : 1.0;
  const double scale_na =
      parser.key_exists("scale-na") ? std::stod(parser.value("scale-na")) : 1.0;
  const double yfx = parser.key_exists("yflip-xf") ? -1.0 : 1.0;
  const double yfn = parser.key_exists("yflip-na") ? -1.0 : 1.0;
  // --residue: use the corrected single-residue contraction (X_f in the A slot
  // against V^{bc} built from the two photon responses; beta's b1 term only)
  // instead of the legacy candidate. --prefactor scales the residue moment.
  const bool residue = parser.key_exists("residue");
  const double prefac =
      parser.key_exists("prefactor") ? std::stod(parser.value("prefactor")) : 1.0;
  const double H2EV = 27.211386245988;

  {
    // ---- read DALTON files -------------------------------------------------
    auto [info, entries] = read_rspvec(rspvec_path);
    DaltonMoldenResult molden = read_molden(molden_path);
    const int n_ao = molden.n_ao, n_mo = molden.n_mo;
    const int n_occ = parser.key_exists("n-occ")
                          ? std::stoi(parser.value("n-occ"))
                          : info.nish[0];
    const int n_vir = n_mo - n_occ;

    // Split entries: eigenvectors vs dipole responses.
    std::vector<size_t> exci;
    for (size_t i = 0; i < entries.size(); ++i)
      if (std::string(entries[i].lab1) == "EXCITLAB") exci.push_back(i);
    auto find_dip = [&](int axis, double freq) -> const RspVecEntry * {
      static const char *labs[3] = {"XDIPLEN", "YDIPLEN", "ZDIPLEN"};
      for (const auto &e : entries)
        if (std::string(e.lab1) == labs[axis] && std::abs(e.freq1 - freq) < 1e-6)
          return &e;
      return nullptr;
    };
    const auto roots = parse_roots_csv(
        parser.key_exists("roots") ? parser.value("roots") : "",
        static_cast<int>(exci.size()));

    print("tpa_from_dalton — contraction-only 2PA from DALTON vectors");
    print("  rspvec =", rspvec_path);
    print("  molden =", molden_path);
    print("  n_ao =", n_ao, " n_mo =", n_mo, " n_occ =", n_occ, " n_vir =", n_vir,
          "  eigenvectors:", static_cast<int>(exci.size()));
    print("  thresh =", thresh, " k =", default_k_for_thresh(thresh), " lo =", lo);
    print("  scale_xf =", scale_xf, " scale_na =", scale_na,
          " yflip_xf =", (yfx < 0), " yflip_na =", (yfn < 0));
    print("  contraction =", residue ? "RESIDUE (X_f | V^{bc}[N,N]), b1-only"
                                     : "legacy candidate (beta_abc, C=X_f)",
          " prefactor =", prefac);

    // ---- MRA setup + ground state ------------------------------------------
    Tensor<double> cell(3L, 2L);
    for (int i = 0; i < 3; i++) { cell(i, 0) = -200.0; cell(i, 1) = 200.0; }
    FunctionDefaults<3>::set_cell(cell);
    FunctionDefaults<3>::set_k(default_k_for_thresh(thresh));
    FunctionDefaults<3>::set_thresh(thresh);

    ResponseGroundState g0;
    for (int i = 0; i < n_occ; ++i) {
      std::vector<double> w(molden.mo_coeffs.begin() + static_cast<ptrdiff_t>(i) * n_ao,
                            molden.mo_coeffs.begin() + static_cast<ptrdiff_t>(i + 1) * n_ao);
      g0.amo.push_back(project_weights(world, molden.basis, std::move(w), thresh));
    }
    truncate(world, g0.amo);
    g0.Qa   = QProjector<double, 3>(g0.amo);
    g0.coulop = poperatorT(CoulombOperatorPtr(world, lo, thresh));
    g0.c_xc = 1.0;   // HF
    g0.lo   = lo;

    // import fidelity: MO overlap vs identity
    auto Smo = matrix_inner(world, g0.amo, g0.amo);
    double dev = 0.0;
    for (int i = 0; i < n_occ; ++i)
      for (int j = 0; j < n_occ; ++j)
        dev = std::max(dev, std::abs(Smo(i, j) - (i == j ? 1.0 : 0.0)));
    print("  MO import: max|<phi_i|phi_j> - delta_ij| =", dev);

    std::array<real_function_3d, 3> mu_op{dipole_operator(world, 0),
                                          dipole_operator(world, 1),
                                          dipole_operator(world, 2)};
    std::array<vecfuncT, 3> mu_amo;   // mu_a * phi_i, reused in diagnostics
    for (int a = 0; a < 3; ++a) mu_amo[a] = mul(world, mu_op[a], g0.amo, true);

    // ---- per-root: import, diagnose, contract ------------------------------
    std::vector<int>                     tbl_f;
    std::vector<double>                  tbl_w;
    std::vector<Tensor<double>>         tbl_S;
    std::vector<tpa::Observables>       tbl_o;

    for (int rsel : roots) {
      if (rsel < 0 || rsel >= static_cast<int>(exci.size())) continue;
      const auto &ee = entries[exci[static_cast<size_t>(rsel)]];
      const double wf = ee.freq1, wf_half = 0.5 * wf;

      // X_f
      auto [Xf_x, Xf_y] = split_ov(ee.vec, n_occ, n_vir);
      ResponseStateXY<ClosedShell> Xf;
      Xf.x_alpha = project_block(world, molden.basis, molden.mo_coeffs, n_ao,
                                 n_mo, n_occ, n_vir, Xf_x, thresh, scale_xf);
      if (!Xf_y.empty())
        Xf.y_alpha = project_block(world, molden.basis, molden.mo_coeffs, n_ao,
                                   n_mo, n_occ, n_vir, Xf_y, thresh,
                                   scale_xf * yfx);
      else
        for (int i = 0; i < n_occ; ++i) {
          auto z = copy(Xf.x_alpha[static_cast<size_t>(i)]); z.scale(0.0);
          Xf.y_alpha.push_back(z);
        }

      const double cX2 = coef_norm2(Xf_x), cY2 = coef_norm2(Xf_y);
      const double fX2 = inner(Xf.x_alpha, Xf.x_alpha) / (scale_xf * scale_xf);
      const double fY2 = inner(Xf.y_alpha, Xf.y_alpha) / (scale_xf * scale_xf);
      printf("\n=== root %d  omega = %.6f au (%.3f eV) ===\n", rsel, wf, wf * H2EV);
      printf("  coeff norms: |X|^2 = %.6f  |Y|^2 = %.6f  X2-Y2 = %.6f\n",
             cX2, cY2, cX2 - cY2);
      printf("  func  norms: |x|^2 = %.6f  |y|^2 = %.6f   (import fidelity)\n",
             fX2, fY2);

      // transition dipole + oscillator-strength pin
      printf("  transition dipole / oscillator pin (DALTON prints per-axis "
             "'Oscillator strength (LENGTH)'):\n");
      for (int a = 0; a < 3; ++a) {
        const double T = inner(mu_amo[a], Xf.x_alpha) + inner(mu_amo[a], Xf.y_alpha);
        const double Tm = inner(mu_amo[a], Xf.x_alpha) - inner(mu_amo[a], Xf.y_alpha);
        printf("    axis %c: T(X+Y) = %+.6f  T(X-Y) = %+.6f  osc(c=1,sqrt2,2) = "
               "%.5e  %.5e  %.5e\n",
               "xyz"[a], T, Tm, (2.0 / 3.0) * wf * T * T,
               (2.0 / 3.0) * wf * 2.0 * T * T, (2.0 / 3.0) * wf * 4.0 * T * T);
      }

      // N^A at wf/2
      std::array<ResponseStateXY<ClosedShell>, 3> mu_resp;
      bool ok = true;
      for (int a = 0; a < 3; ++a) {
        const RspVecEntry *de = find_dip(a, wf_half);
        if (!de) {
          printf("  MISSING %cDIPLEN at %.6f — skip root\n", "xyz"[a], wf_half);
          ok = false;
          break;
        }
        auto [Nx, Ny] = split_ov(de->vec, n_occ, n_vir);
        mu_resp[static_cast<size_t>(a)].x_alpha =
            project_block(world, molden.basis, molden.mo_coeffs, n_ao, n_mo,
                          n_occ, n_vir, Nx, thresh, scale_na);
        mu_resp[static_cast<size_t>(a)].y_alpha =
            project_block(world, molden.basis, molden.mo_coeffs, n_ao, n_mo,
                          n_occ, n_vir, Ny, thresh, scale_na * yfn);
      }
      if (!ok) continue;

      // alpha pin: a_raw = sum_i <x+y|mu*phi>; MADNESS convention alpha=-2*a_raw
      printf("  alpha(w/2) pin (diag):");
      for (int a = 0; a < 3; ++a) {
        const double araw = inner(mu_amo[a], mu_resp[static_cast<size_t>(a)].x_alpha) +
                            inner(mu_amo[a], mu_resp[static_cast<size_t>(a)].y_alpha);
        printf("  %c: raw=%+.4f  -2raw=%+.4f", "xyz"[a], araw, -2.0 * araw);
      }
      printf("\n");

      // ---- the contraction (zero solves) ----
      Tensor<double> S;
      if (parser.key_exists("decompose")) {
        // Term-class decomposition of the residue: S_E3 = the pure two-electron
        // part (V^{bc} with ZERO one-electron operators) and S_1e = full - E3
        // (the mu one-electron/property content of V). The correct residue's
        // term weights (Dalton QRSMO: E[3] - B[2] - 2 A[2]) can then be solved
        // per element from Dalton's reference: Dalton = a*S_E3 + b*S_1e.
        real_function_3d zop = copy(mu_op[0]); zop.scale(0.0);
        std::array<real_function_3d, 3> zops{zop, copy(zop), copy(zop)};
        auto S_e3   = tpa::tpa_moment_residue(world, g0, Xf, mu_resp, zops, 1.0);
        auto S_full = tpa::tpa_moment_residue(world, g0, Xf, mu_resp, mu_op, 1.0);
        if (world.rank() == 0) {
          printf("  DECOMPOSE root %d (Sxx Syy Szz Sxy Sxz Syz):\n", rsel);
          printf("    S_E3  : %9.4f %9.4f %9.4f %9.4f %9.4f %9.4f\n",
                 S_e3(0,0), S_e3(1,1), S_e3(2,2), S_e3(0,1), S_e3(0,2), S_e3(1,2));
          printf("    S_1e  : %9.4f %9.4f %9.4f %9.4f %9.4f %9.4f\n",
                 S_full(0,0)-S_e3(0,0), S_full(1,1)-S_e3(1,1),
                 S_full(2,2)-S_e3(2,2), S_full(0,1)-S_e3(0,1),
                 S_full(0,2)-S_e3(0,2), S_full(1,2)-S_e3(1,2));
        }
        S = S_full;
      } else if (residue) {
        S = tpa::tpa_moment_residue(world, g0, Xf, mu_resp, mu_op, prefac);
      } else {
        auto vbc_b = tpa::tpa_sources(world, g0, Xf, mu_resp, mu_op);
        S = tpa::tpa_moment(world, g0, Xf, mu_resp, mu_op, vbc_b);
      }
      const auto obs = tpa::observables(S, wf);
      tbl_f.push_back(rsel); tbl_w.push_back(wf);
      tbl_S.push_back(S);    tbl_o.push_back(obs);
    }

    // ---- Dalton-format tables ----
    printf("\n                  +--------------------------------+\n");
    printf("                  | Two-photon transition tensor S |  (from DALTON vectors)\n");
    printf("                  +--------------------------------+\n");
    printf("      No  Energy(eV)     Sxx     Syy     Szz     Sxy     Sxz     Syz\n");
    for (size_t r = 0; r < tbl_f.size(); ++r)
      printf("     %3d   %8.2f  %8.3f%8.3f%8.3f%8.3f%8.3f%8.3f\n",
             tbl_f[r] + 1, tbl_w[r] * H2EV, tbl_S[r](0, 0), tbl_S[r](1, 1),
             tbl_S[r](2, 2), tbl_S[r](0, 1), tbl_S[r](0, 2), tbl_S[r](1, 2));
    printf("\n      No Energy(eV)      Df          Dg          D_linear    R\n");
    for (size_t r = 0; r < tbl_f.size(); ++r)
      printf("     %3d   %8.2f   %.4e  %.4e  %.4e  %.2f\n",
             tbl_f[r] + 1, tbl_w[r] * H2EV, tbl_o[r].Df, tbl_o[r].Dg,
             tbl_o[r].D_linear, tbl_o[r].R);
  } // all Functions destruct before finalize()

  finalize();
  return 0;
}
