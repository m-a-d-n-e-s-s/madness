// ===========================================================================
// test_diagonalize_slots.cpp — pins the slot-identity contract of
// rs::diagonalize (kernels/response_space_ops.hpp).
//
// History: an ascending-omega re-sort ("step 5.5") used to run AFTER the
// diagonal-dominance sort. It exactly inverted the dominance permutation
// (dominance_perm collapsed to identity — the [ROT-SLOTS] REORDERED
// diagnostic in es_solver.hpp could never fire), voided the phase fix, and
// rotated the cluster-unmix blocks off their slots, so near-degenerate
// roots swapped slots/flipped signs every iteration at coarse protocols.
// This test pins the fixed contract:
//
//   1. Slot order is DOMINANCE order: with input slots ordered descending
//      in omega, dominance_perm != identity (REORDERED fires) and
//      omega(slot) tracks the slot's input axis, NOT ascending order.
//   2. Phase fix survives to the output: U(i,i) >= 0 for every slot.
//   3. (omega, U) are consistent eigenpairs: A·u_s = omega(s)·S·u_s.
//   4. Rank-deficient overlap: omega / dominance_perm / omega_ascending
//      are padded to the EXPANDED dimension of U (previously returned
//      short — the caller indexed omega out of bounds).
//
//   test_diagonalize_slots   (MPI=1, fast; PASS/FAIL printed, rc reflects it)
// ===========================================================================

#include "../kernels/response_space_ops.hpp"

#include <madness/world/MADworld.h>

#include <cmath>
#include <cstdio>

using namespace madness;
using namespace molresponse_v3;

static int n_fail = 0;
static void check(bool ok, const char *what) {
  std::printf("  %-64s %s\n", what, ok ? "PASS" : "FAIL");
  if (!ok) ++n_fail;
}

int main(int argc, char **argv) {
  [[maybe_unused]] World &world = initialize(argc, argv);
  {
    const double thr = 1.0e-6;

    // ---- Case 1: input slots descending in omega ------------------------
    // A = [[0.5, eps], [eps, 0.3]], S = I. sygvp returns ascending
    // (0.3-ish first); the dominance sort must swap the columns back so
    // slot 0 keeps its 0.5-axis identity. eps is well below the cluster
    // window so unmixing stays out of the way.
    {
      Tensor<double> A(2, 2), S(2, 2);
      A(0, 0) = 0.5;  A(1, 1) = 0.3;
      A(0, 1) = A(1, 0) = 1.0e-3;
      S(0, 0) = S(1, 1) = 1.0;

      auto r = rs::diagonalize(A, S, thr, /*cluster_factor=*/100.0);

      bool ident = true;
      for (long i = 0; i < 2; ++i)
        if (r.dominance_perm[i] != i) ident = false;
      check(!ident, "descending input: dominance_perm != identity (REORDERED)");
      check(r.omega(0L) > r.omega(1L),
            "descending input: omega stays slot-ordered (not re-sorted)");
      check(std::fabs(r.omega(0L) - 0.5) < 1e-4 &&
            std::fabs(r.omega(1L) - 0.3) < 1e-4,
            "descending input: omega values track slot axes");
      check(r.omega_ascending(0L) < r.omega_ascending(1L),
            "omega_ascending is ascending");
      bool phase_ok = true, pairs_ok = true;
      for (long s = 0; s < 2; ++s) {
        if (r.U(s, s) < 0.0) phase_ok = false;
        // ||A u_s - omega_s S u_s||_inf
        for (long i = 0; i < 2; ++i) {
          double lhs = 0.0, rhs = 0.0;
          for (long j = 0; j < 2; ++j) {
            lhs += A(i, j) * r.U(j, s);
            rhs += S(i, j) * r.U(j, s);
          }
          if (std::fabs(lhs - r.omega(s) * rhs) > 1e-10) pairs_ok = false;
        }
      }
      check(phase_ok, "phase fix survives: U(i,i) >= 0");
      check(pairs_ok, "(omega, U) are consistent generalized eigenpairs");
    }

    // ---- Case 2: rank-deficient overlap ---------------------------------
    // Slot 2 duplicates slot 0 (S singular, one zero singular value). The
    // eigenproblem runs at reduced dim 2 but U is expanded back to 3x3;
    // omega / perm / omega_ascending must come back at dim 3 too.
    {
      Tensor<double> A(3, 3), S(3, 3);
      A(0, 0) = 0.4;  A(1, 1) = 0.6;  A(2, 2) = 0.4;
      A(0, 2) = A(2, 0) = 0.4;   // duplicate axis coupling
      S(0, 0) = S(1, 1) = S(2, 2) = 1.0;
      S(0, 2) = S(2, 0) = 1.0;   // slot 2 == slot 0

      auto r = rs::diagonalize(A, S, thr, /*cluster_factor=*/100.0);

      check(r.U.dim(0) == 3 && r.U.dim(1) == 3, "rank-deficient: U is 3x3");
      check(r.omega.dim(0) == 3,
            "rank-deficient: omega padded to expanded dim (was OOB)");
      check(static_cast<long>(r.dominance_perm.size()) == 3,
            "rank-deficient: dominance_perm padded to expanded dim");
      check(r.omega_ascending.dim(0) == 3,
            "rank-deficient: omega_ascending padded to expanded dim");
      bool finite = true;
      for (long i = 0; i < 3; ++i)
        if (!std::isfinite(r.omega(i))) finite = false;
      check(finite, "rank-deficient: all omega finite (bounded placeholder)");
    }

    std::printf("test_diagonalize_slots: %s\n",
                n_fail == 0 ? "ALL PASS" : "FAILURES");
  }
  finalize();
  return n_fail == 0 ? 0 : 1;
}
