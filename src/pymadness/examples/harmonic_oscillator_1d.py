#!/usr/bin/env python3
"""
Solve the 1D quantum harmonic oscillator using MADNESS.

    H psi = (-1/2 d^2/dx^2 + 1/2 x^2) psi = E psi

Ground state energy is 0.5 (in natural units).

This example demonstrates:
  - Potential shift to make the effective energy negative (required for BSH)
  - KAIN (Krylov Accelerated Inexact Newton) subspace solver for convergence
  - Same algorithm as the C++ example src/examples/3dharmonic.cc

The BSH Green's function only works for bound states (negative energy).
Since V(x)=0.5*x^2 has E_0=0.5 > 0, we shift the potential by a constant
DELTA to make the effective energy negative:

    (T + V - DELTA) psi = (E - DELTA) psi

with E_eff = E - DELTA < 0. After convergence, E = E_eff + DELTA.
"""

import numpy as np
import pymadness


class KAIN:
    """Krylov Accelerated Inexact Newton solver.

    Pure Python implementation of the KAIN algorithm from:
        R. J. Harrison, J. Comput. Chem. 25, 328-334 (2004).

    Works with any objects that support inner(), +, -, and scalar *.
    """
    def __init__(self, maxsub=10):
        self.maxsub = maxsub
        self.ulist = []
        self.rlist = []
        self.Q = None

    def update(self, u, r):
        """Given current solution u and residual r, return improved solution.

        Args:
            u: current solution (pymadness Function)
            r: residual = u - G(V*u) or similar
        Returns:
            improved solution
        """
        if self.maxsub == 1:
            return u - r

        m = len(self.ulist)
        self.ulist.append(u.copy())
        self.rlist.append(r.copy())

        # Build/extend the Q matrix: Q[i,j] = <u_i | r_j>
        Qnew = np.zeros((m + 1, m + 1))
        if m > 0:
            Qnew[:m, :m] = self.Q
        for i in range(m + 1):
            Qnew[i, m] = pymadness.inner(self.ulist[i], self.rlist[m])
            Qnew[m, i] = pymadness.inner(self.ulist[m], self.rlist[i])
        self.Q = Qnew

        # Solve KAIN subspace equations
        c = self._solve_kain(self.Q)

        # Form new solution: unew = sum_i c[i] * (u[i] - r[i])
        unew = c[0] * (self.ulist[0] - self.rlist[0])
        for i in range(1, m + 1):
            unew += c[i] * (self.ulist[i] - self.rlist[i])

        # Trim subspace if at maximum size
        if len(self.ulist) == self.maxsub:
            self.ulist.pop(0)
            self.rlist.pop(0)
            self.Q = self.Q[1:, 1:].copy()

        return unew

    @staticmethod
    def _solve_kain(Q, rcond=1e-12):
        """Solve the KAIN subspace equations.

        Given Q[i,j] = <u_i | r_j>, solve for coefficients c such that
        the new solution is sum_i c[i] * (u[i] - r[i]).
        """
        m = Q.shape[0]
        # Build the augmented system
        A = np.zeros((m + 1, m))
        b = np.zeros(m + 1)
        A[:m, :m] = Q
        A[m, :] = 1.0
        b[m] = 1.0

        # Solve via least-squares (handles rank-deficiency)
        c, _, _, _ = np.linalg.lstsq(A, b, rcond=rcond)

        # Safety: clamp if coefficients are too large
        if np.max(np.abs(c)) > 1000.0:
            c = np.zeros(m)
            c[-1] = 1.0

        return c


def run(world):
    # Numerical parameters
    k = 10
    thresh = 1e-8
    L = 10.0
    DELTA = 7.0  # potential shift to make effective energy negative

    pymadness.FunctionDefaults1D.set_k(k)
    pymadness.FunctionDefaults1D.set_thresh(thresh)
    pymadness.FunctionDefaults1D.set_cubic_cell(-L, L)

    # Shifted potential: V(x) = 0.5 * x^2 - DELTA
    V = pymadness.function_1d(world, lambda r: 0.5 * r[0]**2 - DELTA)

    # Initial guess: broad Gaussian (deliberately poor to show KAIN convergence)
    psi = pymadness.function_1d(world, lambda r: np.exp(-0.1 * r[0]**2))
    psi.scale(1.0 / psi.norm2())

    # Compute total energy = kinetic + potential (with shift)
    def compute_energy(psi, V):
        D = pymadness.Derivative1D(world, 0)
        dpsi = D(psi)
        kinetic = 0.5 * dpsi.inner(dpsi)
        potential = psi.inner(V * psi)
        return kinetic + potential

    E = compute_energy(psi, V)
    solver = KAIN(maxsub=10)

    print(f"Potential shift DELTA = {DELTA}")
    print(f"{'iter':>4s}  {'energy':>14s}  {'unshifted':>14s}  {'error':>10s}")

    for iteration in range(50):
        # Construct BSH Green's function with shifted energy
        mu = np.sqrt(-2.0 * E)
        bsh = pymadness.BSHOperator1D(world, mu, 1e-4, thresh)

        Vpsi = V * psi
        Vpsi.truncate()

        # Apply Green's function: G * (V*psi)
        # The residual is: r = psi + 2*G*(V*psi)
        # (would be zero for the exact eigenfunction)
        Gpsi = pymadness.apply(bsh, Vpsi)
        residual = psi + 2.0 * Gpsi

        err = residual.norm2()
        print(f"  {iteration:3d}  {E:14.10f}  {E + DELTA:14.10f}  {err:10.2e}")

        if err < 5.0 * thresh:
            print("\nConverged!")
            break

        # KAIN update: better than simple fixed-point iteration psi = psi - r
        psi = solver.update(psi, residual)

        # Normalize
        psi.scale(1.0 / psi.norm2())

        # Update energy
        E = compute_energy(psi, V)

    print(f"\nFinal energy (shifted):   {E:14.10f}")
    print(f"Final energy (unshifted): {E + DELTA:14.10f}")
    print(f"Exact energy:             {0.5:14.10f}")
    print(f"Error:                    {E + DELTA - 0.5:14.10f}")


def main():
    with pymadness.World(quiet=True) as world:
        run(world)


if __name__ == "__main__":
    main()
