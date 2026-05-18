#!/usr/bin/env python3
"""
Solve the hydrogen atom ground state using MADNESS.

Uses the BSH (bound-state Helmholtz) operator to iteratively solve:
    (-1/2 nabla^2 + V) psi = E psi

where V(r) = -1/r is the nuclear potential.
The exact ground state energy is -0.5 Hartree.


Vectorized callables
--------------------
When pymadness projects a Python callable onto the adaptive grid, MADNESS
evaluates it at many quadrature points.  There are two ways to write the
callable:

**Scalar** (simple, but slow)::

    def V_scalar(r):                          # r has shape (3,)
        return -1.0 / max(np.linalg.norm(r), 1e-6)

**Vectorized** (fast — recommended)::

    def V_vectorized(r):                      # r has shape (npts, 3)
        rr = np.sqrt(np.sum(r**2, axis=1))    # shape (npts,)
        return -1.0 / np.maximum(rr, 1e-6)    # shape (npts,)

The vectorized form receives *all* quadrature points for one box at once as
a numpy array of shape ``(npts, NDIM)`` and must return an array of shape
``(npts,)``.  pymadness auto-detects which convention the callable uses on
the first call:

  - If the callable returns an array of the expected length, it is treated
    as vectorized and all subsequent calls use the batch path.
  - If it raises an exception or returns something else (e.g. a scalar),
    pymadness falls back to calling it once per point.

Either way the GIL is only acquired **once per box** (not once per point),
but the vectorized form is much faster because numpy evaluates the whole
batch in compiled C code rather than looping in Python.  For k=8 in 3D
(512 points per box) the speedup is typically 10–50x.

Writing vectorized callables — tips:

  - Use ``np.sum(r**2, axis=1)`` instead of ``np.dot(r, r)``
    (``np.dot`` on a 2D array does matrix multiplication, not a row-wise
    dot product).
  - Use ``np.maximum(x, c)`` instead of ``max(x, c)`` (the built-in
    ``max`` doesn't work element-wise on arrays).
  - Use ``np.linalg.norm(r, axis=1)`` for row-wise norms.
  - If your function has branches (``if``/``else``), replace them with
    ``np.where(condition, value_if_true, value_if_false)``.
"""

import numpy as np
import pymadness


def run(world):
    # Numerical parameters
    k = 8           # wavelet order (polynomial degree)
    thresh = 1e-6   # truncation threshold
    L = 30.0        # box size [-L, L]^3

    pymadness.FunctionDefaults3D.set_k(k)
    pymadness.FunctionDefaults3D.set_thresh(thresh)
    pymadness.FunctionDefaults3D.set_cubic_cell(-L, L)

    # ----------------------------------------------------------------
    # Nuclear potential: V(r) = -1/|r|, regularized at the origin.
    #
    # Vectorized: r has shape (npts, 3), return shape (npts,).
    # np.maximum replaces the scalar max() so the clamp works on arrays.
    # ----------------------------------------------------------------
    def V_func(r):
        rr = np.sqrt(np.sum(r**2, axis=1))
        return -1.0 / np.maximum(rr, 1e-6)

    # ----------------------------------------------------------------
    # Initial guess: exp(-r^2), a simple Gaussian.
    #
    # Vectorized: np.sum(..., axis=1) computes the row-wise dot product.
    # ----------------------------------------------------------------
    def guess(r):
        return np.exp(-np.sum(r**2, axis=1))

    print("Projecting potential and initial guess...")
    V = pymadness.function_3d(world, V_func)
    psi = pymadness.function_3d(world, guess)

    # Normalize
    norm = psi.norm2()
    psi.scale(1.0 / norm)

    # Initial energy estimate
    E = -0.5
    print(f"Initial energy estimate: {E:.6f}")

    # BSH iteration: solve (T + V) psi = E psi
    # Rearranged: (T - E) psi = -V psi
    # BSH Green's function: G_mu = (-nabla^2 + mu^2)^{-1} with mu = sqrt(-2E)
    # So: psi_new = -2 * G_mu(V * psi)
    for iteration in range(30):
        mu = np.sqrt(-2.0 * E)
        bsh = pymadness.BSHOperator3D(world, mu, 1e-4, thresh)

        Vpsi = V * psi
        Vpsi.truncate()

        # new_psi = -2 * G(V*psi)
        new_psi = pymadness.apply(bsh, Vpsi)
        new_psi.scale(-2.0)
        new_psi.truncate()

        # Compute residual norm for convergence check
        # The BSH residual: r = (T + V - E) psi = new_psi - psi (up to normalization)
        rnorm = (new_psi - psi).norm2()

        # Normalize
        norm = new_psi.norm2()
        new_psi.scale(1.0 / norm)

        # Update energy via expectation value: E = <psi|T+V|psi>
        # Use: <T> = <psi|T|psi> = -0.5 <psi|nabla^2|psi>
        # Or equivalently: E = <psi|V|psi> + <psi|T|psi>
        # Simpler: E = (eps * <psi_old|psi_new>) / <psi_new|psi_new>
        # Best: use the Rayleigh quotient with kinetic energy from derivatives
        Vpsi = V * new_psi
        PE = new_psi.inner(Vpsi)

        # Kinetic energy via gradient: T = 0.5 * sum_i ||d psi/dx_i||^2
        grad = pymadness.gradient(world, new_psi)
        KE = 0.5 * sum(g.inner(g) for g in grad)
        del grad  # free memory

        E_new = PE + KE

        print(f"  iter {iteration:3d}  E = {E_new:14.10f}  "
              f"dE = {E_new - E:10.2e}  residual = {rnorm:.2e}")

        if rnorm < 5 * thresh:
            print("\nConverged!")
            break

        E = E_new
        psi = new_psi

    print(f"\nFinal energy:  {E:14.10f} Hartree")
    print(f"Exact energy:  {-0.5:14.10f} Hartree")
    print(f"Error:         {E - (-0.5):14.10f} Hartree")


def main():
    with pymadness.World(quiet=True) as world:
        run(world)
    # All MADNESS objects (V, psi, bsh, etc.) go out of scope when
    # run() returns, BEFORE the World context manager calls finalize().


if __name__ == "__main__":
    main()
