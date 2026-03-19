#!/usr/bin/env python3
"""
Demonstrate differentiation and operator application in MADNESS.

Creates a 3D Gaussian function and computes its gradient,
Laplacian, and verifies against analytical results.

All callables use the vectorized convention: they accept a numpy array
of shape ``(npts, 3)`` and return an array of shape ``(npts,)``.  This
is typically 10–50x faster than the scalar (per-point) convention because
numpy evaluates the whole batch in compiled C code while the GIL is
acquired only once per box.  See ``hydrogen_atom.py`` for a longer
discussion.
"""

import numpy as np
import pymadness


def run(world):
    pymadness.FunctionDefaults3D.set_k(10)
    pymadness.FunctionDefaults3D.set_thresh(1e-8)
    pymadness.FunctionDefaults3D.set_cubic_cell(-10.0, 10.0)

    alpha = 1.0

    # f(r) = exp(-alpha * r^2)
    # Vectorized: r has shape (npts, 3), np.sum(..., axis=1) gives row-wise r^2.
    f = pymadness.function_3d(world,
        lambda r: np.exp(-alpha * np.sum(r**2, axis=1)))

    print(f"f: norm={f.norm2():.10f}, integral={f.trace():.10f}")
    print(f"   exact integral = (pi/alpha)^(3/2) = {(np.pi/alpha)**1.5:.10f}")

    # Gradient: df/dx_i = -2*alpha*x_i * f
    grad = pymadness.gradient(world, f)
    for i, label in enumerate(["x", "y", "z"]):
        print(f"  df/d{label}: norm={grad[i].norm2():.8f}")

    # Analytical df/dx: -2*alpha*x * exp(-alpha*r^2)
    # Vectorized: r[:, 0] extracts the x-coordinates for all points.
    dfdx_exact = pymadness.function_3d(world,
        lambda r: -2.0 * alpha * r[:, 0] * np.exp(-alpha * np.sum(r**2, axis=1)))

    error = (grad[0] - dfdx_exact).norm2()
    print(f"\n  |df/dx - exact|  = {error:.2e}")

    # Laplacian: sum of second derivatives
    laplacian = pymadness.Function3D()
    for axis in range(3):
        D = pymadness.Derivative3D(world, axis)
        d2f = D(grad[axis])
        if not laplacian.is_initialized():
            laplacian = d2f
        else:
            laplacian += d2f

    # Analytical Laplacian: (4*alpha^2 * r^2 - 6*alpha) * f
    # Vectorized: np.sum(r**2, axis=1) replaces np.dot(r, r).
    lap_exact = pymadness.function_3d(world,
        lambda r: (4*alpha**2 * np.sum(r**2, axis=1) - 6*alpha) *
                   np.exp(-alpha * np.sum(r**2, axis=1)))

    error = (laplacian - lap_exact).norm2()
    print(f"  |Laplacian - exact| = {error:.2e}")

    print("\nDone.")


def main():
    with pymadness.World(quiet=True) as world:
        run(world)


if __name__ == "__main__":
    main()
