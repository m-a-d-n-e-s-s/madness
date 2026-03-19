#!/usr/bin/env python3
"""
Demonstrate differentiation and operator application in MADNESS.

Creates a 3D Gaussian function and computes its gradient,
Laplacian, and verifies against analytical results.
"""

import numpy as np
import pymadness


def run(world):
    pymadness.FunctionDefaults3D.set_k(10)
    pymadness.FunctionDefaults3D.set_thresh(1e-8)
    pymadness.FunctionDefaults3D.set_cubic_cell(-10.0, 10.0)

    alpha = 1.0

    # f(r) = exp(-alpha * r^2)
    f = pymadness.function_3d(world,
        lambda r: np.exp(-alpha * np.dot(r, r)))

    print(f"f: norm={f.norm2():.10f}, integral={f.trace():.10f}")
    print(f"   exact integral = (pi/alpha)^(3/2) = {(np.pi/alpha)**1.5:.10f}")

    # Gradient: df/dx_i = -2*alpha*x_i * f
    grad = pymadness.gradient(world, f)
    for i, label in enumerate(["x", "y", "z"]):
        print(f"  df/d{label}: norm={grad[i].norm2():.8f}")

    # Analytical df/dx: -2*alpha*x * exp(-alpha*r^2)
    dfdx_exact = pymadness.function_3d(world,
        lambda r: -2.0 * alpha * r[0] * np.exp(-alpha * np.dot(r, r)))

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
    lap_exact = pymadness.function_3d(world,
        lambda r: (4*alpha**2 * np.dot(r, r) - 6*alpha) *
                   np.exp(-alpha * np.dot(r, r)))

    error = (laplacian - lap_exact).norm2()
    print(f"  |Laplacian - exact| = {error:.2e}")

    print("\nDone.")


def main():
    with pymadness.World(quiet=True) as world:
        run(world)


if __name__ == "__main__":
    main()
