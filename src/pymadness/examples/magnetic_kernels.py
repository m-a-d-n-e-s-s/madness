#!/usr/bin/env python3
"""The magnetic Green's-function kernels (madness/mra/magnetic_kernels.h) from Python.

Applies (T_k - E)^{-1} of all four kernels (ft, chain_symmetric, chain_landau, mehler) to a Gaussian on the
field axis (Lz = 0, where all four coincide) and prints the relative differences, which are of the order of
the truncation threshold.  The Landau-gauge kernel acts on e^{+i chi} f and is transformed back.

Usage:  python3 magnetic_kernels.py [--B 2] [--E 0] [--k 8] [--thresh 1e-4] [--L 10]
"""
import argparse
import gc

import numpy as np

import pymadness as mad


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--B", type=float, default=2.0)
    ap.add_argument("--E", type=float, default=0.0)
    ap.add_argument("--k", type=int, default=8)
    ap.add_argument("--thresh", type=float, default=1e-4)
    ap.add_argument("--L", type=float, default=10.0)
    a = ap.parse_args()
    if not mad.HAVE_MAGNETIC_KERNEL:
        raise SystemExit("this pymadness build has no magnetic kernels")

    world = mad.World.get(quiet=True)
    mad.FunctionDefaults3D.set_k(a.k)
    mad.FunctionDefaults3D.set_thresh(a.thresh)
    mad.FunctionDefaults3D.set_cubic_cell(-a.L, a.L)     # the kernels set the 1D/2D defaults themselves

    f = mad.ComplexFunctionFactory3D(world).functor(
        lambda r: np.exp(-np.sum(r * r, axis=1)).astype(complex)).create()

    G = {}
    for kind in mad.magnetic_kernel_names():
        K = mad.make_magnetic_kernel(world, mad.MagneticKernelParameters(kind=kind, B=a.B, nodes=50, tmin=1e-5,
                                                                         nbroaden=6))
        g = K.to_symmetric_gauge(K.apply(K.from_symmetric_gauge(f), a.E))
        G[kind] = g
        print(f"{kind:16s} gauge {K.gauge().name:9s} paramagnetic in kernel {K.absorbs_paramagnetic()!s:5s}"
              f"  |G f| = {g.norm2():.8f}")
    ref = G["chain_symmetric"]
    for kind, g in G.items():
        print(f"  |G_{kind} f - G_chain_symmetric f| / |G f| = {(g - ref).norm2() / ref.norm2():.2e}")
    del G, ref, f, g, K          # MADNESS objects must go before the world
    gc.collect()
    mad.World.close()


if __name__ == "__main__":
    main()
