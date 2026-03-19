"""
pymadness - Python interface to MADNESS
(Multiresolution Adaptive Numerical Environment for Scientific Simulation)

Provides adaptive multiresolution representations of functions and operators
for solving integral and differential equations in many dimensions.

Basic usage::

    import pymadness

    with pymadness.World() as world:
        pymadness.FunctionDefaults3D.set_k(8)
        pymadness.FunctionDefaults3D.set_thresh(1e-6)
        pymadness.FunctionDefaults3D.set_cubic_cell(-10.0, 10.0)

        f = pymadness.function_3d(world, lambda r: r[0]**2 + r[1]**2)
        print(f"norm = {f.norm2()}")


Vectorized callables (recommended)
-----------------------------------
For best performance, write callables that operate on batches of points.
A vectorized callable receives a numpy array of shape ``(npts, NDIM)``
and returns an array of shape ``(npts,)``::

    # Scalar (simple, slow — one Python call per quadrature point per box):
    def f_scalar(r):            # r shape: (3,)
        return np.exp(-np.dot(r, r))

    # Vectorized (fast — one Python call per box, numpy does the loop):
    def f_vectorized(r):        # r shape: (npts, 3)
        return np.exp(-np.sum(r**2, axis=1))

pymadness auto-detects which convention the callable uses.  Both signatures
are passed to ``function_3d`` the same way — no flag needed.

For functions that cannot be vectorized with numpy, the scalar form still
works; pymadness falls back to a per-point Python loop but still acquires
the GIL only once per box.


numba @cfunc (zero overhead)
-----------------------------
For maximum performance, compile the function with ``numba.cfunc`` and pass
it via ``function_3d_cfunc``.  This bypasses Python entirely::

    import numba, math
    @numba.cfunc("float64(CPointer(float64))")
    def gaussian_c(r_ptr):
        r = numba.carray(r_ptr, (3,))
        return math.exp(-(r[0]**2 + r[1]**2 + r[2]**2))

    f = pymadness.function_3d_cfunc(world, gaussian_c, k=8, thresh=1e-6)
"""

from _pymadness import (
    # World and runtime
    World as _RawWorld,
    initialize as _initialize,
    finalize as _finalize,
    is_initialized as _is_initialized,

    # FunctionDefaults
    FunctionDefaults1D,
    FunctionDefaults2D,
    FunctionDefaults3D,

    # Function types
    Function1D,
    Function2D,
    Function3D,
    ComplexFunction3D,

    # Operators
    CoulombOperator,
    BSHOperator3D,
    BSHOperator1D,
    SlaterOperator,
    Derivative1D,
    Derivative2D,
    Derivative3D,
    SeparatedConvolution1D,
    SeparatedConvolution3D,

    # Free functions
    apply,
    inner,
    copy,
    gradient,
    diff as _diff,

    # Tensor
    Tensor,
    tensor_to_numpy,
    numpy_to_tensor,
)


class World:
    """Context manager for the MADNESS parallel runtime.

    Usage::

        with pymadness.World() as world:
            # ... use world ...
        # runtime is finalized automatically

    Or without context manager::

        world = pymadness.World.get(quiet=True)
        # ... use world ...
        pymadness.World.close()
    """

    _instance = None

    def __init__(self, quiet=True):
        self._quiet = quiet
        self._world = None

    def __enter__(self):
        self._world = _initialize(self._quiet)
        return self._world

    def __exit__(self, *args):
        import gc
        gc.collect()  # Destroy MADNESS objects before finalizing
        _finalize()
        self._world = None

    @staticmethod
    def get(quiet=True):
        """Initialize MADNESS and return the World object (no context manager)."""
        return _initialize(quiet)

    @staticmethod
    def close():
        """Finalize the MADNESS runtime."""
        _finalize()


def function_1d(world, f, k=None, thresh=None):
    """Create a 1D function from a Python callable.

    The callable can use either calling convention:
      - Scalar: ``f(r)`` where ``r`` has shape ``(1,)``, returns a float.
      - Vectorized: ``f(r)`` where ``r`` has shape ``(npts, 1)``, returns
        an array of shape ``(npts,)``.  This is much faster because the
        GIL is acquired only once per box instead of once per point.

    Args:
        world: MADNESS World object
        f: Python callable (scalar or vectorized)
        k: wavelet order (None = use FunctionDefaults)
        thresh: truncation threshold (None = use FunctionDefaults)

    Returns:
        Function1D
    """
    return Function1D(world, f,
                      k=k if k is not None else -1,
                      thresh=thresh if thresh is not None else -1.0)


def function_2d(world, f, k=None, thresh=None):
    """Create a 2D function from a Python callable.

    The callable can use either calling convention:
      - Scalar: ``f(r)`` where ``r`` has shape ``(2,)``, returns a float.
      - Vectorized: ``f(r)`` where ``r`` has shape ``(npts, 2)``, returns
        an array of shape ``(npts,)``.  This is much faster because the
        GIL is acquired only once per box instead of once per point.

    Args:
        world: MADNESS World object
        f: Python callable (scalar or vectorized)
        k: wavelet order (None = use FunctionDefaults)
        thresh: truncation threshold (None = use FunctionDefaults)

    Returns:
        Function2D
    """
    return Function2D(world, f,
                      k=k if k is not None else -1,
                      thresh=thresh if thresh is not None else -1.0)


def function_3d(world, f, k=None, thresh=None):
    """Create a 3D function from a Python callable.

    The callable can use either calling convention:
      - Scalar: ``f(r)`` where ``r`` has shape ``(3,)``, returns a float.
      - Vectorized: ``f(r)`` where ``r`` has shape ``(npts, 3)``, returns
        an array of shape ``(npts,)``.  This is much faster because the
        GIL is acquired only once per box instead of once per point.

    Args:
        world: MADNESS World object
        f: Python callable (scalar or vectorized)
        k: wavelet order (None = use FunctionDefaults)
        thresh: truncation threshold (None = use FunctionDefaults)

    Returns:
        Function3D
    """
    return Function3D(world, f,
                      k=k if k is not None else -1,
                      thresh=thresh if thresh is not None else -1.0)


def function_1d_cfunc(world, cfunc, k=None, thresh=None):
    """Create a 1D function from a numba @cfunc.

    The cfunc must have signature ``float64(CPointer(float64))`` (a pointer
    to a double[1] coordinate array).  Pass the compiled cfunc object directly;
    its address is extracted automatically.

    Args:
        world: MADNESS World object
        cfunc: a numba ``@cfunc`` compiled function (with ``.address`` attribute)
        k: wavelet order (None = use FunctionDefaults)
        thresh: truncation threshold (None = use FunctionDefaults)

    Returns:
        Function1D
    """
    return Function1D.from_cfunc(world, cfunc.address,
                                 k=k if k is not None else -1,
                                 thresh=thresh if thresh is not None else -1.0)


def function_2d_cfunc(world, cfunc, k=None, thresh=None):
    """Create a 2D function from a numba @cfunc.

    See :func:`function_1d_cfunc` for details on the cfunc convention.

    Args:
        world: MADNESS World object
        cfunc: a numba ``@cfunc`` compiled function
        k: wavelet order (None = use FunctionDefaults)
        thresh: truncation threshold (None = use FunctionDefaults)

    Returns:
        Function2D
    """
    return Function2D.from_cfunc(world, cfunc.address,
                                 k=k if k is not None else -1,
                                 thresh=thresh if thresh is not None else -1.0)


def function_3d_cfunc(world, cfunc, k=None, thresh=None):
    """Create a 3D function from a numba @cfunc.

    See :func:`function_1d_cfunc` for details on the cfunc convention.

    Args:
        world: MADNESS World object
        cfunc: a numba ``@cfunc`` compiled function
        k: wavelet order (None = use FunctionDefaults)
        thresh: truncation threshold (None = use FunctionDefaults)

    Returns:
        Function3D
    """
    return Function3D.from_cfunc(world, cfunc.address,
                                 k=k if k is not None else -1,
                                 thresh=thresh if thresh is not None else -1.0)


def diff(world, f, axis=0):
    """Compute derivative of a function along a given axis.

    Args:
        world: MADNESS World object
        f: Function1D, Function2D, or Function3D
        axis: direction to differentiate (0=x, 1=y, 2=z)

    Returns:
        Derivative of f along the specified axis
    """
    return _diff(world, f, axis)


__all__ = [
    # Context manager
    "World",
    # Function constructors
    "function_1d", "function_2d", "function_3d",
    "function_1d_cfunc", "function_2d_cfunc", "function_3d_cfunc",
    # Function types
    "Function1D", "Function2D", "Function3D", "ComplexFunction3D",
    # Defaults
    "FunctionDefaults1D", "FunctionDefaults2D", "FunctionDefaults3D",
    # Operators
    "CoulombOperator", "BSHOperator3D", "BSHOperator1D", "SlaterOperator",
    "Derivative1D", "Derivative2D", "Derivative3D",
    "SeparatedConvolution1D", "SeparatedConvolution3D",
    # Free functions
    "apply", "inner", "copy", "gradient", "diff",
    # Tensor
    "Tensor", "tensor_to_numpy", "numpy_to_tensor",
]
