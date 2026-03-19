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

    Args:
        world: MADNESS World object
        f: callable taking a numpy array of shape (1,) and returning float
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

    Args:
        world: MADNESS World object
        f: callable taking a numpy array of shape (2,) and returning float
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

    Args:
        world: MADNESS World object
        f: callable taking a numpy array of shape (3,) and returning float
        k: wavelet order (None = use FunctionDefaults)
        thresh: truncation threshold (None = use FunctionDefaults)

    Returns:
        Function3D
    """
    return Function3D(world, f,
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
