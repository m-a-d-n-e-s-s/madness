# pymadness

Python bindings for [MADNESS](https://github.com/m-a-d-n-e-s-s/madness)
(Multiresolution Adaptive Numerical Environment for Scientific Simulation).

`pymadness` exposes `Function<T,NDIM>` (1D–6D), integral operators, and
derivatives to Python, letting you solve PDEs and eigenvalue problems
interactively.

## Installation

pymadness is built as part of MADNESS using CMake.

### Prerequisites

- C++17 compiler (GCC >= 7, Clang >= 5, Apple Clang >= 11)
- CMake >= 3.16
- Python >= 3.8 with development headers
- NumPy
- An MPI implementation (Open MPI, MPICH, etc.)
- LAPACK/BLAS

Optional:

- **matplotlib** — for plotting helpers (`pip install matplotlib`)
- **plotly** — for interactive 3D surface plots (`pip install plotly`)
- **numba** — for zero-overhead `@cfunc` projection (`pip install numba`)

### Build from source

```bash
# Clone MADNESS (if you haven't already)
git clone https://github.com/m-a-d-n-e-s-s/madness.git
cd madness

# Configure with Python bindings enabled
cmake -B build -DENABLE_PYTHON=ON

# Build the extension module
cmake --build build --target _pymadness

# Verify the build
cd build/src/pymadness
PYTHONPATH=. python -c "import pymadness; print('pymadness OK')"
```

pybind11 is fetched automatically via CMake FetchContent if not already
installed.

### Making pymadness importable

After building, the extension module (`_pymadness.*.so`) and the Python
package (`pymadness/`) live under `build/src/pymadness/`.  Choose one of:

**Option A — Set `PYTHONPATH` (simplest for development)**

```bash
export PYTHONPATH=/path/to/madness/build/src/pymadness:$PYTHONPATH
```

Add this to your shell profile (`~/.bashrc`, `~/.zshrc`) to make it
persistent across sessions.

**Option B — Install into a prefix**

```bash
cmake --install build --component python --prefix /path/to/install
# Then add /path/to/install/lib/pymadness to PYTHONPATH.
```

**Option C — Symlink into site-packages**

```bash
SITE=$(python -c "import site; print(site.getsitepackages()[0])")
ln -s /path/to/madness/build/src/pymadness/pymadness "$SITE/pymadness"
ln -s /path/to/madness/build/src/pymadness/_pymadness*.so "$SITE/"
```

### Running the tests

```bash
cd build/src/pymadness
PYTHONPATH=. python -m pytest ../../src/pymadness/tests/ -v
```

## Quick start

```python
import numpy as np
import pymadness

def run(world):
    pymadness.FunctionDefaults3D.set_k(8)
    pymadness.FunctionDefaults3D.set_thresh(1e-6)
    pymadness.FunctionDefaults3D.set_cubic_cell(-10.0, 10.0)

    # Vectorized callable: r has shape (npts, 3), return shape (npts,)
    f = pymadness.function_3d(world,
        lambda r: np.exp(-np.sum(r**2, axis=1)))

    print(f"norm     = {f.norm2():.6f}")
    print(f"integral = {f.trace():.6f}")

    grad = pymadness.gradient(world, f)
    print(f"|grad f| components: {[g.norm2() for g in grad]}")

def main():
    with pymadness.World() as world:
        run(world)

main()
```

## Vectorized callables (recommended)

For best performance, write callables that accept a batch of points as a
NumPy array of shape `(npts, NDIM)` and return an array of shape `(npts,)`:

```python
# Scalar (simple, slower — falls back to per-point Python loop):
def f_scalar(r):             # r shape: (3,)
    return np.exp(-np.dot(r, r))

# Vectorized (fast — NumPy evaluates the whole batch in C):
def f_vectorized(r):         # r shape: (npts, 3)
    return np.exp(-np.sum(r**2, axis=1))
```

pymadness auto-detects which convention the callable uses on the first call.
Both are passed to `function_3d()` the same way — no flag is needed.

A Gaussian centered at an offset `r0`:

```python
r0 = np.array([1.0, 2.0, 0.5])
alpha = 3.0

def gaussian_offset(r):         # r shape: (npts, 3)
    dr = r - r0                  # broadcasting: (npts, 3) - (3,)
    return np.exp(-alpha * np.sum(dr**2, axis=1))

f = pymadness.function_3d(world, gaussian_offset)
```

Tips for writing vectorized callables:

- Use `np.sum(r**2, axis=1)` instead of `np.dot(r, r)` (the latter does
  matrix multiplication on 2D arrays).
- Use `r - r0` with a shape `(NDIM,)` array for offsets — NumPy broadcasts
  automatically.
- Use `np.maximum(x, c)` instead of the built-in `max(x, c)`.
- Use `np.where(cond, a, b)` instead of `if`/`else` branches.
- Use `r[:, 0]`, `r[:, 1]`, `r[:, 2]` to extract coordinate columns.

For maximum performance, compile the function with `numba.cfunc` to bypass
Python entirely:

```python
import numba, math

@numba.cfunc("float64(CPointer(float64))")
def gaussian_c(r_ptr):
    r = numba.carray(r_ptr, (3,))
    return math.exp(-(r[0]**2 + r[1]**2 + r[2]**2))

f = pymadness.function_3d_cfunc(world, gaussian_c, k=8, thresh=1e-6)
```

## Grid evaluation and plotting

Evaluate a function on a regular grid with `eval_cube()`, or use the
plotting helpers from `pymadness.plotting`:

```python
import numpy as np
from pymadness.plotting import plot_1d, plot_line_cut, plot_2d_slice

# 1D function
plot_1d(f1d, lo=-5, hi=5, npt=200, show=True)

# Line cut through a 3D function along the x-axis (y=z=0)
plot_line_cut(f3d, axis=0, npt=200, show=True)

# 2D color map of a 3D function in the z=0 plane
plot_2d_slice(f3d, fixed_axis=2, fixed_value=0.0, npt=200, show=True)

# Raw grid data (returns a NumPy array)
cell = np.array([[-5.0, 5.0], [-5.0, 5.0], [-5.0, 5.0]])
vals = f3d.eval_cube(cell, [100, 100, 100])  # shape (100, 100, 100)
```

### Interactive 2D plots

The interactive plot functions `iplot_2d` and `iplot_2d_slice` re-evaluate
the function whenever you zoom or pan, keeping the number of grid points
(resolution) constant.  This means zooming in reveals finer detail rather
than pixelating.

For 3D function slices, you can also change the viewing plane and slice
position interactively:

- Press **x**, **y**, or **z** to switch which axis is held fixed (the
  slice plane).
- **Scroll** the mouse wheel to move the slice position along the fixed
  axis.

```python
from pymadness.plotting import iplot_2d, iplot_2d_slice

# Interactive 2D function plot — zoom/pan re-evaluates at npt resolution
iplot_2d(f2d, npt=200)

# Interactive slice of a 3D function — zoom, pan, change plane & position
iplot_2d_slice(f3d, fixed_axis=2, fixed_value=0.0, npt=200)
```

**Jupyter notebooks:** Use `%matplotlib widget` (requires
`pip install ipympl`) for interactive zoom/pan inside notebooks.
With `%matplotlib inline`, the static plot functions (`plot_2d_slice`, etc.)
are more appropriate.

### 3D surface plots (plotly)

Display a 2D function (or a 2D slice of a 3D function) as an interactive
3D surface.  Uses plotly's WebGL rendering for fast, smooth rotate/zoom/pan
— works in Jupyter notebooks and standalone scripts (opens a browser).

```bash
pip install plotly
```

```python
from pymadness.plotting import plot_surface

# 3D surface of a Function2D
plot_surface(f2d, npt=100)

# 3D surface of a z=0 slice through a Function3D
plot_surface(f3d, fixed_axis=2, fixed_value=0.0, npt=100)

# Zoom into a region and clamp the function-value axis
plot_surface(f3d, npt=100, xrange=[-2, 2], yrange=[-2, 2], zrange=[-0.5, 1.0])

# Multiple functions overlaid
plot_surface([f1, f2, f3], npt=150, fixed_value=0.1,
             colorscale="Viridis", opacity=0.5)

# Per-function colorscales and labels
plot_surface([rho, V_nuc], npt=100,
             colorscale=["RdBu_r", "Greens"],
             opacity=[0.9, 0.5],
             labels=["density", "nuclear potential"])

# To get the plotly Figure object for further customization:
fig = plot_surface(f, npt=100, show=False)
fig.update_layout(...)
fig.show(renderer="notebook")
```

**Axis ranges:**
- `xrange`, `yrange` — set the spatial evaluation domain per axis.  The
  function is re-evaluated on this region at full `npt` resolution, so
  zooming in reveals finer detail.  `lo`/`hi` set both axes at once.
- `zrange` — clamps the displayed function-value axis (visual only).

When multiple functions are passed, each gets a distinct colorscale
(auto-assigned or explicitly specified) and its own colorbar.

No `%matplotlib widget` or extra Jupyter extensions needed — plotly
renders inline automatically.

## FunctionFactory

For fine-grained control over function construction, use the chainable
`FunctionFactory`:

```python
f = (pymadness.FunctionFactory3D(world)
     .functor(lambda r: np.exp(-np.sum(r**2, axis=1)))
     .k(10)
     .thresh(1e-8)
     .initial_level(3)
     .noautorefine()
     .create())
```

Available methods: `functor()`, `k()`, `thresh()`, `initial_level()`,
`max_refine_level()`, `refine()`, `norefine()`, `autorefine()`,
`noautorefine()`, `truncate_mode()`, `truncate_on_project()`,
`notruncate_on_project()`, `fence()`, `nofence()`, `empty()`, `create()`.

## High-dimensional functions (4D–6D)

pymadness supports functions in up to 6 dimensions:

```python
pymadness.FunctionDefaults4D.set_k(5)
pymadness.FunctionDefaults4D.set_thresh(1e-3)
pymadness.FunctionDefaults4D.set_cubic_cell(-5.0, 5.0)

f4d = pymadness.function_4d(world, lambda r: np.exp(-np.sum(r**2, axis=1)))
```

High-dimensional functions use significantly more memory and computation
time.  Use low wavelet order (k=4–6) and loose thresholds (1e-3 to 1e-4)
for exploration.

## Using pymadness in a Jupyter notebook

pymadness works well in Jupyter notebooks for interactive exploration.

### Setup

1. Make sure pymadness is importable from the notebook kernel's Python.
   The easiest way is to set `PYTHONPATH` **before** launching Jupyter:

   ```bash
   export PYTHONPATH=/path/to/madness/build/src/pymadness:$PYTHONPATH
   jupyter lab
   ```

   Alternatively, add the path at the top of the notebook:

   ```python
   import sys
   sys.path.insert(0, "/path/to/madness/build/src/pymadness")
   ```

2. Install matplotlib for inline plots:

   ```bash
   pip install matplotlib
   ```

### Notebook lifecycle

MADNESS must be initialized exactly once per process and finalized before
exit.  In a notebook the kernel persists across cells, so use `World.get()`
/ `World.close()` instead of a `with` block:

```python
# Cell 1 — initialize (run once)
import numpy as np
import pymadness
from pymadness.plotting import *

world = pymadness.World.get(quiet=True)
pymadness.FunctionDefaults3D.set_k(8)
pymadness.FunctionDefaults3D.set_thresh(1e-6)
pymadness.FunctionDefaults3D.set_cubic_cell(-10.0, 10.0)
```

```python
# Cell 2 — work with functions (run as many times as you like)
f = pymadness.function_3d(world, lambda r: np.exp(-np.sum(r**2, axis=1)))
print(f"norm = {f.norm2():.10f}")
```

```python
# Cell 3 — plot
from pymadness.plotting import plot_line_cut
fig, ax = plot_line_cut(f, axis=0, npt=300)
```

```python
# Cell N — finalize (or just restart the kernel)
pymadness.World.close()
```

If you restart the kernel, MADNESS is finalized automatically — call
`World.get()` again in the first cell.

### Inline plots

With `%matplotlib inline` (the default in most Jupyter environments),
plots from `pymadness.plotting` appear inline automatically.  For
interactive zoom/pan, use `%matplotlib widget` (requires `ipympl`:
`pip install ipympl`).

### Keeping MADNESS objects alive

In a notebook, variables persist across cells.  This means you **do not**
need the `run()` / `main()` pattern described below — just make sure to call
`World.close()` (or restart the kernel) after you `del` any remaining
MADNESS objects.

If you assign a new function to the same variable name, the old function is
automatically destroyed:

```python
# Cell A
f = pymadness.function_3d(world, my_func_v1)

# Cell B — f from Cell A is released when this line runs
f = pymadness.function_3d(world, my_func_v2)
```

## Object lifetime rules

MADNESS objects (`Function`, `SeparatedConvolution`, `Derivative`) are backed
by a parallel runtime (`World`) that is initialized/finalized explicitly.
**All MADNESS objects must be destroyed before `World` is finalized.**

If a `Function` or operator is garbage-collected after `finalize()`, the
program will crash (segfault in `DeferredCleanup::destroy`).

### The pattern: separate `run()` from `main()`

```python
def run(world):
    # All MADNESS objects are local variables here.
    f = pymadness.function_3d(world, my_func)
    g = pymadness.apply(op, f)
    # ... do work ...
    # When run() returns, all locals are destroyed.

def main():
    with pymadness.World() as world:
        run(world)
    # finalize() is called AFTER run()'s locals are gone.

main()
```

**Why this works:** Python destroys local variables when a function returns.
By putting all MADNESS work in `run()`, every `Function`, operator, and
derivative is destroyed when `run()` returns — before the `with` block exits
and calls `finalize()`.

**What breaks:**

```python
# BAD: MADNESS objects live in the same scope as the context manager
with pymadness.World() as world:
    f = pymadness.function_3d(world, my_func)
    g = f * f
# finalize() runs here, but f and g are still alive until the
# enclosing scope exits → crash on cleanup
```

### Rules of thumb

1. Always put MADNESS work in a function called from inside the `with` block.
2. If you must use the non-context-manager API (`World.get()` / `World.close()`),
   call `del` on all MADNESS objects and `gc.collect()` before `World.close()`.
3. Never store MADNESS objects as module-level globals or class attributes that
   outlive `World`.

## API reference

### World

| Function | Description |
|----------|-------------|
| `pymadness.World()` | Context manager — initializes and finalizes the runtime |
| `pymadness.World.get(quiet=True)` | Initialize without context manager |
| `pymadness.World.close()` | Finalize the runtime |

### FunctionDefaults

Static configuration for each dimension (1D–6D):

```python
pymadness.FunctionDefaults3D.set_k(8)            # wavelet order
pymadness.FunctionDefaults3D.set_thresh(1e-6)     # truncation threshold
pymadness.FunctionDefaults3D.set_cubic_cell(-L, L) # simulation domain
```

### Creating functions

```python
# From a vectorized callable (recommended)
f = pymadness.function_3d(world, lambda r: np.exp(-np.sum(r**2, axis=1)))

# With explicit parameters
f = pymadness.function_1d(world, my_func, k=10, thresh=1e-8)

# From a numba @cfunc (zero overhead)
f = pymadness.function_3d_cfunc(world, my_cfunc)

# Via FunctionFactory (fine-grained control)
f = pymadness.FunctionFactory3D(world).functor(my_func).k(10).create()
```

### Arithmetic

```python
h = f + g          # addition
h = f - g          # subtraction
h = f * g          # pointwise multiplication
h = 2.0 * f        # scalar multiplication
f += g              # in-place addition
f.scale(alpha)      # in-place scaling
f.truncate()        # discard small coefficients
```

### Pointwise math functions

```python
g = pymadness.exp(f)     # g(x) = exp(f(x))
g = pymadness.log(f)     # g(x) = log(f(x))
g = pymadness.sqrt(f)    # g(x) = sqrt(f(x))

# In-place variants (modify f directly):
f.unaryop("exp")
f.unaryop("log")
f.unaryop("sqrt")
f.unaryop("abs")
```

### Norms and integrals

```python
f.norm2()                         # L2 norm
f.trace()                         # integral over domain
pymadness.inner(f, g)            # inner product <f|g>
f(np.array([0.0, 0.0, 0.0]))    # point evaluation
```

### Grid evaluation

```python
cell = np.array([[-5, 5], [-5, 5], [-5, 5]])  # (NDIM, 2)
vals = f.eval_cube(cell, [100, 100, 100])       # returns numpy array
```

### Derivatives

```python
D = pymadness.Derivative3D(world, axis)  # axis: 0=x, 1=y, 2=z
df = D(f)

# or as a one-liner:
df = pymadness.diff(world, f, axis=0)

# gradient (3D only):
grad = pymadness.gradient(world, f)  # returns [df/dx, df/dy, df/dz]
```

### Integral operators

```python
# BSH (bound-state Helmholtz): exp(-mu*r)/(4*pi*r)
bsh = pymadness.BSHOperator3D(world, mu, lo, eps)
result = pymadness.apply(bsh, f)

# Coulomb: 1/(4*pi*r)
coulomb = pymadness.CoulombOperator(world, lo, eps)
result = pymadness.apply(coulomb, f)

# Slater: exp(-mu*r)
slater = pymadness.SlaterOperator(world, mu, lo, eps)
```

Parameters: `mu` = exponent, `lo` = smallest length scale, `eps` = precision.

### Save / Load

```python
f.save("function.mad")
g = pymadness.Function3D.load(world, "function.mad")
```

### Plotting

All plotting functions live in `pymadness.plotting`:

```python
from pymadness.plotting import (
    plot_1d, plot_2d, plot_line_cut, plot_2d_slice,  # matplotlib, static
    iplot_2d, iplot_2d_slice,                         # matplotlib, interactive
    plot_surface,                                      # plotly, 3D surface
)
```

| Function | Backend | Description |
|----------|---------|-------------|
| `plot_1d(f, ...)` | matplotlib | Plot a 1D function |
| `plot_2d(f, ...)` | matplotlib | 2D color map of a Function2D |
| `plot_line_cut(f, axis, ...)` | matplotlib | 1D line cut through a 2D/3D function |
| `plot_2d_slice(f, fixed_axis, ...)` | matplotlib | 2D color map slice of a 3D function |
| `iplot_2d(f, ...)` | matplotlib | Interactive 2D plot (re-evaluates on zoom) |
| `iplot_2d_slice(f, ...)` | matplotlib | Interactive 2D slice (zoom + scroll/key to change slice) |
| `plot_surface(f, ...)` | plotly | Interactive 3D surface (rotate/zoom/pan) |

Common parameters:

| Parameter | Description |
|-----------|-------------|
| `npt` | Grid points per axis (constant across zoom levels) |
| `lo`, `hi` | Spatial range for all axes (default: simulation cell) |
| `xrange`, `yrange` | Per-axis spatial range, overrides `lo`/`hi` (`plot_surface` only) |
| `zrange` | Function-value display range (`plot_surface` only) |
| `fixed_axis` | Axis to hold constant for slices (0=x, 1=y, 2=z) |
| `fixed_value` | Value along the fixed axis |
| `colorscale` | Plotly colorscale name or list (`plot_surface` only) |
| `opacity` | Surface opacity 0–1, or list (`plot_surface` only) |
| `labels` | Trace names for legend (`plot_surface` only) |
| `show` | Display immediately (default: True) |
| `ax` | Existing matplotlib Axes to plot into (matplotlib functions only) |

## Examples

- `examples/hydrogen_atom.py` — BSH iteration for the hydrogen ground state
- `examples/harmonic_oscillator_1d.py` — 1D quantum harmonic oscillator with KAIN
- `examples/derivatives.py` — gradient and Laplacian of a 3D Gaussian
