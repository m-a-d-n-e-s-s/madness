# pymadness

Python bindings for [MADNESS](https://github.com/m-a-d-n-e-s-s/madness)
(Multiresolution Adaptive Numerical Environment for Scientific Simulation).

`pymadness` exposes `Function<T,NDIM>`, integral operators, and derivatives to
Python, letting you solve PDEs and eigenvalue problems interactively.

## Building

From the MADNESS build directory:

```bash
cmake .. -DENABLE_PYTHON=ON
cmake --build . --target _pymadness
```

pybind11 is fetched automatically if not found on the system.

## Running

```bash
cd build
PYTHONPATH=src/pymadness python3 ../src/pymadness/examples/hydrogen_atom.py
```

## Quick start

```python
import numpy as np
import pymadness

def run(world):
    pymadness.FunctionDefaults3D.set_k(8)
    pymadness.FunctionDefaults3D.set_thresh(1e-6)
    pymadness.FunctionDefaults3D.set_cubic_cell(-10.0, 10.0)

    f = pymadness.function_3d(world, lambda r: np.exp(-np.dot(r, r)))
    print(f"norm = {f.norm2():.6f}")
    print(f"integral = {f.trace():.6f}")

    grad = pymadness.gradient(world, f)
    print(f"|grad f| components: {[g.norm2() for g in grad]}")

def main():
    with pymadness.World() as world:
        run(world)

main()
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

Static configuration for each dimension (1D, 2D, 3D):

```python
pymadness.FunctionDefaults3D.set_k(8)            # wavelet order
pymadness.FunctionDefaults3D.set_thresh(1e-6)     # truncation threshold
pymadness.FunctionDefaults3D.set_cubic_cell(-L, L) # simulation domain
```

### Creating functions

```python
f = pymadness.function_3d(world, lambda r: np.exp(-np.dot(r, r)))
f = pymadness.function_1d(world, lambda r: np.sin(r[0]), k=10, thresh=1e-8)
```

The callable receives a numpy array of shape `(NDIM,)` and must return a scalar.

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

### Norms and integrals

```python
f.norm2()                    # L2 norm
f.trace()                    # integral over domain
pymadness.inner(f, g)       # inner product <f|g>
f(np.array([0.0, 0.0, 0.0]))  # point evaluation
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

## Examples

- `examples/hydrogen_atom.py` — BSH iteration for the hydrogen ground state
- `examples/harmonic_oscillator_1d.py` — 1D quantum harmonic oscillator
- `examples/derivatives.py` — gradient and Laplacian of a 3D Gaussian
