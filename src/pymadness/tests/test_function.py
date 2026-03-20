"""Tests for pymadness Function bindings."""

import numpy as np
import pytest

import pymadness


@pytest.fixture(scope="session")
def world():
    """Initialize MADNESS once for the entire test session."""
    w = pymadness.World.get(quiet=True)
    pymadness.FunctionDefaults3D.set_k(6)
    pymadness.FunctionDefaults3D.set_thresh(1e-5)
    pymadness.FunctionDefaults3D.set_cubic_cell(-10.0, 10.0)
    pymadness.FunctionDefaults1D.set_k(6)
    pymadness.FunctionDefaults1D.set_thresh(1e-5)
    pymadness.FunctionDefaults1D.set_cubic_cell(-10.0, 10.0)
    yield w
    pymadness.World.close()


class TestFunction3D:
    """Tests for 3D real function."""

    def test_create_from_callable(self, world):
        f = pymadness.function_3d(world, lambda r: np.exp(-np.dot(r, r)))
        assert f.is_initialized()
        assert f.k() == 6

    def test_norm(self, world):
        """norm of exp(-r^2) should be (pi/2)^(3/4)."""
        f = pymadness.function_3d(world, lambda r: np.exp(-np.dot(r, r)))
        expected = (np.pi / 2.0) ** 0.75  # sqrt(integral exp(-2r^2))
        assert abs(f.norm2() - expected) < 1e-3

    def test_trace(self, world):
        """Integral of exp(-r^2) = pi^(3/2)."""
        f = pymadness.function_3d(world, lambda r: np.exp(-np.dot(r, r)))
        expected = np.pi ** 1.5
        assert abs(f.trace() - expected) < 1e-2

    def test_addition(self, world):
        f = pymadness.function_3d(world, lambda r: np.exp(-np.dot(r, r)))
        g = f + f
        # norm of 2f should be 2 * norm(f)
        assert abs(g.norm2() - 2.0 * f.norm2()) < 1e-2

    def test_subtraction(self, world):
        f = pymadness.function_3d(world, lambda r: np.exp(-np.dot(r, r)))
        g = f - f
        assert g.norm2() < 1e-4

    def test_scalar_multiply(self, world):
        f = pymadness.function_3d(world, lambda r: np.exp(-np.dot(r, r)))
        g = 3.0 * f
        assert abs(g.norm2() - 3.0 * f.norm2()) < 1e-3

    def test_pointwise_multiply(self, world):
        """f*f should equal f^2."""
        f = pymadness.function_3d(world, lambda r: np.exp(-np.dot(r, r)))
        f2 = f * f
        # integral of exp(-2*r^2) = (pi/2)^(3/2)
        expected = (np.pi / 2.0) ** 1.5
        assert abs(f2.trace() - expected) < 1e-2

    def test_inner_product(self, world):
        """<f|f> = ||f||^2."""
        f = pymadness.function_3d(world, lambda r: np.exp(-np.dot(r, r)))
        ip = pymadness.inner(f, f)
        assert abs(ip - f.norm2()**2) < 1e-3

    def test_point_evaluation(self, world):
        f = pymadness.function_3d(world, lambda r: np.exp(-np.dot(r, r)))
        val = f(np.array([0.0, 0.0, 0.0]))
        assert abs(val - 1.0) < 1e-2

    def test_copy(self, world):
        f = pymadness.function_3d(world, lambda r: np.exp(-np.dot(r, r)))
        g = f.copy()
        assert abs(f.norm2() - g.norm2()) < 1e-10

    def test_compress_reconstruct(self, world):
        f = pymadness.function_3d(world, lambda r: np.exp(-np.dot(r, r)))
        f.compress()
        assert f.is_compressed()
        f.reconstruct()
        assert f.is_reconstructed()

    def test_truncate(self, world):
        f = pymadness.function_3d(world, lambda r: np.exp(-np.dot(r, r)))
        norm_before = f.norm2()
        f.truncate()
        assert abs(f.norm2() - norm_before) < 1e-4

    def test_repr(self, world):
        f = pymadness.function_3d(world, lambda r: np.exp(-np.dot(r, r)))
        s = repr(f)
        assert "k=6" in s
        assert "nodes=" in s


class TestFunction1D:
    """Tests for 1D real function."""

    def test_create_and_norm(self, world):
        f = pymadness.function_1d(world, lambda r: np.exp(-r[0]**2))
        # ||exp(-x^2)|| = sqrt(integral exp(-2x^2) dx) = (pi/2)^(1/4)
        expected = (np.pi / 2.0) ** 0.25
        assert abs(f.norm2() - expected) < 1e-3

    def test_derivative(self, world):
        """d/dx exp(-x^2) = -2x exp(-x^2)."""
        f = pymadness.function_1d(world, lambda r: np.exp(-r[0]**2))
        D = pymadness.Derivative1D(world, 0)
        df = D(f)
        exact = pymadness.function_1d(world, lambda r: -2.0 * r[0] * np.exp(-r[0]**2))
        error = (df - exact).norm2()
        assert error < 1e-3


class TestVectorized:
    """Tests for vectorized (batch) callable evaluation."""

    def test_vectorized_3d(self, world):
        """Vectorized callable should give same result as scalar."""
        def gaussian_vectorized(r):
            # r has shape (npts, 3) — use numpy broadcasting
            return np.exp(-np.sum(r**2, axis=1))

        f_vec = pymadness.function_3d(world, gaussian_vectorized)
        expected_norm = (np.pi / 2.0) ** 0.75
        assert abs(f_vec.norm2() - expected_norm) < 1e-3

    def test_vectorized_matches_scalar_3d(self, world):
        """Vectorized and scalar callables should produce the same function."""
        def scalar_fn(r):
            return np.exp(-np.dot(r, r))

        def vectorized_fn(r):
            return np.exp(-np.sum(r**2, axis=1))

        f_scalar = pymadness.function_3d(world, scalar_fn)
        f_vec = pymadness.function_3d(world, vectorized_fn)
        diff = (f_scalar - f_vec).norm2()
        assert diff < 1e-4

    def test_vectorized_1d(self, world):
        """Vectorized callable should work for 1D too."""
        def gaussian_1d_vec(r):
            return np.exp(-r[:, 0]**2)

        f = pymadness.function_1d(world, gaussian_1d_vec)
        expected = (np.pi / 2.0) ** 0.25
        assert abs(f.norm2() - expected) < 1e-3

    def test_scalar_still_works(self, world):
        """Old-style scalar callables must still work through vectorized path."""
        # This callable uses np.dot which only works for 1D input
        f = pymadness.function_3d(world, lambda r: np.exp(-np.dot(r, r)))
        expected = (np.pi / 2.0) ** 0.75
        assert abs(f.norm2() - expected) < 1e-3


class TestCfunc:
    """Tests for numba @cfunc path (optional)."""

    @pytest.fixture(autouse=True)
    def _skip_without_numba(self):
        pytest.importorskip("numba")

    def test_cfunc_3d(self, world):
        """numba @cfunc should give same result as Python callable."""
        import numba
        import ctypes

        @numba.cfunc("float64(CPointer(float64))")
        def gaussian_c(r_ptr):
            r = numba.carray(r_ptr, (3,))
            return np.exp(-(r[0]*r[0] + r[1]*r[1] + r[2]*r[2]))

        f = pymadness.function_3d_cfunc(world, gaussian_c)
        expected = (np.pi / 2.0) ** 0.75
        assert abs(f.norm2() - expected) < 1e-3

    def test_cfunc_1d(self, world):
        """numba @cfunc should work for 1D."""
        import numba

        @numba.cfunc("float64(CPointer(float64))")
        def gaussian_1d_c(r_ptr):
            r = numba.carray(r_ptr, (1,))
            return np.exp(-r[0]*r[0])

        f = pymadness.function_1d_cfunc(world, gaussian_1d_c)
        expected = (np.pi / 2.0) ** 0.25
        assert abs(f.norm2() - expected) < 1e-3


class TestEvalCube:
    """Tests for eval_cube grid evaluation."""

    def test_eval_cube_1d(self, world):
        """eval_cube on a 1D Gaussian should match analytical values."""
        f = pymadness.function_1d(world, lambda r: np.exp(-r[:, 0]**2))
        cell = np.array([[-5.0, 5.0]])
        npt = 50
        vals = f.eval_cube(cell, [npt])
        assert vals.shape == (npt,)
        x = np.linspace(-5.0, 5.0, npt)
        expected = np.exp(-x**2)
        assert np.max(np.abs(vals - expected)) < 1e-3

    def test_eval_cube_3d(self, world):
        """eval_cube on a 3D Gaussian at the origin should be ~1."""
        f = pymadness.function_3d(world, lambda r: np.exp(-np.sum(r**2, axis=1)))
        # Small cube centered at origin
        cell = np.array([[-0.1, 0.1], [-0.1, 0.1], [-0.1, 0.1]])
        vals = f.eval_cube(cell, [3, 3, 3])
        assert vals.shape == (3, 3, 3)
        # Center value should be close to 1.0
        assert abs(vals[1, 1, 1] - 1.0) < 1e-2

    def test_eval_cube_2d(self, world):
        """eval_cube on a 2D function."""
        pymadness.FunctionDefaults2D.set_k(6)
        pymadness.FunctionDefaults2D.set_thresh(1e-5)
        pymadness.FunctionDefaults2D.set_cubic_cell(-10.0, 10.0)
        f = pymadness.function_2d(world, lambda r: np.exp(-np.sum(r**2, axis=1)))
        cell = np.array([[-1.0, 1.0], [-1.0, 1.0]])
        vals = f.eval_cube(cell, [5, 5])
        assert vals.shape == (5, 5)
        # Center should be close to 1.0
        assert abs(vals[2, 2] - 1.0) < 1e-2


class TestFunctionFactory:
    """Tests for FunctionFactory Python bindings."""

    def test_factory_basic(self, world):
        """Build a function through FunctionFactory with chained methods."""
        factory = pymadness.FunctionFactory3D(world)
        factory.functor(lambda r: np.exp(-np.sum(r**2, axis=1)))
        factory.k(6)
        factory.thresh(1e-5)
        f = factory.create()
        assert f.is_initialized()
        expected = (np.pi / 2.0) ** 0.75
        assert abs(f.norm2() - expected) < 1e-3

    def test_factory_chaining(self, world):
        """FunctionFactory methods should be chainable."""
        f = (pymadness.FunctionFactory3D(world)
             .functor(lambda r: np.exp(-np.sum(r**2, axis=1)))
             .k(6)
             .thresh(1e-5)
             .create())
        assert f.is_initialized()

    def test_factory_1d(self, world):
        """FunctionFactory1D should work."""
        f = (pymadness.FunctionFactory1D(world)
             .functor(lambda r: np.exp(-r[:, 0]**2))
             .create())
        expected = (np.pi / 2.0) ** 0.25
        assert abs(f.norm2() - expected) < 1e-3

    def test_factory_norefine(self, world):
        """norefine() should produce a shallower tree."""
        f_refined = (pymadness.FunctionFactory3D(world)
                     .functor(lambda r: np.exp(-np.sum(r**2, axis=1)))
                     .create())
        f_shallow = (pymadness.FunctionFactory3D(world)
                     .functor(lambda r: np.exp(-np.sum(r**2, axis=1)))
                     .norefine()
                     .create())
        assert f_shallow.max_depth() <= f_refined.max_depth()


class TestHigherDimensions:
    """Tests for 4D, 5D, 6D functions."""

    def test_function_4d(self, world):
        """4D Gaussian should have correct norm."""
        pymadness.FunctionDefaults4D.set_k(5)
        pymadness.FunctionDefaults4D.set_thresh(1e-3)
        pymadness.FunctionDefaults4D.set_cubic_cell(-5.0, 5.0)

        f = pymadness.function_4d(world,
            lambda r: np.exp(-np.sum(r**2, axis=1)))
        # ||exp(-r^2)||_4D = (pi/2)^(4/4) = pi/2
        expected = (np.pi / 2.0) ** 1.0
        assert abs(f.norm2() - expected) < 0.1  # loose tol for low k

    def test_defaults_4d(self, world):
        """FunctionDefaults4D should be accessible."""
        pymadness.FunctionDefaults4D.set_k(5)
        assert pymadness.FunctionDefaults4D.get_k() == 5


class TestOperators:
    """Tests for operator application."""

    def test_gradient(self, world):
        """Gradient of a Gaussian should have zero x-component integral."""
        f = pymadness.function_3d(world, lambda r: np.exp(-np.dot(r, r)))
        grad = pymadness.gradient(world, f)
        assert len(grad) == 3
        # Each component of gradient of a symmetric Gaussian integrates to zero
        for i in range(3):
            assert abs(grad[i].trace()) < 1e-4
