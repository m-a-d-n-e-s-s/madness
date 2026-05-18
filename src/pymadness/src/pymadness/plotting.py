"""
pymadness.plotting - Visualization helpers for MADNESS functions.

Provides convenience functions for plotting 1D, 2D, and 3D MADNESS functions
using matplotlib.  All functions accept a MADNESS Function object and
evaluate it on a regular grid via ``eval_cube()``.

Requires matplotlib (optional dependency of pymadness).
3D surface plots (``plot_surface``) use plotly for fast WebGL rendering
with interactive rotate/zoom (``pip install plotly``).

Quick examples::

    import pymadness
    from pymadness.plotting import plot_1d, plot_2d_slice, plot_surface

    with pymadness.World() as world:
        pymadness.FunctionDefaults3D.set_k(8)
        pymadness.FunctionDefaults3D.set_thresh(1e-6)
        pymadness.FunctionDefaults3D.set_cubic_cell(-10.0, 10.0)

        f = pymadness.function_3d(world, lambda r: np.exp(-np.sum(r**2, axis=1)))

        # 1D line cut along the x-axis
        plot_line_cut(f, axis=0, npt=200, show=True)

        # Static 2D slice in the z=0 plane
        plot_2d_slice(f, fixed_axis=2, fixed_value=0.0, npt=200, show=True)

        # Interactive 3D surface plot of the same slice (plotly, WebGL)
        plot_surface(f, fixed_axis=2, fixed_value=0.0, npt=100)
"""

import numpy as np


def eval_grid_1d(f, lo=None, hi=None, npt=200):
    """Evaluate a 1D function on a uniform grid.

    Args:
        f: Function1D
        lo: lower bound (default: from FunctionDefaults)
        hi: upper bound (default: from FunctionDefaults)
        npt: number of grid points

    Returns:
        (x, values) — 1D numpy arrays
    """
    import pymadness
    if lo is None or hi is None:
        cell = pymadness.FunctionDefaults1D.get_cell()
        cell_np = pymadness.tensor_to_numpy(cell)
        if lo is None:
            lo = cell_np[0, 0]
        if hi is None:
            hi = cell_np[0, 1]
    cell = np.array([[lo, hi]])
    vals = f.eval_cube(cell, [npt])
    x = np.linspace(lo, hi, npt)
    return x, vals


def plot_1d(f, lo=None, hi=None, npt=200, label=None, ax=None, show=False,
            **kwargs):
    """Plot a 1D MADNESS function.

    Args:
        f: Function1D
        lo, hi: plot range (default: simulation cell)
        npt: number of grid points
        label: line label
        ax: matplotlib Axes (created if None)
        show: call plt.show() when done
        **kwargs: passed to ax.plot()

    Returns:
        (fig, ax)
    """
    import matplotlib.pyplot as plt
    x, vals = eval_grid_1d(f, lo, hi, npt)
    if ax is None:
        fig, ax = plt.subplots()
    else:
        fig = ax.get_figure()
    ax.plot(x, vals, label=label, **kwargs)
    ax.set_xlabel("x")
    if label is not None:
        ax.legend()
    if show:
        plt.show()
    return fig, ax


def plot_line_cut(f, axis=0, fixed=None, npt=200, lo=None, hi=None,
                  label=None, ax=None, show=False, **kwargs):
    """Plot a line cut of a 2D or 3D function along one axis.

    For a 3D function with ``axis=0``, this evaluates f(x, y0, z0) where
    y0, z0 are given by ``fixed`` (default: all zeros).

    Args:
        f: Function2D or Function3D
        axis: which axis to sweep (0=x, 1=y, 2=z)
        fixed: values for the other coordinates (default: zeros).
               For 3D with axis=0: fixed=[y0, z0].
        npt: number of grid points along the swept axis
        lo, hi: range for the swept axis (default: simulation cell)
        label: line label
        ax: matplotlib Axes
        show: call plt.show()
        **kwargs: passed to ax.plot()

    Returns:
        (fig, ax)
    """
    import matplotlib.pyplot as plt
    import pymadness

    ndim = _get_ndim(f)
    if fixed is None:
        fixed = [0.0] * (ndim - 1)
    fixed = list(fixed)
    if len(fixed) != ndim - 1:
        raise ValueError(
            f"fixed must have {ndim - 1} values, got {len(fixed)}")

    # Determine range along swept axis
    defaults = _get_defaults(ndim)
    if lo is None or hi is None:
        cell = pymadness.tensor_to_numpy(defaults.get_cell())
        if lo is None:
            lo = cell[axis, 0]
        if hi is None:
            hi = cell[axis, 1]

    # Build the grid: narrow cell with 1 point in each fixed dimension
    cell_arr = np.zeros((ndim, 2))
    npt_list = [1] * ndim
    for i in range(ndim):
        if i == axis:
            cell_arr[i] = [lo, hi]
            npt_list[i] = npt
        else:
            # Insert the fixed value in the correct position
            idx = i if i < axis else i - 1
            v = fixed[idx]
            cell_arr[i] = [v, v]
            npt_list[i] = 1

    vals = f.eval_cube(cell_arr, npt_list)
    vals = vals.squeeze()  # collapse singleton dimensions

    x = np.linspace(lo, hi, npt)
    if ax is None:
        fig, ax = plt.subplots()
    else:
        fig = ax.get_figure()
    ax.plot(x, vals, label=label, **kwargs)
    ax.set_xlabel("xyz"[axis] if axis < 3 else f"dim {axis}")
    if label is not None:
        ax.legend()
    if show:
        plt.show()
    return fig, ax


def plot_2d_slice(f, fixed_axis=2, fixed_value=0.0, npt=200,
                  lo=None, hi=None, ax=None, show=False, **kwargs):
    """Plot a 2D color map of a 3D function in one plane.

    Args:
        f: Function3D
        fixed_axis: which axis to hold constant (0=x, 1=y, 2=z)
        fixed_value: value of the fixed coordinate
        npt: grid points per free axis
        lo, hi: range for the free axes (default: simulation cell)
        ax: matplotlib Axes
        show: call plt.show()
        **kwargs: passed to ax.pcolormesh()

    Returns:
        (fig, ax, mesh)
    """
    import matplotlib.pyplot as plt
    import pymadness

    ndim = _get_ndim(f)
    if ndim < 3:
        raise ValueError("plot_2d_slice requires a 3D (or higher) function")

    cell = pymadness.tensor_to_numpy(_get_defaults(ndim).get_cell())

    # Build cell and npt arrays for eval_cube
    # Only two axes are sampled; any additional non-fixed axes are pinned to
    # their cell midpoint so that eval_cube always returns a 2D result.
    cell_arr = np.zeros((ndim, 2))
    npt_list = [1] * ndim
    free_axes = []
    for i in range(ndim):
        if i == fixed_axis:
            cell_arr[i] = [fixed_value, fixed_value]
            npt_list[i] = 1
        elif len(free_axes) < 2:
            cell_arr[i, 0] = lo if lo is not None else cell[i, 0]
            cell_arr[i, 1] = hi if hi is not None else cell[i, 1]
            npt_list[i] = npt
            free_axes.append(i)
        else:
            # Pin extra axes to their cell midpoint
            mid = (cell[i, 0] + cell[i, 1]) / 2.0
            cell_arr[i] = [mid, mid]
            npt_list[i] = 1

    vals = f.eval_cube(cell_arr, npt_list)
    vals = vals.squeeze()  # shape (npt, npt)

    # Build coordinate arrays for the two free axes
    x = np.linspace(cell_arr[free_axes[0], 0], cell_arr[free_axes[0], 1], npt)
    y = np.linspace(cell_arr[free_axes[1], 0], cell_arr[free_axes[1], 1], npt)

    if ax is None:
        fig, ax = plt.subplots()
    else:
        fig = ax.get_figure()

    mesh = ax.pcolormesh(x, y, vals.T, shading="auto", **kwargs)
    fig.colorbar(mesh, ax=ax)
    labels = "xyz"
    ax.set_xlabel(labels[free_axes[0]] if free_axes[0] < 3 else f"dim {free_axes[0]}")
    ax.set_ylabel(labels[free_axes[1]] if free_axes[1] < 3 else f"dim {free_axes[1]}")
    ax.set_aspect("equal")
    ax.set_title(f"{labels[fixed_axis]}={fixed_value:.2f}" if fixed_axis < 3
                 else f"dim {fixed_axis}={fixed_value:.2f}")
    if show:
        plt.show()
    return fig, ax, mesh


def plot_2d(f, npt=200, lo=None, hi=None, ax=None, show=False, **kwargs):
    """Plot a 2D MADNESS function as a color map.

    Args:
        f: Function2D
        npt: grid points per axis
        lo, hi: range for both axes (default: simulation cell)
        ax: matplotlib Axes
        show: call plt.show()
        **kwargs: passed to ax.pcolormesh()

    Returns:
        (fig, ax, mesh)
    """
    import matplotlib.pyplot as plt
    import pymadness

    cell = pymadness.tensor_to_numpy(pymadness.FunctionDefaults2D.get_cell())
    cell_arr = np.zeros((2, 2))
    for i in range(2):
        cell_arr[i, 0] = lo if lo is not None else cell[i, 0]
        cell_arr[i, 1] = hi if hi is not None else cell[i, 1]

    vals = f.eval_cube(cell_arr, [npt, npt])
    x = np.linspace(cell_arr[0, 0], cell_arr[0, 1], npt)
    y = np.linspace(cell_arr[1, 0], cell_arr[1, 1], npt)

    if ax is None:
        fig, ax = plt.subplots()
    else:
        fig = ax.get_figure()

    mesh = ax.pcolormesh(x, y, vals.T, shading="auto", **kwargs)
    fig.colorbar(mesh, ax=ax)
    ax.set_xlabel("x")
    ax.set_ylabel("y")
    ax.set_aspect("equal")
    if show:
        plt.show()
    return fig, ax, mesh


# --- 3D surface plots (plotly) ---

_DEFAULT_COLORSCALES = [
    "RdBu_r", "Viridis", "Plasma", "Cividis", "Inferno",
    "Turbo", "Hot", "Electric", "Blues", "Greens",
]


def plot_surface(f, fixed_axis=2, fixed_value=0.0, npt=100,
                 lo=None, hi=None, xrange=None, yrange=None, zrange=None,
                 colorscale="RdBu_r", opacity=1.0,
                 labels=None, show=True, **kwargs):
    """Plot one or more MADNESS functions as interactive 3D surfaces using plotly.

    Accepts a single function or a list/tuple of functions.  When multiple
    functions are given, each is rendered as a separate surface with its own
    colorscale (auto-assigned unless ``colorscale`` is a list).

    For Function2D, plots the function directly as a surface.
    For Function3D (or higher), plots a 2D slice as a surface.

    The plot is rendered via WebGL and supports interactive rotate, zoom,
    and pan out of the box — in Jupyter notebooks, standalone scripts, and
    the browser.

    Requires plotly (``pip install plotly``).

    Args:
        f: a MADNESS Function, or a list/tuple of Functions
        fixed_axis: for 3D+ functions, axis to hold constant (0=x, 1=y, 2=z)
        fixed_value: for 3D+ functions, value along the fixed axis
        npt: grid points per free axis
        lo, hi: range for the free axes (default: simulation cell).
            Sets both spatial axes at once.  Use ``xrange``/``yrange`` to
            override individual axes.
        xrange: [lo, hi] for the first free axis (overrides ``lo``/``hi``)
        yrange: [lo, hi] for the second free axis (overrides ``lo``/``hi``)
        zrange: [lo, hi] for the function-value (z) axis display range
        colorscale: plotly colorscale name, or list of names (one per function).
            Default: auto-cycles through distinct colorscales.
        opacity: surface opacity 0.0–1.0, or list of opacities
        labels: list of trace names (one per function) for the legend
        show: display the figure immediately (default: True)
        **kwargs: passed to each ``go.Surface()``

    Returns:
        plotly ``Figure`` object (only when ``show=False``)

    Examples::

        # Single function
        plot_surface(f, npt=100)

        # Zoom into a region and clamp z-axis
        plot_surface(f, npt=100, xrange=[-2, 2], yrange=[-2, 2],
                     zrange=[-0.5, 1.0])

        # Multiple functions overlaid
        plot_surface([f1, f2, f3], npt=150, fixed_value=0.1,
                     colorscale="Viridis", opacity=0.5)

        # Per-function colorscales and labels
        plot_surface([rho, V_nuc], npt=100,
                     colorscale=["RdBu_r", "Greens"],
                     opacity=[0.9, 0.5],
                     labels=["density", "nuclear potential"])
    """
    try:
        import plotly.graph_objects as go
    except ImportError:
        raise ImportError(
            "plotly is required for plot_surface. Install it with:\n"
            "  pip install plotly") from None

    # Normalize to a list of functions
    if not isinstance(f, (list, tuple)):
        funcs = [f]
    else:
        funcs = list(f)
    nfuncs = len(funcs)

    # Normalize per-trace parameters to lists
    if isinstance(colorscale, str):
        if nfuncs == 1:
            colorscales = [colorscale]
        else:
            # Auto-assign distinct colorscales
            colorscales = [_DEFAULT_COLORSCALES[i % len(_DEFAULT_COLORSCALES)]
                           for i in range(nfuncs)]
    else:
        colorscales = list(colorscale)

    if isinstance(opacity, (int, float)):
        opacities = [opacity] * nfuncs
    else:
        opacities = list(opacity)

    if labels is None:
        trace_labels = [f"f{i+1}" if nfuncs > 1 else None for i in range(nfuncs)]
    else:
        trace_labels = list(labels)
        # Pad or truncate so len(trace_labels) always equals nfuncs
        if len(trace_labels) < nfuncs:
            trace_labels += [f"f{i+1}" for i in range(len(trace_labels), nfuncs)]
        elif len(trace_labels) > nfuncs:
            trace_labels = trace_labels[:nfuncs]

    # Resolve per-axis spatial ranges: xrange/yrange override lo/hi
    eff_xlo = xrange[0] if xrange is not None else lo
    eff_xhi = xrange[1] if xrange is not None else hi
    eff_ylo = yrange[0] if yrange is not None else lo
    eff_yhi = yrange[1] if yrange is not None else hi

    fig = go.Figure()
    axis_labels = None

    for i, func in enumerate(funcs):
        ndim = _get_ndim(func)
        if ndim == 2:
            x, y, vals, ax_labels = _eval_2d_data(
                func, npt=npt, lo=lo, hi=hi,
                xlo=eff_xlo, xhi=eff_xhi, ylo=eff_ylo, yhi=eff_yhi)
        elif ndim >= 3:
            x, y, vals, ax_labels = _eval_2d_slice_data(
                func, fixed_axis=fixed_axis, fixed_value=fixed_value,
                npt=npt, lo=lo, hi=hi,
                xlo=eff_xlo, xhi=eff_xhi, ylo=eff_ylo, yhi=eff_yhi)
        else:
            raise ValueError("plot_surface requires 2D or higher functions")

        if axis_labels is None:
            axis_labels = ax_labels

        trace_kwargs = dict(kwargs)
        if trace_labels[i] is not None:
            trace_kwargs["name"] = trace_labels[i]
        # Show individual colorbars only when there are multiple surfaces
        if nfuncs > 1:
            trace_kwargs.setdefault("showscale", True)
            trace_kwargs.setdefault("colorbar", dict(
                x=1.0 + 0.12 * i, len=0.6, title=trace_labels[i]))

        fig.add_trace(go.Surface(
            x=x, y=y, z=vals.T,
            colorscale=colorscales[i % len(colorscales)],
            opacity=opacities[i % len(opacities)],
            **trace_kwargs
        ))

    title = ""
    ndim0 = _get_ndim(funcs[0])
    if ndim0 >= 3:
        al = "xyzuvw"
        title = f"{al[fixed_axis]} = {fixed_value:.2f}"

    if axis_labels is None:
        axis_labels = ["x", "y"]

    xaxis_opts = dict(title=axis_labels[0])
    yaxis_opts = dict(title=axis_labels[1])
    zaxis_opts = dict(title="f")
    if xrange is not None:
        xaxis_opts["range"] = list(xrange)
    if yrange is not None:
        yaxis_opts["range"] = list(yrange)
    if zrange is not None:
        zaxis_opts["range"] = list(zrange)

    scene = dict(
        xaxis=xaxis_opts,
        yaxis=yaxis_opts,
        zaxis=zaxis_opts,
        aspectmode="manual",
        aspectratio=dict(x=1, y=1, z=0.6),
    )

    fig.update_layout(
        title=title,
        scene=scene,
        margin=dict(l=0, r=0, t=40, b=0),
    )

    if not show:
        return fig
    _show_plotly(fig)
    # Return None so Jupyter doesn't auto-display a second copy.
    # Use show=False to get the Figure object back.


# Keep as an alias for explicit slice calls
plot_surface_slice = plot_surface


def _show_plotly(fig):
    """Display a plotly figure, auto-detecting the best renderer."""
    try:
        from IPython.display import display, HTML
        shell = get_ipython().__class__.__name__  # noqa: F821
        if shell == "ZMQInteractiveShell":
            display(HTML(fig.to_html(include_plotlyjs="cdn", full_html=False)))
        else:
            fig.show()
    except (NameError, ImportError):
        fig.show()


# --- Interactive 2D plots ---

class InteractivePlot2D:
    """Interactive 2D color-map plot that re-evaluates the MADNESS function at
    constant resolution when the view changes (zoom / pan).

    For 3D+ functions displayed as slices, the slice plane and position can be
    changed interactively:

    - Press **x**, **y**, or **z** to switch the fixed (sliced) axis.
    - **Scroll** the mouse wheel to move the slice position along the fixed axis.

    The plot title shows the current slice plane and position.

    Args:
        f: A MADNESS Function (2D or 3D+).
        fixed_axis: For 3D+ functions, the axis to hold constant (0/1/2 = x/y/z).
        fixed_value: Initial value along the fixed axis.
        npt: Grid points per free axis (stays constant on zoom).
        lo, hi: Initial range for the free axes (default: simulation cell).
        scroll_step: Fraction of the axis range to move per scroll tick.
        cmap: Matplotlib colormap name.
        **kwargs: Extra keyword arguments passed to ``pcolormesh()``.

    Returns:
        An ``InteractivePlot2D`` instance (keeps references to fig, ax, etc.).
    """

    def __init__(self, f, fixed_axis=2, fixed_value=0.0, npt=200,
                 lo=None, hi=None, scroll_step=0.02, cmap="RdBu_r",
                 **kwargs):
        import matplotlib.pyplot as plt
        import pymadness

        self.f = f
        self.npt = npt
        self.ndim = _get_ndim(f)
        self.scroll_step = scroll_step
        self.cmap = cmap
        self.mesh_kwargs = kwargs

        if self.ndim < 2:
            raise ValueError(
                f"InteractivePlot2D requires a function of dimension >= 2, "
                f"got {self.ndim}D"
            )

        # For 2D functions there is no slice axis
        if self.ndim == 2:
            self.fixed_axis = None
            self.fixed_value = None
        else:
            if fixed_axis < 0 or fixed_axis >= self.ndim:
                raise ValueError(
                    f"fixed_axis={fixed_axis} is out of range for a {self.ndim}D function "
                    f"(must be in [0, {self.ndim - 1}])"
                )
            self.fixed_axis = fixed_axis
            self.fixed_value = fixed_value

        # Determine free axes and initial ranges from simulation cell
        defaults = _get_defaults(self.ndim)
        cell = pymadness.tensor_to_numpy(defaults.get_cell())
        self.cell = cell  # full simulation cell, for scroll bounds

        if self.ndim == 2:
            self.free_axes = [0, 1]
        else:
            self.free_axes = [i for i in range(self.ndim) if i != self.fixed_axis][:2]

        self.xlim = [lo if lo is not None else float(cell[self.free_axes[0], 0]),
                     hi if hi is not None else float(cell[self.free_axes[0], 1])]
        self.ylim = [lo if lo is not None else float(cell[self.free_axes[1], 0]),
                     hi if hi is not None else float(cell[self.free_axes[1], 1])]

        # Create figure
        self.fig, self.ax = plt.subplots()
        self.mesh = None
        self.colorbar = None
        self._redraw()

        # Connect events — use button_release_event to detect end of zoom/pan
        self.fig.canvas.mpl_connect("button_release_event", self._on_button_release)
        if self.ndim >= 3:
            self.fig.canvas.mpl_connect("scroll_event", self._on_scroll)
            self.fig.canvas.mpl_connect("key_press_event", self._on_key)

        # Display the figure widget in Jupyter
        try:
            from IPython.display import display
            display(self.fig.canvas)
        except ImportError:
            pass

    def _repr_mimebundle_(self, **kwargs):
        """Let Jupyter display the interactive widget when this object is
        the last expression in a cell."""
        if hasattr(self.fig.canvas, '_repr_mimebundle_'):
            return self.fig.canvas._repr_mimebundle_(**kwargs)
        return None

    # -- internal helpers --

    def _evaluate(self):
        """Evaluate the function on the current view."""
        cell_arr = np.zeros((self.ndim, 2))
        npt_list = [1] * self.ndim

        if self.ndim == 2:
            cell_arr[0] = self.xlim
            cell_arr[1] = self.ylim
            npt_list = [self.npt, self.npt]
        else:
            for i in range(self.ndim):
                if i == self.fixed_axis:
                    cell_arr[i] = [self.fixed_value, self.fixed_value]
                    npt_list[i] = 1
                elif i == self.free_axes[0]:
                    cell_arr[i] = self.xlim
                    npt_list[i] = self.npt
                elif i == self.free_axes[1]:
                    cell_arr[i] = self.ylim
                    npt_list[i] = self.npt
                else:
                    mid = (self.cell[i, 0] + self.cell[i, 1]) / 2.0
                    cell_arr[i] = [mid, mid]
                    npt_list[i] = 1

        vals = self.f.eval_cube(cell_arr, npt_list)
        return vals.squeeze()

    def _redraw(self):
        """Full redraw: evaluate the function and update the plot."""
        vals = self._evaluate()
        x = np.linspace(self.xlim[0], self.xlim[1], self.npt)
        y = np.linspace(self.ylim[0], self.ylim[1], self.npt)

        if self.mesh is None:
            # First draw
            self.mesh = self.ax.pcolormesh(x, y, vals.T, shading="auto",
                                           cmap=self.cmap, **self.mesh_kwargs)
            self.colorbar = self.fig.colorbar(self.mesh, ax=self.ax)
        else:
            # Update: remove old mesh, plot new one, update colorbar
            for coll in list(self.ax.collections):
                coll.remove()
            self.mesh = self.ax.pcolormesh(x, y, vals.T, shading="auto",
                                           cmap=self.cmap, **self.mesh_kwargs)
            self.colorbar.update_normal(self.mesh)

        # Force axis limits to match our view (prevent autoscale drift)
        self.ax.set_xlim(self.xlim)
        self.ax.set_ylim(self.ylim)
        self.ax.set_aspect("equal", adjustable="box")

        labels = "xyzuvw"
        self.ax.set_xlabel(labels[self.free_axes[0]])
        self.ax.set_ylabel(labels[self.free_axes[1]])
        self._update_title()
        self.fig.canvas.draw_idle()

    def _update_title(self):
        if self.fixed_axis is not None:
            labels = "xyzuvw"
            self.ax.set_title(
                f"{labels[self.fixed_axis]} = {self.fixed_value:.4f}")

    def _on_button_release(self, event):
        """Re-evaluate after the user finishes a zoom or pan."""
        new_xlim = list(self.ax.get_xlim())
        new_ylim = list(self.ax.get_ylim())
        if new_xlim == self.xlim and new_ylim == self.ylim:
            return
        self.xlim = new_xlim
        self.ylim = new_ylim
        self._redraw()

    def _on_scroll(self, event):
        """Move the slice position on scroll."""
        if self.fixed_axis is None:
            return
        axis_range = (self.cell[self.fixed_axis, 1]
                      - self.cell[self.fixed_axis, 0])
        delta = self.scroll_step * axis_range
        if event.button == "up":
            self.fixed_value += delta
        elif event.button == "down":
            self.fixed_value -= delta
        # Clamp to simulation cell
        self.fixed_value = float(np.clip(
            self.fixed_value,
            self.cell[self.fixed_axis, 0],
            self.cell[self.fixed_axis, 1]))
        self._redraw()

    def _on_key(self, event):
        """Switch slice plane on x / y / z key press."""
        key_map = {"x": 0, "y": 1, "z": 2}
        if event.key not in key_map:
            return
        new_axis = key_map[event.key]
        if new_axis >= self.ndim or new_axis == self.fixed_axis:
            return
        self.fixed_axis = new_axis
        self.fixed_value = 0.0
        self.free_axes = [i for i in range(self.ndim) if i != self.fixed_axis][:2]
        # Reset view to full simulation cell for the new free axes
        self.xlim = [float(self.cell[self.free_axes[0], 0]),
                     float(self.cell[self.free_axes[0], 1])]
        self.ylim = [float(self.cell[self.free_axes[1], 0]),
                     float(self.cell[self.free_axes[1], 1])]
        self._redraw()


def iplot_2d(f, npt=200, lo=None, hi=None, cmap="RdBu_r", **kwargs):
    """Interactive 2D color-map of a Function2D.

    Zoom and pan re-evaluate the function at constant resolution (``npt``
    grid points per axis).

    Args:
        f: Function2D
        npt: grid points per axis (constant across zoom levels)
        lo, hi: initial range for both axes (default: simulation cell)
        cmap: matplotlib colormap name
        **kwargs: passed to ``pcolormesh()``

    Returns:
        InteractivePlot2D instance
    """
    return InteractivePlot2D(f, npt=npt, lo=lo, hi=hi, cmap=cmap, **kwargs)


def iplot_2d_slice(f, fixed_axis=2, fixed_value=0.0, npt=200,
                   lo=None, hi=None, scroll_step=0.02, cmap="RdBu_r",
                   **kwargs):
    """Interactive 2D slice of a 3D (or higher) function.

    Zoom and pan re-evaluate the function at constant resolution.  Additional
    interactive controls:

    - Press **x**, **y**, or **z** to switch the slice plane.
    - **Scroll** the mouse wheel to move the slice position along the
      fixed axis (step size controlled by ``scroll_step``).

    The title displays the current slice axis and position.

    Args:
        f: Function3D (or higher dimension)
        fixed_axis: axis to hold constant (0=x, 1=y, 2=z)
        fixed_value: initial value along the fixed axis
        npt: grid points per free axis (constant across zoom levels)
        lo, hi: initial range for the free axes (default: simulation cell)
        scroll_step: fraction of the axis range per scroll tick (default 0.02)
        cmap: matplotlib colormap name
        **kwargs: passed to ``pcolormesh()``

    Returns:
        InteractivePlot2D instance
    """
    ndim = _get_ndim(f)
    if ndim < 3:
        raise ValueError("iplot_2d_slice requires a 3D (or higher) function")
    return InteractivePlot2D(f, fixed_axis=fixed_axis,
                             fixed_value=fixed_value, npt=npt,
                             lo=lo, hi=hi, scroll_step=scroll_step,
                             cmap=cmap, **kwargs)


# --- Helpers ---

def _get_ndim(f):
    """Infer NDIM from a Function object's type name."""
    if isinstance(f, (list, tuple)):
        raise TypeError(
            f"Expected a single MADNESS Function, got {type(f).__name__}. "
            f"Did you mean to pass a list to plot_surface()?")
    name = type(f).__name__  # e.g. "Function3D"
    for n in range(1, 7):
        if f"{n}D" in name:
            return n
    raise TypeError(f"Cannot determine NDIM from {type(f)}")


def _get_defaults(ndim):
    """Return the FunctionDefaults class for the given NDIM."""
    import pymadness
    return {
        1: pymadness.FunctionDefaults1D,
        2: pymadness.FunctionDefaults2D,
        3: pymadness.FunctionDefaults3D,
        4: pymadness.FunctionDefaults4D,
        5: pymadness.FunctionDefaults5D,
        6: pymadness.FunctionDefaults6D,
    }[ndim]


def _eval_2d_data(f, npt=200, lo=None, hi=None,
                   xlo=None, xhi=None, ylo=None, yhi=None):
    """Evaluate a 2D function on a grid.

    Per-axis bounds (xlo/xhi, ylo/yhi) override lo/hi.

    Returns:
        (x, y, vals, labels) where vals has shape (npt, npt).
    """
    import pymadness
    cell = pymadness.tensor_to_numpy(pymadness.FunctionDefaults2D.get_cell())
    cell_arr = np.zeros((2, 2))
    cell_arr[0, 0] = xlo if xlo is not None else (lo if lo is not None else cell[0, 0])
    cell_arr[0, 1] = xhi if xhi is not None else (hi if hi is not None else cell[0, 1])
    cell_arr[1, 0] = ylo if ylo is not None else (lo if lo is not None else cell[1, 0])
    cell_arr[1, 1] = yhi if yhi is not None else (hi if hi is not None else cell[1, 1])
    vals = f.eval_cube(cell_arr, [npt, npt])
    x = np.linspace(cell_arr[0, 0], cell_arr[0, 1], npt)
    y = np.linspace(cell_arr[1, 0], cell_arr[1, 1], npt)
    return x, y, vals, ["x", "y"]


def _eval_2d_slice_data(f, fixed_axis=2, fixed_value=0.0, npt=200,
                         lo=None, hi=None,
                         xlo=None, xhi=None, ylo=None, yhi=None):
    """Evaluate a 2D slice of a 3D+ function on a grid.

    Per-axis bounds (xlo/xhi, ylo/yhi) override lo/hi for the two free axes.

    Returns:
        (x, y, vals, labels) where vals has shape (npt, npt).
    """
    import pymadness
    ndim = _get_ndim(f)
    defaults = _get_defaults(ndim)
    cell = pymadness.tensor_to_numpy(defaults.get_cell())

    cell_arr = np.zeros((ndim, 2))
    npt_list = [1] * ndim
    free_axes = []
    for i in range(ndim):
        if i == fixed_axis:
            cell_arr[i] = [fixed_value, fixed_value]
            npt_list[i] = 1
        elif len(free_axes) < 2:
            cell_arr[i, 0] = lo if lo is not None else cell[i, 0]
            cell_arr[i, 1] = hi if hi is not None else cell[i, 1]
            npt_list[i] = npt
            free_axes.append(i)
        else:
            # Pin extra axes to their cell midpoint so the result is always 2D
            mid = (cell[i, 0] + cell[i, 1]) / 2.0
            cell_arr[i] = [mid, mid]
            npt_list[i] = 1

    # Apply per-axis overrides to the two free axes
    if xlo is not None:
        cell_arr[free_axes[0], 0] = xlo
    if xhi is not None:
        cell_arr[free_axes[0], 1] = xhi
    if ylo is not None:
        cell_arr[free_axes[1], 0] = ylo
    if yhi is not None:
        cell_arr[free_axes[1], 1] = yhi

    vals = f.eval_cube(cell_arr, npt_list)
    vals = vals.squeeze()

    axis_labels = "xyzuvw"
    x = np.linspace(cell_arr[free_axes[0], 0], cell_arr[free_axes[0], 1], npt)
    y = np.linspace(cell_arr[free_axes[1], 0], cell_arr[free_axes[1], 1], npt)
    labels = [axis_labels[free_axes[0]] if free_axes[0] < 6 else f"dim {free_axes[0]}",
              axis_labels[free_axes[1]] if free_axes[1] < 6 else f"dim {free_axes[1]}"]
    return x, y, vals, labels
