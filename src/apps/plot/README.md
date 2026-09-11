# Plot Functions from Chemistry Libraries

```bash
plot2cube file=mra_orbital_0
plot2plane file=mra_orbital_1
```

needs files
- mra_orbital_0.00000
- input



## Reading the plot files in Python

`read_plots.py` (numpy; scipy/scikit-image/matplotlib for the figures) reads moldft's `.cube` and `.dx`
files and offers isodensity-band statistics of the electrostatic potential, differences between two
runs on a common surface, orbital overlaps between runs, orbital extent measures, isosurface and slice
figures. `python3 read_plots.py --help` lists the modes; the docstring states the sign and frame conventions.
