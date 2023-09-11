# Contents

This directory contains various utilities:

- autocorr.mw ... a Maple worksheet to generate the autocorrelation coefficients
- dump2.py with dependencies ... a Python program to generate the two-scale coefficients. Run with `python dump2.py`
- quadrature.py can be easily modified to generate the full-precision Gauss-Legendre coeffs
- MRAMeshOrbitalPlot3D.wl ... a Wolfram language ("Mathematica", for the old school) package for plotting orbitals stored in the Gaussian Cube format (see `madness::plot_cubefile`) along with the mesh stored a JSON format (see `madness::print_tree_jsonfile`).

## MRAMeshOrbitalPlot3D.wl

- if you have [https://www.wolfram.com/wolframscript/](WolframScript) installed (yes by default on Windows/Linux; on a Mac need to [https://www.wolfram.com/wolframscript/](download/install manually)) to plot from command-line or shell script. See a quick example in `h2-no1.wsl`; to try out execute `./h2-no1.wsl > h2-no1.pdf` 
- more elaborate plotting is best done from a Mathematica notebook. Here's a brief example:
```Wolfram
<< "/path/to/madness/source/dir/src/madness/mra/tools/MRAMeshOrbitalPlot3D.wl"

(* this assumes that data.cube (produced by madness::plot_cubefile) and data.tree.json (produced by madness::print_tree_jsonfile) exist in /path/to/data
MRAMeshOrbitalPlot3D["/path/to/data"]
```
