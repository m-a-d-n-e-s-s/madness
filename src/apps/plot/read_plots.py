#!/usr/bin/env python3
"""Read the plot files moldft writes (.cube or .dx) into numpy and analyse them.

Needs numpy; scipy, scikit-image and matplotlib only for --diff/--surface/--slice.

  read_plots.py total_density.cube                      # header, grid, min/max, grid integral
  read_plots.py --band total_density.cube esp.cube      # ESP statistics on the 8e-5 <= rho <= 1e-4 isodensity band
  read_plots.py --band --surface esp.png total_density.cube esp.cube   # rho = 1e-4 isosurface coloured by the ESP
  read_plots.py --diff ref_rho.cube ref_esp.cube other_rho.cube other_esp.cube [--surface diff.png]
        # ESP(other) - ESP(ref) and rho(other) - rho(ref) on the reference band
  read_plots.py --orbitals dirA dirB [--nocc 22]        # grid overlaps |<a_i|b_j>| between the amo-*.cube of two runs
  read_plots.py --extent dir [--radius 8]               # per orbital: |phi|^2 within R of the nuclear centroid, rms radius, edge fraction
  read_plots.py --slice amo-00022.cube --axis 0 --value 0.0 --png slice.png   # plane cut of a plotted function

Conventions: atomic units on the grid of the calculation frame (moldft recentres, and by default rotates,
the molecule; atom positions are in the cube header). coulomb = vnuc + vcoul (electron convention);
esp = -coulomb = potential felt by a positive test charge, Hartree per e (x 27.2114 eV, x 627.51 kcal/mol).
.dx: ASCII header + little-endian float64, C order (z fastest). .cube: Gaussian cube, z fastest, 5 digits.
For interactive viewing use VMD (.dx and .cube), ParaView or Avogadro (.cube).
"""
import sys, re, argparse
import numpy as np

def read_cube(fn):
    L = open(fn).read().split('\n')
    nat = int(L[2].split()[0]); origin = np.array(L[2].split()[1:4], float)
    n = [int(L[3+i].split()[0]) for i in range(3)]
    delta = np.array([float(L[3+i].split()[1+i]) for i in range(3)])
    atoms = [(int(l.split()[0]), np.array(l.split()[2:5], float)) for l in L[6:6+nat]]
    vals = np.array(' '.join(L[6+nat:]).split(), float).reshape(n)
    return dict(origin=origin, delta=delta, n=np.array(n), atoms=atoms, values=vals, comment=L[1])

def read_dx(fn):
    b = open(fn, 'rb').read()
    i = b.index(b'binary data follows\n') + len(b'binary data follows\n')
    hdr = b[:i].decode()
    n = [int(x) for x in re.search(r'counts\s+(\d+)\s+(\d+)\s+(\d+)', hdr).groups()]
    origin = np.array(re.search(r'origin\s+(\S+)\s+(\S+)\s+(\S+)', hdr).groups(), float)
    deltas = [[float(x) for x in l.split()[1:4]] for l in hdr.splitlines() if l.startswith('delta')]
    delta = np.array([deltas[k][k] for k in range(3)])
    cnt = int(re.search(r'items\s+(\d+)', hdr).group(1))
    vals = np.frombuffer(b[i:i+8*cnt], dtype='<f8').reshape(n)
    return dict(origin=origin, delta=delta, n=np.array(n), atoms=[], values=vals, comment='')

def read(fn):
    return read_cube(fn) if fn.endswith('.cube') else read_dx(fn)

def axes(g):
    return [g['origin'][d] + g['delta'][d]*np.arange(g['n'][d]) for d in range(3)]

def summary(fn, g):
    v = g['values']; dV = np.prod(g['delta'])
    print(f"{fn}: grid {g['n']} origin {g['origin']} spacing {g['delta']} (bohr)")
    if g['comment']: print("  header:", g['comment'])
    for z, p in g['atoms']: print(f"  atom Z={z} at {p}")
    print(f"  min {v.min():.6e} max {v.max():.6e}  grid integral {v.sum()*dV:.6f}  (a density integral is only as good as the grid resolves the cores)")

def band_stats(rho, esp, lo, hi):
    assert rho['values'].shape == esp['values'].shape and np.allclose(rho['origin'], esp['origin'])
    m = (rho['values'] >= lo) & (rho['values'] <= hi)
    e = esp['values'][m]
    print(f"isodensity band {lo:g} <= rho <= {hi:g}: {m.sum()} grid points (spacing {rho['delta'][0]:.3f} bohr)")
    if e.size:
        q = np.percentile(e, [0, 25, 50, 75, 100])
        print("  ESP (Hartree/e): min %.5f  q25 %.5f  median %.5f  q75 %.5f  max %.5f" % tuple(q))
        print("  ESP (kcal/mol/e): min %.2f  q25 %.2f  median %.2f  q75 %.2f  max %.2f" % tuple(q*627.5095))
    return m, e

def surface_png(rho, esp, iso, out, vrange=None, title=None):
    """Marching-cubes isosurface of rho at `iso`, coloured by esp interpolated at the vertices; matplotlib PNG."""
    from skimage import measure
    from scipy.ndimage import map_coordinates
    import matplotlib; matplotlib.use('Agg'); import matplotlib.pyplot as plt
    from mpl_toolkits.mplot3d.art3d import Poly3DCollection
    verts, faces, _, _ = measure.marching_cubes(rho['values'], level=iso)
    e = map_coordinates(esp['values'], verts.T, order=1)
    xyz = rho['origin'] + verts*rho['delta']
    ef = e[faces].mean(axis=1)
    vmax = vrange if vrange else np.percentile(np.abs(ef), 98)
    cmap = plt.get_cmap('coolwarm_r')             # red = negative ESP, blue = positive (Jensen's convention)
    colors = cmap(0.5 + 0.5*np.clip(ef/vmax, -1, 1))
    fig = plt.figure(figsize=(9, 7)); ax = fig.add_subplot(111, projection='3d')
    ax.add_collection3d(Poly3DCollection(xyz[faces], facecolors=colors, edgecolor='none'))
    lo, hi = xyz.min(0), xyz.max(0); c = (lo+hi)/2; r = (hi-lo).max()/2
    ax.set_xlim(c[0]-r, c[0]+r); ax.set_ylim(c[1]-r, c[1]+r); ax.set_zlim(c[2]-r, c[2]+r)
    for z, pos in rho['atoms']: ax.scatter(*pos, s=12, c='k')
    ax.set_xlabel('x / bohr'); ax.set_ylabel('y / bohr'); ax.set_zlabel('z / bohr')
    sm = plt.cm.ScalarMappable(cmap=cmap, norm=plt.Normalize(-vmax, vmax)); sm.set_array([])
    fig.colorbar(sm, ax=ax, shrink=0.6, label='ESP / Hartree per e')
    ax.set_title(title or f'rho = {iso:g} isosurface, {len(faces)} faces')
    fig.savefig(out, dpi=150, bbox_inches='tight'); print(f"wrote {out}: {len(verts)} vertices, ESP on surface min {e.min():.4f} max {e.max():.4f} Hartree/e")

def diff_stats(rref, eref, roth, eoth, lo, hi, surface=None, vrange=None):
    assert rref['values'].shape == roth['values'].shape and np.allclose(rref['origin'], roth['origin'])
    m = (rref['values'] >= lo) & (rref['values'] <= hi)
    de = (eoth['values'] - eref['values'])[m]; dr = (roth['values'] - rref['values'])[m]
    print(f"reference band {lo:g} <= rho_ref <= {hi:g}: {m.sum()} points")
    q = np.percentile(de, [0, 25, 50, 75, 100])
    print("  dESP = other - ref (Hartree/e): min %.5f q25 %.5f median %.5f q75 %.5f max %.5f | MAD %.5f MaxAD %.5f" % (*q, np.abs(de).mean(), np.abs(de).max()))
    print("  drho = other - ref on the band: mean %.2e, max|drho| %.2e (band density 1e-4)" % (dr.mean(), np.abs(dr).max()))
    r = np.corrcoef(eref['values'][m], de)[0, 1]
    print(f"  correlation of dESP with ESP_ref on the band: {r:+.3f} (Jensen: negative slope = charge-transfer signature)")
    if surface:
        d = dict(eoth); d['values'] = eoth['values'] - eref['values']
        surface_png(rref, d, hi, surface, vrange, title=f'ESP(other) - ESP(ref) on the rho_ref = {hi:g} surface')

def orbital_overlaps(dirA, dirB, nocc=None):
    """Pair the orbitals of two runs by grid overlap |<a_i|b_j>| (orbitals grid-normalised); report the diagonal and the swaps."""
    import glob, os
    fa = sorted(glob.glob(os.path.join(dirA, 'amo-*.cube'))); fb = sorted(glob.glob(os.path.join(dirB, 'amo-*.cube')))
    n = min(len(fa), len(fb), nocc or 10**6)
    A = [read_cube(f)['values'].ravel() for f in fa[:n]]; B = [read_cube(f)['values'].ravel() for f in fb[:n]]
    A = np.array([a/np.linalg.norm(a) for a in A]); B = np.array([b/np.linalg.norm(b) for b in B])
    S = np.abs(A @ B.T)
    print(f"orbital overlaps {dirA} vs {dirB}: {n} orbitals; grid-normalised, sign-free")
    for i in range(n):
        j = int(np.argmax(S[i])); tag = "" if j == i else f"   <-- best match is orbital {j} of B"
        print(f"  a{i:3d}: |<a|b_same>| = {S[i,i]:.4f}   best |<a|b_j>| = {S[i,j]:.4f}{tag}")
    return S

def orbital_extent(d, R=8.0):
    """How delocalised is each plotted orbital: |phi|^2 fraction inside a sphere of radius R (bohr) around the nuclear
    centroid, rms radius, and the fraction in the outermost 10% shell of the plot box (a box-continuum state lives there)."""
    import glob, os
    files = sorted(glob.glob(os.path.join(d, 'amo-*.cube')))
    g0 = read_cube(files[0]); ax = axes(g0); X, Y, Z = np.meshgrid(*ax, indexing='ij')
    c = np.mean([p for _, p in g0['atoms']], axis=0); r = np.sqrt((X-c[0])**2 + (Y-c[1])**2 + (Z-c[2])**2)
    lo = g0['origin']; hi = g0['origin'] + g0['delta']*(g0['n']-1); L = (hi-lo)
    edge = (X < lo[0]+0.1*L[0]) | (X > hi[0]-0.1*L[0]) | (Y < lo[1]+0.1*L[1]) | (Y > hi[1]-0.1*L[1]) | (Z < lo[2]+0.1*L[2]) | (Z > hi[2]-0.1*L[2])
    print(f"orbital extent in {d}: box {np.round(lo,1)} .. {np.round(hi,1)} bohr, grid {g0['n']}, centroid {np.round(c,2)}, sphere R = {R} bohr")
    print("  orbital   |phi|^2 in R   rms radius (bohr)   edge-shell fraction   max|phi|")
    for f in files:
        v = read_cube(f)['values']; w = v*v; W = w.sum()
        print(f"  {os.path.basename(f)[4:9]}   {w[r<R].sum()/W:10.4f}   {np.sqrt((w*r*r).sum()/W):14.2f}   {w[edge].sum()/W:16.4f}   {np.abs(v).max():.3e}")

def slice_png(g, axis, value, out, log=False):
    """Colour map of the function on the plane axis = value (nearest grid plane), atoms projected."""
    import matplotlib; matplotlib.use('Agg'); import matplotlib.pyplot as plt
    ax_ = axes(g); i = int(round((value - g['origin'][axis])/g['delta'][axis])); i = min(max(i, 0), g['n'][axis]-1)
    v = np.take(g['values'], i, axis=axis); others = [d for d in range(3) if d != axis]
    A, B = ax_[others[0]], ax_[others[1]]
    fig, axp = plt.subplots(figsize=(7, 6))
    if log:
        im = axp.pcolormesh(A, B, np.log10(np.abs(v.T) + 1e-12), shading='auto', cmap='viridis'); label = 'log10 |f|'
    else:
        m = np.abs(v).max(); im = axp.pcolormesh(A, B, v.T, shading='auto', cmap='RdBu_r', vmin=-m, vmax=m); label = 'f'
    for z, p in g['atoms']: axp.plot(p[others[0]], p[others[1]], 'ko', ms=3)
    axp.set_xlabel(f"{'xyz'[others[0]]} / bohr"); axp.set_ylabel(f"{'xyz'[others[1]]} / bohr"); axp.set_aspect('equal')
    axp.set_title(f"{out}: plane {'xyz'[axis]} = {ax_[axis][i]:.2f} bohr"); fig.colorbar(im, ax=axp, label=label)
    fig.savefig(out, dpi=130, bbox_inches='tight'); print(f"wrote {out}; slice max|f| {np.abs(v).max():.3e}")

if __name__ == '__main__':
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('files', nargs='+'); ap.add_argument('--band', action='store_true', help='files = density, esp: statistics on the band')
    ap.add_argument('--lo', type=float, default=8e-5); ap.add_argument('--hi', type=float, default=1e-4)
    ap.add_argument('--diff', action='store_true', help='files = ref_rho ref_esp other_rho other_esp')
    ap.add_argument('--orbitals', action='store_true', help='files = dirA dirB: overlap matrix of the amo-*.cube orbitals')
    ap.add_argument('--nocc', type=int, default=None)
    ap.add_argument('--slice', action='store_true', help='files = one cube/dx: plane cut, see --axis/--value/--png/--log')
    ap.add_argument('--axis', type=int, default=0); ap.add_argument('--value', type=float, default=0.0)
    ap.add_argument('--png', default='slice.png'); ap.add_argument('--log', action='store_true')
    ap.add_argument('--extent', action='store_true', help='files = dir: delocalisation measures of the amo-*.cube orbitals')
    ap.add_argument('--radius', type=float, default=8.0)
    ap.add_argument('--surface', metavar='PNG', help='with --band: write the rho=hi isosurface coloured by ESP')
    ap.add_argument('--vrange', type=float, default=None, help='colour scale limit in Hartree/e (default: 98th percentile)')
    a = ap.parse_args()
    if a.orbitals:
        orbital_overlaps(a.files[0], a.files[1], a.nocc); sys.exit(0)
    if a.extent:
        orbital_extent(a.files[0], a.radius); sys.exit(0)
    if a.slice:
        slice_png(read(a.files[0]), a.axis, a.value, a.png, a.log); sys.exit(0)
    grids = {f: read(f) for f in a.files}
    for f, g in grids.items(): summary(f, g)
    if a.band:
        band_stats(grids[a.files[0]], grids[a.files[1]], a.lo, a.hi)
        if a.surface: surface_png(grids[a.files[0]], grids[a.files[1]], a.hi, a.surface, a.vrange)
    if a.diff:
        diff_stats(*[grids[f] for f in a.files[:4]], a.lo, a.hi, a.surface, a.vrange)
