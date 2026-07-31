# Example output

What a response calculation actually prints and writes. Every excerpt below is
**real output from a run of this code**, lightly trimmed for length — not
illustrative mock-ups.

Two output surfaces:

- **`<prefix>.out`** — the human-readable summary, one section per workflow task.
- **`response_metadata.json`** (in the calc directory) — the machine-readable
  record: every property tensor keyed by protocol and by its identity fields,
  plus provenance. This is what you script against.

## A polarizability run, end to end

Deck (water, HF, static α):

```text
dft
    xc        hf
    protocol  [1e-4]
    l         200.0
    prefix    "mad"
end
response
    dipole             true
    dipole.directions  xyz
    dipole.frequencies [0.0]
end
molecule
    units atomic
    no_orient true
    O    0.00000000   -0.22490589    0.00000000
    H    1.45235000    0.89962300    0.00000000
    H   -1.45235000    0.89962300    0.00000000
end
```

`madqc --wf=response h2o.in` → `h2o.out`:

```text
======================================================================
  MADQC RESULTS SUMMARY
======================================================================

  Task 0 : SCF  (model = scf)
  ----------------------------------------------------------------------
    Molecule         : H2O
    Total energy     :      -76.065377246 Ha       -2069.8444 eV
    Orbital eps (a)  : -20.5621 -0.9593 -0.9593 -0.7225 -0.5092  (Ha)
    Converged        :  thresh = 0.0001  dconv = 0.0001
    Wall time        :  0.0 s   (1 MPI x 1 threads)

  Task 1 : RESPONSE
  ----------------------------------------------------------------------
    alpha[1e-04_k6](w=0.000)
      x :      9.482637     -0.000000     -0.000000
      y :     -0.000000      8.714917      0.000000
      z :      0.000000      0.000000      7.973547

======================================================================
```

Note the shape: the workflow ran the ground state as task 0 and the response as
task 1, and the α label carries the **protocol key** (`1e-04_k6`) and the
frequency, because both are part of the result's identity. The off-diagonal
elements are zero to ~10⁻⁷ here, which is the symmetry of this geometry showing
up as a numerical check for free.

## Excited states

Per converged root, the solver reports the excitation energy and its transition
properties (real output, water, root 0):

```text
   Response Function 0		8.54365210 eV
   --------------------------------------------

   Transition Dipole Moments
   X: 0.00000000   Y: 0.00000000   Z: -0.45734020

   Dipole Oscillator Strength: 0.04378044

   Transition Quadrupole Moments
                  X                Y                Z
   X       0.00000000      -0.00000000      -0.00000187
   Y      -0.00000000       0.00000000       0.24705616
   Z      -0.00000187       0.24705616      -0.00000000

   Dominant Contributions:
                  x          y
   Occupied 4   0.99991358 0.02164660
   Occupied 3   0.02094647 0.02793969
   Occupied 1   0.01199088 0.00997191
```

**Dominant Contributions** is the useful diagnostic for identifying *which*
state you got: it lists the occupied orbitals carrying the excitation, with
separate `x` (excitation) and `y` (de-excitation) weights. Here the state is
essentially a single excitation out of occupied orbital 4 (weight 0.9999) with
negligible de-excitation character — as expected for a low valence excitation.
Use this, not the root index, when matching states across calculations.

## Two-photon absorption

For a 2PA run the engine prints the transition tensor per state and then the
rotationally-invariant observables (real output, LiH — note roots 2 and 3 are a
degenerate pair):

```text
                  | Two-photon transition tensor S |
                  +--------------------------------+
     -----------------------------------------------------------------
      No  Energy(eV)     Sxx     Syy     Szz     Sxy     Sxz     Syz
     -----------------------------------------------------------------
       1       4.06    25.792  25.792  75.255   0.000  -0.000   0.000
       2       5.07    -0.000   0.000   0.000  -0.000  -7.372 -42.533
       3       5.07    -0.000  -0.000  -0.000   0.000 -42.533   7.372
     -----------------------------------------------------------------

     D = 2*Df+4*Dg (Linear);  D = -2*Df+6*Dg (Circular)
     Df = sum_ij S_ii*S_jj /30;   Dg = sum_ij S_ij^2 /30
     sigma = D*(E/2)^2*AU_TO_GM (GM, 0.1 eV FWHM);  R = (-Df+3Dg)/(Df+2Dg)
```

The definitions are printed with the numbers deliberately, so a table pasted into
a notebook is self-documenting. **Read the invariants, not the components,** for
degenerate states: roots 2 and 3 above distribute the same physics differently
between `Sxz` and `Syz` (the manifold is free to rotate), while their invariants
are identical — which is why any cross-code or cross-basis comparison must be
made on `Df`/`Dg`/`D`/`R`.

## The machine-readable record

`response_metadata.json` in the calc directory holds the same results in a stable
shape. A polarizability row:

```json
"properties": {
  "alpha": {
    "1e-06_k8": [
      { "omega": 0.0, "directions": "z",
        "alpha": [[8.541974761089701]],
        "converged": true, "row_accuracy": ..., "max_bsh_residual": ... }
    ]
  }
}
```

and a two-photon row:

```json
"properties": {
  "tpa": {
    "1e-06_k8": [
      { "es_root_id": 0, "omega": 0.14906966931818236, "omega_ev": 4.056392349378767,
        "S": [[25.79195, 2.2e-12, -2.98e-05],
              [2.2e-12, 25.79195, 2.14e-05],
              [-2.98e-05, 2.14e-05, 75.25462]],
        "Df": 536.267, "Dg": 233.124,
        "D_linear": 2005.028, "D_circular": 326.207, "R": 0.1627,
        "sigma_linear_gm": 24.173, "sigma_circular_gm": 3.933,
        "writer_nproc": 1 }
    ]
  }
}
```

Three conventions worth knowing:

- Results are keyed by **protocol** (`1e-06_k8`), so a run that climbs the ladder
  keeps every rung side by side rather than overwriting.
- Within a protocol, rows are **upserted on their identity fields** (for α that is
  (ω, directions); for 2PA it is `es_root_id`) — so re-running a calculation
  replaces its row instead of appending a duplicate.
- Rows carry provenance (`converged`, accuracy fields, `writer_nproc`), so a
  consumer can tell how trustworthy a number is without re-reading the log.

## Where files land

```
<calc-dir>/
  response_metadata.json          the record above
  es__<protocol>/                 excited-state bundles (one archive per root)
  <pert>__<protocol>__<freq>      converged response states (restartable)
  es_analysis__<protocol>.json    transition properties
<prefix>.out                      the human summary
<prefix>.calc_info.json           the workflow-level record (all tasks)
```

Response states are written per protocol rung and are **restartable**: a re-run
reloads the coarse rung and climbs from it rather than starting over. See
[`madqc/RESPONSE_PROPERTIES.md`](../../../madqc/RESPONSE_PROPERTIES.md) for the
full deck reference and [`madqc/HDF5_IO.md`](../../../madqc/HDF5_IO.md) if you
want these archives as HDF5 instead of the native format.
