# molresponse

A unified MADNESS solver for **molecular response properties** — linear and
non-linear — built on multiresolution analysis (MRA). It is the response engine
behind `madqc --wf=response` (the default and only engine), designed to scale to
larger systems and to export research-quality, MRA-fidelity data. Because it is
basis-set-free, its properties converge in resolution toward the complete-basis
limit rather than being extrapolated in a Gaussian basis.

## What it computes

| Property | Status | Guide |
|---|---|---|
| **Polarizability α** (static + dynamic) | ✅ supported | [guides/polarizability.md](docs/guides/polarizability.md) |
| **First hyperpolarizability β** (static + SHG) | ✅ supported | [guides/hyperpolarizability.md](docs/guides/hyperpolarizability.md) |
| **Vibrational Raman** (single component) | ✅ supported | [guides/raman.md](docs/guides/raman.md) |
| **Excited states** (TDA / Full RPA) | ✅ supported (closed shell) | [guides/excited_states.md](docs/guides/excited_states.md) |
| **Two-photon absorption** | 🟡 preliminary (closed-shell singlets) | [guides/two_photon_absorption.md](docs/guides/two_photon_absorption.md) |
| **DALTON warm-start** (seed ES/SCF from DALTON) | ✅ supported | [guides/dalton_warm_start.md](docs/guides/dalton_warm_start.md) |
| Resonance Raman · full-tensor Raman · open-shell | 🔬 future | — |

**Scope:** closed-shell. Open-shell excited states and full-RPA are out of scope
for this release.

## Documentation

- **Feature guides** — one page per property, in [`docs/guides/`](docs/guides/):
  [polarizability](docs/guides/polarizability.md) ·
  [hyperpolarizability](docs/guides/hyperpolarizability.md) ·
  [raman](docs/guides/raman.md) ·
  [excited states](docs/guides/excited_states.md) ·
  [two-photon absorption](docs/guides/two_photon_absorption.md) ·
  [DALTON warm-start](docs/guides/dalton_warm_start.md).
- **[`docs/guides/example_output.md`](docs/guides/example_output.md)** — what a run
  prints and writes (real output).
- **[`docs/guides/formalism.md`](docs/guides/formalism.md)** — the shared response
  formalism the guides specialize.
- **[`docs/guides/parallelism.md`](docs/guides/parallelism.md)** — subworld /
  state-parallel design.
- **[`docs/FEATURE_STATUS.md`](docs/FEATURE_STATUS.md)** — the capability matrix.
- **[`../madqc/RESPONSE_PROPERTIES.md`](../madqc/RESPONSE_PROPERTIES.md)** — run
  recipes and the deck reference; **[`../madqc/HDF5_IO.md`](../madqc/HDF5_IO.md)** —
  optional HDF5 I/O.

## Using it

**Through madqc (recommended).** Put a `response` block in a madqc deck and run
the response workflow (v3 is the default — no engine line needed):

```text
response
    dipole              true
    dipole.directions   z
    dipole.frequencies  [0.0, 0.04]
    quadratic           true        # also compute β
end
```
```bash
madqc --wf=response h2o.in
```

Outputs land in `<prefix>.calc_info.json` (machine-readable) and `<prefix>.out`
(human-readable), with the property tensors in `response_metadata.json`. See
[`../madqc/RESPONSE_PROPERTIES.md`](../madqc/RESPONSE_PROPERTIES.md) for
per-property knobs and the result schema.

**Standalone binary** (development / low-level testing) — drives the response seam
directly on a ground-state archive:

```bash
molresponse --archive=<path_to_moldft.restartdata>
```

## Relationship to the other response codes

`molresponse_v2` (`MolresponseLib.hpp`) was this engine's predecessor and parity
reference during development; **it has since been removed** — this engine is now
the sole response engine and the default for `madqc --wf=response`, and a deck that
still carries `engine v2` errors out. Migrating a v2 deck is a one-line edit:
delete the `engine v2` line. See [`MIGRATION_FROM_V2.md`](MIGRATION_FROM_V2.md).
