# Unified Persistence Schema (FD + ES)

Status: **design, agreed 2026-05-27.** Implements the concrete schema behind
the conceptual naming in `03_naming_conventions.md` and feeds the (currently
empty) property layer in `04_property_design.md`. Builds on the ES root
identity landed in Increment 9 (commit afe80db14).

## Problem

FD response and excited-state (ES) response are solved over the same protocol
ramp but persist independently. Properties that **combine** them — resonant
Raman (FD at an ES excitation energy), two-photon absorption — are only
physically meaningful when both inputs were computed at the **same protocol
accuracy**. Today there is no shared key to assert that, and FD has no
metadata file at all. We need:

1. one protocol identity both sides agree on, and
2. one metadata file that records FD states, ES bundles, and the property
   links between them, all keyed by that protocol identity.

## Decision 1 — protocol key = physical `(thresh, k)`

The join key is the physical accuracy, **not** a positional ramp index (an
index collides across runs with different ramps and can't assert "same
accuracy"). This reuses the authority already in `ResponseProtocol.hpp` and the
`mad.fock.json` precedent in `GroundState.cpp` (`"thresh: <t> k: <k>"`).

New, in `ResponseProtocol.hpp` (single source of truth):

```cpp
/// Canonical, filename-safe protocol identity. "1e-06_k8".
inline std::string protocol_key(double thresh, int k) {
  char buf[32];
  std::snprintf(buf, sizeof buf, "%.0e_k%d", thresh, k);   // glibc: "1e-06"
  return buf;
}
/// From the active FunctionDefaults<3> (the common funnel both solvers use).
inline std::string protocol_key() {
  return protocol_key(FunctionDefaults<3>::get_thresh(),
                      FunctionDefaults<3>::get_k());
}
```

`index` (position in the ramp) is recorded too, but only as an ordering hint —
never as an identity.

## Decision 2 — one unified `response_metadata.json`

Single top-level file. A shared `protocols` registry, plus `fd_states`,
`excited_states`, and `properties` subtrees, all keyed by `protocol_key`.
Property matching is then a string compare inside one file.

```jsonc
{
  "schema_version": 1,
  "protocols": {
    "1e-04_k6": { "thresh": 1e-4, "k": 6, "index": 0 },
    "1e-06_k8": { "thresh": 1e-6, "k": 8, "index": 1 }
  },

  "fd_states": {
    "dipole_x": {                              // perturbation_description
      "1e-06_k8": {                            // protocol_key
        "f0.05700": {                          // freq_key
          "freq": 0.057,
          "type": "full", "shell": "closed_shell",
          "converged": true, "iter": 8,
          "bsh_residual": 1.2e-7,
          "archive": "dipole_x__1e-06_k8__f0.05700"
        }
      }
    }
  },

  "excited_states": {
    "1e-06_k8": {                              // protocol_key
      "type": "tda", "shell": "closed_shell", "n_roots": 2,
      "bundle_dir": "es_bundle__1e-06_k8",
      "converged": true,
      "slot_permutation": [0, 1],
      "roots": [
        { "stable_index": 0, "root_id": "es_root_0000", "slot": 0,
          "omega": 0.468, "display_name": "S1",
          "bsh_residual": 7e-5, "density_residual": 1e-3 }
      ]
    }
  },

  "properties": {
    "resonant_raman": {
      "1e-06_k8": [
        { "es_root_id": "es_root_0001", "fd_freq": 0.481, "value": 0.0 }
      ]
    }
  }
}
```

### Authority split (avoids dual-write drift)

- **`response_metadata.json` is authoritative for metadata** — protocol
  registry, FD/ES status, root identity, property links, discovery/matching.
- **Binary archives are authoritative for state data.** A FD point is a single
  `ResponseStateXY` archive; an ES bundle is a directory of per-root archives.
- The ES bundle dir keeps a **minimal `bundle.json`** (type/shell/n_roots/k/
  thresh + slot→file map + stable_index) so `load_es_roots` stays
  self-contained and collective without reading the global file. The richer
  copy in `excited_states/<key>` is the aggregate mirror, rebuilt from the
  per-bundle manifest on save. The per-artifact manifest wins on conflict.

## Archive naming (both use `protocol_key`)

```
FD point:    <pert>__<protocol_key>__f<freq>      dipole_x__1e-06_k8__f0.05700
ES bundle:   es_bundle__<protocol_key>/           es_bundle__1e-06_k8/
               ├── bundle.json                       (minimal loader manifest)
               ├── root_0  root_1 ...                (per-root archives)
```

`freq` formatted `f%.5f` (matches the existing FD `SolveKey` ordering in
`main.cpp`). `__` separates fields so a name parses unambiguously.

## Restart precedence (unchanged in spirit, now key-driven)

**FD** target `(pert, protocol_key_i, freq_j)`:
1. exact `(pert, protocol_key_i, freq_j)` archive → reload
2. lower-protocol same `(pert, freq_j)` → project up
3. nearby freq same protocol → initial guess
4. fresh perturbation-seeded guess

**ES** target bundle at `protocol_key_i`:
1. exact `es_bundle__<protocol_key_i>` → reload (matches roots by stable_index)
2. nearest lower protocol → project up, carry stable_index
3. fresh guess

Both read the `protocols` registry to order protocols by `index` for "lower
protocol" lookups.

## Property-matching contract

A property combining ES root `r` with an FD solve at `freq = ω_r` is valid
**iff `fd.protocol_key == es.protocol_key`**. Because ω is an eigenvalue that
shifts per protocol, the FD frequency for resonant Raman is taken from
`excited_states/<protocol_key>/roots[*].omega` for the matching root_id — never
hard-coded. The computed property is recorded under
`properties/<name>/<protocol_key>/` with both the `es_root_id` and the
`fd_freq` used, so the provenance (which protocol, which root) is explicit.

## Implementation increments

- **13a — protocol_key:** add the two functions to `ResponseProtocol.hpp`;
  retrofit ES `roots.json`/`bundle.json` to record `protocol_key` (it already
  records k/thresh). No behaviour change. Cheap, lands first.
- **13b — `ResponseMetadata` writer:** a small JSON aggregator (load-or-create
  `response_metadata.json`, upsert protocol/fd/es entries, atomic rank-0
  write). Both solvers register through it.
- **13c — FD save/load:** `fd_save_load.hpp` + `response_filename(...)` using
  per-state `ResponseStateXY::save/load`; FD restart precedence; `--save`/
  `--load` in `test_v3_fd_skeleton.cpp`.
- **13d — ES into the aggregate:** ES save also upserts `excited_states/<key>`
  into `response_metadata.json` (keep the per-bundle manifest as loader truth).
- **13e — property match (later, with 04):** the assert + `properties/` writer,
  once the first ES-derived property is implemented.

Each increment is independently testable; 13a/13b/13c restore the symmetry the
ES side already has, before any property work.
```
