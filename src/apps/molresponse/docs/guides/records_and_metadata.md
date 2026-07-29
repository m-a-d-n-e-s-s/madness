# Records: metadata, results, and restart

Where a response calculation's results and state live, what the files guarantee,
and the design reasoning behind those guarantees. This is the interface to script
against — if you consume MADNESS response output programmatically, this is the
contract.

## Two files, two audiences

| file | scope | audience |
|---|---|---|
| `<prefix>.calc_info.json` | the whole workflow — one entry per task | the driver's record: what ran, in what order, and whether it succeeded |
| `response_metadata.json` (in the calc dir) | the response calculation | the property and state record: every tensor, every converged state, keyed and provenanced |
| `<prefix>.out` | the whole workflow | humans |

The rest of this page is about `response_metadata.json`, because that is the one
with a schema you build on.

## Structure

```
{
  "schema_version": 1,
  "protocols":      { "<protocol_key>": { thresh, k, index } },
  "ground_state":   { archive, fnv1a64, ... },
  "fd_states":      { "<pert>": { "<protocol_key>": { "<freq_key>":
                      { freq, type, shell, converged, iter, bsh_residual,
                        archive, backend, writer_nproc } } } },
  "excited_states": { "<protocol_key>": { type, shell, n_roots, bundle_dir,
                        converged, slot_permutation, roots[] } },
  "vbc_states":     { "<vbc_id>": { "<protocol_key>": { ... } } },
  "properties":     { "<name>": { "<protocol_key>": [ { ... }, ... ] } },
  "run_summary":    { stop_reason, dropped_work, ... },
  "io":             { backend, hdf5_compiled }
}
```

`<protocol_key>` is the resolution rung, e.g. `1e-04_k6`, `1e-06_k8` — the
threshold and the polynomial order together, because they always move together.

## The three conventions that matter

### 1 · Everything is keyed by protocol

A run that climbs the resolution ladder produces results at each rung, and they
are stored **side by side** rather than overwriting. So `properties.alpha` may hold
both a `1e-04_k6` and a `1e-06_k8` entry for the same frequency.

*Why:* the coarse rung is not garbage — it is the seed for the fine rung, and it is
also the evidence that the fine rung converged to the same place. Overwriting it
would destroy the convergence trail. It also means a consumer can ask "give me the
best available" or "show me the rung-to-rung drift" from one file.

### 2 · Rows are upserted on their identity fields

Within a `(property, protocol)`, a row's identity is the tuple of fields that
physically distinguish it:

| property | identity |
|---|---|
| `alpha` | (`omega`, `directions`) |
| `beta`, `raman` | (`A`, `B`, `C`, `freq_b`, `freq_c`) |
| `tpa`, ES-derived | (`es_root_id`, `fd_freq`) |

Re-running a calculation **replaces** the row with the same identity rather than
appending a second copy.

*Why:* an append-only record accumulated duplicate and stale rows across restarts,
so a consumer could not tell which α at ω = 0 was the current one. Identity-based
upsert makes re-running idempotent, which is what you want when the whole point of
the restart machinery is that runs get re-run.

### 3 · Rows carry provenance, not just values

Every record says how much to trust it:

- `converged` — did the state that produced this actually converge
- accuracy fields (`row_accuracy`, `vbc_accuracy`, `max_bsh_residual`) — the
  quality of the *inputs* that built this row
- `writer_nproc` — how many ranks wrote it
- `backend` — which archive format holds the state (`native` / `hdf5`)

*Why:* a bare number is not a result. A consumer — a plotting script, a regression
harness, a paper table — needs to distinguish a converged production value from a
partially-converged intermediate, without parsing the log. Recording accuracy at
the row level also means a property inherits the honest accuracy of the states that
built it, rather than claiming the accuracy of the rung it was assembled at.

## Restart

Converged response states are written per rung, and the metadata is what makes them
findable: on a re-run the engine looks up the best usable state for the target
(perturbation, frequency, protocol) and resumes from it — a coarse rung seeds the
fine one instead of starting from scratch.

Three properties of that path worth knowing:

- **Honest climb.** A rung that exhausts its iteration budget is recorded with
  `converged: false` and still used to seed the next rung. The record stays honest
  (it does not claim convergence) while the work is not thrown away.
- **Diverged states are never seeded from.** A blown-up state is recorded as such
  and excluded from source selection.
- **Rank-agnostic.** Native archives store their own layout and are redistributed
  on read, so a state written on one process count reloads on another. (Verified:
  a polarizability reloaded across process counts agrees to ~10⁻¹².)

The `backend` field is read back on load, so a state always loads with the format
that wrote it — toggling the I/O backend between runs cannot make a stale twin
shadow the state the metadata describes.

## Atomicity

`response_metadata.json` is written to a temporary file and renamed over the
target. On POSIX that rename is atomic within a filesystem, so a crash during the
write leaves either the previous file intact or the new one complete — never a
half-written record. Native state archives are written the same way.

The writer is rank-0-only and MPI-free by construction; callers guard with a rank
check and then fence, so other ranks observe the file only after it is complete.
Concurrent writers on one file are **not** supported — serialize at the
orchestration layer. (Processes that must share a calc directory read-only can do
so; the tooling that does this suppresses metadata writes entirely.)

## Ground-state fingerprint

`ground_state` records the archive path **and a hash of it**. On a later run the
hash is checked, so a calc directory built against one ground state cannot be
silently reused against a different one — a mixed-ground-state directory is caught
rather than producing quietly wrong properties.

## Consuming it

```python
import json
meta = json.load(open("response_metadata.json"))

# best available alpha at omega=0, preferring the finest rung
for key in sorted(meta["properties"]["alpha"], reverse=True):   # 1e-06_k8 before 1e-04_k6
    rows = [r for r in meta["properties"]["alpha"][key]
            if r["omega"] == 0.0 and r.get("converged")]
    if rows:
        print(key, rows[0]["alpha"]); break
```

Two habits worth adopting: **filter on `converged`** rather than assuming, and
**iterate protocol keys deliberately** rather than taking the first — the file
deliberately keeps every rung, so "which rung" is a choice the consumer makes.

See [example output](example_output.md) for real records, and
[`madqc/HDF5_IO.md`](../../../madqc/HDF5_IO.md) for the optional HDF5 backend.
