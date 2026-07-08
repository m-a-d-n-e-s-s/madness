# Naming Conventions for Response States

> **Concrete implementation schema:** the unified protocol key + the single
> `response_metadata.json` that ties FD and ES together (so properties can
> assert a common protocol) are specified in
> `13_unified_persistence_schema.md`. This doc is the conceptual rationale;
> doc 13 is the authoritative on-disk layout.

## Current FD Response Naming

The FD naming system encodes three things in the state name:

1. **Perturbation identity** — what's perturbing the system (e.g., `Dipole_x`)
2. **Protocol** — the truncation threshold (encoded as a key like `thresh_1e-04`)
3. **Frequency** — the response frequency (encoded as a key like `freq_0.000` or `freq_0.056`)

A `LinearResponseDescriptor` holds the full description of a channel:
perturbation + list of frequencies + list of protocol thresholds.

A `LinearResponsePoint` pins a descriptor to a specific (thresh_index, freq_index),
giving one concrete solve target.

File naming via `response_filename(thresh_index, freq_index)` produces the
archive filename. The metadata record keys by:

```
states / <perturbation_description> / protocols / <protocol_key> / saved|converged / <freq_key>
```

This works because for FD response, the identity of a state is fully
determined by (perturbation, protocol, frequency). Frequency is an input —
you choose it, you solve at it, it doesn't change during iteration.

### Restart Strategy (FD)

For a target point (perturbation, protocol_i, freq_j):
1. Check if exact (perturbation, protocol_i, freq_j) archive exists → reload
2. Try lower protocol for same (perturbation, freq_j) → project up
3. Try nearby frequency at same protocol → use as initial guess
4. Fall back to fresh perturbation-projected guess

This is clean because the frequency is stable — it's the same before
and after the solve, so you can name files by it and find them later.

---

## The ES Naming Problem

For excited states, the "frequency" is the excitation energy omega, which
is an eigenvalue — an OUTPUT of the solve, not an input. It changes every
iteration. At protocol_0 you might get omega = 0.3851, and at protocol_1
after refinement you might get omega = 0.3847. Even within a single protocol,
omega shifts as the subspace rotates.

This means **you cannot use frequency as part of the state identity** the
way FD does. If you did, the file you saved at protocol_0 would have one
frequency in its name, and at protocol_1 you'd be looking for a different
frequency and wouldn't find it.

Additionally, excited states come in bundles. You solve for N roots
simultaneously. The roots can reorder — what was root 3 at protocol_0
might become root 4 at protocol_1 if energies cross. So array position
is also not a stable identity.

### What's Already Built

The current code (from the checklist and ExcitedStateBundleSolver) has
already addressed the identity problem at the metadata level:

**`ExcitedRootDescriptor`** — each root gets:
- `root_id`: stable string key, e.g. `"es_root_0001"` — never changes
- `stable_index`: permanent integer identity assigned at first appearance
- `slot_index`: current position in the response space vector (CAN change)
- `energy`: current excitation energy (CAN change)
- `display_name`: human-readable name assigned at protocol boundaries
  using energy grouping with tolerance `10 * threshold` and suffix `a/b/c`

**`slot_permutation`** — maps current array slot → stable_index, tracking
how roots have reordered across protocols.

**Root naming is done at protocol boundaries, not every iteration.** This
avoids churn during iterative convergence.

**Restart snapshots** are protocol-keyed (like FD), but the content is
the entire bundle — all roots together — not individual root archives.

---

## Proposed Unified Naming Design

### Principle: Identity must be stable across iterations and protocols

For FD: identity = (perturbation, frequency) — both are inputs, both stable.
For ES: identity = root_id (stable_index) — assigned once, never changes.

### FD State Naming (unchanged)

```
Archive:  <perturbation>_p<protocol>_f<frequency>
Metadata: states/<perturbation>/protocols/<protocol_key>/<freq_key>
```

Identity key: (perturbation_description, freq_key)
Protocol key: protocol threshold
Both are known before the solve starts.

### ES State Naming (proposed)

The identity of an excited state is its `root_id` / `stable_index`.
The frequency (excitation energy) is metadata attached to the identity,
not part of the identity itself.

**Bundle-level archive naming:**
```
Archive:  es_bundle_p<protocol>
Metadata: excited_states/protocols/<protocol_key>/
```

The bundle archive contains ALL roots for that protocol. This is
natural because the ES solver produces the full bundle at once —
roots are coupled and cannot be saved independently.

**Per-root identity within the bundle:**
```
Root identity: es_root_<stable_index>  (e.g., es_root_0001)
Slot mapping:  slot_permutation[slot] = stable_index
```

The root_id is the permanent name. The slot_index is the current
array position (which can change due to reordering).

**Why not per-root archives?**
Because excited-state roots are coupled. You can't meaningfully save
root 3 without root 4 — they share a subspace and their identities
are defined relative to each other (through the diagonalization).
The bundle IS the unit of persistence.

### ES Restart Strategy (proposed)

For a target bundle at protocol_i:
1. Check if `es_bundle_p<protocol_i>` exists → reload full bundle
2. Try lower protocol: `es_bundle_p<protocol_{i-1}>` → project up
   - Use the slot_permutation from protocol_{i-1} to match roots
   - The root_ids carry across: es_root_0001 at protocol_0 is still
     es_root_0001 at protocol_1, even if it moved from slot 2 to slot 3
3. Fall back to fresh guess (random or perturbation-seeded subspace)

The slot_permutation is critical for step 2. When you load a lower-protocol
bundle, the roots may be in a different order than what the new protocol
will converge to. The stable_index lets you track which root is which
even as slots shuffle.

---

## How ES Roots Connect to Properties

This is the downstream question: once you have converged excited states,
how do property computations (two-photon absorption, resonant Raman)
reference specific roots?

### Option A: By root_id

Property computation says: "I need es_root_0003 at protocol_2."
It looks up the bundle archive, finds the slot_permutation, and
extracts the vector at the slot corresponding to stable_index 3.

Pro: Completely stable. Works even if roots reorder.
Con: Requires the slot_permutation lookup every time.

### Option B: By energy range

Property computation says: "I need the excited state near omega = 0.385."
It looks up the bundle metadata, finds the root whose energy is closest,
and uses that root's slot.

Pro: More physical — you're asking for a state by its energy.
Con: Ambiguous if roots are near-degenerate. Fragile across protocols
     (energy shifts slightly).

### Option C: By display_name

Property computation says: "I need state S1_a."
It looks up the root whose display_name matches.

Pro: Human-readable.
Con: Display names are assigned with tolerance-based grouping and can
     change across protocols if energies shift enough.

### Recommendation: root_id as primary, energy as diagnostic

Use root_id for all programmatic references. Use energy and display_name
for human-facing output and validation. This matches what's already
built in ExcitedRootDescriptor.

---

## Naming for ES-Derived Properties

When computing properties that depend on specific excited states
(e.g., two-photon absorption cross-section for state S1), the
property record should reference the root_id:

```
properties / two_photon_absorption / es_root_0001 / protocol_key / value
```

This way the property is tied to a stable identity, not an array slot
or an approximate energy.

For resonant Raman, where you need response at a frequency matching
an excitation energy, the connection is:
1. Solve ES bundle → get converged root energies
2. Use root energy as the frequency for FD response solve
3. Record the FD state with its normal naming (perturbation + frequency)
4. Record the link: "this FD solve at freq=0.3847 corresponds to
   es_root_0001 with converged energy 0.3847"

The FD naming system doesn't need to change — it just receives a
frequency that happens to come from an excited-state calculation.
The link between the two is metadata, not naming.

---

## Summary: What Each Naming Layer Needs

| Layer          | FD Response              | ES Response                  |
|----------------|--------------------------|------------------------------|
| State identity | (perturbation, frequency)| root_id (stable_index)       |
| Archive unit   | Single vector per point  | Full bundle per protocol     |
| Archive name   | pert_p{proto}_f{freq}    | es_bundle_p{proto}           |
| Metadata key   | states/pert/proto/freq   | excited_states/proto/root_id |
| Restart match  | Same pert+freq, lower protocol | Same bundle, lower protocol, slot_permutation maps roots |
| Property ref   | (perturbation, frequency)| root_id                      |
| Frequency role | Input (stable, part of identity) | Output (eigenvalue, metadata only) |
