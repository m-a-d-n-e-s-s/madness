<!-- Feature-guide template. Every guide in this directory follows this shape so
     the set reads consistently. One page per property; highlight the feature and
     its validation, point to the recipes and theory rather than duplicating them. -->

# <Property> <symbol>

One sentence: what physical quantity this is and why you'd compute it.

## What it computes

The defining expression (the tensor/quantity), in one or two lines, with the
physical meaning of the indices. Keep it to the essentials — the full derivation
lives in [the formalism guide](formalism.md).

## Status

| configuration | status |
|---|---|
| <e.g. static> | ✅ supported |
| <e.g. dynamic ω> | ✅ / 🚧 in progress |

State the supported configurations only.

## Validation / Preliminary results

Two flavors, depending on the property:

- **Established properties (α, β):** already benchmarked in the literature — cite
  the author's papers and state results are not reproduced here. No tables.
- **New properties (2PA, excited states):** include a small **preliminary**
  results table in the style of the project website (e.g. excitation energies, or
  δ_TPA/σ per state, marked clearly as preliminary), and this is where the DALTON
  warm-start (seeding) discussion lives. The exhaustive/authoritative validation
  campaign stays local for a dedicated results PR.

No graphs in the release docs.

## Run it

The minimal `madqc --wf=response` deck knobs (2–4 lines), then a pointer to the
full recipe in [`../../../madqc/RESPONSE_PROPERTIES.md`](../../../madqc/RESPONSE_PROPERTIES.md).
Say where the result lands (`response_metadata.json` → `properties/<name>`).

## Under the hood

One paragraph: the method (linear response / VBC quadratic source / single
residue / …), what it reuses, and any convention that matters. Link to
[formalism](formalism.md) for the shared (A,B,C) response machinery.

## References

- Method / benchmark papers.
- Companion guides.
