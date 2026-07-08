# Raman State-Parallel Baseline Report (2026-02-26)

## Scope

This baseline captures the latest `raman_paper` state-parallel runs on Seawulf before adaptive scheduling changes.

- Dataset root: `/gpfs/scratch/ahurtado/project_data/data/raman_paper`
- Molecules analyzed: `H2O`, `NH3`, `CH4`
- Configurations compared:
  - `mra-p07-no-mul-sparse`
  - `mra-p07-add-mul-sparse`

## Executive Summary

1. `mra-p07-add-mul-sparse` did not run the new state-parallel implementation and is not a valid performance comparison baseline.
2. `mra-p07-no-mul-sparse` used the newer executable and state-parallel keys correctly.
3. `CH4` and `NH3` completed successfully in parallel.
4. `H2O` reproduced subgroup failure and serial fallback with `Hung queue` timeout symptoms.
5. The `H2O` failure pattern is consistent with straggler-sensitive fence waits plus the default 900 s `await_timeout`.

## Run Matrix

| Molecule | Case | Processors | Git descriptor | State-parallel keys parsed | Outcome |
|---|---|---:|---|---|---|
| CH4 | no-mul-sparse | 32 | `last_svn-6016-g45def148f-dirty` | yes | success |
| NH3 | no-mul-sparse | 32 | `last_svn-6016-g45def148f-dirty` | yes | success |
| H2O | no-mul-sparse | 32 | `last_svn-6016-g45def148f-dirty` | yes | subgroup failure + fallback |
| CH4 | add-mul-sparse | 32 | `last_svn-5978-g3d4f60902-dirty` | no | ran older code path |
| NH3 | add-mul-sparse | 32 | `last_svn-5978-g3d4f60902-dirty` | no | ran older code path |
| H2O | add-mul-sparse | 32 | `last_svn-5978-g3d4f60902-dirty` | no | ran older code path |

## Evidence (Key Log Anchors)

- Old executable + ignored keys (`add-mul-sparse`):
  - `/gpfs/scratch/ahurtado/project_data/data/raman_paper/data/H2O/mra-p07-add-mul-sparse/mad.raman.statepar32.out:35`
  - `/gpfs/scratch/ahurtado/project_data/data/raman_paper/data/H2O/mra-p07-add-mul-sparse/mad.raman.statepar32.out:41`
- New executable + state-parallel plan (`no-mul-sparse`):
  - `/gpfs/scratch/ahurtado/project_data/data/raman_paper/data/H2O/mra-p07-no-mul-sparse/mad.raman.statepar32.out:35`
  - `/gpfs/scratch/ahurtado/project_data/data/raman_paper/data/H2O/mra-p07-no-mul-sparse/mad.raman.statepar32.out:17131`
- H2O subgroup failure:
  - `/gpfs/scratch/ahurtado/project_data/data/raman_paper/data/H2O/mra-p07-no-mul-sparse/mad.raman.statepar32.out:17149`
- Hung queue / timeout:
  - `/gpfs/scratch/ahurtado/project_data/data/raman_paper/data/H2O/mra-p07-no-mul-sparse/mad.raman.statepar32.err:43`
  - `/gpfs/projects/rjh/adrian/development/madness-worktrees/molresponse-feature-next/src/madness/world/thread.h:1470`
  - `/gpfs/projects/rjh/adrian/development/madness-worktrees/molresponse-feature-next/src/madness/world/thread.cc:76`

## H2O Load-Balance Snapshot (no-mul-sparse)

- Planned channels/frequencies at protocol 0:
  - 12 channels total
  - 3 dipole channels with 6 frequencies each
  - 9 nuclear channels with 1 frequency each
- Runtime plan:
  - `requested_groups=32`
  - `channel_owner_groups=12`
  - `effective_point_groups=27`
  - policy: `ti=0 -> channel_series`, `ti>=1 -> channel_point`

Observed subgroup timing spread (from subgroup console logs):

- `thresh=1e-04`: max/min nonzero subgroup wall-sum ratio ~`6.44`
- `thresh=1e-06`: max/min nonzero subgroup wall-sum ratio ~`13.69`
- Longest single subgroup point timing:
  - `Nuc_0x @ 1e-06`: ~`909.8 s`
  - `/gpfs/scratch/ahurtado/project_data/data/raman_paper/data/H2O/mra-p07-no-mul-sparse/mad.raman.statepar32/task_1/molresponse/response_console.group18.log:168`

This leaves idle groups waiting at synchronization fences and increases timeout sensitivity.

## CH4 / NH3 Status (no-mul-sparse)

- CH4 completed derived-state stage:
  - `/gpfs/scratch/ahurtado/project_data/data/raman_paper/data/CH4/mra-p07-no-mul-sparse/mad.raman.statepar32.out:2907`
- NH3 completed derived-state stage:
  - `/gpfs/scratch/ahurtado/project_data/data/raman_paper/data/NH3/mra-p07-no-mul-sparse/mad.raman.statepar32.out:26891`

## Baseline Adaptive Strategy (for next implementation phase)

Use protocol-specific active groups with adaptive halving until estimated imbalance is acceptable.

Proposed runtime heuristic:

1. `active_groups_t = min(requested_groups, pending_channels_t, pending_points_t)`
2. Estimate per-lane cost from persisted point timings (`response_metadata`).
3. While `imbalance_ratio > target` and `active_groups_t > min_groups`, set `active_groups_t = ceil(active_groups_t / 2)`.
4. Rebuild point ownership lanes using the adjusted `active_groups_t`.
5. Keep protocol-0 fallback policy explicit and logged.

Initial practical implication for H2O-like cases:

- `32 -> 16 -> 8` groups can reduce idle-fence pressure when work cardinality is small/skewed.

## Frequency-Coupled Protocol-0 Research Direction

Current restart metadata already confirms nearest-frequency reuse in protocol-0 chain solves (`previous_frequency_memory`).

Next step candidates:

1. Shared-subspace multi-frequency solve for each channel in protocol-0.
2. Shifted-system Krylov/recycling strategy across nearby frequencies.
3. Block KAIN/nonlinear mixing with cross-frequency vectors as candidate subspace enrichments.

These should be evaluated after stabilizing adaptive lane sizing to separate scheduling effects from solver effects.

## New Systems Added For Scaling Tests

Two larger aromatic systems were added to the Raman dataset:

- `Anthracene` (C14H10, 24 atoms)
- `Pyrene` (C16H10, 26 atoms)

Added files:

- `/gpfs/scratch/ahurtado/project_data/data/raman_paper/molecules/Anthracene.mol`
- `/gpfs/scratch/ahurtado/project_data/data/raman_paper/molecules/Pyrene.mol`
- full run scaffolding under:
  - `/gpfs/scratch/ahurtado/project_data/data/raman_paper/data/Anthracene`
  - `/gpfs/scratch/ahurtado/project_data/data/raman_paper/data/Pyrene`

Initial geometry source (3D starting guess):

- Anthracene CID 8418 (PubChem PUG REST 3D SDF)
- Pyrene CID 31423 (PubChem PUG REST 3D SDF)
- URLs:
  - `https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/cid/8418/record/SDF/?record_type=3d`
  - `https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/cid/31423/record/SDF/?record_type=3d`

## Submitted Optimization Jobs (2026-02-26)

MADQC `moldft` geometry optimization jobs:

- `anthracene_moldft` -> job `1722295`
- `pyrene_moldft` -> job `1722296`

DALTON optimize stage jobs (via `seawulf-project-data-ops` wrapper, molecules subset only):

- Anthracene/aug-cc-pVDZ -> `1722297`
- Anthracene/aug-cc-pVQZ -> `1722298`
- Anthracene/aug-cc-pVTZ -> `1722299`
- Pyrene/aug-cc-pVDZ -> `1722300`
- Pyrene/aug-cc-pVQZ -> `1722301`
- Pyrene/aug-cc-pVTZ -> `1722302`

Monitor snapshot created:

- `/gpfs/scratch/ahurtado/project_data/data/raman_paper/reports/new_systems_monitor_2026-02-26.csv`
- Snapshot summary: 6 DALTON optimize jobs pending/running under Raman data scope (molecules subset `Anthracene`, `Pyrene`).
