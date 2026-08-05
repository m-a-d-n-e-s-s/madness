# MADQC Workflow Builder Guide

> This is the **developer** guide for extending madqc. For the **user** guide
> (running calculations, input format, output/schema, visualization), see
> [`README.md`](README.md).

`madqc` builds execution pipelines through `src/madness/chem/WorkflowBuilders.hpp`.

For personal Seawulf interactive build/test commands, see `src/apps/madqc/SEAWULF_INTERACTIVE_WORKFLOW.md`.
For molresponse state-level subgroup-parallel design notes, see `src/apps/molresponse_v2/STATE_PARALLEL_DESIGN.md`.

## Current Architecture

1. `madqc.cpp` parses CLI and creates `Params` + `qcapp::Workflow`.
2. `madqc.cpp` calls:

```cpp
workflow_builders::add_workflow_drivers(world, pm, user_workflow, wf);
```

3. `WorkflowBuilders.hpp` maps workflow names (`scf`, `nemo`, `response`, `mp2/cc2`, `cis`, `oep`) to driver wiring, and `--optimize` to the optimizer variant of a workflow's reference.
4. `wf.run(prefix)` executes the assembled pipeline.

Workflow names are also centralized in `WorkflowBuilders.hpp`:

- `workflow_kind_from_name(...)`
- `runnable_workflows`
- `runnable_workflow_list()`

## How To Add A New Workflow

1. Define the workflow name and expected behavior (single step or multi-step pipeline).
2. Implement a helper in `WorkflowBuilders.hpp`, for example:

```cpp
inline void add_myworkflow_drivers(World& world, Params& pm, qcapp::Workflow& wf) {
  // parameter tweaks (optional)
  // create application(s)
  // add one or more SinglePointDriver entries
}
```

3. Register it in `add_workflow_drivers(...)`:

```cpp
} else if (user_workflow == "myworkflow") {
  add_myworkflow_drivers(world, pm, wf);
}
```

4. Update `help()` text and available workflow list in `src/apps/madqc/madqc.cpp`.
   This should use `runnable_workflow_list()` instead of hard-coded strings.
5. If there is a standalone compatibility executable, keep it thin and call the same builder helper (no duplicated orchestration logic).
6. Add at least one integration test under `src/apps/madqc/` using the new `--wf=<name>` route.
7. Update `src/apps/madqc/test_workflow_builders.cpp` with the new workflow name and expected kind mapping.

## A flag that changes what a workflow's reference does

`--optimize` is the pattern for "same reference engine, different thing done to
it": `add_workflow_drivers` takes a bool and, when set, routes to
`add_optimize_workflow_drivers`, which switches on the *workflow kind* to pick
`OptimizeDriver<moldft_lib>` or `OptimizeDriver<nemo_lib>`. The engine is named
once, by `--wf`. Prefer this over inventing a workflow name per combination, and
prefer a Driver over an Application whenever the step may own sub-runs — the
optimizer will, once numerical gradients arrive as displaced sub-calculations.
Engine-specific code stays behind `if constexpr (std::is_same_v<Calc, SCF>)`
inside the driver: for the optimizer that is only the target
(`MolecularEnergy` wraps an SCF, `Nemo` is already an
`OptimizationTargetInterface`) and the SCF-only protocol preparation.

## Workflow Design Rules

1. Keep orchestration in `WorkflowBuilders.hpp`; keep algorithm code in libraries/applications.
2. Prefer composable driver chains over one large monolithic app.
3. Apply parameter-derived logic in one place before driver creation.
4. Throw explicit errors for unsupported/disabled workflows.

## Minimal Pattern

```cpp
// In builder helper:
auto reference = std::make_shared<SCFApplication<moldft_lib>>(world, pm);
wf.addDriver(std::make_unique<qcapp::SinglePointDriver>(reference));
wf.addDriver(std::make_unique<qcapp::SinglePointDriver>(
    std::make_unique<ResponseApplication<molresponse_lib>>(world, pm, reference->calc())));
```

This pattern keeps workflow composition centralized and reusable across executables.
