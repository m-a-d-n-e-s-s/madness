#!/usr/bin/env python3

import json
import os
import shlex
import subprocess
import sys

sys.path.append("@CMAKE_SOURCE_DIR@/bin")
from test_utilities import cleanup


def run_cmd(cmd, env):
    launcher = env.get("MADQC_LAUNCHER", "").strip()
    full_cmd = shlex.split(launcher) + cmd if launcher else cmd

    print("executing\n ", " ".join(full_cmd))
    p = subprocess.run(
        full_cmd,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
        env=env,
    )
    stdout = p.stdout.decode("utf-8", errors="replace")
    stderr = p.stderr.decode("utf-8", errors="replace")
    print("finished run")
    print(stdout)
    if stderr:
        print(stderr)
    print("exitcode", p.returncode)
    return p


def find_response_task(calc_info):
    for task in calc_info.get("tasks", []):
        if task.get("type") == "response":
            return task
    return None


if __name__ == "__main__":
    print("Testing @BINARY@/@TESTCASE@")
    print("reference files found in directory: @SRCDIR@")

    prefix = "mad_@BINARY@_@TESTCASE@"
    outputfile = prefix + ".calc_info.json"
    protocol_key = f"{1.0e-2:.0e}"
    response_metadata_file = os.path.join(
        prefix, "task_1", "molresponse", "response_metadata.json"
    )

    cleanup(prefix)
    try:
        os.remove(outputfile)
    except FileNotFoundError:
        pass

    env = os.environ.copy()
    cmd = [
        "./@BINARY@",
        "--molecule=he",
        "--wf=response",
        f"--prefix={prefix}",
        "--dft=k=6;maxiter=1;econv=1e-3;dconv=1e-2;protocol=[1e-2];localize=canon",
        "--response=dipole=true;dipole.directions=z;dipole.frequencies=[0.0];"
        "excited.enable=true;excited.num_states=2;excited.guess_max_iter=1;"
        "excited.maxiter=3;excited.maxsub=4;print_level=2",
    ]
    run = run_cmd(cmd, env=env)
    if run.returncode != 0:
        sys.exit(run.returncode)

    if not os.path.exists(outputfile):
        print("Output file not found:", outputfile)
        sys.exit(1)

    with open(outputfile, "r", encoding="utf-8") as f:
        got = json.load(f)
    if not os.path.exists(response_metadata_file):
        print("response_metadata.json not found:", response_metadata_file)
        sys.exit(1)
    with open(response_metadata_file, "r", encoding="utf-8") as f:
        response_metadata = json.load(f)

    success = True

    response_task = find_response_task(got)
    if response_task is None:
        print("Did not find response task in calc_info.json")
        sys.exit(1)

    metadata = response_task.get("metadata", {})
    excited_states = metadata.get("excited_states", {})
    excited_plan = excited_states.get("plan", {})
    excited_protocols = excited_states.get("protocols", {})
    excited_exec = metadata.get("excited_state_planner", {}).get("execution", {})

    if not excited_plan.get("enabled", False):
        print("Expected excited plan to be enabled")
        success = False
    if int(excited_plan.get("num_states", -1)) != 2:
        print("Expected excited num_states=2, got", excited_plan.get("num_states"))
        success = False
    if int(excited_plan.get("guess_max_iter", -1)) != 1:
        print(
            "Expected excited guess_max_iter=1, got",
            excited_plan.get("guess_max_iter"),
        )
        success = False
    if int(excited_plan.get("maxiter", -1)) != 3:
        print("Expected excited maxiter=3, got", excited_plan.get("maxiter"))
        success = False
    if int(excited_plan.get("maxsub", -1)) != 4:
        print("Expected excited maxsub=4, got", excited_plan.get("maxsub"))
        success = False

    if protocol_key not in excited_protocols:
        print("Expected excited protocol key", protocol_key, "in metadata")
        success = False
        protocol_node = {}
    else:
        protocol_node = excited_protocols[protocol_key]

    if not isinstance(protocol_node.get("saved"), bool):
        print("Expected boolean 'saved' in excited protocol node")
        success = False
    if not isinstance(protocol_node.get("converged"), bool):
        print("Expected boolean 'converged' in excited protocol node")
        success = False
    if not isinstance(protocol_node.get("energies"), list):
        print("Expected list 'energies' in excited protocol node")
        success = False
    if not isinstance(protocol_node.get("state_names"), list):
        print("Expected list 'state_names' in excited protocol node")
        success = False
    if not isinstance(protocol_node.get("roots"), list):
        print("Expected list 'roots' in excited protocol node")
        success = False
    if not isinstance(protocol_node.get("slot_permutation"), list):
        print("Expected list 'slot_permutation' in excited protocol node")
        success = False
    if not isinstance(protocol_node.get("residual_norms"), list):
        print("Expected list 'residual_norms' in excited protocol node")
        success = False
    if not isinstance(protocol_node.get("density_change_norms"), list):
        print("Expected list 'density_change_norms' in excited protocol node")
        success = False
    if not isinstance(protocol_node.get("relative_residual_norms"), list):
        print("Expected list 'relative_residual_norms' in excited protocol node")
        success = False
    if not isinstance(protocol_node.get("iteration_max_residuals"), list):
        print("Expected list 'iteration_max_residuals' in excited protocol node")
        success = False
    if not isinstance(protocol_node.get("iteration_max_density_changes"), list):
        print(
            "Expected list 'iteration_max_density_changes' in excited protocol node"
        )
        success = False
    if not isinstance(protocol_node.get("iteration_max_relative_residuals"), list):
        print(
            "Expected list 'iteration_max_relative_residuals' in excited protocol node"
        )
        success = False
    if not isinstance(protocol_node.get("iterations"), int):
        print("Expected integer 'iterations' in excited protocol node")
        success = False
    if not isinstance(protocol_node.get("convergence_mode"), str):
        print("Expected string 'convergence_mode' in excited protocol node")
        success = False
    if not isinstance(protocol_node.get("accelerator_mode"), str):
        print("Expected string 'accelerator_mode' in excited protocol node")
        success = False
    if not isinstance(protocol_node.get("accelerator_subspace"), int):
        print("Expected integer 'accelerator_subspace' in excited protocol node")
        success = False
    if not isinstance(protocol_node.get("density_convergence_target"), (int, float)):
        print(
            "Expected numeric 'density_convergence_target' in excited protocol node"
        )
        success = False
    if not isinstance(protocol_node.get("relative_convergence_target"), (int, float)):
        print(
            "Expected numeric 'relative_convergence_target' in excited protocol node"
        )
        success = False
    if not isinstance(protocol_node.get("max_rotation"), (int, float)):
        print("Expected numeric 'max_rotation' in excited protocol node")
        success = False

    stage_status = protocol_node.get("stage_status", "")
    if not isinstance(stage_status, str) or not stage_status:
        print("Expected non-empty protocol stage_status")
        success = False
    if protocol_node.get("response_variant") != "dynamic_restricted":
        print(
            "Expected response_variant=dynamic_restricted, got",
            protocol_node.get("response_variant"),
        )
        success = False
    if protocol_node.get("restart_support_mode") != "full_bundle_resume":
        print(
            "Expected restart_support_mode=full_bundle_resume, got",
            protocol_node.get("restart_support_mode"),
        )
        success = False
    if protocol_node.get("restart_source") != "fresh_guess":
        print(
            "Expected restart_source=fresh_guess on initial run, got",
            protocol_node.get("restart_source"),
        )
        success = False
    if protocol_node.get("snapshot_kind") != "protocol_bundle":
        print(
            "Expected snapshot_kind=protocol_bundle, got",
            protocol_node.get("snapshot_kind"),
        )
        success = False
    if protocol_node.get("bundle_state_present") is not True:
        print(
            "Expected bundle_state_present=true, got",
            protocol_node.get("bundle_state_present"),
        )
        success = False
    if protocol_node.get("restart_capable") is not True:
        print(
            "Expected restart_capable=true for restricted dynamic run, got",
            protocol_node.get("restart_capable"),
        )
        success = False
    if protocol_node.get("convergence_mode") != "density_relative_dual_gate":
        print(
            "Expected convergence_mode=density_relative_dual_gate, got",
            protocol_node.get("convergence_mode"),
        )
        success = False
    if protocol_node.get("accelerator_mode") != "kain_per_root":
        print(
            "Expected accelerator_mode=kain_per_root, got",
            protocol_node.get("accelerator_mode"),
        )
        success = False
    if int(protocol_node.get("accelerator_subspace", -1)) != 4:
        print(
            "Expected accelerator_subspace=4, got",
            protocol_node.get("accelerator_subspace"),
        )
        success = False
    if abs(float(protocol_node.get("density_convergence_target", -1.0)) - 1.0) > 1.0e-12:
        print(
            "Expected density_convergence_target=1.0, got",
            protocol_node.get("density_convergence_target"),
        )
        success = False
    if abs(float(protocol_node.get("relative_convergence_target", -1.0)) - 0.5) > 1.0e-12:
        print(
            "Expected relative_convergence_target=0.5, got",
            protocol_node.get("relative_convergence_target"),
        )
        success = False
    if abs(float(protocol_node.get("max_rotation", -1.0)) - 2.0) > 1.0e-12:
        print(
            "Expected max_rotation=2.0, got",
            protocol_node.get("max_rotation"),
        )
        success = False
    restart_source_threshold = protocol_node.get("restart_source_threshold")
    if not isinstance(restart_source_threshold, (int, float)) or abs(
        restart_source_threshold - 1.0e-2
    ) > 1.0e-12:
        print(
            "Expected restart_source_threshold=1e-2, got",
            restart_source_threshold,
        )
        success = False

    roots = protocol_node.get("roots", [])
    if isinstance(roots, list):
        seen_root_ids = set()
        slot_roots = {}
        for root in roots:
            if not isinstance(root, dict):
                print("Expected root descriptor objects in 'roots'")
                success = False
                continue
            root_id = root.get("root_id", "")
            if not isinstance(root_id, str) or not root_id:
                print("Expected non-empty root_id in root descriptor", root)
                success = False
            elif root_id in seen_root_ids:
                print("Expected unique root_id values, got duplicate", root_id)
                success = False
            else:
                seen_root_ids.add(root_id)
            slot_index = root.get("slot_index", root.get("root_index"))
            if not isinstance(slot_index, int) or slot_index < 0 or slot_index >= len(roots):
                print("Expected slot_index/root_index to be a valid slot", root)
                success = False
            elif slot_index in slot_roots:
                print("Expected unique slot_index/root_index values across roots", root)
                success = False
            else:
                slot_roots[slot_index] = root
            stable_index = root.get("stable_index", None)
            if not isinstance(stable_index, int):
                print("Expected integer stable_index in root descriptor", root)
                success = False
            display_name = root.get("display_name", root.get("name", ""))
            if not isinstance(display_name, str) or not display_name:
                print("Expected non-empty display name in root descriptor", root)
                success = False
            if not isinstance(root.get("energy"), (int, float)):
                print("Expected numeric energy in root descriptor", root)
                success = False
        if len(protocol_node.get("slot_permutation", [])) != len(roots):
            print("Expected slot_permutation length to match roots length")
            success = False
        elif len(slot_roots) != len(roots):
            print("Expected slot-indexed roots to cover every final slot")
            success = False
        else:
            derived_permutation = [
                slot_roots[slot]["stable_index"] for slot in range(len(roots))
            ]
            derived_names = [
                slot_roots[slot].get("display_name", slot_roots[slot].get("name", ""))
                for slot in range(len(roots))
            ]
            if protocol_node.get("slot_permutation") != derived_permutation:
                print(
                    "Expected slot_permutation to match root stable_index by slot, got",
                    protocol_node.get("slot_permutation"),
                    "expected",
                    derived_permutation,
                )
                success = False
            if protocol_node.get("state_names") != derived_names:
                print(
                    "Expected state_names to match root display names by slot, got",
                    protocol_node.get("state_names"),
                    "expected",
                    derived_names,
                )
                success = False
            derived_energies = [slot_roots[slot].get("energy") for slot in range(len(roots))]
            if protocol_node.get("energies") != derived_energies:
                print(
                    "Expected energies to match root energies by slot, got",
                    protocol_node.get("energies"),
                    "expected",
                    derived_energies,
                )
                success = False
    root_count = len(protocol_node.get("roots", []))
    for field_name in ("residual_norms", "density_change_norms", "relative_residual_norms"):
        field = protocol_node.get(field_name, [])
        if isinstance(field, list) and len(field) != root_count:
            print(f"Expected {field_name} length to match root count")
            success = False
    iter_count = protocol_node.get("iterations", 0)
    for field_name in (
        "iteration_max_residuals",
        "iteration_max_density_changes",
        "iteration_max_relative_residuals",
    ):
        field = protocol_node.get(field_name, [])
        if isinstance(field, list) and len(field) != iter_count:
            print(f"Expected {field_name} length to match iterations")
            success = False

    response_excited_states = response_metadata.get("excited_states", {})
    if response_excited_states != excited_states:
        print(
            "Expected excited_states subtree parity between calc_info metadata and response_metadata.json"
        )
        print("calc_info excited_states:", json.dumps(excited_states, indent=2))
        print(
            "response_metadata excited_states:",
            json.dumps(response_excited_states, indent=2),
        )
        success = False

    if not excited_exec.get("enabled", False):
        print("Expected excited execution.enabled=true")
        success = False
    if int(excited_exec.get("protocol_count", -1)) != 1:
        print(
            "Expected excited execution.protocol_count=1, got",
            excited_exec.get("protocol_count"),
        )
        success = False
    solver_adapter = excited_exec.get("solver_adapter", "")
    if not isinstance(solver_adapter, str) or not solver_adapter:
        print("Expected non-empty excited execution.solver_adapter")
        success = False

    events = excited_exec.get("protocol_events", [])
    if not isinstance(events, list) or len(events) == 0:
        print("Expected at least one excited protocol event")
        success = False
    else:
        event = events[0]
        if event.get("protocol_key") != protocol_key:
            print(
                "Expected first event protocol_key=",
                protocol_key,
                "got",
                event.get("protocol_key"),
            )
            success = False
        event_status = event.get("stage_status", "")
        if not isinstance(event_status, str) or not event_status:
            print("Expected non-empty stage_status in first protocol event")
            success = False
        if event.get("convergence_mode") != protocol_node.get("convergence_mode"):
            print("Expected event convergence_mode to match protocol metadata")
            success = False
        if event.get("accelerator_mode") != protocol_node.get("accelerator_mode"):
            print("Expected event accelerator_mode to match protocol metadata")
            success = False

    print("final success:", success)
    sys.exit(0 if success else 1)
