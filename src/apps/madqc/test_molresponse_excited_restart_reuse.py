#!/usr/bin/env python3

import glob
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
        "excited.enable=true;excited.tda=true;excited.num_states=1;"
        "excited.guess_max_iter=1;excited.maxiter=2;excited.maxsub=4;"
        "print_level=2",
    ]

    run1 = run_cmd(cmd, env=env)
    if run1.returncode != 0:
        sys.exit(run1.returncode)

    if not os.path.exists(response_metadata_file):
        print("Expected response_metadata.json after first run")
        sys.exit(1)
    with open(response_metadata_file, "r", encoding="utf-8") as f:
        response_metadata_first = json.load(f)
    first_protocol_node = (
        response_metadata_first.get("excited_states", {})
        .get("protocols", {})
        .get(protocol_key, {})
    )
    first_roots = first_protocol_node.get("roots", [])
    first_names = first_protocol_node.get("state_names", [])
    first_slot_permutation = first_protocol_node.get("slot_permutation", [])
    first_root_ids = [root.get("root_id") for root in first_roots if isinstance(root, dict)]
    if first_protocol_node.get("response_variant") != "static_restricted":
        print(
            "Expected first-run response_variant=static_restricted, got",
            first_protocol_node.get("response_variant"),
        )
        sys.exit(1)
    if first_protocol_node.get("restart_support_mode") != "full_bundle_resume":
        print(
            "Expected first-run restart_support_mode=full_bundle_resume, got",
            first_protocol_node.get("restart_support_mode"),
        )
        sys.exit(1)
    if first_protocol_node.get("snapshot_kind") != "protocol_bundle":
        print(
            "Expected first-run snapshot_kind=protocol_bundle, got",
            first_protocol_node.get("snapshot_kind"),
        )
        sys.exit(1)
    if first_protocol_node.get("bundle_state_present") is not True:
        print(
            "Expected first-run bundle_state_present=true, got",
            first_protocol_node.get("bundle_state_present"),
        )
        sys.exit(1)
    if first_protocol_node.get("restart_capable") is not True:
        print(
            "Expected first-run restart_capable=true, got",
            first_protocol_node.get("restart_capable"),
        )
        sys.exit(1)

    restart_pattern = os.path.join(
        prefix, "task_1", "molresponse", "*.excited_bundle.*.restartdata"
    )
    restart_files = sorted(glob.glob(restart_pattern))
    if not restart_files:
        print("Expected excited restart snapshot, none found at", restart_pattern)
        sys.exit(1)
    print("Found excited restart snapshots:", restart_files)
    protocol_restart = None
    for candidate in restart_files:
        if f".{protocol_key}.restartdata" in candidate:
            protocol_restart = candidate
            break
    if protocol_restart is None:
        print("Expected protocol restart snapshot for", protocol_key)
        sys.exit(1)
    with open(protocol_restart, "r", encoding="utf-8") as f:
        restart_snapshot = json.load(f)
    if restart_snapshot.get("response_variant") != "static_restricted":
        print(
            "Expected protocol restart snapshot response_variant=static_restricted, got",
            restart_snapshot.get("response_variant"),
        )
        sys.exit(1)
    if restart_snapshot.get("restart_support_mode") != "full_bundle_resume":
        print(
            "Expected protocol restart snapshot restart_support_mode=full_bundle_resume, got",
            restart_snapshot.get("restart_support_mode"),
        )
        sys.exit(1)
    if restart_snapshot.get("snapshot_kind") != "protocol_bundle":
        print(
            "Expected protocol restart snapshot kind=protocol_bundle, got",
            restart_snapshot.get("snapshot_kind"),
        )
        sys.exit(1)
    if restart_snapshot.get("bundle_state_present") is not True:
        print(
            "Expected protocol restart snapshot bundle_state_present=true, got",
            restart_snapshot.get("bundle_state_present"),
        )
        sys.exit(1)
    if restart_snapshot.get("restart_capable") is not True:
        print(
            "Expected protocol restart snapshot restart_capable=true, got",
            restart_snapshot.get("restart_capable"),
        )
        sys.exit(1)
    if restart_snapshot.get("convergence_mode") != "density_relative_dual_gate":
        print(
            "Expected protocol restart snapshot convergence_mode=density_relative_dual_gate, got",
            restart_snapshot.get("convergence_mode"),
        )
        sys.exit(1)
    if restart_snapshot.get("accelerator_mode") != "kain_per_root":
        print(
            "Expected protocol restart snapshot accelerator_mode=kain_per_root, got",
            restart_snapshot.get("accelerator_mode"),
        )
        sys.exit(1)
    if int(restart_snapshot.get("accelerator_subspace", -1)) != 4:
        print(
            "Expected protocol restart snapshot accelerator_subspace=4, got",
            restart_snapshot.get("accelerator_subspace"),
        )
        sys.exit(1)

    run2 = run_cmd(cmd, env=env)
    if run2.returncode != 0:
        sys.exit(run2.returncode)

    if not os.path.exists(outputfile):
        print("Output file not found:", outputfile)
        sys.exit(1)
    if not os.path.exists(response_metadata_file):
        print("response_metadata.json missing after second run")
        sys.exit(1)

    with open(outputfile, "r", encoding="utf-8") as f:
        got = json.load(f)
    with open(response_metadata_file, "r", encoding="utf-8") as f:
        response_metadata_second = json.load(f)

    response_task = find_response_task(got)
    if response_task is None:
        print("Did not find response task in calc_info.json")
        sys.exit(1)

    metadata = response_task.get("metadata", {})
    excited_states = metadata.get("excited_states", {})
    protocol_node = excited_states.get("protocols", {}).get(protocol_key, {})
    if response_metadata_second.get("excited_states", {}) != excited_states:
        print("Expected excited_states subtree parity after restart reuse")
        success = False
    excited_exec = metadata.get("excited_state_planner", {}).get("execution", {})
    events = excited_exec.get("protocol_events", [])

    success = True

    if not isinstance(events, list) or len(events) == 0:
        print("Expected at least one excited protocol event in second run")
        success = False
        event = {}
    else:
        event = events[0]

    restart_reused = bool(event.get("restart_reused", False))
    if not restart_reused:
        print("Expected restart_reused=true in first protocol event")
        success = False
    if event.get("restart_source") != "current_protocol_snapshot":
        print(
            "Expected event restart_source=current_protocol_snapshot, got",
            event.get("restart_source"),
        )
        success = False

    stage_status = protocol_node.get("stage_status", "")
    if not isinstance(stage_status, str) or "restart" not in stage_status.lower():
        print(
            "Expected protocol stage_status to indicate restart reuse, got",
            stage_status,
        )
        success = False
    if protocol_node.get("restart_source") != "current_protocol_snapshot":
        print(
            "Expected protocol restart_source=current_protocol_snapshot, got",
            protocol_node.get("restart_source"),
        )
        success = False
    if protocol_node.get("snapshot_kind") != "protocol_bundle":
        print(
            "Expected protocol snapshot_kind=protocol_bundle after restart reuse, got",
            protocol_node.get("snapshot_kind"),
        )
        success = False
    if protocol_node.get("bundle_state_present") is not True:
        print(
            "Expected protocol bundle_state_present=true after restart reuse, got",
            protocol_node.get("bundle_state_present"),
        )
        success = False
    if protocol_node.get("restart_capable") is not True:
        print(
            "Expected protocol restart_capable=true after restart reuse, got",
            protocol_node.get("restart_capable"),
        )
        success = False
    if protocol_node.get("response_variant") != "static_restricted":
        print(
            "Expected protocol response_variant=static_restricted after restart reuse, got",
            protocol_node.get("response_variant"),
        )
        success = False
    if protocol_node.get("restart_support_mode") != "full_bundle_resume":
        print(
            "Expected protocol restart_support_mode=full_bundle_resume after restart reuse, got",
            protocol_node.get("restart_support_mode"),
        )
        success = False
    if protocol_node.get("convergence_mode") != "density_relative_dual_gate":
        print(
            "Expected protocol convergence_mode=density_relative_dual_gate after restart reuse, got",
            protocol_node.get("convergence_mode"),
        )
        success = False
    if protocol_node.get("accelerator_mode") != "kain_per_root":
        print(
            "Expected protocol accelerator_mode=kain_per_root after restart reuse, got",
            protocol_node.get("accelerator_mode"),
        )
        success = False
    if int(protocol_node.get("accelerator_subspace", -1)) != 4:
        print(
            "Expected protocol accelerator_subspace=4 after restart reuse, got",
            protocol_node.get("accelerator_subspace"),
        )
        success = False
    for field_name in (
        "density_change_norms",
        "relative_residual_norms",
        "iteration_max_density_changes",
        "iteration_max_relative_residuals",
    ):
        if not isinstance(protocol_node.get(field_name), list):
            print(f"Expected list {field_name} after restart reuse")
            success = False
    restart_source_threshold = protocol_node.get("restart_source_threshold")
    if not isinstance(restart_source_threshold, (int, float)) or abs(
        restart_source_threshold - 1.0e-2
    ) > 1.0e-12:
        print(
            "Expected protocol restart_source_threshold=1e-2, got",
            restart_source_threshold,
        )
        success = False

    restart_ready_protocols = int(excited_exec.get("restart_ready_protocols", 0))
    if restart_ready_protocols < 1:
        print(
            "Expected execution.restart_ready_protocols >= 1, got",
            restart_ready_protocols,
        )
        success = False

    second_protocol_node = (
        response_metadata_second.get("excited_states", {})
        .get("protocols", {})
        .get(protocol_key, {})
    )
    second_roots = second_protocol_node.get("roots", [])
    second_names = second_protocol_node.get("state_names", [])
    second_slot_permutation = second_protocol_node.get("slot_permutation", [])
    second_root_ids = [
        root.get("root_id") for root in second_roots if isinstance(root, dict)
    ]
    first_root_name_map = {
        root.get("root_id"): root.get("display_name", root.get("name", ""))
        for root in first_roots
        if isinstance(root, dict) and root.get("root_id")
    }
    second_root_name_map = {
        root.get("root_id"): root.get("display_name", root.get("name", ""))
        for root in second_roots
        if isinstance(root, dict) and root.get("root_id")
    }

    if first_root_ids != second_root_ids:
        print("Expected stable root ids across no-op restart")
        print("first_root_ids:", first_root_ids)
        print("second_root_ids:", second_root_ids)
        success = False

    if first_root_name_map != second_root_name_map:
        print("Expected stable root display-name mapping across no-op restart")
        print("first_root_name_map:", first_root_name_map)
        print("second_root_name_map:", second_root_name_map)
        success = False

    if isinstance(second_roots, list):
        slot_roots = {}
        previous_stable_index = None
        for root in second_roots:
            if not isinstance(root, dict):
                print("Expected dict root descriptors in second run roots")
                success = False
                continue
            stable_index = root.get("stable_index")
            if isinstance(stable_index, int):
                if previous_stable_index is not None and stable_index < previous_stable_index:
                    print("Expected second-run roots manifest to be stable-index sorted")
                    success = False
                previous_stable_index = stable_index
            slot_index = root.get("slot_index", root.get("root_index"))
            if not isinstance(slot_index, int) or slot_index < 0 or slot_index >= len(second_roots):
                print(
                    "Expected root slot_index/root_index to be a valid final slot",
                    root,
                )
                success = False
                continue
            if slot_index in slot_roots:
                print(
                    "Expected unique slot_index/root_index values in second run roots",
                    root,
                )
                success = False
                continue
            slot_roots[slot_index] = root
        if len(slot_roots) != len(second_roots):
            print("Expected second-run roots to cover every final slot")
            success = False
        else:
            derived_permutation = [
                slot_roots[slot]["stable_index"] for slot in range(len(second_roots))
            ]
            derived_names = [
                slot_roots[slot].get("display_name", slot_roots[slot].get("name", ""))
                for slot in range(len(second_roots))
            ]
            if derived_permutation != second_slot_permutation:
                print(
                    "Expected second-run slot_permutation to match roots by slot",
                    derived_permutation,
                    second_slot_permutation,
                )
                success = False
            if derived_names != second_names:
                print(
                    "Expected second-run state_names to match root display names by slot",
                    derived_names,
                    second_names,
                )
                success = False

    print("final success:", success)
    sys.exit(0 if success else 1)
