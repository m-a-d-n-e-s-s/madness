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
    return p, stdout, stderr


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
        "--dft=k=6;maxiter=1;econv=1e-3;dconv=1e-2;protocol=[1e-1,1e-4];localize=canon",
        "--response=dipole=true;dipole.directions=z;dipole.frequencies=[0.0];"
        "excited.enable=true;excited.num_states=2;excited.guess_max_iter=1;"
        "excited.maxiter=2;excited.maxsub=4;print_level=2",
    ]
    run, stdout, stderr = run_cmd(cmd, env=env)
    if run.returncode != 0:
        sys.exit(run.returncode)

    success = True
    combined_output = stdout + "\n" + stderr
    if "EXCITED_BUNDLE_REPROJECT" not in combined_output:
        print("Expected EXCITED_BUNDLE_REPROJECT marker in solver output")
        success = False

    if not os.path.exists(outputfile):
        print("Output file not found:", outputfile)
        sys.exit(1)

    with open(outputfile, "r", encoding="utf-8") as f:
        got = json.load(f)
    response_metadata_file = os.path.join(
        prefix, "task_1", "molresponse", "response_metadata.json"
    )
    if not os.path.exists(response_metadata_file):
        print("Expected response_metadata.json at", response_metadata_file)
        sys.exit(1)
    with open(response_metadata_file, "r", encoding="utf-8") as f:
        response_metadata = json.load(f)

    response_task = find_response_task(got)
    if response_task is None:
        print("Did not find response task in calc_info.json")
        sys.exit(1)

    metadata = response_task.get("metadata", {})
    excited_states = metadata.get("excited_states", {})
    if response_metadata.get("excited_states", {}) != excited_states:
        print("Expected excited_states subtree parity for protocol projection run")
        success = False
    excited_protocols = excited_states.get("protocols", {})
    excited_exec = metadata.get("excited_state_planner", {}).get("execution", {})
    events = excited_exec.get("protocol_events", [])

    expected_protocol_keys = [f"{1.0e-1:.0e}", f"{1.0e-4:.0e}"]
    for key in expected_protocol_keys:
        if key not in excited_protocols:
            print("Expected excited protocol key", key, "in metadata")
            success = False
            continue
        protocol_node = excited_protocols[key]
        if protocol_node.get("response_variant") != "dynamic_restricted":
            print(
                "Expected response_variant=dynamic_restricted for protocol",
                key,
                "got",
                protocol_node.get("response_variant"),
            )
            success = False
        if protocol_node.get("restart_support_mode") != "full_bundle_resume":
            print(
                "Expected restart_support_mode=full_bundle_resume for protocol",
                key,
                "got",
                protocol_node.get("restart_support_mode"),
            )
            success = False
        if protocol_node.get("snapshot_kind") != "protocol_bundle":
            print(
                "Expected snapshot_kind=protocol_bundle for protocol",
                key,
                "got",
                protocol_node.get("snapshot_kind"),
            )
            success = False
        if protocol_node.get("bundle_state_present") is not True:
            print(
                "Expected bundle_state_present=true for protocol",
                key,
                "got",
                protocol_node.get("bundle_state_present"),
            )
            success = False
        if protocol_node.get("restart_capable") is not True:
            print(
                "Expected restart_capable=true for protocol",
                key,
                "got",
                protocol_node.get("restart_capable"),
            )
            success = False
        roots = protocol_node.get("roots", [])
        slot_permutation = protocol_node.get("slot_permutation", [])
        energies = protocol_node.get("energies", [])
        if not isinstance(roots, list) or len(roots) != len(energies):
            print(
                "Expected roots list with same length as energies for protocol",
                key,
            )
            success = False
        if not isinstance(slot_permutation, list) or len(slot_permutation) != len(roots):
            print(
                "Expected slot_permutation length to match roots for protocol",
                key,
            )
            success = False
        slot_roots = {}
        previous_stable_index = None
        for root in roots:
            if not isinstance(root, dict):
                print("Expected dict root descriptors for protocol", key)
                success = False
                continue
            stable_index = root.get("stable_index")
            if not isinstance(stable_index, int):
                print("Expected integer stable_index in root descriptor", root)
                success = False
            elif previous_stable_index is not None and stable_index < previous_stable_index:
                print("Expected roots manifest to be stable-index sorted for protocol", key)
                success = False
            previous_stable_index = stable_index
            slot_index = root.get("slot_index", root.get("root_index"))
            if not isinstance(slot_index, int) or slot_index < 0 or slot_index >= len(roots):
                print(
                    "Expected slot_index/root_index to be a valid final slot",
                    root,
                )
                success = False
            elif slot_index in slot_roots:
                print("Expected unique slot_index/root_index values for protocol", key, root)
                success = False
            else:
                slot_roots[slot_index] = root
            display_name = root.get("display_name", root.get("name", ""))
            if not isinstance(display_name, str) or not display_name:
                print("Expected non-empty display_name/name in root descriptor", root)
                success = False
            energy = root.get("energy")
            if not isinstance(energy, (int, float)):
                print("Expected numeric energy in root descriptor", root)
                success = False
        if len(slot_roots) != len(roots):
            print("Expected roots to cover every final slot for protocol", key)
            success = False
        else:
            derived_permutation = [
                slot_roots[slot]["stable_index"] for slot in range(len(roots))
            ]
            derived_names = [
                slot_roots[slot].get("display_name", slot_roots[slot].get("name", ""))
                for slot in range(len(roots))
            ]
            derived_energies = [slot_roots[slot].get("energy") for slot in range(len(roots))]
            if slot_permutation != derived_permutation:
                print(
                    "Expected slot_permutation to match stable_index by slot for protocol",
                    key,
                )
                print("slot_permutation:", slot_permutation)
                print("derived:", derived_permutation)
                success = False
            if protocol_node.get("state_names", []) != derived_names:
                print("Expected state_names to match root display names by slot for protocol", key)
                success = False
            if energies != derived_energies:
                print("Expected energies to match root energies by slot for protocol", key)
                success = False

    if int(excited_exec.get("protocol_count", -1)) != 2:
        print(
            "Expected excited execution.protocol_count=2, got",
            excited_exec.get("protocol_count"),
        )
        success = False

    if not isinstance(events, list) or len(events) != 2:
        print(
            "Expected exactly two excited protocol events, got",
            len(events) if isinstance(events, list) else "non-list",
        )
        success = False
    else:
        got_keys = [events[0].get("protocol_key", ""), events[1].get("protocol_key", "")]
        if got_keys != expected_protocol_keys:
            print("Expected protocol event keys", expected_protocol_keys, "got", got_keys)
            success = False
        if events[1].get("restart_source") != "lower_protocol_snapshot":
            print(
                "Expected second protocol event restart_source=lower_protocol_snapshot, got",
                events[1].get("restart_source"),
            )
            success = False
        if not events[1].get("restart_reused", False):
            print("Expected second protocol event restart_reused=true")
            success = False

    second_protocol = excited_protocols.get(expected_protocol_keys[1], {})
    if second_protocol.get("restart_source") != "lower_protocol_snapshot":
        print(
            "Expected second protocol restart_source=lower_protocol_snapshot, got",
            second_protocol.get("restart_source"),
        )
        success = False
    restart_source_threshold = second_protocol.get("restart_source_threshold")
    if not isinstance(restart_source_threshold, (int, float)) or abs(
        restart_source_threshold - 1.0e-1
    ) > 1.0e-12:
        print(
            "Expected second protocol restart_source_threshold=1e-1, got",
            restart_source_threshold,
        )
        success = False

    print("final success:", success)
    sys.exit(0 if success else 1)
