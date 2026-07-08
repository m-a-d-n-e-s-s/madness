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


def make_cmd(prefix, protocols, state_parallel_opts):
    protocol_list = ",".join(f"{value:.0e}" for value in protocols)
    return [
        "./@BINARY@",
        "--molecule=he",
        "--wf=response",
        f"--prefix={prefix}",
        f"--dft=k=6;maxiter=1;econv=1e-3;dconv=1e-2;protocol=[{protocol_list}];localize=canon",
        "--response=dipole=true;dipole.directions=z;dipole.frequencies=[0.0];"
        "excited.enable=true;excited.tda=true;excited.num_states=1;"
        "excited.guess_max_iter=1;excited.maxiter=2;excited.maxsub=4;"
        "print_level=2;"
        + state_parallel_opts,
    ]


def run_case(prefix, first_parallel_opts, second_parallel_opts, description, env):
    outputfile = prefix + ".calc_info.json"
    response_metadata_file = os.path.join(
        prefix, "task_1", "molresponse", "response_metadata.json"
    )
    cleanup(prefix)
    try:
        os.remove(outputfile)
    except FileNotFoundError:
        pass

    run1, _, _ = run_cmd(make_cmd(prefix, [1.0e-2], first_parallel_opts), env=env)
    if run1.returncode != 0:
        print("First run failed for", description)
        return False

    run2, stdout2, stderr2 = run_cmd(make_cmd(prefix, [1.0e-2], second_parallel_opts), env=env)
    if run2.returncode != 0:
        print("Second run failed for", description)
        return False

    if "EXCITED_INIT strategy=protocol_restart_guess" not in (stdout2 + "\n" + stderr2):
        print("Expected protocol_restart_guess initialization in second run for", description)
        return False

    if not os.path.exists(outputfile):
        print("Expected calc_info output file for", description)
        return False
    if not os.path.exists(response_metadata_file):
        print("Expected response_metadata.json for", description)
        return False

    with open(outputfile, "r", encoding="utf-8") as f:
        calc_info = json.load(f)
    with open(response_metadata_file, "r", encoding="utf-8") as f:
        response_metadata = json.load(f)

    response_task = find_response_task(calc_info)
    if response_task is None:
        print("Did not find response task in calc_info.json for", description)
        return False

    metadata = response_task.get("metadata", {})
    excited_states = metadata.get("excited_states", {})
    if response_metadata.get("excited_states", {}) != excited_states:
        print("Expected excited_states subtree parity for", description)
        return False

    protocol_key = f"{1.0e-2:.0e}"
    protocol_node = excited_states.get("protocols", {}).get(protocol_key, {})
    if protocol_node.get("restart_source") != "current_protocol_snapshot":
        print(
            "Expected current_protocol_snapshot restart_source for",
            description,
            "got",
            protocol_node.get("restart_source"),
        )
        return False
    if protocol_node.get("response_variant") != "static_restricted":
        print(
            "Expected static_restricted response_variant for",
            description,
            "got",
            protocol_node.get("response_variant"),
        )
        return False
    if protocol_node.get("snapshot_kind") != "protocol_bundle":
        print(
            "Expected protocol_bundle snapshot_kind for",
            description,
            "got",
            protocol_node.get("snapshot_kind"),
        )
        return False
    if protocol_node.get("bundle_state_present") is not True:
        print("Expected bundle_state_present=true for", description)
        return False
    if protocol_node.get("restart_capable") is not True:
        print("Expected restart_capable=true for", description)
        return False
    restart_source_threshold = protocol_node.get("restart_source_threshold")
    if not isinstance(restart_source_threshold, (int, float)) or abs(
        restart_source_threshold - 1.0e-2
    ) > 1.0e-12:
        print(
            "Expected restart_source_threshold=1e-2 for",
            description,
            "got",
            restart_source_threshold,
        )
        return False

    events = metadata.get("excited_state_planner", {}).get("execution", {}).get(
        "protocol_events", []
    )
    if not isinstance(events, list) or len(events) != 1:
        print("Expected exactly one protocol event for", description)
        return False
    if events[0].get("restart_source") != "current_protocol_snapshot":
        print(
            "Expected event restart_source=current_protocol_snapshot for",
            description,
            "got",
            events[0].get("restart_source"),
        )
        return False
    return True


if __name__ == "__main__":
    print("Testing @BINARY@/@TESTCASE@")
    print("reference files found in directory: @SRCDIR@")

    env = os.environ.copy()
    if env.get("SLURM_JOB_ID"):
        host = os.uname().nodename
        env["MADQC_LAUNCHER"] = (
            f"mpirun_rsh -launcher srun -np 2 {host} {host}"
        )
    else:
        env["MADQC_LAUNCHER"] = "mpirun -np 2"

    success = True
    success = run_case(
        "mad_@BINARY@_@TESTCASE@_serial_to_grouped",
        "state_parallel=off",
        "state_parallel=on;state_parallel_groups=2;state_parallel_min_states=1;"
        "state_parallel_point_start_protocol=0",
        "serial_to_grouped",
        env,
    ) and success
    success = run_case(
        "mad_@BINARY@_@TESTCASE@_grouped_to_serial",
        "state_parallel=on;state_parallel_groups=2;state_parallel_min_states=1;"
        "state_parallel_point_start_protocol=0",
        "state_parallel=off",
        "grouped_to_serial",
        env,
    ) and success

    print("final success:", success)
    sys.exit(0 if success else 1)
