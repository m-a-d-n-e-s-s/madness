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
        "excited.enable=true;excited.tda=true;excited.num_states=2;"
        "excited.guess_max_iter=1;excited.maxiter=2;excited.maxsub=4;"
        "print_level=2",
    ]

    run1, stdout1, stderr1 = run_cmd(cmd, env=env)
    if run1.returncode != 0:
        sys.exit(run1.returncode)
    combined_output1 = stdout1 + "\n" + stderr1
    if "EXCITED_FRESH_GUESS_SELECT" not in combined_output1:
        print("Expected explicit fresh-guess selection marker in first run output")
        sys.exit(1)
    if "guess_bundle" not in combined_output1:
        print("Expected fresh-guess selection to preserve a guess_bundle seed")
        sys.exit(1)

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
    first_root_ids = [
        root.get("root_id") for root in first_roots if isinstance(root, dict)
    ]

    guess_archive_file = os.path.join(
        prefix, "task_1", "molresponse", "response.excited_bundle.guess.restartdata"
    )
    guess_archive_pattern = os.path.join(
        prefix, "task_1", "molresponse", "*.excited_bundle.guess.restartdata*"
    )
    guess_archives = sorted(glob.glob(guess_archive_pattern))
    if not guess_archives:
        print("Expected guess archive artifacts, none found at", guess_archive_pattern)
        sys.exit(1)
    print("Found guess archive artifacts:", guess_archives)
    if not os.path.exists(guess_archive_file):
        print("Expected primary guess archive restartdata file at", guess_archive_file)
        sys.exit(1)
    with open(guess_archive_file, "r", encoding="utf-8") as f:
        guess_snapshot = json.load(f)
    if guess_snapshot.get("snapshot_kind") != "guess_bundle":
        print(
            "Expected guess archive snapshot_kind=guess_bundle, got",
            guess_snapshot.get("snapshot_kind"),
        )
        sys.exit(1)
    if guess_snapshot.get("bundle_state_present") is not True:
        print(
            "Expected guess archive bundle_state_present=true, got",
            guess_snapshot.get("bundle_state_present"),
        )
        sys.exit(1)
    if guess_snapshot.get("restart_capable") is not True:
        print(
            "Expected guess archive restart_capable=true, got",
            guess_snapshot.get("restart_capable"),
        )
        sys.exit(1)
    if guess_snapshot.get("response_variant") != "static_restricted":
        print(
            "Expected guess archive response_variant=static_restricted, got",
            guess_snapshot.get("response_variant"),
        )
        sys.exit(1)
    if guess_snapshot.get("restart_support_mode") != "full_bundle_resume":
        print(
            "Expected guess archive restart_support_mode=full_bundle_resume, got",
            guess_snapshot.get("restart_support_mode"),
        )
        sys.exit(1)
    guess_roots = guess_snapshot.get("roots", [])
    guess_names = guess_snapshot.get("state_names", [])
    guess_root_ids = [
        root.get("root_id") for root in guess_roots if isinstance(root, dict)
    ]
    if guess_root_ids != first_root_ids:
        print(
            "Expected guess archive root ids to match first protocol metadata, got",
            guess_root_ids,
            "expected",
            first_root_ids,
        )
        sys.exit(1)
    if guess_names != first_names:
        print(
            "Expected guess archive state_names to match first protocol metadata, got",
            guess_names,
            "expected",
            first_names,
        )
        sys.exit(1)

    protocol_snapshot_pattern = os.path.join(
        prefix,
        "task_1",
        "molresponse",
        f"*.excited_bundle.{protocol_key}.restartdata*",
    )
    protocol_snapshots = sorted(glob.glob(protocol_snapshot_pattern))
    if not protocol_snapshots:
        print(
            "Expected protocol restart artifacts before fallback test, none found at",
            protocol_snapshot_pattern,
        )
        sys.exit(1)
    for path in protocol_snapshots:
        os.remove(path)
    print("Removed protocol restart artifacts:", protocol_snapshots)

    with open(response_metadata_file, "r", encoding="utf-8") as f:
        response_metadata = json.load(f)
    protocol_node = (
        response_metadata.get("excited_states", {})
        .get("protocols", {})
        .get(protocol_key, {})
    )
    protocol_node["saved"] = False
    protocol_node["converged"] = False
    protocol_node["stage_status"] = "forced_guess_archive_fallback"
    with open(response_metadata_file, "w", encoding="utf-8") as f:
        json.dump(response_metadata, f, indent=2)
        f.write("\n")

    run2, stdout2, stderr2 = run_cmd(cmd, env=env)
    if run2.returncode != 0:
        sys.exit(run2.returncode)

    combined_output = stdout2 + "\n" + stderr2
    success = True
    if "EXCITED_INIT strategy=guess_archive" not in combined_output:
        print("Expected guess-archive initialization path in second run output")
        success = False
    if "strategy=protocol_restart_guess" in combined_output:
        print("Did not expect same-protocol restart guess path in fallback test")
        success = False
    if "strategy=lower_protocol_restart_guess" in combined_output:
        print("Did not expect lower-protocol restart guess path in fallback test")
        success = False

    if not os.path.exists(outputfile):
        print("Expected calc_info output file after second run")
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
    if response_metadata_second.get("excited_states", {}) != excited_states:
        print("Expected excited_states subtree parity after guess fallback run")
        success = False

    protocol_node = excited_states.get("protocols", {}).get(protocol_key, {})
    final_roots = protocol_node.get("roots", [])
    final_names = protocol_node.get("state_names", [])
    final_root_ids = [
        root.get("root_id") for root in final_roots if isinstance(root, dict)
    ]
    if protocol_node.get("restart_source") != "guess_archive":
        print(
            "Expected restart_source=guess_archive, got",
            protocol_node.get("restart_source"),
        )
        success = False
    if protocol_node.get("response_variant") != "static_restricted":
        print(
            "Expected response_variant=static_restricted, got",
            protocol_node.get("response_variant"),
        )
        success = False
    if protocol_node.get("restart_support_mode") != "full_bundle_resume":
        print(
            "Expected restart_support_mode=full_bundle_resume, got",
            protocol_node.get("restart_support_mode"),
        )
        success = False
    if protocol_node.get("snapshot_kind") != "protocol_bundle":
        print(
            "Expected snapshot_kind=protocol_bundle after fallback solve, got",
            protocol_node.get("snapshot_kind"),
        )
        success = False
    if protocol_node.get("bundle_state_present") is not True:
        print(
            "Expected bundle_state_present=true after fallback solve, got",
            protocol_node.get("bundle_state_present"),
        )
        success = False
    if protocol_node.get("restart_capable") is not True:
        print(
            "Expected restart_capable=true after fallback solve, got",
            protocol_node.get("restart_capable"),
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
    if final_root_ids != first_root_ids:
        print(
            "Expected guess-archive fallback to preserve root ids, got",
            final_root_ids,
            "expected",
            first_root_ids,
        )
        success = False
    if final_names != first_names:
        print(
            "Expected guess-archive fallback to preserve state_names, got",
            final_names,
            "expected",
            first_names,
        )
        success = False

    events = metadata.get("excited_state_planner", {}).get("execution", {}).get(
        "protocol_events", []
    )
    if not isinstance(events, list) or len(events) != 1:
        print("Expected exactly one excited protocol event after fallback rerun")
        success = False
    elif events[0].get("restart_source") != "guess_archive":
        print(
            "Expected event restart_source=guess_archive, got",
            events[0].get("restart_source"),
        )
        success = False

    print("final success:", success)
    sys.exit(0 if success else 1)
