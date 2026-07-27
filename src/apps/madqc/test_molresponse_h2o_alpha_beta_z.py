#!/usr/bin/env python3

import argparse
import json
import os
import sys
import subprocess

sys.path.append("@CMAKE_SOURCE_DIR@/bin")
from test_utilities import cleanup, skip_on_small_machines


def response_rows_to_map(props):
    """Flatten the v3 response_properties OBJECT (property -> protocol_key ->
    row(s)) into comparable {key: value} pairs. Alpha rows expand per tensor
    element; beta/raman rows are one value each. (The old v2 engine emitted a
    flat array — v2-shaped files are rejected so a stale reference reads as a
    shape error, not a silent pass.)"""
    if not isinstance(props, dict):
        raise ValueError(
            "response_properties is not the v3 object shape — stale v2 "
            "reference? Regenerate the .ref.json from a v3 run."
        )
    mapped = {}
    for prop, by_key in props.items():
        for pkey, rows in by_key.items():
            if not isinstance(rows, list):
                rows = [rows]
            for row in rows:
                if "alpha" in row:  # tensor row
                    dirs = row.get("directions", "")
                    omega = round(float(row.get("omega", 0.0)), 12)
                    m = row["alpha"]
                    for i, r in enumerate(m):
                        for j, v in enumerate(r):
                            di = dirs[i] if i < len(dirs) else str(i)
                            dj = dirs[j] if j < len(dirs) else str(j)
                            mapped[(prop, pkey, di, dj, omega)] = float(v)
                elif "beta" in row:  # beta / raman row
                    key = (
                        prop, pkey,
                        row.get("A"), row.get("B"), row.get("C"),
                        round(float(row.get("freq_b", 0.0)), 12),
                        round(float(row.get("freq_c", 0.0)), 12),
                    )
                    mapped[key] = float(row["beta"])
    return mapped


def get_scf_energy(task):
    # Accept both legacy and nested task schemas.
    if "properties" in task and "energy" in task["properties"]:
        return float(task["properties"]["energy"])
    if "scf_total_energy" in task:
        return float(task["scf_total_energy"])
    if "scf" in task and "scf_total_energy" in task["scf"]:
        return float(task["scf"]["scf_total_energy"])
    raise KeyError("Could not find SCF energy in task JSON")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="H2O response regression test (alpha/beta, z direction)"
    )
    parser.add_argument(
        "--reference_directory",
        action="store",
        default="@SRCDIR@",
        help="directory containing test input and reference output",
    )
    args = parser.parse_args()

    # This response test is expensive on small machines.
    try:
        if skip_on_small_machines():
            print("Skipping this verylong test on small machines")
            sys.exit(77)
    except Exception:
        print("Unable to evaluate machine size from MAD_NUM_THREADS, skipping test")
        sys.exit(77)

    print("Testing @BINARY@/@TESTCASE@")
    print("reference files found in directory:", args.reference_directory)

    prefix = "mad_@BINARY@_@TESTCASE@"
    outputfile = prefix + ".calc_info.json"
    inputfile = os.path.join(args.reference_directory, "test_molresponse_h2o_alpha_beta_z.in")
    referencefile = os.path.join(args.reference_directory, prefix + ".calc_info.ref.json")

    if not os.path.exists(inputfile):
        print("Input file not found:", inputfile)
        sys.exit(1)
    if not os.path.exists(referencefile):
        print("Reference file not found:", referencefile)
        print("Skipping until expected output JSON is added.")
        sys.exit(77)

    cleanup(prefix)

    cmd = f"./@BINARY@ --wf=response --prefix={prefix} {inputfile}"
    print("executing\n", cmd)
    p = subprocess.run(
        cmd,
        shell=True,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
        encoding='utf-8', errors='replace',
    )

    print("finished run")
    print(p.stdout)
    if p.stderr:
        print(p.stderr)
    exitcode = p.returncode
    print("exitcode", exitcode)
    if exitcode != 0:
        sys.exit(exitcode)

    with open(outputfile, "r", encoding="utf-8") as f:
        got = json.load(f)
    with open(referencefile, "r", encoding="utf-8") as f:
        ref = json.load(f)

    success = True

    # Basic structure checks
    if len(got["tasks"]) != len(ref["tasks"]):
        print("Mismatch in number of tasks:", len(got["tasks"]), len(ref["tasks"]))
        success = False

    if len(got["tasks"]) < 2 or len(ref["tasks"]) < 2:
        print("Expected at least two tasks (SCF + response) in output/reference.")
        sys.exit(1)

    # Ground-state checks
    got_scf = got["tasks"][0]
    ref_scf = ref["tasks"][0]
    got_model = got_scf.get("model", "scf")
    ref_model = ref_scf.get("model", "scf")
    if got_model != ref_model:
        print("SCF model mismatch:", got_model, ref_model)
        success = False
    got_energy = get_scf_energy(got_scf)
    ref_energy = get_scf_energy(ref_scf)
    if abs(got_energy - ref_energy) > 1e-4:
        print(
            "SCF energy mismatch:",
            got_energy,
            ref_energy,
        )
        success = False

    # Response checks
    got_rsp = got["tasks"][1]
    ref_rsp = ref["tasks"][1]
    if got_rsp["type"] != ref_rsp["type"]:
        print("Response type mismatch:", got_rsp["type"], ref_rsp["type"])
        success = False

    got_rows = response_rows_to_map(got_rsp["properties"]["response_properties"])
    ref_rows = response_rows_to_map(ref_rsp["properties"]["response_properties"])

    if set(got_rows.keys()) != set(ref_rows.keys()):
        print("Response property key mismatch")
        print("Only in output:", sorted(set(got_rows.keys()) - set(ref_rows.keys())))
        print("Only in reference:", sorted(set(ref_rows.keys()) - set(got_rows.keys())))
        success = False
    else:
        for key in sorted(ref_rows.keys()):
            gval = got_rows[key]
            rval = ref_rows[key]
            if gval is None and rval is None:
                continue
            if gval is None or rval is None:
                print("Missing response value for key:", key, gval, rval)
                success = False
                continue
            # Tolerance: 1e-4 absolute floor + 5e-4 RELATIVE. Solver
            # run-to-run reproducibility at the deck's dconv 1e-4 is ~1e-4
            # relative (cross-node thread-order noise shifts the stopping
            # iteration), so a flat 1e-4 abs on beta ~O(10) false-fails.
            if abs(gval - rval) > max(1e-4, 5e-4 * abs(rval)):
                print("Response value mismatch for key:", key, gval, rval)
                success = False

    print("final success:", success)

    sys.exit(0 if success else 1)
