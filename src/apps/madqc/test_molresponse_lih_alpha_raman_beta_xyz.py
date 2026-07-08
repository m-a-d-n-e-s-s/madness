#!/usr/bin/env python3

import argparse
import json
import os
import shlex
import subprocess
import sys

sys.path.append("@CMAKE_SOURCE_DIR@/bin")
from test_utilities import cleanup


def response_rows_to_map(rows):
    mapped = {}
    for row in rows:
        key = (
            row["property"],
            tuple(row["component"]),
            round(float(row["freqB"]), 12),
            None if "freqC" not in row else round(float(row["freqC"]), 12),
        )
        mapped[key] = None if "value" not in row else float(row["value"])
    return mapped


def get_scf_energy(task):
    if "properties" in task and "energy" in task["properties"]:
        return float(task["properties"]["energy"])
    if "scf_total_energy" in task:
        return float(task["scf_total_energy"])
    if "scf" in task and "scf_total_energy" in task["scf"]:
        return float(task["scf"]["scf_total_energy"])
    raise KeyError("Could not find SCF energy in task JSON")


def run_cmd(cmd, env=None):
    print("executing\n ", " ".join(cmd))
    p = subprocess.run(
        cmd,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
        encoding='utf-8', errors='replace',
        env=env,
    )
    print("finished run")
    print(p.stdout)
    if p.stderr:
        print(p.stderr)
    print("exitcode", p.returncode)
    return p


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="LiH response regression test (alpha/beta/raman, xyz directions)"
    )
    parser.add_argument(
        "--reference_directory",
        action="store",
        default="@SRCDIR@",
        help="directory containing test input and reference output",
    )
    args = parser.parse_args()

    print("Testing @BINARY@/@TESTCASE@")
    print("reference files found in directory:", args.reference_directory)

    prefix = "mad_@BINARY@_@TESTCASE@"
    outputfile = prefix + ".calc_info.json"
    inputfile = os.path.join(
        args.reference_directory, "test_molresponse_lih_alpha_raman_beta_xyz.in"
    )
    referencefile = os.path.join(args.reference_directory, prefix + ".calc_info.ref.json")

    if not os.path.exists(inputfile):
        print("Input file not found:", inputfile)
        sys.exit(1)
    if not os.path.exists(referencefile):
        print("Reference file not found:", referencefile)
        print("Skipping until expected output JSON is added.")
        sys.exit(77)

    cleanup(prefix)
    for stale in (outputfile,):
        try:
            os.remove(stale)
        except FileNotFoundError:
            pass

    env = os.environ.copy()
    launcher = env.get("MADQC_LAUNCHER", "").strip()

    rsp_cmd = [
        "./@BINARY@",
        "--wf=response",
        f"--prefix={prefix}",
        inputfile,
    ]
    if launcher:
        rsp_cmd = shlex.split(launcher) + rsp_cmd
    rsp_run = run_cmd(rsp_cmd, env=env)
    if rsp_run.returncode != 0:
        sys.exit(rsp_run.returncode)

    if not os.path.exists(outputfile):
        print("Response output file not found:", outputfile)
        sys.exit(1)

    with open(outputfile, "r", encoding="utf-8") as f:
        got = json.load(f)
    with open(referencefile, "r", encoding="utf-8") as f:
        ref = json.load(f)

    success = True

    if len(got["tasks"]) != len(ref["tasks"]):
        print("Mismatch in number of tasks:", len(got["tasks"]), len(ref["tasks"]))
        success = False

    if len(got["tasks"]) < 2 or len(ref["tasks"]) < 2:
        print("Expected at least two tasks (SCF + response) in output/reference.")
        sys.exit(1)

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
        print("SCF energy mismatch:", got_energy, ref_energy)
        success = False

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
            if abs(gval - rval) > 1e-4:
                print("Response value mismatch for key:", key, gval, rval)
                success = False

    print("final success:", success)
    sys.exit(0 if success else 1)
