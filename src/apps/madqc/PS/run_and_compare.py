#!/usr/bin/env python3
"""
Run and compare madqc pseudopotential (GTH-LDA and GTH-PBE) benchmark calculations
against reference madqc outputs and NWChem nwpw (pspw) aperiodic plane-wave outputs.

Usage:
  # Compare existing reference outputs in madqc/ and nwchem/:
  python3 run_and_compare.py --compare-only

  # Run madqc on all benchmark inputs and compare against reference outputs:
  python3 run_and_compare.py --madqc /path/to/madqc

  # Run madqc on a subset of tests:
  python3 run_and_compare.py --madqc /path/to/madqc --tests be_lda ne_pbe water_lda water_pbe
"""

import argparse
import os
import re
import shutil
import subprocess
import sys
import tempfile

ALL_TESTS = [
    "be_lda",
    "be_pbe",
    "ne_lda",
    "ne_pbe",
    "mg_lda",
    "mg_pbe",
    "water_lda",
    "water_pbe",
    "water2_lda",
    "water2_pbe",
    "methane_lda",
    "methane_pbe",
    "ethylene_lda",
    "ethylene_pbe",
    "ethylene_allpsp_lda",
    "ethylene_allpsp_pbe",
]


def parse_madqc_output(filepath):
    """Parse the final converged SCF energy components and eigenvalues from a madqc output file."""
    if not os.path.exists(filepath):
        return None
    with open(filepath, "r", encoding="utf-8", errors="replace") as f:
        text = f.read()

    blocks = re.findall(
        r"(              kinetic\s+[-0-9.]+\n(?:.*?\n)*?                total\s+[-0-9.]+)",
        text,
    )
    if not blocks:
        return None

    last_block = blocks[-1]
    data = {}
    patterns = {
        "kinetic": r"kinetic\s+([-0-9.]+)",
        "nonlocal_psp": r"nonlocal psp\s+([-0-9.]+)",
        "nuclear_attraction": r"nuclear attraction\s+([-0-9.]+)",
        "coulomb": r"coulomb\s+([-0-9.]+)",
        "xc": r"exchange-correlation\s+([-0-9.]+)",
        "nuclear_repulsion": r"nuclear-repulsion\s+([-0-9.]+)",
        "total": r"total\s+([-0-9.]+)",
    }
    for key, pat in patterns.items():
        m = re.search(pat, last_block)
        if m:
            data[key] = float(m.group(1))

    evals = re.findall(r"alpha eigenvalues\n\[\*\]\s+([^\n]+)", text)
    if evals:
        data["eigenvalues"] = [float(x) for x in evals[-1].split()]
    return data


def parse_nwchem_output(filepath):
    """Parse the final converged PSPW energy components and eigenvalues from an NWChem output file."""
    if not os.path.exists(filepath):
        return None
    with open(filepath, "r", encoding="utf-8", errors="replace") as f:
        text = f.read()

    summaries = text.split("==  Summary Of Results  ==")
    if len(summaries) < 2:
        return None

    last_summary = summaries[-1]
    data = {}
    patterns = {
        "total": r"total\s+energy\s*:\s*([-+0-9.EdD]+)",
        "coulomb": r"hartree\s+energy\s*:\s*([-+0-9.EdD]+)",
        "xc": r"exc-corr\s+energy\s*:\s*([-+0-9.EdD]+)",
        "nuclear_repulsion": r"ion-ion\s+energy\s*:\s*([-+0-9.EdD]+)",
        "kinetic": r"kinetic \(planewave\)\s*:\s*([-+0-9.EdD]+)",
        "nuclear_attraction": r"V_local \(planewave\)\s*:\s*([-+0-9.EdD]+)",
        "nonlocal_psp": r"V_nl\s+\(planewave\)\s*:\s*([-+0-9.EdD]+)",
    }
    for key, pat in patterns.items():
        m = re.search(pat, last_summary)
        if m:
            data[key] = float(m.group(1).replace("D", "E").replace("d", "e"))
    return data


def find_madqc_binary(user_path=None):
    if user_path:
        if os.path.isfile(user_path) and os.access(user_path, os.X_OK):
            return os.path.abspath(user_path)
        raise FileNotFoundError(f"Specified madqc binary not found or not executable: {user_path}")

    script_dir = os.path.dirname(os.path.abspath(__file__))
    candidates = [
        shutil.which("madqc"),
        "/tmp/madrel-build/src/apps/madqc/madqc",
        os.path.join(script_dir, "..", "madqc"),
        os.path.join(script_dir, "..", "..", "..", "..", "build", "src", "apps", "madqc", "madqc"),
    ]
    for c in candidates:
        if c and os.path.isfile(c) and os.access(c, os.X_OK):
            return os.path.abspath(c)
    return None


def main():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--madqc", help="Path to madqc executable")
    parser.add_argument(
        "--compare-only",
        action="store_true",
        help="Only compare existing reference outputs in madqc/ and nwchem/ without running madqc",
    )
    parser.add_argument(
        "--tests",
        nargs="+",
        default=ALL_TESTS,
        choices=ALL_TESTS,
        help="Subset of tests to run/compare (default: all 14 tests)",
    )
    parser.add_argument(
        "--tol-ref",
        type=float,
        default=1.0e-5,
        help="Maximum allowed |E(run) - E(ref_madqc)| in Ha when running madqc (default: 1e-5)",
    )
    parser.add_argument(
        "--work-dir",
        help="Directory in which to run madqc calculations (default: temporary directory)",
    )
    parser.add_argument(
        "--threads",
        type=int,
        default=int(os.environ.get("MAD_NUM_THREADS", "4")),
        help="Number of MADNESS threads (MAD_NUM_THREADS) to use when running madqc to avoid core oversubscription (default: 4)",
    )
    args = parser.parse_args()

    script_dir = os.path.dirname(os.path.abspath(__file__))
    madqc_dir = os.path.join(script_dir, "madqc")
    nwchem_dir = os.path.join(script_dir, "nwchem")
    chemdata_dir = os.path.abspath(os.path.join(script_dir, "..", "..", "..", "madness", "chem"))

    run_results = {}
    if not args.compare_only:
        madqc_bin = find_madqc_binary(args.madqc)
        if not madqc_bin:
            print(
                "ERROR: Could not locate madqc binary. Pass --madqc /path/to/madqc or use --compare-only.",
                file=sys.stderr,
            )
            return 1

        env = os.environ.copy()
        if "MRA_CHEMDATA_DIR" not in env and os.path.isdir(chemdata_dir):
            env["MRA_CHEMDATA_DIR"] = chemdata_dir
        env["MAD_NUM_THREADS"] = str(args.threads)
        env["OMP_NUM_THREADS"] = "1"
        env["MKL_NUM_THREADS"] = "1"

        work_root = args.work_dir or tempfile.mkdtemp(prefix="madqc_ps_test_")
        os.makedirs(work_root, exist_ok=True)
        print(f"Using madqc binary   : {madqc_bin}")
        print(f"Using MRA_CHEMDATA_DIR: {env.get('MRA_CHEMDATA_DIR', '(default)')}")
        print(f"Using MAD_NUM_THREADS: {env['MAD_NUM_THREADS']} (OMP_NUM_THREADS=1, MKL_NUM_THREADS=1)")
        print(f"Working directory    : {work_root}\n")

        for test in args.tests:
            inp_src = os.path.join(madqc_dir, f"{test}.in")
            job_dir = os.path.join(work_root, test)
            os.makedirs(job_dir, exist_ok=True)
            inp_dst = os.path.join(job_dir, f"{test}.in")
            out_dst = os.path.join(job_dir, f"{test}.out")
            shutil.copyfile(inp_src, inp_dst)

            print(f"Running {test:<14} ... ", end="", flush=True)
            with open(out_dst, "w", encoding="utf-8") as out_f:
                proc = subprocess.run([madqc_bin, f"{test}.in"], cwd=job_dir, env=env, stdout=out_f, stderr=subprocess.STDOUT)
            if proc.returncode != 0:
                print(f"FAILED (exit code {proc.returncode})")
                return 1
            parsed = parse_madqc_output(out_dst)
            if not parsed or "total" not in parsed:
                print("FAILED (could not parse total energy)")
                return 1
            run_results[test] = parsed
            print(f"done (E = {parsed['total']:14.8f} Ha)")
        print()

    # Print comparison table
    header = (
        f"{'Test':<20} {'E(madqc ref) [Ha]':>18} {'E(NWChem pspw) [Ha]':>20} "
        f"{'dE(madqc-NWChem)':>17} {'V_nl(madqc)':>13} {'V_nl(NWChem)':>13}"
    )
    if not args.compare_only:
        header += f" {'dE(run-ref)':>13} {'Status':>8}"
    print(header)
    print("-" * len(header))

    all_passed = True
    for test in args.tests:
        ref_madqc = parse_madqc_output(os.path.join(madqc_dir, f"{test}.out"))
        ref_nwchem = parse_nwchem_output(os.path.join(nwchem_dir, f"{test}.out"))

        e_mad = ref_madqc["total"] if ref_madqc and "total" in ref_madqc else float("nan")
        e_nwc = ref_nwchem["total"] if ref_nwchem and "total" in ref_nwchem else float("nan")
        de_nw = e_mad - e_nwc
        vnl_mad = ref_madqc.get("nonlocal_psp", float("nan")) if ref_madqc else float("nan")
        vnl_nwc = ref_nwchem.get("nonlocal_psp", float("nan")) if ref_nwchem else float("nan")

        row = f"{test:<20} {e_mad:18.8f} {e_nwc:20.8f} {de_nw:+17.8f} {vnl_mad:13.8f} {vnl_nwc:13.8f}"
        if not args.compare_only:
            e_run = run_results[test]["total"]
            de_ref = e_run - e_mad
            passed = abs(de_ref) <= args.tol_ref
            if not passed:
                all_passed = False
            status = "PASS" if passed else "FAIL"
            row += f" {de_ref:+13.8f} {status:>8}"
        print(row)

    print("-" * len(header))
    if not args.compare_only:
        if all_passed:
            print(f"All {len(args.tests)} madqc calculations matched reference outputs within {args.tol_ref:.1e} Ha.")
            return 0
        else:
            print(f"One or more calculations differed from reference outputs by > {args.tol_ref:.1e} Ha.", file=sys.stderr)
            return 1
    return 0


if __name__ == "__main__":
    sys.exit(main())
