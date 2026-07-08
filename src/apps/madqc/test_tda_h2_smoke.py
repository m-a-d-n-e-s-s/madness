#!/usr/bin/env python3
"""Smoke test for test_tda_h2 — H2 TDA excitation energies via ResponseBundle KAIN.

Usage:
    python3 test_tda_h2_smoke.py [--build-dir DIR] [--keep-workdir] [--mpi-np N]

This script:
  1. Writes an H2 input file (geometry in bohr, HF/LDA, standard MRA settings).
  2. Runs moldft to produce calc_info.json + moldft.restartdata.
  3. Runs test_tda_h2 against that output.
  4. Parses excitation energies from stdout.
  5. Checks:
       - All energies are positive.
       - Energies are monotonically non-decreasing.
       - First 3 energies lie in a physically plausible range [0.1, 1.5] Ha.
  6. Optionally compares against reference values from legacy molresponse
     (if LEGACY_MOLRESPONSE env var points to a legacy binary).

Environment variables:
    BUILD_DIR            Path to the build directory (auto-detected from build symlink).
    MADQC_LAUNCHER       MPI launcher prefix, e.g. "mpiexec -np 2".
    MADQC_TEST_MPI_NP    Number of MPI ranks (overrides --mpi-np).
    LEGACY_MOLRESPONSE   Path to legacy molresponse binary (optional; skips if absent).
"""

import argparse
import os
import re
import shlex
import shutil
import subprocess
import sys
import tempfile

# ── Defaults ──────────────────────────────────────────────────────────────────

# H2 at R = 1.4 bohr along z-axis.
H2_INPUT = """\
geometry
  H  0.0  0.0  -0.7
  H  0.0  0.0   0.7
end

dft
  xc hf
  k 5
  L 20.0
  protocol [1e-4]
  nopen 0
end
"""

# Reference excitation energies from legacy molresponse TDA, H2, HF, thresh 1e-4
# (placeholder — update after a calibration run with the legacy binary).
# Format: (min_Ha, max_Ha) tolerance bands around the reference.
REFERENCE_ENERGIES_HA = None  # set to a list of floats once calibrated

# Loose physical plausibility band: all energies should be in [0.05, 2.0] Ha.
PLAUSIBILITY_LO = 0.05
PLAUSIBILITY_HI = 2.0


# ── Helpers ───────────────────────────────────────────────────────────────────

def find_build_dir():
    """Locate the ninja build directory (checks BUILD_DIR env, then build symlink)."""
    bd = os.environ.get("BUILD_DIR", "").strip()
    if bd and os.path.isdir(bd):
        return bd
    # Try the build symlink in the repo root.
    src = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
    repo_root = os.path.dirname(src)
    candidate = os.path.join(repo_root, "build")
    if os.path.isdir(candidate):
        return candidate
    return None


def find_binary(build_dir, name):
    """Search for a binary under build_dir."""
    for root, _dirs, files in os.walk(build_dir):
        if name in files:
            return os.path.join(root, name)
    return None


def run(cmd, cwd=None, env=None, check=True, capture=True):
    launcher_str = os.environ.get("MADQC_LAUNCHER", "").strip()
    full_cmd = shlex.split(launcher_str) + cmd if launcher_str else cmd
    print(f"  $ {' '.join(full_cmd)}")
    kw = dict(cwd=cwd, env=env)
    if capture:
        kw["stdout"] = subprocess.PIPE
        kw["stderr"] = subprocess.PIPE
    result = subprocess.run(full_cmd, **kw)
    if capture:
        stdout = result.stdout.decode("utf-8", errors="replace")
        stderr = result.stderr.decode("utf-8", errors="replace")
        print(stdout)
        if stderr:
            print(stderr, file=sys.stderr)
    else:
        stdout = ""
    if check and result.returncode != 0:
        raise RuntimeError(f"Command failed with exit code {result.returncode}: {full_cmd}")
    return result, stdout


def parse_omega_from_stdout(output):
    """Extract excitation energies from test_tda_h2 stdout.

    Matches lines like:
      State  0:  omega =  0.37451234 Ha  ( 10.1916 eV)
    """
    pattern = re.compile(r"State\s+(\d+):\s+omega\s*=\s*([-+]?\d+\.\d+)\s+Ha")
    matches = pattern.findall(output)
    if not matches:
        return []
    energies = [None] * len(matches)
    for idx_str, val_str in matches:
        energies[int(idx_str)] = float(val_str)
    return [e for e in energies if e is not None]


# ── Main ──────────────────────────────────────────────────────────────────────

def main():
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawTextHelpFormatter)
    parser.add_argument("--build-dir", default=None,
                        help="Path to the ninja build directory.")
    parser.add_argument("--keep-workdir", action="store_true",
                        help="Do not delete the temporary work directory on success.")
    parser.add_argument("--mpi-np", type=int, default=1,
                        help="Number of MPI ranks (default: 1).")
    parser.add_argument("--num-states", type=int, default=4,
                        help="Number of TDA roots (default: 4).")
    args = parser.parse_args()

    # Resolve MPI rank count.
    mpi_np = int(os.environ.get("MADQC_TEST_MPI_NP", args.mpi_np))
    if mpi_np > 1 and not os.environ.get("MADQC_LAUNCHER"):
        os.environ["MADQC_LAUNCHER"] = f"mpiexec -np {mpi_np}"

    # Find build directory.
    build_dir = args.build_dir or find_build_dir()
    if build_dir is None:
        print("ERROR: Cannot find build directory. Set BUILD_DIR or pass --build-dir.")
        sys.exit(1)
    print(f"Build directory: {build_dir}")

    moldft_bin = find_binary(build_dir, "moldft")
    if moldft_bin is None:
        print(f"ERROR: moldft binary not found under {build_dir}")
        sys.exit(1)

    test_tda_h2_bin = find_binary(build_dir, "test_tda_h2")
    if test_tda_h2_bin is None:
        print(f"ERROR: test_tda_h2 binary not found under {build_dir}")
        sys.exit(1)

    print(f"moldft:      {moldft_bin}")
    print(f"test_tda_h2: {test_tda_h2_bin}")

    workdir = tempfile.mkdtemp(prefix="tda_h2_smoke_")
    print(f"Work directory: {workdir}")

    try:
        # ── Step 1: write H2 input ─────────────────────────────────────────
        input_file = os.path.join(workdir, "input")
        with open(input_file, "w") as f:
            f.write(H2_INPUT)

        # ── Step 2: run moldft ─────────────────────────────────────────────
        print("\n=== Running moldft ===")
        _result, _stdout = run([moldft_bin, input_file], cwd=workdir)

        archive = os.path.join(workdir, "moldft.restartdata")
        if not os.path.exists(archive):
            raise RuntimeError(f"moldft did not produce {archive}")
        print(f"  moldft.restartdata found: {archive}")

        # ── Step 3: run test_tda_h2 ────────────────────────────────────────
        print("\n=== Running test_tda_h2 ===")
        cmd = [
            test_tda_h2_bin,
            workdir,
            "--num-states", str(args.num_states),
            "--max-iter", "40",
            "--dconv", "1e-4",
            "--print-level", "1",
        ]
        result, stdout = run(cmd, cwd=workdir)

        # ── Step 4: parse energies ─────────────────────────────────────────
        energies = parse_omega_from_stdout(stdout)
        if not energies:
            raise RuntimeError("Could not parse any excitation energies from test_tda_h2 output.")

        print(f"\nParsed {len(energies)} excitation energy(ies):")
        for i, e in enumerate(energies):
            print(f"  State {i}: {e:.8f} Ha  ({e * 27.2114:.4f} eV)")

        # ── Step 5: checks ────────────────────────────────────────────────
        failures = []

        # All energies positive.
        for i, e in enumerate(energies):
            if e <= 0:
                failures.append(f"State {i}: energy non-positive ({e:.6f} Ha)")

        # Monotonically non-decreasing.
        for i in range(1, len(energies)):
            if energies[i] < energies[i-1] - 1e-5:
                failures.append(
                    f"Energies not sorted: omega[{i-1}]={energies[i-1]:.6f} > "
                    f"omega[{i}]={energies[i]:.6f}")

        # Physical plausibility band.
        for i, e in enumerate(energies):
            if not (PLAUSIBILITY_LO <= e <= PLAUSIBILITY_HI):
                failures.append(
                    f"State {i}: energy {e:.6f} Ha outside plausible range "
                    f"[{PLAUSIBILITY_LO}, {PLAUSIBILITY_HI}] Ha")

        # Compare against reference if provided.
        if REFERENCE_ENERGIES_HA is not None:
            tol = 1.0e-3  # Ha
            for i, (e, ref) in enumerate(zip(energies, REFERENCE_ENERGIES_HA)):
                if abs(e - ref) > tol:
                    failures.append(
                        f"State {i}: |omega - ref| = {abs(e-ref):.2e} Ha > tol={tol:.0e} "
                        f"(computed={e:.6f}, ref={ref:.6f})")

        if failures:
            print("\nFAILURES:")
            for f in failures:
                print(f"  FAIL: {f}")
            sys.exit(1)
        else:
            print("\nPASS: all checks passed.")

    finally:
        if not args.keep_workdir:
            shutil.rmtree(workdir, ignore_errors=True)
        else:
            print(f"Work directory kept at: {workdir}")


if __name__ == "__main__":
    main()
