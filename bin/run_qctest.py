#!/usr/bin/env python3
"""Run one quantum-chemistry example deck and check it against its reference.

A qctest case is a self-contained directory (see src/examples/qc/README.md):

    <case>/<case>.in                        the input deck
    <case>/run.sh                           one-liner invocation
    <case>/check.json                       result keys + tolerances
    <case>/reference/<case>.calc_info.json  compared numerically
    <case>/reference/<case>.out             for humans; never compared

The case is copied into the current working directory and run there, so the
source tree is never written to.  ctest supplies a per-case scratch directory;
for a manual run, cd into an empty directory first.

    bin/run_qctest.py --case src/examples/qc/scf_h2o_hf
    bin/run_qctest.py --case src/examples/qc/scf_h2o_hf --update   # regenerate reference

Exit codes: 0 pass, 1 fail, 77 skipped (ctest's SKIP_RETURN_CODE).
"""

import argparse
import json
import os
import shutil
import subprocess
import sys
from pathlib import Path

sys.path.append(str(Path(__file__).resolve().parent))
from test_utilities import madjsoncompare

SKIP = 77


def load_check(case, updating):
    """Read check.json, failing loudly rather than with a KeyError deep in the run.

    Under --update the checks are not needed yet: generating the reference first
    and then reading it to decide what is worth comparing is the natural order
    when adding a case.
    """
    check_file = case / "check.json"
    if not check_file.is_file():
        if updating:
            print(f"note: no check.json in {case} yet")
            return {}
        sys.exit(f"qctest: no check.json in {case}")
    with open(check_file) as f:
        check = json.load(f)
    if not check.get("checks") and not updating:
        sys.exit(f"qctest: {check_file} lists no checks")
    return check


def should_skip(requires):
    """Resource gates, mirroring test_utilities.skip_on_small_machines but declarative.

    The thread gate reads MAD_NUM_THREADS when it is set and otherwise falls back
    to the detected core count. `skip_on_small_machines` looks at the env var
    alone, so it reports "too small" on any machine that has not exported it --
    which is why those tests skip nearly everywhere.
    """
    min_threads = requires.get("threads")
    if min_threads is not None:
        env = os.environ.get("MAD_NUM_THREADS")
        threads = int(env) if env else (os.cpu_count() or 1)
        if threads < min_threads:
            return (f"needs {min_threads} threads, found {threads} "
                    f"({'MAD_NUM_THREADS' if env else 'detected cores'})")

    if requires.get("mpi"):
        mpiexec = os.environ.get("MPIEXEC", "mpiexec")
        if shutil.which(mpiexec) is None:
            return f"needs MPI, but '{mpiexec}' is not on PATH"

    for var in requires.get("env", []):
        if not os.environ.get(var):
            return f"needs environment variable {var}"

    return None


# Run artifacts that let a warm workdir short-circuit the calculation. ctest's
# per-case scratch directory is created once at configure time and reused by
# every later run (cmake/modules/AddQCTests.cmake), so without this a second
# `ctest -L qctest` reuses the first run's state: madqc finds a valid
# <label>.calc_info.json, returns NextAction::Ok and skips the SCF entirely. The
# case still passes -- against its own previous output -- in a fraction of the
# time, which would hide a regression from anyone iterating locally.
#
# The checkpoint is the actual lever (SCFApplication::run keys the skip on
# has_results); the archives are cleared too so the engine's own restart cannot
# warm-start either.
STATE_GLOBS = ("*.calc_info.json", "*.restartdata*", "*.restartaodata")


def clear_previous_state(workdir):
    """Delete run artifacts left by an earlier run in this workdir.

    Files only, recursively, matching known MADNESS run-artifact names -- never
    whole directories, since a manual run may share the cwd with the user's own
    files.
    """
    removed = 0
    for pattern in STATE_GLOBS:
        for stale in sorted(workdir.rglob(pattern)):
            if stale.is_file():
                stale.unlink()
                removed += 1
    if removed:
        print(f"cleared {removed} artifact(s) from a previous run in {workdir}")


def stage(case, workdir):
    """Copy the case into workdir, skipping the reference (which must stay read-only)."""
    if case == workdir:
        sys.exit(f"qctest: refusing to run inside the case directory {case}; "
                 "cd into an empty scratch directory first")
    for entry in sorted(case.iterdir()):
        if entry.name in ("reference", "check.json", "README.md"):
            continue
        if entry.is_dir():
            shutil.copytree(entry, workdir / entry.name, dirs_exist_ok=True)
        else:
            shutil.copy2(entry, workdir / entry.name)
    return workdir / "run.sh"


def run(runscript, logfile):
    """Execute run.sh through sh, so the checked-in one-liner is what gets tested.

    Output is streamed line by line to our stdout AND to logfile as it arrives,
    rather than captured and printed at exit. A correlated case can run for tens
    of minutes; buffering its output until the end leaves no way to tell a slow
    run from a wedged one except stack samples and file mtimes.
    """
    cmd = ["sh", str(runscript.name)]
    print("executing:", " ".join(cmd), "  (MADQC=%s)" % os.environ.get("MADQC", "madqc"))
    print("streaming output to:", logfile)
    sys.stdout.flush()

    with open(logfile, "w", buffering=1) as log:
        p = subprocess.Popen(cmd, stdout=subprocess.PIPE, stderr=subprocess.STDOUT,
                             encoding="utf-8", errors="replace", bufsize=1)
        for line in p.stdout:
            sys.stdout.write(line)
            sys.stdout.flush()
            log.write(line)
        p.wait()

    print("program exit code:", p.returncode)
    return p.returncode


def lookup(data, keys):
    for k in keys:
        data = data[k]
    return data


def compare(output, reference, checks):
    """Numeric/exact comparison of the declared keys. A missing key is a failure."""
    cmp = madjsoncompare(str(output), str(reference))
    for entry in checks:
        keys, tol = entry["key"], entry.get("tol", 0.0)
        try:
            cmp.compare(keys, tol)
        except (KeyError, IndexError, TypeError) as e:
            print(f"key {keys} not found in both files: {type(e).__name__}: {e}")
            cmp.success = False
            continue

        # A reference value no larger than its own tolerance makes the check
        # vacuous: every conceivable result passes, so it survives any
        # regression. Exactly 0.0 is the common case -- a field the workflow
        # never populated, which is what several of the older hand-written
        # scripted tests compare for `scf_total_energy` -- but a dark state's
        # 1e-26 oscillator strength checked to 1e-3 is just as empty. Require an
        # explicit opt-in where the near-zero is the physics (a
        # symmetry-vanishing dipole, a gradient at a stationary point).
        ref = lookup(cmp.data2, keys)
        if isinstance(ref, float) and abs(ref) <= tol and not entry.get("allow_zero"):
            print(f"key {keys} has reference value {ref!r}, within its own tolerance "
                  f"{tol} -- this check cannot fail. Compare a key that carries a "
                  "value, or set \"allow_zero\": true if the near-zero is physical.")
            cmp.success = False
    return cmp.success


def update_reference(case, output, report):
    refdir = case / "reference"
    refdir.mkdir(exist_ok=True)
    shutil.copy2(output, refdir / output.name)
    print("updated reference:", refdir / output.name)
    if report.is_file():
        shutil.copy2(report, refdir / report.name)
        print("updated reference:", refdir / report.name)


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--case", required=True, help="path to the qctest case directory")
    ap.add_argument("--update", action="store_true",
                    help="overwrite the checked-in reference with this run's results")
    ap.add_argument("--keep-state", action="store_true",
                    help="do not clear checkpoints/restart archives left by an earlier "
                         "run in the work directory; for a case that deliberately "
                         "exercises restart from warm state")
    args = ap.parse_args()

    case = Path(args.case).resolve()
    if not case.is_dir():
        sys.exit(f"qctest: no such case directory: {case}")
    workdir = Path.cwd().resolve()
    check = load_check(case, args.update)

    print(f"qctest case: {case.name}")
    print(f"description: {check.get('description', '(none)')}")
    print(f"case dir   : {case}")
    print(f"work dir   : {workdir}")

    reason = should_skip(check.get("requires", {}))
    if reason:
        print(f"SKIPPED: {reason}")
        return SKIP

    if not args.keep_state:
        clear_previous_state(workdir)

    runscript = stage(case, workdir)
    if not runscript.is_file():
        sys.exit(f"qctest: {case} has no run.sh")

    exitcode = run(runscript, workdir / f"{case.name}.run.log")
    if exitcode == SKIP:
        print("SKIPPED: run.sh returned 77")
        return SKIP

    # madqc names its output after the deck: <case>.in -> <case>.calc_info.json.
    # (A deck literally named `input` would fall back to the default prefix `mad`,
    # which is why cases name the deck after themselves.)
    output = workdir / check.get("output", f"{case.name}.calc_info.json")
    report = workdir / output.name.replace(".calc_info.json", ".out")

    if not output.is_file():
        print(f"FAILED: expected result file {output.name} was not produced")
        return 1

    if args.update:
        update_reference(case, output, report)
        return 0

    reference = case / "reference" / output.name
    if not reference.is_file():
        print(f"FAILED: no reference at {reference}; generate it with --update")
        return 1

    ok = compare(output, reference, check["checks"])
    print("final success:", ok and exitcode == 0)
    return 0 if (ok and exitcode == 0) else 1


if __name__ == "__main__":
    sys.exit(main())
