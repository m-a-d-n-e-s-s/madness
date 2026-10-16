#!/usr/bin/env python3
"""Reruns in one directory must answer the request, not the previous one (#822).

madqc decides whether the results of an earlier run can stand from what is on
disk: the restart planner reads the orbital archive's header, and the results
file is reused only if it was computed from that very archive (its archive_id)
with the same post-SCF requests. Each case below is a rerun that the old
checkpoint logic answered with the previous run's numbers:

  * a plain rerun          -- must reuse the results, without iterating
  * `restart none`         -- must iterate from the atomic guess
  * `nvalpha 1`            -- the archive has no virtual: must iterate
  * `eprec 1e-2` after 1e-4 -- a different Hamiltonian: must re-converge, and
                              the run must keep the requested eprec
  * `dipole true` added    -- the orbitals stand: properties only, no iterations

Helium at one coarse rung, so every run takes seconds. `--optimize` followed by `--wf=scf` is
covered by test_madqc_optimize_results.py.

    ./test_madqc_restart_bookkeeping.py     # after CMake substitution
"""

import glob
import json
import os
import subprocess
import sys

sys.path.append("@CMAKE_SOURCE_DIR@/bin")
from test_utilities import cleanup  # noqa: E402


def run(cmd):
    print("executing\n ", cmd)
    p = subprocess.run(cmd, shell=True, stdout=subprocess.PIPE,
                       stderr=subprocess.STDOUT, encoding="utf-8",
                       errors="replace")
    print(p.stdout)
    print("exitcode ", p.returncode)
    return p.returncode, p.stdout


def results_of(prefix):
    """The task results file of the moldft step, parsed, or None."""
    hits = sorted(glob.glob(os.path.join(prefix, "**", "moldft.results.json"),
                            recursive=True))
    if not hits:
        return None
    with open(hits[0]) as fh:
        return json.load(fh)


def iterated(out):
    return "Iteration 0 at" in out


def check(ok, what):
    print(("  ok   " if ok else "  FAIL ") + what)
    return ok


if __name__ == "__main__":
    print("Testing @BINARY@/@TESTCASE@")
    binary = "./@BINARY@"
    ok = True

    # one coarse rung keeps each run to a few seconds; every case below turns on
    # a decision, not on precision
    def madqc(prefix, dft="", molecule="source_name=he"):
        return run(binary + ' --molecule="' + molecule + '" --prefix=' + prefix +
                   ' --dft="protocol=[1.e-4]; ' + dft + '"')

    # --- a plain rerun reuses -------------------------------------------------
    p = "mad_restart822_rerun"
    cleanup(p)
    rc, _ = madqc(p)
    ok &= check(rc == 0, "rerun: first run exits cleanly")
    first = results_of(p)
    ok &= check(first is not None and first.get("archive_id", "0" * 16) != "0" * 16,
                "rerun: the results file records the archive it came from")
    rc, out = madqc(p)
    ok &= check(rc == 0 and "reusing" in out and not iterated(out),
                "rerun: the second run reuses the results without iterating")

    # --- adding a property recomputes it, not the SCF -------------------------
    rc, out = madqc(p, "dipole=true")
    res = results_of(p)
    ok &= check(rc == 0 and not iterated(out), "dipole: no SCF iterations")
    ok &= check(res is not None and "dipole" in res["results"]["properties"],
                "dipole: the dipole is in the results")

    # --- restart none never reuses --------------------------------------------
    p = "mad_restart822_none"
    cleanup(p)
    madqc(p, "restart=none")
    rc, out = madqc(p, "restart=none")
    ok &= check(rc == 0 and iterated(out) and "reusing" not in out,
                "restart none: the second run iterates from the guess")

    # --- a virtual the archive lacks ------------------------------------------
    p = "mad_restart822_nvalpha"
    cleanup(p)
    madqc(p)
    rc, out = madqc(p, "nvalpha=1")
    res = results_of(p)
    ok &= check(rc == 0 and iterated(out), "nvalpha: the second run iterates")
    ok &= check(res is not None and
                len(res["results"]["scf"]["scf_eigenvalues_a"]["vals"]) == 2,
                "nvalpha: the results hold the occupied and the virtual")

    # --- a different eprec is a different Hamiltonian -------------------------
    p = "mad_restart822_eprec"
    cleanup(p)
    madqc(p, "", "source_name=he; eprec=1.e-4")
    e4 = results_of(p)["results"]["properties"]["energy"]
    rc, out = madqc(p, "", "source_name=he; eprec=1.e-2")
    e2 = results_of(p)["results"]["properties"]["energy"]
    ok &= check(rc == 0 and iterated(out), "eprec: the second run re-converges")
    # eprec 1e-2 raises the He energy by ~1.1e-3 Ha over 1e-4
    ok &= check(abs(e2 - e4) > 1.e-4,
                "eprec: the energy is the 1e-2 one (delta %.3e Ha)" % abs(e2 - e4))
    ok &= check("eprec  1.0000e-04 # defined" not in out,
                "eprec: the requested eprec is not rolled back")

    print("final success: ", bool(ok))
    sys.exit(0 if ok else 1)
