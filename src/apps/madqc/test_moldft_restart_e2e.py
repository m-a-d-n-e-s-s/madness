#!/usr/bin/env python3
"""End-to-end restart / checkpoint-robustness test for the madqc SCF driver.

Covers the three Phase-0 driver-hardening fixes (raman thread brief):

  * defect 3 — NextAction::Restart must actually restart. A coarse converged
    run leaves <prefix>.restartdata + calc_info.json; a rerun at a FINER
    protocol must reload those MOs from restartdata (valid() -> Restart), not
    recompute from an atomic guess.

  * defect 1 — a truncated/corrupt calc_info.json must degrade to Redo, not
    throw past the recovery catch and brick the restart. We truncate the
    checkpoint and assert the rerun still exits cleanly (recomputes).

This test runs madqc several times in the same directory, so it is heavier
than a single-shot scripted test. Run it directly on an allocation:

    ./test_moldft_restart_e2e.py            # after CMake substitution
"""

import os
import sys
import subprocess

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


if __name__ == "__main__":
    print("Testing @BINARY@/@TESTCASE@")

    prefix = "mad_@BINARY@_@TESTCASE@"
    ckpt = prefix + ".calc_info.json"
    cleanup(prefix)
    for f in (ckpt, prefix + ".restartdata.00000"):
        if os.path.exists(f):
            os.remove(f)

    binary = "./@BINARY@"
    mol = ' --molecule="source_name=he; eprec=1.e-6" --prefix=' + prefix
    ok = True

    # --- Run 1: coarse protocol, converge + save restartdata --------------
    rc1, _ = run(binary + mol +
                 ' --dft="protocol=[1.e-4]; maxiter=20; save=1;"')
    ok = ok and (rc1 == 0)
    have_restart = os.path.exists(prefix + ".restartdata.00000")
    print("restartdata written after run 1:", have_restart)
    ok = ok and have_restart

    # --- Run 2: finer protocol -> valid() returns Restart -----------------
    # SCF must load MOs from restartdata (SCF.cc prints "... from restartdata").
    rc2, out2 = run(binary + mol +
                    ' --dft="protocol=[1.e-4,1.e-6]; maxiter=20; save=1;"')
    restart_honored = ("from restartdata" in out2) or \
                      ("loading MOs from restartdata archive" in out2)
    print("restart honored on run 2:", restart_honored)
    ok = ok and (rc2 == 0) and restart_honored

    # --- Run 3: corrupt checkpoint must degrade to Redo, not crash --------
    with open(ckpt, "r+") as fh:
        data = fh.read()
        fh.seek(0)
        fh.write(data[: max(1, len(data) // 2)])  # truncate mid-JSON
        fh.truncate()
    print("truncated checkpoint to simulate a killed mid-write")
    rc3, out3 = run(binary + mol +
                    ' --dft="protocol=[1.e-4,1.e-6]; maxiter=20; save=1;"')
    corrupt_recovered = (rc3 == 0)
    print("recovered from corrupt checkpoint:", corrupt_recovered)
    ok = ok and corrupt_recovered

    print("final success: ", ok)
    sys.exit(0 if ok else 1)
