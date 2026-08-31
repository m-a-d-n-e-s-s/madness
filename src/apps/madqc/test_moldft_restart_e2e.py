#!/usr/bin/env python3
"""End-to-end restart / checkpoint-robustness test for the madqc SCF driver.

Three ways a rerun in a directory that already holds results can go wrong, all
decided in SCFApplication::run (chem/Applications.hpp) BEFORE the engine is
built -- so none of them is reachable from a single-shot test:

  * restart honored -- a coarse converged run leaves an orbital archive and a
    checkpoint; a rerun at a FINER protocol must reload those MOs
    (valid() -> Restart) rather than start over from an atomic guess.

  * a changed Hamiltonian must NOT be answered out of the checkpoint. valid()
    weighs thresholds, requested properties and the archive, every one of which
    a changed `xc` leaves intact, so the cached results pass all of it. Without
    the guard, `xc=lda` in a directory holding a converged `xc=hf` run reported
    the HF energy as the LDA answer and never constructed an SCF. This is the
    one case where being wrong is silent, which is why the assertion here is on
    the energy and not just on an exit code.

  * a truncated/corrupt checkpoint must degrade to Redo, not throw past the
    recovery catch and brick the restart.

NB on paths -- getting these wrong is why this test sat disabled and unrun:
madqc executes each task in <prefix>/task_<n>/<engine>/, so the checkpoint
SCFApplication reads back is <prefix>/task_0/moldft/moldft.calc_info.json and
the orbital archive is <prefix>/task_0/moldft/<prefix>.restartdata.00000.
The <prefix>.calc_info.json sitting in the cwd is the aggregated report, which
nothing ever reads back -- truncating that one tests nothing.

`lda` is used as the second functional deliberately: MADNESS falls back to
xcfunctional_ldaonly.cc when libxc is absent, so this case works in both
configurations.

Runs madqc four times in the same directory, so it is heavier than a
single-shot scripted test.

    ./test_moldft_restart_e2e.py            # after CMake substitution
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


def find_one(*parts):
    """The single file matching <prefix>/**/<name>, or None."""
    hits = sorted(glob.glob(os.path.join(*parts), recursive=True))
    return hits[0] if hits else None


def energy_from(ckpt):
    """The SCF energy recorded in a task checkpoint, or None."""
    try:
        with open(ckpt) as fh:
            return json.load(fh)["properties"]["energy"]
    except Exception as e:                      # noqa: BLE001
        print("could not read an energy from", ckpt, ":", e)
        return None


def check(ok, what):
    print(("  ok   " if ok else "  FAIL ") + what)
    return ok


if __name__ == "__main__":
    print("Testing @BINARY@/@TESTCASE@")

    prefix = "mad_@BINARY@_@TESTCASE@"
    cleanup(prefix)

    binary = "./@BINARY@"
    mol = ' --molecule="source_name=he; eprec=1.e-6" --prefix=' + prefix
    coarse = ' --dft="protocol=[1.e-4]; maxiter=20; save=1;"'
    fine = ' --dft="protocol=[1.e-4,1.e-6]; maxiter=20; save=1;"'
    # identical to `fine` except for the functional, so a recompute here can
    # only be attributed to the change of Hamiltonian
    fine_lda = ' --dft="protocol=[1.e-4,1.e-6]; maxiter=20; save=1; xc=lda;"'

    ok = True

    # --- Run 1: coarse protocol, converge + save the orbital archive --------
    rc, _ = run(binary + mol + coarse)
    ok &= check(rc == 0, "run 1 (coarse) exits cleanly")
    archive = find_one(prefix, "**", prefix + ".restartdata.00000")
    ok &= check(archive is not None, "run 1 wrote an orbital archive")
    print("   archive:", archive)

    ckpt = find_one(prefix, "**", "moldft.calc_info.json")
    ok &= check(ckpt is not None, "run 1 wrote a task checkpoint")
    print("   checkpoint:", ckpt)
    if not ok:
        print("final success: ", False)
        sys.exit(1)

    # --- Run 2: finer protocol -> valid() returns Restart -------------------
    rc, out = run(binary + mol + fine)
    ok &= check(rc == 0, "run 2 (finer) exits cleanly")
    ok &= check("from restartdata" in out,
                "run 2 reloaded the MOs instead of starting from a guess")
    e_hf = energy_from(ckpt)
    ok &= check(e_hf is not None, "run 2 recorded an energy")
    print("   hf energy:", e_hf)

    # --- Run 3: same everything, different functional -> must recompute -----
    rc, out = run(binary + mol + fine_lda)
    ok &= check(rc == 0, "run 3 (xc=lda) exits cleanly")
    ok &= check("different Hamiltonian" in out,
                "run 3 refused the checkpoint written for another functional")
    e_lda = energy_from(ckpt)
    ok &= check(e_lda is not None, "run 3 recorded an energy")
    print("   lda energy:", e_lda)
    # HF and LDA on helium differ by ~2.7e-2 Ha. Any tolerance far above the
    # convergence threshold and far below that gap distinguishes "recomputed"
    # from "handed back the cached HF number", which is the whole point.
    if e_hf is not None and e_lda is not None:
        ok &= check(abs(e_lda - e_hf) > 1.e-3,
                    "run 3 returned the LDA energy, not the cached HF one "
                    "(delta %.3e Ha)" % abs(e_lda - e_hf))

    # --- Run 4: corrupt checkpoint must degrade to Redo, not crash ----------
    with open(ckpt, "r+") as fh:
        data = fh.read()
        fh.seek(0)
        fh.write(data[: max(1, len(data) // 2)])   # truncate mid-JSON
        fh.truncate()
    print("truncated", ckpt, "to simulate a process killed mid-write")
    rc, _ = run(binary + mol + fine_lda)
    ok &= check(rc == 0, "run 4 recovered from a corrupt checkpoint")

    print("final success: ", bool(ok))
    sys.exit(0 if ok else 1)
