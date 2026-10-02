#!/usr/bin/env python3
"""`--optimize` and `--wf=scf` in one directory must not answer for each other (#822).

Both run in <prefix>/task_0/moldft/ and share its orbital archive. What keeps
that safe is the archive_id: every save gets a new one, and each results file
records the id it was computed from.

  1. --wf=scf              the energy at the input geometry
  2. --optimize --wf=scf   optimizes; writes optimize.results.json
  3. --optimize --wf=scf   already finished: no SCF iterations, same energy
  4. --wf=scf              the archive now holds the optimized geometry, so the
                           SCF results of run 1 must NOT be reused: the energy is
                           recomputed at the input geometry. The old checkpoint
                           logic reported run 1's energy while holding the
                           optimized orbitals.

LiH, the scf_lih_optimize deck; four madqc runs, so this is a verylong test.

    ./test_madqc_optimize_results.py        # after CMake substitution
"""

import glob
import json
import os
import subprocess
import sys

sys.path.append("@CMAKE_SOURCE_DIR@/bin")
from test_utilities import cleanup  # noqa: E402

DECK = """
optimization
    maxiter     10
end

dft
    xc          hf
    protocol    [1e-4, 1e-6]
    dconv       1e-4
end

molecule
    units   atomic
    eprec   1e-6
    Li  0.0  0.0  0.0
    H   0.0  0.0  3.05
end
"""


def run(cmd):
    print("executing\n ", cmd)
    p = subprocess.run(cmd, shell=True, stdout=subprocess.PIPE,
                       stderr=subprocess.STDOUT, encoding="utf-8",
                       errors="replace")
    print(p.stdout)
    print("exitcode ", p.returncode)
    return p.returncode, p.stdout


def load(prefix, name):
    hits = sorted(glob.glob(os.path.join(prefix, "**", name), recursive=True))
    if not hits:
        return None
    with open(hits[0]) as fh:
        return json.load(fh)


def check(ok, what):
    print(("  ok   " if ok else "  FAIL ") + what)
    return ok


if __name__ == "__main__":
    print("Testing @BINARY@/@TESTCASE@")
    prefix = "mad_optimize822"
    cleanup(prefix)
    deck = prefix + ".in"
    with open(deck, "w") as fh:
        fh.write(DECK)
    binary = "./@BINARY@"
    ok = True

    rc, _ = run(binary + " --wf=scf " + deck)
    ok &= check(rc == 0, "run 1 (scf) exits cleanly")
    e_input = load(prefix, "moldft.results.json")["results"]["properties"]["energy"]

    rc, _ = run(binary + " --optimize --wf=scf " + deck)
    opt = load(prefix, "optimize.results.json")
    ok &= check(rc == 0 and opt is not None, "run 2 (optimize) writes optimize.results.json")
    e_opt = opt["optimization_results"]["final_energy"]

    rc, out = run(binary + " --optimize --wf=scf " + deck)
    again = load(prefix, "optimize.results.json")["optimization_results"]["final_energy"]
    ok &= check(rc == 0 and "optimization already finished" in out and
                "Iteration 0 at" not in out,
                "run 3 (optimize again) is recognized as finished, no SCF")
    ok &= check(abs(again - e_opt) < 1.e-6, "run 3 reports the optimized energy")

    rc, out = run(binary + " --wf=scf " + deck)
    e_after = load(prefix, "moldft.results.json")["results"]["properties"]["energy"]
    ok &= check(rc == 0 and "reusing" not in out,
                "run 4 (scf) does not reuse results the optimizer made stale")
    # the optimized energy lies ~9e-6 Ha below the input-geometry one; the SCF is
    # converged far tighter than that at protocol 1e-6
    ok &= check(abs(e_after - e_input) < 1.e-6 and abs(e_after - e_opt) > 3.e-6,
                "run 4 reports the input-geometry energy (%.9f vs %.9f, optimized %.9f)"
                % (e_after, e_input, e_opt))

    os.remove(deck)
    print("final success: ", bool(ok))
    sys.exit(0 if ok else 1)
