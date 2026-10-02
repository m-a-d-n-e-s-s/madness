#!/usr/bin/env python3
"""Assert the task-record envelope on a madqc calc_info.json.

Usage: check_task_envelope.py <prefix>.calc_info.json [--expect-response]
Exit 0 iff every assertion holds; prints one line per check.
"""
import json, sys

USAGE = "Usage: check_task_envelope.py <prefix>.calc_info.json [--expect-response]"

REQUIRED_TOP = ("schema_name", "schema_version", "provenance")
REQUIRED_PROVENANCE = ("madness", "workflow", "hostname", "nproc", "threads")
# Always present, even on a reload-only record that never iterated.
REQUIRED_SCF_ALWAYS = ("xc", "precision", "scf_total_energy")
# Only present when the SCF actually iterated; absent -> [SKIP], not [FAIL].
REQUIRED_SCF_ENERGY = ("scf_iterations", "nuclear_repulsion_energy",
                       "scf_one_electron_energy", "scf_two_electron_energy", "scf_kinetic_energy",
                       "scf_nuclear_attraction_energy", "scf_coulomb_energy")
REQUIRED_PRECISION = ("k", "thresh", "protocol", "econv", "dconv", "L", "ncoeff")

def main() -> int:
    if len(sys.argv) < 2:
        print(USAGE)
        return 2
    path = sys.argv[1]
    expect_response = "--expect-response" in sys.argv[2:]
    with open(path) as fh:
        ci = json.load(fh)
    fails = 0
    def check(cond, label):
        nonlocal fails
        print(("  [PASS]  " if cond else "  [FAIL]  ") + label)
        if not cond: fails += 1
    if not isinstance(ci, dict):
        check(False, "top-level JSON is an object")
        print(f"{fails} failure(s)")
        return 1
    for k in REQUIRED_TOP: check(k in ci, f"top-level '{k}' present")
    check("tasks" in ci and isinstance(ci["tasks"], list), "top-level 'tasks' is a list")
    tasks = ci.get("tasks") if isinstance(ci.get("tasks"), list) else []
    for k in REQUIRED_PROVENANCE: check(k in ci.get("provenance", {}), f"provenance.{k} present")
    check("git_commit" in ci.get("provenance", {}).get("madness", {}), "provenance.madness.git_commit present")
    scf = next((t for t in tasks if isinstance(t, dict) and t.get("type") in ("scf", "nemo")), None)
    check(scf is not None, "an SCF task entry has type scf|nemo")
    if scf is not None:
        scf_block = scf.get("scf")
        conv_block = scf.get("convergence")
        scf_ok = isinstance(scf_block, dict)
        check(scf_ok, "scf block present")
        reload_only = False
        if scf_ok:
            for k in REQUIRED_SCF_ALWAYS: check(k in scf_block, f"scf.{k} present")
            for k in REQUIRED_PRECISION: check(k in scf_block.get("precision", {}), f"scf.precision.{k} present")
            check(scf.get("precision") == scf_block.get("precision"), "precision mirrored at task top level")
            # A reload-only run (restart read_only / NextAction::ReloadOnly) never
            # calls e_data.add_data(), so the energy decomposition and iteration
            # count are never filled in (A1). That is a valid record, not a
            # broken one -- skip the checks that only make sense after a solve
            # instead of failing them.
            reload_only = "scf_iterations" not in scf_block
            if reload_only:
                print("  [SKIP]  energies (SCF did not iterate: reload-only record)")
            else:
                for k in REQUIRED_SCF_ENERGY: check(k in scf_block, f"scf.{k} present")
                need = ("scf_one_electron_energy", "scf_two_electron_energy", "nuclear_repulsion_energy", "scf_total_energy")
                if all(k in scf_block for k in need):
                    total = (scf_block["scf_one_electron_energy"] + scf_block["scf_two_electron_energy"]
                             + scf_block["nuclear_repulsion_energy"]
                             + scf_block.get("scf_pcm_energy", 0.0) + scf_block.get("scf_dispersion_correction_energy", 0.0))
                    check(abs(total - scf_block["scf_total_energy"]) < 1e-8,
                          f"energy components sum to scf_total_energy (diff {total - scf_block['scf_total_energy']:.2e})")
                else:
                    check(False, "energy components present for the sum check")
        conv_ok = isinstance(conv_block, dict)
        check(conv_ok, "convergence block present")
        if conv_ok:
            check(conv_block.get("status") in ("converged", "unconverged"), "convergence.status set")
            # Same reload-only exemption as above: conv_res.iterations stays at
            # its -1 default when the SCF never iterated (A1).
            if not reload_only:
                check(conv_block.get("iterations", -1) >= 1, "convergence.iterations >= 1")
        check("wall_s" in scf.get("provenance", {}), "task provenance.wall_s present")
    if expect_response:
        resp = next((t for t in tasks if isinstance(t, dict) and t.get("type") == "response"), None)
        check(resp is not None, "a response task entry exists")
        if resp is not None:
            check("wall_s" in resp.get("provenance", {}), "response task provenance.wall_s present")
    print(f"{fails} failure(s)")
    return 1 if fails else 0

if __name__ == "__main__":
    sys.exit(main())
