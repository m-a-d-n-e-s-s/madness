#!/usr/bin/env python3
"""Comprehensive 4-Tier End-to-End (E2E) Test Suite for MADNESS madqc.

Validates:
  - Calculation restarting and checkpoint integrity across all-electron and GTH pseudopotentials
  - Hamiltonian mismatch detection and stale energy rejection
  - Potential construction integrity (pure_ae / psp_calc parameter precedence)
  - Hellmann-Feynman nuclear force calculation vs external NWChem cc-pVQZ reference gradients
  - Analytical GTH pseudopotential derivatives for pure-PSP and mixed AE/PSP systems
  - Numerical gradient verification vs Central Finite Differences (dE/dR)
  - External validation vs NWChem 7.3.1 pspw plane-wave gradients (task pspw gradient)
  - Translational invariance (vanishing net molecular force) and physical restoring forces

Methodology:
  Tier 1: Feature Coverage (Core Capabilities, AE and PSP single-point forces)
  Tier 2: Boundary & Corner Cases (Hamiltonian Rejection, Finite Differences, NWChem Oracles)
  Tier 3: Cross-Feature Combinations (Mixed AE/PSP, Translational Invariance)
  Tier 4: Real-World Scenarios (Multi-Atom Checkpointing, Physical Restoring Forces)
"""

import argparse
import inspect
import json
import os
import re
import shutil
import subprocess
import sys
import tempfile
from typing import Any, Dict, List, Optional, Tuple

DEFAULT_MADQC = "@BINARY@"
if DEFAULT_MADQC.startswith("@") or not os.path.isabs(DEFAULT_MADQC):
    candidate = os.path.abspath(DEFAULT_MADQC) if not DEFAULT_MADQC.startswith("@") else None
    if candidate and os.path.isfile(candidate):
        DEFAULT_MADQC = candidate
    elif os.path.isfile("/tmp/madrel-build/src/apps/madqc/madqc"):
        DEFAULT_MADQC = "/tmp/madrel-build/src/apps/madqc/madqc"
    elif shutil.which("madqc"):
        DEFAULT_MADQC = shutil.which("madqc")
DEFAULT_NWCHEM = "/home/rjh/Devel/nwchem-7.3.1/bin/LINUX64/nwchem"
DEFAULT_NWCHEM_BASIS = "/home/rjh/Devel/nwchem-7.3.1/src/basis/libraries/"

# Reference NWChem 7.3.1 DFT Gradients (Ha/Bohr)
# All-electron evaluated with: xc slater vwn_5, basis spherical cc-pvqz, noautoz nocenter noautosym
# PSPW evaluated with: xc lda, SC 20.0, ngrid 128 128 128, aperiodic
NWCHEM_REF_FORCES = {
    "h2_lda_1.4bohr": {
        "energy": -1.13726435,
        "gradients": [
            [0.0, 0.0, 0.016597],   # H1 at z = -0.7
            [0.0, 0.0, -0.016597],  # H2 at z = +0.7
        ],
    },
    "h2_gth_lda_1.2bohr": {
        "energy": -1.12938613,
        "gradients": [
            [0.0, 0.0, 0.126737],   # H1 at z = -0.6
            [0.0, 0.0, -0.126737],  # H2 at z = +0.6
        ],
    },
    "h2o_lda_equilibrium": {
        "energy": -75.90704166,
        "gradients": [
            [0.0, 0.0, -0.034823],      # O  at (0.0, 0.0, 0.710024)
            [0.022009, 0.0, 0.017412],   # H1 at (-1.426508, 0.0, -0.355012)
            [-0.022009, 0.0, 0.017412],  # H2 at (1.426508, 0.0, -0.355012)
        ],
    },
    "water_mixed_lda_equilibrium": {
        "energy": -17.18646296,
        "gradients": [
            [0.0, 0.0, -0.034898],      # psO at (0.0, 0.0, 0.710024)
            [0.022884, 0.0, 0.017445],   # H1  at (-1.426508, 0.0, -0.355012)
            [-0.022884, 0.0, 0.017445],  # H2  at (1.426508, 0.0, -0.355012)
        ],
    },
}

# Reference Ground-State Energies (Ha)
REF_ENERGIES = {
    "he_lda_ae": -2.83445,
    "be_lda_ae": -14.4468,
    "be_lda_gth": -0.993398,
    "be_lda_gth_k6": -0.991245,
    "be_pbe_gth": -14.5354,
    "methane_lda_gth": -8.03824,
    "water_lda_mixed": -17.186463,
    "ethylene_allpsp_lda": -13.695742,
    "h2_gth_lda_1.2": -1.123958,
}


def get_environment() -> Dict[str, str]:
    """Return environment variables enforcing strict thread limits."""
    env = os.environ.copy()
    env["MAD_NUM_THREADS"] = "4"
    env["OMP_NUM_THREADS"] = "1"
    env["MKL_NUM_THREADS"] = "1"
    if "NWCHEM_BASIS_LIBRARY" not in env and os.path.exists(DEFAULT_NWCHEM_BASIS):
        env["NWCHEM_BASIS_LIBRARY"] = DEFAULT_NWCHEM_BASIS
    return env


def run_command(
    cmd: List[str], cwd: str, timeout: int = 120, verbose: bool = False
) -> Tuple[int, str]:
    """Execute command in specified directory with thread safety."""
    if verbose:
        print(f"\n[EXEC] (cwd={cwd}): {' '.join(cmd)}")
    try:
        p = subprocess.run(
            cmd,
            cwd=cwd,
            env=get_environment(),
            stdout=subprocess.PIPE,
            stderr=subprocess.STDOUT,
            encoding="utf-8",
            errors="replace",
            timeout=timeout,
        )
        if verbose:
            print(p.stdout)
        return p.returncode, p.stdout
    except subprocess.TimeoutExpired:
        return -1, f"Command timed out after {timeout} seconds"
    except Exception as e:
        return -2, f"Execution failed: {e}"


def read_calc_info(dir_path: str, prefix: str) -> Optional[Dict[str, Any]]:
    """Read the aggregated calc_info.json report."""
    candidates = [
        os.path.join(dir_path, f"{prefix}.calc_info.json"),
        os.path.join(dir_path, f"{prefix.lower()}.calc_info.json"),
        os.path.join(dir_path, "moldft.calc_info.json"),
    ]
    for json_path in candidates:
        if os.path.exists(json_path):
            try:
                with open(json_path, "r", encoding="utf-8") as f:
                    return json.load(f)
            except Exception as e:
                print(f"[WARN] Failed to parse {json_path}: {e}")
    return None


def get_task_energy(info: Optional[Dict[str, Any]], task_idx: int = 0) -> Optional[float]:
    """Extract converged total SCF energy from calc_info JSON."""
    if not info or "tasks" not in info or len(info["tasks"]) <= task_idx:
        return None
    task = info["tasks"][task_idx]
    # Check properties.energy first
    if "properties" in task and "energy" in task["properties"]:
        return float(task["properties"]["energy"])
    # Fallback to scf_total_energy
    if "scf" in task and "scf_total_energy" in task["scf"]:
        return float(task["scf"]["scf_total_energy"])
    return None


def get_task_gradient(
    info: Optional[Dict[str, Any]], task_idx: int = 0
) -> Optional[List[List[float]]]:
    """Extract Cartesian nuclear gradients tensor (natoms x 3) from calc_info JSON."""
    if not info or "tasks" not in info or len(info["tasks"]) <= task_idx:
        return None
    task = info["tasks"][task_idx]
    grad_obj = task.get("properties", {}).get("gradient") or task.get("scf", {}).get(
        "properties", {}
    ).get("gradient")
    if not grad_obj or "vals" not in grad_obj:
        return None
    vals = grad_obj["vals"]
    if len(vals) % 3 != 0:
        return None
    natoms = len(vals) // 3
    return [[vals[i * 3 + j] for j in range(3)] for i in range(natoms)]


def parse_stdout_energies(stdout: str) -> Dict[str, float]:
    """Parse energy components from madqc stdout table."""
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
        matches = re.findall(pat, stdout)
        if matches:
            try:
                data[key] = float(matches[-1])
            except ValueError:
                pass
    return data


def parse_nwchem_ion_forces(stdout: str) -> Optional[List[List[float]]]:
    """Parse Ion Forces: table from NWChem pspw gradient stdout.

    NWChem reports Hellmann-Feynman ion forces F = -dE/dR under 'Ion Forces:'.
    This function converts forces to Cartesian energy gradients dE/dR = -F (Ha/Bohr).
    """
    idx = stdout.find("Ion Forces:")
    if idx == -1:
        return None
    lines = stdout[idx:].splitlines()
    gradients = []
    pattern = re.compile(
        r"^\s*\d+\s+[A-Za-z]+\s*\(\s*([-0-9.Ee+]+)\s+([-0-9.Ee+]+)\s+([-0-9.Ee+]+)\s*\)"
    )
    for line in lines[1:]:
        if "C.O.M." in line or "=====" in line:
            if gradients:
                break
        m = pattern.match(line)
        if m:
            fx, fy, fz = float(m.group(1)), float(m.group(2)), float(m.group(3))
            gradients.append([-fx, -fy, -fz])
    return gradients if gradients else None


class TestCaseResult:
    def __init__(self, test_id: str, tier: int, description: str):
        self.test_id = test_id
        self.tier = tier
        self.description = description
        self.passed = False
        self.skipped = False
        self.error_category: Optional[str] = None
        self.details: str = ""
        self.observed: Dict[str, Any] = {}
        self.expected: Dict[str, Any] = {}


# ============================================================================
# TIER 1: Feature Coverage (Core Capabilities)
# ============================================================================


def test_t1_1_restart_ae_clean(madqc: str, verbose: bool) -> TestCaseResult:
    """T1.1: Small all-electron calculation restart (He, LDA)."""
    res = TestCaseResult("T1.1", 1, "Small AE calculation restart (He, LDA)")
    with tempfile.TemporaryDirectory(prefix="madqc_e2e_t1_1_") as tmpdir:
        # Continuous baseline
        in_cont = os.path.join(tmpdir, "he_cont.in")
        with open(in_cont, "w") as f:
            f.write(
                """dft
  xc lda
  protocol [1e-3, 1e-5]
  dconv 1e-5
  k 6
  maxiter 25
  prefix he_cont
  save 1
end

geometry
  units bohr
  no_orient true
  He 0.0 0.0 0.0
end
"""
            )
        rc_c, out_c = run_command([madqc, in_cont], tmpdir, verbose=verbose)
        if rc_c != 0:
            res.error_category = "BASELINE_FAILURE"
            res.details = f"Continuous baseline failed with exit code {rc_c}"
            return res
        info_c = read_calc_info(tmpdir, "he_cont")
        e_cont = get_task_energy(info_c)
        if e_cont is None:
            res.error_category = "BASELINE_FAILURE"
            res.details = "Continuous baseline did not produce valid energy"
            return res

        # Stage 1: Coarse run
        in_s1 = os.path.join(tmpdir, "he_s1.in")
        with open(in_s1, "w") as f:
            f.write(
                """dft
  xc lda
  protocol [1e-3]
  dconv 1e-4
  k 6
  maxiter 25
  prefix he_stage
  save 1
end

geometry
  units bohr
  no_orient true
  He 0.0 0.0 0.0
end
"""
            )
        rc1, out1 = run_command([madqc, in_s1], tmpdir, verbose=verbose)
        if rc1 != 0:
            res.error_category = "STAGE1_FAILURE"
            res.details = f"Stage 1 coarse run failed with exit code {rc1}"
            return res

        # Stage 2: Finer restart
        in_s2 = os.path.join(tmpdir, "he_s2.in")
        with open(in_s2, "w") as f:
            f.write(
                """dft
  xc lda
  protocol [1e-3, 1e-5]
  dconv 1e-5
  k 6
  maxiter 25
  prefix he_stage
  save 1
end

geometry
  units bohr
  no_orient true
  He 0.0 0.0 0.0
end
"""
            )
        rc2, out2 = run_command([madqc, in_s2], tmpdir, verbose=verbose)
        if rc2 != 0:
            res.error_category = "STAGE2_FAILURE"
            res.details = f"Stage 2 restart run failed with exit code {rc2}"
            return res

        info_s2 = read_calc_info(tmpdir, "he_stage")
        e_restarted = get_task_energy(info_s2)
        if e_restarted is None:
            res.error_category = "STAGE2_FAILURE"
            res.details = "Stage 2 restart did not produce valid energy"
            return res

        res.expected["energy_continuous"] = e_cont
        res.observed["energy_restarted"] = e_restarted
        delta = abs(e_restarted - e_cont)
        res.observed["delta"] = delta

        # Must reload from restartdata
        reloaded = (
            "reading initial orbitals from restartdata" in out2
            or "from restartdata" in out2
        )
        if not reloaded:
            res.error_category = "RELOAD_NOT_TRIGGERED"
            res.details = "Did not reload MOs from restartdata"
            return res

        if delta > 1e-7:
            res.error_category = "ENERGY_MISMATCH"
            res.details = (
                f"Restarted energy {e_restarted:.8f} differs from continuous "
                f"{e_cont:.8f} by {delta:.2e} Ha"
            )
            return res

        res.passed = True
        res.details = f"Clean AE restart converged, delta={delta:.2e} Ha (< 1e-7 Ha)"
        return res


def test_t1_2_restart_psp_clean(madqc: str, verbose: bool) -> TestCaseResult:
    """T1.2: Small GTH pseudopotential calculation restart (Be, LDA)."""
    res = TestCaseResult("T1.2", 1, "Small GTH PSP restart (Be, LDA)")
    with tempfile.TemporaryDirectory(prefix="madqc_e2e_t1_2_") as tmpdir:
        # Continuous baseline
        in_cont = os.path.join(tmpdir, "be_psp_cont.in")
        with open(in_cont, "w") as f:
            f.write(
                """dft
  xc lda
  protocol [1e-3, 1e-5]
  dconv 1e-5
  k 6
  maxiter 25
  psp_calc true
  prefix be_psp_cont
  save 1
end

geometry
  units bohr
  no_orient true
  Be 0.0 0.0 0.0
end
"""
            )
        rc_c, out_c = run_command([madqc, in_cont], tmpdir, verbose=verbose)
        if rc_c != 0:
            res.error_category = "BASELINE_FAILURE"
            res.details = f"Continuous PSP baseline failed with exit code {rc_c}"
            return res
        info_c = read_calc_info(tmpdir, "be_psp_cont")
        e_cont = get_task_energy(info_c)
        if e_cont is None:
            res.error_category = "BASELINE_FAILURE"
            res.details = "Continuous PSP baseline did not produce valid energy"
            return res

        # Stage 1: Coarse run
        in_s1 = os.path.join(tmpdir, "be_psp_s1.in")
        with open(in_s1, "w") as f:
            f.write(
                """dft
  xc lda
  protocol [1e-3]
  dconv 1e-4
  k 6
  maxiter 25
  psp_calc true
  prefix be_psp_stage
  save 1
end

geometry
  units bohr
  no_orient true
  Be 0.0 0.0 0.0
end
"""
            )
        rc1, out1 = run_command([madqc, in_s1], tmpdir, verbose=verbose)
        if rc1 != 0:
            res.error_category = "STAGE1_FAILURE"
            res.details = f"Stage 1 PSP coarse run failed with exit code {rc1}"
            return res

        # Stage 2: Finer restart
        in_s2 = os.path.join(tmpdir, "be_psp_s2.in")
        with open(in_s2, "w") as f:
            f.write(
                """dft
  xc lda
  protocol [1e-3, 1e-5]
  dconv 1e-5
  k 6
  maxiter 25
  psp_calc true
  prefix be_psp_stage
  save 1
end

geometry
  units bohr
  no_orient true
  Be 0.0 0.0 0.0
end
"""
            )
        rc2, out2 = run_command([madqc, in_s2], tmpdir, verbose=verbose)
        if rc2 != 0:
            res.error_category = "STAGE2_FAILURE"
            res.details = f"Stage 2 PSP restart run failed with exit code {rc2}"
            return res

        info_s2 = read_calc_info(tmpdir, "be_psp_stage")
        e_restarted = get_task_energy(info_s2)
        if e_restarted is None:
            res.error_category = "STAGE2_FAILURE"
            res.details = "Stage 2 PSP restart did not produce valid energy"
            return res

        res.expected["energy_continuous"] = e_cont
        res.observed["energy_restarted"] = e_restarted
        delta = abs(e_restarted - e_cont)
        res.observed["delta"] = delta

        reloaded = (
            "reading initial orbitals from restartdata" in out2
            or "from restartdata" in out2
        )
        if not reloaded:
            res.error_category = "RELOAD_NOT_TRIGGERED"
            res.details = "Did not reload PSP MOs from restartdata"
            return res

        if delta > 1e-7:
            res.error_category = "ENERGY_MISMATCH"
            res.details = (
                f"Restarted PSP energy {e_restarted:.8f} differs from continuous "
                f"{e_cont:.8f} by {delta:.2e} Ha"
            )
            return res

        res.passed = True
        res.details = f"Clean PSP restart converged, delta={delta:.2e} Ha (< 1e-7 Ha)"
        return res


def test_t1_3_protocol_reload(madqc: str, verbose: bool) -> TestCaseResult:
    """T1.3: Multi-stage coarse-to-fine protocol restart reload (H2O, AE)."""
    res = TestCaseResult("T1.3", 1, "Multi-stage coarse-to-fine protocol reload (H2O)")
    with tempfile.TemporaryDirectory(prefix="madqc_e2e_t1_3_") as tmpdir:
        in_stage1 = os.path.join(tmpdir, "h2o_s1.in")
        with open(in_stage1, "w") as f:
            f.write(
                """dft
  xc lda
  protocol [1e-3]
  dconv 1e-4
  k 7
  prefix h2o_proto
  save 1
end

geometry
  units bohr
  no_orient true
  O  0.00000000000000  0.00000000000000  0.71002389091700
  H -1.42650808202133  0.00000000000000 -0.35501194545850
  H  1.42650808202133  0.00000000000000 -0.35501194545850
end
"""
            )
        rc1, out1 = run_command([madqc, in_stage1], tmpdir, verbose=verbose)
        if rc1 != 0:
            res.error_category = "STAGE1_FAILURE"
            res.details = f"Stage 1 coarse H2O run failed with exit code {rc1}"
            return res

        in_stage2 = os.path.join(tmpdir, "h2o_s2.in")
        with open(in_stage2, "w") as f:
            f.write(
                """dft
  xc lda
  protocol [1e-3, 1e-5]
  dconv 1e-5
  k 7
  prefix h2o_proto
  save 1
end

geometry
  units bohr
  no_orient true
  O  0.00000000000000  0.00000000000000  0.71002389091700
  H -1.42650808202133  0.00000000000000 -0.35501194545850
  H  1.42650808202133  0.00000000000000 -0.35501194545850
end
"""
            )
        rc2, out2 = run_command([madqc, in_stage2], tmpdir, verbose=verbose)
        if rc2 != 0:
            res.error_category = "STAGE2_FAILURE"
            res.details = f"Stage 2 refined H2O run failed with exit code {rc2}"
            return res

        reloaded = (
            "reading initial orbitals from restartdata" in out2
            or "from restartdata" in out2
        )
        if not reloaded:
            res.error_category = "RELOAD_NOT_TRIGGERED"
            res.details = "Stage 2 did not reload from restartdata"
            return res

        info = read_calc_info(tmpdir, "h2o_proto")
        energy = get_task_energy(info)
        res.observed["energy"] = energy
        res.passed = True
        res.details = f"Protocol reload succeeded with energy {energy:.6f} Ha"
        return res


def test_t1_4_fresh_baseline(madqc: str, verbose: bool) -> TestCaseResult:
    """T1.4: Fresh calculation baseline energy checks (He AE, Be GTH)."""
    res = TestCaseResult("T1.4", 1, "Fresh calculation baseline energies")
    with tempfile.TemporaryDirectory(prefix="madqc_e2e_t1_4_") as tmpdir:
        # He AE
        in_he = os.path.join(tmpdir, "he_base.in")
        with open(in_he, "w") as f:
            f.write(
                """dft
  xc lda
  protocol [1e-4]
  dconv 1e-5
  k 6
  prefix he_baseline
end

geometry
  units bohr
  no_orient true
  He 0.0 0.0 0.0
end
"""
            )
        rc_he, _ = run_command([madqc, in_he], tmpdir, verbose=verbose)
        if rc_he != 0:
            res.error_category = "HE_FAILED"
            res.details = "Fresh He calculation failed"
            return res
        e_he = get_task_energy(read_calc_info(tmpdir, "he_baseline"))

        # Be GTH
        in_be = os.path.join(tmpdir, "be_base.in")
        with open(in_be, "w") as f:
            f.write(
                """dft
  xc lda
  protocol [1e-4]
  dconv 1e-5
  k 6
  psp_calc true
  prefix be_baseline
end

geometry
  units bohr
  no_orient true
  Be 0.0 0.0 0.0
end
"""
            )
        rc_be, _ = run_command([madqc, in_be], tmpdir, verbose=verbose)
        if rc_be != 0:
            res.error_category = "BE_FAILED"
            res.details = "Fresh Be GTH calculation failed"
            return res
        e_be = get_task_energy(read_calc_info(tmpdir, "be_baseline"))

        res.observed["he_energy"] = e_he
        res.observed["be_energy"] = e_be
        res.expected["he_ref"] = REF_ENERGIES["he_lda_ae"]
        res.expected["be_ref"] = REF_ENERGIES["be_lda_gth_k6"]

        if e_he is None or abs(e_he - REF_ENERGIES["he_lda_ae"]) > 5e-3:
            res.error_category = "HE_ENERGY_MISMATCH"
            res.details = f"He energy {e_he} differs from ref {REF_ENERGIES['he_lda_ae']}"
            return res

        if e_be is None or abs(e_be - REF_ENERGIES["be_lda_gth_k6"]) > 1e-3:
            res.error_category = "BE_ENERGY_MISMATCH"
            res.details = f"Be energy {e_be} differs from ref {REF_ENERGIES['be_lda_gth_k6']}"
            return res

        res.passed = True
        res.details = f"Baselines verified: He={e_he:.4f} Ha, Be={e_be:.4f} Ha"
        return res


def test_t1_5_forces_ae_single_point(madqc: str, verbose: bool) -> TestCaseResult:
    """T1.5: Single-point Hellmann-Feynman forces for small AE molecule (H2)."""
    res = TestCaseResult("T1.5", 1, "Single-point Hellmann-Feynman forces (H2)")
    with tempfile.TemporaryDirectory(prefix="madqc_e2e_t1_5_") as tmpdir:
        deck = os.path.join(tmpdir, "h2_forces.in")
        with open(deck, "w") as f:
            f.write(
                """dft
  xc lda
  protocol [1e-4]
  dconv 1e-6
  k 8
  derivatives true
  prefix h2_force
end

geometry
  units bohr
  no_orient true
  H 0.0 0.0 -0.7
  H 0.0 0.0  0.7
end
"""
            )
        rc, out = run_command([madqc, deck], tmpdir, verbose=verbose)
        if rc != 0:
            res.error_category = "MADQC_EXEC_ERROR"
            res.details = f"madqc failed with exit code {rc}"
            return res

        info = read_calc_info(tmpdir, "h2_force")
        grad = get_task_gradient(info)
        if not grad or len(grad) != 2:
            res.error_category = "GRADIENT_MISSING"
            res.details = "properties.gradient missing or invalid in calc_info.json"
            return res

        fz1 = grad[0][2]
        fz2 = grad[1][2]
        res.observed["grad_H1"] = grad[0]
        res.observed["grad_H2"] = grad[1]

        # Non-zero forces along z
        if abs(fz1) < 1e-3 or abs(fz2) < 1e-3:
            res.error_category = "ZERO_FORCES"
            res.details = f"Forces along bond axis are vanishingly small: {fz1}, {fz2}"
            return res

        res.passed = True
        res.details = f"Hellmann-Feynman gradient computed: fz(H1)={fz1:.6f}, fz(H2)={fz2:.6f} Ha/Bohr"
        return res


def test_t1_6_forces_mixed_ae_psp_single_point(
    madqc: str, verbose: bool
) -> TestCaseResult:
    """T1.6: Single-point mixed AE/PSP forces (water_lda: psO with AE H, derivatives true)."""
    res = TestCaseResult(
        "T1.6", 1, "Single-point mixed AE/PSP forces (water_lda with psO and AE H)"
    )
    with tempfile.TemporaryDirectory(prefix="madqc_e2e_t1_6_") as tmpdir:
        deck = os.path.join(tmpdir, "water_mixed_forces.in")
        with open(deck, "w") as f:
            f.write(
                """dft
  xc lda
  protocol [1e-4]
  dconv 1e-5
  k 7
  derivatives true
  prefix water_mixed_force
end

geometry
  units bohr
  no_orient true
  psO  0.00000000000000  0.00000000000000  0.71002389091700
  H   -1.42650808202133  0.00000000000000 -0.35501194545850
  H    1.42650808202133  0.00000000000000 -0.35501194545850
end
"""
            )
        rc, out = run_command([madqc, deck], tmpdir, timeout=120, verbose=verbose)
        if rc != 0:
            res.error_category = "MADQC_EXEC_ERROR"
            res.details = f"madqc failed with exit code {rc}"
            return res

        info = read_calc_info(tmpdir, "water_mixed_force")
        grad = get_task_gradient(info)
        if not grad or len(grad) != 3:
            res.error_category = "GRADIENT_MISSING"
            res.details = "properties.gradient missing or invalid in calc_info.json"
            return res

        fz_O = grad[0][2]
        fx_H1 = grad[1][0]
        fz_H1 = grad[1][2]
        fx_H2 = grad[2][0]
        fz_H2 = grad[2][2]

        res.observed["grad_O"] = grad[0]
        res.observed["grad_H1"] = grad[1]
        res.observed["grad_H2"] = grad[2]

        # Verify non-zero forces
        if abs(fz_O) < 1e-3 or abs(fz_H1) < 1e-3:
            res.error_category = "ZERO_FORCES"
            res.details = f"Forces are vanishingly small: fz(O)={fz_O}, fz(H1)={fz_H1}"
            return res

        # Expected sign: dE/dz(O) < 0 and dE/dz(H) > 0
        if fz_O >= 0 or fz_H1 <= 0 or fz_H2 <= 0:
            res.error_category = "UNEXPECTED_FORCE_SIGN"
            res.details = (
                f"Unexpected force signs: fz(O)={fz_O:.6f} (expected < 0), "
                f"fz(H)={fz_H1:.6f} (expected > 0)"
            )
            return res

        res.passed = True
        res.details = (
            f"Mixed AE/PSP forces computed: fz(O)={fz_O:.6f}, "
            f"fx(H1)={fx_H1:.6f}, fz(H1)={fz_H1:.6f} Ha/Bohr"
        )
        return res


def test_t1_7_forces_pure_psp_single_point(
    madqc: str, verbose: bool
) -> TestCaseResult:
    """T1.7: Single-point pure-PSP forces (GTH-H2)."""
    res = TestCaseResult("T1.7", 1, "Single-point pure-PSP forces (GTH-H2)")
    with tempfile.TemporaryDirectory(prefix="madqc_e2e_t1_7_") as tmpdir:
        deck = os.path.join(tmpdir, "h2_pure_psp_forces.in")
        with open(deck, "w") as f:
            f.write(
                """dft
  xc lda
  protocol [1e-4]
  dconv 1e-5
  k 7
  derivatives true
  prefix h2_pure_psp
end

geometry
  units bohr
  no_orient true
  psH 0.0 0.0 -0.6
  psH 0.0 0.0  0.6
end
"""
            )
        rc, out = run_command([madqc, deck], tmpdir, timeout=60, verbose=verbose)
        if rc != 0:
            res.error_category = "MADQC_EXEC_ERROR"
            res.details = f"madqc failed with exit code {rc}"
            return res

        info = read_calc_info(tmpdir, "h2_pure_psp")
        grad = get_task_gradient(info)
        if not grad or len(grad) != 2:
            res.error_category = "GRADIENT_MISSING"
            res.details = "properties.gradient missing or invalid in calc_info.json"
            return res

        fz1 = grad[0][2]
        fz2 = grad[1][2]
        res.observed["grad_H1"] = grad[0]
        res.observed["grad_H2"] = grad[1]

        if abs(fz1) < 1e-3 or abs(fz2) < 1e-3:
            res.error_category = "ZERO_FORCES"
            res.details = f"Forces along bond are vanishingly small: fz1={fz1}, fz2={fz2}"
            return res

        # Equal and opposite forces along z
        if abs(fz1 + fz2) > 1e-4:
            res.error_category = "ASYMMETRIC_FORCES"
            res.details = f"Forces not equal and opposite: fz1={fz1}, fz2={fz2}"
            return res

        res.passed = True
        res.details = f"Pure-PSP forces computed: fz(H1)={fz1:.6f}, fz(H2)={fz2:.6f} Ha/Bohr"
        return res


# ============================================================================
# TIER 2: Boundary & Corner Cases (Hamiltonian Rejection, Oracles)
# ============================================================================


def test_t2_1_mismatch_ae_to_psp(madqc: str, verbose: bool) -> TestCaseResult:
    """T2.1: Mismatched Hamiltonian rejection (AE checkpoint -> PSP calculation)."""
    res = TestCaseResult(
        "T2.1", 2, "Mismatched Hamiltonian rejection: AE checkpoint -> PSP run"
    )
    with tempfile.TemporaryDirectory(prefix="madqc_e2e_t2_1_") as tmpdir:
        prefix = "be_cross_ham"
        # Step 1: Run All-Electron Be (Z=4, E ~ -14.45 Ha)
        in_ae = os.path.join(tmpdir, "be_ae.in")
        with open(in_ae, "w") as f:
            f.write(
                f"""dft
  xc lda
  protocol [1e-3]
  dconv 1e-4
  k 6
  maxiter 25
  prefix {prefix}
  save 1
end

geometry
  units bohr
  no_orient true
  Be 0.0 0.0 0.0
end
"""
            )
        rc1, out1 = run_command([madqc, in_ae], tmpdir, verbose=verbose)
        if rc1 != 0:
            res.error_category = "AE_SETUP_FAILED"
            res.details = f"Step 1 AE Be failed with exit code {rc1}"
            return res

        info_ae = read_calc_info(tmpdir, prefix)
        e_ae = get_task_energy(info_ae)
        if e_ae is None or abs(e_ae - REF_ENERGIES["be_lda_ae"]) > 0.5:
            res.error_category = "AE_SETUP_FAILED"
            res.details = f"Step 1 AE Be produced unexpected energy: {e_ae}"
            return res

        # Step 2: In the SAME directory, run Pseudopotential Be (Zeff=2, E ~ -0.993 Ha)
        in_psp = os.path.join(tmpdir, "be_psp.in")
        with open(in_psp, "w") as f:
            f.write(
                f"""dft
  xc lda
  protocol [1e-3]
  dconv 1e-4
  k 6
  maxiter 25
  psp_calc true
  prefix {prefix}
  save 1
end

geometry
  units bohr
  no_orient true
  Be 0.0 0.0 0.0
end
"""
            )
        rc2, out2 = run_command([madqc, in_psp], tmpdir, verbose=verbose)
        if rc2 != 0:
            res.error_category = "PSP_RUN_FAILED"
            res.details = f"Step 2 PSP Be failed with exit code {rc2}"
            return res

        info_psp = read_calc_info(tmpdir, prefix)
        e_psp = get_task_energy(info_psp)
        res.observed["e_reported"] = e_psp
        res.expected["e_ae_stale"] = e_ae
        res.expected["e_psp_true"] = REF_ENERGIES["be_lda_gth"]

        # If madqc reused the stale AE energy:
        if e_psp is not None and abs(e_psp - e_ae) < 0.1:
            res.error_category = "RESTART_STALE_ENERGY_REUSED"
            res.details = (
                f"STALE ENERGY REUSED: madqc returned old AE energy {e_psp:.6f} Ha "
                f"for PSP run (expected {REF_ENERGIES['be_lda_gth']:.6f} Ha)"
            )
            return res

        # Must be close to true PSP energy
        if e_psp is None or abs(e_psp - REF_ENERGIES["be_lda_gth"]) > 0.05:
            res.error_category = "PSP_ENERGY_INCORRECT"
            res.details = (
                f"PSP energy {e_psp} differs from expected {REF_ENERGIES['be_lda_gth']}"
            )
            return res

        res.passed = True
        res.details = f"Successfully rejected AE checkpoint and converged to PSP ground state {e_psp:.6f} Ha"
        return res


def test_t2_2_mismatch_functional(madqc: str, verbose: bool) -> TestCaseResult:
    """T2.2: Mismatched functional rejection (LDA -> PBE in same dir)."""
    res = TestCaseResult(
        "T2.2", 2, "Mismatched functional rejection (LDA -> PBE restart)"
    )
    with tempfile.TemporaryDirectory(prefix="madqc_e2e_t2_2_") as tmpdir:
        prefix = "he_func_switch"
        # Step 1: Run LDA
        in_lda = os.path.join(tmpdir, "he_lda.in")
        with open(in_lda, "w") as f:
            f.write(
                f"""dft
  xc lda
  protocol [1e-3]
  dconv 1e-4
  k 6
  maxiter 25
  prefix {prefix}
  save 1
end

geometry
  units bohr
  no_orient true
  He 0.0 0.0 0.0
end
"""
            )
        rc1, out1 = run_command([madqc, in_lda], tmpdir, verbose=verbose)
        if rc1 != 0:
            res.error_category = "LDA_SETUP_FAILED"
            res.details = f"Step 1 LDA run failed with exit code {rc1}"
            return res
        e_lda = get_task_energy(read_calc_info(tmpdir, prefix))

        # Step 2: Run PBE in same dir
        in_pbe = os.path.join(tmpdir, "he_pbe.in")
        with open(in_pbe, "w") as f:
            f.write(
                f"""dft
  xc pbe
  protocol [1e-3]
  dconv 1e-4
  k 6
  maxiter 25
  prefix {prefix}
  save 1
end

geometry
  units bohr
  no_orient true
  He 0.0 0.0 0.0
end
"""
            )
        rc2, out2 = run_command([madqc, in_pbe], tmpdir, verbose=verbose)
        if rc2 != 0:
            res.error_category = "PBE_RUN_FAILED"
            res.details = f"Step 2 PBE run failed with exit code {rc2}"
            return res
        e_pbe = get_task_energy(read_calc_info(tmpdir, prefix))

        res.observed["e_lda"] = e_lda
        res.observed["e_pbe"] = e_pbe

        if e_lda is None or e_pbe is None:
            res.error_category = "ENERGY_MISSING"
            res.details = "Could not extract energies"
            return res

        # LDA and PBE on Helium differ by ~0.05 Ha. Stale energy reuse would yield delta < 1e-4.
        delta = abs(e_pbe - e_lda)
        if delta < 1e-3:
            res.error_category = "FUNCTIONAL_NOT_RECOMPUTED"
            res.details = f"PBE run reused LDA energy (delta={delta:.2e} Ha)"
            return res

        res.passed = True
        res.details = f"Functional switch recomputed correctly: E_lda={e_lda:.4f}, E_pbe={e_pbe:.4f} (delta={delta:.4f} Ha)"
        return res


def test_t2_3_imported_mol_potential(madqc: str, verbose: bool) -> TestCaseResult:
    """T2.3: Imported molecule parameter precedence (nuclear potential not omitted)."""
    res = TestCaseResult(
        "T2.3", 2, "Potential integrity: verify non-zero nuclear attraction in PSP"
    )
    with tempfile.TemporaryDirectory(prefix="madqc_e2e_t2_3_") as tmpdir:
        deck = os.path.join(tmpdir, "be_psp_pot.in")
        with open(deck, "w") as f:
            f.write(
                """dft
  xc lda
  protocol [1e-3]
  dconv 1e-4
  k 6
  psp_calc true
  prefix be_pot
end

geometry
  units bohr
  no_orient true
  Be 0.0 0.0 0.0
end
"""
            )
        rc, out = run_command([madqc, deck], tmpdir, verbose=verbose)
        if rc != 0:
            res.error_category = "MADQC_EXEC_ERROR"
            res.details = f"madqc exited with code {rc}"
            return res

        energies = parse_stdout_energies(out)
        res.observed["energy_components"] = energies

        # If nuclear attraction was omitted, vnuc == 0
        vnuc = energies.get("nuclear_attraction")
        if vnuc is None or abs(vnuc) < 0.1:
            res.error_category = "POTENTIAL_OMISSION"
            res.details = f"Nuclear attraction potential is zero or missing: {vnuc}"
            return res

        res.passed = True
        res.details = f"Nuclear attraction potential present: V_nuc = {vnuc:.6f} Ha"
        return res


def test_t2_4_pbe_psp_parameters(madqc: str, verbose: bool) -> TestCaseResult:
    """T2.4: PBE pseudopotential parameter handling (Be Zeff=4, purely local lmax=-1)."""
    res = TestCaseResult(
        "T2.4", 2, "PBE pseudopotential parameter handling (Be Zeff=4, lmax=-1)"
    )
    with tempfile.TemporaryDirectory(prefix="madqc_e2e_t2_4_") as tmpdir:
        deck = os.path.join(tmpdir, "be_pbe_psp.in")
        with open(deck, "w") as f:
            f.write(
                """dft
  xc pbe
  protocol [1e-3]
  dconv 1e-4
  k 6
  psp_calc true
  prefix be_pbe
end

geometry
  units bohr
  no_orient true
  Be 0.0 0.0 0.0
end
"""
            )
        rc, out = run_command([madqc, deck], tmpdir, verbose=verbose)
        if rc != 0:
            res.error_category = "EXECUTION_FAILURE"
            res.details = f"Be PBE PSP run crashed with exit code {rc}"
            return res

        info = read_calc_info(tmpdir, "be_pbe")
        energy = get_task_energy(info)
        res.observed["energy"] = energy
        res.expected["ref_energy"] = REF_ENERGIES["be_pbe_gth"]

        if energy is None or abs(energy - REF_ENERGIES["be_pbe_gth"]) > 0.05:
            res.error_category = "ENERGY_MISMATCH"
            res.details = (
                f"Be PBE energy {energy} differs from reference {REF_ENERGIES['be_pbe_gth']}"
            )
            return res

        res.passed = True
        res.details = f"Be PBE GTH converged cleanly to {energy:.6f} Ha (lmax=-1 handled)"
        return res


def test_t2_5_forces_nwchem_h2(
    madqc: str, nwchem: Optional[str], verbose: bool
) -> TestCaseResult:
    """T2.5: Hellmann-Feynman force comparison vs NWChem cc-pVQZ reference (H2)."""
    res = TestCaseResult(
        "T2.5", 2, "Hellmann-Feynman force comparison vs NWChem cc-pVQZ (H2)"
    )
    with tempfile.TemporaryDirectory(prefix="madqc_e2e_t2_5_") as tmpdir:
        deck = os.path.join(tmpdir, "h2_force.in")
        with open(deck, "w") as f:
            f.write(
                """dft
  xc lda
  protocol [1e-4, 1e-6]
  dconv 1e-6
  k 8
  derivatives true
  prefix h2_comp
end

geometry
  units bohr
  no_orient true
  H 0.0 0.0 -0.7
  H 0.0 0.0  0.7
end
"""
            )
        rc, out = run_command([madqc, deck], tmpdir, verbose=verbose)
        if rc != 0:
            res.error_category = "MADQC_EXEC_ERROR"
            res.details = f"madqc failed with exit code {rc}"
            return res

        info = read_calc_info(tmpdir, "h2_comp")
        grad = get_task_gradient(info)
        if not grad or len(grad) != 2:
            res.error_category = "GRADIENT_MISSING"
            res.details = "No gradient in madqc output"
            return res

        ref_grad = NWCHEM_REF_FORCES["h2_lda_1.4bohr"]["gradients"]
        res.observed["madqc_gradients"] = grad
        res.expected["nwchem_gradients"] = ref_grad

        # Compare Cartesian components
        max_delta = 0.0
        for at in range(2):
            for c in range(3):
                delta = abs(grad[at][c] - ref_grad[at][c])
                if delta > max_delta:
                    max_delta = delta

        res.observed["max_force_delta"] = max_delta
        # Tolerance: 1e-4 Ha/Bohr
        tol = 1.0e-4
        if max_delta > tol:
            res.error_category = "FORCE_TOLERANCE_EXCEEDED"
            res.details = f"Max force difference {max_delta:.2e} exceeds tolerance {tol:.2e} Ha/Bohr"
            return res

        res.passed = True
        res.details = f"Forces match NWChem cc-pVQZ to within {max_delta:.2e} Ha/Bohr (< 1e-4 Ha/Bohr)"
        return res


def test_t2_6_forces_nwchem_h2o(
    madqc: str, nwchem: Optional[str], verbose: bool
) -> TestCaseResult:
    """T2.6: Hellmann-Feynman force comparison vs NWChem cc-pVQZ reference (H2O AE)."""
    res = TestCaseResult(
        "T2.6", 2, "Hellmann-Feynman force comparison vs NWChem cc-pVQZ (H2O AE)"
    )
    with tempfile.TemporaryDirectory(prefix="madqc_e2e_t2_6_") as tmpdir:
        deck = os.path.join(tmpdir, "h2o_force.in")
        with open(deck, "w") as f:
            f.write(
                """dft
  xc lda
  protocol [1e-4, 1e-6]
  dconv 1e-6
  k 8
  derivatives true
  prefix h2o_comp
end

geometry
  units bohr
  no_orient true
  O  0.00000000000000  0.00000000000000  0.71002389091700
  H -1.42650808202133  0.00000000000000 -0.35501194545850
  H  1.42650808202133  0.00000000000000 -0.35501194545850
end
"""
            )
        rc, out = run_command([madqc, deck], tmpdir, timeout=180, verbose=verbose)
        if rc != 0:
            res.error_category = "MADQC_EXEC_ERROR"
            res.details = f"madqc failed with exit code {rc}"
            return res

        info = read_calc_info(tmpdir, "h2o_comp")
        grad = get_task_gradient(info)
        if not grad or len(grad) != 3:
            res.error_category = "GRADIENT_MISSING"
            res.details = "No gradient in madqc output"
            return res

        ref_grad = NWCHEM_REF_FORCES["h2o_lda_equilibrium"]["gradients"]
        res.observed["madqc_gradients"] = grad
        res.expected["nwchem_gradients"] = ref_grad

        # Check z-component forces on O and H (basis set sensitivity is highest along non-bonded directions)
        delta_O_z = abs(grad[0][2] - ref_grad[0][2])
        delta_H1_z = abs(grad[1][2] - ref_grad[1][2])
        delta_H2_z = abs(grad[2][2] - ref_grad[2][2])
        max_z_delta = max(delta_O_z, delta_H1_z, delta_H2_z)

        res.observed["delta_O_z"] = delta_O_z
        res.observed["delta_H1_z"] = delta_H1_z
        res.observed["delta_H2_z"] = delta_H2_z

        tol = 1.0e-4
        if max_z_delta > tol:
            res.error_category = "FORCE_TOLERANCE_EXCEEDED"
            res.details = f"Max z-force difference {max_z_delta:.2e} exceeds tolerance {tol:.2e} Ha/Bohr"
            return res

        res.passed = True
        res.details = f"H2O z-forces match NWChem cc-pVQZ to {max_z_delta:.2e} Ha/Bohr (< 1e-4 Ha/Bohr)"
        return res


def test_t2_7_forces_mixed_vs_fd(madqc: str, verbose: bool) -> TestCaseResult:
    """T2.7: Mixed AE/PSP force verification vs Central Finite Differences (water_lda)."""
    res = TestCaseResult(
        "T2.7", 2, "Mixed AE/PSP forces vs Central Finite Differences (water_lda)"
    )
    with tempfile.TemporaryDirectory(prefix="madqc_e2e_t2_7_") as tmpdir:
        delta = 1.0e-3
        z_O = 0.71002389091700
        x_H1 = -1.42650808202133
        z_H1 = -0.35501194545850
        x_H2 = 1.42650808202133
        z_H2 = -0.35501194545850

        def eval_water_energy(
            prefix: str,
            cur_z_O: float,
            cur_x_H1: float,
            cur_z_H1: float,
            cur_x_H2: float,
            cur_z_H2: float,
        ) -> Optional[float]:
            prefix_clean = prefix.lower()
            in_file = os.path.join(tmpdir, f"{prefix_clean}.in")
            with open(in_file, "w") as f:
                f.write(
                    f"""dft
  xc lda
  protocol [1e-4, 1e-6]
  dconv 1e-6
  k 8
  maxiter 30
  prefix {prefix_clean}
end

geometry
  units bohr
  no_orient true
  psO  0.00000000000000  0.00000000000000  {cur_z_O:.14f}
  H   {cur_x_H1:.14f}  0.00000000000000 {cur_z_H1:.14f}
  H    {cur_x_H2:.14f}  0.00000000000000 {cur_z_H2:.14f}
end
"""
                )
            rc, out = run_command([madqc, in_file], tmpdir, timeout=120, verbose=verbose)
            if rc != 0:
                return None
            info = read_calc_info(tmpdir, prefix_clean)
            e = get_task_energy(info)
            if e is None:
                e = parse_stdout_energies(out).get("total")
            return e

        # 1. Analytical run at equilibrium
        in_ana = os.path.join(tmpdir, "water_ana.in")
        with open(in_ana, "w") as f:
            f.write(
                f"""dft
  xc lda
  protocol [1e-4, 1e-6]
  dconv 1e-6
  k 8
  maxiter 30
  derivatives true
  prefix water_ana
end

geometry
  units bohr
  no_orient true
  psO  0.00000000000000  0.00000000000000  {z_O:.14f}
  H   {x_H1:.14f}  0.00000000000000 {z_H1:.14f}
  H    {x_H2:.14f}  0.00000000000000 {z_H2:.14f}
end
"""
            )
        rc_ana, out_ana = run_command([madqc, in_ana], tmpdir, timeout=150, verbose=verbose)
        if rc_ana != 0:
            res.error_category = "ANALYTICAL_EXEC_FAILED"
            res.details = f"Analytical water calculation failed with code {rc_ana}"
            return res

        grad_ana = get_task_gradient(read_calc_info(tmpdir, "water_ana"))
        if not grad_ana or len(grad_ana) != 3:
            res.error_category = "GRADIENT_MISSING"
            res.details = "Could not extract analytical gradient"
            return res

        ana_dE_dz_O = grad_ana[0][2]
        ana_dE_dz_H1 = grad_ana[1][2]

        # 2. Central finite difference for psO along z
        ep_O_z = eval_water_energy("w_oz_p", z_O + delta, x_H1, z_H1, x_H2, z_H2)
        em_O_z = eval_water_energy("w_oz_m", z_O - delta, x_H1, z_H1, x_H2, z_H2)

        # 3. Central finite difference for AE H1 along z
        ep_H1_z = eval_water_energy("w_h1z_p", z_O, x_H1, z_H1 + delta, x_H2, z_H2)
        em_H1_z = eval_water_energy("w_h1z_m", z_O, x_H1, z_H1 - delta, x_H2, z_H2)

        if any(e is None for e in [ep_O_z, em_O_z, ep_H1_z, em_H1_z]):
            res.error_category = "FD_ENERGY_MISSING"
            res.details = "One or more finite difference calculations failed to produce energies"
            return res

        fd_dE_dz_O = (ep_O_z - em_O_z) / (2.0 * delta)
        fd_dE_dz_H1 = (ep_H1_z - em_H1_z) / (2.0 * delta)

        diff_O_z = abs(ana_dE_dz_O - fd_dE_dz_O)
        diff_H1_z = abs(ana_dE_dz_H1 - fd_dE_dz_H1)
        max_diff = max(diff_O_z, diff_H1_z)

        res.observed["analytical"] = {"dE_dz_O": ana_dE_dz_O, "dE_dz_H1": ana_dE_dz_H1}
        res.observed["finite_diff"] = {"dE_dz_O": fd_dE_dz_O, "dE_dz_H1": fd_dE_dz_H1}
        res.observed["diffs"] = {"diff_O_z": diff_O_z, "diff_H1_z": diff_H1_z}
        res.observed["max_diff"] = max_diff

        tol = 2.0e-4
        if max_diff > tol:
            res.error_category = "FD_TOLERANCE_EXCEEDED"
            res.details = f"Max force discrepancy {max_diff:.2e} Ha/Bohr exceeds tolerance {tol:.2e} Ha/Bohr"
            return res

        res.passed = True
        res.details = (
            f"Mixed AE/PSP analytical forces match FD within {max_diff:.2e} Ha/Bohr (< 2.0e-4): "
            f"psO(z) ana={ana_dE_dz_O:.6f} vs fd={fd_dE_dz_O:.6f} (diff={diff_O_z:.2e}), "
            f"H1(z) ana={ana_dE_dz_H1:.6f} vs fd={fd_dE_dz_H1:.6f} (diff={diff_H1_z:.2e})"
        )
        return res


def test_t2_8_forces_pure_psp_vs_fd(madqc: str, verbose: bool) -> TestCaseResult:
    """T2.8: Pure-PSP force verification vs Central Finite Differences (GTH-H2)."""
    res = TestCaseResult(
        "T2.8", 2, "Pure-PSP forces vs Central Finite Differences (GTH-H2)"
    )
    with tempfile.TemporaryDirectory(prefix="madqc_e2e_t2_8_") as tmpdir:
        delta = 1.0e-3
        z0 = 0.6  # Compressed bond R = 1.2 Bohr

        def eval_h2_energy(prefix: str, z2_val: float) -> Optional[float]:
            in_file = os.path.join(tmpdir, f"{prefix}.in")
            with open(in_file, "w") as f:
                f.write(
                    f"""dft
  xc lda
  protocol [1e-4, 1e-6]
  dconv 1e-6
  k 8
  maxiter 30
  prefix {prefix}
end

geometry
  units bohr
  no_orient true
  psH 0.0 0.0 -0.6
  psH 0.0 0.0  {z2_val:.14f}
end
"""
                )
            rc, _ = run_command([madqc, in_file], tmpdir, timeout=60, verbose=verbose)
            if rc != 0:
                return None
            return get_task_energy(read_calc_info(tmpdir, prefix))

        # 1. Analytical run
        in_ana = os.path.join(tmpdir, "h2_psp_ana.in")
        with open(in_ana, "w") as f:
            f.write(
                f"""dft
  xc lda
  protocol [1e-4, 1e-6]
  dconv 1e-6
  k 8
  maxiter 30
  derivatives true
  prefix h2_psp_ana
end

geometry
  units bohr
  no_orient true
  psH 0.0 0.0 -0.6
  psH 0.0 0.0  {z0:.14f}
end
"""
            )
        rc_ana, _ = run_command([madqc, in_ana], tmpdir, timeout=60, verbose=verbose)
        if rc_ana != 0:
            res.error_category = "ANALYTICAL_EXEC_FAILED"
            res.details = f"Analytical GTH-H2 calculation failed with code {rc_ana}"
            return res

        grad_ana = get_task_gradient(read_calc_info(tmpdir, "h2_psp_ana"))
        if not grad_ana or len(grad_ana) != 2:
            res.error_category = "GRADIENT_MISSING"
            res.details = "Could not extract analytical gradient"
            return res

        ana_dE_dz2 = grad_ana[1][2]

        # 2. Finite difference: H2 z-displacement
        ep = eval_h2_energy("h2_z_plus", z0 + delta)
        em = eval_h2_energy("h2_z_minus", z0 - delta)

        if ep is None or em is None:
            res.error_category = "FD_ENERGY_MISSING"
            res.details = "Failed to obtain converged energies for finite difference"
            return res

        fd_dE_dz2 = (ep - em) / (2.0 * delta)
        diff = abs(ana_dE_dz2 - fd_dE_dz2)

        res.observed["analytical_dE_dz2"] = ana_dE_dz2
        res.observed["fd_dE_dz2"] = fd_dE_dz2
        res.observed["diff"] = diff
        res.expected["fd_dE_dz2"] = fd_dE_dz2

        tol = 2.0e-4
        if diff > tol:
            res.error_category = "FD_TOLERANCE_EXCEEDED"
            res.details = (
                f"Pure-PSP force discrepancy {diff:.2e} Ha/Bohr "
                f"exceeds tolerance {tol:.2e} Ha/Bohr (ana={ana_dE_dz2:.6f}, fd={fd_dE_dz2:.6f})"
            )
            return res

        res.passed = True
        res.details = (
            f"Pure-PSP analytical force matches FD to {diff:.2e} Ha/Bohr "
            f"(< 2.0e-4 Ha/Bohr; ana={ana_dE_dz2:.6f}, fd={fd_dE_dz2:.6f})"
        )
        return res


def test_t2_9_forces_nwchem_pspw(
    madqc: str, nwchem: Optional[str], verbose: bool
) -> TestCaseResult:
    """T2.9: Pure-PSP force comparison vs NWChem 7.3.1 pspw gradient (GTH-H2)."""
    res = TestCaseResult(
        "T2.9", 2, "Pure-PSP force comparison vs NWChem 7.3.1 pspw gradient (GTH-H2)"
    )
    with tempfile.TemporaryDirectory(prefix="madqc_e2e_t2_9_") as tmpdir:
        if nwchem and os.path.isfile(nwchem):
            # Step 1: Write H.psp matching madness gth.xml parameters
            h_psp_path = os.path.join(tmpdir, "H.psp")
            with open(h_psp_path, "w") as f:
                f.write(
                    """1
H
   1
  -1
   0.20000000  -4.18023680   0.72507482   0.00000000   0.00000000
H GTH-LDA from madness gth.xml
"""
                )

            # Step 2: Write and execute NWChem PSPW deck
            nw_deck = os.path.join(tmpdir, "h2_psp.nw")
            with open(nw_deck, "w") as f:
                f.write(
                    """title "H2 GTH-LDA gradient"
echo
start h2_psp
memory stack 2000 mb heap 2000 mb global 1000 mb

geometry units bohr noautoz nocenter noautosym
  H    0.00000000000000      0.00000000000000     -0.60000000000000
  H    0.00000000000000      0.00000000000000      0.60000000000000
end

nwpw
  simulation_cell
    boundary_conditions aperiodic
    SC 16.0
    ngrid 80 80 80
  end
  xc lda
  mult 1
  lmbfgs
  tolerances 1.0d-7 1.0d-7
end
task pspw gradient
"""
                )
            rc_nwc, out_nwc = run_command([nwchem, "h2_psp.nw"], tmpdir, timeout=180, verbose=verbose)
            if rc_nwc != 0:
                res.error_category = "NWCHEM_EXEC_FAILED"
                res.details = f"NWChem execution failed with exit code {rc_nwc}"
                return res

            nwc_grad = parse_nwchem_ion_forces(out_nwc)
            if not nwc_grad or len(nwc_grad) != 2:
                res.error_category = "NWCHEM_OUTPUT_PARSE_ERROR"
                res.details = "Could not parse Ion Forces from NWChem output"
                return res
        else:
            nwc_grad = NWCHEM_REF_FORCES["h2_gth_lda_1.2bohr"]["gradients"]

        # Step 3: Run madqc analytical gradient on identical system
        mad_deck = os.path.join(tmpdir, "h2_madqc.in")
        with open(mad_deck, "w") as f:
            f.write(
                """dft
  xc lda
  protocol [1e-4, 1e-6]
  dconv 1e-6
  k 8
  derivatives true
  prefix h2_madqc
end

geometry
  units bohr
  no_orient true
  psH 0.0 0.0 -0.6
  psH 0.0 0.0  0.6
end
"""
            )
        rc_mad, out_mad = run_command([madqc, mad_deck], tmpdir, timeout=60, verbose=verbose)
        if rc_mad != 0:
            res.error_category = "MADQC_EXEC_ERROR"
            res.details = f"madqc failed with exit code {rc_mad}"
            return res

        mad_grad = get_task_gradient(read_calc_info(tmpdir, "h2_madqc"))
        if not mad_grad or len(mad_grad) != 2:
            res.error_category = "MADQC_GRADIENT_MISSING"
            res.details = "Could not extract madqc gradient"
            return res

        res.observed["madqc_gradient"] = mad_grad
        res.expected["nwchem_pspw_gradient"] = nwc_grad

        # Compare z-forces
        delta_z0 = abs(mad_grad[0][2] - nwc_grad[0][2])
        delta_z1 = abs(mad_grad[1][2] - nwc_grad[1][2])
        max_delta = max(delta_z0, delta_z1)

        res.observed["max_force_delta"] = max_delta
        tol = 5.0e-4
        if max_delta > tol:
            res.error_category = "FORCE_TOLERANCE_EXCEEDED"
            res.details = (
                f"Max force difference {max_delta:.2e} Ha/Bohr exceeds tolerance {tol:.2e} Ha/Bohr "
                f"(madqc={mad_grad[1][2]:.6f}, nwc={nwc_grad[1][2]:.6f})"
            )
            return res

        res.passed = True
        res.details = (
            f"Pure-PSP forces match NWChem 7.3.1 pspw to {max_delta:.2e} Ha/Bohr "
            f"(< 5.0e-4 Ha/Bohr; madqc={mad_grad[1][2]:.6f}, nwc={nwc_grad[1][2]:.6f})"
        )
        return res


# ============================================================================
# TIER 3: Cross-Feature Interactions (Mixed AE/PSP, Invariance)
# ============================================================================


def test_t3_1_mixed_ae_psp_restart(madqc: str, verbose: bool) -> TestCaseResult:
    """T3.1: Mixed all-electron and pseudopotential system restart (water_lda with psO and AE H)."""
    res = TestCaseResult(
        "T3.1", 3, "Mixed AE/PSP system restart (water_lda with psO, AE H)"
    )
    with tempfile.TemporaryDirectory(prefix="madqc_e2e_t3_1_") as tmpdir:
        in_s1 = os.path.join(tmpdir, "water_mixed_s1.in")
        with open(in_s1, "w") as f:
            f.write(
                """dft
  xc lda
  protocol [1e-3]
  dconv 1e-4
  k 7
  prefix water_mixed
  save 1
end

geometry
  units bohr
  no_orient true
  psO  0.00000000000000  0.00000000000000  0.71002389091700
  H   -1.42650808202133  0.00000000000000 -0.35501194545850
  H    1.42650808202133  0.00000000000000 -0.35501194545850
end
"""
            )
        rc1, out1 = run_command([madqc, in_s1], tmpdir, verbose=verbose)
        if rc1 != 0:
            res.error_category = "STAGE1_FAILED"
            res.details = f"Stage 1 mixed water run failed with exit code {rc1}"
            return res

        energies1 = parse_stdout_energies(out1)
        res.observed["stage1_energies"] = energies1

        # Both nonlocal PSP (from psO) and nuclear attraction (from AE H) must be present
        if (
            energies1.get("nonlocal_psp") is None
            or energies1.get("nuclear_attraction") is None
        ):
            res.error_category = "POTENTIAL_MISSING"
            res.details = "Either nonlocal PSP or nuclear attraction energy is missing in mixed calculation"
            return res

        # Stage 2: Restart finer
        in_s2 = os.path.join(tmpdir, "water_mixed_s2.in")
        with open(in_s2, "w") as f:
            f.write(
                """dft
  xc lda
  protocol [1e-3, 1e-5]
  dconv 1e-5
  k 7
  prefix water_mixed
  save 1
end

geometry
  units bohr
  no_orient true
  psO  0.00000000000000  0.00000000000000  0.71002389091700
  H   -1.42650808202133  0.00000000000000 -0.35501194545850
  H    1.42650808202133  0.00000000000000 -0.35501194545850
end
"""
            )
        rc2, out2 = run_command([madqc, in_s2], tmpdir, verbose=verbose)
        if rc2 != 0:
            res.error_category = "STAGE2_FAILED"
            res.details = f"Stage 2 mixed water restart failed with exit code {rc2}"
            return res

        reloaded = (
            "reading initial orbitals from restartdata" in out2
            or "from restartdata" in out2
        )
        if not reloaded:
            res.error_category = "RELOAD_NOT_TRIGGERED"
            res.details = "Mixed system did not reload from restartdata"
            return res

        res.passed = True
        res.details = "Mixed AE/PSP system successfully assembled potentials and restarted"
        return res


def test_t3_2_refine_cross_functional(madqc: str, verbose: bool) -> TestCaseResult:
    """T3.2: Protocol refinement combined with cross-functional restart."""
    res = TestCaseResult(
        "T3.2", 3, "Refinement combined with functional mismatch detection"
    )
    with tempfile.TemporaryDirectory(prefix="madqc_e2e_t3_2_") as tmpdir:
        prefix = "h2_cross_refine"
        # Run 1: Coarse LDA
        in_1 = os.path.join(tmpdir, "h2_c_lda.in")
        with open(in_1, "w") as f:
            f.write(
                f"""dft
  xc lda
  protocol [1e-3]
  dconv 1e-4
  k 6
  maxiter 25
  prefix {prefix}
  save 1
end

geometry
  units bohr
  no_orient true
  H 0.0 0.0 -0.7
  H 0.0 0.0  0.7
end
"""
            )
        rc1, out1 = run_command([madqc, in_1], tmpdir, verbose=verbose)
        if rc1 != 0:
            res.error_category = "RUN1_FAILED"
            res.details = f"Run 1 coarse LDA failed with exit code {rc1}"
            return res
        e1 = get_task_energy(read_calc_info(tmpdir, prefix))

        # Run 2: Fine PBE in same directory
        in_2 = os.path.join(tmpdir, "h2_f_pbe.in")
        with open(in_2, "w") as f:
            f.write(
                f"""dft
  xc pbe
  protocol [1e-3, 1e-5]
  dconv 1e-5
  k 6
  maxiter 25
  prefix {prefix}
  save 1
end

geometry
  units bohr
  no_orient true
  H 0.0 0.0 -0.7
  H 0.0 0.0  0.7
end
"""
            )
        rc2, out2 = run_command([madqc, in_2], tmpdir, verbose=verbose)
        if rc2 != 0:
            res.error_category = "RUN2_FAILED"
            res.details = f"Run 2 fine PBE failed with exit code {rc2}"
            return res
        e2 = get_task_energy(read_calc_info(tmpdir, prefix))

        res.observed["e_lda_coarse"] = e1
        res.observed["e_pbe_fine"] = e2

        if e1 is None or e2 is None or abs(e1 - e2) < 1e-3:
            res.error_category = "FUNCTIONAL_NOT_RECOMPUTED"
            res.details = "Functional mismatch not respected across protocol change"
            return res

        res.passed = True
        res.details = f"Successfully combined protocol refinement and functional switch: E_lda={e1:.4f}, E_pbe={e2:.4f}"
        return res


def test_t3_3_forces_translational_invariance(
    madqc: str, verbose: bool
) -> TestCaseResult:
    """T3.3: Zero net force / translational invariance check for isolated molecules."""
    res = TestCaseResult(
        "T3.3", 3, "Zero net force translational invariance (H2 and H2O)"
    )
    with tempfile.TemporaryDirectory(prefix="madqc_e2e_t3_3_") as tmpdir:
        # Run H2
        deck_h2 = os.path.join(tmpdir, "h2_trans.in")
        with open(deck_h2, "w") as f:
            f.write(
                """dft
  xc lda
  protocol [1e-4]
  dconv 1e-6
  k 8
  derivatives true
  prefix h2_trans
end

geometry
  units bohr
  no_orient true
  H 0.0 0.0 -0.7
  H 0.0 0.0  0.7
end
"""
            )
        rc_h2, _ = run_command([madqc, deck_h2], tmpdir, verbose=verbose)
        if rc_h2 != 0:
            res.error_category = "H2_EXEC_FAILED"
            res.details = f"H2 calculation failed with code {rc_h2}"
            return res

        info_h2 = read_calc_info(tmpdir, "h2_trans")
        grad_h2 = get_task_gradient(info_h2)
        if not grad_h2:
            res.error_category = "GRADIENT_MISSING"
            res.details = "No gradient found in H2 calculation"
            return res

        # Net force vector
        net_f = [sum(grad_h2[at][c] for at in range(len(grad_h2))) for c in range(3)]
        net_f_norm = sum(x * x for x in net_f) ** 0.5
        res.observed["h2_net_force_norm"] = net_f_norm

        tol = 1.0e-5  # Wavelet basis guarantees vanishing Pulay force
        if net_f_norm > tol:
            res.error_category = "TRANSLATIONAL_INVARIANCE_VIOLATED"
            res.details = (
                f"Net force on isolated H2 {net_f_norm:.2e} exceeds tolerance {tol:.2e}"
            )
            return res

        res.passed = True
        res.details = (
            f"Translational invariance confirmed: net force norm = {net_f_norm:.2e} Ha/Bohr (< 1e-5)"
        )
        return res


def test_t3_4_forces_translational_invariance_psp(
    madqc: str, verbose: bool
) -> TestCaseResult:
    """T3.4: Translational Invariance for GTH and mixed AE/PSP molecules (water_lda and ethylene_allpsp_lda)."""
    res = TestCaseResult(
        "T3.4", 3, "Translational Invariance for mixed and pure-PSP molecules"
    )
    with tempfile.TemporaryDirectory(prefix="madqc_e2e_t3_4_") as tmpdir:
        # Part A: water_lda (mixed psO + AE H)
        deck_water = os.path.join(tmpdir, "water_trans.in")
        with open(deck_water, "w") as f:
            f.write(
                """dft
  xc lda
  protocol [1e-4, 1e-6]
  dconv 1e-6
  k 8
  maxiter 30
  derivatives true
  prefix water_trans
end

geometry
  units bohr
  no_orient true
  psO  0.00000000000000  0.00000000000000  0.71002389091700
  H   -1.42650808202133  0.00000000000000 -0.35501194545850
  H    1.42650808202133  0.00000000000000 -0.35501194545850
end
"""
            )
        rc_w, _ = run_command([madqc, deck_water], tmpdir, timeout=150, verbose=verbose)
        if rc_w != 0:
            res.error_category = "WATER_EXEC_FAILED"
            res.details = f"water_lda calculation failed with code {rc_w}"
            return res

        grad_w = get_task_gradient(read_calc_info(tmpdir, "water_trans"))
        if not grad_w or len(grad_w) != 3:
            res.error_category = "GRADIENT_MISSING"
            res.details = "No gradient in water_lda calculation"
            return res

        net_w = [sum(grad_w[a][c] for a in range(len(grad_w))) for c in range(3)]
        net_w_norm = sum(x * x for x in net_w) ** 0.5
        res.observed["water_net_force_norm"] = net_w_norm

        # Part B: ethylene_allpsp_lda (pure PSP)
        deck_eth = os.path.join(tmpdir, "eth_trans.in")
        with open(deck_eth, "w") as f:
            f.write(
                """dft
  xc lda
  psp_calc true
  L 30.0
  k 7
  protocol [1e-4]
  dconv 1e-5
  maxiter 30
  derivatives true
  prefix eth_trans
end

geometry
  units bohr
  no_orient true
  C  1.25925675384707  0.00000000000000  0.00000000000000
  C -1.25925675384707  0.00000000000000  0.00000000000000
  H  2.34326800978054  1.75410867131933  0.00000000000000
  H  2.34326800978054 -1.75410867131933  0.00000000000000
  H -2.34326800978054  1.75410867131933  0.00000000000000
  H -2.34326800978054 -1.75410867131933  0.00000000000000
end
"""
            )
        rc_e, _ = run_command([madqc, deck_eth], tmpdir, timeout=120, verbose=verbose)
        if rc_e != 0:
            res.error_category = "ETHYLENE_EXEC_FAILED"
            res.details = f"ethylene_allpsp_lda calculation failed with code {rc_e}"
            return res

        grad_e = get_task_gradient(read_calc_info(tmpdir, "eth_trans"))
        if not grad_e or len(grad_e) != 6:
            res.error_category = "GRADIENT_MISSING"
            res.details = "No gradient in ethylene calculation"
            return res

        net_e = [sum(grad_e[a][c] for a in range(len(grad_e))) for c in range(3)]
        net_e_norm = sum(x * x for x in net_e) ** 0.5
        res.observed["ethylene_net_force_norm"] = net_e_norm

        tol = 2.0e-4
        if net_w_norm > tol:
            res.error_category = "TRANSLATIONAL_INVARIANCE_VIOLATED"
            res.details = (
                f"Net force norm on mixed water {net_w_norm:.2e} exceeds tolerance {tol:.2e} Ha/Bohr"
            )
            return res

        if net_e_norm > tol:
            res.error_category = "TRANSLATIONAL_INVARIANCE_VIOLATED"
            res.details = (
                f"Net force norm on pure-PSP ethylene {net_e_norm:.2e} exceeds tolerance {tol:.2e} Ha/Bohr"
            )
            return res

        res.passed = True
        res.details = (
            f"Translational invariance confirmed: water ||sum F|| = {net_w_norm:.2e} Ha/Bohr, "
            f"ethylene ||sum F|| = {net_e_norm:.2e} Ha/Bohr (< 2.0e-4 Ha/Bohr)"
        )
        return res


# ============================================================================
# TIER 4: Real-World Scenarios (Multi-Atom Checkpointing, Physical Restoring Force)
# ============================================================================


def test_t4_1_multi_atom_restart(madqc: str, verbose: bool) -> TestCaseResult:
    """T4.1: Multi-atom molecule (methane) restarted from checkpoint."""
    res = TestCaseResult(
        "T4.1", 4, "Multi-atom molecule checkpoint restart (Methane, GTH-LDA)"
    )
    with tempfile.TemporaryDirectory(prefix="madqc_e2e_t4_1_") as tmpdir:
        in_s1 = os.path.join(tmpdir, "methane_s1.in")
        with open(in_s1, "w") as f:
            f.write(
                """dft
  xc lda
  protocol [1e-3]
  dconv 1e-4
  k 7
  prefix methane_chk
  save 1
end

geometry
  units bohr
  no_orient true
  psC  0.00000000000000  0.00000000000000  0.00000000000000
  H    1.18731518047970  1.18731518047970  1.18731518047970
  H   -1.18731518047970 -1.18731518047970  1.18731518047970
  H   -1.18731518047970  1.18731518047970 -1.18731518047970
  H    1.18731518047970 -1.18731518047970 -1.18731518047970
end
"""
            )
        rc1, out1 = run_command([madqc, in_s1], tmpdir, timeout=180, verbose=verbose)
        if rc1 != 0:
            res.error_category = "STAGE1_FAILED"
            res.details = f"Methane Stage 1 coarse run failed with exit code {rc1}"
            return res

        in_s2 = os.path.join(tmpdir, "methane_s2.in")
        with open(in_s2, "w") as f:
            f.write(
                """dft
  xc lda
  protocol [1e-3, 1e-5]
  dconv 1e-5
  k 7
  prefix methane_chk
  save 1
end

geometry
  units bohr
  no_orient true
  psC  0.00000000000000  0.00000000000000  0.00000000000000
  H    1.18731518047970  1.18731518047970  1.18731518047970
  H   -1.18731518047970 -1.18731518047970  1.18731518047970
  H   -1.18731518047970  1.18731518047970 -1.18731518047970
  H    1.18731518047970 -1.18731518047970 -1.18731518047970
end
"""
            )
        rc2, out2 = run_command([madqc, in_s2], tmpdir, timeout=240, verbose=verbose)
        if rc2 != 0:
            res.error_category = "STAGE2_FAILED"
            res.details = f"Methane Stage 2 restart run failed with exit code {rc2}"
            return res

        info2 = read_calc_info(tmpdir, "methane_chk")
        e2 = get_task_energy(info2)
        res.observed["energy"] = e2
        res.expected["ref_energy"] = REF_ENERGIES["methane_lda_gth"]

        reloaded = (
            "reading initial orbitals from restartdata" in out2
            or "from restartdata" in out2
        )
        if not reloaded:
            res.error_category = "RELOAD_NOT_TRIGGERED"
            res.details = "Methane did not reload MOs from restartdata"
            return res

        if e2 is None or abs(e2 - REF_ENERGIES["methane_lda_gth"]) > 0.05:
            res.error_category = "ENERGY_MISMATCH"
            res.details = (
                f"Methane energy {e2} differs from ref {REF_ENERGIES['methane_lda_gth']}"
            )
            return res

        res.passed = True
        res.details = f"Methane restarted cleanly from checkpoint, E={e2:.6f} Ha"
        return res


def test_t4_2_forces_restoring_sign(madqc: str, verbose: bool) -> TestCaseResult:
    """T4.2: Geometry displacement restoring force sign check (compressed vs stretched H2)."""
    res = TestCaseResult(
        "T4.2", 4, "Geometry displacement restoring force sign consistency"
    )
    with tempfile.TemporaryDirectory(prefix="madqc_e2e_t4_2_") as tmpdir:
        # 1. Compressed bond: R = 1.2 Bohr (Re ~ 1.45 Bohr) -> Repulsive force (dE/dz1 > 0, dE/dz2 < 0)
        in_comp = os.path.join(tmpdir, "h2_comp.in")
        with open(in_comp, "w") as f:
            f.write(
                """dft
  xc lda
  protocol [1e-4]
  dconv 1e-6
  k 8
  derivatives true
  prefix h2_comp
end

geometry
  units bohr
  no_orient true
  H 0.0 0.0 -0.6
  H 0.0 0.0  0.6
end
"""
            )
        rc_c, _ = run_command([madqc, in_comp], tmpdir, verbose=verbose)
        if rc_c != 0:
            res.error_category = "COMP_FAILED"
            res.details = "Compressed H2 calculation failed"
            return res

        grad_comp = get_task_gradient(read_calc_info(tmpdir, "h2_comp"))
        if not grad_comp:
            res.error_category = "GRADIENT_MISSING"
            res.details = "No gradient in compressed calculation"
            return res

        # 2. Stretched bond: R = 1.8 Bohr -> Attractive restoring force (dE/dz1 < 0, dE/dz2 > 0)
        in_stretch = os.path.join(tmpdir, "h2_stretch.in")
        with open(in_stretch, "w") as f:
            f.write(
                """dft
  xc lda
  protocol [1e-4]
  dconv 1e-6
  k 8
  derivatives true
  prefix h2_stretch
end

geometry
  units bohr
  no_orient true
  H 0.0 0.0 -0.9
  H 0.0 0.0  0.9
end
"""
            )
        rc_s, _ = run_command([madqc, in_stretch], tmpdir, verbose=verbose)
        if rc_s != 0:
            res.error_category = "STRETCH_FAILED"
            res.details = "Stretched H2 calculation failed"
            return res

        grad_stretch = get_task_gradient(read_calc_info(tmpdir, "h2_stretch"))
        if not grad_stretch:
            res.error_category = "GRADIENT_MISSING"
            res.details = "No gradient in stretched calculation"
            return res

        fz_comp_H2 = grad_comp[1][2]  # Atom 2 at +z
        fz_stretch_H2 = grad_stretch[1][2]

        res.observed["fz_comp_H2"] = fz_comp_H2
        res.observed["fz_stretch_H2"] = fz_stretch_H2

        # Physics:
        # At R=1.2 (< Re), repulsive force pushes atom 2 outwards (+z), so -Fz > 0 => dE/dz2 < 0.
        # At R=1.8 (> Re), attractive restoring force pulls atom 2 inwards (-z), so -Fz < 0 => dE/dz2 > 0.
        if fz_comp_H2 >= 0:
            res.error_category = "UNPHYSICAL_FORCE_SIGN"
            res.details = (
                f"Compressed bond has non-repulsive force: dE/dz2 = {fz_comp_H2}"
            )
            return res

        if fz_stretch_H2 <= 0:
            res.error_category = "UNPHYSICAL_FORCE_SIGN"
            res.details = (
                f"Stretched bond has non-attractive force: dE/dz2 = {fz_stretch_H2}"
            )
            return res

        res.passed = True
        res.details = (
            f"Physical restoring forces verified: compressed dE/dz={fz_comp_H2:.4f} < 0 (repulsive), "
            f"stretched dE/dz={fz_stretch_H2:.4f} > 0 (attractive)"
        )
        return res


def test_t4_3_forces_gth_restoring_sign(
    madqc: str, verbose: bool
) -> TestCaseResult:
    """T4.3: Physical restoring force on GTH bond: Compressed (R=1.1) vs Stretched (R=1.8) GTH-H2."""
    res = TestCaseResult(
        "T4.3", 4, "Physical restoring force on GTH bond (compressed vs stretched GTH-H2)"
    )
    with tempfile.TemporaryDirectory(prefix="madqc_e2e_t4_3_") as tmpdir:
        # 1. Compressed bond: R = 1.1 Bohr (Re ~ 1.45 Bohr) -> Repulsive force (dE/dz2 < 0)
        in_comp = os.path.join(tmpdir, "h2_gth_comp.in")
        with open(in_comp, "w") as f:
            f.write(
                """dft
  xc lda
  protocol [1e-4]
  dconv 1e-5
  k 7
  derivatives true
  prefix h2_gth_comp
end

geometry
  units bohr
  no_orient true
  psH 0.0 0.0 -0.55
  psH 0.0 0.0  0.55
end
"""
            )
        rc_c, _ = run_command([madqc, in_comp], tmpdir, timeout=60, verbose=verbose)
        if rc_c != 0:
            res.error_category = "COMP_FAILED"
            res.details = f"Compressed GTH-H2 calculation failed with code {rc_c}"
            return res

        grad_comp = get_task_gradient(read_calc_info(tmpdir, "h2_gth_comp"))
        if not grad_comp:
            res.error_category = "GRADIENT_MISSING"
            res.details = "No gradient in compressed calculation"
            return res

        # 2. Stretched bond: R = 1.8 Bohr -> Attractive restoring force (dE/dz2 > 0)
        in_stretch = os.path.join(tmpdir, "h2_gth_stretch.in")
        with open(in_stretch, "w") as f:
            f.write(
                """dft
  xc lda
  protocol [1e-4]
  dconv 1e-5
  k 7
  derivatives true
  prefix h2_gth_stretch
end

geometry
  units bohr
  no_orient true
  psH 0.0 0.0 -0.9
  psH 0.0 0.0  0.9
end
"""
            )
        rc_s, _ = run_command([madqc, in_stretch], tmpdir, timeout=60, verbose=verbose)
        if rc_s != 0:
            res.error_category = "STRETCH_FAILED"
            res.details = f"Stretched GTH-H2 calculation failed with code {rc_s}"
            return res

        grad_stretch = get_task_gradient(read_calc_info(tmpdir, "h2_gth_stretch"))
        if not grad_stretch:
            res.error_category = "GRADIENT_MISSING"
            res.details = "No gradient in stretched calculation"
            return res

        fz_comp_H2 = grad_comp[1][2]    # Atom 2 at +z
        fz_stretch_H2 = grad_stretch[1][2]

        res.observed["fz_comp_H2"] = fz_comp_H2
        res.observed["fz_stretch_H2"] = fz_stretch_H2

        # Physics:
        # At R=1.1 (< Re), repulsive force pushes atom 2 outwards (+z), so dE/dz2 < 0.
        # At R=1.8 (> Re), attractive restoring force pulls atom 2 inwards (-z), so dE/dz2 > 0.
        if fz_comp_H2 >= 0:
            res.error_category = "UNPHYSICAL_FORCE_SIGN"
            res.details = f"Compressed GTH bond has non-repulsive force: dE/dz2 = {fz_comp_H2}"
            return res

        if fz_stretch_H2 <= 0:
            res.error_category = "UNPHYSICAL_FORCE_SIGN"
            res.details = f"Stretched GTH bond has non-attractive force: dE/dz2 = {fz_stretch_H2}"
            return res

        res.passed = True
        res.details = (
            f"Physical GTH restoring forces verified: compressed dE/dz={fz_comp_H2:.4f} < 0 (repulsive), "
            f"stretched dE/dz={fz_stretch_H2:.4f} > 0 (attractive)"
        )
        return res


ALL_TESTS = [
    # Tier 1
    ("T1.1", test_t1_1_restart_ae_clean),
    ("T1.2", test_t1_2_restart_psp_clean),
    ("T1.3", test_t1_3_protocol_reload),
    ("T1.4", test_t1_4_fresh_baseline),
    ("T1.5", test_t1_5_forces_ae_single_point),
    ("T1.6", test_t1_6_forces_mixed_ae_psp_single_point),
    ("T1.7", test_t1_7_forces_pure_psp_single_point),
    # Tier 2
    ("T2.1", test_t2_1_mismatch_ae_to_psp),
    ("T2.2", test_t2_2_mismatch_functional),
    ("T2.3", test_t2_3_imported_mol_potential),
    ("T2.4", test_t2_4_pbe_psp_parameters),
    ("T2.5", lambda m, v, nw=None: test_t2_5_forces_nwchem_h2(m, nw or DEFAULT_NWCHEM, v)),
    ("T2.6", lambda m, v, nw=None: test_t2_6_forces_nwchem_h2o(m, nw or DEFAULT_NWCHEM, v)),
    ("T2.7", test_t2_7_forces_mixed_vs_fd),
    ("T2.8", test_t2_8_forces_pure_psp_vs_fd),
    ("T2.9", lambda m, v, nw=None: test_t2_9_forces_nwchem_pspw(m, nw or DEFAULT_NWCHEM, v)),
    # Tier 3
    ("T3.1", test_t3_1_mixed_ae_psp_restart),
    ("T3.2", test_t3_2_refine_cross_functional),
    ("T3.3", test_t3_3_forces_translational_invariance),
    ("T3.4", test_t3_4_forces_translational_invariance_psp),
    # Tier 4
    ("T4.1", test_t4_1_multi_atom_restart),
    ("T4.2", test_t4_2_forces_restoring_sign),
    ("T4.3", test_t4_3_forces_gth_restoring_sign),
]


def run_all_tests(
    madqc: str,
    nwchem: Optional[str] = None,
    filter_tier: Optional[int] = None,
    filter_id: Optional[str] = None,
    verbose: bool = False,
) -> List[TestCaseResult]:
    """Execute selected E2E tests and return results."""
    results = []
    nwchem_bin = nwchem or DEFAULT_NWCHEM
    print("=" * 80)
    print("MADNESS madqc 4-Tier E2E Verification Suite")
    print(f"madqc binary:  {madqc}")
    print(f"nwchem binary: {nwchem_bin or 'N/A'}")
    print("=" * 80)

    for test_id, func in ALL_TESTS:
        if filter_id and test_id != filter_id:
            continue
        tier_num = int(test_id[1])
        if filter_tier and tier_num != filter_tier:
            continue

        print(f"\nRunning [{test_id}] ...", flush=True)
        sig = inspect.signature(func)
        if "nw" in sig.parameters:
            res = func(madqc, verbose, nw=nwchem_bin)
        else:
            res = func(madqc, verbose)
        status = "PASS" if res.passed else ("SKIP" if res.skipped else "FAIL")
        color = "\033[92m" if res.passed else "\033[91m"
        reset = "\033[0m"
        print(f"  Result: {color}[{status}]{reset} {res.description}")
        print(f"  Details: {res.details}")
        results.append(res)

    print("\n" + "=" * 80)
    print("Summary of E2E Test Suite Results:")
    print("=" * 80)
    total = len(results)
    passed = sum(1 for r in results if r.passed)
    failed = sum(1 for r in results if not r.passed and not r.skipped)
    skipped = sum(1 for r in results if r.skipped)

    for r in results:
        status_sym = "[OK]  " if r.passed else "[FAIL]"
        print(f"  {status_sym} [{r.test_id}] (Tier {r.tier}): {r.description}")
        if not r.passed:
            print(f"         Reason: {r.error_category} - {r.details}")

    print("-" * 80)
    print(
        f"Total: {total} | Passed: {passed} | Failed: {failed} | Skipped: {skipped}"
    )
    print("=" * 80)
    return results


def main():
    parser = argparse.ArgumentParser(
        description="MADNESS madqc 4-Tier E2E Verification Suite"
    )
    parser.add_argument(
        "--madqc",
        default=DEFAULT_MADQC,
        help="Path to madqc binary (default: %(default)s)",
    )
    parser.add_argument(
        "--nwchem",
        default=DEFAULT_NWCHEM,
        help="Path to nwchem binary (default: %(default)s)",
    )
    parser.add_argument(
        "--tier", type=int, choices=[1, 2, 3, 4], help="Run only tests from specified tier"
    )
    parser.add_argument(
        "--test", dest="test_id", help="Run a specific test by ID (e.g. T1.1, T2.1)"
    )
    parser.add_argument(
        "--smoke", action="store_true", help="Run a fast smoke test subset (T1.1, T1.5, T2.1)"
    )
    parser.add_argument(
        "-v", "--verbose", action="store_true", help="Verbose command execution output"
    )
    parser.add_argument(
        "--reference_directory", help="Path to reference directory (provided by CTest)"
    )
    parser.add_argument(
        "--reference_dir", help="Path to reference directory (alias)"
    )
    args = parser.parse_args()

    args.madqc = os.path.abspath(args.madqc)
    if not os.path.isfile(args.madqc):
        print(f"Error: madqc binary not found at '{args.madqc}'")
        sys.exit(2)

    if args.smoke:
        results = []
        for tid in ["T1.1", "T1.5", "T2.1"]:
            results.extend(
                run_all_tests(
                    args.madqc,
                    args.nwchem,
                    filter_id=tid,
                    verbose=args.verbose,
                )
            )
    else:
        results = run_all_tests(
            args.madqc,
            args.nwchem,
            filter_tier=args.tier,
            filter_id=args.test_id,
            verbose=args.verbose,
        )

    all_passed = all(r.passed for r in results)
    sys.exit(0 if all_passed else 1)


if __name__ == "__main__":
    main()
