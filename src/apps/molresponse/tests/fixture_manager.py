"""
fixture_manager.py — Orchestrate test fixture generation, execution, and validation.

Manages the bridge between git-tracked fixture definitions (systems, properties,
numerical settings) and ephemeral scratch data (checkpoints, response outputs).

Design principles:
  - Fixture definitions in repo are the source of truth
  - Scratch data is generated on demand and cached
  - Every operation is idempotent (skip if valid output exists)
  - Numerical settings are controlled by tier (low/medium/high)
  - Properties are composable: individual fixtures for focused tests,
    combined runs for production (raman→alpha+nuclear, beta→alpha, etc.)

Usage:
    fm = FixtureManager()
    fm.ensure_ground_state("h2o_hf", "medium")
    fm.ensure_response("h2o_alpha_static", "medium")
    result = fm.validate("h2o_alpha_static", "medium")

Scratch layout:
    $GECKO_SCRATCH/molresponse_fixtures/
    ├── systems/<system_name>/<tier>/
    │   ├── moldft.in
    │   ├── moldft.restartdata.00000
    │   ├── moldft.calc_info.json
    │   └── moldft.fock.json
    ├── response/<run_name>/<tier>/
    │   ├── response.in
    │   ├── response_metadata.json
    │   └── *.calc_info.json
    └── reference/dalton/<system>/<basis>/
        ├── *.dal
        ├── *.mol
        └── *.out
"""

from __future__ import annotations

import json
import os
import shutil
import subprocess
from dataclasses import dataclass, field
from datetime import datetime
from pathlib import Path
from typing import Any, Optional


# ---------------------------------------------------------------------------
# Paths
# ---------------------------------------------------------------------------

def _fixtures_dir() -> Path:
    """Return the repo-local fixtures directory."""
    return Path(__file__).parent / "fixtures"


def _scratch_root() -> Path:
    """Return the scratch root for generated data."""
    scratch = os.environ.get("GECKO_SCRATCH")
    if not scratch:
        raise EnvironmentError(
            "GECKO_SCRATCH not set. Source setenv.sh or export it."
        )
    root = Path(scratch) / "molresponse_fixtures"
    root.mkdir(parents=True, exist_ok=True)
    return root


# ---------------------------------------------------------------------------
# Data classes
# ---------------------------------------------------------------------------

@dataclass
class NumericalTier:
    """Numerical settings for one accuracy tier."""

    name: str
    dft: dict[str, Any]
    molecule: dict[str, Any]
    auto_k: int
    effective_dconv: float
    use_case: str

    def to_gecko_set_params(self, xc: str = "hf") -> list[str]:
        """Convert tier settings to gecko create_input set_params list."""
        params = [f"dft.xc={xc}"]
        for key, val in self.dft.items():
            if isinstance(val, list):
                params.append(f"dft.{key}=[{','.join(str(v) for v in val)}]")
            elif isinstance(val, bool):
                params.append(f"dft.{key}={'true' if val else 'false'}")
            else:
                params.append(f"dft.{key}={val}")
        for key, val in self.molecule.items():
            params.append(f"molecule.{key}={val}")
        return params


@dataclass
class SystemFixture:
    """Metadata for a molecular system."""

    name: str
    molecule: str
    method: str
    spin: str
    n_atoms: int
    n_electrons: int
    n_occupied: int
    geometry_xyz: Path
    scf_settings: dict[str, Any]
    # Open-shell
    nopen: int = 0
    n_alpha: Optional[int] = None
    n_beta: Optional[int] = None


@dataclass
class PropertyFixture:
    """A single property test configuration."""

    name: str
    system: str
    property_type: str
    tier: str = "medium"
    frequencies: list[float] = field(default_factory=lambda: [0.0])
    directions: str = "xyz"
    perturbation: str = "dipole"
    response_settings: dict[str, Any] = field(default_factory=dict)
    expected: dict[str, Any] = field(default_factory=dict)
    reference_sources: dict[str, Any] = field(default_factory=dict)

    # For combined runs
    includes_properties: list[str] = field(default_factory=list)


@dataclass
class CombinedRun:
    """A multi-property response run combining several fixtures.

    Example: a "production_h2o" run that computes alpha + beta + raman
    in a single madqc invocation, using the same ground state.
    """

    name: str
    system: str
    tier: str
    property_fixtures: list[str]  # names of PropertyFixture entries
    description: str = ""


# ---------------------------------------------------------------------------
# Loading helpers
# ---------------------------------------------------------------------------

def load_numerical_settings() -> dict[str, NumericalTier]:
    """Load all tiers from numerical_settings.json."""
    path = _fixtures_dir() / "numerical_settings.json"
    with open(path) as f:
        data = json.load(f)

    tiers = {}
    for tier_name, tier_data in data["tiers"].items():
        tiers[tier_name] = NumericalTier(
            name=tier_name,
            dft=tier_data["dft"],
            molecule=tier_data["molecule"],
            auto_k=tier_data["auto_k_at_final_protocol"],
            effective_dconv=tier_data["effective_dconv_at_final_step"],
            use_case=tier_data["use_case"],
        )
    return tiers


def load_system(name: str) -> SystemFixture:
    """Load a system fixture definition from the repo."""
    sys_dir = _fixtures_dir() / "systems" / name
    with open(sys_dir / "system.json") as f:
        data = json.load(f)

    return SystemFixture(
        name=data["name"],
        molecule=data["molecule"],
        method=data["method"],
        spin=data["spin"],
        n_atoms=data["n_atoms"],
        n_electrons=data["n_electrons"],
        n_occupied=data["n_occupied"],
        geometry_xyz=sys_dir / "geometry.xyz",
        scf_settings=data.get("scf_settings", {}),
        nopen=data.get("nopen", 0),
        n_alpha=data.get("n_alpha"),
        n_beta=data.get("n_beta"),
    )


def load_property_fixture(category: str, name: str) -> PropertyFixture:
    """Load a property fixture from the repo.

    Args:
        category: Subdirectory under property/ (e.g. "alpha_static")
        name: JSON filename without extension (e.g. "h2o")
    """
    path = _fixtures_dir() / "property" / category / f"{name}.json"
    with open(path) as f:
        data = json.load(f)

    # Merge top-level response-relevant keys into response_settings
    # so that _build_response_input can find them in one place.
    resp_settings = dict(data.get("response_settings", {}))
    if "num_states" in data and "num_states" not in resp_settings:
        resp_settings["num_states"] = data["num_states"]
    if "method" in data:
        # "tda" → tda=True; "full"/"full_tddft" → tda=False
        if data["method"] == "tda" and "tda" not in resp_settings:
            resp_settings["tda"] = True
        elif data["method"] in ("full", "full_tddft") and "tda" not in resp_settings:
            resp_settings["tda"] = False

    return PropertyFixture(
        name=data["name"],
        system=data["system"],
        property_type=data["property"],
        tier=data.get("tier", "medium"),
        frequencies=data.get("frequencies", [0.0]),
        directions=data.get("directions", "xyz"),
        perturbation=data.get("perturbation", "dipole"),
        response_settings=resp_settings,
        expected=data.get("expected", {}),
        reference_sources=data.get("reference_sources", {}),
    )


# ---------------------------------------------------------------------------
# SLURM profiles
# ---------------------------------------------------------------------------

@dataclass
class ClusterConfig:
    """Configuration for a specific HPC cluster/partition type."""
    name: str
    env_script: Optional[str]
    build_dir: str
    moldft_bin: str
    madqc_bin: str
    molresponse_v3_bin: str
    partitions: dict[str, str]
    default_partition: str
    mem_per_node_gb: int
    cores_per_node: int
    numa_mask: Optional[str]
    account: Optional[str]

    def resolve_bin(self, binary: str) -> str:
        """Resolve a binary name to its full path."""
        return str(Path(self.build_dir) / binary)


@dataclass
class ResourceProfile:
    """SLURM resource requirements for a specific (size_class, tier)."""
    nodes: int
    ntasks_per_node: int
    cpus_per_task: int
    mad_num_threads: int
    wall_time: str


def load_slurm_profiles() -> dict[str, Any]:
    """Load SLURM profiles from slurm_profiles.json."""
    path = _fixtures_dir() / "slurm_profiles.json"
    with open(path) as f:
        return json.load(f)


def load_cluster_config(cluster: str) -> ClusterConfig:
    """Load configuration for a specific cluster."""
    profiles = load_slurm_profiles()
    c = profiles["clusters"][cluster]
    return ClusterConfig(
        name=cluster,
        env_script=c["env_script"],
        build_dir=c["build_dir"],
        moldft_bin=c["moldft_bin"],
        madqc_bin=c["madqc_bin"],
        molresponse_v3_bin=c["molresponse_v3_bin"],
        partitions=c["partitions"],
        default_partition=c["default_partition"],
        mem_per_node_gb=c["mem_per_node_gb"],
        cores_per_node=c["cores_per_node"],
        numa_mask=c.get("numa_mask"),
        account=c.get("account"),
    )


def get_resource_profile(system: str, tier: str) -> ResourceProfile:
    """Look up SLURM resource requirements for a system + tier."""
    profiles = load_slurm_profiles()
    size_class = profiles["size_class_map"].get(system)
    if not size_class:
        # Infer from n_occupied
        sys_fixture = load_system(system)
        n_occ = sys_fixture.n_occupied
        if n_occ <= 2:
            size_class = "tiny"
        elif n_occ <= 5:
            size_class = "small"
        elif n_occ <= 15:
            size_class = "medium_size"
        else:
            size_class = "large"

    rp = profiles["resource_profiles"][size_class][tier]
    return ResourceProfile(**rp)


def build_slurm_script(
    job_name: str,
    work_dir: Path,
    cluster: ClusterConfig,
    resources: ResourceProfile,
    executable: str,
    exec_args: str = "",
    partition_class: str = "short",
) -> str:
    """Generate a SLURM job script.

    All output (stdout, stderr) goes to work_dir so everything is co-located.
    """
    partition = cluster.partitions.get(partition_class, cluster.default_partition)
    total_tasks = resources.nodes * resources.ntasks_per_node
    bin_path = cluster.resolve_bin(executable)

    lines = [
        "#!/usr/bin/env bash",
        f"#SBATCH --job-name={job_name}",
        f"#SBATCH --partition={partition}",
        f"#SBATCH --nodes={resources.nodes}",
        f"#SBATCH --ntasks-per-node={resources.ntasks_per_node}",
        f"#SBATCH --cpus-per-task={resources.cpus_per_task}",
        f"#SBATCH --time={resources.wall_time}",
        f"#SBATCH --output={work_dir}/slurm-%j.out",
        f"#SBATCH --error={work_dir}/slurm-%j.err",
    ]

    if cluster.account:
        lines.append(f"#SBATCH --account={cluster.account}")

    lines.append("")
    lines.append("set -euo pipefail")
    lines.append("")

    # Environment setup
    if cluster.env_script:
        lines.append(f"source {cluster.env_script}")
    lines.append("")

    # Runtime settings
    lines.append(f"export MAD_NUM_THREADS={resources.mad_num_threads}")
    lines.append(f"export OMP_NUM_THREADS={resources.cpus_per_task}")
    lines.append("")

    # Change to work directory
    lines.append(f"cd {work_dir}")
    lines.append("")

    # Build the mpirun command
    mpi_args = ["mpirun"]
    if resources.nodes > 1 or resources.ntasks_per_node > 1:
        mpi_args.append("--map-by numa")
    if cluster.numa_mask:
        mpi_args.extend(["numactl", f"--preferred-many={cluster.numa_mask}"])

    mpi_args.append(bin_path)
    if exec_args:
        mpi_args.append(exec_args)

    lines.append("echo \"SLURM_JOB_START $(date -Iseconds)\"")
    lines.append(f"echo \"EXECUTABLE {bin_path}\"")
    lines.append(f"echo \"WORK_DIR {work_dir}\"")
    lines.append(f"echo \"NTASKS {total_tasks}\"")
    lines.append("")
    lines.append(" ".join(mpi_args))
    lines.append("")
    lines.append("EXIT_CODE=$?")
    lines.append("echo \"SLURM_JOB_END $(date -Iseconds) exit_code=$EXIT_CODE\"")
    lines.append("exit $EXIT_CODE")
    lines.append("")

    return "\n".join(lines)


def _select_partition_class(resources: ResourceProfile) -> str:
    """Pick a partition time class based on wall time."""
    h, m, s = resources.wall_time.split(":")
    total_hours = int(h) + int(m) / 60
    if total_hours <= 0.5:
        return "short"
    elif total_hours <= 4:
        return "medium"
    elif total_hours <= 12:
        return "long"
    else:
        return "extended"


# ---------------------------------------------------------------------------
# Property composition
# ---------------------------------------------------------------------------

# Map of which properties subsume others.  If you run raman, you
# automatically get alpha (dipole) and nuclear responses.
PROPERTY_DEPENDENCIES = {
    "polarizability": {"dipole"},
    "hyperpolarizability": {"dipole", "quadratic"},
    "raman": {"dipole", "nuclear", "quadratic"},
    "excitation_energies": {"dipole", "excited"},
}

# Composite property bundles for production runs
PRODUCTION_BUNDLES = {
    "alpha_beta_raman": {
        "description": "Full linear+quadratic property sweep",
        "response_flags": {
            "dipole": True,
            "dipole.directions": "xyz",
            "nuclear": True,
            "nuclear.directions": "xyz",
            "quadratic": True,
            "kain": True,
        },
    },
    "alpha_only": {
        "description": "Linear dipole polarizability only",
        "response_flags": {
            "dipole": True,
            "dipole.directions": "xyz",
            "kain": True,
        },
    },
    "excited_tda": {
        "description": "TDA excited states",
        "response_flags": {
            "dipole": True,
            "dipole.directions": "z",
            "excited.enable": True,
            "excited.tda": True,
        },
    },
}


def _parse_geometry_xyz(xyz_path: Path) -> list[tuple[str, float, float, float]]:
    """Read a geometry.xyz file and return a list of (symbol, x, y, z) tuples.

    XYZ format:
        Line 1: atom count
        Line 2: comment (ignored)
        Lines 3+: symbol x y z
    """
    lines = xyz_path.read_text().strip().splitlines()
    if len(lines) < 3:
        raise ValueError(f"XYZ file {xyz_path} has fewer than 3 lines")

    n_atoms = int(lines[0].strip())
    atoms: list[tuple[str, float, float, float]] = []
    for line in lines[2:]:
        parts = line.split()
        if len(parts) < 4:
            continue
        symbol = parts[0]
        x, y, z = float(parts[1]), float(parts[2]), float(parts[3])
        atoms.append((symbol, x, y, z))

    if len(atoms) != n_atoms:
        raise ValueError(
            f"XYZ file {xyz_path} declares {n_atoms} atoms but found {len(atoms)}"
        )
    return atoms


def _format_molecule_block(atoms: list[tuple[str, float, float, float]], eprec: float) -> str:
    """Format a MADNESS molecule block from parsed atoms."""
    lines = ["molecule"]
    lines.append(f"    eprec {eprec}")
    lines.append("    units angstrom")
    for symbol, x, y, z in atoms:
        lines.append(f"    {symbol} {x:.12f} {y:.12f} {z:.12f}")
    lines.append("end")
    return "\n".join(lines)


def _build_moldft_input(
    tier_settings: NumericalTier,
    system: SystemFixture,
    xc: str,
) -> str:
    """Build a complete moldft.in file as a string.

    Reads the system's geometry.xyz and combines it with the tier's DFT
    settings to produce a valid MADNESS input file.
    """
    atoms = _parse_geometry_xyz(system.geometry_xyz)
    dft_params = tier_settings.dft

    # System-specific L override (e.g. small systems use smaller boxes)
    box_l = system.scf_settings.get("L", dft_params.get("l"))

    # Build dft block
    dft_lines = ["dft"]
    dft_lines.append(f'    xc "{xc}"')
    if "dconv" in dft_params:
        dft_lines.append(f"    dconv {dft_params['dconv']}")
    if "econv" in dft_params:
        dft_lines.append(f"    econv {dft_params['econv']}")
    if "protocol" in dft_params:
        proto_str = ",".join(str(v) for v in dft_params["protocol"])
        dft_lines.append(f"    protocol [{proto_str}]")
    if box_l is not None:
        dft_lines.append(f"    l {box_l}")
    if "maxiter" in dft_params:
        dft_lines.append(f"    maxiter {dft_params['maxiter']}")
    if "localize" in dft_params:
        dft_lines.append(f'    localize "{dft_params["localize"]}"')
    if dft_params.get("dipole", False):
        dft_lines.append("    dipole true")

    # Open-shell
    if system.nopen > 0:
        dft_lines.append(f"    nopen {system.nopen}")

    dft_lines.append("end")

    # Build molecule block
    eprec = tier_settings.molecule.get("eprec", 1e-4)
    mol_block = _format_molecule_block(atoms, eprec)

    return "\n".join(dft_lines) + "\n" + mol_block + "\n"


def _build_response_input(
    fixture: PropertyFixture,
    tier: NumericalTier,
    archive_path: Path,
    system: SystemFixture,
) -> str:
    """Build a complete response .in file as a string.

    Includes:
    - A dft block with no_compute true and tier protocol settings
    - A response block built from the property fixture
    - A molecule block from the system's geometry.xyz
    """
    atoms = _parse_geometry_xyz(system.geometry_xyz)
    dft_params = tier.dft

    # -- DFT block (no_compute mode) --
    dft_lines = ["dft"]
    dft_lines.append("    no_compute true")

    xc = system.scf_settings.get("xc", "hf")
    dft_lines.append(f'    xc "{xc}"')

    if "protocol" in dft_params:
        proto_str = ",".join(str(v) for v in dft_params["protocol"])
        dft_lines.append(f"    protocol [{proto_str}]")
    if "l" in dft_params:
        dft_lines.append(f"    l {dft_params['l']}")
    if "localize" in dft_params:
        dft_lines.append(f'    localize "{dft_params["localize"]}"')

    dft_lines.append("end")

    # -- Response block --
    resp_lines = ["response"]
    resp_lines.append(f"    archive {archive_path}")

    prop_type = fixture.property_type
    resp_settings = fixture.response_settings

    # Map property type to response flags
    if prop_type in ("polarizability", "hyperpolarizability", "raman"):
        resp_lines.append("    dipole true")
        resp_lines.append(f"    dipole.directions {fixture.directions}")
        freq_str = ",".join(str(f) for f in fixture.frequencies)
        resp_lines.append(f"    dipole.frequencies [{freq_str}]")

    if prop_type == "hyperpolarizability":
        resp_lines.append("    quadratic true")

    if prop_type == "raman":
        resp_lines.append("    nuclear true")
        resp_lines.append(f"    nuclear.directions {fixture.directions}")
        freq_str = ",".join(str(f) for f in fixture.frequencies)
        resp_lines.append(f"    nuclear.frequencies [{freq_str}]")

    if prop_type == "excitation_energies":
        resp_lines.append("    dipole true")
        resp_lines.append("    excited.enable true")
        num_states = resp_settings.get("num_states", 4)
        resp_lines.append(f"    excited.num_states {num_states}")
        # TDA vs full: check method or response_settings
        tda = resp_settings.get("tda", True)
        resp_lines.append(f"    excited.tda {'true' if tda else 'false'}")

    # Common response settings from the fixture
    if resp_settings.get("kain"):
        resp_lines.append("    kain true")
    if "maxiter" in resp_settings:
        resp_lines.append(f"    maxiter {resp_settings['maxiter']}")
    if "dconv" in resp_settings:
        resp_lines.append(f"    dconv {resp_settings['dconv']}")
    if "maxsub" in resp_settings:
        resp_lines.append(f"    maxsub {resp_settings['maxsub']}")

    resp_lines.append("end")

    # -- Molecule block --
    eprec = tier.molecule.get("eprec", 1e-4)
    mol_block = _format_molecule_block(atoms, eprec)

    return "\n".join(dft_lines) + "\n" + "\n".join(resp_lines) + "\n" + mol_block + "\n"


def merge_response_flags(*fixtures: PropertyFixture) -> dict[str, Any]:
    """Merge response settings from multiple property fixtures into a
    single response block for a combined run.

    Later fixtures override earlier ones for conflicting keys.
    """
    merged: dict[str, Any] = {}
    all_frequencies: set[float] = set()

    for fix in fixtures:
        # Accumulate frequencies
        all_frequencies.update(fix.frequencies)

        # Merge response settings
        merged.update(fix.response_settings)

        # Set perturbation flags based on property type
        deps = PROPERTY_DEPENDENCIES.get(fix.property_type, {"dipole"})
        for dep in deps:
            merged[dep] = True

        # Directions: widen to xyz if any fixture wants it
        if fix.directions == "xyz":
            merged["dipole.directions"] = "xyz"
            if "nuclear" in deps:
                merged["nuclear.directions"] = "xyz"

    merged["dipole.frequencies"] = sorted(all_frequencies)
    if "nuclear" in merged:
        merged["nuclear.frequencies"] = sorted(all_frequencies)

    return merged


# ---------------------------------------------------------------------------
# FixtureManager
# ---------------------------------------------------------------------------

class FixtureManager:
    """Central orchestrator for fixture generation and validation."""

    def __init__(self):
        self.fixtures_dir = _fixtures_dir()
        self.scratch = _scratch_root()
        self.tiers = load_numerical_settings()

    # -- Path helpers -------------------------------------------------------

    def ground_state_dir(self, system: str, tier: str) -> Path:
        return self.scratch / "systems" / system / tier

    def response_dir(self, run_name: str, tier: str) -> Path:
        return self.scratch / "response" / run_name / tier

    def reference_dir(self, system: str, basis: str) -> Path:
        return self.scratch / "reference" / "dalton" / system / basis

    # -- Ground state -------------------------------------------------------

    def generate_moldft_input(
        self,
        system: str,
        tier: str,
        xc: str = "hf",
    ) -> Path:
        """Generate a moldft.in using the tier settings and system geometry.

        Returns the path to the generated input file.
        """
        tier_settings = self.tiers[tier]
        sys_fixture = load_system(system)
        out_dir = self.ground_state_dir(system, tier)
        out_dir.mkdir(parents=True, exist_ok=True)

        input_path = out_dir / "moldft.in"

        # Build the actual moldft.in content
        input_content = _build_moldft_input(tier_settings, sys_fixture, xc)
        input_path.write_text(input_content)

        # Also write a manifest for debugging / provenance tracking
        params = tier_settings.to_gecko_set_params(xc=xc)
        if sys_fixture.nopen > 0:
            params.append(f"dft.nopen={sys_fixture.nopen}")

        manifest = {
            "system": system,
            "tier": tier,
            "xc": xc,
            "gecko_set_params": params,
            "geometry_xyz": str(sys_fixture.geometry_xyz),
            "generated_input": str(input_path),
        }
        with open(out_dir / "generation_manifest.json", "w") as f:
            json.dump(manifest, f, indent=2)

        return input_path

    def has_valid_checkpoint(self, system: str, tier: str) -> bool:
        """Check if a valid moldft checkpoint exists.

        Checks for both 'moldft' and 'mad' prefixes since the default
        MADNESS prefix is 'mad'.
        """
        gs_dir = self.ground_state_dir(system, tier)
        for prefix in ("moldft", "mad"):
            restart = gs_dir / f"{prefix}.restartdata.00000"
            calc_info = gs_dir / f"{prefix}.calc_info.json"
            if restart.exists() and calc_info.exists():
                return True
        return False

    def ensure_ground_state(
        self,
        system: str,
        tier: str,
        xc: str = "hf",
        force: bool = False,
    ) -> Path:
        """Ensure a ground-state checkpoint exists, generating if needed.

        Returns the ground-state directory path.
        """
        gs_dir = self.ground_state_dir(system, tier)

        if self.has_valid_checkpoint(system, tier) and not force:
            return gs_dir

        input_path = self.generate_moldft_input(system, tier, xc=xc)

        # TODO: Call madness MCP run_madqc or subprocess moldft here.
        # self._run_moldft(input_path, gs_dir)

        return gs_dir

    # -- Response -----------------------------------------------------------

    def generate_response_input(
        self,
        fixture: PropertyFixture,
        tier: Optional[str] = None,
    ) -> Path:
        """Generate a response input file from a property fixture.

        The response input references the ground-state archive from the
        corresponding system + tier checkpoint.
        """
        tier = tier or fixture.tier
        run_name = fixture.name
        out_dir = self.response_dir(run_name, tier)
        out_dir.mkdir(parents=True, exist_ok=True)

        gs_dir = self.ground_state_dir(fixture.system, tier)
        # Use relative path from the response run dir to the ground state archive
        archive_path = Path(
            f"../../../systems/{fixture.system}/{tier}/moldft.restartdata"
        )

        tier_settings = self.tiers[tier]
        sys_fixture = load_system(fixture.system)

        # Build the actual response.in content
        input_content = _build_response_input(
            fixture, tier_settings, archive_path, sys_fixture
        )
        input_path = out_dir / "response.in"
        input_path.write_text(input_content)

        # Also write a manifest for debugging / provenance tracking
        response_flags = merge_response_flags(fixture)
        response_flags["archive"] = str(archive_path)

        manifest = {
            "fixture": fixture.name,
            "system": fixture.system,
            "tier": tier,
            "response_flags": response_flags,
            "ground_state_dir": str(gs_dir),
        }
        with open(out_dir / "generation_manifest.json", "w") as f:
            json.dump(manifest, f, indent=2)

        return input_path

    def generate_combined_run(
        self,
        run: CombinedRun,
    ) -> Path:
        """Generate a single response input that computes multiple properties.

        This is the production path: one madqc invocation computes alpha +
        beta + raman (or whatever the bundle specifies) from a single
        ground-state checkpoint.
        """
        out_dir = self.response_dir(run.name, run.tier)
        out_dir.mkdir(parents=True, exist_ok=True)

        # Load all constituent property fixtures
        fixtures = []
        for pf_name in run.property_fixtures:
            # Parse "category/name" format
            parts = pf_name.split("/")
            if len(parts) == 2:
                fixtures.append(load_property_fixture(parts[0], parts[1]))
            else:
                raise ValueError(
                    f"Property fixture name must be 'category/name', got: {pf_name}"
                )

        # Merge all response flags
        combined_flags = merge_response_flags(*fixtures)

        gs_dir = self.ground_state_dir(run.system, run.tier)
        combined_flags["archive"] = str(gs_dir / "moldft.restartdata")

        manifest = {
            "run_name": run.name,
            "system": run.system,
            "tier": run.tier,
            "description": run.description,
            "constituent_fixtures": run.property_fixtures,
            "combined_response_flags": combined_flags,
            "ground_state_dir": str(gs_dir),
        }
        with open(out_dir / "generation_manifest.json", "w") as f:
            json.dump(manifest, f, indent=2)

        return out_dir / "response.in"

    def ensure_response(
        self,
        fixture_name: str,
        tier: Optional[str] = None,
        force: bool = False,
    ) -> Path:
        """Ensure a response calculation is complete.

        Generates ground state if needed, then response input, then runs.
        """
        # Parse fixture_name as "category/name"
        parts = fixture_name.split("/")
        fixture = load_property_fixture(parts[0], parts[1])
        tier = tier or fixture.tier

        # Ensure ground state first
        self.ensure_ground_state(fixture.system, tier)

        run_dir = self.response_dir(fixture.name, tier)
        calc_info = run_dir / "response.calc_info.json"

        if calc_info.exists() and not force:
            return run_dir

        self.generate_response_input(fixture, tier=tier)
        # TODO: Run madqc --wf=response here
        return run_dir

    # -- Validation ---------------------------------------------------------

    def validate(
        self,
        fixture_name: str,
        tier: Optional[str] = None,
        tolerance: Optional[float] = None,
    ) -> dict[str, Any]:
        """Compare computed results against expected values in the fixture.

        Returns a dict with pass/fail status and per-property comparisons.
        """
        parts = fixture_name.split("/")
        fixture = load_property_fixture(parts[0], parts[1])
        tier = tier or fixture.tier

        if tolerance is None:
            # Default tolerance based on tier
            tolerance = self.tiers[tier].effective_dconv * 10

        run_dir = self.response_dir(fixture.name, tier)
        calc_info_path = run_dir / "response.calc_info.json"

        if not calc_info_path.exists():
            return {
                "status": "missing",
                "message": f"No results found at {calc_info_path}",
            }

        with open(calc_info_path) as f:
            results = json.load(f)

        # TODO: Extract computed property values from results
        # and compare against fixture.expected
        return {
            "status": "not_implemented",
            "fixture": fixture.name,
            "tier": tier,
            "tolerance": tolerance,
            "expected": fixture.expected,
        }

    # -- Dalton references --------------------------------------------------

    def generate_dalton_reference(
        self,
        system: str,
        basis: str = "aug-cc-pVDZ",
        property_type: str = "alpha",
    ) -> Path:
        """Generate Dalton input files for a reference calculation.

        Uses gecko generate_calc_inputs with code="dalton".
        """
        ref_dir = self.reference_dir(system, basis)
        ref_dir.mkdir(parents=True, exist_ok=True)

        # TODO: Call gecko MCP generate_calc_inputs with code="dalton"
        return ref_dir

    def load_dalton_reference(
        self,
        system: str,
        basis: str,
        property_type: str = "alpha",
    ) -> Optional[dict[str, Any]]:
        """Load parsed Dalton results for a system/basis/property.

        Uses gecko's Dalton output parser if results exist.
        Returns None if no results found.
        """
        ref_dir = self.reference_dir(system, basis)
        # Look for Dalton output files
        out_files = list(ref_dir.glob("*.out"))
        if not out_files:
            return None

        # TODO: Call gecko MCP to parse Dalton output
        return None

    # -- SLURM --------------------------------------------------------------

    def generate_job_script(
        self,
        system: str,
        tier: str,
        cluster: str = "xeonmax",
        calc_type: str = "moldft",
    ) -> Path:
        """Generate a SLURM job script for a ground-state or response calc.

        Returns the path to the generated job.sh.
        """
        cluster_cfg = load_cluster_config(cluster)
        resources = get_resource_profile(system, tier)
        partition_class = _select_partition_class(resources)

        if calc_type == "moldft":
            work_dir = self.ground_state_dir(system, tier)
            job_name = f"moldft_{system}_{tier}"
            executable = cluster_cfg.moldft_bin
            exec_args = "--input=moldft.in"
        elif calc_type == "response":
            # For response, work_dir would be the response run dir
            work_dir = self.ground_state_dir(system, tier)  # placeholder
            job_name = f"response_{system}_{tier}"
            executable = cluster_cfg.madqc_bin
            exec_args = "--wf=response --input=response.in"
        else:
            raise ValueError(f"Unknown calc_type: {calc_type}")

        work_dir.mkdir(parents=True, exist_ok=True)

        script = build_slurm_script(
            job_name=job_name,
            work_dir=work_dir,
            cluster=cluster_cfg,
            resources=resources,
            executable=executable,
            exec_args=exec_args,
            partition_class=partition_class,
        )

        job_path = work_dir / "job.sh"
        job_path.write_text(script)
        return job_path

    def submit_ground_state(
        self,
        system: str,
        tier: str,
        cluster: str = "xeonmax",
        force: bool = False,
    ) -> Optional[str]:
        """Generate inputs + job script, then sbatch.

        Returns the SLURM job ID, or None if skipped.
        """
        gs_dir = self.ground_state_dir(system, tier)

        if self.has_valid_checkpoint(system, tier) and not force:
            return None

        # Ensure moldft.in exists
        self.generate_moldft_input(system, tier)

        # Generate job script
        job_path = self.generate_job_script(system, tier, cluster, "moldft")

        # Submit
        result = subprocess.run(
            ["sbatch", str(job_path)],
            capture_output=True,
            text=True,
            cwd=str(gs_dir),
        )

        if result.returncode != 0:
            print(f"FAIL  sbatch {system}/{tier}: {result.stderr.strip()}")
            return None

        # Parse job ID from "Submitted batch job 12345"
        job_id = result.stdout.strip().split()[-1]

        # Record status
        status = {
            "job_id": job_id,
            "system": system,
            "tier": tier,
            "cluster": cluster,
            "submitted": datetime.now().isoformat(),
            "status": "SUBMITTED",
        }
        with open(gs_dir / "job_status.json", "w") as f:
            json.dump(status, f, indent=2)

        return job_id

    def check_jobs(self, tier: str = "medium") -> list[dict[str, Any]]:
        """Check status of all submitted jobs for a tier.

        Queries sacct/squeue and updates job_status.json files.
        """
        results = []
        for system in self.list_systems():
            gs_dir = self.ground_state_dir(system, tier)
            status_file = gs_dir / "job_status.json"

            if not status_file.exists():
                results.append({
                    "system": system,
                    "tier": tier,
                    "status": "NO_JOB",
                    "checkpoint": self.has_valid_checkpoint(system, tier),
                })
                continue

            with open(status_file) as f:
                status = json.load(f)

            job_id = status["job_id"]

            # Query SLURM for current state
            query = subprocess.run(
                ["sacct", "-j", job_id, "--format=State", "--noheader", "-P"],
                capture_output=True,
                text=True,
            )

            if query.returncode == 0 and query.stdout.strip():
                # sacct may return multiple lines (job + job.batch); take first
                slurm_state = query.stdout.strip().split("\n")[0].strip()
                status["status"] = slurm_state
                with open(status_file, "w") as f:
                    json.dump(status, f, indent=2)
            else:
                slurm_state = status.get("status", "UNKNOWN")

            results.append({
                "system": system,
                "tier": tier,
                "job_id": job_id,
                "status": slurm_state,
                "checkpoint": self.has_valid_checkpoint(system, tier),
            })

        return results

    # -- Inventory ----------------------------------------------------------

    def list_systems(self) -> list[str]:
        """List all defined systems."""
        systems_dir = self.fixtures_dir / "systems"
        return sorted(
            d.name for d in systems_dir.iterdir()
            if d.is_dir() and (d / "system.json").exists()
        )

    def list_property_fixtures(self) -> dict[str, list[str]]:
        """List all property fixtures by category."""
        prop_dir = self.fixtures_dir / "property"
        result = {}
        for cat_dir in sorted(prop_dir.iterdir()):
            if cat_dir.is_dir():
                result[cat_dir.name] = sorted(
                    p.stem for p in cat_dir.glob("*.json")
                )
        return result

    def status(self) -> dict[str, Any]:
        """Report what exists in scratch vs what's defined in fixtures."""
        systems = self.list_systems()
        props = self.list_property_fixtures()

        ground_states = {}
        for sys in systems:
            ground_states[sys] = {}
            for tier in self.tiers:
                ground_states[sys][tier] = self.has_valid_checkpoint(sys, tier)

        return {
            "scratch_root": str(self.scratch),
            "systems": systems,
            "property_fixtures": props,
            "ground_state_status": ground_states,
            "tiers": list(self.tiers.keys()),
        }


# ---------------------------------------------------------------------------
# CLI entry point
# ---------------------------------------------------------------------------

def main():
    """Simple CLI for fixture management."""
    import argparse

    parser = argparse.ArgumentParser(
        description="Manage molresponse_v3 test fixtures"
    )
    sub = parser.add_subparsers(dest="command")

    # status
    sub.add_parser("status", help="Show fixture inventory and scratch status")

    # generate-ground-states
    gen_gs = sub.add_parser(
        "generate-ground-states",
        help="Generate moldft inputs for all systems at a given tier",
    )
    gen_gs.add_argument("--tier", default="medium", choices=["low", "medium", "high"])
    gen_gs.add_argument("--system", default=None, help="Single system (default: all)")
    gen_gs.add_argument("--xc", default="hf")

    # generate-job-scripts
    gen_jobs = sub.add_parser(
        "generate-job-scripts",
        help="Generate SLURM job scripts for all systems at a given tier",
    )
    gen_jobs.add_argument("--tier", default="medium", choices=["low", "medium", "high"])
    gen_jobs.add_argument("--system", default=None, help="Single system (default: all)")
    gen_jobs.add_argument("--cluster", default="xeonmax", choices=["xeonmax", "40core", "96core"])

    # submit
    submit = sub.add_parser(
        "submit",
        help="Generate inputs + job scripts and submit via sbatch",
    )
    submit.add_argument("--tier", default="medium", choices=["low", "medium", "high"])
    submit.add_argument("--system", default=None, help="Single system (default: all)")
    submit.add_argument("--cluster", default="xeonmax", choices=["xeonmax", "40core", "96core"])
    submit.add_argument("--force", action="store_true", help="Re-submit even if checkpoint exists")

    # check
    check = sub.add_parser("check", help="Check status of submitted SLURM jobs")
    check.add_argument("--tier", default="medium", choices=["low", "medium", "high"])

    # validate
    val = sub.add_parser("validate", help="Validate computed results against fixtures")
    val.add_argument("fixture", help="category/name (e.g. alpha_static/h2o)")
    val.add_argument("--tier", default=None)

    args = parser.parse_args()

    fm = FixtureManager()

    if args.command == "status":
        import pprint
        pprint.pprint(fm.status())

    elif args.command == "generate-ground-states":
        systems = [args.system] if args.system else fm.list_systems()
        for sys in systems:
            path = fm.generate_moldft_input(sys, args.tier, xc=args.xc)
            print(f"  {sys}/{args.tier} -> {path}")

    elif args.command == "generate-job-scripts":
        systems = [args.system] if args.system else fm.list_systems()
        for sys in systems:
            # Ensure moldft.in exists first
            fm.generate_moldft_input(sys, args.tier)
            path = fm.generate_job_script(sys, args.tier, args.cluster)
            print(f"  {sys}/{args.tier} -> {path}")

    elif args.command == "submit":
        systems = [args.system] if args.system else fm.list_systems()
        submitted = 0
        skipped = 0
        for sys in systems:
            job_id = fm.submit_ground_state(sys, args.tier, args.cluster, args.force)
            if job_id:
                print(f"  SUBMIT  {sys}/{args.tier}  job_id={job_id}")
                submitted += 1
            else:
                if fm.has_valid_checkpoint(sys, args.tier):
                    print(f"  SKIP    {sys}/{args.tier}  (checkpoint exists)")
                else:
                    print(f"  FAIL    {sys}/{args.tier}")
                skipped += 1
        print(f"\nSubmitted: {submitted}, Skipped: {skipped}")

    elif args.command == "check":
        results = fm.check_jobs(args.tier)
        # Print as a table
        print(f"{'System':<20} {'Status':<15} {'Job ID':<12} {'Checkpoint'}")
        print("-" * 60)
        for r in results:
            chk = "yes" if r["checkpoint"] else "no"
            job_id = r.get("job_id", "-")
            print(f"{r['system']:<20} {r['status']:<15} {job_id:<12} {chk}")

    elif args.command == "validate":
        result = fm.validate(args.fixture, tier=args.tier)
        import pprint
        pprint.pprint(result)

    else:
        parser.print_help()


if __name__ == "__main__":
    main()
