# molresponse Usage Guide

.. toctree:: \:maxdepth: 2

## Introduction

`molresponse` is a command\u2011line application for computing molecular response properties (e.g. polarizability, hyperpolarizability, Raman) using multiresolution analysis (MRA). It reads a structured input file (INI or JSON), sets up the ground-state calculation, and generates response states based on user-defined perturbations or high-level property requests.

## Usage

.. code-block:: console

```
molresponse <input_file>
```

Where `<input_file>` can be:

* A plain-text INI-style file
* A JSON file with pre-defined keys

## Examples

---

Minimal polarizability calculation (INI):

.. code-block:: ini

```
[task]
model      = moldft
actions    = energy, response

[response]
requested_properties = polarizability
dipole               = true
dipole.directions    = xyz
dipole.frequencies   = 0.0
```

JSON equivalent:

.. code-block:: json

```
{
  "task": {
    "model": "moldft",
    "actions": ["energy","response"]
  },
  "response": {
    "requested_properties": ["polarizability"],
    "dipole": true,
    "dipole": {
      "directions": ["x","y","z"],
      "frequencies": [0.0]
    }
  }
}
```

MADQC application example (plain-text INI):

This first `task` block is only used by the MADQC driver, which orchestrates
calls to lower-level MADNESS applications (DFT, response, geometry, etc.).

.. code-block:: ini

```
task
    model moldft
    driver energy
end

dft
    dconv    1e-04
    l        200.0
    maxiter  30
    protocol [0.0001,1e-06]
end

molecule
    eprec      1e-06
    field      [0.0,0.0,0.0]
    no_orient  true
    psp_calc   false
    pure_ae    true
    symtol     -0.01
    O  0  0      0.221665
    H  0  1.4309 -0.88666
    H  0 -1.4309 -0.88666
end

response
    quadratic           true
    dipole.directions   xyz
    dipole.frequencies  [0.0, 0.0198666875]
    kain                true
    step_restriction    false
    maxiter             15
end
```

Parameter Reference

.. glossary::

task Top-level block specifying the electronic-structure method (`model`) and which tasks (`actions`) to run (e.g. `energy`, `response`, `optimize`).

numeric Solver-level parameters controlling convergence thresholds, grid tolerances, SCF settings, and maximum iterations.

geometry Molecule geometry options (source\_type, source\_name, units, symtol, etc.).

response Controls response-state generation:

```
  - **requested_properties** (*string list*)
    High\u2011level alias: `polarizability`, `hyperpolarizability`, `raman`.
  - **dipole** (*bool*)
  - **dipole.directions** (*string* or *char list*)
  - **dipole.frequencies** (*float list*)
  - **nuclear** (*bool*)
  - **nuclear.directions** (*string* or *char list*)
  - **nuclear.frequencies** (*float list*)
  - **nuclear.atom_indices** (*int list*)
  - **quadratic** (*bool*): treat dipole kicks as quadratic perturbations.
```

optimize Geometry optimization parameters (maxiter, algopt, value\_precision, etc.).

## Derived Parameters & Validation

* If `requested_properties` is omitted, the code will derive it from `dipole`, `nuclear`, and `quadratic` flags:

  * *linear* + *dipole* \u2192 `polarizability`
  * *quadratic*         \u2192 `polarizability` & `hyperpolarizability`
  * *quadratic* + *dipole* + *nuclear* \u2192 `raman`

* When `requested_properties` is explicitly set, you **must** supply matching perturbation details:

  * For `polarizability` or `hyperpolarizability`: `dipole.directions` & `dipole.frequencies`
  * For `raman`: both `dipole.*` and `nuclear.*` blocks

## Output

* **Restart data:** writes response-state vectors to disk
* **JSON summary:** merged input+derived settings printed to stdout or `<prefix>.json`
* **Logs:** based on `print_level` parameter

## See also

* \:doc:`/api/response` for the C++ API of response- and state-generation classes
* \:doc:`/tutorials/raman` for a step-by-step Raman example
