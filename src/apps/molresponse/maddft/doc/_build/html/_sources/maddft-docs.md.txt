# MADNESS-QCEngine Documentation

Purpose: To provide a single executable for running chemistry applications. This will be 
the main executable used within QCEngine. 

## Current Functionality

- hf/dft 
    - geometry optimization
- linear response
    - excited states
    - dipole
    - nuclear
- response properties
    - polarizability
    - hyperpolarizability
    - gradients of excited states and response states - optimized for excited states


The executable is designed to be a single executable that can be used to run all of the above from either 
an input file in the MADNESS format or JSON input.

### Examples:

H2 Dipole frequency response calculation with hyperpolarizability calculation


```bash
dft
    econv 0.01
    protocol [0.0001]
end

response
    dipole true
    first_order true
    freq_range [0.0,0.056,0.1]
    kain true
    maxiter 10
    maxsub 10
    omega 0.0
    protocol [0.0001]
    quadratic true
end

geometry
    eprec 0.0001
    field [0.0,0.0,0.0]
    no_orient false
    psp_calc false
    pure_ae true
    symtol -0.01
    H 0 0 0
    H 0 0 0.7414
end

```
Or in JSON format:

```json
{
    "dft": {
        "econv": 0.01,
        "protocol": [
            0.0001
        ]
    },
    "molecule": {
        "geometry": [
            [
                0.0,
                0.0,
                0.0
            ],
            [
                0.0,
                0.0,
                0.7414
            ]
        ],
        "parameters": {
            "eprec": 0.0001,
            "field": [
                0.0,
                0.0,
                0.0
            ],
            "no_orient": false,
            "psp_calc": false,
            "pure_ae": true,
            "symtol": -0.01
        },
        "symbols": [
            "H",
            "H"
        ]
    },
    "response": {
        "dipole": true,
        "first_order": true,
        "freq_range": [
            0.0,
            0.056,
            0.1
        ],
        "kain": true,
        "maxiter": 10,
        "maxsub": 10,
        "omega": 0.0,
        "protocol": [
            0.0001
        ],
        "quadratic": true
    }
}

```

To run the calculation, use the following command:

```bash
maddft input 
```
or

```bash
maddft input.json
```
as well as in parallel:

```bash
mpirun -np 4 maddft input
mpirun -np 4 maddft input.json
```

The outputs will be stored in a directory called `output` in the current working directory.

#### Linear Response Directories

Since each response calculation is independent they are placed in separate directories with
the following naming convention:

```[perturbation]_[frequency]``` so for example the dipole response at 0.056 would be in the directory `dipole_0-056000`









