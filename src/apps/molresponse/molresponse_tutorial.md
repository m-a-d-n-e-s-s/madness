---
theme: gaia
_class: lead
paginate: true
backgroundColor: #fff
backgroundImage: url('https://marp.app/assets/hero-background.svg')
style: |
  .columns {
    display: grid;
    grid-template-columns: repeat(2, minmax(0, 1fr));
    gap: 1rem;
  }
marp: true
math: mathjax
---

# What is molresponse?

`molresponse` is one of the response codes to compute frequency-dependent response properties as well as excited-states.  We will focus on the frequency-dependent response in the tutorial.  In the most basic form, the solver solves the fruency-dependent coupled resonse equations.

---

## The solver class `FrequencyResponse`

Within `molresponse` the solver class is `FrequencyResponse` which solves,
$$\boldsymbol{A}\boldsymbol{X}=\boldsymbol{b}$$
in the form,
$$\boldsymbol{X}=\boldsymbol{G(\omega_a)}\star[\boldsymbol{VX}+\boldsymbol{b}]$$
Where $\boldsymbol{X}$ is the response vector containing  $\boldsymbol{x}$ and $\boldsymbol{y}$ transition functions.  And $\boldsymbol{b}$ is the perturbation vector. and $\boldsymbol{G(\omega_a)}$ is the Green's function.

---

## Solving $\boldsymbol{X}=\boldsymbol{G(\omega_a)}\star[\boldsymbol{VX}+\boldsymbol{b}]$

- Constructor for the `FrequencyResponse` class.  

```cpp
 FrequencyResponse calc(World &world, const CalcParams &params, double frequency,
                      RHS_Generator rhs);
```

- The `World` object is the MADNESS world.
- The `CalcParams` object contains the parameters for the calculation.
- The `frequency` is the frequency of the perturbation.
- The `RHS_Generator` is a function that generates the RHS vector $\boldsymbol{b}$.

---

## X_Space class

'X_Space' is the main struct that FrequencyResponse operates on.  It holds the response vector $\boldsymbol{X}$.

```cpp
    typedef std::vector<vector_real_function_3d> response_matrix;
    struct X_space {
        size_t n_states;  // Num. of resp. states
        size_t n_orbitals;// Num. of ground states
        response_space x, y;// x and y transition functions    };
```

- The number of states is defined by the perturbation type.  
- The number of orbitals is the number of ground state orbitals.
- The `response_space` is a struct holding **x** and **y** transition functions.

---

## Two Current Options for the `RHS_Generator`

```cpp
X_space dipole_generator(World &world, FrequencyResponse &calc);
X_space nuclear_generator(World &world, FrequencyResponse &calc);
```

- `dipole = True`
  - Generates the RHS vector for the dipole operator.
  - x,y,z components of the dipole operator are generated making dimension of  `X_space` 3 x num_orbitals .
- nuclear
    - Generates the RHS vector for the nuclear derivative operators.


---

## Input File

```input
response
             maxiter  25         # defined   maximum number of iterations
               dconv  1.0000e-04 # defined   recommended values: 1.e-4 < dconv < 1.e-8
           guess_xyz  false      # defined   TODO : check what this is for
            protocol  [1.0000e-04, 1.0000e-06] # defined   calculation protocol
             restart  false      # defined   Flag to restart scf loop from file
                kain  true       # defined   Turn on Krylov Accelarated Inexact Newton Solver
              maxsub  5          # defined   size of iterative subspace ... set to 0 or 1 to disable
                  xc  hf         # defined   XC input line
                save  true       # defined   if true save orbitals to disk
           save_file  restart_dipole_hf_0-000000 # defined   File name to save orbitals for restart
         first_order  true       # defined   Flag to turn on first order response calc
              dipole  true       # defined   Flag to turn on frequency dependent property calc
               omega  0.0000e+00 # defined   Incident energy for dynamic response
end
```

---

# Prerequisites
<!-- Describe the necessary theoretical knowledge, practical skills, and software requirements -->

---

# Installation and Setup
<!-- Provide a brief guide on how to install and set up MADNESS -->

---

# Key Concepts
<!-- Discuss the key concepts and methods used by molresponse -->

---

# Code Walkthrough
<!-- Introduce the basics of how to use molresponse -->

## Input Files
<!-- Describe how to set up input files -->

## Running Calculations
<!-- Describe how to run calculations -->

## Interpreting Output
<!-- Describe how to interpret the output -->

---

# Example Calculations
<!-- Provide hands-on examples of calculations -->

## Simple System
<!-- Show a calculation on a simple system -->

## Complex System
<!-- Show a calculation on a more complex system -->

---

# Common Issues and Troubleshooting
<!-- Discuss common issues and how to solve them -->

---

# Conclusion and Further Resources
<!-- Wrap up the tutorial and provide links to further resources -->

---

# Exercises
<!-- Provide exercises for the participants -->