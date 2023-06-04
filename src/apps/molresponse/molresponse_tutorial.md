# Tutorial: Understanding and using the `molresponse` code

 In this tutorial, we will explore the computation of frequency-dependent response properties using the `FrequencyResponse` class within `molresponse`.

## Overview of `molresponse`

`molresponse` is a multiresolution solver designed to compute response properties and excited states of molecular systems.  Solver computes response states starting with a `moldft` calculation, which provides the ground state of the system

## Computing Frequency-Dependent Response with `FrequencyResponse` class

In `molresponse`, the computation of frequency-dependent response properties is achieved by solving the coupled response equations in integral form:

$$\boldsymbol{X}=\boldsymbol{G(\omega_a)}\star[\boldsymbol{VX}+\boldsymbol{P}]$$

Here's a breakdown of the variables involved:

- $\boldsymbol{X}$: The response vector containing the transition functions.
- $\boldsymbol{G(\omega_a)}$: The Green's function.
- $\boldsymbol{V}$: The potential.
- $\boldsymbol{P}$: The perturbation operator (e.g., dipole, nuclear, second-order).
- $\omega_a$: The frequency of the perturbation.

Therefore, response calculation is parameterized by:

1. The ground state calculation with molecular geometry, ground state orbitals, and orbital energies.
2. The perturbation operator.
3. The frequency of the perturbation.

The solver iteratively computes the response vector, starting from an initial $\boldsymbol{X}$, until convergence is reached. Convergence is determined by the bsh residuals of the response vectors and the change in the response density, both falling below defined thresholds.
With convergence, the response vectors are used to computed the frequency-dependent response properties.

### Understanding the `X_space` class

The `X_Space` class is a fundamental component in molresponse as it encapsulates all response vectors computed in a response calculation. The solver is capable of solvin multiple perturbations at the same time.

```cpp
    struct X_space {
        size_t n_states;  // Num. of resp. states
        size_t n_orbitals;// Num. of ground states
        response_space x, y;// x and y transition functions    };
```

- `n_states`: This represents the number of response states. The actual number is defined by the type of perturbation used.

- `n_orbitals`: This refers to the number of ground state orbitals present in the system.

- `response_space x, y`: These are structures that hold the x and y transition functions, respectively. These functions are central to how the system responds to external perturbations.

___

## Solving for **X**

To find X, we iterate the X_Space object through the FrequencyResponse class.

Below is the constructor for the FrequencyResponse class. It takes the frequency and the RHS_Generator (representing the perturbation operator) as inputs:

```cpp
 FrequencyResponse calc(World &world, const CalcParams &params, double frequency,RHS_Generator rhs);
```

Key parameters:

- 'const CalcParams &params': Encapsulates the parameters for the calculation. This include molecular geometry and ground state orbitals.
- 'double frequency': The frequency of the perturbation.
- 'RHS_Generator rhs': The perturbation operator.

---

## Choosing the Perturbation  `RHS_Generator`

Currently, `molresponse` supports two options for the Right Hand Side (RHS) vector generation:

```cpp
X_space dipole_generator(World &world, FrequencyResponse &calc);
X_space nuclear_generator(World &world, FrequencyResponse &calc);
```

**Option 1: Dipole Generator**

- The `dipole_generator` creates the RHS vector for the dipole operator.
- It generates the x, y, and z components of the dipole operator.
- Consequently, the dimension of `X_space` becomes `3 x num_orbitals`.
- To use it, set the `dipole` flag to `True` in the input file

**Option 2: Nuclear Generator**

- The `nuclear_generator` generates the RHS vector for the nuclear derivative operators.
- It is still under development and hasn't been thoroughly tested, so use with caution.

## Prerequisites: Preparing for the `molresponse` Calculation

Before we dive into the `FrequencyResponse` calculation, let's make sure we have everything set up correctly.

1. **MADNESS Environment**: Make sure that MADNESS is properly installed and set up in your environment.

2. **Restart Files**: This tutorial assumes that you already have `moldft` restart files. These files contain the molecular orbitals and are required to start the `molresponse` calculation.

    If you don't have these files yet, you can generate them by running a `moldft` calculation with the appropriate input parameters for your system.

3. **Directory Structure**: It is helpful if your working directory is set up as shown.  Here `dipole_hf_0-000000` is the directory for the response calculation.  The `moldft` restart files are in the parent directory.

    ```makrkdown
    Be/
    ├── dipole_hf_0-000000
    │   ├── response_base.json
    │   ├── response.in
    │   └── restart_dipole_hf_0-000000.00000
    ├── dipole_hf_0-011116
    ├── dipole_hf_0-022233
    ├── dipole_hf_0-033349
    ├── dipole_hf_0-044465
    ├── dipole_hf_0-055582
    ├── dipole_hf_0-066698
    ├── moldft.energies.dat
    ├── moldft.calc_info.json
    ├── moldft.in
    ├── moldft.restartaodata
    └── moldft.restartdata.00000
    ```

Now that we have everything set up, let's dive into the `FrequencyResponse` calculation!

---

## Essential Parameters in the Response Input File

In `molresponse`, the input file for the `FrequencyResponse` calculation is based on the `moldft` input file with an additional `response` section.

Key parameters you need to define:

- **Perturbation Operator**: Set by `dipole`. When `True`, the dipole operator is used.
- **Perturbation Frequency**: Defined by the `omega` parameter.
- **Ground-State Restart File**: By default, `molresponse` looks for `../moldft.restartdata`. You can specify a different path using the `archive` parameter.

Below, we have an example response input file named `response.in`, which calculates the static response for the dipole operator.

```input
response
         print_level  1         
             maxiter  25        
               dconv  1.0000e-04 
            protocol  [1.0000e-04, 1.0000e-06] 
             restart  false      
                kain  true       
              maxsub  5          
                  xc  hf         
                save  true       
           save_file  restart_dipole_hf_0-000000 
         first_order  true       
              dipole  true       
               omega  0.0000e+00 
end
```

In this example we have specfied a restart_file in the `save_file` parameter.  This will save the restart file in the `dipole_hf_0-000000` directory.  This is useful if you want to restart the response calculation from the static solution.  We will do this in the next example.

Notes:

- The default input file name for the `FrequencyResponse` calculation is `response.in`.
- For further details on the parameters, refer to the `molresponse` documentation.
- The output data will be saved in the `response_base.json` file.

---

## Computing Dynamic Response: Restarting from Static Solution

In this example, we compute the frequency-dependent response, starting from a precomputed static solution. This approach can often be more efficient than starting from scratch.

To achieve this, we:

1. Set the `restart` parameter to `true`.
2. Specify the static solution restart file in the `restart_file` parameter.
3. Update the `omega` parameter for the dynamic response frequency.

Here's the revised input file for the `FrequencyResponse` calculation:

```input
response
         print_level  1          # 0: no output; 1: final energy; 2: iterations; 3: timings; 10: debug
   plot_all_orbitals  true       # Turn on 2D plotting of response orbitals
             maxiter  25         # Maximum number of iterations
               dconv  1.0000e-04 # Recommended values: 1.e-4 < dconv < 1.e-8
            protocol  [1.0000e-04, 1.0000e-06] # Calculation protocol
             restart  true       # Flag to restart calculation from file
        restart_file  ../dipole_hf_0-000000/restart_dipole_hf_0-000000 # File to read ground parameters from
                kain  true       # Turn on Krylov Accelerated Inexact Newton Solver
              maxsub  5          # Size of iterative subspace ... set to 0 or 1 to disable
                  xc  hf         # XC input line
                save  true       # If true, save orbitals to disk
           save_file  restart_dipole_hf_0-011116 # File name to save orbitals for restart
         first_order  true       # Flag to turn on first order response calculation
              dipole  true       # Sets RHS to dipole operator 3 x num_orbitals
               omega  1.1116e-02 # Incident energy for dynamic response
end
```

This configuration computes the frequency-dependent response at 0.011116 a.u. frequency, starting from the static solution saved in the `../dipole_hf_0-000000/restart_dipole_hf_0-000000` file.

---

## Conclusion and Next Steps

We've explored the inner workings of `molresponse` and seen it in action through practical examples. By now, you should have a clear understanding of:

1. The theoretical framework that underpins `molresponse`.
2. The core components of the code and their roles.
3. How to set up and perform static and dynamic response calculations.
