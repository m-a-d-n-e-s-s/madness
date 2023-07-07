
# `molresponse` Tutorial

In the following example we show the basic usage of the `molresponse` application.
We will compute the frequency-dependent response of a Be atom.  To compute response 
properties with `molresponse` we first need to compute the ground-state solution with `moldft`.
We will assume that you have already computed the ground-state solution and have the restart files

## Prerequisites: Preparing for the `molresponse` Calculation

1. **MADNESS Environment**: Make sure that MADNESS is properly installed and set up in your environment.

2. **Restart Files**: This tutorial assumes that you already have `moldft` restart files. These files contain the molecular orbitals and are required to start the `molresponse` calculation.

    If you don't have these files yet, you can generate them by running a `moldft` calculation with the appropriate input parameters for your system.

3. **Directory Structure**:  Below is an example of running a series of response calculations for Be atom.
in a range of frequencies.  The directory structure is as follows:

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

### molrespones input and output

In the base directory you have input and output files associated to the `moldft` calculation.
In this example, `dipole_hf_0*` directories contain the input and output files for the `molresponse` calculations.
For the ground-state, `molresponse` reads in the `moldft` restart files.
To define the response specific calculation parameters, we use the `response.in` file.

Below is an example input file running a dipole response calculation
at zero frequency.

```input
response
         print_level  20         # defined   0: no output; 1: final energy; 2: iterations; 3: timings; 10: debug
                plot  false      # defined   turn on plotting of final orbitals. Output format is .vts
   plot_all_orbitals  true       # defined   Turn on 2D plotting of response orbitals
             maxiter  25         # defined   maximum number of iterations
               dconv  1.0000e-04 # defined   recommended values: 1.e-4 < dconv < 1.e-8
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
Above you can read all the basic parameters for the `molresponse` calculation.
The key response specific parameters are:

- dipole: Flag to turn on frequency dependent property calculation
- omega: Incident energy for dynamic response
- xc: XC input line

In this example we have specified a restart_file in the `save_file` parameter.  This will save the restart file in the `dipole_hf_0-000000` directory.  This is useful if you want to restart the response calculation from the static solution.  We will do this in the next example.



### Computing Dynamic Response: Restarting from Static Solution

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

### Nuclear Response

To run a nuclear response calculation instead of adding a `dipole` parameter, add a `nuclear` parameter.  This will compute the response of the system to a nuclear perturbation.

### Next steps

For more details checkout the [molresponse documentation](/src/apps/molresponse/details.md).

