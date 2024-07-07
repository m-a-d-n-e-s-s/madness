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
