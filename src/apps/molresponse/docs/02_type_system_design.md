# molresponse Type System and Algorithm Design

Design principle: **clarity over abstraction.** When reading the solver,
you should immediately know what solver you're in (FD or ES), what type
definitions are being used (static/full/TDA), and what each algorithmic
step does. The code should read like the math.

> **Authoritative formulas:** this document describes the *type-system
> architecture*. For the exact operator formulas, metric/sign
> conventions, density prefactors, and the legacy source line that is the
> ground truth for each, see **`operator_contracts.md`** — that file is
> the numeric reference and is kept in sync with the validated code.
> Where this design sketch and `operator_contracts.md` differ on a
> number, `operator_contracts.md` wins.

---

## Two Solvers

### FD Response (Frequency-Dependent)

Solves for one response vector at a time. Each (perturbation, frequency)
pair is independent. This is what enables state-parallel execution —
independent vectors can be solved on independent subworlds.

The algorithm for one vector:
1. Construct the right-hand side (perturbation-dependent)
2. Iterate until convergence:
   a. Compute response density from current response functions
   b. Apply potential operators (V0, gamma) to build theta
   c. Apply BSH Green's function update
   d. Check convergence (residual, density change)
3. Save converged state

### ES Response (Excited-State)

Solves for a coupled bundle of response vectors (X_space). The vectors
are coupled through the eigenvalue problem — they share a subspace and
are rotated together.

The algorithm for the bundle:
1. Generate initial guess (X_space)
2. Iterate until convergence:
   a. Compute Lambda * X (potential energy matrix on response space)
   b. Rotate / diagonalize in the subspace
   c. Construct Theta_X (residual-like quantity for BSH update)
   d. Apply BSH Green's function update to each vector
   e. Check convergence (collective across all vectors)
3. Save converged states and eigenvalues

### What the Two Solvers Share

Both solvers need the same categories of building blocks — they just
call them differently (one vector vs. bundle). The building block
definitions are determined by the response type (static/full/TDA),
not by the solver.

---

## Three Response Types

The response type determines how the building blocks are defined.
The key distinction is **the density definition**, which propagates
into every operator.

### Static (x-only, y = x)

- **Storage:** x functions only — {x_i, i = 1..n_occ}
- **Implicit relation:** y_i = x_i
- **Response density:** rho^(1) = 4 * sum_i [phi_i * x_i]
  (spin factor 2 × y=x factor 2)
- **BSH shift:** mu_i = sqrt(-2 * epsilon_i), no frequency shift
- **Used by:** FD solver at omega = 0

### Full (x + y, independent)

- **Storage:** both x and y — {x_i, y_i, i = 1..n_occ}
- **No implicit relation:** x and y are independent
- **Response density:** rho^(1) = 2 * sum_i [phi_i * (x_i + y_i)]
  (spin factor 2; reduces to Static's 4*sum(phi*x) when y=x)
- **BSH shift:** mu+_i = sqrt(-2*(epsilon_i + omega)) for x,
                  mu-_i = sqrt(-2*(epsilon_i - omega)) for y
- **Used by:** FD solver at omega != 0, ES solver (full TDDFT/TDHF)

### TDA (x-only, y = 0)

- **Storage:** x functions only — {x_i, i = 1..n_occ}
- **Implicit relation:** y_i = 0
- **Response density:** rho^(1) = 2 * sum_i [phi_i * x_i]
  (spin factor 2; half of Static, which carries the extra y=x factor)
- **BSH shift:** mu_i = sqrt(-2 * (epsilon_i + omega)) for x only
- **Used by:** ES solver (Tamm-Dancoff approximation)

### Why TDA != Static

Both store only x functions, and both carry the spin factor of 2. But:
- Static assumes y = x → density picks up an *extra* factor of 2
  (`4*sum(phi*x)`)
- TDA assumes y = 0 → density has only the spin factor (`2*sum(phi*x)`)

So Static density is exactly twice TDA density. This changes every
operator that depends on the density (Coulomb, exchange, XC). They are
genuinely different types, not the same type with a flag.

### Why Full FD == Full ES (at the building block level)

Both store x and y. Both use the same density definition. Both use the
same operator definitions. The difference is in the solver: FD processes
one vector at a given omega, ES processes a bundle and solves for omega.
The building blocks don't care which solver is calling them.

---

## Building Blocks by Type

These are the foundational pieces that each type must define. The solver
(FD or ES) calls these — it doesn't know or care which type is active.

### Building Block Catalog

#### 1. Response Density Construction

Given ground-state orbitals {phi_i} and response functions, compute rho^(1).

| Type   | Definition                                    |
|--------|-----------------------------------------------|
| Static | rho^(1) = 4 * sum_i [phi_i * x_i]            |
| Full   | rho^(1) = 2 * sum_i [phi_i * (x_i + y_i)]    |
| TDA    | rho^(1) = 2 * sum_i [phi_i * x_i]            |

(Restricted/closed-shell prefactors. Unrestricted halves the spin
factor: Static 2, Full 1. See `operator_contracts.md` density row.)

#### 2. Coulomb Potential (gamma)

Given the response density, compute the Coulomb potential.

| Type   | Definition                                    |
|--------|-----------------------------------------------|
| Static | gamma = J[4 * sum_i(phi_i * x_i)]             |
| Full   | gamma = J[2 * sum_i(phi_i * (x_i + y_i))]     |
| TDA    | gamma = J[2 * sum_i(phi_i * x_i)]             |

All three call the same Coulomb operator J, but on the type's response
density rho^(1) from #1 above. The Coulomb operator itself is
type-independent.

#### 3. Exchange Potential (HF)

Given ground orbitals and response functions, compute exchange.

| Type   | Acts on                                       |
|--------|-----------------------------------------------|
| Static | K(x) using density with y = x assumption      |
| Full   | K(x, y) using both x and y                    |
| TDA    | K(x) using density with y = 0 assumption      |

This is where the type distinction matters most. The exchange operator
must know the density definition to be correct.

#### 4. XC Kernel (DFT only)

Given the response density, apply the XC kernel f_xc.

| Type   | Definition                                    |
|--------|-----------------------------------------------|
| Static | f_xc * rho^(1)_static                         |
| Full   | f_xc * rho^(1)_full                           |
| TDA    | f_xc * rho^(1)_tda                            |

Same kernel f_xc, different density. For HF, this is a no-op.

#### 5. Potential Assembly (V0 + gamma → theta)

Combine the ground-state potential, Coulomb, exchange, and XC
contributions into the total potential acting on the response functions.

FD solver: theta_x = V0 * x + gamma * phi + K(response) + f_xc(response)
           (and theta_y for full type)

ES solver: same assembly, but applied to the full X_space bundle.

The components are the same; the assembly is type-specific because
which pieces exist (theta_x only vs theta_x + theta_y) depends on
whether we have y functions.

#### 6. BSH Green's Function Application

Apply the bound-state Helmholtz operator to update response functions.

| Type   | Update                                        |
|--------|-----------------------------------------------|
| Static | x_new = -2 * G(mu) * theta_x                 |
| Full   | x_new = -2 * G(mu+) * theta_x                |
|        | y_new = -2 * G(mu-) * theta_y                |
| TDA    | x_new = -2 * G(mu+) * theta_x                |

Where:
- mu   = sqrt(-2 * epsilon_i)              (static)
- mu+  = sqrt(-2 * (epsilon_i + omega))    (full x, TDA)
- mu-  = sqrt(-2 * (epsilon_i - omega))    (full y)

#### 7. Orthogonality Projection

Project out occupied-space components from response functions.

| Type   | Apply to                                      |
|--------|-----------------------------------------------|
| Static | Q * x  (Q = 1 - sum_i |phi_i><phi_i|)        |
| Full   | Q * x, Q * y                                  |
| TDA    | Q * x                                         |

Same projector Q, but applied to the right set of functions.

#### 8. Residual / Convergence

Compute the residual (difference between old and new response functions)
and check against convergence threshold.

| Type   | Measure                                       |
|--------|-----------------------------------------------|
| Static | ||x_new - x_old||                             |
| Full   | max(||x_new - x_old||, ||y_new - y_old||)     |
| TDA    | ||x_new - x_old||                             |

---

## Additional ES-Specific Building Blocks

The ES solver has steps that FD does not:

#### 9. Lambda Matrix Construction (ES only)

Compute the potential energy matrix Lambda on the response subspace:
  Lambda_ij = <X_i | H_response | X_j>

This uses the same operator building blocks (Coulomb, exchange, XC)
but contracts them across the bundle.

#### 10. Subspace Rotation / Diagonalization (ES only)

Diagonalize Lambda in the current subspace to get eigenvalues (excitation
energies) and rotate X_space to the eigenbasis.

#### 11. Theta Construction for Bundle (ES only)

After rotation, construct the residual-like quantity theta for BSH:
  theta_i = Lambda_i * X_i - omega_i * X_i  (schematic)

The operator pieces are the same as FD; the assembly accounts for
the eigenvalue structure.

---

## Perturbation Type (Orthogonal to Response Type)

The perturbation determines the RHS of the response equations and how
many independent channels exist. It does NOT affect operator definitions.

| Perturbation         | Channels | RHS construction              |
|----------------------|----------|-------------------------------|
| Dipole (electric)    | 3 (x,y,z) | mu_a * phi_i (dipole operator)|
| Nuclear displacement | 3*N_atoms | dV_nuc/dR_a * phi_i          |
| [Future types]       | varies   | [type-specific]               |

For FD solver: each channel is solved independently.
For ES solver: perturbation defines the initial guess space.

---

## Spin Model (Orthogonal to Response Type)

| Model       | Ground orbitals  | Response functions per channel |
|-------------|------------------|---------------------------------|
| Closed-shell| {phi_i} doubly occupied | {x_i} (and y_i if full)    |
| Open-shell  | {phi_i^alpha, phi_i^beta} | {x_i^alpha, x_i^beta} etc.|

The operator building blocks need to account for spin (exchange is
spin-dependent), but the solver skeleton iterates the same way.

---

## Solver × Type Combinations

Valid combinations and their use cases:

| Solver | Type   | Use case                               |
|--------|--------|----------------------------------------|
| FD     | Static | Polarizability at omega=0, static Raman|
| FD     | Full   | Dynamic polarizability, hyperpolarizability|
| ES     | Full   | Full TDDFT/TDHF excited states (direct (X,Y) iteration, `ESSolver<Full>`)|
| ES     | TDA    | Tamm-Dancoff excited states            |
| FD     | TDA    | (not standard, but type system allows) |

Full ES is the generic `ESSolver<Full, ClosedShell>` — direct
paired-(X,Y) iteration with the indefinite-metric subspace step.
(A symmetric-reduction `ESSolverFullRPA` (Davidson on u = X+Y via
(A−B)(A+B)u = ω²u) existed as a cross-check but was removed 2026-06 —
not a direction we are pursuing.)

---

## What the Algorithm Looks Like (Pseudocode for Clarity)

### Single-Vector Building Blocks (shared by both solvers)

All building block definitions operate on a single response vector.
This is the atomic unit of work. Both solvers use these same definitions.

```
# These are what type_defs provides — one vector in, one vector out.

rho_i       = type_defs.compute_density_i(phi_i, x_i, y_i)
gamma       = apply_coulomb(rho)           # type-independent operator
K_i         = type_defs.apply_exchange_i(phi, x_i, y_i)
Vxc_i       = type_defs.apply_xc_i(rho)   # no-op for HF
theta_x_i, theta_y_i = type_defs.assemble_potential_i(...)
x_new_i, y_new_i     = type_defs.apply_bsh_i(theta_x_i, theta_y_i, eps_i, omega)
```

Note: the response density rho is a sum over all orbitals, so it
requires a loop (or reduction) over i. But the per-orbital contributions
to rho are defined per-vector by the type. The Coulomb operator
then acts on the total density once.

### FD Solver (one response vector at a time)

Calls the single-vector building blocks directly.

```
solve_fd_response(ground, perturbation, omega, type_defs):
    x, y = initial_guess(ground, perturbation, type_defs)
    for iteration in 1..max_iter:
        # density is a sum over orbitals — loop the per-orbital definition
        rho = sum over i: type_defs.compute_density_i(phi_i, x_i, y_i)
        gamma = apply_coulomb(rho)
        for i in 1..n_occ:
            K_i   = type_defs.apply_exchange_i(phi, x_i, y_i)
            Vxc_i = type_defs.apply_xc_i(rho)
            theta_i = type_defs.assemble_potential_i(x_i, y_i, phi_i, gamma, K_i, Vxc_i)
            x_new_i, y_new_i = type_defs.apply_bsh_i(theta_i, eps_i, omega)
            x_new_i, y_new_i = project_occupied_i(x_new_i, y_new_i, phi)
        residual = type_defs.compute_residual(x, y, x_new, y_new)
        if converged(residual, threshold):
            break
        x, y = x_new, y_new
    return x, y
```

### ES Solver (coupled bundle, same building blocks)

Loops the same single-vector building blocks over each state in the
bundle. The bundle-level operations (Lambda, rotation) sit on top.

```
solve_es_response(ground, num_states, type_defs):
    X = initial_guess_bundle(ground, num_states, type_defs)
    for iteration in 1..max_iter:

        # --- Bundle-level: Lambda matrix construction ---
        # For each state s in the bundle, apply operators to X_s
        # using the SAME single-vector building blocks:
        for s in 1..num_states:
            rho_s = sum over i: type_defs.compute_density_i(phi_i, X_s.x_i, X_s.y_i)
            gamma_s = apply_coulomb(rho_s)
            for i in 1..n_occ:
                K_s_i   = type_defs.apply_exchange_i(phi, X_s.x_i, X_s.y_i)
                Vxc_s_i = type_defs.apply_xc_i(rho_s)
                Lambda_X_s_i = type_defs.assemble_potential_i(X_s.x_i, X_s.y_i,
                                   phi_i, gamma_s, K_s_i, Vxc_s_i)
            # ^ This loop over i is the same code path as FD.
            #   Future: macrotask this loop over states s.

        # --- Bundle-level: rotation and diagonalization ---
        Lambda_matrix = compute_inner_products(X, Lambda_X)
        X, omega = rotate_and_diagonalize(X, Lambda_matrix)

        # --- Bundle-level: BSH update (loop same building blocks) ---
        for s in 1..num_states:
            theta_s = construct_theta_from_rotated(Lambda_X_s, omega_s, X_s)
            for i in 1..n_occ:
                X_new_s.x_i, X_new_s.y_i = type_defs.apply_bsh_i(
                    theta_s_i, eps_i, omega_s)
                X_new_s.x_i, X_new_s.y_i = project_occupied_i(
                    X_new_s.x_i, X_new_s.y_i, phi)
            # ^ Again, same single-vector BSH as FD.

        residual = compute_bundle_residual(X, X_new)
        if converged(residual, threshold):
            break
        X = X_new
    return X, omega
```

The key point: the loops over states (s) in the ES solver call the
exact same `type_defs` functions as the FD solver. The ES solver just
wraps them in a loop over the bundle and adds the rotation step.

When macrotasking is introduced later, those `for s in 1..num_states`
loops become parallel dispatch over the same single-vector building
blocks — no new operator implementations needed.

Both solvers read like the math. Both use the same building blocks.
The difference is structural (one vector vs bundle + rotation), not
in the operator definitions.

---

## What Lives Where

| Concern                | Where it belongs                          |
|------------------------|-------------------------------------------|
| Solver skeleton (FD)   | Solver layer — algorithm reads like math  |
| Solver skeleton (ES)   | Solver layer — algorithm reads like math  |
| Type definitions       | Type layer — static/full/TDA building blocks|
| Perturbation RHS       | Perturbation layer — channel definitions  |
| Naming conventions     | Naming layer — encodes state/freq/protocol|
| Restart logic          | Persistence layer — save/load/check       |
| Protocol ramping       | Protocol layer — threshold management     |
| State-parallel         | Orchestration layer — subworld management |
| Property computation   | Property layer — combine states → property|

The current molresponseLib is almost entirely orchestration + persistence.
The solver skeleton and type definitions should be in their own files,
readable without navigating scheduling or restart logic.
