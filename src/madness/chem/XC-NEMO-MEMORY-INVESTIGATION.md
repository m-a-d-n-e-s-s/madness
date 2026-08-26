# nemo meta-GGA memory runaway: findings

Investigation log, 2026-08-25/26, on branch `pr-nemo-dft` (base `d8a438201`).
Companion to `XC-MADNESS-ASSESSMENT.md` §5.4, whose first version of this
section is superseded by what follows.

**Status: cause identified, cure measured, not yet committed.**

The open question of the first round — *which term makes the potential large on
the first pass of a protocol rung* — is answered: the semilocal divergence
`div(2 de/dsigma grad rho)`, and the reason it is large is the **refinement
depth** at which the flux is represented, not anything about the flux itself.
nemo's meta-GGA carries every XC intermediate at depth 20–22 where moldft
carries them at 10–13, because `refine_to_common_level` lifts all of them to
τ's level and τ's level is set by U1. A derivative operator turns that extra
depth into pointwise excursions of ±8·10⁴ in the potential (§3, §4).

The cure is to stop representing a multiplicative XC potential at all: the
**weak form**, with the divergence commuted through the Green's function so that
the flux is only ever convolved, never differentiated (§8). Measured on the
repro, same machine, same threads: the baseline is killed at 5.8 GB after 608 s,
the weak form **converges in 447 s** to −2.90966370 against the moldft
reference −2.909663718, and the potential's pointwise range goes from
±10⁵ to O(1).

Three cheaper reformulations were tried first and all three failed. Every one is
numerically correct — each reproduces −2.9096636x — and none of them touches the
runaway. They are recorded in §7 with the numbers that killed them, because two
of them are the obvious things to try and each costs a day.

Read §6 (the methodological trap) before adding a probe of your own. The first
round drew three wrong conclusions from node counts; this round only got
anywhere by measuring **pointwise min/max** alongside them.

---

## 1. Symptom

He, TPSS, restricted, `--wf=nemo`, `protocol [1e-4, 1e-6, 1e-8]`. node26,
`MAD_NUM_THREADS=20`.

| path | wall | outcome |
|---|---|---|
| moldft (`--wf=scf`) | 2 m 52 s | converged, -2.909663718 |
| nemo | 56 m | **killed** at 167 GB of 187 GB, no progress after a fock-matrix build |

Under a `ulimit -v` the run reliably walks up to whatever ceiling it is given —
40 GB in ~130 s, 60 GB in ~372 s, 12/20/30 GB proportionately sooner. It does
not converge on a size. That is the signature of a runaway, not of a fixed
over-refinement, and it is the single most important qualitative fact.

**nemo GGA is unaffected.** The same system and protocol with `xc pbe` converges
in 32 s. So τ is necessary to trigger it — and §4 shows why: τ sets the
refinement depth, and PBE's depth is 14–17 where TPSS's is 20–22.

Nothing registered goes near this. `nemo_he_tpss` runs in 24 s at thresh 1e-6.

---

## 2. Established: the chain into the orbital update

Instrumented in `Nemo::solve` around the BSH/KAIN block, first iteration of the
1e-8 rung:

```
UPD iter  0  nemo  (in)          tree     457   depth  8
UPD iter  0  Vnemo (in)          tree   60185   depth 22      <-- 132x the orbital
UPD iter  0  after BSH apply     tree   78745   depth 15
UPD iter  0  after KAIN          tree   78745
UPD iter  0  after truncate      tree   29825                 <-- next iteration's orbital
UPD iter  0  after normalize     tree   29825
UPD iter  0  after orthonormalize tree  29825
```

So the orbital inherits the potential's tree. `truncate` at `thresh` cuts 78745
to 29825, still 65x the 457-node orbital it started from, and that 29825-node
nemo is what arrives at `prep_xc_args` on the next pass:

```
STEP ==== first over-limit object: IN arho ====
STEP IN amo (orbitals)   tree 29681      (was 457 on the previous pass)
STEP IN arho             tree 29697      (was 1033)
```

Everything inside the XC code downstream of that is a faithful consequence.
The loop is: large `Vnemo` -> BSH -> truncate -> larger orbital -> larger
density and τ -> larger `Vnemo`. Each cycle multiplies.

The same disparity is present one rung earlier and merely survivable: at thresh
1e-6, `Vnemo` is 16409-24497 nodes against a 393-881-node orbital.

---

## 3. Answered: which term, and why it is large on the first pass

A per-term census inside `make_xc_potential` and after
`compute_nemo_potentials`, printing norm / tree / depth / **pointwise min and
max**, on every pass. First pass of the 1e-8 rung, baseline nemo TPSS:

```
XCPROBE tau_a                norm 3.18e-01  tree  8265  depth 22
XCPROBE local(de/drho)       norm 2.65e+00  tree  8265  depth 22  min -1.70e+00  max  1.6e-17
XCPROBE vtau                 norm 3.86e-02  tree  8265  depth 22  min -7.2e-20   max  7.00e-02
XCPROBE flux_same[x]         norm 8.15e-02  tree  8265  depth 22  min -6.52e-02  max  6.52e-02
XCPROBE div(flux_same)       norm 8.80e-01  tree 61121  depth 22  min -6.49e+04  max  8.23e+04
XCPROBE xc_pot               norm 2.89e+00  tree 61121  depth 22  min -8.23e+04  max  6.49e+04

VNEMO nemo    norm 2.52e+00  tree   457  depth  8
VNEMO Unemo   norm 3.13e+00  tree  1801  depth 20
VNEMO Jnemo   norm 1.89e+00  tree  2185  depth 11
VNEMO Knemo   norm 0.00e+00  tree     9  depth  1
VNEMO xcnemo  norm 6.68e-01  tree 60217  depth 22
```

**It is the semilocal divergence, and nothing else.** `Unemo`, `Jnemo` and
`Knemo` are all ≤ 2185 nodes; `xcnemo` alone is 60217 against a 457-node nemo,
which is the 60185/457 disparity of §2 accounted for in full. The local term
`de/drho` is perfectly behaved (min −1.70, max ~0). `vtau` is perfectly behaved
(max 0.07). The **flux** is perfectly behaved: norm 0.0815, pointwise ±0.065.
Its **divergence** is ±8·10⁴.

The next iteration takes `local(de/drho)` to tree 51473 and the run dies in
`new failed nd=3 type=3 size=1000`.

---

## 4. Established: the amplification is set by the refinement depth

The same probes on the two paths that survive this calculation:

| variant | xc_args depth (per rung) | `div(flux)` min … max | `xc_pot` tree | outcome |
|---|---|---|---|---|
| moldft TPSS | 10 / 10–11 / 12–13 | −2.3 … +17 (1e-4, settled) | 4809–21265 | converges |
| nemo PBE (no τ) | – / 14 / 17 | +132 (d14) → +1596 (d17) | 7369–22409 | converges, 32 s |
| **nemo TPSS** | 20 / 21 / 22 | ±684 (d20) → ±1.9e4 (d21) → ±8.2e4 (d22) | 8841 → 17481 → 61121 | **runaway** |
| nemo TPSS, τ without U1 (§5) | 11 / 14 / 17 | ±52 (d11) → ±596 (d14) → ±7.8e3 (d17) | 6281 → 17481 → 79041 | reaches 1e-8, still exhausts 12 GB |

Three things to take from this table.

**The flux is the same object on both paths.** moldft: norm 0.0816, pointwise
±0.056. nemo: norm 0.0815, pointwise ±0.065. Same physics, same magnitude,
same everything — except that moldft represents it on a depth-10 tree and nemo
on a depth-22 one. moldft is not more correct here; it is less refined, and it
gets away with it.

**The extra depth carries no information.** For nemo PBE the divergence's norm
is **0.198284 at depth 14 and 0.198286 at depth 17** — identical to six digits
— while its pointwise max grows from 132 to 1596, a factor 12. A real feature
would move the norm. This is the resolution of a measure-zero structure, and
resolving it further only makes the excursion bigger.

**Removing τ's depth removes most of the excursion.** Dropping the two U1 terms
from `set_tau` (`MAD_XC_TAU_NO_U1`, a diagnostic that gives a wrong τ) takes τ
from depth 20–22 to 11/14/17 — exactly nemo PBE's depths — and the divergence's
excursion falls with it, by a factor 10 at the 1e-8 rung. That closes the chain
end to end:

```
U1's per-axis products (x/r non-smooth componentwise, |U1|^2 not)  §5
  -> tau at depth 20-22
  -> refine_to_common_level lifts every xc_arg to that depth
  -> div() of a depth-22 flux gives +-8e4 pointwise
  -> xc_pot tree 61121, xcnemo 60217 against a 457-node nemo
  -> BSH -> truncate -> larger orbital -> larger tau -> repeat
```

### What the deep-level content actually is

Not sub-threshold noise that truncation forgot. Truncating the flux at 1×, 10×
and 100× `thresh` immediately before `div()` changes essentially nothing —
tree 61121 → 60097 → 59705 → 58617, excursion 8.2e4 → 7.3e4 → 7.3e4 → 7.3e4.
The content is genuinely there in the represented function.

What is there is a **jump discontinuity at the nucleus**. `zeta = grad log(rho)`
tends to `-2Z r̂` as `r -> 0`, and the Cartesian components of `r̂` flip sign
across the origin, so each flux component `2 (de/dsigma) rho zeta_x` has a jump
of order `4 Z rho(0) de/dsigma` at `r = 0`. The divergence of a jump is
delta-like: its resolved pointwise amplitude grows without bound as the mesh
refines, which is exactly the measured behaviour. Whether one calls that
"amplified noise" or "an unphysical functional singularity being resolved" does
not matter — both have the same cure, which is not to resolve it.

### This is the hazard `XC-IMPLEMENTATION-NOTES.md` §4.2 already names

That section quotes Balbas/Martins/Soler verbatim: *"Since |grad rho| has cusps
at extrema of rho, grad g is discontinuous at these points, what causes problems
for its numerical representation"*, and adds that the sigma convention "hides
this but does not remove the underlying cusp". That is the measurement above,
arrived at from the other end. What this investigation adds is the quantity: in
an adaptive basis the problem is not that the representation is *bad*, it is that
it gets **worse the more you refine**, without bound, and that the refinement
level is set by a completely unrelated part of the code (tau, via
`refine_to_common_level`).

Two corrections to those notes follow from it.

**Section 4.3's table mislocates the hazard.** The row for form (a) reads
"noise: X ~ rho^{-4/3} grad rho in the tail, then differentiated". The tail is
not where it bites. `munge` zeroes the flux below `rhotol` and the flux carries a
factor rho above it, so the asymptotic region is quiet -- measured, the flux is
pointwise +-0.065 everywhere. The excursion is at the **nucleus**. That
mislocation is what sent this investigation at 7.1 first, which targets the tail
and changes nothing.

**Section 4.4's identity is right but its justification is not.** The remedy

    v_xc psi = (df/drho) psi - div(X psi) + X . grad(psi)

is offered there because "since psi decays, X psi is far better conditioned in
the asymptotic region". True, and irrelevant: at the nucleus psi is O(1), so
`div(X psi)` inherits the same jump-derived excursion that `div(X)` has --- the
same trade 7.1 makes with rho in place of psi, and 7.1 was measured to change the
excursion by nothing. Evaluated as written, this identity buys nothing.

It is nevertheless the *correct decomposition*, and §8 uses exactly it. What
makes it pay is not the tail but the observation that the Green's function
commutes with the divergence, so `div(X psi)` never has to be formed at all.

Two things in section 4 that this work *confirms*, incidentally. The rule "use
the same derivative operator for E_xc and v_xc" (4.3) was being violated before
`div_dft_deriv` (section 11): `dft_deriv` selected the operator for
`grad(log rho)`, which feeds both the energy and the potential, while `div()`
hardcoded ABGV. And the ban on expanding `div X` analytically (4.2, because it
drags in the Hessian) is not violated by 7.1 -- the rho-factored form keeps
exactly one divergence and puts the remainder in a pointwise term.

---

## 5. Established: τ's depth, and where it comes from

At the same threshold and for the same physics, nemo's τ is ten levels deeper
than moldft's while carrying the same information:

| | nemo | moldft |
|---|---|---|
| `rho_a` | tree 649, depth 11 | tree 521, depth 9 |
| `zeta_ax` | tree 1161, depth 11 | tree 905, depth 9 |
| `tau_a` | tree 2377, depth **20** | tree 969, depth **10** |
| `tau_a` norm | 6.84702e-01 | 6.84740e-01 |

Norms agree to four digits, so the extra depth carries nothing.

Localized to one step. Within `set_tau`, `tau after |dF|^2` is depth 17 and
`tau after cross` is depth 22 — the U1 cross term is what adds it. And the reason
is a property of U1 worth recording, because it is not obvious:

```
U1_x       norm 4.07e-01   tree 2889   depth 18
U1dot      norm 5.93e-01   tree 2441   depth 11        (= |U1|^2)
R_square   norm 1.00e+03   tree 2441   depth 11
F          norm 9.30e-01   tree  457   depth  8
dF_x       norm 7.64e-01   tree  969   depth  9
U1_x * F   norm 2.06e-01   tree 3401   depth 18        (inherits it)
```

**The Cartesian components of U1 need depth 18; their square needs 11.** U1 is
proportional to `Sr_div_S` times `smoothed_unitvec`, and the components of
`x/r` are individually non-smooth at the origin — direction is undefined at
r = 0 — whereas `|x/r|^2 = 1` is trivially smooth. Any product that uses a
single component inherits depth 18; anything that uses only the contraction or
the square does not.

`refine_to_common_level` then promotes every `xc_arg` to τ's depth. Measured:
`rho_a` goes from tree 1033 to 8265 at depth 22, an 8x inflation for nothing, and
`vtau` inherits the same. §4 shows that this promotion is the mechanism, not
merely a memory tax.

**τ also goes negative**, which it cannot: min -9.67e-08 against a peak of 23.4
in the bracket, -4.37e-08 in the returned τ. The excursion is absolute noise set
by the near-nucleus peak, so it lands wherever the true τ is smaller than that,
i.e. in the outer region. `make_libxc_args` floors those points to `tautol`,
which replaces an outer shell with a floor value. Whether that matters is
untested; the sign violation itself is real and is the same defect class as
`XC-IMPLEMENTATION-NOTES.md` §6.9.

---

## 6. The trap, twice over

**A node count measured on a late pass tells you nothing about causation.** Once
the orbital has grown, *every* object downstream of it grows in proportion, at
fixed threshold and fixed depth. A ratio like "24795 in, 61217 out" then looks
like amplification when it is inheritance. Only the **first pass at a protocol
rung** is informative. Three of the refutations in §9 (1, 9, and the first
version of §5.4 in the assessment) are corrections of conclusions drawn from
late-pass data.

**A node count is the wrong observable anyway.** The first round measured trees
and depths only, and on that evidence refutation 9 concluded that `div` is not
an amplifier because on first passes it *shrinks* its input (10011 -> 8841
nodes). That is true and it is beside the point: the tree barely moves while the
pointwise range moves by four orders of magnitude. Nothing in §4 is visible
without min/max. `Function` carries no min/max, so this needs a `unary_op`
functor with a mutex-guarded shared accumulator — that is what
`MAD_XC_PROBE=2` does, and it is the single highest-value probe in this
investigation.

Related: the growth is a **single-step jump at fixed depth with constant norm** —
`rho_a` 1033 -> 29753 nodes and `tau_a` 6537 -> 50281 in one iteration, depth
17 and 22 respectively before and after, norms constant to six digits. So it is
lateral fill-in, not deepening, and it represents no new information. "Over-
refinement in depth" was the wrong frame for the *growth*; depth is the wrong
frame for the growth and the right frame for the *amplification*.

---

## 7. Reformulations tried and measured, none of them a cure

All three are exact identities or strict smoothings, all three reproduce the
converged energy (−2.9096636x against the moldft reference −2.909663718), and
all three leave the runaway in place. Each sits behind an env switch, off by
default, so they can be A/B'd on one binary. The one that works is in §8; what
separates it from these is that it stops differentiating the flux altogether,
rather than rearranging what the derivative acts on.

### 7.1 `MAD_XC_FACTORED_GGA` — factor ρ out of the flux, not just the gradient

The existing `zeta` representation already uses `grad rho = rho grad zeta`
*inside* the flux. Push it one level up: with `A = rho B`, `B = 2 (de/dsigma) zeta`,

```
div A = rho div B + grad(rho).B = rho div B + rho (2 de/dsigma) |zeta|^2
```

The second term is purely pointwise — `|zeta|^2` is already contracted
pointwise in `make_libxc_args` — so it needs no derivative at all and folds
into the local `de/drho` term for free. Only `div B` is differentiated, and `B`
carries no factor `rho`, so the noise of the divergence is damped by `rho`
where the term is small instead of standing on its own.

**Refuted by measurement.** `rho*div(B)` has min/max ±5.8e3 at the 1e-4 rung —
*the same* as the unfactored ±5.8e3, iteration for iteration. The excursions
sit at the nucleus, where `rho` is O(1), so factoring `rho` out damps nothing
where it matters. It is also worse-conditioned overall: `B` has norm 86 against
the flux's 0.082, a factor 1000, and `xc_pot`'s tree goes 8841 -> 16777.

One real cost of the factored form, recorded because it is not obvious: `munge`
zeroes the density hard below `rhotol`, which is a harmless O(1e-7) jump in `A`
but an O(10^2) jump in `B`, since `de/dsigma` grows like `rho^{-1/3}` in the
tail. `XCfunctional::gga_ramp` replaces it with a C² switch over
`rhotol..ggatol`. That ramp is sound and reusable whatever happens to the rest.

### 7.2 `MAD_XC_NEMO_ZETA` — build ζ from the regularized density

nemo hands `prep_xc_args` the *physical* density `rho = R^2 rho_reg` and the
code then takes `grad(log rho)` numerically, i.e. the nuclear cusp goes under
the derivative operator — exactly what the regularization exists to avoid.
`zeta = grad log(rho_reg) - 2 U1` is exact, with `U1 = -grad(R)/R` analytic and
already cached, and it is the same decomposition `set_tau` uses.

**Refuted, and it makes things worse.** Trees do shrink slightly (3337 -> 2953,
8841 -> 8073) but the divergence's excursion grows: ±6.7e3 against the
baseline's ±684 at the first pass of the 1e-4 rung. The reason is §5: the cusp
is traded for U1's own direction discontinuity at the origin, which is the
non-smoothness that forced depth 18 in the first place. Both objects are
non-smooth at `r = 0`; only which one is not changes.

### 7.3 `MAD_XC_FLUX_TRUNC` — drop the deep levels before differentiating

Truncate the flux at `f * thresh` immediately before `div()`, on the theory that
the depth-11-to-22 coefficients are O(thresh) noise the common refinement level
added.

**Refuted.** `f = 1, 10, 100` gives tree 60097 / 59705 / 58617 against the
baseline's 61121, and excursion 7.3e4 against 8.2e4. The content is not below
100x thresh; it is genuinely represented. That is what pointed at a jump
discontinuity rather than at noise (§4).

### 7.4 Untested

`MAD_XC_SMOOTH_VTAU` replaces the hard `binary_munge` screen on `de/dtau` — a
jump on the `rho = ggatol` iso-surface, and `apply_tau_term` differentiates
exactly that product — with `gga_ramp`. It is the same defect class as the rest
of this section and it is meta-GGA-only, which is suggestive. It has not been
run, and the measured `vtau` (max 0.07 on settled passes, 0.83 and 24.3 on
early ones) does not make it look like the dominant term.

---

## 8. Cure: the weak form, with the divergence commuted through G

`MAD_XC_WEAK_GGA`. This is the one that works.

### 8.1 Why nothing multiplicative can work

`div(X)` *is* a derivative of `X`, and `X` has a jump at every nucleus (§4). No
algebraic rearrangement of a multiplicative potential avoids differentiating it;
§7.1 and §7.2 only change *which* discontinuous object goes under the operator,
and both were measured to change the excursion by nothing. The weak form is
different in kind, because there `X` appears undifferentiated:

    <phi|v|psi> = int (df/drho) phi psi + int X . grad(phi psi)

A jump is perfectly harmless under an integral against a smooth `grad(phi psi)`.
Only differentiating it is fatal.

### 8.2 The commutation, which is what makes it usable in this code

The weak form gives matrix elements, and the Fock matrix is happy with that. The
orbital update is not: BSH needs `V psi` as a *function*. The way through is that
`G` is a radial convolution, so it commutes with the gradient:

    G * (psi div X) = div(G * (X psi)) - G * (X . grad psi)

Evaluate the left-hand side via the right, and the divergence acts on
`G * (X psi)`, which is C^1 for bounded `X psi`. The jump is only ever convolved.
Concretely, BSHApply computes `G_i * (-2 V psi_i)`, and the piece of `V psi_i`
with no multiplicative representation is `-div(Y_i)`, so

    G_i * (-2 (-div Y_i)) = 2 div(G_i * Y_i)

is added to the update. Note that this is exactly the decomposition
`XC-IMPLEMENTATION-NOTES.md` §4.4 writes down --- the identity was right, it just
does not pay off until the divergence is pushed through `G` rather than evaluated
directly. Same treatment for the meta-gga term `-1/2 div(v_tau grad psi)`, whose
`v_tau` is hard-screened at `ggatol` and therefore also jumps.

The split, with nemo kets `F_i` and `W_i = v_tau (grad F_i - U1 F_i)`:

    mult_i = X . grad F_i + 1/2 U1 . W_i        -> joins V psi as usual
    Y_i    = X F_i        + 1/2 W_i             -> pushed through G

and `v_xc F_i = (df/drho) F_i + mult_i - div(Y_i)` identically. The R factors
cancel in the tau part exactly as in `apply_tau_term`, so nothing is divided by R
and the cusp stays in the analytic U1.

The Fock block has to be built separately, because what is missing from the
multiplicative part is precisely the term with no multiplicative representation:

    F_ij = <psi_i|df/drho|psi_j> + int X.grad(psi_i psi_j)
                                 + 1/2 int v_tau grad psi_i . grad psi_j

with `grad(R^2 F_i F_j) = R^2 [F_j grad F_i + F_i grad F_j - 2 U1 F_i F_j]`. It is
symmetric by construction, which the divergence form is not --- `Nemo::solve`
still carries the comment "not symmetric actually" about the old one.

### 8.3 Measured

He/TPSS, `protocol [1e-4, 1e-6, 1e-8]`, `econv 1e-8`, `dconv 1e-6`. Both runs on
the same machine with the same thread count, under a 6 GB RSS watchdog:

| | baseline (divergence form) | weak form + commutation |
|---|---|---|
| outcome | **killed** at 5.8 GB, 608 s | **converged**, 447 s, peak 5.4 GB |
| energy | never finished | **-2.90966370** (moldft ref -2.909663718) |
| `xc_pot` pointwise | min -5.26e4, max +6.69e4 | min -1.70, max ~0 |
| `xcnemo` tree | 57649, growing | 6297 -> 12001..18057, bounded |
| nemo tree | 457 -> 29825 -> ... | 457 -> 5193..8145, bounded |

The excursion is not reduced, it is never formed: `xc_pot` is `df/drho` alone and
the flux stays at its honest +-0.065. `xcnemo` and the orbital both oscillate
within a band instead of growing, which is the runaway signature of §1 gone.
Convergence is monotone from -2.90966325 on the first 1e-8 iteration.

Correctness, three independent checks --- the last two exercise the tau terms, so
the R cancellation and the 1/2 factors are covered:

| deck | divergence form | weak form | agreement |
|---|---|---|---|
| He/PBE, protocol 1e-6 | -2.89293474 | -2.89293462 | 1.2e-7 |
| He/TPSS, protocol 1e-6 | -2.90966357 | -2.90966316 | 4.1e-7 |
| He/TPSS, protocol 1e-8 | (dies) | -2.90966370 | 1.8e-8 vs moldft |

**Cost.** Three extra Green's-function applications per orbital. Measured on the
first 1e-8 iteration: `BSH apply (xc flux)` 6.89 s against 0.01 s for the
ordinary `BSH apply`, and `compute XCnemo` 5.33 s. On this case it is still a net
win because the baseline is dragging 57649-node trees through the same machinery,
but on a system where the divergence form is *not* pathological the weak form
will be slower.

**What it does not fix.** Peak memory is 5.4 GB against the baseline's 5.8 GB
kill point. The weak form makes the calculation terminate, not cheap: tau still
forces depth 22 on every xc_arg (`tau_a` and `local(de/drho)` both tree 15113 at
depth 22), so §5's depth inflation is untouched. Fixing it (§10) is now a pure
efficiency question rather than a correctness one.

**Untested.** Spin-polarized (the flux sums the same-spin and cross-spin
contributions, which is right for either case, but no open-shell deck has been
run); the surface term that the weak form drops and the divergence form keeps
(`XC-IMPLEMENTATION-NOTES.md` §4.1 --- irrelevant for He in a 50-bohr box where
rho at the walls is ~1e-40, not irrelevant in a tight cell); `ac_data` and OEP,
which both want a potential *function* and would have to keep building the
divergence-form one.

---

## 9. Refuted in the first round, with the measurement that killed each

Kept as-is; §6 explains why 1 and 9 are weaker than they read.

| # | hypothesis | how it died |
|---|---|---|
| 1 | `refine()` traffic in the τ machinery | every τ object is <= 0.05 GB against a 167 GB blowup: `U1dot` 0.018, `R_square` 0.018, `tau` 0.043, `W` 0.018 |
| 2 | τ truncated at `0.01*thresh` (= 1e-10 at the 1e-8 rung) rather than `thresh` | real but partial: depth 22 -> 21, τ 0.043 -> 0.029 GB, `rho_a` 0.054 -> 0.042. Growth unchanged (7753 -> 52377 nodes), still capped at 40 GB in 2:10 |
| 3 | the three-term expansion of the square (which cannot stay non-negative) vs squaring the field directly | both abort identically. Squaring gives τ min -6.3e-08, still negative, because `dot(G,G)` is itself projected and truncated. 38.7 GB, 2:09 |
| 4 | contracting U1 over the axis index prunes the depth each per-axis product needs | it does not. `U1_x*dF_x` tree 3081 depth 18; `dot(U1,dF)` tree 3465 depth **18**; truncated 1289 depth **18**. Truncation cuts nodes, never depth |
| 5 | `smoothed_unitvec`'s cutoff is `eprec`, used as a radius in bohr, so U1's depth ~ log2(L/eprec) | numerically seductive — log2(50/1e-4) ~ 19 matches the measured 18 — and wrong. `eprec` 1e-8 / 1e-6 / default all abort at 59.5 GB in 373/371/375 s. Identical to three digits |
| 6 | the loop closes through `apply_tau_term`, i.e. the non-multiplicative term poisons the orbital update | opening the loop (skip `xcnemo += apply_tau_term`) only delays it: 29.0 GB / 80 s closed, 29.5 GB / 122 s open. Consistent with §4: τ's role is to set the depth, not to carry the term |
| 7 | the orbitals are the source, upstream of all XC code | no: on the first pass of the rung the orbital is 457 nodes while `Vnemo` is already 60185 |
| 8 | `dft_deriv bspline`/`ble` as shipped | all abort: 18.4 GB/121 s, 19.5 GB/134 s, vs abgv 19.4 GB/69 s. Because `div()` never saw the setting — see 10 |
| 9 | `div(semilocal, true)` is the amplifier, turning a 24795-node flux into 61217 nodes | that ratio was measured on a *late* pass, and node counts are the wrong observable: `div` **is** the amplifier, in pointwise range rather than in tree size. See §4 and §6 |
| 10 | routing the divergence through `dft_deriv` (new `div_dft_deriv`, since `vmra.h`'s `div()` hardcodes plain ABGV) | helps locally, does not cure: ble 19.5 GB/79 s, bspline 18.9 GB/84 s, abgv 19.4 GB/71 s. BLE gives ~10% fewer nodes and depth 14 vs 20 on settled passes |

---

## 10. What is left

The runaway is cured (§8). What is left is validation and efficiency.

**Validation, before any of this is committed.**

1. **node26 at 20 threads.** Everything in §8.3 was measured locally at 6
   threads under an RSS watchdog, because node26 was unreachable for the whole
   session (`bra08` timing out). The comparison is internally consistent — both
   runs on the same machine — but the absolute numbers are not the ones the rest
   of this document quotes.
2. **The full qctest sweep with the switch off.** The non-weak path is
   textually identical (the new code is all inside `if (weak_xc)` branches and
   `compute_nemo_potentials`'s new out-params stay null), and seven cases passed
   after the `div_dft_deriv` change, but not after the weak-form edits.
3. **An open-shell deck.** §8 sums the same-spin and cross-spin flux into one
   vector field, which is right for either case, but nothing polarized has been
   run.
4. **`ac_data` and OEP**, which want a potential function and currently would
   silently get `df/drho` only. They need either the divergence-form potential
   rebuilt for their own use or an explicit refusal.

**Efficiency, now that correctness no longer depends on it.** τ still forces
depth 22 on every xc_arg, which is what keeps peak memory at 5.4 GB where moldft
needs a fraction of that. τ needs that depth only because `U1_x * (...)` products
appear in `set_tau`, and refutation 4 rules out fixing it by contracting after
the fact. What is left is to never represent U1 per-axis as a `Function` at all:
evaluate the whole cross term `-2 F (U1 . grad F)` pointwise, with U1 supplied by
its analytic functor, the way `make_libxc_args` already contracts ζ pointwise
rather than carrying `chi` as its own function. Then τ lands at ρ's depth and the
depth-22 tax on every intermediate disappears. On the evidence of §4's last table
that is worth roughly a factor 10 in the size of everything downstream.

Capping the refinement level of the flux, which was the leading candidate before
§8, is no longer worth doing: it addressed the same mechanism but by refusing to
resolve a feature rather than by never forming it, and it would have introduced
an arbitrary level with no principled value and a silent effect on every GGA
user.

Two other things worth measuring, neither on the critical path:

- **`||de/dtau||_inf`.** The BSH conditioning concern of
  `XC-IMPLEMENTATION-NOTES.md` §7 predicts an order-0, non-compact composition
  whose norm is set by that quantity and which refinement does not improve. The
  probe now reports it: `vtau` max is 0.07 on settled passes, 0.83 and 24.3 on
  early ones. Nothing alarming so far.
- **`MAD_XC_SMOOTH_VTAU`** (§7.4), never run.

---

## 11. What is in the tree

`div_dft_deriv` is a real fix and behaviour-neutral at the default
`dft_deriv abgv`: `vmra.h`'s `div()` hardcodes `gradient_operator`, so
`dft_deriv` reached `grad(log rho)` — the well-conditioned derivative — and
never the flux divergence, which is the badly-conditioned one. This is
`XC-MADNESS-ASSESSMENT.md` §4.2 exactly. Measured effect: BLE gives ~10% fewer
nodes and depth 14 instead of 20 on settled passes. Two traps in writing it,
both of which produce a `TENSOR ASSERTION FAILED ... use dot` from apparently
unrelated code: the `Derivative` objects must outlive the un-fenced `apply`
calls, and the per-axis results come back *reconstructed*, so they have to go
through `sum()` (which compresses) rather than `operator+=`.

The weak form (§8) is the piece meant to survive. It is `MAD_XC_WEAK_GGA` for
now, and it consists of:

- `XCOperator::is_weak_form()`, and a `make_xc_potential()` that returns
  `df/drho` alone and stashes the summed flux in `semilocal_flux`;
- `XCOperator::weak_xc_terms()` — the `mult_i` / `Y_i` split of §8.2;
- `XCOperator::weak_xc_matrix()` — the Fock block, symmetric by construction;
- `flux_bsh_term()` in `nemo.cc` — `2 div(G_i * Y_i)`, with the Green's function
  built *exactly* as `BSHApply` builds it (same `eps_in_green` clamp, same `lo`,
  same `bshtol`), or the two halves of the update belong to different operators;
- `compute_nemo_potentials()` gains two optional out-params, because in weak form
  the Fock block can no longer be read off `xcnemo`;
- the flux is rotated with the orbitals on canonicalisation — `Y_i` is linear in
  `F_i`, so the same unitary applies componentwise.

**Trap worth knowing before touching any of this.** The *scalar*
`Function::gaxpy` checks tree states and throws if they disagree; the *vector*
`gaxpy` in `vmra.h` coerces them instead. So `flux[i][axis] += Xket[i]` on a
freshly `compressed()` accumulator explodes while `Y[axis] += Xket` does not, and
the failure surfaces as `TENSOR ASSERTION FAILED ... use dot` from code that has
nothing to do with the mistake. Accumulate over vectors, not over scalars. The
same trap ate two build cycles in `div_dft_deriv`.

Everything else is investigation scaffolding, all off by default:

| switch | what it does | §
|---|---|---|
| `MAD_XC_WEAK_GGA` | **the cure**: weak form + divergence through G | 8 |
| `MAD_XC_PROBE` | 1 = tree/depth/norm census per term, 2 = also pointwise min/max | 6 |
| `MAD_XC_FACTORED_GGA` | ρ-factored semilocal potential | 7.1 |
| `MAD_XC_NEMO_ZETA` | ζ from the regularized density and analytic U1 | 7.2 |
| `MAD_XC_FLUX_TRUNC` | truncate the flux at N×thresh before `div()` | 7.3 |
| `MAD_XC_SMOOTH_VTAU` | C² ramp instead of `binary_munge` on `de/dtau` | 7.4 |
| `MAD_XC_TAU_NO_U1` | **diagnostic, wrong τ**: drop U1 from `set_tau` | 4 |

`XCfunctional::gga_ramp` (C² switch over `rhotol..ggatol`, smooth in log ρ) is
the one piece of §7 worth keeping regardless of what happens to the rest.

`ncf->square()` is still not cached: it re-projects a functor on every call,
unlike `U1vec()` which returns members, and `set_tau` runs at two nemo sites —
so four functor projections per SCF iteration (2x R^2, 2x |U1|^2). `Nemo`
already caches `R_square` as a member; `XCOperator` holds only `ncf` and cannot
reach it.

Regression status with all switches off: `nemo_he_tpss`, `nemo_he_pbe`,
`nemo_he_hf`, `scf_he_tpss`, `scf_li_tpss`, `scf_h2o_hf`, `nemo_h2o_canon` all
pass — but that sweep was run *before* the weak-form edits, and has to be
repeated (§10). Note `MAD_NUM_THREADS` must be set: unset on node26 it defaults
to 96 and `madqc` aborts with *"When configured with
MADNESS_TASK_BACKEND=Pthreads MAD_NUM_THREADS cannot exceed 64"*, which through
ctest reads as seven simultaneous deck regressions.

---

## 12. Reproduction

```
protocol [1e-4, 1e-6, 1e-8], econv 1e-8, dconv 1e-6, He, xc tpss, --wf=nemo
```

Always under `ulimit -v` — unbounded it will take the whole machine. 12 GB
reaches the 1e-8 rung and dies within ~90 s, which is enough for first-pass
instrumentation and is the setting to use. `MAD_NUM_THREADS=20`; leaving it
unset on node26 gives 96 and `madqc` aborts with *"When configured with
MADNESS_TASK_BACKEND=Pthreads MAD_NUM_THREADS cannot exceed 64"*, which through
ctest looks exactly like seven simultaneous deck regressions.

`MAD_XC_PROBE=2` gives the census of §3 directly. The pointwise min/max needs a
`unary_op` functor with a mutex-guarded shared accumulator — `Function` carries
no min/max — and per §6 it is the observable without which none of this is
visible. Run the same deck with `--wf=scf` and with `xc pbe` for the two
controls in §4's table; they cost 76 s and 32 s and they are what turn a number
into a comparison.

`ulimit -v` does not exist on Darwin (`setrlimit failed: invalid argument`), so
the local runs of §8.3 were capped by sampling RSS from a watchdog loop and
killing at 6 GB. That is what makes it safe to run the *baseline* on a
workstation at all: unbounded it takes the whole machine.

The `Nemo::solve` BSH/KAIN trace of §2 and the step-by-step census through
`prep_xc_args`/`set_tau` (with an auto-abort on the first object over a node
limit, which saved most of the run time) were temporary and are not in the tree.
Rebuild them from this description rather than guessing at new probes.
