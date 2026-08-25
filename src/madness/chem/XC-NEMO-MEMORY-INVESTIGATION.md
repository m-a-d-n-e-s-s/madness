# nemo meta-GGA memory runaway: findings

Investigation log, 2026-08-25, on branch `pr-nemo-dft` (base `d8a438201`).
Companion to `XC-MADNESS-ASSESSMENT.md` §5.4, whose first version of this
section is superseded by what follows.

**Status: not resolved.** Ten candidate causes tested, ten refuted by
measurement. What is established is a causal chain from the potential into the
orbital update; what is not established is which term makes that potential large
in the first place. The negatives are the useful part — and so is the
methodological trap in §5, which produced three of my own wrong conclusions and
will produce more if the next person does not read it first.

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
in 32 s. So τ is necessary to trigger it. That is not the same as τ being the
amplifier — see §4.

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

Everything inside the XC code downstream of that is a faithful consequence:
`grad log(arho)` 5979 -> 120299, `xcargs post-truncate` 6820 -> 149804,
`tau post-R2` 7753 -> 53481. The XC code is handed an already-grown orbital.

The same disparity is present one rung earlier and merely survivable: at thresh
1e-6, `Vnemo` is 16409-24497 nodes against a 393-881-node orbital.

**What this establishes.** The loop is: large `Vnemo` -> BSH -> truncate -> larger
orbital -> larger density and τ -> larger `Vnemo`. Each cycle multiplies. The
question the chain does *not* answer is why `Vnemo` is 60185 nodes on the *first*
pass of the rung, when the orbital is still 457 nodes and no feedback has
happened yet. That is the open question.

---

## 3. Established: τ's depth, and where it comes from

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
`vtau` inherits the same.

**τ also goes negative**, which it cannot: min -9.67e-08 against a peak of 23.4
in the bracket, -4.37e-08 in the returned τ. The excursion is absolute noise set
by the near-nucleus peak, so it lands wherever the true τ is smaller than that,
i.e. in the outer region. `make_libxc_args` floors those points to `tautol`,
which replaces an outer shell with a floor value. Whether that matters is
untested; the sign violation itself is real and is the same defect class as
`XC-IMPLEMENTATION-NOTES.md` §6.9.

---

## 4. Refuted, with the measurement that killed each

| # | hypothesis | how it died |
|---|---|---|
| 1 | `refine()` traffic in the τ machinery | every τ object is <= 0.05 GB against a 167 GB blowup: `U1dot` 0.018, `R_square` 0.018, `tau` 0.043, `W` 0.018 |
| 2 | τ truncated at `0.01*thresh` (= 1e-10 at the 1e-8 rung) rather than `thresh` | real but partial: depth 22 -> 21, τ 0.043 -> 0.029 GB, `rho_a` 0.054 -> 0.042. Growth unchanged (7753 -> 52377 nodes), still capped at 40 GB in 2:10 |
| 3 | the three-term expansion of the square (which cannot stay non-negative) vs squaring the field directly | both abort identically. Squaring gives τ min -6.3e-08, still negative, because `dot(G,G)` is itself projected and truncated. 38.7 GB, 2:09 |
| 4 | contracting U1 over the axis index prunes the depth each per-axis product needs | it does not. `U1_x*dF_x` tree 3081 depth 18; `dot(U1,dF)` tree 3465 depth **18**; truncated 1289 depth **18**. Truncation cuts nodes, never depth |
| 5 | `smoothed_unitvec`'s cutoff is `eprec`, used as a radius in bohr, so U1's depth ~ log2(L/eprec) | numerically seductive — log2(50/1e-4) ~ 19 matches the measured 18 — and wrong. `eprec` 1e-8 / 1e-6 / default all abort at 59.5 GB in 373/371/375 s. Identical to three digits |
| 6 | the loop closes through `apply_tau_term`, i.e. the non-multiplicative term poisons the orbital update | opening the loop (skip `xcnemo += apply_tau_term`) only delays it: 29.0 GB / 80 s closed, 29.5 GB / 122 s open |
| 7 | the orbitals are the source, upstream of all XC code | no: on the first pass of the rung the orbital is 457 nodes while `Vnemo` is already 60185 |
| 8 | `dft_deriv bspline`/`ble` as shipped | all abort: 18.4 GB/121 s, 19.5 GB/134 s, vs abgv 19.4 GB/69 s. Because `div()` never saw the setting — see 10 |
| 9 | `div(semilocal, true)` is the amplifier, turning a 24795-node flux into 61217 nodes | that ratio was measured on a *late* pass. On first passes `div` **shrinks** its input: abgv 10011 -> 8841, ble 8859 -> 7945. See §5 |
| 10 | routing the divergence through `dft_deriv` (new `div_dft_deriv`, since `vmra.h`'s `div()` hardcodes plain ABGV) | helps locally, does not cure: ble 19.5 GB/79 s, bspline 18.9 GB/84 s, abgv 19.4 GB/71 s. BLE gives ~10% fewer nodes and depth 14 vs 20 on settled passes |

---

## 5. The trap that produced three of these

**A node count measured on a late pass tells you nothing about causation.** Once
the orbital has grown, *every* object downstream of it grows in proportion, at
fixed threshold and fixed depth. A ratio like "24795 in, 61217 out" then looks
like amplification when it is inheritance.

Only the **first pass at a protocol rung** is informative, before any feedback
has occurred. Three of the refutations above (1, 9, and the first version of
§5.4 in the assessment) are corrections of conclusions I drew from late-pass
data.

Related: the growth is a **single-step jump at fixed depth with constant norm** —
`rho_a` 1033 -> 29753 nodes and `tau_a` 6537 -> 50281 in one iteration, depth
17 and 22 respectively before and after, norms constant to six digits. So it is
lateral fill-in, not deepening, and it represents no new information. "Over-
refinement in depth" was the wrong frame throughout.

---

## 6. Open, and how to attack it

The one question the chain does not answer: **why is `Vnemo` 60185 nodes on the
first pass of the 1e-8 rung, against a 457-node orbital?**

`Vnemo = Unemo + Jnemo - Knemo [+ pcmnemo] [+ xcnemo]` (`nemo.cc:523`), and

```
xcnemo  = truncate(xc_pot * nemo);              // multiplicative
xcnemo += xcoperator.apply_tau_term(nemo);      // non-multiplicative
```

with `xc_pot` = `de/drho - div(2 (de/dsigma) grad rho)` for restricted meta-GGA
(the cross-spin flux term is polarized-only).

**Do this first, and only on first passes:** report the tree of `Unemo`,
`Jnemo`, `Knemo` and `xcnemo` separately, plus the three constituents of
`xc_pot`, on the first `compute_nemo_potentials` call of each rung. My own
per-term numbers (`int[0]` 8265, flux 24795, `div(flux)` 61217,
`xc_pot*nemo` 59105, `tau_term` 41545, `xcnemo` 60265) were taken on a late pass
and are worthless for this question. It is entirely possible that a
non-XC term dominates and that τ's role is only to raise the common refinement
level.

Then, in rough order of expected value:

1. **Cap the potential's tree before it multiplies the orbitals.** A 60185-node
   potential acting on a 457-node orbital is wrong on its face, whatever the
   term. `truncate` at `thresh` demonstrably does not do it (78745 -> 29825, still
   65x). This is the most direct interruption of the loop and does not require
   knowing which term is guilty.
2. **Keep U1 out of per-axis products.** Only the contraction and `|U1|^2` are
   shallow. This is a real inefficiency (τ depth 22 vs moldft's 10) even though
   it is not the runaway, and it also affects `apply_tau_term`, where
   `mul(world, U1[axis], ...)` appears twice.
3. **Measure `||de/dtau||_inf`.** The BSH conditioning concern of
   `XC-IMPLEMENTATION-NOTES.md` §7 predicts an order-0, non-compact composition
   whose norm is set by that quantity and which refinement does not improve. It
   has never been measured on either path, and a stalling-but-not-diverging SCF
   at the tightest rung is exactly its predicted signature. nemo does not set the
   operator's print level, so it currently logs nothing.

---

## 7. Two fixes found on the way, neither a cure

Both are real, both are independent of the runaway, and neither is committed.

**`div_dft_deriv`.** `vmra.h`'s `div()` hardcodes
`gradient_operator<T,NDIM>(world)` with no variant parameter, so `dft_deriv`
reached `grad(log rho)` — the well-conditioned derivative — and never the flux
divergence, which is the badly-conditioned one. This is
`XC-MADNESS-ASSESSMENT.md` §4.2 exactly, and it is a five-line fix: build the
three derivative operators with `set_bspline1()`/`set_ble1()` per `dft_deriv` and
sum `apply(D_axis, v[axis])`, as `set_tau` already does. Measured effect: BLE
gives ~10% fewer nodes and depth 14 instead of 20 on settled passes. It changes
the GGA path for every user, so it needs its own qctest validation.

**`ncf->square()` is not cached.** It re-projects a functor on every call, unlike
`U1vec()` which returns members, and `set_tau` runs at two nemo sites — so four
functor projections per SCF iteration (2x R^2, 2x |U1|^2). `Nemo` already caches
`R_square` as a member; `XCOperator` holds only `ncf` and cannot reach it.

---

## 8. Reproduction

```
protocol [1e-4, 1e-6, 1e-8], econv 1e-8, dconv 1e-6, He, xc tpss, --wf=nemo
```

Always under `ulimit -v` — unbounded it will take the whole machine. 12 GB
reaches the 1e-8 rung and dies within ~90 s, which is enough for first-pass
instrumentation and is the setting to use.

The instrumentation used here was temporary and is not in the tree. It was four
env-var-guarded probes: a per-`xc_arg` census before `refine_to_common_level`
(norm, pointwise min/max, `tree_size()`, `max_depth()`, `get_size()`); a
step-by-step census through `prep_xc_args` and `set_tau` with an auto-abort on
the first object over a node limit; a per-term census inside
`make_xc_potential`; and a trace of `Nemo::solve`'s BSH/KAIN block. The
pointwise min/max needs a `unary_op` functor with a mutex-guarded accumulator —
`Function` carries no min/max — and that is what makes the sign violations
visible. Rebuild them from this description rather than guessing at new probes;
the auto-abort in particular saved most of the run time.
