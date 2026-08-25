# MADNESS XC implementation: assessment

Companion to `XC-IMPLEMENTATION-NOTES.md`, which collects the literature working
equations. This document measures MADNESS against them, with an emphasis on
derivative stability and on the numerical hazards specific to a real-space code.
Response is out of scope.

Originally assessed at `399acaf9e` (`pr-pass-libxc-syntax`); **revised
2026-08-24 against `pr-implement-mgga`**, where ground-state meta-GGA landed and
most of the original defect list was closed. Files:
`src/madness/chem/xcfunctional*.{h,cc}`, `SCFOperators.{h,cc}`, with
`nemo.{h,cc}`, `correlationfactor.h`, `SCF.cc`, `oep.h` and the tests.

Every claim carries a `file:line` or a measurement. Literature citations are in
the companion document; **[M]** marks something measured in this repo.

---

## 1. Architecture as built

```
SCF::apply_potential                                        SCF.cc
  └─ XCOperator ctor  (fresh, every iteration, per spin)
       └─ prep_xc_args(arho, brho)                          SCFOperators.cc:758
            ├─ logdensa = unary_op(arho, logme())
            ├─ grada = grad{,_bspline_one,_ble_one}(logdensa)
            └─ store rhoa[, rhob], zeta_{a,b}{x,y,z}
  ├─ set_tau(amo, aocc, bmo, bocc)      (meta-GGA only)     SCFOperators.cc:477
  ├─ compute_xc_energy()
  └─ vloc += make_xc_potential()
       ├─ refine_to_common_level(xc_args)
       ├─ multi_to_multi_op_values(xc_potential)
       ├─ dft_pot  = intermediates[0]
       ├─ dft_pot -= div(intermediates[1..3])       same-spin flux
       ├─ dft_pot -= div(intermediates[4..6])       cross-spin flux (polarized)
       └─ vtau     = intermediates[7 | 4]           (meta-GGA only)
  then  Vpsi = vloc*psi  +  apply_tau_term(psi)             (meta-GGA only)
```

`make_libxc_args` (`xcfunctional_libxc.cc:284`) reconstructs libxc's inputs
pointwise:

$$
\rho=\mathrm{munge}(2\rho_\alpha),\qquad
\nabla\rho_s=\rho_s\zeta_s,\qquad
\sigma_{st}=\rho_s\rho_t\,(\zeta_s\!\cdot\!\zeta_t),\qquad
\tau_s\ \ge\ \rho_s(\zeta_s\!\cdot\!\zeta_s)/8 ,
$$

and `vxc` returns the flux vector componentwise. The non-multiplicative term is
$-\tfrac12\sum_x D_x(v_\tau D_x\psi_i)$ (`SCFOperators.cc:548`), applied to the
orbitals rather than folded into `vloc`.

---

## 2. What the implementation gets right

### 2.1 The ζ parametrization is better than the literature default

MADNESS never stores $\nabla\rho$ as a `Function`. It stores
$\zeta_\sigma=\nabla\ln\rho_\sigma$ and reconstructs
$\nabla\rho_\sigma=\rho_\sigma\zeta_\sigma$ pointwise, applying the gradient to
`log(rho)` directly rather than dividing $\nabla\rho$ by $\rho$.

This is a genuinely good answer to the tail hazard of
`XC-IMPLEMENTATION-NOTES.md` §6.5–6.7. Where $\rho\sim e^{-2\kappa r}$,
$\zeta\to-2\kappa\hat r$ — **bounded** — while $\nabla\rho\to0$ and
$\partial f/\partial\sigma\sim\rho^{-4/3}\to\infty$. The delicate product of a
vanishing and a diverging quantity is replaced by a bounded one. Structural
advantage over grid codes that differentiate $\rho$.

> **⚠ The earlier version of this section claimed more than the code delivered,
> and the gap was a wrong-answer bug.** It said that
> $\sigma_{ss}=\rho_s^2|\zeta_s|^2\ge0$ holds "identically", that Cauchy–Schwarz
> holds "automatically", and that "MADNESS's construction cannot violate it".
> Those statements are true of the *algebra* and were false of the *code*,
> because $\chi_{st}=\zeta_s\!\cdot\!\zeta_t$ was carried as its own
> multiwavelet function while $\nabla\rho_s$ was built from the $\zeta$
> components. Two independently projected representations of one product are not
> pointwise equal, and near the nuclear cusp they differ by O(1). See §4.5.
>
> **The advantage is real, but it is conditional on contracting $\chi$ pointwise
> from $\zeta$ at the quadrature point** — which is what the code now does
> (`xcfunctional_libxc.cc:272`, `chi_of()`). Then all three Gram relations hold
> exactly and the positivity floors become no-ops. A general statement of the
> principle is in `XC-IMPLEMENTATION-NOTES.md` §6.9.

### 2.2 The GGA potential and its factors

`make_xc_potential` builds the multiplicative divergence form with same-spin and
cross-spin terms as separate `div` applications, citing Yanai 2005 Eq. (12) by
name. Correct choice for this code: form (b) would forfeit the multiplicative
potential the BSH iteration needs, and form (c) (White–Bird) is ill-defined under
adaptive refinement (companion §4.3).

The factors are right: `2.0*vs[2*ispin]*ddens[ispin]` same-spin, and
`vs[1]*ddens[1-ispin]` cross-spin — factor 1, multiplying the *opposite* spin's
gradient. The cross-spin term, the one most often omitted, is present in both the
functional and the operator.

### 2.3 The response-kernel convention is self-consistent

The companion §8.1 flags a factor-of-2 hazard: the $\{2,2,2,4\}$ pattern is
correct only under $\sigma_{\rm pt}=\nabla\rho\!\cdot\!\nabla\rho_{\rm pt}$
without the 2. MADNESS stores
`enum_sigma_pta_div_rho`$=\zeta_\alpha\!\cdot\!\nabla\rho_{\rm pt}$ and multiplies
by the full density, giving exactly
$\sigma_{\rm pt}=\nabla\rho\!\cdot\!\nabla\rho_{\rm pt}$ — no factor 2 — and
`fxc_apply`'s coefficients are `{2,2,4,2}`, matching. **Correct**, and worth
recording explicitly because the pairing is fragile and undocumented where it
matters.

### 2.4 meta-GGA is implemented the way the literature says to

- $\tau$ built from orbital gradients (route A), occupation-weighted, with the
  libxc $\tfrac12$ (`SCFOperators.cc:477`). Virtuals carry occupation zero and
  are skipped; fractional occupations are honoured.
- The non-multiplicative term uses the **nested** form
  $-\tfrac12\sum_xD_x(v_\tau D_x\psi)$, so $v_\tau$ is only ever multiplied,
  never differentiated — the single most important trick (companion §6.2) — and
  the form is self-adjoint by construction rather than up to discretization
  error.
- It is applied to the orbitals alongside exact exchange, not folded into `vloc`.
- TPSS first, with r²SCAN and SCAN as aliases; the deck comments record why that
  order (companion §6.3).
- The ellipticity diagnostic $1+\min v_\tau>0$ is logged at `print_level>=2`
  (moldft only — nemo does not set the operator's print level).
- It works on the regularized path too: with a nuclear correlation factor both
  $\tau$ and the non-multiplicative term are built by the $\psi=RF$ product rule,
  with the cusp confined to the analytic $\mathbf U_1$ and the $R$ factors of the
  operator cancelling identically. See §4.4.
- The tau floor is $\tau_s\ge\rho_s\chi_{ss}/8$ — the von Weizsäcker bound as a
  *product*, so no division by the density and no contact with $\sigma$'s
  positivity floor, which fed into $\sigma/(8\rho)$ would grow like $1/\rho$
  exactly where the true $\tau_W$ vanishes.

### 2.5 Smoothed derivative operators exist; nemo has an asymptotic correction

`dft_deriv` selects ABGV (default), B-spline or BLE, and is now plumbed into all
six `XCOperator` constructors including the nemo ones. `nemo.cc` applies the
Tozer–Handy correction from `AC.h` to the finished XC potential — the requirement
of companion §6.8, since the BSH kernel needs bound orbitals. **The plain moldft
path has no equivalent.**

---

## 3. Defects

### 3.1 Closed since the original assessment

| # | Defect | Consequence then | Closed by |
|---|---|---|---|
| 3.1 | `∇ρ_β` built from `enum_zetaa_*` | wrong potential, every open-shell GGA | `zetab_x = xc_args[enum_zetab_x]` |
| 3.2 | `σ_ab` clamped to `max(1e-14, …)` | sign of a bilinear quantity destroyed | unfloored; CS bound re-imposed after the diagonal floors |
| 3.4 | `prep_xc_args` truncated the member, not the local | `extra_truncation` was a no-op everywhere | `truncate(world, xcargs, extra_truncation)` |
| 3.5 | `dft_deriv` absent from the nemo ctors | `dft_deriv bspline` silently ignored in nemo | plumbed into all ctors |
| 3.6 | `hf_coeff` only for hardcoded aliases | `xc HYB_GGA_XC_B3LYP` ran with **zero** exact exchange | `xc_hyb_exx_coef` query |
| 3.7 | `XC_FAMILY_HYB_MGGA` missing | TPSSh/M06-2X threw from the wrong place | added throughout the dispatch |
| 3.9 | `test_dft.cc` dead kernel assertion (stale `err`) | the kernel check could not fail | `err` recomputed |
| — | `fxc_apply` weighted the accumulated result | wrong for any non-unit multi-functional weight except by luck of ordering | weight applied per contribution |
| §6 | no meta-GGA at all | — | implemented; see §2.4 |

**[M]** The two wrong-answer bugs (3.1, 3.2) plus the Gram-matrix defect of §4.5
were all invisible to the test suite for the same structural reason: the only
polarized quantity exercised was the *energy*, fed a spin-symmetric density
(`arho, arho`), where $\zeta_\alpha\equiv\zeta_\beta$ and
$\sigma_{\alpha\beta}>0$. A perfect blind spot. That is now covered end-to-end by
the `scf_li_tpss` qctest deck (open shell, unequal spin densities, asserting on
energies *and* eigenvalues) and at tensor level by
`test_meta_gga_dedtau_polarized`. The **unit-level** gap remains: see §3.3.

### 3.2 Still live

- **`has_fxc()` and `has_kxc()` return `false` unconditionally**
  (`xcfunctional_libxc.cc:235`, `:240`, comment "not thought about this yet")
  while `fxc_apply` is fully implemented and used. Called nowhere, so
  dead-but-wrong.
- **`is_dft()` has divergent semantics between build configurations:**
  `funcs.size()>0` with libxc (`:230`, with the alternative commented out) versus
  `is_lda()||is_gga()||is_meta()` without (`xcfunctional_ldaonly.cc:77`).
- **`fxc_apply` calls `xc_gga_fxc` and `xc_gga_vxc` back-to-back** where
  `xc_gga_vxc_fxc` exists — same duplication as §5.2.
- `expme` (`SCFOperators.h:910`) is defined and never used.
- `XCOperator`'s matrix-element overloads throw, closing off the weak-form route
  entirely.
- `SCF::apply_potential` uses the density-**recomputing** ctor although
  `arho`/`brho` are in scope at the call site and a ctor taking them exists.
- `SCF.cc` carries its own open question in a comment: `??RJH?? Won't this
  incorrectly exclude hybrid DFT with coeff=1.0?`
- **Slots 10–13 of `xc_arg`** (`enum_saa/sab/sbb/sigtot`) are written by nothing
  and read by nothing; `XCfunctional::plot()` fills `enum_saa`, which
  `make_libxc_args` has never read, then dereferences a null pointer — broken for
  anything above LDA, before and after this work, and unreferenced. Delete it or
  fix it deliberately, and sweep the enum against its actual consumers.

### 3.3 Test coverage

| Test | Label | Covers |
|---|---|---|
| `test_dft.cc` | **short** | `LDA_X`, hybrid coefficients, meta-GGA identities, polarized `de/dtau` (tensor level) |
| `test_SCFOperators.cc` | **long** | closed-shell `lda`, `LDA_X`, `pbe`, `bp` against hardwired references |
| `testxc.cc` | **built, but not a test** | — |
| `scf_he_pbe0`, `scf_he_tpss`, `scf_li_tpss` (qctest) | **medium** | closed-shell hybrid + meta-GGA, open-shell meta-GGA |
| `nemo_he_pbe`, `nemo_he_tpss` (qctest) | **medium** | GGA and meta-GGA on the regularized path, each cross-checked against its moldft counterpart |

Both C++ tests are gated on `TARGET Libxc::xc`.

Remaining gaps:

- **No *reference-valued* polarized test, and no `nspin=1` vs `nspin=2`
  consistency test.** `test_XCOperator` does assert on the polarized
  `Function`-level potential, and on a genuinely spin-asymmetric density: it
  checks the exact symmetry $v_\alpha[\rho_a,\rho_b]=v_\beta[\rho_b,\rho_a]$
  pointwise, needing no reference values. What is missing is a polarized case
  pinned to a reference number, and the companion §10 test — the same
  closed-shell system through `nspin=1` and `nspin=2`, agreeing to roundoff on
  energy, potential *and* kernel. Note two things that test provably cannot
  catch: a Gram-matrix violation (companion §6.9), and anything symmetric under
  the spin swap, which is why §4.5's defect survived `test_XCOperator` even
  though that test asserts on exactly the right object. Its densities are also
  smooth Gaussians, so the O(1) cusp projection error never arises there at all.
- `testxc.cc` is not in `CHEM_TEST_SOURCES_{SHORT,LONG}` but *is* built, as
  `CHEM_OTHER_TESTS` via `add_mad_executable` with the comment "not included in
  the unit tests". So it compiles on every build and produces no ctest entry —
  a manually-run executable, like `plotxc`. `testlda.cc` and
  `test_exchangeoperator.cc` are genuinely absent from `CMakeLists.txt`.

---

## 4. Numerical stability

The substance of this assessment. Items 4.1–4.4 are unchanged from the original;
4.5 is the case study that came out of the open-shell meta-GGA divergence and
supersedes a claim in §2.1.

### 4.1 The flux vector is differentiated directly

`div(semilocal, true)` differentiates

$$
\mathbf X_\sigma=2\frac{\partial f}{\partial\sigma_{\sigma\sigma}}\nabla\rho_\sigma+\frac{\partial f}{\partial\sigma_{\alpha\beta}}\nabla\rho_{\sigma'}
$$

itself. The companion §4.4 identity,
$-\psi\nabla\!\cdot\!\mathbf X=-\nabla\!\cdot\!(\mathbf X\psi)+\mathbf X\!\cdot\!\nabla\psi$,
would move the differentiation onto the better-conditioned $\mathbf X\psi$ — and
cannot be used, because `make_xc_potential()` must return an orbital-independent
`real_function_3d` for `vloc`.

The ζ parametrization substantially mitigates this ($\mathbf X$ is built from a
bounded $\zeta$ and a munged $\rho$ rather than from a vanishing $\nabla\rho$
times a diverging $\partial f/\partial\sigma$), but $\mathbf X$ still contains
$\partial f/\partial\sigma$ with its $\rho^{-4/3}$ growth, and that is the object
differentiated.

### 4.2 The smoothing is applied to the easy derivative, not the hard one

**The most consequential stability finding, and still unaddressed.**

- `grad log(rho)` — the **well-conditioned** derivative, of a function smooth
  except at the nuclear cusp and the munge isosurface — uses the user-selectable
  smoothed operator.
- `div(semilocal, true)` — the **badly-conditioned** derivative, of a function
  containing $\partial f/\partial\sigma\sim\rho^{-4/3}$ — always uses plain ABGV.
  `div` at `vmra.h` hardcodes `gradient_operator<T,NDIM>(world)` with no variant
  parameter.

So `dft_deriv` has **no effect whatsoever on the divergence step**, and whatever
noise suppression B-spline/BLE provide is spent on the derivative that needed it
least. Given the SCAN/r²SCAN finding that grid sensitivity in $\tau$-dependent
functionals traces specifically to oscillations in
$\nabla(\partial e/\partial\tau)$ (companion §6.3), the asymmetry matters more now
that meta-GGA has landed, not less.

`apply_tau_term` does honour `dft_deriv`, via its own `make_derivative(axis)`
(`SCFOperators.cc:602`) rather than the free `div` — so the meta-GGA term is
already doing the right thing and the GGA flux divergence is the outlier. A
variant-selectable `div` (or the same local construction) would close it.

### 4.3 The munge is a hard clamp, and it sits under a derivative

`logme` (`SCFOperators.h:894`): `log(std::max(1.e-14,val))+14.0`. The `max` makes
$\log\rho$ exactly constant below $\rho=10^{-14}$, so $\zeta\equiv0$ there and the
GGA correctly degrades to LDA-like behaviour in the far tail rather than
diverging. Sound, and the `+14` shift is inert under differentiation.

But the clamp leaves $\log\rho$ **$C^0$ and not $C^1$** across the isosurface
$\rho=10^{-14}$, and a derivative operator is applied straight across that kink.
In an adaptive basis a kink on a large-area surface is expensive (refinement
chases it) and noisy — companion §6.7 item 1.

The codebase *has* a smooth alternative and does not use it:
`XCfunctional::polyn` (`xcfunctional.h`) is the unique quintic joining a constant
to the identity with matched first and second derivatives, exposed as
`munge_old` — and referenced nowhere. The live `munge` is the hard cutoff
`if (rho <= rhotol) rho=rhomin;`. Switching `logme` to a $C^2$ ramp is a small,
self-contained experiment with a clear predicted benefit.

Live defaults: `rhotol = 1e-7`, `rhomin = 0.0`, `ggatol = 1e-4`,
`tautol = 1e-12`, all overridable from the input line
(`RHOTOL`/`RHOMIN`/`GGATOL`/`TAUTOL`). Note **three** screening mechanisms engage
on different isosurfaces: `logme`'s $10^{-14}$, `munge`'s `rhotol` $=10^{-7}$, and
the absolute $10^{-14}$ floor on $\sigma$ in `make_libxc_args` (unrelated to
either and not user-adjustable). `rhomin = 0.0` means `munge` maps the tail to
*exactly zero*, which is what makes the $\sigma$ floor the only thing standing
between libxc and a zero density — see §4.5.

### 4.4 The nuclear cusp: half fixed, half still open

**Status.** The $\tau$ half is implemented (`nemo: meta-gga on the regularized
path`). The $\zeta$/$\sigma$ half — the GGA reduced gradient — is still open, and
is now the clearest actionable finding here.

`XCOperator`'s nemo constructor forms the **full physical density**, cusp
included, and hands it to the same code path as moldft:
`arho = (arhonemo * nemo->R_square).truncate(...)`. `prep_xc_args` then takes
`log(arho)` and differentiates it numerically. So the $R^2$ multiplication
reintroduces the cusp and a numerical gradient is taken across it — discarding
the entire point of the nemo factorization for this quantity.

Meanwhile `XCOperator` holds `ncf`, documented verbatim as "the nuclear
correlation factor, if it exists, **for computing derivatives for GGA**". It is
now read by `set_tau` and `apply_tau_term`, which use it to apply the product rule
below — but *not* by `prep_xc_args`, which is the one place the comment actually
names. The GGA path still differentiates `log(arho)` numerically across the cusp.

**The fix is an identity, not an approximation.** With $\rho=R^2\rho_{\rm nemo}$
and $\mathbf U_1=-\nabla R/R$ (`correlationfactor.h`),

$$
\boxed{\;\zeta=\nabla\ln\rho=-2\mathbf U_1+\nabla\ln\rho_{\rm nemo}\;}
$$

The cusp lives entirely in the first term, which is **precomputed, stored and
analytic**, built on a `smoothed_unitvec` so it is well-behaved at $r=0$. Only
$\nabla\ln\rho_{\rm nemo}$ — cusp-free by construction — is differentiated
numerically.

Why $F$ is cusp-free is what makes this exact rather than approximate: every NCF
in this family satisfies $S_A'(0)/S_A(0)=-Z_A$ exactly — for the default Slater
factor $S_A=\frac{1}{a-1}e^{-aZ_Ar}+1$, `Sr_div_S` gives $-aZ/a=-Z$ at $r=0$,
**independent of $a$** — so $R$ carries precisely the same leading cusp as $\psi$
and $F=\psi/R$ has a continuous gradient at the nucleus.

**And the analytic version was already written, then never called.**
`Nemo::make_sigma` (`nemo.h`, `nemo.cc`) computes the GGA reduced gradient fully
analytically, with exactly this product-rule decomposition, per its own doc
comment. It is **dead code** — `grep -rn make_sigma src/` finds only its
declaration and definition, and `git log -S make_sigma` shows it untouched since
the `chem` directory move. Its sibling `Nemo::make_ddensity` is live, but only for
the nuclear Hessian. Reviving `make_sigma` is smaller than writing the
decomposition fresh, but check it against the current ζ interface rather than
assuming it correct after years without a caller.

**The same applies to $\tau$, and that half is now done.** With $\psi_i=RF_i$,
$\nabla\psi_i=R(\nabla F_i-\mathbf U_1F_i)$, so

$$
\boxed{\;\tau_\sigma=\tfrac12R^2\sum_iw_i\bigl|\nabla F_i-\mathbf U_1F_i\bigr|^2\;}
$$

which is the decomposition `OEP::compute_total_kinetic_density` (`oep.h`) already
used — and which already honoured `dft_deriv`, unlike the GGA flux path (§4.2).
`set_tau` now applies it directly, branching on its own `ncf`, so callers just
pass what they have.

The non-multiplicative term needed a second identity. nemo's equations are for
$F$, so what `apply_tau_term` must return is
$R^{-1}[-\tfrac12\nabla\!\cdot\!(v_\tau\nabla(RF))]$. With
$\mathbf W=v_\tau(\nabla F-\mathbf U_1F)$, so that $v_\tau\nabla(RF)=R\mathbf W$,

$$
\nabla\!\cdot\!(R\mathbf W)=R\,\nabla\!\cdot\!\mathbf W+\nabla R\!\cdot\!\mathbf W
 = R\bigl(\nabla\!\cdot\!\mathbf W-\mathbf U_1\!\cdot\!\mathbf W\bigr)
$$

so the $R$ factors cancel identically, leaving
$-\tfrac12(\nabla\!\cdot\!\mathbf W-\mathbf U_1\!\cdot\!\mathbf W)$. **Nothing is
divided by $R$**, $v_\tau$ is still only multiplied, and with no ncf $\mathbf U_1$
drops out and it reduces to the plain nested form — so the moldft path is
untouched. Wired into `compute_nemo_potentials` and
`compute_energy_regularized`; `Nemo::make_fock_operator` refuses a meta-GGA
instead, since a `LocalPotentialOperator` cannot carry a differential operator.

**Measured, and this is the check worth keeping.** moldft builds $\tau$ by
differentiating the physical orbitals across the cusp with no $\mathbf U_1$
anywhere — a wholly different code path to the same physics. On He/TPSS at
`protocol [1e-4, 1e-6, 1e-8]`, nemo gives $-2.90966373$ (one iteration short of
its 1e-06 residual target) against moldft's converged $-2.909663718$, agreeing to
1.2e-08, both 14.4 uHa below TURBOMOLE TPSS/aug-cc-pV6Z. A sign error or a
dropped term in either identity shows up there and nowhere else. `nemo_he_tpss`
registers it at cheap settings.

Note the earlier expectation in this section — that the nemo $\tau$ path would be
*more accurate* than moldft's — is not what was observed: the two agree, and it is
nemo that is far more expensive (see §5.4). The accuracy argument for the
regularized path shows up instead in how fast it converges in `thresh` (§5.4).

Corroborating evidence that second derivatives of the density are the wrong road
here, from the tree itself (`nemo.h`):

> The Laplacian should currently only be used for subsequent convolution with a
> Green's function (which is reasonably stable), but not on its own! … **the
> Laplacian of the regularized density is still very noisy**

i.e. second derivatives are hard *even after* nemo regularization. That is a
strong argument for the $\tau$ route over $\nabla^2\rho$, and for never rewriting
$\tau$ as $\tau_L+\tfrac14\nabla^2\rho$ in this code. `nemo` also carries **four**
competing kinetic-energy routines differing only in how much of $\nabla\psi$ is
handled analytically, the last flagged *"probably imprecise"* — and that is the
one which differentiates $\psi=RF$ directly. The existence and labelling of those
variants is itself the numerical evidence for doing the product rule analytically.

### 4.5 Pointwise consistency of derived quantities — a wrong-answer class

**[M]** 2026-08-24. The general principle is in `XC-IMPLEMENTATION-NOTES.md`
§6.9; the mechanism and the measurements are in §6.10 there. What belongs here is
what it means for this codebase.

The bug: $\chi_{st}=\zeta_s\!\cdot\!\zeta_t$ was formed by `prep_xc_args` as its
own multiwavelet function (`dot(world, grada, gradb)`, truncated, then
`refine_to_common_level`'d), while `make_libxc_args` built
$\nabla\rho_s=\rho_s\zeta_s$ from the $\zeta$ components. At the nucleus of an
open-shell Li atom the two disagree by O(1): stored $\chi_{aa}=-2.07$ against
$+0.196$ from contracting the stored $\zeta$. A sum of squares came out negative,
both diagonal $\sigma$'s hit their $10^{-14}$ floor while $\sigma_{ab}$ kept
$-98.6$, and the total $\sigma$ — which *is* $|\nabla\rho|^2$ — came to $-197$.

Effect: **the correlation kernels return NaN for `vsigma` as soon as the total
goes negative**, in plain GGA as much as meta-GGA. Under libxc 7.1.2 that aborts;
under 7.0.0 it returns finite garbage and the SCF diverges instead — which is why
this presented as a platform-dependent, nondeterministic failure. The *energy*
stayed right to three digits while $\partial e/\partial\tau_\beta$ was wrong by
492×.

Fixed by contracting $\sigma$ pointwise from the same $\zeta$ handed to libxc as
$\nabla\rho$, so it is the exact Gram matrix at every quadrature point, with the
Cauchy–Schwarz bound re-imposed after the diagonal floors (raising a diagonal can
break the bound on its own). The three $\chi$ products and enum slots 22–24 are
gone.

**Four things this codebase should carry forward.**

1. **The ζ/χ scheme's structural guarantees are conditional on pointwise
   contraction.** §2.1's claim was about the algebra; the code has to actually
   realize it. Any future derived quantity with a constraint —
   $\tau\ge\tau_W$, $|\zeta|\le1$, a Gram relation — belongs in
   `make_libxc_args` as a contraction of primitives, not in `prep_xc_args` as a
   projected `Function`.
2. **`rhomin = 0.0` plus absolute floors is a dangerous combination.** `munge`
   maps the tail to exactly zero while `max(1e-14, …)` keeps $\sigma$ strictly
   positive, so the two screening mechanisms can hand libxc a state that is
   individually plausible and jointly impossible. The floors are now applied
   consistently (diagonals, then CS on the cross term), but the general hazard —
   independent absolute floors on quantities that constrain each other — remains
   wherever else it appears. §4.3's three-isosurface list is where to look.
3. **A bound derived from admissible inputs proves nothing about a domain
   error.** A 1440-point scan had bounded the kernel by $|\partial e/\partial\rho|\le1.04$
   "for any input" and was used to rule the kernel out; it had sampled only
   positive-semidefinite $\sigma$. When a measurement contradicts a bound, check
   the bound's premises first.
4. **Watch the potential, not only the energy.** Every energy-based diagnostic
   was blind here. The residual per spin channel was the only thing that showed
   it, and it showed it in the *minority* channel — which is where an open-shell
   meta-GGA is least constrained and most worth instrumenting.

Validation record for the fix, protocol `[1e-4, 1e-6]`, `dft_deriv abgv`, node26
(GNU-13, libxc 7.1.2) with cross-checks on macOS (clang/arm64, libxc 7.0.0):
restricted-vs-unrestricted agreement Be TPSS 5.1e-07, He TPSS 1.9e-07, Li⁺ TPSS
2.7e-07 against the Be PBE control's 6.2e-07; closed-shell references move by
≤3.8e-08 (PBE0 not at all); H2O⁺ PBE/PBE0 reproduce earlier converging-configuration
values to 2e-09; H2O⁺ TPSSh $=-76.004565565$ on both platforms with the default
`abgv`, where it previously needed `bspline` on one; Li TPSS $=-7.489068143$, 286 µHa
below TURBOMOLE UKS/aug-cc-pV5Z. All 17 qctest cases and the 23 chem/qc
short+medium ctest entries pass.

---

## 5. Performance

### 5.1 Derivative applications per SCF iteration

| | ζ (in ctor) | `div` | τ + operator | total |
|---|---|---|---|---|
| LDA | 0 | 0 | 0 | **0** |
| GGA, closed shell | 3 | 3 | 0 | **6** |
| GGA, spin-polarized | 6 × 2 operators = 12 | 6 × 2 = 12 | 0 | **24** |
| meta-GGA | as GGA | as GGA | $3N$ (τ) + $6N$ (operator) | **+9N** |

The polarized count is 24 because `SCF::apply_potential` is called once per spin
and each call constructs a fresh `XCOperator` whose `prep_xc_args` computes
**both** spins' log-gradients — byte-identical output, twice per iteration.
Building the intermediates once and sharing them across the two spin operators is
a pure win and still not done.

Each `grad`/`div` also constructs three fresh `Derivative` objects, which are
`WorldObject`s — collectively constructed and destroyed 3 at a time, twice per
build, every iteration. Caching them on the operator is straightforward.

**The meta-GGA saving worth taking:** the inner $D_x\psi_i$ in `apply_tau_term`
is **the same object** as the $\nabla\psi_i$ needed for $\tau$ in `set_tau`.
Computing the orbital gradients once per iteration and sharing them drops the
marginal cost of the operator over $\tau$ to the $3N$ outer divergence. As built,
they are computed twice. This is the largest remaining meta-GGA win.

So meta-GGA moves the XC cost from $O(1)$ to $O(N)$ derivative applications. The
APW literature reports 2–3× a GGA iteration (companion §6.4); for MADNESS the
factor will be larger, because GGA here is unusually cheap (6 applications total)
and the $\tau$ machinery is unavoidably per-orbital.

### 5.2 The functional is evaluated twice per iteration

`compute_xc_energy()` and `make_xc_potential()` each do their own
`refine_to_common_level` → pointwise evaluation → `truncate` (which
*compresses*). The same `xc_args` are refined, evaluated and compressed **twice**
back-to-back, and libxc is called twice (`xc_gga_exc` then `xc_gga_vxc`;
`xc_mgga_exc` then `xc_mgga_vxc`).

libxc provides fused `xc_{gga,mgga}_exc_vxc`, which share the internal kernel
evaluation. Using them, plus a single refine/evaluate/truncate pass returning both
the energy density and the potential components, roughly halves the XC cost per
iteration. **The single largest easy performance win, and untouched.**
`fxc_apply` has the same pattern (`xc_gga_fxc` + `xc_gga_vxc` where
`xc_gga_vxc_fxc` exists).

### 5.3 Other costs

- `div(semilocal, true)` passes `refine=true`, so the three flux components are
  refined (roughly doubling tree size) before differentiation, every call.
- Removing the $\chi$ products (§4.5) took 3 multiwavelet multiplications plus a
  sum out of the closed-shell path and 9 out of the polarized one, along with
  their `redundant` tree-state changes, their truncations, and their contribution
  to the `refine_to_common_level` level. **[M]** Note the last of these means the
  common level is now set by $\rho$, $\zeta$ and $\tau$ alone; the measured
  effect on converged energies was ≤3.8e-08, well inside threshold.
- Roughly 8–10 fences in `prep_xc_args` and 8–11 in `make_xc_potential` for
  closed-shell GGA, of which only two are explicit; the rest are implicit in
  `reconstruct`/`refine`/`mul`/`sum`/`truncate`/`compress`.
- No caching anywhere: a fresh `XCOperator` per use site per iteration, in
  `SCF.cc`, `nemo.cc` (three sites) and `TDHF.cc` (two).

### 5.4 Measured: nemo meta-GGA is expensive, and the memory is the problem

**[M]** Unlike everything above, these are measurements. He/TPSS at
`protocol [1e-4, 1e-6, 1e-8]`, node26, `MAD_NUM_THREADS=20`:

| path | wall time | outcome |
|---|---|---|
| moldft | 2 m 52 s | converged, -2.909663718 |
| nemo | 56 m 16 s | **killed** at 167 GB resident (of 187 GB), no progress after a fock-matrix build |

So ~20x the wall time for the same number on the same system, and a memory
profile that does not survive the tightest rung. Per-rung timings show where it
goes: ~5 s per energy evaluation at 1e-4, ~190 s at 1e-8 with $k=10$. The
$O(1)\to O(N)$ shift in derivative applications (§5.1) accounts for the time; it
does not by itself account for 167 GB.

Prime suspect is the `refine()` traffic in the $\tau$ machinery. `set_tau` refines
a copy of the orbital vector per axis, and `apply_tau_term` refines both the ket
copy and $\mathbf W$ per axis, so at $k=10$ several refined orbital vectors are
live simultaneously. That makes §6's gradient-sharing item a memory fix as much
as a speed one: the inner $D_x\psi_i$ in `apply_tau_term` is **the same object**
as the $\nabla\psi_i$ `set_tau` needs, and computing it once would cut both the
time and the peak.

Nothing registered goes near this — `nemo_he_tpss` is 24 s — but a molecule at
1e-8 on this path is not currently feasible.

**The regularized path converges faster in `thresh`, and that is where its
accuracy argument lives.** He/PBE, same settings both paths: at thresh 1e-6 nemo
gives -2.892934590 and moldft -2.892864, 70 uHa apart; tightened to 1e-8 they
agree to 6e-09 at -2.89293476. The 70 uHa is moldft's remaining discretization
error, which nemo has already shed two rungs earlier. That is the property
`nemo_he_pbe` asserts on, and the one §4.4's $\zeta$ work should preserve.

---

## 6. Recommended order of work

Done since the first revision: a nemo GGA regression deck (`nemo_he_pbe`), and
meta-GGA on the nemo path with its deck (`nemo_he_tpss`) — the former item 1, the
latter item 8 plus the nemo half of the operator work.

| # | Item | § | Effort | Value |
|---|---|---|---|---|
| 1 | Share the orbital gradients between `set_tau` and `apply_tau_term` | 5.1, 5.4 | small | **now a memory fix, not just 2×** — 167 GB at $k=10$ |
| 2 | Fuse `exc`/`vxc` via `xc_{gga,mgga}_exc_vxc` | 5.2 | small | ~2× on the XC step |
| 3 | Share intermediates across the two spin operators | 5.1 | small | 2× on polarized GGA |
| 4 | Variant-selectable `div`, so `dft_deriv` reaches the flux divergence | 4.2 | small | smoothing where it is needed |
| 5 | `nspin=1` vs `nspin=2` consistency test, and a reference-valued polarized case | 3.3 | small | companion §9 rows 1,2,3,5,8,9 |
| 6 | Wire up analytic ζ/σ in nemo — revive `make_sigma` | 4.4 | medium | the remaining half of the cusp work; `nemo_he_pbe` now measures it |
| 7 | `logme` on a $C^2$ ramp instead of a hard clamp | 4.3 | small | tree size and noise in the tail |
| 8 | A unit-level $\int\tau=T$ check *with* an ncf | 3.3 | small | localizes a nemo τ regression instead of failing an SCF |
| 9 | Housekeeping: `has_fxc`/`has_kxc`, `is_dft` semantics, dead `plot()`/`expme`, unread enum slots 10–13 | 3.2 | small | — |
| 10 | `Nemo::make_fock_operator` refuses a meta-GGA; register the `XCOperator` itself, or add a second operator for the τ term | 4.4 | medium | only TDHF consumes it, and TDHF rejects meta-GGAs anyway |
| 11 | Decide `testxc.cc`'s fate — promote to a test or delete | 3.3 | trivial | — |

Item 1 is now the highest-value entry: it is the difference between a molecule at
thresh 1e-8 being feasible on the nemo meta-GGA path and not (§5.4).

**Explicitly deferred.** Response/TDDFT with $\tau$ — `fxc_apply` refuses
meta-GGA up front, since there is no perturbed kinetic energy density and $\tau$'s
variation is orbital-dependent rather than a functional of $\rho_{\rm pt}$.
Laplacian-level functionals — refused in `initialize()`; slots reserved, nothing
implemented, and the in-tree evidence of §4.4 argues against ever needing them.
OEP/multiplicative meta-GGA potentials — gKS is what every production code uses,
and the OEP's meta-GGA gaps are close to GGA's anyway (companion §5.5).

---

## 7. Verification status

Everything in §1–§3 and §5 was read from the source at the cited locations.
Literature claims are traceable to `XC-IMPLEMENTATION-NOTES.md` and carry its
`[V]`/`[D]`/`[U]` marking. §4.5's measurements are reproducible: the libxc
behaviour in ~40 lines against libxc alone, the SCF results from the decks named
there.

Not verified:

- **Performance numbers in §5 are operation counts from reading the code, not
  measurements**, except where marked **[M]**. The claim that fusing `exc`/`vxc`
  roughly halves the XC step is an inference from the duplicated
  refine/evaluate/compress cycle, not a benchmark.
- **`Nemo::make_sigma` (§4.4) was read, not tested.** No caller since the `chem`
  directory move, so it has ridden through every subsequent refactor of the
  surrounding interfaces without a compiler check on its *semantics*. Treat it as
  a correct derivation and a suspect implementation.
- **The accuracy benefit of the analytic-ζ fix (§4.4) is predicted, not
  demonstrated.** The identity is exact; whether it measurably improves converged
  energies at production thresholds needs a test. `nemo_he_pbe` is now that test's
  baseline (§5.4), so the measurement is available but has not been made.
- **The τ half of §4.4 *is* demonstrated** — two independent code paths agreeing
  to 1.2e-08 (§4.4) — but only on one system, He, closed shell. Nothing has
  exercised the regularized τ on more than one atom, and `set_tau`'s open-shell
  branch has never run with an ncf at all.
- **The BSH conditioning concern** (companion §7) is unsupported by any published
  source and still unmeasured, on either path. $\|\partial e/\partial\tau\|_\infty$
  is logged at `print_level>=2` in moldft; nemo does not set the print level, so
  it logs nothing. Worth wiring up and reading before designing around it — and
  note §5.4's stall at the 1e-8 rung is a candidate symptom, though memory is the
  more likely cause there.
- **`clang-tidy` has not been run** on this subsystem: it is installed on neither
  development machine and neither build tree exports `compile_commands.json`.
- **The sign convention $\mathbf U_1=-\nabla R/R$ is confirmed** from three places
  (`correlationfactor.h` *"NOTE THE SIGN … U1 = -S'/S"*, the functor body, and
  `nemo.h` *"note: U1=-grad(R)/R"*). The file header's summary writes
  $\vec U_1=R^{-1}(\vec\nabla R)$ without the minus, and the comment describing
  `U2` as *"-S''/S - Z/r"* is missing a factor $\tfrac12$ (the small-$r$ limit of
  `Spp_div_S` is $Z^2(1-1.5a)$, matching $-\tfrac12\Delta S/S-Z/r$). **The code is
  authoritative; that header block is not.** Worth fixing while in the area.
