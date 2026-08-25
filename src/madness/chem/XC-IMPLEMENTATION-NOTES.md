# Working equations for LDA, GGA and meta-GGA in a real-space code

Literature notes, assembled 2026-08-22, condensed 2026-08-24. Scope: the
exchange–correlation working equations needed for rungs 1–3 of Jacob's ladder in
a real-space / wavelet code, the libxc conventions they must be expressed in,
and the numerical traps between the two. Companion:
`XC-MADNESS-ASSESSMENT.md`, which measures MADNESS against this.

Atomic units. Spin labels $\sigma,\sigma'\in\{\alpha,\beta\}$,
$\rho=\rho_\alpha+\rho_\beta$.

**Verification legend.** **[V]** read verbatim from a primary source (§12);
**[D]** derived here from **[V]** definitions; **[M]** measured against libxc or
MADNESS in this repo; **[U]** unverified — standard citation, text not opened.
*Do not propagate [U] into a paper without checking.* The split is load-bearing:
several widely-repeated "facts" here are convention-, version- or
implementation-dependent.

---

## 1. Common setup

$$
E_{xc} = \int f\bigl(\rho_\sigma,\nabla\rho_\sigma,\tau_\sigma,\nabla^2\rho_\sigma\bigr)\,d^3r
$$

with the contracted gradient invariants

$$
\sigma_{\alpha\alpha}=\nabla\rho_\alpha\!\cdot\!\nabla\rho_\alpha,\quad
\sigma_{\alpha\beta}=\nabla\rho_\alpha\!\cdot\!\nabla\rho_\beta,\quad
\sigma_{\beta\beta}=\nabla\rho_\beta\!\cdot\!\nabla\rho_\beta .
$$

$\sigma_{\alpha\beta}$ is the **plain** dot product, no factor 2 (§2.1).

| Rung | Ingredients | Potential |
|---|---|---|
| LDA/LSDA | $\rho_\sigma$ | multiplicative, pointwise |
| GGA | $+\,\sigma_{\sigma\sigma'}$ | multiplicative, one divergence |
| meta-GGA | $+\,\tau_\sigma$ (rarely $\nabla^2\rho_\sigma$) | **non-multiplicative** — a differential operator |

That last row is the qualitative break; §5 is about it.

**The $\sigma$ matrix is a Gram matrix.** With $\mathbf g_\sigma=\nabla\rho_\sigma$,
$\sigma_{st}=\mathbf g_s\!\cdot\!\mathbf g_t$, so

$$
\sigma_{ss}\ge0,\qquad
|\sigma_{\alpha\beta}|\le\sqrt{\sigma_{\alpha\alpha}\sigma_{\beta\beta}},\qquad
\sigma_{\alpha\alpha}+2\sigma_{\alpha\beta}+\sigma_{\beta\beta}=|\nabla\rho|^2\ge0 .
$$

All three are consequences of one object, not three independent constraints.
libxc depends on all three; violating any of them is a domain error, not an
accuracy loss. §6.9–6.10 are about how a real-space code violates them by
accident, and what it costs.

---

## 2. libxc interface conventions

Verified against the libxc 5.1.x and 3.0.x manuals and the `devel` source
(`src/xc.h`, `functionals.c`, `work_gga_inc.c`, `work_mgga_inc.c`, `gga.c`,
`scripts/maple2c_lib/gga.py`, `ChangeLog.md`).

### 2.1 Inputs

**[V]**
```
xc_lda_*  (p, np, rho,                    <outputs>)
xc_gga_*  (p, np, rho, sigma,             <outputs>)
xc_mgga_* (p, np, rho, sigma, lapl, tau,  <outputs>)
```
`rho[0,1]`$=\rho_{\uparrow,\downarrow}$; `sigma[0,1,2]`$=(\uparrow\uparrow,\uparrow\downarrow,\downarrow\downarrow)$;
`lapl[0,1]`, `tau[0,1]` per spin. `XC_UNPOLARIZED==1`, `XC_POLARIZED==2` — the
*number of spin channels*, not a bit flag. "The routines expect a nonnegative
density."

**No factor 2 in `sigma[1]`, confirmed twice independently.** The codegen
(`maple2c_lib/gga.py`) builds the polarized total as
`sqrt(sigma0 + 2*sigma1 + sigma2)`; that explicit `2*` is correct only if
`sigma[1]` itself carries none. **[V]**

### 2.2 `zk` is per particle

**[V]** $\texttt{zk}\equiv\texttt{exc}=e=\varepsilon/\rho$, so
$E_{xc}=\int\rho\,\texttt{zk}$. There is exactly **one** `exc` component per
point in both spin cases. Every `v*` output differentiates $\varepsilon$ (per
volume), not $e$ — for LDA libxc has already applied the product rule:

$$
\texttt{vrho}_\sigma=\frac{\partial(\rho\,\texttt{zk})}{\partial\rho_\sigma}
=\texttt{zk}+\rho\frac{\partial\texttt{zk}}{\partial\rho_\sigma}
\neq\frac{\partial\texttt{zk}}{\partial\rho_\sigma}.
$$

**Never multiply `vrho` by $\rho$.** A code accumulating $E_{xc}=\int\texttt{zk}$
while using $v_{xc}=\texttt{vrho}$ fails every finite-difference and virial check
while looking plausible. **[U]** No libxc release used a per-volume `zk`; treat
per-volume as a *cross-code* hazard (older non-libxc layers), not a version
change.

### 2.3 Output orderings

**[V]** `nspin=2`:

| Array | Index list | Count |
|---|---|---|
| `vrho` | $(0),(1)$ | 2 |
| `vsigma` | $(0),(1),(2)$ | 3 |
| `v2rho2` | $(0,0),(0,1),(1,1)$ | 3 |
| `v2rhosigma` | $(0,0)\ldots(1,2)$ | 6 |
| `v2sigma2` | $(0,0),(0,1),(0,2),(1,1),(1,2),(2,2)$ | 6 |

**The two 6-element arrays are packed differently.** `v2rhosigma` is the full
$2\times3$ **rectangle**, row-major (offset $3\tau+st$) — different variables, no
symmetry. `v2sigma2` is the **upper triangle of a symmetric $3\times3$** —
off-diagonals stored once, used twice in any contraction. Both are `6*np`
doubles, so confusing them crashes nothing. **[V]**

**[V]** Size rule: a variable with $s$ spin components differentiated $k$ times
contributes $\binom{s+k-1}{k}$; the array size is the product over variables in
the name, ordered `[rho, sigma, lapl, tau]`. $s=2$ for `rho`/`lapl`/`tau`, $s=3$
for `sigma`. Worked: `v4rho2sigma2`$=3\times6=18$, `v4sigma2lapltau`$=6\times3\times2=36$.

**[D]** meta-GGA second derivatives, from that rule (the manual does not
tabulate them): `v2rholapl`, `v2rhotau`, `v2lapltau` $=4$ (rectangle);
`v2sigmalapl`, `v2sigmatau` $=6$ (rectangle); `v2lapl2`, `v2tau2` $=3$
(triangle). The *counts* are rigorous; row-major for the rectangles is inferred
by analogy with the **[V]** `v2rhosigma` listing. Confirm by finite difference
before depending on it.

**[V]** Unwanted outputs may be `NULL`. Check
`XC_FLAGS_HAVE_{EXC,VXC,FXC,KXC,LXC}` before requesting order $n$. Fused entry
points (`xc_gga_exc_vxc`, `xc_gga_exc_vxc_fxc`, `xc_mgga_exc_vxc`, …) share one
kernel evaluation — worth using for both potential-and-energy and
potential-and-kernel.

### 2.4 The $\tau$ convention

**[V]** manual: "the kinetic energy density is defined with the 1/2 factor
(until Libxc version 2.0.0, tau was defined as twice the kinetic energy
density): $\tau_\alpha=\tfrac12\sum_i|\nabla\psi_{\alpha i}|^2$" — positive
definite, no Laplacian folded in (`lapl` is separate).

A documentation conflict, adjudicated: the 5.1.x manual and `ChangeLog.md`
(`[2.0.0] 2012-11-28`, "Definition of tau in the metaGGAs. Now tau is the exact
kinetic energy density (with the factor 1/2 included)") agree on **2.0.0**; the
3.0.x page says 1.2 and is stale. Both state the identical formula, so only the
cutoff was ever in doubt. **Any libxc ≥ 2.0.0 is the $\tfrac12$ convention.**

**Independent confirmation from the source [V]+[D].** `work_mgga_inc.c` enforces
the von Weizsäcker / Fermi-hole-curvature bound
`my_sigma[0] = m_min(my_sigma[0], 8.0*my_rho[0]*my_tau[0])`, i.e.
$\sigma_{ss}\le8\rho_s\tau_s$. For one occupied real orbital
$|\nabla\rho_s|^2=4\rho_s|\nabla\phi|^2$, and with $\tau_s=\tfrac12|\nabla\phi|^2$
that is exactly $8\rho_s\tau_s$ (equality, as it must be). With
$\tau_s=|\nabla\phi|^2$ the constant would be 4. The `8.0` pins the $\tfrac12$.

**[V]** The wider literature diverges: Ye (PRB **91**, 075101) writes $\tau$
*without* the $\tfrac12$; TPSS/SCAN/r²SCAN and Doumont et al. carry it; the
Minnesota lineage commonly omits it **[U]**. *Read the prefactor in each paper;
never infer it.*

### 2.5 Unpolarized → polarized mapping

**[V]** from `maple2c_lib/gga.py`, which *generates* the `nspin=1` kernel by
substitution `mf := (r0,r1,s0,s1,s2) -> eval(2*f_chan(r0/2, s0/4))`:

$$
\rho_\alpha=\rho_\beta=\rho/2,\qquad
\sigma_{\alpha\alpha}=\sigma_{\alpha\beta}=\sigma_{\beta\beta}=\sigma/4 .
$$

Consistent with §2.1 **[D]**: $\tfrac\sigma4+2\tfrac\sigma4+\tfrac\sigma4=\sigma$ ✓.

**[D]** The derivative mapping — the part codes get wrong:

$$
\texttt{vrho}_{\rm unpol}=\texttt{vrho}_\alpha,\qquad
\texttt{vsigma}_{\rm unpol}=\tfrac14\bigl(2\,\texttt{vsigma}_{\alpha\alpha}+\texttt{vsigma}_{\alpha\beta}\bigr).
$$

`vrho` is unchanged between conventions; `vsigma` picks up the $\tfrac14$ and the
channel sum.

---

## 3. LDA / LSDA

$$
v_{xc,\sigma}=\frac{\partial f}{\partial\rho_\sigma}=\texttt{vrho}[\sigma],
$$
purely pointwise, no derivative operators.

**Dirac/Slater exchange — the reference case.** **[V]**
$f_x=-\tfrac32(3/4\pi)^{1/3}(\rho_\alpha^{4/3}+\rho_\beta^{4/3})$, hence **[D]**
$v_{x,\sigma}=-2(3/4\pi)^{1/3}\rho_\sigma^{1/3}$, reducing at
$\rho_\alpha=\rho_\beta=\rho/2$ to $v_x=-(3/\pi)^{1/3}\rho^{1/3}$. **This is the
consistency check for the whole spin-convention chain — run it first.** It is
also the origin of the $\rho^{1/3}$ non-smoothness: $v_x$ is finite but has
infinite derivative wherever $\rho\to0$, and a kink at every nuclear cusp.

Correlation parametrizations (PW92, VWN) were **not** read from the original
papers and their formulas are omitted here rather than reproduced from memory.
Verify against libxc's `maple/lda_exc/lda_c_pw.mpl` and `lda_c_vwn*.mpl`. Note
VWN3 vs VWN5 differ in the spin-stiffness fit and quantum-chemistry codes
disagree about the default.

---

## 4. GGA

### 4.1 Variation, and where the factors come from

$E_{xc}=\int f(\rho_\alpha,\rho_\beta,\sigma_{\alpha\alpha},\sigma_{\alpha\beta},\sigma_{\beta\beta})$
— Pople/Gill/Johnson (PGJ) Eq. (7) **[V]**. **[D]**

$$
\delta\sigma_{\alpha\alpha}=2\nabla\rho_\alpha\!\cdot\!\nabla\delta\rho_\alpha,\qquad
\delta\sigma_{\alpha\beta}=\nabla\rho_\beta\!\cdot\!\nabla\delta\rho_\alpha+\nabla\rho_\alpha\!\cdot\!\nabla\delta\rho_\beta .
$$

The **factor 2** on the same-spin term is just
$\partial\sigma_{\sigma\sigma}/\partial\nabla\rho_\sigma=2\nabla\rho_\sigma$ — that
invariant is *quadratic*. The cross term carries **no** 2 because
$\sigma_{\alpha\beta}$ is *bilinear*; its coefficient is the gradient of the
**opposite** spin density. That is the term most often forgotten, and the reason
open-shell GGA cannot be assembled channel-by-channel.

**XC gradient flux vector:**

$$
\boxed{\;\mathbf X_\sigma=2\frac{\partial f}{\partial\sigma_{\sigma\sigma}}\nabla\rho_\sigma
+\frac{\partial f}{\partial\sigma_{\alpha\beta}}\nabla\rho_{\sigma'}\;}
$$

so $\delta E_{xc}=\sum_\sigma\int[\,\partial f/\partial\rho_\sigma\,\delta\rho_\sigma+\mathbf X_\sigma\!\cdot\!\nabla\delta\rho_\sigma]$,
and one integration by parts gives §4.2.

**The boundary term is the one that gets silently dropped.** It vanishes when
$\delta\rho\to0$ on $\partial\Omega$ (a box where $\rho$ is numerically zero at
the walls) or under PBC — *not* for a tightly truncated cell. There the
divergence and weak forms are **different operators**, the potential is no longer
the derivative of the energy being computed, and the SCF loses variationality.

### 4.2 Form (a): divergence form / multiplicative potential

$$
\boxed{\;v_{xc,\sigma}=\frac{\partial f}{\partial\rho_\sigma}-\nabla\!\cdot\!\mathbf X_\sigma\;}
$$

i.e. **[V]** (PGJ Eq. 9; Yanai/Harrison/Handy Eq. 12)

$$
v_{xc,\alpha}=\texttt{vrho}_\alpha-\nabla\!\cdot\!\bigl(2\,\texttt{vsigma}_{\alpha\alpha}\nabla\rho_\alpha+\texttt{vsigma}_{\alpha\beta}\nabla\rho_\beta\bigr),
$$
unpolarized $v_{xc}=\texttt{vrho}-2\nabla\!\cdot\!(\texttt{vsigma}\,\nabla\rho)$.

**Closed-shell cross-check [D].** Inserting the §2.5 mapping into the unpolarized
formula reproduces the polarized one exactly, with no leftover factor:

$$
\texttt{vrho}_\alpha-\nabla\!\cdot\!\Bigl(2\cdot\tfrac14(2\texttt{vsigma}_{\alpha\alpha}+\texttt{vsigma}_{\alpha\beta})\cdot2\nabla\rho_\alpha\Bigr)
=\texttt{vrho}_\alpha-\nabla\!\cdot\!\bigl(2\texttt{vsigma}_{\alpha\alpha}\nabla\rho_\alpha+\texttt{vsigma}_{\alpha\beta}\nabla\rho_\beta\bigr)\;\checkmark
$$

This validates the 2/1 pattern, the $\sigma/4$ mapping and the
$\tfrac14(2\texttt{vsigma}_{\alpha\alpha}+\texttt{vsigma}_{\alpha\beta})$ relation
at once. **Put it in a unit test.**

**Never expand $\nabla\!\cdot\!\mathbf X$ analytically.** Doing so brings in the
full Hessian $\nabla\nabla\rho$ and the functional's second derivatives, even for
a plain SCF step **[D]** (PGJ note only that "Hessian matrices … are involved"
**[V]**). Instead build $\mathbf X_\sigma$ as a real-space *function* from
quantities already at hand (`vsigma` from libxc, $\nabla\rho$ from the code's
derivative operator) and apply **one** divergence. Only first derivatives ever
touch the density.

**In the $|\nabla\rho|$ convention it is explicitly singular.** With
$g=|\nabla\rho|$, Balbás/Martins/Soler Eq. (2) **[V]** has $1/g$ and $1/g^2$
terms, so form (a) diverges wherever $\nabla\rho=0$ — every nucleus, every bond
critical point, all of vacuum. Verbatim **[V]**: "Since $|\nabla\rho|$ has cusps
at extrema of $\rho$, $\nabla g$ is discontinuous at these points, what causes
problems for its numerical representation." The $\sigma$ convention hides this
but does not remove the underlying cusp. Conversion:
$\partial f/\partial\sigma=(1/2g)\,\partial f/\partial g$.

### 4.3 Forms (b) and (c), and the choice

**(b) Weak / matrix-element form** — PGJ Eq. (16), reproduced verbatim by
Lehtola et al. Eqs. (51)–(52) **[V]**:

$$
\langle\phi|v_{xc,\sigma}|\psi\rangle=\int\Bigl[\frac{\partial f}{\partial\rho_\sigma}\phi\psi+\mathbf X_\sigma\!\cdot\!\nabla(\phi\psi)\Bigr].
$$

No Hessians — PGJ call that "a major computational advantage" **[V]**. The cost
in a real-space code is that $v_{xc}$ is no longer a *function*: BSH/Green's
iteration, KAIN updates on $\psi$, potential plotting and asymptotic corrections
all break.

**(c) White–Bird** — the derivative of the *discretized* energy,
$\tilde v_i=(1/w_i)\partial\tilde E_{xc}/\partial\rho_i$
(Balbás/Martins/Soler Eqs. 3–5 **[V]**), structurally the discrete adjoint of the
§4.1 integration by parts **[D]**. Verbatim **[V]**: Eqs. (2) and (5) are
"asymptotically equivalent … but different in practice", and (5) "and not
Eq. (2), is the 'correct' definition of $v_{xc}$" if the discretized form is what
is minimized. This is what GPAW does **[V]**, reporting that the discretization
"tend[s] to smooth out the large pathological wiggles developed by the PW91-GGA
functional in regions of zero density or zero density gradient" — and its own
footnote, "probably the most complicated part of the code".

**[D] Why (c) does not transfer to an adaptive multiwavelet basis:**
$\partial\tilde E/\partial\rho_i$ presupposes a fixed, refinement-independent
coefficient set. Under adaptive refinement the set changes, so the quantity is
not a well-defined local potential. Hence the multiwavelet lineage uses form (a).

| | needs from $\rho$ | needs from $\psi$ | main hazard |
|---|---|---|---|
| (a) divergence | $\nabla\rho$, one div of $\mathbf X$ | nothing | noise: $\mathbf X\sim\rho^{-4/3}\nabla\rho$ in the tail, then differentiated |
| (b) weak form | $\nabla\rho$ only | $\nabla\phi,\nabla\psi$ | no potential *function* — breaks BSH |
| (c) White–Bird | $\nabla\rho$ from your own operator | nothing | needs a fixed coefficient set |

**One rule spanning all three:** use the same derivative operator for $E_{xc}$
and $v_{xc}$. Mixing an analytic $\nabla\rho$ in the energy with a
finite-difference one in the potential is the inconsistency White–Bird exists to
fix.

### 4.4 The identity worth keeping

If flux-vector noise in form (a) bites, this **[D]** moves the differentiation
off $\mathbf X$:

$$
\hat v_{xc}\psi=\frac{\partial f}{\partial\rho_\sigma}\psi-\nabla\!\cdot\!(\mathbf X_\sigma\psi)+\mathbf X_\sigma\!\cdot\!\nabla\psi .
$$

Since $\psi$ decays, $\mathbf X\psi$ is far better conditioned in the asymptotic
region than $\mathbf X$ alone, and $\nabla\psi$ is smooth. It keeps a
multiplicative *action* while improving the tail — but not a multiplicative
*potential*.

---

## 5. meta-GGA

### 5.1 The two kinetic energy densities, and the gauge trap

**[V]** $\tau=\tfrac12\sum_i|\nabla\psi_i|^2$ (positive definite),
$\tau_L=-\tfrac12\sum_i\psi_i\nabla^2\psi_i$, and
$\tau_L=\tau-\tfrac14\nabla^2\rho$. Zahariev/Leang/Gordon **[V]**: the two
"differ by a full differential … that integrates to zero", so both give the same
kinetic energy — *but*

> Although all of these expressions upon integration give the same kinetic
> energy value, the use of different forms … would lead to different
> exchange-correlation energies. … the additional full differential expression
> does not simply integrate away.

**$\tau$ and $\tau_L$ are interchangeable inside a *linear* integral; a meta-GGA
is a nonlinear function of its argument.** Substituting $\tau_L$ defines a
**different functional**. You cannot "simplify" TPSS or SCAN by switching to the
Laplacian form and claim the same energy.

### 5.2 Why $\delta E_{xc}/\delta\rho$ has no closed form

$\tau$ is built from orbitals, which depend on $\rho$ only implicitly.
Zahariev et al. Eq. (8) **[V]**: the chain rule needs
$\delta\psi_i[\rho]/\delta\rho$, i.e. the KS response function and its
**inverse** — an OEP. Yang/Peng/Sun/Perdew Eqs. (4)–(6) **[V]** split it as
$v_{xc,\sigma}=v^{\rm GGA}_{xc,\sigma}+\int\!d^3r'\,(\partial e_{xc}/\partial\tau_\sigma)(\mathbf r')\,\delta\tau_\sigma(\mathbf r')/\delta n_\sigma(\mathbf r)$,
verbatim: the first "is a multiplicative potential, but $v^{\tau-\rm dep}$ does
not have a closed form."

### 5.3 The generalized Kohn–Sham working equation

Differentiate with respect to the **orbital**, not the density:

$$
\boxed{\;\bigl(\hat v_{xc,\sigma}\psi_i\bigr)\supset-\frac12\nabla\!\cdot\!\left[\frac{\partial e_{xc}}{\partial\tau_\sigma}\nabla\psi_i\right]\;}
$$

**Verified verbatim in two independent primary sources** — Yang et al. Eq. (7)
**[V]** and Doumont/Tran/Blaha Eq. (6) **[V]**, the latter giving the complete
operator. Weak form, Yang et al. Eq. (16) **[V]**:
$-\tfrac12\int\psi_i^*\nabla\!\cdot\![v_\tau\nabla\psi_j]=\tfrac12\int v_\tau\,\nabla\psi_i^*\!\cdot\!\nabla\psi_j$.

**Normalization check [D]:** set $e_{xc}=\tau$, so $v_\tau=1$; the operator gives
$-\tfrac12\nabla^2\psi_i$, exactly the derivative of
$T_s=\tfrac12\int|\nabla\psi|^2$ ✓. The $\tfrac12$ is inherited from libxc's
$\tau$; a code using $\tau=\sum|\nabla\phi|^2$ needs
$-\nabla\!\cdot\!(v_\tau\nabla\psi)$ instead.

**Laplacian contribution [D].** Two integrations by parts give sign $(+1)^2$:
$v_{xc,\sigma}\supset+\nabla^2(\partial e_{xc}/\partial\nabla^2\rho_\sigma)$ —
**multiplicative**, unlike the $\tau$ term, but putting *four* net derivatives of
$\rho$ into the potential. libxc gates it behind `XC_FLAGS_NEEDS_LAPLACIAN`
**[V]**: the library itself treats it as exceptional.

### 5.4 The complete operator, and the effective-mass reading

$$
\boxed{\begin{aligned}
\hat v_{xc,\sigma}\psi&=V_{\rm mult}\psi-\tfrac12\nabla\!\cdot\!(v_{\tau,\sigma}\nabla\psi)\\
V_{\rm mult}&=\frac{\partial e}{\partial\rho_\sigma}-\nabla\!\cdot\!\Bigl[2\frac{\partial e}{\partial\sigma_{\sigma\sigma}}\nabla\rho_\sigma+\frac{\partial e}{\partial\sigma_{\sigma\sigma'}}\nabla\rho_{\sigma'}\Bigr]\;\bigl[+\nabla^2(\partial e/\partial\nabla^2\rho_\sigma)\bigr]\\
v_{\tau,\sigma}&=\texttt{vtau}[\sigma]
\end{aligned}}
$$

Adding the kinetic operator:

$$
\hat T+\hat v_\tau=-\frac12\nabla\!\cdot\!\bigl[(1+v_{\tau,\sigma})\nabla\bigr],
$$

a position-dependent inverse effective mass $1/m^*=1+\partial e_{xc}/\partial\tau_\sigma$.
Yang et al. **[V]** put the design constraint in one sentence: the gKS potential
is "a non-multiplicative self-adjoint operator … **a differential operator for a
meta-GGA and an integral operator for a hybrid**". *meta-GGA is to the Laplacian
what exact exchange is to the Coulomb operator.*

**Ellipticity condition [D]** — a cheap and essential runtime assertion:

$$
1+\frac{\partial e_{xc}}{\partial\tau_\sigma}(\mathbf r)>0\quad\forall\mathbf r .
$$

Below $-1$ the local effective mass changes sign, the quadratic form
$\tfrac12\int(1+v_\tau)|\nabla\psi|^2$ loses positive-definiteness, and a
high-frequency mode's "kinetic energy" goes negative — a variational collapse
that a real-space code expresses as uncontrolled refinement at the finest scale.
**[U]** as a literature claim; the algebra is elementary. **[M]** caution: a
lattice-sampled $\min v_\tau$ is an indicator, not a bound, and sampling the cell
centre of an atom samples the nuclear cusp — where one point out of 241 is
genuinely non-elliptic while the neighbourhood is smooth. Do not act on a
single-point warning.

### 5.5 Lineage, and a factor-of-2 warning about the canonical derivation

Both sources read attribute the practical scheme to Neumann–Nobes–Handy (1996)
**[V]**; Zahariev et al. rederive it as the ODDM and expose its hidden
approximation — that
$\delta E_{xc}[\tau[\rho]]/\delta\psi_i\approx\delta E_{xc}[\tau]/\delta\psi_i^{\rm KS}$
equates functional derivatives that "differ in meaning" **[V]**. Their
conclusion: the ODDM gives "non-unique derivatives … and **non-local**
exchange-correlation potentials". Yang et al. instead treat gKS as a
self-contained variational scheme whose energies lie *below* the OEP's **[V]**.
Both are right about different questions: gKS is variationally sound as an
*orbital* theory but does not deliver $\delta E_{xc}/\delta\rho$.

> ### ⚠ Factor-of-2 inconsistency *inside* Zahariev et al. (2013)
>
> Their Eq. (21) has $\delta\tau/\delta\psi_i^{\rm KS}=2\nabla\psi_i\nabla\delta$.
> With their own $\tau=\tfrac12\sum|\nabla\psi_i|^2$ and real orbitals the correct
> prefactor is **1**. The stray 2 propagates into Eqs. (22), (25), (26), (36), so
> their Eq. (26) lacks the $\tfrac12$ that Yang et al. Eq. (16) carries. Their own
> Appendix (A1)/(A2) *is* consistent with the $\tfrac12$ convention
> ($\tfrac14$ / $\tfrac12$ / $1$ prefactors, one $\tfrac12$ per $\tau$-derivative).
> **The Appendix is right; the main text carries a typo. Use Yang et al. Eq. (16)
> and Doumont et al. Eq. (6) as the authority.** This matters because Zahariev
> et al. is the most-cited canonical derivation.

**What "improved meta-GGA gaps" really means.** Yang et al. **[V]**: the OEP
meta-GGA $v_{xc}$'s "are close to the GGA $v_{XC}$'s … Consequently the OEP
meta-GGA gaps … are not improved significantly over the GGA gaps." The celebrated
gap improvement is an artifact of the **non-multiplicative operator**, not of the
density functional. Papers rarely say which they mean.

---

## 6. Real-space numerics

### 6.1 Obtaining $\nabla\rho$ and $\tau$

$\nabla\rho_\sigma=\mathbf D\rho_\sigma$ (A) or
$\sum_i2\phi_i\nabla\phi_i$ (B); both are in production use. Balbás et al. argue
the *discrete* gradient can be the better object **[V]**: "in practice $g_i$ is
frequently a better approximation than $g(\mathbf r_i)$ to the average value …
within the spatial 'pixel'", provided the density grid has twice the cutoff of
the wave functions. GPAW takes (A) with a finite-difference operator **[V]**. The
conditioning argument for (B) — $\rho$ needs roughly twice the bandwidth of
$\phi$, so differentiating at the smoother orbital level keeps amplification down
— is **[U]** as stated.

For $\tau$, **route A is what every verified source uses** — libxc, Yang et al.
Eq. (3), Doumont et al. Eq. (5), r²SCAN, Ye **[V]**:
$\tau_\sigma=\tfrac12\sum_i w_i|\nabla\psi_{i\sigma}|^2$. One gradient per
occupied orbital, positive by construction (which is what makes $\alpha$ well
defined), and **the same gradients are reused in the operator application**.

Route B eliminates orbital gradients via the KS equation
($\tau_L=\sum_i\epsilon_i|\psi_i|^2-v_{\rm eff}\rho$, $\tau=\tau_L+\tfrac14\nabla^2\rho$).
**[D]** Its drawbacks are disqualifying: exact only at self-consistency with
exact eigenpairs, so away from convergence $\tau$ becomes a function of the SCF
error — dangerous inside a variational loop and fatal for fixed-point analysis;
it trades well-conditioned $\nabla\psi$ for badly-conditioned $\nabla^2\rho$ plus
a near-cancelling difference of two large core terms; and it does not guarantee
$\tau\ge0$, so $\alpha$ can go negative.

> **Correction of a common misattribution.** Ye, PRB **91**, 075101 (2015) is
> frequently cited for a direct-gradients-vs-KS-equation comparison. Its abstract
> presents *one* scheme and contains no such comparison. Its verified value here
> is the missing-$\tfrac12$ convention divergence, and the demonstration that a
> better $\tau$ evaluation fixes a convergence pathology of the mBJ potential.

### 6.2 The single most important implementation trick

Yang et al. **[V]**: "the grid sensitivity of gKS is much lower than that of OEP.
There is no contradiction since $\nabla\partial e_{XC}/\partial\tau$ does not have
to be evaluated directly in gKS meta-GGA. **Using integration by parts** …"

**Never form $\nabla(\partial e_{xc}/\partial\tau)$:**

$$
\boxed{\;\nabla\!\cdot\!(v_\tau\nabla\psi)=\sum_{x=1}^{3}D_x\bigl(v_\tau\,D_x\psi\bigr)\;}
$$

One gradient of $\psi$, three multiplications, one divergence — $v_\tau$ is only
ever *multiplied*. **[D]** This form is also self-adjoint **by construction**,
since $\langle\phi|-\tfrac12\nabla\!\cdot\!(v\nabla)|\psi\rangle=\tfrac12\int v\nabla\phi\!\cdot\!\nabla\psi$
is exactly symmetric, whereas the expanded
$-\tfrac12(v\nabla^2\psi+\nabla v\!\cdot\!\nabla\psi)$ is symmetric only up to
discretization error *and* needs $\nabla^2\psi$, which has a cusp at every
nucleus. The same reasoning is why the OEP is grid-hungry: its Slater term
contains $\nabla\partial e_{XC}/\partial\tau$ explicitly and cannot integrate it
away **[V]**.

**With a regularizing prefactor, the same trick applies twice. [M]** If the code
solves for $F$ with $\psi=RF$ for some smooth positive $R$ (a nuclear correlation
factor, a Jastrow prefactor, an augmentation), then $\tau$ and the operator both
want the product rule rather than a numerical derivative of $\psi$:

$$
\nabla\psi=R(\nabla F-\mathbf U_1F),\qquad \mathbf U_1=-\nabla R/R
$$
$$
\tau_\sigma=\tfrac12R^2\sum_iw_i|\nabla F_i-\mathbf U_1F_i|^2
$$

and for the operator, writing $\mathbf W=v_\tau(\nabla F-\mathbf U_1F)$ so that
$v_\tau\nabla\psi=R\mathbf W$,

$$
\nabla\!\cdot\!(R\mathbf W)=R\bigl(\nabla\!\cdot\!\mathbf W-\mathbf U_1\!\cdot\!\mathbf W\bigr)
\;\Longrightarrow\;
R^{-1}\Bigl[-\tfrac12\nabla\!\cdot\!(v_\tau\nabla\psi)\Bigr]
 = -\tfrac12\bigl(\nabla\!\cdot\!\mathbf W-\mathbf U_1\!\cdot\!\mathbf W\bigr).
$$

**The $R$ factors cancel identically — nothing is divided by $R$**, and if
$\mathbf U_1$ is analytic (as it is when $R$ is chosen to carry the nuclear cusp)
then every remaining derivative acts on the smooth $F$. Forming $\psi=RF$ and
differentiating that instead puts the cusp back under the derivative operator,
which defeats the regularization. Verified in MADNESS by agreement to 1.2e-08
between this route and the direct one on He/TPSS (`XC-MADNESS-ASSESSMENT.md`
§4.4) — but note the regularized route was ~20x more expensive there, and its
memory profile did not survive the tightest threshold (§5.4 of the same).

### 6.3 Where the noise comes from — the $\alpha\approx1$ mechanism

The most useful diagnostic result in the literature surveyed. Yang et al.
**[V]**: "**SCAN needs the biggest grid to converge**", and the explanation —

> $f$ decreases monotonically with $\alpha$ from 1 at $\alpha=0$ to 0 at
> $\alpha=1$ … the $f(\alpha)$ of SCAN and MS2 are constructed to have
> $d^2f/d\alpha^2=0$ there … Therefore **$df/d\alpha$ has a flat- or
> linearly-topped peak at $\alpha=1$, which explains the oscillation.** … **TPSS
> does not have this feature since it uses $z$ instead of $\alpha$.**

with $\alpha=(\tau-\tau_W)/\tau_{\rm unif}$, $\tau_W=|\nabla n|^2/8n$,
$\tau_{\rm unif}=\tfrac{3}{10}(3\pi^2)^{2/3}n^{5/3}$ **[V]**.

**Two actionable conclusions.** The noise is localized near $\alpha\approx1$ —
spatially, wherever the system crosses between single-orbital and uniform-gas
character, i.e. bonding regions and the density tail. And **TPSS is structurally
easier on a grid than SCAN**, being built on $z=\tau_W/\tau$. Bring up TPSS
first, then r²SCAN, then SCAN. That ordering is derived, not folklore.

r²SCAN's fix is *explicitly numerical* **[V]**: "SCAN's utility for large scale
projects is limited by its sensitivity to the density of the numerical
integration grid"; "constraining the interpolation function to have zero
derivatives at $\alpha=1$ results in a twisted function that harms numerical
performance"; "It would be difficult to assert that any of the grid settings
present a converged SCAN energy, with SCAN errors varying unpredictably by a
factor of 2." The regularization is
$\bar\alpha=(\tau-\tau_W)/(\tau_{\rm unif}+\eta\tau_W)$, $\eta=10^{-3}$, whose
$\eta\tau_W$ specifically fixes the $\alpha\to0$ ($\tau\to\tau_W$) $0/0$.

### 6.4 Measured costs (literature)

- **[V]** Doumont/Tran/Blaha (APW): SCAN needs $G_{\max}=20$–24 vs 12 for PBE;
  MGGA potential evaluation "2–3× more expensive" than GGA, "primarily from
  doubling the Fourier grid". MoS₂ with TASK took 65 SCF iterations vs 9–15.
- **[V]** Yang et al.: for Ne the best built-in BAND grid has 120 radial points;
  the SCAN **OEP** needed ≥266, and 504 for a tighter criterion. **gKS was fine
  on the built-in grids.**
- **[V]** Yang et al. on robustness: "the OEP meta-GGA in general needs more SCF
  cycles… For small gap materials, the OEPs of some meta-GGA functionals do not
  converge, while the corresponding gKS calculations converge normally. … SCAN
  has the least convergence problem." **Note this is the opposite ranking from
  the grid-size finding** — one is integration-grid size, the other OEP SCF
  convergence. Quote with context.

### 6.5 Non-smoothness of $e_{xc}$

All **[V]**: $\rho^{4/3},\rho^{1/3},r_s^{1/2},r_s^{3/2}$ are branch points at
$\rho=0$, so $e_{xc}$ is nowhere analytic where the density vanishes and a
polynomial (multiwavelet) basis converges only algebraically there; $r_s\to\infty$
overflows PW92's polynomial long before $\epsilon_c\to0$ analytically;
$\zeta\to0/0$ in vacuum and $f(\zeta)$'s $(1\pm\zeta)^{4/3}$ has infinite second
derivative at $|\zeta|=1$, where a fully polarized tail sits (hence
`xc_func_set_zeta_threshold`); the reduced gradient
$x_\sigma=|\nabla\rho_\sigma|/\rho_\sigma^{4/3}\sim\rho^{-1/3}$ **diverges in the
tail**, so every GGA enhancement factor is evaluated at an argument growing
without bound exactly where the energy content is least — the origin of the
"pathological wiggles"; and $\rho$ has a Kato cusp at each nucleus while
$|\nabla\rho|$ has cusps at *every extremum* of $\rho$.

### 6.6 libxc's screening, exactly

**[V]** manual: below a density threshold "the library will skip the evaluation
of the functional and instead set all the output quantities to zero", with
per-functional defaults that "might depend on the computer architecture or the
compiler". Setters `xc_func_set_{dens,zeta,sigma,tau}_threshold` all propagate
into auxiliary functionals and **all silently ignore non-positive arguments** —
passing 0 to "disable" screening is a no-op.

**[V]** `functionals.c`, with its own comment:

```c
func->dens_threshold  = func->info->dens_threshold;
/* the density and sigma cutoffs should be connected ... This is the correct scaling */
func->sigma_threshold = pow(func->info->dens_threshold, 4.0/3.0);
func->zeta_threshold  = DBL_EPSILON;
func->tau_threshold   = 1e-20;
```

$\sigma_{\rm thr}$ is **derived** from $\rho_{\rm thr}$ at init. **[D]
Consequence: calling `xc_func_set_dens_threshold` after `xc_func_init` does NOT
update `sigma_threshold`.** Set both explicitly.

**[V]** `work_gga_inc.c`:

```c
dens = (p->nspin == XC_POLARIZED) ? VAR(rho,ip,0) + VAR(rho,ip,1) : VAR(rho,ip,0);
if(dens < p->dens_threshold) return;
my_rho[0]   = m_max(p->dens_threshold, VAR(rho, ip, 0));
my_sigma[0] = m_max(p->sigma_threshold * p->sigma_threshold, VAR(sigma, ip, 0));
if(p->nspin == XC_POLARIZED){
  s_ave = 0.5*(my_sigma[0] + my_sigma[2]);
  /* | grad n |^2 = |grad n_up + grad n_down|^2 > 0 */
  my_sigma[1] = (my_sigma[1] >= -s_ave ? my_sigma[1] : -s_ave);
  my_sigma[1] = (my_sigma[1] <= +s_ave ? my_sigma[1] : +s_ave);
}
```

Six things, **[V]** unless marked:

1. The screen is on the **total** density. A point with $\rho_\alpha=10^{-20}$,
   $\rho_\beta=1$ is evaluated, not skipped.
2. On skip it `return`s **without writing**; outputs are zero only because
   `gga.c` memsets them up front. **[D]** Drive the kernels below that layer and
   you own the zeroing.
3. Above the screen, inputs are **clamped, not rejected**. **[D] So in the
   clamped region the returned `vsigma` is not the derivative of the returned
   `zk` with respect to the `sigma` you passed** — it is the derivative at the
   clamped point. Finite-difference validation must sit well above these floors.
4. The clamp floor is `sigma_threshold * sigma_threshold` — the threshold
   **squared** ($\sim10^{-40}$), not the threshold. Whether intentional is
   unclear; it is what the code does.
5. The cross-spin clamp $|\sigma_{\alpha\beta}|\le(\sigma_{\alpha\alpha}+\sigma_{\beta\beta})/2$
   is **looser than Cauchy–Schwarz** (AM–GM) **[D]**. Its actual purpose, per the
   source comment, is to keep the *total* non-negative — and it does that
   exactly, since $\sigma_{aa}+2\sigma_{ab}+\sigma_{bb}\ge0$ follows from
   $\sigma_{ab}\ge-(\sigma_{aa}+\sigma_{bb})/2$. **[M]** But it protects the total
   only while the *diagonals* are physical; see §6.10.
6. $\sigma_{\alpha\beta}$ has **no lower magnitude floor** — only that two-sided
   clamp. Only same-spin channels get floored.

**meta-GGA additions [V]** (`work_mgga_inc.c`):

```c
/* Many functionals shamelessly divide by tau, so we set a reasonable threshold */
if(p->info->flags & XC_FLAGS_NEEDS_TAU){
  my_tau[0] = m_max(p->tau_threshold, VAR(tau, ip, 0));
  if(p->info->flags & XC_FLAGS_ENFORCE_FHC)
    my_sigma[0] = m_min(my_sigma[0], 8.0*my_rho[0]*my_tau[0]);
}
/* lapl can have any values */
```

**`lapl` is not clamped, floored or checked at all** — while being the noisiest
quantity you can hand it. There is also a documented pointer hazard: when the
functional does not use the Laplacian the caller must pass a genuine `NULL`,
because forming `&VAR(lapl,ip,0)` on a null base yields a bogus *non-*`NULL`
address that the kernel's own guard cannot detect. **Never pass a zero-filled
buffer instead.**

**[V]** `XC_FLAGS_NEEDS_TAU` is **new in 7.0.0**, so a code supporting ≤6.x must
assume every MGGA needs $\tau$. The ChangeLog records these flags still being
*corrected* recently (`MGGA_X_2D_PRP10: added the missing XC_NEEDS_TAU flag`), so
do not treat them as infallible. **[V]** The `xc-threshold` utility helps pick
thresholds per platform — the ChangeLog documents a case where on 32-bit x86 with
the 387 unit the regression suite missed `GGA_X_HJS_PBE`'s `v2sigma2` by twice
its tolerance and `MGGA_X_2D_PRHG07_PRP10`'s `v2rho2` by a factor of 2000.

### 6.7 Why a real-space code needs its own screening on top

**[D]** throughout, reasoning from the **[V]** facts above.

1. **libxc's screening creates a discontinuity, and adaptive bases resolve
   discontinuities by refining.** Exact zeros below the threshold and clamped
   values just above make a kink — and across the skip boundary a jump — on a
   large-area surface out in the tail. A uniform grid merely samples it; an
   adaptive representation *tries to represent it*, refining to the finest level
   in vacuum. The cost appears as tree growth, not as a wrong number, which makes
   it easy to misdiagnose.
2. **The divergence form differentiates a libxc output**, amplifying all of
   point 1. Either apply your own smooth floor to $\rho$ well above
   `dens_threshold`, or use the weak form.
3. **`vsigma` diverges while $\nabla\rho$ vanishes.** For typical GGA exchange
   $\texttt{vsigma}\sim\rho^{-4/3}$; the meaningful object
   $\texttt{vsigma}\cdot\nabla\rho$ is a product of two quantities both heading
   for their noise floors, and it is then *differentiated*. libxc's $10^{-15}$
   floor is far too permissive; real-space codes typically need
   $\rho_{\min}\in[10^{-8},10^{-12}]$, basis- and threshold-dependent.
4. **Clamping breaks the consistency your convergence tests rely on** (§6.6
   item 3): virial checks and finite-difference gradients show small unexplained
   residuals sourced entirely from the tail.
5. **meta-GGA is worse in two independent ways.** $\tau$ is floored only when
   `XC_FLAGS_NEEDS_TAU` is set — absent before 7.0.0 — and `lapl` is not checked
   at all while being a second derivative of a small number.
6. **For response the exposure is squared.** The kernel needs `v2sigma2`, whose
   tail divergence is worse than `vsigma`'s, and $v_{\rm pt}$ takes a divergence
   of a term containing $\texttt{v2sigma2}\cdot\sigma_{\rm pt}\cdot\nabla\rho$ —
   three separately noisy tail quantities multiplied, then differentiated. Screen
   the response at least as aggressively as the ground state, and screen on the
   *ground-state* $\rho$, not on $\rho_{\rm pt}$ (which legitimately has nodes,
   and whose smallness carries no information about conditioning).

### 6.8 Asymptotic correction

The bare LDA/GGA potential decays exponentially rather than as $-1/r$. For an
integral-equation code this is not cosmetic: the bound-state Helmholtz kernel
needs real parameters, hence $\epsilon_{\rm HOMO}\le\mu\le-\epsilon_{\rm HOMO}$,
so orbitals and response functions must actually be bound. Yanai/Harrison/Handy
Eq. (15) **[V]** use the Tozer–Handy form
$v_{xc}^{\rm AC}\to-1/r+I+\epsilon_h$, shifting by $I+\epsilon_h$ (no explicit
virtuals available) and switching to a truncated multipole form outside
$r_A<X\sigma_A$; reported "numerically stable and efficient" **[V]**. Their
remark on the alternative is a pure real-space hazard **[V]**: Hirata's AC "has a
numerical difficulty in computing the Slater potential numerically, which
involves the density function in the denominator" — the same class as the
$1/|\nabla\rho|$ factors of §4.2.

### 6.9 Consistency of derived quantities — a real-space-specific trap

**[M]** This one is not in the literature and cost a full investigation to find
(see §6.10 for the case study). It generalizes beyond XC.

A functional's inputs are usually **not independent**: $\sigma$ is the Gram
matrix of the gradients (§1); $\tau\ge\tau_W=|\nabla\rho|^2/8\rho$; $|\zeta|\le1$.
Each such relation is an *identity* among quantities the host code computes.

In a real-space code with an adaptive basis, a quantity computed as a
**separately projected function** is not pointwise equal to the same quantity
**contracted at the quadrature point** from its constituents. The two agree to
the representation threshold *on average*, and disagree pointwise — worst
wherever the constituents are non-smooth, i.e. at nuclear cusps, where the
discrepancy is O(1) rather than O(thresh).

Consequences, in increasing order of nastiness:

1. An identity that holds analytically fails numerically.
2. A quantity that is a **sum of squares can come out negative**, so any floor
   guarding it (`max(eps, ...)`) fires — silently replacing a physical value by
   an unphysical one, with no diagnostic.
3. A **constraint linking several such quantities is then violated**, because
   one has been floored and its partners have not. This is where a library
   evaluates outside its domain.
4. The library's own guards do **not** save you: they are written assuming the
   inputs are mutually consistent, and typically bound one variable using the
   others (§6.6 item 5).

**The rule.** Derive constrained quantities by **pointwise contraction of the
primitive fields at the quadrature point**, never by reading a second
multiwavelet representation of the same product. Then the constraints hold by
construction, exactly, and every positivity floor becomes a no-op that only fires
where the true value really is below it. This is also cheaper — the products
disappear from the intermediate set entirely.

A corollary for validation: a functional-layer unit test that supplies mutually
consistent inputs **cannot** catch a violation of this kind, because the bug lives
in the host code's construction of the inputs, not in the functional. Only an
end-to-end test on a system whose densities differ per spin will show it (§10).

### 6.10 Case study: the $\sigma$ matrix that was not a Gram matrix

**[M]** MADNESS, branch `pr-implement-mgga`, 2026-08-24. Full account in
`XC-MADNESS-ASSESSMENT.md` §4.5; kept here because the mechanism is generic to
real-space codes and the numbers pin every step.

MADNESS stores $\zeta_\sigma=\nabla\ln\rho_\sigma$ and formed
$\chi_{st}=\zeta_s\!\cdot\!\zeta_t$ as its **own** multiwavelet function, then
built $\nabla\rho_s=\rho_s\zeta_s$ from the $\zeta$ components and
$\sigma_{st}=\rho_s\rho_t\chi_{st}$ from the $\chi$ functions. At the nucleus of
an open-shell Li atom ($\rho_a=6.57$, $\rho_b=6.31$):

| | stored $\chi$ function | contracted from stored $\zeta$ |
|---|---|---|
| $\chi_{aa}$ | **−2.0715** | +0.196488 |
| $\chi_{bb}$ | **−2.4844** | +0.000259 |
| $\chi_{ab}$ | −2.3804 | −0.007138 |

$\chi_{aa}$ and $\chi_{bb}$ are sums of squares. Both diagonals therefore hit the
$10^{-14}$ positivity floor while $\sigma_{ab}$ — correctly unfloored, being
bilinear — kept its full $-98.6$, giving
$\sigma_{aa}+2\sigma_{ab}+\sigma_{bb}=-197$: a negative $|\nabla\rho|^2$.

Measured against libxc alone, sweeping $\sigma_{ab}$ with both diagonals on the
floor: **`exc` stays finite and independent of $\sigma_{ab}$, while `vsigma` of
`GGA_C_PBE` and `MGGA_C_TPSS` turns NaN exactly when the total crosses zero**
($+0.01$ fine, $-0.01$ NaN). The exchange kernels never fail. So the damage is in
the *potential*, not the energy, and it is the **correlation** kernels that break
— which is why plain open-shell PBE and PBE0 aborted alongside TPSS. With
*physical* diagonals a negative total is unreachable: libxc's two-sided clamp
(§6.6 item 5) bounds the total below by zero. A floored diagonal is a
precondition.

Three lessons worth carrying:

- **libxc 7.1.2 returns NaN here; 7.0.0 returns finite garbage.** The two are
  digit-identical on every *admissible* input, and differ only on the
  inadmissible one — which presented as a platform-dependent, seemingly
  nondeterministic bug. A domain error that one version catches and another does
  not is worse than one neither catches.
- **The energy stayed plausible while the potential was destroyed.** At the same
  point `exc` differed in the third digit while $\partial e/\partial\tau_\beta$
  came out $2.49\times10^{-3}$ against a true $1.23$ — wrong by 492×, the
  minority channel losing its whole effective-mass term. Symptom: a believable
  energy next to a diverging minority-spin residual. Any diagnostic that watches
  only the energy is blind to this.
- **A bound derived from admissible inputs tells you nothing here.** A
  1440-point scan had established $|\partial e/\partial\rho_\beta|\le1.04$ for
  "any input" and was used to *rule out* the kernel — but it had sampled only
  positive-semidefinite $\sigma$. Outside the domain there is no bound. When a
  measurement contradicts a bound, check the bound's premises before the
  measurement.

---

## 7. Consequences for a Green's-function / integral-equation code

This is where a MADNESS-class code differs from every published meta-GGA
implementation. The standard iteration replaces the differential eigenproblem by
an integral one:

$$
\psi_i=-2\hat G_{\mu_i}[V\psi_i],\quad\mu_i=\sqrt{-2\epsilon_i},\quad
\hat G_\mu=(-\nabla^2+\mu^2)^{-1},\;\text{kernel }\frac{e^{-\mu|\mathbf r-\mathbf r'|}}{4\pi|\mathbf r-\mathbf r'|}.
$$

Inverting the kinetic operator **analytically** is the entire point: it removes
the ill-conditioning and gives mesh-independent convergence. Exact exchange fits
because $\hat K\psi$ is a *function*; LDA and GGA fit because their potentials
are multiplicative — precisely why the multiwavelet lineage adopted the
divergence form.

**meta-GGA does not.** Formally

$$
\psi_i=-2\hat G_{\mu_i}\Bigl[V_{\rm mult}\psi_i-\tfrac12\sum_xD_x(v_\tau D_x\psi_i)\Bigr]
$$

and the term is cheap to build. The problem is conditioning, not cost.

**[D]** $\hat G_\mu$ smooths by order $-2$; $\nabla\!\cdot\!(v_\tau\nabla\,\cdot)$
differentiates by order $+2$, so the composition is of **order 0** — bounded but
**not compact**, with norm set by $\|v_\tau\|_\infty$ rather than by the mesh or
$\mu$. The fixed-point map is contractive in that channel only while
$\|v_\tau\|_\infty$ is comfortably below 1, and **refining does not help**. In
practice $|v_\tau|$ is small for well-behaved functionals, so this should be
benign — but it needs *measuring*. **[U]** as a literature claim; this analysis
was found in no source read, and it is the thing most worth checking against
multiwavelet convergence theory.

**Options, least to most invasive.**

- **(a) Right-hand side.** Use the equation above directly; monitor
  $\max|v_\tau|$ per iteration. Add damping or KAIN on the $\tau$ channel only if
  the residual stalls. *Start here.*
- **(b) Absorb the mean into $\hat G$.** Write $1+v_\tau=(1+\bar v)(1+\delta)$,
  rescale to $\mu'=\sqrt{-2\epsilon/(1+\bar v)}$ with prefactor $1/(1+\bar v)$,
  leaving only $\delta$ on the right-hand side. Shrinks the order-0 norm to
  $\|\delta\|_\infty$. **[U]**, a construction rather than a citation.
- **(c) Deorbitalize** — replace $\tau$ by $\tau[\rho,\nabla\rho,\nabla^2\rho]$,
  restoring a multiplicative potential (SCAN-L, Mejía-Rodríguez & Trickey
  **[U]**). Least invasive for a GF code, at the cost of a different and less
  accurate functional, and it reintroduces $\nabla^2\rho$ noise.
- **(d) Non-variational.** Evaluate $E_{xc}^{\rm mGGA}$ on GGA orbitals. Yang
  et al. use this as a reference point **[V]**: "the gKS total energy should
  always be lower than the OEP total energy, and both should be lower than the
  meta-GGA total energy evaluated with non-variational orbitals." Defensible for
  energetics; wrong for forces and anything needing self-consistency.

**Forces.** Once the operator is non-multiplicative the Pulay-type terms change.
ONETEP (Womack, Mardirossian, Head-Gordon & Skylaris, JCP **145**, 204114 (2016))
is **[C]** the closest published analogue to the multiwavelet problem and the
first thing to read.

**State of the art.** **[V]** implementations from papers read: APW/all-electron
(Doumont/Tran/Blaha), periodic OEP+gKS in BAND (Yang et al.), GAMESS (Zahariev
et al.), FLAPW $\tau$ (Ye). **[C]** ONETEP. **[U]** GPAW/Octopus/BigDFT status
unchecked. **No published meta-GGA implementation in any multiwavelet code was
located**, and the search was not thorough enough to assert absence — state it as
"none located", not "none exists". Worth checking whether MRChem's XCFun-based
layer exposes $\tau$.

---

## 8. Linear response — TDDFT / CPKS

### 8.1 Unpolarized GGA kernel, and a factor-of-2 hazard

**[D]** Perturb $\rho\to\rho+\lambda\rho_{\rm pt}$ and differentiate
$v_{xc}=\varepsilon_\rho-\nabla\!\cdot\!(2\varepsilon_\sigma\nabla\rho)$ at
$\lambda=0$. With $\sigma_{\rm pt}=2\nabla\rho\!\cdot\!\nabla\rho_{\rm pt}$ (the
actual first-order change of libxc's $\sigma$):

$$
\boxed{\begin{aligned}
v_{\rm pt}=\;&\texttt{v2rho2}\,\rho_{\rm pt}+\texttt{v2rhosigma}\,\sigma_{\rm pt}\\
&-\nabla\!\cdot\!\bigl[2\,\texttt{vsigma}\,\nabla\rho_{\rm pt}+2\,\texttt{v2rhosigma}\,\rho_{\rm pt}\nabla\rho+2\,\texttt{v2sigma2}\,\sigma_{\rm pt}\nabla\rho\bigr]
\end{aligned}}
$$

> ### ⚠ The hazard
>
> The coefficient pattern $\{2,2,2,4\}$ is correct **only** under
> $\sigma_{\rm pt}=\nabla\rho\!\cdot\!\nabla\rho_{\rm pt}$, *without* the 2.
> Pairing it with $\sigma_{\rm pt}=2\nabla\rho\!\cdot\!\nabla\rho_{\rm pt}$
> double-counts in **exactly two of five terms** — the local
> $\varepsilon_{\rho\sigma}\sigma_{\rm pt}$ (2 instead of 1) and the
> $\varepsilon_{\sigma\sigma}\sigma_{\rm pt}\nabla\rho$ divergence (4 instead of
> 2). Two of five wrong by 2× is the worst failure mode: excitation energies
> shift by a few tenths of an eV rather than blowing up, so it survives casual
> testing.
>
> **Use the form above.** Then $\sigma_{\rm pt}$ *means* "the first-order change
> in the quantity libxc calls `sigma`", the same variable feeds `v2rhosigma` and
> `v2sigma2`, and there is exactly one place where the factor 2 lives.

**Sanity checks [D].** Pure LDA gives $v_{\rm pt}=\texttt{v2rho2}\,\rho_{\rm pt}$
— the textbook ALDA kernel ✓. For $\varepsilon=a(\rho)+b(\rho)\sigma$ it
reproduces the direct derivative of that model's $v_{xc}$ ✓.

### 8.2 Spin-polarized GGA kernel

**[D]** First-order changes — note the cross channel gets **two** terms with
coefficient 1 each, the spin-polarized image of the §4.1 asymmetry and the second
place a factor of 2 goes missing:

$$
\sigma_{{\rm pt},\alpha\alpha}=2\nabla\rho_\alpha\!\cdot\!\nabla\rho_{{\rm pt},\alpha},\qquad
\sigma_{{\rm pt},\alpha\beta}=\nabla\rho_\alpha\!\cdot\!\nabla\rho_{{\rm pt},\beta}+\nabla\rho_\beta\!\cdot\!\nabla\rho_{{\rm pt},\alpha}
$$

$$
\begin{aligned}
v_{{\rm pt},\alpha}=\;&\sum_\tau\texttt{v2rho2}[\alpha,\tau]\rho_{{\rm pt},\tau}+\sum_{st}\texttt{v2rhosigma}[\alpha,st]\sigma_{{\rm pt},st}\\
-\nabla\!\cdot\!\Bigl[\;&2\texttt{vsigma}_{\alpha\alpha}\nabla\rho_{{\rm pt},\alpha}+\texttt{vsigma}_{\alpha\beta}\nabla\rho_{{\rm pt},\beta}\\
&+\Bigl(2\sum_\tau\texttt{v2rhosigma}[\tau,\alpha\alpha]\rho_{{\rm pt},\tau}+2\sum_{st}\texttt{v2sigma2}[\alpha\alpha,st]\sigma_{{\rm pt},st}\Bigr)\nabla\rho_\alpha\\
&+\Bigl(\sum_\tau\texttt{v2rhosigma}[\tau,\alpha\beta]\rho_{{\rm pt},\tau}+\sum_{st}\texttt{v2sigma2}[\alpha\beta,st]\sigma_{{\rm pt},st}\Bigr)\nabla\rho_\beta\Bigr]
\end{aligned}
$$

and $\beta$ by $\alpha\!\leftrightarrow\!\beta$,
$\alpha\alpha\!\leftrightarrow\!\beta\beta$ ($\alpha\beta$ maps to itself).
Indexing: $\sum_\tau$ over $\{\alpha,\beta\}$, $\sum_{st}$ over the **three**
sigma channels; `v2rhosigma` is the $2\times3$ rectangle (offset $3\tau+st$)
**[D]**; **`v2sigma2` is triangle-packed — unpack it into a symmetric $3\times3$
before contracting**, or the off-diagonals get counted once instead of twice.
That is the third factor-of-2 trap.

### 8.3 meta-GGA kernel

Zahariev/Leang/Gordon claim priority **[V]** ("presented here for the first
time"). Their Eqs. (27)–(28) replace the four-index $f_{xc}$ contraction with one
over $[\nabla\psi_p][\nabla\psi_q]\,\delta^2E_{xc}/\delta\tau\delta\tau\,[\nabla\psi_s][\nabla\psi_t]$;
third order follows with prefactor $1/8$ and three $\nabla\psi$ pairs **[V]**.
Their singlet/triplet equations (A1)/(A2) carry the prefactor pattern
$\tfrac14/\tfrac12/1$ — one $\tfrac12$ per $\tau$-derivative, with the
$\gamma$-derivative's 2 cancelling one — which is the internal consistency check
that pins the $\tfrac12$ and exposes the main text's stray 2 (§5.5). Triplet
takes the spin-cross terms with a **minus** sign.

**[D]** Structurally the $\tau$ sector contributes both a local term (through
`v2rhotau`, `v2sigmatau`) and a *perturbation of the non-multiplicative
operator*, $-\tfrac12\nabla\!\cdot\!(v_{\tau,{\rm pt}}\nabla\psi_i)$ with
$v_{\tau,{\rm pt}}=\sum\texttt{v2rhotau}\rho_{\rm pt}+\sum\texttt{v2sigmatau}\sigma_{\rm pt}+\sum\texttt{v2tau2}\tau_{\rm pt}$.

**Response-side grid warning [V]** (Yang et al.): "TDDFT with meta-GGA xc kernels
is also expected to show this grid sensitivity." The kernel needs
$\partial^2e/\partial\tau^2$ — *second* derivatives of functions whose *first*
derivatives already oscillate (§6.3) — so expect a materially finer grid than the
ground state.

### 8.4 How the multiwavelet lineage does it

Yanai/Harrison/Handy compute the second derivatives fully numerically. Verbatim
**[V]**:

> The numerical computation is carried out literally according to the arithmetic
> expressions of equations (13) and (14), which involves the multiplications of
> functions and the applications of the differentiation operator to functions. …
> Therefore, the computational cost required is O(N).

That is the real-space idiom in one sentence: **never form an integral kernel;
multiply functions and apply the derivative operator pointwise in 3D.**

> **Missing reference.** For the GGA implementation details their text defers to
> "T. Yanai, R. J. Harrison, *J. Chem. Phys.*, submitted (2004)", which does not
> appear to have been published. If so, the detailed multiwavelet GGA write-up
> does not exist in the literature and **the code is the source of truth.**

---

## 9. Factor-of-2, packing and consistency traps, collected

| # | Location | Correct | Wrong-way failure | Status |
|---|---|---|---|---|
| 1 | $\sigma_{\alpha\beta}$ definition | $\nabla\rho_\alpha\!\cdot\!\nabla\rho_\beta$, **no** 2 | 2× in all cross-spin terms | **[V]** manual + codegen |
| 2 | GGA potential, same-spin | $2\,\texttt{vsigma}_{\alpha\alpha}\nabla\rho_\alpha$ | 2× too small | **[D]** §4.1 |
| 3 | GGA potential, cross-spin | $1\times\texttt{vsigma}_{\alpha\beta}\nabla\rho_\beta$ | 2× too large | **[D]** §4.1 |
| 4 | `zk` vs $\varepsilon$ | $\varepsilon=\rho\,\texttt{zk}$; all `v*` differentiate $\varepsilon$ | energy/potential inconsistent by $\rho$ | **[V]** manual |
| 5 | unpol → pol mapping | $\rho/2$, $\sigma/4$; $\texttt{vsigma}_{\rm unpol}=\tfrac14(2\texttt{vsigma}_{\alpha\alpha}+\texttt{vsigma}_{\alpha\beta})$ | closed-shell `nspin=1` and `2` disagree | **[V]**+**[D]** §2.5 |
| 6 | $\tau$ normalization | $\tfrac12\sum|\nabla\phi|^2$; operator $-\tfrac12\nabla\!\cdot\!(v_\tau\nabla\psi)$ | 2× in the $\tau$ term | **[V]** manual + FHC clamp |
| 7 | $\sigma_{\rm pt}$ in the kernel | one convention throughout, never mixed | 2× in the `v2rhosigma` local term and the `v2sigma2` div term *only* | **[D]** §8.1 |
| 8 | `v2sigma2` packing | triangle; unpack before contracting | off-diagonals counted once | **[V]** |
| 9 | `v2rhosigma` packing | full $2\times3$ rectangle | wrong element; same size, no crash | **[V]** |
| 10 | `lapl` when unused | pass genuine `NULL` | bogus non-`NULL` address defeats the guard | **[V]** source comment |
| 11 | `sigma_threshold` after init | set explicitly; derived from `dens_threshold` at init only | screening silently inconsistent | **[V]** |
| 12 | **$\sigma$ as a Gram matrix** | contract $\sigma_{st}$ pointwise from the *same* gradients handed to libxc | diagonals hit their floor, cross term does not, total goes negative → correlation `vsigma` NaN (7.1.2) or finite garbage (7.0.0) | **[M]** §6.9–6.10 |

---

## 10. Test strategy

**The single highest-value test.** Run a closed-shell system through both the
`nspin=1` and `nspin=2` code paths and require agreement to roundoff for
**energy, potential *and* kernel**. Rows 1, 2, 3, 5, 8, 9 of §9 all fail it.

**But it does not catch row 12** — and this is worth stating explicitly, because
it was the blind spot that let the §6.10 defect survive. A spin-compensated
system run through the polarized path has $\rho_\alpha=\rho_\beta$, so every
cross term is degenerate with its diagonal and every consistency violation
cancels. **Add an end-to-end test on a system whose two spin densities genuinely
differ** — the cheapest is a lithium atom — and assert on the *potential*, not
only the energy (§6.10: the energy stayed right to three digits while the
potential was wrong by 492×).

**Then** finite-difference the kernel against the potential
($v_{\rm pt}\stackrel{?}{=}[v_{xc}(\rho+h\rho_{\rm pt})-v_{xc}(\rho-h\rho_{\rm pt})]/2h$)
for row 7, and the potential against the energy (or a virial check) for rows 4
and 6. **Do all of it well above the screening thresholds of §6.6** — inside the
clamped region `vsigma` is not the derivative of the returned `zk`, so
finite-difference tests there report artifacts.

**Higher-rung specifics:**

- Dirac exchange closed-shell reduction $v_x=-(3/\pi)^{1/3}\rho^{1/3}$ (§3).
- $e_{xc}=\tau$ must reproduce $-\tfrac12\nabla^2\psi$ (§5.3) — pins the
  $\tfrac12$ and the libxc $\tau$ convention in one test.
- $\int\tau\,d^3r$ must equal the kinetic energy from the code's own operator.
- The operator must satisfy its own weak form,
  $\langle\phi|v_\tau|\psi\rangle=\tfrac12\int v_\tau\nabla\phi\!\cdot\!\nabla\psi$
  — the integration by parts *and* the self-adjointness the Fock build relies on.
- Finite-difference $\partial e/\partial\tau$ against $e_{xc}$ **in the polarized
  path specifically**: it uses a different libxc entry point, a different result
  index and a different stride from the closed-shell one, so an index/stride/spin
  error there is invisible to closed-shell tests.
- Log $\min/\max v_\tau$ and check $1+\min v_\tau>0$ (§5.4), remembering the
  cusp-sampling caveat.
- Assert the Gram-matrix relations on the $\sigma$ actually handed to libxc, or
  better, construct it so they cannot fail (§6.9).

---

## 11. Condensed meta-GGA recipe for a Green's-function code

```
Per SCF iteration, per spin sigma:
  1.  g_i     = grad psi_i                        (reuse in step 7)
  2.  tau     = (1/2) sum_i w_i |g_i|^2           # the 1/2 -- libxc >= 2.0.0
  3.  rho_s;  sigma_{st} CONTRACTED POINTWISE from the same grad rho_s
      handed to libxc, so the Gram relations hold exactly (§6.9)
      lapl ONLY if XC_FLAGS_NEEDS_LAPLACIAN; else pass a genuine NULL
  4.  libxc: zk (per particle!), vrho, vsigma, vtau (per volume);
      prefer the fused xc_{gga,mgga}_exc_vxc
  5.  V_mult  = vrho - div( 2*vsigma_ss grad rho_s + vsigma_ss' grad rho_s' )
                [+ lap(vlapl) only at Laplacian level]
      v_tau   = vtau
  6.  log min/max v_tau ; check 1 + min v_tau > 0
  7.  W psi_i = V_mult psi_i - (1/2) sum_x D_x( v_tau * D_x psi_i )
      # NEVER form grad(v_tau)
      # with psi = R F, instead:  -(1/2)( div W - U1.W ),
      #   W = v_tau (grad F - U1 F)      -- the R factors cancel (6.2)
  8.  psi_i <- -2 G_{mu_i}[ W psi_i + (other potentials) psi_i ]
      optionally rescale G by a mean (1 + v_bar) to shrink the order-0 term
  9.  Functional order of attack: TPSS (z-based, no alpha=1 oscillation),
      then r2SCAN, and only then SCAN.
```

Expect 2–3× the cost of a GGA iteration **[V]** and a materially tighter
threshold requirement, concentrated where $\alpha$ crosses 1 **[V]**. Measure
$\|\partial e_{xc}/\partial\tau\|_\infty$ early — it decides between options (a),
(b) and (c) of §7.

---

## 12. Bibliography

**Read in full and quoted verbatim — [V]**

- libxc **5.1.x** manual — <https://libxc.gitlab.io/manual/libxc-5.1.x/> — and
  **3.0.x** — <https://libxc.gitlab.io/manual/previous/libxc-3.0.x/>
- libxc `devel` source — `src/xc.h`, `functionals.c`, `work_gga_inc.c`,
  `work_mgga_inc.c`, `gga.c`, `gga_x_pbe.c`, `lda_x.c`,
  `scripts/maple2c_lib/gga.py`, `ChangeLog.md` —
  <https://gitlab.com/libxc/libxc/-/tree/devel>
- **Pople, Gill & Johnson**, Chem. Phys. Lett. **199**, 557 (1992) —
  <http://rsc.anu.edu.au/~pgill/papers/031KSLCAO.pdf>. Eqs. (7), (8) and prose
  read directly; Eq. (9) is OCR-degraded and (16) a bitmap, both verified through
  Lehtola et al.'s verbatim reproduction.
- **Lehtola, Blockhuys & Van Alsenoy**, arXiv:1912.12029, Eqs. (46)–(55).
- **Balbás, Martins & Soler**, arXiv:cond-mat/0104171 (published PRB **64**,
  165110 (2001); the *journal* reference is itself **[U]** — the preprint was
  read). Eqs. (1)–(6), §I, §V.
- **Mortensen, Hansen & Jacobsen**, PRB **71**, 035109 (2005),
  arXiv:cond-mat/0411218. §VI, Eq. (72), footnote 37.
- **Yanai, Harrison & Handy**, Mol. Phys. **103**, 413 (2005),
  doi:10.1080/00268970412331319236, OA at <https://zenodo.org/record/1234395>.
  Eqs. (10)–(16), §2.2–2.3. **The primary reference for LDA/GGA in a multiwavelet
  basis.**
- **Yang, Peng, Sun & Perdew**, PRB **93**, 205205 (2016), arXiv:1603.00512.
  Eqs. (1)–(7), (16), §III.
- **Zahariev, Leang & Gordon**, JCP **138**, 244108 (2013),
  doi:10.1063/1.4811270. Eqs. (1)–(3), (8), (17)–(30), (36), Appendix.
- **Doumont, Tran & Blaha**, arXiv:2112.13145. Eqs. (5), (6), (A.1), Table 1.
- **Furness, Kaplan, Ning, Perdew & Sun** [r²SCAN], JPCL **11**, 8208 (2020),
  arXiv:2008.03374. Abstract and numerical-motivation statements.
- **Ye**, PRB **91**, 075101 (2015), arXiv:1501.01390. **Abstract only.**

**Citation and purpose confirmed by a source that was read — [C]**

Seidl/Görling/Vogl/Majewski/Levy PRB **53**, 3764 (1996) [generalized KS];
**Womack/Mardirossian/Head-Gordon/Skylaris JCP 145, 204114 (2016)** [ONETEP
meta-GGA — the closest published analogue to the multiwavelet problem];
Neumann/Nobes/Handy Mol. Phys. **87**, 1 (1996) and Neumann/Handy CPL **266**, 16
(1997) [the NNH scheme]; Adamo/Ernzerhof/Scuseria JCP **112**, 2643 (2000);
Arbuznikov/Kaupp CPL **381**, 495 (2003), JCP **131**, 084103 (2009), PCCP **4**,
5467 (2002); Gritsenko/Baerends PRA **64**, 042506 (2001);
Grüning/Gritsenko/Baerends JCP **116**, 6435 (2002); Ayers/Parr/Nagy IJQC **90**,
309 (2002) [the $\tau/\tau_L$ gauge identity]; Eich/Hellgren JCP **141**, 224107
(2014); Leang/Zahariev/Gordon JCP **136**, 104101 (2012); JCP **127**, 214103
(2007) *Avoiding singularity problems associated with meta-GGA … kinetic energy
density* (**authors not verified**; the most on-point older reference for
$\tau$-related singularities — worth chasing); JCP **137**, 164105 (2012);
White & Bird PRB **50**, 4954 (1994) (known only through Balbás et al. and
GPAW's citation); Hamann PRB **54**, 1568 (1996); Bylander & Kleinman
J. Comput. Phys. **136**, 599 (1997); Johnson/Gill/Pople JCP **98**, 5612 (1993)
(abstract only); Lehtola/Steigemann/Oliveira/Marques SoftwareX **7**, 1 (2018)
(the citation libxc requests; every OA mirror tried returned 403/404, so all §2
conventions are anchored to the manual and source instead).

**Standard citations, not verified — [U]**

SCAN (Sun/Ruzsinszky/Perdew PRL **115**, 036402 (2015)); TPSS (Tao/Perdew/
Staroverov/Scuseria PRL **91**, 146401 (2003)); rSCAN (Bartók & Yates JCP **150**,
161101 (2019)); Laplacian level (Perdew & Constantin PRB **75**, 155109 (2007));
PW92; VWN; Furche & Ahlrichs JCP **117**, 7433 (2002) **+ erratum JCP 121, 12772
(2004)**; Van Caillie & Amos CPL **291**, 71 (1998), **317**, 159 (2000);
Lehtola & Marques JCP **157**, 174114 (2022); the MADNESS multiwavelet papers
(Harrison/Fann/Yanai/Gan/Beylkin JCP **121**, 11587 (2004); Yanai et al. JCP
**121**, 2866 and 6680 (2004); Alpert/Beylkin/Gines/Vozovoi JCP **182**, 149
(2002)); GPAW (Enkovaara et al. JPCM **22**, 253202 (2010)); Octopus; BigDFT;
Mejía-Rodríguez & Trickey [SCAN-L]; MRChem JCTC 2023.

---

## 13. Known gaps

1. **§8's coefficients rest on derivation from verified libxc definitions, not on
   a fetched paper.** Furche & Ahlrichs (JCP **117**, 7433 *and its erratum*) is
   the cross-check that did not happen. When obtaining it, check (a) whether it
   uses $\sigma=|\nabla\rho|^2$, $g=|\nabla\rho|$ or the reduced
   $x=|\nabla\rho|/\rho^{4/3}$ — each rescales the $\sigma$-derivative
   coefficients differently — and (b) whether its $\sigma_{\alpha\beta}$ carries
   a factor 2, which several older papers do and libxc does not.
2. **No published meta-GGA implementation in a multiwavelet code was located**,
   and the search was not exhaustive enough to assert absence.
3. **The order-0 / non-compactness analysis in §7 is unsupported by any source
   read** — the claim most worth checking against multiwavelet convergence theory.
4. **The GPAW / Octopus / BigDFT meta-GGA survey did not happen.**
5. **PW92 and VWN were not read from the originals**; their formulas are omitted
   here rather than reproduced from memory. Verify against libxc's maple sources.
6. **Eggbox amplitude for meta-GGA is unquantified.** The verified grid data
   ($G_{\max}$ 20–24 vs 12; 120→266→504 radial points) concerns integration
   grids, not eggbox.
7. **The `v2*` meta-GGA packing orderings in §2.3 are [D]**, not [V]. Confirm by
   finite difference.
8. **§6.10's "correlation kernels NaN on a negative total $\sigma$" is [M], not
   traced to libxc source.** The behaviour is reproducible in ~40 lines against
   libxc alone and is version-dependent (7.1.2 NaN, 7.0.0 finite garbage); *which*
   internal expression fails was not identified.
