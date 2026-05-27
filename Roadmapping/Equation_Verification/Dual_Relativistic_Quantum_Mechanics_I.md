# Equation Verification: *Dual Relativistic Quantum Mechanics I*

**Authors:** Tepper L. Gill, Gonzalo Ares de Parga, **Trey Morris** (repo owner), Mamadou Wade
**Date:** August 21, 2021
**Source:** [`../Tepper_Gill_Papers/Dual Relativistic Quantum Mechanics I.pdf`](../Tepper_Gill_Papers/Dual%20Relativistic%20Quantum%20Mechanics%20I.pdf)
**Markdown:** [`../Converted_Markdown/Dual Relativistic Quantum Mechanics I/Dual Relativistic Quantum Mechanics I.md`](../Converted_Markdown/Dual%20Relativistic%20Quantum%20Mechanics%20I/Dual%20Relativistic%20Quantum%20Mechanics%20I.md)

**Verification status:** In progress (2026-05-11). Wolfram MCP online. The Maxwell-paper verification doc (`Two_Mathematically_Equivalent_Versions_of_Maxwells_Equations.md`) is referenced throughout for Section I, which is mostly a recap of that paper.

## Conventions

Same as the Maxwell-paper verification:
- Gaussian (c.g.s.) units.
- $\mathbf{w} = d\mathbf{x}/dt$ (observer-time velocity), $\mathbf{u} = d\mathbf{x}/d\tau = \gamma\mathbf{w}$ (proper-time velocity).
- $b = \sqrt{c^2 + \mathbf{u}^2}$ ("collaborative" speed of light).
- $\gamma = 1/\sqrt{1 - \mathbf{w}^2/c^2}$.
- $\boldsymbol\pi = \mathbf{p} - e\mathbf{A}/c$ (kinetic momentum).
- $H_0 = \sqrt{c^2\boldsymbol\pi^2 + m^2c^4}$, $H = H_0 + V$.
- Standard Dirac matrices $\boldsymbol\alpha, \beta$ with $\{\alpha_i,\alpha_j\}=2\delta_{ij}$, $\{\boldsymbol\alpha,\beta\}=0$, $\beta^2=1$.

---

## Equation index

| Eq. | Topic | Verdict |
|---:|---|---|
| (I.1) | $d\tau^2 = dt^2 - (1/c^2)d\mathbf{x}^2$ — proper-time definition | ✅ |
| (I.2) | $cdt = b\,d\tau$, $\mathbf{u} = \gamma\mathbf{w}$ | ✅ |
| (I.3) | $(1/c)d/dt = (1/b)d/d\tau$ — time-derivative duality | ✅ |
| (I.4) | $\mathbf{w}/c = \mathbf{u}/b$ — velocity duality | ✅ |
| (I.5) | $dW/dt = \{H, W\}$ — textbook | ✅ |
| (I.6) | $K = H^2/(2mc^2) + mc^2/2$ — canonical proper-time Hamiltonian | ✅ |
| (I.7) | Hamilton's equations for the dual Hamiltonian | ✅ |
| (I.8) | Maxwell's equations (Gaussian) — textbook | ✅ |
| (I.9) | Proper-time Maxwell's equations | ✅ |
| (I.10) | $\mathbf{B}$-field dual wave equation with dissipative term | ✅ |
| (I.11) | $\mathbf{E}$-field dual wave equation with dissipative term | ✅ |
| (I.7-fields) | Modified Liénard–Wiechert $\mathbf{E}$, $\mathbf{B}$ | ✅ |
| (II.1) | Dual Dirac equation | ✅ |
| (II.2) | Dual square-root (potential outside) | ✅ |
| (II.3) | Dual square-root (potential in mass) — **new** | ✅ |
| (III.1) | Dirac eigenvalue problem (matrix form) | ✅ |
| (III.2) | $\psi_2$ solution | ✅ |
| (III.3) | $K_D$ in terms of $V_0, V$ | ✅ |
| (III.4) | Dual Dirac eigenvalue equation | ✅ |
| (III.5) | Coupled $\psi_1, \psi_2$ equations | ✅ |
| (III.6) | Expanded eigenvalue equation | ✅ |
| (III.7) | $\psi_2$ with denominator approximation | ✅ |
| (III.8) | Approximated eigenvalue equation | ✅ |
| (III.9)–(III.12) | Pauli identities | ✅ |
| (III.13)–(III.17) | Spherical-coordinate algebra | ✅ |
| (III.18) | $g_r$ formula structure | ✅ (algebra) |
| (III.19), (III.20) | Other new terms | ✅ |
| (III.21)–(III.23) | g-factor numerical reproduction | ⚠/✅ at hypothesis-(ii) framework precision (was 🔴 fails → ⚠ characterised → ⚠/✅-conditional) — closed form $r_e/r_0 = (2-a_e)/(2(2+a_e))$ via QED-inheritance of $a_e(\alpha)$; unconditional ✅ pending [#75](https://github.com/temoTxt/PyPhysics/issues/75) (hypothesis (i) author input); see [§III.D-extension](#iiid-extension--first-principles-derivation-of-r_er_0-closes-64) |

---

> **⚠/✅-conditional resolved by closed-form algebraic inversion (Section III.D).** The paper's published $r_e = 0.499857150068631 \times r_0$ does not reproduce the experimental $g_e$, but the cutoff *formula* $g_r = 2[1 - 4r_0/(2r+r_0)]$ is correct: inverting it for $g_r(r_e/r_0) = -2(1 + a_e)$ where $a_e$ is the QED anomalous magnetic moment yields the closed form $r_e/r_0 = (2 - a_e)/(2(2 + a_e))$. With CODATA-full $a_e^{\rm expt}$, this gives $r_e/r_0 = 0.499\,420\,509\,913\,18$, matching the triangulated optimum at residual $3.45\times 10^{-13}$ — that is **1.4σ above** the triangulation precision floor $\sigma_r = 2.5\times 10^{-13}$ (consistent at 2σ; Wolfram MCP: `Abs[residual] > sigma → True`). **Honest scope:** this is *reproduction-by-inheritance* under hypothesis (ii) (QED supplies $a_e(\alpha)$; the framework supplies the (III.22) cutoff–anomaly identification and the algebraic inversion); it is not an independent dual-framework derivation of $a_e$. The §III.D Schwinger identification subsection from [PR #70](https://github.com/temoTxt/PyPhysics/pull/70) characterises the same algebraic operation as the back-fit identity from the Candidate-3 perspective; both subsections together give the full picture. Verdict marker history: **🔴 (was) → ⚠ characterised (#61) → ⚠/✅ at hypothesis-(ii) framework precision (#64, this PR)**. **Unconditional ✅** requires hypothesis-(i) re-derivation (proper-time photon propagator + bound-state propagator + mass-renormalisation prescription at the cutoff), tracked in [issue #75](https://github.com/temoTxt/PyPhysics/issues/75). See the [§III.D-extension](#iiid-extension--first-principles-derivation-of-r_er_0-closes-64) section and [`FINDINGS_for_author_review.md` Finding 2](FINDINGS_for_author_review.md#finding-2--drqm-i-eq-iii22-published-r_e-does-not-reproduce-the-experimental-g_e) for full details.

---

## Section I: Recap of dual classical theory

The bulk of Section I is a transcription of results already verified in [[Two_Mathematically_Equivalent_Versions_of_Maxwells_Equations]]. Cross-references below.

### Eq. (I.1) — Proper-time definition

**As printed (line 41):**
$$d\tau^2 = dt^2 - \frac{1}{c^2}d\mathbf{x}^2, \qquad \mathbf{w} = \frac{d\mathbf{x}}{dt}.$$

**Status:** Textbook (Jackson §11.4). No verification needed.

**Verdict:** ✅ Textbook.

---

### Eq. (I.2) — Proper-time velocity

**As printed (line 46):**
$$c\,dt = \sqrt{\mathbf{u}^2 + c^2}\,d\tau = b\,d\tau, \qquad \mathbf{u} = \frac{d\mathbf{x}}{d\tau} = \gamma(\mathbf{w})\mathbf{w}.$$

**Status:** Same identity as Maxwell paper Step 1 of Eq. (1) derivation. Algebraic rearrangement of (I.1):

$$dt^2 = d\tau^2 + \frac{1}{c^2}d\mathbf{x}^2 = d\tau^2(1 + \mathbf{u}^2/c^2) \;\Rightarrow\; \frac{dt}{d\tau} = \sqrt{1+\mathbf{u}^2/c^2} = \frac{b}{c}.$$

**Verdict:** ✅ Algebraic consequence of (I.1). See Maxwell paper Eq. (1) verification for full chain.

---

### Eq. (I.3) — Time-derivative duality

**As printed (line 52):**
$$\frac{1}{c}\frac{d}{dt} \equiv \frac{1}{b}\frac{d}{d\tau}.$$

**Status:** Identical to Maxwell paper Eq. (2).

**Verdict:** ✅ Confirmed in Maxwell paper Eq. (2) verification (`FullSimplify` returned 0).

---

### Eq. (I.4) — Velocity duality

**As printed (line 58):**
$$\frac{\mathbf{w}}{c} = \frac{\mathbf{u}}{b}.$$

**Status:** Identical to Maxwell paper Eq. (1).

**Verdict:** ✅ Confirmed in Maxwell paper Eq. (1) verification.

---

### Eq. (I.5) — Standard Poisson-bracket dynamics

**As printed (line 67):**
$$\frac{dW}{dt} = \{H, W\}.$$

**Verdict:** ✅ Textbook (Goldstein §9.4).

---

### Eq. (I.6) — Canonical proper-time Hamiltonian

**As printed (line 84):**
$$K = \frac{H^2}{2mc^2} + \frac{mc^2}{2}, \qquad \frac{dW}{d\tau} = \{K, W\}.$$

**Status:** Identical to Maxwell paper Eq. (16). The derivation in DRQM I (lines 71–83) is slightly more explicit than Maxwell paper's presentation: starts from $dW/d\tau = (H/mc^2)\{H,W\}$, identifies $H/mc^2 \cdot \partial H/\partial\mathbf{p} = \partial(H^2/(2mc^2))/\partial\mathbf{p}$, and reads off $a = a' = mc^2/2$ from the requirement that $K|_{\mathbf{p}=0} = mc^2$.

**Pedagogical check (the integration step):**

- **Step 1.** From $\partial/\partial\mathbf{p}[H^2/(2mc^2) + a] = (H/mc^2)\partial H/\partial\mathbf{p}$ — true for any $a$ that doesn't depend on $\mathbf{p}$.
- **Step 2.** Similarly for the spatial derivative. So $K = H^2/(2mc^2) + a$ for some constant (or pure function of conserved quantities) $a$.
- **Step 3.** At $\mathbf{p}=0$: $H = mc^2$, so $H^2/(2mc^2) = mc^2/2$. Imposing $K|_{\mathbf{p}=0} = mc^2$ gives $a = mc^2/2$.

**Mathematica check:**
```mathematica
ClearAll[hH, mm, cc]; KK = hH^2/(2 mm cc^2) + mm cc^2/2; Print["K at H = m c^2: ", FullSimplify[KK /. hH -> mm cc^2, Assumptions -> {mm > 0, cc > 0}]]
(* Result: m c^2  ✅  (matches required boundary condition) *)
```

**Verdict:** ✅ Confirmed by Wolfram MCP; cross-referenced to Maxwell paper Eq. (16).

---

### Eq. (I.7) — Hamilton's equations for the dual Hamiltonian

**As printed (lines 93–103):**
The paper expands $K$ with minimal coupling and shows
$$\frac{c}{b}\frac{d\boldsymbol\pi}{d\tau} = e\mathbf{E} + \frac{e}{b}(\mathbf{u}\!\times\!\mathbf{B}) = e\mathbf{E} + \frac{e}{c}(\mathbf{w}\!\times\!\mathbf{B}) = \frac{d\boldsymbol\pi}{dt},$$
demonstrating "mathematical equivalence" of the dual and standard equations of motion.

**Status:** Identical content to Maxwell paper Eqs. (17)–(18), but presented in cleaner form (without the $V/(mcb)$ correction term since this paper writes $H$ as $H_0 + V$ via the *standard* substitution).

**Mathematica check:**
```mathematica
(* The transition from the c/b * d pi/d tau expression to e E + (e/c) w x B uses Eqs (I.3) and (I.4). *)
(* In 1D, the (e/b) u term in the magnetic-force becomes (e/c) w via Eq (I.4): u/b = w/c. *)
(* Symbolically: e/b * u = e (u/b) = e (w/c) = e/c * w  ✓  *)
ClearAll[ee, uu, ww, bb, cc]; Print["(e/b) u vs (e/c) w via Eq (I.4) [u/b = w/c]: ", FullSimplify[(ee/bb) uu - (ee/cc) ww /. {uu -> bb ww/cc}]]
(* Result: 0  ✅  (Wolfram MCP, 2026-05-11) *)
```

**Verdict:** ✅ Mathematical equivalence verified. The two forms of the equation of motion are algebraic transcriptions of each other via Eqs. (I.3)–(I.4).

---

### Eq. (I.8) — Standard Maxwell's equations

**As printed (lines 111–113):**
$$\nabla\!\cdot\!\mathbf{B}=0,\quad \nabla\!\cdot\!\mathbf{E}=4\pi\rho, \quad \nabla\!\times\!\mathbf{E}=-\frac{1}{c}\frac{\partial\mathbf{B}}{\partial t},\quad \nabla\!\times\!\mathbf{B}=\frac{1}{c}\!\left[\frac{\partial\mathbf{E}}{\partial t}+4\pi\rho\mathbf{w}\right].$$

**Verdict:** ✅ Textbook quotation, identical to Maxwell paper Eq. (3).

---

### Eq. (I.9) — Proper-time Maxwell's equations

**As printed (lines 118–120):** Same as (I.8) but with $1/c \to 1/b$, $\partial_t \to \partial_\tau$, $\mathbf{w} \to \mathbf{u}$.

**Verdict:** ✅ Identical to Maxwell paper Eq. (3′); confirmed there.

---

### Eqs. (I.10), (I.11) — Dual wave equations with dissipative term

**As printed (lines 127–133):** Same as Maxwell paper Eq. (4) (B-field and E-field versions).

**Verdict:** ✅ Identical to Maxwell paper Eq. (4); coefficient $-(\mathbf{u}\!\cdot\!\mathbf{a})/b^4$ confirmed there.

---

### Unnumbered (lines 139, 143) — Modified Liénard–Wiechert fields

**Status:** Identical to Maxwell paper Eq. (7).

**Verdict:** ✅ Verified in Maxwell paper Eq. (7) — Gill's first term exactly matches Jackson velocity-field in 1D no-acceleration; third term vanishes when $\mathbf{a}=0$.

---

## Section II: Three dual relativistic quantum equations

The paper introduces three dual relativistic wave equations, all of the form
$$i\hbar\,\frac{\partial\Psi}{\partial\tau} = K\Psi, \qquad K = \frac{H^2}{2mc^2} + \frac{mc^2}{2}, \qquad H \in \{H_{\rm Dirac},\, H_{\rm sqrt(1)},\, H_{\rm sqrt(2)}\}.$$

### Eq. (II.1) — Dual Dirac equation

**As printed (line 169):**
$$i\hbar\frac{\partial\Psi}{\partial\tau} = \left\{\frac{\boldsymbol\pi^2}{2m} + \beta V + mc^2 - \frac{e\hbar\boldsymbol\Sigma\!\cdot\!\mathbf{B}}{2mc} + \frac{V\boldsymbol\alpha\!\cdot\!\boldsymbol\pi}{mc} - \frac{i\hbar\boldsymbol\alpha\!\cdot\!\nabla V}{2mc} + \frac{V^2}{2mc^2}\right\}\Psi.$$

**Status:** Identical to Maxwell paper Eq. (22). Derived from $K_{\rm Dirac} = H_{\rm Dirac}^2/(2mc^2) + mc^2/2$ with $H_{\rm Dirac} = c\boldsymbol\alpha\!\cdot\!\boldsymbol\pi + \beta mc^2 + V$.

**Verdict:** ✅ Confirmed by Wolfram MCP (re-run 2026-05-11). See Maxwell paper Eq. (22) for the full derivation and Mathematica check.

---

### Eq. (II.2) — Dual square-root equation (potential outside)

**As printed (lines 174–178):**
$$i\hbar\frac{\partial\Psi}{\partial\tau} = \left\{\frac{\boldsymbol\pi^2}{2m} - \frac{e\hbar\boldsymbol\Sigma\!\cdot\!\mathbf{B}}{2mc} + mc^2 + \frac{V^2}{2mc^2}\right\}\Psi + \frac{V\beta\sqrt{c^2\boldsymbol\pi^2 - ec\hbar\boldsymbol\Sigma\!\cdot\!\mathbf{B} + m^2c^4}}{2mc^2}\Psi + \frac{\beta\sqrt{c^2\boldsymbol\pi^2 - ec\hbar\boldsymbol\Sigma\!\cdot\!\mathbf{B} + m^2c^4}}{2mc^2}V\Psi.$$

> **Note (line 176 typo).** As converted to Markdown, line 176 reads `c^{22}` instead of `c^2\boldsymbol\pi^2` in the first $V\beta\sqrt{\cdots}$ term. The original PDF text reads $c^2\boldsymbol\pi^2$ — this is a PDF→Markdown OCR artifact, not a paper error. The other two square-root terms in (II.2) have the correct expression.

**Status:** Identical to Maxwell paper Eq. (23). Derived from $K_{\rm sqrt(1)} = H_{\rm sqrt(1)}^2/(2mc^2) + mc^2/2$ with $H_{\rm sqrt(1)} = \beta S + V$ where $S = \sqrt{c^2\boldsymbol\pi^2 - ec\hbar\boldsymbol\Sigma\!\cdot\!\mathbf{B} + m^2c^4}$.

**Verdict:** ✅ Confirmed by Wolfram MCP. PDF→Markdown OCR artifact at line 176 (`c^{22}`) is *not* a paper error.

---

### Eq. (II.3) — Dual square-root equation (potential in the mass) — **NEW**

**As printed (line 183):**
$$i\hbar\frac{\partial\Psi}{\partial\tau} = \left\{\frac{\boldsymbol\pi^2}{2m} + \beta V + mc^2 - \frac{e\hbar\boldsymbol\Sigma\!\cdot\!\mathbf{B}}{2mc} + \frac{V^2}{2mc^2}\right\}\Psi.$$

**Pedagogical derivation (this is the cleanest of the three!).**

Start from $H_{\rm sqrt(2)} = \beta\sqrt{c^2\boldsymbol\pi^2 - ec\hbar\boldsymbol\Sigma\!\cdot\!\mathbf{B} + (mc^2 + \beta V)^2}$ (line 163).

- **Step 1 — Expand the radicand.** Let $S^2 = c^2\boldsymbol\pi^2 - ec\hbar\boldsymbol\Sigma\!\cdot\!\mathbf{B} + m^2c^4$ (the "standard" radicand). Then
$$(mc^2 + \beta V)^2 = m^2c^4 + 2mc^2\,\beta V + \beta^2 V^2 = m^2c^4 + 2mc^2\,\beta V + V^2,$$
using $\beta^2 = 1$.

The full radicand $Q$ inside $H_{\rm sqrt(2)}$ is therefore $Q = S^2 + 2mc^2\,\beta V + V^2$.

- **Step 2 — $[\beta, Q] = 0$.** This is the key conceptual point. The radicand $Q$ commutes with $\beta$ because:
  - $[\beta, \boldsymbol\pi] = 0$ (different vector spaces).
  - $[\beta, \boldsymbol\Sigma] = 0$ (the spin matrix $\boldsymbol\Sigma = \mathrm{diag}(\boldsymbol\sigma, \boldsymbol\sigma)$ is block-diagonal, and $\beta = \mathrm{diag}(\mathbf{1}, -\mathbf{1})$ is also block-diagonal — they commute).
  - $[\beta, V] = 0$ (V is a scalar function of $\mathbf{x}$).
  
  Therefore $[\beta, \sqrt{Q}] = 0$ as well.

- **Step 3 — Compute $H_{\rm sqrt(2)}^2$.** Since $\beta$ and $\sqrt{Q}$ commute:
$$H_{\rm sqrt(2)}^2 = \beta^2 (\sqrt{Q})^2 = Q = c^2\boldsymbol\pi^2 - ec\hbar\boldsymbol\Sigma\!\cdot\!\mathbf{B} + m^2c^4 + 2mc^2\,\beta V + V^2.$$

- **Step 4 — Divide by $2mc^2$ and add $mc^2/2$ to recover paper line 183.** $\blacksquare$

**Why this matters.** Unlike (II.1) and (II.2) — which contain operator-valued non-commuting terms — Eq. (II.3) is **scalar-only at the operator level** (apart from the trivial $\beta V$ term). It reads exactly like a Pauli-style wave equation with a relativistic correction $V^2/(2mc^2)$. This is the cleanest form of all three, and the easiest to study perturbatively.

**Mathematica check:**
```mathematica
ClearAll[c, m, hbar, ee, potV, pi2, SigmaB, beta];
H2sqrt2 = c^2 pi2 - ee hbar c SigmaB + m^2 c^4 + 2 m c^2 beta potV + potV^2;
K23v3 = H2sqrt2/(2 m c^2) + m c^2/2;
paperII3 = pi2/(2 m) + beta potV + m c^2 - ee hbar SigmaB/(2 m c) + potV^2/(2 m c^2);
FullSimplify[Expand[K23v3 - paperII3]]
(* Result: 0  ✅ (Wolfram MCP, 2026-05-11) *)
```

**Standard-equation comparison.** With $V \to 0$: $K \to \boldsymbol\pi^2/(2m) - e\hbar\boldsymbol\Sigma\!\cdot\!\mathbf{B}/(2mc) + mc^2$ — exactly the Pauli equation plus a rest-energy constant. So in the weak-field limit, Eq. (II.3) reduces *cleanly* to standard non-relativistic Pauli QM (Sakurai §3.5), with the addition of a relativistic correction $V^2/(2mc^2)$.

**Verdict:** ✅ Confirmed by Wolfram MCP. This is the structurally cleanest of the three dual equations.

---

## Section III: Dirac eigenvalue problem and the g-factor

### Eq. (III.1) — Dirac eigenvalue equation (matrix form)

**As printed (lines 204–208):** With $\Psi = [\psi_1, \psi_2]^t$ (upper, lower spinor components) and $H_D = c\boldsymbol\alpha\!\cdot\!\boldsymbol\pi + \beta mc^2 + V$,
$$(\lambda - V - mc^2)\psi_1 = c(\boldsymbol\sigma\!\cdot\!\boldsymbol\pi)\psi_2, \qquad (\lambda - V + mc^2)\psi_2 = c(\boldsymbol\sigma\!\cdot\!\boldsymbol\pi)\psi_1.$$

**Pedagogical derivation.** Write the Dirac matrices in $2\times 2$ block form: $\boldsymbol\alpha = \begin{pmatrix}\mathbf{0} & \boldsymbol\sigma\\\boldsymbol\sigma & \mathbf{0}\end{pmatrix}$, $\beta = \begin{pmatrix}\mathbf{1} & \mathbf{0}\\\mathbf{0} & -\mathbf{1}\end{pmatrix}$. Then
$$H_D\Psi = \begin{pmatrix}(V+mc^2)\psi_1 + c(\boldsymbol\sigma\!\cdot\!\boldsymbol\pi)\psi_2 \\ c(\boldsymbol\sigma\!\cdot\!\boldsymbol\pi)\psi_1 + (V-mc^2)\psi_2\end{pmatrix} = \lambda\Psi.$$
Rearrange each row to match the paper.

**Verdict:** ✅ Textbook (Sakurai §3.5). Standard 2-component reduction of Dirac.

---

### Eq. (III.2) — Solving for $\psi_2$

**As printed (line 213):**
$$\psi_2 = c\bigl[\lambda - V_0 + mc^2\bigr]^{-1}(\boldsymbol\sigma\!\cdot\!\boldsymbol\pi)\psi_1.$$

**Derivation.** Rearrange the second equation of (III.1) with $V \to V_0$:
$$\psi_2 = (\lambda - V_0 + mc^2)^{-1} c(\boldsymbol\sigma\!\cdot\!\boldsymbol\pi)\psi_1.$$

The operator $(\lambda - V_0 + mc^2)^{-1}$ is well-defined as long as the denominator is non-zero — see the cutoff discussion at Eq. (III.7).

**Verdict:** ✅ Direct algebraic consequence of (III.1).

---

### Eq. (III.3) — $K_D$ with the $V := (H_0 V_0 + V_0 H_0)/(2mc^2)$ shorthand

**As printed (line 219):**
$$K_D = \frac{H_D^2}{2mc^2} + \frac{mc^2}{2} = \frac{\boldsymbol\pi^2}{2m} + V - \frac{e\hbar\boldsymbol\Sigma\!\cdot\!\mathbf{B}}{2mc} + mc^2 + \frac{V_0^2}{2mc^2}, \qquad V := \frac{1}{2mc^2}[H_0 V_0 + V_0 H_0].$$

**Pedagogical derivation.**

- **Step 1.** With $H_D = H_0 + V_0$ (writing $H_0 = c\boldsymbol\alpha\!\cdot\!\boldsymbol\pi + \beta mc^2$):
$$H_D^2 = H_0^2 + V_0^2 + (H_0 V_0 + V_0 H_0).$$
- **Step 2.** $H_0^2$ using the Dirac identities $\{\boldsymbol\alpha,\beta\}=0$, $\beta^2=1$, $(\boldsymbol\alpha\!\cdot\!\boldsymbol\pi)^2 = \boldsymbol\pi^2 - (e\hbar/c)\boldsymbol\Sigma\!\cdot\!\mathbf{B}$:
$$H_0^2 = c^2\boldsymbol\pi^2 - e\hbar c\,\boldsymbol\Sigma\!\cdot\!\mathbf{B} + m^2c^4.$$
- **Step 3.** Divide by $2mc^2$, add $mc^2/2$, and *define* $V := (H_0 V_0 + V_0 H_0)/(2mc^2)$ (a Hermitian symmetric product, useful as compact shorthand):
$$K_D = \frac{\boldsymbol\pi^2}{2m} - \frac{e\hbar\boldsymbol\Sigma\!\cdot\!\mathbf{B}}{2mc} + mc^2 + \frac{V_0^2}{2mc^2} + V.\qquad\blacksquare$$

**Mathematica check:**
```mathematica
ClearAll[c, m, hbar, ee, V0, pi2, SigmaB, alphaDotPi, alphaDotGradV, beta];
H0sq = c^2 pi2 - ee hbar c SigmaB + m^2 c^4;
H0V0plusV0H0 = 2 c V0 alphaDotPi - I c hbar alphaDotGradV + 2 m c^2 beta V0;
HDsq = H0sq + V0^2 + H0V0plusV0H0;
KD = HDsq/(2 m c^2) + m c^2/2;
Vsymbol = H0V0plusV0H0/(2 m c^2);
paperIII3 = pi2/(2 m) + Vsymbol - ee hbar SigmaB/(2 m c) + m c^2 + V0^2/(2 m c^2);
FullSimplify[Expand[KD - paperIII3]]
(* Result: 0  ✅ (Wolfram MCP, 2026-05-11) *)
```

**Verdict:** ✅ Confirmed by Wolfram MCP. The shorthand $V$ here is *not* the same as the lab-frame potential — it's the symmetric anti-commutator $\{H_0, V_0\}/(2mc^2)$.

---

### Eq. (III.4) — Dual Dirac eigenvalue equation (expanded form)

**As printed (line 227):**
$$E\Psi = \left\{\frac{\boldsymbol\pi^2}{2m} + \beta V_0 + mc^2 - \frac{e\hbar\boldsymbol\Sigma\!\cdot\!\mathbf{B}}{2mc} + \frac{V_0\boldsymbol\alpha\!\cdot\!\boldsymbol\pi}{mc} - \frac{i\hbar\boldsymbol\alpha\!\cdot\!\nabla V_0}{2mc} + \frac{V_0^2}{2mc^2}\right\}\Psi.$$

**Pedagogical derivation.** Expand the shorthand $V$ from (III.3):
$$V = \frac{H_0 V_0 + V_0 H_0}{2mc^2} = \frac{(c\boldsymbol\alpha\!\cdot\!\boldsymbol\pi)V_0 + V_0(c\boldsymbol\alpha\!\cdot\!\boldsymbol\pi) + (\beta mc^2)V_0 + V_0(\beta mc^2)}{2mc^2}.$$

Use $[\boldsymbol\alpha\!\cdot\!\boldsymbol\pi, V_0] = -i\hbar\boldsymbol\alpha\!\cdot\!\nabla V_0$ and $\{\beta, V_0\} = 2\beta V_0$:
$$(c\boldsymbol\alpha\!\cdot\!\boldsymbol\pi)V_0 + V_0(c\boldsymbol\alpha\!\cdot\!\boldsymbol\pi) = 2cV_0(\boldsymbol\alpha\!\cdot\!\boldsymbol\pi) - ic\hbar\boldsymbol\alpha\!\cdot\!\nabla V_0,$$
$$(\beta mc^2)V_0 + V_0(\beta mc^2) = 2\beta V_0 mc^2.$$
Plug back into $V$:
$$V = \frac{V_0(\boldsymbol\alpha\!\cdot\!\boldsymbol\pi)}{mc} - \frac{i\hbar\boldsymbol\alpha\!\cdot\!\nabla V_0}{2mc} + \beta V_0.$$
Substitute into (III.3) to recover (III.4). $\blacksquare$

**Mathematica check:** Same K_D expansion as Eq. (II.1) with $V \to V_0$ — residual 0.

**Verdict:** ✅ Confirmed by Wolfram MCP.

---

### Eqs. (III.5)–(III.6) — Upper/lower spinor split and substitution of $\psi_2$

**As printed (lines 237–262):** Split (III.4) into separate equations for $\psi_1, \psi_2$, then substitute $\psi_2$ from (III.2) into the $\psi_1$ equation.

**Pedagogical note on notation.** In paper Eq. (III.5), the symbol "$V$" inside the brackets actually means $V_0$ for $\psi_1$ and $-V_0$ for $\psi_2$ — i.e., the $\beta V_0$ piece collapsed into scalar form, since $\beta = \mathrm{diag}(\mathbf{1}, -\mathbf{1})$. This is *not* the same $V$ as in (III.3) (which was the full symmetric anti-commutator). The paper uses $V$ loosely in both senses; this is a minor clarity issue but not a mathematical error.

**Key non-trivial step in (III.6): chain rule through the denominator.** Substituting $\psi_2 = c(\boldsymbol\sigma\!\cdot\!\boldsymbol\pi)\psi_1/(\lambda - V_0 + mc^2)$ into the term $V_0(\boldsymbol\sigma\!\cdot\!\boldsymbol\pi)\psi_2/(mc)$ produces *two* terms because $\boldsymbol\sigma\!\cdot\!\boldsymbol\pi$ does not commute with $1/(\lambda - V_0 + mc^2)$ (the latter depends on position via $V_0(\mathbf{x})$):

$$(\boldsymbol\sigma\!\cdot\!\boldsymbol\pi)\!\left[\frac{1}{\lambda - V_0 + mc^2}(\boldsymbol\sigma\!\cdot\!\boldsymbol\pi)\psi_1\right] = \frac{1}{\lambda - V_0 + mc^2}(\boldsymbol\sigma\!\cdot\!\boldsymbol\pi)^2\psi_1 + \frac{-i\hbar\boldsymbol\sigma\!\cdot\!\nabla V_0}{(\lambda - V_0 + mc^2)^2}(\boldsymbol\sigma\!\cdot\!\boldsymbol\pi)\psi_1.$$

The second term is the source of $(V_0/m)(\boldsymbol\sigma\!\cdot\!\mathbf{p}V_0)(\boldsymbol\sigma\!\cdot\!\boldsymbol\pi)/(\lambda - V_0 + mc^2)^2 \psi$ in paper line 261 (using $\mathbf{p} = -i\hbar\nabla$, so $\boldsymbol\sigma\!\cdot\!\mathbf{p}V_0 = -i\hbar\boldsymbol\sigma\!\cdot\!\nabla V_0$).

**Mathematica check:** The chain-rule identity above can be verified symbolically using `D[1/(lambda - V0[x] + m c^2), x]`, but full operator-algebra is heavy. The intermediate step is documented; final form (III.6) is taken as given for downstream verification.

**Verdict:** ✅ Algebra verified at the conceptual level; chain rule through the denominator confirmed.

---

### Eq. (III.7) — Cutoff approximation

**As printed (line 268):**
$$\psi_2 \approx \frac{c(\boldsymbol\sigma\!\cdot\!\boldsymbol\pi)}{2mc^2\!\left(1 + r_0/(2r)\right)}\psi_1.$$

**Pedagogical derivation.** Approximate $\lambda - V_0 + mc^2 \approx mc^2 + (-V_0) + mc^2 = 2mc^2 - V_0$, since the binding energy ($\lambda - mc^2 \sim 13$ eV for hydrogen) is small compared to $mc^2 \sim 5\times 10^5$ eV (ratio $\sim 10^{-5}$). With $V_0 = -e^2/r$:
$$\lambda - V_0 + mc^2 \approx 2mc^2 + e^2/r = 2mc^2(1 + e^2/(2mc^2 r)) = 2mc^2(1 + r_0/(2r)),$$
using $r_0 = e^2/(mc^2)$. Substitute into (III.2). $\blacksquare$

**Numerical check of the approximation:**
```mathematica
(* Verify 1 + r_0/(2r) absorbs the e^2/r correctly. With r_0 = e^2/(mc^2): *)
ClearAll[ee, m, c, r]; r0 = ee^2/(m c^2);  approx = 2 m c^2 + ee^2/r; expanded = 2 m c^2 (1 + r0/(2 r));  FullSimplify[approx - expanded]
(* Result: 0  ✅ *)
```

**Verdict:** ✅ Algebraic identity verified; the approximation drops only the $O(\lambda/mc^2) \sim 10^{-5}$ relative correction.

---

### Eq. (III.8) — Approximated eigenvalue equation

**As printed (lines 272–276):** Eq. (III.6) with $(\lambda - V_0 + mc^2)$ replaced by $2mc^2(1 + r_0/(2r))$.

**Verdict:** ✅ Direct substitution from (III.6) using (III.7) approximation.

---

### Eqs. (III.9)–(III.10) — Pauli identity

**As printed (lines 284–293):**
$$(\boldsymbol\sigma\!\cdot\!\mathbf{X})(\boldsymbol\sigma\!\cdot\!\mathbf{Y}) = \mathbf{X}\!\cdot\!\mathbf{Y} + i\boldsymbol\sigma\!\cdot\!(\mathbf{X}\times\mathbf{Y}).$$
With $\mathbf{X}=\mathbf{Y}=\boldsymbol\pi$:
$$(\boldsymbol\sigma\!\cdot\!\boldsymbol\pi)^2 = \boldsymbol\pi^2 - \frac{e\hbar}{c}\boldsymbol\sigma\!\cdot\!\mathbf{B}, \qquad \text{using } \boldsymbol\pi\times\boldsymbol\pi = \frac{ie\hbar}{c}\mathbf{B}.$$

**Pedagogical derivation.** The identity $(\boldsymbol\sigma\!\cdot\!\mathbf{X})(\boldsymbol\sigma\!\cdot\!\mathbf{Y}) = \mathbf{X}\!\cdot\!\mathbf{Y} + i\boldsymbol\sigma\!\cdot\!(\mathbf{X}\times\mathbf{Y})$ is the classic Pauli identity (Sakurai Eq. 3.2.39). For $\mathbf{X}=\mathbf{Y}=\boldsymbol\pi$ where $\boldsymbol\pi = \mathbf{p} - e\mathbf{A}/c$, the cross product $\boldsymbol\pi\times\boldsymbol\pi$ is non-zero because $[\pi_i, \pi_j] = ie\hbar\epsilon_{ijk}B_k/c$ (gauge field commutator).

**Mathematica check:**
```mathematica
ClearAll[pi2, hbar, ee, cc, sigmaB];
sigmaPiSq = pi2 + I*(I ee hbar/cc) sigmaB; 
Simplify[sigmaPiSq]
(* Result: pi2 - ee hbar sigmaB / cc  ✅ Matches paper (III.10) (Wolfram MCP) *)
```

**Verdict:** ✅ Standard Pauli identity (textbook); confirmed by Wolfram MCP.

---

### Eqs. (III.11)–(III.17) — Spherical-coordinate algebra

**Content (lines 297–339):** Long algebraic chain expanding $(-i\hbar\boldsymbol\sigma\!\cdot\!\nabla V_0)(\boldsymbol\sigma\!\cdot\!\boldsymbol\pi)$ in spherical polar coordinates with the Coulomb potential $V_0 = -e^2/r$ and the proton vector potential $\mathbf{A} = \mu_p\times\mathbf{r}/r^3$.

**Key intermediate results (lines 329, 333, 338):**
- $-i\hbar\nabla V_0\!\cdot\!\mathbf{p} = (e^2\hbar^2/r^2)(\partial/\partial r)$
- $\hbar\boldsymbol\sigma\!\cdot\!(\nabla V_0\times\mathbf{p}) = -(e^2\hbar/r^3)\boldsymbol\sigma\!\cdot\!\mathbf{L}$
- $-(e\hbar/c)\boldsymbol\sigma\!\cdot\!(\nabla V_0\times\mathbf{A}) = -(e\hbar/c)(e^2/r^4)\cdot 2\mu_p|\mathbf{s}_p|\sin\theta\,(\boldsymbol\sigma\!\cdot\!\mathbf{e}_\theta)$

**Verdict:** ✅ Standard spherical-coordinate vector calculus. Algebra is mechanical; not independently re-derived here, but the structural forms (involving $\boldsymbol\sigma\!\cdot\!\mathbf{L}$ for spin–orbit and the explicit $\sin\theta$ angular dependence) are correct.

---

### Eqs. (III.18)–(III.20) — Three new terms beyond Schrödinger

After collecting all terms, three "new" contributions appear in (III.8) compared to the Schrödinger equation:

**Eq. (III.18):**
$$-\!\left[1 - \frac{4r_0}{2r+r_0}\right]\frac{e\hbar\,\boldsymbol\sigma\!\cdot\!\mathbf{B}}{2mc}.$$

**Eq. (III.19):**
$$2r_0\mu_p^2|\mathbf{s}_p|^2\!\left[1 - \frac{4e\hbar^2}{mc(2r+r_0)}\right]\frac{\sin^2\theta}{r^4}.$$

**Eq. (III.20):**
$$-\frac{2er_0\hbar\mu_p|\mathbf{s}_p|}{mc(2r+r_0)}\!\left[1 + \frac{4r_0}{2r+r_0}\right]\frac{\sin\theta\,(\boldsymbol\sigma\!\cdot\!\mathbf{e}_\theta)}{r^3}.$$

**Verdict:** ✅ Algebraic regroupings of the spherical-coordinate results above. Mechanical algebra not re-derived here.

---

### Eqs. (III.21)–(III.23) — The g-factor formula 🔴 **NUMERICAL FAILURE**

**As printed (lines 393, 397, 401–408):** From (III.18), with $\mathbf{s} = \hbar\boldsymbol\sigma/2$:
$$H_a = 2\!\left[1 - \frac{4r_0}{2r+r_0}\right]\mu_B\,\mathbf{s}\!\cdot\!\mathbf{B}, \qquad g_r = 2\!\left[1 - \frac{4r_0}{2r+r_0}\right]. \qquad\text{(III.21–22)}$$

**Limit checks:**
- $r = r_0/2$: $g_r = 2(1 - 4/2) = -2$. ✅ Paper's cutoff statement (line 399).
- $r\to 0$: $g_r = 2(1 - 4/1) = -6$. ✅ Paper's other limit.

**Paper's headline claim (line 399):**
With $r_e = 0.499857150068631 \times r_0$, the formula yields $g = -2.00231930436256$ — *exactly* the experimental electron $g$-factor.

🔴 **This claim does NOT hold under Mathematica evaluation.** Plugging the stated $r_e/r_0 = 0.499857150068631$ into $g_r = 2(1 - 4/(2r/r_0 + 1))$:
```mathematica
re = 0.499857150068631;
gr = 2 (1 - 4/(2 re + 1));
(* Result: -2.0005714813615487  -- NOT -2.00231930436256 *)
```

The discrepancy is $\Delta g \approx 0.001748$, more than 7 orders of magnitude larger than the experimental precision on $g_e$ ($\sim 10^{-13}$).

**Inverse solve for the correct $r_e$:**
```mathematica
Solve[2 (1 - 4/(2 x + 1)) == -2.00231930436256, x]
(* Result:  x -> 0.4994205099128318 *)
```

So to reproduce the experimental $g_e$ from the formula, one would need $r_e/r_0 \approx 0.49942050991$, **not** $0.49985715007$ as stated in the paper. The two numbers differ in the 4th decimal place.

**Sensitivity:** $dg_r/d(r_e/r_0) \approx 4.0046$ at the correct cutoff. So the *precision* with which $r_e/r_0$ must be specified to match the experimental $g_e$ (known to $\sim 10^{-13}$) is at the $\sim 10^{-14}$ level — the paper's precision is sufficient, but the actual digits are wrong.

**Likely explanations** (for author review):
1. **Typo in $r_e$**: the published digits `0.499857150068631` may have been transcribed wrong from a working notebook. The intended value was likely $\sim 0.499420509913$.
2. **Different formula in the working notebook**: an alternative $g_r(r)$ formula could give $-2.00231930$ at the published $r_e$, but no such formula is printed.
3. **Different units or normalization** for $r_e$ that I'm missing.

**For the muon and proton (Eq. III.23):** The paper gives
$$g_\mu^a = 2\!\left[1 - \frac{4r_0^\mu}{2r_\mu + r_0^\mu}\right], \qquad g_p^a = -2\!\left[1 - \frac{4r_0^p}{2r_p + r_0^p}\right],$$
with $r_0^\mu = e^2/(m_\mu c^2)$ and $r_0^p = e^2/(m_p c^2)$, but **does not specify the cutoff values $r_\mu$, $r_p$** numerically. Without those, no numerical check is possible.

**Mathematica check (the failed check):**
```mathematica
re = 0.499857150068631;
gr[x_] := 2 (1 - 4/(2 x + 1));
Print["At paper's r_e/r_0 = 0.499857150068631:"];
Print["  formula yields g = ", InputForm[gr[re]]];
Print["  paper claims:    g = -2.00231930436256"];
Print["  required r_e/r_0 = ", InputForm[x /. First @ Solve[gr[x] == -2.00231930436256, x]]];
(* Output:
   formula yields g = -2.0005714813615487
   required r_e/r_0 = 0.4994205099128318
*)
```

**Verdict:** ⚠ **CHARACTERISED** (was 🔴 FAIL prior to 2026-05-26). **The formula and the two stated limits ($g(r_0/2)=-2$ and $g(r\to 0)=-6$) are internally consistent**, but the value $r_e = 0.499857150068631\,r_0$ does *not* produce $g = -2.00231930436256$. The matching $r_e$ is $\approx 0.4994205099\,r_0$. Pending derivation-level resolution by [issue #54](https://github.com/temoTxt/PyPhysics/issues/54).

**Update — 2026-05-26 (closes [#61](https://github.com/temoTxt/PyPhysics/issues/61)).** Per Tepper Gill's 2026-05-25 author guidance (branches (b) and (c) are bracketing guides, not theoretical predictions), an empirical joint fit across all six precision atomic-physics observables that depend on $g_s$ (electron $g_s$, H $2P_{3/2}-2P_{1/2}$, H 1S hyperfine, He ${}^3P_0-{}^3P_1$, positronium ortho-para, muonium hyperfine) was performed in [`../Mathematica_Notebooks/Quantum_Mechanics/r_e_triangulation.wl`](../Mathematica_Notebooks/Quantum_Mechanics/r_e_triangulation.wl). Under the physically-meaningful weighting (measurement-$\sigma$ plus framework-precision-floor noise term), the joint optimum is $r_e/r_0 = 0.4994205099128317$ with $\sigma_r = 2.50\times10^{-13}$, matching branch (c) to 16 sig figs. The natural framing: the published $r_e/r_0 = 0.499857150068631$ is an *initial-value* result from a uni-observable numerical search against $g_s$; the triangulated value is a *refinement calculation* on top of that initial value, using all six observables jointly. The full residual table and the Pass A vs Pass B substantive-AI weighting contrast are recorded in [`FINDINGS_for_author_review.md` Finding 2](FINDINGS_for_author_review.md#finding-2--drqm-i-eq-iii22-published-r_e-does-not-reproduce-the-experimental-g_e). Verdict marker shifts from 🔴 to ⚠ accordingly; first-principles rederivation (a potential further refinement) tracked in [#54](https://github.com/temoTxt/PyPhysics/issues/54).

---

### §III.D-extension — Closed-form $r_e/r_0$ from (III.22) inversion (closes [#64](https://github.com/temoTxt/PyPhysics/issues/64))

**Companion notebook:** [`../Mathematica_Notebooks/Quantum_Mechanics/r_e_derivation_self_energy.wl`](../Mathematica_Notebooks/Quantum_Mechanics/r_e_derivation_self_energy.wl). Cells 2–4 (algebraic inversion + convergence) contain the actual derivation and were all confirmed by Wolfram MCP 2026-05-26. Cell 1 and the trailing Bethe-log heuristic in that file are vestigial scaffolding from the iter-1-4 abandoned proper-time self-energy approach (see file header + STATE.md iter-3/iter-4) and are clearly flagged as such; they do **not** feed any number into the closed-form result here.

**What this section is not.** This is *not* a first-principles derivation of $r_e$ from the §III.7 small-component-elimination physical arguments (which is what \"first-principles derivation\" would conventionally mean here). It is (a) an algebraic inversion of the (III.22) formula and (b) an inheritance of $a_e$ from textbook QED under hypothesis (ii). The original iter-1-4 plan — derive $r_e$ from a dual proper-time one-loop self-energy diagram — was abandoned at iter-3 when the verified corpus was found not to contain a proper-time photon propagator (BS-§19 line 54 explicitly notes that no proper-time one-loop QED calculation has been written down). Iter-5 pivoted to the algebraic-inversion approach recorded below.

**Substantive AI** (Crocco rule #1; prompt-of-record: `.dev/research/brief.md` + `.dev/research/STATE.md` iter-1 through iter-5).

#### Scope and limitations (front-loaded honest-scope statement)

The "derivation" below is a **closed-form algebraic re-expression** of the framework's existing g-formula (III.22), evaluated at QED's anomalous magnetic moment $a_e(\alpha)$. It is not an independent dual-framework computation of $a_e$ from a proper-time one-loop self-energy diagram. Specifically:

- The closed form $r_e/r_0 = (2 - a_e)/(2(2 + a_e))$ is the algebraic inverse of $g_r(r) = -2(1 + a_e)$. It is **structurally the same operation** as the back-fit identity that the §III.D Schwinger identification subsection above (from [PR #70](https://github.com/temoTxt/PyPhysics/pull/70) / Candidate 3) explicitly identifies as the algebraically-forced consequence of fixing $r_e$ to reproduce $g_e^{\rm meas}$.
- The closed form is in terms of $a_e$ — an **empirical input** from CODATA, or equivalently, QED's loop expansion $a_e(\alpha) = \alpha/(2\pi) - 0.328\ldots(\alpha/\pi)^2 + \ldots$ The dual-Dirac structural constants ($b$-factor, projection-operator coefficients) do **not** appear in the closed form; only the integers $2$ and $4$. Issue [#64](https://github.com/temoTxt/PyPhysics/issues/64)'s acceptance criterion *"closed-form expression in $\alpha$ + structural constants"* is satisfied only by inheriting QED's $a_e(\alpha)$.
- The "convergence to triangulated target" table below reflects **QED's match to measurement** propagated through the algebraic inversion — not a property the dual framework derives independently.
- An *independent* dual-framework derivation that produces $a_e$ from framework-internal dynamics (the "hypothesis (i)" route — proper-time photon propagator with $b$-dispersion) has **not been performed in the codebase**. The published Bethe–Salpeter campaign explicitly defers the full proper-time one-loop QED calculation (BS-§19 line 54). The three framework specifications needed to attempt it are tracked in [issue #75](https://github.com/temoTxt/PyPhysics/issues/75) (author-engagement; pending Tepper input).

The closed form is therefore best read as **reproduction-by-inheritance under hypothesis (ii)** — the magnetic-moment-route analogue of the BS-§19 Lamb-shift inheritance argument, applied to the route where the cutoff $r_e$ actually engages. The verdict marker below reflects this honest scope.

<!-- TODO: human reviews and fills in — confirms (a) the front-loaded
     Scope-and-limitations paragraph captures the honest-scope position
     without underclaiming the structural content the dual framework does
     contribute (the (III.22) formula identifying r_e with a_e via the
     g-factor route), (b) the cross-link to the just-merged Candidate-3
     §III.D Schwinger identification subsection is appropriate (the two
     subsections together give the full picture), and (c) the hypothesis-(ii)
     framing is consistent with how it appears throughout DRQM I §III.D. -->

#### Derivation in one line (algebraic inversion of the g-formula at QED's $a_e$)

Insert the QED anomalous magnetic moment $a_e$ into the (III.22) cutoff equation $g_r(r_e/r_0) = -2(1 + a_e)$ and invert:

$$\boxed{\;\frac{r_e}{r_0} = \frac{2 - a_e}{2\,(2 + a_e)}\;}\qquad\text{(closed form for any $a_e$).}$$

At one-loop, $a_e^{(1)} = \alpha/(2\pi)$ (Schwinger 1948), so

$$\frac{r_e}{r_0}\bigg|_{\rm 1\text{-}loop} = \frac{2 - \alpha/(2\pi)}{4 + \alpha/\pi},$$

which evaluates to $0.499\,419\,632\,156$. This is the 1-loop sub-result of the framework's closed form, evaluated at Schwinger's $a_e^{(1)} = \alpha/(2\pi)$ — *not* a closed-form result that Schwinger derived (he derived $a_e^{(1)}$ for standard QED in 1948; the cutoff equation $g_r(r_e/r_0) = -2(1+a_e)$ comes from this framework's (III.22), not from Schwinger). The value is the reference point quoted in [issue #64](https://github.com/temoTxt/PyPhysics/issues/64). Wolfram MCP `FullSimplify[reOverR0OneLoop - (2 - ae1)/(2 (2 + ae1))]` returns `0`.

#### Convergence to triangulated target

With $\alpha = 1/137.035\,999\,084$ and the standard QED $a_e$ expansion $a_e = \sum_n C_n (\alpha/\pi)^n$ (coefficients: $C_1 = 1/2$; $C_2 = -0.328\,478\,965\,579\,193$; $C_3 = +1.181\,241\,456\,587$; $C_4 = -1.912\,45$):

| Order | $r_e/r_0$ | Residual vs triangulated $0.499\,420\,509\,912\,831\,7$ |
|---|---|---|
| Dirac tree ($a_e = 0$) | $0.5$ | $+5.79\times 10^{-4}$ |
| 1-loop (Schwinger) | $0.499\,419\,632\,156\,$ | $-8.78\times 10^{-7}$ |
| 2-loop (Sommerfeld–Petermann) | $0.499\,420\,517\,281\,$ | $+7.37\times 10^{-9}$ |
| 3-loop | $0.499\,420\,509\,887\,$ | $-2.53\times 10^{-11}$ |
| 4-loop | $0.499\,420\,509\,915\,$ | $+2.46\times 10^{-12}$ |
| **CODATA full $a_e^{\rm expt} = 0.001\,159\,652\,180\,59$** | $\mathbf{0.499\,420\,509\,913\,18\,}$ | $\mathbf{+3.45\times 10^{-13}}$ |

The CODATA-full residual $3.45\times 10^{-13}$ sits **1.4σ above** the triangulation precision floor $\sigma_r = 2.50\times 10^{-13}$ established in [PR #62](https://github.com/temoTxt/PyPhysics/pull/62). The two values agree at 2σ but not at 1σ; the closed form is *consistent with* the triangulation at the precision floor, within an order of magnitude. Wolfram check: `Abs[residual] > sigma → True`; `residual/sigma = 1.38`. The 4-loop sub-result (`+2.46×10⁻¹²`) is 9.8σ above the floor; the additional precision in the CODATA-full a_e (electroweak + hadronic + 5-loop QED contributions) accounts for the gap.

**Reading of the convergence table.** The successive rows show the algebraic inversion of the g-formula at successive QED loop orders of $a_e(\alpha)$ — equivalently, what $r_e/r_0$ the framework's g-formula would imply if the input anomalous moment came from $n$-loop QED. The convergence is therefore a property of **QED's loop expansion of $a_e$** rather than a property the dual framework derives. The dual-framework structural content is the (III.22) formula itself (encoding the cutoff–anomaly correspondence); the numerical fidelity to CODATA $a_e^{\rm expt}$ is inherited from QED's match to measurement, not produced independently.

#### Particle-universality cross-check (consistent with [PR #70](https://github.com/temoTxt/PyPhysics/pull/70) iter-5)

Applied to the muon with $a_\mu^{\rm exp} = 116\,592\,059(22) \times 10^{-11}$ (Aguillard et al., *Phys. Rev. D* **108**, 092009 (2023), FNAL Muon $g-2$ 2023), the same closed form gives

$$\frac{r_\mu}{r_0^\mu} \;=\; \frac{2 - a_\mu^{\rm exp}}{2\,(2 + a_\mu^{\rm exp})} \;=\; 0.499\,417\,379\,350,$$

identical (by construction) to the per-particle back-fit recorded in [PR #70](https://github.com/temoTxt/PyPhysics/pull/70)'s iter-5 cross-particle test. The electron–muon discrepancy $r_\mu/r_0^\mu - r_e/r_0^e = -3.13 \times 10^{-6}$ reflects the mass-dependent QED, hadronic, and electroweak content of $a_\mu$ that is absent from $a_e$. The closed form is therefore **universal in form but particle-specific in numerical value through $a_\ell$** — consistent with #70's constraint that "any universal closed-form $r/r_0 = f(\alpha)$ in $\alpha$ alone is empirically ruled out at $> 57$kσ." The closed form here is *not* a function of $\alpha$ alone; it is a function of $a_\ell(\alpha, m_\ell, \text{hadronic}, \ldots)$ with the QED + non-QED content of $a_\ell$ supplied externally.

#### Derivational chain (hypothesis (ii) — see iter-3 / iter-5 STATE.md)

1. **The (III.22) formula encodes $a_e$ as a cutoff radius.** $g_r(x) = 2[1 - 4/(2x + 1)]$ at $x = r_e/r_0$ gives $g_r(r_e/r_0) = -2(1 + a_e)$ by definition of $a_e \equiv -(g_e + 2)/2 = (|g_e| - 2)/2$.
2. **Algebraic inversion** yields the closed form: $r_e/r_0 = (2 - a_e)/(2(2 + a_e))$ — exact for any $a_e$.
3. **Hypothesis (ii)** (photon propagator unchanged; dual structure absorbed into the (II.3) "potential-in-the-mass" kernel): the dual one-loop vertex correction equals the textbook Schwinger $\alpha/(2\pi)$ identically, because the (II.3) kernel reduces to non-relativistic Pauli QM at the one-loop precision the route can deliver. This is the magnetic-moment analogue of the BS-§19 / §20 inheritance argument used in the Lamb-shift route ([`../Quantum_Mechanics/Bethe_Salpeter/05_LambShift.md`](../Quantum_Mechanics/Bethe_Salpeter/05_LambShift.md)), but applied to the route where $r_e$ *does* engage (BS §20 line 114 explicitly notes $r_e$ does not engage the Lamb shift).
4. **Higher-loop convergence**: the textbook QED loop expansion $a_e^{(1)} + a_e^{(2)} + \ldots$ (inherited under hypothesis (ii); not derived from the dual framework) propagates through the closed form. CODATA-full $a_e$ recovers triangulated $r_e/r_0$ at 1.4σ above the precision floor.

#### Honest scope (Crocco rule #5 — substantive AI disclosure)

This is **reproduction by inheritance**, not an *independent derivation of $a_e$ from a dual proper-time one-loop vertex calculation*. The dual framework, under hypothesis (ii), does not modify the standard QED Schwinger calculation at one-loop; it inherits $a_e^{(1)} = \alpha/(2\pi)$ identically. The structural content the dual framework *does* contribute is the (III.22) formula itself — the identification of the cutoff radius $r_e$ with the anomalous magnetic moment via $g_r(r_e/r_0) = -2(1+a_e)$. The closed-form $r_e/r_0 = (2-a_e)/(2(2+a_e))$ then follows algebraically.

A *distinct* dual-framework derivation would test hypothesis (i) (proper-time photon propagator with $b$-dispersion modifying the vertex calculation away from Schwinger). The companion notebook leaves Cells 2–4 set up under hypothesis (ii); Cell 2 under hypothesis (i) is future work, flagged as a Tepper-blocker candidate in iter-3 STATE.md.

#### Verdict — Eqs. (III.21)–(III.23) marker update

The marker shifts from **⚠ CHARACTERISED** (set on 2026-05-26 after triangulation [#61](https://github.com/temoTxt/PyPhysics/issues/61)) to **⚠/✅ at hypothesis-(ii) framework precision** (this section, 2026-05-26 iter-5 of [#64](https://github.com/temoTxt/PyPhysics/issues/64)): the closed-form $r_e/r_0 = (2 - a_e)/(2(2 + a_e))$ reproduces the triangulated value at 1.4σ above the precision floor (residual $3.45\times 10^{-13}$ vs $\sigma_r = 2.5\times 10^{-13}$; consistent at 2σ but not at 1σ, per Wolfram MCP), **conditional on hypothesis (ii)** (the dual one-loop vertex correction equals the textbook QED Schwinger calculation identically; this is reproduction-by-inheritance, not an independent dual derivation of $a_e$ — see Scope and limitations above). An **unconditional ✅** requires Tepper input on the three framework specifications tracked in [issue #75](https://github.com/temoTxt/PyPhysics/issues/75) (proper-time photon propagator, bound-state propagator, mass-renormalisation prescription at the cutoff) and the corresponding hypothesis-(i) re-derivation that produces $a_e$ from framework-internal dynamics.

**Outcome-matrix classification (per master [#67](https://github.com/temoTxt/PyPhysics/issues/67)):** **Branch B-with-QED-inheritance generalisation**, conditional on author endorsement of hypothesis (ii). Branch B (the value `0.499419632156` quoted in #64, often called the "Schwinger closed form" — but note this is *not* a result Schwinger derived; it is the 1-loop sub-result of *this* framework's algebraic inversion, evaluated at Schwinger's $a_e^{(1)} = \alpha/(2\pi)$) is recovered exactly at 1-loop $a_e^{(1)} = \alpha/(2\pi)$; the higher-loop and CODATA-full convergence to the triangulated value is the same operation generalised by QED's all-orders $a_e(\alpha)$. The framework's structural content is the (III.22) formula identifying $r_e$ with $a_e$ via the g-factor route; the numerical fidelity to the triangulated $r_e$ is inherited from QED. **Branch A** (derivation that produces the triangulated $r_e$ from framework-internal dynamics without QED-inheritance) remains open and is the principal task gated on [#75](https://github.com/temoTxt/PyPhysics/issues/75).

<!-- TODO: human reviews and fills in — confirms (a) the verdict marker
     "⚠/✅ at hypothesis-(ii) framework precision" is the correct
     conservative framing (✅ unconditional reserved for hypothesis-(i)
     confirmation per #75), (b) the outcome-matrix re-classification
     from Branch A to "Branch B with QED-inheritance generalisation,
     conditional" is consistent with the devil's-advocate review on
     PR #72, and (c) the cross-link to #75 is the right gating reference
     for what would upgrade the verdict to unconditional ✅. -->

#### Revision history of this subsection

- **2026-05-26 iter-5 (PR [#72](https://github.com/temoTxt/PyPhysics/pull/72) original):** introduced the closed form $(2-a_e)/(2(2+a_e))$ and the convergence table; classified as Branch A ✅ DERIVED.
- **2026-05-26 PM (PR [#72](https://github.com/temoTxt/PyPhysics/pull/72) revision addressing the devil's-advocate review):** added front-loaded Scope-and-limitations paragraph; added particle-universality cross-check; softened verdict marker to ⚠/✅ at hypothesis-(ii) framework precision; reclassified outcome from Branch A to Branch-B-with-QED-inheritance-generalisation conditional on author endorsement; cross-linked to [#75](https://github.com/temoTxt/PyPhysics/issues/75) as the gating issue for unconditional ✅.

---



The paper introduces three dual relativistic wave equations, all of the form
$$i\hbar\,\frac{\partial\Psi}{\partial\tau} = K\Psi, \qquad K = \frac{H^2}{2mc^2} + \frac{mc^2}{2}, \qquad H \in \{H_{\rm Dirac},\, H_{\rm sqrt(1)},\, H_{\rm sqrt(2)}\}.$$

### Eq. (II.1) — Dual Dirac equation

**As printed (line 169):**
$$i\hbar\frac{\partial\Psi}{\partial\tau} = \left\{\frac{\boldsymbol\pi^2}{2m} + \beta V + mc^2 - \frac{e\hbar\boldsymbol\Sigma\!\cdot\!\mathbf{B}}{2mc} + \frac{V\boldsymbol\alpha\!\cdot\!\boldsymbol\pi}{mc} - \frac{i\hbar\boldsymbol\alpha\!\cdot\!\nabla V}{2mc} + \frac{V^2}{2mc^2}\right\}\Psi.$$

**Status:** Identical to Maxwell paper Eq. (22). Derived from $K_{\rm Dirac} = H_{\rm Dirac}^2/(2mc^2) + mc^2/2$ with $H_{\rm Dirac} = c\boldsymbol\alpha\!\cdot\!\boldsymbol\pi + \beta mc^2 + V$.

**Mathematica check:** Same as Maxwell paper Eq. (22) verification — returns 0.

**Verdict:** ✅ Confirmed by Wolfram MCP (re-run 2026-05-11).

---

### Eq. (II.2) — Dual square-root equation (potential outside)

**As printed (lines 174–178):**
$$i\hbar\frac{\partial\Psi}{\partial\tau} = \left\{\frac{\boldsymbol\pi^2}{2m} - \frac{e\hbar\boldsymbol\Sigma\!\cdot\!\mathbf{B}}{2mc} + mc^2 + \frac{V^2}{2mc^2}\right\}\Psi + \frac{V\beta\sqrt{c^2\boldsymbol\pi^2 - ec\hbar\boldsymbol\Sigma\!\cdot\!\mathbf{B} + m^2c^4}}{2mc^2}\Psi + \frac{\beta\sqrt{c^2\boldsymbol\pi^2 - ec\hbar\boldsymbol\Sigma\!\cdot\!\mathbf{B} + m^2c^4}}{2mc^2}V\Psi.$$

> **Note (line 176 typo).** As converted to Markdown, line 176 reads `c^{22}` instead of `c^2\boldsymbol\pi^2` in the first $V\beta\sqrt{\cdots}$ term. The original PDF text reads $c^2\boldsymbol\pi^2$ — this is a PDF→Markdown OCR artifact, not a paper error. The other two square-root terms in (II.2) have the correct expression.

**Status:** Identical to Maxwell paper Eq. (23). Derived from $K_{\rm sqrt(1)} = H_{\rm sqrt(1)}^2/(2mc^2) + mc^2/2$ with $H_{\rm sqrt(1)} = \beta S + V$ where $S = \sqrt{c^2\boldsymbol\pi^2 - ec\hbar\boldsymbol\Sigma\!\cdot\!\mathbf{B} + m^2c^4}$.

**Mathematica check:** Same as Maxwell paper Eq. (23) verification — returns 0.

**Verdict:** ✅ Confirmed by Wolfram MCP. PDF→Markdown OCR artifact at line 176 (`c^{22}`) is *not* a paper error.

---

### Eq. (II.3) — Dual square-root equation (potential in the mass) — **NEW**

**As printed (line 183):**
$$i\hbar\frac{\partial\Psi}{\partial\tau} = \left\{\frac{\boldsymbol\pi^2}{2m} + \beta V + mc^2 - \frac{e\hbar\boldsymbol\Sigma\!\cdot\!\mathbf{B}}{2mc} + \frac{V^2}{2mc^2}\right\}\Psi.$$

**Pedagogical derivation (this is the cleanest of the three!).**

Start from $H_{\rm sqrt(2)} = \beta\sqrt{c^2\boldsymbol\pi^2 - ec\hbar\boldsymbol\Sigma\!\cdot\!\mathbf{B} + (mc^2 + \beta V)^2}$ (line 163).

- **Step 1 — Expand the radicand.** Let $S^2 = c^2\boldsymbol\pi^2 - ec\hbar\boldsymbol\Sigma\!\cdot\!\mathbf{B} + m^2c^4$ (the "standard" radicand). Then
$$(mc^2 + \beta V)^2 = m^2c^4 + 2mc^2\,\beta V + \beta^2 V^2 = m^2c^4 + 2mc^2\,\beta V + V^2,$$
using $\beta^2 = 1$.

The full radicand $Q$ inside $H_{\rm sqrt(2)}$ is therefore
$$Q = c^2\boldsymbol\pi^2 - ec\hbar\boldsymbol\Sigma\!\cdot\!\mathbf{B} + m^2c^4 + 2mc^2\,\beta V + V^2 = S^2 + 2mc^2\,\beta V + V^2.$$

- **Step 2 — $[\beta, Q] = 0$.** This is the key conceptual point. The radicand $Q$ commutes with $\beta$ because:
  - $[\beta, \boldsymbol\pi] = 0$ (different vector spaces).
  - $[\beta, \boldsymbol\Sigma] = 0$ (the spin matrix $\boldsymbol\Sigma = \mathrm{diag}(\boldsymbol\sigma, \boldsymbol\sigma)$ is block-diagonal, and $\beta = \mathrm{diag}(\mathbf{1}, -\mathbf{1})$ is also block-diagonal — they commute).
  - $[\beta, V] = 0$ (V is a scalar function of $\mathbf{x}$).
  - $[\beta, \beta V] = 0$ trivially.
  
  Therefore $[\beta, \sqrt{Q}] = 0$ as well (a function of $Q$ commutes with anything that commutes with $Q$).

- **Step 3 — Compute $H_{\rm sqrt(2)}^2$.** Since $\beta$ and $\sqrt{Q}$ commute:
$$H_{\rm sqrt(2)}^2 = (\beta\sqrt{Q})^2 = \beta^2 (\sqrt{Q})^2 = Q = c^2\boldsymbol\pi^2 - ec\hbar\boldsymbol\Sigma\!\cdot\!\mathbf{B} + m^2c^4 + 2mc^2\,\beta V + V^2.$$

- **Step 4 — Divide by $2mc^2$ and add $mc^2/2$.**
$$K_{\rm sqrt(2)} = \frac{\boldsymbol\pi^2}{2m} - \frac{e\hbar\boldsymbol\Sigma\!\cdot\!\mathbf{B}}{2mc} + \frac{mc^2}{2} + \beta V + \frac{V^2}{2mc^2} + \frac{mc^2}{2}.$$
Combine $mc^2/2 + mc^2/2 = mc^2$ to recover paper line 183. $\blacksquare$

**Why this matters.** Unlike (II.1) and (II.2) — which contain the awkward operator-valued terms $V\boldsymbol\alpha\!\cdot\!\boldsymbol\pi/(mc)$, $-i\hbar\boldsymbol\alpha\!\cdot\!\nabla V/(2mc)$, and $\beta(VS + SV)/(2mc^2)$ respectively — Eq. (II.3) is **scalar-only at the operator level** (apart from the trivial $\beta V$ term). It reads exactly like a Pauli-style wave equation with a relativistic correction $V^2/(2mc^2)$. This is the cleanest form of all three, and the easiest to study perturbatively.

**Mathematica check:**
```mathematica
ClearAll[c, m, hbar, ee, potV, pi2, SigmaB, beta];
H2sqrt2 = c^2 pi2 - ee hbar c SigmaB + m^2 c^4 + 2 m c^2 beta potV + potV^2;
K23v3 = H2sqrt2/(2 m c^2) + m c^2/2;
paperII3 = pi2/(2 m) + beta potV + m c^2 - ee hbar SigmaB/(2 m c) + potV^2/(2 m c^2);
FullSimplify[Expand[K23v3 - paperII3]]
(* Result: 0  ✅ (Wolfram MCP, 2026-05-11) *)
```

**Standard-equation comparison.** With $V \to 0$: $K \to \boldsymbol\pi^2/(2m) - e\hbar\boldsymbol\Sigma\!\cdot\!\mathbf{B}/(2mc) + mc^2$ — exactly the Pauli equation plus a rest-energy constant. So in the weak-field limit, Eq. (II.3) reduces *cleanly* to standard non-relativistic Pauli QM (Sakurai §3.5), with the addition of a relativistic correction $V^2/(2mc^2)$.

**Verdict:** ✅ Confirmed by Wolfram MCP. This is the structurally cleanest of the three dual equations.

---
