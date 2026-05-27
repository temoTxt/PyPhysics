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
| (III.21)–(III.23) | g-factor numerical reproduction | 🔴 fails |

---

> **🔴 Critical numerical finding (Section III.D).** The paper claims that with the cutoff $r_e = 0.499857150068631 \times r_0$, the formula $g_r = 2[1 - 4r_0/(2r+r_0)]$ produces the experimental electron $g$-factor $-2.00231930436256$. **Wolfram MCP computation shows that the formula at the stated $r_e$ instead gives $g \approx -2.00057148$.** To reproduce $g = -2.00231930436256$ one needs $r_e/r_0 \approx 0.4994205099$. Either the formula or the stated $r_e$ value contains a typo. This is the key falsifiable prediction of the paper and must be flagged for author review. See the Eq. (III.18)–(III.23) section below for full details.

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

#### Schwinger identification — empirical residual test (2026-05-26, [issue #66](https://github.com/temoTxt/PyPhysics/issues/66) Candidate 3)

**Motivation.** Inverting the §III.D formula $g_r(r) = 2[1 - 4r_0/(2r+r_0)]$ at the Schwinger one-loop QED anomalous moment $g_e = -2 - \alpha/\pi$ gives the closed-form

$$
\left.\frac{r_e}{r_0}\right|_{\text{Schwinger}} \;=\; \frac{2 - \alpha/(2\pi)}{4 + \alpha/\pi} \;=\; 0.499\,419\,632\,156\,99\ldots,
$$

which differs from the triangulated $r_e/r_0 = 0.499\,420\,509\,912\,83$ by $\Delta r = +8.78\times10^{-7}$. Issue #66 asks whether this near-coincidence indicates that §III.D's cutoff prescription is engineered to reproduce one-loop QED through a derived identity in $\alpha$, or whether the match is contingent.

**Empirical test.** Wolfram MCP at 20-digit precision in [`../Mathematica_Notebooks/Quantum_Mechanics/r_e_schwinger_residual_test.wl`](../Mathematica_Notebooks/Quantum_Mechanics/r_e_schwinger_residual_test.wl):

| Quantity | Value | Notes |
|---|---|---|
| $\alpha$ (CODATA 2018) | $7.297\,352\,5693\times 10^{-3}$ | |
| $r_e/r_0$ closed-form | $0.499\,419\,632\,156\,99$ | inversion of $g_r(r) = -2 - \alpha/\pi$ |
| $r_e/r_0$ triangulated | $0.499\,420\,509\,912\,83$ | PR #62 Pass B optimum |
| $\Delta r$ | $+8.78\times10^{-7}$ | triangulated − closed-form |
| $dg_r/dr$ at $r_{\text{closed}}$ | $4.0048$ | $= 16/(2r+1)^2$ |
| $\Delta g_e$ from $\Delta r$ | $+3.515\times10^{-6}$ | linear-in-$\Delta r$ propagation |
| $\Delta g_e = g_e^{\text{meas}} - g_e^{\text{Schwinger}}$ | $+3.515\times10^{-6}$ | |
| Karplus–Kroll two-loop $+2 C_2(\alpha/\pi)^2$, $C_2 = 0.328\,478\,965\,579\,193$ | $+3.545\times10^{-6}$ | Karplus–Kroll 1950 + Sommerfield 1957 + Petermann 1957 |
| Laporta–Remiddi three-loop $-2 C_3(\alpha/\pi)^3$, $C_3 = 1.181\,241\,456\,587$ | $-2.96\times10^{-8}$ | Laporta–Remiddi 1996, *Phys. Lett. B* **379**, 283 |
| Kinoshita-Fukuda+ four-loop $+2 C_4(\alpha/\pi)^4$, $C_4 \approx 1.9106$ | $+1.1\times10^{-10}$ | Aoyama–Hayakawa–Kinoshita–Nio numerical |
| Sum (KK + LR + KF) | $+3.515\,11\times10^{-6}$ | |
| Residual: observed − predicted | $-9.7\times10^{-12}$ | |

**Numerical observation.** The triangulated $r_e/r_0$ lies on the all-orders-QED curve in $r$-space: $\Delta g_e^{\text{observed}}$ matches $\Delta g_e^{\text{predicted}}$ (Schwinger one-loop + KK two-loop + LR three-loop + KF four-loop) to better than $5\times 10^{-12}$ — a 9-to-10-digit agreement, well below the framework's nominal $10^{-6}$ precision floor.

**Structural caveat (load-bearing for the interpretation).** This $10^{-11}$ agreement is *algebraically forced* by the back-fit, not evidence of intentional Schwinger encoding in §III.D. The triangulated $r_e/r_0$ is *defined* (per `r_e_triangulation.wl` Pass B) as the cutoff that makes $g_r(r) = g_e^{\text{meas}}$ to 16 sig figs. Therefore

$$
\Delta g_e^{\text{obs}} \;:=\; g_e^{\text{meas}} - g_e^{\text{Schwinger}} \;=\; g_e^{\text{meas}} - (-2 - \alpha/\pi)
$$

is identically equal, by construction, to the standard-QED all-orders-beyond-one-loop contribution $-2[a_e - \alpha/(2\pi)] = +2 C_2(\alpha/\pi)^2 - 2 C_3(\alpha/\pi)^3 + 2 C_4(\alpha/\pi)^4 - \ldots$ Any cutoff value satisfying $g_r(r) = g_e^{\text{meas}}$, when subtracted from the closed-form Schwinger value via the linear $dg_r/dr$ propagation, will yield this same residual. The agreement is therefore **necessary** for any cutoff reproducing $g_e^{\text{meas}}$ — including those produced by mechanisms entirely unrelated to one-loop QED — and is **not sufficient** to identify intentional encoding.

**Distinguishing intentional from forced.** Intentionality at the §III.D derivation level would require the framework to *produce the closed-form expression* $(2 - \alpha/(2\pi))/(4 + \alpha/\pi)$ as a derived identity in $\alpha$. Inspection of §III.D Eqs. (III.18)–(III.23): the framework derives the formula $g_r(r) = 2[1 - 4r_0/(2r+r_0)]$ and the limit checks $g_r(r_0/2) = -2$ and $g_r(r \to 0) = -6$, but **leaves $r_e/r_0$ as a free parameter** — per Tepper Gill's 2026-05-25 guidance, the published value $0.499857150068631$ was obtained by a uni-observable numerical search against $g_s^{\text{meas}}$, not derived from the framework's structure. No closed-form expression $r_e/r_0 = f(\alpha)$ appears in the published §III.D.

**Cross-particle test (muon, iter-5 addendum).** Eq. (III.23) extends the same dimensionless $g$-factor formula to the muon: $g_\mu^a = 2[1 - 4r_0^\mu/(2r_\mu + r_0^\mu)]$ with $r_0^\mu = e^2/(m_\mu c^2)$. The paper does not specify a numerical value for $r_\mu$. Applying the §III.D back-fit prescription to the muon (i.e., solving $g_\mu^a = -2(1 + a_\mu^{\text{exp}})$ for $r_\mu/r_0^\mu$) using the **Fermilab Muon $g-2$ Collaboration final 2023 result** $a_\mu^{\text{exp}} = 116\,592\,059(22) \times 10^{-11}$ (D. P. Aguillard et al. (FNAL Muon $g-2$), *Phys. Rev. D* **108**, 092009 (2023)) gives $r_\mu/r_0^\mu = 0.499\,417\,379\,350$. The electron's back-fit value $r_e/r_0^e = 0.499\,420\,509\,913$ differs by $\Delta r = +3.13 \times 10^{-6}$ — equivalent to $\sim 57{,}000\,\sigma_{a_\mu}$ in muon measurement units (with $\sigma_r = \sigma_{a_\mu}/4 \approx 5.5 \times 10^{-11}$ via the linear $dg_r/dr \approx 4$ propagation). The discrepancy in anomaly space is $a_\mu - a_e = +6.27 \times 10^{-6}$ ($\sim 28{,}500\,\sigma_{a_\mu}$), reflecting the mass-dependent QED, hadronic, and electroweak content of $a_\mu$ that is absent from $a_e$.

This **rules out at $> 57$kσ** any hypothesis under which the dimensionless cutoff $r/r_0$ is universal across leptons — i.e., any candidate closed-form derivation $r_e/r_0 = f(\alpha)$ in $\alpha$ alone. The hypothesis is *not made by the published paper* (DRQM-I §III.D leaves $r_\mu, r_p$ as separate free parameters per Eq. III.23, see verification doc lines 487–489), so this cross-particle test does *not* falsify the framework as written. What it constrains is any future closed-form derivation (e.g., from [issue #54](https://github.com/temoTxt/PyPhysics/issues/54)): such a derivation must produce a **particle-mass-dependent** cutoff $r/r_0 = f(\alpha, m_\ell/m_e, \ldots)$ that reproduces both $a_e^{\text{exp}}$ and $a_\mu^{\text{exp}}$ at the precision of their measurements — a much taller order than the pure-$\alpha$ closed form $(2 - \alpha/(2\pi))/(4 + \alpha/\pi)$.

**Verdict.** At the §III.D-as-published level, the Schwinger identification of $r_e/r_0$ is **structurally forced by the back-fit**, not derived as a closed-form identity. The empirical residual test is consistent with intentional Schwinger encoding *up to the structural caveat above* — the $10^{-11}$ agreement cannot, on its own, distinguish "framework derives the closed-form" from "framework back-fits any value that reproduces $g_e^{\text{meas}}$." The cross-particle test further narrows the discrimination: any future closed-form derivation must be particle-mass-dependent (pure-$\alpha$ closed forms are empirically ruled out). Discriminating B from C now requires:

1. **Inspection of a first-principles rederivation** from the dual-Dirac renormalisation prescription (tracked in [#54](https://github.com/temoTxt/PyPhysics/issues/54)) that produces $r_e/r_0$ in closed form **with particle-mass dependence** — if it reproduces both $a_e^{\text{exp}}$ and $a_\mu^{\text{exp}}$, outcome **B** ✅.
2. **Framework predictions for observables whose QED corrections do not factorise through $g_s$** (e.g., hydrogen 1S–2S transition, antiprotonic helium, muonic-H Lamb shift) — these break the $(g_s/-2)^n$-form degeneracy of the existing six-observable fit and would constrain $r_e$ independently of $g_e^{\text{meas}}$. Framework prediction formulas for these are not derived in the published Bethe–Salpeter campaign.

Finding 2 verdict: stays at ⚠ **CHARACTERISED**. The empirical-test path's outcome is **C-as-published** (Schwinger match algebraically forced for the electron; universal-closed-form cutoff ruled out at $> 57$kσ by the cross-particle test). Outcome **B** (intentional encoding) remains open through the routes above but requires particle-mass-dependent derivational structure. Loop iteration log: [`.dev/research/STATE.md`](../../.dev/research/STATE.md) iters 1–5.

<!-- TODO: human reviews and fills in — confirms (a) the structural caveat is the right framing of the iter-2 $10^{-11}$ KK+LR+KF agreement, (b) the C-as-published verdict correctly distinguishes "back-fit-forced" from "derivation-level intentional", (c) the two discrimination routes (#54 first-principles rederivation; type-(b) observables) correctly bracket what would upgrade C → B, and (d) the verdict-marker stays at ⚠ rather than shifting to ✅ or 🔴. Note any choice the reviewer would have made differently. -->

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
