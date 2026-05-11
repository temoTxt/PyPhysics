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
| (III.1) | Dirac eigenvalue problem (matrix form) | ⬜ |
| (III.2) | $\psi_2$ solution | ⬜ |
| (III.3) | $K_D$ in terms of $V_0, V$ | ⬜ |
| (III.4) | Dual Dirac eigenvalue equation | ⬜ |
| (III.5) | Coupled $\psi_1, \psi_2$ equations | ⬜ |
| (III.6) | Expanded eigenvalue equation | ⬜ |
| (III.7) | $\psi_2$ with denominator approximation | ⬜ |
| (III.8) | Approximated eigenvalue equation | ⬜ |
| (III.9)–(III.12) | Pauli identities | ⬜ |
| (III.13)–(III.17) | Spherical-coordinate algebra | ⬜ |
| (III.18) | $g_r = 2[1 - 4r_0/(2r+r_0)]$ result | 🔴 numerical |
| (III.19), (III.20) | Other new terms | ⬜ |
| (III.21)–(III.23) | g-factor formulas for $e^-$, $\mu$, $p$ | 🔴 numerical |

---

> **🔴 Critical numerical finding (Section III.D).** The paper claims that with the cutoff $r_e = 0.499857150068631 \times r_0$, the formula $g_r = 2[1 - 4r_0/(2r+r_0)]$ produces the experimental electron $g$-factor $-2.00231930436256$. **Wolfram MCP computation shows that the formula at the stated $r_e$ instead gives $g \approx -2.00057148$.** To reproduce $g = -2.00231930436256$ one needs $r_e/r_0 \approx 0.4994205099$. Either the formula or the stated $r_e$ value contains a typo. This is the key falsifiable prediction of the paper and must be flagged for author review. See the Eq. (III.18)–(III.23) section below for full details.

---

## Section I: Recap of dual classical theory

The bulk of Section I is a transcription of results already verified in [`Two_Mathematically_Equivalent_Versions_of_Maxwells_Equations.md`](./Two_Mathematically_Equivalent_Versions_of_Maxwells_Equations.md). Cross-references below.

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
