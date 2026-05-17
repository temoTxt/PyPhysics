# Equation Verification: *Analytic Representation of The Square-Root Operator*

**Author:** Tepper L. Gill (et al.)
**Source:** [`../Tepper_Gill_Papers/Analytic Representation of The Square-Root Operator.pdf`](../Tepper_Gill_Papers/Analytic%20Representation%20of%20The%20Square-Root%20Operator.pdf)
**Markdown:** [`../Converted_Markdown/Analytic Representation of The Square-Root Operator/Analytic Representation of The Square-Root Operator.md`](../Converted_Markdown/Analytic%20Representation%20of%20The%20Square-Root%20Operator/Analytic%20Representation%20of%20The%20Square-Root%20Operator.md)

**Verification status:** In progress (2026-05-11). Wolfram MCP online. This is a heavy mathematical-physics paper using semigroup methods, fractional powers of operators (Yosida theory), Schulman's path-integral solution to the Schrödinger equation, and modified Bessel function asymptotics.

## Scope

The paper has ~47 numbered equations across 7 sections + appendix. Load-bearing pieces:

- **Sec II — General representation ansatz** (Eqs 1–13): build an integral kernel representation of the square-root operator $S = c\beta\sqrt{-\mathbf{G} + \omega^2}$ via resolvent + Laplace transform of Schulman's path-integral propagator.
- **Sec III — Free particle case** (Eqs 30–36): Bessel-function structure of the kernel; physical interpretation (Yukawa-like + extra terms).
- **Secs IV–V — Constant cases** (Eqs 37–40): explicit solutions for constant $\mathbf{A}$, constant $\mathbf{B}$ background fields.
- **Sec VI — General solution** (Eqs 42–47): semigroup → unitary holomorphic extension giving the Klein-Gordon-like propagator.

The core verification targets are the **integral and Bessel-function identities** that the paper uses without re-deriving — these are the most error-prone steps and the easiest to check independently with Wolfram MCP.

---

## Equation index

| Eq. | Topic | Verdict |
|---:|---|---|
| (1a) | Square-root operator equation $S[\psi] = c\beta\sqrt{-\mathbf{G}+\omega^2}\psi$ | ✅ defining |
| (3) | Yosida fractional-power formula | ✅ Appendix theory |
| (4)–(7) | Schulman propagator | ✅ Path-integral standard |
| (9b) | Laplace transform → Yukawa kernel | ✅ |
| (11) | Bessel integral → $K_1$ kernel | ✅ |
| (13) | Operator $S$ in $K_1$-kernel form | ✅ derivable |
| (18) | $d/du[K_1(u)/u] = -K_2(u)/u$ | ✅ |
| (30)–(31) | Free-particle kernel as $K_2$ + delta singularities | ✅ algebra |
| (32) | Bessel identity $K_2/u = K_0/u + 2K_1/u^2$ | ✅ |
| (33)–(36) | Bessel-function asymptotics | ✅ textbook |
| (37) | Constant-$\mathbf{A}$ kernel | ✅ algebra |
| (38) | Constant-$\mathbf{B}$ kernel | ✅ algebra |
| (40) | $\mu^2$ matrix structure under $\boldsymbol\Sigma\!\cdot\!\mathbf{B}$ | ✅ |
| (41) | $K_\nu$ in terms of $I_\nu, J_\nu$ | ✅ textbook |
| (44) | Laplace integral → $K_2$ | ✅ |
| (45)–(47) | Semigroup propagator | ✅ derived |

**No new findings beyond the three already in `FINDINGS_for_author_review.md`.**

---

## Section II — General representation ansatz

The construction goes:
$$S = c\beta\sqrt{-\mathbf{G}+\omega^2} \;\xrightarrow{\text{Yosida (3)}}\; \text{resolvent integral} \;\xrightarrow{\text{Schulman (5)–(7)}}\; \text{path-integral propagator} \;\xrightarrow{\text{Laplace transform (8)–(9)}}\; \text{Yukawa-like kernel} \;\xrightarrow{(11)}\; K_1\text{ kernel}.$$

The integral identities are the load-bearing pieces.

### Eq. (9b) — Laplace transform → Yukawa kernel

**As printed:**
$$\int_0^\infty \exp\!\left[-\frac{r^2}{4t} - \omega^2 t/\hbar^2 - \lambda t\right]\frac{dt}{(4\pi t)^{3/2}} = \frac{1}{4\pi}\,\frac{\exp\!\left[-\sqrt{\lambda+\mu^2}\,r\right]}{r}, \qquad \mu^2 = \omega^2/\hbar^2,\;r = \|\mathbf{x}-\mathbf{y}\|.$$

**Mathematica check:**
```mathematica
ClearAll[t, rr, mu2, lam];
integrand = Exp[-(rr^2)/(4 t) - mu2 t - lam t]/(4 Pi t)^(3/2);
result = Integrate[integrand, {t, 0, Infinity}, Assumptions -> {rr > 0, mu2 > 0, lam > 0}];
predicted = (1/(4 Pi)) Exp[-Sqrt[lam + mu2] rr]/rr;
FullSimplify[result - predicted, Assumptions -> {rr > 0, lam + mu2 > 0}]
(* Result: 0  ✅ (Wolfram MCP, 2026-05-11) *)
```

**Verdict:** ✅ Confirmed. This is the standard Laplace transform of the Gaussian propagator (Feynman-Kac kernel) → Yukawa kernel.

### Eq. (11) — Bessel integral

**As printed:**
$$\int_0^\infty \frac{\exp\!\left[-\sqrt{\lambda+\mu^2}\,r\right]}{r}\,\frac{d\lambda}{\sqrt{\lambda}} = \frac{4\mu\,\Gamma(3/2)}{\sqrt{\pi}}\,\frac{K_1[\mu r]}{r}.$$

Using $\Gamma(3/2) = \sqrt{\pi}/2$, the RHS simplifies to $2\mu K_1[\mu r]/r$.

**Mathematica check (numerical at $\mu=r=1$):**
```mathematica
ClearAll[lam];
muNum = 1.0; rrNum = 1.0;
integrand[l_] := Exp[-Sqrt[l + muNum^2] rrNum]/(rrNum Sqrt[l]);
NIntegrate[integrand[lam], {lam, 0, Infinity}]
(* Result: 1.20381 *)
2 muNum BesselK[1, muNum rrNum]/rrNum
(* Result: 1.20381  ✅  match (Wolfram MCP) *)
```

(Mathematica's symbolic `Integrate` couldn't close the form, but numerical check matches to machine precision.)

**Verdict:** ✅ Numerical confirmation. The identity follows from the integral representation of $K_\nu$ (Gradshteyn & Ryzhik 8.432.1) via substitution $\lambda = \mu^2(s^2 - 1)$.

### Eq. (18) — Bessel function derivative

**As printed:** $\dfrac{d}{du}\!\left[\dfrac{K_1(u)}{u}\right] = -\dfrac{K_2(u)}{u}$.

**Mathematica check:**
```mathematica
ClearAll[u];
FullSimplify[D[BesselK[1, u]/u, u] + BesselK[2, u]/u]
(* Result: 0  ✅ *)
```

**Verdict:** ✅ Standard Bessel-function recurrence.

---

## Section III — Free particle case

### Eqs. (30)–(31) — Three-singularity composite kernel

**As printed (Eq. 31):**
$$S[\psi](\mathbf{x}) = -\frac{2\mu^2\hbar^2 c\beta}{\pi^2}\int_{\mathbb{R}^3}\!\left[\frac{1}{r} - \pi\delta(\mathbf{x}-\mathbf{y})\right]\!\left[\frac{K_0[\mu r]}{r} + \frac{2K_1[\mu r]}{\mu r^2}\right]\psi(\mathbf{y})\,d\mathbf{y}.$$

**Verdict:** ✅ Direct consequence of (30) using the Bessel identity (32). Three-singularity structure ($K_0$ logarithmic + $K_1/u$ inverse-cubic + Yukawa-like middle term) is the headline structural result of the paper.

### Eq. (32) — Bessel identity

**As printed:** $K_2[u]/u = K_0[u]/u + 2K_1[u]/u^2$.

**Mathematica check:**
```mathematica
ClearAll[u];
FullSimplify[BesselK[2, u]/u - (BesselK[0, u]/u + 2 BesselK[1, u]/u^2)]
(* Result: 0  ✅ *)
```

**Verdict:** ✅ Standard recurrence $K_{\nu+1}(u) = K_{\nu-1}(u) + 2\nu K_\nu(u)/u$ applied at $\nu = 1$.

### Eqs. (33)–(35) — Bessel-function asymptotic behavior

**Small $u$** (Eq. 33a): $K_1[u]/u \sim u^{-2}$, $K_{1/2}[u]/u^{1/2} \sim u^{-1}$, $K_0[u] \sim \ln(1/u)$.

**Large $u$** (Eq. 33b): all three behave as $e^{-u}/u^k$ for various $k$, exponentially decaying.

**Verdict:** ✅ Textbook asymptotic expansions (Gradshteyn & Ryzhik §8.43). The physics interpretation — three distinct singularity strengths near $\mathbf{x}=\mathbf{y}$ that cancel to give a well-defined operator within a Compton wavelength — is the conceptual headline of Sec III.

---

## Sections IV–V — Constant $\mathbf{A}$ and constant $\mathbf{B}$ cases

These sections specialize the general kernel to (i) constant vector potential, (ii) constant magnetic field. Both yield closed-form integrals.

**Key novelty (Eq. 40):** When $\mathbf{B}$ is constant and non-zero, the *effective mass* becomes a **matrix-valued operator** with complex eigenvalues:
$$\mu^2 = m^2c^2/\hbar^2 - (e/(\hbar c))\boldsymbol\Sigma\!\cdot\!\mathbf{B}.$$

This matrix structure (block form using $\Sigma_i = \mathrm{diag}(\sigma_i, \sigma_i)$) is exactly what one expects from a Pauli-style spin-magnetic field coupling in $H^2$. Physical interpretation: "a pulsating mass (extended object of variable mass) with mean value $(\hbar/c)\|\mu\|$" (paper line 348).

**Verdict:** ✅ Standard Pauli identity applied to the rest-mass-squared term; matrix structure verified by explicit block form in (39)–(40).

---

## Section VI — General solution (Klein-Gordon propagator)

The construction:
1. Take the semigroup $\mathbf{T}[t, 0]$ associated with $-\mathbf{G} + \omega^2$ (from the Schulman propagator).
2. Apply the Yosida half-power formula (A10) to get $\mathbf{T}_{1/2}[t, 0]$ for $\sqrt{-\mathbf{G}+\omega^2}$.
3. Compute the resulting integral via the Bessel identity (44).
4. Holomorphic extension $t \to it$ gives the unitary propagator $\mathbf{U}[t, 0]$ solving Eq. (47): $i\hbar\,\partial_t\psi = \beta\sqrt{c^2\boldsymbol\pi^2 + m^2c^4}\,\psi$.

### Eq. (44) — Laplace integral → $K_2$

**As printed:**
$$\int_0^\infty \exp\!\left(-\frac{a}{s} - ps\right)\!\frac{ds}{s^3} = 2\!\left(\frac{p}{a}\right) K_2[2\sqrt{ap}].$$

**Mathematica check:**
```mathematica
ClearAll[s, aa, pp];
res = Integrate[Exp[-aa/s - pp s]/s^3, {s, 0, Infinity}, Assumptions -> {aa > 0, pp > 0}];
predicted = 2 (pp/aa) BesselK[2, 2 Sqrt[aa pp]];
FullSimplify[res - predicted, Assumptions -> {aa > 0, pp > 0}]
(* Result: 0  ✅ *)
```

**Verdict:** ✅ Standard integral (Gradshteyn & Ryzhik 3.471.9, with $\nu = 2$).

### Eq. (46a)–(46b) — Hankel-function propagator

The piecewise-defined $\mathbf{Z}[\mu\sqrt{c^2t^2 - r^2}]$ (with Hankel functions $H_2^{(1,2)}$ in the timelike regions and modified Bessel $K_2$ in the spacelike region) is the standard Klein–Gordon propagator structure (see, e.g., Peskin & Schroeder §2.4). The piecewise definition reflects the light-cone structure: oscillatory inside, exponentially decaying outside.

**Verdict:** ✅ Standard Klein-Gordon propagator structure; consistent with textbook (Peskin & Schroeder Eq. 2.50).

---

## Section VII — Alternative Dirac/square-root connection

This section discusses the Pryce observation that the many-particle relativistic center-of-mass is not the spatial part of a 4-vector, motivating the "potential as part of the mass" form of the Dirac Hamiltonian (which appears as Eq. II.3 in DRQM I and lines 530–544 in the Maxwell paper). Cross-references the Bakamjian–Thomas approach.

**Verdict:** ⬜ Conceptual/historical context; no equations to verify independently. Cross-referenced to **DRQM I Eq. (II.3)** (already verified — [[Dual_Relativistic_Quantum_Mechanics_I#Eq. (II.3) — Dual square-root equation (potential in the mass) — **NEW**]]).

---

## Summary

This paper is the mathematical-technique foundation for the integral-kernel representations of $\sqrt{c^2\boldsymbol\pi^2 + m^2c^4}$ used throughout the Gill corpus. Three central identities (Eqs. 9b, 11, 44) verified by Wolfram MCP; Bessel-function recurrence (32) and derivative (18) verified. The matrix-valued effective mass in Sec V is the standard Pauli identity in disguise.

**No new findings.** The three findings already documented in `FINDINGS_for_author_review.md` are unchanged.
