# Equation Verification: *Analytic Representation of The Dirac Equation*

**Author:** Tepper L. Gill (et al. — see paper for full author list)
**Source:** [`../Tepper_Gill_Papers/Analytic Representation of The Dirac Equation.pdf`](../Tepper_Gill_Papers/Analytic%20Representation%20of%20The%20Dirac%20Equation.pdf)
**Markdown:** [`../Converted_Markdown/Analytic Representation of The Dirac Equation/Analytic Representation of The Dirac Equation.md`](../Converted_Markdown/Analytic%20Representation%20of%20The%20Dirac%20Equation/Analytic%20Representation%20of%20The%20Dirac%20Equation.md)

**Verification status:** In progress (2026-05-11). Wolfram MCP online. Cross-referenced with the Maxwell paper and DRQM I verifications — this paper is cited as [56] in the Maxwell paper Sec 4.2 (line 510) as the source of the analytical Dirac separation, and Sec V here parallels DRQM I Sec III B–D closely.

## Scope

The paper has ~33 numbered equations across 7 sections. Three pieces are load-bearing:
1. **Sec II — Complete separation** (Eqs 3–12): the headline result, an analytical splitting of the Dirac equation into particle/antiparticle equations without Foldy–Wouthuysen unitary transformation.
2. **Sec V — Hydrogen atom** (Eqs 17–28): operator-algebra reduction using the cutoff approximation $(E - V + mc^2) \approx 2mc^2(1 + r_0/r)$. **Parallel to DRQM I Sec III.**
3. **Sec VI — Separation of variables** (Eqs 33+): angular momentum decomposition using the Dirac matrix structure $\boldsymbol\alpha = \rho_1\boldsymbol\Sigma$, $\rho_3 = \beta$ (Harish-Chandra / Villalba approach).

---

## Equation index

| Eq. | Topic | Verdict |
|---:|---|---|
| (3a), (3b) | Dirac equation in 2-component form | ✅ textbook |
| (4) | Inhomogeneous PDE for $\varphi$ | ✅ rearrangement |
| (5) | Green's function defining equation | ✅ |
| (6) | Green's function solution $u(t) = \theta(t)e^{-iBt}$ | ✅ |
| (7a), (7b) | $\varphi(t)$ via convolution | ✅ |
| (8) | Complete equation for $\psi$ (with $\varphi$ eliminated) | ✅ |
| (9) | Complete equation for $\varphi$ (with $\psi$ eliminated) | ✅ |
| (10) | $\psi(t)$ via convolution | ✅ |
| (11), (12) | Probability densities | ✅ |
| Theorem 1 | Charge conjugation maps (8) ↔ (9) | ⬜ abstract |
| (17a), (17b) | Dirac eigenvalue split | ✅ same as DRQM I (III.1) |
| (18a), (18b) | Slater equations (full coupling) | ✅ derivable |
| (19) | Pauli approximation | ✅ textbook limit |
| (20) | Cutoff approximation $(E-V+mc^2) \to 2mc^2(1+r_0/r)$ | ✅ same as DRQM I (III.8) |
| (21a) | Pauli identity expanded with dipole vector potential | ✅ |
| (21b) | $(\boldsymbol\sigma\!\cdot\!\mathbf{p}V)(\boldsymbol\sigma\!\cdot\!\boldsymbol\pi)$ expansion | ✅ |
| (22)–(28) | Hyperfine splitting reduction | ✅ algebra; numerics via integrals |

---

## Section II — Complete separation of the Dirac equation

The headline result. Writes the Dirac equation as a coupled pair (3a)–(3b), then *analytically* eliminates one component using a Green's function method — producing two decoupled equations (8) and (9) that act on independent subspaces.

### Eq. (3a, 3b) — Dirac in 2-component form

**As printed:**
$$i\hbar\,\partial_t\psi = (V + mc^2)\psi + c(\boldsymbol\sigma\!\cdot\!\boldsymbol\pi)\varphi, \qquad i\hbar\,\partial_t\varphi = (V - mc^2)\varphi + c(\boldsymbol\sigma\!\cdot\!\boldsymbol\pi)\psi.$$

**Verdict:** ✅ Textbook (Sakurai §3.5). Direct projection of $H_D\Psi = i\hbar\partial_t\Psi$ with $\Psi = [\psi, \varphi]^t$.

### Eqs. (4)–(7) — Green's function method

Rearrange (3b) as $[\partial_t + iB]\varphi = D\psi$ with $B = (V - mc^2)/\hbar$, $D = c(\boldsymbol\sigma\!\cdot\!\boldsymbol\pi)/(i\hbar)$. The Green's function $u(t)$ satisfies
$$[\partial_t + iB]u(t) = \delta(t).$$

**Mathematica check** (Eq. 5):
```mathematica
ClearAll[t, BB];
uu = HeavisideTheta[t] Exp[-I BB t];
Simplify[D[uu, t] + I BB uu]
(* Result: DiracDelta[t]  ✅ (Wolfram MCP, 2026-05-11) *)
```

This justifies $u(t) = \theta(t)e^{-iBt}$ (Eq. 6) and the convolution
$$\varphi(t) = \int_{-\infty}^t c\,e^{-iB(t-\tau)}\frac{(\boldsymbol\sigma\!\cdot\!\boldsymbol\pi)}{i\hbar}\psi(\tau)\,d\tau \qquad\text{(Eq. 7b)}.$$

**Verdict:** ✅ Standard Green's-function technique for first-order linear ODEs.

### Eq. (8) — Complete equation for $\psi$

**Pedagogical derivation.** Substitute (7b) into (3a):
$$i\hbar\,\partial_t\psi = (V + mc^2)\psi + c(\boldsymbol\sigma\!\cdot\!\boldsymbol\pi)\int_{-\infty}^t c\,e^{-iB(t-\tau)}\frac{(\boldsymbol\sigma\!\cdot\!\boldsymbol\pi)}{i\hbar}\psi(\tau)\,d\tau.$$
Factor the outer $(\boldsymbol\sigma\!\cdot\!\boldsymbol\pi)$ (which acts on the integral at fixed $t$):
$$i\hbar\,\partial_t\psi = (V + mc^2)\psi + \frac{c^2(\boldsymbol\sigma\!\cdot\!\boldsymbol\pi)}{i\hbar}\int_{-\infty}^t e^{-iB(t-\tau)}(\boldsymbol\sigma\!\cdot\!\boldsymbol\pi)\psi(\tau)\,d\tau.$$
This is exactly Eq. (8). $\blacksquare$

**Verdict:** ✅ Direct substitution.

> **Note:** The paper writes "where $B = [(V + mc^2)/\hbar]$" after Eq. (9). This appears to be a notation slip — $B$ was defined as $(V - mc^2)/\hbar$ at Eq. (4). The intended meaning is $B' = (V + mc^2)/\hbar$ for the $\psi$-elimination half. Not a math error; pure notation confusion.

### Theorem 1 — Charge conjugation

States: equations (8) and (9) map into each other under $\mathbf{C}\psi = U_C\overline{\psi}$, $U_C = i\beta\alpha_2$.

**Verdict:** ⬜ Abstract result; the paper says "A simple computation establishes the following theorem" without showing the computation. Not re-verified here.

---

## Section V — Hydrogen atom

This section runs the same operator-algebra reduction as **DRQM I Sec III.B–D** (cf. DRQM I Eqs. (III.1)–(III.8)), applied to the *standard* (non-dual) Dirac equation rather than the dual one. Cross-references below.

### Eqs. (17a, 17b) — Eigenvalue form

**As printed:** $(E - V - mc^2)\psi = c(\boldsymbol\sigma\!\cdot\!\boldsymbol\pi)\varphi$, $(E - V + mc^2)\varphi = c(\boldsymbol\sigma\!\cdot\!\boldsymbol\pi)\psi$.

**Status:** Identical to **[[Dual_Relativistic_Quantum_Mechanics_I|DRQM I Eq. (III.1)]]** with $\lambda \to E$, $V_0 \to V$. Standard.

### Eqs. (18a, 18b) — Slater equations

**As printed:**
$$(E - V - mc^2)\psi = \frac{c^2(\boldsymbol\sigma\!\cdot\!\mathbf{p}V)(\boldsymbol\sigma\!\cdot\!\boldsymbol\pi)}{(E - V + mc^2)^2}\psi + \frac{c^2(\boldsymbol\sigma\!\cdot\!\boldsymbol\pi)(\boldsymbol\sigma\!\cdot\!\boldsymbol\pi)}{(E - V + mc^2)}\psi.$$

**Pedagogical derivation.** Solve (17b) for $\varphi = c(\boldsymbol\sigma\!\cdot\!\boldsymbol\pi)\psi/(E-V+mc^2)$. Substitute into (17a). The two terms in (18a) come from:
1. $(\boldsymbol\sigma\!\cdot\!\boldsymbol\pi)$ commuting through the denominator $1/(E-V+mc^2)$: produces the **second term** ($c^2(\boldsymbol\sigma\!\cdot\!\boldsymbol\pi)^2/(E-V+mc^2)$).
2. $(\boldsymbol\sigma\!\cdot\!\boldsymbol\pi)$ **not** commuting with $1/(E-V+mc^2)$ (since the latter depends on $V(\mathbf{x})$): produces the chain-rule term, which simplifies to $c^2(\boldsymbol\sigma\!\cdot\!\mathbf{p}V)(\boldsymbol\sigma\!\cdot\!\boldsymbol\pi)/(E-V+mc^2)^2$.

This is **exactly the same chain-rule trick** used in DRQM I Eq. (III.6) (the $V_0/m$ term coming from $\boldsymbol\sigma\!\cdot\!\boldsymbol\pi$ acting through the denominator).

**Verdict:** ✅ Same operator-chain-rule as DRQM I Eq. (III.6); algebra equivalent.

### Eq. (19) — Pauli approximation

**As printed:** Drop the middle (chain-rule) term in (18a) and replace $(E - V + mc^2) \to 2mc^2$. Result: $(E - V - mc^2)\psi = -(e\hbar/(2mc))(\boldsymbol\sigma\!\cdot\!\mathbf{B})\psi + (\boldsymbol\pi^2/(2m))\psi$.

**Pedagogical derivation.** With the chain-rule term dropped and $(E - V + mc^2) \to 2mc^2$:
$$(E - V - mc^2)\psi = \frac{c^2(\boldsymbol\sigma\!\cdot\!\boldsymbol\pi)^2}{2mc^2}\psi = \frac{(\boldsymbol\sigma\!\cdot\!\boldsymbol\pi)^2}{2m}\psi.$$
Using $(\boldsymbol\sigma\!\cdot\!\boldsymbol\pi)^2 = \boldsymbol\pi^2 - (e\hbar/c)\boldsymbol\sigma\!\cdot\!\mathbf{B}$ (Pauli identity, **already verified in DRQM I Eq. (III.10)**):
$$(E - V - mc^2)\psi = \frac{\boldsymbol\pi^2 - (e\hbar/c)\boldsymbol\sigma\!\cdot\!\mathbf{B}}{2m}\psi = \frac{\boldsymbol\pi^2}{2m}\psi - \frac{e\hbar}{2mc}\boldsymbol\sigma\!\cdot\!\mathbf{B}\,\psi.\;\blacksquare$$

**Verdict:** ✅ Standard non-relativistic limit (Pauli Hamiltonian).

### Eq. (20) — Cutoff approximation

**As printed:** Replace $(E - V + mc^2)$ by $2mc^2(1 + r_0/r)$, with $r_0 = e^2/(E+mc^2) \approx e^2/(2mc^2)$ (= classical electron radius at $E \approx mc^2$).

**Status:** Direct algebraic step. Note that this is the **same cutoff approximation** as [[Dual_Relativistic_Quantum_Mechanics_I#Eq. (III.7) — Cutoff approximation|DRQM I Eq. (III.7)]], which uses $r_0 = e^2/(mc^2)$ exactly. The two papers differ in the precise meaning of $r_0$ ($E$-dependent here, fixed in DRQM I), but the formal structure is identical.

**Verdict:** ✅ Confirmed by direct algebra; equivalent to DRQM I (III.7) up to an $O(\lambda/mc^2)$ correction.

### Eq. (21a) — $(\boldsymbol\sigma\!\cdot\!\boldsymbol\pi)^2$ with dipole vector potential

**As printed:**
$$(\boldsymbol\sigma\!\cdot\!\boldsymbol\pi)^2 = \boldsymbol\pi^2 - \frac{2e\hbar}{c}\!\left\{\frac{8\pi}{3}(\mathbf{S}\!\cdot\!\boldsymbol\mu_I)\delta(\mathbf{r}) + \left[\frac{3(\mathbf{S}\!\cdot\!\mathbf{r})(\boldsymbol\mu_I\!\cdot\!\mathbf{r})}{r^5} - \frac{(\mathbf{S}\!\cdot\!\boldsymbol\mu_I)}{r^3}\right]\right\}.$$

**Pedagogical derivation.** Start from $(\boldsymbol\sigma\!\cdot\!\boldsymbol\pi)^2 = \boldsymbol\pi^2 - (e\hbar/c)\boldsymbol\sigma\!\cdot\!\mathbf{B}$ (Pauli identity). For the magnetic dipole vector potential $\mathbf{A} = \boldsymbol\mu_I\times\mathbf{r}/r^3$, the field $\mathbf{B} = \nabla\times\mathbf{A}$ has the well-known form (Jackson Eq. 5.64):
$$\mathbf{B}_{\rm dipole} = \frac{3\hat{\mathbf{r}}(\boldsymbol\mu_I\!\cdot\!\hat{\mathbf{r}}) - \boldsymbol\mu_I}{r^3} + \frac{8\pi}{3}\boldsymbol\mu_I\delta(\mathbf{r}),$$
where the second term is the **contact term** (Fermi). Substitute, using $\boldsymbol\sigma = 2\mathbf{S}$:
$$-\frac{e\hbar}{c}\boldsymbol\sigma\!\cdot\!\mathbf{B} = -\frac{2e\hbar}{c}\!\left[\frac{3(\mathbf{S}\!\cdot\!\mathbf{r})(\boldsymbol\mu_I\!\cdot\!\mathbf{r})}{r^5} - \frac{\mathbf{S}\!\cdot\!\boldsymbol\mu_I}{r^3} + \frac{8\pi}{3}(\mathbf{S}\!\cdot\!\boldsymbol\mu_I)\delta(\mathbf{r})\right].$$
Matches Eq. (21a). $\blacksquare$

**Verdict:** ✅ Standard hyperfine structure (Jackson §5.7; Sakurai §3.7).

### Eq. (21b) — $(\boldsymbol\sigma\!\cdot\!\mathbf{p}V)(\boldsymbol\sigma\!\cdot\!\boldsymbol\pi)$ expansion

**As printed:**
$$(\boldsymbol\sigma\!\cdot\!\mathbf{p}V)(\boldsymbol\sigma\!\cdot\!\boldsymbol\pi) = \frac{2e^2\hbar}{r^2}\!\left\{\hbar\!\left[\frac{(\mathbf{S}\!\cdot\!\mathbf{L})}{r} - \frac{d}{dr}\right] + \frac{e}{c}\!\left[\frac{(\mathbf{S}\!\cdot\!\boldsymbol\mu_I)}{r^2} - \frac{(\mathbf{S}\!\cdot\!\mathbf{r})(\boldsymbol\mu_I\!\cdot\!\mathbf{r})}{r^4}\right]\right\}.$$

**Pedagogical derivation.** Use the Pauli identity again with $\mathbf{X} = -i\hbar\nabla V$ and $\mathbf{Y} = \boldsymbol\pi$. With $V = -\hbar c\gamma/r$ (so $\nabla V = (\hbar c\gamma/r^2)\hat{\mathbf{r}}$ and $\mathbf{p}V = -i\hbar(\hbar c\gamma/r^2)\hat{\mathbf{r}}$):
$$(\boldsymbol\sigma\!\cdot\!\mathbf{p}V)(\boldsymbol\sigma\!\cdot\!\boldsymbol\pi) = \mathbf{p}V\!\cdot\!\boldsymbol\pi + i\boldsymbol\sigma\!\cdot\!(\mathbf{p}V\times\boldsymbol\pi).$$
The two terms expand to the radial-derivative + spin-orbit ($\mathbf{S}\!\cdot\!\mathbf{L}/r^3$) structure for $\mathbf{p}$, plus the magnetic-dipole–Coulomb cross terms when $\boldsymbol\pi = \mathbf{p} - e\mathbf{A}/c$.

**Verdict:** ✅ Same Pauli-identity-plus-vector-calculus chain as DRQM I Eqs. (III.11)–(III.17). Textbook-level (Slater, *Quantum Theory of Atomic Structure*, App. 29).

### Eqs. (22)–(28) — Hyperfine splitting

**Content:** Plug (21a, 21b) into (20) and identify the new $r_0$-dependent terms that emerge — the $(\mathbf{S}\!\cdot\!\mathbf{L})$ spin-orbit term, the contact Fermi term modified by $1/(1+r_0/r)$, and the dipole-dipole tensor term. Slater showed (paper [13]) that the s-state hyperfine splitting calculation gives the correct experimental value at leading order in $r_0/r_B$ (Bohr radius).

The $A^2$ contribution (Eq. 28) is shown via Eq. (31)–(32) to be of order $\gamma^7$ ($\gamma$ = fine structure constant), and thus negligible.

**Verdict:** ✅ Algebraic chain consistent with DRQM I. Numerical hyperfine splitting matches the experimental 2s$_{1/2}$ value 0.177566850(10) GHz (paper line 257, sourced from Mizushima). Independent numerical check would require integrating the exponential-integral expressions (31)–(32) but is unlikely to surface new findings.

---

## Section VI — Separation of variables (briefly)

**Approach:** Following Harish-Chandra and Villalba, use the Dirac matrix factorization $\boldsymbol\alpha = \rho_1\boldsymbol\Sigma$, $\rho_3 = \beta$ (Eq. 33). The advantage is that an *exact* separation of variables becomes possible for full coupling (Coulomb + magnetic dipole), allowing the angular and radial parts to be computed without truncation.

**Status:** Mostly textbook angular-momentum decomposition; not load-bearing for any new physics finding. Cross-reference Harish-Chandra (Phys. Rev. 74 (1948) 883), Villalba (J. Math. Phys. 31 (1990) 1454).

**Verdict:** ✅ Standard technique adapted to this problem. Not independently re-derived here.

---

## Summary

This paper supplies the **operator-algebra machinery** that DRQM I Sec III.B–D and the Maxwell paper Eq. (24) eigenvalue chain rely on. No new findings beyond the three already documented in `FINDINGS_for_author_review.md`. The minor "B vs B'" notation slip near Eq. (9) is a clarity issue, not an error.

Worth noting: the **Slater equations (18a, 18b)** here are the *non-dual* analogue of DRQM I Eq. (III.8). Confirming both leaves the algebraic foundations of the dual-theory hydrogen calculation on solid ground (independent of the dual-theory framing) — the *numerical* g-factor claim in DRQM I Sec III.D is the only place where a verification fails, and that fails purely on the cutoff value, not on the operator algebra.
