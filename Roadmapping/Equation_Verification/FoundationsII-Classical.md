# Equation Verification: *FOUNDATIONS FOR QED II: CLASSICAL THEORY*

**Authors:** Tepper L. Gill, Gonzalo Ares de Parga
**Source:** [`../Tepper_Gill_Papers/FoundationsII-Classical.pdf`](../Tepper_Gill_Papers/FoundationsII-Classical.pdf)
**Markdown:** [`../Converted_Markdown/FoundationsII-Classical/FoundationsII-Classical.md`](../Converted_Markdown/FoundationsII-Classical/FoundationsII-Classical.md)

**Verification status:** In progress (2026-05-11). Wolfram MCP online.

## Scope

Companion to **FOUNDATIONS FOR QED I** focusing on the classical theory. Sections:
- **Sec 2** — Einstein-Newtonian extension; three K forms (overlap with TCEP 5.4, 5.6, 5.7); **Einstein-Newtonian particle** with the load-bearing critical-point analysis at $r = r_0$.
- **Sec 3** — Einstein and dual theory (overlap with Maxwell paper and DRQM I).
- **Sec 4** — Light: nature, mass, speed (philosophical/historical).
- **Sec 5** — Historical time, CMBR, symmetric big bang (speculative cosmology).

Sec 2.2 contains a result that **resolves Open Question 3** from the Maxwell paper.

---

## Equation index

| Eq. | Topic | Verdict |
|---:|---|---|
| (2.1)–(2.4) | Proper time, dual identities | ✅ (= Maxwell Eqs 1–2) |
| Thm 1 | Minkowski Incompatible Theorem | ✅ (= QED I Thm 1.1) |
| Thm 2 | Einstein Compatibility Theorem (1) | ✅ |
| Thm 3 | Einstein Dual Compatibility Theorem | ✅ |
| (2.5) | Chain rule for $dW/d\tau$ | ✅ (= TCEP Eq 5.2) |
| (2.6) | $K = H^2/(2mc^2) + mc^2/2 = \mathbf{p}^2/(2m) + mc^2$ | ✅ |
| (2.7) | General solution | ✅ |
| (2.8) | Fixed-Lorentz-frame K | ✅ (= TCEP 5.6) |
| (2.9) | Fixed-momentum K | ✅ (= TCEP 5.7) |
| **(2.10)** | **K with H = √(c²p² + (mc² + V)²) — V as part of mass** | **✅ + resolves Open Q3** |
| **(2.11)** | **F_K = −∇V(1 + V/mc²), critical point at $r = r_0$** | **✅ resolves Open Q3** |
| (2.12)–(2.13) | Hamilton equations for "V as additive" case | ✅ (= DRQM I I.7) |
| Sec 3.x | Einstein + dual Maxwell theory | ✅ (cross-ref Maxwell paper) |

**No new findings beyond the three already documented.** Significant cross-reference: this paper provides the **explicit derivation** of the $r_0$ critical-point claim that was flagged as Open Question 3 in the Maxwell paper.

---

## Section 2.2 — Einstein-Newtonian Particle (the load-bearing novel content)

### Eq. (2.10) — K when V appears in the mass

**As printed:**
$$H = \sqrt{c^2\mathbf{p}^2 + (mc^2 + V)^2}, \qquad K = \frac{\mathbf{p}^2}{2m} + V + \frac{V^2}{2mc^2} + mc^2.$$

**Pedagogical derivation.**
- **Step 1.** Compute $H^2 = c^2\mathbf{p}^2 + (mc^2 + V)^2 = c^2\mathbf{p}^2 + m^2c^4 + 2mc^2 V + V^2$.
- **Step 2.** $K = H^2/(2mc^2) + mc^2/2$:
$$K = \frac{c^2\mathbf{p}^2}{2mc^2} + \frac{m^2c^4}{2mc^2} + \frac{2mc^2 V}{2mc^2} + \frac{V^2}{2mc^2} + \frac{mc^2}{2} = \frac{\mathbf{p}^2}{2m} + \frac{mc^2}{2} + V + \frac{V^2}{2mc^2} + \frac{mc^2}{2}.$$
- **Step 3.** Combine $\frac{mc^2}{2} + \frac{mc^2}{2} = mc^2$. Result matches paper. $\blacksquare$

**Mathematica check:**
```mathematica
ClearAll[p, V, m, cc];
H = Sqrt[cc^2 p^2 + (m cc^2 + V)^2];
K = H^2/(2 m cc^2) + m cc^2/2;
paperK = p^2/(2 m) + V + V^2/(2 m cc^2) + m cc^2;
FullSimplify[K - paperK]
(* Result: 0  ✅ (Wolfram MCP, 2026-05-11) *)
```

**Verdict:** ✅ Confirmed.

### Eq. (2.11) — Force with critical point at $r_0$ ⭐ resolves Maxwell Open Q3

**As printed:**
$$\mathbf{u} = \frac{\partial K}{\partial \mathbf{p}} = \frac{\mathbf{p}}{m}, \qquad \mathbf{F}_K = -\frac{\partial K}{\partial \mathbf{x}} = -\nabla V\!\left(1 + \frac{V}{mc^2}\right).$$

With $V = -mc^2 r_0/r$ (Coulomb-like), so $V/(mc^2) = -r_0/r$:
- **(1) $r > r_0$:** $1 + V/(mc^2) = 1 - r_0/r > 0$, and $-\nabla V = -mc^2 r_0/r^2 \cdot \hat{\mathbf{r}}$ pointing toward origin. ⟹ $\mathbf{F}_K$ **attractive**.
- **(2) $r = r_0$:** $1 - r_0/r = 0$. ⟹ $\mathbf{F}_K = 0$. **Critical point.**
- **(3) $r < r_0$:** $1 + V/(mc^2) < 0$, sign flip. ⟹ $\mathbf{F}_K$ **repulsive**.

**Mathematica check:**
```mathematica
ClearAll[r, r0, m, cc];
V[r_] := -m cc^2 r0/r;
Fnabla = -D[V[r], r] (1 + V[r]/(m cc^2));
(* r-hat component of F_K *)
Print["At r = r0/2 (should be > 0, repulsive): ", Simplify[Fnabla /. r -> r0/2, Assumptions -> {m > 0, cc > 0, r0 > 0}]];
Print["At r = r0   (should be = 0): ",          Simplify[Fnabla /. r -> r0,   Assumptions -> {m > 0, cc > 0, r0 > 0}]];
Print["At r = 2r0  (should be < 0, attractive): ", Simplify[Fnabla /. r -> 2 r0, Assumptions -> {m > 0, cc > 0, r0 > 0}]];
(* Results:
   At r = r0/2: +4 cc^2 m / r0    ✅ repulsive
   At r = r0:   0                  ✅ critical
   At r = 2 r0: -cc^2 m / (8 r0)   ✅ attractive
*)
```

**Verdict:** ✅ **All three sign cases confirmed. This is the explicit derivation of the "$r_0$ is a critical point" claim from Maxwell paper line 256.**

**Cross-reference to Maxwell paper Open Question 3:** The Maxwell paper made this critical-point claim *without* the derivation. This paper supplies it. With $V = -e^2/r$ and identifying $r_0 = e^2/(mc^2)$ (classical electron radius), $V = -mc^2$ at $r = r_0$, which is exactly the condition $1 + V/(mc^2) = 0$. The Maxwell paper's claim is **algebraically correct in the framework where the Hamiltonian uses "$V$ as part of the mass" form $H = \sqrt{c^2\mathbf{p}^2 + (mc^2 + V)^2}$**, *not* in the standard "$V$ as additive" form $H = \sqrt{c^2\mathbf{p}^2 + m^2c^4} + V$.

This resolves the units/sign concern in Open Question 3.

---

## Section 2.2.2 — Second case ("V as additive")

**As printed:** $H = \sqrt{c^2\mathbf{p}^2 + m^2c^4} + V$ ⟹ $K = \boldsymbol\pi^2/(2m) + mc^2 + V^2/(2mc^2) + V H_0/(mc^2)$.

**Status:** Identical to DRQM I Eq. (I.7) and Maxwell paper Eq. (16) expansion. The Hamilton equations (2.12, 2.13) reproduce DRQM I (I.7) exactly.

**Verdict:** ✅ Already verified.

---

## Section 3 — Einstein and dual Maxwell theory

This section re-derives results from the Maxwell paper and DRQM I in the language of dual classical theory. Cross-references:
- Sec 3.x: dual proper-time framework — Maxwell paper Eqs (1)–(2), (3)–(3'), DRQM I Eqs (I.1)–(I.9).
- Sec 3.6: dual Maxwell theory — Maxwell paper Eqs (3)–(11), DRQM I Eqs (I.8)–(I.11).

**Verdict:** ✅ All already verified in their respective primary docs.

---

## Sections 4–5 — Light/mass and historical cosmology

**Content:** Sec 4 discusses the photon mass (cross-references Maxwell paper Eq. 6). Sec 5 presents the speculative "symmetric big bang" / CMBR ideas in the dual-theory framework.

**Verdict:** ⬜ Conceptual / no new verifiable equations. Sec 4 photon mass formula matches Maxwell paper Eq. (6) (already verified).

---

## Summary

This paper has **important novel content**: Sec 2.2 supplies the explicit derivation of the "$r_0$ critical-point" claim that was Open Question 3 in the Maxwell paper. The mathematical result is correct as derived.

**Open Question 3 from the Maxwell paper is now RESOLVED:**
> The Maxwell paper claim that "$r_0$ is a critical point of the force $-\nabla\Phi - \nabla\Phi\,V/(mc^2)$" requires the Hamiltonian framework with $H = \sqrt{c^2\mathbf{p}^2 + (mc^2 + V)^2}$ (V as part of the mass), not the standard "V as additive" framework. Within that framework, the algebra is correct: $V = -e^2/r$, $r_0 = e^2/(mc^2)$, and $\mathbf{F}_K = -\nabla V(1 + V/(mc^2))$ has a critical point at $V = -mc^2$, which occurs at $r = r_0$. Wolfram MCP confirms attractive/zero/repulsive behavior for $r > r_0$, $r = r_0$, $r < r_0$.

**No new findings beyond the three already documented.** I'll update the Maxwell paper Open Questions section to note that #3 is resolved.

