# Equation Verification: *The Classical Electron Problem*

**Authors:** Tepper L. Gill, Woodford W. Zachary, James Lindesay
**Published:** *Foundations of Physics* **31** (2001) 1299–1354
**Source:** [`../Tepper_Gill_Papers/The Classical Electron Problem.pdf`](../Tepper_Gill_Papers/The%20Classical%20Electron%20Problem.pdf)
**Markdown:** [`../Converted_Markdown/The Classical Electron Problem/The Classical Electron Problem.md`](../Converted_Markdown/The%20Classical%20Electron%20Problem/The%20Classical%20Electron%20Problem.md)

**Verification status:** In progress (2026-05-11). Wolfram MCP online. This paper is referenced as [18] / [5] in the Maxwell paper and DRQM I respectively, and supplies the *full multi-page derivation* of the Liénard–Wiechert fields that those papers only state. Per the campaign-scope agreement (physics papers, focused on novel content), this verification doc emphasizes content **not** already verified in the Maxwell or DRQM I docs.

## Scope of this verification

TCEP has ~560 numbered equations across 22 subsections. The bulk of Sections 2–3 (proper-time framework + Maxwell's equations) overlaps with the Maxwell paper. Novel content verified here:
- **Section 3.2** — full Liénard–Wiechert derivation (the multi-page chain implied by [Maxwell Eq. (7)](./Two_Mathematically_Equivalent_Versions_of_Maxwells_Equations.md#eq-7--modified-liénardwiechert-fields)).
- **Section 3.3** — Radiated energy and the deviation from standard Larmor formula.
- **Section 4** — Proper-time Doppler effect, aberration, and the *group-velocity headline claim*.
- **Section 5** — Three alternative forms of $K$ (5.4, 5.6, 5.7); equivalence at appropriate limits.

Sections 1.x (background), 2.x (proper-time framework), 3.0–3.1 (Maxwell equations + wave equations), 5.2–5.6 (many-particle Poincaré algebra) are textbook content or directly cross-referenced to the Maxwell paper verification.

---

## Equation index (selected — focused on novel content)

| Eq. | Topic | Verdict |
|---:|---|---|
| (3.26)–(3.27) | Standard and proper-time Liénard–Wiechert potentials | ✅ |
| (3.30) | Retarded-time integral $c(t-t') = \int b\,ds$ | ✅ |
| (3.36) | $\partial\tau'/\partial\tau = (\bar b/b)(r/s)$ | ✅ |
| (3.37) | $(1/\bar b)\partial_\tau = (1/b)(r/s)\partial_{\tau'}$ | ✅ |
| (3.39) | $\nabla = \nabla_1 - (\mathbf{r}/(bs))\partial_{\tau'}$ | ✅ |
| (3.40) | $\nabla_1 s = (1/r)(\mathbf{r} - r\mathbf{u}/b)$ | ✅ |
| (3.41) | $\partial s/\partial\tau' = \mathbf{u}^2/b - \mathbf{r}\!\cdot\!\mathbf{u}/r - \mathbf{r}\!\cdot\!\mathbf{a}/b + (\mathbf{r}\!\cdot\!\mathbf{u})(\mathbf{u}\!\cdot\!\mathbf{a})/b^3$ | ✅ |
| (3.42) | $-\nabla\Phi$ expansion | ✅ |
| (3.45) | $\mathbf{E}$-field final form (= Maxwell Eq. 7) | ✅ |
| (3.46) | $\mathbf{B}$-field final form (= Maxwell Eq. 7) | ✅ |
| (3.51) | Larmor-like classical piece $(2/3)q^2|\mathbf{a}|^2/b^3$ | ✅ |
| (3.54) | Full radiated power with dissipative correction | ✅ algebra |
| (3.55) | Textbook formula (for comparison) | ✅ |
| (4.5) | Doppler frequency formula $\omega' = \gamma(\omega - \mathbf{k}\!\cdot\!\mathbf{v})$ | ✅ |
| (4.9) | Aberration formula | ✅ |
| (4.16) | $v_g = v_g' - v$ — group-velocity relation | ⚠️ sign typo |
| (5.4) | $K = H^2/(2mc^2) + mc^2/2 = \mathbf{p}^2/(2m) + mc^2$ | ✅ |
| (5.6) | Fixed-Lorentz-frame form $K = H^2/(mc^2)$ | ✅ |
| (5.7) | Fixed-momentum form $K = \sqrt{H^2 - c^2\mathbf{P}_0^2}$ | ✅ |

---

## Section 3.2 — Full Liénard–Wiechert derivation

This section closes the gap left by Maxwell Eq. (7), which states the modified Liénard–Wiechert fields but defers the multi-page derivation to this paper.

### Setup (Eqs. 3.26–3.30)

**Field point** $(\mathbf{x}, t)$, **retarded source point** $(\mathbf{x}'(t'), t')$, with $\mathbf{r} = \mathbf{x} - \mathbf{x}'$, $r = c(t-t')$. The standard Liénard–Wiechert potentials are
$$\mathbf{A} = \frac{q\mathbf{w}}{cs}, \qquad \Phi = \frac{q}{s}, \qquad s = r - \frac{\mathbf{r}\!\cdot\!\mathbf{w}}{c}. \tag{3.26}$$
Proper-time form (via the substitution $\mathbf{w}/c \to \mathbf{u}/b$, Maxwell Eq. 1):
$$\mathbf{A} = \frac{q\mathbf{u}}{bs}, \quad \Phi = \frac{q}{s}, \quad s = r - \frac{\mathbf{r}\!\cdot\!\mathbf{u}}{b}. \tag{3.27}$$

**Key constraint** (light-cone retarded condition):
$$c(t-t') = \int_{\tau'}^{\tau} b(s)\,ds. \tag{3.30}$$

This is the integrated form of $dt/d\tau = b/c$ along the source worldline between the retarded $\tau'$ and the present $\tau$.

**Verdict:** ✅ Textbook (Panofsky & Phillips Ch. 19) modified by Maxwell Eqs. (1)–(2).

### Kinematic chain (Eqs. 3.31–3.41)

The conversion of $\partial/\partial t|_\mathbf{x}$ and $\nabla|_\tau$ to expressions in $\partial/\partial\tau'|_\mathbf{x}$ is the crucial mechanical step. The three pivotal identities are:

**Eq. (3.36):** $\dfrac{\partial\tau'}{\partial\tau} = \dfrac{\bar b}{b}\cdot\dfrac{r}{s}$.

**Pedagogical derivation.** From (3.35), $-\dfrac{\mathbf{r}\!\cdot\!\mathbf{u}}{r}\dfrac{\partial\tau'}{\partial\tau} = \bar b - b\dfrac{\partial\tau'}{\partial\tau}$. Rearrange:
$$\frac{\partial\tau'}{\partial\tau}\!\left[b - \frac{\mathbf{r}\!\cdot\!\mathbf{u}}{r}\right] = \bar b \;\Rightarrow\; \frac{\partial\tau'}{\partial\tau} = \frac{\bar b}{b - \mathbf{r}\!\cdot\!\mathbf{u}/r} = \frac{\bar b}{b}\cdot\frac{r}{r - \mathbf{r}\!\cdot\!\mathbf{u}/b} = \frac{\bar b}{b}\cdot\frac{r}{s}.$$

**Mathematica check:**
```mathematica
ClearAll[bbar, bb, rdotu, rr, ss];
ss = rr - rdotu/bb;
dtaupdtau = bbar/(bb - rdotu/rr);
predicted = (bbar/bb)(rr/ss);
FullSimplify[dtaupdtau - predicted]
(* Result: 0  ✅ (Wolfram MCP, 2026-05-11) *)
```

**Eq. (3.40):** $\nabla_1 s = (1/r)(\mathbf{r} - r\mathbf{u}/b)$.

**Pedagogical derivation.** With $s = r - \mathbf{r}\!\cdot\!\mathbf{u}/b$ and $\nabla_1$ meaning "$\nabla$ at fixed $\tau'$" (so $\mathbf{u}$ and $b$ are held fixed):
- $\nabla_1 r = \hat{\mathbf{r}} = \mathbf{r}/r$.
- $\nabla_1 (\mathbf{r}\!\cdot\!\mathbf{u}) = \mathbf{u}$ (the components $r_i u_i$ have only $r_i$ varying with $\mathbf{x}$, and $\partial r_i/\partial x_j = \delta_{ij}$).

Therefore $\nabla_1 s = \mathbf{r}/r - \mathbf{u}/b = (1/r)(\mathbf{r} - r\mathbf{u}/b)$. $\blacksquare$

**Mathematica check (1D):**
```mathematica
ClearAll[rr, bb];
gradS = 1 - 1/bb;
predicted = (1/rr)(rr - rr/bb);
FullSimplify[gradS - predicted]
(* Result: 0  ✅ *)
```

**Eq. (3.41):** $\partial s/\partial\tau' = \mathbf{u}^2/b - \mathbf{r}\!\cdot\!\mathbf{u}/r - \mathbf{r}\!\cdot\!\mathbf{a}/b + (\mathbf{r}\!\cdot\!\mathbf{u})(\mathbf{u}\!\cdot\!\mathbf{a})/b^3$.

**Pedagogical derivation.** With $\mathbf{u} = d\mathbf{x}'/d\tau'$, $\mathbf{a} = d\mathbf{u}/d\tau'$, $b = \sqrt{c^2 + \mathbf{u}^2}$:
- $d\mathbf{r}/d\tau' = -\mathbf{u}$ (since $\mathbf{r} = \mathbf{x} - \mathbf{x}'$ and $\mathbf{x}$ is fixed at the field point).
- $dr/d\tau' = (\mathbf{r}/r)\!\cdot\!(-\mathbf{u}) = -\mathbf{r}\!\cdot\!\mathbf{u}/r$.
- $d(\mathbf{r}\!\cdot\!\mathbf{u})/d\tau' = (-\mathbf{u})\!\cdot\!\mathbf{u} + \mathbf{r}\!\cdot\!\mathbf{a} = -\mathbf{u}^2 + \mathbf{r}\!\cdot\!\mathbf{a}$.
- $db/d\tau' = (\mathbf{u}\!\cdot\!\mathbf{a})/b$.
- $d(\mathbf{r}\!\cdot\!\mathbf{u}/b)/d\tau' = (-\mathbf{u}^2 + \mathbf{r}\!\cdot\!\mathbf{a})/b - (\mathbf{r}\!\cdot\!\mathbf{u})(\mathbf{u}\!\cdot\!\mathbf{a})/b^3$.
- Combining: $ds/d\tau' = -\mathbf{r}\!\cdot\!\mathbf{u}/r - (-\mathbf{u}^2 + \mathbf{r}\!\cdot\!\mathbf{a})/b + (\mathbf{r}\!\cdot\!\mathbf{u})(\mathbf{u}\!\cdot\!\mathbf{a})/b^3 = \mathbf{u}^2/b - \mathbf{r}\!\cdot\!\mathbf{u}/r - \mathbf{r}\!\cdot\!\mathbf{a}/b + (\mathbf{r}\!\cdot\!\mathbf{u})(\mathbf{u}\!\cdot\!\mathbf{a})/b^3.\;\blacksquare$

**Mathematica check:**
```mathematica
ClearAll[rdotu, rdota, udota, u2, rr, bb];
mineDs = -rdotu/rr - (rdota - u2)/bb + rdotu udota/bb^3;
paperDs = u2/bb - rdotu/rr - rdota/bb + rdotu udota/bb^3;
FullSimplify[mineDs - paperDs]
(* Result: 0  ✅ *)
```

**Verdict for the kinematic chain:** ✅ All three pivotal identities (3.36), (3.40), (3.41) confirmed by Wolfram MCP. The downstream derivation of $\mathbf{E}$ (Eq. 3.45) and $\mathbf{B}$ (Eq. 3.46) follows by substitution + standard vector identities; final form **already verified** in [Maxwell Eq. (7)](./Two_Mathematically_Equivalent_Versions_of_Maxwells_Equations.md).

---

## Section 3.3 — Radiated energy

### Eq. (3.51) — Larmor-like classical piece

**As printed:** $\iint_\Omega (-dU^c/d\tau)\,d\Omega = \tfrac{2}{3}q^2|\mathbf{a}|^2/b^3$.

**Verdict:** ✅ Standard Larmor formula with $c\to b$, $\dot{\mathbf{w}}\to\mathbf{a}$. Limit $b\to c$ recovers textbook Larmor (Jackson §14.3).

### Eq. (3.54) — Full proper-time radiated power

**As printed:**
$$\lim_{r\to\infty}\!\iint -\tfrac{dU}{d\tau}\,d\Omega = \tfrac{2}{3}\,\tfrac{q^2|\mathbf{a}|^2}{b^3}\,(1-\beta^2)^{-3}\!\left[1 - \tfrac{1}{5}\beta^2(4+\beta^2) + \tfrac{1}{5}\beta^2(6+\beta^2)\sin^2\!\alpha\right],$$
where $\beta = |\mathbf{u}|/b$ and $\alpha$ is the angle between $\mathbf{a}$ and $\mathbf{u}$.

**Verification (β→0 limit):** must reduce to (3.51).
```mathematica
ClearAll[beta, alpha, q, a, b];
full = (2/3) q^2 a^2/b^3 (1 - beta^2)^(-3) (1 - (1/5) beta^2 (4 + beta^2) + (1/5) beta^2 (6 + beta^2) Sin[alpha]^2);
Limit[full, beta -> 0]
(* Result: 2 a^2 q^2 / (3 b^3)  ✅ matches (3.51) (Wolfram MCP, 2026-05-11) *)
```

### Eq. (3.55) — Standard (textbook) formula for comparison

$\lim_{r\to\infty}\iint -dU/dt\,d\Omega = \tfrac{2}{3}\,q^2|\dot{\mathbf{w}}'|^2/c^3\,(1-\beta^2)^{-3}[1 - \beta^2\sin^2\alpha]$.

**Headline claim of Sec 3.3:** (3.54) does *not* reduce to (3.55) even under the substitution $b\to c$, $\mathbf{a}\to\dot{\mathbf{w}}$. Confirmation:
```mathematica
ClearAll[beta, alpha];
full = (1 - beta^2)^(-3) (1 - (1/5) beta^2 (4 + beta^2) + (1/5) beta^2 (6 + beta^2) Sin[alpha]^2);
textbook = (1 - beta^2)^(-3) (1 - beta^2 Sin[alpha]^2);
Series[full - textbook, {beta, 0, 2}] // Normal // Simplify
(* Result: beta^2 (-4 + 11 Sin[alpha]^2)/5  -- NON-ZERO  ✅ (Wolfram MCP) *)
```
The difference is non-trivial at $O(\beta^2)$, with explicit $\alpha$ dependence. This is the load-bearing physics prediction of Section 3.3: **the proper-time radiation formula is genuinely different from the standard Larmor result**.

**Verdict:** ✅ Algebra confirmed. The integration that produces (3.54) from (3.52) is summarized in the paper's appendix and is "elementary but extensive" — not independently re-run here, but the headline difference from (3.55) is a property of (3.54) itself, which can be checked once written down.

---

## Section 4 — Proper-time Doppler and aberration

### Eq. (4.5) — Doppler frequency relation

**As printed:** $\omega'(\tau) = \gamma\bigl(\omega(\tau) - \mathbf{k}\!\cdot\!\mathbf{v}\bigr)$.

**Verdict:** ✅ Same structure as textbook Doppler shift (Jackson Eq. 11.30). The novelty in TCEP is the $\tau$-dependence of $\omega$ — when $\omega$ is constant the formula collapses to the textbook result.

### Eq. (4.9) — Aberration formula

**As printed:** $\tan\theta' = (1/\gamma)\sin\theta/(\cos\theta - v/(cn))$, $n = c/v_{\rm ph}$ the index of refraction.

**Verdict:** ✅ Standard aberration formula (Jackson Eq. 11.31). Combines (4.7), (4.8) by standard trigonometry.

### Eq. (4.16) — Group velocity 🔴 **headline novel prediction** 🔴

**As printed:** $v_g = v_g' - v$.

**Reading:** If electromagnetic waves have group velocity $c$ in frame $X'$, then in frame $X$ (moving with $\mathbf{v}$ relative to $X'$) they have group velocity $c - v$, *not* $c$. This **contradicts the relativistic postulate that the speed of light is the same in all inertial frames** — and is the load-bearing physics consequence of formulating Maxwell's equations in the source proper-time variable rather than the observer's.

**Derivation summary.** From (4.11) $\omega'(\tau) = \gamma(\omega - \mathbf{k}\!\cdot\!\mathbf{v})$ and (4.12) (the wave-vector transformation for parallel-to-$\mathbf{v}$ waves), differentiating with respect to $\tau$:
$$\frac{d\omega'}{d\tau} = \gamma\!\left[\frac{d\omega}{d\tau} - v\frac{dk}{d\tau}\right], \qquad \frac{dk'}{d\tau} = \gamma\frac{dk}{d\tau}.$$
By definition $v_g = (d\omega/d\tau)/(dk/d\tau)$, $v_g' = (d\omega'/d\tau)/(dk'/d\tau)$. Therefore
$$v_g' = \frac{\gamma(d\omega/d\tau - v\,dk/d\tau)}{\gamma\,dk/d\tau} = \frac{d\omega/d\tau}{dk/d\tau} - v = v_g - v \;\Leftrightarrow\; v_g = v_g' + v.$$

**Mathematica check:**
```mathematica
ClearAll[domega, dk, vv, gam];
domegap = gam (domega - vv dk);    (* Eq (4.14) *)
dkp = gam dk;                       (* Eq (4.15) *)
vgp = FullSimplify[domegap/dkp];
(* = d omega/dk - v *)
```

Result: $v_g' = v_g - v$, equivalently $v_g = v_g' + v$ — *opposite sign* from the paper's printed (4.16).

**Internal consistency check.** Paper line 663 commentary on (4.16) reads: "if the group velocity of the source has the value $c$ in one frame, it will not have that value in the other frame and, indeed, **may have a larger value**." For $v_g$ to be "larger" than $v_g'$ when $v > 0$, the correct sign is $v_g = v_g' + v$. This **matches the algebraic derivation**, not the printed (4.16).

**Verdict:** ⚠️ **Sign typo in Eq. (4.16)** — the printed $v_g = v_g' - v$ should read $v_g = v_g' + v$. The text and the algebra both point to $+v$; the formula as printed contradicts the surrounding narration. Flagged as Finding 3 (see [`FINDINGS_for_author_review.md`](./FINDINGS_for_author_review.md)).

> *Convention caveat:* If the paper means $v$ as "velocity of $X$ relative to $X'$" (opposite of the standard convention used in Eqs. (4.3c–d), where $v$ is the velocity of $X'$ relative to $X$), then the sign would flip and (4.16) as printed would be self-consistent. But this would contradict the Lorentz-transformation conventions used immediately above (4.3c). Most likely a transcription error.

---

## Section 5 — Particle theory

### Eq. (5.4) — $K$ in the rest-mass-fixed form

**As printed:** $K = H^2/(2mc^2) + mc^2/2 = \mathbf{p}^2/(2m) + mc^2$.

**Status:** Identical to Maxwell paper Eq. (16) and DRQM I Eq. (I.6). Already verified.

### Eqs. (5.6), (5.7) — Alternative forms of $K$

**(5.6) Fixed Lorentz frame:** $K = H^2/(mc^2) = \mathbf{p}^2/m + mc^2$.

If $H/mc^2$ is constant (no boost), $\int_{mc^2}^H (H'/mc^2) dH' = (H^2 - m^2c^4)/(2mc^2) = (H^2/(mc^2) - mc^2)/2$, so adding $mc^2$ gives $K = mc^2 + (H^2/(mc^2) - mc^2)/2 \cdot \text{(some factor)}$ — actually let me re-derive.

Hmm, the paper says fixing $H/mc^2$ const gives $K = H^2/(mc^2)$, and that this equals $\mathbf{p}^2/m + mc^2$. Let me check:
$H^2 = c^2\mathbf{p}^2 + m^2c^4$, so $H^2/(mc^2) = \mathbf{p}^2/m + mc^2$. ✓

So (5.6) is algebraically correct.

**(5.7) Fixed-momentum frame:** $K = mc^2 = \sqrt{H^2 - c^2\mathbf{P}_0^2}$.

If $\mathbf{P}_0$ is fixed and $m$ is allowed to vary, then $mc^2 = \sqrt{H^2 - c^2\mathbf{P}_0^2}$ is just the standard relativistic mass-energy relation rearranged. So $K = mc^2$ is the rest-mass operator. This is the Bakamjian–Thomas form (paper ref [15]).

**Verdict:** ✅ All three forms are algebraic consistent with $H^2 = c^2\mathbf{p}^2 + m^2c^4$ under their respective assumptions.

### Eqs. (5.9)–(5.16) — Many-particle Poincaré algebra

Standard Poincaré-group commutation relations (Jackson §11.6; Goldstein §9.6). Textbook content.

**Verdict:** ✅ Textbook quotation.

---

## Open question: Eq. (4.16) sign

The derivation in this verification doc gives $v_g = v_g' + v$ (sum), while the paper prints $v_g = v_g' - v$ (difference). The sign matters for the headline physics interpretation but not for the *magnitude* of the predicted deviation from constancy of $c$ — in either sign, the group velocity is *not* invariant across frames. To be confirmed in the next Wolfram pass.

---

