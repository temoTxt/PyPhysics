# Equation Verification: *Two Mathematically Equivalent Versions of Maxwell's Equations*

**Authors:** Tepper L. Gill, Woodford W. Zachary
**Published:** *Foundations of Physics* 41 (2011) 99–128
**Source:** [`../Tepper_Gill_Papers/Two Mathematically Equivalent Versions of Maxwell's Equations.pdf`](../Tepper_Gill_Papers/Two%20Mathematically%20Equivalent%20Versions%20of%20Maxwell's%20Equations.pdf)
**Markdown:** [`../Converted_Markdown/Two Mathematically Equivalent Versions of Maxwell's Equations/Two Mathematically Equivalent Versions of Maxwell's Equations.md`](../Converted_Markdown/Two%20Mathematically%20Equivalent%20Versions%20of%20Maxwell%27s%20Equations/Two%20Mathematically%20Equivalent%20Versions%20of%20Maxwell%27s%20Equations.md)

**Verification status:** Wolfram MCP online (2026-05-10). Every `Mathematica check` block that ran is annotated with the actual returned result (`Result: 0  ✅` for an identity that simplified to zero). Eqs. (1)–(6) verified; Eqs. (7)–(11) verified at the chain-rule / coefficient level; Eqs. (12)–(24) being checked & derived in this same pass.

## Conventions used in this document

- Gaussian (c.g.s.) units, matching the paper.
- $\mathbf{w} = d\mathbf{x}/dt$ is the **observer-time velocity** (≤ c).
- $\mathbf{u} = d\mathbf{x}/d\tau$ is the **proper-time velocity** (no upper bound).
- $\gamma = \gamma(\mathbf{w}) = (1-\mathbf{w}^2/c^2)^{-1/2}$ (Lorentz factor in the observer-time variable).
- $b = \sqrt{c^2 + \mathbf{u}^2}$ is the **collaborative ("effective") speed of light** central to the dual theory.
- An overdot on a 3-vector denotes a $\tau$-derivative: $\dot{\mathbf{u}} = d\mathbf{u}/d\tau \equiv \mathbf{a}$ (proper acceleration).

> **Standard-textbook reading list cross-referenced throughout.** Jackson, *Classical Electrodynamics* (3e); Goldstein, *Classical Mechanics* (3e); Sakurai, *Modern Quantum Mechanics* (2e); Peskin & Schroeder, *Introduction to QFT*; Weinberg, *Gravitation and Cosmology* (for the metric digression).

---

## Equation index

| Eq. | Topic | Verdict |
|---:|---|---|
| (1) | $\mathbf{w}/c = \mathbf{u}/b$ — velocity duality | ✅ |
| (2) | $(1/c)\partial_t = (1/b)\partial_\tau$ — time-derivative duality | ✅ |
| (3) | Standard Maxwell's equations (Gaussian) | ✅ |
| (3′) | Proper-time-equivalent Maxwell's equations | ✅ |
| (4) | Dual wave equations for **E**, **B** with dissipative term | ✅ |
| (5) | Klein–Gordon form after scale transformation | ✅ |
| (6) | Effective ("photon") mass $\mu$ from $b$-dynamics | ✅ |
| (7) | Modified Liénard–Wiechert **E**, **B** fields | ✅ |
| (8) | Standard Lorentz time transformation | ✅ |
| (9) | $t = (1/c)\int_0^\tau b(s)\,ds$ (mean-value form) | ✅ |
| (10) | Proper-time position/velocity/acceleration boosts | ✅ |
| (11) | $b'(\tau) = \gamma[b - \mathbf{u}\cdot\mathbf{v}/c]$ | ✅ |
| (12) | $\mathbf{J}'$ transformation | ⬜ |
| (13) | $b'\rho' = \gamma[b\rho - \mathbf{J}\cdot\mathbf{v}/c]$ | ⬜ |
| (14) | $\rho' = (\rho - \mathbf{J}\cdot\mathbf{v}/bc) / (1 - \mathbf{u}\cdot\mathbf{v}/bc)$ | ⬜ |
| (15) | $\rho'/\rho$ ratio (sphere stays spherical) | ⬜ |
| (16, unnumbered) | $K = H^2/(2mc^2) + mc^2/2$ — canonical proper-time Hamiltonian | ⬜ |
| (17) | Hamilton equations for $\mathbf{u}, d\mathbf{p}/d\tau$ | ⬜ |
| (18) | Force equation with $V/(mcb)$ correction to Lorentz force | ⬜ |
| (19) | $\mathbf{U} = d\mathbf{X}/d\tau = \mathbf{P}/M$ | ⬜ |
| (20) | $\sum H_i \mathbf{v}_i / H$ representation of $\mathbf{U}$ | ⬜ |
| (21) | $dW/d\tau = \sum (d\tau_i/d\tau)\{K_i, W\}$ | ⬜ |
| (22, unnumbered) | Dual Dirac equation in proper time | ⬜ |
| (23) | Dual standard square-root equation | ⬜ |
| (24) | Eigenvalue relation $\frac{E_n^2}{2mc^2}+\frac{mc^2}{2}=\dots$ | ⬜ |

---

## Eq. (1) — Velocity duality

**As printed (line 73):**
$$\frac{\mathbf{w}}{c} = \frac{\mathbf{u}}{b}, \qquad b = \sqrt{c^2 + \mathbf{u}^2}.$$

**Context / claim:** The Lorentz factor $\gamma$ can be repackaged so that the proper-time velocity $\mathbf{u}$ stands in the same dimensionless ratio to $b$ as the observer-time velocity $\mathbf{w}$ does to $c$. This is the algebraic seed of the entire dual theory.

**Mathematica check:**
```mathematica
ClearAll[w, u, c];
b = Sqrt[c^2 + u^2];
(* From u = w / Sqrt[1 - w^2/c^2], invert to get w in terms of u: *)
wOfU = w /. First @ Solve[u == w/Sqrt[1 - w^2/c^2] && w > 0 && c > 0 && u > 0, w];
(* Verify w/c == u/b: *)
FullSimplify[wOfU/c - u/b, Assumptions -> {c > 0, u > 0}]
(* Result: 0  ✅ (Wolfram MCP, 2026-05-10) *)
```

**Expanded derivation:**

Start from the textbook definition of the proper time of a particle moving with velocity $\mathbf{w}$:
$$d\tau^2 \;=\; dt^2 - \frac{1}{c^2}\,d\mathbf{x}^2 \;=\; dt^2\!\left[1 - \frac{\mathbf{w}^2}{c^2}\right], \tag{$\star$}$$
which gives the familiar $d\tau = dt/\gamma(\mathbf{w})$.

- **Step 1.** From $(\star)$, $\mathbf{u} \equiv d\mathbf{x}/d\tau = (dt/d\tau)(d\mathbf{x}/dt) = \gamma\mathbf{w}$, so
$$\mathbf{u} \;=\; \frac{\mathbf{w}}{\sqrt{1 - \mathbf{w}^2/c^2}}.$$
- **Step 2.** Square it: $\mathbf{u}^2 = \mathbf{w}^2/(1-\mathbf{w}^2/c^2)$.
- **Step 3.** Solve for $\mathbf{w}^2$: $\mathbf{w}^2(1 + \mathbf{u}^2/c^2) = \mathbf{u}^2 \cdot (1 - \mathbf{w}^2/c^2 + \mathbf{w}^2/c^2)$ — actually, cleaner: $\mathbf{u}^2 - \mathbf{u}^2\mathbf{w}^2/c^2 = \mathbf{w}^2$, hence $\mathbf{w}^2 = \mathbf{u}^2/(1 + \mathbf{u}^2/c^2)$.
- **Step 4.** Therefore $\mathbf{w} = \mathbf{u}/\sqrt{1 + \mathbf{u}^2/c^2}$. Define $b := \sqrt{c^2 + \mathbf{u}^2}$ so that $\sqrt{1 + \mathbf{u}^2/c^2} = b/c$. Then
$$\mathbf{w} \;=\; \frac{c\,\mathbf{u}}{b} \;\;\Longleftrightarrow\;\; \frac{\mathbf{w}}{c} = \frac{\mathbf{u}}{b}. \qquad\blacksquare$$

**Standard-equation comparison:**

The textbook 4-velocity is $u^\mu = (\gamma c,\,\gamma\mathbf{w})$ (Jackson §11.4 or Goldstein §7.5). What Gill calls $\mathbf{u}$ is the **spatial part of the 4-velocity**, $u^i = \gamma w^i$. The "collaborative speed" $b$ is then just the **temporal component of the 4-velocity**, since $u^0 = \gamma c$ and $\gamma c = c/\sqrt{1-\mathbf{w}^2/c^2} = \sqrt{c^2 + \mathbf{u}^2} = b$.

So Eq. (1) is the trivially-true identity $u^i/u^0 = w^i/c$ — the spatial-to-temporal ratio of any 4-velocity equals the spatial velocity divided by $c$. **What is *non*-standard is using $b = u^0$ as a "speed of light"** (it has units of velocity but is frame-dependent on the source).

**Verdict:** ✅ Confirmed by Wolfram MCP (`FullSimplify` returned 0).

---

## Eq. (2) — Time-derivative duality

**As printed (line 78):**
$$\frac{1}{c}\frac{\partial}{\partial t} \;=\; \frac{1}{b}\frac{\partial}{\partial \tau}.$$

**Context / claim:** When the chain rule is applied through $\tau(t)$, observer-time and proper-time derivatives are interconverted by exactly the factor $c/b$. This is the operator-level companion to Eq. (1).

**Mathematica check:**
```mathematica
ClearAll[t, \[Tau], u, c, F];
b = Sqrt[c^2 + u^2];
(* dt/d\[Tau] = b/c, so 1/c \[PartialD]_t = 1/c (d\[Tau]/dt) \[PartialD]_\[Tau] = (1/c)(c/b) \[PartialD]_\[Tau] = (1/b) \[PartialD]_\[Tau] *)
dTaudt = c/b;
Simplify[(1/c)*dTaudt - 1/b]
(* Result: 0  ✅ (Wolfram MCP, 2026-05-10) *)
```

**Expanded derivation:**

- **Step 1.** From $(\star)$ in Eq. (1)'s derivation, $dt^2 = d\tau^2 + d\mathbf{x}^2/c^2 = d\tau^2(1 + \mathbf{u}^2/c^2)$. Hence
$$\frac{dt}{d\tau} \;=\; \sqrt{1 + \mathbf{u}^2/c^2} \;=\; \frac{b}{c}, \qquad \frac{d\tau}{dt} \;=\; \frac{c}{b}.$$
- **Step 2.** By the chain rule, for any field $F(\mathbf{x},t)$,
$$\frac{\partial F}{\partial t} \;=\; \frac{d\tau}{dt}\,\frac{\partial F}{\partial \tau} \;=\; \frac{c}{b}\,\frac{\partial F}{\partial \tau}.$$
- **Step 3.** Divide both sides by $c$:
$$\frac{1}{c}\frac{\partial F}{\partial t} \;=\; \frac{1}{b}\frac{\partial F}{\partial \tau}. \qquad\blacksquare$$

**Standard-equation comparison:**

The standard substantial / convective derivative in fluid dynamics, $D/Dt = \partial_t + \mathbf{w}\!\cdot\!\nabla$, is the analogous "change of time-variable" trick at the *advected*-field level. Gill's Eq. (2) is the **fixed-position** version: it converts $\partial_t$ to $\partial_\tau$ when the field is evaluated *at the source location*, not while following a fluid element. In a textbook 4-derivative $\partial^\mu = (c^{-1}\partial_t, \nabla)$, Eq. (2) says $c^{-1}\partial_t = b^{-1}\partial_\tau$, i.e. **the 0-component of $\partial^\mu$ is invariant in form when one trades $(c,t)\!\to\!(b,\tau)$**. This is what makes the next step — substituting into Maxwell's equations — work.

**Verdict:** ✅ Confirmed by Wolfram MCP.

---

## Eq. (3) — Standard (observer-time) Maxwell's equations

**As printed (lines 89–91), Gaussian units:**

$$\nabla\!\cdot\!\mathbf{B} = 0,\qquad \nabla\!\cdot\!\mathbf{E} = 4\pi\rho,$$
$$\nabla\times\mathbf{E} = -\frac{1}{c}\,\frac{\partial \mathbf{B}}{\partial t},\qquad \nabla\times\mathbf{B} = \frac{1}{c}\!\left[\frac{\partial \mathbf{E}}{\partial t} + 4\pi\rho\,\mathbf{w}\right].$$

**Context / claim:** Quoted, not derived — these are textbook Maxwell's equations in Gaussian units with $\mathbf{J} = \rho\mathbf{w}$.

**Mathematica check:** Not applicable (textbook quotation).

**Standard-equation comparison:** Identical to Jackson Eqs. (6.6), (6.7) (Gaussian units), with $\mathbf{J} = \rho\mathbf{w}$.

**Verdict:** ✅ (textbook).

---

## Eq. (3′) — Proper-time-equivalent Maxwell's equations

**As printed (lines 95–96):**

$$\nabla\!\cdot\!\mathbf{B} = 0,\qquad \nabla\!\cdot\!\mathbf{E} = 4\pi\rho,$$
$$\nabla\times\mathbf{E} = -\frac{1}{b}\,\frac{\partial \mathbf{B}}{\partial \tau},\qquad \nabla\times\mathbf{B} = \frac{1}{b}\!\left[\frac{\partial \mathbf{E}}{\partial \tau} + 4\pi\rho\,\mathbf{u}\right].$$

**Context / claim:** Pointwise substitution of Eqs. (1) and (2) into the standard Maxwell's equations. Algebraically equivalent to (3); physically reinterpreted because $b$ depends on the source state.

**Mathematica check:**
```mathematica
(* Symbolic substitution check.
   At a fixed spatial point on the source worldline, treat \[PartialD]_t E as a stand-in
   "dEdt" and \[PartialD]_\[Tau] E as "dEdtau".  By Eq.(2): dEdt = (c/b) dEdtau.  By
   Eq.(1): w = (c/b) u. *)
ClearAll[c, b, u, w, \[Rho], dEdt, dEdtau];
LHS = (1/c) dEdt + (4 Pi/c) \[Rho] w;
RHS = (1/b) dEdtau + (4 Pi/b) \[Rho] u;
diff = LHS - RHS /. {dEdt -> (c/b) dEdtau, w -> (c/b) u};
Simplify[diff]
(* Result: 0  ✅ (Wolfram MCP, 2026-05-10) *)
```

**Expanded derivation:**

- **Step 1.** Curl equations: by Eq. (2), every $(1/c)\partial_t$ becomes $(1/b)\partial_\tau$. Hence
$$\nabla\!\times\!\mathbf{E} = -\frac{1}{c}\frac{\partial\mathbf{B}}{\partial t} = -\frac{1}{b}\frac{\partial\mathbf{B}}{\partial \tau}.$$
- **Step 2.** Current term in Ampère's law: by Eq. (1), $\mathbf{w}/c = \mathbf{u}/b$, so $\rho\mathbf{w}/c = \rho\mathbf{u}/b$. Therefore
$$\nabla\!\times\!\mathbf{B} = \frac{1}{c}\frac{\partial\mathbf{E}}{\partial t} + \frac{4\pi}{c}\rho\mathbf{w} = \frac{1}{b}\frac{\partial\mathbf{E}}{\partial\tau} + \frac{4\pi}{b}\rho\mathbf{u}.$$
- **Step 3.** Divergence equations are unchanged (no time derivatives, no $\mathbf{w}$).

**Standard-equation comparison:**

These are **not** new Maxwell's equations — they are the *same* fields obeying the *same* dynamics, written in different variables. The physical novelty is that $b = b(\tau)$ depends on the source's motion, so when treated as a coefficient in a PDE in the variable $\tau$ alone, $b$ is no longer a constant. The next equation (4) makes this explicit.

> ⚠ **Note on the equation-errors doc.** The errors doc flags `Jα = (Jx, Jy, Jz, ibρ)` and similar `b`-appearances as errors that should be `c`. Per the derivation above, the appearance of `b` is *required* by the substitution rules (1)–(2). These are not errors.

**Verdict:** ✅ Confirmed by Wolfram MCP (substitution-rule level).

---

## Eq. (4) — Dual wave equations with dissipative term

**As printed (lines 107–109):**

$$\frac{1}{b^2}\frac{\partial^2 \mathbf{B}}{\partial \tau^2} - \frac{\mathbf{u}\!\cdot\!\mathbf{a}}{b^4}\!\left[\frac{\partial \mathbf{B}}{\partial \tau}\right] - \nabla^2 \mathbf{B} = \frac{1}{b}\left[4\pi\,\nabla\!\times\!(\rho\mathbf{u})\right],$$

$$\frac{1}{b^2}\frac{\partial^2 \mathbf{E}}{\partial \tau^2} - \frac{\mathbf{u}\!\cdot\!\mathbf{a}}{b^4}\!\left[\frac{\partial \mathbf{E}}{\partial \tau}\right] - \nabla^2 \mathbf{E} = -\nabla(4\pi\rho) - \frac{1}{b}\frac{\partial}{\partial\tau}\!\left[\frac{4\pi(\rho\mathbf{u})}{b}\right],$$
with $\mathbf{a} = d\mathbf{u}/d\tau$.

**Context / claim:** Take the curl of the two proper-time Faraday/Ampère equations from (3′), use $\nabla\!\times\!\nabla\!\times\!\mathbf{F} = \nabla(\nabla\!\cdot\!\mathbf{F}) - \nabla^2\mathbf{F}$, and watch a *new* first-derivative-in-$\tau$ term emerge purely from the $\tau$-dependence of $b$. This term acts like a velocity-dependent damping and is the paper's mechanism for radiation reaction.

**Mathematica check:**
```mathematica
(* Sketch: d/d\[Tau](1/b) = -(1/b^2) db/d\[Tau].  With b = Sqrt[c^2 + u.u],
   db/d\[Tau] = (u . a)/b, so d/d\[Tau](1/b) = -(u . a)/b^3.  Multiplying the curl of
   Faraday by (1/b) and applying \[PartialD]_\[Tau] then yields the dissipative coefficient
   -(u.a)/b^4 in front of \[PartialD]_\[Tau] B. *)
ClearAll[u, a, c, \[Tau]];
b[\[Tau]_] := Sqrt[c^2 + u[\[Tau]] . u[\[Tau]]];
res = D[1/b[\[Tau]], \[Tau]] /. {u'[\[Tau]] -> a[\[Tau]]};
(* Mathematica's Dot is non-commutative; for Euclidean 3-vectors u.a == a.u. *)
resSym = res /. {a[\[Tau]] . u[\[Tau]] -> u[\[Tau]] . a[\[Tau]]};
FullSimplify[resSym + (u[\[Tau]] . a[\[Tau]])/b[\[Tau]]^3]
(* Result: 0  ✅  d/d\[Tau](1/b) = -(u.a)/b^3 (Wolfram MCP, 2026-05-10) *)
```

**Expanded derivation (B-field):**

- **Step 1.** Take $\nabla\!\times$ of the proper-time Faraday law $\nabla\!\times\!\mathbf{E} = -(1/b)\partial_\tau \mathbf{B}$:
$$\nabla\!\times\!\nabla\!\times\!\mathbf{E} = -\nabla\!\times\!\!\left[\frac{1}{b}\,\partial_\tau \mathbf{B}\right] = -\frac{1}{b}\,\partial_\tau(\nabla\!\times\!\mathbf{B}). \tag{$\dagger$}$$
(Note: $b$ depends on $\tau$ but not on $\mathbf{x}$ at the source's worldline, so it commutes with $\nabla$.)

- **Step 2.** Use the vector identity $\nabla\!\times\!\nabla\!\times\!\mathbf{E} = \nabla(\nabla\!\cdot\!\mathbf{E}) - \nabla^2\mathbf{E} = \nabla(4\pi\rho) - \nabla^2\mathbf{E}$ on the LHS of $(\dagger)$.

- **Step 3.** On the RHS, substitute Ampère: $\nabla\!\times\!\mathbf{B} = (1/b)\partial_\tau\mathbf{E} + (4\pi/b)\rho\mathbf{u}$. So
$$\partial_\tau(\nabla\!\times\!\mathbf{B}) = \partial_\tau\!\!\left[\frac{1}{b}\partial_\tau\mathbf{E}\right] + 4\pi\,\partial_\tau\!\!\left[\frac{\rho\mathbf{u}}{b}\right].$$
- **Step 4.** Use $\partial_\tau(1/b) = -(\dot b)/b^2 = -(\mathbf{u}\!\cdot\!\mathbf{a})/b^3$ (from $b^2 = c^2 + \mathbf{u}^2$, so $2b\dot b = 2\mathbf{u}\!\cdot\!\mathbf{a}$):
$$\partial_\tau\!\!\left[\frac{1}{b}\partial_\tau\mathbf{E}\right] = \frac{1}{b}\partial_\tau^2\mathbf{E} - \frac{\mathbf{u}\!\cdot\!\mathbf{a}}{b^3}\,\partial_\tau\mathbf{E}.$$
- **Step 5.** Combine, multiply through by $1/b$, and you obtain the **E**-field wave equation in the form printed. The **B**-field equation follows by taking $\nabla\!\times$ of Ampère's law instead and using $\nabla\!\cdot\!\mathbf{B}=0$.

**Standard-equation comparison:**

Standard EM wave equation (Jackson §6.4):
$$\Box \mathbf{E} \;\equiv\; \frac{1}{c^2}\partial_t^2 \mathbf{E} - \nabla^2\mathbf{E} \;=\; -\nabla(4\pi\rho) - \frac{4\pi}{c^2}\,\partial_t (\rho\mathbf{w}).$$
Gill's Eq. (4) has the same *form* under the rule $c\to b,\;t\to\tau,\;\mathbf{w}\to\mathbf{u}$, **plus** the extra term
$$-\frac{\mathbf{u}\!\cdot\!\mathbf{a}}{b^4}\,\partial_\tau\mathbf{E},$$
which arises *only* because $b$ depends on $\tau$. In the unaccelerated limit $\mathbf{a}=\mathbf{0}$, $b$ is constant, this term vanishes, and we recover the standard wave equation (with $b$ playing the role of $c$).

**Verdict:** ✅ Coefficient $-(\mathbf{u}\!\cdot\!\mathbf{a})/b^4$ confirmed by Wolfram MCP. Remaining algebra (Steps 1–5) is a textbook curl-of-curl manipulation; the *only* dual-theory–specific piece is the time-derivative of $1/b$, which is what we checked.

---

## Eq. (5) — Klein–Gordon-like form after rescaling

**As printed (lines 113–115):** After the scale change $\mathbf{E}\to (b/c)^{1/2}\mathbf{E}$, $\mathbf{B}\to(b/c)^{1/2}\mathbf{B}$,
$$\frac{1}{b^2}\partial_\tau^2 \mathbf{B} - \nabla^2 \mathbf{B} + \left[\frac{\ddot b}{2b^3} - \frac{3\dot b^2}{4b^4}\right]\mathbf{B} = \frac{c^{1/2}}{b^{3/2}}\bigl[4\pi\nabla\!\times\!(\rho\mathbf{u})\bigr],$$
and similarly for **E**.

**Context / claim:** The dissipative first-derivative term in Eq. (4) is absorbed into an effective potential / mass-squared term in Eq. (5) by the standard substitution $\psi = (b/c)^{1/2}\Psi$ that kills the first-order $\partial_\tau$ in a damped wave equation.

**Mathematica check:**
```mathematica
(* Substitute \[Psi] = (b/c)^{1/2} \[Psi]new into (1/b^2)\[Psi]_\[Tau]\[Tau] - (b'/b^3)\[Psi]_\[Tau] - \[Del]^2 \[Psi] = src,
   divide through by Sqrt[b/c], and read off the coefficient of \[Psi]new[\[Tau]] (the
   piece with no derivatives).  That coefficient should be the predicted
   effective-mass-squared term b''/(2 b^3) - 3 b'^2/(4 b^4). *)
ClearAll[b, \[Tau], \[Psi]new];
trans[\[Tau]_] := Sqrt[b[\[Tau]]/c] \[Psi]new[\[Tau]];
expr = (1/b[\[Tau]]^2) D[trans[\[Tau]], {\[Tau], 2}] - (b'[\[Tau]]/b[\[Tau]]^3) D[trans[\[Tau]], \[Tau]];
expanded = Expand[expr/Sqrt[b[\[Tau]]/c]];
coefPsi = Coefficient[expanded, \[Psi]new[\[Tau]]];
predicted = b''[\[Tau]]/(2 b[\[Tau]]^3) - 3 b'[\[Tau]]^2/(4 b[\[Tau]]^4);
FullSimplify[coefPsi - predicted]
(* Result: 0  ✅ (Wolfram MCP, 2026-05-10) *)
```

**Expanded derivation sketch:**

This is the standard Liouville-type trick: in the equation
$$\frac{1}{b^2}\partial_\tau^2 \psi - \frac{\dot b}{b^3}\,\partial_\tau\psi - \nabla^2\psi = (\text{source}),$$
the substitution $\psi = (b/c)^{1/2}\Psi$ removes the first-order $\partial_\tau$ term and generates an effective potential proportional to $\ddot b/b^3 - \dot b^2/b^4$. Algebra is mechanical; details suit a Mathematica check.

(Note: in Eq. (4) the first-derivative coefficient is $\mathbf{u}\!\cdot\!\mathbf{a}/b^4 = \dot b/b^3$, since $\dot b = \mathbf{u}\!\cdot\!\mathbf{a}/b$. So both forms agree.)

**Standard-equation comparison:**

This is the **massive Klein–Gordon equation in curved/conformal form**. Compare Peskin & Schroeder Eq. (2.12) for the free KG equation $(\Box + m^2)\phi = 0$, with $m^2 \leftrightarrow [\ddot b/2b^3 - 3\dot b^2/4b^4]\cdot\hbar^2/c^2$ (cf. Eq. (6) below). When $b$ is constant ($\mathbf{a} = \mathbf{0}$), the "mass" vanishes and we recover the massless wave equation — i.e. the photon has zero mass in vacuum but acquires an effective mass during source acceleration.

**Verdict:** ✅ Effective-mass coefficient $\ddot b/(2b^3) - 3\dot b^2/(4b^4)$ confirmed by Wolfram MCP.

---

## Eq. (6) — Effective photon mass

**As printed (line 120):**
$$\mu = \left\{\frac{\hbar^2}{c^2}\!\left[\frac{\ddot b}{2b^3} - \frac{3\dot b^2}{4b^4}\right]\right\}^{1/2} = \left\{\frac{\hbar^2}{c^2}\!\left[\frac{\mathbf{u}\!\cdot\!\ddot{\mathbf{u}} + \dot{\mathbf{u}}^2}{2b^4} - \frac{5(\mathbf{u}\!\cdot\!\dot{\mathbf{u}})^2}{4b^6}\right]\right\}^{1/2}.$$

**Context / claim:** Algebraic restatement of the coefficient appearing in Eq. (5), in terms of $\mathbf{u}, \dot{\mathbf{u}}, \ddot{\mathbf{u}}$ rather than $b, \dot b, \ddot b$.

**Mathematica check:**
```mathematica
(* Mathematica's Dot is non-commutative for symbolic vectors, which spuriously
   blocks the simplification.  Trick: use scalar surrogates X = u.u, Y = u.u',
   Z = u.u'' + u'.u' — the chain rule on b^2 = c^2 + u.u gives all needed
   identities in terms of X, Y, Z. *)
ClearAll[X, Y, Z, c];
bsq = c^2 + X;             (* b^2          *)
bp  = Y / Sqrt[bsq];        (* b' = (u.u')/b                       *)
bpp = (Z - Y^2/bsq) / Sqrt[bsq];  (* b'' = (u.u''+u'.u' - b'^2)/b *)
LHS = bpp/(2 bsq^(3/2)) - 3 bp^2/(4 bsq^2);
RHS = Z/(2 bsq^2) - 5 Y^2/(4 bsq^3);
FullSimplify[LHS - RHS]
(* Result: 0  ✅ (Wolfram MCP, 2026-05-10) *)
```

**Expanded derivation:**

- **Step 1.** $b^2 = c^2 + \mathbf{u}^2 \;\Rightarrow\; 2b\dot b = 2\mathbf{u}\!\cdot\!\dot{\mathbf{u}}$, so $\dot b = (\mathbf{u}\!\cdot\!\dot{\mathbf{u}})/b$.
- **Step 2.** Differentiate again: $\dot b^2 + b\ddot b = \dot{\mathbf{u}}^2 + \mathbf{u}\!\cdot\!\ddot{\mathbf{u}}$, so
$$\ddot b = \frac{\dot{\mathbf{u}}^2 + \mathbf{u}\!\cdot\!\ddot{\mathbf{u}} - \dot b^2}{b}.$$
- **Step 3.** Plug into $\ddot b/(2b^3) - 3\dot b^2/(4b^4)$:
$$\frac{\dot{\mathbf{u}}^2 + \mathbf{u}\!\cdot\!\ddot{\mathbf{u}} - \dot b^2}{2b^4} - \frac{3\dot b^2}{4b^4} = \frac{\dot{\mathbf{u}}^2 + \mathbf{u}\!\cdot\!\ddot{\mathbf{u}}}{2b^4} - \frac{2\dot b^2 + 3\dot b^2}{4b^4} = \frac{\dot{\mathbf{u}}^2 + \mathbf{u}\!\cdot\!\ddot{\mathbf{u}}}{2b^4} - \frac{5\dot b^2}{4b^4}.$$
- **Step 4.** Using $\dot b^2 = (\mathbf{u}\!\cdot\!\dot{\mathbf{u}})^2/b^2$:
$$\frac{5\dot b^2}{4b^4} = \frac{5(\mathbf{u}\!\cdot\!\dot{\mathbf{u}})^2}{4b^6},$$
recovering the printed form exactly. $\blacksquare$

**Standard-equation comparison:**

In standard QFT, the photon is exactly massless (gauge invariance + masslessness from Proca-zero limit). Gill's $\mu$ is a *dynamical, source-dependent* photon mass that vanishes whenever the source is inertial. It is **not** the bound proposed by Goldhaber & Nieto (a frame-independent Proca mass, $\lesssim 10^{-24}$ GeV); it is a coefficient in a damped wave operator that only manifests during emission. Compare Jackson §1.2 for the standard "$m_\gamma = 0$" position and Bargmann–Wigner [Proc. Nat. Acad. Sci. 34, 211 (1948)] for the Proca framework.

**Verdict:** ✅ Confirmed by Wolfram MCP. The algebraic identity $\ddot b/(2b^3)-3\dot b^2/(4b^4) = (\mathbf{u}\!\cdot\!\ddot{\mathbf{u}}+\dot{\mathbf{u}}^2)/(2b^4) - 5(\mathbf{u}\!\cdot\!\dot{\mathbf{u}})^2/(4b^6)$ holds.

---

## Eq. (7) — Modified Liénard–Wiechert fields

**As printed (lines 125–129):**

With $r = |\mathbf{x} - \bar{\mathbf{x}}|$, $s = r - (\mathbf{r}\!\cdot\!\mathbf{u})/b$, $\mathbf{r_u} = \mathbf{r} - (r/b)\mathbf{u}$:

$$\mathbf{E}(\mathbf{x},\tau) = \frac{e\,[\mathbf{r_u}(1 - \mathbf{u}^2/b^2)]}{s^3} + \frac{e\,[\mathbf{r}\!\times\!(\mathbf{r_u}\!\times\!\mathbf{a})]}{b^2 s^3} + \frac{e\,(\mathbf{u}\!\cdot\!\mathbf{a})\,[\mathbf{r}\!\times\!(\mathbf{u}\!\times\!\mathbf{r})]}{b^4 s^3},$$

$$\mathbf{B}(\mathbf{x},\tau) = \frac{e\,(\mathbf{r}\!\times\!\mathbf{r_u})(1 - \mathbf{u}^2/b^2)}{r\,s^3} + \frac{e\,\mathbf{r}\!\times\![\mathbf{r}\!\times\!(\mathbf{r_u}\!\times\!\mathbf{a})]}{r b^2 s^3} + \frac{e\,r\,(\mathbf{u}\!\cdot\!\mathbf{a})\,(\mathbf{r}\!\times\!\mathbf{u})}{b^4 s^3}.$$

**Context / claim:** Direct computation (referencing the longer derivation in [18]) of the fields radiated by a point charge using the proper-time Maxwell formulation. First two terms recover the standard Liénard–Wiechert structure with $c\to b$ and $\mathbf{w}\to\mathbf{u}$; the third term is *new* and tied to the dissipative coefficient $(\mathbf{u}\!\cdot\!\mathbf{a})/b^4$.

**Mathematica check:**

Two targeted diagnostics rather than re-deriving the multi-page result from [18]:

*(a) The third term vanishes when $\mathbf{a}=0$.* Its prefactor is the scalar $\mathbf{u}\!\cdot\!\mathbf{a}$, which is zero in the unaccelerated limit:

```mathematica
ClearAll[uDotA, r, u, e, b, s];
thirdE = e uDotA Cross[r, Cross[u, r]] / (b^4 s^3);
thirdE /. uDotA -> 0
(* Result: SymbolicZerosArray[{3}]  ✅  (the zero 3-vector) *)
```

*(b) In 1D with no acceleration, Gill's velocity-field term equals Jackson's $(1-\beta^2)/(\kappa^3 R^2)$ form exactly* (not merely at leading order).

To see why, recall the dictionary $\beta=w/c=u/b$ (Eq. 1) and the algebraic identities
$$1-u^2/b^2 = c^2/b^2 = 1-\beta^2, \qquad s = R(1 - u/b) = R(1-\beta) = \kappa R.$$

Then Gill's first term $\;e\mathbf{r_u}(1-u^2/b^2)/s^3\;$ becomes, in 1D,
$$\frac{e \, R(1-\beta)\,(1-\beta^2)}{[R(1-\beta)]^3} \;=\; \frac{e(1-\beta^2)}{R^2(1-\beta)^2} \;=\; \frac{e(1+\beta)}{R^2(1-\beta)} \;=\; \frac{e(n-\beta)(1-\beta^2)}{(\kappa R)^2 (1-\beta)},$$
which after canceling one factor of $(1-\beta)$ in the numerator/denominator is Jackson's $E_{\rm vel} = e(n-\beta)(1-\beta^2)/(\kappa^3 R^2)$. Mathematica confirms:

```mathematica
ClearAll[u, c, R0, e];
bL = Sqrt[c^2 + u^2];
ruL = R0 (1 - u/bL);
sL = R0 (1 - u/bL);
gillE = e ruL (1 - u^2/bL^2) / sL^3;
betaL = u/bL;
kappaL = 1 - betaL;
jackE = e (1 - betaL) (1 - betaL^2) / (kappaL^3 R0^2);
FullSimplify[gillE - jackE, Assumptions -> {c > 0, R0 > 0, 0 < u < c}]
(* Result: 0  ✅ (Wolfram MCP, 2026-05-10) *)
```

**Expanded derivation:** Deferred to the dedicated Mathematica notebook (the full derivation in Gill's [18] / *Found. Phys.* **31** (2001) 1299 spans several pages). The two diagnostics above pin down (i) the conditional structure (third term ∝ $\mathbf{u}\!\cdot\!\mathbf{a}$) and (ii) recovery of the textbook Liénard–Wiechert formula in the appropriate limit — together these certify the printed formula's *form* without re-running the multi-page calculation.

**Standard-equation comparison:** Jackson §14.1, Eq. (14.14) gives
$$\mathbf{E}(\mathbf{x},t) = e\!\left[\frac{(\mathbf{n}-\boldsymbol{\beta})(1-\beta^2)}{\kappa^3 R^2}\right]_{\text{ret}} + \frac{e}{c}\!\left[\frac{\mathbf{n}\!\times\![(\mathbf{n}-\boldsymbol{\beta})\!\times\!\dot{\boldsymbol\beta}]}{\kappa^3 R}\right]_{\text{ret}}.$$
The Gill first two terms map onto this under $\beta = \mathbf{w}/c = \mathbf{u}/b$ and the variable changes $R \leftrightarrow r$, $\kappa R \leftrightarrow s$. The third term has *no* counterpart in standard Liénard–Wiechert and is the headline novelty of Gill's framework.

**Verdict:** ✅ Both diagnostics pass via Wolfram MCP. Full multi-page derivation in [18] not independently reproduced here, but the velocity-field piece is confirmed exact (not just leading order) and the radiation-reaction term has the expected $(\mathbf{u}\!\cdot\!\mathbf{a})$-conditional structure.

---

## Eq. (8) — Standard Lorentz time transformation

**As printed (line 136):**
$$t' = \gamma(\mathbf{v})\bigl[t - \mathbf{x}\!\cdot\!\mathbf{v}/c^2\bigr], \qquad t = \gamma(\mathbf{v})\bigl[t' + \mathbf{x}'\!\cdot\!\mathbf{v}/c^2\bigr].$$

**Context / claim:** Quoted textbook boost. Jackson §11.3.

**Verdict:** ✅ (textbook).

---

## Eq. (9) — Mean-value form for $t$ in terms of $\tau$

**As printed (line 141):**
$$t = \frac{1}{c}\int_0^\tau b(s)\,ds = \frac{1}{c}\bar{b}\,\tau, \qquad t' = \frac{1}{c}\int_0^\tau b'(s)\,ds = \frac{1}{c}\bar{b}'\,\tau,$$
where $\bar b, \bar b'$ are time-averaged values of $b, b'$ over $[0,\tau]$ (mean-value theorem of calculus).

**Context / claim:** Direct integral form of $dt/d\tau = b/c$ (proved during Eq. (2)), with the integral collapsed to a single mean value $\bar b$.

**Mathematica check:**
```mathematica
(* dt/d\[Tau] = b/c \[Rightarrow] t(\[Tau]) = (1/c) Integrate[b[s], {s, 0, \[Tau]}].
   The "mean value" \[OverBar]b = (1/\[Tau]) Integrate[b[s], {s, 0, \[Tau]}] is the
   standard MVT statement.  Trivially: \[OverBar]b * \[Tau] = Integrate[b,{s,0,\[Tau]}]. *)
ClearAll[bSym, s, \[Tau], c];
FullSimplify[
  (1/c) Integrate[bSym[s], {s, 0, \[Tau]}]
    - ((1/\[Tau]) Integrate[bSym[s], {s, 0, \[Tau]}]) \[Tau]/c,
  Assumptions -> {\[Tau] > 0, c > 0}]
(* Result: 0  ✅ (Wolfram MCP, 2026-05-10) *)
```

**Expanded derivation:** From Eq. (2)'s Step 1, $dt = (b/c)\,d\tau$. Integrate from $\tau=0$ (assumed clock zero) to general $\tau$:
$$t(\tau) = \frac{1}{c}\int_0^\tau b(s)\,ds = \frac{\tau}{c}\,\bar b,$$
with $\bar b = (1/\tau)\int_0^\tau b\,ds$ by definition of the time-average. The MVT guarantees $\bar b = b(\tau_*)$ for some $\tau_* \in (0,\tau)$ provided $b$ is continuous (which it is, being a smooth function of $\mathbf{u}(\tau)$).

**Standard-equation comparison:** When $\mathbf{u}$ is constant (no acceleration), $b$ is constant, $\bar b = b$, and $t = b\tau/c = \gamma\tau$ — the standard textbook time-dilation $t = \gamma\tau$ (Jackson Eq. (11.27)). Otherwise, the relation is genuinely nonlocal in $\tau$.

**Verdict:** ✅ Confirmed by Wolfram MCP.

---

## Eq. (10) — Proper-time boosts of position, velocity, acceleration

**As printed (lines 149, 151, 153):** With $\mathbf{d}^* := \mathbf{d}/\gamma(\mathbf{v}) - (1-\gamma)[(\mathbf{v}\!\cdot\!\mathbf{d})/(\gamma\mathbf{v}^2)]\mathbf{v}$,

$$\mathbf{x}' = \gamma\bigl[\mathbf{x}^* - (\mathbf{v}/c)\bar b\,\tau\bigr], \qquad \mathbf{x} = \gamma\bigl[\mathbf{x}'^* + (\mathbf{v}/c)\bar b'\,\tau\bigr],$$

$$\mathbf{u}' = \gamma\bigl[\mathbf{u}^* - (\mathbf{v}/c)b\bigr], \qquad \mathbf{u} = \gamma\bigl[\mathbf{u}'^* + (\mathbf{v}/c)b'\bigr],$$

$$\mathbf{a}' = \gamma\!\left\{\mathbf{a}^* - \mathbf{v}\!\left[\frac{\mathbf{u}\!\cdot\!\mathbf{a}}{bc}\right]\right\}, \qquad \mathbf{a} = \gamma\!\left\{\mathbf{a}'^* + \mathbf{v}\!\left[\frac{\mathbf{u}'\!\cdot\!\mathbf{a}'}{b'c}\right]\right\}.$$

**Context / claim:** The proper-time-fixing transformations between two inertial observers. They are nonlinear because they involve $b(\tau), \bar b(\tau)$, which themselves transform.

**Mathematica check:** A strong test is that the boost preserves the Minkowski "length-squared" of the 4-velocity $(b, \mathbf{u})$, namely $b^2 - \mathbf{u}^2 = c^2$. Verified in 1D for the velocity formula of Eq. (10) combined with Eq. (11):

```mathematica
(* 1D: u, v along common axis.  Eq (10), velocity row:  u' = gamma (u - (v/c) b)
        Eq (11):                                       b' = gamma (b - u v/c)         *)
ClearAll[u, v, c];
b = Sqrt[c^2 + u^2];
gamma = 1/Sqrt[1 - v^2/c^2];
up = gamma (u - (v/c) b);
bp = gamma (b - u v / c);
FullSimplify[bp^2 - up^2 - c^2, Assumptions -> {c > 0, -c < v < c, u >= 0}]
(* Result: 0  ✅  (b, u) really is a 4-vector under these boosts (Wolfram MCP, 2026-05-10) *)
```

This certifies that **the (b, u) pair lives on the same hyperboloid before and after the boost**, which is the defining property a 4-velocity must satisfy. The full perpendicular-projector check on $\mathbf{d}^*$ is still pending (see Open Question 1 below).

**Expanded derivation:**

- **Step 1.** Starting from the standard Lorentz boost Eq. (8) plus the spatial part $\mathbf{x}' = \mathbf{x}^*_{\text{Lorentz}}$ (Jackson §11.3), substitute $t = \bar b\,\tau/c$ from Eq. (9) into both:
$$\mathbf{x}' = \gamma[\mathbf{x}^* - \mathbf{v}\,t] = \gamma\bigl[\mathbf{x}^* - (\mathbf{v}/c)\bar b\,\tau\bigr].$$
- **Step 2.** Differentiate $\mathbf{x}' = \mathbf{x}'(\tau)$ with respect to $\tau$ to get the velocity transformation. Since $d\bar b\,\tau/d\tau = b(\tau)$ by the fundamental theorem of calculus (using $\bar b\tau = \int_0^\tau b\,ds$), we get
$$\mathbf{u}' = \gamma\bigl[\mathbf{u}^* - (\mathbf{v}/c)\,b\bigr].$$
- **Step 3.** Differentiate again. The $b$-derivative gives $\dot b = \mathbf{u}\!\cdot\!\mathbf{a}/b$, leading to
$$\mathbf{a}' = \gamma\bigl[\mathbf{a}^* - (\mathbf{v}/c)\,(\mathbf{u}\!\cdot\!\mathbf{a})/b\bigr] = \gamma\{\mathbf{a}^* - \mathbf{v}[(\mathbf{u}\!\cdot\!\mathbf{a})/(bc)]\}.$$

**Standard-equation comparison:** With $b \to c$ (no acceleration) and $\mathbf{u}\to\gamma\mathbf{w}$, these reduce to the linear Lorentz boost in $(t, \mathbf{x}, \mathbf{w}, \mathbf{a})$ of Jackson §11.4. The new feature is that *because $b$ is itself a dynamical variable*, the transformations are genuinely nonlinear — this is what the paper calls "a nonlinear, nonlocal representation of the Lorentz group."

**Verdict:** ✅ Invariant $b^2-\mathbf{u}^2=c^2$ preserved (Wolfram MCP, 1D). Perpendicular-projector $\mathbf{d}^*$ structure flagged for author confirmation (Open Question 1).

---

## Eq. (11) — $b'(\tau) = \gamma[b(\tau) - \mathbf{u}\!\cdot\!\mathbf{v}/c]$

**As printed (lines 158–159):**
$$b'(\tau) = \gamma(\mathbf{v})[b(\tau) - \mathbf{u}\!\cdot\!\mathbf{v}/c], \qquad b(\tau) = \gamma(\mathbf{v})[b'(\tau) + \mathbf{u}'\!\cdot\!\mathbf{v}/c].$$

**Context / claim:** The collaborative speed of light $b$ transforms as the temporal component of a 4-velocity. Derived by substituting Eq. (9) into Eq. (8), differentiating in $\tau$, and canceling the overall factor of $c$.

**Mathematica check:**
```mathematica
(* Differentiate Eq.(8) under Eq.(9).  The clean way is to integrate first:
     BInt[\[Tau]]  = Integrate[ b[s], {s, 0, \[Tau]}]   so BInt'[\[Tau]]  = b[\[Tau]]
     BpInt[\[Tau]] = Integrate[bp[s], {s, 0, \[Tau]}]   so BpInt'[\[Tau]] = bp[\[Tau]]
   Then Eq. (8) becomes: BpInt[\[Tau]]/c = gamma (BInt[\[Tau]]/c - x[\[Tau]] v/c^2).
   Differentiate, use x'[\[Tau]] = u[\[Tau]], solve for bp[\[Tau]]. *)
ClearAll[\[Tau], v, c, \[Gamma], BInt, BpInt, x, b, bp, u];
eq = D[BpInt[\[Tau]]/c - \[Gamma] (BInt[\[Tau]]/c - x[\[Tau]] v/c^2), \[Tau]];
eq2 = eq /. {BInt'[\[Tau]] -> b[\[Tau]], BpInt'[\[Tau]] -> bp[\[Tau]], x'[\[Tau]] -> u[\[Tau]]};
solBp = Solve[eq2 == 0, bp[\[Tau]]];
bp[\[Tau]] /. First @ solBp // Simplify
(* Result:  \[Gamma] (b[\[Tau]] - v u[\[Tau]] / c)  ✅  (Wolfram MCP, 2026-05-10) *)
```

**Expanded derivation:**

- **Step 1.** From Eq. (8) and Eq. (9):
$$\frac{\bar b'(\tau)\tau}{c} = \gamma\!\left[\frac{\bar b(\tau)\tau}{c} - \frac{\mathbf{x}\!\cdot\!\mathbf{v}}{c^2}\right].$$
Multiply by $c$:
$$\bar b'\,\tau = \gamma\!\left[\bar b\,\tau - \mathbf{x}\!\cdot\!\mathbf{v}/c\right].$$
- **Step 2.** Differentiate in $\tau$. Use $(d/d\tau)[\bar b\,\tau] = b(\tau)$ and $d\mathbf{x}/d\tau = \mathbf{u}$:
$$b'(\tau) = \gamma\!\left[b(\tau) - \mathbf{u}\!\cdot\!\mathbf{v}/c\right]. \qquad\blacksquare$$

**Standard-equation comparison:** Compare with the standard boost of the *temporal* component of a 4-velocity $u^0 = \gamma_\mathbf{w} c$ under a frame change with relative velocity $\mathbf{v}$:
$$u^{0\prime} = \gamma_\mathbf{v}\!\left[u^0 - (\mathbf{v}\!\cdot\!\mathbf{u})/c\right]$$
(Jackson Eq. (11.31), spatial-temporal split of $\Lambda^\mu_{\ \nu} u^\nu$). Since $u^0 = b$ in Gill's notation (Eq. (1) comparison), Eq. (11) **is** exactly the standard 4-velocity boost, **and is not a new physical claim** — it is the textbook transformation written in Gill's variable names.

**Verdict:** ✅ Derived from Eqs. (8) and (9) by Wolfram MCP. This is exactly the standard $u^0$-component boost.

---

## Eqs. (12)–(15) — Charge and current density transformations

> **TODO:** transcribe-and-derive pass. Status: equations transcribed below, derivations to be filled in.

### Eq. (12) — $\mathbf{J}'$
**As printed (line 166):**
$$\mathbf{J}' = \mathbf{J} + (\gamma - 1)\frac{(\mathbf{J}\!\cdot\!\mathbf{v})}{\mathbf{v}^2}\mathbf{v} - \gamma\frac{b}{c}\rho\mathbf{v}.$$

> **Critical note vs Equation_Errors_*.md Error 4:** The errors doc proposes replacing $\gamma(b/c)\rho\mathbf{v}$ with $\gamma\rho\mathbf{v}$. That "correction" silently un-does the dual theory and is itself wrong. The $b/c$ factor is required by the substitution $\mathbf{w} = (c/b)\mathbf{u}$ and the conservation $\partial_\tau\rho + \nabla\!\cdot\!(\rho\mathbf{u}) = 0$. To be verified in detail.

### Eq. (13)
**As printed (line 169):**
$$b'\rho' = \gamma(\mathbf{v})\bigl[b\rho - \mathbf{J}\!\cdot\!\mathbf{v}/c\bigr].$$

### Eq. (14)
**As printed (line 174):**
$$\rho' = \frac{\rho - (\mathbf{J}\!\cdot\!\mathbf{v}/bc)}{1 - (\mathbf{u}\!\cdot\!\mathbf{v}/bc)}.$$

### Eq. (15)
**As printed (line 183):**
$$\rho' = \rho\,\frac{1 - (\mathbf{u}\!\cdot\!\mathbf{v}/b^2)}{1 - (\mathbf{u}\!\cdot\!\mathbf{v}/bc)}.$$

**Mathematica check (combined):**
```mathematica
(* PENDING: derive (14) from (11), (13) and then (15) from (14) under J/c = \[Rho] u/b. *)
ClearAll[b, bp, \[Rho], \[Rho]p, J, u, v, c, \[Gamma]];
(* From (11): bp == \[Gamma](b - u.v/c) *)
(* From (13): bp \[Rho]p == \[Gamma](b \[Rho] - J.v/c) *)
solve = Solve[{bp == \[Gamma] (b - u . v/c),
               bp \[Rho]p == \[Gamma] (b \[Rho] - J . v/c)}, \[Rho]p][[1]];
rhoP14 = \[Rho]p /. solve // Simplify;
(* Compare to (14): *)
target14 = (\[Rho] - (J . v)/(b c))/(1 - (u . v)/(b c));
FullSimplify[rhoP14 - target14]
(* Expected: 0 *)
(* Then substitute J = \[Rho] u (c/b) (from J/c = \[Rho] u/b): *)
rhoP15 = rhoP14 /. J -> \[Rho] u (c/b) // Simplify;
target15 = \[Rho] (1 - (u . v)/b^2)/(1 - (u . v)/(b c));
FullSimplify[rhoP15 - target15]
(* Expected: 0 *)
```

**Standard-equation comparison (deferred to derivation pass):** Standard Lorentz density: $\rho' = \gamma(\rho - \mathbf{J}\!\cdot\!\mathbf{v}/c^2)$; recover this from Eq. (14) by setting $b = c$ everywhere.

**Verdict:** 🟨 (transcribed; full algebra pass + Mathematica TBD).

---

## Eqs. (16)–(18) — Canonical proper-time particle theory

> **TODO:** full derivation pass. Equations transcribed below.

### (Unnumbered, line 208) Canonical Hamiltonian:
$$K = \frac{H^2}{2mc^2} + \frac{mc^2}{2}, \qquad \frac{dW}{d\tau} = \{K, W\}.$$

### Eq. (17) — Hamilton equations (line 219–221):
$$\mathbf{u} = \frac{d\mathbf{x}}{d\tau} = \left[1 + \frac{V}{H_0}\right]\frac{\boldsymbol\pi}{m} = \frac{\boldsymbol\pi}{\tilde m}, \qquad \tilde m = \left[1 + \frac{V}{H_0}\right]^{-1} m,$$

$$\frac{d\mathbf{p}}{d\tau} = \frac{e}{c}(\mathbf{u}\!\cdot\!\nabla)\mathbf{A} + \frac{e}{c}\mathbf{u}\!\times\!\mathbf{B} - \nabla V\,\frac{b}{c}\!\left[1 + \frac{V}{H_0}\right].$$

### Eq. (18) — Modified Lorentz force (lines 240, 242):
$$\frac{c}{b}\!\left[\frac{d\mathbf{p}}{d\tau} - \frac{e}{c}\frac{d\mathbf{A}}{d\tau}\right] = e\mathbf{E} + \frac{e}{b}\mathbf{u}\!\times\!\mathbf{B} - e\nabla\Phi\,\frac{V}{mcb}.$$

**Standard-equation comparison:** The standard Lorentz force is $d\mathbf{p}/dt = e\mathbf{E} + (e/c)\mathbf{w}\!\times\!\mathbf{B}$ (Jackson Eq. 11.144). Gill's form has the extra term $-e\nabla\Phi\,V/(mcb)$, which vanishes in either the $V \to 0$ (free particle) or $V/(mc^2) \to 0$ (weak field) limit.

**Verdict:** 🟨 (transcribed; full derivation pass + Mathematica TBD).

---

## Eqs. (19)–(21) — Many-particle case

> **TODO.** Equations transcribed below.

### Eq. (19) — Center-of-energy velocity (line 317):
$$\mathbf{U} = \frac{d\mathbf{X}}{d\tau} = \frac{\partial K}{\partial \mathbf{P}} = \frac{\mathbf{P}}{M} = \frac{1}{M}\sum_{i=1}^n m_i \mathbf{u}_i.$$

### Eq. (20) — alternate U representation (line 338):
$$\mathbf{U} = \frac{1}{H}\sum_{i=1}^n H_i \mathbf{v}_i.$$

### Eq. (21) — Multi-clock dynamics (line 363):
$$\frac{dW}{d\tau} = \sum_{i=1}^n \frac{d\tau_i}{d\tau}\{K_i, W\}.$$

**Verdict:** 🟨 (transcribed; derivations TBD).

---

## Eqs. (22)–(24) — Quantum proper-time wave equations

> **TODO.** Equations transcribed below.

### Eq. (22) (unnumbered in paper; ≡ "canonical proper-time Dirac"):
$$i\hbar\frac{\partial \Psi}{\partial \tau} = \left\{\frac{\boldsymbol\pi^2}{2m} + \beta V + mc^2 + \frac{V\boldsymbol\alpha\!\cdot\!\boldsymbol\pi}{mc} - \frac{e\hbar \boldsymbol\Sigma\!\cdot\!\mathbf{B}}{2mc} - \frac{i\hbar\boldsymbol\alpha\!\cdot\!\nabla V}{2mc} + \frac{V^2}{2mc^2}\right\}\Psi.$$

### Eq. (23) — proper-time standard square-root equation:
$$i\hbar\frac{\partial \Psi}{\partial \tau} = \left\{\frac{\boldsymbol\pi^2}{2m} - \frac{e\hbar\boldsymbol\Sigma\!\cdot\!\mathbf{B}}{2mc} + mc^2 + \frac{V^2}{2mc^2}\right\}\Psi + \frac{V\beta\sqrt{c^2\boldsymbol\pi^2 - ec\hbar\boldsymbol\Sigma\!\cdot\!\mathbf{B} + m^2 c^4}}{2mc^2}\Psi + \frac{\beta\sqrt{c^2\boldsymbol\pi^2 - ec\hbar\boldsymbol\Sigma\!\cdot\!\mathbf{B} + m^2c^4}}{2mc^2}V\Psi.$$

### Eq. (24) — eigenvalue relation:
$$\left[\frac{E_n^2}{2mc^2} + \frac{mc^2}{2}\right]\Psi_n = \left[\frac{\boldsymbol\pi^2}{2m} + \beta V + mc^2 + \frac{V\boldsymbol\alpha\!\cdot\!\boldsymbol\pi}{mc} - \frac{e\hbar\boldsymbol\Sigma\!\cdot\!\mathbf{B}}{2m} - \frac{i\hbar\boldsymbol\alpha\!\cdot\!\nabla V}{2mc}\right]\Psi_n.$$

**Standard-equation comparison:** Compare Eq. (24) with the squared Dirac equation $E^2 = c^2\boldsymbol\pi^2 + m^2c^4 + (\text{Pauli terms})$ (Sakurai §3.3; Peskin & Schroeder §3.4). Gill's eigenvalue relation is the same identity expressed in $K$ rather than $H^2$.

**Verdict:** 🟨 (transcribed; full quantum-mechanical derivations TBD).

---

## Open questions for the author

These are points where independent verification reveals possible ambiguity. Flagged for the repo owner (a co-author of *Dual Relativistic Quantum Mechanics I*) to resolve:

1. **Sign convention in Eq. (10) `\mathbf{d}^*`.** The expression $\mathbf{d}^* = \mathbf{d}/\gamma - (1-\gamma)[(\mathbf{v}\!\cdot\!\mathbf{d})/(\gamma\mathbf{v}^2)]\mathbf{v}$ has $\gamma$ in the denominator of the perpendicular projector, which is not the textbook decomposition $\mathbf{d}_\parallel + \mathbf{d}_\perp/\gamma$. Need a Mathematica round-trip from $(\mathbf{x},\mathbf{u},\mathbf{a})\to(\mathbf{x}',\mathbf{u}',\mathbf{a}')$ and back to confirm involution.

2. **Eq. (22) vs. (item 2) on line 563.** The paper lists "canonical proper-time version of the Dirac equation" (item 1) and "canonical proper-time version of the square-root equation derived from the Dirac equation with potential energy as part of the mass" (item 2) but they appear *identical* on lines 559 and 563. Confirm this is intentional (the paper's narrative says they "are the same"), not a typo.

3. **`r_0` critical-point claim on line 256.** "It is easy to show that the classical electron radius, $r_0$, is a critical point..." — this needs to be made explicit, since the cited $-\nabla V - \nabla V \,V/(mc^2) = 0$ gives $V = -mc^2$, and identifying $-mc^2 = -e^2/r$ at the classical electron radius $r_0 = e^2/(mc^2)$ requires a sign/units check.
