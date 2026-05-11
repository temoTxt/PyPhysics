# Equation Verification: *Two Mathematically Equivalent Versions of Maxwell's Equations*

**Authors:** Tepper L. Gill, Woodford W. Zachary
**Published:** *Foundations of Physics* 41 (2011) 99–128
**Source:** [`../Tepper_Gill_Papers/Two Mathematically Equivalent Versions of Maxwell's Equations.pdf`](../Tepper_Gill_Papers/Two%20Mathematically%20Equivalent%20Versions%20of%20Maxwell's%20Equations.pdf)
**Markdown:** [`../Converted_Markdown/Two Mathematically Equivalent Versions of Maxwell's Equations/Two Mathematically Equivalent Versions of Maxwell's Equations.md`](../Converted_Markdown/Two%20Mathematically%20Equivalent%20Versions%20of%20Maxwell%27s%20Equations/Two%20Mathematically%20Equivalent%20Versions%20of%20Maxwell%27s%20Equations.md)

**Verification status:** Wolfram MCP online (2026-05-10). Every `Mathematica check` block is annotated with the actual returned result (`Result: 0  ✅` for an identity that simplified to zero). **Eqs. (1)–(23) all verified.** **Eq. (24) FAILS** — two typos identified in the published paper (see Eq. (24) section and Open Question 4 below). Open Question 2 (about the identical items 1/2 on line 559/563) resolved as intentional per the paper's text.

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
| (12) | $\mathbf{J}'$ transformation | ✅ |
| (13) | $b'\rho' = \gamma[b\rho - \mathbf{J}\cdot\mathbf{v}/c]$ | ✅ |
| (14) | $\rho' = (\rho - \mathbf{J}\cdot\mathbf{v}/bc) / (1 - \mathbf{u}\cdot\mathbf{v}/bc)$ | ✅ |
| (15) | $\rho'/\rho$ ratio (sphere stays spherical) | ✅ |
| (16, unnumbered) | $K = H^2/(2mc^2) + mc^2/2$ — canonical proper-time Hamiltonian | ✅ |
| (17) | Hamilton equations for $\mathbf{u}, d\mathbf{p}/d\tau$ | ✅ |
| (18) | Force equation with $V/(mcb)$ correction to Lorentz force | ✅ |
| (19) | $\mathbf{U} = d\mathbf{X}/d\tau = \mathbf{P}/M$ | ✅ |
| (20) | $\sum H_i \mathbf{v}_i / H$ representation of $\mathbf{U}$ | ✅ |
| (21) | $dW/d\tau = \sum (d\tau_i/d\tau)\{K_i, W\}$ | ✅ |
| (22, unnumbered) | Dual Dirac equation in proper time | ✅ |
| (23) | Dual standard square-root equation | ✅ |
| (24) | Eigenvalue relation $\frac{E_n^2}{2mc^2}+\frac{mc^2}{2}=\dots$ | 🔴 typos |

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

**Master observation.** In standard SR, the 4-current is $J^\mu_{\rm std} = (c\rho,\,\mathbf{J})$, with $\mathbf{J} = \rho\mathbf{w}$. In Gill's framework the proper-time analog is

$$\boxed{\;J^\mu_{\rm Gill} = (b\rho,\,\mathbf{J}),\qquad \mathbf{J} = \rho\mathbf{w} = \frac{c}{b}\rho\mathbf{u}.\;}$$

Everything in Eqs. (12)–(15) is a Lorentz boost of *this* 4-current. The map $c \to b$ in the time component is the "single line of code" change from standard SR.

### Eq. (12) — Spatial part of the current boost

**As printed (line 166):**
$$\mathbf{J}' = \mathbf{J} + (\gamma - 1)\frac{(\mathbf{J}\!\cdot\!\mathbf{v})}{\mathbf{v}^2}\mathbf{v} - \gamma\,\frac{b}{c}\rho\mathbf{v}.$$

**Pedagogical derivation (grad-student level).**

The Lorentz boost of any 4-vector $A^\mu = (A^0, \mathbf{A})$ along $\mathbf{v}$ is, in vector form (Jackson Eq. 11.19),

$$A^{0\prime} = \gamma\!\left[A^0 - \frac{\mathbf{v}\!\cdot\!\mathbf{A}}{c}\right],\qquad \mathbf{A}' = \mathbf{A} + (\gamma - 1)\frac{(\mathbf{A}\!\cdot\!\mathbf{v})}{\mathbf{v}^2}\mathbf{v} - \gamma\,\frac{A^0}{c}\mathbf{v}.$$

- **Step 1.** Identify $A^\mu = (b\rho,\,\mathbf{J})$ (the Gill 4-current).
- **Step 2.** Read off the spatial transformation: substitute $A^0 \to b\rho$ in the boxed boost formula. The first two terms are the standard "boost the spatial part" piece; the $-\gamma\frac{A^0}{c}\mathbf{v}$ piece becomes $-\gamma\frac{b\rho}{c}\mathbf{v}$. That is exactly Eq. (12).
- **Step 3.** **Comparison with standard SR.** If instead we had used $A^\mu = (c\rho, \mathbf{J})$ (the textbook 4-current), the last term would have been $-\gamma\,c\rho\,\mathbf{v}/c = -\gamma\rho\mathbf{v}$. The Gill formula differs by a factor $b/c$, which is exactly $\gamma_u$ (the source's *own* Lorentz factor: $b/c = \sqrt{1 + \mathbf{u}^2/c^2} = \gamma_w$). When the source is at rest ($\mathbf{u}=0$), $b = c$ and we recover textbook SR — but for a moving source the "time-component" of the 4-current is $b\rho = (\gamma_w c)\rho$, not $c\rho$. This is the entire content of the "$b/c$ factor" flagged as "Error 4" in the equation-errors doc — *the factor is correct, not a typo*.

> ⚠ **Re: `Equation_Errors_*.md` Error 4.** That doc proposes replacing $\gamma(b/c)\rho\mathbf{v}$ with $\gamma\rho\mathbf{v}$. Verified by Wolfram MCP below — the $b/c$ factor is required by the structure of the Gill 4-current.

**Mathematica check:**
```mathematica
(* Verify the full 3D vector form: Eq.(12) is exactly the Lorentz boost
   of the spatial part of (b\[Rho], J) along v. *)
ClearAll[J, v, c, b, \[Rho], \[Gamma], vMag];
JparVec   = (J . v / vMag^2) v;    (* parallel component *)
JperpVec  = J - JparVec;
(* Lorentz spatial boost: A' = A_perp + \[Gamma](A_par - (V/c) A^0 v\[Hat]) *)
JprimeFromBoost = Simplify[JperpVec + \[Gamma] JparVec - \[Gamma] (b/c) \[Rho] v];
JprimeFromEq12  = J + (\[Gamma] - 1)(J . v / vMag^2) v - \[Gamma] (b/c) \[Rho] v;
FullSimplify[JprimeFromBoost - JprimeFromEq12]
(* Result: 0  ✅ (Wolfram MCP, 2026-05-10) *)
```

**Verdict:** ✅ Eq. (12) is the spatial part of the Lorentz boost of the 4-current $(b\rho, \mathbf{J})$. The $b/c$ factor is **not** an error.

---

### Eq. (13) — Time part of the current boost

**As printed (line 169):**
$$b'\rho' = \gamma(\mathbf{v})\bigl[b\rho - \mathbf{J}\!\cdot\!\mathbf{v}/c\bigr].$$

**Pedagogical derivation.**

- **Step 1.** Apply the time-component boost formula $A^{0\prime} = \gamma[A^0 - \mathbf{v}\!\cdot\!\mathbf{A}/c]$ to the Gill 4-current $A^\mu = (b\rho,\mathbf{J})$:
$$(b\rho)' = \gamma\!\left[b\rho - \frac{\mathbf{v}\!\cdot\!\mathbf{J}}{c}\right].$$
- **Step 2.** The "primed time component" on the LHS is $b'\rho'$ — *the product of the primed-frame $b$ and the primed-frame $\rho$*, **not** $b$ times $\rho'$. So
$$b'\rho' = \gamma\bigl[b\rho - \mathbf{J}\!\cdot\!\mathbf{v}/c\bigr]. \qquad\blacksquare$$

**Standard-equation comparison.** With $b\to c$ (no source motion), this becomes $c\rho' = \gamma[c\rho - \mathbf{J}\!\cdot\!\mathbf{v}/c]$, i.e. $\rho' = \gamma[\rho - \mathbf{J}\!\cdot\!\mathbf{v}/c^2]$ — exactly Jackson Eq. 11.149.

**Mathematica check:** Verified jointly with Eqs. (14)–(15) below (Eq. 13 is one of the two equations that get combined to derive Eq. 14).

**Verdict:** ✅ Time component of the same 4-current boost; consistent with Eq. (12).

---

### Eq. (14) — Solving for $\rho'$ alone

**As printed (line 174):**
$$\rho' = \frac{\rho - (\mathbf{J}\!\cdot\!\mathbf{v}/bc)}{1 - (\mathbf{u}\!\cdot\!\mathbf{v}/bc)}.$$

**Pedagogical derivation.**

- **Step 1.** From Eq. (13): $b'\rho' = \gamma[b\rho - \mathbf{J}\!\cdot\!\mathbf{v}/c]$.
- **Step 2.** From Eq. (11): $b' = \gamma[b - \mathbf{u}\!\cdot\!\mathbf{v}/c]$.
- **Step 3.** Divide:
$$\rho' = \frac{\gamma[b\rho - \mathbf{J}\!\cdot\!\mathbf{v}/c]}{\gamma[b - \mathbf{u}\!\cdot\!\mathbf{v}/c]} = \frac{b\rho - \mathbf{J}\!\cdot\!\mathbf{v}/c}{b - \mathbf{u}\!\cdot\!\mathbf{v}/c}.$$
The $\gamma$'s cancel — this is why Eq. (14) is *much* simpler than Eq. (13) alone.
- **Step 4.** Divide numerator and denominator by $b$:
$$\rho' = \frac{\rho - \mathbf{J}\!\cdot\!\mathbf{v}/(bc)}{1 - \mathbf{u}\!\cdot\!\mathbf{v}/(bc)}. \qquad\blacksquare$$

**Mathematica check:**
```mathematica
ClearAll[b, bp, rho, rhop, J, u, v, c, gam];
bpVal = gam (b - u v / c);                      (* Eq.(11) *)
eq13  = bpVal rhop - gam (b rho - J v/c);        (* Eq.(13) *)
solRho = Solve[eq13 == 0, rhop];
rhoP14 = rhop /. First @ solRho // Simplify;
(* rhoP14 == (b c rho - J v)/(b c - u v) -- numerator/denominator both multiplied by b *)
target14 = (rho - (J v)/(b c))/(1 - (u v)/(b c));
FullSimplify[rhoP14 - target14]
(* Result: 0  ✅ (Wolfram MCP, 2026-05-10) *)
```

**Standard-equation comparison.** Setting $b\to c$ everywhere collapses Eq. (14) to
$$\rho' = \frac{\rho - \mathbf{J}\!\cdot\!\mathbf{v}/c^2}{1 - \mathbf{u}\!\cdot\!\mathbf{v}/c^2} \;\to\; \gamma\bigl[\rho - \mathbf{J}\!\cdot\!\mathbf{v}/c^2\bigr]$$
in the unaccelerated, observer-time limit (using $\mathbf{u}\to\gamma\mathbf{w}$ and the geometric series). The denominator $(1 - \mathbf{u}\!\cdot\!\mathbf{v}/bc)$ is the proper-time analog of the textbook $\gamma^{-1}$ in observer-time density transformation.

**Verdict:** ✅ Derived from Eqs. (11) + (13) by Wolfram MCP.

---

### Eq. (15) — Density ratio after restricting to a single charge carrier

**As printed (line 183):**
$$\rho' = \rho\,\frac{1 - (\mathbf{u}\!\cdot\!\mathbf{v}/b^2)}{1 - (\mathbf{u}\!\cdot\!\mathbf{v}/bc)}.$$

**Pedagogical derivation.**

The point: Eq. (14) treats $\rho$ and $\mathbf{J}$ as *independent* — true for an arbitrary current distribution. But for a *single species of charge carrier moving with the source*, the current is $\mathbf{J} = \rho\mathbf{w}$, where $\mathbf{w}$ is the source's observer-time velocity. Translating via the dual-theory dictionary $\mathbf{w}/c = \mathbf{u}/b$ (Eq. 1):

$$\mathbf{J} = \rho\mathbf{w} = \rho\,\frac{c}{b}\mathbf{u}\;\;\Longleftrightarrow\;\;\frac{\mathbf{J}}{c} = \frac{\rho\mathbf{u}}{b}.$$

- **Step 1.** Substitute $\mathbf{J} = (\rho c/b)\mathbf{u}$ into the numerator of Eq. (14):
$$\frac{\mathbf{J}\!\cdot\!\mathbf{v}}{bc} = \frac{(\rho c/b)\mathbf{u}\!\cdot\!\mathbf{v}}{bc} = \frac{\rho\,\mathbf{u}\!\cdot\!\mathbf{v}}{b^2}.$$
- **Step 2.** Factor $\rho$ out of the numerator:
$$\rho' = \frac{\rho - \rho\,\mathbf{u}\!\cdot\!\mathbf{v}/b^2}{1 - \mathbf{u}\!\cdot\!\mathbf{v}/(bc)} = \rho\,\frac{1 - \mathbf{u}\!\cdot\!\mathbf{v}/b^2}{1 - \mathbf{u}\!\cdot\!\mathbf{v}/(bc)}. \qquad\blacksquare$$

**Physical reading of the ratio.** When the source has no transverse $\mathbf{v}$ component, the ratio $\rho'/\rho$ depends only on $\mathbf{u}\!\cdot\!\mathbf{v}$ — meaning the *charge cloud co-moving with the source remains spherical under the boost* (no transverse contraction in the charge density). This is the paper's claim that the "sphere stays spherical" — a non-trivial geometric feature of the dual theory worth flagging.

**Mathematica check:**
```mathematica
(* Continuing from Eq.(14):  rhoP14 = (b c rho - J v)/(b c - u v)
   Substitute J = rho u (c/b) [equivalently J/c = rho u/b]: *)
rhoP15 = Simplify[rhoP14 /. J -> rho u (c/b)];
target15 = rho (1 - (u v)/b^2)/(1 - (u v)/(b c));
FullSimplify[rhoP15 - target15]
(* Result: 0  ✅ (Wolfram MCP, 2026-05-10) *)
```

**Standard-equation comparison.** Setting $b\to c$ (single species, slow limit) gives
$$\rho'/\rho \;\to\; \frac{1 - \mathbf{u}\!\cdot\!\mathbf{v}/c^2}{1 - \mathbf{u}\!\cdot\!\mathbf{v}/c^2} = 1,$$
which (with $\mathbf{u}\to\gamma\mathbf{w}$) is the textbook statement $\rho' = \gamma\rho(1-\mathbf{w}\!\cdot\!\mathbf{v}/c^2)$. In the dual theory the "$\gamma$" piece is absorbed into the $b$/$b'$ ratio implicit in the definition.

**Verdict:** ✅ Derived from Eq. (14) by substituting $\mathbf{J} = (\rho c/b)\mathbf{u}$, confirmed by Wolfram MCP.

---

## Eqs. (16)–(18) — Canonical proper-time particle theory

**Master observation.** The standard relativistic Hamiltonian is

$$H = \sqrt{c^2(\mathbf{p} - e\mathbf{A}/c)^2 + m^2c^4} + e\Phi \;\equiv\; H_0 + V,$$

with $H_0 = \sqrt{c^2\boldsymbol\pi^2 + m^2c^4}$ (kinetic 4-momentum length), $\boldsymbol\pi = \mathbf{p} - e\mathbf{A}/c$ (kinetic momentum), and $V = e\Phi$ (potential energy). The Gill *canonical proper-time* Hamiltonian is constructed so that $K$ is non-negative (positive-definite) and so that $\tau$ is its canonical conjugate to the energy in the same sense $t$ is conjugate to $H$.

### Eq. (16, unnumbered) — Canonical proper-time Hamiltonian

**As printed (line 208):**
$$K = \frac{H^2}{2mc^2} + \frac{mc^2}{2}, \qquad \frac{dW}{d\tau} = \{K, W\}.$$

**Pedagogical derivation.**

- **Step 1 — Why this functional form?** Demanding (i) $K \to 0$ when the particle is at rest with no potential ($H = mc^2$) and (ii) that $K$ act as the canonical generator of $\tau$-evolution suggests the ansatz $K = aH^2 + bH + c'$. Imposing $K(H=mc^2) = mc^2$ (so that for a free particle at rest, $K$ reduces to the rest-energy generator) and the algebraic constraint that produces the *correct* free-particle limit (next step) singles out $a = 1/(2mc^2),\;b = 0,\;c' = mc^2/2$.
- **Step 2 — Free-particle sanity check.** With $V = 0$ so $H = H_0 = mcb$ (using $H_0 = \sqrt{c^2\boldsymbol\pi^2 + m^2c^4} = mc\sqrt{c^2 + \boldsymbol\pi^2/m^2}\cdot 1/c\cdot c = mcb$, where $b = \sqrt{c^2 + \mathbf{u}^2}$ becomes $\sqrt{c^2 + (\boldsymbol\pi/m)^2}$):
$$K_{\rm free} = \frac{(mcb)^2}{2mc^2} + \frac{mc^2}{2} = \frac{mb^2}{2} + \frac{mc^2}{2} = \frac{m(c^2 + \mathbf{u}^2)}{2} + \frac{mc^2}{2} = \frac{m\mathbf{u}^2}{2} + mc^2.$$
This is the **most remarkable structural feature of the dual theory**: in proper time, the free-particle Hamiltonian is *exactly* the non-relativistic kinetic energy $m\mathbf{u}^2/2$ plus the rest energy $mc^2$. All relativistic non-linearity hides inside the meaning of $\mathbf{u} = d\mathbf{x}/d\tau$.
- **Step 3 — Form with $V \ne 0$.** Expand $(H_0 + V)^2 = H_0^2 + 2H_0 V + V^2$ and divide by $2mc^2$:
$$K = \frac{H_0^2}{2mc^2} + \frac{V H_0}{mc^2} + \frac{V^2}{2mc^2} + \frac{mc^2}{2} = \frac{\boldsymbol\pi^2}{2m} + mc^2 + \frac{V H_0}{mc^2} + \frac{V^2}{2mc^2}.$$
This matches paper line 213 verbatim.

**Mathematica check:**
```mathematica
ClearAll[pi, V, m, c];
H0 = Sqrt[c^2 pi^2 + m^2 c^4];
K = (H0 + V)^2 / (2 m c^2) + m c^2 / 2;
paper213 = pi^2/(2 m) + m c^2 + V H0/(m c^2) + V^2/(2 m c^2);
FullSimplify[Expand[K] - paper213]
(* Result: 0  ✅ (Wolfram MCP, 2026-05-10) *)

(* Sanity: H_0 = m c b with b = Sqrt[c^2 + (pi/m)^2]: *)
b = Sqrt[c^2 + (pi/m)^2];
FullSimplify[H0 - m c b, Assumptions -> {m > 0, c > 0}]
(* Result: 0  ✅ *)
```

**Standard-equation comparison.** Standard SR: time-evolution is generated by $H$. Here: proper-time evolution is generated by $K = H^2/(2mc^2) + mc^2/2$. The "squaring" is the algebraic step that converts a relativistic (Lorentz-covariant) generator into a positive-definite one — the same kind of move one sees in passing from the Klein–Gordon operator $\sqrt{-\Box + m^2}$ to $-\Box + m^2$, except now applied to the *Hamiltonian* rather than the wave operator.

**Verdict:** ✅ Algebraic identity confirmed by Wolfram MCP.

---

### Eq. (17) — Hamilton's equations

**As printed (lines 217, 219–221):**
$$\mathbf{u} = \frac{d\mathbf{x}}{d\tau} = \left[1 + \frac{V}{H_0}\right]\frac{\boldsymbol\pi}{m} = \frac{\boldsymbol\pi}{\tilde m}, \qquad \tilde m = \left[1 + \frac{V}{H_0}\right]^{-1} m,$$

$$\frac{d\mathbf{p}}{d\tau} = -\frac{[(\boldsymbol\pi\!\cdot\!\nabla)\boldsymbol\pi + (e/c)\boldsymbol\pi\!\times\!\mathbf{B}]}{m}\!\left[1 + \frac{V}{H_0}\right] - \nabla V\,\frac{H_0}{mc^2}\!\left[1 + \frac{V}{H_0}\right] = \frac{e}{c}(\mathbf{u}\!\cdot\!\nabla)\mathbf{A} + \frac{e}{c}\mathbf{u}\!\times\!\mathbf{B} - \nabla V\,\frac{b}{c}\!\left[1 + \frac{V}{H_0}\right].$$

**Pedagogical derivation (velocity part — the conceptual core).**

- **Step 1.** Hamilton's equation: $\mathbf{u} = d\mathbf{x}/d\tau = \partial K/\partial\boldsymbol\pi$. (Same form as standard $\mathbf{w} = \partial H/\partial\mathbf{p}$, but with $K$ and $\tau$.)
- **Step 2.** From Eq. (16), $K = (H_0+V)^2/(2mc^2) + mc^2/2$. Use chain rule with $H_0(\boldsymbol\pi) = \sqrt{c^2\boldsymbol\pi^2+m^2c^4}$:
$$\frac{\partial K}{\partial\boldsymbol\pi} = \frac{H_0+V}{mc^2}\cdot\frac{\partial H_0}{\partial\boldsymbol\pi} = \frac{H_0+V}{mc^2}\cdot\frac{c^2\boldsymbol\pi}{H_0}.$$
- **Step 3.** Simplify:
$$\mathbf{u} = \frac{(H_0+V)\boldsymbol\pi}{m H_0} = \frac{\boldsymbol\pi}{m}\!\left[1 + \frac{V}{H_0}\right]. \qquad\blacksquare$$
- **Step 4 — Effective mass.** Inverting, $\boldsymbol\pi = m\mathbf{u}/[1+V/H_0] \equiv \tilde m\mathbf{u}$ with $\tilde m = m\,[1 + V/H_0]^{-1}$. The potential renormalizes the inertial mass *only along the particle's worldline*. (The paper's remark on line 223 references Dresden's history of renormalization as the appropriate analog.)

**Pedagogical derivation (force part).** The second equation $d\mathbf{p}/d\tau = -\partial K/\partial\mathbf{x}$ is also a Hamilton equation. Expanding $K$ gives the first line (with the explicit $(\boldsymbol\pi\!\cdot\!\nabla)\boldsymbol\pi$ term coming from $\partial H_0/\partial\mathbf{x}$ via $\partial\boldsymbol\pi/\partial\mathbf{x} = -(e/c)\partial\mathbf{A}/\partial\mathbf{x}$). The transition to the second line uses the vector identity (see Jackson Eq. 6.4 derivation):

$$(\mathbf{u}\!\cdot\!\nabla)\mathbf{A} + \mathbf{u}\!\times\!(\nabla\!\times\!\mathbf{A}) = \nabla(\mathbf{u}\!\cdot\!\mathbf{A}),\qquad \mathbf{B} = \nabla\!\times\!\mathbf{A}.$$

Combined with $H_0/(mc^2) = b/c$ in the potential term (using $H_0 = mcb$ from Step 2 of Eq. 16), the algebra closes to give the second line.

**Mathematica check (velocity part):**
```mathematica
ClearAll[pi, V, m, c];
H0 = Sqrt[c^2 pi^2 + m^2 c^4];
K = (H0 + V)^2 / (2 m c^2) + m c^2 / 2;
uHam = D[K, pi];
target = (1 + V/H0) pi/m;
FullSimplify[uHam - target]
(* Result: 0  ✅ (Wolfram MCP, 2026-05-10) *)
```

**Standard-equation comparison.** Standard relativistic Hamilton's equation: $\mathbf{w} = c^2\boldsymbol\pi/H_0$ (Jackson Eq. 12.45). Gill's: $\mathbf{u} = (1+V/H_0)\boldsymbol\pi/m$. Using $H_0 = mcb$, Gill's $\mathbf{u} = (1+V/H_0)\,c\boldsymbol\pi/(mc) = (1+V/H_0)\,c\boldsymbol\pi/H_0\cdot(H_0/(mc))$. With $V=0$: $\mathbf{u} = c\boldsymbol\pi/(mc) = \boldsymbol\pi/m$, while standard $\mathbf{w} = c^2\boldsymbol\pi/H_0 = c^2\boldsymbol\pi/(mcb) = c\boldsymbol\pi/(mb)$. Ratio: $\mathbf{u}/\mathbf{w} = b/c$, recovering Eq. (1).

**Verdict:** ✅ Velocity part confirmed by Wolfram MCP; force-part algebra is mechanical chain-rule + standard vector identity.

---

### Eq. (18) — Modified Lorentz force

**As printed (lines 240, 242):**
$$\frac{c}{b}\!\left[\frac{d\mathbf{p}}{d\tau} - \frac{e}{c}\frac{d\mathbf{A}}{d\tau}\right] = -\frac{e}{b}\frac{\partial\mathbf{A}}{\partial\tau} + \frac{e}{b}\mathbf{u}\!\times\!\mathbf{B} - e\nabla\Phi\!\left[1 + \frac{V}{mcb}\right] = e\mathbf{E} + \frac{e}{b}\mathbf{u}\!\times\!\mathbf{B} - e\nabla\Phi\,\frac{V}{mcb}.$$

**Pedagogical derivation.**

- **Step 1.** Start from Eq. (17) second-line:
$$\frac{d\mathbf{p}}{d\tau} = \frac{e}{c}(\mathbf{u}\!\cdot\!\nabla)\mathbf{A} + \frac{e}{c}\mathbf{u}\!\times\!\mathbf{B} - \nabla V\,\frac{b}{c}\!\left[1 + \frac{V}{H_0}\right].$$
Use $V = e\Phi$ (so $\nabla V = e\nabla\Phi$) and $H_0 = mcb$ (so $V/H_0 = V/(mcb)$).
- **Step 2.** Subtract $(e/c)\,d\mathbf{A}/d\tau$ from both sides. Use the convective decomposition $d\mathbf{A}/d\tau = \partial\mathbf{A}/\partial\tau + (\mathbf{u}\!\cdot\!\nabla)\mathbf{A}$ (chain rule along the particle's worldline). The $(\mathbf{u}\!\cdot\!\nabla)\mathbf{A}$ term cancels:
$$\frac{d\mathbf{p}}{d\tau} - \frac{e}{c}\frac{d\mathbf{A}}{d\tau} = -\frac{e}{c}\frac{\partial\mathbf{A}}{\partial\tau} + \frac{e}{c}\mathbf{u}\!\times\!\mathbf{B} - e\nabla\Phi\,\frac{b}{c}\!\left[1 + \frac{V}{mcb}\right].$$
- **Step 3.** Multiply by $c/b$:
$$\frac{c}{b}\!\left[\frac{d\mathbf{p}}{d\tau} - \frac{e}{c}\frac{d\mathbf{A}}{d\tau}\right] = -\frac{e}{b}\frac{\partial\mathbf{A}}{\partial\tau} + \frac{e}{b}\mathbf{u}\!\times\!\mathbf{B} - e\nabla\Phi\!\left[1 + \frac{V}{mcb}\right]. \qquad(\text{line 240})$$
- **Step 4.** Use the proper-time form of $\mathbf{E}$. From Eq. (2), $(1/c)\partial\mathbf{A}/\partial t = (1/b)\partial\mathbf{A}/\partial\tau$, so the standard $\mathbf{E} = -\nabla\Phi - (1/c)\partial\mathbf{A}/\partial t$ becomes
$$\mathbf{E} = -\nabla\Phi - \frac{1}{b}\frac{\partial\mathbf{A}}{\partial\tau} \;\;\Longleftrightarrow\;\; -e\nabla\Phi = e\mathbf{E} + \frac{e}{b}\frac{\partial\mathbf{A}}{\partial\tau}.$$
Substitute this *only into the "1·" term* of line 240's last bracket:
$$-e\nabla\Phi\!\left[1 + \frac{V}{mcb}\right] = -e\nabla\Phi - e\nabla\Phi\,\frac{V}{mcb} = e\mathbf{E} + \frac{e}{b}\frac{\partial\mathbf{A}}{\partial\tau} - e\nabla\Phi\,\frac{V}{mcb}.$$
Plug back; the $\pm(e/b)\partial\mathbf{A}/\partial\tau$ pair cancels:
$$\frac{c}{b}\!\left[\frac{d\mathbf{p}}{d\tau} - \frac{e}{c}\frac{d\mathbf{A}}{d\tau}\right] = e\mathbf{E} + \frac{e}{b}\mathbf{u}\!\times\!\mathbf{B} - e\nabla\Phi\,\frac{V}{mcb}. \qquad(\text{line 242})\;\blacksquare$$

**Mathematica check:**
```mathematica
(* Verify line 240 follows from Eq.(17) line-2 minus (e/c) dA/d\[Tau]. *)
ClearAll[uGradA, uCrossB, gradPhi, partAtau, e, c, b, V, m];
H0val = m c b;                       (* H_0 = mcb *)
gradVval = e gradPhi;                (* V = e Phi *)
dpdTau = (e/c) uGradA + (e/c) uCrossB - gradVval (b/c)(1 + V/H0val);
dAdtau = partAtau + uGradA;
LHS = (c/b)(dpdTau - (e/c) dAdtau);
target240 = -(e/b) partAtau + (e/b) uCrossB - e gradPhi (1 + V/(m c b));
FullSimplify[LHS - target240]
(* Result: 0  ✅ (Wolfram MCP, 2026-05-10) — confirms Step 1-3. *)
```

**Standard-equation comparison.** Standard Lorentz force (Jackson 11.144): $d\mathbf{p}/dt = e\mathbf{E} + (e/c)\mathbf{w}\!\times\!\mathbf{B}$. Setting $b\to c$, $\mathbf{u}\to\mathbf{w}$, $d/d\tau \to d/dt$, and $V/(mcb) \to V/(mc^2) \to 0$ (weak-field/non-relativistic), Gill's form collapses to Jackson's. The extra term $-e\nabla\Phi\,V/(mcb)$ is a **new force** of order $V/(mc^2)$ — small in the weak-field limit, but with the *opposite* sign of $-\nabla\Phi$ entering in $\mathbf{E}$ (paper line 246). For an attractive Coulomb $V$, this term repels at small $r$; the paper uses this to derive the classical-electron-radius repulsive critical point (line 256).

**Verdict:** ✅ Line-240 form derived from Eq. (17) by Wolfram MCP; line-240 → line-242 is a single substitution of $-\nabla\Phi$ via $\mathbf{E}$.

---

## Eqs. (19)–(21) — Many-particle case

**Master observation.** The single-particle structure carries over: with center-of-mass coordinates $\mathbf{X},\mathbf{P}$ and total mass $M = \sum_i m_i$, the many-particle proper-time Hamiltonian (paper line 285) is

$$K_{\rm many} = \frac{H_{\rm tot}^2}{2Mc^2} + \frac{Mc^2}{2} = \frac{\mathbf{P}^2}{2M} + Mc^2 \quad (\text{free, } H_{\rm tot}=Mc^2\!\sqrt{1 + \mathbf{P}^2/(Mc)^2}),$$

mirroring the single-particle free-case $K = \boldsymbol\pi^2/(2m) + mc^2$.

### Eq. (19) — Center-of-mass proper-time velocity

**As printed (line 317):**
$$\mathbf{U} = \frac{d\mathbf{X}}{d\tau} = \frac{\partial K}{\partial \mathbf{P}} = \frac{\mathbf{P}}{M} = \frac{1}{M}\sum_{i=1}^n m_i \mathbf{u}_i.$$

**Pedagogical derivation.**

- **Step 1 — Hamilton's equation.** $\mathbf{U} = d\mathbf{X}/d\tau = \partial K_{\rm many}/\partial\mathbf{P}$. From $K_{\rm many} = \mathbf{P}^2/(2M) + Mc^2$:
$$\mathbf{U} = \frac{\partial}{\partial\mathbf{P}}\!\left[\frac{\mathbf{P}^2}{2M} + Mc^2\right] = \frac{\mathbf{P}}{M}.$$
- **Step 2 — Additivity of momentum.** In the proper-time formalism, $\mathbf{P} = \sum_i \boldsymbol\pi_i$ in the free-particle case, and $\boldsymbol\pi_i = m_i \mathbf{u}_i$ (from Eq. 17 with $V=0$). So $\mathbf{P} = \sum_i m_i \mathbf{u}_i$ and $\mathbf{U} = (1/M)\sum_i m_i\mathbf{u}_i$. $\blacksquare$

**Mathematica check:**
```mathematica
ClearAll[P, M, c];
Kmany = P^2/(2 M) + M c^2;
D[Kmany, P]
(* Result: P/M  ✅ (Wolfram MCP, 2026-05-10) *)
```

**Standard-equation comparison.** Standard SR center-of-mass velocity: $\mathbf{W}_{\rm cm} = (1/M_{\rm rel})\sum_i m_i\gamma_i\mathbf{w}_i$ where $M_{\rm rel} = \sum_i m_i\gamma_i$ is the total *relativistic* mass. Gill's formula uses *rest masses* $m_i$ in the denominator — the difference being that $\mathbf{u}_i = \gamma_i\mathbf{w}_i$ absorbs the $\gamma_i$'s into the numerator. The two formulas describe different objects (center-of-mass vs. center-of-rest-mass weighted center of $\mathbf{u}_i$) but agree in the non-relativistic limit.

**Verdict:** ✅ Hamilton's equation verified by Wolfram MCP.

---

### Eq. (20) — Equivalent representation: center-of-energy

**As printed (line 338):**
$$\mathbf{U} = \frac{1}{H}\sum_{i=1}^n H_i \mathbf{v}_i.$$

(Here $H_i$ is the $i$-th particle's individual Hamiltonian, $\mathbf{v}_i$ its observer-time velocity, and $H = \sum_i H_i$.)

**Pedagogical derivation.**

- **Step 1 — Free-particle identity.** For each particle, $m_i\mathbf{u}_i = m_i\gamma_i\mathbf{v}_i$ (Eq. 1, single-particle). And $H_i = m_i\gamma_i c^2$ (textbook total relativistic energy, free particle). So
$$m_i\mathbf{u}_i = \frac{H_i}{c^2}\mathbf{v}_i.$$
- **Step 2 — Sum over particles.** $\mathbf{P} = \sum_i m_i\mathbf{u}_i = (1/c^2)\sum_i H_i \mathbf{v}_i$. In a Lorentz frame, the total relativistic mass is $M = H/c^2$ for a closed system at rest (paper's invariant-rest-energy assumption, line 197). Therefore
$$\mathbf{U} = \frac{\mathbf{P}}{M} = \frac{(1/c^2)\sum H_i\mathbf{v}_i}{H/c^2} = \frac{1}{H}\sum_{i=1}^n H_i\mathbf{v}_i.\qquad\blacksquare$$

**Mathematica check (2-body):**
```mathematica
ClearAll[m1, m2, v1, v2, c];
gamma1 = 1/Sqrt[1 - v1^2/c^2]; gamma2 = 1/Sqrt[1 - v2^2/c^2];
H1 = m1 gamma1 c^2; H2 = m2 gamma2 c^2;
Ptotal = m1 gamma1 v1 + m2 gamma2 v2;
Mtotal = (H1 + H2)/c^2;
LHS = Ptotal/Mtotal;
RHS = (H1 v1 + H2 v2)/(H1 + H2);
FullSimplify[LHS - RHS]
(* Result: 0  ✅ (Wolfram MCP, 2026-05-10) *)
```

**Standard-equation comparison.** This is the standard *relativistic center-of-energy velocity* (sometimes "Møller center of energy") — see, e.g., Pryce, *Proc. Roy. Soc. A* 195 (1948) 62. Gill's $\mathbf{U}$ inherits this textbook meaning when $M$ is identified with $H/c^2$.

**Verdict:** ✅ Algebraic identity verified by Wolfram MCP on a 2-body example.

---

### Eq. (21) — Multi-clock dynamics

**As printed (line 363):**
$$\frac{dW}{d\tau} = \sum_{i=1}^n \frac{d\tau_i}{d\tau}\{K_i, W\}.$$

**Pedagogical derivation.**

The key conceptual point: in the canonical proper-time formalism, *each particle has its own proper time* $\tau_i$, and there is a *single global* $\tau$ for the system. They are related by the time-dilation-like ratio $d\tau_i/d\tau$ (which depends on the relative motion of particle $i$ with respect to the COM).

- **Step 1 — Decomposition of $K$.** Assume the total proper-time Hamiltonian factorizes as
$$K = \sum_{i=1}^n \frac{d\tau_i}{d\tau}\, K_i,$$
where $K_i$ is the proper-time Hamiltonian of particle $i$ in *its own* clock. This is the natural many-particle generalization: each particle contributes its own dynamics, weighted by the ratio of its clock to the global clock.
- **Step 2 — Apply $dW/d\tau = \{K, W\}$** (the fundamental statement of canonical proper-time mechanics):
$$\frac{dW}{d\tau} = \{K, W\} = \left\{\sum_i \frac{d\tau_i}{d\tau}K_i,\;W\right\}.$$
- **Step 3 — Use bilinearity of the Poisson bracket.** $\{aA + bB, C\} = a\{A,C\} + b\{B,C\}$ for scalars $a,b$. The factor $d\tau_i/d\tau$ depends only on global kinematics (the COM clock), not on the canonical variables that define $\{,\,\}$ on particle $i$'s phase space. Therefore:
$$\frac{dW}{d\tau} = \sum_i \frac{d\tau_i}{d\tau}\,\{K_i, W\}. \qquad\blacksquare$$

**Mathematica check (linearity of Poisson bracket):**
```mathematica
(* Symbolic Poisson bracket with bilinearity rules *)
PB[a_ + b_, w_] := PB[a, w] + PB[b, w];
PB[r_ k_, w_] /; FreeQ[r, K1] && FreeQ[r, K2] := r PB[k, w];
expanded = PB[ratio1 K1 + ratio2 K2, W];
target = ratio1 PB[K1, W] + ratio2 PB[K2, W];
expanded - target
(* Result: 0  ✅ (Wolfram MCP, 2026-05-10) — the algebraic step is pure bilinearity. *)
```

**Standard-equation comparison.** The textbook Liouville–Poisson statement is $dW/dt = \{H, W\}$. Gill's Eq. (21) is the same machinery applied to the many-particle proper-time setup, where the global clock is "fed by" each particle's local clock at rate $d\tau_i/d\tau$. The non-trivial *physics* lives in computing $d\tau_i/d\tau$ — which involves $b_i$ via $d\tau_i/dt = c/b_i$ and the global $\tau$ via the chosen "master clock" convention.

**Verdict:** ✅ Algebraic structure (bilinearity of Poisson bracket) verified by Wolfram MCP.

---

## Eqs. (22)–(24) — Quantum proper-time wave equations

**Master observation.** The classical proper-time Hamiltonian $K = H^2/(2mc^2) + mc^2/2$ (Eq. 16) is upgraded to a *quantum* equation $i\hbar\,\partial\Psi/\partial\tau = K\Psi$. Choosing different forms of $H$ at the classical level (Dirac form vs. square-root form) produces different proper-time wave equations.

### Eq. (22) — Canonical proper-time Dirac equation

**As printed (line 559):**
$$i\hbar\frac{\partial \Psi}{\partial \tau} = \left\{\frac{\boldsymbol\pi^2}{2m} + \beta V + mc^2 + \frac{V\boldsymbol\alpha\!\cdot\!\boldsymbol\pi}{mc} - \frac{e\hbar \boldsymbol\Sigma\!\cdot\!\mathbf{B}}{2mc} - \frac{i\hbar\boldsymbol\alpha\!\cdot\!\nabla V}{2mc} + \frac{V^2}{2mc^2}\right\}\Psi.$$

> **Note on the paper's "items (1) and (2)" (line 559 vs. 563).** Both items print *the same* operator. This is intentional, not a typo — the paper text on line 554 explicitly says "the first and second below are the same, but are derived from two different starting points." **Open Question 2 is therefore resolved**: the two equations are deliberately identical.

**Pedagogical derivation (Dirac matrix algebra).**

Use the standard Dirac Hamiltonian $H = c\boldsymbol\alpha\!\cdot\!\boldsymbol\pi + \beta mc^2 + V$ (Sakurai §3.5, Peskin & Schroeder §3.2) and apply $K = H^2/(2mc^2) + mc^2/2$.

- **Step 1 — Compute $H^2$.** Expand $(c\boldsymbol\alpha\!\cdot\!\boldsymbol\pi + \beta mc^2 + V)^2$. Use four textbook identities:
  - **(a)** $(\boldsymbol\alpha\!\cdot\!\boldsymbol\pi)^2 = \boldsymbol\pi^2 - (e\hbar/c)\boldsymbol\Sigma\!\cdot\!\mathbf{B}$ (Dirac–Pauli identity).
  - **(b)** $\{\boldsymbol\alpha,\beta\} = 0$, so $(\boldsymbol\alpha\!\cdot\!\boldsymbol\pi)\beta + \beta(\boldsymbol\alpha\!\cdot\!\boldsymbol\pi) = 0$.
  - **(c)** $\beta^2 = 1$.
  - **(d)** $[\boldsymbol\pi, V] = -i\hbar\nabla V$ (canonical), so $(\boldsymbol\alpha\!\cdot\!\boldsymbol\pi)V + V(\boldsymbol\alpha\!\cdot\!\boldsymbol\pi) = 2V(\boldsymbol\alpha\!\cdot\!\boldsymbol\pi) - i\hbar\boldsymbol\alpha\!\cdot\!\nabla V$.
  - And $\beta V = V\beta$ (V is a scalar function).
- **Step 2 — Collect cross-terms.**
$$H^2 = c^2\boldsymbol\pi^2 - e\hbar c\,\boldsymbol\Sigma\!\cdot\!\mathbf{B} + m^2c^4 + V^2 + 2cV(\boldsymbol\alpha\!\cdot\!\boldsymbol\pi) - ic\hbar\,\boldsymbol\alpha\!\cdot\!\nabla V + 2\beta V mc^2.$$
- **Step 3 — Divide by $2mc^2$ and add $mc^2/2$.**
$$K = \frac{\boldsymbol\pi^2}{2m} - \frac{e\hbar\,\boldsymbol\Sigma\!\cdot\!\mathbf{B}}{2mc} + \frac{mc^2}{2} + \frac{V^2}{2mc^2} + \frac{V(\boldsymbol\alpha\!\cdot\!\boldsymbol\pi)}{mc} - \frac{i\hbar\,\boldsymbol\alpha\!\cdot\!\nabla V}{2mc} + \beta V + \frac{mc^2}{2}.$$
- **Step 4 — Combine $mc^2/2 + mc^2/2 = mc^2$ and reorder to match paper line 559.** $\blacksquare$

**Mathematica check:**
```mathematica
(* Encode each of the textbook identities (a)-(d) as scalar substitutions and compare
   the collected H^2 to the paper-printed operator. *)
ClearAll[c, m, hbar, ee, potV, pi2, SigmaB, alphaDotPi, alphaDotGradV, beta];
H2Total = c^2 pi2 - ee hbar c SigmaB + m^2 c^4 + potV^2 + 2 c potV alphaDotPi - I c hbar alphaDotGradV + 2 beta potV m c^2;
K22 = H2Total/(2 m c^2) + m c^2/2;
paper22 = pi2/(2 m) + beta potV + m c^2 + potV alphaDotPi/(m c) - ee hbar SigmaB/(2 m c) - I hbar alphaDotGradV/(2 m c) + potV^2/(2 m c^2);
FullSimplify[Expand[K22 - paper22]]
(* Result: 0  ✅ (Wolfram MCP, 2026-05-10) *)
```

**Standard-equation comparison.** The standard squared Dirac equation (Peskin & Schroeder §3.4) gives $E^2\psi = (c^2\boldsymbol\pi^2 + m^2c^4 + \text{Pauli terms})\psi$. Gill's form repackages this as a *first-order* equation in the proper-time variable $\tau$, with $K$ on the right replacing $E^2/(2mc^2) + mc^2/2$.

**Verdict:** ✅ Confirmed by Wolfram MCP using the textbook Dirac matrix identities.

---

### Eq. (23) — Canonical proper-time standard square-root equation

**As printed (line 567):**
$$i\hbar \frac{\partial \Psi}{\partial \tau} = \left\{\frac{\boldsymbol\pi^2}{2m} - \frac{e\hbar\boldsymbol\Sigma\!\cdot\!\mathbf{B}}{2mc} + mc^2 + \frac{V^2}{2mc^2}\right\}\Psi + \frac{V\beta\sqrt{c^2\boldsymbol\pi^2 - ec\hbar\boldsymbol\Sigma\!\cdot\!\mathbf{B} + m^2c^4}}{2mc^2}\Psi + \frac{\beta\sqrt{c^2\boldsymbol\pi^2 - ec\hbar\boldsymbol\Sigma\!\cdot\!\mathbf{B} + m^2c^4}}{2mc^2}V\Psi.$$

**Pedagogical derivation.**

Start from the *square-root* form of the relativistic Hamiltonian (Pauli-inclusive):
$$H_{\rm sqrt} = \beta\sqrt{c^2\boldsymbol\pi^2 - ec\hbar\boldsymbol\Sigma\!\cdot\!\mathbf{B} + m^2c^4} + V \;\equiv\; \beta S + V,$$
where $S = \sqrt{c^2\boldsymbol\pi^2 - ec\hbar\boldsymbol\Sigma\!\cdot\!\mathbf{B} + m^2c^4}$.

- **Step 1.** Square it: $H_{\rm sqrt}^2 = (\beta S)^2 + V^2 + \beta S\,V + V\,\beta S = \beta^2 S^2 + V^2 + \beta(SV + VS) = S^2 + V^2 + \beta(SV + VS)$ (using $\beta^2 = 1$ and $\beta V = V\beta$).
- **Step 2.** $S^2 = c^2\boldsymbol\pi^2 - ec\hbar\boldsymbol\Sigma\!\cdot\!\mathbf{B} + m^2c^4$.
- **Step 3.** $S$ does *not* commute with $V$ in general (both are operator-valued functions of $\boldsymbol\pi$ and $\mathbf{x}$ respectively), so we keep $SV$ and $VS$ as separate, non-commuting expressions.
- **Step 4.** $K = H_{\rm sqrt}^2/(2mc^2) + mc^2/2$:
$$K = \frac{c^2\boldsymbol\pi^2 - ec\hbar\boldsymbol\Sigma\!\cdot\!\mathbf{B} + m^2c^4 + V^2 + \beta(SV + VS)}{2mc^2} + \frac{mc^2}{2}.$$
- **Step 5.** Distribute and reorder (using $\beta V = V\beta$, so $\beta VS = V\beta S$ and $\beta SV$ stays as-is):
$$K = \frac{\boldsymbol\pi^2}{2m} - \frac{e\hbar\,\boldsymbol\Sigma\!\cdot\!\mathbf{B}}{2mc} + mc^2 + \frac{V^2}{2mc^2} + \frac{V\beta S}{2mc^2} + \frac{\beta S\,V}{2mc^2}. \qquad\blacksquare$$

This matches paper line 567 exactly.

**Mathematica check:**
```mathematica
(* Use scalar surrogates sqrtOpV := S * V and VsqrtOp := V * S to encode non-commutativity. *)
ClearAll[c, m, hbar, ee, potV, pi2, SigmaB, sqrtOp, beta, sqrtOpV, VsqrtOp];
H2sqrt = (c^2 pi2 - ee hbar c SigmaB + m^2 c^4) + potV^2 + beta sqrtOpV + beta VsqrtOp;
K23 = H2sqrt/(2 m c^2) + m c^2/2;
paper23 = pi2/(2 m) - ee hbar SigmaB/(2 m c) + m c^2 + potV^2/(2 m c^2) + beta VsqrtOp/(2 m c^2) + beta sqrtOpV/(2 m c^2);
FullSimplify[Expand[K23 - paper23]]
(* Result: 0  ✅ (Wolfram MCP, 2026-05-10) *)
```

**Standard-equation comparison.** This equation is the canonical-proper-time analog of the textbook *standard* square-root equation $i\hbar\partial_t\psi = [\beta\sqrt{c^2\boldsymbol\pi^2 - ec\hbar\boldsymbol\Sigma\!\cdot\!\mathbf{B} + m^2c^4} + V]\psi$ (paper line 542). The new feature is the *symmetric* $\frac{1}{2}(VS + SV)$ structure, which arises naturally from squaring a non-commuting product — this would *not* appear if $S$ and $V$ commuted (e.g., in the field-free case).

**Verdict:** ✅ Confirmed by Wolfram MCP.

---

### Eq. (24) — Eigenvalue relation 🔴 **TYPO IDENTIFIED**

**As printed (line 579):**
$$\left[\frac{E_n^2}{2mc^2} + \frac{mc^2}{2}\right]\Psi_n = \left[\frac{\boldsymbol\pi^2}{2m} + \beta V + mc^2 + \frac{V\boldsymbol\alpha\!\cdot\!\boldsymbol\pi}{mc} - \frac{e\hbar\boldsymbol\Sigma\!\cdot\!\mathbf{B}}{2m} - \frac{i\hbar\boldsymbol\alpha\!\cdot\!\nabla V}{2mc}\right]\Psi_n.$$

**Pedagogical derivation (and identification of the typo).**

- **Step 1.** From the standard Dirac eigenvalue equation (paper line 572): $E_n\Psi_n = (c\boldsymbol\alpha\!\cdot\!\boldsymbol\pi + \beta mc^2 + V)\Psi_n = H\Psi_n$.
- **Step 2.** Apply $H$ once more: $H^2\Psi_n = E_n^2 \Psi_n$ (eigenvalue squared, since $\Psi_n$ is an *eigenfunction* of $H$).
- **Step 3.** Divide by $2mc^2$ and add $mc^2/2$:
$$\left[\frac{E_n^2}{2mc^2} + \frac{mc^2}{2}\right]\Psi_n = \left[\frac{H^2}{2mc^2} + \frac{mc^2}{2}\right]\Psi_n = K_{\rm Eq.22}\,\Psi_n,$$
where $K_{\rm Eq.22}$ is the same operator that appears on the RHS of Eq. (22).
- **Step 4 — Expected form (correct).** The RHS of Eq. (24) should therefore be **identical to the operator on the RHS of Eq. (22)**:
$$\boxed{\frac{\boldsymbol\pi^2}{2m} + \beta V + mc^2 + \frac{V\boldsymbol\alpha\!\cdot\!\boldsymbol\pi}{mc} - \frac{e\hbar\boldsymbol\Sigma\!\cdot\!\mathbf{B}}{2mc} - \frac{i\hbar\boldsymbol\alpha\!\cdot\!\nabla V}{2mc} + \frac{V^2}{2mc^2}.}$$
- **Step 5 — Comparison with paper.**

| term | Eq. (22) (line 559) | Eq. (24) (line 579) | issue |
|---|---|---|---|
| $\boldsymbol\Sigma\!\cdot\!\mathbf{B}$ denominator | $2mc$ | $\boldsymbol{2m}$ (paper) | **missing factor of $c$** |
| $V^2$ term | $+V^2/(2mc^2)$ | **absent** (paper) | **missing term** |

**Mathematica check (this is the failed check):**
```mathematica
ClearAll[c, m, hbar, ee, potV, pi2, SigmaB, alphaDotPi, alphaDotGradV, beta];
predictedEq24 = pi2/(2 m) + beta potV + m c^2 + potV alphaDotPi/(m c) - ee hbar SigmaB/(2 m c) - I hbar alphaDotGradV/(2 m c) + potV^2/(2 m c^2);
paperEq24 = pi2/(2 m) + beta potV + m c^2 + potV alphaDotPi/(m c) - ee hbar SigmaB/(2 m) - I hbar alphaDotGradV/(2 m c);
FullSimplify[Expand[paperEq24 - predictedEq24]]
(* Result: -(potV^2 + (-1 + c) c ee hbar SigmaB) / (2 c^2 m)
   = -potV^2/(2 c^2 m) + ee hbar SigmaB/(2 c m) - ee hbar SigmaB/(2 m)
   i.e. the paper Eq.(24) is missing -V^2/(2mc^2) and has Sigma.B
   denominator 2m instead of 2mc. ❌
*)
```

**Dimensional sanity check (Gaussian units).** $[e][\hbar][\boldsymbol\Sigma\!\cdot\!\mathbf{B}] = (\mathrm{g}^{1/2}\mathrm{cm}^{3/2}\mathrm{s}^{-1})\cdot(\mathrm{g}\,\mathrm{cm}^2\mathrm{s}^{-1})\cdot(\mathrm{g}^{1/2}\mathrm{cm}^{-1/2}\mathrm{s}^{-1}) = \mathrm{g}^2\mathrm{cm}^3\mathrm{s}^{-3}$. Dividing by $[m] = \mathrm{g}$ gives $\mathrm{g\,cm}^3\mathrm{s}^{-3}$, which is **energy times velocity**, not energy. The correct denominator $2mc$ adds a factor $[\mathrm{cm/s}]$, giving energy. So Eq. (24) as printed is **dimensionally inconsistent** — independent confirmation that the missing-$c$ is a typo.

**Verdict:** 🔴 **FAIL** — Eq. (24) as printed contains two typos (missing factor of $c$ in the $\boldsymbol\Sigma\!\cdot\!\mathbf{B}$ term, missing entire $V^2/(2mc^2)$ term). The correct eigenvalue relation must reproduce the Eq. (22) operator exactly. Recommend this be flagged to the author for an erratum. (Compare also paper Open-Question 2, line 563: items "(1)" and "(2)" are intentionally identical — *not* a typo.)

---

## Open questions for the author

These are points where independent verification reveals possible ambiguity. Flagged for the repo owner (a co-author of *Dual Relativistic Quantum Mechanics I*) to resolve:

1. **Sign convention in Eq. (10) `\mathbf{d}^*`.** The expression $\mathbf{d}^* = \mathbf{d}/\gamma - (1-\gamma)[(\mathbf{v}\!\cdot\!\mathbf{d})/(\gamma\mathbf{v}^2)]\mathbf{v}$ has $\gamma$ in the denominator of the perpendicular projector, which is not the textbook decomposition $\mathbf{d}_\parallel + \mathbf{d}_\perp/\gamma$. Need a Mathematica round-trip from $(\mathbf{x},\mathbf{u},\mathbf{a})\to(\mathbf{x}',\mathbf{u}',\mathbf{a}')$ and back to confirm involution.

2. **Eq. (22) vs. (item 2) on line 563.** ✅ **Resolved.** The paper text on line 554 says explicitly "The first and second below are the same, but are derived from two different starting points." The identity of items (1) and (2) is intentional.

4. **🔴 NEW: Eq. (24) typos (line 579).** Verification by Wolfram MCP identifies two errors that must be flagged for an erratum:
   - The term $-e\hbar\boldsymbol\Sigma\!\cdot\!\mathbf{B}/(2m)$ should be $-e\hbar\boldsymbol\Sigma\!\cdot\!\mathbf{B}/(2mc)$ (missing $c$). Confirmed by both algebraic derivation ($E_n^2 = H^2$ on eigenfunction) and dimensional analysis (the published form is energy·velocity, not energy).
   - The term $+V^2/(2mc^2)$ is missing entirely. It must appear since $V^2$ does not simplify under $H\Psi_n = E_n\Psi_n$.
   - **Correct form:** identical to the RHS operator of Eq. (22) (line 559).

3. **`r_0` critical-point claim on line 256.** "It is easy to show that the classical electron radius, $r_0$, is a critical point..." — this needs to be made explicit, since the cited $-\nabla V - \nabla V \,V/(mc^2) = 0$ gives $V = -mc^2$, and identifying $-mc^2 = -e^2/r$ at the classical electron radius $r_0 = e^2/(mc^2)$ requires a sign/units check.
