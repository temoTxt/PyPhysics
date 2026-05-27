# 01 — Setup and verification of the modified-Newtonian acceleration

**Status:** ✅ all Wolfram-MCP verified.
**Source paper:** [`../Tepper_Gill_Papers/Dual Newtonian Theory.tex`](../Tepper_Gill_Papers/Dual%20Newtonian%20Theory.tex), §1 *Dual Classical Theory*.

## 1. The framework's interaction Hamiltonian

The paper begins with the dual proper-time Hamiltonian (Eq. 2.4, identical to DRQM I Eq. I.6 verified ✅ at [`../Equation_Verification/Dual_Relativistic_Quantum_Mechanics_I.md` line 124](../Equation_Verification/Dual_Relativistic_Quantum_Mechanics_I.md)):

$$K = \frac{H^2}{2mc^2} + \frac{mc^2}{2}.$$

For interaction with a potential `V`, the paper chooses `H₁ = √(c²π² + (mc² + V)²)` (the "potential-in-the-mass" form, identical to DRQM I §II.3 verified ✅). Substituting and dropping `O(c⁻²)` cross-terms beyond `1/c²`:

$$K \;=\; \frac{\boldsymbol{\pi}^2}{2m} + mc^2 + V + \frac{V^2}{2mc^2}.$$

This is the load-bearing structural input to the rest of the paper.

## 2. Eq. (h3) — modified-Newtonian acceleration

Starting from `K` above with vanishing vector potential (`π = p`), Hamilton's equation `ṗ = −∂K/∂x` gives:

$$\frac{\partial K}{\partial x} \;=\; \nabla V + \frac{V \nabla V}{mc^2} \;=\; \nabla V\!\left(1 + \frac{V}{mc^2}\right).$$

Dividing by `m`:

$$\boxed{\;\mathbf{a} \;=\; -\frac{\nabla V}{m}\!\left(1 + \frac{V}{mc^2}\right).\;}\quad\text{(paper Eq. h3)}$$

**Wolfram MCP check:**

```wolfram
ClearAll[r, cc, mm, GG, MM, potV, KK, gradV, accel]; KK[r_] = potV[r] + potV[r]^2/(2 mm cc^2); dKdr = D[KK[r], r]; gradV = D[potV[r], r]; Print["dK/dr = ", FullSimplify[dKdr]]; Print["FullSimplify[dKdr - ∇V·(1 + V/(mc²))] = ", FullSimplify[dKdr - gradV (1 + potV[r]/(mm cc^2))]]
```

**Result:** `dK/dr = (1 + potV[r]/(cc² mm)) potV′[r]`; `FullSimplify[dKdr - ∇V·(1+V/(mc²))] = 0`. ✅

## 3. Eq. (h4) — two-body equations of motion under gravity

Apply Eq. (h3) with `V = −GMm/r` (gravitational potential between Sun mass `M` and test particle Mercury mass `m`):

$$\nabla V \;=\; \frac{\partial V}{\partial r}\,\hat{\mathbf{e}}_r \;=\; +\frac{GMm}{r^2}\,\hat{\mathbf{e}}_r, \qquad \frac{V}{mc^2} \;=\; -\frac{GM}{c^2 r}.$$

Substituting:

$$\boxed{\;\mathbf{a}_m \;=\; -\frac{GM}{r^2}\!\left(1 - \frac{GM}{c^2 r}\right)\hat{\mathbf{e}}_r.\;}\quad\text{(paper Eq. h4, Mercury-side)}$$

The Sun-side equation `a_s = (Gm/r²)(1 − Gm/(c²r))·ê_r` follows by m ↔ M symmetry.

**Wolfram MCP check:**

```wolfram
ClearAll[r, cc, mm, GG, MM, potV, accel]; potV[r_] := -GG MM mm/r; gradV = D[potV[r], r]; accel = -gradV/mm (1 + potV[r]/(mm cc^2)); paperh4 = -(GG MM/r^2) (1 - GG MM/(cc^2 r)); Print["a (derived) = ", FullSimplify[accel]]; Print["FullSimplify[a - paper_h4] = ", FullSimplify[accel - paperh4]]
```

**Result:** `a (derived) = GG MM (GG MM - cc² r)/(cc² r³)`; `FullSimplify[a - paper_h4] = 0`. ✅

## 4. Physical interpretation

The `(1 − GM/(c²r))` factor is *less than 1* for finite `r`. Two consequences:

1. **Gravitational attraction is *weaker* than Newtonian** by a fraction `GM/(c²r)`. For Mercury at the Sun: `GM/(c²r) ≈ 2.55 × 10⁻⁸` — small, but present.

2. **At `r = r₀ ≡ GM/c²`** (the Sun's gravitational radius, ≈ `1.48 km` — half its Schwarzschild radius) the framework predicts zero gravity, and for `r < r₀` gravity *repels*. This is the paper's "classical principle of impenetrability"; the framework predicts no singularity at `r = 0`.

The first consequence is the load-bearing input to the perihelion calculation; the second is structurally interesting but not Mercury-relevant.

**Sign-of-correction note.** Standard GR's perihelion advance comes from an extra *attractive* `1/r³` term in the effective potential (the `3GM·L²/(c²r³)` term from the Schwarzschild geodesic). The framework's Eq. (h4) gives the *opposite* sign of relativistic correction — the `(1 − GM/(c²r))` factor *reduces* the `1/r²` attraction. This sign difference becomes load-bearing in [`03_numerical_predictions.md` §3](03_numerical_predictions.md) where the framework's full dual under-predicts the perihelion advance and (under the natural interpretation) predicts the wrong sign.

<!-- TODO: human reviews and fills in — confirms the "wrong sign of relativistic correction" reading is the correct interpretation of the (1 − GM/(c²r)) factor's mechanism. The paper does not explicitly compare its prediction's sign convention to GR's perihelion-advance sign; the comparison is the campaign's contribution. -->

## 5. Verdict

✅ Paper's algebra for Eqs. (2.4), (h3), (h4) reproduces under Wolfram-MCP verification. The load-bearing input to the subsequent orbital-dynamics calculation is the modified gravitational acceleration `a = −(GM/r²)(1 − GM/(c²r))·ê_r`.
