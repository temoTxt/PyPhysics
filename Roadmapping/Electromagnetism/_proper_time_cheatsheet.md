# Proper-time substitution rules — cheat sheet

The proper-time reformulation of Maxwell's equations, established by Gill and Zachary in [[Two_Mathematically_Equivalent_Versions_of_Maxwells_Equations]] and verified at [`Roadmapping/Equation_Verification/Two_Mathematically_Equivalent_Versions_of_Maxwells_Equations.md`](../Equation_Verification/Two_Mathematically_Equivalent_Versions_of_Maxwells_Equations.md), is governed by a finite list of substitution rules. Every per-problem reformulation reduces to applying them. This document is the one-page reference; the per-problem documents under [`Jackson/`](Jackson/) cite it rather than restating the rules.

## 1. Definitions

Let `w = dx/dt` denote the velocity of a particle as measured by an observer's clock, and let `u = dx/dτ` denote the velocity as measured by the particle's local clock (its proper time). The effective speed of light in the proper-time frame is

```
b² = c² + u²,
```

so `b ≥ c` always, and `b → c` when the particle is at rest in the observer's frame. The proper-time acceleration is `a = du/dτ`. It follows from `b² = c² + u²` that

```
2 b ḃ = 2 u · a,    so    ḃ = (u · a) / b,
```

a result that arises whenever a τ-derivative meets `b`'s τ-dependence.

## 2. Substitution rules

| Rule | Form | Origin | Application |
|---|---|---|---|
| Velocity duality | `w / c = u / b` | Eq. (1) of the Maxwell paper | Replace `w/c` by `u/b` in any kinematic expression. |
| Time-derivative duality | `(1/c) ∂_t = (1/b) ∂_τ` | Eq. (2) | Replace any classical time derivative by the proper-time analogue. |
| Gauss for **E**, divergenceless **B** | `∇·E = 4πρ`, `∇·B = 0` | Eq. (3′) | Unchanged from standard Maxwell. |
| Faraday | `∇×E = −(1/b) ∂_τ B` | Eq. (3′) | The classical `1/c` becomes `1/b`; the time variable becomes `τ`. |
| Ampère–Maxwell | `∇×B = (1/b) ∂_τ E + (4π/b) ρ u` | Eq. (3′) | The current density picks up an explicit `(b/c)` factor when one identifies `J = (b/c) ρ w = ρ u`. |
| Bare chain rule | `∂_τ(1/b) = −(u·a)/b³` | Derived from Eq. (4) | Applies whenever a τ-derivative acts through a `1/b`. |
| Dissipative coefficient (wave equation) | `−(u·a)/b⁴ · ∂_τ` | Eq. (4) | The same `−(u·a)/b³` factor, multiplied by an additional `1/b` once it has been collected as a coefficient of `∂_τ` in a second-order wave equation. Algebraically `(u·a)/b⁴ = ḃ/b³`; the two forms are equivalent. |
| Boost of `b` | `b' = γ(v) [b − u · v / c]` | Eq. (11) | Lorentz boost in the proper-time formulation. |
| Boost of current density | `J' = J + (γ − 1)(J · v) v / v² − γ (b/c) ρ v` | Eq. (12) | Spatial part of the 4-current `(b ρ, J)` boost. |
| Modified Lorentz force | `F = e [E + (u/b) × B] + ...` | Eq. (18) | Includes the proper-time correction; the trailing terms are the dissipative back-reaction. |

The same result expressed in two forms (the bare chain rule `−(u·a)/b³` and the wave-equation coefficient `−(u·a)/b⁴ · ∂_τ`) is one of the most common sources of confusion in early proper-time derivations. The two are equivalent under `ḃ = (u·a)/b`; pick whichever form is natural for the equation at hand.

## 3. The decision tree for a Jackson problem

In order to apply the substitution rules to a Jackson problem, we follow four steps:

1. Identify the classical equation(s) the problem uses.
2. Determine whether the equation is a Maxwell equation, a force law, or a derived quantity (potential, energy, momentum, etc.).
3. Apply the substitution rules to the underlying Maxwell / force law, and then re-derive the derived quantity *from scratch* in the proper-time frame. We do not simply substitute `c → b` into the classical answer — derived quantities pick up additional terms from `b(τ)`'s τ-dependence, as the structure of Maxwell-paper Eq. (4) demonstrates.
4. Where the proper-time answer differs from the `c → b` redressing of the classical answer, document the extra term and its physical interpretation (radiation reaction, longitudinal radiation component, etc.).

## 4. The Eq. 24 branching rule

[`FINDINGS_for_author_review.md`](../Equation_Verification/FINDINGS_for_author_review.md) records that Maxwell-paper Eq. (24) — the proper-time Pauli-form Hamiltonian — carries a flagged finding: it appears to be missing a factor of `c` in the spin-magnetic-moment term and the entire `V²/(2mc²)` term. The finding is unresolved as of campaign start.

For any Jackson problem whose proper-time derivation invokes Eq. (24) — likely a subset of Ch. 12 (relativistic dynamics) and any spin-orbit / hyperfine problem — the per-problem document carries a **branched treatment** per [§5.1 of the plan](../../.dev/tasks/42-electromagnetism-jackson-proper-time.md#51-branched-treatment-for-eq-24-touching-problems). The proper-time section is worked twice (`(c)` as-published and `(c′)` with the proposed correction), and the verdict is recorded per branch.

A problem touches Eq. (24) if its derivation needs the Pauli-form Hamiltonian

```
H = (mc² + p²/2m) − eϕ + (e/c) A · p + e ℏ Σ · B / (2m) + V² / (2mc²)
```

or its non-relativistic reduction. Most pure-EM Jackson problems (Chs. 1–11, 14, 16) do not touch Eq. (24); they operate at the field-equation level (Eqs. (3′), (4), (7)), not the Hamiltonian level.

## 4a. Honest open questions (flagged here so per-problem documents can cite back)

The following are *not* substitution rules but open questions about the framework that the per-problem documents have encountered and not resolved. They are gathered here so that any per-problem document can flag dependence on them via a single cross-reference rather than restating the issue:

- **Covariance of `\mathbf{u}\cdot\mathbf{a}`.** The dissipative coefficient `(\mathbf{u}\cdot\mathbf{a})/b^{4}` is a 3-vector dot product between `\mathbf{u} = d\mathbf{x}/d\tau` and `\mathbf{a} = d\mathbf{u}/d\tau`. It is *not* the 4-vector contraction `u^{\mu}a_{\mu}` of standard SR (which vanishes identically by 4-velocity / 4-acceleration orthogonality). How `\mathbf{u}\cdot\mathbf{a}` transforms under the proper-time boost of Eq. (11) of the Maxwell paper has **not been derived in the campaign**. Per-problem documents that rely on `\mathbf{u}\cdot\mathbf{a}` having a frame-independent operational meaning (notably [`Ch14_P14.2`](Jackson/Ch14_Radiation_by_Moving_Charges.md#problem-j3e-p142--li%C3%A9nard-wiechert-fields-with-the-proper-time-third-term) and [`Ch16_P16.5`](Jackson/Ch16_Radiation_Damping.md#problem-j3e-p165--proper-time-rr-prediction-for-the-colepoder-geometry)) carry an inline reminder pointing back to this entry.

- **Far-field `\hat k`-decomposition of the Eq. (7) third term.** Whether the third term `e(\mathbf{u}\cdot\mathbf{a})\mathbf{r}\times(\mathbf{u}\times\mathbf{r})/(b^{4}s^{3})` survives the radiation-zone limit `r\to\infty` as a true outgoing radiation component (consistent with `\nabla\!\cdot\!\mathbf{E} = 0` in vacuum), or instead reorganises into a near-field plus a transverse-radiation correction, has **not been computed in this PR**. The associated `dP/d\Omega` integration over a sphere at `r\to\infty` is the cleanest test and is referred forward to the issue [#43](https://github.com/temoTxt/PyPhysics/issues/43) numerical setup.

- **Explicit proper-time equation of motion under radiation reaction.** The proper-time radiation-reaction force `\mathbf{F}_{\text{RR, PT}}(\mathbf{u}, \mathbf{a})` — needed for any explicit step-force IVP, and ultimately for the dissolution-of-pathologies claim of [`Ch16_P16.1`](Jackson/Ch16_Radiation_Damping.md#problem-j3e-p161--abrahamlorentz-radiation-reaction-and-the-proper-time-dissolution-claim) — has not been derived from first principles in the campaign. The dissolution argument is currently a structural / operator-order statement, not a derived result. See P16.1's "Honest gap" subsection for the minimum-bar derivation steps.

These three open questions are the load-bearing follow-on items that turn the PR D / PR E headline content from "structural prediction of the framework" into "derived consequence of the framework" — the further work that issue #43 and any successor issue will need to do.

## 5. A pedagogical note on Eqs. (3′), (4), and (7)

It is useful to keep in mind that the three equations cited most often in this campaign — (3′) for the proper-time Maxwell equations, (4) for the wave equation with its dissipative term, and (7) for the modified Liénard–Wiechert fields — describe the same physics at three levels of derivation. Eq. (3′) is the underlying field theory. Eq. (4) is what one obtains by taking the curl of Eq. (3′) and using standard vector identities; the dissipative term `−(u·a)/b⁴ · ∂_τ` arises from `∂_τ(1/b)` acting on the inhomogeneous term during the curl-of-curl manipulation. Eq. (7) is the explicit field of a point charge, computed by solving Eq. (4) with the appropriate retarded conditions; the third term, proportional to `(u·a)/b⁴`, is the same dissipative effect now expressed in the field's spatial structure rather than as a wave-equation coefficient.

The three levels are mathematically equivalent, but each is the natural level at which to attack a particular class of Jackson problem. Field-configuration problems (Chs. 1–5) live at Eq. (3′). Wave-equation problems (Chs. 7–10) live at Eq. (4). Radiation problems (Chs. 14, 16) live at Eq. (7). When in doubt about which level to start from, identify the Jackson chapter and the corresponding equation should follow.
