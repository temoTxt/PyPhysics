# Ch. 1 — Introduction to Electrostatics

This chapter contains Jackson canonical problems worked in the proper-time reformulation alongside their classical CGS and SI solutions. Per [§4 of the plan](../../../.dev/tasks/42-electromagnetism-jackson-proper-time.md#4-unit-system-handling), Ch. 1 of Jackson is presented in the three-system regime: both CGS and SI versions of each problem appear, followed by the proper-time reformulation.

The Ch. 1 problems below are part of PR 0 (the fluency-warm-up prequel — see [§7.3 of the plan](../../../.dev/tasks/42-electromagnetism-jackson-proper-time.md#73-pr-a-prequel-adopted)). They share the structural feature that the source charges are at rest in the observer's frame, so `u = 0` and `b² = c² + u² = c²`. The proper-time formulation reduces to the classical formulation identically in this limit; the work of each problem is in the CGS↔SI translation and in the documentation that the reduction is indeed identical.

## Problems in this chapter

| Problem | Status | Role |
|---|---|---|
| [Problem J3e-P1.5 — Electrostatic self-energy of a uniformly charged sphere](#problem-j3e-p15--electrostatic-self-energy-of-a-uniformly-charged-sphere) | drafted (PR 0) | fluency-builder |

---

### Problem J3e-P1.5 — Electrostatic self-energy of a uniformly charged sphere

**Selection provenance** (Crocco §5 substantive-AI note):
- *Chosen because:* it is the simplest electrostatic problem with a non-trivial volumetric integral, and its solution in either unit system is a textbook result. As the first per-problem document in the campaign, it serves as a canary for the [`_template_problem.md`](../_template_problem.md) structure, the [`VOICE_MATCH_GILL.md`](../../Tooling/VOICE_MATCH_GILL.md) read-aloud test, and the Mathematica MCP workflow.
- *Alternatives considered:* J3e-P1.1 (Gauss's law for a closed surface — too short to exercise the template) and J3e-P1.7 (energy of charge distributions in general — less concrete than the specific charged-sphere result).
- *Role in this PR:* fluency-builder. The proper-time formulation reduces to the classical formulation identically; the work of the problem is in (a)↔(b) and in confirming the reduction in (c).

<!-- TODO: human reviews and fills in — confirms this Selection-provenance reasoning before the problem is rolled into PR 0's commit -->

**Source:** Jackson, *Classical Electrodynamics*, 3e Problem 1.5 (and 2e Problem 1.5, equivalent at this number). *Paraphrased; consult the textbook for the precise statement.*

**Paraphrased statement:** Compute the total electrostatic energy stored in a solid sphere of radius `R` carrying a uniform volume charge density, with total charge `Q`. Express the result in Gaussian units and, separately, in SI units.

**Setup:** A solid sphere of radius `R` carries a uniform volume charge density `ρ = 3Q/(4πR³)`. The sphere is at rest in the observer's frame, so the source velocity satisfies `w = 0`, and consequently the proper-time velocity is `u = 0` and `b² = c² + u² = c²`. There are no currents and no time-varying fields; the problem is purely electrostatic.

#### (a) Classical solution — Gaussian (CGS)

In Gaussian units, Gauss's law applied to a spherical surface of radius `r` centred on the sphere gives the electric field

$$
\mathbf{E}(r) = \begin{cases} \dfrac{Q\,r}{R^{3}}\hat r & r \le R, \\[1ex] \dfrac{Q}{r^{2}}\hat r & r \ge R. \end{cases}
$$

The electrostatic energy density in Gaussian units is `u_E = E²/(8π)`. The total energy is the volume integral of `u_E` over all space; in the spherical symmetry of the present configuration, it decomposes into an inside-the-sphere contribution and an outside-the-sphere contribution,

$$
W = \int_{0}^{R}\!\!\frac{1}{8\pi}\!\left(\frac{Q\,r}{R^{3}}\right)^{\!2}\! 4\pi r^{2}\, dr + \int_{R}^{\infty}\!\!\frac{1}{8\pi}\!\left(\frac{Q}{r^{2}}\right)^{\!2}\! 4\pi r^{2}\, dr.
$$

Evaluating the two integrals separately,

$$
W = \frac{Q^{2}}{2 R^{6}}\!\int_{0}^{R}\! r^{4}\, dr + \frac{Q^{2}}{2}\!\int_{R}^{\infty}\!\frac{dr}{r^{2}} = \frac{Q^{2}}{10 R} + \frac{Q^{2}}{2 R} = \frac{3}{5}\,\frac{Q^{2}}{R}.
$$

**Mathematica check** (Wolfram MCP, 2026-05-24):

```mathematica
ClearAll[r, capR, capQ];
wInside = Integrate[(capQ^2 r^2/capR^6)/(8 Pi) 4 Pi r^2, {r, 0, capR}, Assumptions -> capR > 0];
wOutside = Integrate[(capQ^2/r^4)/(8 Pi) 4 Pi r^2, {r, capR, Infinity}, Assumptions -> capR > 0];
FullSimplify[wInside + wOutside]
(* Result: (3 capQ^2)/(5 capR)  ✅ *)
```

(The `capQ` and `capR` variable names are stand-ins for `Q` and `R`; we use the longer forms to avoid any ambiguity with single-letter Wolfram protected symbols, in keeping with the Mathematica MCP discipline recorded in [`CLAUDE.md`](../../../CLAUDE.md#how-equation-verification-works).)

We observe that the inside-the-sphere contribution accounts for one-sixth of the total and the outside-the-sphere contribution for five-sixths. This is the well-known classical result for the electrostatic self-energy of a uniformly charged solid sphere, recorded for example in Jackson 3e Eq. (1.55) and in equivalent form in 2e.

#### (b) Classical solution — SI

In SI units the electrostatic energy density is `u_E = (ε₀/2) E²`, and the electric field carries an additional `1/(4πε₀)` factor:

$$
\mathbf{E}(r) = \begin{cases} \dfrac{Q\,r}{4\pi\varepsilon_{0} R^{3}}\hat r & r \le R, \\[1ex] \dfrac{Q}{4\pi\varepsilon_{0} r^{2}}\hat r & r \ge R. \end{cases}
$$

Substituting into `W = ∫ (ε₀/2) E² dV`, the volume integrals are structurally identical to those of (a), and one obtains

$$
W = \frac{3}{5}\,\frac{Q^{2}}{4\pi\varepsilon_{0} R}.
$$

The translation between (a) and (b) consists of one substitution at the level of the answer: `Q²_{\text{Gaussian}}/R \to Q²_{\text{SI}}/(4\pi\varepsilon_{0} R)`. The geometric `3/5` factor is unchanged by the choice of unit system, as one expects, since it arises from the volume integrals alone and is therefore insensitive to whether the field is expressed in Gaussian or SI form.

#### (c) Proper-time reformulation

The source charges are at rest in the observer's frame, so `w = 0`. It follows from the velocity duality of Eq. (1) of [[Two_Mathematically_Equivalent_Versions_of_Maxwells_Equations]] that `u = 0`, and from `b² = c² + u²` that `b = c` identically throughout the configuration. Each of the substitution rules summarised in [`_proper_time_cheatsheet.md`](../_proper_time_cheatsheet.md) reduces in this limit as follows:

- `c → b` is the identity, because `b = c` everywhere.
- `(1/c)\, \partial_t = (1/b)\, \partial_\tau` is the trivial statement that two zero quantities are equal; the configuration is static, so no time derivative appears in either formulation.
- The dissipative chain-rule term `\partial_\tau(1/b) = -(u\cdot a)/b^3` vanishes by `u = 0`.
- The current-density rescaling `J → (b/c) J` is the identity, because there is no current and `b/c = 1`.

Each substitution rule that the proper-time formulation could introduce is either trivially the identity or annihilated by the `u = 0` condition. It is clear that the proper-time prediction for the total electrostatic energy is identical to the classical result of (a):

$$
W_{\text{proper-time}} = \frac{3}{5}\,\frac{Q^{2}}{R}.
$$

The two formulations are *mathematically equivalent and physically equivalent* in this case — a stronger statement than the general "mathematically equivalent but not physically equivalent" relation of the Gill–Zachary paper, because the configuration's static character closes the gap between the two clocks. The local clock of the source agrees with the observer's clock when the source is at rest, and so the information that the local clock would otherwise encode (per Eq. (4) of the Maxwell paper) is identically zero here.

**Comparison:**

| Quantity | Classical (CGS) | Classical (SI) | Proper-time |
|---|---|---|---|
| `W` (total energy) | `(3/5)\,Q^{2}/R` | `(3/5)\,Q^{2}/(4\pi\varepsilon_{0} R)` | identical to CGS |

**Does the proper-time answer differ from a pure `c → b` redressing?** ✅ no. With `u = 0` we have `b = c`, so the `c → b` redressing is the identity, and the answer is unchanged.

**Verdict:** ✅ all three solutions consistent. The proper-time formulation contains no machinery that engages in the static case; it reproduces the classical answer exactly.

**Notes for author review:** none. This problem is the cleanest possible reduction of the proper-time formulation to the classical one, and serves as the campaign's first canary — for the [`_template_problem.md`](../_template_problem.md) structure, for the [`VOICE_MATCH_GILL.md`](../../Tooling/VOICE_MATCH_GILL.md) read-aloud test, and for the Mathematica MCP workflow. Any drift in the template or in the voice surfaces here first.

**Companion notebook:** [`Roadmapping/Mathematica_Notebooks/Electromagnetism/JacksonCh01_P1_5.wl`](../../Mathematica_Notebooks/Electromagnetism/JacksonCh01_P1_5.wl) — runnable independent of the Mathematica MCP.

---

*Note on chapter scope.* The remaining Ch. 1 problems planned for PR 0 reside in Ch. 2 (image-method problems), Ch. 5 (magnetostatic problems), and the optional fourth problem also in Ch. 5. Per the campaign's per-chapter document convention, the image-method problem appears under [`Ch02_Boundary_Value_Problems_Electrostatics_I.md`](Ch02_Boundary_Value_Problems_Electrostatics_I.md) rather than here.
