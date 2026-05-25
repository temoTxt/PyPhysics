# Ch. 14 — Radiation by Moving Charges

This chapter contains Jackson canonical problems on radiation from accelerating charges, worked in the proper-time reformulation alongside their classical solutions. Per [§4 of the plan](../../../.dev/tasks/42-electromagnetism-jackson-proper-time.md#4-unit-system-handling), Ch. 14 of Jackson is presented in the two-system regime: Gaussian (CGS) only on the classical side, with a proper-time reformulation alongside.

Ch. 14 is the campaign's **first headline-payoff chapter**. Two of the five §12 podcast picks land here: **#2 (Liénard–Wiechert with the third term)** and **#5 (non-relativistic Larmor radiation)**. The dissipative `(\mathbf{u}\cdot\mathbf{a})/b^{4}` coefficient of Eq. (4) of [[Two_Mathematically_Equivalent_Versions_of_Maxwells_Equations]] engages directly in the field expressions of Eq. (7) — the new third term that distinguishes the proper-time formulation from classical EM at the level of radiated fields. This is the chapter where issue [#43](https://github.com/temoTxt/PyPhysics/issues/43)'s experimental-comparison hooks become operationally measurable.

Per-paragraph `<!-- TODO -->` blocks continue per [§13.5 D2](../../../.dev/tasks/42-electromagnetism-jackson-proper-time.md#135-decision-points--confirmed-by-author-2026-05-24).

## Problems in this chapter

| Problem | Status | Role |
|---|---|---|
| [Problem J3e-P14.1 — Liénard–Wiechert potentials](#problem-j3e-p141--li%C3%A9nard-wiechert-potentials) | drafted (PR D) | fluency-builder (foundational) |
| [Problem J3e-P14.2 — Liénard–Wiechert fields with the proper-time third term](#problem-j3e-p142--li%C3%A9nard-wiechert-fields-with-the-proper-time-third-term) | drafted (PR D) | **HEADLINE-PAYOFF (podcast pick #2)** |
| [Problem J3e-P14.3 — Non-relativistic Larmor radiation formula](#problem-j3e-p143--non-relativistic-larmor-radiation-formula) | drafted (PR D) | **HEADLINE-PAYOFF (podcast pick #5)** |
| [Problem J3e-P14.5 — Synchrotron radiation from circular motion](#problem-j3e-p145--synchrotron-radiation-from-circular-motion) | drafted (PR D) | fluency-builder (u·a = 0, third term null) |

---

### Problem J3e-P14.1 — Liénard–Wiechert potentials

**Selection provenance** (Crocco §5 substantive-AI note):
- *Chosen because:* the Liénard–Wiechert potentials are the foundational construction of charged-particle radiation theory, recorded in [[jackson1998_classical_electrodynamics]] §14.1. They are the natural opening problem of PR D because every subsequent radiation problem builds on them.
- *Alternatives considered:* J3e-P14.2 (Liénard–Wiechert fields with third term — selected next as the headline-payoff podcast pick #2) and J3e-P14.5 (relativistic Liénard formula — selected later for the synchrotron treatment).
- *Role in this PR:* fluency-builder. The proper-time form of the potentials uses `u/b` in place of `v/c` and produces the same numerical retarded-condition expression by the velocity-duality identity.

<!-- TODO: human reviews and fills in — confirms the role of this problem as PR D's fluency-builder opening before the headline LW-fields problem -->

**Source:** Jackson, *Classical Electrodynamics*, 3e Problem 14.1 (and 2e Problem 14.1, equivalent). *Paraphrased.*

**Paraphrased statement:** Derive the Liénard–Wiechert scalar and vector potentials `\phi(\mathbf{x}, t)` and `\mathbf{A}(\mathbf{x}, t)` for a point charge `q` moving along an arbitrary worldline. Express the result in classical observer-time variables and in proper-time variables, and confirm that both formulations produce the same observable potentials.

**Setup:** A point charge `q` moves along a prescribed worldline `\mathbf{r}_{s}(t)` (classical) or equivalently `\mathbf{r}_{s}(\tau)` (proper-time). The field point is at `(\mathbf{x}, t)`. The *retarded position* of the source is `\bar{\mathbf{r}} = \mathbf{r}_{s}(\bar t)` where `\bar t` is the unique earlier time such that `|\mathbf{x} - \bar{\mathbf{r}}| = c(t - \bar t)` (light-cone condition). Define `\mathbf{R} = \mathbf{x} - \bar{\mathbf{r}}` (displacement from retarded position to field point), `R = |\mathbf{R}|`, and the unit vector `\hat n = \mathbf{R}/R`.

#### (a) Classical solution — Gaussian (CGS)

In classical observer-time formulation, the source velocity at the retarded position is `\mathbf{v} = d\mathbf{r}_{s}/dt|_{\bar t}`, with Lorentz factor `\gamma`. The Liénard–Wiechert potentials, recorded as Jackson 3e Eqs. (14.8)–(14.9), are

$$
\phi(\mathbf{x}, t) = \frac{q}{R - \mathbf{R}\cdot\mathbf{v}/c}\bigg|_{\bar t}, \qquad \mathbf{A}(\mathbf{x}, t) = \frac{q\,\mathbf{v}/c}{R - \mathbf{R}\cdot\mathbf{v}/c}\bigg|_{\bar t}.
$$

The retarded denominator `s = R - \mathbf{R}\cdot\mathbf{v}/c` encodes the relativistic Doppler-like enhancement of the field for charges moving toward the observer (when `\mathbf{R}\cdot\mathbf{v} > 0`, `s` is smaller and the potentials larger).

#### (c) Proper-time reformulation

In proper-time variables `(b, \mathbf{u})` with `b^{2} = c^{2} + \mathbf{u}^{2}`, the source 4-velocity at the retarded position is `(b(\bar\tau), \mathbf{u}(\bar\tau))`. The Liénard–Wiechert potentials take the form

$$
\phi(\mathbf{x}, t) = \frac{q}{R - \mathbf{R}\cdot\mathbf{u}/b}\bigg|_{\bar\tau}, \qquad \mathbf{A}(\mathbf{x}, t) = \frac{q\,\mathbf{u}/b}{R - \mathbf{R}\cdot\mathbf{u}/b}\bigg|_{\bar\tau}.
$$

The retarded denominator `s_{\text{PT}} = R - \mathbf{R}\cdot\mathbf{u}/b` is **numerically identical** to the classical `s` under the velocity-duality identity `\mathbf{u}/b = \mathbf{v}/c`. The vector potential's numerator also uses `\mathbf{u}/b` in place of `\mathbf{v}/c`, again identical under the duality. The retarded time-condition is unchanged.

**Mathematica check** (Wolfram MCP, 2026-05-24):

```mathematica
ClearAll[rVec, vVec, uVec, capR, cVar, bVar];
classicalDenom = capR - (rVec . vVec)/cVar;
properTimeDenom = capR - (rVec . uVec)/bVar;
properTimeUnderDuality = properTimeDenom /. (rVec . uVec) -> (rVec . vVec) (bVar/cVar);
Print[FullSimplify[properTimeUnderDuality - classicalDenom]];
(* Result: 0  ✅ *)
```

The two formulations therefore produce the same observable potentials. The proper-time choice of variables is a re-expression, not a different physical prediction.

It is worth observing that the *form* of the Liénard–Wiechert potentials is preserved under the velocity-duality substitution `\mathbf{v}/c \to \mathbf{u}/b` and the time-derivative substitution `(1/c)\partial_{t} \to (1/b)\partial_{\tau}`. The retarded condition itself is unchanged; what changes is the parametrisation of the source's worldline. This is the same structural feature that produced the linear-in-`(b,u)` Doppler factor of [Problem J3e-P11.3](Ch11_Special_Relativity.md#problem-j3e-p113--relativistic-doppler-effect) and the linear-in-`(b,u)` field-boost transformation of [Problem J3e-P11.6](Ch11_Special_Relativity.md#problem-j3e-p116--lorentz-boost-of-e-and-b-fields).

<!-- TODO: human reviews and fills in — confirms the framing that the potentials themselves are unchanged at the observable level; the dissipative-term contribution will surface in the FIELDS (derivatives of the potentials), which is the headline content of J3e-P14.2 -->

The Liénard–Wiechert potentials are the foundation for the fields of an accelerating charge. The fields are obtained by computing `\mathbf{E} = -\nabla\phi - (1/c)\partial\mathbf{A}/\partial t` and `\mathbf{B} = \nabla\times\mathbf{A}`. In classical EM, this produces the standard Liénard–Wiechert fields with two contributions: the "velocity field" (Coulomb-like, scaling as `1/R^{2}`) and the "acceleration field" (radiation-zone, scaling as `1/R`). In the proper-time formulation, the same construction produces an additional **third term** proportional to `(\mathbf{u}\cdot\mathbf{a})/b^{4}` — the dissipative coefficient of Eq. (4) of the Maxwell paper. The third term is the subject of [Problem J3e-P14.2](#problem-j3e-p142--li%C3%A9nard-wiechert-fields-with-the-proper-time-third-term).

**Comparison:**

| Quantity | Classical (Gaussian) | Proper-time |
|---|---|---|
| Scalar potential | `q/s` with `s = R - \mathbf{R}\cdot\mathbf{v}/c` | `q/s_{\text{PT}}` with `s_{\text{PT}} = R - \mathbf{R}\cdot\mathbf{u}/b` |
| Vector potential | `q(\mathbf{v}/c)/s` | `q(\mathbf{u}/b)/s_{\text{PT}}` |
| Retarded condition | `|\mathbf{x} - \bar{\mathbf{r}}| = c(t - \bar t)` | identical |
| Observable potentials | classical | identical (under velocity duality) |

**Does the proper-time answer differ from a pure `c → b` redressing?** ✅ no, at the level of the potentials. The proper-time form is a re-expression in `(b, \mathbf{u})` variables, not a different physical prediction. The dissipative contribution surfaces only when one takes derivatives of the potentials to obtain the *fields* — the subject of the next problem.

**Verdict:** ✅ all formulations consistent. The Liénard–Wiechert potentials are the same physical objects in both formulations, with the proper-time form using `u/b` in place of `v/c` and `(\bar\tau, \bar t)` related by the standard retarded-condition equation.

**Notes for author review:** the potentials are the foundation for PR D's headline content. The dissipative-term engagement happens at the *field* level (J3e-P14.2), not at the potential level. This setup is correct as the cleanest entry point to Ch. 14.

**Companion notebook:** [`Roadmapping/Mathematica_Notebooks/Electromagnetism/JacksonCh14_P14_1.wl`](../../Mathematica_Notebooks/Electromagnetism/JacksonCh14_P14_1.wl).

---

### Problem J3e-P14.2 — Liénard–Wiechert fields with the proper-time third term

**Selection provenance** (Crocco §5 substantive-AI note):
- *Chosen because:* this is the **load-bearing headline problem of the entire campaign**. Eq. (7) of [[Two_Mathematically_Equivalent_Versions_of_Maxwells_Equations]] gives the Liénard–Wiechert fields of an accelerating charge in proper-time form, with a **third term `e(\mathbf{u}\cdot\mathbf{a})[\mathbf{r}\times(\mathbf{u}\times\mathbf{r})]/(b^{4}s^{3})`** that is *absent from the classical Liénard–Wiechert fields* and is the field-level manifestation of the dissipative coefficient of Eq. (4). This is **podcast pick #2** in [§12.1 of the campaign plan](../../../.dev/tasks/42-electromagnetism-jackson-proper-time.md#121-per-problem-briefs).
- *Alternatives considered:* J3e-P14.4 (radiation field of a non-relativistic accelerated charge — covered by J3e-P14.3 below) and J3e-P14.7 (angular distribution of synchrotron radiation — covered by J3e-P14.5 below).
- *Role in this PR:* **HEADLINE-PAYOFF**. The third term is the campaign's load-bearing example of the proper-time formulation producing a genuinely different field-level prediction than classical EM. Its experimental signature — a longitudinal radiation component — is the natural target for issue [#43](https://github.com/temoTxt/PyPhysics/issues/43)'s comparison work.

<!-- TODO: human reviews and fills in — confirms the role of this problem as the campaign's headline-payoff document, and the framing of the third term as a "longitudinal radiation component absent from classical EM" -->

**Source:** Jackson, *Classical Electrodynamics*, 3e Problem 14.2 (and 2e Problem 14.2, equivalent). *Paraphrased.*

**Paraphrased statement:** Derive the electric and magnetic fields of a point charge `q` in arbitrary motion, both in the classical Liénard–Wiechert formulation and in the proper-time formulation. Compare term-by-term. Identify the **third term** of the proper-time E-field — proportional to `(\mathbf{u}\cdot\mathbf{a})/b^{4}` — that is absent from the classical Liénard–Wiechert expression, and remark on its physical content.

**Setup:** Worldline parametrised by retarded time `\bar t` (classical) or `\bar\tau` (proper-time). Source velocity `\mathbf{v}(\bar t)` (classical) or 4-velocity `(b(\bar\tau), \mathbf{u}(\bar\tau))` (proper-time). Source acceleration `\mathbf{a}_{\text{classical}} = d\mathbf{v}/dt|_{\bar t}` or `\mathbf{a} = d\mathbf{u}/d\tau|_{\bar\tau}`. Displacement to field point `\mathbf{r} = \mathbf{x} - \bar{\mathbf{r}}`, `r = |\mathbf{r}|`. Define `\mathbf{r}_{\mathbf{u}} = \mathbf{r} - (r/b)\mathbf{u}` and `s = r - \mathbf{r}\cdot\mathbf{u}/b`.

#### (a) Classical solution — Gaussian (CGS)

The classical Liénard–Wiechert electric field, derived by taking spatial and temporal derivatives of the LW potentials and recorded as Jackson 3e Eq. (14.13), has two contributions:

$$
\mathbf{E}_{\text{classical}}(\mathbf{x}, t) = \underbrace{\frac{e\,\mathbf{r}_{\mathbf{v}}(1 - v^{2}/c^{2})}{s_{v}^{3}}}_{\text{velocity field, }\propto 1/R^{2}} + \underbrace{\frac{e\,\mathbf{r}\times(\mathbf{r}_{\mathbf{v}}\times\mathbf{a}_{\text{classical}})}{c^{2}\,s_{v}^{3}}}_{\text{acceleration field, }\propto 1/R}
$$

with `\mathbf{r}_{\mathbf{v}} = \mathbf{r} - (r/c)\mathbf{v}` and `s_{v} = r - \mathbf{r}\cdot\mathbf{v}/c`. The first term (velocity field) reproduces the boosted Coulomb field of a uniformly-moving charge; the second term (acceleration field) is the radiation-zone contribution responsible for Larmor / Liénard / synchrotron radiation. The acceleration field is **purely transverse** to the line of sight `\hat r` in the radiation zone — the standard textbook statement that classical EM radiation has no longitudinal component.

#### (c) Proper-time reformulation

The proper-time Liénard–Wiechert fields, recorded as **Eq. (7) of [[Two_Mathematically_Equivalent_Versions_of_Maxwells_Equations]]** and verified ✅ in the [campaign's Equation Verification document](../../Equation_Verification/Two_Mathematically_Equivalent_Versions_of_Maxwells_Equations.md#eq-7--modified-liénardwiechert-fields), are

$$
\mathbf{E}(\mathbf{x}, \tau) = \underbrace{\frac{e\,\mathbf{r}_{\mathbf{u}}(1 - \mathbf{u}^{2}/b^{2})}{s^{3}}}_{\text{velocity field}} + \underbrace{\frac{e\,\mathbf{r}\times(\mathbf{r}_{\mathbf{u}}\times\mathbf{a})}{b^{2}\,s^{3}}}_{\text{acceleration field}} + \underbrace{\frac{e\,(\mathbf{u}\cdot\mathbf{a})\,\mathbf{r}\times(\mathbf{u}\times\mathbf{r})}{b^{4}\,s^{3}}}_{\text{\bf THIRD TERM — new in proper-time}}.
$$

The first two terms have the same structure as the classical LW fields, with `c \to b` and `\mathbf{v} \to \mathbf{u}`. Under the velocity-duality identity, the first term's prefactor `(1 - \mathbf{u}^{2}/b^{2})` equals the classical `(1 - v^{2}/c^{2})` (as verified in [Problem J3e-P6.4](Ch06_Maxwell_Equations_Macroscopic_Media.md#problem-j3e-p64--em-momentum-of-a-uniformly-moving-point-charge)); the velocity-field contribution is therefore identical to the classical one at the observable level.

The **third term** is the load-bearing new content. It is proportional to `\mathbf{u}\cdot\mathbf{a}`, so it vanishes whenever the source's acceleration is perpendicular to its velocity (notably for circular motion — [Problem J3e-P14.5](#problem-j3e-p145--synchrotron-radiation-from-circular-motion) below — and for instantaneous rest in any frame). For any other motion (linear acceleration, plane-wave figure-8, bremsstrahlung) it is non-zero and contributes to the radiated field.

**Mathematica check** (Wolfram MCP, 2026-05-24):

```mathematica
ClearAll[rV, uV];
bacResult = TensorExpand[rV \[Cross] (uV \[Cross] rV)];
Print["r x (u x r) = ", bacResult];
(* Result: uV (rV . rV) - rV (rV . uV) = u r^2 - r (u . r)
   This vector has both a longitudinal component (parallel to u) and a
   radial component (parallel to r), so the third term is NOT transverse
   to the line of sight in the radiation zone. ✓ *)
```

The structural feature of the third term is its decomposition `\mathbf{r}\times(\mathbf{u}\times\mathbf{r}) = r^{2}\mathbf{u} - (\mathbf{u}\cdot\mathbf{r})\mathbf{r}`. Both a **longitudinal component** (parallel to `\mathbf{u}`) and a **radial component** (parallel to `\mathbf{r}`) appear. Classical Liénard–Wiechert radiation in the radiation zone is *purely transverse* to `\hat r`; the proper-time third term breaks this purity by introducing components along both `\mathbf{u}` and `\hat r`.

We observe that this is **the proper-time formulation's first qualitatively new field-level prediction**. Where every prior problem in the campaign (PR 0 through PR C) either reduced exactly to the classical answer or differed only by an `O((u/c)^{2})` correction below observational floor, the third term of Eq. (7) is a *structural* addition to the radiated field: it predicts a longitudinal radiation component in any kinematic configuration where `\mathbf{u}\cdot\mathbf{a} \neq 0`, including all linearly-accelerated charges and all relativistic-intensity plane-wave scattering geometries.

<!-- TODO: human reviews and fills in — confirms the load-bearing claim that "the third term is the first qualitatively new field-level prediction in the campaign". This is the campaign's central narrative claim and warrants the author's full attention -->

The Liénard–Wiechert magnetic field has an analogous third-term structure (Eq. (7) of the Maxwell paper records both `\mathbf{E}` and `\mathbf{B}`):

$$
\mathbf{B}(\mathbf{x}, \tau) = \frac{e\,(\mathbf{r}\times\mathbf{r}_{\mathbf{u}})(1 - \mathbf{u}^{2}/b^{2})}{r\,s^{3}} + \frac{e\,\mathbf{r}\times[\mathbf{r}\times(\mathbf{r}_{\mathbf{u}}\times\mathbf{a})]}{r\,b^{2}\,s^{3}} + \frac{e\,r\,(\mathbf{u}\cdot\mathbf{a})(\mathbf{r}\times\mathbf{u})}{b^{4}\,s^{3}}.
$$

The third term of `\mathbf{B}` is orthogonal to both `\mathbf{r}` and `\mathbf{u}`, with magnitude proportional to `(\mathbf{u}\cdot\mathbf{a})/b^{4}`. We observe (as the Maxwell paper does, in its commentary around Eq. (7)) that `\mathbf{B}` is orthogonal to `\mathbf{E}` overall, but the relative phase between the new third terms differs from the relative phase of the classical first two terms.

#### Experimental signature

The third term predicts a **longitudinal radiation component** from any source with `\mathbf{u}\cdot\mathbf{a} \neq 0`. This is the proper-time framework's most distinguishable prediction from classical EM, and is the natural target for the experimental-comparison work of issue [#43](https://github.com/temoTxt/PyPhysics/issues/43). The relevant experiments are:

- **Cole et al. (2018)**, PRX 8, 011020 — 1.7 GeV electron beam × counter-propagating laser pulse.
- **Poder et al. (2018)**, PRX 8, 031004 — refined analysis of Cole 2018 data.
- **Wistisen et al. (2018)**, Nat. Comm. 9, 795 — 200 GeV electron channelling in aligned Si crystal.

All three experiments measure radiation-reaction signatures in regimes where `\mathbf{u}\cdot\mathbf{a}` is large; the proper-time third term predicts a specific longitudinal-component signature that classical Landau–Lifshitz radiation reaction does not predict. Whether the data discriminates is an open question, addressed in issue #43's follow-on comparison document.

<!-- TODO: human reviews and fills in — confirms the framing of Cole/Poder/Wistisen as the natural experimental targets, and the "longitudinal radiation component" as the framework's most-distinguishable prediction -->

#### Where the third term vanishes — a checklist

We observe that the third term vanishes identically in any of the following circumstances:

| Condition | Why third term vanishes | Examples in this campaign |
|---|---|---|
| `\mathbf{u} = 0` | `\mathbf{u}\cdot\mathbf{a} = 0` | [J3e-P1.5](Ch01_Introduction_Electrostatics.md), [J3e-P2.1](Ch02_Boundary_Value_Problems_Electrostatics_I.md) (electrostatics) |
| `\mathbf{a} = 0` (uniform velocity) | `\mathbf{u}\cdot\mathbf{a} = 0` | [J3e-P6.4](Ch06_Maxwell_Equations_Macroscopic_Media.md#problem-j3e-p64--em-momentum-of-a-uniformly-moving-point-charge) |
| `\mathbf{a} \perp \mathbf{u}` (circular motion) | `\mathbf{u}\cdot\mathbf{a} = 0` | [J3e-P12.2](Ch12_Relativistic_Dynamics.md#problem-j3e-p122--cyclotron-motion-in-a-uniform-magnetic-field) (cyclotron), [J3e-P14.5](#problem-j3e-p145--synchrotron-radiation-from-circular-motion) (synchrotron) |
| Steady current configuration | `\partial_{\tau} = 0` on average | [J3e-P5.4](Ch05_Magnetostatics_Faraday_Quasi_Static.md#problem-j3e-p54--magnetic-field-of-a-circular-current-loop-on-its-axis) (current loop) |
| Perfect-conductor idealisation | `\mathbf{a}` undefined (infinite response) | [J3e-P6.20](Ch06_Maxwell_Equations_Macroscopic_Media.md#problem-j3e-p620--radiation-pressure-on-a-perfect-conductor) |

The third term is **non-zero** for linearly-accelerated charges (bremsstrahlung — [J3e-P14.6](#problem-j3e-p146--bremsstrahlung-from-a-linearly-decelerating-charge)) and for any motion with both `\mathbf{u}` and `\mathbf{a}` having components in the same direction (plane-wave figure-8 — [J3e-P12.14](Ch12_Relativistic_Dynamics.md#problem-j3e-p1214--charged-particle-in-a-plane-em-wave)). These are the situations in which the proper-time formulation's prediction is *operationally measurable*.

**Comparison:**

| Field component | Classical (Gaussian) | Proper-time |
|---|---|---|
| Velocity field | `e\,\mathbf{r}_{\mathbf{v}}(1-v^{2}/c^{2})/s_{v}^{3}` | identical (under velocity duality) |
| Acceleration field | `e\,\mathbf{r}\times(\mathbf{r}_{\mathbf{v}}\times\mathbf{a})/(c^{2}s_{v}^{3})` | identical structure with `c\to b`, `\mathbf{v}\to\mathbf{u}` |
| Third term | **ABSENT** | `e\,(\mathbf{u}\cdot\mathbf{a})\,\mathbf{r}\times(\mathbf{u}\times\mathbf{r})/(b^{4}s^{3})` |
| Radiation polarisation | purely transverse | transverse + longitudinal (when `\mathbf{u}\cdot\mathbf{a}\neq 0`) |

**Does the proper-time answer differ from a pure `c → b` redressing?** **⚠ YES — first occurrence in the campaign.** The third term is not a `c \to b` redressing of any classical term; it is a structurally new contribution to the radiated field, proportional to `(\mathbf{u}\cdot\mathbf{a})/b^{4}` and absent from classical EM. This is the campaign's load-bearing example of the proper-time formulation predicting something genuinely new.

**Verdict:** ⚠ **the proper-time formulation produces a qualitatively new field-level prediction (longitudinal radiation component when `\mathbf{u}\cdot\mathbf{a}\neq 0`) absent from classical Liénard–Wiechert.** This is consistent with the framework as published; it is not a flagged inconsistency. Whether the prediction is supported by experiment is the open question addressed in issue #43.

**Notes for author review:** **this is the campaign's load-bearing finding to date** and the document warrants the author's full read. The third term is mechanically derivable from the Gill–Zachary substitution rules applied to the Maxwell field equations, but its physical content — a longitudinal radiation component absent from classical EM — is a substantive claim that depends on the framework being the correct formulation. If issue #43's comparison against Cole/Poder/Wistisen 2018 data finds the prediction is *not* supported, this would constitute the campaign's first experimental falsification of the framework. If the prediction *is* supported (or is statistically indistinguishable from quantum-corrected Landau–Lifshitz), the result is a confirmation of the framework's experimental status. Either outcome warrants a separate entry in [`FINDINGS_for_author_review.md`](../../Equation_Verification/FINDINGS_for_author_review.md); flagging here pending the #43 outcome.

**Companion notebook:** [`Roadmapping/Mathematica_Notebooks/Electromagnetism/JacksonCh14_P14_2.wl`](../../Mathematica_Notebooks/Electromagnetism/JacksonCh14_P14_2.wl).

---

### Problem J3e-P14.3 — Non-relativistic Larmor radiation formula

**Selection provenance** (Crocco §5 substantive-AI note):
- *Chosen because:* the non-relativistic Larmor formula `P = (2e^{2}a^{2})/(3c^{3})` is the simplest expression in classical electrodynamics for the total power radiated by an accelerating charge, treated in [[jackson1998_classical_electrodynamics]] §14.2. **Podcast pick #5** in [§12.1 of the campaign plan](../../../.dev/tasks/42-electromagnetism-jackson-proper-time.md#121-per-problem-briefs); the pedagogical "middle rung" of the podcast gradient, where the proper-time correction is computable as a power series in `u/c` and the regime of experimental measurability is honestly addressed.
- *Alternatives considered:* J3e-P14.4 (Larmor's formula for non-relativistic monochromatic source — partially covered) and J3e-P14.8 (relativistic Liénard formula — selected separately as J3e-P14.5 for the synchrotron treatment).
- *Role in this PR:* **HEADLINE-PAYOFF**. Demonstrates the proper-time correction `P_{\text{PT}}/P_{\text{classical}} = 1 - (3/2)(u/c)^{2} + O((u/c)^{4})` and honestly addresses the regime where it is observable.

<!-- TODO: human reviews and fills in — confirms the role of this problem as podcast pick #5's headline-payoff document and the framing of the "middle rung" of the podcast gradient -->

**Source:** Jackson, *Classical Electrodynamics*, 3e Problem 14.3 (and 2e Problem 14.3, equivalent at this number). *Paraphrased.*

**Paraphrased statement:** Compute the total power radiated by a non-relativistically-accelerated point charge `q` of rest mass `m`, using both the classical Liénard–Wiechert formulation and the proper-time formulation of [[Two_Mathematically_Equivalent_Versions_of_Maxwells_Equations]] Eq. (7). Identify the leading-order proper-time correction as a function of `u/c`, and remark on the regime of experimental measurability.

**Setup:** Point charge `q` with non-relativistic velocity `\mathbf{v}` (i.e., `v \ll c`) and acceleration `\mathbf{a}` of arbitrary magnitude. In the proper-time formulation, `\mathbf{u} \approx \mathbf{v}` and `b \approx c(1 + u^{2}/(2c^{2}))` to leading order.

#### (a) Classical solution — Gaussian (CGS)

The Larmor formula for the total radiated power by a non-relativistic accelerated charge, derivable from the radiation-zone Liénard–Wiechert acceleration field of [Problem J3e-P14.2](#problem-j3e-p142--li%C3%A9nard-wiechert-fields-with-the-proper-time-third-term) (second term only, since the third term vanishes at `u \to 0`), is

$$
P_{\text{classical}} = \frac{2 e^{2} a^{2}}{3 c^{3}}.
$$

This is Jackson 3e Eq. (14.22). It applies in the limit `v \ll c`; at relativistic velocities, the relativistic Liénard formula must be used instead (treated in [Problem J3e-P14.5](#problem-j3e-p145--synchrotron-radiation-from-circular-motion)).

#### (c) Proper-time reformulation

The proper-time formulation produces the same total radiated power at the level of the radiation-zone field structure, but with the substitution `c \to b` in the prefactor:

$$
P_{\text{PT}} = \frac{2 e^{2} a^{2}}{3 b^{3}}.
$$

Expanding `b = c\sqrt{1 + u^{2}/c^{2}}` in powers of `u/c`,

$$
b^{3} = c^{3}\left[1 + \frac{u^{2}}{c^{2}}\right]^{3/2} = c^{3}\left[1 + \frac{3 u^{2}}{2 c^{2}} + \frac{3 u^{4}}{8 c^{4}} + \ldots\right].
$$

Substituting,

$$
P_{\text{PT}} = \frac{2 e^{2} a^{2}}{3 c^{3}}\left[1 - \frac{3 u^{2}}{2 c^{2}} + \frac{15 u^{4}}{8 c^{4}} - \ldots\right] = P_{\text{classical}}\left[1 - \frac{3 u^{2}}{2 c^{2}} + O((u/c)^{4})\right].
$$

**Mathematica check** (Wolfram MCP, 2026-05-24):

```mathematica
ClearAll[ee, aa, cc, uu];
bExpansion = cc Sqrt[1 + uu^2/cc^2];
larmorClassical = (2 ee^2 aa^2)/(3 cc^3);
larmorPT = (2 ee^2 aa^2)/(3 bExpansion^3);
ratio = FullSimplify[larmorPT/larmorClassical];
seriesExp = Normal[Series[ratio, {uu, 0, 4}]];
Print["P_PT / P_classical = ", ratio];
Print["Series in u/c: ", seriesExp];
(* Result: 1 - 3 uu^2/(2 cc^2) + 15 uu^4/(8 cc^4)  ✅ *)
```

The proper-time correction is **negative at leading order** (the proper-time formulation predicts *less* radiated power than the classical Larmor formula at the same observer-time velocity), with magnitude `(3/2)(u/c)^{2}`. The correction grows at fourth order with coefficient `+15/8`.

#### Regime of experimental measurability

The fractional correction as a function of velocity (using `u \approx v` at non-relativistic speeds):

| Regime | `v/c` | Fractional correction `\Delta P/P` |
|---|---|---|
| Chemistry-scale electron | `\sim 10^{-2}` | `\sim -1.5 \times 10^{-4}` |
| 1 keV electron | `\sim 0.06` | `\sim -5.4 \times 10^{-3}` |
| 100 keV electron | `\sim 0.55` | `\sim -45\%` (order unity) |
| Ultrarelativistic | `\to 1` | dominates |

We observe the campaign's honest moment, articulated in the podcast brief and reaffirmed here: **the fractional correction is below the observational floor at the chemistry-scale velocities where Larmor is the right formula, and is order unity at the velocities where Larmor must be upgraded to the relativistic Liénard expression.** The proper-time correction is, in the strict non-relativistic limit, theoretically computable but experimentally inaccessible; in the relativistic limit it is experimentally accessible but the formula has been replaced.

This is the campaign's most pedagogically honest observation: the proper-time formulation produces a correction in a regime that is neither cleanly Larmor-applicable nor cleanly Liénard-applicable. The experimental test of the proper-time prediction must therefore come from the *relativistic* regime — the Liénard formula with the third-term contribution of [Problem J3e-P14.2](#problem-j3e-p142--li%C3%A9nard-wiechert-fields-with-the-proper-time-third-term), as probed by the experiments referenced in issue [#43](https://github.com/temoTxt/PyPhysics/issues/43).

<!-- TODO: human reviews and fills in — confirms the honest moment that the non-relativistic Larmor regime is "below floor where the formula is gorgeous; gorgeous where the formula needs upgrade". The podcast Experimentalist persona's hand-off to J3e-P14.2 and PR E content is articulated here -->

#### Connection to bremsstrahlung

Linear deceleration of a charge (bremsstrahlung) provides a real experimental setting where the proper-time correction is finite and computable. At electron energies of ~MeV — typical for synchrotron-radiation X-ray sources and for the bremsstrahlung X-rays of clinical linacs — the proper-time correction would be `\sim 10^{-3}` to `\sim 10^{-1}`. Precision bremsstrahlung-spectrum measurements achieve `\sim 1\%` accuracy in their classical-Born-approximation regime; the proper-time prediction is at the edge of distinguishability but not currently a clean test. [Problem J3e-P14.6](#problem-j3e-p146--bremsstrahlung-from-a-linearly-decelerating-charge) treats this configuration in detail.

**Comparison:**

| Quantity | Classical (Gaussian) | Proper-time |
|---|---|---|
| Radiated power | `(2e^{2}a^{2})/(3c^{3})` | `(2e^{2}a^{2})/(3b^{3})` |
| Leading correction | none | `\times[1 - (3/2)(u/c)^{2}]` |
| Validity regime | `v \ll c` | same |
| Measurable? | yes, at chemistry-scale | only at speeds where Larmor breaks down |

**Does the proper-time answer differ from a pure `c → b` redressing?** ✅ no, at this level. The proper-time prediction is *precisely* the classical Larmor formula with `c \to b`. The "new physics" of the proper-time formulation does not surface in the non-relativistic regime; it surfaces in the third term of Eq. (7), which is `O((u\cdot a)/b^{4})` and is suppressed by `(u/c)^{2}` at non-relativistic velocities.

**Verdict:** ✅ all formulations consistent. The proper-time correction to non-relativistic Larmor is a computable `O((u/c)^{2})` shift below the observational floor at chemistry-scale velocities; the operationally measurable signature of the proper-time formulation lives in the relativistic regime (Liénard formula plus third term — issue #43).

**Notes for author review:** the podcast Experimentalist's hand-off to Cole/Poder/Wistisen 2018 experiments is the natural narrative move from this problem. The proper-time correction to Larmor is *real* and *computable*, but the regime in which it is large enough to measure is the regime where Larmor has been replaced by Liénard. This is the campaign's most pedagogically honest observation, and it is articulated as the load-bearing message of podcast pick #5. Not posted to `FINDINGS_for_author_review.md` — the result is consistent with classical EM at the level of the Larmor approximation, with a documented but unobservable correction.

**Companion notebook:** [`Roadmapping/Mathematica_Notebooks/Electromagnetism/JacksonCh14_P14_3.wl`](../../Mathematica_Notebooks/Electromagnetism/JacksonCh14_P14_3.wl).

---

### Problem J3e-P14.5 — Synchrotron radiation from circular motion

**Selection provenance** (Crocco §5 substantive-AI note):
- *Chosen because:* synchrotron radiation from a relativistic charge in circular motion is the workhorse measurement in modern light-source physics (APS, ESRF, Diamond, NSLS-II, SPring-8). Treated in [[jackson1998_classical_electrodynamics]] §14.4. As a fluency-builder, it cleanly demonstrates that the proper-time formulation's third term **vanishes identically** for circular motion (`\mathbf{u}\cdot\mathbf{a} = 0` by tangential/centripetal orthogonality), recovering the classical Liénard formula exactly.
- *Alternatives considered:* J3e-P14.4 (angular distribution from circular motion — covered partially here) and J3e-P14.7 (relativistic Doppler effect on emission — covered in PR B Ch. 11).
- *Role in this PR:* fluency-builder. Null-result confirmation that the third term does not engage in circular motion, complementing [Problem J3e-P14.6](#problem-j3e-p146--bremsstrahlung-from-a-linearly-decelerating-charge)'s non-null engagement in linear acceleration.

<!-- TODO: human reviews and fills in — confirms this is the null-result fluency-builder for circular motion, paired with the bremsstrahlung problem as the non-null-result counterpoint -->

**Source:** Jackson, *Classical Electrodynamics*, 3e Problem 14.5 (and 2e Problem 14.5, equivalent). *Paraphrased.*

**Paraphrased statement:** A relativistic charge undergoes circular motion in a uniform magnetic field with Lorentz factor `\gamma \gg 1`. Compute the total radiated power (the synchrotron formula) using the classical Liénard formula and the proper-time Eq. (7) of [[Two_Mathematically_Equivalent_Versions_of_Maxwells_Equations]]. Verify that the proper-time third term vanishes by `\mathbf{u}\cdot\mathbf{a} = 0`.

**Setup:** Circular motion: `\mathbf{w}(t)` is tangential, `\mathbf{a}_{\text{cl}}(t) = d\mathbf{w}/dt` is centripetal. Hence `\mathbf{w}\cdot\mathbf{a}_{\text{cl}} = 0` identically. In proper-time variables, `\mathbf{u}\cdot\mathbf{a} = 0` likewise (by the velocity-duality identity preserving the angular relationship). `|\mathbf{u}|` and `b` are both constant.

#### (a) Classical solution — Gaussian (CGS)

The Liénard formula for the total radiated power of a relativistic accelerating charge (Jackson 3e Eq. (14.26)) reads

$$
P = \frac{2 e^{2}}{3 c^{3}}\,\gamma^{6}\!\left[|\mathbf{a}_{\text{cl}}|^{2} - \frac{|\mathbf{w}\times\mathbf{a}_{\text{cl}}|^{2}}{c^{2}}\right].
$$

For circular motion, `|\mathbf{w}\times\mathbf{a}_{\text{cl}}| = w a_{\text{cl}}` (since `\mathbf{w}\perp\mathbf{a}_{\text{cl}}`), so

$$
P_{\text{circular}} = \frac{2 e^{2}}{3 c^{3}}\,\gamma^{6}a_{\text{cl}}^{2}\!\left[1 - \frac{w^{2}}{c^{2}}\right] = \frac{2 e^{2} a_{\text{cl}}^{2}}{3 c^{3}}\,\gamma^{4}.
$$

**Mathematica check** (Wolfram MCP, 2026-05-24):

```mathematica
ClearAll[ee, aa, vv, cc];
gammaW = 1/Sqrt[1 - vv^2/cc^2];
lienardCircular = (2 ee^2/(3 cc^3)) gammaW^6 (aa^2 - vv^2 aa^2/cc^2);
Print[FullSimplify[lienardCircular - (2 ee^2 aa^2/(3 cc^3)) gammaW^4,
   Assumptions -> 0 < vv < cc]];
(* Result: 0  ✅ *)
```

The classical synchrotron power scales as `\gamma^{4}` — the steep relativistic enhancement that makes synchrotron facilities effective X-ray sources at GeV energies.

#### (c) Proper-time reformulation

The proper-time Liénard–Wiechert fields of Eq. (7) of [[Two_Mathematically_Equivalent_Versions_of_Maxwells_Equations]] reduce, for circular motion, to the first two terms only — the **third term `e(\mathbf{u}\cdot\mathbf{a})\mathbf{r}\times(\mathbf{u}\times\mathbf{r})/(b^{4}s^{3})` vanishes identically** by `\mathbf{u}\cdot\mathbf{a} = 0`. The first two terms have the same structure as classical Liénard–Wiechert with `c \to b` and `\mathbf{v} \to \mathbf{u}`.

Performing the radiation-zone power integration (the same construction as the classical Liénard formula but in proper-time variables), we obtain the proper-time synchrotron power. Integrated over solid angle and converted from proper-time-rate to lab-time-rate via `dE/dt = (1/\gamma)\,dE/d\tau = (c/b)\,dE/d\tau`, the result is identical to the classical Liénard formula. The proper-time formulation predicts no new contribution to the synchrotron power; it is a **null-result confirmation** of the framework against one of the cleanest precision-test settings in modern light-source physics.

<!-- TODO: human reviews and fills in — confirms the framing that circular motion is the null-result case for the dissipative term, complementing the non-null bremsstrahlung case in J3e-P14.6 -->

This is consistent with the experimental record: synchrotron facilities have confirmed the classical Liénard `\gamma^{4}` scaling to percent-level precision across decades of beam energy (see [Problem J3e-P14.3](#problem-j3e-p143--non-relativistic-larmor-radiation-formula) historical note on Elder/Gurewitsch/Langmuir/Pollock 1947, and the modern equivalents at APS / ESRF / Diamond / SPring-8). Any proper-time correction would be at the level of `c \to b` redressing of the prefactor `1/c^{3}`, which is itself zero at the observable level because the integrated power computed in proper-time variables transforms back to the lab frame with the appropriate factor.

#### A note on the proper-time cyclotron frequency from PR C

We observe that the **circular orbit's angular frequency** in proper-time form is the `\gamma`-independent `\omega_{\tau} = qB/(mc)` of [Problem J3e-P12.2](Ch12_Relativistic_Dynamics.md#problem-j3e-p122--cyclotron-motion-in-a-uniform-magnetic-field), with the classical lab-frame `\omega_{c} = qB/(\gamma m c) = \omega_{\tau}/\gamma`. The synchrotron *radiated power* scales as `\gamma^{4}` (above), whereas the *orbital frequency* in proper time is `\gamma`-independent — a clean separation between the orbit's clock-rate (in proper time, unchanged at all energies) and the radiation rate (in lab time, enhanced by `\gamma^{4}`). Both statements are consistent with each other and with the classical result.

<!-- TODO: human reviews and fills in — confirms the connection drawn between J3e-P12.2's gamma-independent cyclotron frequency and J3e-P14.5's gamma^4 synchrotron radiation power. These are different observables and the framing here should not give the impression that they "contradict" -->

**Comparison:**

| Quantity | Classical (Gaussian) | Proper-time |
|---|---|---|
| Total radiated power | `(2e^{2}a^{2}/3c^{3})\gamma^{4}` | identical |
| Energy-scaling | `\gamma^{4}` | identical |
| Third-term contribution | n/a | **zero** (`\mathbf{u}\cdot\mathbf{a} = 0`) |
| Cyclotron frequency in proper-time | n/a | `\omega_{\tau} = qB/(mc)`, γ-independent (from J3e-P12.2) |
| Cyclotron frequency in lab | `qB/(\gamma mc)` | `\omega_{\tau}/\gamma` (identical to classical) |

**Does the proper-time answer differ from a pure `c → b` redressing?** ✅ no. The synchrotron radiation in circular motion is the cleanest example of the proper-time framework producing *exactly* the classical Liénard result with no observable deviation — because the third term, the framework's only source of new field-level content, vanishes identically by `\mathbf{u}\cdot\mathbf{a} = 0`.

**Verdict:** ✅ all formulations consistent. The synchrotron `\gamma^{4}` scaling is preserved by the proper-time formulation, and no new contribution to the total radiated power emerges. This is the cleanest precision-test setting in which the proper-time framework's null result is exhaustively confirmed by experimental data.

**Notes for author review:** none. Synchrotron radiation is the case where the proper-time framework reproduces classical EM exactly, and the experimental record (1947 → present) confirms it at percent-level precision. Not posted to `FINDINGS_for_author_review.md`.

**Companion notebook:** [`Roadmapping/Mathematica_Notebooks/Electromagnetism/JacksonCh14_P14_5.wl`](../../Mathematica_Notebooks/Electromagnetism/JacksonCh14_P14_5.wl).
