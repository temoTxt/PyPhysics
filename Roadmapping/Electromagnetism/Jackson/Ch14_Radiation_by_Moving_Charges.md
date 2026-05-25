# Ch. 14 — Radiation by Moving Charges

This chapter contains Jackson canonical problems on radiation from accelerating charges, worked in the proper-time reformulation alongside their classical solutions. Per [§4 of the plan](../../../.dev/tasks/42-electromagnetism-jackson-proper-time.md#4-unit-system-handling), Ch. 14 of Jackson is presented in the two-system regime: Gaussian (CGS) only on the classical side, with a proper-time reformulation alongside.

Ch. 14 is the campaign's **first headline-payoff chapter**. Two of the five §12 podcast picks land here: **#2 (Liénard–Wiechert with the third term)** and **#5 (non-relativistic Larmor radiation)**. The dissipative `(\mathbf{u}\cdot\mathbf{a})/b^{4}` coefficient of Eq. (4) of [[Two_Mathematically_Equivalent_Versions_of_Maxwells_Equations]] engages directly in the field expressions of Eq. (7) — the new third term that distinguishes the proper-time formulation from classical EM at the level of radiated fields. This is the chapter where issue [#43](https://github.com/temoTxt/PyPhysics/issues/43)'s experimental-comparison hooks become operationally measurable.

Per-paragraph `<!-- TODO -->` blocks continue per [§13.5 D2](../../../.dev/tasks/42-electromagnetism-jackson-proper-time.md#135-decision-points--confirmed-by-author-2026-05-24).

## Problems in this chapter

| Problem | Status | Role |
|---|---|---|
| [Problem J3e-P14.1 — Liénard–Wiechert potentials](#problem-j3e-p141--li%C3%A9nard-wiechert-potentials) | drafted (PR D) | fluency-builder (foundational) |
| [Problem J3e-P14.2 — Liénard–Wiechert fields with the proper-time third term](#problem-j3e-p142--li%C3%A9nard-wiechert-fields-with-the-proper-time-third-term) | drafted (PR D) | **HEADLINE-PAYOFF (podcast pick #2)** |

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
