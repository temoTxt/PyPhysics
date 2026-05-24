# Ch. 14 — Radiation by Moving Charges

This chapter contains Jackson canonical problems on radiation from accelerating charges, worked in the proper-time reformulation alongside their classical solutions. Per [§4 of the plan](../../../.dev/tasks/42-electromagnetism-jackson-proper-time.md#4-unit-system-handling), Ch. 14 of Jackson is presented in the two-system regime: Gaussian (CGS) only on the classical side, with a proper-time reformulation alongside.

Ch. 14 is the campaign's **first headline-payoff chapter**. Two of the five §12 podcast picks land here: **#2 (Liénard–Wiechert with the third term)** and **#5 (non-relativistic Larmor radiation)**. The dissipative `(\mathbf{u}\cdot\mathbf{a})/b^{4}` coefficient of Eq. (4) of [[Two_Mathematically_Equivalent_Versions_of_Maxwells_Equations]] engages directly in the field expressions of Eq. (7) — the new third term that distinguishes the proper-time formulation from classical EM at the level of radiated fields. This is the chapter where issue [#43](https://github.com/temoTxt/PyPhysics/issues/43)'s experimental-comparison hooks become operationally measurable.

Per-paragraph `<!-- TODO -->` blocks continue per [§13.5 D2](../../../.dev/tasks/42-electromagnetism-jackson-proper-time.md#135-decision-points--confirmed-by-author-2026-05-24).

## Problems in this chapter

| Problem | Status | Role |
|---|---|---|
| [Problem J3e-P14.1 — Liénard–Wiechert potentials](#problem-j3e-p141--li%C3%A9nard-wiechert-potentials) | drafted (PR D) | fluency-builder (foundational) |

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
