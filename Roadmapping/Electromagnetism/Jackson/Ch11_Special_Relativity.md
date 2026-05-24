# Ch. 11 — Special Theory of Relativity

This chapter contains Jackson canonical problems on special relativity, worked in the proper-time reformulation alongside their classical solutions. Per [§4 of the plan](../../../.dev/tasks/42-electromagnetism-jackson-proper-time.md#4-unit-system-handling), Ch. 11 of Jackson is presented in the two-system regime: Gaussian (CGS) only on the classical side, with a proper-time reformulation alongside. The 3rd edition reverts to Gaussian for Chs. 11+ for exactly the reason these problems make manifest: the Lorentz-covariant structure of relativistic electrodynamics is most cleanly expressed in Gaussian units, without the `c`-factor bookkeeping of SI.

Ch. 11 is where the **proper-time group** of [[Two_Mathematically_Equivalent_Versions_of_Maxwells_Equations]] §1.3 is exercised in detail. The structural observation surfaced in [Problem J3e-P6.11](Ch06_Maxwell_Equations_Macroscopic_Media.md#problem-j3e-p611--symmetric-stress-tensor-and-lorentz-behaviour) — that the proper-time group is distinct from the standard Lorentz group, with the source's local clock held fixed across observers rather than the observer's clock held fixed across sources — becomes the load-bearing context for every problem in this chapter. Per-paragraph `<!-- TODO -->` blocks continue per [§13.5 D2 of the plan](../../../.dev/tasks/42-electromagnetism-jackson-proper-time.md#135-decision-points--confirmed-by-author-2026-05-24).

## Problems in this chapter

| Problem | Status | Role |
|---|---|---|
| [Problem J3e-P11.1 — Relativistic velocity addition](#problem-j3e-p111--relativistic-velocity-addition) | drafted (PR B) | fluency-builder (velocity-duality exercise) |
| [Problem J3e-P11.3 — Relativistic Doppler effect](#problem-j3e-p113--relativistic-doppler-effect) | drafted (PR B) | fluency-builder (time-derivative duality) |
| [Problem J3e-P11.5 — 4-vector position transformation](#problem-j3e-p115--4-vector-position-transformation) | drafted (PR B) | fluency-builder (bookkeeping) |
| [Problem J3e-P11.6 — Lorentz boost of **E** and **B** fields](#problem-j3e-p116--lorentz-boost-of-e-and-b-fields) | drafted (PR B) | headline-adjacent (group-distinction sharpened) |

---

### Problem J3e-P11.1 — Relativistic velocity addition

**Selection provenance** (Crocco §5 substantive-AI note):
- *Chosen because:* relativistic velocity addition is the textbook entry point to special relativity and is recorded in [[jackson1998_classical_electrodynamics]] §11.3. As the first per-problem document in Ch. 11, it exercises the velocity-duality identity `w/c = u/b` of Eq. (1) of [[Two_Mathematically_Equivalent_Versions_of_Maxwells_Equations]] in the simplest possible setting — composition of two boosts.
- *Alternatives considered:* J3e-P11.2 (the relativistic Doppler shift — selected separately as the next PR B problem) and J3e-P11.4 (relativistic mass — deferred to PR C because of its closer connection to dynamics).
- *Role in this PR:* fluency-builder. The proper-time formulation gives the same observable lab-frame velocity composition as classical SR; the work of the problem is in verifying the algebraic equivalence and in articulating what the proper-time group adds to the standard Lorentz picture.

<!-- TODO: human reviews and fills in — confirms the selection of velocity addition as PR B's opening problem and the role of "fluency-builder" rather than "headline-payoff" given that proper-time SR predicts the same kinematics as classical SR -->

**Source:** Jackson, *Classical Electrodynamics*, 3e Problem 11.1 (and 2e Problem 11.1, equivalent). *Paraphrased; consult the textbook for the precise statement.*

**Paraphrased statement:** Frame `S'` moves with velocity `\mathbf{v}` (along `+\hat x`) relative to frame `S`. A particle has velocity `\mathbf{w}'` (along `+\hat x`) as observed in `S'`. Compute the particle's velocity `\mathbf{w}` as observed in `S`, using (i) the standard Lorentz velocity-addition formula, and (ii) the proper-time 4-velocity transformation. Verify that the two formulations give identical observable predictions.

**Setup:** All velocities are along `\hat x`. We use `w` for an observer-time velocity (`dx/dt`) and `u` for a proper-time velocity (`dx/d\tau`), with `u/b = w/c` and `b^{2} = c^{2} + u^{2}` per Eqs. (1)–(2) of the Maxwell paper. The Lorentz factor of the boost between frames is `\gamma(v) = 1/\sqrt{1 - v^{2}/c^{2}}`.

#### (a) Classical solution — Gaussian (CGS)

The Lorentz velocity-addition formula for observer-time velocities reads

$$
w = \frac{w' + v}{1 + w' v/c^{2}}.
$$

This is the standard result of Jackson 3e Eq. (11.31), derivable from the Lorentz transformation of the 4-vector position. The formula is symmetric in `w'` and `v` and reduces to the Galilean `w = w' + v` in the non-relativistic limit `w', v \ll c`. Its load-bearing feature is the denominator `1 + w'v/c^{2}`, which ensures that the addition of two velocities each less than `c` produces a result less than `c`.

#### (b) Classical solution — SI

The SI form of the addition formula is dimensionally identical to (a), since velocities have the same units in both systems and `c` is a universal constant. No translation is needed; the formula is `w = (w' + v)/(1 + w'v/c^{2})` in either unit system.

#### (c) Proper-time reformulation

In the proper-time formulation, the natural object is the 4-velocity `(b, \mathbf{u})` with invariant Minkowski-length-squared `b^{2} - \mathbf{u}^{2} = c^{2}` (deriving from `b^{2} = c^{2} + u^{2}`). Under a Lorentz boost with velocity `v` along `\hat x`, the components transform as

$$
b' = \gamma(v)\!\left[b - \frac{u_{x}\,v}{c}\right], \qquad u_{x}' = \gamma(v)\!\left[u_{x} - \frac{v\,b}{c}\right].
$$

These are Eq. (11) of [[Two_Mathematically_Equivalent_Versions_of_Maxwells_Equations]] and its complement, derived in [`Two_Mathematically_Equivalent_Versions_of_Maxwells_Equations.md` §Eq. (11)](../../Equation_Verification/Two_Mathematically_Equivalent_Versions_of_Maxwells_Equations.md) and verified ✅.

**Mathematica check** (Wolfram MCP, 2026-05-24): we verify three things — that the Minkowski-length invariant `b^{2} - u^{2} = c^{2}` is preserved by the boost; that the boost reproduces the lab-frame observer-time velocity `w_{\text{lab}}` from the classical addition formula; and that the velocity-duality identity `w/c = u/b` is consistent under composition.

```mathematica
ClearAll[uMoving, vRel, cc];
bMovingExpr = Sqrt[cc^2 + uMoving^2];
uLab = (1/Sqrt[1 - vRel^2/cc^2]) (uMoving + vRel bMovingExpr/cc);
bLab = (1/Sqrt[1 - vRel^2/cc^2]) (bMovingExpr + vRel uMoving/cc);
labInvariant = FullSimplify[bLab^2 - uLab^2 - cc^2,
   Assumptions -> 0 < vRel < cc && uMoving > 0];
wLab = cc uLab/bLab;
wMovingFromProper = cc uMoving/bMovingExpr;
wLabClassical = (wMovingFromProper + vRel)/(1 + wMovingFromProper vRel/cc^2);
Print["Lab-frame invariant b^2 - u^2 - c^2 = ", labInvariant];
Print["w_lab (proper-velocity composition) = w_lab (classical)?  ",
   FullSimplify[wLab - wLabClassical,
      Assumptions -> 0 < vRel < cc && uMoving > 0] === 0];
(* Lab-frame invariant: 0  ✅
   w_lab match: True  ✅ *)
```

The 4-velocity Minkowski invariant `b^{2} - u^{2} = c^{2}` is preserved by any Lorentz boost — a structural feature of the proper-time formulation that follows from the velocity duality alone. The lab-frame observer-time velocity `w_{\text{lab}} = c\,u_{\text{lab}}/b_{\text{lab}}` reproduces the classical velocity-addition formula exactly. We observe that the proper-time and classical formulations predict identical observable kinematics for any composition of boosts.

It is worth observing what the proper-time formulation adds, structurally, to the classical Lorentz picture. The Lorentz group, as it appears in the classical formulation, exchanges observer-time coordinates `(t, \mathbf{x})` between inertial frames; the proper-time `\tau` of a particle is computed from `dt` and `\mathbf{w}` by the standard `d\tau = dt\sqrt{1 - w^{2}/c^{2}}`. In the Gill–Zachary formulation, by contrast, the proper-time of the *source* is the primary coordinate, and the local-time of each particle is held fixed across observer frames per [§1.3 of the Maxwell paper](../../Equation_Verification/Two_Mathematically_Equivalent_Versions_of_Maxwells_Equations.md). The boost transformation Eq. (11) above is therefore not a standard Lorentz boost on `(t, \mathbf{x})`; it is the **proper-time-group boost on `(b, \mathbf{u})`**, related to but not identical with the Lorentz boost.

<!-- TODO: human reviews and fills in — confirms the framing that the boost transformation Eq. (11) is the proper-time-group's analogue of the Lorentz boost, rather than the Lorentz boost itself. This is the load-bearing structural claim for Ch. 11 and was first surfaced in J3e-P6.11; sharpening here -->

The observable consequence of this structural distinction is *not* a different velocity-addition formula at the kinematic level — the addition rule for the lab-time `w` is the standard one above, and the rule for the proper-time `u` is mechanically derivable from it by the velocity duality. The distinction shows up only when one asks about the *covariance* of further constructions: the stress tensor of [Problem J3e-P6.11](Ch06_Maxwell_Equations_Macroscopic_Media.md#problem-j3e-p611--symmetric-stress-tensor-and-lorentz-behaviour), the radiation reaction of PR D, and the dual-theory Hamiltonian of [Problem J3e-P12.1](Ch12_Relativistic_Dynamics.md) (in PR C) all live in spaces where the proper-time group acts, and their transformation properties cannot be read off from the standard Lorentz group alone.

<!-- TODO: human reviews and fills in — confirms the claim that the proper-time group's distinctness from the Lorentz group becomes operationally relevant only in the stress-tensor / radiation-reaction / Hamiltonian contexts. This is a load-bearing forward-looking claim about PR C/D/E content -->

**Comparison:**

| Quantity | Classical (Gaussian) | Classical (SI) | Proper-time |
|---|---|---|---|
| Velocity addition formula | `w = (w' + v)/(1 + w'v/c^{2})` | same | same observable `w`, derived from Eq. (11) 4-velocity boost |
| 4-velocity invariant | `\gamma^{2}c^{2}(1 - w^{2}/c^{2}) = c^{2}` | same | `b^{2} - u^{2} = c^{2}` (equivalent under `\gamma c = b, \gamma w = u`) |
| Covariance group | Lorentz | Lorentz | proper-time group of Maxwell paper §1.3 |
| Observable kinematics | Lorentz-invariant | same | identical to classical |

**Does the proper-time answer differ from a pure `c → b` redressing?** ✅ no, at the level of observable kinematics. The proper-time formulation produces the same lab-frame velocity composition as classical SR, by the velocity-duality identity. The structural difference — the covariance group — is not observable at the kinematic level; it surfaces only in derived structures (the stress tensor, the radiation reaction force, the Hamiltonian).

**Verdict:** ✅ all three solutions consistent at the kinematic level. The proper-time group's distinctness from the Lorentz group is not observable in velocity-addition problems; it becomes operationally relevant in PR C (dynamics), PR D (radiation), and PR E (radiation damping). Ch. 11 will exercise the structural distinction at the level of field transformations (subsequent problems in this chapter).

**Notes for author review:** the framing of "proper-time group ≠ Lorentz group" is sharpened in this problem relative to the seed observation in [Problem J3e-P6.11](Ch06_Maxwell_Equations_Macroscopic_Media.md#problem-j3e-p611--symmetric-stress-tensor-and-lorentz-behaviour). I have not posted it to [`FINDINGS_for_author_review.md`](../../Equation_Verification/FINDINGS_for_author_review.md) because it is a structural feature of the framework's geometric content, not an unresolved inconsistency. If the subsequent Ch. 11 problems (Lorentz boost of fields, Lorentz invariants) sharpen the observation further, the findings document may be updated then.

**Companion notebook:** [`Roadmapping/Mathematica_Notebooks/Electromagnetism/JacksonCh11_P11_1.wl`](../../Mathematica_Notebooks/Electromagnetism/JacksonCh11_P11_1.wl) — runnable independent of the Mathematica MCP.

---

### Problem J3e-P11.3 — Relativistic Doppler effect

**Selection provenance** (Crocco §5 substantive-AI note):
- *Chosen because:* the relativistic Doppler shift is the canonical "moving source / stationary observer" problem and the simplest exercise of the time-derivative duality `(1/c)\partial_{t} = (1/b)\partial_{\tau}` of Eq. (2) of [[Two_Mathematically_Equivalent_Versions_of_Maxwells_Equations]]. Treated in [[jackson1998_classical_electrodynamics]] §11.4.
- *Alternatives considered:* J3e-P11.4 (transverse Doppler — defer to PR F+ as a finer fluency-builder).
- *Role in this PR:* fluency-builder.

<!-- TODO: human reviews and fills in — confirms the role of this problem as a time-derivative-duality exercise -->

**Source:** Jackson, *Classical Electrodynamics*, 3e Problem 11.3 (and 2e Problem 11.3, equivalent). *Paraphrased.*

**Paraphrased statement:** A monochromatic source emits radiation at angular frequency `\omega'` in its rest frame. The source moves with velocity `\mathbf{v} = v\hat z` directly toward an observer at rest. Compute the angular frequency `\omega` observed in the lab frame, using both the standard Lorentz boost of the photon 4-vector and the proper-time formulation.

**Setup:** Photon 4-vector in the source's rest frame: `k'^{\mu} = (\omega'/c, k'\hat z)` with `k' = \omega'/c` (null 4-vector). The source moves with lab-velocity `v`, so the proper-velocity of a source fluid element is `u = \gamma(v) v` and the corresponding `b = \gamma(v) c = \sqrt{c^{2} + u^{2}}`.

#### (a) Classical solution — Gaussian (CGS)

The lab-frame angular frequency follows from the Lorentz boost of the photon 4-vector. For longitudinal motion (source approaching observer along `+\hat z`),

$$
\omega = \omega' \sqrt{\frac{1 + \beta}{1 - \beta}}, \qquad \beta = \frac{v}{c}.
$$

This is the standard Doppler shift formula of Jackson 3e Eq. (11.30). It is the product `\gamma(1+\beta)` of the time-dilation factor `\gamma` and the geometric "leading-edge" factor `(1 + \beta)`, reduced to a single square-root form.

#### (c) Proper-time reformulation

In the proper-time variables `(b, u)` with `b^{2} = c^{2} + u^{2}`, the Doppler formula takes a particularly clean form. Using `\beta = w/c = u/b` and `\gamma = b/c`, one finds `\gamma(1 - \beta) = (b - u)/c`, and the lab-frame frequency reads

$$
\omega = \omega'\,\frac{c}{b - u} = \omega'\,\frac{b + u}{c}.
$$

The two forms are equivalent under `b^{2} - u^{2} = c^{2}`, since `c^{2}/(b-u) = (b+u)`. The proper-time formulation exchanges the classical `\sqrt{(1+\beta)/(1-\beta)}` for `(b+u)/c`, which is a linear function of the proper velocity rather than a square-root function of the observer-time velocity. It is worth observing that the proper-time form makes manifest the photon's linear coupling to the source's 4-velocity component along the line of sight.

**Mathematica check** (Wolfram MCP, 2026-05-24):

```mathematica
ClearAll[uu, cc, omegaPrime];
bSubstituted = Sqrt[cc^2 + uu^2];
classicalDoppler = omegaPrime Sqrt[(1 + uu/bSubstituted)/(1 - uu/bSubstituted)];
properTimeDoppler = omegaPrime cc/(bSubstituted - uu);
altForm = omegaPrime (bSubstituted + uu)/cc;
Print[FullSimplify[classicalDoppler^2 - altForm^2,
   Assumptions -> uu > 0 && cc > 0]];
Print[FullSimplify[properTimeDoppler^2 - altForm^2,
   Assumptions -> uu > 0 && cc > 0]];
(* Both: 0  ✅ *)
```

All three forms agree algebraically (after squaring) and numerically (sample `u = c`, `c = 1`, `\omega' = 1`: all three give `\omega = 1 + \sqrt 2 \approx 2.414`).

<!-- TODO: human reviews and fills in — confirms the framing that the proper-time form `(b+u)/c` is "linear in the source 4-velocity" and that this is worth highlighting as a structural feature of the formulation -->

**Comparison:**

| Quantity | Classical (Gaussian) | Proper-time |
|---|---|---|
| Doppler factor | `\sqrt{(1+\beta)/(1-\beta)}` | `(b+u)/c` |
| Observable `\omega` | identical | identical |

**Does the proper-time answer differ from a pure `c → b` redressing?** ✅ no, at the kinematic level. The proper-time form is a re-expression of the classical Doppler shift in proper-velocity variables, with no observable difference.

**Verdict:** ✅ all formulations consistent. The proper-time form `(b+u)/c` is a structurally cleaner re-expression of the classical Doppler factor; the observable frequency is unchanged.

**Notes for author review:** the observation that the proper-time Doppler formula is *linear* in `(b, u)` rather than square-root in `\beta` is worth recording — it is the simplest example in the campaign of the proper-time variables producing a cleaner algebraic form than the observer-time variables. Not posted to `FINDINGS_for_author_review.md`; structural feature, not a finding.

**Companion notebook:** [`Roadmapping/Mathematica_Notebooks/Electromagnetism/JacksonCh11_P11_3.wl`](../../Mathematica_Notebooks/Electromagnetism/JacksonCh11_P11_3.wl).

---

### Problem J3e-P11.5 — 4-vector position transformation

**Selection provenance** (Crocco §5 substantive-AI note):
- *Chosen because:* the Lorentz transformation of the position 4-vector is the foundational kinematic statement of special relativity, treated in [[jackson1998_classical_electrodynamics]] §11.3 and §11.6. As the kinematic backbone of every relativistic derivation in Jackson, it is the cleanest place to articulate the relationship between the standard Lorentz group acting on `x^{\mu} = (ct, \mathbf{x})` and the proper-time group of [[Two_Mathematically_Equivalent_Versions_of_Maxwells_Equations]] §1.3 acting on the 4-velocity `(b, \mathbf{u})`.
- *Alternatives considered:* J3e-P11.6 (Lorentz boost of `\mathbf{E}` and `\mathbf{B}` — selected next as the headline-adjacent problem for PR B) and J3e-P11.7 (Lorentz invariants — selected as the closing problem of PR B).
- *Role in this PR:* fluency-builder. The position 4-vector transforms identically under both the Lorentz group and the proper-time group at the level of components; the structural distinction surfaces only when one passes from positions to 4-velocities and to derived quantities.

<!-- TODO: human reviews and fills in — confirms the framing that the 4-position transforms identically under both group structures, with the distinction appearing only at the 4-velocity level -->

**Source:** Jackson, *Classical Electrodynamics*, 3e Problem 11.5 (and 2e Problem 11.5, equivalent). *Paraphrased.*

**Paraphrased statement:** Frame `S'` moves with velocity `\mathbf{v} = v\hat x` relative to frame `S`. Write down the Lorentz transformation of the 4-vector position `x^{\mu} = (ct, x, y, z)`, verify that the Minkowski interval `c^{2}t^{2} - x^{2} - y^{2} - z^{2}` is invariant, and remark on how the same transformation is understood in the proper-time formulation.

**Setup:** The Lorentz boost along `\hat x` with velocity `v` has parameter `\gamma(v) = 1/\sqrt{1 - v^{2}/c^{2}}`. Position 4-vector in frame `S`: `(ct, x, y, z)`. Position 4-vector in frame `S'`: `(ct', x', y', z')`. The transverse components `y` and `z` are unchanged by the boost.

#### (a) Classical solution — Gaussian (CGS)

The standard Lorentz transformation, derivable from the postulates of special relativity (Jackson 3e Eqs. (11.16)–(11.17)), reads

$$
ct' = \gamma(v)\!\left[ct - \beta x\right], \qquad x' = \gamma(v)\!\left[x - \beta\,ct\right], \qquad y' = y, \quad z' = z,
$$

with `\beta = v/c`. The Minkowski interval `s^{2} = c^{2}t^{2} - x^{2} - y^{2} - z^{2}` is invariant under this transformation, a property recorded as the kinematic statement of the relativity postulate.

**Mathematica check** (Wolfram MCP, 2026-05-24):

```mathematica
ClearAll[tt, xx, yy, zz, vv, cc];
gammaV = 1/Sqrt[1 - vv^2/cc^2];
tPrime = gammaV (tt - vv xx/cc^2);
xPrime = gammaV (xx - vv tt);
interval = FullSimplify[
   cc^2 tPrime^2 - xPrime^2 - yy^2 - zz^2
      - (cc^2 tt^2 - xx^2 - yy^2 - zz^2),
   Assumptions -> 0 < vv < cc];
(* Result: 0  ✅ *)
```

The interval is exactly invariant; the Lorentz boost is therefore a proper isometry of Minkowski spacetime.

#### (c) Proper-time reformulation

The position 4-vector `x^{\mu}` is defined by an observer in a particular inertial frame; its components have no dependence on which clock parametrises a particle's worldline. We observe that the Lorentz transformation above acts identically in the classical and proper-time formulations: the components of `x^{\mu}` transform the same way regardless of whether we use the observer's `t` or the source's `\tau` as the affine parameter along any particular worldline.

What does change between the two formulations is the **affine parameter of the worldline itself**. Classical SR uses the proper time of a particle, computed from observer-frame quantities as `d\tau = dt\sqrt{1 - w^{2}/c^{2}}` where `w = |d\mathbf{x}/dt|`. In the Gill–Zachary formulation, `\tau` is the source's local clock, taken as primary, and the relation `(1/c)\partial_{t} = (1/b)\partial_{\tau}` of Eq. (2) of [[Two_Mathematically_Equivalent_Versions_of_Maxwells_Equations]] expresses the connection between the two clocks at the level of derivatives.

At the level of the 4-position `x^{\mu}`, both formulations agree on the components and on the Minkowski-interval invariance. The structural distinction — Lorentz group vs proper-time group — surfaces when one passes from `x^{\mu}` to derived 4-vectors. The 4-velocity `(b, \mathbf{u})` transforms under the proper-time group via Eq. (11) (verified in [Problem J3e-P11.1](#problem-j3e-p111--relativistic-velocity-addition)), which is related to but distinct from the standard Lorentz boost on `\dot{x}^{\mu} = (\gamma c, \gamma\mathbf{w})`. The identity `\gamma c = b` and `\gamma\mathbf{w} = \mathbf{u}` reconciles the two at the level of components but does *not* reconcile them at the level of which group acts as the fundamental symmetry of the dynamics.

<!-- TODO: human reviews and fills in — confirms the framing that "Lorentz boost on position acts identically in both formulations; the group distinction lives at the 4-velocity level and beyond". This is the load-bearing kinematic claim for Ch. 11 -->

**Comparison:**

| Quantity | Classical (Gaussian) | Proper-time |
|---|---|---|
| Position transformation | `ct' = \gamma(ct - \beta x), x' = \gamma(x - \beta ct)` | identical |
| Minkowski interval | invariant | invariant |
| Affine parameter | observer-time `t`, with `\tau` derived | source-time `\tau` primary, with `t` related via Eq. (2) |
| Covariance group acting | Lorentz | proper-time group of Maxwell paper §1.3 |

**Does the proper-time answer differ from a pure `c → b` redressing?** ✅ no, at the position level. The Lorentz boost on `x^{\mu}` is the same operation in either formulation. The choice of affine parameter (and hence of fundamental symmetry group) is independent of the position transformation itself.

**Verdict:** ✅ all formulations consistent at the kinematic 4-position level. The proper-time group's distinctness from the Lorentz group, established in [Problem J3e-P6.11](Ch06_Maxwell_Equations_Macroscopic_Media.md#problem-j3e-p611--symmetric-stress-tensor-and-lorentz-behaviour) and sharpened in [Problem J3e-P11.1](#problem-j3e-p111--relativistic-velocity-addition), does not surface here because positions are 4-vectors irrespective of which clock parametrises the worldline.

**Notes for author review:** none. The 4-position transformation is the cleanest kinematic statement in which the two formulations agree exactly, with no daylight between them at the component level.

**Companion notebook:** [`Roadmapping/Mathematica_Notebooks/Electromagnetism/JacksonCh11_P11_5.wl`](../../Mathematica_Notebooks/Electromagnetism/JacksonCh11_P11_5.wl).

---

### Problem J3e-P11.6 — Lorentz boost of **E** and **B** fields

**Selection provenance** (Crocco §5 substantive-AI note):
- *Chosen because:* the Lorentz transformation of the electromagnetic field tensor is the standard test of field-theoretic Lorentz covariance, treated in [[jackson1998_classical_electrodynamics]] §11.10. It is PR B's headline-adjacent problem because it exercises the proper-time group on the *fields* (rather than on positions or velocities) and surfaces the cleanest algebraic re-expression of the boost in `(b, u)` variables.
- *Alternatives considered:* J3e-P11.7 (Lorentz invariants only — selected as the closing problem of PR B because it is a structural confirmation of what J3e-P11.6 already establishes) and J3e-P11.13 (field of a uniformly moving charge — defer to PR F+ as a derived application).
- *Role in this PR:* headline-adjacent. The classical `(\gamma, \gamma\beta)` form of the field-boost equations rewrites in `(b/c, u/c)` form via the velocity-duality identity; the resulting expressions are *linear* in `(b, u)`, paralleling the structural observation in [Problem J3e-P11.3](#problem-j3e-p113--relativistic-doppler-effect) about the Doppler factor.

<!-- TODO: human reviews and fills in — confirms the role of this problem as headline-adjacent for PR B and the framing of the (b, u) form as "structurally linear" -->

**Source:** Jackson, *Classical Electrodynamics*, 3e Problem 11.6 (and 2e Problem 11.6, equivalent). *Paraphrased.*

**Paraphrased statement:** Frame `S'` moves with velocity `\mathbf{v} = v\hat x` relative to frame `S`. Given the electric field `\mathbf{E}` and magnetic field `\mathbf{B}` measured in `S`, derive the corresponding fields `\mathbf{E}'` and `\mathbf{B}'` measured in `S'`. Verify that the scalar invariant `\mathbf{E}^{2} - \mathbf{B}^{2}` and the pseudoscalar invariant `\mathbf{E}\cdot\mathbf{B}` are preserved by the transformation. Carry out the analysis in both standard `(\gamma, \gamma\beta)` variables and the proper-time `(b, u)` variables.

**Setup:** As in [Problem J3e-P11.5](#problem-j3e-p115--4-vector-position-transformation), `\gamma = \gamma(v)` and `\beta = v/c`. The fields are evaluated at a single spacetime event; we work component-wise.

#### (a) Classical solution — Gaussian (CGS)

The standard Lorentz transformation of the electromagnetic field, recorded as Jackson 3e Eqs. (11.149), reads (with the boost along `+\hat x`):

$$
\begin{aligned}
E_{x}' &= E_{x}, & E_{y}' &= \gamma\!\left(E_{y} - \beta\,B_{z}\right), & E_{z}' &= \gamma\!\left(E_{z} + \beta\,B_{y}\right),\\
B_{x}' &= B_{x}, & B_{y}' &= \gamma\!\left(B_{y} + \beta\,E_{z}\right), & B_{z}' &= \gamma\!\left(B_{z} - \beta\,E_{y}\right).
\end{aligned}
$$

The longitudinal components `E_{x}` and `B_{x}` are unchanged; the transverse components mix `\mathbf{E}` and `\mathbf{B}` via the boost factor. The two Lorentz invariants — the scalar `I_{1} = \mathbf{E}^{2} - \mathbf{B}^{2}` and the pseudoscalar `I_{2} = \mathbf{E}\cdot\mathbf{B}` — are preserved by direct algebraic check.

**Mathematica check** (Wolfram MCP, 2026-05-24):

```mathematica
ClearAll[capE1, capE2, capE3, capB1, capB2, capB3, vv, cc];
gammaV = 1/Sqrt[1 - vv^2/cc^2]; beta = vv/cc;
eP = {capE1, gammaV (capE2 - beta capB3), gammaV (capE3 + beta capB2)};
bP = {capB1, gammaV (capB2 + beta capE3), gammaV (capB3 - beta capE2)};
diffI1 = FullSimplify[(eP.eP - bP.bP)
   - ((capE1^2 + capE2^2 + capE3^2) - (capB1^2 + capB2^2 + capB3^2)),
   Assumptions -> 0 < vv < cc];
diffI2 = FullSimplify[eP.bP
   - (capE1 capB1 + capE2 capB2 + capE3 capB3),
   Assumptions -> 0 < vv < cc];
(* diffI1 = 0, diffI2 = 0  ✅ — both invariants preserved *)
```

#### (c) Proper-time reformulation

Using the velocity-duality identities `u/b = w/c = \beta` and `b/c = \gamma` (with `b^{2} = c^{2} + u^{2}`), the field-transformation equations rewrite in the proper-time variables `(b, u)` as

$$
\begin{aligned}
E_{x}' &= E_{x}, & E_{y}' &= \frac{1}{c}\!\left(b\,E_{y} - u\,B_{z}\right), & E_{z}' &= \frac{1}{c}\!\left(b\,E_{z} + u\,B_{y}\right),\\
B_{x}' &= B_{x}, & B_{y}' &= \frac{1}{c}\!\left(b\,B_{y} + u\,E_{z}\right), & B_{z}' &= \frac{1}{c}\!\left(b\,B_{z} - u\,E_{y}\right).
\end{aligned}
$$

We observe that the proper-time form is **linear** in the boost variables `(b, u)`, in contrast to the classical form's non-linear `(\gamma, \gamma\beta)` dependence on `v`. This parallels the structural observation noted in [Problem J3e-P11.3](#problem-j3e-p113--relativistic-doppler-effect): the proper-time variables produce algebraically cleaner expressions because they remove the explicit square-root that the classical `\gamma` factor introduces.

The transformation is, of course, the same Lorentz transformation as in (a) — it acts on the same field tensor `F^{\mu\nu}` and produces the same boosted fields. The two forms differ only in how the boost parameter is expressed; the underlying physical content is unchanged.

**Mathematica check** (Wolfram MCP, 2026-05-24):

```mathematica
ClearAll[uu, bb, cc];
ePProperTime = {capE1, (1/cc) (bb capE2 - uu capB3), (1/cc) (bb capE3 + uu capB2)};
bPProperTime = {capB1, (1/cc) (bb capB2 + uu capE3), (1/cc) (bb capB3 - uu capE2)};
diffI1PT = FullSimplify[
   (ePProperTime.ePProperTime - bPProperTime.bPProperTime)
      - ((capE1^2 + capE2^2 + capE3^2) - (capB1^2 + capB2^2 + capB3^2))
   /. bb -> Sqrt[cc^2 + uu^2],
   Assumptions -> uu > 0 && cc > 0];
diffI2PT = FullSimplify[
   ePProperTime.bPProperTime
      - (capE1 capB1 + capE2 capB2 + capE3 capB3)
   /. bb -> Sqrt[cc^2 + uu^2],
   Assumptions -> uu > 0 && cc > 0];
(* diffI1PT = 0, diffI2PT = 0  ✅ — both invariants preserved in (b, u) form *)
```

The two Lorentz invariants are preserved in either parametrisation of the boost, as the consistency of the velocity-duality identity requires.

<!-- TODO: human reviews and fills in — confirms the framing that the (b, u) form of the field-boost is "structurally linear" and that this is worth highlighting alongside the analogous Doppler observation -->

It is worth observing what this means for the proper-time group's action on the field tensor. The fields `\mathbf{E}` and `\mathbf{B}` are observer-frame quantities: they are measured at fixed `\mathbf{x}` over an interval of observer time `t`. Under a boost to a frame moving with velocity `v`, the same fields are re-expressed in the moving frame's coordinates. This transformation is the same operation whether one parameterises by `(v, \gamma)` or by `(u, b)`; the proper-time group's distinctness from the Lorentz group, established at the 4-velocity level in [Problem J3e-P11.1](#problem-j3e-p111--relativistic-velocity-addition), does *not* produce a different transformation of the fields. Both groups act on `F^{\mu\nu}` identically when restricted to inertial boosts.

The proper-time group's distinct action surfaces only when one considers boosts between frames in which **the source is not co-moving with one of the two observers**. In that more general setting, the proper-time group preserves the source's local clock as the fundamental kinematic invariant, whereas the Lorentz group preserves the Minkowski interval of observers. For boosts between two observer frames (both treating the same source as external), the two groups agree. For boosts that bring an observer into the source's instantaneous rest frame, they coincide. For boosts between two source frames (rare in classical EM, but conceivable for multi-source problems), they differ — and this is the multi-source ambiguity surfaced in [Problem J3e-P6.1](Ch06_Maxwell_Equations_Macroscopic_Media.md#problem-j3e-p61--maxwell-with-magnetic-monopoles-and-electricmagnetic-duality).

<!-- TODO: human reviews and fills in — flags the multi-source-boost case as the place where the proper-time and Lorentz groups truly diverge, connecting the present problem to the earlier multi-species observation -->

**Comparison:**

| Quantity | Classical `(\gamma, \beta)` | Proper-time `(b, u)` |
|---|---|---|
| `E_y'` | `\gamma(E_y - \beta B_z)` | `(b E_y - u B_z)/c` |
| `B_y'` | `\gamma(B_y + \beta E_z)` | `(b B_y + u E_z)/c` |
| Functional form | non-linear in `v` (via `\gamma`) | linear in `(b, u)` |
| `I_1 = E^2 - B^2` | preserved | preserved |
| `I_2 = E\cdot B` | preserved | preserved |

**Does the proper-time answer differ from a pure `c → b` redressing?** ✅ no, at the level of single-source inertial boosts. The proper-time form `(b E - u B)/c` is a re-expression of the classical `\gamma(E - \beta B)` using the velocity-duality identity; the underlying transformation is the same.

**Verdict:** ✅ all formulations consistent. The proper-time form is algebraically cleaner (linear vs square-root) but represents the same Lorentz boost as the classical form. The proper-time group's distinct action on the field tensor surfaces only in multi-source problems, where the source's local clock is not unique — as noted in [Problem J3e-P6.1](Ch06_Maxwell_Equations_Macroscopic_Media.md#problem-j3e-p61--maxwell-with-magnetic-monopoles-and-electricmagnetic-duality).

**Notes for author review:** the connection between the field-boost linearity in `(b, u)` and the structurally linear Doppler factor of [Problem J3e-P11.3](#problem-j3e-p113--relativistic-doppler-effect) is worth recording as a pattern. The proper-time variables consistently produce linear expressions where the classical `(\gamma, \beta)` variables produce square-roots, suggesting that `(b, u)` are the natural coordinates for relativistic kinematics in the Gill–Zachary framework. Not posted to `FINDINGS_for_author_review.md`; this is a structural observation about variable choice, not a finding about physics.

**Companion notebook:** [`Roadmapping/Mathematica_Notebooks/Electromagnetism/JacksonCh11_P11_6.wl`](../../Mathematica_Notebooks/Electromagnetism/JacksonCh11_P11_6.wl).
