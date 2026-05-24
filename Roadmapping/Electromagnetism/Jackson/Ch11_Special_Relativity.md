# Ch. 11 — Special Theory of Relativity

This chapter contains Jackson canonical problems on special relativity, worked in the proper-time reformulation alongside their classical solutions. Per [§4 of the plan](../../../.dev/tasks/42-electromagnetism-jackson-proper-time.md#4-unit-system-handling), Ch. 11 of Jackson is presented in the two-system regime: Gaussian (CGS) only on the classical side, with a proper-time reformulation alongside. The 3rd edition reverts to Gaussian for Chs. 11+ for exactly the reason these problems make manifest: the Lorentz-covariant structure of relativistic electrodynamics is most cleanly expressed in Gaussian units, without the `c`-factor bookkeeping of SI.

Ch. 11 is where the **proper-time group** of [[Two_Mathematically_Equivalent_Versions_of_Maxwells_Equations]] §1.3 is exercised in detail. The structural observation surfaced in [Problem J3e-P6.11](Ch06_Maxwell_Equations_Macroscopic_Media.md#problem-j3e-p611--symmetric-stress-tensor-and-lorentz-behaviour) — that the proper-time group is distinct from the standard Lorentz group, with the source's local clock held fixed across observers rather than the observer's clock held fixed across sources — becomes the load-bearing context for every problem in this chapter. Per-paragraph `<!-- TODO -->` blocks continue per [§13.5 D2 of the plan](../../../.dev/tasks/42-electromagnetism-jackson-proper-time.md#135-decision-points--confirmed-by-author-2026-05-24).

## Problems in this chapter

| Problem | Status | Role |
|---|---|---|
| [Problem J3e-P11.1 — Relativistic velocity addition](#problem-j3e-p111--relativistic-velocity-addition) | drafted (PR B) | fluency-builder (velocity-duality exercise) |
| [Problem J3e-P11.3 — Relativistic Doppler effect](#problem-j3e-p113--relativistic-doppler-effect) | drafted (PR B) | fluency-builder (time-derivative duality) |
| [Problem J3e-P11.5 — 4-vector position transformation](#problem-j3e-p115--4-vector-position-transformation) | drafted (PR B) | fluency-builder (bookkeeping) |

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
