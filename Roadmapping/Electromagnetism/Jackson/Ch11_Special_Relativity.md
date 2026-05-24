# Ch. 11 — Special Theory of Relativity

This chapter contains Jackson canonical problems on special relativity, worked in the proper-time reformulation alongside their classical solutions. Per [§4 of the plan](../../../.dev/tasks/42-electromagnetism-jackson-proper-time.md#4-unit-system-handling), Ch. 11 of Jackson is presented in the two-system regime: Gaussian (CGS) only on the classical side, with a proper-time reformulation alongside. The 3rd edition reverts to Gaussian for Chs. 11+ for exactly the reason these problems make manifest: the Lorentz-covariant structure of relativistic electrodynamics is most cleanly expressed in Gaussian units, without the `c`-factor bookkeeping of SI.

Ch. 11 is where the **proper-time group** of [[Two_Mathematically_Equivalent_Versions_of_Maxwells_Equations]] §1.3 is exercised in detail. The structural observation surfaced in [Problem J3e-P6.11](Ch06_Maxwell_Equations_Macroscopic_Media.md#problem-j3e-p611--symmetric-stress-tensor-and-lorentz-behaviour) — that the proper-time group is distinct from the standard Lorentz group, with the source's local clock held fixed across observers rather than the observer's clock held fixed across sources — becomes the load-bearing context for every problem in this chapter. Per-paragraph `<!-- TODO -->` blocks continue per [§13.5 D2 of the plan](../../../.dev/tasks/42-electromagnetism-jackson-proper-time.md#135-decision-points--confirmed-by-author-2026-05-24).

## Problems in this chapter

| Problem | Status | Role |
|---|---|---|
| [Problem J3e-P11.1 — Relativistic velocity addition](#problem-j3e-p111--relativistic-velocity-addition) | drafted (PR B) | fluency-builder (velocity-duality exercise) |

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
