# Ch. 5 — Magnetostatics, Faraday's Law, Quasi-Static Fields

This chapter contains Jackson canonical problems on magnetostatics, worked in the proper-time reformulation alongside their classical CGS and SI solutions. Per [§4 of the plan](../../../.dev/tasks/42-electromagnetism-jackson-proper-time.md#4-unit-system-handling), Ch. 5 of Jackson is presented in the three-system regime.

The Ch. 5 problems below are part of PR 0 (the fluency-warm-up prequel — see [§7.3 of the plan](../../../.dev/tasks/42-electromagnetism-jackson-proper-time.md#73-pr-a-prequel-adopted)). They differ from the Ch. 1 and Ch. 2 problems in that the source carries a steady current, so the proper-time formulation engages the current-density rescaling rule `J → (b/c)\,J = \rho\mathbf{u}`. As we shall see, this rescaling cancels exactly with the `1/b` prefactor of the proper-time Ampère–Maxwell law, leaving the magnetic field identical to the classical result. The fluency lesson is therefore not that the rescaling produces a small correction, but that the rescaling is *algebraically* annihilated for any time-independent current — a stronger statement than the static-case lesson of Chs. 1–2.

## Problems in this chapter

| Problem | Status | Role |
|---|---|---|
| [Problem J3e-P5.4 — Magnetic field of a circular current loop on its axis](#problem-j3e-p54--magnetic-field-of-a-circular-current-loop-on-its-axis) | drafted (PR 0) | fluency-builder |

---

### Problem J3e-P5.4 — Magnetic field of a circular current loop on its axis

**Selection provenance** (Crocco §5 substantive-AI note):
- *Chosen because:* the circular current loop is the simplest non-trivial magnetostatic configuration, treated in [[jackson1998_classical_electrodynamics]] §5.4 and recorded in every graduate textbook of classical electrodynamics. It is the first per-problem document in the campaign in which a current density appears, and so the first occasion on which the proper-time substitution `J → (b/c)\,J` engages.
- *Alternatives considered:* J3e-P5.6 (B-field of a long straight wire — too symmetric to test the substitution rule with any geometric content) and J3e-P5.13 (B inside a uniformly rotating charged sphere — selected separately as PR 0's optional fourth problem).
- *Role in this PR:* fluency-builder. The configuration is time-independent, so `∂_\tau = 0`; the only proper-time machinery that could engage is the `(b/c)` current-density rescaling, which we show is algebraically cancelled by the `1/b` prefactor of the proper-time Ampère–Maxwell law.

<!-- TODO: human reviews and fills in — confirms the selection of the current-loop problem (and the rejection of the straight-wire problem) before this problem is rolled into PR 0's commit -->

**Source:** Jackson, *Classical Electrodynamics*, 3e Problem 5.4 (and 2e Problem 5.4, equivalent at this number). *Paraphrased; consult the textbook for the precise statement.*

**Paraphrased statement:** A circular loop of radius `R` carries a steady current `I`. Compute the magnetic field on the axis of the loop, at perpendicular distance `z` from the loop's centre.

**Setup:** Place the loop in the `xy`-plane, centred at the origin, with current `I` flowing in the `+\hat\phi` direction. The field point lies on the symmetry axis at `(0, 0, z)`. The configuration is time-independent in the observer's frame; the current carriers (typically conduction electrons) drift along the loop with some velocity `\mathbf{w}`, but for the macroscopic current density we are interested in only the steady value `\mathbf{J}`, not the carrier dynamics. By the standard fluid-element interpretation of the proper-time formulation in [[Two_Mathematically_Equivalent_Versions_of_Maxwells_Equations]], the corresponding proper-time current density is `\rho\mathbf{u} = (b/c)\mathbf{J}`, with `b² = c² + u²` and `\mathbf{u}` the proper velocity of the fluid element.

#### (a) Classical solution — Gaussian (CGS)

In Gaussian units, the Biot–Savart law gives the magnetic field of a current element `I\,d\boldsymbol\ell` as `d\mathbf{B} = (I/c)\,d\boldsymbol\ell \times \hat r / r^{2}`. Parameterising the loop by the azimuthal angle `\phi \in [0, 2\pi]`, the loop element is `d\boldsymbol\ell = (-R\sin\phi, R\cos\phi, 0)\,d\phi`, the displacement from the loop to the field point is `\mathbf{r} = (-R\cos\phi, -R\sin\phi, z)`, and `|\mathbf{r}| = \sqrt{R^{2} + z^{2}}`. The cross product evaluates to

$$
d\boldsymbol\ell \times \mathbf{r} = (R\,z\cos\phi,\; R\,z\sin\phi,\; R^{2})\,d\phi.
$$

The off-axis components integrate to zero by azimuthal symmetry. The axial component integrates to `2\pi R^{2}`, and so the total field on the axis is

$$
\mathbf{B}(0, 0, z) = \frac{I}{c}\,\frac{2\pi R^{2}}{(R^{2} + z^{2})^{3/2}}\,\hat z = \frac{2\pi I R^{2}}{c\,(R^{2} + z^{2})^{3/2}}\,\hat z.
$$

**Mathematica check** (Wolfram MCP, 2026-05-24):

```mathematica
ClearAll[phi, capR, capZ];
dlCrossR = {capR capZ Cos[phi], capR capZ Sin[phi], capR^2};
loopIntegral = Integrate[dlCrossR, {phi, 0, 2 Pi}];
(* loopIntegral = {0, 0, 2 Pi capR^2}  ✅ *)
```

The transverse components vanish identically; the axial component yields the prefactor `2\pi R^{2}` quoted above. Dividing by `(R^{2} + z^{2})^{3/2}` and multiplying by `I/c` reproduces the textbook result, recorded as Jackson 3e Eq. (5.41) and in equivalent form in 2e.

#### (b) Classical solution — SI

In SI units, the Biot–Savart law carries an additional `\mu_{0}/(4\pi)` factor rather than `1/c`. The geometric content of (a) is unchanged, and one obtains

$$
\mathbf{B}(0, 0, z) = \frac{\mu_{0} I R^{2}}{2\,(R^{2} + z^{2})^{3/2}}\,\hat z.
$$

The translation between (a) and (b) is the standard `1/c \to \mu_{0}/(4\pi)` substitution for the Biot–Savart prefactor, equivalent to `4\pi/c \to \mu_{0}` for the prefactor of Ampère's law and consistent with `\varepsilon_{0}\mu_{0} c^{2} = 1` of the SI system.

#### (c) Proper-time reformulation

The configuration is time-independent in the observer's frame, so the proper-time derivative `\partial_\tau` of any field quantity vanishes. The proper-time Ampère–Maxwell law of Eq. (3′) of [[Two_Mathematically_Equivalent_Versions_of_Maxwells_Equations]] therefore reduces to

$$
\nabla \times \mathbf{B} = \frac{4\pi}{b}\,\rho\mathbf{u}.
$$

The source term `(4\pi/b)\rho\mathbf{u}` looks superficially different from the classical Ampère source `(4\pi/c)\mathbf{J}`. We observe that the two are in fact identical, by direct application of the velocity duality `\mathbf{w}/c = \mathbf{u}/b` and the identification `\mathbf{J} = \rho\mathbf{w}`: from these two facts it follows that `\rho\mathbf{u} = (b/c)\mathbf{J}`, and hence

$$
\frac{4\pi}{b}\,\rho\mathbf{u} = \frac{4\pi}{b}\cdot\frac{b}{c}\,\mathbf{J} = \frac{4\pi}{c}\,\mathbf{J}.
$$

**Mathematica check** (Wolfram MCP, 2026-05-24):

```mathematica
ClearAll[bb, cc, capJ];
ampMaxProperTime = (4 Pi/bb) (bb/cc) capJ;
ampMaxClassical = (4 Pi/cc) capJ;
FullSimplify[ampMaxProperTime - ampMaxClassical]
(* Result: 0  ✅ *)
```

The proper-time and classical source terms agree identically for any steady current. The `(b/c)` rescaling of the current density is exactly cancelled by the `1/b` prefactor of Ampère's law in proper-time form. We emphasise that this is an *exact* algebraic cancellation, not an approximate one valid only in the non-relativistic limit `u \ll c`; the proper-time and classical formulations of magnetostatics produce identical field equations and, by uniqueness of solutions to the curl/divergence boundary-value problem, identical magnetic fields.

The proper-time prediction for the field of the current loop is therefore identical to the classical Gaussian result of (a):

$$
\mathbf{B}_{\text{proper-time}}(0, 0, z) = \frac{2\pi I R^{2}}{c\,(R^{2} + z^{2})^{3/2}}\,\hat z.
$$

It is worth observing that the cancellation found here propagates to *every* time-independent magnetostatic problem in Jackson Chs. 5–8: as long as `\partial_\tau \mathbf{E} = 0` (no displacement current in proper-time form) and the current is steady, the proper-time and classical formulations are operationally equivalent at the level of the field equations. The two formulations are mathematically equivalent and physically equivalent in this case, a stronger statement than the general "mathematically equivalent but not physically equivalent" relation of the Gill–Zachary paper, because the time-independence of the configuration closes the gap between the two clocks just as it did for the electrostatic cases of [Problem J3e-P1.5](Ch01_Introduction_Electrostatics.md#problem-j3e-p15--electrostatic-self-energy-of-a-uniformly-charged-sphere) and [Problem J3e-P2.1](Ch02_Boundary_Value_Problems_Electrostatics_I.md#problem-j3e-p21--point-charge-near-a-grounded-conducting-plane).

**Comparison:**

| Quantity | Classical (CGS) | Classical (SI) | Proper-time |
|---|---|---|---|
| `B_z(0, 0, z)` | `\dfrac{2\pi I R^{2}}{c\,(R^{2} + z^{2})^{3/2}}` | `\dfrac{\mu_{0} I R^{2}}{2\,(R^{2} + z^{2})^{3/2}}` | identical to CGS |

**Does the proper-time answer differ from a pure `c → b` redressing?** ✅ no. The naïve `c → b` redressing (which would predict `B_z \propto 1/b` rather than `1/c`) is *not* what the proper-time formulation actually gives; the cancellation between the `(b/c)` current-density rescaling and the `1/b` prefactor of Ampère's law restores the classical `1/c`. This is a sharper conclusion than that of the static-case problems, where no rescaling was active at all.

**Verdict:** ✅ all three solutions consistent. The proper-time formulation reduces to the classical formulation for any time-independent magnetostatic problem, by the exact algebraic cancellation documented above.

**Notes for author review:** the result that proper-time and classical magnetostatics produce identical field equations is a structural feature of the Gill–Zachary substitution rules, not a fortuitous coincidence of the present problem. As the campaign progresses into Chs. 5–8, every steady-current problem will exhibit the same cancellation; we record it once here and cite this entry as the canonical reference when subsequent problems engage the rescaling rule.

**Companion notebook:** [`Roadmapping/Mathematica_Notebooks/Electromagnetism/JacksonCh05_P5_4.wl`](../../Mathematica_Notebooks/Electromagnetism/JacksonCh05_P5_4.wl) — runnable independent of the Mathematica MCP.
