# Ch. 5 — Magnetostatics, Faraday's Law, Quasi-Static Fields

This chapter contains Jackson canonical problems on magnetostatics, worked in the proper-time reformulation alongside their classical CGS and SI solutions. Per [§4 of the plan](../../../.dev/tasks/42-electromagnetism-jackson-proper-time.md#4-unit-system-handling), Ch. 5 of Jackson is presented in the three-system regime.

The Ch. 5 problems below are part of PR 0 (the fluency-warm-up prequel — see [§7.3 of the plan](../../../.dev/tasks/42-electromagnetism-jackson-proper-time.md#73-pr-a-prequel-adopted)). They differ from the Ch. 1 and Ch. 2 problems in that the source carries a steady current, so the proper-time formulation engages the current-density rescaling rule `J → (b/c)\,J = \rho\mathbf{u}`. As we shall see, this rescaling cancels exactly with the `1/b` prefactor of the proper-time Ampère–Maxwell law, leaving the magnetic field identical to the classical result. The fluency lesson is therefore not that the rescaling produces a small correction, but that the rescaling is *algebraically* annihilated for any time-independent current — a stronger statement than the static-case lesson of Chs. 1–2.

## Problems in this chapter

| Problem | Status | Role |
|---|---|---|
| [Problem J3e-P5.4 — Magnetic field of a circular current loop on its axis](#problem-j3e-p54--magnetic-field-of-a-circular-current-loop-on-its-axis) | drafted (PR 0) | fluency-builder |
| [Problem J3e-P5.13 — Magnetic dipole moment of a uniformly rotating charged sphere](#problem-j3e-p513--magnetic-dipole-moment-of-a-uniformly-rotating-charged-sphere) | drafted (PR 0) | fluency-builder (optional 4th) |

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

---

### Problem J3e-P5.13 — Magnetic dipole moment of a uniformly rotating charged sphere

**Selection provenance** (Crocco §5 substantive-AI note):
- *Chosen because:* this is the simplest non-trivial extension of [Problem J3e-P5.4](#problem-j3e-p54--magnetic-field-of-a-circular-current-loop-on-its-axis) — the current distribution now has full three-dimensional structure (a continuous distribution of azimuthal currents at every depth and latitude of the sphere) rather than the planar circular geometry of the current loop. The magnetic dipole moment of a rigidly rotating uniformly charged sphere is a textbook result, treated in [[jackson1998_classical_electrodynamics]] §5.6 and recorded in equivalent form in [[jackson1975_classical_electrodynamics]].
- *Alternatives considered:* the surface-charged rotating sphere (also Jackson §5, slightly easier because the current is confined to a single shell and reduces to a uniformly magnetised sphere by an exact equivalence), and the field on the axis inside the rotating sphere (a substantially more involved integral). The solid-charged case for the dipole moment alone is the optimal fluency target for PR 0: rich enough to exercise the spherical-coordinate machinery, simple enough to land in a single document.
- *Role in this PR:* fluency-builder (optional fourth). Together with [Problem J3e-P5.4](#problem-j3e-p54--magnetic-field-of-a-circular-current-loop-on-its-axis), it confirms that the proper-time and classical formulations of magnetostatics agree on the dipole moment of a steady current distribution.

<!-- TODO: human reviews and fills in — confirms that the dipole-moment treatment (rather than the full field-inside computation) is the appropriate PR 0 scope before this problem is rolled into PR 0's commit -->

**Source:** Jackson, *Classical Electrodynamics*, 3e Problem 5.13 (and 2e Problem 5.13, equivalent). *Paraphrased; consult the textbook for the precise statement.*

**Paraphrased statement:** A solid sphere of radius `R` is uniformly charged with total charge `Q` (volume charge density `\rho = 3Q/(4\pi R^{3})`) and rotates rigidly about a diameter with angular velocity `\boldsymbol\omega`. Compute the magnetic dipole moment of the resulting steady current distribution.

**Setup:** Place the centre of the sphere at the origin and align `\boldsymbol\omega = \omega\hat z`. The current density at a point `\mathbf{r}` inside the sphere is `\mathbf{J}(\mathbf{r}) = \rho\,(\boldsymbol\omega\times\mathbf{r})`. The configuration is time-independent in the observer's frame; as in [Problem J3e-P5.4](#problem-j3e-p54--magnetic-field-of-a-circular-current-loop-on-its-axis), `\partial_\tau = 0` for both fields and current in the proper-time formulation, and the only proper-time machinery that engages is the `\rho\mathbf{u} = (b/c)\mathbf{J}` rescaling rule.

#### (a) Classical solution — Gaussian (CGS)

The magnetic dipole moment of a current distribution in Gaussian units is

$$
\mathbf{m} = \frac{1}{2c}\int \mathbf{r} \times \mathbf{J}(\mathbf{r})\, dV.
$$

Substituting `\mathbf{J} = \rho\,(\boldsymbol\omega \times \mathbf{r})` and applying the vector identity `\mathbf{r} \times (\boldsymbol\omega \times \mathbf{r}) = \boldsymbol\omega (\mathbf{r}\cdot\mathbf{r}) - \mathbf{r}(\mathbf{r}\cdot\boldsymbol\omega)`,

$$
\mathbf{m} = \frac{\rho}{2c}\int\!\left[\boldsymbol\omega\,r^{2} - \mathbf{r}(\mathbf{r}\cdot\boldsymbol\omega)\right] dV.
$$

By the azimuthal symmetry of the configuration about the rotation axis, only the component of `\mathbf{m}` along `\hat z` survives. With `\boldsymbol\omega = \omega\hat z` and `\mathbf{r}\cdot\boldsymbol\omega = \omega r\cos\theta` in spherical coordinates,

$$
m_{z} = \frac{\rho\omega}{2c}\int\!\left[r^{2} - r^{2}\cos^{2}\theta\right] dV = \frac{\rho\omega}{2c}\int r^{2}\sin^{2}\theta\, dV.
$$

We evaluate the integral over the volume of the sphere:

$$
\int r^{2}\sin^{2}\theta\, dV = \int_{0}^{R}\!\!\int_{0}^{\pi}\!\!\int_{0}^{2\pi} r^{2}\sin^{2}\theta\,\cdot r^{2}\sin\theta\, dr\, d\theta\, d\phi = 2\pi \cdot \frac{R^{5}}{5} \cdot \frac{4}{3} = \frac{8\pi R^{5}}{15}.
$$

Substituting back and using `\rho = 3Q/(4\pi R^{3})`,

$$
m_{z} = \frac{\rho\omega}{2c}\cdot\frac{8\pi R^{5}}{15} = \frac{4\pi\rho\omega R^{5}}{15\,c} = \frac{Q\omega R^{2}}{5\,c}.
$$

**Mathematica check** (Wolfram MCP, 2026-05-24):

```mathematica
ClearAll[r, theta, phi, capR, capRho, omega, capQ];
mZ = (capRho omega/2) Integrate[
   r^2 Sin[theta]^2 r^2 Sin[theta],
   {r, 0, capR}, {theta, 0, Pi}, {phi, 0, 2 Pi},
   Assumptions -> capR > 0
   ];
mZSimplified = FullSimplify[mZ /. capRho -> 3 capQ/(4 Pi capR^3)];
(* Result: capQ omega capR^2 / 5  ✅  (SI form; Gaussian carries an extra 1/c) *)
```

We observe that the result is independent of the carrier mass — the magnetic dipole moment depends only on the charge distribution and the angular velocity, not on the dynamical properties of the material carrying the current.

#### (b) Classical solution — SI

The magnetic dipole moment of a current distribution in SI units is

$$
\mathbf{m} = \frac{1}{2}\int \mathbf{r} \times \mathbf{J}(\mathbf{r})\, dV,
$$

with no `1/c` factor. The integrals of (a) are unchanged, and we obtain

$$
m_{z} = \frac{Q\omega R^{2}}{5}.
$$

The translation between (a) and (b) is the single substitution `m_{\text{Gaussian}} = m_{\text{SI}}/c`.

#### (c) Proper-time reformulation

The configuration is time-independent in the observer's frame, so `\partial_\tau` of every field quantity vanishes. The proper-time current density is `\rho\mathbf{u} = (b/c)\mathbf{J}` by the velocity duality `\mathbf{w}/c = \mathbf{u}/b` and the identification `\mathbf{J} = \rho\mathbf{w}`, as recorded in [`_proper_time_cheatsheet.md`](../_proper_time_cheatsheet.md).

The magnetic dipole moment in the proper-time formulation is defined by the analogous integral over `\rho\mathbf{u}` rather than `\mathbf{J}/c`. Writing `\mathbf{m}_{\text{proper-time}} = (1/2)\int \mathbf{r} \times \rho\mathbf{u}\, dV` and substituting `\rho\mathbf{u} = (b/c)\mathbf{J}`,

$$
\mathbf{m}_{\text{proper-time}} = \frac{1}{2}\int \mathbf{r} \times \frac{b}{c}\mathbf{J}\, dV = \frac{b}{2c}\int \mathbf{r} \times \mathbf{J}\, dV = \frac{b}{c}\,\mathbf{m}_{\text{SI}}\bigg|_{c\text{-factor restored}}.
$$

Note that this is *not* the same as the Gaussian moment, because the prefactor is `b/c`, not `1/c`. The proper-time moment differs from the Gaussian moment by a factor of `b/1 = b`, that is, by exactly the factor `\sqrt{c^{2} + u^{2}}` divided by `c`.

We must therefore ask: what is `u` for the rotating sphere? The proper velocity `u` of a fluid element is the magnitude of `du/d\tau` for the world-line of that element, which for a rigid rotation is `u = \omega R \sin\theta\,(b/c)` — implicit in `b`. Solving `u^{2} = c^{2}(b/c)^{2} - c^{2}` is straightforward in the non-relativistic limit `\omega R \ll c`, where one obtains `u \approx \omega r \sin\theta` and `b \approx c[1 + (1/2)(\omega r\sin\theta/c)^{2} + \cdots]`. The leading-order correction to the dipole moment is therefore `O(u^{2}/c^{2}) = O((\omega R/c)^{2})`, which is far below the observational floor for any ordinary rotating charged sphere.

In the strict limit `\omega R/c \to 0`, the `b/c` factor in the proper-time moment is unity and the proper-time and classical moments agree:

$$
\mathbf{m}_{\text{proper-time}} = \mathbf{m}_{\text{Gaussian}}\cdot c = \mathbf{m}_{\text{SI}} = \frac{Q\omega R^{2}}{5}\,\hat z\quad\text{(non-relativistic limit)}.
$$

This is consistent with the structural cancellation noted in [Problem J3e-P5.4](#problem-j3e-p54--magnetic-field-of-a-circular-current-loop-on-its-axis) at the level of the field equations: as long as the source velocity is much less than `c`, proper-time and classical magnetostatics produce the same observables to leading order, and a small `O((u/c)^{2})` correction at the next order. The correction is the proper-time analogue of the standard special-relativistic correction to magnetic moments, expressed in the proper-time variables rather than the observer-time ones.

**Comparison:**

| Quantity | Classical (CGS) | Classical (SI) | Proper-time (non-relativistic limit) |
|---|---|---|---|
| `m_z` | `Q\omega R^{2}/(5c)` | `Q\omega R^{2}/5` | `Q\omega R^{2}/5` (matches SI to leading order; `O((u/c)^{2})` correction at next order) |

**Does the proper-time answer differ from a pure `c → b` redressing?** ⚠ yes, at second order in `u/c`. The strict `c → b` redressing of the Gaussian result would replace `1/c` by `1/b`, giving an apparent `(c/b)`-suppression of the dipole moment. The proper-time formulation introduces a *separate* `b/c` factor from the current-density rescaling, and the two factors combine to leave the moment equal (to leading order) to the Gaussian one with no `b`-dependence at all. This is the same cancellation noted in [Problem J3e-P5.4](#problem-j3e-p54--magnetic-field-of-a-circular-current-loop-on-its-axis), but applied to a derived quantity (the dipole moment) rather than directly to the field equation. The residual `O((u/c)^{2})` correction reflects the fact that `b` itself depends on position through `u(\mathbf{r}) = \omega r\sin\theta\,(b/c)`, and the position-dependence is not captured by a simple constant prefactor.

**Verdict:** ✅ all three solutions consistent at leading order, with a documented `O((u/c)^{2})` proper-time correction that is below the observational floor for non-relativistic configurations. The campaign records this entry as the first per-problem document in which the proper-time formulation produces a non-zero deviation from the classical answer, however small.

**Notes for author review:** the leading-order proper-time correction to the magnetic dipole moment of a rotating charged sphere is, to my knowledge, not recorded in the published Gill–Zachary corpus or in any standard textbook treatment. It is mechanically derivable from the velocity-duality rule alone and does not require Eq. (24) or any other flagged finding. The magnitude `(\omega R/c)^{2}` is far below the precision of any classical-mechanics-scale experiment on a rotating macroscopic charged sphere. For a Penning-trap electron in cyclotron motion with `\omega R/c \sim 10^{-2}`, the analogous correction would be `\sim 10^{-4}`, which is large compared to the `10^{-13}` precision of g-2 measurements; however, the relevant geometry (single electron, not a rotating macroscopic charge distribution) is sufficiently different that the present result does not transfer directly. This is a candidate observation for inclusion in `Roadmapping/Equation_Verification/FINDINGS_for_author_review.md` if PR D or later work confirms a similar `(u/c)^{2}` correction in the radiation-reaction regime.

**Companion notebook:** [`Roadmapping/Mathematica_Notebooks/Electromagnetism/JacksonCh05_P5_13.wl`](../../Mathematica_Notebooks/Electromagnetism/JacksonCh05_P5_13.wl) — runnable independent of the Mathematica MCP.
