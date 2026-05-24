# Ch. 6 — Maxwell Equations, Macroscopic Electromagnetism, Conservation Laws

This chapter contains Jackson canonical problems on Maxwell's equations and their conservation laws, worked in the proper-time reformulation alongside their classical CGS and SI solutions. Per [§4 of the plan](../../../.dev/tasks/42-electromagnetism-jackson-proper-time.md#4-unit-system-handling), Ch. 6 of Jackson is presented in the three-system regime.

Ch. 6 is the campaign's first chapter in which the proper-time machinery engages on every problem. The configurations are no longer static or purely-magnetostatic; time-varying fields, displacement currents, and (in [Problem J3e-P6.20](#problem-j3e-p620--radiation-pressure-on-a-perfect-conductor)) the dissipative `−(\mathbf{u}\cdot\mathbf{a})/b^{4}` coefficient of Eq. (4) of [[Two_Mathematically_Equivalent_Versions_of_Maxwells_Equations]] all enter the discussion. Per [§13.5 D2 of the plan](../../../.dev/tasks/42-electromagnetism-jackson-proper-time.md#135-decision-points--confirmed-by-author-2026-05-24), per-paragraph `<!-- TODO: human reviews and fills in -->` blocks accompany each substantive interpretive claim throughout this chapter.

## Problems in this chapter

| Problem | Status | Role |
|---|---|---|
| [Problem J3e-P6.1 — Maxwell with magnetic monopoles + electric–magnetic duality](#problem-j3e-p61--maxwell-with-magnetic-monopoles-and-electricmagnetic-duality) | drafted (PR A) | headline-adjacent |
| [Problem J3e-P6.4 — EM momentum of a uniformly-moving point charge](#problem-j3e-p64--em-momentum-of-a-uniformly-moving-point-charge) | drafted (PR A) | headline-payoff (podcast pick #1) |
| [Problem J3e-P6.5 — Poynting theorem in macroscopic media](#problem-j3e-p65--poynting-theorem-in-macroscopic-media) | planned (PR A) | headline-adjacent |
| [Problem J3e-P6.11 — Symmetric stress tensor and Lorentz behaviour](#problem-j3e-p611--symmetric-stress-tensor-and-lorentz-behaviour) | planned (PR A) | headline-adjacent |
| [Problem J3e-P6.20 — Radiation pressure on a perfect conductor](#problem-j3e-p620--radiation-pressure-on-a-perfect-conductor) | planned (PR A) | headline-payoff (first `(u·a)/b⁴` test) |

---

### Problem J3e-P6.1 — Maxwell with magnetic monopoles and electric–magnetic duality

**Selection provenance** (Crocco §5 substantive-AI note):
- *Chosen because:* the magnetic-monopole extension of Maxwell's equations is a textbook test of the structural symmetries of the theory, treated in [[jackson1998_classical_electrodynamics]] §6.11 and §6.12 and recorded in equivalent form in [[jackson1975_classical_electrodynamics]]. It is the natural opening problem for Ch. 6 because it exercises *all four* Maxwell equations on equal footing and forces the campaign to confront how the proper-time substitution rules act on the source terms in both the curl and divergence equations.
- *Alternatives considered:* J3e-P6.2 (gauge transformations and the Lorenz gauge — interesting but operates only on the potentials, not the field equations) and J3e-P6.3 (energy in the Lorenz-gauge formulation — too closely related to J3e-P6.5 to merit a separate entry).
- *Role in this PR:* headline-adjacent. The classical derivation is short, but the proper-time analysis surfaces a structural observation about how the framework's `b² = c² + u²` depends on the choice of source species — a non-trivial issue when both electric and magnetic monopoles are present.

<!-- TODO: human reviews and fills in — confirms the choice of monopoles as PR A's opening problem (rather than J3e-P6.2 gauge transformations) before this problem is rolled into the PR -->

**Source:** Jackson, *Classical Electrodynamics*, 3e Problem 6.1 (and 2e Problem 6.1, equivalent). *Paraphrased; consult the textbook for the precise statement.*

**Paraphrased statement:** Generalise Maxwell's equations to admit magnetic-monopole sources (a magnetic charge density `\rho_{m}` and a magnetic current density `\mathbf{J}_{m}`). Demonstrate that the resulting set of equations is invariant under the electric–magnetic duality transformation `\mathbf{E} \to \mathbf{B}`, `\mathbf{B} \to -\mathbf{E}`, `\rho_{e} \to \rho_{m}`, `\rho_{m} \to -\rho_{e}`, `\mathbf{J}_{e} \to \mathbf{J}_{m}`, `\mathbf{J}_{m} \to -\mathbf{J}_{e}`. Carry out the analysis in CGS, in SI, and in the proper-time reformulation.

**Setup:** Following [[jackson1998_classical_electrodynamics]] §6.11, we adopt the symmetric form of Maxwell's equations in which a magnetic-charge density and a magnetic-current density appear on the right-hand sides of `\nabla\cdot\mathbf{B}` and `\nabla\times\mathbf{E}` respectively, mirroring the electric sources on `\nabla\cdot\mathbf{E}` and `\nabla\times\mathbf{B}`. The electric monopoles have charge density `\rho_{e}` and current density `\mathbf{J}_{e} = \rho_{e}\mathbf{w}_{e}` (with `\mathbf{w}_{e}` their observer-time velocity); similarly `\rho_{m}` and `\mathbf{J}_{m} = \rho_{m}\mathbf{w}_{m}` for the magnetic monopoles. In the proper-time formulation, the corresponding proper velocities are `\mathbf{u}_{e}` and `\mathbf{u}_{m}`, with `b_{e}^{2} = c^{2} + \mathbf{u}_{e}^{2}` and `b_{m}^{2} = c^{2} + \mathbf{u}_{m}^{2}`.

#### (a) Classical solution — Gaussian (CGS)

The Maxwell equations extended with magnetic-monopole sources, in Gaussian units, read

$$
\begin{aligned}
\nabla\cdot\mathbf{E} &= 4\pi\rho_{e},\\
\nabla\cdot\mathbf{B} &= 4\pi\rho_{m},\\
\nabla\times\mathbf{E} &= -\frac{1}{c}\frac{\partial\mathbf{B}}{\partial t} - \frac{4\pi}{c}\mathbf{J}_{m},\\
\nabla\times\mathbf{B} &= \frac{1}{c}\frac{\partial\mathbf{E}}{\partial t} + \frac{4\pi}{c}\mathbf{J}_{e}.
\end{aligned}
$$

The Faraday law (third equation) acquires a magnetic-current source `−(4\pi/c)\mathbf{J}_{m}` in symmetry with the electric-current source of the Ampère–Maxwell law (fourth equation); the divergence equations carry the charge-density sources `4\pi\rho_{e}` and `4\pi\rho_{m}` respectively.

It is easy to show that this set is invariant under the duality transformation. Applying `\mathbf{E}\to\mathbf{B}`, `\mathbf{B}\to-\mathbf{E}`, `\rho_{e}\to\rho_{m}`, `\rho_{m}\to-\rho_{e}`, `\mathbf{J}_{e}\to\mathbf{J}_{m}`, `\mathbf{J}_{m}\to-\mathbf{J}_{e}` to the first equation gives `\nabla\cdot\mathbf{B} = 4\pi\rho_{m}`, which is the second equation; conversely the second yields `\nabla\cdot(-\mathbf{E}) = 4\pi(-\rho_{e})`, equivalent to the first. The third equation transforms to `\nabla\times\mathbf{B} = -(1/c)\partial(-\mathbf{E})/\partial t - (4\pi/c)(-\mathbf{J}_{e}) = (1/c)\partial\mathbf{E}/\partial t + (4\pi/c)\mathbf{J}_{e}`, which is the fourth; and the fourth transforms to the third. The four equations therefore form a closed orbit of size two under the duality, and the transformation is a symmetry of the field equations.

The duality is, as recorded by Jackson 3e §6.11, the simplest of a one-parameter family of continuous rotations of the `(\mathbf{E},\mathbf{B})` doublet, with the discrete `\pi/2`-rotation above as one limit. We restrict attention to the discrete case for the present problem.

<!-- TODO: human reviews and fills in — confirms the framing of the duality as a "structural symmetry of the field equations" rather than a deeper claim about physical equivalence between electric and magnetic charges -->

#### (b) Classical solution — SI

In SI units, Maxwell's equations with magnetic-monopole sources read

$$
\begin{aligned}
\nabla\cdot\mathbf{E} &= \frac{\rho_{e}}{\varepsilon_{0}},\\
\nabla\cdot\mathbf{B} &= \mu_{0}\rho_{m},\\
\nabla\times\mathbf{E} &= -\frac{\partial\mathbf{B}}{\partial t} - \mu_{0}\mathbf{J}_{m},\\
\nabla\times\mathbf{B} &= \mu_{0}\varepsilon_{0}\frac{\partial\mathbf{E}}{\partial t} + \mu_{0}\mathbf{J}_{e}.
\end{aligned}
$$

The duality transformation in SI is conventionally written `\mathbf{E}\to c\mathbf{B}`, `c\mathbf{B}\to-\mathbf{E}`, `\rho_{e}\to\rho_{m}/c`, `c\rho_{m}\to-\rho_{e}`, with similar `c`-factors for the currents. The `c`-factors arise because `\mathbf{E}` and `c\mathbf{B}` share a dimension in SI but not the same numerical magnitude, so the duality cannot exchange them without an explicit unit-conversion factor. We observe that this dimensional asymmetry, absent in Gaussian units, is one of the reasons that Jackson 3e Chs. 11+ revert to Gaussian when treating field-theoretic dualities and Lorentz covariance.

<!-- TODO: human reviews and fills in — confirms the framing of the c-factors as "unit-system bookkeeping that is asymmetric in SI" rather than a deeper claim about field-theoretic structure -->

#### (c) Proper-time reformulation

The proper-time Maxwell equations of Eq. (3′) of [[Two_Mathematically_Equivalent_Versions_of_Maxwells_Equations]], extended with magnetic-monopole sources, read

$$
\begin{aligned}
\nabla\cdot\mathbf{E} &= 4\pi\rho_{e},\\
\nabla\cdot\mathbf{B} &= 4\pi\rho_{m},\\
\nabla\times\mathbf{E} &= -\frac{1}{b_{e}}\frac{\partial\mathbf{B}}{\partial\tau_{e}} - \frac{4\pi}{b_{m}}\rho_{m}\mathbf{u}_{m},\\
\nabla\times\mathbf{B} &= \frac{1}{b_{e}}\frac{\partial\mathbf{E}}{\partial\tau_{e}} + \frac{4\pi}{b_{e}}\rho_{e}\mathbf{u}_{e}.
\end{aligned}
$$

We have written `b_{e}` and `\tau_{e}` for the proper-time variables associated with the electric source in the Ampère–Maxwell equation, and `b_{m}` for the variable associated with the magnetic source in the Faraday equation. **The question of which proper-time variable to use in the displacement-current term of the Faraday equation does not have a unique answer in the framework as published**, because the Gill–Zachary paper formulates the proper-time substitution for a single source species. In the single-species case, this ambiguity does not arise; the magnetic-monopole extension is the first place in the campaign where it becomes operationally relevant.

<!-- TODO: human reviews and fills in — flags the multi-species ambiguity as a structural observation about the framework's domain of validity; needs author judgement on whether this is an extension worth pursuing in DRQM-II or a non-issue because the displacement-current term is unique in the lab frame anyway -->

The source-term reductions to the classical formulation are exact, by the same algebraic cancellation noted in [Problem J3e-P5.4](Ch05_Magnetostatics_Faraday_Quasi_Static.md#problem-j3e-p54--magnetic-field-of-a-circular-current-loop-on-its-axis): using `\rho_{e}\mathbf{u}_{e} = (b_{e}/c)\mathbf{J}_{e}` and `\rho_{m}\mathbf{u}_{m} = (b_{m}/c)\mathbf{J}_{m}`,

$$
\frac{4\pi}{b_{e}}\rho_{e}\mathbf{u}_{e} = \frac{4\pi}{c}\mathbf{J}_{e}, \qquad \frac{4\pi}{b_{m}}\rho_{m}\mathbf{u}_{m} = \frac{4\pi}{c}\mathbf{J}_{m}.
$$

**Mathematica check** (Wolfram MCP, 2026-05-24):

```mathematica
ClearAll[bb, cc, capJe, capJm];
electricSource = (4 Pi/bb) (bb/cc) capJe;
electricSourceClassical = (4 Pi/cc) capJe;
magneticSource = (4 Pi/bb) (bb/cc) capJm;
magneticSourceClassical = (4 Pi/cc) capJm;
FullSimplify[electricSource - electricSourceClassical]   (* 0  ✅ *)
FullSimplify[magneticSource - magneticSourceClassical]   (* 0  ✅ *)
```

Both source-term reductions vanish identically. The cancellation found at the field-equation level in [Problem J3e-P5.4](Ch05_Magnetostatics_Faraday_Quasi_Static.md#problem-j3e-p54--magnetic-field-of-a-circular-current-loop-on-its-axis) for electric sources extends symmetrically to magnetic sources, as the duality between the two species would require.

It follows that, restricting to either pure electric monopoles (`\rho_{m} = 0`, `\mathbf{J}_{m} = 0`) or pure magnetic monopoles (`\rho_{e} = 0`, `\mathbf{J}_{e} = 0`), the proper-time formulation reduces to the classical formulation by the same exact-cancellation mechanism as in PR 0. The duality transformation acts on the *classical* equations and on the *proper-time* equations symmetrically; the proper-time formulation does not break the duality symmetry of the classical theory.

<!-- TODO: human reviews and fills in — confirms the framing that "the duality symmetry is preserved by the proper-time substitution" is an algebraic claim about the substitution rules, not a deeper physical statement about whether magnetic monopoles exist or are observable -->

**Comparison:**

| Quantity | Classical (CGS) | Classical (SI) | Proper-time |
|---|---|---|---|
| `\nabla\cdot\mathbf{E}` source | `4\pi\rho_{e}` | `\rho_{e}/\varepsilon_{0}` | `4\pi\rho_{e}` |
| `\nabla\cdot\mathbf{B}` source | `4\pi\rho_{m}` | `\mu_{0}\rho_{m}` | `4\pi\rho_{m}` |
| `\nabla\times\mathbf{E}` source term | `-(4\pi/c)\mathbf{J}_{m}` | `-\mu_{0}\mathbf{J}_{m}` | `-(4\pi/b_{m})\rho_{m}\mathbf{u}_{m} = -(4\pi/c)\mathbf{J}_{m}` |
| `\nabla\times\mathbf{B}` source term | `(4\pi/c)\mathbf{J}_{e}` | `\mu_{0}\mathbf{J}_{e}` | `(4\pi/b_{e})\rho_{e}\mathbf{u}_{e} = (4\pi/c)\mathbf{J}_{e}` |
| Duality symmetry under `\pi/2` rotation | preserved | preserved (with `c`-factor bookkeeping) | preserved (single-species); structurally ambiguous for mixed sources |

**Does the proper-time answer differ from a pure `c → b` redressing?** ✅ no, for the source-term content of each equation individually. The proper-time formulation produces the same source structure as the classical one for any single-species sub-problem. The displacement-current `(1/b)\partial_\tau` term has the same content as classical `(1/c)\partial_t` under the time-derivative duality of Eq. (2) of the Maxwell paper.

**Verdict:** ✅ all three solutions are consistent at the level of the field equations. ⚠ A structural observation: the proper-time framework as published does not specify how to handle the displacement-current term when both electric and magnetic monopoles are present with distinct proper velocities. The ambiguity does not affect the present problem's verification of the duality symmetry, but it is the first place in the campaign where the framework's single-species formulation becomes operationally visible.

**Notes for author review:** the multi-species ambiguity in the proper-time formulation (which `\tau` governs the displacement-current term when both electric and magnetic monopoles are present?) is a structural observation about the framework's domain of validity, not a flagged inconsistency. It does not bear on the single-species derivations of [[Two_Mathematically_Equivalent_Versions_of_Maxwells_Equations]] Eqs. (1)–(23), all of which are verified. It is worth recording as a candidate question for any future Gill–Zachary follow-up paper that addresses multi-species sources (the kinetic-theory regime, plasma physics with both electrons and positrons treated separately, etc.). I have flagged it in this document but not posted it to [`FINDINGS_for_author_review.md`](../../Equation_Verification/FINDINGS_for_author_review.md) because it is an extension question, not an inconsistency in the published framework.

**Companion notebook:** [`Roadmapping/Mathematica_Notebooks/Electromagnetism/JacksonCh06_P6_1.wl`](../../Mathematica_Notebooks/Electromagnetism/JacksonCh06_P6_1.wl) — runnable independent of the Mathematica MCP.

---

### Problem J3e-P6.4 — EM momentum of a uniformly-moving point charge

**Selection provenance** (Crocco §5 substantive-AI note):
- *Chosen because:* the electromagnetic momentum of a slowly-moving point charge is the textbook entry to the Abraham–Lorentz "4/3 problem" and to the broader question of where electromagnetic momentum lives — in the field, in the matter, or in some mix. It is recorded as a standard problem in [[jackson1998_classical_electrodynamics]] and is podcast pick #1 in [§12.1 of the campaign plan](../../../.dev/tasks/42-electromagnetism-jackson-proper-time.md#121-per-problem-briefs).
- *Alternatives considered:* J3e-P6.3 (energy and momentum of a Lorenz-gauge formulation — too closely related to J3e-P6.5 to merit separate treatment) and J3e-P12.6 (radiation pressure on a charged particle in a plane wave — a Ch. 12 problem, deferred to PR C).
- *Role in this PR:* headline-payoff. The Abraham–Lorentz puzzle is one of the campaign's load-bearing narrative threads, and the proper-time analysis must address it transparently: the framework reproduces the classical `(4/3)` factor exactly for uniform motion, neither dissolving nor exacerbating the puzzle.

<!-- TODO: human reviews and fills in — confirms the role of this problem as a headline-payoff for PR A, and the framing that the proper-time framework does not dissolve the 4/3 puzzle (vs the framing that it does, which would be a stronger claim) -->

**Source:** Jackson, *Classical Electrodynamics*, 3e Problem 6.4 (and 2e Problem 6.4, equivalent). *Paraphrased; consult the textbook for the precise statement.*

**Paraphrased statement:** A point charge `q` (regularised as a uniformly charged solid sphere of radius `R` to make the field-energy integrals finite) moves with constant velocity `\mathbf{v}` of magnitude `v \ll c` in the observer's frame. Compute the linear momentum carried by the electromagnetic field of the moving charge. Carry out the analysis in CGS, in SI, and in the proper-time reformulation, and remark on the famous "4/3 factor" that the result exhibits.

**Setup:** The charge moves with velocity `\mathbf{v} = v\hat z` in the observer's frame. In its rest frame the charge is a uniformly charged solid sphere of radius `R` carrying total charge `q`, with the rest-frame electric field given by the result of [Problem J3e-P1.5](Ch01_Introduction_Electrostatics.md#problem-j3e-p15--electrostatic-self-energy-of-a-uniformly-charged-sphere). To leading order in `v/c`, the lab-frame electric field is the rest-frame field translated rigidly with the source, and the lab-frame magnetic field follows from the relativistic field transformation `\mathbf{B} = (\mathbf{v}/c)\times\mathbf{E}` (Gaussian). The proper-time velocity of the source fluid element is `\mathbf{u} = (b/c)\mathbf{v}` with `b^{2} = c^{2} + u^{2}`, so `u \to v` and `b \to c` in the non-relativistic limit.

#### (a) Classical solution — Gaussian (CGS)

The electromagnetic momentum density in Gaussian units is `\mathbf{g} = (1/(4\pi c))\,\mathbf{E}\times\mathbf{B}`. Substituting `\mathbf{B} = (\mathbf{v}/c)\times\mathbf{E}` and applying the BAC–CAB identity,

$$
\mathbf{g} = \frac{1}{4\pi c^{2}}\left[\mathbf{v}\,(\mathbf{E}\cdot\mathbf{E}) - \mathbf{E}\,(\mathbf{v}\cdot\mathbf{E})\right].
$$

For `\mathbf{v} = v\hat z`, the `\hat z`-component of `\mathbf{g}` is

$$
g_{z} = \frac{1}{4\pi c^{2}}\left[v\,(E_{x}^{2} + E_{y}^{2} + E_{z}^{2}) - v\,E_{z}^{2}\right] = \frac{v}{4\pi c^{2}}(E_{x}^{2} + E_{y}^{2}).
$$

Integrating over all space and using the spherical-symmetry identity `\int E_{x}^{2}\, dV = \int E_{y}^{2}\, dV = \int E_{z}^{2}\, dV = (1/3)\int E^{2}\, dV`,

$$
P_{z} = \int g_{z}\, dV = \frac{v}{4\pi c^{2}}\cdot\frac{2}{3}\int E^{2}\, dV.
$$

The total field-squared integral was computed in [Problem J3e-P1.5](Ch01_Introduction_Electrostatics.md#problem-j3e-p15--electrostatic-self-energy-of-a-uniformly-charged-sphere); evaluating it on the present configuration gives `\int E^{2}\, dV = 24\pi q^{2}/(5 R)`. Substituting,

$$
P_{z} = \frac{v}{4\pi c^{2}}\cdot\frac{2}{3}\cdot\frac{24\pi q^{2}}{5 R} = \frac{4}{3}\cdot\frac{q^{2}}{5 R c^{2}}\cdot v \cdot 3 = \frac{4\,q^{2}\,v}{5\,R\,c^{2}}.
$$

(The factor of three in the second line is an algebraic cancellation; the result reduces cleanly to `4 q^{2} v / (5 R c^{2})`.)

**Mathematica check** (Wolfram MCP, 2026-05-24):

```mathematica
ClearAll[r, qq, capR, vv, cc];
eSquaredIn = Integrate[(qq r/capR^3)^2 4 Pi r^2, {r, 0, capR}, Assumptions -> capR > 0];
eSquaredOut = Integrate[(qq/r^2)^2 4 Pi r^2, {r, capR, Infinity}, Assumptions -> capR > 0];
totalESquared = FullSimplify[eSquaredIn + eSquaredOut];
pField = (vv/(4 Pi cc^2)) (2/3) totalESquared;
FullSimplify[pField]
(* Result: 4 qq^2 vv / (5 capR cc^2)  ✅ *)
```

Comparing with the field energy `U = (3/5)(q^{2}/R)` from [Problem J3e-P1.5](Ch01_Introduction_Electrostatics.md#problem-j3e-p15--electrostatic-self-energy-of-a-uniformly-charged-sphere), the momentum can also be written

$$
\mathbf{P}_{\text{field}} = \frac{4}{3}\,\frac{U}{c^{2}}\,\mathbf{v}.
$$

We observe that the coefficient is `4/3`, not `1`. This is the celebrated "4/3 problem" of classical electron theory, recorded in detail by Abraham (1903), Lorentz (1904), and Poincaré (1906): if the electromagnetic mass of the electron is identified with `U/c^{2}`, then the momentum of a moving electron exceeds the naïve `m v` by a factor of `4/3`. The puzzle is not that the electromagnetic momentum is wrong — the calculation above is straightforward — but that the identification of `U/c^{2}` as a Lorentz-covariant rest mass is in tension with the momentum scaling. Poincaré resolved the tension by introducing internal stresses (the "Poincaré stresses") that contribute an additional `-(1/3) U/c^{2} \mathbf{v}` of momentum, restoring `\mathbf{P}_{\text{total}} = (U/c^{2})\mathbf{v}` as a Lorentz-covariant relation.

<!-- TODO: human reviews and fills in — confirms the Abraham–Lorentz–Poincaré historical framing and its scope: the puzzle is about the relation between the electromagnetic mass and the Lorentz-covariant momentum, not about whether the (4/3) factor is "correct" in isolation -->

#### (b) Classical solution — SI

The same computation in SI units, using `\mathbf{g} = \varepsilon_{0}\,\mathbf{E}\times\mathbf{B}` and `\mathbf{B} = (1/c^{2})\,\mathbf{v}\times\mathbf{E}` (the SI form of the boost), gives

$$
\mathbf{P}_{\text{field}} = \frac{4}{3}\,\frac{U_{\text{SI}}}{c^{2}}\,\mathbf{v} = \frac{q^{2}\,v}{6\pi\varepsilon_{0}\,R\,c^{2}}\,\hat z,
$$

with `U_{\text{SI}} = q^{2}/(8\pi\varepsilon_{0} R) \cdot (3/5)` per [Problem J3e-P1.5](Ch01_Introduction_Electrostatics.md#problem-j3e-p15--electrostatic-self-energy-of-a-uniformly-charged-sphere). The `4/3` factor is unchanged from CGS to SI, as one expects for a dimensionless geometric coefficient; the unit-system substitution affects only the prefactor `U/c^{2}`.

#### (c) Proper-time reformulation

The configuration is uniform motion in the observer's frame: `\mathbf{w} = \mathbf{v} = v\hat z` is constant, and the proper-time acceleration `\mathbf{a} = d\mathbf{u}/d\tau = 0`. With `\mathbf{u}\cdot\mathbf{a} = 0`, the dissipative coefficient `-(\mathbf{u}\cdot\mathbf{a})/b^{4}` of Eq. (4) of [[Two_Mathematically_Equivalent_Versions_of_Maxwells_Equations]] vanishes; the new third term of the Liénard–Wiechert fields of Eq. (7) likewise vanishes; the modified Lorentz force law of Eq. (18) reduces to `\mathbf{F} = q[\mathbf{E} + (\mathbf{u}/b)\times\mathbf{B}]`, which by the velocity-duality identity `\mathbf{u}/b = \mathbf{w}/c` of Eq. (1) is numerically identical to the classical Lorentz force.

The proper-time Liénard–Wiechert formula for the field of a uniformly-moving charge (the first term of Eq. (7), with the third term suppressed by `\mathbf{u}\cdot\mathbf{a} = 0`) reads

$$
\mathbf{E}(\mathbf{x},\tau) = \frac{e\,\mathbf{r}_{\mathbf{u}}\,(1 - \mathbf{u}^{2}/b^{2})}{s^{3}}.
$$

We observe that the prefactor `(1 - \mathbf{u}^{2}/b^{2})` looks superficially different from the classical Liénard–Wiechert prefactor `(1 - \mathbf{w}^{2}/c^{2})`. It is easy to show that the two are equal under the velocity-duality identity. From `\mathbf{w}/c = \mathbf{u}/b` and `b^{2} = c^{2} + \mathbf{u}^{2}` it follows that `b = c/\sqrt{1 - \mathbf{w}^{2}/c^{2}}` and `\mathbf{u} = \mathbf{w}/\sqrt{1 - \mathbf{w}^{2}/c^{2}}`. Substituting,

$$
1 - \frac{\mathbf{u}^{2}}{b^{2}} = 1 - \frac{\mathbf{w}^{2}/(1 - \mathbf{w}^{2}/c^{2})}{c^{2}/(1 - \mathbf{w}^{2}/c^{2})} = 1 - \frac{\mathbf{w}^{2}}{c^{2}}.
$$

**Mathematica check** (Wolfram MCP, 2026-05-24):

```mathematica
ClearAll[ww, cc];
bSolved = cc/Sqrt[1 - ww^2/cc^2];
uSolved = ww/Sqrt[1 - ww^2/cc^2];
identityCheck = FullSimplify[
   1 - uSolved^2/bSolved^2 - (1 - ww^2/cc^2),
   Assumptions -> 0 < ww < cc
   ];
(* Result: 0  ✅ *)
```

The proper-time Liénard–Wiechert formula is therefore numerically identical to the classical formula for any uniformly-moving charge, expressed in different variables. The electromagnetic field of the moving charge is the same field in either formulation, and the integral that yields `\mathbf{P}_{\text{field}} = (4/3)(U/c^{2})\mathbf{v}` is the same integral.

It follows that the proper-time formulation reproduces the classical `4/3` factor exactly. The framework does *not* dissolve the Abraham–Lorentz puzzle through some new mechanism of EM-momentum bookkeeping; the puzzle's resolution still requires either the Poincaré stresses (or their modern analogue, the covariant stress-energy tensor of a self-consistent classical electron), or the abandonment of the point-charge regularisation in favour of a fully field-theoretic treatment. The proper-time framework is consistent with the classical resolution; it does not provide an alternative one.

<!-- TODO: human reviews and fills in — confirms the framing that the proper-time framework is "consistent with" the Poincaré-stress resolution rather than offering a competing one. This is the load-bearing claim of the (c) section and is the subtlest piece of the per-problem document; needs author judgement on whether to leave it as conditional ("the framework is consistent with, but does not require, the Poincaré stresses") or stronger ("the framework requires the same Poincaré stresses as classical EM") -->

In the language of [§12.1's podcast brief](../../../.dev/tasks/42-electromagnetism-jackson-proper-time.md#121-per-problem-briefs), the proper-time prediction for the field momentum of a uniformly-moving point charge is identical to the classical prediction, with no observable deviation. This is the cleanest demonstration in the campaign that the velocity-duality substitution rule is *consistent* with the classical EM mass calculation, neither weakening nor strengthening the Abraham–Lorentz puzzle. The framework's most distinguishable predictions live in problems where `\mathbf{u}\cdot\mathbf{a} \neq 0` — specifically in [Problem J3e-P6.20](#problem-j3e-p620--radiation-pressure-on-a-perfect-conductor) of this chapter and in PRs D–E.

**Comparison:**

| Quantity | Classical (CGS) | Classical (SI) | Proper-time |
|---|---|---|---|
| `\mathbf{P}_{\text{field}}` | `(4/3)(U/c^{2})\mathbf{v}` | `(4/3)(U_{\text{SI}}/c^{2})\mathbf{v}` | identical to CGS / SI |
| Coefficient | `4/3` | `4/3` | `4/3` |
| Resolution of 4/3 puzzle | Poincaré stresses (1906) | same | same — framework reproduces but does not dissolve the puzzle |

**Does the proper-time answer differ from a pure `c → b` redressing?** ✅ no, for uniform motion. The velocity-duality identity `1 - u^{2}/b^{2} = 1 - w^{2}/c^{2}` makes the proper-time and classical Liénard–Wiechert formulas numerically equivalent for any uniformly-moving charge.

**Verdict:** ✅ all three solutions consistent. The proper-time formulation reproduces the `4/3` factor exactly and is consistent with the Poincaré-stress resolution of the Abraham–Lorentz puzzle. It does not offer an alternative resolution; the puzzle remains where the classical theory left it.

**Notes for author review:** the framing in (c) — that the proper-time framework is *consistent with* but does not *require* the Poincaré stresses — is the subtlest interpretive claim in this document and is flagged with a per-paragraph `<!-- TODO -->`. It would be defensible to make a stronger claim (that the framework requires the same Poincaré stresses as classical EM, since the field equations and the field momentum are identical), but I have left the claim conditional pending the author's judgement on whether a stronger statement is warranted. No `FINDINGS_for_author_review.md` entry is recommended for this problem; the result is a confirmation of the classical Abraham–Lorentz puzzle, not a new finding.

**Companion notebook:** [`Roadmapping/Mathematica_Notebooks/Electromagnetism/JacksonCh06_P6_4.wl`](../../Mathematica_Notebooks/Electromagnetism/JacksonCh06_P6_4.wl) — runnable independent of the Mathematica MCP.
