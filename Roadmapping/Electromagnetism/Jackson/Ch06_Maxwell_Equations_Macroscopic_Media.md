# Ch. 6 — Maxwell Equations, Macroscopic Electromagnetism, Conservation Laws

This chapter contains Jackson canonical problems on Maxwell's equations and their conservation laws, worked in the proper-time reformulation alongside their classical CGS and SI solutions. Per [§4 of the plan](../../../.dev/tasks/42-electromagnetism-jackson-proper-time.md#4-unit-system-handling), Ch. 6 of Jackson is presented in the three-system regime.

Ch. 6 is the campaign's first chapter in which the proper-time machinery engages on every problem. The configurations are no longer static or purely-magnetostatic; time-varying fields, displacement currents, and (in [Problem J3e-P6.20](#problem-j3e-p620--radiation-pressure-on-a-perfect-conductor)) the dissipative `−(\mathbf{u}\cdot\mathbf{a})/b^{4}` coefficient of Eq. (4) of [[Two_Mathematically_Equivalent_Versions_of_Maxwells_Equations]] all enter the discussion. Per [§13.5 D2 of the plan](../../../.dev/tasks/42-electromagnetism-jackson-proper-time.md#135-decision-points--confirmed-by-author-2026-05-24), per-paragraph `<!-- TODO: human reviews and fills in -->` blocks accompany each substantive interpretive claim throughout this chapter.

## Problems in this chapter

| Problem | Status | Role |
|---|---|---|
| [Problem J3e-P6.1 — Maxwell with magnetic monopoles + electric–magnetic duality](#problem-j3e-p61--maxwell-with-magnetic-monopoles-and-electricmagnetic-duality) | drafted (PR A) | headline-adjacent |
| [Problem J3e-P6.4 — EM momentum of a uniformly-moving point charge](#problem-j3e-p64--em-momentum-of-a-uniformly-moving-point-charge) | planned (PR A) | headline-payoff (podcast pick #1) |
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
