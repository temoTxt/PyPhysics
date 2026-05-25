# Ch. 16 — Radiation Damping, Classical Models of Charged Particles

This chapter contains Jackson canonical problems on radiation damping — the back-reaction of an accelerating charge's own radiated field on its motion — worked in the proper-time reformulation alongside their classical solutions. Per [§4 of the plan](../../../.dev/tasks/42-electromagnetism-jackson-proper-time.md#4-unit-system-handling), Ch. 16 of Jackson is presented in the two-system regime: Gaussian (CGS) only on the classical side, with a proper-time reformulation alongside.

Ch. 16 is the campaign's **second headline-payoff chapter** (alongside Ch. 14). The classical Abraham–Lorentz–Dirac theory of radiation reaction carries notorious pathologies — runaway solutions, pre-acceleration, and 120 years of unresolved foundational debate. The proper-time formulation's claim is that the dissipative `(\mathbf{u}\cdot\mathbf{a})/b^{4}` coefficient of Eq. (4) of [[Two_Mathematically_Equivalent_Versions_of_Maxwells_Equations]] **dissolves these pathologies** by producing radiation-reaction drag at the level of the field equations, with first-order (not third-order) time-derivative structure, eliminating the unphysical solution branches that classical Abraham–Lorentz–Dirac admits.

**Per [§13 of the campaign plan](../../../.dev/tasks/42-electromagnetism-jackson-proper-time.md#13-devils-advocate-review-and-what-we-cannot-honestly-fix), every claim in this chapter is conditional on the Gill–Zachary framework being the correct formulation.** Whether the framework's prediction matches measured radiation-reaction signatures is the open question of issue [#43](https://github.com/temoTxt/PyPhysics/issues/43); whether it provides a consistent classical-electron model is the open question that runs back to Abraham (1903). Per-paragraph `<!-- TODO -->` blocks continue per [§13.5 D2](../../../.dev/tasks/42-electromagnetism-jackson-proper-time.md#135-decision-points--confirmed-by-author-2026-05-24).

## Problems in this chapter

| Problem | Status | Role |
|---|---|---|
| [Problem J3e-P16.1 — Abraham–Lorentz radiation reaction and the proper-time dissolution claim](#problem-j3e-p161--abrahamlorentz-radiation-reaction-and-the-proper-time-dissolution-claim) | drafted (PR E) | **HEADLINE-PAYOFF (podcast pick #3)** |

---

### Problem J3e-P16.1 — Abraham–Lorentz radiation reaction and the proper-time dissolution claim

**Selection provenance** (Crocco §5 substantive-AI note):
- *Chosen because:* the Abraham–Lorentz radiation-reaction force is the foundational result of classical radiation-damping theory, treated in [[jackson1998_classical_electrodynamics]] §16.x and recorded historically across Abraham (1903) → Lorentz (1904) → Dirac (1938) → Wheeler–Feynman (1945) → Landau–Lifshitz (1962) → Rohrlich / Spohn (2000s). It is **podcast pick #3** in [§12.1 of the campaign plan](../../../.dev/tasks/42-electromagnetism-jackson-proper-time.md#121-per-problem-briefs) and the campaign's most rhetorically loaded problem: the proper-time framework's claim to *dissolve* the classical pathologies by producing radiation-reaction drag intrinsically in the field equations rather than as an auxiliary force.
- *Alternatives considered:* J3e-P16.2 (radiation-reaction damping of an oscillator — selected next as the concrete-example fluency-builder) and J3e-P16.5 (Lorentz model of an electron — defer to PR F+ as classical-electron-model territory).
- *Role in this PR:* **HEADLINE-PAYOFF**. This is the campaign's most consequential single per-problem document.

<!-- TODO: human reviews and fills in — confirms the role of this problem as the campaign's most rhetorically loaded document, and the framing of the proper-time "dissolution claim" as conditional on the framework being the correct formulation -->

**Source:** Jackson, *Classical Electrodynamics*, 3e Problem 16.1 (and 2e Problem 16.1, equivalent). *Paraphrased.*

**Paraphrased statement:** Derive the Abraham–Lorentz radiation-reaction force `\mathbf{F}_{\text{RR}}` on a non-relativistic point charge `q` of mass `m`. Identify the classical pathologies of the resulting equation of motion (runaway solutions and pre-acceleration). Carry out the analogous derivation in the proper-time framework of [[Two_Mathematically_Equivalent_Versions_of_Maxwells_Equations]], and remark on the framework's claim to dissolve the classical pathologies.

**Setup:** Non-relativistic point charge of mass `m` and charge `q`. External force `\mathbf{F}_{\text{ext}}(t)`. The radiated power is computed from the Larmor formula (treated in [Problem J3e-P14.3](Ch14_Radiation_by_Moving_Charges.md#problem-j3e-p143--non-relativistic-larmor-radiation-formula)) and the back-reaction on the particle is derived by an energy-conservation argument.

#### (a) Classical solution — Gaussian (CGS)

The Larmor power radiated by the accelerating charge is `P = (2e^{2}a^{2})/(3c^{3})`. By energy conservation, the particle must lose this energy to the radiation field. Writing `dE/dt = -P + \mathbf{F}_{\text{RR}}\cdot\mathbf{v}` and assuming `\mathbf{F}_{\text{RR}}` is proportional to `da/dt` (the only choice that makes `\mathbf{F}_{\text{RR}}\cdot\mathbf{v} = -P` time-averaged for a bounded motion), one obtains the **Abraham–Lorentz formula**

$$
\mathbf{F}_{\text{RR}} = \frac{2 e^{2}}{3 c^{3}}\,\frac{d\mathbf{a}}{dt} = m\,\tau_{0}\,\frac{d\mathbf{a}}{dt},
$$

with the **radiation-reaction time scale**

$$
\tau_{0} = \frac{2 e^{2}}{3 m c^{3}}.
$$

For an electron, `\tau_{0} \approx 6.3 \times 10^{-24}` s — comparable to the Compton time `\hbar/(m_{e}c^{2})` and many orders of magnitude shorter than any classical EM time scale. The equation of motion is

$$
m\,\mathbf{a} = \mathbf{F}_{\text{ext}} + m\,\tau_{0}\,\frac{d\mathbf{a}}{dt}.
$$

##### The runaway pathology

When `\mathbf{F}_{\text{ext}} = 0`, the equation reduces to `d\mathbf{a}/dt = \mathbf{a}/\tau_{0}`, with the solution

$$
\mathbf{a}(t) = \mathbf{a}(0)\,\exp(t/\tau_{0}).
$$

Any initial acceleration grows exponentially with timescale `\tau_{0}`. For an electron with even an infinitesimal initial acceleration, this predicts unbounded acceleration on a `\sim 10^{-23}` s timescale — **unphysical**. This is the "runaway solution" pathology, identified by Abraham (1903) and never satisfactorily resolved in classical Abraham–Lorentz–Dirac theory.

**Mathematica check** (Wolfram MCP, 2026-05-24):

```mathematica
ClearAll[ee, mm, cc, tt, tau0, alpha];
tau0 = 2 ee^2/(3 mm cc^3);
runawaySol = DSolve[{D[a[tt], tt] == a[tt]/tau0, a[0] == alpha}, a[tt], tt];
Print[runawaySol];
(* Result: a(t) = alpha Exp[t/tau_0]  ✅ exponential runaway *)
```

##### The pre-acceleration pathology

To eliminate the runaway, one imposes the *asymptotic boundary condition* `\mathbf{a}(t \to \infty) \to \mathbf{F}_{\text{ext}}(\infty)/m`. Solving the resulting integral equation gives, for a step force `\mathbf{F}_{\text{ext}} = F_{0}\,\theta(t)\,\hat x` turned on at `t = 0`,

$$
\mathbf{a}(t) = \begin{cases}\dfrac{F_{0}}{m}\exp(t/\tau_{0})\,\hat x & t < 0, \\ \dfrac{F_{0}}{m}\,\hat x & t \ge 0.\end{cases}
$$

The particle *accelerates before the force is applied*, over the `\tau_{0} \sim 10^{-24}` s timescale. This is the "pre-acceleration" pathology — equally unphysical, equally unresolved in classical Abraham–Lorentz–Dirac theory.

The relativistic generalisation (the Lorentz–Abraham–Dirac equation, Dirac 1938) inherits both pathologies in 4-vector form. Wheeler–Feynman absorber theory (1945) attempts a resolution by invoking advanced potentials; Landau–Lifshitz (1962) gives a perturbative "reduced" equation that avoids runaway but has its own foundational issues; Spohn (2000s) and others have given rigorous existence proofs under restrictive conditions but the pathologies remain in the unrestricted theory. After 120 years, the consensus position in classical Abraham–Lorentz–Dirac theory is that the pathologies are **artefacts of the point-particle limit** and the equations are not literally physical at the timescale `\tau_{0}`.

#### (c) Proper-time reformulation

The Gill–Zachary proper-time formulation produces radiation-reaction drag intrinsically in the field equations rather than as an auxiliary force. Per Eq. (4) of [[Two_Mathematically_Equivalent_Versions_of_Maxwells_Equations]], the proper-time wave equation for the electric field of an accelerating source carries a **dissipative term**

$$
-\frac{\mathbf{u}\cdot\mathbf{a}}{b^{4}}\,\frac{\partial \mathbf{E}}{\partial \tau},
$$

which acts to oppose changes in the field consistent with energy loss to radiation. The radiation-reaction force on the particle is then computed by integrating the back-reaction of the self-field on the source — but the key structural distinction from classical Abraham–Lorentz–Dirac is that this dissipative term is **first-order in `\partial/\partial\tau`**, not third-order.

##### The proper-time dissolution claim

The claim of the Gill–Zachary framework, articulated around Eq. (4) of the Maxwell paper:

> *"The new term in equation (4) is dissipative, acts to oppose the acceleration, is zero when `\mathbf{a} = 0` and arises instantaneously with the action of forces on the particle. … Thus, in this approach, there is no need to assume advanced potentials, self-interaction, mass renormalization along with the Lorentz–Dirac equation in order to account for it (radiation reaction)…"*

Because the dissipative term is first-order in `\partial_{\tau}`, the equation of motion in the proper-time formulation does **not** carry the third-time-derivative `da/dt` term that produces the classical Abraham–Lorentz runaway. The integral equation for `\mathbf{a}(t)` is structurally different; the exponential-runaway solution branch does not exist. By the same argument, the pre-acceleration pathology — which arises from imposing causal boundary conditions on a third-order equation — does not arise in the proper-time formulation because the underlying equation is first-order in `\partial_{\tau}` and the standard causal initial-value problem is well-posed.

<!-- TODO: human reviews and fills in — confirms the framing that "first-order in partial_tau, therefore no third-order pathologies" is the load-bearing structural claim of the dissolution. This is the campaign's most rhetorically loaded paragraph and warrants the author's full read -->

##### What this dissolution claim does NOT establish

To be honest about the scope of the dissolution claim:

1. **It does not prove that the proper-time radiation-reaction force matches experimental data.** Whether the proper-time prediction agrees with the Cole/Poder/Wistisen 2018 measurements is the subject of issue #43, which is open. The dissolution is a *structural* claim about the absence of pathologies in the proper-time formulation, not a *predictive* claim that the proper-time RR force is numerically correct at all energies.

2. **It does not address the classical-electron self-energy problem.** The classical electron's electromagnetic self-energy diverges in the point-particle limit, independent of the radiation-reaction question. The proper-time formulation has not been shown to dissolve this distinct pathology; the framework concerns the radiation back-reaction on the particle's motion, not its self-energy.

3. **It is conditional on the framework being the correct formulation.** Per [§13 of the campaign plan](../../../.dev/tasks/42-electromagnetism-jackson-proper-time.md#13-devils-advocate-review-and-what-we-cannot-honestly-fix), every claim in this campaign is "*if the Gill–Zachary framework is the correct formulation, then X*". The dissolution claim is meaningful only conditional on that premise. Whether the premise is correct is the open question of the framework itself.

<!-- TODO: human reviews and fills in — confirms the three honest caveats: not proven experimentally, not addressing self-energy, conditional on framework correctness.  These caveats are load-bearing for the campaign's framing per §13, and must not be softened -->

##### Comparison of structural form

| Equation | Classical Abraham–Lorentz | Proper-time (Gill–Zachary) |
|---|---|---|
| Form | `m\mathbf{a} = \mathbf{F}_{\text{ext}} + m\tau_{0}\,d\mathbf{a}/dt` | `m\,d\mathbf{u}/d\tau = \mathbf{F}_{\text{ext}} + \text{(first-order in }\partial_{\tau}\text{ dissipation)}` |
| Time-derivative order | third (via `da/dt = d^{3}x/dt^{3}`) | first (via the `\partial_{\tau}` dissipation) |
| Runaway solutions | exist, with timescale `\tau_{0}` | **not present** by structural absence |
| Pre-acceleration | required to eliminate runaway | **not required** (initial-value problem well-posed) |
| Mass renormalization | required (point-particle limit) | not addressed (separate question) |
| Advanced potentials | invoked by Wheeler–Feynman | not invoked |
| Experimental validation | partial (Cole/Poder/Wistisen 2018; debated) | conditional on framework correctness; #43 is testing |

**Does the proper-time answer differ from a pure `c → b` redressing?** **⚠ YES — and structurally, in a load-bearing way.** The proper-time formulation does not just rescale the Abraham–Lorentz coefficient by `c \to b`; it produces a structurally different equation of motion with no third-time-derivative term. The dissolution of the classical pathologies is a consequence of this structural difference, not of a coefficient rescaling.

**Verdict:** ⚠ **the proper-time formulation predicts the absence of the runaway and pre-acceleration pathologies of classical Abraham–Lorentz–Dirac theory.** This is consistent with the framework as published and is the campaign's most rhetorically significant result. **Whether this constitutes a correct physical description of radiation reaction is the open question of issue #43.**

**Notes for author review:** **this is the campaign's load-bearing rhetorical claim** and warrants the author's most careful read. Three operational caveats are flagged in the document body and should not be softened: (1) the dissolution is a structural claim, not a predictive validation against experiment; (2) it does not address the classical self-energy divergence; (3) it is conditional on the Gill–Zachary framework being the correct formulation. **If issue #43's comparison against Cole/Poder/Wistisen 2018 data supports the proper-time RR prediction, this together with the structural dissolution claim of this problem would warrant a full entry in [`FINDINGS_for_author_review.md`](../../Equation_Verification/FINDINGS_for_author_review.md) — the campaign's most significant author-review item.** If the #43 comparison is inconclusive or unfavourable, the dissolution claim is still structurally interesting (no pathologies in the framework) but its physical significance is downgraded accordingly. The author should consider this the campaign's primary deliverable for review.

**Companion notebook:** [`Roadmapping/Mathematica_Notebooks/Electromagnetism/JacksonCh16_P16_1.wl`](../../Mathematica_Notebooks/Electromagnetism/JacksonCh16_P16_1.wl).
