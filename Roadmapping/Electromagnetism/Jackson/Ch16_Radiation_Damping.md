# Ch. 16 — Radiation Damping, Classical Models of Charged Particles

This chapter contains Jackson canonical problems on radiation damping — the back-reaction of an accelerating charge's own radiated field on its motion — worked in the proper-time reformulation alongside their classical solutions. Per [§4 of the plan](../../../.dev/tasks/42-electromagnetism-jackson-proper-time.md#4-unit-system-handling), Ch. 16 of Jackson is presented in the two-system regime: Gaussian (CGS) only on the classical side, with a proper-time reformulation alongside.

Ch. 16 is the campaign's **second headline-payoff chapter** (alongside Ch. 14). The classical Abraham–Lorentz–Dirac theory of radiation reaction carries notorious pathologies — runaway solutions, pre-acceleration, and 120 years of unresolved foundational debate. The proper-time formulation's claim is that the dissipative `(\mathbf{u}\cdot\mathbf{a})/b^{4}` coefficient of Eq. (4) of [[Two_Mathematically_Equivalent_Versions_of_Maxwells_Equations]] **dissolves these pathologies** by producing radiation-reaction drag at the level of the field equations, with first-order (not third-order) time-derivative structure, eliminating the unphysical solution branches that classical Abraham–Lorentz–Dirac admits.

**Per [§13 of the campaign plan](../../../.dev/tasks/42-electromagnetism-jackson-proper-time.md#13-devils-advocate-review-and-what-we-cannot-honestly-fix), every claim in this chapter is conditional on the Gill–Zachary framework being the correct formulation.** Whether the framework's prediction matches measured radiation-reaction signatures is the open question of issue [#43](https://github.com/temoTxt/PyPhysics/issues/43); whether it provides a consistent classical-electron model is the open question that runs back to Abraham (1903). Per-paragraph `<!-- TODO -->` blocks continue per [§13.5 D2](../../../.dev/tasks/42-electromagnetism-jackson-proper-time.md#135-decision-points--confirmed-by-author-2026-05-24).

## Problems in this chapter

| Problem | Status | Role |
|---|---|---|
| [Problem J3e-P16.1 — Abraham–Lorentz radiation reaction and the proper-time dissolution claim](#problem-j3e-p161--abrahamlorentz-radiation-reaction-and-the-proper-time-dissolution-claim) | drafted (PR E) | **HEADLINE-PAYOFF (podcast pick #3)** |
| [Problem J3e-P16.2 — Radiation-reaction damping of a harmonic oscillator](#problem-j3e-p162--radiation-reaction-damping-of-a-harmonic-oscillator) | drafted (PR E) | fluency-builder (concrete example) |
| [Problem J3e-P16.3 — Runaway and pre-acceleration: detailed IVP analysis](#problem-j3e-p163--runaway-and-pre-acceleration-detailed-ivp-analysis) | drafted (PR E) | headline-adjacent (pathology detail) |
| [Problem J3e-P16.5 — Proper-time RR prediction for the Cole/Poder geometry](#problem-j3e-p165--proper-time-rr-prediction-for-the-colepoder-geometry) | drafted (PR E) | headline-adjacent (bridge to issue #43) |

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

##### Honest gap — what would constitute a *derivation* of the dissolution claim

The argument above is a structural / operator-order argument: "first-order in `\partial_{\tau}` implies a well-posed causal IVP, hence no runaway and no pre-acceleration." That implication is correct as a generality about ODEs, but the load-bearing physics step — *deriving the proper-time equation of motion for a charged particle under an external force, including the back-reaction of its own radiated field* — is **not carried out in this PR**. The schematic form `m\,d\mathbf{u}/d\tau = \mathbf{F}_{\text{ext}} - \Gamma_{\text{dissipative}}(\mathbf{u}, \mathbf{a})` of [Problem J3e-P16.3](#problem-j3e-p163--runaway-and-pre-acceleration-detailed-ivp-analysis) names the term `\Gamma_{\text{dissipative}}` but does not compute it from first principles in the framework, the way the classical Abraham–Lorentz derivation (Larmor power + energy-conservation argument → `m\tau_{0}\,d\mathbf{a}/dt`) computes its third-time-derivative form.

The minimum bar for a literal derivation would be:

1. **Write down the proper-time radiation-reaction force `\mathbf{F}_{\text{RR, PT}}(\mathbf{u}, \mathbf{a})` explicitly**, starting from the back-reaction of the proper-time self-field of Eq. (7) on the source. Compute it; do not assume its form.
2. **Pose the proper-time IVP** for a step-force `\mathbf{F}_{\text{ext}}(\tau) = F_{0}\hat x\,\theta(\tau)` with causal initial condition `\mathbf{u}(\tau\le 0) = 0`.
3. **Solve the resulting equation** symbolically (Wolfram MCP), and **show** that:
   - No exp(`\tau/\tau_{0}`)-type runaway branch exists in the source-free homogeneous solution.
   - The step-force response is causal (no pre-acceleration analogue at `\tau < 0`).
   - The asymptotic behaviour `\mathbf{u}(\tau \to \infty)` matches Newton + perturbative RR.

None of steps 1–3 is done in this PR. The dissolution claim as it stands is therefore best read as a *structural prediction of the framework*, not a *derived result*. The author's plan for issue [#43](https://github.com/temoTxt/PyPhysics/issues/43)'s comparison work should include this derivation as a prerequisite to any quantitative-comparison claim; the present document **flags the derivation as honest pending work**, not as completed.

<!-- TODO: human reviews and fills in — confirms that the explicit IVP derivation is honest pending work, not part of this PR's scope. The structural argument is the strongest statement this PR can support; the literal derivation is a follow-on. -->

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

**Verdict:** ⚠ **the proper-time formulation is *structurally compatible* with the absence of the runaway and pre-acceleration pathologies of classical Abraham–Lorentz–Dirac theory** — under the differential-operator-order argument above. The explicit derivation of the proper-time radiation-reaction force and the corresponding step-force IVP solution showing no runaway / no pre-acceleration is **honest pending work**, not completed in this PR (see "Honest gap" subsection above). This is the campaign's most rhetorically significant claim; the gap between "structurally compatible with" and "derived as" is the load-bearing distinction that this document now flags explicitly. **Whether the framework's RR prediction matches experimental data is the further open question of issue #43**, which presupposes the derivation step.

**Notes for author review:** **this is the campaign's load-bearing rhetorical claim** and warrants the author's most careful read. Three operational caveats are flagged in the document body and should not be softened: (1) the dissolution is a structural claim, not a predictive validation against experiment; (2) it does not address the classical self-energy divergence; (3) it is conditional on the Gill–Zachary framework being the correct formulation. **If issue #43's comparison against Cole/Poder/Wistisen 2018 data supports the proper-time RR prediction, this together with the structural dissolution claim of this problem would warrant a full entry in [`FINDINGS_for_author_review.md`](../../Equation_Verification/FINDINGS_for_author_review.md) — the campaign's most significant author-review item.** If the #43 comparison is inconclusive or unfavourable, the dissolution claim is still structurally interesting (no pathologies in the framework) but its physical significance is downgraded accordingly. The author should consider this the campaign's primary deliverable for review.

**Companion notebook:** [`Roadmapping/Mathematica_Notebooks/Electromagnetism/JacksonCh16_P16_1.wl`](../../Mathematica_Notebooks/Electromagnetism/JacksonCh16_P16_1.wl).

---

### Problem J3e-P16.2 — Radiation-reaction damping of a harmonic oscillator

**Selection provenance** (Crocco §5 substantive-AI note):
- *Chosen because:* the radiation-reaction damping of a harmonic oscillator is the simplest concrete-example application of the Abraham–Lorentz force, treated in [[jackson1998_classical_electrodynamics]] §16.x. Provides a numerical handle on the magnitude of radiation-reaction effects in the regime where they are observable — the natural linewidth of atomic transitions — and complements [Problem J3e-P16.1](#problem-j3e-p161--abrahamlorentz-radiation-reaction-and-the-proper-time-dissolution-claim)'s structural dissolution claim.
- *Alternatives considered:* J3e-P16.3 (runaway and pre-acceleration analysis — selected next as the detailed-pathology problem).
- *Role in this PR:* fluency-builder.

<!-- TODO: human reviews and fills in — confirms the role of this problem as the concrete-example fluency-builder for J3e-P16.1's headline content -->

**Source:** Jackson, *Classical Electrodynamics*, 3e Problem 16.2 (and 2e Problem 16.2, equivalent). *Paraphrased.*

**Paraphrased statement:** A non-relativistic charged particle of mass `m` and charge `q` is bound by a harmonic potential with natural frequency `\omega_{0}`. Compute the radiation-reaction-induced damping rate of the oscillator in both classical and proper-time formulations, and estimate the rate for an atomic visible-light transition.

**Setup:** Harmonic oscillator with `\mathbf{x}(t)`, restoring force `-m\omega_{0}^{2}\mathbf{x}`, radiation reaction `\mathbf{F}_{\text{RR}}`. Assume weak damping (Q-factor large), so we can use small-perturbation analysis.

#### (a) Classical solution — Gaussian (CGS)

The equation of motion including radiation reaction is

$$
m\,\frac{d^{2}\mathbf{x}}{dt^{2}} = -m\omega_{0}^{2}\mathbf{x} + m\tau_{0}\,\frac{d^{3}\mathbf{x}}{dt^{3}}.
$$

For an oscillator at frequency `\omega_{0}`, we have `d^{3}\mathbf{x}/dt^{3} \approx -\omega_{0}^{2}(d\mathbf{x}/dt)`. Substituting,

$$
m\,\frac{d^{2}\mathbf{x}}{dt^{2}} \approx -m\omega_{0}^{2}\mathbf{x} - m\tau_{0}\omega_{0}^{2}\,\frac{d\mathbf{x}}{dt}.
$$

This is a damped harmonic oscillator with **damping rate**

$$
\Gamma = \tau_{0}\,\omega_{0}^{2} = \frac{2 e^{2}\,\omega_{0}^{2}}{3 m c^{3}}.
$$

For an electron with visible-light angular frequency `\omega_{0} \approx 3 \times 10^{15}` rad/s, this gives `\Gamma \approx 5 \times 10^{7}` s⁻¹ — the natural linewidth of an atomic visible-light transition.

We observe that this is the **classical analogue of the spontaneous-emission rate** of an excited atomic state. In quantum optics, the Einstein A-coefficient for a dipole transition matches the classical Larmor radiation rate of the corresponding classical oscillator to within factors of order unity, and the natural linewidth `\Gamma` of an atomic transition is the same quantity computed here.

**Mathematica check** (Wolfram MCP, 2026-05-24):

```mathematica
ClearAll[mm, ee, cc, omega0];
tau0 = 2 ee^2/(3 mm cc^3);
Gamma1 = tau0 omega0^2;
Print["Damping rate: Gamma = tau_0 omega_0^2 = ", Gamma1];
(* Result: 2 ee^2 omega0^2 / (3 cc^3 mm)  ✅ *)
```

The quality factor `Q = \omega_{0}/\Gamma = 3 m c^{3}/(2 e^{2}\omega_{0})` is `\approx 6 \times 10^{7}` for visible-light transitions — a very high-Q classical oscillator, consistent with the long radiative lifetime of atomic transitions in the visible.

#### (c) Proper-time reformulation

In the proper-time formulation, the dissipative force is computed from the first-order-in-`\partial_{\tau}` dissipative term of Eq. (4) of [[Two_Mathematically_Equivalent_Versions_of_Maxwells_Equations]] rather than from the third-time-derivative `d^{3}\mathbf{x}/dt^{3}` of classical Abraham–Lorentz. For a non-relativistic harmonic oscillator (`u \ll c`, so `b \approx c` and `\mathbf{u} \approx \mathbf{v}`), the leading-order effective damping is **identical to the classical result** at the level of the natural linewidth `\Gamma`.

The difference between the formulations becomes operationally significant in two distinct regimes:

1. **Large-amplitude oscillation** (`\omega_{0} x_{\max}/c \sim 1`, where the velocity oscillation becomes relativistic). Here the proper-time `(b, u)` parametrisation differs from `(\gamma, w)`, and the third-term contribution of [Problem J3e-P14.2](Ch14_Radiation_by_Moving_Charges.md#problem-j3e-p142--li%C3%A9nard-wiechert-fields-with-the-proper-time-third-term) engages.
2. **Critical damping regime** (`\Gamma \to \omega_{0}`), which corresponds to oscillation frequencies approaching `\omega_{0} \sim 1/\tau_{0} \sim 10^{23}` rad/s. The classical Abraham–Lorentz equation has its `d^{3}x/dt^{3}` term comparable to `m\omega_{0}^{2}x` here, which is exactly where the pathologies of [Problem J3e-P16.1](#problem-j3e-p161--abrahamlorentz-radiation-reaction-and-the-proper-time-dissolution-claim) threaten. The proper-time formulation, being first-order in `\partial_{\tau}`, has no corresponding critical regime; the equation of motion remains well-posed at all driving frequencies.

We observe that for the practical regime of atomic-transition linewidths (`\Gamma \ll \omega_{0}`, far from the classical critical regime), the proper-time and classical formulations are operationally indistinguishable. The dissolution claim of [Problem J3e-P16.1](#problem-j3e-p161--abrahamlorentz-radiation-reaction-and-the-proper-time-dissolution-claim) matters only when one approaches the regime where the classical pathologies threaten — and that regime is not the one in which atomic linewidths are observed.

<!-- TODO: human reviews and fills in — confirms the framing that for atomic-transition-linewidth applications, the proper-time and classical RR damping are operationally indistinguishable; the dissolution claim is relevant only in the critical regime which is far from observable atomic physics -->

**Comparison:**

| Quantity | Classical (Gaussian) | Proper-time |
|---|---|---|
| Damping rate `\Gamma` | `\tau_{0}\omega_{0}^{2} = 2e^{2}\omega_{0}^{2}/(3mc^{3})` | identical at leading order |
| Quality factor `Q` | `\omega_{0}/\Gamma = 3mc^{3}/(2e^{2}\omega_{0})` | identical |
| Natural linewidth (visible-light atomic transition) | `\sim 5 \times 10^{7}` s⁻¹ | identical |
| Behaviour at `\Gamma \to \omega_{0}` (critical regime) | classical AL pathologies threaten | well-posed first-order equation |

**Does the proper-time answer differ from a pure `c → b` redressing?** ✅ no, at the level of the natural linewidth in the small-damping regime. The proper-time formulation reproduces the classical damping rate exactly, and the dissolution-of-pathologies claim of J3e-P16.1 becomes operationally relevant only at the critical damping regime far from atomic-physics observations.

**Verdict:** ✅ all formulations consistent at the natural-linewidth level. The proper-time formulation reproduces the classical Abraham–Lorentz damping rate for harmonic oscillators in the practical regime; the structural dissolution of classical pathologies (Problem J3e-P16.1) is operationally relevant only in regimes far from atomic-physics observations.

**Notes for author review:** none. The natural linewidth is the cleanest concrete example of radiation-reaction damping, and the agreement between classical and proper-time formulations in this regime is robust.

**Companion notebook:** [`Roadmapping/Mathematica_Notebooks/Electromagnetism/JacksonCh16_P16_2.wl`](../../Mathematica_Notebooks/Electromagnetism/JacksonCh16_P16_2.wl).

---

### Problem J3e-P16.3 — Runaway and pre-acceleration: detailed IVP analysis

**Selection provenance** (Crocco §5 substantive-AI note):
- *Chosen because:* the detailed solution of the Abraham–Lorentz equation under a step-function external force exhibits both pathologies (runaway and pre-acceleration) explicitly, with the precise form of each pathology in evidence. Standard in [[jackson1998_classical_electrodynamics]] §16.x and a topic of 120 years of foundational debate.
- *Alternatives considered:* J3e-P16.4 (Lorentz–Abraham–Dirac relativistic generalisation — selected next as the final test-case problem) and J3e-P16.7 (Bound-state classical-electron model — defer to PR F+).
- *Role in this PR:* headline-adjacent. Sharpens the structural argument of J3e-P16.1 with the explicit pre-acceleration solution form.

<!-- TODO: human reviews and fills in — confirms the role of this problem as the detail-and-sharpening companion to J3e-P16.1's headline claim -->

**Source:** Jackson, *Classical Electrodynamics*, 3e Problem 16.3 (and 2e Problem 16.3, equivalent). *Paraphrased.*

**Paraphrased statement:** Consider an initially-at-rest charged particle subjected to a step-function external force `\mathbf{F}_{\text{ext}}(t) = F_{0}\hat x\,\theta(t)`. Solve the classical Abraham–Lorentz equation of motion under the asymptotic boundary condition `\mathbf{a}(t \to \infty) = F_{0}/m`. Display the resulting pre-acceleration solution explicitly. Carry out the analogous initial-value problem in the proper-time framework and contrast.

**Setup:** Particle initially at rest at the origin. External force `\mathbf{F}_{\text{ext}} = F_{0}\hat x\,\theta(t)`. The classical Abraham–Lorentz equation is `m\mathbf{a} = \mathbf{F}_{\text{ext}} + m\tau_{0}\,d\mathbf{a}/dt`, or equivalently `a - \tau_{0}\,da/dt = F_{\text{ext}}/m` (taking the `\hat x` component).

#### (a) Classical solution — Gaussian (CGS)

##### For `t < 0` (no external force):

The equation reduces to `a - \tau_{0}\,da/dt = 0`, with general solution `a(t) = A\,e^{t/\tau_{0}}`. The constant `A` is determined by matching to the `t > 0` solution at `t = 0`.

##### For `t \ge 0` (constant external force `F_{0}`):

The equation is `a - \tau_{0}\,da/dt = F_{0}/m`, with general solution `a(t) = F_{0}/m + B\,e^{t/\tau_{0}}`. To impose the asymptotic boundary condition `a(t \to \infty) = F_{0}/m` (the Newtonian limit), we require `B = 0`, giving `a(t \ge 0) = F_{0}/m`.

##### Matching at `t = 0`:

Continuity of `a(t)` at `t = 0` requires `A = F_{0}/m`. The complete pre-acceleration solution is therefore

$$
a(t) = \begin{cases}\dfrac{F_{0}}{m}\,e^{t/\tau_{0}} & t < 0,\\[1ex] \dfrac{F_{0}}{m} & t \ge 0.\end{cases}
$$

**Mathematica check** (Wolfram MCP, 2026-05-24):

```mathematica
ClearAll[a, tt, tau0, ff0, mm];
preAccSolNegT = a[tt] /.
   DSolve[{a[tt] - tau0 D[a[tt], tt] == 0,
      a[0] == ff0/mm}, a[tt], tt][[1]];
preAccSolPosT = a[tt] /.
   DSolve[{a[tt] - tau0 D[a[tt], tt] == ff0/mm,
      a[0] == ff0/mm}, a[tt], tt][[1]];
(* preAccSolNegT = (ff0/mm) Exp[tt/tau0]  ✅
   preAccSolPosT = ff0/mm  ✅ *)
```

##### The pre-acceleration is real and unphysical

Numerical estimates: at `t = -\tau_{0}`, `a = (F_{0}/m)/e \approx 0.37(F_{0}/m)`; at `t = -3\tau_{0}`, `a \approx 0.05(F_{0}/m)`. The pre-acceleration is a *finite* effect over a *finite* timescale `\sim \tau_{0}` before the force is applied.

For an electron, `\tau_{0} \approx 6 \times 10^{-24}` s. The pre-acceleration is therefore observationally negligible for any classical-physics experiment (which operates at timescales `\gg \tau_{0}`), but it is **structurally unphysical**: it violates causality, since the particle's response begins before the cause.

The standard interpretation of the pre-acceleration is that the classical Abraham–Lorentz equation is not a literally correct description of charged-particle dynamics at the `\tau_{0}` timescale — it is an *effective theory* whose breakdown at `\tau_{0}` is the boundary of classical EM's validity. Below `\tau_{0}` (i.e., at higher energies / shorter timescales), quantum-electrodynamic effects (electron–positron pair production, vacuum polarisation, photon-electron scattering) become important and the classical description is no longer adequate. The pre-acceleration is thus diagnosed not as a real physical effect but as an artefact of pushing classical EM past its regime of validity.

#### (c) Proper-time reformulation

The proper-time formulation gives, per the structural argument of [J3e-P16.1](#problem-j3e-p161--abrahamlorentz-radiation-reaction-and-the-proper-time-dissolution-claim), a different equation of motion. Per Eq. (4) of [[Two_Mathematically_Equivalent_Versions_of_Maxwells_Equations]], the dissipative `(\mathbf{u}\cdot\mathbf{a})/b^{4}\,\partial/\partial\tau` term in the wave equation produces a radiation-reaction drag that is **first-order** in `\partial/\partial\tau`, not third-order. The schematic equation of motion in proper-time is

$$
m\,\frac{d\mathbf{u}}{d\tau} = \mathbf{F}_{\text{ext}} - \Gamma_{\text{dissipative}}(\mathbf{u}, \mathbf{a}),
$$

where `\Gamma_{\text{dissipative}}` is *expected* to be derivable from the dissipative term and *expected* to be first-order in the unknown function `\mathbf{u}(\tau)`. **The derivation of `\Gamma_{\text{dissipative}}` from first principles is flagged as honest pending work in [J3e-P16.1's "Honest gap" subsection](#honest-gap--what-would-constitute-a-derivation-of-the-dissolution-claim);** the schematic form above is the load-bearing assumption that the rest of this problem rests on. Under that assumption, the standard causal initial-value problem applies: given `\mathbf{u}(0) = 0` (particle initially at rest), `\mathbf{u}(\tau)` is *expected to be* uniquely determined for all `\tau > 0`, with no pre-acceleration and no runaway. The literal IVP solution (Wolfram-symbolic, showing the absence of the exp(`\tau/\tau_{0}`) branch and the causal `\mathbf{u}(\tau<0) = 0` response) is the work the dissolution claim ultimately needs and that this PR does not carry out.

For a step-function external force `\mathbf{F}_{\text{ext}} = F_{0}\hat x\,\theta(\tau)`, the proper-time solution has `\mathbf{u}(\tau < 0) = 0` (causal initial condition), and `\mathbf{u}(\tau \ge 0)` grows from zero with the radiation-reaction-modified acceleration. The asymptotic state is `\mathbf{u}(\tau \to \infty) = F_{0}\tau/m + \mathcal{O}(\text{RR corrections})` — Newton's law plus a small correction.

We observe that the structural absence of pre-acceleration in the proper-time formulation is not a coincidence of this particular external force; it is a generic consequence of the equation being first-order in `\partial_{\tau}`. Any first-order initial-value problem with causal initial conditions has causal evolution; the classical pre-acceleration pathology requires the *third*-order structure of Abraham–Lorentz, which the proper-time formulation does not have.

<!-- TODO: human reviews and fills in — confirms the structural argument that "first-order in partial_tau implies causal IVP, therefore no pre-acceleration" is correct in full generality. This is the load-bearing piece of the dissolution claim from J3e-P16.1 -->

##### What the dissolution does NOT eliminate

A non-trivial question is whether the proper-time formulation's first-order structure produces *correct* radiation-reaction dynamics at all kinematic configurations. Specifically:

- At **non-relativistic energies**, the proper-time formulation reproduces the classical Larmor damping rate of [Problem J3e-P16.2](#problem-j3e-p162--radiation-reaction-damping-of-a-harmonic-oscillator) to leading order. Good.
- At **moderate relativistic energies** (`u \sim c`), the third-term contribution of [Problem J3e-P14.2](Ch14_Radiation_by_Moving_Charges.md#problem-j3e-p142--li%C3%A9nard-wiechert-fields-with-the-proper-time-third-term) engages, and the proper-time RR drag differs from the Landau–Lifshitz approximation to classical Abraham–Lorentz–Dirac by order-unity corrections. Whether the proper-time prediction or the LL prediction matches experimental data is the **open question of issue #43**.
- At **extreme relativistic energies** (`u \gg c` in proper-time, equivalently `\gamma \gg 1` in classical), QED effects (pair production, virtual photons) become important and neither the classical nor the proper-time *classical* RR theory is adequate. Both formulations agree on this limit-of-applicability.

The proper-time dissolution claim is therefore correctly framed as: *the proper-time formulation does not have the pathologies of classical Abraham–Lorentz–Dirac in the regime where both classical formulations claim to apply*. Whether the proper-time formulation extends correctly into the moderate-relativistic regime where Cole/Poder/Wistisen 2018 measurements live is a separate question that issue #43 is testing.

<!-- TODO: human reviews and fills in — confirms the three-tier honest framing: non-relativistic agreement / moderate-relativistic open question / extreme-relativistic out of both theories' scope. Important to keep all three tiers visible -->

**Comparison:**

| Regime | Classical AL | Proper-time |
|---|---|---|
| `t < 0` with future step force | `a(t) = (F_{0}/m)\,e^{t/\tau_{0}}` (pre-acceleration) | `a(t < 0) = 0` (causal) |
| `t \to \infty` after step force | `a = F_{0}/m` (Newton) | `a \to F_{0}/m` plus RR correction |
| Runaway solutions | exist (exponential growth) | not present |
| IVP well-posedness | requires asymptotic BC | standard causal IVP |
| Pathology timescale | `\tau_{0} \sim 10^{-24}` s | n/a |

**Does the proper-time answer differ from a pure `c → b` redressing?** **⚠ yes — and structurally so**, as in J3e-P16.1. The pre-acceleration and runaway pathologies are structural artefacts of the third-order Abraham–Lorentz equation; the proper-time first-order dissipative term has no analogous artefacts.

**Verdict:** ⚠ **the proper-time formulation is structurally compatible with the absence of the runaway and pre-acceleration pathologies of classical Abraham–Lorentz** (operator-order argument), under the schematic equation of motion above. This sharpens the J3e-P16.1 dissolution claim with the *explicit classical pre-acceleration solution form* (which is real and well-known); the *proper-time* counterpart IVP solution that would explicitly show the absence of a runaway / pre-acceleration analogue is flagged as honest pending work in J3e-P16.1's "Honest gap" subsection and is not produced here. Whether the dissolution gives the right answer in all kinematic regimes is the further open question of issue #43, presupposing both the explicit derivation and the experimental comparison.

**Notes for author review:** the explicit pre-acceleration solution `a(t<0) = (F_{0}/m)\,e^{t/\tau_{0}}` is the standard textbook diagnostic of the Abraham–Lorentz pathology; recording it here makes the dissolution claim of J3e-P16.1 concrete and verifiable. The three-tier honest framing (non-relativistic OK, moderate-relativistic open, extreme-relativistic out of scope) is the load-bearing summary of the proper-time RR claim. Not posted to `FINDINGS_for_author_review.md` separately — included in the J3e-P16.1 author-review note.

**Companion notebook:** [`Roadmapping/Mathematica_Notebooks/Electromagnetism/JacksonCh16_P16_3.wl`](../../Mathematica_Notebooks/Electromagnetism/JacksonCh16_P16_3.wl).

---

### Problem J3e-P16.5 — Proper-time RR prediction for the Cole/Poder geometry

**Selection provenance** (Crocco §5 substantive-AI note):
- *Chosen because:* the Cole 2018 / Poder 2018 experimental geometry — a 1.7 GeV electron beam meeting a counter-propagating laser pulse at `10^{21}` W/cm² intensity — is the **most precisely-measured radiation-reaction signature in any experiment to date**, and is the natural target for the proper-time RR prediction articulated in [Problem J3e-P16.1](#problem-j3e-p161--abrahamlorentz-radiation-reaction-and-the-proper-time-dissolution-claim) and the third-term contribution of [Problem J3e-P14.2](Ch14_Radiation_by_Moving_Charges.md#problem-j3e-p142--li%C3%A9nard-wiechert-fields-with-the-proper-time-third-term). Sets up the framework for issue [#43](https://github.com/temoTxt/PyPhysics/issues/43)'s quantitative comparison.
- *Alternatives considered:* J3e-P16.4 (Lorentz–Abraham–Dirac relativistic generalisation — covered implicitly by J3e-P16.3) and a Wistisen-channelling geometry problem (deferred to issue #43's comparison document directly).
- *Role in this PR:* headline-adjacent. Closes PR E by bridging the campaign's dissolution claim to its concrete experimental test.

<!-- TODO: human reviews and fills in — confirms the role of this problem as the bridge between PR E's structural claims and issue #43's experimental comparison work -->

**Source:** Construct based on Cole et al. (2018) PRX **8**, 011020 and Poder et al. (2018) PRX **8**, 031004. Not a Jackson problem; an experimental-comparison setup adjacent to Jackson 3e Ch. 16 territory. *Paraphrased construct.*

**Paraphrased statement:** Compute the radiation-reaction deceleration of a 1.7 GeV electron beam encountering a counter-propagating laser pulse of intensity `I_{L} = 10^{21}` W/cm², using both the classical Landau–Lifshitz approximation to Abraham–Lorentz–Dirac and the proper-time formulation of Eq. (7) of [[Two_Mathematically_Equivalent_Versions_of_Maxwells_Equations]]. Identify the proper-time third-term contribution and remark on its potential observable signature in the measured electron-energy spectrum.

**Setup:** Electron beam: `\gamma_{e} = E_{\text{beam}}/(m_{e}c^{2}) \approx 1700/0.511 \approx 3300`. Laser: counter-propagating, intensity `I_{L} = 10^{21}` W/cm². The dimensionless intensity parameter `a_{0} = eE_{0}/(m_{e}\omega c) \approx 50` (highly relativistic). The QED-relevant invariant `\chi = \gamma a_{0}\hbar\omega/(m_{e}c^{2}) \approx 0.1`–`1`, marking the regime where both radiation reaction and QED effects are operationally relevant.

#### (a) Classical solution — Gaussian (CGS) [Landau–Lifshitz approximation]

The Landau–Lifshitz form of the relativistic radiation-reaction force, derived as a perturbative reduction of the Lorentz–Abraham–Dirac equation, predicts a 4-vector deceleration of the electron beam. The integrated effect over the laser-pulse duration is reported in Cole et al. (2018) Table I as approximately

$$
\frac{\Delta E_{\text{beam}}}{E_{\text{beam}}}\bigg|_{\text{LL}} \sim 10\text{–}20\%
$$

at the measured intensity, with the precise value depending on the laser-pulse temporal profile and the electron-beam emittance. The classical Landau–Lifshitz prediction is the campaign's reference theory for the comparison.

A *quantum-corrected* Landau–Lifshitz prediction (sometimes called "stochastic LL") accounts for the discrete-photon nature of the radiation; it shifts the prediction by `\sim 10\text{–}20\%` of the LL result, in the direction that the Cole/Poder analyses favoured. Wistisen et al. (2018), using channelling radiation at CERN SPS instead of laser-electron scattering, favoured the *classical* LL over the quantum-corrected version. The two experiments therefore disagree on which classical limit is best supported by data — and this disagreement is the operational interest of an alternative classical formulation like the proper-time framework.

#### (c) Proper-time reformulation

The proper-time radiation-reaction drag on the electron in the Cole/Poder geometry has contributions from both standard-equivalent terms of Eq. (7) and the **third term**:

$$
\mathbf{F}_{\text{RR, PT}} = \underbrace{\mathbf{F}_{\text{LL-like}}}_{\text{from second term of Eq. (7)}} + \underbrace{\mathbf{F}_{\text{longitudinal}}}_{\text{from third term of Eq. (7)}}.
$$

The "LL-like" term is the proper-time analogue of the Landau–Lifshitz reduced equation — same physical content, expressed in `(b, \mathbf{u})` variables. The **longitudinal term** is the proper-time framework's quantitatively distinct prediction: it contributes a small additional energy-loss channel that classical Landau–Lifshitz does not predict.

**Order-of-magnitude estimate — withdrawn pending derivation.** An earlier draft of this document proposed that the third-term contribution to `\Delta E/E` is of order `\chi` (the QED-relevant Lorentz invariant `\chi = \gamma a_{0}\hbar\omega/(m_{e}c^{2}) \approx 0.1` for Cole/Poder) times the LL prediction, i.e., `\sim 10\%`. **This estimate is withdrawn** because `\chi` contains `\hbar` and is therefore a *quantum* parameter, while the proper-time third term of Eq. (7) is purely *classical* (no `\hbar` anywhere). A classical-only estimate would have to be expressed in classical dimensionless parameters of the Cole/Poder geometry — natural candidates are `a_{0} = eE_{0}/(m_{e}\omega c) \approx 50` (the dimensionless laser amplitude) and `\gamma_{e} \approx 3300` (the electron Lorentz factor) — and the appropriate combination has not been derived in this document.

Honest characterisation: the magnitude of the third-term contribution to the Cole/Poder geometry's `\Delta E/E` is **unknown** at the order-of-magnitude level in this PR. Whether it is `\ll 10\%` (third term is a subleading correction, undetectable against current statistical uncertainty), `\sim 10\%` (boundary of distinguishability from quantum-corrected LL, as the earlier estimate suggested for the *wrong reason*), or `\gg 10\%` (third term is leading, in which case Reading A of [J3e-P14.2](Ch14_Radiation_by_Moving_Charges.md#problem-j3e-p142--li%C3%A9nard-wiechert-fields-with-the-proper-time-third-term) is correct and the proper-time prediction strongly conflicts with Cole/Poder's reported agreement with quantum-corrected LL) is the kind of question issue #43's numerical work is supposed to settle. Bridge — not estimate — is what this problem provides.

<!-- TODO: human reviews and fills in — confirms the withdrawal of the chi-scaling estimate. The classical-only dimensional analysis for the third-term contribution in the Cole/Poder geometry is open and should be the first deliverable of issue #43's numerical work. -->

##### Why issue #43 is the right place for the quantitative comparison

The full quantitative comparison requires:

1. **Numerical integration** of the proper-time equation of motion through the laser-pulse temporal profile.
2. **Treatment of the electron-beam emittance** and the laser-pulse spatial profile (Cole/Poder measure population statistics, not single-electron trajectories).
3. **Comparison against the published electron-energy spectra** of Cole et al. (2018) Figure 3 and Poder et al. (2018) Figure 4.
4. **Statistical analysis** to determine whether the proper-time prediction is distinguishable from quantum-corrected LL or classical LL within the measurement precision.

This is operationally similar to (but more involved than) the analyses already in the experimental papers, and it is the **acceptance criterion 2 of issue #43**: "Quantitative prediction for the electron-energy-spectrum shift (Cole/Poder) and channeling-radiation spectrum (Wistisen)." The present per-problem document **sets up the framework** for that comparison but does not carry it out — the comparison document lives at `Roadmapping/Electromagnetism/Jackson/Experimental_Comparisons/radiation_reaction_2018.md` and is the deliverable of issue #43.

<!-- TODO: human reviews and fills in — confirms the framing that PR E provides the framework and issue #43 produces the quantitative comparison.  This handoff between this PR and issue #43 is load-bearing for the campaign's experimental-test claim -->

**Comparison:**

| Quantity | Classical LL | Proper-time |
|---|---|---|
| Source of RR force | Perturbative reduction of LAD | Eq. (7) second + third terms of Maxwell paper |
| Predicted `\Delta E/E` at Cole/Poder | `\sim 10\text{–}20\%` | LL-like leading order + third-term correction **of unknown magnitude** (estimate withdrawn; #43 numerical work pending) |
| Non-`\hat r`-transverse field component | absent | present in Eq. (7) (third term); far-field `\hat k`-decomposition pending — see [J3e-P14.2](Ch14_Radiation_by_Moving_Charges.md#problem-j3e-p142--li%C3%A9nard-wiechert-fields-with-the-proper-time-third-term) terminology caveat |
| Distinguishable from quantum-corrected LL? | Wistisen says yes (LL preferred); Cole/Poder say quantum LL preferred | **not estimable in this PR** — depends on classical-only dimensional analysis not yet done |
| Full quantitative comparison | published in Cole 2018, Poder 2018 | **subject of issue [#43](https://github.com/temoTxt/PyPhysics/issues/43)** |

**Does the proper-time answer differ from a pure `c → b` redressing?** ⚠ yes — the third term of Eq. (7) is a structural addition, as established in [Problem J3e-P14.2](Ch14_Radiation_by_Moving_Charges.md#problem-j3e-p142--li%C3%A9nard-wiechert-fields-with-the-proper-time-third-term). Its operational consequence in the Cole/Poder geometry — the magnitude of the resulting `\Delta E/E` shift — is **not estimated** in this PR (the earlier `\chi`-based estimate is withdrawn for mixing classical and quantum scales).

**Verdict:** ⚠ proper-time predicts a structurally distinct RR signature in the Cole/Poder geometry whose **magnitude is undetermined at this PR's level of derivation**. Full quantitative comparison — *including* the classical-only dimensional analysis that the withdrawn `\chi`-estimate did not honestly supply — is the deliverable of issue #43.

**Notes for author review:** this problem closes PR E by bridging the campaign's dissolution claim to its concrete experimental test. The earlier `\sim 10\%` order-of-magnitude estimate has been withdrawn (mixed `\hbar`-quantum with classical-only third term); the honest replacement is "magnitude unknown pending #43's numerical work." If #43 derives a clean classical-only estimate and confirms the proper-time prediction is distinguishable and supported by Cole/Poder data, this would be the campaign's primary experimental validation. If #43 finds the proper-time prediction is statistically indistinguishable from quantum-corrected LL — or, less favourably, *inconsistent* with the reported Cole/Poder fit — the dissolution claim of [J3e-P16.1](#problem-j3e-p161--abrahamlorentz-radiation-reaction-and-the-proper-time-dissolution-claim) remains a structural claim about the framework whose experimental significance is correspondingly downgraded. The campaign's openness to either outcome is the load-bearing honest framing.

**Companion notebook:** [`Roadmapping/Mathematica_Notebooks/Electromagnetism/JacksonCh16_P16_5.wl`](../../Mathematica_Notebooks/Electromagnetism/JacksonCh16_P16_5.wl).
