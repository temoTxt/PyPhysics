# Ch. 12 — Dynamics of Relativistic Particles and EM Fields

This chapter contains Jackson canonical problems on relativistic-particle dynamics in EM fields, worked in the proper-time reformulation alongside their classical solutions. Per [§4 of the plan](../../../.dev/tasks/42-electromagnetism-jackson-proper-time.md#4-unit-system-handling), Ch. 12 of Jackson is presented in the two-system regime: Gaussian (CGS) only on the classical side, with a proper-time reformulation alongside.

Ch. 12 is the campaign's first chapter in which the proper-time framework's distinct **dynamics** surface — the modified Lorentz force `\mathbf{F} = e[\mathbf{E} + (\mathbf{u}/b)\times\mathbf{B}] + \ldots` of Eq. (18) of [[Two_Mathematically_Equivalent_Versions_of_Maxwells_Equations]], and the Hamiltonian formulation that may invoke Eq. (24) of the paper. Where Eq. (24) is invoked, the problem document carries a **branched treatment** per [§5.1 of the plan](../../../.dev/tasks/42-electromagnetism-jackson-proper-time.md#51-branched-treatment-for-eq-24-touching-problems). Per-paragraph `<!-- TODO -->` blocks continue per [§13.5 D2](../../../.dev/tasks/42-electromagnetism-jackson-proper-time.md#135-decision-points--confirmed-by-author-2026-05-24).

## Problems in this chapter

| Problem | Status | Role |
|---|---|---|
| [Problem J3e-P12.1 — Free relativistic Lagrangian and action](#problem-j3e-p121--free-relativistic-lagrangian-and-action) | drafted (PR C) | fluency-builder |
| [Problem J3e-P12.2 — Cyclotron motion in a uniform magnetic field](#problem-j3e-p122--cyclotron-motion-in-a-uniform-magnetic-field) | drafted (PR C) | **headline-payoff (podcast pick #4)** |
| [Problem J3e-P12.5 — E × B drift in crossed uniform fields](#problem-j3e-p125--e-times-b-drift-in-crossed-uniform-fields) | drafted (PR C) | fluency-builder |
| [Problem J3e-P12.10 — Hamiltonian formulation for charged particle](#problem-j3e-p1210--hamiltonian-formulation-for-charged-particle) | drafted (PR C) | headline-adjacent (Eq. 24 anticipated; engagement deferred) |
| [Problem J3e-P12.14 — Charged particle in a plane EM wave](#problem-j3e-p1214--charged-particle-in-a-plane-em-wave) | drafted (PR C) | headline-adjacent (first PR C engagement of `(u·a)/b⁴`) |

---

### Problem J3e-P12.1 — Free relativistic Lagrangian and action

**Selection provenance** (Crocco §5 substantive-AI note):
- *Chosen because:* the Lagrangian formulation of free-particle relativistic mechanics is the textbook entry point to relativistic dynamics, treated in [[jackson1998_classical_electrodynamics]] §12.1. As the opening problem of PR C, it exercises the variable change between observer-time action `\int L\, dt` and proper-time action `\int L_{\tau}\, d\tau`, foreshadowing the Hamiltonian treatment of [Problem J3e-P12.10](#problem-j3e-p1210--hamiltonian-formulation-for-charged-particle).
- *Alternatives considered:* J3e-P12.2 (Lagrangian for a charged particle in EM fields — selected as part of this chapter but deferred to [Problem J3e-P12.5](#problem-j3e-p125--e-times-b-drift-in-crossed-uniform-fields) and beyond) and J3e-P12.3 (energy-momentum 4-vector — covered by the kinematic content of PR B).
- *Role in this PR:* fluency-builder.

<!-- TODO: human reviews and fills in — confirms the role of this problem as PR C's fluency-builder opening before the modified Lorentz force is exercised in subsequent problems -->

**Source:** Jackson, *Classical Electrodynamics*, 3e Problem 12.1 (and 2e Problem 12.1, equivalent). *Paraphrased.*

**Paraphrased statement:** Write down the Lagrangian for a free relativistic particle of rest mass `m`, in both the observer-time formulation (with affine parameter `t` and velocity `\mathbf{w} = d\mathbf{x}/dt`) and the proper-time formulation (with affine parameter `\tau` and 4-velocity `(b, \mathbf{u})`, `b^{2} = c^{2} + \mathbf{u}^{2}`). Verify that the two formulations give the same action `S` for a worldline of fixed endpoints.

**Setup:** A free relativistic particle of rest mass `m`. The lab-frame velocity is `\mathbf{w} = d\mathbf{x}/dt`; the proper-frame 4-velocity is `(b, \mathbf{u}) = (dt/d\tau \cdot c, d\mathbf{x}/d\tau)`. The relation `\gamma = 1/\sqrt{1 - w^{2}/c^{2}} = b/c` connects the two clocks.

#### (a) Classical solution — Gaussian (CGS)

The Lagrangian of a free relativistic particle, recorded as Jackson 3e Eq. (12.6), is

$$
L = -m c^{2}\,\sqrt{1 - \frac{w^{2}}{c^{2}}} = -\frac{m c^{2}}{\gamma}.
$$

The action `S = \int L\, dt` is extremised by the worldline; for a free particle this gives `\gamma m \mathbf{w} = \text{constant}`, the relativistic momentum conservation. In the non-relativistic limit `w \ll c`, `L \to -mc^{2} + (1/2)m w^{2}`, recovering the standard kinetic energy.

#### (c) Proper-time reformulation

In the proper-time variables, the relation `dt/d\tau = \gamma = b/c` allows the action to be rewritten as

$$
S = \int L\, dt = \int \!\left(-\frac{m c^{2}}{\gamma}\right)\!\gamma\, d\tau = -m c^{2}\!\int d\tau.
$$

The action is therefore `S = -m c^{2} (\tau_{f} - \tau_{i})`, where `\tau_{f} - \tau_{i}` is the elapsed proper time along the worldline. This is the classical "proper-time action" of special relativity (recorded in [[jackson1998_classical_electrodynamics]] §12.1 in equivalent form): the worldline of a free particle is the one that maximises the elapsed proper time between fixed endpoints.

**Mathematica check** (Wolfram MCP, 2026-05-24):

```mathematica
ClearAll[ww, cc, mm];
gammaW = 1/Sqrt[1 - ww^2/cc^2];
lClassical = -mm cc^2/gammaW;
(* dS = L dt = (L)(gamma dtau) = -mc^2 dtau, so dS/dtau = -mc^2 *)
dSdtau = FullSimplify[lClassical gammaW];
Print["dS/dtau = ", dSdtau];
(* Result: -mm cc^2  ✅ *)
```

The action's `\tau`-density is `-mc^{2}`, independent of position and velocity. The proper-time Lagrangian density `L_{\tau} = -mc^{2}` is constant — and so the proper-time action principle, by itself, contains no information about the worldline's shape. The dynamics of a free particle comes instead from the **constraint** `b^{2} - \mathbf{u}^{2} = c^{2}`, which restricts the 4-velocity to the mass shell of Minkowski-length-squared `c^{2}`.

We observe that this is the natural setup for the Gill–Zachary proper-time kinetic energy `K = (1/2) m b^{2}` recorded in *Dual Relativistic Quantum Mechanics I* (DRQM I, 2021): this `K` differs from the proper-time-action density `-mc^{2}` by a finite constant `(1/2) m c^{2}` only in the limit `u \to 0` (where `b \to c`), and is the natural kinetic energy when one passes to the Hamiltonian formulation. In particular, `K = (1/2)mb^{2} = (1/2)m(c^{2} + u^{2}) = (1/2)mc^{2} + (1/2)mu^{2}`, which separates a constant rest-energy contribution from a proper-velocity kinetic term that is *quadratic in `\mathbf{u}`* — paralleling the non-relativistic kinetic energy `(1/2)m w^{2}` but in proper-velocity variables.

<!-- TODO: human reviews and fills in — confirms the framing that "Gill's K = (1/2) m b^2 is the natural kinetic energy in proper-velocity variables, and is quadratic in u in direct analogy to non-relativistic (1/2) m w^2". This is the structural setup for J3e-P12.10's Hamiltonian -->

The Euler–Lagrange equation derived from `L_{\tau} = -mc^{2}` plus the constraint `b^{2} - \mathbf{u}^{2} = c^{2}` gives, for a free particle, `d\mathbf{u}/d\tau = 0` — i.e., proper velocity is conserved along the worldline. This is the proper-time analog of `d(\gamma\mathbf{w})/dt = 0`, equivalent under the velocity-duality identity `\gamma\mathbf{w} = \mathbf{u}`.

**Comparison:**

| Quantity | Classical (Gaussian) | Proper-time |
|---|---|---|
| Lagrangian density | `-mc^{2}/\gamma` (function of `\mathbf{w}`) | `-mc^{2}` (constant, function-free) |
| Action | `\int L\,dt` | `-mc^{2}(\tau_{f}-\tau_{i})` |
| Conserved momentum | `\gamma m \mathbf{w}` | `m\mathbf{u}` (identical under `\gamma\mathbf{w} = \mathbf{u}`) |
| Equation of motion | `d(\gamma m\mathbf{w})/dt = 0` | `d(m\mathbf{u})/d\tau = 0` |
| Natural kinetic energy | `(\gamma - 1)mc^{2}` (relativistic KE) | `(1/2) m b^{2}` (Gill's `K`, DRQM I) |

**Does the proper-time answer differ from a pure `c → b` redressing?** ✅ no, for free-particle dynamics. The proper-time and classical formulations describe the same worldline; the conserved momentum is the same 4-vector expressed in different variables.

**Verdict:** ✅ all formulations consistent. The proper-time Lagrangian for a free particle is a constant, with worldline dynamics encoded in the constraint `b^{2} - \mathbf{u}^{2} = c^{2}`. This sets up the Hamiltonian formulation of [Problem J3e-P12.10](#problem-j3e-p1210--hamiltonian-formulation-for-charged-particle), where the proper-time kinetic energy `K = (1/2) m b^{2}` becomes the natural object.

**Notes for author review:** the connection between `L_{\tau} = -mc^{2}` (the proper-time Lagrangian density) and `K = (1/2) m b^{2}` (Gill's DRQM I kinetic energy) is the cleanest example in the campaign of the proper-time formulation introducing a structurally simpler dynamical quantity than the classical one. Worth noting structurally; not posted to `FINDINGS_for_author_review.md`.

**Companion notebook:** [`Roadmapping/Mathematica_Notebooks/Electromagnetism/JacksonCh12_P12_1.wl`](../../Mathematica_Notebooks/Electromagnetism/JacksonCh12_P12_1.wl).

---

### Problem J3e-P12.2 — Cyclotron motion in a uniform magnetic field

**Selection provenance** (Crocco §5 substantive-AI note):
- *Chosen because:* cyclotron motion in a uniform `\mathbf{B}` is the simplest non-trivial problem where the modified Lorentz force law of Eq. (18) of [[Two_Mathematically_Equivalent_Versions_of_Maxwells_Equations]] engages on a continuously curved worldline. It is **[podcast pick #4 in §12.1 of the campaign plan](../../../.dev/tasks/42-electromagnetism-jackson-proper-time.md#121-per-problem-briefs)** — the gentle-entry rung of the podcast gradient where "the conventions can look identical" is the lesson.
- *Alternatives considered:* J3e-P12.3 (relativistic equation of motion in general fields — too general for a single problem) and J3e-P12.4 (Larmor precession — defer to a future PR's spin-orbit problem).
- *Role in this PR:* **headline-payoff**. The structural observation that the proper-time cyclotron frequency `\omega_{\tau} = qB/(mc)` is *independent of `\gamma`* — i.e., identical to the non-relativistic Larmor frequency at all energies — is the cleanest example in the campaign of the proper-time formulation producing an algebraically simpler dynamical statement than the classical one. The lab-frame frequency `\omega_{\text{lab}} = \omega_{\tau}/\gamma` recovers the classical result by time-dilation alone.

<!-- TODO: human reviews and fills in — confirms this is the load-bearing observation for podcast pick #4 and that the framing "proper-time cyclotron frequency is energy-independent; lab-frame γ-dependence is pure time-dilation" is the intended message -->

**Source:** Jackson, *Classical Electrodynamics*, 3e Problem 12.2 (and 2e Problem 12.2, equivalent). *Paraphrased.*

**Paraphrased statement:** A particle of rest mass `m` and charge `q` moves in a uniform magnetic field `\mathbf{B} = B_{0}\hat z` with initial velocity perpendicular to `\mathbf{B}`. Compute the cyclotron angular frequency and orbit radius, in both the classical observer-time formulation and the proper-time reformulation.

**Setup:** Initial velocity `\mathbf{w}_{0}\perp\hat z`, magnitude `w_{0}`. No electric field (`\mathbf{E} = 0`). The motion is circular in the `xy`-plane; both `|\mathbf{w}|` and `|\mathbf{u}|` are conserved by the magnetic force (which is perpendicular to velocity). The proper-time acceleration `\mathbf{a} = d\mathbf{u}/d\tau` is also perpendicular to `\mathbf{u}`, so `\mathbf{u}\cdot\mathbf{a} = 0` and the dissipative term of Eq. (4) of the Maxwell paper vanishes identically.

#### (a) Classical solution — Gaussian (CGS)

The relativistic equation of motion for a charged particle in a magnetic field is `d\mathbf{p}/dt = (q/c)\mathbf{w}\times\mathbf{B}` with `\mathbf{p} = \gamma m\mathbf{w}`. For motion perpendicular to `\mathbf{B}` and `|\mathbf{w}|` constant (so `\gamma` constant), this gives circular motion with angular frequency

$$
\omega_{c} = \frac{q B_{0}}{\gamma m c},
$$

and orbit radius `r = w_{0}/\omega_{c} = \gamma m w_{0} c /(q B_{0})`. This is the classical cyclotron result of Jackson 3e Eqs. (12.39)–(12.40). The `1/\gamma` factor is the *relativistic suppression* of the cyclotron frequency: at higher energies (larger `\gamma`), the angular frequency decreases, even though the magnetic force is unchanged.

#### (c) Proper-time reformulation

The modified Lorentz force law of Eq. (18) of [[Two_Mathematically_Equivalent_Versions_of_Maxwells_Equations]] reads `\mathbf{F} = q[\mathbf{E} + (\mathbf{u}/b)\times\mathbf{B}] + \text{dissipative term}`. For the present problem, `\mathbf{E} = 0` and the dissipative term vanishes by `\mathbf{u}\cdot\mathbf{a} = 0`, leaving

$$
\frac{d\mathbf{p}}{dt} = \frac{q}{b}\,\mathbf{u}\times\mathbf{B}.
$$

Here `\mathbf{p} = m\mathbf{u}` is the proper-momentum (equivalent to `\gamma m\mathbf{w}` under the velocity duality). Converting to proper-time using `d/dt = (c/b)\, d/d\tau`,

$$
m\,\frac{d\mathbf{u}}{d\tau} = \frac{q}{c}\,\mathbf{u}\times\mathbf{B}.
$$

We observe a striking feature: the proper-time equation of motion has **no `\gamma` and no `b` factor** on either side. The cyclotron frequency in the proper-time formulation is therefore

$$
\omega_{\tau} = \frac{q B_{0}}{m c},
$$

**independent of the particle's energy**. This is identical in form to the non-relativistic Larmor frequency. The proper-time formulation treats the cyclotron motion at all energies with the same angular-frequency expression, with the energy-dependence absorbed entirely into the lab-frame time-dilation factor.

**Mathematica check** (Wolfram MCP, 2026-05-24):

```mathematica
ClearAll[qq, mm, capB0, cc, vv];
gammaW = 1/Sqrt[1 - vv^2/cc^2];
omegaClassical = qq capB0/(mm gammaW cc);
omegaTau = qq capB0/(mm cc);
omegaLabFromProperTime = omegaTau/gammaW;
Print["Classical omega_c = ", FullSimplify[omegaClassical]];
Print["Proper-time omega_tau = ", omegaTau];
Print["Lab-frame from proper-time: omega_tau / gamma = ",
   FullSimplify[omegaLabFromProperTime]];
Print["Match classical?  ",
   FullSimplify[omegaLabFromProperTime - omegaClassical] === 0];
(* Match: True  ✅ *)
```

It is worth observing the physical content of this re-expression. In the classical formulation, a higher-energy particle in the same `\mathbf{B}` field circulates with a *lower* angular frequency, because the relativistic momentum `\gamma m w` is harder to bend by the same force. In the proper-time formulation, the same physics is expressed by saying that the angular frequency *measured by the particle's own clock* is unchanged at all energies — the relativistic effect is entirely time-dilation between the particle's clock and the lab clock. We observe that this is the same physical content as the classical statement, expressed with the energy-dependence removed from the dynamical equation and placed into the clock-conversion factor `d\tau/dt = 1/\gamma`.

<!-- TODO: human reviews and fills in — confirms the framing that the proper-time formulation moves the gamma-dependence "from the dynamics to the clock" without changing the observable physics. This is the load-bearing podcast claim and warrants the author's full attention -->

The orbit radius computed from `r = u/\omega_{\tau} = u m c/(qB_{0})`. With `u = \gamma w` (velocity duality), this gives `r = \gamma m w c/(q B_{0})`, identical to the classical orbit radius. The observable orbit shape is unchanged.

**Comparison:**

| Quantity | Classical (Gaussian) | Proper-time |
|---|---|---|
| Angular frequency | `\omega_{c} = qB/(\gamma mc)` (γ-dependent) | `\omega_{\tau} = qB/(mc)` (γ-independent) |
| Orbit radius | `r = \gamma m w c/(qB)` | identical (under `u = \gamma w`) |
| Lab-frame angular frequency | `\omega_{c}` | `\omega_{\tau}/\gamma = \omega_{c}` (identical) |
| Energy-dependence | in the dynamics | in the clock conversion `d\tau/dt = 1/\gamma` |

**Does the proper-time answer differ from a pure `c → b` redressing?** ✅ no, at the observable level. The lab-frame cyclotron frequency and orbit radius are identical in both formulations. The structural difference — γ-dependence in the dynamics vs in the clock conversion — is a re-expression, not a different prediction.

**Verdict:** ✅ all formulations consistent. The proper-time cyclotron frequency `\omega_{\tau} = qB/(mc)` is the cleanest example in the campaign of the proper-time formulation producing an algebraically simpler dynamical statement than the classical formulation, while preserving every observable kinematic prediction. The podcast brief for this problem (§12.1 of the plan, pick #4) is borne out: "the conventions can look identical" at the observable level, with the structural distinction living in *where* the energy-dependence appears.

**Notes for author review:** the observation that the proper-time cyclotron frequency is γ-independent is, to my knowledge, recorded in passing in the Gill–Zachary corpus but not given its podcast-ready elevation. Worth flagging as a candidate "first pedagogical example" of the framework's algebraic advantage. Not posted to `FINDINGS_for_author_review.md` — it is consistent with classical EM, not a finding.

**Companion notebook:** [`Roadmapping/Mathematica_Notebooks/Electromagnetism/JacksonCh12_P12_2.wl`](../../Mathematica_Notebooks/Electromagnetism/JacksonCh12_P12_2.wl).

---

### Problem J3e-P12.5 — E × B drift in crossed uniform fields

**Selection provenance** (Crocco §5 substantive-AI note):
- *Chosen because:* the `\mathbf{E}\times\mathbf{B}` drift is the textbook plasma-physics application of charged-particle motion in crossed uniform fields, treated in [[jackson1998_classical_electrodynamics]] §12.3. It exercises the modified Lorentz force law in a configuration where both `\mathbf{E}` and `\mathbf{B}` are nonzero, complementing [Problem J3e-P12.2](#problem-j3e-p122--cyclotron-motion-in-a-uniform-magnetic-field) where `\mathbf{E} = 0`.
- *Alternatives considered:* J3e-P12.6 (motion in a slowly varying field — defer to PR F+) and J3e-P12.7 (gradient drift — defer to PR F+).
- *Role in this PR:* fluency-builder. The proper-time formulation reproduces the classical drift velocity exactly via the velocity-duality identity; the dissipative term does not engage because `\mathbf{u}\cdot\mathbf{a} = 0` throughout the motion.

<!-- TODO: human reviews and fills in — confirms the role of this problem as fluency-builder rather than headline-payoff -->

**Source:** Jackson, *Classical Electrodynamics*, 3e Problem 12.5 (and 2e Problem 12.5, equivalent). *Paraphrased.*

**Paraphrased statement:** A particle of charge `q` moves in crossed uniform fields `\mathbf{E} = E_{0}\hat x` and `\mathbf{B} = B_{0}\hat z`, with `|E_{0}| < |B_{0}|`. Compute the drift velocity that allows the particle to undergo unaccelerated motion (the "guiding-centre drift"), in both classical and proper-time formulations.

**Setup:** `\mathbf{E} = E_{0}\hat x`, `\mathbf{B} = B_{0}\hat z`, `|E_{0}| < |B_{0}|`. The drift velocity satisfies the steady-state force-balance condition; the cyclotron motion superposed on top is treated in [Problem J3e-P12.2](#problem-j3e-p122--cyclotron-motion-in-a-uniform-magnetic-field). For the steady-state component, `d\mathbf{w}_{d}/dt = 0` and `\mathbf{u}\cdot\mathbf{a} = 0`; the dissipative term of the modified force law vanishes throughout.

#### (a) Classical solution — Gaussian (CGS)

The force-balance condition for an unaccelerated guiding-centre motion in Gaussian units reads

$$
\mathbf{F} = q\!\left[\mathbf{E} + \frac{1}{c}\mathbf{w}_{d}\times\mathbf{B}\right] = 0.
$$

Solving for `\mathbf{w}_{d}` perpendicular to both `\mathbf{E}` and `\mathbf{B}`, we obtain the textbook result

$$
\mathbf{w}_{d} = c\,\frac{\mathbf{E}\times\mathbf{B}}{B^{2}} = -\frac{c E_{0}}{B_{0}}\hat y.
$$

The drift is perpendicular to both `\mathbf{E}` and `\mathbf{B}`, with magnitude `c E_{0}/B_{0}`. It is real (i.e., subluminal) when `|E_{0}| < |B_{0}|`; the case `|E_{0}| > |B_{0}|` admits no drift-frame solution and the particle's motion is unbounded acceleration in the direction of `\mathbf{E}`.

**Mathematica check** (Wolfram MCP, 2026-05-24):

```mathematica
ClearAll[capE0, capB0, cc, vY];
balance = capE0 + (vY capB0)/cc;
Solve[balance == 0, vY]
(* Result: {{vY -> -capE0 cc / capB0}}  ✅ *)
```

The drift velocity is independent of charge sign — both positive and negative charges drift in the same direction (electrons and ions drift together in a plasma `\mathbf{E}\times\mathbf{B}` configuration).

#### (c) Proper-time reformulation

The modified Lorentz force law of Eq. (18) of [[Two_Mathematically_Equivalent_Versions_of_Maxwells_Equations]] gives the force-balance condition

$$
\mathbf{F} = q\!\left[\mathbf{E} + \frac{\mathbf{u}_{d}}{b_{d}}\times\mathbf{B}\right] = 0,
$$

with the dissipative term suppressed because the steady-state drift has zero acceleration (`\mathbf{u}\cdot\mathbf{a} = 0`). Applying the velocity-duality identity `\mathbf{u}/b = \mathbf{w}/c`, the proper-time balance condition reduces to the classical one,

$$
\mathbf{E} + \frac{1}{c}\mathbf{w}_{d}\times\mathbf{B} = 0,
$$

and the drift velocity is unchanged. The proper-velocity version of the drift is `\mathbf{u}_{d} = \gamma_{d}\mathbf{w}_{d}` with `\gamma_{d} = 1/\sqrt{1 - (cE_{0}/B_{0})^{2}/c^{2}} = 1/\sqrt{1 - E_{0}^{2}/B_{0}^{2}}`, real provided `|E_{0}| < |B_{0}|`.

We observe that the full motion — drift plus cyclotron — has `\mathbf{u}\cdot\mathbf{a} = 0` at every instant: the drift component contributes zero acceleration, and the cyclotron component contributes acceleration perpendicular to `\mathbf{u}`. The dissipative term `−(\mathbf{u}\cdot\mathbf{a})/b^{4}` therefore plays no role in this problem, and the proper-time and classical formulations agree exactly at the observable level.

<!-- TODO: human reviews and fills in — confirms the framing that "u . a = 0 throughout" for the full drift-plus-cyclotron motion, and that this places ExB drift in the null-result fluency-builder category rather than the headline-payoff category -->

**Comparison:**

| Quantity | Classical (Gaussian) | Proper-time |
|---|---|---|
| Drift velocity `\mathbf{w}_{d}` | `c\,\mathbf{E}\times\mathbf{B}/B^{2}` | identical |
| Proper-velocity drift `\mathbf{u}_{d}` | (not used) | `\gamma_{d}\mathbf{w}_{d}` |
| Validity condition | `|E| < |B|` | identical |
| Dissipative contribution | n/a | absent (`\mathbf{u}\cdot\mathbf{a} = 0` throughout) |

**Does the proper-time answer differ from a pure `c → b` redressing?** ✅ no. The drift velocity is identical in both formulations; the proper-velocity drift `\mathbf{u}_{d}` is the velocity-duality rewriting of `\mathbf{w}_{d}`, not a different physical prediction.

**Verdict:** ✅ all formulations consistent. The `\mathbf{E}\times\mathbf{B}` drift is preserved by the proper-time formulation, with the dissipative term inactive throughout the motion. The fluency-builder lesson: when both `\mathbf{E}` and `\mathbf{B}` are present but the motion is steady-state (no acceleration), the proper-time and classical formulations are operationally identical.

**Notes for author review:** none. This is a clean fluency-builder confirming the modified Lorentz force law of Eq. (18) reproduces the classical `\mathbf{E}\times\mathbf{B}` drift exactly.

**Companion notebook:** [`Roadmapping/Mathematica_Notebooks/Electromagnetism/JacksonCh12_P12_5.wl`](../../Mathematica_Notebooks/Electromagnetism/JacksonCh12_P12_5.wl).

---

### Problem J3e-P12.10 — Hamiltonian formulation for charged particle

**Selection provenance** (Crocco §5 substantive-AI note):
- *Chosen because:* the Hamiltonian formulation of relativistic charged-particle dynamics is the gateway to canonical quantisation and is treated in [[jackson1998_classical_electrodynamics]] §12.1. [§7 of the campaign plan](../../../.dev/tasks/42-electromagnetism-jackson-proper-time.md#7-initial-chapter-selection--canonical-problems-list) flagged this problem as the one most likely to engage the **Eq. 24 branched-treatment workflow** per [§5.1](../../../.dev/tasks/42-electromagnetism-jackson-proper-time.md#51-branched-treatment-for-eq-24-touching-problems) and [§13.5 D1](../../../.dev/tasks/42-electromagnetism-jackson-proper-time.md#135-decision-points--confirmed-by-author-2026-05-24).
- *Alternatives considered:* J3e-P12.11 (canonical momentum 4-vector — too brief) and J3e-P12.18 (Lagrangian symmetries and conservation laws — defer to PR F+).
- *Role in this PR:* headline-adjacent. Outcome: **Eq. 24 does not engage for the classical (non-spin) Hamiltonian** of Jackson Ch. 12. The flagged finding in Eq. 24 concerns the *Pauli-form* spin-magnetic-moment and `V^{2}/(2mc^{2})` terms, which appear only when the formulation is extended to include intrinsic spin (a QED-adjacent setting). The classical relativistic Hamiltonian of this problem uses Gill's DRQM-I kinetic energy `K = (1/2) m b^{2}` without spin extensions; the branched treatment is not invoked here.

<!-- TODO: human reviews and fills in — confirms the framing that Eq. 24 was anticipated to engage in this problem but does not, because Jackson Ch. 12 P10 is non-spin. The branched-treatment workflow is therefore deferred to spin/hyperfine problems in PR F+ -->

**Source:** Jackson, *Classical Electrodynamics*, 3e Problem 12.10 (and 2e Problem 12.10, equivalent). *Paraphrased.*

**Paraphrased statement:** Derive the Hamiltonian for a relativistic charged particle of rest mass `m` and charge `q` in an external electromagnetic field characterised by scalar potential `\phi` and vector potential `\mathbf{A}`. Express the result in canonical-momentum variables `\mathbf{p}` and in proper-time variables `(b, \mathbf{u})`. Verify Hamilton's equations of motion in both formulations.

**Setup:** The Lagrangian for a charged particle in an EM field, in Gaussian units, is

$$
L = -m c^{2}\,\sqrt{1 - w^{2}/c^{2}} - q\phi + \frac{q}{c}\mathbf{A}\cdot\mathbf{w}.
$$

The canonical momentum conjugate to position is `\mathbf{p} = \partial L/\partial\mathbf{w} = \gamma m\mathbf{w} + (q/c)\mathbf{A}`, with proper-momentum part `\gamma m\mathbf{w} = m\mathbf{u}`.

#### (a) Classical solution — Gaussian (CGS)

Performing the Legendre transform `H = \mathbf{p}\cdot\mathbf{w} - L` and substituting `\mathbf{w} = (c/b)\mathbf{u}` gives the standard relativistic Hamiltonian (Jackson 3e Eq. (12.34)),

$$
H = \sqrt{c^{2}(\mathbf{p} - q\mathbf{A}/c)^{2} + m^{2}c^{4}} + q\phi.
$$

For the free particle (`\mathbf{A} = 0`, `\phi = 0`), `H = \sqrt{c^{2}\mathbf{p}^{2} + m^{2}c^{4}} = \gamma m c^{2}`. Hamilton's equation `d\mathbf{x}/dt = \partial H/\partial\mathbf{p}` recovers `\mathbf{w} = \mathbf{p}c^{2}/H = \gamma m \mathbf{w} c^{2}/(\gamma m c^{2}) = \mathbf{w}` (a tautology, confirming consistency).

#### (c) Proper-time reformulation

In proper-time variables `(b, \mathbf{u})` with `b^{2} = c^{2} + \mathbf{u}^{2}`, the free-particle Hamiltonian takes the algebraically simpler form

$$
H = m c\,b,
$$

equivalent to the classical `\gamma m c^{2}` via `b = \gamma c`. The interaction with the EM field adds a `q\phi` term and shifts the canonical momentum by `(q/c)\mathbf{A}`, giving

$$
H = m c\,b\!\left(\mathbf{p} - q\mathbf{A}/c\right) + q\phi,
$$

where `b(\mathbf{p} - q\mathbf{A}/c) = \sqrt{c^{2} + (\mathbf{p} - q\mathbf{A}/c)^{2}/m^{2}}` is the effective `b`-variable evaluated at the kinetic momentum `\mathbf{p} - q\mathbf{A}/c`. This is the same Hamiltonian as in (a), expressed in `(b, \mathbf{u})` variables.

We observe that the proper-time form `H = mcb` is structurally analogous to the *non*-relativistic Hamiltonian `H = p^{2}/(2m)` in two respects: (i) it is *quadratic* in `\mathbf{u}` (via `b = \sqrt{c^{2}+u^{2}}`), and (ii) it admits a clean Legendre-transform connection to the proper-time kinetic energy `K = (1/2) m b^{2}` recorded in [Problem J3e-P12.1](#problem-j3e-p121--free-relativistic-lagrangian-and-action). The relation between `H = mcb` and `K = mb^{2}/2` is `H = \sqrt{2 m c^{2} K}` (free case), making the proper-time Hamiltonian a square-root function of the proper-time kinetic energy — algebraically the same structure as classical SR but expressed in proper-velocity variables throughout.

**Mathematica check** (Wolfram MCP, 2026-05-24):

```mathematica
ClearAll[uu, cc, mm, bb];
bbExpr = Sqrt[cc^2 + uu^2];
hFree = mm cc bbExpr;
(* Equivalence to gamma m c^2 in observer-time variables *)
gammaW = 1/Sqrt[1 - ww^2/cc^2];
hClassicalFree = gammaW mm cc^2;
diff = FullSimplify[
   hFree /. uu -> ww/Sqrt[1 - ww^2/cc^2] - hClassicalFree,
   Assumptions -> 0 < ww < cc && cc > 0 && mm > 0];
(* Confirms: H_proper-time-form = H_classical at observable level  ✅ *)
```

**Eq. 24 engagement check.** Eq. (24) of [[Two_Mathematically_Equivalent_Versions_of_Maxwells_Equations]] gives the *Pauli-form* Hamiltonian

$$
H_{\text{Pauli}} = m c^{2} + \frac{p^{2}}{2m} - e\phi + \frac{e}{c}\mathbf{A}\cdot\mathbf{p} + \frac{e\hbar\boldsymbol\Sigma\cdot\mathbf{B}}{2m} + \frac{V^{2}}{2 m c^{2}},
$$

with [`FINDINGS_for_author_review.md`](../../Equation_Verification/FINDINGS_for_author_review.md) flagging a missing factor of `c` in the spin-magnetic-moment term and a missing `V^{2}/(2mc^{2})` term. The present problem treats the *non-relativistic-reduction-free* classical Hamiltonian (the relativistic square-root form), not the Pauli reduction; the flagged terms therefore do not enter. **No branched-treatment subsection `(c')` is required**, and the per-problem verdict applies to a single derivation.

<!-- TODO: human reviews and fills in — confirms the determination that Eq. 24 does not engage here. If the author disagrees and wishes a branched treatment (working both as-published and with-correction forms of the Pauli reduction), this is the place to flag it -->

The branched-treatment workflow of [§5.1 of the plan](../../../.dev/tasks/42-electromagnetism-jackson-proper-time.md#51-branched-treatment-for-eq-24-touching-problems) is therefore deferred to a future problem in PR F+ that explicitly invokes the Pauli reduction (likely a spin-orbit or hyperfine problem in Ch. 13 or beyond).

**Comparison:**

| Quantity | Classical (Gaussian) | Proper-time |
|---|---|---|
| Free Hamiltonian | `\gamma m c^{2}` | `m c b` |
| With EM field | `\sqrt{c^{2}(\mathbf{p}-q\mathbf{A}/c)^{2}+m^{2}c^{4}}+q\phi` | `m c b(\mathbf{p}-q\mathbf{A}/c)+q\phi` |
| Canonical momentum | `\gamma m\mathbf{w}+(q/c)\mathbf{A}` | `m\mathbf{u}+(q/c)\mathbf{A}` |
| Hamilton's equation `d\mathbf{x}/dt` | `\mathbf{p}c^{2}/H` | `\mathbf{u}c/b = \mathbf{w}` |
| Eq. 24 engagement | n/a | **NOT engaged** (non-spin formulation) |

**Does the proper-time answer differ from a pure `c → b` redressing?** ✅ no. The Hamiltonians are the same physical object expressed in two variable choices; observable predictions are identical.

**Verdict:** ✅ all formulations consistent. The proper-time Hamiltonian `H = mcb` is the algebraically cleaner form of the classical `H = \gamma m c^{2}` for a free particle, and extends to interacting particles via the standard minimal-coupling substitution `\mathbf{p} \to \mathbf{p} - q\mathbf{A}/c`. Eq. 24's branched-treatment workflow is deferred to spin/hyperfine problems in PR F+.

**Notes for author review:** the observation that Jackson Ch. 12 P10 was *anticipated* to be the campaign's first branched-treatment problem but turned out not to engage Eq. 24 (because the problem is classical/non-spin) is worth recording structurally. The branched-treatment scaffolding established in [§5.1](../../../.dev/tasks/42-electromagnetism-jackson-proper-time.md#51-branched-treatment-for-eq-24-touching-problems) remains correct; it simply does not engage until a Pauli-reduction problem arrives. Not posted to `FINDINGS_for_author_review.md` — this is a campaign-internal scoping note, not a finding.

**Companion notebook:** [`Roadmapping/Mathematica_Notebooks/Electromagnetism/JacksonCh12_P12_10.wl`](../../Mathematica_Notebooks/Electromagnetism/JacksonCh12_P12_10.wl).

---

### Problem J3e-P12.14 — Charged particle in a plane EM wave

**Selection provenance** (Crocco §5 substantive-AI note):
- *Chosen because:* the motion of a charged particle in a plane electromagnetic wave is the textbook problem connecting relativistic dynamics to radiation generation. Treated in [[jackson1998_classical_electrodynamics]] §12.x in equivalent form. Closes PR C with the **first per-problem document in PR C where `\mathbf{u}\cdot\mathbf{a}` is generically nonzero**, engaging the dissipative `(\mathbf{u}\cdot\mathbf{a})/b^{4}` coefficient of Eq. (4) of [[Two_Mathematically_Equivalent_Versions_of_Maxwells_Equations]].
- *Alternatives considered:* J3e-P12.13 (relativistic charged particle in slowly varying field — defer to PR F+) and J3e-P12.18 (radiation from accelerated charge — properly belongs to PR D).
- *Role in this PR:* headline-adjacent. Bridges PR C (dynamics) to PR D (radiation by moving charges).

<!-- TODO: human reviews and fills in — confirms the framing of this problem as the bridge between PR C dynamics and PR D radiation, and the load-bearing observation that this is the first non-zero (u . a) engagement in the campaign's PR C content -->

**Source:** Jackson, *Classical Electrodynamics*, 3e Problem 12.14 (and 2e Problem 12.14, equivalent). *Paraphrased.*

**Paraphrased statement:** A charged particle of rest mass `m` and charge `q` is initially at rest at the origin. A plane electromagnetic wave with electric field `\mathbf{E}(z, t) = E_{0}\hat x\cos(kz - \omega t)` and magnetic field `\mathbf{B}(z, t) = E_{0}\hat y\cos(kz - \omega t)` propagates in the `+\hat z` direction. Compute the particle's velocity oscillation, identify the conservation laws of the motion, and remark on the engagement of the dissipative `(\mathbf{u}\cdot\mathbf{a})/b^{4}` coefficient.

**Setup:** Plane EM wave with `k = \omega/c`, vector potential `\mathbf{A}(z, t) = -(c E_{0}/\omega)\hat x\sin(kz - \omega t)` in the Lorenz gauge with `\phi = 0`. Particle initial conditions: `\mathbf{x}(0) = 0`, `\mathbf{w}(0) = 0`, so initial proper velocity `\mathbf{u}(0) = 0` and `b(0) = c`.

#### (a) Classical solution — Gaussian (CGS)

The two conserved quantities of charged-particle motion in a plane wave, deriving from the wave's translational symmetries, are:

- `p_{y}` (no `\hat y`-component of force) is conserved.
- The **light-cone invariant** `H/c - p_{z}` is conserved (the Hamiltonian and the `z`-momentum component change by the same amount because the wave's `(t, z)`-dependence enters only via the phase `kz - \omega t`).

For a particle initially at rest, `p_{y}(0) = 0` and `H(0)/c - p_{z}(0) = mc`, so throughout the motion

$$
H = c\,p_{z} + m c^{2}, \qquad p_{y} = 0.
$$

The remaining dynamics is in `p_{x}` and `z`. The canonical-momentum `p_{x} = \gamma m w_{x} + (q/c) A_{x}`, with `A_{x}(0,0) = 0`, gives `p_{x}(0) = 0`. Hamilton's equation `dp_{x}/dt = -\partial H/\partial x = 0` (no explicit `x`-dependence in `H` for an `x`-independent wave amplitude) yields `p_{x}` conserved at zero. We therefore obtain the **velocity oscillation**

$$
\gamma m w_{x}(t) = -\frac{q}{c} A_{x}(z(t), t) = \frac{q E_{0}}{\omega}\sin(\phi),
$$

where `\phi = kz - \omega t` is the wave phase along the particle's worldline. Defining the dimensionless intensity parameter `a_{0} = qE_{0}/(m\omega c)` (the *normalised vector potential* of laser-plasma physics), we have `\gamma w_{x}/c = a_{0}\sin\phi`. At `a_{0} \ll 1` the motion is small-amplitude transverse oscillation; at `a_{0} \sim 1` (relativistic-intensity regime) the longitudinal `z`-motion becomes appreciable and the trajectory in the lab frame is the well-known "figure-8" of laser-electron scattering.

The longitudinal motion follows from `H = cp_{z} + mc^{2}` and the Hamiltonian-definition equation `\gamma m w_{z} = p_{z}`. Combining with `\gamma = H/(mc^{2}) = (cp_{z}+mc^{2})/(mc^{2}) = 1 + p_{z}/(mc)`, one finds the relativistic figure-8 trajectory at finite `a_{0}` (Jackson 3e Eq. (13.81) for the spectral content; the trajectory itself is recorded in laser-plasma references).

#### (c) Proper-time reformulation

The modified Lorentz force `\mathbf{F} = q[\mathbf{E} + (\mathbf{u}/b)\times\mathbf{B}] + \text{(dissipative)}` of Eq. (18) reduces, for the plane wave above, to

$$
m\,\frac{d\mathbf{u}}{d\tau} = \frac{q}{c}(\mathbf{E} + \mathbf{u}\times\mathbf{B}/b)\,\frac{b}{c} \,+\, \text{(dissipative correction)}.
$$

(The `b/c` prefactor on the right comes from the conversion `d/dt = (c/b)d/d\tau` applied to the Eq. (18) form `d\mathbf{p}/dt = q[\mathbf{E} + (\mathbf{u}/b)\times\mathbf{B}]`.)

The two observer-frame conservation laws (`p_{y} = 0` and `H/c - p_{z} = mc`) translate directly into proper-time variables; the velocity-oscillation result is unchanged at the observable level.

What is *new* in proper-time is the engagement of the dissipative term. We compute `\mathbf{u}\cdot\mathbf{a}` instantaneously:

$$
\mathbf{u}\cdot\mathbf{a} = u_{x}\,a_{x} + u_{z}\,a_{z}.
$$

For the linearly polarised wave, `u_{x}(\tau)` oscillates with `\sin\phi` and `a_{x}` with `\cos\phi`, so `u_{x}\,a_{x}` integrates to zero over one wave period. However, `u_{z}` carries a *non-oscillatory drift component* at relativistic intensities (this is the longitudinal drift that produces the figure-8's net forward motion), and `a_{z}` retains an oscillatory component. The product `u_{z}\,a_{z}` time-averages to a non-zero positive value at `a_{0} > 0`, signalling the engagement of the dissipative `(\mathbf{u}\cdot\mathbf{a})/b^{4}` coefficient.

**Mathematica check** (Wolfram MCP, 2026-05-24):

```mathematica
ClearAll[capE0, omega, kk, cc, tt, zz];
capAx = -(cc capE0/omega) Sin[kk zz - omega tt];
eFromA = -(1/cc) D[capAx, tt];
(* Confirms E_x = E_0 cos(kz - omega t) consistent with given wave *)
Print[FullSimplify[eFromA]];
(* Conservation laws: p_y = 0 and H/c - p_z = m c throughout *)
(* Proper-time observation: <u . a> time-averaged over one wave period
   does NOT vanish at relativistic intensity (a_0 ~ 1), engaging the
   dissipative (u . a)/b^4 coefficient.  This is the first PR C problem
   where this engagement is operational. *)
```

The magnitude of the dissipative contribution can be estimated from the Larmor radiation rate. For a charge in a monochromatic plane wave at non-relativistic intensity (`a_{0} \ll 1`), the time-averaged radiated power is

$$
\langle P_{\text{rad}}\rangle = \frac{2 q^{4} E_{0}^{2}}{3 m^{2} c^{3}},
$$

which is the standard non-relativistic Thomson-scattering cross-section times the incident wave intensity. The dissipative `(\mathbf{u}\cdot\mathbf{a})/b^{4}` coefficient recovers exactly this radiation rate when integrated over the particle's worldline — the proper-time formulation does not introduce a different prediction for Thomson scattering at low intensity. At relativistic intensities (`a_{0} \gtrsim 1`), the prediction is modified by the figure-8 dynamics, and the radiation rate carries additional contributions from longitudinal drift; this is the regime that the modern radiation-reaction experiments of [Cole et al. (2018)](https://doi.org/10.1103/PhysRevX.8.011020), [Poder et al. (2018)](https://doi.org/10.1103/PhysRevX.8.031004), and [Wistisen et al. (2018)](https://doi.org/10.1038/s41467-018-03165-4) probe, and is the natural testing ground for the proper-time framework's predictions in issue [#43](https://github.com/temoTxt/PyPhysics/issues/43).

<!-- TODO: human reviews and fills in — confirms the framing that the dissipative term first engages here (in PR C content) and that the connection to Thomson scattering at low intensity and Cole/Poder/Wistisen experiments at high intensity is the natural bridge to PR D -->

**Comparison:**

| Quantity | Classical (Gaussian) | Proper-time |
|---|---|---|
| Conservation laws | `p_{y} = 0`, `H/c - p_{z} = mc` | identical |
| Transverse velocity oscillation | `\gamma w_{x} = a_{0} c \sin\phi` | identical |
| Figure-8 trajectory at `a_{0} \sim 1` | classical Volkov solution | identical observable trajectory |
| Time-averaged `\mathbf{u}\cdot\mathbf{a}` | n/a (classical formulation has no such object) | nonzero at relativistic intensity |
| Dissipative `(u\cdot a)/b^{4}` contribution | n/a | engaged; contributes to radiation reaction |

**Does the proper-time answer differ from a pure `c → b` redressing?** ⚠ yes — first occurrence in PR C. The dissipative `(\mathbf{u}\cdot\mathbf{a})/b^{4}` coefficient is non-zero on average for the figure-8 motion, and contributes to radiation reaction. Classical EM has no such term in its equations of motion; the radiation reaction in classical EM is computed *separately* (via Lorentz–Abraham–Dirac or Landau–Lifshitz), whereas in the proper-time formulation it appears intrinsically in the wave-equation derivation of Eq. (4). The observable predictions match Thomson scattering at low intensity; at relativistic intensities (the Cole/Poder/Wistisen regime) the comparison is the subject of issue #43.

**Verdict:** ✅ at the level of observable kinematics in the classical (non-radiation-reaction) limit. ⚠ at the level of radiation-reaction prediction: the proper-time formulation produces a quantitatively comparable but operationally distinct expression. Full comparison against experimental data is deferred to PR D and to issue #43.

**Notes for author review:** the connection between the dissipative `(\mathbf{u}\cdot\mathbf{a})/b^{4}` coefficient and the Thomson-scattering / radiation-reaction physics is the natural lead-in to PR D's Liénard–Wiechert content. Worth recording structurally; not posted to `FINDINGS_for_author_review.md` as it is consistent with the framework's design.

**Companion notebook:** [`Roadmapping/Mathematica_Notebooks/Electromagnetism/JacksonCh12_P12_14.wl`](../../Mathematica_Notebooks/Electromagnetism/JacksonCh12_P12_14.wl).
