# Ch. 12 — Dynamics of Relativistic Particles and EM Fields

This chapter contains Jackson canonical problems on relativistic-particle dynamics in EM fields, worked in the proper-time reformulation alongside their classical solutions. Per [§4 of the plan](../../../.dev/tasks/42-electromagnetism-jackson-proper-time.md#4-unit-system-handling), Ch. 12 of Jackson is presented in the two-system regime: Gaussian (CGS) only on the classical side, with a proper-time reformulation alongside.

Ch. 12 is the campaign's first chapter in which the proper-time framework's distinct **dynamics** surface — the modified Lorentz force `\mathbf{F} = e[\mathbf{E} + (\mathbf{u}/b)\times\mathbf{B}] + \ldots` of Eq. (18) of [[Two_Mathematically_Equivalent_Versions_of_Maxwells_Equations]], and the Hamiltonian formulation that may invoke Eq. (24) of the paper. Where Eq. (24) is invoked, the problem document carries a **branched treatment** per [§5.1 of the plan](../../../.dev/tasks/42-electromagnetism-jackson-proper-time.md#51-branched-treatment-for-eq-24-touching-problems). Per-paragraph `<!-- TODO -->` blocks continue per [§13.5 D2](../../../.dev/tasks/42-electromagnetism-jackson-proper-time.md#135-decision-points--confirmed-by-author-2026-05-24).

## Problems in this chapter

| Problem | Status | Role |
|---|---|---|
| [Problem J3e-P12.1 — Free relativistic Lagrangian and action](#problem-j3e-p121--free-relativistic-lagrangian-and-action) | drafted (PR C) | fluency-builder |
| [Problem J3e-P12.2 — Cyclotron motion in a uniform magnetic field](#problem-j3e-p122--cyclotron-motion-in-a-uniform-magnetic-field) | drafted (PR C) | **headline-payoff (podcast pick #4)** |

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
