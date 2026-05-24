# Ch. 12 — Dynamics of Relativistic Particles and EM Fields

This chapter contains Jackson canonical problems on relativistic-particle dynamics in EM fields, worked in the proper-time reformulation alongside their classical solutions. Per [§4 of the plan](../../../.dev/tasks/42-electromagnetism-jackson-proper-time.md#4-unit-system-handling), Ch. 12 of Jackson is presented in the two-system regime: Gaussian (CGS) only on the classical side, with a proper-time reformulation alongside.

Ch. 12 is the campaign's first chapter in which the proper-time framework's distinct **dynamics** surface — the modified Lorentz force `\mathbf{F} = e[\mathbf{E} + (\mathbf{u}/b)\times\mathbf{B}] + \ldots` of Eq. (18) of [[Two_Mathematically_Equivalent_Versions_of_Maxwells_Equations]], and the Hamiltonian formulation that may invoke Eq. (24) of the paper. Where Eq. (24) is invoked, the problem document carries a **branched treatment** per [§5.1 of the plan](../../../.dev/tasks/42-electromagnetism-jackson-proper-time.md#51-branched-treatment-for-eq-24-touching-problems). Per-paragraph `<!-- TODO -->` blocks continue per [§13.5 D2](../../../.dev/tasks/42-electromagnetism-jackson-proper-time.md#135-decision-points--confirmed-by-author-2026-05-24).

## Problems in this chapter

| Problem | Status | Role |
|---|---|---|
| [Problem J3e-P12.1 — Free relativistic Lagrangian and action](#problem-j3e-p121--free-relativistic-lagrangian-and-action) | drafted (PR C) | fluency-builder |

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
