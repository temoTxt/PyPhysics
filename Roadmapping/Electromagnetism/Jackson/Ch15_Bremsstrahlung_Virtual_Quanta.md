# Ch. 15 — Bremsstrahlung, Method of Virtual Quanta

This chapter contains Jackson canonical problems on bremsstrahlung radiation theory and the Weizsäcker–Williams method of virtual quanta, worked in the proper-time reformulation alongside their classical solutions. Part of **PR J (backfill)** per [§7 of the plan](../../../.dev/tasks/42-electromagnetism-jackson-proper-time.md#7-initial-chapter-selection--canonical-problems-list). Compact format.

**Note:** Ch. 15 is closely related to issue [#48](https://github.com/temoTxt/PyPhysics/issues/48) (precision MeV-bremsstrahlung experimental comparison). The problems here treat bremsstrahlung from the *theoretical*-formulation side; issue #48 addresses the experimental-comparison side.

## Problems in this chapter

| Problem | Status | Role |
|---|---|---|
| [Problem J3e-P15.1 — Bremsstrahlung in classical Born approximation](#problem-j3e-p151--bremsstrahlung-in-classical-born-approximation) | drafted (PR J) | fluency-builder; **issue #48 reference** |
| [Problem J3e-P15.5 — Method of virtual quanta](#problem-j3e-p155--method-of-virtual-quanta) | drafted (PR J) | fluency-builder |
| [Problem J3e-P15.9 — Soft-photon limit and infrared divergence](#problem-j3e-p159--soft-photon-limit-and-infrared-divergence) | drafted (PR J) | fluency-builder |

---

### Problem J3e-P15.1 — Bremsstrahlung in classical Born approximation

**Selection provenance:** the classical bremsstrahlung cross-section; foundational for the QED Born-approximation calculation. Directly relevant to issue #48's experimental-comparison work. *Pragmatic AI.*

**Paraphrased statement:** Compute the classical-mechanics differential bremsstrahlung cross-section for an electron decelerating in the Coulomb field of a nucleus, in the non-relativistic Born approximation.

#### Classical solution

The non-relativistic Bethe–Heitler-style cross-section (Gaussian) takes the form

$$
\frac{d\sigma}{d\hbar\omega} = \frac{16}{3}\,\alpha\,r_{e}^{2}\,Z^{2}\,\frac{1}{\hbar\omega}\,\ln\!\left(\frac{E_{i} + p_{i}c}{E_{f} + p_{f}c}\right) \approx \text{(constant)}\cdot\frac{1}{\omega}
$$

with the famous `1/\omega` low-frequency behaviour (the *infrared catastrophe* signature, regularised in QED by soft-photon resummation).

#### Proper-time reformulation

For the bremsstrahlung geometry (electron decelerating in nuclear Coulomb field), the electron has `u\cdot a` finite (linear deceleration along the trajectory). The third term of Eq. (7) engages with magnitude `(u\cdot a)/b^{4}` per [Problem J3e-P14.6](Ch14_Radiation_by_Moving_Charges.md#problem-j3e-p146--bremsstrahlung-from-a-linearly-decelerating-charge). At MeV electron energies the correction is order unity relative to classical Bethe–Heitler; at non-relativistic energies it is `(u/c)^{2}` suppressed.

**Verdict:** ⚠ third term engages at MeV energies; experimental significance is the subject of **issue [#48](https://github.com/temoTxt/PyPhysics/issues/48)**.

**Companion notebook:** [`Roadmapping/Mathematica_Notebooks/Electromagnetism/JacksonCh15_P15_1.wl`](../../Mathematica_Notebooks/Electromagnetism/JacksonCh15_P15_1.wl).

---

### Problem J3e-P15.5 — Method of virtual quanta

**Selection provenance:** the Weizsäcker–Williams method — treating the Coulomb field of a fast-moving charged particle as a flux of "virtual photons". *Pragmatic AI.*

**Paraphrased statement:** Express the Coulomb field of a relativistic charged particle as an equivalent photon flux in a frame where the particle is at rest.

#### Classical solution

A relativistic particle's Coulomb field is Lorentz-contracted into a transverse-EM-pulse spectrum. The number of equivalent photons per unit frequency per unit impact parameter is

$$
\frac{dN}{d\omega} \approx \frac{2 Z^{2}\alpha}{\pi}\,\frac{1}{\omega}\,\ln\!\left(\frac{b_{\max}}{b_{\min}}\right),
$$

logarithmic dependence on impact parameter cut-offs. This is the equivalent-photon method that underlies many high-energy QED calculations.

#### Proper-time reformulation

Same Lorentz-contraction argument applies; the proper-time formulation's `(b/c)` factor on the current density and `(1/b)` factor on the field gradients combine to leave the equivalent-photon flux unchanged at the level of total number per `\omega`.

**Verdict:** ✅ identical at leading order.

**Companion notebook:** [`Roadmapping/Mathematica_Notebooks/Electromagnetism/JacksonCh15_P15_5.wl`](../../Mathematica_Notebooks/Electromagnetism/JacksonCh15_P15_5.wl).

---

### Problem J3e-P15.9 — Soft-photon limit and infrared divergence

**Selection provenance:** the soft-photon emission limit; classical analogue of the QED infrared divergence. *Pragmatic AI.*

**Paraphrased statement:** Show that the classical bremsstrahlung cross-section has a `1/\omega` low-frequency divergence (the classical-mechanics infrared catastrophe). Identify the physical resolution.

#### Classical solution

In the soft-photon limit `\hbar\omega \to 0`, the bremsstrahlung cross-section scales as `d\sigma/d\omega \to (\text{const})/\omega`. Integrated over all frequencies, this is logarithmically divergent at the soft end.

The physical resolution is that *real* (measurable) cross-sections always include a soft-photon-emission resolution `\omega_{\min}`, which regularises the divergence. The number of total emitted photons is infinite but each carries vanishing energy; the total energy is finite. In QED, the Bloch–Nordsieck mechanism resums infrared divergences and produces finite measurable cross-sections.

#### Proper-time reformulation

Same low-frequency divergence structure (the third term contributes at fixed `\omega`, not specifically at the soft limit). Same physical resolution via soft-photon-resolution regulator.

**Verdict:** ✅ identical at the level of the infrared structure.

**Companion notebook:** [`Roadmapping/Mathematica_Notebooks/Electromagnetism/JacksonCh15_P15_9.wl`](../../Mathematica_Notebooks/Electromagnetism/JacksonCh15_P15_9.wl).
