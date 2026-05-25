# Ch. 13 — Collisions, Energy Loss, Scattering of Charged Particles

This chapter contains Jackson canonical problems on charged-particle energy loss in matter, worked in the proper-time reformulation alongside their classical solutions. Part of **PR J (backfill)** per [§7 of the plan](../../../.dev/tasks/42-electromagnetism-jackson-proper-time.md#7-initial-chapter-selection--canonical-problems-list). Compact format.

## Problems in this chapter

| Problem | Status | Role |
|---|---|---|
| [Problem J3e-P13.1 — Bohr classical stopping-power formula](#problem-j3e-p131--bohr-classical-stopping-power-formula) | drafted (PR J) | fluency-builder |
| [Problem J3e-P13.5 — Cherenkov radiation](#problem-j3e-p135--cherenkov-radiation) | drafted (PR J) | fluency-builder |
| [Problem J3e-P13.6 — Bethe stopping power for relativistic particle](#problem-j3e-p136--bethe-stopping-power-for-relativistic-particle) | drafted (PR J) | fluency-builder |

---

### Problem J3e-P13.1 — Bohr classical stopping-power formula

**Selection provenance:** the original Bohr (1913) calculation of energy loss; foundational for the Bethe formula. *Pragmatic AI.*

**Paraphrased statement:** A fast charged particle traverses a medium of atomic-electron density `n`. Compute the classical stopping power `-dE/dx` due to small-angle Coulomb scattering off the atomic electrons.

#### Classical solution

Bohr's formula (Gaussian):

$$
-\frac{dE}{dx} = \frac{4\pi e^{4} n}{m v^{2}}\ln\!\left(\frac{b_{\max}}{b_{\min}}\right),
$$

with the impact-parameter cutoffs determined by atomic-binding and kinematic constraints. At non-relativistic speeds, this is the textbook stopping power.

#### Proper-time reformulation

The incident charged particle has `u/c \sim v/c`; at MeV-scale energies `v/c \sim 0.9` and `u/c` is finite. Third-term correction is `O((u/c)^{2})` but suppressed further by the small-angle Coulomb-scattering kinematics. Negligible at MeV scales; potentially measurable at GeV scales (the regime of issue #43).

**Verdict:** ✅ identical at leading Bohr order; small-angle scattering kinematics suppress the third-term contribution.

**Companion notebook:** [`Roadmapping/Mathematica_Notebooks/Electromagnetism/JacksonCh13_P13_1.wl`](../../Mathematica_Notebooks/Electromagnetism/JacksonCh13_P13_1.wl).

---

### Problem J3e-P13.5 — Cherenkov radiation

**Selection provenance:** the radiation produced when a charged particle's velocity exceeds the local phase velocity of light in a medium. *Pragmatic AI.*

**Paraphrased statement:** A charged particle moves with velocity `v > c/n` through a medium of refractive index `n`. Compute the Cherenkov radiation cone angle and the radiated energy spectrum.

#### Classical solution

Cherenkov cone angle: `\cos\theta_{c} = c/(n v) = 1/(n\beta)`. Frank–Tamm spectrum:

$$
\frac{d^{2}E}{dx\, d\omega} = \frac{e^{2}}{c^{2}}\,\omega\!\left[1 - \frac{1}{n^{2}\beta^{2}}\right].
$$

#### Proper-time reformulation

The particle moves with `u/c = v/(c\,\gamma^{-1}) > c/n` in the medium frame, so `u/b = v/c > 1/n` and the Cherenkov condition `b\cos\theta_{c} = c/n` reproduces the classical formula. Same Frank–Tamm spectrum at leading order; same `O((u/c)^{2})` third-term correction structure as elsewhere in PR I/J.

**Verdict:** ✅ identical at leading Frank–Tamm order.

**Companion notebook:** [`Roadmapping/Mathematica_Notebooks/Electromagnetism/JacksonCh13_P13_5.wl`](../../Mathematica_Notebooks/Electromagnetism/JacksonCh13_P13_5.wl).

---

### Problem J3e-P13.6 — Bethe stopping power for relativistic particle

**Selection provenance:** the relativistic generalisation of Bohr's stopping formula; the standard medical-physics stopping-power formula. *Pragmatic AI.*

**Paraphrased statement:** Generalise the Bohr stopping-power formula to a relativistic charged particle. Recover the Bethe formula.

#### Classical solution

Bethe's formula (Gaussian, classical limit):

$$
-\frac{dE}{dx} = \frac{4\pi e^{4} n Z}{m_{e} c^{2}\beta^{2}}\!\left[\ln\!\left(\frac{2 m_{e} c^{2}\beta^{2}\gamma^{2}}{I}\right) - \beta^{2}\right],
$$

with `I` the mean ionisation energy of the medium. At ultra-relativistic energies (`\gamma \gg 1`), the formula scales as `\ln(2 m_{e} c^{2}\gamma^{2}/I)`. This is the foundation of charged-particle dosimetry in medical physics.

#### Proper-time reformulation

The incident particle has `u/c` finite at MeV-GeV energies; third-term correction is `O((u/c)^{2})` but with logarithmic suppression from the impact-parameter integration. The proper-time formulation reproduces Bethe at leading order with the same small-angle-scattering simplifications.

**Verdict:** ✅ identical at leading Bethe order.

**Notes for author review:** This problem is the natural connection point between the campaign's RR work (issues #43, #48) and stopping-power physics in medical-linac calibration. The third-term contribution to bremsstrahlung (issue #48) is in the same energy regime as Bethe stopping; a comprehensive comparison would relate the two.

**Companion notebook:** [`Roadmapping/Mathematica_Notebooks/Electromagnetism/JacksonCh13_P13_6.wl`](../../Mathematica_Notebooks/Electromagnetism/JacksonCh13_P13_6.wl).
