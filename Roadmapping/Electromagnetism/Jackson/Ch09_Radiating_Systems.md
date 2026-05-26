# Ch. 9 — Radiating Systems, Multipole Fields and Radiation

This chapter contains Jackson canonical problems on multipole radiation, worked in the proper-time reformulation alongside their classical solutions. Per [§4 of the plan](../../../.dev/tasks/42-electromagnetism-jackson-proper-time.md#4-unit-system-handling), Ch. 9 of Jackson is presented in the three-system regime.

Ch. 9 is part of **PR I (backfill)** per [§7 of the plan](../../../.dev/tasks/42-electromagnetism-jackson-proper-time.md#7-initial-chapter-selection--canonical-problems-list). Unlike PR F (Ch. 7 vacuum waves) and PR G (Ch. 4 statics) where `\mathbf{u} = 0` throughout, **Ch. 9 radiation problems involve sources with finite `\mathbf{u}\cdot\mathbf{a}`** — the same regime where the proper-time third term of Eq. (7) engages, as established in [Problem J3e-P14.2](Ch14_Radiation_by_Moving_Charges.md#problem-j3e-p142--li%C3%A9nard-wiechert-fields-with-the-proper-time-third-term). For oscillating sources at non-relativistic amplitudes, the third-term contribution is `O((u/c)^{2})` below the classical result and below typical observational floors, but its structural presence is real. Per-problem documents flag this engagement explicitly.

## Problems in this chapter

| Problem | Status | Role |
|---|---|---|
| [Problem J3e-P9.1 — Electric dipole radiation](#problem-j3e-p91--electric-dipole-radiation) | drafted (PR I) | fluency-builder (dissipative-term `O((u/c)^{2})` correction) |
| [Problem J3e-P9.2 — Magnetic dipole radiation](#problem-j3e-p92--magnetic-dipole-radiation) | drafted (PR I) | fluency-builder |
| [Problem J3e-P9.3 — Electric quadrupole radiation](#problem-j3e-p93--electric-quadrupole-radiation) | drafted (PR I) | fluency-builder |
| [Problem J3e-P9.4 — Center-fed linear antenna](#problem-j3e-p94--center-fed-linear-antenna) | drafted (PR I) | fluency-builder |

---

### Problem J3e-P9.1 — Electric dipole radiation

**Selection provenance:** the simplest non-trivial radiation source; foundational for all of Ch. 9. *Pragmatic AI use.*

**Paraphrased statement:** An electric dipole oscillates with time-dependence `\mathbf{p}(t) = \mathbf{p}_{0}\cos(\omega t)`. Compute the radiation field, the angular distribution, and the total radiated power.

#### Classical solution

In the far-zone (radiation zone), the electric dipole field is

$$
\mathbf{E}_{\text{rad}}(\mathbf{r}, t) = \frac{1}{c^{2} r}\,[\hat r \times (\ddot{\mathbf{p}}\times\hat r)]\bigg|_{\text{ret}},
$$

with the standard "transverse" structure (perpendicular to `\hat r`). The angular distribution is `\sin^{2}\theta` (the classic doughnut pattern). The total radiated power is the Larmor-style expression

$$
P = \frac{2 |\ddot{\mathbf{p}}|^{2}}{3 c^{3}} = \frac{2 p_{0}^{2}\omega^{4}}{3 c^{3}}.
$$

#### Proper-time reformulation

The oscillating dipole has `\mathbf{u}\cdot\mathbf{a} \neq 0` instantaneously (the charge velocity and acceleration are 90° out of phase, so `\mathbf{u}\cdot\mathbf{a}` oscillates with twice the dipole frequency). Per the third term of [Problem J3e-P14.2](Ch14_Radiation_by_Moving_Charges.md#problem-j3e-p142--li%C3%A9nard-wiechert-fields-with-the-proper-time-third-term), this engages the dissipative coefficient `(\mathbf{u}\cdot\mathbf{a})/b^{4}` and introduces a small longitudinal-radiation-component correction.

For oscillation amplitude `x_{\max}` and frequency `\omega`, the proper-time velocity is `u \sim \omega x_{\max}`, and the fractional correction to the radiated power is `\sim (u/c)^{2} = (\omega x_{\max}/c)^{2}`. For atomic-scale oscillations (`x_{\max} \sim 10^{-10}` m, optical `\omega \sim 10^{15}` Hz), `\omega x_{\max}/c \sim 10^{-5}/3 \times 10^{8} \sim 10^{-12}` — far below any observational floor.

**Verdict:** ✅ identical to classical at the level of the leading dipole-radiation formula. ⚠ structural `O((u/c)^{2})` correction from the third term, below observational floor for atomic-scale dipole oscillations.

**Companion notebook:** [`Roadmapping/Mathematica_Notebooks/Electromagnetism/JacksonCh09_P9_1.wl`](../../Mathematica_Notebooks/Electromagnetism/JacksonCh09_P9_1.wl).

---

### Problem J3e-P9.2 — Magnetic dipole radiation

**Selection provenance:** the next-order radiation source after electric dipole; suppressed by `(\omega/c)^{2}` relative to electric dipole. *Pragmatic AI use.*

**Paraphrased statement:** A magnetic dipole oscillates with time-dependence `\mathbf{m}(t) = \mathbf{m}_{0}\cos(\omega t)`. Compute the radiation pattern and the total radiated power.

#### Classical solution

The magnetic-dipole radiation power is suppressed by `(v/c)^{2}` relative to electric-dipole, since the charge source motion is solenoidal rather than translational:

$$
P = \frac{2 |\ddot{\mathbf{m}}|^{2}}{3 c^{5}} = \frac{2 m_{0}^{2}\omega^{4}}{3 c^{5}}.
$$

The angular distribution is `\sin^{2}\theta` (same as electric dipole), but the electric field is rotated 90° in the plane perpendicular to `\hat r`.

#### Proper-time reformulation

Magnetic dipole oscillation involves rotating currents rather than oscillating charge separation. The mean `\mathbf{u}\cdot\mathbf{a}` over a cycle is non-zero (similar to the electric-dipole case), with the same `(u/c)^{2}` suppression of the third-term contribution.

**Verdict:** ✅ identical at leading order, with the same structural `O((u/c)^{2})` correction.

**Companion notebook:** [`Roadmapping/Mathematica_Notebooks/Electromagnetism/JacksonCh09_P9_2.wl`](../../Mathematica_Notebooks/Electromagnetism/JacksonCh09_P9_2.wl).

---

### Problem J3e-P9.3 — Electric quadrupole radiation

**Selection provenance:** the third radiation channel after electric and magnetic dipoles; further suppressed by `(\omega/c)^{2}`. *Pragmatic AI use.*

**Paraphrased statement:** Compute the radiation from a time-varying electric quadrupole moment `Q_{ij}(t)`.

#### Classical solution

The total radiated power is

$$
P = \frac{1}{180 c^{5}}\sum_{ij}|\dddot{Q}_{ij}|^{2},
$$

with the third-time-derivative of the quadrupole moment characteristic of the radiation. The angular distribution has more complex `\theta`-dependence than the dipole `\sin^{2}\theta`.

#### Proper-time reformulation

Same structural `O((u/c)^{2})` third-term correction as electric dipole; below observational floor for atomic-scale quadrupole oscillations.

**Verdict:** ✅ identical at leading order.

**Companion notebook:** [`Roadmapping/Mathematica_Notebooks/Electromagnetism/JacksonCh09_P9_3.wl`](../../Mathematica_Notebooks/Electromagnetism/JacksonCh09_P9_3.wl).

---

### Problem J3e-P9.4 — Center-fed linear antenna

**Selection provenance:** the canonical practical-radiation-system problem; bridges Ch. 9 theory to antenna engineering. *Pragmatic AI use.*

**Paraphrased statement:** A center-fed linear antenna of half-length `L/2` carries a sinusoidal current distribution. Compute the radiation pattern and the radiation resistance.

#### Classical solution

For a thin linear antenna with current `I(z) = I_{0}\sin(k|L/2 - |z||)`, the far-zone radiation pattern is

$$
\frac{dP}{d\Omega} = \frac{I_{0}^{2}}{2\pi c}\left|\frac{\cos(kL\cos\theta/2) - \cos(kL/2)}{\sin\theta}\right|^{2}.
$$

The half-wave dipole (`kL = \pi`) has radiation resistance `R_{\text{rad}} \approx 73`Ω.

#### Proper-time reformulation

Antenna current is a non-relativistic drift of conduction electrons; `\mathbf{u} \sim 10^{-3}` m/s `\ll c`. The third-term correction is `\sim (10^{-3}/3\times 10^{8})^{2} \sim 10^{-22}` — astronomically below observational floor. The classical radiation resistance and pattern are unchanged.

**Verdict:** ✅ identical.

**Companion notebook:** [`Roadmapping/Mathematica_Notebooks/Electromagnetism/JacksonCh09_P9_4.wl`](../../Mathematica_Notebooks/Electromagnetism/JacksonCh09_P9_4.wl).
