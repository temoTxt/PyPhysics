# Ch. 8 — Waveguides, Resonant Cavities, Optical Fibers

This chapter contains Jackson canonical problems on EM-wave propagation in waveguides and resonant cavities, worked in the proper-time reformulation alongside their classical solutions. Per [§4 of the plan](../../../.dev/tasks/42-electromagnetism-jackson-proper-time.md#4-unit-system-handling), Ch. 8 of Jackson is presented in the three-system regime.

Ch. 8 is part of **PR H (backfill)** per [§7 of the plan](../../../.dev/tasks/42-electromagnetism-jackson-proper-time.md#7-initial-chapter-selection--canonical-problems-list). All Ch. 8 configurations are source-free wave-equation problems with boundary conditions; the proper-time formulation reduces identically since `\mathbf{u} = 0` everywhere. Compact per-problem format per PR F precedent.

## Problems in this chapter

| Problem | Status | Role |
|---|---|---|
| [Problem J3e-P8.1 — TE and TM modes in rectangular waveguide](#problem-j3e-p81--te-and-tm-modes-in-rectangular-waveguide) | drafted (PR H) | fluency-builder |
| [Problem J3e-P8.2 — Cylindrical waveguide modes](#problem-j3e-p82--cylindrical-waveguide-modes) | drafted (PR H) | fluency-builder |
| [Problem J3e-P8.5 — Resonant cavity modes](#problem-j3e-p85--resonant-cavity-modes) | drafted (PR H) | fluency-builder |
| [Problem J3e-P8.7 — Power loss and Q factor in waveguide](#problem-j3e-p87--power-loss-and-q-factor-in-waveguide) | drafted (PR H) | fluency-builder |

---

### Problem J3e-P8.1 — TE and TM modes in rectangular waveguide

**Selection provenance:** the canonical waveguide problem; Jackson 3e §8.4. *Pragmatic AI use.*

**Paraphrased statement:** A rectangular waveguide of cross-section `a \times b` has perfectly conducting walls. Find the TE and TM mode solutions and the dispersion relation.

#### Classical solution

For TE modes (`E_{z} = 0`, `B_{z}` characteristic), the dispersion relation is

$$
k_{z}^{2} = \frac{\omega^{2}}{c^{2}} - \left(\frac{m\pi}{a}\right)^{2} - \left(\frac{n\pi}{b}\right)^{2}.
$$

Modes propagate when `\omega > \omega_{mn,\text{cutoff}} = c\sqrt{(m\pi/a)^{2} + (n\pi/b)^{2}}`. TM modes have analogous structure with `B_{z} = 0` and `E_{z}` characteristic. Lowest mode is typically TE_{10}.

#### Proper-time reformulation

Wave equation inside waveguide, no source: `\mathbf{u} = 0`, `b = c`. The proper-time formulation reduces identically; dispersion and cutoff frequencies are unchanged.

**Verdict:** ✅ identical. Waveguide modes are unaffected by the proper-time formulation.

**Companion notebook:** [`Roadmapping/Mathematica_Notebooks/Electromagnetism/JacksonCh08_P8_1.wl`](../../Mathematica_Notebooks/Electromagnetism/JacksonCh08_P8_1.wl).

---

### Problem J3e-P8.2 — Cylindrical waveguide modes

**Selection provenance:** waveguide with circular cross-section; Bessel-function radial structure. *Pragmatic AI use.*

**Paraphrased statement:** A cylindrical waveguide of radius `R` has perfectly conducting walls. Find the TE and TM mode solutions in terms of Bessel functions.

#### Classical solution

TM modes have `E_{z}(r, \phi, z, t) \propto J_{m}(k_{c}r)\cos(m\phi)e^{i(k_{z}z - \omega t)}` with `J_{m}(k_{c}R) = 0` (boundary condition). TE modes have `B_{z}` characteristic with `J_{m}'(k_{c}R) = 0`. The cutoff transverse wavenumber `k_{c}` is determined by zeroes of Bessel functions (or their derivatives), and the dispersion relation is `k_{z}^{2} = \omega^{2}/c^{2} - k_{c}^{2}`.

#### Proper-time reformulation

Source-free wave equation in cylindrical geometry: identical structure under `\mathbf{u} = 0`. Bessel-function modes and cutoff frequencies unchanged.

**Verdict:** ✅ identical.

**Companion notebook:** [`Roadmapping/Mathematica_Notebooks/Electromagnetism/JacksonCh08_P8_2.wl`](../../Mathematica_Notebooks/Electromagnetism/JacksonCh08_P8_2.wl).

---

### Problem J3e-P8.5 — Resonant cavity modes

**Selection provenance:** standing-wave modes in a closed resonant cavity; the eigenfrequency spectrum is discrete. Jackson 3e §8.7. *Pragmatic AI use.*

**Paraphrased statement:** Compute the resonant-frequency spectrum of a closed rectangular cavity of dimensions `a \times b \times d`. Identify the lowest mode.

#### Classical solution

Boundary conditions at all six walls produce a discrete spectrum

$$
\omega_{lmn} = c\sqrt{\left(\frac{l\pi}{a}\right)^{2} + \left(\frac{m\pi}{b}\right)^{2} + \left(\frac{n\pi}{d}\right)^{2}}.
$$

The lowest mode is TE_{101} (for `a > d > b` ordering convention). The Q-factor is determined by ohmic losses in the cavity walls — see [Problem J3e-P8.7](#problem-j3e-p87--power-loss-and-q-factor-in-waveguide).

#### Proper-time reformulation

Source-free standing waves; `\mathbf{u} = 0`. Frequency spectrum unchanged.

**Verdict:** ✅ identical.

**Companion notebook:** [`Roadmapping/Mathematica_Notebooks/Electromagnetism/JacksonCh08_P8_5.wl`](../../Mathematica_Notebooks/Electromagnetism/JacksonCh08_P8_5.wl).

---

### Problem J3e-P8.7 — Power loss and Q factor in waveguide

**Selection provenance:** the standard treatment of finite-conductivity wall losses in waveguides; Jackson 3e §8.5. *Pragmatic AI use.*

**Paraphrased statement:** Compute the attenuation rate of a TE_{10} mode in a rectangular waveguide with finite-conductivity walls. Express the result in terms of the skin depth `\delta`.

#### Classical solution

The wall surface currents dissipate energy at rate `dP/dz \propto 1/\delta`, where the skin depth `\delta = c/\sqrt{2\pi\omega\sigma\mu}` (Gaussian) was computed in [Problem J3e-P7.2](Ch07_Plane_EM_Waves.md#problem-j3e-p72--dispersion-in-a-conducting-medium). The attenuation coefficient `\alpha` is then proportional to the wall-dissipation rate divided by the propagating power, giving a frequency-dependent expression that increases with `1/\sqrt\omega` at high frequencies. Q-factor `Q = \omega/(2\alpha v_{g})` characterises the cavity's quality.

#### Proper-time reformulation

Wall-dissipation analysis depends on the skin depth of the wall material (from Ch. 7) and on the field configuration in the waveguide (which is unchanged). The proper-time formulation reduces identically since both inputs are unaffected.

**Verdict:** ✅ identical.

**Companion notebook:** [`Roadmapping/Mathematica_Notebooks/Electromagnetism/JacksonCh08_P8_7.wl`](../../Mathematica_Notebooks/Electromagnetism/JacksonCh08_P8_7.wl).
