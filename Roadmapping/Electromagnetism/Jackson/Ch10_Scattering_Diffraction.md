# Ch. 10 — Scattering and Diffraction

This chapter contains Jackson canonical problems on EM scattering and diffraction, worked in the proper-time reformulation alongside their classical solutions. Part of **PR J (backfill)** per [§7 of the plan](../../../.dev/tasks/42-electromagnetism-jackson-proper-time.md#7-initial-chapter-selection--canonical-problems-list). Compact format.

## Problems in this chapter

| Problem | Status | Role |
|---|---|---|
| [Problem J3e-P10.1 — Thomson scattering from a free electron](#problem-j3e-p101--thomson-scattering-from-a-free-electron) | drafted (PR J) | fluency-builder |
| [Problem J3e-P10.5 — Rayleigh scattering by small dielectric sphere](#problem-j3e-p105--rayleigh-scattering-by-small-dielectric-sphere) | drafted (PR J) | fluency-builder |
| [Problem J3e-P10.7 — Kirchhoff diffraction from circular aperture](#problem-j3e-p107--kirchhoff-diffraction-from-circular-aperture) | drafted (PR J) | fluency-builder |

---

### Problem J3e-P10.1 — Thomson scattering from a free electron

**Selection provenance:** the simplest scattering process; non-relativistic limit. *Pragmatic AI.*

**Paraphrased statement:** Compute the differential and total cross-section for Thomson scattering — elastic scattering of low-energy photons off a free electron.

#### Classical solution

Thomson differential cross-section:

$$
\frac{d\sigma}{d\Omega} = r_{e}^{2}\,\frac{1 + \cos^{2}\theta}{2}, \quad r_{e} = \frac{e^{2}}{m c^{2}} \approx 2.82 \times 10^{-13}\,\text{cm}.
$$

Total cross-section: `\sigma_{T} = (8\pi/3)r_{e}^{2} \approx 6.65 \times 10^{-25}` cm² (the Thomson cross-section).

#### Proper-time reformulation

Electron oscillates at the incident photon frequency with amplitude `x_{0} \sim eE_{0}/(m\omega^{2})`. For non-relativistic incident-photon energies (`\hbar\omega \ll m c^{2}`), `u/c \ll 1` and the third-term correction is `O((u/c)^{2}) \ll 1`. Thomson cross-section reproduces exactly at leading order.

**Verdict:** ✅ identical at leading order. At high incident-photon energies (`\hbar\omega \sim m c^{2}`, Compton regime), the classical Thomson formula breaks down and QED takes over; the proper-time formulation's structural validity in that regime is the same open question as in PR D.

**Companion notebook:** [`Roadmapping/Mathematica_Notebooks/Electromagnetism/JacksonCh10_P10_1.wl`](../../Mathematica_Notebooks/Electromagnetism/JacksonCh10_P10_1.wl).

---

### Problem J3e-P10.5 — Rayleigh scattering by small dielectric sphere

**Selection provenance:** the canonical small-particle scattering problem; explains "why the sky is blue". *Pragmatic AI.*

**Paraphrased statement:** A small dielectric sphere of radius `a \ll \lambda` and refractive index `n` scatters incident EM radiation. Compute the differential cross-section in the Rayleigh limit.

#### Classical solution

The induced dipole moment of the sphere is `\mathbf{p} = (n^{2}-1)/(n^{2}+2) \cdot a^{3}\mathbf{E}_{0}`. The radiated power follows from the dipole-radiation formula (J3e-P9.1), giving differential cross-section

$$
\frac{d\sigma}{d\Omega} = k^{4} a^{6}\left(\frac{n^{2}-1}{n^{2}+2}\right)^{2}\,\frac{1 + \cos^{2}\theta}{2}.
$$

Total: `\sigma \propto k^{4} a^{6}` — the famous `1/\lambda^{4}` dependence that produces blue-sky scattering.

#### Proper-time reformulation

Bound-electron oscillation in the dielectric sphere; `u/c \ll 1`. Same `O((u/c)^{2})` third-term structure as J3e-P9.1, below observational floor.

**Verdict:** ✅ identical at leading order.

**Companion notebook:** [`Roadmapping/Mathematica_Notebooks/Electromagnetism/JacksonCh10_P10_5.wl`](../../Mathematica_Notebooks/Electromagnetism/JacksonCh10_P10_5.wl).

---

### Problem J3e-P10.7 — Kirchhoff diffraction from circular aperture

**Selection provenance:** the classical Fraunhofer-diffraction result; standard in optics. *Pragmatic AI.*

**Paraphrased statement:** Compute the Fraunhofer-zone diffraction pattern of a plane EM wave passing through a circular aperture of radius `a` in an opaque screen.

#### Classical solution

The Kirchhoff integral evaluates in the Fraunhofer limit to the Airy pattern:

$$
\frac{dP}{d\Omega} \propto \left|\frac{2 J_{1}(k a \sin\theta)}{k a \sin\theta}\right|^{2},
$$

with the first minimum at `\sin\theta = 1.22\,\lambda/(2a)` — the Rayleigh resolution criterion.

#### Proper-time reformulation

Diffraction is a wave-propagation phenomenon with no source motion at the aperture; `\mathbf{u} = 0`, identical to classical.

**Verdict:** ✅ identical.

**Companion notebook:** [`Roadmapping/Mathematica_Notebooks/Electromagnetism/JacksonCh10_P10_7.wl`](../../Mathematica_Notebooks/Electromagnetism/JacksonCh10_P10_7.wl).
