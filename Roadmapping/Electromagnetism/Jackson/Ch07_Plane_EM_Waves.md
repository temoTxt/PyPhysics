# Ch. 7 — Plane Electromagnetic Waves and Wave Propagation

This chapter contains Jackson canonical problems on plane EM wave propagation in vacuum and linear media, worked in the proper-time reformulation alongside their classical solutions. Per [§4 of the plan](../../../.dev/tasks/42-electromagnetism-jackson-proper-time.md#4-unit-system-handling), Ch. 7 of Jackson is presented in the three-system regime.

Ch. 7 is part of **PR F (backfill chapters)** per [§7 of the plan](../../../.dev/tasks/42-electromagnetism-jackson-proper-time.md#7-initial-chapter-selection--canonical-problems-list). The problems below are fluency-builders that confirm the proper-time formulation reduces to classical EM in vacuum and linear-media wave propagation — there is no source motion in any of these problems, so `\mathbf{u} = 0` and `b = c`, making the dissipative `(\mathbf{u}\cdot\mathbf{a})/b^{4}` coefficient inactive throughout.

## Problems in this chapter

| Problem | Status | Role |
|---|---|---|
| [Problem J3e-P7.1 — Plane EM wave in vacuum](#problem-j3e-p71--plane-em-wave-in-vacuum) | drafted (PR F) | fluency-builder (foundational) |
| [Problem J3e-P7.2 — Dispersion in a conducting medium](#problem-j3e-p72--dispersion-in-a-conducting-medium) | drafted (PR F) | fluency-builder |
| [Problem J3e-P7.3 — Reflection and refraction at a plane boundary](#problem-j3e-p73--reflection-and-refraction-at-a-plane-boundary) | drafted (PR F) | fluency-builder |
| [Problem J3e-P7.5 — Polarisation states of a plane wave](#problem-j3e-p75--polarisation-states-of-a-plane-wave) | drafted (PR F) | fluency-builder |

---

### Problem J3e-P7.1 — Plane EM wave in vacuum

**Selection provenance:** the simplest possible EM wave solution; foundational for all of Ch. 7. *Pragmatic AI use only; no substantive interpretive content beyond the standard textbook derivation.*

**Source:** Jackson, *Classical Electrodynamics*, 3e Problem 7.1 (and 2e equivalent). *Paraphrased.*

**Paraphrased statement:** Derive the plane-wave solution to Maxwell's equations in vacuum, identify the dispersion relation `\omega = kc`, and verify that the proper-time formulation produces the identical result.

**Setup:** Free-space Maxwell equations, no source. Trial solution `\mathbf{E} = E_{0}\hat x\cos(kz - \omega t)`, `\mathbf{B} = E_{0}\hat y\cos(kz - \omega t)`.

#### (a)/(b) Classical solution (Gaussian and SI)

The wave equation `(1/c^{2})\partial^{2}\mathbf{E}/\partial t^{2} = \nabla^{2}\mathbf{E}` is satisfied for the trial solution provided `\omega = kc`. The transverse-wave condition `\mathbf{E}\cdot\hat k = 0` and `\mathbf{B} = \hat k\times\mathbf{E}` (Gaussian; with `c` factors in SI) determine the polarisation structure.

**Mathematica check** (Wolfram MCP, 2026-05-24):

```mathematica
ClearAll[capE0, kk, omega, cc, zz, tt];
eField = capE0 Cos[kk zz - omega tt];
residual = FullSimplify[(1/cc^2) D[eField, {tt, 2}] - D[eField, {zz, 2}]];
(* residual = (capE0 (cc^2 kk^2 - omega^2) Cos[omega tt - kk zz]) / cc^2
   Vanishes iff omega = k c  ✅ *)
```

#### (c) Proper-time reformulation

With no source, `\mathbf{u} = 0` ⇒ `b = c` throughout. The proper-time wave equation `(1/b^{2})\partial^{2}\mathbf{E}/\partial\tau^{2} = \nabla^{2}\mathbf{E}` reduces to the classical form. Dispersion `\omega = kc` is unchanged.

**Verdict:** ✅ identical. Vacuum plane-wave propagation is a null-result for the proper-time framework.

**Companion notebook:** [`Roadmapping/Mathematica_Notebooks/Electromagnetism/JacksonCh07_P7_1.wl`](../../Mathematica_Notebooks/Electromagnetism/JacksonCh07_P7_1.wl).

---

### Problem J3e-P7.2 — Dispersion in a conducting medium

**Selection provenance:** treats the modification of the dispersion relation by finite conductivity, producing complex wavenumber and skin-depth attenuation. *Pragmatic AI use.*

**Paraphrased statement:** Compute the dispersion relation for a plane EM wave in a conducting medium with conductivity `\sigma` and permittivity `\varepsilon`. Show that the wavenumber is complex, with imaginary part determining the skin depth.

#### Classical solution

For a wave `\mathbf{E} \propto e^{i(kz - \omega t)}` in a conductor obeying Ohm's law `\mathbf{J} = \sigma\mathbf{E}`, Maxwell's equations give the dispersion

$$
k^{2} = \frac{\omega^{2}\varepsilon\mu}{c^{2}} + i\frac{4\pi\omega\sigma\mu}{c^{2}}.
$$

For a good conductor (`\sigma \gg \omega\varepsilon/(4\pi)`), `k \approx (1+i)/\delta` with skin depth `\delta = c/\sqrt{2\pi\omega\sigma\mu}`. The wave decays exponentially over `\delta`.

#### Proper-time reformulation

The medium is at rest in the observer's frame; the conduction electrons have drift velocities far below `c`, so the macroscopic `\mathbf{J}` has `\mathbf{u} \approx \mathbf{v}_{\text{drift}} \ll c` and the `\rho\mathbf{u} = (b/c)\mathbf{J}` substitution is approximately the identity. The dispersion relation is unchanged.

**Verdict:** ✅ identical. The skin depth and complex wavenumber are the same in both formulations.

**Companion notebook:** [`Roadmapping/Mathematica_Notebooks/Electromagnetism/JacksonCh07_P7_2.wl`](../../Mathematica_Notebooks/Electromagnetism/JacksonCh07_P7_2.wl).

---

### Problem J3e-P7.3 — Reflection and refraction at a plane boundary

**Selection provenance:** Fresnel equations for EM-wave reflection / transmission at a planar interface between two linear media. *Pragmatic AI use.*

**Paraphrased statement:** Derive the Fresnel reflection and transmission coefficients for a plane EM wave incident at angle `\theta_{i}` on a planar boundary between two media with refractive indices `n_{1}` and `n_{2}`. Confirm that the proper-time formulation reproduces the standard Fresnel result.

#### Classical solution

For s-polarisation (`\mathbf{E}` perpendicular to plane of incidence):

$$
r_{s} = \frac{n_{1}\cos\theta_{i} - n_{2}\cos\theta_{t}}{n_{1}\cos\theta_{i} + n_{2}\cos\theta_{t}}, \qquad t_{s} = \frac{2 n_{1}\cos\theta_{i}}{n_{1}\cos\theta_{i} + n_{2}\cos\theta_{t}}.
$$

For p-polarisation (`\mathbf{E}` in plane of incidence):

$$
r_{p} = \frac{n_{2}\cos\theta_{i} - n_{1}\cos\theta_{t}}{n_{2}\cos\theta_{i} + n_{1}\cos\theta_{t}}, \qquad t_{p} = \frac{2 n_{1}\cos\theta_{i}}{n_{2}\cos\theta_{i} + n_{1}\cos\theta_{t}}.
$$

Snell's law `n_{1}\sin\theta_{i} = n_{2}\sin\theta_{t}` relates the angles.

#### Proper-time reformulation

Boundary conditions on `\mathbf{E}_{\parallel}` and `\mathbf{B}_{\parallel}` are the same in both formulations (derived from Maxwell's equations with no source velocity at the boundary). The Fresnel coefficients are unchanged.

**Verdict:** ✅ identical. Reflection / refraction at a planar boundary is unaffected by the proper-time formulation.

**Companion notebook:** [`Roadmapping/Mathematica_Notebooks/Electromagnetism/JacksonCh07_P7_3.wl`](../../Mathematica_Notebooks/Electromagnetism/JacksonCh07_P7_3.wl).

---

### Problem J3e-P7.5 — Polarisation states of a plane wave

**Selection provenance:** the standard linear / circular / elliptical polarisation decomposition for a monochromatic plane wave. *Pragmatic AI use.*

**Paraphrased statement:** Decompose a general plane wave's polarisation state into linear, circular, and elliptical components. Verify Stokes parameters and the structure of the Poincaré sphere.

#### Classical solution

A general plane wave can be written

$$
\mathbf{E}(\mathbf{r}, t) = (E_{0x}\hat x + E_{0y}\hat y\,e^{i\delta})\,e^{i(kz - \omega t)},
$$

with `E_{0x}`, `E_{0y}` real and `\delta` the relative phase. Stokes parameters `(I, Q, U, V)` characterise the polarisation state, with `I^{2} = Q^{2} + U^{2} + V^{2}` for fully-polarised light. The Poincaré sphere encodes the state geometrically.

#### Proper-time reformulation

Polarisation states are properties of the field configuration, independent of source motion. The proper-time formulation produces the same polarisation analysis as classical EM for any plane-wave configuration in vacuum or a stationary linear medium.

**Verdict:** ✅ identical. Polarisation analysis is unaffected by the proper-time formulation in source-free wave propagation.

**Companion notebook:** [`Roadmapping/Mathematica_Notebooks/Electromagnetism/JacksonCh07_P7_5.wl`](../../Mathematica_Notebooks/Electromagnetism/JacksonCh07_P7_5.wl).
