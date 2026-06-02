# Ch. 4 — Multipoles, Electrostatics of Macroscopic Media

This chapter contains Jackson canonical problems on the multipole expansion and electrostatics in linear dielectric media, worked in the proper-time reformulation alongside their classical CGS and SI solutions. Per [§4 of the plan](../../../.dev/tasks/42-electromagnetism-jackson-proper-time.md#4-unit-system-handling), Ch. 4 of Jackson is presented in the three-system regime.

Ch. 4 is part of **PR G (backfill)** per [§7 of the plan](../../../.dev/tasks/42-electromagnetism-jackson-proper-time.md#7-initial-chapter-selection--canonical-problems-list). Like PR F's Ch. 7, these problems are pure-fluency null-result confirmations: every Ch. 4 configuration is static or quasi-static, so `u = 0` and `b = c`, and the proper-time formulation reduces identically to classical EM. Compact per-problem format per the PR F precedent.

## Problems in this chapter

| Problem | Status | Role |
|---|---|---|
| [Problem J3e-P4.1 — Multipole expansion of potential](#problem-j3e-p41--multipole-expansion-of-potential) | drafted (PR G) | fluency-builder (foundational) |
| [Problem J3e-P4.3 — Polarization in a linear dielectric](#problem-j3e-p43--polarization-in-a-linear-dielectric) | drafted (PR G) | fluency-builder |
| [Problem J3e-P4.5 — Energy in dielectric medium](#problem-j3e-p45--energy-in-dielectric-medium) | drafted (PR G) | fluency-builder |
| [Problem J3e-P4.7 — Boundary conditions at dielectric interface](#problem-j3e-p47--boundary-conditions-at-dielectric-interface) | drafted (PR G) | fluency-builder |

---

### Problem J3e-P4.1 — Multipole expansion of potential

**Selection provenance:** the multipole expansion is the foundational tool for treating localised charge distributions at large distance; Jackson 3e §4.1. *Pragmatic AI use.*

**Paraphrased statement:** Derive the multipole expansion of the electrostatic potential due to a localised charge distribution, identifying the monopole, dipole, and quadrupole contributions.

#### Classical solution (Gaussian)

For a localised charge density `\rho(\mathbf{r}')` centred at the origin, the potential at field point `\mathbf{r}` with `r \gg r'_{\max}` is

$$
\phi(\mathbf{r}) = \frac{q}{r} + \frac{\mathbf{p}\cdot\hat r}{r^{2}} + \frac{1}{2 r^{3}}\sum_{i,j}Q_{ij}\hat r_{i}\hat r_{j} + O(r^{-4}),
$$

with `q = \int\rho\, dV`, `\mathbf{p} = \int\mathbf{r}'\rho\, dV`, `Q_{ij} = \int(3 r'_{i} r'_{j} - r'^{2}\delta_{ij})\rho\, dV`. SI form differs by an overall `1/(4\pi\varepsilon_{0})` factor.

#### Proper-time reformulation

Static configuration: `\mathbf{u} = 0`, `b = c`. The proper-time potential reduces to the classical expression exactly. The multipole moments are the same.

**Verdict:** ✅ identical. The multipole expansion is unaffected by the proper-time formulation for static charge distributions.

**Companion notebook:** [`Roadmapping/Mathematica_Notebooks/Electromagnetism/JacksonCh04_P4_1.wl`](../../Mathematica_Notebooks/Electromagnetism/JacksonCh04_P4_1.wl).

---

### Problem J3e-P4.3 — Polarization in a linear dielectric

**Selection provenance:** macroscopic-medium response of a linear dielectric to an external field; Jackson 3e §4.3. *Pragmatic AI use.*

**Paraphrased statement:** A linear dielectric of susceptibility `\chi_{e}` is placed in a uniform external field `\mathbf{E}_{0}`. Compute the polarization `\mathbf{P}`, the bound charge density `\rho_{b}`, and the displacement field `\mathbf{D}`.

#### Classical solution

For a linear dielectric, `\mathbf{P} = \chi_{e}\mathbf{E}` (Gaussian, where `\chi_{e}` is dimensionless). The displacement is `\mathbf{D} = \mathbf{E} + 4\pi\mathbf{P} = (1 + 4\pi\chi_{e})\mathbf{E} = \varepsilon\mathbf{E}` with `\varepsilon = 1 + 4\pi\chi_{e}`. Bound charges: volume `\rho_{b} = -\nabla\cdot\mathbf{P}`, surface `\sigma_{b} = \mathbf{P}\cdot\hat n`.

SI form: `\mathbf{P} = \chi_{e}\varepsilon_{0}\mathbf{E}`, `\mathbf{D} = \varepsilon_{0}\mathbf{E} + \mathbf{P} = \varepsilon_{0}(1+\chi_{e})\mathbf{E}`.

#### Proper-time reformulation

The bound charges in a stationary linear dielectric have zero drift velocity at the macroscopic level (the polarisation establishes itself essentially instantaneously and then doesn't move). `\mathbf{u} = 0` for the macroscopic bound-charge distribution, so `b = c` and the proper-time formulation reduces identically.

**Verdict:** ✅ identical. Linear-dielectric polarisation is unaffected by the proper-time formulation in static configurations.

**Companion notebook:** [`Roadmapping/Mathematica_Notebooks/Electromagnetism/JacksonCh04_P4_3.wl`](../../Mathematica_Notebooks/Electromagnetism/JacksonCh04_P4_3.wl).

---

### Problem J3e-P4.5 — Energy in dielectric medium

**Selection provenance:** the electrostatic energy stored in the field of a linear dielectric; Jackson 3e §4.7. *Pragmatic AI use.*

**Paraphrased statement:** Compute the total electrostatic energy stored in the field of a linear dielectric with permittivity `\varepsilon`. Compare to the vacuum case.

#### Classical solution

In a linear dielectric (Gaussian), the energy density is `u = (1/(8\pi))\mathbf{E}\cdot\mathbf{D} = (\varepsilon/(8\pi))E^{2}`. The total energy is `W = \int u\, dV`. The presence of the dielectric reduces the energy stored at given free-charge configuration by a factor `\varepsilon`, because the polarisation partially screens the field.

In SI: `u = (1/2)\mathbf{E}\cdot\mathbf{D} = (\varepsilon\varepsilon_{0}/2)E^{2}` (where `\varepsilon` here is the dimensionless relative permittivity).

#### Proper-time reformulation

Static dielectric: `\mathbf{u} = 0`, `b = c`. The energy expression is unchanged.

**Verdict:** ✅ identical. Energy in linear dielectrics is unaffected by the proper-time formulation.

**Companion notebook:** [`Roadmapping/Mathematica_Notebooks/Electromagnetism/JacksonCh04_P4_5.wl`](../../Mathematica_Notebooks/Electromagnetism/JacksonCh04_P4_5.wl).

---

### Problem J3e-P4.7 — Boundary conditions at dielectric interface

**Selection provenance:** the standard boundary conditions for the EM field at the interface between two linear dielectrics; Jackson 3e §4.4. *Pragmatic AI use.*

**Paraphrased statement:** Derive the boundary conditions on `\mathbf{E}` and `\mathbf{D}` at the interface between two linear dielectrics in the absence of free surface charges.

#### Classical solution

From Maxwell's equations and the divergence theorem (for `\mathbf{D}`) and Stokes' theorem (for `\mathbf{E}`):

$$
D_{1n} = D_{2n} \quad (\text{normal component of } \mathbf{D} \text{ continuous when no free surface charge}),
$$

$$
\mathbf{E}_{1t} = \mathbf{E}_{2t} \quad (\text{tangential component of } \mathbf{E} \text{ continuous}).
$$

Combining with the constitutive relations `\mathbf{D}_{i} = \varepsilon_{i}\mathbf{E}_{i}` gives `\varepsilon_{1}E_{1n} = \varepsilon_{2}E_{2n}`. The field bends at the interface — the analogue of Snell's law for static fields.

#### Proper-time reformulation

The boundary conditions on `\mathbf{E}` and `\mathbf{D}` derive directly from Maxwell's equations applied at the interface; they depend only on the field equations and not on source-velocity content. The proper-time formulation produces the same boundary conditions.

**Verdict:** ✅ identical. Dielectric-interface boundary conditions are unaffected by the proper-time formulation.

**Companion notebook:** [`Roadmapping/Mathematica_Notebooks/Electromagnetism/JacksonCh04_P4_7.wl`](../../Mathematica_Notebooks/Electromagnetism/JacksonCh04_P4_7.wl).
