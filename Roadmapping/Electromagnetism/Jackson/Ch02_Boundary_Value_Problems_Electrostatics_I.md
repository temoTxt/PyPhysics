# Ch. 2 — Boundary-Value Problems in Electrostatics I

This chapter contains Jackson canonical problems on boundary-value electrostatics, worked in the proper-time reformulation alongside their classical CGS and SI solutions. Per [§4 of the plan](../../../.dev/tasks/42-electromagnetism-jackson-proper-time.md#4-unit-system-handling), Ch. 2 of Jackson is presented in the three-system regime: both CGS and SI versions of each problem appear, followed by the proper-time reformulation.

The Ch. 2 problems below are part of PR 0 (the fluency-warm-up prequel — see [§7.3 of the plan](../../../.dev/tasks/42-electromagnetism-jackson-proper-time.md#73-pr-a-prequel-adopted)). As in Ch. 1, the source charges are at rest in the observer's frame, so `u = 0` and `b = c` identically throughout the configuration. The proper-time formulation reduces to the classical formulation; the work of each problem is in the CGS↔SI translation, in the geometric content of the image method, and in the documentation that the reduction is identical.

## Problems in this chapter

| Problem | Status | Role |
|---|---|---|
| [Problem J3e-P2.1 — Point charge near a grounded conducting plane](#problem-j3e-p21--point-charge-near-a-grounded-conducting-plane) | drafted (PR 0) | fluency-builder |

---

### Problem J3e-P2.1 — Point charge near a grounded conducting plane

**Selection provenance** (Crocco §5 substantive-AI note):
- *Chosen because:* the method of images, established by Kelvin in the mid-nineteenth century and treated systematically in [[jackson1998_classical_electrodynamics]] Ch. 2, is the standard geometric technique for electrostatics with conductor boundaries. The grounded-plane case is the simplest configuration where the method applies non-trivially, and the resulting expressions (induced surface charge, image force, electrostatic energy) are textbook quantities recorded in every graduate-level treatment.
- *Alternatives considered:* J3e-P2.2 (point charge above an infinite grounded plane with a hemispherical boss — a richer image configuration but more complex geometry) and J3e-P2.6 (point charge near a grounded conducting sphere — an alternative image-method exercise). The simpler grounded-plane case is the more appropriate fluency-builder for PR 0; the sphere problem is a natural extension once the plane case is settled.
- *Role in this PR:* fluency-builder. The configuration is static, so the proper-time formulation reduces to the classical formulation; this problem exercises the image method (a geometric technique) and the surface-charge integration (a non-trivial volumetric integral over the plane).

<!-- TODO: human reviews and fills in — confirms the Kelvin attribution and the alternative-problem reasoning before this problem is rolled into PR 0's commit -->

**Source:** Jackson, *Classical Electrodynamics*, 3e Problem 2.1 (and 2e Problem 2.1, equivalent). *Paraphrased; consult the textbook for the precise statement.*

**Paraphrased statement:** A point charge `q` is located at a perpendicular distance `d` above an infinite grounded conducting plane. Find (i) the electrostatic potential everywhere above the plane, (ii) the induced surface charge density on the plane, (iii) the force on the point charge, and (iv) the total electrostatic energy of the configuration.

**Setup:** Place the conducting plane at `z = 0`, with conductor occupying `z ≤ 0` and the charge `q` at `(0, 0, d)` with `d > 0`. Both the charge and the plane are at rest in the observer's frame, so `w = 0`, `u = 0`, and `b² = c² + u² = c²` throughout. By the symmetry of the configuration we use cylindrical coordinates `(ρ, ϕ, z)` with `ρ² = x² + y²`. Write `r₁ = \sqrt{ρ² + (z - d)²}` for the distance from the field point `(ρ, ϕ, z)` to the source charge, and `r₂ = \sqrt{ρ² + (z + d)²}` for the distance to the image position `(0, 0, -d)`.

#### (a) Classical solution — Gaussian (CGS)

The method of images, in the form recorded by Jackson 3e Eqs. (2.1)–(2.4), replaces the conducting plane by an image charge `-q` placed at the mirror position `(0, 0, -d)`. In the half-space `z > 0`, the potential of the original configuration (charge `q` plus grounded plane) is identical to that of the two-point-charge configuration (`q` at `(0, 0, d)`, `-q` at `(0, 0, -d)`) — both satisfy Laplace's equation with the same boundary condition `Φ(z = 0) = 0` and the same source at `(0, 0, d)`. By the uniqueness theorem for Laplace's equation with Dirichlet boundary conditions, the two potentials must agree in `z > 0`. We therefore write, for `z > 0`,

$$
\Phi(\rho, z) = \frac{q}{r_{1}} - \frac{q}{r_{2}} = \frac{q}{\sqrt{\rho^{2} + (z - d)^{2}}} - \frac{q}{\sqrt{\rho^{2} + (z + d)^{2}}}.
$$

It is easy to show that `Φ(z = 0) = 0` identically in `ρ`, since `r₁ = r₂ = \sqrt{ρ² + d²}` on the plane.

The induced surface charge density on the conductor is recovered from the discontinuity in the normal component of the electric field, which in Gaussian units satisfies `σ = E_n/(4π)` with `E_n` measured immediately outside the conductor in the outward direction. Differentiating the potential,

$$
E_{z}(\rho, 0^{+}) = -\,\frac{\partial \Phi}{\partial z}\bigg|_{z = 0} = -\,\frac{2 q d}{(\rho^{2} + d^{2})^{3/2}},
$$

so that

$$
\sigma(\rho) = \frac{E_{z}(\rho, 0^{+})}{4\pi} = -\,\frac{q\,d}{2\pi (\rho^{2} + d^{2})^{3/2}}.
$$

The induced charge is negative everywhere on the plane (provided `q > 0`), as one expects from electrostatic attraction.

The total induced charge is obtained by integrating `σ` over the plane. This is the cleanest verification of the image method's bookkeeping, since the total induced charge must equal `-q` by global charge balance:

**Mathematica check** (Wolfram MCP, 2026-05-24):

```mathematica
ClearAll[rho, capD, capQ];
sigma = -capQ capD/(2 Pi (rho^2 + capD^2)^(3/2));
totalInduced = Integrate[sigma 2 Pi rho, {rho, 0, Infinity}, Assumptions -> capD > 0];
FullSimplify[totalInduced]
(* Result: -capQ  ✅ *)
```

The induced surface charge integrates to exactly `-q`, confirming that the image method reproduces the global charge balance.

The force on the point charge follows directly from the image: the field at `(0, 0, d)` produced by the image charge `-q` at `(0, 0, -d)` is a Coulomb field of magnitude `q/(2d)²` directed toward the plane. The force on the real charge is therefore

$$
\mathbf{F} = -\,\frac{q^{2}}{4 d^{2}}\,\hat z,
$$

attractive, directed toward the plane.

An independent computation, using the Maxwell stress tensor evaluated on the conductor surface, gives the force on the conductor itself: the electrostatic pressure pulling the conductor toward `q` is `σ²/(8\pi)` per unit area (in Gaussian units; equivalently, `E_n²/(8\pi)` outward from the conductor). Integrating over the plane:

```mathematica
ClearAll[rho, capD, capQ];
eField = -2 capQ capD/(rho^2 + capD^2)^(3/2);
pressure = eField^2/(8 Pi);
forceOnConductor = Integrate[pressure 2 Pi rho, {rho, 0, Infinity}, Assumptions -> capD > 0];
FullSimplify[forceOnConductor]
(* Result: capQ^2/(4 capD^2)  ✅ *)
```

The force on the conductor is `q²/(4d²)` in the `+\hat z` direction (toward `q`), and by Newton's third law the force on `q` is `-q²/(4d²)\hat z` (toward the conductor). The Maxwell-stress computation and the image-charge computation agree, as the consistency of the method of images requires.

The total electrostatic energy stored in the field above the conducting plane is

$$
W = -\,\frac{q^{2}}{4 d},
$$

which one obtains either by direct volume integration of `E²/(8\pi)` over the half-space `z > 0` (this integral is structurally identical to the integral for a dipole, restricted to one half-space), or, more economically, by integrating the image-method force `-q²/(4r²)` from `r = \infty` to `r = d`. We observe that this is *half* of the energy of the dipole `(q, -q)` with separation `2d`, namely `-q²/(2d)`, since only the half-space above the conductor contributes; the energy stored "inside" the conductor — that is, what one would have to associate with the image charge in a fully two-charge picture — is not physical, because no field exists below `z = 0`.

#### (b) Classical solution — SI

In SI units, Coulomb's law carries an additional `1/(4\pi\varepsilon_{0})` factor, and the energy density is `(\varepsilon_{0}/2) E²` rather than `E²/(8\pi)`. The image-method geometry is unchanged; the boundary condition `Φ(z = 0) = 0` and the uniqueness theorem apply equally well. We therefore obtain, by direct substitution into the expressions of (a):

$$
\Phi(\rho, z) = \frac{q}{4\pi\varepsilon_{0}}\!\left[\frac{1}{\sqrt{\rho^{2} + (z - d)^{2}}} - \frac{1}{\sqrt{\rho^{2} + (z + d)^{2}}}\right],
$$

$$
\sigma(\rho) = -\,\frac{q\,d}{2\pi (\rho^{2} + d^{2})^{3/2}},
$$

$$
\mathbf{F} = -\,\frac{1}{4\pi\varepsilon_{0}}\,\frac{q^{2}}{4 d^{2}}\,\hat z,
$$

$$
W = -\,\frac{1}{4\pi\varepsilon_{0}}\,\frac{q^{2}}{4 d}.
$$

The induced surface charge density is identical to the Gaussian expression — a notable feature of Gaussian electrostatics is that quantities expressible as `1/(4\pi)` times a length-cubed-inverse density (such as `σ`) take the same form in both systems, because the `4\pi` factor in Gauss's law of one system is absorbed by the `1/(4\pi\varepsilon_{0})` of the other. The total induced charge remains exactly `-q` in either system.

#### (c) Proper-time reformulation

As in [Problem J3e-P1.5](Ch01_Introduction_Electrostatics.md#problem-j3e-p15--electrostatic-self-energy-of-a-uniformly-charged-sphere), the source velocity is zero and the configuration is static. It follows from the velocity duality of Eq. (1) of [[Two_Mathematically_Equivalent_Versions_of_Maxwells_Equations]] that `u = 0`, and from `b² = c² + u²` that `b = c` throughout the configuration. Each substitution rule of [`_proper_time_cheatsheet.md`](../_proper_time_cheatsheet.md) reduces trivially:

- `c → b` is the identity, since `b = c`.
- The proper-time derivative `(1/b) ∂_τ` and the observer-time derivative `(1/c) ∂_t` are both zero, since the configuration is static; the substitution `(1/c) ∂_t \to (1/b) ∂_τ` therefore equates two vanishing quantities.
- The dissipative term `∂_τ(1/b) = -(\mathbf{u}·\mathbf{a})/b^{3}` vanishes by `\mathbf{u} = 0`.
- The current-density rescaling `\mathbf{J} → (b/c)\mathbf{J}` is the identity, since there is no current and `b/c = 1`.

The proper-time predictions for the potential, the induced surface charge density, the force on `q`, and the electrostatic energy are therefore identical to the classical Gaussian results of (a). We do not repeat them here.

It is worth observing that this identical reduction holds for *every* purely electrostatic problem in Jackson Chs. 1–3: the static character of the configuration forecloses any difference between the two clocks of [[Two_Mathematically_Equivalent_Versions_of_Maxwells_Equations]]. The two formulations are mathematically equivalent and physically equivalent in this case. The information encoded by the local clock that distinguishes the two formulations in general (Eq. (4) of the Maxwell paper) is zero whenever `\mathbf{u} = 0`; one cannot make a measurement on a static configuration that distinguishes the local clock of the source from the observer's clock, because the two clocks tick at the same rate when the source is at rest.

**Comparison:**

| Quantity | Classical (CGS) | Classical (SI) | Proper-time |
|---|---|---|---|
| Total induced charge | `-q` | `-q` | identical to CGS |
| Force on `q` | `-q^{2}/(4 d^{2})\hat z` | `-q^{2}/(4\pi\varepsilon_{0} \cdot 4 d^{2})\hat z` | identical to CGS |
| Energy | `-q^{2}/(4 d)` | `-q^{2}/(16 \pi\varepsilon_{0} d)` | identical to CGS |

**Does the proper-time answer differ from a pure `c → b` redressing?** ✅ no. With `\mathbf{u} = 0` we have `b = c`, and the `c → b` redressing is the identity.

**Verdict:** ✅ all three solutions consistent. The proper-time formulation contains no machinery that engages in the static case; it reproduces the classical answer exactly. The Maxwell-stress and image-method computations of the force on `q` agree, as the consistency of the image method requires.

**Notes for author review:** none. This is the second canary problem in PR 0; the template and voice have held up under a problem with more derivational content than [Problem J3e-P1.5](Ch01_Introduction_Electrostatics.md#problem-j3e-p15--electrostatic-self-energy-of-a-uniformly-charged-sphere). Two independent Mathematica MCP checks (induced charge and Maxwell-stress force) are recorded inline, as the verification convention requires.

**Companion notebook:** [`Roadmapping/Mathematica_Notebooks/Electromagnetism/JacksonCh02_P2_1.wl`](../../Mathematica_Notebooks/Electromagnetism/JacksonCh02_P2_1.wl) — runnable independent of the Mathematica MCP.

---

## PR K supplementary problems (added 2026-05-24)

Three Ch. 2 supplementary problems added in PR K's final-backfill batch. All static, `u = 0`, `b = c`; proper-time reduces identically to classical.

### Problem J3e-P2.3 — Point charge above grounded sphere

**Statement:** Point charge `q` at distance `d > R` from centre of grounded conducting sphere of radius `R`. Compute the image charge and the induced surface-charge distribution.

**Classical:** Image charge `q' = -qR/d` at distance `R^{2}/d` from centre. Surface charge density and force on `q` follow from the image construction.

**Proper-time:** Static, identical. **Verdict:** ✅. **Companion:** [`JacksonCh02_P2_3.wl`](../../Mathematica_Notebooks/Electromagnetism/JacksonCh02_P2_3.wl).

### Problem J3e-P2.5 — Two grounded planes meeting at angle

**Statement:** Two grounded conducting half-planes meet at a right angle. A point charge `q` lies in the quadrant. Compute the image-charge configuration.

**Classical:** Three image charges at the reflections of `q` through each plane and through the corner. Sum-of-Coulomb-fields potential.

**Proper-time:** Static, identical. **Verdict:** ✅. **Companion:** [`JacksonCh02_P2_5.wl`](../../Mathematica_Notebooks/Electromagnetism/JacksonCh02_P2_5.wl).

### Problem J3e-P2.7 — Conducting cylinder in uniform field

**Statement:** A long conducting cylinder of radius `R` is placed in a previously-uniform external field perpendicular to its axis. Compute the induced surface charge.

**Classical:** Two-dimensional separation of variables in cylindrical coordinates. Induced linear dipole per unit length `\lambda_{p} = 2E_{0}R^{2}`; field outside is `\mathbf{E}_{0}` plus pure-2D-dipole field.

**Proper-time:** Static, identical. **Verdict:** ✅. **Companion:** [`JacksonCh02_P2_7.wl`](../../Mathematica_Notebooks/Electromagnetism/JacksonCh02_P2_7.wl).
