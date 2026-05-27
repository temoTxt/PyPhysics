# Effect 02 — Gravitational time dilation

**One-line summary:** A clock in the weaker gravitational potential of GPS orbit ticks *faster* than a clock on the geoid by `+45.65 μs/day`. This is the larger of the two terms in the headline ±38 μs/day.

**Status:** ✅ drafted.
**Ashby (2003) cross-reference:** §3.1, Eq. (19); MTW §40.4.

## 1. Effect statement

A clock at higher gravitational potential ticks faster than a clock at lower potential. The Earth's gravitational potential is more negative near the surface (deeper in the well) than at GPS altitude (`≈ 20 182 km` above sea level), so a satellite clock — *purely from this effect* — runs faster than a ground clock. The frequency excess accumulates to `+45.65 μs` per day.

The effect is computed from the gravitational part of the master equation in effect 01, dropping the velocity term entirely (that term is the subject of effect 03; the two are combined in effect 04). The "ground" reference is the geoid — Earth's equipotential surface coinciding with mean sea level — because all clocks on the geoid tick at the same rate, regardless of latitude, by definition of the geoid.

## 2. Setup

Two clocks at rest in the Schwarzschild geometry: one at the geoid potential `Φ_geoid` (effective, including the centrifugal term from Earth's rotation), one at GPS orbital radius `a_GPS = 26 560 km` with the gravitational potential `Φ(a_GPS) = −GM⊕/a_GPS` (centrifugal terms not applicable in the non-rotating ECI frame in which we work).

| Symbol | Value | Meaning |
|---|---|---|
| `Φ_geoid/c²` | `−6.969 290 134 × 10⁻¹⁰` | IAU 2000 / W₀ value (combined gravitational + centrifugal) |
| `Φ(a_GPS)/c²` | `−1.669 814 × 10⁻¹⁰` | gravitational potential at GPS orbit |

Note that the geoid value `Φ_geoid` *already includes* the centrifugal contribution `−(1/2)ω²R²cos²λ` from Earth's rotation — this is what makes the geoid potential constant on the geoid surface despite the surface itself being non-spherical. So when we compare the GPS satellite's gravitational potential to `Φ_geoid`, we are correctly comparing potentials in their respective rest frames.

## 3. Derivation

From effect 01, the master equation for a clock at rest in the Schwarzschild geometry (velocity term dropped):

$$\frac{d\tau}{dt} \;\approx\; 1 \;-\; \frac{GM}{rc^2}.$$

For the geoid, the relevant potential is the *effective* potential `Φ_geoid` that combines the gravitational and centrifugal contributions; for the satellite in ECI (non-rotating) coordinates, the relevant potential is purely gravitational. The rate-difference *from the gravitational term alone* is:

$$\left.\frac{d\tau_{\rm sat}}{dt}\right|_{\rm grav} - \left.\frac{d\tau_{\rm grnd}}{dt}\right|_{\rm grav} \;=\; -\frac{GM}{a_{\rm GPS}\, c^2} \;-\; \frac{\Phi_{\rm geoid}}{c^2}\bigg|_{\rm grav-only}.$$

The "grav-only" piece of the geoid potential is `−GM/R⊕ + (J₂ correction)`; the centrifugal piece of `Φ_geoid` is part of effect 03's analysis. Decomposing:

$$\Phi_{\rm geoid} \;=\; \underbrace{-\,\frac{GM}{R_\oplus}}_{\rm pointmass} \;+\; \underbrace{\text{(J}_2\text{ correction)}}_{\text{effect 08}} \;-\; \underbrace{\frac{1}{2}\, \omega_\oplus^2 R_\oplus^2 \cos^2\!\lambda}_{\rm centrifugal,\ part\ of\ effect\ 03}.$$

Effect 02 (this document) takes the *pure gravitational* piece, deferring J₂ to effect 08 and the centrifugal contribution to effect 03. To leading order in the point-mass model:

$$\boxed{\;\left.\frac{\Delta f}{f}\right|_{\rm grav} \;=\; \frac{GM_\oplus}{R_\oplus c^2} \;-\; \frac{GM_\oplus}{a_{\rm GPS}\, c^2}\;}$$

The sign is positive: the satellite ticks faster than the ground clock by this fractional rate.

## 4. Numerical evaluation

Substituting:

$$\frac{GM_\oplus}{R_\oplus c^2} = 6.953\,489\,3 \times 10^{-10}, \qquad \frac{GM_\oplus}{a_{\rm GPS}\, c^2} = 1.669\,814\,2 \times 10^{-10}.$$

Difference:

$$\left.\frac{\Delta f}{f}\right|_{\rm grav} = 5.283\,675\,1 \times 10^{-10}.$$

Per day (`86 400 s`):

$$\left.\Delta t\right|_{\rm grav,\,daily} = +5.2837 \times 10^{-10} \times 86\,400 \;\approx\; +4.565 \times 10^{-5}\ {\rm s} \;=\; \boxed{\;+45.65\ \mu{\rm s/day}\;}.$$

The satellite gains `45.65 μs` per day relative to the ground from gravitational redshift alone, before any velocity-dilation correction.

## 5. Wolfram MCP check

```wolfram
cc=299792458; GMe=3.986004418*^14; RR=6378137; aGPS=26560000; day=86400; gravRate = GMe/(RR cc^2) - GMe/(aGPS cc^2); Print["grav rate = ", ScientificForm[N[gravRate,8]]]; Print["grav microsec/day = ", N[gravRate*day*1*^6,6]]
```

**Result:**
```
grav rate = 5.2836751×10⁻¹⁰
grav microsec/day = 45.6509
```
✅ matches §4.

## 6. Comparison with Ashby (2003)

Ashby §3.1 derives the gravitational redshift in his Eq. (19), giving the same `GM/c² × (1/R⊕ − 1/a)` form. His value `+45.7 μs/day` matches ours to within the rounding-precision of his quoted parameters (Ashby uses `GM⊕ = 3.9860044 × 10¹⁴`, a 7-figure value that rounds our 10-figure GM to the same 5-figure final answer).

A common alternative decomposition in operational GPS documents combines the geoid potential `Φ_geoid` directly:

$$\left.\frac{\Delta f}{f}\right|_{\rm grav} = \frac{\Phi_{\rm geoid}}{c^2} - \frac{\Phi(a_{\rm GPS})}{c^2}.$$

With `Φ_geoid/c² = −6.969 290 × 10⁻¹⁰` (the IAU 2000 value with J₂ and centrifugal absorbed) and `Φ(a_GPS)/c² = −1.669 814 × 10⁻¹⁰`, this gives `5.299 476 × 10⁻¹⁰` → `+45.787 μs/day`. The `0.14 μs/day` difference between this and our `45.65` decomposes into three contributions (using effect 08's numbers): centrifugal at equator (`+0.10 μs/day`, carried into effect 03's `v_eq²/(2c²)` term); J₂ at equator (`+0.033 μs/day`, accounted for in effect 08); and small residuals (`~0.01 μs/day`) from higher zonal harmonics and tidal averages folded into the empirical IAU 2000 W₀ value. **The two decompositions agree on the total — they differ in where each sub-contribution is bookkept.**

## 7. Verdict

✅ reproduces Ashby (2003) §3.1 Eq. (19) at quoted precision. The `+45.65 μs/day` figure assumes the point-mass `−GM/R⊕` geoid; the operational `+45.79 μs/day` figure absorbs the centrifugal contribution into the geoid value. The two are reconciled in effect 04.
