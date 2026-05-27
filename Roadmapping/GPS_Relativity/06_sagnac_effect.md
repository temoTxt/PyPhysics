# Effect 06 — Sagnac effect (rotating-frame light propagation)

**One-line summary:** When the GPS signal is processed in the Earth-Centered Earth-Fixed (ECEF) rotating frame, the speed of light is *not* `c` along a generic path — there is a path-dependent correction `Δt_Sagnac = (2ω⊕/c²) · A_E`, where `A_E` is the projected area swept by the photon path in the equatorial plane. For GPS this contributes up to `±137 ns` per signal, with the sign set by whether the signal travels in the same or opposite direction as Earth's rotation.

**Status:** ✅ drafted.
**Ashby (2003) cross-reference:** §3.4, Eqs. (29)–(35).

## 1. Effect statement

GPS pseudoranges are computed in the ECEF (rotating) frame, the natural frame for a ground-based receiver. Maxwell's equations have a constant speed of light `c` only in inertial frames; in a rotating frame, the line element of the Minkowski metric picks up an off-diagonal term `g_{0i} ∝ (ω × r)_i / c` (the GR analogue of a Coriolis effect on light), and the coordinate time of flight for a photon depends on its path through that off-diagonal field.

The correction is the GR-language version of the **Sagnac effect** measured by Sagnac (1913), the same effect used in fiber-optic gyroscopes. For Earth, the magnitude scales as `2ω⊕/c² ≈ 1.62 × 10⁻²¹ s/m²` × projected area. For a GPS-to-ground signal, the maximum projected area is `~a_GPS × R⊕/2 ≈ 8.5 × 10¹³ m²`, giving a maximum per-signal Sagnac correction of `~137 ns`.

The correction is *path-dependent* — it cannot be absorbed into either a constant frequency offset or a per-satellite ephemeris term. It must be computed by the receiver from the signal geometry.

## 2. Setup

The ECEF metric, to order `1/c²`, is obtained from the ECI Schwarzschild metric by the coordinate transformation `φ_ECEF = φ_ECI − ω⊕ t`:

$$ds^2_{\rm ECEF} = -\left(1 - \frac{2\Phi(r)}{c^2}\right) c^2 dt^2 + 2\, c\, (\boldsymbol{\omega}_\oplus \times \mathbf{r}) \cdot d\mathbf{x}\, dt + |d\mathbf{x}|^2,$$

where `(ω⊕ × r)` is the velocity field of an ECEF-rest point as seen from ECI. The new off-diagonal term `g_{0i} = (ω⊕ × r)_i / c` is the source of the Sagnac correction; the diagonal terms recover the ECI line element on the diagonal in the no-rotation limit.

| Symbol | Value | Meaning |
|---|---|---|
| `ω⊕` | `7.292 × 10⁻⁵ rad/s` | Earth sidereal rotation rate |
| `2ω⊕/c²` | `1.623 × 10⁻²¹ s/m²` | Sagnac coefficient |
| `A_E` | path-dependent (m²) | projected equatorial area swept by photon |

## 3. Derivation

A photon travels on a null geodesic, `ds² = 0`. From the ECEF metric:

$$0 = -\left(1 - \frac{2\Phi}{c^2}\right) c^2 dt^2 + 2 c\, (\boldsymbol{\omega} \times \mathbf{r}) \cdot d\mathbf{x}\, dt + |d\mathbf{x}|^2.$$

Solving the quadratic for `c\,dt` to leading order in `1/c²`:

$$c\, dt \;\approx\; |d\mathbf{x}| \;-\; \frac{1}{c}\, (\boldsymbol{\omega}_\oplus \times \mathbf{r}) \cdot d\mathbf{x} \;+\; \mathcal{O}(c^{-3}).$$

Integrating along the photon's path from transmitter at `r_T` to receiver at `r_R`:

$$c\, \Delta t \;=\; \int |d\mathbf{x}| \;-\; \frac{1}{c} \int (\boldsymbol{\omega}_\oplus \times \mathbf{r}) \cdot d\mathbf{x} \;+\; \mathcal{O}(c^{-3}).$$

The first integral is the Euclidean path length `L`; the second integral, using the identity `(ω × r) · dr = ω · (r × dr)`, becomes `ω · (∮(1/2) r × dr) × 2 = ω · 2A`, where `A` is the vector area swept by the position vector along the path. Hence:

$$\Delta t \;=\; \frac{L}{c} \;-\; \frac{2 \boldsymbol{\omega}_\oplus \cdot \mathbf{A}}{c^2}.$$

The first term is the naive light-travel time at speed `c`; the second is the Sagnac correction. Since `ω⊕` points along the Earth's rotation axis (the `ẑ`-axis in ECEF), `ω · A = ω⊕ A_E`, where `A_E` is the projection of the vector area onto the equatorial plane:

$$\boxed{\;\Delta t_{\rm Sagnac} \;=\; -\,\frac{2 \omega_\oplus}{c^2}\, A_E\;}$$

The sign: a westbound signal (against Earth's rotation) has `A_E < 0` and `Δt_Sagnac > 0` (the signal takes *longer*); an eastbound signal has `A_E > 0` and `Δt_Sagnac < 0`. This is the same direction-dependent asymmetry as the original Sagnac interferometer.

## 4. Numerical evaluation

The Sagnac coefficient is:

$$\frac{2 \omega_\oplus}{c^2} \;=\; \frac{2 \times 7.292\,115 \times 10^{-5}}{(2.998 \times 10^8)^2} \;=\; 1.623 \times 10^{-21}\ {\rm s/m^2}.$$

The maximum projected area for a GPS signal: the triangle whose vertices are the Earth's center, the GPS satellite, and the ground receiver. When this triangle is in the equatorial plane and at maximum extent, `A_E ≈ a_GPS · R⊕ / 2 ≈ 8.47 × 10¹³ m²`. Substituting:

$$\big|\Delta t_{\rm Sagnac}\big|_{\rm max} \;=\; 1.623 \times 10^{-21} \times 8.47 \times 10^{13} \;\approx\; 1.374 \times 10^{-7}\ {\rm s} \;=\; \boxed{\;\approx 137\ {\rm ns}.\;}$$

Operational GPS literature commonly quotes "`~133 ns` maximum"; the small discrepancy with our `137` is because the geometric maximum is constrained by the visibility horizon (the satellite must be above the receiver's local horizon to be received), not the full half-disk.

A `137 ns` per-signal correction corresponds to a pseudorange error of `c × 137.4 ns ≈ 41.2 m` if uncorrected — substantially larger than the meter-scale civilian GPS accuracy. Every GPS receiver computes the Sagnac correction from the broadcast satellite ephemeris and the receiver's known ECEF position; this is one of the standard receiver-side corrections in any GPS firmware.

## 5. Wolfram MCP check

```wolfram
cc=299792458; om=7.2921151467*^-5; aGPS=26560000; RR=6378137; Amax = aGPS RR/2; Print["Sagnac coefficient (s/m^2) = ", ScientificForm[N[2 om/cc^2,8]]]; Print["max projected area (m^2) = ", ScientificForm[N[Amax]]]; Print["max Sagnac (ns) = ", N[2 om/cc^2 Amax 1*^9,6]]; (* Pseudorange error *) Print["pseudorange error (m) = ", N[cc 2 om/cc^2 Amax,6]]
```

**Result:**
```
Sagnac coefficient (s/m^2) = 1.6227145×10⁻²¹
max projected area (m^2) = 8.4701659×10¹³
max Sagnac (ns) = 137.4466
pseudorange error (m) = 41.2055
```
✅ matches §4. (Note: an earlier version of this block quoted `1.624×10⁻²¹` / `137.587 ns` / `41.249 m` — those values do not reproduce from a clean Wolfram session with the parameter values in this campaign's README, and have been corrected.)

## 6. Comparison with Ashby (2003)

Ashby §3.4 derives the Sagnac correction in the same form (his Eq. (29) for the metric, Eq. (35) for the integrated correction). His maximum value `133 ns` is the visibility-horizon-limited value for a receiver at `~45°` latitude; our `137 ns` is the geometric maximum without the visibility constraint. The two agree once the visibility constraint is folded in.

The Sagnac correction is structurally distinct from effects 02–05: those are *clock-rate* corrections (modifying `dτ/dt`); the Sagnac is a *signal-propagation* correction (modifying `Δt_received − Δt_emitted` for fixed clock rates). It belongs to the same family as effect 07 (Shapiro) — both are corrections to the propagation time, not the clock rate.

Note also that the Sagnac correction is *frame-dependent*: in ECI, it is absent (light travels at `c` everywhere); it appears only when one chooses to work in the rotating ECEF frame. ECEF is operationally easier (the receiver's coordinates are stationary in it), which is why the correction is computed. There is no physical disagreement between ECI and ECEF analyses — they give the same observable arrival time.

## 7. Verdict

✅ reproduces Ashby (2003) §3.4 Eq. (35) at quoted precision. The geometric maximum `~137 ns` matches the operational `~133 ns` visibility-limited maximum to within the geometric refinement. Sagnac correction is implemented in every GPS receiver.
