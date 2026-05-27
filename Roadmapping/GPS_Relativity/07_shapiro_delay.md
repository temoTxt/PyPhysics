# Effect 07 — Shapiro propagation delay (gravitational time delay)

**One-line summary:** Light signals travel slightly *slower* through Earth's gravitational potential well than they would in flat space — the Shapiro delay. For GPS, the maximum per-signal correction is `~62 ps` for a satellite near the horizon, `~42 ps` for a satellite at zenith. The effect is negligible for civilian GPS positioning but is included in geodetic-precision GPS post-processing (e.g., IGS final products).

**Status:** ✅ drafted.
**Ashby (2003) cross-reference:** §6, Eqs. (53)–(54).

## 1. Effect statement

The Shapiro delay (Shapiro 1964) is the general-relativistic prediction that light passing through a gravitational potential takes longer than the Euclidean path length divided by `c`. The extra delay arises because the spatial part of the Schwarzschild metric is curved: at radius `r`, the effective spatial coordinate `dr_proper = dr/√(1 − 2GM/rc²) > dr`, so the photon traverses more proper distance per coordinate `dr`.

For Earth's mass, the prefactor `2GM⊕/c³ ≈ 30 ps` sets the scale. Multiplied by a logarithm of the geometric ratio of distances, the maximum GPS-link Shapiro delay is at the `~60 ps` level — small compared to all other corrections in this campaign, but not negligible for precise geodetic applications.

## 2. Setup

A null geodesic connects transmitter at `r_T` to receiver at `r_R`, with the photon's path of Euclidean length `L = |r_R − r_T|`. In the Schwarzschild geometry, the coordinate time of flight is `L/c` plus a positive correction `Δt_Shapiro` that depends on the geometry only through `r_T`, `r_R`, and `L`.

| Symbol | Value | Meaning |
|---|---|---|
| `2GM⊕/c³` | `2.957 × 10⁻¹¹ s` | Shapiro prefactor |
| `r_T = a_GPS` | `26 560 000 m` | satellite radial distance |
| `r_R = R⊕` | `6 378 137 m` | ground-receiver radial distance |
| `L` (slant range) | `~ 2.0–2.6 × 10⁷ m` | direct Euclidean distance |

The maximum slant range occurs when the satellite is just on the horizon: `L_max = √(r_T² − r_R²) ≈ 2.578 × 10⁷ m`. The minimum slant range occurs at zenith: `L_min = r_T − r_R ≈ 2.018 × 10⁷ m`.

## 3. Derivation

Starting from the Schwarzschild line element with `ds² = 0` for a photon, restricted to the orbital plane (`θ = π/2`) for clarity:

$$0 = -\left(1 - \frac{2GM}{rc^2}\right) c^2 dt^2 + \left(1 - \frac{2GM}{rc^2}\right)^{-1} dr^2 + r^2 d\varphi^2.$$

The coordinate time elapsed along the null geodesic between `r_T` and `r_R` is:

$$c\, \Delta t \;=\; \int_{\rm path} \frac{1}{\sqrt{1 - \frac{2GM}{rc^2}}} \cdot \frac{1}{1 - \frac{2GM}{rc^2}}\; dr \;+\; \text{(angular terms)}.$$

For weak field, expand to first order in `2GM/(rc²)`:

$$c\, \Delta t \;\approx\; \int_{\rm path} \!\left(1 \;+\; \frac{2GM}{rc^2}\right) ds_{\rm flat},$$

where `ds_flat` is the flat-space line element. The first term gives the Euclidean light-travel time `L/c`; the second term is the Shapiro correction.

For a straight-line path from `r_T` to `r_R` with closest-approach distance `b` (impact parameter), the integral evaluates to the standard result:

$$\boxed{\;\Delta t_{\rm Shapiro} \;=\; \frac{2 GM_\oplus}{c^3}\, \ln\!\left(\frac{r_T + r_R + L}{r_T + r_R - L}\right)\;}$$

where `L` is the Euclidean distance between transmitter and receiver. This is Ashby's Eq. (53) and the classical Shapiro formula (Misner-Thorne-Wheeler §40.4).

## 4. Numerical evaluation

For a GPS satellite at zenith over the receiver:

- `L_zenith = a_GPS − R⊕ = 2.018 × 10⁷ m`
- `r_T + r_R = 3.294 × 10⁷ m`
- `(r_T + r_R + L)/(r_T + r_R − L) = 5.312/1.276 = 4.162`
- `ln(4.162) = 1.427`
- `Δt_Shapiro = 2.957 × 10⁻¹¹ × 1.427 = 4.22 × 10⁻¹¹ s = 42.2 ps`

For a GPS satellite near the horizon (maximum slant range):

- `L_max = √(a_GPS² − R⊕²) = 2.578 × 10⁷ m`
- `(r_T + r_R + L)/(r_T + r_R − L) = 5.872/0.716 = 8.20`
- `ln(8.20) = 2.105`
- `Δt_Shapiro = 2.957 × 10⁻¹¹ × 2.105 = 6.23 × 10⁻¹¹ s = 62.3 ps`

The total maximum Shapiro delay per signal is `~62 ps`, corresponding to a pseudorange error of `c × 62 × 10⁻¹² ≈ 1.87 cm` if uncorrected. This is below the meter-scale civilian GPS noise floor but is included in geodetic post-processing (IGS final products) and in two-way satellite-to-satellite GPS time-transfer.

## 5. Wolfram MCP check

```wolfram
cc=299792458; GMe=3.986004418*^14; aGPS=26560000; RR=6378137; pre = 2 GMe/cc^3; Lz = aGPS - RR; Lh = Sqrt[aGPS^2 - RR^2]; ratz = (aGPS+RR+Lz)/(aGPS+RR-Lz); rath = (aGPS+RR+Lh)/(aGPS+RR-Lh); Print["2 GM/c^3 (s) = ", ScientificForm[N[pre,8]]]; Print["zenith L (m) = ", N[Lz,6]]; Print["zenith log ratio = ", N[Log[ratz],6]]; Print["zenith Shapiro (ps) = ", N[pre Log[ratz] 1*^12,6]]; Print["horizon L (m) = ", N[Lh,6]]; Print["horizon log ratio = ", N[Log[rath],6]]; Print["horizon Shapiro (ps) = ", N[pre Log[rath] 1*^12,6]]; Print["horizon pseudorange error (cm) = ", N[cc pre Log[rath] 100,6]]
```

**Result:**
```
2 GM/c^3 (s) = 2.9573902×10⁻¹¹
zenith L (m) = 2.018186×10⁷
zenith log ratio = 1.42653
zenith Shapiro (ps) = 42.2072
horizon L (m) = 2.578282×10⁷
horizon log ratio = 2.10494
horizon Shapiro (ps) = 62.2795
horizon pseudorange error (cm) = 1.86694
```
✅ matches §4.

## 6. Comparison with Ashby (2003)

Ashby §6 derives the Shapiro delay starting from the same null-geodesic integral. His Eq. (53) is identical to the boxed result above. He notes (Ashby §6, paragraph after Eq. (54)) that the Shapiro delay for an Earth-orbiting receiver "rarely exceeds 2 cm in equivalent path length"; our `~1.87 cm` for the horizon case matches.

A common point of confusion: the Shapiro delay is sometimes quoted at `~19 ps` for GPS. That smaller value refers to the *receiver-applied* correction in standard single-frequency GPS, which uses an approximate, geometry-averaged value rather than the per-signal exact formula. For high-precision applications (IGS, satellite laser ranging, geodesy), the exact per-signal value (`~42–62 ps`) is used.

The Shapiro delay differs structurally from the Sagnac (effect 06) in that it is *frame-independent* (depends only on the invariant Schwarzschild radial coordinates and the slant length) and *sign-definite* (always positive — light always takes longer). It is part of the same propagation-correction family but cannot be cancelled by choosing a different frame.

## 7. Verdict

✅ reproduces Ashby (2003) §6 Eq. (53) at quoted precision. The maximum per-signal value `~62 ps` matches the Ashby-quoted "rarely exceeds 2 cm" path-length equivalent. Implemented in geodetic-precision GPS post-processing; routinely omitted from civilian receiver firmware where its contribution is below the noise floor.
