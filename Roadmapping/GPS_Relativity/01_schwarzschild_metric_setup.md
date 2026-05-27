# Effect 01 — Schwarzschild metric setup (post-Newtonian to order 1/c²)

**One-line summary:** The foundational line element from which every clock-rate effect in this campaign is derived. The Schwarzschild solution outside a spherical Earth, expanded to order `1/c²`, gives a single master expression `dτ/dt ≈ 1 − GM/(rc²) − v²/(2c²)` that is the parent of effects 02–04 and an essential input to effects 05–08.

**Status:** ✅ drafted.
**Ashby (2003) cross-reference:** §3, Eqs. (12)–(18); MTW *Gravitation* §40.

## 1. Effect statement

The geometry outside the Earth is modeled as the Schwarzschild solution of Einstein's equations for a non-rotating, spherically symmetric mass. Two independent corrections are flagged but not folded into the leading-order analysis: (a) the Kerr / Lense-Thirring rotational gravitomagnetic correction from Earth's angular momentum `J⊕` (entering `g_{0φ}`) is `~ 10⁻¹⁶` per orbit and is flagged in effect 08 §6 as below the operational floor; (b) the static `J₂` quadrupole correction to the Earth's mass distribution (entering `g_{00}`) is the subject of effect 08 itself. These are independent effects with different physical origins (rotation vs mass distribution) and they are not "absorbed into" each other. The setup here is the spherically symmetric weak-field metric, expanded to the order at which all subsequent clock-rate calculations are carried out.

The post-Newtonian expansion truncates at `1/c²` — corrections at `1/c⁴` are at the `10⁻²⁰` level for GPS (products of two `1/c²` factors, each `~ 10⁻¹⁰` at GPS) and are below the operational precision floor.

## 2. Setup

Coordinates `(t, r, θ, φ)` are Schwarzschild coordinates, with the Earth's center at `r = 0`. The metric is:

$$ds^2 = -\left(1 - \frac{2GM}{rc^2}\right) c^2 \, dt^2 + \left(1 - \frac{2GM}{rc^2}\right)^{-1} dr^2 + r^2 \left(d\theta^2 + \sin^2\!\theta \, d\varphi^2\right).$$

The coordinate time `t` is **Geocentric Coordinate Time** (TCG): the time of a hypothetical clock at rest at infinity, after rescaling to remove the constant rate offset relative to clocks on the geoid. Operationally, GPS uses a different time scale (GPS time, derived from TAI), but the rate analysis is invariant under that choice.

| Symbol | Value | Meaning |
|---|---|---|
| `c` | `299 792 458 m/s` | speed of light (exact) |
| `GM⊕` | `3.986 004 418 × 10¹⁴ m³/s²` | geocentric gravitational constant |
| `r` | variable | radial coordinate |
| `v` | variable | coordinate-frame 3-velocity magnitude |
| `M⊕` | `5.972 × 10²⁴ kg` | Earth mass (informational; `GM⊕` is the directly-measured quantity) |

## 3. Derivation

The proper time `τ` of a clock following a worldline is

$$d\tau^2 = -\frac{1}{c^2} ds^2 = \left(1 - \frac{2GM}{rc^2}\right) dt^2 - \frac{1}{c^2}\left(1 - \frac{2GM}{rc^2}\right)^{-1} dr^2 - \frac{r^2}{c^2}\left(d\theta^2 + \sin^2\!\theta \, d\varphi^2\right).$$

Dividing by `dt²` and writing `v²` for the coordinate 3-velocity squared,

$$v^2 = \left(1 - \frac{2GM}{rc^2}\right)^{-1}\!\!\left(\frac{dr}{dt}\right)^{\!2} + r^2\!\left[\left(\frac{d\theta}{dt}\right)^{\!2} + \sin^2\!\theta \left(\frac{d\varphi}{dt}\right)^{\!2}\right],$$

we obtain the exact form

$$\left(\frac{d\tau}{dt}\right)^{\!2} = 1 - \frac{2GM}{rc^2} - \frac{v^2}{c^2}\left(1 + \mathcal{O}\!\left(\frac{2GM}{rc^2}\right)\right).$$

The first parenthesis on the right is `1 − 2GM/(rc²) ≈ 1 − 1.4 × 10⁻⁹` at the Earth's surface and `1 − 3.3 × 10⁻¹⁰` at GPS altitude. The second parenthesis, when expanded, contributes corrections of order `(GM/rc²)(v²/c²) ≈ 10⁻²⁰`, which we drop. To order `1/c²`:

$$\boxed{\;\frac{d\tau}{dt} \;\approx\; 1 \;-\; \frac{GM}{rc^2} \;-\; \frac{v^2}{2c^2}\;}$$

This is the master equation. The `−GM/(rc²)` term yields the gravitational time dilation of effect 02; the `−v²/(2c²)` term yields the velocity time dilation of effect 03; their combination at GPS-orbit parameters yields the net headline of effect 04. Effects 05 and 08 are corrections to the same expression for non-circular orbits and oblate gravitational fields; effects 06 and 07 are propagation-side corrections that do not appear in `dτ/dt` directly.

## 4. Numerical evaluation

The order-of-magnitude estimates for each term at GPS-relevant locations:

| Location | `GM/(rc²)` | `v²/(2c²)` |
|---|---|---|
| Geoid (sea level) | `6.953 × 10⁻¹⁰` | `1.20 × 10⁻¹²` (at equator) |
| GPS orbit (`r = 26 560 km`) | `1.670 × 10⁻¹⁰` | `8.35 × 10⁻¹¹` (circular orbit) |

Note that at GPS altitude `GM/(rc²) ≈ 2·v²/(2c²)` — the gravitational and kinematic effects are of the same order of magnitude, and they nearly cancel for low-Earth orbit; at GPS altitude gravity dominates by a factor of two, giving the net positive offset of effect 04. At the GPS orbital radius these two ratios stand in the algebraic relation `GM/(rc²) = 2·v²/(2c²)` exactly, because for a circular orbit `v² = GM/r`.

## 5. Wolfram MCP check

```wolfram
cc=299792458; GMe=3.986004418*^14; RR=6378137; aGPS=26560000; om=7.2921151467*^-5; Print["GM/(R c^2) = ", ScientificForm[N[GMe/(RR cc^2),8]]]; Print["GM/(a c^2) = ", ScientificForm[N[GMe/(aGPS cc^2),8]]]; Print["v_GPS = ", N[Sqrt[GMe/aGPS],8], " m/s"]; Print["v_eq = ", N[om RR,6], " m/s"]; Print["v_GPS^2/(2 c^2) = ", ScientificForm[N[GMe/(2 aGPS cc^2),8]]]; Print["v_eq^2/(2 c^2) = ", ScientificForm[N[(om RR)^2/(2 cc^2),8]]]
```

**Result:**
```
GM/(R c^2) = 6.9534893×10⁻¹⁰
GM/(a c^2) = 1.6698142×10⁻¹⁰
v_GPS = 3873.9621 m/s
v_eq = 465.101 m/s
v_GPS^2/(2 c^2) = 8.3490712×10⁻¹¹
v_eq^2/(2 c^2) = 1.2034088×10⁻¹²
```
✅ matches §4.

## 6. Comparison with Ashby (2003)

Ashby §3 develops the post-Newtonian expansion in the same form. His Eq. (12) is the line element above; his Eq. (16) is the master expression we boxed. Ashby retains an additional `+U·v²/c⁴` cross-term in Eq. (16) for completeness; that term is `~10⁻²⁰` for GPS and is dropped from all subsequent operational equations. Our truncation at `1/c²` matches Ashby's operational treatment in §4 onwards.

The sign convention `g_{00} = −(1 − 2GM/rc²)` (mostly-plus metric) matches Ashby and MTW; Weinberg uses the opposite sign convention but his clock-rate formulas, after the convention conversion, reduce to the same boxed result.

## 7. Verdict

✅ reproduces Ashby (2003) §3 Eq. (16) at the order at which all subsequent effects are calculated. The master equation is the parent of effects 02–04 (directly) and 05–08 (as the unperturbed starting point).
