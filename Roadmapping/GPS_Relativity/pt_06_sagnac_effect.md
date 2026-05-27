# Effect pt_06 — proper-time companion to Sagnac effect

**One-line summary:** The Sagnac effect derivation runs through the flat-space Maxwell equations in a rotating frame — exactly the framework's verified domain. Under the PT substitution rules, the photon-null condition gives `Δt_Sagnac = −2ω·A/b²` (where `b` is evaluated for the photon's path in the rotating frame). At GPS-relevant velocities `b ≈ c`, the prediction reduces to the standard `Δt_Sagnac = −2ω·A/c²` with relative error `~ 10⁻¹⁰`, matching the standard `~137 ns` maximum.

**Status:** ✅ matches standard at GPS precision.
**Framework applicability:** SR + flat-space Maxwell (direct framework application).
**Standard-method companion:** [`06_sagnac_effect.md`](06_sagnac_effect.md).

## 0. Framework applicability

The Sagnac effect is derived from the null condition `ds² = 0` for a photon in the rotating ECEF frame, with the off-diagonal Coriolis-like term `g_{0i} = (ω⊕ × r)_i / c`. The off-diagonal term comes from the coordinate transformation `φ_ECEF = φ_ECI − ω⊕ t` applied to the *flat-space* Minkowski metric. This is fully within the framework's verified domain — no GR-extension is needed.

The PT framework's substitution `c → b`, `(1/c)∂_t → (1/b)∂_τ` applies directly to the photon-null condition and the rotating-frame metric. The derivation parallels the standard one with `c → b`; the result reduces to the standard expression at `b → c`.

This document is in the **pragmatic** Crocco category — direct framework application of verified rules.

## 1. Effect statement

In the rotating ECEF frame, photons satisfy a null condition modified by the Coriolis-like off-diagonal metric term. Integrating along the photon's path gives a coordinate-time correction `Δt = L/c − 2ω·A/c²` (standard), where `A` is the area swept by the position vector and `L` is the Euclidean path length.

Under the PT framework, the photon-null condition uses `b` instead of `c` in the off-diagonal term, giving `Δt = L/b − 2ω·A/b²`. At GPS, where the satellite's velocity gives `b/c − 1 ≈ 8.3 × 10⁻¹¹`, the PT prediction matches the standard to relative precision `~10⁻¹⁰` — yielding the same `~137 ns` maximum.

A subtle point: the "satellite's `b`" is the framework's effective speed of light *as measured in the satellite's rest frame for its own motion*. For the photon (which is massless), the framework's `b` is the speed of light along the photon's worldline; for a photon in vacuum, this is just `c` (no proper-velocity correction). The PT correction to the Sagnac formula thus only appears when the *receiver* (which is rotating with Earth, so has `u_recv ≈ 465 m/s`) interprets the integrated signal — and that correction is at the `(u_recv/c)² ~ 10⁻¹²` level.

## 2. Setup

| Symbol | Value | PT-relevant note |
|---|---|---|
| `ω⊕` | `7.292 × 10⁻⁵ rad/s` | Earth sidereal rotation |
| `c` | `299 792 458 m/s` | speed of light in vacuum |
| `b_recv` (ground receiver) | `c (1 + 1.20 × 10⁻¹²)` | rotation gives `u_recv ≤ 465 m/s` |
| `A_E` (max projected area) | `8.47 × 10¹³ m²` | satellite-ground triangle, max |
| `2ω⊕/c²` (standard coefficient) | `1.623 × 10⁻²¹ s/m²` | matches `2ω⊕/b²` to `~10⁻¹²` |

The PT correction relative to standard is `(c²/b²) − 1 = −u²/b² ≈ −1.2 × 10⁻¹²` at the receiver — meaning the PT Sagnac is *smaller* than standard by 0.000 000 17 ns at maximum — far below GPS noise floor.

## 3. Derivation

The ECEF metric, under the framework's `c → b` substitution applied to the time-time and off-diagonal components, becomes:

$$ds^2 \;=\; -\!\left(1 - \frac{2\Phi}{b^2}\right) b^2 dt^2 \;+\; 2\, b\, (\boldsymbol{\omega}_\oplus \times \mathbf{r}) \cdot d\mathbf{x}\, dt \;+\; |d\mathbf{x}|^2.$$

For the photon (null condition `ds² = 0`), solving the quadratic for `b\, dt`:

$$b\, dt \;\approx\; |d\mathbf{x}| \;-\; \frac{1}{b}\, (\boldsymbol{\omega}_\oplus \times \mathbf{r}) \cdot d\mathbf{x}.$$

Integrating along the photon's path and using the vector-area identity `∮(ω×r)·dr = ω·2A`:

$$\boxed{\;\Delta t_{\rm Sagnac,\,PT} \;=\; -\,\frac{2\, \boldsymbol{\omega}_\oplus \cdot \mathbf{A}}{b^2}.\;}$$

This is the PT form. The structure differs from standard only in `b² → c²` substitution: in the standard formulation, `c` is the universal speed; in PT, `b` is the local effective speed. For the *receiver's frame* (where the Sagnac correction is observed), `b² = c² + u_recv² ≈ c²` for any Earth-rotation velocity.

## 4. Numerical evaluation

The PT and standard predictions for the maximum Sagnac correction:

$$|\Delta t_{\rm Sagnac}|_{\rm max} \;=\; \frac{2 \omega_\oplus A_{E,\,\max}}{b^2_{\rm recv}}.$$

With `b_recv² = c²(1 + 1.20 × 10⁻¹²)` and `A_{E,max} = 8.47 × 10¹³ m²`:

$$|\Delta t_{\rm Sagnac,\,PT}|_{\rm max} \;=\; \frac{1.623 \times 10^{-21}}{1 + 1.20 \times 10^{-12}} \times 8.47 \times 10^{13} \;\approx\; 1.374 \times 10^{-7}\ {\rm s} \;\approx\; \boxed{\;137\ {\rm ns}.\;}$$

The PT prediction differs from the standard `137.447 ns` by `~3.3 × 10⁻¹⁰ ns` (the `u_recv²/c²` correction); operationally invisible at any conceivable GPS clock-stability floor.

## 5. Wolfram MCP check

```wolfram
ClearAll[cc,om,RR,aGPS,Amax,vRecv,bRecv,PT,std]; cc=299792458; om=7.2921151467*^-5; RR=6378137; aGPS=26560000; Amax=aGPS RR/2; vRecv=om RR; bRecv=Sqrt[cc^2 + vRecv^2]; PT=2 om Amax/bRecv^2; std=2 om Amax/cc^2; Print["std Sagnac max (s) = ", ScientificForm[N[std, 14]]]; Print["PT Sagnac max (s) = ", ScientificForm[N[PT, 14]]]; Print["PT - std (s) = ", ScientificForm[N[PT - std, 8]]]; Print["PT - std (ns) = ", ScientificForm[N[(PT - std) 1*^9, 8]]]; Print["max Sagnac (ns) = ", N[PT 1*^9, 8]]
```

**Result:**
```
std Sagnac max (s) = 1.3744661×10⁻⁷
PT Sagnac max (s) = 1.3744661×10⁻⁷
PT - std (s) = -3.30793×10⁻¹⁹
PT - std (ns) = -3.30793×10⁻¹⁰
max Sagnac (ns) = 137.4466
```
✅ matches standard. PT correction relative to standard is `~3 × 10⁻¹⁰ ns` (operationally invisible).

## 5b. GPS-precision-limit equivalence to standard

The PT-to-standard reduction is the simplest in the campaign:

$$\frac{1}{b_{\rm recv}^2} \;=\; \frac{1}{c^2 + u_{\rm recv}^2} \;=\; \frac{1}{c^2}\,\frac{1}{1 + u_{\rm recv}^2/c^2} \;\approx\; \frac{1}{c^2}\!\left(1 - \frac{u_{\rm recv}^2}{c^2}\right).$$

The leading term gives the standard `2ω·A/c²`. The PT correction is at `𝒪(u²/c⁴)` — for the GPS receiver (rotating at `≤ 465 m/s`), this is `~ 10⁻²² s/m²` × `A`, giving sub-picosecond corrections to the Sagnac.

For the satellite (which experiences the *full* Sagnac through its own non-zero `u_sat`), the contribution to ECI-to-ECEF transformations is `~ 10⁻¹⁰` of standard — operationally just below the GPS clock-stability floor. *In principle*, this is a measurement-accessible regime for the next generation of optical-clock GPS receivers — but it's not a current measurement.

## 6. Comparison with standard-method companion

Standard [`06_sagnac_effect.md`](06_sagnac_effect.md) §3 derives the Sagnac correction from `g_{0i} = (ω×r)_i/c` in the ECEF metric. The PT framework replaces `c → b` in the same off-diagonal term; the algebraic structure of the derivation is identical, and the numerical result agrees to `~10⁻¹⁰` at GPS precision.

The framework's prediction is structurally distinct from standard in that `b` is *local-frame-dependent* (different `b` at satellite vs receiver). For GPS, both frames have `b ≈ c` so the difference is invisible. For a clock comparison involving *different rotating frames* (e.g., one rotating with Earth, one rotating with Jupiter), the framework's prediction would explicitly carry the local `b` of each frame.

## 6b. Where the framework would diverge

The PT prediction for Sagnac would differ from standard at observable levels when:

1. **The receiver moves at relativistic speed in its rotating frame.** A spacecraft Sagnac interferometer rotating at near-light-speed would carry a PT correction of order `(u/c)²` to the receiver-frame Sagnac formula.
2. **The light source has non-trivial proper motion.** For atomic-clock signals from a uniformly-moving source, the framework's `b` for the source frame contributes a small additional Sagnac term. For GPS, this is the satellite's `b_sat ≈ c(1 + 8.3 × 10⁻¹¹)` — giving a correction at the `10⁻¹⁰` level relative to standard.

Neither regime is currently a GPS measurement, but the second one is at the edge of next-generation optical-clock GPS systems.

## 7. Verdict

✅ matches standard-method companion [`06_sagnac_effect.md`](06_sagnac_effect.md) §7 at GPS precision. PT framework directly applies (no GR extension needed); the structure differs only in `c → b` at the off-diagonal term, with `b ≈ c` for GPS receivers reducing the prediction to the standard form. PT correction relative to standard is `~ 10⁻¹⁰ ns` — operationally invisible.
