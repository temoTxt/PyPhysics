# Effect pt_05 — proper-time companion to eccentricity correction

**One-line summary:** The eccentricity correction's main term `−(2√(GMa)/c²) e sin E` carries over to the PT framework with `c → b` substitution; the reduction to standard at GPS precision is the same `~10⁻²⁰` story as effects pt_01–pt_04. The proper-time framework's *third term* `(u·a)/b⁴ · ∂_τ` (Maxwell-paper Eq. (4) wave-equation dissipative coefficient) makes a separate contribution to the eccentricity correction; for circular orbits it vanishes identically (`u ⊥ a`); for elliptical GPS orbits at `e = 0.02` the contribution is `~ 10⁻³⁴` per signal — twenty orders of magnitude below GPS clock stability.

**Status:** ⚠ partial — main eccentricity term is standard; third-term contribution is direct framework application and is operationally negligible.
**Framework applicability:** mixed (orbital dynamics ✅; gravitational potential under speculative extension; third-term contribution direct).
**Standard-method companion:** [`05_eccentricity_correction.md`](05_eccentricity_correction.md).

## 0. Framework applicability

The main eccentricity correction is a Kepler-orbit integral over `dτ/dt = 1 − GM/(rc²) − v²/(2c²)`. The PT formulation replaces this with `1 − GM/(rb²) − u²/(2b²)`, which (per pt_01 §5b) reduces to the standard formula at GPS precision.

The PT framework additionally predicts a *third-term contribution* to any dynamical quantity that involves the satellite's velocity and acceleration: the dissipative coefficient `−(u·a)/b⁴ · ∂_τ` in the wave equation (Maxwell-paper Eq. (4)). For a satellite in orbit, `u` is the orbital velocity and `a` is the gravitational acceleration. **For circular orbits, `u ⊥ a` (velocity is tangential, acceleration is radial), so `u·a = 0` exactly and the third term vanishes.** For elliptical orbits, `u` has a radial component when the satellite is between perigee and apogee, and `u·a ≠ 0`. We compute this contribution explicitly and show it is operationally negligible for GPS.

This document is in the **pragmatic** Crocco category for the third-term derivation (direct framework application of a verified rule), and **substantive** for the speculative extension to the gravitational potential (inherited from pt_01).

## 1. Effect statement

The standard eccentricity correction `Δτ_ecc = −(2√(GMa)/c²) e sin E` arises from the periodic variation of `r` and `v` over an elliptical Kepler orbit. The PT framework reproduces this main term at GPS precision via the same `c → b` reduction as in pt_01–pt_04.

The framework additionally predicts a *new* contribution: the third-term `(u·a)/b⁴` coefficient in the wave equation governing EM signal propagation in the satellite's local frame. This is a framework-specific prediction, not present in standard SR + GR. The contribution is operationally `~10⁻³⁴ s` per signal at the maximum GPS eccentricity of `e = 0.02`, twenty orders of magnitude below GPS clock stability.

## 2. Setup

For an elliptical Keplerian orbit (parameters as in standard effect 05):

| Symbol | Value at GPS, `e = 0.02` | Meaning |
|---|---|---|
| `e` | `0.02` (max spec) | orbital eccentricity |
| `p = a(1−e²)` | `2.655 × 10⁷ m` | semi-latus rectum |
| `u_r` (max radial proper velocity) | `√(GM/p) × e ≈ 77.5 m/s` | at `f = π/2`, true anomaly |
| `a_r` (radial gravitational acc.) | `GM/r² ≈ 0.565 m/s²` | always radially inward |
| `u·a` (max) | `u_r × a_r ≈ 43.8 m²/s³` | only nonzero between perigee/apogee |
| `(u·a)/b⁴` | `≈ 5.42 × 10⁻³³ /s` | wave-equation dissipative coefficient |

The maximum `u·a` occurs at the orbital phase where the satellite is moving fastest radially (typically at true anomaly `f ≈ ±π/2`), and is zero at perigee and apogee.

## 3. Derivation

### 3a. Main eccentricity term under PT framework

Substituting `c → b` in the standard derivation:

$$\frac{d\tau}{dt}\bigg|_{\rm PT,\,ecc} \;=\; 1 \;-\; \frac{GM}{r b^2} \;-\; \frac{u^2}{2 b^2} \;=\; 1 \;-\; \frac{2GM}{r b^2} \;+\; \frac{GM}{2 a b^2}\quad\text{(vis-viva, with } u^2 = GM(2/r - 1/a)).$$

Reducing `b → c` at GPS precision (relative error `~10⁻¹⁰`):

$$\Delta\tau_{\rm ecc,\,PT}(t) \;\approx\; -\,\frac{2\sqrt{GM_\oplus a_{\rm GPS}}}{c^2}\, e\, \sin E(t) \;+\; \mathcal{O}(10^{-20})\, e\, \sin E.$$

**Result:** identical to standard at GPS precision. Peak `~46 ns` at `e = 0.02`, same as effect 05.

### 3b. Third-term contribution (framework-specific)

The wave equation Eq. (4) of the Maxwell paper, applied to a signal propagating in the satellite's local frame, carries a dissipative coefficient `−(u·a)/b⁴ · ∂_τ`. For a signal of duration `τ_sig` (the typical signal-emission timescale, `~ 1/f_GPS ≈ 7 × 10⁻⁸ s` for the L1 carrier; the relevant accumulated time for the satellite-to-ground signal is the full propagation time `~ 0.07 s`):

$$\Delta\tau_{\rm 3rd\,term} \;\sim\; \frac{(\boldsymbol{u} \cdot \boldsymbol{a})}{b^4}\, \tau_{\rm sig}.$$

For GPS at `e = 0.02`, maximum `|u·a| ≈ 43.8 m²/s³`, `b ≈ c`, `τ_sig ≈ 0.07 s`:

$$|\Delta\tau_{\rm 3rd\,term}|_{\rm max} \;\approx\; \frac{43.8}{(3 \times 10^8)^4} \times 0.07 \;\approx\; 3.8 \times 10^{-34}\ {\rm s}.$$

Twenty orders of magnitude below GPS clock stability (`~10⁻¹⁴`); operationally invisible.

## 4. Numerical evaluation

| Component | Magnitude (peak) | Relative to GPS clock stability |
|---|---|---|
| Main eccentricity term (PT = standard) | `±46 ns` | `~10⁻⁵` (operationally large) |
| Third-term contribution (PT-specific) | `~ 4 × 10⁻³⁴ s` per signal | `~10⁻²⁰` (operationally invisible) |

The main term is what GPS receivers compute and apply (every receiver firmware on Earth implements it). The third term is the framework's *additional* prediction; for GPS it is unobservable, but it would matter in regimes where `u·a` is much larger.

## 5. Wolfram MCP check

```wolfram
cc=299792458; GMe=3.986004418*^14; aGPS=26560000; ecc=0.02; p=aGPS(1-ecc^2); urMax=Sqrt[GMe/p] ecc; arMax=GMe/aGPS^2; uDotAmax=urMax arMax; tauSig=0.07; (* full satellite-to-ground propagation time *) Print["max u_r (m/s) = ", N[urMax, 6]]; Print["a_r (m/s^2) = ", N[arMax, 6]]; Print["max u·a (m^2/s^3) = ", N[uDotAmax, 6]]; Print["(u·a)/b^4 (per s) = ", ScientificForm[N[uDotAmax/cc^4, 8]]]; Print["max third-term per signal (s) = ", ScientificForm[N[uDotAmax/cc^4 tauSig, 8]]]; (* compare to main term *) mainTerm = 2 Sqrt[GMe aGPS]/cc^2 ecc; Print["main eccentricity term peak (s) = ", ScientificForm[N[mainTerm, 8]]]; Print["ratio third / main = ", ScientificForm[N[(uDotAmax/cc^4 tauSig)/mainTerm, 8]]]
```

**Result:**
```
max u_r (m/s) = 77.4947
a_r (m/s^2) = 0.565043
max u·a (m^2/s^3) = 43.7878
(u·a)/b^4 (per s) = 5.42089×10⁻³³
max third-term per signal (s) = 3.79462×10⁻³⁴
main eccentricity term peak (s) = 4.57933×10⁻⁸
ratio third / main = 8.29×10⁻²⁷
```
✅ confirms framework-specific third-term contribution is `~10⁻²⁶` of the main eccentricity term (operationally invisible).

## 5b. GPS-precision-limit equivalence to standard

The main term is identical to standard at `𝒪(c⁻²)` (reduction `b → c` at GPS precision is `~10⁻¹⁰` relative; the main term is `~10⁻⁸` so the relative error in PT prediction is `~10⁻¹⁸` of the main result — well below the noise floor).

The third-term contribution is operationally invisible (`~10⁻²⁶` relative to the main term). At GPS the framework reproduces the standard receiver-implemented correction without any observable addition.

## 6. Comparison with standard-method companion

Standard [`05_eccentricity_correction.md`](05_eccentricity_correction.md) §3 derives the peak `±46 ns` correction from `Δτ_ecc = −(2/c²)√(GMa) e sin E`. The PT framework reproduces this main term identically at GPS precision. The framework *additionally* predicts a third-term contribution from `(u·a)/b⁴`; for GPS this is `~10⁻²⁶` of the main term and is operationally invisible.

The PT framework's third-term contribution is the same dissipative coefficient that drives the radiation-reaction predictions in [issue #43](https://github.com/temoTxt/PyPhysics/issues/43) and in [Ch16_Radiation_Damping](../Electromagnetism/Jackson/Ch16_Radiation_Damping.md) of the Electromagnetism campaign. For GPS, the third term is operationally invisible; for radiation-reaction experiments with high acceleration (the Cole/Poder/Wistisen 2018 regime, `a ~ 10³⁰ m/s²`), the third term contributes at order unity and is potentially distinguishable. The framework is internally consistent — the same coefficient appears in both regimes.

## 6b. Where the framework would diverge

The third-term contribution to eccentricity-like periodic corrections would become operationally observable when `(u·a) τ_sig / b⁴ ~ 10⁻¹⁴` (the GPS clock-stability floor). For `τ_sig ~ 0.07 s` (one propagation time), this requires `u·a ~ 10⁻¹⁴ × 0.07⁻¹ × (3×10⁸)⁴ ≈ 1.2 × 10³³ m²/s³`. With `u ~ c` (relativistic), this means `a ~ 4 × 10²⁴ m/s²` — which is the acceleration scale of an electron in a `10²⁰ W/cm²` laser pulse. **This is exactly the [#43](https://github.com/temoTxt/PyPhysics/issues/43) experimental regime** — not GPS.

For lower accelerations, the framework's third-term contribution is below operational precision and the framework reproduces the standard eccentricity correction.

## 7. Verdict

⚠ matches standard-method companion [`05_eccentricity_correction.md`](05_eccentricity_correction.md) §7 at GPS precision for the main `±46 ns` term. The framework's third-term contribution `(u·a)/b⁴` is operationally invisible for GPS (`~10⁻³⁴ s` per signal at `e = 0.02`). The third term *would* become observable in high-acceleration regimes — the same regime as the [#43](https://github.com/temoTxt/PyPhysics/issues/43) radiation-reaction sub-investigation.
