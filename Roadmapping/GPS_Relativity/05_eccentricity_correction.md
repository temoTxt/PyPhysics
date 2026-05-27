# Effect 05 — Eccentricity correction (periodic `e·sin E` term)

**One-line summary:** GPS orbits are not perfectly circular — typical eccentricity `e ≈ 0.005–0.02`. The departure from circular orbit makes both `r` and `v` vary periodically over an orbit, producing a periodic time-dilation correction of amplitude `±2√(GMa)/c² × e` that reaches `±46 ns` at the maximum allowed eccentricity. This correction must be applied by the *receiver* (not absorbed into the satellite's pre-launch frequency offset) because it depends on the satellite's instantaneous orbital phase.

**Status:** ✅ drafted.
**Ashby (2003) cross-reference:** §4, Eqs. (40)–(45).

## 1. Effect statement

In an elliptical Keplerian orbit, the radial distance varies as `r = a(1 − e cos E)` and the orbital speed varies through the vis-viva relation `v² = GM(2/r − 1/a)`. The headline calculation in effect 04 used the *circular-orbit* values `r = a`, `v² = GM/a`; both vary periodically along an elliptical orbit. The deviation from those circular values, integrated over time, gives a periodic correction `Δt_ecc(E) = −(2/c²)√(GMa) · e · sin E`, where `E` is the eccentric anomaly.

This is **not** absorbed into the pre-launch frequency offset for two reasons:
1. It is **periodic**, not secular — its average over an orbit is zero, so a constant frequency offset cannot absorb it.
2. It depends on the satellite's orbital phase (`E`), which changes from satellite to satellite and orbit to orbit. The GPS ground segment uploads orbital ephemerides; the receiver computes `e sin E` from the broadcast ephemeris and applies the correction locally.

## 2. Setup

| Symbol | Value | Meaning |
|---|---|---|
| `a` | `26 560 000 m` | semi-major axis (Kepler) |
| `e` | `≤ 0.02` (operational) | orbital eccentricity (mission spec) |
| `E` | variable | eccentric anomaly |
| `M` | variable | mean anomaly, `M = n(t − t_p)` |
| `n` | `√(GM/a³)` | mean motion |
| `t_p` | epoch | time of perigee passage |
| `T_orb` | `2π/n ≈ 11.97 h` | orbital period (≈ half a sidereal day) |

Kepler's equation `M = E − e sin E` connects mean anomaly to eccentric anomaly. The true anomaly `f` (angle from perigee to satellite in the orbit plane) is related to `E` by `tan(f/2) = √((1+e)/(1−e)) tan(E/2)`.

## 3. Derivation

Starting from the master equation (effect 01), now with `r` and `v` allowed to vary:

$$\frac{d\tau}{dt} \;=\; 1 \;-\; \frac{GM}{r c^2} \;-\; \frac{v^2}{2c^2}.$$

Substituting the vis-viva relation `v² = GM(2/r − 1/a)`:

$$\frac{d\tau}{dt} \;=\; 1 \;-\; \frac{GM}{r c^2} \;-\; \frac{GM}{c^2}\left(\frac{1}{r} - \frac{1}{2a}\right) \;=\; 1 \;-\; \frac{2 GM}{r c^2} \;+\; \frac{GM}{2 a c^2}.$$

The mean rate over the orbit uses the time-average `⟨1/r⟩_t = 1/a` (a standard Kepler result, distinct from the angle-average which is `1/(a(1−e²)^{1/2})`):

$$\left\langle\frac{d\tau}{dt}\right\rangle_t \;=\; 1 \;-\; \frac{2GM}{a c^2} \;+\; \frac{GM}{2 a c^2} \;=\; 1 \;-\; \frac{3 GM}{2 a c^2}.$$

This matches the circular-orbit result of effect 04 — the mean rate of an elliptical orbit equals the circular-orbit rate at semi-major axis `a`, **independent of eccentricity**. This is what allows GPS to apply a single constant pre-launch frequency offset for all satellites in the constellation: the offset depends only on `a`, which is nominally the same.

The departure from the mean is the periodic correction:

$$\Delta\!\left(\frac{d\tau}{dt}\right) \;=\; -\frac{2 GM}{c^2}\!\left(\frac{1}{r} - \frac{1}{a}\right).$$

To first order in `e`: `1/r − 1/a ≈ (e/a) cos E + 𝒪(e²)`, so:

$$\Delta\!\left(\frac{d\tau}{dt}\right) \;=\; -\frac{2 GM \, e}{a c^2} \cos E + \mathcal{O}(e^2).$$

Integrating in time using `dt = (1 − e cos E)/n · dE`:

$$\Delta\tau(t) - \Delta t \;=\; \int_0^t \Delta\!\left(\frac{d\tau}{dt'}\right) dt' \;=\; -\frac{2 GM \, e}{a c^2 n} \int_0^E \cos E'\,(1 - e\cos E')\, dE'.$$

To first order in `e`, the integral is `sin E`. Using `n² = GM/a³`, so that `GM/(a²n) = GM/(a²) · √(a³/GM) = √(GM/a) = √(GMa)/a`, the prefactor becomes:

$$\frac{2 GM \, e}{a c^2 n} \;=\; \frac{2 e}{c^2} \sqrt{GM \, a}.$$

Hence the periodic eccentricity correction is:

$$\boxed{\;\Delta\tau_{\rm ecc}(t) \;=\; -\,\frac{2\sqrt{GM_\oplus \, a}}{c^2}\, e \, \sin E(t)\;}$$

This is identical to Ashby (2003) Eq. (43). The minus sign means: when the satellite is *near perigee* (small `r`, large `v`, both effects making the clock tick slow), `sin E > 0` and `Δτ < 0` (clock runs slow); when near apogee, the opposite. The correction is bounded by `±(2/c²)√(GMa) · e`.

## 4. Numerical evaluation

The amplitude coefficient:

$$\frac{2\sqrt{GM_\oplus \, a_{\rm GPS}}}{c^2} \;=\; \frac{2 \sqrt{3.986 \times 10^{14} \times 2.656 \times 10^7}}{(2.998 \times 10^8)^2} \;\approx\; 2.290 \times 10^{-6}\ {\rm s}.$$

Peak amplitudes for the GPS eccentricity range:

| `e` (typical) | Peak `|Δτ_ecc|` (ns) |
|---|---|
| `0.005` (low) | `11.45 ns` |
| `0.01` (typical) | `22.90 ns` |
| `0.02` (max spec) | `45.79 ns` |

The peak occurs at `E = π/2` (quarter-orbit from perigee), with period equal to the orbital period (`11.97 h`).

A peak of `46 ns` corresponds to a pseudorange error of `c × 46 × 10⁻⁹ ≈ 14 m` if uncorrected — comparable to the other dominant GPS error budget items (ionospheric, tropospheric). The eccentricity correction is therefore *operationally* essential, even though it is `~10⁻³` smaller than the headline `38 μs/day`.

## 5. Wolfram MCP check

```wolfram
cc=299792458; GMe=3.986004418*^14; aGPS=26560000; coeff = 2 Sqrt[GMe aGPS]/cc^2; Print["amplitude coefficient (s) = ", ScientificForm[N[coeff,8]]]; Print["e=0.005 peak (ns) = ", N[coeff 0.005 1*^9,6]]; Print["e=0.01 peak (ns) = ", N[coeff 0.01 1*^9,6]]; Print["e=0.02 peak (ns) = ", N[coeff 0.02 1*^9,6]]; Print["orbital period (h) = ", N[2 Pi Sqrt[aGPS^3/GMe]/3600,6]]
```

**Result:**
```
amplitude coefficient (s) = 2.2896629×10⁻⁶
e=0.005 peak (ns) = 11.4483
e=0.01 peak (ns) = 22.8966
e=0.02 peak (ns) = 45.7933
orbital period (h) = 11.9660
```
✅ matches §4.

## 6. Comparison with Ashby (2003)

Ashby §4 develops the eccentricity correction starting from the same vis-viva-substituted master equation. His Eq. (43) gives the boxed result above with identical sign and prefactor. The GPS Interface Control Document (ICD-200) cites Ashby's Eq. (43) directly as the receiver-implemented correction; the formula is built into every consumer GPS chip on Earth. The implementation uses the broadcast `e` and the eccentric anomaly `E` computed from the broadcast mean anomaly via Kepler's equation.

A subtle point: the eccentricity correction *exists* even though the constant frequency offset (effect 04) has been applied. The frequency offset is calibrated to the *time-averaged* satellite-orbit rate, which equals the circular-orbit rate at `r = a` (the result derived in §3 above). The periodic departure from that mean is what `e sin E` parametrises.

## 7. Verdict

✅ reproduces Ashby (2003) §4 Eq. (43) at quoted precision. The peak `≤46 ns` figure for `e = 0.02` matches the operational ICD-200 cited maximum. This correction is implemented in every GPS receiver.
