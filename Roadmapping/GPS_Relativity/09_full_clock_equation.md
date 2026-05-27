# Effect 09 — Full clock equation (synthesis of all corrections)

**One-line summary:** The operational GPS timing equation, combining all eight effects from the campaign into a single expression that maps satellite-clock proper time `τ_sat` to ground-clock proper time `τ_grnd` as a function of orbital phase and signal geometry. The headline secular rate `+38.57 μs/day` is built into the pre-launch frequency offset; the remaining periodic, path-dependent, and geometry-dependent terms are applied by the GPS ground segment (for orbit corrections) and the receiver (for Sagnac and Shapiro).

**Status:** ✅ drafted.
**Ashby (2003) cross-reference:** §8, summary equations.

<!-- TODO: human reviews and fills in — confirms the synthesis below correctly combines the eight per-effect documents and identifies the operational responsibility (pre-launch offset / ephemeris broadcast / receiver firmware) for each term -->

## 1. Effect statement

GPS uses a layered correction architecture. Some corrections are *secular* (constant in time) and can be absorbed into the pre-launch frequency offset; others are *periodic* (vary with orbital phase) and require ephemeris-broadcast parameters; others are *path-dependent* (vary with signal geometry) and must be computed by the receiver. This document collects the corrections of effects 02–08 into a single operational equation and labels each term by where it is applied in the GPS architecture.

The total transformation from coordinate time `t` (TCG) to a clock's proper time `τ` is the integral of `dτ/dt` over the clock's worldline plus the propagation-correction terms applied to inter-clock signals.

## 2. Setup

We follow the convention of Ashby (2003) §8: separate the corrections into

1. **Secular rate offset** (effects 02 + 03 + 08): fixed multiplicative rate, absorbed into the pre-launch frequency setting.
2. **Periodic eccentricity term** (effect 05): time-varying, broadcast in ephemeris, applied by receiver.
3. **Propagation corrections** (effects 06 + 07): signal-path-dependent, applied by receiver from broadcast satellite position.

## 3. The full clock equation

The proper time of a GPS satellite clock, expressed in terms of GPS system time, is:

$$\boxed{
\begin{aligned}
\tau_{\rm sat}(t) \;=\;& (1 + \delta_{\rm offset})\, t \\
                 & \;-\; \frac{2 \sqrt{GM_\oplus \, a}}{c^2}\, e \sin E(t) \\
                 & \;+\; \delta_{\rm clock}(t) \;+\; \delta_{\rm relativistic\ extras}.
\end{aligned}
}$$

where:

| Term | Magnitude | Effect | Source |
|---|---|---|---|
| `(1 + δ_offset) t` | `δ_offset = +4.4647 × 10⁻¹⁰` (secular) | 02 + 03 + 08 | pre-launch frequency offset (`10.229 999 995 43 MHz`) |
| `−(2√(GMa)/c²) e sin E` | `≤ ±46 ns` periodic (`e ≤ 0.02`) | 05 | computed by receiver from broadcast `e`, `E` |
| `δ_clock(t)` | `~ ns` random walk | atomic clock noise | corrected via ground-segment monitoring + broadcast clock corrections |
| `δ_relativistic extras` | `< 1 ns` | residuals | gravitational-wave / tidal / next-order PN |

The signal-arrival-time correction, applied by the receiver, is:

$$\boxed{
\begin{aligned}
t_{\rm arrival\,ECEF} \;=\;& \frac{L}{c} \\
                       & \;-\; \frac{2 \omega_\oplus A_E}{c^2} \\
                       & \;+\; \frac{2 GM_\oplus}{c^3} \ln\!\left(\frac{r_T + r_R + L}{r_T + r_R - L}\right).
\end{aligned}
}$$

| Term | Magnitude | Effect | Source |
|---|---|---|---|
| `L/c` | `~ 0.07 s` (geometric) | classical light-time | direct geometric computation |
| `−(2ω⊕/c²) A_E` | `≤ ±137 ns` | 06 Sagnac | receiver-computed from satellite and receiver ECEF positions |
| `(2GM/c³) ln(...)` | `≤ ~62 ps` | 07 Shapiro | typically omitted (below noise floor); included in geodetic post-processing |

## 4. Numerical budget summary

The full GPS relativistic error budget, in canonical units:

| Effect | Magnitude | Type | Where applied |
|---|---|---|---|
| 02 Gravitational time dilation | `+45.65 μs/day` | secular | pre-launch offset |
| 03 Velocity time dilation | `−7.11 μs/day` | secular | pre-launch offset |
| 04 Net headline (02 + 03) | `+38.54 μs/day` (pointmass) | secular | pre-launch offset |
| 08a J₂ correction to geoid | `+0.033 μs/day` | secular | pre-launch offset |
| 08b J₂ correction to satellite orbit | `−0.001 μs/day` (avg) | secular | pre-launch offset |
| **04 + 08 operational net** | **`+38.57 μs/day`** | **secular** | **pre-launch (`Δf/f = −4.4647 × 10⁻¹⁰`)** |
| 05 Eccentricity periodic | `≤ ±46 ns` (`e=0.02`) | periodic | broadcast ephemeris + receiver |
| 06 Sagnac per signal | `≤ ±137 ns` | path-dependent | receiver firmware |
| 07 Shapiro per signal | `≤ +62 ps` | path-dependent | geodetic post-processing only |
| 08c Lense-Thirring | `~ 10⁻¹⁶` | sub-ps | (not currently corrected) |
| 08d Higher PN | `~ 10⁻¹⁹` | below noise floor | (not currently corrected) |

If any of effects 02 + 03 + 04 + 05 + 06 were ignored, the resulting pseudorange error per day:

| Ignored effect | Pseudorange drift |
|---|---|
| 02 + 03 + 04 (headline) | `c × 38.57 μs ≈ 11.56 km/day` |
| 05 (eccentricity max) | `c × 46 ns ≈ 14 m` periodic |
| 06 (Sagnac max) | `c × 137 ns ≈ 41 m` per signal |
| 07 (Shapiro max) | `c × 62 ps ≈ 1.87 cm` per signal |
| 08 (J₂ corrections) | `c × 0.033 μs ≈ 10 m/day` |

The popular `~10 km/day` figure cited in the Space Daily article ([issue #57](https://github.com/temoTxt/PyPhysics/issues/57)) is the dominant headline `38 μs/day` × `c`, exactly the `11.6 km/day` we compute.

## 5. Wolfram MCP check

```wolfram
cc=299792458; secOffsetDay = 38.57*^-6; Print["headline km/day drift = ", N[cc secOffsetDay/1000, 6]]; ecc05Drift = cc 46*^-9; Print["eccentricity m drift = ", N[ecc05Drift, 6]]; sagDrift = cc 137*^-9; Print["Sagnac m drift = ", N[sagDrift, 6]]; shaDrift = cc 62*^-12; Print["Shapiro m drift = ", N[shaDrift, 6]]; J2Drift = cc 0.033*^-6; Print["J2 m/day drift = ", N[J2Drift, 6]]
```

**Result:**
```
headline km/day drift = 11.5630
eccentricity m drift = 13.7904
Sagnac m drift = 41.0716
Shapiro m drift = 0.0185872
J2 m/day drift = 9.89315
```
✅ matches §4.

## 6. Comparison with Ashby (2003) — the operational picture

Ashby §8 lays out the full operational equation as a closed-form expression. The principal difference between this campaign's synthesis and Ashby's is purely organizational: we have separated the corrections into nine per-effect documents, each with its own derivation and Wolfram MCP check, whereas Ashby gives the synthesis directly. The agreement on every quoted number is exact (within rounding precision of the input parameters); the disagreement is only in the bookkeeping of which decomposition is used for which contribution.

The campaign reproduces:
- the headline `+38.57 μs/day` (Ashby Eq. (39)) ✅
- the pre-launch frequency offset `−4.4647 × 10⁻¹⁰` (Ashby Eq. (39)) ✅
- the periodic eccentricity amplitude `±(2/c²)√(GMa) · e` (Ashby Eq. (43)) ✅
- the Sagnac correction `(2ω⊕/c²) · A_E` (Ashby Eq. (35)) ✅
- the Shapiro delay `(2GM/c³) ln(...)` (Ashby Eq. (53)) ✅
- the J₂ correction to the geoid potential (Ashby Eq. (24)) ✅

with no structural disagreement.

## 7. Verdict

✅ **The standard-method derivation of all eight effects reproduces Ashby (2003) §8's operational clock equation at the precision quoted by Ashby for every effect.** The campaign's value is not the production of new numbers; it is the consolidation of the full derivation in one place, with each effect numerically reproduced via Wolfram MCP from named parameter values.

## Next: the proper-time companion campaign

The companion campaign will repeat each of effects 01–08 under the Gill–Zachary proper-time substitution rules `c → b = √(c² + u²)`, `t → τ`, `∂_t → (b/c) ∂_τ`. Three structural questions to confront in that campaign:

1. **Effect 01 master equation.** Does the proper-time master equation reduce to `1 − GM/(rc²) − u²/(2b²)` (with `b` not `c`), and if so, what does that change?
2. **Effects 02–04 headline.** Does the proper-time formulation reproduce `+38.57 μs/day`, and if not, where does the deviation appear? At GPS velocities `u ≪ c`, `b ≈ c` to 11 figures, so the proper-time prediction should agree at the GPS-relevant precision; any structural deviation would be a flag for the proper-time framework's interpretation of clock rates.
3. **Effect 05 eccentricity.** The eccentricity correction's `2 sin E / n` integrating-factor structure assumes a specific orbital-motion frame. Does this survive the substitution to proper-time orbital motion?

These questions are deferred to the proper-time companion campaign and are not addressed in this PR.
