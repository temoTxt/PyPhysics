# Effect pt_09 — proper-time companion to full clock equation (synthesis + open questions)

**One-line summary:** The proper-time framework's predictions for the eight GPS effects match the standard-method predictions at GPS-relevant precision (relative agreement `~10⁻¹⁰` or better for every effect). The framework's full clock equation has the same operational structure — pre-launch frequency offset `Δf/f = −4.4647 × 10⁻¹⁰`, periodic eccentricity broadcast in ephemeris, Sagnac and Shapiro applied by the receiver. The framework adds *no operationally visible terms* for GPS; the campaign records four open questions where a derived (not speculative) GR-extended framework would clarify the structural claims behind the agreement.

**Status:** (synthesis).
**Standard-method companion:** [`09_full_clock_equation.md`](09_full_clock_equation.md).

<!-- TODO: human reviews and fills in — confirms that the synthesis honestly summarises the per-effect results, and that the four open questions in §3 capture the right load-bearing items for a future GR-extension paper. -->

## 1. Summary of per-effect verdicts

The eight per-effect documents reach the following verdicts:

| Effect | Standard prediction | PT prediction | Verdict |
|---|---|---|---|
| pt_01 metric setup | `dτ/dt = 1 − GM/(rc²) − v²/(2c²)` | `dτ/dt = 1 − GM/(rb²) − u²/(2b²)` | ⚠ structural difference; reduces to standard at GPS |
| pt_02 gravitational dilation | `+45.65 μs/day` | matches standard (relative `~10⁻²⁰`) | ❌ framework GR-extended speculatively |
| pt_03 velocity dilation | `−7.11 μs/day` | `−7.11 μs/day` *exactly* by algebraic identity | ✅ exact framework equivalence |
| pt_04 net headline | `+38.57 μs/day` | `+38.57 μs/day` (matches to `~10⁻¹⁸`) | ⚠ partial; combined PT = standard at GPS |
| pt_05 eccentricity | `±46 ns` periodic | matches; third term `~10⁻³⁴ s/signal` | ⚠ framework adds invisible third term |
| pt_06 Sagnac | `≤ 137 ns` | `≤ 137 ns` (matches to `~10⁻¹⁰ ns`) | ✅ direct framework application |
| pt_07 Shapiro | `≤ 62 ps` | identical to standard (vacuum photon: `b = c`) | ❌ framework GR-extended; vacuum case trivial |
| pt_08 J₂ oblateness | `+0.033 μs/day` | matches standard (relative `~10⁻¹²`) | ❌ framework GR-extended speculatively |

**Operational headline:** the framework predicts the same pre-launch frequency offset `Δf/f = −4.4647 × 10⁻¹⁰` and the same `~38.6 μs/day` net rate as standard. There is *no observable difference* between the two formulations at GPS-relevant precision.

## 2. The full PT clock equation

The PT operational equation, expressed in the same form as standard [`09_full_clock_equation.md`](09_full_clock_equation.md) §3:

$$\boxed{
\begin{aligned}
\tau_{\rm sat}(t) \;=\;& (1 + \delta_{\rm offset, PT})\, t \\
                 & \;-\; \frac{2 \sqrt{GM_\oplus a}}{b_{\rm sat}^2}\, e \sin E(t) \\
                 & \;+\; \delta_{\rm clock}(t) \;+\; \delta_{\rm 3rd\text{-}term} \;+\; \delta_{\rm relativistic\ extras}.
\end{aligned}
}$$

with

$$\delta_{\rm offset,\,PT} \;=\; \frac{GM_\oplus}{R_\oplus b_{\rm grnd}^2} - \frac{3 GM_\oplus}{2 a_{\rm GPS} b_{\rm sat}^2} + \frac{\omega_\oplus^2 R_\oplus^2}{2 b_{\rm grnd}^2} \;+\; \delta_{J_2,\,{\rm PT}}.$$

The two innovations relative to standard:
- All `1/c²` factors become `1/b²` (different `b` for satellite vs receiver). At GPS, the difference is invisible.
- `δ_3rd-term` is the framework's third-term contribution to the eccentricity correction. For circular orbits it vanishes; for `e = 0.02` GPS orbits, `~10⁻³⁴ s` per signal — invisible.

The signal-arrival-time correction:

$$\boxed{
\begin{aligned}
t_{\rm arrival,\,ECEF} \;=\;& \frac{L}{c} \\
                       & \;-\; \frac{2 \boldsymbol{\omega}_\oplus \cdot \mathbf{A}}{b_{\rm recv}^2} \\
                       & \;+\; \frac{2 GM_\oplus}{b_{\rm photon}^3} \ln\!\left(\frac{r_T + r_R + L}{r_T + r_R - L}\right).
\end{aligned}
}$$

For vacuum photons, `b_photon = c` exactly (pt_07 §3), so the Shapiro term is identical to standard. The Sagnac term carries `b_recv ≈ c` for an Earth-rotating receiver — invisible correction.

## 3. Four open questions for a future GR-extended Gill framework

The campaign's PT predictions are conditional on speculative extensions in four places. A derived (not speculative) GR-extended Gill framework would resolve these.

### Q1. What is the correct curved-spacetime extension of the framework?

The minimal-extension hypothesis used in pt_01, pt_02, pt_04, pt_07, pt_08 is `c → b` in the Schwarzschild metric. This is *one* of several possible extensions; alternatives include (a) `c → b` only in the time-time component, (b) treating null vs timelike geodesics differently, (c) replacing the metric structure entirely with a `b`-dependent connection. The minimal extension reduces to standard GR at low `u/c`; non-minimal extensions might predict observable corrections.

For GPS, any reasonable extension reduces to the same `+38.6 μs/day` at the operational precision floor. The framework's GR extension cannot be tested by GPS — it would require either much-larger-velocity clock comparisons, or much-stronger-gravity environments, or both.

### Q2. How does `u · a` transform under the proper-time boost?

This question is flagged in [`_proper_time_cheatsheet.md` §4a](../Electromagnetism/_proper_time_cheatsheet.md) for the radiation-reaction sub-investigation. It also affects pt_05's third-term contribution: whether `u · a` is frame-independent matters for predicting the third-term's observable size in non-inertial frames. For GPS circular orbits this is moot (`u · a = 0`); for elliptical GPS orbits the third term is `~10⁻³⁴ s/signal` regardless of how `u · a` transforms.

### Q3. Does the framework predict a Sagnac-like contribution from the satellite's own motion?

In pt_06, the standard Sagnac formula carries `b_recv` (the receiver's local light speed). A symmetric reading of the framework might suggest that the *satellite*'s `b_sat` should also appear — perhaps as a satellite-side Sagnac-like correction. For GPS this is `~10⁻¹⁰` of the standard Sagnac, at the edge of next-generation optical-clock receivers but invisible to current GPS.

### Q4. Is there a framework-specific J₂-like correction beyond the speculative extension?

pt_08 records the J₂ correction under the speculative extension. A first-principles derivation of multipole-induced clock corrections in the framework would clarify whether the standard `+0.033 μs/day` is the leading-order term or whether additional `b`-dependent contributions appear at higher order.

## 4. Honest summary

The proper-time framework reproduces standard GPS predictions at all operationally relevant precisions. This is the expected outcome of the framework's claim of mathematical equivalence to standard SR + flat-space Maxwell in the appropriate limit, and is consistent with — but does not independently validate — the framework. **The campaign found no observable disagreement between PT and standard at GPS precision.** A non-zero finding (PT predicting something different at GPS precision that disagrees with measurement) would have been a serious problem for the framework; none was found.

The campaign's value to the research program is twofold:

1. **A consolidated point-for-point reference** for what the framework predicts about GPS, with explicit per-effect derivations and honest scoping of where the framework requires (speculative) extension.
2. **A set of four open questions** flagged for a future GR-extended Gill framework. Each is identified by which per-effect document surfaces it; each is operationally invisible at GPS precision but structurally load-bearing for the framework's eventual extension to curved spacetime.

## 5. Suggested next work

- **Resolve the GR extension question** in a future Gill–Zachary paper. The campaign's speculative `c → b` hypothesis is the *minimal* extension; deriving the right extension from first principles would close pt_01's `<!-- TODO -->` block.
- **Pick up the `u · a` covariance question** in the Electromagnetism sub-investigation ([Ch14_P14.2](../Electromagnetism/Jackson/Ch14_Radiation_by_Moving_Charges.md), [Ch16](../Electromagnetism/Jackson/Ch16_Radiation_Damping.md), [issue #43](https://github.com/temoTxt/PyPhysics/issues/43)). GPS is not a sensitive probe of this; high-acceleration experiments are.
- **Optical-clock GPS** (next-generation, `~10⁻¹⁷` clock stability) is the GPS-side regime where the framework's `~10⁻¹⁰` PT-vs-standard correction could *just* start to matter. Re-running this campaign at optical-clock precision would be useful when the relevant satellite missions come online.

## 6. Verdict

✅ **The proper-time framework reproduces the full standard GPS relativistic clock correction at all operationally relevant precisions.** The framework's predictions are conditional on a speculative GR extension (where applicable); under the minimal-extension hypothesis, all eight effects match standard. The framework adds *no operationally visible terms* for current-generation GPS atomic clocks. Four open questions for a future GR-extended Gill framework are flagged in §3.
