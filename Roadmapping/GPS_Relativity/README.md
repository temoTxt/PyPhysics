# GPS Relativity — standard-method clock corrections

This research thread reproduces, in one consolidated derivation, the full relativistic clock correction applied to the GPS satellite constellation. The campaign covers nine separately-named effects, from the Schwarzschild line element through to the operational clock equation programmed into the GPS ground segment. Each per-effect document is paired with a Wolfram-MCP numerical sanity check and a comparison against the canonical [Ashby (2003) *Living Reviews* treatment](https://link.springer.com/article/10.12942/lrr-2003-1).

The thread is governed by the plan at [`.dev/tasks/57-gps-standard-derivation.md`](../../.dev/tasks/57-gps-standard-derivation.md) (issue [#57](https://github.com/temoTxt/PyPhysics/issues/57)).

<!-- TODO: human reviews and fills in — confirms the framing above accurately states the thread's scope and conditional posture before the proper-time companion campaign begins -->

## Scope

The "standard method" is the textbook Schwarzschild + post-Newtonian treatment expanded to order `1/c²`, as worked out by Ashby (1975, 2003), Allan & Ashby (1985), and Hodge (1992). It is *not* the proper-time formulation of the Gill–Zachary substitution rules — that companion campaign is deferred to a future PR.

The campaign is **not** a critique of the standard answer; it is a careful consolidation of it. Where this document reproduces a result already in Ashby (2003), the verdict line names the section + equation in Ashby for direct cross-reference.

## Status

| # | Effect | File | Magnitude (per day) | Status |
|---|---|---|---|---|
| 01 | Schwarzschild metric setup | [`01_schwarzschild_metric_setup.md`](01_schwarzschild_metric_setup.md) | (foundation) | ✅ drafted |
| 02 | Gravitational time dilation | [`02_gravitational_time_dilation.md`](02_gravitational_time_dilation.md) | **+45.7 μs** | ✅ drafted |
| 03 | Velocity time dilation | [`03_velocity_time_dilation.md`](03_velocity_time_dilation.md) | **−7.1 μs** | ✅ drafted |
| 04 | Net headline ±38 μs/day | [`04_net_38us_per_day.md`](04_net_38us_per_day.md) | **+38.5 μs** | ✅ drafted |
| 05 | Eccentricity correction | [`05_eccentricity_correction.md`](05_eccentricity_correction.md) | ≤ 46 ns (periodic, e=0.02) | ✅ drafted |
| 06 | Sagnac effect | [`06_sagnac_effect.md`](06_sagnac_effect.md) | ≤ 137 ns (path-dependent) | ✅ drafted |
| 07 | Shapiro delay | [`07_shapiro_delay.md`](07_shapiro_delay.md) | ≤ 62 ps (path-dependent) | ✅ drafted |
| 08 | J₂ oblateness + higher PN | [`08_J2_oblateness.md`](08_J2_oblateness.md) | ≈ sub-ps secular | ✅ drafted |
| 09 | Full clock equation | [`09_full_clock_equation.md`](09_full_clock_equation.md) | (synthesis) | ✅ drafted |

The headline result reproduces Ashby (2003) Eq. (39): the satellite-clock frequency is offset by `−4.4647 × 10⁻¹⁰` *before launch* so that, once in orbit, it ticks at the same rate as a clock on the geoid.

## Conventions

Three meta-documents govern every per-effect write-up:

- [`_template_effect.md`](_template_effect.md) — the per-effect document template. Every `0N_*.md` mirrors this structure.
- [`../Tooling/CROCCO_COMPLIANCE.md`](../Tooling/CROCCO_COMPLIANCE.md) — AI-use disclosure. The derivations themselves are *pragmatic* per [§5](../Tooling/CROCCO_COMPLIANCE.md); framing paragraphs in the README and §9 synthesis carry the standard substantive `<!-- TODO -->` block.
- [`../Tooling/VOICE_MATCH_GILL.md`](../Tooling/VOICE_MATCH_GILL.md) — voice guide. Applies to interpretive prose only; derivation steps are operational.

### Units and parameters

All derivations use SI units. Parameter values are quoted from the IERS Conventions (2010) and the WGS-84 geodetic standard, taken as exact for the purpose of numerical reproduction:

| Symbol | Value | Source |
|---|---|---|
| `c` | `299 792 458 m/s` (exact) | SI definition |
| `GM⊕` | `3.986 004 418 × 10¹⁴ m³/s²` | IERS / WGS-84 (defined) |
| `R⊕` (equatorial) | `6 378 137 m` | WGS-84 (defined) |
| `ω⊕` (sidereal rotation) | `7.292 115 1467 × 10⁻⁵ rad/s` | IERS 2010 |
| `J₂` | `1.082 626 7 × 10⁻³` | EGM2008 |
| `a_GPS` (semi-major axis) | `26 560 000 m` (nominal) | GPS ICD-200 |
| `v_GPS` (circular-orbit speed) | `3 874 m/s` (≡ √(GM⊕/a)) | derived |
| `Φ_geoid/c²` (geoid potential, relative) | `−6.969 290 134 × 10⁻¹⁰` | IAU 2000 / W₀ |

The GPS L-band carrier is generated from a fundamental clock frequency of `10.23 MHz`. The pre-launch frequency offset programmed into the satellite is `Δf/f = −4.4647 × 10⁻¹⁰`, giving an as-launched fundamental of `10.229 999 995 43 MHz`.

## Source material

- **Ashby, N. (2003). *Relativity in the Global Positioning System*.** Living Reviews in Relativity, 6(1). Bibliography stub: [`Roadmapping/History/Bibliography/Retrospective/ashby2003_relativity_in_gps.md`](../History/Bibliography/Retrospective/ashby2003_relativity_in_gps.md). The canonical reference; each per-effect verdict cites the relevant section.
- **Allan, D. W., Ashby, N., & Hodge, C. C. (1997).** *The Science of Timekeeping* (HP Application Note 1289). The operational engineering perspective.
- **Misner, Thorne, Wheeler, *Gravitation* (1973) §§25, 40.** Foundational treatment of the Schwarzschild geodesic equation and weak-field expansion.
- **Weinberg, S. (1972) *Gravitation and Cosmology* §9.** The post-Newtonian expansion in the form used throughout this campaign.

## Related work

- [`Roadmapping/Equation_Verification/`](../Equation_Verification/) — verification thread for the Gill papers. The proper-time companion campaign will depend on the verified Maxwell-paper substitution rules.
- [`Roadmapping/Electromagnetism/`](../Electromagnetism/) — Jackson canonical problems × proper-time. Methodologically parallel: textbook standard answer, then proper-time reformulation as a companion.

## Honest limits

This campaign reproduces a result that has been worked out by many authors and is operationally embedded in a deployed satellite constellation. Its value to the repository is consolidation in one place, with each effect numerically reproduced via Wolfram MCP from named parameter values. The campaign does **not** make a novel physics claim; it does not validate or critique the dual-theory framework. The companion proper-time campaign — to be opened as a sibling PR — will revisit each effect and record any structural difference.
