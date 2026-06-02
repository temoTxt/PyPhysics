# GPS Relativity — proper-time companion campaign

This thread is the **proper-time companion** to the standard-method GPS relativity campaign. It applies the [Gill–Zachary proper-time substitution rules](_proper_time_GPS_cheatsheet.md) — verified at [`Roadmapping/Equation_Verification/Two_Mathematically_Equivalent_Versions_of_Maxwells_Equations.md`](../Equation_Verification/Two_Mathematically_Equivalent_Versions_of_Maxwells_Equations.md) — to each of the nine effects in the standard derivation, recording for each effect whether the framework applies, what it predicts, and where (if anywhere) it would diverge from standard.

The thread is governed by the plan at [`.dev/tasks/57-gps-proper-time-companion.md`](../../.dev/tasks/57-gps-proper-time-companion.md) (issue [#57](https://github.com/temoTxt/PyPhysics/issues/57), sibling to PR [#71](https://github.com/temoTxt/PyPhysics/pull/71)).

<!-- TODO: human reviews and fills in — confirms the framing below honestly states the framework's applicability boundary at the SR/GR transition before the per-effect documents go into substantive detail -->

## Scope and the SR / GR boundary

The Gill–Zachary substitution rules govern (a) Maxwell's equations and (b) flat-space relativistic dynamics. They have **not** been extended to curved spacetime in any verified Gill paper. For GPS, this divides the nine effects into three categories:

| Category | Effects | Framework applies as-published? |
|---|---|---|
| **SR + flat-space EM** | 03 (velocity dilation), 06 (Sagnac) | ✅ yes — derivation runs through `b² = c² + u²` Maxwell rules |
| **Schwarzschild + GR** | 02 (gravitational redshift), 07 (Shapiro), 08 (J₂) | ❌ no — framework would need a curved-spacetime extension |
| **Combined / kinematic + gravitational** | 01 (metric setup), 04 (headline), 05 (eccentricity) | ⚠ partial — kinematic piece applies, gravitational piece needs extension |

For the categories ❌ and ⚠ the per-effect documents record a *speculative minimal extension* (`c → b` substituted into the Schwarzschild metric, with `u² → 0` reducing it to standard GR), clearly flagged with substantive-AI disclosure per [Crocco §1](../Tooling/CROCCO_COMPLIANCE.md). This is not a derivation of GR-extended proper-time gravity — it is a "what would the framework predict *if* it extended this way" exercise, included so that the side-by-side comparison with the standard derivation is complete.

The headline operational result is **unchanged from the standard derivation**: at GPS-relevant precision (`u² / c² ≈ 1.67 × 10⁻¹⁰`, gravitational potential `GM⊕/(rc²) ≈ 10⁻¹⁰`), the proper-time formulation reduces to standard SR + GR with no observable difference.

## Status

| # | Effect | Standard PR file | PT companion file | Framework applies? | PT predicts |
|---|---|---|---|---|---|
| pt_01 | PT Schwarzschild setup | [`01_...`](01_schwarzschild_metric_setup.md) | [`pt_01_...`](pt_01_proper_time_metric_setup.md) | ⚠ partial | `dτ/dt = 1 − GM/(rc²) − u²/(2b²)` (speculative extension) |
| pt_02 | PT gravitational time dilation | [`02_...`](02_gravitational_time_dilation.md) | [`pt_02_...`](pt_02_gravitational_time_dilation.md) | ❌ GR-only | matches standard (speculative extension reduces) |
| pt_03 | PT velocity time dilation | [`03_...`](03_velocity_time_dilation.md) | [`pt_03_...`](pt_03_velocity_time_dilation.md) | ✅ yes | `dτ/dt = c/b`, exactly equivalent to `1/γ` |
| pt_04 | PT net headline | [`04_...`](04_net_38us_per_day.md) | [`pt_04_...`](pt_04_net_38us_per_day.md) | ⚠ partial | `+38.57 μs/day` (same number; bookkeeping reorganised) |
| pt_05 | PT eccentricity correction | [`05_...`](05_eccentricity_correction.md) | [`pt_05_...`](pt_05_eccentricity_correction.md) | ⚠ partial | matches standard; third term `(u·a)/b⁴` contribution vanishes for circular orbit |
| pt_06 | PT Sagnac effect | [`06_...`](06_sagnac_effect.md) | [`pt_06_...`](pt_06_sagnac_effect.md) | ✅ yes | `Δt = −2ω · A / b²` — same algebra with `c → b`, matches standard at GPS precision |
| pt_07 | PT Shapiro delay | [`07_...`](07_shapiro_delay.md) | [`pt_07_...`](pt_07_shapiro_delay.md) | ❌ GR-only | matches standard (speculative extension reduces) |
| pt_08 | PT J₂ oblateness | [`08_...`](08_j2_oblateness.md) | [`pt_08_...`](pt_08_j2_oblateness.md) | ❌ GR-only | matches standard (speculative extension reduces) |
| pt_09 | PT synthesis | [`09_...`](09_full_clock_equation.md) | [`pt_09_...`](pt_09_synthesis_and_open_questions.md) | (synthesis) | full proper-time clock equation + open extension questions |

## Conventions

The PT per-effect template at [`_template_pt_effect.md`](_template_pt_effect.md) mirrors the standard template with two additions: §0 (Framework applicability) and §5b (GPS-precision-limit equivalence to standard).

- **Pragmatic / substantive labelling.** pt_03, pt_05, pt_06 are *pragmatic* — direct framework application. pt_01, pt_02, pt_04, pt_07, pt_08 are *substantive* (carry a speculative GR-extension) — each carries `<!-- TODO: human reviews and fills in -->` on the extension paragraph per [Crocco §1](../Tooling/CROCCO_COMPLIANCE.md).
- **Wolfram MCP checks** confirm the numerical equivalence to standard at the operational order. They use the same parameter values as the sibling PR — see [`README.md`](README.md) §Units and parameters.

## Source material

- **Gill, T. L. & Zachary, M. W. — proper-time substitution rules.** Verified at [`Roadmapping/Equation_Verification/Two_Mathematically_Equivalent_Versions_of_Maxwells_Equations.md`](../Equation_Verification/Two_Mathematically_Equivalent_Versions_of_Maxwells_Equations.md). Eqs. (1)–(23) verified ✅; Eq. (24) carries an unresolved finding (not invoked in this campaign).
- **Ashby, N. (2003).** *Relativity in the Global Positioning System*. Living Reviews in Relativity, 6(1). The standard-method reference. Bibliography stub: [`Roadmapping/History/Bibliography/Retrospective/ashby2003_relativity_in_gps.md`](../History/Bibliography/Retrospective/ashby2003_relativity_in_gps.md).
- **Proper-time cheatsheets.** [`_proper_time_GPS_cheatsheet.md`](_proper_time_GPS_cheatsheet.md) (this PR, the GPS-specific applicability matrix), and the EM-focused [`Roadmapping/Electromagnetism/_proper_time_cheatsheet.md`](../Electromagnetism/_proper_time_cheatsheet.md) (verified rules).

## Honest limits

The framework's prediction at GPS-relevant precision is the same as standard SR + GR by construction (mathematical equivalence at this order). A *negative* finding — the framework predicting something at GPS precision that disagrees with measurement — would be a serious problem for the framework. The campaign is structured to be sensitive to such a finding, not to be biased against it. None was found. This is consistent with the framework's claim of mathematical equivalence in the appropriate limit and does not constitute a novel validation of either formulation.

The author's prior involvement in the framework (co-author of *Dual Relativistic Quantum Mechanics I*, 2021) is acknowledged: this is exploratory work in the Gill–Zachary framework, not a neutral evaluation of it. Per-effect verdicts are conditional predictions of the form "*if* the framework extends to curved spacetime as the minimal-extension hypothesis suggests, *then* the prediction at GPS precision is X." The "if" is the load-bearing caveat.

## Related work

- [PR #71 (sibling)](https://github.com/temoTxt/PyPhysics/pull/71) — the standard-method derivation this PR companions.
- [`Roadmapping/Electromagnetism/`](../Electromagnetism/) — Jackson canonical problems × proper-time. Methodologically parallel campaign; the same SR + flat-space PT Maxwell substitution rules apply.
- [`Roadmapping/Equation_Verification/Two_Mathematically_Equivalent_Versions_of_Maxwells_Equations.md`](../Equation_Verification/Two_Mathematically_Equivalent_Versions_of_Maxwells_Equations.md) — the verified Gill–Zachary paper. The rules of this campaign come from there.
