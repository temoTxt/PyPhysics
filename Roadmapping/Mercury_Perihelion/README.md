# Mercury perihelion — dual-Newtonian framework check

This campaign verifies the math of Gill's [*Relativistic Newtonian Theory*](../Tepper_Gill_Papers/Dual%20Newtonian%20Theory.tex) (user-supplied filename: *Dual Newtonian Theory*) and computes the Mercury perihelion advance under the three framework predictions the paper offers — Corda, full dual, approximate dual — comparing each against standard GR and the observed value.

The campaign is governed by the plan at [`.dev/tasks/81_Perhelion_of_Mercury.md`](../../.dev/tasks/81_Perhelion_of_Mercury.md) (issue [#81](https://github.com/temoTxt/PyPhysics/issues/81)).

<!-- TODO: human reviews and fills in — confirms the framing of the campaign as a verification of the paper's *gravitational* extension of the dual framework, distinct from but cross-referenced to the GPS author report's Q1 (curved-spacetime extension of the framework). -->

## Headline result

The paper's algebra reproduces ✅ via Wolfram MCP; the numerical perihelion advance is computed below.

| Prediction | Δφ per century | vs observed (43″/century) | vs GR (42.99″/century) |
|---|---|---|---|
| Corda `Δφ_c = πm/M` | **44.66″** | +3.9% | +3.9% |
| Full dual `Δφ_d` | **37.79″** | **−12.1%** | **−12.1%** |
| Approximate dual `Δφ_{d₁}` | **37.79″** | **−12.1%** | **−12.1%** |
| Standard GR `6πGM/[c²a(1−e²)]` | 42.99″ | +0.0% (reference) | (reference) |
| **Observed** (Le Verrier residual) | **≈ 43″** | (reference) | +0.02% |

**Headline finding.** The paper highlights Corda's `44.39″/century` as the framework's match to observation, but Corda's value is the *reduced-mass-only* Newtonian effect `πm/M`, **not** the framework's structural gravitational prediction. The framework's own full dual prediction is 37.79″/century — about 12% below observed and with the *opposite sign* of relativistic correction compared to GR.

Decomposition of the framework's full dual:

- Reduced-mass (Corda) contribution: `+44.66″/century`
- Relativistic contribution from the `(1 − GM/(c²r))` factor in Eq. (h4): **`−6.87″/century`** (note: *negative* — the framework's modified gravity is *weaker* at small `r`, producing perihelion *recession*, not advance)
- Net: `+44.66 − 6.87 = +37.79″/century`

This is a substantive finding for the GPS author report's open question Q1; see [`04_findings_and_GPS_Q1.md`](04_findings_and_GPS_Q1.md).

## Status

| # | Document | Content | Status |
|---|---|---|---|
| 01 | [`01_setup_and_verification.md`](01_setup_and_verification.md) | Paper's `K = π²/(2m) + mc² + V + V²/(2mc²)` → Eq. (h3) → Eq. (h4) chain | ✅ all Wolfram-MCP verified |
| 02 | [`02_orbital_dynamics.md`](02_orbital_dynamics.md) | Eqs. (h5)–(h8) (orbital speed, period, angular velocity) | ✅ all Wolfram-MCP verified |
| 03 | [`03_numerical_predictions.md`](03_numerical_predictions.md) | Numerical Δφ for Mercury under all three predictions + GR + observed | ✅ computed |
| 04 | [`04_findings_and_GPS_Q1.md`](04_findings_and_GPS_Q1.md) | Honest scoping + cross-reference to GPS Q1 | ✅ flagged |

Verification doc with per-equation entries: [`../Equation_Verification/Dual_Newtonian_Theory.md`](../Equation_Verification/Dual_Newtonian_Theory.md).

## Conventions

All Wolfram MCP blocks follow the `ClearAll` discipline established in the GPS campaign's effect 08 fix (see [PR #71 commit `210750d`](https://github.com/temoTxt/PyPhysics/pull/71/commits/210750d)) to prevent stale-variable contamination. Symbol conventions per CLAUDE.md: `ee` for charge, `potV` for potential, `GG` / `MM` / `mm` for `G` / `M` / `m` (avoid built-in collisions).

### Mercury orbital parameters (constants throughout)

| Quantity | Value | Source |
|---|---|---|
| `M_⊙` (solar mass) | `1.98892 × 10³⁰ kg` | IAU 2015 |
| `m_Mercury` | `3.3011 × 10²³ kg` | NASA fact sheet |
| `a_Mercury` (semi-major axis) | `5.7909 × 10¹⁰ m` | NASA fact sheet |
| `e_Mercury` (eccentricity) | `0.20563` | NASA fact sheet |
| `T_orb_Mercury` (sidereal period) | `87.969 days` | NASA fact sheet |
| `G` | `6.6743 × 10⁻¹¹ m³ kg⁻¹ s⁻²` | CODATA 2018 |
| `c` | `299 792 458 m/s` | SI definition |
| `m/M` (Mercury/Sun) | `1.6597 × 10⁻⁷` | derived |
| Orbits per tropical century | `415.20` | `100·365.25·86400 / T_orb` |
| `GM_⊙/(c² a_Mercury)` | `2.55 × 10⁻⁸` | derived |

## Source material

- **Gill, T. L. — *Relativistic Newtonian Theory*** (Howard Univ., undated; user-supplied .tex). Filename in issue #81: `Dual Newtonian Theory.pdf`. Paper title per `\title{}` field: *Relativistic Newtonian Theory*. Source: [`Roadmapping/Tepper_Gill_Papers/Dual Newtonian Theory.tex`](../Tepper_Gill_Papers/Dual%20Newtonian%20Theory.tex).
- **Corda, C.** — *Secret of perihelion of Mercury* (PDF attached to issue #81; not yet committed; will be ingested under `Roadmapping/Historical_Papers/Retrospective/` per the History campaign's PDF-acquisition policy).
- **Van Flandern, T. (1976)** — *A determination of the rate of change of G* and related perihelion analysis (PDF attached to issue #81; same disposition).
- **Standard GR reference:** Misner-Thorne-Wheeler *Gravitation* §40.5 for the Schwarzschild geodesic-derived `Δφ_GR = 6πGM/[c²a(1−e²)]`.

## Related work

- [GPS author report §5 Q1](../Author_Reports/2026-05_gps_relativity_summary_for_gill.md) — the open question this campaign provides a (negative) answer to.
- [Equation_Verification/Dual_Newtonian_Theory.md](../Equation_Verification/Dual_Newtonian_Theory.md) — per-equation verification entries.
- [GPS_Relativity/_proper_time_GPS_cheatsheet.md §3](../GPS_Relativity/_proper_time_GPS_cheatsheet.md) — applicability matrix for the GPS campaign; the Mercury campaign's finding informs the matrix's "❌ Framework as-published doesn't extend" entries.

## Honest limits

This campaign reproduces standard Mercury perihelion math (a 19th-century calculation refined by Einstein 1915) and compares it against the framework's predictions; the value to the repository is the explicit Wolfram-MCP verification of each step and the per-prediction breakdown. **The negative finding** — that the framework's full dual prediction under-predicts the observed Mercury perihelion advance by 12% and has the wrong sign of relativistic correction — is reported as a structural feature of the framework's modified-gravity proposal in this paper, *not* as a critique of the framework as a whole. The framework may have other proposals for the gravitational extension (Q1 lists four hypotheses considered in the GPS campaign); this paper tests *one* such extension, and the result is informative regardless of direction.
