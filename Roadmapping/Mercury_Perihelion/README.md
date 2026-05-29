# Mercury perihelion — dual-Newtonian framework check

This campaign verifies the math of Gill's [*Relativistic Newtonian Theory*](../Tepper_Gill_Papers/Dual%20Newtonian%20Theory.tex) (user-supplied filename: *Dual Newtonian Theory*) and computes the Mercury perihelion advance under the three framework predictions the paper offers — Corda, full dual, approximate dual — comparing each against standard GR and the observed value.

The campaign is governed by the plan at [`.dev/tasks/81_Perhelion_of_Mercury.md`](../../.dev/tasks/81_Perhelion_of_Mercury.md) (issue [#81](https://github.com/temoTxt/PyPhysics/issues/81)).

<!-- TODO: human reviews and fills in — confirms the framing of the campaign as a verification of the paper's *gravitational* extension of the dual framework, distinct from but cross-referenced to the GPS author report's Q1 (curved-spacetime extension of the framework). -->

## Headline result

The paper's algebra reproduces ✅ via Wolfram MCP; the numerical perihelion advance is computed below.

| Prediction | Δφ per century | vs observed (43″/century) | vs GR (42.99″/century) |
|---|---|---|---|
| Corda `Δφ_c = πm/M` (reduced-mass; *not* a real precession) | 44.66″ | +3.9% | +3.9% |
| Paper full dual `Δφ_d` (circular-orbit heuristic) | 37.79″ | −12.1% | −12.1% |
| **Framework, proper eccentric-orbit calculation** ([doc 05](05_nbody_orbital_mechanics.md)) | **−7.17″** | **wrong sign** | **−1/6 exactly** |
| Standard GR `6πGM/[c²a(1−e²)]` | 42.99″ | +0.0% (reference) | (reference) |
| **Observed** (modern residual; Newcomb–Clemence reduction) | **≈ 43″** | (reference) | +0.02% |

The **proper orbital-mechanics row** is the load-bearing one: the paper's `+37.79″` comes from a circular-orbit heuristic that conflates orbital-frequency change with ellipse precession. The genuine perihelion advance of the framework's force law, computed via the apsidal angle for an eccentric orbit, is `−7.17″/century` — exactly `−1/6` of GR. See [doc 05](05_nbody_orbital_mechanics.md).

**Headline finding (sharpened in [doc 05](05_nbody_orbital_mechanics.md)).** A *proper orbital-mechanics* calculation (eccentric orbit, apsidal angle — not the paper's circular-orbit heuristic) shows the framework's force law `a = −(GM/r²)(1 − GM/(c²r))·ê_r` gives a genuine perihelion advance of **exactly `−1/6` of GR** — opposite sign, one-sixth magnitude: `−7.17″/century` vs GR's `+42.99″`. This is the classic "1/6 factor" of a force-law-only modification that lacks GR's spatial-metric curvature. **The N-body refinement cannot close this**: planetary perturbations contribute framework relativistic corrections `~10⁻⁴–10⁻⁵` of the Sun-Mercury term (the relativistic factor `GM/(c²r)` scales with the attracting mass, and the planets are light). The earlier "12% gap" framing (doc 03, from the paper's circular-orbit `+37.79″`) is superseded by this structural result: the discrepancy is not a numerical near-miss but a sign-and-factor-of-6 structural difference.

The paper's headline `44.39″/century` (Corda's `πm/M`) is the *reduced-mass period correction*, which causes **no** perihelion precession for a closed two-body ellipse — it is not a genuine perihelion advance at all (see [doc 05 §1](05_nbody_orbital_mechanics.md)).

Decomposition of the framework's full dual (one of several choices; the framework itself does not separate the prediction this way):

- Reduced-mass (Corda) contribution: `+44.66″/century`
- Contribution from the `(1 − GM/(c²r))` factor in Eq. (h4): **`−6.87″/century`** (this *individual contribution* is negative because the factor *reduces* the gravitational attraction — but it adds to the larger positive reduced-mass contribution to give a positive net)
- Net: `+44.66 − 6.87 = +37.79″/century` (positive forward advance)

The `(1 − GM/(c²r))`-induced contribution has the *opposite sign* from GR's `+42.99″/century` *purely-relativistic* contribution — meaning the framework's relativistic correction *subtracts* from the reduced-mass advance rather than adding to it as GR does. But the *total* framework prediction (`+37.79″`) is positive perihelion advance, not recession.

This is a substantive finding for the GPS author report's open question Q1; see [`04_findings_and_GPS_Q1.md`](04_findings_and_GPS_Q1.md).

## Status

| # | Document | Content | Status |
|---|---|---|---|
| 01 | [`01_setup_and_verification.md`](01_setup_and_verification.md) | Paper's `K = π²/(2m) + mc² + V + V²/(2mc²)` → Eq. (h3) → Eq. (h4) chain | ✅ all Wolfram-MCP verified |
| 02 | [`02_orbital_dynamics.md`](02_orbital_dynamics.md) | Eqs. (h5)–(h8) (orbital speed, period, angular velocity) | ✅ all Wolfram-MCP verified |
| 03 | [`03_numerical_predictions.md`](03_numerical_predictions.md) | Numerical Δφ for Mercury under all three paper predictions + GR + observed (paper's circular-orbit heuristic) | ✅ computed |
| 04 | [`04_findings_and_GPS_Q1.md`](04_findings_and_GPS_Q1.md) | Honest scoping + cross-reference to GPS Q1 | ✅ flagged |
| 05 | [`05_nbody_orbital_mechanics.md`](05_nbody_orbital_mechanics.md) | Proper eccentric-orbit apsidal-angle calc (`−1/6` of GR, exact) + N-body quantification + "1/6 factor" / Van Flandern connection | ✅ **the decisive result** |

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

This campaign reproduces standard Mercury perihelion math (a 19th-century calculation refined by Einstein 1915) and compares it against the framework's predictions; the value to the repository is the explicit Wolfram-MCP verification of each step and the per-prediction breakdown.

**The finding** is that the framework's full dual prediction (`+37.79″/century`) under-predicts the observed Mercury perihelion advance (`≈43″/century`) by `~12%`. This is *not* a precision-level rule-out: the paper's text claims only that Corda's value "well approximates" GR and the observed value (i.e., few-percent agreement), and a 12% disagreement between the framework's structural prediction and observation is consistent with that level of claimed precision but does not improve on it.

**The result is also subject to a two-body vs N-body caveat.** The 43″ observational residual is computed by *subtracting* Newtonian planetary perturbations (~531″/century from other planets on Mercury) from the total observed precession (~5600″/century). The subtraction assumes Newtonian gravity for the perturbing planets. If the framework's modified Newtonian gravity (with the `(1 − GM/(c²r))` factor) is applied consistently to *all* planetary perturbations on Mercury, the residual would change by a comparable amount to the framework's two-body correction itself. A consistent test of the framework would require running its modified gravity through the full N-body calculation — which this campaign has not done.

**Interpretive caveat.** The paper's `V²/(2mc²)` term was originally verified in the DRQM I (II.3) "potential-in-the-mass" kernel where `V` is an *EM* potential (between charged particles). Whether the same kernel is intended to apply when `V` is a *gravitational* potential is an interpretive choice the paper makes but does not justify from first principles. The author may intend the EM-only reading and consider the Mercury application a *separate* proposal; the result tests this specific substitution, not the kernel itself.

The framework may have other proposals for the gravitational extension (Q1 of the GPS author report lists four metric-based hypotheses; this paper proposes a categorically different kind of answer — modified Newtonian dynamics on flat space, with no curved-spacetime apparatus). This paper tests *one* such extension and finds it does not reach precision-level agreement with Mercury observation; the result is informative for the author's decision about which gravitational extension to develop further.
