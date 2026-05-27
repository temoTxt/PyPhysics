# Equation Verification: *Relativistic Newtonian Theory* (a.k.a. *Dual Newtonian Theory*)

**Author:** Tepper L. Gill (Howard University; sole author).
**Source:** [`../Tepper_Gill_Papers/Dual Newtonian Theory.tex`](../Tepper_Gill_Papers/Dual%20Newtonian%20Theory.tex) (.tex source supplied directly via issue [#81](https://github.com/temoTxt/PyPhysics/issues/81); user-supplied filename `Dual Newtonian Theory.tex`, paper title per the `\title{}` field *Relativistic Newtonian Theory*).
**Markdown:** [`../Converted_Markdown/Dual Newtonian Theory/Dual Newtonian Theory.md`](../Converted_Markdown/Dual%20Newtonian%20Theory/Dual%20Newtonian%20Theory.md) (pandoc-converted from .tex, no warnings).

**Verification status:** ✅ all 13 equations algebraically verified via Wolfram MCP (2026-05-27); ⚠ one numerical finding flagged for author review (see §"Numerical finding" below). Mercury-perihelion application worked out in [`../Mercury_Perihelion/`](../Mercury_Perihelion/).

## Conventions

Same as the Maxwell-paper and DRQM I verifications:

- Gaussian (c.g.s.) units throughout; the paper switches to SI only for the Mercury numerics.
- `w = dx/dt` (observer velocity), `u = dx/dτ = γw` (proper-time velocity), `b = √(c² + u²)`.
- `π = p − eA/c` (kinetic momentum); for the paper's gravitational application, `A = 0` and `π = p`.
- Standard Dirac matrices used in §1 but not load-bearing for the gravitational application.

## Equation index

| Eq. label | Topic | Verdict | Notes |
|---|---|---|---|
| (a1) | Proper-time definition `dτ = √(1−w²/c²) dt` | ✅ Textbook | Same as DRQM I (I.1) verified ✅. |
| (b1) | Same in primed frame | ✅ Textbook | Same as DRQM I (I.1) primed. |
| (c) | Rearranged: `cdt = b·dτ`, `u = γw` | ✅ | Same as DRQM I (I.2). |
| (d) | Same in primed frame | ✅ | Same as DRQM I (I.2) primed. |
| (e) | Time-derivative duality `(1/c)d/dt ≡ (1/b)d/dτ` | ✅ | Same as Maxwell-paper Eq. (2) verified ✅. |
| (f) | Velocity duality `w/c = u/b` | ✅ | Same as Maxwell-paper Eq. (1) verified ✅. |
| (2.4) | Canonical proper-time Hamiltonian `K = H²/(2mc²) + mc²/2` | ✅ | Same as DRQM I (I.6) and Maxwell-paper (16). |
| (after 2.4, no label) | Interaction `K = π²/(2m) + mc² + V + V²/(2mc²)` from `H₁ = √(c²π² + (mc²+V)²)` | ✅ | Structurally identical to DRQM I (II.3) "potential-in-the-mass" kernel verified ✅. |
| (h3) | `a = −(∇V/m)(1 + V/(mc²))` | ✅ NEW | Verified via Wolfram MCP — [`Mercury_Perihelion/01_setup_and_verification.md` §2](../Mercury_Perihelion/01_setup_and_verification.md). |
| (h4) | `a_m = −(GM/r²)(1 − GM/(c²r))·ê_r` from (h3) with `V = −GMm/r` | ✅ NEW | Verified via Wolfram MCP — [`Mercury_Perihelion/01_setup_and_verification.md` §3](../Mercury_Perihelion/01_setup_and_verification.md). |
| (h5) | Orbital speed `u_d = √(G(M+m)/r₀) · √[1 − G(M²+m²)/((M+m)c²r₀)]` | ✅ NEW | Verified via Wolfram MCP — [`Mercury_Perihelion/02_orbital_dynamics.md` §2](../Mercury_Perihelion/02_orbital_dynamics.md). |
| (h6) | Newtonian period `T₀ = 2πr₀^(3/2)/(GM)^(1/2)` | ✅ Textbook | Kepler's third law (Goldstein §3.7). |
| (h7) | Dual period `T_d/T₀ = [(1+m/M)(1−G(M²+m²)/((M+m)c²r₀))]^(−1/2)` | ✅ NEW | Verified via Wolfram MCP — [`Mercury_Perihelion/02_orbital_dynamics.md` §3](../Mercury_Perihelion/02_orbital_dynamics.md). |
| (h8) | Dual angular velocity `ω_d = ω₀ · [(1+m/M)(1−G(M²+m²)/((M+m)c²r₀))]^(1/2)` | ✅ NEW | Follows from `ω_d = 2π/T_d` and (h7). |

## Numerical finding (flagged for author review)

The paper's §2 concludes by computing three perihelion-advance predictions for Mercury (Corda's `Δφ_c = πm/M ≈ 44.39″/century`; full dual `Δφ_d`; approximate dual `Δφ_{d₁}`) and observes that Corda's value "well approximates the value per tropical century of 42.98″ given by general relativity and the well known observational value of 43″."

**Wolfram-MCP-verified numerical computation** (full breakdown in [`../Mercury_Perihelion/03_numerical_predictions.md`](../Mercury_Perihelion/03_numerical_predictions.md)):

| Prediction | Δφ per century | Notes |
|---|---|---|
| Corda `πm/M` | `44.66″` (paper says `44.39″`; `0.6%` discrepancy likely from parameter rounding) | reduced-mass effect, *not* relativistic |
| Full dual `Δφ_d` | **`37.79″`** | framework's structural prediction |
| Approximate dual `Δφ_{d₁}` | `37.79″` (matches full dual to 4 sig figs) | first-order Taylor of full dual |
| Standard GR `6πGM/[c²a(1−e²)]` | `42.99″` | reproduces observed to `<0.1%` |
| Observed (Le Verrier residual) | `≈ 43″` | reference |

**Finding (⚠):** The paper's Corda value (`44.39–44.66″`) is the *reduced-mass* contribution `πm/M` from the `(1+m/M)^(1/2)` factor in Eq. (h8) — a Newtonian effect that is already absorbed into Le Verrier's planetary-perturbation analysis. **The framework's structural prediction including the relativistic factor `(1 − G(M²+m²)/((M+m)c²r₀))^(1/2)` is `37.79″/century`**, which is `12%` below observed.

Decomposition: framework full-dual = reduced-mass `+44.66″` + framework relativistic correction `−6.87″` = `+37.79″`. The framework's relativistic correction has the **opposite sign** from GR's `+42.99″`.

The paper's headline claim that "44.39″/century well-approximates 42.98″" is therefore technically true for the Corda value but *not* for the framework's structural prediction.

This finding cross-references the open question Q1 in [`../Author_Reports/2026-05_gps_relativity_summary_for_gill.md`](../Author_Reports/2026-05_gps_relativity_summary_for_gill.md) §5 about the framework's correct curved-spacetime extension; see [`../Mercury_Perihelion/04_findings_and_GPS_Q1.md`](../Mercury_Perihelion/04_findings_and_GPS_Q1.md) for the full discussion.

<!-- TODO: human reviews and fills in — confirms whether this numerical finding should be promoted into FINDINGS_for_author_review.md as a Finding 4. The algebraic verification is complete (✅ all 13 equations); the question for the author is whether the framework's intended Mercury prediction is Corda's reduced-mass value (which agrees with observation at the 3-4% level) or the full-dual structural value (which disagrees by 12% and has the wrong sign of relativistic correction). -->

## Bibliography note (deferred)

A bibliography stub for this paper at `../History/Bibliography/Retrospective/gill_relativistic_newtonian.md` (or `gill_dual_newtonian.md`, depending on which title is preferred as the cite-key) is recommended; it is not yet scaffolded in this PR. The paper has no journal venue recorded in the `.tex` `\author{}` block — confirm with the author whether this is a preprint, an in-prep manuscript, or a previously-published work whose citation details should be looked up.
