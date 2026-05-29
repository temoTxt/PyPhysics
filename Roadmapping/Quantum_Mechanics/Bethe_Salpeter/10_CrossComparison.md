# PR J — Cross-comparison summary

**PR J — campaign-closing chapter.** This is the closing per-PR document of the Bethe–Salpeter precision-predictions campaign. PR J does not introduce new results; it consolidates the verdicts across PRs A–I into a single cross-comparison table and articulates the campaign's combined experimental status.

The campaign's load-bearing observation, articulated across all 28 results across the 10 PRs:

> *For non-relativistic atomic-physics observables, the proper-time canonical Hamiltonian `K = H_{0}^{2}/(2 m c^{2}) + m c^{2}/2` reduces — by algebraic definition — to `m c^{2} + p^{2}/(2 m)`, so every Schrödinger-level result of Bethe–Salpeter passes through unchanged. For precision-spectroscopy observables that depend linearly or quadratically on the electron magnetic moment `g_{s}`, the dual-theory prediction is `(g_{s}/-2)^{n}\,\times\,(\text{textbook QED leading term})` in each case. Six such `g_{s}`-dependent observables (hydrogen fine structure, 21-cm hyperfine, helium ³P fine structure, M1 transition rates, positronium ortho-para, muonium hyperfine) all share this single scaling structure; the campaign's "six-observable" `r_{e}` evaluation is therefore one structural fact applied to six observables, not six independent tests. The DRQM I §III.D cutoff value has, as of [PR #62](https://github.com/temoTxt/PyPhysics/pull/62), been refined from its as-published `r_{e}/r_{0} = 0.499857150068631` to the triangulated `r_{e}/r_{0} = 0.499420509912831\,7` (the joint-best-fit across all six observables, per Tepper Gill's 2026-05-25 author guidance that the as-published value was a uni-observable numerical search, not a derived quantity). At the triangulated `r_{e}`, the `(g_{s}/-2)^{n}\,\times\,\text{textbook}` formula reproduces measurement by construction at the precision the textbook leading-`g_{s}` route delivers — because the triangulated `r_{e}` is defined as the value that gives measured `g_{s}` across the joint fit. For the Lamb shift, the proper-time formulation's matrix elements, energy denominators, log-Bethe cutoff, and mass-renormalisation subtraction are all formulation-independent (per PRs A, B, E), so the framework's Lamb prediction reproduces the textbook Bethe-1947 calculation by construction and inherits its `\sim 42` MHz residual against measurement.*

The campaign's experimental status, honestly framed: *of 28 recorded results, 20 are formulation-independent (the framework agrees with textbook because it reduces to textbook), 6 are `g_{s}`-dependent observables evaluated at the triangulated `r_{e}` that test the back-fit self-consistency of the cutoff against measured `g_{s}` (jointly across all six), and 2 are at-Bethe-precision results that reproduce textbook predictions (Lamb shift + helium ground state). Zero of the 28 results constitute independent corroborations of the dual theory's content distinct from textbook QM.* This is the campaign's honest scope. The empirical resolution of the `r_{e}` finding is settled at the triangulated value per [PR #62](https://github.com/temoTxt/PyPhysics/pull/62); the remaining open question (tracked in [issue #54](https://github.com/temoTxt/PyPhysics/issues/54)) is whether a first-principles derivation from the dual-Dirac renormalisation prescription reproduces that cutoff. The question is about parameter-fitting consistency, not about experimental discrimination from standard QED.

---

## §1. PR-by-PR result inventory

| PR | Section | Result | Verdict | Notes |
|---|---|---|---|---|
| **A** | BS-§3 | Schrödinger equation, operator identity | ✅ | `K + V₀ = mc² + p²/(2m) + V₀` exactly (Wolfram MCP-verified) |
| **A** | BS-§4 | Hydrogen spectrum `E_{n}` | ✅ | Matches Rydberg `R_∞` to all orders in non-rel |
| **A** | BS-§5 | Bohr radius + radial eigenfunctions | ✅ | Identical to textbook |
| **A** | BS-§6 | `O(4)` symmetry + degeneracy | ✅ | Inherited from Schrödinger algebra |
| **B** | BS-§8 | Dipole matrix elements | ✅ | Operator + wavefunction identity |
| **B** | BS-§10 | Oscillator strengths + TRK sum rule | ✅ | `[r,K]=iℏp/m` survives rest-energy offset |
| **B** | BS-§13 | Radial integrals + selection rules | ✅ | Formulation-independent |
| **C** ⭐ | BS-§14.1 | Sommerfeld–Dirac spin–orbit | ✅ | Dual-Dirac FW reproduces leading `(Zα)²` |
| **C** ⭐ | BS-§14.2 | `2P₃/₂–2P₁/₂` splitting | ✅ at triangulated `r_e` | `-7` MHz residual at Bethe-estimate precision |
| **C** ⭐ | BS-§14.3 | Relativistic kinetic + Darwin | ✅ | Same coefficient in dual Dirac FW |
| **D** | BS-§15 | Nuclear recoil `μ/m_e` | ✅ | Two-body kinematic identity |
| **D** | BS-§17 | `(Zα)⁴` Pauli/FW expansion | ✅ | Weyl ordering reproduces textbook |
| **D** | BS-§18 | Bethe–Salpeter equation | ✅ | Apparatus inherited (ladder approximation) |
| **E** ⭐ | BS-§19 | Self-energy (Bethe 1947) | ✅ | Matrix elements + energy denominators formulation-independent |
| **E** ⭐ | BS-§20 | **Lamb shift `2S₁/₂–2P₁/₂`** | **✅ at Bethe precision** | `~1016` MHz vs `1057.845(9)` measured; `~42` MHz residual is standard Bethe-estimate gap |
| **E** ⭐ | BS-§21 | Uehling vacuum polarisation | ✅ | Identical one-loop |
| **F** ⭐ | BS-§22.1 | **Hydrogen 1S hyperfine (21-cm)** | ✅ at triangulated `r_e` | `-0.4` MHz residual at Bethe-estimate precision |
| **F** ⭐ | BS-§22.2 | Muonium + positronium (deferred) | (PR I) | — |
| **G** | BS-§24 | Photoionisation K-edge | ✅ | Cross-section formulation-independent |
| **G** | BS-§30 | M1 + E2 multipole transitions | ✅ at triangulated `r_e` | M1 inherits `(g_s/-2)^2 × textbook` via `g_s`; E2 formulation-independent |
| **G** | BS-§35 | Third-term effect on dipole approximation | ✅ | Doubly suppressed; below precision floor |
| **H** | BS-§47 | Helium HF `Z'=27/16` | ✅ | Variational identity |
| **H** | BS-§52 | Hylleraas correlation expansion | ✅ | `^{4}He` spin-singlet escapes `r_e` |
| **H** | BS-§60 | Full helium ground-state energy | ✅ at Bethe precision | `~6` meV residual |
| **I** | BS-§64 | Helium singlet-triplet exchange | ✅ | Formulation-independent (Fermi statistics) |
| **I** | BS-§72 | **Helium `^{3}P_J` fine structure** | ✅ at triangulated `r_e` *provisional* | `~kHz` residual after spin-spin/spin-orbit recombination (not derived in PR I) |
| **I** | BS-§80 | **Positronium ortho-para + muonium hyperfine** | ✅ at triangulated `r_e` (positronium *provisional*) | Muonium `~0.1` MHz residual; positronium requires positronium-specific QED absorbing correction (not derived in PR I) |

**Total: 28 results across 10 PRs**, of which:

- **20 ✅** unconditional — *but read this honestly:* every one of these is annotated "formulation-independent", "operator + wavefunction identity", "identical to textbook", "inherited from Schrödinger algebra", "Fermi statistics", etc. They are agreement-by-structure: the dual-theory framework's prediction *equals* the textbook prediction because the framework reduces to textbook QM at the level the observable engages. These results confirm the framework is *consistent* with textbook QM where it overlaps; they do **not** independently corroborate any distinctive dual-theory content.
- **6 ✅ at the triangulated `r_e`** (precision `g_{s}`-dependent observables) — *all six* are computed as `(g_{s}/-2)^{n}\,\times\,(\text{textbook leading term})` with `n \in \{1, 2\}`. The triangulated `r_{e}/r_{0} = 0.499420509912831\,7` (per [PR #62](https://github.com/temoTxt/PyPhysics/pull/62)) substitutes `g_{s} = -2.00231930\ldots` (the measured value, by which the joint fit defines `r_{e}`) and gets textbook agreement by construction. The "six independent precision discriminators" framing of earlier campaign drafts collapses to one structural fact (single-parameter joint fit of `g_{s}` via `r_{e}`) applied to six `g_{s}`-dependent observables.
- **2 ✅ at Bethe-estimate precision** (Lamb shift + helium ground state) — the dual-theory framework, where its non-relativistic-equivalent piece is exercised, reproduces the textbook Bethe-estimate calculation and inherits its residual against measurement. These are not independent corroborations either; they are recognitions that the framework, when reduced to textbook QM, gives textbook predictions.

**Zero of 28 results constitute independent corroboration of the dual theory's content distinct from textbook QM.** The campaign's value is in (i) tracing where the dual theory does and does not reduce to textbook, (ii) flagging the `r_{e}` finding's downstream self-consistency at high precision, and (iii) identifying the future-work direction (proper-time one-loop dual-QED) where a *distinct* dual-theory prediction could be derived and tested. The campaign does **not** itself produce such a distinct prediction.

---

## §2. The `r_{e}` self-consistency across six `g_{s}`-dependent observables (at the triangulated value)

**Reframe note.** Earlier drafts of this chapter called the table below "six *independent* precision discriminators all pointing to corrected `r_{e}`." That framing is withdrawn. Each row is a `g_{s}`-dependent observable for which textbook QED gives the predicted value as `f_{i}(g_{s}) \times E_{i,0}` with `f_{i}` a known polynomial in `(g_{s}/-2)` (linear for spin–orbit and Fermi contact contributions; quadratic for two-fermion spin–spin contributions like positronium ortho-para). The framework's prediction at the triangulated `r_{e}/r_{0} = 0.499420509912831\,7` (per [PR #62](https://github.com/temoTxt/PyPhysics/pull/62)) substitutes `g_{s} = -2.00231930\ldots` (the experimentally measured value, which is also by construction the value the joint fit defines `r_{e}` to give) and produces the textbook prediction per observable, which reproduces measurement at textbook leading-`g_{s}` precision.

The table therefore is **one calculation** (a `(g_{s}/-2)^{n}\,\times\,\text{textbook}` substitution at the triangulated `r_{e}`) repeated for six different `g_{s}`-dependent observables, not six independent tests of an underlying parameter. The campaign's load-bearing quantitative result, honestly stated, is:

| Observable | Framework prediction at triangulated `r_e` | Measurement | Residual | Source of residual | NIST MCP cross-check (issue #90) |
|---|---|---|---|---|---|
| Electron `g_{s}` (Finding 2 itself) | `-2.00231930` | `-2.00231930…` (CODATA) | `0` | matches by construction | `-2.00231930436256(35)` — CODATA 2018, `get_constant` |
| Hydrogen `2P₃/₂–2P₁/₂` (PR C) | `~10\,962` MHz | `10\,969.13(10)` MHz | `~7` MHz | Bethe-estimate floor (recoil + 2-loop) | `10\,969.049(15)` MHz — ASD levels, **theoretical/Ritz**, `get_levels` |
| Hydrogen 1S hyperfine (PR F) | `~1\,420.04` MHz | `1\,420.405\,751\,768(2)` MHz | `~0.4` MHz | sub-leading QED (α/π corrections + recoil + nuclear-structure) | not in ASD at this precision (measured value) |
| Helium `^{3}P₀-^{3}P_{1}` (PR I) | `~29\,616.95` MHz | `29\,616.952` MHz | `~kHz` | full QED residual *(provisional — spin-spin/spin-orbit recombination not derived in PR I)* | not derived via MCP (ASD He I levels exist; interval not auto-extracted) |
| Positronium ortho-para (PR I) | `~203\,389` MHz | `203\,389(2)` MHz | matches at Bethe precision | positronium-specific QED *(provisional — absorbing correction not derived in PR I)* | not on NIST (exotic atom; PDG / primary literature) |
| Muonium hyperfine (PR I) | `~4\,463.4` MHz | `4\,463.302\,776(51)` MHz | `~0.1` MHz | sub-leading QED | not on NIST (exotic atom; PDG / primary literature) |

**On the `NIST MCP cross-check` column (provenance, pragmatic per Crocco §1).** This column was populated programmatically via the `nist` MCP server ([issue #90](https://github.com/temoTxt/PyPhysics/issues/90), [`mcp_servers/nist_mcp/`](../../../mcp_servers/nist_mcp/)) on 2026-05-29; it is a retrieval cross-check, **not** a replacement for the `Measurement` column. Two honest caveats are load-bearing: (1) only the electron `g_{s}` is sourced from an authoritative independent reference (CODATA-2018 via `scipy.constants`); (2) the hydrogen `2P₃/₂–2P₁/₂` entry is the **theoretical/Ritz** level difference tabulated in the NIST ASD (`82\,259.2850014 − 82258.9191133` cm⁻¹, converted at `c`), not an independent measurement — it must not be read as corroborating the framework, and it differs from the measured `10\,969.13(10)` MHz by `~0.08` MHz. The MHz/sub-MHz splittings for hyperfine, helium `^{3}P`, positronium, and muonium are **not** carried by ASD at this precision and are deliberately left blank rather than back-filled with theoretical surrogates. <!-- TODO: human reviews and fills in — confirms the MCP cross-check column is correctly framed as a retrieval aid (g_s from CODATA; H fine structure as ASD-theoretical, not measurement), and that leaving the non-ASD splittings blank is the right honest choice rather than substituting theoretical values into a measurement comparison -->

**The pattern is consistent — *and consistent by construction.*** All six observables share the same `(g_{s}/-2)^{n}` scaling, and the triangulated `r_{e}` is by construction the joint-best-fit value across all six. The "six observables consistent" pattern is therefore a single-fact pattern (single-parameter joint fit) repeated six times, not six independent corroborations. Per-observable residual table including the full Pass A vs Pass B objective-function contrast: [`FINDINGS_for_author_review.md` Finding 2 update](../../Equation_Verification/FINDINGS_for_author_review.md#finding-2--drqm-i-eq-iii22-published-r_e-does-not-reproduce-the-experimental-g_e).

**An honest reading of what the campaign establishes.** *If* the dual-theory framework is taken to require a single-parameter `r_{e}` calibrated against `g_{s}`, *then* the triangulated `r_{e} \approx 0.4994205099\,r_{0}` is the joint-best-fit value across the `g_{s}`-dependent precision spectroscopy of hydrogen, helium, positronium, and muonium. This is a self-consistency statement about parameter-fitting within the framework, not an independent experimental discrimination of the dual-theory framework from standard QED. The framework uses the same `(g_{s}/-2)^{n}\,\times\,\text{textbook}` formula at the cutoff that textbook QED uses with measured `g_{s}`; the campaign's content is the cutoff identification, not a distinct prediction.

**The flagged finding's resolution path** is now: (1) the *empirical* resolution is settled at the triangulated `r_{e} \approx 0.4994205099\,r_{0}` per [PR #62](https://github.com/temoTxt/PyPhysics/pull/62) — the framework's `g_{s}`-dependent precision predictions are self-consistent at textbook leading-`g_{s}` precision when this cutoff is used; (2) the *theoretical* resolution — whether a first-principles derivation from the dual-Dirac renormalisation prescription reproduces this cutoff — remains open and is tracked in [issue #54](https://github.com/temoTxt/PyPhysics/issues/54). Finding 2's verdict shifts 🔴 → ⚠ (was "flagged", now "characterised") in [`FINDINGS_for_author_review.md`](../../Equation_Verification/FINDINGS_for_author_review.md); it would shift ⚠ → ✅ only if the first-principles derivation in #54 agrees with the triangulated value.

The campaign's *experimental status* with respect to the dual-theory framework's content distinct from standard QED is — at the level the campaign's apparatus exercises — *untested*, because every prediction it produces is either textbook-equivalent (20/28 results) or textbook-equivalent up to a single `g_{s}` substitution at the triangulated `r_{e}` (6/28 results) or textbook-equivalent at Bethe-estimate precision (2/28 results).

<!-- TODO: human reviews and fills in — confirms (a) the reframe from "six independent corroborations" to "one back-fit applied to six g_s-dependent observables at the triangulated `r_e`" is the correct honest framing post-[PR #62](https://github.com/temoTxt/PyPhysics/pull/62), (b) the resolution-path articulation (empirical settled at triangulated value; theoretical open in #54) is the right way to point the author's next step, and (c) the campaign's overall experimental status as "untested with respect to distinctive dual-theory content" is the right scope statement, not overselling -->

---

## §3. The Lamb shift result — *reproduction*, not endorsement

Earlier drafts of this chapter called the Lamb-shift result "the framework's strongest independent experimental endorsement." That phrasing is withdrawn. The honest reading of PR E ([`05_LambShift.md`](05_LambShift.md)) is:

- The proper-time framework reproduces the Bethe (1947) self-energy estimate **because** matrix elements, energy denominators, the Bethe-log UV cutoff, and the mass-renormalisation subtraction are all formulation-independent (per BS-§19 lines 47–50). The proper-time `\sim 1\,016` MHz prediction *is* the textbook Bethe-1947 calculation with a different label, not a distinct derivation.
- The `\sim 42` MHz residual against the measured `1\,057.845(9)` MHz is the standard Bethe-estimate-route residual that classical Bethe-1947 has had since 1947; it is the gap to the full one-loop QED calculation. The dual-theory framework inherits this residual exactly because it inherits the calculation exactly.
- An "endorsement" requires the framework to make a *distinguishable* prediction that is then *confirmed*. The Lamb-shift result here distinguishes the dual-theory framework from nothing — it reproduces the textbook calculation by document construction. The honest framing is therefore "the framework, where it reduces to the textbook Bethe-1947 calculation, gives the textbook Bethe-1947 prediction and inherits its textbook residual," not "the framework agrees with the Lamb-shift measurement."
- A full proper-time one-loop dual-QED calculation — required to push below the Bethe-estimate residual and to produce a *distinguishable* prediction from the standard one-loop result — is **out of scope** for this campaign. That calculation is where a Lamb-shift endorsement of the dual-theory framework could in principle live; this campaign does not produce it.

The Lamb shift result is `r_{e}`-independent because the leading log-Bethe contribution is `g_{s}`-symmetric; that observation is correct and continues to hold under the reframe. But it does not upgrade reproduction-of-textbook into endorsement-of-framework. The honest summary is: PR E shows the framework's non-relativistic-equivalent piece, when applied to the Lamb shift, recovers Bethe-1947 — exactly as the framework's algebraic non-rel reduction (BS-§3) guarantees it will.

---

## §4. Campaign-closing observation and outlook

The Bethe–Salpeter precision-predictions campaign — opened by issue [#50](https://github.com/temoTxt/PyPhysics/issues/50) as the experimental-discrimination counterpart to the [Griffiths pedagogical campaign (#49)](../Griffiths/) — has delivered:

1. **All four acceptance criteria of #50** met:
   - PR A merged — `Bethe_Salpeter/` scaffolded + §§1–7 complete (4 results), with `K → p²/2m + V₀` reduction verified ✅
   - PR C merged — §14 fine structure rederived; `2P₃/₂–2P₁/₂` splitting verdict recorded (now ✅ at the triangulated `r_e` per [PR #62](https://github.com/temoTxt/PyPhysics/pull/62)) ✅
   - PR E merged — §§19–21 Lamb shift treatment; explicit numerical prediction vs `1\,057.845(9)` MHz measured ✅
   - Flagged finding cross-posted to [`FINDINGS_for_author_review.md`](../../Equation_Verification/FINDINGS_for_author_review.md) ✅ — the `r_{e}` finding's consequences across all six precision observables, now evaluated at the triangulated value

2. **One consolidated honest-scope statement**: the dual-theory framework reproduces non-relativistic atomic-physics observables to all measured precision *because the framework reduces to non-relativistic QM by algebraic definition of `K`* (see [BS-§3](01_NonRelHydrogen.md#result-bs-3--schrödinger-equation-for-hydrogen)); it reproduces `g_{s}`-dependent precision observables at the triangulated `r_{e}` *because the triangulated `r_{e}` is the joint-best-fit value across the six observables, defined to give the measured `g_{s}`, and the predictions are `(g_{s}/-2)^{n} \times \text{textbook}`*; and it reproduces the Lamb shift at Bethe-estimate precision *because the framework's matrix elements and energy denominators are formulation-independent and the Bethe-1947 calculation passes through unchanged*. **The campaign does not produce a single prediction that is both distinct from textbook QM and tested against experiment.** That clause is the honest scope.

3. **One identified future-work direction**: a full proper-time one-loop dual-QED calculation. This is the natural place where a *distinct* dual-theory prediction (beyond the textbook Bethe-1947 leading-log) could be derived and meaningfully tested against the Lamb-shift measurement. The campaign identifies but does not produce this calculation. (Separately, a *first-principles derivation* of the cutoff `r_{e}` itself from the dual-Dirac renormalisation prescription is tracked in [issue #54](https://github.com/temoTxt/PyPhysics/issues/54); the empirical triangulation in [PR #62](https://github.com/temoTxt/PyPhysics/pull/62) settled the value, but a derivation of why the framework specifies that value remains open.)

4. **A noted reviewability gap** (Crocco §5 substantive-AI exposure): the ≥25 BS-section-number citations across the per-chapter documents rely on the Bethe & Salpeter (1977) text being available for verification. Per the campaign's bibliography policy, the Bethe–Salpeter reference is held local-only (`pdf_status: acquired`); reviewers without local access cannot independently audit the section-number chain. This is a known gap, recorded here.

**Honest closing.** The campaign's experimental status with respect to the dual-theory framework's content *distinct from* standard QED is — at the level the campaign's apparatus exercises — *untested*. At the triangulated `r_{e}` per [PR #62](https://github.com/temoTxt/PyPhysics/pull/62), the framework's `g_{s}`-dependent precision predictions are self-consistent at textbook leading-`g_{s}` precision. Whether the framework actually outperforms or differs from standard QED at higher precision is the question a future one-loop dual-QED calculation would have to answer; whether a first-principles derivation reproduces the triangulated `r_{e}` is the question [issue #54](https://github.com/temoTxt/PyPhysics/issues/54) tracks.

The campaign is closed at this honest scope, not at a stronger "framework passes the precision-experiment test" framing.

<!-- TODO: human reviews and fills in — confirms (a) the closing observation accurately summarises the campaign's experimental outcome, (b) the campaign's acceptance criteria from #50 are correctly identified as satisfied, and (c) the future-work direction (proper-time one-loop dual-QED) is the correct natural next step. The "honest framing" of the campaign — conditional on r_e resolution + Bethe-estimate precision floor — is faithfully recorded throughout this closing chapter -->
