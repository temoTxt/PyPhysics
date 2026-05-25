# Plan: issue #48 — proper-time third-term predictions for MeV bremsstrahlung

**Tracks:** [#48 — proper-time third-term predictions for precision MeV bremsstrahlung measurements](https://github.com/temoTxt/PyPhysics/issues/48). Sub-issue of [#42](https://github.com/temoTxt/PyPhysics/issues/42); companion to [#43](https://github.com/temoTxt/PyPhysics/issues/43) (different energy/intensity regime — see issue's *"Why this is worth a separate issue from #43"* table).
**Status:** plan only; no code or content committed yet. Branch `48-electromagnetism-proper-time-third-term-predictions-for-precision-mev-bremsstrahlung-measurements` is checked out.

## 1. Four findings from the pre-plan triage that materially shape execution

### 1a. **PR D dependency is now satisfied; the #42 campaign has shipped through PR K**

(Revised 2026-05-25 after Trey confirmed the push completed and the #42 branch advanced.) The full #42 campaign now lives on this branch via merge `f2a6010`: PR B (Ch. 11), PR C (Ch. 12), **PR D (Ch. 14)**, PR E (Ch. 16), PR F (Ch. 7), PR G (Ch. 4), PR H (Ch. 8), PR I (Ch. 9), PR J (Chs. 10 + 13 + 15), PR K (Chs. 1, 2, 3, 5 supplementaries). The expected commits `31c7779` (PR D problem 1 of 5: P14.1) and `8793f90` (PR D problem 5 of 5: P14.6) are present.

**The canonical bremsstrahlung derivation #48's prediction document is built on top of lives at [[Ch14_Radiation_by_Moving_Charges#Problem J3e-P14.6 — Bremsstrahlung from a linearly decelerating charge]].** That document already establishes:

- The classical Liénard `γ^6 a^2` formula for linear motion (contrasts with the `γ^4` for circular motion in P14.5).
- The proper-time third-term scaling `|E_3rd|/|E_acc,classical| ~ (u/c)^2` at non-relativistic intensity → **order unity at MeV** electron energies → **dominates** in the ultra-relativistic limit.
- A table of (`u/c`) → ratio at chemistry-scale, 10 keV, 100 keV, MeV, and ultrarelativistic.
- The verdict `⚠` — proper-time produces a quantitatively distinct radiated field for any non-circular acceleration; the longitudinal polarisation signature is the experimental discriminator.
- The closing "Notes for author review" already explicitly identifies MeV bremsstrahlung at clinical-linac energies as a candidate follow-on to #43 — i.e., this issue #48 was anticipated in PR D's own commentary.

**Implication for #48's prediction document.** §2 reduces from "self-contained restatement of Eq. (7)" (the #43 pattern) to a wikilink pointer to the Ch14 P14.6 canonical derivation. The substantive work of #48 is **per-energy specialization** to clinical-linac 1 / 6 / 18 MeV against the specific bibliography of bremsstrahlung-spectrum measurements scaffolded in Stage 1, plus the cross-experiment summary. The §6 "PR-D reconciliation list" originally planned as a deferred-work item is now satisfied at the structural level by virtue of the canonical derivation existing; what remains is the per-experiment specialisation.

### 1b. **The third-term contribution is genuinely large at clinical-linac energies**

The third term of Eq. (7) carries a `(u·a)/b⁴` prefactor. At a typical clinical-linac electron energy of 6 MeV, the electron velocity is `u/c ≈ 0.997`, the Lorentz factor `γ ≈ 12.7`, and `b² = c² + u² ≈ 2 c²`. The third term's contribution relative to the standard acceleration field (term II of Eq. 7) scales as roughly `(u/c)² · |a_∥|/|a_⊥| · (s/b)⁻¹`, which for linear deceleration in a bremsstrahlung target — where the acceleration *is* longitudinal by construction — is order unity, not perturbative. This is the regime the issue's headline observation calls out: *"order unity at MeV electron energies (`u/c ~ 0.94`)"*. The calculation in §3 needs to be done carefully because the small-parameter expansions that hide the third term in low-energy or weak-field regimes are not available here.

### 1c. **The candidate experimental signature is a longitudinal polarisation component**

For linear deceleration the standard textbook Liénard–Wiechert radiation is purely transverse in the far-field (radiation-zone) limit; the polarisation vector lies in the plane spanned by `n̂ × a` where `n̂` is the direction from the source to the observer. The proper-time third term's geometric factor is `(r × (u × r))`, which for `u ∥ a` (linear deceleration) produces a non-zero projection along `u` itself — a *longitudinal* component that classical EM does not predict. This is the experimentally testable signature; the question is whether published bremsstrahlung-spectrum measurements at clinical energies are *polarisation-resolved* enough to discriminate it, or whether the comparison must be made indirectly (e.g., via the angular distribution or via the total spectrum's normalisation).

### 1d. **Bibliography candidates resolved via Crossref; three behind paywalls and three open**

The issue identifies five candidate references by name. After Crossref lookup (2026-05-24 in plan-revision pass) the actual reference set sharpens to six papers — the issue's "Pratt & Tseng review" is most usefully realised as *two* papers (Koch-Motz 1959 review + Tseng-Pratt 1971 cross-section), since Pratt and Tseng's seminal bremsstrahlung work is the 1971 cross-section paper rather than a *Rev. Mod. Phys.* review per se:

| Cite-key (proposed) | DOI | Title | License | Acquisition route |
|---|---|---|---|---|
| `koch1959_bremsstrahlung_review` | `10.1103/RevModPhys.31.920` | Bremsstrahlung Cross-Section Formulas and Related Data | APS default (paywall) | **Trey downloads** |
| `tseng1971_screened_bremsstrahlung` | `10.1103/PhysRevA.3.100` | Exact Screened Calculations of Atomic-Field Bremsstrahlung | APS default (paywall) | **Trey downloads** |
| `seltzer1986_bremsstrahlung_tables` | `10.1016/0092-640X(86)90014-8` | Bremsstrahlung energy spectra from electrons 1 keV–10 GeV | Elsevier (paywall) | **Trey downloads** |
| `aapm_tg51_1999_clinical_dosimetry` | `10.1118/1.598691` | AAPM TG-51 protocol for clinical reference dosimetry | **CC-BY 3.0** ✅ | open — automated |
| `salvat2019_penelope` (or similar) | (OECD/NEA report; not in Crossref) | PENELOPE: A Code System for Monte Carlo Simulation | OECD/NEA — usually open | automated; fall back to local if NEA blocks |
| `nist_egsnrc_reference` | (NIST tech report) | EGSnrc Monte Carlo reference data | NIST — open | automated; specific tech-report DOI resolved during Stage 1 |

Three of six papers require Trey to download (institutional access via APS / Elsevier). If repository-public commit of the markdown extracts triggers copyright concerns — even fair-use academic quotation — the per-paper resolution per the campaign's existing two-tier policy is `pdf_status: acquired`, PDFs local-only, marker-pdf-converted markdown committed. *If the marker-pdf markdown of any specific paper exceeds the fair-use envelope (e.g. extensive equation reproduction from Seltzer-Berger 1986's tables), the conservative fallback is to commit only a verdict-level summary of the paper's relevant content rather than the full conversion.* Trey raised the option of making the repository private if needed; defer to him on per-paper basis as the markdown conversions land.

## 2. What this plan does NOT do

(Revised 2026-05-25 — PR D dependency now satisfied per §1a above.) The prediction document does **not** re-derive the proper-time third term self-contained; it wikilinks to [[Ch14_Radiation_by_Moving_Charges#Problem J3e-P14.6 — Bremsstrahlung from a linearly decelerating charge]] and specializes from there. The previously-planned §6 PR-D reconciliation list is folded into the document's structure from the outset rather than being a deferred-work item.

Per the issue's own scope guards: this plan does NOT engage clinical-physics interpretation (dosimetry), Born-approximation refinements that already account for the proper-time third term (we compare to *current* tabulated cross-sections), or radiation-reaction effects (negligible at clinical-linac energies; this is a pure Liénard–Wiechert third-term test, distinct from #43).

Per the campaign convention: this plan does **not** modify the per-problem template or per-chapter conventions established in PR 0 + PR A; #48's prediction document follows the same voice and Crocco-substantive-gating conventions but uses the `Experimental_Comparisons/` structure (cross-experiment analysis document, not per-problem) introduced by #43.

This plan does **not** propose acquiring PDFs that are behind paywalls without Trey's per-paper authorisation. The Crossref DOI lookups are non-controversial; PDF download routes (publisher direct, arXiv preprint, sci-hub-or-similar, library access) need a per-source decision.

## 3. Staged execution

### Stage 1 — Bibliography scaffolding (pragmatic)

For each of the 3–5 reference candidates in §1d:

1. Resolve DOI + Crossref metadata.
2. Scaffold bib stub via `uv run python Roadmapping/History/Bibliography/scaffold_bib_note.py --cite-key <key> --type retrospective --from-doi <doi>`.
3. Hand-edit to add tags (`#era/forward`, `#thread/electromagnetism`, `#gill-silent`, `bremsstrahlung`, `clinical-physics`, plus per-paper-specific topic tags) and `gill_corpus_overlap: [Two_Mathematically_Equivalent_Versions_of_Maxwells_Equations]`.
4. PDF acquisition: per Trey's resolved §6 Q2 below. Default — PDFs local-only, marker-pdf-converted markdown committed; `pdf_status: acquired`. PDFs that *cannot* be accessed even locally get `pdf_status: unavailable` and the bib stub records the access status.
5. Convert all acquired PDFs in one CPU run (Stage 1's slow leg; ~10–30 min/paper on CPU; the local `llama-server` continues to hold the A4500 VRAM).
6. Regenerate `Historical_Papers/Acquisition_Tracker.md`.

Proposed initial cite-keys (subject to confirmation in Stage 1 after Crossref lookup):

| Cite-key | Paper |
|---|---|
| `pratt1973_bremsstrahlung_review` | Pratt & Tseng, *Rev. Mod. Phys.* 1973 |
| `seltzer1986_bremsstrahlung_cross_sections` | Seltzer & Berger, *At. Data Nucl. Data Tables* 1986 |
| `salvat2019_penelope` (or similar) | Salvat PENELOPE OECD/NEA report (most recent edition) |
| `nist1985_egsnrc` | NIST EGSnrc reference data |
| `aapm_tg_NN_<topic>` | One AAPM Task Group report covering clinical-linac bremsstrahlung-spectrum measurements (specific TG number selected at scaffold time) |

### Stage 2 — Prediction document scaffold

Create `Roadmapping/Electromagnetism/Jackson/Experimental_Comparisons/bremsstrahlung_MeV_spectra.md` following the same structural template established by #43's [`radiation_reaction_2018.md`](../../Roadmapping/Electromagnetism/Jackson/Experimental_Comparisons/radiation_reaction_2018.md):

- §1 Selection provenance (Crocco-substantive — `<!-- TODO -->`).
- §2 The proper-time Liénard–Wiechert third term (pragmatic): §2.1 field expressions (quoted from Eq. 7); §2.2 when the third term contributes (linear deceleration → `u · a ≠ 0` instantaneously, not just cycle-averaged); §2.3 longitudinal polarisation prediction (substantive interpretation `<!-- TODO -->`).
- §3 Per-energy comparison: subsections for 1 MeV, 6 MeV, and 18 MeV clinical-linac energies (the issue's three target energies). Each subsection has a *geometry* paragraph (pragmatic — restated from the bib-stub paper geometries plus the linear-deceleration kinematics of bremsstrahlung in a high-Z target), a *proper-time prediction setup* paragraph (mix), a *quantitative prediction* block (Crocco-substantive `<!-- TODO -->`), and a *comparison + verdict* block (Crocco-substantive `<!-- TODO -->`).
- §4 Cross-energy summary table (Crocco-substantive verdict column `<!-- TODO -->`).
- §5 Crocco compliance (pragmatic).
- §6 PR-D reconciliation list (pragmatic; fires when PR D lands).
- §7 FINDINGS cross-post (conditional; deferred).

Voice conforms to [`VOICE_MATCH_GILL.md`](../../Roadmapping/Tooling/VOICE_MATCH_GILL.md). The "mathematically equivalent but not physically equivalent" signature is the natural anchor when contrasting the proper-time prediction against the Born-approximation cross-section.

### Stage 3 — Quantitative third-term calculation (substantive — gated)

The §3 quantitative-prediction blocks compute:

1. The proper-time third-term contribution to the bremsstrahlung **angular distribution** `dW/dΩ` at each of 1, 6, 18 MeV.
2. The predicted **longitudinal polarisation component** as a function of emission angle.
3. The corresponding **modification to the integrated cross-section** (likely small; the dominant signature is angular / polarisational).
4. Comparison against the Born-approximation tabulated values (Seltzer-Berger 1986 + PENELOPE).

This calculation is doable analytically for the linear-deceleration limit; full numerical evaluation across the bremsstrahlung spectrum requires a Wolfram Language notebook (or equivalent). Companion notebook: `Roadmapping/Mathematica_Notebooks/Electromagnetism/Bremsstrahlung_MeV_proper_time.wl`.

**Tooling note.** The Wolfram MCP server (`mcp__wolfram__WolframLanguageEvaluator`) appears disconnected in the current session per the [system-reminder posted 2026-05-24]; the campaign's standard "MCP verification inline" pattern from PR 0 / PR A may need to fall back to manual Mathematica runs from the .wl notebook against a local Wolfram installation (or wait for MCP to come back). This is not a blocker — the .wl notebook is the source-of-record by convention, and the inline MCP block is a redundant verification — but worth flagging.

Per Crocco §5: the calculation steps themselves are *pragmatic* (mechanical application of Eq. 7 with substitution into bremsstrahlung kinematics); the *physical interpretation* of the longitudinal polarisation as a discriminator between proper-time and classical EM is *substantive* and goes behind the `<!-- TODO -->` block. The verdicts (✅/⚠/❌) are substantive.

### Stage 4 — FINDINGS cross-post (conditional, deferred)

If any §3 verdict comes back ⚠ or ❌, add to [`Roadmapping/Equation_Verification/FINDINGS_for_author_review.md`](../../Roadmapping/Equation_Verification/FINDINGS_for_author_review.md). Same format-extension caveat as #43's §7: prediction-vs-experiment findings are structurally different from internal-derivation findings, and the extension is best worked out in PR D's discussion before the first such finding lands.

## 4. What this PR ships

Realistic scope when Stages 1 + 2 + 3-pragmatic-content complete:

- [ ] This plan doc.
- [ ] 3–5 bib stubs under `Bibliography/Retrospective/`.
- [ ] Per-paper converted markdown (Stage 1's slow leg).
- [ ] Acquisition tracker regenerated.
- [ ] Prediction document at `Roadmapping/Electromagnetism/Jackson/Experimental_Comparisons/bremsstrahlung_MeV_spectra.md` — §1 (TODO), §2.1 + §2.2 filled, §2.3 (TODO), §3 per-energy structure with geometry + setup paragraphs filled, quantitative-prediction + verdict blocks (TODO).
- [ ] (Optional, possibly) Companion `.wl` notebook scaffold — empty cells with section headers matching §3's calculations.
- [ ] Jackson README updated (add row to the Experimental_Comparisons table).

Substantive blocks (§1 selection provenance, §2.3 longitudinal-polarisation interpretation, §3 quantitative predictions per energy, §3 verdicts per energy, §4 verdict column + narrative) remain `<!-- TODO -->` for the substantive reviewer to fill — see §5 below.

**Issue #48 stays open** after this PR merges until the substantive blocks land.

## 5. Substantive reviewer — resolved: @Zorns-Lemmon

Per Trey (2026-05-24): route Stage-3-substantive to `@Zorns-Lemmon`, the same substantive reviewer engaged on #43. The MeV-clinical-bremsstrahlung regime is outside Zorns's cosmic-ray-physics primary specialty, but bremsstrahlung physics generally — including the angular and polarisation structure of radiation from accelerating charges — is foundational to high-energy and cosmic-ray work. Zorns is well-positioned to evaluate the longitudinal-polarisation prediction's plausibility and the comparison against tabulated Born-approximation cross-sections, even if specific clinical-physics measurement-protocol details (TG-51 calibration nuances, etc.) require additional consultation.

The handoff to Zorns mirrors the #43 pattern: a PR-comment summary explicitly enumerates the Crocco-compliance-gating tasks (read primary sources, validate geometry blocks, fill substantive blocks for §1 / §2.3 / §3 quantitative predictions / §3 verdicts / §4 verdict column, or flag where additional theoretical work is needed). Authoring of that comment is Stage 5 below.

If a clinical-medical-physics collaborator becomes reachable independently of this PR's timeline, a second-pass review can be requested then; the substantive blocks are version-controlled and a subsequent reviewer can append.

## 6. Resolved + open decisions

**Resolved as of 2026-05-24:**

1. ✅ Branch already created (`48-electromagnetism-proper-time-third-term-predictions-for-precision-mev-bremsstrahlung-measurements`).
2. ✅ Prediction document location follows the issue spec verbatim.
3. ✅ **Substantive reviewer: @Zorns-Lemmon** (Trey, 2026-05-24). See §5 above.
4. ✅ **Wolfram MCP reconnected** (system-reminder, 2026-05-24). Stage 3's MCP-verification pattern from PR 0 / PR A is available again; no fallback to local Mathematica needed.
5. ✅ **PDF acquisition split: Trey downloads paywalled (3 papers), automated for open-access (3 papers).** Per-paper assignments in §1d's table. Markdown-only commitment by default; if any marker-pdf extract exceeds fair-use envelope, fall back to verdict-level summary only. Repository-private route reserved as a last-resort fallback per Trey's offer.

**Open (need Trey's input before Stage 2/3 ramp-up):**

1. ~~PR-D commits not yet visible~~ **RESOLVED 2026-05-25:** Trey pushed; the full #42 campaign (PR B through PR K) merged into this branch as `f2a6010`. PR D's Ch. 14 P14.6 derivation is now the canonical anchor for #48 — see §1a above.
2. **Three target energies.** The issue calls out 1, 6, 18 MeV as "typical of clinical linacs." Lock to those three, or add (e.g.) 25 MeV (high-end clinical) or 0.5 MeV (sub-MeV) to expand the signature curve? *Recommendation: lock to the three the issue specifies; expansion happens in a follow-on PR if the comparison surfaces interesting energy-dependence.*
3. **`.wl` notebook scope.** Empty-cell scaffold ships in this PR, or full Wolfram computation deferred until Stage 3? *Recommendation: empty scaffold ships in this PR (mirrors the campaign convention; per-problem .wl notebooks are now plentiful on the merged branch); full Wolfram MCP-augmented computation in Stage 3 once Zorns engages.*
4. **Relationship to Ch15 P15.x backfill.** PR J landed [`Ch15_Bremsstrahlung_Virtual_Quanta.md`](../Roadmapping/Electromagnetism/Jackson/Ch15_Bremsstrahlung_Virtual_Quanta.md) — should #48's prediction doc cross-reference it, and conversely should Ch15 grow a forward-pointer to the Experimental_Comparisons document? *Recommendation: bidirectional cross-link in Stage 2.*
