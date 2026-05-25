# Implementation plan: Interim report for Tepper Gill (author review)

**Tracks:** issue [#58](https://github.com/temoTxt/PyPhysics/issues/58).
**Status:** plan; report drafting + LaTeX/PDF build pipeline to follow.
**Dependencies:** [#42 PR #51](https://github.com/temoTxt/PyPhysics/pull/51), [#49 PR #52](https://github.com/temoTxt/PyPhysics/pull/52), [#50 PR #53](https://github.com/temoTxt/PyPhysics/pull/53) — all three campaign PRs delivered as of 2026-05-25; this work consolidates them into an author-engageable summary.

This is **not** a multi-PR campaign. It is a single report-authoring task with three deliverables (markdown source + LaTeX + PDF) plus a small amount of supporting tooling (build script + `Author_Reports/` folder convention). The plan inherits the [§13 honest-framing discipline of #42](42-electromagnetism-jackson-proper-time.md#13-devils-advocate-review-and-what-we-cannot-honestly-fix) verbatim — the report's value to Tepper is that he can trust its scoping.

---

## 1. Goal and scope

Produce an interim report for Tepper Gill summarising the repository's verification + campaign work since 2026-05-11, with explicit author-side questions and decision points. The report's audience is a single known collaborator (a co-author of *DRQM I* and the foundational dual-theory papers); its purpose is to surface honest scoping, flagged findings, and open questions in a form Tepper can act on in one sitting.

**Honest framing of the request** (revised 2026-05-25 after the phone call). The report is **solicited**: Tepper requested an interim summary via a phone call with Trey. The campaign documents (#51, #52, #53) record the verification and per-PR findings in detail, but no consolidated author-engageable summary exists, and that summary is what Tepper has asked for. This is a deliverable to a known collaborator who has explicitly invited it — not an unsolicited peer note.

The implication for the report's framing: §1's cover note opens by acknowledging the phone call as the report's origin (no apology for the imposition, no peer-note-asymmetry hedge), then states what we have for him. The honest-scoping discipline is preserved (everything in §10 still applies); what changes is that we are not also apologising for the existence of the report.

### Why this earns its own thread

- **Author-engageable format.** The verification campaign's [`FINDINGS_for_author_review.md`](../../Roadmapping/Equation_Verification/FINDINGS_for_author_review.md) covers three substantive errata findings; the per-PR campaign documents (50+ files across #42 / #49 / #50) cover the load-bearing structural results. Neither is an author-engageable summary.
- **PDF deliverable.** Tepper receives a PDF, not a GitHub URL trail. The LaTeX/PDF build pipeline is the bridge between the markdown source-of-truth (which keeps `<!-- TODO -->` blocks for human-acceptance) and the deliverable (which does not).
- **Folder structure for any successor reports.** `Roadmapping/Author_Reports/` becomes the location *if* further author reports are written; the README documents naming conventions so the folder remains navigable. *Not* a claim that this report's build pipeline is the standard convention for hypothetical future reports — that generalisation is deferred until a second use-case actually appears (see §5 scope note).

### Explicit non-goals

- Not a journal manuscript. Working document for a known collaborator.
- Not a re-summary of any per-PR campaign document — cross-reference, do not duplicate.
- Not a pre-judgement of Tepper's answers. The report frames the questions; it does not propose Tepper's preferred answer for him.
- Not a new-campaign plan. Suggested next steps in the report's §7 are tactical (resolve specific questions), not strategic (open new campaigns).

---

## 2. Repository structure

New top-level folder, sibling to `Equation_Verification/` / `Electromagnetism/` / `Quantum_Mechanics/` / `History/`:

```
Roadmapping/Author_Reports/
├── README.md                              # convention for future reports
├── build_report.sh (or .py)               # single-command markdown → LaTeX → PDF
├── 2026-05_interim_for_gill.md            # markdown source-of-truth
├── 2026-05_interim_for_gill.tex           # LaTeX source (derived)
└── 2026-05_interim_for_gill.pdf           # built PDF (committed)
```

The folder's `README.md` documents the convention so subsequent author reports follow the same naming + build pipeline.

---

## 3. Report structure (the markdown source)

**Length budget: 3–7 pages in the final PDF, soft target 5.** Earlier drafts of this plan set a hard `[3, 5]` cap as a load-bearing constraint; that cap is relaxed because the campaign's honest-scoping caveats (back-fit caveat on branch (c), K-identity-is-by-definition caveat, Lamb-shift-is-reproduction-not-endorsement caveat, asymmetric absorbing-correction caveat) are *paragraph-scale* and the steel-man edits to PRs #51/#53 confirmed that compressing them costs fidelity. The honest-scoping discipline is load-bearing for the report; the page count is not. If the natural drafted length comes in at 6–7 pages because the caveats demand the space, the report is the right length and the caveats stay.

Soft section targets (at the 5-page nominal): §§1, 2, 7 ≈ ½ page each; §§3, 4, 5 ≈ ¾ page each; §6 the load-bearing Q1–Q8 section ≈ 1 page. If the report goes long, the page budget yields *before* the caveats do.

The build-step defensive check (§5 below) fails loudly **above 7** (overshoot is real) or **below 3** (undershoot suggests missing substance, most likely §6 Q1–Q8 needs more context per question or §4 conditional predictions are under-cited). Inside `[3, 7]` the check passes; the human-acceptance pass owns the final judgement on length. **Honest scoping never gets sacrificed to hit a length target.**

`Roadmapping/Author_Reports/2026-05_interim_for_gill.md` has the following sections, in order. **Emphasis revised 2026-05-25**: lean *into* process and initial results — especially the Jackson (#42) findings and the candidate experiments — and *away from* future-roadmap content. The author's primary value-add to Tepper is *what we found while doing the work*, not *what we propose to do next*; the §7 next-steps section is reduced to a brief closing list rather than a load-bearing section.

1. **§1 Cover note** (≈ ½ page). One paragraph — opens by acknowledging the phone call that prompted the report (Tepper's request), states what we have for him in the report, and notes what response would help us most. Neutral co-author-to-co-author register per §4 voice override. *No peer-note hedging — the report is solicited.*
2. **§2 How we got here — the process** (≈ ¾ page, *expanded*). The verification campaign's posture (systematic Wolfram-Mathematica re-derivation, Crocco compliance, branched-treatment workflow for flagged findings). What discipline we adopted, what we did not adopt, why. This section's load-bearing function is to establish that the work was done *carefully* — every claim downstream rests on whether Tepper trusts the process. References to [`CROCCO_COMPLIANCE.md`](../../Roadmapping/Tooling/CROCCO_COMPLIANCE.md), [`VOICE_MATCH_GILL.md`](../../Roadmapping/Tooling/VOICE_MATCH_GILL.md), [`Equation_Verification/README.md`](../../Roadmapping/Equation_Verification/README.md), and the three campaign plan files. Cover the four threads (verification + #42 + #49 + #50) in process terms — what each set out to do, not what each found (that's §3–§4).
3. **§3 Jackson highlights — initial results** (≈ 1 page, *expanded; promoted to lead-section*). The interesting Jackson (#42) results, in order of what we believe Tepper will care about most:
   - **J3e-P14.2 LW third term** (longitudinal radiation component absent from classical EM; podcast pick #2). Reference [Ch14_Radiation_by_Moving_Charges.md](../../Roadmapping/Electromagnetism/Jackson/Ch14_Radiation_by_Moving_Charges.md). This is the campaign's most novel structural prediction outside the precision-spectroscopy sector.
   - **J3e-P16.1 Abraham-Lorentz dissolution claim** (proper-time first-order `∂_τ` structure removes runaway + pre-acceleration pathologies; podcast pick #3). Reference [Ch16_Radiation_Damping.md](../../Roadmapping/Electromagnetism/Jackson/Ch16_Radiation_Damping.md). The claim is structurally clean; whether it *also* fixes the empirical defects (out of scope for #42) is what motivates the candidate experiments below.
   - **J3e-P5.13 Zorn surprise finding** (the non-zero result flagged via #42's [Zorn comment](https://github.com/temoTxt/PyPhysics/issues/42)). What the calculation showed and what we asked Zorn to verify.
   - **Honest scoping for each**: every Jackson result is labelled with what it predicts vs what the dual-theory framework *predicts unconditionally* (e.g., the Abraham-Lorentz dissolution is conditional on author intent — see §6 Q6).
4. **§4 Candidate experiments — what we found to look at** (≈ ¾ page, *new emphasis*). The two experimental sub-investigations the #42 campaign surfaced:
   - **#43 Cole / Poder / Wistisen 2018 radiation-reaction experiments** ([issue #43](https://github.com/temoTxt/PyPhysics/issues/43)). Three published experiments (Cole *PRX* 8, Poder *PRX* 8, Wistisen *Nat. Comm.* 9) in the GeV / extreme-intensity regime where the proper-time third term and dissipative coefficient `(u·a)/b⁴` make predictions distinguishable from classical Landau-Lifshitz. The two papers disagree on which classical limit fits, which itself is a hint that the proper-time formulation is a candidate alternative worth working out. **Status: scoped, not executed.**
   - **#48 MeV bremsstrahlung in clinical-linac regime** ([issue #48](https://github.com/temoTxt/PyPhysics/issues/48)). Precision-experiment regime distinct from #43: medical-physics calibration measurements at MeV electron energies, where the third term contributes at order unity and would predict a longitudinal polarisation component absent from classical Liénard-Wiechert. **Status: scoped, not executed.**
   - **Why these matter to Tepper**: the third term he and Zachary derived in *Two Mathematically Equivalent Versions of Maxwell's Equations* (Eq. 7, ✅ verified) makes operationally distinguishable predictions in two experimentally-accessible regimes. The campaign has done the scoping; we are flagging the candidates for his awareness, not asking him to execute them.
5. **§5 What we found unconditionally vs conditionally** (≈ ¾ page). Brief recap of what the campaigns showed *without* and *with* conditioning on the `r_e` finding's resolution. Three load-bearing structural results (verification ✅: Maxwell-paper structure, DRQM I §II FW reduction, dual-Dirac eqs. II.1–II.3) plus the six `r_e`-discriminator observables (the table from #50 PR J, **printed in full in §5** — see §6 below for why). The Bethe-estimate Lamb-shift result is recorded here as the framework's strongest independent endorsement (independent of `r_e`).
6. **§6 The `r_e` branch question** (≈ 1¼ pages, **lead question of the report; expanded 2026-05-25**). The single most load-bearing question we have for Tepper. Detailed treatment, not a one-line entry:

   - **§6.1 What the two branches are.** Restate the finding from [`FINDINGS_for_author_review.md`](../../Roadmapping/Equation_Verification/FINDINGS_for_author_review.md) Finding 2 in two paragraphs. Branch **(b)** is the as-published `r_e/r_0 = 0.499857150068631` from DRQM I §III.D, which yields `g_e = -2.0005714` (a ≈`1.75 × 10⁻³` disagreement from measured `-2.00231930…`). Branch **(c)** is the corrected `r_e/r_0 ≈ 0.499420510`, which reproduces the measured `g_e` and the textbook anomalous moment. The numerical reduction is given explicitly — both formula evaluation (`g_r = 2(1 - 4r_0/(2r + r_0))`) and the table of branch-(b) vs branch-(c) predictions for the six precision observables.
   - **§6.2 Why the question is load-bearing.** The same `r_e` value propagates through every `g_s`-linear (and `g_s²`-linear) observable in the campaign. The #50 PR J cross-comparison table — six independent precision atomic-physics observables (electron `g_s`, hydrogen 2P₃/₂–2P₁/₂, hydrogen 1S hyperfine 21-cm line, helium ³P₀–³P₁, positronium ortho-para, muonium hyperfine) — *consistently* exhibits the same pattern: branch (b) is in ~10⁻³ disagreement; branch (c) agrees at the Bethe-estimate precision floor. On the 21-cm hyperfine line (the most precisely measured atomic-physics frequency, 12 sig figs), branch (b)'s disagreement is ~6 orders of magnitude beyond measurement uncertainty. **The campaign's experimental status depends on which branch is the framework's intended choice.** We are not asking Tepper to pick the answer that best fits the data — we are asking whether the published value is the intended one, or whether it was a transcription error from a working notebook.
   - **§6.3 What we are asking for, concretely.** Q1 is decomposed into three sub-questions so the report makes the resolution actionable:
     - Q1a — Is `r_e/r_0 = 0.499857150068631` (branch (b), as published) the intended value of the framework? *If yes, the framework predicts a measurable ~10⁻³ disagreement with measured hyperfine, fine structure, and `g`; the campaign's experimental verdicts in #50 should be recorded as 🔴 (ruled out at current precision).*
     - Q1b — Is `r_e/r_0 ≈ 0.499420510` (branch (c), the value that reproduces measured `g_e`) the intended value? *If yes, the campaign's six precision verdicts collapse to ✅ at Bethe-estimate precision, and the as-published number is a transcription error to be erratum'd.*
     - Q1c — Is there a third possibility we have not considered (a different formula, a different derivation, a context where one value is correct and another emerges)? *If yes, we will rework the campaign's verdicts with the third value in mind.*
   - **§6.4 How the framework supports answering this.** Two routes: (i) retrieve the original working notebook in which `r_e` was computed and confirm the digit string; (ii) rederive `r_e` from the dual-Dirac equation's renormalisation prescription. The campaign cannot do (i) (no access to working notebooks); the campaign could do (ii) but only if the renormalisation prescription is documented (this is closely related to Q5 operator-ordering — see §6.5). Tepper is the source-of-record for both routes. Even a one-line confirmation ("the intended value is (b)" or "(c)" or "I'll check") is sufficient for the campaign to move forward.

   The structural emphasis in §6 is intentional and analogous to the campaign-wide weight assigned to the Lamb shift and the hyperfine line in #50: this is *the* question whose resolution most efficiently determines the dual-theory framework's experimental status.

7. **§7 Secondary findings + questions** (≈ ½ page, *condensed from previous §6*). The remaining errata findings and questions, after §6 has handled `r_e`. Three to four items:
   - Q2 — Confirm or correct the Maxwell Eq. 24 erratum (Finding 1: missing `c`, missing `V²/(2mc²)`).
   - Q3 — Confirm or correct the TCEP Eq. (4.16) sign typo (Finding 3).
   - Q4 — `r_μ` and `r_p` numerical values for DRQM I Eq. III.23 (so muonium/proton magnetic-moment predictions are derivable).
   - Q5 *(optional, only if space permits within the [3, 7] page bound)* — Operator-ordering choice in DRQM I §II.D FW reduction (sub-leading `(Zα)⁴` results depend on it; this is the question §6.4 route (ii) needs to make progress on `r_e`).
   - **Confirm the radiation-reaction dissolution interpretation (J3e-P16.1)** — moved into §3 Jackson highlights inline rather than asked as a standalone question, since it is structurally about a Jackson result rather than an erratum.
   - **Removed from this report**: Q7/Q8 roadmap-and-apparatus questions (proper-time one-loop dual-QED). Documented in [#55](https://github.com/temoTxt/PyPhysics/issues/55) for the record; not surfaced here. Report emphasis is process + initial results.
8. **§8 Closing note** (≈ ¼ page, *reduced*). One short paragraph closing the report — thanks for his time and the phone call, mention what would help us hear from him (Q1 first; Q2–Q4 secondary), and an honest line about whether to keep the campaigns running on the current trajectory. **Not a "next steps" plan.** The campaigns continue if the resolution of Q1 makes that worthwhile.

**Section-length budget summary** (soft targets, total ≈ 5½ pages at the nominal; within the `[3, 7]` page bound):

| § | Topic | Length |
|---|---|---|
| §1 | Cover note (solicited; acknowledges phone call) | ½ page |
| §2 | Process (campaign discipline) | ¾ page |
| §3 | Jackson highlights (Jackson initial-results emphasis) | 1 page |
| §4 | Candidate experiments (#43, #48) | ¾ page |
| §5 | Unconditional vs conditional results (incl. 6-observable `r_e` table) | ¾ page |
| §6 | **The `r_e` branch question — lead** (Q1a/b/c + framework support) | **1¼ pages** |
| §7 | Secondary findings + Q2–Q5 | ½ page |
| §8 | Closing | ¼ page |
| **Total** | | **~5½ pages** |

The report's structural weight is now distributed across three load-bearing emphases: §3 + §4 (process/results/experiments, ~1¾ pages); §6 (`r_e` branch question, 1¼ pages); §5 + §7 (unconditional results + secondary findings, ~1¼ pages). §1, §2, §8 are framing material (~1½ pages combined). §6's promotion to dedicated lead-question status is the key change in this revision: the report actively *seeks* the `r_e` resolution rather than listing it as one of many questions.

---

## 4. Per-section authoring discipline

The report is **substantive AI** end-to-end. Per the Crocco compliance protocol, every paragraph carries a per-paragraph `<!-- TODO: human reviews and fills in -->` block in the markdown source. The committed `.md` retains these for the human-acceptance pass; the built PDF does not (per the build-script's TODO-stripping step in §5 below).

Voice for the **AI draft (Phase 2): neutral co-author-to-co-author register** — measured, deferential, peer note. *Not* Gill's published-paper voice. Rationale (confirmed 2026-05-25 via plan-acceptance question): addressing a senior co-author in their own writing style risks reading as parody or pandering more easily than as deference, especially in an unsolicited peer note. The default is the safer register; the human-acceptance pass in Phase 4 can pull individual paragraphs toward [`VOICE_MATCH_GILL.md`](../../Roadmapping/Tooling/VOICE_MATCH_GILL.md) if a specific section reads better there, but the AI draft does not start in that voice.

This is an **explicit override of the campaign-wide voice discipline** (which the #42, #49, #50 per-PR documents use uniformly). The override is scoped to the report and is documented here so the Phase 3.5 devil's-advocate self-review can verify the AI draft does not drift back toward the campaign voice by inertia.

Honest scoping in every section: algebraic identities labelled as such; conditional predictions labelled as such; out-of-scope claims labelled as such. The model is the [`BetheSalpeter_S3.wl` honest-scope linter note](../../Roadmapping/Mathematica_Notebooks/Quantum_Mechanics/BetheSalpeter_S3.wl) — verbatim discipline.

---

## 5. Build pipeline

A single-command build from the markdown source to the deliverable PDF.

**Scope note (revised).** Earlier drafts justified the six-step pipeline as "convention-setting for future reports" (cf. §1.2). That justification is downgraded: there is no second author report planned, no second recipient, and no second use case. Per the CLAUDE.md *don't design for hypothetical future requirements* discipline, the build script is scoped to *this* report only. If a second report appears later, the script will be generalised then with the benefit of an actual second use-case. The pipeline is still six steps because each one solves a concrete problem (TODO-strip, math passthrough, page-count guard, etc.), not because it is forward-engineered for hypothetical reports.

Pipeline:

1. **Markdown → LaTeX (pandoc)**. `pandoc 2026-05_interim_for_gill.md -o 2026-05_interim_for_gill.tex` with appropriate template + filter flags. Math passes through unchanged (`$...$`, `$$...$$`). Tables convert via pandoc's reliable markdown-table → LaTeX-table conversion. Code blocks via `listings` or `verbatim`.
2. **TODO-comment stripping**. Markdown HTML comments `<!-- TODO: human reviews and fills in -->` must not leak into the PDF. Either:
   - Use a pandoc filter (Lua filter) that drops HTML comments at conversion time, or
   - Pre-process the `.md` to strip `<!-- ... -->` lines before pandoc, or
   - Use pandoc's `--strip-comments` if available.
   The build script must **fail loudly** if any `<!-- TODO -->` is present at PDF-build time after the strip step (defensive check).
3. **LaTeX → PDF**. `pdflatex` or `xelatex`, run twice (for cross-references). Document class `article` with `geometry`, `hyperref`, `amsmath`, `booktabs`. PDF metadata set via pandoc `--metadata` or LaTeX `\hypersetup`: Title, Author = "Trey Morris with Claude Opus 4.7", Date, Subject = "Interim report for author review".
4. **Cross-reference handling** (revised 2026-05-25 via plan-acceptance question). GitHub-relative `../blob/main/...` links cannot be followed from a PDF. The PDF keeps the ~dozen most-cited paths as `\footnote{See \texttt{path} in the repository.}` references; lower-frequency cross-references in the markdown source are *dropped from the PDF* while remaining in the markdown source-of-truth (where GitHub renders them as live links). The drop-set is decided in Phase 3.5 (devil's-advocate pass) by counting citation frequency in the draft; the human-acceptance pass in Phase 4 confirms.

   The rule for "is this safe to drop": if the report's claim *survives without the cross-reference* (i.e., the claim is stated and the reader can act on it without following the link), the cross-reference is decorative and can be dropped from the PDF. If the claim *requires* the link to be verifiable (e.g., the six branched verdict numbers in §4 require the PR #53 table), the cross-reference is load-bearing and stays. This is the same scoping discipline the plan applies to the report's substantive content: keep what is load-bearing, drop what is decorative.
5. **Bibliography integration** (if cited). Reuse the existing [`Roadmapping/History/Bibliography/build_bibtex.py`](../../Roadmapping/History/Bibliography/build_bibtex.py) pipeline to produce a `bibliography.bib` and reference it from the LaTeX via `\bibliography{}`. Cite only existing bib stubs; do not fabricate citations.
6. **Page-count defensive check.** After PDF build completes, the script runs `pdfinfo` (or `pdftk … dump_data`, or `qpdf --show-npages`) to extract the page count and fails loudly if outside `[3, 7]` (per the relaxed §3 budget). Sample check:

   ```bash
   pages=$(pdfinfo "$PDF" | awk '/^Pages:/ {print $2}')
   if [ "$pages" -lt 3 ] || [ "$pages" -gt 7 ]; then
     echo "ERROR: PDF page count $pages is outside the [3, 7] budget per plan §3." >&2
     exit 1
   fi
   ```

   This check is symmetric with the TODO-comment defensive check: both fail the build rather than silently shipping a non-conforming PDF. The window is `[3, 7]` rather than the tighter `[3, 5]` an earlier draft proposed; honest scoping is load-bearing and is allowed to push the page count into the 6–7 range when needed (see §3).

System dependencies (reused from Manim sub-project per [`CLAUDE.md`](../../CLAUDE.md)): `pandoc`, `texlive-latex-base`, `texlive-latex-extra`, `texlive-pictures`, plus `poppler-utils` (for `pdfinfo`). Document install command in the script's docstring.

---

## 6. PR sequencing

This task is delivered as a **single PR** (not a multi-PR campaign). The phase ordering is **deliberate**: no PDF is built and committed before the human-acceptance pass completes, so Tepper never sees a finished-looking PDF that is actually a pre-review AI draft (which would materially violate Crocco §5 — see §10 honest framing).

| Phase | Scope | Status |
|---|---|---|
| Phase 1 | `Roadmapping/Author_Reports/` scaffold + `README.md` convention + build-script skeleton (script tested on a stub `.md`, but PDF *not* committed yet) | pending |
| Phase 2 | `2026-05_interim_for_gill.md` draft — §1 through §7 with per-paragraph TODO blocks | pending |
| Phase 3 | Build pipeline *dry-runs* against the draft `.md` (verify pandoc + LaTeX + page-count + TODO-strip mechanics on the actual content) but does **not** commit a PDF yet | pending |
| Phase 3.5 | **Devil's-advocate self-review** of the Phase 2 draft against the campaign's §13 honest-framing discipline + **steel-man patches** applied to the `.md`. Analogous to the workflow that pulled PRs #51 and #53's per-chapter docs back into alignment with §13 after they drifted above it. The self-review's output (the critique + the steel-man diff) is written to a *transient* file `2026-05_interim_for_gill_self_review.md` for the Phase 4 human reviewer to consult; the file is **deleted before PR open** so the `Author_Reports/` folder stays clean. The critique-and-patch content survives in git history via the Phase 3.5 commit message | pending |
| Phase 4 | Human-acceptance pass: review remaining TODO blocks (and the transient self-review file), edit prose, confirm honest scoping after Phase 3.5 patches, **build and commit final PDF**, **delete the transient self-review file** | requires human review |
| Phase 5 | PR opened (closes #58) | pending |

Phases 1–3 are AI-authored work that does not ship a PDF; Phase 3.5 is the AI-authored devil's-advocate / steel-man gate that catches honest-framing drift before human review; Phase 4 is human-authored (the load-bearing review pass — the committed `.md` is a draft until this pass completes, and the PDF is built and committed only at the end of Phase 4); Phase 5 closes the loop.

---

## 7. Acceptance criteria

Per [issue #58](https://github.com/temoTxt/PyPhysics/issues/58):

1. Markdown source at `Roadmapping/Author_Reports/2026-05_interim_for_gill.md` (criterion 1 of issue).
2. Report structure per §3 above (criterion 2 of issue): §§1–7 with the eight Q1–Q8 open questions explicit.
3. Honest scoping discipline throughout (criterion 3): per-paragraph TODO blocks, no overselling, scoping labels on every identity / prediction / out-of-scope claim.
4. Voice match per [`VOICE_MATCH_GILL.md`](../../Roadmapping/Tooling/VOICE_MATCH_GILL.md) (criterion 4).
5. Crocco compliance (criterion 5): substantive AI end-to-end, per-paragraph TODO blocks, no AI-fabricated citations.
6. Cross-references present (criterion 6): the seven specific cross-references named in #58's body.
7. LaTeX/PDF deliverable (criterion 7): `.tex`, `.pdf`, and `build_report.sh` (or `.py`) all committed; build is single-command and rebuildable (see §13 reproducibility scope; *not* byte-identical determinism); PDF carries correct metadata; no `<!-- TODO -->` leaks into the PDF.
8. **PDF page count is in `[3, 7]`** per §3 length budget (relaxed from an earlier `[3, 5]` cap to protect honest-scoping caveats); the build's `pdfinfo`-based defensive check (§5) passes.
9. **The PDF is committed only at the end of Phase 4 (human-acceptance pass)**, not during Phases 1–3. The build script's dry-run capability is exercised in Phase 3 to verify the pipeline mechanics, but the deliverable PDF — the artifact Tepper would actually receive — does not enter the repository until a human has reviewed the AI-drafted `.md` (and the Phase 3.5 devil's-advocate self-review's steel-man patches) and accepted them.
10. **Phase 3.5 devil's-advocate self-review artifact** is created as a *transient* file (one to two pages) during the Phase 3.5 → 4 transition and **deleted before PR open** (revised 2026-05-25 via plan-acceptance question). The critique-and-patch content survives in the Phase 3.5 commit message body so git history retains the record without cluttering `Author_Reports/`. The Phase 4 human reviewer consults the transient file during their pass; once acceptance is signed off, the file is removed.

The PR is mergeable when all ten hold and a human-acceptance pass has reviewed (and signed off on) the report's prose.

---

## 8. Files to modify / create

| File | Change |
|---|---|
| `Roadmapping/Author_Reports/README.md` | Create — convention for future author reports (naming, build pipeline, expectations) |
| `Roadmapping/Author_Reports/build_report.sh` (or `.py`) | Create — single-command markdown → LaTeX → PDF build script with TODO-stripping + defensive check |
| `Roadmapping/Author_Reports/2026-05_interim_for_gill.md` | Create — report markdown source, §§1–7 per §3 of this plan |
| `Roadmapping/Author_Reports/2026-05_interim_for_gill.tex` | Create (derived) — LaTeX from pandoc, human-edited as needed; *built and committed at end of Phase 4, after human-acceptance* |
| `Roadmapping/Author_Reports/2026-05_interim_for_gill.pdf` | Create (derived) — built PDF; *committed at end of Phase 4, after human-acceptance* (NOT during Phases 1–3, per §6 ordering + §7 acceptance criterion 9) |
| `Roadmapping/Author_Reports/2026-05_interim_for_gill_self_review.md` | Create (transient) — Phase 3.5 devil's-advocate self-review artifact + steel-man patch summary (one to two pages). Lives only during Phase 3.5 → 4 transition for the human reviewer to consult; **deleted before PR open** per §6 phase ordering + §7 acceptance criterion 10. Content survives in git history via the Phase 3.5 commit message |

No source-code or campaign-document changes. The report only *consolidates* existing campaign documents; it does not modify them.

---

## 9. Out-of-scope

- **Re-authoring or correcting any campaign document.** The report references the campaign documents; it does not edit them. If a reviewer identifies an error in a per-PR document while drafting the report, that error is filed as a separate issue, not fixed inline.
- **Building infrastructure for arbitrary report types.** The build script is scoped to *this* report — academic article style, with the specific TODO-stripping discipline. Earlier drafts of this plan justified the script as "convention-setting for future reports"; that framing is downgraded (see §5 scope note) because no second report is planned. A more general report-pipeline (with options for different document classes, journal templates, etc.) is a future thread that will be designed when a second use-case actually appears.
- **Author-side response to the report.** Tepper's answers to Q1–Q8 are recorded in follow-on threads (#54 already opened for Q1; #55 for Q7; further as needed); they are out of scope for this PR.
- **Translation of the campaign verdicts.** The report records what the campaigns said; it does not re-derive their verdicts or change them. Disagreements with the recorded verdicts are filed as separate issues against the relevant campaign.

---

## 10. Honest framing

Per [§13 of #42](42-electromagnetism-jackson-proper-time.md#13-devils-advocate-review-and-what-we-cannot-honestly-fix), inherited verbatim. The report is the campaign's most direct test of the honest-framing discipline: it is the document Tepper will read first, and if it oversells the campaigns' results, the rest of the work is worth less to him.

**The report is solicited** (revised 2026-05-25 after a phone call with Tepper). The earlier draft of this plan framed the report as unsolicited and made apologising-for-imposition a load-bearing discipline; that framing is withdrawn because Tepper has requested the report. The implication: §1's cover note opens directly with what the report is and what we have for him, no peer-note hedge.

The honest-scoping discipline is still load-bearing, but for a different reason. With an unsolicited report, the discipline justified the imposition (we are only sending what we are confident in); with a solicited report, the discipline serves the recipient's ability to act on the material (Tepper needs to know what is load-bearing vs decorative, what is verified vs conditional, what is in scope vs out). The same per-paragraph TODO markers, the same scoping labels, the same `[3, 7]` page bound — they are now justified by *usefulness to the recipient* rather than by *non-imposition on the recipient*.

**Q7/Q8 specifically — removed from the report.** Earlier drafts of this plan listed Q7 as *"framework's intent to develop proper-time one-loop dual-QED."* That phrasing reads as asking Tepper to commit to a multi-year research programme via a checkbox response. The §3 listing now **removes Q7 (and Q8 propagator-sufficiency) from the report entirely** (revised 2026-05-25, in line with the user's "focus more on process and initial results, not next steps" direction). Both remain documented in [#55](https://github.com/temoTxt/PyPhysics/issues/55) for the record. The report's emphasis is what the campaigns *found*, not what they propose for future work.

**Crocco §5 substantive-AI compliance, on the report-as-deliverable.** The report is substantive AI end-to-end. The PDF must not be built and committed before the human-acceptance pass (Phase 4) is complete (see §6's revised phase ordering). The earlier draft of this plan ordered Phase 3 (PDF build) before Phase 4 (human review); that ordering would have stripped per-paragraph `<!-- TODO -->` markers from the deliverable before a human reviewed the draft, materially violating Crocco §5's *"the human-acceptance section is explicitly blank for the human to fill in"* rule. The revised ordering — dry-run the build pipeline in Phase 3, run the devil's-advocate self-review + steel-man in Phase 3.5, then human review + final PDF in Phase 4 — keeps the AI markers visible to the human reviewer and ships the PDF only after acceptance.

The three honest caveats explicit in #50 carry over to the report:

1. **The Lamb-shift route is the Bethe (1947) estimate, not a one-loop QED calculation.** The report's §5 records this explicitly when summarising the Lamb-shift verdict. *Q7 (apparatus development) is no longer surfaced in the report* per the §3 emphasis revision; the caveat lives in §5 as a scoping label, not as a question for Tepper to answer.
2. **`r_e`-dependent observables are conditional on the `r_e` finding's resolution.** The report's §5 records both branches (b) and (c) and prints the six-observable cross-comparison table in full. §6 is the dedicated section asking Tepper to resolve which branch is intended — *not* by best-fit-to-data, but by retrieving the working-notebook digit string or rederiving from the framework's renormalisation prescription. The report's framing of Q1 is decomposed into Q1a/b/c sub-questions to make the resolution actionable.
3. **The campaign's "passes the precision-spectroscopy test" claim is conditional on branch (c).** The report does not assert that the framework agrees with measurement; it asserts that *if* branch (c) is the intended `r_e`, the framework agrees at the Bethe-estimate precision floor across six independent observables. §6 elevates the conditionality from a parenthetical caveat (as it appeared in earlier drafts) to the report's lead question.
4. **Candidate experiments (§4) are scoped, not executed.** The two experimental sub-investigations (#43, #48) are presented as work the campaign has surfaced for Tepper's awareness, *not* as work the campaign has done or is committing to do. The honest framing is explicit: we found these experimental regimes where the framework's predictions are operationally distinguishable; we have not yet executed the comparison.

These caveats are load-bearing in the report. They are not buried in a closing-disclaimer paragraph; they are recorded at the point in §3 / §4 / §5 where the relevant claim appears.

---

## 11. Decision points

Inherited from [§13.5 of #42](42-electromagnetism-jackson-proper-time.md#135-decision-points--confirmed-by-author-2026-05-24) and [§8 of #49](49-quantum-mechanics-griffiths-proper-time.md#8-decision-points):

| Inherited | Override for #58 | Status |
|---|---|---|
| Per-paragraph TODO for substantive content | same — applies end-to-end (the report is substantive end-to-end) | confirmed |
| Voice = Gill's published-paper voice | same — and load-bearing because the audience IS Gill | confirmed |
| Honest scoping discipline | same — every identity/prediction/out-of-scope claim labelled | confirmed |
| Mechanical deferral on flagged findings | does not engage — the report does not modify or re-derive findings | confirmed |
| AI is never an author | same — Tepper's "author response" recipient is human; report's "Author = Trey Morris with Claude Opus 4.7" is honest pragmatic-AI tooling attribution, not journal authorship | confirmed |

One #58-specific decision: **how many cross-references in §2 to cite by URL vs by relative path.** Relative paths assume the reader has the repository checked out; URL paths assume PDF-only reading. Default: footnote-style URL references for the dozen most-cited paths (so the PDF stands alone), relative paths in the markdown source-of-truth (so internal navigation works on GitHub).

---

## 12. Linked PRs / branches

- Parent issue: [#58](https://github.com/temoTxt/PyPhysics/issues/58)
- Branch: `58-author-review-interim-report-for-tepper-gill-goals-progress-and-open-questions-across-the-verification-campaign-work` (current)
- PR (to be opened after Phase 3 build is operational + Phase 4 human-acceptance pass): closes #58

The branch base is `origin/main`. Implementation requires the campaign content (Bethe-Salpeter, Griffiths, Jackson docs) to be present on this branch for the report to reference local files. By the time this PR is being implemented, [#51](https://github.com/temoTxt/PyPhysics/pull/51) / [#52](https://github.com/temoTxt/PyPhysics/pull/52) / [#53](https://github.com/temoTxt/PyPhysics/pull/53) should have merged to main; if not, the implementer rebases this branch onto the most-recent campaign branch (stacked, as #53 stacks on #52).

---

## 13. Working notes

- The build script is scoped to *this* report's format; generalising it to handle multiple document classes / journal templates is a future thread.
- The `Author_Reports/` folder may host multiple reports over time (`2026-05_interim_for_gill.md`, `2026-XX_followup_for_gill.md`, etc.); the `README.md` documents naming convention (`{YYYY}-{MM}_{report_type}_for_{recipient}.{ext}`) so the folder remains navigable.
- The committed PDF is rebuilt from the markdown source-of-truth on every meaningful change to the `.md`. The build script's **reproducibility scope** (not the same as byte-identical determinism — pandoc + pdflatex output is rarely byte-identical across pandoc/texlive versions, font availability, or `/tmp` paths embedded by some LaTeX packages) is: *rebuildable from the same `.md` on a clean checkout without manual intervention*. The build script pins the PDF metadata date via `--metadata date=YYYY-MM-DD` (no timestamp leakage), documents the pandoc + texlive versions it was developed against in its docstring, and accepts that byte-identical output across machines requires a pinned container that is out of scope for this script. Earlier drafts overclaimed byte-identical determinism; the relaxed requirement is the honest one.

---

## 14. Disclosure: this plan is substantive AI

Per Crocco §5 substantive-AI compliance, this plan document — making architectural decisions about the report's structure, build pipeline, phase ordering, honest-framing discipline, and acceptance criteria — is **substantive AI** (Trey Morris with Claude Opus 4.7, 1M context), not pragmatic AI. The decisions encoded here (Phase 3.5 self-review gate, relaxed `[3, 7]` page budget, PDF-after-human-acceptance ordering, voice-match-as-flag-not-assertion, Q7 downgrade, byte-identical-determinism withdrawal, cross-reference no-drops policy, build-pipeline scoped to this report only) are AI-proposed and require a human-acceptance pass on the plan itself before Phase 1 work begins. The plan's substantive-AI status is symmetric with the report's: both require human acceptance before the deliverable ships.

This disclosure was added after a devil's-advocate self-review (analogous to the workflow that produced commits `393efb4` on #51 and `913802d` on #53). The earlier draft of this plan carried no substantive-AI marker; this section corrects that omission. The committed self-review artifact for the *plan* (the critiques + the steel-man patches) is the commit that introduces this disclosure plus the §3 / §5.4 / §5.6 / §6 / §7 / §10 / §13 edits — they are the load-bearing patches that pulled the plan back into alignment with the §13 honest-framing discipline of #42 that the plan claims (and now actually) inherits verbatim.
