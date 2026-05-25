# Phase 3.5 devil's-advocate self-review — `2026-05_interim_for_gill.md`

**Transient artifact.** This file is created in Phase 3.5 for the Phase 4 human reviewer to consult. It is **deleted before the PR opens** per the campaign plan; the critique-and-patch content survives in the Phase 3.5 commit message body.

Devil's-advocate pass against the §13 honest-framing discipline of [#42](https://github.com/temoTxt/PyPhysics/issues/42), applied to the AI-drafted report. Pattern follows the steel-man passes on PRs #51 and #53 (commits `393efb4`, `913802d`).

## Critiques caught and steel-man patches applied

### Critique 1 — §5 oversold the Lamb-shift result as "the framework's strongest independent endorsement"

**Original wording:** *"The Bethe-estimate Lamb shift (above) is the framework's strongest independent endorsement: it does not depend on `r_e`, because the leading log-Bethe self-energy contribution is `g_s = 2`-symmetric…"*

**Problem.** "Endorsement" implies the framework was *distinguished from standard QED* by this result. It was not. The Bethe-estimate route in the dual-theory formulation gives the same `~1,016` MHz number as the standard QED Bethe-estimate route. The `~42` MHz gap to measurement is shared between both formulations — it is what neither captures without a full one-loop calculation. Claiming an "endorsement" overstates the result by ignoring that it does not differentiate the framework from textbook QED at the precision the apparatus can presently access.

**Patch applied.** Rewrote the paragraph to characterise the Lamb-shift result as the campaign's *clearest measurement-level non-falsification* at the Bethe-estimate precision floor — and explicitly noted that the same number is what standard QED's Bethe-estimate route gives. Added a sentence that the framework's apparatus cannot push below this precision floor without the one-loop dual-QED work documented in [#55](https://github.com/temoTxt/PyPhysics/issues/55). The point that *the framework does not fail* at Bethe-estimate precision is preserved (it is the strongest empirical claim the apparatus supports); the overclaim that this constitutes an *endorsement* is withdrawn.

**Why this matters for the recipient.** Tepper is a senior physicist who will read "strongest independent endorsement" and immediately ask: endorsement *relative to what alternative?* If the answer is "the same Bethe estimate that standard QED gives," then the language is positioning rather than reporting. The patched language reports.

### Critique 2 — §6.3 Q1a "the framework owes an explanation" read as accusatory

**Original wording:** *"…the campaign's verdicts on the six observables in §5 should be recorded as 🔴 (ruled out at current measurement precision), and the framework owes an explanation for the ~10⁻³ disagreement across hyperfine, fine structure, and `g`."*

**Problem.** "Owes an explanation" carries a faintly accusatory register that does not belong in a peer note to a senior co-author about a finding that is plausibly a transcription error. The branched verdict structure is honest, but the framing of the consequence should be honest *and* diplomatic: a `~10⁻³` disagreement is a real consequence, but the appropriate framing is *"this would prompt a question about reconciliation"* rather than *"the framework owes an explanation."*

**Patch applied.** Rewrote the conditional consequence as *"the open question becomes how the framework reconciles the ~10⁻³ disagreement across hyperfine, fine structure, and `g` — either through a feature of the derivation we have not captured, or through a regime distinction we have missed."* This makes the consequence concrete (the question that would arise) without imputing fault to the framework or the author. Also softened "should be recorded as 🔴" to "would be recorded as 🔴" (matches the conditional tone).

**Why this matters for the recipient.** Tepper is also being asked, in the same paragraph, to confirm that the published `r_e` value *was* a transcription error (Q1b). If Q1a's consequence is framed as "the framework owes an explanation," the report reads as putting him in a corner: either confess to the transcription error or accept the accusatory framing. The patched language is diplomatic; the substance is the same.

## Critiques considered but not patched

### Critique 3 — §3 J3e-P14.2 "longitudinal radiation component"

**Possible issue.** The Maxwell paper Eq. (7) third term `e(u·a)[r×(u×r)]/(b⁴s³)` contains the vector `r × (u × r) = r² u - r(u·r)`, which is the projection of `u` perpendicular to `r`. In the standard radiation-zone picture, the propagation direction is `n = r/|r|`, and *transverse* means perpendicular to `n`. The third term is perpendicular to `n`, so by the standard categorisation it is *transverse*, not longitudinal.

**Reasoning for not patching here.** The Jackson Ch. 14 campaign document (PR #51) claims the third term predicts a longitudinal component, and the report should be consistent with the campaign work it summarises. If the campaign work itself has misframed the result, that is a separate issue (file an issue against the campaign document, not against the report). The report's TODO block on this paragraph already flags the question explicitly: *"confirm the longitudinal-component framing is accurate; flag any alternative interpretation that should be foregrounded."*

**Recommendation for Phase 4 human reviewer.** Before sending, verify the longitudinal-vs-transverse claim against the campaign's actual derivation in `Roadmapping/Electromagnetism/Jackson/Ch14_Radiation_by_Moving_Charges.md`. If the campaign's claim is "transverse but with a different angular distribution than classical LW," the report's §3 sentence should be rewritten to match. If the campaign's claim is genuinely "longitudinal in some sense the report is correctly summarising," leave as is.

### Critique 4 — §3 J3e-P5.13 included despite being unsettled

**Possible issue.** The Zorn surprise finding has not been confirmed by a second pair of eyes (the @Zorns-Lemmon tag on issue #42 is awaiting review). The report promotes the result to a "Jackson highlight" while flagging the un-settled-ness in the TODO. Some readers will argue that the report should not promote un-settled findings to a co-author summary.

**Reasoning for not patching here.** The result is presented honestly: "we have not yet had a second pair of eyes on the derivation." The reason for including it in the report is exactly that — getting a senior co-author's view is one valid way to settle it. Demoting it would lose the opportunity. The Phase 4 human reviewer can decide whether to drop it before sending.

**Recommendation for Phase 4 human reviewer.** If Tepper is the right reviewer for the J3e-P5.13 result, keep it as-is. If a Zorn review is expected before the report is sent and would replace the need for Tepper's view, demote J3e-P5.13 to a footnote in §3 or move it to §7 as Q6.

### Critique 5 — Voice register

**Possible issue.** The plan §4 required *neutral co-author-to-co-author* register for the AI draft (override of campaign-wide Gill voice). The draft uses "we" pervasively and connective openers occasionally, which both read like Gill's voice. Has the override actually been applied?

**Reasoning for not patching here.** Pluralisation ("we did X") is standard academic prose and is not exclusive to Gill's voice; ditto connective openers. The draft does *not* use the load-bearing Gill-voice fingerprints from `VOICE_MATCH_GILL.md` — no "mathematically equivalent but not physically equivalent" phrases, no "reformative framing" stylistic moves, no narrative attribution. The register reads as neutral academic prose. If the Phase 4 human reviewer disagrees, individual paragraphs can be reworded; the global register seems correct.

**Recommendation for Phase 4 human reviewer.** Read §1 (cover note) aloud and judge whether it reads as a peer note from one co-author to another. If it reads as performatively Gill-voiced, soften the connective openers in §1; otherwise leave the voice as-is across the whole report.

## Verification of build pipeline

The Phase 3 dry-run after these patches were applied — verified in the Phase 3.5 commit's terminal output — reports 5 pages, within the `[3, 7]` budget, with all defensive checks passing (no `<!-- TODO` markers leaked to LaTeX; no `human reviews and fills in` text leaked; page-count check pass). The patched draft is ready for the Phase 4 human-acceptance pass.

## Notes for Phase 4

1. The two load-bearing patches (Critiques 1 + 2) have been applied. Read the patched paragraphs in §5 and §6.3 specifically; if you disagree with the rephrasing, revert and pick a different wording.

2. Critique 3 (J3e-P14.2 transverse-vs-longitudinal) is *not* patched because the report should be consistent with the campaign document. Resolve the question in the campaign document first if needed.

3. Critique 4 (J3e-P5.13 un-settled) is *not* patched because the result is presented honestly. Drop it from the report if you would prefer to settle it via Zorn first.

4. Critique 5 (voice) is judged to be already correct (neutral register, not Gill-voiced). Confirm or override per your read.

5. The per-paragraph `<!-- TODO -->` markers in the markdown source are the load-bearing review prompts. Resolve each one before building the final PDF in Phase 4.

6. **The final PDF is built and committed at the end of Phase 4**, not before. After you resolve the TODO markers, run `./build_report.sh 2026-05_interim_for_gill` (without `--dry-run`) to produce `2026-05_interim_for_gill.tex` and `2026-05_interim_for_gill.pdf` in this folder. Then **delete this `_self_review.md` file** before opening the PR (the critique-and-patch content survives in the Phase 3.5 commit message).
