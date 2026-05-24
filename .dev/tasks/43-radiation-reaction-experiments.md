# Plan: issue #43 — proper-time predictions for the 2018 radiation-reaction experiments

**Tracks:** [#43 — proper-time predictions for radiation-reaction experiments (Cole 2018, Poder 2018, Wistisen 2018)](https://github.com/temoTxt/PyPhysics/issues/43). Sub-issue of [#42](https://github.com/temoTxt/PyPhysics/issues/42) (Jackson canonical problems × CGS/SI × proper-time Maxwell).
**Status:** plan only; no code or content committed yet. Branch `43-electromagnetism-proper-time-predictions-for-radiation-reaction-experiments-cole-2018-poder-2018-wistisen-2018` is checked out.

## 1. Three findings from the pre-plan triage that materially shape execution

### 1a. **All three papers are CC-BY 4.0 — but PDFs stay local per the resolved decision**

Confirmed via Crossref:

| DOI | License | Authors |
|---|---|---|
| `10.1103/PhysRevX.8.011020` (Cole) | CC-BY 4.0 | 25 |
| `10.1103/PhysRevX.8.031004` (Poder) | CC-BY 4.0 | 24 |
| `10.1038/s41467-018-03165-4` (Wistisen) | CC-BY 4.0 | 4 |

The CC-BY 4.0 licenses *would* allow committing the PDFs publicly. **Per the resolved decision (Trey, 2026-05-24): keep PDFs local-only and commit only the marker-pdf-converted markdown.** This means `pdf_status: acquired` (not `out_of_copyright_public`) — matching the issue's original framing on the *operational* axis (markdown only) even though the underlying licensing axis is more permissive. Consistency over discoverability for this round.

### 1b. **The substantive prediction work is blocked on parent #42 PR D**

#42's plan ([`.dev/tasks/42-electromagnetism-jackson-proper-time.md`](42-electromagnetism-jackson-proper-time.md) §7) sequences a 12-PR campaign:

> PR A (Ch. 6) → PR B (Ch. 11) → PR C (Ch. 12) → **PR D (Ch. 14 Radiation by Moving Charges)** → PR E (Ch. 16 Radiation Damping) → PR F+ (backfill)

PR D is described as "First chapter where the dissipative `(u·a)/b⁴` term *should* contribute non-trivially. Liénard–Wiechert problems are the pivot." Issue #43 explicitly cites PR D as the source of the per-problem derivation template.

**Implication.** Without PR D's per-problem template + the Liénard–Wiechert proper-time third-term derivation in canonical-problem form, the prediction document at `Electromagnetism/Jackson/Experimental_Comparisons/radiation_reaction_2018.md` would need to **re-derive the third term from scratch in its own document** — which duplicates work that PR D will do once, generically, for all Ch. 14 problems. The issue itself acknowledges this: "Depends on #42 PR D ... ; *can be worked in parallel once the three bib stubs land*."

### 1c. **No `Electromagnetism/` directory exists yet**

#42 PR A hasn't shipped either. The campaign's repository scaffolding (`Electromagnetism/Jackson/...`, `Electromagnetism/_template_problem.md`) is not in place. The acceptance criterion "(2) Proper-time prediction document at `Electromagnetism/Jackson/Experimental_Comparisons/radiation_reaction_2018.md`" requires that directory structure first — which is part of PR A's scope.

## 2. What this plan does NOT do

Per §1b: this plan does **not** attempt to re-derive the proper-time Liénard–Wiechert third term in the prediction document. That derivation lives in PR D's Ch. 14 work; #43's prediction document *references* it. If #43 had to do the derivation independently, the work would later need to be reconciled with PR D's canonical form. Better to wait.

Per §1c: this plan does **not** scaffold the `Electromagnetism/` directory tree. That is PR A's scope; doing it here would pre-empt the per-problem template debate that PR A is intended to settle.

Per §1a: this plan does **not** commit the PDFs (only the converted markdown), matching Trey's resolved decision.

Per the issue's own scope guards: this plan does NOT introduce new theoretical predictions in regimes outside 2018 data, re-analyse raw experimental data, or attempt QED corrections to the proper-time framework.

## 3. Staged execution

### Stage 1 — Acquire + convert the three papers (PDFs local-only)

1. Download all three PDFs from publisher (PRX at journals.aps.org; Nat Commun at nature.com). Place under `Roadmapping/Historical_Papers/Retrospective/` — the directory is gitignored by default, so the PDFs stay local without needing to opt them out.
2. Run `uv run python Roadmapping/parse_papers.py --input ./Roadmapping/Historical_Papers/Retrospective --output ./Roadmapping/Historical_Converted_Markdown/Retrospective` to produce per-paper markdown. The `--skip-existing` default protects the already-converted issue-#41 papers from re-conversion.
3. CPU run (no GPU available — same constraint as the issue-#41 workflow because the local `llama-server` holds GPU memory). Wall-clock estimate: 3 papers × ~10–30 min/paper = 30–90 min total. Run in background.

### Stage 2 — Scaffold the three bib stubs

```bash
uv run python Roadmapping/History/Bibliography/scaffold_bib_note.py \
  --cite-key cole2018_radiation_reaction \
  --type retrospective \
  --from-doi 10.1103/PhysRevX.8.011020

uv run python Roadmapping/History/Bibliography/scaffold_bib_note.py \
  --cite-key poder2018_radiation_reaction \
  --type retrospective \
  --from-doi 10.1103/PhysRevX.8.031004

uv run python Roadmapping/History/Bibliography/scaffold_bib_note.py \
  --cite-key wistisen2018_radiation_reaction \
  --type retrospective \
  --from-doi 10.1038/s41467-018-03165-4
```

Then hand-edit each stub to:

- `era: forward` (all three are 2018 experiments, post the "1948–1965" tail of the historical periodization — `forward` matches the convention).
- `tags`: `['#era/forward', '#thread/electromagnetism', '#gill-silent', 'radiation-reaction', 'strong-field-qed']` plus per-paper: Cole/Poder add `'laser-electron-collision'`, Wistisen adds `'crystal-channeling'`.
- `pdf_status: acquired` per finding §1a.
- `pdf_path` to the local (gitignored) location for traceability.
- `converted_md` to the committed markdown path.
- `gill_corpus_overlap: ["Two_Mathematically_Equivalent_Versions_of_Maxwells_Equations"]` since the proper-time framework that's being tested against these experiments is verified in that doc.
- `chapters_citing: []` — no chapter cites them yet; #43's prediction document will be the first, but it does not exist as a chapter.
- `human_reviewed: false` — flip only when Trey has personally read each paper.

Regenerate `Historical_Papers/Acquisition_Tracker.md` via `update_acquisition_tracker.py`.

### Stage 3 — Defer the prediction document until PR D lands

Issue #43's acceptance criterion #2 (the `radiation_reaction_2018.md` prediction document) requires:

(a) The `Electromagnetism/Jackson/Experimental_Comparisons/` directory scaffolded (PR A's scope).
(b) The per-problem template at `Electromagnetism/_template_problem.md` (PR A's scope).
(c) The proper-time Liénard–Wiechert third term re-derived in canonical-problem form (PR D's scope).

Without all three, the prediction document either duplicates work or has to be rewritten when PR D lands. **Recommendation:** comment on #43 noting that the prediction document is staged for after PR D, and explicitly time-order it: PR A → PR D → #43 prediction doc. This is consistent with the issue's own "can be worked in parallel once the three bib stubs land" phrasing — meaning: bib stubs now, prediction document later.

### Stage 4 — Crocco compliance pre-work (deferred to Stage 3 timeframe)

Per #42's compliance §6 + this issue's AC#4:

- **Pragmatic** (no separate disclosure section): Mathematica MCP for symbolic verification of the third-term contribution to the energy-spectrum shift; running the standard Landau–Lifshitz comparison; formatting the per-problem template; translating between CGS and SI as needed.
- **Substantive** (`<!-- TODO: human reviews and fills in -->` block required): any prose that proposes a *physical interpretation* of why the proper-time prediction fits or misses; any narrative framing ("this is the headline payoff", "this dissolves the Abraham–Lorentz pathology"); the actual verdict (✅ / ⚠ / ❌) when it requires judgement about which classical comparator is appropriate; the prediction document's top-of-document **Selection provenance** note.

The committed prompt-of-record for any substantive AI use in Stage 3 would land under `Roadmapping/Tooling/_prompts/` before the prediction document commits. Likely candidate: a new `derive_proper_time_prediction.md` prompt, or reuse of `chapter_qa_review.md` for the final QA pass. **Defer the prompt-template choice to when Stage 3 is unblocked.**

### Stage 5 — FINDINGS cross-post (conditional, deferred)

If Stage 3 produces a ⚠ or ❌ verdict on any experiment, add a section to [`Roadmapping/Equation_Verification/FINDINGS_for_author_review.md`](../../Roadmapping/Equation_Verification/FINDINGS_for_author_review.md) following the existing pattern (currently 3 findings: Maxwell Eq. 24, DRQM I Eq. III.22, TCEP Eq. 4.16). The new finding(s) would be structurally different — they're predictions-vs-experiment comparisons rather than internal-derivation errors — so the FINDINGS format may need a sub-template extension. Worth a small section in PR D's discussion to settle the format before the first prediction-finding goes in.

## 4. What this PR ships (if §3 Stages 1–2 are executed)

Realistic scope **today**:

- [ ] Three PDFs in `Roadmapping/Historical_Papers/Retrospective/` (NOT committed — kept local under the gitignore default).
- [ ] Three converted markdown trees in `Roadmapping/Historical_Converted_Markdown/Retrospective/` (committed).
- [ ] Three hand-edited bib stubs in `Roadmapping/History/Bibliography/Retrospective/` (committed).
- [ ] Updated `Roadmapping/Historical_Papers/Acquisition_Tracker.md` (committed).
- [ ] This plan doc committed.
- [ ] A comment on issue #43 explaining the deferral of the prediction document until PR D, and that the bib-stub work is complete.

**Issue #43 stays open** until the prediction document lands. The PR closes only the bib-stub portion of AC#1; AC#2, #3, #4 remain pending PR D.

## 5. Alternative: pivot to PR A

If the intent of "work on issue #43" is really "make progress on the radiation-reaction comparison," the highest-leverage move may actually be **PR A of #42** (Ch. 6, the campaign's foundation chapter), since that unblocks PR D, which in turn unblocks #43's prediction document. PR A has its own 5 problems (Maxwell equations with magnetic monopoles, EM momentum of point charge, etc.) — different work, but it's the path to #43's headline payoff.

Worth raising as an option but not recommended without explicit user direction, since PR A is a substantial effort (5 problems × full three-way template, plus the campaign-foundation scaffolding).

## 6. Resolved decisions

1. ✅ **PDF acquisition tier:** keep PDFs local, commit only converted markdown. `pdf_status: acquired` on all three stubs. (Resolved Trey, 2026-05-24.)
2. ✅ **Scope of this PR:** override default recommendation — **attempt the prediction document anyway and accept the duplication** that PR D will eventually reconcile. (Resolved Trey, 2026-05-24.) The prediction document at `Electromagnetism/Jackson/Experimental_Comparisons/radiation_reaction_2018.md` ships in this PR. The directory tree it creates (`Electromagnetism/Jackson/Experimental_Comparisons/`) is minimal and will be subsumed by PR A's full scaffolding when that lands. Crocco-substantive sections (verdicts, physical interpretations) carry the mandatory `<!-- TODO: human reviews and fills in -->` blocks rather than being filled in by the AI.
3. ✅ **No pivot to PR A** — work the prediction document on the current branch per Q2 above.
4. ⏸ **Prediction-finding FINDINGS format (Stage 5):** defer to PR D discussion. *No action required from Trey now.* If this PR's prediction document produces a ⚠ or ❌ verdict (filled in by Trey post-PR), the FINDINGS cross-post is a separate small PR after.
