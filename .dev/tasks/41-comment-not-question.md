# Plan: respond to issue #41 (comment, not question)

**Tracks:** [#41 — "Something for your (or Claude's) consideration"](https://github.com/temoTxt/PyPhysics/issues/41)
**Status:** plan only; no code or content committed yet. No GitHub response posted yet.

## 1. What was actually shared

The thread is a *forwarded lead*, not a specific question — Zorns-Lemmon explicitly disclaims having read either paper ("I can't say whether it is irrelevant, mundane, brilliant, or crackpot"). Hence the task-file name: **comment, not question**. The job is to triage two unrelated papers and produce a substantive reply, not to "answer" anything narrowly.

| # | Source | What | Where |
|---|---|---|---|
| A | Issue body | "Quantum-Mechanical Physics as Invariant Geometric Structure" | [zenodo.org/records/20109773](https://zenodo.org/records/20109773) + PDF attached to issue |
| B | First comment (Zorns-Lemmon) | Bootstrap → string-theory amplitudes (Veneziano, Virasoro-Shapiro) from ultrasoft Regge zeros | [arXiv:2508.09246](https://arxiv.org/abs/2508.09246) + PDF attached |
| — | Second comment (Trey) | Already acknowledged; flagged that Discussions tab is now enabled | — |

Both PDFs are GitHub-attachment-hosted and not yet in the repo. Verified by `find`: no local copy of either.

## 2. Constraints this plan must honor

- **Crocco substantive vs pragmatic.** Any *evaluation* of either paper that gets published (in a GitHub reply, in chapter prose, or in the bibliography body) is **substantive**. The prompt-of-record must live under [`Roadmapping/Tooling/_prompts/`](Roadmapping/Tooling/_prompts/) and the output must include a `<!-- TODO: human reviews and fills in -->` acceptance section per CLAUDE.md rule 5.
- **No invented citations.** If the reply or any synth report references either paper by cite-key, scaffold the bib stub first via [`scaffold_bib_note.py`](Roadmapping/History/Bibliography/scaffold_bib_note.py). Both have DOIs/arXiv IDs, so `--from-doi` should work for B; A's Zenodo record needs the DOI extracted from the page (10.5281/zenodo.20109773 if that's the canonical form — confirm before scaffolding).
- **`human_reviewed: false`** on both bib stubs until Trey reads the primary sources.
- **AI is never an author** of the eventual GitHub reply. Claude drafts; Trey edits and posts.
- **Don't auto-post to GitHub.** The final reply lands when Trey decides, not when Claude finishes drafting.

## 3. Staged execution (one PR per stage, or one bundled PR — see §6)

### Stage 1 — Acquire and convert both PDFs

1. Download the two GitHub-attached PDFs into a local working directory. Zenodo paper looks like a preprint (likely public-domain or CC-licensed — confirm before force-adding); arXiv preprint is public by default.
2. Per the PDF acquisition policy: if both license-clear, place under `Historical_Papers/Retrospective/` and `git add -f`. Otherwise keep PDFs local, commit only the converted markdown.
3. Run `uv run python Roadmapping/parse_papers.py --input <dir> --output Roadmapping/Historical_Converted_Markdown/Retrospective/`.

### Stage 2 — Scaffold bib stubs

```bash
# Paper B — has arXiv ID, easy
uv run python Roadmapping/History/Bibliography/scaffold_bib_note.py \
    --cite-key <author>2025_bootstrap_string_amplitudes \
    --from-doi 10.48550/arXiv.2508.09246

# Paper A — Zenodo, may need manual metadata if Crossref lookup fails
uv run python Roadmapping/History/Bibliography/scaffold_bib_note.py \
    --cite-key <author>2025_invariant_geometric_qm
# then hand-edit YAML
```

Both go under `Roadmapping/History/Bibliography/Retrospective/` (third-party works, not Gill corpus). Tag both with `#thread/external-suggestion` (new tag — adopt or rename) so future triage can find collaborator-forwarded papers. Set `chapters_citing: []` until a chapter actually cites them. Regenerate the acquisition tracker with `update_acquisition_tracker.py`.

### Stage 3 — Substantive evaluation per paper

Write `Roadmapping/Tooling/_prompts/evaluate_external_suggestion.md` (new) — a versioned prompt template that asks for:

1. One-paragraph plain-English summary of the paper's central claim.
2. Verdict on Zorns-Lemmon's four-way axis: *irrelevant / mundane / brilliant / crackpot* — with an honest "uncertain" allowed.
3. Specific intersections with the Gill program:
   - **Paper A:** does an "invariant geometric structure" formulation parallel Gill's dual-proper-time geometry (the `b`/`c` modifications)? Cite verification-doc anchors if so.
   - **Paper B:** is the bootstrap method (unitarity + locality + Lorentz invariance + finite-energy constraint → unique amplitude) viable as a derivational path for Gill's `b/c` factors? This is the *methodological* question Zorns-Lemmon actually raised; treat it as the central question of the reply.
4. Any open questions worth opening as follow-up issues or discussion threads (matches Trey's "branch threads" remark in comment 2).

Output goes to `Roadmapping/Tooling/synth_reports/2026-05-23_issue41_external_suggestions.md` with the mandatory human-acceptance section blank.

### Stage 4 — Draft the GitHub reply

A single comment on issue #41 (or a Discussion thread, see §5) that:

- Thanks Zorns-Lemmon and confirms both papers were read, not skimmed.
- Gives a one-paragraph verdict on each, framed in Zorns-Lemmon's own four-way axis.
- Engages substantively with the bootstrap-method suggestion: is it a plausible derivational route for `b/c`, or do the Gill modifications fail one of the bootstrap input axioms (e.g., is the dual-mass geometry compatible with standard locality)?
- Links any spun-off discussion threads / sub-issues by URL.
- **Does not** auto-merge either paper into chapter prose or claim verification anchors that don't yet exist.

Draft lives in the synth report; Trey copies-edits-posts.

### Stage 5 — Wire spin-offs (only if §3 surfaces them)

If the bootstrap-method suggestion looks viable, open a new issue: *"Investigate bootstrap-style derivation of Gill `b/c` modifications from {unitarity, locality, Lorentz invariance, finite energy}"*. Label `research-direction`. Cross-reference from a Forward chapter (7–9) under a `#speculative` claim — never `#established` on the strength of one paper.

If the Zenodo paper's geometric framing connects to Foundations II Sec 2.2.1 (the `r_0` critical-point thread that already resolved Maxwell-paper open-Q3), add a *retrospective* note to `FoundationsII-Classical.md` cross-referencing the bib stub — but only after Trey reads the Zenodo paper personally and flips `human_reviewed: true`.

## 4. What this plan deliberately is NOT

- **Not a verdict.** This plan does not pre-judge either paper. The verdict is produced in Stage 3 against a committed prompt; pre-judging here would smuggle a substantive AI claim into the planning doc.
- **Not a chapter-prose change.** Nothing in Stages 1–4 edits a History chapter, an Equation_Verification doc, or a Manim scene. Spin-off chapter edits, if any, are Stage 5's territory and require a separate PR.
- **Not an automatic Discussions migration.** Trey opened the Discussions tab; this plan does not move the thread on its own.

## 5. Decisions (resolved 2026-05-23)

1. **Reply venue:** comment on issue #41 (not Discussions). Stage 4 draft targets this issue.
2. **License check on Paper A (Zenodo):** skipped this round. Treat both PDFs as `pdf_status: acquired` and keep them local — commit only the marker-pdf-converted markdown under `Historical_Converted_Markdown/Retrospective/`. Do not `git add -f` the PDFs. Revisit license posture later if either paper is promoted into the corpus proper.
3. **Tag adoption:** tag forwarded papers with `#contributor/zorns-lemmon` (new convention — credits the human who surfaced the lead, distinct from `#era/...` and `#thread/...` axes). Future collaborator-forwarded papers get `#contributor/<github-handle>`.
4. **Scope of the reply:** engage the bootstrap-as-derivational-tool methodological angle in the Stage 4 reply. Treat it as the substantive core of the response — but stay bounded ("here's why it's worth/not worth pursuing, candidate follow-up issue linked" rather than a multi-paragraph speculation). Per-paper verdicts come first; methodological engagement second.

## 6. PR shape

Default: **one bundled PR** containing Stages 1–3 (PDFs in / converted / bib stubs scaffolded / synth report drafted) plus the draft reply text in the synth report. Stage 4 is Trey posting; Stage 5 is conditional and earns its own PR if triggered. Reasoning: Stages 1–3 are tightly coupled and meaningless in isolation, and the synth report carries the human-acceptance gate that already separates Claude's draft from Trey's published reply.
