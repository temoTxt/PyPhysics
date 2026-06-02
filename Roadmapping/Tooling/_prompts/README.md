# Prompt templates

Implements Crocco et al. (2026) reporting expectation #6 — **"Provide prompts or model specifications when AI played a substantive role"** — by committing every substantive-AI prompt as a version-controlled artefact in this directory.

## What lives here

Each `.md` file is a single prompt template plus its compliance metadata. The header frontmatter encodes:

- **`ai_use_class`** — `pragmatic` or `substantive` (Crocco's central distinction; see the [compliance doc](../CROCCO_COMPLIANCE.md)).
- **`touchpoint`** — one of Crocco's six touchpoints: `literature_search`, `screening`, `data_extraction`, `synthesis`, `conceptual_dev`, `writing`.
- **`invoked_by`** — what tool / slash command consumes this prompt.
- **`model`** — recommended Claude model (or generic family).
- **`status`** — `draft` (not yet used in committed output) or `active` (used by current tooling).

The body is the prompt itself.

## Why version-controlled?

If the prompt that produced a clustered claim or a suggested cross-reference changes, the reproducibility of the synthesis output changes. Crocco emphasises: *"AI-assisted analyses are interpretable only when readers know what corpus was used, how topics were defined and validated, and what prompts generated the output"* (§"Reporting Expectations").

When a synth report under [`Roadmapping/Tooling/synth_reports/`](../synth_reports/) cites a prompt from this directory, it includes the **commit SHA** at the time of use — so anyone reviewing the report can `git show <SHA>:Roadmapping/Tooling/_prompts/<file>.md` and see the exact prompt that ran.

## Current templates

| File | Touchpoint | Class | Used by | Status |
|---|---|---|---|---|
| [`triage_summary.md`](triage_summary.md) | screening | pragmatic | `/triage-papers` (Phase 5) | draft |
| [`synth_cluster_claims.md`](synth_cluster_claims.md) | synthesis | **substantive** | `cluster_claims.py` (Phase 4) | draft |
| [`synth_suggest_cross_refs.md`](synth_suggest_cross_refs.md) | synthesis | **substantive** | `suggest_cross_refs.py` (Phase 4) | draft |
| [`chapter_qa_review.md`](chapter_qa_review.md) | writing | pragmatic | ad-hoc Claude Code chats | draft |
| [`evaluate_external_suggestion.md`](evaluate_external_suggestion.md) | synthesis | **substantive** | ad-hoc Claude Code chats responding to collaborator-forwarded papers (issues, discussions) | active |

## Adding a new prompt

1. Decide whether the use is pragmatic or substantive (see the compliance doc for the rubric).
2. Create `Roadmapping/Tooling/_prompts/<short_name>.md` with the frontmatter above + the prompt body.
3. Update the table in this README.
4. When a tool starts invoking the prompt, change `status: draft` to `status: active`.
5. **Substantive prompts (only)**: subsequent edits to the body must be traceable in git history — never amend; always commit a new revision.
