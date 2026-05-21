# Crocco et al. (2026) compliance — design commitments

This repository's tool suite — the bibliography pipeline, the crawl + triage queue, the synth tools, the slash commands — is designed to satisfy the eight reporting expectations articulated by:

> Crocco, O. S., Rasdi, R. M., & Garavan, T. N. (2026). Responsible AI in Non-Empirical Research: Guidance for Authors, Reviewers, and Editors in HRD. *Human Resource Development Review*. https://doi.org/10.1177/15344843261445531

This document records (1) a short summary of the article's framework, (2) the explicit mapping from each expectation to a concrete commitment in this repo, and (3) a verification checklist a future reviewer (Crocco himself or anyone else) can run.

## 1. Why this document exists

Issue [#34](https://github.com/temoTxt/PyPhysics/issues/34) tagged **@olivercrocco** (the article's first author) as the human-in-the-loop reviewer of the dissertation-tooling survey. Crocco's article is the canonical statement of how AI should be used in *non-empirical* scholarship — exactly the genre this repo produces (literature reviews, conceptual chapters, theory-comparison work, bibliometric analysis). Aligning our tooling to his expectations *before* substantive AI use happens means:

- His review of our plan operates against a target he already authored.
- We avoid retroactive disclosure churn (rewriting reports after-the-fact to add provenance).
- The repo doubles as a *worked example* of Crocco's framework applied in a different field (physics-history rather than HRD).

## 2. The article's framework, condensed

### Six AI touchpoints

Crocco identifies six points in non-empirical research where AI may be invoked:

(a) **Literature search and retrieval** — finding papers.
(b) **Screening and selection** — triaging which to include.
(c) **Data extraction and coding from existing literature** — pulling structured facts out of papers.
(d) **Synthesis and thematic identification** — finding patterns across papers.
(e) **Conceptual development and framework building** — proposing new analytical structures.
(f) **Writing, structuring, and refining manuscripts** — composing the final text.

### Pragmatic vs. substantive distinction (load-bearing)

- **Pragmatic use**: AI as efficiency tool that does *not* shape the intellectual contribution. Examples: grammar, formatting, translation, reference management, bibliometric visualisation, post-hoc claim-citation auditing.
- **Substantive use**: AI in activities that shape *what the manuscript argues, synthesizes, or proposes*. Examples: generating thematic categories, proposing conceptual frameworks, producing draft sections that convey core arguments.

The distinction is **central to disclosure** — substantive uses require fuller transparency than pragmatic ones.

### Three-tier typology

- **Legitimate**: researcher's intellectual contribution remains entirely their own. Efficiency, grammar, translation, ref-mgmt, claim-source auditing, bibliometric viz.
- **Contested**: middle ground; defensible if the researcher demonstrates ownership. AI-assisted thematic coding the researcher verifies, AI-suggested frameworks substantially developed, AI-generated drafts fundamentally rewritten.
- **Illegitimate**: submitting AI-generated literature reviews as original synthesis; fabricating or failing to verify AI-suggested citations; presenting AI conceptual models without substantial development; producing manuscript sections via AI without disclosure.

## 3. The eight reporting expectations — and our commitments

| # | Crocco's expectation | Our commitment in this repo |
|---|---|---|
| **1** | Disclose all AI tools used and their specific role in each phase | The 9 chapters' commits cite Claude Opus 4.7 via `Co-Authored-By` trailer. The plan doc [`.dev/tasks/34-dissertation-tooling.md`](../../.dev/tasks/34-dissertation-tooling.md) records which AI is used where (Phase 3 crawl = no AI, Phase 4 synth = substantive Claude calls, Phase 5 triage = pragmatic Claude calls). Every synth report under `synth_reports/` (Phase 4) will include a per-section `ai_use_class` tag. |
| **2** | Specify whether AI was used for pragmatic support or substantive contribution | Every prompt template in [`_prompts/`](_prompts/) declares `ai_use_class: pragmatic\|substantive` in its frontmatter. Synth reports inherit this from the prompt they invoked. |
| **3** | Describe verification procedures for AI-generated content, particularly citations and source accuracy | The bibliography pipeline's three-step verification chain: (a) crawl candidates come pre-tagged with verified DOIs from Crossref / S2 / ArXiv; (b) `sync_from_zotero.py` flags drift on existing entries rather than overwriting; (c) the YAML field `human_reviewed` (per expectation 7) confirms the human read the source. The `validate_wikilinks.py` tool catches broken citations before commit. |
| **4** | Clarify how the argument was developed | Synth reports (Phase 4) have a mandatory "human acceptance" section blank-filled by the human after reading the AI's suggestions — recording which were accepted, which rejected, and why. Chapter PRs that act on synth-report suggestions cite the suggestion + the human's reasoning in the PR description. |
| **5** | Address potential biases introduced by AI tools | [`crawl/BIASES.md`](crawl/BIASES.md) documents each crawler's known systematic skews (S2 English-leaning, ArXiv discipline-restricted, Crossref publication-channel-biased, Zotero curation-biased) and mitigations. The bias section is referenced from any synth report whose inputs depend on crawler-surfaced candidates. |
| **6** | Provide prompts or model specifications when AI played a substantive role | Every substantive prompt is committed under [`_prompts/`](_prompts/) with its `ai_use_class` frontmatter, intended invoker, and recommended model. Synth reports cite the prompt's path + git commit SHA so the *exact* prompt that ran can be retrieved. |
| **7** | Confirm engagement with primary sources | The bib YAML schema includes `human_reviewed: bool`. New stubs default to `false`. The audit tool [`audit_human_reviewed.py`](../History/Bibliography/audit_human_reviewed.py) scans the bib + sets the field via a conservative body-word-count heuristic. Stubs flagged `false` should be flipped manually only after the human reads the source and writes a 2-3-paragraph body summary. |
| **8** | Take full responsibility for all content (AI cannot be listed as an author) | Every commit's `Co-Authored-By: Claude Opus 4.7 (1M context) <noreply@anthropic.com>` trailer lists Claude as a *tool*, not an author. The PR author of record on every chapter is Trey. AI is not credited as an author of any committed artefact. |

## 4. Touchpoint coverage in this repo

| Crocco touchpoint | Tool / pipeline | Use class |
|---|---|---|
| (a) Literature search + retrieval | `crawl/from_s2.py`, `crawl/from_arxiv.py`, `crawl/from_crossref.py`, `crawl/from_zotero.py` | **No AI invoked** — pure API queries against verified citation graphs |
| (b) Screening + selection | Phase 5 `/triage-papers` slash command | pragmatic (Claude summarises abstracts; human decides) |
| (c) Data extraction + coding | (Not yet implemented; would land in a future PR) | n/a |
| (d) Synthesis + thematic identification | Phase 4 `cluster_claims.py` (via `_prompts/synth_cluster_claims.md`) | **substantive** |
| (e) Conceptual development + framework building | Manual; Claude Code chats during chapter authoring | varies (declared per chat) |
| (f) Writing, structuring, refining manuscripts | Manual; Claude Code chats; `_prompts/chapter_qa_review.md` for QA pass | pragmatic (QA); substantive (drafting) — declared per use |

## 5. Verification checklist (for a future reviewer)

A reviewer (Crocco or anyone else) can audit this repo's compliance by running:

1. **Tool disclosure (expectation 1).** Read the most recent 5 chapter PRs' descriptions; confirm each cites the AI tools used in writing.
   ```bash
   gh pr list --state merged --limit 5 --json title,body,number
   ```

2. **Pragmatic vs substantive classification (expectation 2).** List every prompt template and confirm each has `ai_use_class`:
   ```bash
   for f in Roadmapping/Tooling/_prompts/*.md; do head -5 "$f" | grep ai_use_class; done
   ```

3. **Citation verification (expectation 3).** Run the wikilink validator + the human_reviewed audit:
   ```bash
   uv run python Roadmapping/History/_tools/validate_wikilinks.py
   uv run python Roadmapping/History/Bibliography/audit_human_reviewed.py
   ```

4. **Argument-development transparency (expectation 4).** Inspect any committed synth report under `synth_reports/` (Phase 4 will populate). Each should have non-empty "human acceptance" section.

5. **Bias documentation (expectation 5).** Confirm `crawl/BIASES.md` exists and lists every active crawler.

6. **Prompt-template availability (expectation 6).** For each synth report citing a prompt, confirm the prompt exists in `_prompts/` and the cited commit SHA resolves.

7. **Primary-source engagement (expectation 7).** Count bib stubs with `human_reviewed: true` vs `false`. Sample a few of each and confirm the body summary matches the flag.
   ```bash
   uv run python Roadmapping/History/Bibliography/audit_human_reviewed.py --quiet
   ```

8. **No AI as author (expectation 8).** Confirm no `Author:` field of any commit lists Claude / GPT / Gemini / etc. — only `Co-Authored-By:` trailers.
   ```bash
   git log --all --format='%an <%ae>' | grep -i 'claude\|chatgpt\|gpt\|gemini' || echo "no AI authors found"
   ```

## 6. Items deferred

Things Crocco's framework suggests we *could* do but haven't yet committed to:

- **Per-section provenance in chapters.** Crocco hints at section-by-section disclosure ("which parts of this chapter were drafted vs. revised by AI?"). We currently disclose at the PR level. Section-level disclosure is more granular; could be added if any chapter has an unusually mixed AI/human composition.
- **Bias mitigations beyond crawler-level.** Crocco's expectation 5 also covers *synthesis bias* — that AI's choice of which themes to cluster reflects its training data. Our Phase 4 prompts (`synth_cluster_claims.md`) include "quote, don't paraphrase" rules to make biased inferences traceable, but a dedicated synthesis-bias audit tool would go further.
- **Reproducibility prompts as code-block-only artefacts.** Crocco notes prompts should be reproducible. We commit them in markdown with prose framing — a stricter discipline would commit them as plain text with no surrounding commentary so they can be programmatically piped to a model. Open question whether the simpler form helps or hurts review.

These are open for follow-up if Crocco's review surfaces a specific need.

## 7. Provenance of this document

This compliance doc and its enabling tools (the YAML field, the audit tool, the `_prompts/` directory, the BIASES doc) were drafted by Claude Opus 4.7 in conversation with Trey Morris, based on a reading of Crocco et al. 2026 attached to issue #34. The reading was performed by the AI; the design choices encoded here were proposed by the AI in conversation and accepted by the human author. Per Crocco's expectation 8, Trey holds final responsibility for what's in this repo.

The article's prompt for this compliance doc maps onto Crocco's own disclosure section: items (1), (4), and (5) of his disclosure ("searching for and cross-referencing", "auditing references", "copyediting") are pragmatic uses analogous to how `_prompts/chapter_qa_review.md` is intended; items (2) and (3) ("developing and refining the manuscript's structural outline", "generating preliminary draft text") are substantive uses analogous to how `_prompts/synth_cluster_claims.md` will be used in Phase 4. The structural symmetry is intentional — Crocco's own paper is the worked example we're patterning on.
