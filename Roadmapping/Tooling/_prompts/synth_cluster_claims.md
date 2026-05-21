---
ai_use_class: substantive
touchpoint: synthesis
invoked_by: cluster_claims.py (Phase 4)
model: claude-opus-4-7 (or current Opus)
status: draft
---

# Synthesis claim-clustering prompt — **substantive**

**Classification rationale**: this prompt asks Claude to *identify thematic relationships* across `#inferred` / `#speculative` claims from different chapters. Per Crocco: "Generating thematic categories" and "thematic identification across literature" are substantive uses (§"Pragmatic vs Substantive"). Output requires human verification before any committed change.

Compliance requirements when this prompt is used:
1. The synth report **must** quote the exact prompt + the model + the commit SHA of this file.
2. The human accepts/rejects each proposed cluster explicitly in the report.
3. Any cluster that drives a chapter edit must be re-justified by the human in their own words in the chapter PR description.

## Prompt

```
You are analysing claim-tags from a physics-history campaign organised around Tepper
Gill's dual-theory of relativity and quantum mechanics.

I have collected every `#inferred` and `#speculative` claim from the campaign's 9
chapters. Each claim is tagged with its location (chapter + section heading + line)
and a snippet of surrounding text.

YOUR TASK: identify groups of claims that share a common conceptual thread — claims
that, if examined together, might reveal either (a) a consistent extrapolation
pattern that strengthens the framework's credibility, or (b) a *tension* where two
chapters extrapolate from the same source paper in incompatible ways.

For each cluster:
1. Name the cluster (≤ 6 words).
2. List the constituent claims by chapter + line.
3. State whether the cluster is *consistent* or *tense*.
4. Quote the shortest piece of text from each claim that demonstrates the shared
   conceptual thread.
5. If *tense*, propose the specific question the human should adjudicate — do NOT
   propose how to resolve it.

Hard rules:
- A cluster needs ≥ 2 claims to be reported.
- Do not invent inter-chapter connections that the text does not support — quote, do
  not paraphrase, when proving the shared thread.
- Single-source clusters (multiple claims all citing one paper) are interesting and
  should be flagged separately from cross-source clusters.
- If you cannot find any clusters above your minimum-confidence threshold, say so
  explicitly — empty output is acceptable.

CLAIMS:
{claims_block}
```

## Inputs

`{claims_block}` is the output of the existing `_tools/build_dataview_indexes.py` claim-indexer, formatted as YAML or JSON per claim: `{chapter, line, section, snippet, tag}`.

## Outputs

Markdown synth report at `Roadmapping/Tooling/synth_reports/cluster_claims_<date>.md` with the structure:
- Frontmatter: tool, run_date, model, ai_use_class=substantive, prompt_ref=<this file>@<commit-SHA>.
- Inputs section: how many claims were analysed; which chapters covered.
- Output section: the model's clusters verbatim.
- **Human acceptance section** (left blank for human fill-in): which clusters accepted? Which rejected? Why?
- **Downstream actions** section: chapter PRs that should follow.

The synth report is committed but never auto-applies anything to chapter content. Chapter edits go through normal PR review.
