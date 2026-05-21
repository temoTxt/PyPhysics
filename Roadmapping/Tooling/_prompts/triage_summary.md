---
ai_use_class: pragmatic
touchpoint: screening
invoked_by: /triage-papers slash command (Phase 5)
model: claude-opus-4-7 (or current Opus)
status: draft
---

# Triage-summary prompt

**Classification rationale**: this prompt asks Claude to *summarise* an abstract for human triage. The summary is read-only-by-human and does not enter the committed manuscript — therefore **pragmatic** per Crocco's distinction (analogous to grammar checking or formatting). The human still makes the keep/skip decision after reading the abstract directly + the summary.

## Prompt

```
You are helping triage a candidate research paper for inclusion in a physics-history
campaign organised around Tepper Gill's dual-theory of relativity and quantum mechanics.

For each candidate below, produce a SINGLE PARAGRAPH (~3-4 sentences) explaining why
this paper might or might not be worth including, in the specific context of the
chapter the queue.md entry suggests. Cite specific aspects of the abstract — do not
hedge with generic statements like "this could be relevant".

Hard rules:
1. NEVER invent claims that aren't in the abstract.
2. If the abstract is missing or empty, say so explicitly — do not infer content.
3. If the paper looks low-quality (no clear method, no concrete results), say so.
4. End your paragraph with one of:
   - "Recommend KEEP." (only if you're confident it's strong + on-topic)
   - "Recommend LATER." (worth reading but not urgent for this chapter)
   - "Recommend SKIP." (low quality or off-topic)
   - "Insufficient information to recommend." (abstract too thin)

The human reading your output makes the final keep/skip/later decision. Your role is to
surface the signal in the abstract for them, not to decide.

CANDIDATES:
{candidates_block}
```

## Inputs

`{candidates_block}` is the top-N candidates from `queue.md` as a YAML or JSON block, each with:
- `title`, `authors`, `year`, `journal`
- `doi` or `arxiv_id`
- `abstract` (raw text, no AI rewriting)
- `suggested_chapter` (e.g., "08_quantum_computing_open_questions")
- `why_surfaced` (the crawler's reason)
- `score` (the crawler's score)

## Outputs

One paragraph per candidate, ending with one of the four recommendation tokens. The human edits the `decision:` field in `queue.md` and may override the recommendation freely.

## What this prompt does NOT do

- Does not write into `queue.md` directly.
- Does not push decisions to Zotero.
- Does not invent citations or supplementary context not present in the abstract.

All of those are downstream of the human's edits to `queue.md`.
