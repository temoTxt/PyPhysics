---
ai_use_class: substantive
touchpoint: synthesis
invoked_by: suggest_cross_refs.py (Phase 4)
model: claude-opus-4-7 (or current Opus)
status: draft
---

# Cross-reference suggestion prompt — **substantive**

**Classification rationale**: this prompt asks Claude to *propose new wikilinks* between chapters and bibliography entries that the author has not yet linked. The connections proposed shape what the manuscript argues. Per Crocco: this is substantive use, requiring full disclosure + human verification.

## Prompt

```
You are reviewing a physics-history chapter for possible cross-references it should add.

I will give you:
- One chapter's body text (markdown).
- A list of bibliography entries (cite_key, title, abstract) that exist in the repo
  but are NOT yet wikilinked from this chapter.

YOUR TASK: identify the top 5 bibliography entries that this chapter ought to cite,
ranked by relevance to specific passages in the chapter.

For each suggestion:
1. The cite_key of the proposed reference.
2. The exact line / section in the chapter where it should be inserted.
3. A one-sentence rationale grounded in concrete content overlap — quote a fragment
   from the chapter AND a fragment from the bib entry to demonstrate the link.
4. Confidence: high / medium / low.

Hard rules:
- ONLY suggest entries that already exist in the bib I provide. NEVER invent a cite_key.
- The chapter is the AUTHORITY; suggest only where the chapter would be improved, not
  where the bib entry would benefit from being cited.
- Prefer few high-confidence suggestions over many low-confidence ones — empty output
  is fine if nothing crosses your bar.
- A suggestion is invalid if you can't quote concrete content from both sides.

Chapter:
{chapter_body}

Unlinked bibliography (cite_key | title | abstract):
{bib_block}
```

## Inputs

- `{chapter_body}`: the full markdown body of one chapter (provided directly; not paraphrased).
- `{bib_block}`: every entry from `Bibliography/{Primary,Retrospective}/` whose cite-key is NOT already linked from this chapter, formatted as a table.

## Outputs

Markdown synth report at `Roadmapping/Tooling/synth_reports/suggest_cross_refs_<chapter>_<date>.md`:
- Frontmatter (model, commit SHA, ai_use_class=substantive).
- Chapter: the target chapter file path.
- Inputs: how many candidate bib entries were considered.
- Suggestions section: the model's ranked list verbatim.
- **Human acceptance section** (blank): which suggestions accepted? Which rejected? Reasons?
- Action items: text-level edits Trey commits if he accepts a suggestion.

Suggestions never auto-apply. The human writes the wikilink edits themselves and runs `validate_wikilinks.py` before commit.
