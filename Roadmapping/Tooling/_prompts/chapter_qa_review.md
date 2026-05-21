---
ai_use_class: pragmatic
touchpoint: writing
invoked_by: ad-hoc Claude Code chats during QA pass
model: any (Sonnet/Opus/Haiku all fine)
status: draft
---

# Chapter QA-review prompt

**Classification rationale**: this prompt asks Claude to *check* a chapter for surface-level quality issues (grammar, internal consistency of claims, broken wikilinks, missing references). Output is editorial feedback the human accepts or rejects — it does not generate manuscript content. **Pragmatic** per Crocco's distinction (analogous to grammar/copy-editing).

## Prompt

```
You are doing a copy-edit + internal-consistency review of one chapter of a
physics-history campaign. The chapter is already drafted; your job is to catch
errors and inconsistencies, NOT to rewrite or restructure.

Check, in order:

1. **Grammar, typos, sentence-level prose.** List concrete fixes; no rewriting beyond
   the level of swapping a single word or punctuation mark.

2. **Internal-claim consistency.** Are there places where the chapter contradicts
   itself? (E.g., paragraph 3 says X, paragraph 8 says ~X.) Quote both and ask
   the human to reconcile.

3. **Wikilink targets.** Are any `[[citekey]]` references that look like bib stubs
   missing from `Roadmapping/History/Bibliography/{Primary,Retrospective}/`?
   (You may not have access to that directory — if so, say so and skip this step.)

4. **Framing-principle compliance.** Per the campaign's load-bearing rule:
   `#inferred` and `#speculative` claims must be qualified with a
   precision-regime statement. List any tagged claims that lack this
   qualification.

5. **Surface-level claims you cannot verify.** Flag specific assertions of fact
   (dates, numerical values, paper attributions) where you couldn't independently
   confirm them from the chapter text alone. Quote them — don't invent corrections.

Hard rules:
- Do NOT rewrite paragraphs. Suggest targeted edits only.
- Do NOT invent citations to fill perceived gaps.
- Do NOT claim a fact is wrong without quoting both the chapter's version and
  the basis for your concern.
- If a section is good, say so and move on — empty findings are valid output.

CHAPTER:
{chapter_body}
```

## Inputs

`{chapter_body}` is the full markdown of one chapter.

## Outputs

Markdown report (typically pasted into a Claude Code conversation rather than committed). The human selects fixes worth applying and commits them as a normal chapter edit.

No synth-report file is created for this prompt — it's an interactive QA aid, not a synthesis step.
