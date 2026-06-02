---
ai_use_class: substantive
touchpoint: synthesis
invoked_by: ad-hoc Claude Code chats responding to collaborator-forwarded papers (issues, discussions)
model: claude-opus-4-7 (or current Opus)
status: draft
---

# External-suggestion evaluation prompt — **substantive**

**Classification rationale**: when a collaborator forwards a paper through an issue or discussion and asks for an evaluation, the response shapes whether the paper enters the bibliography corpus, drives a chapter edit, or seeds a follow-up research direction. That is squarely Crocco's "synthesis + thematic identification" touchpoint and is **substantive** (§"Pragmatic vs Substantive"). Verdicts produced by this prompt require human verification before any committed change beyond the synth report itself.

Compliance requirements when this prompt is used:

1. The synth report **must** quote the exact prompt + the model + the commit SHA of this file.
2. The human accepts/rejects each verdict explicitly in the report's human-acceptance section.
3. The bib stub for each evaluated paper must already exist (scaffolded via `scaffold_bib_note.py`) before the evaluation runs — Crocco rule 2 (never invent citations).
4. The synth report does not commit chapter or bibliography body changes; chapter edits and `human_reviewed: true` flips happen in separate PRs after the human reads the primary source.
5. Honest "uncertain" verdicts are required when the paper falls outside the evaluator's competence or when the substantive claim cannot be verified without primary-source effort beyond a single read.

## Prompt

```
You are evaluating a paper that a collaborator surfaced via a GitHub issue or discussion
thread in the PyPhysics repository. The repo is built around Tepper Gill's dual theory
of relativity and quantum mechanics — see CLAUDE.md and the FoundationsII-Classical /
Equation_Verification corpus for the framework's specific claims (b/c modifications to
proper time, dual relativistic dynamics, the r_e/r_0 critical-point thread).

INPUTS (one block per paper):
- cite_key: the bib stub already scaffolded under Bibliography/Retrospective/
- converted_md: the marker-pdf-converted markdown path; read the full file
- issue_url: the GitHub thread that surfaced this paper
- forwarder_handle: the collaborator's GitHub handle, for attribution

YOUR TASK: produce a structured evaluation per paper, plus (if the forwarding comment
raised a *methodological* suggestion) a separate engagement section.

For EACH paper:

1. **Plain-English summary (one paragraph).** What does the paper actually claim? Quote
   the central claim verbatim and translate it into accessible language. Note explicitly
   if the author disclaims novel quantitative predictions.

2. **Technical fidelity check (one paragraph).** Are the cited identities / theorems /
   derivations correct? If the central claim reduces to a textbook identity, say so. If
   the paper extends published mathematics in a new direction, quote the new step.

3. **Verdict on the forwarder's four-way axis** (irrelevant / mundane / brilliant /
   crackpot — Zorns-Lemmon's axis from issue #41; substitute the actual forwarder's
   wording if different). Pick ONE primary verdict + at most one secondary qualifier.
   Acceptable to say "mundane core + speculative extensions" or similar; not acceptable
   to refuse to commit to a verdict.

4. **Relevance to the Gill program (one paragraph).** Does the paper engage with — or
   bear on — any specific Gill claim? Cite anchors in the verification docs or
   FoundationsII-Classical when concrete overlap exists. Say "tangential" or "no
   overlap" honestly when that is the case; do NOT manufacture a Gill connection.

5. **Provenance signals.** Author affiliations, peer-review status (journal vs. preprint
   vs. Zenodo / personal-archive), references to author's own prior work vs. external
   citation, self-described scope. These inform whether the verdict can be drawn from
   the paper alone or needs follow-up consultation.

If the forwarding comment raised a methodological suggestion (e.g., "could method X
apply to the Gill program?"), produce a SEPARATE section:

**Methodological engagement.** State the suggestion as forwarded. Then in one or two
paragraphs: where does the analogy hold? Where does it break? What would a Gill-
specific adaptation actually look like? Be willing to say "the spirit transfers but
the specific machinery does not."

CLOSING SECTIONS for the synth report (the evaluator fills these; the human reviews):

- **Candidate follow-up actions.** List as bullet points. Examples: open issue
  "investigate X", scaffold chapter section for Forward chapter N, no action needed.
- **Draft GitHub reply.** Write the actual text the human can copy-edit-post on the
  forwarding thread. Frame per-paper verdicts in the forwarder's own axis; engage with
  any methodological suggestion. Cap at ~400 words. Sign off as the human (NOT as
  Claude); the human posts.
- **Human acceptance section** (MUST be left blank for the human): per-verdict
  accept/reject, methodological-engagement accept/reject, which candidate follow-up
  actions are approved.

Hard rules:
- Quote when proving a claim — do not paraphrase the paper into stronger or weaker
  language than the original.
- If you cannot read past a section because of equation rendering damage or PDF
  conversion artifacts, say so explicitly and bound the verdict accordingly.
- Do NOT scaffold new bib stubs from inside this prompt — referenced papers already
  have stubs by the time this prompt runs (Crocco rule 2).
- The draft reply must NOT post itself: it lives in the synth report until the human
  edits and posts.

PAPER INPUTS:
{paper_blocks}

FORWARDING COMMENT (if any methodological suggestion):
{forwarding_comment}
```

## Inputs

- `{paper_blocks}`: one block per paper with `cite_key`, `converted_md` path, `issue_url`, and `forwarder_handle`. The evaluator (Claude Code in an ad-hoc chat) reads the full markdown — there is no chunked summarization step.
- `{forwarding_comment}`: verbatim text of the comment that surfaced the methodological angle, if any. For issue #41 this is Zorns-Lemmon's bootstrap suggestion.

## Outputs

Markdown synth report at `Roadmapping/Tooling/synth_reports/<date>_issue<NN>_external_suggestions.md` with the structure:

- **Frontmatter**: tool, run_date, model, ai_use_class=substantive, prompt_ref=`evaluate_external_suggestion.md@<commit-SHA>`, issue_ref, forwarder_handle, evaluated_cite_keys.
- **Inputs section**: what was read; which bib stubs are tied to this report.
- **Per-paper evaluation** (one block per cite-key): the five required subsections (summary / technical fidelity / verdict / Gill relevance / provenance).
- **Methodological engagement section** (if applicable).
- **Candidate follow-up actions**.
- **Draft GitHub reply** (the text the human will post).
- **Human acceptance section** — left blank, fenced as `<!-- TODO: human reviews and fills in -->`.

The synth report commits as-is; the draft reply does not post itself; chapter edits or `human_reviewed: true` flips on the cited stubs happen in subsequent PRs after the human reads the primary source.
