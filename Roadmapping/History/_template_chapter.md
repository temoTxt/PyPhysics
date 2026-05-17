---
chapter: NN
title: "<Era title>"
era: "<YYYY-YYYY>"
threads: [electromagnetism, quantum, solid-state]   # any subset
animations: []                                       # filenames under ../Animations/manim_scenes/
verification_anchors: []                             # wikilinks into Equation_Verification/
status: draft                                        # draft | reviewed | shipped
---

# Chapter NN — <Era title> (<YYYY>–<YYYY>)

> **Framing principle (load-bearing).** We are exploring differences in mathematical conventions and their consequences for the physical interpretation of established experiments. We are not inventing new physics. Every "Gill would predict X" claim below is qualified by what experimental precision would distinguish it from the standard prediction; the dual-theory framework reproduces all experimentally confirmed predictions of standard SR + QM + QED within current measurement precision.

## 1. Overview

One paragraph framing the era: what's happening in physics, who the central figures are, what experimental capabilities are coming online, and why this era matters for the standard-vs-proper-time contrast. End with one sentence pointing forward to the next chapter.

## 2. Historical narrative

Standard chronological account of the period. Subheadings by thread or by decade as appropriate. Cite primary sources inline as Obsidian wikilinks:

- [[firstauthorYYYY_slug]] — short context on what the source contributes here.

Quote directly from the converted markdown of a primary source when a passage is load-bearing; cite the file + line number from `Historical_Converted_Markdown/`.

### 2.1 Electromagnetism thread

…

### 2.2 Quantum thread

…

### 2.3 Solid-state thread (if applicable)

…

## 3. Proper-time commentary

The dual-theory reframing, claim by claim. **Every claim is tagged** with one of `#verified` / `#inferred` / `#speculative` / `#gill-silent`, and every `#inferred` or `#speculative` claim is qualified with what experimental regime would distinguish it from the standard prediction.

### 3.1 What's directly verified

For each claim drawn from the verification campaign:

- **<claim>.** `#verified` — see [[<verification_doc>#<anchor>]]. {1–3 sentences of context for the historical narrative above.}

### 3.2 What's mechanically inferred

For each claim that's a direct extrapolation of a published Gill equation:

- **<claim>.** `#inferred` from [[<paper>#<eq>]]. {Inference walked through in 3–6 sentences.} **Experimental distinguishability:** {what precision regime would tell the two framings apart, with order-of-magnitude estimate.}

### 3.3 What Gill is silent on

- **<topic>.** `#gill-silent`. {One sentence acknowledging the gap.}

## 4. Key derivations worth animating

| Manim scene | Status | What it shows |
|---|---|---|
| [`<scene_file>.py`](../Animations/manim_scenes/<scene_file>.py) | rendered \| proposed | one-line description |

Scenes proposed but not yet rendered are tracked as TODOs in the chapter's PR description.

## 5. Primary sources cited

- [[firstauthorYYYY_slug]] — <one-line context>
- …

## 6. Retrospective reviews drawn on

- [[firstauthorYYYY_slug]] — <one-line context>
- …

## 7. Cross-references

- Forward: [[<next_chapter>]]
- Verification anchors: [[<verification_doc>]], …
- Findings: [[FINDINGS_for_author_review]] (if any era-relevant)
