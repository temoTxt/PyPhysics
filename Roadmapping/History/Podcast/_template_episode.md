---
episode: NN
title: "<Episode title>"
era: "<YYYY-YYYY or 'forward'>"
chapter: <NN_chapter_file>
speakers: [Historian, Physicist, Experimentalist]    # drop Experimentalist for 2-voice forward episodes
target_runtime_min: 35                                # 30–45 min typical (≈5,000–7,000 words)
animations_cued: []                                   # scene-file basenames; auto-resolved by lint_episode.py
status: draft                                         # draft | reviewed | tts-rendered
---

# Episode NN — <Episode title>

> Companion dialogue script for [[<NN_chapter_file>]]. Same primary research, conversational form.

## Cold open

**Historian:** {scene-setting: a famous figure, an experiment in progress, a journal entry — sets the era}.

**Experimentalist:** {hook: "and here's what was puzzling about that result" / "the apparatus they were using was…"}.

**Physicist:** {within the first ~60 seconds, deliver the framing-principle reminder in their own voice — something like:} "Before we go further: what we're doing in this series is exploring how the same experimental record can be re-expressed in a different mathematical convention — Gill's proper-time framework. We're not proposing new physics. Where the dual-theory framing makes a different numerical prediction from standard SR or QM, we'll flag the precision regime that could tell them apart. Most of the time, that regime is beyond current measurement."

`[cue: title card]`

## Historical sweep

Historian leads; others interject with questions and elaborations. Wikilinked primary sources can be read aloud as "as Maxwell put it in his 1865 paper…":

**Historian:** …

**Experimentalist:** {questions about what was actually measured}.

**Physicist:** {questions about the equations being introduced}.

`[cue: animation <scene_name> ]`

## Physics deep dive

Physicist takes the floor; works through one or two key derivations at chalk-talk pace. Stage directions cue Manim scenes for the produced audio version.

**Physicist:** …

`[cue: animation <scene_name> ]`

**Experimentalist:** {challenges with "and what does that predict for the experiment?"}

## Proper-time interlude

The "now let's contrast that with how Gill's framework would frame this" beat. Physicist leads.

**Physicist:** {introduces the dual-theory framing, cross-referencing [[<verification_doc>#<anchor>]]}.

**Experimentalist:** {"but how would we tell the difference experimentally?"}

**Physicist:** {answers honestly — usually "at current precision, you can't; here's the regime where you could"}.

**Historian:** {connects this back to a historical figure who would have appreciated the move, or to a present-day stake}.

## Closing

**Historian:** {forward pointer to the next era}.

**Experimentalist:** {drops a hook for unresolved questions to be picked up next episode}.

**Physicist:** {one-sentence reminder of the framing principle as the sign-off}.

`[cue: end card with bibliography wikilinks for show notes]`

## Show notes

Auto-generated from the chapter's bibliography by `lint_episode.py`; primary + retrospective sources cited above.
