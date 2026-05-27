# Podcast — 3-voice dialogue scripts

Companion dialogue versions of every chapter. Same primary research, same bibliography wikilinks, same animation cross-references — conversational form. Each script is publishable as written dialogue (a radio play) without ever running TTS; an optional audio pipeline lands in PR L (Phase 10) if/when wanted.

## Cast (fixed across all episodes)

A 3-voice cast keeps the narrative coherent. For the speculative forward chapters (Ch 8 QC, Ch 9 fusion) the Experimentalist may drop out — there's less measured data to interrogate — and the format becomes 2-voice.

- **The Historian** — pushes the chronological narrative forward, cites primary sources, supplies biographical and institutional context. Scholarly but accessible; doesn't lean on jargon.
- **The Physicist** — explains the equations and the proper-time framework; takes the listener through algebra at chalk-talk pace. Slight enthusiasm for the dual-theory contrasts; always returns to the framing principle when a "Gill would predict" sentence comes up.
- **The Experimentalist** — focuses on what was actually measured, what wasn't, what the apparatus looked like, what the data ruled out. Skeptical voice; asks the others to defend their claims with evidence. Particularly important in the Wolfram-verified chapters (Ch 2, 4, 5) where the dual-theory predictions are concrete.

The same three personas carry across all nine episodes. They have continuity of voice — a callback in episode 5 to something the Historian said in episode 2 is legitimate and worth using.

## Persona → voice mapping (placeholder; finalized in PR L)

| Persona | Voice profile | TTS voice ID |
|---|---|---|
| Historian | warm, mid-range, measured pacing | (pending PR L) |
| Physicist | brighter, faster, slight rising inflection on equations | (pending PR L) |
| Experimentalist | lower, drier, occasional pause before challenging a claim | (pending PR L) |

Final voice IDs are written into this table by PR L when the TTS backend is chosen (`pyttsx3` offline, `coqui-tts` neural, ElevenLabs / OpenAI / Anthropic API). Until then, the personas exist only as written style, and the scripts ship as markdown.

## Episode structure

Mirrors the chapter structure but in dialogue form. Target ~30–45 min spoken (≈5,000–7,000 words script):

1. **Cold open + framing-principle disclaimer.** Historian sets the scene with a famous figure or experiment of the era; Experimentalist hooks with "and here's what was puzzling about that result". **Within the first ~60 seconds**, whoever has the floor (usually the Physicist) delivers the framing-principle reminder in their own voice — "we're exploring a different mathematical convention for known physics, not proposing new physics".
2. **Historical sweep.** Historian leads; others interject with questions. Wikilinked primary sources read aloud as "as Maxwell put it in his 1865 paper…".
3. **Physics deep dive.** Physicist takes the floor; works through one or two key derivations. Stage directions `[cue: animation <scene>]` mark where the Manim scenes drop in for the produced audio/video version.
4. **Proper-time interlude.** Explicit "now let's contrast that with how Gill's framework would frame this" moment. Physicist leads; Experimentalist asks "but how would we tell the difference experimentally?". The honest answer is usually "you can't at current precision; here's the regime where you could".
5. **Closing.** Historian closes with a forward pointer to the next era; Experimentalist drops a hook for unresolved questions; Physicist signs off with one-sentence reminder of the framing principle.

## Episode YAML frontmatter

```yaml
---
episode: 02
title: "Classical Synthesis 1860-1900"
era: "1860-1900"
chapter: 02_classical_synthesis_1860_1900
speakers: [Historian, Physicist, Experimentalist]
target_runtime_min: 35
animations_cued:
  - hist_maxwell_synthesis
  - hist_michelson_morley
status: draft   # draft | reviewed | tts-rendered
---
```

The `animations_cued` field lets the optional Phase 10 build pipeline auto-build a video version of the episode with the animations dropped in at the right cue points.

## Dialogue format

Each speaker's lines begin with a bold persona name followed by a colon. Stage directions are inline in italics or set off with `[brackets]`.

```markdown
**Historian:** In the spring of 1864, Maxwell sat down to write…

**Experimentalist:** And the apparatus he was thinking about — Faraday's induction-ring setup — was already a decade old at that point.

**Physicist:** Right. Let me take us into the equation that ties this all together.

`[cue: animation hist_maxwell_synthesis]`

**Physicist:** Start with Ampère's law…
```

`lint_episode.py` (PR C) validates: YAML schema, speakers ∈ canonical cast, every `animations_cued` resolves to a real Manim scene file under `../../Animations/manim_scenes/`, every wikilinked source resolves, rough word-count → runtime cross-check.

## Preshow prep documents

`preshow_*.md` files are a distinct doc type from the `episode_NN_*.md` scripts: a **prompt bank** that supports a recorded conversation with a guest, rather than a written-dialogue script. The 3-voice cast (Historian / Physicist / Experimentalist) still anchors the format, but each menu item carries a voice tag (`[H]` / `[P]` / `[E]`) indicating the natural opener and a citation back to the source verification doc or primary bibliography note. The guest is a fourth voice; the regulars draw the guest into dialogue rather than scripting them.

Frontmatter variant for preshow prep docs:

```yaml
---
doc_type: preshow_prep
title: "<Pre-show title>"
guest: "<Guest name>"
speakers: [Historian, Physicist, Experimentalist]
target_runtime_min: 120                              # 60–180 min typical for a pre-show
buckets: [history, society_culture, experiments, thought_experiments]
status: draft                                        # draft | reviewed | recorded
issue: <NN>                                          # tracking issue
---
```

Naming convention: `preshow_<short_topic>.md` (e.g., `preshow_tepper_proper_time.md`). Preshow docs are **not** validated by `lint_episode.py` — the linter targets the episode-script schema; the preshow schema is intentionally lighter (no `animations_cued`, no chapter binding, voice tags instead of speaker turns).

## Layout

```
Podcast/
├── README.md                            # this file
├── _template_episode.md
├── episode_01_early_electromagnetism.md
├── episode_02_classical_synthesis.md
├── episode_03_old_quantum_theory.md
├── episode_04_quantum_mechanics.md
├── episode_05_QED_renormalization_solid_state.md
├── episode_06_synthesis_divergence_map.md
├── episode_07_PNT_GPS_SLR_QKD.md
├── episode_08_quantum_computing.md
├── episode_09_fusion.md
├── preshow_tepper_proper_time.md        # preshow prep doc (issue #74)
├── lint_episode.py                      # introduced in PR C
├── build_audio.py                       # introduced in PR L (optional)
├── build_episode_video.py               # introduced in PR L (optional)
└── audio/                               # gitignored; TTS output
```
