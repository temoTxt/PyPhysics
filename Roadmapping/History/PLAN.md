# History of Physics 1800–1965 — Multi-Step Plan

**Tracks:** GitHub issue #7
**Status:** plan-only (this document). No execution PRs opened yet.

A multi-step plan to recreate the history of physics from 1800 to 1965 — through the birth of particle physics — comparing the **standard development** with **Tepper Gill's proper-time / dual-theory framework**. Five chronological chapters covering three parallel threads (electromagnetism, quantum mechanics, solid-state / superconductivity / transistors), each ending in 1965 at the cusp of the standard model and Moore's law.

**Plus three forward-looking chapters** (`Forward/`) treating quantum computing, fusion, and **Position-Navigation-Timing (PNT)** — GPS, Satellite Laser Ranging (SLR), Quantum Key Distribution (QKD) — as *consequences* of the 1800–1965 arc.

The QC and fusion chapters are **open-question roadmaps** (tag-dominant `#speculative`) — what *might* the proper-time framework say about these fields where Gill hasn't published. The PNT chapter is **derivationally heavier** (tag-dominant `#inferred`) because the relevant proper-time content is in already-verified Gill papers (Maxwell Eqs 1–2, 9, 10–11; TCEP Eqs 4.5, 4.16); the chapter derives PNT basics from first principles, then walks GPS / SLR / QKD as applications, framing each through both standard SR+GR and Gill's framework.

**Each chapter ships in three forms**: a reference markdown document (the "chapter"), a Manim animated video for key derivations, and a **podcast-episode script** — a 2–3-voice dialogue walking through the same material conversationally. Same primary research, three audience-facing artifacts: reader, viewer, listener.

## Framing principle (load-bearing — applies to every chapter, scene, and episode)

> **We are exploring differences in mathematical conventions and their consequences for the physical interpretation of established experiments. We are not inventing new physics.**

Every chapter narrative, every Manim scene, every podcast episode must lead with — or in the closing of — an explicit version of this statement. The dual-theory framework reproduces all experimentally confirmed predictions of standard SR + QM + QED within current measurement precision; the contrast with standard physics is at the level of mathematical conventions (proper time as the natural variable, `b = √(c² + u²)` as the collaborative speed, positive-definite Hamiltonians via squaring) and the interpretive consequences that follow. Every concrete "Gill's framework would predict X" claim must be qualified with what experimental regime would (or would not) distinguish that prediction from the standard one. Findings flagged for author review (e.g., `[[FINDINGS_for_author_review]]`) are *internal-consistency* failures of the cited papers — algebra that doesn't reproduce — not "Gill is wrong about physics".

This principle protects the project from being misread as a fringe-physics venture. The whole point is the *opposite*: an honest examination of how the same experimental record can be coherently re-expressed in a different mathematical idiom, with attention to where that idiom produces sharper or more natural derivations (e.g., dissolution of Klein-Gordon's negative-probability problem, radiation reaction without Lorentz-Dirac, GPS clock synchronization via Maxwell Eq. 9).

## Scoping decisions (already settled)

| Question | Decision |
|---|---|
| Primary output format | Markdown chapters + companion Manim scenes (mirrors verification campaign #3 + animations campaign #5). |
| What when Gill has no published paper covering a topic | Extrapolate cautiously; mark every extrapolation `#inferred`. |
| Solid-state / transistor track | Include with speculative proper-time commentary, also flagged `#inferred`. |
| Existing verification + animations docs | Retroactively converted to Obsidian wikilinks in a small preparatory PR (PR A below). |
| BibTeX bridge | Built in Phase 0 (PR B), so chapters can cite via wikilinks from day one. |
| PDF storage | `Historical_Papers/` gitignored by default. Allow-list small public-domain PDFs via `git add -f`. In-copyright PDFs are never committed; only the converted markdown is, under fair-use academic quotation. |

## Output directory layout

```
.obsidian/                                       # vault config (PR A)
Roadmapping/
├── Tepper_Gill_Papers/                          # existing
├── Converted_Markdown/                          # existing
├── Equation_Verification/                       # existing — wikilinked in PR A
├── Animations/                                  # existing — wikilinked in PR A
├── Mathematica_Notebooks/                       # existing
├── Historical_Papers/                           # NEW (PR B) — gitignored allow-list
│   ├── Primary/
│   │   ├── maxwell1865_dynamical_theory.pdf     # public-domain → committed
│   │   ├── michelson1887_relative_motion.pdf    # public-domain → committed
│   │   └── ...
│   ├── Retrospective/
│   │   ├── pais1982_subtle_is_the_lord.pdf      # in-copyright → NOT committed
│   │   └── ...
│   └── Acquisition_Tracker.md                   # NEW — master status table
├── Historical_Converted_Markdown/               # NEW (PR B) — output of parse_papers.py
│   ├── Primary/
│   │   └── maxwell1865_dynamical_theory/
│   │       ├── maxwell1865_dynamical_theory.md
│   │       └── images/
│   └── Retrospective/...
└── History/                                     # NEW
    ├── README.md                                # methodology + Obsidian conventions
    ├── PLAN.md                                  # this file
    ├── _template_chapter.md                     # copyable scaffold
    ├── 01_early_electromagnetism_1800_1860.md
    ├── 02_classical_synthesis_1860_1900.md
    ├── 03_old_quantum_theory_1900_1925.md
    ├── 04_quantum_mechanics_1925_1948.md
    ├── 05_QED_renormalization_solid_state_1948_1965.md
    ├── 06_synthesis_divergence_map.md
    ├── Forward/                                 # NEW — forward-looking applications + speculative roadmaps
    │   ├── README.md                            # framing per chapter (derivational vs speculative)
    │   ├── 07_PNT_GPS_SLR_QKD.md                # derivational; tag-dominant #inferred
    │   ├── 08_quantum_computing_open_questions.md
    │   └── 09_fusion_open_questions.md          # both tag-dominant #speculative
    ├── Podcast/                                 # NEW — 2-3 voice dialogue scripts per chapter
    │   ├── README.md                            # speaker personas, conventions
    │   ├── _template_episode.md
    │   ├── episode_01_early_electromagnetism.md
    │   ├── episode_02_classical_synthesis.md
    │   ├── episode_03_old_quantum_theory.md
    │   ├── episode_04_quantum_mechanics.md
    │   ├── episode_05_QED_renormalization_solid_state.md
    │   ├── episode_06_synthesis_divergence_map.md
    │   ├── episode_07_PNT_GPS_SLR_QKD.md
    │   ├── episode_08_quantum_computing.md
    │   ├── episode_09_fusion.md
    │   └── audio/                               # TTS output; gitignored
    └── Bibliography/
        ├── README.md                            # cite-key conventions, YAML schema, tags
        ├── build_bibtex.py                      # YAML → bibliography.bib generator
        ├── bibliography.bib                     # generated; gitignored
        ├── Primary/
        │   └── *.md                             # citation cards
        └── Retrospective/
            └── *.md
```

## Five eras + divergence axis

| Ch | Era | EM / QM thread | Solid-state thread | Where proper-time diverges |
|---|---|---|---|---|
| 1 | **1800–1860** | Volta → Ørsted → Ampère → Faraday | (none) | `#inferred`: quasi-static sources, `b ≈ c` throughout. **Null comparison.** |
| 2 | **1860–1900** | Maxwell (1865) → Hertz (1888) → Michelson-Morley (1887) → Thomson electron (1897) | Hall effect (1879), photoconductivity (1873), Braun rectification (1874) | `#verified`: Eqs (3)→(3′) verified in PR #4. **Michelson-Morley** is the natural launch point for the dual-theory contrast. |
| 3 | **1900–1925** | Planck → photoelectric → Bohr → Compton → de Broglie → matrix/wave mechanics | Onnes He liquefaction (1908) → superconductivity (1911) | `#inferred`: proper-time Bohr levels (K replaces H); Compton on free electron unchanged (μ_photon=0 for inertial sources). |
| 4 | **1925–1948** | Schrödinger → Klein-Gordon → Dirac → Anderson positron | Bloch (1928), band theory (Wilson 1931), Meissner (1933) | `#verified`: DRQM I verified in PR #4. **Klein-Gordon's negative-probability problem dissolved** in proper-time. Antimatter via Santilli isoduals. |
| 5 | **1948–1965** | Lamb shift → Schwinger (1947) → Feynman/Tomonaga (1948–50) → Dyson → parity violation (Wu 1957) | Transistor (1947), BCS (1957), Josephson (1962) | `#verified`: DRQM I **Finding 2** — the **headline payoff** of the whole project. `#inferred`: Dyson conjecture dissolved by Gill's time-ordered operator calculus. |
| 6 | **Synthesis** | Divergence map across all eras; reuses + extends [`synthesis_tour.py`](../Animations/manim_scenes/synthesis_tour.py). |

## Obsidian + citation conventions

### Cite-key format

`firstauthor + year + slug`, lowercase, snake-case:
- `[[maxwell1865_dynamical_theory]]`
- `[[lamb1947_fine_structure]]`
- `[[bcs1957_superconductivity]]` (BCS as conventional alias)
- `[[schweber1994_qed_and_men]]`

### Bibliography note YAML schema

```yaml
---
title: "A Dynamical Theory of the Electromagnetic Field"
authors: ["James Clerk Maxwell"]
year: 1865
type: primary                    # primary | retrospective
era: "1860-1900"
tags: [electromagnetism, maxwell-equations, foundational]
journal: "Philosophical Transactions of the Royal Society"
volume: 155
pages: "459-512"
doi: "10.1098/rstl.1865.0008"
pdf_status: out_of_copyright_public   # acquired | pending | unavailable | out_of_copyright_public
pdf_path: "../../../Historical_Papers/Primary/maxwell1865_dynamical_theory.pdf"
converted_md: "../../../Historical_Converted_Markdown/Primary/maxwell1865_dynamical_theory/maxwell1865_dynamical_theory.md"
gill_corpus_overlap: ["Two_Mathematically_Equivalent_Versions_of_Maxwells_Equations"]
---
```

### Tag taxonomy

- **By framework claim type** (confidence tiers, strongest → weakest):
  - `#verified` — cross-ref to a Wolfram-confirmed verification doc.
  - `#inferred` — direct/mechanical extrapolation of a published Gill equation to a topic he didn't cover (e.g., Bohr levels with K replacing H).
  - `#speculative` — beyond mechanical extrapolation; proposes where the proper-time framework *might* bear on an open question, without doing the derivation. Used heavily in the `Forward/` chapters.
  - `#gill-silent` — no claim about proper-time relevance at all.
- **By thread**: `#thread/electromagnetism`, `#thread/quantum`, `#thread/solid-state`.
- **By era**: `#era/1800-1860`, `#era/1860-1900`, `#era/1900-1925`, `#era/1925-1948`, `#era/1948-1965`.
- **By topic**: free-form, e.g. `#hyperfine-splitting`, `#renormalization`, `#superconductivity`.

### Historical chapter template (Ch 1–6)

Each historical chapter has these sections (full template in `_template_chapter.md`):
1. **Overview** — 1-paragraph framing of the era. **Closes with the framing-principle disclaimer** (see "Framing principle" section above) — verbatim or in the author's words.
2. **Historical narrative** — standard physics, key experiments highlighted; wikilinks to primary sources inline.
3. **Proper-time commentary** — section per thread (EM / QM / solid-state); each claim tagged `#verified`, `#inferred`, or `#gill-silent`. **Every "Gill would predict X" claim is qualified by what experimental precision would distinguish it from the standard prediction.**
4. **Key derivations worth animating** — list of which animations exist (cross-ref) or are proposed for this chapter.
5. **Primary sources cited** — wikilink list to `Bibliography/Primary/*.md`.
6. **Retrospective reviews drawn on** — wikilink list to `Bibliography/Retrospective/*.md`.

### Forward chapter templates — two variants, *different epistemic status*

**Variant A — open-question roadmap (Ch 8 QC, Ch 9 fusion).** Template (`Forward/_template_forward_open_questions.md`):
1. **Overview** — what the field is and why it's a consequence of the 1800–1965 arc. **Closes with the framing-principle disclaimer.**
2. **Historical roots** — wikilinks back into Ch 1–5 identifying the moments in the historical arc that seeded this field.
3. **Current state (2026 perspective)** — brief landscape paragraph, *not* a literature review.
4. **Major open questions** — numbered list; each with a one-sentence statement + current status.
5. **Speculative proper-time implications** — for each open question, where the dual-theory framework *might* bear on it. **Every claim tagged `#speculative`**; explicit disclaimer "Gill has not published on this — extrapolated from `[[…]]`".
6. **Experimental tests that could distinguish frameworks** — if any exist or are buildable.
7. **Bibliography** — mostly retrospective reviews + modern literature (the field is ongoing); primary sources sparser than historical chapters.

Animations are sparser: 0–1 per chapter. Natural format is a roadmap diagram (conceptual flowchart of open questions and their interrelationships) rather than a derivation walkthrough.

**Variant B — derivational applications (Ch 7 PNT).** Template (`Forward/_template_forward_derivational.md`). Used when the field has *already-verified* Gill content directly bearing on it — derivations are concrete, not speculative. Sections:
1. **Overview** — what the field is and why it's a consequence of the 1800–1965 arc. **Closes with the framing-principle disclaimer.**
2. **Historical roots** — wikilinks back into Ch 1–5 (e.g., for PNT: time-keeping → atomic clocks; Einstein 1905 SR → time dilation; Einstein 1916 GR → gravitational time dilation; Doppler 1842; Sagnac 1913).
3. **First-principles derivation of the field basics** — equations, step by step, both standard and Gill-framed. Closely mirrors historical-chapter style.
4. **Applications walk-through** — each application gets its own subsection deriving the headline equations and showing where standard vs proper-time framings agree/diverge. For Ch 9 PNT: §A GPS, §B SLR, §C QKD.
5. **What proper-time changes (and what it doesn't)** — concrete `#inferred` or `#verified` claims with numerical estimates. Most often the proper-time framework reproduces standard predictions at current precision; flagging *where* it would diverge sets up future experimental tests.
6. **Bibliography** — mix of primary engineering sources (e.g., for PNT: Vessot 1980 Gravity Probe A, original GPS papers) and modern reviews (Ashby 2003 *Living Reviews in Relativity*).

Animations richer than Variant A: 3–4 per chapter, closer to the historical-chapter Manim density.

## Podcast episodes — conventions

Each chapter ships with a companion podcast-episode script in `Roadmapping/History/Podcast/`. The script is a markdown file with dialogue between 2–3 named speakers, plus stage directions ("`[pause]`", "`[cue: animation maxwell_eq04_dissipative.py]`", etc.). Same primary research as the chapter; same wikilinks to bibliography; conversational form.

### Speaker personas (3 voices)

A fixed cast across all episodes to keep the narrative coherent:

- **The Historian** — pushes the chronological narrative forward, cites primary sources, supplies context. Reads as scholarly but accessible.
- **The Physicist** — explains the equations and the proper-time framework; takes the reader through the algebra at chalk-talk pace. Slight enthusiasm for the dual-theory contrasts.
- **The Experimentalist** — focuses on what was actually measured, what wasn't, what the apparatus looked like, and what the data ruled out. Skeptical voice; asks the others to defend their claims with evidence.

The 3-voice format is the default; for episodes where the third voice would be padding (e.g., the speculative forward chapters where there's less experimental data), the Experimentalist may drop out and it becomes a 2-voice dialogue.

### Episode structure

Mirrors the chapter structure but in dialogue form:

1. **Cold open + framing-principle disclaimer** — Historian sets the scene with a famous figure or experiment of the era; Experimentalist hooks the listener with "and here's what was puzzling about that result". **Then the Physicist (or whoever has the floor) delivers the framing-principle reminder** — "Quick reminder: we're exploring a different mathematical convention for known physics, not proposing new physics" — in their own voice, within the first ~60 seconds.
2. **Historical sweep** — Historian leads, others interject with questions and elaborations. Wikilinked primary sources read aloud as "as Maxwell put it in his 1865 paper…".
3. **Physics deep dive** — Physicist takes the floor; works through one or two key derivations. Stage directions cue the corresponding Manim scenes for the produced audio version.
4. **Proper-time interlude** — explicit "now let's contrast that with how Gill's framework would frame this" moment. Physicist leads; Experimentalist asks "but how would we tell the difference experimentally?"
5. **Closing** — Historian closes with a forward pointer to the next era; Experimentalist drops a hook for unresolved questions.

Target length: ~30–45 minutes spoken (≈ 5,000–7,000 words script per episode).

### Episode YAML frontmatter

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

The `animations_cued` field lets a future agent (or the TTS pipeline) auto-build a video version of the podcast with the animations dropped in at the right cue points.

### Production pipeline (deferred — Phase 9, optional)

The committed artifact is **always the script** (markdown). Audio production is optional and lives in a separate Phase 9 PR if/when we want it:

- TTS via a Python tool (candidates: `pyttsx3` offline; ElevenLabs / OpenAI / Anthropic TTS via API; coqui-tts for open-source neural). Decision deferred to Phase 9.
- Per-speaker voice selection: each persona gets a distinct voice, set in `Roadmapping/History/Podcast/README.md`.
- Output: `Podcast/audio/episode_NN_*.mp3`; gitignored, regenerable from the script.
- Optional: a `build_episode_video.py` that composites the audio with the cued animations into a finished video.

If we never produce audio, the scripts are still publishable as written dialogue — like a series of physics-themed radio plays.

## Twelve-PR execution plan (11 required + 1 optional)

### PR A — Obsidian retrofit *(prep, small)*

- Add `.obsidian/` vault config; enable plugins: **Tag Pane**, **Backlinks**, **Outgoing Links**, **Dataview**, **Templates**.
- `.gitignore` for `.obsidian/`: commit `app.json`, `core-plugins.json`, `community-plugins.json`, `hotkeys.json`, `templates.json`; ignore `workspace.json`, `workspace-mobile.json`, `cache/`.
- Convert relative-URL cross-references in:
  - `Roadmapping/Equation_Verification/*.md` (11 files)
  - `Roadmapping/Equation_Verification/FINDINGS_for_author_review.md`
  - `Roadmapping/Animations/README.md`
- Manual spot-check of the graph view (Trey verifies that the existing 12 docs form an interconnected cluster).
- **Out of scope:** rewriting any content; just link syntax.

### PR B — Phase 0 infrastructure

**Documentation + skeleton:**
- `Roadmapping/History/README.md` — methodology + Obsidian conventions + comparison template.
- `Roadmapping/History/_template_chapter.md` — copyable scaffold for historical chapters.
- `Roadmapping/History/Forward/_template_forward_chapter.md` — copyable scaffold for forward chapters.
- `Roadmapping/History/Podcast/_template_episode.md` — copyable scaffold for podcast episodes.
- `Roadmapping/History/Bibliography/README.md` — cite-key conventions, YAML schema, tag taxonomy.
- `Roadmapping/History/Podcast/README.md` — speaker personas + canonical cast + persona → voice mapping (voices empty until PR L).
- Empty chapter files (`01_…md` through `06_…md`, `07_…md`, `08_…md`) with their respective templates + section headings.
- Update root `CLAUDE.md` with paragraphs on Obsidian conventions, PDF-acquisition workflow, and tooling locations.

**Tools (see "Required tools and scripts" section above):**
- **Refactor** `Roadmapping/parse_papers.py` to accept `--input` / `--output` / `--skip-existing` CLI args; backwards-compatible.
- **New**: `Roadmapping/History/Bibliography/build_bibtex.py`.
- **New**: `Roadmapping/History/Bibliography/scaffold_bib_note.py`.
- **New**: `Roadmapping/History/Bibliography/update_acquisition_tracker.py`.

**Directory + tracker setup:**
- Create `Roadmapping/Historical_Papers/Primary/` + `Retrospective/` + `Roadmapping/Historical_Converted_Markdown/Primary/` + `Retrospective/` (with `.gitkeep`).
- Add `Roadmapping/Historical_Papers/Acquisition_Tracker.md` skeleton listing all ~123 cite-keys with `pdf_status: pending` (use `update_acquisition_tracker.py` to bootstrap).
- Add `Historical_Papers/` to `.gitignore` with allow-list pattern for `*.pdf` matching public-domain entries (revisited per-chapter).

**Dependencies:**
- Add `pyyaml` to root `pyproject.toml` (for bibliography frontmatter parsing).
- Add `requests` to root `pyproject.toml` (for `scaffold_bib_note.py --from-doi` Crossref lookup).

### PR C–G — Phases 1–5 (one chapter per PR)

**PR C also introduces** four per-chapter helper tools (used by every subsequent chapter PR; see "Required tools and scripts" above): `_tools/fetch_pdf.py`, `_tools/validate_wikilinks.py`, `_tools/qa_converted_markdown.py`, `_tools/chapter_status.py`, plus `Podcast/lint_episode.py`. They land in PR C so PRs D–J can use them; their development is part of the PR C scope.

Each chapter PR follows the same **six-step** pattern:

1. **Bibliography stubs** — YAML frontmatter only, for ~10–24 era sources.
2. **PDF acquisition pass** — fetch what's freely available (Royal Society archive, Wikisource, arxiv, archive.org); update `Acquisition_Tracker.md`; commit only public-domain PDFs (allow-listed).
3. **Conversion pass** — `uv run python Roadmapping/parse_papers.py --input Roadmapping/Historical_Papers/Primary --output Roadmapping/Historical_Converted_Markdown/Primary`. Commit converted markdown + images.
4. **Narrative pass** — write the chapter. Inline wikilinks resolve to the now-filled bibliography notes; quote directly from the converted markdown where useful.
5. **QA + retrospective summaries + animation** — eyeball converted equations for OCR errors and fix the load-bearing ones; fill in 2–3-paragraph summaries on retrospective notes; render any new Manim scenes.
6. **Podcast-episode script** — write the dialogue version of the chapter for `Roadmapping/History/Podcast/episode_NN_*.md`. Reuses the same bibliography wikilinks and animation cross-references; converts the narrative into 2–3-voice dialogue (target ~30–45 min runtime). Audio production deferred to optional Phase 9.

Per-chapter sizing (each PR includes chapter narrative + bibliography stubs + animation + podcast script):

| PR | Phase | Chapter | New bib notes | New Manim scenes | Podcast script lines | Est. total lines (excl. converted PDFs) |
|---|---|---|---|---|---|---|
| C | 1 | 1: 1800–1860 | ~13 | 1 | ~600 | ~1,200 |
| D | 2 | 2: 1860–1900 (load-bearing) | ~16 | 2 | ~900 | ~2,000 |
| E | 3 | 3: 1900–1925 | ~20 | 2 | ~800 | ~1,800 |
| F | 4 | 4: 1925–1948 | ~17 | 2 | ~900 | ~2,000 |
| G | 5 | 5: 1948–1965 (headline) | ~24 | 2 | ~1,000 | ~2,300 |
| H | 6 | 6: Synthesis + indexes | 0 | 1 (extends `synthesis_tour`) | ~600 | ~1,100 |
| I | 7 | 7: PNT — GPS / SLR / QKD *(forward, derivational)* | ~20 (incl. Le Verrier 1859, Einstein 1915 Mercury, Vessot 1980, Ashby 2003, BB84 1984) | ~5 (Mercury perihelion, GPS clocks, SLR, BB84, basics) | ~1,000 | ~2,400 |
| J | 8 | 8: Quantum computing open questions *(forward, speculative)* | ~15 (mostly retrospective + modern review) | 0–1 (optional roadmap diagram) | ~700 | ~1,500 |
| K | 9 | 9: Fusion open questions *(forward, speculative)* | ~13 | 0–1 | ~700 | ~1,500 |

**Totals (incl. all three forward chapters + podcast scripts):** ~138 bibliography notes, ~15–17 new Manim scenes, 9 podcast episodes (~7,200 lines of dialogue scripts), 9 reference chapters, ~15,800 total lines of narrative + bibliography + dialogue (PDFs and converted markdown not counted).

### PR H — Phase 6 synthesis

- `06_synthesis_divergence_map.md` — closing chapter: where standard and proper-time agree, where they diverge sharply.
- **New tools** (see "Required tools and scripts"): `_tools/build_dataview_indexes.py` (generates `_index_by_year.md`, `_index_by_tag.md`, `_index_inferred_claims.md` covering every `#inferred` + `#speculative` claim across all chapters), and *optional* `_tools/citation_graph.py` (emits a static citation-network PNG for embedding in the chapter).
- Extend `Roadmapping/Animations/manim_scenes/synthesis_tour.py` to span the full 165-year arc, or add `hist_full_synthesis.py` as a complement.

### PR I — Phase 7: Position-Navigation-Timing (PNT) — GPS, SLR, QKD *(forward, derivational)*

`07_PNT_GPS_SLR_QKD.md` — derivational chapter walking PNT basics → three modern applications. Tag-dominant `#inferred` (not `#speculative`); the verified Gill papers directly bear on PNT through Maxwell Eqs (1–2, 9, 10–11) + TCEP Eqs (4.5, 4.16).

**Historical roots** (back-links into Ch 1–5):
- Doppler 1842 (Ch 1–2) → frequency shift from relative motion.
- Le Verrier 1859 (Ch 2) → discovers 43″/century anomaly in Mercury's perihelion advance; Newtonian gravity + planetary perturbations cannot account for it.
- Maxwell 1865 (Ch 2) → EM signal propagation; light as the metric of distance.
- Einstein SR 1905 (Ch 3) → time dilation `v²/2c²`.
- Sagnac 1913 (Ch 3) → rotation-induced fringe shift; baseline for rotating-frame corrections.
- Einstein GR 1915–16 (Ch 4) → derives Mercury's 43″/century *exactly* from Schwarzschild metric — the **first experimental confirmation of GR**. Gravitational time dilation falls out of the same metric.
- Essen-Parry cesium clock 1955 (Ch 5) → atomic time standard makes the GR clock effects observable.

**Chapter sections**:

**§A — PNT basics from first principles.** Derives distance-from-time-of-flight `d = c·τ`; required precision; clock standards; transmitter-position knowledge. Both standard SR + GR and Gill's framework introduced; **Maxwell Eq. (9) `t = (1/c) ∫b(s) ds` is the natural launching point** for the proper-time framing of clock synchronization across an observer–satellite pair.

**§B — Mercury's perihelion precession — the precursor that motivates GPS's GR correction.** Walks through Einstein's 1915 derivation step by step at the grad-student level:
1. Schwarzschild metric `ds² = (1 − 2GM/(c²r))c²dt² − (1 − 2GM/(c²r))⁻¹dr² − r²dΩ²`.
2. Geodesic equation in the equatorial plane → effective radial potential with an extra `−GM L²/(c²r³)` term beyond Newton.
3. Orbit integral → per-revolution perihelion advance `Δφ = 6πGM/(c²a(1−e²))`.
4. Numerical plug-in for Mercury: 43″/century — exactly Le Verrier's residual.

Frames the result as: **GR is not optional even for slow planetary motion** — at-rest clocks differing in gravitational potential accumulate different proper times. The same physics that produces Mercury's 43″/century produces the GPS satellite's +45 μs/day. `#gill-silent` from the proper-time side (Gill hasn't published on extending the dual theory to GR), but it's the *pedagogical bridge* that makes the GPS section make sense.

**§C — GPS as the most precise everyday test of relativity.** Standard derivation: satellite-frame SR correction ≈ −7 μs/day; GR correction ≈ +45 μs/day; net +38 μs/day. *Without* both corrections, GPS position error grows ≈ 10 km/day. Gill-framed derivation: each clock has its own τ_i; the GPS-Time global τ is related via `H_i/mc²`. At satellite orbital velocities `u/c ≈ 1.3×10⁻⁵`, the collaborative speed `b = √(c²+u²) ≈ c(1 + 8.5×10⁻¹¹)` — so the velocity-correction prediction agrees with standard SR to better than current GPS precision. The GR piece inherits the `#gill-silent` framing from §B.

**§D — Satellite Laser Ranging (SLR).** Round-trip time `2d/c` between ground station and a corner-cube retroreflector on a satellite (e.g., LAGEOS); millimeter-precision distance. Corrections: atmospheric refraction, Earth tides, Shapiro delay (a §B-related GR effect). Gill-framed: TCEP Eq. (4.16) group velocity `v_g = v_g' + v` predicts no detectable deviation at SLR precision (`b ≈ c` for sub-orbital optical photons in vacuum).

**§E — Quantum Key Distribution (QKD).** BB84 protocol: Alice sends polarized photons in random basis (rectilinear/diagonal); Bob measures in random basis; classical sifting + privacy amplification. Security from no-cloning. Eavesdropper detection via QBER threshold. Gill-framed: QKD photons live in standard QM Hilbert space, so the protocol is the same. **Speculative `#speculative` hook**: extending to Gill's KS-Hilbert space (`[[FOUNDATIONS_FOR_QED_I_MATHEMATICAL]]` Sec 3) could redefine the no-cloning argument; effective photon mass μ (`[[Two_Mathematically_Equivalent_Versions_of_Maxwells_Equations]]` Eq. 6) is non-zero only during source acceleration so unaffected for fiber-based QKD; satellite-QKD (e.g., Micius 2017) gets brief mention.

**Bibliography (~20 sources)**: Le Verrier 1859 *Comptes Rendus* "Théorie du mouvement de Mercure", Einstein 1915 *Sitz. König. Preuss. Akad. Wiss.* "Erklärung der Perihelbewegung des Merkur aus der allgemeinen Relativitätstheorie", Schwarzschild 1916, Hofmann-Wellenhof *GPS Theory and Practice*, Ashby 2003 *Living Reviews in Relativity* "Relativity in the Global Positioning System", Vessot et al. 1980 Gravity Probe A (`Phys. Rev. Lett.`), Parkinson et al. 1996 *Global Positioning System: Theory and Applications* Vols I–II, original GPS spec papers; Smith & Christodoulidis 1985 SLR; Bennett & Brassard 1984 BB84, Ekert 1991, Pirandola et al. 2020 *Advances in Quantum Cryptography*, Liao et al. 2018 Micius.

**Scenes (~5)**: see per-chapter scene table.

### PR J — Phase 8: Quantum Computing open questions *(forward, speculative)*

`08_quantum_computing_open_questions.md` — open-question roadmap for QC as a consequence of the historical arc.

**Historical roots** (back-links into Ch 1–5):
- Bohr atom (Ch 3) → quantization of states (qubit prototype).
- Heisenberg matrix mechanics (Ch 4) → linear algebra of QM (the lingua franca of QC).
- Dirac equation (Ch 4) → spin-½ systems (natural qubit).
- BCS Cooper pairs (Ch 5) → superconducting qubits.
- Josephson junctions (Ch 5) → transmon qubits.
- Lamb shift + cavity QED (Ch 5) → circuit QED.

**Major open questions, with speculative proper-time bearing**:
- **Decoherence channels**: does Gill's KS-Hilbert space + time-ordered operator calculus give a different account of decoherence than standard QFT? `#speculative` (extrapolated from `[[FOUNDATIONS_FOR_QED_I_MATHEMATICAL]]` Sec 3–4).
- **Cavity QED coherence**: does the dynamical "effective photon mass" μ in dual EM affect coherence times in superconducting circuit QED? `#speculative` (extrapolated from `[[Two_Mathematically_Equivalent_Versions_of_Maxwells_Equations]]` Eq. 6).
- **QED error-budget verification**: if Dyson's conjecture is dissolved in Gill's framework, does the QED error budget for high-precision QC verification experiments change? `#speculative`.
- **g-factor measurements as QC tests**: could a high-precision superconducting qubit array be designed to distinguish standard QED from Gill's dual QED predictions (cf. `[[FINDINGS_for_author_review]]` Finding 2)? `#speculative` — most concrete experimental hook.
- **Topological QC / anyons**: Gill is silent; flag `#gill-silent`.

**Bibliography**: mostly modern QC reviews (Nielsen-Chuang, Preskill notes, surface-code papers) + a few foundational classics (Feynman 1982 "Simulating Physics with Computers", Deutsch 1985, Shor 1994, Kitaev 2003). Most sources via arxiv (open access).

**Scenes**: optionally one — a roadmap diagram showing the open questions and their interrelationships, possibly extending `synthesis_tour.py`'s style.

### PR K — Phase 9: Fusion open questions *(forward, speculative)*

`09_fusion_open_questions.md` — open-question roadmap for fusion as a consequence of the historical arc.

**Historical roots** (back-links into Ch 1–5):
- Atomic theory + Rutherford nucleus (Ch 3) → atoms exist and have a nucleus.
- Einstein E=mc² (Ch 3) → mass-energy equivalence → basis for fusion energy release.
- Aston mass defect (Ch 3) → binding-energy curve.
- Eddington stellar fusion conjecture (Ch 4) → stars powered by H→He fusion.
- Bethe CNO cycle / pp-chain (Ch 4) → quantitative astrophysical fusion.
- Teller-Ulam thermonuclear (Ch 5) → military fusion.
- ZETA / Tokamak / Lawson criterion (Ch 5) → controlled fusion engineering.

**Major open questions, with speculative proper-time bearing**:
- **Relativistic plasma corrections**: in dense plasmas (e.g., inertial confinement at peak compression) where ions move at significant fractions of c, does the "collaborative speed" `b = √(c²+u²)` differ measurably from c in cross-section calculations? `#speculative`.
- **Modified Lorentz force in plasma**: Gill's Maxwell paper Eq. (18) has a `V/(mcb)` correction term. For typical fusion plasma `V/(mc²) ~ 10⁻⁶` (eV vs MeV rest mass), so the correction is below current detection — quantify exactly how much below. `#speculative`.
- **Coulomb barrier tunneling**: standard nuclear physics; firmly `#gill-silent`.
- **Confinement criteria** (Lawson, triple product, ignition): mostly standard plasma physics; `#gill-silent`.

The forward chapter is honest that fusion is largely `#gill-silent` — most of the interesting connection is *historical-arc context*, not proper-time-frames-the-prediction. The chapter's payoff is showing how nuclear physics emerged from the same 1800–1965 development as the dual-theory thread.

**Bibliography**: Aston classics, Bethe 1939, Lawson 1957, ITER design docs, modern fusion reviews. Mix of pre-1929 PD (Aston), mid-century in-copyright (Bethe), and modern arxiv.

**Scenes**: optionally one — binding-energy curve animated with the historical experimental data points overlaid, or a roadmap diagram.

### PR L — Phase 10 *(optional)*: Podcast audio production

Decoupled from chapter PRs. Tackles end-to-end audio production for the 8 episode scripts written in PRs C–J:

- **TTS-engine decision**: `pyttsx3` offline / `coqui-tts` open-source neural / commercial API (ElevenLabs, OpenAI). Per-speaker voice IDs codified in `Roadmapping/History/Podcast/README.md`. Decision deferred to the start of PR L based on quality/cost tradeoff at that time.
- **New tools** (see "Required tools and scripts"): `Podcast/build_audio.py` (parses script → TTS per speaker → concatenated MP3 per episode) and optional companion `Podcast/build_episode_video.py` (composites audio with cued animations into MP4).
- **TTS-backend dependencies** added to root `pyproject.toml` only in this PR (not earlier).
- All audio + video output gitignored (regenerable from the scripts).

This is an "if we want it" PR — the scripts themselves are publishable as written dialogue without ever running TTS.

## Required tools and scripts

The campaign needs a suite of supporting scripts. They land in the PR that first needs them (most in Phase 0); each is small (<300 lines) and follows the project convention `uv run python path/to/script.py …`. All Python paths assume the root project's 3.14 env unless explicitly noted.

### Phase 0 essentials (block PR B)

These must exist before any chapter PR can start.

| Script | Path | Purpose | Lines |
|---|---|---|---|
| `parse_papers.py` *(refactor)* | `Roadmapping/parse_papers.py` | PDF → markdown via marker-pdf. Existing script; gains `--input` / `--output` / `--skip-existing` CLI args so it can serve both the Gill corpus and the new historical-papers tree. Backwards-compatible when called with no args. | ~30 line diff |
| `build_bibtex.py` | `Roadmapping/History/Bibliography/build_bibtex.py` | Walks `Bibliography/**/*.md` YAML frontmatter; emits `Bibliography/bibliography.bib`. Output gitignored — YAML is canonical. YAML → BibTeX field mapping documented in the script's docstring. | ~150 |
| `scaffold_bib_note.py` | `Roadmapping/History/Bibliography/scaffold_bib_note.py` | Given `--cite-key`, optional `--doi`, optional `--year/--author/--title`, emits a skeleton `Bibliography/Primary/<cite-key>.md` with valid YAML frontmatter. Optional `--from-doi` flag auto-fills metadata via the Crossref REST API. Speeds up "bibliography stubs" step in chapter PRs by ~10×. | ~200 |
| `update_acquisition_tracker.py` | `Roadmapping/History/Bibliography/update_acquisition_tracker.py` | Reads all `Bibliography/**/*.md` frontmatter; regenerates the status table in `Historical_Papers/Acquisition_Tracker.md`. Idempotent. Run after any batch of new bib stubs. | ~120 |

### Per-chapter helpers (introduced in chapter PRs)

These land in the earliest chapter PR that needs them (typically PR C) and persist for the rest of the campaign.

| Script | Path | Introduced in | Purpose | Lines |
|---|---|---|---|---|
| `fetch_pdf.py` | `Roadmapping/History/_tools/fetch_pdf.py` | PR C | Given `--cite-key` and `--url` (or `--doi` + Crossref lookup), downloads the PDF into `Historical_Papers/Primary/` or `Retrospective/`. Updates the bibliography note's `pdf_status` and `pdf_path` fields. Logs to `Acquisition_Tracker.md`. | ~250 |
| `validate_wikilinks.py` | `Roadmapping/History/_tools/validate_wikilinks.py` | PR C | Walks every `.md` in `History/`, `Equation_Verification/`, `Animations/`. Checks each `[[wikilink]]` (including section-anchor forms `[[file#heading]]`) resolves to a real file/heading. Outputs a broken-link report; non-zero exit on failure. **Pre-commit hook candidate.** Crucial once the bibliography grows past ~30 notes. | ~180 |
| `qa_converted_markdown.py` | `Roadmapping/History/_tools/qa_converted_markdown.py` | PR C | Scans `Historical_Converted_Markdown/**/*.md` for marker-pdf OCR failure patterns: lone `V` near equations (Vanadium artifact), exponent malformations like `c^{22}` adjacent to `\pi`, missing `\hbar` in known QM contexts, page-break artifacts, running-header bleed-through. **Flags only; doesn't auto-fix.** Output is a per-paper checklist for the chapter's QA step. | ~280 |
| `chapter_status.py` | `Roadmapping/History/_tools/chapter_status.py` | PR C | Compact dashboard: reads all chapter files + bibliography YAML + Acquisition_Tracker; prints a per-chapter table (bib stubs filled / PDFs acquired / scenes rendered / podcast script status / wikilinks valid). Helps a future agent (or me, mid-campaign) know what's outstanding without spelunking. | ~200 |

### Synthesis tools (PR H)

| Script | Path | Purpose | Lines |
|---|---|---|---|
| `build_dataview_indexes.py` | `Roadmapping/History/_tools/build_dataview_indexes.py` | Generates static markdown index pages mimicking Dataview output for users without the plugin: `_index_by_year.md`, `_index_by_tag.md`, `_index_inferred_claims.md` (every `#inferred` + `#speculative` claim across all chapters with anchors). Idempotent; re-runs as new chapters land. | ~250 |
| `citation_graph.py` | `Roadmapping/History/_tools/citation_graph.py` | *Optional.* Builds a directed graph of citations (chapter → bib note → verification doc → animation) and emits a PNG via graphviz for embedding in the synthesis chapter. Obsidian's built-in graph view already does this interactively; this is for the static publication. | ~180 |

### Podcast tooling (PR C introduces the linter; PR L optional for audio)

| Script | Path | Introduced in | Purpose | Lines |
|---|---|---|---|---|
| `lint_episode.py` | `Roadmapping/History/Podcast/lint_episode.py` | PR C *(first episode)* | Validates a podcast episode script: YAML schema, speakers ∈ canonical cast, every `animations_cued` resolves to a real Manim scene file, every wikilinked source exists, rough word-count → runtime cross-check (5,000–7,000 words ≈ 30–45 min). Run per episode in the chapter PR. | ~200 |
| `build_audio.py` | `Roadmapping/History/Podcast/build_audio.py` | PR L *(optional, Phase 10)* | Parses episode markdown (YAML + dialogue lines like `Historian: …`); calls a TTS engine per speaker with persona → voice mapping from `Podcast/README.md`; concatenates with appropriate pauses; emits `Podcast/audio/episode_NN_*.mp3`. **Configurable TTS backend**: `pyttsx3` (offline), `coqui-tts` (open-source neural), or commercial API. | ~400 |
| `build_episode_video.py` | `Roadmapping/History/Podcast/build_episode_video.py` | PR L *(optional)* | Composites the audio (from `build_audio.py`) with the cued animations (from `animations_cued:` frontmatter) into a finished video using ffmpeg. Output: `Podcast/video/episode_NN_*.mp4`. Only relevant if we want YouTube-ready episodes. | ~250 |

### Nice-to-have / deferred

These can wait until a concrete need arises:

- **Manim scene template generator** — emits a new scene file from the canonical structure (docstring + render command + Scene subclass skeleton).
- **Bidirectional cross-reference checker** — extends `validate_wikilinks.py` to flag asymmetric links (A wikilinks B but B doesn't wikilink A).
- **PDF allow-list manager** — given the gitignore-by-default policy for `Historical_Papers/`, reads `Acquisition_Tracker.md` and ensures every `pdf_status: out_of_copyright_public` PDF is force-added while every `pdf_status: acquired` (in-copyright) PDF is *not* in the index.
- **OCR-fixer for pre-1900 scans** — beyond flagging, attempt automated repair for known patterns in Maxwell-era scans.

### Tooling-conventions checklist

Every new script lands with:

1. **Top-of-file docstring** — purpose + usage (`uv run python …`) + I/O contract.
2. **`if __name__ == "__main__": main()`** entrypoint with `argparse` CLI.
3. **No state writes outside the directory it's in** unless explicitly path-arg'd (e.g., bibliography scripts can write into `Bibliography/`, but `fetch_pdf.py` writes into `Historical_Papers/` only via the configured `--output-dir`).
4. **`--dry-run`** flag on any script that writes / mutates state.
5. **No new top-level dependencies** without explicit note in the PR. The root pyproject already has `arxiv`, `marker-pdf`, `torch`; bibliography tooling shouldn't need anything beyond `pyyaml` + `requests` (Crossref). TTS dependencies land only in the Phase 9 PR.

## Per-chapter scene proposals

| Chapter | Proposed new scenes |
|---|---|
| 1 | `hist_faraday_induction.py` — field lines through a moving loop; induced EMF. |
| 2 | `hist_michelson_morley.py` — predicted fringe shift vs measured null. `hist_maxwell_synthesis.py` — Ampère + Faraday → wave equation → `c` emerges. |
| 3 | `hist_bohr_proper_time.py` — Bohr levels rederived with `K = H²/(2mc²) + mc²/2`. `hist_compton_null.py` — Gill's framework recovers standard Compton for inertial electrons. |
| 4 | `hist_klein_gordon_vs_dual.py` — negative-probability problem and how the dual K dissolves it. `hist_positron_isodual.py` — Anderson + Santilli isodual reinterpretation. |
| 5 | `hist_lamb_shift_contrast.py`. Reuses [`drqm_eq22_g_factor_finding.py`](../Animations/manim_scenes/drqm_eq22_g_factor_finding.py) as the chapter payoff. |
| 6 | Extends [`synthesis_tour.py`](../Animations/manim_scenes/synthesis_tour.py). |
| 7 *(forward, derivational)* | `pnt_basics_time_of_flight.py` — distance = c·t; introduces atomic-clock precision. `pnt_mercury_perihelion.py` — animated orbit + Einstein 1915 derivation landing on 43″/century. `pnt_gps_relativity.py` — SR (−7 μs/day) + GR (+45 μs/day) clock corrections, framed via Maxwell Eq (9) `t = (1/c) ∫b ds`. `pnt_slr_geometry.py` — round-trip laser ranging + Sagnac. `pnt_qkd_bb84.py` — BB84 protocol + no-cloning. |
| 8 *(forward)* | Optional `forward_qc_roadmap.py` — open-question roadmap diagram; flowchart-style rather than derivational. |
| 9 *(forward)* | Optional `forward_fusion_binding_curve.py` — binding-energy curve with historical data overlay. |

## OCR-quality + copyright caveats

1. **Modern PDFs (post-1980, arxiv-style)** convert cleanly via `marker-pdf`.
2. **Mid-century journal scans (1947 Lamb-Retherford, 1957 BCS)** convert OK with some equation cleanup needed. Plan a per-chapter QA pass.
3. **Pre-1900 scans (Maxwell 1865, Hertz 1893)** require significant manual cleanup of equations and figures. Budget extra QA time.
4. **Copyright reality** — two-tier policy:
   - `pdf_status: out_of_copyright_public` → commit the PDF.
   - `pdf_status: acquired` (in copyright) → DO NOT commit the PDF. Convert locally, commit *only* the markdown if licensing allows (likely fair-use for academic quotation), and note source in YAML.
   - For books we cite but can't redistribute, convert and commit only the chapters relevant to our citations, with an explicit `pages_quoted` field in YAML.

## Cross-referencing the verification campaign

The 11 papers verified in PR #4 are the "anchor texts" for proper-time commentary. Each chapter cross-references them via wikilinks (post-PR-A):

| Chapter | Verification docs cross-referenced |
|---|---|
| 2 | `[[Two_Mathematically_Equivalent_Versions_of_Maxwells_Equations]]` |
| 3 | `[[FoundationsII-Classical]]`, `[[Analytic_Representation_of_The_Square-Root_Operator]]` |
| 4 | `[[Dual_Relativistic_Quantum_Mechanics_I]]`, `[[Analytic_Representation_of_The_Dirac_Equation]]`, `[[FoundationsII-Classical]]` |
| 5 | `[[Dual_Relativistic_Quantum_Mechanics_I]]`, `[[FOUNDATIONS_FOR_QED_I_MATHEMATICAL]]`, `[[Feynman_Operator_Calculus_Papers]]`, `[[Relativistic_Transformations_of_Thermodynamics]]`, `[[FINDINGS_for_author_review]]` |
| 7 *(forward, PNT)* | `[[Two_Mathematically_Equivalent_Versions_of_Maxwells_Equations]]` (Eqs. 1–2 velocity/time duality, Eq. 9 t-τ integral, Eqs. 10–11 boost), `[[The_Classical_Electron_Problem]]` (Eqs. 4.5 Doppler, 4.16 group velocity), `[[Dual_Relativistic_Quantum_Mechanics_I]]` (Eq. II.3 + KS-Hilbert space — QKD anchor) |
| 8 *(forward, QC)* | `[[FOUNDATIONS_FOR_QED_I_MATHEMATICAL]]` (KS-Hilbert space, time-ordered operator calculus), `[[Two_Mathematically_Equivalent_Versions_of_Maxwells_Equations]]` (μ photon mass), `[[FINDINGS_for_author_review]]` (Finding 2 as QC test hook) |
| 9 *(forward, fusion)* | `[[Two_Mathematically_Equivalent_Versions_of_Maxwells_Equations]]` (Eq. 18 modified Lorentz force; collaborative speed `b` in plasma) |

## Open items left for execution time

- Final per-chapter bibliography lists (specific cite-keys with paper titles) — drafted during the bibliography-stub step of each chapter PR.
- Specific PDF source URLs for each cite-key — populated into `Acquisition_Tracker.md` as PDFs are sourced.
- The exact YAML schema for the `gill_corpus_overlap` field — finalized in PR B once the Obsidian retrofit (PR A) sets the canonical filenames.
- **TTS backend choice** for Phase 10 — deferred until PR L starts; depends on quality/cost tradeoff at that time.
- Whether to LaTeX-typeset a final bound document (e.g. via pandoc) at the end — deferred; not in this plan.

## Tool development summary

For quick orientation, every script the campaign needs, by phase:

| Phase | Tool | Status | New deps |
|---|---|---|---|
| 0 (PR B) | `parse_papers.py` *(CLI refactor)* | exists; refactor | none |
| 0 (PR B) | `Bibliography/build_bibtex.py` | new | pyyaml |
| 0 (PR B) | `Bibliography/scaffold_bib_note.py` | new | pyyaml, requests |
| 0 (PR B) | `Bibliography/update_acquisition_tracker.py` | new | pyyaml |
| 1 (PR C) | `_tools/fetch_pdf.py` | new | requests |
| 1 (PR C) | `_tools/validate_wikilinks.py` | new | none |
| 1 (PR C) | `_tools/qa_converted_markdown.py` | new | none |
| 1 (PR C) | `_tools/chapter_status.py` | new | pyyaml |
| 1 (PR C) | `Podcast/lint_episode.py` | new | pyyaml |
| 6 (PR H) | `_tools/build_dataview_indexes.py` | new | pyyaml |
| 6 (PR H) | `_tools/citation_graph.py` *(optional)* | new | graphviz, pygraphviz |
| 10 (PR L, optional) | `Podcast/build_audio.py` | new | TTS-engine of choice |
| 10 (PR L, optional) | `Podcast/build_episode_video.py` | new | ffmpeg-python |

Total new tools: **12 essential + 2 optional**, ~2,700 lines of Python across the whole campaign. All single-file, single-purpose, `uv run python …` invocable. Net new root-level dependencies through PR K: **2** (`pyyaml`, `requests`); through PR L: **+TTS backend + ffmpeg-python**.
