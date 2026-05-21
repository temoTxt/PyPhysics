# Crawler biases

Implements Crocco et al. (2026) reporting expectation #5 — **"Address potential biases introduced by AI tools in literature selection and synthesis. Acknowledge limitations and describe steps taken to mitigate them."**

Each crawler in this directory pulls candidates from a specific external source. None of them is unbiased. This document records the systematic skews per source so the triage human can compensate.

## `from_s2.py` — Semantic Scholar Recommendations

**Coverage:** the [Semantic Scholar corpus](https://www.semanticscholar.org/about) (≈ 200 million papers, AI-extracted via Allen Institute pipelines).

**Known biases:**
- **English-language over-representation.** S2's pipeline produces best metadata for English papers; non-English work is indexed but the recommendation engine's relevance signal weakens.
- **CS / biomedical over-representation.** S2's recommendation engine was trained heavily on CS-conference + biomedical-journal citation graphs. Physics is well-covered, but lower-citation-count physics papers may surface less reliably than equivalent-quality biomedical papers.
- **Recency bias.** Recently published (high-citation-velocity) papers are over-surfaced relative to older equally-influential papers, because the recommendation engine partly uses citation velocity.
- **Author-overlap bonus reinforces the bias.** Our `from_s2.py` adds a small score bonus when authors overlap with our existing bib — this amplifies the *existing* composition of our library, which leans Gill-collaboration-heavy. Over time this could narrow the literature we discover via S2.

**Mitigations:**
- Per-source cap (10/day in `build_queue.py`) limits S2's voice in the queue.
- Run alongside ArXiv + Crossref so single-source surfacings get flagged in the merged `sources` field.
- Periodically re-seed by deliberately adding bib entries from outside the current cluster (chapter scaffolding should explicitly cite work the existing bib does not).

## `from_arxiv.py` — ArXiv API daily listings

**Coverage:** the [ArXiv preprint server](https://arxiv.org/), specifically the physics + quantum + condensed-matter categories selected in the crawler.

**Known biases:**
- **Discipline-restricted by construction.** Default categories are `physics.hist-ph`, `quant-ph`, `cond-mat.supr-con`, `physics.gen-ph`. This intentionally excludes most experimental high-energy, astrophysics, biophysics, and entire non-physics disciplines (history of science as practised by humanities scholars, philosophy of physics, etc.).
- **English-language bias** (same as S2; ArXiv is overwhelmingly English).
- **Geographic bias.** ArXiv has high penetration in North America, Europe, India, and East Asia; uneven in South America, Africa, the Middle East.
- **Self-selection bias.** Researchers who don't post preprints (more common in some sub-fields than others) are systematically under-surfaced.
- **Keyword-density scoring amplifies preexisting tag composition.** A paper that uses our chapter-tag vocabulary verbatim scores higher than one that addresses the same physics with different terminology — under-rewarding cross-disciplinary work.

**Mitigations:**
- Per-source cap (10/day).
- `from_crossref.py` partially compensates: it pulls papers that *cite our existing bib* regardless of whether they're on ArXiv. Highly-cited published papers are very likely to be Crossref-indexed even if they skipped ArXiv.
- Trey's manual Zotero captures should deliberately include papers found *outside* ArXiv (e.g., via Google Scholar, journal alerts).

## `from_crossref.py` — Crossref event data

**Coverage:** every DOI registered with Crossref (≈ 150 million works), which is the vast majority of journal-published academic literature.

**Known biases:**
- **Publication-channel bias.** Papers in journals that don't pay for Crossref membership are missing. Most major Western journals are members; many Global South journals are not. (Open-access publishers in OASPA are usually members; that mitigates some of the gap.)
- **Event-data lag.** The event-data API surfaces citation events with a multi-week delay after publication. Very recent papers (< 4 weeks) won't surface here.
- **Self-citation amplification.** A new paper that cites multiple of our bib entries gets a high `from_crossref.py` score. That's intentional — but it does mean papers from the same author network as our existing bib over-surface.
- **No abstract.** Crossref doesn't reliably surface abstracts. Triage decisions on Crossref-only candidates rely on title + journal + author + citation overlap alone, which can mislead.

**Mitigations:**
- Per-source cap (10/day).
- `from_crossref.py` enforces `--min-overlap 1` by default — only surfaces papers that cite at least one of our entries. We can raise this if the queue gets noisy.
- When a Crossref-only candidate looks promising, the human should look up its abstract via Semantic Scholar or the publisher's page before keeping it.

## `from_zotero.py` — Zotero local-API delta

**Coverage:** whatever Trey captures via the Zotero Connector browser extension.

**Known biases:**
- **100% selection bias.** This is the *most biased* source by design — it reflects exactly what Trey was reading on the day(s) he captured. This is fine for an individual's research practice but means the surfacings here are not in any sense a sample of the field.
- **Capture-context bias.** Trey captures more when he's reading (journals, ArXiv listings, Twitter/Mastodon physics-community shares) than when he's not. Field events (conferences, anniversaries of publications) drive bursts that don't reflect long-term importance.
- **No bias-mitigation needed** — this source is meant to encode Trey's curation choices. The other three sources collectively counteract its narrowness by reaching outside his reading list.

**Mitigations (collective):** running all four crawlers in parallel means Zotero's narrowness is balanced by the breadth of the other three. The `sources: [...]` field in deduplicated candidates makes single-source Zotero surfacings legible.

## Overall composition health-check

A queue dominated by Zotero entries indicates Trey is reading heavily — *good*. A queue dominated by S2 entries indicates the library is becoming insular (S2 keeps suggesting near-neighbours of what's already in our bib). When that happens, the corrective is **deliberate counter-bias action**:

- Run `from_arxiv.py` with a *different* set of categories than the default for one week.
- Manually capture (via Zotero Connector) papers from review articles, dissertations, or anthologies that aren't via the usual channels.
- Search for recent work in adjacent disciplines (philosophy of physics, history of science, science studies) that cites our corpus from outside the physics literature.

## Cross-reference

- Crocco's expectation #5 in [`Roadmapping/Tooling/CROCCO_COMPLIANCE.md`](../CROCCO_COMPLIANCE.md).
- Per-crawler implementations: `from_s2.py`, `from_arxiv.py`, `from_crossref.py`, `from_zotero.py` in this directory.
- Cap enforcement: `Roadmapping/Tooling/triage/build_queue.py`.
