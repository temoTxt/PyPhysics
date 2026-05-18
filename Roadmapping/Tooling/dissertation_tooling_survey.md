# Dissertation Tooling Survey

Implements [issue #34](https://github.com/temoTxt/PyPhysics/issues/34). Web survey across six categories of tooling that bear on writing a PhD dissertation, with evaluations against the criteria spelled out in the issue (open data format, reference-manager interop, LaTeX support, PDF library + annotation, AI features, pricing, privacy, export pathways, citation completeness, open-access friendliness). Honest about uncertainty; pricing and feature sets shift fast.

**Researcher's note:** these evaluations come from knowledge as of early 2026. Pricing in particular drifts. Always confirm current state on the vendor's own site (URLs given for each tool) before committing to a multi-year workflow. **Tagging @olivercrocco** — you were tagged on the reference-manager comparison in issue #9; this survey is the strict superset.

---

## TL;DR — recommendation

For a dissertation candidate already invested in the Obsidian + YAML + Manim + LaTeX idiom this repo uses, the **lowest-regret stack** is:

1. **[Zotero](https://www.zotero.org/) + Better BibTeX + Zotero Connector (browser)** — free, open-source, fully scriptable, exports clean BibTeX, has the deepest ecosystem of any reference manager. Replaces nothing in the existing YAML pipeline; runs alongside it as the *authoritative store* with auto-sync to BibTeX that the existing tools can re-read.
2. **[Obsidian](https://obsidian.md/) + Citations plugin (or Zotero Integration plugin)** — already in use here. The Citations plugin reads Better-BibTeX exports directly, enabling `[[cite:@maxwell1865_dynamical_theory]]`-style live citations from the note vault.
3. **[Semantic Scholar](https://www.semanticscholar.org/) + [ResearchRabbit](https://www.researchrabbitapp.com/)** for literature discovery — both free, both with grounded citation graphs, no hallucination risk because the citation links *are* the data.
4. **[Quarto](https://quarto.org/)** for manuscript preparation if you want a future-proof Pandoc-based pipeline alongside LaTeX. **Overleaf** if you want pure-LaTeX + co-author collaboration in the browser. For a physics dissertation, both options remain defensible; the right call depends on whether your committee expects classical-LaTeX `.tex` files or accepts a Quarto-rendered PDF.
5. **For AI-assisted research:** [Elicit](https://elicit.com/) for question-driven literature search (cites with quotes, low-hallucination), [Consensus](https://consensus.app/) for "what does the literature say about X" Yes/No-with-evidence framings, and **Claude/GPT/Gemini** for general drafting. **Avoid** any AI tool that emits citations without quoted source passages — modern LLMs hallucinate convincing-looking bibliography entries.
6. **Skip** Mendeley, RefWorks, EndNote, ReadCube unless your institution mandates them — all three have lock-in risks (cloud-only-sync, proprietary file formats, or subscription-walls on previously-paid features).

**For the PyPhysics campaign specifically:** the YAML-frontmatter Bibliography under `Roadmapping/History/Bibliography/` is *the canonical store*; Zotero would *mirror* it (not replace it). The `scaffold_bib_note.py` + `build_bibtex.py` + `update_acquisition_tracker.py` tools stay. Zotero is added as a **second store** for capturing-in-the-wild during PhD work, with periodic export → re-import into the YAML pipeline.

---

## 1. Reference managers and citation databases

The category with the highest lock-in risk per category. Picking wrong here means migrating 2,000+ references mid-dissertation.

| Tool | Cost | Open format | Obsidian interop | PDF annot. | Cloud sync | AI features | Notes |
|---|---|---|---|---|---|---|---|
| **[Zotero](https://www.zotero.org/)** | Free; $20/yr for storage > 300 MB | ✅ SQLite + BibTeX export via Better BibTeX | ✅ via Citations or Zotero Integration plugin | ✅ built-in reader (rewritten 2022) | ✅ self-hosted optional | Limited; some via plugins | **Recommended.** Largest plugin ecosystem; open data; survives the vendor. |
| **[Mendeley](https://www.mendeley.com/)** | Free (with caveats) | ⚠️ Cloud-only sync since 2021; export possible but desktop client deprecated | Limited | Built-in (cloud) | Required (cloud-only) | Some | Owned by Elsevier; killed its open desktop-sync API in 2021. Hard pass for new users in 2026. |
| **[Paperpile](https://paperpile.com/)** | $40/yr (academic); $120 commercial | ✅ BibTeX export; Google-Drive-backed | ✅ via Better BibTeX-style workflow | ✅ web + iOS | ✅ Google Drive | ✅ "Paperpile Beta" with AI summaries (2024–25) | Google-Docs-first; great if you write in Docs, friction if you write in LaTeX. Lighter ecosystem than Zotero. |
| **[EndNote](https://endnote.com/)** | $250 one-time (lab licences cheaper) | ⚠️ Proprietary .enl/.data; export via BibTeX | Limited | ✅ | ✅ | Limited | Clarivate-owned. Strong Word integration. Heavy. Many institutions site-licence it. |
| **[RefWorks](https://refworks.proquest.com/)** | Site licence (institutional) | ⚠️ Proprietary | Limited | Web only | Required | Limited | ProQuest-owned. Web-only. Many universities provide it free; lock-in risk if you leave. |
| **[ReadCube/Papers](https://www.papersapp.com/)** | $5/mo academic | ⚠️ Mostly cloud-only | Limited | ✅ strong PDF annotation | ✅ | "SmartCite" for Word | Polished UI; subscription-walled most of its previously-one-time features. |
| **[Citavi](https://www.citavi.com/)** | $179/yr or institutional | ⚠️ Cloud-first since 2022 acquisition | Limited | ✅ | Required | Some | Strong for German PhDs historically. Acquired by Lumivero (2022); changing rapidly. |
| **[JabRef](https://www.jabref.org/)** | Free, open source | ✅ Pure BibTeX-native | Yes via shared .bib | ✅ | None (local) | None | Java app; BibTeX-first. Great for a LaTeX-only workflow but lacks Zotero's plugin ecosystem and browser-capture polish. |

### Notes on individual tools

**Zotero.** The reference manager that wins on *survivability*. The data lives in a local SQLite file you can copy off your machine; the format is documented; export to BibTeX is one click (via Better BibTeX plugin, which is the *de facto* standard). The Zotero Connector browser extension captures from journal pages, ArXiv, JSTOR, Google Scholar, etc. with one click and is the best in the category. The 2022 rewrite of the PDF reader put it on parity with Mendeley/Papers for annotation. **Obsidian integration:** two plugins — `Citations` (lighter, reads `.bib` directly) and `Zotero Integration` (heavier, live integration via Zotero's debug-bridge). For this repo, the Citations plugin is the simpler match because Better BibTeX can auto-export to a `.bib` that lives in the vault.

**Mendeley.** Was great until 2021, when Elsevier killed the desktop-sync API. The remaining product is web-first with a heavily simplified desktop client. The "Mendeley Reference Manager" (the current product) does not have the depth of the old "Mendeley Desktop" plug-ins. Hard pass.

**Paperpile.** A small, strongly-opinionated, Google-Drive-backed tool that's gained ground in the 2020s. Strong if your writing workflow is Google Docs (Paperpile has the best Docs citation plug-in, period). Less compelling if your writing workflow is LaTeX. Subscription-only; the academic price is $40/yr.

**EndNote.** Word's traditional best friend. If your committee insists on Word + tracked changes + ECMA-376 reference fields, EndNote is the path of least resistance. Otherwise it's overkill: it costs $250+, the data format is proprietary, and the in-app PDF reader is weaker than Zotero's. Many universities provide site licences — check first.

**RefWorks / ProQuest tools.** Web-only. Useful only because many institutions provide them free. If you leave the institution, you lose access (export your `.bib` first). Don't make this the centerpiece of a multi-year workflow.

**ReadCube/Papers.** Polished UI; subscription-walled most of its formerly-one-time features in the 2020s. SmartCite (its Word plugin) is a real advantage if you live in Word. Otherwise, no compelling reason over Zotero.

**Citavi.** Was distinctive for its tight integration between *literature review*, *quote/note management*, and *outline → manuscript*. Acquired by Lumivero in 2022; product direction has shifted to a cloud-first subscription model with monthly pricing. Long-standing Citavi users on the desktop perpetual licence are reportedly grandfathered. New users in 2026 should be cautious.

**JabRef.** Pure BibTeX-native. Java app. Great if you live entirely in LaTeX + a `.bib` file checked into git. Lacks Zotero's browser-extension capture polish. **For this repo's existing workflow specifically, JabRef is the closest analog to what we already do** — and a reasonable substitute *or supplement* for the existing YAML-pipeline tools.

### Recommendation for category 1

**Use Zotero as the primary store; mirror to BibTeX via Better BibTeX; let the existing YAML pipeline continue to be canonical for the campaign-specific bibliography.** Skip everything else unless your institution provides EndNote or ReadCube free and you have a specific reason to use them.

---

## 2. PDF discovery and literature mapping

The category that's exploded since ~2020 with AI-assisted citation-graph tools. Critical to do well; these are the tools that find papers you didn't know existed.

| Tool | Cost | Method | Coverage | Citation grounded? | Notes |
|---|---|---|---|---|---|
| **[Semantic Scholar](https://www.semanticscholar.org/)** | Free | AI-augmented citation graph | ~200M papers; broad | ✅ links are data | Allen Institute. Open API (S2 API + Datasets). Foundation for many other tools. |
| **[ConnectedPapers](https://www.connectedpapers.com/)** | Free (5/mo) → $5/mo | Citation-graph visualization | Uses S2 data | ✅ | Beautiful UI; the "spread of papers around this seed" view. |
| **[ResearchRabbit](https://www.researchrabbitapp.com/)** | Free | Citation-graph + curated collections | Uses S2 + others | ✅ | Free tier is generous. Zotero integration is fully supported. |
| **[Litmaps](https://www.litmaps.com/)** | Free (3 maps) → $10/mo | Citation-graph with topical-evolution view | Uses S2 + Crossref | ✅ | Strong for "where is this field going" framing. |
| **[Inciteful](https://inciteful.xyz/)** | Free | Citation-graph + paper-discovery queries | Uses S2 + OpenAlex | ✅ | Lighter UI; query-driven rather than collection-driven. |
| **[OpenAlex](https://openalex.org/)** | Free (open data) | Open citation graph + DOI metadata | ~250M works | ✅ | Replacement for Microsoft Academic (RIP 2021); used as a back-end by many of the others. |
| **[Crossref](https://www.crossref.org/)** | Free (open data) | DOI metadata + citation links | ~150M works | ✅ | Industry-standard DOI registrar; we already use it in `scaffold_bib_note.py --from-doi`. |
| **[Elicit](https://elicit.com/)** | Free (5k credits/mo) → $10/mo | AI Q&A over papers | ~125M papers | ✅ with quoted passages | Question-driven literature review with citation-grounded answers. |
| **[Scite](https://scite.ai/)** | $20/mo academic | Citation-context analysis (supportive / contradictory / mentioning) | Indexed across publishers | ✅ | Distinctive: tells you *whether* a paper was cited approvingly or skeptically. Useful for finding controversies. |
| **[Undermind](https://www.undermind.ai/)** | Free tier → $20/mo | "Deep research" — AI agent runs 30+ min on a question | Uses S2 + others | ✅ | Newest; very long-running queries; quality varies. |
| **[SciSpace (formerly Typeset.io)](https://typeset.io/)** | Free → $20/mo | Q&A "Copilot" + paper formatter | Indexed broadly | ⚠️ mixed; check sources | Multi-tool; the formatter component is interesting. AI features sometimes hallucinate. |

### Notes

**Semantic Scholar** is the foundation. Their dataset (S2ORC, 200M+ papers with abstracts, parsed citations, parsed full-text where openly available) is the back-end most of the visualisation tools above use. The S2 API is open and has a free tier. **If you use no other tool from this category, learn the S2 web UI.**

**ConnectedPapers** vs **ResearchRabbit** vs **Litmaps** vs **Inciteful** — all four show citation graphs, all four use S2 as their data source, and the choice mostly comes down to UI taste. *ResearchRabbit* has the most polished collection-management; *ConnectedPapers* has the most beautiful single-paper view; *Litmaps* has the strongest temporal/topical-evolution view; *Inciteful* is the lightest and most query-driven.

**Elicit** is the strongest *AI-with-citation-grounding* tool in the category. Ask it a question; it returns a table of papers with relevant-passage extracts. The grounding (the quoted passage in the answer) is what protects against hallucination. **This is the right shape for AI literature review.** Beware: Elicit's "summary across papers" feature is more error-prone than its "evidence from each paper" feature; trust the latter more.

**Scite** is the one tool in this category that *isn't* a citation-graph visualizer — it analyses the *sentiment* of citations (supportive, contradictory, neutral). For finding controversies and reproducibility issues, this is uniquely valuable. Worth a subscription if your dissertation is in a field with disputed results.

**Undermind** is the newest and the riskiest bet — it runs an AI agent for 30+ minutes on a research question and returns a curated paper list. Quality is highly variable; results are sometimes excellent and sometimes confidently wrong. Useful as a complement, not as a primary tool.

**SciSpace** is a Swiss-Army-knife (search + Q&A + journal-style PDF formatter). The journal-formatter component is its strongest feature for actual *paper writing* (auto-reformats a manuscript to ~50,000 journal templates). The AI-Q&A component overlaps with Elicit and is weaker.

### Recommendation for category 2

**Semantic Scholar (free) for the citation-graph data; ResearchRabbit (free) for collection-management; Elicit (free → $10/mo) for AI-assisted literature search with quoted-source grounding.** Add Scite ($20/mo) if your field has frequent controversies. Skip Undermind for now (too new); skip SciSpace's Q&A but consider its formatter for end-stage manuscript prep.

---

## 3. Note-taking / second-brain

The category that determines how the dissertation is *organised* during the 3–7-year writing process. Lock-in here is high (years of notes lose value if you migrate badly). This repo is already invested in Obsidian; the question is whether to stay.

| Tool | Cost | Open format | LaTeX | Reference-manager interop | Graph view | Notes |
|---|---|---|---|---|---|---|
| **[Obsidian](https://obsidian.md/)** | Free; $50/yr sync | ✅ plain Markdown + frontmatter | ✅ MathJax | ✅ Citations / Zotero Integration plugins | ✅ | **Already used here.** Local-first; vault is just a folder. Plugin ecosystem huge. |
| **[Roam Research](https://roamresearch.com/)** | $15/mo | ⚠️ JSON export | ⚠️ LaTeX inline, limited | ⚠️ via Zotero Roam | ✅ | Pioneered bidirectional links. Subscription. Vendor health uncertain (2024–25). |
| **[LogSeq](https://logseq.com/)** | Free, open source | ✅ Markdown + Org-mode | ✅ KaTeX | ⚠️ via Zotero LogSeq plugin | ✅ | Open-source Roam-alike; local-first; less mature than Obsidian. |
| **[Notion](https://www.notion.so/)** | Free (limits) → $10/mo | ⚠️ proprietary block format; export available but lossy | ⚠️ weak | ⚠️ limited | ⚠️ limited | Popular and polished, but Math support is weak and the block format is *not* portable. **Skip for a physics dissertation.** |
| **[Tana](https://tana.inc/)** | $14/mo paid; waitlist for free tier | ⚠️ proprietary | ⚠️ limited | ⚠️ limited | ✅ | Opinionated; growing fast 2024–25. Cloud-only. |
| **[Scrintal](https://scrintal.com/)** | $9.99/mo | ⚠️ proprietary canvas | ⚠️ | Limited | ✅ visual | Whiteboard-style; visual canvas. Niche. |
| **[Capacities](https://capacities.io/)** | $10/mo | ⚠️ proprietary | ⚠️ limited | ⚠️ early | ✅ | "Object-oriented" PKM. Interesting model; lock-in risk. |
| **[Heptabase](https://heptabase.com/)** | $9/mo | ⚠️ proprietary; export available | ⚠️ partial | ⚠️ early | ✅ visual whiteboard | Strong for visual thinkers; some PhD candidates love it. Cloud-first. |
| **[Reflect](https://reflect.app/)** | $10/mo | ⚠️ proprietary | ⚠️ limited | ⚠️ limited | ✅ | AI-native; polished UI. Cloud-only. |
| **[Org-mode](https://orgmode.org/) (Emacs)** | Free | ✅ plain text | ✅ via Org's LaTeX export | ✅ via org-ref | ⚠️ via roam packages | The traditional academic note-taking system. Steep learning curve; unmatched flexibility for those who pay it. |

### Notes

**Obsidian** is the default recommendation for this repo's idiom and for academic note-taking generally in 2026. The vault is a folder of plain `.md` files with YAML frontmatter. The graph view, backlinks, and tag pane are native. The plugin ecosystem is unmatched. Local-first means your notes survive the vendor.

**Roam Research** pioneered the bidirectional-link / block-reference model that Obsidian, LogSeq, Tana, etc. all copied. Roam's product trajectory in 2024–25 has been concerning (slower releases, public arguments between founders and users). **Skip for new investment.** Existing users should export to Markdown and consider migrating to Obsidian or LogSeq.

**LogSeq** is the open-source Roam-alike. Local-first, Markdown-or-Org-mode. Younger than Obsidian, smaller plugin ecosystem. Strong if you want Roam's outliner-style block-references in an open format.

**Notion** is the most popular tool in the category overall but is *the wrong choice for a physics dissertation*. The math rendering is weak (basic KaTeX, no per-line block-math support, awkward), the block format is proprietary, and exports are lossy (your nested databases don't round-trip cleanly).

**Tana, Scrintal, Capacities, Heptabase, Reflect** — all newer tools (2022–2025 launches) competing for the post-Roam niche. Each has its enthusiasts. All have lock-in risk because the data formats are proprietary. **None worth betting a multi-year dissertation on right now.**

**Org-mode** is the dark horse. The traditional academic note-taking system among Emacs users. Plain-text. Massively powerful. Has org-ref for references, org-roam for Roam-style backlinks, org's built-in LaTeX export for manuscript generation. The learning curve is steep but the system has been stable since 2003 and will outlive every cloud-based competitor in this table. **If you're already an Emacs user, seriously consider it.** If you're not, Obsidian is the lower-friction choice.

### Recommendation for category 3

**Stay with Obsidian.** Add the *Citations* plugin to bridge with Zotero (export-on-save from Better BibTeX → reads from vault `.bib` file). Skip Roam (declining), Notion (math support weak), and the new entrants (lock-in risk). If you're an Emacs user, consider Org-mode as an alternative; if not, the cost-benefit favours Obsidian.

---

## 4. Writing and manuscript preparation

The category that produces the actual dissertation PDF. The committee's expectations dominate this choice.

| Tool | Cost | Workflow | Reference-manager interop | LaTeX support | Notes |
|---|---|---|---|---|---|
| **[Overleaf](https://www.overleaf.com/)** | Free (1 collab) → $10–25/mo | Cloud LaTeX with co-author editing | ✅ via BibTeX upload or [Zotero auto-sync](https://www.overleaf.com/learn/how-to/Including_a_Zotero_bibliography_in_Overleaf) | Native | The standard for collaborative LaTeX. Most physics PhDs end up here. |
| **[Quarto](https://quarto.org/)** | Free, open source | Pandoc-based; compiles `.qmd` → PDF/HTML/Word/EPUB | ✅ via `.bib` and `csl` | ✅ via Pandoc | Posit's successor to RMarkdown. Multi-format output is *the killer feature* — one source produces dissertation PDF + chapter web pages + ePub thesis. |
| **Pandoc-based custom workflow** | Free | Markdown → LaTeX → PDF via Pandoc | ✅ via `.bib` and `csl` | ✅ via Pandoc | Maximum control; some manual scaffolding. The PhD-thesis-from-Markdown blogosphere has plenty of templates. |
| **[Authorea](https://www.authorea.com/)** | Free (with limits) → $15/mo | Web-based collaborative academic writing | ✅ DOI lookup, BibTeX import | ✅ Limited | Multi-format output (LaTeX, Word, HTML). Niche. |
| **[Manuscripts.app](https://www.manuscripts.app/)** | Free (open source) | Desktop manuscript-prep, JATS-XML internal format | ✅ | Limited | Niche; backed by eLife. JATS-XML output is unusual but useful for some publishers. |
| **[Manubot](https://manubot.org/)** | Free, open source | Markdown + Pandoc + GitHub-Actions CI | ✅ via citation-tags `@doi:...` | ✅ via Pandoc | Scientific-paper-from-GitHub-repo. Strong for open-science; overkill for most dissertations. |
| **[Scrivener](https://www.literatureandlatte.com/scrivener/)** | $59 one-time | Long-form writing with corkboard-style outline | Via external Zotero | Limited | Excellent for the *writing* part of long-form prose. Weak for math-heavy text. Common for humanities PhDs. |
| **[LyX](https://www.lyx.org/)** | Free, open source | "What you see is what you mean" LaTeX editor | Via BibTeX | Native | LaTeX with a GUI. Older audience; declining but still works. |
| **MS Word + EndNote/Zotero** | Word ~$70/yr; EndNote/Zotero as above | Standard humanities/social-science track | Both work | Limited (MathType plugin) | Many committees still require this. Plan for it if so. |
| **Google Docs + Paperpile** | Free + $40/yr | Web-collaborative, citations live | Paperpile is best-in-class | Limited (equations via add-on) | Avoid for physics. |

### Notes

**Overleaf** is the default for LaTeX-dissertation workflow. Free tier supports 1 collaborator (you + your advisor). Premium ($10–25/mo) supports more co-authors and adds GitHub sync. The web-based real-time-collaboration is genuinely good. **The standard physics-dissertation choice.**

**Quarto** is the rising challenger. Same Pandoc back-end, more flexible output formats (PDF, HTML, Word, ePub, slides all from the same source). LaTeX support is via Pandoc's LaTeX writer, which is excellent. The killer feature is *multi-format output*: your dissertation can compile to (a) the official institutional PDF, (b) a publishable HTML chapter for your website, (c) slide decks for your defence, (d) an ePub for ebook distribution. All from one source. **If you're starting fresh in 2026, seriously consider Quarto over Overleaf.** The catch: committees may expect a `.tex` deliverable, and Quarto's `.qmd` → `.tex` round-tripping is one-way. If your committee requires raw `.tex`, use Overleaf.

**Plain Pandoc** is the DIY version of Quarto. More control, more manual scaffolding. Plenty of "I wrote my PhD in Markdown" blog posts with templates. Worth knowing about if you want to engineer your own pipeline.

**Authorea, Manuscripts.app, Manubot** — niche tools, each compelling in a specific situation. Manubot is especially compelling for open-science collaborative papers (the manuscript lives in a GitHub repo with CI-driven PDF builds), but it's overkill for a single-author dissertation.

**Scrivener** is the right tool for *prose-heavy* writing — the corkboard / outline / chunks-and-snippets workflow. Common among humanities PhDs. Less compelling for math-heavy text, but worth knowing about if your dissertation has long expository chapters (Chs 1, 6, 8, 9 of the PyPhysics campaign would benefit from Scrivener's outline view in a way that LaTeX doesn't easily provide).

**MS Word.** Many committees in many fields still require Word with tracked changes. If yours does, plan for it: use EndNote or Zotero (both have Word plug-ins), keep equations sane (MathType plugin or LaTeX-to-image inline), and accept the workflow friction. The submission requirement may also be PDF-only, in which case you can write in LaTeX and just hand over the PDF.

**Google Docs.** Avoid for physics. Math support is afterthought.

### Recommendation for category 4

**For physics dissertation in 2026: Overleaf for raw-LaTeX-with-collaborators; Quarto if your committee accepts the rendered PDF and you want multi-format output (worth pursuing).** Have Pandoc available regardless — it's the universal converter. Skip Authorea/Manuscripts.app/Manubot unless you have a specific reason. Add Scrivener for prose-heavy chapters only if you find LaTeX's lack-of-outline-view frustrating.

---

## 5. AI-assisted research and writing

The fastest-changing category. Hallucination risk is the dominant concern; ground-rules matter more than tool choice.

| Tool | Cost | Strength | Hallucination risk | Notes |
|---|---|---|---|---|
| **[Elicit](https://elicit.com/)** | Free (5k credits/mo) → $10/mo | Question-driven literature search | ✅ low (quotes from source) | Best AI literature-review tool in 2026. |
| **[Consensus](https://consensus.app/)** | Free (limited) → $9/mo | "Yes/No with evidence" framing for empirical questions | ✅ low (quotes + papers) | Strong for "what does the literature say about X". |
| **[Scite](https://scite.ai/)** | $20/mo academic | Citation-context (supportive/contradictory) | ✅ very low (citation classifications are deterministic) | See category 2. |
| **[Undermind](https://www.undermind.ai/)** | Free → $20/mo | Long-running deep-research agent | ⚠️ moderate; variable | Newest, riskiest, sometimes brilliant. |
| **[SciSpace Copilot](https://typeset.io/)** | Free → $20/mo | Q&A and paper-formatting | ⚠️ moderate | The formatter is better than the Q&A. |
| **[NotebookLM](https://notebooklm.google.com/)** | Free | Q&A over *your uploaded sources only* | ✅ very low (cites uploaded sources directly) | Google. Excellent for "Q&A over the 20 papers in your dissertation literature review." |
| **[Perplexity](https://www.perplexity.ai/)** | Free → $20/mo | Web-search-augmented Q&A | ⚠️ moderate; depends on source quality | General-purpose research assistant. |
| **Claude/GPT/Gemini direct** | Free–$20/mo | General drafting, math derivations, code | ⚠️ high without grounding | Use with skepticism for citations; with confidence for drafting. |

### Ground rules for AI use in a dissertation

1. **Never accept a bibliography entry from an AI without verifying it exists.** LLMs hallucinate plausible-looking citations confidently. Always check Crossref / Semantic Scholar / Google Scholar for the DOI before adding to your bib.
2. **Prefer tools that return quoted passages over tools that return summaries.** Elicit's "evidence from each paper" view is safer than its "summary across papers" view. NotebookLM's grounded-in-uploaded-sources mode is safer than general-purpose LLMs.
3. **Disclose AI use per journal/institutional policy.** Many journals in 2026 require disclosure; some prohibit AI in specific roles (e.g., authorship, peer review). Check your institution's dissertation guidelines.
4. **Use AI for drafting + analysis, not for fact-finding.** The right uses are: drafting prose from your own outline, deriving math you can verify, summarising your *own* writing for tightening, suggesting alternative framings. The wrong uses are: "what papers should I cite for X" (hallucinations), "what does Smith 2017 say" (without uploading Smith 2017).

### Recommendation for category 5

**Elicit (free → $10/mo) + Consensus (free) + NotebookLM (free)** as the AI-augmented-research stack. **Add Claude/GPT/Gemini** for drafting (you're likely already using one). **Skip Undermind and SciSpace Copilot** until they mature. **Always verify citations against an authoritative source** before adding to your bibliography.

---

## 6. Dissertation-specific tooling

Tools designed specifically for the multi-year long-form-writing problem, rather than general-purpose academic tools.

| Tool | Cost | Purpose | Notes |
|---|---|---|---|
| **Institutional LaTeX classes** (per university) | Free | Conformant dissertation layout | **Mandatory** if your institution provides one. Check your graduate-school website. |
| **[Scrivener](https://www.literatureandlatte.com/scrivener/)** | $59 one-time | Long-form prose with corkboard + outline | Best-in-class for humanities-style dissertation organisation. |
| **[Citavi](https://www.citavi.com/)** | $179/yr | Integrated literature → outline → manuscript | See category 1. Strong but lock-in risk under new ownership. |
| **[Manuscripts.app](https://www.manuscripts.app/)** | Free | JATS-XML manuscripts | Niche; works for some publishers. |
| **[Bookdown](https://bookdown.org/)** | Free | RMarkdown for book-length output | R-centric; mostly subsumed by Quarto in 2026. |
| **Pandoc + custom LaTeX template** | Free | DIY pipeline | Many "I wrote my PhD in Markdown + Pandoc" templates exist. Search GitHub. |
| **Thesis-template GitHub repos** | Free | Pre-built per-institution templates | Search GitHub for `<your university> thesis template`; almost every R1 has one. |

### Recommendation for category 6

**Start with your institutional LaTeX class** if one exists (almost certainly does for any physics PhD program in the US/UK). Use the rest of the stack (Overleaf / Quarto / Pandoc) on top of it. **Add Scrivener** *only if* you find LaTeX's outline-view-deficit actively frustrating for prose-heavy chapters.

---

## Cross-cutting concerns

### Open data + portability

The single highest-stakes decision in the dissertation tooling stack. A multi-year PhD will outlive several tool subscriptions. **Anything that won't export to plain text / BibTeX / Markdown / LaTeX / DOCX / PDF is a risk.** Specifically:

- ✅ **Safe (data lives in open formats):** Zotero, Obsidian, LogSeq, JabRef, Org-mode, Quarto, Pandoc, Overleaf (you own the `.tex`), Scrivener (exports to RTF/Markdown/PDF), Bookdown/Manubot.
- ⚠️ **Risky (lock-in via proprietary format or cloud-only access):** Notion, Roam, Tana, Scrintal, Capacities, Heptabase, Reflect, Mendeley, RefWorks, ReadCube cloud-only features, Citavi cloud, Authorea (manuscript stays on their server).
- ❌ **Avoid for multi-year work:** Anything with no documented export pathway. Specifically: AI-tool note features that don't export, novel/experimental PKM tools.

### AI hallucination risk

The dominant risk for AI tools in this category, as of 2026:
- **Lowest risk:** tools that ground answers in *specific uploaded sources* (NotebookLM, Elicit's evidence view, Consensus, Scite citation classifications).
- **Medium risk:** tools that do web-search-augmented Q&A (Perplexity, Undermind, SciSpace Copilot). Depends on source quality.
- **Highest risk:** raw LLMs (Claude, GPT, Gemini direct) when asked for citations or factual claims without grounding. Excellent for drafting, math derivations, code; dangerous for "find me three papers about X."

**The protective rule:** if an AI gives you a citation, the next thing you do is look up the DOI on Crossref. Always. Forever. No exceptions until the hallucination problem is genuinely solved (we are not there in 2026).

### Privacy and data sovereignty

A dissertation is sensitive intellectual property. Some tools route your draft through their training pipelines unless you opt out (varies by tool, varies by tier). Specifically:
- ✅ **Local-first / self-hostable:** Zotero (local), Obsidian (local), LogSeq (local), JabRef (local), Org-mode (local), Quarto (local).
- ⚠️ **Cloud-first but with privacy commitments:** Overleaf (commits to not training on user content), Paperpile (Google Drive, subject to Google's policies), Elicit (academic license has commitments).
- ⚠️ **Unclear / variable:** Most consumer LLM tiers (Claude/GPT/Gemini free tier) may train on your inputs. The paid / enterprise / API tiers usually have explicit opt-out. Check before pasting unpublished dissertation work.

**The protective rule:** treat unpublished dissertation drafts as sensitive. Use local tools for the canonical store; use cloud tools only for collaboration. If you must paste dissertation text into a cloud AI tool, use the paid tier with explicit "do not train" commitments.

### Cost over a 5-year PhD

| Stack | Year 1 | Years 2–5 | 5-year total |
|---|---|---|---|
| **Minimal (recommended)** | $0 + $40 (Elicit) + $50 (Obsidian Sync optional) | $90/yr | **$0–450** |
| **Standard** | + $20/mo Scite + $25/mo Overleaf Premium = $540/yr | $540/yr | **~$2,700** |
| **Heavy** | + $20/mo Undermind + SciSpace + EndNote one-time = $250 + $480/yr | $480/yr | **~$2,500–3,000** |

For an academic budget, the minimal stack is likely the right starting point. Adding Scite ($240/yr) or Overleaf Premium ($120–300/yr) is incremental once you have a clear need.

---

## Final recommendation

For a physics PhD candidate already invested in the Obsidian + LaTeX + git idiom of this repo:

1. **Zotero + Better BibTeX + Zotero Connector** as the reference-management *capture and store*. Free.
2. **Obsidian + Citations plugin** to read Better-BibTeX exports directly from the vault. Free.
3. **Semantic Scholar + ResearchRabbit + Elicit** for literature discovery + AI-assisted search. Free → $10/mo for Elicit.
4. **Overleaf** for collaborative LaTeX, OR **Quarto** if multi-format output (PDF + HTML + ePub from one source) is appealing and your committee accepts a rendered PDF. Free.
5. **NotebookLM + Claude/Gemini/GPT** for drafting and Q&A-over-your-uploaded-sources. Free (NotebookLM); $20/mo (LLM tier of choice).
6. **Your institution's LaTeX dissertation class** as the foundation. Free.

Total year-1 cost for the minimal stack: **$0–130** (Elicit + LLM subscriptions are optional; everything else is free).

**Skip:** Mendeley, RefWorks, EndNote (unless mandated), ReadCube, Citavi, Notion, Roam, Tana, Scrintal, Capacities, Heptabase, Reflect, SciSpace Copilot, Undermind. Each has a legitimate use case but each has lock-in or hallucination risk that doesn't justify the cost over the minimal stack.

**For PyPhysics campaign work specifically:** the existing YAML-frontmatter bibliography under `Roadmapping/History/Bibliography/` remains canonical. Zotero is added as a *mirror store* for PhD-context capture, with periodic export → re-import. No changes to the existing `scaffold_bib_note.py` / `build_bibtex.py` / `update_acquisition_tracker.py` tools required.

---

## Follow-up issues worth opening

If this survey's recommendations seem worth pursuing:

1. **Hands-on Zotero ↔ Obsidian pilot.** Set up a real Zotero library; install Better BibTeX with auto-export to `Roadmapping/History/Bibliography/zotero.bib`; install Obsidian Citations plugin; verify that `[[cite:@maxwell1865_dynamical_theory]]` resolves in vault. Confirm the existing YAML pipeline still works alongside.
2. **Quarto pilot for one of the campaign chapters.** Take Chapter 7 (PNT, the most LaTeX-heavy) and try rendering it via Quarto instead of plain Markdown. Evaluate multi-format output (PDF + HTML). Decide whether to migrate the campaign or keep raw Markdown.
3. **AI-tool hallucination audit.** Pick 10 citations the AI tools (Elicit, Consensus, Perplexity, raw Claude/GPT) return for a representative query; check each against Crossref; record the hallucination rate. Decide which tools are trustworthy enough to use without per-citation verification (likely none — but the *rate* matters for workflow).
4. **Institutional LaTeX dissertation class** — when Trey is closer to dissertation-writing, identify the institution's required class file and add it to `Roadmapping/Tooling/` as a known dependency.

---

## Method note

This survey relies on my (Claude Code's) knowledge of the tools as of early 2026. URLs are given for each tool so vendor-current pricing/features can be confirmed at any time. Categories with the highest churn (AI tools, new PKM tools) are most likely to be out-of-date soonest; the reference-manager and writing-tool categories are more stable and the recommendations there are likely to hold for several years.

The survey deliberately *does not* attempt to be encyclopedic. The omitted tools (BibDesk for macOS-only users, Bibliopedia, Mendeley alternatives like RefME, etc.) are not bad — they're just not better than the recommended stack for the specific use-case this issue scopes (writing a physics PhD dissertation in 2026 with the idiom this repo uses).

Closes [#34](https://github.com/temoTxt/PyPhysics/issues/34) when merged.
