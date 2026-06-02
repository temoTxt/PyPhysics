# Per-result template — Bethe–Salpeter precision predictions

Each result in the Bethe–Salpeter campaign uses the following structure. This template is the source-of-record for the format; per-PR chapter documents instantiate it for each individual result.

---

### Bethe–Salpeter §N (or Eq. (N.M)) — short name

**Source:** Bethe & Salpeter §N (Springer 1977 reprint, p. X — local PDF only per [`pdf_status: acquired`](../../History/Bibliography/Retrospective/bethe1977_one_two_electron_atoms.md)). *Pragmatic / Substantive AI.*

**As printed in Bethe–Salpeter:** (LaTeX exactly as printed, fair-use quotation; one or two equations max)

**Modern measurement / CODATA value:** (with year + reference; e.g., CODATA-2018 hydrogen Lamb shift `1057.845(9)` MHz)

**Proper-time / dual-theory derivation:** step-by-step under `K` (Eq. I.6 of DRQM I) and/or the dual Dirac equation. Explicit treatment of the `(u/c)^k` terms that distinguish the proper-time formulation from the textbook. Cite the proper-time `K` cheatsheet [`../_proper_time_K_cheatsheet.md`](../_proper_time_K_cheatsheet.md) rather than restating the canonical relations.

**Wolfram MCP check:** single-line WL verifying the non-rel reduction or the algebraic identity used in the derivation. Record `Result: 0 ✅` inline.

**Numerical comparison:**

| Source | Prediction | Residual vs measurement |
|---|---|---|
| Bethe–Salpeter (textbook) | (value with units) | (Δ in MHz / Hz / ppm) |
| Proper-time / dual-theory | (value with units) | (Δ in MHz / Hz / ppm) |
| Modern measurement | (value ± uncertainty) | — |

**Verdict:**
- ✅ proper-time agrees with measurement to within current precision
- ⚠ predicts a distinguishable shift at the current measurement precision floor (record the shift's size and the precision floor in the comparison row)
- 🔴 disagrees with measurement (cross-post to [`FINDINGS_for_author_review.md`](../../Equation_Verification/FINDINGS_for_author_review.md))

**Companion notebook:** [`BetheSalpeter_§N.wl`](../../Mathematica_Notebooks/Quantum_Mechanics/BetheSalpeter_§N.wl) (if non-trivial algebra; for one-line MCP checks, the inline record suffices).

**Branched treatment (when `r_e` is invoked):** PR C (fine structure) and PR F (hyperfine) carry both branches (a) as-published `r_e ≈ 0.499857`, gives `g = -2.0005714`, and (b) corrected `r_e ≈ 0.499420510`, gives measured `g = -2.00231930…`. Both branches are recorded; the verdict block reports both numerical comparisons.

---

## Voice + Crocco notes

- Substantive content (interpretation of any divergence; physical-significance prose) tagged `*Substantive AI.*` at the source line, with per-paragraph `<!-- TODO: human reviews and fills in -->` blocks following the paragraph.
- Pragmatic content (Mathematica verification, transcription, CODATA lookup) tagged `*Pragmatic AI.*`.
- Voice follows [`VOICE_MATCH_GILL.md`](../../Tooling/VOICE_MATCH_GILL.md). "Mathematically equivalent but not physically equivalent" is the anchor phrase for any verdict that records a measurable difference.
