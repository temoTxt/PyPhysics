# Quantum Mechanics — Wolfram Notebooks

Companion Wolfram Language (`.wl`) notebooks for the QM-side verification work. Each notebook is a single-file script runnable under the Wolfram MCP server (or `wolframscript`) per the conventions in [`../../../CLAUDE.md`](../../../CLAUDE.md) "How equation verification works." The three MCP gotchas (single-line code; avoid `V` reserved for Vanadium; avoid `e` resolving to Euler's number; non-commutative `Dot`) apply.

Notebooks here are companions to:

- **Griffiths problem set** (`Roadmapping/Quantum_Mechanics/Griffiths/`): per-problem cells of the form `GriffithsCh{N}_P{N}_{M}.wl`.
- **Bethe–Salpeter scaffold result** (`Roadmapping/Quantum_Mechanics/Bethe_Salpeter/`): currently `BetheSalpeter_S3.wl` for the PR A `K = mc² + p²/(2m)` identity.
- **DRQM I §III.D triangulation** (this issue): [`r_e_triangulation.wl`](r_e_triangulation.wl), the joint-fit notebook for $r_e/r_0$ across the six precision atomic-physics observables, following Tepper Gill's 2026-05-25 author guidance that branches (b) and (c) are bracketing guides rather than theoretical predictions. Companion to [`../../Equation_Verification/FINDINGS_for_author_review.md`](../../Equation_Verification/FINDINGS_for_author_review.md) Finding 2 and [`../../Author_Reports/2026-05_followup_for_gill_re_triangulation.md`](../../Author_Reports/2026-05_followup_for_gill_re_triangulation.md).

Each notebook carries a docstring block stating its scope, the companion markdown doc, and the Crocco compliance posture (pragmatic vs substantive AI per [`../../../CLAUDE.md`](../../../CLAUDE.md)).
