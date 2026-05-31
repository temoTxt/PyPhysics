# Task 103: [White paper] Research arc — dual/proper-time theory applied to IBM quantum computer coherence (microwave-coupling angle/shape focus)

## Objective

Draft a **3-page white paper** committed at `Roadmapping/Author_Reports/2026-05_quantum_coherence_whitepaper.md` (+ built `.tex` + `.pdf` via the `build_report.sh` pipeline, page count strict at 3 within the `[3, 7]` budget) that proposes a research arc applying the dual / proper-time theory framework to **IBM quantum computer coherence**. The white paper must center on the user's triage-time hook: the proper-time Maxwell equation's **extra dissipative term** ([`Two_Mathematically_Equivalent_Versions_of_Maxwells_Equations.md:31`](../../Roadmapping/Equation_Verification/Two_Mathematically_Equivalent_Versions_of_Maxwells_Equations.md#L31), [line 188 onward](../../Roadmapping/Equation_Verification/Two_Mathematically_Equivalent_Versions_of_Maxwells_Equations.md#L188)) and whether its $\mathbf{u}\!\cdot\!\mathbf{a}/b^4$ structure contributes additional noise via the **microwave coupling** to the qubit — a contribution that depends on the angle between the microwave drive vector and the qubit's effective dipole / charge-displacement axis. The white paper proposes a phased research arc whose Phase B prediction is: **changing the angle or shape of the microwave coupling reduces this framework-predicted noise contribution**, falsifiable on IBM Quantum hardware.

## Background

**The extra term, located precisely.** The verification of Gill's dual Maxwell equations records the dissipative-term equation at [`Two_Mathematically_Equivalent_Versions_of_Maxwells_Equations.md:188`](../../Roadmapping/Equation_Verification/Two_Mathematically_Equivalent_Versions_of_Maxwells_Equations.md#L188) (Eq. (4) — *Dual wave equations for **E**, **B** with dissipative term*, ✅ verified). The discriminating structural claim is at [line 232 of that file](../../Roadmapping/Equation_Verification/Two_Mathematically_Equivalent_Versions_of_Maxwells_Equations.md#L232): *"Gill's Eq. (4) has the same form under the rule $c\to b$, $t\to\tau$, $\mathbf{w}\to\mathbf{u}$, **plus** the extra term"*, with the next line ([234](../../Roadmapping/Equation_Verification/Two_Mathematically_Equivalent_Versions_of_Maxwells_Equations.md#L234)) noting *"which arises only because $b$ depends on $\tau$. In the unaccelerated limit $\mathbf{a}=\mathbf{0}$, $b$ is constant, this term vanishes."* The third-term partner in the proper-time Liénard–Wiechert is recorded at [line 332](../../Roadmapping/Equation_Verification/Two_Mathematically_Equivalent_Versions_of_Maxwells_Equations.md#L332) of the same file: *"the third term is new and tied to the dissipative coefficient $(\mathbf{u}\!\cdot\!\mathbf{a})/b^4$"*. **This is the apparatus-side hook the white paper must lead with**: the term's angle dependence (via $\mathbf{u}\!\cdot\!\mathbf{a}$) and shape dependence (via the source's acceleration profile $\mathbf{a}(\tau)$) is exactly the variable an experimenter controls when shaping the microwave drive of a transmon qubit.

**Where the campaign currently sits.** The dual / proper-time framework has been validated against precision atomic-physics observables (Bethe-Salpeter campaign [#50](https://github.com/temoTxt/PyPhysics/issues/50); Li$^{2+}$ Z-extension [#78](https://github.com/temoTxt/PyPhysics/issues/78); hydrogenic-ion g-factor Z-scan [#82](https://github.com/temoTxt/PyPhysics/issues/82) with Outcome C — the framework's cutoff has no Z-dependence and inherits QED's bound-state structure). At the precision the apparatus currently delivers (leading-$g_s$ / Bethe-estimate), the framework reproduces textbook QED *by construction* — not because it has been independently corroborated by experiment, but because its predictions reduce to the textbook formulas at that precision floor. Applying the framework to superconducting transmon qubits is therefore a substantive conceptual leap: the natural domain is bound electrons in atoms, while transmon qubits are condensed-matter quantum devices. The white paper must front-load this honestly.

**Why the microwave-coupling hook is the right Phase A entry point.** The dissipative term in Eq. (4) is the framework's *only* structural distinction from textbook Maxwell at the classical-field level (the other dual modifications are $c \to b$ substitutions or kinematic re-parametrisations). It vanishes when the source is unaccelerated. Transmon-qubit gate operations drive the Cooper-pair condensate with shaped microwave pulses — i.e., the source IS accelerated, putting the framework's extra term in a regime where it can contribute. The angle/shape dependence ($\mathbf{u}\!\cdot\!\mathbf{a}/b^4$) is the variable the experimenter controls. This makes Phase A concrete: derive the dissipative-term contribution to the effective qubit-drive Hamiltonian and compute its angle/shape dependence under realistic drive parameters.

**Plan-file format precedent.** The shape and depth of this plan follow [`.dev/tasks/61_Triangulate_optimal_r_e_r_0_across_6_precision_observables.md`](../../.dev/tasks/61_Triangulate_optimal_r_e_r_0_across_6_precision_observables.md) (≈128 lines; in-line code refs to specific files+lines; explicit Crocco-compliance posture; per-step substantive-AI flagging). The author-report format precedent is [`Roadmapping/Author_Reports/2026-05_collaborator_intro_spectroscopy.md`](../../Roadmapping/Author_Reports/2026-05_collaborator_intro_spectroscopy.md) (the most recent author report, in [PR #89](https://github.com/temoTxt/PyPhysics/pull/89)): YAML frontmatter, a header block (`**Date:**` / `**From:**` / `**Re:**`), substantive body, optional backing-table at end, terminal `<!-- TODO: human reviews and fills in -->` block.

**Crocco compliance posture.** Choice of the three-phase research-arc framing, the Phase-A null-result-most-likely position, the microwave-coupling angle/shape as the specific predictive hook, and the wording of the honest-scope paragraph are **substantive-AI choices**. Each requires per-section `<!-- TODO: human reviews and fills in -->` markers per [`CLAUDE.md`](../../CLAUDE.md) Crocco §1. The build pipeline (`build_report.sh`) is pragmatic.

**Build-pipeline gotchas already encountered.** The `build_report.sh` Unicode substitution table does **not** cover `⇒` (U+21D2) or `⁺` (U+207A); both have bitten prior reports (PR #68 and PR #89 respectively). The implementer must either avoid those characters or use LaTeX math mode (e.g., `$\Rightarrow$`, `$^{2+}$`) for inline ion-like notation. The page-count check is strict — exactly 3 pages within the `[3, 7]` budget (failure: rebuild with content trimmed or expanded).

## Implementation Plan

1. **Create the author-report file.** Path: `Roadmapping/Author_Reports/2026-05_quantum_coherence_whitepaper.md`. YAML frontmatter (title / author / date / subject) following the [`2026-05_collaborator_intro_spectroscopy.md`](../../Roadmapping/Author_Reports/2026-05_collaborator_intro_spectroscopy.md) precedent. Header block (`**Date:** 2026-05-31. **From:** Trey Morris (with Claude Opus 4.7). **Re:** ...`) names the audience choice at draft-time (issue #103's *Audience options* section lists four; default to *internal-campaign planning* unless Trey overrides at draft start).

2. **Page 1 — motivation + the dissipative-term hook.** Two paragraphs:
   - *What's known apparatus-side:* the dual Maxwell Eq. (4) extra term with $(\mathbf{u}\!\cdot\!\mathbf{a})/b^4$ structure (cite [`Two_Mathematically_Equivalent_Versions_of_Maxwells_Equations.md` line 232](../../Roadmapping/Equation_Verification/Two_Mathematically_Equivalent_Versions_of_Maxwells_Equations.md#L232) for the equation and line 332 for the corresponding third-term in the modified Liénard–Wiechert).
   - *Why it might bear on transmon-qubit coherence:* the term is nonzero only under source acceleration; microwave-driven Cooper-pair motion accelerates the source; therefore the term contributes to the radiation reaction the qubit sees, with magnitude depending on the angle between the drive direction and the Cooper-pair velocity. Conclude with a one-sentence statement of the falsifiable hypothesis: **shaping the microwave coupling to minimise $\mathbf{u}\!\cdot\!\mathbf{a}/b^4$ should reduce framework-predicted noise — and standard QED predicts no such angle/shape dependence at this order**.

3. **Page 2 — three-phase research arc.** Replicate the structure already drafted in the issue body (Phase A / B / C), tightened to fit one page, with each phase having explicit acceptance criteria:
   - **Phase A** — apparatus mapping. Derive the dissipative-term contribution to the effective qubit-drive Hamiltonian (or equivalent Lindblad operator) starting from Eq. (4). Acceptance: a closed-form (or numerical) expression for the noise contribution as a function of the microwave drive angle $\theta$ and pulse-shape envelope $A(\tau)$. **Null-result-most-likely condition**: if the contribution is identically zero in the non-relativistic transmon-qubit limit ($\beta = |\mathbf{u}|/b \to 0$ where $b \to c$), the program closes honestly at Phase A.
   - **Phase B** — concrete predictions. Specific angles and pulse-shapes where the dissipative-term contribution is maximised or minimised. Acceptance: at least one explicit prediction of the form "varying $\theta$ from $0$ to $\pi/2$ at fixed drive amplitude changes the coherence-loss contribution by $\Delta T_2 = \ldots$" — quantitative + experimentally measurable on IBM Quantum hardware.
   - **Phase C** — IBM Quantum test. Specify the IBM Quantum Network resources required (qubit count, gate calibration access, microwave-pulse-shape control). Comparison with standard CPMG / Uhrig dynamical-decoupling sequences and with robust DCG variants. Acceptance: ✅ framework predicts measured coherence dependence on $\theta$/shape; ⚠ framework matches at precision floor only; 🔴 framework refuted.

4. **Page 3 — methods + asks + honest scope.**
   - Methods table: standard coherence-extension techniques (charge-noise filtering; optimised pulse shaping; tunable couplers) alongside what the framework adds *if* Phase A succeeds (angle/shape-modulated drive geometry).
   - Encoding-sequence direction: specific candidate sequences for further development (proper-time-spaced CPMG analogues; dissipative-Maxwell-motivated robust DCG variants).
   - Asks (audience-dependent — see issue #103 Audience options): IBM Quantum compute time; a superconducting-qubit theorist for Phase-A vetting; possibly a student or theorist familiar with the dual-theory framework (Tepper Gill or one of his students).
   - **Honest-scope paragraph (load-bearing — non-negotiable):** the white paper proposes a falsifiable research arc; it does not claim a derived prediction. The most likely outcome at Phase A is null (the dissipative term reduces to a textbook-equivalent contribution in the transmon limit), in which case the program closes honestly. The paper's value is **specifying a falsifiable test that does not currently exist in the campaign's repertoire**, not asserting that the framework will win that test. Mirror the structure of the [`2026-05_collaborator_intro_spectroscopy.md`](../../Roadmapping/Author_Reports/2026-05_collaborator_intro_spectroscopy.md) honest-scope paragraph.

5. **Crocco compliance markers.** Embed at minimum the following `<!-- TODO: human reviews and fills in -->` blocks:
   - One on the apparatus-mapping paragraph (Page 1) — substantive AI move connecting Eq. (4) to transmon-qubit drive physics.
   - One on the Phase B prediction paragraph (Page 2) — substantive AI move proposing a specific quantitative experimental signature.
   - One on the audience choice (Page 3) — substantive AI default of *internal-campaign planning*; Trey must confirm or override.
   - One terminal block summarising the open substantive-AI items for human sign-off before sharing the paper externally.

6. **Build the report.** From `Roadmapping/Author_Reports/`:

    ```bash
    ./build_report.sh 2026-05_quantum_coherence_whitepaper
    ```

    Verify build success and `pages: 3` exactly. If the page count is `2` (under) or `4` (over), trim or expand a specific paragraph; do not pad with filler. If the build fails on Unicode (likely culprits: `⁺`, `⇒`, em-dash `—`), substitute LaTeX math-mode equivalents (`$^{2+}$`, `$\Rightarrow$`) — these have already broken prior reports (PR #68, PR #89) and the same pattern applies here.

7. **Cross-link to the campaign.** Add explicit issue/PR cross-references in the body to: [#50](https://github.com/temoTxt/PyPhysics/issues/50) (Bethe-Salpeter source-of-record), [#78](https://github.com/temoTxt/PyPhysics/issues/78) (Li$^{2+}$ Z-extension), [#82](https://github.com/temoTxt/PyPhysics/issues/82) (Z-scan Outcome C), [#55](https://github.com/temoTxt/PyPhysics/issues/55) (proper-time one-loop, open frontier), [#75](https://github.com/temoTxt/PyPhysics/issues/75) (framework-specs, Tepper engagement), and Eq. (4) of the dual Maxwell verification doc at line 188.

8. **Commit and PR.** New commit on `103-ibm-quantum-coherence-white-paper`. PR title `103: white paper — proper-time theory on IBM quantum computer coherence (microwave-coupling angle/shape focus)`. PR body summarises the three pages, the load-bearing honest-scope paragraph, and the open audience-confirmation TODO. Does *not* close #103 (closes when the white paper is sent / acted on, per the audience choice).

## Files to Modify

| File | Change |
|---|---|
| `Roadmapping/Author_Reports/2026-05_quantum_coherence_whitepaper.md` | New file — 3-page white paper, microwave-coupling angle/shape focus, three-phase research arc, honest-scope paragraph, Crocco TODO blocks. |
| `Roadmapping/Author_Reports/2026-05_quantum_coherence_whitepaper.tex` | New file — auto-generated by `./build_report.sh 2026-05_quantum_coherence_whitepaper`. Implementer commits the generated file (matches the precedent of prior author reports in `Author_Reports/`). |
| `Roadmapping/Author_Reports/2026-05_quantum_coherence_whitepaper.pdf` | New file — auto-generated by `./build_report.sh`; **must be exactly 3 pages**. Committed as a binary artifact (matches prior author-report precedent). |

Nothing else is in scope. In particular: no edits to `Two_Mathematically_Equivalent_Versions_of_Maxwells_Equations.md`, `FINDINGS_for_author_review.md`, `10_CrossComparison.md`, `Dual_Relativistic_Quantum_Mechanics_I.md`, or any of the campaign's existing apparatus / verification docs. The white paper *references* those files but does not modify them — they are the source-of-record this proposal builds on, not the implementation surface.

## Dependencies

- `pandoc` (≥ 3.x), `texlive-latex-base`, `texlive-latex-extra`, `texlive-pictures`, `poppler-utils` — all already in the toolchain per `Roadmapping/Author_Reports/build_report.sh`'s dependency comment. No new system or Python dependencies.
- No new Wolfram MCP calls required at white-paper-writing time. The paper *describes* a Wolfram-MCP-derivable Phase A computation; deriving it is the follow-on work the white paper proposes, not part of this ticket.
- No new bibliography entries (`Roadmapping/History/Bibliography/`). External citations stay as inline DOI references in the report body, matching the precedent of prior author reports.

## Acceptance Criteria

- [ ] `Roadmapping/Author_Reports/2026-05_quantum_coherence_whitepaper.md` exists with YAML frontmatter (title / author / date / subject) matching the precedent of prior author reports in that directory.
- [ ] Page 1 names Eq. (4) of the dual Maxwell paper, cites `Two_Mathematically_Equivalent_Versions_of_Maxwells_Equations.md` line 232 for the extra-term statement and line 332 for the corresponding third-term in the modified Liénard–Wiechert, and states the falsifiable hypothesis about microwave-coupling angle/shape modulation of the framework-predicted noise.
- [ ] Page 2 contains all three phases (A apparatus mapping / B concrete prediction / C IBM Quantum test) with explicit acceptance criteria each; Phase A's null-result-most-likely closure condition is named explicitly.
- [ ] Page 3 covers methods + encoding-sequence directions + asks + honest-scope paragraph; honest-scope paragraph explicitly states *proposal, not derived prediction* and *null result at Phase A most likely*.
- [ ] PDF built via `./build_report.sh 2026-05_quantum_coherence_whitepaper` reports `[BUILD SUCCESS]` and `pages: 3` exactly. Re-build until the page count is exactly 3 (trim or expand a specific paragraph, never pad with filler).
- [ ] `.tex` file passes the `<!-- TODO` defensive check (the build script strips TODO comments before LaTeX compilation; verify no leakage).
- [ ] At least four `<!-- TODO: human reviews and fills in -->` blocks present in the `.md` per the Crocco-compliance-markers list above.
- [ ] Cross-references present in the body to issues #50 / #78 / #82 / #55 / #75 and to the dual Maxwell verification doc line 188.
- [ ] No edits to any file outside the three rows of the *Files to Modify* table.
- [ ] PR opened with title `103: white paper — proper-time theory on IBM quantum computer coherence (microwave-coupling angle/shape focus)`; PR does *not* close #103 (closes when the white paper is acted on per the audience choice).

## Testing

Reviewer-side verification commands:

```bash
# 1. Build the PDF from a fresh checkout, verify the page count
cd Roadmapping/Author_Reports
./build_report.sh 2026-05_quantum_coherence_whitepaper

# 2. Verify the page count is exactly 3
pdfinfo 2026-05_quantum_coherence_whitepaper.pdf | grep "^Pages:" | awk '{print $2}'
# Expected: 3

# 3. Verify the .tex has no TODO markers (the build script strips them; this is a regression check)
grep -c "TODO" 2026-05_quantum_coherence_whitepaper.tex
# Expected: 0

# 4. Verify the body cites the dual Maxwell extra term
grep -E "(u.*a.*b\^4|dissipative|extra term|Eq\\. ?\\(4\\))" 2026-05_quantum_coherence_whitepaper.md
# Expected: at least one match per pattern

# 5. Verify the four canonical cross-references are present
for ref in "#50" "#78" "#82" "#55" "#75"; do grep -c "$ref" 2026-05_quantum_coherence_whitepaper.md; done
# Expected: each grep returns ≥ 1
```

No new automated tests are added (the repo has no test framework per [`CLAUDE.md`](../../CLAUDE.md); verification is by the PDF build pipeline and grep checks above). The Phase A theoretical derivation — the substantive Wolfram-MCP computation the paper *proposes* — is **out of scope** for this ticket; the white paper specifies the work, it does not perform it.
