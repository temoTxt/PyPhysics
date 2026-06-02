# Task 81: Perhelion of Mercury

## Objective

Issue [#81](https://github.com/temoTxt/PyPhysics/issues/81) attaches three PDFs (Van Flandern 1976, Corda's "Secret", and Gill's *Dual Newtonian Theory*) and the issue body is empty otherwise. The implementer's task: (1) verify the math in Gill's paper end-to-end via Wolfram MCP from the equation chain at lines `\lb{2.4}` → `\lb{h4}` → `\lb{h7}` → `\lb{h8}` of [`/home/tmorris/Downloads/Dual Newtonian Theory.tex`](file:///home/tmorris/Downloads/Dual%20Newtonian%20Theory.tex); (2) finish the numerical computation of Mercury's perihelion advance under the three predictions the paper offers (Corda, full dual, approximate dual); (3) compare against the standard GR prediction `6πGM/[c²a(1−e²)]` and the observed value 43″/century; (4) flag any insights the paper's gravitational-side dynamics provide for the four open questions in [`Roadmapping/Author_Reports/2026-05_gps_relativity_summary_for_gill.md` §5](../../Roadmapping/Author_Reports/2026-05_gps_relativity_summary_for_gill.md) — particularly Q1 (curved-spacetime extension of the framework).

The paper's Eq. (h4) `a_m = −(GM/r²)(1 − GM/(c²r))·ê_r` is a candidate framework prediction for the gravitational dynamics that GPS PR #73 explicitly flagged as missing from the verified corpus. If the math reproduces the observed perihelion advance, this paper is the framework's GR-extension answer — directly answering GPS Q1 and potentially closing several follow-on threads.

## Background

The Mercury .tex source is currently *stashed* (see `git stash list`): pre-triage, the implementer committed it to `Roadmapping/Tepper_Gill_Papers/Dual Newtonian Theory.tex` and ran `pandoc -f latex -t markdown` to produce `Roadmapping/Converted_Markdown/Dual Newtonian Theory/Dual Newtonian Theory.md` (279 lines, no warnings, all 13 `\label{...}` equation IDs preserved). The triage step requires a clean working tree, so the files were stashed; the implementer should `git stash pop` (or `git stash apply`) on the new branch to restore them.

Key equations in the paper (from inspection of the .tex source):

- **Eq. (2.4)** (line 227 in the .tex): the proper-time Hamiltonian `K = H²/(2mc²) + mc²/2` — identical to DRQM I Eq. (I.6) and Maxwell-paper Eq. (16); verified ✅ in [`Roadmapping/Equation_Verification/Dual_Relativistic_Quantum_Mechanics_I.md` line 124](../../Roadmapping/Equation_Verification/Dual_Relativistic_Quantum_Mechanics_I.md).
- **Eq. (after 2.4)** (line 241): interaction Hamiltonian `H₁ = √(c²π² + (mc² + V)²)` chosen over `H₂` for the Newtonian-Maxwell unification; gives `K = π²/(2m) + mc² + V + V²/(2mc²)`. This is structurally identical to the DRQM I (II.3) "potential-in-the-mass" kernel verified ✅ at [`Roadmapping/Equation_Verification/Dual_Relativistic_Quantum_Mechanics_I.md` line 41](../../Roadmapping/Equation_Verification/Dual_Relativistic_Quantum_Mechanics_I.md).
- **Eq. (h3)** (line 249): `a = −(∇V/m)(1 + V/(mc²))` — the relativistic-Newtonian acceleration with the framework's `V²` correction.
- **Eq. (h4)** (line 257): two-body equations `a_m = −(GM/r²)(1 − GM/(c²r))·ê_r`, `a_s = (Gm/r²)(1 − Gm/(c²r))·ê_r`. Applies the Eq. (h3) form with gravitational `V = −GMm/r`.
- **Eq. (h7)** (line 286): dual orbital period `T_d = T₀·[(1+m/M)(1 − G(M²+m²)/((M+m)c²r₀))]^(−1/2)`.
- **Eq. (h8)** (line 293): dual angular frequency `ω_d = ω₀·[(1+m/M)(1 − G(M²+m²)/((M+m)c²r₀))]^(1/2)`.
- **Three predictions** (lines 310–333): Corda `Δφ_c = πm/M ≈ 44.39″/century`; full dual `Δφ_d = 2π·[...]^(1/2) − 2π`; approximate dual `Δφ_{d₁}`.
- **Headline** (line 336): paper claims `Δφ_c = 44.39″/century` ≈ GR's 42.98″/century ≈ observed 43″/century.

The standard GR prediction comes from the Schwarzschild geodesic via:

```
Δφ_GR = 6πGM_⊙ / [c² a (1 − e²)]   per orbit
```

For Mercury: `a = 5.7909 × 10¹⁰ m`, `e = 0.20563`, orbital period `T = 87.969` days. Over a tropical century (`100 yr ÷ T_orb`) this gives 42.98″/century — the standard GR figure.

The paper has **no Wolfram-MCP verification anywhere in the repository yet** — the corpus has 22 verified Gill papers but this one was added in this issue. The implementer is establishing the verification baseline.

**GPS cross-reference (load-bearing).** The proper-time companion campaign [PR #73 (merged at `ac19d44`)](https://github.com/temoTxt/PyPhysics/pull/73) and the follow-on author report at [`Roadmapping/Author_Reports/2026-05_gps_relativity_summary_for_gill.md`](../../Roadmapping/Author_Reports/2026-05_gps_relativity_summary_for_gill.md) flag four open questions about curved-spacetime extension of the dual framework. The Mercury paper proposes exactly such an extension: a *modified Newtonian* gravitational acceleration `a = −(GM/r²)(1 − GM/(c²r))·ê_r` derived from the framework's `K = π²/(2m) + V + V²/(2mc²)` kernel with `V = −GMm/r`. This is **not** `c → b` in the Schwarzschild metric (which was the speculative hypothesis adopted in GPS PR #73 §pt_01 §0); it is a different extension entirely (the framework's `V²/(2mc²)` term *applied to gravity*). The implementer must record whether this paper's extension reproduces GR predictions (Mercury perihelion is the canonical test) and what that says about the four GPS open questions.

The two reference PDFs attached to #81 (Van Flandern 1976, Corda's "Secret") supply non-Gill context. Van Flandern (1976) gives an independent classical-mechanics derivation of perihelion advance; Corda's paper produces the 44.39″/century value the Gill paper cites at line 315. The implementer ingests these as supporting bibliography but does not re-verify them in the verification campaign.

## Implementation Plan

1. **Restore the Mercury files** that were stashed during triage: `git stash list` to find the relevant stash, `git stash apply stash@{N}` (or `git stash pop` if it's at top). Confirm `Roadmapping/Tepper_Gill_Papers/Dual Newtonian Theory.tex` and `Roadmapping/Converted_Markdown/Dual Newtonian Theory/Dual Newtonian Theory.md` are back in the working tree.

2. **Ingest the two reference PDFs.** Download `Van-Flandern-1.pdf` and `Secret.C.Corda.pdf` from issue #81's attachments to `Roadmapping/Historical_Papers/Retrospective/`. Force-add via `git add -f` if needed per the History campaign's PDF-acquisition policy (see [`Roadmapping/History/Bibliography/README.md`](../../Roadmapping/History/Bibliography/README.md)).

3. **Scaffold bibliography stubs** for both reference papers under `Roadmapping/History/Bibliography/Retrospective/` using `Roadmapping/History/Bibliography/scaffold_bib_note.py` — `vanflandern1976_perihelion.md` and `corda_secret_perihelion.md`. Cite-keys follow the firstauthor+year+slug convention per [`Roadmapping/History/Bibliography/README.md`](../../Roadmapping/History/Bibliography/README.md). Also scaffold `gill_relativistic_newtonian.md` for the Gill paper itself (note: the .tex title is "Relativistic Newtonian Theory"; the user-supplied filename is "Dual Newtonian Theory" — confirm which name to use as cite-key).

4. **Create the verification document** at `Roadmapping/Equation_Verification/Dual_Newtonian_Theory.md` following the per-paper template established by the other 12 Equation_Verification docs. Per-equation entries for the 13 `\label{...}` equation IDs in the paper. Each entry: "As printed" (LaTeX + line ref), Wolfram MCP check (single-line code per CLAUDE.md gotchas), expanded derivation, standard-equation comparison (Goldstein §3.10 for the Kepler problem; MTW §40.5 for GR Mercury), verdict (✅/⚠/🔴).

5. **Create the campaign folder** `Roadmapping/Mercury_Perihelion/` with the structure below. Each per-effect doc mirrors the campaign-template convention from [`Roadmapping/GPS_Relativity/_template_effect.md`](../../Roadmapping/GPS_Relativity/_template_effect.md). The README has the standard status table + scope + parameter table (Mercury orbital constants + solar mass).

6. **Per-effect derivations and Wolfram MCP verifications** (the math work):

   - `01_dual_newtonian_setup.md` — derive the chain `H = √(c²p² + m²c⁴)` → `K = H²/(2mc²) + mc²/2` → with `H₁` interaction → `K = π²/(2m) + mc² + V + V²/(2mc²)`. Verify the `V²` term comes from the `(mc² + V)²` cross-term per the paper's `H₁` choice. Reference DRQM I §I.6 + §II.3 verified ✅ for the parent identities.
   - `02_two_body_equations_motion.md` — derive Eq. (h3) `a = −(∇V/m)(1 + V/(mc²))` from `K`; apply with `V = −GMm/r` to recover Eq. (h4) `a_m = −(GM/r²)(1 − GM/(c²r))·ê_r`. Wolfram MCP check via symbolic differentiation.
   - `03_orbital_period_and_omega.md` — derive Eq. (h5)–(h8) from the centripetal force balance for circular orbit at radius `r₀`. Wolfram MCP: substitute `M = M_⊙ = 1.989e30 kg`, `m = m_Mercury = 3.301e23 kg`, `r₀ = 5.791e10 m` (Mercury semi-major axis), `c = 2.998e8 m/s`, `G = 6.674e-11 N·m²/kg²`. Compute `T_d`, `ω_d` numerically and compare to `T₀`, `ω₀`.
   - `04_perihelion_advance.md` — compute the three predicted `Δφ` values: Corda `Δφ_c = πm/M`, full dual `Δφ_d`, approximate dual `Δφ_{d₁}`. Convert from radians-per-revolution to arcsec-per-tropical-century (multiply by (180/π)·3600·(centuries_per_tropical_century)·(orbits_per_tropical_century)). Mercury makes ≈ 414.93 orbits/century. Verify the paper's claimed 44.39″/century for Corda's value.
   - `05_comparison_with_GR.md` — compute the standard GR prediction `Δφ_GR = 6πGM/[c²a(1−e²)]` for Mercury (account for eccentricity `e = 0.20563`); produce `42.98″/century`. Tabulate: Corda (paper-quoted 44.39); full dual (compute); approximate dual (compute); GR (42.98); observed (43″). State which is closest and within what precision.

7. **Update the GPS author report (or create follow-on document) to flag GPS Q1 resolution.** If the Mercury paper's gravitational extension reproduces the observed perihelion advance within the framework's precision floor, that is a *concrete* answer to GPS Q1 (the curved-spacetime extension is "the framework's `V²/(2mc²)` kernel applied to gravity via `V = −GMm/r`", not `c → b` in the Schwarzschild metric). Add a closing-pointer paragraph to [`Roadmapping/Author_Reports/2026-05_gps_relativity_summary_for_gill.md` §6](../../Roadmapping/Author_Reports/2026-05_gps_relativity_summary_for_gill.md) noting the cross-reference, and re-build the PDF via `Roadmapping/Author_Reports/build_report.sh`.

8. **FINDINGS update.** If the math reveals any discrepancies in the paper (sign typos, factor errors, missing terms), record in `Roadmapping/Equation_Verification/FINDINGS_for_author_review.md` as a new Finding 4. If everything reproduces cleanly, record a "no findings" entry noting the paper passed verification.

## Files to Modify

| File | Change |
|---|---|
| `Roadmapping/Tepper_Gill_Papers/Dual Newtonian Theory.tex` | restore from stash; commit |
| `Roadmapping/Converted_Markdown/Dual Newtonian Theory/Dual Newtonian Theory.md` | restore from stash; commit |
| `Roadmapping/Historical_Papers/Retrospective/Van-Flandern-1.pdf` | download from issue #81; `git add -f` |
| `Roadmapping/Historical_Papers/Retrospective/Secret-Corda.pdf` | download from issue #81; `git add -f` |
| `Roadmapping/Historical_Converted_Markdown/Retrospective/Van-Flandern-1.md` | marker-pdf convert |
| `Roadmapping/Historical_Converted_Markdown/Retrospective/Secret-Corda.md` | marker-pdf convert |
| `Roadmapping/History/Bibliography/Retrospective/vanflandern1976_perihelion.md` | scaffold via `scaffold_bib_note.py` |
| `Roadmapping/History/Bibliography/Retrospective/corda_secret_perihelion.md` | scaffold |
| `Roadmapping/History/Bibliography/Retrospective/gill_relativistic_newtonian.md` | scaffold |
| `Roadmapping/Equation_Verification/Dual_Newtonian_Theory.md` | new — per-equation verification doc |
| `Roadmapping/Equation_Verification/FINDINGS_for_author_review.md` | append Finding 4 (or no-findings note) |
| `Roadmapping/Mercury_Perihelion/README.md` | new — campaign README |
| `Roadmapping/Mercury_Perihelion/01_dual_newtonian_setup.md` | new |
| `Roadmapping/Mercury_Perihelion/02_two_body_equations_motion.md` | new |
| `Roadmapping/Mercury_Perihelion/03_orbital_period_and_omega.md` | new |
| `Roadmapping/Mercury_Perihelion/04_perihelion_advance.md` | new |
| `Roadmapping/Mercury_Perihelion/05_comparison_with_GR.md` | new |
| `Roadmapping/Author_Reports/2026-05_gps_relativity_summary_for_gill.md` | append closing pointer if Q1 is answered; re-build PDF |
| `Roadmapping/Author_Reports/2026-05_gps_relativity_summary_for_gill.tex` | re-built artifact |
| `Roadmapping/Author_Reports/2026-05_gps_relativity_summary_for_gill.pdf` | re-built artifact |

## Dependencies

- **Wolfram MCP** (`mcp__wolfram__WolframLanguageEvaluator`) — for symbolic and numerical verification per CLAUDE.md gotchas (single-line code; avoid `V`/`e` as variable names — use `potV`/`ee`).
- **pandoc 3.1+** — for the .tex → .md conversion (already used to convert the Mercury paper).
- **marker-pdf** — for the two reference PDFs (`uv run python Roadmapping/parse_papers.py` with `--input` pointing at `Roadmapping/Historical_Papers/Retrospective/`).
- **build_report.sh** — for the author-report PDF rebuild after the §6 closing pointer is added.
- **Mercury orbital parameters** (constants, not dependencies): `a = 5.7909 × 10¹⁰ m`, `e = 0.20563`, `T_orb = 87.969 days`, `m_Mercury = 3.301 × 10²³ kg`. Solar mass `M_⊙ = 1.989 × 10³⁰ kg`. `G = 6.674 × 10⁻¹¹`. `c = 2.998 × 10⁸`. Centuries-per-tropical-century: `100·365.25 / 87.969 ≈ 415.20` orbits/century (use the paper's convention).

## Acceptance Criteria

- [ ] Mercury .tex + converted .md restored from stash and committed to the Tepper_Gill_Papers/ and Converted_Markdown/ folders.
- [ ] Both reference PDFs (Van Flandern, Corda) downloaded, marker-pdf converted, and bib stubs scaffolded (three bib stubs total: the two reference papers + the Gill paper itself).
- [ ] `Roadmapping/Equation_Verification/Dual_Newtonian_Theory.md` covers all 13 numbered equations from the .tex with Wolfram MCP checks.
- [ ] `Roadmapping/Mercury_Perihelion/` campaign folder created with README + 5 per-effect documents, mirroring the GPS_Relativity/ template structure.
- [ ] Eq. (h4) algebraic derivation from Eq. (h3) reproduced and Wolfram-MCP verified.
- [ ] Numerical computation of `Δφ_c`, `Δφ_d`, `Δφ_{d₁}` for Mercury with explicit arcsec-per-tropical-century conversion; paper's claimed 44.39″/century for Corda reproduced (or discrepancy flagged).
- [ ] Standard GR `Δφ_GR = 6πGM/[c²a(1−e²)]` computed for Mercury (target 42.98″/century).
- [ ] Comparison table: Corda / full dual / approximate dual / GR / observed.
- [ ] GPS author report §6 updated with a cross-reference paragraph if (and only if) the math reveals that the paper's gravitational extension is a concrete answer to GPS Q1. PDF rebuilt; page count remains in `[3, 7]`.
- [ ] FINDINGS Finding 4 (or "no findings" note) appended after the verification doc completes.
- [ ] All Wolfram MCP code blocks use the `ClearAll` discipline established in the GPS campaign's effect 08 fix (per the devil's advocate review of PR #71, see [PR #71 commit `210750d`](https://github.com/temoTxt/PyPhysics/pull/71/commits/210750d)).

## Testing

```bash
# Verify the Mercury paper's .tex still converts cleanly:
pandoc -f latex -t markdown \
  "Roadmapping/Tepper_Gill_Papers/Dual Newtonian Theory.tex" \
  -o /tmp/mercury-check.md
diff /tmp/mercury-check.md \
  "Roadmapping/Converted_Markdown/Dual Newtonian Theory/Dual Newtonian Theory.md"
# (should be byte-identical; mismatches reveal stale conversion)

# Re-run the Wolfram MCP checks in the verification doc end-to-end from a
# clean Wolfram session — each ClearAll-headed block should reproduce its
# claimed result.

# Re-build the author report after the §6 closing pointer is added:
cd Roadmapping/Author_Reports && ./build_report.sh 2026-05_gps_relativity_summary_for_gill --dry-run
# (defensive checks: no TODO leakage; page count in [3, 7])

# Confirm Mercury perihelion numerical match:
#   Paper claim:   Δφ_c (Corda)        = 44.39″/century
#   GR prediction: Δφ_GR               = 42.98″/century
#   Observed:                          = 43″   /century
#   Acceptable verdict for the campaign: dual prediction within ~5% of observed.
```

Reviewer should also confirm that the GPS Q1 cross-reference in the author report (if added) is **honest** — the Mercury paper's gravitational extension is a candidate answer to Q1, not a verified one, until §III review is performed by the author per the existing Crocco-compliance discipline.
