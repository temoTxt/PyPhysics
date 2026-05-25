# Proper-time third-term predictions for MeV bremsstrahlung

This document specialises the canonical proper-time bremsstrahlung derivation of [[Ch14_Radiation_by_Moving_Charges#Problem J3e-P14.6 — Bremsstrahlung from a linearly decelerating charge]] to the clinical-linac MeV regime (1, 6, 18 MeV) and compares the proper-time prediction against the published Born-approximation bremsstrahlung literature: the cross-section reviews of [[koch1959_bremsstrahlung_review]] and [[tseng1971_screened_bremsstrahlung]], the canonical tabulation of [[seltzer1986_bremsstrahlung_tables]], the Monte Carlo reference implementations [[kawrakow2000_egsnrc]] and [[salvat2018_penelope]], and the clinical-dosimetry calibration protocol [[aapm_tg51_1999_clinical_dosimetry]]. The point of the document is to determine whether the order-unity proper-time correction the Ch. 14 P14.6 problem predicts in the MeV regime is detectable against the percent-level precision floor of the published bremsstrahlung literature.

The document is structured as follows. §1 records the selection provenance and the framework's scope. §2 wikilinks to the canonical Ch. 14 P14.6 derivation and summarises the load-bearing scaling result. §3 takes 1 MeV, 6 MeV, and 18 MeV in turn, recording the kinematic setup, the per-energy proper-time prediction, and the verdict against the Born-approximation comparator. §4 collects the cross-energy summary. §5–§7 record Crocco compliance, the relationship to the now-shipped PR D + the Ch. 15 bremsstrahlung problems, and the conditional FINDINGS cross-post.

**Status note (operational).** This document is a draft. The substantive content — selection provenance in §1, per-energy quantitative predictions and verdicts in §3, and the cross-energy summary verdict column in §4 — is gated by `<!-- TODO: human reviews and fills in -->` blocks per the campaign-wide convention recorded in [§13.5 D2 of the #42 plan](../../../../.dev/tasks/42-electromagnetism-jackson-proper-time.md#135-decision-points--confirmed-by-author-2026-05-24). The pragmatic content — §2 wikilink-and-summary block, §3 kinematic setup paragraphs, §5/§6/§7 operational sections — is filled in here without gating, per [`CROCCO_COMPLIANCE.md`](../../../Tooling/CROCCO_COMPLIANCE.md) §3 (expectations 2 and 6).

**Predecessor work.** The full proper-time Liénard–Wiechert derivation, including the explicit field expressions of Eq. (7) of [[Two_Mathematically_Equivalent_Versions_of_Maxwells_Equations]] and the verification that the third term engages non-trivially in linear-acceleration geometries, lives at [[Ch14_Radiation_by_Moving_Charges#Problem J3e-P14.6 — Bremsstrahlung from a linearly decelerating charge]]. The Ch. 14 work uniquely establishes the `|E_3rd|/|E_acc,classical| ~ (u/c)^2` scaling at non-relativistic intensity and the order-unity contribution at MeV — both of which §2 below restates compactly.

---

## 1. Selection provenance

<!-- TODO: human reviews and fills in -->

*Required content for this section (per #42 plan §6 substantive-block convention):*

- Why the MeV-clinical regime specifically, vs the alternative of (i) staying at GeV with the #43 experiments, (ii) targeting non-relativistic atomic-physics bremsstrahlung (sub-keV), or (iii) waiting for next-generation precision-XFEL bremsstrahlung benchmarks. The Ch. 14 P14.6 author-review note ("Worth flagging as a candidate for issue #43's follow-on work") records the *intuition*; this block records the specific reasoning for prioritising MeV bremsstrahlung in *this* document.
- Why the three energies 1, 6, 18 MeV specifically. The issue calls these "typical of clinical linacs"; flesh out whether they span the relevant clinical-physics range adequately, whether 25 MeV (high-end clinical) or 0.5 MeV (sub-MeV) extensions are worth adding, and what the experimental-comparator precision floor is at each energy.
- Why the six bibliography candidates (cross-section reviews + tabulation + Monte Carlo implementations + clinical-dosimetry protocol) and not others. The Stage-1 selection was made under triage pressure; this block records the selection rationale that survives Trey's primary-source read.

---

## 2. The canonical proper-time bremsstrahlung framework

**See:** [[Ch14_Radiation_by_Moving_Charges#Problem J3e-P14.6 — Bremsstrahlung from a linearly decelerating charge]] for the full derivation.

For self-contained reference, the load-bearing result the Ch. 14 work establishes is the following. For a charge undergoing linear deceleration — `\mathbf{u}` and `\mathbf{a}` along the same line, so `\mathbf{u}\cdot\mathbf{a} = -|u||a|` — the proper-time Liénard–Wiechert third term of Eq. (7) of [[Two_Mathematically_Equivalent_Versions_of_Maxwells_Equations]],

$$\mathbf{E}_{3}(\mathbf{x},\tau) = \frac{e\,(\mathbf{u}\!\cdot\!\mathbf{a})\,[\mathbf{r}\!\times\!(\mathbf{u}\!\times\!\mathbf{r})]}{b^{4}\,s^{3}},$$

engages non-trivially. We observe that the ratio of the third-term radiated-field magnitude to the classical Liénard acceleration-field magnitude scales as `(u/c)^2` at non-relativistic intensity, growing to order unity at MeV electron energies, and dominating in the ultra-relativistic limit. The numerical table established in Ch. 14 P14.6 is reproduced compactly for orientation:

| Electron kinetic energy | `u/c` | `|\mathbf{E}_3|/|\mathbf{E}_{\text{acc,classical}}|` |
|---|---|---|
| Chemistry-scale | `\sim 10^{-2}` | `\sim 10^{-4}` |
| 10 keV | `\sim 0.20` | `\sim 4 \times 10^{-2}` |
| 100 keV | `\sim 0.55` | `\sim 0.3` |
| **1 MeV** | `\sim 0.94` | order unity |
| **6 MeV** | `\sim 0.997` | order unity (saturated) |
| **18 MeV** | `\sim 0.9996` | order unity (saturated) |
| Ultrarelativistic | `\to 1` | dominates |

The 1, 6, 18 MeV rows are the targets of §3 below. The classical-EM prediction is *purely transverse* radiation in the radiation zone; the proper-time third term carries an `\mathbf{r}\times(\mathbf{u}\times\mathbf{r})` geometric factor that, for `\mathbf{u}` parallel to `\mathbf{a}` (the linear-acceleration geometry), produces a non-zero projection along `\mathbf{u}` itself — a **longitudinal polarisation component**. This is the experimentally testable signature: classical EM forbids it, proper-time predicts it.

---

## 3. Per-energy comparison

Each subsection records the kinematic setup (pragmatic — derived from the Ch. 14 P14.6 scaling table and the bibliography geometry), a per-energy proper-time prediction setup (mixed), and a comparison + verdict against the Born-approximation comparators (substantive — gated TODO block).

### 3.1 1 MeV electron beam

**Kinematic setup (pragmatic):** at electron kinetic energy `T = 1 \text{ MeV}`, the Lorentz factor is `\gamma = 1 + T/(m_e c^2) \approx 2.96`, giving `\beta = v/c \approx 0.941` and `u/c = \gamma\beta \approx 2.79`. The proper-time scalar is `b = \sqrt{c^2 + u^2} = \gamma c \approx 2.96 c`. The scaling `|\mathbf{E}_3|/|\mathbf{E}_{\text{acc,classical}}| \sim (u/c)^2 / (\text{order-unity geometric factors})` evaluates to roughly `\sim 1` — the third term is comparable in magnitude to the classical acceleration field but not yet dominant.

**Proper-time prediction setup:** the linear-deceleration geometry of an electron entering a high-Z target follows the Ch. 14 P14.6 derivation directly. The decelerating force is approximately constant within a screening length of the nucleus (Bethe-Heitler regime), and the resulting `\mathbf{u}\cdot\mathbf{a}` is large and negative throughout the deceleration. We observe that for 1 MeV — the lower end of the clinical-linac range — the proper-time longitudinal-polarisation prediction should be most easily distinguishable from the classical purely-transverse prediction because the third-term contribution is comparable in magnitude to (rather than dominated by or dominating) the classical contribution. The cleanest experimental signature would be a polarisation-resolved measurement of bremsstrahlung emitted at angles `\theta \sim 1/\gamma \approx 20°` from the electron direction, where the standard relativistic-beaming pattern peaks.

**Quantitative prediction:**

<!-- TODO: human reviews and fills in -->

*Required content (substantive — Crocco §5):* compute the proper-time-augmented angular distribution `dW/dΩ` at 1 MeV, including the third-term longitudinal contribution. Compare numerically against the Born-approximation result tabulated by [[seltzer1986_bremsstrahlung_tables]] and implemented in [[kawrakow2000_egsnrc]] / [[salvat2018_penelope]]. The Mathematica MCP work lives in `Roadmapping/Mathematica_Notebooks/Electromagnetism/Bremsstrahlung_1MeV_proper_time.wl` (companion notebook to be created).

**Comparison and verdict:**

<!-- TODO: human reviews and fills in -->

*Required content (substantive — Crocco §5):* place the proper-time prediction on the same plot as the Born-approximation curve. Quote the chi-squared against any polarisation-resolved 1 MeV bremsstrahlung data in the literature (or note the absence of such data). Verdict: ✅ (consistent with current data) / ⚠ (distinguishable but within precision floor) / ❌ (ruled out by current data).

### 3.2 6 MeV electron beam

**Kinematic setup (pragmatic):** at electron kinetic energy `T = 6 \text{ MeV}` — the canonical clinical-linac low-energy mode — the Lorentz factor is `\gamma \approx 12.7`, with `\beta \approx 0.997`. The third-term-to-classical-acceleration ratio is solidly in the order-unity saturated regime (per the Ch. 14 P14.6 scaling table). The 6 MeV mode is the most extensively-characterised clinical-linac configuration in the medical-physics literature, with TG-51 dosimetry calibration ([[aapm_tg51_1999_clinical_dosimetry]]) operating directly at this energy.

**Proper-time prediction setup:** identical linear-deceleration geometry to §3.1; the difference is purely in `(u/c)` and `\gamma`. We observe that at 6 MeV the third-term and classical contributions are *both* in the radiation-zone regime where they can be compared at the level of the angular distribution. The relativistic-beaming half-angle has tightened to `1/\gamma \approx 4.5°`, so polarisation-resolved measurements need finer angular resolution than at 1 MeV, but the integrated cross-section measurements from clinical-dosimetry calibration data offer the highest precision benchmark.

**Quantitative prediction:**

<!-- TODO: human reviews and fills in -->

*Required content (substantive — Crocco §5):* same template as §3.1, evaluated at 6 MeV. Special attention to whether the integrated cross-section is modified at the percent level (TG-51 precision floor) or only at the angular-distribution / polarisation level. If the integrated cross-section is unchanged at the percent level — which the scaling argument suggests, since polarisation differences typically cancel in integration — then 6 MeV is the *consistency* test rather than the *discriminating* test. Companion notebook: `Bremsstrahlung_6MeV_proper_time.wl`.

**Comparison and verdict:**

<!-- TODO: human reviews and fills in -->

*Required content (substantive — Crocco §5):* same template as §3.1.

### 3.3 18 MeV electron beam

**Kinematic setup (pragmatic):** at electron kinetic energy `T = 18 \text{ MeV}` — the upper end of the standard clinical-linac range — the Lorentz factor is `\gamma \approx 36.2`, with `\beta \approx 0.9996`. The third-term-to-classical-acceleration ratio is deep in the order-unity saturated regime; the relativistic-beaming half-angle is `1/\gamma \approx 1.6°`. The 18 MeV configuration is used clinically for deep-seated tumour treatments and has historically been the upper-end calibration target for the AAPM dosimetry protocols.

**Proper-time prediction setup:** same as §3.1 and §3.2 in geometric structure; the third-term contribution is now significant but the classical Liénard `\gamma^6` scaling has also grown sharply, so the *ratio* is no longer the most sensitive signature. The polarisation-resolved signature is the cleaner discriminator at this energy, because the integrated power is dominated by the classical contribution's `\gamma^6` enhancement and any third-term correction at the order-unity ratio level is masked.

**Quantitative prediction:**

<!-- TODO: human reviews and fills in -->

*Required content (substantive — Crocco §5):* same template, evaluated at 18 MeV. Note that 18 MeV is approaching the regime where the precision Born-approximation literature thins out — most modern bremsstrahlung cross-section work is either at lower energies (clinical dosimetry) or at GeV (synchrotron / radiation-reaction experiments of #43). The 18 MeV verdict is therefore expected to lean toward "⚠ consistent within the published precision but at the edge of the experimental coverage." Companion notebook: `Bremsstrahlung_18MeV_proper_time.wl`.

**Comparison and verdict:**

<!-- TODO: human reviews and fills in -->

*Required content (substantive — Crocco §5):* same template as §3.1.

---

## 4. Cross-energy summary

<!-- TODO: human reviews and fills in -->

*Required content (substantive — Crocco §5):* fill in the verdict column after §3 verdicts are in. The published-comparator columns are pragmatic restatements of what each cited tabulation predicts at the indicated energy.

| Energy | Born-approximation (Koch-Motz / Tseng-Pratt / Seltzer-Berger) | Monte Carlo (EGSnrc / PENELOPE) | Proper-time (this document) |
|---|---|---|---|
| 1 MeV | tabulated; ~1–2% precision | implemented; matches tabulation to ~0.5% | ? |
| 6 MeV | tabulated; ~1% precision (TG-51 calibration energy) | implemented; matches to ~0.5% | ? |
| 18 MeV | tabulated; ~1–2% precision | implemented; matches to ~0.5% | ? |

The three energies span the clinical-linac range. The interesting question is whether the proper-time prediction's order-unity correction to the radiation field translates into a *measurable* deviation in the integrated cross-section (in which case the comparison is sharp and the precision floor decisive) or only into a polarisation-resolved angular-distribution signature (in which case the existing integrated-cross-section literature does not discriminate, and a dedicated polarisation experiment is needed). The §3 quantitative-prediction blocks need to settle this for each energy.

---

## 5. Crocco compliance

This document mixes pragmatic and substantive AI use per [`CROCCO_COMPLIANCE.md`](../../../Tooling/CROCCO_COMPLIANCE.md) and [#42's plan §6](../../../../.dev/tasks/42-electromagnetism-jackson-proper-time.md#6-source-material-treatment-copyright--crocco-compliance):

- **Pragmatic** (no separate disclosure block in the body): the §2 wikilink-and-summary block (restated compactly from the verified Ch. 14 P14.6 derivation); the §3 kinematic-setup paragraphs (mechanical Lorentz-factor arithmetic at the three indicated energies); the §3 proper-time-prediction-setup paragraphs (which framework features engage at each energy, derived from the Ch. 14 scaling table); the §4 table headers and published-comparator columns; the §1 introductory framing; this §5; the §6 PR-D-and-Ch.-15 cross-link section; the §7 FINDINGS pointer.
- **Substantive** (`<!-- TODO: human reviews and fills in -->` blocks throughout): §1 selection provenance; every §3 per-energy *quantitative prediction* and *verdict*; the §4 cross-energy summary table's proper-time-verdict column and the discriminating-vs-consistency-test narrative paragraph.

The substantive blocks are deliberately blank in this draft. They are gated by Zorns-Lemmon's review per the #43 routing convention (see §6 of [`.dev/tasks/48-bremsstrahlung-mev-spectra.md`](../../../../.dev/tasks/48-bremsstrahlung-mev-spectra.md)); Trey holds final responsibility for the posted artifact per Crocco rule 8.

No prompt-of-record was used to generate any of the pragmatic content above; it is all derived from the canonical Ch. 14 P14.6 work and the six bib-stub navigational summaries by direct re-statement and elementary kinematic arithmetic. If the substantive blocks are filled in by a future AI pass (rather than by Zorns or Trey directly), the pass MUST use a committed prompt under [`Roadmapping/Tooling/_prompts/`](../../../Tooling/_prompts/) — candidate name `derive_proper_time_prediction.md` (new; the same prompt #43's §5 floats).

---

## 6. Relationship to Ch. 14 P14.6 and Ch. 15

The proper-time bremsstrahlung framework this document specialises lives at [[Ch14_Radiation_by_Moving_Charges#Problem J3e-P14.6 — Bremsstrahlung from a linearly decelerating charge]] (PR D, shipped 2026-05-24). The closely-related QED-adjacent bremsstrahlung problems — full Born-approximation cross-section, Bethe-Heitler formula, method-of-virtual-quanta treatments — live in [[Ch15_Bremsstrahlung_Virtual_Quanta]] (PR J, shipped). When this document's §3 verdicts land, both Ch. 14 P14.6 and Ch. 15 should be edited to add forward-pointers; the `chapters_citing` field of each of the six bib stubs scaffolded in Stage 1 should also be updated to record that this prediction document cites them. None of those bidirectional cross-links land in this PR; they happen in follow-up commits after Zorns engages.

---

## 7. FINDINGS cross-post (conditional, deferred)

If any §3 verdict comes back ⚠ or ❌ when Zorns fills the substantive blocks, the corresponding finding is recorded in [`Roadmapping/Equation_Verification/FINDINGS_for_author_review.md`](../../../Equation_Verification/FINDINGS_for_author_review.md) following the pattern established by [`Roadmapping/Electromagnetism/Jackson/Experimental_Comparisons/radiation_reaction_2018.md`](radiation_reaction_2018.md)'s §7 — predictions-vs-experiment findings are structurally different from the FINDINGS doc's existing three internal-derivation entries, and the format extension referenced there applies equally here. No cross-post happens in this PR; the conditional cross-post is a small follow-up PR after the §3 verdicts are recorded.
