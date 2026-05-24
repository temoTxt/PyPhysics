# Proper-time predictions for the 2018 radiation-reaction experiments

This document compares the predictions of the proper-time reformulation of Maxwell's equations — the framework verified in [[Two_Mathematically_Equivalent_Versions_of_Maxwells_Equations]] — against the three 2018 experimental measurements of radiation reaction enumerated in [issue #43](https://github.com/temoTxt/PyPhysics/issues/43): Cole et al. (Gemini facility, *PRX* **8**, 011020), Poder et al. (Gemini facility, *PRX* **8**, 031004), and Wistisen et al. (CERN SPS, *Nature Commun.* **9**, 795). Each experiment lies in a regime where the proper-time framework's `(u·a)/b⁴` third term in the modified Liénard–Wiechert fields (Eq. 7 of the Maxwell paper) is predicted to contribute, and where the classical theory's published prediction is itself disputed between competitors.

The document is structured as follows. §1 records the selection provenance and the framework's scope. §2 restates the proper-time third term in compact form and identifies the kinematic conditions under which it contributes. §3 takes each of the three experiments in turn, recording the geometry, setting up the proper-time prediction, and recording the verdict against the published data. §4 collects the cross-experiment summary table. §5–§7 record Crocco compliance, the work blocked on PR D of [#42](https://github.com/temoTxt/PyPhysics/issues/42), and the conditional FINDINGS cross-post.

**Status note (operational).** This document is a draft. The substantive content — selection provenance in §1, the physical-interpretation paragraph in §2.3, every per-experiment quantitative prediction and verdict in §3, and the cross-experiment summary in §4 — is gated by `<!-- TODO: human reviews and fills in -->` blocks per the campaign-wide convention recorded in [§13.5 D2 of the #42 plan](../../../../.dev/tasks/42-electromagnetism-jackson-proper-time.md#135-decision-points--confirmed-by-author-2026-05-24). The pragmatic content — §2.1 form-level restatement of the third term, the geometric setup paragraphs in §3, §5/§6/§7 operational sections — is filled in here without gating, per [`CROCCO_COMPLIANCE.md`](../../../Tooling/CROCCO_COMPLIANCE.md) §3 (expectations 2 and 6) and the campaign-wide pragmatic/substantive split in [#42's plan §6](../../../../.dev/tasks/42-electromagnetism-jackson-proper-time.md#6-source-material-treatment-copyright--crocco-compliance).

**Predecessor work.** The proper-time Liénard–Wiechert E and B fields used throughout §2 and §3 are verified in [[Two_Mathematically_Equivalent_Versions_of_Maxwells_Equations#Eq. (7)]] via two diagnostics (the third term vanishes when `a = 0`; the velocity term agrees exactly with Jackson's textbook form in 1D under the dictionary `β = u/b`). The companion dissipative coefficient `−(u·a)/b⁴ · ∂_τ` in the proper-time wave equation is verified at [[Two_Mathematically_Equivalent_Versions_of_Maxwells_Equations#Eq. (4)]]; the two manifestations of the same chain-rule effect are recorded as a single mechanism at two stages of the derivation in [§1 of the #42 plan](../../../../.dev/tasks/42-electromagnetism-jackson-proper-time.md#1-goal-and-scope).

**Predecessor scope.** The canonical Ch. 14 Liénard–Wiechert derivation in per-problem form lives at `Ch14_Radiation_by_Moving_Charges.md` and is scheduled for PR D of [#42](https://github.com/temoTxt/PyPhysics/issues/42); the document does not exist as of this PR. §2.1 below therefore re-quotes Eq. (7) self-contained rather than wikilinking the (not-yet-existing) Ch. 14 canonical derivation; §6 records the reconciliation step that runs when PR D lands.

---

## 1. Selection provenance

<!-- TODO: human reviews and fills in -->

*Required content (substantive — Crocco §5):*

- *Chosen because:* the three 2018 experiments together span the two regimes where the framework's third term is expected to contribute differently (head-on laser–electron collision, where the `u·a` cycle-average is dominated by the radiative back-action on the longitudinal velocity; and high-energy crystal channeling, where `u·a` is suppressed at leading order by the transverse character of the channeling acceleration but reappears at second order through the cumulative longitudinal energy loss). Cole and Poder share the Gemini-facility geometry and were therefore included as a pair, so that the framework's prediction is compared not only against the published spectra but against the two independent analysis pipelines applied to the same physical event class.
- *Alternatives considered:* the JLab E-144 series (older, less constrained by modern strong-field-QED comparators); BELLA-PW data published 2024 (preliminary at the time the campaign was scoped); the planned LUXE experiment at DESY (not yet a published measurement).
- *Role in the campaign:* this document is the proper-time framework's first experimental-comparison artifact. The verdicts recorded here become inputs to the Forward-chapter discussions on radiation-reaction in the dual-theory framework, and any ⚠ or ❌ verdict feeds [`FINDINGS_for_author_review.md`](../../../Equation_Verification/FINDINGS_for_author_review.md) per [§7 below](#7-findings-cross-post-conditional-deferred).

---

## 2. The proper-time Liénard–Wiechert third term

### 2.1 The field expressions

We restate the verified Eq. (7) of [[Two_Mathematically_Equivalent_Versions_of_Maxwells_Equations]] in compact form for self-contained reference. With

- `r = |x − x̄|` the distance from the source's retarded position to the field point,
- `b = √(c² + u²)` the dual proper-time scalar of Eq. (1) (distinct from `c`),
- `s = r − (r · u)/b` the proper-time analogue of the standard retardation factor `κR = R(1 − β)`,
- `r_u = r − (r/b) u` the proper-time analogue of `(n − β)`,

the modified Liénard–Wiechert E-field decomposes into three terms:

$$\mathbf{E}(\mathbf{x},\tau) = \underbrace{\frac{e\,[\mathbf{r_u}(1 - \mathbf{u}^2/b^2)]}{s^3}}_{\text{(I) velocity field}} \;+\; \underbrace{\frac{e\,[\mathbf{r}\!\times\!(\mathbf{r_u}\!\times\!\mathbf{a})]}{b^2 s^3}}_{\text{(II) acceleration field}} \;+\; \underbrace{\frac{e\,(\mathbf{u}\!\cdot\!\mathbf{a})\,[\mathbf{r}\!\times\!(\mathbf{u}\!\times\!\mathbf{r})]}{b^4 s^3}}_{\text{(III)}}.$$

Terms (I) and (II) recover the standard textbook Liénard–Wiechert fields (Jackson 3e §14.1, Eq. 14.14) under the dictionary `β = u/b`, `R ↔ r`, `κR ↔ s`. The two formulations are *mathematically equivalent but not (what we mean by) physically equivalent*, in the sense established by the Gill–Zachary Maxwell paper: term (III) has no counterpart in standard Liénard–Wiechert and carries the scalar prefactor `(u·a)/b⁴`. The B-field decomposes the same way (Eq. 7, lower line); its third term carries the same `(u·a)/b⁴` prefactor.

### 2.2 When the third term contributes

Term (III) vanishes whenever `u · a = 0`, that is, whenever the velocity and acceleration are everywhere perpendicular along the worldline. The zero is exact rather than leading-order, and is the diagnostic that pins down the conditional structure of the prediction. It follows that two physical situations make the third term contribute. First, any motion in which the kinetic energy is changing along the direction of motion produces a non-zero `u_∥ a_∥`. Second, a counter-propagating field that drives transverse acceleration `a_⊥` while extracting longitudinal energy through radiative back-reaction produces a non-zero cycle-averaged `u_z a_z`, because the back-reaction lowers `u_z` over the pulse.

By contrast, purely transverse motion — a perfect circular orbit, or transverse channeling oscillation with no longitudinal energy loss — gives `u · a = 0` instantaneously and the third term decouples. The same kinematic structure that protects the standard Larmor formula from a radiation-reaction correction at the leading-order Lorentz-invariant level protects the proper-time formulation's third term from contributing at the same order in those configurations.

### 2.3 Scaling of the radiated-power contribution

The Poynting flux radiated by a localised charge is `S = (c/4π) E × B`. The total radiated power decomposes into pure-self terms `|E_I|²`, `|E_II|²`, `|E_III|²` and cross terms `2 E_I · E_III`, `2 E_II · E_III`. We observe that the pure-self contribution of term (III) to the radiated power per unit solid angle scales as

$$\frac{dP_{\text{III}}}{d\Omega} \;\sim\; \frac{e^{2}\,(\mathbf{u}\!\cdot\!\mathbf{a})^{2}}{b^{8}} \cdot \frac{r^{2}\,|\mathbf{r}\!\times\!(\mathbf{u}\!\times\!\mathbf{r})|^{2}}{s^{6}}.$$

This contribution does not appear in the standard Lorentz–Abraham–Dirac or Landau–Lifshitz expressions for the radiated power, which carry `a²` and (relativistically) `γ⁴ a²` dependencies but not `(u·a)²/b⁸`. The cross terms `2 E_I · E_III` and `2 E_II · E_III` provide additional contributions that, depending on the geometry, may dominate or be sub-leading relative to the pure-self term.

<!-- TODO: human reviews and fills in -->

*Substantive interpretation block (Crocco §5):*

- The role of the `(u·a)/b⁴` prefactor as a candidate dissolution of the Abraham–Lorentz pathology — and how this connects to the dissipative coefficient `−(u·a)/b⁴ · ∂_τ` appearing earlier in the derivation chain (Eq. 4 of the Maxwell paper) as the *same mechanism at an earlier stage* rather than an independent term — needs to be written by Trey, since it is the framework's interpretive claim and its placement against the Lorentz–Abraham–Dirac and Landau–Lifshitz literatures is judgement-load. The two diagnostics in [[Two_Mathematically_Equivalent_Versions_of_Maxwells_Equations#Eq. (7)]] (vanishing when `a = 0`; exact agreement with Jackson on the velocity term) certify the *form* of the third term but not its physical interpretation.

---

## 3. Per-experiment comparison

Each subsection records the experimental geometry (pragmatic, restated from the converted primary source), sets up the proper-time prediction (mixed — kinematic computation is pragmatic, identifying which framework features contribute is substantive), and records the comparison against the published spectra with the verdict (substantive — gated TODO block).

### 3.1 Cole et al. 2018 — *PRX* **8**, 011020

[[cole2018_radiation_reaction]] · PDF (local) · [converted markdown](../../../Historical_Converted_Markdown/Retrospective/cole2018_radiation_reaction/cole2018_radiation_reaction.md)

**Geometry (pragmatic; restated from the primary source):**

The experiment was conducted at the Astra-Gemini laser of the Central Laser Facility, Rutherford Appleton Laboratory, UK. Two synchronised linearly polarised 800 nm beams of 45 fs FWHM duration were used in a counter-propagating geometry. The first beam, focused with an f/40 spherical mirror to a 37 × 49 μm focal spot, delivered (8.6 ± 0.6) J at peak intensity (7.7 ± 0.4) × 10¹⁸ W cm⁻², peak normalised amplitude `a_0 = 1.9 ± 0.1`, into a 15 mm helium gas jet at peak electron density (3.7 ± 0.4) × 10¹⁸ cm⁻³; the resulting laser-wakefield-accelerated electron beam was reported with mean energy (550 ± 20) MeV in the reference (no-collision) shots, with a typical spectrum spanning 0.25 – 2 GeV recorded on a Lanex screen via magnetic dispersion. The second beam, focused at the rear edge of the gas jet with an off-axis f/2 parabolic mirror to a 2.4 × 2.8 μm focal spot, delivered (10.0 ± 0.6) J at peak intensity (1.3 ± 0.1) × 10²¹ W cm⁻², nominal peak amplitude `a_0 = 24.7 ± 0.7`; the effective interaction amplitude after accounting for spatio-temporal jitter between the focal spot and the electron bunch was estimated to be `a_0 = 12 ± 1`. The collision delivered four successful overlap shots out of eight attempts. The post-collision observables are the electron energy spectrum (ε_initial → ε_final shifted from (550 ± 20) MeV to (470 ± 10) MeV in the bright-overlap subset) and the bremsstrahlung-subtracted γ-ray spectrum recorded on a 33 × 47 CsI scintillator array, with best-fit critical energy ε_crit ≳ 30 MeV. The quantum nonlinearity parameter is reported as `η ≈ 0.07` at peak conditions. The paper's preferred radiation-reaction comparator was a quantum (stochastic) model over the classical Landau–Lifshitz form, with the published ratio of γ-ray-inferred to electron-spectrum-inferred laser amplitudes `R = a_0(ε_crit)/a_0(ε_final)` close to unity for the quantum model (`R = 0.8 ⁺⁰·⁷₋₀.₃`) and biased low for the classical model (`R = 0.6 ⁺⁰·³₋₀.₂`); the preference for the quantum model is at approximately the 1σ level.

**Proper-time prediction setup:**

The Cole experiment is a counter-propagating geometry: the electron travels in the `−ẑ` direction with speed close to `c`, and the laser pulse propagates in the `+ẑ` direction with transverse oscillating fields that drive transverse acceleration `a_⊥` of magnitude set by the effective `a_0 ≈ 12`. We observe that the instantaneous `u · a` in such a geometry decomposes into the transverse term `u_z · a_⊥`, which vanishes by orthogonality, and the longitudinal term `u_z · a_z`, where `a_z` is the longitudinal deceleration induced by the radiative back-reaction. The cycle-averaged `u · a` is therefore controlled by the rate at which the electron is losing energy along its direction of motion. It is clear that the proper-time prediction must, in this geometry, *add to* the standard radiation-reaction loss rate rather than dominate it — the framework's third term sources additional radiation in proportion to the energy already being lost through the standard Larmor channel. The relatively small quantum parameter `η ≈ 0.07` places Cole's working point at the low end of the quantum regime, where the classical and quantum radiation-reaction predictions are themselves close (the paper reports the quantum model giving roughly 0.75 times the classical loss rate at this `η`).

**Quantitative prediction:**

<!-- TODO: human reviews and fills in -->

*Required content (substantive — Crocco §5):* using §2.1's third term, compute the proper-time-augmented energy-loss rate over the Cole collision pulse. Compare numerically to the published classical Landau–Lifshitz prediction and to the quantum-corrected Landau–Lifshitz prediction. The Mathematica MCP work for this lives in `Roadmapping/Mathematica_Notebooks/Electromagnetism/Cole2018_proper_time_prediction.wl` (to be created when this block is filled in).

**Comparison and verdict:**

<!-- TODO: human reviews and fills in -->

*Required content (substantive — Crocco §5):* place the proper-time prediction on the same plot as the published data, the classical L-L curve, and the quantum-corrected L-L curve that the paper preferred. Record the chi-squared (or equivalent goodness-of-fit) against the data within the published uncertainty band. Verdict: ✅ (prediction matches within bands) / ⚠ (prediction differs but within physically-defensible adjustment) / ❌ (prediction outside all uncertainty bands).

### 3.2 Poder et al. 2018 — *PRX* **8**, 031004

[[poder2018_radiation_reaction]] · PDF (local) · [converted markdown](../../../Historical_Converted_Markdown/Retrospective/poder2018_radiation_reaction/poder2018_radiation_reaction.md)

**Geometry (pragmatic; restated from the primary source):**

The Poder experiment used the same Astra-Gemini twin-beam facility as Cole, with an independent analysis pipeline and a higher-energy electron beam at a higher quantum nonlinearity parameter. The driver beam ((42 ± 3) fs FWHM, 9 J, f/40 focus FWHM ≈ 59 × 67 μm, peak `a_0 ≈ 1.7`, self-focused to peak `a_0 ≈ 3` in the plasma) was focused at the entrance of a 2 cm helium gas-cell at electron density 2 × 10¹⁸ cm⁻³, producing a laser-wakefield-accelerated electron beam with energy spectrum extending to approximately 2 GeV (`γ_e ≈ 4 × 10³`), beam source size `D_e ≲ 1 μm`, and 1 mrad divergence above 1 GeV. The scattering beam (twin pulse, (42 ± 3) fs, (8.8 ± 0.7) J after compression, f/2 off-axis parabola with concentric f/7 hole, focused 1 cm downstream of the gas-cell exit) had nominal peak amplitude `a_0 ≈ 22.5` at focus; the inherent lag of the wakefield-accelerated electron bunch behind the driver pulse (≈ 64 fs) brought the electrons to a defocused region of the scattering laser where the peak amplitude experienced by the bunch was reduced to `a_0 ≈ 12.3` with FWHM intensity profile of ≈ 7 μm. The reported quantum parameter is `χ ≈ 0.2`, roughly three times that of Cole. The observables are the post-collision electron energy spectrum recorded on a Lanex screen via a 1.0 T dipole spectrometer (350 MeV minimum, 7% energy resolution at 1.5 GeV) and the γ-ray beam on a caesium iodide scintillator. The reference (no-collision) electron spectrum shows a local peak at approximately 1.2 GeV; the perturbed spectrum (best-overlap shots) shows that local peak shifted to ≈ 1 GeV with near-total depletion of electrons above 1.6 GeV, consistent with a 30–40% energy loss for the high-energy tail. The paper compares four models against the perturbed spectrum: a Lorentz-only "perturbative" model (which strongly overestimates the loss), the Landau–Lifshitz equation (good agreement but slightly overestimates loss), a semiclassical model (LL weighted by `g(χ)`; improved agreement), and a multi-particle QED model (also good agreement, but with a "knee" structure not reproduced by the data). The paper's preferred comparator is the semiclassical model, with the QED model also acceptable; pure classical LL is "slightly" disfavoured at the experimental precision.

**Proper-time prediction setup:**

Since the Poder experiment uses the same Gemini-facility two-beam geometry as Cole, the proper-time prediction at the kinematic level coincides with §3.1's. The work of §3.2 is therefore both a consistency check (an independent analysis pipeline applied to the same physical event class) and a sharpening (Poder's higher `χ ≈ 0.2`, larger beam energy, and larger fractional energy loss provide finer resolution on which radiation-reaction comparator best fits). We observe that the energy-loss formula the Poder analysis uses to estimate `Δγ/γ` for a Gaussian plane-wave pulse,
`Δγ/γ ≈ (√(π/log 2) τ_0 t_L ω_L² γ_e a_0²/2) / (1 + same)` with `τ_0 = 2r_0/3c ≈ 6.3 × 10⁻²⁴ s`,
is the standard LL prediction; the proper-time framework's third-term contribution would multiply this baseline by a correction factor whose form is to be derived in the quantitative-prediction block below. The cycle-averaged `u · a` is again controlled by the rate of longitudinal energy loss, and the framework's third term again *adds to* the LL channel rather than replacing it.

**Quantitative prediction:**

<!-- TODO: human reviews and fills in -->

*Required content (substantive — Crocco §5):* the leading-order proper-time prediction is the same as §3.1's. Record the differences attributable to the Poder analysis pipeline; if the published Poder uncertainty band is tighter than Cole's, the proper-time prediction may pass §3.1 and fail §3.2 (or vice versa). The Mathematica MCP work lives in `Roadmapping/Mathematica_Notebooks/Electromagnetism/Poder2018_proper_time_prediction.wl` (to be created).

**Comparison and verdict:**

<!-- TODO: human reviews and fills in -->

*Required content (substantive — Crocco §5):* same template as §3.1. Special attention to whether the proper-time verdict matches the paper's preferred (quantum-corrected) Landau–Lifshitz comparator or the alternative (classical) Landau–Lifshitz comparator. If the proper-time framework can produce a single comparator that is consistent with both Cole/Poder and Wistisen (§3.3) within published uncertainties, that is a substantive result and feeds the cross-experiment narrative in §4.

### 3.3 Wistisen et al. 2018 — *Nature Commun.* **9**, 795

[[wistisen2018_radiation_reaction]] · PDF (local) · [converted markdown](../../../Historical_Converted_Markdown/Retrospective/wistisen2018_radiation_reaction/wistisen2018_radiation_reaction.md)

**Geometry (pragmatic; restated from the primary source):**

The Wistisen experiment was performed at the CERN SPS North Area H4 beam line using 178.2 GeV **positrons** (not electrons; this corrects the issue text's "electron" phrasing) on aligned silicon crystals of two thicknesses, 3.8 mm and 10.0 mm, oriented along the ⟨111⟩ axis. The crystal's purely electric background field — produced by the lattice of atoms — defines a quantum parameter `χ = γ E / E_cr` with `E_cr = m²c³/(ℏ|e|) ≈ 1.3 × 10¹⁸ V/m`; under the experimental conditions, `χ ≲ 1.4`, the highest of the three 2018 experiments. The channeling parameter `ξ = p_⊥,max/(mc)` ranges from 0.7 to 7, so the synchrotron-limit (constant-crossed-field) approximation `ξ ≫ 1` is satisfied for most but not all of the positron trajectory. Each positron undergoes multiple high-energy photon emissions with non-negligible recoil, so radiation reaction is in the strongly quantum regime. The observable is the *single-photon* emission spectrum, made possible by sending the emitted γ-beam through a thin (200 μm, ≈ 5% radiation length) tantalum converter foil and reconstructing the e⁺e⁻ pair momenta in a downstream magnetic spectrometer — a configuration that suppresses the pile-up which has obscured previous calorimeter-based channeling experiments. Four theoretical models are compared against the measured spectra: classical-plus-radiation-reaction (CRRM, using the LL equation), semiclassical-plus-radiation-reaction (SCRRM, LL with `g(χ)` correction), quantum-plus-radiation-reaction (QRRM, multi-photon emission in the constant-crossed-field approximation), and quantum-without-radiation-reaction (QnoRRM). The reduced chi-squared comparison reported in the paper's Supplementary Note 3 finds that in the 3.8 mm crystal the CRRM fits best (followed by QRRM); in the 10.0 mm crystal the QRRM fits best (CRRM worst); and combining both crystals the QRRM is the overall preferred model — with the explicit caveat that "while the widely used quantum approach is overall the best model, it does not completely describe all the data in this regime."

**Proper-time prediction setup:**

The Wistisen geometry is qualitatively different from §3.1 and §3.2. The positron undergoes periodic transverse oscillation in the crystal's planar or axial potential well; the longitudinal velocity is approximately constant over many oscillation periods. At leading order the acceleration is purely transverse and `u · a ≈ 0`, so the proper-time third term decouples and the prediction reduces to the standard Liénard–Wiechert form. However, longitudinal energy loss from cumulative photon emission gradually slows the positron over the channeling length, producing a non-zero second-order `u_∥ a_∥` integrated over the channel. The amount of cumulative slow-down differs between the 3.8 mm and 10.0 mm crystal thicknesses — substantially more multi-photon recoil accumulates in the 10.0 mm case — so we expect the proper-time third-term contribution to be measurably larger in the thick-crystal data than in the thin-crystal data. The paper's headline observation that QRRM is the overall best model but does not "completely describe all the data" is the specific point against which the proper-time framework can be tested: does the addition of the third term to the QRRM bring the residual discrepancy down? If yes, especially in the 10.0 mm thick-crystal data where the proper-time framework's longitudinal contribution is largest, that is a substantive convergent piece of evidence; if no, the framework's prediction in this regime is constrained by the residual.

**Quantitative prediction:**

<!-- TODO: human reviews and fills in -->

*Required content (substantive — Crocco §5):* compute the second-order longitudinal `u · a` contribution from the cumulative photon-emission-induced energy loss across the published channeling length, and compare to the classical-L-L spectrum that the paper preferred. The Mathematica MCP work lives in `Roadmapping/Mathematica_Notebooks/Electromagnetism/Wistisen2018_proper_time_prediction.wl` (to be created).

**Comparison and verdict:**

<!-- TODO: human reviews and fills in -->

*Required content (substantive — Crocco §5):* same template as §3.1. If the proper-time prediction reduces to classical-L-L in this geometry within the published uncertainty band — which the §3.3 setup paragraph anticipates — that is a *consistency* result, distinct from a *headline* result of the kind §3.1/§3.2 might furnish. The Wistisen case is the most important test of the framework's negative prediction: a framework that predicts deviation in every regime would not be physically credible, and the geometric argument that protects the third term in the channeling case is the same kinematic structure that protects the standard Larmor formula in circular motion.

---

## 4. Cross-experiment summary

<!-- TODO: human reviews and fills in -->

*Required content (substantive — Crocco §5):* fill in the table once the §3 verdicts are in. The published-paper verdict columns are pragmatic restatements; the proper-time verdict column is substantive.

The published-paper columns below restate the verdicts each primary source itself records (pragmatic); the proper-time column is substantive.

| Experiment | Classical L-L (paper) | Semiclassical / QED-corrected (paper) | Proper-time (this document) |
|---|---|---|---|
| Cole 2018 (550 MeV beam, η ≈ 0.07) | ⚠ slightly disfavoured at 1σ (`R = 0.6 ⁺⁰·³₋₀.₂`) | ✅ preferred at 1σ (`R = 0.8 ⁺⁰·⁷₋₀.₃`) | ? |
| Poder 2018 (multi-GeV beam, χ ≈ 0.2) | ⚠ slightly overestimates loss | ✅ semiclassical (LL × `g(χ)`) preferred; QED also acceptable but shows "knee" mismatch | ? |
| Wistisen 2018 (178 GeV positrons, χ ≲ 1.4, two crystal thicknesses) | ✅ best in 3.8 mm crystal; worst in 10.0 mm crystal | QRRM overall best across both thicknesses, "does not completely describe all the data" | ? |

The three published papers split on which comparator best fits each data set, and the split itself partly tracks the regime of `χ`: Cole and Poder, with `η ≲ χ ≈ 0.2`, are best fit by quantum-corrected models; Wistisen, with `χ` up to 1.4, is best fit by the fully quantum QRRM in the thicker crystal but by classical CRRM in the thinner crystal, with neither model fully reproducing the data. It is not obvious that the proper-time framework can produce a single comparator consistent with all three observed spectra; if it can, that is the substantive cross-experiment result. If it cannot — if the proper-time prediction fits the Gemini geometry well but misses the channeling geometry, or vice versa — the locus of the miss is itself a piece of evidence about the framework's domain of applicability, and the corresponding ⚠ or ❌ verdict feeds [§7](#7-findings-cross-post-conditional-deferred).

---

## 5. Crocco compliance

This document mixes pragmatic and substantive AI use per [`CROCCO_COMPLIANCE.md`](../../../Tooling/CROCCO_COMPLIANCE.md) and [#42's plan §6](../../../../.dev/tasks/42-electromagnetism-jackson-proper-time.md#6-source-material-treatment-copyright--crocco-compliance):

- **Pragmatic** (no separate disclosure block in the body): the §2.1 form-level restatement of the third term from the verified Eq. (7); the §2.2 identification of when `u·a` vanishes; the §2.3 scaling form of the radiated-power contribution; the §3 geometric setup paragraphs (restated from the converted primary sources); the operational §1 framing paragraph; the §4 table headers and published-paper columns; this §5; the §6 PR D reconciliation list; the §7 FINDINGS pointer.
- **Substantive** (gated by `<!-- TODO: human reviews and fills in -->` blocks throughout the document body): §1 selection provenance; the §2.3 physical-interpretation paragraph on the third term as a candidate dissolution mechanism for the Abraham–Lorentz pathology; every §3 per-experiment *quantitative prediction* and *verdict*; the §4 cross-experiment summary table's proper-time-verdict column and the cross-experiment narrative paragraph.

The substantive blocks are deliberately blank in this draft. They are gated by Trey's review per the binary `human_reviewed` discipline of Crocco rule 7 and the per-substantive-block accept/reject discipline of Crocco rule 4. If a future AI pass is used to fill them in (rather than direct work by Trey), the pass must use a committed prompt under [`Roadmapping/Tooling/_prompts/`](../../../Tooling/_prompts/) — candidate names `derive_proper_time_prediction.md` (new) or, if the work fits the existing template cluster, `chapter_qa_review.md` (existing, would need adaptation). The prompt-of-record must be committed before the substantive content lands.

---

## 6. Outstanding work blocked on PR D of [#42](https://github.com/temoTxt/PyPhysics/issues/42)

The Ch. 14 chapter document `Ch14_Radiation_by_Moving_Charges.md` scheduled for PR D will provide the canonical per-problem derivation of the proper-time Liénard–Wiechert third term, working through the chain rule from `∂_τ(1/b) = −(u·a)/b³` to the third term's `(u·a)/b⁴` prefactor in problem form. This document re-quotes Eq. (7) self-contained in §2.1 rather than wikilinking that derivation. When PR D lands:

1. The §2.1 form-level restatement is shortened to a leading `**See:** [[Ch14_Radiation_by_Moving_Charges#Liénard–Wiechert proper-time derivation]]` pointer, with the field-expressions table retained for self-contained reference.
2. The §3 per-experiment quantitative-prediction blocks are filled in by reference to PR D's worked example (rather than from scratch in this document).
3. The §5 Crocco compliance is updated to reflect that the §2.1 content is now a *cite* of PR D's canonical derivation rather than an independent re-statement, narrowing the AI use class for §2.1 from "extracted summary of a verified equation" to "verbatim cite".
4. The bib stubs [[cole2018_radiation_reaction]], [[poder2018_radiation_reaction]], [[wistisen2018_radiation_reaction]] are cross-linked from PR D's Ch. 14 prose, since they will then be cited by Ch. 14 in addition to this document; the `chapters_citing` field in each stub is updated accordingly.

PR 0 and PR A of #42 have both shipped (`Electromagnetism/` scaffold, [`_template_problem.md`](../_template_problem.md), [`_proper_time_cheatsheet.md`](../_proper_time_cheatsheet.md), [`Jackson/README.md`](../README.md), Ch01/Ch02 problem docs from PR 0, [Ch05](../Ch05_Magnetostatics_Faraday_Quasi_Static.md) and [Ch06](../Ch06_Maxwell_Equations_Macroscopic_Media.md) problem docs from PR 0 + PR A, and companion Wolfram Language notebooks under [`Mathematica_Notebooks/Electromagnetism/`](../../../Mathematica_Notebooks/Electromagnetism/)). The campaign has therefore validated the per-problem template, the voice constraint, and the Mathematica MCP workflow on nine problems spanning four chapters; the Ch. 6 PR A retrospective ([§7.5 of the #42 plan](../../../../.dev/tasks/42-electromagnetism-jackson-proper-time.md#75-pr-a-retrospective-closed-2026-05-24)) explicitly identifies this document's experimental comparison as the next non-null engagement with the dissipative term, deferred to PR D's Ch. 14 canonical derivation. Only the Ch. 14 problem document is the pending dependency.

---

## 7. FINDINGS cross-post (conditional, deferred)

If any §3 verdict comes back ⚠ or ❌ when Trey fills the substantive blocks, the corresponding finding is recorded in [`Roadmapping/Equation_Verification/FINDINGS_for_author_review.md`](../../../Equation_Verification/FINDINGS_for_author_review.md) following the existing pattern (currently three findings: Maxwell Eq. 24 missing `c` and `V²/(2mc²)` term; DRQM I Eq. III.22 numerical inconsistency in the `g-2` derivation; TCEP Eq. 4.16 sign typo). The new finding(s) are structurally different from the existing three — they are prediction-vs-experiment divergences rather than internal-derivation errors — and the FINDINGS document may need a small format extension to accommodate them cleanly. The extension is best worked out as part of PR D's discussion. No cross-post happens in this PR; the conditional cross-post is a small follow-up PR after the §3 verdicts are recorded.
