---
title: "Interim report — verification + campaign work to date"
author: "Trey Morris with Claude Opus 4.7"
date: "2026-05-25"
---

> **Pre-standard notice (per [#115](https://github.com/temoTxt/PyPhysics/issues/115), [#117](https://github.com/temoTxt/PyPhysics/issues/117), and [#120](https://github.com/temoTxt/PyPhysics/issues/120)).** This report was written before the writing standard was codified. It does not yet comply with [`VOICE_MATCH_GILL.md`](../Tooling/VOICE_MATCH_GILL.md) §3.bis (numbered display equations only), §3.ter (short complete sentences; math does the talking), §3.quater (numbered footnotes plus a bibliography section), or §3.quinquies (no emojis or pseudo-emoji symbols). If reissued to a new recipient, it is to be rewritten to the new standard first.

# §1 Cover note

Following our phone call, this report consolidates the verification work and three downstream campaigns since 2026-05-11 into a form that can be read in one sitting. The structure is: how we did the work (§2), what we found in Jackson that we believe will interest you most (§3), the two experimental regimes the Jackson work surfaced as candidates for further attention (§4), the broader unconditional and conditional results across the campaigns (§5), and the single most load-bearing question we have for you (§6 — the resolution of the published $r_e$ value). Secondary errata findings and questions follow in §7. The closing in §8 notes what would help us most to hear back on.

A note on numbering: tags of the form `#NN` throughout this report refer to **GitHub issues or pull requests in the public [`temoTxt/PyPhysics`](https://github.com/temoTxt/PyPhysics) repository**. Each such reference is a clickable link to the corresponding issue or PR; the URLs follow the pattern `https://github.com/temoTxt/PyPhysics/issues/NN` for issues and `…/pull/NN` for PRs. The verification documents, per-PR campaign chapters, Wolfram-Mathematica notebooks, and bibliography stubs cited in this report all live in that repository as well, accessible at the GitHub URLs listed in their respective issues.

## Papers and source texts referenced

The table below lists every paper and textbook this report cites by abbreviation, so the abbreviations below are unambiguous. Dual-theory papers from your corpus appear in the top group; the next group covers other primary sources (Bethe 1947 supplies the Lamb-shift route used in §5); the bottom group covers the campaign textbooks (Jackson, Griffiths, Bethe-Salpeter) whose canonical problems and precision-spectroscopy results the campaigns work through.

| Abbreviation in report | Full title | Authors | Year | Venue |
|---|---|---|---|---|
| **Maxwell paper** | *Two Mathematically Equivalent Versions of Maxwell's Equations* | Gill, Zachary | 2011 | *Foundations of Physics* **41**, 99–128 |
| **Analytic Dirac** | *Analytic Representation of the Dirac Equation* | Gill, Zachary | ca. 2002 | *Foundations of Physics Letters* (subject to author confirmation of citation details) |
| **TCEP** | *The Classical Electron Problem* | Gill, Zachary, Lindesay | 2001 | *Foundations of Physics* **31**, 1299–1354 |
| **DRQM I** | *Dual Relativistic Quantum Mechanics I* | Gill, Ares de Parga, Morris, Wade | 2021 | Manuscript dated 21 August 2021 (Howard University) |
| **Bethe (1947)** | *The Electromagnetic Shift of Energy Levels* | Bethe | 1947 | *Physical Review* **72**, 339 |
| **Jackson** | *Classical Electrodynamics* (3rd ed.) | Jackson | 1998 | Wiley |
| **Griffiths 3e** | *Introduction to Quantum Mechanics* (3rd ed.) | Griffiths, Schroeter | 2018 | Cambridge University Press |
| **Bethe-Salpeter** | *Quantum Mechanics of One- and Two-Electron Atoms* | Bethe, Salpeter | 1957 / 1977 reprint | Plenum 1957; Springer 1977 reprint |

The "Analytic Dirac" venue is left tentative because the converted-PDF metadata does not record the publication line cleanly; one of the Q-items below is whether the citation we should use matches the conference / journal record you intend.

# §2 How we got here — the process

The work began in May 2026 with a systematic Wolfram-Mathematica re-derivation of every numbered equation in the dual-theory paper corpus (Maxwell paper, Analytic Dirac, DRQM I, TCEP — see the table above for full citation details — plus the Banach-space and Navier-Stokes mathematics papers at the load-bearing-results level). The output was [`Roadmapping/Equation_Verification/`](../Equation_Verification/), a per-paper verification document recording each equation as it appears in the paper, a single-line Wolfram check, a step-by-step expanded derivation, a cross-reference to the standard textbook equivalent (Jackson / Sakurai / Goldstein / Weinberg as applicable), and a verdict (✅ / ⚠ / 🔴).


That campaign produced three substantive findings flagged for author review (§7 below) and three load-bearing structural results that survived the re-derivation cleanly. The verification work then opened three downstream campaigns: *[#42 Electromagnetism / Jackson](https://github.com/temoTxt/PyPhysics/issues/42)* (the canonical Jackson problems × CGS/SI × proper-time Maxwell, 66 problems across Chs. 1–16), *[#49 Quantum Mechanics / Griffiths](https://github.com/temoTxt/PyPhysics/issues/49)* (canonical Griffiths problems × 2e/3e × proper-time / dual-theory QM, 41 problems across Chs. 1–12), and *[#50 Bethe-Salpeter precision predictions](https://github.com/temoTxt/PyPhysics/issues/50)* (the precision-experiment counterpart to Griffiths, 28 results across Chs. on hydrogen / helium / Lamb shift / hyperfine / positronium / muonium). The three campaigns are recorded in Pull Requests [PR #51](https://github.com/temoTxt/PyPhysics/pull/51), [PR #52](https://github.com/temoTxt/PyPhysics/pull/52), [PR #53](https://github.com/temoTxt/PyPhysics/pull/53) respectively.


Two discipline choices we want to flag explicitly. First, all AI-assisted work in the repository carries Crocco compliance markers per the *Responsible AI in Non-Empirical Research* framework (Crocco, Rasdi, Garavan 2026): substantive AI work (the kind that shapes what a manuscript argues) requires per-paragraph TODO blocks and a documented human-acceptance pass; pragmatic AI work (transcription, Mathematica verification, formatting) does not. This report is substantive AI end-to-end. Second, every per-PR document carries an explicit honest-framing discipline (the §13 review introduced in the [#42](https://github.com/temoTxt/PyPhysics/issues/42) plan) that distinguishes algebraic identities from physical claims, conditional predictions from unconditional ones, and within-scope verdicts from out-of-scope ones. The Phase 3.5 devil's-advocate pass on this report is the same discipline applied to the report itself.


# §3 Jackson highlights — initial results

The [#42](https://github.com/temoTxt/PyPhysics/issues/42) Electromagnetism campaign worked through 66 canonical Jackson problems under the proper-time / dual-theory formulation of the Maxwell paper. Three results merit foregrounding.


**J3e-P14.2 — the proper-time Liénard-Wiechert third term.** The Maxwell paper's Eq. (7) third term $e(u\cdot a)[r\times(u\times r)]/(b^4 s^3)$ enters the radiation-zone fields with a coefficient that is $(u\cdot a)/b^4$-suppressed at non-relativistic intensity but contributes at order unity for $u/c \approx 1$. The Jackson Ch. 14 problem worked out the explicit angular distribution of the radiation under this term and identified that the prediction is a **longitudinal radiation component** in the far field, absent from classical Liénard-Wiechert (which is purely transverse in the radiation zone). The result is a structural prediction of the formulation, not an empirical claim — but it is the most natural place to look for an operational signature of the dual-theory radiation theory. The Jackson per-PR document under [PR #51](https://github.com/temoTxt/PyPhysics/pull/51) records the algebra and the cross-reference to the Maxwell paper.


**J3e-P16.1 — Abraham-Lorentz dissolution.** The proper-time formulation supplies the radiation-reaction force as a *first-order* $\partial_\tau$ correction (dissipative coefficient $(u\cdot a)/b^4$ from Maxwell paper Eq. 4) rather than the textbook *third-order* $\partial_t^3$ Abraham-Lorentz form. The Jackson Ch. 16 problem worked through the dissolution claim explicitly: the first-order structure removes the runaway-solution pathology and the pre-acceleration pathology that mar the classical theory, because the equation of motion is first-order in proper-time and so does not admit the higher-order spurious solutions. The claim is structural; whether it *also* resolves the empirical defects that motivated the Cole / Poder / Wistisen 2018 experiments is what we believe deserves a separate look — see §4.


**J3e-P5.13 — a non-zero result we did not anticipate.** Working through Jackson's magnetic-moment-of-a-current-loop problem in the proper-time formulation gave a non-zero result where the classical calculation gives zero. The result is flagged with the @Zorns-Lemmon-tagged comment on [issue #42](https://github.com/temoTxt/PyPhysics/issues/42) for cross-verification; we have not yet had a second pair of eyes on the derivation. The campaign's per-PR document records the calculation in full so the result is reproducible.


# §4 Candidate experiments — what we found to look at

The Jackson work surfaced two experimental regimes where the dual-theory formulation's predictions are operationally distinguishable from textbook radiation theory at current measurement precision. We have scoped both as sibling issues to [#42](https://github.com/temoTxt/PyPhysics/issues/42) but have not yet executed the comparison; we flag them here for your awareness.

**Cole 2018 / Poder 2018 / Wistisen 2018 ([#43](https://github.com/temoTxt/PyPhysics/issues/43)).** Three radiation-reaction experiments in the GeV / extreme-intensity regime (Gemini laser-wakefield + counter-propagating laser pulse at $10^{21}$ W/cm² for Cole and Poder; CERN SPS 200 GeV electron channeling through aligned silicon for Wistisen). The three papers themselves disagree on which classical limit fits their data (Cole and Poder favour quantum-corrected Landau-Lifshitz; Wistisen prefers classical Landau-Lifshitz). The disagreement is itself a signal that the regime is ambiguous; the proper-time third term and dissipative coefficient make distinguishable predictions in both geometries. [Issue #43](https://github.com/temoTxt/PyPhysics/issues/43) records the scoping.

**MeV bremsstrahlung at clinical-linac energies ([#48](https://github.com/temoTxt/PyPhysics/issues/48)).** A precision-experiment regime distinct from [#43](https://github.com/temoTxt/PyPhysics/issues/43). Medical-physics calibration measurements at clinical-linac electron energies (1, 6, 18 MeV typical) deliver $\sim 1$–$2\%$ precision on the bremsstrahlung spectrum across decades of data (NIST EGSnrc + AAPM Task Group standards). At MeV electron energies the third term contributes at order unity ($u/c \approx 0.94$), and the prediction is a longitudinal polarisation component absent from classical Born-approximation QED. [Issue #48](https://github.com/temoTxt/PyPhysics/issues/48) records the scoping; the bibliography stubs for the target references (Seltzer-Berger tables, Pratt-Tseng review, NIST EGSnrc) are not yet in place.


# §5 What the campaigns found — unconditional and conditional

Three structural results across the four campaigns are unconditional in the sense that they do not depend on the resolution of the $r_e$ branch question (§6). First, the algebraic structure of the Maxwell paper (Eqs. 22–24 modulo the typos in §7 Finding 1) reproduces under Wolfram verification. Second, the DRQM I §II Foldy-Wouthuysen reduction of the dual-Dirac equation produces the standard Pauli Hamiltonian at leading-$(Z\alpha)^2$ order, with the standard spin-orbit, relativistic-kinetic, and Darwin coefficients. Third, Bethe's 1947 mass-renormalisation estimate of the 2S₁/₂–2P₁/₂ Lamb shift reproduces under the proper-time matrix elements; the campaign's predicted Lamb shift is $\sim 1{,}016$ MHz, the same number the standard Bethe-estimate route gives, with the $\sim 42$ MHz residual to the measured $1{,}057.845(9)$ MHz attributed to dropped sub-leading one-loop QED contributions (out of scope for the campaign's apparatus).


The remaining results across [#50](https://github.com/temoTxt/PyPhysics/issues/50) (Bethe-Salpeter precision predictions) are *conditional* on the resolution of the $r_e$ finding. Six independent precision atomic-physics observables exhibit a consistent branched verdict:

| Observable | Branch (b) prediction | Branch (c) prediction | Measurement | Uncertainty |
|---|---|---|---|---|
| Electron $g_s$ | $-2.0005714$ | $-2.00231930$ | $-2.00231930\ldots$ | $\sim 10^{-12}$ |
| H $2P_{3/2} - 2P_{1/2}$ | $\sim 10{,}952$ MHz | $\sim 10{,}962$ MHz | $10{,}969.13(10)$ MHz | $\sim 10^{-8}$ |
| H 1S hyperfine (21-cm) | $\sim 1{,}418.81$ MHz | $\sim 1{,}420.04$ MHz | $1{,}420.405\,751\,768(2)$ MHz | $\sim 10^{-12}$ |
| He ${}^3P_0 - {}^3P_1$ | $\sim 29{,}617.4$ MHz | $\sim 29{,}616.95$ MHz | $29{,}616.952$ MHz | $\sim 10^{-9}$ |
| Positronium ortho-para | $\sim 203{,}505$ MHz | $\sim 203{,}389$ MHz | $203{,}389(2)$ MHz | $\sim 10^{-5}$ |
| Muonium hyperfine | $\sim 4{,}464.6$ MHz | $\sim 4{,}463.4$ MHz | $4{,}463.302\,776(51)$ MHz | $\sim 10^{-8}$ |

Across all six, branch (b) is in $\sim 10^{-3}$ fractional disagreement with measurement; branch (c) agrees at the campaign's Bethe-estimate precision floor. The pattern is the consequence of the $r_e$ finding propagating through every $g_s$-linear and $g_s^2$-linear observable in the same way.


The Bethe-estimate Lamb shift above is the campaign's **clearest measurement-level non-falsification** at the precision the route can deliver: the framework reproduces the same $\sim 1{,}016$ MHz number as the standard QED Bethe-estimate route gives, and the $\sim 42$ MHz residual is shared between the two routes (it is what neither captures without a full one-loop calculation). It is *not* an "endorsement" in the sense of distinguishing the framework from standard QED — the two formulations give the same Bethe-estimate number — and the campaign's apparatus cannot push below this precision floor (a proper-time one-loop dual-QED calculation is documented as future work in [#55](https://github.com/temoTxt/PyPhysics/issues/55)). The relevant point is that *the framework does not fail* at the Bethe-estimate precision, which is the strongest empirical claim the present apparatus supports. The six observables in the table, by contrast, do depend on $r_e$ at leading order.


# §6 The $r_e$ branch question — the load-bearing ask

This is the single most consequential question we have for you. We have carried two values of $r_e/r_0$ through the campaigns because the published value and the value required to reproduce the measured $g_e$ disagree at the fourth decimal place.


## §6.1 What the two branches are

DRQM I §III.D records the anomalous-$g$-factor formula $g_r = 2(1 - 4 r_0/(2r + r_0))$ together with a stated $r_e/r_0 = 0.499857150068631$ that the paper claims yields $g_e = -2.00231930436256$ (the experimental value). Wolfram evaluation of the formula at the stated $r_e/r_0$ instead gives $g_e = -2.0005714813615487$. The value of $r_e/r_0$ that *does* reproduce $g_e = -2.00231930436256$ is $\approx 0.4994205099128318$. The two values are very close — they differ in the fourth decimal place — but the formula $g_r = 2(1 - 4 r_0/(2r + r_0))$ has a sensitivity $dg_r / d(r_e/r_0) \approx 4$ at the cutoff, so even small differences in $r_e$ correspond to substantial differences in $g_e$. We label the as-published value **branch (b)** and the corrected value that reproduces measured $g_e$ **branch (c)**.


## §6.2 Why the question is load-bearing

The $r_e$ value enters every observable that depends on the electron's anomalous magnetic moment, which means every spin-dependent precision observable in atomic and molecular spectroscopy. The six-observable table in §5 above shows the pattern: branch (b) consistently predicts a $\sim 10^{-3}$ disagreement; branch (c) consistently agrees with measurement at the Bethe-estimate precision floor. On the hydrogen 1S hyperfine line (the 21-cm transition, the most precisely measured atomic-physics frequency at $\sim 10^{-12}$ relative uncertainty), branch (b)'s disagreement is $\sim 6$ orders of magnitude beyond the measurement uncertainty. We are *not* asking you to pick the value that best fits the data — the campaign's verdict-recording is not a fitting exercise. We are asking which value the framework intends, so the campaign can record verdicts ✅ or 🔴 unconditionally rather than carrying both branches indefinitely.


## §6.3 What we are asking, concretely

The question decomposes into three sub-questions, each of which has a clear consequence for the campaign.

- **Q1a.** Is $r_e/r_0 = 0.499857150068631$ (the as-published value, branch (b)) the value the framework intends? If yes, the campaign's verdicts on the six observables in §5 would be recorded as 🔴 (ruled out at current measurement precision), and the open question becomes how the framework reconciles the $\sim 10^{-3}$ disagreement across hyperfine, fine structure, and $g$ — either through a feature of the derivation we have not captured, or through a regime distinction we have missed.

- **Q1b.** Is $r_e/r_0 \approx 0.499420510$ (the corrected value that reproduces measured $g_e$, branch (c)) the value the framework intends? If yes, the published number in DRQM I §III.D is a transcription error to be erratum'd, and the campaign's verdicts collapse to ✅ at the Bethe-estimate precision floor.

- **Q1c.** Is there a third possibility we have not considered — a different formula, a different derivation, or a context-dependent choice where one value is correct in one regime and another in another? If yes, we will rework the campaign's verdicts with the third value in mind.


## §6.4 How the framework supports answering this

Two routes are available. The first is to retrieve the original working notebook in which $r_e$ was computed and confirm the digit string. The campaign does not have access to working notebooks; you do. Even a one-line confirmation that the published value is the intended one, or that the corrected value is the intended one, is sufficient for the campaign to move forward. The second route is to rederive $r_e$ from the dual-Dirac equation's renormalisation prescription — but this requires the prescription itself to be documented (the operator-ordering choice in DRQM I §II.D, which is also Q5 below). The campaign could in principle do the rederivation under the Weyl ordering we have assumed, but the result is only as good as the assumption about the prescription.


# §7 Secondary findings and questions

Three additional errata findings and supporting questions, all of which are pragmatic-level confirmations rather than load-bearing structural decisions.


**Finding 1 — Maxwell paper Eq. (24): two typos.** The published Eq. (24) is missing a factor of $c$ in the $\Sigma \cdot \mathbf{B}$ term (should be $-e \hbar \Sigma \cdot \mathbf{B} / (2mc)$, the value already in Eq. (22)), and is missing the $+V^2 / (2mc^2)$ cross-term that arises in computing $H^2 \to K = H^2/(2mc^2) + mc^2/2$. The campaign's Wolfram check is reproducible; the algebraic argument is straightforward.

**Finding 3 — TCEP Eq. (4.16): sign typo.** Eq. (4.16) writes $v_g = v_g' - v$, but the algebra from the paper's own (4.14)–(4.15) gives $v_g = v_g' + v$. The paper's commentary surrounding the equation also reads as the $+$ sign; the sign in the equation itself appears to be a typesetting error.

**Q4 — $r_\mu$ and $r_p$ numerical values for DRQM I Eq. III.23.** The formulas $g_\mu^a$ and $g_p^a$ for the muon and proton anomalous moments are stated in DRQM I §III.D alongside the electron formula, but the corresponding cutoff values $r_\mu$ and $r_p$ are not given numerically. Without them, muonium and proton magnetic-moment predictions cannot be computed from the framework. If you have these values, they would close out the analogous discriminator analyses we are presently unable to do for the muon and proton.

**Q5 (optional) — operator-ordering choice in DRQM I §II.D.** Sub-leading $(Z\alpha)^4$ results in the dual-Dirac Foldy-Wouthuysen reduction depend on the operator-ordering choice. We have assumed Weyl ordering throughout the campaign, which reproduces the standard Pauli Hamiltonian at leading order. If the framework intends a different ordering, the sub-leading results would shift and the §6.4 route-(ii) rederivation of $r_e$ would proceed under the different prescription.

**Q6 (citation hygiene) — *Analytic Dirac* citation details.** The Analytic Dirac paper's publication venue is recorded as *Foundations of Physics Letters* ca. 2002 in the table on page 1, but the converted PDF's metadata does not record the publication line cleanly. The cleanest fix would be a one-line confirmation of the journal, volume, and year as you intend the citation to read — so the table can be locked in and any future bibliographic export from the repository carries the correct line.


# §8 Closing

Thank you for the phone call and for the opening to send this. The single most useful response would be a resolution of §6 Q1a/b/c — even a one-line answer determines which of the campaign's six precision verdicts collapse to ✅ and which become 🔴. Confirmations or corrections on Finding 1, Finding 3, and Q4 ($r_\mu$, $r_p$) would close out the verification-level items. The campaigns continue under the current trajectory if you find the work useful in this form; if you have a different direction you would prefer, we will redirect.

