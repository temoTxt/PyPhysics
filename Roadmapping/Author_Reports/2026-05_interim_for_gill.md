---
title: "Interim report — verification and campaign work to date"
author: "Trey Morris with Claude Opus 4.7"
date: "2026-05-25"
---

# §1 Cover note

Following our phone call, we consolidate here the verification work and the three downstream campaigns undertaken since 2026-05-11. The intent is a document that can be read in one sitting and responded to in another.

The structure is as follows. §2 describes the process. §3 records the three Jackson results we believe will interest you most. §4 lists the two experimental regimes the Jackson work surfaced as candidates for further attention. §5 collects the broader unconditional and conditional results across the four campaigns. §6 carries the single most load-bearing question we have for you, namely the resolution of the published value of $r_{e}$. §7 records the secondary errata findings and citation questions. §8 closes.

A note on numbering. The original draft of this report carried inline issue and pull-request tags drawn from the public PyPhysics repository[^13]. Every such tag has been moved into a numbered footnote. The file paths and external links cited by the campaign documents have been moved likewise. The bibliography section at the end resolves each footnote.

## Papers and source texts referenced

The table below lists every paper and textbook cited by abbreviation, so that the abbreviations are unambiguous in what follows. The first group covers the dual-theory papers from your corpus. The second group covers other primary sources used in §5 (Bethe 1947 supplies the Lamb-shift route). The third group covers the campaign textbooks (Jackson, Griffiths, Bethe-Salpeter).

| Abbreviation | Full title | Authors | Year | Venue |
|---|---|---|---|---|
| Maxwell paper | *Two Mathematically Equivalent Versions of Maxwell's Equations* | Gill, Zachary | 2011 | *Foundations of Physics* **41**, 99–128 |
| Analytic Dirac | *Analytic Representation of the Dirac Equation* | Gill, Zachary | ca. 2002 | *Foundations of Physics Letters* (subject to author confirmation) |
| TCEP | *The Classical Electron Problem* | Gill, Zachary, Lindesay | 2001 | *Foundations of Physics* **31**, 1299–1354 |
| DRQM I | *Dual Relativistic Quantum Mechanics I* | Gill, Ares de Parga, Morris, Wade | 2021 | Manuscript dated 21 August 2021 (Howard University) |
| Bethe (1947) | *The Electromagnetic Shift of Energy Levels* | Bethe | 1947 | *Physical Review* **72**, 339 |
| Jackson | *Classical Electrodynamics* (3rd ed.) | Jackson | 1998 | Wiley |
| Griffiths 3e | *Introduction to Quantum Mechanics* (3rd ed.) | Griffiths, Schroeter | 2018 | Cambridge University Press |
| Bethe-Salpeter | *Quantum Mechanics of One- and Two-Electron Atoms* | Bethe, Salpeter | 1957 / 1977 reprint | Plenum 1957; Springer 1977 reprint |

The Analytic Dirac venue is left tentative because the converted-PDF metadata does not record the publication line cleanly. One of the questions in §7 (Q6) is whether the citation we have matches the journal record you intend.

# §2 How we got here — the process

The work began in May 2026 with a systematic Wolfram-Mathematica re-derivation of every numbered equation in the dual-theory paper corpus. The corpus covers the Maxwell paper, the Analytic Dirac paper, DRQM I, and TCEP. The Banach-space and Navier-Stokes mathematics papers were verified at the load-bearing-results level only. The output is a per-paper verification document[^1] recording each equation as it appears in the paper, a single-line Wolfram check, a step-by-step expanded derivation, a cross-reference to the standard textbook equivalent (Jackson, Sakurai, Goldstein, or Weinberg as applicable), and a verdict.

That verification campaign produced three substantive findings flagged for author review (§7) and three load-bearing structural results that survived the re-derivation cleanly. The verification work then opened three downstream campaigns.

The first is the Electromagnetism / Jackson campaign[^2]. It works through 66 canonical Jackson problems in CGS and SI under the proper-time Maxwell formulation. The problems span Chapters 1–16.

The second is the Quantum Mechanics / Griffiths campaign[^3]. It works through 41 canonical Griffiths problems across editions 2 and 3, under the proper-time / dual-theory QM formulation. The problems span Chapters 1–12.

The third is the Bethe-Salpeter precision-prediction campaign[^4]. It is the precision-experiment counterpart to Griffiths. It covers 28 results across the hydrogen, helium, Lamb-shift, hyperfine, positronium, and muonium chapters.

The three campaigns are recorded in three pull requests[^5][^6][^7].

Two discipline choices we wish to flag. First, all AI-assisted work in the repository carries Crocco compliance markers per the *Responsible AI in Non-Empirical Research* framework[^8]. Under that framework, substantive AI work (the kind that shapes what a manuscript argues) requires per-paragraph review blocks and a documented human-acceptance pass; pragmatic AI work (transcription, Mathematica verification, formatting) does not. The present report is substantive AI end-to-end.

Second, every per-PR document carries an explicit honest-framing discipline, introduced as the §13 review in the Jackson plan. The discipline separates algebraic identities from physical claims. It separates conditional predictions from unconditional ones. It separates within-scope verdicts from out-of-scope ones. The Phase 3.5 devil's-advocate pass on the present report is the same discipline applied to the report itself.

# §3 Jackson highlights — initial results

The Electromagnetism campaign[^2] worked through 66 canonical Jackson problems under the proper-time dual-theory formulation of the Maxwell paper. Three results merit foregrounding.

## §3.1 Problem 14.2 — the proper-time Liénard-Wiechert third term

We observe that Eq. (7) of the Maxwell paper contains a third term in the radiation-zone fields,
$$
\frac{e\,(\mathbf{u}\cdot\mathbf{a})\,[\mathbf{r}\times(\mathbf{u}\times\mathbf{r})]}{b^{4}\,s^{3}},
\qquad b = \sqrt{c^{2} + \mathbf{u}^{2}},
\tag{1}
$$
which is absent from the classical Liénard-Wiechert expression. The coefficient is suppressed by $(\mathbf{u}\cdot\mathbf{a})/b^{4}$ at non-relativistic intensity. It contributes at order unity for $\mathbf{u}/c \approx 1$.

Working through Jackson Problem 14.2 under the dual-theory formulation, one obtains the explicit angular distribution of the radiation under (1). It follows that the formulation predicts a longitudinal radiation component in the far field. The classical Liénard-Wiechert expression is purely transverse in the radiation zone.

The result is a structural prediction, not yet an empirical claim. Indeed, it is (to the extent we have looked) the most natural place to seek an operational signature of the dual-theory radiation theory. The Jackson per-PR document under the Electromagnetism PR[^5] records the algebra and the cross-reference to Eq. (7) of the Maxwell paper.

## §3.2 Problem 16.1 — Abraham-Lorentz dissolution

The proper-time formulation supplies the radiation-reaction force as a first-order $\partial_{\tau}$ correction. The dissipative coefficient is the $(\mathbf{u}\cdot\mathbf{a})/b^{4}$ factor of Eq. (4) of the Maxwell paper. By contrast, the textbook Abraham-Lorentz force is third-order in $\partial_{t}$,
$$
\mathbf{F}_{\rm AL} \;=\; \tfrac{2}{3}\,\frac{e^{2}}{c^{3}}\,\dddot{\mathbf{x}}.
\tag{2}
$$
The third-order time derivative in (2) is the source of the runaway-solution and pre-acceleration pathologies of the classical theory. The proper-time and standard equations of motion are mathematically equivalent at the level of particle trajectories, but not (what we mean by) physically equivalent at the level of radiation-reaction structure.

We observe that the first-order proper-time structure removes both pathologies. The equation of motion is first-order in $\tau$ and therefore does not admit the higher-order spurious solutions. The Jackson Problem 16.1 derivation makes the dissolution claim explicit.

The claim is structural. Whether it also resolves the empirical defects that motivated the Cole, Poder, and Wistisen 2018 experiments is what we believe deserves a separate look. The candidate experiments are recorded in §4.

## §3.3 Problem 5.13 — a non-zero result we did not anticipate

Working through Jackson's magnetic-moment-of-a-current-loop problem in the proper-time formulation gave a non-zero result, where the classical calculation gives zero. The result is flagged with a tagged comment on the Electromagnetism issue[^9] for cross-verification. An independent re-derivation has not yet been performed. The campaign's per-PR document records the calculation in full so that the result is reproducible.

# §4 Candidate experiments — what we found to look at

The Jackson work surfaced two experimental regimes in which the dual-theory predictions appear to be operationally distinguishable from textbook radiation theory at current measurement precision. We have scoped both as sibling issues to the Electromagnetism master, but have not yet executed the comparison. We flag them here for your awareness.

## §4.1 Cole 2018 / Poder 2018 / Wistisen 2018

Three radiation-reaction experiments operate in the GeV and extreme-intensity regime. Cole and Poder use the Gemini laser-wakefield accelerator with a counter-propagating laser pulse at intensity
$$
I \;\sim\; 10^{21}\;{\rm W/cm^{2}}.
\tag{3}
$$
Wistisen uses the CERN SPS 200 GeV electron beam channeling through aligned silicon. The three papers themselves disagree on which classical limit fits their data: Cole and Poder favour quantum-corrected Landau-Lifshitz, while Wistisen prefers classical Landau-Lifshitz. The disagreement is itself a signal that the regime is ambiguous. The proper-time third term of (1) and the dissipative coefficient of (2) make distinguishable predictions in both geometries. The scoping is recorded on the dedicated sibling issue[^10].

## §4.2 MeV bremsstrahlung at clinical-linac energies

A precision-experiment regime distinct from §4.1 emerges from medical-physics calibration. Calibration measurements at clinical-linac electron energies (typically 1, 6, and 18 MeV) deliver $\sim 1$–$2\%$ precision on the bremsstrahlung spectrum across decades of data, against the NIST EGSnrc and AAPM Task Group standards. At MeV electron energies the third term of (1) contributes at order unity, with
$$
\frac{u}{c} \;\approx\; 0.94 \quad\text{at}\quad E_{e} = 1\;{\rm MeV}.
\tag{4}
$$
The prediction is a longitudinal polarisation component absent from classical Born-approximation QED. The scoping is recorded on the dedicated sibling issue[^11]. The bibliography stubs for the target references (Seltzer-Berger tables, Pratt-Tseng review, NIST EGSnrc) are not yet in place.

# §5 What the campaigns found — unconditional and conditional

Three structural results across the four campaigns are unconditional, in the sense that they do not depend on the resolution of the $r_{e}$ branch question of §6.

First, the algebraic structure of the Maxwell paper (Eqs. (22)–(24), modulo the typos recorded under Finding 1 of §7) reproduces under Wolfram verification.

Second, the §II Foldy-Wouthuysen reduction of the dual-Dirac equation in DRQM I produces the standard Pauli Hamiltonian at leading $(Z\alpha)^{2}$ order, with the standard spin-orbit, relativistic-kinetic, and Darwin coefficients.

Third, Bethe's 1947 mass-renormalisation estimate of the 2S$_{1/2}$–2P$_{1/2}$ Lamb shift reproduces under the proper-time matrix elements. We obtain the campaign's predicted Lamb shift,
$$
\Delta\nu_{\rm Lamb}^{\rm framework}
\;\approx\; 1{,}016\;{\rm MHz},
\tag{5}
$$
which is the same number that the standard Bethe-estimate route gives. The residual to the measured value,
$$
\Delta\nu_{\rm Lamb}^{\rm meas} \;=\; 1{,}057.845(9)\;{\rm MHz},
\tag{6}
$$
is $\sim 42$ MHz. The residual is attributable to the sub-leading one-loop QED contributions dropped by the Bethe estimate, which are out of scope for the campaign's apparatus.

The remaining results across the Bethe-Salpeter campaign[^4] are *conditional* on the resolution of the $r_{e}$ finding. Six independent precision atomic-physics observables exhibit a consistent branched verdict, set out in the table below.

| Observable | Branch (b) prediction | Branch (c) prediction | Measurement | Uncertainty |
|---|---|---|---|---|
| Electron $g_{s}$ | $-2.000\,571\,4$ | $-2.002\,319\,30$ | $-2.002\,319\,30\ldots$ | $\sim 10^{-12}$ |
| H $2P_{3/2} - 2P_{1/2}$ | $\sim 10{,}952$ MHz | $\sim 10{,}962$ MHz | $10{,}969.13(10)$ MHz | $\sim 10^{-8}$ |
| H 1S hyperfine (21 cm) | $\sim 1{,}418.81$ MHz | $\sim 1{,}420.04$ MHz | $1{,}420.405\,751\,768(2)$ MHz | $\sim 10^{-12}$ |
| He ${}^{3}P_{0} - {}^{3}P_{1}$ | $\sim 29{,}617.4$ MHz | $\sim 29{,}616.95$ MHz | $29{,}616.952$ MHz | $\sim 10^{-9}$ |
| Positronium ortho-para | $\sim 203{,}505$ MHz | $\sim 203{,}389$ MHz | $203{,}389(2)$ MHz | $\sim 10^{-5}$ |
| Muonium hyperfine | $\sim 4{,}464.6$ MHz | $\sim 4{,}463.4$ MHz | $4{,}463.302\,776(51)$ MHz | $\sim 10^{-8}$ |

Across all six observables, branch (b) is in $\sim 10^{-3}$ fractional disagreement with measurement, while branch (c) agrees with measurement at the campaign's Bethe-estimate precision floor. The pattern is the consequence of the $r_{e}$ finding propagating through every $g_{s}$-linear and $g_{s}^{2}$-linear observable in the same way.

The Bethe-estimate Lamb shift of (5)–(6) is the campaign's clearest measurement-level non-falsification at the precision the route can deliver. The framework reproduces the same $\sim 1{,}016$ MHz number that the standard QED Bethe-estimate route gives. The $\sim 42$ MHz residual is shared between the two routes; it is what neither captures without a full one-loop calculation. The result is not, in the strict sense, a discriminator between the framework and standard QED. The two formulations give the same Bethe-estimate number. The campaign's apparatus cannot push below this precision floor. A proper-time one-loop dual-QED calculation is documented as future work[^12]. The relevant point is that the framework does not fail at the Bethe-estimate precision. That is the strongest empirical claim the present apparatus supports. The six observables in the table, by contrast, do depend on $r_{e}$ at leading order.

# §6 The $r_{e}$ branch question — the load-bearing ask

This is the single most consequential question we have for you. We have carried two values of $r_{e}/r_{0}$ through the campaigns, because the published value and the value required to reproduce the measured $g_{e}$ disagree at the fourth decimal place.

## §6.1 What the two branches are

DRQM I §III.D records the anomalous-$g$-factor formula
$$
g_{r} \;=\; 2\!\left(1 - \frac{4\,r_{0}}{2r + r_{0}}\right),
\tag{7}
$$
together with a stated value
$$
\left(\frac{r_{e}}{r_{0}}\right)_{\!\rm published} \;=\; 0.499\,857\,150\,068\,631,
\tag{8}
$$
that the paper claims yields
$$
g_{e}^{\rm claimed} \;=\; -2.002\,319\,304\,362\,56,
\tag{9}
$$
which is the experimental value of $g_{e}$. We have evaluated (7) at the value (8) in Wolfram and instead obtain
$$
g_{e}^{\rm Wolfram} \;=\; -2.000\,571\,481\,361\,548\,7.
\tag{10}
$$
The value of $r_{e}/r_{0}$ that does reproduce (9) is
$$
\left(\frac{r_{e}}{r_{0}}\right)_{\!\rm corrected} \;\approx\; 0.499\,420\,509\,912\,831\,8.
\tag{11}
$$
The two values (8) and (11) are very close — they differ in the fourth decimal place — but the formula (7) has a sensitivity
$$
\left.\frac{d g_{r}}{d(r_{e}/r_{0})}\right|_{r_{e}/r_{0}\approx 1/2} \;\approx\; 4,
\tag{12}
$$
so that even small differences in $r_{e}$ correspond to substantial differences in $g_{e}$. We label the as-published value (8) **branch (b)**, and the corrected value (11) that reproduces the measured $g_{e}$ **branch (c)**.

## §6.2 Why the question is load-bearing

The value of $r_{e}$ enters every observable that depends on the electron's anomalous magnetic moment. That is every spin-dependent precision observable in atomic and molecular spectroscopy. The six-observable table in §5 shows the pattern. Branch (b) consistently predicts a $\sim 10^{-3}$ disagreement; branch (c) consistently agrees with measurement at the Bethe-estimate precision floor.

Consider the hydrogen 1S hyperfine line, the 21-cm transition. It is the most precisely measured atomic-physics frequency, at $\sim 10^{-12}$ relative uncertainty. On that line, branch (b)'s disagreement is of order $\sim 10^{6}$ standard deviations beyond the measurement uncertainty.

We are not asking you to pick the value that best fits the data. The campaign's verdict-recording is not a fitting exercise. We are asking which value the framework intends, so that the campaign can record verdicts as **Pass** or **Refuted** unconditionally rather than carry both branches indefinitely.

## §6.3 What we are asking, concretely

The question decomposes into three sub-questions, each of which has a clear consequence for the campaign.

**Q1a.** Is the as-published value (8) the value the framework intends, i.e. branch (b)? If yes, the campaign's verdicts on the six observables of §5 would be recorded as **Refuted** (ruled out at current measurement precision). The open question would then become how the framework reconciles the $\sim 10^{-3}$ disagreement across hyperfine, fine structure, and $g$ — either through a feature of the derivation we have not captured, or through a regime distinction we have missed.

**Q1b.** Is the corrected value (11) the value the framework intends, i.e. branch (c)? If yes, the published number (8) in DRQM I §III.D is a transcription error to be erratum'd, and the campaign's verdicts collapse to **Pass** at the Bethe-estimate precision floor.

**Q1c.** Is there a third possibility we have not considered — a different formula, a different derivation, or a context-dependent choice where one value is correct in one regime and another in another? If yes, we will rework the campaign's verdicts with the third value in mind.

## §6.4 How the framework supports answering this

Two routes are available. The first is to retrieve the original working notebook in which $r_{e}$ was computed and confirm the digit string. The campaign does not have access to the working notebooks; you do. Even a one-line confirmation that the published or the corrected value is the intended one is sufficient for the campaign to move forward.

The second route is to rederive $r_{e}$ from the renormalisation prescription of the dual-Dirac equation. This requires the prescription itself to be documented (the operator-ordering choice in DRQM I §II.D, which is also Q5 of §7). The campaign could in principle perform the rederivation under the Weyl ordering we have assumed, but the result is only as good as the assumption about the prescription.

# §7 Secondary findings and questions

Three additional errata findings and three supporting questions, all of which are pragmatic-level confirmations rather than load-bearing structural decisions.

**Finding 1 — Maxwell paper Eq. (24): two typos.** The published Eq. (24) is missing a factor of $c$ in the $\Sigma\cdot\mathbf{B}$ term. The term should read
$$
-\,\frac{e\,\hbar\,\Sigma\cdot\mathbf{B}}{2 m c},
\tag{13}
$$
which is the value already in Eq. (22). Furthermore, the published Eq. (24) is missing the cross-term
$$
+\,\frac{V^{2}}{2 m c^{2}},
\tag{14}
$$
which arises in computing $H^{2} \to K = H^{2}/(2 m c^{2}) + m c^{2}/2$. The campaign's Wolfram check is reproducible, and the algebraic argument is straightforward.

**Finding 3 — TCEP Eq. (4.16): sign typo.** Equation (4.16) is written as
$$
v_{g} \;=\; v_{g}' - v,
\tag{15}
$$
but the algebra from the paper's own (4.14)–(4.15) gives
$$
v_{g} \;=\; v_{g}' + v.
\tag{16}
$$
The paper's commentary surrounding the equation also reads as the $+$ sign. The sign in (15) appears to be a typesetting error.

**Q4 — $r_{\mu}$ and $r_{p}$ numerical values for DRQM I Eq. (III.23).** The formulas for the muon and proton anomalous moments,
$$
g_{\mu}^{a}, \qquad g_{p}^{a},
\tag{17}
$$
are stated in DRQM I §III.D alongside the electron formula. The corresponding cutoff values $r_{\mu}$ and $r_{p}$ are not given numerically. Without them, the muonium and proton magnetic-moment predictions cannot be computed from the framework. If you have these values, they would close out the analogous discriminator analyses we are presently unable to do for the muon and proton.

**Q5 (optional) — operator-ordering choice in DRQM I §II.D.** The sub-leading $(Z\alpha)^{4}$ results in the dual-Dirac Foldy-Wouthuysen reduction depend on the operator-ordering choice. We have assumed Weyl ordering throughout the campaign, which reproduces the standard Pauli Hamiltonian at leading order. If the framework intends a different ordering, the sub-leading results would shift, and the §6.4 route-(ii) rederivation of $r_{e}$ would proceed under the different prescription.

**Q6 (citation hygiene) — Analytic Dirac citation details.** The Analytic Dirac paper's publication venue is recorded in the table at the head of this report as *Foundations of Physics Letters*, ca. 2002. The converted-PDF metadata does not record the publication line cleanly. The cleanest fix would be a one-line confirmation of the journal, volume, and year as you intend the citation to read, so that the table can be locked in and any future bibliographic export from the repository carries the correct line.

# §8 Closing

Thank you for the phone call and for the opening to send this. The single most useful response would be a resolution of §6 (Q1a, Q1b, Q1c). Even a one-line answer determines which of the campaign's six precision verdicts collapse to **Pass** and which become **Refuted**. Confirmations or corrections on Finding 1, Finding 3, and Q4 ($r_{\mu}$, $r_{p}$) would close out the verification-level items. The campaigns continue under the current trajectory if you find the work useful in this form. If you have a different direction you would prefer, we will redirect.

## References

[^1]: Per-paper verification documents, `Roadmapping/Equation_Verification/`, repository `temoTxt/PyPhysics`.

[^2]: Electromagnetism / Jackson canonical-problems campaign. Issue #42, `github.com/temoTxt/PyPhysics/issues/42`.

[^3]: Quantum Mechanics / Griffiths canonical-problems campaign. Issue #49, `github.com/temoTxt/PyPhysics/issues/49`.

[^4]: Bethe-Salpeter precision-prediction campaign. Issue #50, `github.com/temoTxt/PyPhysics/issues/50`.

[^5]: Electromagnetism per-PR campaign record. PR #51, `github.com/temoTxt/PyPhysics/pull/51`.

[^6]: Quantum Mechanics per-PR campaign record. PR #52, `github.com/temoTxt/PyPhysics/pull/52`.

[^7]: Bethe-Salpeter per-PR campaign record. PR #53, `github.com/temoTxt/PyPhysics/pull/53`.

[^8]: Crocco, Rasdi, and Garavan, *Responsible AI in Non-Empirical Research* (2026). Repository compliance mapping: `Roadmapping/Tooling/CROCCO_COMPLIANCE.md`.

[^9]: Cross-verification tag (`@Zorns-Lemmon`) on the Electromagnetism master issue. Issue #42, `github.com/temoTxt/PyPhysics/issues/42`.

[^10]: Cole 2018 / Poder 2018 / Wistisen 2018 radiation-reaction scoping. Issue #43, `github.com/temoTxt/PyPhysics/issues/43`.

[^11]: MeV-bremsstrahlung clinical-linac scoping. Issue #48, `github.com/temoTxt/PyPhysics/issues/48`.

[^12]: Proper-time one-loop dual-QED Lamb-shift calculation (future work). Issue #55, `github.com/temoTxt/PyPhysics/issues/55`.

[^13]: Public PyPhysics repository, `github.com/temoTxt/PyPhysics`.
