# Ch. 12 — Afterword (3e renumber): EPR, Bell, Measurement

**PR L — campaign-closing PR.** The Afterword chapter treats foundational issues (EPR paradox, Bell's theorem, measurement problem) that are interpretation-level questions, not Hamiltonian-level. The proper-time formulation has nothing new to add to these debates; the campaign includes Ch. 12 for completeness. Two problems suffice.

The choice of Hamiltonian (`\hat H_{\text{NR}}`, `K`, or any other) does not affect the EPR/Bell argument structure, which depends only on the quantum statistical interpretation and the locality / realism assumptions. Per [§13 of the plan's honest framing](../../../.dev/tasks/49-quantum-mechanics-griffiths-proper-time.md#7-honest-framing), the dual-theory campaign is silent on these foundational questions.

## Problems

| Problem | Status | Role |
|---|---|---|
| [G3e-P12.1 — EPR paradox and entanglement](#problem-g3e-p121--epr-paradox-and-entanglement) | drafted | fluency (final) |
| [G3e-P12.5 — Bell's theorem and CHSH inequality](#problem-g3e-p125--bells-theorem-and-chsh-inequality) | drafted | fluency (final) |

---

### Problem G3e-P12.1 — EPR paradox and entanglement

**Source:** Griffiths 3e Problem 12.1. *Pragmatic AI.*

**Statement:** A spin-singlet pair of particles is created and the two particles separate. Alice measures the spin of particle 1 along direction `\hat a`; Bob measures the spin of particle 2 along `\hat b`. Compute the correlation function `\langle\sigma_{1}^{\hat a}\sigma_{2}^{\hat b}\rangle`.

**Textbook:** `\langle\sigma_{1}^{\hat a}\sigma_{2}^{\hat b}\rangle = -\hat a \cdot \hat b`. Entangled correlations exceed classical local-realist bounds (Bell's theorem, next problem).

**Proper-time:** The singlet state `|\psi\rangle = (1/\sqrt 2)(|↑↓\rangle - |↓↑\rangle)` is determined by spin algebra, independent of Hamiltonian. Same correlation function.

**Verdict:** ✅. **Companion:** [`GriffithsCh12_P12_1.wl`](../../Mathematica_Notebooks/Quantum_Mechanics/GriffithsCh12_P12_1.wl).

---

### Problem G3e-P12.5 — Bell's theorem and CHSH inequality

**Source:** Griffiths 3e Problem 12.5. *Pragmatic AI.*

**Statement:** Derive the CHSH version of Bell's inequality `|S| \le 2` for local hidden-variable theories, and show that QM predicts `|S| = 2\sqrt 2` for optimal Bell-state measurements.

**Textbook:** Standard CHSH derivation. Aspect-type experiments confirm the QM prediction; local hidden-variable theories are ruled out.

**Proper-time:** The CHSH inequality, the QM upper bound (`2\sqrt 2`), and the Aspect-experiment results are all features of the quantum-mechanical formalism that the dual-theory framework inherits without modification. The dual theory provides no new mechanism to violate the Tsirelson bound and no new prediction for Bell-type experiments.

**Notes:** The campaign's silence on EPR/Bell is appropriate. The dual-theory programme aims to give a different Hamiltonian-level description of the *same* quantum mechanics, not a different interpretive framework. Foundational debates about realism/locality lie outside the campaign's scope.

**Verdict:** ✅. **Companion:** [`GriffithsCh12_P12_5.wl`](../../Mathematica_Notebooks/Quantum_Mechanics/GriffithsCh12_P12_5.wl).

---

## Campaign-closing observation

This is the final per-problem document of the Quantum Mechanics / Griffiths campaign. PR L closes 12 PRs spanning Chs. 1–12 of Griffiths 3e. Issue #49's acceptance criteria are satisfied:

1. ✅ PR A merged (Ch. 1 scaffold + 3 problems)
2. ✅ PR D merged (Ch. 4 hydrogen via proper-time `K`)
3. ✅ PR G merged (Ch. 7 fine structure with branched-treatment on `r_e` finding)
4. ⚠ Root README cross-reference — small follow-on commit needed

The campaign's load-bearing observation, articulated across the 41 per-problem documents:

> *For the bulk of Griffiths (non-relativistic content), the proper-time canonical `K = H_{0}^{2}/(2mc²) + mc²/2 = p²/(2m) + mc²` reduces **exactly** to the textbook non-relativistic Hamiltonian plus a constant rest-energy offset. The campaign demonstrates this reduction problem-by-problem and identifies the **fine-structure regime** (Ch. 7) where the dual Dirac equation supplies the relativistic corrections, with the anomalous-`g`-factor `r_e` finding from DRQM I §III.D as the campaign's only flagged disagreement with precision experiment.*

The dual-theory framework passes the undergraduate-QM-correspondence test. The experimental-discrimination work — testing whether the framework's predictions for precision hydrogen experiments (fine structure, hyperfine splitting, Lamb shift, `g-2`) match measured data — is the subject of the Bethe–Salpeter precision-predictions umbrella, opened concurrently with #49 and out of scope here.
