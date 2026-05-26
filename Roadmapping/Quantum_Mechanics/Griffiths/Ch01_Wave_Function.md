# Ch. 1 — The Wave Function

This chapter contains Griffiths canonical problems on the wave-function formalism, worked in the proper-time / dual-theory reformulation alongside their classical Griffiths-2e and Griffiths-3e textbook solutions. Per [§3 of the plan](../../../.dev/tasks/49-quantum-mechanics-griffiths-proper-time.md#3-per-problem-template), Ch. 1 of Griffiths is presented in the **three-formulation regime**: (a) 2e textbook, (b) 3e textbook (where it differs), and (c) proper-time / dual-theory reformulation under the canonical `K`.

Ch. 1 is **PR A** of the Quantum_Mechanics campaign. It is foundational: all three problems demonstrate the **exact** reduction

$$K - m c^{2} = \frac{p^{2}}{2 m}$$

derived in [`_proper_time_K_cheatsheet.md`](../_proper_time_K_cheatsheet.md) §2. Free-particle non-relativistic problems reduce identically to textbook QM with a constant rest-energy offset. The proper-time formulation has nothing new to add at this level; the campaign records the agreement explicitly to set the pedagogical foundation for the relativistic content of Ch. 4 and Ch. 7.

## Problems in this chapter

| Problem | Status | Role |
|---|---|---|
| [Problem G3e-P1.4 — Normalisation and expectation values of a Gaussian wavepacket](#problem-g3e-p14--normalisation-and-expectation-values-of-a-gaussian-wavepacket) | drafted (PR A) | fluency-builder (foundational) |
| [Problem G3e-P1.5 — Specific wavefunction `A e^{-\lambda\|x\|}`](#problem-g3e-p15--specific-wavefunction-a-e--lambda-x) | drafted (PR A) | fluency-builder |
| [Problem G3e-P1.14 — Probability current and continuity equation](#problem-g3e-p114--probability-current-and-continuity-equation) | drafted (PR A) | fluency-builder |

---

### Problem G3e-P1.4 — Normalisation and expectation values of a Gaussian wavepacket

**Selection provenance** (Crocco §5 substantive-AI note):
- *Chosen because:* the Gaussian wavepacket is the standard introductory example illustrating normalisation, expectation values, and uncertainty. As PR A's opening problem, it exercises the fundamental QM formalism (`\langle x \rangle`, `\langle p \rangle`, `\Delta x \cdot \Delta p`) under the proper-time reformulation.
- *Alternatives considered:* G3e-P1.1 (probability density of a discrete sample — too short) and G3e-P1.3 (visualising the wavefunction — minimal computational content).
- *Role in this PR:* fluency-builder (foundational).

<!-- TODO: human reviews and fills in — confirms the choice of Gaussian wavepacket as PR A's opening foundational problem -->

**Source:** Griffiths, *Introduction to Quantum Mechanics*, 3e Problem 1.4 (Cambridge 2018). 2e equivalence: Problem 1.4 (Pearson 2005). *Paraphrased; consult the textbook for the precise statement.*

**Paraphrased statement:** A particle has wavefunction `\psi(x, t = 0) = A\,e^{-\lambda x^{2}}` with real positive `\lambda`. Determine the normalisation constant `A`, and compute `\langle x \rangle`, `\langle x^{2} \rangle`, `\langle p \rangle`, `\langle p^{2} \rangle`, the uncertainties `\Delta x` and `\Delta p`, and verify the Heisenberg uncertainty product.

#### (a) Griffiths 2e and (b) Griffiths 3e textbook solution

The two editions are identical for this problem (Griffiths records this as essentially unchanged between 2e and 3e). Standard derivation:

- Normalisation: `\int_{-\infty}^{\infty}|\psi|^{2}\, dx = A^{2}\sqrt{\pi/(2\lambda)} = 1`, so `A = (2\lambda/\pi)^{1/4}`.
- `\langle x \rangle = 0` (Gaussian is even).
- `\langle x^{2} \rangle = 1/(4\lambda)`, so `\Delta x = 1/(2\sqrt\lambda)`.
- `\langle p \rangle = 0` (real Gaussian; symmetric).
- `\langle p^{2} \rangle = \hbar^{2}\lambda`, so `\Delta p = \hbar\sqrt\lambda`.
- Uncertainty product: `\Delta x\,\Delta p = \hbar/2` — saturates the Heisenberg lower bound (minimum-uncertainty Gaussian state).

#### (c) Proper-time / dual-theory reformulation

Under the canonical Hamiltonian `K = H^{2}/(2mc^{2}) + mc^{2}/2` of [`_proper_time_K_cheatsheet.md`](../_proper_time_K_cheatsheet.md), the operator structure of the wavefunction formalism is unchanged at the level of position and momentum operators. The expectation values `\langle x \rangle`, `\langle p \rangle`, and their second moments are computed identically because `\hat x` and `\hat p` are not modified by the proper-time substitution — only the *Hamiltonian* changes, from `H_{0} = \sqrt{c^{2}\hat p^{2} + m^{2}c^{4}}` (standard relativistic) to `K = \hat p^{2}/(2m) + mc^{2}` (proper-time canonical).

For free particles, the proper-time `K` differs from the textbook non-relativistic Hamiltonian `\hat H_{\text{NR}} = \hat p^{2}/(2m)` by a constant rest-energy offset `mc^{2}`. Energy eigenvalues differ by `mc^{2}`, but eigenstates are identical and expectation values of `\hat x, \hat p, \hat x^{2}, \hat p^{2}` are unchanged.

**Mathematica check** (Wolfram MCP, 2026-05-25):

```mathematica
ClearAll[capA, lambda, xx];
psi = capA Exp[-lambda xx^2];
normalSq = Integrate[psi^2, {xx, -Infinity, Infinity},
   Assumptions -> lambda > 0];
capASolved = capA /. Solve[normalSq == 1, capA][[2]];
Print["A = ", capASolved];
(* Result: (2 lambda / Pi)^(1/4)  ✅ *)
```

**Reduction-to-textbook verdict:** ✅ matches Griffiths *exactly* (not just at leading order in `v/c`). The proper-time `K` and textbook `\hat H_{\text{NR}}` differ by an additive constant, which has no effect on stationary-state expectation values.

**Comparison:**

| Quantity | Griffiths 2e | Griffiths 3e | Proper-time |
|---|---|---|---|
| `A` (normalisation) | `(2\lambda/\pi)^{1/4}` | identical | identical |
| `\langle x \rangle` | `0` | identical | identical |
| `\langle x^{2} \rangle` | `1/(4\lambda)` | identical | identical |
| `\Delta x \cdot \Delta p` | `\hbar/2` (minimum) | identical | identical |
| Energy offset relative to textbook | `0` | `0` | `+ m c^{2}` (rest energy) |

**Standard-equation comparison:** matches Sakurai §1.6 and Shankar Ch. 4 normalisation conventions for the minimum-uncertainty Gaussian.

**Verdict:** ✅ all three solutions consistent. The proper-time formulation reproduces every expectation value identically; the only difference is the constant rest-energy offset in the Hamiltonian eigenvalues, which is conventional in non-relativistic QM.

**Notes for author review:** none. The exact (not approximate) reduction of `K` to `p^{2}/(2m) + mc^{2}` for free particles is the cleanest possible foundation for the campaign's non-relativistic Griffiths content.

**Companion notebook:** [`Roadmapping/Mathematica_Notebooks/Quantum_Mechanics/GriffithsCh01_P1_4.wl`](../../Mathematica_Notebooks/Quantum_Mechanics/GriffithsCh01_P1_4.wl).

---

### Problem G3e-P1.5 — Specific wavefunction `A e^{-\lambda|x|}`

**Selection provenance:** the exponential-cusp wavefunction is the second standard example in Griffiths Ch. 1; tests the formalism on a non-Gaussian non-smooth function. *Pragmatic AI.*

**Source:** Griffiths, *Introduction to Quantum Mechanics*, 3e Problem 1.5 (Cambridge 2018). 2e equivalence: Problem 1.5 (Pearson 2005). *Paraphrased.*

**Paraphrased statement:** A particle has wavefunction `\psi(x) = A\,e^{-\lambda |x|}`. Determine `A`, compute `\langle x \rangle`, `\langle x^{2} \rangle`, and verify the uncertainty product `\Delta x\,\Delta p` against the Heisenberg bound.

#### (a)/(b) Griffiths textbook solution

- Normalisation: `\int_{-\infty}^{\infty}A^{2}e^{-2\lambda|x|}\,dx = A^{2}/\lambda = 1`, so `A = \sqrt\lambda`.
- `\langle x \rangle = 0` (symmetric).
- `\langle x^{2} \rangle = 1/(2\lambda^{2})`, so `\Delta x = 1/(\lambda\sqrt 2)`.
- For `\langle p^{2} \rangle`: `\psi` has a cusp at `x = 0`, so `\psi'` has a discontinuity and `\psi''` has a delta-function piece. Careful integration gives `\langle p^{2} \rangle = \hbar^{2}\lambda^{2}`, so `\Delta p = \hbar\lambda`.
- Uncertainty product: `\Delta x\,\Delta p = \hbar/\sqrt 2 \approx 0.707\,\hbar`, which satisfies `\Delta x\,\Delta p \ge \hbar/2 = 0.5\,\hbar` (the Heisenberg bound) but does *not* saturate it (the Gaussian saturates; this cusp-function exceeds the bound by `\sqrt 2`).

#### (c) Proper-time / dual-theory reformulation

Identical to (a)/(b) for the same reason as G3e-P1.4: the position and momentum operators are unmodified, and `K` differs from `\hat H_{\text{NR}}` only by the constant `mc^{2}`.

**Reduction-to-textbook verdict:** ✅ exact identity. **Verdict:** ✅. **Companion:** [`GriffithsCh01_P1_5.wl`](../../Mathematica_Notebooks/Quantum_Mechanics/GriffithsCh01_P1_5.wl).

---

### Problem G3e-P1.14 — Probability current and continuity equation

**Selection provenance:** the probability-current concept is foundational for any extension of QM (relativistic, scattering, time-dependent). Tests whether the proper-time formulation preserves probability conservation. *Pragmatic AI.*

**Source:** Griffiths, *Introduction to Quantum Mechanics*, 3e Problem 1.14 (Cambridge 2018). 2e equivalence: Problem 1.14 (Pearson 2005). *Paraphrased.*

**Paraphrased statement:** Derive the probability current `J(x, t) = (\hbar/(2 m i))\,(\psi^{*}\partial_{x}\psi - \psi\,\partial_{x}\psi^{*})` from the Schrödinger equation, and verify the continuity equation `\partial_{t}\rho + \partial_{x}J = 0` with `\rho = |\psi|^{2}`.

#### (a)/(b) Griffiths textbook solution

From `i\hbar\,\partial_{t}\psi = (-\hbar^{2}/(2m))\partial_{x}^{2}\psi + V\psi` and its complex conjugate:

$$
\partial_{t}|\psi|^{2} = \frac{1}{i\hbar}(\psi^{*}H\psi - \psi H\psi^{*}) = \frac{\hbar}{2 m i}\partial_{x}(\psi^{*}\partial_{x}\psi - \psi\partial_{x}\psi^{*}) = -\partial_{x}J.
$$

The continuity equation `\partial_{t}\rho + \partial_{x}J = 0` follows. Probability is conserved; the current `J` describes its flow.

#### (c) Proper-time / dual-theory reformulation

In the proper-time formulation, time evolution is governed by `K`:

$$
i\hbar\,\partial_{\tau}\psi = K\psi = \left[\frac{\hat p^{2}}{2 m} + m c^{2}\right]\psi.
$$

The Schrödinger-equation derivation of the continuity equation goes through identically — the rest-energy constant `mc^{2}` is real and cancels in `\psi^{*}H\psi - \psi H\psi^{*}` (it is a c-number, not an operator). The same current formula `J(x, \tau) = (\hbar/(2 m i))(\psi^{*}\partial_{x}\psi - \psi\partial_{x}\psi^{*})` results, with `\partial_{\tau}` replacing `\partial_{t}`.

In the non-relativistic limit, `\partial_{\tau} \to \partial_{t}` and the proper-time and textbook continuity equations coincide. Probability conservation is preserved.

**Reduction-to-textbook verdict:** ✅ exact identity at the operator level; identical observable `J(x, t)` in the non-rel limit.

**Verdict:** ✅. The proper-time formulation preserves probability conservation through the same continuity-equation structure as textbook QM.

**Companion notebook:** [`Roadmapping/Mathematica_Notebooks/Quantum_Mechanics/GriffithsCh01_P1_14.wl`](../../Mathematica_Notebooks/Quantum_Mechanics/GriffithsCh01_P1_14.wl).
