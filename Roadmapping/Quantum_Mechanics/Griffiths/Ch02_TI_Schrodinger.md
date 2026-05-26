# Ch. 2 — Time-Independent Schrödinger Equation

Part of **PR B** per [§4 of the plan](../../../.dev/tasks/49-quantum-mechanics-griffiths-proper-time.md#4-pr-sequencing-12-prs-al-per-issue-body). Bound-state problems with non-relativistic spectra; proper-time `K` reduces *exactly* to `p²/(2m) + mc²` per [Ch. 1's foundational reduction](Ch01_Wave_Function.md). All four problems are pure null-result fluency-builders.

## Problems

| Problem | Status | Role |
|---|---|---|
| [G3e-P2.4 — Infinite square well](#problem-g3e-p24--infinite-square-well) | drafted | fluency |
| [G3e-P2.10 — Harmonic oscillator](#problem-g3e-p210--harmonic-oscillator) | drafted | fluency |
| [G3e-P2.23 — Delta-function well](#problem-g3e-p223--delta-function-well) | drafted | fluency |
| [G3e-P2.32 — Finite square well](#problem-g3e-p232--finite-square-well) | drafted | fluency |

---

### Problem G3e-P2.4 — Infinite square well

**Source:** Griffiths 3e Problem 2.4 (2e equivalent). *Pragmatic AI.*

**Statement:** Particle of mass `m` confined to `0 < x < a` by infinite potential walls. Find the energy eigenvalues and eigenfunctions.

**Textbook:** `\psi_{n}(x) = \sqrt{2/a}\sin(n\pi x/a)`, `E_{n} = n^{2}\pi^{2}\hbar^{2}/(2 m a^{2})` for `n = 1, 2, 3, \ldots`.

**Proper-time:** `K = p^{2}/(2m) + mc^{2}` exactly; eigenfunctions unchanged; eigenvalues shifted by `+mc^{2}` (rest-energy offset).

**Verdict:** ✅ identical eigenfunctions; eigenvalues differ by constant `mc^{2}`. **Companion:** [`GriffithsCh02_P2_4.wl`](../../Mathematica_Notebooks/Quantum_Mechanics/GriffithsCh02_P2_4.wl).

---

### Problem G3e-P2.10 — Harmonic oscillator

**Source:** Griffiths 3e Problem 2.10. *Pragmatic AI.*

**Statement:** Particle in potential `V(x) = (1/2)m\omega^{2}x^{2}`. Find the energy spectrum.

**Textbook:** `E_{n} = \hbar\omega(n + 1/2)` for `n = 0, 1, 2, \ldots`; eigenfunctions are Hermite-polynomial weighted Gaussians.

**Proper-time:** Same eigenfunctions; eigenvalues `E_{n} + mc^{2}`. The harmonic oscillator's ladder structure is preserved by the proper-time formulation.

**Verdict:** ✅. **Companion:** [`GriffithsCh02_P2_10.wl`](../../Mathematica_Notebooks/Quantum_Mechanics/GriffithsCh02_P2_10.wl).

---

### Problem G3e-P2.23 — Delta-function well

**Source:** Griffiths 3e Problem 2.23. *Pragmatic AI.*

**Statement:** Particle in `V(x) = -\alpha\delta(x)`. Find the bound-state energy and wavefunction.

**Textbook:** Single bound state with `E = -m\alpha^{2}/(2\hbar^{2})`, `\psi(x) = \sqrt{m\alpha}/\hbar \cdot e^{-m\alpha|x|/\hbar^{2}}`.

**Proper-time:** Same wavefunction; energy shifted by `+mc^{2}`. The delta-well's single bound state is preserved.

**Verdict:** ✅. **Companion:** [`GriffithsCh02_P2_23.wl`](../../Mathematica_Notebooks/Quantum_Mechanics/GriffithsCh02_P2_23.wl).

---

### Problem G3e-P2.32 — Finite square well

**Source:** Griffiths 3e Problem 2.32. *Pragmatic AI.*

**Statement:** Particle in `V(x) = -V_{0}` for `|x| < a`, `V = 0` outside. Find the bound-state spectrum (transcendental equation).

**Textbook:** Bound states determined by `\tan(ka) = \sqrt{(V_{0} - E)/E}` (even) or `\cot(ka) = -\sqrt{(V_{0} - E)/E}` (odd), with `k = \sqrt{2m(V_{0} + E)}/\hbar`. Number of bound states grows with `V_{0}a^{2}`.

**Proper-time:** Same transcendental equation; energies shifted by `+mc^{2}`.

**Verdict:** ✅. **Companion:** [`GriffithsCh02_P2_32.wl`](../../Mathematica_Notebooks/Quantum_Mechanics/GriffithsCh02_P2_32.wl).
