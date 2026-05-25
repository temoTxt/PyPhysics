# Ch. 9 — WKB Approximation

**PR I.** Semiclassical (WKB) approximation: tunnelling, classical-turning-point matching, Bohr–Sommerfeld quantisation. Non-relativistic; proper-time `K` reduces exactly. Three fluency-builders.

## Problems

| Problem | Status | Role |
|---|---|---|
| [G3e-P9.1 — WKB connection formula](#problem-g3e-p91--wkb-connection-formula) | drafted | fluency |
| [G3e-P9.7 — Tunnelling through a square barrier](#problem-g3e-p97--tunnelling-through-a-square-barrier) | drafted | fluency |
| [G3e-P9.14 — Bohr–Sommerfeld quantisation](#problem-g3e-p914--bohrsommerfeld-quantisation) | drafted | fluency |

---

### Problem G3e-P9.1 — WKB connection formula

**Source:** Griffiths 3e Problem 9.1. *Pragmatic AI.*

**Statement:** Derive the WKB wavefunction in classically-allowed and classically-forbidden regions, and the connection formula at classical turning points.

**Textbook:** Classically allowed: `\psi \approx (1/\sqrt{p(x)})[A\cos(\int p/\hbar\,dx) + B\sin(\int p/\hbar\,dx)]`. Classically forbidden: exponentially decaying / growing.

**Proper-time:** WKB analysis is a semiclassical expansion in `\hbar`; the proper-time `K - mc^{2} = p^{2}/(2m)` is identical to the non-relativistic textbook kinetic energy, so the WKB expansion is unchanged.

**Verdict:** ✅. **Companion:** [`GriffithsCh09_P9_1.wl`](../../Mathematica_Notebooks/Quantum_Mechanics/GriffithsCh09_P9_1.wl).

---

### Problem G3e-P9.7 — Tunnelling through a square barrier

**Source:** Griffiths 3e Problem 9.7. *Pragmatic AI.*

**Statement:** Compute the WKB tunnelling probability through a square barrier `V_{0}` of width `a`.

**Textbook:** Transmission coefficient `T \approx e^{-2\kappa a}` with `\kappa = \sqrt{2m(V_{0} - E)}/\hbar`. Exponential suppression of tunnelling probability with barrier width.

**Proper-time:** Same `\kappa`, same `T`. The relativistic correction to `\kappa` from the next-order `K` expansion would be `O((v/c)^{2})`-small and is well below typical observational precision for tunnelling experiments.

**Verdict:** ✅. **Companion:** [`GriffithsCh09_P9_7.wl`](../../Mathematica_Notebooks/Quantum_Mechanics/GriffithsCh09_P9_7.wl).

---

### Problem G3e-P9.14 — Bohr–Sommerfeld quantisation

**Source:** Griffiths 3e Problem 9.14. *Pragmatic AI.*

**Statement:** Derive the Bohr–Sommerfeld quantisation rule `\oint p\,dx = (n + 1/2)\,2\pi\hbar` from WKB connection formulae.

**Textbook:** Standard semiclassical derivation; recovers harmonic-oscillator spectrum `E_{n} = \hbar\omega(n + 1/2)`.

**Proper-time:** Same semiclassical derivation; `p^{2}/(2m)` unchanged ⇒ same quantisation rule. `K`-eigenvalue offset is `mc^{2}` (constant), as throughout the campaign.

**Verdict:** ✅. **Companion:** [`GriffithsCh09_P9_14.wl`](../../Mathematica_Notebooks/Quantum_Mechanics/GriffithsCh09_P9_14.wl).
