# Ch. 3 — Formalism

Part of **PR C** per [§4 of the plan](../../../.dev/tasks/49-quantum-mechanics-griffiths-proper-time.md#4-pr-sequencing-12-prs-al-per-issue-body). Hilbert-space recap and generalised statistical interpretation. Proper-time formulation preserves all operator algebra; the only difference is `K = p²/(2m) + mc²` as the Hamiltonian. Three null-result fluency-builders.

## Problems

| Problem | Status | Role |
|---|---|---|
| [G3e-P3.13 — Hermitian operators and observables](#problem-g3e-p313--hermitian-operators-and-observables) | drafted | fluency |
| [G3e-P3.15 — Position-momentum commutator and uncertainty](#problem-g3e-p315--position-momentum-commutator-and-uncertainty) | drafted | fluency |
| [G3e-P3.21 — Dirac notation and projection operators](#problem-g3e-p321--dirac-notation-and-projection-operators) | drafted | fluency |

---

### Problem G3e-P3.13 — Hermitian operators and observables

**Source:** Griffiths 3e Problem 3.13. *Pragmatic AI.*

**Statement:** Show that physical observables are represented by Hermitian operators (eigenvalues real, eigenvectors of different eigenvalues orthogonal).

**Textbook:** Standard `\langle\phi|A\psi\rangle = \langle A\phi|\psi\rangle` for Hermitian `A`; eigenvalue/eigenvector orthogonality follows.

**Proper-time:** Hermiticity is a property of operators, not Hamiltonian-dependent. `\hat x`, `\hat p`, and `K` are all Hermitian; same theorem applies.

**Verdict:** ✅. **Companion:** [`GriffithsCh03_P3_13.wl`](../../Mathematica_Notebooks/Quantum_Mechanics/GriffithsCh03_P3_13.wl).

---

### Problem G3e-P3.15 — Position-momentum commutator and uncertainty

**Source:** Griffiths 3e Problem 3.15. *Pragmatic AI.*

**Statement:** Derive `[\hat x, \hat p] = i\hbar` and apply the generalised uncertainty principle `\Delta A\,\Delta B \ge (1/2)|\langle[A, B]\rangle|`.

**Textbook:** Canonical commutation relation; minimum uncertainty `\Delta x\,\Delta p \ge \hbar/2`.

**Proper-time:** Position and momentum operators unmodified by the proper-time substitution. `[\hat x, \hat p] = i\hbar` unchanged. Generalised uncertainty principle holds identically.

**Verdict:** ✅. **Companion:** [`GriffithsCh03_P3_15.wl`](../../Mathematica_Notebooks/Quantum_Mechanics/GriffithsCh03_P3_15.wl).

---

### Problem G3e-P3.21 — Dirac notation and projection operators

**Source:** Griffiths 3e Problem 3.21. *Pragmatic AI.*

**Statement:** Express the completeness relation `\sum_{n}|n\rangle\langle n| = \hat 1` and verify it for the harmonic-oscillator eigenstates.

**Textbook:** Standard Dirac formalism; completeness gives the resolution of identity.

**Proper-time:** Eigenstates of `K` differ from eigenstates of `\hat H_{\text{NR}}` only by an irrelevant constant (`K = \hat H_{\text{NR}} + mc^{2}\hat 1`), so the eigenstate basis is identical and the completeness relation is unchanged.

**Verdict:** ✅. **Companion:** [`GriffithsCh03_P3_21.wl`](../../Mathematica_Notebooks/Quantum_Mechanics/GriffithsCh03_P3_21.wl).
