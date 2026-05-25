# Ch. 6 — Symmetries and Conservation Laws

**PR F.** Symmetries and conserved quantities under the proper-time formulation. Continuous symmetries (translations, rotations) give conserved momenta exactly as in textbook QM; discrete symmetries (parity, time reversal) are unaffected by the choice of Hamiltonian. Three fluency-builders.

## Problems

| Problem | Status | Role |
|---|---|---|
| [G3e-P6.3 — Translation symmetry and momentum conservation](#problem-g3e-p63--translation-symmetry-and-momentum-conservation) | drafted | fluency |
| [G3e-P6.9 — Rotation symmetry and angular momentum](#problem-g3e-p69--rotation-symmetry-and-angular-momentum) | drafted | fluency |
| [G3e-P6.20 — Parity and time reversal](#problem-g3e-p620--parity-and-time-reversal) | drafted | fluency |

---

### Problem G3e-P6.3 — Translation symmetry and momentum conservation

**Source:** Griffiths 3e Problem 6.3. *Pragmatic AI.*

**Statement:** Show that translation invariance of the Hamiltonian implies momentum conservation.

**Textbook:** `[\hat H, \hat p] = 0` ⇒ `\langle\hat p\rangle = const`. Translation operator `e^{i\hat p a/\hbar}` commutes with `\hat H` for translationally-invariant systems.

**Proper-time:** `[\hat K, \hat p] = [\hat p^{2}/(2m) + mc^{2}, \hat p] = 0` (free particle). Same conservation law.

**Verdict:** ✅. **Companion:** [`GriffithsCh06_P6_3.wl`](../../Mathematica_Notebooks/Quantum_Mechanics/GriffithsCh06_P6_3.wl).

---

### Problem G3e-P6.9 — Rotation symmetry and angular momentum

**Source:** Griffiths 3e Problem 6.9. *Pragmatic AI.*

**Statement:** Show rotation invariance of the Hamiltonian implies angular momentum conservation.

**Textbook:** `[\hat H, \hat L_{z}] = 0` for central potential. Rotation operator `e^{i\hat L_{z}\phi/\hbar}` commutes with central-potential `\hat H`.

**Proper-time:** Angular momentum operator unchanged; `\hat K` for central potential commutes with `\hat L_{z}` (the central nature of the potential is unaffected by the proper-time `K`).

**Verdict:** ✅. **Companion:** [`GriffithsCh06_P6_9.wl`](../../Mathematica_Notebooks/Quantum_Mechanics/GriffithsCh06_P6_9.wl).

---

### Problem G3e-P6.20 — Parity and time reversal

**Source:** Griffiths 3e Problem 6.20. *Pragmatic AI.*

**Statement:** Discuss parity and time-reversal symmetries of the Hamiltonian. Find which eigenstates have definite parity.

**Textbook:** Eigenstates of `\hat H = \hat p^{2}/(2m) + V(\hat x)` have definite parity if `V(-x) = V(x)` (even potential).

**Proper-time:** The proper-time `K = \hat p^{2}/(2m) + V(\hat x) + mc^{2}` has the same parity property as textbook `\hat H`. Time-reversal symmetry preserved because the additive constant `mc^{2}` does not affect anti-unitary structure.

**Verdict:** ✅. **Companion:** [`GriffithsCh06_P6_20.wl`](../../Mathematica_Notebooks/Quantum_Mechanics/GriffithsCh06_P6_20.wl).
