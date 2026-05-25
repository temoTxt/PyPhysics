# Ch. 8 — Variational Principle

**PR H.** Variational principle and applications to many-body systems. Non-relativistic; proper-time `K` reduces exactly. Three fluency-builders.

## Problems

| Problem | Status | Role |
|---|---|---|
| [G3e-P8.1 — Variational principle and the ground-state bound](#problem-g3e-p81--variational-principle-and-the-ground-state-bound) | drafted | fluency |
| [G3e-P8.10 — Helium ground-state energy](#problem-g3e-p810--helium-ground-state-energy) | drafted | fluency |
| [G3e-P8.14 — Hydrogen-molecule ion `H_{2}^{+}`](#problem-g3e-p814--hydrogen-molecule-ion-h_2) | drafted | fluency |

---

### Problem G3e-P8.1 — Variational principle and the ground-state bound

**Source:** Griffiths 3e Problem 8.1. *Pragmatic AI.*

**Statement:** Show that for any trial wavefunction `\psi_{\text{trial}}`, the expectation value `\langle\psi_{\text{trial}}|H|\psi_{\text{trial}}\rangle \ge E_{0}` (the true ground-state energy).

**Textbook:** Proof via expansion of `\psi_{\text{trial}}` in `H`-eigenstates.

**Proper-time:** Proof depends only on Hermiticity of `H` and orthonormality of eigenstates. `K` is Hermitian with orthonormal eigenstates (identical to textbook's, modulo `mc^{2}` energy shift). Same theorem.

**Verdict:** ✅. **Companion:** [`GriffithsCh08_P8_1.wl`](../../Mathematica_Notebooks/Quantum_Mechanics/GriffithsCh08_P8_1.wl).

---

### Problem G3e-P8.10 — Helium ground-state energy

**Source:** Griffiths 3e Problem 8.10. *Pragmatic AI.*

**Statement:** Estimate the helium ground-state energy using a trial wavefunction with effective nuclear charge `Z' < 2`.

**Textbook:** Trial: `\psi_{\text{trial}}(r_{1}, r_{2}) = (Z'^{3}/\pi a^{3}) e^{-Z'(r_{1}+r_{2})/a}`. Minimising `\langle H \rangle` with respect to `Z'` gives `Z'_{\min} = 27/16 \approx 1.69` and `E_{\text{trial}} \approx -77.5` eV (experimentally `-79.0` eV; 2% error).

**Proper-time:** Same trial functions; same expectation values (operators unchanged, integrals identical). Energy offset by `2 m_{e} c^{2}` (two electrons), which cancels in the variational ratio. `Z'_{\min}` and binding energy identical.

**Verdict:** ✅. **Companion:** [`GriffithsCh08_P8_10.wl`](../../Mathematica_Notebooks/Quantum_Mechanics/GriffithsCh08_P8_10.wl).

---

### Problem G3e-P8.14 — Hydrogen-molecule ion `H_{2}^{+}`

**Source:** Griffiths 3e Problem 8.14. *Pragmatic AI.*

**Statement:** Estimate the binding energy and bond length of `H_{2}^{+}` using LCAO trial functions.

**Textbook:** Linear combination of atomic orbitals `\psi_{\pm}(\mathbf{r}) = N[\psi_{1s}(\mathbf{r} - \mathbf{R}/2) \pm \psi_{1s}(\mathbf{r} + \mathbf{R}/2)]`. Bonding (`+`) gives binding energy `\sim 1.76` eV at bond length `\sim 2.5\,a_{0}` (experimentally `2.65` eV at `2.0\,a_{0}`).

**Proper-time:** Same LCAO structure; same integrals; same binding-energy result.

**Verdict:** ✅. **Companion:** [`GriffithsCh08_P8_14.wl`](../../Mathematica_Notebooks/Quantum_Mechanics/GriffithsCh08_P8_14.wl).
