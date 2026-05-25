# Ch. 5 — Identical Particles

**PR E.** Bose / Fermi statistics under the proper-time formulation: no expected divergence from textbook, since (anti)symmetrisation is intrinsic to the wavefunction structure and not affected by which Hamiltonian governs time evolution.

## Problems

| Problem | Status | Role |
|---|---|---|
| [G3e-P5.1 — Two non-interacting identical particles in a box](#problem-g3e-p51--two-non-interacting-identical-particles-in-a-box) | drafted | fluency |
| [G3e-P5.5 — Periodic table (electron shell-filling)](#problem-g3e-p55--periodic-table-electron-shell-filling) | drafted | fluency |
| [G3e-P5.18 — Fermi sea and free-electron gas](#problem-g3e-p518--fermi-sea-and-free-electron-gas) | drafted | fluency |

---

### Problem G3e-P5.1 — Two non-interacting identical particles in a box

**Source:** Griffiths 3e Problem 5.1. *Pragmatic AI.*

**Statement:** Two identical particles (bosons or fermions) in a 1D infinite square well. Compute the ground-state energy and wavefunction.

**Textbook:** Bosons: `\psi_{B} = \psi_{1}(x_{1})\psi_{1}(x_{2})`, `E = 2E_{1}`. Fermions: anti-symmetrised, `\psi_{F} = (1/\sqrt 2)[\psi_{1}(x_{1})\psi_{2}(x_{2}) - \psi_{2}(x_{1})\psi_{1}(x_{2})]`, `E = E_{1} + E_{2}`. The fermion ground state is higher than the boson's because Pauli forbids both fermions occupying `\psi_{1}`.

**Proper-time:** (Anti)symmetrisation is intrinsic to multi-particle wavefunctions; the proper-time `K` does not modify the (anti)symmetry requirement. Same ground-state structure and energies (offset by `2 m c^{2}` rest energy of the two particles).

**Verdict:** ✅. **Companion:** [`GriffithsCh05_P5_1.wl`](../../Mathematica_Notebooks/Quantum_Mechanics/GriffithsCh05_P5_1.wl).

---

### Problem G3e-P5.5 — Periodic table (electron shell-filling)

**Source:** Griffiths 3e Problem 5.5. *Pragmatic AI.*

**Statement:** Apply the Pauli exclusion principle to fill electron shells in hydrogen-like atoms. Identify the noble-gas configurations.

**Textbook:** Aufbau principle; shells fill in order `1s^{2}, 2s^{2} 2p^{6}, 3s^{2} 3p^{6}, \ldots` giving noble-gas configurations at `Z = 2, 10, 18, 36, \ldots`.

**Proper-time:** Shell-filling depends on single-electron energy ordering (Bohr + spin–orbit corrections) and Pauli exclusion. Bohr ordering is unchanged in proper-time per [Problem G3e-P4.11](Ch04_QM_3D.md#problem-g3e-p411--bohr-spectrum-from-proper-time-k). Same shell-filling sequence; same periodic table.

**Verdict:** ✅. **Companion:** [`GriffithsCh05_P5_5.wl`](../../Mathematica_Notebooks/Quantum_Mechanics/GriffithsCh05_P5_5.wl).

---

### Problem G3e-P5.18 — Fermi sea and free-electron gas

**Source:** Griffiths 3e Problem 5.18. *Pragmatic AI.*

**Statement:** Compute the Fermi energy of a free-electron gas at density `n`.

**Textbook:** `E_{F} = (\hbar^{2}/2m)(3\pi^{2}n)^{2/3}`. Density of states `g(E) \propto \sqrt E` for non-relativistic 3D free electrons.

**Proper-time:** Same single-electron kinetic energy `p^{2}/(2m)`; same Fermi sphere; same Fermi energy formula.

**Verdict:** ✅. **Companion:** [`GriffithsCh05_P5_18.wl`](../../Mathematica_Notebooks/Quantum_Mechanics/GriffithsCh05_P5_18.wl).
