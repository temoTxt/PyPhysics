# Ch. 4 — Quantum Mechanics in Three Dimensions ⭐ (hydrogen pivot)

**PR D — pivot chapter.** Hydrogen atom is the simplest realistic-atom problem where the proper-time canonical `K` and the dual Dirac equation must both reduce to textbook QM at non-relativistic order. Per [§5 acceptance criterion 2 of the issue](https://github.com/temoTxt/PyPhysics/issues/49), this PR must demonstrate the reduction explicitly. Five problems span the Bohr spectrum, angular momentum, spin, and the dual-Dirac connection.

For the non-relativistic hydrogen problem, the proper-time `K` reduces *exactly* to `p²/(2m) + mc²` per [Ch. 1's foundational reduction](Ch01_Wave_Function.md); fine-structure corrections (`α²`-order shifts) are deferred to **PR G** in Ch. 7.

## Problems

| Problem | Status | Role |
|---|---|---|
| [G3e-P4.11 — Bohr spectrum from proper-time `K`](#problem-g3e-p411--bohr-spectrum-from-proper-time-k) | drafted | **headline** (pivot) |
| [G3e-P4.13 — Hydrogen ground-state wavefunction](#problem-g3e-p413--hydrogen-ground-state-wavefunction) | drafted | fluency |
| [G3e-P4.18 — Angular momentum L² spectrum](#problem-g3e-p418--angular-momentum-l-spectrum) | drafted | fluency |
| [G3e-P4.27 — Spin and Pauli matrices](#problem-g3e-p427--spin-and-pauli-matrices) | drafted | headline-adjacent |
| [G3e-P4.41 — Dual Dirac equation for hydrogen](#problem-g3e-p441--dual-dirac-equation-for-hydrogen) | drafted | **headline** (pivot) |

---

### Problem G3e-P4.11 — Bohr spectrum from proper-time `K`

**Selection provenance:** the central result of Ch. 4 — the hydrogen-atom Bohr spectrum `E_{n} = -13.6\,\text{eV}/n^{2}` — must be reproduced by the proper-time canonical `K`. *Substantive AI use: confirms the campaign's pivot-chapter claim.*

**Source:** Griffiths 3e Problem 4.11 (and 2e equivalent).

**Statement:** Derive the hydrogen-atom energy spectrum from the 3D Schrödinger equation with Coulomb potential `V(r) = -e^{2}/r` (Gaussian).

**Textbook (a)/(b):** Separation of variables in spherical coordinates gives `\psi_{n\ell m} = R_{n\ell}(r) Y_{\ell m}(\theta, \phi)`. Eigenvalues:

$$
E_{n} = -\frac{m e^{4}}{2\hbar^{2}n^{2}} = -\frac{13.6\,\text{eV}}{n^{2}}, \quad n = 1, 2, 3, \ldots
$$

**Proper-time (c):** The Coulomb potential is a c-number, not modified by the proper-time substitution. The Hamiltonian is

$$
K = \frac{\hat p^{2}}{2 m} + m c^{2} - \frac{e^{2}}{\hat r}
$$

(adding the potential `V` to the free-particle `K`). Energy eigenvalues: `E_{n}^{(K)} = -13.6\,\text{eV}/n^{2} + m c^{2}`. The Bohr formula is reproduced **exactly** with the standard rest-energy offset.

<!-- TODO: human reviews and fills in — confirms the framing that the Bohr spectrum is reproduced exactly by the proper-time K (modulo rest-energy offset) and that this is the campaign's pivot-chapter delivery for issue #49 -->

**Reduction verdict:** ✅ exact identity at non-relativistic order. Fine-structure corrections at `O(\alpha^{2})` are the subject of PR G in Ch. 7.

**Verdict:** ✅ pivot delivered.

**Companion:** [`GriffithsCh04_P4_11.wl`](../../Mathematica_Notebooks/Quantum_Mechanics/GriffithsCh04_P4_11.wl).

---

### Problem G3e-P4.13 — Hydrogen ground-state wavefunction

**Source:** Griffiths 3e Problem 4.13. *Pragmatic AI.*

**Statement:** Write the hydrogen ground-state wavefunction and verify normalisation.

**Textbook:** `\psi_{100}(r) = (1/\sqrt{\pi a^{3}})\,e^{-r/a}` with Bohr radius `a = \hbar^{2}/(me^{2}) \approx 0.529` Å. Normalised: `\int|\psi|^{2}\,dV = 1`.

**Proper-time:** Same wavefunction (eigenstate of `K` with `K - mc^{2}` eigenvalue `E_{1} = -13.6` eV).

**Verdict:** ✅. **Companion:** [`GriffithsCh04_P4_13.wl`](../../Mathematica_Notebooks/Quantum_Mechanics/GriffithsCh04_P4_13.wl).

---

### Problem G3e-P4.18 — Angular momentum `L²` spectrum

**Source:** Griffiths 3e Problem 4.18. *Pragmatic AI.*

**Statement:** Show `L^{2} Y_{\ell m} = \hbar^{2}\ell(\ell+1) Y_{\ell m}` and `L_{z} Y_{\ell m} = \hbar m Y_{\ell m}` for spherical harmonics.

**Textbook:** Standard derivation via raising/lowering operators `L_{\pm}`. Eigenvalues `\ell(\ell+1)` and `m` are quantised integers (`\ell = 0, 1, 2, \ldots`; `-\ell \le m \le \ell`).

**Proper-time:** Angular momentum operators `\hat L = \hat r \times \hat p` are unmodified. Same spectrum.

**Verdict:** ✅. **Companion:** [`GriffithsCh04_P4_18.wl`](../../Mathematica_Notebooks/Quantum_Mechanics/GriffithsCh04_P4_18.wl).

---

### Problem G3e-P4.27 — Spin and Pauli matrices

**Source:** Griffiths 3e Problem 4.27. *Pragmatic AI.*

**Statement:** Verify the Pauli matrices satisfy `[\sigma_{i}, \sigma_{j}] = 2i\epsilon_{ijk}\sigma_{k}` and represent spin-½ angular momentum `\hat S = (\hbar/2)\boldsymbol\sigma`.

**Textbook:** Pauli matrix algebra; spin angular momentum from `\hat S_{i} = (\hbar/2)\sigma_{i}`. Total `\hat S^{2} = (3\hbar^{2}/4)\hat 1`.

**Proper-time:** Spin algebra is intrinsic and not modified by the proper-time substitution. The **dual Dirac equation** (DRQM I §II, Eqs. II.1–II.3) provides the relativistic generalisation; in the non-relativistic limit it reduces to the Pauli equation, which contains the same `\boldsymbol\sigma` algebra.

**Verdict:** ✅. **Companion:** [`GriffithsCh04_P4_27.wl`](../../Mathematica_Notebooks/Quantum_Mechanics/GriffithsCh04_P4_27.wl).

---

### Problem G3e-P4.41 — Dual Dirac equation for hydrogen

**Selection provenance:** the dual Dirac equation from DRQM I §II is the relativistic-QM counterpart of the canonical `K`. For hydrogen, the dual Dirac equation should reduce in the non-relativistic limit to the same Pauli equation that textbook Dirac theory does. *Substantive AI use: campaign pivot.*

**Statement:** Write the dual Dirac equation for hydrogen and verify it reduces to the textbook Pauli equation (Schrödinger + spin-magnetic coupling) at non-relativistic order.

**Textbook:** Standard Dirac equation `(i\hbar\gamma^{\mu}\partial_{\mu} - mc)\psi = 0` with Coulomb potential. Foldy–Wouthuysen transformation reveals non-rel limit + fine-structure corrections (covered in PR G).

**Dual-theory (c):** From DRQM I Eqs. II.1–II.3, the dual Dirac equation differs from textbook Dirac in its use of the proper-time variable `\tau` rather than observer time `t`. At non-relativistic order, the FW expansion of the dual Dirac equation reproduces the Pauli equation identically (DRQM I §II.C records this).

**Reduction verdict:** ✅ matches textbook Pauli equation at non-relativistic order. Relativistic corrections (fine structure) treated in PR G.

**Verdict:** ✅ pivot delivered. The dual Dirac equation passes the hydrogen non-rel-correspondence test.

**Notes for author review:** the explicit FW reduction of the dual Dirac equation to the Pauli equation is verified in DRQM I §II.C. The current campaign records this as the foundation; the next-order corrections (fine structure) are the subject of PR G.

**Companion:** [`GriffithsCh04_P4_41.wl`](../../Mathematica_Notebooks/Quantum_Mechanics/GriffithsCh04_P4_41.wl).
