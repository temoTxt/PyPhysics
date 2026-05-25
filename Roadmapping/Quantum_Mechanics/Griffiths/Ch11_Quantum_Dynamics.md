# Ch. 11 — Quantum Dynamics (3e renumber)

**PR K.** Time-dependent perturbation theory; emission/absorption; sudden and adiabatic approximations. Most problems are non-relativistic; emission/absorption couples to EM radiation and bridges to PR D of #42 (Liénard–Wiechert) at the relativistic edge. Four problems.

Note: Ch. 11 (3e) was renumbered between editions; the 2e equivalent is approximately Ch. 9 in some printings. Per-problem source lines record both edition numbers.

## Problems

| Problem | Status | Role |
|---|---|---|
| [G3e-P11.5 — Time-dependent perturbation theory (TDPT)](#problem-g3e-p115--time-dependent-perturbation-theory-tdpt) | drafted | fluency |
| [G3e-P11.10 — Stimulated emission and absorption](#problem-g3e-p1110--stimulated-emission-and-absorption) | drafted | fluency |
| [G3e-P11.15 — Sudden approximation](#problem-g3e-p1115--sudden-approximation) | drafted | fluency |
| [G3e-P11.20 — Spontaneous emission rate (Einstein A)](#problem-g3e-p1120--spontaneous-emission-rate-einstein-a) | drafted | headline-adjacent |

---

### Problem G3e-P11.5 — Time-dependent perturbation theory (TDPT)

**Source:** Griffiths 3e Problem 11.5. *Pragmatic AI.*

**Statement:** Derive the first-order TDPT transition probability `P_{i \to f}(t) = |\langle f|H'(t)|i\rangle/\hbar|^{2}\,\sin^{2}((\omega_{f i} - \omega)t/2)/((\omega_{f i} - \omega)/2)^{2}` for a sinusoidal perturbation.

**Textbook:** Standard Fermi golden rule derivation.

**Proper-time:** TDPT formalism depends on the time-evolution structure `i\hbar\,\partial_{t}|\psi\rangle = H|\psi\rangle`. In proper-time, evolution is `i\hbar\,\partial_{\tau}|\psi\rangle = K|\psi\rangle`. In the non-relativistic limit `\partial_{\tau} \to \partial_{t}` and `K \to \hat H_{\text{NR}} + mc^{2}`; TDPT structure preserved.

**Verdict:** ✅. **Companion:** [`GriffithsCh11_P11_5.wl`](../../Mathematica_Notebooks/Quantum_Mechanics/GriffithsCh11_P11_5.wl).

---

### Problem G3e-P11.10 — Stimulated emission and absorption

**Source:** Griffiths 3e Problem 11.10. *Pragmatic AI.*

**Statement:** Derive the Einstein B coefficient for stimulated emission/absorption of an atom in a thermal radiation field.

**Textbook:** `B_{f i} = (\pi e^{2}/(3\varepsilon_{0}\hbar^{2}))|\langle f|\mathbf{r}|i\rangle|^{2}` (electric dipole). Connects to A coefficient via `A = (\hbar\omega^{3}/(\pi^{2}c^{3})) B`.

**Proper-time:** Electric-dipole approximation involves the position operator (unchanged) and the incident EM field (treated semiclassically; same as textbook). Result identical.

**Verdict:** ✅. **Companion:** [`GriffithsCh11_P11_10.wl`](../../Mathematica_Notebooks/Quantum_Mechanics/GriffithsCh11_P11_10.wl).

---

### Problem G3e-P11.15 — Sudden approximation

**Source:** Griffiths 3e Problem 11.15. *Pragmatic AI.*

**Statement:** Apply the sudden approximation to compute transition probabilities under a rapid change of the Hamiltonian.

**Textbook:** `|\psi_{\text{after}}\rangle = |\psi_{\text{before}}\rangle` (state unchanged in sudden limit); probability of finding in new eigenstate `|n_{\text{new}}\rangle` is `|\langle n_{\text{new}}|\psi_{\text{before}}\rangle|^{2}`.

**Proper-time:** Sudden approximation is a kinematic argument about the state vector under rapid Hamiltonian change; independent of which Hamiltonian (NR vs proper-time) governs the evolution. Result identical.

**Verdict:** ✅. **Companion:** [`GriffithsCh11_P11_15.wl`](../../Mathematica_Notebooks/Quantum_Mechanics/GriffithsCh11_P11_15.wl).

---

### Problem G3e-P11.20 — Spontaneous emission rate (Einstein A)

**Selection provenance:** spontaneous emission connects QM directly to classical radiation (via Larmor's formula at the classical analogue limit, [J3e-P14.3 in #42](../../Electromagnetism/Jackson/Ch14_Radiation_by_Moving_Charges.md#problem-j3e-p143--non-relativistic-larmor-radiation-formula)). The Einstein A coefficient should be reproduced identically.

**Source:** Griffiths 3e Problem 11.20.

**Statement:** Compute the spontaneous-emission rate (Einstein A coefficient) for the `2p \to 1s` transition in hydrogen.

**Textbook:** `A_{2p \to 1s} = (\omega^{3}/(\pi\varepsilon_{0}\hbar c^{3}))\,|\langle 1s|e\mathbf{r}|2p\rangle|^{2} \approx 6.27 \times 10^{8}` s⁻¹. Lifetime `\tau \approx 1.6` ns.

**Proper-time:** Spontaneous emission rate is computed by quantising the EM radiation field (semiclassically: combination of stimulated emission + vacuum-fluctuation argument). The QED computation involves the dual Dirac equation if pursued from first principles, but at the level of Griffiths' textbook derivation (matrix elements of `e\mathbf{r}`), the result is unchanged.

**Cross-reference to #42**: this Einstein A coefficient is the classical analogue of the radiation-damping rate `\Gamma = \tau_{0}\omega_{0}^{2}` computed in [Problem J3e-P16.2](../../Electromagnetism/Jackson/Ch16_Radiation_Damping.md#problem-j3e-p162--radiation-reaction-damping-of-a-harmonic-oscillator). The two campaigns produce consistent estimates for the natural linewidth of atomic transitions.

**Verdict:** ✅. **Companion:** [`GriffithsCh11_P11_20.wl`](../../Mathematica_Notebooks/Quantum_Mechanics/GriffithsCh11_P11_20.wl).
