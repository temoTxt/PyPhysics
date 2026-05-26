# Ch. 10 — Scattering

**PR J.** Scattering theory: partial waves, Born approximation, relativistic-kinematic scattering. The Born-approximation problem is non-relativistic; the relativistic-kinematics problems engage the proper-time formulation more substantively. Four problems.

## Problems

| Problem | Status | Role |
|---|---|---|
| [G3e-P10.1 — Partial-wave expansion](#problem-g3e-p101--partial-wave-expansion) | drafted | fluency |
| [G3e-P10.7 — Born approximation for Coulomb scattering](#problem-g3e-p107--born-approximation-for-coulomb-scattering) | drafted | fluency |
| [G3e-P10.10 — Hard-sphere scattering and s-wave phase shift](#problem-g3e-p1010--hard-sphere-scattering-and-s-wave-phase-shift) | drafted | fluency |
| [G3e-P10.18 — Relativistic-kinematics scattering](#problem-g3e-p1018--relativistic-kinematics-scattering) | drafted | headline-adjacent |

---

### Problem G3e-P10.1 — Partial-wave expansion

**Source:** Griffiths 3e Problem 10.1. *Pragmatic AI.*

**Statement:** Expand a scattered wave in partial waves and derive the relation between phase shifts `\delta_{\ell}` and the scattering amplitude `f(\theta)`.

**Textbook:** `f(\theta) = (1/k)\sum_{\ell}(2\ell+1)\sin\delta_{\ell}\,e^{i\delta_{\ell}}\,P_{\ell}(\cos\theta)`. Total cross-section `\sigma = (4\pi/k^{2})\sum(2\ell+1)\sin^{2}\delta_{\ell}` (optical theorem).

**Proper-time:** Partial-wave expansion is a kinematic property of the scattering wavefunction; the dispersion `k = \sqrt{2mE}/\hbar` is unchanged by the proper-time `K - mc^{2} = p^{2}/(2m)`.

**Verdict:** ✅. **Companion:** [`GriffithsCh10_P10_1.wl`](../../Mathematica_Notebooks/Quantum_Mechanics/GriffithsCh10_P10_1.wl).

---

### Problem G3e-P10.7 — Born approximation for Coulomb scattering

**Source:** Griffiths 3e Problem 10.7. *Pragmatic AI.*

**Statement:** Compute the Born-approximation scattering amplitude for the Coulomb potential `V(r) = q_{1}q_{2}/r`.

**Textbook:** `f(\theta) = -2 m q_{1} q_{2}/(\hbar^{2} q^{2})` with `q = 2 k \sin(\theta/2)`. Rutherford cross-section `d\sigma/d\Omega = (q_{1} q_{2})^{2}/(16 E^{2}\sin^{4}(\theta/2))`.

**Proper-time:** Born approximation uses non-rel `k`; result identical at leading order. Relativistic corrections (`Z\alpha`-order) are beyond Born approximation and require the dual Dirac equation; not engaged at this level.

**Verdict:** ✅. **Companion:** [`GriffithsCh10_P10_7.wl`](../../Mathematica_Notebooks/Quantum_Mechanics/GriffithsCh10_P10_7.wl).

---

### Problem G3e-P10.10 — Hard-sphere scattering and s-wave phase shift

**Source:** Griffiths 3e Problem 10.10. *Pragmatic AI.*

**Statement:** Compute the s-wave (`\ell = 0`) phase shift for hard-sphere scattering of radius `a`.

**Textbook:** `\delta_{0} = -k a`. Low-energy cross-section `\sigma \to 4\pi a^{2}` (four times geometric cross-section).

**Proper-time:** Phase shift depends on boundary condition `\psi(r = a) = 0` and free-wave kinematics outside; both unchanged. Result identical.

**Verdict:** ✅. **Companion:** [`GriffithsCh10_P10_10.wl`](../../Mathematica_Notebooks/Quantum_Mechanics/GriffithsCh10_P10_10.wl).

---

### Problem G3e-P10.18 — Relativistic-kinematics scattering

**Selection provenance:** at high energies, scattering kinematics involves the relativistic dispersion `E = \sqrt{p^{2}c^{2} + m^{2}c^{4}}`. The proper-time `K = p²/(2m) + mc²` does NOT reproduce this directly; the dual Dirac equation is required for relativistic scattering. *Substantive AI: notes a regime where the campaign's K-only treatment is insufficient.*

**Source:** Griffiths 3e Problem 10.18 (relativistic-kinematics scattering).

**Statement:** Compute the Compton-scattering kinematics in the photon energy range where `\hbar\omega \sim m_{e}c^{2}`.

**Textbook:** Compton formula: `\Delta\lambda = (\hbar/m_{e}c)(1 - \cos\theta)`. Derived from energy-momentum conservation with relativistic kinematics.

**Proper-time:** Compton scattering's kinematics requires the full relativistic dispersion `E = \sqrt{p^{2}c^{2} + m^{2}c^{4}}` for the electron's final-state energy — which is `H_{0}`, not `K`. The dual Dirac equation provides this; classical kinematics is unchanged. The dual-Dirac Compton calculation reproduces the textbook formula.

**Verdict:** ✅ at non-rel limit; ⚠ at relativistic kinematic precision, the framework's predictions depend on the dual Dirac equation (which DRQM I §II argues reproduces textbook Compton). Not the same as the `K`-only treatment of earlier Griffiths chapters.

**Notes for author review:** Compton scattering is the first Griffiths problem in the campaign where the campaign's `K`-centric framework is insufficient and the full dual Dirac equation must be invoked. The transition from `K`-treatable (non-relativistic) to dual-Dirac-treatable (relativistic) is a structural distinction worth recording.

**Companion:** [`GriffithsCh10_P10_18.wl`](../../Mathematica_Notebooks/Quantum_Mechanics/GriffithsCh10_P10_18.wl).
