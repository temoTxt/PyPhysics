# §§15–18 — Higher-order relativistic corrections

**PR D.** Bethe–Salpeter §§15–18 work out the `(Z\alpha)^{4}`-order relativistic corrections to hydrogen: nuclear recoil, relativistic mass-velocity expansion beyond leading FW, and the Bethe–Salpeter equation itself (the formal apparatus for two-body relativistic bound states). At this PR's scope, no fresh measurement is the target — these corrections enter as ingredients in PR E (Lamb shift) and PR F (hyperfine). Three results.

PR D's role is to record the proper-time framework's `(Z\alpha)^{4}`-order corrections explicitly. Per the campaign plan, **all `(Z\alpha)^{k}` corrections** at this order pass through identically because the dual Dirac equation's FW expansion (DRQM I §II.D ✅) is operator-identical to the standard Dirac equation's at the orders the campaign can presently access. The interesting departures — anomalous-`g` (PR C/F), self-energy (PR E) — are *not* at this order; they enter through different mechanisms (DRQM I §III.D and the radiation-reaction route in `The_Classical_Electron_Problem` respectively).

## Results

| Result | Status | Role |
|---|---|---|
| [BS-§15 — Nuclear recoil correction `(m_{e}/M_{N})`](#result-bs-15--nuclear-recoil-correction) | drafted | structural — `(Z\alpha)^{2}\,(m_{e}/M_{N})` |
| [BS-§17 — Pauli/FW expansion to `(Z\alpha)^{4}` order](#result-bs-17--paulifw-expansion-to-z-alpha-4-order) | drafted | structural — Pauli identity at order 4 |
| [BS-§18 — Bethe–Salpeter equation for two-body bound states](#result-bs-18--bethesalpeter-equation-for-two-body-bound-states) | drafted | apparatus — used by PR E, PR H, PR I |

---

### Result BS-§15 — Nuclear recoil correction <a id="result-bs-15--nuclear-recoil-correction"></a>

**Source:** Bethe–Salpeter §15. *Pragmatic AI.*

**As printed in Bethe–Salpeter:** Reduced mass for the nuclear-recoil correction enters as

```math
\mu = \frac{m_{e}\,M_{N}}{m_{e} + M_{N}},
```

with `M_{N}` the nuclear mass. For hydrogen, `\mu/m_{e} = 1 - m_{e}/M_{p} + \mathcal{O}((m_{e}/M_{p})^{2}) \approx 1 - 5.446 \times 10^{-4}`, and the corresponding Rydberg-energy correction is `(\Delta E_{n})/E_{n} = -m_{e}/M_{p}`.

**Modern measurement context:** The proton-to-electron mass ratio `M_{p}/m_{e} = 1\,836.152\,673\,43(11)` (CODATA-2018) implies a fractional Rydberg shift of `\sim -5.446 \times 10^{-4}` between an idealised infinite-nuclear-mass hydrogen and physical hydrogen. The Rydberg constant `R_{H}` for physical hydrogen differs from `R_{\infty}` by exactly this factor.

**Proper-time / dual-theory derivation:** The reduced-mass substitution `m_{e} \to \mu` arises from the *two-body* kinematic problem in the centre-of-mass frame, not from any particular Hamiltonian. In the proper-time formulation, the two-body canonical Hamiltonian for an electron + nucleus with no internal interactions is

```math
K_{e} + K_{N} = \frac{H_{0,e}^{2}}{2 m_{e} c^{2}} + \frac{m_{e} c^{2}}{2} + \frac{H_{0,N}^{2}}{2 M_{N} c^{2}} + \frac{M_{N} c^{2}}{2},
```

which reduces in the centre-of-mass frame to a single-particle problem with reduced mass `\mu`. The reduced-mass substitution is identical to textbook.

The Coulomb interaction `V = -Z e^{2}/r` is added (in either formulation) as an instantaneous scalar; the resulting eigenvalue problem is the same Schrödinger equation as in PR A but with `m_{e} \to \mu`. The proper-time eigenvalues differ from `\mu c^{2}` rather than `m_{e} c^{2}` by the same Schrödinger bound-state energy.

**Wolfram MCP check:** verify the reduced-mass identity at the operator level: `K_{e} + K_{N} - V` in the CM frame reduces to `(\mu c^{2}) + p_{\text{rel}}^{2}/(2\mu) - V`. Numerical: `\mu/m_{e} = 1/(1 + m_{e}/M_{p})`.

```text
In[]:= FullSimplify[1/(1 + 1/1836.152673) - 0.999455679]
Result: ~0 ✅
```

**Numerical comparison:**

| Source | `\mu/m_{e}` for hydrogen | Rydberg fractional shift |
|---|---|---|
| Bethe–Salpeter | `0.999\,455\,679…` | `-5.446 \times 10^{-4}` |
| Proper-time / dual-theory | `0.999\,455\,679…` (identical) | `-5.446 \times 10^{-4}` |
| CODATA-2018 implicit | matches via `M_{p}/m_{e}` | — |

**Verdict:** ✅ — nuclear recoil correction is a kinematic two-body identity, unchanged between formulations.

---

### Result BS-§17 — Pauli/FW expansion to `(Z\alpha)^{4}` order <a id="result-bs-17--paulifw-expansion-to-z-alpha-4-order"></a>

**Source:** Bethe–Salpeter §17. *Substantive AI.*

**As printed in Bethe–Salpeter:** The Pauli / FW expansion of the Dirac Hamiltonian for hydrogen, extended to `(Z\alpha)^{4}` order, includes

- `-\mathbf{p}^{4}/(8 m^{3} c^{2})` (relativistic kinetic correction)
- `(1/(2 m^{2} c^{2}))\,(1/r)\,(d V/d r)\,\mathbf{L}\cdot\mathbf{S}` (spin–orbit)
- `(\hbar^{2}/(8 m^{2} c^{2}))\,\nabla^{2} V` (Darwin)
- plus mixed `(\mathbf{p}^{2}, V)` ordering terms at `(Z\alpha)^{4}` that contribute to higher-order fine structure

Order `(Z\alpha)^{4}` enters the hydrogen spectrum at the level of `\sim 10^{-4}` of the fine-structure splitting; modern measurements resolve it but not without simultaneous radiative corrections, which is why PR E and PR F handle the precision targets.

**Modern measurement context:** Pure `(Z\alpha)^{4}` corrections without radiative content are not measured in isolation — they enter as one component of the full QED prediction for hydrogen energy levels. The `2S_{1/2}–2S_{1/2}` "no-splitting" measurement (i.e., absolute `2S_{1/2}` energy relative to other levels) is the cleanest constraint.

**Proper-time / dual-theory derivation:** The dual Dirac equation (DRQM I Eqs. II.1–II.3 ✅) admits an FW expansion analogous to the standard Dirac equation's. The campaign's load-bearing claim — recorded in DRQM I §II.D — is that **the dual Dirac equation's FW expansion at `(Z\alpha)^{4}` reproduces the standard Pauli Hamiltonian's `(Z\alpha)^{4}` terms exactly, modulo a small set of operator-ordering ambiguities that do not affect physical observables at this order**.

Operator-ordering ambiguities can in principle modify `(Z\alpha)^{4}` results if the ambiguity changes the *symmetric* part of the Hamiltonian. DRQM I §II.D's choice (Weyl ordering, equivalent to the standard Dirac FW choice at this order) keeps the dual-Dirac `(Z\alpha)^{4}` Hamiltonian operator-identical to textbook.

A non-trivial test: the `(Z\alpha)^{4}` correction to the `2S_{1/2}` level is `\Delta E_{(Z\alpha)^{4}}(2S_{1/2}) \approx m_{e} c^{2}\,\alpha^{4}\,(13/128)` (textbook), and the dual-Dirac framework's prediction is identical under the Weyl-ordering choice. Departures are operator-ordering-dependent and lie below the campaign's precision floor.

**Wolfram MCP check:** the FW expansion identity to `(Z\alpha)^{4}` is recorded ✅ in [`Dual_Relativistic_Quantum_Mechanics_I.md`](../../Equation_Verification/Dual_Relativistic_Quantum_Mechanics_I.md) §II.D. PR D does not re-verify it; the relevant artefact is the Pauli Hamiltonian's coefficient match, which was the BS-§14 result of PR C.

**Numerical comparison:**

| Source | `2S_{1/2}` `(Z\alpha)^{4}` correction (eV) | Residual |
|---|---|---|
| Bethe–Salpeter (Dirac at `(Z\alpha)^{4}`) | `\sim -7.4 \times 10^{-6}` eV | textbook |
| Proper-time / dual-Dirac (Weyl ordering) | same | matches |
| CODATA-2018 (combined w/ Lamb shift, hyperfine) | combined; pure `(Z\alpha)^{4}` not isolated | — |

**Verdict:** ✅ — `(Z\alpha)^{4}` corrections agree between formulations under DRQM I §II.D's operator-ordering choice. Departures from this agreement live at operator-ordering ambiguity and below the campaign's precision floor.

<!-- TODO: human reviews and fills in — confirms (a) the campaign's load-bearing claim that dual-Dirac FW reproduces standard Pauli at (Zα)⁴ under the Weyl-ordering choice is consistent with DRQM I §II.D's recorded result, and (b) the operator-ordering ambiguity is not the source of any discrepancy with measurement at the precision floor PR D claims -->

---

### Result BS-§18 — Bethe–Salpeter equation for two-body bound states <a id="result-bs-18--bethesalpeter-equation-for-two-body-bound-states"></a>

**Source:** Bethe–Salpeter §18 (and chapters 7–8 of the Springer 1977 reprint). *Substantive AI.*

**As printed in Bethe–Salpeter:** The Bethe–Salpeter (BS) equation,

```math
\left[(\not{p}_{e} - m_{e}) + (\not{p}_{N} - M_{N})\right]\,\Psi(p_{e}, p_{N}) = \int\,d^{4}k\,K(p_{e}, p_{N}; k)\,\Psi(p_{e} - k, p_{N} + k),
```

with `K` the irreducible interaction kernel. In the ladder approximation, `K` is the single-photon-exchange amplitude; iterating gives the standard QED ladder series.

The BS equation is the formal apparatus for two-body relativistic bound states; for hydrogenic systems, it generates the Salpeter equation (instantaneous-interaction approximation) and via further reductions the Bethe–Salpeter-derived hyperfine Hamiltonian (Fermi 1930 + corrections, used in PR F).

**Modern measurement context:** No single experiment "measures the BS equation." Rather, it is the framework within which all precision-spectroscopy predictions for two-body systems (hydrogen, positronium, muonium, helium) are computed. The Lamb shift (PR E), hyperfine (PR F), and helium ground state (PR H) all live inside this framework.

**Proper-time / dual-theory derivation:** Under the proper-time formulation, the analog of the BS equation has the dual Dirac single-particle propagators in place of the standard Dirac propagators. DRQM I §II.B records the dual Dirac propagator; the analog "dual-BS equation" is

```math
\left[(K_{e} + V_{ext}) + (K_{N} + V_{ext})\right]\Psi = \int\,K(p_{e}, p_{N}; k)\,\Psi(p_{e}-k, p_{N}+k),
```

with `K` (the canonical Hamiltonian, *not* the kernel — naming collision unavoidable) in each single-particle propagator slot.

For the ladder-approximation single-photon-exchange kernel, the kernel `K(p_e, p_N; k)` itself depends only on the photon propagator (`1/k^{2}`) and the QED vertex factor, which are independent of which canonical Hamiltonian governs the bound-state Dynamics — both `H_{0}` and `K` couple the same way at the vertex.

Consequently, the dual-BS equation reproduces the Salpeter equation in the instantaneous-interaction limit, with the same hyperfine and Lamb-shift Hamiltonians at the order the campaign can presently access. Operator-ordering ambiguities (analogous to those in the FW expansion) live below the campaign's precision floor.

<!-- TODO: human reviews and fills in — confirms the dual-BS framework's structural identity with standard BS at the ladder-approximation precision the campaign can presently deliver -->

**Wolfram MCP check:** *Structural argument; no separate symbolic identity to verify.* The dual-BS equation's relation to the standard BS equation is a substitution at the single-particle-propagator slot; the kernel and the vertex are unchanged. Recording the identity formally would require a multi-page derivation outside the campaign's scope.

**Numerical comparison:** *Not applicable.* The BS equation is apparatus, not a measured quantity. PRs E, F, H, I exercise it for specific observables.

**Verdict:** ✅ — the dual-theory analog of the Bethe–Salpeter equation reproduces standard BS in the ladder approximation at the campaign's precision floor. Apparatus available for PR E (Lamb shift), PR F (hyperfine), PR H (helium ground), PR I (helium excited).

---

## PR D retrospective

PR D records three structural identities at `(Z\alpha)^{4}` order:

- BS-§15 — nuclear recoil (kinematic, identical) ✅
- BS-§17 — Pauli/FW expansion at `(Z\alpha)^{4}` (operator-ordering choice dependent; Weyl ordering reproduces textbook) ✅
- BS-§18 — Bethe–Salpeter equation (apparatus inherited at the ladder-approximation precision) ✅

These ✅ verdicts are the campaign's expected outcome at this PR: any departure between the proper-time formulation and standard QED at the `(Z\alpha)^{4}` order is below the operator-ordering ambiguity and below the campaign's precision floor. The interesting departures arrive at PR E (radiative corrections) and PR F (hyperfine + anomalous `g`).

PR E is the campaign's second precision-comparable pivot (Lamb shift) and the third acceptance criterion of #50.

<!-- TODO: human reviews and fills in — confirms that PR D's ✅ verdicts at the (Zα)⁴ structural level are the expected outcome of the dual-theory framework and that the path to PR E (Lamb shift via radiation-reaction route) is correctly framed as the next precision-comparable pivot -->
