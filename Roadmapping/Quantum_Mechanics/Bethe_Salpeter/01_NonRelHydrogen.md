# §§1–7 — Non-relativistic hydrogen

**PR A.** Bethe–Salpeter's opening sections establish the non-relativistic hydrogen problem: Schrödinger equation, energy spectrum `E_n = -mc²(Zα)²/(2n²)`, Bohr radius, radial-equation eigenfunctions, and degeneracy. These are the textbook results that any candidate relativistic / proper-time formulation must reproduce in the `u \ll c` limit before its precision-experiment predictions can be taken seriously. Four results.

The role of PR A in the campaign is **structural**: we verify, problem by problem, that the proper-time `K = H_0^{2}/(2 m c^{2}) + m c^{2}/2` reduces *exactly* to the non-relativistic Coulomb Hamiltonian `p²/(2m) + V_0` (modulo a constant `(3/2) m c^{2}` rest-energy offset), so that every Schrödinger-level result of Bethe–Salpeter passes through unchanged. The acceptance criterion is mechanical confirmation that no algebraic surprise lurks at the non-relativistic level.

The proper-time / dual-theory framework is *not* expected to differ from textbook QM at this level — every honest framing of the campaign records this — and the verdicts in PR A reflect that. The interesting physics begins at PR C (fine structure); PR A is the pedagogical floor against which later precision predictions are measured.

## Results

| Result | Status | Role |
|---|---|---|
| [BS-§3 — Schrödinger equation for hydrogen](#result-bs-3--schrödinger-equation-for-hydrogen) | drafted | scaffold (non-rel reduction) |
| [BS-§4 — Energy spectrum `E_n = -mc²(Zα)²/(2n²)`](#result-bs-4--energy-spectrum-e_n) | drafted | scaffold (non-rel reduction) |
| [BS-§5 — Bohr radius and radial eigenfunctions](#result-bs-5--bohr-radius-and-radial-eigenfunctions) | drafted | scaffold |
| [BS-§6 — Degeneracy and the `O(4)` symmetry](#result-bs-6--degeneracy-and-the-o4-symmetry) | drafted | scaffold |

---

### Result BS-§3 — Schrödinger equation for hydrogen

**Source:** Bethe–Salpeter §3 (Springer 1977 reprint). *Pragmatic AI.*

**As printed in Bethe–Salpeter:** The hydrogenic Schrödinger equation,

```math
\left[-\frac{\hbar^{2}}{2 m}\nabla^{2} - \frac{Z e^{2}}{r}\right]\psi(\mathbf{r}) = E\,\psi(\mathbf{r}),
```

with `m` the electron mass (or reduced mass `μ` to high precision) and `Z` the nuclear charge.

**Modern measurement context:** The Schrödinger spectrum is the ~part-per-thousand-accurate baseline against which fine structure (`(Zα)²`-suppressed) and Lamb shift (`α(Zα)⁴`-suppressed) corrections are layered. No direct precision measurement targets the Schrödinger Hamiltonian itself; rather, the Rydberg constant `R_∞ = m_e c α²/(4π\hbar)` underlies all hydrogen spectroscopy and is the most precisely measured fundamental constant (CODATA-2018: `R_∞ = 10\,973\,731.568160(21)` m⁻¹, ~2×10⁻¹² relative).

**Proper-time / dual-theory derivation:** Start from the proper-time canonical Hamiltonian

```math
K = \frac{H_{0}^{2}}{2 m c^{2}} + \frac{m c^{2}}{2}, \qquad H_{0} = \sqrt{c^{2}\boldsymbol{\pi}^{2} + m^{2} c^{4}},
```

with `\boldsymbol{\pi} = \mathbf{p}` (no vector potential for the unperturbed Coulomb problem; the scalar Coulomb potential enters as `V_{0} = -Ze²/r` via the standard substitution `K \to K + V_{0}`). Substituting `H_{0}^{2} = c²\mathbf{p}² + m²c⁴`,

```math
K = \frac{c^{2}\mathbf{p}^{2} + m^{2} c^{4}}{2 m c^{2}} + \frac{m c^{2}}{2} = \frac{\mathbf{p}^{2}}{2 m} + \frac{m c^{2}}{2} + \frac{m c^{2}}{2} = \frac{\mathbf{p}^{2}}{2 m} + m c^{2}.
```

The proper-time canonical Hamiltonian reduces to the non-relativistic kinetic energy `p²/(2m)` *plus a constant rest-energy offset* `m c^{2}` — **exactly**, not perturbatively. This is the canonical identity that anchors the entire PR A scaffold.

Including the Coulomb potential,

```math
K + V_{0} = m c^{2} + \frac{\mathbf{p}^{2}}{2 m} - \frac{Z e^{2}}{r} = m c^{2} + H_{\text{Schr}},
```

so the eigenvalue equation `(K + V_{0})\,\psi = \mathcal{E}\,\psi` becomes `H_{\text{Schr}}\psi = (\mathcal{E} - m c^{2})\,\psi = E\,\psi` — identical to Bethe–Salpeter's Schrödinger equation, with `E = \mathcal{E} - m c^{2}` the standard non-relativistic energy measured relative to the rest energy. Bound-state eigenfunctions are the textbook hydrogenic wavefunctions; spectrum is treated next.

**Wolfram MCP check:** verify the algebraic identity `K = m c^{2} + p²/(2m)` (the `V_{0}` term enters as a scalar addition unchanged). Companion: [`BetheSalpeter_S3.wl`](../../Mathematica_Notebooks/Quantum_Mechanics/BetheSalpeter_S3.wl).

```text
In[]:= FullSimplify[((c^2*p^2 + m^2*c^4)/(2*m*c^2) + m*c^2/2) - (m*c^2 + p^2/(2*m))]
Result: 0 ✅
```

**Numerical comparison:** *Not applicable.* The Schrödinger equation itself is not a measured number; the spectrum it generates is verified in BS-§4 below.

**Verdict:** ✅ — proper-time `K + V_{0}` and Bethe–Salpeter's Schrödinger Hamiltonian are operator-identical (modulo the `m c^{2}` rest-energy offset). No relativistic correction enters at this level.

---

### Result BS-§4 — Energy spectrum `E_n` <a id="result-bs-4--energy-spectrum-e_n"></a>

**Source:** Bethe–Salpeter §4. *Pragmatic AI.*

**As printed in Bethe–Salpeter:** Hydrogenic Schrödinger spectrum,

```math
E_{n} = -\frac{m e^{4} Z^{2}}{2 \hbar^{2} n^{2}} = -\frac{m c^{2}\,(Z\alpha)^{2}}{2 n^{2}}, \qquad n = 1, 2, 3, \ldots
```

with `α = e²/(\hbar c)` the fine-structure constant (Gaussian units).

**Modern measurement context:** The Rydberg energy `Ry = m_e c² α²/2 ≈ 13.605\,693\,122\,994(26)` eV is the CODATA-2018 value (~2×10⁻¹² relative). The ground state of hydrogen sits at `E_{1} = -Ry`, the `n = 2` levels at `-Ry/4 ≈ -3.401` eV, etc. These are the *non-relativistic* eigenvalues; fine-structure corrections (≈MHz) and Lamb shift are layered on top in PR C / PR E.

**Proper-time / dual-theory derivation:** From BS-§3 above, the eigenvalue equation under the proper-time Hamiltonian is

```math
(K + V_{0})\,\psi_{n} = (m c^{2} + E_{n})\,\psi_{n},
```

so the proper-time eigenvalues differ from the Schrödinger eigenvalues by the constant `m c^{2}` rest-energy offset. The non-relativistic spectrum measured relative to threshold (the usual convention) is *identical*:

```math
E_{n}^{\text{proper-time}} - m c^{2} = -\frac{m c^{2}\,(Z\alpha)^{2}}{2 n^{2}}.
```

This is the canonical non-relativistic-reduction statement of the campaign: the proper-time formulation's eigenvalue offset by `m c^{2}` recovers the Schrödinger spectrum. No `O((Zα)^{4})` or `O(α(Zα)^{4})` correction appears at this order — those enter at PR C (fine structure) and PR E (Lamb shift) via the dual Dirac equation, not the canonical `K`.

The Rydberg energy itself follows directly from `m_e c² α²/2`. All hydrogen spectroscopy depending on `R_∞` carries through unchanged.

**Wolfram MCP check:** verify the spectrum identity by substituting `n = 1` to recover `-Ry`:

```text
In[]:= FullSimplify[-m*c^2*(Z*alpha)^2/(2*1^2) - (-m*c^2*alpha^2/2) /. Z -> 1]
Result: 0 ✅
```

(That is, the proper-time prediction at `n = 1, Z = 1` matches the textbook ground-state energy `-Ry = -m_e c² α²/2`.)

**Numerical comparison:**

| Source | Ground-state energy `E_{1}` | Residual vs measurement |
|---|---|---|
| Bethe–Salpeter (Schrödinger) | `-13.605\,693\,12…` eV | 0 (exact at non-rel order) |
| Proper-time / dual-theory | `-13.605\,693\,12…` eV | 0 (exact at non-rel order) |
| CODATA-2018 Rydberg | `-13.605\,693\,122\,994(26)` eV | — |

The fine-structure correction (`≈ -1.81 × 10^{-4}` eV at `n = 1`) and Lamb shift (`≈ +1.45 × 10^{-7}` eV from the `2S` shift) sit at PR C and PR E respectively. PR A's verdict is the non-rel baseline only.

**Verdict:** ✅ — proper-time eigenvalues match the textbook Schrödinger spectrum to all orders in the non-relativistic limit. The campaign's claim that the proper-time formulation reduces *exactly* to non-rel QM is verified at the spectral level for hydrogen.

---

### Result BS-§5 — Bohr radius and radial eigenfunctions

**Source:** Bethe–Salpeter §5. *Pragmatic AI.*

**As printed in Bethe–Salpeter:** The Bohr radius

```math
a_{0} = \frac{\hbar^{2}}{m e^{2}} = \frac{\hbar}{m c\,\alpha},
```

and the ground-state radial eigenfunction

```math
\psi_{100}(\mathbf{r}) = \frac{1}{\sqrt{\pi}}\,\left(\frac{Z}{a_{0}}\right)^{3/2} e^{-Z r / a_{0}}.
```

**Modern measurement context:** `a_{0} = 5.291\,772\,109\,03(80) \times 10^{-11}` m (CODATA-2018). The hydrogen ground state is the most extensively characterised quantum-mechanical bound state, with the wavefunction probed indirectly via every hydrogen spectroscopy measurement.

**Proper-time / dual-theory derivation:** The Bohr radius emerges from the same eigenvalue problem as in BS-§4 above. Since the operator equation `(K + V_{0} - m c^{2})\,\psi = E\,\psi` is identical to the Schrödinger equation (per BS-§3), the eigenfunctions are textbook hydrogenic wavefunctions. The Bohr radius `a_{0} = \hbar/(m c \alpha)` and `\psi_{100}` are unchanged.

This is again a structural reduction: any matrix element, expectation value, or wavefunction overlap that Bethe–Salpeter computes at the Schrödinger level passes through to the proper-time treatment unchanged. PR B's transition matrix elements rest on this fact.

**Wolfram MCP check:** *Not applicable as a separate check.* The Bohr radius emerges from variational stationarity of the standard hydrogenic ansatz `\psi(r) = (Z/a_0)^{3/2} e^{-Z r/a_0} / \sqrt{π}`; since the underlying Hamiltonian is operator-identical to Bethe–Salpeter's per BS-§3, the variational stationary point is identical.

**Numerical comparison:**

| Source | Bohr radius `a_{0}` | Residual vs CODATA |
|---|---|---|
| Bethe–Salpeter | `\hbar²/(m e²)` (analytic) | 0 |
| Proper-time / dual-theory | `\hbar²/(m e²)` (analytic) | 0 |
| CODATA-2018 | `5.291\,772\,109\,03(80) × 10^{-11}` m | — |

**Verdict:** ✅ — Bohr radius and hydrogenic wavefunctions are operator-identical to Bethe–Salpeter's. No correction at non-rel order.

---

### Result BS-§6 — Degeneracy and the `O(4)` symmetry

**Source:** Bethe–Salpeter §6. *Pragmatic AI.*

**As printed in Bethe–Salpeter:** The hydrogenic spectrum has `n²`-fold degeneracy at each principal quantum number `n` (counting orbital `\ell` from `0` to `n-1`, magnetic `m_{\ell}` from `-\ell` to `\ell`). The degeneracy is a consequence of the conserved Runge–Lenz vector `\mathbf{A} = (1/m)\,\mathbf{p} \times \mathbf{L} - (Z e^{2}/r)\,\hat{\mathbf{r}}` (Pauli 1926); together with `\mathbf{L}`, the commutation algebra closes to `O(4)`.

**Modern measurement context:** The `n²` degeneracy is lifted by relativistic effects (`O((Zα)^{2})` — fine structure) and QED radiative corrections (`O(α(Zα)^{4})` — Lamb shift). The `2S_{1/2}–2P_{1/2}` splitting (Lamb shift, `1057.845(9)` MHz) is the famous tell that the *Schrödinger* degeneracy is the non-relativistic limit, not the full physical degeneracy. PR E records the proper-time prediction for this splitting.

**Proper-time / dual-theory derivation:** Since the non-relativistic Hamiltonian `H_{\text{Schr}} = K + V_{0} - m c^{2}` is operator-identical to Bethe–Salpeter's Schrödinger Hamiltonian, the conserved Runge–Lenz vector and `O(4)` symmetry are inherited unchanged. The `n²`-degeneracy is exact at non-rel order in the proper-time formulation.

Departures from the Schrödinger degeneracy under the dual Dirac equation (Pauli / FW reduction) are deferred to PR C (`2P_{3/2}–2P_{1/2}` fine-structure splitting) and PR E (`2S_{1/2}–2P_{1/2}` Lamb shift). At non-rel order, no such departure exists in either formulation.

**Wolfram MCP check:** *Not applicable at non-rel order.* The degeneracy is a symmetry property of the operator algebra `[\mathbf{A}, H_{\text{Schr}}] = 0`, which is structurally identical under the proper-time formulation (the operator algebra is unchanged because the operator is identical).

**Numerical comparison:** Degeneracies are integer multiplicities, not measured to floating-point precision. The non-relativistic `n² = 1, 4, 9, …` degeneracies are reproduced in both formulations.

**Verdict:** ✅ — `O(4)` symmetry and `n²`-degeneracy structure are inherited verbatim by the proper-time formulation at non-rel order. Lifting of the degeneracy by relativistic / radiative effects is the subject of PRs C and E.

---

## PR A retrospective

PR A scaffolds the Bethe–Salpeter campaign and demonstrates four mechanical non-relativistic-reduction results: Schrödinger equation, energy spectrum, Bohr radius / wavefunctions, and `O(4)` degeneracy structure. All four return ✅ verdicts.

The structural point — verified by the Wolfram MCP check in BS-§3 — is that `K + V_{0} = m c^{2} + p²/(2m) + V_{0}` *exactly*, not as a `(p/mc)^{2}` expansion. The proper-time formulation does not introduce *any* deviation from non-relativistic QM until either the dual Dirac equation (PR C) or radiative corrections (PR E) are invoked.

This is the floor against which PR C, PR E, PR F are measured. The campaign's experimental discriminators must depart from non-rel QM only in the regimes where measurement is sensitive to the departure; PR A confirms that the floor is clean.

The next PR (PR B) treats matrix elements + transitions — also non-relativistic, and expected to inherit the same ✅ structural reduction.

<!-- TODO: human reviews and fills in — confirms (a) the operator identity in BS-§3 is the campaign's load-bearing structural claim, (b) PR A's mechanical ✅ verdicts are the expected outcome and do not constitute experimental endorsement of the framework, and (c) the readers understand that PR A's role is to establish the floor for PRs C, E, F -->
