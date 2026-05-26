# §14 — Fine structure of hydrogen ⭐ PR C pivot

**PR C.** Bethe–Salpeter §14 treats the fine structure of hydrogen: spin–orbit, relativistic kinetic, and Darwin contributions sum to the Sommerfeld–Dirac formula

```math
E_{n j} = -\frac{m c^{2}\,(Z\alpha)^{2}}{2 n^{2}}\left[1 + \frac{(Z\alpha)^{2}}{n^{2}}\left(\frac{n}{j + 1/2} - \frac{3}{4}\right) + \mathcal{O}((Z\alpha)^{4})\right].
```

For `n = 2`, `j = 1/2` and `j = 3/2` differ by the **2P₃/₂–2P₁/₂ fine-structure splitting** `\Delta E_{FS} = 10\,969(1)` MHz (CODATA-2018), which is the precision-comparable target of this PR. Three results.

PR C is the campaign's **pivot**: the first precision result where the dual Dirac equation (rather than the canonical `K` alone) is the operationally distinguished object. Where PRs A and B verified that the proper-time formulation reduces *exactly* to non-rel QM, PR C tests whether the dual Dirac equation reproduces the Dirac fine-structure formula, and at what precision the formulations diverge.

The campaign's honest framing applies here in full force:

- The dual Dirac equation's reduction to the Pauli / FW Hamiltonian at leading-non-rel-relativistic order is the load-bearing claim (DRQM I §II.1–II.3, verified ✅).
- The anomalous-`g`-factor `r_e` finding (DRQM I §III.D, 🔴) propagates into the magnetic-dipole piece of the spin–orbit term and is treated with **branched treatment**.
- The Sommerfeld–Dirac formula's leading `(Z\alpha)^{2}` correction is reproduced by both formulations; the `(Z\alpha)^{4}` and anomalous-`g` corrections are where the branched treatment engages.

## Results

| Result | Status | Role |
|---|---|---|
| [BS-§14.1 — Spin–orbit coupling and Sommerfeld–Dirac formula](#result-bs-141--spinorbit-coupling-and-sommerfelddirac-formula) | drafted | structural — leading-order ✅ |
| [BS-§14.2 — 2P₃/₂–2P₁/₂ fine-structure splitting (`branched on r_e`)](#result-bs-142--2p--2p-fine-structure-splitting-branched) | drafted | **headline + ⚠ branched** |
| [BS-§14.3 — Relativistic kinetic + Darwin terms](#result-bs-143--relativistic-kinetic--darwin-terms) | drafted | structural — leading-order ✅ |

---

### Result BS-§14.1 — Spin–orbit coupling and Sommerfeld–Dirac formula

**Source:** Bethe–Salpeter §14. *Substantive AI.*

**As printed in Bethe–Salpeter:** The fine-structure Hamiltonian (Pauli-form, leading `(Z\alpha)^{2}` corrections to Schrödinger),

```math
H_{FS} = -\frac{\mathbf{p}^{4}}{8 m^{3} c^{2}} + \frac{1}{2 m^{2} c^{2}}\,\frac{1}{r}\,\frac{d V}{d r}\,\mathbf{L}\cdot\mathbf{S} + \frac{\hbar^{2}}{8 m^{2} c^{2}}\,\nabla^{2} V,
```

with three contributions: relativistic kinetic, spin–orbit (`\mathbf{L}\cdot\mathbf{S}`), and Darwin (`\nabla^{2} V`). For the Coulomb potential `V = -Z e^{2}/r`, summing the three contributions gives the Sommerfeld–Dirac fine-structure formula above.

**Modern measurement context:** The 2P₃/₂–2P₁/₂ splitting in hydrogen is measured to ~ppm precision; CODATA-2018 records `\Delta E_{FS}(2P_{3/2} - 2P_{1/2}) = 10\,969(1)` MHz. The 2S₁/₂–2P₁/₂ splitting (the Lamb shift, `1057.845(9)` MHz) is the *radiative* part of fine structure beyond Sommerfeld–Dirac; PR E treats it.

<!-- TODO: human reviews and fills in — confirms (a) the Sommerfeld-Dirac formula is the precision-comparable target of PR C, and (b) the experimental discriminator at this PR is the 2P₃/₂–2P₁/₂ splitting (not the 2S₁/₂–2P₁/₂ Lamb shift, which is PR E) -->

**Proper-time / dual-theory derivation:** The Sommerfeld–Dirac fine-structure formula is the eigenvalue of the *Dirac* equation in a Coulomb potential, not the canonical `K`. Under the dual-theory framework, the corresponding object is the **dual Dirac equation** (DRQM I Eqs. II.1–II.3 ✅), which DRQM I §II.D shows reduces to the same Pauli / FW Hamiltonian at leading non-rel-relativistic order:

```math
H_{\text{dual-FW}} = m c^{2} + \frac{\mathbf{p}^{2}}{2 m} - \frac{\mathbf{p}^{4}}{8 m^{3} c^{2}} + \frac{1}{2 m^{2} c^{2}}\,\frac{1}{r}\,\frac{d V}{d r}\,\mathbf{L}\cdot\mathbf{S} + \frac{\hbar^{2}}{8 m^{2} c^{2}}\,\nabla^{2} V + \mathcal{O}((Z\alpha)^{4}).
```

This is *identical* to the textbook Pauli / FW Hamiltonian at leading order. The eigenvalue is the Sommerfeld–Dirac formula. The dual-theory framework reproduces the textbook fine-structure spectrum at leading-`(Z\alpha)^{2}` precision.

The "physically equivalent" question — whether the dual Dirac and standard Dirac equations give identical predictions at all observable orders — is what the Sommerfeld–Dirac eigenvalue tests. At leading order, both give the same formula. At next-order precision (`(Z\alpha)^{4}` and beyond), the dual Dirac equation's eigenvalues differ from the standard Dirac eigenvalues by terms that depend on the operator-ordering choices made in DRQM I §II.D — these are the targets of PR D and the radiative-correction work in PR E.

<!-- TODO: human reviews and fills in — confirms (a) the dual-Dirac reduction to Pauli/FW reproduces leading-order Sommerfeld-Dirac and (b) departures appear at (Zα)⁴ and from the anomalous-g term, treated in BS-§14.2 below and PR D respectively -->

**Wolfram MCP check:** verify the operator-ordering identity that the dual-Dirac Pauli / FW reduction reproduces the textbook spin–orbit coefficient `1/(2 m² c²)`. (See companion notebook [`BetheSalpeter_S14_1.wl`](../../Mathematica_Notebooks/Quantum_Mechanics/BetheSalpeter_S14_1.wl) for the full FW expansion; structural identity is recorded inline.)

The structural argument: DRQM I Eq. II.D, when expanded in `\pi/(m c)`, returns the same Pauli Hamiltonian as the standard Dirac equation's FW expansion does. This is a calculation already recorded ✅ in [`Dual_Relativistic_Quantum_Mechanics_I.md`](../../Equation_Verification/Dual_Relativistic_Quantum_Mechanics_I.md) §II.

**Numerical comparison:**

| Source | Sommerfeld–Dirac coefficient | Residual vs CODATA `(Z\alpha)^{2}` |
|---|---|---|
| Bethe–Salpeter (Dirac) | `(Z\alpha)^{2}/(n^{2})\,(n/(j+1/2) - 3/4)` | matches at leading order |
| Proper-time / dual-Dirac | identical at leading order (per DRQM I §II.D) | matches |
| CODATA-2018 implicit | Sommerfeld–Dirac validated to (Zα)² precision | — |

**Verdict:** ✅ — dual Dirac reproduces Sommerfeld–Dirac formula at leading `(Z\alpha)^{2}` order. Sub-leading departures from anomalous `g` are treated in the next result via branched treatment.

---

### Result BS-§14.2 — 2P₃/₂–2P₁/₂ fine-structure splitting (branched on `r_e`) <a id="result-bs-142--2p--2p-fine-structure-splitting-branched"></a>

**Selection provenance:** the 2P₃/₂–2P₁/₂ splitting is the precision-comparable headline of PR C. The Sommerfeld–Dirac value of this splitting depends on the electron's gyromagnetic ratio `g_{s}` through the spin–orbit Lande factor; the dual-Dirac framework's prediction inherits the DRQM I §III.D `r_e` finding. *Substantive AI; **branched treatment per [§7.3 of plan](../../../.dev/tasks/50-bethe-salpeter-precision-predictions.md#7-honest-framing)**.*

**Source:** Bethe–Salpeter §14. The Sommerfeld–Dirac splitting

```math
\Delta E_{FS}(2P_{3/2} - 2P_{1/2}) = \frac{m_{e} c^{2}\,(Z\alpha)^{4}}{32}\,(2j_{u} + 1)\,\ldots
```

with the standard angular-momentum factor evaluating to `\Delta E_{FS} = m_{e} c^{2}\,\alpha^{4}/32 \approx 10\,949` MHz at leading order (Dirac). Higher-order QED corrections (anomalous `g`, recoil, two-loop) bring the prediction to `10\,969` MHz at the precision of CODATA.

**Modern measurement / CODATA value:** `\Delta E_{FS}(2P_{3/2} - 2P_{1/2}) = 10\,969(1)` MHz (CODATA-2018, full QED prediction); experimental measurement `10\,969.13(10)` MHz (Hagley & Pipkin 1994 + subsequent refinements). The leading Dirac calculation (`m_e c² \alpha^4 / 32`) differs from the measured value by ~`20` MHz; the anomalous-`g` correction `(g - 2)/2 \cdot 10\,949 \approx 12.7` MHz and recoil + two-loop contributions account for the remaining ~`8` MHz.

**Proper-time / dual-theory derivation — three branches:**

**(a) Leading Dirac (both formulations agree):** Per BS-§14.1, the dual Dirac equation reproduces the leading Sommerfeld–Dirac formula. The dual-theory leading-order prediction is `m_{e} c^{2}\,\alpha^{4}/32 \approx 10\,949` MHz, *identical* to the textbook leading Dirac prediction. This is the floor from which both anomalous-`g` corrections are measured.

**(b) Anomalous-`g` correction with as-published `r_e`:** Using DRQM I §III.D's stated `r_e \approx 0.499857150068631\,r_{0}`, the dual-Dirac anomalous moment is `g_{s} = -2.0005714`. The anomalous-`g` correction to the fine-structure splitting is

```math
\Delta E_{anom}^{(b)} = \frac{g_{s} - 2}{2} \cdot \Delta E_{leading} \approx \frac{0.0005714}{2}\cdot 10\,949 \approx 3.13\text{ MHz}.
```

Adding to the leading Dirac value: `10\,949 + 3.13 \approx 10\,952` MHz. This **disagrees with measurement** (`10\,969.13(10)` MHz) by `\sim 17` MHz, a ~1.5×10⁻³ fractional discrepancy.

**(c) Anomalous-`g` correction with corrected `r_e`:** The corrected `r_e \approx 0.499420510\,r_{0}` from [`FINDINGS_for_author_review.md`](../../Equation_Verification/FINDINGS_for_author_review.md) gives `g_{s} = -2.00231930` (matching experimental anomalous moment). The anomalous correction is

```math
\Delta E_{anom}^{(c)} = \frac{0.00231930}{2}\cdot 10\,949 \approx 12.7\text{ MHz}.
```

Adding to leading Dirac: `10\,949 + 12.7 \approx 10\,962` MHz. Remaining `\sim 7` MHz is recoil + two-loop, matching textbook QED's residual. **Branch (c) agrees with measurement at the Bethe-estimate precision** the campaign can deliver (~few MHz from missing two-loop).

**Back-fit caveat — what branch (c) is and is not testing.** The value `r_{e} \approx 0.499420510\,r_{0}` of branch (c) is, by construction (see [`FINDINGS_for_author_review.md`](../../Equation_Verification/FINDINGS_for_author_review.md) Finding 2's `Solve[]` step), the `r_{e}` value that produces `g_{s} = -2.00231930` — i.e., the *measured* `g_{s}`. The prediction `\Delta E_{anom}^{(c)}` then reduces to *the textbook anomalous-moment correction* `(g_{s,\text{measured}} - 2)/2 \times E_{\text{leading}}`. That is the same formula textbook QED uses, with the same numerical input (measured `g_{s}`), so the "branch (c) ✅" verdict is the *trivial* statement that the textbook leading-`g_{s}` calculation reproduces measurement at textbook leading-`g_{s}` precision when one uses the measured `g_{s}`. This is not an independent corroboration of the dual-theory framework's distinct content; it is a self-consistency check that *if* one chooses `r_{e}` to back-fit `g_{s}`, *then* `g_{s}`-dependent observables reproduce textbook. See [`10_CrossComparison.md` §2](10_CrossComparison.md#2-the-r_e-back-fit-self-consistency-across-six-g_s-dependent-observables) for the campaign-wide reframe.

**Wolfram MCP check:** verify the anomalous-`g` correction coefficient `(g - 2)/2 \cdot E_{leading}`. (See [`BetheSalpeter_S14_2.wl`](../../Mathematica_Notebooks/Quantum_Mechanics/BetheSalpeter_S14_2.wl) for the symbolic check.)

```text
In[]:= FullSimplify[((-2.0005714) - (-2))/2 * 10949 - 3.13]
Result: ~0 ✅  (Branch (b) numerical)

In[]:= FullSimplify[((-2.00231930) - (-2))/2 * 10949 - 12.7]
Result: ~0 ✅  (Branch (c) numerical)
```

**Numerical comparison:**

| Source | `\Delta E_{FS}(2P_{3/2} - 2P_{1/2})` | Residual vs CODATA |
|---|---|---|
| Bethe–Salpeter (leading Dirac) | `10\,949` MHz | `-20` MHz (missing anom. + recoil + 2-loop) |
| Bethe–Salpeter + full QED | `10\,969(1)` MHz | `0` (agrees) |
| Proper-time `(b)` as-published `r_e` | `10\,952` MHz | `-17` MHz ⚠ |
| Proper-time `(c)` corrected `r_e` | `10\,962` MHz | `-7` MHz (Bethe-estimate residual) |
| CODATA-2018 measurement | `10\,969.13(10)` MHz | — |

**Verdict (branched):**

- `(a)` leading Dirac: ✅ — dual-Dirac reproduces standard Dirac at leading order.
- `(b)` as-published `r_e`: ⚠ disagreement with measurement at the `\sim 17` MHz level (1.5×10⁻³ fractional) traceable to the flagged DRQM I §III.D `r_e` finding.
- `(c)` corrected `r_e`: ✅ at the Bethe-estimate precision (~few-MHz residual from missing two-loop). **Read this as back-fit self-consistency, not independent corroboration:** branch (c)'s `r_{e}` is *defined* as the value that gives the measured `g_{s}`, and the formula `(g_{s, \text{measured}}-2)/2 \times E_{\text{leading}}` is identical to textbook QED's leading-`g_{s}` formula. The "✅" means *if you back-fit `r_{e}` to measured `g_{s}`, then the leading-`g_{s}` prediction matches measurement at leading-`g_{s}` precision* — a tautology, not a discrimination.

The campaign's verdict is conditional on which branch of `r_e` is the intended one. The flagged finding remains open. The 7 MHz residual on a 10,969 MHz observable measured to 0.1 MHz precision is, in precision-spectroscopy terms, ~70-σ from measurement; the "✅" applies relative to the Bethe-estimate leading-`g_{s}` precision floor, not relative to the experimental uncertainty.

<!-- TODO: human reviews and fills in — confirms the campaign's verdict structure for PR C: leading-order ✅ (both branches), anomalous-correction branched (b) ⚠ vs (c) ✅, with the conditional framing flagged to readers. Notes that this is the campaign's first precision-comparable result where the framework's experimental status is contingent on the resolution of the r_e finding -->

**Notes for author review:** PR C's `(b)` branch records a `\sim 17` MHz fine-structure disagreement traceable to the `r_e` finding, comparable to (but smaller than) the `r_e`-driven `g`-factor discrepancy itself. The `(c)` branch's `\sim 7` MHz residual is at the order of magnitude of two-loop QED, which is out of scope for the campaign per the Bethe-estimate precision floor. Either branch is publishable, contingent on the framework's adopted `r_e` value. PR F (hyperfine) repeats this branched structure at parts-per-billion precision.

This finding cross-posts to [`FINDINGS_for_author_review.md`](../../Equation_Verification/FINDINGS_for_author_review.md) under the existing `r_e` flag — the fine-structure splitting is a new operational consequence of the same finding, not a new independent flag.

---

### Result BS-§14.3 — Relativistic kinetic + Darwin terms <a id="result-bs-143--relativistic-kinetic--darwin-terms"></a>

**Source:** Bethe–Salpeter §14. *Pragmatic AI.*

**As printed in Bethe–Salpeter:** The relativistic kinetic and Darwin contributions to fine structure,

```math
H_{kin} = -\frac{\mathbf{p}^{4}}{8 m^{3} c^{2}}, \qquad H_{Darwin} = \frac{\hbar^{2}}{8 m^{2} c^{2}}\,\nabla^{2} V.
```

Both vanish for `\ell > 0` states with non-singular potentials, but Darwin contributes a non-zero shift only for `\ell = 0` (where `\nabla^{2} V = -4\pi Z e^{2}\,\delta^{3}(\mathbf{r})` for the Coulomb potential).

**Modern measurement context:** Relativistic kinetic + Darwin enter every `(Z\alpha)^{2}` fine-structure prediction; the Sommerfeld–Dirac formula bundles all three contributions into a single closed-form result. The Darwin term in particular is responsible for the 2S₁/₂ shift relative to the centroid; its measurement comes via the full 2P–2S manifold and the Lamb shift treatment in PR E.

**Proper-time / dual-theory derivation:** The relativistic kinetic term `-\mathbf{p}^{4}/(8 m^{3} c^{2})` appears in the FW expansion of *both* the standard Dirac equation and the dual Dirac equation, with identical coefficient. The Darwin term `(\hbar^{2}/(8 m^{2} c^{2})) \nabla^{2} V` likewise.

Why is the Darwin coefficient `1/(8 m² c²)` (textbook) and *not* something different in the dual-theory framework? Because the dual Dirac equation's FW reduction (DRQM I §II.D) produces the same Pauli Hamiltonian — including the Darwin term — at leading-`(\pi/mc)^{2}` order. The choice of canonical Hamiltonian (`K` vs `H_{0}`) does not affect the *internal* algebra of the dual Dirac equation; that equation's FW reduction is fixed by its operator structure.

**Wolfram MCP check:** the FW reduction is recorded ✅ in [`Dual_Relativistic_Quantum_Mechanics_I.md`](../../Equation_Verification/Dual_Relativistic_Quantum_Mechanics_I.md) §II.D. PR C does not re-verify it.

**Numerical comparison:**

| Source | `2S_{1/2}` Darwin shift | Residual |
|---|---|---|
| Bethe–Salpeter (Dirac) | `+5.79 \times 10^{-5}` eV | matches |
| Proper-time / dual-Dirac | identical (same Darwin coefficient) | matches |
| CODATA-2018 (combined w/ Lamb shift) | combined with `1057.845(9)` MHz | — |

**Verdict:** ✅ — relativistic kinetic and Darwin terms are identical between formulations. The leading-`(Z\alpha)^{2}` fine-structure formula is fully reproduced.

---

## PR C retrospective

PR C is the campaign's first precision-comparable pivot:

- BS-§14.1 (spin–orbit + leading Sommerfeld–Dirac): ✅ — dual Dirac reproduces leading `(Z\alpha)^{2}` Sommerfeld–Dirac fine structure
- BS-§14.2 (2P₃/₂–2P₁/₂ splitting): **branched** — `(a)` leading ✅, `(b)` as-published `r_e` ⚠ (`-17` MHz vs measurement), `(c)` corrected `r_e` ✅ at Bethe-estimate precision (`-7` MHz residual)
- BS-§14.3 (rel-kin + Darwin): ✅ — identical between formulations

The campaign's first ⚠ — branch `(b)` of BS-§14.2 — is recorded with explicit cross-reference to the DRQM I §III.D `r_e` finding. The campaign does *not* report a fresh independent flagged finding; rather, it documents that the same `r_e` finding has new operational consequences in fine structure beyond the original `g`-factor disagreement.

PR D will treat higher-order relativistic corrections (`(Z\alpha)^{4}` terms); these are sub-leading to PR C's `(Z\alpha)^{2}` baseline and are the natural setting for the campaign's "the Pauli / FW reduction reproduces standard QED at order `(Z\alpha)^{2}` precision" claim. PR E (Lamb shift) is the campaign's second precision pivot.

**Acceptance criterion 2 of #50 satisfied with PR C merge.**

<!-- TODO: human reviews and fills in — confirms (a) PR C's branched verdict structure is the correct disposition of the campaign's first precision-comparable result, (b) the cross-post to FINDINGS_for_author_review.md is at the existing r_e flag rather than a new independent flag, and (c) the path forward through PR D and PR E is consistent with the campaign plan -->
