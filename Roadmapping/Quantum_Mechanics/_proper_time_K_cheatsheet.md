# Proper-time canonical `K` and dual Dirac — cheat sheet

The proper-time formulation of relativistic quantum mechanics, established by Gill and co-authors in [[Dual_Relativistic_Quantum_Mechanics_I]] and verified at [`Roadmapping/Equation_Verification/Dual_Relativistic_Quantum_Mechanics_I.md`](../Equation_Verification/Dual_Relativistic_Quantum_Mechanics_I.md), is governed by a small set of canonical relations. Every Griffiths problem's proper-time treatment reduces to applying them. This document is the one-page reference; the per-problem documents under [`Griffiths/`](Griffiths/) cite it rather than restating the relations.

## 1. The canonical proper-time Hamiltonian (DRQM I Eq. I.6)

$$
K = \frac{H_{0}^{2}}{2 m c^{2}} + \frac{m c^{2}}{2}, \qquad H_{0} = \sqrt{c^{2}\boldsymbol\pi^{2} + m^{2} c^{4}}.
$$

`H_{0}` is the relativistic free-particle Hamiltonian. The canonical `K` is the proper-time conjugate energy, taking the place of the standard relativistic `H` in the proper-time formulation.

For a charged particle in an external EM field (Gaussian units), `\boldsymbol\pi = \mathbf{p} - q\mathbf{A}/c` is the minimal-coupling momentum.

## 2. Non-relativistic reduction

For `\boldsymbol\pi \ll m c`, expand `H_{0}` in powers of `\pi/(mc)`:

$$
H_{0} = m c^{2}\sqrt{1 + \frac{\boldsymbol\pi^{2}}{m^{2} c^{2}}} = m c^{2}\left[1 + \frac{\boldsymbol\pi^{2}}{2 m^{2} c^{2}} - \frac{\boldsymbol\pi^{4}}{8 m^{4} c^{4}} + \ldots\right] = m c^{2} + \frac{\boldsymbol\pi^{2}}{2 m} - \frac{\boldsymbol\pi^{4}}{8 m^{3} c^{2}} + \ldots
$$

Substituting into `K = H_{0}^{2}/(2 m c^{2}) + m c^{2}/2` and dropping the constant `(m c^{2}) + m c^{2}/2 = (3/2) m c^{2}` rest-energy offset (which is conventional in non-relativistic QM),

$$
K - \text{const} = \frac{\boldsymbol\pi^{2}}{2 m} + \mathcal{O}\!\left(\frac{\pi^{4}}{m^{3} c^{2}}\right).
$$

This is the **standard Schrödinger kinetic energy** at leading order. The proper-time formulation reduces to ordinary non-relativistic quantum mechanics in the `u \ll c` limit.

## 3. Substitution rules summary

| Relation | Form | Origin | Application |
|---|---|---|---|
| Canonical `K` | `K = H_{0}^{2}/(2 m c^{2}) + m c^{2}/2` | DRQM I Eq. I.6 ✅ | Proper-time Hamiltonian; replaces standard `H` |
| Free `H_{0}` | `H_{0} = \sqrt{c^{2}\boldsymbol\pi^{2} + m^{2} c^{4}}` | DRQM I Eq. I.5 | Relativistic kinetic + rest energy |
| Minimal coupling | `\boldsymbol\pi = \mathbf{p} - q\mathbf{A}/c` | standard | EM coupling unchanged from textbook |
| Non-rel reduction | `K - (3/2)mc^{2} \to \boldsymbol\pi^{2}/(2m)` | from §2 above | Recovers Schrödinger kinetic energy |
| Velocity duality | `b² = c² + u²`, `u/b = w/c` | Gill–Zachary Maxwell paper | Same as #42 campaign |
| Dual Dirac equation | DRQM I Eqs. II.1–II.3 ✅ | DRQM I §II | Spin-½ relativistic equation; reduces to Pauli for `u \ll c` |

## 4. Flagged finding — DRQM I §III.D anomalous-`g`-factor (`r_e`) issue

[`FINDINGS_for_author_review.md`](../Equation_Verification/FINDINGS_for_author_review.md) records that DRQM I §III.D's stated `r_e \approx 0.499857150068631 r_0` gives `g = -2.0005714`, not the claimed experimental `-2.00231930436256`. The required `r_e` for the experimental value would be `r_e \approx 0.499420510 r_0`.

For Griffiths problems whose proper-time derivation invokes the anomalous-`g`-factor numerics — **expected primarily in Ch. 7 fine-structure problems and possibly in Ch. 4 hydrogen** — a **branched treatment** applies per the per-problem template: derive once with the as-published `r_e` and once with the corrected value, record both verdicts.

Most non-relativistic Griffiths problems do not invoke `r_e` and the branched treatment does not engage.

## 5. Pedagogical note on `K` vs `H`

The proper-time `K` is *not* the same operator as the standard relativistic Hamiltonian `H_{0}`. The relation `K = H_{0}^{2}/(2mc²) + mc²/2` is a quadratic-in-`H_{0}` definition, which has both technical and physical motivations recorded in DRQM I §I.B. For the campaign's purposes:

- **Spectra match in the non-relativistic limit.** Eigenvalues of `K` and `H_{0}^{2}/(2mc²) + mc²/2` agree, and at leading order in `\pi/(mc)` both reduce to non-rel Schrödinger eigenvalues plus rest-energy offsets.
- **Spectra differ at relativistic order.** The eigenvalue spectrum of `K` differs from `H_{0}` itself at `\mathcal{O}((v/c)^{2})` — this is the operational signature in Ch. 7 fine structure.
- **Time evolution is governed by `K`** in the proper-time formulation: `i\hbar\,\partial_{\tau}\psi = K\psi`. In the non-relativistic limit `\partial_{\tau} \to \partial_{t}` and standard Schrödinger evolution recovers.

The reduction is mechanical and the per-problem documents document it explicitly for each Ch. 1–3 problem (and where relevant in Ch. 4+).
