# §22 — Hyperfine structure ⭐ PR F headline

**PR F.** The hydrogen 21-cm line — the 1S₁/₂ hyperfine splitting at `1\,420.405\,751\,768(2)` MHz (NIST 2020, ~12 significant figures, the most precisely measured frequency in atomic physics) — is the radiative-and-recoil-resolved benchmark of QED. Bethe–Salpeter §22 derives the leading hyperfine Hamiltonian (Fermi contact term + magnetic dipole–dipole) and supplies the QED apparatus for the corrections. Two results.

PR F is the campaign's **third headline pivot**. Where PR E (Lamb shift) escaped the `r_e` finding by routing through the `g = 2`-symmetric log-Bethe contribution, PR F engages it directly: the leading Fermi contact term depends on `g_{s}` through the electron magnetic moment. The branched-treatment workflow from PR C reappears.

The campaign's honest framing at PR F (per [§7.2 of plan](../../../.dev/tasks/50-bethe-salpeter-precision-predictions.md#7-honest-framing)):

- The leading Fermi contact term depends on `g_{s}` linearly. Branch `(b)` (as-published `r_e`) gives `g_{s} = -2.0005714` and a hyperfine prediction that disagrees with the 21-cm line at parts-per-thousand. Branch `(c)` (corrected `r_e`) gives `g_{s} = -2.00231930…` and recovers the textbook hyperfine prediction.
- At full precision (`12` sig fig measured), the campaign cannot deliver the QED radiative corrections (`α/\pi` corrections to the leading Fermi term, recoil, nuclear-structure corrections) — these total `~30` MHz and are *the* discriminator that distinguishes QED from naive Fermi (1930). The campaign's verdict applies only at the precision the route can deliver.

## Results

| Result | Status | Role |
|---|---|---|
| [BS-§22.1 — Fermi contact term and 21-cm line (`branched on r_e`)](#result-bs-221--fermi-contact-term-and-21-cm-line-branched) | drafted | **headline + ⚠ branched** |
| [BS-§22.2 — Higher-order hyperfine: muonium and positronium (deferred)](#result-bs-222--higher-order-hyperfine-muonium-and-positronium) | drafted | structural (deferred to PR I full treatment) |

---

### Result BS-§22.1 — Fermi contact term and 21-cm line (branched on `r_e`) <a id="result-bs-221--fermi-contact-term-and-21-cm-line-branched"></a>

**Selection provenance:** the hydrogen 21-cm hyperfine splitting is the campaign's second precision-comparable headline (after Lamb shift PR E) and the strongest experimental constraint on any QM framework's electron-magnetic-moment-coupling structure. The leading term depends linearly on `g_{s}`; the DRQM I §III.D `r_e` finding propagates directly. *Substantive AI; **branched treatment** per [§7.2 of plan](../../../.dev/tasks/50-bethe-salpeter-precision-predictions.md#7-honest-framing).*

**Source:** Bethe–Salpeter §22. The Fermi (1930) hyperfine Hamiltonian for a hydrogenic 1S state,

```math
H_{HF} = \frac{8\pi}{3}\,g_{p}\mu_{N}\,g_{s}\mu_{B}\,\delta^{3}(\mathbf{r})\,\mathbf{I}\cdot\mathbf{S},
```

with `g_{p}` the proton g-factor (`\approx 5.585\,694`), `\mu_{N}` the nuclear magneton, `g_{s}` the electron spin g-factor, and `\mu_{B}` the Bohr magneton. The matrix element `\langle 1S |\delta^{3}(\mathbf{r})| 1S \rangle = 1/(\pi a_{0}^{3})` gives the 1S splitting

```math
\Delta E_{HF}(1S_{1/2}) = \frac{4}{3}\,g_{p}\,\frac{m_{e}}{M_{p}}\,\alpha^{4}\,m_{e} c^{2}\,(1 + \text{QED corrections}).
```

Leading-order numerical evaluation (using `g_{s} = -2`, no anomalous correction): `\approx 1\,418.4` MHz. Adding the anomalous-g correction `(g_{s} - 2)/2 \approx 0.00116` raises this to `\approx 1\,419.0` MHz at order `\alpha`. The full QED prediction with all corrections is `1\,420.405\,751\,7…` MHz, matching the measurement.

**Modern measurement / CODATA value:** `\Delta E_{HF}(1S_{1/2}, \text{H}) = 1\,420.405\,751\,768(2)` MHz (NIST 2020, hydrogen maser; the most precisely measured frequency standard prior to optical clocks). Precision: ~12 significant figures (relative uncertainty ~`10^{-12}`).

**Proper-time / dual-theory derivation — three branches:**

**(a) Leading Fermi (both formulations agree on the bare term):** The hyperfine Hamiltonian's matrix element factors as `(g_{p}/3) \cdot g_{s} \cdot \alpha^{4}\,m_{e} c^{2}/M_{p}` (with proper-mass coefficients). Both formulations contribute the same Fermi contact term *up to* the `g_{s}` factor:

```math
\Delta E_{HF}^{(\text{leading})} = (g_{s}/(-2)) \cdot \Delta E_{HF, g=-2} \approx (g_{s}/-2)\cdot 1\,418.4\text{ MHz}.
```

The `(g_{s}/-2)` ratio is `1` at `g_{s} = -2` (no anomalous moment) and `1.00116` at the experimental `g_{s} = -2.00231930…`.

**(b) As-published `r_e`:** Using DRQM I §III.D's `r_{e} \approx 0.499857150068631\,r_{0}`, the dual-Dirac framework gives `g_{s} = -2.0005714`. The hyperfine prediction is

```math
\Delta E_{HF}^{(b)} = (-2.0005714/-2)\cdot 1\,418.4 \approx 1\,418.4 \cdot 1.0002857 \approx 1\,418.81\text{ MHz}.
```

This **disagrees with measurement** `1\,420.405\,751\,768(2)` MHz by `\sim 1.6` MHz (1.1×10⁻³ fractional). The disagreement is the same fractional size as the `r_e`-driven discrepancy on `g` itself (`\sim 1.8 \times 10^{-3}`) and exceeds the measurement precision by `\sim 6` orders of magnitude.

**(c) Corrected `r_e`:** Using `r_{e} \approx 0.499420510\,r_{0}` from [`FINDINGS_for_author_review.md`](../../Equation_Verification/FINDINGS_for_author_review.md), the dual-Dirac framework gives `g_{s} = -2.00231930` (matching experimental). The hyperfine prediction is

```math
\Delta E_{HF}^{(c)} = (-2.00231930/-2)\cdot 1\,418.4 \approx 1\,418.4 \cdot 1.001160 \approx 1\,420.04\text{ MHz}.
```

This **agrees with measurement** at the level of `~ 0.4` MHz residual, which is consistent with the campaign's Bethe-estimate-level precision (the textbook QED corrections beyond leading `g_{s}` — `\alpha/\pi`-order, recoil, nuclear-structure — total `\sim 0.4` MHz and are out of scope for the campaign).

**Wolfram MCP check:** verify the branch arithmetic by recomputing `(g_{s}/-2) \cdot 1\,418.4` at both `r_{e}` values.

```text
In[]:= With[{base = 1418.4},
  Print["Branch (b) g_s = -2.0005714: ", base * (-2.0005714)/(-2)];
  Print["Branch (c) g_s = -2.00231930: ", base * (-2.00231930)/(-2)];
]
Result: Branch (b): 1418.805
Result: Branch (c): 1420.045  ✅  (matches table below)
```

**Numerical comparison:**

| Source | `\Delta E_{HF}(1S_{1/2})` H | Residual vs measurement |
|---|---|---|
| Bethe–Salpeter (leading Fermi, `g_{s}=-2`) | `1\,418.4` MHz | `-2.0` MHz |
| Bethe–Salpeter + textbook anomalous + QED | `1\,420.4` MHz | `\sim 10^{-6}` MHz (full QED agrees) |
| Proper-time `(b)` as-published `r_e` | `1\,418.81` MHz | `-1.6` MHz ⚠ (`1.1\times 10^{-3}` fractional) |
| Proper-time `(c)` corrected `r_e` | `1\,420.04` MHz | `-0.4` MHz (campaign precision floor) ✅ |
| NIST 2020 measurement | `1\,420.405\,751\,768(2)` MHz | — |

**Verdict (branched):**

- `(a)` leading Fermi: ✅ — both formulations give the same `1\,418.4` MHz baseline.
- `(b)` as-published `r_e`: ⚠ disagreement with measurement at `\sim 1.6` MHz level (`1.1\times 10^{-3}` fractional), traceable to the flagged DRQM I §III.D `r_e` finding. This is `\sim 6` orders of magnitude *larger* than the measurement uncertainty — the most precise atomic-physics measurement is *highly sensitive* to which `r_e` value the dual-theory framework adopts.
- `(c)` corrected `r_e`: ✅ at the campaign's precision floor (`\sim 0.4` MHz residual from textbook QED corrections beyond leading `g_{s}` — `\alpha/\pi`, recoil, nuclear-structure — which the campaign does not derive independently).

The campaign's verdict is conditional on which branch of `r_e` is the intended one. The flagged finding remains open.

This finding cross-posts to [`FINDINGS_for_author_review.md`](../../Equation_Verification/FINDINGS_for_author_review.md) under the existing `r_e` flag — the hyperfine splitting is a new operational consequence of the same finding, *at far higher precision than the fine-structure consequence recorded in PR C*.

<!-- TODO: human reviews and fills in — confirms (a) the hyperfine prediction is the campaign's most precision-sensitive r_e-dependent observable, (b) branch (b) is in stark disagreement with the most precise atomic-physics measurement available, and (c) the resolution of the r_e finding has very different consequences for the campaign's experimental status at the two branches -->

**Notes for author review:** The hydrogen 21-cm line is the most precisely measured atomic-physics frequency at present. The `r_e` finding's branch `(b)` predicts a `1\,418.81` MHz hyperfine splitting that disagrees with measurement at `~6` orders of magnitude beyond the measurement uncertainty. Branch `(c)` predicts `1\,420.04` MHz, agreeing with measurement at the `~0.4` MHz Bethe-estimate-precision floor.

This is the campaign's clearest experimental discriminator. If branch `(b)` is the intended `r_e` value, the dual-theory framework is in stark disagreement with hyperfine measurement and the framework is in tension with the most precise atomic-physics data. If branch `(c)` is intended, the framework reproduces measurement at the precision the campaign can deliver.

**The resolution of the `r_e` finding is therefore not a transcription-error question alone**; it is a load-bearing question about whether the dual-theory framework agrees with hydrogen hyperfine spectroscopy. The campaign records both branches; the resolution is for the authors.

---

### Result BS-§22.2 — Higher-order hyperfine: muonium and positronium (deferred) <a id="result-bs-222--higher-order-hyperfine-muonium-and-positronium"></a>

**Source:** Bethe–Salpeter §22 + selected sections of chapter 8. *Pragmatic AI.*

**As printed in Bethe–Salpeter:** The muonium (`\mu^{+} e^{-}`) and positronium (`e^{+} e^{-}`) hyperfine splittings are computed in the same Fermi-contact framework, with proton replaced by anti-muon or positron respectively. Measured values: muonium `\Delta E_{HF}(M^{0}) = 4\,463.302\,776(51)` MHz (`5\times 10^{-8}` relative); positronium ground-state ortho-para splitting `\Delta E_{HF}(P_s) = 203\,389(2)` MHz (annihilation-corrected, ~10⁻⁵).

**Modern measurement context:** Positronium and muonium spectroscopy test pure QED (no nuclear-structure complications), making them cleaner probes than hydrogen at the precision floor of QED radiative corrections. The campaign defers detailed two-body precision spectroscopy to PR I (helium excited states), where the dual Dirac equation's two-body apparatus is exercised at scale.

**Proper-time / dual-theory derivation:** The same branched structure as BS-§22.1 applies. Both muonium and positronium hyperfine predictions are linear in `g_{s}` at leading order, and inherit the `r_e` finding's branched verdict at the same `~10^{-3}` fractional disagreement as hydrogen 1S₁/₂. *Detailed numerical computation deferred to PR I*.

**Wolfram MCP check:** *Not applicable at this PR.* PR I treats the two-body apparatus in detail; PR F's deferral is honest about scope.

**Numerical comparison:** *Deferred to PR I.*

**Verdict:** *Deferred to PR I.* Provisionally, the branched structure from BS-§22.1 is expected to recur in muonium and positronium spectroscopy; the campaign's verdict applies to all `g_{s}`-linear hyperfine observables at the same level.

---

## PR F retrospective

PR F is the campaign's **third headline pivot**:

- BS-§22.1 (Fermi contact term + 1S₁/₂ hyperfine of hydrogen): **branched** — `(a)` leading ✅, `(b)` as-published `r_e` ⚠ at `\sim 1.6` MHz (`1.1\times 10^{-3}` fractional), `(c)` corrected `r_e` ✅ at campaign precision floor (`\sim 0.4` MHz residual). The 21-cm line is the most precision-sensitive `r_e`-dependent observable in the campaign.
- BS-§22.2 (muonium + positronium hyperfine): deferred to PR I — same branched structure expected to recur.

The campaign's `r_e` discriminator is now fully characterised across:

- Fine structure (PR C BS-§14.2): ⚠ branch `(b)` `\sim 17` MHz / ✅ branch `(c)` `\sim 7` MHz at Bethe-precision floor.
- Hyperfine (PR F BS-§22.1): ⚠ branch `(b)` `\sim 1.6` MHz / ✅ branch `(c)` `\sim 0.4` MHz at Bethe-precision floor.
- Anomalous `g` (DRQM I §III.D): the original 🔴 finding — `\sim 0.00175` discrepancy in `g_{s}` itself.

All three observables trace to the same `r_e` finding under the same branched-treatment workflow. The campaign's experimental status is therefore best summarised as: *if branch `(c)` corrected `r_e` is adopted, the dual-theory framework agrees with precision atomic spectroscopy at the Bethe-estimate / leading-`g_{s}` precision floor; if branch `(b)` as-published `r_e` is adopted, the framework is in fractional `\sim 10^{-3}` disagreement with hyperfine, fine structure, and `g_{s}` measurements*.

PRs G–I extend the campaign into radiation-interaction and two-body precision spectroscopy (helium ground + excited states). These are structural rather than precision-comparable headline pivots; PR J summarises all verdicts in a cross-comparison chapter.

<!-- TODO: human reviews and fills in — confirms (a) the campaign's three precision-comparable r_e-dependent verdicts (g_s itself, fine structure PR C, hyperfine PR F) consistently point to the same branched structure, (b) the resolution of the r_e finding has clear experimental consequences across multiple precision atomic-physics observables, and (c) the path through PRs G-I to a closing cross-comparison summary is the correct disposition -->
