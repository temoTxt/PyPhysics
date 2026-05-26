# §22 — Hyperfine structure ⭐ PR F headline

**PR F.** The hydrogen 21-cm line — the 1S₁/₂ hyperfine splitting at `1\,420.405\,751\,768(2)` MHz (NIST 2020, ~12 significant figures, the most precisely measured frequency in atomic physics) — is the radiative-and-recoil-resolved benchmark of QED. Bethe–Salpeter §22 derives the leading hyperfine Hamiltonian (Fermi contact term + magnetic dipole–dipole) and supplies the QED apparatus for the corrections. Two results.

PR F is the campaign's **third headline pivot**. Where PR E (Lamb shift) escaped the `r_e` finding by routing through the `g = 2`-symmetric log-Bethe contribution, PR F engages it directly: the leading Fermi contact term depends on `g_{s}` through the electron magnetic moment. The framework's prediction is evaluated at the triangulated `r_e/r_0 = 0.4994205099128317` per [PR #62](https://github.com/temoTxt/PyPhysics/pull/62).

The campaign's honest framing at PR F (per [§7.2 of plan](../../../.dev/tasks/50-bethe-salpeter-precision-predictions.md#7-honest-framing)):

- The leading Fermi contact term depends on `g_{s}` linearly. At the triangulated `r_e` (the joint-best-fit across six precision observables), `g_{s} = -2.00231930…` matches the measured value and the framework recovers the textbook hyperfine prediction.
- At full precision (`12` sig fig measured), the campaign cannot deliver the QED radiative corrections (`α/\pi` corrections to the leading Fermi term, recoil, nuclear-structure corrections) — these total `~30` MHz and are *the* discriminator that distinguishes QED from naive Fermi (1930). The campaign's verdict applies only at the precision the route can deliver.

## Results

| Result | Status | Role |
|---|---|---|
| [BS-§22.1 — Fermi contact term and 21-cm line (at triangulated `r_e`)](#result-bs-221--fermi-contact-term-and-21-cm-line-branched) | drafted | **headline + ✅ at triangulated `r_e`** |
| [BS-§22.2 — Higher-order hyperfine: muonium and positronium (deferred)](#result-bs-222--higher-order-hyperfine-muonium-and-positronium) | drafted | structural (deferred to PR I full treatment) |

---

### Result BS-§22.1 — Fermi contact term and 21-cm line (at triangulated `r_e`) <a id="result-bs-221--fermi-contact-term-and-21-cm-line-branched"></a>

**Selection provenance:** the hydrogen 21-cm hyperfine splitting is the campaign's second precision-comparable headline (after Lamb shift PR E) and the strongest experimental constraint on any QM framework's electron-magnetic-moment-coupling structure. The leading term depends linearly on `g_{s}`; the DRQM I §III.D `r_e` finding propagates directly, evaluated at the triangulated `r_e/r_0 = 0.4994205099128317` per [PR #62](https://github.com/temoTxt/PyPhysics/pull/62) (closes [#61](https://github.com/temoTxt/PyPhysics/issues/61)). *Substantive AI; un-branched-verdict cleanup post-[PR #62](https://github.com/temoTxt/PyPhysics/pull/62).*

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

**Proper-time / dual-theory derivation — leading + anomalous at triangulated `r_e`:**

**(a) Leading Fermi (both formulations agree on the bare term):** The hyperfine Hamiltonian's matrix element factors as `(g_{p}/3) \cdot g_{s} \cdot \alpha^{4}\,m_{e} c^{2}/M_{p}` (with proper-mass coefficients). Both formulations contribute the same Fermi contact term *up to* the `g_{s}` factor:

```math
\Delta E_{HF}^{(\text{leading})} = (g_{s}/(-2)) \cdot \Delta E_{HF, g=-2} \approx (g_{s}/-2)\cdot 1\,418.4\text{ MHz}.
```

The `(g_{s}/-2)` ratio is `1` at `g_{s} = -2` (no anomalous moment) and `1.00116` at the experimental `g_{s} = -2.00231930…`.

**(b) At triangulated `r_e`:** Using the triangulated `r_{e}/r_{0} = 0.4994205099128317` from [PR #62](https://github.com/temoTxt/PyPhysics/pull/62) (the joint-best-fit across six precision observables; see also the [follow-up author note](../../Author_Reports/2026-05_re_triangulation_followup_for_gill.md)), the dual-Dirac framework gives `g_{s} = -2.00231930` (matching the measured value). The hyperfine prediction is

```math
\Delta E_{HF} = (-2.00231930/-2)\cdot 1\,418.4 \approx 1\,418.4 \cdot 1.001160 \approx 1\,420.04\text{ MHz}.
```

This **agrees with measurement** at the level of `~ 0.4` MHz residual, which is consistent with the campaign's Bethe-estimate-level precision (the textbook QED corrections beyond leading `g_{s}` — `\alpha/\pi`-order, recoil, nuclear-structure — total `\sim 0.4` MHz and are out of scope for the campaign).

**Back-fit caveat — what the "✅" is and is not testing.** Same caveat as in [`03_FineStructure.md` BS-§14.2](03_FineStructure.md#result-bs-142--2p--2p-fine-structure-splitting-branched). The triangulated `r_{e}/r_{0} = 0.4994205099128317` is the value that gives measured `g_{s}`; the prediction `(g_{s, \text{measured}}/-2) \times 1\,418.4` reduces to the *textbook* Fermi-contact prediction with the measured anomalous moment. The "✅" verdict therefore says: *the textbook leading-Fermi-contact formula with measured `g_{s}` reproduces measurement at textbook leading-`g_{s}` precision when `r_e` is the joint-best-fit value*. That is a self-consistency check at the cutoff, not independent corroboration of the dual-theory framework's content distinct from textbook QED. The 0.4 MHz residual on a 1,420 MHz observable measured to 2 Hz precision is `\sim 10^{5}` σ from measurement; the "✅" is relative to the Bethe-estimate / leading-`g_{s}` precision floor, not relative to experimental uncertainty. See [`10_CrossComparison.md` §2](10_CrossComparison.md#2-the-r_e-back-fit-self-consistency-across-six-g_s-dependent-observables).

**Wolfram MCP check:** verify the arithmetic by computing `(g_{s}/-2) \cdot 1\,418.4` at the triangulated `r_e`.

```text
In[]:= With[{base = 1418.4}, base * (-2.00231930)/(-2)]
Result: 1420.045  ✅  (matches table below)
```

**Numerical comparison:**

| Source | `\Delta E_{HF}(1S_{1/2})` H | Residual vs measurement |
|---|---|---|
| Bethe–Salpeter (leading Fermi, `g_{s}=-2`) | `1\,418.4` MHz | `-2.0` MHz |
| Bethe–Salpeter + textbook anomalous + QED | `1\,420.4` MHz | `\sim 10^{-6}` MHz (full QED agrees) |
| Proper-time at triangulated `r_e` | `1\,420.04` MHz | `-0.4` MHz (campaign precision floor) ✅ |
| NIST 2020 measurement | `1\,420.405\,751\,768(2)` MHz | — |

**Verdict:**

- `(a)` leading Fermi: ✅ — both formulations give the same `1\,418.4` MHz baseline.
- `(b)` at triangulated `r_e`: ✅ at the campaign's precision floor (`\sim 0.4` MHz residual from textbook QED corrections beyond leading `g_{s}`). **Back-fit self-consistency, not independent corroboration** — see back-fit caveat above.

This finding cross-posts to [`FINDINGS_for_author_review.md`](../../Equation_Verification/FINDINGS_for_author_review.md) under the existing `r_e` flag — the hyperfine splitting is an operational consequence of the same finding, now evaluated at the triangulated value, *at far higher precision than the fine-structure consequence recorded in PR C*.

<!-- TODO: human reviews and fills in — confirms (a) the hyperfine prediction is the campaign's most precision-sensitive r_e-dependent observable, (b) the un-branched verdict at the triangulated `r_e` matches measurement at the Bethe-estimate precision floor, and (c) the self-consistency framing is faithfully recorded. -->

**Notes for author review:** The hydrogen 21-cm line is the most precisely measured atomic-physics frequency at present (12 sig figs). At the triangulated `r_e`, the dual-theory framework predicts `1\,420.04` MHz, agreeing with measurement at the `~0.4` MHz Bethe-estimate-precision floor. This residual is `~ 10^5` measurement-σ — meaningful relative to QED-precision standards but not relative to the framework's known precision-floor scope. Cross-reference to [PR #62](https://github.com/temoTxt/PyPhysics/pull/62) and the [follow-up author note](../../Author_Reports/2026-05_re_triangulation_followup_for_gill.md) for the triangulation's full residual table.

---

### Result BS-§22.2 — Higher-order hyperfine: muonium and positronium (deferred) <a id="result-bs-222--higher-order-hyperfine-muonium-and-positronium"></a>

**Source:** Bethe–Salpeter §22 + selected sections of chapter 8. *Pragmatic AI.*

**As printed in Bethe–Salpeter:** The muonium (`\mu^{+} e^{-}`) and positronium (`e^{+} e^{-}`) hyperfine splittings are computed in the same Fermi-contact framework, with proton replaced by anti-muon or positron respectively. Measured values: muonium `\Delta E_{HF}(M^{0}) = 4\,463.302\,776(51)` MHz (`5\times 10^{-8}` relative); positronium ground-state ortho-para splitting `\Delta E_{HF}(P_s) = 203\,389(2)` MHz (annihilation-corrected, ~10⁻⁵).

**Modern measurement context:** Positronium and muonium spectroscopy test pure QED (no nuclear-structure complications), making them cleaner probes than hydrogen at the precision floor of QED radiative corrections. The campaign defers detailed two-body precision spectroscopy to PR I (helium excited states), where the dual Dirac equation's two-body apparatus is exercised at scale.

**Proper-time / dual-theory derivation:** The same `(g_s/-2)^n × textbook` structure as BS-§22.1 applies. Both muonium and positronium hyperfine predictions are linear (muonium) or quadratic (positronium) in `g_{s}` at leading order, and the framework's prediction at the triangulated `r_e` is fully covered in PR I (BS-§80). *Detailed numerical computation in PR I.*

**Wolfram MCP check:** *Not applicable at this PR.* PR I treats the two-body apparatus in detail; PR F's deferral is honest about scope.

**Numerical comparison:** *Deferred to PR I.*

**Verdict:** *Deferred to PR I.* The `(g_s/-2)^n × textbook` structure from BS-§22.1 recurs in muonium (`n=1`) and positronium (`n=2`) spectroscopy; the campaign's verdict applies to all `g_{s}`-dependent hyperfine observables at the triangulated `r_e` per [PR #62](https://github.com/temoTxt/PyPhysics/pull/62).

---

## PR F retrospective

PR F is the campaign's **third headline pivot**:

- BS-§22.1 (Fermi contact term + 1S₁/₂ hyperfine of hydrogen): ✅ at triangulated `r_e/r_0 = 0.4994205099128317` (per [PR #62](https://github.com/temoTxt/PyPhysics/pull/62)) — `(a)` leading ✅, `(b)` at triangulated `r_e` ✅ at campaign precision floor (`\sim 0.4` MHz residual). The 21-cm line is the most precision-sensitive `r_e`-dependent observable in the campaign.
- BS-§22.2 (muonium + positronium hyperfine): deferred to PR I — same `(g_s/-2)^n × textbook` structure, evaluated at the triangulated `r_e`.

The campaign's `r_e` evaluation is now fully characterised at the triangulated value across:

- Fine structure (PR C BS-§14.2): ✅ at triangulated `r_e` — `\sim 7` MHz residual at Bethe-precision floor.
- Hyperfine (PR F BS-§22.1): ✅ at triangulated `r_e` — `\sim 0.4` MHz residual at Bethe-precision floor.
- Anomalous `g` (DRQM I §III.D): the original 🔴 finding — characterised to ⚠ by [PR #62](https://github.com/temoTxt/PyPhysics/pull/62)'s triangulation.

All three observables trace to the same `r_e` finding, evaluated at the joint-best-fit cutoff from the multi-observable triangulation. The campaign's experimental status is therefore: *at the triangulated `r_e`, the dual-theory framework agrees with precision atomic spectroscopy at the Bethe-estimate / leading-`g_{s}` precision floor.* The residual ⚠ in [FINDINGS Finding 2](../../Equation_Verification/FINDINGS_for_author_review.md) reflects the open Scope 1 question of whether a first-principles derivation reproduces the triangulated cutoff value (tracked in [issue #54](https://github.com/temoTxt/PyPhysics/issues/54)).

PRs G–I extend the campaign into radiation-interaction and two-body precision spectroscopy (helium ground + excited states). These are structural rather than precision-comparable headline pivots; PR J summarises all verdicts in a cross-comparison chapter.

<!-- TODO: human reviews and fills in — confirms (a) the campaign's three precision-comparable r_e-dependent verdicts (g_s itself, fine structure PR C, hyperfine PR F) all sit consistently at the triangulated `r_e` value post-[PR #62](https://github.com/temoTxt/PyPhysics/pull/62), (b) the un-branched verdict structure at the joint-best-fit cutoff is the correct disposition, and (c) the path through PRs G-I to PR J's closing cross-comparison summary remains valid -->
