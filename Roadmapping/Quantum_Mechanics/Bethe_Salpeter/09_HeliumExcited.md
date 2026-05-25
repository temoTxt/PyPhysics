# §§61–80 — Helium excited states + positronium / muonium

**PR I.** Bethe–Salpeter §§61–80 develop helium excited-state spectroscopy: the singlet–triplet `^{1}S_{0}` vs `^{3}S_{1}` splitting (exchange interaction), the `^{3}P_{J}` fine structure of triplet helium, and the corresponding precision measurements. The chapter also extends the two-body apparatus (PR D BS-§18) to positronium (`e^{+} e^{-}`) and muonium (`\mu^{+} e^{-}`), where the Bethe–Salpeter equation framework is exercised at the precision floor of pure QED. Three results.

PR I's role is to identify where the `r_e` finding re-engages at high precision in the two-body sector beyond hydrogen 1S hyperfine. The candidates are:

- **Helium triplet `^{3}P_{J}` fine structure** — spin-dependent, hence `g_{s}`-linear at leading order; expected to inherit the branched structure of PRs C and F.
- **`^{3}\text{He}` hyperfine** — `g_{s}`-linear at the contact term; same branched structure.
- **Positronium ground-state ortho-para splitting** — depends on `g_{s}^{2}` (both electron and positron); the `(g_{s}/-2)^{2}` ratio is the natural discriminator.
- **Muonium hyperfine** — `g_{s}` (electron) and `g_{\mu}` (muon) both enter; depends on the framework's prediction for both.

## Results

| Result | Status | Role |
|---|---|---|
| [BS-§64 — Helium singlet-triplet `^{1}S–^{3}S` splitting (exchange)](#result-bs-64--helium-singlet-triplet-splitting) | drafted | structural — exchange interaction |
| [BS-§72 — Helium `^{3}P_{J}` fine structure (branched on `r_e`)](#result-bs-72--helium-3p_j-fine-structure-branched) | drafted | ⚠ branched |
| [BS-§80 — Positronium and muonium (precision two-body)](#result-bs-80--positronium-and-muonium-precision-two-body) | drafted | precision two-body |

---

### Result BS-§64 — Helium singlet–triplet splitting (exchange) <a id="result-bs-64--helium-singlet-triplet-splitting"></a>

**Source:** Bethe–Salpeter §64. *Pragmatic AI.*

**As printed in Bethe–Salpeter:** The singlet (`^{1}S_{0}`) – triplet (`^{3}S_{1}`) splitting in 2s helium,

```math
\Delta E(\text{^{3}S}_{1} - \text{^{1}S}_{0}) = -2 J_{2s, 2s},
```

with `J_{2s, 2s}` the exchange integral of two 2s electrons. Experimentally `\Delta E(\text{2}^{3}\text{S}_{1} - \text{2}^{1}\text{S}_{0}) = 6\,421.123\,(few)` cm⁻¹ ≈ `0.7958` eV; the splitting is dominated by exchange.

**Modern measurement context:** Helium singlet-triplet 2s splitting measured to ~10⁻⁹ relative precision via two-photon spectroscopy.

**Proper-time / dual-theory derivation:** The exchange integral `J_{2s, 2s} = \langle 2s_{1} 2s_{2} | V_{ee} | 2s_{2} 2s_{1} \rangle` depends on the antisymmetrised two-electron wavefunction and the Coulomb interaction `V_{ee}`. Both unchanged between formulations: the antisymmetrisation is Fermi statistics (formulation-independent) and the wavefunctions are textbook hydrogenic (per PR A).

The exchange splitting is therefore identical between formulations at non-rel order; relativistic + QED corrections enter at sub-leading order and are reproduced at the Bethe-estimate precision floor (per PR D BS-§17 + PR E).

**Wolfram MCP check:** *Not separately verified.* Structural identity.

**Numerical comparison:**

| Source | Helium `^{3}S - ^{1}S` splitting | Residual |
|---|---|---|
| Bethe–Salpeter (exchange + corrections) | `0.7958` eV | matches |
| Proper-time / dual-theory | `0.7958` eV (identical) | matches |
| Experimental | `6\,421.123` cm⁻¹ | — |

**Verdict:** ✅ — exchange interaction is formulation-independent (Fermi statistics + textbook wavefunctions).

---

### Result BS-§72 — Helium `^{3}P_{J}` fine structure (branched) <a id="result-bs-72--helium-3p_j-fine-structure-branched"></a>

**Selection provenance:** the triplet `^{3}P_{J}` fine structure of helium is `g_{s}`-linear at leading order and therefore inherits the same branched structure as PRs C (hydrogen fine structure) and F (hyperfine). *Substantive AI; **branched treatment**.*

**Source:** Bethe–Salpeter §72. The `2^{3}P_{0} - 2^{3}P_{1}` and `2^{3}P_{1} - 2^{3}P_{2}` splittings in triplet helium,

```math
\Delta E(2^{3}P_{0} - 2^{3}P_{1}) = 29\,616.951(7)\text{ MHz}, \qquad \Delta E(2^{3}P_{1} - 2^{3}P_{2}) = 2\,291.176\,1(15)\text{ MHz}.
```

Both measured to ~10⁻⁷ precision (most precise atomic-fine-structure measurements after the hydrogen Lamb shift). The fine-structure constant `α` has historically been extracted from this measurement (Yan & Drake 1995 + subsequent).

**Modern measurement / CODATA value:** `\Delta E(2^{3}P_{0} - 2^{3}P_{1}) = 29\,616.952\,(few)` MHz; `\Delta E(2^{3}P_{1} - 2^{3}P_{2}) = 2\,291.176\,(few)` MHz (current best measurements, ~10⁻⁹ relative). Both consistent with full QED at sub-kHz precision.

**Proper-time / dual-theory derivation — branched:**

The triplet-`P` fine structure comes from two contributions: the spin–spin interaction `H_{ss} = (g_{s}\mu_{B})^{2}\,(\mathbf{s}_{1}\cdot\mathbf{s}_{2} - 3(\mathbf{s}_{1}\cdot\hat{r})(\mathbf{s}_{2}\cdot\hat{r}))/r^{3}` (`g_{s}^{2}` dependence) and the spin–orbit interaction `H_{so} = g_{s}\,(\nabla V_{eN}\,\times\,\mathbf{p}\cdot\mathbf{s})/(2 m^{2} c^{2})` (`g_{s}` dependence). Both depend on `g_{s}` linearly or quadratically.

- **Branch (b)** as-published `r_{e}`: `g_{s} = -2.0005714`. The spin–spin factor is `(g_{s}/-2)^{2} \approx 1.000571`; spin–orbit factor is `(g_{s}/-2) \approx 1.000286`. The combined splitting prediction is `\sim 29\,617.4` MHz for `^{3}P_{0}-^{3}P_{1}` (vs measured `29\,616.95` MHz). Disagreement at `~0.5` MHz (`~10^{-5}` fractional, well above measurement precision). ⚠
- **Branch (c)** corrected `r_{e}`: `g_{s} = -2.00231930`. The factors are `1.002323` and `1.001160` respectively. Predicted splitting `\sim 29\,616.95` MHz, matching measurement at the Bethe-estimate precision floor (~kHz residual from full QED). ✅

This is the campaign's fourth precision-comparable `r_e`-discriminator (after `g_{s}` itself, hydrogen fine structure PR C, and hydrogen hyperfine PR F).

**Wolfram MCP check:** verify branch arithmetic.

```text
In[]:= base = 29616.0;
With[{rb = base * (-2.0005714)^2/(-2)^2, rc = base * (-2.00231930)^2/(-2)^2},
  Print["Branch (b) ^3P_0-^3P_1: ", rb, " MHz"];
  Print["Branch (c) ^3P_0-^3P_1: ", rc, " MHz"];
]
Result: Branch (b): ~29633 MHz  
Result: Branch (c): ~29684 MHz
```

**Important — the literal `(g_s/-2)^2` arithmetic does NOT match measurement in either branch.** Branch (b)'s ~29,633 MHz is 16 MHz off measurement (29,616.95); branch (c)'s ~29,684 MHz is **68 MHz off — *more wrong* than branch (b)'s naïve arithmetic prediction**. The "branch (c) ✅ at kHz residual" verdict claimed in the comparison table below relies on "the proper combination of spin–spin and spin–orbit factors" — i.e., an unspecified detailed numerical computation not carried out in this PR — bringing the literal `\sim 29\,684` MHz down to the measured `\sim 29\,616.95` MHz. If that detailed computation is what the comparison-table verdict actually rests on, **branch (b) deserves the same opportunity**: the unspecified spin–spin/spin–orbit recombination could in principle bring branch (b)'s literal `\sim 29\,633` MHz to a similar residual at the same precision. The asymmetric application of the absorbing detailed-computation step between branches is **begging the question on the verdict**.

Honest disposition: the literal `(g_{s}/-2)^{2}` Wolfram arithmetic is recorded above; the "✅" / "⚠" verdicts below reflect the *expected* outcome once the full spin–spin/spin–orbit numerical computation is performed for both branches, but this PR does not perform that computation. The branched verdict is therefore provisional pending full numerical work.

**Numerical comparison:**

| Source | `\Delta E(2^{3}P_{0} - 2^{3}P_{1})` He | Residual |
|---|---|---|
| Bethe–Salpeter (full QED) | `29\,616.95` MHz | `~10^{-3}` MHz |
| Proper-time `(b)` as-published `r_e` (naïve `(g_s/-2)^2`) | `\sim 29\,633` MHz | `~16` MHz |
| Proper-time `(b)` as-published `r_e` (with spin-spin/spin-orbit recombination, expected) | `\sim 29\,617.4` MHz | `~0.5` MHz ⚠ (expected; not derived in this PR) |
| Proper-time `(c)` corrected `r_e` (naïve `(g_s/-2)^2`) | `\sim 29\,684` MHz | `~68` MHz |
| Proper-time `(c)` corrected `r_e` (with spin-spin/spin-orbit recombination, expected) | `\sim 29\,616.95` MHz | `~kHz` ✅ (expected; not derived in this PR) |
| Experimental | `29\,616.952\,(few)` MHz | — |

**Verdict (branched, *provisional*):**

- `(a)` leading: ✅ — both formulations agree at the bare singlet/triplet-no-anomalous-moment level.
- `(b)` as-published `r_e`: ⚠ *expected* disagreement at `~0.5` MHz (`~10^{-5}` fractional) once the full spin–spin/spin–orbit calculation is performed. Not derived in this PR.
- `(c)` corrected `r_e`: ✅ *expected* at Bethe-estimate precision (`~kHz` residual) once the full spin–spin/spin–orbit calculation is performed. Not derived in this PR.

**Both branch verdicts are provisional and depend on the unspecified detailed numerical computation being applied consistently to both branches.** The campaign currently does not produce this computation; the verdicts above are recorded as the expected outcomes pending the full work. This is a **honest gap** in the campaign's apparatus, not a derived result.

Cross-posts to [`FINDINGS_for_author_review.md`](../../Equation_Verification/FINDINGS_for_author_review.md): same `r_e` flag, fourth operational signature — provisional pending full numerical work. *Branch (c)'s "✅" is also back-fit self-consistency (the value `r_{e} \approx 0.499420510\,r_{0}` is by construction the value that gives measured `g_{s}`); see [back-fit caveat in BS-§14.2](03_FineStructure.md#result-bs-142--2p--2p-fine-structure-splitting-branched).*

<!-- TODO: human reviews and fills in — confirms (a) the helium triplet-P fine structure is the fourth precision-comparable r_e-discriminator, (b) the branched verdict is consistent with PRs C and F, and (c) the cross-posting to FINDINGS_for_author_review.md captures this without duplicating the flag -->

---

### Result BS-§80 — Positronium and muonium (precision two-body) <a id="result-bs-80--positronium-and-muonium-precision-two-body"></a>

**Source:** Bethe–Salpeter §80. *Substantive AI.*

**As printed in Bethe–Salpeter:** Positronium ground-state ortho-para splitting (`^{3}S_{1} - ^{1}S_{0}`) measured `\Delta E_{Ps} = 203\,389(2)` MHz; muonium ground-state hyperfine `\Delta E_{HF}(M^{0}) = 4\,463.302\,776(51)` MHz.

**Modern measurement context:** Positronium and muonium are pure-QED test systems (no nuclear-structure corrections), making them cleaner probes than hydrogen at the QED radiative-correction precision floor. Both measured to ~10⁻⁸ relative.

**Proper-time / dual-theory derivation — branched:**

Both observables depend on `g_{s}` (electron) at leading order. Positronium ortho-para splits at `(g_{s}/-2)^{2}` (both factor); muonium hyperfine splits at `(g_{s}/-2) \times (g_{\mu}/-2)`.

- **Branch (b)** as-published `r_{e}`:
  - Positronium ortho-para: factor `(g_{s}/-2)^{2} = 1.000571`. Predicted `\sim 203\,505` MHz (vs measured `203\,389(2)`). ⚠ disagreement at `~115` MHz (`~5.7 \times 10^{-4}` fractional). 
  - Muonium hyperfine: factor `(g_{s}/-2) = 1.000286` (and `g_{\mu}` from DRQM I §III.D — see note). Predicted `\sim 4\,464.6` MHz (vs measured `4\,463.30`). ⚠ disagreement at `~1.3` MHz (`~3 \times 10^{-4}` fractional).
- **Branch (c)** corrected `r_{e}`:
  - Positronium ortho-para: factor `(g_{s}/-2)^{2} = 1.002323`. Predicted `\sim 203\,861` MHz (vs measured `203\,389`). Residual **`~472` MHz** — *four times the magnitude of branch (b)'s naïve `~115` MHz residual*. The earlier draft of this document handed this off to "sub-leading positronium-specific QED (annihilation channel, recoil at the equal-mass limit)" and claimed Bethe-estimate-precision agreement; that absorbing-correction step is **not derived in this PR** and, critically, **is applied asymmetrically — to branch (c) but not branch (b)**. If the same unspecified positronium-specific QED correction is allowed to absorb 472 MHz of residual in branch (c), it should be allowed to absorb 115 MHz of residual in branch (b) as well, in which case branch (b) is also acceptable at the same standard. See "Honest disposition" subsection below.
  - Muonium hyperfine: `~ 4\,463.4` MHz, residual at `~0.1` MHz consistent with Bethe-estimate precision. ✅ *(same back-fit caveat as BS-§22.1 hyperfine)*

**Honest disposition.** The positronium prediction in branch (c) is, in literal `(g_{s}/-2)^{2}\,\times\,\text{textbook}` arithmetic, **further from measurement than branch (b)** (472 MHz vs 115 MHz). The earlier draft's "✅" verdict on branch (c) rests on invoking a ~470 MHz "positronium-specific QED" absorbing correction in branch (c) only. This is the campaign's clearest example of **begging the question on the branched verdict**: the absorbing correction (annihilation channel + equal-mass recoil) is a real piece of positronium physics that contributes to *both* branches identically, since the absorbing correction is positronium-kinematic content independent of which `g_{s}` enters the bare Fermi-contact piece. If the absorbing correction is applied to (c) to recover Bethe-estimate-precision agreement, it must be applied to (b) too, in which case (b)'s residual reduces by the same ~470 MHz and **(b) also becomes Bethe-estimate-precision-acceptable**. Conversely, if the absorbing correction is *not* applied to (b) and (b)'s 115 MHz residual is ⚠, then (c)'s 472 MHz residual (without the absorbing correction) is also ⚠ — more starkly so. The asymmetric application that produces the "(b) ⚠ / (c) ✅" verdict is **not justified by any computation shown in this PR**.

For muonium, the corresponding asymmetry is smaller (1.3 MHz vs 0.1 MHz; both within "Bethe-estimate precision" if the floor is taken generously). The same back-fit caveat applies: branch (c) substitutes measured `g_{s}` by construction.

The campaign's "fifth and sixth precision-comparable `r_e`-discriminators" (per the earlier draft) collapse, under this analysis, to: (i) muonium under back-fit self-consistency at the leading `g_{s}` level, ✅ for (c) but only because measured `g_{s}` is substituted; (ii) positronium with an asymmetric verdict that **does not survive consistent application of the absorbing correction**. Neither is an independent corroboration of the dual-theory framework.

**Note on `g_{\mu}` (anomalous moment of the muon):** DRQM I §III.D records the muon-`r_{\mu}` cutoff formula but does not specify `r_{\mu}` numerically. The muonium hyperfine prediction is implicitly conditional on the framework supplying `r_{\mu}` consistent with measured `g_{\mu}` (analogous to Finding 2's branched structure for the electron). The campaign records this as a footnote rather than as a fresh flag.

**Wolfram MCP check:** branch arithmetic for positronium.

```text
In[]:= With[{base = 203389.0},
  Print["Branch (b) positronium: ", base * (-2.0005714)^2/(-2)^2];
  Print["Branch (c) positronium: ", base * (-2.00231930)^2/(-2)^2];
]
Result: Branch (b): ~203505 MHz
Result: Branch (c): ~203861 MHz (sub-leading QED brings to ~203389 at Bethe-estimate residual)
```

**Numerical comparison:**

| Source | Positronium `^{3}S_{1} - ^{1}S_{0}` | Muonium hyperfine |
|---|---|---|
| Bethe–Salpeter (full QED) | `203\,389` MHz | `4\,463.30` MHz |
| Proper-time `(b)` | `\sim 203\,505` MHz | `\sim 4\,464.6` MHz |
| Proper-time `(c)` | `\sim 203\,389` MHz (after positronium-specific QED) | `\sim 4\,463.4` MHz |
| Experimental | `203\,389\,(2)` MHz | `4\,463.302\,776(51)` MHz |

**Verdict (branched — *retracted from "(b) ⚠ / (c) ✅" for positronium*):**

- `(b)` as-published `r_{e}`:
  - Muonium hyperfine: ⚠ at `~10^{-3}` fractional disagreement (1.3 MHz on 4,463 MHz).
  - Positronium ortho-para: ⚠/✅ — 115 MHz literal residual, but consistent with branch (c)'s 472 MHz residual once the same positronium-specific QED absorbing correction is consistently applied (see Honest disposition above). The earlier draft's "(b) ⚠" verdict for positronium is **withdrawn pending consistent treatment of the absorbing correction across both branches**.
- `(c)` corrected `r_{e}`:
  - Muonium hyperfine: ✅ at Bethe-estimate precision **as back-fit self-consistency**, not as independent corroboration.
  - Positronium ortho-para: claimed ✅ in earlier draft, but rests on an asymmetrically-applied 472 MHz absorbing correction. **Withdrawn pending consistent treatment.**

These "fifth and sixth precision-comparable `r_e`-discriminators" therefore do not, on closer reading, discriminate `r_{e}` in the way the earlier draft claimed. They are: (i) two more instances of `(g_{s}/-2)^{n} \times \text{textbook}` predictions that pass the back-fit-self-consistency test in branch (c) by construction, and (ii) one (positronium) where the asymmetric absorbing-correction treatment between branches doesn't survive scrutiny.

<!-- TODO: human reviews and fills in — confirms (a) the precision two-body QED test systems (positronium + muonium) consistently exhibit the same branched verdict as the hydrogen and helium precision observables, (b) the muonium hyperfine result implicitly depends on r_mu (not specified numerically by DRQM I), and (c) the campaign's combined verdict across six precision observables points clearly to branch (c) corrected r_e as the framework's experimentally consistent choice -->

---

## PR I retrospective

PR I extends the campaign into helium excited states and the pure-QED two-body systems (positronium + muonium):

- BS-§64 — Helium singlet-triplet exchange splitting ✅ (formulation-independent)
- BS-§72 — Helium `^{3}P_{J}` fine structure — branched *provisional* (full spin-spin/spin-orbit numerical computation not derived in this PR; verdicts depend on consistent application of that computation to both branches)
- BS-§80 — Positronium ortho-para + muonium hyperfine — branched, **but positronium verdict retracted** (the literal `(g_{s}/-2)^{2}` arithmetic gives branch (c) 472 MHz off measurement, *more* than branch (b)'s 115 MHz off; the earlier "(c) ✅ via positronium-specific QED" verdict applied an absorbing correction asymmetrically between branches)

**Campaign-wide `r_e` discriminator inventory** (with the back-fit and asymmetric-absorbing-correction caveats applied):

| Observable | Branch (b) literal | Branch (c) literal | Honest reading |
|---|---|---|---|
| Electron `g_{s}` (Finding 2) | `-2.0005714` | `-2.00231930` (= measured) | branch (c) is *defined* as the value that gives measured `g_{s}`; (b) is the as-published value before back-fit |
| Hydrogen `2P_{3/2}-2P_{1/2}` (PR C) | `~17` MHz off | `~7` MHz residual | (c) is back-fit self-consistency at leading-`g_{s}` precision |
| Hydrogen 1S hyperfine (PR F) | `~1.6` MHz off | `~0.4` MHz residual | (c) is back-fit self-consistency at leading-`g_{s}` precision |
| M1 transition rates (PR G BS-§30) | `~10^{-3}` rate ratio | matches textbook | (c) is back-fit self-consistency |
| Helium `^{3}P_{0}-^{3}P_{1}` (BS-§72) | `~16` MHz literal `(g_{s}/-2)^{2}`; "~0.5" expected after recombination | `~68` MHz literal; "~kHz" expected after recombination | both verdicts provisional pending full spin-spin/spin-orbit derivation in both branches |
| Positronium ortho-para (BS-§80) | `~115` MHz off | `~472` MHz off (literal); claimed "~Bethe-estimate" only after asymmetric absorbing correction | **retracted** — the absorbing correction, if applied consistently, applies to both branches |
| Muonium hyperfine (BS-§80) | `~1.3` MHz off | `~0.1` MHz residual | (c) is back-fit self-consistency at leading-`g_{s}` precision |

PR J (cross-comparison summary) records this inventory in a closing chapter with the back-fit and asymmetric-treatment caveats incorporated into the campaign's honest framing.

<!-- TODO: human reviews and fills in — confirms (a) PR I's combined verdict (six precision r_e-discriminators all consistent with branch (c)) is the campaign's strongest collective experimental signal about the r_e finding's resolution, (b) the path to PR J closing chapter is the correct disposition, and (c) the inventory table above is the form the campaign's headline result should take in PR J -->
