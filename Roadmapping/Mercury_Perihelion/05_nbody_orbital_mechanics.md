# 05 — Proper orbital mechanics + the N-body question: can the framework reach the full 43″?

**Status:** ✅ Wolfram-MCP verified (exact + numerical).
**Companion notebook:** [`../Mathematica_Notebooks/Mercury_Perihelion/perihelion_advance.wl`](../Mathematica_Notebooks/Mercury_Perihelion/perihelion_advance.wl) — 6 cells (Eq. h3, Eq. h4, exact −1/6 apsidal-angle ratio, Mercury numerics, independent orbit integration confirming retrograde sign, N-body scaling). Runs in Mathematica or `wolframscript`.
**Motivating question** (from the devil's advocate review of [PR #83](https://github.com/temoTxt/PyPhysics/pull/83)): the observed `43″/century` is an *N-body-derived* residual; the paper's two-body calculation isn't directly comparable. Does a proper N-body / eccentric-orbit treatment of the framework's force law close the 12% gap to observation?

**Answer: no.** Two independent reasons, both decisive:

1. The framework's force law gives a perihelion advance of **exactly `−1/6` of GR** (opposite sign, one-sixth magnitude) when computed properly for an *eccentric* orbit — `−7.17″/century` vs GR's `+42.99″`. This is a structural feature of the `(1 − GM/(c²r))` correction, not an artifact of the paper's circular-orbit approximation.
2. The N-body planetary perturbations contribute framework relativistic corrections `~10⁻⁴–10⁻⁵` of the Sun-Mercury term (because the relativistic factor `GM/(c²r)` scales with the *attracting* mass, and the planets are light). Even summing all planets gives `~0.01″/century` — three orders of magnitude too small to bridge the `−7 → +43` gap.

## 1. Why the paper's circular-orbit treatment is a heuristic

Perihelion advance is intrinsically about *eccentric* orbits — it is the precession of the orbital ellipse's major axis. **A circular orbit has no perihelion.** The paper's Eqs. (h5)–(h8) compute the dual *angular frequency* `ω_d` for a circular orbit at radius `r₀` and interpret `ω_d·T₀ − 2π` as a perihelion advance. That is a heuristic: it conflates a change in orbital *frequency* with a precession of the orbital *ellipse*. To compute the genuine perihelion advance, one must use the orbit's apsidal angle for the actual eccentric trajectory.

## 2. The framework's force law as a potential perturbation

The framework's relative two-body acceleration (paper Eq. h4, combined):

$$\mathbf{a} \;=\; -\frac{G(M+m)}{r^2}\!\left[1 - \frac{G(M^2+m^2)}{(M+m)c^2 r}\right]\hat{\mathbf{e}}_r \;=\; -\frac{G(M+m)}{r^2}\hat{\mathbf{e}}_r \;+\; \frac{G^2(M^2+m^2)}{c^2 r^3}\hat{\mathbf{e}}_r.$$

The correction term `+G²(M²+m²)/(c²r³)` is a *repulsive* `1/r³` force. The corresponding potential (per unit reduced mass):

$$\Phi(r) \;=\; -\frac{G(M+m)}{r} \;+\; \frac{B}{r^2}, \qquad B \;=\; \frac{G^2(M^2+m^2)}{2c^2} \;>\; 0.$$

The framework adds a **`1/r²` potential perturbation** (repulsive). This is structurally different from GR, whose Schwarzschild geodesic adds a **`1/r³` potential perturbation** (the attractive `−GMℓ²/(c²r³)` term in the effective potential).

## 3. Exact perihelion advance — the `1/r²` perturbation is exactly solvable

A `1/r²` potential perturbation merges with the centrifugal term: `ℓ²/(2r²) + B/r² = (ℓ² + 2B)/(2r²)`, i.e., the orbit is a Kepler problem with a *rescaled* angular momentum `ℓ' = √(ℓ² + 2B)`. The apsidal angle is then exactly `π·ℓ/ℓ'`, and the precession per orbit is:

$$\Delta\varphi_{\rm framework} \;=\; 2\pi\!\left(\frac{\ell}{\sqrt{\ell^2 + 2B}} - 1\right) \;\approx\; -\frac{2\pi B}{\ell^2}\quad(\text{small }B).$$

Substituting `B = G²(M²+m²)/(2c²)` and `ℓ² = G(M+m)a(1−e²)`:

$$\boxed{\;\Delta\varphi_{\rm framework} \;=\; -\frac{\pi G(M^2+m^2)}{c^2 a (1-e^2)(M+m)} \;\xrightarrow{m\ll M}\; -\frac{\pi GM}{c^2 a (1-e^2)}.\;}$$

Compare GR (MTW §40.5):

$$\Delta\varphi_{\rm GR} \;=\; +\frac{6\pi GM}{c^2 a (1-e^2)}.$$

**The ratio is exactly:**

$$\boxed{\;\frac{\Delta\varphi_{\rm framework}}{\Delta\varphi_{\rm GR}} \;=\; -\frac{1}{6}\cdot\frac{M^2+m^2}{M(M+m)} \;\xrightarrow{m\ll M}\; -\frac{1}{6}.\;}$$

The framework's perihelion advance is **−1/6 of GR's: opposite sign, one-sixth magnitude.**

### Wolfram MCP check

```wolfram
ClearAll[BB, ellSq, aa, ecc, GG, MM, mm, cc]; deltaExact = -2 Pi BB/ellSq; dF = deltaExact /. {BB -> GG^2 (MM^2 + mm^2)/(2 cc^2), ellSq -> GG (MM + mm) aa (1 - ecc^2)}; dGR = 6 Pi GG MM/(cc^2 aa (1 - ecc^2)); Print["Δφ_framework = ", FullSimplify[dF]]; Print["ratio framework/GR (full) = ", FullSimplify[dF/dGR]]; Print["ratio (m→0) = ", FullSimplify[dF/dGR /. mm -> 0]]
```

**Result:**
```
Δφ_framework = GG (mm² + MM²) Pi / (aa cc² (-1 + ecc²) (mm + MM))
ratio framework/GR (full) = -(1/6) (mm² + MM²)/(MM (mm + MM))
ratio (m→0) = -1/6
```
✅ exact `−1/6`.

## 4. Numerical values for Mercury

```wolfram
ClearAll[Msun, mMerc, aMerc, eMerc, Torb, GG, cc, orbitsPerCentury, radToArcsec]; Msun = 1.98892*^30; mMerc = 3.3011*^23; aMerc = 5.7909*^10; eMerc = 0.20563; Torb = 87.969*86400.; GG = 6.6743*^-11; cc = 299792458.; orbitsPerCentury = 100*365.25*86400./Torb; radToArcsec = 180/Pi*3600.; dGR = 6 Pi GG Msun/(cc^2 aMerc (1 - eMerc^2)); dFrame = -Pi GG (Msun^2 + mMerc^2)/(cc^2 aMerc (1 - eMerc^2) (Msun + mMerc)); Print["GR per century = ", N[dGR orbitsPerCentury radToArcsec, 6]]; Print["framework (eccentric) per century = ", N[dFrame orbitsPerCentury radToArcsec, 6]]; Print["ratio = ", N[dFrame/dGR, 8]]
```

**Result:**
```
GR per century = 42.9918 arcsec
framework (eccentric) per century = -7.1653 arcsec
ratio = -0.166667
```

| Quantity | Value |
|---|---|
| GR perihelion advance | `+42.99″/century` (matches observed `≈43″` to <0.1%) |
| Framework (proper, eccentric orbit) | **`−7.17″/century`** (`= −1/6 × GR`) |
| Framework (paper's circular approximation) | `−6.86″/century` (the `(1−e²)`-suppressed version; matches the paper's heuristic to within the eccentricity factor `1/(1−e²) = 1.044`) |

The paper's circular-orbit heuristic (`−6.86″`) and the proper eccentric-orbit calculation (`−7.17″`) agree on sign and magnitude; the proper calculation is the eccentricity-enhanced refinement.

## 5. The "1/6 factor" — a classic diagnostic

A perihelion advance of `1/6` the GR value is one of the most recognizable diagnostics in the history of gravitation. Pre-GR and alternative-gravity treatments that capture only *part* of the GR effect — for example, special-relativistic velocity-dependent-mass arguments, or simple potential modifications that lack the spatial-metric-curvature contribution — characteristically land at a fraction of the GR value, frequently `1/6`. Einstein's full `6/6` in 1915 came precisely from the spatial curvature that flat-space modifications do not reproduce.

The framework's modified-Newtonian gravity is exactly such a flat-space modification: it adjusts the *force law* (via the `V²/(2mc²)` kernel applied to `V = −GMm/r`) but introduces no spatial-metric curvature. It therefore lands at `−1/6` — the `1/6` fraction characteristic of force-law-only modifications, with the *negative* sign coming from the framework's correction being a repulsive `1/r²` potential perturbation (weaker gravity) rather than GR's attractive `1/r³` term.

**This connects directly to the Van Flandern reference attached to issue [#81](https://github.com/temoTxt/PyPhysics/issues/81).** Van Flandern analyzed how alternative gravity formulations produce fractions of the GR Mercury perihelion; the framework's `−1/6` places it squarely in the class of force-law modifications that Van Flandern's analysis covers. *(A full reading of the Van Flandern PDF is deferred to the bibliography-ingestion step of the [plan](../../.dev/tasks/81_Perhelion_of_Mercury.md); this cross-reference is provisional until that reading is done.)*

<!-- TODO: human reviews and fills in — confirms (a) the −1/6 ratio is the correct exact result for the framework's force law; (b) the "1/6 is a classic diagnostic of force-law-only modifications" framing is accurate and not overstated; (c) the Van Flandern connection is correctly characterized (provisional pending a reading of the actual PDF). -->

## 6. The N-body question, quantified

The 43″ residual is computed by subtracting Newtonian planetary perturbations (~531″/century from Venus, Earth, Jupiter, etc. on Mercury) from the total observed precession (~5600″/century). Would the framework's modified gravity, applied consistently to *all* planet-on-Mercury forces, change the residual enough to matter?

**No.** The framework's relativistic correction for a body of mass `M_body` at distance `r` carries the factor `GM_body/(c²r)`. This scales with the *attracting* mass:

```wolfram
ClearAll[GG, cc, Msun, aMerc, Mjup, aJup]; GG = 6.6743*^-11; cc = 299792458.; Msun = 1.98892*^30; aMerc = 5.7909*^10; Mjup = 1.898*^27; aJup = 7.785*^11; Print["Sun-Mercury factor GM/(c²r) = ", ScientificForm[N[GG Msun/(cc^2 aMerc), 4]]]; Print["Jupiter-Mercury factor = ", ScientificForm[N[GG Mjup/(cc^2 aJup), 4]]]; Print["ratio Jupiter/Sun = ", ScientificForm[N[(Mjup/aJup)/(Msun/aMerc), 4]]]
```

**Result:**
```
Sun-Mercury factor GM/(c²r) = 2.551×10⁻⁸
Jupiter-Mercury factor = 1.811×10⁻¹²
ratio Jupiter/Sun = 7.099×10⁻⁵
```

Jupiter — the largest planetary perturber — contributes a framework relativistic correction `~7×10⁻⁵` of the Sun-Mercury term. Summing all planets gives at most `~10⁻⁴` of the Sun-Mercury `−7.17″`, i.e., `~10⁻³ arcsec/century`. **The N-body planetary corrections are three to four orders of magnitude too small to bridge the `−7 → +43` gap.** The dominant framework relativistic effect is, and remains, the Sun-Mercury two-body term — which gives `−7.17″/century`.

## 7. Verdict

❌ **The N-body / proper-orbital-mechanics formulation does not reach the observed `+43″/century`.** The framework's force law gives a perihelion advance of exactly `−1/6` of GR (`−7.17″/century`: opposite sign, one-sixth magnitude), and the N-body planetary corrections are `~10⁻³`–`10⁻⁴` too small to matter. To reach the full measurement, the framework would need a *fundamentally different* gravitational correction — specifically the spatial-metric-curvature contribution (the `1/r³` effective-potential term) that GR has and the flat-space `V²/(2mc²)` modification does not produce. This is a clean, decisive answer to the open caveat from the [PR #83 devil's advocate review](https://github.com/temoTxt/PyPhysics/pull/83), and it sharpens the GPS author report's Q1: the flat-space modified-Newtonian extension is ruled out for Mercury *not* by a 12% numerical gap but by the structural `−1/6` factor, which no N-body refinement can close.

<!-- TODO: human reviews and fills in — confirms the §7 verdict is the right honest conclusion: the N-body exploration does not save the framework's gravitational extension, and the −1/6 structural result (rather than the earlier 12% numerical gap) is the load-bearing finding. Confirms this should be promoted into FINDINGS_for_author_review.md as the substance of the Mercury finding. -->
