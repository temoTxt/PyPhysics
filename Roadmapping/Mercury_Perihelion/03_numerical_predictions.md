# 03 — Numerical predictions: Mercury perihelion under three framework forms + GR + observed

**Status:** ✅ all Wolfram-MCP computed.
**Source paper:** [`../Tepper_Gill_Papers/Dual Newtonian Theory.tex`](../Tepper_Gill_Papers/Dual%20Newtonian%20Theory.tex), §2 *The Sun-Mercury System* (lines 309–338 list the three predictions).

## 1. Parameters used

Per [README §Conventions](README.md):

```
M_⊙ = 1.98892 × 10³⁰ kg
m_Mercury = 3.3011 × 10²³ kg
a_Mercury = 5.7909 × 10¹⁰ m
e_Mercury = 0.20563
T_orb = 87.969 days
G = 6.6743 × 10⁻¹¹ m³ kg⁻¹ s⁻²
c = 299 792 458 m/s
```

Derived: `m/M = 1.6597 × 10⁻⁷`; `GM/(c²a) = 2.55 × 10⁻⁸`; orbits per tropical century = `415.20`.

## 2. The three paper predictions and standard GR

The paper offers three predictions (lines 309–333 of the .tex):

**(1) Corda's value** — uses only the reduced-mass factor:

$$\omega_c \;=\; \omega_0\,(1 + m/M)^{1/2} \;\approx\; \omega_0\,(1 + m/(2M)).$$

$$\Delta\varphi_c \;=\; 2\pi\,(1 + m/(2M)) - 2\pi \;=\; \frac{\pi m}{M}.$$

**(2) Full dual** — the framework's structural prediction including the relativistic factor:

$$\Delta\varphi_d \;=\; 2\pi\left[(1 + m/M)\!\left(1 - \tfrac{G(M^2 + m^2)}{(M+m)c^2 r_0}\right)\right]^{1/2} - 2\pi.$$

**(3) Approximate dual** — first-order Taylor expansion of (2):

$$\Delta\varphi_{d_1} \;\approx\; 2\pi\!\left[\frac{m}{2M} - \frac{GM(1+m/M)}{2c^2 r_0}\right].$$

**Standard GR reference** (MTW §40.5 / Schwarzschild perihelion):

$$\Delta\varphi_{\rm GR} \;=\; \frac{6\pi GM}{c^2 a (1 - e^2)}\quad\text{per orbit}.$$

## 3. Wolfram MCP computation

```wolfram
ClearAll[Msun, mMerc, aMerc, eMerc, Torb, GG, cc, orbitsPerCentury, radToArcsec]; Msun = 1.98892*^30; mMerc = 3.3011*^23; aMerc = 5.7909*^10; eMerc = 0.20563; Torb = 87.969 * 86400.; GG = 6.6743*^-11; cc = 299792458.; orbitsPerCentury = 100*365.25*86400./Torb; radToArcsec = 180/Pi*3600.; (* Corda *) dphiCorda = N[Pi mMerc/Msun]; Print["Δφ_c per orbit rad = ", dphiCorda]; Print["Δφ_c per century arcsec = ", dphiCorda * orbitsPerCentury * radToArcsec]; (* Full dual *) factor = (1 + mMerc/Msun) (1 - GG (Msun^2 + mMerc^2)/((Msun + mMerc) cc^2 aMerc)); dphiDual = 2 Pi (Sqrt[factor] - 1); Print["Δφ_d per orbit rad = ", N[dphiDual, 12]]; Print["Δφ_d per century arcsec = ", N[dphiDual * orbitsPerCentury * radToArcsec, 8]]; (* Approximate dual *) dphiD1 = 2 Pi (mMerc/(2 Msun) - GG Msun (1 + mMerc/Msun)/(2 cc^2 aMerc)); Print["Δφ_d1 per orbit rad = ", N[dphiD1, 12]]; Print["Δφ_d1 per century arcsec = ", N[dphiD1 * orbitsPerCentury * radToArcsec, 8]]; (* Standard GR *) dphiGR = 6 Pi GG Msun/(cc^2 aMerc (1 - eMerc^2)); Print["Δφ_GR per orbit rad = ", N[dphiGR, 12]]; Print["Δφ_GR per century arcsec = ", N[dphiGR * orbitsPerCentury * radToArcsec, 8]]
```

**Result:**

```
Δφ_c per orbit rad = 5.21424×10⁻⁷
Δφ_c per century arcsec = 44.6557
Δφ_d per orbit rad = 4.41296×10⁻⁷
Δφ_d per century arcsec = 37.7934
Δφ_d1 per orbit rad = 4.41296×10⁻⁷
Δφ_d1 per century arcsec = 37.7934
Δφ_GR per orbit rad = 5.01995×10⁻⁷
Δφ_GR per century arcsec = 42.9918
```

## 4. Summary table

| Prediction | Δφ per orbit (rad) | Δφ per century (arcsec) | vs observed 43″ | vs GR 42.99″ |
|---|---|---|---|---|
| Corda `πm/M` | `5.21 × 10⁻⁷` | **44.66″** | `+3.9%` | `+3.9%` |
| Full dual `Δφ_d` | `4.41 × 10⁻⁷` | **37.79″** | **`−12.1%`** | **`−12.1%`** |
| Approximate dual `Δφ_{d₁}` | `4.41 × 10⁻⁷` | **37.79″** | **`−12.1%`** | **`−12.1%`** |
| Standard GR | `5.02 × 10⁻⁷` | **42.99″** | `+0.0%` | (reference) |
| **Observed** (Le Verrier residual) | — | **`≈ 43″`** | (reference) | `+0.02%` |

Standard GR reproduces observed to `≈ 0.05%` — the canonical Einstein 1915 result.

## 5. Decomposition of the framework's full-dual prediction

$$\Delta\varphi_d \;\approx\; \underbrace{\pi m/M}_{\rm reduced\text{-}mass\ (Corda)} \;+\; \underbrace{(-\pi GM/(c^2 r_0))}_{\rm framework\ relativistic\ correction}$$

For Mercury:

| Contribution | Δφ per century (arcsec) |
|---|---|
| Reduced-mass (Corda's `πm/M`) | `+44.66″` |
| Framework's relativistic correction (`−πGM/(c²r₀)`) | **`−6.87″`** |
| **Net framework prediction `Δφ_d`** | **`+37.79″`** |

**The framework's relativistic correction has the *opposite sign* from GR's `+42.99″/century`** and is about ~6× smaller in magnitude.

This is the load-bearing numerical finding of the campaign. The paper's headline claim that "44.39″/century well-approximates 42.98″ and observed 43″" is true only for the Corda value, which is the *reduced-mass-only* Newtonian effect — **not** the framework's structural prediction including the `(1 − GM/(c²r))` correction from Eq. (h4).

## 6. Discrepancy with the paper's quoted Corda value

The paper quotes `Δφ_c = 44.39″/century`; our computation gives `44.66″`. The `0.6%` discrepancy is consistent with the paper's likely use of slightly different orbital parameters (`m/M ratio`, sidereal vs tropical century, mean vs semi-major orbital radius). Both values exceed observed `43″` by `3-4%` and both exceed GR by the same.

## 7. Honest scoping notes

- **The Corda value** is mathematically a Newtonian reduced-mass effect (the orbital angular frequency scales as `√(G(M+m)/r³) = √(GM/r³) · √(1+m/M)`). It is not a relativistic prediction. The 19th-century Le Verrier analysis of Mercury's perihelion already accounted for reduced-mass effects in the planetary perturbation calculations; the `43″/century` residual is what remains *after* such Newtonian corrections.
- **The framework's structural prediction** (full dual or approximate dual, identical to leading order in this regime) is `37.79″/century`. This under-predicts observation by `12%` and disagrees with the sign convention of GR's relativistic correction (the framework's correction is *negative*, reducing the reduced-mass-driven advance; GR's is *positive*).
- The honest reading: the paper's modified-Newtonian extension does **not** reproduce observed Mercury perihelion at the precision the observation supports (`~10⁻⁴` relative). What "matches" is the reduced-mass-only Corda value, which is not the framework's relativistic prediction.

## 8. Verdict

✅ Numerical computation reproduces the paper's claimed values for all three predictions (Corda value within `0.6%`; full and approximate dual exact to printed precision). ⚠ The framework's full-dual prediction (`37.79″/century`) under-predicts observed Mercury perihelion advance by `12%` and has the opposite sign of relativistic correction compared to GR. Headline finding flagged for the GPS author report's open question Q1 — see [`04_findings_and_GPS_Q1.md`](04_findings_and_GPS_Q1.md).
