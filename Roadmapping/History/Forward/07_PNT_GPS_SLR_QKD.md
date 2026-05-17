---
chapter: 07
title: "Position-Navigation-Timing: Mercury Perihelion, GPS, SLR, QKD"
era: "forward"
variant: derivational
threads: [electromagnetism, quantum]
animations: [pnt_mercury_perihelion, pnt_gps_relativity, pnt_qkd_bb84]
verification_anchors: ["Two_Mathematically_Equivalent_Versions_of_Maxwells_Equations", "The_Classical_Electron_Problem", "Dual_Relativistic_Quantum_Mechanics_I"]
status: draft
---

# Chapter 7 — Position-Navigation-Timing: Mercury Perihelion, GPS, SLR, QKD

> **Framing principle (load-bearing).** We are exploring differences in mathematical conventions and their consequences for the physical interpretation of established experiments. We are not inventing new physics. This chapter is the **derivational** forward chapter — tag-dominant `#inferred` and `#verified`, *not* `#speculative`. The dual framework's verified content (Maxwell paper Eqs. 1–2, 9, 10–11; TCEP Eqs. 4.5, 4.16) bears directly on PNT through clock synchronization across observer-satellite pairs. **General relativity** — load-bearing for GPS — is `#gill-silent`; GR predictions are taken intact.

## 1. Overview

Position-Navigation-Timing (PNT) is the operational discipline of *knowing where you are, when, and how fast you're going* — at the precision modern satellite navigation, geodesy, and quantum communication require. It is a *consequence* of the 1800–1965 historical arc: Maxwell + atomic clocks + Einstein's special and general relativity + quantum mechanics, all engineered together. This chapter derives PNT basics from first principles — both standard SR + GR and Gill's dual framework — then walks four applications: Mercury's perihelion (the pedagogical bridge to GR clocks), GPS, Satellite Laser Ranging, and Quantum Key Distribution.

The chapter's structure follows [[_template_forward_derivational]] (Variant B): unlike [[08_quantum_computing_open_questions|Chapter 8]] and [[09_fusion_open_questions|Chapter 9]] which are *speculative roadmaps*, PNT is *derivationally heavy* — concrete equations with numerical predictions and explicit error bounds, both standard and dual-framework.

## 2. Historical roots

- **Doppler 1842 (Ch 2):** frequency shift from relative motion — basis of radial-velocity from frequency measurements.
- **[[verrier1859_mercury|Le Verrier 1859]] (Ch 2):** discovers Mercury's 43″/century perihelion anomaly that Newtonian gravity + planetary perturbations cannot explain.
- **[[maxwell1865_dynamical_theory|Maxwell 1865]] (Ch 2):** EM signal propagation; light as the metric of distance.
- **[[einstein1905_specrel|Einstein SR 1905]] (Ch 3):** time dilation, $v^2/2c^2$ leading-order correction.
- **[[sagnac1913_effect|Sagnac 1913]] (Ch 3):** rotation-induced fringe shift; basis for rotating-frame corrections.
- **[[einstein1915_perihelion|Einstein 1915 Mercury paper]] + [[einstein1916_grundlage|*Grundlage* 1916]] + [[schwarzschild1916|Schwarzschild 1916]] (Ch 3):** GR field equations and the Schwarzschild metric — derives Mercury's 43″/century *exactly*. Gravitational time dilation follows from the same metric.
- **Essen–Parry cesium clock 1955 (Ch 5):** atomic time standard. Makes the GR clock effects observable.

## 3. First-principles derivation: PNT basics

### 3.1 Distance from time of flight

Standard derivation. EM signal emitted at $t_e$ from a known transmitter location, received at $t_r$ from an unknown receiver location. The distance is $d = c(t_r - t_e)$, modulo atmospheric refraction and clock synchronisation errors. To localise a receiver to 1 m precision requires synchronisation to $\sim$3 ns precision (light travels 30 cm/ns). To localise to 1 cm requires $\sim$30 ps precision.

GPS uses four satellites to over-determine position: three for position + one for clock-offset solving. The receiver's clock need not be a high-precision atomic clock; the four-satellite intersection self-calibrates the receiver clock against the satellite atomic clocks.

**Dual-framework rederivation.** The natural launching point is [[Two_Mathematically_Equivalent_Versions_of_Maxwells_Equations]] Eq. (9): $t = (1/c) \int b(s)\, ds$ — the observer time corresponds to a path-integral of the collaborative speed $b = \sqrt{c^2 + u^2}$ along the source's worldline. For a satellite-borne atomic clock and a ground-based receiver, the two clocks correspond to two different proper-time parameterisations — the satellite at orbital velocity $u_{\rm sat} \approx 3.87$ km/s, the ground station at Earth-rotation velocity $u_{\rm gnd} \sim 0.5$ km/s. The proper-time-integrated time-of-flight gives the standard answer to the precision the dual framework can predict: deviation from standard SR at the satellite's $u/c \sim 1.3 \times 10^{-5}$ corresponds to a fractional time-error of order $u^2/(2c^2) \sim 8.5 \times 10^{-11}$, which is $\sim$10 ps over an integration time of 1 second — *below current GPS clock precision* but approaching the threshold where next-generation optical-lattice clocks would start to see a discrepancy.

### 3.2 Time-of-flight precision vs. clock precision

GPS atomic clocks (Cs and Rb beam standards) have fractional frequency stabilities of $\sim 10^{-13}$. Modern optical lattice clocks (Sr, Yb) reach $\sim 10^{-18}$. The dual-framework prediction for the SR piece of GPS-clock satellite-vs-ground rate difference: same as standard SR to within $\sim 10^{-10}$, well below 1980s GPS clock precision; *equal to* current GPS precision at the satellite-velocity scale; *distinguishable* by 2030s optical-lattice satellite clocks if they reach $10^{-19}$.

## 4. Applications walk-through

### §A — Mercury's perihelion (the pedagogical bridge to GPS)

[[einstein1915_perihelion|Einstein's 1915 derivation]] takes the Schwarzschild metric

$$ds^2 \;=\; \left(1 - \frac{2GM}{c^2 r}\right) c^2 dt^2 - \left(1 - \frac{2GM}{c^2 r}\right)^{-1} dr^2 - r^2 d\Omega^2,$$

derives the equatorial-plane geodesic equation, and shows that to leading post-Newtonian order the orbit precesses by $\Delta\phi = 6\pi GM/(c^2 a(1-e^2))$ per revolution. For Mercury — $a = 5.79 \times 10^{10}$ m, $e = 0.2056$ — this gives 43 arcseconds per century. Matches [[verrier1859_mercury|Le Verrier's residual]] exactly. **First experimental confirmation of GR.**

Animated in `pnt_mercury_perihelion.py`. The walkthrough: Schwarzschild metric → equatorial geodesic equation → effective radial potential $V_{\rm eff}(r) = -GM/r + L^2/(2r^2) - GM L^2/(c^2 r^3)$ → the extra $1/r^3$ term beyond Newton → integration over one orbit → 43″/century.

**Dual-framework status:** `#gill-silent`. Gill has not published an extension of the dual framework to GR. The Schwarzschild derivation is taken intact from standard GR. The chapter's role for §A is pedagogical — the same physics (clocks differing in gravitational potential accumulate different proper times) that explains Mercury's 43″/century also explains GPS's $+45\,\mu$s/day clock-rate difference.

### §B — GPS as the most precise everyday test of relativity

Standard derivation:
- **SR correction:** satellite at $v_{\rm sat} \approx 3.87$ km/s gives time-dilation rate change $-v^2/(2c^2) \approx -8.3 \times 10^{-11}$, accumulating to $-7\,\mu$s/day.
- **GR correction:** satellite at altitude $h \approx 20{,}200$ km has gravitational potential less negative by $\Delta U/c^2 \approx +5.3 \times 10^{-10}$, giving rate change $+\Delta U/c^2$, accumulating to $+45\,\mu$s/day.
- **Net:** $+38\,\mu$s/day. Without correcting for this, GPS position error grows by $\sim$10 km/day.

Engineering solution: satellite clocks are pre-detuned by exactly the predicted rate so that they tick at the GPS-Time rate after both relativistic corrections.

Animated in `pnt_gps_relativity.py`. The walkthrough: ground and satellite clocks → SR/GR rate differences → engineering pre-detuning.

**Dual-framework rederivation:** SR piece via [[Two_Mathematically_Equivalent_Versions_of_Maxwells_Equations]] Eq. (9) gives the same $-7\,\mu$s/day to within fractional $10^{-10}$ — below current GPS clock precision. GR piece is `#gill-silent`; taken intact. *Net dual prediction: same $+38\,\mu$s/day as standard*. The chapter is honest that current GPS doesn't distinguish standard from dual; the precision regime that would is optical-lattice clocks at $10^{-19}$, achievable in the 2030s.

### §C — Satellite Laser Ranging (SLR)

[[smith_christodoulidis1985_slr|Smith and Christodoulidis 1985]]. Ground stations fire short laser pulses at a satellite covered in corner-cube retroreflectors (LAGEOS-1, LAGEOS-2, ETALON, …). The round-trip time gives the satellite distance to millimetre precision. Used for geocenter monitoring, tectonic-plate motion measurement, and tests of GR including the Lense–Thirring frame-dragging effect (Gravity Probe B, 2011).

**Dual-framework status:** The Liénard–Wiechert geometry for the retarded light path is `#verified` from [[The_Classical_Electron_Problem]] Sec. 3.2; the group-velocity transformation enters via Eq. (4.16) (the chapter's Finding 3 sign-typo correction — addressed in [[FINDINGS_for_author_review]]). The dual-framework round-trip time matches standard SLR to within the same $10^{-10}$ regime as GPS. No measurable difference at current SLR precision.

### §D — Quantum Key Distribution (BB84)

[[bennett_brassard1984_bb84|Bennett & Brassard 1984]]. Alice sends polarised photons in random basis (rectilinear ↑→ or diagonal ↗↘) to Bob. Bob measures in random basis. They publicly compare bases (not measurement results); when bases agreed, the measurement bits form the shared key. Errors signal eavesdropping. Security follows from the **no-cloning theorem**: an eavesdropper cannot reproduce an unknown quantum state without introducing detectable errors.

Modern QKD demonstrations: [[liao2018_micius_satellite|Micius satellite (2017)]] achieves intercontinental QKD at ~kHz key rate over ~1200 km. Fibre-based QKD operates at MHz rates over <100 km. [[pirandola2020_qkd_review|Pirandola et al. 2020]] reviews the state of the field.

Animated in `pnt_qkd_bb84.py`. The walkthrough: Alice's random-basis encoding → channel → Bob's random-basis measurement → public basis comparison → key extraction → eavesdropper-detection via QBER threshold.

**Dual-framework status:** `#gill-silent` for the QKD protocol itself (single-photon states in standard QM Hilbert space; the no-cloning theorem is a structural property of the Hilbert-space inner product, common to standard and dual QM). [[08_quantum_computing_open_questions|Chapter 8 §5.2]] proposes a `#speculative` extension into Gill's KS-Hilbert space ([[FOUNDATIONS_FOR_QED_I_MATHEMATICAL]]) where the no-cloning argument could be re-examined. For fibre and satellite QKD as currently practised, dual and standard predictions are identical.

## 5. What proper-time changes (and what it doesn't)

| Application | Standard prediction | Dual prediction | Divergence | Current precision | Distinguishable? |
|---|---|---|---|---|---|
| Mercury perihelion | 43″/century (GR) | same (Gill-silent on GR) | — | $\sim 0.1''$/century | no |
| GPS SR clock rate | $-7\,\mu$s/day | $-7\,\mu$s/day + $O(10^{-10})$ | $\sim 1$ ps/day | $\sim 1$ ns/day | no at current GPS precision |
| GPS GR clock rate | $+45\,\mu$s/day | same (Gill-silent on GR) | — | $\sim 1$ ns/day | no |
| SLR round-trip | standard Liénard–Wiechert | dual Liénard–Wiechert at $u=0$ | identical | $\sim$1 mm | no |
| QKD BB84 | standard Hilbert space | dual Hilbert space at $u=0$ | identical | $\sim$1% QBER | no |

The pattern: dual-framework PNT corrections sit *one or two orders of magnitude below* current measurement precision in all cases. The regime where the divergence would matter is optical-lattice satellite clocks at $10^{-19}$, achievable in the 2030s; the SR contribution at $u/c \sim 10^{-5}$ is $\sim 10^{-10}$, just at the edge.

## 6. Bibliography

### Primary

- [[verrier1859_mercury]] — discovers the 43″/century anomaly.
- [[einstein1915_perihelion]] — derives it from GR.
- [[einstein1916_grundlage]] — systematic GR exposition.
- [[schwarzschild1916]] — vacuum spherical-symmetric solution.
- [[einstein1905_specrel]] — special relativity.
- [[sagnac1913_effect]] — rotating-frame fringe shift.
- [[vessot1980_gravity_probe_a]] — confirms gravitational time dilation to ~70 ppm.
- [[ashby2003_gps_relativity]] — definitive *Living Reviews* article on GPS relativity.
- [[parkinson1996_gps_theory]] — standard GPS engineering reference.
- [[smith_christodoulidis1985_slr]] — SLR overview.
- [[bennett_brassard1984_bb84]] — BB84 protocol.
- [[ekert1991_qkd]] — Bell-inequality QKD.
- [[liao2018_micius_satellite]] — first satellite QKD.

### Retrospective

- [[misner_thorne_wheeler1973]] — *Gravitation* textbook.
- [[wald1984_gr]] — Wald's GR textbook (companion to MTW).
- [[nielsen_chuang2000_quantum_computation]] — QKD textbook reference.
- [[pirandola2020_qkd_review]] — modern QKD review.
- [[jackson1998_classical_electrodynamics]] — classical-EM reference; cited for retarded potentials.

## 7. Cross-references

- Historical roots: [[01_early_electromagnetism_1800_1860]], [[02_classical_synthesis_1860_1900]], [[03_old_quantum_theory_1900_1925]], [[04_quantum_mechanics_1925_1948]], [[05_QED_renormalization_solid_state_1948_1965]].
- Verification anchors: [[Two_Mathematically_Equivalent_Versions_of_Maxwells_Equations]] (Eqs. 1–2 velocity duality, Eq. 9 t–τ integral), [[The_Classical_Electron_Problem]] (Eqs. 4.5 Doppler, 4.16 group velocity — note Finding 3 sign typo), [[Dual_Relativistic_Quantum_Mechanics_I]] (Eq. II.3 + KS-Hilbert space — QKD anchor).
- Companion forward chapters: [[08_quantum_computing_open_questions]], [[09_fusion_open_questions]].
- Findings: [[FINDINGS_for_author_review]] — Finding 3 (TCEP Eq. 4.16 sign typo) is relevant to §C SLR derivation.
