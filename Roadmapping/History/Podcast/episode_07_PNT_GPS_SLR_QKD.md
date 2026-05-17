---
episode: 07
title: "Position-Navigation-Timing: Mercury, GPS, SLR, QKD"
era: "forward"
chapter: 07_PNT_GPS_SLR_QKD
speakers: [Historian, Physicist, Experimentalist]
target_runtime_min: 12
animations_cued:
  - pnt_mercury_perihelion
  - pnt_gps_relativity
  - pnt_qkd_bb84
status: draft
---

# Episode 7 — Position-Navigation-Timing

> Companion dialogue script for [[07_PNT_GPS_SLR_QKD]]. First of the three forward chapters; the *derivational* one with concrete numerical content.

## Cold open

**Historian:** You unlock your phone, open a map app, and wait three seconds. The dot appears on the screen — your location to ten metres. Behind that three seconds: four satellites orbiting Earth at 20,200 kilometres up, broadcasting precisely-timed pseudorandom codes that travel through the ionosphere at roughly the speed of light; an atomic clock on each satellite that is pre-detuned by exactly 38 microseconds per day to account for general-relativistic and special-relativistic time-dilation effects; and a receiver in your phone that solves a four-equation, four-unknown system for your position plus its own clock offset, in software, sixty times per second.

**Experimentalist:** Take away the GR correction and your position drifts by ten kilometres a day. Take away the SR correction and you still drift, less far but still unusable. This is the most precise everyday test of relativity ever built. Hundreds of millions of phones simultaneously.

**Physicist:** And in this chapter — Chapter 7 of the history-of-physics campaign — we work through PNT, position-navigation-timing, as a *derivational* application of the dual-theory framework we've been tracking. Four applications: Mercury perihelion as the pedagogical bridge to GR clocks, GPS, satellite laser ranging, and quantum key distribution.

**Experimentalist:** And the dual framework's bearing on PNT?

**Physicist:** Three pieces. First: the SR contribution to GPS clock rates — that $-7$ microseconds per day — is rederivable in the dual framework via Gill's Maxwell paper equation 9, the $t = (1/c) \int b(s) ds$ relation between observer time and a path-integral of the collaborative speed along the source worldline. Same numerical answer as standard SR to within $10^{-10}$, well below current GPS clock precision but approaching optical-lattice satellite-clock precision in the 2030s. Second: the GR contribution — the $+45$ microseconds per day — is `#gill-silent`. Gill hasn't published an extension to GR. We take it intact from standard GR. Third: QKD — the BB84 protocol — works in both standard and dual Hilbert space because the no-cloning argument is structural. Same QBER predictions.

**Historian:** And the chapter's pedagogical bridge from history to engineering: Mercury's perihelion. The 1859 Le Verrier residual of 43 arcseconds per century is the same physics — clocks differing in gravitational potential — that drives the GPS clock correction.

**Physicist:** Let's start with Mercury.

## Historical sweep + first derivation

**Historian:** 1859. Urbain Le Verrier in Paris, the same astronomer who in 1846 had used Newtonian perturbation theory to predict Neptune from anomalies in Uranus's orbit, examines the data on Mercury. He finds that Mercury's perihelion — the closest point of its orbit to the Sun — advances by an amount that Newtonian gravity plus the perturbations of all other known planets cannot account for. The residual is 43 arcseconds per century. Small but unambiguous. Le Verrier speculates about an unseen planet inside Mercury's orbit, provisionally named Vulcan. None is found.

**Experimentalist:** The residual sits as an open astronomical puzzle for 56 years.

**Historian:** And in November 1915, in Berlin, Einstein completes the general relativity field equations and applies them to Mercury's orbit. The paper, [[einstein1915_perihelion]], is *Erklärung der Perihelbewegung des Merkur* — explanation of Mercury's perihelion motion. The Schwarzschild metric, derived by [[schwarzschild1916|Karl Schwarzschild]] within weeks while serving on the Russian front, gives the spacetime around a spherically symmetric mass. Apply the geodesic equation in the equatorial plane. Get an effective radial potential that's identical to Newton plus an extra $-GM L^2/(c^2 r^3)$ term.

`[cue: animation pnt_mercury_perihelion]`

**Physicist:** The animation walks through the calculation. Schwarzschild metric → equatorial geodesic equation → effective radial potential → integration over one orbit → perihelion shift per revolution $\Delta\phi = 6\pi G M / (c^2 a (1 - e^2))$. Plug in Mercury's semi-major axis $a = 5.79 \times 10^{10}$ meters and eccentricity $e = 0.2056$ and the gravitational parameter $GM$ of the Sun. Out comes 43 arcseconds per century. Matches Le Verrier's residual exactly.

**Experimentalist:** First experimental confirmation of GR.

**Historian:** And the same physics — gravitational time dilation, clocks deeper in a gravitational potential tick slower — drives GPS. The transition from Mercury's orbit to a GPS satellite is just a question of swapping the masses and altitudes.

## GPS as the most precise everyday test of relativity

**Physicist:** GPS satellites orbit at 20,200 kilometres altitude, with orbital velocity 3.87 kilometres per second. Two relativistic clock effects.

`[cue: animation pnt_gps_relativity]`

**Physicist:** The SR correction first. The satellite is moving at $v \approx 3.87$ km/s relative to the ground; time dilation makes its clock tick slower by a fractional rate $v^2/(2c^2) \approx 8.3 \times 10^{-11}$. Over a day, that accumulates to about $-7$ microseconds.

**Experimentalist:** And the GR correction?

**Physicist:** The satellite is *higher* in the Earth's gravitational potential than the ground station. Clocks at higher gravitational potential tick *faster*. The rate-increase is $\Delta U/c^2 \approx +5.3 \times 10^{-10}$ at 20,200 km altitude — about six times larger in magnitude than the SR effect, and with the opposite sign. Over a day, accumulates to $+45$ microseconds.

**Experimentalist:** Net?

**Physicist:** $+38$ microseconds per day. If you don't correct for it, your position error grows by *ten kilometres per day*. GPS is unusable without relativistic corrections.

**Historian:** And the engineering solution?

**Physicist:** Brilliant in its simplicity. The satellite atomic clocks are pre-detuned at manufacture by exactly the predicted rate. They run at the GPS-Time rate *after* both relativistic corrections. The receiver doesn't have to apply any relativistic correction; the clocks are calibrated to compensate from the moment they're switched on in orbit.

**Experimentalist:** Confirmed experimentally to high precision by [[vessot1980_gravity_probe_a|Vessot's 1976 Gravity Probe A]] — a hydrogen-maser-on-a-rocket sub-orbital flight that measured the predicted altitude-vs-time-dilation curve to ~70 ppm. Also verified continuously by the GPS system itself: any deviation from the predicted relativistic offset would show up as a system-wide ranging error.

**Physicist:** And the dual-framework rederivation. The SR piece — the $-7$ microseconds per day — falls out of [[Two_Mathematically_Equivalent_Versions_of_Maxwells_Equations]] equation 9: the observer time corresponds to a path-integral of the collaborative speed $b = \sqrt{c^2 + u^2}$ along the satellite's worldline. At $u/c \sim 1.3 \times 10^{-5}$ for the satellite, $b/c - 1 \sim 8.5 \times 10^{-11}$. Multiply by 86,400 seconds per day, get $\sim 7$ microseconds. Same answer as standard SR to within fractional $10^{-10}$.

**Experimentalist:** And the GR piece?

**Physicist:** `#gill-silent`. Gill hasn't published an extension to GR. We take the standard $+45$ microseconds per day intact. The dual framework's net prediction for GPS is *the same $+38$ microseconds per day as standard*.

**Experimentalist:** So at current GPS precision, the standard and dual framings are indistinguishable.

**Physicist:** Right. The regime where the SR-piece divergence would become measurable is *optical-lattice satellite clocks at $10^{-19}$ fractional stability* — achievable in the 2030s. Until then, GPS is a beautiful demonstration of GR + SR working together, with the dual-framework's reframing of the SR piece being conceptually interesting but predictively equivalent.

## Satellite Laser Ranging

**Historian:** Satellite laser ranging. Ground stations fire short laser pulses — picosecond duration — at satellites covered in corner-cube retroreflectors. The pulses return with a measurable round-trip time. Distance from $d = c(t_r - t_e)/2$. Precision: millimetres.

**Experimentalist:** Used for geocenter monitoring (where exactly is the center of the Earth's mass?), tectonic-plate motion measurement (centimetres per year), tides, polar motion, and tests of GR — including the Lense–Thirring frame-dragging effect detected by Gravity Probe B in 2011.

**Physicist:** The Liénard–Wiechert geometry for the retarded light path — the standard classical-EM tool for tracking a moving source-and-receiver pair — was verified in the original equation-verification campaign as part of [[The_Classical_Electron_Problem]] section 3.2. The group-velocity relation that enters at the photon-energy level is the TCEP equation 4.16 that's flagged for author review as Finding 3 — a sign typo in the paper, with the algebra and the paper's own commentary agreeing on the correct sign. SLR predictions match between standard and dual frameworks at SLR precision.

## Quantum Key Distribution

**Historian:** And QKD. [[bennett_brassard1984_bb84|Bennett and Brassard's 1984 BB84 protocol]]. Alice sends single photons polarised in random basis — rectilinear up-down or diagonal — to Bob. Bob measures each photon in random basis. They publicly compare bases, not measurement results. When the bases agreed, the measurement bits form a shared key.

`[cue: animation pnt_qkd_bb84]`

**Experimentalist:** The animation walks through the protocol. Alice's random encoding, Bob's random measurement, basis comparison, key extraction. Plus the eavesdropper-detection step: any in-channel measurement disturbs the quantum state, introducing detectable errors at a quantum bit error rate above some threshold.

**Physicist:** Security follows from the no-cloning theorem — Wootters and Zurek 1982, also Dieks 1982. An eavesdropper cannot make an arbitrary copy of an unknown quantum state. The dual framework's structural properties of Hilbert space — inner product, norm, completeness — are the same as standard QM's. The no-cloning argument is structural; the BB84 protocol works identically in both framings.

**Experimentalist:** Modern QKD: [[liao2018_micius_satellite|the Chinese Micius satellite (2017)]] achieves intercontinental QKD at kilohertz rates over 1200-kilometre satellite-to-ground links. Fibre QKD operates at megahertz rates over hundred-kilometre links. [[pirandola2020_qkd_review|Pirandola et al. 2020 review]] covers the state of the field.

**Physicist:** And the speculative bridge to [[08_quantum_computing_open_questions|Chapter 8]]: Gill's [[FOUNDATIONS_FOR_QED_I_MATHEMATICAL|Foundations I Mathematical]] paper organises QFT on a Kuelbs–Steadman Hilbert space. The no-cloning argument in *that* representation could in principle be re-examined. Whether the security guarantee changes — *tagged `#speculative`* in Chapter 8 — is a question deferred.

## Closing

**Historian:** Chapter 7 in one episode. Mercury perihelion as the GR pedagogical bridge. GPS as the most precise everyday test of relativity. SLR for geodesy and GR tests. QKD for cryptographic key distribution. All four sit on the 1800–1965 historical arc, all four have engineering precision now that lets us tell the difference between standard and dual-framework predictions in some cases — though for GPS the relevant regime is a generation in the future.

**Experimentalist:** The dual framework's contribution to PNT is `#inferred` on the SR pieces — same predictions as standard to within $10^{-10}$ — and `#gill-silent` on the GR pieces. Net: identical engineering predictions at current precision; structurally interesting reframing of the SR-piece derivation.

**Physicist:** Next episode is Chapter 8 — quantum-computing open questions, the first of the two *speculative* forward chapters. Roadmap rather than derivation. Honest about what we don't know. Thanks for listening.

`[cue: end card with bibliography wikilinks for show notes]`

## Show notes

Auto-generated from the chapter's bibliography by `lint_episode.py`.
