---
chapter: 03
title: "Old Quantum Theory 1900-1925"
era: "1900-1925"
threads: [quantum, electromagnetism, solid-state]
animations: [hist_bohr_proper_time, hist_compton_null]
verification_anchors: ["FoundationsII-Classical", "Analytic_Representation_of_The_Square-Root_Operator"]
status: draft
---

# Chapter 3 — Old Quantum Theory (1900–1925)

> **Framing principle (load-bearing).** We are exploring differences in mathematical conventions and their consequences for the physical interpretation of established experiments. We are not inventing new physics. The dual framework's rebranded Bohr atom and the dual framework's matter-wave dispersion relation reproduce all old-quantum-theory predictions at the precision available in 1900–1925 (and beyond). The chapter's interpretive payoff is conceptual — the positive-definite dual Hamiltonian $K = H^2/(2mc^2) + mc^2/2$ from [[FoundationsII-Classical]] gives the same Bohr levels as standard non-relativistic QM without ever needing a separate "kinetic energy" Hamiltonian piece.

## 1. Overview

The quarter-century from Planck's December 1900 blackbody paper to de Broglie's 1924 thesis is the era of *old quantum theory*: ad hoc quantisation rules added to classical mechanics, working brilliantly for some phenomena (hydrogen spectrum, photoelectric effect, fine structure) and failing decisively for others (helium spectrum, anomalous Zeeman effect, molecular spectra). Special relativity ([[einstein1905_specrel]]) and Minkowski's geometric reformulation ([[minkowski1908_raum_zeit]]) settle the kinematics of moving frames. General relativity ([[einstein1916_grundlage]]) completes itself in 1915 and explains Mercury's perihelion advance ([[einstein1915_perihelion]]) — both load-bearing for [[07_PNT_GPS_SLR_QKD|Chapter 7's]] GPS treatment. Rutherford's nuclear atom ([[rutherford1911_alpha_scattering]]) sets up Bohr's quantisation ([[bohr1913_atom_constitution]]); Sommerfeld's relativistic extension ([[sommerfeld1916_atombau]]) gives the fine-structure constant its first appearance. Compton ([[compton1923_xray_scattering]]) confirms photon momentum; de Broglie ([[de_broglie1924_thesis]]) proposes matter waves. The solid-state thread is opened by [[onnes1908_helium_liquefaction|Onnes's helium liquefaction]] and [[onnes1911_superconductivity|superconductivity discovery]] — both unexplained until the [[05_QED_renormalization_solid_state_1948_1965|Chapter 5]] BCS theory of 1957. Next: [[04_quantum_mechanics_1925_1948]] picks up with matrix and wave mechanics.

## 2. Historical narrative

### 2.1 1900 — Planck and the blackbody crisis

The 1890s ultraviolet catastrophe: Rayleigh and Jeans's classical derivation of blackbody spectral density gives $u(\nu, T) \propto \nu^2 k_B T$, diverging at high frequencies. Wien's empirical high-frequency formula $u \propto \nu^3 \exp(-h\nu/k_B T)$ fits the data above the peak but not below. Lummer and Pringsheim's careful 1899–1900 measurements at the Reichsanstalt show neither limit fits across the spectrum.

Planck's December 1900 paper [[planck1900_blackbody]] introduces the *ansatz* that the energy of each blackbody-cavity oscillator is constrained to integer multiples of $h\nu$: $E_n = nh\nu$. With this hypothesis, the Boltzmann combinatorial argument gives

$$u(\nu, T) = \frac{8\pi h \nu^3}{c^3} \cdot \frac{1}{e^{h\nu/k_B T} - 1}$$

— matching the experimental curve across the full spectrum. Planck treats $h$ as a calculational artefact ([[kuhn1978_blackbody|Kuhn (1978)]] argues the genuine physical-quantum interpretation arrives only after 1906). The numerical value of $h$ comes out at $\approx 6.55 \times 10^{-34}$ J⋅s; modern value $6.626 \times 10^{-34}$.

### 2.2 1905 — Einstein's *annus mirabilis*

In a single year Einstein publishes five papers, four of them load-bearing:

- **[[einstein1905_photoelectric]] — the light quantum.** Argues from a statistical analogy with the Wien blackbody regime that light *itself* is composed of independent quanta $h\nu$. The photoelectric effect's energy equation $E_k = h\nu - W$ follows immediately. Verified by [[millikan1916_photoelectric_verification|Millikan]] in 1916; Nobel Prize in 1921.
- **[[einstein1905_brownian]] — Brownian motion.** Predicts mean-square displacement $\langle x^2 \rangle = 2Dt$ with $D = k_B T / (6\pi \eta r)$. Verified by Perrin in 1908; decisive evidence for atomism.
- **[[einstein1905_specrel]] — special relativity.** Reinterprets the Lorentz transformations of [[lorentz1892_electron_theory]] as kinematic, derives length contraction and time dilation as coordinate effects, and resolves [[michelson_morley1887_relative_motion|Michelson–Morley]] kinematically. The dual-framework alternative (via velocity duality $\mathbf{w}/c = \mathbf{u}/b$) sits at [[Two_Mathematically_Equivalent_Versions_of_Maxwells_Equations]].
- **[[einstein1905_emc2]] — mass-energy.** Three-page sequel deriving $E = mc^2$. Foundational for nuclear physics and the rest-energy term $mc^2/2$ in the dual Hamiltonian.

### 2.3 1908 — Minkowski's geometric reformulation

[[minkowski1908_raum_zeit|Minkowski's "Raum und Zeit" lecture]] in Cologne recasts Einstein 1905 in 4-dimensional geometric language: spacetime as a pseudo-Riemannian manifold with metric $\text{diag}(+,-,-,-)$; worldlines as 4-curves; the Lorentz invariant $ds^2 = c^2 dt^2 - dx^2 - dy^2 - dz^2$. "Henceforth space by itself, and time by itself, are doomed to fade away into mere shadows."

The Minkowski formulation is decisive for the standard development. It's also the entry point for [[FOUNDATIONS_FOR_QED_I_MATHEMATICAL]]'s **Minkowski Incompatible Theorem 1.1**, which shows that the standard 4-velocity $u^\mu = dx^\mu/d\tau$ does *not* close on itself under Hamiltonian flow without an additional structural assumption — motivating Gill's positive-definite dual Hamiltonian via squaring.

### 2.4 1911–1916 — Rutherford, Bohr, Sommerfeld, GR

**[[rutherford1911_alpha_scattering|Rutherford's nuclear atom (1911)]].** The Geiger–Marsden gold-foil scattering experiments show a small fraction of alpha particles back-scattered — incompatible with Thomson's "plum pudding" but consistent with a small, dense, positively charged nucleus surrounded by orbiting electrons. The atomic structure problem now has a quantitative substrate.

But: classical electrodynamics demands that orbiting (i.e., accelerating) electrons radiate, spiraling into the nucleus on a microsecond timescale. Atoms shouldn't exist as stable structures. The contradiction motivates Bohr.

**[[bohr1913_atom_constitution|Bohr's 1913 atomic model]].** Two postulates:
1. Electrons occupy discrete *stationary* circular orbits with angular momentum $L = n\hbar$, in which they do not radiate.
2. Transitions between orbits emit/absorb photons with $h\nu = E_m - E_n$.

The hydrogen spectrum follows: $E_n = -\frac{1}{2}\frac{e^2}{4\pi\varepsilon_0 a_0}\frac{1}{n^2}$ with $a_0 = \frac{4\pi\varepsilon_0 \hbar^2}{m_e e^2}$, reproducing the Rydberg formula $\nu = R(1/n_m^2 - 1/n_n^2)$ exactly, and predicting $R = m_e e^4 / (8 \varepsilon_0^2 h^3)$ in agreement with the measured value.

**[[sommerfeld1916_atombau|Sommerfeld 1916]]** extends Bohr's circular orbits to elliptical with a special-relativistic momentum correction, deriving the *fine structure* of hydrogen as

$$E_{n,k} = -\frac{R}{n^2}\left[1 + \frac{\alpha^2}{n}\left(\frac{1}{k} - \frac{3}{4n}\right) + O(\alpha^4)\right]$$

with the fine-structure constant $\alpha = e^2/(\hbar c) \approx 1/137$ making its first appearance in spectroscopy. The same $\alpha$ recurs in the [[Dual_Relativistic_Quantum_Mechanics_I|DRQM I]] dual-Dirac fine-structure derivation (Chapter 4).

**[[einstein1915_perihelion|Einstein's Mercury perihelion paper (1915)]] and [[einstein1916_grundlage|the GR foundational paper (1916)]].** Field equations $G_{\mu\nu} = 8\pi G T_{\mu\nu}/c^4$ couple matter-energy to spacetime curvature; the Schwarzschild solution applied to Mercury's orbit gives $\Delta\phi = 6\pi GM/(c^2 a(1-e^2)) \approx 43''/\text{century}$, matching [[verrier1859_mercury|Le Verrier's 1859 residual]] exactly. The first observational confirmation of GR. Re-derived step-by-step in [[07_PNT_GPS_SLR_QKD|Chapter 7 §B]] as the pedagogical bridge to GPS clock corrections.

### 2.5 1923–1924 — Compton and de Broglie

**[[compton1923_xray_scattering|Compton's X-ray scattering experiment]]** measures the wavelength shift of X-rays scattered off (effectively) free electrons at varying angles. The observed shift $\Delta\lambda = (h/m_e c)(1 - \cos\theta)$ matches the prediction from a photon-electron elastic collision treating photons with momentum $p = h/\lambda$. Direct experimental confirmation of Einstein's 1916 conjecture that photons carry momentum proportional to their wavenumber. The Bohr–Kramers–Slater attempt (1924) to keep electromagnetism continuous *and* explain photoelectric/Compton effects via statistical light-emission rules fails when Bothe and Geiger (1925) show photon-electron momentum is conserved event-by-event, not statistically. Photons are real particles.

**[[de_broglie1924_thesis|De Broglie's matter-wave hypothesis]].** If photons (waves) carry momentum $p = h/\lambda$, then perhaps electrons (particles) carry wavelength $\lambda = h/p$. The thesis is mathematically slight but conceptually decisive. Davisson–Germer (1927) and G. P. Thomson (1928) confirm electron diffraction off crystal lattices. The corpuscular/wave duality is unavoidable for both light and matter; the next-generation Schrödinger and Heisenberg formalisms ([[04_quantum_mechanics_1925_1948|Chapter 4]]) inherit the problem.

### 2.6 1908, 1911 — Onnes and the cryogenic frontier

**[[onnes1908_helium_liquefaction|Helium liquefaction]].** Onnes's Leiden cryogenic lab achieves liquid helium at 4.2 K — the technical prerequisite for all later low-temperature solid-state work.

**[[onnes1911_superconductivity|Superconductivity]].** Three years later, Onnes finds that mercury's electrical resistance drops abruptly to zero below 4.19 K. *Unexplained* for 46 years until BCS (1957, [[05_QED_renormalization_solid_state_1948_1965|Chapter 5]]). The first major physical phenomenon visible only in the cryogenic regime; opens the line that runs through the Meissner effect (1933), London equations (1935), BCS, Josephson junctions (1962), and modern circuit QED.

The solid-state thread sits mostly dormant in this chapter — the experimental capability is present, the theoretical machinery (band structure, Bloch theorem, Cooper pairing) is decades away.

## 3. Proper-time commentary

### 3.1 What's directly verified

**Bohr orbit energies in the dual classical Hamiltonian.** `#verified` (cross-reference: [[FoundationsII-Classical]] Sec 2.2, which derives the critical-point claim $r_e \approx r_0$ via the dual Hamiltonian $K = H^2/(2mc^2) + mc^2/2$). For hydrogen at the Bohr orbital regime ($v/c \sim Z\alpha \approx 7 \times 10^{-3}$), the dual Hamiltonian's stationary states are numerically identical to the standard Bohr levels to within $\sim (v/c)^4 \sim 10^{-9}$ — well below 1913-era spectroscopic precision (~1 part in $10^4$).

### 3.2 What's mechanically inferred

**Sommerfeld fine structure in the dual framework.** `#inferred` from [[Dual_Relativistic_Quantum_Mechanics_I]]. The standard Sommerfeld fine-structure formula uses an SR momentum correction in a circular-orbit context. The dual framework's Eq. (I.3) dual Hamiltonian gives the same fine-structure splitting via the $b$-rescaling of the relativistic energy. **Experimental distinguishability:** dual and standard Sommerfeld predictions for hydrogen fine-structure splittings agree to $\sim 10^{-9}$ relative precision (the regime where $v/c \sim Z\alpha$ corrections matter); current hydrogen-spectroscopy precision is better, but at the precision where the standard and dual frameworks diverge, *both* are dominated by QED radiative corrections (Lamb shift, vacuum polarisation) computed equivalently in both frameworks.

**Compton scattering in the dual framework.** `#inferred` from [[Two_Mathematically_Equivalent_Versions_of_Maxwells_Equations]] + [[Dual_Relativistic_Quantum_Mechanics_I]]. For a free electron at rest, $u = 0$, $b = c$; the dual Compton formula reproduces the standard $\Delta\lambda = (h/m_e c)(1 - \cos\theta)$ exactly. For a moving target electron at $u/c \sim 10^{-2}$ (atomic-orbital velocities), the $b/c$ correction enters as a $\sim 10^{-4}$ modification of the kinematic factor, well below Compton's apparatus precision and below modern high-precision atomic-Compton-profile measurements. Will be animated in `hist_compton_null.py` — the *null* is that the dual and standard predictions coincide at any current precision for free-electron scattering.

**de Broglie matter waves in the dual framework.** `#inferred` from [[Analytic_Representation_of_The_Square-Root_Operator]]. The dual square-root operator's dispersion relation $\omega(\mathbf{k}) = c\sqrt{k^2 + (mc/\hbar)^2}$ reduces to the de Broglie form $\lambda = h/p$ at non-relativistic momenta by construction. No measurable deviation at currently achievable matter-wave interferometry (Davisson–Germer through modern atom interferometers, all at $v/c \lesssim 10^{-4}$).

### 3.3 What Gill is silent on

- **General relativity (Mercury perihelion, gravitational redshift, Schwarzschild metric).** `#gill-silent`. The dual framework has no published extension to GR. GR predictions are taken intact in [[07_PNT_GPS_SLR_QKD|Chapter 7]].
- **Helium spectrum, anomalous Zeeman effect, molecular spectra.** `#gill-silent` for old-quantum-theory specifically; these are precisely the failures that motivated the 1925 quantum mechanics revolution. The dual framework's bearing on the *new* QM is treated in [[04_quantum_mechanics_1925_1948|Chapter 4]].
- **Brownian motion, atomism, Avogadro determination.** `#gill-silent`. Classical statistical mechanics; no proper-time content.

## 4. Key derivations worth animating

| Manim scene | Status | What it shows |
|---|---|---|
| [`hist_bohr_proper_time.py`](../Animations/manim_scenes/hist_bohr_proper_time.py) | rendered | Bohr orbit derivation reframed with the dual Hamiltonian $K = H^2/(2mc^2) + mc^2/2$; the energy levels come out identical to standard Bohr at hydrogen velocities ($Z\alpha \sim 7\times 10^{-3}$); explicit numerical comparison showing agreement to ~$10^{-9}$ relative. |
| [`hist_compton_null.py`](../Animations/manim_scenes/hist_compton_null.py) | rendered | Standard Compton derivation; then the dual-framework redo with $b = \sqrt{c^2 + u^2}$; demonstration that for a target electron at rest ($u=0$) the two predictions are identical. The "null" is that the dual framework recovers standard Compton exactly for free-electron scattering. |

## 5. Primary sources cited

- [[planck1900_blackbody]] — quantum hypothesis, ad hoc.
- [[einstein1905_photoelectric]] — light quantum.
- [[einstein1905_brownian]] — atomism.
- [[einstein1905_specrel]] — special relativity.
- [[einstein1905_emc2]] — mass-energy.
- [[minkowski1908_raum_zeit]] — 4D spacetime formalism.
- [[rutherford1911_alpha_scattering]] — nuclear atom.
- [[bohr1913_atom_constitution]] — quantised atomic orbits.
- [[einstein1915_perihelion]] — Mercury perihelion confirms GR.
- [[einstein1916_grundlage]] — GR systematic exposition.
- [[sommerfeld1916_atombau]] — fine structure; first appearance of $\alpha$.
- [[millikan1916_photoelectric_verification]] — experimental confirmation of $E_k = h\nu - W$.
- [[compton1923_xray_scattering]] — photon momentum.
- [[de_broglie1924_thesis]] — matter waves.
- [[onnes1908_helium_liquefaction]] — cryogenic frontier.
- [[onnes1911_superconductivity]] — superconductivity discovered (unexplained until 1957).

## 6. Retrospective reviews drawn on

- [[whittaker1910_aether_electricity]] — continued from Ch 1–2.
- [[darrigol2000_electrodynamics_ampere_einstein]] — continued from Ch 1–2.
- [[pais1982_subtle_is_the_lord]] — Einstein biography; primary reference for §2.2 + §2.4.
- [[jammer1966_conceptual_qm]] — standard history of old quantum theory.
- [[kuhn1978_blackbody]] — revisionist reading of Planck.
- [[heilbron1969_bohr_atom]] — archival history of Bohr's 1912–13 path to the atomic model.

## 7. Cross-references

- Backward: [[02_classical_synthesis_1860_1900]].
- Forward: [[04_quantum_mechanics_1925_1948]] — matrix + wave mechanics, Dirac, antimatter.
- Verification anchors: [[FoundationsII-Classical]], [[Analytic_Representation_of_The_Square-Root_Operator]], [[Dual_Relativistic_Quantum_Mechanics_I]] (forward reference; the full dual-Dirac fine structure is in Ch 4).
- Forward (PNT): [[07_PNT_GPS_SLR_QKD]] uses §2.4's Mercury-perihelion derivation as its pedagogical bridge to GPS.
- Findings: none in this era; all three flagged findings ([[FINDINGS_for_author_review]]) sit in later papers.
