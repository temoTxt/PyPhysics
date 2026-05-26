# Candidate 2 — overnight state log

## Iteration 0 — 2026-05-26 — initialized

- Branch `65-theory-candidate-2-variational-determination-of-r_e-via-renormalised-dual-dirac-equation` checked out.
- `.dev/research/brief.md` written.
- `.dev/research/loop_prompt.md` written.
- No derivation work yet.
- **Next:** read source-of-record §1 (`Dual_Relativistic_Quantum_Mechanics_I.md` §II and §III), record the eigenvalue-problem structure with $r_e$ as a free parameter in this STATE log.
- **Outcome-matrix:** not yet determinable.
- **Status:** READY.

---

## Iteration 1 — 2026-05-26T17:30Z — eigenvalue-problem structure recorded

**Advanced:** Read `Roadmapping/Equation_Verification/Dual_Relativistic_Quantum_Mechanics_I.md` §II (Eqs. II.1–II.3) and §III (Eqs. III.1–III.23). Recording here the eigenvalue-problem structure with $r_e$ as a free parameter; this is the algebraic substrate the variational determination must operate on.

### Substrate equations (from DRQM I §III)

The dual Dirac eigenvalue equation reads (Eq. III.4):
$$E\,\Psi \;=\; \left\{\,\frac{\boldsymbol\pi^{2}}{2m} + \beta V_{0} + mc^{2} - \frac{e\hbar\boldsymbol\Sigma\!\cdot\!\mathbf{B}}{2mc} + \frac{V_{0}\boldsymbol\alpha\!\cdot\!\boldsymbol\pi}{mc} - \frac{i\hbar\boldsymbol\alpha\!\cdot\!\nabla V_{0}}{2mc} + \frac{V_{0}^{2}}{2mc^{2}}\,\right\}\Psi.$$
With Coulomb $V_{0} = -e^{2}/r$ and the substitution $\psi_{2} = c(\lambda - V_{0} + mc^{2})^{-1}(\boldsymbol\sigma\!\cdot\!\boldsymbol\pi)\psi_{1}$ (Eq. III.2), the 2-component reduction yields a single equation for $\psi_{1}$ in which $(\lambda - V_{0} + mc^{2})^{-1}$ is the *only* place the cutoff parameter $r_e$ can enter.

### Where $r_e$ enters

**Eq. (III.7) cutoff approximation.** Using $\lambda - mc^{2} \ll mc^{2}$ (binding energy $\sim 13$ eV vs $mc^{2} \sim 5\times 10^{5}$ eV; relative error $\sim 10^{-5}$):
$$\lambda - V_{0} + mc^{2} \;\approx\; 2mc^{2} + \frac{e^{2}}{r} \;=\; 2mc^{2}\!\left(1 + \frac{r_{0}}{2r}\right),\qquad r_{0} \equiv \frac{e^{2}}{mc^{2}}.$$
This identifies $r_{0}$ as the natural radial scale (classical electron radius in Gaussian units, modulo a factor of $4\pi\varepsilon_0$). The denominator $(1 + r_{0}/(2r))$ is what carries the $r$-dependence through the rest of the derivation.

**Eqs. (III.18)–(III.20).** After spherical-coordinate expansion of $-i\hbar\boldsymbol\alpha\!\cdot\!\nabla V_{0}$ and the $V_{0}\boldsymbol\alpha\!\cdot\!\boldsymbol\pi/(mc)$ chain-rule term, the three new contributions all carry an explicit $1/(2r + r_{0})$ structure (the denominator from III.7 after pulling out $2mc^{2}$). The radial cutoff $r_{e}$ enters when one *evaluates these operators at* $r = r_{e}$ rather than integrating across all $r$, i.e. interprets $r_{e}$ as the hard lower cutoff on the radial integration domain.

**Eq. (III.22) — the g-factor formula.** Collecting the spin–field term from (III.18) at $r = r_{e}$:
$$g_{r}(r_{e}) \;=\; 2\!\left[1 - \frac{4r_{0}}{2r_{e} + r_{0}}\right].$$
This is the *one* numerical constraint already in the paper that determines $r_{e}/r_{0}$ from observable data. Inverting against the measured $g_{e} = -2.00231930436256$ gives $r_{e}/r_{0} = 0.4994205099128318$ (triangulated value from PR #62, confirmed independently by joint fit across six observables to $\sigma_r = 2.5\times 10^{-13}$).

### Variational-route framing

The candidate-2 question is whether $r_{e}/r_{0}$ can be fixed *without* reference to $g_{e}$, by demanding the dual-Dirac equation's eigenvalue itself reproduce $m_{e}c^{2}$ (plus framework-internal binding contributions) at the cutoff. Concretely, the variational functional is
$$E[r_{e}] \;\equiv\; \frac{\langle\Psi_{r_{e}} | K_{D} | \Psi_{r_{e}}\rangle}{\langle\Psi_{r_{e}} | \Psi_{r_{e}}\rangle}$$
where $\Psi_{r_{e}}$ is the ground-state spinor evaluated with the radial integration domain $[r_{e}, \infty)$ (or equivalently with $1/r \to r/(r^{2} + r_{e}^{2})$ or similar regularisation; the choice of regulator is itself a closure question to be flagged). The renormalisation condition is
$$E(r_{e}) \;=\; m_{e}c^{2} + \Delta E_{\rm bind}^{\rm framework} + \Delta E_{\rm SE}^{\rm framework},$$
where the framework-internal binding and self-energy contributions $\Delta E_{\rm bind}^{\rm framework}$, $\Delta E_{\rm SE}^{\rm framework}$ are *not yet specified* in the paper's published apparatus. **This is the load-bearing gap.**

### Closure-condition inventory (preview, to be expanded next iteration)

Candidate closure conditions, to be classified `framework-internal` vs `ad-hoc`:

1. **Stationarity of $E[r_{e}]$ w.r.t. $r_{e}$** — `framework-internal` if and only if the renormalisation prescription identifies the cutoff as a variational parameter (analogous to mass renormalisation $\partial E/\partial m_{\rm bare} = 0$ in QED); `ad-hoc` if just imposed.
2. **Gauge invariance at the cutoff** — current conservation $\partial_{\mu} J^{\mu}|_{r=r_{e}} = 0$ at the radial boundary. Framework-internally motivated *if* DRQM I requires the radial boundary to preserve $U(1)$. **Status: not yet checked in paper text.**
3. **Magnetic-moment relation closure** — demanding the spin-$\boldsymbol\sigma\!\cdot\!\mathbf{B}$ coefficient at $r=r_{e}$ match an independent framework-internal prediction (e.g. from proper-time dynamics rather than from $g_{e}$ data). `framework-internal` if such an independent prediction exists.
4. **Normalisation closure** — $\int_{r_{e}}^{\infty} |\Psi|^{2}\, d^{3}r = 1$ as a non-trivial condition fixing $r_{e}$. *Likely ad-hoc* — normalisation alone has no preferred $r_{e}$.

**Next:** Enumerate candidate closure conditions in full and start the framework-internal-vs-ad-hoc classification, beginning with #2 (gauge invariance at the cutoff). Read DRQM I §II.A and the companion verification doc for the proper-time Maxwell equations to see whether the framework specifies current conservation at a radial boundary.

**Outcome-matrix:** not yet determinable (need at least one framework-internal closure condition before any branch is reachable).

**Status:** READY for next iteration.

---

## Iteration 2 — 2026-05-26T17:55Z — closure-condition inventory + framework-internal/ad-hoc classification

**Advanced:** Reviewed companion verification docs for any framework-specified boundary/gauge structure at the cutoff. Key finding: **FoundationsII-Classical Sec 2.2** (verified in `Roadmapping/Equation_Verification/FoundationsII-Classical.md` Eq. 2.11) establishes that **$r = r_0$ is a critical point of the dual Hamiltonian force** $\mathbf{F}_K = -\nabla V(1 + V/(mc^2))$, *within the framework's H-as-part-of-mass form* $H = \sqrt{c^2\mathbf{p}^2 + (mc^2 + V)^2}$. This is the only framework-internal preferred radial scale beyond $r_0$ itself that I have located. The Maxwell paper §I.D 4-current $J^\mu_{\rm Gill} = (b\rho, \mathbf{J})$ (Maxwell paper Eq. 12, with $\mathbf{J} = (c/b)\rho\mathbf{u}$) implies a *bulk* continuity equation in proper-time, but I found no published radial-boundary continuity condition.

### Closure-condition inventory (full enumeration)

| # | Condition | Framework-internal evidence | Classification | Notes |
|---|---|---|---|---|
| 1 | **Variational stationarity** $\partial E[r_e]/\partial r_e = 0$ | None — DRQM I does *not* state $r_e$ is determined by stationarity of any functional. The Dresden-renormalisation-analog remark on Maxwell paper line 223 ([Two_Mathematically_Equivalent_Versions_of_Maxwells_Equations.md:696](../Roadmapping/Equation_Verification/Two_Mathematically_Equivalent_Versions_of_Maxwells_Equations.md)) hints at mass-renormalisation but does not pin a variational principle on $r_e$. | **ad-hoc** (unless author confirms) | Tractable as scaffolding for the Mathematica route; cannot terminally fix an outcome-A/B result without author endorsement. |
| 2 | **Radial-boundary current conservation** $\mathbf{J}\!\cdot\!\hat{\mathbf{r}}\,\big|_{r=r_e} = 0$ | Bulk continuity is guaranteed by the dual-current structure (Maxwell paper Eqs. 12–15), but no boundary form is published. | **ad-hoc** pending author input | Would force a node in $\Psi$ at $r=r_e$; numerically over-constrains. |
| 3 | **$g$-factor closure** $g_r(r_e) = g_e^{\rm exp}$ | This is what the published paper does — uses experimental input. | **framework-external** | Already executed in PR #62 (triangulated $r_e/r_0 = 0.4994205099128317$). Not eligible as a *first-principles* derivation. |
| 4 | **Critical-point locking** $r_e = r_0/2$ | $r_0$ is the critical point of $\mathbf{F}_K$ (FoundationsII-Classical Eq. 2.11). A "midpoint" cutoff $r_e = r_0/2$ would lock into the dual-Hamiltonian's structural geometry. | **framework-internal but wrong-precision** | Direct plug-in: $g_r(r_0/2) = 2[1 - 4r_0/(r_0+r_0)] = -2$. Misses the Schwinger $\alpha/\pi$ correction. So this condition yields the Dirac-tree-level $g=-2$ exactly; it cannot reproduce $0.4994$. |
| 5 | **Normalisation closure** $\int_{r_e}^{\infty} |\Psi|^2 d^3 r = 1$ | None — trivially achievable for any $r_e$ via rescaling. | **ad-hoc / no info** | Discard. |
| 6 | **Schwinger one-loop closure** $r_e/r_0 = (2 - \alpha/(2\pi))/(4 + \alpha/\pi)$ | Inverts $g_r$ against the Schwinger $g_e^{(1\text{-loop})} = -2 - \alpha/\pi$ analytically; $\sim 10^{-6}$ from triangulated. | **framework-external** | This is Candidate 3's route. Requires the framework's renormalisation prescription to produce an $\alpha$-dependent ratio internally — not established. |
| 7 | **Energy-eigenvalue mass-renormalisation** $\langle K_D \rangle_{r_e} = m_e c^2 + \Delta E_{\rm bind}^{\rm framework} + \Delta E_{\rm SE}^{\rm framework}$ | The renormalisation-condition language in the candidate-2 brief; analogous to QED's $\alpha(\mu_0) = \alpha_{\rm phys}$ at a chosen reference scale. | **framework-internal IFF** the framework specifies $\Delta E_{\rm bind}^{\rm framework}$ and $\Delta E_{\rm SE}^{\rm framework}$ — currently *unspecified* in DRQM I as published. | **Load-bearing.** This is the only framework-internal candidate with non-trivial dynamics. Requires author input or independent identification of the framework's binding/SE prescription. |

### Classification summary

- **Framework-internal & sufficient-precision:** *none yet* — #4 is internal but yields only $g=-2$; #7 is internal but contains an unspecified term.
- **Framework-internal & insufficient-precision:** #4 (critical-point locking → $g=-2$ exactly, off by $\sim 10^{-3}$).
- **Framework-external:** #3, #6 (use measured $g_e$ or Schwinger one-loop respectively).
- **Ad-hoc:** #1, #2, #5.

### Forward path

Two routes to break the impasse without author input:

**Route X — Push #7 with a *minimal* framework-internal binding/SE specification.** Try $\Delta E_{\rm bind} = $ Coulomb expectation value $\langle V_0 \rangle = -e^2 \langle 1/r \rangle$ evaluated on the cutoff-truncated ground state; $\Delta E_{\rm SE} = 0$ (no separate SE at this stage). If $\langle K_D \rangle_{r_e} = m_e c^2 + \langle V_0\rangle_{r_e}$ produces a tractable equation, solve symbolically and compare against $0.4994205099$. The classification of this result depends on whether the framework endorses "$\Delta E_{\rm bind} = \langle V_0\rangle$" as the renormalisation prescription (probably yes — it's the standard textbook reading, and DRQM I's perturbation structure does not contradict it).

**Route Y — Sharpen the #4 critical-point analysis.** $r_0$ is a critical point of $\mathbf{F}_K$. Where is the critical point of *the dual-Dirac eigenvalue* $\lambda$ as a function of $r$? If $\partial\lambda/\partial r|_{r=r_e} = 0$ has a definite solution, that's a fully framework-internal closure condition with no external input. The eigenvalue $\lambda$ comes from solving Eq. (III.4)–(III.8) under the radial-cutoff convention; this requires actually setting up the radial Dirac equation with the cutoff regulator. Tractable but lengthy.

**Next:** Begin Route X. Create `Roadmapping/Mathematica_Notebooks/Quantum_Mechanics/r_e_derivation_variational.wl` (`.wl` style template per CLAUDE.md, single-line Wolfram cells, avoid `V`/`e`/`Dot` pitfalls). First cell: set up the radial Dirac equation with the cutoff regulator (radial domain $[r_e, \infty)$, hydrogenic trial $\psi_1 = N e^{-r/a}/r$ with $a$ to be optimised). Second cell: compute $\langle K_D \rangle_{r_e}$ symbolically for this trial. **Goal of next iteration is the file scaffolding + first cell only** — full eigenvalue evaluation is a 2-3 iteration arc.

**Outcome-matrix:** still not yet determinable. Framework-internal closure exists in principle (#7 with minimal binding spec, or #4 sharpened to dual-Dirac eigenvalue criticality), but no definite $r_e/r_0$ has been computed. **No BLOCKED state yet** — Routes X and Y are both pursuable without author input; if Route X produces an $r_e/r_0$ that matches triangulated, that's outcome A; if it matches Schwinger, outcome B; if it produces a different value, outcome C; if intractable, then BLOCKED on author input for the binding/SE prescription.

**Status:** READY for next iteration (Route X start).

---

## Iteration 3 — 2026-05-26T18:18Z — Route X start: `r_e_derivation_variational.wl` scaffold + Section 1 (norm)

**Advanced:** Created `Roadmapping/Mathematica_Notebooks/Quantum_Mechanics/r_e_derivation_variational.wl` with:

1. **Header docstring** documenting source-of-record (DRQM I §III, FoundationsII-Classical §2.2, STATE.md iter 1-2), Route X framing, honest scope (substantive AI), trial-wavefunction choice (exp(-r/aa) soft cutoff), and Wolfram-MCP gotchas per CLAUDE.md.
2. **Section 1: Symbol setup + trial wavefunction.** Defines $r_0 = e^2/(mc^2)$ (classical electron radius, Gaussian units), trial $\psi_1(r; aa) = \exp(-r/aa)$ on radial domain $[r_e, \infty)$ (soft cutoff — amplitude truncated below $r_e$, not forced to zero).
3. **Section 2: Cutoff-restricted normalisation.** Computed the radial integral $\int_{r_e}^\infty r^2 e^{-2r/aa}\, dr$ via Wolfram MCP. **Closed form:**
$$\int_{r_e}^{\infty} r^2 e^{-2r/aa}\, dr \;=\; \frac{aa\,(aa^2 + 2\,aa\,r_e + 2\,r_e^2)}{4\,e^{2r_e/aa}}.$$
**Sanity check $r_e \to 0$:** $\to aa^3/4$, matching the standard 1s hydrogen norm. ✅ (Wolfram MCP 2026-05-26).
4. **Section 3: PLACEHOLDER for $\langle K_D\rangle$.** Documented the plan: for the field-free s-state, only scalar terms in (III.4) survive (no B, no spin–orbit), so
$$\langle K_D\rangle_{r_e, aa} = \langle \pi^2/(2m)\rangle + mc^2 + \langle V_0\rangle + \langle V_0^2/(2mc^2)\rangle.$$
Closure #7 collapsed to its scalar reading: with $\Delta E_{\rm bind} = \langle V_0\rangle$, the mass-renormalisation condition becomes
$$\langle \pi^2/(2m)\rangle_{r_e, aa} + \langle V_0^2/(2mc^2)\rangle_{r_e, aa} = 0.$$
This is one equation for two unknowns $(aa, r_e)$; the second condition is variational stationarity $\partial\langle K_D\rangle/\partial aa = 0$ on the trial.
5. **Human-acceptance stub** (Crocco): three substantive-AI choices flagged — trial-wavefunction form, the reading of $\Delta E_{\rm bind}$ as $\langle V_0\rangle$, and the choice of $aa$-stationarity as the companion condition.

**Closure-condition classification update:** Condition #7 collapsed to the form $\langle\pi^2/(2m)\rangle + \langle V_0^2/(2mc^2)\rangle = 0$. Whether this is solvable for finite $(aa, r_e)$ is the question for next iteration. Currently *tentatively framework-internal* — the reading $\Delta E_{\rm bind} = \langle V_0\rangle$ is the textbook default for non-relativistic mass renormalisation and consistent with the framework's "$V$ as part of the mass" form (FoundationsII-Classical §2.2), but the dual framework's published apparatus does not pin this prescription explicitly. **Will become BLOCKED on author input** if the calculation yields an $r_e/r_0$ that disagrees with both the triangulated and Schwinger values — then we'll need to determine whether $\Delta E_{\rm bind} \ne \langle V_0\rangle$ in this framework.

**Next:** Compute the three remaining expectation values symbolically in Wolfram MCP:
- $\langle V_0\rangle_{r_e,aa} = -e^2 \langle 1/r\rangle_{r_e,aa}$ (one radial integral with $1/r$)
- $\langle V_0^2/(2mc^2)\rangle_{r_e,aa} = e^4/(2mc^2) \langle 1/r^2\rangle_{r_e,aa}$ (one radial integral with $1/r^2$)
- $\langle \pi^2/(2m)\rangle_{r_e,aa}$ — kinetic energy expectation: the trial $\psi_1 = N e^{-r/aa}$ has $\nabla^2 \psi_1 = (1/aa^2 - 2/(aa\, r))\psi_1$, so $\langle -\hbar^2 \nabla^2/(2m)\rangle = (\hbar^2/(2m))[1/aa^2 - (2/aa)\langle 1/r\rangle]$. Two integrals reused from above.

Then assemble closure equation $\langle\pi^2/(2m)\rangle + \langle V_0^2/(2mc^2)\rangle = 0$ in symbolic form. **Goal of next iteration: write Section 3 of the notebook with these three integrals + the assembled closure equation; no numerical solve yet.**

**Outcome-matrix:** still not yet determinable — full closure equation not yet assembled. **No BLOCKED state.**

**Status:** READY for next iteration (Section 3 expectation-value integrals).
