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
