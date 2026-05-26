# Candidate 1 — overnight state log

## Iteration 0 — 2026-05-26 — initialized

- Branch `64-theory-candidate-1-proper-time-self-energy-integral-derivation-of-r_e` checked out.
- `.dev/research/brief.md` written.
- `.dev/research/loop_prompt.md` written.
- No derivation work yet.
- **Next:** read source-of-record §1 (`Dual_Relativistic_Quantum_Mechanics_I.md`), record key equations from §II and §III.D in this STATE log.
- **Outcome-matrix:** not yet determinable.
- **Status:** READY.

## Iteration 1 — 2026-05-26 — DRQM I §II + §III.D ingested

**Advanced.** Read `Roadmapping/Equation_Verification/Dual_Relativistic_Quantum_Mechanics_I.md` (591 lines) end-to-end. Recorded below the structural identities needed to set up a one-loop self-energy diagram in the dual framework.

### Conventions confirmed (Gaussian, c.g.s.)

- $\mathbf{u} = d\mathbf{x}/d\tau = \gamma\mathbf{w}$; $b = \sqrt{c^2+\mathbf{u}^2}$; $cdt = b\,d\tau$; $(1/c)\partial_t = (1/b)\partial_\tau$.
- $\boldsymbol\pi = \mathbf{p} - e\mathbf{A}/c$.
- $r_0 := e^2/(mc^2)$ (the classical electron-radius scale used as Coulomb-cutoff yardstick).

### §II — three dual relativistic wave equations, common form $i\hbar\,\partial_\tau\Psi = K\Psi$ with $K = H^2/(2mc^2) + mc^2/2$.

- **(II.1) Dual Dirac**, with $H_{\rm D} = c\boldsymbol\alpha\!\cdot\!\boldsymbol\pi + \beta mc^2 + V$:
  $$K_{\rm D} = \frac{\boldsymbol\pi^2}{2m} + \beta V + mc^2 - \frac{e\hbar\boldsymbol\Sigma\!\cdot\!\mathbf{B}}{2mc} + \frac{V\boldsymbol\alpha\!\cdot\!\boldsymbol\pi}{mc} - \frac{i\hbar\boldsymbol\alpha\!\cdot\!\nabla V}{2mc} + \frac{V^2}{2mc^2}.$$
- **(II.2) Sqrt(1)**, $H_{\rm s1} = \beta S + V$ with $S = \sqrt{c^2\boldsymbol\pi^2 - ec\hbar\boldsymbol\Sigma\!\cdot\!\mathbf{B} + m^2c^4}$ — operator-valued $\{V,S\}$ structure.
- **(II.3) Sqrt(2) "potential-in-the-mass"** — *cleanest*, $H_{\rm s2} = \beta\sqrt{c^2\boldsymbol\pi^2 - ec\hbar\boldsymbol\Sigma\!\cdot\!\mathbf{B} + (mc^2+\beta V)^2}$, gives
  $$K_{\rm s2} = \frac{\boldsymbol\pi^2}{2m} + \beta V + mc^2 - \frac{e\hbar\boldsymbol\Sigma\!\cdot\!\mathbf{B}}{2mc} + \frac{V^2}{2mc^2}.$$

  *Implication for self-energy.* (II.3) has *no* $V\,\boldsymbol\alpha\!\cdot\!\boldsymbol\pi$ or $\boldsymbol\alpha\!\cdot\!\nabla V$ pieces — the dual proper-time kinetic kernel reduces to a Pauli kinetic term + scalar $V^2/(2mc^2)$ relativistic correction. This is the kernel the proper-time integral should propagate (one-loop self-energy on the (II.3) Hamiltonian, since the (II.1) form double-counts the spin–orbit/Darwin terms that arise from operator non-commutativity).

### §III.D — where $r_e$ enters

- **(III.7) cutoff approximation.** $\lambda - V_0 + mc^2 \approx 2mc^2(1 + r_0/(2r))$ with $V_0 = -e^2/r$. The denominator therefore becomes $2mc^2(1+r_0/(2r))$ everywhere $\psi_2$ appears in the small-component elimination.
- **(III.18) "new" $g$-factor contribution.** After collecting all expanded terms:
  $$\Delta_{g}\,H = -\Bigl[1 - \frac{4r_0}{2r+r_0}\Bigr]\frac{e\hbar\,\boldsymbol\sigma\!\cdot\!\mathbf{B}}{2mc}.$$
- **(III.21–22) $g$-factor.** $g_r(r) = 2[1 - 4r_0/(2r+r_0)]$. Limits: $g_r(r_0/2) = -2$; $g_r(0) = -6$. Cutoff distance $r_e$ is the spatial scale where the $\psi_2$-substitution stops being trusted — physically a radiative-correction cutoff.
- **(III.22) target.** Triangulated empirical value: $r_e/r_0 = 0.499\,420\,509\,912\,831\,7$ (PR #62, joint fit across 6 observables; $\sigma_r = 2.5\times10^{-13}$).
- **(III.23) muon and proton analogues** use $r_0^{\mu,p} = e^2/(m_{\mu,p}c^2)$ — same functional form, different mass scale.

### Schwinger closed-form near-match

$$g_S(\alpha) := \frac{2 - \alpha/(2\pi)}{4 + \alpha/\pi}\Big|_{r=r_e\,{\rm root}}\!\!\!\Rightarrow r_e/r_0 = 0.499\,419\,632\,156\ldots,$$

i.e. agrees with the triangulated value at the $\sim 10^{-6}$ level — the residual ($\sim 8.8\times 10^{-7}$) is the gap a higher-loop / framework-modified self-energy must explain.

### Self-energy strategy (working hypothesis for next iterations)

A one-loop self-energy in the (II.3) kernel, written as a Schwinger proper-time integral $\Sigma(p) = \int_0^\infty ds\,e^{-s\,\mathcal{O}(p,m)}\,\ldots$, with the dual proper-time relation $cdt = b\,d\tau$ providing a *natural* lower-cutoff at $s_{\min} \sim r_e/c$ (in proper-time units). Mass-renormalisation condition: $m_{\rm phys} = m + \Sigma(p)\!\restriction_{p^2 = m^2c^2}$, fixing the cutoff $s_{\min}$ (equivalently $r_e$) by matching to the framework's tree-level mass parameter. This is the open thread to develop iteratively.

**Next.** Read source-of-record §2 — `Roadmapping/Equation_Verification/The_Classical_Electron_Problem.md` — and extract the proper-time radiation-reaction structure (per #55, a candidate ingredient at radiative-correction order). Record the proper-time photon-propagator form and any framework-internal mass-renormalisation conditions Tepper has previously stated.

**Outcome-matrix:** not yet determinable (still scoping the kernel; no numerical $r_e/r_0$ yet).

**Status:** READY.
