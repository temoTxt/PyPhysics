# Dual Relativistic Quantum Mechanics I

Tepper L. Gill[∗](#page-0-0) *EECS, Mathematics and Computational Physics Laboratory, Howard University Washington DC 20059 USA*

Gonzalo Ares de Parga[†](#page-0-1)

*Departmento de F´ısica, Escuela Superior de F´ısica y Matem´aticas, Instituto Polit´ecnico Nacional COFAA; Edif 9, U. P. Adolfo Løpez Mateos, ´ Zacatenco, Lindavista, 07738, M´exico D.F., M´exico*

Trey Morris[‡](#page-0-2) and Mamadou Wade[§](#page-0-3) *Department of EECS, Howard University Washington DC 20059 USA* (Dated: August 21, 2021)

It was shown in [\[1\]](#page-6-0) that the ultra-violet divergence in quantum electrodynamics (QED) is caused by a violation of the time-energy uncertainly relationship, due to the implicit assumption of infinitesimal time information (Dyson's conjecture). In [\[2\]](#page-7-0) it was shown that Einstein's special theory of relativity and Maxwell's field theory have mathematically equivalent dual versions. The dual versions arise from an identity relating observer time to proper time as a contact transformation on configuration space, which leaves phase space invariant. The special theory has a dual version in the sense that, for any set of n particles, every observer has two unique sets of global variables (X, t) and (X, τ ) to describe the dynamics, where X is the (unique) canonical center of mass. In the (X, t) variables, time is relative and the speed of light is unique, while in the (X, τ ) variables, time is unique and the speed of light is relative with no upper bound. In the Maxwell case, the two sets of particle wave equations are not equivalent. The dual version contains an additional longitudinal (dissipative) radiation term that appears instantaneously with acceleration, leading to the prediction that radiation from a betatron (of any frequency) will not produce photoelectrons. A major outcome is the dual unification of Newtonian mechanics and classical electrodynamics with Einstein's special theory of relativity, without a self-energy divergency, or need of the problematic Lorentz-Dirac equation or any assumptions about the size or structure of a particle. The purpose of this paper is to introduce and develop the dual theory of relativistic quantum mechanics. We obtain three distinct dual relativistic wave equations that reduce to the Schr¨odinger equation when minimal coupling is turned off. We show that the dual Dirac equation provides a new formula for the anomalous magnetic moment of a charged particle. We can obtain the exact value for the electron g-factor and phenomenological values for the muon and proton g-factors.

## INTRODUCTION

In classical electrodynamics, Dirac partially by-passed many of the problems left over from the nineteenth century by replacing particles by fields (see [\[3\]](#page-7-1)). This approach led to the first example of a divergent theory (infinite self-energy). Dirac showed that, by using both advanced and retarded fields and a limiting procedure, one obtained a dissipative term, which accounted for the radiation reaction problem as an addition to the Lorentz equation (Lorentz-Dirac equation). This self-energy divergency was the main motivation for the Wheeler-Feynman approach to classical electrodynamics (see [\[4\]](#page-7-2)). Their theory gave the same dissipative term while avoiding the self-energy divergency. This approach could not be quantized, but provided insight for Feynman's approach to QED.

The failure to directly solve the classical problems forced researchers to use the Dirac theory as the basis for relativistic quantum mechanics and QED. This program maintained the self-energy divergence and introduced a few others. These divergences were later by-passed by Feynman, Schwinger and Tomonaga in the late 1940's leading to the great successes of that era. Neither Feynman, Schwinger or Tomonaga considered their work a complete theory or final solution. Their methods did not account for the full spectra Hydrogen, but required the solution of the eigenvalue problem from the Dirac equation as initial input.

The predominant belief at the time was that they were on the right track, and the remaining difficulties would eventually be cleaned up by the mathematicians. However, the major mathematical investigations were restricted to the limited task of providing justification for the subtraction methods used to handle the divergencies. This justification never came and by the early 1980's, it became clear that it never would. The development of the electro-weak theory and the standard model each added additional problems, caused by extensions of QED methods to higher energy scales.

<span id="page-0-0"></span><sup>∗</sup> [tgill@howard.edu](mailto:tgill@howard.edu)

<span id="page-0-1"></span><sup>†</sup> [gadpau@hotmail.com](mailto:gadpau@hotmail.com)

<span id="page-0-2"></span><sup>‡</sup> [Morris.Trey.J@gmail.com](mailto:Morris.Trey.J@gmail.com)

<span id="page-0-3"></span><sup>§</sup> [mamadou.wade@howard.edu](mailto:mamadou.wade@howard.edu)

Another (less known) line of investigation sought to directly deal with the physical cause for these problems based on a number of suggestions from Dirac, Dyson and Feynman (see [1]). A major outcome was that the ultra-violet divergency came from a violation of the time-energy uncertainty relationship (as suggested by Dyson) and was not a hint of some (unknown) deeper problems as many believed.

The success of this and other efforts suggested that an investigation into the physical justification for time as a forth coordinate was in order. The lack of justification resulted in the discovery of the dual theory of special relativity [2] and the dual Maxwell theory [6, 7]. The dual Maxwell theory identifies radiation reaction as a dissipative term in the E field equation. This led to the elimination of three major problems with the Dirac version of classical electrodynamics: the self-energy divergence disappeared, the need for point particles disappeared and, the need for the (problematic) Lorentz-Dirac equation disappeared. In addition, it was shown that the dissipative term is equivalent to a small (dynamical) mass for the photon. This latter property implies that a quantum field based on the dual Maxwell theory will not lead to an infrared-divergence (see [6]).

The purpose of this paper is to introduce the dual relativistic quantum theory. After a brief review of the dual single particle and Maxwell theory in the second section, we introduce the dual relativistic quantum theory in the third section. The fourth section is devoted to the dual Dirac equation, which provides a new formula for the g-factor of the electron and can also be used to obtain exact (phenomenological) values for the muon and proton g-factors.

#### I. DUAL CLASSICAL THEORY

### A. Particle Clock

To develop the dual classical theory, we assume a classical interacting particle defined on phase space with variables  $(t, \mathbf{x}, \mathbf{p})$  and Hamiltonian H as seen by an observer O in an inertial frame (in the standard setup). If  $\mathbf{w}$  is the particle velocity, let  $\gamma^{-1}(\mathbf{w}) = \sqrt{1 - \mathbf{w}^2/c^2}$ . The classical proper time is defined by  $d\tau = \sqrt{1 - \mathbf{w}^2/c^2}dt$ ,

$$\mathbf{w} = \frac{d\mathbf{x}}{dt}, \quad d\tau^2 = dt^2 - \frac{1}{c^2}d\mathbf{x}^2 \tag{I.1}$$

Rearranging the last term, we get  $dt^2 = d\tau^2 + \frac{1}{c^2}d\mathbf{x}^2$ , so

<span id="page-1-0"></span>
$$cdt = \left(\sqrt{\mathbf{u}^2 + c^2}\right) d\tau, \quad \mathbf{u} = \frac{d\mathbf{x}}{d\tau} = \gamma(\mathbf{w})\mathbf{w}.$$
 (I.2)

If we let  $b = \sqrt{\mathbf{u}^2 + c^2}$ , the first term in equation (I.2) becomes  $cdt = bd\tau$ . This leads to our first identity:

<span id="page-1-1"></span>
$$\frac{1}{c}\frac{d}{dt} \equiv \frac{1}{b}\frac{d}{d\tau}.$$
 (I.3)

This identity provides the correct way to define the relationship between the proper time and observer time for the particle. If we apply the identity to  $\mathbf{x}$ , we obtain our second new identity, showing that the transformation leaves the configuration (or tangent) space invariant:

<span id="page-1-2"></span>
$$\frac{\mathbf{w}}{c} = \frac{1}{c} \frac{d\mathbf{x}}{dt} \equiv \frac{1}{b} \frac{d\mathbf{x}}{d\tau} = \frac{\mathbf{u}}{b}.$$
 (I.4)

The new particle coordinates are  $(\mathbf{x}, \tau)$ . In this representation, the position  $\mathbf{x}$  is uniquely defined relative to O, while  $\tau$  is uniquely defined by the particle. The particle momentum can be represented as  $\mathbf{p} = m\gamma(\mathbf{w})\mathbf{w} = m\mathbf{u}$ , where m is the particle rest mass. Thus, the phase space variables  $(\mathbf{x}, \mathbf{p})$ , are left invariant. For later use, we also have  $\gamma(\mathbf{w}) = H/mc^2$ . This allows us to also write  $d\tau = (mc^2/H) dt$ .

#### B. Dual Particle Theory

For compatibility with quantum theory, we require that any change of clocks be canonical. The key concept to our approach is seen by examining the time evolution of a dynamical parameter  $W(\mathbf{x}, \mathbf{p})$ , via the standard formulation in terms of the Poisson brackets:

$$\frac{dW}{dt} = \{H, W\}. \tag{I.5}$$

To represent the dynamics via the proper time of the particle, we use the representation  $d\tau = (mc^2/H)dt$ , so that:

$$\frac{dW}{d\tau} = \frac{dt}{d\tau} \frac{dW}{dt} = \frac{H}{mc^2} \left\{ H, W \right\}.$$

Using the invariant rest energy  $mc^2$ , we determine the canonical proper-time Hamiltonian K such that:

$$\{K, W\} = \frac{H}{mc^2} \{H, W\}, \quad K|_{\mathbf{p}=0} = H|_{\mathbf{p}=0} = mc^2.$$

From

$$\begin{split} \{K,W\} &= \left[\frac{H}{mc^2}\frac{\partial H}{\partial \mathbf{p}}\right]\frac{\partial W}{\partial \mathbf{x}} - \left[\frac{H}{mc^2}\frac{\partial H}{\partial \mathbf{x}}\right]\frac{\partial W}{\partial \mathbf{p}} \\ &= \frac{\partial}{\partial \mathbf{p}}\left[\frac{H^2}{2mc^2} + a\right]\frac{\partial W}{\partial \mathbf{x}} - \frac{\partial}{\partial \mathbf{x}}\left[\frac{H^2}{2mc^2} + a'\right]\frac{\partial W}{\partial \mathbf{p}}, \end{split}$$

we see that  $a = a' = \frac{1}{2}mc^2$ . Thus, assuming no explicit time dependence, we have:

<span id="page-1-3"></span>
$$K = \frac{H^2}{2mc^2} + \frac{mc^2}{2}$$
, and  $\frac{dW}{d\tau} = \{K, W\}$ . (I.6)

Since  $\tau$  remains invariant during interaction (minimal coupling), we assume K also remains invariant. Thus, if  $\sqrt{c^2\mathbf{p}^2 + m^2c^4} \to \sqrt{c^2\boldsymbol{\pi}^2 + m^2c^4} + V$ , where  $\boldsymbol{\pi} = \mathbf{p} - \frac{e}{c}\mathbf{A}$ , with  $\mathbf{A}$  the vector potential and V the potential energy. In this case, K becomes:

$$K = \frac{\pi^2}{2m} + mc^2 + \frac{V^2}{2mc^2} + \frac{V\sqrt{c^2\pi^2 + m^2c^4}}{mc^2}$$

If we set  $H_0 = \sqrt{c^2 \pi^2 + m^2 c^4}$ , use standard vector identities with  $\nabla \times \pi = -\frac{e}{c} \mathbf{B}$ , and compute Hamilton's equations, we get:

$$\frac{d\mathbf{x}}{d\tau} = \frac{\partial K}{\partial \mathbf{p}} = \frac{H}{mc^2} \left( \frac{c^2 \boldsymbol{\pi}}{H_0} \right) = \frac{b}{c} \left( \frac{c^2 \boldsymbol{\pi}}{H_0} \right) \Rightarrow \frac{d\mathbf{x}}{d\tau} = \frac{b}{c} \frac{d\mathbf{x}}{dt}$$

and

$$\frac{d\mathbf{p}}{d\tau} = \frac{b}{c} \frac{\left[ \left( c^{2}\boldsymbol{\pi} \cdot \boldsymbol{\nabla} \right) \mathbf{A} + \frac{e}{b} \left( c^{2}\boldsymbol{\pi} \times \mathbf{B} \right) \right]}{H_{0}} - \frac{b}{c} \boldsymbol{\nabla} V$$

$$= \frac{b}{c} \left[ \left( \mathbf{u} \cdot \boldsymbol{\nabla} \right) \mathbf{A} + \frac{e}{b} \left( \mathbf{u} \times \mathbf{B} \right) \right] - \frac{b}{c} \boldsymbol{\nabla} V$$

$$= \frac{b}{c} \left[ e\mathbf{E} + \frac{e}{b} \left( \mathbf{u} \times \mathbf{B} \right) + \frac{e}{b} \frac{d\mathbf{A}}{d\tau} \right] \quad \Rightarrow \quad (I.7)$$

$$\frac{c}{b} \frac{d\boldsymbol{\pi}}{d\tau} = \left[ e\mathbf{E} + \frac{e}{b} \left( \mathbf{u} \times \mathbf{B} \right) \right] = \left[ e\mathbf{E} + \frac{e}{c} \left( \mathbf{w} \times \mathbf{B} \right) \right] = \frac{d\boldsymbol{\pi}}{dt}.$$

The above shows that the standard and dual equations of motion are mathematically equivalent. (They are clearly not physically equivalent.)

#### C. Dual Maxwell Theory

To study the field of a particle, we write Maxwell's equations (in c.g.s. units):

$$\nabla \cdot \mathbf{B} = 0, \qquad \nabla \cdot \mathbf{E} = 4\pi \rho,$$

$$\nabla \times \mathbf{E} = -\frac{1}{c} \frac{\partial \mathbf{B}}{\partial t}, \quad \nabla \times \mathbf{B} = \frac{1}{c} \left[ \frac{\partial \mathbf{E}}{\partial t} + 4\pi \rho \mathbf{w} \right].$$
(I.8)

Using equations (I.3) and (I.4), we have (the mathematically identical representation):

$$\nabla \cdot \mathbf{B} = 0, \qquad \nabla \cdot \mathbf{E} = 4\pi\rho,$$

$$\nabla \times \mathbf{E} = -\frac{1}{b} \frac{\partial \mathbf{B}}{\partial \tau}, \quad \nabla \times \mathbf{B} = \frac{1}{b} \left[ \frac{\partial \mathbf{E}}{\partial \tau} + 4\pi\rho \mathbf{u} \right].$$
(I.9)

Thus, we obtain a mathematically equivalent set of Maxwell's equations using the local time of the particle to describe its fields.

To derive the corresponding wave equations, we take the curl of the last two equations in (I.9), and use standard vector identities, to get:

$$\frac{1}{b^2} \frac{\partial^2 \mathbf{B}}{\partial \tau^2} - \frac{\mathbf{u} \cdot \mathbf{a}}{b^4} \left[ \frac{\partial \mathbf{B}}{\partial \tau} \right] - \nabla^2 \cdot \mathbf{B} = \frac{1}{b} \left[ 4\pi \nabla \times (\rho \mathbf{u}) \right],$$
(I.10)

$$\frac{1}{b^2} \frac{\partial^2 \mathbf{E}}{\partial \tau^2} - \frac{\mathbf{u} \cdot \mathbf{a}}{b^4} \left[ \frac{\partial \mathbf{E}}{\partial \tau} \right] - \nabla^2 \cdot \mathbf{E} =$$

$$- \nabla (4\pi\rho) - \frac{1}{b} \frac{\partial}{\partial \tau} \left[ \frac{4\pi(\rho \mathbf{u})}{b} \right],$$
(I.11)

where  $\mathbf{a} = d\mathbf{u}/d\tau$  is the particle acceleration. The new term in equation (I.10) is dissipative, acts to oppose the acceleration, is zero when  $\mathbf{a} = 0$  or perpendicular to  $\mathbf{u}$  and arises instantaneously with the force. This makes it clear that the local clock encodes information about the particle's interaction that is unavailable when the clock of the observer or co-moving observer is used to describe the fields. Furthermore, this term does not depend on the nature of the force. This is exactly what one expects of the back reaction caused by inertial resistance of the particle to accelerated motion and, according to Wheeler and Feynman [4], is precisely what is meant by radiation reaction. It follows that no consideration of the action of a particle on itself or the problematic Lorentz-Dirac equation is required to account for radiation reaction.

The **E** and **B** fields can be computed in the standard manner (using only retarded potentials) to get: (see [5])

$$\begin{split} \mathbf{E}\left(\mathbf{x},\tau\right) &= \frac{q\mathbf{r}_{u}\left(1-\mathbf{u}^{2}/b^{2}\right)}{s^{3}} + \frac{q\left[\mathbf{r}\times\left(\mathbf{r}_{u}\times\mathbf{a}\right)\right]}{b^{2}s^{3}} \\ &+ \frac{q\left(\mathbf{u}\cdot\mathbf{a}\right)\left[\mathbf{r}\times\left(\mathbf{u}\times\mathbf{r}\right)\right]}{b^{4}s^{3}} \end{split}$$

and

$$\mathbf{B}\left(\mathbf{x},\tau\right) = \frac{q\left(\mathbf{r}_{u} \times \mathbf{r}\right)\left(1 - \mathbf{u}^{2}/b^{2}\right)}{rs^{3}} + \frac{q\mathbf{r} \times \left[\mathbf{r} \times \left(\mathbf{r}_{u} \times \mathbf{a}\right)\right]}{rb^{2}s^{3}} + \frac{qr\left(\mathbf{u} \cdot \mathbf{a}\right)\left(\mathbf{r} \times \mathbf{u}\right)}{b^{4}s^{3}}.$$

<span id="page-2-0"></span>It is easy to see that **B** is orthogonal to **E**. The last term in each case arises because of the dissipative terms in the respective equation. These terms are zero if a is zero or orthogonal to u. In the first case, there is no radiation and the particle moves with constant velocity so that the field is massless. The second case depends on the creation of motion which keeps  $\mathbf{a}$  orthogonal to  $\mathbf{u}$  (for example a betatron). Since  $\mathbf{r} \times (\mathbf{u} \times \mathbf{r}) = r^2 \mathbf{u} - (\mathbf{u} \cdot \mathbf{r}) \mathbf{r}$ , we see that there is a component along the direction of propagation (longitudinal). (Thus, the E field has a longitudinal part.) This shows that the new dissipative term is equivalent to an effective mass, meaning that the cause for radiation reaction comes directly from the use of the local clock to formulate Maxwell's equations. Thus, there is no need to assume advanced potentials, self-interaction or mass renormalization along with the Lorentz-Dirac equation in order to account for radiation reaction as is done in Dirac's theory. Furthermore, no assumptions about the structure of the source are required.

<span id="page-2-1"></span>Remark I.1. We conjecture that the above effective mass is the actual source of the photoelectric effect and that the photon is a real particle of non-zero (dynamical) mass, which travels with the fields. If this conjecture is correct, radiation from a betatron (of any frequency) exposed to a metal surface will not produce photo electrons.

# II. DUAL RELATIVISTIC QUANTUM THEORY

The Klein-Gordon and Dirac equations were first discovered in early attempts to make quantum mechanics compatible with the Minkowski formulation of the special theory of relativity. Both were partially successful but could no longer be interpreted as particle equations. A complete solution required quantum field theory and its associated problems. In this section we introduce the dual relativistic quantum theory, which always has a single particle theory.

Using equation (I.6), we follow the standard procedure to quantize leading to:

$$i\hbar\frac{\partial\Phi}{\partial\tau}=K\Phi=\left[\frac{H^2}{2mc^2}+\frac{mc^2}{2}\right]\Phi.$$

In addition to the Dirac Hamiltonian, there are two other possible Hamiltonians, depending on the way the potential appears with the square-root operator:

$$\beta \sqrt{c^2 \pi^2 - ec\hbar \Sigma \cdot \mathbf{B} + m^2 c^4} + V$$

and

$$\beta \sqrt{c^2 \pi^2 - ec\hbar \Sigma \cdot \mathbf{B} + (mc^2 + \beta V)^2}$$

This gives us three possible dual relativistic particle equations for spin- $\frac{1}{2}$  particles (see also [7]).

1. The dual Dirac equation:

$$i\hbar \frac{\partial \Psi}{\partial \tau} = \left\{ \frac{\boldsymbol{\pi}^2}{2m} + \boldsymbol{\beta} V + mc^2 - \frac{e\hbar \boldsymbol{\Sigma} \cdot \mathbf{B}}{2mc} + \frac{V\boldsymbol{\alpha} \cdot \boldsymbol{\pi}}{mc} - \frac{i\hbar \boldsymbol{\alpha} \cdot \nabla V}{2mc} + \frac{V^2}{2mc^2} \right\} \Psi.$$
(II.1)

2. The dual version of the square-root equation, using the first possibility:

$$i\hbar \frac{\partial \Psi}{\partial \tau} = \left\{ \frac{\pi^2}{2m} - \frac{e\hbar \mathbf{\Sigma} \cdot \mathbf{B}}{2mc} + mc^2 + \frac{V^2}{2mc^2} \right\} \Psi$$

$$+ \frac{V\beta \sqrt{c^{22} - ec\hbar \mathbf{\Sigma} \cdot \mathbf{B} + m^2c^4}}{2mc^2} \Psi$$

$$+ \frac{\beta \sqrt{c^2 \pi^2 - ec\hbar \mathbf{\Sigma} \cdot \mathbf{B} + m^2c^4}}{2mc^2} V\Psi.$$
(II.2)

3. The dual version of the square-root equation, using the second possibility:

$$i\hbar \frac{\partial \Psi}{\partial \tau} = \left\{ \frac{\pi^2}{2m} + \beta V + mc^2 - \frac{e\hbar \mathbf{\Sigma} \cdot \mathbf{B}}{2mc} + \frac{V^2}{2mc^2} \right\} \Psi.$$
 (II.3)

If  $\mathbf{A}$  and V are zero, all equations reduce to:

$$i\hbar\frac{\partial\Psi}{\partial\tau}=\left\{\frac{\mathbf{p^2}}{2m}+mc^2\right\}\Psi,$$

which is the Schrödinger equation with an added mass term. This makes it easy to see that, in all cases, K is positive definite. In mathematical terms, the lower order terms are relatively bounded with respect to  $\mathbf{p}^2/2m$ . It follows that, unlike the Dirac and Klein-Gordon approach, we can interpret these equations as representations for actual particles. In the above equations, we have assumed that V is time independent. (However, since  $\mathbf{A}(\mathbf{x},\tau)$  can have general time-dependence,  $\sqrt{c^2\pi^2-ec\hbar\Sigma\cdot\mathbf{B}+m^2c^4}$  need not be related to the Dirac operator by a Foldy-Wouthuysen type transformation.)

### III. THE DUAL DIRAC THEORY

We restrict our investigation to the dual Dirac equation. Let  $\mathbf{s}_p$  and  $\boldsymbol{\mu}_p = 2\mu_p \mathbf{s}_p$  be the proton spin and magnetic moment operators respectively. Let  $r_0 = e^2/mc^2$  be the classical electron radius,  $\alpha = \frac{e^2}{\hbar c}$  be the fine structure constant and let  $\boldsymbol{\alpha} = (\alpha_1, \alpha_2, \alpha_3)$  be the standard Dirac matrix, where  $\alpha_i = [\mathbf{0}, \sigma_i, \sigma_i, \mathbf{0}]$ ,

$$\sigma_1 = \begin{pmatrix} \mathbf{0} & 1 \\ 1 & \mathbf{0} \end{pmatrix}, \ \sigma_2 = \begin{pmatrix} \mathbf{0} & -i \\ i & \mathbf{0} \end{pmatrix}, \ \sigma_3 = \begin{pmatrix} \mathbf{1} & \mathbf{0} \\ \mathbf{0} & -\mathbf{1} \end{pmatrix}$$

and  $\sigma = [\sigma_1, \sigma_2, \sigma_3]^t$ . The potentials can be written as  $V_0 = -mc^2r_0/r$ ,  $\mathbf{A} = \boldsymbol{\mu}_p \times \mathbf{r}/r^3$ , where the spin orientation is along the z-axis (i.e.,  $A_r = A_\theta = 0$  and  $A_\phi = \frac{2\mu_p s_p \sin \theta}{r^2}$ ). In what follows,  $\boldsymbol{\pi} = \mathbf{p} - \frac{\mathbf{e}}{\mathbf{c}} \mathbf{A}$  and  $\boldsymbol{\pi}$  is the area of the unit circle.

### A. The Dirac Equation

The eigenvalue problem for the Dirac equation  $\lambda \Psi = H_D \Psi$ , with  $\Psi = [\psi_1, \psi_2]$ , can be written as:

$$(\lambda - V - mc^2)\psi_1 = c(\sigma \cdot \boldsymbol{\pi})\psi_2$$
  

$$(\lambda - V + mc^2)\psi_2 = c(\sigma \cdot \boldsymbol{\pi})\psi_1.$$
(III.1)

Solving the second equation for  $\psi_2$  we have:

<span id="page-3-0"></span>
$$\psi_2 = c \left[ \lambda - V_0 + mc^2 \right]^{-1} (\sigma \cdot \boldsymbol{\pi}) \,\psi_1 \qquad \text{(III.2)}$$

## B. The Dual Dirac Equation

With  $H_D = c\boldsymbol{\alpha} \cdot \boldsymbol{\pi} + mc^2\beta + V_0 = H_0 + V_0$ , let  $V = \frac{1}{2mc^2}[H_0V_0 + V_0H_0]$ . Then, we can write the dual Dirac Hamiltonian as:

$$K_D = \frac{H_D^2}{2mc^2} + \frac{mc^2}{2} = \frac{\pi^2}{2m} + V - \frac{e\hbar\Sigma \cdot \mathbf{B}}{2mc} + mc^2 + \frac{V_0^2}{2mc^2},$$
 (III.3)

## C. The Eigenvalue Problem

The general eigenvalue problem is:

<span id="page-4-0"></span>
$$E\Psi = \left\{ \frac{\pi^2}{2m} + \beta V_0 + mc^2 - \frac{e\hbar\Sigma \cdot \mathbf{B}}{2mc} + \frac{V_0\alpha \cdot \pi}{mc} - \frac{i\hbar\alpha \cdot \nabla V_0}{2mc} + \frac{V_0^2}{2mc^2} \right\} \Psi.$$
(III.4)

Where, as before  $\Psi = [\psi_1, \psi_2]^t$ , with  $\psi_1$ ,  $\psi_2$  the upper and lower spinor components. With  $\mathbf{A} = \mathbf{0}$ , and the exact eigenvalues for  $\lambda \Psi = H_D \Psi$ , we can use  $\left[\frac{\lambda^2}{2mc^2} + \frac{mc^2}{2}\right] \Psi = K_D \Psi$  to find the exact eigenvalues for:

$$\begin{split} E\Psi &= \left\{ \frac{\mathbf{p}^2}{2m} + \beta V_0 + mc^2 + \frac{V_0^2}{2mc^2} \right\} \Psi \\ &+ \left\{ \frac{V_0 \alpha \cdot \mathbf{p}}{2m} - \frac{i\hbar\alpha \cdot \nabla V_0}{2mc} + \frac{V_0 \alpha \cdot \mathbf{p}}{2m} - \frac{i\hbar\alpha \cdot \nabla V_0}{2mc} \right\} \Psi. \end{split}$$

For further analysis, it is convenient to split (III.4) into two equations:

<span id="page-4-1"></span>
$$E\psi_{1} = \left\{ \frac{\boldsymbol{\pi}^{2}}{2m} + V + mc^{2} - \frac{e\hbar\boldsymbol{\sigma}\cdot\mathbf{B}}{2mc} + \frac{V_{0}^{2}}{2mc^{2}} \right\} \psi_{1}$$

$$+ \left\{ \frac{V_{0}\boldsymbol{\sigma}\cdot\boldsymbol{\pi}}{mc} - \frac{i\hbar\boldsymbol{\sigma}\cdot\nabla V_{0}}{2mc} \right\} \psi_{2}$$

$$E\psi_{2} = \left\{ \frac{\boldsymbol{\pi}^{2}}{2m} - V + mc^{2} - \frac{e\hbar\boldsymbol{\sigma}\cdot\mathbf{B}}{2mc} + \frac{V_{0}^{2}}{2mc^{2}} \right\} \psi_{2}$$

$$+ \left\{ \frac{V_{0}\boldsymbol{\sigma}\cdot\boldsymbol{\pi}}{mc} - \frac{i\hbar\boldsymbol{\sigma}\cdot\nabla V_{0}}{2mc} \right\} \psi_{1}.$$
(III.5)

If we now use equation (III.2), with the denominator to the left, we get:

$$\psi_2 = \frac{c\boldsymbol{\sigma} \cdot \boldsymbol{\pi}}{\lambda - V_0 + mc^2} \psi_1.$$

We can now drop the second equation in (III.5) and convert the first to the stationary case and, (using  $\psi$ ) to get the eigenvalue equation:

$$E\psi = \left\{ \frac{\pi^2}{2m} + V_0 + mc^2 - \frac{e\hbar\boldsymbol{\sigma}\cdot\mathbf{B}}{2mc} + \frac{V_0^2}{2mc^2} \right\} \psi + \left\{ \frac{V_0\boldsymbol{\sigma}\cdot\boldsymbol{\pi}}{mc} - \frac{i\hbar\boldsymbol{\sigma}\cdot\nabla V_0}{2mc} \right\} \frac{c\boldsymbol{\sigma}\cdot\boldsymbol{\pi}}{(\lambda - V_0 + mc^2)} \psi.$$

Expanding, we have:

<span id="page-4-2"></span>
$$E\psi = \left\{ \frac{\pi^2}{2m} + V_0 + mc^2 - \frac{e\hbar\boldsymbol{\sigma}\cdot\mathbf{B}}{2mc} + \frac{V_0^2}{2mc^2} \right\} \psi$$

$$-\frac{i\hbar\left(\boldsymbol{\sigma}\cdot\nabla V_0\right)\left(\boldsymbol{\sigma}\cdot\boldsymbol{\pi}\right)}{2m\left(\lambda - V_0 + mc^2\right)} \psi$$

$$+\frac{V_0}{m} \frac{\left(\boldsymbol{\sigma}\cdot\mathbf{p}V_0\right)\left(\boldsymbol{\sigma}\cdot\boldsymbol{\pi}\right)}{\left(\lambda - V_0 + mc^2\right)^2} \psi + \frac{V_0}{m} \frac{\left(\boldsymbol{\sigma}\cdot\boldsymbol{\pi}\right)\left(\boldsymbol{\sigma}\cdot\boldsymbol{\pi}\right)}{\left(\lambda - V_0 + mc^2\right)} \psi.$$
(III.6)

Since the binding energy in Hydrogen is 13ev and the rest mass of the electron is  $5\times10^5{\rm ev}$ , the ratio is  $2.6\times10^5{\rm ev}$ 

 $10^{-5}$ . Thus, there is little loss if we replace  $\lambda - V + mc^2$  by  $2mc^2 + \frac{e^2}{r}$  in equation (III.2). This allows us to bypass the non-linear eigenvalue problem but we must still impose a cut-off since since the denominator is undefined at r=0. With  $r_0=\frac{e^2}{mc^2}$ , we can write (III.2) as:

$$\psi_2 = \frac{c(\boldsymbol{\sigma} \cdot \boldsymbol{\pi})}{2mc^2 \left(1 + \frac{r_0}{2r}\right)} \psi_1. \tag{III.7}$$

Using this, equation (III.6) becomes:

$$E\psi = \left\{ \frac{\pi^2}{2m} + V_0 + mc^2 - \frac{e\hbar\boldsymbol{\sigma}\cdot\mathbf{B}}{2mc} + \frac{V_0^2}{mc^2} \right\} \psi$$

$$-\frac{i\hbar\left(\boldsymbol{\sigma}\cdot\nabla V_0\right)\left(\boldsymbol{\sigma}\cdot\boldsymbol{\pi}\right)}{4m^2c^2\left(1 + \frac{r_0}{2r}\right)} \psi \qquad (III.8)$$

$$+\frac{V_0\left(\boldsymbol{\sigma}\cdot\mathbf{p}V_0\right)\left(\boldsymbol{\sigma}\cdot\boldsymbol{\pi}\right)}{4m^3c^4\left(1 + \frac{r_0}{2r}\right)^2} \psi + \frac{V_0\left(\boldsymbol{\sigma}\cdot\boldsymbol{\pi}\right)\left(\boldsymbol{\sigma}\cdot\boldsymbol{\pi}\right)\left(\boldsymbol{\sigma}\cdot\boldsymbol{\pi}\right)}{2m^2c^2\left(1 + \frac{r_0}{2r}\right)} \psi.$$

The terms inside the first brace are essentially the leading terms for the Schrödinger equation (when  $\mathbf{A}=\mathbf{0}$ ). For proof of concept, we will treat the remaining terms as a first order perturbation.

# 1. The S-state Problem

Our main interest is in the s-state spectra, but before proceeding, we need to calculate the terms which contain  $(\boldsymbol{\sigma} \cdot \boldsymbol{\pi}) (\boldsymbol{\sigma} \cdot \boldsymbol{\pi})$  and  $-i\hbar (\boldsymbol{\sigma} \cdot \boldsymbol{\nabla} V_0) (\boldsymbol{\sigma} \cdot \boldsymbol{\pi})$ . For this, we use the relations

$$(\boldsymbol{\sigma} \cdot \mathbf{X}) (\boldsymbol{\sigma} \cdot \mathbf{Y}) = \mathbf{X} \cdot \mathbf{Y} + i \, \boldsymbol{\sigma} \cdot (\mathbf{X} \times \mathbf{Y}).$$
 (III.9)

If  $\mathbf{X} = \mathbf{Y} = \pi$ , we have

$$(\boldsymbol{\sigma} \cdot \boldsymbol{\pi}) (\boldsymbol{\sigma} \cdot \boldsymbol{\pi}) = \boldsymbol{\pi}^2 + i \, \boldsymbol{\sigma} \cdot (\boldsymbol{\pi} \times \boldsymbol{\pi})$$
$$\boldsymbol{\pi} \times \boldsymbol{\pi} = \frac{ie\hbar}{c} \mathbf{B}$$

$$(\boldsymbol{\sigma} \cdot \boldsymbol{\pi}) (\boldsymbol{\sigma} \cdot \boldsymbol{\pi}) = \boldsymbol{\pi}^2 - \frac{e\hbar}{c} \boldsymbol{\sigma} \cdot \boldsymbol{B}.$$
 (III.10)

If  $\mathbf{X} = -i\hbar \nabla V_0$  and  $\mathbf{Y} = \boldsymbol{\pi}$ , we have:

$$(-i\hbar\boldsymbol{\sigma}\cdot\boldsymbol{\nabla}V_0)(\boldsymbol{\sigma}\cdot\boldsymbol{\pi}) = -i\hbar\nabla V_0\cdot\boldsymbol{\pi} + i\,\boldsymbol{\sigma}\cdot(-i\hbar\boldsymbol{\nabla}V_0\times\boldsymbol{\pi}).$$
(III.11)

By using  $\boldsymbol{\pi} = (\mathbf{p} - \frac{\mathbf{e}}{\mathbf{c}} \mathbf{A})$ , we arrive at

$$(-i\hbar\boldsymbol{\sigma}\cdot\boldsymbol{\sigma}V_{0})(\boldsymbol{\sigma}\cdot\boldsymbol{\pi})$$

$$=-i\hbar\boldsymbol{\sigma}V_{0}\cdot\boldsymbol{\pi}+\hbar\boldsymbol{\sigma}\cdot(\boldsymbol{\sigma}V_{0}\times\boldsymbol{\pi})$$

$$=-i\hbar\boldsymbol{\sigma}V_{0}\cdot\mathbf{p}+\hbar\boldsymbol{\sigma}\cdot(\boldsymbol{\sigma}V_{0})\times\mathbf{p}$$

$$+\frac{ie\hbar}{c}\left[(\boldsymbol{\sigma}V_{0}\cdot\mathbf{A})+i\boldsymbol{\sigma}\cdot(\boldsymbol{\sigma}V_{0}\times\mathbf{A})\right].$$

Since  $\mathbf{A} \propto \mathbf{e}_{\varphi}$ , we see that  $\nabla V_0 \cdot \mathbf{A} \propto \mathbf{e}_r \cdot \mathbf{e}_{\varphi} = 0$ , so that

$$(-i\hbar\boldsymbol{\sigma}\cdot\boldsymbol{\nabla}V_{0})(\boldsymbol{\sigma}\cdot\boldsymbol{\pi}) =$$

$$-i\hbar(\boldsymbol{\nabla}V_{0}\cdot\mathbf{p}) + \hbar\boldsymbol{\sigma}\cdot(\boldsymbol{\nabla}V_{0}\times\mathbf{p}) - \frac{\mathbf{e}\hbar}{c}\boldsymbol{\sigma}\cdot(\boldsymbol{\nabla}V_{0}\times\mathbf{A})$$
(III.12)

If we write  $\mathbf{p}$  in spherical polar coordinates, we have:

$$\mathbf{p} = -i\hbar \nabla = -i\hbar \left( \mathbf{e}_r \frac{\partial}{\partial r} + \frac{1}{r} \mathbf{e}_{\varphi} \frac{\partial}{\partial \varphi} + \frac{1}{\sin \theta} \mathbf{e}_{\theta} \frac{\partial}{\partial \theta} \right). \tag{III.13}$$

It then follows that:

$$\nabla V_0 = \left( \mathbf{e}_r \frac{\partial V_0}{\partial r} + \frac{1}{r} \mathbf{e}_\theta \frac{\partial V_0}{\partial \theta} + \frac{1}{r \sin \theta} \mathbf{e}_\varphi \frac{\partial V_0}{\partial \varphi} \right). \quad \text{(III.14)}$$

Then since  $V_0 = -\frac{e^2}{r}$ , we have

$$-i\hbar\nabla V_0 \cdot \mathbf{p} =$$

$$-i\hbar \left\{ \frac{\partial V_0}{\partial r} \mathbf{e}_r \cdot \left[ -i\hbar \left( \mathbf{e}_r \frac{\partial}{\partial r} \right) \right] \right\} = \frac{e^2 \hbar^2}{r^2} \frac{\partial}{\partial r} \cdot ^{\text{(III.15)}}$$

We also have

$$\hbar \sigma \cdot (\nabla V_0 \times \mathbf{p}) = -\hbar \frac{e^2}{r^3} \sigma \cdot \mathbf{L}$$
 (III.16)

Finnaly, with  $\mathbf{A} = \frac{2\mu_p |s_p| \sin \theta}{r^2} \mathbf{e}_{\varphi}$ , we have

$$-\frac{\mathbf{e}\hbar}{c}\sigma\cdot(\nabla V_0\times\mathbf{A}) = -\frac{\mathbf{e}\hbar}{c}\frac{e^2}{r^4}2\mu_p |s_p|\sin\theta (\sigma\cdot\mathbf{e}_\theta).$$
(III.17)

From these results, we can write the last three terms in equation (III.8) as:

$$\begin{split} &(\mathbf{a}) \quad -\frac{i\hbar\left(\sigma\cdot\nabla V_{0}\right)\left(\sigma\cdot\boldsymbol{\pi}\right)}{2m^{2}c^{2}\left(1+\frac{r_{0}}{2r}\right)} = -\frac{e^{2}\hbar^{2}}{2m^{2}c^{2}\left(1+\frac{r_{0}}{2r}\right)r^{2}}\frac{\partial}{\partial r} \\ &-\frac{e^{2}\hbar}{2m^{2}c^{2}\left(1+\frac{r_{0}}{2r}\right)r^{3}}\boldsymbol{\sigma}\cdot\mathbf{L} - \frac{e^{3}\hbar\mu_{p}\left|\mathbf{s}_{p}\right|}{m^{2}c^{3}\left(1+\frac{r_{0}}{2r}\right)r^{4}}\sin\theta\left(\boldsymbol{\sigma}\cdot\mathbf{e}_{\theta}\right) \\ &(\mathbf{b}) \quad +\frac{V_{0}\left(\boldsymbol{\sigma}\cdot\boldsymbol{\pi}\right)\left(\boldsymbol{\sigma}\cdot\boldsymbol{\pi}\right)}{m^{2}c^{2}\left(1+\frac{r_{0}}{2r}\right)} = -\frac{e^{2}\hbar^{2}\mathbf{p}^{2}}{m^{2}c^{2}r\left(1+\frac{r_{0}}{2r}\right)} \\ &-\frac{e^{3}\hbar^{2}\mathbf{A}^{2}}{m^{2}c^{3}r\left(1+\frac{r_{0}}{2r}\right)} + \frac{e^{3}\hbar\boldsymbol{\sigma}\cdot\mathbf{B}}{m^{2}c^{3}r\left(1+\frac{r_{0}}{2r}\right)} \\ &(\mathbf{c}) \quad +\frac{V_{0}\left(\boldsymbol{\sigma}\cdot\mathbf{p}V_{0}\right)\left(\boldsymbol{\sigma}\cdot\boldsymbol{\pi}\right)}{m^{3}c^{4}\left(1+\frac{r_{0}}{2r}\right)^{2}} = -\frac{e^{4}\hbar^{2}}{m^{3}c^{4}\left(1+\frac{r_{0}}{2r}\right)^{2}r^{3}}\frac{\partial}{\partial r} \\ &-\frac{e^{4}\hbar}{m^{3}c^{4}\left(1+\frac{r_{0}}{2r}\right)^{2}r^{4}}\boldsymbol{\sigma}\cdot\mathbf{L} - \frac{2e^{5}\hbar\mu_{p}\left|\mathbf{s}_{p}\right|}{m^{3}c^{5}\left(1+\frac{r_{0}}{2r}\right)^{2}r^{5}}\sin\theta\left(\boldsymbol{\sigma}\cdot\mathbf{e}_{\theta}\right). \end{split}$$

When  $\mathbf{A} = \mathbf{0}$ , these terms become:

$$(\mathbf{a}') \quad -\frac{i\hbar\left(\sigma\cdot\nabla V_{0}\right)\left(\sigma\cdot\mathbf{p}\right)}{2m^{2}c^{2}\left(1+\frac{r_{0}}{2r}\right)} = -\frac{e^{2}\hbar^{2}}{2m^{2}c^{2}\left(1+\frac{r_{0}}{2r}\right)r^{2}}\frac{\partial}{\partial r}$$

$$-\frac{e^{2}\hbar}{2m^{2}c^{2}\left(1+\frac{r_{0}}{2r}\right)r^{3}}\sigma\cdot\mathbf{L}$$

$$(\mathbf{b}') \quad +\frac{V_{0}\left(\sigma\cdot\mathbf{p}\right)\left(\sigma\cdot\mathbf{p}\right)}{m^{2}c^{2}\left(1+\frac{r_{0}}{2r}\right)} = -\frac{e^{2}\hbar^{2}\mathbf{p}^{2}}{m^{2}c^{2}r\left(1+\frac{r_{0}}{2r}\right)}$$

$$(\mathbf{c}') \quad +\frac{V_{0}\left(\sigma\cdot\mathbf{p}V_{0}\right)\left(\sigma\cdot\mathbf{p}\right)}{m^{3}c^{4}\left(1+\frac{r_{0}}{2r}\right)^{2}} = -\frac{e^{4}\hbar^{2}}{m^{3}c^{4}\left(1+\frac{r_{0}}{2r}\right)^{2}r^{3}}\frac{\partial}{\partial r}$$

$$-\frac{e^{4}\hbar}{m^{3}c^{4}\left(1+\frac{r_{0}}{2r}\right)^{2}r^{4}}\sigma\cdot\mathbf{L}.$$

The new terms that arise, separating equation (III.8) from the Schrödinger equation, and when  $\mathbf{A} \neq \mathbf{0}$  are the two terms from inside first brace of equation (III.8):

$$\frac{2e^2\mu_p^2|\mathbf{s}_p|^2\sin^2\theta}{mc^2r^4} - \frac{e\hbar\sigma \cdot \mathbf{B}}{2mc}$$

and the following three terms from above:

$$-\frac{e^{3}\hbar\mu_{p}|\mathbf{s}_{p}|}{m^{2}c^{3}\left(1+\frac{r_{0}}{2r}\right)r^{4}}\sin\theta\left(\sigma\cdot\mathbf{e}_{\theta}\right)-\frac{4e^{3}\hbar^{2}\mu_{p}^{2}|\mathbf{s}_{p}|^{2}\sin^{2}\theta}{m^{2}c^{3}r^{5}\left(1+\frac{r_{0}}{2r}\right)} + \frac{e^{3}\hbar\sigma\cdot\mathbf{B}}{m^{2}c^{3}r\left(1+\frac{r_{0}}{2r}\right)}-\frac{2e^{5}\hbar\mu_{p}|\mathbf{s}_{p}|}{m^{3}c^{5}\left(1+\frac{r_{0}}{2r}\right)^{2}r^{5}}\sin\theta\left(\sigma\cdot\mathbf{e}_{\theta}\right).$$

Grouping and rearranging the terms, we have:

<span id="page-5-0"></span>
$$\frac{4er_0\hbar\sigma\cdot\mathbf{B}}{2mc(2r+r_0)}-\frac{e\hbar\sigma\cdot\mathbf{B}}{2mc}=-\left[1-\frac{4r_0}{(2r+r_0)}\right]\frac{e\hbar\sigma\cdot\mathbf{B}}{2mc}(\text{III}.18)$$

<span id="page-5-1"></span>
$$\frac{2r_0\mu_p^2|\mathbf{s}_p|^2\sin^2\theta}{r^4} - \frac{4er_0\hbar^2\mu_p^2|\mathbf{s}_p|^2\sin^2\theta}{mcr^5\left(1 + \frac{r_0}{2r}\right)} 
= 2r_0\mu_p^2|\mathbf{s}_p|^2\left[1 - \frac{4e\hbar^2}{mc\left(2r + r_0\right)}\right]\frac{\sin^2\theta}{r^4}$$
(III.19)

and

<span id="page-5-2"></span>
$$-\left[\frac{er_0\hbar\mu_p|\mathbf{s}_p|}{mc(1+\frac{r_0}{2r})r^4} + \frac{2er_0^2\hbar\mu_p|\mathbf{s}_p|}{mc(1+\frac{r_0}{2r})^2r^5}\right]\sin\theta\left(\sigma\cdot\mathbf{e}_\theta\right)$$

$$= -\frac{2er_0\hbar\mu_p|\mathbf{s}_p|}{mc(2r+r_0)}\left[1 + \frac{4r_0}{(2r+r_0)}\right]\frac{\sin\theta}{r^3}\left(\sigma\cdot\mathbf{e}_\theta\right)$$
(III.20)

In the remainder of the paper, we focus on the implications of equation (III.18) and the anomalous magnetic moment. The implications of (III.19) and (III.20) will be a part of another study.

#### D. Anomalous Magnetic Moment

In this section, we investigate the equation (III.18) under the assumption that the charged, spin-1/2 particle does not possess any internal structure (a Dirac particle). In this case, the spin magnetic moment is given by:

$$\mu = g \frac{e}{2mc} \mathbf{s} = g\mu_B \mathbf{s},$$

where  $\mathbf{s} = \frac{\hbar \sigma}{2}$  is the intrinsic spin operator. We can also write the above as

$$H_a = 2\left[1 - \frac{4r_0}{(2r + r_0)}\right] \mu_B \mathbf{s} \cdot \mathbf{B} \qquad \text{(III.21)}$$

Thus, we have that:

$$g_r = 2 \left[ 1 - \frac{4r_0}{(2r + r_0)} \right] \tag{III.22}$$

If we take the cutoff at  $r = \frac{r_0}{2}$ , then g = -2, while if we take the cut off at  $g = \lim_{r \to 0} g_r$ , we obtain g = -6. Taking  $r_e = 0.499857150068631 \times r_0$ , we obtain the correct experimental result:

$$q = -2.00231930436256.$$

If we treat the muon and proton phenomenologically we can also obtain their q-factors:

$$g_{\mu}^{a} = 2 \left[ 1 - \frac{4r_{0}^{\mu}}{(2r_{\mu} + r_{0}^{\mu})} \right]$$

$$g_{p}^{a} = -2 \left[ 1 - \frac{4r_{0}^{p}}{(2r_{p} + r_{0}^{p})} \right],$$

$$e^{2} \qquad 1 \quad r \quad e^{2}$$
(III.23)

where  $r_0^{\mu} = \frac{e^2}{m_{\mu}c^2}$  and  $r_0^p = \frac{e^2}{m_pc^2}$ .

### IV. DISCUSSION

At the classical level we find that the standard and dual theories are mathematically equivalent. At the quantum level, the dual Dirac equation is not mathematically equivalent to the Dirac equation. The dual Dirac equation is strictly positive definite, so that there are no problems with using it as a particle equation. However, we must now directly face the existence of antiparticles.

In order to do this, let us first revisit our conceptual view of the real numbers and their representation. Recall that a field is a set  $\mathbb{A}$  that has two binary operations  $\oplus$ and  $\odot$  that satisfies all our common experience with real numbers. Formally:

**Definition IV.1.** The real numbers is a triplet  $(\mathbb{R}, +, \cdot)$ , which is a field, with 0 as the additive identity (i.e., a +  $0 = a \text{ for all } a \in \mathbb{R}$ ) and 1 as the multiplicative identity (i.e.,  $a \cdot 1 = a$  for all  $a \in \mathbb{R}$ ).

This structure was designed by mathematicians without regard to its possible use in physics. Santilli [8] defined the isodual number field for use in physics and that is what we need.

**Definition IV.2.** The isodual real numbers  $(\hat{\mathbb{R}}, +, *)$  is a field, with  $0 = \hat{0}$  as the additive identity (i.e.,  $\hat{a} + \hat{0} = \hat{a}$ for all  $-a = \hat{a} \in \hat{\mathbb{R}}$ ) and  $\hat{1} = -1$  as the multiplicative identity (i.e.,  $\hat{a} * \hat{1} = (-a)(-1)(-1) = \hat{a}$  for all  $\hat{a} \in \mathbb{R}$ ).

We note that we can obtain the isodual of any physical quantity A from the equation A + A = 0.

<span id="page-6-0"></span>In our theory, the evolution of a particle is formally defined on a Hilbert space  $\mathcal{H}$  over the complex numbers  $\mathbb{C} = \mathbb{R} + i\mathbb{R}$ , with Hamiltonian K by the equation

$$i\hbar\frac{\partial\psi}{\partial\tau}=K\psi.$$
 The conjugate equation is:

$$-i\hbar\frac{\partial\psi^*}{\partial\tau} = K\psi^*.$$

If we use  $\hat{\mathcal{C}}$  as our number field, we can write the above equation as:

$$\hat{i} * \hat{\hbar} * \frac{\partial \psi^*}{\partial \hat{\tau}} = \hat{K} * \psi^*$$

This approach allows us to naturally view anti-particles as particles with their proper time reversed and their evolution defined on  $\mathcal{H}^*$  over  $\mathcal{C}$ . (This does not imply that the time of the observer is reversed.)

Remark IV.3. Santilli [8] has shown that charge conjugation and isoduallity are equivalent for the particleantiparticle symmetry operation.

### CONCLUSION

In this paper we have introduced the dual relativistic quantum theory corresponding to Einstein's special theory of relativity and Maxwell's field theory [2]. The dual classical theory was shown to be mathematically equivalent, but the dual quantum theory is not. We have found three distinct dual relativistic wave equations that reduce to the Schrödinger equation when minimal coupling is turned off. We have focused on the dual Dirac equation and used it to derive a new formula for the gfactor of a spin-1/2 particle. This allowed us to obtain the exact value for the electron g-factor. The formula can also be applied to the muon and the proton. Using the isodual numbers of Santilli [8], we have shown that our theory naturally interprets antiparticles as particles moving backwards in their proper time (and not the time of the observer).

# ACKNOWLEDGMENTS

We would like to thank Professors Netsivi Ben-Amots, Alexander Gersten, Larry Horwitz, Martin C. Land, Elliot Leib and Ruggero M. Santilli for their continued interest, support and suggestions.

- (2002), 69-93.
- <span id="page-7-0"></span>[2] T.L. Gill, and G. Ares de Parga,*The Einstein Dual Theory of Relativity*, Advanced Studies in Theoretical Physics 13(8) (2019) 337-377.
- <span id="page-7-1"></span>[3] P. A. M. Dirac, Proc. Roy. Soc (London) A117 (1928) 610, A118 (1928) 351.
- <span id="page-7-2"></span>[4] J.A. Wheeler and R.P. Feynman, Interaction with the absorber as the mechanism of radiation, *Rev. Mod. Phys.*, 17 (1949), 157-181.
- <span id="page-7-5"></span>[5] T.L. Gill, W.W. Zachary, and J. Lindesay, *The Classical Electron Problem*, Foundations of Phys. 31 (2001) 1299-

- 1354.
- <span id="page-7-3"></span>[6] T. L. Gill and W. W. Zachary, *Two Mathematically Equivalent Versions of Maxwell's Equations*, Foundations of Phys. 41 (2011) 99-128.
- <span id="page-7-4"></span>[7] T.L. Gill, T. Morris and S. K. Kurtz Universal Journal of Physics and Application 3(1)(2015) 24-40.
- <span id="page-7-6"></span>[8] R. M. Santilli, *Isonumbers and genonumbers of dimension 1, 2, 4, 8, their isoduals, and pseudoduals, and "hidden numbers" of dimension 3, 5, 6, 7*, Algebras, Groups and Geometries, Vol. 10, 273-322, 1993.