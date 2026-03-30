## FOUNDATIONS FOR QED I: MATHEMATICAL

## TEPPER L. GILL AND GONZALO ARES DE PARGA

Abstract. This is the first of a series on the foundations for relativistic quantum theory and quantum electrodynamics (QED). In this paper we devote ourselves to the mathematical foundations. We provide proofs of the following:

- (1) The Minkowski Incompatibility Theorem: Minkowski' postulate is incompatible with Einstein's first postulate for two or more particles.
- (2) The Pauli approximation is invalid for s-states and cut-offs are already required before QED.
- (3) The Dirac equation is not physically equivalent to the square root equation but is related by a unitary transformation.
- (4) The Hilbert space designed for Feynman's path integral formulation of quantum theory KS<sup>2</sup> [R n ], is not physically equivalent to L 2 [R n ], the Hilbert space for the Sch¨odinger formulation, but is related by a unitary transformation.

A major contribution of this paper is that mathematical equivalence is not the same as physical equivalence and that unbiased physical analysis and justification is required before acceptance of all mathematical concepts into physical theory. Another major outcome is the construction of Feynman's time ordered operator calculus, and its use to prove the last two remaining conjectures of Dyson concerning the foundations of QED. We also provide the mathematical foundations for Feynman's path integral and use his operator calculus to extend the path integral to allow a larger class of interactions.

## Introduction and Problem

We are ending the first quarter of the 21st century and find ourselves still facing difficult problems left over from the 19th century (see [24, 25]). Of particular concern, we have yet to obtained a relativistic version of Newtonian mechanics compatible with Maxwell's theory and we can only partially claim a one-particle relativistic quantum theory (Dirac's equation requires holes for completion). The replacement of particles by fields has led to results in excellent agreement with experiment. However, the methods used have led to an additional set of problems, which have defied resolution for almost five generations. The purpose of this paper is to provide the mathematical foundations for a final solution. A major conclusion is that mathematical

1

Key words and phrases. Dyson conjectures; Dirac theory, Pauli equation, Feynman calculus, Dyson conjectures;

assumptions must receive careful vetting before acceptance as physically relevant and correct. The second paper will address the classical foundations, while the third will provide the quantum foundations.

- 0.1. Framework and Historical Background. There are five postulates that provide the foundations for all of physical science:
- (1) There is a real (external) physical world, which exists and is independent of our existence.
  - (2) This physical world is as it appears to us in our consciousness.
- (3) This physical world is knowable to us via our senses (including our surrogate instruments).

All of physical science is rooted in these postulates. Our scientific forefathers no longer exist and yet the laws they left us still exist (in some form), so the first postulate is clear. Without the second we can never be sure that our senses give truth, making the third meaningless. Thus, we take these postulates to be obvious.

- 0.1.1. Experimental Physics. Distortions and even illusions about the world can appear in the consciousness of one or more us and not appear to others (the collective). Thus, we need a collective physical filter that allows us to distinguish physical truth from individual or group illusion. We call this filter a scientific experiment.
- Definition 0.1. A scientific experiment is a controlled reproducible experience. The basic postulate of experimental physics is that:
- (4) We can obtain objective reproducible empirical information about the physical world by experimental investigation.
- Definition 0.2. We say that a physical subject is known when all possible cause effect relations (of physical interest) can be a priori predicted under well-defined control conditions.
- 0.1.2. Theoretical Physics. The objective of theoretical physics is to design and construct mathematical representations of the physical world. These representations must describe the cause effect relationships observed in experiment, must be physically and mathematically consistent, using a minimal number of variables and parameters. The basic postulate of theoretical physics is that:
- (5) Mathematics is the correct tool for the design, analysis and certification of the consistency of representations of physical reality.

The Newton-Maxwell problem was faced directly by Abraham, Einstein, Lorentz, Poincar´e, Ritz and other major thinkers of the 1900's. Starting with a complete analysis of the microscopic behavior of electrons and ions, Lorentz was able to explain all the macroscopic properties of optics and electrodynamics (see [37], [38] ). He was the first to obtain the transformations between different observers, showing how their relative velocities

affected their results when in different inertial frames. In his investigations, Poincar´e discovered an error in the derivation and, after correction, named them for Lorentz. In addition, Poincar´e observed that if time is treated as an imaginary coordinate, these Lorentz transformations formed a group. (He also derived the proper time, which was later claimed by Minkowski.)

Einstein's work on the photo-electric effect, give him a different perspective. As noted by Brown [6], he did not believe that Maxwell's theory would survive the existence of photons, but considered the Lorentz transformations fundamental. Thus, he derived them from kinematical arguments. (At that time, Maxwell's theory was still not accepted by all, see [31] and references therein.)

Einstein [14, 15] observed that the constant c appears in Maxwell's equations for all inertial observers. At that time experimental information about the speed of light was meager, being restricted to macroscopic and astronomical studies. Einstein was the first to realized that a formal postulate on the velocity of light was necessary. His proposal was that all physical theories should satisfy the (now well-known) postulates of special relativity:

- (1) The physical laws of nature and the results of all experiments are independent of the particular inertial frame of the observer (in which the experiment is performed).
- (2) The speed of light in empty space is constant and is independent of the motion of the source or receiver (in an inertial frame).

Like others before him, Einstein realized that the idea of absolute space was neither knowable nor necessary. The first postulate abandons this idea. Einstein's second postulate appears to abandon of the idea of absolute time.

Unlike Minkowski, Poincar´e's insight and understanding of the difference between physics and mathematics helped him to resist the temptation to use the "physically unjustified" mathematical observation that time be made a forth geometric coordinate as a (necessary) tool for the representation of physical reality. Minkowski made this leap later with much philosophical fanfare, but lacking any physical justification. Minkowski was a number theorist with few accomplishments of note in physics. In the work of Walters, we find a clear discussion of Minkowski's physics background and his knowledge of Poincar´e's work (see [60]).

Thus, we make explicit Minkowski's unacknowledged third postulate for Einstein's special theory of relativity:

(3) The correct implementation of the first two postulates requires that time be treated as a fourth coordinate, and the relationship between components so constrained as to satisfy the natural invariance induced by the Lorentz group of electrodynamics, (Minkowski space).

The physics community voted in favor of the Lorentz group as the proper transformation theory. Today, many believe that Newtonian theory is covered by Minkowski's version of the special theory. As will be shown later, the failure to conduct a complete (physical) analysis of the third postulate is at the root of all major problems in the classical foundations.

0.2. **Summary.** In the first section we show that there are three representations for the proper time which follow from Einstein's special theory. The first representation is well-known and is used to prove that Minkowski's postulate is incompatible with the first postulate of Einstein for two or more particles. We then discuss the center of mass problem due to Pryce [51], the many particle quantum problem due to the Bakamjian and Thomas [7] and the no-interaction theorem first proven by Currie, Jordan, and Sudarshan [9]. All of these foundational problems are a direct consequence Minkowski's postulate.

In the second section, we used the second representation of proper time and a number of important examples, to explicitly show the danger in equating mathematical equivalence to physical equivalence, and the explicit need for mathematics designed for physics. We prove that the Dirac and square root equations are not physically equivalent, and yet (mathematically) unitary equivalent. As an aside, we provide the correct physical interpretation of the well-known zitterbewegung first observed by Schödinger in the Dirac equation. We end the section with a proof that the Pauli approximation is invalid for the study of s-states and that cut-offs are already required before QED.

In the third section, we introduce the Henstock-Kurzweil integral (HK). It is very close to the one learned in calculus but extends the Lebesgue integral to include nonabsolutely integrable functions (like  $e^{ix^2/2}$ ). We then design the natural Hilbert space for HK integrals, which allows us to construct the Feynman path integral in the manner he intended, while retaining both the physical intuition and computational advantage of his approach. This space  $KS^2[\mathbb{R}^n]$  contains  $L^2[\mathbb{R}^n]$  and the Schwartz distributions or generalized functions as continuous dense embeddings. From a physical point of view, use of the HK-integral and  $KS^2[\mathbb{R}^n]$  eliminates the extra time and effort required to first learn Lebesgue measure theory and then the theory of distributions.

In section four we use the third representation of proper time to prove the existence of a universal clock first predicted by Horwitz and Fanchi, which we call the historical time (see [16] and [34]). We use it to construct Feynman's time ordered operator calculus, construct his path integral, and to prove the last two remaining conjectures of Dyson for the foundations of QED.

## 1. Proper time I and Minkowski's Postulate

Let m be the mass of a particle or the effective mass for a system of particles. We assume the particle or system is defined on phase space with variables  $(t, \mathbf{x}, \mathbf{p})$  as seen by observer O, and variables  $(t', \mathbf{x}', \mathbf{p}')$  as seen by observer O' (both inertial frames). Where  $\mathbf{x}, \mathbf{x}'$  is the position of the particle (or center of mass) and  $(\mathbf{p}, \mathbf{p}')$  is the particle momentum (or center of mass momentum). If  $\gamma^{-1}(\mathbf{v}) = \sqrt{1 - \frac{\mathbf{v}^2}{c^2}}$ ,  $\mathbf{v} = \frac{d\mathbf{x}}{dt}$ , where  $\mathbf{v}(\mathbf{v}')$  is the velocity of the particle, the first representation of the proper time is defined by our observers as:

$$d\tau = \gamma^{-1}(\mathbf{v}) dt$$
,  $d\tau = \gamma^{-1}(\mathbf{v}') dt'$  proper time 1. (1.1)

For the second representation, first square (1.1) to get:

$$d\tau^2 = dt^2 - \frac{1}{c^2} d\mathbf{x}^2, \quad d\tau^2 = dt'^2 - \frac{1}{c^2} d\mathbf{x}'^2.$$

Reaarrange terms to obtain:

$$d\tau^2 + \frac{1}{c^2}d\mathbf{x}^2 = dt^2, \quad d\tau^2 + \frac{1}{c^2}d\mathbf{x'}^2 = dt'^2 \Rightarrow$$

$$\sqrt{1 + \frac{\mathbf{u}^2}{c^2}}d\tau = dt, \quad \mathbf{u} = \frac{d\mathbf{x}}{d\tau}, \quad \sqrt{1 + \frac{\mathbf{u'}^2}{c^2}}d\tau = dt', \quad \mathbf{u'} = \frac{d\mathbf{x'}}{d\tau}.$$

If we let  $b = \sqrt{c^2 + \mathbf{u^2}}$ ,  $(b' = \sqrt{c^2 + \mathbf{u'}^2})$ , then we have our second representation

$$cdt = bd\tau$$
 and  $cdt' = b'd\tau$  proper time 2. (1.2)

For the third representation, use the relations  $H = \gamma(\mathbf{v}) mc^2$  for observer O and  $H' = \gamma(\mathbf{v}') mc^2$  for observer O', to obtain:

$$d\tau = \frac{mc^2}{H}dt$$
,  $d\tau = \frac{mc^2}{H'}dt'$  proper time 3. (1.3)

1.1. **The Minkowski Problem.** Minkowski used the first representation of proper time as the metric for his four-spacetime theory. Einstein was the first to point out that it is not the differential of an exact form as the value of  $\tau$  depended on the forces acting on the particle. In order to use his postulate, Minkowski introduced the clock of a co-moving observer as a substitute (see the notes in Sommerfeld in [45]). This additional postulate which is mathematically convenient but should have alarmed all physics thinkers, since no co-moving observer can actually sit on the tangent plane of a bullet.

The idea of a four dimensional space was very attractive to mathematicians who entered and dominated the field for at least a decade. In the eyes of everyone, this also made Minkowski's postulate a natural part of Einstein's special theory. As a consequence, the intense physical investigation and justification normally given for new postulates did not happen. By the time problems in attempts to develop a many particle relativistic quantum theory appeared, Minkowski's postulate was already sacred. For the general

theory, Einstein only considered an approach that extended Minkowski's postulate (see Pais [43]). The following theorem shows how important an early investigation would have been.

Theorem 1.1. (Minkowski Incompatible Theorem) The addition of Minkowski's postulate to the postulates of Einstein is incompatible for two or more particles,

*Proof.* Let O and O' be two inertial observers. Without loss, we can assume that both clocks begin when their origins coincide and O' is moving with uniform velocity  $\mathbf{v}$  as seen by O. Let two particles, each the source of an electromagnetic field, move with velocities  $\mathbf{w}_i$  (i = 1, 2), as seen by O, and  $\mathbf{w}'_i$  (i = 1, 2), as seen by O', such that:

$$\mathbf{x}'_{i} = \mathbf{x}_{i} - \gamma(\mathbf{v})\mathbf{v}t + \left[\gamma(\mathbf{v}) - 1\right]\left(\mathbf{x}_{i} \cdot \mathbf{v} / \|\mathbf{v}\|^{2}\right)\mathbf{v},$$

$$\mathbf{x}_{i} = \mathbf{x}'_{i} + \gamma(\mathbf{v})\mathbf{v}t' + \left[\gamma(\mathbf{v}) - 1\right]\left(\mathbf{x}'_{i} \cdot \mathbf{v} / \|\mathbf{v}\|^{2}\right)\mathbf{v} \quad \text{and},$$
(1.4)

$$\mathbf{w'}_{i} = \mathbf{w}_{i} - \gamma(\mathbf{v})\mathbf{v} + [\gamma(\mathbf{v}) - 1](\mathbf{w}_{i} \cdot \mathbf{v} / ||\mathbf{v}||^{2})\mathbf{v},$$

$$\mathbf{w}_{i} = \mathbf{w'}_{i} + \gamma(\mathbf{v})\mathbf{v} + [\gamma(\mathbf{v}) - 1](\mathbf{w'}_{i} \cdot \mathbf{v} / ||\mathbf{v}||^{2})\mathbf{v}.$$
(1.5)

Thus, there is no problem in requiring the positions and velocities to transform as expected. However, when we try to transform the clocks, we see the problem at once because

$$t' = \gamma(\mathbf{v}) \left( t - \mathbf{x}_1 \cdot \mathbf{v}/c^2 \right)$$
 and  $t' = \gamma(\mathbf{v}) \left( t - \mathbf{x}_2 \cdot \mathbf{v}/c^2 \right)$ . (1.6)

This is impossible unless  $\mathbf{x}_1 \cdot \mathbf{v} = \mathbf{x}_2 \cdot \mathbf{v}$  (for all time). This means that the two observers cannot use their clocks to share information (with other observers) about two or more particles. It follows that, if we attempt to replace  $\mathbf{x}_i$  and  $\mathbf{x}'_i$  with four vectors using t and t', the first postulate fails. Thus, these three postulates are both mathematically and physically incompatible.  $\square$ 

Remark 1.2. The Lorentz transformations (equations (1.4) (1.5) and (1.6)) contain information about each observer and each particle, so it is not clear that they can be converted into the tensor notation that came into fashion after Sommerfeld simplified the (truly) complicated notation of Minkowski.

1.1.1. The *n*-particle position problem. To construct the many-particle clock for the O frame observer, assume that n interacting particles have Hamiltonians  $H_i$  and a total Hamiltonian  $H = \sum_{i=1}^n H_i$ . The effective mass energy  $Mc^2$ , and total momentum  $\mathbf{P}$ , are defined as follows:

$$Mc^2 = \sqrt{H^2 - c^2 \mathbf{P}^2}, \quad \mathbf{P} = \sum_{i=1}^n \mathbf{p}_i.$$

We can now also represent the Hamiltonian by  $H = \sqrt{c^2 \mathbf{P}^2 + M^2 c^4}$ . In 1948 Pryce conducted the first study of the relativistic center-of-mass for

two or more particles (see [51]). He found only one was canonical and available to all the observers. He observed that the canonical representation of the center-of-mass cannot be the three-vector part of a four-vector. This problem arises because the canonical center of mass position  $\mathbf{X}$  is defined by: (For a modern representation see Longhi, Lusanna, and Pons [40].)

$$\mathbf{X} = \frac{1}{H} \sum_{i=1}^{n} H_i \mathbf{x}_i + \frac{c^2 \left(\mathbf{S} \times \mathbf{P}\right)}{H \left(Mc^2 + H\right)},\tag{1.7}$$

where S is the global spin of the system of particles relative to O. (It is clear that (1.7) cannot represent the vector part of a four-vector.) However, X is the canonical conjugate of P and this is precisely what is needed for a many-particle relativistic quantum theory. (It is clear that the problem is the third postulate.)

Five years after Pryce, Bakamjian and Thomas constructed a many-particle relativistic quantum theory using the canonical center of mass X. They explicitly showed that their theory satisfied the two postulates of Einstein (see [7]). They conjectured that, any attempt to enforce the third postulate would only be compatible with free particles (BT-conjecture).

In 1965 Currie, Jordan and Sudarshan [9], gave a proof of the BT-conjecture for two particles, which they renamed "The No-Interaction Theorem". It has since been extended to an arbitrary number of particles by Leutwyler [39].

**Theorem 1.3.** (No-Interaction Theorem) Suppose we have a system of n particles with phase space variables  $\{(\mathbf{x}_i, \mathbf{p}_i)\}_{i=1}^n$  defined on  $\mathbb{R}^{3n} \times \mathbb{R}^{3n}$  with the following properties:

- (1) The system has a Hamiltonian representation.
- (2) The system has a canonical representation of the Poincaré group.
- (3) Each  $\mathbf{x}_i$  is the vector part of a four-vector.

Then these assumptions are only compatible with free particles.

Thus, we see that Minkowski's postulate is the culprit, which has hindered progress in both classical and quantum mechanics for over a 120 years. The first hint at a solution was given by Feynman when he showed that, by letting time recover its role in physics as it appears in our consciousness, he could construct his operator calculus, his formulation of quantum mechanics and QED.

## 2. Proper Time II and Quantum Problems

In this section, we first show how the second definition of proper time avoids the problems created by Minkowski's postulate. We then analyze the Dirac equation and show why the Pauli approximation cannot be used to study s-states in Hydrogen. We then show that although related by a unitary transformation, the Dirac operator and the square root operator are not physically equivalent.

2.0.1. One-Particle Clock. The proper time for each particle is uniquely defined for each observer. Assume a system of *n*-particles observed by O, who is able to identify each particle and attach a vector  $\mathbf{x}_i$  to the  $i^{\text{th}}$  particle, denoting its spacial distance to the origin. If  $\mathbf{w}_i$  is the velocity of particle  $i^{\text{th}}$  as seen by O, let  $\gamma^{-1}(\mathbf{w}_i) = \sqrt{1 - \mathbf{w}_i^2/c^2}$ . The  $i^{\text{th}}$  particle proper time is defined as before by:

$$d\tau_i = \gamma^{-1}(\mathbf{w}_i)dt, \quad \mathbf{w}_i = \frac{d\mathbf{x}_i}{dt}, \quad d\tau_i^2 = dt^2 - \frac{1}{c^2}d\mathbf{x}_i^2.$$
 (2.1)

Rewrite the last term to get:

$$dt^{2} = d\tau_{i}^{2} + \frac{1}{c^{2}}d\mathbf{x}_{i}^{2}, \Rightarrow cdt = \left(\sqrt{\mathbf{u}_{i}^{2} + c^{2}}\right)d\tau_{i}, \quad \mathbf{u}_{i} = \frac{d\mathbf{x}_{i}}{d\tau_{i}} = \gamma(\mathbf{w}_{i})\mathbf{w}_{i}.$$
(2.2)

If we let  $b_i = \sqrt{\mathbf{u}_i^2 + c^2}$ , the second term in equation (2.2) becomes  $cdt = b_i d\tau_i$ . This leads to a new identity:

$$\frac{1}{c}\frac{d}{dt} \equiv \frac{1}{b_i}\frac{d}{d\tau_i}, i = 1, \dots, n.$$
(2.3)

Thus, we obtain an interesting set of relationships between the proper time of each particle and and that of the observer. Applying it to  $\mathbf{x}_i$ , we obtain another set of relationships showing that the tangent space is left invariant in all cases.

$$\frac{\mathbf{w}_i}{c} = \frac{1}{c} \frac{d\mathbf{x}_i}{dt} \equiv \frac{1}{b_i} \frac{d\mathbf{x}_i}{d\tau_i} = \frac{\mathbf{u}_i}{b_i}, i = 1, \dots, n.$$
 (2.4)

The new particle coordinates are  $(\mathbf{x}_i, \tau_i)$ . In this representation,  $\mathbf{x}_i$  is uniquely defined relative to O, while  $\tau_i$  is uniquely defined by each particle.

Remark 2.1. Equation (2.4) implies that the use of  $\beta = \mathbf{w}_i/c \equiv \mathbf{u}_i/b_i$  as a parameter for measurements is no-longer reliable and that time dilation is a physically incorrect concept.

Since the momentum is  $\mathbf{p}_i = m_i \gamma(\mathbf{w}_i) \mathbf{w}_i$  and  $\mathbf{u}_i = \gamma(\mathbf{w}_i) \mathbf{w}_i$ , we can also represent the momentum as  $\mathbf{p}_i = m_i \mathbf{u}_i$ . It follows that the particle phase space is also left invariant. As  $m_i$  is the rest mass, we see that the concept of relativistic mass increase is also is a physically incorrect concept.

**Theorem 2.2.** If the observer time is replaced by each particle's proper time, the Minkowski incompatibility theorem no-longer holds and the transformations that preserve the first postulate are:

$$b_i' = \gamma (\mathbf{v}) (b_i - \mathbf{u}_i \cdot \mathbf{v}/c)$$

and

$$b_i = \gamma \left( \mathbf{v} \right) \left( b_i' + \mathbf{u'}_i \cdot \mathbf{v}/c \right).$$

*Proof.* If we use each particle proper time then, for observer O'

$$\frac{\mathbf{w}_i'}{c} = \frac{1}{c} \frac{d\mathbf{x}_i'}{dt} \equiv \frac{1}{b_i'} \frac{d\mathbf{x}_i'}{d\tau_i} = \frac{\mathbf{u}_i'}{b_i'}, i = 1, \dots, n.$$

From  $t' = \gamma(\mathbf{v}) (t - \mathbf{x}_i \cdot \mathbf{v}/c^2)$  we have  $cdt' = \gamma(\mathbf{v}) (cdt - d\mathbf{x}_i \cdot \mathbf{v}/c)$ . Using the two identities  $cdt' = b_i'd\tau_i$  and  $cdt = b_id\tau_i$ , we obtain:

$$b_i' d\tau_i = \gamma(\mathbf{v}) \left( b_i d\tau_i - d\mathbf{x}_i \cdot \mathbf{v}/c \right)$$

Dividing by  $d\tau_i$ , we see that the transformations between  $b_i$  and  $b'_i$  are

$$b_i' = \gamma(\mathbf{v}) (b_i - \mathbf{u}_i \cdot \mathbf{v}/c)$$
.

By a similar calculation, we obtain the reverse relationship:

$$b_i = \gamma \left( \mathbf{v} \right) \left( b_i' + \mathbf{u'}_i \cdot \mathbf{v}/c \right).$$

It now follows that the two postulates are valid for all observers.

Each observer can now use their clock and the canonical center of mass to quantize the n-particle system (following [7]).

**Remark 2.3.** We have three important comments:

- (1) If O and O' are at rest relative to each other then  $\mathbf{v} = 0$ ,  $\gamma(\mathbf{v}) = 1$  and  $b_i = b'_i$ .
- (2) If particle i is at rest in O, then  $\tau_i = t$  and  $b_i = c$  but  $\tau_i \neq t'$  and  $b_i' > c$ .
- (3) Since the particle frames need not be inertial,  $b_i > c$  does not violate the second postulate.
- 2.1. Quantum Problems. This section is devoted to the discussion a number of missed and misunderstood mathematical issues in the foundations of relativistic quantum theory.
- 2.1.1. The Dirac Equation. After Palui proposed his exclusion principle in 1925 ([48]), Dirac was able to find the first physically viable way around the square-root problem in 1926 (see [11]). Raised in the British tradition of quaternions, Dirac knew that quaternions could be used to generate a spin-1/2 representation for the rotation group. He used the Pauli matrices to write the square-root operator as:

$$\sqrt{c^2\mathbf{p}^2 + m^2c^4} = \sqrt{\left[c\boldsymbol{\alpha}\cdot\mathbf{p} + \boldsymbol{\beta}mc^2\right]^2}, \quad \text{where,} \quad \boldsymbol{\beta} = \begin{bmatrix} I_2 & 0 \\ 0 & -I_2. \end{bmatrix} \quad (2.5)$$

and the matrix  $\alpha$  is defined by  $\alpha = (\alpha_1, \alpha_2, \alpha_3)$ , with

$$\alpha_i = \begin{pmatrix} \mathbf{0} & \sigma_i \\ \sigma_i & \mathbf{0} \end{pmatrix}, \ \sigma_1 = \begin{pmatrix} \mathbf{0} & 1 \\ 1 & \mathbf{0} \end{pmatrix}, \ \sigma_2 = \begin{pmatrix} \mathbf{0} & -i \\ i & \mathbf{0} \end{pmatrix}, \ \sigma_3 = \begin{pmatrix} 1 & \mathbf{0} \\ \mathbf{0} & -1 \end{pmatrix}$$

He thus obtained an alternative representation of equation (2.5), now known as the Dirac equation:

$$i\hbar \frac{\partial \Psi}{\partial t} = \left[ c\boldsymbol{\alpha} \cdot \mathbf{p} + \boldsymbol{\beta} mc^2 \right] \Psi$$
 (2.6)

We must now view  $\Psi$  as a vector-valued function or spinor. That is,  $\Phi \in L^2(\mathbb{R}^3, C^4) = L^2(\mathbb{R}^3) \otimes C^4$  is a four-component column vector  $\Psi = (\psi_1, \psi_2, \phi_1, \phi_2)^t$ . In this representation,  $\psi = (\psi_1, \psi_2)^t$  is the particle component and  $\phi = (\phi_1, \phi_2)^t$  is the antiparticle component. Dirac's equation gave the correct spin value for the electron and the then known (at that time) spectrum for Hydrogen. His equation also predicted negative energy solutions, which led Dirac to propose his hole theory of the vacuum. In this view, a hole is an electron with negative energy and positive charge. This

means that Dirac's theory cannot be single particle one. However, the later discovery of the positron confirmed his approach. Schweber [59] provides a

good discussion of the difficulties, which eventually led to QED.

A complete understanding of the Dirac equation only became possible after the development of geometric algebra. We now see that Dirac used a Lorentz-invariant Clifford algebra to represent the complex algebra of observables for the electron (see Hestenes [33]). However, it is apparently believed that, when interaction is present a complete separation of the particle and antiparticle components of his equation without approximation is not possible. The algebraic researches of Foldy-Wouthuysen [20], Feynman and Gell-Mann [18] and Pauli [49], along with the various approximations approaches in the literature may have supported this idea.

The first exact analytical separation of Dirac's equation with minimal coupling was found by Gill, Zachary and Alfred in 2005. Since this is important for our physical interpretation of Dirac's equation, the correct physical interpretation of the zitterbewegung, the correct physical relation between Dirac's equation and the square root equation and, to understand the physical failure of the Palui approximation for the study of s-states in Hydrogen, we provide a brief outline of this method. (See [30] for details.)

If  $\mathbf{A}(\mathbf{x},t)$  and  $V(\mathbf{x})$  are vector and scalar potentials make the standard identifications  $\mathbf{p} \to \boldsymbol{\pi} = \mathbf{p} - (e/c)\mathbf{A}$  and write (2.6) in two-component form:

$$i\hbar \frac{\partial \psi}{\partial t} - (mc^2 + V) \psi = c (\boldsymbol{\sigma} \cdot \boldsymbol{\pi}) \varphi,$$
 (2.7)

and

$$i\hbar \frac{\partial \varphi}{\partial t} + (mc^2 - V) \varphi = c (\sigma \cdot \pi) \psi.$$
 (2.8)

Where the vectors  $\boldsymbol{\psi} = [\psi_1, \psi_2]^t$  and  $\boldsymbol{\varphi} = [\varphi_1, \varphi_2]^t$  represent the particle and anti-particle components  $\boldsymbol{\Psi} = [\boldsymbol{\psi}, \boldsymbol{\varphi}]^t$ . We can separate these two equations using the method of integrating factors as in ODE's. We solve for  $\boldsymbol{\varphi}$  in the second equation and use the result in the first to obtain an equation which only depends on  $\boldsymbol{\psi}$ . Thus, we need to determine  $\delta$  such that:

$$e^{-\delta t}i\hbar\frac{\partial}{\partial t}\left[\varphi e^{\delta t}\right]=i\hbar\frac{\partial\varphi}{\partial t}+\left(mc^2-V\right)\varphi.$$

Simplify the left-hand side and compare to obtain that  $\delta = \frac{1}{i\hbar} (mc^2 - V)$ . We now assume that  $\psi$  and  $\varphi$  vanish at  $-\infty$  and integrate to obtain:

$$i\hbar \frac{\partial \psi}{\partial t} - (mc^2 + V) \psi = c^2 (\boldsymbol{\sigma} \cdot \boldsymbol{\pi}) \left\{ \int_{-\infty}^t \exp\left\{\frac{i}{\hbar} (mc^2 - V) (t - \tau)\right\} (\boldsymbol{\sigma} \cdot \boldsymbol{\pi}) \psi(\mathbf{x}, \tau) d\tau \right\}.$$
(2.9)

Using the same method, we obtain our second equation for  $\varphi$ :

$$i\hbar \frac{\partial \varphi}{\partial t} + (mc^2 - V) \varphi = c^2 (\sigma \cdot \pi) \left\{ \int_{-\infty}^t \exp\left\{ -\frac{i}{\hbar} (mc^2 + V) (t - \tau) \right\} (\sigma \cdot \pi) \varphi(\mathbf{x}, \tau) d\tau \right\}.$$
(2.10)

Equations (2.9) and (2.10) now decomposed  $L^2(\mathbb{R}^3, C^4)$  into a direct sum  $L^2(\mathbb{R}^3, C^4) = L^2(\mathbb{R}^3, C^2) \oplus L^2(\mathbb{R}^3, C^2)$ . For later use, the probability density for  $\psi$  is:

$$\rho_{\boldsymbol{\psi}}(t) = \left| \boldsymbol{\psi}\left(\mathbf{x}, t\right) \right|^{2} + \left| c \int_{-\infty}^{t} \exp\left\{ \frac{i}{\hbar} \left( mc^{2} - V \right) (t - \tau) \right\} (\boldsymbol{\sigma} \cdot \boldsymbol{\pi}) \, \boldsymbol{\psi}\left(\mathbf{x}, \tau\right) d\tau \right|^{2} (2.11)$$

Thus, we obtain a complete analytic diagonalization of the particle, antiparticle wave functions (see [30]). The decomposition is analytic, compared to Foldy and Wouthuysen's group theoretical framework in that the wave function is not transformed (see [20]). We will discuss this in further detail, after we study the square root operator.

2.1.2. The Pauli Approximation. For the eigenvalue problem, equations (2.7) and (2.8) can be written as:

$$(E - V - mc^{2}) \psi = c (\sigma \cdot \pi) \varphi,$$
  

$$(E - V + mc^{2}) \varphi = c (\sigma \cdot \pi) \psi.$$

Using

$$\varphi = c(E - V + mc^2)^{-1} (\boldsymbol{\sigma} \cdot \boldsymbol{\pi}) \psi$$

in the first equation and

$$\psi = c(E - V - mc^2)^{-1} (\boldsymbol{\sigma} \cdot \boldsymbol{\pi}) \varphi$$

in the second equation. After simplification, we obtain:

$$(E - V - mc^{2}) \psi = \frac{c^{2} (\boldsymbol{\sigma} \cdot \mathbf{p}V) (\boldsymbol{\sigma} \cdot \boldsymbol{\pi})}{(E - V + mc^{2})^{2}} \psi + \frac{c^{2} (\boldsymbol{\sigma} \cdot \boldsymbol{\pi}) (\boldsymbol{\sigma} \cdot \boldsymbol{\pi})}{(E - V + mc^{2})} \psi$$
(2.12)

and

$$(E - V + mc^{2})\varphi = \frac{c^{2}(\sigma \cdot \mathbf{p}V)(\sigma \cdot \pi)}{(E - V - mc^{2})^{2}}\varphi + \frac{c^{2}(\sigma \cdot \pi)(\sigma \cdot \pi)}{(E - V - mc^{2})}\varphi. \quad (2.13)$$

Equations (2.12) and (2.13) first appeared in a paper by Williams in 1940 (see [62]). He was a student of Slater, who included them in the appendix of his book published in 1960 (see [58]). We call them the Williams-Slater

equations. As will be seen, they are important for understanding the mathematical foundations of QDE. In equation (2.12), replace  $(E - V + mc^2)^{-2}$ by 0 and approximate  $E-V+mc^2$  by  $2mc^2$ , to obtain the Pauli equation:

$$(E - V - mc^{2}) \psi = \frac{(\boldsymbol{\sigma} \cdot \boldsymbol{\pi}) (\boldsymbol{\sigma} \cdot \boldsymbol{\pi})}{2m} \psi.$$

The probability density for this equation now becomes:

$$\rho_{\boldsymbol{\psi}} = \left| \boldsymbol{\psi} \left( \mathbf{x} \right) \right|^2 + \left| \frac{\left( \boldsymbol{\sigma} \cdot \boldsymbol{\pi} \right)}{2mc} \boldsymbol{\psi} \left( \mathbf{x} \right) \right|^2$$

If we use

$$(\boldsymbol{\sigma} \cdot \boldsymbol{\pi}) (\boldsymbol{\sigma} \cdot \boldsymbol{\pi}) = \boldsymbol{\pi}^2 - \frac{e\hbar}{c} \boldsymbol{\sigma} \cdot \mathbf{B},$$

We obtain the standard form:

$$E\psi = \left[\frac{\pi^2}{2m} - \frac{e\hbar (\boldsymbol{\sigma} \cdot \mathbf{B})}{2mc} + V + mc^2\right]\psi$$
 (2.14)

In the case of interest,  $V_0 = -mc^2r_0/r$ ,  $\mathbf{A} = \boldsymbol{\mu}_p \times \mathbf{r}/r^3$ , with the spin orientated along the z-axis, so that  $A_r = A_\theta = 0$  and  $A_\phi = \frac{2\mu_p s_p \sin \theta}{r^2}$ . Thus, we see that  $\pi^2 = \mathbf{p^2} + \mathbf{A^2}$ . We recall that QED provides corrections to the Dirac equation for s-states in Hydrogen based on the Pauli approximation. Thus, it is important to know the extent that we can believe any conclusions from the Pauli equation about the behavior of the Dirac equation in sstates. Recalling that a Dirac electron has a finite probability of being at the origin in the s-state, we se that the following make the Pauli approximation invalided:

- (1)  $(E-V+mc^2)$  cannot be approximated by  $2mc^2$  for s-states. (2) As **A** has a  $r^{-2}$  term and **A**<sup>2</sup> has a  $r^{-4}$  term they cannot omitted for s-state calculations.
- (3) The use of  $\rho_{\psi} = |\psi(\mathbf{x})|^2$  as the density is incorrect and the corrected density  $\rho_{\psi} = |\psi(\mathbf{x})|^2 + \left| \frac{(\boldsymbol{\sigma} \cdot \boldsymbol{\pi})}{2mc} \psi(\mathbf{x}) \right|^2$  cannot be used because of (2).

To resolve these difficulties, we note that the binding energy in Hydrogen is 13ev, while the rest mass of the electron is  $5 \times 10^5$  ev, so the ratio is  $2.6 \times 10^{-5}$ . Thus, we can replace  $E - V + mc^2$  by  $2mc^2 + \frac{e^2}{r}$  with little error, but a cut-off is required at the classical electron radius (see [30] for details). The important conclusion is that any analysis of the Dirac Hydrogen atom problem for s-states requires a cut-off, prior to any QED calculations.

2.1.3. The Square Root Equation. In 2005 Gill and Zachary [29] used the theory of fractional powers of linear operators to construct an (analytic) representation for the square-root energy operator valid for all spin values.

A better representation was found by Samko (see [55] pg. 270):

$$\beta \hbar c \left[ \sqrt{\omega^{2} I - \Delta} \right] \psi (\mathbf{x}, t) = \beta m c^{2} \psi (\mathbf{x}, t)$$

$$- \frac{\beta \hbar c}{2\pi^{2}} \left\{ \int_{\mathbb{R}^{3}} \frac{K_{2} \left[ \omega \left| \mathbf{x} - \mathbf{y} \right| \right]}{\left| \mathbf{x} - \mathbf{y} \right|^{2}} \left[ \psi (\mathbf{x}, t) - \psi (\mathbf{y}, t) \right] d\mathbf{y} \right\}.$$
(2.15)

The integral is over  $\mathbb{R}^3$  showing in what sense the square root operator is nonlocal in space. For a better understanding of this extended object, we need more information about the Bessel functions  $K_2[\omega|\mathbf{y}|]$ ,  $K_1[\omega|\mathbf{y}|]$ ,  $K_{\frac{1}{2}}[\omega|\mathbf{y}|]$  and  $K_0[\omega|\mathbf{y}|]$  (see Gradshteyn and Ryzhik [26]).

$$K_{2} \left[\omega \left|\mathbf{y}\right|\right] = K_{0} \left[\omega \left|\mathbf{y}\right|\right] + \frac{K_{1} \left[\omega \left|\mathbf{y}\right|\right]}{\left[\omega \left|\mathbf{y}\right|\right]} \quad \text{and}$$

$$\frac{K_{1/2} \left[\omega \left|\mathbf{y}\right|\right]}{\left[\omega \left|\mathbf{y}\right|\right]^{1/2}} = \sqrt{\frac{\pi}{2}} \frac{exp\left\{-\omega \left|\mathbf{y}\right|\right\}}{\left[\omega \left|\mathbf{y}\right|\right]}.$$

Since the integral of  $|\mathbf{x} - \mathbf{y}|^{-2}$  is finite, the effective singularity properties for  $0 < z \ll 1$  are:

$$\frac{K_{1}[z]}{z} = C_{1} [1 + \theta_{1}(z)] z^{-2} 
\frac{K_{1/2}[z]}{z^{1/2}} = \sqrt{\frac{\pi}{2}} [1 + \theta_{2}(z)] z^{-1} 
K_{0}[z] = C_{0} [1 + \theta_{0}(z)] \ln z^{-1}$$

The  $K_1$  diverges at 0 like  $z^{-2}$  (strong interaction?) and the  $K_0$  is actually an integrable singularity. We include the  $K_{1/2}$  for comparison as it is essentially the well-known Yakawa potential (expected to account for the short range nuclear interaction). The effective properties for  $z \gg 1$  are:

$$\frac{K_{1}[z]}{z} = C_{1} [1 + \Theta_{1}(z)] z^{-3/2} e^{-z} 
\frac{K_{1/2}[z]}{z^{1/2}} = \sqrt{\frac{\pi}{2}} [1 + \Theta_{1/2}(z)] z^{-1} e^{-z} 
K_{0} [z] = C_{0} [1 + \Theta_{0}(z)] z^{-1/2} e^{-z}$$

In this region, we have a cut-off range and the strength of interaction reverses, with the  $K_0$  having the longest tail (weak interaction?). It should be clear that the square root operator is nonlocal in space and can be treated like the identity outside a few Compton wavelengths.

2.1.4. Discussion I. The important conclusion from the last two sections is that the Dirac and square root operators are not physically equivalent. The Dirac operator is nonlocal in time while the square root operator is nonlocal in space. It is shown in [30] that the correct physical interpretation of the zitterbewegung and the fact that the expected instantaneous value of a velocity measurement of a Dirac particle is  $\pm c$  is that it instantaneously jumps between +c and -c at each point in time to make it appear as a point in space.

We must now account for the known fact that the Foldy-Wouthuysen is a unitary transformation mapping the coupled particle-antiparticle Dirac operator to the uncoupled particle-antiparticle square root operator. From equations (2.9) and (2.10), we see that the Foldy-Wouthuysen unitary operator maps a nonlocal in time operator to a nonlocal in space operator, so that the unitary property of an operator is not sufficient to give confidence in the physics. We should also remember that the concept of unitary equivalence is a purely mathematical one, taken over to physics without any physical justification.

Remark 2.4. The Schödinger representation (on  $L^2[\mathbb{R}^3]$ ) and the Heisenberg representation  $\ell_2[\mathbb{R}^3]$  were known to be physically equivalent before they were proven to also be mathematical equivalent (i.e, related by a unitary mapping) by von Neumann. However, there are an infinite number of Hilbert spaces related to  $\ell_2[\mathbb{R}^3]$  by a unitary mapping, which have nothing to due with quantum theory. In a physically justified sense,  $L^2[\mathbb{R}^3]$  is not the correct Hilbert space for Feynman's path integral formulation of quantum mechanics. We will construct it in the next section. This space will also serve as another counter example to the notion that mathematical equivalence means physical equivalence.

# 3. Hilbert Space for Feynman's Formulation of Quantum Mechanics

In this section, we provide an introduction to the Henstock-Kurzweil integral (HK). The integral is well defined for operator-valued functions that may not be separably valued (where both the Bochner and Pettis integrals are undefined). Intuitively, one uses a version of the Riemann integral (of calculus) with the interior points chosen first, while the size of the base rectangle around any interior point is determined by an arbitrary positive function defined at that point. This integral was discovered independently by Henstock [HS] and Kurzweil [KW]. Our aim is to provide a sense of the conceptual and technical simplicity of the HK-integral. Those interested in more detail and proofs are directed to Gill and Zachary [28]. Let  $\mathcal{H}$  be a separable Hilbert space and let  $L(\mathcal{H})$  be the algebra of bounded linear operators on  $\mathcal{H}$ . Let  $[a,b] \subset \mathbb{R}$  and, for each  $t \in [a,b]$ , let  $H(t) \in L(\mathcal{H})$  be a given family of operators.

**Definition 3.1.** Let  $\delta(t)$  map  $[a,b] \rightarrow (0,\infty)$ , and let  $\mathbf{P} = \{t_0, \tau_1, t_1, \tau_2, \cdots, \tau_n, t_n\}$ , where  $a = t_0 \leqslant \tau_1 \leqslant t_1 \leqslant \cdots \leqslant \tau_n \leqslant t_n = b$ . We call  $\mathbf{P}$  a HK-partition for  $\delta$  (or HK-partition when  $\delta$  is understood) provided that for  $0 \leqslant i \leqslant n-1$ ,  $t_i, t_{i+1} \in (\tau_{i+1} - \delta(\tau_{i+1}), \tau_{i+1} + \delta(\tau_{i+1}))$ .

**Definition 3.2.** The family H(t),  $t \in [a,b]$ , is said to have a (uniform) HK-integral if there is an operator Q[a,b] in  $L(\mathcal{H})$  such that, for each  $\varepsilon > 0$ , there exists a function  $\delta$  from  $[a,b] \to (0,\infty)$  such that, whenever  $\mathbf{P}$  is a HK-partition for  $\delta$ , then

$$\left\| \sum_{i=1}^{n} \Delta t_i H(\tau_i) - Q[a, b] \right\| < \varepsilon.$$

In this case, we write

$$Q[a,b] = (HK) \int_{a}^{b} H(t)dt.$$

**Example 3.3.** The following example shows there are HK-integrable functions, which are not Lebesgue integrable. If  $F'(t) = 2t\sin(1/t^2) - (2/t)\cos(1/t^2)$ , for  $t \in (0,1)$  and F'(0) = 0. It is easy to see that  $F(t) = t^2\sin(1/t^2)$ , so that

$$(HK)$$
  $\int_0^1 \left(2t\sin\frac{1}{t^2} - 2\frac{1}{t}\cos\frac{1}{t^2}\right) dt = \sin 1.$ 

The following theorem explains the difference between the two integrals.

**Theorem 3.4.** Let  $f(t): [a,b] \to \mathbb{R}$ .

- (1) If f(t) is L-integrable on [a,b], then it is HK-integrable on [a,b] and:  $HK-\int_a^b f(t)dt = L-\int_a^b f(t)dt$ .
- (2) If f(t) is HK-integrable on [a,b], then  $\sup_{t>a} \left| \int_a^t f(s)ds \right| < \infty$ .
- (3) If f(t) is HK-integrable and bounded on [a, b], then it is L-integrable on [a, b].

## 3.1. The KS-Hilbert Space.

In this section, our objective is to construct a particular (separable) Hilbert space  $KS^2[\mathbb{R}^n]$ . This space is of special interest because it contains the class of HK-integrable functions, the space of the test functions  $\mathcal{D}[\mathbb{R}^n]$ , the space  $\mathfrak{M}[\mathbb{R}^n]$  of finitely additive measures on  $\mathbb{R}^n$  and  $L^p[\mathbb{R}^n]$  for  $1 \leq p \leq \infty$ . Furthermore, each of the latter class spaces is contained in  $KS^2[\mathbb{R}^n]$  as continuous dense and compact embeddings (e.g., weakly convergent sequences in each of the above spaces are strongly convergent in  $KS^2[\mathbb{R}^n]$ ). In particular,  $KS^2[\mathbb{R}^n]$  is perfect for the construction of Feynman's path integral formulation of quantum mechanics (in the form he suggested).

To construct  $KS^2[\mathbb{R}^n]$ , fix n, and let  $\mathbb{Q}^n$  be the set all vectors  $\{\mathbf{x} = (x_1, x_2 \cdots, x_n) \in \mathbb{R}^n\}$  such that  $x_i$  is rational for each i. This is a countable dense set in  $\mathbb{R}^n$ , so we can arrange it as  $\mathbb{Q}^n = \{\mathbf{x}_1, \mathbf{x}_2, \mathbf{x}_3 \cdots\}$ . For each l and i, let  $\mathbf{B}_l(\mathbf{x}_i)$  be the closed cube centered at  $\mathbf{x}_i$  of side  $r_l = 2^{-l}, l \in \mathbb{N}$  parallel to the coordinate axis. Now choose an order so that the set  $\{\mathbf{B}_k(\mathbf{x}_k), k \in \mathbb{N}\}$  contains all closed cubes  $\{\mathbf{B}_l(\mathbf{x}_i) \mid (l,i) \in \mathbb{N} \times \mathbb{N}\}$  centered at a point in  $\mathbb{Q}^n$ . Let  $\mathcal{E}_k(\mathbf{x}) = 1$  if  $\mathbf{x} \in \mathbf{B}_k(\mathbf{x}_k)$  and be zero otherwise, so that  $\mathcal{E}_k(\mathbf{x})$  is in  $L^p[\mathbb{R}^n]$  for  $1 \leq p \leq \infty$ . Define  $F_k(\cdot)$  on  $L^1[\mathbb{R}^n]$  by

$$F_k(f) = \int_{\mathbb{R}^n} \mathcal{E}_k(\mathbf{x}) f(\mathbf{x}) d\mathbf{x}.$$
 (3.1)

It is clear that  $F_k(\cdot)$  is a bounded linear functional on  $L^p[\mathbb{R}^n]$  for each k,  $||F_k||_{\infty} \leq 1$  and if  $F_k(f) = 0$  for all k, f = 0 so that  $\{F_k\}$  is fundamental

on  $L^p[\mathbb{R}^n]$  for  $1 \leqslant p \leqslant \infty$ . Fix  $\lambda$ , set  $t_{\lambda}^k = \lambda^{k-1} e^{-\lambda} / (k-1)!$  and define a measure  $d\mathbf{P}_{\lambda}(\mathbf{x}, \mathbf{y})$  on  $\mathbf{R}^n \times \mathbf{R}^n$  by:

$$d\mathbf{P}_{\lambda}(\mathbf{x}, \mathbf{y}) = \left[ \sum\nolimits_{k=1}^{\infty} t_{\lambda}^{k} \mathcal{E}_{k}(\mathbf{x}) \mathcal{E}_{k}(\mathbf{y}) \right] d\mathbf{x} d\mathbf{y}.$$

We can now define an inner product  $(\cdot)$  on  $\mathbf{L}^1[\mathbf{R}^n]$  by

$$(f,g) = \int_{\mathbb{R}^n \times \mathbb{R}^n} f(\mathbf{x}) g(\mathbf{y})^* d\mathbf{P}_{\lambda}(\mathbf{x}, \mathbf{y})$$

$$= \sum_{k=1}^{\infty} t_{\lambda}^k \left[ \int_{\mathbb{R}^n} \mathcal{E}_k(\mathbf{x}) f(\mathbf{x}) d\mathbf{x} \right] \left[ \int_{\mathbb{R}^n} \mathcal{E}_k(\mathbf{y}) g(\mathbf{y}) d\mathbf{y} \right]^*.$$
(3.2)

We call the completion of  $L^1[\mathbb{R}^n]$ , with the above inner product, the Kuelbs-Steadman space  $KS^2[\mathbb{R}^n]$ . To see that it contains the HK-integrable functions, let f be HK-integrable on  $\mathbb{R}^n$  then:

$$||f||_{KS}^2 = \sum_{k=1}^{\infty} t_{\lambda}^k \left| \int_{\mathbb{R}^n} \mathcal{E}_k(\mathbf{x}) f(\mathbf{x}) d\mathbf{x} \right|^2 \leqslant \sup_{k} \left| \int_{\mathbb{R}^n} \mathcal{E}_k(\mathbf{x}) f(\mathbf{x}) d\mathbf{x} \right|^2 \leqslant \left| \int_{\mathbb{R}^n} f(\mathbf{x}) d\mathbf{x} \right|^2 < \infty,$$
 so  $f \in KS^2[\mathbb{R}^n]$ .

**Definition 3.5.** We say that  $\mathcal{B} \subset \mathcal{H}$  is a continuous dense embedding if:

- (1)  $\mathcal{B}$  is dense in  $\mathcal{H}$ .
- (2) For each  $f \in \mathcal{B}$ ,  $||f||_{\mathcal{H}} \leq ||f||_{\mathcal{B}}$ .

**Theorem 3.6.** For each  $p, 1 \leq p \leq \infty$ ,  $L^p[\mathbb{R}^n] \subset KS^2[\mathbb{R}^n]$  is a continuous dense embedding.

*Proof.* By construction,  $KS^2[\mathbb{R}^n]$  contains  $L^1[\mathbb{R}^n]$  a continuous dense embedding, so we need only show that  $KS^2[\mathbb{R}^n] \supset L^p[\mathbb{R}^n]$  for  $p \neq 1$ . If  $p < \infty, L^p[\mathbb{R}^n] \cap L^1[\mathbb{R}^n]$  is dense in  $L^1[\mathbb{R}^n]$  and hence dense in  $KS^2[\mathbb{R}^n]$ . If  $f \in L^p[\mathbb{R}^n]$  then:

$$||f||_{KS^{2}} = \left[ \sum_{k=1}^{\infty} t_{\lambda}^{k} \left| \int_{\mathbb{R}^{n}} \mathcal{E}_{k}(\mathbf{x}) f(\mathbf{x}) d\mathbf{x} \right|^{\frac{2p}{p}} \right]^{1/2}$$

$$\leq \left[ \sum_{k=1}^{\infty} t_{\lambda}^{k} \left( \int_{\mathbb{R}^{n}} \mathcal{E}_{k}(\mathbf{x}) |f(\mathbf{x})|^{p} d\mathbf{x} \right)^{\frac{2}{p}} \right]^{1/2}$$

$$\leq \sup_{k} \left[ \int_{\mathbf{R}^{n}} \mathcal{E}_{k}(\mathbf{x}) |f(\mathbf{x})|^{p} d\mathbf{x} \right]^{\frac{1}{p}} \leq ||f||_{p}.$$

Hence,  $f \in KS^2[\mathbb{R}^n]$  and the embedding is continuous. For  $p = \infty$ , we have

$$||f||_{KS^2} = \left[ \sum_{k=1}^{\infty} t_{\lambda}^k \left| \int_{\mathbb{R}^n} \mathcal{E}_k(\mathbf{x}) f(\mathbf{x}) d\mathbf{x} \right|^2 \right]^{1/2}$$

$$\leq \left[ \left[ \sum_{k=1}^{\infty} t_{\lambda}^k [vol(\mathbf{B}_k)]^2 \right] [ess \sup |f|]^2 \right]^{1/2} \leq M ||f||_{\infty}.$$

Thus  $f \in KS^2[\mathbb{R}^n]$ , and  $L^{\infty}[\mathbb{R}^n] \subset KS^2[\mathbb{R}^n]$  as a continuous dense embedding.

We can also construct  $KS^p[\mathbb{R}^n]$  for all p by defining

$$||f||_{KS^p} = \begin{cases} \left\{ \sum_{k=1}^{\infty} t_k \left| \int_{\mathbb{R}^n} \mathcal{E}_k(\mathbf{x}) f(\mathbf{x}) d\mathbf{x} \right|^p \right\}^{1/p}, 1 \leq p < \infty, \\ \sup_{k \geq 1} \left| \int_{\mathbb{R}^n} \mathcal{E}_k(\mathbf{x}) f(\mathbf{x}) d\mathbf{x} \right|, p = \infty. \end{cases}$$

It is easy to see that  $\|\cdot\|_{KS^p}$  defines a norm on  $L^p[\mathbb{R}^n]$ . If  $KS^p[\mathbb{R}^n]$  is the completion of  $L^p[\mathbb{R}^n]$  with respect to this norm, we have:

**Theorem 3.7.** For each  $q, 1 \leq q \leq \infty$ ,  $KS^p[\mathbb{R}^n] \supset L^q[\mathbb{R}^n]$  as a continuous dense embedding and,  $KS^{\infty}[\mathbb{R}^n] \subset KS^p[\mathbb{R}^n]$  for each p.

Remark 3.8. We observe that  $KS^2[\mathbb{R}^n] = KS^2[\mathbb{R}^n]^* \supset L^1[\mathbb{R}^n]^{**} = \mathfrak{M}[\mathbb{R}^n]$ , the space of finitely additive measures on  $\mathbb{R}^n$ , and in particular, the Dirac measure  $\delta(\mathbf{x}) \in \mathfrak{M}[\mathbb{R}^n]$ . The theory of distributions was developed via topological vector space methods because of the belief that there was no Hilbert space which contained the test functions  $\mathcal{D}[\mathbb{R}^n]$ . The next result shows their mistake was due to a lack of interest in non-absolutely integrable functions.

**Definition 3.9.** The space of test functions  $\mathcal{D}[\mathbb{R}^n] = C_c^{\infty}[\mathbb{R}^n]$  (the set of continuous functions with the sequential compact support topology). A sequence  $\phi_j \in \mathcal{D}[\mathbb{R}^n]$  converges to  $\phi \in \mathcal{D}[\mathbb{R}^n]$  in the sequential compact support topology, if there exists a compact set  $K \subset \mathbb{R}^n$  containing the support of  $\phi_j - \phi$  for all j and, for every multi-index  $\alpha$ ,  $D^{\alpha}\phi_j$  converges to  $D^{\alpha}\phi$  uniformly on K.

**Theorem 3.10.** The test functions  $\mathcal{D}[\mathbb{R}^n] \subset KS^2[\mathbb{R}^n]$  as a continuous dense embedding.

Proof. Since  $KS^{\infty}(\mathbb{R}^n)$  is continuously embedded in  $KS^2[\mathbb{R}^n]$ , it suffices to prove the result for  $KS^{\infty}(\mathbb{R}^n)$ . Suppose that  $\phi_j \to \phi$  in  $\mathcal{D}[\mathbb{R}^n]$ , so that there exists a compact set  $K \subset \mathbb{R}^n$ , containing the support of  $\phi_j - \phi$  and  $D^{\alpha}\phi_j$  converges to  $D^{\alpha}\phi$  uniformly on K for every multi-index  $\alpha$ . Let  $L = \{l \in \mathbb{N} : \text{the support of } \mathcal{E}_l, stp\{\mathcal{E}_l\} \cap K \neq \emptyset\}$  and  $M = \sup_{l \in L} [vol(\mathbf{B}_l)]$ , then

$$\lim_{j \to \infty} \left\| D^{\alpha} \phi - D^{\alpha} \phi_{j} \right\|_{KS^{\infty}} = \lim_{j \to \infty} \sup_{l \in L} \left| \int_{\mathbb{R}^{n}} \left[ D^{\alpha} \phi\left(x\right) - D^{\alpha} \phi_{j}\left(x\right) \right] \mathcal{E}_{l}\left(x\right) dx \right|$$

$$\leq M \lim_{j \to \infty} \sup_{x \in K} \left| D^{\alpha} \phi\left(x\right) - D^{\alpha} \phi_{j}\left(x\right) \right| = 0.$$

Thus, as the convergence is uniform, it follows that  $\mathcal{D}[\mathbb{R}^n] \subset KS^2[\mathbb{R}^n]$  as a continuous dense embedding. Since  $KS^2[\mathbb{R}^n]$  is self dual, we see that the Schwartz distributions,  $\mathcal{D}^*[\mathbb{R}^n] \subset KS^2[\mathbb{R}^n]$ .

If  $f, g \in L^1[\mathbb{R}^n]$ , we denote the Fourier transform of f and the convolution of g with respect to f by  $\mathfrak{F}(f)$  and  $\mathfrak{C}_f(g)$  respectively.

**Theorem 3.11.** The Fourier transform  $\mathfrak{F}(\cdot)$  and, for each  $f \in L^1[\mathbb{R}^n]$ , the convolution with respect to f,  $\mathfrak{C}_f(\cdot)$  both extend to  $KS^2[\mathbb{R}^n]$  as bounded linear operators.

- **Remark 3.12.** The possibility that  $KS^2[\mathbb{R}^n]$  could exist was noticed by Gill and Zachary after reading the paper by Kuelbs [KB]. It was actually constructed by V. Steadman and first appeared as a part of her dissertation (see [ST]). As will be seen in the next section, this space is perfect for the Feynman path integral formulation of quantum mechanics.
- 3.2. **Feynman Path Integral I.** The properties of  $KS^2[\mathbb{R}^n]$  derived earlier shows that it is a better replacement for  $L^2[\mathbb{R}^n]$  for the study of the path integral formulation of quantum theory developed by Feynman. It is easy to prove that both the position and momentum operators have closed, densely defined extensions to  $KS^2[\mathbb{R}^n]$ . Furthermore, the extension of  $\mathfrak{F}(\cdot)$  and  $\mathfrak{C}_f(\cdot)$  by Theorem 3.11 insures that both the Schrödinger and Heisenberg theories have faithful representations on  $KS^2[\mathbb{R}^n]$ .

Since  $KS^2[\mathbb{R}^n]$  contains  $\mathfrak{M}[\mathbb{R}^n]$  the space of all finitely and countably additive measures as compact embeddings, it follows that all the approximating sequences for the Dirac measure converge strongly to it in the  $KS^2[\mathbb{R}^n]$  topology. Thus, the finitely additive set function defined on the Borel sets (Feynman kernel):

$$\mathbb{K}_{\mathbf{F}}[t, \mathbf{x}; s, B] = \int_{B} (2\pi i (t - s))^{-1/2} \exp\{i |\mathbf{x} - \mathbf{y}|^{2} / 2(t - s)\} d\mathbf{y}$$

is in  $KS^2[\mathbb{R}^n]$  and  $\|\mathbb{K}_{\mathbf{F}}[t,\mathbf{x}\,;\,s,B]\|_{KS^2}\leqslant 1$  and

$$\mathbb{K}_{\mathbf{F}}[t, \mathbf{x}; s, B] = \int_{\mathbb{R}^n} \mathbb{K}_{\mathbf{F}}[t, \mathbf{x}; \tau, d\mathbf{z}] \mathbb{K}_{\mathbf{F}}[\tau, \mathbf{z}; s, B], \text{ (HK-integral)}.$$

**Definition 3.13.** Let  $\mathbf{P}_n = \{t_0, \tau_1, t_1, \tau_2, \cdots, \tau_n, t_n\}$  be a HK-partition of the interval [0,t] for each n, with  $\lim_{n\to\infty} \Delta \mu_n = 0$  (mesh). Set  $\Delta t_j = t_j - t_{j-1}, \tau_0 = 0$  and for  $\psi \in KS^2[\mathbb{R}^n]$  define

$$\int_{\mathbb{R}^{n[0,t]}} \mathbb{K}_{\mathbf{F}}[\mathcal{D}_{\lambda} \mathbf{x}(\tau) ; \mathbf{x}(0)] = e^{-\lambda t} \sum_{k=0}^{[|\lambda t|]} \frac{(\lambda t)^k}{k!} \left\{ \prod_{j=1}^k \int_{\mathbb{R}^n} \mathbb{K}_{\mathbf{F}}[t_j, \mathbf{x}(\tau_j); t_{j-1}, d\mathbf{x}(\tau_{j-1})] \right\},$$

and

$$\int_{\mathbb{R}^{[0,t]}} \mathbb{K}_{\mathbf{F}}[\mathcal{D}\mathbf{x}(\tau); \mathbf{x}(0)] \psi[\mathbf{x}(0)]$$

$$= \lim_{\lambda \to \infty} \int_{\mathbb{R}^{n[0,t]}} \mathbb{K}_{\mathbf{F}}[\mathcal{D}_{\lambda}\mathbf{x}(\tau); \mathbf{x}(0)] \psi[\mathbf{x}(0)]$$
(3.3)

whenever the limit exists.

The next result is now elementary but a more general (sum over paths) result covering most areas of application will be given after the next section.

Theorem 3.14. The function ψ(x) ≡ 1 ∈ KS<sup>2</sup> [R n ] and

$$\int_{\mathbb{R}^{3[s,t]}} \mathbb{K}_{\mathbf{F}}[\mathcal{D}\mathbf{x}(\tau) ; \mathbf{x}(s)] = \mathbb{K}_{\mathbf{F}}[t, \mathbf{x}; s, \mathbf{y}] = \frac{1}{\sqrt{2\pi i(t-s)}} \exp\{i|\mathbf{x} - \mathbf{y}|^2 / 2(t-s)\}.$$

The above result is exactly what Feynman expected to obtain.

Remark 3.15. The above results also hold with no changes if we replace KF[t, x ; s, B] by the kernel KW[t, x ; s, B], for the Wiener measure.

Since L 2 [R n ] is dense, KS<sup>2</sup> [R n ] is separable and every orthonormal basis for L 2 [R n ] may be normalized to make an orthonormal basis for KS<sup>2</sup> [R n ]. Thus they are unitarily equivalent but:

- (1) KS<sup>2</sup> [R n ] contains both D[R n ] and D<sup>∗</sup> [R n ], while L 2 [R n ] does not.
- (2) The integral of e ix<sup>2</sup> is finite in KS<sup>2</sup> [R n ] but infinite in L 2 [R n ].

Thus, L 2 [R n ] and KS<sup>2</sup> [R n ] are mathematically equivalent but not physically equivalent.

We will show that interactions of a general nature become possible when we consider the path integral as a part of Feynman's time ordered operator calculus in the next section.

## 4. Feynman Operator Calculus and the Dyson Conjectures

The purpose of this section is to provide an introduction to the Feynman operator calculus. This section closely follows the work in [24]. Those interested in more detail along with proofs can consult [27]. We first need a universal order parameter (see [24, 25]).

4.1. Horwitz-Fanchi Time. Newton suggested that the universe had an absolute definition of time which existed independent of any observer. He assume that this clock and provides a unique definition of simultaneity and ordering for events, but may be only understood as a mathematical object. (In Newton's era mathematics and personal philosophy still played a dominant role in one's physics). Einstein's special theory became problematic for Newton's absolute time, while Minkowski's postulate eliminated any sense of ordering. Feynman later suggested that we view a physical event as a 3 dimensional motion picture exposing more of the outcome as the film unfolds (reintroduced the ordering property).

The first to use a fifth parameter to account for the missing order property Fock, and later Stueckelberg followed by Feynman and Schwinger. They each considered it a global clock as an addition to Minkowski space time. Horwitz and Piron [34], and Fanchi and Collins gave it the name historical time. They were the first to realize that there should be some physical justification for this clock (see [16]).

The purpose of this section is to provide a proof of their conjecture. Our proof uses the third definition of proper time and is based on a relatively weak assumption about the properties of the universe.

**Theorem 4.1.** If the universe is representable in the sense that it is independent of any observer's particular portion of the universe, then the universe has a unique clock.

*Proof.* If  $Mc^2$  is the total mass energy and H is the total energy Hamiltonian for the universe as seen by O, under the stated conditions  $H/Mc^2$  is constant. By assumption this view is observer independent, so that O' will also obtain the same ratio. It follows from (1.3) that

$$d\tau = \frac{Mc^2}{H}dt = \frac{Mc^2}{H}dt', \quad \Rightarrow t = t' \equiv \tau_h$$

and  $\tau_h$  uniquely defines a clock for the universe.

The principle of homogeneity states that the distribution of matter does not vary significantly with position on large scale. Several recent studies suggest that the principle of homogeneity is not as secure as believed in the past (see [42]). Studies have also recently led to questions about the principle of isotropy (see [5, 57]). Thus, if future research support these results, the cosmological principle may not hold. However, the assumption of Theorem 4.1 can remain valid.

- 4.2. **Feynman-Dyson Space.** Let  $\mathcal{H} = KS^2[\mathbb{R}^n]$ , with a fixed orthonormal basis  $\{e^i\}$ . Let J = [-T, T] be an interval of HF-time and, for each  $t \in J$ , define  $\mathcal{H}(t) = \mathcal{H}$  and let  $\mathcal{H}_{\otimes} = \hat{\otimes}_t \mathcal{H}(t)$  be the continuous tensor product Hilbert space of von Neumann over J. For each  $i \in \mathbb{N}$ , set  $e^i_t = e^i$  and  $E^i = \otimes_t e^i_t$ . We define  $\mathcal{FD}^i \subset \mathcal{H}_{\otimes}$  to be the smallest Hilbert space containing  $E^i$  and call  $\mathcal{FD} = \bigoplus_{i=1}^{\infty} \mathcal{FD}^i$  the Feynman-Dyson space over J for  $\mathcal{H}$ . (This is the film for Feynman's space-time events.)
- 4.3. **Time Ordered Operators.** If  $\mathcal{C}(\mathcal{H}_{\otimes})$  is a set of closed densely defined linear operators on  $\mathcal{H}_{\otimes}$ , and  $\{H_I(t), t \in J\}$  is a family of generators for unitary groups, we define  $\mathcal{C}(\mathcal{H}(t)) \subset \mathcal{C}(\mathcal{H}_{\otimes})$  by:

$$\mathcal{C}(\mathcal{H}(t)) = \left\{ \mathbf{H}_I(t) \mid \mathbf{H}_I(t) \right\} = \widehat{\bigotimes}_{b \geqslant s > t} \mathbf{I}_s \otimes H_I(t) \otimes (\bigotimes_{t > s \geqslant -T} I_s) ,$$

where  $I_s$  is the identity operator at time s. It follows that

$$\mathbf{H}_I(t)\mathbf{H}_I(s) = \mathbf{H}_I(s)\mathbf{H}_I(t), t \neq s.$$

Thus, the operators are ordered in time, commute when acting at different times and still maintain their mathematically defined position on paper.

4.4. **Time Ordered Integrals.** If  $\{H_I(t)|t \in J\}$  is a family of time dependent Hamiltonian generators of unitary groups on  $\mathcal{H}$  and  $\{\mathbf{H}_I(t)|t \in J\}$  is the time ordered version defined on  $\mathcal{FD}$  then:

**Theorem 4.2.** The time ordered integral  $\mathbf{Q}[t, -T] = \int_{-T}^{t} \mathbf{H}_{I}(s) ds$  exists (a.e), has a dense domain and is the generator of a strongly continuous unitary group  $\mathbf{U}[t, -T]$ . If  $\Psi_{0}$  is in the domain of  $\mathbf{Q}[t, -T]$  then

(1)
$$\Psi(t) = \mathbf{U}[t, -T] \Psi_0 = \exp\left\{-\frac{i}{\hbar}\mathbf{Q}[t, -T]\right\} \Psi_0, \quad satisfies:$$
(2)
$$i\hbar \frac{\partial \Psi(t)}{\partial t} = \mathbf{H}_I(t) \Psi(t), \quad \Psi(-T) = \Psi_0.$$

4.5. Interaction Representation. A theorem of Haag shows that the equal time commutation relations of an interacting field are equivalent to those of a free field (see Streater and Wightman [STW], pg. 161). A recent experiment of Lindner et. al. shows there is some quantum interference in time for the wave function of a particle (see Horwitz [HW]). In this section we show the interaction representation is well-defined if time smearing exists.

Let us assume that  $\boldsymbol{H}_0(t)$  and  $\boldsymbol{H}_1(t)$  are generators of unitary groups for each  $t \in J$ ,  $\boldsymbol{H}_I(t) = \boldsymbol{H}_0(t) \oplus \boldsymbol{H}_1(t)$  is densely defined,  $\boldsymbol{H}_1^n(t) = n\boldsymbol{H}_1(t)R(n,\boldsymbol{H}_1(t))$  is the Yosida approximator for  $\boldsymbol{H}_1(t)$  and Theorem (4.2) is satisfied. Define  $\mathbf{U}_n[t,a]$ ,  $\mathbf{U}_0[t,a]$  and  $\mathbf{U}_0^n[t,a]$  by:

$$\begin{aligned} \mathbf{U}_n[t,a] &= \exp\{(-i/\hbar)\int\limits_a^t [\bm{H}_0(s)\oplus \bm{H}_1^n(s)]ds\}, \ \mathbf{U}_0[t,a] &= \exp\{(-i/\hbar)\int\limits_a^t \bm{H}_0(s)ds\}, \ \bm{U}_0^\kappa[t,a] &= \exp\{(-i/\hbar)\int\limits_a^t \bm{H}_0^\kappa(s)ds\}, \ \bm{H}_0^\kappa(t) &= \int\limits_{-\infty}^\infty \rho_\kappa(t,s)\mathbf{E}[t,s]\bm{H}_0(s)ds, \end{aligned}$$

and  $\rho_{\kappa}(t,s)$  is our smearing density, which we take to depend on the Planck length  $\kappa$  and  $\int_{-\infty}^{\infty} \rho_{\kappa}(t,s) ds = 1$  (for example,  $\rho_{\kappa}(t,s) = [1/\sqrt{2i\pi\kappa^2}] \exp\{i(t-s)^2/2\kappa^2\}$ ).

We now obtain:

$$\boldsymbol{H}_{I}^{n}(t) = \mathbf{U}_{0}^{\kappa}[a, t]\boldsymbol{H}_{1}^{n}(t)\mathbf{U}_{0}^{\kappa}[t, a],$$

and the terms do not commute. If we set  $\Psi_n(t) = \mathbf{U}_0^{\kappa}[a,t]\mathbf{U}_n[t,a]\Phi$ , we have

$$\begin{split} &\frac{\partial}{\partial t}\Psi_n(t)=\frac{i}{\hbar}\mathbf{U}_0^{\kappa}[a,t]\boldsymbol{H}_0(t)\mathbf{U}_n[t,a]\Phi-\frac{i}{\hbar}\mathbf{U}_0^{\kappa}[a,t]\left[\boldsymbol{H}_0(t)+\boldsymbol{H}_1^n(t)\right]\mathbf{U}_n[t,a]\Phi\\ &\text{so that}\quad \frac{\partial}{\partial t}\Psi_n(t)=\frac{i}{\hbar}\{\mathbf{U}_0^{\kappa}[a,t]\boldsymbol{H}_1^n(t)\mathbf{U}_0^{\kappa}[t,a]\}\mathbf{U}_0^{\kappa}[a,t]\mathbf{U}_n[t,a]\Phi\\ &\text{and}\quad i\hbar\frac{\partial}{\partial t}\Psi_n(t)=\boldsymbol{H}_I^n(t)\Psi_n(t),\ \Psi_n(a)=\Phi. \end{split}$$

With the same conditions as Theorem (4.2), we have

**Theorem 4.3.** If  $Q_1[t,a] = \int_a^t H_1(s)ds$  generates a unitary group on  $\mathcal{H}$ , then the time-ordered integral  $\mathbf{Q}_I[t,a] = \int_a^t \mathbf{H}_I(s)ds$ , where  $\mathbf{H}_I(t) = \mathbf{U}_0^{\kappa}[a,t]\mathbf{H}_1(t)\mathbf{U}_0^{\kappa}[t,a]$  generates a unitary group on  $\mathcal{FD}_{\otimes}^2$ , and

$$\exp\{(-i/\hbar)\mathbf{Q}_I^n[t,a]\} \to \exp\{(-i/\hbar)\mathbf{Q}_I[t,a]\},$$

where  $\mathbf{Q}_I^n[t,a] = \int_a^t \mathbf{H}_I^n(s) ds$ , and:

$$i\hbar \frac{\partial}{\partial t} \Psi(t) = \boldsymbol{H}_I(t)\Psi(t), \ \Psi(a) = \Phi.$$

- 4.6. **The Dyson Conjectures.** Freeman Dyson was the first to understand the different approaches of both Feynman and Schwinger-Tomonaga. He showed that the real difference represented two ways of looking at Heisenberg's S-matrix. This insight allowed him to both simplify and unify the Feynman-Schwinger-Tomonaga theories. Following this work, Dyson made his four famous conjectures about QED.
  - (1) He first conjectured that the divergences in QED were due to an idealized conception of measurability resulting from the assumed infinitely precise knowledge of the space-time positions of particles implied by the Hamiltonian formulation, which leads to a violation of the Heisenberg uncertainty principle.
  - (2) His secondly conjectured that, in general we could only expect the renormalized perturbation expansion of QED to be asymptotic.
  - (3) His third conjecture was that a particular set of Feynman graphs, which gave overlapping divergencies may be made to systematically cancel out (proven by Salam [55]).
  - (4) His forth conjecture was that, a proof that the ultraviolet divergences cancel to all orders demanded that that a certain Feynman integral must converge (proven by Weinberg [63]).

His first two conjectures remained open due to a lack of mathematical meaning for Feynman's time operator calculus that Dyson (following Feynman) used in his analysis to provide a physical interpretation. The first conjecture also required a formulation of the S-matrix that provided a physically meaningful reason for the ultraviolet divergency as a violation of the uncertainty principle. In this section we provide proofs of these two remaining conjectures (see [27]).

4.6.1. The Second Conjecture. The second conjecture is essentially a mathematical problem that could not be considered without meaning for Feynman's operator calculus, so we prove it first.

**Definition 4.4.** The operator  $\mathbf{U}^w[t,a] = \exp\{w\mathbf{Q}[t,a]\}$  is asymptotic in the sense of Poincaré if, for each n and each  $\Psi_a \in D\left[(\mathbf{Q}[t,a])^{n+1}\right]$  (domain),

$$\lim_{w \to 0} w^{-(n+1)} \left\{ \mathbf{U}^w[t, a] - \sum_{k=1}^n \frac{(w\mathbf{Q}[t, a])^k}{k!} \right\} \Psi_a = \frac{\mathbf{Q}[t, a]^{n+1}}{(n+1)!} \Psi_a.$$
 (4.1)

This is the operator version of an asymptotic expansion corresponds to the classical sense of asymptotic expansion for an infinite series of functions as given by Poincaré (see [35]).

If w is a parameter (charge), and  $\mathbf{Q}[t,a]$  generates a unitary group on  $\mathcal{FD}$ , we have:

**Theorem 4.5.** For -T < a < t < T,

- (1)  $\mathbf{U}^w[t,a] = exp\{w\mathbf{Q}[t,a]\}$  is asymptotic in the sense of Poincaré.
- (2) If  $\Psi_a \in D\left[ (\mathbf{Q}[t,a])^{n+1} \right]$  and  $\Psi(t) = \mathbf{U}^w[t,a]\Psi_a$ , for each n we have:

$$\Psi(t) = \Psi_a + \sum_{k=1}^n w^k \int_a^t ds_1 \int_a^{s_1} ds_2 \cdots \int_a^{s_{k-1}} ds_k \mathbf{H}_I(s_1) \cdots \mathbf{H}_I(s_k) \Psi_a$$

$$+ \int_0^w (w - \xi)^n d\xi \int_a^t ds_1 \cdots \int_a^{s_n} ds_{n+1} \mathbf{H}_I(s_1) \cdots \mathbf{H}_I(s_{n+1}) \mathbf{U}^{\xi}[s_{n+1}, a] \Psi_a.$$

$$(4.2)$$

The above theorem provides a precise statement of Dyson's first conjecture ([27]). Equation (4.2) is the corrected Dyson expansion, in that it provides the remainder term, making the expansion exact for all n (see ([27]).

Remark 4.6. There are rare cases in which the perturbation series converge to the solution. There are also cases where the renormalized series may diverge, but still respond to some summability method (usually Borel). A good discussion, with references, can be found in Glimm and Jaffe [21].

4.6.2. Dyson's first Conjecture. In order to provide a basis for Dyson's first conjecture, we define an exchange operator.

**Definition 4.7.** An exchange operator E[t,t'] on  $C[\mathcal{H}_{\otimes}]$  is a linear map defined for pairs t, t' such that:

- (1)  $E[t,t']: \mathcal{C}[\mathcal{H}(t)] \to \mathcal{C}[\mathcal{H}(t')], (bijective mapping),$
- (2) E[s, t']E[t, s] = E[t, t'],(3)  $E[t, t']E[t', t] = \mathbf{I}_{\otimes}, (identity)$

(4) for 
$$s \neq t$$
,  $t'$ ,  $E[t, t']\mathbf{H}_I(s) = \mathbf{H}_I(s)$ , for all  $\mathbf{H}_I(s) \in \mathcal{C}[\mathcal{H}(s)]$ .

This operator will allow the exchange of time positions between a pair of operators in an expression. We define the experimental time ordered integral with information concentrated at time points  $\{\tau_1, \tau_2, \cdots, \tau_n\}$  by:

$$\mathbf{Q}_{e}\left[\tau_{1}, \tau_{2}, \dots \tau_{n}\right] = \sum_{j=1}^{n} \int_{t_{j-1}}^{t_{j}} E\left[\tau_{j}, s\right] \mathbf{H}_{I}(s) ds$$
$$= \sum_{j=1}^{n} \Delta t_{j} \left[\frac{1}{\Delta t_{j}} \int_{t_{j-1}}^{t_{j}} E\left[\tau_{j}, s\right] \mathbf{H}_{I}(s) ds\right].$$

For each interval  $[t_{j-1}, t_j]$ , the average is concentrated at the point  $\tau_j$ . Let  $\Psi \in \mathcal{FD}$  and let N(t) be a Poisson distributed random variable with parameter  $\lambda$ , we define the function  $\mathbf{U}[N(t), 0]\Psi$  by:

$$\mathbf{U}[N(t), 0]\Psi = \exp\left\{\mathbf{Q}_e\left[\tau_1, \tau_2 \cdots, \tau_{N(t)}\right]\right\}\Psi. \tag{4.3}$$

 $\mathbf{U}[N(t),0]\Psi$  is a  $\mathcal{FD}$ -valued random variable which distributes of the number of and time placements of information points on a film between times  $a \to t$ . To relate  $\mathbf{U}[N(t),0]\Psi$  to experiment, we need to compute its expected value.

$$\bar{\mathbf{U}}_{\lambda}[t,0]\Psi = \mathcal{E}\left\{\mathbf{U}[N(t),0]\Psi\right\}$$

$$= \sum_{n=0}^{\infty} \mathcal{E}\left\{\mathbf{U}[N(t),0]\Psi \mid N(t)=n\right\} Prob[N(t)=n],$$
(4.4)

$$\mathcal{E}\left\{\mathbf{U}[N(t),0]\Psi\mid N(t)=n\right\}$$

$$=\int_{0}^{t}\frac{d\tau_{1}}{t-0}\int_{\tau_{1}}^{t}\frac{d\tau_{2}}{t-\tau_{1}}\cdots\int_{\tau_{n-1}}^{t}\frac{d\tau_{n}}{t-\tau_{n-1}}\mathbf{U}[\tau_{n},\tau_{n-1}\cdots,\tau_{1}]\Psi=\bar{\mathbf{U}}_{n}[t,0]\Psi, \tag{4.5}$$

with  $Prob[N(t) = n] = \frac{(\lambda t)^n}{n!} \exp\{-\lambda t\}$ . Since we are only interested in asymptotic results, when  $\lambda \to \infty$ , we can take  $\tau_j = \frac{jt}{n}$ ,  $1 \le j \le n$ ,  $\Delta t_j = \frac{t}{n}$  and replace  $\bar{\mathbf{U}}_n[t,0]\Psi$  by  $\mathbf{U}_n[t,0]\Psi$ , so that

$$\mathbf{U}_n[t,0]\Psi = \exp\left\{\sum_{j=1}^n \int_{t_{j-1}}^{t_j} E[\tau_j,s]\mathbf{H}_I(s)ds\right\}\Psi.$$

We define the experimental evolution operator  $\mathbf{U}_{\lambda}[t,0]\Psi$  by

$$\mathbf{U}_{\lambda}[t,0]\Psi = \sum_{n=0}^{\infty} \frac{(\lambda t)^n}{n!} \exp\{-\lambda t\} \mathbf{U}_n[t,0]\Psi. \tag{4.6}$$

Borel summability is regular,  $\mathbf{U}_n[t,0]\Psi \to \mathbf{U}[t,0]\Psi$  implies  $\mathbf{U}_{\lambda}[t,0]\Psi \to \mathbf{U}[t,0]\Psi$ . Thus, by Theorem 4.2:

$$\lim_{\lambda \to \infty} \mathbf{U}_{\lambda}[t, 0] \Psi = \mathbf{U}[t, 0] \Psi. \tag{4.7}$$

As  $\lambda \to \infty \Rightarrow \lambda^{-1} \to 0$ , we obtain a continuous (mean) path. This continuous path is the result of averaging over an infinite number of paths (discrete).

We now replace  $\mathbf{H}_I(t)$  by  $\frac{-i}{\hbar}\mathbf{H}_I(t)$ , and  $\mathbf{U}_{\lambda}[T, -T]\Psi$  by the experimental S-matrix  $\mathbf{S}_{\lambda}[T, -T]\Psi$ , such that

$$\mathbf{S}_{\lambda}[T, -T]\Psi = \sum_{n=0}^{\infty} \frac{(2\lambda T)^n}{n!} \exp\left[-2\lambda T\right] \mathbf{S}_n[T, -T]\Psi. \tag{4.8}$$

Combining  $\exp[-2\lambda T]$  and  $\mathbf{S}_n[T, -T]\Psi$ , we obtain:

$$\exp\left[-2\lambda T\right]\mathbf{S}_n[T, -T]\Psi = e^{-\frac{i}{\hbar}\left\{\sum_{j=1}^n \int_{t_{j-1}}^{t_j} \left[E[\tau_j, s]\mathbf{H}_I(s) - i\lambda\hbar\mathbf{I}_{\otimes}\right]ds\right\}}\Psi. \tag{4.9}$$

To obtain a physical interpretation of this equation, we take a closer look at the term in the exponential in equation (4.9):

$$\sum_{j=1}^{n} \left\{ \frac{1}{\Delta t_j} \int_{t_{j-1}}^{t_j} \left[ E\left[\tau_j, s\right] \mathbf{H}_I(s) - i\lambda \hbar I_{\otimes} \right] ds \right\} \Delta t_j. \tag{4.10}$$

The term  $-i\lambda\hbar I_{\otimes}$  corresponds to absorption of photon energy of amount  $\lambda\hbar$  in each interval  $[t_{j-1},t_j]$  and concentrated at the point  $\tau_j$ . Taking the limit in equation (4.9), we obtain the S-matrix. This now means that an infinite amount of photon energy is deposited at each point of the interval [-T,T]. This is clearly the ultra-violet divergence. To relate it to the time-energy uncertainty relation, we recall that it is a variance, so it can only appear in the second order of the perturbation expansion and all higher order terms. This is what we see in QED. In addition, the method of elimination for the second order terms was shown to work for the higher order ones by Weinberg (as predicted by Dyson).

- Remark 4.8. There are two other divergencies left to accounted for in QED. If we let  $T \to \infty$ , we obtain the infrared divergence, while the self-energy divergence is resolved by introducing an infinite mass counter term in the Hamiltonian density  $\mathbf{H}_I$ . These will be addressed in the second paper of this series (but see [24]).
- 4.7. Path Integrals II. The purpose of this section is show that the path integral is a special case of the time-ordered operator calculus (as conjectured by Feynman).

Let U[t, a] be an evolution operator on  $KS^2(\mathbb{R}^3)$ , with time-dependent generator H(t), and reproducing kernel  $K[\mathbf{x}(t), t; \mathbf{x}(s), s]$ :

$$\begin{split} K\left[\mathbf{x}(t),\,t;\,\mathbf{x}(s),\,s\right] &= \int_{\mathbb{R}^3} K\left[\mathbf{x}(t),\,t;\,d\mathbf{x}(\tau),\,\tau\right] K\left[\mathbf{x}(\tau),\,\tau;\,\mathbf{x}(s),\,s\right],\\ U[t,s]\varphi(s) &= \int_{\mathbb{R}^3} K\left[\mathbf{x}(t),\,t;\,d\mathbf{x}(s),\,s\right]\varphi(s). \end{split}$$

If  $\mathbf{U}[t,s]$  is the corresponding time-ordered version defined on  $\mathcal{FD}^2_{\otimes} \subset \mathcal{H}^2_{\otimes}$ , with kernel  $\mathbb{K}_{\mathbf{f}}[\mathbf{x}(t), t; \mathbf{x}(s), s]$ . Since  $\mathbf{U}[t,\tau]\mathbf{U}[\tau,s] = \mathbf{U}[t,s]$ , we have:

$$\mathbb{K}_{\mathbf{f}}\left[\mathbf{x}(t), t; \mathbf{x}(s), s\right] = \int_{\mathbb{R}^3} \mathbb{K}_{\mathbf{f}}\left[\mathbf{x}(t), t; d\mathbf{x}(\tau), \tau\right] \mathbb{K}_{\mathbf{f}}\left[\mathbf{x}(\tau), \tau; \mathbf{x}(s), s\right].$$

From our sum over paths representation for  $\mathbf{U}[t,s]$ , we have:

$$\mathbf{U}[t, s]\Phi(s) = \lim_{\lambda \to \infty} \mathbf{U}_{\lambda}[t, s]\Phi(s)$$
$$= \lim_{\lambda \to \infty} e^{-\lambda(t-s)} \sum_{k=0}^{n} \frac{[\lambda(t-s)]^{k}}{k!} \mathbf{U}_{k}[t, s]\Phi(s),$$

where  $n = [|\lambda(t-s)|]$  (the greatest integer in  $\lambda(t-s)$ ) and

$$\mathbf{U}_k[t,s]\Phi(s) = \exp\left\{(-i/\hbar)\sum_{j=1}^k \int_{t_{j-1}}^{t_j} \mathbf{E}[\tau_j,\tau]\mathbf{H}(\tau)d\tau\right\}\Phi(s).$$

As in Section 3.2, define  $\mathbb{K}_{\mathbf{f}}[\mathcal{D}_{\lambda}\mathbf{x}(t) ; \mathbf{x}(s)]$  by:

$$\mathbb{K}_{\mathbf{f}}[\mathcal{D}_{\lambda}\mathbf{x}(t) ; \mathbf{x}(s)]$$

$$=:e^{-\lambda(t-s)}\sum_{k=1}^{n}\frac{\left[\lambda(t-s)\right]}{k!}\left\{\prod_{j=1}^{k}\int_{\mathbb{R}^{3}}\mathbb{K}_{f}\left[t_{j},\mathbf{x}\left(t_{j}\right);d\mathbf{x}\left(t_{j-1}\right),t_{j-1}\right]|^{\tau_{j}}\right\}$$

and  $|\tau_j|$  denotes that the integration is performed at the time  $\tau_j$ .

**Definition 4.9.** We define the Feynman path integral associated with U[t, s] by:

$$\mathbf{U}[t,s]\Phi(s) = \int_{\mathbb{R}^{3[t,s]}} \mathbb{K}_{\mathbf{f}}[\mathcal{D}\mathbf{x}(\tau) \; ; \; \mathbf{x}(s)]\Phi(s) = \lim_{\lambda \to \infty} \int_{\mathbb{R}^{3[t,s]}} \mathbb{K}_{\mathbf{f}}[\mathcal{D}_{\lambda}\mathbf{x}(\tau) \; ; \; \mathbf{x}(s)]\Phi(s).$$

**Theorem 4.10.** For the Feynman time-ordered theory, whenever a reproducing kernel exists on  $KS^2[\mathbb{R}^3]$ , we have:

$$\lim_{\lambda \to \infty} \mathbf{U}_{\lambda}[t, s] \Phi(s) = \mathbf{U}[t, s] \Phi(s) = \int_{\mathbb{R}^{3[t, s]}} \mathbb{K}_{\mathbf{f}}[\mathcal{D}^{\lambda} \mathbf{x}(\tau) ; \mathbf{x}(s)] \Phi[\mathbf{x}(s)],$$

and the limit is independent of the space of continuous functions.

Let us assume that  $H_0(t)$  and  $H_1(t)$  are strongly continuous generators of unitary groups, with a common dense domain D(t), for each  $t \in J = [a, b]$  and, let  $\mathbf{H}_{1,\rho}(t) = \rho \mathbf{H}_1(t) R(\rho, \mathbf{H}_1(t))$  be the Yosida approximator for the time-ordered version of  $H_1(t)$ , with dense domain  $D = \bigotimes_{t \in I} D(t)$ . Define  $\mathbf{U}^{\rho}[t, a]$  and  $\mathbf{U}^0[t, a]$  by:

$$\begin{split} \mathbf{U}^{\rho}[t,a] &= \exp\{(-i/\hbar) \int\limits_{a}^{t} [\boldsymbol{H}_{0}(s) + \boldsymbol{H}_{1,\rho}(s)] ds\}, \\ \mathbf{U}^{0}[t,a] &= \exp\{(-i/\hbar) \int\limits_{a}^{t} \boldsymbol{H}_{0}(s) ds\}. \end{split}$$

Since  $H_{1,\rho}(s)$  is bounded,  $H_0(s) + H_{1,\rho}(s)$  is a generator of a unitary group for each  $s \in J$  and finite  $\rho$ . Now assume that  $\mathbf{U}^0[t,a]$  has an associated reproducing kernel, so that  $\mathbf{U}^0[t,a] = \int_{\mathbb{R}^3[t,s]} \mathbb{K}_{\mathbf{f}}[\mathcal{D}\mathbf{x}(\tau);\mathbf{x}(a)]$ . We now have

the following general result, which is independent of the space of continuous functions.

**Theorem 4.11.** (Extended Feynman-Kac) If  $\mathbf{H}_0(s) \oplus \mathbf{H}_1(s)$  is a generator of a unitary group, then

$$\begin{split} &\lim_{\rho \to \infty} \mathbf{U}^{\rho}[t,a] \Phi(a) = \mathbf{U}[t,a] \Phi(a) \\ &= \int_{\mathbb{R}^{3[t,a]}} \mathbb{K}_{\mathbf{f}}[\mathcal{D}\mathbf{x}(\tau) \; ; \; \mathbf{x}(a)] \exp\{(-i/\hbar) \int_{a}^{\tau} \boldsymbol{H}_{1}(s) ds]\} \Phi[\mathbf{x}(a)]. \end{split}$$

*Proof.* The fact that  $\mathbf{U}^{\rho}[t,a]\Phi(a) \to \mathbf{U}[t,a]\Phi(a)$  is clear. To prove that

$$\mathbf{U}[t,a]\Phi(a) = \int_{\mathbb{R}^{3[t,a]}} \mathbb{K}_{\mathbf{f}}[\mathcal{D}\mathbf{x}(\tau);\mathbf{x}(a)] \exp\{(-i/\hbar) \int_{a}^{t} \mathbf{H}_{1}(s) ds\} \Phi(a),$$

first note that, since the time-ordered integral exists and we are only interested in the limit, we can write for each  $\boldsymbol{k}$ 

$$U_k^{\rho}[t,a]\Phi(a) = \exp\left\{\left(-i/\hbar\right)\sum_{j=1}^k \int_{t_{j-1}}^{t_j} \left[\mathbf{E}[\tau_j,s]\boldsymbol{H}_0(s) + \mathbf{E}[\tau_j',s]\boldsymbol{H}_{1,\rho}(s)\right] ds\right\}\Phi(a),$$

where  $\tau_j$  and  $\tau'_j$  are distinct points in the interval  $(t_{j-1}, t_j)$ . Thus, we can also write  $U_k^{\rho}[t, a]\Phi(a)$  as

$$\begin{aligned} &\mathbf{U}_{k}^{\rho}[t,a] = exp\left\{\frac{-i}{\hbar}\sum_{j=1}^{k}\int_{t_{j-1}}^{t_{j}}\mathbf{E}[\tau_{j},s]\boldsymbol{H}_{0}(s)ds\right\}exp\left\{\frac{-i}{\hbar}\sum_{j=1}^{k}\int_{t_{j-1}}^{t_{j}}\mathbf{E}[\tau'_{j},s]\boldsymbol{H}_{1,\rho}(s)ds\right\} \\ &= \prod_{j=1}^{k}exp\left\{\frac{-i}{\hbar}\int_{t_{j-1}}^{t_{j}}\mathbf{E}[\tau_{j},s]\boldsymbol{H}_{0}(s)ds\right\}exp\left\{\frac{-i}{\hbar}\sum_{j=1}^{k}\int_{t_{j-1}}^{t_{j}}\mathbf{E}[\tau'_{j},s]\boldsymbol{H}_{1,\rho}(s)ds\right\} \\ &= \prod_{j=1}^{k}\int_{\mathbb{R}^{3}}\mathbb{K}_{\mathbf{f}}\left[t_{j},\mathbf{x}(t_{j})\;;\;t_{j-1},d\mathbf{x}(t_{j-1})\right]exp\left\{\frac{-i}{\hbar}\sum_{j=1}^{k}\int_{t_{j-1}}^{t_{j}}\mathbf{E}[\tau'_{j},s]\boldsymbol{H}_{1,\rho}(s)ds\right\}. \end{aligned}$$

If we put this in our experimental evolution operator  $\mathbf{U}_{\lambda}^{\rho}[t,a]\Phi(a)$  and compute the limit, we have:

$$\mathbf{U}^{\rho}[t, a]\Phi(a)$$

$$= \int_{\mathbb{D}^{3[t, a]}} \mathbb{K}_{\mathbf{f}}[\mathcal{D}\mathbf{x}(t); \mathbf{x}(a)] \exp\left\{ (-i/\hbar) \int_{a}^{t} \boldsymbol{H}_{1, \rho}(s) ds \right\} \Phi(a).$$

Since the limit as  $\rho \to \infty$  on the left exists, it defines the limit on the right.

4.8. **Examples.** In this section, we pause to discuss a few examples. Theorem 4.11 is somewhat abstract, so our first example covers almost all of non-relativistic quantum theory.

Let  $\Delta$  be the Laplacian on  $\mathbb{R}^n$  and let V be any potential such that  $H = (-\hbar^2/2)\Delta + V$  generates a unitary group. Then, using time as an index, the problem:

$$(i\hbar)\partial\psi(\mathbf{x},t)/\partial t = \mathbf{H}(t)\psi(\mathbf{x},t), \ \psi(\mathbf{x},0) = \psi_0(\mathbf{x}),$$

has a solution with the extended Feynman-Kac representation.

Our second example is due to Albeverio and Mazzucchi [1]. Let  $\mathbb{C}$  be a completely symmetric positive definite fourth-order covariant tensor on  $\mathbb{R}^n$ , let  $\Omega$  be a symmetric positive definite  $n \times n$  matrix and let  $\lambda$  be a nonnegative constant. Then:

$$H = -\frac{\hbar^2}{2} \mathbf{\Delta} + \frac{1}{2} \mathbf{x} \Omega^2 \mathbf{x} + \lambda \mathbb{C} [\mathbf{x}, \mathbf{x}, \mathbf{x}, \mathbf{x}]$$

is known to be the generator for a unitary group on  $L^2[\mathbb{R}^n]$ . Albeverio and Mazzucchi [1] prove that  $\bar{A}$  has a path integral representation as the analytic continuation (in the parameter  $\lambda$ ) of an infinite dimensional generalized oscillatory integral. (Their version of Feynman's path integral.)

Using the results of the previous sections, we can extension H to  $KS^2[\mathbb{R}^n]$ , which still generates a unitary group. Let  $V = \frac{1}{2}\mathbf{x}\Omega^2\mathbf{x} + \lambda\mathbb{C}[\mathbf{x},\mathbf{x},\mathbf{x},\mathbf{x}]$  and  $V_{\rho} = V(I + \rho V^2)^{-1/2}$ ,  $\rho > 0$ . We can prove that  $V_{\rho}$  is a bounded generator which converges to V. Since  $-\frac{\hbar^2}{2}\mathbf{\Delta}$  generates a unitary group,  $H_{\rho} = -\frac{\hbar^2}{2}\mathbf{\Delta} + V_{\rho}$  also generates one and converges to H. Let

$$\boldsymbol{H}(\tau) = (\mathop{\hat{\otimes}}_{t \geq s > \tau} \mathbf{I}_s) \otimes \boldsymbol{H} \otimes (\mathop{\otimes}_{\tau > s \geq 0} \mathbf{I}_s),$$

then  $\boldsymbol{H}(t)$  generates a unitary group for each t and  $\boldsymbol{H}_{\rho}(t)$  converges to  $\boldsymbol{H}(t)$  on  $\mathcal{FD}_{\otimes}^{2}$ . We can now obtain:

$$\mathbf{U}[t, a]\Phi = \int_{\mathbb{R}^{3[t, a]}} \mathbb{K}_{\mathbf{F}}[\mathcal{D}\mathbf{x}(\tau); \mathbf{x}(a)] \exp\left\{-(i/\hbar) \int_{a}^{\tau} V(s) ds\right\} \Phi$$
$$= \lim_{\rho \to 0} \int_{\mathbb{R}^{3[t, a]}} \mathbb{K}_{\mathbf{F}}[\mathcal{D}\mathbf{x}(\tau); \mathbf{x}(a)] \exp\left\{-(i/\hbar) \int_{a}^{\tau} V_{\rho}(s) ds\right\} \Phi.$$

We refer to [28] for additional examples, including path integrals for kernels that are not perturbations of the Laplacian.

#### **Discussion and Conclusion**

In this paper we have faced directly a number of mathematical problems that have clouded our ability to see how to correctly address some of the physical problems left over from the last century.

- 4.9. Classical Level. At the classical level, we have made Minkowski's unacknowledged third (mathematical) postulate for Einstein's special theory of relativity explicit and provided a direct proof that it is incompatible with the first postulate for two or more particles. We then use one of three representations for the proper time to show how to resolve this incompatibility. This then allowed us to resolve the center of mass problem first noticed by Pryce [51] and the no Interaction problem first conjectured by Bakamjian and Thomas [7] and later proved by Currie, Jordan and Sudarshan [9]. This clears the way for a direct extension of Newton's theory to one compatible with that of Maxwell (addressed in the second paper).
- 4.10. **Quantum Level.** At the quantum level, we have first provided an analytical separation (diagonalization) of the full (minimal coupling) Dirac equation into particle and antiparticle components. The diagonalization is analytic in that it was achieved without transforming the wave functions and revealed the nonlocal time behavior of the particle-antiparticle relationship. We also derived the probability densities, providing a complete decoupling of the spinor space  $L^2[\mathbb{R}^3, \mathbb{C}^4] = L^2[\mathbb{R}^3, \mathbb{C}^2] \oplus L^2[\mathbb{R}^3, \mathbb{C}^2]$ . We have given the correct physical interpretation of the zitterbewegung and the fact that the expected instantaneous value of a velocity measurement of a Dirac particle  $= \pm c$  is that the Dirac equation makes a spatially extended particle appear as a point in the present by forcing it to oscillate instantaneously between the past and future with speed c (at each point in time).
- 4.10.1. Pauli Approximation. It is generally agreed that QED is in excellent agreement with experiments. The fact that it is very successful is not in doubt. However, there are still some serious technical and foundational issues which must be addressed. Historically, when Lamb and Retherford confirmed suspicions that the  $2s_{1/2}$  state of hydrogen was shifted above the  $2p_{1/2}$  state, the Pauli approximation to the Dirac equation was used to decide that the Dirac equation was not sufficient. Thus, from a foundational view, it is important to understand the extent that we can believe the ineffectiveness of physical statements about the s-state behavior of the Dirac equation based on the Pauli equation (or any of its variants). In this paper, we have first recalled that a Dirac particle has a finite probability of being at the origin in the s-state. This is not true for p-states or higher so that the Pauli approximation applies. However, we show explicitly that the Pauli approximation is not valid for any calculation in s-states for two reasons:
  - (1) We cannot replace  $(E-V+mc^2)$  by  $2mc^2$  as V is undefined at r=0.
  - (2) The correct probability density for the Pauli approximation is  $\rho_{\psi} = |\psi(\mathbf{x})|^2 + \left|\frac{(\boldsymbol{\sigma} \cdot \boldsymbol{\pi})}{2mc} \psi(\mathbf{x})\right|^2$  but it cannot be used to compute s-state perturbations because it contains a  $r^{-2}$  term making it undefined at r = 0.

We then show that these difficulties can be avoided, with and error of ≈ 10−<sup>5</sup> ev, provided a cut-off at the classical electron radius is used. The important conclusion is that a cut-off is already required before QED. This implies that the extra covariant formulations and renormalization procedures may not be necessary to obtain and justify the excellent results of QED.

- 4.11. Square Root Operator. It is well-known that the square root operator is nonlocal in space but and explicit representation and it physical properties has not been available. In this paper we have provided such a representation, (due to Samko [55]). We have discussed its singularity properties and shown that it looks like point particle outside a few Compton wave lengths. The fact that it is also directly related to the Dirac operator by a unitary transformation (Foldy-Wouthuysen [20]) shows that mathematical equivalence does not always imply physical equivalence.
- 4.12. Feynman Path Integral. Historically, the mathematics community has had two responses to the introduction of a new mathematical idea or method into physics. The first response has been to fit the idea or method into an existing framework. The second and more exciting response is when such an idea or method leads to the development of a new branch of mathematics. The most prevalent and successful response has been in finding an existing mathematical structure that will reasonably accommodate the physical theory and provide (at least) the framework for mathematical rigor. In some rare but important instances, there is no obvious mathematical structure which can completely accommodate the theory in the manner presented by physicists. In this case, mathematicians have extended and/or adapted an existing mathematical theory, developed new mathematical structures or suggested (in frustration) that any conclusions derived from the use of these ideas or methods are at least suspect. Over the last seventy years, all of the above positions have appeared in response to Feynman?s introduction of his path integral into quantum theory.

In this paper, we have introduced a more general integral (HK) and designed a new Hilbert space (KS<sup>2</sup> [R n ]) for the Feynman theory. This approach provides the correct framework for mathematical rigor without removing the physically intuitive and computationally effective approach suggested by Feynman. This space is also unitarily equivalent to L 2 [R n ] showing that unitary equivalence is not always indicative of the same physics.

4.13. Feynman's Operator Calculus and Dyson' Conjectures. In his seminal paper [13], Dyson gave four conjectures concerning the foundations for QED. The last two were proven by Salam [2] and Weinberg [63] respectively. The first two remained open because Dyson used Feynman's time ordered operator calculus, which had no foundations. In this paper, we used the theory developed by Gill and Zachary to provide solutions showing that

Dyson was correct. In particular, we showed that the ultra-violet divergence is caused by a violation of the time-energy uncertainty relation at each point in time. Thus, a cut-off is also required by QED.

## Acknowledments

We would like to thank Professors Netsivi Ben-Amots, Alexander Gersten, Lawrence P. Horwitz, Martin C. Land, Elliot Leib and Ruggero M. Santilli for their continued interest and support. This paper is dedicated in honor of

Professor Lawrence P. Horwitz (Larry) on the occasion of his ninety-second birthday. The first author would like to register his deep admiration and affection for the many kind ways the Larry has freely given support, helpful suggestions and advice over the last 26 years. He has always conducted himself in a manner that represents the best model of what a scientific professional should be.

## Declaration. The authors certify that:

- (1) they have no relevant financial or non-financial interests to disclose;
- (2) they have no conflicts of interest to declare that are relevant to the content of this manuscript;
- (3) they have no financial or proprietary interests in any material discussed in this manuscript; and
- (4) they have no affiliations with or involvement in any organization or entity with any financial interest or non-financial interest in the subject matter or materials discussed in this manuscript.

## References

## References

- [1] S. Albeverio and S. Mazzucchi Feynman path integrals for polynomially growing potentials, J. Funct. Anal. 221 (2005), 83-121.
- [2] A. Salam, Overlapping divergence and the S-matrix, Phys. Rev. 82, (1951) 217-227.
- [3] L. C. M. Benjamin, R. M. Antonio and A. P. Gonzalo, The 2.7 ◦K blackbody radiation background reference frame, Chin. Phys. B 19 No. 4 (2010) 04XXXX1-5.
- [4] E. T. Bell, Men of Mathematics , Simon & Schuster, New York, (1937).
- [5] B. Javanmardi; C. Porciani; P. Kroupa; J. Pflamm-Altenburg Probing the Isotropy of Cosmic Acceleration Traced By Type Ia Supernovae The Astrophysical Journal Letters. 810 (1) (2015) 47.
- [6] H.R Brown, Eur. J. Phys. 26, S85 (2005).
- [7] B. Bakamjian and L. H. Thomas, Phys. Rev. 92, 1300 (1953).
- [8] C. D. Anderson, The Positive Electron, Physical Review 43, 491 (1943).
- [9] D. G. Currie, T. F. Jordan, and E. C. G. Sudarshan, Rev. Mod. Phys. 35, 350 (1963).
- [10] P. A. M. Dirac, Classical theory of radiating electrons, Proceedings of the Royal Soc. of London A 167 (1938) 148-168.
- [11] P. A. M. Dirac, Proc. Roy. Soc (London) A117 (1928) 610, A118 (1928) 351.
- [12] Einstein Centennial Symposium, October 1-5, 1979, University of Illinois at Carbondale, USA.

- [13] F. Dyson, The S-matrix in quantum electrodynamics, Phys. Rev. 75 (1949), 1736- 1755.
- [14] A. Einstein, Ann. d. Phys. 17, 891 (1905).
- [15] A. Einstein, Jahrbuch Radioaktivitat V, 422 (1907) (Berichtigungen).
- [16] J. R. Fanchi, Confronting the ENIGMA of TIME, World Scientific, Singapore, (2023).
- [17] R. P. Feynman, Phys. Rev. 81, 108 (1951).
- [18] R. P. Feynman and M. Gell-Mann, Phys. Rev. 109, 193 (1958).
- [19] M. Frisch, Inconsistency Asymmetry and Non-Locality, Oxford University Press N. Y. (2005).
- [20] L. L. Foldy and S. A. Wouthuysen, Phys. Rev. 78, 29 (1950).
- [21] J. Glimm and A. Jaffe, Quantum Physics. A Functional Integral Point of View, Springer, N. Y., (1987).
- [22] T. L. Gill and J. Lindesay, Canonical Proper Time Formulation of Relativistic Particle Dynamics, Int. J. of Theor. Phys. 32 (1993), 2087-2098.
- [23] T. L. Gill, The Square-Root Operator, Proper-Time and a Particle Interpretation for the Klein-Gordon Equation, FERMILAB-Pub-82/60-THY (1982).
- [24] T. L. Gill, and G. Ares de Parga https://iopscience.iop.org/article/10.1088/1742- 6596/2482/1/012015
- [25] T. L. Gill, and G. Ares de Parga https://iopscience.iop.org/article/10.1088/1742- 6596/1239/1/012013
- [26] I. S. Gradshteyn and I. M. Ryzhik Tables of Integrals, Series and Products, Academic Press, NewYork, (1980).
- [27] T. L. Gill and W. W. Zachary, Foundations for relativistic quantum theory I: Feynman's operator calculus and the Dyson conjectures, Journal of Mathematical Physics 43 (2002), 69-93.
- [28] T. L. Gill and W. W. Zachary, Functional Analysis and the Feynman operator Calculus, Springer New York, (2016).
- [29] T.L. Gill, W.W. Zachary, Analytic representation of the square-root operator, J. Phys. A: Math. Gen. 38 (2005) 1-18.
- [30] T.L. Gill, W.W. Zachary, and M. Alfred, Analytic representation of the Dirac equation, J. Phys. A: Math. Gen. 38 (2005) 1-22.
- [31] T.L. Gill, W.W. Zachary and J. Lindesay, The Classical Electron Problem, Found. Phys. 31, 1299 (2001).
- [32] https://www.spiedigitallibrary.org/conference-proceedings-of-spie/7421/1/Twoversions-of-Maxwells-equations-and-the-nature-of-light/10.1117/12.824191.
- [33] D. Hestenes, Space Time Algebra, Gordon and Breach Pub. New York, (1966)
- [34] L. P. Horwitz and C. Piron, Helv. Phys. Acta 46, 316 (1981).
- [35] E. Hille, and R. S. Phillips, Functional Analysis and Semigroups, American Mathematical Society Colloquium Publication No. 31 American Mathematical Society, Providence, RI, (1957).
- [HS] R. Henstock, The General Theory of Integration, Clarendon Press, Oxford, (1991).
- [HW] L. P. Horwitz, On the significance of a recent experiment demonstrating quantum interference in time, (to appear in Phys. Rev. Letters, see arXiv:quant-ph/0507044).
- [36] T. Ichinose, Essential selfadjointness of the Weyl quantized relativistic hamiltonian, Ann. Inst. Henri Poincar´e (Physique th´eorique) 51, (1989) 265-298.
- [KB] J. Kuelbs, Gaussian measures on a Banach space, Journal of Functional Analysis 5 (1970), 354–367.
- [KW] J. Kurzweil, Nichtabsolut konvergente Integrale, Teubner-Texte z¨ur Mathematik, Band 26, Teubner Verlagsgesellschaft, Leipzig, (1980).
- [37] H. A. Lorentz, Archives Neerlandaises des Sciences Exactes et Naturelles, 25, 353 (1892).
- [38] H. A. Lorentz, The Theory of Electrons B. G. Teubner, Leipzig, 1906; (reprinted by Dover, New York, 1952).

- [39] H. Leutwyler. A no-interaction theorem in classical relativistic hamiltonian particle mechanics, Nuovo Cim., 37, 556 (1965).
- [40] G. Longhi, L. Lusanna, and J. M. Pons, J. Math. Phys. 30, 1893 (1989).
- [41] H. Minkowski, Physikalische Zeitschrift 10,104 (1909).
- [42] Nathan J. Secrest; Sebastian von Hausegger; Mohamed Rameez; Roya Mohayaee; Subir Sarkar; Jacques Colin. A Test of the Cosmological Principle with Quasars, The Astrophysical Journal Letters, 908 (2) 2021): L51
- [43] A. Pais, Subtle is the lord: the science and life of Albert Einstein, Oxford University Press (1982).
- [44] P.J.E. Peebles, Principles of Physical Cosmology, Princeton University Press, London (1993).
- [45] W. Perret and G. B. Jeffery (translators, with additional notes by A. Sommerfeld), The Principle of Relativity by H. A. Lorentz, A. Einstein, H. Minkowski and H. Weyl, Dover, New York, 1952.
- [46] H. Poincar´e, Sur la dynamique de l' ´electron, Rendiconti del Circolo matematico Rendiconti del Circolo di Palermo, 21, 129-176 (1906).
- [47] W. K. H. Panofsky and M. Phillips, Classical Electricity and Magnetism, (Second Edition), Addision-Wesley, Reading, MA. (1962).
- [48] W. Pauli, Uber den Zusammenhang des Abschlusses der Elektronengruppen im Atom ¨ mit der Komplexstruktur der Spektren. Zeitschrift f¨ur Physik. 31 (1): 765-783 (1925).
- [49] W. Pauli, Handbuch der Physik 2nd ed. 24, (1933).
- [50] A. A. Penzias and R. W. Wilson, Ap. J. 142, 419 (1965).
- [51] M.H.L. Pryce, Proc. Roy. Soc. London A 195, 400 (1948).
- [52] G. Ares de Parga, T.L. Gill and W.W. Zachary, The Thomas program and the canonical proper-time theory J. of Comp. Methods in Sci. and Eng. 3 (2013) 117-134.
- [53] M. Radovan On the Nature of Time, arXiv:1509.01498
- [54] F. Rohrlich, Classical Charged Particles, Advanced Book Program, Addison-Wesley, New York (1990).
- [55] S. G. Samko, Hypersingular Integrals and Their Applications, Analytical Methods and Special Functions, Vol. 5 Taylor & Frances, Pub., (2002).
- [56] J. Schwinger, Found. Phys. 13, 2573 (1998).
- [57] Saadeh D, Feeney SM, Pontzen A, Peiris HV, McEwen, JD . How Isotropic is the Universe, Physical Review Letters. 117 (13) (2016) 131302.
- [58] J. C. Slater, Quantum Theory of Atomic Structure, Vol. II McGraw-Hill, New York (1960).
- [ST] V. Steadman, Theory of operators on Banach spaces, Ph.D thesis, Howard University, 1988.
- [STW] R. F. Streater and A. S. Wightman, PCT, Spin and Statistics and All That, Benjamin, New York, (1964).
- [59] S. Schweber, QED And The Men Who Made It , Princeton University Press, London (1994).
- [60] S. Walters, In: Gonner, H., Renn, J., Ritter, J. (eds.) The Expanding Worlds of General Relativity. Einstein Studies, vol. 7, pp. 45-86. Birkh¨auser, Boston (1999).
- [61] J.A. Wheeler and R.P. Feynman, Interaction with the absorber as the mechanism of radiation Rev. Mod. Phys. 17 (1949) 157-181.
- [62] A. O. Williams, Phys. Rev. 58 (1940), 723.
- [63] S. Weinberg, High energy behavior in quantum field theory, Phys. Rev. 118 (1960), 838-849.

Department of EECS, Mathematics and Computational Physics Laboratory, Howard University, Washington DC 20059 USA, tgill@howard.edu

Departmento de F´ısica, Escuela Superior de F´ısica y Matematicas, Insti- ´ tuto Politecnico Nacional COFAA, Edif 9, U. P. Adolfo L ´ øpez Mateos, Zaca- ´ tenco, Lindavista, 07738, Mexico D.F., M ´ exico, ´ E-mail : gadpau@hotmail.com