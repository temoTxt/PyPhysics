![](_page_0_Picture_1.jpeg)

### **PAPER • OPEN ACCESS**

# Mathematical Concepts in Physics

To cite this article: Tepper L. Gill 2025 J. Phys.: Conf. Ser. 2987 012003

View the [article online](https://doi.org/10.1088/1742-6596/2987/1/012003) for updates and enhancements.

## You may also like

- [Analysis of student's mathematical](/article/10.1088/1742-6596/1211/1/012107) [connection ability in linear equation system](/article/10.1088/1742-6596/1211/1/012107) [with two variables](/article/10.1088/1742-6596/1211/1/012107) D Rahmawati, Budiyono and D R S Saputro -
- [Exploration Ethnomathematics on](/article/10.1088/1742-6596/1776/1/012016) [Traditional House Ume Kbubu in North](/article/10.1088/1742-6596/1776/1/012016) [Central Timor Districts](/article/10.1088/1742-6596/1776/1/012016) Maritce Alfrida Tlonaen and Yohanis Ndapa Deda -
- [The effect of contextual teaching and](/article/10.1088/1755-1315/314/1/012064) [learning approach and motivation of](/article/10.1088/1755-1315/314/1/012064) [learning on the ability of understanding the](/article/10.1088/1755-1315/314/1/012064) [mathematics concepts of grade V student](/article/10.1088/1755-1315/314/1/012064) Yuvita Rama Yeni, Hendra Syarifuddin and Riska Ahmad -

doi:10.1088/1742-6596/2987/1/012003

# Mathematical Concepts in Physics

## Tepper L. Gill<sup>1</sup>

 $^1\mathrm{Department}$  of EECS, Mathematics and Computational Physics Laboratory, Howard University, Washington DC 20059 USA

E-mail: tgill@howard.edu

**Abstract.** The purpose of this paper show how mathematical concepts not designed for physics can muddle understanding and/or prevent progress. We provide the following examples to justify our position:

- 1. A proof that Minkowski's (third) postulate is incompatible with Einstein's first postulate for two or more particles.
- 2. A proof that the Dirac equation is not physically equivalent to the square root equation but is related by a unitary transformation.
- 3. The design of a Hilbert space for Feynman's path integral formulation of quantum theory  $KS^2[\mathbb{R}^n]$ , which is not physically equivalent to  $L^2[\mathbb{R}^n]$ , the Hilbert space for the Schödinger formulation, but is related by a unitary transformation.

A major contribution of this paper is that mathematical equivalence is not the same as physical equivalence and that unbiased physical analysis and justification of mathematical concepts should be required before acceptance into physical theory.

### 0 Introduction and Problem

The most important occurrence that led to the scientific revolution was not a fact but an idea. That is the general acceptance of the idea that we can understand the world and its workings by experiment (controlled experiences) and mathematical modeling. Following the important studies of Galileo Galilei, the first major outcome was the publication of the Principia by Isaac Newton in 1687. Newton is rightfully considered the father of theoretical physics. However, an equally important contribution was another idea, Newton's realization that mathematics is a self-consistent intellectual tool that can be molded (or designed) to provide a faithful representation of the world as he understood it. This unacknowledged contribution has had as much (if not more) impact on mankind's mathematical and scientific progress. From this perspective, mathematics is amazingly effective in science because much of it was designed for that purpose.

Before World War II, theoretical physicists worked on mathematical models that would explain experiments. Since then, the order has been reversed and theoretical physicists consider their task the construction of conceptual frameworks that can later be tested by experiment. This reversal has made it possible to ignore and even downplay difficult foundational problems. One outcome of this approach is that there are now 16 (mathematically beautiful) conceptual frameworks for the merge of quantum field theory with general relativity, with few, if any realizable experiments to test them (see Esposito). On the other hand, the first physically and mathematically consistent relativistic unification of Newton's mechanics with Maxwell's electrodynamics has only recently appeared (2024). This is over a hundred years after the problem was first identified (see [26]).

Content from this work may be used under the terms of the Creative Commons Attribution 4.0 licence. Any further distribution of this work must maintain attribution to the author(s) and the title of the work, journal citation and DOI.

doi:10.1088/1742-6596/2987/1/012003

## 0.1 Framework and Historical Background

There are five postulates that provide the foundations for all of physical science:

- (1) There is a real (external) physical world, which exists and is independent of our existence.
- (2) This physical world is as it appears to us in our consciousness.
- (3) This physical world is knowable to us via our senses (including our surrogate instruments).

All of physical science is rooted in these postulates. Our scientific forefathers no longer exist and yet the laws they left us still exist (in some form), so the first postulate is clear. Without the second we can never be sure that our senses give truth, making the third meaningless. Thus, we take these postulates to be obvious.

0.1.1 Experimental Physics Distortions and even illusions about the world can appear in the consciousness of one or more us and not appear to others (the collective). Thus, we need a collective physical filter that allows us to distinguish physical truth from individual or group illusion. We call this filter a scientific experiment.

Definition 0.1. A scientific experiment is a controlled reproducible experience. The basic postulate of experimental physics is that:

(4) We can obtain objective reproducible empirical information about the physical world by experimental investigation.

Definition 0.2. We say that a physical subject is known when all possible cause effect relations (of physical interest) can be a priori predicted under well-defined control conditions.

- 0.1.2 Theoretical Physics The objective of theoretical physics is to design and construct mathematical representations of the physical world. These representations must describe the cause effect relationships observed in experiment, must be physically and mathematically consistent, using a minimal number of variables and parameters. The basic postulate of theoretical physics is that:
- (5) Mathematics is the correct tool for the design, analysis and certification of the consistency of representations of physical reality. (The most distinct difference between mathematics and physics is that mathematics need not be concerned with real world objects.)

The Newton-Maxwell problem was faced directly by Abraham, Einstein, Lorentz, Poincar´e, Ritz and other major thinkers of the 1900's. Starting with a complete analysis of the microscopic behavior of electrons and ions, Lorentz was able to explain all the macroscopic properties of optics and electrodynamics (see [39], [40] ). He was the first to obtain the transformations between different observers, showing how their relative velocities affected their results when in different inertial frames. In his investigations, Poincar´e discovered an error in the derivation and, after correction, named them for Lorentz. In addition, Poincar´e observed that if time is treated as an imaginary coordinate, these Lorentz transformations formed a group. (He also derived the proper time, which was later claimed by Minkowski.)

Einstein's work on the photo-electric effect, give him a different perspective. As noted by Brown [6], he did not believe that Maxwell's theory would survive the existence of photons, but considered the Lorentz transformations fundamental. Thus, he derived them from kinematical arguments. (At that time, Maxwell's theory was still not accepted by all, see [33] and references therein.)

Einstein [14, 15] observed that the constant c appears in Maxwell's equations for all inertial observers. At that time experimental information about the speed of light was meager, being restricted to macroscopic and astronomical studies. Einstein was the first to realized that a formal postulate on the velocity of light was necessary. His proposal was that all physical theories should satisfy the (now well-known) postulates of special relativity:

- (1) The physical laws of nature and the results of all experiments are independent of the particular inertial frame of the observer (in which the experiment is performed).
- (2) The speed of light in empty space is constant and is independent of the motion of the source or receiver (in an inertial frame).

Like others before him, Einstein realized that the idea of absolute space was neither knowable nor necessary. The first postulate abandons this idea. Einstein's second postulate appears to abandon of the idea of absolute time.

doi:10.1088/1742-6596/2987/1/012003

Unlike Minkowski, Poincaré's insight and understanding of the difference between physics and mathematics helped him to resist the temptation to use the "physically unjustified" mathematical observation that time be made a forth geometric coordinate as a (necessary) tool for the representation of physical reality. Minkowski made this leap later with much philosophical fanfare, but lacking any physical justification. Minkowski was a number theorist with few accomplishments of note in physics. In the work of Walters, we find a clear discussion of Minkowski's physics background and his knowledge of Poincaré's work (see [62]).

Thus, we make explicit Minkowski's unacknowledged third postulate for Einstein's special theory of relativity:

(3) The correct implementation of the first two postulates requires that time be treated as a fourth coordinate, and the relationship between components so constrained as to satisfy the natural invariance induced by the Lorentz group of electrodynamics, (Minkowski space).

The physics community voted in favor of the Lorentz group as the proper transformation theory. Today, many believe that Newtonian theory is covered by Minkowski's version of the special theory. As will be shown later, the failure to conduct a complete (physical) analysis of the third postulate is at the root of all major problems in the classical foundations.

#### 0.2 Summary

In the first section we show that there are three representations for the proper time which follow from Einstein's special theory. The first representation is well-known and is used to prove that Minkowski's postulate is incompatible with the first postulate of Einstein for two or more particles. We then discuss the center of mass problem due to Pryce [53], the many particle quantum problem due to the Bakamjian and Thomas [7] and the no-interaction theorem first proven by Currie, Jordan, and Sudarshan [9]. All of these foundational problems are a direct consequence Minkowski's postulate.

In the second section, we used the second representation of proper time and a number of important examples, to explicitly show the danger in equating mathematical equivalence to physical equivalence, and the explicit need for mathematics designed for physics. We prove that the Dirac and square root equations are not physically equivalent, and yet (mathematically) unitary equivalent. As an aside, we provide the correct physical interpretation of the well-known zitterbewegung first observed by Schödinger in the Dirac equation.

In the third section, we introduce the Henstock-Kurzweil integral (HK). It is very close to the one learned in calculus but extends the Lebesgue integral to include nonabsolutely integrable functions (like  $e^{ix^2/2}$ ). We then design the natural Hilbert space for HK integrals, which allows us to construct the Feynman path integral in the manner he intended, while retaining both the physical intuition and computational advantage of his approach. This space  $KS^2[\mathbb{R}^n]$  contains  $L^2[\mathbb{R}^n]$  and the Schwartz distributions or generalized functions as continuous dense embeddings. From a physical point of view, use of the HK-integral and  $KS^2[\mathbb{R}^n]$  eliminates the extra time and effort required to first learn Lebesgue measure theory and then the theory of distributions.

In section four we use the third representation of proper time to prove the existence of a universal clock first predicted by Horwitz and Fanchi, which we call the historical time (see [18] and [36]). We use it to construct Feynman's time ordered operator calculus, construct his path integral, and to prove the last two remaining conjectures of Dyson for the foundations of QED.

## 1 Proper time I and Minkowski's Postulate

Let m be the mass of a particle or the effective mass for a system of particles. We assume the particle or system is defined on phase space with variables  $(t, \mathbf{x}, \mathbf{p})$  as seen by observer O, and variables  $(t', \mathbf{x}', \mathbf{p}')$  as seen by observer O' (both inertial frames). Where  $\mathbf{x}, \mathbf{x}'$  is the position of the particle (or center of mass) and  $(\mathbf{p}, \mathbf{p}')$  is the particle momentum (or center of mass momentum). If  $\gamma^{-1}(\mathbf{v}) = \sqrt{1 - \frac{\mathbf{v}^2}{c^2}}$ ,  $\mathbf{v} = \frac{d\mathbf{x}}{dt}$ , where  $\mathbf{v}$  ( $\mathbf{v}'$ ) is the velocity of the particle, the first representation of the proper time is defined by our observers as:

$$d\tau = \gamma^{-1}(\mathbf{v}) dt, \quad d\tau = \gamma^{-1}(\mathbf{v}') dt' \quad \text{proper time 1.}$$
 (1.1)

doi:10.1088/1742-6596/2987/1/012003

For the second representation, first square (1.1) to get:

$$d\tau^2 = dt^2 - \frac{1}{c^2}d\mathbf{x}^2, \quad d\tau^2 = dt'^2 - \frac{1}{c^2}d\mathbf{x}'^2.$$

Reaarrange terms to obtain:

$$d\tau^2 + \frac{1}{c^2}d\mathbf{x}^2 = dt^2, \quad d\tau^2 + \frac{1}{c^2}d\mathbf{x'}^2 = dt'^2 \Rightarrow$$

$$\sqrt{1 + \frac{\mathbf{u}^2}{c^2}}d\tau = dt, \quad \mathbf{u} = \frac{d\mathbf{x}}{d\tau}, \quad \sqrt{1 + \frac{\mathbf{u'}^2}{c^2}}d\tau = dt', \quad \mathbf{u'} = \frac{d\mathbf{x'}}{d\tau}.$$

If we let  $b = \sqrt{c^2 + \mathbf{u}^2}$ ,  $(b' = \sqrt{c^2 + \mathbf{u'}^2})$ , then we have our second representation

$$cdt = bd\tau$$
 and  $cdt' = b'd\tau$  proper time 2. (1.2)

For the third representation, use the relations  $H = \gamma(\mathbf{v}) mc^2$  for observer O and  $H' = \gamma(\mathbf{v}') mc^2$  for observer O', to obtain:

$$d\tau = \frac{mc^2}{H}dt$$
,  $d\tau = \frac{mc^2}{H'}dt'$  proper time 3. (1.3)

#### 1.1 The Minkowski Problem

Minkowski used the first representation of proper time as the metric for his four-spacetime theory. In order to use his postulate, Minkowski introduced the clock of a co-moving observer as a substitute metric (see the notes in Sommerfeld in [47]). This substitute is mathematically convenient but should have alarmed all physics thinkers, since no co-moving observer can actually sit on the tangent plane of a bullet.

Einstein, Lorentz, Poincar?e, and Ritz, regarded it as a mathematical abstraction lacking physical content and maintained that space and time have distinct physical properties. In fact, Einstein was the first (and only one) to oppose the Minkowski openly. Einstein and Jakob Laub wrote two papers that introduced a different approach to electromagnetic theory, which was simpler and did not depend on the spacetime formalism (see [16, 17]). They argued that the spacetime formalism was too complicated and did not add any new physics.

The idea of a four dimensional space was very attractive to mathematicians who entered and dominated the field for at least a decade. In the eyes of everyone, this also made Minkowski's postulate a natural part of Einstein's special theory. As a consequence, the intense physical investigation and justification normally given for new postulates did not happen. By the time problems in attempts to develop a many particle relativistic quantum theory appeared, Minkowski's postulate was already sacred. For the general theory, Einstein only considered an approach that extended Minkowski's postulate (see Pais [45]). The following theorem shows how important an early investigation would have been.

Theorem 1.1. (Minkowski Incompatible Theorem) The addition of Minkowski's postulate to the postulates of Einstein is incompatible for two or more particles,

*Proof.* Let O and O' be two inertial observers. Without loss, we can assume that both clocks begin when their origins coincide and O' is moving with uniform velocity  $\mathbf{v}$  as seen by O. Let two particles, each the source of an electromagnetic field, move with velocities  $\mathbf{w}_i$  (i = 1, 2), as seen by O, and  $\mathbf{w}'_i$  (i = 1, 2), as seen by O', such that:

$$\mathbf{x}'_{i} = \mathbf{x}_{i} - \gamma(\mathbf{v}) \mathbf{v}t + [\gamma(\mathbf{v}) - 1] \left(\mathbf{x}_{i} \cdot \mathbf{v} / \|\mathbf{v}\|^{2}\right) \mathbf{v},$$

$$\mathbf{x}_{i} = \mathbf{x}'_{i} + \gamma(\mathbf{v}) \mathbf{v}t' + [\gamma(\mathbf{v}) - 1] \left(\mathbf{x}'_{i} \cdot \mathbf{v} / \|\mathbf{v}\|^{2}\right) \mathbf{v} \quad \text{and},$$

$$(1.4)$$

$$\mathbf{w'}_{i} = \mathbf{w}_{i} - \gamma(\mathbf{v})\mathbf{v} + \left[\gamma(\mathbf{v}) - 1\right] \left(\mathbf{w}_{i} \cdot \mathbf{v} / \|\mathbf{v}\|^{2}\right)\mathbf{v},$$

$$\mathbf{w}_{i} = \mathbf{w'}_{i} + \gamma(\mathbf{v})\mathbf{v} + \left[\gamma(\mathbf{v}) - 1\right] \left(\mathbf{w'}_{i} \cdot \mathbf{v} / \|\mathbf{v}\|^{2}\right)\mathbf{v}.$$
(1.5)

doi:10.1088/1742-6596/2987/1/012003

Thus, there is no problem in requiring the positions and velocities to transform as expected. However, when we try to transform the clocks, we see the problem at once because

$$t' = \gamma(\mathbf{v}) \left( t - \mathbf{x}_1 \cdot \mathbf{v}/c^2 \right)$$
 and  $t' = \gamma(\mathbf{v}) \left( t - \mathbf{x}_2 \cdot \mathbf{v}/c^2 \right)$ . (1.6)

This is impossible unless  $\mathbf{x}_1 \cdot \mathbf{v} = \mathbf{x}_2 \cdot \mathbf{v}$  (for all time). This means that the two observers cannot use their clocks to share information (with other observers) about two or more particles. It follows that, if we attempt to replace  $\mathbf{x}_i$  and  $\mathbf{x}'_i$  with four vectors using t and t', the first postulate fails. Thus, these three postulates are both mathematically and physically incompatible.

**Remark 1.2.** The Lorentz transformations (equations (1.4) (1.5) and (1.6)) contain information about each observer and each particle, so it is not clear that they can be converted into the tensor notation that came into fashion after Sommerfeld simplified the (truly) complicated notation of Minkowski.

1.1.1 The *n*-particle position problem To construct the many-particle clock for the O frame observer, assume that n interacting particles have Hamiltonians  $H_i$  and a total Hamiltonian  $H = \sum_{i=1}^n H_i$ . The effective mass energy  $Mc^2$ , and total momentum  $\mathbf{P}$ , are defined as follows:

$$Mc^2 = \sqrt{H^2 - c^2 \mathbf{P}^2}, \quad \mathbf{P} = \sum_{i=1}^n \mathbf{p}_i.$$

We can now also represent the Hamiltonian by  $H = \sqrt{c^2 \mathbf{P}^2 + M^2 c^4}$ . In 1948 Pryce conducted the first study of the relativistic center-of-mass for two or more particles (see [53]). He found only one was canonical and available to all the observers. He observed that the canonical representation of the center-of-mass cannot be the three-vector part of a four-vector. This problem arises because the canonical center of mass position  $\mathbf{X}$  is defined by: (For a modern representation see Longhi, Lusanna, and Pons [42].)

$$\mathbf{X} = \frac{1}{H} \sum_{i=1}^{n} H_i \mathbf{x}_i + \frac{c^2 \left( \mathbf{S} \times \mathbf{P} \right)}{H \left( Mc^2 + H \right)}, \tag{1.7}$$

where S is the global spin of the system of particles relative to O. (It is clear that (1.7) cannot represent the vector part of a four-vector.) However, X is the canonical conjugate of P and this is precisely what is needed for a many-particle relativistic quantum theory. (It is clear that the problem is the third postulate.)

Five years after Pryce, Bakamjian and Thomas constructed a many-particle relativistic quantum theory using the canonical center of mass X. They explicitly showed that their theory satisfied the two postulates of Einstein (see [7]). They conjectured that, any attempt to enforce the third postulate would only be compatible with free particles (BT-conjecture).

In 1965 Currie, Jordan and Sudarshan [9], gave a proof of the BT-conjecture for two particles, which they renamed "The No-Interaction Theorem". It has since been extended to an arbitrary number of particles by Leutwyler [41].

**Theorem 1.3.** (No-Interaction Theorem) Suppose we have a system of n particles with phase space variables  $\{(\mathbf{x}_i, \mathbf{p}_i)\}_{i=1}^n$  defined on  $\mathbb{R}^{3n} \times \mathbb{R}^{3n}$  with the following properties:

- 1. The system has a Hamiltonian representation.
- 2. The system has a canonical representation of the Poincaré group.
- 3. Each  $\mathbf{x}_i$  is the vector part of a four-vector.

Then these assumptions are only compatible with free particles.

Thus, we see that Minkowski's postulate is the culprit, which has hindered progress in both classical and quantum mechanics for over a 120 years. The first hint at a solution was given by Feynman when he showed that, by letting time recover its role in physics as it appears in our consciousness, he could construct his operator calculus, his formulation of quantum mechanics and QED.

doi:10.1088/1742-6596/2987/1/012003

#### 2 Proper Time II and Quantum Problems

In this section, we first show how the second definition of proper time avoids the problems created by Minkowski's postulate. We the analyze the Dirac equation and show why the Pauli approximation cannot be used to study s-states in Hydrogen. We then show that although related by a unitary transformation, the Dirac operator and the square root operator are not physically equivalent.

#### 2.1 One-Particle Clock

The proper time for each particle is uniquely defined for each observer. Assume a system of *n*-particles observed by O, who is able to identify each particle and attach a vector  $\mathbf{x}_i$  to the  $i^{\text{th}}$  particle, denoting its spacial distance to the origin. If  $\mathbf{w}_i$  is the velocity of particle i as seen by O, let  $\gamma^{-1}(\mathbf{w}_i) = \sqrt{1 - \mathbf{w}_i^2/c^2}$ . The  $i^{\text{th}}$  particle proper time is defined as before by:

$$d\tau_i = \gamma^{-1}(\mathbf{w}_i)dt, \quad \mathbf{w}_i = \frac{d\mathbf{x}_i}{dt}, \quad d\tau_i^2 = dt^2 - \frac{1}{c^2}d\mathbf{x}_i^2.$$
 (2.1)

Rewrite the last term to get:

$$dt^{2} = d\tau_{i}^{2} + \frac{1}{c^{2}}d\mathbf{x}_{i}^{2}, \Rightarrow cdt = \left(\sqrt{\mathbf{u}_{i}^{2} + c^{2}}\right)d\tau_{i}, \quad \mathbf{u}_{i} = \frac{d\mathbf{x}_{i}}{d\tau_{i}} = \gamma(\mathbf{w}_{i})\mathbf{w}_{i}. \tag{2.2}$$

If we let  $b_i = \sqrt{\mathbf{u}_i^2 + c^2}$ , the second term in equation (2.2) becomes  $cdt = b_i d\tau_i$ . This leads to a new identity:

$$\frac{1}{c}\frac{d}{dt} \equiv \frac{1}{b_i}\frac{d}{d\tau_i}, i = 1, \dots, n. \tag{2.3}$$

Thus, we obtain an interesting set of relationships between the proper time of each particle and and that of the observer. Applying it to  $\mathbf{x}_i$ , we obtain another set of relationships showing that the tangent space is left invariant in all cases.

$$\frac{\mathbf{w}_i}{c} = \frac{1}{c} \frac{d\mathbf{x}_i}{dt} \equiv \frac{1}{b_i} \frac{d\mathbf{x}_i}{d\tau_i} = \frac{\mathbf{u}_i}{b_i}, i = 1, \dots, n.$$
(2.4)

The new particle coordinates are  $(\mathbf{x}_i, \tau_i)$ . In this representation,  $\mathbf{x}_i$  is uniquely defined relative to O, while  $\tau_i$  is uniquely defined by each particle.

**Remark 2.1.** Equation (2.4) implies that the use of  $\beta = \mathbf{w}_i/c \equiv \mathbf{u}_i/b_i$  as a parameter for measurements is no-longer reliable and that time dilation is a physically incorrect concept.

Since the momentum is  $\mathbf{p}_i = m_i \gamma(\mathbf{w}_i) \mathbf{w}_i$  and  $\mathbf{u}_i = \gamma(\mathbf{w}_i) \mathbf{w}_i$ , we can also represent the momentum as  $\mathbf{p}_i = m_i \mathbf{u}_i$ . It follows that the particle phase space is also left invariant. As  $m_i$  is the rest mass, we see that the concept of relativistic mass increase is also is a physically incorrect concept.

**Theorem 2.2.** If the observer time is replaced by each particle's proper time, the Minkowski incompatibility theorem no-longer holds and the transformations that preserve the first postulate are:

$$b_i' = \gamma(\mathbf{v}) (b_i - \mathbf{u}_i \cdot \mathbf{v}/c)$$

and

$$b_i = \gamma(\mathbf{v}) (b'_i + \mathbf{u}'_i \cdot \mathbf{v}/c)$$

*Proof.* If we use each particle proper time then, for observer O'

$$\frac{\mathbf{w}_i'}{c} = \frac{1}{c} \frac{d\mathbf{x}_i'}{dt} \equiv \frac{1}{b_i'} \frac{d\mathbf{x}_i'}{d\tau_i} = \frac{\mathbf{u}_i'}{b_i'}, i = 1, \dots, n.$$

From  $t' = \gamma(\mathbf{v}) \left( t - \mathbf{x}_i \cdot \mathbf{v}/c^2 \right)$  we have  $cdt' = \gamma(\mathbf{v}) \left( cdt - d\mathbf{x}_i \cdot \mathbf{v}/c \right)$ . Using the two identities  $cdt' = b_i' d\tau_i$  and  $cdt = b_i d\tau_i$ , we obtain:

$$b_i' d\tau_i = \gamma (\mathbf{v}) (b_i d\tau_i - d\mathbf{x}_i \cdot \mathbf{v}/c)$$

Dividing by  $d\tau_i$ , we see that the transformations between  $b_i$  and  $b'_i$  are

$$b_i' = \gamma(\mathbf{v}) (b_i - \mathbf{u}_i \cdot \mathbf{v}/c)$$
.

doi:10.1088/1742-6596/2987/1/012003

By a similar calculation, we obtain the reverse relationship:

$$b_i = \gamma (\mathbf{v}) (b'_i + \mathbf{u'}_i \cdot \mathbf{v}/c).$$

It now follows that the two postulates are valid for all observers.

Each observer can now use their clock and the canonical center of mass to quantize the n-particle system (following [7]).

## Remark 2.3. We have three important comments:

- 1. If O and O′ are at rest relative to each other then v = 0, γ(v) = 1 and b<sup>i</sup> = b ′ i .
- 2. If particle i is at rest in O, then τ<sup>i</sup> = t and b<sup>i</sup> = c but τ<sup>i</sup> ̸= t ′ and b ′ <sup>i</sup> > c.
- 3. Since the particle frames need not be inertial, b<sup>i</sup> > c does not violate the second postulate.

## 2.2 Quantum Problems

This section is devoted to the discussion a number of missed and misunderstood mathematical issues in the foundations of relativistic quantum theory.

2.2.1 The Dirac Equation After Palui proposed his exclusion principle in 1925 ([50] ), Dirac was able to find the first physically viable way around the square-root problem in 1926 (see [11]). Raised in the British tradition of quaternions, Dirac knew that quaternions could be used to generate a spin-1/2 representation for the rotation group. He used the Pauli matrices to write the square-root operator as:

$$\sqrt{c^2 \mathbf{p}^2 + m^2 c^4} = \sqrt{\left[c\alpha \cdot \mathbf{p} + \beta mc^2\right]^2}, \quad \text{where,} \quad \beta = \begin{bmatrix} I_2 & 0\\ 0 & -I_2. \end{bmatrix}$$
 (2.5)

and the matrix α is defined by α = (α1, α2, α3), with

$$\alpha_i = \begin{pmatrix} \mathbf{0} & \sigma_i \\ \sigma_i & \mathbf{0} \end{pmatrix}, \ \ \sigma_1 = \begin{pmatrix} \mathbf{0} & 1 \\ 1 & \mathbf{0} \end{pmatrix}, \ \ \sigma_2 = \begin{pmatrix} \mathbf{0} & -i \\ i & \mathbf{0} \end{pmatrix}, \ \ \sigma_3 = \begin{pmatrix} 1 & \mathbf{0} \\ \mathbf{0} & -1 \end{pmatrix}$$

He thus obtained an alternative representation of equation (2.5), now known as the Dirac equation:

$$i\hbar \frac{\partial \Psi}{\partial t} = \left[c\boldsymbol{\alpha} \cdot \mathbf{p} + \boldsymbol{\beta} mc^2\right] \Psi$$
 (2.6)

We must now view Ψ as a vector-valued function or spinor. That is, Φ ∈ L 2 (R 3 , C<sup>4</sup> ) = L 2 (R 3 ) ⊗ C 4 is a four-component column vector Ψ = (ψ1, ψ2, ϕ1, ϕ2) t . In this representation, ψ = (ψ1, ψ2) t is the particle component and ϕ = (ϕ1, ϕ2) t is the antiparticle component. Dirac's equation gave the correct spin value for the electron and the then known (at that time) spectrum for Hydrogen. His equation also predicted negative energy solutions, which led Dirac to propose his hole theory of the vacuum. In this view, a hole is an electron with negative energy and positive charge. This means that Dirac's theory cannot be single particle one. However, the later discovery of the positron confirmed his approach. Schweber [61] provides a good discussion of the difficulties, which eventually led to QED.

A complete understanding of the Dirac equation only became possible after the development of geometric algebra. We now see that Dirac used a Lorentz-invariant Clifford algebra to represent the complex algebra of observables for the electron (see Hestenes [35]). However, it is apparently believed that, when interaction is present a complete separation of the particle and antiparticle components of his equation without approximation is not possible. The algebraic researches of Foldy-Wouthuysen [22], Feynman and Gell-Mann [20] and Pauli [51], along with the various approximations approaches in the literature may have supported this idea.

The first exact analytical separation of Dirac's equation with minimal coupling was found by Gill, Zachary and Alfred in 2005. Since this is important for our physical interpretation of Dirac's equation, the correct physical interpretation of the zitterbewegung, the correct physical relation between Dirac's equation and the square root equation and, to understand the physical failure of the Palui approximation for the study of s-states in Hydrogen, we provide a brief outline of this method. (See [32] for details.)

doi:10.1088/1742-6596/2987/1/012003

If  $\mathbf{A}(\mathbf{x},t)$  and  $V(\mathbf{x})$  are vector and scalar potentials make the standard identifications  $\mathbf{p} \to \boldsymbol{\pi} = \mathbf{p} - (e/c)\mathbf{A}$  and write (2.6) in two-component form:

$$i\hbar \frac{\partial \psi}{\partial t} - (mc^2 + V)\psi = c(\sigma \cdot \pi)\varphi,$$
 (2.7)

and

$$i\hbar \frac{\partial \varphi}{\partial t} + (mc^2 - V)\varphi = c(\boldsymbol{\sigma} \cdot \boldsymbol{\pi})\psi.$$
 (2.8)

Where the vectors  $\boldsymbol{\psi} = [\psi_1, \psi_2]^t$  and  $\boldsymbol{\varphi} = [\varphi_1, \varphi_2]^t$  represent the particle and anti-particle components  $\Psi = [\boldsymbol{\psi}, \boldsymbol{\varphi}]^t$ . We can separate these two equations using the method of integrating factors as in ODE's. We solve for  $\boldsymbol{\varphi}$  in the second equation and use the result in the first to obtain an equation which only depends on  $\boldsymbol{\psi}$ . Thus, we need to determine  $\delta$  such that:

$$e^{-\delta t}i\hbar \frac{\partial}{\partial t} \left[ \varphi e^{\delta t} \right] = i\hbar \frac{\partial \varphi}{\partial t} + \left( mc^2 - V \right) \varphi$$

Simplify the left-hand side and compare to obtain that  $\delta = \frac{1}{i\hbar} (mc^2 - V)$ . We now assume that  $\psi$  and  $\varphi$  vanish at  $-\infty$  and integrate to obtain:

$$i\hbar \frac{\partial \psi}{\partial t} - (mc^2 + V) \psi = c^2 (\boldsymbol{\sigma} \cdot \boldsymbol{\pi}) \left\{ \int_{-\infty}^t \exp\left\{\frac{i}{\hbar} (mc^2 - V) (t - \tau)\right\} (\boldsymbol{\sigma} \cdot \boldsymbol{\pi}) \psi(\mathbf{x}, \tau) d\tau \right\}.$$
(2.9)

Using the same method, we obtain our second equation for  $\varphi$ :

$$i\hbar \frac{\partial \varphi}{\partial t} + (mc^2 - V) \varphi = c^2 (\sigma \cdot \pi) \left\{ \int_{-\infty}^t \exp\left\{ -\frac{i}{\hbar} (mc^2 + V) (t - \tau) \right\} (\sigma \cdot \pi) \varphi(\mathbf{x}, \tau) d\tau \right\}.$$
(2.10)

Equations (2.9) and (2.10) now decomposed  $L^2(\mathbb{R}^3, C^4)$  into a direct sum  $L^2(\mathbb{R}^3, C^4) = L^2(\mathbb{R}^3, C^2) \oplus L^2(\mathbb{R}^3, C^2)$ . For later use, the probability density for  $\psi$  is:

$$\rho_{\psi}(t) = \left| \psi\left(\mathbf{x}, t\right) \right|^{2} + \left| c \int_{-\infty}^{t} \exp\left\{ \frac{i}{\hbar} \left( mc^{2} - V \right) (t - \tau) \right\} (\boldsymbol{\sigma} \cdot \boldsymbol{\pi}) \psi\left(\mathbf{x}, \tau\right) d\tau \right|^{2}.$$
 (2.11)

Thus, we obtain a complete analytic diagonalization of the particle, antiparticle wave functions (see [32]). The decomposition is analytic, compared to Foldy and Wouthuysen's group theoretical framework in that the wave function is not transformed (see [22]). We will discuss this in further detail, after we study the square root operator.

2.2.2 The Square Root Equation In 2005 Gill and Zachary [31] used the theory of fractional powers of linear operators to construct an (analytic) representation for the square-root energy operator valid for all spin values. A better representation was found by Samko (see [57] pg. 270):

$$\beta \hbar c \left[ \sqrt{\omega^{2} I - \Delta} \right] \psi (\mathbf{x}, t) = \beta m c^{2} \psi (\mathbf{x}, t)$$

$$- \frac{\beta \hbar c}{2\pi^{2}} \left\{ \int_{\mathbb{R}^{3}} \frac{K_{2} \left[ \omega \left| \mathbf{x} - \mathbf{y} \right| \right]}{\left| \mathbf{x} - \mathbf{y} \right|^{2}} \left[ \psi (\mathbf{x}, t) - \psi (\mathbf{y}, t) \right] d\mathbf{y} \right\}.$$
(2.12)

The integral is over  $\mathbb{R}^3$  showing in what sense the square root operator is nonlocal in space. For a better understanding of this extended object, we need more information about the Bessel functions  $K_2[\omega|\mathbf{y}|]$ ,  $K_1[\omega|\mathbf{y}|]$ ,  $K_{\frac{1}{2}}[\omega|\mathbf{y}|]$  and  $K_0[\omega|\mathbf{y}|]$  (see Gradshteyn and Ryzhik [28]).

$$\begin{split} K_2\left[\omega\left|\mathbf{y}\right|\right] &= K_0\left[\omega\left|\mathbf{y}\right|\right] + \frac{K_1\left[\omega\left|\mathbf{y}\right|\right]}{\left[\omega\left|\mathbf{y}\right|\right]} \quad \text{and} \\ \frac{K_{1/2}\left[\omega\left|\mathbf{y}\right|\right]}{\left[\omega\left|\mathbf{y}\right|\right]^{1/2}} &= \sqrt{\frac{\pi}{2}} \frac{\exp\left\{-\omega\left|\mathbf{y}\right|\right\}}{\left[\omega\left|\mathbf{y}\right|\right]}. \end{split}$$

doi:10.1088/1742-6596/2987/1/012003

Since the integral of  $|\mathbf{x} - \mathbf{y}|^{-2}$  is finite, the effective singularity properties for  $0 < z \ll 1$  are:

$$\frac{\frac{K_{1}[z]}{z}}{\frac{z}{z^{1/2}}} = C_{1} \left[ 1 + \theta_{1} \left( z \right) \right] z^{-2}$$

$$\frac{\frac{K_{1/2}[z]}{z^{1/2}}}{z^{1/2}} = \sqrt{\frac{\pi}{2}} \left[ 1 + \theta_{2} \left( z \right) \right] z^{-1}$$

$$K_{0}[z] = C_{0} \left[ 1 + \theta_{0} \left( z \right) \right] \ln z^{-1}$$

The  $K_1$  diverges at 0 like  $z^{-2}$  (strong interaction?) and the  $K_0$  is actually an integrable singularity. We include the  $K_{1/2}$  for comparison as it is essentially the well-known Yakawa potential (expected to account for the short range nuclear interaction). The effective properties for  $z \gg 1$  are:

$$\frac{K_{1}[z]}{z} = C_{1} \left[ 1 + \Theta_{1} \left( z \right) \right] z^{-3/2} e^{-z} \\ \frac{K_{1/2}[z]}{z^{1/2}} = \sqrt{\frac{\pi}{2}} \left[ 1 + \Theta_{1/2} \left( z \right) \right] z^{-1} e^{-z} \\ K_{0} \left[ z \right] = C_{0} \left[ 1 + \Theta_{0} \left( z \right) \right] z^{-1/2} e^{-z} \\ \right\}.$$

In this region, we have a cut-off range and the strength of interaction reverses, with the  $K_0$  having the longest tail (weak interaction?). It should be clear that the square root operator is nonlocal in space and can be treated like the identity outside a few Compton wavelengths.

2.2.3 Discussion I The important conclusion from the last two sections is that the Dirac and square root operators are not physically equivalent. The Dirac operator is nonlocal in time while the square root operator is nonlocal in space. It is shown in [32] that the correct physical interpretation of the zitterbewegung and the fact that the expected instantaneous value of a velocity measurement of a Dirac particle is  $\pm c$  is that it instantaneously jumps between +c and -c at each point in time to make it appear as a point in space.

We must now account for the known fact that the Foldy-Wouthuysen is a unitary transformation mapping the coupled particle-antiparticle Dirac operator to the uncoupled particle-antiparticle square root operator. From equations (2.9) and (2.10), we see that the Foldy-Wouthuysen unitary operator maps a nonlocal in time operator to a nonlocal in space operator, so that the unitary property of an operator is not sufficient to give confidence in the physics. We should also remember that the concept of unitary equivalence is a purely mathematical one, taken over to physics without any physical justification.

**Remark 2.4.** The Schödinger representation (on  $L^2[\mathbb{R}^3]$ ) and the Heisenberg representation  $\ell_2[\mathbb{R}^3]$  were known to be physically equivalent before they were proven to also be mathematical equivalent (i.e, related by a unitary mapping) by von Neumann. However, there are an infinite number of Hilbert spaces related to  $\ell_2[\mathbb{R}^3]$  by a unitary mapping, which have nothing to due with quantum theory. In a physically justified sense,  $L^2[\mathbb{R}^3]$  is not the correct Hilbert space for Feynman's path integral formulation of quantum mechanics. We will construct it in the next section. This space will also serve as another counter example to the notion that mathematical equivalence means physical equivalence.

## 3 Hilbert Space for Feynman's Formulation of Quantum Mechanics

In this section, we provide an introduction to the Henstock-Kurzweil integral (HK). The integral is well defined for operator-valued functions that may not be separably valued (where both the Bochner and Pettis integrals are undefined). Intuitively, one uses a version of the Riemann integral (of calculus) with the interior points chosen first, while the size of the base rectangle around any interior point is determined by an arbitrary positive function defined at that point. This integral was discovered independently by Henstock and Kurzweil. Our aim is to provide a sense of the conceptual and technical simplicity of the HK-integral. Those interested in more detail and proofs are directed to Gill and Zachary [30]. Let  $\mathcal{H}$  be a separable Hilbert space and let  $L(\mathcal{H})$  be the algebra of bounded linear operators on  $\mathcal{H}$ . Let  $[a,b] \subset \mathbb{R}$  and, for each  $t \in [a,b]$ , let  $H(t) \in L(\mathcal{H})$  be a given family of operators.

**Definition 3.1.** Let  $\delta(t)$  map  $[a,b] \to (0,\infty)$ , and let  $\mathbf{P} = \{t_0, \tau_1, t_1, \tau_2, \cdots, \tau_n, t_n\}$ , where  $a = t_0 \leqslant \tau_1 \leqslant t_1 \leqslant \cdots \leqslant \tau_n \leqslant t_n = b$ . We call  $\mathbf{P}$  a HK-partition for  $\delta$  (or HK-partition when  $\delta$  is understood) provided that for  $0 \leqslant i \leqslant n-1$ ,  $t_i, t_{i+1} \in (\tau_{i+1} - \delta(\tau_{i+1}), \tau_{i+1} + \delta(\tau_{i+1}))$ .

**Definition 3.2.** The family H(t),  $t \in [a,b]$ , is said to have a (uniform) HK-integral if there is an operator Q[a,b] in  $L(\mathcal{H})$  such that, for each  $\varepsilon > 0$ , there exists a function  $\delta$  from  $[a,b] \to (0,\infty)$  such that, whenever  $\mathbf{P}$  is a HK-partition for  $\delta$ , then

$$\left\| \sum_{i=1}^{n} \Delta t_i H(\tau_i) - Q[a, b] \right\| < \varepsilon.$$

doi:10.1088/1742-6596/2987/1/012003

In this case, we write

$$Q[a,b] = (HK) \int_{a}^{b} H(t)dt.$$

**Example 3.3.** The following example shows there are HK-integrable functions, which are not Lebesgue integrable. If  $F'(t) = 2t\sin(1/t^2) - (2/t)\cos(1/t^2)$ , for  $t \in (0,1)$  and F'(0) = 0. It is easy to see that  $F(t) = t^2\sin(1/t^2)$ , so that

$$(HK)$$
  $\int_{0}^{1} \left( 2t \sin \frac{1}{t^2} - 2\frac{1}{t} \cos \frac{1}{t^2} \right) dt = \sin 1.$ 

The following theorem explains the difference between the two integrals.

**Theorem 3.4.** Let  $f(t): [a,b] \to \mathbb{R}$ .

- 1. If f(t) is L-integrable on [a,b], then it is HK-integrable on [a,b] and:  $HK-\int_a^b f(t)dt = L-\int_a^b f(t)dt$ .
- 2. If f(t) is HK-integrable on [a,b], then  $\sup_{t>a}\left|\int_a^t f(s)ds\right|<\infty$ .
- 3. If f(t) is HK-integrable and bounded on [a, b], then it is L-integrable on [a, b].

#### 3.1 The KS-Hilbert Space

In this section, our objective is to construct a particular (separable) Hilbert space  $KS^2[\mathbb{R}^n]$ . This space is of special interest because it contains the class of HK-integrable functions, the space of the test functions  $\mathcal{D}[\mathbb{R}^n]$ , the space  $\mathfrak{M}[\mathbb{R}^n]$  of finitely additive measures on  $\mathbb{R}^n$  and  $L^p[\mathbb{R}^n]$  for  $1 \leq p \leq \infty$ . Furthermore, each of the latter class spaces is contained in  $KS^2[\mathbb{R}^n]$  as continuous dense and compact embeddings (e.g., weakly convergent sequences in each of the above spaces are strongly convergent in  $KS^2[\mathbb{R}^n]$ ). In particular,  $KS^2[\mathbb{R}^n]$  is perfect for the construction of Feynman's path integral formulation of quantum mechanics (in the form he suggested).

To construct  $KS^2[\mathbb{R}^n]$ , fix n, and let  $\mathbb{Q}^n$  be the set all vectors  $\{\mathbf{x}=(x_1,x_2\cdots,x_n)\in\mathbb{R}^n\}$  such that  $x_i$  is rational for each i. This is a countable dense set in  $\mathbb{R}^n$ , so we can arrange it as  $\mathbb{Q}^n=\{\mathbf{x}_1,\mathbf{x}_2,\mathbf{x}_3\cdots\}$ . For each l and i, let  $\mathbf{B}_l(\mathbf{x}_i)$  be the closed cube centered at  $\mathbf{x}_i$  of side  $r_l=2^{-l}, l\in\mathbb{N}$  parallel to the coordinate axis. Now choose an order so that the set  $\{\mathbf{B}_k(\mathbf{x}_k), k\in\mathbb{N}\}$  contains all closed cubes  $\{\mathbf{B}_l(\mathbf{x}_i) \mid (l,i)\in\mathbb{N}\times\mathbb{N}\}$  centered at a point in  $\mathbb{Q}^n$ . Let  $\mathcal{E}_k(\mathbf{x})=1$  if  $\mathbf{x}\in\mathbf{B}_k(\mathbf{x}_k)$  and be zero otherwise, so that  $\mathcal{E}_k(\mathbf{x})$  is in  $L^p[\mathbb{R}^n]$  for  $1\leqslant p\leqslant\infty$ . Define  $F_k(\cdot)$  on  $L^1[\mathbb{R}^n]$  by

$$F_k(f) = \int_{\mathbb{R}^n} \mathcal{E}_k(\mathbf{x}) f(\mathbf{x}) d\mathbf{x}.$$
 (3.1)

It is clear that  $F_k(\,\cdot\,)$  is a bounded linear functional on  $L^p[\mathbb{R}^n]$  for each k,  $\|F_k\|_{\infty} \leqslant 1$  and if  $F_k(f) = 0$  for all k, f = 0 so that  $\{F_k\}$  is fundamental on  $L^p[\mathbb{R}^n]$  for  $1 \leqslant p \leqslant \infty$ . Fix  $\lambda$ , set  $t_{\lambda}^k = \lambda^{k-1} e^{-\lambda} / (k-1)!$  and define a measure  $d\mathbf{P}_{\lambda}(\mathbf{x}, \mathbf{y})$  on  $\mathbf{R}^n \times \mathbf{R}^n$  by:

$$d\mathbf{P}_{\lambda}(\mathbf{x}, \mathbf{y}) = \left[ \sum_{k=1}^{\infty} t_{\lambda}^{k} \mathcal{E}_{k}(\mathbf{x}) \mathcal{E}_{k}(\mathbf{y}) \right] d\mathbf{x} d\mathbf{y}.$$

We can now define an inner product  $(\cdot)$  on  $\mathbf{L}^1[\mathbf{R}^n]$  by

$$(f,g) = \int_{\mathbb{R}^n \times \mathbb{R}^n} f(\mathbf{x}) g(\mathbf{y})^* d\mathbf{P}_{\lambda}(\mathbf{x}, \mathbf{y})$$

$$= \sum_{k=1}^{\infty} t_{\lambda}^k \left[ \int_{\mathbb{R}^n} \mathcal{E}_k(\mathbf{x}) f(\mathbf{x}) d\mathbf{x} \right] \left[ \int_{\mathbb{R}^n} \mathcal{E}_k(\mathbf{y}) g(\mathbf{y}) d\mathbf{y} \right]^*.$$
(3.2)

We call the completion of  $L^1[\mathbb{R}^n]$ , with the above inner product, the Kuelbs-Steadman space  $KS^2[\mathbb{R}^n]$ . To see that it contains the HK-integrable functions, let f be HK-integrable on  $\mathbb{R}^n$  then:

$$\|f\|_{KS}^2 = \sum_{k=1}^{\infty} t_{\lambda}^k \left| \int_{\mathbb{D}^n} \mathcal{E}_k(\mathbf{x}) f(\mathbf{x}) d\mathbf{x} \right|^2 \leqslant \sup_{k} \left| \int_{\mathbb{D}^n} \mathcal{E}_k(\mathbf{x}) f(\mathbf{x}) d\mathbf{x} \right|^2 \leqslant \left| \int_{\mathbb{D}^n} f(\mathbf{x}) d\mathbf{x} \right|^2 < \infty,$$

so  $f \in KS^2[\mathbb{R}^n]$ .

Journal of Physics: Conference Series 2987 (2025) 012003

doi:10.1088/1742-6596/2987/1/012003

**Definition 3.5.** We say that  $\mathcal{B} \subset \mathcal{H}$  is a continuous dense embedding if:

- 1.  $\mathcal{B}$  is dense in  $\mathcal{H}$ .
- 2. For each  $f \in \mathcal{B}$ ,  $||f||_{\mathcal{H}} < ||f||_{\mathcal{B}}$ .

**Theorem 3.6.** For each  $p, 1 \leq p \leq \infty$ ,  $L^p[\mathbb{R}^n] \subset KS^2[\mathbb{R}^n]$  is a continuous dense embedding.

*Proof.* By construction,  $KS^2[\mathbb{R}^n]$  contains  $L^1[\mathbb{R}^n]$  a continuous dense embedding, so we need only show that  $KS^2[\mathbb{R}^n] \supset L^p[\mathbb{R}^n]$  for  $p \neq 1$ . If  $p < \infty, L^p[\mathbb{R}^n] \cap L^1[\mathbb{R}^n]$  is dense in  $L^1[\mathbb{R}^n]$  and hence dense in  $KS^2[\mathbb{R}^n]$ . If  $f \in L^p[\mathbb{R}^n]$  then:

$$||f||_{KS^{2}} = \left[ \sum_{k=1}^{\infty} t_{\lambda}^{k} \left| \int_{\mathbb{R}^{n}} \mathcal{E}_{k}(\mathbf{x}) f(\mathbf{x}) d\mathbf{x} \right|^{\frac{2p}{p}} \right]^{1/2}$$

$$\leq \left[ \sum_{k=1}^{\infty} t_{\lambda}^{k} \left( \int_{\mathbb{R}^{n}} \mathcal{E}_{k}(\mathbf{x}) |f(\mathbf{x})|^{p} d\mathbf{x} \right)^{\frac{2}{p}} \right]^{1/2}$$

$$\leq \sup_{k} \left[ \int_{\mathbf{R}^{n}} \mathcal{E}_{k}(\mathbf{x}) |f(\mathbf{x})|^{p} d\mathbf{x} \right]^{\frac{1}{p}} \leq ||f||_{p}.$$

Hence,  $f \in KS^2[\mathbb{R}^n]$  and the embedding is continuous. For  $p = \infty$ , we have

$$||f||_{KS^2} = \left[\sum_{k=1}^{\infty} t_{\lambda}^k \left| \int_{\mathbb{R}^n} \mathcal{E}_k(\mathbf{x}) f(\mathbf{x}) d\mathbf{x} \right|^2 \right]^{1/2}$$

$$\leq \left[ \left[\sum_{k=1}^{\infty} t_{\lambda}^k [vol(\mathbf{B}_k)]^2 \right] [ess \sup |f|]^2 \right]^{1/2} \leq M ||f||_{\infty}.$$

Thus  $f \in KS^2[\mathbb{R}^n]$ , and  $L^{\infty}[\mathbb{R}^n] \subset KS^2[\mathbb{R}^n]$  as a continuous dense embedding.

We can also construct  $KS^p[\mathbb{R}^n]$  for all p by defining

$$||f||_{KS^p} = \begin{cases} \left\{ \sum_{k=1}^{\infty} t_k \left| \int_{\mathbb{R}^n} \mathcal{E}_k(\mathbf{x}) f(\mathbf{x}) d\mathbf{x} \right|^p \right\}^{1/p}, 1 \leqslant p < \infty, \\ \sup_{k \geqslant 1} \left| \int_{\mathbb{R}^n} \mathcal{E}_k(\mathbf{x}) f(\mathbf{x}) d\mathbf{x} \right|, p = \infty. \end{cases}$$

It is easy to see that  $\|\cdot\|_{KS^p}$  defines a norm on  $L^p[\mathbb{R}^n]$ . If  $KS^p[\mathbb{R}^n]$  is the completion of  $L^p[\mathbb{R}^n]$  with respect to this norm, we have:

**Theorem 3.7.** For each  $q, 1 \leq q \leq \infty$ ,  $KS^p[\mathbb{R}^n] \supset L^q[\mathbb{R}^n]$  as a continuous dense embedding and,  $KS^{\infty}[\mathbb{R}^n] \subset KS^p[\mathbb{R}^n]$  for each p.

Remark 3.8. We observe that  $KS^2[\mathbb{R}^n] = [KS^2]^*[\mathbb{R}^n] \supset L^1[\mathbb{R}^n]^{**} = \mathfrak{M}[\mathbb{R}^n]$ , the space of finitely additive measures on  $\mathbb{R}^n$ , and in particular, the Dirac measure  $\delta(\mathbf{x}) \in \mathfrak{M}[\mathbb{R}^n]$ . The theory of distributions was developed via topological vector space methods because of the belief that there was no Hilbert space which contained the test functions  $\mathcal{D}[\mathbb{R}^n]$ . The next result shows their mistake was due to a lack of interest in non-absolutely integrable functions.

**Definition 3.9.** The space of test functions  $\mathcal{D}[\mathbb{R}^n] = C_c^{\infty}[\mathbb{R}^n]$  (the set of continuous functions with the sequential compact support topology). A sequence  $\phi_j \in \mathcal{D}[\mathbb{R}^n]$  converges to  $\phi \in \mathcal{D}[\mathbb{R}^n]$  in the sequential compact support topology, if there exists a compact set  $K \subset \mathbb{R}^n$  containing the support of  $\phi_j - \phi$  for all j and, for every multi-index  $\alpha$ ,  $D^{\alpha}\phi_j$  converges to  $D^{\alpha}\phi$  uniformly on K.

**Theorem 3.10.** The test functions  $\mathcal{D}[\mathbb{R}^n] \subset KS^2[\mathbb{R}^n]$  as a continuous dense embedding.

*Proof.* Since  $KS^{\infty}(\mathbb{R}^n)$  is continuously embedded in  $KS^2[\mathbb{R}^n]$ , it suffices to prove the result for  $KS^{\infty}(\mathbb{R}^n)$ . Suppose that  $\phi_j \to \phi$  in  $\mathcal{D}[\mathbb{R}^n]$ , so that there exists a compact set  $K \subset \mathbb{R}^n$ , containing the support

doi:10.1088/1742-6596/2987/1/012003

of  $\phi_j - \phi$  and  $D^{\alpha}\phi_j$  converges to  $D^{\alpha}\phi$  uniformly on K for every multi-index  $\alpha$ . Let  $L = \{l \in \mathbb{N} : \text{the support of } \mathcal{E}_l, \ stp\{\mathcal{E}_l\} \cap K \neq \emptyset\}$  and  $M = \sup_{l \in L} [vol(\mathbf{B}_l)]$ , then

$$\lim_{j \to \infty} \|D^{\alpha} \phi - D^{\alpha} \phi_{j}\|_{KS^{\infty}} = \lim_{j \to \infty} \sup_{l \in L} \left| \int_{\mathbb{R}^{n}} \left[ D^{\alpha} \phi\left(x\right) - D^{\alpha} \phi_{j}\left(x\right) \right] \mathcal{E}_{l}\left(x\right) dx \right|$$

$$\leq M \lim_{j \to \infty} \sup_{x \in K} \left| D^{\alpha} \phi\left(x\right) - D^{\alpha} \phi_{j}\left(x\right) \right| = 0.$$

Thus, as the convergence is uniform, it follows that  $\mathcal{D}[\mathbb{R}^n] \subset KS^2[\mathbb{R}^n]$  as a continuous dense embedding. Since  $KS^2[\mathbb{R}^n]$  is self dual, we see that the Schwartz distributions,  $\mathcal{D}^*[\mathbb{R}^n] \subset KS^2[\mathbb{R}^n]$ .

If  $f, g \in L^1[\mathbb{R}^n]$ , we denote the Fourier transform of f and the convolution of g with respect to f by  $\mathfrak{F}(f)$  and  $\mathfrak{C}_f(g)$  respectively.

**Theorem 3.11.** The Fourier transform  $\mathfrak{F}(\cdot)$  and, for each  $f \in L^1[\mathbb{R}^n]$ , the convolution with respect to f,  $\mathfrak{C}_f(\cdot)$  both extend to  $KS^2[\mathbb{R}^n]$  as bounded linear operators.

**Remark 3.12.** The possibility that  $KS^2[\mathbb{R}^n]$  could exist was noticed by Gill and Zachary after reading the paper by Kuelbs [KB]. It was actually constructed by V. Steadman and first appeared as a part of her dissertation (see [ST]). As will be seen in the next section, this space is perfect for the Feynman path integral formulation of quantum mechanics.

### 3.2 Feynman Path Integral I

The properties of  $KS^2[\mathbb{R}^n]$  derived earlier shows that it is a better replacement for  $L^2[\mathbb{R}^n]$  for the study of the path integral formulation of quantum theory developed by Feynman. It is easy to prove that both the position and momentum operators have closed, densely defined extensions to  $KS^2[\mathbb{R}^n]$ . Furthermore, the extension of  $\mathfrak{F}(\cdot)$  and  $\mathfrak{C}_f(\cdot)$  by Theorem 3.11 insures that both the Schrödinger and Heisenberg theories have faithful representations on  $KS^2[\mathbb{R}^n]$ .

Since  $KS^2[\mathbb{R}^n]$  contains  $\mathfrak{M}[\mathbb{R}^n]$  the space of all finitely and countably additive measures as compact embeddings, it follows that all the approximating sequences for the Dirac measure converge strongly to it in the  $KS^2[\mathbb{R}^n]$  topology. Thus, the finitely additive set function defined on the Borel sets (Feynman kernel):

$$\mathbb{K}_{\mathbf{F}}[t, \mathbf{x}; s, B] = \int_{B} \left(2\pi i(t-s)\right)^{-1/2} \exp\{i|\mathbf{x} - \mathbf{y}|^{2} / 2(t-s)\} d\mathbf{y}$$

is in  $KS^2[\mathbb{R}^n]$  and  $\|\mathbb{K}_{\mathbf{F}}[t,\mathbf{x};s,B]\|_{KS^2} \leqslant 1$  and

$$\mathbb{K}_{\mathbf{F}}[t,\mathbf{x}\,;\,s,B] = \int_{\mathbb{R}^n} \mathbb{K}_{\mathbf{F}}[t,\mathbf{x}\,;\,\tau,d\mathbf{z}] \mathbb{K}_{\mathbf{F}}[\tau,\mathbf{z}\,;\,s,B], \ \ (\text{HK-integral}).$$

**Definition 3.13.** Let  $\mathbf{P}_n = \{t_0, \tau_1, t_1, \tau_2, \cdots, \tau_n, t_n\}$  be a HK-partition of the interval [0, t] for each n, with  $\lim_{n\to\infty} \Delta\mu_n = 0$  (mesh). Set  $\Delta t_j = t_j - t_{j-1}, \tau_0 = 0$  and for  $\psi \in KS^2[\mathbb{R}^n]$  define

$$\int_{\mathbb{R}^{n[0,t]}} \mathbb{K}_{\mathbf{F}}[\mathcal{D}_{\lambda} \mathbf{x}(\tau) ; \mathbf{x}(0)] = e^{-\lambda t} \sum_{k=0}^{[|\lambda t|]} \frac{(\lambda t)^k}{k!} \left\{ \prod_{j=1}^k \int_{\mathbb{R}^n} \mathbb{K}_{\mathbf{F}}[t_j, \mathbf{x}(\tau_j); t_{j-1}, d\mathbf{x}(\tau_{j-1})] \right\},$$

and

$$\int_{\mathbb{R}^{[0,t]}} \mathbb{K}_{\mathbf{F}}[\mathcal{D}\mathbf{x}(\tau); \mathbf{x}(0)] \psi[\mathbf{x}(0)]$$

$$= \lim_{\lambda \to \infty} \int_{\mathbb{R}^{n[0,t]}} \mathbb{K}_{\mathbf{F}}[\mathcal{D}_{\lambda}\mathbf{x}(\tau); \mathbf{x}(0)] \psi[\mathbf{x}(0)] \tag{3.3}$$

whenever the limit exists.

The next result is now elementary but a more general (sum over paths) result covering most areas of application will be given after the next section.

doi:10.1088/1742-6596/2987/1/012003

**Theorem 3.14.** The function  $\psi(\mathbf{x}) \equiv 1 \in KS^2[\mathbb{R}^n]$  and

$$\int_{\mathbb{R}^{3[s,t]}} \mathbb{K}_{\mathbf{F}}[\mathcal{D}\mathbf{x}(\tau) ; \mathbf{x}(s)] = \mathbb{K}_{\mathbf{F}}[t, \mathbf{x} ; s, \mathbf{y}] = \frac{1}{\sqrt{2\pi i(t-s)}} \exp\{i|\mathbf{x} - \mathbf{y}|^2 / 2(t-s)\}.$$

The above result is exactly what Feynman expected to obtain.

**Remark 3.15.** The above results also hold with no changes if we replace  $\mathbb{K}_{\mathbf{F}}[t, \mathbf{x}; s, B]$  by the kernel  $\mathbb{K}_{\mathbf{W}}[t, \mathbf{x}; s, B]$ , for the Wiener measure.

Since  $L^2[\mathbb{R}^n]$  is dense,  $KS^2[\mathbb{R}^n]$  is separable and every orthonormal basis for  $L^2[\mathbb{R}^n]$  may be normalized to make an orthonormal basis for  $KS^2[\mathbb{R}^n]$ . Thus they are unitarily equivalent but:

- 1.  $KS^2[\mathbb{R}^n]$  contains both  $\mathcal{D}[\mathbb{R}^n]$  and  $\mathcal{D}^*[\mathbb{R}^n]$ , while  $L^2[\mathbb{R}^n]$  does not.
- 2. The integral of  $e^{ix^2}$  is finite in  $KS^2[\mathbb{R}^n]$  but infinite in  $L^2[\mathbb{R}^n]$ .

Thus,  $L^2[\mathbb{R}^n]$  and  $KS^2[\mathbb{R}^n]$  are mathematically equivalent but not physically equivalent.

We show that interactions of a general nature become possible when the path integral is extended as part of Feynman's time ordered operator calculus (see [30]).

### Discussion and Conclusion

In this paper we have faced directly a number of mathematical problems that have clouded our ability to see how to correctly address some of the physical problems left over from the last century.

#### 3.3 Classical Level

At the classical level, we have made Minkowski's unacknowledged third (mathematical) postulate for Einstein's special theory of relativity explicit and provided a direct proof that it is incompatible with the first postulate for two or more particles. We then use one of three representations for the proper time to show how to resolve this incompatibility. This then allowed us to resolve the center of mass problem first noticed by Pryce [53] and the no Interaction problem first conjectured by Bakamjian and Thomas [7] and later proved by Currie, Jordan and Sudarshan [9]. This clears the way for a direct extension of Newton's theory to one compatible with that of Maxwell (addressed in the second paper).

## 3.4 Dirac and Square Root Operator

We have first provided an analytical separation (diagonalization) of the full (minimal coupling) Dirac equation into particle and antiparticle components. The diagonalization is analytic in that it was achieved without transforming the wave functions and reveals the nonlocal time behavior of the particle-antiparticle relationship.

It is well-known that the square root operator is nonlocal in space but and explicit representation and it physical properties has not been available. We have provided such a representation, (due to Samko [57]). We have discussed its singularity properties and shown that it looks like point particle outside a few Compton wave lengths. The fact that it is also directly related to the Dirac operator by a unitary transformation (Foldy-Wouthuysen [22]) shows that mathematical equivalence does not always imply physical equivalence.

## 3.5 Feynman Path Integral

Historically, the mathematics community has had two responses to the introduction of a new mathematical idea or method into physics. The first response has been to fit the idea or method into an existing framework. The second and more exciting response is when such an idea or method leads to the development of a new branch of mathematics. The most prevalent and successful response has been in finding an existing mathematical structure that will reasonably accommodate the physical theory and provide (at least) the framework for mathematical rigor. In some rare but important instances, there is no obvious mathematical structure which can completely accommodate the theory in the manner presented by physicists. In this case, mathematicians have extended and/or adapted an existing mathematical theory,

doi:10.1088/1742-6596/2987/1/012003

developed new mathematical structures or suggested (in frustration) that any conclusions derived from the use of these ideas or methods are at least suspect. Over the last seventy years, all of the above positions have appeared in response to Feynman?s introduction of his path integral into quantum theory.

In this paper, we have introduced a more general integral (HK) and designed a new Hilbert space  $(KS^2[\mathbb{R}^n])$  for the Feynman theory. This approach provides the correct framework for mathematical rigor without removing the physically intuitive and computationally effective approach suggested by Feynman. This space is also unitarily equivalent to  $L^2[\mathbb{R}^n]$  showing that unitary equivalence is not always indicative of the same physics.

#### Acknowledments

We would like to thank Professors Netsivi Ben-Amots, Alexander Gersten, Lawrence P. Horwitz, Martin C. Land, Elliot Leib and Ruggero M. Santilli for their continued interest and support.

This paper is dedicated in honor of Professor Lawrence P. Horwitz (Larry) on the occasion of his ninety-second birthday. The first author would like to register his deep admiration and affection for the many kind ways the Larry has freely given support, helpful suggestions and advice over the last 26 years. He has always conducted himself in a manner that represents the best model of what a scientific professional should be.

#### Declaration

The authors certify that:

- 1. they have no relevant financial or non-financial interests to disclose;
- 2. they have no conflicts of interest to declare that are relevant to the content of this manuscript;
- 3. they have no financial or proprietary interests in any material discussed in this manuscript; and
- 4. they have no affiliations with or involvement in any organization or entity with any financial interest or non-financial interest in the subject matter or materials discussed in this manuscript.

#### References

- [1] S. Albeverio and S. Mazzucchi Feynman path integrals for polynomially growing potentials, J. Funct. Anal. 221 (2005), 83-121.
- [2] A. Salam, Overlapping divergence and the S-matrix, Phys. Rev. 82, (1951) 217-227.
- [3] L. C. M. Benjamin, R. M. Antonio and A. P. Gonzalo, *The 2.7 °K blackbody radiation background reference frame*, Chin. Phys. B **19** No. 4 (2010) 04XXXX1-5.
- [4] E. T. Bell, Men of Mathematics, Simon & Schuster, New York, (1937).
- [5] B. Javanmardi; C. Porciani; P. Kroupa; J. Pflamm-Altenburg *Probing the Isotropy of Cosmic Acceleration Traced By Type Ia Supernovae* The Astrophysical Journal Letters. **810** (1) (2015) 47.
- [6] H.R Brown, Eur. J. Phys. 26, S85 (2005).
- [7] B. Bakamjian and L. H. Thomas, Phys. Rev. **92**, 1300 (1953).
- [8] C. D. Anderson, The Positive Electron, Physical Review 43, 491 (1943).
- [9] D. G. Currie, T. F. Jordan, and E. C. G. Sudarshan, Rev. Mod. Phys. 35, 350 (1963).
- [10] P. A. M. Dirac, Classical theory of radiating electrons, Proceedings of the Royal Soc. of London A 167 (1938) 148-168.
- [11] P. A. M. Dirac, Proc. Roy. Soc (London) A117 (1928) 610, A118 (1928) 351.
- [12] Einstein Centennial Symposium, October 1-5, 1979, University of Illinois at Carbondale, USA.

- [13] F. Dyson, The S-matrix in quantum electrodynamics, Phys. Rev. 75 (1949), 1736-1755.
- [14] A. Einstein, Ann. d. Phys. 17, 891 (1905).
- [15] A. Einstein, Jahrbuch Radioaktivitat V, 422 (1907) (Berichtigungen).
- [16] A. Einstein and J Laub, Uber die elektromagnetischen Grundgleichugen f¨ur bewegte K¨oper, Ann. d. ¨ Phys. 33(8), 532-540 (1908).
- [17] A. Einstein and J Laub, Uber die im elektromagnetischen Felde auf ruhende K¨oper ausge¨ubten ¨ ponderomotorischen K¨rafte, Ann. d. Phys. 331(8), 551-550 (1908).
- [18] J. R. Fanchi, Confronting the ENIGMA of TIME, World Scientific, Singapore, (2023).
- [19] R. P. Feynman, Phys. Rev. 81, 108 (1951).
- [20] R. P. Feynman and M. Gell-Mann, Phys. Rev. 109, 193 (1958).
- [21] M. Frisch, Inconsistency Asymmetry and Non-Locality, Oxford University Press N. Y. (2005).
- [22] L. L. Foldy and S. A. Wouthuysen, Phys. Rev. 78, 29 (1950).
- [23] J. Glimm and A. Jaffe, Quantum Physics. A Functional Integral Point of View, Springer, N. Y., (1987).
- [24] T. L. Gill and J. Lindesay, Canonical Proper Time Formulation of Relativistic Particle Dynamics, Int. J. of Theor. Phys. 32 (1993), 2087-2098.
- [25] T. L. Gill, The Square-Root Operator, Proper-Time and a Particle Interpretation for the Klein-Gordon Equation, FERMILAB-Pub-82/60-THY (1982).
- [26] T. L. Gill, and G. Ares de Parga https://iopscience.iop.org/article/10.1088/1742- 6596/2482/1/012015
- [27] T. L. Gill, and G. Ares de Parga https://iopscience.iop.org/article/10.1088/1742- 6596/1239/1/012013
- [28] I. S. Gradshteyn and I. M. Ryzhik Tables of Integrals, Series and Products, Academic Press, NewYork, (1980).
- [29] T. L. Gill and W. W. Zachary, Foundations for relativistic quantum theory I: Feynman's operator calculus and the Dyson conjectures, Journal of Mathematical Physics 43 (2002), 69-93.
- [30] T. L. Gill and W. W. Zachary, Functional Analysis and the Feynman operator Calculus, Springer New York, (2016).
- [31] T.L. Gill, W.W. Zachary, Analytic representation of the square-root operator, J. Phys. A: Math. Gen. 38 (2005) 1-18.
- [32] T.L. Gill, W.W. Zachary, and M. Alfred, Analytic representation of the Dirac equation, J. Phys. A: Math. Gen. 38 (2005) 1-22.
- [33] T.L. Gill, W.W. Zachary and J. Lindesay, The Classical Electron Problem, Found. Phys. 31, 1299 (2001).
- [34] https://www.spiedigitallibrary.org/conference-proceedings-of-spie/7421/1/Two-versions-of-Maxwells-equations-and-the-nature-of-light/10.1117/12.824191.
- [35] D. Hestenes, Space Time Algebra, Gordon and Breach Pub. New York, (1966)
- [36] L. P. Horwitz and C. Piron, Helv. Phys. Acta 46, 316 (1981).
- [37] E. Hille, and R. S. Phillips, Functional Analysis and Semigroups, American Mathematical Society Colloquium Publication No. 31 American Mathematical Society, Providence, RI, (1957).
- [HW] L. P. Horwitz, On the significance of a recent experiment demonstrating quantum interference in time, (to appear in Phys. Rev. Letters, see arXiv:quant-ph/0507044).

- [38] T. Ichinose, Essential selfadjointness of the Weyl quantized relativistic hamiltonian, Ann. Inst. Henri Poincar´e (Physique th´eorique) 51, (1989) 265-298.
- [KB] J. Kuelbs, Gaussian measures on a Banach space, Journal of Functional Analysis 5 (1970), 354–367.
- [39] H. A. Lorentz, Archives Neerlandaises des Sciences Exactes et Naturelles, 25, 353 (1892).
- [40] H. A. Lorentz, The Theory of Electrons B. G. Teubner, Leipzig, 1906; (reprinted by Dover, New York, 1952).
- [41] H. Leutwyler. A no-interaction theorem in classical relativistic hamiltonian particle mechanics, Nuovo Cim., 37, 556 (1965).
- [42] G. Longhi, L. Lusanna, and J. M. Pons, J. Math. Phys. 30, 1893 (1989).
- [43] H. Minkowski, Physikalische Zeitschrift 10,104 (1909).
- [44] Nathan J. Secrest; Sebastian von Hausegger; Mohamed Rameez; Roya Mohayaee; Subir Sarkar; Jacques Colin. A Test of the Cosmological Principle with Quasars, The Astrophysical Journal Letters, 908 (2) 2021): L51
- [45] A. Pais, Subtle is the lord: the science and life of Albert Einstein, Oxford University Press (1982).
- [46] P.J.E. Peebles, Principles of Physical Cosmology, Princeton University Press, London (1993).
- [47] W. Perret and G. B. Jeffery (translators, with additional notes by A. Sommerfeld), The Principle of Relativity by H. A. Lorentz, A. Einstein, H. Minkowski and H. Weyl, Dover, New York, 1952.
- [48] H. Poincar´e, Sur la dynamique de l' ´electron, Rendiconti del Circolo matematico Rendiconti del Circolo di Palermo, 21, 129-176 (1906).
- [49] W. K. H. Panofsky and M. Phillips, Classical Electricity and Magnetism, (Second Edition), Addision-Wesley, Reading, MA. (1962).
- [50] W. Pauli, Uber den Zusammenhang des Abschlusses der Elektronengruppen im Atom mit der Kom- ¨ plexstruktur der Spektren. Zeitschrift f¨ur Physik. 31 (1): 765-783 (1925).
- [51] W. Pauli, Handbuch der Physik 2nd ed. 24, (1933).
- [52] A. A. Penzias and R. W. Wilson, Ap. J. 142, 419 (1965).
- [53] M.H.L. Pryce, Proc. Roy. Soc. London A 195, 400 (1948).
- [54] G. Ares de Parga, T.L. Gill and W.W. Zachary, The Thomas program and the canonical proper-time theory J. of Comp. Methods in Sci. and Eng. 3 (2013) 117-134.
- [55] M. Radovan On the Nature of Time, arXiv:1509.01498
- [56] F. Rohrlich, Classical Charged Particles, Advanced Book Program, Addison-Wesley, New York (1990).
- [57] S. G. Samko, Hypersingular Integrals and Their Applications, Analytical Methods and Special Functions, Vol. 5 Taylor & Frances, Pub., (2002).
- [58] J. Schwinger, Found. Phys. 13, 2573 (1998).
- [59] Saadeh D, Feeney SM, Pontzen A, Peiris HV, McEwen, JD . How Isotropic is the Universe, Physical Review Letters. 117 (13) (2016) 131302.
- [60] J. C. Slater, Quantum Theory of Atomic Structure, Vol. II McGraw-Hill, New York (1960).
- [ST] V. Steadman, Theory of operators on Banach spaces, Ph.D thesis, Howard University, 1988.
- [STW] R. F. Streater and A. S. Wightman, PCT, Spin and Statistics and All That, Benjamin, New York, (1964).
- [61] S. Schweber, QED And The Men Who Made It , Princeton University Press, London (1994).

doi:10.1088/1742-6596/2987/1/012003

- [62] S. Walters, In: Gonner, H., Renn, J., Ritter, J. (eds.) The Expanding Worlds of General Relativity. Einstein Studies, vol. 7, pp. 45-86. Birkh¨auser, Boston (1999).
- [63] J.A. Wheeler and R.P. Feynman, Interaction with the absorber as the mechanism of radiation Rev. Mod. Phys. 17 (1949) 157-181.
- [64] A. O. Williams, Phys. Rev. 58 (1940), 723.
- [65] S. Weinberg, High energy behavior in quantum field theory, Phys. Rev. 118 (1960), 838-849.