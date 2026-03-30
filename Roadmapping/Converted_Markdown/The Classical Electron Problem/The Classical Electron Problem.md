## The Classical Electron Problem

Tepper L. Gill1,2,<sup>4</sup> , W. W. Zachary1,<sup>5</sup> and J.Lindesay<sup>3</sup>

> <sup>1</sup>Department of Electrical Engineering <sup>2</sup>Department of Mathematics <sup>3</sup>Department of Physics Howard University,Washington, DC 20059 E-mail tgill@howard.edu

> > <sup>4</sup>Department of Physics University of Michigan Ann Arbor, Mich. 48109

<sup>5</sup>Department of Mathematics and Statistics University of Maryland University College College Park, Maryland 20742 E-mail: wwzachary@earthlink.net

# Abstract

In this paper, we construct a parallel image of the conventional Maxwell theory by replacing the observer-time by the proper-time of the source. This formulation is mathematically, but not physically, equivalent to the conventional form. The change induces a new symmetry group which is distinct from, but closely related to the Lorentz group, and fixes the clock of the source for all observers. The new wave equation contains an additional term (dissipative), which arises instantaneously with acceleration. This shows that the origin of radiation reaction is not the action of a "charge" on itself but arises from inertial resistance to changes in motion. This dissipative term is equivalent to an effective mass so that classical radiation has both a massless and a massive part. Hence, at the local level the theory is one of particles and fields but there is no self-energy divergence (nor any of the other problems). We also show that, for any closed system of particles, there is a global inertial frame and unique (invariant) global proper-time (for each observer) from which to observe the system. This global clock is intrinsically related to the proper clocks of the individual particles and provides a unique definition of simultaneity for all events associated with the system. We suggest that this clock is the historical clock of Horwitz, Piron, and Fanchi. At this level, the theory is of the action-at-a-distance type and the absorption hypothesis of Wheeler and Feynman follows from global conservation of energy.

PACS classification codes: 03.30.+,03.50.De

Keywords: generalized Maxwell theory, special relativity, radiation reaction

# 1.0 Introduction

It was 1865 when James Clark Maxwell published his theory of electrodynamics. The slow but steady progress made by our understanding and use of mechanics and thermodynamics was given a major boost by Maxwell's theory made practical. For example, starting from 1866, a continuous communications link has existed between Europe and the US ( due in no small part to the efforts of Lord Kelvin). By 1883, Edison had a workable light bulb, while Bell invented the telephone in 1886. The radio waves predicted by Maxwell were discovered by Hertz in 1887, and electricity, producing new inventions weekly, was well on the way to providing what we now consider normal.

In the intervening 41 years between Maxwell and the introduction of the special theory of relativity in 1905, a scientific and technological revolution had taken firm roots. Indeed, it has been suggested by Feynman<sup>1</sup> that, " From the long view of the history of mankindseen from, say ten thousand years... there can be little doubt that the most significant event of the 19th century will be judged as Maxwell's discovery of the laws of electrodynamics."

When the founding fathers, Lorentz, Poincar´e, Einstein, and their contemporaries began to study the issues associated with the foundations of electrodynamics; they had a number of options open to them in addressing the fact that the Newtonian theory and the Maxwell theory were invariant under different transformation groups: (see Jackson<sup>2</sup> )

- 1. *Both theories are incorrect and a correct theory is yet to be found.*
- 2. *The "proper" Maxwell theory will be invariant under the Galilean group.*
- 3. *The "proper" Newtonian theory will be invariant under the Lorentz group.*
- 4. *The assumption of an ether for electromagnetic propagation is correct so that*

*Galilean relativity applies to mechanics while electromagnetism has a preferred reference frame.*

At the time, it was unthinkable that the Maxwell theory had any serious flaws. Lorentz3,<sup>4</sup> had recently shown that all of the macroscopic phenomena of electrodynamics and optics could be accounted for based on an analysis of the microscopic behavior of electrons and ions.

Einstein<sup>5</sup> rejected the fourth possibility and, as noted by Spencer and Shama<sup>6</sup> , was the " first scientist with the foresight to realize that a formal postulate on the velocity of light was necessary." He proposed that all physical theories should satisfy the (now well-known) postulates of special relativity:

- 1. *The physical laws of nature and the results of all experiments are independent of the particular inertial frame of the observer (in which the experiment is performed).*
- 2. *The speed of light in empty space is constant and is independent of the motion of the source or receiver.*

The first postulate abandons the notion of an absolute space , while the second abandons absolute time. In a later paper, Einstein<sup>7</sup> modified the second postulate to make it explicit that he always referred to observers in inertial frames:

2 ′ . *The speed of light in empty space is constant and independent of the motion of the source or receiver in any inertial frame.*

Einstein formulated his theory in the usual three-dimensional notation, making a

distinction between time and space. It was noted by Poincar´e 8 that the transformations of Lorentz could be treated as rotations if time is made an imaginary coordinate. Poincar´e had also introduced the metric now attributed to Minkowski<sup>9</sup> .

Although Poincar´e discovered the proper-time, it was Minkowski who recognized its importance in physical theory and showed that it is the only unique variable associated with the source and available to all observers. Motivated by philosophical concerns, he further proposed that space and time should not be treated separately, but should be unified in the now well-known fashion leading to Minkowski space. Given the tremendous impact of the then-recent work in geometry on science, it was natural for him to think along these lines. Once he accepted this approach, it was also natural to assume that the proper-time of the source be used to parameterize the motion, acting as the metric for the underlying geometrization of the special theory of relativity, thus implicitly requiring that another postulate be added:

3. *The correct implementation of the first two postulates requires that time be treated as a fourth coordinate, and the relationship between components so constrained to satisfy the natural invariance induced by the Lorentz group of electrodynamics, (Minkowski space).*

The four-geometry postulate was very popular at the time and was embraced by many; but other important physical thinkers, including Einstein, Lorentz, Poincar´e, and Ritz, regarded it as a mathematical abstraction lacking physical content and maintained that space and time have distinct physical properties. Although Einstein demurred, the feeling among many of the leading physicists at that time was that an alternative implementation should be possible which preserves some remnant of an absolute time variable (true time), while still allowing for the constancy of the speed of light. It was noted by Whittaker<sup>10</sup> that a few weeks before he died, Lorentz is reported to have maintained his belief in the existence of this "true time". Dresden<sup>11</sup> reports that "· · · He retained his beliefs in a Euclidean, Newtonian space time, and in absolute simultaneity · · ·."

## 1.1 Perspective

The general focus on, and excitement about, the four-geometry left little room for serious alternative investigations (separated from philosophical debates). This is unfortunate since the diversion is part of the reason that the physical foundations of classical electrodynamics did not receive the early intense investigation accorded mechanics. Possibly because of the apparent completeness of the special theory, interest in statistical mechanics, quantum theory, and the problem of accelerated motion (the general theory), Einstein was preoccupied with these other important areas. On the other hand, the physics community lost three important thinkers on the subject by 1912. Ritz died in 1908, Minkowski died (shortly after his paper appeared) in 1909, and Poincar´e died in 1912. The First World War began in 1914 and within four years decimated a whole generation. Furthermore, by 1913 interests had already shifted from electrodynamics to the new quantum theory. The longer this investigation into classical electrodynamics was delayed, the more Minkowski's approach became embedded in the culture of physics, permeating the foundations for all future theories. By the time problems in attempts to merge the special theory of relativity with quantum theory forced researchers to take a new look at the foundations of classical electrodynamics, the Minkowski approach to the implementation of the special theory was considered almost sacred.

We are now taking our first steps into the twenty first century, one hundred and fortyfive years later. Electromagnetism is now in the hands of the engineers, mathematicians, and philosophers, and much of it is not considered mainstream physics. For those who learned physics in the sixties and seventies, "electrodynamics seems as old as mechanics". The continued success of quantum mechanics and the "apparent" successes of quantum electrodynamics and the standard model has made the subject pass´e. Today, students study the subject as an introduction to the special theory, preparation for advanced quantum theory, and as a simple example of a gauge theory. From this perspective, there is no real reason to believe that the first possibility should be rejected out of hand (i.e., that both the Newtonian and Maxwell theories could in some way be incorrect). Such a possibility is even more likely in light of the fact that the problems facing the early workers are still with us in one form or another. Furthermore, additional problems have arisen from both theory and experiment.

# 1.2 Problems

## Newtonian Mechanics

Once it was accepted that the "proper" Newtonian theory should be invariant under the Lorentz group, work on this problem was generally ignored until after World War Two when everyone realized that the quantum theory did not solve the problems left open by the classical theory. In particular, it was first noticed that (at the classical level) Minkowski's approach only works (as expected) in the one-particle case. It was 1948 when Pryce<sup>12</sup> showed that the canonical center-of-mass is not the three-vector part of a four-vector.

This variable is required for any "natural" relativistic many-particle theory. Virtually all research since then has focused on attempts to avoid this problem while maintaining use of the proper-time of the observer as the fourth coordinate for Minkowski geometry.

In order to provide a simple approach to the problem encountered by Pryce, let us consider two inertial observers X and X' with the same orientation. Assume that the (proper) clocks of X and X' both begin when their origins coincide and X' is moving with uniform velocity  $\mathbf{v}$  as seen by X. Let two particles, each the source of an electromagnetic field, move with velocities  $\mathbf{w}_i$  (i = 1, 2), as seen by X, and  $\mathbf{w}'_i$  (i = 1, 2), as seen by X', so that:

$$\mathbf{x}_{i}' = \mathbf{x}_{i} - \gamma(\mathbf{v})\mathbf{v}t + (\gamma(\mathbf{v}) - 1)(\mathbf{x}_{i} \cdot \mathbf{v} / \|\mathbf{v}\|^{2})\mathbf{v}, \tag{1.1a}$$

$$\mathbf{x}_{i} = \mathbf{x}_{i}' + \gamma(\mathbf{v})\mathbf{v}t' + (\gamma(\mathbf{v}) - 1)\left(\mathbf{x}_{i}' \cdot \mathbf{v} / \|\mathbf{v}\|^{2}\right)\mathbf{v}, \tag{1.1b}$$

with  $\gamma(\mathbf{v}) = 1/\left[1 - (\mathbf{v}/c)^2\right]^{1/2}$ , represent the spacial Lorentz transformations between the corresponding observers. Thus, there is clearly no problem in requiring that the positions transform as expected. However, when we try to transform the clocks, we see the problem at once since we must have, for example,

$$t' = \gamma(\mathbf{v}) \left( t - \mathbf{x}_1 \cdot \mathbf{v}/c^2 \right), \qquad t' = \gamma(\mathbf{v}) \left( t - \mathbf{x}_2 \cdot \mathbf{v}/c^2 \right).$$
 (1.2a)

This is clearly impossible except under very special conditions on all other observers. Furthermore, if we write down the center-of-mass position  $\mathbf{X}$  and require that it transform as above, we add another (impossible) constraint on the clock of any other observer. Pryce's approach is more abstract (and complicated), but leads to the same result.

In his 1949 paper, Dirac<sup>13</sup> observed that we must choose a particular realization of the Poincaré algebra in order to identify the appropriate variables for theory formulation. He showed that there are three possible choices of distinct three-dimensional hypersurfaces that are invariant under subgroups of the Poincar´e group and intersect every particle worldline once; the instant form, the point form, and the front form. (It was later shown by Leutwyler and Stern<sup>14</sup> that there are five choices. However, the other two are not especially interesting.) The instant form is best known. It is based on normal time-evolution and uses spacelike hyperplanes in Minkowski space; the point form is based on mass hyperboloids; while the front form is based on null hyperplanes.

Following Dirac's work, Bakamjian and Thomas<sup>15</sup> showed that one can construct a quantizable many-particle theory that satisfies the first two postulates of Einstein. However, they suggested that when interaction is introduced, their approach would not permit both a global theory and provide an invariant particle world-line description (satisfy the third postulate). This conjecture was generalized and later proved by Currie et al<sup>16</sup> to the effect that the requirements of Hamiltonian formulation, (canonical) independent-particle variables, and relativistic covariance (i.e the canonical positions transform as geometric coordinates), are only compatible with noninteracting particles (The No-Interaction Theorem). There are many references on the subject, but the book by Sudarshan and Mukunda<sup>17</sup> gives a comprehensive review of the problems and attempts to solve them (up to 1974). All attempts have ended in failure for one or more reasons which usually include the inability to quantize.

The No-Interaction Theorem led many to suspend the requirement that canonical positions transform as geometric coordinates and to focus on the construction of the "correct many-particle representation for the Poincar´e algebra". However, a very important (but

not well-known) theorem was proved by Fong and Sucher<sup>18</sup> in 1964 for the quantum case, and by Peres<sup>19</sup> in 1971 for the classical case:

Theorem 1.0 *(Fong-Sucher-Peres) Suppose that no restriction is put on the transformation law of the canonical variables of a many-particle system. Then given any Hamiltonian* H*, total momentun* P*, and angular momentum* J *satisfying:*

$$\begin{split} dH/dt &= 0, \quad [H,P_m] = 0, \quad [H,J_m] = 0, \\ [P_m,P_n] &= 0, \quad [J_m,P_n] = \varepsilon_{mns}P_s, \quad [J_m,J_n] = \varepsilon_{mns}J_s, \end{split}$$

*it is always possible to find a boost generator* L *so that the full set of commutation relations of the Poincar´e algebra for the inhomogeneous Lorentz group will be satisfied.*

In order to underscore the importance of this theorem, Peres showed explicitly how to construct a "clearly" nonrelativistic Hamiltonian and appropriate boost generator (along with canonical center-of-mass, total momentum, and angular momentum). Thus, this theorem implies that a relativistic classical (or quantum) many-particle theory requires something else besides the commutation relations for the inhomogeneous Lorentz group. On the other hand, this is the only requirement imposed on us by Maxwell's equations! It follows that, contrary to common belief, our historical (intellectual) state of affairs is not dictated by the Maxwell theory. We conclude that the Minkowski postulate imposes an additional condition on the special theory (not required by Maxwell's equations), but we are still unable to correctly account for Newtonian mechanics (after almost a hundred years). Those willing to dismiss the issue as arcane should be aware that the same problem also exists for the general theory. *Thus, the major problem facing us in the twenty first century is to construct a quantizable classical theory which satisfies the first two postulates of Einstein in some reasonable form and includes Newtonian Mechanics*.

# Interpretation

There are interpretation problems with the Minkowski approach that are not well-known. First, it should be noted that the conventional use of the words coordinate time tends to obscure the fact that this is the proper-time of the observer. This makes physical interpretation complicated and strange because one is required to refer back to the propertime of the source (or the postulated clock of a co-moving observer) in order to acquire a complete interpretation and analysis of experiments. Thus, the "parameter" (used to define the four-geometry) must also be viewed as a physically real measurable quantity when the theory is used for experimental analysis. At the classical level this asymmetrical relationship may be vexing, but it is not contradictory. However, at the quantum level this same problem becomes more fundamental. At this level, the observer proper-time is a c-number that transforms to an operator under the Lorentz group, while the proper-time of the source is an operator that remains invariant (see Wigner<sup>20</sup>).

## Radiation Reaction and the Lorentz-Dirac Equation

The problems associated with the radiation of accelerated charged particles, and those of the Lorentz-Dirac equation are old and well-known. Two books that have contributed to a clearer understanding of these basic problems are those of Rohrlich<sup>21</sup> and Parrott<sup>22</sup> . Rohrlich provides a comprehensive study of the classical theory up to 1965, which includes a nice review of the history. (Those unaware of the continuing effort to solve the classical electron problem should also see Rohrlich<sup>23</sup>.) Parrott's book is both clear and insightful. (His chapter on the Lorentz-Dirac equation is unbiased, well done, and should be required reading for any serious student of the subject.) The classics, Panofsky and Phillips<sup>24</sup>, and

Jackson<sup>2</sup> are also important sources of insight and history. The elementary (but correct) account by Feynman<sup>1</sup> in volume II of his famous lecture series has done much to educate those with little or no concern with the foundations.

The radiation of accelerated charged particles is known to occur instantaneously with acceleration and its nature has been the object of much speculation (see Wheeler and Feynman25). The great success of Lorentz in using the Maxwell (field) theory, along with his aether, to show that all of the macroscopic electrodynamics and optics could be derived from a microscopic analysis has done much to foster our faith in the correctness of the theory. This success carried with it our first introduction to the divergences of a field theory. He found that the energy density and the field momentum for each particle diverges unless the particle has a finite radius. In addition, the derived (Lorentz) force law did not provide the appropriate dissipation to account for the observed radiation. It was also known that the electromagnetic mass defined by the electrostatic energy divided by c <sup>2</sup> and that defined via the electomagnetic momentum did not agree, giving the well known 4/3's problem (see Schwinger<sup>26</sup>).

These problems led to the study of various finite-size models for charged particles and, in turn, forced serious consideration of the action of one part of a charge on itself (self-energy) and also required the introduction of extra forces to hold the particle together (Poincar´e stresses).

The appearance of the classical divergence difficulties in the quantized theory (along with a few new ones) led many to hope that the successful construction of a consistent classical theory would help to solve the corresponding problems in quantum electrodynamics.

For this reason, many attempts were made to formulate such a theory. The most wellknown early attempts are due to Born and Infield27, Dirac28, Bopp29, and Wheeler and Feynman25. (Less well-known other attempts are due to Rosen30, Podolsky and Schwed<sup>31</sup> , and Feynman32.) Each ran into problems with quantization and are a part of the history. However, the point particle reduction theory of Dirac and the Wheeler-Feynman approach have special importance.

The use of particles of finite radius causes serious problems with Lorentz invariance, so a major advance was made when Dirac constructed a point particle reduction theory for the Lorentz model. To do this, he used Maxwell's equations to find the retarded field of the particle, assuming that at large distances the field only contains outgoing waves, and then calculated the advanced field assuming that at large distances the field only contains converging waves. He then defined half the difference between the retarded and the advanced fields evaluated at the particle position, multiplied by the charge, as the force of radiation reaction. This term was added to the Lorentz force to provide the appropriate dissipation term (the Lorentz-Dirac equation). This provided the same dissipation term obtained by Lorentz (in a nonrelativistic calculation), but was independent of the particle radius. Thus, Dirac produced a point particle theory while all the other problems remained unchanged, and this is essentially what we have today. It should also be noted that point particles of finite mass imply infinite density. This was a real problem during Newton's time, but does not appear to cause problems today.

Wheeler and Feynman took a different ploy. They showed that we could use point particles, obtain the same radiation reaction term as above, and eliminate the self-energy divergence. Their approach assumes that the field which acts on a given particle arises only from other particles (adjunct field). They used half the sum of the retarded and the advanced fields, and assumed that there are sufficiently many particles in the system to completely absorb all radiation given off from any one of them (absorption hypothesis). The theory is of the action-at-a-distance type and also eliminates the divergences associated with the energy and momentum densities. Unfortunately, the theory could not be quantized, but this work made it clear that the action-at-a-distance and field theory approaches are much closer than was generally expected. (Indeed, Wheeler and Feynman argued that the two theories are complimentary views.)

The two best-known problems with the Lorentz-Dirac equation are runaway solutions and preacceleration. The equation has solutions for a free particle (with no force) that can self-accelerate off to infinity. It was conjectured that these solutions were eliminated by the asymptotic condition proposed by Haag<sup>33</sup>. However, Parrott<sup>22</sup> (pg. 196) notes that the asymptotic condition is necessary to ensure conservation of energy-momentum, but may not be sufficient to eliminate all strange solutions. Furthermore, the recent paper of Parrott and Endres<sup>34</sup> makes this conjecture doubtful. It has been recently shown by Low<sup>35</sup> that this problem also shows up at the nonrelativistic quantum level. Things are better for quantum electrodynamics (they don't appear), but caution is required as the possible existence of a Landau-like anomalous pole in the photon propagator or the electron-massive photon forward scattering amplitude could produce the runaway effect.

The preacceleration problem arises because the equation is nonlocal in time. This means that the particle can accelerate prior to the action of a force. The problem is generally ignored with the observation that the natural time interval for this effect (say for an electron) is of the order of 6.<sup>2</sup> <sup>×</sup> <sup>10</sup>−<sup>24</sup> sec, so that no classical particle can enter from a free state into interaction over such a small time interval.

These problems have existed for sometime now and no solution seems to be in sight. It is clear that the first problem is based on the assumption that the dissipation should be in the Lorentz force and, since this term is of third order in the position variable, the difficulty follows. The second problem can be traced back to the use of advanced fields which are necessary for the theory (Dirac and Wheeler and Feynman), and to get the correct dissipation term.

# Mach's Principle and the 2.7 ◦K MBR

Today, we know that a unique preferred frame of rest exists throughout the universe and is available to all observers. This is the 2.7 ◦K microwave background radiation (MBR) which was discovered by Penzias and Wilson<sup>36</sup> in 1965 using basic microwave equipment (by today's standards). This radiation is now known to be highly isotropic with anisotropy limits set at 0.001%. Futhermore, direct measurements have been made of the velocity of both our Solar System and Galaxy through this radiation (370 and 600 km/sec respectively, see Peebles<sup>37</sup> ). One can only speculate as to what impact this information would have had on the thinking of Einstein, Lorentz, Minkowski, Poincar´e, Ritz and the many other investigators of the early 1900's who were concerned with the foundations of electrodynamics and mechanics. The importance of this discovery for the foundations of electrodynamics in our view is that this frame is caused by *radiation from accelerated charged particles* (independent of the various cosmological suggestions).

As noted by Peebles, the MBR does not violate the special theory. However, general relativity predicts that at each point we can adjust our acceleration locally to find a freely falling frame where the special theory holds. In this frame, all observers with constant velocity are equivalent. Thus, according to the general theory we have an infinite family of freely falling frames. Within this context, the Penzias and Wilson findings show that there is a *unique frame in which both the acceleration and velocity can be set equal to zero at each point in the universe.*

As suggested by Rohrlich21, "Mach's principle was originally designed to ensure that there is no difference between the rotation of the earth with repect to the fixed stars or the fixed stars with respect to the earth." It now appears that the fixed stars are not needed and the earth really does rotate. Our concern with this principle is associated with the fact that an accelerated charged particle experiences a damping force simultaneously with the moment of acceleration (relative to any inertial frame). Thus, it appears that a charged particle can be used to identify accelerating frames and raises the question: what is a charged particle accelerating with respect to? Put another way, charged particles appear to know when they experience a force. Furthermore, even if the force is constant, the effect cannot be transformed away. This is a problem for any theory that seeks to unify electromagnetism with gravity.

# 1.3 Purpose

Dirac<sup>41</sup> was critical of the use of Minkowski geometry as fundamental. As late as 1963, he noted that "...the picture with four-dimensional symmetry does not give us the whole situation... Quantum theory has taught us that we must take a three-dimensional section of what appears to our consciousness at one time (an observation), and relate it to another three-dimensional section at another time." In reviewing attempts to merge gravitation with quantum theory, Dirac goes on to question the fundamental nature of the fourdimensional requirement in physics and notes that, in some cases, physical descriptions are simplified when one departs from it. The real question is: What do we replace it with that solves the outstanding problems and has some contact with the physics we know?

A major part of our strong belief in the fundamental nature of the covariant Minkowski approach to theory construction is based on the Feynman-Schwinger-Tomonaga formulation of QED and their great computational success in accounting for the Lamb shift and the anomalous magnetic moment. The correct history is at variance with this belief (see Schweber<sup>38</sup>). It should first be noted that, using noncovariant methods, French and Weisskopf<sup>39</sup>, and Kroll and Lamb<sup>40</sup> were the first to get the correct results. The history of the French and Weisskopf paper can be found in Schweber and is well worth reading. Both Schwinger and Feynman initially got incorrect results using their covariant formulation and only after the work of French and Weisskopf was circulated did they find their mistakes. Later, Tomonaga got the correct results but used noncovariant methods in the middle of the calculation (see Schweber<sup>38</sup>, pg. 270).

In attempting to solve the problems of the classical electron, almost every possible change has been explored except the Minkowski four-geometry requirement. Our purpose in this paper is to carefully study the mathematical and physical implications that arise when we replace the observer proper-time by the source proper-time in Maxwell's equations. In order to see how this is possible, we first recall Minkowski's definition of the proper-time

of a source:

$$d\tau^2 = dt^2 - \frac{1}{c^2} d\mathbf{x}^2 = dt^2 \left[ 1 - \left( \frac{\mathbf{w}}{c} \right)^2 \right], \ \mathbf{w} = \frac{d\mathbf{x}}{dt}, \tag{1.3a}$$

$$d\tau^2 = dt'^2 - \frac{1}{c^2} d\mathbf{x}'^2 = dt^2 \left[ 1 - \left( \frac{\mathbf{w}'}{c} \right)^2 \right], \ \mathbf{w}' = \frac{d\mathbf{x}'}{dt'}. \tag{1.3b}$$

Minkowski was aware that  $d\tau$  is not an exact one-form and this observation may have affected his decision to restrict its use to being a parameter for the four-geometry. However, there is an important physical reason why it is not an exact (mathematical) one-form. Physically, a particle can traverse many different paths (in space) during any given  $\tau$  interval. This reflects the fact that the distance a particle can travel in a given time interval depends on the forces acting on it. This implies that the clock of the source carries physical information, and there is no a priori physical reason to believe that this information is properly encoded when  $\tau$  is used as a parameter. We rewrite (1.3) as

$$dt^{2} = d\tau^{2} + \frac{1}{c^{2}}d\mathbf{x}^{2} = (d\tau)^{2} \left[ 1 + \left(\frac{\mathbf{u}}{c}\right)^{2} \right], \ \mathbf{u} = \frac{d\mathbf{x}}{d\tau}, \tag{1.4a}$$

$$dt'^{2} = d\tau^{2} + \frac{1}{c^{2}}d\mathbf{x}'^{2} = (d\tau)^{2}\left[1 + \left(\frac{\mathbf{u}'}{c}\right)^{2}\right], \ \mathbf{u}' = \frac{d\mathbf{x}'}{d\tau}.$$
 (1.4b)

Thus, another possibility appears (which does give an exact one-form). In case we have two or more particles, our new time transformations are replaced by (in the simplest case)

$$a_i'\tau_i = \gamma(\mathbf{v})[a_i\tau_i - \mathbf{x}_i \cdot \mathbf{v}/c^2],$$
 (1.2b)

where  $\tau_i$  is the proper-time of the i-th particle and  $a_i$  and  $a'_i$  are terms which depend only on  $\tau_i$ .

In Section 2 we construct the invariance group which fixes the proper-time of the source in the single particle case and then explore some of the physical implications and

interpretations of this approach. At this level we see that the speed of particles may be faster than the speed of light. The physical interpretation is that the mass and the mean lifetime of unstable particles are now both constant, while the velocity computed using the clock of the source replaces the velocity computed using the observer clock. Thus, as will be seen, there is no contradiction with the second postulate, only a change in conventions. The second postulate is shown to always hold for experiments conducted with the source at rest in the frame of the observer, as is the case for the Michelson-Morley experiment.

In Section 3 we show explicitly that Maxwell's equations have an equivalent representation which fixes the proper-time of the source for all observers. We then prove that this formulation is left covariant under the action of the proper-time group. Although the fields have the same transformation properties as the conventional formulation, both the current and charge densities transform differently. In particular, we prove that if the charge density is at rest in any inertial frame then it is invariant (not just covariant) for all observers. By example, even in the accelerating case when the proper velocity is 2c, the relative velocity of our observers must be a subtantial fraction of c for them to detect any difference in their measured properties of the charge distribution.

In this Section we also derive the corresponding wave equations and show that they contain an additional dissipative term,which arises instantaneously with acceleration. By a change of variables, we show that the dissipative term is equivalent to an effective mass for electromagnetic radiation. We validate this interpretation by directly calculating the energy radiated by an accelerated charge in the proper-time formulation. The radiation formulas obtained are close in form but differ from those computed via the conventional

formulation, but agree in the low-velocity limit. In particular, the proper-time theory predicts an additional term for the E-field which acts along the direction of motion (longitudinal), proving the validity of our interpretation of the wave equation. This result means that in the proper-time formulation, there is no need to to require that the charge act back on itself in order to account for radiation reaction. When we couple this result with the the invariance of the charge density, we are able to prove that the proper-time theory is independent of the particle size, structure and geometry.

In Section 4 we derive the related versions of the optical Doppler effect and the aberration of wave vectors. These two phenomena are both well-known and ubiquitous. However, the general forms are usually derived using Lorentz transformations<sup>2</sup>,<sup>42</sup>. Here, we derive them from the proper-time theory, using the new invariance group. In addition to the usual terms, we obtain new results because of the nonlocal frequency effects implied by our theory. These effects play an important role in our derivation of the group velocity for electromagnetic waves. Here we show that the group velocity is c only when measured in the (rest) frame of the observer, but will not be c for any other observer moving relative to that frame. The new value (in the simplest case) will be either c+v or c−v, depending on the direction of the relative motion. However, as will be shown in Section 6, this effect is in the noise for experiments conducted up to now because of theory interpretation.

In Section 5 we formulate a global interacting many-particle theory. With an eye towards the quantum theory, we require that the change from observer proper-time to source proper-time be canonical. This leads to the Hamiltonian which generates τ translations. To accomplish this, we use a representation of the proper-time that is independent of the number of particles. We derive our many-particle theory via the commutation relations for the Poincar´e algebra. As a side benefit, we show that the global system has a (unique) proper-time (avaliable for all observers). This clock provides a unique definition of simultaneity for all events associated with the system and is (shown to be) intrinsically related to the proper-times of the particles (in the system). From these results, it follows that at the local level, during interaction, the proper-time group is a nonlinear and nonlocal representation of the Lorentz group. On the other hand, at the global level, the proper-time group differs from the Lorentz group by a scale transformation. It follows from the work in this Section and in Section 2 that the group representation space is Euclidean.

In Section 6 we explore the ramifications and implications of our formulation and discuss some apparent disadvantages.

## 2.0 Proper-Time Transformations

In this section, we derive the transformations that fix the proper-time of the source for all observers. If we set b <sup>2</sup> = u <sup>2</sup> + c 2 , then from (1.1) and (1.4) we have that

$$t = (1/c) \int_{0}^{\tau} b(s)ds$$
 and  $t' = (1/c) \int_{0}^{\tau} b(s)'ds$ . (2.0)

It follows that t and t ′ are nonlocal as functions of τ in the sense that their values depend on the particular physical history (proper-time path) of the source. By the mean value property for integrals, we can find a unique s(τ ) for each τ , 0 < s(τ ) < τ , such that u<sup>τ</sup> = u(τ − s(τ )), and

$$t = (1/c) \int_{0}^{\tau} b(s)ds = (\bar{b}_{\tau}/c)\tau,$$
 (2.1a)

$$t' = (1/c) \int_{0}^{\tau} b'(s)ds = (\bar{b}'_{\tau}/c)\tau.$$
 (2.1b)

It is clear that this property is observer-independent since

$$t' = \gamma(\mathbf{v})(t - \mathbf{x} \cdot \mathbf{v}/c^2) \Rightarrow (\bar{b}_{\tau}'/c)\tau = \gamma(\mathbf{v})[(\bar{b}_{\tau}/c)\tau - (\mathbf{x} \cdot \mathbf{v}/c^2)]. \tag{2.2}$$

With a fixed clock for all observers, we can now develop a theory in which only the spatial coordinates are transformed. Using (2.2), the required transformations are

$$\mathbf{x}' = \mathbf{x} - \gamma(\mathbf{v})(\bar{b}_{\tau}/c)\mathbf{v}\tau + (\gamma(\mathbf{v}) - 1)(\mathbf{x} \cdot \mathbf{v}/||\mathbf{v}||^2)\mathbf{v}, \tag{2.3a}$$

$$\mathbf{x} = \mathbf{x}' + \gamma(\mathbf{v})(\bar{b}_{\tau}'/c)\mathbf{v}\tau + (\gamma(\mathbf{v}) - 1)(\mathbf{x}' \cdot \mathbf{v}/||\mathbf{v}||^2)\mathbf{v}.$$
(2.3b)

From a physical point of view, (2.3) tells us (explicitly) that observers can only share information about the past position of a given physical system. The above approach also gives us the only (presently known) rational solution to the problem of distant simultaneity. It is clear that all observers have the option of using their proper clocks with no hope of agreeing on the time occurrence of any event associated with the source. On the other hand, if each observer agrees to use the proper clock of the source, we see that they will always agree on the time occurrence of any event associated with the source.

We now see that  $a_i = (\bar{b}_i/c)$  and  $a'_i = (\bar{b}'_i/c)$  in equation (1.2b). The unit for b and b' is velocity so that physical interpretation is very important. It will arise naturally when we represent Maxwell's equations using the proper-time of the source. For now, we note that they are related by

$$b' = \gamma(\mathbf{v}) \left[ b - \frac{\mathbf{u} \cdot \mathbf{v}}{c} \right], \qquad b = \gamma(\mathbf{v}) \left[ b' + \frac{\mathbf{u}' \cdot \mathbf{v}}{c} \right].$$
 (2.4)

For any vector  $\mathbf{d}$ , set

$$\mathbf{d}^{\dagger} = \mathbf{d}/\gamma(\mathbf{v}) - (1 - \gamma(\mathbf{v})) \left[ \mathbf{v} \cdot \mathbf{d}/(\gamma(\mathbf{v})\mathbf{v}^2) \right] \mathbf{v}. \tag{2.5}$$

Then the full set of transformations between observers that fix the proper-time of the source take the (almost) familiar form

$$\mathbf{x}' = \gamma(\mathbf{v}) \left[ \mathbf{x}^{\dagger} - (\mathbf{v}/c)\bar{b}_{\tau}\tau \right], \qquad \mathbf{x} = \gamma(\mathbf{v}) \left[ \mathbf{x}'^{\dagger} + (\mathbf{v}/c)\bar{b}'_{\tau}\tau \right],$$
 (2.6)

$$\mathbf{u}' = \gamma(\mathbf{v}) \left[ \mathbf{u}^{\dagger} - (\mathbf{v}/c)b \right], \qquad \mathbf{u} = \gamma(\mathbf{v}) \left[ \mathbf{u}'^{\dagger} + (\mathbf{v}/c)b' \right],$$
 (2.7)

$$\mathbf{a}' = \gamma(\mathbf{v}) \left\{ \mathbf{a}^{\dagger} - \mathbf{v} \left[ \mathbf{u} \cdot \mathbf{a}/(bc) \right] \right\}, \quad \mathbf{a} = \gamma(\mathbf{v}) \left\{ \mathbf{a}'^{\dagger} + \mathbf{v} \left[ \mathbf{u}' \cdot \mathbf{a}'/(b'c) \right] \right\},$$
 (2.8)

where  $\mathbf{a}$  ( $\mathbf{a}'$ ) is the particle proper-(three) acceleration. The above transformations (along with (2.4)) form the proper-time group. In this formulation, we now have only one clock as an intrinsic part of the theory.

The above transformations are so close to Lorentz transformations that one might wonder if any new physics is possible. Not only is there new physics, as will be seen later, but just as importantly, there are new physical interpretations of old ideas. For example, relativistic momentum increase is attributed to relativistic mass increase so that

$$\mathbf{p} = m\mathbf{w}, \quad m = m_0[1 - w^2/c^2]^{-1/2}.$$
 (2.9a)

In the new interpretation,

$$\mathbf{p} = m_0 \mathbf{u}, \quad \mathbf{u} = \mathbf{w} [1 - w^2/c^2]^{-1/2},$$
 (2.9b)

so there is no mass increase, the (proper) velocity increases. Thus, in particle experiments the particle will have a fixed mass and decay constant, independent of its velocity. On the other hand, the particle can have (proper) speeds larger than the speed of light since its velocity is now interpreted to be dx/dτ . All cases where time dilation is discussed in the standard approach are replaced by statements about u in the new approach.

Note that the relationship between u and w can be viewed as dual in the sense that

$$\mathbf{u} = \mathbf{w}[1 - w^2/c^2]^{-1/2},\tag{2.10}$$

$$\mathbf{w} = \mathbf{u}[1 + u^2/c^2]^{-1/2}. (2.11a)$$

This relationship was first derived by Schott<sup>43</sup> in the famous 1915 paper in which he also derived the well-known Schott term of classical electrodynamics. Dividing by c in (2.11a), we get

$$\frac{\mathbf{w}}{c} = \frac{\mathbf{u}}{b}.\tag{2.11b}$$

It is easy to show that [1 + u 2 c 2 ] <sup>1</sup>/<sup>2</sup> = [1 <sup>−</sup> <sup>w</sup> 2 c 2 ] −1/2 . Expanding both sides and using (2.11b), we have

$$[1 + u^2/c^2]^{1/2} = 1 + \frac{1}{2} \frac{u^2}{c^2} - \frac{1}{8} \frac{u^4}{c^4} + \cdots,$$
 (2.12)

$$[1 - w^2/c^2]^{-1/2} = 1 + \frac{1}{2} \frac{w^2}{c^2} + \frac{3}{8} \frac{w^4}{c^4} + \cdots,$$
 (2.13a)

$$[1 - w^2/c^2]^{-1/2} = [1 - u^2/b^2]^{-1/2} = 1 + \frac{1}{2} \frac{u^2}{b^2} + \frac{3}{8} \frac{u^4}{b^4} + \cdots$$
 (2.13b)

Thus, all three expressions agree in the low-velocity region. It follows that all the results derived from the standard implementation of special relativity using w/c can also be consistently derived using u/b. This result will be repeatedly exploited in this paper to provide an alternative interpretation of much of classical electrodynamics. The real question that arises is which of these definitions of velocity is appropriate in the construction of faithful representations of physical reality. (see Section 6).

# 3.0 Proper-Time Maxwell Equations

In order to formulate the corresponding Maxwell theory, we need the following theorem which is derived from (2.0) and (2.6):

Theorem 3.1 *The transformation properties of the derivatives when the observers use the clock of the source are*:

$$\frac{1}{c}\frac{\partial}{\partial t} = \frac{1}{b}\frac{\partial}{\partial \tau}, \ \frac{1}{c}\frac{\partial}{\partial t'} = \frac{1}{b'}\frac{\partial}{\partial \tau}, \tag{3.1}$$

$$\nabla = \gamma(\mathbf{v}) \left[ \nabla' - (\mathbf{v}/cb')(\partial/\partial\tau) \right]. \quad \nabla' = \gamma(\mathbf{v}) \left[ \nabla + (\mathbf{v}/cb)(\partial/\partial\tau) \right], \tag{3.2}$$

Proof: For each case, we prove the first result. For the first case, we use the chain rule so that (1/c)∂/∂t = (1/c)(∂τ /∂t)(∂/∂τ ). Using equation (1.3a) and the fact that [1 <sup>−</sup> <sup>w</sup><sup>2</sup>/c<sup>2</sup> ] <sup>1</sup>/<sup>2</sup> = [1 + u <sup>2</sup>/c<sup>2</sup> ] −1/2 , we have

$$(1/c)(\partial \tau/\partial t) = (1/c)[1 - \mathbf{w}^2/c^2]^{1/2} = (1/c)[1 + \mathbf{u}^2/c^2]^{-1/2} = (1/b).$$
(3.3a)

This gives the first part of (3.1). To prove the first part of (3.2), we use equation (1.1a) to get that (with an obvious abuse of notation)

$$\frac{\partial}{\partial \mathbf{x}} = \frac{\partial \mathbf{x}'}{\partial \mathbf{x}} \frac{\partial}{\partial \mathbf{x}'} + \frac{\partial t'}{\partial \mathbf{x}} \frac{\partial \tau}{\partial t'} \frac{\partial}{\partial \tau}.$$
(3.3b)

Now note that (∂x ′/∂x) = <sup>γ</sup>(v), (∂t′/∂x) = <sup>−</sup>γ(v)v/c<sup>2</sup> , and (∂τ /∂t′ ) = (c/b′ ). Putting these terms in equation (3.3b) gives our result.

We can now formulate the proper-time version of Maxwell's equations. The conventional form of these equations for two observers is (in Gaussian units):

$$\nabla \cdot \mathbf{B} = 0, \qquad \nabla \times \mathbf{E} + \frac{1}{c} \frac{\partial \mathbf{B}}{\partial t} = 0,$$
 (3.4a)

$$\nabla \cdot \mathbf{E} = 4\pi \rho, \quad \nabla \times \mathbf{B} = \frac{1}{c} \left[ \frac{\partial \mathbf{E}}{\partial t} + 4\pi \rho \mathbf{w} \right],$$
 (3.4b)

$$\nabla' \cdot \mathbf{B}' = 0, \qquad \nabla' \times \mathbf{E}' + \frac{1}{c} \frac{\partial \mathbf{B}'}{\partial t'} = 0,$$
 (3.5a)

$$\nabla' \cdot \mathbf{E}' = 4\pi \rho', \quad \nabla' \times \mathbf{B}' = \frac{1}{c} \left[ \frac{\partial \mathbf{E}'}{\partial t'} + 4\pi \rho' \mathbf{w}' \right].$$
 (3.5b)

Using (2.11) and (3.1) - (3.2), the above equations can be rewritten using the proper-time of the source to get

$$\nabla \cdot \mathbf{B} = 0, \qquad \nabla \times \mathbf{E} + \frac{1}{b} \frac{\partial \mathbf{B}}{\partial \tau} = 0,$$
 (3.6a)

$$\nabla \cdot \mathbf{E} = 4\pi \rho, \quad \nabla \times \mathbf{B} = \frac{1}{b} \left[ \frac{\partial \mathbf{E}}{\partial \tau} + 4\pi \rho \mathbf{u} \right],$$
 (3.6b)

$$\nabla' \cdot \mathbf{B}' = 0, \qquad \nabla' \times \mathbf{E}' + \frac{1}{b'} \frac{\partial \mathbf{B}'}{\partial \tau} = 0,$$
 (3.7a)

$$\nabla' \cdot \mathbf{E}' = 4\pi \rho', \quad \nabla' \times \mathbf{B}' = \frac{1}{b'} \left[ \frac{\partial \mathbf{E}'}{\partial \tau} + 4\pi \rho' \mathbf{u}' \right].$$
 (3.7b)

We see that when observers use the proper-time of the source, the velocity of electromagnetic waves depends on the motion (of the source), and has magnitude larger than c.

This may seem strange and even contradictory to the second postulate: "The speed of
light in any inertial frame is constant and is independent of the motion of the source or
receiver." This is not the case. On closer inspection, it is clear that the second postulate
assumes that the observer's proper-clock is being used to measure time. Thus, there is no
contradiction, just a change in conventions.

In the Michelson-Morley experiment, the source is at rest in the frame of the observer so that  $\mathbf{u} = \mathbf{0}$  and b = c. It follows that this approach (also) explains the Michelson-Morley null result. It also provides agreement with the conceptual (but not technical) framework proposed by Ritz<sup>44</sup>; namely, that the speed of light does depend on the (proper) motion of the source. In this sense, both Einstein and Ritz were correct.

We could follow Einstein's method<sup>5</sup> in proving the covariance of the proper-time equations (using (3.1)-(3.2)). However, we use the four-vector approach first, to emphasize the fact that our theory is compatible with four-vectors (in the one-particle case) and second, because it will be convenient for our derivation of the proper-time transformation of plane waves in Section 4. (The plane waves will be used to derive formulas for the Doppler effect and aberration of wave vectors.) Writing our equations in four-dimensional form as

$$F = \begin{bmatrix} 0 & B_z & -B_y & -iE_x \\ -B_z & 0 & B_x & -iE_y \\ B_y & -B_x & 0 & -iE_z \\ iE_x & iE_y & iE_z & 0 \end{bmatrix}, \quad \frac{\partial}{\partial x_4} = -\frac{i}{b}\frac{\partial}{\partial \tau}, \tag{3.8}$$

it follows that

$$\frac{\partial F_{\alpha\beta}}{\partial x_{\gamma}} + \frac{\partial F_{\beta\gamma}}{\partial x_{\alpha}} + \frac{\partial F_{\gamma\alpha}}{\partial x_{\beta}} = 0, \quad (\alpha, \beta, \gamma = 1, 2, 3, 4), \tag{3.9}$$

is equivalent to the sourceless equations (3.4a) and

$$\frac{\partial F_{\alpha\beta}}{\partial x_{\beta}} = \frac{4\pi}{b} J_{\alpha}, \quad J_{\alpha} = (J_x, J_y, J_z, ib\rho), \tag{3.10}$$

is equivalent to the proper-time equations with sources (3.4b). It should be noted that, in (3.9) and (3.10) and in the sequel, the summation convention is in force for repeated indices. If we now define  $[a_{\mu\nu}]$  by

$$[a_{\mu\nu}] = \begin{bmatrix} 1 + (\gamma - 1)(v_x^2/v^2) & (\gamma - 1)[(v_x v_y)/v^2] & (\gamma - 1)[(v_x v_z)/v^2] & i\gamma \frac{v_x}{c} \\ (\gamma - 1)[(v_x v_y)/v^2] & 1 + (\gamma - 1)(v_y^2/v^2) & (\gamma - 1)[(v_y v_z)/v^2] & i\gamma \frac{v_y}{c} \\ (\gamma - 1)[(v_x v_z)/v^2] & (\gamma - 1)[(v_y v_z)/v^2] & 1 + (\gamma - 1)(v_z^2/v^2) & i\gamma \frac{v_z}{c} \\ -i\gamma \frac{v_x}{c} & -i\gamma \frac{v_y}{c} & \gamma \end{bmatrix}, \quad (3.11)$$

with  $\gamma = [1 - (v/c)^2]^{-1/2}$ ; then the transformations

$$x'_{\mu} = a_{\mu\nu}x_{\nu} \quad (\mu, \nu = 1, 2, 3, 4),$$
 (3.12)

correspond for  $\mu = 1, 2, 3$  to the first set of equations in (2.6) with  $x_4 = i\bar{b}_{\tau}\tau = i\int_0^{\tau} b(s)ds$ . Integrating the first equation in (2.4), we have

$$\int_0^\tau b'(s)ds = \gamma(\mathbf{v}) \left[ \int_0^\tau b(s)ds - \frac{\mathbf{x} \cdot \mathbf{v}}{c} \right]. \tag{3.13}$$

Since the transformations (3.12) are equivalent to our proper-time transformations, we can transform the fields between observers using the four-vector approach just as is commonly done using Lorentz transformations<sup>24,45,46</sup>. Thus, we see that the transformations  $F'_{\mu\nu} = a_{\mu\alpha}a_{\nu\beta}F_{\alpha\beta}$  ( $\mu, \nu, \alpha, \beta = 1, 2, 3, 4$ ) are equivalent to

$$\mathbf{E}' = \gamma \left[ \mathbf{E} + \frac{1}{c} \left( \mathbf{v} \times \mathbf{B} \right) \right] - (\gamma - 1) \frac{(\mathbf{E} \cdot \mathbf{v})}{\mathbf{v}^2} \mathbf{v}, \tag{3.14}$$

$$\mathbf{B'} = \gamma \left[ \mathbf{B} - \frac{1}{c} \left( \mathbf{v} \times \mathbf{E} \right) \right] - (\gamma - 1) \frac{(\mathbf{B} \cdot \mathbf{v})}{\mathbf{v}^2} \mathbf{v}.$$
 (3.15)

It should not be surprising that equations (3.14) and (3.15) are the same as would be obtained if our observers used their own clocks. This is because the transformation coefficient matrix (3.11) is the same as is used for Lorentz transformations between fields. On the other hand, when we look at the current and charge densities, the transformations  $J'_{\mu} = a_{\mu\alpha}J_{\alpha} \quad (\mu, \alpha = 1, 2, 3, 4) \text{ are equivalent to}$ 

$$\mathbf{J}' = \mathbf{J} + (\gamma - 1) \frac{(\mathbf{J} \cdot \mathbf{v})}{\mathbf{v}^2} v - \gamma \frac{b}{c} \rho \mathbf{v}, \tag{3.16a}$$

$$b'\rho' = \gamma(\mathbf{v}) \left[ b\rho - (\mathbf{J} \cdot \mathbf{v}/c) \right]. \tag{3.16b}$$

Using the first equation of (2.4) in (3.16b), we get:

$$\rho' = \frac{\rho - (\mathbf{J} \cdot \mathbf{v}/bc)}{1 - (\mathbf{u} \cdot \mathbf{v}/bc)}.$$
(3.16c)

This result is different from the standard one, (which we obtain if we set b' = b = c in (3.16b)),

$$\rho' = \gamma(\mathbf{v}) \left[ \rho - (\mathbf{J} \cdot \mathbf{v}/c^2) \right]. \tag{3.16d}$$

To see a further difference, if we insert the expression  $\mathbf{J}/c = \rho(\mathbf{u}/b)$  for the current density in (3.16c) and  $\mathbf{J} = \rho \mathbf{w}$  in (3.16d); we obtain

$$\rho' = \rho \frac{1 - (\mathbf{u} \cdot \mathbf{v}/b^2)}{1 - (\mathbf{u} \cdot \mathbf{v}/bc)},\tag{3.17a}$$

$$\rho' = \rho \gamma(\mathbf{v}) \left[ 1 - (\mathbf{w} \cdot \mathbf{v}/c^2) \right]. \tag{3.17b}$$

In order to obtain a sense of the difference between  $\rho$  and  $\rho'$ , assume that

$$u = 2c \approx u', \ b = \sqrt{5}c, \ \Rightarrow w = \frac{2}{\sqrt{5}}c \approx c, \Rightarrow$$
$$\rho' = \rho \left[ \frac{1 - \frac{2v}{5c}}{1 - \frac{2v}{\sqrt{5}c}} \right], \ \& \quad \rho = \rho' \left[ \frac{1 + \frac{2v}{5c}}{1 + \frac{2v}{\sqrt{5}c}} \right].$$

It follows that, unless the relative speed of our two observers is a substantial fraction of c, they will decide that  $\rho = \rho'$ . In fact, we obtain the following remarkable result from equation (3.17a):

**Theorem 3.2** If the source is at rest in the X frame then  $\rho = \rho'$  for all other observers.

**Proof:** The proof is easy, just note that if  $\mathbf{u} = \mathbf{0}$  in X then b = c and, from equation (2.17a),  $\rho = \rho'$ . Since X' is arbitrary, the result is true for all observers.

The above theorem means that, in the proper-time formulation, a spherical charge distribution at rest in any inertial frame will appear spherical to all other inertial observers. As will be shown in the next section, the radiation from an accelerated charged particle appears as a dissipative term in the wave equations for the fields (i.e., neither self-interaction or advanced fields are required). From these two results, we see that the proper-time formulation is independent of particle size or structure.

#### 3.1 Proper-Time Wave Equations

If in equations (3.6), we set

$$\mathbf{B} = \nabla \times \mathbf{A}, \qquad \mathbf{E} = -\frac{1}{b} \frac{\partial \mathbf{A}}{\partial \tau} - \nabla \Phi, \tag{3.18}$$

then we obtain

$$\nabla \left[ \nabla \cdot \mathbf{A} + \frac{1}{b} \frac{\partial \Phi}{\partial \tau} \right] + \frac{1}{b} \frac{\partial}{\partial \tau} \left[ \frac{1}{b} \frac{\partial \mathbf{A}}{\partial \tau} \right] - \nabla^2 \mathbf{A} = \frac{1}{b} \left( 4\pi \rho \mathbf{u} \right), \tag{3.19}$$

and

$$-\nabla^2 \Phi - \frac{1}{b} \frac{\partial}{\partial \tau} \left[ \nabla \cdot \mathbf{A} \right] = 4\pi \rho. \tag{3.20}$$

Imposing the (proper-time) Lorentz gauge

$$\nabla \cdot \mathbf{A} + \frac{1}{b} \frac{\partial \Phi}{\partial \tau} = 0, \tag{3.21}$$

we get the wave equations

$$\frac{1}{b^2} \frac{\partial^2 \mathbf{A}}{\partial \tau^2} - \frac{1}{b^4} (\mathbf{u} \cdot \mathbf{a}) \frac{\partial \mathbf{A}}{\partial \tau} - \nabla^2 \mathbf{A} = \frac{1}{b} [4\pi \rho \mathbf{u}], \qquad (3.22a)$$

$$\frac{1}{b^2} \frac{\partial^2 \Phi}{\partial \tau^2} - \frac{1}{b^4} (\mathbf{u} \cdot \mathbf{a}) \frac{\partial \Phi}{\partial \tau} - \nabla^2 \Phi = 4\pi \rho. \tag{3.22b}$$

We thus obtain a new term that arises because the proper-time of the source carries information about the interaction that is not available when the proper-time of the observer is used in formulating theory. In Section 5 the wave equations will be derived for the fields directly to get (no gauge required):

$$\frac{1}{b^2} \frac{\partial^2 \mathbf{E}}{\partial \tau^2} - \frac{1}{b^4} (\mathbf{u} \cdot \mathbf{a}) \frac{\partial \mathbf{E}}{\partial \tau} - \nabla^2 \mathbf{E} = -\nabla \left[ 4\pi \rho \mathbf{u} \right] - \frac{1}{b} \frac{\partial}{\partial \tau} \left[ \frac{4\pi \mathbf{J}}{b} \right], \tag{3.23a}$$

$$\frac{1}{b^2} \frac{\partial^2 \mathbf{B}}{\partial \tau^2} - \frac{1}{b^4} (\mathbf{u} \cdot \mathbf{a}) \frac{\partial \mathbf{B}}{\partial \tau} - \nabla^2 \mathbf{B} = \frac{1}{b} \frac{\partial}{\partial \tau} \left[ \frac{4\pi \nabla \times \mathbf{J}}{b} \right]. \tag{3.23b}$$

Thus, the new term is independent of the gauge. The physical interpretation is clear, this is a dissipative term which is zero if a is zero or orthogonal to u. Furthermore, it arises instantaneously with the acceleration of the source. This is exactly what one expects of the radiation caused by the inertial resistance of the source to accelerated motion and is precisely what one means by radiation reaction (see Wheeler and Feynman<sup>25</sup>). It should be noted that the creation of real physical conditions which will make a orthogonal to u is almost impossible since a arises because of an external force and has no relationship to u. In order to get some insight into the meaning of the new dissipative terms, let us focus on equation (3.22b). If we use p = m0u from equation (2.9b), we see that the external force Fext satisfies (This is only approximate as will be seen in Section 5.4, equation (5.58).)

$$\mathbf{F}_{ext} = \frac{d\mathbf{p}}{d\tau} = m_0 \mathbf{a},\tag{3.24a}$$

so that equation (3.22b) becomes

$$\frac{1}{b^2} \frac{\partial^2 \Phi}{\partial \tau^2} - \left(\frac{\mathbf{u}}{b}\right) \cdot \left(\frac{\mathbf{F}_{ext}}{m_0 b^2}\right) \left(\frac{1}{b} \frac{\partial \Phi}{\partial \tau}\right) - \nabla^2 \Phi = 4\pi \rho. \tag{3.24b}$$

If we identify m0b <sup>2</sup> with the effective interaction energy of the particle, then the middle term can be interpreted as the reactive power loss per unit interaction energy of the particle due to its resistance to Fext. To see this additional term in another physically important way, use the change of variables Φ = (b/c) 1/2 g in (3.22b) to get (see Courant and Hilbert<sup>47</sup>)

$$\frac{1}{b^2} \frac{\partial^2 g}{\partial \tau^2} - \nabla^2 g + \left[ \frac{\ddot{b}}{2b^3} - \frac{5\dot{b}^2}{4b^4} \right] g = 4\pi\rho \left( \frac{c}{b} \right)^{1/2}. \tag{3.24c}$$

This is the Klein-Gordon equation with an effective mass µ given by

$$\mu = \left\{ \frac{\hbar^2}{b^2} \left[ \frac{\ddot{b}}{2b^3} - \frac{5\dot{b}^2}{4b^4} \right] \right\}^{1/2}.$$
 (3.25)

Hence, the reactive power loss per unit interaction energy in (3.24b) is equivalent to an effective mass for the photon that depends on the external force acting on the particle.

We have only considered our equations at the source. If we look at them in a region outside the source, there is a major change. The dissipative term is now constant with its value fixed at the time the radiation left the source. Thus, a new picture emerges. Every

accelerated charged particle emits a continuous stream of (very) small particles (photons) in all directions. The energy and the velocity of the particles depend on the velocity of the source at the moment of emission. The velocity of the particles remains constant until they are scattered or absorbed.

### 3.2 Radiation From An Accelerated Charge

In this section, we compute the radiation from an accelerated charge using the propertime theory. We can solve equation (3.24b) directly, but a better approach is to first find the solution using the proper-time of the observer and then transform the result to the proper-time of the source. This makes the computations easier to follow and gives the result quicker. We follow closely the approach in Panofsky and Phillips<sup>24</sup>. In this section,  $(\mathbf{x}(t),t)$  represents the field position and  $(\mathbf{x}'(t'),t')$  represents the retarded position of a point charge source q, with  $\mathbf{r} = \mathbf{x} - \mathbf{x}'$ ,  $d\mathbf{r}/dt' = -\mathbf{w}$ , and  $d^2\mathbf{r}/dt'^2 = \dot{\mathbf{w}}$ . The field solutions using the standard Lienard-Wiechert potentials are given by

$$\mathbf{A} = \frac{q\mathbf{w}}{cs}, \qquad \Phi = \frac{q}{s}, \qquad s = r - \left(\frac{\mathbf{r} \cdot \mathbf{w}}{c}\right). \tag{3.26}$$

The proper-time form is obtained by replacing  $\mathbf{w}/c$  by  $\mathbf{u}/b$  to get

$$\mathbf{A} = \frac{q\mathbf{u}}{bs}, \quad \Phi = \frac{q}{s}, \quad s = r - \left(\frac{\mathbf{r} \cdot \mathbf{u}}{b}\right).$$
 (3.27)

The field and source-point variables are related by the condition

$$r = |\mathbf{x} - \mathbf{x}'| = c(t - t'). \tag{3.28}$$

Here,  $d\mathbf{r}/d\tau' = -\mathbf{u} = -d\mathbf{x}'/d\tau'$ , where  $\tau'$  denotes the retarded proper-time of the source. The corresponding  $\mathbf{E}$  and  $\mathbf{B}$  fields can be computed using equation (3.18) in the form

$$\mathbf{E}(\mathbf{x},\tau) = -\frac{1}{\overline{b}} \frac{\partial \mathbf{A}(\mathbf{x},\tau)}{\partial \tau} - \nabla \Phi(\mathbf{x},\tau), \ \mathbf{B}(\mathbf{x},\tau) = \nabla \times \mathbf{A}(\mathbf{x},\tau)$$
(3.29)

with  $\bar{\mathbf{u}} = d\mathbf{x}/d\tau$ , where  $\tau$  denotes the proper-time of the present position of the source and  $\bar{b} = (\bar{\mathbf{u}}^2 + c^2)^{1/2}$ . In order to compute the fields from the potentials, we note that the components of the  $\nabla$  operator are partials at constant time  $\tau$ , and therefore are not at constant  $\tau'$ . Also, the partial derivatives with respect to  $\tau$  imply constant  $\mathbf{x}$  and hence refer to the comparison of potentials at a given point over an interval in which the coordinates of the source will have changed. Since only time variations with respect to  $\tau'$  are given, we must transform  $(\partial/\partial\tau)|_{\mathbf{x}}$  and  $\nabla|_{\tau}$  to expressions in terms of  $\partial/\partial\tau'|_{\mathbf{x}}$ . To do this, we must first transform (3.28) into a relationship between  $\tau$  and  $\tau'$ . The required correspondence is

$$c(t - t') = \int_{\tau'}^{\tau} b(s)ds. \tag{3.30}$$

It is easier to first relate  $\partial/\partial t|_{\mathbf{x}}$  to  $\partial/\partial t'|_{\mathbf{x}}$  and then convert them to relationships between  $\partial/\partial \tau|_{\mathbf{x}}$  and  $\partial/\partial \tau'|_{\mathbf{x}}$ . The following are in reference 24, pg. 298:

$$\frac{\partial r}{\partial t'} = -\frac{\mathbf{r} \cdot \mathbf{w}}{r}, \quad \frac{\partial r}{\partial t} = c \left( 1 - \frac{\partial t'}{\partial t} \right) = \frac{\partial r}{\partial t'} \cdot \frac{\partial t'}{\partial t} = -\frac{\mathbf{r} \cdot \mathbf{w}}{r} \frac{\partial t'}{\partial t}. \tag{3.31}$$

Since  $\partial \tau / \partial t = c/b$ , we have

$$\frac{\partial r}{\partial t} = c \frac{\partial}{\partial t} (t - t') = \frac{\partial \tau}{\partial t} \frac{\partial}{\partial \tau} \int_{\tau'}^{\tau} b(s) ds = \frac{c}{\bar{b}} \left[ \bar{b} - b \frac{\partial \tau'}{\partial \tau} \right]. \tag{3.32}$$

We also have, using  $\partial \tau'/\partial t'=c/b$  , that

$$\frac{\partial r}{\partial t'} = \frac{\partial r}{\partial \tau'} \frac{\partial \tau'}{\partial t'} = \frac{c}{b} \frac{\partial r}{\partial \tau'} \Rightarrow \frac{1}{b} \frac{\partial r}{\partial \tau'} = -\frac{\mathbf{r} \cdot \mathbf{w}}{rc} = -\frac{\mathbf{r} \cdot \mathbf{u}}{rb}, \tag{3.33}$$

so  $\partial r/\partial \tau' = -\mathbf{r} \cdot \mathbf{u}/r$  and hence

$$\frac{\partial r}{\partial t} = \frac{\partial r}{\partial \tau} \frac{c}{\bar{b}} = \frac{c}{\bar{b}} \left[ \bar{b} - b \frac{\partial \tau'}{\partial \tau} \right] \Rightarrow \frac{\partial r}{\partial \tau} = \left[ \bar{b} - b \frac{\partial \tau'}{\partial \tau} \right], \tag{3.34}$$

$$\frac{\partial r}{\partial \tau} = \frac{\partial r}{\partial \tau'} \frac{\partial \tau'}{\partial \tau} = -\frac{\mathbf{r} \cdot \mathbf{u}}{r} \frac{\partial \tau'}{\partial \tau} \Rightarrow -\frac{\mathbf{r} \cdot \mathbf{u}}{r} \frac{\partial \tau'}{\partial \tau} = \left[ \bar{b} - b \frac{\partial \tau'}{\partial \tau} \right]. \tag{3.35}$$

Solving (3.35) for  $\partial \tau'/\partial \tau$ , we get

$$\frac{\partial \tau'}{\partial \tau} = \frac{\bar{b}}{b} \frac{r}{s}, \quad s = r - \frac{\mathbf{r} \cdot \mathbf{u}}{b}. \tag{3.36}$$

Using this, we see that

$$\frac{1}{\bar{b}}\frac{\partial}{\partial \tau} = \frac{1}{b} \cdot \frac{r}{s}\frac{\partial}{\partial \tau'}.$$
(3.37)

From  $\nabla r = -c\nabla t' = \nabla_1 r + (\partial r/\partial t')\nabla t'$ , we see that

$$\nabla r = \frac{\mathbf{r}}{r} - \frac{c}{b} \cdot \frac{\mathbf{r} \cdot \mathbf{u}}{r} \nabla t' \quad \Rightarrow -c \nabla t' = \frac{\mathbf{r}}{r} - \frac{c}{b} \cdot \frac{\mathbf{r} \cdot \mathbf{u}}{r} \nabla t'. \tag{3.38}$$

Using  $c\nabla t' = b\nabla \tau'$  and solving for  $\nabla \tau'$ , we get  $\nabla \tau' = -(\mathbf{r}/bs)$ , so that

$$\nabla = \nabla_1 - \frac{\mathbf{r}}{bs} \cdot \frac{\partial}{\partial \tau'}.$$
 (3.39)

We now compute  $\nabla_1 s$  and  $\partial s/\partial \tau'$ . The calculations are easy, so we simply state the results:

$$\nabla_1 s = \frac{\mathbf{r}}{r} - \frac{\mathbf{u}}{b} = \frac{1}{r} \left( \mathbf{r} - \frac{r\mathbf{u}}{b} \right), \tag{3.40}$$

$$\frac{\partial s}{\partial \tau'} = \frac{\mathbf{u}^2}{b} - \frac{\mathbf{r} \cdot \mathbf{u}}{r} - \frac{\mathbf{r} \cdot \mathbf{a}}{b} + \frac{(\mathbf{r} \cdot \mathbf{u}) (\mathbf{u} \cdot \mathbf{a})}{b^3}.$$
 (3.41)

We can now calculate the fields. The computations are long but follow those of reference 24, so we only record a few selected results. We obtain

$$-\nabla\Phi = \frac{q}{s^{2}}\nabla s = \frac{q}{s^{2}}\left(\nabla_{1}s - \frac{\mathbf{r}}{bs} \cdot \frac{\partial s}{\partial \tau}\right) \Rightarrow$$

$$-\nabla\Phi = \frac{q\left[\mathbf{r}\left(1 - \mathbf{u}^{2}/b^{2}\right) - \mathbf{u}s/b\right]}{s^{3}} + \frac{q\mathbf{r}\left(\mathbf{r} \cdot \mathbf{a}\right)}{b^{2}s^{3}} - \frac{q\mathbf{r}\left(\mathbf{r} \cdot \mathbf{u}\right)\left(\mathbf{u} \cdot \mathbf{a}\right)}{b^{4}s^{3}}.$$
(3.42)

Now use equation (3.37) to get

$$-\frac{1}{\overline{b}}\frac{\partial \mathbf{A}}{\partial \tau} = \left(-\frac{1}{b}\right) \left(\frac{r}{s}\right) \frac{\partial \mathbf{A}}{\partial \tau'} \Rightarrow$$

$$-\frac{1}{\overline{b}}\frac{\partial A}{\partial \tau} = \frac{-\left(qr\mathbf{u}/b\right)\left\{\left(\mathbf{u}/b\right)\cdot\left[\left(\mathbf{r}/r\right)-\left(\mathbf{u}/b\right)\right]\right\}}{s^{3}} + \frac{-qr^{2}\mathbf{a} + qr\left\{\mathbf{r}\times\left[\mathbf{a}\times\left(\mathbf{u}/b\right)\right]\right\}}{b^{2}s^{3}} + \frac{q\mathbf{u}\left[\left(\mathbf{r}\cdot\mathbf{r}\right)\left(\mathbf{u}\cdot\mathbf{a}\right)\right]}{b^{4}s^{3}}.$$
(3.43)

Combining (3.42) and (3.43), we get

$$\mathbf{E}(\mathbf{x},\tau) = -\frac{1}{\bar{b}} \frac{\partial \mathbf{A}(\mathbf{x},\tau)}{\partial \tau} - \nabla \Phi(x,\tau) \Rightarrow$$

$$\mathbf{E}(\mathbf{x},\tau) = \frac{q \left[ \mathbf{r} \left( 1 - \mathbf{u}^2 / b^2 \right) - \mathbf{u} s / b \right]}{s^3} - \frac{\left( -q r \mathbf{u} / b \right) \left[ (\mathbf{u} / b) \cdot (\mathbf{r} / r - \mathbf{u} / b) \right]}{s^3}$$

$$+ \frac{-q \left[ r^2 \mathbf{a} - \mathbf{r} \left( \mathbf{r} \cdot \mathbf{a} \right) \right] + q r \left[ \mathbf{r} \times (\mathbf{a} \times \mathbf{u} / b) \right]}{b^2 s^3} + \frac{q \left( \mathbf{u} \cdot \mathbf{a} \right) \left[ \mathbf{u} r^2 - \mathbf{r} \left( \mathbf{r} \cdot \mathbf{u} \right) \right]}{b^4 s^3}.$$
(3.44)

Finally, using standard vector identities and combining terms, we get (with  $\mathbf{r_u} = \mathbf{r} - \mathbf{u}r/b$ )

$$\mathbf{E}(\mathbf{x}, \tau) = \frac{q \left[ \mathbf{r}_{\mathbf{u}} \left( 1 - \mathbf{u}^2 / b^2 \right) \right]}{s^3} + \frac{q \left\{ \mathbf{r} \times \left[ \mathbf{r}_{\mathbf{u}} \times \mathbf{a} \right] \right\}}{b^2 s^3} + \frac{q \left( \mathbf{u} \cdot \mathbf{a} \right) \left[ \mathbf{r} \times \left( \mathbf{u} \times \mathbf{r} \right) \right]}{b^4 s^3}.$$
(3.45)

The computation of  $\mathbf{B}$  is similar:

$$\mathbf{B}(\mathbf{x},\tau) = \frac{q\left[(\mathbf{r} \times \mathbf{r_u})(1 - \mathbf{u}^2/b^2)\right]}{rs^3} + \frac{q\mathbf{r} \times \{\mathbf{r} \times [\mathbf{r_u} \times \mathbf{a}]\}}{rb^2s^3} + \frac{qr(\mathbf{u} \cdot \mathbf{a})(\mathbf{r} \times \mathbf{u})}{b^4s^3}.$$
(3.46)

It is easy to see that we have  $\mathbf{B} = (\mathbf{r}/r) \times \mathbf{E}$  so that  $\mathbf{B}$  is orthogonal to  $\mathbf{E}$ . The first two terms in (3.45) and (3.46) are the same as (19-13) and (19-14) in reference 24 (pg. 299). The last term in each case arises because of the dissipative terms in equations (3.22) and (3.23).

The last terms in (3.45) and (3.46) are zero if **a** is zero or orthogonal to **u**. In the first case, there is no radiation and the particle moves with constant velocity so that the field is massless. As noted earlier, the second case depends on conditions that are impossible in practice, namely the creation of motion which keeps **a** orthogonal to **u**. Since  $\mathbf{r} \times (\mathbf{u} \times \mathbf{r}) = r^2\mathbf{u} - (\mathbf{u} \cdot \mathbf{r})\mathbf{r}$ , we see that there is a component along the direction

of propagation (longitudinal). Hence, in all other cases, there is a small mass associated with electromagnetic radiation which varies with the acceleration of the particle.

# 3.3 Radiated Energy

In light of the difference in the calculated fields, it becomes important to also compute the radiated energy for the proper-time theory and compare it with the Minkowski formulation. It is well-known that the radiated energy is determined by the Poynting vector, which is defined by  $\mathbf{P} = (c/4\pi) (\mathbf{E} \times \mathbf{B})$ .

To calculate the angular distribution of the radiated energy, we must be careful to note that the rate of radiation is the amount of energy lost by the charge in a time interval  $d\tau'$  during the emission of the signal  $(-dU/d\tau')$ . However (at a field point), the Poynting vector  $\mathbf{P}$  represents the energy flow per unit time measured at the present time  $(\tau)$ . With this understanding, the same approach that leads to the above formula gives  $\mathbf{P} = (\bar{b}/4\pi) (\mathbf{E} \times \mathbf{B})$  in the proper-time formulation. We thus obtain the rate of energy loss of a charged particle into a given infinitesimal solid angle  $d\Omega$  as

$$-\frac{dU}{d\tau'}(\Omega)d\Omega = (\bar{b}/4\pi) \left[\mathbf{n} \cdot (\mathbf{E} \times \mathbf{B})\right] \mathbf{r}^2 \frac{d\tau}{d\tau'} d\Omega. \tag{3.47}$$

Using equation (3.36), we get that  $(d\tau/d\tau') = bs/\bar{b}r$ , so that (3.47) becomes

$$-\frac{dU}{d\tau'}(\Omega)d\Omega = (b/4\pi)\left[\mathbf{n}\cdot(\mathbf{E}\times\mathbf{B})\right]rsd\Omega. \tag{3.48}$$

As is well-known, only those terms that fall off as (1/r) (the radiation terms) in (3.45) and (3.46) contribute to the integral of (3.48). It is easy to see that our theory gives the following radiation terms:

$$\mathbf{E}_{rad} = \frac{q\left\{\mathbf{r} \times \left[\mathbf{r}_{\mathbf{u}} \times \mathbf{a}\right]\right\}}{b^2 s^3} + \frac{q\left(\mathbf{u} \cdot \mathbf{a}\right)\left[\mathbf{r} \times \left(\mathbf{u} \times \mathbf{r}\right)\right]}{b^4 s^3} = \mathbf{E}_{rad}^c + \mathbf{E}_{rad}^d, \tag{3.49}$$

$$\mathbf{B}_{rad} = \frac{q\mathbf{r} \times \{\mathbf{r} \times [\mathbf{r}_{\mathbf{u}} \times \mathbf{a}]\}}{rb^2s^3} + \frac{qr(\mathbf{u} \cdot \mathbf{a})(\mathbf{r} \times \mathbf{u})}{b^4s^3} = \mathbf{B}_{rad}^c + \mathbf{B}_{rad}^d, \tag{3.50}$$

where  $\mathbf{E}_{rad}^c$ ,  $\mathbf{B}_{rad}^c$  are of the same form as the classical terms with c replaced by b,  $\mathbf{w}'$  by  $\mathbf{u}$ , and  $\dot{\mathbf{w}}'$  by  $\mathbf{a}$ . The two terms  $\mathbf{E}_{rad}^d$ ,  $\mathbf{B}_{rad}^d$ , are new and come directly from the dissipation term in the wave equations. (Note the characteristic  $(\mathbf{u} \cdot \mathbf{a})/b^4$ .) We can easily integrate the classical terms to see that

$$\iint_{\Omega} (-dU^{c}/d\tau)d\Omega$$

$$= (b/4\pi) \iint_{\Omega} \left[ \mathbf{n} \cdot (\mathbf{E}_{rad}^{c} \times \mathbf{B}_{rad}^{c}) \right] rsd\Omega = \frac{2}{3} \frac{q^{2} |\mathbf{a}|^{2}}{b^{3}}.$$
(3.51)

This agrees with the standard result for small proper-velocity and proper-acceleration of the charge when  $b \approx c$  and  $\mathbf{a} \approx d\mathbf{w}/dt$ .

In the general case, our theory gives additional effects because of the dissipative terms. To compute the integral of (3.48), we use spherical coordinates with the proper-velocity  $\mathbf{u}$  directed along the positive z-axis. Without loss of generality, we orient the coordinate system so that the proper-acceleration  $\mathbf{a}$  lies in the xz-plane. Let  $\alpha$  denote the acute angle between  $\mathbf{a}$  and  $\mathbf{u}$ , and substitute (3.49) and (3.50) in (3.48) to obtain

$$-\frac{dU}{d\tau}(\Omega)d\Omega =$$

$$= \frac{q^2 |\mathbf{a}|^2}{4\pi b^3} \left\{ (1 - \beta \cos \theta)^{-4} \left[ 1 - \sin^2 \theta \sin^2 \alpha \cos \phi \right] - \cos^2 \theta \cos^2 \alpha - (1/2) \sin 2\theta \sin 2\alpha \cos \phi \right]$$

$$-2\beta (1 - \beta \cos \theta)^{-5} \left( \sin^2 \theta \cos \alpha - (1/2) \sin 2\theta \sin \alpha \cos \phi \right) \chi$$

$$+\beta^2 \sin^2 \theta (1 - \beta \cos \theta)^{-6} \chi^2 \right\},$$
(3.52)

where

$$\chi = \frac{b^2}{r|\mathbf{a}|} \left( 1 - \beta^2 \right) + \beta \cos \alpha \left( 1 - \frac{1}{\beta} \cos \theta \right) - \sin \theta \sin \alpha \cos \phi, \tag{3.53}$$

and  $\beta = (|\mathbf{u}|/b)$ .

The integration of (3.52) over the surface of the sphere is elementary, and we obtain, after some extensive but easy computations (which are summarized in the appendix):

$$\lim_{r \to \infty} \iint -\frac{dU}{d\tau}(\Omega) d\Omega$$

$$= \frac{2}{3} \frac{q^2 |\mathbf{a}|^2}{b^3} (1 - \beta^2)^{-3} \left[ 1 - \frac{1}{5} \beta^2 (4 + \beta^2) + \frac{1}{5} \beta^2 (6 + \beta^2) \sin^2 \alpha \right].$$
(3.54)

As can be seen, this result agrees with (3.51) at the lowest order. For comparison, the same calculation using the observer's clock for the case of general orientation of velocity  $d\mathbf{x}'/dt'$  and acceleration  $d\mathbf{w}'/dt'$  is

$$\lim_{r \to \infty} \iint -\frac{dU}{dt}(\Omega) d\Omega = \frac{2}{3} \frac{q^2 |\dot{\mathbf{w}}'|^2}{c^3} (1 - \beta^2)^{-3} [1 - \beta^2 \sin^2 \alpha], \qquad (3.55)$$

where  $\beta = (|\mathbf{w}'|/c)$ .

We observe that, in general, for an arbitrary angle  $\alpha$  with  $0 \le \alpha \le \pi/2$  and arbitrary  $\beta$  between 0 and 1, our result does not agree with (3.55) even if we replace b with c and a with  $d\mathbf{w}'/dt'$ . This shows, along with our other results, that the apparently small change in clocks induces large changes in the physical predictions. We will return to this point in the conclusion of the paper.

### 4.0 Proper-Time Doppler Effect and Aberration

In this section, we apply our proper-time theory to compute the optical Doppler effect and aberration. To do this, we first consider the transformation properties of plane wave solutions to Maxwell's equations. Assuming that our observers are in the far-field of the source so that, to a good approximation, the waves are plane when they arrive at the observers' positions, we want solutions of the form ( $\mathbf{E}_0 = const$ ,  $\mathbf{B}_0 = const$ )

$$\mathbf{E} = \Re \left\{ \mathbf{E}_0 \exp \left[ i \left( \mathbf{k} \cdot \mathbf{x} - \frac{1}{c} \int_0^\tau \omega(s) b(s) ds \right) \right] \right\}, \tag{4.1a}$$

$$\mathbf{B} = \Re \left\{ \mathbf{B}_0 \exp \left[ i \left( \mathbf{k} \cdot \mathbf{x} - \frac{1}{c} \int_0^\tau \omega(s) b(s) ds \right) \right] \right\}, \tag{4.1b}$$

where, in accordance with equations (2.0), we have modified the plane wave representations to allow for proper-time (nonlocal) dependence of the frequency. Assuming that the frequency is a differentiable function of time, we get that the above plane wave representations of the fields are solutions of the wave equations in the far-field region (where the charge and current densities are zero),

$$\frac{1}{b^2} \frac{\partial^2 \mathbf{E}}{\partial \tau^2} - \frac{1}{b^4} (\mathbf{u} \cdot \mathbf{a}) \frac{\partial \mathbf{E}}{\partial \tau} - \nabla^2 \mathbf{E} = 0, \tag{3.23a}$$

$$\frac{1}{b^2} \frac{\partial^2 \mathbf{B}}{\partial \tau^2} - \frac{1}{b^4} (\mathbf{u} \cdot \mathbf{a}) \frac{\partial \mathbf{B}}{\partial \tau} - \nabla^2 \mathbf{B} = 0, \tag{3.23b}$$

provided that

$$\mathbf{k}^2 = \frac{\omega(\tau)^2}{c^2} \left[ 1 + i \frac{c\dot{\omega}(\tau)}{b\omega(\tau)^2} \right]. \tag{4.2a}$$

In addition, from (3.4) we have

$$\mathbf{k} \cdot \mathbf{B}_0 = 0, \quad \mathbf{k} \cdot \mathbf{E}_0 = 0, \tag{4.2b}$$

$$\mathbf{k} \times \mathbf{E}_0 = \frac{\omega(\tau)}{c} \mathbf{B}_0. \tag{4.2c}$$

It follows from (4.2a) that the wave vector  $\mathbf{k}$  depends on  $\omega(\tau)$  and its derivative  $\dot{\omega}(\tau)$ .

To obtain the transformation properties of the plane waves, we use (3.14) and (3.15) along with (4.1) to get

$$\mathbf{E}' = \Re \left\{ \mathbf{E}'_0 \exp \left[ i \left( \mathbf{k} \cdot \mathbf{x} - \frac{1}{c} \int_0^\tau \omega(s) b(s) ds \right) \right] \right\}, \tag{4.3a}$$

$$\mathbf{B}' = \Re \left\{ \mathbf{B}'_0 \exp \left[ i \left( \mathbf{k} \cdot \mathbf{x} - \frac{1}{c} \int_0^\tau \omega(s) b(s) ds \right) \right] \right\}, \tag{4.3b}$$

with

$$\mathbf{E}_0' = \gamma \left[ \mathbf{E}_0 + \frac{1}{c} \left( \mathbf{v} \times \mathbf{B}_0 \right) \right] - (\gamma - 1) \frac{\left( \mathbf{E}_0 \cdot \mathbf{v} \right)}{\mathbf{v}^2} \mathbf{v}. \tag{4.3c}$$

$$\mathbf{B}_0' = \gamma \left[ \mathbf{B}_0 - \frac{1}{c} \left( \mathbf{v} \times \mathbf{E}_0 \right) \right] - (\gamma - 1) \frac{(\mathbf{B}_0 \cdot \mathbf{v})}{\mathbf{v}^2} \mathbf{v}. \tag{4.3d}$$

We now use the inverse transformations (2.1a), (2.3b), and (2.4) to transform the phase

$$\Phi = i \left( \mathbf{k} \cdot \mathbf{x} - (1/c) \int_0^\tau \omega(s) b(s) ds \right)$$
 (4.3*e*)

in (4.3a) and (4.3b) to the corresponding expression in the primed variables:

$$\Phi' = i \left( \mathbf{k}' \cdot \mathbf{x}' - (1/c) \int_0^\tau \omega'(s) b'(s) ds \right), \tag{4.3}f$$

where the wave number  $\mathbf{k}'$  and the frequency  $\omega'(s)$  are to be determined by the requirement that the transformed phase  $\Phi'$  has the indicated form. Substituting (2.1b) and (2.4) into (4.3e), we get

$$\Phi = i \left[ \left( \mathbf{k} + (\gamma(\mathbf{v}) - 1) \left( \frac{\mathbf{k} \cdot \mathbf{v}}{||\mathbf{v}||^{2}} \right) \mathbf{v} \right) \cdot \mathbf{x}' + \gamma(\mathbf{v}) \frac{\mathbf{k} \cdot \mathbf{v}}{c} \int_{0}^{\tau} b'(s) ds - \frac{1}{c} \int_{0}^{\tau} \omega(s) b'(s) ds \right] \\
= i \left[ \left( \mathbf{k} + (\gamma - 1) \left( \frac{\mathbf{k} \cdot \mathbf{v}}{||\mathbf{v}||^{2}} \right) \mathbf{v} \right) \cdot \mathbf{x}' - \frac{\gamma}{c} \int_{0}^{\tau} (\omega(s) - \mathbf{k} \cdot \mathbf{v}) b'(s) ds - \frac{\gamma}{c^{2}} \int_{0}^{\tau} \omega(s) \mathbf{u}' \cdot \mathbf{v} ds \right].$$
(4.4)

Integrating the last term in (4.4) by parts, we obtain the desired form for  $\Phi'$ , where the frequency relation is given by

$$\omega'(\tau) = \gamma(\omega(\tau) - \mathbf{k} \cdot \mathbf{v}), \tag{4.5}$$

and the wave number relation (contributing the nonlocal part to  $\Phi'$ ) is given by:

$$\mathbf{k}' \cdot \mathbf{x}'(\tau) = \mathbf{k} \cdot \mathbf{x}'(\tau) + (\gamma - 1) \left[ \frac{(\mathbf{k} \cdot \mathbf{v})(\mathbf{v} \cdot \mathbf{x}'(\tau))}{||\mathbf{v}||^2} \right] - \frac{\gamma \omega(\tau)}{c^2} (\mathbf{v} \cdot \mathbf{x}'(\tau)) + \frac{\gamma \omega(0)}{c^2} (\mathbf{v} \cdot \mathbf{x}'(0)) + \frac{\gamma}{c^2} \int_0^{\tau} \frac{d\omega(s)}{ds} \left[ \mathbf{v} \cdot \mathbf{x}'(s) \right] ds.$$

$$(4.6)$$

The wave vectors in our two frames differ by an extra nonlocal term compared to the standard result, while the transformations of the frequencies (4.5) agree with the normal case except for the  $\tau$  dependence. This nonlocal term occurs because we allowed the frequency of the wave to vary. It is easy to check that, if  $\omega$  is constant (and the source passes though the origin(s) at  $\tau = 0$ ), we get the standard result.

We now consider the planar representation (4.6) with the velocity  $\mathbf{v}$  taken along the  $\mathbf{x} = \mathbf{x}'$  axes with angle  $\theta$  defined as that between  $\mathbf{k}$  and  $\mathbf{v}$ , and  $\theta'$  the angle between  $\mathbf{k}'$  and  $\mathbf{v}$ ,  $\omega$  constant, and assume that the source passes though the origin(s) at  $\tau = 0$ . We then obtain from (4.6) the following relations between the angles  $\theta$  and  $\theta'$ :

$$k'\cos\theta' = \gamma k\cos\theta - \gamma \frac{\omega}{c^2}v,\tag{4.7}$$

$$k'\sin\theta' = k\sin\theta. \tag{4.8}$$

They combine in the standard manner to give

$$\tan \theta' = \frac{1}{\gamma} \frac{\sin \theta}{\cos \theta - \frac{v}{c} \frac{\omega}{kc}}.$$
 (4.9a)

This is the standard result for the aberration of wave vectors due to the relative motion of the two reference frames. It should be noted that we have not assumed that the X frame is at rest relative to the medium. Furthermore, we see from (4.2a) that  $kc = \omega$  in free space (under the above assumptions). In general,  $kc = \omega(\tau) \left[1 + i\left(c\dot{\omega}(\tau)/b\omega(\tau)^2\right)\right]^{1/2}$  so that our theory allows for nonlocal effects.

For any homogeneous medium,  $\omega/\Re k$  is equal to the phase velocity,  $v_{ph}$ , of the wave,

$$v_{ph} = c\Re\left[1 + i\left(c\dot{\omega}(\tau)/b\omega(\tau)^2\right)\right]^{-1/2},\tag{4.10}$$

and  $c/v_{ph}$  is defined to be the index of refraction, n, of the medium. Thus, (4.9a) becomes:

$$\tan \theta' = \frac{1}{\gamma} \frac{\sin \theta}{\cos \theta - \frac{v}{cn}}.$$
 (4.9b)

This is what we would normally expect from the standard theory. However, the importance of (4.10) becomes clear when we consider the group velocity, rather than the phase velocity, of electromagnetic waves. As is well-known, the group velocity represents the rate of energy transmission, and is defined by  $v_g = \Re(d\omega/dk)$ . We know that use of observer clocks (proper-times) gives  $v_g = v_g' = c$ . The question is, what is this relationship in the source proper-time theory?

To determine how  $v_g$  is related to  $v'_g$ , we restrict ourselves to the case when the waves are moving parallel to the motion of the X' frame relative to the X frame, so that the wave vectors  $\Re \mathbf{k}$  and  $\Re \mathbf{k}'$  are parallel to the velocity  $\mathbf{v}$ . Then the frequency and wave number relations (4.5) and (4.6) become (under these conditions)

$$\omega'(\tau) = \gamma(\omega(\tau) - \mathbf{k} \cdot \mathbf{v}), \tag{4.11}$$

$$k'x' = \gamma \left(k - \frac{\gamma v\omega(\tau)}{c^2}\right)x'(\tau) + \frac{\gamma v\omega(0)}{c^2}x'(0) + \frac{\gamma v}{c^2} \int_0^{\tau} \frac{d\omega(s)}{ds} \left[x'(s)\right]ds, \tag{4.12}$$

where, in the last equation, we have replaced the vector  $\mathbf{x}'(\tau)$  by the scalar  $x'(\tau)$  because we are only interested in the  $\tau$  dependence of the frequencies and wave numbers.

Defining the group velocity in the X, X' frames by

$$v_g \equiv \Re(\frac{d\omega}{dk}) = \Re(\frac{d\omega}{d\tau} / \frac{dk}{d\tau}), \qquad v_g' \equiv \Re(\frac{d\omega'}{dk'}) = \Re(\frac{d\omega'}{d\tau} / \frac{dk'}{d\tau}),$$
 (4.13)

we obtain from (4.11) the equation

$$\frac{d\omega'}{d\tau} = \gamma \left(\frac{d\omega}{d\tau} - v\frac{dk}{d\tau}\right),\tag{4.14}$$

and from (4.12) (after canceling terms),

$$\frac{dk'}{d\tau}x'(\tau) = \gamma \frac{dk}{d\tau}x'(\tau). \tag{4.15}$$

Substitution of (4.14) and (4.15) into (4.13) gives the relation

$$v_g = v_g' - v \tag{4.16}$$

between the group velocities in the X and X' frames respectively. It is clear that, if the group velocity of the source has the value c in one frame, it will not have that value in the other frame and, indeed, may have a larger value. Furthermore, the Doppler formula (4.11) can be written as

$$\omega'(\tau) = \gamma \omega(\tau) (1 - \beta n \left[\omega(\tau)\right] \cos \theta), \tag{4.17}$$

where we have used  $\beta = v/c$ ,  $k = \omega/v_{ph}$ , and  $n = c/v_{ph}$ . Because of (4.10), this is a nonlinear relationship.

#### 5.0 Particle Theory

### 5.1 One-Particle Theory

In order to understand the additional changes implied by fixing the proper-time of the source for all observers, we need only consider the question of particle dynamics. Since our motivation is quantum theory, any change of variables must be canonical. (We focus

on the X-frame equation, but the same results can also be derived for the X′ -frame.) In the conventional formulation of quantum theory, the Hamiltonian H is the generator of observer proper-time translations. We now seek to identify the Hamiltonian K which will generate source proper-time translations. To see how this may be done, let W be any classical observable so that the Poisson bracket defines Hamilton's equations in the X frame by: (here, H = p c <sup>2</sup>p<sup>2</sup> + m2c 4)

$$\frac{dW}{dt} = \frac{\partial H}{\partial \mathbf{p}} \frac{\partial W}{\partial \mathbf{x}} - \frac{\partial H}{\partial \mathbf{x}} \frac{\partial W}{\partial \mathbf{p}} = \{H, W\}.$$
 (5.1)

Now use the fact that the Hamiltonian for a free particle of mass m can be represented as H = mc<sup>2</sup>γ(w), so that γ(w) = H/mc<sup>2</sup> . This implies that

$$d\tau = (mc^2/H) dt.$$

The time evolution of the functional W is given by the chain rule:

$$\frac{dW}{d\tau} = \frac{dt}{d\tau}\frac{dW}{dt} = \frac{H}{mc^2}\left\{H, W\right\}. \tag{5.2}$$

The energy functional K conjugate to the proper-time τ must satisfy {K, W} = (H/mc<sup>2</sup> ){H, W}. The direct solution is obtained by rewriting the Poisson bracket relation in (5.2) as

$$\frac{dW}{d\tau} = \left[ \frac{H}{mc^2} \frac{\partial H}{\partial \mathbf{p}} \right] \frac{\partial W}{\partial \mathbf{x}} - \left[ \frac{H}{mc^2} \frac{\partial H}{\partial \mathbf{x}} \right] \frac{\partial W}{\partial \mathbf{p}} 
= \frac{\partial}{\partial \mathbf{p}} \left[ \frac{H^2}{2mc^2} + a \right] \frac{\partial W}{\partial \mathbf{x}} - \frac{\partial}{\partial \mathbf{x}} \left[ \frac{H^2}{2mc^2} + a \right] \frac{\partial W}{\partial \mathbf{p}}.$$
(5.3)

Now impose the condition that <sup>p</sup> = 0 <sup>⇒</sup> <sup>K</sup> <sup>=</sup> <sup>H</sup> <sup>=</sup> mc<sup>2</sup> . This gives a = a ′ = mc<sup>2</sup>/2, and

$$K = \frac{H^2}{2mc^2} + \frac{mc^2}{2} = \frac{\mathbf{p}^2}{2m} + mc^2.$$
 (5.4)

This equation was derived by Gill and Lindesay48. It looks like the nonrelativistic case but is fully relativistic and (partially) eliminates the problems associated with the square root in the conventional implementation. The most general solution is

$$K = mc^{2} + \int_{mc^{2}}^{H} (dt/d\tau)d\bar{H} = mc^{2} + \int_{mc^{2}}^{H} (\bar{H}/mc^{2})d\bar{H}.$$
 (5.5)

There are three possible solutions to this equation depending on the assumptions made.

1. If we fix the Lorentz frame, then H/mc<sup>2</sup> is constant and we get

$$K = \frac{H^2}{mc^2} = \frac{\mathbf{p}^2}{m} + mc^2. \tag{5.6}$$

This form was first derived by Gill<sup>49</sup>, and used to give a particle representation for the Klein-Gordon equation with positive probability density and with the source proper-time as an operator.

- 2. If we keep the mass fixed and allow the Lorentz frame to vary (boost), we get equation (5.4).
- 3. If we keep the momentum P = P<sup>0</sup> fixed and allow the Lorentz frame H and the mass m to vary, we get

$$K = mc^2 = \sqrt{H^2 - c^2 \mathbf{P}_0^2}. (5.7)$$

This is the appropriate Hamiltonian in the constant momentum frame. This form has received the most attention, having been used to associate the source proper-time with the (off-shell) mass operator in parametrized relativistic quantum theories. See Aparicio et al<sup>50</sup> for a recent discussion of this case. The book by Fanchi<sup>51</sup> surveys all work up to 1993 (see also Fanchi52). In all three cases, a generator can be constructed proving that they are true canonical transformations. For the first two cases, the generators are constructed in references 48 and 49 respectively. The construction of the generator for the third case was done in the seminal work of Bakamjian and Thomas <sup>15</sup> .

We plan to use equation (5.4) in our work for a number of interesting reasons. First, it is simple, directly related to the nonrelativistic case, and the quantized version is (will be) positive definite. Furthermore, since the mass is fixed, it, along with the spin, are natural choices to label the irreducible representations of the (proper-time) Poincar´e algebra describing elementary particles (see equations (5.24)-(5.32) and Wigner<sup>53</sup>). In addition, it should be noted that some of the best models for quark dynamics within nucleons "appear" to be nonrelativistic (see, for example, Strobel<sup>54</sup> and references therein).

The following theorem provides an explicit representation of the generator for the canonical change of variables for (5.4). (The result can be proved by direct computation<sup>55</sup>.)

Theorem 5.1. *If* <sup>S</sup> = (mc<sup>2</sup> <sup>−</sup> <sup>K</sup>)τ, *then* <sup>S</sup> *is the generator for the canonical change of variables from* (x, p, t, H) *to* (x, p, τ, K) *(by our X-frame observer) and*:

$$\mathbf{p} \cdot d\mathbf{x} - Hdt = \mathbf{p} \cdot d\mathbf{x} - Kd\tau + dS. \tag{5.8}$$

*It follows that the proper-time (free particle) equations will be form invariant (covariant) for all observers*.

# 5.2 Many-Particle Theory

Suppose we have a closed system of n particles with individual Hamiltonians H<sup>i</sup> and total Hamiltonian H (in the X-frame). We assume that H is of the form

$$H = \sum_{i=1}^{n} H_i. \tag{5.9}$$

If we define the effective mass M and total momentum P by

$$Mc^2 = \sqrt{H^2 - c^2 \mathbf{P}^2}, \quad \mathbf{P} = \sum_{i=1}^n \mathbf{p}_i,$$
 (5.10)

H also has the representation

$$H = \sqrt{c^2 \mathbf{P}^2 + M^2 c^4}. (5.11)$$

To construct the many-particle theory, we observe that the representation

$$d\tau = (Mc^2/H)dt (5.12)$$

does not depend on the number of particles in the system. Thus, we can uniquely define the proper-time of the system for all observers. (In the primed frame, we have a similar representation.) If we let L be the boost (generator of pure Lorentz transformations) and define the total angular momentum J by

$$\mathbf{J} = \sum_{i=1}^{n} \mathbf{x}_i \times \mathbf{p}_i, \tag{5.13}$$

we then have the following Poisson Bracket relations characteristic of the algebra for the Poincar´e group (when we use the observer proper-time):

$$\frac{d\mathbf{P}}{dt} = \{H, \mathbf{P}\} = \mathbf{0} \qquad \frac{d\mathbf{J}}{dt} = \{H, \mathbf{J}\} = \mathbf{0} \qquad \{P_i, P_j\} = 0$$
 (5.14)

$$\{J_i, P_j\} = \varepsilon_{ijk} P_k \qquad \{J_i, J_j\} = \varepsilon_{ijk} J_k \qquad \{J_i, L_j\} = \varepsilon_{ijk} L_k$$
 (5.15)

$$\frac{d\mathbf{L}}{dt} = \{H, \mathbf{L}\} = -\mathbf{P} \qquad \{P_i, L_j\} = -\delta_{ij}H/c^2, \qquad \{L_i, L_j\} = -\varepsilon_{ijk}J_k/c^2. \tag{5.16}$$

It is easy to see that M commutes with H, P, and J, and to show that M commutes with L. Constructing K as in the one-particle case, we have

$$K = \frac{H^2}{2Mc^2} + \frac{Mc^2}{2} = \frac{\mathbf{P}^2}{2M} + Mc^2.$$

Thus, we can use the same definitions for P, J, and L to obtain our new commutation relations:

$$\frac{d\mathbf{P}}{d\tau} = \{K, \mathbf{P}\} = \mathbf{0}, \quad \frac{d\mathbf{J}}{d\tau} = \{K, \mathbf{J}\} = \mathbf{0}, \quad \{P_i, P_j\} = 0, \tag{5.17}$$

$$\{J_i, P_j\} = \varepsilon_{ijk} P_k, \quad \{J_i, J_j\} = \varepsilon_{ijk} J_k, \quad \{J_i, L_j\} = \varepsilon_{ijk} L_k,$$
 (5.18)

$$\frac{d\mathbf{L}}{d\tau} = \{K, \mathbf{L}\} = \frac{-H}{Mc^2} \mathbf{P}, \quad \{P_i, L_j\} = -\delta_{ij} H/c^2, \quad \{L_i, L_j\} = -\varepsilon_{ijk} J_k/c^2. \tag{5.19}$$

It follows that, except for a constant scale change, the proper-time group is generated by the same algebra as the Lorentz group. This result is not surprising given the close relation between the two groups. It also proves our earlier statement that the form of K is fully relativistic.

Let the map from (x<sup>i</sup> , t) → (x<sup>i</sup> , τ ) be denoted by C[t, τ ], and let P(X′ , X) be the Poincar´e map from X → X′ .

Theorem 5.2 *The proper-time coordinates of the system as seen by an observer at* X *are related to those of an observer at* X′ *by the transformation*:

$$\mathbf{R}_{M}[\tau] = \mathbf{C}[t', \tau] \mathbf{P}(X', X) \mathbf{C}^{-1}[t, \tau]. \tag{5.20}$$

**Proof:** The proof follows since the diagram below is commutative.

$$X(\{\mathbf{x}_i\}, t) \longrightarrow X'(\{\mathbf{x}'_i\}, t')$$

$$\mathbf{C}^{-1}[t,\tau] \qquad \qquad \qquad \mathbf{C}[t',\tau] \tag{5.21}$$

$$X(\{\mathbf{x}_i\}, \tau) \longleftarrow X'(\{\mathbf{x}'_i\}, \tau)$$

The top diagram is the Poincaré map from  $X \to X'$ . It is important to note that this map is between the coordinates of observers. In this sense, our approach may be viewed as a direct generalization of the conventional theory. In the global case, when **U** is constant, t is related to  $\tau$  by a scale transformation so that we have a group with the same algebra as the Poincaré group (up to a constant scale), but it has an Euclidean metric! In this case, Theorem 5.2 proves that  $\mathbf{R}_M$  is in the proper-time group, formed by a similarity action on the Poincaré group by the canonical group  $\mathbf{C}_{\tau}$ . On the other hand, Theorem 5.2 is true in general. This means that in both the local and global cases (when the acceleration is nonzero) t is related to t and t via nonlocal (nonlinear) transformations. It follows that, in general, the group action is not linear, and hence is not covered by the Cartan classification.

Since K does not depend on the center-of-mass position X, it is easy to see that

$$\mathbf{U} = \frac{d\mathbf{X}}{d\tau} = \frac{\partial K}{\partial \mathbf{P}} = \frac{\mathbf{P}}{M} = \frac{1}{M} \sum_{i=1}^{n} m_i \mathbf{u}_i, \tag{5.22}$$

where  $\mathbf{u}_i = d\mathbf{x}_i/d\tau_i$ . We can now define b by

$$b = \sqrt{\mathbf{U}^2 + c^2} \Rightarrow H = Mcb. \tag{5.23}$$

Thus, equation (5.12) can also be represented as

$$d\tau = (c/b)dt. (5.24)$$

If we set  $\mathbf{v}_i = d\mathbf{x}_i/d\tau$ , an easy calculation shows that

$$\mathbf{u}_{i} = \frac{d\mathbf{x}_{i}}{d\tau_{i}} = \frac{d\tau}{d\tau_{i}} \frac{d\mathbf{x}_{i}}{d\tau} = \frac{b_{i}}{b} \mathbf{v}_{i} \Rightarrow \frac{\mathbf{u}_{i}}{b_{i}} = \frac{\mathbf{v}_{i}}{b}.$$
 (5.25)

The velocity  $\mathbf{v}_i$  is the one our observer sees when he uses the global proper-clock of the system to compute the particle velocity, while  $\mathbf{u}_i$  is the one seen when he uses the local proper clock of the particle to compute its velocity. Solving for  $\mathbf{u}_i$  and  $b_i$  in terms of  $\mathbf{v}_i$  and b, we get

$$\mathbf{u}_i = \frac{c\mathbf{v}_i}{\sqrt{b^2 - \mathbf{v}_i^2}}, b_i = \frac{cb}{\sqrt{b^2 - \mathbf{v}_i^2}} \text{ or } \frac{b_i}{b} = \frac{c}{\sqrt{b^2 - \mathbf{v}_i^2}}.$$
 (5.26)

Note that, since  $b^2 = \mathbf{U}^2 + c^2$ , if **U** is not zero, then any  $\mathbf{v}_i$  can be larger than c. On the other hand, if **U** is zero, b = c and, from the global perspective, our theory looks like the conventional one. Using (5.26), we can rewrite **U** as

$$\mathbf{U} = \frac{1}{M} \sum_{i=1}^{n} m_i \mathbf{u}_i = \frac{1}{M} \sum_{i=1}^{n} \frac{m_i c \mathbf{v}_i}{\sqrt{b^2 - \mathbf{v}_i^2}} = \frac{1}{M} \sum_{i=1}^{n} \frac{b_i m_i \mathbf{v}_i}{b} = \frac{1}{H} \sum_{i=1}^{n} H_i \mathbf{v}_i.$$
 (5.27)

It follows that the position of the center-of-mass (energy) satisfies

$$\mathbf{X} = \frac{1}{H} \sum_{i=1}^{n} H_i \mathbf{x}_i + \mathbf{Y}, \quad \frac{d\mathbf{Y}}{d\tau} = \mathbf{0}.$$
 (5.28)

It is natural to choose Y so that X is the canonical center of mass:

$$\mathbf{X} = \frac{1}{H} \sum_{i=1}^{n} H_i \mathbf{x}_i + \frac{c^2 (\mathbf{S} \times \mathbf{P})}{H(Mc^2 + H)},\tag{5.29}$$

where **S** is the (conserved) spin of the system. The important point is that  $(\mathbf{X}, \mathbf{P}, \tau, K)$  is the new set of (global) variables for the system.

**Theorem 5.3** If  $S = (mc^2 - K)\tau$ , then S is the generator for the change of variables from  $(\{\mathbf{x}_i\}, \{\mathbf{p}_i\}, \mathbf{t}, H) \to (\{\mathbf{x}_i\}, \{\mathbf{p}_i\}, \tau, K), \text{ from } (\mathbf{X}, \mathbf{P}, \mathbf{t}, H) \to (\mathbf{X}, \mathbf{P}, \tau, K), \text{ and:}$ 

$$\sum_{i=1}^{n} \mathbf{p}_i d\mathbf{x}_i - H dt = \sum_{i=1}^{n} \mathbf{p}_i d\mathbf{x}_i - K d\tau + dS,$$
(5.30)

$$\mathbf{P} \cdot d\mathbf{X} - Hdt = \mathbf{P} \cdot d\mathbf{X} - Kd\tau + dS. \tag{5.31}$$

We can now write down the transformations that fix the proper-time of the system of particles for any observer. If V is the relative velocity between two observers, we have

$$b' = \gamma(\mathbf{V}) \left[ b - \mathbf{U} \cdot \mathbf{V}/c \right], \ b = \gamma(\mathbf{V}) \left[ b' + \mathbf{U}' \cdot \mathbf{V}/c \right], \tag{5.32}$$

$$\mathbf{X}' = \gamma(\mathbf{V}) \left[ \mathbf{X}^{\dagger} - (\mathbf{V}/c)b\tau \right], \quad \mathbf{X} = \gamma(\mathbf{V}) \left[ \mathbf{X}'^{\dagger} + (\mathbf{V}/c)b'\tau \right], \quad (5.33)$$

$$\mathbf{U}' = \gamma(\mathbf{V}) \left[ \mathbf{U}^{\dagger} - (\mathbf{V}/c)b \right], \quad \mathbf{U} = \gamma(\mathbf{V}) \left[ \mathbf{U}'^{\dagger} + (\mathbf{V}/c)b' \right]. \tag{5.34}$$

As our system is closed, U is constant and τ is linearly related to t. Yet, the physical interpretation is different in the extreme if U is not zero. Furthermore, we see from equation (5.34) that, even if U is zero in one frame, it will not be zero in any other frame that is in relative motion. It is clear that τ is uniquely determined by the particles in the system and is available to all observers. Just as important is the fact that there is a very basic relationship between the global system clock and the clocks of the individual particles. In order to derive this relationship, we return to our definition of the global Hamiltonian K and let W be any observable. Then

$$\frac{dW}{d\tau} = \{K, W\} = \frac{H}{Mc^2} \{H, W\} = \frac{H}{Mc^2} \sum_{i=1}^n \{H_i, W\} 
= \frac{H}{Mc^2} \sum_{i=1}^n \frac{m_i c^2}{H_i} \left[ \frac{H_i}{m_i c^2} \{H_i, W\} \right] = \sum_{i=1}^n \frac{Hm_i}{MH_i} \{K_i, W\}$$
(5.35)

Using the (easily derived) fact that dτi/dτ = Hmi/MH<sup>i</sup> = bi/b, we get

$$\frac{dW}{d\tau} = \sum_{i=1}^{n} \frac{d\tau_i}{d\tau} \left\{ K_i, W \right\}. \tag{5.36}$$

Equation (5.36) is very important because it relates the global systems dynamics to the local systems dynamics and provides the basis for a direct approach to the quantum relativistic many-body problem using one (universal) wave function. The use of a many-times approach is not new and dates back to the early work of Dirac et al56. Our many-times approach is like that of Rohrlich and Horwitz<sup>57</sup> (see also Longhi et al58). Our approach is distinct, as is clear from (5.36) and the fact that all our times are unique and invariant for all observers.

## 5.3 Interaction (Global External)

In this section, we follow convention (in the simplest fashion) and introduce an external global interaction via minimal coupling in the free Hamiltonian. This means that we fix the position X, the momentum P, and the mass M. It is still possible for the angular momentum J to be conserved but, in general, it need not be equal to the angular momentum in the noninteracting case. Our interaction Hamiltonian becomes

$$K = \frac{\Pi^2}{2M} + Mc^2 + V(\mathbf{X}),\tag{5.37}$$

where A = A(X, τ ), V = V (X, τ ) are the vector and scalar potentials of the external field, and Π = P − (q/c)A. (In the next section, we derive an alternative equation appropriate when the cause of the external field is included in the theory to form a closed system.) Using (5.37) and Hamilton's equations, we get

$$\dot{\mathbf{X}} = \mathbf{U} = \frac{\Pi}{M}, \quad \dot{\mathbf{P}} = -\frac{\nabla \Pi^2}{2M} - \nabla V. \tag{5.38}$$

Using standard vector identities, elementary calculations give the (proper-time) Lorentz force

$$\frac{Mc}{b}\frac{d\mathbf{U}}{d\tau} = q\mathbf{E} + \frac{q}{b}\mathbf{U} \times \mathbf{B},\tag{5.39}$$

$$\mathbf{E} = -\frac{1}{b} \frac{\partial \mathbf{A}}{\partial \tau} - \nabla V, \quad \mathbf{B} = \nabla \times \mathbf{A}. \tag{5.40}$$

The fact that we can derive (a generalized form of) the Lorentz force from a (apparently) nonrelativistic Hamiltonian is well-known (see Hughes<sup>59</sup>). However, in order to see how the nonuniqueness of the Maxwell-Lorentz theory shows up here, we need only recall that  $\mathbf{W}/c = \mathbf{U}/b$  and  $(1/b)\partial/\partial\tau = (1/c)\partial/\partial t$ , so we can also write equations (5.39) and (5.40) as  $(\mathbf{W} = d\mathbf{X}/dt)$ 

$$M\frac{d\mathbf{U}}{dt} = q\mathbf{E} + \frac{q}{c}\mathbf{W} \times \mathbf{B},\tag{5.41}$$

$$\mathbf{E} = -\frac{1}{c} \frac{\partial \mathbf{A}}{\partial t} - \nabla V, \quad \mathbf{B} = \nabla \times \mathbf{A}. \tag{5.42}$$

This is the "original" force derived by Lorentz<sup>3</sup> (in 1892) and used as a part of his theory of the electrodynamics and optics of macroscopic phenomena. What is truly remarkable is the fact that the two equations (5.39) and (5.41) are mathematically equivalent, but clearly not physically equivalent, with radically different physical interpretations.

#### Global Field Theory

We can now discuss the fields of our global system of particles in a given external field. Using  $(1/c)(\partial/\partial t) = (1/b)(\partial/\partial \tau)$  (as in the one-particle case), we can write Maxwell's equations for the global system of particles as:

$$\nabla \cdot \mathbf{B} = 0, \qquad \nabla \times \mathbf{E} + \frac{1}{b} \frac{\partial \mathbf{B}}{\partial \tau} = 0,$$
 (5.43a)

$$\nabla \cdot \mathbf{E} = 4\pi \rho, \quad \nabla \times \mathbf{B} = \frac{1}{b} \left[ \frac{\partial \mathbf{E}}{\partial \tau} + 4\pi \mathbf{J} \right],$$
 (5.43b)

where  $\rho$  and  $\mathbf{J}$  represent the charge and current density of the system (as a whole) relative to its external environment. Taking the curl of the last equations of (5.43a) and (5.43b), using the standard vector identity (for any sufficiently differentiable  $\mathbf{W}$ )

$$\nabla \times (\nabla \times \mathbf{W}) = \nabla(\nabla \cdot \mathbf{W}) - \nabla^2 \mathbf{W},$$

and the first equations of (5.43a) and (5.43b), we get the corresponding global wave equations

$$\frac{1}{b} \frac{\partial}{\partial \tau} \left[ \frac{1}{b} \frac{\partial \mathbf{E}}{\partial \tau} \right] - \nabla^2 \cdot \mathbf{E} = -\nabla (4\pi\rho) - \frac{1}{b} \frac{\partial}{\partial \tau} \left[ \frac{4\pi \mathbf{J}}{b} \right],$$

$$\frac{1}{b} \frac{\partial}{\partial \tau} \left[ \frac{1}{b} \frac{\partial \mathbf{B}}{\partial \tau} \right] - \nabla^2 \cdot \mathbf{B} = \frac{1}{b} \frac{\partial}{\partial \tau} \left[ \frac{4\pi \nabla \times \mathbf{J}}{b} \right].$$
(5.44)

Computing the derivatives, these equations may also be written as

$$\frac{1}{b^{2}} \frac{\partial^{2} \mathbf{E}}{\partial \tau^{2}} - \left[ \frac{\mathbf{U}}{b^{4}} \cdot \frac{d\mathbf{U}}{\partial \tau} \right] \left[ \frac{\partial \mathbf{E}}{\partial \tau} \right] - \nabla^{2} \mathbf{E} = -\nabla (4\pi\rho) - \frac{1}{b} \frac{\partial}{\partial \tau} \left[ \frac{4\pi \mathbf{J}}{b} \right],$$

$$\frac{1}{b^{2}} \frac{\partial^{2} \mathbf{B}}{\partial \tau^{2}} - \left[ \frac{\mathbf{U}}{b^{4}} \cdot \frac{d\mathbf{U}}{\partial \tau} \right] \left[ \frac{\partial \mathbf{B}}{\partial \tau} \right] - \nabla^{2} \mathbf{B} = \frac{1}{b} \frac{\partial}{\partial \tau} \left[ \frac{4\pi \nabla \times \mathbf{J}}{b} \right].$$
(5.45)

From (5.45), we see directly that the dissipative term does not depend on the gauge. These equations imply that the field of the global system dissipates energy (radiation) throughout the enclosing domain. Since  $\mathbf{U} = (1/\mathrm{M}) \sum_{i=1}^{n} \mathbf{m_i u_i}$ , this radiation depends on the average of the (local proper) motion of all the particles in the system (e.g.,  $\mathbf{u}_i = d\mathbf{x}_i/d\tau_i$ ). This suggests that the particles live in a heat bath of radiation created by the global system's (inertial) reaction to the external field. This heat bath will fill out any domain enclosing the system of particles.

When **U** is constant,  $\dot{\mathbf{U}} = \mathbf{0}$  so that there are only velocity fields (and no radiation fields). This is necessarily the case if energy is conserved on the global level and implies the following theorem:

**Theorem 5.4** If **U** is constant then all radiation generated by internal interactions must be absorbed by the particles in the system.

The above theorem was a (required) postulate for the Wheeler-Feynman formulation. It should be noted that our formulation does not require advanced fields. As will be seen in the next Section, the individual particle interaction from the local point of view (using the particle proper-time), is of the local field type. In Section 5.5, we will see that the individual particle interaction, from the global point of view (using the global proper-time), is of the action-at-a-distance type. This confirms and refines the Wheeler-Feynman conjecture concerning the relationship between these two views.

It is clear that, in general, the above theorem is only approximately true and it is more reasonable to consider conservation of energy in a statistical sense. For example, our galaxy is clearly not a conserved system in the absolute sense, but may be considered conserved in the mean. Thus, the radiation we receive from the other galaxies is, on the average, equal to the radiation leakage from our galaxy.

### 5.4 Interaction (Internal)

In this section we assume that the system of n interacting particles can be represented via:

$$H = \sum_{i=1}^{n} H_i = H_0 + V, \quad H_i = H_{oi} + V_i,$$

$$H_{0i} = \sqrt{c^2 \pi_i^2 + m_i^2 c^4}, \quad \pi_i = \mathbf{p}_i - \frac{e_i}{c} \mathbf{A}_i,$$
(5.46)

$$H_{0} = \sum_{i=1}^{n} H_{0i}, \quad \mathbf{A}_{i} = \sum_{i \neq j} \mathbf{A}_{ji}, \quad e_{i} \mathbf{A}_{ji} = \frac{e_{i} e_{j} \left(\mathbf{w}_{j} - \mathbf{w}_{i}\right)}{2s_{ji}},$$

$$V = \sum_{i=1}^{n} V_{i}, \quad V_{i} = \sum_{i \neq j} \frac{e_{i} e_{j}}{2s_{ij}}, \quad s_{ji} = s_{ij}, \quad \frac{\partial}{\partial \mathbf{x}_{i}} (s_{ij}) = -\frac{\partial}{\partial \mathbf{x}_{j}} (s_{ij}).$$

$$(5.47)$$

Since we have specified the internal interactions, it is not a priori clear that the system is closed. Under the stated conditions, the following results can be proven by direct computation.

**Lemma 5.1** Set 
$$\mathbf{P} = \sum_{i=1}^{n} \mathbf{p}_i$$
,  $\Pi = \sum_{i=1}^{n} \pi_i$ , then  $\Pi = \mathbf{P}$ .

**Theorem 5.5** 
$$\{H, \mathbf{P}\} = 0$$
,  $\{H, V\} = 0$ ,  $\{\mathbf{P}, V\} = 0$ .

It follows that, as in Section 5.2, we can define the total effective mass M by M c<sup>2</sup> = p H<sup>2</sup> − c <sup>2</sup>P2, so that H = √ c <sup>2</sup>P<sup>2</sup> + M2c 4.

Lemma 5.2 {H, M} = 0, {P, M} = 0.

Using the above results, it now follows that the set {H<sup>i</sup> | 1 ≤ i ≤ n}, forms a closed system satisfying all the conditions of Section 5.2.

# 5.5 Particle Interaction (Local View)

We are now ready to investigate the nature of the dynamics of the ith-particle (say) caused by the action of the other particles on it. Since there are two possible clocks, τ and τ<sup>i</sup> , there are two different views, or answers, to our question. Let W<sup>i</sup> be any observable of the ith-particle, then

$$\frac{dW_{i}}{d\tau} = \{K, W_{i}\} = \sum_{j=1}^{n} \frac{\partial K}{\partial \mathbf{p}_{j}} \frac{\partial W_{i}}{\partial \mathbf{x}_{j}} - \frac{\partial K}{\partial \mathbf{x}_{j}} \frac{\partial W_{i}}{\partial \mathbf{p}_{j}}, 
\frac{dW_{i}}{d\tau_{i}} = \{K_{i}, W_{i}\} = \frac{\partial K_{i}}{\partial \mathbf{p}_{i}} \frac{\partial W_{i}}{\partial \mathbf{x}_{i}} - \frac{\partial K_{i}}{\partial \mathbf{x}_{i}} \frac{\partial W_{i}}{\partial \mathbf{p}_{i}},$$
(5.48)

$$K = \frac{H^2}{2Mc^2} + \frac{Mc^2}{2}, \quad K_i = \frac{H_i^2}{2m_ic^2} + \frac{m_ic^2}{2}.$$
 (5.49)

The equations of motion can be computed rather easily in the second case. The Hamiltonian has an explicit representation as

$$K_{i} = \frac{\pi_{i}^{2}}{2m_{i}} + m_{i}c^{2} + \frac{V_{i}^{2}}{2m_{i}c^{2}} + \frac{H_{i0}V_{i}}{m_{i}c^{2}},$$

$$\Rightarrow \frac{d\mathbf{x}_{i}}{d\tau_{i}} = \frac{\partial K_{i}}{\partial \mathbf{p}_{i}} = \frac{\pi_{i}}{m_{i}} \left(\frac{H_{i}}{H_{i0}}\right),$$
(5.50)

$$\frac{d\mathbf{p}_i}{d\tau_i} = -\frac{\partial K_i}{\partial \mathbf{x}_i} = -\frac{\nabla_i \pi_i^2}{2m_i} \left(\frac{H_i}{H_{i0}}\right) - \nabla_i V_i \left(\frac{H_i}{m_i c^2}\right). \tag{5.51}$$

Using (Hi/Hi<sup>0</sup>)∇iπ 2 <sup>i</sup> = −2(ei/c) [(u<sup>i</sup> · ∇i)A<sup>i</sup> + u<sup>i</sup> × (∇<sup>i</sup> × Ai)] , B<sup>i</sup> = (∇<sup>i</sup> × Ai), and  $(H_i/m_ic^2) = (b_i/c)$ , we have

$$\frac{d\mathbf{p}_i}{d\tau_i} = \frac{e_i}{c} \left[ (\mathbf{u}_i \cdot \nabla_i) \mathbf{A}_i + \mathbf{u}_i \times \mathbf{B}_i \right] - \frac{b_i}{c} \nabla_i V_i.$$
 (5.52)

Finally, using  $(\mathbf{u}_i \cdot \nabla_i) \mathbf{A}_i = (d\mathbf{A}_i/d\tau_i) - (\partial \mathbf{A}_i/\partial \tau_i)$ , and  $V_i = e_i \Phi_i$ , we have

$$\frac{c}{b_i}\frac{d\pi_i}{d\tau_i} = e_i \mathbf{E}_i + \frac{e_i}{b_i} \left( \mathbf{u}_i \times \mathbf{B}_i \right), \tag{5.53}$$

$$\mathbf{E}_{i} = -\frac{1}{b_{i}} \frac{\partial \mathbf{A}_{i}}{\partial \tau_{i}} - \nabla_{i} \Phi_{i}. \tag{5.54}$$

We call this the local view since it gives information about the action of the external field on the particle but provides no information about the action of the particle on the source of the external force. Equation (5.53) is of the same form as (5.39), so if we use  $(1/b_i)(\partial/\partial \tau_i) = (1/c)(\partial/\partial t)$  and  $(\mathbf{u}_i/b_i) = (\mathbf{w}_i/c)$ , we have

$$\frac{d\pi_i}{dt} = e_i \mathbf{E}_i + \frac{e_i}{c} \left( \mathbf{w}_i \times \mathbf{B}_i \right), \tag{5.55}$$

$$\mathbf{E}_{i} = -\frac{1}{c} \frac{\partial \mathbf{A}_{i}}{\partial t} - \nabla_{i} \Phi_{i}. \tag{5.56}$$

This is the same result we found in Section 5.3 when we used minimal coupling directly for the global case. For later reference we return to equation (5.50), solve for  $\pi_i$ , and differentiate, to get

$$\dot{\boldsymbol{\pi}}_i = \bar{m}_i \dot{\mathbf{u}}_i - \bar{m}_i \mathbf{u}_i \left[ \frac{(\mathbf{u}_i \cdot \nabla_i) V_i}{H_i} \right], \quad \bar{m}_i = m_i \left( 1 - \frac{V_i}{H_i} \right). \tag{5.57}$$

Putting this term in (5.53), and taking the dot product, we have

$$(\mathbf{u}_i \cdot \dot{\mathbf{u}}_i) = \frac{1}{2} \frac{d}{d\tau_i} \|\mathbf{u}_i\|^2 = \|\mathbf{u}_i\|^2 \left[ \frac{(\mathbf{u}_i \cdot \nabla_i) V_i}{H_i} \right] + \frac{e_i}{\hat{m}_i} (\mathbf{u}_i \cdot \mathbf{E}_i,).$$

$$\hat{m}_i = \frac{c}{b_i} \bar{m}_i = m_i \frac{c}{b_i} \left[ 1 - \frac{V_i}{H_i} \right].$$

$$(5.58)$$

### 5.6 Particle Interaction (Global View)

Let us now see what changes occur when we focus on the motion of the same particle as seen from the global point of view. In this case, we have

$$\frac{d\mathbf{x}_i}{d\tau} = \mathbf{v}_i = \frac{\partial K}{\partial \mathbf{p}_i} = \left(\frac{H}{M}\right) \frac{\pi_i}{H_{i0}},\tag{5.59}$$

$$\frac{d\mathbf{p}_i}{d\tau} = -\frac{\partial K}{\partial \mathbf{x}_i} = -\left(\frac{H}{Mc^2}\right) \sum_{k=1}^n \left[ \frac{c^2 \nabla_i \pi_k^2}{H_{k0}} - \nabla_i V_k \right]. \tag{5.60}$$

Using standard calculations as in the local view, and  $(H/Mc^2) = (b/c)$ , we have

$$\frac{d\mathbf{p}_i}{d\tau} = \sum_{k=1}^n \left\{ \frac{e_k}{c} \left[ (\mathbf{v}_k \cdot \nabla_i) \mathbf{A}_k + \mathbf{v}_k \times (\nabla_i \times \mathbf{A}_k) \right] - \frac{b}{c} \nabla_i V_k \right\}.$$
 (5.61)

Now use

$$(\mathbf{v}_i \cdot \nabla_i) \mathbf{A}_i = (d\mathbf{A}_i/d\tau) - (\partial \mathbf{A}_i/\partial \tau)$$

to get

$$\frac{d\mathbf{p}_{i}}{d\tau} - \frac{e_{i}}{c} \frac{d\mathbf{A}_{i}}{d\tau} = \frac{e_{i}}{c} \left[ \mathbf{v}_{i} \times \mathbf{B}_{i} \right] - \frac{e_{i}}{c} \frac{\partial \mathbf{A}_{i}}{\partial \tau} - \frac{b}{c} \nabla_{i} V_{i} 
+ \sum_{k \neq i}^{n} \left\{ \frac{e_{k}}{c} \left[ (\mathbf{v}_{k} \cdot \nabla_{i}) \mathbf{A}_{k} + \mathbf{v}_{k} \times (\nabla_{i} \times \mathbf{A}_{k}) \right] - \frac{b}{c} \nabla_{i} V_{k} \right\}.$$
(5.62)

From, (5.46) and (5.47) we see that  $(\mathbf{v}_k \cdot \nabla_i) \mathbf{A}_k = -(\mathbf{v}_k \cdot \nabla_k) \mathbf{A}_{ik}$ , etc, so we may write (5.62) in the form (using  $\mathbf{E}_i = -(1/b)(\partial \mathbf{A}_i/\partial \tau) - \nabla_i \Phi_i, \mathbf{B}_i = \mathbf{v}_i \times \mathbf{A}_i$ )

$$\frac{c}{b}\frac{d\pi_{i}}{d\tau} = e_{i}\mathbf{E}_{i} + \frac{e_{i}}{b}\left[\mathbf{v}_{i} \times \mathbf{B}_{i}\right] 
-\sum_{k \neq i}^{n} \left\{ \frac{e_{k}}{b}\left[\left(\mathbf{v}_{k} \cdot \nabla_{k}\right)\mathbf{A}_{ik} + \mathbf{v}_{k} \times \left(\nabla_{k} \times \mathbf{A}_{ik}\right)\right] - e_{k}\nabla_{k}\Phi_{ik} \right\}.$$
(5.63)

If we now set  $(\mathbf{v}_k \cdot \nabla_k) \mathbf{A}_{ik} = (d\mathbf{A}_{ik}/d\tau) - (\partial \mathbf{A}_{ik}/\partial\tau), \quad \mathbf{B}_{ik} = \nabla_k \times \mathbf{A}_{ik}, \quad \mathbf{E}_{ik} = -(1/b)(\partial \mathbf{A}_{ik}/\partial\tau) - \nabla_k \Phi_{ik}, \quad \text{and} \quad \mathbf{F}_{ik} = e_k \mathbf{E}_{ik} + (e_k/b)\mathbf{v}_k \times \mathbf{B}_{ik}, \text{ we have}$ 

$$\frac{c}{b}\frac{d\pi_i}{d\tau} = \mathbf{F}_i - \sum_{k \neq i}^n \left\{ \mathbf{F}_{ik} + \frac{e_k}{b} \frac{d\mathbf{A}_{ik}}{d\tau} \right\}. \tag{5.64}$$

If we use  $\pi_i = \bar{m}_i \mathbf{u}_i$ ,  $\bar{m}_i = m_i [1 - V_i/H_i]$ ,  $\mathbf{u}_i = [c \mathbf{v}_i/(b^2 - \mathbf{v}_i^2)^{1/2}]$ , we get (b is constant)

$$\frac{d}{d\tau} \left( \frac{\tilde{m}_i \mathbf{v}_i}{\sqrt{1 - (\mathbf{v}_i^2 / b^2)}} \right) = \mathbf{F}_i - \sum_{k \neq i}^n \left\{ \mathbf{F}_{ik} + \frac{e_k}{b} \frac{d\mathbf{A}_{ik}}{d\tau} \right\}, \tilde{m}_i = \left( \frac{c}{b} \right)^2 \bar{m}_i.$$
(5.65)

In order to interpret equation (5.65), we return to equation (5.54) and use the fact that  $(1/b_i)(\partial/\partial\tau_i) = (1/b)(\partial/\partial\tau)$  and  $(\mathbf{u}_i/b_i) = (\mathbf{v}_i/b)$  to get

$$-\frac{1}{b_i}\frac{\partial \mathbf{A}_i}{\partial \tau_i} - \nabla_i \Phi_i = -\frac{1}{b}\frac{\partial \mathbf{A}_i}{\partial \tau} - \nabla_i \Phi_i, \ \frac{e_i}{b_i} \left( \mathbf{u}_i \times \mathbf{B}_i \right) = \frac{e_i}{b} \left( \mathbf{v}_i \times \mathbf{B}_i \right).$$

This means that our force  $\mathbf{F}_i$  in (5.65) is identical to the right-hand side of (5.53) (the local Lorentz force). Equation (5.65) is our replacement for the Lorentz-Dirac equation. The second term on the right-hand side is the necessary dissipative term required to satisfy Newton's third law, and represents the action of the i-th particle on all the other particles in the system. It is important to note that this equation contains no third-order derivatives, so that it will satisfy the standard conditions for existence and uniqueness of solutions for initial value problems. It will not contain runaway solutions, nor advanced actions, etc. Furthermore, the equation does not depend on the structure of the particles in the system.

We now see that the global view of particle interactions is a pure action-at a- distance theory while from the local point of view particle interactions are mediated by the fields (a field theory).

For future reference, we assume that the global system is interacting with an external force, so that  $\dot{\mathbf{U}}$  is not zero. If we differentiate the left-hand side of (5.65), we get (using  $(\beta_i^2 = \mathbf{v}_i^2/b^2)$ ),

$$\frac{c}{b}\frac{d\pi_{i}}{d\tau} = \frac{\tilde{m}_{i}\dot{\mathbf{v}}_{i}}{\left[1 - \beta_{i}^{2}\right]^{1/2}} + \frac{\tilde{m}_{i}\mathbf{v}_{i}\left[\mathbf{v}_{i}\cdot\dot{\mathbf{v}}_{i} - \mathbf{U}\cdot\dot{\mathbf{U}}\right]}{b^{2}\left[1 - \beta_{i}^{2}\right]^{3/2}} - \frac{\tilde{m}_{i}\mathbf{v}_{i}}{\left[1 - \beta_{i}^{2}\right]^{1/2}}\frac{d}{d\tau}\left[\ln\left(1 - \frac{V_{i}}{H_{i}}\right)\right], \tilde{m}_{i} = m_{i}\frac{c^{2}}{b^{2}}\left(1 - \frac{V_{i}}{H_{i}}\right).$$
(5.66)

Taking the dot product with  $\mathbf{v}_i$ , we obtain the effective power transfer for the i-th particle

$$\frac{\tilde{m}_{i}}{2\left[1-\beta_{i}^{2}\right]^{3/2}} \frac{d\left\|\mathbf{v}_{i}\right\|^{2}}{d\tau} - \frac{\tilde{m}_{i}\left\|\mathbf{v}_{i}\right\|^{2}\left[\mathbf{U}\cdot\dot{\mathbf{U}}\right]}{b^{2}\left[1-\beta_{i}^{2}\right]^{3/2}} - \frac{\tilde{m}_{i}\left\|\mathbf{v}_{i}\right\|^{2}}{\left[1-\beta_{i}^{2}\right]^{1/2}} \frac{d}{d\tau} \left[\ln\left(1-\frac{V_{i}}{H_{i}}\right)\right]$$

$$= \mathbf{v}_{i}\cdot\mathbf{F}_{i} - \sum_{k\neq i}^{n} \left\{\mathbf{v}_{i}\cdot\mathbf{F}_{ik} + \frac{e_{k}}{b}\left(\mathbf{v}_{i}\cdot\frac{d\mathbf{A}_{ik}}{d\tau}\right)\right\}.$$
(5.67)

If U = 0, (5.64) and (5.67) become  $(\beta_i^2 = \mathbf{v}_i^2/c^2, \tilde{m}_i = \bar{m}_i, \text{ and } \tau = t)$ 

$$\frac{d}{dt} \left( \frac{\bar{m}_i \mathbf{v}_i}{\sqrt{1 - \beta_i^2}} \right) = \mathbf{F}_i - \sum_{k \neq i}^n \left\{ \mathbf{F}_{ik} + \frac{e_k}{c} \frac{d\mathbf{A}_{ik}}{dt} \right\}, \tag{5.68}$$

$$\frac{\bar{m}_{i}}{2\left[1-\beta_{i}^{2}\right]^{3/2}} \frac{d\left\|\mathbf{v}_{i}\right\|^{2}}{dt} - \frac{\bar{m}_{i}\left\|\mathbf{v}_{i}\right\|^{2}}{\left[1-\beta_{i}^{2}\right]^{1/2}} \frac{d}{dt} \left[\ln\left(1-\frac{V_{i}}{H_{i}}\right)\right]$$

$$= \mathbf{v}_{i} \cdot \mathbf{F}_{i} - \sum_{k \neq i}^{n} \left\{\mathbf{v}_{i} \cdot \mathbf{F}_{ik} + \frac{e_{k}}{c} \left(\mathbf{v}_{i} \cdot \frac{d\mathbf{A}_{ik}}{dt}\right)\right\}.$$
(5.69)

It follows that, even when the global system is at rest in the frame of the observer, our theory is distinct. In closing this section we note that summing equation (5.65) or (5.68) on i gives zero as expected, reflecting conservation of the global momentum.

## 6.0 Discussion

#### 6.1 Proper-time of the Source

In this paper, we have shown that Maxwell's equations have a mathematically equivalent formulation and additional symmetry group that fixes the proper-time of the source for all observers. The new group is closely related to the Lorentz group and, in fact, at the local level, is a nonlinear and nonlocal representation. We have constructed a dual theory using the proper-time of the source and have shown that it is covariant with respect to this group. However, the speed of light now depends on the motion of the source and the new

group replaces time transformations between observers by transformations of the velocity of light with respect to the source for different observers. This implies that the speed of light can be greater than its value in any fixed inertial frame. In the new formulation, the second postulate of the special relativity is only true when the source is in the rest frame of the observer. We have further shown that, for any closed system of particles, there is a global inertial frame and unique (invariant) global proper-clock (for each observer) from which to observe the system. In this case, the corresponding group differs from the Lorentz group by a scale transformation. This global proper-clock is intrinsically related to the proper-clocks of the individual particles in the system and provides a unique definition of simultaneity for all events associated with the system. Hence, at the global level, we can always choose a unique observer-independent measure of time for the study of physical systems. One important consequence of this result can be stated as a theorem.

Theorem 6.1 *Suppose that the observable universe is representable in the sense that the observed ratio of mass to total energy is constant and independent of our observed portion of the universe. Then the universe has a unique clock that is available to all observers*.

The above assumptions are equivalent to the homogeneity and isotropy of the energy and mass density of the universe.

The use of a global variable without attaching physical meaning to it dates back to the early work of Tetrode and Fock (for a review, see Fanchi<sup>51</sup>). However, starting in the 1970's, Horwitz and Piron<sup>60</sup> and later Fanchi<sup>51</sup> began to suggest the use of a special clock for global systems which they called the historical time. They predicted that such a variable should exist as a real physical parameter and Fanchi <sup>52</sup> suggested experiments to

detect this clock. In our approach we treat the transformation from observer proper-time to global system proper-time as a canonical (contact) transformation on extended phase space. This approach allows us to identify the canonical Hamiltonian and the associated Lie algebra (Poisson) bracket. Hence, we suggest that this global proper-time is the one sought by the above researchers. From an operational point of view, all observers can identify the time according to this (global) clock by recording the time on their clock, use the experimentally determine value for the velocity W, of the center of mass of the system, and then use equation (1.3a).

Rohrlich<sup>61</sup> has recently conducted a very interesting study of the classical self-force for the dynamics of finite-sized particles with both electromagnetic and gravitational selfinteractions (using the Lorentz-Dirac equation). He posits his model as a replacement for the point-particle model which is beyond the validity of the classical theory. His approximations neglect the nonlinear terms in the derivatives of the acceleration and leads to more reasonable equations of motion, but violates time-reversal invariance. This suggests that a successful classical theory which does not require the point-particle concept may help to explain time-reversal noninvariance at the macro-level.

As noted earlier, the proper-time theory does not depend on the size, structure, or geometry of the charge distribution. Furthermore, the global fields of any system of radiating particles in a closed domain will quickly leak radiation into every part of the domain. Since the field equations carry intrinsic information about the velocity and acceleration of each particle at the moment of dissipation, any observer will only receive information about the past behavior of the particles in the system. Since the observed radiation is an average over all the particles, this provides an explanation for the arrow of time as a statistical effect as suggested by Einstein. Also, since we only use the retarded solutions of Maxwell's equations, we may follow the suggestion of Feynman <sup>32</sup> and St¨uckelberg<sup>62</sup> and treat antimatter as matter with its proper-time reversed.

The above approach also provides us with a simple answer for questions about conservation laws during the big bang. If we assume that the big bang created two separate universes, one with matter (moving forward in proper-time), and one with antimatter (moving backward in proper-time). Then all global (physical) quantities in our universe will be conserved while providing us with a nice explanation for the lack of large concentrations of antimatter in our universe.

# 6.2 Equivalent Theories and Convention

It is no doubt more unsettling to many that the two theories could be mathematically equivalent but not physically equivalent. It is more natural to expect that two mathematically equivalent theories would also be physically equivalent, and there are a number of historical examples to support such expectations; the Lagrange-Hamiltonian formulation of classical mechanics, the Heisenberg-Schr¨odinger formulation of quantum mechanics, and the Feynman-Schwinger-Tomonaga formulation of quantum electrodynamics. In the first case, both formulations have proved equally valuable depending on the purpose. However, the latter two cases raise interesting questions.

After Feynman constructed a path integral formulation of quantum mechanics, it was shown to physically include the Heisenberg-Schrodinger formulation. However, it has never been shown to be mathematically equivalent since there are well-known serious

foundational problems with the mathematical notion of a path integral for quantum theory. On the other hand, it has not been shown that the Heisenberg-Schr¨odinger formulation is physically equivalent to the Feynman path integral approach. (There are theories where the path integral approach is easy, while the other two approaches are difficult to construct.)

In order to prove that the Feynman formulation of QED was physically equivalent to the Schwinger-Tomonaga formulation, Dyson<sup>63</sup> assumed that time had the additional property of an index which kept track of the time an operator operates (time-ordering). This represents a new physical input to theory formulation and has only recently received any mathematical attention<sup>64</sup>,65,<sup>66</sup>. Thus, mathematical equivalence has not been shown, and although some progress has been made, we are far from a solution.

In our opinion, the Feynman path integral approach is physically more general than that of Heisenberg and Schr¨odinger, and his formulation of QED is physically more general than that of Schwinger and Tomonaga. In both cases, he introduces new concepts that make it physically easier to think about and solve problems. What Feynman did was to show that it is still possible to formulate theories which more closely represent the way the world appears to us in our consciousness.

It was Poincar´e<sup>67</sup> who first noticed that some hypotheses (assumptions), which are made for theory construction, arise because of empirical data, while others occur because they are convenient. The convenient hypotheses are generally imposed by the mathematical structures we use to represent physical theories. These hypotheses are called conventions by Poincar´e in order to point out the fact that different conventions could lead to different theories which would be mathematically equivalent. He was not sure that the theories would be physically different, but he seems to have left open that possibility. The work of this paper shows that different conventions can lead to different physical theories. Since all inerital reference frames are equivalent, the one chosen by any observer is a convention. If we seek simplicity, we can all attach our frames to the MBR and use the proper-time of the universe for our global clock. In this case, we could satisfy the two postulates of the special theory, while the field and particle equations of any system would be invariant under the action of the Lorentz group (for all observers).

# 6.3 Velocity of Light

The price paid for the results of this paper will certainly seem high to many. We have rejected the third postulate of Minkowski that time be put on an equal footing with position and made a coordinate for four-geometry. We have also rejected the assumption (convention) that the observer proper-time be used to define the dynamics of an observed system. Thus, in our approach, the time is a (intrinsic) dynamical variable which must be determined by experiment along with other properties (of the observed system). This leads to a new interpretive framework in which the second postulate is only true when the source is at rest in the frame of the observer. Thus, we have reduced the observer reference frame to the prerelativistic three-geometry of Euclidean space. The observer's clock is now a part of the measuring equipment which is used to determined the proper-time of the source.

The proper-time formulation has an obvious disadvantage since, it is generally believed that, all the available experimental evidence supports the second postulate of special relativity (that the velocity of light is constant). Einstein<sup>68</sup> pointed out in a footnote to his second paper: "The principle of the constancy of the velocity of light is of course contained in Maxwell's equations." What he meant by this was that the second postulate follows from the fact that the constant c in Maxwell's equations is an invariant for all (inertial) observers. Since that time, many experiments have been done to verify that assumption. However, in 1965, Fox<sup>69</sup> wrote a very important paper which reviewed the evidence for constant c and against the emission theory of Ritz44. His conclusion was that all previous experiments were flawed for a number of reasons. In many cases, analysis of the experimental data failed to take into account the (now well-known) extinction theorem of Ewald and Oseen (see Jackson<sup>2</sup> ). The only data found that firmly supported the second postulate came from experiments on the lifetime of fast mesons and the velocity of γ rays and light from moving sources. In his conclusion Fox states that " · · · Unless something has been overlooked, these seem to be the only pieces of experimental evidence we have. This is surprising in light of the long history and importance of the problem." These "pieces of experimental evidence" have another interpretation in the proper-time theory. As noted in Section 1, the lifetime of fast mesons is the fixed value measured when they are at rest while their velocity is now computed using the proper-time of the meson which is derived from the experiment. The same interpretation applies to γ rays and light from moving sources. Thus, the same experiments that support c as constant when we assume that the observer proper-time should be used to formulate the theory also supports the result that the speed of light depends on the motion of the source when we assume that the source proper-time should be used to formulate the theory.

## 6.4 Photon Mass

Work on the question of photon mass has focused on the addition of a mass term to the

Lagrangian density for Maxwell's equations and generally leads to the Proca equation ( see Bargmann and Wigner70). Early work in this direction can be traced back from the paper of Schr¨odinger and Bass71. As in our approach, the speed of light is no longer constant in all reference frames. In this case, the fields are distorted by the mass term and experiments of Goldhaber and Nieto<sup>72</sup> use geomagnetic data to set an upper bound of 3 <sup>×</sup> <sup>10</sup>−<sup>24</sup> GeV for the mass term (see Jackiw<sup>73</sup> ). This approach causes gauge problems, and has not found favor at the classical level. The proper-time theory is fully gauge invariant and the (photon) mass is dynamical, appearing only during acceleration of the source.

It should be recalled that Maxwell's equations are (spin 1) relativistic wave equations (see Akhiezer and Berestetskii<sup>74</sup>). On the other hand, the experiments of Pound and Snider<sup>75</sup> show directly that photons have an apparent weight (as one would expect of any material object). These experiments do not depend on either the special or general theory of relativity and are not directly dependent on frequency or wavelength measurements. The existence of a small mass for the photon has important implications for QED. It is wellknown that a small photon mass can eliminate the infrared catastrophe (see Feynman<sup>76</sup>).

# Acknowledments

Work for this paper was begun while the first author was supported as a member of the School of Mathematics in the Institute for Advanced Study, Princeton, N. J., and completed during a visting appointment in the physics department at the University of Michigan. The authors would like to acknowledge important discussions, comments and encouragement from Professors G. Wienreich and H. Winful of the University of Michigan, and Professor Horwitz from the University of Tel Aviv, Isreal. We would like to give special thanks to Professor M. Wegener of the University of Aarhus, Denmark for an introduction to the work of Poincar´e.

## Appendix

In this appendix, we outline the derivation of (3.54) from the angular distribution (3.52) by taking the limit as r → ∞ after integrating over a sphere of radius r. The integrations over the azimuthal angle φ are easily done. Then, for the integrations over the polar angle θ, it is convenient to make the change of variable µ = cos θ and for a = 2, 3, . . ., b = 0, 1, 2 . . ., define the following sequence of integrals:

$$I_{a,b} \equiv \int_{-1}^{1} (1 - \beta \mu)^{-a} \mu^{b} d\mu. \tag{A1}$$

We then obtain from (3.52) that

$$\lim_{r \to \infty} \iint -\frac{dU}{dt}(\Omega) d\Omega = \frac{bq^2 |\bar{\mathbf{a}}|^2}{\bar{b}^4} \left\{ \left( 1 - \frac{1}{2} \sin^2 \alpha \right) I_{4,0} + \left( \frac{1}{2} \sin^2 \alpha - \cos^2 \alpha \right) I_{4,2} - 2\beta \left[ \beta \cos^2 \alpha \left( I_{5,0} - I_{5,2} \right) + \left( -\cos^2 \alpha + \frac{1}{2} \sin^2 \alpha \right) \left( I_{5,1} - I_{5,3} \right) \right] + \beta^2 \left[ \left( \beta^2 \cos^2 \alpha + \frac{1}{2} \sin^2 \alpha \right) \left( I_{6,0} - I_{6,2} \right) + 2\beta \cos^2 \alpha \left( I_{6,3} - I_{6,1} \right) + \left( \cos^2 \alpha - \frac{1}{2} \sin^2 \alpha \right) \left( I_{6,2} - I_{6,4} \right) \right] \right\}.$$
(A2)

Relations among the integrals (A1) for different integer values of a and b are easily obtained by integration by parts:

$$I_{a,b} = \frac{1}{\beta (a-1)} \left[ (1-\beta)^{-(a-1)} - (-1)^b (1+\beta) \right] - \frac{b}{\beta (a-1)} I_{a-1,b-1}, \tag{A3}$$

for  $a \ge 2$ ,  $b \ge 1$ ; and for b = 0, only the first term contributes:

$$I_{a,0} = \frac{1}{\beta (a-1)} \left[ (1-\beta)^{-(a-1)} - (-1)^b (1+\beta) \right], \ a \ge 2.$$
 (A4)

We note that the differences of the integrals (A1) that occur in (A2) are of the type  $I_{a,b} - I_{a,b+2}$  for given values of a and b. For differences of this type, the term in brackets in (A3) does not contribute and we have:

$$I_{a,b} - I_{a,b+2} = \frac{1}{\beta (a-1)} \left[ -b \left( I_{a-1,b-1} - I_{a-1,b+1} \right) + 2I_{a-1,b+1} \right], \tag{A5}$$

for integer values of a and b such that  $a \ge 3$ ,  $b \ge 1$ . For b = 0 the difference term on the right hand side of (A5) is missing and we have:

$$I_{a,0} - I_{a,2} = \frac{2}{\beta (a-1)} I_{a-1,1}, \quad a \ge 3.$$
 (A6)

To use the above results to evaluate (A2), we start with differences of the form (A5) and (A6) with a = 6 and b = 1, 2. The terms which arise from the difference term on the right-hand side of (A5) combine with the terms with a = 5 which are already present in (A2). After combining the coefficients of similar terms, we can then apply the process again to the integrals (A1) with a = 5. Now we have a difference from the situation with the integrals involving a = 6 that , in addition to having differences of the form (A5) with a = 5 and b = 1 and of (A6) with a = 5, we also have the integrals I5,<sup>1</sup> and I5,<sup>0</sup> which are not differences. However, these are easily evaluated by use of (A3) (giving a term involving I<sup>4</sup>,<sup>0</sup>) and (A4), respectively. We can continue this procedure to successively lower values of a, terminating at the value a = 2. The substitution of the various values of Ia,b and elimination of cos<sup>2</sup> θ by use of the identity cos<sup>2</sup> <sup>θ</sup> = 1 <sup>−</sup> sin<sup>2</sup> θ leads to the result (3.54).

# References

- [1] R. P. Feynman, R. B. Leighton, and M. Sands, The Feynman Lectures on Physics, Vol II. Addison-Wesley, New York (1974).
- [2] J. D. Jackson, Classical Electrodynamics, John Wiley & Sons, New York (1975).
- [3] H. A. Lorentz, Archives Neerlandaises des Sciences Exactes et Naturelles, 25, 353 (1892).
- [4] H. A. Lorentz, The Theory of Electrons B. G. Teubner, Leipzig, 1906; (reprinted by Dover, New York, 1952).
- [5] A. Einstein, Ann. d. Phys. 17, 891 (1905).
- [6] D. E. Spencer and U. Y. Shama, Physics Essays 9, 476 (1996).
- [7] A. Einstein, Jahrbuch Radioaktivitat V, 422 (1907) (Berichtigungen).
- [8] H. Poincar´e, C.R. Acad. Sci. Paris 140, 1504 (1905).
- [9] H. Minkowski, Physikalische Zeitschrift 10,104 (1909).
- [10] E. Whittaker, A History of Aether and Electricity, Vol. I, Thomas Nelson and Sons, London, (1951).
- [11] M. Dresden, in Renormalization: From Lorentz to Landau (and Beyond), L. M. Brown (ed.), Springer-Verlag, New York, (1993).
- [12] M.H.L. Pryce, Proc. Roy. Soc. London A 195, 400 (1948).
- [13] P. A. M. Dirac, Rev. Mod. Phys. 21, 392 (1949).
- [14] H. Leutwyler and J. Stern, Ann. Phys. (N. Y.) 112, 94 (1978).
- [15] B. Bakamjian and L. H. Thomas, Phys. Rev. 92, 1300 (1953).

- [16] D. G. Currie, T. F. Jordan, and E. C. G. Sudarshan, Rev. Mod. Phys. 35, 350 (1963).
- [17] E. C. G. Sudarshan and N. Mukunda, Classical Dynamics: A Modern Perspective. John Wiley & Sons, New York (1974).
- [18] R. Fong and J. Sucher, J. Math. Phys. 5, 456 (1964).
- [19] A. Peres, Symposia Mathematica 12, 61 (1973).
- [20] E. P. Wigner, in Aspects of Quantum Theory, in Honor of P. A. M. Dirac's 70th Birthday, edited by A. Salam and E. P. Wigner (Cambridge Univ. Press, London, 1972).
- [21] F. Rohrlich,Classical Charged Particles: Foundations of Their Theory , Addison-Wesley, Reading, Massachusetts (1965).
- [22] S. Parrott, Relativistic Electrodynamics and Differential Geometry, Springer, New York (1987).
- [23] F. Rohrlich, Phys. Rev. D 58, 116002 (1999).
- [24] W. K. H. Panofsky and M. Phillips, Classical Electricity and Magnetism, (Second Edition), Addision-Wesley, Reading, MA. (1962).
- [25] J. A. Wheeler and R. P. Feynman, Rev. of Mod. Phys. 21, 425 (1949).
- [26] J. Schwinger, Found. Phys. 13, 2573 (1998).
- [27] M. Born and L. Infield, Proc. R. Soc. London A144, 425 (1934).
- [28] P. A. M. Dirac, Proc. R. Soc. London A167, 148 (1938).
- [29] F. Bopp, Ann. d. Phys. 42, 573 (1942).
- [30] N. Rosen, Phys. Rev. 72, 298 (1947).

- [31] B. Podolsky and P. Schwed, Rev. Mod. Phys. 20, 40 (1948).
- [32] R. P. Feynman, Phys. Rev. 74, 939 (1948).
- [33] R. Haag, Z. f. Naturf. 10a, 752 (1955).
- [34] S. Parrott and D. J. Endres, Found. Phys. 25, 441 (1995).
- [35] F. E. Low, Ann. Phys. 266, 274 (1998).
- [36] A. A. Penzias and R. W. Wilson, Ap. J. 142, 419 (1965).
- [37] P.J.E. Peebles, Principles of Physical Cosmology, Princeton University Press, London (1993).
- [38] S. Schweber, QED And The Men Who Made It , Princeton University Press, London (1994).
- [39] B. French and V. Weisskopf, Phys. Rev. 75, 1240 (1949).
- [40] N. Kroll and W. Lamb, Phys. Rev. 75, 388 (1949).
- [41] P. A. M. Dirac, Sci. Amer. 208, 45 (1963).
- [42] T. P. Gill, The Doppler Effect , Logos Press, London (1965).
- [43] G. A. Schott, Phil. Mag. 29, 49 (1915).
- [44] W. Ritz, Archives des Sciences Physiques et Naturelles 16, 209 (1908).
- [45] C. Moller, The Theory of Relativity, Oxford Clarendon Press, London (1960).
- [46] C. H. Papas, Theory of Electromagnetic Wave Propagation, Dover Press, New York (1988).
- [47] R. Courant and D. Hilbert, Methods of Mathematical Physics, vol. II , Wiley-Interscience, New York (1965).

- [48] T. L. Gill and J. Lindesay, Inter. J. Theor. Phys. 32, 2087 (1993).
- [49] T. L. Gill, Fermilab-Pub-82/60-THY.
- [50] J. P. Aparicio, F. H. Gaioli, and E. T. Garcia-Alvarez, Phys. Rev. A 51, 96 (1995).
- [51] J. R. Fanchi, Parametrized Relativistic Quantum Theory, Kluwer, Dordrecht, (1993).
- [52] J.R. Fanchi, Found. Phys. 23, 487 (1993).
- [53] E. P. Wigner, Ann. Math. 40, 149 (1939).
- [54] G. L. Strobel, Inter. J. Theor. Phys. 37, 2087 (1998).
- [55] T. L. Gill, W.W. Zachary, and J. Lindesay, Inter. J. Theor. Physics 37, 2573 (1998).
- [56] P. A. M. Dirac, V. A. Fock, and B. Podolsky, Phys. Z. Sowj. Un. 2, 6 (1932) [Reprinted in J. Schwinger, ed. Selected Papers in Quantum Electrodynamics, Dover, New York (1958)].
- [57] F. Rohrlich and L. P. Horwitz, Phys. Rev. D 24, 1528 (1981).
- [58] G. Longhi, L. Lusanna, and J. M. Pons, J. Math. Phys. 30, 1893 (1989).
- [59] R. J. Hughes, Amer. J. Phys. 60, 301 (1992).
- [60] L. P. Horwitz and C. Piron, Helv. Phys. Acta 46, 316 (1981).
- [61] F. Rohrlich, Phys. Rev. D 60, 084017 (1999).
- [62] E. C. G. St¨uckelberg, Helv. Phys. Acta 15, 23 (1942).
- [63] F. J. Dyson, Phys. Rev. D 75, 486, 1736 (1949).
- [64] G. W. Johnson and M. L. Lapidus, Mem. Amer. Math. Soc. 62, 1 (1986).
- [65] T. L. Gill and W.W. Zachary, J. Math. Phys. 28, 1459 (1987).

- [66] G. W. Johnson and M. L. Lapidus, *The Feynman Integral and Feynman's Operational Calculus*, Oxford U. Press, New York, (2000).
- [67] H. Poincar´e, *Science and Hypothesis*, Dover Press, New York, (1952).
- [68] A. Einstein, Ann. d. Phys. 18, 639 (1905).
- [69] J. G. Fox, Amer. J. Phys. 33, 1 (1965).
- [70] V. Bargmann and E. P. Wigner, Proc. Nat. Acad. Sci. 34, 211 (1948).
- [71] E. Schr¨odinger and L. Bass, Proc. R. Soc. London A232, 1 (1938).
- [72] A. Goldhaber and M. Nieto, Rev. Mod. Phys. 43, 277 (1971).
- [73] R. Jackiw, Comments Mod. Phys. 1A, 1 (1999).
- [74] A. I. Akhiezer and V. D. Berestetskii, Quantum Electrodynamics, Wiley-Interscience, New York (1965)
- [75] R. V. Pound and J. L. Snider, Phys. Rev. 140, B788 (1965).
- [76] R. P. Feynman, Quantum Electrodynamics, W. A. Benjamin, New York (1964)