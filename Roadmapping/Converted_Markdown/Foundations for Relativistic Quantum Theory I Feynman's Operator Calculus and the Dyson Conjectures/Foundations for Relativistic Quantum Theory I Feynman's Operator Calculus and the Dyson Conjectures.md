# **Foundations for Relativistic Quantum Theory I: Feynman's Operator Calculus and the Dyson Conjectures**

**Tepper L. Gill1,2,3 and W. W. Zachary1,4**

**<sup>1</sup>**Department of Electrical Engineering **<sup>2</sup>**Department of Mathematics Howard University Washington, DC 20059 E-mail: tgill@howard.edu

> **<sup>3</sup>**Department of Physics University of Michigan Ann Arbor, Mich. 48109

**<sup>4</sup>**Department of Mathematics and Statistics University of Maryland University College College Park, Maryland 20742 E-mail: wwzachary@earthlink.net

#### **Abstract**

In this paper, we provide a representation theory for the Feynman operator calculus. This allows us to solve the general initial-value problem and construct the Dyson series. We show that the series is asymptotic, thus proving Dyson's second conjecture for QED. In addition, we show that the expansion may be considered exact to any finite order by producing the remainder term. This implies that every nonperturbative solution has a perturbative expansion. Using a physical analysis of information from experiment versus that implied by our models, we reformulate our theory as a sum over paths. This allows us to relate our theory to Feynman's path integral, and to prove Dyson's first conjecture that the divergences are in part due to a violation of Heisenberg's uncertainly relations.

PACS classification codes: 02.30.Sa, 02.30.Tb, 03.65.Bz, 12.20.-m, 11.80.-m Keywords: Feynman operator calculus, Dyson conjecture, divergences in QED

### **1.0 Introduction**

Following Dirac's quantization of the electromagnetic field in 19271 , and his relativistic electron theory in 19282 , the equations for quantum electrodynamics (QED) were developed by Heisenberg and Pauli3,4 in the years 1929-30 (see Miller5 and Schweber6 ). From the beginning, when researchers attempted to use the straightforward and physically intuitive timedependent perturbation expansion to compute physical observerables, a number of divergent expressions appeared. Although it was known that the same problems also existed in classical electrodynamics, it was noted by Oppenheimer7 that there was a fundamental difference in the quantum problem as compared to the classical one. (Dirac9 had shown that, in the classical case, one could account for the problem of radiation reaction without directly dealing with the self-energy divergence by using both advanced and retarded fields and a particular limiting procedure.)

Early attempts to develop subtraction procedures for the divergent expressions were very discouraging because they depended on both the gauge and the Lorentz frame, making them appear ambiguous. Although the equations of QED were both Lorentz and gauge covariant, it was generally believed that, in a strict sense, they had no solutions expandable in powers of the charge. The thinking of the times was clearly expressed by Oppenheimer8 in his 1948 report to the Solvay Conference, "If one wishes to explore these solutions, bearing in mind that certain infinite terms will, in a later theory, no longer be infinite, one needs a covariant way of identifying these terms; and for that, not merely the field equations themselves, but the whole method of approximation and solution must at all stages preserve covariance."

The solution to the problem posed by Oppenheimer was made (independently) by Tomonaga10, Schwinger11 and Feynman12,13. (These papers may be found in Schwinger14.) Tomonaga introduced what is now known as the interaction representation and showed how the approximation process could be carried out in a covariant manner. Schwinger developed the general theory and applied it to many of the important problems. Feynman took a holistic view of physical reality in his development. He suggested that we view a physical event as occurring on a film which exposes more and more of the outcome as the film unfolds. His idea was to deal directly with the solutions to the equations describing the physical system, rather than the equations themselves. In addition to solving the problem posed by Oppenheimer, Feynman's approach led to a new perturbation series, which provided an easy, intuitive, and computationally simple method to study interacting particles while giving physical meaning to each term in his expansion.

Since Feynman's method and approach was so different, it was not clear how it related to that of Schwinger and Tomonaga. Dyson15,16, made a major contribution. Dyson realized that Feynman and Schwinger were both dealing with different versions of Heisenberg's S-matrix. He then formally introduced time-ordering and provided a unified approach by demonstrating the equivalence of the Feynman and Schwinger-Tomonaga theories. This approach also allowed him to show how the Schwinger theory could be greatly simplified and extended to all orders of the perturbation expansion. Dyson's time-ordering idea was actually obtained from discussions with Feynman, who later explored and fully developed it into his time-ordered operator calculus17.

#### **1.1 Background**

After the problem proposed by Oppenheimer was resolved, attitudes toward the renormalization program and quantum field theory could be classified into three basic groups. The first group consisted of those who were totally dissatisfied with the renormalization program. The second group considered the renormalization program an interim step and believed that the divergences were an indication of additional physics, which could not be reached by present formulations. The first two groups will not be extensively discussed in this paper. However, we can associate the names of Dirac and Landau with the first group, and Sakata and Schwinger with the second. (See Dirac18, Sakata19, Schwinger20 and also Schweber6 .)

The third group was more positive, and directed its attention towards investigating the mathematical foundations of quantum field theory with the hope of providing a more orderly approach to the renormalization program (assuming that the theory proved consistent). This direction was clearly justified since part of the problem had been consistently blamed on a mathematical issue, the perturbation expansion. Indeed, the whole renormalization program critically depended on the expansion of the Smatrix in powers of the coupling constant. This concern was further supported since attempts to use the expansion when the coupling constant was large led to meaningless results. Additional unease could be attributed to the fact that, at that time, not much was actually known about the physically important cases where one was dealing with unbounded operatorvalued functions (distributions).

Researchers working on the mathematical foundations of quantum electrodynamics, and quantum field theory adopted the name axiomatic field theory starting in the fifties. These researchers focused on trying to find out what could be learned about the existence of local relativistic quantum field theories based on certain natural assumptions which included the postulates of quantum mechanics, locality, Poincaré invariance, and a reasonable spectrum. This approach was initiated by the work of Wightman21, and Lehmann, Symanzik and Zimmermann22,23. Here, the quantized field is interpreted mathematically as an operator-valued Schwartz distribution. Explicit use of the theory of distributions was a major step, which helped to partially make the theory (mathematically) sound by smoothing out the fields locally. (The recent paper by Wightman24 provides an inspired introduction to the history of Heisenberg's early observations on the latter concept and its relationship to the divergences25.)

The axiomatic approach proved very fruitful, providing the first rigorous proofs of a number of important general results, and attracted many able researchers. The favored name today is Algebraic Quantum Field Theory. The books by Jost26, Streater and Wightman27 and Bogolubov and Shirkov28 are the classics, while more recent work can be found in Haag29. (See also the book by Bogolubov, Logunov and Todorov30, and the recent review paper by Buchholz31.)

For a number of reasons, most notably a lack of nontrivial examples, the axiomatic approach evolved in a number of directions. One major direction is called "constructive" quantum field theory. Here, one focuses on attempts to directly construct solutions of various model field theories, which either have exact (nonperturbative) solutions, or have an asymptotic perturbative expansion which can be summed to the exact solution. In this approach, instead of formulating the theory in Minkowski spacetime, one passes to imaginary time and formulates it in Euclidean space (an idea which first appeared in Dyson15). This leads to a formulation in terms of "Schwinger functions", also known as Euclidean Green's functions. The advantage of this approach is that hyperbolic equations are transformed to elliptic ones, and Gaussian kernels, for which a very rich set of analytic tools has been developed, replace Feynman kernels. The output of this enterprise is truly impressive. Constructive solutions have been obtained for a number of important models. Furthermore, this approach has given us a clearer picture of the problems associated with the rigorous construction of a relativistic quantum field theory and provided new mathematical methods. An early summary of this approach may be found in the lecture notes32, while more recent progress is contained in the lecture notes33, both edited by Velo and Wightman (see also references 41 and 6). The books by Glimm and Jaffe34 and Simon35 give a different flavor and point of departure.

Although a great deal of work has been done in constructive field theory over the last thirty years, many difficult problems still remain. For example, the appearance of difficulties with the constructive approach to polynomial types of field theories is discussed in the paper by Sokal36. He conjectured that the lj≥<sup>4</sup> <sup>4</sup> theory (lj<sup>4</sup> in four or more spacetime dimensions) is a generalized free field, where l is the coupling constant. This theory represents a self-interacting boson field. The conjecture was proven by Aizenman and Graham37 and Fröhlich38. Three years later, Gawedzki and Kupiainen39 proved that, if we change the sign of the coupling constant, the solution exists (as a tempered distribution) and the perturbation expansion is asymptotic to the solution. This state of affairs led Wightman (reference 33, pg. 1) to lament that "We do not know whether the lack of an existence theorem for solutions with the "right" sign reflects the non-existence of solutions or merely the lack of a technique to construct them." Things are further complicated by the fact that the lj<sup>4</sup> <sup>4</sup> theory has a perturbative solution! This led Gallavotti45 to suggest that constructive approaches other than the ferromagnetic lattice approximation, used by Aizenman and Graham, and Fröhlich, may be required.

The most well known method for quantum field theory calculations is perturbative renormalization theory. This approach is discussed in most standard texts on quantum field theory and has an interesting history that is best told by Wightman40. (The first book to include Dyson's reformulation of the Feynman-Schwinger-Tomonaga theory is the classic by Jauch and Rohrlich42.) Early work in the perturbative approach focused on the development of different renormalization methods with the hope of identifying those for which rigorous mathematical methods could be used. The methods generally consisted of two parts. First, the Green's functions were regularized in a relativistically and gauge invariant manner28,40,41,46 to yield well-defined tempered distributions, even on the light cone. Then appropriate counter-terms were introduced so that, in the limit, when the regularization was removed, the various divergences of the S-matrix were also removed. It was found that all renormalization procedures are equivalent up to a finite renormalization (cf. references 40, 41). Today, theories are classified as "renormalizable" or "unrenormalizable" according as the number of renormalizable constants is finite or infinite, respectively.

Some model theories in less than four spacetime dimensions considered in constructive field theory belong to a special subclass of renormalizable theories called "super renormalizable", for which the renormalization process can be carried out without using perturbation theory32-35. For these theories, the renormalized perturbation series can be shown to be Borel summable to the exact nonperturbative solution. A nice summary of these developments was given by Glimm and Jaffe34. On the other hand, constructive models of the Gross-Neveu type are renormalizable but not super renormalizable (see reference 33).

Feldman et al43 have studied the mathematical foundations of quantum electrodynamics from the perturbative point of view (see also Rosen in reference 33, pg. 201). Here, a renormalized formal power series (renormalized tree expansion) is obtained for a measure on the space of fields within the Euclidean formulation of QED. (The tree expansion method is an outgrowth of Wilson's44 renormalization group approach as distilled by Gallavotti45 and co-workers.) It is then shown that QED in four (Euclidean) dimensions is locally Borel summable. Their work is truly remarkable and represents the first (formal) proof that (Euclidean) quantum electrodynamics can be renormalized using gauge invariant counterterms. However, in general, it is a nontrivial problem to return from the Euclidean regime to Minkowski space. The return trip requires application of the Osterwalder-Schhrader reconstruction theorem (see reference 32). This theorem places conditions on the Euclidean Green's functions which guarantees analytic continuation back to the real-time vacuum expectation values. When these conditions are fulfilled, the Lehmann, Symanzik, and Zimmermann (LSZ)22,23,32 reduction formulae may then be used to obtain the S-matrix. For technical reasons, they were not able to directly apply the Osterwalder-Schhrader theorem. They could still get back to QED in Minkowski spacetime by following the methods of Hepp46 and Lowenstein and Speer47. However, nothing could be said about the convergence properties of their series.

## **1.2 Purpose**

It is clear that Dyson's use of time-ordering was the fundamental conceptual tool which allowed him to relate the Feynman and Schwinger-Tomonaga theories. This tool has now become a natural part of almost every branch of physics and is even used in parts of engineering. Its importance to the foundations of quantum field theory led Segal53 to suggest that the identification of mathematical meaning for Feynman's time-ordered operator calculus is one of the major problems. A number of investigators have attempted to solve this problem. Miranker and Weiss54 showed how the Feynman ordering process could be done formally using the theory of Banach algebras. Nelson55 used Banach algebras to developed a theory of "operants" as an alternate (formal) approach. Araki56, motivated by the work of Fujiwara, used Banach algebras to develop yet another formal approach. (Fujiwara57 had earlier suggested that the Feynman program could be implemented if one used a sheet of unit operators at every point except at time t, where the true operator should be placed.) Maslov58 used the idea of a T-product to formally order operators and developed an operational theory. Another important approach to this problem via the idea of an index may be found in the works of Johnson and Lapidus59-61, see also Johnson, Lapidus, and DeFacio62.

This paper is a part of a new investigation into the physical and mathematical foundations of relativistic quantum theory. Our overall goal is to construct a self-consistent relativistic quantum theory of particles and fields. For this paper, we have two specific objectives. Our first (and major) objective is to construct a physically simple and computationally useful representation theory for the Feynman time-ordered operator calculus.

A correct formulation and representation theory for the Feynman time-ordered operator calculus should at least have the following desirable features:

- **1.** It should provide a transparent generalization of current analytic methods without sacrificing the physically intuitive and computationally useful ideas of Feynman.
- **2.** It should provide a clear approach to some of the mathematical problems of relativistic quantum theory.
- **3.** It should explain the connection with path integrals.

In the course of his analysis, unification, and simplification of the Feynman-Schwinger-Tomonaga theory, Dyson made two important suggestions (conjectures). The first conjecture concerned the divergences in QED, while the second was concerned with the convergence of the renormalized perturbation series. In addressing the problem of divergences, Dyson conjectured that they may be due to an idealized conception of measurability resulting from the infinitely precise knowledge of the spacetime positions of particles (implied by our Hamiltonian formulation) which leads to a violation of the Heisenberg uncertainty principle. This point of view can be traced directly to the Bohr-Rosenfeld theory of measurability for field operators and, according to Schweber6 , is an outgrowth of Dyson's discussions with Oppenheimer.

In addressing the renormalized S-matrix16, Dyson suggested that it might be more reasonable to expect the expansion to be asymptotic rather than convergent and gave physical arguments to support his claim. The lack of a clear mathematical framework made it impossible to formulate and investigate his suggestions.

Schweber6 notes that Dyson made two other well-known conjectures. The "overlapping divergences" conjecture was proved by Salam48, Ward49, Mills and Yang50, and Hepp51. Dyson's conjecture that a certain Feynman integral converges, necessary for showing that the ultraviolet divergences cancel to all orders, was proved by Weinberg52.

Our second objective is to provide proofs of the above two conjectures under general conditions that should apply to any formulation of quantum field theory which does not abandon Hamiltonian generators for unitary solution operators. The proof of the first conjecture is, to some extent expected, and is a partial vindication of our belief in the consistency of quantum electrodynamics in the sense that the ultraviolet problem is caused by an effect that is basically "simple". Such a result is partly anticipated since the effect can be made to disappear via appropriate cutoffs. We also identify (special) conditions under which the renormalized perturbation series may actually converge. A proof of the above conjectures is implicit in, and is one of the major achievements of constructive field theory for the models studied. In fact, these theories verify a stronger version of the second conjecture since, as noted earlier, the renormalized perturbation series is summable to the true solution.

## **1.3 Summary**

The work in this paper is both a generalization and simplification of earlier work63-65 that is easier and requires the weakest known conditions. We construct a new representation Hilbert space and von Neumann algebra for the Feynman (time-ordered) operator calculus. In order to make the theory applicable to other areas, we develop it using semigroups of contractions and the Riemann integral. A contraction semigroup on a Hilbert space H can always be extended to a unitary group on a larger space H¢. Thus, for quantum theory we may replace the semigroups by unitary groups and assume that our space is H¢ without any loss in understanding.

The Riemann integral can be easily replaced by the operator-valued Riemann-complete integral of Henstock66 and Kurzweil67, which generalizes the Bochner and Pettis integrals (see Gill63). This integral is easier to understand (and learn) compared to the Lebesgue or Bochner integrals, and provides useful variants of the same theorems that have made those integrals so important. Furthermore, it arises from a simple (transparent) generalization of the Riemann integral that was taught in elementary calculus. Its usefulness in the construction of Feynman path integrals was first shown by Henstock68, and has been further explored in the recent book by Muldowney69.

In Section 1.4 we provide a brief review of the necessary operator theory in order to make the paper self-contained. In Section 2 we construct an infinite tensor product Hilbert space and define what we mean by time ordering. In Section 3 we construct time-ordered integrals and evolution operators and prove that they have the expected properties. In Section 4 we define what is meant by the phase "asymptotic in the sense of Poincaré" for operators, and use it to prove Dyson's second conjecture for contraction semigroups. We then discuss conditions under which the perturbation series may be expected to converge.

In Section 5 we take a photograph of a track left by an elementary particle in a bubble chamber as a prototype to conduct a physical analysis of what is actually known from experiment. This approach is used to rederive our time-ordered evolution operator as the limit of a probabilistic sum over paths. We use it to briefly discuss our theory in relationship to the Feynman path integral, and show that it provides a general and natural definition for the path integral that is independent of measure theory and the space of continuous paths.

The results from Section 5 are applied to the S-matrix expansion in Section 6 to provide a formulation and proof of Dyson's first conjecture. In particular, we show that, within our formulation, the assumption of precise time information over a particle's trajectory introduces an infinite amount of energy into the system at each point in time. We use Dyson's original notation partly for reasons of nostalgia, but also to point out what we are not able to explain within our framework. Also, since all renormalization procedures are equivalent, there is no loss.

# **1.4 Operator Theory**

In this subsection we establish notation and quote some results from operator theory used in the paper. Let H denote a separable Hilbert space over **C** (complex numbers), **B**(H) the set of bounded linear operators, and **C**(H) the set of closed densely defined linear operators on H.

**Definition 1.0** A family of bounded linear operators {, } U( ), < *t t* 0 0 £ • defined on H is a *strongly continuous semigroup (or* C0 - *semigroup)* if

1. U(0,0) = I, 2. U(t+s,0) = U(t,0)U(s,0), 3.  $\lim_{t\to 0} U(t,0)\varphi = \varphi$ ,  $\forall \varphi \in \mathcal{H}$  U(t,0) is a contraction semigroup in case  $||U(t,0)|| \le 1$ . If we replace 2 by 2'.  $U(t,\tau) = U(t,s)U(s,\tau)$ ,  $0 \le \tau \le s \le t < \infty$ , then we call  $U(t,\tau)$  a *strongly continuous evolution family*.

**Definition 1.2** A densely defined operator H is said to be *maximal dissipative* if  $Re\langle H\varphi,\varphi\rangle \leq 0$ ,  $\forall \varphi \in D(H)$ , and  $\mathcal{R}an(I-H) = \mathcal{H}$  (range of (I-H)).

The following results may be found in Goldstein<sup>70</sup> or Pazy<sup>71</sup>.

**Theorem 1.2** Let U(t,0) be a  $C_0$ -semigroup of contraction operators on  $\mathcal{H}$ .

Then

1) 
$$H\varphi = \lim_{t\to 0} \frac{\mathrm{U}(t,0)\varphi - \varphi}{t}$$
 exists for  $\varphi$  in a dense set.

2) 
$$R(z, H) = (zI - H)^{-1}$$
 exists for  $z > 0$  and  $||R(z, H)|| \le \frac{1}{z}$ .

**Theorem 1.3** Suppose H is a maximal dissipative operator. Then H generates a unique  $C_0$ -semigroup  $\{U(t,0) | 0 \le t < \infty\}$  of contraction operators on  $\mathcal{H}$ .

**Theorem 1.4** If H is densely defined with both H and  $H^*$  dissipative, then H is maximal dissipative.

# 2.0 Infinite Tensor Product von Neumann Algebras

In this section we define time-ordered operators and construct the representation space which will be used in the next section to develop our theory of time-ordered integrals and evolution operators. Much of the

material in this section was developed by von Neumann72 for other purposes, but is perfectly suited for our program. In order to see how natural our approach is, let H<sup>ƒ</sup> = ƒˆ <sup>s</sup>H*(s)* denote the infinite tensor product Hilbert space of von Neumann, where H(*s*) =H for s Œ [a,b] and ƒˆ denotes closure. If **B**(Hƒ) is the set of bounded operators on Hƒ, define **B**(H(*t*))Ã**B**(Hƒ) by

$$\mathbf{B}(\mathcal{H}(t)) = \{\mathbf{H}(t) | \mathbf{H}(t) = \hat{\otimes}_{a \ge s > t} \mathbf{I}_{s} \otimes H(t) \otimes (\otimes_{t > s \ge -a} \mathbf{I}_{s}), \forall H(t) \in \mathbf{B}(\mathcal{H})\}, \quad (2.1a)$$

where Is denotes an identity operator, and let **B#** (Hƒ) be the uniform closure of the von Neumann algebra generated by the family {**B**(H(*t*)), E *t* Œ }. If the family { ( ) | E} *Ht t* Œ is in **B**(H), then the corresponding operators { ( ) | E} **H** *t t* Œ Œ**B#** (Hƒ) commute when acting at different times: *t s* π fi

$$\mathbf{H}(t)\mathbf{H}(s) = \mathbf{H}(s)\mathbf{H}(t). \tag{2.1b}$$

**Definition 2.0** The smallest space F D<sup>ƒ</sup> Õ Hƒwhich, leaves the family { ( ) **H** *t t*| E} Œ invariant is called a Feynman-Dyson space for the family. (This is the film.)

We need the following results about operators on Hƒ.

**Theorem 2.1:** (von Neumann72 *The mapping* **T**<sup>q</sup> *t* **: B**(H) Æ**B**(H(*t*)) *is an isometric isomorphism of algebras.* (*We call* **T**<sup>q</sup> *<sup>t</sup> the time-ordering morphism.*)

**Definition 2.2** *The vector* F=ƒ*<sup>s</sup>* <sup>f</sup>*<sup>s</sup> is said to be equivalent to* Y=ƒ*<sup>s</sup>* <sup>y</sup> *<sup>s</sup> and we write* F Yª , *if and only if*

$$\sum_{s} \left| \left\langle \phi_{s}, \psi_{s} \right\rangle_{s} - 1 \right| < \infty. \tag{2.2}$$

Here,  $\langle \cdot, \cdot \rangle_s$  is the inner product on  $\mathcal{H}(s)$ , and it is understood that the sum is meaningful only if at most a countable number of terms are different from zero.

Let  $\mathcal{H}_{\Phi} = cl \bigg\{ \Psi \, \bigg| \, \Psi = \sum_{i=1}^n \Psi_i, \, \Psi_i \approx \Phi, \, n \in \mathbf{N} \bigg\}$  (closure),  $\Phi \in \mathcal{H}_{\otimes}$ , and let  $\mathbf{P}_{\Phi}$  denote the projection from  $\mathcal{H}_{\otimes}$  onto  $\mathcal{H}_{\Phi}$ . The space  $\mathcal{H}_{\Phi}$  is known as the incomplete tensor product generated by  $\Phi$ . The details on incomplete tensor product spaces as well as proofs of the next two theorems may be found in von Neumann<sup>72</sup>.

**Theorem 2.3** The relation defined above is an equivalence relation on  $\mathcal{H}_{\otimes}$  and 1) if  $\Psi$  is not equivalent to  $\Phi$ , then  $\mathcal{H}_{\Phi} \cap \mathcal{H}_{\Psi} = \{0\}$  (i.e.,  $\mathcal{H}_{\Phi} \perp \mathcal{H}_{\Psi}$ );

- 2) if  $\psi_s \neq \phi_s$  occurs for at most a finite number of s, then  $\Phi = \bigotimes \phi_s \approx \Psi = \bigotimes \psi_s$ ;
- 3) if  $T \in B^{\#}(\mathcal{H}_{\otimes})$ , then  $P_{\Phi}T = TP_{\Phi}$  so that  $P_{\Phi}T \in B^{\#}(\mathcal{H}_{\Phi})$ .

The second condition in Theorem 2.3 implies that, for each fixed  $\Phi = \underset{s}{\otimes} \phi_{s}$ , there is an uncountable number of  $\Psi = \underset{s}{\otimes} \psi_{s}$  equivalent to  $\Phi$ , while the third condition implies that every bounded linear operator on  $\mathcal{H}_{\otimes}$  restricts to a bounded linear operator on  $\mathcal{H}_{\Phi}$  for each  $\Phi$ .

We can now construct our film  $\mathcal{FD}_{\otimes}$ . Let  $\left\{e^{i} \middle| i \in N\right\}$  denote an arbitrary ordered complete orthonormal basis (c.o.b) for  $\mathcal{H}$ . For each  $t \in E, i \in N$ , let  $e^{i}_{t} = e^{i}$ ,  $E^{i} = \underset{t \in E}{\otimes} e^{i}_{t}$ , and define  $\mathcal{FD}^{i}$  to be the incomplete tensor product generated by the vector  $E^{i}$ . Setting  $\mathcal{FD}_{\otimes} = \underset{i=1}{\overset{\infty}{\oplus}} \mathcal{FD}^{i}$ , it will be clear in the next section that  $\mathcal{FD}_{\otimes}$  is (one of an infinite number of) the natural representation

space(s) for Feynman's time-ordered operator theory. It should be noted that  $\mathcal{FD}_{\otimes}$  is a nonseparable Hilbert (space) bundle over [a, b]. However, it is not hard to see that each fiber is isomorphic to  $\mathcal{H}$ .

In order to facilitate the proofs in the next section, we need an explicit basis for each  $\mathcal{FD}^i$ . To construct it, fix i and let  $f^i$  denote the set of all functions  $\{j(t)|t\in E\}$  mapping  $E\to N\cup\{0\}$  such that j(t) is zero for all but a finite number of t. Let  $I(j)=\{j(t)|t\in E\}$  denote the function j and set  $E^i_{I(j)}=\bigotimes_{t\in F}e^i_{t,j(t)}$  with  $e^i_{t,0}=e^i$ , and  $j(t)=k\Rightarrow e^i_{t,k}=e^k$ .

**Theorem 2.4** The set  $\left\{E_{I(j)}^{i}\middle|I(j)\in f^{i}\right\}$  is a (c.o.b) for each  $\mathcal{FD}^{i}$ .

For each 
$$\Phi^{i}$$
,  $\Psi^{i} \in F^{i}$ , set  $a^{i}_{I(j)} = \left\langle \Phi^{i}, E^{i}_{I(j)} \right\rangle$ ,  $b^{i}_{I(j)} = \left\langle \Psi^{i}, E^{i}_{I(j)} \right\rangle$ , so that  $\Phi^{i} = \sum_{I(j) \in F^{i}} a^{i}_{I(j)} E^{i}_{I(j)}$ ,  $\Psi^{i} = \sum_{I(j) \in F^{i}} b^{i}_{I(j)} E^{i}_{I(j)}$  and  $\left\langle \Phi^{i}, \Psi^{i} \right\rangle = \sum_{I(j) \in F^{i}} a^{i}_{I(j)} \overline{b}^{i}_{I(k)} \left\langle E^{i}_{I(j)}, E^{i}_{I(k)} \right\rangle$ . Now,  $\left\langle E^{i}_{I(j)}, E^{i}_{I(k)} \right\rangle = \prod_{I(j) \in F^{i}} \left\langle e^{i}_{t, I(j)}, e^{i}_{t, I(k)} \right\rangle = 0$ , unless  $j(t) = k(t)$ ,  $\forall t \in E$ , so that  $\left\langle \Phi^{i}, \Psi^{i} \right\rangle = \sum_{I(j) \in F^{i}} a^{i}_{I(j)} \overline{b}^{i}_{I(j)}$ .

We need the notion of an exchange operator. (Theorem 2.6 is in reference 63.)

**Definition 2.5** An exchange operator  $\mathbf{E}[t,t']$  is a linear map defined for pairs  $t, t' \in [a, b]$  such that :

- 1.  $\mathbf{E}[t,t']: \mathbf{B}(\mathcal{H}(t)) \to \mathbf{B}(\mathcal{H}(t'))$  onto,
- 2.  $\mathbf{E}[t,s]\mathbf{E}[s,t'] = \mathbf{E}[t,t'],$
- 3.  $\mathbf{E}[t,t']\mathbf{E}[t',t] = 1$ ,
- 4. if  $s \neq t,t'$ , then  $\mathbf{E}[t,t']\mathbf{H}(s) = \mathbf{H}(s) \ \forall \mathbf{H}(s) \in \mathbf{B}(\mathcal{H}(s))$ .

#### Theorem 2.6

- 1)  $\mathbf{E}[\cdot,\cdot]$  exists and is a Banach algebra isomorphism on  $\mathbf{B}^{\#}(\mathcal{H}_{\otimes})$ .
- 2)  $\mathbf{E}[s,s']\mathbf{E}[t,t'] = \mathbf{E}[t,t']\mathbf{E}[s,s']$  for distinct pairs (s,s') and (t,t') in  $\mathbf{E}$ .

#### 3.0 Time-Ordered Integrals

In this section we construct time-ordered integrals and evolution operators for a fixed family  $\{H(t)|t\in E\}\subset C(\mathcal{H})$  of generators of contraction semigroups on  $\mathcal{H}$ . We assume that, for each t, H(t) and  $H^*(t)$  are dissipative (so that the family is maximal dissipative for each t). In the following discussion we adopt the notation:

- 1). (e.o.v): "except for at most one s value";
- 2). (e.f.n.v): "except for an at most finite number of s values"; and
- 3). (a.s.c): "almost surely and the exceptional set is at most countable".

The s value referred to is in our fixed interval E.

For the given family  $\{H(t)|t \in E\} \subset C(\mathcal{H})$ , define  $\exp\{\tau \mathbf{H}(t)\}$  by

$$\exp\{\tau \mathbf{H}(t)\} = \hat{\bigotimes}_{s \in [b,t)} \mathbf{I}_s \otimes \left(\exp\{\tau \mathbf{H}(t)\}\right) \otimes \left(\bigotimes_{s \in (t,a]} \mathbf{I}_s\right), \tag{3.1}$$

and set  $\mathbf{H}_z(t) = z\mathbf{H}(t)\mathbf{R}(z,\mathbf{H}(t))$ , z > 0, where  $\mathbf{R}(z,\mathbf{H}(t)) = (z\mathbf{I}_{\otimes} - \mathbf{H}(t))^{-1}$  is the resolvent of  $\mathbf{H}(t)$ . It is known that  $H_z(t)$  generates a uniformly bounded contraction semigroup and  $\lim_{z \to \infty} H_z(t)\phi = H(t)\phi$  for  $\phi \in D(H(t))$ .

**Theorem 3.1** Suppose for each t,  $\{H(t)|t\in E\}\subset C(\mathcal{H})$  generates a strongly continuous contraction semigroup on  $\mathcal{H}$ . Then  $\mathbf{H}(t)\mathbf{H}_z(t)\Phi = \mathbf{H}_z(t)\mathbf{H}(t)\Phi$ ,  $\Phi\in D$ , (where D denotes the domain of the family  $\{\mathbf{H}(t)|t\in E\}$ ), and

- 1. The family  $\{\mathbf{H}_z(t)|t\in E\}$  generates a uniformly bounded contraction semigroup on  $\mathcal{FD}_{\otimes}$  for each t and  $\lim_{z\to\infty}\mathbf{H}_z(t)\Phi=\mathbf{H}(t)\Phi$ ,  $\Phi\in D$ .
- 2. The family  $\{\mathbf{H}(t)|t\in E\}\subset \mathbf{C}(\mathcal{H}_{\otimes})$  generates a strongly continuous contraction semigroup on  $\mathcal{FD}_{\otimes}$  (so that  $\{\mathbf{H}(t)|t\in E\}\subset \mathbf{C}(\mathcal{FD}_{\otimes})$ ).

**Proof:** The proof of 1. is standard. Note that  $\mathbf{H}_{z}(t) = z^{2}\mathbf{R}(z,\mathbf{H}(t)) - z\mathbf{I}_{\otimes}$  and  $\|\mathbf{R}(z,\mathbf{H}(t))\|_{\otimes} \leq 1/z$ , so  $\|\exp\{s\mathbf{H}_{z}(t)\}\|_{\otimes} = \|\exp\{-sz\}\exp\{sz^{2}\mathbf{R}(z,\mathbf{H}(t))\}\|_{\otimes} \leq 1$ . Now recall that  $\lim_{z\to\infty} \{z\mathbf{R}(z,\mathbf{H}(t))\Phi\} = \Phi$ ,  $\Phi \in \mathcal{FD}_{\otimes}$ , so that, for  $\Phi \in D$ , we have that  $\lim_{z\to\infty} \mathbf{H}_{z}(t)\Phi = \lim_{z\to\infty} \{z\mathbf{H}(t)\mathbf{R}(z,\mathbf{H}(t))\Phi\} = \lim_{z\to\infty} \{z\mathbf{R}(z,\mathbf{H}(t))\}\mathbf{H}(t)\Phi = \mathbf{H}(t)\Phi$ .

To prove 2., first recall (Gill<sup>73</sup>) that a tensor product norm,  $\|\cdot\|_{\otimes}$ , is uniform if, for  $\hat{\otimes}_{s\in E} T_s \in \mathbf{B}(\mathcal{H}_{\otimes})$ ,

$$\left\| \hat{\otimes}_{s \in \mathcal{E}} T_{s} \right\|_{\otimes} \le \prod_{s \in \mathcal{E}} \| T_{s} \|. \tag{3.2}$$

Using the uniform property of the (Hilbert space) tensor product norm, it is easy to see that  $\exp{\{\tau \mathbf{H}(t)\}}$  is a contraction semigroup.

To prove strong continuity, we need to identify a dense core for the family  $\{\mathbf{H}(t)|t\in E\}\subset \mathbf{C}(\mathcal{FD}_{\otimes})$ . Let  $D_1$  denote the ordered tensor product of the domains of the family  $\{H(t)|t\in E\}\subset \mathbf{C}(\mathcal{H})$ , (so that  $D_1\subset D$ )

$$D_{1} = \underset{s \in E}{\otimes} D(H(s)) = \left\{ \sum_{i=1}^{n} \underset{s}{\otimes} \varphi_{s}^{i} \middle| \varphi_{s}^{i} \in D(H(s)), s \in E \right\}.$$
 (3.3)

It is clear that  $D_1$  is a dense core in  $\mathcal{H}_{\otimes}$ , so  $D_0 = D_1 \cap \mathcal{FD}_{\otimes}$  is a dense core in  $\mathcal{FD}_{\otimes}$ . Using our standard basis, if  $\Phi, \Psi \in D_0$ ,  $\Phi = \sum_i \sum_{I(i)} a^i_{I(j)} E^i_{I(j)}$ ,  $\Psi = \sum_i \sum_{I(k)} b^i_{I(k)} E^i_{I(k)}$ ;

then, since  $\left(\exp\{\tau \mathbf{H}(t)\} - I_{\otimes}\right)$  is invariant on  $\mathcal{FD}^{i}$  and  $I_{\otimes}$  is the identify on  $\mathcal{FD}_{\otimes}$ , we have

$$\left\langle \left( \exp\{\tau \mathbf{H}(t)\} - I_{\otimes} \right) \Phi, \Psi \right\rangle = \sum_{i} \sum_{\mathbf{I}(i)} \sum_{\mathbf{I}(k)} a^{i}_{\mathbf{I}(j)} \overline{b}^{i}_{\mathbf{I}(k)} \left\langle \left( \exp\{\tau \mathbf{H}(t)\} - I_{\otimes} \right) E^{i}_{\mathbf{I}(j)}, E^{i}_{\mathbf{I}(k)} \right\rangle, \quad (3.4a)$$

and

$$\left\langle \left( \exp\{\tau \mathbf{H}(t)\} - \mathbf{I}_{\otimes} \right) E_{\mathbf{I}(j)}^{i}, E_{\mathbf{I}(k)}^{i} \right\rangle = \prod_{s \neq t} \left\langle e_{s,j(s)}^{i}, e_{s,k(s)}^{i} \right\rangle \left\langle \left( \exp\{\tau H(t)\} - \mathbf{I} \right) e_{t,j(t)}^{i}, e_{t,k(t)}^{i} \right\rangle (3.4b)$$

$$= \left\langle \left( \exp\{\tau H(t)\} - \mathbf{I} \right) e_{t,j(t)}^{i}, e_{t,j(t)}^{i} \right\rangle (e.o.v),$$

$$= \left\langle \left( \exp\{\tau H(t)\} - \mathbf{I} \right) e^{i}, e^{i} \right\rangle (e.f.n.v.),$$

$$\Rightarrow \left\langle \left( \exp\{\tau H(t)\} - \mathbf{I}_{\otimes} \right) \Phi, \Psi \right\rangle = \sum_{i} \sum_{\mathbf{I}(i)} a_{\mathbf{I}(i)}^{i} \overline{b}_{\mathbf{I}(j)}^{i} \left\langle \left( \exp\{\tau H(t)\} - \mathbf{I} \right) e^{i}, e^{i} \right\rangle (a.s.c). (3.4c)$$

Since all sums are finite, we have

$$\lim_{\tau \to 0} \left\langle \left( \exp\{\tau \mathbf{H}(t)\} - \mathbf{I}_{\otimes} \right) \Phi, \Psi \right\rangle = \sum_{i} \sum_{\mathbf{I}(j)} a^{i}_{\mathbf{I}(j)} \overline{b}^{i}_{\mathbf{I}(j)} \left\{ \lim_{\tau \to 0} \left\langle \left( \exp\{\tau \mathbf{H}(t)\} - \mathbf{I} \right) e^{i}, e^{i} \right\rangle \right\} = 0 \text{ (a.s.c)}. (3.4d)$$
The if and only if part is now clear. Since  $\exp\{\tau \mathbf{H}(t)\}$  is bounded on  $\mathcal{H}_{\otimes}$  and

the above limit exists on  $D_0$  (which is dense in  $\mathcal{F}\mathcal{D}_{\otimes}$ ), we see that  $\exp\{\tau \mathbf{H}(t)\}$  extends to a contraction semigroup on  $\mathcal{F}\mathcal{D}_{\otimes}$ . Now use the fact that, if a bounded semigroup converges weakly to the identity, it converges strongly (see Pazy<sup>71</sup>, pg. 44).

We now assume that the family  $\{H(t)| t \in E\} \subset C(\mathcal{H})$  has a weak Riemann integral  $Q = \int_a^b H(t)dt \in C(\mathcal{H})$ . It follows that the family  $\{H_z(t)| t \in E\} \subset B(\mathcal{H})$  also has a weak Riemann integral  $Q_z = \int_a^b H_z(t)dt \in B(\mathcal{H})$ . Let  $P_n$  be a sequence of

partitions (of E) so that the mesh  $\mu(P_n) \to 0$  as  $n \to \infty$ . Set  $Q_{z,n} = \sum_{l=1}^n H_z(\bar{t}_l) \Delta t_l, \ Q_{z,m} = \sum_{q=1}^m H_z(\bar{s}_q) \Delta s_q; \ \mathbf{Q}_{z,n} = \sum_{l=1}^n \mathbf{H}_z(\bar{t}_l) \Delta t_l, \ \mathbf{Q}_{z,m} = \sum_{q=1}^m \mathbf{H}_z(\bar{s}_q) \Delta s_q; \ \text{and}$   $\Delta Q_z = Q_{z,n} - Q_{z,m}, \ \Delta \mathbf{Q}_z = \mathbf{Q}_{z,n} - \mathbf{Q}_{z,m}. \ \text{Let } \Phi, \Psi \in \mathbf{D}_0; \ \Phi = \sum_{i}^J \Phi^i = \sum_{i}^J \sum_{\mathbf{I}(j)}^K a^i_{\mathbf{I}(j)} E^i_{\mathbf{I}(j)},$   $\Psi = \sum_{i}^L \Psi^i = \sum_{i}^L \sum_{\mathbf{I}(k)}^M b^i_{\mathbf{I}(k)} E^i_{\mathbf{I}(k)}, \ \text{and set } \phi = \sum_{i}^J \sum_{\mathbf{I}(j)}^K a^i_{\mathbf{I}(j)} e^i \ \text{and} \ \Psi = \sum_{i}^J \sum_{\mathbf{I}(j)}^K b^i_{\mathbf{I}(j)} e^i. \ \text{Then we}$ 

have:

**Theorem 3.2** (First Fundamental Theorem for Time-Ordered Integrals)

$$\langle \Delta \mathbf{Q}_{z} \Phi, \Psi \rangle = \sum_{i}^{J} \sum_{I(j)}^{K} a_{I(j)}^{i} \overline{b}_{I(j)}^{i} \langle \Delta Q_{z} e^{i}, e^{i} \rangle \text{ (a.s.c)}.$$
 (3.5)

**Note** The form of (3.5) is quite general since  $\Delta \mathbf{Q}_z$  can be replaced by other terms which also give a true relationship. For example, it is easy to show that the family  $\{\mathbf{H}_z(t)|t\in E\}$  is weakly measurable, weakly continuous, weakly differentiable, etc if and only if the same is true for the family  $\{H_z(t)|t\in E\}$ .

 $\begin{aligned} \textbf{Proof:} & \left\langle \Delta \mathbf{Q}_{z} \Phi, \Psi \right\rangle = \sum_{i} \sum_{\mathbf{I}(j)} \sum_{\mathbf{I}(k)} a_{\mathbf{I}(j)}^{\mathbf{i}} \bar{b}_{\mathbf{I}(k)}^{\mathbf{i}} \left\langle \Delta \mathbf{Q}_{z} E_{\mathbf{I}(j)}^{\mathbf{i}}, E_{\mathbf{I}(k)}^{\mathbf{i}} \right\rangle \text{ (we omit the upper limit). Now} \\ & \left\langle \Delta \mathbf{Q}_{z} E_{\mathbf{I}(j)}^{\mathbf{i}}, E_{\mathbf{I}(k)}^{\mathbf{i}} \right\rangle = \sum_{l=1}^{n} \Delta t_{l} \left\langle \mathbf{H}_{z}(\bar{t}_{l}) E_{\mathbf{I}(j)}^{\mathbf{i}}, E_{\mathbf{I}(k)}^{\mathbf{i}} \right\rangle - \sum_{q=1}^{m} \Delta s_{q} \left\langle \mathbf{H}_{z}(\bar{s}_{q}) E_{\mathbf{I}(j)}^{\mathbf{i}}, E_{\mathbf{I}(k)}^{\mathbf{i}} \right\rangle \\ & = \sum_{l=1}^{n} \Delta t_{l} \prod_{\mathbf{t} \neq \bar{t}_{l}} \left\langle e_{\mathbf{t}, \mathbf{j}(\mathbf{t})}^{\mathbf{i}}, e_{\mathbf{t}, \mathbf{k}(\mathbf{t})}^{\mathbf{i}} \right\rangle \left\langle H_{z}(\bar{t}_{l}) e_{\bar{t}_{l}, \mathbf{j}(\bar{t}_{l})}^{\mathbf{i}}, e_{\bar{t}_{l}, \mathbf{k}(\bar{t}_{l})}^{\mathbf{i}} \right\rangle - \sum_{q=1}^{m} \Delta s_{q} \prod_{\mathbf{t} \neq \bar{s}_{q}} \left\langle e_{\mathbf{t}, \mathbf{j}(\mathbf{t})}^{\mathbf{i}}, e_{\mathbf{t}, \mathbf{k}(\mathbf{t})}^{\mathbf{i}} \right\rangle \left\langle H_{z}(\bar{s}_{q}) e_{\bar{s}_{q}, \mathbf{j}(\bar{s}_{q})}^{\mathbf{i}}, e_{\bar{s}_{q}, \mathbf{j}(\bar{s}_{q})}^{\mathbf{i}} \right\rangle \\ & = \sum_{l=1}^{n} \Delta t_{l} \left\langle H_{z}(\bar{t}_{l}) e_{\bar{t}_{l}, \mathbf{j}(\bar{t}_{l})}^{\mathbf{i}}, e_{\bar{t}_{l}, \mathbf{j}(\bar{t}_{l})}^{\mathbf{i}} \right\rangle - \sum_{q=1}^{m} \Delta s_{q} \left\langle H_{z}(\bar{s}_{q}) e_{\bar{s}_{q}, \mathbf{j}(\bar{s}_{q})}^{\mathbf{i}}, e_{\bar{s}_{q}, \mathbf{j}(\bar{s}_{q})}^{\mathbf{i}} \right\rangle \\ & = \left\langle \Delta Q_{z} e^{\mathbf{i}}, e^{\mathbf{i}} \right\rangle \text{ (e.f.n.v). This result leads to (3.5). \end{aligned}$ 

**Theorem 3.3** (Second Fundamental Theorem for Time-Ordered Integrals) *If the family*  $\{H_z(t)|t\in E\}$  *has a weak Riemann (Riemann-Complete) integral, then* 

- 1. the family  $\{\mathbf{H}_z(t)|t\in E\}\subset \mathbf{B}^\#(\mathcal{FD}_\otimes)$  has a weak Riemann (Riemann-Complete) integral.
- 2. *If, in addition, we assume that for each*  $\Phi$  *with*  $\|\Phi\| = 1$ ,

$$\sup_{t \in \mathbb{E}} \left| \int_{a}^{t} (\left\| \mathbf{H}_{z}(s)\Phi \right\|^{2} - \left| \langle \mathbf{H}_{z}(s)\Phi, \Phi \rangle \right|^{2}) ds \right| < \infty, \tag{3.6}$$

then the family  $\{\mathbf{H}_z(t)|t\in E\}$  has a strong integral  $\mathbf{Q}_z[t,a] = \int_a^t \mathbf{H}_z(s)ds$  which generates a uniformly continuous contraction semigroup on  $\mathcal{FD}_{\otimes}$ .

#### **Notes:**

- **1.** It is sufficient that  $\sup_{t \in E} \left| \int_{a}^{t} (\|\mathbf{H}_{z}(s)E^{i}\|^{2} |\langle \mathbf{H}_{z}(s)E^{i}, E^{i} \rangle|^{2}) ds \right| < \infty$  for each i.
- **2.** Condition (3.6) is satisfied if  $\|\mathbf{H}_{z}(s)E^{i}\|^{2}$  is Lebesgue integrable for each i. In this case, we replace the Riemann integral by the Riemann-Complete integral.
- **3.** In general, the family  $\{\mathbf{H}_z(t)|t\in E\}$  need not be a Bochner or Pettis integral, as it is not required that  $\|\mathbf{H}_z(t)\Phi\|, \langle \mathbf{H}_z(t)\Phi, \Phi \rangle$  be (square) Lebesgue integrable. It is possible that  $\int_a^b \|\mathbf{H}_z(t)\Phi\|^2 dt = \infty$  &  $\int_a^b |\langle \mathbf{H}_z(t)\Phi, \Phi \rangle|^2 dt = \infty$ , while (3.6) is zero.

For example, let f(t) be any nonabsolutely (square) integrable function and set  $\mathbf{H}_z(t) = f(t)\mathbf{I}_{\otimes}$ . Then the above possibility holds while  $\int_a^t (\|\mathbf{H}_z(s)\Phi\|^2 - |\langle \mathbf{H}_z(s)\Phi, \Phi \rangle|^2) ds = 0 \text{ for all } t \text{ in E.}$ 

**Proof:** The proof of 1. is easy and follows from (3.5). To see that (3.6) makes  $Q_z$  a strong limit, let  $\Phi \in D_0$ . Then

$$\begin{split} \left\langle \mathbf{Q}_{z,n} \mathbf{\Phi}, \mathbf{Q}_{z,n} \mathbf{\Phi} \right\rangle &= \sum_{i=I(j),I(h)}^{J} \sum_{I(j)}^{K} a_{I(h)}^{i} \left( \sum_{k,m}^{n} \Delta t_{k} \Delta t_{m} \left\langle H_{z}(s_{k}) E_{I(j)}^{i}, H_{z}(s_{m}) E_{I(h)}^{i} \right\rangle \right) \\ &= \sum_{i}^{J} \sum_{I(j)}^{K} \left| a_{I(j)}^{i} \right|^{2} \left( \sum_{k \neq m}^{n} \Delta t_{k} \Delta t_{m} \left\langle H_{z}(s_{k}) e_{s_{k},j(s_{k})}^{i}, e_{s_{k},j(s_{k})}^{i} \right\rangle \left\langle e_{s_{m},j(s_{m})}^{i}, H_{z}(s_{m}) e_{s_{m},j(s_{k})}^{i} \right\rangle \right) \\ &+ \sum_{i}^{J} \sum_{I(j)}^{K} \left| a_{I(j)}^{i} \right|^{2} \left( \sum_{k}^{n} (\Delta t_{k})^{2} \left\langle H_{z}(s_{k}) e_{s_{k},j(s_{k})}^{i}, H_{z}(s_{k}) e_{s_{k},j(s_{k})}^{i} \right\rangle \right). \end{split} \tag{3.7}$$

This can be rewritten as

$$\begin{split} \left\| \mathbf{Q}_{z,n} \Phi \right\|_{\otimes}^{2} &= \sum_{i=1}^{J} \sum_{I(j)}^{K} \left| a_{I(j)}^{i} \right|^{2} \left\{ \left| \left\langle Q_{z,n} e^{i}, e^{i} \right\rangle \right|^{2} \\ &+ \sum_{k=1}^{L} (\Delta t_{k})^{2} \left( \left\| H_{z}(s_{k}) e^{i} \right\|^{2} - \left| \left\langle H_{z}(s_{k}) e^{i}, e^{i} \right\rangle \right|^{2} \right) \right\}, (a.s.c). \end{split}$$
(3.8)

The last term can be written as

$$\left|\sum_{k}^{n} (\Delta t_{k})^{2} \left( \left\| H_{z}(s_{k}) e^{i} \right\|^{2} - \left| \left\langle H_{z}(s_{k}) e^{i}, e^{i} \right\rangle \right|^{2} \right) \right| \leq \mu_{n} M \sup_{t \in \mathbb{E}} \left| \int_{a}^{t} \left( \left\| H_{z}(s) e^{i} \right\|^{2} - \left| \left\langle H_{z}(s) e^{i}, e^{i} \right\rangle \right|^{2} \right) ds \right|,$$

where M is a constant and  $\mu_n$  is the mesh of  $P_n$ , with  $\mu_n \to 0$  as  $n \to \infty$ . Now note that  $\|\mathbf{H}_z(t)E^i\|_{\otimes} = \|H_z(t)e^i\|$  and  $\langle \mathbf{H}_z(t)E^i, E^i\rangle = \langle H_z(t)e^i, e^i\rangle$  (e.o.v) so that

$$\sup_{t \in E} \left| \int_{a}^{t} \left( \left\| H_{z}(s) e^{i} \right\|^{2} - \left| \left\langle H_{z}(s) e^{i}, e^{i} \right\rangle \right|^{2} \right) ds \right| = \sup_{t \in E} \left| \int_{a}^{t} \left( \left\| \mathbf{H}_{z}(s) E^{i} \right\|^{2} - \left| \left\langle \mathbf{H}_{z}(s) E^{i}, E^{i} \right\rangle \right|^{2} \right) ds \right|$$
 (a.s.c).

We can now use (3.6) to get

$$\left\|\mathbf{Q}_{z,n}\mathbf{\Phi}\right\|_{\otimes}^{2} \leq \sum_{i=1}^{J} \sum_{I(j)}^{K} \left|a_{I(j)}^{i}\right|^{2} \left\{ \left|\left\langle Q_{z,n}e^{i}, e^{i}\right\rangle\right|^{2} + \mu_{n}M \sup_{t} \left|\int_{a}^{t} \left(\left\|\mathbf{H}_{z}(t)E^{i}\right\|^{2} - \left|\left\langle\mathbf{H}_{z}(t)E^{i}, E^{i}\right\rangle\right|^{2}\right) ds\right|_{s}\right\}, (a.s.c).$$

Thus,  $\mathbf{Q}_{\mathbf{z},\mathbf{n}}\Phi$  converges strongly to  $\mathbf{Q}_{\mathbf{z}}\Phi$  on  $D_{0}$  and hence has a strong limit on  $\mathcal{F}\mathcal{D}_{\otimes}$ . To show that  $Q_{\mathbf{z}}[t,a]$  generates a uniformly continuous contraction, it suffices to show that  $Q_{\mathbf{z}}[t,a]$  and  $Q_{\mathbf{z}}^{*}[t,a]$  are dissipative. Let  $\Phi$  be in  $D_{0}$ , then  $\langle \mathbf{Q}_{\mathbf{z}}[t,a]\Phi,\Phi\rangle = \sum_{i}^{J}\sum_{\mathbf{l}(j)}^{K}a_{\mathbf{l}(j)}^{i}\overline{b}_{\mathbf{l}(j)}^{i}\langle Q_{\mathbf{z}}e^{i},e^{i}\rangle$  (a.s.c) and, since  $Q_{\mathbf{z},\mathbf{n}}[t,a]$  is disspative for

each n. we have

$$\langle Q_z[t,a]e^i,e^i\rangle = \langle Q_{z,n}[t,a]e^i,e^i\rangle + \langle [Q_z[t,a]-Q_{z,n}[t,a]]e^i,e^i\rangle \leq \langle [Q_z[t,a]-Q_{z,n}[t,a]]e^i,e^i\rangle$$
  
Letting  $n \to \infty$ , we get  $\langle Q_z[t,a]e^i,e^i\rangle \leq 0$ , so that  $\langle \mathbf{Q}_z[t,a]\Phi,\Phi\rangle \leq 0$ . The same argument applies to  $\mathbf{Q}_z^*[t,a]$ . Since  $\mathbf{Q}_z[t,a]$  is dissipative and densely defined, it has a (bounded) dissipative closure on  $\mathcal{FD}_{\otimes}$ .

It should be noted that the theorem is still true if we allow the approximating sums for condition (3.6) to diverge but at an order less than  $\mu_n^{-1+\delta}, \ 0 < \delta < 1, \text{ that is, } \sup_t \left| \int_a^t \left( \left\| \mathbf{H}_z(t) E^i \right\|^2 - \left| \left\langle \mathbf{H}_z(t) E^i, E^i \right\rangle \right|^2 \right) ds \right| = \infty, \text{ with }$   $\left| \sum_{k}^n (\Delta t_k)^2 \left( \left\| H_z(s_k) e^i \right\|^2 - \left| \left\langle H_z(s_k) e^i, e^i \right\rangle \right|^2 \right) \right| \le M \mu_n^{\delta}.$ 

We also note that:

$$\|\mathbf{Q}_{z}[t,a]\Phi\|_{\otimes}^{2} = \sum_{i}^{J} \sum_{I(j)}^{K} |a_{I(j)}^{i}|^{2} |\langle Q_{z}e^{i}, e^{i} \rangle|^{2} \quad (a.s.c),$$
 (3.9)

in either of the above cases. This representation makes it easy to prove the next theorem.

#### Theorem 3.4

1. 
$$\mathbf{Q}_{z}[t,s] + \mathbf{Q}_{z}[s,a] = \mathbf{Q}_{z}[t,a]$$
 (a.s.c),

2. 
$$s - \lim_{h \to 0} \frac{\mathbf{Q}_z[t+h,a] - \mathbf{Q}_z[t,a]}{h} = s - \lim_{h \to 0} \frac{\mathbf{Q}_z[t+h,t]}{h} = \mathbf{H}_z(t)$$
 (a.s.c),

3. 
$$s - \lim_{h \to 0} \mathbf{Q}_z[t+h,t] = 0$$
 (a.s.c),

4. 
$$s - \lim_{h \to 0} \exp \{ \tau \mathbf{Q}_z[t+h,t] \} = I_{\otimes} \text{ (a.s.c)}, \tau \ge 0.$$

**Proof:** In each case, it suffices to prove the result for  $\Phi \in D_0$ . To prove 1., use

$$\begin{split} \left\| \left[ \mathbf{Q}_{z}[t,s] + \mathbf{Q}_{z}[s,a] \right] \Phi \right\|_{\otimes}^{2} &= \sum_{i}^{J} \sum_{I(j)}^{K} \left| a_{I(j)}^{i} \right|^{2} \left| \left[ Q_{z}[t,s] + Q_{z}[s,a] \right] e^{i}, e^{i} \right|^{2} (\text{a.s.c}) \\ &= \sum_{i}^{J} \sum_{I(j)}^{K} \left| a_{I(j)}^{i} \right|^{2} \left| \left\langle Q_{z}[t,a] e^{i}, e^{i} \right\rangle \right|^{2} = \left\| \mathbf{Q}_{z}[t,a] \Phi \right\|_{\otimes}^{2} (\text{a.s.c}). \end{split}$$

To prove 2., use 1 to get 
$$\mathbf{Q}_{z}[t+h,a] - \mathbf{Q}_{z}[t,a] = \mathbf{Q}_{z}[t+h,t]$$
 (a.s.), so that 
$$\lim_{h \to 0} \left\| \frac{\mathbf{Q}_{z}[t+h,t]}{h} \mathbf{\Phi} \right\|_{\otimes}^{2} = \sum_{i=1}^{J} \sum_{I(j)}^{K} \left| a_{I(j)}^{i} \right|^{2} \lim_{h \to 0} \left| \left\langle \frac{Q_{z}[t+h,t]}{h} e^{i}, e^{i} \right\rangle \right|^{2} = \left\| \mathbf{H}_{z}(t) \mathbf{\Phi} \right\|_{\otimes}^{2}$$
 (a.s.c.).

The proof of 3., follows from 2., and the proof of 4. follows from 3.

**Theorem 3.5** Suppose that  $\lim_{z\to\infty} \langle Q_z[t,a]\phi,\psi\rangle = \langle Q[t,a]\phi,\psi\rangle$  exists for  $\phi$  in a dense set  $\forall \psi \in \mathcal{H}$  (weak convergence). Then:

- 1. Q[t,a] generates a strongly continuous contraction semigroup on  $\mathcal{H}$ ,
- 2.  $\lim_{z \to \infty} \mathbf{Q}_z[t,a] \Phi = \mathbf{Q}[t,a] \Phi$  for  $\Phi \in D_0$  and  $\mathbf{Q}[t,a]$  is the generator of a strongly continuous contraction semigroup on  $\mathcal{FD}_{\otimes}$ ,
- 3.  $\mathbf{Q}[t,s] + \mathbf{Q}[s,a] = \mathbf{Q}[t,a]$  (a.s.c.),
- 4.  $\lim_{h\to 0} \frac{\mathbf{Q}[t+h,a] \mathbf{Q}[t,a]}{h} \Phi = \lim_{h\to 0} \frac{\mathbf{Q}[t+h,t]}{h} \Phi = \mathbf{H}(t) \Phi \text{ (a.s.c.)},$
- 5.  $\lim_{h\to 0} \mathbf{Q}[t+h,t]\Phi = 0$  (a.s.c.), and
- 6.  $\lim_{h\to 0} \exp\{\tau \mathbf{Q}[t+h,t]\}\Phi = \Phi \text{ (a.s.c.)}, \tau \ge 0.$

**Proof:** The proofs are easy. For 1., first note that Q[t,a] is closable and use  $\langle Q[t,a]\phi,\phi\rangle = \langle Q_z[t,a]\phi,\phi\rangle + \langle [Q[t,a]-Q_z[t,a]]\phi,\phi\rangle \leq \langle [Q[t,a]-Q_z[t,a]]\phi,\phi\rangle$  and let

 $z \to \infty$ . Then do likewise for  $\langle \phi, Q^*[t, a] \phi \rangle$  to get that Q[t, a] is maximal dissipative. To prove 2., use (3.9) in the form

$$\left\| \left[ \mathbf{Q}_{z}[t,a] - \mathbf{Q}_{z'}[t,a] \right] \Phi \right\|_{\otimes}^{2} = \sum_{i}^{J} \sum_{I(i)}^{K} \left| a_{I(i)}^{i} \right|^{2} \left| \left[ Q_{z}[t,a] - Q_{z'}[t,a] \right] e^{i}, e^{i} \right|^{2}, \text{ (a.s.c.)}.$$

This proves that  $\mathbf{Q}_{z}[t,a] \xrightarrow{s} \mathbf{Q}[t,a]$ . Since  $\mathbf{Q}[t,a]$  is densely defined, it is closable. The same method as above shows that it is maximal dissipative. Proofs of the other results follow the methods of the previous theorem.

Since  $\mathbf{Q}[t,a]$  and  $\mathbf{Q}_z[t,a]$  generate contraction semigroups, set  $\mathbf{U}[t,a] = \exp{\{\mathbf{Q}[t,a]\}}$ ,  $\mathbf{U}_z[t,a] = \exp{\{\mathbf{Q}[t,a]\}}$ , for  $t \in E$ . They are evolution operators and the following theorem is a slight modification of a result due to Hille and Phillips<sup>74</sup>, known as the second exponential formula.

**Theorem 3.6** If  $\mathbf{Q}^{\iota}[t,a] = w\mathbf{Q}[t,a]$  is the generator of a strongly continuous contraction semigroup, and  $\mathbf{U}^{w}[t,a] = \exp\{w\mathbf{Q}[t,a]\}$ , then, for each n and  $\Phi \in D[(\mathbf{Q}[t,a])^{n+1}]$ , we have (where w is a parameter)

$$\mathbf{U}^{w}[t,a]\Phi = \left\{ I_{\otimes} + \sum_{k=1}^{n} \frac{\left(w\mathbf{Q}[t,a]\right)^{n}}{n!} + \frac{1}{n!} \int_{0}^{w} (w-\xi)^{n} \mathbf{Q}[t,a]^{n+1} \mathbf{U}^{\xi}[t,a] d\xi \right\} \Phi.$$
 (3.10)

**Proof:** The proof is easy. Start with  $\left[\mathbf{U}_{z}^{w}[t,a]\Phi - I_{\otimes}\right]\Phi = \int_{0}^{w} \mathbf{Q}_{z}[t,a]\mathbf{U}_{z}^{\xi}[t,a]d\xi\Phi$  and use integration by parts to get that

$$\left[\mathbf{U}_{z}^{w}[t,a]\Phi - I_{\otimes}\right]\Phi = w\mathbf{Q}_{z}[t,a]\Phi + \int_{0}^{w} (w-\xi)\left[\mathbf{Q}_{z}[t,a]\right]^{2}\mathbf{U}_{z}^{\xi}[t,a]d\xi\Phi.$$

It is clear how to get the n-th term. Finally, let  $z \to \infty$  to get (3.10).

**Theorem 3.7** *If* a < t < b,

1. 
$$\lim_{z\to\infty} \mathbf{U}_z[t,a]\Phi = \mathbf{U}[t,a]\Phi, \ \Phi \in \mathcal{FD}_{\otimes},$$

2. 
$$\frac{\partial}{\partial t} \mathbf{U}_{z}[t,a]\Phi = \mathbf{H}_{z}(t)\mathbf{U}_{z}[t,a]\Phi = \mathbf{U}_{z}[t,a]\mathbf{H}_{z}(t)\Phi, \ \Phi \in \mathcal{FD}_{\otimes},$$
 and

3. 
$$\frac{\partial}{\partial t}\mathbf{U}[t,a]\Phi = \mathbf{H}(t)\mathbf{U}[t,a]\Phi = \mathbf{U}[t,a]\mathbf{H}(t)\Phi, \ \Phi \in \mathbf{D}(\mathbf{Q}[b,a]) \supset \mathbf{D}_0.$$

**Proof:** To prove 1., use the fact that  $\mathbf{H}_{x}(t)$  and  $\mathbf{H}(t)$  commute along with

$$\mathbf{U}[t,a]\Phi - \mathbf{U}_{z}[t,a]\Phi = \int_{0}^{1} (d/ds) \left(e^{s\mathbf{Q}[t,a]}e^{(1-s)\mathbf{Q}_{z}[t,a]}\right) \Phi ds$$
$$= \int_{0}^{1} s \left(e^{s\mathbf{Q}[t,a]}e^{(1-s)\mathbf{Q}_{z}[t,a]}\right) \left(\mathbf{Q}[t,a] - \mathbf{Q}_{z}[t,a]\right) \Phi ds, \text{ so that}$$

$$\|\mathbf{U}[t,a]\Phi - \mathbf{U}_{z}[t,a]\Phi\| \le \|\mathbf{Q}[t,a]\Phi - \mathbf{Q}_{z}[t,a]\Phi\|$$

To prove 2., use

$$\mathbf{U}_{z}[t+h,a] - \mathbf{U}_{z}[t,a] = \mathbf{U}_{z}[t,a] \left(\mathbf{U}_{z}[t+h,t] - \mathbf{I}\right) = \left(\mathbf{U}_{z}[t+h,t] - \mathbf{I}\right)\mathbf{U}_{z}[t,a], \text{ so that,}$$

$$\frac{\left(\mathbf{U}_{z}[t+h,a] - \mathbf{U}_{z}[t,a]\right)}{h} = \mathbf{U}_{z}[t,a] \frac{\left(\mathbf{U}_{z}[t+h,t] - \mathbf{I}\right)}{h}.$$

Now set  $\Phi_z^t = \mathbf{U}_z[t,a]\Phi$  and use (3.10) with n = 1 and w = 1 to get:

$$\mathbf{U}_{z}[t+h,t]\boldsymbol{\Phi}_{z}^{t} = \left\{ I_{\otimes} + \mathbf{Q}_{z}[t+h,t] + \int_{0}^{1} (1-\xi) \mathbf{U}_{z}^{\xi}[t+h,t]\mathbf{Q}_{z}[t+h,t]^{2} d\xi \right\} \boldsymbol{\Phi}_{z}^{t}, \text{ so that}$$

$$\frac{\left(\mathbf{U}_{z}[t+h,t]-\mathbf{I}\right)}{h} \boldsymbol{\Phi}_{z}^{t} - \mathbf{H}_{z}(t)\boldsymbol{\Phi}_{z}^{t} = \frac{\mathbf{Q}_{z}[t+h,t]}{h} \boldsymbol{\Phi}_{z}^{t} - \mathbf{H}_{z}(t)\boldsymbol{\Phi}_{z}^{t}$$

$$+ \int_{0}^{1} (1-\xi) \mathbf{U}_{z}^{\xi}[t+h,t] \frac{\mathbf{Q}_{z}[t+h,t]^{2}}{h} \boldsymbol{\Phi}_{z}^{t} d\xi$$

It follows that

$$\left\| \frac{\left(\mathbf{U}_{z}[t+h,t]-\mathbf{I}\right)}{h} \boldsymbol{\Phi}_{z}^{t} - \mathbf{H}_{z}(t) \boldsymbol{\Phi}_{z}^{t} \right\|_{\otimes} \leq \left\| \frac{\mathbf{Q}_{z}[t+h,t]}{h} \boldsymbol{\Phi}_{z}^{t} - \mathbf{H}_{z}(t) \boldsymbol{\Phi}_{z}^{t} \right\|_{\otimes} + \frac{1}{2} \left\| \frac{\mathbf{Q}_{z}[t+h,t]^{2}}{h} \boldsymbol{\Phi}_{z}^{t} \right\|_{\otimes}.$$

The result now follows from Theorem (3.4)-2 and 3.

To prove 3., note that  $\mathbf{H}_z(t)\Phi = \mathbf{H}(t)\{z\mathbf{R}(z,\mathbf{H}(t))\}\Phi = \{z\mathbf{R}(z,\mathbf{H}(t))\}\mathbf{H}(t)\Phi$ , so that  $\{z\mathbf{R}(z,\mathbf{H}(t))\}$  commutes with  $\mathbf{U}[t,a]$  and  $\mathbf{H}(t)$ . Now show that

$$\begin{split} \left\| \mathbf{H}_{\mathbf{z}}(t)\mathbf{U}_{\mathbf{z}}[t,a]\Phi - \mathbf{H}_{\mathbf{z}'}(t)\mathbf{U}_{\mathbf{z}'}[t,a]\Phi \right\| &\leq \left\| \left[ \mathbf{U}_{\mathbf{z}}[t,a]\Phi - \mathbf{U}_{\mathbf{z}'}[t,a] \right] \mathbf{H}(t)\Phi \right\| \\ &+ \left\| \left[ \mathbf{z}\mathbf{R}(\mathbf{z},\mathbf{H}(t))\Phi - \mathbf{z}' \, \mathbf{R}(\mathbf{z}',\mathbf{H}(t)) \right] \mathbf{H}(t)\Phi \right\| \to 0, \ \mathbf{z},\mathbf{z}' \to \infty, \end{split}$$

so that, for 
$$\Phi \in D(\mathbf{Q}[b,a])$$
,  $\mathbf{H}_z(t)\mathbf{U}_z[t,a]\Phi \to \mathbf{H}(t)\mathbf{U}[t,a]\Phi = \frac{\partial}{\partial t}\mathbf{U}[t,a]\Phi$ .

The previous theorems form the core of our approach to the Feynman operator calculus. Our theory applies to both hyperbolic and parabolic equations. In the conventional approach, these two cases require different methods (see Pazy<sup>71</sup>). It is not hard to show that the requirements imposed in these cases are stronger than (our condition of) weak integral. This will be discussed in a later paper devoted to the general problem on Banach spaces.

## 4. Perturbation Theory

**Definition 4.1** The evolution operator  $\mathbf{U}^{w}[t,a] = \exp\{w\mathbf{Q}[t,a]\}$ , is said to be asymptotic in the sense of Poincaré, if for each n and each  $\Phi_a \in D[(\mathbf{Q}[t,a])^{n+1}]$ , we have

$$\lim_{w \to 0} w^{-(n+1)} \left\{ \mathbf{U}^{w}[t,a] - \sum_{k=1}^{n} \frac{\left(w\mathbf{Q}[t,a]\right)^{k}}{k!} \right\} \Phi_{a} = \frac{\mathbf{Q}[t,a]^{n+1}}{(n+1)!} \Phi_{a}. \tag{4.1}$$

This is the operator version of an asymptotic expansion in the classical sense, but here  $\mathbf{Q}[t,a]$  is an unbounded operator.

As noted earlier, Dyson16 analyzed the (renormalized) perturbation expansion for quantum electrodynamics and suggested that it actually diverges. He concluded that we could, at best, hope that the series is asymptotic. His arguments were based on (not completely convincing) physical considerations, but no precise formulation of the problem was possible at that time. However, the calculations of Hurst75, Thirring76, Peterman77, and Jaffe78 for specific models all support Dyson's contention that the renormalized perturbation series diverges. In his recent book91 (pg. 13-16), Dyson's views on the perturbation series and renormalization are reiterated: " … in spite of all the successes of the new physics, the two questions that defeated me in 1951 remain unsolved." Here, he is referring to the question of mathematical consistency for the whole renormalization program, and our ability to (reliably) calculate nuclear processes in quantum chromodynamics. (For other details and references to additional works, see Schweber6,80, Wightman84 and Zinn-Justin79.)

The general construction of a physically simple and mathematically satisfactory formulation of quantum electrodynamics is still an open problem. The next theorem establishes Dyson's (second) conjecture under conditions that would apply to any (future) theory that does not require a radical departure from the present foundations of quantum theory (unitary solution operators). It also applies to the renormalized expansions in some areas of condensed matter physics where the solution operators are contraction semigroups.

**Theorem 4.2** *Suppose the conditions for Theorem 3.5 are satisfied. Then:*

1. **U Q** *<sup>w</sup>*[ ]= *ta w ta* , [, ] exp{ } *is asymptotic in the sense of Poincaré.*

2. For each n and each  $\Phi_a \in D[(\mathbf{Q}[t,a])^{n+1}]$ , we have

$$\Phi(t) = \Phi_{a} + \sum_{k=1}^{n} w^{k} \int_{a}^{t} ds_{1} \int_{a}^{s_{1}} ds_{2} \cdots \int_{a}^{s_{k-1}} ds_{k} \mathbf{H}(s_{1}) \mathbf{H}(s_{2}) \cdots \mathbf{H}(s_{k}) \Phi_{a} 
+ \int_{0}^{w} (w - \xi)^{n} d\xi \int_{a}^{t} ds_{1} \int_{a}^{s_{1}} ds_{2} \cdots \int_{a}^{s_{n}} ds_{n+1} \mathbf{H}(s_{1}) \mathbf{H}(s_{2}) \cdots \mathbf{H}(s_{n+1}) \mathbf{U}^{\xi}[s_{n+1}, a] \Phi_{a},$$
(4.2)

where  $\Phi(t) = \mathbf{U}^w[t,a]\Phi_a$ .

**Proof:** From (3.10), we have

$$\mathbf{U}^{w}[t,a]\Phi = \left\{ \sum_{k=0}^{n} \frac{\left(w\mathbf{Q}[t,a]\right)^{n}}{n!} + \frac{1}{n!} \int_{0}^{w} (w-\xi)^{n} \mathbf{Q}[t,a]^{n+1} \mathbf{U}^{\xi}[t,a] d\xi \right\} \Phi,$$

so that

$$w^{-(n+1)} \left\{ \mathbf{U}^{w}[t,a] \Phi_{a} - \sum_{k=0}^{n} \frac{\left(w \mathbf{Q}[t,a]\right)^{k}}{k!} \Phi_{a} \right\} = + \frac{(n+1)}{(n+1)!} w^{-(n+1)} \int_{0}^{w} (w - \xi)^{n} d\xi \mathbf{U}^{\xi}[t,a] \mathbf{Q}[t,a]^{n+1} \Phi_{a}.$$

Replace the right hand side by

$$I = \frac{(n+1)}{(n+1)!} w^{-(n+1)} \int_{0}^{w} (w-\xi)^{n} d\xi \Big\{ \mathbf{U}_{z}^{\xi}[t,a] + \Big[ \mathbf{U}^{\xi}[t,a] - \mathbf{U}_{z}^{\xi}[t,a] \Big] \Big\} \mathbf{Q}[t,a]^{n+1} \Phi_{a}.$$

Now, expand the term  $\mathbf{U}_z^{\xi}[t,a]$  in a two-term Taylor series about zero to get

$$\mathbf{U}_{z}^{\xi}[t,a] = I_{\otimes} + \xi \mathbf{Q}_{z}[t,a] + R_{z}^{\xi}.$$

Put the above in I, compute the elementary integrals showing that only the  $I_{\otimes}$  term gives a nonzero value (of 1/(n+1)) when  $w \to 0$ . Then let  $z \to \infty$  to get

$$\lim_{w\to 0} (n+1)w^{-(n+1)} \int_{0}^{w} d\xi (w-\xi)^{n} \mathbf{U}^{\xi}[t,a] \mathbf{Q}[t,a]^{n+1} \Phi_{a} = \mathbf{Q}[t,a]^{n+1} \Phi_{a}.$$

This proves that  $\mathbf{U}[t,a] = \exp{\{\mathbf{Q}[t,a]\}}$  is asymptotic in the sense of Poincaré. To prove (4.2), let  $\Phi_a \in D[(\mathbf{Q}[t,a])^{n+1}]$  for each  $k \le n+1$ , and use the fact that (Dollard and Friedman<sup>81</sup>)

$$\left(\mathbf{Q}_{z}[t,a]\right)^{k}\boldsymbol{\Phi}_{a} = \left(\int_{a}^{t}\mathbf{H}_{z}(s)ds\right)^{k}\boldsymbol{\Phi}_{a} = (k!)\int_{a}^{t}ds_{1}\int_{a}^{s_{1}}ds_{2}\cdots\int_{a}^{s_{k-1}}ds_{n}\mathbf{H}_{z}(s_{1})\mathbf{H}_{z}(s_{2})\cdots\mathbf{H}_{z}(s_{k})\boldsymbol{\Phi}_{a}.$$
(4.3)

Letting  $z \rightarrow \infty$  gives the result.

Our conditions are very weak. For example, the recent work of Tang and Li<sup>82</sup> required that ||H(t)|| be Lebesgue integrable.

There are well known special cases in which the perturbation series may actually converge to the solution. This can happen, for example, if the generator is bounded or if it is analytic in some sector. More generally, when the generator is of the form  $\mathbf{H}(t) = \mathbf{H}_0(t) + \mathbf{H}_i(t)$ , where  $\mathbf{H}_0(t)$  is analytic and  $\mathbf{H}_i(t)$  is some reasonable perturbation, which need not be bounded, there are conditions that allow the interaction representation to have a convergent Dyson expansion. These results can be formulated and proven in our formalism. However, the proofs are essentially the same as in the standard case so we will present them in a later paper devoted to the operator calculus on Banach spaces. The recent book by Engel and Nagel<sup>83</sup> provides some new results in this general area.

There are also cases where the (renormalized) series may diverge, but still respond to some summability method. This phenomenon is well known in classical analysis. In field theory, things can be much more complicated. A good discussion, with references, can be found in the review by Wightman84 and the book by Glimm and Jaffe34.

#### **5. Sum Over Paths**

In this section we first review and make a distinction between what is actually known and what we think we know about the foundations for our physical view of the micro-world. The objective is to provide the background for a number of physically motivated postulates that will be used to develop a theory of measurement for the micro-world (sufficient for our purposes). This will allow us to relate the theory of Sections 3 and 4 to Feynman's sum over paths approach and prove Dyson's second conjecture. This section differs from the previous ones in that we shift the orientation and perspective from that of mathematical physics to that of theoretical physics.

 In spite of the enormous successes of the physical sciences in the past century, our information and understanding about the micro-world is still rather meager. In the macro-world we are quite comfortable with the view that physical systems evolve continuously in time and our results justify this view. Indeed, the success of continuum physics is the basis for a large part of our technical advances in the twentieth century. On the other hand, the same view is also held at the micro-level and, in this case, our position is not very secure. The ability to measure physical events continuously in time at the micro-level must be considered a belief which, although convenient, has no place in science as an a priori constraint.

In order to establish perspective, let us consider this belief within the context of a satisfactory, and well-justified theory, Brownian motion. This theory lies at the interface between the macro- and the micro-worlds. Some presentations of this theory (the careful ones) make a distinction between the mathematical and the physical foundations of Brownian motion and that distinction is important for our discussion.

When Einstein85 began his investigation of the physical issues associated with this phenomenon, he was forced to assume that physical information about the state of a Brownian particle (position, velocity, etc) can only be known in time intervals that are large compared with the mean time between molecular collisions. (It is known that, under normal physical conditions, a Brownian particle receives about 1021 collisions per second.) Wiener took the mathematical step and assumed that this mean time (between collisions) could be made zero, thus providing a mathematical Brownian particle. This corresponds physically to the assumption that the ratio of the mass of the particle to the friction of the fluid is zero in the limit (see Wiener et al86).

From the physical point of view, use of Wiener's idealization of the Einstein model was not satisfactory since it led to problems of unbounded path length and nondifferentiability at all points. The first problem is physically impossible while the second is physically unreasonable. Of course, the idealization has turned out to be quite satisfactory in areas where the information required need not be detailed, such as large parts of electrical engineering, chemistry, and the biological sciences. Ornstein and Uhlenbeck87 later constructed a model that, gives the Einstein view asymptotically, but in small-time regions, is equivalent to the assumption that the particle travels a linear path between collisions. This model

provides finite path length and differentiability. (The theory was later idealized by Doob88.) What we do know is that the very nature of the liquid state implies collective behavior among the molecules. *This means that we do not know what path the particle travels in between collisions*. However, since the tools and methods of analysis require some form of continuity, some such (in between observation) assumptions must be made.It is clear that the need for these assumptions is imposed by the available mathematical structures within which we must represent physical reality as a model.

Theoretical science concerns itself with the construction of mathematical representations of certain restricted portions of physical reality. Various trends and philosophies that are prevalent at the time temper these constructs. A consistent theme has been the quest for simplicity. This requirement is born out of the natural need to restrict models to the minimum number of variables, relationships, constraints, etc, which give a satisfactory account of known experimental results and possibly allow the prediction of heretofore unknown consequences. One important outcome of this approach has been to implicitly eliminate all reference to the background within which physical systems evolve. In the micro-world, such an action cannot be justified without prior investigation. We propose to replace the use of mathematical coordinate systems by "physical coordinate systems" in order to (partially) remedy this problem.

We denote a physical coordinate system at time *t* by **R**<sup>p</sup> 3 ( )*t* . This coordinate system is attached to an observer (including measuring devices) and is envisioned as **R**<sup>3</sup> plus any background effects, either local or distant, which affect the observer's ability to obtain precise (ideal) experimental information about physical reality. This in turn affects our observer's ability to construct precise (ideal) representations and make precise predictions about physical reality (in the micro-world).

More specifically, consider the evolution of some micro-system on the interval E = [a,b]. Physically this evolution manifests itself as a curve on  $\mathbf{X}$ , where

$$\prod_{t\in E} \mathbf{R}_{p}^{3}(t) = \mathbf{X}.$$

Thus, true physical events occur on X where actual experimental information is modified by fluctuations in  $\mathbb{R}^3_p(t)$ , and by the interaction of the micro-system with the measuring equipment. Based on the success of our models, we know that such small changes are in the noise region, and they have no effect on our predictions for macro-systems. However, there is no a priori reason to believe that the effects will be small on micro-systems.

In terms of our theoretical representations, we are forced to model the evolution of physical systems in terms of wave functions, amplitudes, and/or operator-valued distributions, etc. There are thus two spaces, the physical space of evolution for the micro-system and the observer's space of obtainable information concerning this evolution. The lack of distinction between these two spaces seems to be the cause for some of the confusion and lack of physical clarity. For example, it may be perfectly correct to assume that a particle travels a continuous path on **X**. However, the assumption that the observer's space of obtainable information includes infinitesimal spacetime knowledge of this path is completely unfounded. This leads to our first postulate:

**Postulate 1.** Physical reality is a continuous process in time.

We thus take this view, fully recognizing that experiment does not provide continuous information about physical reality, and that there is no reason to believe that our mathematical representations contain precise information about the continuous spacetime behavior of physical processes at this level.

Since the advent of the special theory of relativity, there is much discussion about events, which generally means a point in **R**<sup>4</sup> with the Minkowski metric. In terms of real physics, this is a fiction which is frequently useful for reasons of presentation but so widely used that, to avoid confusion, it is appropriate to define what we mean by a *physical event.*

**Definition 5.1** *A physical event is a set of physical changes in a given system that can be verified directly by experiment or indirectly via subsequent changes, where conclusions are based on an a priori agreed-upon model of the physical process.*

This definition corresponds more closely to what is meant by physical events. It explicitly recognizes the evolution of scientific inference and the need for general agreement about what is being observed (based on specific models).

Before continuing, it will be helpful to have a particular physical picture in mind that makes the above discussion explicit. For this purpose, we take this picture to be a photograph showing the track left by a p-meson in a bubble chamber (and take seriously the amount of information available). In particular, we assume that the following reaction occurs:

$$\pi^+ \to \mu^+ + \nu$$
.

 We further assume that the orientation of our photograph is such that the pmeson enters on the left at time *t=0* and the tracks left by the µ-meson disappear on the right at time *t=T*, where *T* is of the order of 10-<sup>3</sup> sec, the time exposure for photographic film. Although the neutrino does not appear in the photograph, we also include a track for it. In Figure I we present a simplified picture of this photograph.

![](_page_36_Picture_2.jpeg)

**Figure I**

We have drawn the photograph as if we continuously see the particles in the picture. However, experiment only provides us individual bubbles, which do not necessarily overlap, from which we must extract physical information. A more accurate (though still not realistic) depiction is given in Figure II.

![](_page_37_Picture_0.jpeg)

**Figure II**

Let us assume that we have magnified a portion of our photograph to the extent that we may distinguish the individual bubbles created by the pmeson as it passes through the chamber. In Figure III, we present a simplified model of adjacent bubbles.

![](_page_37_Picture_3.jpeg)

**Figure III**

**Postulate 2.** *We assume that the center of each bubble represents the average knowable effect of the particle in a symmetric time interval about the center.*

By average knowable effect, we mean the average of the physical observables. In Fig III, we consider the existence of a bubble at time <sup>t</sup> *<sup>j</sup>* to be caused by the average of the physical observables over the time interval *t t* [ ] *j j* -1, , where *tj jj* - - 1 1 = + ( / )[ 1 2 <sup>t</sup> <sup>t</sup> ] and *tj jj* = + <sup>+</sup> ( / )[ 1 2 ] <sup>t</sup> <sup>t</sup> <sup>1</sup> . This postulate requires some justification. In general, the resolution of film and the relaxation time for distinct bubbles in the chamber vapor are limited. This means that if the p-meson creates two bubbles that are closely spaced in time, the bubbles may coalesce and appear as one. If this does not occur, it is still possible that the film will record the event as one bubble because of its inability to resolve events is such small time intervals.

Let us now recognize that we are dealing with one photograph so that, in order to obtain all available information, we must analyze a large number of photographs of the same reaction obtained under similar conditions (preprepared states). It is clear that the number of bubbles and the time placement of the bubbles will vary (independently of each other) from photograph to photograph. Let <sup>l</sup> -1 denote the average time for the appearance of a bubble in the film.

**Postulate 3.** *We assume that the number of bubbles in any film is a random variable.*

**Postulate 4.** *We assume that, given that n bubbles have appeared on a film, the time positions of the centers of the bubbles are uniformly distributed.*

**Postulate 5.** *We assume that N(t), the number of bubbles up to time t in a given film, is a Poisson-distributed random variable with parameter* <sup>l</sup> *.*

To motivate Postulate 5, recall that  $\tau_j$  is the time center of the j-th bubble and  $\lambda^{-1}$  is the average (experimentally determined) time between bubbles. The following results can be found in Ross<sup>89</sup>.

**Theorem 5.1** The random variables  $\Delta \tau_j = \tau_j - \tau_{j-1}$  ( $\tau_0 = 0$ ) are independent identically distributed random variables of exponential type with mean  $\lambda^{-1}$ , for  $1 \le j \le n$ .

The arrival times  $\tau_1, \tau_2, \dots, \tau_n$  are not independent, but their density function can be computed from

$$\operatorname{Prob}[\tau_{1}, \tau_{2}, \dots, \tau_{n}] = \operatorname{Prob}[\tau_{1}] \operatorname{Prob}[\tau_{2} \mid \tau_{1}] \dots \operatorname{Prob}[\tau_{n} \mid \tau_{1}, \tau_{2}, \dots, \tau_{n-1}]. \tag{5.1a}$$

We now use Theorem 5.1 to conclude that, for  $k \ge 1$ ,

$$\operatorname{Prob}[\tau_{k} \mid \tau_{1}, \tau_{2}, \dots, \tau_{k-1}] = \operatorname{Prob}[\tau_{k} \mid \tau_{k-1}]. \tag{5.1b}$$

We don't know this conditional probability however, the natural assumption is that given n bubbles appear, they are equally (uniformly) distributed on the interval. We can now construct what we call the experimental evolution operator. Assume that the conditions for Theorem 3.5 are satisfied and that the family  $\{\tau_1, \tau_2, \dots, \tau_n\}$  represents the time positions of the centers of n bubbles in our film of Fig III. Set a = 0 and define  $\mathbf{Q}_E[\tau_1, \tau_2, \dots, \tau_n]$  by

$$\mathbf{Q}_{E}[\tau_{1}, \tau_{2}, \dots, \tau_{n}] = \sum_{i=1}^{n} \int_{t_{j-1}}^{t_{j}} E[\tau_{j}, s] \mathbf{H}(s) ds.$$
 (5.2a)

Here,  $t_0 = \tau_0 = 0$ ,  $t_j = (1/2)[\tau_j + \tau_{j+1}]$  (for  $1 \le j \le n$ ), and  $E[\tau_j, s]$  is the exchange operator defined in Section 2. The effect of our exchange operator

 $E[\tau_j, s]$  is to concentrate all information contained in  $[t_{j-1}, t_j]$  at  $\tau_j$ . This is how we implement our postulate that the known physical event of the bubble at time  $\tau_j$  is due to an average of physical effects over  $[t_{j-1}, t_j]$  with information concentrated at  $\tau_j$ . We can rewrite  $\mathbf{Q}_E[\tau_1, \tau_2, \dots, \tau_n]$  as

$$\mathbf{Q}_{E}[\boldsymbol{\tau}_{1}, \boldsymbol{\tau}_{2}, \cdots, \boldsymbol{\tau}_{n}] = \sum_{j=1}^{n} \Delta t_{j} \left[ \frac{1}{\Delta t_{j}} \int_{t_{j-1}}^{t_{j}} E[\boldsymbol{\tau}_{j}, s] \mathbf{H}(s) ds \right].$$
 (5.2b)

Thus, we indeed have an average as required by Postulate 2. The evolution operator is given by

$$U[\tau_1, \tau_2, \dots, \tau_n] = \exp\left\{\sum_{j=1}^n \Delta t_j \left[\frac{1}{\Delta t_j} \int_{t_{j-1}}^{t_j} E[\tau_j, s] \mathbf{H}(s) ds\right]\right\}.$$
 (5.3a)

For  $\Phi \in \mathcal{FD}_{\otimes}$ , we define the function  $U[N(t), 0]\Phi$  by:

$$\mathbf{U}[N(t),0]\Phi = U\left[\tau_1, \tau_2, \dots, \tau_{N(t)}\right]\Phi. \tag{5.3b}$$

The function  $U[N(t),0]\Phi$  is a  $\mathcal{FD}_{\otimes}$ -valued random variable, which represents the distribution of the number of bubbles that may appear on our film up to time t. In order to relate  $U[N(t),0]\Phi$  to actual experimental results, we must compute its expected value. Using Postulates 3, 4, and 5, we have

$$\overline{\mathbf{U}}_{\lambda}[t,0]\Phi = \mathcal{E}[\mathbf{U}[N(t),0]\Phi] = \sum_{n=0}^{\infty} \mathcal{E}\{\mathbf{U}[N(t),0]\Phi | N(t) = n\} \operatorname{Prob}[N(t) = n], \tag{5.4a}$$

$$\mathcal{E}\left\{\mathbf{U}[N(t),0]\Phi | N(t) = n\right\} = \int_{0}^{t} \frac{d\tau_{1}}{t} \int_{\tau_{1}}^{t} \frac{d\tau_{2}}{t-\tau_{1}} \cdots \int_{\tau_{n-1}}^{t} \frac{d\tau_{n}}{t-\tau_{n-1}} \mathbf{U}[\tau_{n},\cdots,\tau_{1}]\Phi = \overline{\mathbf{U}}_{n}[t,0]\Phi, \quad (5.5a)$$

and

$$\operatorname{Prob}[N(t) = n] = \frac{\left(\lambda t\right)^n}{n!} \exp\{-\lambda t\}. \tag{5.6}$$

The integral in (5.4a) acts to distribute uniformly the time positions  $\tau_j$  over the successive intervals  $[t, \tau_{j-1}]$ ,  $1 \le j \le n$ , given that  $\tau_{j-1}$  has been determined. This is a natural result given our lack of knowledge.

The integral (5.4a) is of theoretical value but is not easy to compute. Since we are only interested in what happens when  $\lambda \to \infty$ , and as the mean number of bubbles in the film at time t is  $\lambda t$ , we can take  $\tau_j = (jt/n)$ ,  $1 \le j \le n$ ,  $(\Delta t_j = t/n \text{ for each n})$ . We can now replace  $\overline{\mathbf{U}}_n[t,0]\Phi$  by  $\mathbf{U}_n[t,0]\Phi$ , and with this understanding, we continue to use  $\tau_j$ , so that

$$\mathbf{U}_{n}[t,0]\Phi = \exp\left\{\sum_{j=1}^{n} \int_{t_{j-1}}^{t_{j}} E[\tau_{j},s]\mathbf{H}(s)ds\right\}\Phi.$$
 (5.5b)

We define our experimental evolution operator  $\mathbf{U}_{_{1}}[t,0]\Phi$  by

$$\mathbf{U}_{\lambda}[t,0]\Phi = \sum_{n=0}^{\infty} \frac{(\lambda t)^n}{n!} \exp\{-\lambda t\} \mathbf{U}_n[t,0]\Phi.$$
 (5.4b)

We now have the following result, which is a consequence of the fact that Borel summability is regular.

**Theorem 5.4** Assume that the conditions for Theorem 3.5 are satisfied. Then

$$\lim_{\lambda \to \infty} \overline{\mathbf{U}}_{\lambda}[t, 0] \Phi = \lim_{\lambda \to \infty} \mathbf{U}_{\lambda}[t, 0] \Phi = \mathbf{U}[t, 0] \Phi. \tag{5.7}$$

Since  $\lambda \to \infty \Rightarrow \lambda^{-1} \to 0$ , this means that the average time between bubbles is zero (in the limit) so that we get a continuous path. It should be observed that

this continuous path arises from averaging the sum over an infinite number of (discrete) paths. The first term in (5.4b) corresponds to the path of a  $\pi$ -meson that created no bubbles (i.e., the photograph is blank). This event has probability  $\exp\{-\lambda t\}$  (which approaches zero as  $\lambda \to \infty$ ). The n-th term corresponds to the path of a  $\pi$ -meson that created n bubbles, (with probability  $[(\lambda t)^n/n!]\exp\{-\lambda t\}$ ) etc. Before deriving a physical relationship, let  $P[t;s,\lambda]=0$  if  $s \le 0$  and, for  $0 < s < \infty$ , define it as:

$$P[t;s,\lambda] = e^{-\lambda t} \sum_{k=0}^{\lceil \lambda s \rceil} \frac{(\lambda t)^k}{k!},$$
(5.8)

where  $n = \lceil \lambda s \rceil$  is the greatest integer  $\leq \lambda s$ . We can now write  $\mathbf{U}[t,0]\Phi$  as

$$\mathbf{U}[t,0]\Phi = \lim_{\lambda \to \infty} \int_0^\infty d_s P[t;s,\lambda] \mathbf{U}_{\lceil \lambda s \rceil}[s,0]\Phi,$$

$$\mathbf{U}_{\lceil \lambda s \rceil}[s,0]\Phi = \exp\left\{\sum_{j=1}^{\lceil \lambda s \rceil} \int_{t_{j-1}}^{t_j} E[\tau_j,u] \mathbf{H}(u) du\right\}\Phi$$
(5.9)

Equation (5.9) means that we get both a sum over paths and a probability interpretation for our formalism. This allows us to give a new definition for path integrals.

Suppose the evolution operator  $\mathbf{U}[t,0]$  has a kernel,  $\mathbf{K}[\mathbf{x}(t), t; \mathbf{x}(0), 0]$ , such that

1. 
$$\mathbf{K}[\mathbf{x}(t), t; \mathbf{x}(s), s] = \int_{\mathbf{R}^3} \mathbf{K}[\mathbf{x}(t), t; \mathbf{x}(s), s] \mathbf{K}[\mathbf{x}(s), s; \mathbf{x}(0), 0] d\mathbf{x}(s)$$
, and

2. 
$$\mathbf{U}[t,0]\Phi = \int_{\mathbf{R}^3} \mathbf{K}[\mathbf{x}(t), t; \mathbf{x}(0), 0]d\mathbf{x}(0).$$

Then, from equation (5.9), we have that

$$\mathbf{U}[t,0]\Phi = \lim_{\lambda \to \infty} \int_0^\infty d_s P[t;s,\lambda] \left\{ \prod_{j=1}^{\lambda s} \int_{\mathbf{R}^3} \mathbf{K} \left[ \mathbf{x}(t_j), \ t_j; \ \mathbf{x}(t_{j-1}), \ t_{j-1} \right]_j \prod_{j=1}^{\lambda s} d\mathbf{x}(t_{j-1}) \Phi(0) \right\}.$$

Thus, whenever we can associate a kernel with our evolution operator, the time-ordered version always provides a well-defined path-integral as a sum over paths. The definition does not (directly) depend on the space of continuous paths and is independent of a theory of measure on infinite dimensional spaces. Feynman suggested that the operator calculus was more general, in his book with Hibbs90 (see pg. 355-6).

#### **6. The S-Matrix**

The objective of this section is to provide a formulation of the S-matrix that will allow us to investigate the sense in which we can believe Dyson's first conjecture. At the end of his second paper on the relationship between the Feynman and Schwinger-Tomonaga theories, he explored the difference between the divergent Hamiltonian formalism that one must begin with and the finite S-matrix that results from renormalization. He takes the view that it is a contrast between a real observer and a fictitious (ideal) observer. The real observer can only determine particle positions with limited accuracy and always gets finite results from his measurements. Dyson then suggests that "... The ideal observer, however, using non-atomic apparatus whose location in space and time is known with infinite precision, is imagined to be able to disentangle a single field from its interactions with others, and to measure the interaction. In conformity with the Heisenberg uncertainty principle, it can perhaps be considered a physical consequence of the infinitely precise knowledge of (particle) location allowed to the ideal observer, that the value obtained when he measures (the interaction) is infinite." He goes on to

remark that if his analysis is correct, the problem of divergences is attributable to an idealized concept of measurability.

In order to explore this idea, we work in the interaction representation with obvious notation. Replace the interval [t, 0] by [T, -T],  $\mathbf{H}(t)$  by  $(-i/\hbar)\mathbf{H}_{\mathrm{I}}(t)$ , and our experimental evolution operator  $\mathbf{U}_{\lambda}[T, -T]\Phi$  by the experimental scattering operator  $\mathbf{S}_{\lambda}[T, -T]\Phi$ , where

$$\mathbf{S}_{\lambda}[T,-T]\Phi = \sum_{n=0}^{\infty} \frac{(2\lambda T)^n}{n!} \exp[-2\lambda T] \mathbf{S}_n[T,-T]\Phi,$$
(6.1)

$$\mathbf{S}_{n}[T,-T]\Phi = \exp\left\{ (-i/\hbar) \sum_{j=1}^{n} \int_{t_{j-1}}^{t_{j}} E[\tau_{j},s] \mathbf{H}_{I}(s) ds \right\} \Phi, \tag{6.2}$$

and  $\mathbf{H}_{\mathrm{I}}(t) = \int_{\mathbf{R}^3} \mathbf{H}_{\mathrm{I}}(\mathbf{x}(t), t) d\mathbf{x}(t)$  is the interaction energy. We follow Dyson for consistency (see also the discussion), so that  $\delta mc^2$  is the mass counter-term designed to cancel the self-energy divergence, and

$$\mathbf{H}_{\mathrm{I}}(\mathbf{x}(t),t) = -ie\mathbf{A}_{\mu}(\mathbf{x}(t),t)\overline{\psi}(\mathbf{x}(t),t)\gamma_{\mu}\psi(\mathbf{x}(t),t) - \delta mc^{2}\overline{\psi}(\mathbf{x}(t),t)\psi(\mathbf{x}(t),t). \quad (6.3)$$

We now give a physical interpretation of our formalism. Rewrite equation (6.1) as

$$\mathbf{S}_{\lambda}[T,-T]\Phi = \sum_{n=0}^{\infty} \frac{(2\lambda T)^{n}}{n!} \exp\left\{ (-i/\hbar) \sum_{j=1}^{n} \int_{t_{j-1}}^{t_{j}} \left[ E[\tau_{j},s] \mathbf{H}_{I}(s) - i\lambda \hbar \mathbf{I}_{\otimes} \right] ds \right\} \Phi.$$
 (6.4)

In this form, it is clear that the term  $-i\lambda\hbar I_{\otimes}$  has a physical interpretation as the absorption of photon energy of amount  $\lambda\hbar$  in each subinterval  $[t_j,t_{j-1}]$  (cf. Mott and Massey<sup>92</sup>). When we compute the limit, we get the standard S-matrix (on

[T, -T]). It follows that we must add an infinite amount of photon energy to the mathematical description of the experimental picture (at each point in time) in order to obtain the standard scattering operator. This is the ultraviolet divergence and shows explicitly that the transition from the experimental to the ideal scattering operator requires that we illuminate the particle throughout its entire path. Thus, it appears that we have, indeed, violated the uncertainty relation. This is further supported if we look at the form of the standard Smatrix:

$$\mathbf{S}[T, -T]\Phi = \exp\left\{ (-i/\hbar) \int_{-T}^{T} \mathbf{H}_{I}(s) ds \right\} \Phi, \tag{6.5}$$

and note that the differential ds in the exponent implies perfect infinitesimal time knowledge at each point, strongly suggesting that the energy should be totally undetermined. If violation of the Heisenberg uncertainty relation is the cause for the ultraviolet divergence, then as it is a variance relation, it will not appear in first order (perturbation) but should show up in all higher-order terms. On the other hand, if we eliminate the divergent terms in second order, we would expect our method to prevent them from appearing in any higher order term of the expansion. The fact that this is precisely the case in quantum electrodynamics is a clear verification of Dyson's conjecture.

If we allow T to become infinite, we once again introduce an infinite amount of energy into the mathematical description of the experimental picture, as this is also equivalent to precise time knowledge (at infinity). Of course, this is the well-known infrared divergence and can be eliminated by keeping T finite (see Dahmen et al<sup>93</sup>) or introducing a small mass for the photon (see Feynman<sup>12</sup>, pg. 769). If we hold  $\lambda$  fixed while letting T become infinite, the experimental S-matrix takes the form:

$$\mathbf{S}_{\lambda}[\infty, -\infty] \Phi = \exp \left\{ (-i/\hbar) \sum_{j=1}^{\infty} \int_{t_{j-1}}^{t_{j}} E[\tau_{j}, s] \mathbf{H}_{I}(s) ds \right\} \Phi,$$

$$\bigcup_{j=1}^{\infty} \left[ t_{j-1}, t_{j} \right] = (-\infty, \infty), \quad \& \quad \Delta t_{j} = \lambda^{-1}.$$

$$(6.6)$$

This form is interesting since it shows how a minimal time eliminates the ultraviolet divergence. Of course, this is not unexpected, and has been known at least since Heisenberg<sup>94</sup> introduced his fundamental length as a way around the divergences. This was a prelude to the various lattice approximation methods. The review by Lee<sup>95</sup> is interesting in this regard.

In closing this section, we record our exact expansion for the S-matrix to any finite order. With  $\Phi(-\infty) \in D[(\mathbf{Q}[\infty, -\infty])^{n+1}]$ , we have

$$\mathbf{S}[\infty, -\infty]\Phi(-\infty) = \sum_{k=0}^{n} \left(\frac{-i}{\hbar}\right)^{k} \int_{-\infty}^{\infty} ds_{1} \int_{-\infty}^{s_{1}} ds_{2} \cdots \int_{-\infty}^{s_{k-1}} ds_{k} \mathbf{H}_{1}(s_{1}) \mathbf{H}_{1}(s_{2}) \cdots \mathbf{H}_{1}(s_{k}) \Phi(-\infty)$$

$$+ \left(\frac{-i}{\hbar}\right)^{n+1} \int_{0}^{1} (1-\xi)^{n} d\xi \int_{-\infty}^{\infty} ds_{1} \int_{-\infty}^{s_{1}} ds_{2} \cdots \int_{-\infty}^{s_{n}} ds_{n+1} \mathbf{H}_{1}(s_{1}) \mathbf{H}_{1}(s_{2}) \cdots \mathbf{H}_{1}(s_{n+1}) \mathbf{S}^{\xi}[s_{n+1}, -\infty] \Phi(-\infty).$$

$$(6.7)$$

It follows that (in a theoretical sense) we can consider the standard S-matrix expansion to be exact, when truncated at any order, by adding the last term of equation (6.7) to give the remainder. This result also means that whenever we can construct an exact nonperturbative solution, it always implies the

existence of a perturbative solution valid to any order. However, in general, only in particular cases can we know if the series at some n (without the remainder) approximates the solution.

#### **Discussion**

In this paper we have shown how to construct a natural representation Hilbert space for Feynman's time-ordered operator calculus. This space allows us to construct the time-ordered integral and evolution operator (propagator) under the weakest known conditions. Using the theory, we have shown that the perturbation expansion relevant to quantum theory is asymptotic in the sense of Poincaré. This provides a precise formulation and proof of Dyson's second conjecture16 that, in general, we can only expect the expansion to be asymptotic.

Our investigation into the extent that our continuous models for the micro-world faithfully represent the amount of information available from experiment has led to a derivation of the time-ordered evolution operator in a more physical way. This approach made it possible to prove that the ultraviolet divergence is caused by a violation of the Heisenberg uncertainty relation at each point in time, thus partially confirming Dyson's first conjecture.

We used Dyson's original notation so as to explicitly exhibit the counter-term necessary to eliminate the self-energy divergence that occurs in QED. This divergence is not accounted for and is outside the scope of the current investigation. Thus, within our present framework, we cannot say that all the divergences arise from our disregard of some simple physics, and are not the result of deeper problems. Thus, Dyson's concerns about the mathematical consistency of quantum electrodynamics, and quantum field theory in general, is still an open problem.

Although we are not working in the framework of axiomatic field theory, our approach may make some uneasy since Haag's theorem suggests that the interaction representation does not exist (see Streater and Wightman27 pg. 161). (Haag's theorem assumes, among other things, that the equal time commutation relations for the canonical variables of a interacting field are equivalent to those of a free field.) In trying to explain this unfortunate result, these authors point out that ( see pg. 168) "… What is even more likely in physically interesting quantum field theories is that equal time commutation relations will make no sense at all; the field might not be an operator unless smeared in time as well as space. " The work in Sections 5 and 6 of this paper strongly suggests that there is no physical basis to assume that we know anything about canonical variables at one instant in time (see postulate 2 and the following paragraph). Thus, our approach actually confirms the above comments of Streater and Wightman.

## **Acknowledgments**

Work for this paper was begun while the first author was supported as a member of the School of Mathematics in the Institute for Advanced Study, Princeton, NJ, and completed while visiting in the physics department of the University of Michigan.

#### **References**

1 P. A. M. Dirac, *Proc. Roy. Soc. London* **A114** (1927), 243.

2 P. A. M. Dirac, *Proc. Roy. Soc. London* **A117** (1928), 610.

3 W. Heisenberg and W. Pauli, *Z. Phys.* **56** (1929), 1.

4 W. Heisenberg and W. Pauli, *Z. Phys.* **59** (1930), 168.

5 A. I. Miller, *Early quantum electrodynamics*, Cambridge University Press, Cambridge, UK, 1994.

6 S. S. Schweber, *QED and the Men Who Made It: Dyson, Feynman, Schwinger, and Tomonaga*, Princeton University Press, Princeton, NJ 1994.

7 J. R. Oppenheimer, *Physical Review* **35** (1930), 461.

8 J. R. Oppenheimer, Rapports du 8e Conseil de Physique, Solvay, (1950) 269.

9 P. A. M. Dirac, *Proc. Roy. Soc. London* **A167** (1938), 148.

10S. Tomonaga, *Prog. Theor. Phys.* **1** (1949), 27.

11J. Schwinger, *Physical Review* **73** (1948), 416.

- 12R. P. Feynman, *Physical Review* **76** (1949), 749, 769.
- 13R. P. Feynman, *Physical Review* **80** (1950), 440.
- 14J. Schwinger, *Selected Papers on Quantum Electrodynamics*, Dover Publications, New York 1958.
- 15F. J. Dyson, *Physical Review* **75** (1949), 486, 1736.
- 16F. J. Dyson, *Physical Review* **85** (1952), 631.
- 17R. P. Feynman, *Physical Review* **81** (1951), 108.
- 18P. A. M. Dirac, in *The Birth of Particle Physics*, edited by L. M. Brown and L. Hodderson, Cambridge University Press, Cambridge, UK, 1983, 39.
- 19S. Sakata, *Prog. Theor. Phys.* **16** (1956), 686.
- 20J. Schwinger, in *The Birth of Particle Physics*, edited by L. M. Brown and L. Hodderson, Cambridge University Press, Cambridge, UK, 1983, 329.
- 21A. S. Wightman, Phys. Rev. **101** (1956), 860.
- 22H. Lehmann, K. Symanzik, and W. Zimmermann, *Nuovo Cim.* **1** (1955), 205.
- 23H. Lehmann, K. Symanzik, and W. Zimmermann, *Nuovo Cim.* **6** (1957), 319.
- 24A. S. Wightman, *Fortschr. Phys.* **44** (1996), 143.
- 25W. Heisenberg, Sächs. Akad. Wiss. (Leipzig), Berichte d. math.-phys. KI. **83** (1931), 3.

- 26R. Jost, *The General Theory of Quantized Fields*, Amer. Math. Soc. Providence, RI, 1965.
- 27R. F. Streater and A. S. Wightman, *PCT, Spin and statistics and all that*, Benjamin, New York, 1964.
- 28N. N. Bogolubov and D. V. Shirkov, *Introduction to the theory of quantized fields*, Interscience, London, 1958.
- 29R. Haag, *Local Quantum Physics*, (2nd edition) Springer, New York, 1996.
- 30N. N. Bogolubov, A. A. Logunov and I. T. Todorov, *Introduction to axiomatic quantum field theory*, Benjamin, New York, 1975.
- 31D. Buchholz, *Current Trends in Axiomatic Quantum Field Theory*, (preprint) hep-th/9811233, 1998.
- 32G. Velo and A. S. Wightman (editors), *Constructive Quantum Field Theory*, Lecture Notes in Phys. **25**, Springer, New York, 1974.
- 33G. Velo and A. S. Wightman (editors), *Constructive Quantum Field Theory II*, NATO ASI Series B, Vol. 234, 1988.
- 34J. Glimm and A. Jaffe, *Quantum Physics. A functional integral point of view*, Springer, New York, 1987.
- 35B. Simon, *Functional Integration and Quantum Physics*, Academic Press, New York, 1979.
- 36A. D. Sokal, *Ann. Inst. Henri Poincaré*, **A37** (1982), 13.
- 37M. Aizenman and R. Graham, *Nucl. Phys.* **B225** (1983), 261.

- 38J. Fröhlich, *Nucl. Phys.* **B200** (1982), 281.
- 39K. Gawedzki and A. Kupiainen, *Commun. Math. Phys.* **102** (1985), 1.
- 40A. S. Wightman in *Renormalization Theory*, NATO ASI Series C, Vol. 23, edited by G. Velo and A. S. Wightman, D. Reidel Publishers, Dordrecht, Holland, 1976.
- 41O. I. Zavialov, *Renormalized Quantum Field Theory*, Kluwer, Dordrecht, Holland, 1990.
- 42J. M. Jauch and F. Rohrlich, *The Theory of Photons and Electrons,* Addison-Wesley, Reading, MA 1955.
- 4 3 J. Feldman, T. Hurd, L. Rosen and J. Wright, *QED: A Proof of Renormalizability*, Springer Lecture Notes in Physics, **312**, Springer, New York, 1988.
- 44K. G. Wilson, *Phys. Rev.* **D6** (1972), 419.
- 45G. Gallavotti, *Rev. Mod. Phys.* **57** (1985), 471.
- 46K. Hepp, *Commun. Math. Phys.* **2** (1966), 301.
- 47J. Lowenstein and E. Speer, *Commun. Math. Phys.* **47** (1976), 4.
- 48A. Salam, *Phys. Rev.* **84** (1951), 426.
- 49J. Ward, *Proc. Phys. Soc. London*, **A64** (1951), 54.
- 50R. L. Mills and C. N. Yang, *Prog. Theoret. Phys. Supp.* **37** (1966), 507.

- 51K. Hepp, *Théorie de la Rénormalization,* Springer Lecture Notes in Physics, Vol. **2** Springer, New York, 1969.
- 52S. Weinberg, *Phys. Rev.* **118** (1960), 838.
- 53I. E. Segal, in *Lectures in Modern Analysis and Applications II*, (G. T. Taam, editor.) Lectures Notes in Math. **140**, Springer, New York, 1970.
- 54W. L. Miranker, and B. Weiss, *SIAM Rev.* **6** (1966), 104.
- 55E. Nelson, in *Functional Analysis and Related Fields* (F. Browder, ed.) Springer-Verlag New York, 1970.
- 56H. Araki, *Ann. Sci. Ecole Norm. Sup.* **6** (1973), 67.
- 57I. Fujiwara, *Prog. Theor. Phys.* **7** (1952), 433.
- 58V. P. Maslov, *Operational Methods*, Mir, Moscow, 1976.
- 59G. W. Johnson and M. L. Lapidus, *Mem. Amer. Math. Soc.* **62** (1986), 1.
- 60G. W. Johnson and M. L. Lapidus, *J. Funct. Anal.* **81** (1988), 74.
- 61G. W. Johnson and M. L. Lapidus, *The Feynman Integral and Feynman's Operational Calculus,* Oxford U. Press, New York, 2000.
- 62G. W. Johnson, M. L. Lapidus, and B. DeFacio, in *Stochastics Processes: A Festschrift in Honour of Gopinath Kallianpur*, (S. Cambanis, J. K. Ghosh, R. L. Karandikar and P. K. Sen, editors) Springer-Verlag, New York, 1993.
- 63T. L. Gill, *Transactions of the American Mathematical Society* **266** (1981), 161.

- 64T. L. Gill, *Transactions of the American Mathematical Society* **279** (1983), 617.
- 65T. L. Gill, and W.W. Zachary, *J. of Math. Phys.* **28** (1987), 1459.
- 66R. Henstock, *Theory of Integration*, Butterworth, London, 1963.
- 67J. Kurzweil, *Czechoslovak Math. J.* **7** (1957), 418.
- 68R. Henstock, *Proc. London Math. Soc.* **27** (1973), 317.
- 69P. Muldowney, *A General Theory of Integration in Function Spaces*, Pitman Research Notes in Mathematics, John Wiley & Sons, New York, 1987.
- 70J. A. Goldstein, *Semigroups of Linear Operators and Applications,* Oxford U. Press, New York, 1985.
- 71A. Pazy, *Semigroups of Linear Operators and Applications to Partial Differential Equations,* Appl. Math. Sci. **44**, Springer-Verlag, New York, 1983.
- 72J. von Neumann, *Compositio Math.* **6** (1938), 1.
- 73T. L. Gill, *J. Funct. Anal.* **30** (1978), 17.
- 74E. Hille, and R. S. Phillips, *Functional Analysis and Semigroups,* Amer. Math. Soc. Colloq. Pub. **31**, Amer. Math. Soc. Providence, RI, 1957.
- 75C. A. Hurst, *Proc. Camb. Phil*. Soc. **48** (1952), 625.
- 76W. Thirring, *Helv. Phys. Acta* **26** (1953), 33.
- 77A. Petermann, *Helv. Phys. Acta* **26** (1953), 291.

- 78A. Jaffe, *Commun. Math. Phys.* **1** (1965), 127.
- 79J. Zinn-Justin, *Quantum Field Theory and Critical Phenomena,* 2nd Ed. Clarendon Press, Oxford, 1993.
- 80S. S. Schweber, *An Introduction To Relativistic Quantum Field Theory,* Harper & Row, New York, 1961.
- 81J. D. Dollard and C. N. Friedman, *Product Integration with Applications to Differential Equations,* Encyclopedia of Math. **10**, Addison-Wesley, Reading Mass.,1979.
- 82T. Tang, and Z. Li, *Lett. Math. Phys.* **43** (1998), 55.
- 83K.–J. Engel and R. Nagel, *One-Parameter Semigroups for Linear Evolution Equations*, Grad. Texts in Math., Springer, New York, 2000.
- 84A. S. Wightman, in *The Lesson of Quantum Theory*, (edited by J. de Boer, E. Dal and O. Ulfbeck) Elsevier, Amsterdam, 1986, pg. 201.
- 85A. Einstein, *Ann. Phys.* (Leipzig) **17** (1905), 549.
- 86N. Wiener, A. Siegel, B. Rankin, and W. T. Martin, *Differential Space, Quantum Systems, and Prediction,* M. I. T. Press Cambridge, MA, 1966.
- 87L. S. Ornstein and G. E. Uhlenbeck, *Physical Review* **36** (1930), 1103.
- 88J. L. Doob, *Ann. Math.* **43** (1954), 62.
- 89S. M. Ross, *Introduction To Probability Models*, (4th. edition) Academic Press, New York, 1989.

- 90R. P. Feynman and A. R. Hibbs, *Quantum Mechanics and Path Integrals*, McGraw-Hill, New York, 1965.
- 91F. Dyson, *Selected Papers of Freeman Dyson with Commentary*, Amer. Math. Soc. Providence, RI, 1996.
- 92N. F. Mott, H. S. Massey, *Theory of Atomic Collisions*, Clarendon Press, Oxford, U. K, 1965.
- 93H. D. Dahmen, B. Scholz, and F. Steiner, *Nuclear Phys.* **B202** (1982), 365.
- 94W. Heisenberg, *Ann. Phys.* (Leipzig) **32** (1938), 20.
- 95T. D. Lee, in *The Lesson of Quantum Theory*, (edited by J. de Boer, E. Dal, and O. Ulfbeck) Elsevier, Amsterdam, 1986, pg. 181.