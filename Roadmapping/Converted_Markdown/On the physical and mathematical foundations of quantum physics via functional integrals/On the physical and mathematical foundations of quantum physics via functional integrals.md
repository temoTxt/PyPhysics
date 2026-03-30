## ON THE PHYSICAL AND MATHEMATICAL FOUNDATIONS OF QUANTUM PHYSICS VIA FUNCTIONAL INTEGRALS

#### GIAMPIERO ESPOSITO AND TEPPER L. GILL

Abstract. In order to preserve the leading role of the action principle in formulating all field theories one needs quantum field theory, with the associated BRST symmetry, and Feynman-DeWitt-Faddeev-Popov ghost fields. Such fields result from the fibre-bundle structure of the space of histories, but the physics-oriented literature used them formally because a rigorous theory of measure and integration was lacking. Motivated by this framework, this paper exploits previous work of Gill and Zachary, where the use of Banach spaces for the Feynman integral was proposed. The Henstock-Kurzweil integral is first introduced, because it makes it possible to integrate functions like exp{ix2/2}. The Lebesgue measure on R<sup>∞</sup> is then built and used to define the measure on every separable Hilbert space. The subsequent step is the construction of a new Hilbert space KS<sup>2</sup> [Rn], which contains L<sup>2</sup> [Rn] as a continuous dense embedding, and contains both the test functions D[Rn] and their dual D∗[Rn], the Schwartz space of distributions, as continuous embeddings. This space allows us to construct the Feynman path integral in a manner that maintains its intuitive and computational advantages. We also extend this space to KS<sup>2</sup> [H], where H is any separable Hilbert space. Last, the existence of a unique universal definition of time, τh, that we call historical time, is proven. We use τ<sup>h</sup> as the order parameter for our construction of Feynman's time ordered operator calculus, which in turn is used to extend the path integral in order to include all time dependent groups and semigroups with a reproducing kernel representation.

## 1. Introduction

The present paper is devoted to the physical and mathematical foundations of Feynman's functional integrals. On the physical side, one should first acknowledge the early work of Dirac, who studied, first, the Lagrangian in quantum mechanics both in a paper [DI1] and in his famous book on the principles of quantum mechanics [DI2]. Later on Feynman, who had studied the action at a distance with Wheeler, became interested in a model where two particles can interact with each other [WH], but cannot interact with themselves [FY1]. Feynman realized that he could not write down Hamilton's equations of motion for these particles, and hence was led to study a quantization of the two-particle system not relying upon Hamiltonian methods, but rather using Lagrangian methods along the seminal work by Dirac. Over more than two decades Feynman developed first what he called a spacetime approach to nonrelativistic quantum mechanics [FY2], and then extended formally his approach to gauge field theories [FY3, FY5], discovering eventually the fields [FY4] that came to be known as Feynman-DeWitt-Faddeev-Popov ghost fields. After the pioneering work of Misner [MI], Feynman presented his idea at a conference in Poland [FY4]. Within a few years, Faddeev and Popov had developed a formal derivation for quantum Yang-Mills theory [FAD], whereas DeWitt obtained a systematic derivation for all gauge field theories [DW1], showing that ghost fields arise from the fibre-bundle structure of the space of field histories, when a spacetime [DW3], global approach to quantum field theory [DW4] is developed. With hindsight, the need for a quantum theory of gauge fields might be described as follows.

When we impose a supplementary condition in field theory [DW4] in order to obtain an invertible operator on gauge fields, e.g., the Lorenz condition in electrodynamics

(1) 
$$\Phi_L(A) = \sum_{\mu=1}^4 \nabla^{\mu} A_{\mu} = \tau,$$

or the de Donder condition for gravitation:

(2) 
$$\Phi_{\mu}(h) = \sum_{\nu=1}^{4} \nabla^{\nu} \left( h_{\mu\nu} - \frac{1}{2} g_{\mu\nu} \operatorname{tr}(h) \right) = \tau_{\mu},$$

where hµν are the components of perturbations of the metric g, we are dealing with functionals that associate to physical fields some equations which can be denoted by P <sup>α</sup>, the α being Lie-algebra indices. At classical level one can remark for example that the Maxwell action functional in Minkowski spacetime with metric η provides the noninvertible operator

$$Q_{\mu\nu} = -\eta_{\mu\nu}\Box + \partial_{\mu}\partial_{\nu}.$$

Such an operator acts on smooth sections of the bundle of 1-form fields on Minkowski spacetime, as is clearer by writing it in the form

$$Q_{\mu}^{\ \nu} = -\delta_{\mu}^{\ \nu} \Box + \partial_{\mu} \partial^{\nu}.$$

We note now that an additional term in the action such as the square of the Lorenz gauge would provide a contribution −∂µ∂<sup>ν</sup> to such an operator, by virtue of the identity

(3) 
$$\Phi_L^2(A) = \sum_{\mu,\nu=1}^4 \left[ \partial^{\mu} (A_{\mu} \partial^{\nu} A_{\nu}) - A_{\mu} \partial^{\mu} \partial^{\nu} A_{\nu} \right].$$

The resulting operator on A<sup>ν</sup> would then be Pµν = −ηµν✷, which is instead invertible. More generally, any operator of the form (here ρ ∈ R − {0})

(4) 
$$P_{\mu\nu} = -\eta_{\mu\nu}\Box + (1 - 1/\rho)\,\partial_{\mu}\partial_{\nu}$$

is invertible because its symbol is the matrix (here k <sup>2</sup> = P<sup>4</sup> µ,ν=1 ηµνk µk ν )

(5) 
$$\sigma_{\mu\nu} = \eta_{\mu\nu}k^2 - (1 - 1/\rho)k_{\mu}k_{\nu},$$

whose inverse is [E]

(6) 
$$\Sigma^{\nu\lambda} = \frac{1}{k^2} \eta^{\nu\lambda} + (\rho - 1) \frac{k^{\nu} k^{\lambda}}{k^4}.$$

The expression of the additional term in the action for arbitrary choice of supplementary condition is therefore

$$\frac{1}{2} \sum_{\alpha,\beta} \int d^4x \int d^4x' \ P^{\alpha}(x) \omega_{\alpha\beta}(x,x') P^{\beta}(x') = \frac{1}{2} P^{\alpha} \omega_{\alpha\beta'} P^{\beta'}.$$

With this notation, repeated indices without summation symbol represent actually summation as well as integration [DW3, DW4].

Moreover, from the identity (Latin lower case indices being used for fields)

(7) 
$$\delta P^{\alpha} = P^{\alpha}_{,i} \delta \varphi^{i} = P^{\alpha}_{,i} Q^{i}_{\beta} \delta \xi^{\beta} = \widehat{F}^{\alpha}_{\beta} \delta \xi^{\beta}$$

we discover the existence of the operator defined as follows:

(8) 
$$\widehat{F}^{\alpha}_{\beta}[\varphi] \equiv P^{\alpha}_{,i}[\varphi]Q^{i}_{\beta}[\varphi].$$

Such an operator should act upon fields χ<sup>α</sup> and ψ β , in general independent of each other. For example, in the case of electrodynamics, <sup>F</sup><sup>b</sup> is a map from the space of smooth functions on spacetime into itself:

$$\widehat{F}: C^{\infty}(M,g) \to C^{\infty}(M,g).$$

ON THE PHYSICAL AND MATHEMATICAL FOUNDATIONS OF QUANTUM PHYSICS VIA FUNCTIONAL INTEGRALS5 For gravitation, <sup>F</sup><sup>b</sup> maps smooth vector fields <sup>σ</sup><sup>1</sup> over spacetime into smooth vector fields σ<sup>2</sup> over spacetime (T(M, g) being our notation for the tangent bundle of the spacetime manifold (M, g)):

$$\sigma_1:(M,g)\longrightarrow T(M,g),\ \sigma_2:(M,g)\longrightarrow T(M,g),\ \widehat{F}:\sigma_1\longrightarrow \sigma_2.$$

In light of Eq. (8) we are led to assume that one can build the action functional

(9) 
$$\widetilde{S}[\varphi, \chi, \psi] = S[\varphi] + \frac{1}{2} P^{\alpha}[\varphi] \omega_{\alpha\beta'} P^{\beta'}[\varphi] + \chi_{\alpha} \widehat{F}^{\alpha}_{\beta'}[\varphi] \psi^{\beta'}.$$

The term quadratic in the functionals P <sup>α</sup> spoils gauge invariance of the classical action, whereas the sum of the three terms in Eq. (9) has a new invariance property (see theorem below) under a class of transformations which are called [BE] Becchi-Rouet-Stora-Tyutin (hereafter BRST) transformation. They can be written in the form [DW2]

(10) 
$$\delta \varphi^i = Q^i_{\alpha}[\varphi] \psi^{\alpha} \delta \lambda,$$

(11) 
$$\delta \chi_{\alpha} = \omega_{\alpha\beta} P^{\beta}[\varphi] \delta \lambda,$$

(12) 
$$\delta \psi^{\alpha} = -\frac{1}{2} C^{\alpha}_{\beta \gamma} \psi^{\beta} \psi^{\gamma} \delta \lambda,$$

where δλ is a constant which commutes with ϕ <sup>i</sup> and anticommutes with χ<sup>α</sup> and ψ <sup>α</sup>. The Q<sup>i</sup> <sup>α</sup>[ϕ] are the generators of infinitesimal gauge transformations (we write δ<sup>G</sup> in order to avoid confusion with δ used for BRST in Eqs. (10)-(12), and we write δGξ <sup>α</sup> to denote a set of linearly independent group parameters):

(13) 
$$\delta_G \varphi^i = Q^i_{\alpha}[\varphi] \delta_G \xi^{\alpha}.$$

Such Q<sup>i</sup> <sup>α</sup> are restricted from the identity of group-theoretical nature [DW4]

(14) 
$$Q^{i}_{\alpha,j}[\varphi]Q^{j}_{\beta}[\varphi] - Q^{i}_{\beta,j}[\varphi]Q^{j}_{\alpha}[\varphi] = Q^{i}_{\gamma}[\varphi]C^{\gamma}_{\alpha\beta},$$

where C α βγ are the structure constants of the Lie algebra of the gauge group. For local theories, the Q<sup>i</sup> <sup>α</sup> are linear combinations of Dirac's δ and its derivatives. The ghost operator is discovered from Eq. (7), which shows how the functionals P α used to fix the supplementary conditions are varying under gauge transformations. The classical roots of such a statement can be understood by remarking for example that, under a gauge transformation of the potential, the 1-form A, the functional Φ<sup>L</sup> of the Lorenz gauge changes by the amount (on denoting by ε a function of class C 2 in the gauge transformation of the components of A)

(15) 
$$\Phi_L(A) - \Phi_L(A + \nabla \varepsilon) = -\Box \varepsilon.$$

Similarly, the functional (2) for the de Donder gauge varies under infinitesimal diffeomorphisms according to the relation

(16) 
$$\Phi_{\mu}(h) - \Phi_{\mu}(h + L_X g) = -\frac{1}{2} \sum_{\nu=1}^{4} (g_{\mu\nu} \Box + R_{\mu\nu}) X^{\nu}.$$

The second-order differential operators in such equations are two different realizations of one and the same concept, i.e., the ghost operator [DW3, DW4] defined in Eq. (8). Now we can state and prove a well-known key theorem, as follows.

Theorem on the existence of a BRST-invariant action functional. The action functional defined in Eq. (9) is invariant under the BRST transformations defined in Eqs. (10), (11) and (12).

**Proof.** First of all let us point out that, by virtue of Eqs. (10)-(12), the infinitesimal BRST variation of  $\widetilde{S}[\varphi, \chi, \psi]$  takes the form

$$\delta \widetilde{S}[\varphi, \chi, \psi] = S_{,i}[\varphi]Q^{i}_{\alpha}[\varphi]\psi^{\alpha}\delta\lambda + \frac{1}{2}P^{\alpha}_{,i}[\varphi]Q^{i}_{\gamma}[\varphi]\omega_{\alpha\beta}P^{\beta}[\varphi]\psi^{\gamma}\delta\lambda$$

$$+ \frac{1}{2}P^{\beta}[\varphi]\omega_{\beta\alpha}P^{\alpha}_{,i}[\varphi]Q^{i}_{\gamma}[\varphi]\psi^{\gamma}\delta\lambda + \omega_{\alpha\beta}P^{\beta}[\varphi]\delta\lambda P^{\alpha}_{,i}[\varphi]Q^{i}_{\gamma}[\varphi]\psi^{\gamma}$$

$$+ \chi_{\alpha}P^{\alpha}_{,ij}[\varphi]Q^{j}_{\gamma}[\varphi]\psi^{\gamma}\delta\lambda Q^{i}_{\beta}[\varphi]\psi^{\beta}$$

$$+ \chi_{\alpha}P^{\alpha}_{,i}[\varphi]Q^{i}_{\beta,j}[\varphi]Q^{j}_{\gamma}[\varphi]\psi^{\gamma}\delta\lambda\psi^{\beta}$$

$$(17) \qquad - \frac{1}{2}\chi_{\alpha}P^{\alpha}_{,i}[\varphi]Q^{i}_{\beta}[\varphi]C^{\beta}_{\gamma\delta}\psi^{\gamma}\psi^{\delta}\delta\lambda.$$

From the gauge invariance of the classical action one finds that

$$S_{i}[\varphi]Q^{i}_{\gamma}[\varphi] = 0,$$

because

$$\delta_G S = S_{,i}[\varphi] \delta_G \varphi^i = S_{,i}[\varphi] Q^i_{\alpha}[\varphi] \delta_G \xi^{\alpha} = 0.$$

Moreover, the sum of second, third and fourth term on the right-hand side of Eq. (17) vanishes as well, because  $(\delta\lambda)\psi^{\alpha} = -\psi^{\alpha}(\delta\lambda)$ , and exploiting the symmetric nature of  $\omega_{\alpha\beta}$ . The fifth term on the right-hand side of Eq. (17) reduces to

$$-\chi_{\alpha} P^{\alpha}_{,ij}[\varphi] Q^{j}_{(\gamma}[\varphi] Q^{i}_{\beta)}[\varphi] \psi^{[\gamma} \psi^{\beta]} \delta \lambda$$

and hence vanishes as well. Last, upon using the identity (14) in order to express  $Q^{i}_{\beta}[\varphi]C^{\beta}_{\gamma\delta}$ , the sum of sixth and seventh term on the right-hand side of (17) is given by

$$\chi_{\alpha}P^{\alpha}_{,i}[\varphi]Q^{i}_{\beta,j}[\varphi]Q^{j}_{\gamma}[\varphi]\psi^{\beta}\psi^{\gamma}\delta\lambda - \frac{1}{2}\chi_{\alpha}P^{\alpha}_{,i}[\varphi]Q^{i}_{\gamma,j}[\varphi]Q^{j}_{\delta}[\varphi]\psi^{\gamma}\psi^{\delta}\delta\lambda$$
$$+ \frac{1}{2}\chi_{\alpha}P^{\alpha}_{\ i}[\varphi]Q^{i}_{\ \gamma\ j}[\varphi]Q^{j}_{\beta}[\varphi]\psi^{\beta}\psi^{\gamma}\delta\lambda$$

which is found to vanish after relabelling indices and exploiting the identity ψ <sup>β</sup>ψ <sup>γ</sup> = −ψ <sup>γ</sup>ψ β . Q.E.D.

Thus, the BRST transformations involve also anticommuting fields, and hence the BRST symmetry is not classical. Conceptually, an important conclusion is found to emerge [DW4, E]:

- (1) Relativity suggests using the action principle as a foundation for all field theories.
- (2) The gauge principle leads to a gauge-invariant action functional S.
- (3) The resulting operator S,ij on gauge fields (it maps indeed gauge fields into gauge fields) is not invertible, so that no Green function can be defined for S,ij .
- (4) Such a drawback is amended by adding to the action the second term on the right-hand side of Eq. (9).
- (5) The full action becomes the functional <sup>S</sup><sup>e</sup> of Eq. (9), which is no longer gaugeinvariant, but rather BRST-invariant. This is a quantum symmetry, since it needs anticommuting fields associated to bosonic fields.
- (6) With hindsight, quantum field theory turns out to be the branch of physics which is needed in order to preserve the guiding role of the action principle, and the driving force is provided by relativity, since it puts the emphasis on constructing the action upon completing the original gauge-invariant action by the addition of the other terms in Eq. (9).

In quantum field theory, the anticommuting fields discussed before are accommodated within a global approach that relies upon functional integrals. It is therefore

ON THE PHYSICAL AND MATHEMATICAL FOUNDATIONS OF QUANTUM PHYSICS VIA FUNCTIONAL INTEGRALS9 compelling to build a mathematically consistent theory of such integrals [GZ09], which is the topic of our paper.

1.1. Outline of this paper. After some preliminaries, in section 2, we introduce the Henstock-Kurzweil integral (HK). This integral extends those of Bochner and Pettis and is defined for operator-valued functions that are not separately valued (where both the Bochner and Pettis integrals are not defined). Intuitively, one uses a version of the Riemann integral (of calculus) with the interior points chosen first, while the size of the base rectangle around any interior point is determined by an arbitrary positive function defined at that point (see also Appendix). It is important to note that it integrates functions such as exp{ix2/2}. This integral was discovered by Henstock [HS1], [HS2], and Kurzweil [KW].

In Section 3 we construct the Lebesgue measure on R<sup>∞</sup> by using its equivalent space R<sup>∞</sup> I , which is defined in the section. This approach makes it possible for us to define Lebesgue measure on every separable Hilbert space in Section 4.

In Section 5 we construct a new Hilbert space KS<sup>2</sup> [R n I ], which contains L 2 [R n I as a continuous dense embedding, and contains both the test functions D[R n I ] and their dual D<sup>∗</sup> [R n I ], the Schwartz space of distributions, as continuous embeddings. This space allows us to construct the Feynman path integral in a manner that maintains its intuitive [GMG] and computational advantages. We also extend this space to KS<sup>2</sup> [H], where H is a separable Hilbert space.

In Section 6 we prove that a unique universal definition of time, τh, exists (which we call historical time). We use τ<sup>h</sup> as the order parameter for the construction of Feynman's time ordered operator calculus, which in turn is used to extend the path integral to include all time dependent groups and semigroups with a kernel.

Section 7 is devoted to concluding remarks, and some basic concepts of real and functional analysis are further described in the Appendix. It has been our choice to describe modern analysis concepts which, despite being widespead in mathematics, are often ignored in the physics-oriented literature. We think that their potentialities deserve greater attention in theoretical physics.

### 1.2. Preliminaries.

1.2.1. Extensions and Cubes. Before we discuss the Henstock-Kurzweil integral on R <sup>n</sup>, we need an extension theorem for functions defined on a domain of R <sup>n</sup> and a result that shows that a domain in R <sup>n</sup> can be written as a union of nonoverlapping closed cubes. (proofs of these results may be found in Evans [EV] and Stein [STE], respectively.) Let D be a bounded open connected set of R <sup>n</sup> (a domain) with topological boundary ∂D and closure D¯.

Definition 1.1. Let k be a positive integer: we say that ∂D is of class C<sup>k</sup> if, for every point x ∈ ∂D, there is a homeomorphism ϕ of a neighborhood U of x into R n such that both ϕ and ϕ <sup>−</sup><sup>1</sup> have k continuous derivatives with

$$\varphi\left(\mathbb{D}\cap U\right)\subset\left\{\left(x_{1},\ldots,x_{n}\right)\in\mathbb{R}^{n}:x_{n}>0\right\}$$

and

$$\varphi(\partial \mathbb{D} \cap U) \subset \{(x_1, \dots, x_n) \in \mathbb{R}^n : x_n = 0\}.$$

Theorem 1.2. Let D be a domain in R <sup>n</sup> with ∂D of class C<sup>1</sup> . Let U be any bounded open set such that the closure of D is a compact subset of U. Then there exists a linear operator E mapping functions on D to functions on R <sup>n</sup> such that:

- (1) C maps W1,<sup>2</sup> (D) (see below) continuously into W1,<sup>2</sup> (R <sup>n</sup>).
- (2) C(f)|<sup>D</sup> = f (e. g., E(·) is an extension operator).
- (3) E(f)(x) = 0 for x ∈ U c (e.g., E(f) has support inside U).

Theorem 1.3. Let D be a domain in R <sup>n</sup>. Then D is the union of a sequence of closed cubes {Dk} whose sides are parallel to the coordinate axes and whose interiors are mutually disjoint.

Remark 1.4. Thus, if a function f is defined on a domain in R <sup>n</sup>, by Theorem 1.2 it can be extended to the whole space. By Theorem 1.3, without loss of generality, we can assume that the domain is a cube with sides parallel to the coordinate axes. In either case, the HK-integral can be constructed under these conditions.

1.2.2. Distributions. Let D(R <sup>n</sup>) = C<sup>∞</sup> c (R <sup>n</sup>) denote the space of infinitely differentiable functions φ : R <sup>n</sup> → R with compact support. We say that a sequence of functions {φn} ⊂ D converges to φ ∈ D if there is a fixed compact set U such that all functions φ<sup>n</sup> have their support in U and, for each k ≥ 0, the sequence of k-derivatives of φn, φ(k) <sup>n</sup> , converges uniformly [FMS] to φ (k) on U. We call a function φ belonging to D(R <sup>n</sup>) a test function.

Let u ∈ C 1 (R <sup>n</sup>). Then, if φ ∈ C<sup>∞</sup> c (R <sup>n</sup>), integration by parts gives:

$$\int_{\mathbb{R}^n} (u\phi_{y_i}) d\lambda = \int_{\partial U} (u\phi)\nu_i d\mathbf{S} - \int_{\mathbb{R}^n} (\phi u_{y_i}) d\lambda, \ 1 \le i \le n,$$

where ν is the unit outward normal to U. Since φ vanishes on the boundary, we see that the above reduces to:

$$\int_{\mathbb{R}^n} (u\phi_{y_i}) d\lambda = -\int_{\mathbb{R}^n} (\phi u_{y_i}) \ d\lambda, \ 1 \le i \le n.$$

In the general case, for any u ∈ C <sup>m</sup>(R <sup>n</sup>) and any multi-index α = (α1, . . . , αn), |α| = P<sup>n</sup> <sup>α</sup>=1 α<sup>i</sup> = m, we have

$$\int_{\mathbb{R}^n} u(D^{\alpha}\phi)d\lambda = (-1)^m \int_{\mathbb{R}^n} \phi(D^{\alpha}u)d\lambda.$$

Definition 1.5. If α is a multi-index and u, v ∈ L 1 loc(<sup>R</sup> <sup>n</sup>), we say that v is the α th-weak (or distributional) partial derivative of u and write Dαu = v provided that

$$\int_{\mathbb{R}^n} u(D^{\alpha}\phi) d\lambda = (-1)^{|\alpha|} \int_{\mathbb{R}^n} \phi v \ d\lambda$$

for all functions φ ∈ C<sup>∞</sup> c (R <sup>n</sup>). Thus, v is in the dual space D<sup>∗</sup> (R <sup>n</sup>) of D(R <sup>n</sup>).

The next result is easy.

Lemma 1.6. If a weak α th-partial derivative exists for u, then it is a unique λ-a.e. (i.e., except for a set of measure zero).

Definition 1.7. If m ≥ 0 is fixed and 1 ≤ p ≤ ∞, we define the Sobolev space W m,<sup>2</sup> (R <sup>n</sup>) as the set of all locally summable functions u : R <sup>n</sup> → R such that, for each multi-index α with |α| 6 m, D<sup>α</sup>u exists in the weak sense and belongs to L 2 (R <sup>n</sup>). At a deeper level, Sobolev spaces can be defined in at least two different ways [ACM]. On the one hand, once an open domain Ω of R <sup>n</sup> is given and p ∈ [1, ∞) is fixed, we can start from the space C 1 (Ω; ¯ R) which is the subset of C 1 (Ω; R) consisting of functions u such that both u and its gradient admit a continuous ON THE PHYSICAL AND MATHEMATICAL FOUNDATIONS OF QUANTUM PHYSICS VIA FUNCTIONAL INTEGRALS13 extension to the closure of Ω. We then consider the subspace of regular functions u ∈ C 1 (Ω; ¯ R) such that the norm

(18) 
$$||u||_{W^{1,p}(\Omega;\mathbb{R})} = \left[ \left( ||u||_{L^p(\Omega;\mathbb{R})} \right)^p + \left( ||\nabla u||_{L^p(\Omega;\mathbb{R})} \right)^p \right]^{\frac{1}{p}}$$

is finite. The Sobolev space H1,p(Ω; R) is, by definition, the completion of such a subspace of C 1 (Ω; ¯ R) with respect to the W1,p norm. Alternatively, for p ∈ [1, ∞], the Sobolev space W1,p(Ω; R) is the subset of L p (Ω; R) whose elements are weakly differentiable with corresponding derivatives also belonging to L p (Ω; R).

Definition 1.8. If D is a domain in R <sup>n</sup>, we define W m,2 0 (D) as the closure of C<sup>∞</sup> c (D) in W m,<sup>2</sup> (D).

Remark 1.9. Thus, W m,2 0 (D) contains the functions u ∈ W m,<sup>2</sup> (D) such that, for all |α| ≤ m − 1, Dαu = 0 on the boundary of D, ∂D.

We further note that it is also standard to use Hm(D) = W m,<sup>2</sup> (D) and H<sup>m</sup> 0 (D) = W m,2 0 (D).

# 2. The Henstock-Kurzweil Integral

In this section we begin with an elementary introduction to integrals for functions defined on an interval based on Riemann, Darboux, Kurzweil and Henstock. We then discuss the Henstock-Kurzweil (HK) integral for functions defined on R <sup>n</sup> and for operator-valued functions (see [HS1] and [RUD]).

2.1. Riemann's integral. Riemann provided the first rigorous definition of an integral for a real-valued function in 1868.

Definition 2.1. (Riemann) Let f be a bounded real-valued function defined in the interval −∞ < a < b < ∞. For each partition P = {a = x<sup>0</sup> < x<sup>1</sup> · · · < x<sup>n</sup> = b} and each choice of t<sup>j</sup> ∈ [xj−1, x<sup>j</sup> ], define the corresponding Riemann sum by

$$S(f, P, \{t_j\}) = \sum_{j=1}^{n} f(t_j) \Delta x_j,$$

where ∆x<sup>j</sup> = x<sup>j</sup> − xj−1. Define the norm of P by by

$$||P|| = \max \{\Delta x_j, \ 1 \leqslant j \leqslant n\}.$$

We say that f is Riemann integrable over [a, b] if there exists a number I such that, for each ε > 0, there exists a δ > 0 such that, whenever kPk < δ, we have

$$|I - S(f, P, \{t_j\})| < \varepsilon.$$

In this case, we write:

$$I = R \int_{a}^{b} f(x) dx.$$

Remark 2.2. The Riemann integral is used for many applications and numerical approximation purposes. We can also allow f to be complex-valued on [a, b].

2.2. Darboux integral. When f is a real-valued function on [a, b], an integral equivalent to the Riemann integral was defined by Darboux, which is often discussed in elementary calculus.

## Definition 2.3. Let

$$a = x_0 < x_1 < x_2 < \dots < x_{n-1} < x_n = b$$

be a partition P of [a, b] and, set ∆x<sup>j</sup> = x<sup>j</sup> − xj−1. Define

$$M_j = \sup_{x_{j-1} \le t \le x_j} \{ f(t) \}$$
 and  $m_j = \inf_{x_{j-1} \le t \le x_j} \{ f(t) \}$ 

ON THE PHYSICAL AND MATHEMATICAL FOUNDATIONS OF QUANTUM PHYSICS VIA FUNCTIONAL INTEGRALS15 and define the upper and lower Darboux sums by:

$$I^{u} = \sum\nolimits_{j=1}^{n} M_{j} \Delta x_{j} \quad and \quad I_{l} = \sum\nolimits_{j=1}^{n} m_{j} \Delta x_{j}.$$

If the following limits exist

$$\lim_{\|P\| \to 0} I^u = \lim_{\|P\| \to 0} I_l.$$

We say that the Darboux integral exists for the function f(x) on [a, b] and write:

$$I = D \int_a^b f(x) dx = \lim_{\|P\| \to 0} I^u.$$

The most important result of elementary calculus is the Fundamental Theorem of Calculus, which relates differentiation to integration.

Theorem 2.4. (Fundamental Theorem of Calculus, Riemann)Suppose the derivative F 0 (x) = f(x) is Riemann integrable on [a, b], then R <sup>x</sup> a f(t)dt = F(x) − F(a), for all x ∈ [a, b].

The work of Baire, Borel, Cantor and others showed that the concepts of length and area were, like that of number, more delicate than expected, and required a deeper investigation. These issues have led to a closer study of differentiation and integration. In the early 1900s, Lebesgue weakened the conditions for the Fundamental Theorem. This led to the Lebesgue integral and the following:

Theorem 2.5. (Fundamental Theorem of Calculus, Lebesgue)Suppose the derivative F 0 (x) = f(x) exists on [a, b] and F 0 (x) is bounded, then the Lebesgue integral exists and L R x a f(t)dt = F(x) − F(a), for all x ∈ [a, b].

One is naturally led to ask the following question: Is it possible to define an integration process (?) for which the following holds.

Theorem 2.6. If F is differentiable in [a, b], then F 0 is (?)-integrable in [a, b] and the Fundamental Theorem is valid.

Three integration processes exist that accomplish this task, and they were developed by Denjoy, Perron and later by Kurzweil and Henstock (HK) (discussed below). The HK approach is closest to the Riemann integration process, while retaining all the advantages of Lebesgue's theory. Teaching in basic courses still focuses (mainly) on the Riemann and Lebesgue theories. However, the concept Lebesgue measure is very useful in analysis and geometry. In addition, it can be easily extended to more general settings, whereas the corresponding idea associated with the HK integral does not extend as well [GO].

- 2.3. The Henstock-Kurzweil Integral. While studying a generalized approach to ordinary differential equations Kurzweil noticed that, if the base intervals were chosen to be small where the function was steep, and large where it was flat, this was sufficient to extend the Riemann process to obtain a stronger integral than that of Lebesgue (see [KW1]). Independently, Henstock discovered the same approach and provided the first systematic study of the subject (see [HS1]).
- 2.3.1. One-dimensional HK-integral. Let δ(t) > 0 be a function defined on the compact interval [a, b] and let t<sup>1</sup> < t<sup>2</sup> · · · < t<sup>n</sup> be points in the open interval (a, b). The base intervals [xj−1, x<sup>j</sup> ] used to approximate the integral are chosen such that, [xj−1, x<sup>j</sup> ] ⊆ (t<sup>j</sup> − δ(t<sup>j</sup> ), t<sup>j</sup> + δ(t<sup>j</sup> )). We call each (t<sup>j</sup> , [xj−1, x<sup>j</sup> ]) a tagged interval subordinate to δ and the collection:

$$\mathcal{P} = \{ (t_j, [x_{j-1}, x_j]), 1 \le j \le n \}$$

a HK-δ partition of [a, b], if [a, b] = ∪ n <sup>j</sup>=1 [xj−1, x<sup>j</sup> ]. **Definition 2.7.** We say that a function f(t) defined on [a,b] is HK integrable with value I, if for each  $\varepsilon > 0$ , there exists a function  $\delta(t)$  and a HK- $\delta$  partition of [a,b] such that:

$$\left|I - \sum_{j=1}^{n} f(t_j) \Delta x_j\right| < \varepsilon$$

and we write  $I = HK \int_a^b f(x) dx$ .

The following theorem provides additional information.

**Theorem 2.8.** Let  $f(t) : [a, b] \to \mathbb{R}$ .

- (1) If f(t) is Lebesgue integrable in [a,b], then it is HK-integrable in [a,b] and  $HK\int_a^b f(t)dt = L\int_a^b f(t)dt$ .
- (2) If f(t) is HK-integrable in [a,b], then  $\left| \int_a^b f(x) dx \right| < \infty$ .
- (3) If f(t) is HK-integrable and nonnegative in [a,b], then it is Lebesgue integrable in [a,b].
- (4) If f(t) is HK-integrable on every measurable subset of [a,b], then it is Lebesgue integrable in [a,b].

Remark 2.9. Item (4) for the HK-integral is equivalent to the statement that, if every rearrangement of a series converges, it is an absolutely convergent series.

A function has a property (n.e.) if it has the property except for a countable number of points.

Corollary 2.10. Let  $F:[a,b] \to \mathbb{R}$  be continuous. If F is differentiable (n.e.) on [a,b], then F' is HK-integrable and  $HK\int_a^t F'(s)ds = F(t) - F(a)$  for each  $t \in [a,b]$ .

This last result shows the sense in which the HK-integral is the reverse of the derivative. This result is not true for Lebesgue integrals. The standard counterexample [GO] is F 0 (t) = 2tsin(π/t<sup>2</sup> ) − (2π/t)cos(π/t<sup>2</sup> ) for t irrational in (0, 1) and equal to 0 for t rational in (0, 1). (It can be seen that F(t) = t 2 sin(π/t<sup>2</sup> ).)

2.3.2. The n-dimensional HK-integral. There are a number of ways to approach the HK-integral for R <sup>n</sup>. We follow the approach of Lee Tuo-Yeong (see [TY] and [TY1]). All norms are equivalent on R <sup>n</sup>, however, for the HK-integral the maximal norm kxk = max 16k6n |xk|, is natural. With this norm, the open ball B0 (x, r), is a cube centered at x with sides parallel to the coordinate axis of length 2r. (open interval when n = 1.) If the open interval for side i about x<sup>i</sup> is (a<sup>i</sup> , bi), we can represent B<sup>0</sup> (x, r) = (J 0 , x), where J <sup>0</sup> = Q<sup>n</sup> <sup>i</sup>=1(a<sup>i</sup> , bi). (For a closed ball B(x, r) about x, we can represent it as B(x, r) = (J, x), with J = Q<sup>n</sup> <sup>i</sup>=1[a<sup>i</sup> , b<sup>i</sup> ].) Let λn[ · ] be Lebesgue measure on R n.

Definition 2.11. If E is a compact ball in R <sup>n</sup>, a partition P, of E is a collection {(J<sup>i</sup> , xi) : x<sup>i</sup> ∈ J<sup>i</sup> , 1 6 i 6 m}, where J1, J<sup>2</sup> . . . J<sup>m</sup> are non-overlapping intervals (i.e, λ<sup>n</sup> [J<sup>i</sup> ∩ J<sup>j</sup> ] = 0, i 6= j ) and S<sup>m</sup> <sup>i</sup>=1 J<sup>i</sup> = E.

Definition 2.12. If δ is a positive function on E, we say that P is a HK-δ partition for E if for each i, J<sup>i</sup> ⊂ B<sup>0</sup> (x<sup>i</sup> , δ(xi)).

Remark 2.13. The function δ is also called a gauge<sup>1</sup> on E and a HK-δ partition is known as a δ-fine partition of E.

<sup>1</sup>This has nothing to do with the word gauge used for physical theories of fundamental interactions.

**Lemma 2.14.** Cousin's Lemma If  $\delta(\cdot)$  is a positive function in E, then a HK- $\delta$  partition exists for E.

**Definition 2.15.** A function  $f: E \to \mathbb{R}$  is said to be HK-integrable on E, if there exists a number I such that for any  $\varepsilon > 0$  there is a HK- $\delta$  partition on E such that

(19) 
$$\left| \sum_{i=1}^{m} f(\mathbf{x}_i) \lambda_n[J_i] - I \right| < \varepsilon.$$

In this case, we write

$$I = HK \int_{\Gamma} f(\mathbf{x}) d\lambda_n(\mathbf{x}).$$

The following theorem provides a constructive definition of what it means to say that a function is absolutely continuous: (A proof of the following can be found in [GZ].)

**Theorem 2.16.** Let  $f \in L^1[\mathbb{R}^n]$ . If  $\varepsilon > 0$ , then there is a  $\delta > 0$  such that, whenever E is a measurable set with  $\lambda_n[E] < \delta$ ,

$$\left| \int_{E} f(\mathbf{x}) d\lambda_{n}(\mathbf{x}) \right| < \varepsilon.$$

The next result shows that the Lebesgue integral is a special case of the HK-integral, first proven by Davies and Schuss [DZ], but see [GZ].

**Theorem 2.17.** If E is a measurable subset of  $\mathbb{R}^n$  and  $f: E \to \mathbb{R}$ , has a finite Lebesgue integral on E, then

$$HK \int_{E} f(\mathbf{x}) d\lambda_n(\mathbf{x}) = L \int_{E} f(\mathbf{x}) d\lambda_n(x).$$

2.4. Integration of Operator-valued Functions. In this section, we extend the HK-integral to operator-valued functions<sup>2</sup> . Let Ω be a subset of R, B[Ω] be the set of Borel measurable subsets of Ω and let λ be the Lebesgue measure on R. A function A(t), t ∈ Ω, defined on the measure space (Ω, B[Ω], λ) with values in B = L[H], the bounded linear operators on a Hilbert space H, are called operatorvalued functions.

The problem with integration for operator-valued functions is that the functions may not have a Lebesgue integral (see [HP], pg. 71-80). To understand this problem, we begin with:

### Definition 2.18. The function A(t) is said to be:

- (1) almost surely separably valued (or essentially separably valued) if there exists a subset N ⊂ Ω with λ(N) = 0 such that A(Ω \ N) ⊂ B is separable,
- (2) countably-valued if it assumes at most a countable number of values in B, assuming each value 6= 0 on a measurable subset of Ω, and
- (3) strongly measurable if there exists a sequence {An(t))} of countably-valued functions converging (a.s.) to A(t).
- (4) Bochner integrable if kA(t)k<sup>B</sup> is Lebesgue integrable.
- (5) Gelfand-Pettis integrable if hA(t),Li<sup>B</sup> is Lebesgue integrable for each bounded linear functional L ∈ B<sup>∗</sup> .

In order to define constructively the integral of A(t), we must be able to approximate it with simple operator-valued functions in either a strong (Bochner) or weak sense (Gelfand-Pettis). However, the function must be countably-valued and

<sup>2</sup>We refer to a review by Wightman [W] for the related concept of quantum fields as operatorvalued distributions

strongly measurable in the first case or weakly measurable in the second case. In the case of current interest,  $\mathcal{B}$  is not separable, and the family of operator-valued functions  $A(t):\Omega\to\mathcal{B}$  need not be almost separably valued, and hence need not be strongly or weakly measurable. In this section, we extend the HK integral to operator-valued functions on  $\mathbb{R}$ . Because this version will be used for the Feynman operator calculus, we provide some detail.

Let  $[a,b] \subset \mathbb{R}$  and, for each  $t \in [a,b]$ , let  $A(t) \in L(\mathcal{H})$  be a given family of operators.

Recall that, if  $\delta(t)$  maps  $[a,b] \to (0,\infty)$ , and  $\mathcal{P} = \{t_0, \tau_1, t_1, \tau_2, \cdots, \tau_n, t_n\}$ , where  $a = t_0 \leqslant \tau_1 \leqslant t_1 \leqslant \cdots \leqslant \tau_n \leqslant t_n = b$ , we say it is a HK- $\delta$  partition provided that, for  $0 \leqslant i \leqslant n-1$ ,  $t_i, t_{i+1} \in (\tau_{i+1} - \delta(\tau_{i+1}), \tau_{i+1} + \delta(\tau_{i+1}))$ .

**Lemma 2.19.** Let  $\delta_1(t)$  and  $\delta_2(t)$  map  $[a,b] \to (0,\infty)$ , and assume that  $\delta_1(t) \le \delta_2(t)$ . Thus, if  $\mathcal{P}_1$  is a HK- $\delta_1$  partition, it is also a HK- $\delta_2$  partition.

**Definition 2.20.** The family A(t),  $t \in [a,b]$ , is said to have a (uniform) HK-integral if there is an operator Q[a,b] in  $L(\mathcal{H})$  such that, for each  $\varepsilon > 0$ , there exists a  $\delta$  and a HK- $\delta$  partition such that

$$\left\| \sum_{i=1}^{n} \Delta t_i A(\tau_i) - Q[a, b] \right\| < \varepsilon.$$

In this case, we write

$$Q[a,b] = (HK) \int_a^b A(t)dt.$$

**Theorem 2.21.** For  $t \in [a,b]$ , suppose that the operators  $A_1(t)$  and  $A_2(t)$  both have HK-integrals, then their sum has as well and

$$(HK)\int_{a}^{b} [A_{1}(t) + A_{2}(t)]dt = (HK)\int_{a}^{b} A_{1}(t)dt + (HK)\int_{a}^{b} A_{2}(t)dt.$$

**Theorem 2.22.** Suppose  $\{A_k(t) \mid k \in \mathbb{N}\}$  is a family of operator-valued functions in  $L[\mathcal{H}]$ , converging uniformly to A(t) on [a,b], and  $A_k(t)$  has a HK-integral  $Q_k[a,b]$  for each k; then, A(t) has a HK-integral Q[a,b] and  $Q_k[a,b] \to Q[a,b]$  uniformly.

**Theorem 2.23.** Suppose A(t) is Bochner integrable on [a,b], then A(t) has a HK-integral Q[a,b] and:

(20) 
$$(B) \int_{a}^{b} A(t)dt = (HK) \int_{a}^{b} A(t)dt.$$

Proof. First, let E be a measurable subset of [a,b] and assume that  $A(t) = A\chi_E(t)$ , where  $\chi_E(t)$  is the characteristic function of E. In this case, we show that  $Q[a,b] = A\lambda(E)$ , where  $\lambda(E)$  is the Lebesgue measure of E. Let  $\varepsilon > 0$  be given and let D be a compact subset of E. Let  $F \subset [a,b]$  be an open set containing E such that  $\lambda(F\backslash D) < \varepsilon/\|A\|$ ; and define  $\delta: [a,b] \to (0,\infty)$  such that:

$$\delta(t) = \begin{cases} d(t, [a, b] \backslash F), \ t \in E \\ d(t, D), \ t \in [a, b] \backslash E, \end{cases}$$

where d(x,y)=|x-y| denotes the distance function. Let  $\mathcal{P}=\{t_0,\tau_1,t_1,\tau_2,\cdots,\tau_n,t_n\}$  be a HK- $\delta$  partition. If  $\tau_i\in E$  for  $1\leqslant i\leqslant n$ , then  $(t_{i-1},t_i)\subset F$  so that

(21) 
$$\left\| \sum_{i=1}^{n} \Delta t_i A(\tau_i) - A\lambda(F) \right\| = \|A\| \left[ \lambda(F) - \sum_{\tau_i \in E} \Delta t_i \right].$$

On the other hand, if  $\tau_i \notin E$  then  $(t_{i-1}, t_i) \cap D = \emptyset$  (empty set), then it follows that:

(22) 
$$\left\| \sum_{i=1}^{n} \Delta t_i A(\tau_i) - A\lambda(D) \right\| = \|A\| \left[ \sum_{\tau_i, d, E} \Delta t_i - \lambda(D) \right].$$

ON THE PHYSICAL AND MATHEMATICAL FOUNDATIONS OF QUANTUM PHYSICS VIA FUNCTIONAL INTEGRAES

Combining equations (21) and (22), we have that

(23) 
$$\left\| \sum_{i=1}^{n} \Delta t_{i} A(\tau_{i}) - A \lambda(E) \right\| = \|A\| \left[ \sum_{\tau_{i} \in E} \Delta t_{i} - \lambda(E) \right]$$
$$\leq \|A\| \left[ \lambda(F) - \lambda(E) \right] \leq \|A\| \left[ \lambda(F) - \lambda(D) \right] \leq \|A\| \lambda(F \setminus D) < \varepsilon.$$

Suppose that  $A(t) = \sum_{k=1}^{\infty} A_k \chi_{E_k}(t)$ . By definition, A(t) is Bochner integrable if and only if ||A(t)|| is Lebesgue integrable with:

$$(B)\int_{a}^{b} A(t)dt = \sum_{k=1}^{\infty} A_k \lambda(E_k),$$

and (cf. Hille and Phillips [HP])

$$(L) \int_{a}^{b} ||A(t)|| dt = \sum_{k=1}^{\infty} ||A_{k}|| \lambda(E_{k}).$$

As the partial sums converge uniformly, Q[a, b] exists and

$$Q[a,b] \equiv (HK) \int_a^b A(t)dt = (B) \int_a^b A(t)dt.$$

Let A(t) be an arbitrary Bochner integrable operator-valued function in  $L(\mathcal{H})$ , uniformly measurable and defined on [a, b]. By definition, there exists a sequence  $\{A_k(t)\}$  of countably valued operator-valued functions in  $L(\mathcal{H})$  that converges to A(t) in the uniform operator topology such that:

$$\lim_{k \to \infty} (L) \int_{a}^{b} ||A_{k}(t) - A(t)|| dt = 0,$$

and

$$(B) \int_a^b A(t)dt = \lim_{k \to \infty} (B) \int_a^b A_k(t)dt.$$

Since the  $A_k(t)$  are countably-valued,

$$(HK)\int_a^b A_k(t)dt = (B)\int_a^b A_k(t)dt,$$

so

$$(B) \int_a^b A(t)dt = \lim_{k \to \infty} (HK) \int_a^b A_k(t)dt.$$

We are done if we show that Q[a,b] exists. Since every L-integral is a HK integral,  $f_k(t) = ||A_k(t) - A(t)||$  has a HK integral. This means that  $\lim_{k \to \infty} (HK) \int_a^b f_k(t) dt = 0$ . Let  $\varepsilon > 0$  and let m be so large that

$$\left\| (B) \int_{a}^{b} A(t)dt - (HK) \int_{a}^{b} A_{m}(t)dt \right\| < \varepsilon/4$$

and

$$(HK)\int_{a}^{b} f_{k}(t)dt < \varepsilon/4.$$

Choose  $\delta_1$  so that, if  $\{t_0, \tau_1, t_1, \tau_2, \cdots, \tau_n, t_n\}$  is a HK- $\delta_1$  partition, then

$$\left\| (HK) \int_a^b A_m(t) dt - \sum_{i=1}^n \Delta t_i A_m(\tau_i) \right\| < \varepsilon/4.$$

Now choose  $\delta_2$  so that, whenever  $\{t_0, \tau_1, t_1, \tau_2, \cdots, \tau_n, t_n\}$  is a HK- $\delta_2$  partition,

$$\left\| (HK) \int_a^b f_m(t) dt - \sum_{i=1}^n \Delta t_i f_m(\tau_i) \right\| < \varepsilon/4.$$

Set  $\delta = \delta_1 \wedge \delta_2$  so that, by Lemma 2.19,  $\{t_0, \tau_1, t_1, \tau_2, \cdots, \tau_n, t_n\}$  is a HK- $\delta_1$  and HK- $\delta_2$  partition so that:

$$\left\| (B) \int_{a}^{b} A(t)dt - \sum_{i=1}^{n} \Delta t_{i} A(\tau_{i}) \right\| \leq \left\| (B) \int_{a}^{b} A(t)dt - (HK) \int_{a}^{b} A_{m}(t)dt \right\|$$

$$+ \left\| (HK) \int_{a}^{b} A_{m}(t)dt - \sum_{i=1}^{n} \Delta t_{i} A_{m}(\tau_{i}) \right\| + \left| (HK) \int_{a}^{b} f_{m}(t)dt - \sum_{i=1}^{n} \Delta t_{i} f_{m}(\tau_{i}) \right|$$

$$+ (HK) \int_{a}^{b} f_{m}(t)dt < \varepsilon.$$

**Definition 2.24.** A function  $g:[a,b] \subset \mathbb{R} \to \mathcal{H}$  is said to be of bounded variation or BV, if (cf. Appendix and [AFP])

$$\sup_{\mathbb{P}} \left\| \sum_{i=1}^{n} \left[ g(b_i) - g(a_i) \right] \right\| < \infty,$$

for all partitions  $\mathbb{P} = \{(a_1, b_1), \dots, (a_n, b_n)\}$  of the non-overlapping subintervals of [a, b]. In this case, we set

$$\sup_{\mathbb{P}} \left\| \sum_{i=1}^{n} [g(b_i) - g(a_i)] \right\| = BV_a^b(g).$$

**Theorem 2.25.** Let  $g:[a,b] \to \mathcal{H}$ , be of BV.

(1) If h is continuous in [a, b], then

$$I = HK \int_{a}^{b} h(s) dg(s)$$

exists.

(2) If, in addition A is a closed densely defined linear operator on  $\mathcal{H}$ ,  $g \in D(A)$  and Ag(s) = f(s) is of BV, then

(24) 
$$AI = A \int_a^b h(s)dg(s) = \int_a^b h(s)df(s).$$

Let  $(\Omega, \mathfrak{B}(\Omega), \lambda)$  be a measure space, where  $\Omega$  is a subset of  $\mathbb{R}^n$  and  $\lambda$  is the Lebesgue measure. We know that  $L^{\infty}(\Omega, \mathfrak{B}(\Omega), \lambda) = L^1(\Omega, \mathfrak{B}(\Omega), \lambda)^*$ , is the dual space of  $L^1(\Omega, \mathfrak{B}(\Omega), \lambda)$ . For a later study of the Feynman integral, we must describe the dual space  $L^{\infty}(\Omega, \mathfrak{B}(\Omega), \lambda)^*$  of  $L^{\infty}(\Omega, \mathfrak{B}(\Omega), \lambda)$  in more detail. (Note that the indicator function of a set B,  $I_B(x) = 1, x \in B$  and 0 otherwise.)

**Theorem 2.26.** If  $\ell \in L^{\infty}(\Omega, \mathfrak{B}(\Omega), \lambda)^*$ , there is a finitely additive complex-valued signed measure  $\mu_{\ell}$  of bounded total variation and absolutely continuous with respect to  $\lambda$ , such that

$$\ell(\phi) = \int_{\Omega} \phi(x) d\mu_{\ell}(x), \quad \phi \in L^{\infty}[\Omega, \mathfrak{B}[\Omega], \lambda].$$

Thus,  $L^{\infty}(\Omega, \mathfrak{B}(\Omega), \lambda)^* = \mathfrak{M}(\Omega, \mathfrak{B}(\Omega), \lambda)$ , the finitely additive measures on  $\Omega$ .

Proof. According to the Jordan Decomposition Theorem, every complex measure  $\mu$  can be written as  $\mu = \mu_1 + \mu_2 + i(\mu_3 + \mu_4)$ , where  $\mu_1, \mu_3$  are positive measures and  $\mu_2, \mu_4$  are negative measures (see [GZ]). Thus, it is sufficient to prove the theorem when  $\mu$  is real. Let  $\ell \in L^{\infty}(\Omega, \mathfrak{B}(\Omega), \lambda)^*$  and, for each  $B \in \mathfrak{B}[\Omega]$  set  $\mu_{\ell}(B) = \ell(I_B)$ , where  $I_B$  is the indicator function of B. If  $B_1, B_2 \in \mathfrak{B}[\Omega]$ ,  $B_1 \cap B_2 = \emptyset$ , then  $I_{B_1+B_2} = I_{B_1} + I_{B_2}$  so that

$$\ell(I_{B_1} + I_{B_2}) = \ell(I_{B_1}) + \ell(I_{B_2}) \Rightarrow \mu_\ell(B_1 \cup B_2) = \mu_\ell(B_1) + \mu_\ell(B_2).$$

Since

$$\sup_{B \in \mathfrak{B}} |\mu_{\ell}(B)| = \sup_{B \in \mathfrak{B}} |\ell(I_B)| \leqslant ||\ell|| \, ||I_B|| < \infty,$$

we see that  $\mu$  is of bounded variation.

Let  $\phi \in L^{\infty}(\Omega, \mathfrak{B}(\Omega), \lambda)$  be arbitrary. For any  $\varepsilon > 0$ , a simple function  $s_{\varepsilon}$  exists such that

$$s_{\varepsilon} = \sum_{i=1}^{N} a_{i} I_{B_{i}}, \quad \lambda \left( B_{i} \cap B_{j} \right) = 0, \ i \neq j, \quad \bigcup_{i=1}^{N} B_{i} = \Omega$$

and

$$\left\|\phi - \sum_{i=1}^{N} a_i I_{B_i}\right\| < \varepsilon$$
, so that  $\left\|\ell(\phi) - \sum_{i=1}^{N} a_i \mu_{\ell}(B_i)\right\| < \varepsilon \left\|\ell\right\|$ .

ON THE PHYSICAL AND MATHEMATICAL FOUNDATIONS OF QUANTUM PHYSICS VIA FUNCTIONAL INTEGRALS.

It follows that

$$\ell(\phi) = \int_{\Omega} \phi(x) d\mu_{\ell}(x), \quad \text{and} \quad \|\ell\| = \sup_{ess \sup |\phi| \leqslant 1} \left| \int_{\Omega} \phi(x) d\mu_{\ell}(x) \right|.$$

Finally, because  $\lambda(B) = 0$  implies that  $I_B = 0$  (a.s.), it follows that  $\mu_{\ell}(B) = 0$  so that  $\mu_{\ell}$  is absolutely continuous with respect to  $\lambda$  (by definition, see [GZ]).

From here, we see that  $L^1(\Omega, \mathfrak{B}(\Omega), \lambda)^{**} = \mathfrak{M}(\Omega, \mathfrak{B}(\Omega), \lambda)$  and, the injection of  $L^1(\Omega, \mathfrak{B}(\Omega), \lambda) \to \mathfrak{M}(\Omega, \mathfrak{B}(\Omega), \lambda)$  is dense. Because  $L^1(\mathbb{R}^n, \mathfrak{B}(\mathbb{R}^n), \lambda)$  is a Banach algebra under convolution, it is easy to prove that property.

Corollary 2.27.  $\mathfrak{M}(\mathbb{R}^n, \mathfrak{B}(\mathbb{R}^n), \lambda)$  is a Banach algebra under convolution.

Recall that, if  $\Omega$  is an open subset of  $\mathbb{R}^n$ , then  $\mathbb{C}_c(\Omega)$  is the set of all continuous functions defined on  $\Omega$ , that vanish outside a compact set.

Corollary 2.28. If  $\phi \in \mathbb{C}_c(\Omega)$ , then for each  $\ell \in \mathbb{C}_c(\Omega)^*$ , there is a countably additive complex measure  $\mu_\ell$  such that.

$$\ell(\phi) = \int_{\Omega} \phi(x) d\mu_{\ell}(x).$$

### 3. Lebesgue Measure on $\mathbb{R}^{\infty}$

In this section we briefly review the theory of Lebesgue's measure on infinitedimensional space (see [GZ]). We used a basic requirement for any mapping that would serve as an acceptable version of the Lebesgue measure on  $\mathbb{R}^{\infty}$ . If  $I_0 =$  $\left[-\frac{1}{2},\frac{1}{2}\right]^{\aleph_0}$ , we know that any definition  $\lambda_{\infty}(\cdot)$  of Lebesgue measure must satisfy  $\lambda_{\infty}[I_0] = 1$ . This requirement is the centerfold of our approach.

Let  $I = [-\frac{1}{2}, \frac{1}{2}]$ ,  $I_n = \times_{i=n+1}^{\infty} I$ , so that  $I_0 = \times_{i=1}^{\infty} I$  and let  $\mathbb{R}_I^n = \mathbb{R}^n \times I_n$ . It is evident that  $\mathbb{R}_I^n \subset \mathbb{R}_I^{n+1}$ . Since this is an increasing sequence, we can define  $\mathbb{R}_I^{n}$ 

28

by:

$$\mathbb{R'}_I^{\infty} = \lim_{n \to \infty} \mathbb{R}_I^n = \bigcup_{k=1}^{\infty} \mathbb{R}_I^k.$$

Let  $\tau_1$  be the topology on  $\mathbb{R}'_I^{\infty} = \mathfrak{X}_1$  induced by the class of open sets  $\mathfrak{O}$  defined by:

$$\mathfrak{O} = \bigcup_{n=1}^{\infty} \mathfrak{O}_n = \bigcup_{n=1}^{\infty} \{ U_n : U_n = U \times I_n, \ U \text{ open in } \mathbb{R}^n \},$$

and let  $\tau_2$  be the discrete topology on  $\mathbb{R}^{\infty} \setminus \mathbb{R}'_I^{\infty} = \mathfrak{X}_2$  induced by the metric  $d_2$ , for which  $d_2(x,y) = 1$ ,  $x \neq y$  and  $d_2(x,y) = 0$ , x = y, for all  $x,y \in \mathfrak{X}_2$ .

**Definition 3.1.** We define  $(\mathbb{R}_I^{\infty}, \tau) = (\mathfrak{X}_1, \tau_1) \oplus (\mathfrak{X}_2, \tau_2)$ , the topological direct sum of  $(\mathfrak{X}_1, \tau_1)$  and  $(\mathfrak{X}_2, \tau_2)$ , so that every open set in  $(\mathbb{R}_I^{\infty}, \tau)$  is the disjoint union of two sets  $G_1 \cup G_2$ , where  $G_1$  is open in  $(\mathfrak{X}_1, \tau_1)$  and  $G_2$  is open in  $(\mathfrak{X}_2, \tau_2)$ .

**Definition 3.2.** If  $A_n = A \times I_n$ ,  $B_n = B \times I_n \in \mathbb{R}^n_I$ , we define:

- $(1) A_n \cup B_n = A \cup B \times I_n,$
- (2)  $A_n \cap B_n = A \cap B \times I_n$ , and
- (3)  $B_n^c = B^c \times I_n$ .

**Definition 3.3.** We define  $\mathbb{R}_I^n = \mathbb{R}^n \times I_n \subset R^\infty$ . If T is a linear transformation on  $\mathbb{R}^n$  and  $A_n = A \times I_n$ , we define the extension operator  $T_e$  on  $\mathbb{R}_I^n$  by  $T_e[A_n] = T[A] \times I_n$ .

We define the topology on  $\mathbb{R}^n_I$  via the following class of open sets:

$$\mathfrak{O}_n = \{U \times I_n :, U \text{ open in } \mathbb{R}^n \}$$

and let  $\mathfrak{B}[\mathbb{R}^n_I]$  be the natural Borel  $\sigma$ -algebra.

It follows from our construction that  $\mathbb{R}_I^{\infty} = \mathbb{R}^{\infty}$  as sets, but not as topological spaces. It is easy to prove the following result, which shows that convergence in the  $\tau$ -topology always implies convergence in the Tychonoff topology.

**Theorem 3.4.** If  $y_k$  converges to x in the  $\tau$ -topology, then  $y_k$  converges to x in the Tychonoff topology.

In a similar manner, if  $\mathfrak{B}(\mathbb{R}^n_I)$  is the Borel  $\sigma$ -algebra for  $\mathbb{R}^n_I$  (i.e., the smallest  $\sigma$ -algebra generated by the  $\mathfrak{O}_n$ ), then  $\mathfrak{B}(\mathbb{R}^n_I) \subset \mathfrak{B}(\mathbb{R}^{n+1}_I)$ , hence we can define  $\mathfrak{B}'(\mathbb{R}^\infty_I)$  by:

$$\mathfrak{B}'(\mathbb{R}_I^\infty) = \mathrm{lim}_{n \to \infty} \mathfrak{B}(\mathbb{R}_I^n) = \mathop{\cup}\limits_{k=1}^\infty \mathfrak{B}(\mathbb{R}_I^k).$$

If  $\mathcal{P}(A)$  denotes the power set of A, let  $\mathfrak{B}(\mathbb{R}_I^{\infty})$  be the smallest  $\sigma$ -algebra containing  $\mathfrak{B}'(\mathbb{R}_I^{\infty}) \cup \mathcal{P}(R_I^{\infty} \setminus \bigcup_{n=1}^{\infty} \mathbb{R}_I^n)$ . It is evident that the class  $\mathfrak{B}(\mathbb{R}_I^{\infty})$  coincides with the Borel  $\sigma$ -algebra generated by the  $\tau$ -topology on  $\mathbb{R}^{\infty}$ . For any  $A \in \mathfrak{B}[\mathbb{R}^n]$ , we define  $\lambda_{\infty}(A_n)$  on  $\mathbb{R}_I^n$  by product measure:

$$\lambda_{\infty}(A_n) = \lambda_n(A) \times \prod_{i=n+1}^{\infty} \lambda_1(I) = \lambda_n(A),$$

**Theorem 3.5.**  $\lambda_{\infty}$  is a measure on  $\mathfrak{B}(\mathbb{R}^n_I)$ , equivalent to n-dimensional Lebesgue measure on  $\mathbb{R}^n$ .

*Proof.* If  $A = \underset{i=1}{\overset{\infty}{\times}} A_i \in \mathfrak{B}(\mathbb{R}^n_I)$ , then  $\lambda(A_i) = 1$  for i > n such that the series  $\lambda_{\infty}(A) = \prod_{i=1}^{\infty} \lambda(A_i)$  always converges. Furthermore,

(25) 
$$0 < \lambda_{\infty}(A) = \prod_{i=1}^{\infty} \lambda(A_i) = \prod_{i=1}^{n} \lambda(A_i) = \lambda_n (\underset{i=1}{\overset{n}{\times}} A_i).$$

Bearing in mind that sets of the type  $A = \underset{i=1}{\overset{n}{\times}} A_i$  generate  $\mathfrak{B}(\mathbb{R}^n)$ , we see that  $\lambda_{\infty}(\cdot)$ , restricted to  $\mathbb{R}^n_I$ , is equivalent to  $\lambda_n(\cdot)$ .

Corollary 3.6. The measure  $\lambda_{\infty}(\cdot)$  is translationally and rotationally invariant on  $(\mathbb{R}^n_I, \mathfrak{B}[\mathbb{R}^n_I])$ .

- 3.1. Extension of  $\lambda_{\infty}(\cdot)$  to  $\mathbb{R}_{I}^{\infty}$ . From the previous section, we see that the extension of  $\lambda_{\infty}(\cdot)$  to  $\mathfrak{B}(\mathbb{R}_{I}^{\infty})$  provides the following important consequences:
  - (1)  $\mathbb{R}_I^n \in \mathfrak{B}(\mathbb{R}_I^{\infty})$  for all n and  $\lambda_{\infty}(\cdot)$  is equivalent to  $\lambda_n(\cdot)$  on  $\mathbb{R}_I^n$ ,
  - (2)  $\mathfrak{B}(\mathbb{R}^\infty_I)$  has a large number of open sets of finite measure and
  - (3)  $\lambda_{\infty}(I_0) = 1$ .

It is shown in [GZ] that  $\lambda_{\infty}(\cdot)$  can be extended to a countably additive measure on  $\mathfrak{B}(\mathbb{R}_{I}^{\infty})$  in a constructive manner that is closely related to the same approach used for Lebesgue measure on  $\mathbb{R}^{n}$ . We now consider an equivalent definition of  $\lambda_{\infty}(\cdot)$  that proves useful, in that it is direct and allows us to relate our measure to a number of other approaches. We first recall the following characterization of a measure.

**Theorem 3.7.** If A is a  $\sigma$ -algebra, the mapping  $\mu : A \to [0, \infty]$  is a measure if and only if

- (1)  $\mu(\emptyset) = 0$ .
- (2) If  $A, B \in \mathcal{A}$ , and  $\mu(A \cap B) = 0$ , then  $\mu(A \cup B) = \mu(A) + \mu(B)$ .
- (3) If  $\{B_n\} \subset \mathcal{A} \text{ and } B_n \subset B_{n+1}, \text{ then }$

$$\mu\left(\bigcup_{k=1}^{\infty} B_k\right) = \lim_{n \to \infty} \mu\left(B_n\right).$$

*Proof.* If  $\mu$  is a measure, these conditions are necessary. Thus, we need only prove that these conditions are sufficient. Since  $\mu$  is positive and finitely additive, it

ON THE PHYSICAL AND MATHEMATICAL FOUNDATIONS OF QUANTUM PHYSICS VIA FUNCTIONAL INTEGRALS31

suffices to show that it is countably additive. Let 
$$\{A_n\} \subset \mathcal{A}$$
 be disjoint. From

$$\bigcup_{k=1}^{n} A_k \subset \bigcup_{k=1}^{n+1} A_k, \quad \text{and} \quad \mu\left(\bigcup_{k=1}^{n} A_k\right) = \sum_{k=1}^{n} \mu(A_k),$$

if we let

$$B_n = \bigcup_{k=1}^n A_k,$$

we have B<sup>n</sup> ⊂ Bn+1, so we can apply condition (3) to get

(26) 
$$\lim_{n \to \infty} \mu(B_n) = \mu\left(\bigcup_{k=1}^{\infty} B_k\right) = \mu\left(\bigcup_{k=1}^{\infty} A_k\right)$$

and

(27) 
$$\lim_{n \to \infty} \mu(B_n) = \lim_{n \to \infty} \mu\left(\bigcup_{k=1}^n A_k\right) = \lim_{n \to \infty} \sum_{k=1}^n \mu(A_k) = \sum_{k=1}^\infty \mu(A_k).$$

Combining conditions (2) and (3) proves sufficiency.

Definition 3.8. We define a measure m<sup>k</sup> on B(R<sup>∞</sup>) by

$$m_k(A) = \lambda_{\infty}(A \cap \mathbb{R}^k_I), \text{ for each } A \in \mathfrak{B}(\mathbb{R}^{\infty}),$$

and set

(28) 
$$m(A) = \lim_{k \to \infty} m_k(A).$$

Theorem 3.9. The mapping m : B(R<sup>∞</sup>) → [0, ∞] is a measure.

Proof. It is clear that conditions (1) and (2) of Theorem 3.7 are satisfied by m, thus we need only check (3). Since A<sup>i</sup> ∩ A<sup>j</sup> = ∅ unless i = j, we see that the same is true for A<sup>i</sup> ∩ R k I and A<sup>j</sup> ∩ R k I . Let N ∈ N and note that

$$\lambda_{\infty}\left[\left(\bigcup_{i=1}^{N}A_{i}\right)\cap\mathbb{R}_{I}^{k}\right]=\sum_{i=1}^{N}\lambda_{\infty}\left(A_{i}\cap\mathbb{R}_{I}^{k}\right).$$

Since A<sup>i</sup> ∩ R k <sup>I</sup> ⊂ A<sup>i</sup> ∩ R k+1 I , all the terms are increasing, and therefore

$$m\left(\bigcup_{i=1}^{N} A_{i}\right) = \lim_{k \to \infty} m_{k}\left(\bigcup_{i=1}^{N} A_{i}\right) = \sum_{i=1}^{N} \lim_{k \to \infty} m_{k}\left(A_{i}\right) = \sum_{i=1}^{N} m\left(A_{i}\right).$$

Letting N → ∞, completes the proof.

Corollary 3.10. The completion mˆ (·), of m(·) is equivalent to λ∞(·) defined on B(R<sup>∞</sup> I ).

Theorem 3.11. λ<sup>∞</sup> [R<sup>∞</sup> \ (∪<sup>∞</sup> <sup>n</sup>=1R n I )] = 0.

Proof. Using the extension of m<sup>k</sup> to B(R k I ), we have

$$\lambda_{\infty} \left[ \left( \mathbb{R}^{\infty} \setminus \bigcup_{k=1}^{\infty} \mathbb{R}_{I}^{k} \right) \right] = \lim_{n \to \infty} \hat{m}_{n} \left[ \left( \mathbb{R}^{\infty} \setminus \bigcup_{k=1}^{\infty} \mathbb{R}_{I}^{k} \right) \cap \mathbb{R}_{I}^{n} \right] = 0.$$

Corollary 3.12. λ∞(X2) = 0.

Theorem 3.13. There exists a family of sets {Ak} ⊂ B(R<sup>∞</sup> I ) with λ∞(Ak) < ∞, and a set N with λ∞(N) = 0 such that:

(29) 
$$\mathbb{R}_I^{\infty} = \bigcup_{k=1}^{\infty} A_k \cup N.$$

Proof. Because λ∞(·) is concentrated on X1, we set N = X2. To show that (28) holds, let {xk} be the set of vectors in R<sup>∞</sup> with rational coordinates and let B<sup>k</sup> be the unit cube with center xk, so that λ∞(Bk) = 1 for all k and

$$\mathbb{R}_I^{\infty} = \cup_{k=1}^{\infty} B_k.$$

Let 
$$A_k = B_k \setminus \mathfrak{X}_2$$
.

#### 4. Lebesgue Measure on Hilbert Space

In this section we construct the Lebesgue measure on a separable Hilbert space. Let  $\mathcal{H}$  be a separable Hilbert space and let  $\{e_n\}$  be an orthogonal basis for  $\mathcal{H}$  such that for any  $x \in \mathcal{H}$  there is a sequence of scalars  $(x_n)_{n=1}^{\infty}$ , such that  $x = \sum_{n=1}^{\infty} x_n e_n$  and  $\sum_{n=1}^{\infty} |x_n|^2 < \infty$ 

**Definition 4.1.** We define  $\mathcal{H}_I^n$  and  $\mathcal{H}_I$  by:

$$\mathcal{H}_{I}^{n} = \left\{ (x_k)_{k=1}^{n} : x = \sum_{k=1}^{n} x_k e_k \in \mathcal{H} \right\} \times I_n$$

and

$$\mathcal{H}_I = \left\{ (x_k)_{k=1}^{\infty} : x = \sum_{n=1}^{\infty} x_k e_k \in \mathcal{H} \right\}.$$

It can be seen that  $\mathcal{H}_I^n = \mathbb{R}_I^n$ . This observation will prove important when we discuss integration in the next section.

**Definition 4.2.** Define  $T_n: \mathcal{H} \to \mathcal{H}_I^n$  by  $T_n(x) = (x_k)_{k=1}^n$  and  $T: \mathcal{H} \to \mathcal{H}_I$ ,  $T(x) = (x_k)_{k=1}^{\infty}$ . Define

$$\mathfrak{B}(\mathcal{H}_I) = \mathcal{H}_I \cap \mathfrak{B}(\mathbb{R}_I^{\infty})$$

and define

$$\mathfrak{B}_{I}[\mathcal{H}] = \left\{ T^{-1}(A) \mid A \in \mathfrak{B}[\mathcal{H}_{I}] \right\} =: T^{-1} \left\{ \mathfrak{B}\left[\mathcal{H}_{I}\right] \right\}.$$

**Definition 4.3.** Let  $J_k = \left[ -\frac{1}{2k^{3/2}}, \frac{1}{2k^{3/2}} \right]$  and  $g_k(x_k) = k^{3/2} \chi_{J_k}(x_k)$ . For each n, we define  $\nu_n$  by:

$$\nu_n = \underset{k=1}{\overset{n}{\otimes}} \lambda \otimes \left( \underset{k=n+1}{\overset{\infty}{\otimes}} \mu_k \right) = \lambda_n \otimes \left( \underset{k=n+1}{\overset{\infty}{\otimes}} \mu_k \right), \text{ where } d\mu_k = g_k(x_k) d\lambda(x_k)$$

and for each  $A \in \mathfrak{B}_I(\mathcal{H})$ , we define  $\lambda_{\mathcal{H}}(A)$  by:

$$\lambda_{\mathcal{H}}(A) = \lim \sup_{n \to \infty} \nu_n[T(A) \cap \mathbb{R}_I^n].$$

**Theorem 4.4.** The measure  $\lambda_{\mathcal{H}}$  is a nonvanishing  $\sigma$ -finite Borel measure on  $\mathcal{H}$ .

*Proof.* The result follows if we show that  $\lambda_{\mathcal{H}}$  is nonzero and countably additive. Since

$$\left\| \sum_{k=1}^{\infty} \frac{e_k}{k^{3/2}} \right\|_{\mathcal{H}} \leqslant \sum_{k=1}^{\infty} \frac{1}{k^{3/2}} < \infty,$$

we see that the set  $A = \left\{ \sum_{k=1}^{\infty} x_k e_k : |x_k| \leqslant \frac{1}{k^{3/2}} \right\} \in \mathfrak{B}_I(\mathcal{H})$ , so that

$$\lambda_{\mathcal{H}}(A) \geqslant \nu_1(T[A]) = \lambda(I_1) \prod_{k=2}^{\infty} \mu_k(J_k) = 1.$$

We now note that  $\nu_1[T(A)] \geqslant \nu_n(T[A])$  for all n and, because  $\{\nu_n\}$  is an increasing sequence, we see that  $\lambda_{\mathcal{H}}(A) = 1$ .

To see that  $\lambda_{\mathcal{H}}$  is countably additive, let  $\{A_k\}$  be disjoint and let  $B_k = \bigcup_{i=1}^k A_i$ . Because the family  $\{B_k\}$  increases for each n,

$$\nu_n \left[ T \left( \bigcup_{i=1}^{\infty} A_i \right) \cap \mathbb{R}_I^n \right] = \lim_{k \to \infty} \nu_n \left[ T(B_k) \cap \mathbb{R}_I^n \right] = \sum_{i=1}^{\infty} \nu_n \left[ T \left( A_i \right) \cap \mathbb{R}_I^n \right],$$

so that

$$\lambda_{\mathcal{H}}\left[\bigcup_{i=1}^{\infty} A_{i}\right] = \sum_{i=1}^{\infty} \lim \sup_{n \to \infty} \nu_{n} \left[T\left(A_{i}\right) \cap \mathbb{R}_{I}^{n}\right] = \sum_{i=1}^{\infty} \lambda_{\mathcal{H}} \left[A_{i}\right].$$

**Definition 4.5.** We call  $\mathcal{H}_I$  the canonical representation of  $\mathcal{H}$  in  $\mathbb{R}_I^{\infty}$ .

The following theorem (proved in [GZ]), characterizes the properties of the Lebesgue measure on  $\mathbb{R}_I^{\infty}$  (or  $\mathcal{H}$ ). We state the results for  $\mathbb{R}_I^{\infty}$ .

**Theorem 4.6.** The measure space  $(\mathbb{R}_I^{\infty}, \mathfrak{B}[\mathbb{R}_I^{\infty}], \lambda_{\infty})$  has the following properties:

ON THE PHYSICAL AND MATHEMATICAL FOUNDATIONS OF QUANTUM PHYSICS VIA FUNCTIONAL INTEGRALS35

- (1) λ∞(X2) = 0.
- (2) For every A ∈ L[R<sup>∞</sup> I ] (Lebesgue sets) and ε > 0, there exists a compact set F ⊂ A and an open set G ⊃ A such that λ∞(G \ F) < ε, so that λ∞(·) is regular.
- (3) There exists a family of compact sets {An} ⊂ B[R<sup>∞</sup> I ], with λ∞[An] < ∞ and a set N with λ∞[N] = 0, such that R<sup>∞</sup> <sup>I</sup> = S<sup>∞</sup> <sup>n</sup>=1 A<sup>n</sup> ∪ N (i.e., λ∞(·) is σ-finite).
- (4) For A ∈ B[R<sup>∞</sup> I ], λ∞(A − x) = λ∞(A) if and only if x ∈ `1.
- 4.1. Measurable Functions. In this subsection we discuss measurable functions on R<sup>∞</sup> I (or H). Let x = (x1, x2, x3, . . .) ∈ R<sup>∞</sup> I . Fixing n with I<sup>n</sup> = Q<sup>∞</sup> <sup>k</sup>=n+1 -− 1 2 , 1 2 , we set hn(ˆx) = χI<sup>n</sup> (ˆx) (indicator function of In), where ˆx = (xi)<sup>∞</sup> <sup>i</sup>=n+1.

Definition 4.7. Let M<sup>n</sup> represent the class of Lebesgue measurable functions on R <sup>n</sup>. If x ∈ R<sup>∞</sup> I and f <sup>n</sup> ∈ Mn, let x¯ = (xi) n <sup>i</sup>=1 and define an essentially tame measurable function of order n (or en-tame) on R<sup>∞</sup> I by f(x) = f <sup>n</sup>(¯x) ⊗ hn(ˆx). We let

$$\mathcal{M}_I^n = \{ f(x) : f(x) = f^n(\bar{x}) \otimes h_n(\hat{x}), \ x \in \mathbb{R}_I^\infty \}$$

be the class of all en-tame functions.

Definition 4.8. A function f : R<sup>∞</sup> <sup>I</sup> → R is said to be measurable and we write f ∈ M<sup>I</sup> , if there is a sequence {f<sup>n</sup> ∈ M<sup>n</sup> I } of en-tame functions, such that limn→∞ fn(x) → f(x) λ∞-(a.e).

The existence of functions satisfying Definition 4.8 is not obvious, so we have [GM]:

Theorem 4.9. (Existence) Suppose that f : R<sup>∞</sup> <sup>I</sup> → (−∞, ∞) and f −1 (A) ∈ B[R<sup>∞</sup> I for all A ∈ B[R]. Then there exists a family of functions {fn}, f<sup>n</sup> ∈ M<sup>n</sup> I , such that fn(x) → f(x), λ∞-(a.e).

Remark 4.10. From Theorem 4.6 (1), we can see that any set A, of nonzero measure is concentrated in X<sup>1</sup> (i.e., λ∞(A) = λ∞(A ∩ X1)). It also follows that the essential support of the limit function f(x) in definition 2.14 (i.e., {x |f(x) 6= 0}), is concentrated in R N I , for some N.

4.2. Integration Theory on R<sup>∞</sup> <sup>I</sup> or H. In this section all theorems remain true when R<sup>∞</sup> I is replaced by H, hence it is sufficient to provide all results for R<sup>∞</sup> I .

We provide a constructive theory of integration on R<sup>∞</sup> <sup>I</sup> using the known properties of integration on R n I . This approach has the advantage that all standard theorems for the Lebesgue measure apply. (The proofs are the same as for integration on R <sup>n</sup>.) Let L 1 [R n I ] be the class of integrable functions on R n I . Since L 1 [R n I ] ⊂ L 1 [R n+1 I ], we define L 1 [Rˆ<sup>∞</sup> I ] = S<sup>∞</sup> <sup>n</sup>=1 L 1 [R n I ].

Definition 4.11. We say that a measurable function f ∈ L 1 [R<sup>∞</sup> I ], if there is a Cauchy sequence {fn} ⊂ L 1 [Rˆ<sup>∞</sup> I ], with f<sup>n</sup> ∈ L 1 [R n I ] and limn→∞fn(x) = f(x), λ∞- (a.e).

Theorem 4.12. L 1 [Rˆ<sup>∞</sup> I ] = L 1 [R<sup>∞</sup> I ].

Proof. We know that L 1 [Rˆ<sup>∞</sup> I ] ⊃ L 1 [R n I ], for all n, hence it suffices to prove that L 1 [Rˆ<sup>∞</sup> I ] is closed. Let f be a limit point of L 1 [Rˆ<sup>∞</sup> I ] (f ∈ L 1 [R<sup>∞</sup> I ]). If f = 0, we are done, thus assume f 6= 0. From the remarks above, we know that if A<sup>f</sup> is the support of f, then λ∞(A<sup>f</sup> ) = λ∞(A<sup>f</sup> ∩X1). Thus, A<sup>f</sup> ∩X<sup>1</sup> ⊂ R N I for some N. This ON THE PHYSICAL AND MATHEMATICAL FOUNDATIONS OF QUANTUM PHYSICS VIA FUNCTIONAL INTEGRALS37 implies that there is a function f <sup>0</sup> ∈ L 1 [R N+1 I ], with λ<sup>∞</sup> ({x : f(x) 6= f 0 (x)}) = 0. It follows that, f(x) = f 0 (x) (a.e). Recalling that L 1 [R n I ] is a set of equivalence classes, we observe that L 1 [Rˆ<sup>∞</sup> I ] = L 1 [R<sup>∞</sup> I ].

Definition 4.13. If f ∈ L 1 [R<sup>∞</sup> I ], we define the integral of f by

$$\int_{\mathbb{R}_{I}^{\infty}} f(x)d\lambda_{\infty}(x) = \lim_{n \to \infty} \int_{\mathbb{R}_{I}^{\infty}} f_{n}(x)d\lambda_{\infty}(x),$$

where {fn} ⊂ L 1 [R<sup>∞</sup> I ] is any Cauchy sequence converging to f(x)-(a.e).

Theorem 4.14. If f ∈ L 1 [R<sup>∞</sup> I ], then the above integral exists and all theorems that are true for f ∈ L 1 [R n I ], also hold for f ∈ L 1 [R<sup>∞</sup> I ].

## 5. The Kuelbs-Steadman Spaces

If we want a class of spaces designed to include the HK-integrable functions on R n I , for a fixed n, let Q<sup>n</sup> <sup>I</sup> be the set {x = (x1, x<sup>2</sup> · · ·) ∈ R n I } such that the first n coordinates (x1, x<sup>2</sup> · · · xn) are rational. Since every countable set is isomorphic to the set of natural numbers, Q<sup>n</sup> I is a countable set which is also dense in R n I . Therefore, we can arrange it as Q<sup>n</sup> <sup>I</sup> = {x1, x2, x3, · · ·}. For each l and i, let Bl(xi) be the closed cube centered at x<sup>i</sup> , with sides parallel to the coordinate axes and edge e<sup>l</sup> = 1 2 <sup>l</sup>−1√ n , l ∈ N. Now choose the natural order that maps N × N bijectively to N, and let {Bk, k ∈ N} be the resulting set of (all) closed cubes {Bl(xi) |(l, i) ∈ N × N} centered at a point in Q<sup>n</sup> I . Let Ek(x) be the characteristic function of Bk, so that Ek(x) is in L p [R n I ] ∩ L<sup>∞</sup>[R n I ] for 1 ≤ p < ∞. We define Fk( · ) on L 1 [R n I ] as

(30) 
$$F_k(f) = \int_{\mathbb{R}_I^n} \mathcal{E}_k(\mathbf{x}) f(\mathbf{x}) d\lambda_{\infty}(\mathbf{x}) = \int_{\mathbf{B}_k} f(\mathbf{x}) d\lambda_{\infty}(\mathbf{x}).$$

Since each  $\mathbf{B}_k$  is a cube with sides parallel to the coordinate axes,  $F_k(\cdot)$  is defined for all HK-integrable functions, and is a bounded linear functional on  $L^p[\mathbb{R}_I^n]$  for each k, with  $||F_k||_{\infty} \leq 1$ . Moreover, if  $F_k(f) = 0$  for all k, then f = 0, thus  $\{F_k\}$  is fundamental on  $L^p[\mathbb{R}_I^n]$  for  $1 \leq p < \infty$ . Fix  $t_k > 0$  such that  $\sum_{k=1}^{\infty} t_k = 1$  and define a new inner product  $(\cdot)$  on  $L^2[\mathbb{R}_I^n]$  by

(31) 
$$(f,g) = \sum_{k=1}^{\infty} t_k \left[ \int_{\mathbb{R}_I^n} \mathcal{E}_k(\mathbf{x}) f(\mathbf{x}) d\lambda_{\infty}(\mathbf{x}) \right] \overline{\left[ \int_{\mathbb{R}_I^n} \mathcal{E}_k(\mathbf{y}) g(\mathbf{y}) d\lambda_{\infty}(\mathbf{y}) \right]}.$$

The completion of  $L^2[\mathbb{R}_I^n]$  in this inner product is the Kuelbs-Steadman space,  $KS^2[\mathbb{R}_I^n]$ . To see directly that  $KS^2[\mathbb{R}_I^n]$  contains the HK integrable functions, let f be HK integrable,

$$\|f\|_{KS^2}^2 = \sum_{k=1}^{\infty} t_k \left| \int_{\mathbb{R}_I^n} \mathcal{E}_k(\mathbf{x}) f(\mathbf{x}) d\lambda_{\infty}(\mathbf{x}) \right|^2 \leqslant \sup_k \left| \int_{\mathbf{B}_k} f(\mathbf{x}) d\lambda_{\infty}(\mathbf{x}) \right|^2 < \infty,$$

hence  $f \in KS^2[\mathbb{R}^n_I]$ .

**Theorem 5.1.** For each  $p, 1 \leq p \leq \infty$ ,  $KS^2[\mathbb{R}_I^n] \supset L^p[\mathbb{R}_I^n]$  as a continuous dense subspace.

*Proof.* By construction,  $KS^2[\mathbb{R}_I^n]$  contains  $L^2[\mathbb{R}_I^n]$  densely thus, we only need to show that  $KS^2[\mathbb{R}_I^n] \supset L^q[\mathbb{R}_I^n]$  for  $q \neq 2$ . If  $f \in L^q[\mathbb{R}_I^n]$  and  $q < \infty$ , then we have

$$\begin{aligned} \|f\|_{KS^2} &= \left[ \sum_{k=1}^{\infty} t_k \left| \int_{\mathbb{R}_I^n} \mathcal{E}_k(\mathbf{x}) f(\mathbf{x}) d\lambda_{\infty}(\mathbf{x}) \right|^{\frac{2q}{q}} \right]^{1/2} \\ &\leqslant \left[ \sum_{k=1}^{\infty} t_k \left( \int_{\mathbb{R}_I^n} \mathcal{E}_k(\mathbf{x}) \left| f(\mathbf{x}) \right|^q d\lambda_{\infty}(\mathbf{x}) \right)^{\frac{2}{q}} \right]^{1/2} \\ &\leqslant \sup_k \left( \int_{\mathbb{R}_I^n} \mathcal{E}_k(\mathbf{x}) \left| f(\mathbf{x}) \right|^q d\lambda_{\infty}(\mathbf{x}) \right)^{\frac{1}{q}} \leqslant \|f\|_q \,. \end{aligned}$$

Thus,  $f \in KS^2[\mathbb{R}_I^n]$ . For  $q = \infty$ , first note that  $vol(\mathbf{B}_k)^2 \leq \left[\frac{1}{\sqrt{n}}\right]^{2n} \leq 1$ , hence we have

$$||f||_{KS^2} = \left[ \sum_{k=1}^{\infty} t_k \left| \int_{\mathbf{R}_I^n} \mathcal{E}_k(\mathbf{x}) f(\mathbf{x}) d\lambda_{\infty}(\mathbf{x}) \right|^2 \right]^{1/2}$$

$$\leq \left[ \left[ \sum_{k=1}^{\infty} t_k [vol(\mathbf{B}_k)]^2 \right] [ess \sup |f|]^2 \right]^{1/2} \leq ||f||_{\infty}.$$

Therefore,  $f \in KS^2[\mathbb{R}_I^n]$  and  $L^{\infty}[\mathbb{R}_I^n] \subset KS^2[\mathbb{R}_I^n]$ .

Before proceeding to an additional discussion, we construct the  $KS^p[\mathbb{R}^n_I]$  spaces, for  $1 \leq p \leq \infty$ . Define:

$$||f||_{KS^p} = \begin{cases} \left\{ \sum_{k=1}^{\infty} t_k \left| \int_{\mathbb{R}_I^n} \mathcal{E}_k(\mathbf{x}) f(\mathbf{x}) d\lambda_{\infty}(\mathbf{x}) \right|^p \right\}^{1/p}, 1 \leqslant p < \infty, \\ \sup_{k \geqslant 1} \left| \int_{\mathbb{R}_I^n} \mathcal{E}_k(\mathbf{x}) f(\mathbf{x}) d\lambda_{\infty}(\mathbf{x}) \right|, p = \infty. \end{cases}$$

It is easy to see that  $\|\cdot\|_{KS^p}$  defines a norm on  $L^p$ . If  $KS^p$  is the completion of  $L^p$  with respect to this norm, we have:

**Theorem 5.2.** For each  $q, 1 \leq q \leq \infty, KS^p[\mathbb{R}_I^n] \supset L^q[\mathbb{R}_I^n]$  as a dense continuous embedding.

*Proof.* As in the previous theorem, by construction  $KS^p[\mathbb{R}^n_I]$  contains  $L^p[\mathbb{R}^n_I]$  densely, hence we only need to show that  $KS^p[\mathbb{R}^n_I] \supset L^q[\mathbb{R}^n_I]$  for  $q \neq p$ . First, suppose that  $p < \infty$ . If  $f \in L^q[\mathbb{R}^n_I]$  and  $q < \infty$ , then we have

$$||f||_{KS^{p}} = \left[ \sum_{k=1}^{\infty} t_{k} \left| \int_{\mathbb{R}_{I}^{n}} \mathcal{E}_{k}(\mathbf{x}) f(\mathbf{x}) d\lambda_{\infty}(\mathbf{x}) \right|^{\frac{qp}{q}} \right]^{1/p}$$

$$\leq \left[ \sum_{k=1}^{\infty} t_{k} \left( \int_{\mathbb{R}_{I}^{n}} \mathcal{E}_{k}(\mathbf{x}) |f(\mathbf{x})|^{q} d\lambda_{\infty}(\mathbf{x}) \right)^{\frac{p}{q}} \right]^{1/p}$$

$$\leq \sup_{k} \left( \int_{\mathbb{R}_{I}^{n}} \mathcal{E}_{k}(\mathbf{x}) |f(\mathbf{x})|^{q} d\lambda_{\infty}(\mathbf{x}) \right)^{\frac{1}{q}} \leq ||f||_{q}.$$

Hence,  $f \in KS^p[\mathbb{R}^n]$ . For  $q = \infty$ , we have

$$||f||_{KS^p} = \left[ \sum_{k=1}^{\infty} t_k \left| \int_{\mathbb{R}_I^n} \mathcal{E}_k(\mathbf{x}) f(\mathbf{x}) d\lambda_{\infty}(\mathbf{x}) \right|^p \right]^{1/p}$$

$$\leq \left[ \left[ \sum_{k=1}^{\infty} t_k [vol(\mathbf{B}_k)]^p \right] [ess \sup |f|]^p \right]^{1/p} \leq M ||f||_{\infty}.$$

Thus,  $f \in KS^p[\mathbb{R}^n_I]$  and  $L^\infty[\mathbb{R}^n_I] \subset KS^p[\mathbb{R}^n_I]$ . The case  $p = \infty$  is obvious.  $\square$ 

**Theorem 5.3.** For  $KS^p$ ,  $1 \le p \le \infty$ , we have:

- (1) If  $f, g \in KS^p$ , then  $||f + g||_{KS^p} \le ||f||_{KS^p} + ||g||_{KS^p}$  (Minkowski inequality).
- (2) If K is a weakly compact subset of  $L^p$ , then it is a compact subset of  $KS^p$ .
- (3) If  $1 , then <math>KS^p$  is uniformly convex.
- (4) If  $1 and <math>p^{-1} + q^{-1} = 1$ , the dual space of  $KS^p$  is  $KS^q$ .
- (5)  $KS^{\infty} \subset KS^p$ , for  $1 \le p < \infty$ .

**Theorem 5.4.** For each p,  $1 \leq p \leq \infty$ , the test functions  $\mathcal{D} \subset KS^p(\mathbb{R}^n_I)$  as a continuous embedding.

Proof. Because  $KS^{\infty}(\mathbb{R}^n_I)$  is continuously embedded in  $KS^p(\mathbb{R}^n_I)$ ,  $1 \leq q < \infty$ , it suffices to prove the result for  $KS^{\infty}(\mathbb{R}^n_I)$ . Suppose that  $\phi_j \to \phi$  in  $\mathcal{D}[\mathbb{R}^n_I]$ , so that there exists a compact set  $K \subset \mathbb{R}^n_I$ , containing the support of  $\phi_j - \phi$  and  $D^{\alpha}\phi_j$  converges to  $D^{\alpha}\phi$  uniformly on K for every multi-index  $\alpha$ . Let  $L = \{l \in \mathbb{N} : \text{the support of } \mathcal{E}_l, stp\{\mathcal{E}_l\} \subset K\}$ , then

$$\lim_{j \to \infty} \left\| D^{\alpha} \phi - D^{\alpha} \phi_{j} \right\|_{KS} = \lim_{j \to \infty} \sup_{l \in L} \left| \int_{\mathbb{R}^{n}_{I}} \left[ D^{\alpha} \phi\left(x\right) - D^{\alpha} \phi_{j}\left(x\right) \right] \mathcal{E}_{l}\left(x\right) d\lambda_{\infty}\left(x\right) \right|$$

$$\leq vol(\mathbf{B}_{l}) \lim_{j \to \infty} \sup_{x \in K} |D^{\alpha} \phi(x) - D^{\alpha} \phi_{j}(x)| \leqslant \lim_{j \to \infty} \sup_{x \in K} |D^{\alpha} \phi(x) - D^{\alpha} \phi_{j}(x)| = 0.$$

It follows that  $\mathcal{D}[\mathbb{R}^n_I] \subset KS^p[\mathbb{R}^n_I]$  as a continuous embedding, for  $1 \leq p \leq \infty$ . Thus, using the Hahn-Banach theorem, we see that the Schwartz distributions,  $\mathcal{D}^*[\mathbb{R}^n_I] \subset [KS^p(\mathbb{R}^n_I)]^*$ , for  $1 \leq p \leq \infty$ . ON THE PHYSICAL AND MATHEMATICAL FOUNDATIONS OF QUANTUM PHYSICS VIA FUNCTIONAL INTEGRALS41

5.1. The Family KS<sup>p</sup> [R<sup>∞</sup> I ]. We can now construct the spaces KS<sup>p</sup> [R<sup>∞</sup> I ], 1 ≤ p ≤ ∞, using the same approach that led to L 1 [R<sup>∞</sup> I ]. Since KS<sup>p</sup> [R n I ] ⊂ KS<sup>p</sup> [R n+1 I ], we define KS<sup>p</sup> [Rˆ<sup>∞</sup> I ] = S<sup>∞</sup> <sup>n</sup>=1 KS<sup>p</sup> [R n I ].

Definition 5.5. We say that a measurable function f ∈ KS<sup>p</sup> [R<sup>∞</sup> I ], if there is a Cauchy sequence {fn} ⊂ KS<sup>p</sup> [Rˆ<sup>∞</sup> I ], with f<sup>n</sup> ∈ KS<sup>p</sup> [R n I ] and limn→∞fn(x) = f(x), λ∞-(a.e).

The same proof as Theorem 1.18 shows that functions in KS<sup>p</sup> [Rˆ<sup>∞</sup> I ] differ from functions in its closure KS<sup>p</sup> [R<sup>∞</sup> I ], by sets of measure zero.

Theorem 5.6. KS<sup>p</sup> [Rˆ<sup>∞</sup> I ] = KS<sup>p</sup> [R<sup>∞</sup> I ].

Definition 5.7. If f ∈ KS<sup>p</sup> [R<sup>∞</sup> I ], we define the integral of f by:

$$\int_{\mathbb{R}_{I}^{\infty}} f(x)d\lambda_{\infty}(x) = \lim_{n \to \infty} \int_{\mathbb{R}_{I}^{\infty}} f_{n}(x)d\lambda_{\infty}(x),$$

where f<sup>n</sup> ∈ KS<sup>p</sup> [R n I ] is any Cauchy sequence converging to f(x).

Theorem 5.8. If f ∈ KS<sup>p</sup> [R<sup>∞</sup> I ], then the above integral exists and all theorems that are true for f ∈ KS<sup>p</sup> [R n I ], also hold for f ∈ KS<sup>p</sup> [R<sup>∞</sup> I ].

Remark 5.9. Theorem 5.8 is true if (R<sup>∞</sup> I , λ∞) is replaced by (H, λH).

## 6. Feynman Operator Calculus

In this section we describe the construction of Feynman's time ordered operator calculus. We first prove that a unique definition of time exists based on an assumption that is weaker than the cosmological principle (i.e., the universe is homogenous and isotropic in the large). If O and O<sup>0</sup> are any two observers of a particle with velocity  $\mathbf{v}$  ( $\mathbf{v}'$ ) and  $\gamma^{-1}(\mathbf{v}) = \sqrt{1 - \frac{\mathbf{v}^2}{c^2}}$ ,  $\mathbf{v} = \frac{d\mathbf{x}}{dt}$ , the standard definition of proper time is:

(32) 
$$d\tau = \gamma^{-1}(\mathbf{v}) dt, \quad d\tau = \gamma^{-1}(\mathbf{v}') dt'.$$

Using the relations  $H = \gamma(\mathbf{v}) mc^2$  for observer O and  $H' = \gamma(\mathbf{v}') mc^2$  for observer O', we obtain

(33) 
$$d\tau = \frac{mc^2}{H}dt, \quad d\tau = \frac{mc^2}{H'}dt'.$$

**Theorem 6.1.** If the universe is representable, in the sense that the ratio of the observed total mass energy to the total energy is independent of any observer's particular portion of the universe, then the universe has a unique clock.

*Proof.* Let O and O' be two observers viewing the universe. If  $Mc^2$  is the total mass energy and H is the total energy Hamiltonian for the universe as seen by O, under the stated conditions  $H/Mc^2$  is constant. By assumption this view is observer independent, so that O' will also obtain the same ratio. It follows from (33) that

$$d\tau = \frac{Mc^2}{H}dt = \frac{Mc^2}{H}dt' = d\tau_h.$$

Thus,  $\tau_h = \frac{Mc^2}{H}t$  uniquely defines a universal clock for any observer.

The fifth parameter was first introduced by Fock, and later by Stueckelberg followed by Feynman and Schwinger. They each considered it a global clock in addition to the Minkowski spacetime. Horwitz and Piron [HP], and Fanchi and Collins named it historical time. They were the first to realize that there should be physical justification for this clock (see [FA]).

- 6.1. **Feynman-Dyson Space.** Let  $\mathcal{H} = KS^2[\mathbb{R}_I^n]$ , with a fixed orthonormal basis  $\{e^i\}$ . Let J = [-T, T] be an interval of historical time and, for each  $t \in J$ , define  $\mathcal{H}(t) = \mathcal{H}$  and let  $\mathcal{H}_{\otimes} = \hat{\otimes}_t \mathcal{H}(t)$  be the continuous tensor product Hilbert space of von Neumann over J. For each  $i \in \mathbb{N}$ , set  $e^i_t = e^i$  and  $E^i = \otimes_t e^i_t$ . We define  $\mathcal{FD}^i \subset \mathcal{H}_{\otimes}$  to be the smallest Hilbert space containing  $E^i$ . We call  $\mathcal{FD} = \bigoplus_{i=1}^{\infty} \mathcal{FD}^i$  the Feynman-Dyson space over J for  $\mathcal{H}$  (the film for Feynman's spacetime events).
- 6.2. **Time Ordered Operators.** If  $\mathcal{C}(\mathcal{H}_{\otimes})$  is a set of closed densely defined linear operators on  $\mathcal{H}_{\otimes}$ , and  $\{H_I(t), t \in J\}$  is a family of generators for unitary groups, we define  $\mathcal{C}(\mathcal{H}(t)) \subset \mathcal{C}(\mathcal{H}_{\otimes})$  by:

$$\mathcal{C}(\mathcal{H}(t)) = \left\{ \mathbf{H}_I(t) \mid \mathbf{H}_I(t) = \widehat{\bigotimes}_{T \geqslant s > t} \mathbf{I}_s \otimes H_I(t) \otimes (\bigotimes_{t > s \geqslant -T} I_s) \right\},\,$$

where  $I_s$  is the identity operator at time s. It follows that

$$\mathbf{H}_{I}(t)\mathbf{H}_{I}(s) = \mathbf{H}_{I}(s)\mathbf{H}_{I}(t), t \neq s.$$

Thus, the operators are ordered in time, commute when acting at different times and maintain their mathematically defined positions.

6.3. **Time Ordered Integrals.** We now state the fundamental theorem for time-ordered operators (see [GZ]):

**Theorem 6.2.** (Fundamental Theorem for Time-Ordered Operators) If  $\{H_I(t)|t \in J\}$  is a family of weakly continuous Hamiltonian generators of unitary groups on  $\mathcal{H}$  and  $\{\mathbf{H}_I(t)|t \in J\}$  is the time ordered version defined on  $\mathcal{FD}^2_{\otimes}$  then:

(1) The family  $\{\mathbf{H}_I(t)|t\in J\}$  is strongly continuous and the time-ordered HK-integral  $\mathbf{Q}[t, -T] = \int_{-T}^t \mathbf{H}_I(s) ds$  exists (a.e), has a dense domain and is

the generator of a strongly continuous unitary group U[t, −T] on FD<sup>2</sup> ⊗. Moreover:

(2) If Ψ<sup>0</sup> is in the domain of Q [t, −T] then Ψ (t) = U [t, −T] Ψ<sup>0</sup> = exp − i <sup>~</sup>Q [t, −T] Ψ<sup>0</sup> satisfies:

(34) 
$$i\hbar \frac{\partial \Psi(t)}{\partial t} = \mathbf{H}_{I}(t) \Psi(t), \ \Psi(-T) = \Psi_{0}.$$

6.4. Interaction Representation. Haag's theorem shows that the (sharp) equaltime commutation relations of an interacting field are equivalent to those of a free field (Streater and Wightman [SW]). Lindner et. al. have shown that there is experimental support for interference in time of the wave function for a particle (see Horwitz [HO2]). In this section we show that if time smearing exists, the interaction representation is defined. Before proving this result, we consider the concept of an exchange operator.

Definition 6.3. An exchange operator E[t, t<sup>0</sup> ], is defined for pairs t, t<sup>0</sup> such that:

- (1) E[t, t<sup>0</sup> ] : C[H(t)] → C[H(t 0 )], (bijective mapping),
- (2) E[s, t<sup>0</sup> ]E[t, s] = E[t, t<sup>0</sup> ],
- (3) E[t, t<sup>0</sup> ]E[t 0 , t] = I⊗, (identity)
- (4) for s 6= t, t<sup>0</sup> , E[t, t<sup>0</sup> ]H<sup>I</sup> (s) = H<sup>I</sup> (s), for all H<sup>I</sup> (s) ∈ C[H(s)].

Assume that H0(t) and H1(t) are generators of unitary groups for each t ∈ J, H<sup>I</sup> (t) = H0(t) ⊕ H1(t) is densely defined, H<sup>n</sup> 1 (t) = nH1(t)R (n, H1(t)) is the Yosida approximator for H1(t) and Theorem (6.2) is satisfied. Define Un[t, a], ON THE PHYSICAL AND MATHEMATICAL FOUNDATIONS OF QUANTUM PHYSICS VIA FUNCTIONAL INTEGRALS

 $\mathbf{U}_0[t,a]$  and  $\mathbf{U}_0^{\kappa}[t,a]$  as

$$\begin{aligned} \mathbf{U}_n[t,a] &= \exp\{(-i/\hbar) \int_a^t [\boldsymbol{H}_0(s) \oplus \boldsymbol{H}_1^n(s)] ds\}, \\ \mathbf{U}_0[t,a] &= \exp\{(-i/\hbar) \int_a^t \boldsymbol{H}_0(s) ds\}, \\ \boldsymbol{U}_0^{\kappa}[t,a] &= \exp\{(-i/\hbar) \int_a^t \boldsymbol{H}_0^{\kappa}(s) ds\}, \\ \boldsymbol{H}_0^{\kappa}(t) &= \int_{-\infty}^{\infty} \rho_{\kappa}(t,s) \mathbf{E}[t,s] \boldsymbol{H}_0(s) ds, \end{aligned}$$

where  $\rho_{\kappa}(t,s)$  is the smearing density, which depends on the Planck length  $\kappa$  and  $\int_{-\infty}^{\infty} \rho_{\kappa}(t,s) ds = 1 \text{ (for example, } \rho_{\kappa}(t,s) = \left[1 \middle/ \sqrt{2i\pi\kappa^2}\right] \exp\{i(t-s)^2 \middle/ 2\kappa^2\}).$ 

We now obtain:

$$\boldsymbol{H}_{I}^{n}(t) = \mathbf{U}_{0}^{\kappa}[a, t]\boldsymbol{H}_{1}^{n}(t)\mathbf{U}_{0}^{\kappa}[t, a],$$

and the terms do not commute. If we set  $\Psi_n(t) = \mathbf{U}_0^{\kappa}[a,t]\mathbf{U}_n[t,a]\Phi$ , we have

$$\begin{split} &\frac{\partial}{\partial t}\Psi_n(t) = \frac{i}{\hbar}\mathbf{U}_0^{\kappa}[a,t]\boldsymbol{H}_0(t)\mathbf{U}_n[t,a]\Phi - \frac{i}{\hbar}\mathbf{U}_0^{\kappa}[a,t]\left[\boldsymbol{H}_0(t) + \boldsymbol{H}_1^n(t)\right]\mathbf{U}_n[t,a]\Phi \\ &\text{so that} \quad \frac{\partial}{\partial t}\Psi_n(t) = \frac{i}{\hbar}\{\mathbf{U}_0^{\kappa}[a,t]\boldsymbol{H}_1^n(t)\mathbf{U}_0^{\kappa}[t,a]\}\mathbf{U}_0^{\kappa}[a,t]\mathbf{U}_n[t,a]\Phi \\ &\text{and} \quad i\hbar\frac{\partial}{\partial t}\Psi_n(t) = \boldsymbol{H}_I^n(t)\Psi_n(t), \ \Psi_n(a) = \Phi. \end{split}$$

With the same conditions as Theorem (6.2), we have

**Theorem 6.4.** If  $Q_1[t,a] = \int_a^t H_1(s)ds$  generates a unitary group on  $\mathcal{H}$ , then the time-ordered integral  $\mathbf{Q}_I[t,a] = \int_a^t \mathbf{H}_I(s)ds$ , where  $\mathbf{H}_I(t) = \mathbf{U}_0^{\kappa}[a,t]\mathbf{H}_1(t)\mathbf{U}_0^{\kappa}[t,a]$  generates a unitary group on  $\mathcal{FD}_{\otimes}^2$ , and

$$\exp\{(-i/\hbar)\mathbf{Q}_I^n[t,a]\} \to \exp\{(-i/\hbar)\mathbf{Q}_I[t,a]\},$$

where  $\mathbf{Q}_{I}^{n}[t,a] = \int_{a}^{t} \mathbf{H}_{I}^{n}(s)ds$ , and:

$$i\hbar \frac{\partial}{\partial t} \Psi(t) = \boldsymbol{H}_I(t)\Psi(t), \ \Psi(a) = \Phi.$$

- 6.5. **Feynman Path Integral I.** Recall that  $KS^2[\mathbb{R}_I^n] = KS^2[\mathbb{R}_I^n]^* \supset L^1[\mathbb{R}_I^n]^{**} = \mathfrak{M}[\mathbb{R}_I^n]$ , the space of finitely (and countable) additive measures on  $\mathbb{R}_I^n$ . (It also contains  $\mathcal{D}^*[\mathbb{R}_I^n]$ , the space of distributions.) Thus, the Schrödinger equation with Dirac measure initial data is well posed on  $KS^2[\mathbb{R}_I^n]$ . In this section we construct the Feynman path integral in a manner that preserves both intuitive and computational advantages. We begin with a few operator extensions to  $KS^2[\mathbb{R}_I^n]$ .
- 6.5.1. Extensions to  $KS^2[\mathbb{R}_I^n]$ . Bearing in mind that the position operator  $\mathbf{x}$  and momentum operator  $\mathbf{p}$  are closed and densely defined on  $L^2[\mathbb{R}^n]$  and, our extension operator  $T_e$  extends this property to  $L^2[\mathbb{R}_I^n]$ , it is easy to see that both have closed, densely defined extensions to  $KS^2[\mathbb{R}_I^n]$ . If  $f, g \in L^1[\mathbb{R}_I^n]$ , we denote the Fourier transform of f and the convolution of g with respect to f by  $\mathfrak{F}(f)$  and  $\mathfrak{C}_f(g)$ , respectively. The following theorem extends them to bounded linear operators on  $KS^2[\mathbb{R}_I^n]$ , ensuring that both the Schrödinger and Heisenberg theories have faithful representations on  $KS^2[\mathbb{R}_I^n]$ .

**Theorem 6.5.** For each  $f, g \in L^1[\mathbb{R}^n_I] \vee L^2[\mathbb{R}^n_I]$ , the Fourier transform  $\mathfrak{F}(f)$  and the convolution of g with respect to f,  $\mathfrak{C}_f(g)$ , both extend to  $KS^2[\mathbb{R}^n_I]$  as bounded linear operators.

If  $B \in \mathcal{B}(\mathbb{R}^n_I)$  is bounded,  $\mathbf{x} \in \mathbb{R}^n_I$  fixed and  $t, s \in \mathbb{R}$ , then the kernel:

$$\mathbb{K}_{\mathbf{f}}[t, \mathbf{x}; s, B] = \int_{B} (2\pi i (t - s))^{-1/2} \exp\{i |\mathbf{x} - \mathbf{y}|^{2} / 2(t - s)\} d\mathbf{y}$$

ON THE PHYSICAL AND MATHEMATICAL FOUNDATIONS OF QUANTUM PHYSICS VIA FUNCTIONAL INTEGRALS47 is in KS<sup>2</sup> [R n I ], it is a finitely additive measure with kK<sup>f</sup> [t, <sup>x</sup> ; s, B]kKS<sup>2</sup> <sup>6</sup> 1 and

$$\mathbb{K}_{\mathbf{f}}[t,\mathbf{x}\,;\,s,B] = \int_{\mathbb{R}^n_I} \mathbb{K}_{\mathbf{f}}[t,\mathbf{x}\,;\,\tau,d\mathbf{z}] \mathbb{K}_{\mathbf{f}}[\tau,\mathbf{z}\,;\,s,B], \ \ (\text{HK-integral}).$$

Definition 6.6. Let P<sup>n</sup> = {t0, τ1, t1, τ2, · · · , τn, tn} one of a family of HK-δ<sup>n</sup> partitions of the interval [0, t] for each n, with lim supn→∞ δn(τ ) = 0. Set ∆t<sup>j</sup> = t<sup>j</sup> − tj−1, τ<sup>0</sup> = 0 and for ψ ∈ KS<sup>2</sup> [R n I ] define

(35)
$$\int_{\mathbb{R}_{I}^{n[0,t]}} \mathbb{K}_{\mathbf{f}}[\mathcal{D}_{\lambda}\mathbf{x}(\tau) ; \mathbf{x}(0)] = e^{-\lambda t} \sum_{k=0}^{[|\lambda t|]} \frac{(\lambda t)^{k}}{k!} \left\{ \prod_{j=1}^{k} \int_{\mathbb{R}_{I}^{n}} \mathbb{K}_{\mathbf{f}}[t_{j}, \mathbf{x}(\tau_{j}) ; t_{j-1}, d\mathbf{x}(\tau_{j-1})] \right\},$$

and

(36) 
$$\int_{\mathbb{R}_{I}^{n[0,t]}} \mathbb{K}_{\mathbf{f}}[\mathcal{D}\mathbf{x}(\tau);\mathbf{x}(0)]\psi[\mathbf{x}(0)]$$
$$= \lim_{\lambda \to \infty} \int_{\mathbb{R}_{I}^{n[0,t]}} \mathbb{K}_{\mathbf{f}}[\mathcal{D}_{\lambda}\mathbf{x}(\tau);\mathbf{x}(0)]\psi[\mathbf{x}(0)]$$

whenever the limit exists.

The following is clear. A general result is presented in the next section.

Theorem 6.7. The function ψ(x) ≡ 1 ∈ KS<sup>2</sup> [R n I ] and

(37)
$$\int_{\mathbb{R}_{I}^{n[s,t]}} \mathbb{K}_{\mathbf{f}}[\mathcal{D}\mathbf{x}(\tau) ; \mathbf{x}(s)] = \mathbb{K}_{\mathbf{f}}[t, \mathbf{x} ; s, \mathbf{y}] = \frac{1}{\sqrt{2\pi i(t-s)}} \exp\{i|\mathbf{x} - \mathbf{y}|^{2} / 2(t-s)\}.$$

The above result is exactly what Feynman expected.

6.6. Path Integral II. Feynman suggested that his path integral is a special case of the time-ordered operator calculus. In this section we consider the advantages of time ordering.

**Theorem 6.8.** If U[t,a] is any well defined evolution operator on  $KS^2(\mathbb{R}^n_I)$ , with a time-dependent generator H(t), and reproducing kernel  $K[\mathbf{x}(t),t\,;\,\mathbf{x}(s),s]$  such that:

$$\begin{split} K\left[\mathbf{x}(t),\,t;\,\mathbf{x}(s),\,s\right] &= \int_{\mathbb{R}^n_I} K\left[\mathbf{x}(t),\,t;\,d\mathbf{x}(\tau),\,\tau\right] K\left[\mathbf{x}(\tau),\,\tau;\,\mathbf{x}(s),\,s\right],\\ U[t,a]\varphi(a) &= \int_{\mathbb{R}^n_I} K\left[\mathbf{x}(t),\,t;\,d\mathbf{x}(s),\,s\right]\varphi(s),\,\,such\,\,that\\ &\frac{\partial}{\partial t} U\left[t,a\right]\varphi\left(a\right) = \dot{u}\left(t\right) = H\left(t\right)u\left(t\right),\,\,u\left(a\right) = \varphi\left(a\right). \end{split}$$

Then the corresponding time-ordered version  $\mathbf{U}[t,s]$  defined on  $\mathcal{FD}^2_{\otimes} \subset \mathcal{H}^2_{\otimes}$ , with  $\text{kernel } \mathbb{K}_{\mathbf{f}}[\mathbf{x}(t),\,t;\,\mathbf{x}(s),\,s]$  satisfies the conditions of our fundamental theorem.

Since 
$$\mathbf{U}[t,\tau]\mathbf{U}[\tau,s] = \mathbf{U}[t,s],$$

we have:

$$\mathbb{K}_{\mathbf{f}}\left[\mathbf{x}(t),\,t;\,\mathbf{x}(s),\,s\right] = \int_{\mathbb{R}_{I}^{n}} \mathbb{K}_{\mathbf{f}}\left[\mathbf{x}(t),\,t;\,d\mathbf{x}(\tau),\,\tau\right] \mathbb{K}_{\mathbf{f}}\left[\mathbf{x}(\tau),\,\tau;\,\mathbf{x}(s),\,s\right].$$

From our sum over paths representation for  $\mathbf{U}[t,s]$ , we have:

$$\mathbf{U}[t, s]\Phi(s) = \lim_{\lambda \to \infty} \mathbf{U}_{\lambda}[t, s]\Phi(s)$$
$$= \lim_{\lambda \to \infty} e^{-\lambda(t-s)} \sum_{k=0}^{n} \frac{\left[\lambda(t-s)\right]^{k}}{k!} \mathbf{U}_{k}[t, s]\Phi(s),$$

where  $n = [|\lambda(t-s)|]$  (the greatest integer in  $\lambda(t-s)$ ) and

$$\mathbf{U}_{k}[t,s]\Phi(s) = \exp\left\{ (-i/\hbar) \sum_{j=1}^{k} \int_{t_{j-1}}^{t_{j}} \mathbf{E}[\tau_{j},\tau] \mathbf{H}(\tau) d\tau \right\} \Phi(s).$$

As before, define  $\mathbb{K}_{\mathbf{f}}[\mathcal{D}_{\lambda}\mathbf{x}(t) ; \mathbf{x}(s)]$  by

$$\mathbb{K}_{\mathbf{f}}[\mathcal{D}_{\lambda}\mathbf{x}(t) ; \mathbf{x}(s)]$$

$$=:e^{-\lambda(t-s)}\sum_{k=1}^{n}\frac{\left[\lambda(t-s)\right]}{k!}\left\{\prod_{j=1}^{k}\int_{\mathbb{R}_{I}^{n}}\mathbb{K}_{f}\left[t_{j},\mathbf{x}\left(t_{j}\right);d\mathbf{x}\left(t_{j-1}\right),t_{j-1}\right]|^{\tau_{j}}\right\}$$

and  $|\tau_j|$  denotes that the integration is performed at the time  $\tau_i$ .

Definition 6.9. We define the Feynman path integral associated with U[t, s] by:

$$\mathbf{U}[t,s]\Phi(s) = \int_{\mathbb{R}_I^{n[t,s]}} \mathbb{K}_{\mathbf{f}}[\mathcal{D}\mathbf{x}(\tau) \ ; \ \mathbf{x}(s)]\Phi(s) = \lim_{\lambda \to \infty} \int_{\mathbb{R}_I^{n[t,s]}} \mathbb{K}_{\mathbf{f}}[\mathcal{D}_{\lambda}\mathbf{x}(\tau) \ ; \ \mathbf{x}(s)]\Phi(s).$$

Theorem 6.10. For the Feynman time-ordered theory, whenever a reproducing kernel exists on KS<sup>2</sup> [R n I ], we have

$$\lim_{\lambda \to \infty} \mathbf{U}_{\lambda}[t, s] \Phi(s) = \mathbf{U}[t, s] \Phi(s) = \int_{\mathbb{R}^{n[t, s]}_{I}} \mathbb{K}_{\mathbf{f}}[\mathcal{D}^{\lambda} \mathbf{x}(\tau) ; \mathbf{x}(s)] \Phi[\mathbf{x}(s)],$$

and the limit is independent of the space of continuous functions that we choose.

Remark 6.11. This result includes the Wiener path integral, which requires additional effort to restrict it to the space of continuous paths. We also note that the base space H is allowed to continuously change in time and the family of spaces need not be all Hilbert. In addition, the intersection of the corresponding domains of the generators can even be the empty set (see Goldstein [GS]).

Let us assume that H0(t) and H1(t) are strongly continuous generators of unitary groups, with a common dense domain D(t), for each t ∈ J = [a, b], where H1,ρ(t) = ρH1(t)R(ρ, H1(t)) is the Yosida approximator for the time-ordered version of H1(t), with dense domain D = ⊗t∈ID(t). Define U<sup>ρ</sup> [t, a] and U<sup>0</sup> [t, a] as follows:

$$\mathbf{U}^{\rho}[t,a] = \exp\{(-i/\hbar) \int_{a}^{t} [\boldsymbol{H}_{0}(s) + \boldsymbol{H}_{1,\rho}(s)]ds\},$$
  
$$\mathbf{U}^{0}[t,a] = \exp\{(-i/\hbar) \int_{a}^{t} \boldsymbol{H}_{0}(s)ds\}.$$

Because H1,ρ(s) is bounded, H0(s) +H1,ρ(s) is a generator of a unitary group for each s ∈ J and a finite ρ. Now assume that U<sup>0</sup> [t, a] has an associated reproducing kernel such that  $\mathbf{U}^0[t,a] = \int_{\mathbb{R}_I^{n[t,s]}} \mathbb{K}_{\mathbf{f}}[\mathcal{D}\mathbf{x}(\tau);\mathbf{x}(a)]$ . We now have the following general result.

**Theorem 6.12.** (Extended Feynman-Kac) If  $\mathbf{H}_0(s) \oplus \mathbf{H}_1(s)$  is a generator of a unitary group, then

$$\lim_{\rho \to \infty} \mathbf{U}^{\rho}[t, a] \Phi(a) = \mathbf{U}[t, a] \Phi(a)$$

$$= \int_{\mathbb{R}_{I}^{n[t, a]}} \mathbb{K}_{\mathbf{f}}[\mathcal{D}\mathbf{x}(\tau) ; \mathbf{x}(a)] \exp\{(-i/\hbar) \int_{a}^{\tau} \boldsymbol{H}_{1}(s) ds]\} \Phi[\mathbf{x}(a)].$$

*Proof.* The fact that  $\mathbf{U}^{\rho}[t,a]\Phi(a) \to \mathbf{U}[t,a]\Phi(a)$  is clear. To prove that

$$\mathbf{U}[t,a]\Phi(a) = \int_{\mathbb{R}_I^{n[t,a]}} \mathbb{K}_{\mathbf{f}}[\mathcal{D}\mathbf{x}(\tau);\mathbf{x}(a)] \exp\{(-i/\hbar) \int_a^t \boldsymbol{H}_1(s) ds\} \Phi(a).$$

First, note that because the time-ordered integral exists and we are only interested in the limit, we can write for each k,

$$U_k^{\rho}[t,a]\Phi(a) = \exp\left\{\left(-i/\hbar\right)\sum_{j=1}^k \int_{t_{j-1}}^{t_j} \left[\mathbf{E}[\tau_j,s]\boldsymbol{H}_0(s) + \mathbf{E}[\tau_j',s]\boldsymbol{H}_{1,\rho}(s)\right] ds\right\}\Phi(a),$$

where  $\tau_j$  and  $\tau'_j$  are distinct points in  $(t_{j-1}, t_j)$ . Thus, we can also write  $U_k^{\rho}[t, a]\Phi(a)$  as

$$\begin{aligned} &\mathbf{U}_{k}^{\rho}[t,a] = exp\left\{\frac{-i}{\hbar}\sum_{j=1}^{k}\int_{t_{j-1}}^{t_{j}}\mathbf{E}[\tau_{j},s]\boldsymbol{H}_{0}(s)ds\right\}exp\left\{\frac{-i}{\hbar}\sum_{j=1}^{k}\int_{t_{j-1}}^{t_{j}}\mathbf{E}[\tau'_{j},s]\boldsymbol{H}_{1,\rho}(s)ds\right\} \\ &= \prod_{j=1}^{k}exp\left\{\frac{-i}{\hbar}\int_{t_{j-1}}^{t_{j}}\mathbf{E}[\tau_{j},s]\boldsymbol{H}_{0}(s)ds\right\}exp\left\{\frac{-i}{\hbar}\sum_{j=1}^{k}\int_{t_{j-1}}^{t_{j}}\mathbf{E}[\tau'_{j},s]\boldsymbol{H}_{1,\rho}(s)ds\right\} \\ &= \prod_{j=1}^{k}\int_{\mathbb{R}^{n}_{I}}\mathbb{K}_{\mathbf{f}}\left[t_{j},\mathbf{x}(t_{j})\;;\;t_{j-1},d\mathbf{x}(t_{j-1})\right]exp\left\{\frac{-i}{\hbar}\sum_{j=1}^{k}\int_{t_{j-1}}^{t_{j}}\mathbf{E}[\tau'_{j},s]\boldsymbol{H}_{1,\rho}(s)ds\right\}. \end{aligned}$$

ON THE PHYSICAL AND MATHEMATICAL FOUNDATIONS OF QUANTUM PHYSICS VIA FUNCTIONAL INTEGRAES. If we include this in our candidate evolution operator  $\mathbf{U}^{\rho}_{\lambda}[t,a]\Phi(a)$  and compute the limit, we have

$$\mathbf{U}^{\rho}[t,a]\Phi(a)$$

$$= \int_{\mathbb{R}^{n[t,a]}} \mathbb{K}_{\mathbf{f}}[\mathcal{D}\mathbf{x}(t);\mathbf{x}(a)] \exp\left\{ (-i/\hbar) \int_{a}^{t} \boldsymbol{H}_{1,\rho}(s) ds \right\} \Phi(a).$$

Since the limit as  $\rho \to \infty$  exists on the left, it defines the limit on the right.

6.7. **Examples.** Theorem 6.12 is somewhat abstract. The following example covers most of non-relativistic quantum theory.

**Theorem 6.13.** Let  $\Delta$  be the Laplacian on  $L^2[\mathbb{R}_I^n]$ , or some other Hilbert space,  $\mathcal{H}$  and let V be any potential such that  $H=(-\hbar^2/2)\Delta \oplus V$  generates a unitary group on  $\mathcal{H}$  (see remarks below). Using time as an index, the problem

$$(i\hbar)\partial\psi(\mathbf{x},t)/\partial t = \mathbf{H}(t)\psi(\mathbf{x},t), \ \psi(\mathbf{x},0) = \psi_0(\mathbf{x}),$$

has a a unique solution on  $\mathcal{FD}^2_{\otimes}$  with the extended Feynman-Kac representation.

Remark 6.14. We have used  $\oplus$  to allow a generalized definition of addition (i.e, Trotter-Kato). In fact, Kato has shown that V can be any self-adjoint generator and Goldstein has called it a generalized Lie sum (see [KA])

Our second example is due to Albeverio and Mazzucchi [AM]. It is like the first but we provide a different approach. Let  $\mathbb C$  be a completely symmetric positive definite fourth-order covariant tensor on  $L^2[\mathbb R_I^n]$ , let  $\Omega$  be a symmetric positive-definite  $n \times n$  matrix and let  $\lambda$  be a nonnegative constant. Then:

$$H = -\frac{\hbar^2}{2} \mathbf{\Delta} + \frac{1}{2} \mathbf{x} \Omega^2 \mathbf{x} + \lambda \mathbb{C}[\mathbf{x}, \mathbf{x}, \mathbf{x}, \mathbf{x}]$$

is known to be the generator for a unitary group on  $L^2[\mathbb{R}_I^n]$ . Albeverio and Mazzucchi [AM] prove that  $\bar{H}$  (closure) has a path integral representation as the analytic continuation (in the parameter  $\lambda$ ) of an infinite dimensional generalized oscillatory integral. (Their version of Feynman's path integral.)

Using the results of the previous sections, we can extend H to  $KS^2[\mathbb{R}_I^n]$ , which generates a unitary group. Let  $V = \frac{1}{2}\mathbf{x}\Omega^2\mathbf{x} + \lambda\mathbb{C}[\mathbf{x},\mathbf{x},\mathbf{x},\mathbf{x}]$  and  $V_\rho = V(I + \rho V^*V)^{-1/2}$ ,  $\rho > 0$ . We can prove that  $V_\rho$  is a bounded generator that converges to V. Since  $-\frac{\hbar^2}{2}\mathbf{\Delta}$  generates a unitary group,  $H_\rho = -\frac{\hbar^2}{2}\mathbf{\Delta} + V_\rho$  also generates one and converges to H. Let

$$\boldsymbol{H}(\tau) = (\mathop{\hat{\otimes}}_{t \geqslant s > \tau} \mathbf{I}_s) \otimes H \otimes (\mathop{\otimes}_{\tau > s \geqslant 0} \mathbf{I}_s),$$

then  $\boldsymbol{H}(t)$  generates a unitary group for each t and  $\boldsymbol{H}_{\rho}(t)$  converges to  $\boldsymbol{H}(t)$  on  $\mathcal{FD}_{\otimes}^{2}$ . We can now obtain:

$$\mathbf{U}[t, a]\Phi = \int_{\mathbb{R}_{I}^{n[t, a]}} \mathbb{K}_{\mathbf{F}}[\mathcal{D}\mathbf{x}(\tau); \mathbf{x}(a)] \exp\left\{-(i/\hbar) \int_{a}^{\tau} V(s) ds\right\} \Phi$$
$$= \lim_{\rho \to 0} \int_{\mathbb{R}_{I}^{n[t, a]}} \mathbb{K}_{\mathbf{F}}[\mathcal{D}\mathbf{x}(\tau); \mathbf{x}(a)] \exp\left\{-(i/\hbar) \int_{a}^{\tau} V_{\rho}(s) ds\right\} \Phi.$$

We refer to [GZ] for additional examples, including path integrals for kernels that are not perturbations of the Laplacian.

### 7. Concluding remarks

Modern physical theories of fundamental interactions need the space  $\Phi$  of all field histories over spacetime [DW4]. This is an infinite-dimensional manifold, and the framework of field theory makes it therefore compelling to look for an appropriate mathematical language. This can be obtained by the choice of a separable Banach space, following the beautiful and profound presentation of Geroch [GER], with a subsequent application of the results presented in the main body of our paper. The infinite-dimensional setting is still an unknown land for the majority of the physics community. For example, contraction of tensor arguments is not defined therein [GER], and hence all geometric invariants which contribute to ultraviolet divergences of gravity [THO, GOR] cannot be defined.

7.1. Manifolds modelled on Banach spaces. Within this framework, the idea is to pick out, from the rich structure of a Banach space, a particular type of structure called by Geroch local smoothness [GER]. As a first step, one has to introduce a mechanism by means of which structure can be carried from Banach spaces to other mathematical entities. If M is a set and E is a Banach space, an E-chart on M consists of a subset U of M jointly with a map ψ from U to E, such that ψ is one-to-one, and the image ψ(U) of U by ψ is open in E. A chart establishes therefore a one-to-one correspondence between a certain subset U of M and a certain open subset ψ(U) of E. It is this correspondence that carries structure from E to M.

As a second step, one needs the concept of agreement between two charts as regards their induced smoothness structures on M. Let (U, ψ) and (U 0 , ψ<sup>0</sup> ) be two E-charts on the set M. On the intersection V = U ∩U 0 , two smoothness structures are induced, one from ψ and the other from ψ 0 . The former defines a correspondence between V and the subset ψ(V ) of E, while the latter defines a correspondence between V and the subset ψ 0 (V ) of E.

In order to compare these smoothness structures, let us consider the map ψ <sup>0</sup>ψ −1 from ψ(V ) to ψ 0 (V ), and its inverse ψ  ψ 0−1 , from ψ 0 (V ) to ψ(V ). These maps describe the interaction between the E-charts (U, ψ) and (U 0 , ψ<sup>0</sup> ). At this stage, we have got rid of the manifold M, and we deal with maps between subsets of Banach spaces. Now a mathematical symbol p is fixed, either a non-negative integer or the ∞ symbol, and the charts (U, ψ) and (U 0 , ψ<sup>0</sup> ) on M are said to be C p -compatible if the images ψ(V ) and ψ 0 (V ) are both open subsets of E, and if the maps

$$\psi' \odot \psi^{-1} : \psi(V) \to E, \ \psi \odot {\psi'}^{-1} : \psi'(V) \to E$$

are both C <sup>p</sup> maps. The key role is played by the second condition. Instead of requiring that our maps preserve vector space or norm structure, we just require that they preserve C <sup>p</sup> differential structure. In such a way, a single type of structure is isolated.

A manifold modelled on a Banach space consists of a non-empty set M, a Banach space E, a symbol p, and a collection ζ of E-charts on M, in such a way that the following conditions are satisfied:

- (i) Any two charts in ζ are C p compatible.
- (ii) The charts in ζ cover the set M, i.e., every point of M lies in at least one of the U's-
- (iii) Any chart on M which is compatible with all the charts in ζ is itself an element of ζ.
- (iv) Given distinct points p and p <sup>0</sup> of M, there exist charts (U, ψ) and (U 0 , ψ<sup>0</sup> ) in ζ such that p lies in U and p 0 lies in U 0 , and such that there exists a ball B centred at ψ(p) in ψ(U) and a ball B<sup>0</sup> centred at ψ 0 (p 0 ) in ψ 0 (U 0 ), with the inverse images ψ −1 (B) and ψ 0−1 (B<sup>0</sup> ) having empty intersection in M. The charts are then said to separate points of M.

Condition (i) means that, whenever two charts in the collection ζ induce smoothness structures in the same region of M, these structures agree. Condition (ii) requires that smoothness structure has been induced over all of M. Condition (iii) makes sure that no additional structure has been induced on M. Last, condition (iv) rules out non-Hausdorff manifolds. These four conditions define the concept of C<sup>∞</sup> manifold M based on a Banach space E. The charts in ζ are said to be the admissible charts. A subset O of the manifold M is said to be open if, for every admissible chart (U, ψ), the image ψ(U ∩ O) is open in E. These are just the open sets for a topology on M.

7.2. Further perspectives. In the physics-oriented literature, the object of interest is the in-out amplitude with its functional integral formal representation

$$\langle \text{out} | \text{in} \rangle = \int e^{i(\dot{S}[\varphi] + \frac{1}{2}\omega_{\alpha\beta}K^{\alpha}[\varphi]K^{\beta}[\varphi])} (\det \widehat{G}(\varphi))^{-1} \dot{\mu}[\varphi][d\varphi],$$

where S˙ denotes the classical action S supplemented by all counterterms that are needed to render the amplitude finite [DW4]. It is the factor (detGb[ϕ])<sup>−</sup><sup>1</sup> that gives rise to all ghost loops in the loop expansion, and the exponent −1 is what makes the ghost field a fermion. Such a ghost arises entirely from the fibre-bundle structure of the space Φ of field histories, from the Jacobian of the transformation from fibre-adapted coordinates to local fields ϕ i [DW4].

If the action functional for the in-out amplitude can be expanded to quadratic order in gauge and ghost fields, it can be given a rigorous meaning as a Henstock integral, by virtue of the material presented from Sec. 2 to Sec. 6. However, as far as we can see, a functional integral for nonperturbative gauge theories is still beyond the present capabilities of the framework presented in our paper.

If the choice of infinite-dimensional manifolds modeled on separable Banach spaces is the right thing to do, we expect that this will be the starting point for further progress on functional integrals for quantum gauge theories (see also the remarkable monograph by Glimm and Jaffe [GL]).

For further literature on the framework that we have outlined, we refer the readers to the work in all References not cited so far.

### Acknowledgment

Giampiero Esposito is grateful to INDAM for membership.

## 8. Appendix

At the risk of some repetitions, this appendix is aimed at physics-oriented readers who would appreciate a self-contained description of concepts of real and functional analysis.

8.1. Completion of a normed space. For any normed vector space S it is always possible to build a Banach space S<sup>1</sup> in such a way that S contains a linear space Σ dense upon S<sup>1</sup> and equivalent to S. If S is not complete, every space S<sup>1</sup> as above is said to be completion of S. In order to prove the theorem [MIR], let E be the space having as elements the sequences {fk} of elements of S which fulfill the Cauchy condition, and let us introduce in E an equivalence relation R upon assuming by definition that {fk} and {f 0 k } are equivalent if

(38) 
$$\lim_{k \to \infty} ||f_k - f_k'|| = 0.$$

Let us set S<sup>1</sup> = E/R, with the letter X used to denote its elements. If {fk} and {f 0 k } are two representatives of the element X of S1, the sequences {cfk} and {cf<sup>0</sup> k } ON THE PHYSICAL AND MATHEMATICAL FOUNDATIONS OF QUANTUM PHYSICS VIA FUNCTIONAL INTEGRALS57 are, for all c ∈ K, equivalent to each other and hence represent the same element of S<sup>1</sup> that we denote by cX. Thus, if we consider another element Y of S1, upon denoting by h<sup>k</sup> and h 0 k two of its representatives, also the sequences {f<sup>k</sup> + hk} and {f 0 <sup>k</sup> + h 0 k } are equivalent to each other and hence represent the same element of S1, that we denote by X + Y .

By virtue of such definitions S<sup>1</sup> is a vector space. In particular, the origin of S<sup>1</sup> is represented by any sequence {ωk} such that

$$\lim_{k\to\infty}\|\omega_k\|=0.$$

Let us now show that in S<sup>1</sup> one can define a norm by setting

(39) 
$$||X||_{S_1} = \lim_{k \to \infty} ||f_k||_S,$$

where {fk} is an arbitrary representative of X. Since {fk} satisfies the Cauchy condition, and by virtue of the inequality

$$(40)$$

it follows that the sequence of norms {kfkkS} is convergent and hence it is legitimate to define the norm of X as we have done. Furthermore, from the definition of the equivalence relation, it follows that the S1-norm of X is independent of the particular representative chosen to define it. The properties of a norm

$$||X|| \ge 0$$
,  $||\lambda X|| = |\lambda| ||X||$ ,  $||X + Y|| \le ||X|| + ||Y||$ 

are proved immediately, including the fact that only the zero vector has vanishing norm.

Let us now prove the completeness of S1. For this purpose, we consider a sequence {Xn} of elements of S<sup>1</sup> which fulfills the Cauchy condition. Given ε > 0 one can thus find an index ν<sup>ε</sup> such that

$$n > \nu_{\varepsilon} \Longrightarrow \|X_{n+p} - X_n\|_{S_1} < \frac{\varepsilon}{2}, \ \forall p \in \mathbb{N}.$$

If we then denote by kf (n) k k a representative of Xn, one can write that

(41) 
$$n > \nu_{\varepsilon} \Longrightarrow \lim_{k \to \infty} \|f_k^{(n+p)} - f_k^{(n)}\|_S < \frac{\varepsilon}{2}, \ \forall p \in \mathbb{N}.$$

On the other hand, since the sequence {f (n) k } verifies the Cauchy condition, it is possible to associate to every n an index k<sup>n</sup> such that

(42) 
$$k > k_n \Longrightarrow ||f_k^{(n)} - f_{k_n}^{(n)}||_S < \frac{1}{n}.$$

One can also assume, without loss of generality, that kn+1 > kn. Let us then consider the sequence

$$(43) f_{k_1}^{(1)}, f_{k_2}^{(2)}, ..., f_{k_n}^{(n)}, ...$$

and let us show that it verifies the Cauchy condition. For this purpose we point out that, for all n, one finds for all values of k

$$||f_{k_{n+p}}^{(n+p)} - f_{k_n}^{(n)}||_S \le ||f_{k_{n+p}}^{(n+p)} - f_k^{(n+p)}||_S + ||f_k^{(n+p)} - f_k^{(n)}||_S + ||f_k^{(n)} - f_{k_n}^{(n)}||_S.$$

Thus, if one takes k > kn+<sup>p</sup> > kn, one obtains by virtue of (42)

$$||f_{k_{n+p}}^{(n+p)} - f_{k_n}^{(n)}||_S \le \frac{1}{(n+p)} + \frac{1}{n} + ||f_k^{(n+p)} - f_k^{(n)}||_S.$$

Furthermore, upon taking n > <sup>4</sup> ε , one finds

$$||f_{k_{n+p}}^{(n+p)} - f_{k_n}^{(n)}||_S \le \frac{\varepsilon}{2} + ||f_k^{(n+p)} - f_k^{(n)}||_S.$$

ON THE PHYSICAL AND MATHEMATICAL FOUNDATIONS OF QUANTUM PHYSICS VIA FUNCTIONAL INTEGRALS59

Since k can be chosen as large as we please, it follows from (41) that

$$n > \max\left(\frac{4}{\varepsilon}, \nu_{\varepsilon}\right) \Longrightarrow \|f_{k_{n+p}}^{(n+p)} - f_{k_n}^{(n)}\|_{S} \le \varepsilon, \ \forall p \in N.$$

We have therefore proved that the sequence (43) fulfills the Cauchy condition and hence it represents an element X of S1. We now want to prove that

$$\lim_{n \to \infty} X_n = X.$$

Indeed, one has

$$||X_n - X||_{S_1} = \lim_{r \to \infty} ||f_r^{(n)} - f_{k_r}^{(r)}||_S.$$

Upon fixing for the moment n > νε, let us choose ν anf s in such a way that

$$r > \max(n, k_n), \ s > k_r \ge k_n.$$

One finds therefore

$$||f_r^{(n)} - f_{k_r}^{(r)}||_S \le ||f_r^{(n)} - f_{k_n}^{(n)}||_S + ||f_{k_n}^{(n)} - f_s^{(n)}||_S$$
$$+ ||f_s^{(n)} - f_s^{(r)}||_S + ||f_s^{(r)} - f_{k_r}^{(r)}||_S.$$

k<sup>r</sup>

By virtue of this condition and of the majorization (42), it follows that

$$||f_r^{(n)} - f_{k_r}^{(r)}||_S \le \frac{2}{n} + ||f_s^{(n)} - f_s^{(r)}||_S + \frac{1}{r}.$$

On passing to the limit as s → ∞ one obtains, by virtue of the limit (41),

$$||f_r^{(n)} - f_{k_r}^{(r)}||_S \le \frac{2}{n} + \frac{1}{r} + \frac{\varepsilon}{2}.$$

The subsequent limit as r → ∞ yields

$$n > \nu_{\varepsilon} \Longrightarrow \|X_n - X\|_{S_1} \le \frac{2}{n} + \frac{\varepsilon}{2},$$

which in turn implies the limit (44). We have thus proved the completeness of S1.

If we now consider the space Σ of the X ∈ S<sup>1</sup> for which the sequence {fk} that represents X is convergent, we can consider the map

(45) 
$$x = S \cdot \lim_{k \to \infty} x_k \in S \longrightarrow X = \{f_k\} \in \Sigma,$$

that we agree to denote by X = T(x). Since by virtue of the limit (39) one can write that

$$||T(x)||_{S_1} = ||x||_S,$$

the equivalence of S and Σ has been proved.

Let us now prove that Σ is dense upon S1. For this purpose, let X be an element of S<sup>1</sup> while {fk} is one of its representatives. Let us set

$$\forall n \in N, \ h_n^{(k)} = f_k, \ Y^{(k)} = \{h_n^{(k)}\}_{n \in N}.$$

Since the sequences {h (k) <sup>n</sup> }n∈<sup>N</sup> are all constant, one has Y (k) ∈ Σ. On the other hand, by virtue of the Cauchy condition satisfied by X, one has

$$\forall \varepsilon > 0 \; \exists \nu_{\varepsilon} : n, k > \nu_{\varepsilon} \Longrightarrow \|f_n - h_n^{(k)}\|_S < \varepsilon$$

and hence, by passage to the limit as n → ∞, one finds

$$k > \nu_{\varepsilon} \Longrightarrow ||X - Y^{(k)}||_{S_1} > \varepsilon$$

i.e.

$$X = \lim_{k \to \infty} Y^{(k)}.$$

The last remaining task is to prove that, if S 0 1 is another Banach space containing a space Σ<sup>0</sup> equivalent to S and dense upon S 0 1 , the two spaces S<sup>1</sup> and S 0 <sup>1</sup> are equivalent. Indeed, if one defines an isomorphism between Σ<sup>0</sup> and Σ by establishing ON THE PHYSICAL AND MATHEMATICAL FOUNDATIONS OF QUANTUM PHYSICS VIA FUNCTIONAL INTEGRALS61 a correspondence between X<sup>0</sup> ∈ Σ and X ∈ Σ when they have the same image in S, one has

$$||X'||_{\Sigma'} = ||X||_{\Sigma}$$

and hence such an isomorphism is an equivalence between Σ and Σ<sup>0</sup> . Bearing in mind that Σ is dense upon S<sup>1</sup> and Σ<sup>0</sup> is dense upon S 0 1 , such an equivalence is extended by continuity to an equivalence between S<sup>1</sup> and S 0 1 .

- 8.2. Abstract measure theory. Following [AFP], let X be a nonempty set and let E be a collection of subsets of X.
- (i) The set E is said to be an algebra if the empty set ∅ ∈ E, the union E<sup>1</sup> ∪ E<sup>2</sup> ∈ E and the set-theoretic difference X − E<sup>1</sup> ∈ E whenever E<sup>1</sup> and E<sup>2</sup> belong to E.
- (ii) An algebra E is said to be a σ-algebra if, for any sequence (Eh) ⊂ E, its union ∪hE<sup>h</sup> belongs to E.
- (iii) For any collection G of subsets of X, the σ-algebra generated by G is the smallest σ-algebra containing G. If (X, τ ) is a topological space, one denotes by B(X) the σ-algebra of Borel subsets of X, i.e., the σ-algebra generated by the open subsets of X.
- (iv) If E is a σ-algebra in X, the pair (X, E) is said to be a measure space.

By virtue of De Morgan laws, algebras are closed under finite intersections, and σ-algebras are closed under countable intersections. Furthermore, since the intersection of any family of σ-algebras is a σ-algebra, the concept of generated σalgebra is meaningful. Sets endowed with a σ-algebra are the appropriate framework to introduce measures.

If (X, E) is a measure space, let m be a natural number ≥ 1.

(a) The function µ : E → R <sup>m</sup> is a measure if µ(∅) = 0 and, for any sequence (Eh) of pairwise disjoint elements of E, countable additivity holds, i.e.,

$$\mu\left(\bigcup_{h=0}^{\infty} E_h\right) = \sum_{h=0}^{\infty} \mu(E_h).$$

If m = 1, µ is said to be a real measure, whereas, if m > 1, µ is said to be a vector measure.

(b) If µ is a measure, its total variation |µ| for every E ∈ E is defined according to

$$|\mu|(E) \equiv \sup \left\{ \sum_{h=0}^{\infty} |\mu(E_h)| : \ E_h \in \mathcal{E} \ \text{pairwise disjoint}, \ E = \bigcup_{h=0}^{\infty} E_h \right\}.$$

(c) If µ is a real measure, its positive part µ <sup>+</sup> and negative part µ <sup>−</sup> are defined as follows:

$$\mu^+ = \frac{|\mu| + \mu}{2}, \ \mu^- = \frac{|\mu| - \mu}{2}.$$

In particular, if the set X is nonempty and E is the σ-algebra of all its subsets, one can define the Dirac measure on (X, E) as follows. To each x ∈ X we associate the measure δ<sup>x</sup> defined by δx(E) = 1 if x ∈ E, δx(E) = 0 otherwise. If (xh) is a sequence in X and if (ch) is a sequence in R <sup>m</sup> such that the series P h |ch| is convergent, we can set

$$\left(\sum_{h=0}^{\infty} c_h \delta_{x_h}\right)(E) = \sum_{x_h \in E} c_h$$

and obtain a R <sup>m</sup>-valued measure. Measures of this kind are said to be purely atomic. The set S<sup>µ</sup> of atoms of a measure µ in a measure space (X, E) is defined by

$$S_{\mu} = \{ x \in X : \ \mu(\{x\}) \neq 0 \},$$

ON THE PHYSICAL AND MATHEMATICAL FOUNDATIONS OF QUANTUM PHYSICS VIA FUNCTIONAL INTEGRALS63 provided that the singletons {x} are elements of E. If µ is finite or σ-finite, the set of atoms is at most countable.

8.3. Functions of bounded variation. If Ω is a generic open set in R <sup>N</sup> , and if u ∈ L 1 (Ω), we say that u is a function of bounded variation in Ω if the distributional derivative of u is representable by a measure in Ω, i.e., if [AFP]

(46) 
$$\int_{\Omega} u \frac{\partial \phi}{\partial x_i} dx = -\int_{\Omega} \phi dD_i u, \ \forall \phi \in C_c^{\infty}(\Omega), \ i = 1, ..., N$$

for some R <sup>N</sup> − valued measure Du = (D1u, ..., D<sup>N</sup> u) in Ω. The vector space of all functions of bounded variation in Ω is denoted by BV (Ω).

A smoothing argument shows that the integration by parts just written is still true for any φ ∈ C 1 c (Ω), or even for Lipschitz functions φ with compact support in Ω. All these formulae can be written concisely in the form

(47) 
$$\int_{\Omega} u \operatorname{div} \varphi \, dx = -\sum_{i=1}^{N} \int_{\Omega} \varphi_i dD_i u \, \forall \varphi \in [C_c^1(\Omega)]^N$$

The same notation can be used also for functions u ∈ [BV (Ω)]m. In such a case, Du is a m × N matrix of measures Diu <sup>α</sup> in Ω satisfying

(48) 
$$\sum_{\alpha=1}^{m} \int_{\Omega} u^{\alpha} \operatorname{div} \varphi^{\alpha} dx = -\sum_{\alpha=1}^{m} \sum_{i=1}^{n} \int_{\Omega} \varphi_{i}^{\alpha} dD_{i} u^{\alpha}, \ \forall \varphi \in [C_{c}^{1}(\Omega)]^{mN}.$$

The Sobolev space W<sup>1</sup>,<sup>1</sup> (Ω) is contained in BV (Ω), and the inclusion is strict. If u belongs to [BVloc(Ω)]<sup>m</sup>, one can prove the following properties:

(a) If the distributional derivative Du vanishes, u is equivalent to a constant in any connected component of Ω.

(b) For any bounded Lipschitz function ψ : Ω → R, the function uψ belongs to [BVloc(Ω)]<sup>m</sup> and

$$D(u\psi) = \psi Du + (u \otimes \nabla \psi) \mathcal{L}^N,$$

where ∇ is the approximate pointwise differential and L <sup>N</sup> denotes N-dimensional Lebesgue measure.

(c) If ρ is any convolution kernel and

$$\Omega_{\varepsilon} = \{ x \in \Omega : \operatorname{dist}(x, \partial \Omega) > \varepsilon \},$$

then

$$\nabla(u\star\rho_{\varepsilon})=Du\star\rho_{\varepsilon}\text{ in }\Omega.$$

### References

- [AM] S. Albeverio and S. Mazzucchi Feynman path integrals for polynomially growing potentials, J. Funct. Anal. 221 (2005), 83-121.
- [AFP] L. Ambrosio, N. Fusco and D. Pallara, Functions of Bounded Variation and Free Discontinuity Problems, Clarendon Press, Oxford, (2000).
- [ACM] L. Ambrosio, A. Carlotto and A. Massaccesi, Lectures on Elliptic Partial Differential Equations, Scuola Normale Superiore, Pisa, (2018).
- [BE] C. Becchi, A. Rouet and R. Stora, Renormalization of the Abelian Higgs-Kibble model, Commun. Math. Phys. 42 (1975), 127-162.
- [DZ] R. O. Davies and Z. Schuss, A proof that Henstock's integral includes Lebesgue's, London Math. Soc. 2 (1970), 561-562.
- [DW1] B. S. DeWitt, Quantum theory of gravity. II. The manifestly covariant theory, Phys. Rev. 162 (1967), 1195-1239.
- [DW2] B. S. DeWitt, Quantum gravity: the new synthesis, in General Relativity, an Einstein Centenary Survey, eds. S. W. Hawking and W. Israel, Cambridge University Press, Cambridge (1979), 680-745.

- [DW3] B. S. DeWitt, The Space-time approach to quantum field theory, in Relativity, Groups and Topology II, eds. B. S. DeWitt and R. Stora, North Holland, Amsterdam (1984), 381-738.
- [DW4] B. S. DeWitt, The Global Approach to Quantum Field Theory, Clarendon Press, Oxford, (2003).
- [DI1] P. A. M. Dirac, The Lagrangian in quantum mechanics, Physikalische Zeitschrift der Sowjetunion 3 (1933), 64-72.
- [DI2] P. A. M. Dirac, The Principles of Quantum Mechanics, Clarendon Press, Oxford, (1958).
- [E] G. Esposito, Fondamenti di Teoria Classica dei Campi, Amazon, (2023).
- [EV] L. C. Evans, Partial Differential Equations, AMS Graduate Studies in Math. 18, Providence, R.I, (1998).
- [FAD] L. D. Faddeev and V. Popov Feynman diagrams for the Yang-Mills field, Phys. Lett. 25 B (1967), 29-30.
- [FA] J. R. Fanchi, Confronting the ENIGMA of TIME, World Scientific, Singapore, (2023).
- [FY1] R. P. Feynman, The Principle of Least Action in Quantum Mechanics, PhD. Dissertation, Physics. Princeton University, Princeton, NJ, (1942). (Available from University Microfilms Publications, No. 2948, Ann Arbor, Michigan.)
- [FY2] R. P. Feynman, Space-time approach to non-relativistic quantum mechanics, Rev. Mod. Phys. 20 (1948), 367-387.
- [FY3] R. P. Feynman, An operator calculus having applications in quantum electrodynamics, Phys. Rev. 84 (1951), 108-128.
- [FY4] R. P. Feynman, Quantum theory of gravitation, Acta Physica Polonica 24 (1963) 697-722.
- [FY5] R. P. Feynman, Quantum Electrodynamics, Benjamin, New York, (1964).
- [FH] R. P. Feynman and A. R. Hibbs, Quantum Mechanics and Path Integrals, McGraw-Hill, New York, (1965).
- [FMS] N. Fusco, P. Marcellini and C. Sbordone, Mathematical Analysis: Functions of Several Real Variables and Applications, Springer Nature, Berlin, (2023).
- [GER] R. Geroch, Infinite-Dimensional Manifolds, Minkowski Institute Press, Montreal, (2013).
- [GL] T. L. Gill, Banach spaces for the Schwartz distributions, Real Anal. Exch. 43(1) (2017), 1-36.

- [GM] T. L. Gill and T. Myers, Constructive analysis on Banach spaces, Real Anal. Exch. 44(1) (2019), 1-36.
- [GS] J . A. Goldstein, Semigroups of Linear Operators and Applications, Oxford University Press, New York, 1985
- [GZ09] T. L. Gill and W. W. Zachary, Banach spaces for the Feynman integral, Real Anal. Exch. 34 (2009), 1-44.
- [GZ] T. L. Gill and W. W. Zachary, Functional Analysis and the Feynman Operator Calculus, Springer, New York, (2016).
- [GL] J. Glimm and A. Jaffe, Quantum Physics, A Functional Integral Point of View, Springer, New York, (1987).
- [GO] R. A. Gordon, The Integrals of Lebesgue, Denjoy, Perron, and Henstock, AMS, Providence, (1994).
- [GOR] M. H. Goroff and A. Sagnotti, The ultraviolet behavior of Einstein gravity, Nicl. Phus. B 266 (1986), 709-736.
- [G] L. Gross, Abstract Wiener spaces, Proc. Fifth Berkeley Symposium on Mathematics Statistics and Probability, 1965, pp. 31-42. MR 35:3027.
- [GMG] M. Guadalupe Morales and R. Gait´an, Feynman's path integral seen as a Henstock integral, J. Phys.: Conf. Ser. 912 (2017) 012014.
- [HS1] R. Henstock, A Riemann-type integral of Lebesgue power, Can. J. Math. 20 (1968), 79-87.
- [HS2] R. Henstock, The General Theory of Integration, Clarendon Press, Oxford, (1991).
- [HP] E. Hille and R. S. Phillips, Functional Analysis and Semigroups, Amer. Math. Soc. Colloq. Pub. 31, Amer. Math. Soc. Providence, RI, (1957).
- [HP] L. P. Horwitz, W. C. Schieve and C. Piron, Gibbs ensembles in relativistic classical and quantum mechanics, Ann. Phys. (N.Y.) 137 (1981), 306-340.
- [HO2] L. P. Horwitz, On the significance of a recent experiment demonstrating quantum interference in time, Phys. Lett. A 355 (2006), 1-6.
- [J] F. Jones, Lebesgue Integration on Euclidean Space, Revised Edition, Jones and Bartlett Publishers, Boston, (2001).

- [KA] T. Kato Trotter's product formula for an arbitrary pair of self-adjoint contraction semigroups, in Topics in Functional Analysis. Advances in Mathematics: Supplementary Studies, vol. 3 (Academic, Press New York, 1978), pp. 185–195
- [KB] J. Kuelbs, Gaussian measures on a Banach space, Journal of Functional Analysis 5 (1970), 354-367.
- [KW] J. Kurzweil, Nichtabsolut konvergente Integrale, Teubner-Texte z¨ur Mathematik, Band 26, Teubner Verlagsgesellschaft, Leipzig, (1980).
- [KW1] J. Kurzweil, Generalized ordinary differential equations and continuous dependence on a parameter, Czech. Math. J. 7 (1957), 418-446.
- [L] P. D. Lax, Symmetrizable linear tranformations, Comm. Pure Appl. Math. 7 (1954), 633-647.
- [LL] E. H. Lieb and M. Loss, Analysis, AMS Graduate Studies in Math. 14, Providence, R.I, (1997).
- [MIR] C. Miranda, Istituzioni di Analisi Funzionale Lineare, Unione Matematica Italiana, Bologna, (1978).
- [MI] C. W. Misner, Feynman quantization of general relativity, Rev. Mod. Phys. 29 (1957), 497- 509.
- [RO] H. L. Royden, Real Analysis, (2nd Ed.) Macmillan Press, New York, (1968).
- [RU] W. Rudin, Functional Analysis, McGraw-Hill Press, New York, (1973).
- [RUD] W. Rudin, Principles of Mathematical Analysis, Third edition, McGraw-Hill Press, New York, (1976).
- [SK] S. Saks, Theory of the Integral, Dover Publications, New York, (1964).
- [SH] I. A. Shishmarev, On the Cauchy problem and T-products for hypoelliptic systems, Math. USSR Izvestiya 20 (1983), 577-609.
- [ST] V. Steadman, Theory of operators on Banach spaces, Ph.D thesis, Howard University, (1988).
- [STE] E. Stein, Singular Integrals and Differentiability Properties of Functions, Princeton University Press, Princeton, New Jersey, (1970).
- [SW] R. F. Streater and A. S. Wightman, PCT, Spin and Statistics, and All That, Benjamin, New York, (1964).

- [SU] V. N. Sudakov, Linear sets with quasi-invariant measure, Dokl.Akad.Nauk SSSR; 127 (1959), 524-525 (in Russian).
- [THO] G. 't Hooft and M. Veltman, One-loop divergences in the theory of gravitation, Ann. Inst. Henri Poincar´e 20 (1974) 69-94.
- [TY] L. Tuo-Yeong, Some full descriptive charaterizations of the Henstock-Kurweil integral in Euclidean space, Czechoslovak Math. Journal 55 (2005), 625-637.
- [TY1] L. Tuo-Yeong, Henstock-Kurweil Integration on Euclidean Spaces, Series in Real Analysis-Vol 12 World Scientific, New Jersey, (2011).
- [VN1] J. von Neumann, The uniqueness of Haar's measure, Rec. Math. (Mat. Sbornik) N.S., 1 (1936), 721-734.
- [VN2] J. von Neumann, On infinite direct products, Compositio Mathematica, 6 (1938), 1-77.
- [WE] A. Weil, L'int´egration dans les groupes topologiques et ses applications, Actualit´es Scientifiques et Industrielles, no. 869, Paris, 1940.
- [WH] J. A. Wheeler and R. P. Feynman, Classical electrodynamics in terms of direct interparticle action, Rev. Mod. Phys. 21 (1949), 425-433.
- [W] A. S. Wightman, How it was learned that quantized fields are operator-valued distributions, Fortschr. Phys. 44 (1996), 143-178.
- [YA] Y. Yamasaki, Measures on Infinite-Dimensional Spaces, World Scientific, Singapore, (1985).
- [YH] K. Yosida and E. Hewitt, Finitely additive measures, Trans. Amer. Math. Soc. 72 (1952), 46-66.

### Declaration

The authors certify that:

- (1) they have no relevant financial or non-financial interests to disclose;
- (2) they have no conflicts of interest to declare that are relevant to the content of this manuscript;
- (3) they have no financial or proprietary interests in any material discussed in this manuscript; and

ON THE PHYSICAL AND MATHEMATICAL FOUNDATIONS OF QUANTUM PHYSICS VIA FUNCTIONAL INTEGRALS69

(4) they have no affiliations with or involvement in any organization or entity with any financial interest or non-financial interest in the subject matter or materials discussed in this manuscript.

(Giampiero Esposito) Dipartimento di Fisica "Ettore Pancini", Complesso Universitario di Monte S. Angelo, Via Cintia Edificio 6, 80126 Napoli, Italy, INFN Sezione di Napoli, Complesso Universitario di Monte S. Angelo, Via Cintia Edificio 6, 80126 Napoli, Italy

(Tepper L. Gill) Department of EECS, Mathematics and Computational Physics Laboratory, Howard University, Washington DC 20059, USA