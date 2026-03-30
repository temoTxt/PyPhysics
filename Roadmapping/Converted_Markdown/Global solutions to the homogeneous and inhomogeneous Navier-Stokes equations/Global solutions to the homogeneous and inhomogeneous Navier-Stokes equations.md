T. L. GILL, D. WILLIAMS, AND W. W. ZACHARY\*

Abstract. In this paper we take a new approach to a proof of existence and uniqueness of solutions for the 3D-Navier-Stokes equations, which leads to essentially the same proof for both bounded and unbounded domains and for homogeneous or inhomogeneous incompressible fluids. Our approach is to construct the largest separable Hilbert space SD<sup>2</sup> [R<sup>3</sup> ], for which the Leray-Hopf (type) solutions in L<sup>2</sup> [R<sup>3</sup> ] are strong solutions in SD<sup>2</sup> [R<sup>3</sup> ]. We say Leray-Hopf type because our solutions are weak in the spatial sense but not in time.

When the body force is zero, we prove that, there exists a positive constant u+, such that, for all divergence-free vector fields in a dense set D contained in the closed ball <sup>B</sup> of radius (1−ε) 2 u+, 0 < ε < 1, the initial value problem has unique global weak solutions in C 1 ((0, ∞), B). When the body force is nonzero, we obtain the same result for vector fields in a dense set D contained in the annulus bounded by constants u<sup>−</sup> and <sup>1</sup> 2 u+. In either case, we obtain existence and uniqueness for the Leray-Hopf weak solutions on R 3 . Moreover, with mild conditions on the decay properties of the initial data, we obtain pointwise and time-decay of the solutions.

<sup>1991</sup> Mathematics Subject Classification. Primary (35Q30) Secondary(47H20), (76DO3) .

Key words and phrases. Global, 3D-Navier-Stokes Equations, homogeneous, inhomogeneous .

<sup>\*</sup>deceased.

## Introduction

Let [L 2 (R 3 )]<sup>3</sup> be the Hilbert space of square integrable functions on R 3 , let H[R 3 ] be the completion of the set of functions in <sup>u</sup> <sup>∈</sup> <sup>C</sup><sup>∞</sup> 0 [R 3 3 | ∇ · u = 0 which vanish at infinity with respect to the inner product of [L 2 (R 3 )]3 , and let V[R 3 ] be the completion of the above functions which vanish at infinity with respect to the inner product of H<sup>1</sup> [R 3 ], the functions in H[R 3 ] with weak derivatives in [L 2 (R 3 )]3 . The classical Navier-Stokes initial-value problem (on R <sup>3</sup> and all T > 0) is to find a function u : [0, T ] × R <sup>3</sup> <sup>→</sup> <sup>R</sup> <sup>3</sup> and <sup>p</sup> : [0, T ] <sup>×</sup> <sup>R</sup> <sup>3</sup> <sup>→</sup> <sup>R</sup> such that

$$\partial_t \mathbf{u} + (\mathbf{u} \cdot \nabla) \mathbf{u} - \nu \Delta \mathbf{u} + \nabla p = \mathbf{f}(t) \text{ in } (0, T) \times \mathbb{R}^3,$$

(1) 
$$\nabla \cdot \mathbf{u} = 0 \text{ in } (0,T) \times \mathbb{R}^3 \text{ (in the weak sense)},$$
 
$$\mathbf{u}(0,\mathbf{x}) = \mathbf{u}_0(\mathbf{x}) \text{ in } \mathbb{R}^3.$$

The equations describe the time evolution of the fluid velocity u(x, t) and the pressure p of an incompressible viscous homogeneous Newtonian fluid with constant viscosity coefficient ν in terms of a given initial velocity u0(x) and given external body forces f(x, t).

Let P be the (Leray) orthogonal projection of (L 2 [R 3 ])<sup>3</sup> onto H[R 3 ] and define the Stokes operator by: Au =: <sup>−</sup>P∆u, for <sup>u</sup> <sup>∈</sup> <sup>D</sup>(A) <sup>⊂</sup> <sup>H</sup><sup>2</sup> [R 3 ], the domain of A. If we apply P to equation (1), with B(u, u) = P(u · ∇)u, we can recast equation (1) into the standard form:

(2) 
$$\partial_t \mathbf{u} = -\nu \mathbf{A} \mathbf{u} - B(\mathbf{u}, \mathbf{u}) + \mathbb{P} \mathbf{f}(t) \text{ in } (0, T) \times \mathbb{R}^3,$$
$$\mathbf{u}(0, \mathbf{x}) = \mathbf{u}_0(\mathbf{x}) \text{ in } \mathbb{R}^3,$$

where the orthogonal complement of H relative to {L 2 (R 3 )} 3 , {v : v = ∇q, q ∈ H<sup>1</sup> [R 3 ]}, is used to eliminate the pressure term (see Galdi [GA] or [SY, T1,T2]).

**Background.** The existence of global weak solutions for (2) was proved by Leray [Le] for all divergence-free initial data  $\mathbf{u}_0 \in \mathbb{H}(\mathbb{R}^3)$ . (Hopf [Ho] solved the same problem for a bounded open domain  $\Omega \subset \mathbb{R}^n, n \geq 2$  (see also [Li1, T1, vW]).)

Leray used

(3) 
$$\mathbf{b}(\mathbf{u}, \mathbf{u}, \mathbf{u}) = \langle B(\mathbf{u}, \mathbf{u}), \mathbf{u} \rangle_{\mathbb{H}} = \int_{\mathbb{R}^3} [\mathbf{u}(\mathbf{x}) \cdot \nabla \mathbf{u}(\mathbf{x})] \cdot \mathbf{u}(\mathbf{x}) d\mathbf{x} = 0$$

to show that, for such initial data, the global solution  $\mathbf{u}(t, \mathbf{x})$  satisfies the well-known energy inequality:

$$\|\mathbf{u}(t)\|_{\mathbb{H}}^2 + 2\nu \int_0^t \|\mathbf{A}^{1/2}\mathbf{u}(s)\|_{\mathbb{H}}^2 ds \leqslant \|\mathbf{u}_0\|_{\mathbb{H}}^2, \quad \forall t \geqslant 0.$$

These solutions are called Leray-Hopf solutions. There are two open questions in this case. The first is whether or not all Leray-Hopf solutions are unique and the second is whether or not those solutions with smooth initial data are regular. (A weak solution  $\mathbf{u}(t, \mathbf{x})$  is regular if  $\|\mathbf{u}(t)\|_{\mathbb{V}}$  is continuous.)

Until 1964, another open problem was the existence of global-in-time strong solutions (in the  $\mathbb{H}$  norm) for the three-dimensional Navier-Stokes initial value problem. In that year, Fujita and Kato [FK] proved that strong, global-in-time, smooth three-dimensional solutions exist in the Sobolev space  $\mathbb{H}^{1/2}$  provided that the body forces are small and the initial data is small (compared to the viscosity term  $-\nu \mathbf{A}\mathbf{u}$ , see Section 3 and also [KF], [CH] and [T3]). Their work was extended to  $L^p$  spaces by F. Weissler [WE] and considered in Besov spaces (of negative index of regularity) by M. Cannone, Y. Meyer and F. Planchon [CMP].

Since then, a number of papers have appeared proving existence of global-intime, solutions for small initial data, which are strong in a particular norm. These results will be discussed briefly in Section 3 but, the interested reader is directed to [\[LE\]](#page-33-3) (see also [\[KA1\]](#page-33-4), [\[FRT\]](#page-32-1), [\[CA\]](#page-31-3), [\[PL\]](#page-34-2) and [\[KT\]](#page-33-5)). The authors of [\[GIP\]](#page-32-2) observe that, in Leray's theory the nonlinear term becomes zero because of equation (3). They note that, "For global existence to hold, the point in those "strong solution" theorems is that the smallest assumption enables one to get rid of the nonlinear term, which can be absorbed by the Laplacian."

The interesting paper of Chemin and Gallagher [\[CG\]](#page-31-4) discusses all the relevant spaces, provides their own approach to the problem, and gives a nice picture of the kind of results one can expect from the methods used. Many of the recent approaches exploit the interesting invariance properties of the Navier-Stokes equations (with zero body forces) to construct their spaces. However, it appears that these spaces do not maintain their invariance properties when body forces are present. Furthermore, as first observed by Kato [\[KA1\]](#page-33-4), strong global solutions in L <sup>n</sup>[R <sup>n</sup>] for example, are not necessarily weak solutions in the sense of Leray-Hopf (i.e., they need not have finite energy norm).

Purpose. Our approach differs sharply from other attempts. We first construct the largest separable Hilbert space, for which the Leray-Hopf (type) solutions on R 3 are strong ones. We then obtain the smallest viscosity and largest body forces that balance each other in such a way as to allow global strong solutions for reasonable velocities (see below). A major advantage is that, the methods developed also apply to bounded domains and inhomogeneous fluids (with minor adjustments).

Asymptotic Properties. A number of studies have been conducted on the asymptotic and stability behavior of solutions, u(t, x), of the Navier-Stokes equations. The paper by Brandolese and Vigneron [\[BV\]](#page-31-5) provides a comprehensive analysis of

the problem and a clear presentation of the latest results in this direction (see also [GIP] and Miyakawa [MI]).

The Problem. The general problem in a bounded and unbounded domain is closely related. However, there is one major difference in the two cases. In order to understand the nature of an additional difficulty for the unbounded domain, it is important to discuss a problem that occurs when the domain  $\Omega \subset \mathbb{R}^3$  is bounded.

In the bounded domain case, it has been shown by Brušlinskaja [BR] that, if the viscosity coefficient  $\nu$  is sufficiently small, then stationary solutions of (2) lose stability and at least one eigenvalue of the linear Stokes operator passes from the left halfplane to the right halfplane. This means that the stationary solution becomes a limit cycle, which may be either stable or unstable. In either case, these bifurcations result in the existence of nonunique solutions of the three-dimensional Navier-Stokes equations. This problem has been discussed in a more general way by Foias and Temam [FT]. From this, it's clear that the viscosity coefficient is not a passive constant, but plays an important role in determining the physical properties of the solution(s). In fact, if  $\lambda_1$  is the first eigenvalue of the Stokes operator, the quantity

$$G = \frac{\left\| \mathbf{f} \right\|_{\infty}^2}{\nu^2 \lambda_1^{3/4}},$$

known as the Grashof number, appears naturally and provides a measure of the dynamical complexity of solutions. The Grashof number is similar to the Reynolds number and the dynamical complexity increases with increasing G (see Foias et al [FMTT]). For both physical and mathematical reasons, the corresponding unbounded domain problem is more difficult to study and, to our knowledge, has not received attention in the literature.

0.0.1. *Stability of flow.* The natural requirement of stability of flow in R 3 is usually implemented with the requirement that, under reasonable conditions, the velocity vector field u(t) approaches zero for large t. However, this implies that the total energy of the fluid also approaches zero. Physically this means that the boundary at infinity is kept (at least) at the solid state phase point (i.e., zero degrees for water). Since the physical properties of our fluid change radically at this point, we must be slightly more precise in this case.

On the other hand, in turbulent flow, there is a very important difference between the two and three-dimensional case. In the two-dimensional case, the fluid kinetic energy is transferred from both large to small and small to large scales by nonlinear interactions between different scales of motion (see Fj¨ortoft [\[FJ\]](#page-31-7) and Thompson [\[TH\]](#page-35-2)). However, in the three-dimensional case, there is a one way nonlinear cascade from large to small scales of motion. It follows that, as the kinetic energy becomes large, the fluid velocity becomes more erratic at smaller and smaller scales. It is generally assumed that the condition for a nonturbulent flow is captured by the size of the viscosity coefficient. Physically, this is a good measure, but also implies a number of other well-defined conditions, usually related to body forces, temperature, pressure and relative domain configuration (i.e., obstacles, constrictions, etc). Thus, any reasonable solution should also lead to an upper bound on the total energy (i.e., the velocity in L 2 -norm). The current discussion implies some bounds on the body forces. We will be more precise later (see Theorem 25) but for now, it suffices to assume that f = supt∈R<sup>+</sup> kPf(t)k<sup>H</sup> < ∞.

Definition 1. *We say that a velocity vector field in* R 3 *is* reasonable *if for* 0 ≤ t < ∞*, there is a continuous function* m(t) > 0*, depending only on* t *and a constant*

 $M_0$ , which may depend on  $\mathbf{u}(0)$  and f, such that

$$0 < m(t) \leqslant \|\mathbf{u}(t)\|_{\mathbb{H}} \le M_0.$$

The above definition formalizes the requirement that the fluid has bounded positive definite energy, However, this condition still allows the velocity to approach zero at infinity in a weaker norm.

0.1. **Statement of Results.** Let  $\mathbf{SD}^2[\mathbb{R}^3]$  be our separable Hilbert space, which is constructed in Section 1.1 and contains  $[L^2(\mathbb{R}^3)]^3$  as a compact dense embedding. Let  $\mathbb{H}_{sd}$  be the completion in the  $\mathbf{SD}^2[\mathbb{R}^3]$  norm of the set of functions in  $\{\mathbf{u} \in \mathbb{C}_0^{\infty}[\mathbb{R}^3]^3 \mid \nabla \cdot \mathbf{u} = 0\}$  and let  $\mathbb{P}$  be the orthogonal projection of  $\mathbf{SD}^2[\mathbb{R}^3]$  onto  $\mathbb{H}_{sd}$ .

Rewrite the first equation in (2) in the form:

(4) 
$$\partial_t \mathbf{u} = \mathcal{A}(\mathbf{u}, t) \text{ in } (0, T) \times \mathbb{R}^3,$$
$$\mathcal{A}(\mathbf{u}, t) = -\nu \mathbf{A} \mathbf{u} - B(\mathbf{u}, \mathbf{u}) + \mathbb{P} \mathbf{f}(t).$$

Let  $\mathbb{B}$  be a closed convex subset of  $\mathbb{H}_{sd}$ , which will be identified during the proof of the following theorem.

Theorem 2. If  $\mathbf{f} \neq 0$ , for each  $t \in [0, \infty)$  there exist positive constants  $u_+$ ,  $u_-$ , depending on f and  $\nu$  such that, for all initial data  $\mathbf{u}_0 \in \mathbb{B} \cap D(\mathbf{A}) \subset \mathbb{H}_{sd}$  with  $0 \leq u_- < \|\mathbf{u}_0\|_{sd} \leq \frac{1}{2}u_+$ , the operator  $\mathcal{A}(\cdot,t)$  is the generator of a strongly continuous nonlinear contraction semigroup on  $\mathbb{B}$ . If  $\mathbf{f} = 0$ , we replace  $\frac{1}{2}u_+$  by  $\frac{(1-\varepsilon)}{2}u_+$ , where  $0 < \varepsilon < 1$ , and use the ball  $\mathbb{B}$  of radius  $\frac{(1-\varepsilon)}{2}u_+$  centered at the origin.

If T(t) is the nonlinear semigroup generated by  $\mathcal{A}(\cdot,t)$ , then  $\mathbf{u}(t,\mathbf{x}) = T(t)\mathbf{u}_0(\mathbf{x})$  solves the initial value problem (2). We now have:

**Theorem 3.** For each  $T \in \mathbf{R}^+$ ,  $t \in (0,T)$  and  $\mathbf{u}_0 \in \mathbb{B} \cap D(\mathbf{A})$ , the global-in-time Navier-Stokes initial-value problem in  $\mathbb{R}^3$ :

$$\partial_t \mathbf{u} + (\mathbf{u} \cdot \nabla) \mathbf{u} - \nu \Delta \mathbf{u} + \nabla p = \mathbf{0} \text{ in } (0, T) \times \mathbb{R}^3,$$

(5) 
$$\nabla \cdot \mathbf{u} = 0 \ in \ (0, T) \times \mathbb{R}^3,$$

$$\mathbf{u}(0,\mathbf{x}) = \mathbf{u}_0(\mathbf{x}) \text{ in } \mathbb{R}^3,$$

has a unique strong solution  $\mathbf{u}(t,\mathbf{x})$ , which is in  $\mathbf{SD}^2[[0,\infty);\mathbb{H}_{sd}]$ .

**Theorem 4.** If  $\mathbf{u}(t, \mathbf{x})$  is a solution in the sense of Kato [KA1], in any one of the  $\mathbf{L}^p[\mathbb{R}^3]$  spaces, or any space  $\mathcal{B}$ , which is continuously embedded in  $\mathbf{L}^p[\mathbb{R}^3]$ , then  $\mathbf{u}(t, \mathbf{x})$  is a solution in  $\mathbf{SD}^2[\mathbb{R}^3]$ .

Let S(t) be the semigroup generated by the Stokes operator. It is well-known that  $\lim_{t\to\infty} S(t)\mathbf{u}_0 = 0$ . Since any strong solution is a mild solution, we also have that

$$\mathbf{u}(t, \mathbf{x}) = S(t)\mathbf{u}_0(\mathbf{x}) + \int_0^t S(t - s) \left[ -B\left(\mathbf{u}(s, \mathbf{x}), \mathbf{u}(s, \mathbf{x})\right) + \mathbb{P}\mathbf{f}(s) \right] ds.$$

If we introduce the following energy matrices due to Brandolese and Vigneron,

$$\mathcal{E}_{h,k}(t) = \int_{\mathbb{R}^3} (u_h u_k)(\mathbf{x}, t) d\mathbf{x} \quad \text{and} \quad \mathcal{K}_{h,k}(t) = \int_0^t \int_{\mathbb{R}^3} (u_h u_k)(\mathbf{x}, s) d\mathbf{x} ds,$$

then we have:

**Theorem 5.** Let  $\mathbf{u}(t, \mathbf{x}) = T(t)\mathbf{u}_0(\mathbf{x})$  be the solution to the Navier-Stokes initial value problem (2). If  $\mathbf{f}(t) = \mathbf{0}$  and there is a  $\delta > 0$  such that

(6) 
$$ess \sup_{\mathbf{x} \in \mathbb{R}^3} (1 + |\mathbf{x}|)^{2+\delta} |\mathbf{A}\mathbf{u}_0(\mathbf{x})| < \infty,$$

then

(1) *There is a constant* γ *such that*

$$\mathbf{u}(t, \mathbf{x}) = S(t)\mathbf{u}_0(\mathbf{x}) + \gamma \nabla \left( \sum_{h,k} \frac{\delta_{h,k} |\mathbf{x}|^2 - 3x_h x_k}{3 |\mathbf{x}|^5} K_{h,k}(t) \right) + 0 \left( \frac{1}{|\mathbf{x}|^4} \right).$$

(2) *There is a constant* p<sup>0</sup> *such that*

$$p(t, \mathbf{x}) = p_0 - \gamma \sum_{h,k} \left( \frac{\delta_{h,k}}{3 |\mathbf{x}|^3} - \frac{x_h x_k}{|\mathbf{x}|^5} \right) E_{h,k}(t) + O_t \left( \frac{1}{|\mathbf{x}|^4} \right).$$

It follows from this that, assuming Theorem 4, equation (6) is sufficient to ensure stability of global solutions to the Navier-Stokes equations. Furthermore, Picard's iterative scheme clearly applies in this case.

It is known that, if <sup>u</sup><sup>0</sup> <sup>∈</sup> <sup>V</sup> and <sup>f</sup>(t) <sup>∈</sup> <sup>L</sup><sup>∞</sup>[(0, <sup>∞</sup>), <sup>H</sup>], then there is a time T > <sup>0</sup> such that a Leray-Hopf (weak) solution of the Navier-Stokes equation is uniquely determined on any subinterval of [0, T ) (see Sell and You, [\[SY\]](#page-34-4) p. 396). Thus, we also have that:

Corollary 6. *For each* <sup>t</sup> <sup>∈</sup> <sup>R</sup><sup>+</sup> *and* <sup>u</sup><sup>0</sup> <sup>∈</sup> <sup>B</sup><sup>∩</sup> <sup>D</sup>(A)*, the Navier-Stokes initial-value problem in* R 3 :

$$\partial_t \mathbf{u} + (\mathbf{u} \cdot \nabla) \mathbf{u} - \nu \Delta \mathbf{u} + \nabla p = \mathbf{f}(t) \text{ in } (0, T) \times \mathbb{R}^3,$$

(7) 
$$\nabla \cdot \mathbf{u} = 0 \ in \ (0, T) \times \mathbb{R}^3,$$

$$\mathbf{u}(0,\mathbf{x}) = \mathbf{u}_0(\mathbf{x}),$$

*has a unique weak solution* u(t, x)*, which is in* L 2 loc[[0, ∞); Hsd] *and in* L<sup>∞</sup> loc[[0, ∞); Vsd] ∩ C 1 [(0, ∞); Hsd]*. Moreover, in this case, we also have that* limt→∞ ku(t)k<sup>H</sup>sd = 0*. To be more precise,*

$$\lim_{t \to \infty} \|\mathbf{u}(t) - S(t)\mathbf{u}_0\|_{\mathbb{H}_{sd}} = O(t^{-\alpha/2}),$$

*where* 0 < α < 1/2*.*

It follows from here that existence in the sense of Leray-Hopf is suffcient to ensure asymptotic decay of global solutions to the Navier-Stokes equations.

0.2. **Summary.** In the first section, we bring together a number of basic analytic tools that we use to prove our main results. We then construct our Hilbert space and obtain strong a priori bounds for the nonlinear term in the Navier-Stokes equations. The second section is devoted to proofs of our main results. In the third section, we discuss how our approach allows us to solve the inhomogeneous problem (on  $\mathbb{R}^3$ ). Finally, we note that, with minor changes, our results also apply to the bounded domain case for both homogeneous and inhomogeneous fluids.

## 1. Basic Tools

We now establish a number of results that will be used in the sequel.

**Definition 7.** We say that a (generally nonlinear) operator  $A(\cdot,t)$  is (for each t)

- (1)  $\theta$ -Dissipative if  $\langle \mathcal{A}(\mathbf{u},t), \mathbf{u} \rangle_{\mathbb{H}_{sd}} \leq 0$ ,
- (2) Dissipative if  $\langle \mathcal{A}(\mathbf{u},t) \mathcal{A}(\mathbf{v},t), \mathbf{u} \mathbf{v} \rangle_{\mathbb{H}_{sd}} \leq 0$ ,
- (3) Strongly dissipative if there exists a  $\beta > 0$ , which may depend on t, such that

$$\langle \mathcal{A}(\mathbf{u},t) - \mathcal{A}(\mathbf{v},t), \mathbf{u} - \mathbf{v} \rangle_{\mathbb{H}_{sd}} \leq -\beta \|\mathbf{u} - \mathbf{v}\|_{\mathbb{H}_{sd}}^2$$
.

Note that, if  $\mathcal{A}(\cdot,t)$  is a linear operator, definitions (1) and (2) coincide. Theorem 9 below is essentially due to Browder [B], while Theorem 10 is a slight extension of one from Miyadera [M, p. 185, Theorem 6.20]. The extension follows from Theorem A1 and Theorem A2 of Crandall and Pazy [CP] along with the time-dependent version of the Crandall-Liggett Theorem [CL] (see the appendix to the

GLOBAL SOLUTIONS TO THE HOMOGENEOUS AND INHOMOGENEOUS NAVIER-STOKES EQUATIONS11 first section of [\[CL\]](#page-31-9)). Taken together, this is an extension of Theorems I and II in Kato [\[KA2\]](#page-33-7).

Theorem 8. *Let* B *be a closed, bounded, convex subset of* Hsd*. If* A(·, t) : D(A(·, t)) ∩ B → Hsd *is a densely defined strongly dissipative mapping for each fixed* t ∈ [0, ∞)*, then for each* λ > 0*,* Ran*[*I − λA(·, t)] ⊃ B*).*

Theorem 9. *Let* B *is a closed convex set and let* A(·, t), t ∈ I = [0, ∞) *be a densely defined family of operators on* Hsd *with domains* D(A(·, t)) ∩ B) = D*, independent of* t*, such that:*

- (1) *The operator* A(·, t) *is the generator of a contraction semigroup on* D *for each* t ∈ I*.*
- (2) *The function* A(u, t) *is continuous in both variables on* D × I*.*

*Then* A(·, t) *extends uniquely to the generator of a contraction semigroup on* B *and, for every* u<sup>0</sup> ∈ D ∩ B*, the problem* ∂tu(t, x) = A(u(t, x), t)*,* u(0, x) = u0(x)*, has a unique solution* u(t, x) ∈ C 1 (I; B)*.*

1.1. The Hilbert Space SD<sup>2</sup> . The purpose of this section is to construct a special class of functions in C<sup>∞</sup> c [R <sup>n</sup>] (i.e., functions that are infinitely differentiable with compact support). We will use this class to construct a separable HIlbert space, SD<sup>2</sup> [R n ], which contains Wk,p[R n ], for all k ∈ N and 1 ≤ p ≤ ∞.

**Definition 10.** For  $x \in \mathbb{R}$ ,  $0 \le y < \infty$  and  $1 < a < \infty$ , we define the Jones functions by g(x,y), h(x) by (see Jones [JO], p. 249):

$$g(x,y) = \exp\left\{-y^a e^{iax}\right\},$$

$$h(x) = \begin{cases} \int_0^\infty g(x,y)dy, & x \in \left[-\frac{\pi}{2a}, \frac{\pi}{2a}\right], \\ 0 & \text{otherwise}. \end{cases}$$

The following properties of g are easy to check:

(1)

$$\frac{\partial g(x,y)}{\partial x} = -iay^a e^{iax} g(x,y),$$

(2)

$$\frac{\partial g(x,y)}{\partial y} = -ay^{a-1}e^{iax}g(x,y),$$

so that

(3)

$$iy \frac{\partial g(x,y)}{\partial y} = \frac{\partial g(x,y)}{\partial x}.$$

It is also easy to see that h(x) is in  $L^1[-\frac{\pi}{2a}, \frac{\pi}{2a}]$  and,

(8) 
$$\frac{dh(x)}{dx} = \int_0^\infty \frac{\partial g(x,y)}{\partial x} dy = \int_0^\infty iy \frac{\partial g(x,y)}{\partial y} dy.$$

Integration by parts in the last expression of equation (8) shows that h'(x) = -ih(x), so that  $h(x) = h(0)e^{-ix}$  for  $x \in [-\frac{\pi}{2a}, \frac{\pi}{2a}]$ . Since  $h(0) = \int_0^\infty \exp\{-y^a\}dy$ , an additional integration by parts shows that  $h(0) = \Gamma(\frac{1}{a} + 1)$ . For each  $l \in \mathbb{N}$  let  $a = a_l = 3 \times 2^{l-1}$ ,  $h(x) = h_l(x)$ ,  $x \in [-\frac{\pi}{2a_l}, \frac{\pi}{2a_l}]$  and set  $\varepsilon_l = \frac{\pi}{4a_l}$ .

Let  $\mathbb{Q}$  be the set of rational numbers in  $\mathbb{R}$  and for each  $x^i \in \mathbb{Q}$ , define

$$f_l^i(x) = f_l(x - x^i) = \begin{cases} c_l \exp\left\{\frac{\varepsilon_l^2}{|x - x^i|^2 - \varepsilon_l^2}\right\}, & |x - x^i| < \varepsilon_l, \\ 0, & |x - x^i| \ge \varepsilon_l, \end{cases}$$

where  $c_l$  is the standard normalizing constant. It is easy to check that  $f_l^i \neq 0$  for  $-\varepsilon_l < x - x^i < \varepsilon_l$ , so that the support,  $\operatorname{spt}(f_l^i) \subset [-\varepsilon_l, \varepsilon_l] = [-\frac{\pi}{4a_l}, \frac{\pi}{4a_l}]$ .

Now set  $\chi_l^k(x) = (f_l^k * h_l)(x)$ , so that  $\operatorname{spt}(\chi_l^k) \subset [-\frac{\pi}{2^{l+1}}, \frac{\pi}{2^{l+1}}]$ . For  $x \in \operatorname{spt}(\chi_l^k)$ , we can also write  $\chi_l^k(x) = \chi_l(x - x^k)$  as:

$$\int_{-\infty}^{\infty} f_l[(x-x^k)-z)h_l(z)dz = \int_{-\infty}^{\infty} h_l[(x-x^k)-z]f_l(z)dz = e^{-i(x-x^k)}\int_{-\infty}^{\infty} e^{iz}f_l(z)dz.$$

It is easy to see that  $-\pi - |x^k| < x < \pi + |x^k|$ . Thus, if  $\alpha_l = \int_{-\infty}^{\infty} e^{iz} f_l(z) dz$  and  $I_k = \operatorname{spt}(\chi_l^k)$ , we can now define:

$$\xi_l^k(x) = \frac{1}{n\alpha_l} \frac{\chi_l^k[i(x)]}{3^{\pi + |x^k|}} = \frac{1}{n\alpha_l} \frac{\chi_l[i(x - x^k)]}{3^{\pi + |x^k|}} = \begin{cases} \frac{1}{n} \frac{e^{(x - x^k)}}{3^{\pi + |x^k|}}, & x \in I_k \\ 0, & x \notin I_k, \end{cases}$$

so that  $\left|\xi_l^k(x)\right| < \frac{1}{n}$ .

1.2. **The Space.** To construct our space on  $\mathbb{R}^n$ , let  $\mathbb{Q}^n$  be the set  $\{\mathbf{x} = (x_1, x_2 \cdots, x_n) \in \mathbb{R}^n\}$  such that  $x_j$  is rational for each j. Since this is a countable dense set in  $\mathbb{R}^n$ , we can arrange it as  $\mathbb{Q}^n = \{\mathbf{x}^1, \mathbf{x}^2, \mathbf{x}^3, \cdots\}$ . For each l and i, let  $\mathbf{B}_l(\mathbf{x}^i)$  be the closed cube centered at  $\mathbf{x}^i$  with edge  $\frac{\pi}{a_l}$  and diagonal of length  $r_l = \frac{\pi}{a_l} \sqrt{n}$ .

We now choose the natural order which maps  $\mathbb{N} \times \mathbb{N}$  bijectively to  $\mathbb{N}$ :

$$\{(1,1), (2,1), (1,2), (1,3), (2,2), (3,1), (3,2), (2,3), \dots\}.$$

Let  $\{\mathbf{B}_k, k \in \mathbb{N}\}$  be the resulting set of (all) closed cubes  $\{\mathbf{B}_l(\mathbf{x}^i) \mid (l,i) \in \mathbb{N} \times \mathbb{N}\}$  centered at a point in  $\mathbb{Q}^n$ . For  $\mathbf{x} \in \mathbf{B}_k$ , let

(9) 
$$\mathcal{E}_k(\mathbf{x}) \triangleq \left(\xi_l^i(x_1), \xi_l^i(x_2), \dots \xi_l^i(x_n)\right).$$

It is easy to see that  $\mathcal{E}_k(\mathbf{x})$  is in  $L^p[\mathbb{R}^n]^n \cap L^\infty[\mathbb{R}^n]^n$  for  $1 \leq p < \infty$ . Let  $L^p[\mathbb{R}^n]^n = \mathbf{L}^p[\mathbb{R}^n]$  and define  $F_k(\cdot)$  on  $\mathbf{L}^p[\mathbb{R}^n]$  by

(10) 
$$F_k(f) = \int_{\mathbb{R}^n} \mathcal{E}_k(\mathbf{x}) \cdot f(\mathbf{x}) d\mathbf{x}.$$

It is clear that  $F_k(\ \cdot\ )$  is a bounded linear functional on  $\mathbf{L}^p[\mathbb{R}^n]$  for each  $k,\ \|F_k\|_{\infty} \leq 1$ . Furthermore, if  $F_k(f)=0$  for all  $k,\ f=0$  (a.s.), so that  $\{F_k\}$  is fundamental on  $\mathbf{L}^p[\mathbb{R}^n]$  for  $1\leq p\leq \infty$ .

Set  $t_k = \frac{1}{2^k}$  so that  $\sum_{k=1}^{\infty} t_k = 1$  and define a new inner product  $(\cdot)$  on  $\mathbf{L}^2[\mathbb{R}^n]$  by

(11) 
$$(f,g) = \sum_{k=1}^{\infty} t_k \left[ \int_{\mathbb{R}^n} \mathcal{E}_k(\mathbf{x}) \cdot f(\mathbf{x}) d\mathbf{x} \right] \left[ \int_{\mathbb{R}^n} \mathcal{E}_k(\mathbf{y}) \cdot g(\mathbf{y}) d\mathbf{y} \right].$$

The completion of  $\mathbf{L}^2[\mathbb{R}^n]$ , with the above inner product, is also a Hilbert space,  $\mathbf{SD}^2[\mathbb{R}^n]$ .

Remark 11. This approach is related to the one used in [GZ2], to construct another Hilbert space. Here, one wanted to show that  $\mathbf{L}^1[\mathbb{R}^n]$  can be embedded in a Hilbert space which contains the Denjoy-integrable functions (i.e., functions for which  $|\int_{\mathbb{R}^n} f(\mathbf{x}) d\mathbf{x}| < \infty$ ) and was used to provide the first rigorous mathematical foundations for the Feynman path integral formulation of quantum mechanics (see also [GZ2]).

We recall that Alexiewicz [AL] has shown that the class,  $D(\mathbb{R})$ , of Denjoy integrable functions (restricted and wide sense) can be normed in the following manner: for  $f \in D(\mathbb{R})$ , define  $||f||_D$  by

(12) 
$$||f||_D = \sup_s \left| \int_{-\infty}^s f(r) dr \right|.$$

The restricted Denjoy integral is equivalent to the Henstock-Kurzweil integral (see [HS] and [KW]).

Replacing  $\mathbb{R}$  by  $\mathbb{R}^n$  in (11), for  $f \in D(\mathbb{R}^n)$ , we can also define a norm on  $D(\mathbb{R}^n)$ :

$$(13) \qquad \|f\|_D = \sup_{r>0} \left| \int_{\mathbf{B}_r} f(\mathbf{x}) d\mathbf{x} \right| = \sup_{r>0} \left| \int_{\mathbf{R}^n} \mathbf{I}_{\mathbf{B}_r}(\mathbf{x}) f(\mathbf{x}) d\mathbf{x} \right| < \infty,$$

where  $\mathbf{B}_r$  is any closed cube of diagonal r centered at the origin in  $\mathbb{R}^n$ , with sides parallel to the coordinate axes and  $\mathbf{I}_{\mathbf{B}_r}(\mathbf{x})$  is the indicator function of  $\mathbf{B}_r$ .

1.2.1. Functions of Bounded Variation. The objective of this section is to show that every HK-integrable function is in  $SD^2[\mathbb{R}^n]$ . To do this, we need to discuss a certain class of functions of bounded variation. For functions defined on  $\mathbb{R}$ , the definition of bounded variation is unique. However, for functions on  $\mathbb{R}^n$ ,  $n \geq 2$ , there are a number of distinct definitions.

**Definition 12.** A function  $f \in L^1[\mathbb{R}^n]$  is said to be of bounded variation in the sense of Cesari or  $f \in BV_c[\mathbb{R}^n]$ , if  $f \in L^1[\mathbb{R}^n]$  and each  $i, 1 \le i \le n$ , there exists a signed Radon measure  $\mu_i$ , such that

$$\int_{\mathbb{R}^n} f(\mathbf{x}) \frac{\partial \phi(\mathbf{x})}{\partial x_i} d\lambda_n(\mathbf{x}) = -\int_{\mathbb{R}^n} \phi(\mathbf{x}) d\mu_i(\mathbf{x}),$$

for all  $\phi \in \mathbb{C}_0^{\infty}[\mathbb{R}^n]$ .

This is the definition known to most analysts and is the standard one used in geometric measure theory and partial differential equations.

The class of functions of bounded variation in the sense of Vitali [YE], is well known to applied mathematicians and engineers interested in error estimates associated with research in financial derivatives, control theory, robotics, high speed networks and in the calculation of certain integrals. (See, for example [KAA], [NI], [PT] or [PTR] and references therein.)

For the general definition, see Yeong ([YE], p. 175). We present a definition that is sufficient for continuously differentiable functions.

**Definition 13.** A function f with continuous partials is said to be of bounded variation in the sense of Vitali or  $f \in BV_v[\mathbb{R}^n]$  if for all intervals  $[a_i, b_i]$ ,  $1 \le i \le n$ ,

$$V(f) = \int_{a_1}^{b_1} \cdots \int_{a_n}^{b_n} \left| \frac{\partial^n f(\mathbf{x})}{\partial x_1 \partial x_2 \cdots \partial x_n} \right| d\lambda_n(\mathbf{x}) < \infty.$$

**Definition 14.** We define  $BV_{v,0}[\mathbb{R}^n]$  by:

$$BV_{v,0}[\mathbb{R}^n] = \{ f(\mathbf{x}) \in BV_v[\mathbb{R}^n] : f(\mathbf{x}) \to 0, \text{ as } x_i \to -\infty \},$$

where  $x_i$  is any component of  $\mathbf{x}$ .

The following two theorems may be found in [YE]. (See p. 184 and 187, where the first is used to prove the second.) If  $[a_i, b_i] \subset \mathbb{R}$ , we define  $[\mathbf{a}, \mathbf{b}] \in \mathbb{R}^n$  by  $[\mathbf{a}, \mathbf{b}] = \prod_{k=1}^n [a_i, b_i]$ . (The notation (RS) means Riemann-Stieltjes.)

**Theorem 15.** Let f be HK-integrable on  $[\mathbf{a}, \mathbf{b}]$  and let  $g \in BV_{v,0}[\mathbb{R}^n]$ , then fg is HK-integrable and

$$(HK) \int_{[\mathbf{a}, \mathbf{b}]} f(\mathbf{x}) g(\mathbf{x}) d\lambda_n(\mathbf{x}) = (RS) \int_{[\mathbf{a}, \mathbf{b}]} \left\{ (HK) \int_{[\mathbf{a}, \mathbf{x}]} f(\mathbf{y}) d\lambda_n(\mathbf{y}) \right\} dg(\mathbf{x})$$

**Theorem 16.** Let f be HK-integrable on  $[\mathbf{a}, \mathbf{b}]$  and let  $g \in BV_{v,0}[\mathbb{R}^n]$ , then fg is HK-integrable and

$$\left| (HK) \int_{[\mathbf{a}, \mathbf{b}]} f(\mathbf{x}) g(\mathbf{x}) d\lambda_n(\mathbf{x}) \right| \le \|f\|_D V_{[\mathbf{a}, \mathbf{b}]}(g)$$

**Lemma 17.** The space  $D[\mathbb{R}^n]$ , of all HK-integrable functions is contained in  $SD^2[\mathbb{R}^n]$ .

*Proof.* Since each  $\mathcal{E}_k[\mathbf{x}]$  is a continuous and differentiable on its domain, for  $f \in D[\mathbb{R}^n]$ , from Theorem 16, we have:

$$||f||_{\mathbf{SD}^2}^2 = \sum_{k=1}^{\infty} t_k \left| \int_{\mathbb{R}^n} \mathcal{E}_k(\mathbf{x}) \cdot f(\mathbf{x}) d\mathbf{x} \right|^2 \leqslant \sup_k \left| \int_{\mathbb{R}^n} \mathcal{E}_k(\mathbf{x}) \cdot f(\mathbf{x}) d\mathbf{x} \right|^2$$
  
$$\leqslant ||f||_D^2 \left[ \sup_k V(\mathcal{E}_k) \right]^2 < \infty,$$

so that 
$$f \in \mathbf{SD}^2[\mathbb{R}^n]$$
.

1.2.2. Properties of  $\mathbf{SD}^2$ . We now discuss the general properties of  $\mathbf{L}^p[\mathbb{R}^n]$ . The first two parts of the following theorem are natural, but the last part is an unexpected benefit. It means that a weakly convergent sequence in any of the  $\mathbf{L}^p[\mathbb{R}^n]$  spaces is strongly convergent in  $\mathbf{SD}^2[\mathbb{R}^n]$ .

**Theorem 18.** For each  $p, 1 \leq p \leq \infty$ ,  $SD^2[\mathbb{R}^n] \supset L^p[\mathbb{R}^n]$  as a dense, continuous and compact embedding.

*Proof.* First, by construction,  $\mathbf{SD}^2[\mathbb{R}^n]$  contains  $\mathbf{L}^2[\mathbb{R}^n]$  densely, so we need only show that  $\mathbf{SD}^2[\mathbb{R}^n] \supset \mathbf{L}^q[\mathbb{R}^n]$  for  $q \neq 2$ . If  $f \in \mathbf{L}^q[\mathbb{R}^n]$  and  $q < \infty$ , we have

$$\begin{aligned} & \|f\|_{\mathbf{SD}^{2}} = \left\{ \sum_{k=1}^{\infty} t_{k} \left| \int_{\mathbb{R}^{n}} \mathcal{E}_{k}(\mathbf{x}) \cdot f(\mathbf{x}) d\mathbf{x} \right|^{2} \right\}^{1/2} \\ & \leqslant \left\{ \sum_{k=1}^{\infty} t_{k} \left[ \int_{\mathbb{R}^{n}} \left| \mathcal{E}_{k}(\mathbf{x}) \right|^{q} \cdot \left| f(\mathbf{x}) \right|^{q} d\mathbf{x} \right]^{\frac{2}{q}} \right\}^{1/2} \\ & \leqslant \sup_{k} \left\{ \left[ \int_{\mathbb{R}^{n}} \left| \mathcal{E}_{k}(\mathbf{x}) \right|^{q} \cdot \left| f(\mathbf{x}) \right|^{q} d\mathbf{x} \right]^{\frac{1}{q}} \right\} \leqslant \sup_{k} \left\| \mathcal{E}_{k} \right\|_{q} \left\| f \right\|_{q} \leqslant \left\| f \right\|_{q}. \end{aligned}$$

In the last term, we used  $\sup_{k} \|\mathcal{E}_{k}\|_{q} < 1$ , so that  $f \in \mathbf{SD}^{2}[\mathbb{R}^{n}]$ . if  $q = \infty$ , we have

$$||f||_{\mathbf{SD}^{2}} = \left\{ \sum_{k=1}^{\infty} t_{k} \left| \int_{\mathbb{R}^{n}} \mathcal{E}_{k}(\mathbf{x}) \cdot f(\mathbf{x}) d\mathbf{x} \right|^{2} \right\}^{1/2}$$

$$\leq \sup_{k} \left\{ \left( \int_{\mathbb{R}^{n}} |\mathcal{E}_{k}(\mathbf{x})| |f(\mathbf{x})| d\mathbf{x} \right) \right\} \leq \sup_{k} ||\mathcal{E}_{k}||_{1} ||f||_{\infty} \leq ||f||_{\infty},$$

since  $\|\mathcal{E}_k\|_{L^1} < 1$ .

The proof of compactness follows from the fact that, if  $\{f_n\}$  is any weakly convergent sequence in  $\mathbf{L}^p[\mathbb{R}^n]$ ,  $1 \leq p < \infty$  with limit f, then since  $\mathcal{E}_k(\mathbf{x}) \in \mathbf{L}^q[\mathbb{R}^n]$ ,  $1 \leq q \leq \infty$ ,

$$\int_{\mathbb{R}^n} \mathcal{E}_k(\mathbf{x}) \cdot [f_n(\mathbf{x}) - f(\mathbf{x})] d\mathbf{x} \to 0$$

for each k. Thus,  $\{f_n\}$  converges strongly to f in  $SD^2[\mathbb{R}^n]$ .

Finally, we note that  $\mathbf{SD}^2[\mathbb{R}^n] \supset \mathbf{L}^1[\mathbb{R}^n]^{**} = \mathfrak{M}[\mathbb{R}^n]$ , the space of finitely additive measures on  $\mathbb{R}^n$ . It follows that  $d\mu_k(\mathbf{x}) = \mathcal{E}_k(\mathbf{x})d\mathbf{x}$  defines an element in  $\mathfrak{M}[\mathbb{R}^n]$  (the dual space of  $\mathbf{L}^{\infty}[\mathbb{R}^n]$ ). Thus, if  $\{f_n\}$  is any weakly convergent sequence to f in  $\mathbf{L}^{\infty}[\mathbb{R}^n]$ ,  $\{f_n\}$  converges strongly to  $f \in \mathbf{SD}^2[\mathbb{R}^n]$ .

**Remark 19.** Since  $\mathbf{L}^{\infty}[\mathbb{R}^n] \subset \mathbf{SD}^2[\mathbb{R}^n]$ , while  $\mathbf{SD}^2[\mathbb{R}^n]$  is separable, we see in a clear and forceful manner that separability is not an inherited property.

#### GLOBAL SOLUTIONS TO THE HOMOGENEOUS AND INHOMOGENEOUS NAVIER-STOKES EQUATIONS19 Definition 20. *We call* SD<sup>2</sup> [R <sup>n</sup>] *the strong distribution Hilbert space for* R n*.*

In order to justify our definition, let α be a multi-index of nonnegative integers, α = (α1, α2, · · · αn), with |α| = P<sup>n</sup> <sup>j</sup>=1 α<sup>j</sup> . If D denotes the standard partial differential operator, let

$$D^{\alpha} = D^{\alpha_1} D^{\alpha_2} \cdots D^{\alpha_k}.$$

Theorem 21. *If* <sup>u</sup> <sup>∈</sup> SD<sup>2</sup> [R <sup>n</sup>] *and* <sup>D</sup>α<sup>u</sup> <sup>=</sup> <sup>v</sup><sup>α</sup> *in the weak sense, then* <sup>v</sup><sup>α</sup> <sup>∈</sup> SD<sup>2</sup> [R n]*.*

*Proof.* From our construction, each <sup>E</sup><sup>k</sup> <sup>∈</sup> <sup>C</sup><sup>∞</sup> c [R <sup>n</sup>], so that

$$\int_{\mathbb{R}^n} \mathcal{E}_k(\mathbf{x}) \cdot D^{\alpha} \mathbf{u}(\mathbf{x}) d\mathbf{x} = (-1)^{|\alpha|} \int_{\mathbb{R}^n} D^{\alpha} \mathcal{E}_k(\mathbf{x}) \cdot \mathbf{v}_{\alpha}(\mathbf{x}) d\mathbf{x}.$$

An easy calculation shows that, for any j, R <sup>R</sup><sup>n</sup> ∂jEk(x) · u(x)dx = R <sup>R</sup><sup>n</sup> Ek(x) · uα(x)dx, so that

$$\int_{\mathbb{R}^n} \mathcal{E}_k(\mathbf{x}) \cdot D^{\alpha} \mathbf{u}(\mathbf{x}) d\mathbf{x} = (-1)^{|\alpha|} \int_{\mathbb{R}^n} \mathcal{E}_k(\mathbf{x}) \cdot \mathbf{v}_{\alpha}(\mathbf{x}) d\mathbf{x}.$$

It now follows that, for any <sup>g</sup> <sup>∈</sup> SD<sup>2</sup> [R <sup>n</sup>], (D<sup>α</sup>u, g)SD<sup>2</sup> = (−1)<sup>|</sup>α<sup>|</sup> (vα, g)SD<sup>2</sup> , so that <sup>v</sup><sup>α</sup> <sup>∈</sup> SD<sup>2</sup> [R <sup>n</sup>].

The next result follows from Theorem 21 and explains our use of the term strong distribution in describing SD<sup>2</sup> [R n].

Corollary 22. *If* <sup>u</sup> *is in the domain of* <sup>D</sup>α*, then for any* <sup>g</sup> <sup>∈</sup> SD<sup>2</sup> [R <sup>n</sup>], (D<sup>α</sup>u, g)SD<sup>2</sup> = (−1)<sup>|</sup>α<sup>|</sup> (u, g)SD<sup>2</sup> *so that, in particular,* <sup>k</sup>D<sup>α</sup>ukSD<sup>2</sup> <sup>=</sup> kukSD<sup>2</sup> *.*

Recall that a function **u** is said to be in  $\mathbf{W}^{k,p}[\mathbb{R}^n]$ ,  $k \in \mathbb{N}$ ,  $1 \le p \le \infty$ , if

$$\|\mathbf{u}\|_{k,p} = \left\{ \sum_{0 \leqslant |\alpha| \leqslant k} \|D^{\alpha}\mathbf{u}\|_{\mathbf{L}^p}^p \right\}^{1/p} < \infty, \text{ if } 1 \leqslant p < \infty,$$

$$\|\mathbf{u}\|_{k,\infty} = \max_{0 \leqslant |\alpha| \leqslant k} \|D^{\alpha}\mathbf{u}\|_{\mathbf{L}^{\infty}} < \infty, \quad \text{if} \quad p = \infty.$$

**Lemma 23.** For any  $p, 1 \leq p \leq \infty$  and all  $k \in \mathbb{N}, \mathbf{W}^{k, p}[\mathbb{R}^n] \subset \mathbf{SD}^2[\mathbb{R}^n]$ .

*Proof.* If  $u \in \mathbf{W}^{k, p}[\mathbb{R}^n]$ , then, for any  $\alpha, 0 \leq |\alpha| \leq k$ , we have

$$||u||_{\mathbf{SD}^2}^p = ||D^{\alpha}u||_{\mathbf{SD}^2}^p \le \sum_{0 \le |\alpha| \le k} ||D^{\alpha}u||_{\mathbf{SD}^2}^p \le \sum_{0 \le |\alpha| \le k} ||D^{\alpha}u||_{\mathbf{L}^p}^p.$$

It follows that  $u \in \mathbf{SD}^2[\mathbb{R}^n]$ . The case of  $p = \infty$  is clear.

1.3. The Nonlinear Term: A Priori Estimates. The difficulty in proving the existence of global-in-time strong solutions for equation (4) is directly linked to the problem of getting good a priori estimates for the nonlinear term  $B(\mathbf{u}, \mathbf{u})$ .

**Theorem 24.** If **A** is the Stokes operator and  $\mathbf{u}(\mathbf{x},t) \in D(\mathbf{A})$  is a reasonable vector field, then

- (1)  $\langle -\nu \mathbf{A} \mathbf{u}, \mathbf{u} \rangle_{\mathbb{H}_{sd}} = -\nu \|\mathbf{A} \mathbf{u}\|_{\mathbb{H}_{sd}}^{2}$ .
- (2) For  $\mathbf{u}(\mathbf{x},t) \in \mathbf{SD}^2 \cap D(\mathbf{A})$  and each  $t \in [0,\infty)$ , there exists a constant  $M = M(\mathbf{u}(\mathbf{x},0)) > 0$ , such that

(14) 
$$\left| \langle B(\mathbf{u}, \mathbf{u}), \mathbf{u} \rangle_{\mathbb{H}_{sd}} \right| \le M \|\mathbf{u}\|_{\mathbb{H}_{sd}}^{3}.$$

(3)

(15) 
$$\left| \langle B(\mathbf{u}, \mathbf{v}), \mathbf{w} \rangle_{\mathbb{H}_{sd}} \right| \le M \|\mathbf{u}\|_{\mathbb{H}_{sd}} \|\mathbf{w}\|_{\mathbb{H}_{sd}} \|\mathbf{v}\|_{\mathbb{H}_{sd}}.$$

(4)

$$(16) \qquad \max\{\|B(\mathbf{u}, \mathbf{v})\|_{\mathbb{H}_{sd}}, \ \|B(\mathbf{v}, \mathbf{u})\|_{\mathbb{H}_{sd}}\} \leqslant M \|\mathbf{u}\|_{\mathbb{H}_{sd}} \|\mathbf{v}\|_{\mathbb{H}_{sd}}.$$

*Proof.* From equation (12), we have

$$\left\langle -\nu \mathbf{A}\mathbf{u}, \mathbf{u} \right\rangle_{\mathbb{H}_{sd}} = -\nu \sum\nolimits_{k=1}^{\infty} t_k \left[ \int_{\mathbb{R}^n} \mathcal{E}_k(\mathbf{x}) \cdot \mathbf{A}\mathbf{u}(\mathbf{x}) d\mathbf{x} \right] \left[ \int_{\mathbb{R}^n} \mathcal{E}_k(\mathbf{y}) \cdot \mathbf{u}(\mathbf{y}) d\mathbf{y} \right].$$

Using the fact that  $\mathbf{u} \in D(\mathbf{A})$  and that k = (l, i) (see equation (9)), so that

$$\int_{\mathbb{R}^n} \mathcal{E}_k(\mathbf{y}) \cdot \partial_{y_j}^2 \mathbf{u}(\mathbf{y}) d\mathbf{y} = \int_{\mathbb{R}^n} \partial_{y_j}^2 \mathcal{E}_k(\mathbf{y}) \cdot \mathbf{u}(\mathbf{y}) d\mathbf{y} 
= \int_{I_i} \partial_{y_j}^2 \left( \xi_l^i(y_1), \xi_l^i(y_2), \dots, \xi_l^i(y_n) \right) \cdot \mathbf{u}(\mathbf{y}) d\mathbf{y} = \int_{\mathbb{R}^n} \mathcal{E}_k(\mathbf{y}) \cdot \mathbf{u}(\mathbf{y}) d\mathbf{y}.$$

Using this in the above equation and summing on j, we have

$$\int_{\mathbb{R}^n} \mathcal{E}_k(\mathbf{y}) \cdot \mathbf{A} \mathbf{u}(\mathbf{y}) d\mathbf{y} = \int_{\mathbb{R}^n} \mathcal{E}_k(\mathbf{y}) \cdot \mathbf{u}(\mathbf{y}) d\mathbf{y}.$$

It follows that

$$\langle \mathbf{A}\mathbf{u}, \mathbf{u} \rangle_{\mathbb{H}_{sd}} = \sum\nolimits_{k=1}^{\infty} t_k \left[ \int_{\mathbb{R}^n} \mathcal{E}_k(\mathbf{x}) \cdot \mathbf{A}\mathbf{u}(\mathbf{x}) d\mathbf{x} \right] \left[ \int_{\mathbb{R}^n} \mathcal{E}_k(\mathbf{y}) \cdot \mathbf{A}\mathbf{u}(\mathbf{y}) d\mathbf{y} \right] = \|\mathbf{A}\mathbf{u}\|_{\mathbb{H}_{sd}}^2.$$

This proves (1). To prove (2), let  $\vec{\delta}(\mathbf{x}) = (\delta(x_1), \dots \delta_k(x_3))$ , the *n*-dimensional Dirac delta function and set  $\hat{\varepsilon} = \|\vec{\delta}(\mathbf{x})\|_{\mathbb{H}_{sd}}$ . We start with

$$b(\mathbf{u}, \mathbf{u}, \mathcal{E}_k) = \left| \langle B(\mathbf{u}, \mathbf{u}), \mathcal{E}_k \rangle_{\mathbb{H}_{sd}} \right| = \left| \int_{\mathbb{R}^3} \left( \mathbf{u}(\mathbf{x}) \cdot \nabla \right) \mathbf{u}(\mathbf{x}) \cdot \mathcal{E}_k(\mathbf{x}) d\mathbf{x} \right|$$

and integrate by parts, to get

$$\left| \int_{\mathbb{R}^3} \left\{ \sum_{i=1}^3 u_i(\mathbf{x})^2 \mathcal{E}_k^i(\mathbf{x}) d\mathbf{x} \right\} \right| \leqslant \sup_k \|\mathcal{E}_k\|_{\infty} \|\mathbf{u}\|_{\mathbb{H}}^2 \leq \|\mathbf{u}\|_{\mathbb{H}}^2.$$

Since **u** is reasonable, there is a constant  $\bar{M}$  depending on **u**(0) and f, such that  $\|\mathbf{u}\|_2^2 \leq \bar{M} \|\mathbf{u}\|_{\mathbb{H}_{sd}}^2$ . We now have

$$\begin{aligned} & \left| \langle B(\mathbf{u}, \mathbf{u}), \mathbf{u} \rangle_{\mathbb{H}_{sd}} \right| = \left| \sum_{k=1}^{\infty} t_k \left[ \int_{\mathbb{R}^3} \left( \mathbf{u}(\mathbf{x}) \cdot \nabla \right) \mathbf{u}(\mathbf{x}) \cdot \mathcal{E}_k(\mathbf{x}) d\mathbf{x} \right] \left[ \int_{\mathbb{R}^3} \mathbf{u}(\mathbf{y}) \cdot \mathcal{E}_k(\mathbf{y}) d\mathbf{y} \right] \right| \\ & \leq \bar{M} \hat{\varepsilon}^{-2} \left\| \mathbf{u} \right\|_{\mathbb{H}_{sd}}^2 \left| \sum_{k=1}^{\infty} t_k \left[ \int_{\mathbb{R}^3} \vec{\delta}(\mathbf{x}) \cdot \mathcal{E}_k(\mathbf{x}) d\mathbf{x} \right] \left[ \int_{\mathbb{R}^3} \mathbf{u}(\mathbf{y}) \cdot \mathcal{E}_k(\mathbf{y}) d\mathbf{y} \right] \right| \\ & \leq M \left\| \mathbf{u} \right\|_{\mathbb{H}_{sd}}^3, \end{aligned}$$

where  $M = \bar{M}\hat{\varepsilon}^{-1}$  and the third line above follows from Schwartz's inequality. The proofs of (3) and (4) are easy.

1.4. **Generation Theorem.** We now begin with a study of the operator  $\mathcal{A}(\cdot,t)$ , for fixed t, and establish conditions depending on  $\mathbf{A}$ ,  $\nu$ , and  $\mathbf{f}(t)$  which guarantee that  $\mathcal{A}(\cdot,t)$  generates a contraction semigroup. Clearly  $\mathcal{A}(\cdot,t)$  is defined on  $D(\mathbf{A})$  and, since  $\nu \mathbf{A}$  is a closed positive (m-accretive) operator,  $-\nu \mathbf{A}$  generates a linear contraction semigroup. Thus, we need to ensure that  $\mathcal{A}(\cdot,t)$  will be m-dissipative for each t. We assume that  $\mathbf{f}(t) \in L^{\infty}[[0,\infty); \mathbb{H}_{sd}]$  and is Hölder continuous in t, with  $\|\mathbf{f}(t) - \mathbf{f}(\tau)\|_{\mathbb{H}_{sd}} \leq a |t - \tau|^{\theta}$ , a > 0,  $0 < \theta < 1$ .

**Theorem 25.** If  $0 \neq f = \sup_{t \in \mathbf{R}^+} \| \mathbb{P} \mathbf{f}(t) \|_{\mathbb{H}_{sd}} < \infty$ , there exist positive constants  $u_+, u_-$ , depending only on f,  $\mathbf{A}$  and  $\nu$  such that, for all  $\mathbf{u}$  with  $0 < u_- \le \| \mathbf{u} \|_{\mathbb{H}_{sd}} \le \frac{1}{2} u_+$ ,  $\mathcal{A}(\cdot, t)$  is strongly dissipative.

If f = 0 this implies that  $u_{-} = 0$ . In this case, we replace  $\frac{1}{2}u_{+}$  by  $\frac{(1-\varepsilon)}{2}u_{+}$ , so that  $\mathcal{A}(\cdot,t)$  is strongly dissipative on  $0 < \|\mathbf{u}\|_{\mathbb{H}_{sd}} \leq \frac{(1-\varepsilon)}{2}u_{+}$ .

*Proof.* The proof of our assertion has two parts. First, for  $f \neq 0$ , we require that the nonlinear operator  $\mathcal{A}(\cdot,t)$  be 0-dissipative, which gives us an upper bound  $u_+$  and lower bound  $u_-$  in terms of the norm (i.e.,  $\|\mathbf{u}\|_{\mathbb{H}_{sd}} \leq u_+$ ). We then use this

part to show that  $\mathcal{A}(\cdot,t)$  is strongly dissipative on  $D(\mathbf{A}) \cap \mathbb{B}$ , for any closed convex set,  $\mathbb{B}$ , inside the annulus defined by  $\{\mathbf{u} \in D(\mathbf{A}) : 0 \le u_- \le \|\mathbf{u}\|_{\mathbb{H}_{sd}} \le \frac{1}{2}u_+\}$ . We then consider adjustments when f = 0.

Part 1) From equation (4), we consider the expression

$$\begin{split} & \langle \mathcal{A}(\mathbf{u},t),\mathbf{u} \rangle_{\mathbb{H}_{sd}} = -\nu \, \langle \mathbf{A}\mathbf{u},\mathbf{u} \rangle_{\mathbb{H}_{sd}} + \langle \left[ -B(\mathbf{u},\mathbf{u}) + \mathbb{P}\mathbf{f} \right],\mathbf{u} \rangle_{\mathbb{H}_{sd}} \\ & = -\nu \, \|\mathbf{A}\mathbf{u}\|_{\mathbb{H}_{sd}}^2 - \langle B(\mathbf{u},\mathbf{u}),\mathbf{u} \rangle_{\mathbb{H}_{sd}} + \langle \mathbb{P}\mathbf{f},\mathbf{u} \rangle_{\mathbb{H}_{sd}} \,. \end{split}$$

It follows that

$$\langle \mathcal{A}(\mathbf{u},t),\mathbf{u} \rangle_{\mathbb{H}_{sd}} \leq -\nu \|\mathbf{u}\|_{\mathbb{H}_{sd}}^2 + \bar{M} \|\mathbf{u}\|_{\mathbb{H}_{sd}}^3 + f \|\mathbf{u}\|_{\mathbb{H}_{sd}}.$$

Since  $\|\mathbf{u}\|_{\mathbb{H}_{sd}} > 0$ , we have that  $\mathcal{A}(\cdot, t)$  is 0-dissipative if

$$(17) \qquad -\nu \|\mathbf{u}\|_{\mathbb{H}_{sd}} + \bar{M} \|\mathbf{u}\|_{\mathbb{H}_{sd}}^2 + f \leqslant 0$$

If  $\mathbf{v} \in D(\mathbf{A})$  is a reasonable vector field satisfying (17), we let  $M = \sup \{\bar{M}_{\mathbf{v}}\}$ . If we solve inequality (17) using M, we get

$$u_{\pm} = \frac{\nu}{2M} \left\{ 1 \pm \sqrt{1 - (4fM)/(\nu)^2} \right\} = \frac{\nu}{2M} \left\{ 1 \pm \sqrt{1 - \gamma} \right\},$$

where  $\gamma = (4fM)/\nu^2$ . Since we want real distinct solutions, we must require that

(18) 
$$\gamma = \frac{4fM}{\nu^2} < 1 \implies 2\sqrt{fM} < \nu.$$

It is clear that, if  $\mathbb{P}\mathbf{f} \neq \mathbf{0}$ , then  $u_- < u_+$ , and our requirement that  $\mathcal{A}(\mathbf{u}, t)$  is 0-dissipative implies that, since our solution factors as  $(\|\mathbf{u}\|_{\mathbb{H}_{sd}} - u_+)(\|\mathbf{u}\|_{\mathbb{H}_{sd}} - u_-) \leq$  0, we must have that:

$$\|\mathbf{u}\|_{\mathbb{H}_{sd}} - u_+ \le 0, \|\mathbf{u}\|_{\mathbb{H}_{sd}} - u_- \ge 0.$$

It follows that, for  $u_{-} \leq \|\mathbf{u}\|_{\mathbb{H}_{sd}} \leq u_{+}$ ,  $\langle \mathcal{A}(\mathbf{u},t), \mathbf{u} \rangle_{\mathbb{H}_{sd}} \leq 0$ . (It is clear that, when  $\mathbb{P}\mathbf{f}(t) = \mathbf{0}$ ,  $u_{-} = 0$ , and  $u_{+} = \frac{\nu}{M}$ .)

Part 2): Now, for any  $\mathbf{u}, \mathbf{v} \in D(\mathbf{A})$  with  $\mathbf{u} - \mathbf{v} \in D(\mathbf{A})$  and

$$u_{-} \le \min(\|\mathbf{u}\|_{\mathbb{H}_{sd}}, \|\mathbf{v}\|_{\mathbb{H}_{sd}}) \le \max(\|\mathbf{u}\|_{\mathbb{H}_{sd}}, \|\mathbf{v}\|_{\mathbb{H}_{sd}}) \le (1/2)u_{+},$$

we have that

$$\langle \mathcal{A}(\mathbf{u},t) - \mathcal{A}(\mathbf{v},t), (\mathbf{u} - \mathbf{v}) \rangle_{\mathbb{H}_{sd}} = -\nu \| \mathbf{A}(\mathbf{u} - \mathbf{v}) \|_{\mathbb{H}_{sd}}^{2}$$

$$- \langle [B(\mathbf{u}, \mathbf{u} - \mathbf{v}) + B(\mathbf{v}, \mathbf{u} - \mathbf{v})], (\mathbf{u} - \mathbf{v}) \rangle_{\mathbb{H}_{sd}}$$

$$\leq -\nu \| \mathbf{u} - \mathbf{v} \|_{\mathbb{H}_{sd}}^{2} + M \| \mathbf{u} - \mathbf{v} \|_{\mathbb{H}_{sd}}^{2} (\| \mathbf{u} \|_{\mathbb{H}_{sd}} + \| \mathbf{v} \|_{\mathbb{H}_{sd}})$$

$$\leq -\nu \| \mathbf{u} - \mathbf{v} \|_{\mathbb{H}_{sd}}^{2} + M \| \mathbf{u} - \mathbf{v} \|_{\mathbb{H}_{sd}}^{2} u_{+}$$

$$= -\nu \| \mathbf{u} - \mathbf{v} \|_{\mathbb{H}_{sd}}^{2} + M \| \mathbf{u} - \mathbf{v} \|_{\mathbb{H}_{sd}}^{2} \left( \frac{\nu}{2M} \left\{ 1 + \sqrt{1 - \gamma} \right\} \right)$$

$$= -\frac{\nu}{2} \| \mathbf{u} - \mathbf{v} \|_{\mathbb{H}_{sd}}^{2} \left\{ 1 - \sqrt{1 - \gamma} \right\}$$

$$= -\sigma \| \mathbf{u} - \mathbf{v} \|_{\mathbb{H}_{sd}}^{2}, \ \sigma = \frac{\nu}{2} \left\{ 1 - \sqrt{1 - \gamma} \right\}.$$

If f = 0, in the computation above, we see that  $\sigma = 0$ . To obtain our result in this case, we replace  $\frac{1}{2}u_+$  by  $\frac{(1-\varepsilon)}{2}u_+$ . The same computation shows that  $\sigma = \nu\varepsilon$ 

Let  $\mathbb{B}$  be any closed convex set inside the annulus bounded by  $\frac{1}{2}u_+$  and  $u_-$ . The first part of the next Lemma follows easily from the properties of  $\mathbf{f}(t)$ , the second part follows from Part 2) of Theorem 24 (see equation 18), while the third part is trivial.

**Lemma 26.** Let  $t, \tau \in I = [0, \infty)$  and  $\mathbf{u}, \mathbf{v} \in D(\mathbf{A})$ . Then

- (1) The mapping  $\mathcal{A}(\mathbf{u},t)$  is Hölder continuous in t, with  $\|\mathcal{A}(\mathbf{u},t) \mathcal{A}(\mathbf{u},\tau)\|_{\mathbb{H}_{sd}} \leq a|t-\tau|^{\theta}$ , where a is the Hölder constant for the function  $\mathbf{f}(t)$ .
- (2) The mapping  $A(\mathbf{u},t)$  is a Lipschitz continuous function in  $\mathbf{u}$ , with

$$\|\mathcal{A}(\mathbf{u},t) - \mathcal{A}(\mathbf{v},t)\|_{\mathbb{H}_{sd}} \le b \|\mathbf{u} - \mathbf{v}\|_{\mathbb{H}_{sd}}$$
.

(3) The mapping  $A(\mathbf{u},t) - \mathbb{P}\mathbf{f}(t)$  is coercive:

$$\lim_{\left\|\mathbf{u}\right\|_{\mathbb{H}_{sd}}\rightarrow\infty}\frac{\left\langle \left(\mathcal{A}(\mathbf{u},t)-\mathbb{P}\mathbf{f}(t)\right),\mathbf{u}\right\rangle _{\mathbb{H}_{sd}}}{\left\|\mathbf{u}\right\|_{\mathbb{H}_{sd}}}=\infty.$$

**Theorem 27.** The operator  $\mathcal{A}(\cdot,t)$  is strongly dissipative and jointly continuous in  $\mathbf{u}$  and t. Furthermore, for each  $t \in \mathbf{R}^+$  and  $\beta > 0$ ,  $Ran[I - \beta \mathcal{A}(t)] \supset \mathbb{B}$ , so that  $\mathcal{A}(t)$  is m-dissipative on  $\mathbb{B}$ .

*Proof.* From Theorem 25,  $\mathcal{A}(\cdot,t)$  is strongly dissipative. A strongly dissipative operator is maximal dissipative, so that  $Ran[I-\beta\mathcal{A}(\cdot,t)]\supset \mathbb{B}$ . It follows from [CP] that, since  $\mathbb{H}_{\mathbb{H}_{sd}}$  is a Hilbert space,  $\mathcal{A}(\cdot,t)$  is m-dissipative on  $\mathbb{B}$  for each  $t\in \mathbf{R}^+$ .

To see that  $\mathcal{A}(\mathbf{u},t)$  is continuous in both variables, let  $\mathbf{u}_n, \mathbf{u} \in \mathbb{D}$ ,  $\|(\mathbf{u}_n - \mathbf{u})\|_{\mathbb{H}_{sd}} \to 0$ , with  $t_n, t \in I$  and  $t_n \to t$ . Using  $\|\mathbf{A}\mathbf{u}\|_{\mathbb{H}_{sd}} = \|\mathbf{u}\|_{\mathbb{H}_{sd}}$ , we have

$$\|\mathcal{A}(\mathbf{u}_{n}, t_{n}) - \mathcal{A}(\mathbf{u}, t)\|_{\mathbb{H}_{sd}} \leq \|\mathcal{A}(\mathbf{u}, t_{n}) - \mathcal{A}(\mathbf{u}, t)\|_{\mathbb{H}_{sd}} + \|\mathcal{A}(\mathbf{u}_{n}, t_{n}) - \mathcal{A}(\mathbf{u}, t_{n})\|_{\mathbb{H}_{sd}}$$

$$= \|[\mathbb{P}\mathbf{f}(t_{n}) - \mathbb{P}\mathbf{f}(t)]\|_{\mathbb{H}_{sd}} + \|\nu\mathbf{A}(\mathbf{u}_{n} - \mathbf{u}) + [B(\mathbf{u}_{n} - \mathbf{u}, \mathbf{u}_{n}) + B(\mathbf{u}, \mathbf{u}_{n} - \mathbf{u})]\|_{\mathbb{H}_{sd}}$$

$$\leq d \|t_{n} - t\|^{\theta} + \nu \|\mathbf{A}(\mathbf{u}_{n} - \mathbf{u})\|_{\mathbb{H}_{sd}} + \|B(\mathbf{u}_{n} - \mathbf{u}, \mathbf{u}_{n}) + B(\mathbf{u}, \mathbf{u}_{n} - \mathbf{u})\|_{\mathbb{H}_{sd}}$$

$$\leq d \|t_{n} - t\|^{\theta} + \nu \|(\mathbf{u}_{n} - \mathbf{u})\|_{\mathbb{H}_{sd}} + M \|(\mathbf{u}_{n} - \mathbf{u})\|_{\mathbb{H}_{sd}} \left\{ \|\mathbf{u}_{n}\|_{\mathbb{H}_{sd}} + \|\mathbf{u}\|_{\mathbb{H}_{sd}} \right\}$$

$$\leq d \|t_{n} - t\|^{\theta} + \nu \|(\mathbf{u}_{n} - \mathbf{u})\|_{\mathbb{H}_{sd}} + M \|(\mathbf{u}_{n} - \mathbf{u})\|_{\mathbb{H}_{sd}} \cdot H$$

It follows that A(u, t) is continuous in both variables.

When f = 0, B is the ball of radius (1−ε) 2 u+. We see from Theorem 9 that this completes the proof of Theorem 3.

Proof of Theorem 4: Theorem 3 allows us to conclude that, when u<sup>0</sup> ∈ D(A)∩ B, the initial value problem is solved and the solution u(t, x) is in C 1 [(0, ∞); B]. Since <sup>D</sup>(A) <sup>∩</sup> <sup>B</sup> <sup>⊂</sup> <sup>H</sup><sup>2</sup> , it follows that u(t, x) is also in V for each t > 0. However, we can only conclude that, for any T > 0,

$$\int_0^T \left\| \mathbf{u}(t) \right\|_{\mathbb{H}_{sd}}^2 dt < \infty, \text{ and } \sup_{0 < t < T} \left\| \mathbf{u}(t) \right\|_{\mathbb{V}_{sd}}^2 < \infty.$$

This condition is not strong enough to ensure that ku(t)k<sup>V</sup> is continuous (which is required to resolve the singularity question).

Proof of Theorem 5: The proof of Theorem 5 follows from the fact that H1/<sup>2</sup> and L p , <sup>1</sup> <sup>≤</sup> <sup>p</sup> ≤ ∞ are continuous and densely embedded in SD<sup>2</sup> . Furthermore, the statement concerning B is obvious.

Proof of Theorem 6: The assertions (1) and (2), follow from Theorem 1.2 in Brandolese and Vigneron [\[BV\]](#page-31-5).

The proof of the first part of Corollary 7 follows, since every weak (distributional) solution in L 2 [R 3 ] is a strong (norm) solution in SD<sup>2</sup> [R 3 ]. The second assertion is a special case of a result due to Kato [KA1] (see Remark 1.2 and Theorem 40 ).

## 2. The Inhomogeneous Problem

In the inhomogeneous case, equation (1) becomes (see [\[Li2\]](#page-33-11) and [\[GR\]](#page-32-6)):

$$\rho[\partial_t \mathbf{u} + (\mathbf{u} \cdot \nabla)\mathbf{u}] - \mu \Delta \mathbf{u} + \nabla p = \rho \mathbf{f}(t) \text{ in } (0, T) \times \mathbb{R}^3,$$

∇ · u = 0 in (0, T ) × R 3 (in the weak sense),

(20) 
$$\mathbf{u}(0, \mathbf{x}) = \mathbf{u}_0(\mathbf{x}) \text{ in } \mathbb{R}^3.$$
$$\frac{\partial \rho}{\partial t} + \mathbf{u} \cdot \nabla \rho = 0 \text{ in } (0, T) \times \mathbb{R}$$

$$\rho(0, \mathbf{x}) = \rho_0(\mathbf{x}) \text{ in } \mathbb{R}^3.$$

We assume that the initial density satisfies 0 ≤ ρ0(x) ≤ β for some constant 0 < β.

3 ,

If we use the Leray projection to eliminate the pressure and divide the first equation by ρ, we get

$$\partial_{t}\mathbf{u} + \mathbb{P}(\mathbf{u} \cdot \nabla)\mathbf{u} - \frac{\mu}{\rho}\mathbb{P}\Delta\mathbf{u} = \mathbb{P}\mathbf{f}(t) \text{ in } (0, T) \times \mathbb{R}^{3},$$

$$\mathbf{u}(0, \mathbf{x}) = \mathbf{u}_{0}(\mathbf{x}) \text{ in } \mathbb{R}^{3}.$$

$$(21)$$

$$\frac{\partial \rho}{\partial t} + \mathbf{u} \cdot \nabla \rho = 0 \text{ in } (0, T) \times \mathbb{R}^{3},$$

$$\rho(0, \mathbf{x}) = \rho_{0}(\mathbf{x}) \text{ in } \mathbb{R}^{3}.$$

The solution of the density equation is well-known, ρ(t, x) = U[t, 0]ρ0(x), where U[t, 0] is an isometry (which depends on u(t, x)). Using the Feynman operator calculus (see [\[GZ2\]](#page-32-4)), we can write it symbolically as U[t, 0]ρ0(x) = ρ0(x− R t 0 u(τ, x)dτ ). It follows that 0 ≤ ρ(t, x) ≤ β for all t ∈ [0, T ].

Equation (21) now becomes:

(22) 
$$\langle \mathcal{A}(\mathbf{u}, t), \mathbf{u} \rangle_{\mathbb{H}_{sd}} = -\mu \left\langle \frac{\mathbf{A}\mathbf{u}}{\rho}, \mathbf{u} \right\rangle_{\mathbb{H}_{sd}} + \langle [-B(\mathbf{u}, \mathbf{u}) + \mathbb{P}\mathbf{f}], \mathbf{u} \rangle_{\mathbb{H}_{sd}}$$

$$\leq -\frac{\mu}{\beta} \|\mathbf{A}\mathbf{u}\|_{\mathbb{H}_{sd}}^2 - \langle B(\mathbf{u}, \mathbf{u}), \mathbf{u} \rangle_{\mathbb{H}_{sd}} + \langle \mathbb{P}\mathbf{f}, \mathbf{u} \rangle_{\mathbb{H}_{sd}}.$$

It follows that, on setting  $\nu = \frac{\mu}{\beta}$ , we see that Theorems 4 and 5 also hold for the inhomogeneous problem with minor adjustments.

2.1. **Bounded Domains.** Finally, a close review of the results of this paper show that they also hold for both homogeneous and inhomogeneous flows on bounded domains, with only minor changes.

## 3. Relationship to Other Spaces

As noted earlier, Kato and fujita [KF] proved that the Navier-Stokes equations have global solutions for small initial data in  $\dot{H}^{1/2}$ . If  $\hat{h}(\mathbf{x}) = \int_{\mathbb{R}^3} e^{-2\pi i \mathbf{x} \cdot \mathbf{y}} h(\mathbf{y}) d\mathbf{y}$  is the Fourier transform of a function  $h \in L^2[\mathbb{R}^3]$ , then  $h \in \dot{H}^{1/2}[\mathbb{R}^3]$  if and only if

$$||h||_{\dot{H}^{1/2}} = \left[ \int_{\mathbb{R}^3} 2\pi \left| \mathbf{x} \right| \left| \hat{h}(\mathbf{x}) \right|^2 d\mathbf{x} \right]^{1/2} < \infty.$$

This is one of the spaces where the solutions (with zero body forces) are invariant under translation and scaling transformations. If  $\mathbf{u}(t, \mathbf{x})$  is a solution then  $\mathbf{u}_{\lambda}(t, \mathbf{x}) = \lambda \mathbf{u}(\lambda^2 t, \lambda \mathbf{x})$  is also one and  $\mathbf{u}_{0,\lambda}(\mathbf{x}) = \lambda \mathbf{u}_0(\lambda \mathbf{x})$ . (Kato proved a similar result in [KA1].)

The paper by Koch and Tataru [KT] is of special interest. They define the smallness of their solution in terms of the norm in the space  $BMO^{-1}$ . They consider this space natural for the problem, since it is an invariant space for translations and the scaling transformations,  $\mathbf{u}_{\lambda}(t,\mathbf{x}) = \lambda \mathbf{u}(\lambda^2 t, \lambda \mathbf{x})$ , that leave the Navier-Stokes equations invariant.

In their approach, Koch and Tataru [KT] use the subspace  $BMO^{-1}[\mathbb{R}^n]$  of  $BMO[\mathbb{R}^n]$ , the functions of bounded mean oscillation to construct strong solutions for the Navier-Stokes equations. In this section we study the relationship of their work to ours. (The main result is that  $BMO^{-1}[\mathbb{R}^n] \subset \mathbf{SD}^2[\mathbb{R}^n]$ .)

In order to do this, we first need to construct  $\mathbf{SD}^p[\mathbb{R}^n]$  for all p. For  $f \in \mathbf{L}^p[\mathbb{R}^n]$ , define:

$$||f||_{\mathbf{SD}^{p}} = \begin{cases} \left\{ \sum_{k=1}^{\infty} t_{k} \left| \int_{\mathbb{R}^{n}} \mathcal{E}_{k}(\mathbf{x}) \cdot f(\mathbf{x}) d\mathbf{x} \right|^{p} \right\}^{1/p}, 1 \leq p < \infty, \\ \sup_{k \geq 1} \left| \int_{\mathbb{R}^{n}} \mathcal{E}_{k}(\mathbf{x}) \cdot f(\mathbf{x}) d\mathbf{x} \right|, p = \infty. \end{cases}$$

It is easy to see that  $\|\cdot\|_{\mathbf{SD}^p}$  defines a norm on  $\mathbf{L}^p[\mathbb{R}^n]$ . If  $\mathbf{SD}^p[\mathbb{R}^n]$  is the completion of  $\mathbf{L}^p[\mathbb{R}^n]$  with respect to this norm, we have:

**Theorem 28.** For each  $q, 1 \leq q \leq \infty$ ,  $SD^p[\mathbb{R}^n] \supset L^q[\mathbb{R}^n]$  as continuous, dense and compact embeddings.

The proof of this result as well as the next is essentially the same as in Section 3 of [GZ1], so we omit them.

**Theorem 29.** For  $SD^p[\mathbb{R}^n]$ ,  $1 \le p \le \infty$ , we have:

- (1) If  $f, g \in \mathbf{SD}^p[\mathbb{R}^n]$ , then  $||f + g||_{\mathbf{SD}^p} \leqslant ||f||_{\mathbf{SD}^p} + ||g||_{\mathbf{SD}^p}$  (Minkowski inequality).
- (2) If  $1 and <math>p^{-1} + q^{-1} = 1$ , then the dual space of  $\mathbf{SD}^p[\mathbb{R}^n]$  is  $\mathbf{SD}^q[\mathbb{R}^n]$ .
- (3) The dual space of  $SD^1[\mathbb{R}^n]$  is  $SD^{\infty}[\mathbb{R}^n]$ .
- (4)  $\mathbf{SD}^{\infty}[\mathbb{R}^n] \subset \mathbf{SD}^p[\mathbb{R}^n]$ , as a continuous embedding, for  $1 \leq p < \infty$ .

### 4. Functions of Bounded Mean Oscillation

If Q is a cube in  $\mathbb{R}^3$  and g is locally  $\mathbf{L}^1[\mathbb{R}^3]$ , define  $m_{g,Q}$  by:

$$m_{g,Q} = \frac{1}{|Q|} \int_Q g(\mathbf{x}) d\mathbf{x}.$$

We call it the average or mean of g over Q.

**Definition 30.** If g is locally  $\mathbf{L}^1[\mathbb{R}^3]$ , we say that g is of bounded mean oscillation, and write  $g \in BMO[\mathbb{R}^3]$ , provided that:

$$\|g\|_{BMO} = \sup_{Q} \frac{1}{|Q|} \int_{Q} |g(\mathbf{x}) - m_{g,Q}| d\mathbf{x} < \infty,$$

where the sup is over all cubes in  $\mathbb{R}^3$  (see Grafakos [GRA], chapter 7).

It is well-known that  $\mathbf{L}^{\infty}[\mathbb{R}^3] \subset BMO[\mathbb{R}^3]$ , with the inclusion proper. For example,  $BMO[\mathbb{R}^3]$  also contains  $\ln |p(\mathbf{x})|$ , for all polynomials  $p(\mathbf{x})$ .

 $BMO^{-1}[\mathbb{R}^n]$  is a subspace of tempered distributions and the following theorem provides a nice characterization. (For a proof, see Koch and Tataru [KT].)

**Theorem 31.** Let the vector field u be a tempered distribution. Then  $u \in BMO^{-1}[\mathbb{R}^n]$  if and only if there exist  $f^i \in BMO[\mathbb{R}^n]$ , such that,  $u = \sum_{i=1}^n \partial_i f^i$ .

The following theorem shows that  $BMO^{-1}[\mathbb{R}^n] \subset \mathbf{SD}^{\infty}[\mathbb{R}^n]$ , so that necessarily,  $BMO^{-1}[\mathbb{R}^n] \subset \mathbf{SD}^2[\mathbb{R}^n]$ .

**Theorem 32.** If  $\mathbf{u} \in \mathcal{S}'[\mathbb{R}^n]$ , the space of tempered distributions, then  $\mathbf{u} \in \mathbf{SD}^{\infty}[\mathbb{R}^3]$ .

*Proof.* Since  $\mathcal{E}_k(\mathbf{x}) \in \mathbb{C}_c^{\infty}[\mathbb{R}^n]$  and of slow growth at infinity, it is in  $\mathcal{S}[\mathbb{R}^n]$ , so that,

$$\sup_{k} |\langle \mathbf{u}, \mathcal{E}_k \rangle| = \sup_{k} \left| \int_{\mathbb{R}^n} \mathbf{u}(\mathbf{x}) \cdot \mathcal{E}_k(\mathbf{x}) d\mathbf{x} \right| < \infty, \quad k \in \mathbb{N}.$$

It follows that <sup>u</sup> <sup>∈</sup> SD∞[<sup>R</sup> <sup>n</sup>].

4.1. Conclusion. In this paper we have introduced a new Hilbert space, which allows us to obtain uniqueness for the Leray-Hopf solutions on R 3 , with or without body forces. We also prove global-in-time strong solutions for the three-dimensional Navier-Stokes equations for both bounded and unbounded domains and for a homogeneous or inhomogeneous incompressible fluid. In addition, with mild conditions on the decay properties of the initial data, we obtain pointwise and time-decay of the solutions. However, our methods do not allow us to resolve the singularity question. Our space also contains the Kato solution and those in L p spaces. Although the space used by Koch-Tataru [\[KT\]](#page-33-5), BMO<sup>−</sup><sup>1</sup> <sup>⊂</sup> SD<sup>2</sup> , we are unable to ensure that the embedding is continuous. Thus, we are not able to show that solutions in their sense are solutions in SD<sup>2</sup> .

This paper replaces an earlier one, which contained a fatal error that could not be fixed in the manner we had hoped (see [\[GZ3\]](#page-32-9)).

## References

- <span id="page-30-2"></span>[AD] R. A. Adams, Sobolev Spaces, Academic Press, New York, 1975.
- [AL] A. Alexiewicz, Linear functionals on Denjoy-integrable functions, Colloq. Math. 1 (1948), 289-293.
- <span id="page-30-1"></span>[B] F. Browder, Nonlinear operators and nonlinear equations of evolution in Banach spaces, Proc. Sympos. Pure Math., Vol. 18 part II, Amer. Math. Soc., Providence, RI, 1970.
- <span id="page-30-0"></span>[BR] N. N. Bru˘slinskaja, The behavior of solutions of the equations of hydrodynamics when the Reynolds number passes through a critical point, Dokl. Akad. Nauk SSSR

- 162 (1965), 731-734 (Russian); English translation: Sovet Math. Dokl. 6 (1965), 724-728.
- <span id="page-31-5"></span>[BV] L. Brandolese and F. Vigneron, New asymptotic profiles of nonstationary solutions of the Navier-Stokes system, Journal de Math´ematiques Pure et Appliqu´es 88 (2007), 64-86.
- <span id="page-31-3"></span>[CA] M. Cannone, Ondelettes, paraproduits et NavierStokes, Diderot ´editeur, Arts et Sciences, (1995).
- [CF] P. Constantin and C. Foias, Navier-Stokes Equations, University of Chicago Press, Chicago, ILL, 1988.
- <span id="page-31-4"></span>[CG] J. Y. Chemin and I. Gallagher, Wellposedness and stability results for the Navier-Stokes equations in R 3 , Annales de l'Institut H. Poincar´e Analyse non Lin´eaire, 26 (2009), 599-624.
- <span id="page-31-1"></span>[CH] J. Y. Chemin, Remarques sur l'existence globale pour le systeme de Navier-Stokes incompressible, Siam. J. Math. Anal. 23 (1992), 20-28.
- <span id="page-31-9"></span>[CL] M. Crandall and T. Liggett, Generation of semigroups of nonlinear transformations on general Banach spaces, Amer. J. Math. 93 (1971), 265-293.
- <span id="page-31-2"></span>[CMP] M. Cannone, Y. Meyer et F. Planchon, Solutions autosimilaires des ´equations de Navier-Stokes, S´eminaire "Equations aux D´eriv´ees Partielles" de l' ´ Ecole polytech- ´ nique, Expos´e VIII, 1993-1994.
- <span id="page-31-8"></span>[CP] M. Crandall and A. Pazy, Semigroups of nonlinear contractions and dissipative sets, J. Functional Analysis 3 (1969), 376-418.
- <span id="page-31-7"></span>[FJ] R. Fj¨ortoft, On the change in the spectral distribution of kinetic energy for twodimensional non-divergent flow, Tellus, 5 (1953), 225-230.
- <span id="page-31-0"></span>[FK] H. Fujita and T. Kato, On the Navier-Stokes initial value problem I, Archive for Rational Mechanics and Analysis, 16 (1964), 269-315.
- <span id="page-31-6"></span>[FMTT] C. Foias, O. T. Manley, R. Temam and Y. M. Treve, Asymptotic analysis of the Navier-Stokes equations, Physica D 9 (1983), 157-188.

- <span id="page-32-1"></span>[FRT] G. Furioli, P. Lemari´e-Rieusset, and E. Terraneo. Unicit´e dans L 3 (R 3 ) et d'autres espaces fonctionnels limites pour Navier-Stokes, Rev. Mat. Iberoamericana, 16(3) (2000), 605-667.
- <span id="page-32-3"></span>[FT] C. Foias et R. Temam, Remaques sur les ´equations de Navier-Stokes stationnaires et les ph´enom`enes successifs de bifurcation, Ann. Scuola Norm Sup Pisa 5 (1978), 29-63.
- [G] Y. Giga, Solutions for semilinear parabolic equations in L<sup>p</sup> and regularity of weak solutions of the Navier-Stokes system, J. Diff. Eq. 62 (1986), 186-212.
- [GA] G. P. Galdi, An introduction to the mathematical theory of the Navier-Stokes equations, 2nd Edition, Vol. II, Springer Tracts in Natural Philosophy, Vol. 39 Springer, New York, 1997.
- <span id="page-32-2"></span>[GIP] I. Gallagher, D. Iftimie and F. Planchon, Asymptotics and stability for global solutions to the Navier-Stokes equations, Annales de l'Institut Fourier 53 (2003), 1387- 1424.
- <span id="page-32-6"></span>[GR] P. Germain, Strong solutions and weak-strong uniqueness for the nonhomogeneous Navier-Stokes equation, J. Anal. Math. 105 (2008), 169-196.
- <span id="page-32-8"></span>[GRA] L. Grafakos, Classical and Modern Fourier Analysis, Person Prentice-Hall, New Jersey, 2004.
- <span id="page-32-7"></span>[GZ1] T. L. Gill and W. W. Zachary, Banach Spaces for the Feynman integral, Real Analysis Exchange 34(2) (2008)/(2009), 267-310.
- <span id="page-32-4"></span>[GZ2] T. L. Gill and W. W. Zachary, Feynman operator calculus: The constructive theory, Expositiones Mathematicae 29 (2011), 165-203.
- <span id="page-32-9"></span>[GZ3] T. L. Gill and W. W. Zachary, A sufficiency class for global (in time) solutions to the 3D Navier-Stokes equations, Nonlinear. Anal. 73 (2010), 3116-3122.
- <span id="page-32-0"></span>[Ho] E. Hopf, Uber die Anfangswertaufgabe f¨ur die hydrodynamischen Gru ¨ ndgleichungen, Math. Nachr. 4 (1951), 213-231.
- <span id="page-32-5"></span>[HS] R. Henstock, The General Theory of Integration, Clarendon Press, Oxford, (1991).

- <span id="page-33-8"></span>[JO] F. Jones , Lebesgue Integration on Euclidean Space, Revised Edition, Jones and Bartlett Publishers, Boston (2001).
- <span id="page-33-10"></span>[KAA] V. Koltchinskii, C.T. Abdallah, M. Ariola, P. Dorato and D. Panchenko, Improved Sample Complexity Estimates for Statistical Learning Control of Uncertain Systems, IEEE Trans. Automatic Control 45 (2000), 2383-2388.
- <span id="page-33-4"></span>[KA1] T. Kato Strong L p solutions of the Navier-Stokes equations in R <sup>m</sup> with applications to weak solutions, Math. Zeit. 187 (1984), 471-480.
- <span id="page-33-7"></span>[KA2] T. Kato Nonlinear semigroups and evolution equations, J. Math. Soc. Japan 19 (1967), 508-520.
- <span id="page-33-2"></span>[KF] T. Kato and H. Fujita, On the nonstationary Navier-Stokes system, Rend. Sem. Mat. Univ. Padova 32 (1966), 243-260.
- <span id="page-33-5"></span>[KT] H. Koch and D. Tataru. Well-posedness for the Navier-Stokes equations, Adv. Math. 157(1) (2001), 22-35.
- <span id="page-33-9"></span>[KW] J. Kurzweil, Nichtabsolut konvergente Integrale, Teubner-Texte z¨ur Mathematik, Band 26, Teubner Verlagsgesellschaft, Leipzig, (1980).
- <span id="page-33-3"></span>[LE] P. G. Lemarie-Rieuss´et, Recent developments in the Navier-Stokes problem, CHAP-MAN & HALL/CRC Press, Boca Raton, Florida 2002.
- <span id="page-33-0"></span>[Le] J. Leray, Sur le mouvement d'un liquide visqueux emplissant l'espace, Acta. Math. 63 (1934), 193–248.
- <span id="page-33-1"></span>[Li1] J. L. Lions, Quelques m´ethodes de resolution des proble`mes aux limites non lin´eaires, Dunod, Paris, 1969.
- <span id="page-33-11"></span>[Li2] J. L. Lions, On some problems connected with Navier-Stokes equations, in "Nonlinear Evolution Equations" (Proc. Sympos. at the University of Wisconsin-Madison,1970) M.C. Crandall, ed. Academic Press, New York 1978.
- <span id="page-33-6"></span>[M] I. Miyadera, Nonlinear semigroups, Translations of Mathematical Monographs, Vol. 109, Amer. Math. Soc., Providence, RI, 1977.

- <span id="page-34-3"></span>[MI] T. Miyakawa On the space-time decay properties of nonstationary incompressible Navier-Stokes flows in R n , Funkcialaj Ekvacioj 43 (2000), 541-557.
- <span id="page-34-5"></span>[NI] H. Niederreiter, Random Number Generation and Quasi-Monte Carlo Methods, SIAM, (1992).
- <span id="page-34-2"></span>[PL] F. Planchon. Global strong solutions in Sobolev or Lebesgue spaces to the incompressible Navier Stokes equations in R 3 , Annales de l'Institut Henri Poincar´e, 13 (1996), 319-336.
- <span id="page-34-6"></span>[PT] A. Papageorgiou and J.G. Traub, Faster Evaluation of Multidimesional Integrals, Computers in Physics, Nov., (1997), 574-578.
- <span id="page-34-7"></span>[PTR] S. Paskov and J.G. Traub, Faster Valuation of Financial Derivatives, Journal of Portfolio Management, 22 (1995), 113-120.
- [PZ] A. Pazy, Semigroups of Linear Operators and Applications to Partial Differential Equations Applied Mathematical Sciences, 44, Springer New York, (1983).
- <span id="page-34-4"></span>[SY] G. R. Sell and Y. You, Dynamics of evolutionary equations, Applied Mathematical Sciences, Vol. 143, Springer, New York, 2002.
- [ST] V. Steadman, Theory of operators on Banach spaces, Ph.D thesis, Howard University, 1988.
- <span id="page-34-0"></span>[T1] R. Temam, Navier-Stokes Equations, Theory and Numerical Analysis, AMS Chelsea Pub., Providence, RI, 2001.
- [T2] R. Temam, Infinite dimensional dynamical systems in mechanics and physics, Applied Mathematical Sciences, Vol. 68, Springer, New York, 1988.
- <span id="page-34-1"></span>[T3] R. Temam, Some developments on the Navier-Stokes equations in the second half of the 20th century, in D´evelopment des Math´ematiques au cours de la seconde half du XX´eme si´ecle, J. P. Pier, Managing Editor, Birkhauser, Boston, 1999.
- [TA] E. Talvila, The distributional Denjoy integral, Real Analysis Exchange, 33 (2008), 51-82.

- <span id="page-35-2"></span>[TH] P. D. Thompson Some exact statistics of two-dimensional viscous flow with random forcing, J. Fluid Mech. 55 (1972), 711-717.
- <span id="page-35-0"></span>[vW] W. von Wahl, The equations of Navier-Stokes and abstract parabolic equations, F. Vieweg und Sohn, Braunschweig, 1985.
- <span id="page-35-1"></span>[WE] F. Weissler, The Navier-Stokes Initial Value Problem in L p , Archiv for Rational Mechanics and Analysis, 74 (1980), 219-230.
- <span id="page-35-3"></span>[YE] L. T. Yeong, Henstock-Kurzweil Integration on Euclidean Spaces, Series in Real Analysis, Vol. 12, World Scientific, New Jersey, 2011.

(Tepper L. Gill) Department of Mathematics, E&CE and Computational Physics Laboratory, Howard University, Washington DC 20059, USA, E-mail : tgill@howard.edu

(Daniel Williams) Department of Mathematics, Howard University, Washington DC 20059, USA, E-mail : dwilliams@howard.edu

(Woodford W. Zachary) Department of Mathematics, E&CE and Computational Physics Laboratory,, Howard University, Washington DC 20059, USA, E-mail : wwzachary@earthlink.net