# A SUFFICIENCY CLASS FOR GLOBAL (IN TIME) SOLUTIONS TO THE 3D NAVIERSTOKES EQUATIONS II

T. L. GILL AND W. W. ZACHARY

Abstract. In this paper, we simplify and extend the results of [\[GZ\]](#page-13-0) to include the case in which Ω = R 3 . Let [L 2 (R 3 )]<sup>3</sup> be the Hilbert space of square integrable functions on R<sup>3</sup> and let H[R<sup>3</sup> <sup>3</sup> =: H be the completion of the set, <sup>u</sup> <sup>∈</sup> (C∞<sup>0</sup> [R<sup>3</sup> ])3 | ∇ · u = 0 , with respect to the inner product of [L<sup>2</sup> (R<sup>3</sup> )]3 . In this paper, we consider sufficiency conditions on a class of functions in H which allow global-in-time strong solutions to the three-dimensional Navier-Stokes equations on R<sup>3</sup> . These equations describe the time evolution of the fluid velocity and pressure of an incompressible viscous homogeneous Newtonian fluid in terms of a given initial velocity and given external body forces. Our approach uses the analytic nature of the Stokes semigroup to construct an equivalent norm for H which allows us to prove a reverse of the Poincar´e inequality. This result allows us to provide strong bounds on the nonlinear term. We then prove that, under appropriate conditions, there exists a positive constant u+, depending only on the domain, the viscosity and the body forces such that, for all functions in a dense set D contained in the closed ball B(R 3 ) =: B of radius (1/2)u<sup>+</sup> in H, the Navier-Stokes equations have unique strong solutions in C 1 ((0, ∞), H).

<sup>1991</sup> Mathematics Subject Classification. Primary (35Q30) Secondary(47H20), (76DO3) .

Key words and phrases. Global (in time), 3D-Navier-Stokes Equations.

#### Introduction

Let  $[L^2(\mathbb{R}^3)]^3$  be the Hilbert space of square integrable functions on  $\mathbb{R}^3$  and let  $\mathbb{H}_0[\mathbb{R}^3]$  be the completion of the set of functions in  $\{\mathbf{u} \in \mathbb{C}_0^\infty[\mathbb{R}^3]^3 \mid \nabla \cdot \mathbf{u} = 0\}$  which vanish at infinity with respect to the inner product of  $[L^2(\mathbb{R}^3)]^3$ , and let  $\mathbb{V}_0[\mathbb{R}^3]$  be the completion of the above functions which vanish at infinity with respect to the inner product of  $\mathbb{H}_0^1[\mathbb{R}^3]$ , the functions in  $\mathbb{H}_0[\mathbb{R}^3]$  with weak derivatives in  $[L^2(\mathbb{R}^3)]^3$ . The global-in-time classical Navier-Stokes initial-value problem (on  $\mathbb{R}^3$  and all T > 0) is to find functions  $\mathbf{u} : [0,T] \times \mathbb{R}^3 \to \mathbb{R}^3$  and  $p : [0,T] \times \mathbb{R}^3 \to \mathbb{R}$  such that

$$\partial_t \mathbf{u} + (\mathbf{u} \cdot \nabla) \mathbf{u} - \nu \Delta \mathbf{u} + \nabla p = \mathbf{f}(t) \text{ in } (0, T) \times \mathbb{R}^3,$$

$$\nabla \cdot \mathbf{u} = 0 \text{ in } (0, T) \times \mathbb{R}^3 \text{ (in the weak sense)},$$

$$\lim_{\|\mathbf{x}\| \to \infty} \mathbf{u}(t, \mathbf{x}) = 0 \text{ on } (0, T) \times \mathbb{R}^3,$$

$$\mathbf{u}(0, \mathbf{x}) = \mathbf{u}_0(\mathbf{x}) \text{ in } \mathbb{R}^3.$$

The equations describe the time evolution of the fluid velocity  $\mathbf{u}(\mathbf{x},t)$  and the pressure p of an incompressible viscous homogeneous Newtonian fluid with constant viscosity coefficient  $\nu$  in terms of a given initial velocity  $\mathbf{u}_0(\mathbf{x})$  and given external body forces  $\mathbf{f}(\mathbf{x},t)$ . (Note that our third condition,  $\lim_{\|\mathbf{x}\|\to\infty}\mathbf{u}(t,\mathbf{x})=0$  on  $(0,T)\times\mathbb{R}^3$ , is natural in this case since it is well-known that  $\mathbb{H}_0^k[\mathbb{R}^3]^3=\mathbb{H}^k[\mathbb{R}^3]^3$  (see Stein [S] or [SY].)

### Purpose

Let  $\mathbb{P}$  be the (Leray) orthogonal projection of  $(L^2[\mathbb{R}^3])^3$  onto  $\mathbb{H}_0[\mathbb{R}^3]$  and define the Stokes operator by:  $\mathbf{A}\mathbf{u} =: -\mathbb{P}\Delta\mathbf{u}$ , for  $\mathbf{u} \in D(\mathbf{A}) \subset \mathbb{H}_0^2[\mathbb{R}^3]$ , the domain of  $\mathbf{A}$ . The purpose of this paper is to prove that there exists a number  $u_+$ , depending A SUFFICIENCY CLASS FOR GLOBAL (IN TIME) SOLUTIONS TO THE 3D NAVIERSTOKES EQUATIONS IS

only on  $\mathbf{A}$ , f and  $\nu$  such that, for all functions in a certain subset (defined in the paper) of  $\mathbb{D} = D(\mathbf{A}) \cap \mathbb{B}$ , where  $\mathbb{B}$  is the closed ball of radius  $\frac{1}{2}u_+$  in  $\mathbb{H}_0(\mathbb{R}^3)$ , the Navier-Stokes equations have unique strong solutions in  $\mathbf{u} \in L^{\infty}_{loc}[[0, \infty); \mathbb{V}_0(\mathbb{R}^3)] \cap \mathbb{C}^1[(0, \infty); \mathbb{H}_0(\mathbb{R}^3)]$ .

#### Preliminaries

Applying the Leray projection to equation (1), with  $C(\mathbf{u}, \mathbf{u}) = \mathbb{P}(\mathbf{u} \cdot \nabla)\mathbf{u}$ , we can recast equation (1) in the standard form:

$$\partial_t \mathbf{u} = -\nu \mathbf{A} \mathbf{u} - C(\mathbf{u}, \mathbf{u}) + \mathbb{P} \mathbf{f}(t) \text{ in } (0, T) \times \mathbb{R}^3,$$

(2) 
$$\lim_{\|\mathbf{x}\| \to \infty} \mathbf{u}(t, \mathbf{x}) = 0 \text{ on } (0, T) \times \mathbb{R}^3,$$

$$\mathbf{u}(0,\mathbf{x}) = \mathbf{u}_0(\mathbf{x}) \text{ in } \mathbb{R}^3,$$

where we have used the fact that the orthogonal complement of  $\mathbb{H}_0$  relative to  $\{L^2(\mathbb{R}^3)\}^3$  is  $\{\mathbf{v}: \mathbf{v} = \nabla q, q \in (H^1)^3\}$  to eliminate the pressure term (see Galdi [GA] or [SY, T1,T2]).

**Definition 1.** We say that the operator  $A(\cdot,t)$  is (for each t)

- (1) 0-Dissipative if  $\langle \mathcal{A}(\mathbf{u},t), \mathbf{u} \rangle_{\mathbb{H}} \leq 0$ .
- (2) Dissipative if  $\langle \mathcal{A}(\mathbf{u},t) \mathcal{A}(\mathbf{v},t), \mathbf{u} \mathbf{v} \rangle_{\mathbb{H}} \leq 0$ .
- (3) Strongly dissipative if there exists an  $\delta > 0$  such that

$$\langle \mathcal{A}(\mathbf{u}, t) - \mathcal{A}(\mathbf{v}, t), \mathbf{u} - \mathbf{v} \rangle_{\mathbb{H}} \le -\delta \|\mathbf{u} - \mathbf{v}\|_{\mathbb{H}}^{2}.$$

Note that, if  $\mathcal{A}(\cdot, t)$  is a linear operator, definitions (1) and (2) coincide. Theorem 2 below is essentially due to Browder [B], see Zeidler [Z, Corollary 32.27, page 868 and Corollary 32.35 page 887, in Vol. IIB], while Theorem 3 is from Miyadera [M,

p. 185, Theorem 6.20], and is a modification of the Crandall-Liggett Theorem [\[CL\]](#page-12-1) (see the appendix to the first section of [\[CL\]](#page-12-1)) .

Theorem 2. *Let* B *be a closed, bounded, convex subset of* H*. If* A(·, t) : B → H *is a strongly dissipative mapping for each fixed* t ≥ 0*, then for each* b ∈ B*, there is a* u ∈ B *with* A(u, t) = b *(i.e., the range,* Ran*[*A(·, t)] ⊃ B*).*

Theorem 3. *Let* A(·, t), t ∈ I = [0, ∞)} *be a family of operators defined on* H *with domains* D(A(·, t)) = D*, independent of* t*. We assume that* D = D(A) ∩ B *is a closed convex set (in an appropriate topology):*

- (1) *The operator* A(·, t) *is the generator of a contraction semigroup for each* t ∈ I*.*
- (2) *The function* A(u, t) *is continuous in both variables on* D × I*.*

*Then, for every* u<sup>0</sup> ∈ D*, the problem* ∂tu(t, x) = A(u(t, x), t)*,* u(0, x) = u0(x)*, has a unique solution* u(t, x) ∈ C 1 (I; D)*.*

Stokes Equation. The difficulty in proving the existence of global-in-time strong solutions for equation (2) is directly linked to the problem of getting good estimates for the nonlinear term C(u, u). In [\[GZ\]](#page-13-0), we obtained an extension of the important result due to Constantin and Foias [\[CF\]](#page-12-2)). This result, see below, is one of the major estimates used to study this equation. In what follows, we assume that u, v ∈ D(A).

Theorem 4. *Let* 0 ≤ α<sup>i</sup> , 1 ≤ i ≤ 3*, satisfy* α<sup>1</sup> + α<sup>2</sup> + α<sup>3</sup> = 3/2 *and*

$$(\alpha_1, \alpha_2, \alpha_3) \notin \{(3/2, 0, 0), (0, 3/2, 0), (0, 0, 3/2)\}.$$

A SUFFICIENCY CLASS FOR GLOBAL (IN TIME) SOLUTIONS TO THE 3D NAVIERSTOKES EQUATIONS IN Then there is a positive constant  $c=c(\alpha_i,\mathbb{R}^3)$  such that

(3) 
$$|\langle C(\mathbf{u}, \mathbf{v}), \mathbf{w} \rangle_{\mathbb{H}}| \le c \|\mathbf{A}^{\alpha_1/2} \mathbf{u}\|_{\mathbb{H}} \|\mathbf{A}^{(1+\alpha_2)/2} \mathbf{v}\|_{\mathbb{H}} \|\mathbf{A}^{\alpha_3/2} \mathbf{w}\|_{\mathbb{H}}.$$

In [GZ] we showed that, by renorming  $\mathbb{H}$ , we could prove a very strong inequality for equation (3). Since this result is not well-known and important for this paper, we give a proof. First we need to review the Stokes equation.

If we drop the nonlinear term, we get the well-known Stokes equation ( $\mathbb{P}\mathbf{f}(t) = \mathbf{0}$ ):

$$\partial_t \mathbf{u} = -\nu \mathbf{A} \mathbf{u} \text{ in } (0, T) \times \mathbb{R}^3,$$

$$\lim_{\|\mathbf{x}\| \to \infty} \mathbf{u}(t, \mathbf{x}) = 0 \text{ on } (0, T) \times \mathbb{R}^3,$$

$$\mathbf{u}(0, \mathbf{x}) = \mathbf{u}_0(\mathbf{x}) \text{ in } \mathbb{R}^3.$$

A proof of the next theorem may be found in Sell and You [SY] (page 114):

**Theorem 5.** Let **A** be the Stokes operator on  $\mathbb{R}^3$ . Then the following holds:

- (1) The operator  $\mathbf{A}$  is a positive selfadjoint generator of a contraction semigroup T(t).
- (2) The operator **A** is sectorial and T(t) is analytic.

**Equivalent Norms.** Recall that an equivalent norm,  $\|\cdot\|_{\mathcal{H},1}$ , on a Hilbert space  $\mathcal{H}$ , with norm  $\|\cdot\|_{\mathcal{H}}$ , is one that satisfies: (for positive constants  $M, M_1$ )

$$||u||_{\mathcal{H}} \leqslant M ||u||_{\mathcal{H}_1} \leqslant M_1 ||u||_{\mathcal{H}}, \ u \in \mathcal{H}.$$

It is easy to show that any equivalent norm on  $\mathcal{H}$  can be identified with a transformation of  $\mathcal{H}$  which preserves the topology. In order to see how an equivalent norm can help us, let  $T(t) = exp\{-t\mathbf{A}\}$  be the analytic contraction semigroup generated by the Stokes operator  $\mathbf{A}$ , with  $\|T(t)\mathbf{u}\|_{\mathbb{H}} \leqslant e^{-\omega t} \|\mathbf{u}\|_{\mathbb{H}}$ . Let  $S(t) = e^{\omega t}T(t)$  and

choose M so that  $\|\mathbf{u}\|_{\mathbb{H},1} = \|S(r)\mathbf{u}\|_{\mathbb{H}}$  is an equivalent norm, where r is a fixed value, to be determined. Since  $\mathbf{A}$  is analytic, there is a constant  $c_z$  such that, for  $\mathbf{u} \in D(\mathbf{A}^z)$ ,

$$\|\mathbf{A}^z \mathbf{u}\|_{\mathbb{H},1} = e^{\omega r} \|\mathbf{A}^z T(r) \mathbf{u}\|_{\mathbb{H}} \leqslant e^{\omega r} e^{-\omega r} \frac{c_z}{(r)^z} \|\mathbf{u}\|_{\mathbb{H}} \leqslant \frac{Mc_z}{(r)^z} \|\mathbf{u}\|_{\mathbb{H},1}.$$

Since the norms are equivalent, we also have that

(4) 
$$\|\mathbf{A}^{z}\mathbf{u}\|_{\mathbb{H}} \leq M \|\mathbf{A}^{z}\mathbf{u}\|_{\mathbb{H},1} \leqslant \frac{M^{2}c_{z}}{(r)^{z}} \|\mathbf{u}\|_{\mathbb{H},1}.$$

From Theorem 4, we have the following result:

Theorem 6. Let  $\mathbf{u} \in D(\mathbf{A})$ , set  $\mathbf{S} = S(r)$  and renorm  $\mathbb{H}$  so that  $\|\mathbf{u}\|_{\mathbb{H},1} = \|\mathbf{S}\mathbf{u}\|_{\mathbb{H}}$ . We define  $\mathbf{b}(\mathbf{u}, \mathbf{v}, \mathbf{w})_{\mathbb{H},1} = \langle \mathbf{S}C(\mathbf{u}, \mathbf{v}), \mathbf{S}\mathbf{w} \rangle_{\mathbb{H}}$ . Then:

(1) If we let  $\alpha_1 = 0$ ,  $\alpha_2 = 1$  and  $\alpha_3 = 1/2$ , there are positive constants  $c = c(\alpha_i, \mathbb{R}^3)$ ,  $c_1$  and  $c_2$  such that

(5) 
$$\left| \langle C(\mathbf{u}, \mathbf{v}), \mathbf{w} \rangle_{\mathbb{H}, 1} \right| \leq \frac{M^4 c c_1 c_2}{r^{5/4}} \|\mathbf{u}\|_{\mathbb{H}, 1} \|\mathbf{w}\|_{\mathbb{H}, 1} \|\mathbf{v}\|_{\mathbb{H}, 1} \text{ and} \right| \\ \left| \langle C(\mathbf{v}, \mathbf{u}), \mathbf{w} \rangle_{\mathbb{H}, 1} \right| \leq \frac{M^4 c c_1 c_2}{r^{5/4}} \|\mathbf{u}\|_{\mathbb{H}, 1} \|\mathbf{w}\|_{\mathbb{H}, 1} \|\mathbf{v}\|_{\mathbb{H}, 1}.$$

$$(2)$$

(6) 
$$\max\{\|C(\mathbf{u}, \mathbf{v})\|_{\mathbb{H}, 1}, \|C(\mathbf{v}, \mathbf{u})\|_{\mathbb{H}, 1}\} \leqslant \frac{M^4 c c_1 c_2}{r^{5/4}} \|\mathbf{u}\|_{\mathbb{H}, 1} \|\mathbf{v}\|_{\mathbb{H}, 1}.$$

*Proof.* We prove the first equation of (5), the proof of the second is similar. Set  $S(r) = \mathbf{S}$  and  $\mathbf{S}^2 \mathbf{w} = \mathbf{w}_1$ , then we have:

$$\mathbf{b}(\mathbf{u}, \mathbf{v}, \mathbf{w})_{\mathbb{H}, 1} = \langle \mathbf{S}C(\mathbf{u}, \mathbf{v}), \mathbf{S}\mathbf{w} \rangle_{\mathbb{H}} = \mathbf{b}(\mathbf{u}, \mathbf{v}, \mathbf{w}_1)_{\mathbb{H}}.$$

Using the selfadjoint property of **A**, and integration by parts, we have

$$\mathbf{b}(\mathbf{u}, \mathbf{v}, \mathbf{w}_1)_{\mathbb{H}} = -\mathbf{b}(\mathbf{u}, \mathbf{w_1}, \mathbf{v})_{\mathbb{H}}.$$

A SUFFICIENCY CLASS FOR GLOBAL (IN TIME) SOLUTIONS TO THE 3D NAVIERSTOKES EQUATIONS IT It follows that:

$$\left| \langle C(\mathbf{u}, \mathbf{v}), \mathbf{w} \rangle_{\mathbb{H}, 1} \right| \leq c \left\| \mathbf{A}^{\alpha_1/2} \mathbf{u} \right\|_{\mathbb{H}} \left\| \mathbf{A}^{(1+\alpha_2)/2} \mathbf{w}_1 \right\|_{\mathbb{H}} \left\| \mathbf{A}^{\alpha_3/2} \mathbf{v} \right\|_{\mathbb{H}}.$$

Setting  $\alpha_1 = 0$ ,  $\alpha_2 = 1$   $\alpha_3 = 1/2$  and, using equation (4), we have:

$$\begin{split} & \left| \langle C(\mathbf{u}, \mathbf{v}), \mathbf{w} \rangle_{\mathbb{H}, 1} \right| \leq c \, \|\mathbf{u}\|_{\mathbb{H}} \, \|\mathbf{A}\mathbf{w}_{1}\|_{\mathbb{H}} \, \left\| \mathbf{A}^{1/4} \mathbf{v} \right\|_{\mathbb{H}} \\ & \leq \frac{M c c_{1} c_{2}}{r^{5/4}} \, \|\mathbf{u}\|_{\mathbb{H}} \, \|\mathbf{w}\|_{\mathbb{H}} \, \|\mathbf{v}\|_{\mathbb{H}} \\ & \leq \frac{M^{4} c c_{1} c_{2}}{r^{5/4}} \, \|\mathbf{u}\|_{\mathbb{H}, 1} \, \|\mathbf{w}\|_{\mathbb{H}, 1} \, \|\mathbf{v}\|_{\mathbb{H}, 1} \, . \end{split}$$

The proof of (6) is clear.

The following extension of the Poincaré inequality will also prove useful.

**Lemma 7.** Let  $\mathbf{A}^{1/2}$  generate an analytic contraction semigroup S(t). If r > 0 is any fixed number, then there exists an  $\alpha = \alpha(r) > 0$  such that for any  $\mathbf{u} \in D(\mathbf{A}^{1/2})$ ,

$$\alpha^{1/2} \|\mathbf{u}\|_{H,1} \leqslant r^{1/2} \|\mathbf{A}^{1/2}\mathbf{u}\|_{H,1}.$$

*Proof.* First observe that

$$\int_0^{r^{1/2}} \mathbf{A}^{1/2} S(t) \mathbf{u} dt = \int_0^{r^{1/2}} \frac{d}{dt} S(t) \mathbf{u} dt = S(r^{1/2}) \mathbf{u} - \mathbf{u}.$$

Now choose  $\alpha^{1/2}$  so that  $\|S(r^{1/2})\mathbf{u} - \mathbf{u}\|_{H,1} \ge \alpha^{1/2} \|\mathbf{u}\|_{H,1}$ . It follows that

$$\alpha^{1/2} \|\mathbf{u}\|_{H,1} \leqslant \int_{0}^{r^{1/2}} \left\| S(t) \mathbf{A}^{1/2} \mathbf{u} \right\|_{H,1} dt \leqslant r^{1/2} \left\| \mathbf{A}^{1/2} \mathbf{u} \right\|_{H,1}.$$

## M-Dissipative Conditions

Let us assume that f(t) ∈ L<sup>∞</sup>[[0, ∞); H] and is H¨older continuous in t, with kf(t) − f(τ)kH,<sup>1</sup> ≤ d |t − τ | θ , d > 0, 0 < θ < 1. We can now rewrite equation (2) in the form:

(7) 
$$\partial_t \mathbf{u} = \mathcal{A}(\mathbf{u}, t) \text{ in } (0, T) \times \mathbb{R}^3,$$
$$\mathcal{A}(\mathbf{u}, t) = -\nu \mathbf{A} \mathbf{u} - C(\mathbf{u}, \mathbf{u}) + \mathbb{P} \mathbf{f}(t).$$

We begin with a study of the operator A(·, t), for fixed t, and seek conditions depending on A, ν, and f(t) which guarantee that A(·, t) is m-dissipative for each t. Clearly A(·, t) is defined on D(A) and, since νA is a closed positive (m-accretive) operator, −ν(A) generates a linear contraction semigroup. Thus, we need to ensure that νA(·, t) will be m-dissipative for each t. The following Lemma follows from the properties of f(t).

Lemma 8. *For* t ∈ I = [0, ∞) *and, for each fixed* u ∈ D(A)*,* A(u, t) *is H¨older continuous, with* kA(u, t) − A(u, τ)k<sup>H</sup>,<sup>1</sup> ≤ d |t − τ| θ *, where* d *is the H¨older constant for the function* f(t)*.*

### Main Results

Theorem 9. *Let* <sup>f</sup> = supt∈R<sup>+</sup> <sup>k</sup>Pf(t)k<sup>H</sup>,<sup>1</sup> <sup>&</sup>lt; <sup>∞</sup>*, then there exists a positive constants* u+, u−*, depending only on* f*,* A *and* ν *such that, for all* u *with* 0 ≤ u<sup>−</sup> ≤ kuk<sup>H</sup>,<sup>1</sup> ≤ u+*,* A(·, t) *is strongly dissipative.*

*Proof.* The proof of our first assertion has two parts. First, we require that the nonlinear operator A(·, t) be 0-dissipative, which gives us an upper bound u<sup>+</sup> and lower bound u<sup>−</sup> in terms of the norm (i.e., kukH,<sup>1</sup> 6 u<sup>+</sup> ). We then use this part

#### A SUFFICIENCY CLASS FOR GLOBAL (IN TIME) SOLUTIONS TO THE 3D NAVIERSTOKES EQUATIONS 19

to show that  $\mathcal{A}(\cdot,t)$  is strongly dissipative on any closed convex ball,  $\mathbb{D}$  inside the annulus defined by  $\left\{\mathbf{u} \in D(\mathbf{A}) : 0 \le u_- \le \|\mathbf{u}\|_{\mathbb{H},1} \le \frac{1}{2}u_+\right\}$ .

Part 1) From equation (7), we consider the expression

$$\begin{split} & \left\langle \mathcal{A}(\mathbf{u},t),\mathbf{u} \right\rangle_{\mathbb{H},1} = -\nu \left\langle \mathbf{A}\mathbf{u},\mathbf{u} \right\rangle_{\mathbb{H},1} + \left\langle \left[ -C(\mathbf{u},\mathbf{u}) + \mathbb{P}\mathbf{f} \right],\mathbf{u} \right\rangle_{\mathbb{H},1} \\ & = -\nu \left\| \mathbf{A}^{1/2}\mathbf{u} \right\|_{\mathbb{H}_{1}}^{2} - \left\langle C(\mathbf{u},\mathbf{u}),\mathbf{u} \right\rangle_{\mathbb{H},1} + \left\langle \mathbb{P}\mathbf{f},\mathbf{u} \right\rangle_{\mathbb{H},1}. \end{split}$$

We now use the fact that

$$\frac{\alpha}{r} \left\| \mathbf{u} \right\|_{\mathbb{H},1}^{2} \leqslant \left\| \mathbf{A}^{1/2} \mathbf{u} \right\|_{\mathbb{H},1}^{2} \Rightarrow -\nu \frac{\alpha}{r} \left\| \mathbf{u} \right\|_{\mathbb{H},1}^{2} \geqslant -\nu \left\| \mathbf{A}^{1/2} \mathbf{u} \right\|_{\mathbb{H},1}^{2}$$

to get that

(8) 
$$\langle \mathcal{A}(\mathbf{u}, t), \mathbf{u} \rangle_{\mathbb{H}, 1} \leq -\nu \frac{\alpha}{r} \|\mathbf{u}\|_{\mathbb{H}, 1}^2 + \frac{M^4 c c_1 c_2}{r^{5/4}} \|\mathbf{u}\|_{\mathbb{H}, 1}^3 + f \|\mathbf{u}\|_{\mathbb{H}, 1}.$$

In the last line, we used our estimate from Theorem 6.

Since  $\|\mathbf{u}\|_{\mathbb{H},1} > 0$ , we have that  $\mathcal{A}(\cdot,t)$  is 0-dissipative if

$$-\frac{\nu\alpha}{r} \|\mathbf{u}\|_{H,1} + \frac{M^4cc_1c_2}{r^{5/4}} \|\mathbf{u}\|_{H,1}^2 + f \leqslant 0$$

Solving the equality, we get that

$$u_{\pm} = \frac{\nu \alpha r^{1/4}}{2M^4 c c_1 c_2} \left\{ 1 \pm \sqrt{1 - (4f M^4 c c_1 c_2) / (\nu^2 r^{1/2} \alpha^2)} \right\} = \frac{\nu \alpha r^{1/4}}{2M^4 c c_1 c_2} \left\{ 1 \pm \sqrt{1 - \gamma} \right\},$$

where  $\gamma = (4r^{3/4}fM^4cc_1c_2)/(\nu^2\alpha^2)$ . Since we want real distinct solutions, we must require that

$$\gamma = \frac{4r^{3/4}fM^4cc_1c_2}{\nu^2\alpha^2} < 1 \quad \Rightarrow \frac{2M^2r^{3/8}}{\alpha} \left[fcc_1c_2\right]^{1/2} < \nu.$$

It follows that, if  $\mathbb{P}\mathbf{f} \neq \mathbf{0}$ , then  $u_- < u_+$ , and our requirement that  $\mathcal{A}(\mathbf{u}, t)$  is 0-dissipative implies that, since our solution factors as  $(\|\mathbf{u}\|_{\mathbb{H},1} - u_+)(\|\mathbf{u}\|_{\mathbb{H},1} - u_-) \leq 0$ ,

we must have that:

$$\|\mathbf{u}\|_{\mathbb{H}_1} - u_+ \le 0, \|\mathbf{u}\|_{\mathbb{H}_1} - u_- \ge 0.$$

It follows that, for  $u_- \leq \|\mathbf{u}\|_{\mathbb{H},1} \leq u_+$ ,  $\langle \mathcal{A}(\mathbf{u},t), \mathbf{u} \rangle_{\mathbb{H},1} \leq 0$ . (It is clear that, when  $\mathbb{P}\mathbf{f}(t) = \mathbf{0}$ ,  $u_- = 0$ , and  $u_+ = \frac{\nu \alpha r^{1/4}}{M^4 c c_1 c_2}$ .)

Part 2): Now, for any  $\mathbf{u}, \mathbf{v} \in D(\mathbf{A})$  with  $\mathbf{u} - \mathbf{v} \in D(\mathbf{A})$  and

$$u_{-} \le \min(\|\mathbf{u}\|_{\mathbb{H},1}, \|\mathbf{v}\|_{\mathbb{H},1}) \le \max(\|\mathbf{u}\|_{\mathbb{H},1}, \|\mathbf{v}\|_{\mathbb{H},1}) \le (1/2)u_{+},$$

we have that

$$\begin{split} & \langle \mathcal{A}(\mathbf{u},t) - \mathcal{A}(\mathbf{v},t), (\mathbf{u} - \mathbf{v}) \rangle_{\mathbb{H},1} = -\nu \left\| \mathbf{A}^{1/2} (\mathbf{u} - \mathbf{v}) \right\|_{\mathbb{H},1}^{2} \\ & - \langle [C(\mathbf{u}, \mathbf{u} - \mathbf{v}) + C(\mathbf{v}, \mathbf{u} - \mathbf{v})], (\mathbf{u} - \mathbf{v}) \rangle_{\mathbb{H},1} \\ & \leqslant -\frac{\nu\alpha}{r} \left\| \mathbf{u} - \mathbf{v} \right\|_{\mathbb{H},1}^{2} + [1/(r^{5/4})] M^{4} c c_{1} c_{2} \left\| \mathbf{u} - \mathbf{v} \right\|_{\mathbb{H},1}^{2} \left( \left\| \mathbf{u} \right\|_{\mathbb{H},1} + \left\| \mathbf{v} \right\|_{\mathbb{H},1} \right) \\ & \leq -\frac{\nu\alpha}{r} \left\| \mathbf{u} - \mathbf{v} \right\|_{\mathbb{H},1}^{2} + [1/(r^{5/4})] M^{4} c c_{1} c_{2} \left\| \mathbf{u} - \mathbf{v} \right\|_{\mathbb{H},1}^{2} u_{+} \\ & = -\frac{\nu\alpha}{r} \left\| \mathbf{u} - \mathbf{v} \right\|_{\mathbb{H},1}^{2} + [1/(r^{5/4})] M^{4} c c_{1} c_{2} \left\| \mathbf{u} - \mathbf{v} \right\|_{\mathbb{H},1}^{2} \left( \frac{\nu\alpha r^{1/4}}{2M^{4} c c_{1} c_{2}} \left\{ 1 + \sqrt{1 - \gamma} \right\} \right) \\ & = -\frac{\nu\alpha}{2r} \left\| \mathbf{u} - \mathbf{v} \right\|_{\mathbb{H},1}^{2} \left\{ 1 - \sqrt{1 - \gamma} \right\} . \end{split}$$

Let  $\mathbb D$  be any closed convex set (in the graph norm of  $\mathbf A$ ) inside the annulus bounded by  $\frac{1}{2}u_+$  and  $u_-$ .

**Theorem 10.** The operator  $\mathcal{A}(\cdot,t)$  is closed, strongly dissipative and jointly continuous in  $\mathbf{u}$  and t. Furthermore, for each  $t \in \mathbf{R}^+$  and  $\beta > 0$ ,  $Ran[I - \beta \mathcal{A}(t)] \supset \mathbb{D}$ , so that  $\mathcal{A}(t)$  is m-dissipative on  $\mathbb{D}$ .

### A SUFFICIENCY CLASS FOR GLOBAL (IN TIME) SOLUTIONS TO THE 3D NAVIERSTOKES EQUATIONS II 11

*Proof.* It is easy to see that A(·, t) is closed. Since A(·, t) is strongly dissipative, it is maximal dissipative, so that Ran[I − βA(·, t)] ⊃ D. It follows that A(·, t) is m-dissipative on D for each t ∈ R<sup>+</sup> (since H is a Hilbert space). To see that A(u, t) is continuous in both variables, let <sup>u</sup>n, <sup>u</sup> <sup>∈</sup> <sup>B</sup>, <sup>k</sup>(u<sup>n</sup> <sup>−</sup> <sup>u</sup>)kH,<sup>1</sup> <sup>→</sup> 0, with <sup>t</sup>n, t <sup>∈</sup> <sup>I</sup> and t<sup>n</sup> → t. Then, if kAukH,<sup>1</sup> 6 Mc<sup>3</sup> r kukH,<sup>1</sup> , we have

$$\begin{split} &\|\mathcal{A}(\mathbf{u}_{n},t_{n}) - \mathcal{A}(\mathbf{u},t)\|_{\mathbb{H},1} \leq \|\mathcal{A}(\mathbf{u},t_{n}) - \mathcal{A}(\mathbf{u},t)\|_{\mathbb{H},1} + \|\mathcal{A}(\mathbf{u}_{n},t_{n}) - \mathcal{A}(\mathbf{u},t_{n})\|_{\mathbb{H},1} \\ &= \|[\mathbb{P}\mathbf{f}(t_{n}) - \mathbb{P}\mathbf{f}(t)]\|_{\mathbb{H},1} + \|\nu\mathbf{A}(\mathbf{u}_{n} - \mathbf{u}) + [C(\mathbf{u}_{n} - \mathbf{u}, \mathbf{u}_{n}) + C(\mathbf{u}, \mathbf{u}_{n} - \mathbf{u})]\|_{\mathbb{H},1} \\ &\leq d \|t_{n} - t\|^{\theta} + \nu \|\mathbf{A}(\mathbf{u}_{n} - \mathbf{u})\|_{\mathbb{H},1} + \|C(\mathbf{u}_{n} - \mathbf{u}, \mathbf{u}_{n}) + C(\mathbf{u}, \mathbf{u}_{n} - \mathbf{u})\|_{\mathbb{H},1} \\ &\leq d \|t_{n} - t\|^{\theta} + \nu \frac{Mc_{3}}{r} \|(\mathbf{u}_{n} - \mathbf{u})\|_{\mathbb{H},1} + \frac{M^{4}cc_{1}c_{2}}{r^{5/4}} \|(\mathbf{u}_{n} - \mathbf{u})\|_{\mathbb{H},1} \left\{ \|\mathbf{u}_{n}\|_{\mathbb{H},1} + \|\mathbf{u}\|_{\mathbb{H},1} \right\} \\ &\leq d \|t_{n} - t\|^{\theta} + \nu \frac{Mc_{3}}{r} \|(\mathbf{u}_{n} - \mathbf{u})\|_{\mathbb{H},1} + \frac{M^{4}cc_{1}c_{2}}{r^{5/4}} \|(\mathbf{u}_{n} - \mathbf{u})\|_{\mathbb{H},1} u_{+}. \end{split}$$

When f = 0, D is the graph closure of D(A) ∩ B<sup>+</sup> in the H norm, where B<sup>+</sup> is the ball of radius <sup>1</sup> 2 u+. In this case, it follows that D is a closed, bounded, convex set. We now have:

It follows that A(u, t) is continuous in both variables.

Theorem 11. *For each* T ∈ R<sup>+</sup>*,* t ∈ (0, T ) *and* u<sup>0</sup> ∈ D*, the global-in-time Navier-Stokes initial-value problem in* R 3 :

$$\partial_{t}\mathbf{u} + (\mathbf{u} \cdot \nabla)\mathbf{u} - \nu \Delta \mathbf{u} + \nabla p = \mathbf{0} \ in \ (0, T) \times \mathbb{R}^{3},$$

$$\nabla \cdot \mathbf{u} = 0 \ in \ (0, T) \times \mathbb{R}^{3},$$

$$\lim_{\|\mathbf{x}\| \to \infty} \mathbf{u}(t, \mathbf{x}) = \mathbf{0} \ on \ (0, T) \times \mathbb{R}^{3},$$

$$\mathbf{u}(0, \mathbf{x}) = \mathbf{u}_{0}(\mathbf{x}) \ in \ \mathbb{R}^{3},$$

*has a unique strong solution* u(t, x)*, which is in* L 2 loc[[0, ∞); H] *and in* L<sup>∞</sup> loc[[0, ∞); V] ∩ C 1 [(0, ∞); H]*.*

*Proof.* Theorem 3 allows us to conclude that, when  $\mathbf{u}_0 \in \mathbb{D}$ , the initial value problem is solved and the solution  $\mathbf{u}(t, \mathbf{x})$  is in  $\mathbb{C}^1[(0, \infty); \mathbb{D}]$ . Since  $\mathbb{D} \subset \mathbb{H}^2$ , it follows that  $\mathbf{u}(t, \mathbf{x})$  is also in  $\mathbb{V}$  for each t > 0. It is now clear that, for any T > 0,

$$\int_0^T \left\| \mathbf{u}(t,\mathbf{x}) \right\|_{\mathbb{H}}^2 dt < \infty, \text{ and } \sup_{0 < t < T} \left\| \mathbf{u}(t,\mathbf{x}) \right\|_{\mathbb{V}}^2 < \infty.$$

This gives our conclusion.

When  $\mathbf{f} \neq \mathbf{0}$ ,  $u_{-} \neq 0$ . Let  $\mathbb{k} = \left\{ \mathbf{u} : \|\mathbf{u}\|_{\mathbb{H},1} < u_{-} \&, \|\mathbf{u}\|_{\mathbb{H},1} > \frac{1}{2}u_{+} \right\}$  and set  $\mathbb{B}_{-} = \mathbb{B} \cap \mathbb{k}^{c}$ , where  $\mathbb{k}^{c}$  is the complement of  $\mathbb{k}$ . We can now take the graph closure of  $\mathbb{B}_{-} \cap D(\mathbf{A})$  and use the largest closed convex set containing the initial data inside this set.

#### DISCUSSION

It is known that, if  $\mathbf{u}_0 \in \mathbb{V}$  and  $\mathbf{f}(t) \in L^{\infty}[(0, \infty), \mathbb{H}]$ , then there is a time T > 0 such that a weak solution with this data is uniquely determined on any subinterval of [0, T) (see Sell and You, page 396, [SY]). Thus, we also have that:

Corollary 12. For each  $t \in \mathbf{R}^+$  and  $\mathbf{u}_0 \in \mathbb{D}$  the Navier-Stokes initial-value problem on  $\mathbb{R}^3$ :

$$\partial_{t}\mathbf{u} + (\mathbf{u} \cdot \nabla)\mathbf{u} - \nu\Delta\mathbf{u} + \nabla p = \mathbf{f}(t) \ in \ (0, T) \times \mathbb{R}^{3},$$

$$\nabla \cdot \mathbf{u} = 0 \ in \ (0, T) \times \mathbb{R}^{3},$$

$$\lim_{\|\mathbf{x}\| \to \infty} \mathbf{u}(t, \mathbf{x}) = \mathbf{0} \ on \ (0, T) \times \mathbb{R}^{3},$$

$$\mathbf{u}(0, \mathbf{x}) = \mathbf{u}_{0}(\mathbf{x}) \ in \ \mathbb{R}^{3},$$

has a unique weak solution  $\mathbf{u}(t,\mathbf{x})$  which is in  $L^2_{loc}[[0,\infty);\mathbb{H}^2]$  and in  $L^\infty_{loc}[[0,\infty);\mathbb{V}]\cap \mathbb{C}^1[(0,\infty);\mathbb{H}]$ .

### A SUFFICIENCY CLASS FOR GLOBAL (IN TIME) SOLUTIONS TO THE 3D NAVIERSTOKES EQUATIONS II 13

As in [\[GZ\]](#page-13-0), our results show that the Leray-Hopf weak solutions do not develop singularities if u0(x) ∈ H<sup>2</sup> (see Giga [\[G\]](#page-12-3) and references therein).

We should note that the constant α in Lemma 7 depends on r so we can't change r without affecting α. This means that the size of u<sup>+</sup> need not increase with large values of r.

A close review of the results of this paper show that all theorems hold for the bounded domain case. This provides an improvement of the results in [\[GZ\]](#page-13-0). Furthermore, in that case, we can take α = λ1, the first eigenvalue of A, which is independent of r. This means that choosing larger values for r could increase the possible size of D for bounded domains. However, the inequality for ν must be maintained, so that increasing D is not certain.

# References

- <span id="page-12-0"></span>[B] F. Browder, Nonlinear operators and nonlinear equations of evolution in Banach spaces, Proc. Sympos. Pure Math., Vol. 18 part II, Amer. Math. Soc., Providence, RI, 1970.
- <span id="page-12-2"></span>[CF] P. Constantin and C. Foi¸as, Navier-Stokes Equations, University of Chicago Press, Chicago, IL, 1988.
- <span id="page-12-1"></span>[CL] M. Crandall and T. Liggett, Generation of semigroups of nonlinear transformations on general Banach spaces, Amer. J. Math. 93 (1971), 265-293.
- [GA] G. P. Galdi, An introduction to the mathematical theory of the Navier-Stokes equations, 2nd Edition, Vol. II, Springer Tracts in Natural Philosophy, Vol. 39, Springer, New York, 1997.
- <span id="page-12-3"></span>[G] Y. Giga, Solutions for semilinear parabolic equations in L<sup>p</sup> and regularity of weak solutions of the Navier-Stokes system, J. Diff. Eq. 62 (1986), 186-212.

- <span id="page-13-0"></span>[GZ] T. L. Gill and W. W. Zachary, Sufficiency Class for Global (in time) Solutions to The 3D-Navier-Stokes Equations, Nonlinear Analysis A: Theory, Methods & Applications (2010) DOI: 10.1016/j.na.2010.06.083.
- <span id="page-13-4"></span>[M] I. Miyadera, Nonlinear semigroups, Translations of Mathematical Monographs, Vol. 109, Amer. Math. Soc., Providence, RI, 1977.
- [PZ] A. Pazy, Semigroups of Linear Operators and Applications to Partial Differential Equations Applied Mathematical Sciences, 44, Springer New York, (1983).
- <span id="page-13-2"></span>[SY] G. R. Sell and Y. You, Dynamics of evolutionary equations, Applied Mathematical Sciences, Vol. 143, Springer, New York, 2002.
- <span id="page-13-1"></span>[S] E. M. Stein, Singular integrals and differentiability properties of functions, Princeton University Press, Princeton, NJ, 1970.
- [T1] R. Temam, Navier-Stokes Equations, Theory and Numerical Analysis, AMS Chelsea Pub., Providence, RI, 2001.
- [T2] R. Temam, Infinite dimensional dynamical systems in mechanics and physics, Applied Mathematical Sciences, Vol. 68, Springer, New York, 1988.
- <span id="page-13-3"></span>[Z] E. Zeidler, Nonlinear functional analysis and its applications, Vol. IIB, Springer, New York, 1985.
- (Tepper L. Gill) Department of Electrical Engineering, Howard University, Washington DC 20059, USA, E-mail : tgill@howard.edu

(Woodford W. Zachary) Department of Electrical Engineering, Howard University, Washington DC 20059, USA, E-mail : wwzachary@earthlink.net