#### GLOBAL (IN TIME) SOLUTIONS TO THE 3D-NAVIER-STOKES EQUATIONS ON R 3

## T. L. GILL AND W. W. ZACHARY

Abstract. In two recent papers ([\[GZ1\]](#page-15-0) [\[GZ2\]](#page-15-1)), we provided solutions to the well-known unsolved problem of constructing sufficiency classes of functions in H[R<sup>3</sup> <sup>3</sup> and V[R<sup>3</sup> 3 , which would allow global, in time, strong solutions to the three-dimensional Navier-Stokes equations. These equations describe the time evolution of the fluid velocity and pressure of an incompressible viscous homogeneous Newtonian fluid in terms of a given initial velocity and given external body forces. In both previous papers, our solution was restricted to functions defined on a bounded open domain of class C<sup>3</sup> contained in R<sup>3</sup> . In this paper, we study this problem for functions defined on all of R<sup>3</sup> . We prove that, under appropriate conditions, there exists a positive constant a and a number u+, depending only on the domain, the viscosity, the body forces and the eigenvalues of the "Hermite" Stokes operator (defined below) such that, for all functions in a dense set D contained in the closed ball B(R<sup>3</sup> ) of radius (1/2)u<sup>+</sup> in H[R<sup>3</sup> 3 , the Navier-Stokes equations have unique strong solutions in C<sup>1</sup> ` (0, ∞), H[R<sup>3</sup> 3 ´ .

## Introduction

Let L 2 [R 3 <sup>3</sup> be the real Hilbert space of square integrable functions on R 3 with values in R 3 , and let H0[R 3 <sup>3</sup> be the completion of the set of functions in

<sup>1991</sup> Mathematics Subject Classification. Primary (35Q30) Secondary(47H20, 76DO3) .

Key words and phrases. Global (in time), 3D-Navier-Stokes Equations.

 $\{\mathbf{u} \in \mathbb{C}_0^{\infty}[\mathbb{R}^3]^3 \mid \nabla \cdot \mathbf{u} = 0\}$  which vanish at infinity with respect to the inner product of  $\mathbb{L}^2[\mathbb{R}^3]^3$ , and let  $\mathbb{V}_0[\mathbb{R}^3]^3$  be the completion of the above functions which vanish at infinity with respect to the inner product of  $\mathbb{H}_0^1[\mathbb{R}^3]$ , the functions in  $\mathbb{H}_0[\mathbb{R}^3]^3$  with weak derivatives in  $(\mathbb{L}^2[\mathbb{R}^3])^3$ . The global in time classical Navier-Stokes initial-value problem (on  $\mathbb{R}^3$  and all T > 0) is to find functions  $\mathbf{u} : [0, T] \times \mathbb{R}^3 \to \mathbb{R}^3$  and  $p : [0, T] \times \mathbb{R}^3 \to \mathbb{R}$  such that

$$\partial_t \mathbf{u} + (\mathbf{u} \cdot \nabla) \mathbf{u} - \nu \Delta \mathbf{u} + \nabla p = \mathbf{f}(t) \text{ in } (0, T) \times \mathbb{R}^3,$$

$$\nabla \cdot \mathbf{u} = 0 \text{ in } (0, T) \times \mathbb{R}^3 \text{ (in the weak sense)},$$

$$\lim_{\|\mathbf{x}\| \to \infty} \mathbf{u}(t, \mathbf{x}) = 0 \text{ on } (0, T) \times \mathbb{R}^3,$$

$$\mathbf{u}(0, \mathbf{x}) = \mathbf{u}_0(\mathbf{x}) \text{ in } \mathbb{R}^3.$$

The equations describe the time evolution of the fluid velocity  $\mathbf{u}(\mathbf{x},t)$  and the pressure p of an incompressible viscous homogeneous Newtonian fluid with constant viscosity coefficient  $\nu$  in terms of a given initial velocity  $\mathbf{u}_0(\mathbf{x})$  and given external body forces  $\mathbf{f}(\mathbf{x},t)$ . (Note that our third condition,  $\lim_{\|\mathbf{x}\|\to\infty}\mathbf{u}(t,\mathbf{x})=0$  on  $(0,T)\times\mathbb{R}^3$ , is natural in this case since it is well-known that  $\mathbb{H}_0^k[\mathbb{R}^3]^3=\mathbb{H}^k[\mathbb{R}^3]^3$  (see Stein [S] or [SY].)

### Purpose

Let  $\mathbb{P}$  be the (Leray) orthogonal projection of  $(\mathbb{L}^2[\mathbb{R}^3])^3$  onto  $\mathbb{H}_0[\mathbb{R}^3]^3$  and define the Stokes operator by:  $\mathbf{A}\mathbf{u} =: -\mathbb{P}\Delta\mathbf{u}$ , for  $\mathbf{u} \in D(\mathbf{A}) \subset \mathbb{H}_0^2[\mathbb{R}^3]^3$ , the domain of  $\mathbf{A}$ . Let  $\mathbf{B}\mathbf{u} =: 1/2\mathbb{P}(-\Delta + |\mathbf{x}|^2)\mathbf{u}$  for  $\mathbf{u} \in D(\mathbf{B})$ . We call  $\mathbf{B}$  the Hermite-Stokes operator. The purpose of this paper is to prove that there exists a number  $\mathbf{u}_+$ , depending only on  $\mathbf{A}$ ,  $\mathbf{B}$ , f,  $\nu$  and  $\mathbb{R}^3$ , such that, for all functions in  $\mathbb{D} = D(\mathbf{A}) \cap \mathbb{B}(\mathbb{R}^3)$ , where GLOBAL (IN TIME) SOLUTIONS TO THE 3D-NAVIER-STOKES EQUATIONS ON  $\mathbb{R}^3$   $\mathbb{B}(\mathbb{R}^3)$  is the closed ball of radius  $\mathbf{u}_+$  in  $\mathbb{H}_0(\mathbb{R}^3)^3$ , the Navier-Stokes equations have unique strong solutions in  $\mathbf{u} \in L^\infty_{\mathrm{loc}}[[0,\infty); \mathbb{V}_0(\mathbb{R}^3)^3] \cap \mathbb{C}^1[(0,\infty); \mathbb{H}_0(\mathbb{R}^3)^3$ .

#### **PRELIMINARIES**

In terms of notation and convention, we follow Sell and You [SY]. In order to simplify notation, we let  $\mathbb{H}$  denote  $\mathbb{H}_0[\mathbb{R}^3]^3$  and  $\mathbb{V}$  denote  $\mathbb{V}_0[\mathbb{R}^3]^3$ . Our use of the Fourier transform follows the definition of Rudin [RU]:  $\mathfrak{F}(h) = \frac{1}{[2\pi]^{3/2}} \int_{\mathbb{R}^3} e^{i\mathbf{x}\cdot\mathbf{y}} h(\mathbf{y}) d\mathbf{y}$ , so that no factors of  $2\pi$  appear in the transform pairs. In order to simplify our proofs, we always assume that all functions  $\mathbf{u}$ ,  $\mathbf{v}$  are in  $D(\mathbf{A})$  and, as in [GZ2], we take  $c = max\{c_i\}$ , where  $c_i$  is one of the nine positive constants that appear on pages 363-367 in [SY]. It will also be convenient to use the fact that the norms of  $\mathbb{V}$  and  $\mathbb{V}^{-1}$  are equivalent in their respective graph norms relative to  $\mathbb{H}$ .

#### THE STOKES OPERATOR

It is known that  $\mathbf{A}$  is a nonnegative linear operator which generates an analytic contraction semigroup. It follows that the fractional powers  $\mathbf{A}^{1/2}$  and  $\mathbf{A}^{-1/2}$  are well defined. Moreover, it is also known (cf., [SY], [T1]) that the norms  $\|\mathbf{A}^{1/2}\mathbf{u}\|_{\mathbb{H}}$  and  $\|\mathbf{A}^{-1/2}\mathbf{u}\|_{\mathbb{H}}$  are equivalent to the corresponding norms induced by the Sobolev space  $(H^1[\mathbb{R}^3])^3$ , so that:

(2) 
$$\|\mathbf{u}\|_{\mathbb{V}} \equiv \|\mathbf{A}^{1/2}\mathbf{u}\|_{\mathbb{H}} \text{ and } \|\mathbf{u}\|_{\mathbb{V}^{-1}} \equiv \|\mathbf{A}^{-1/2}\mathbf{u}\|_{\mathbb{H}}.$$

In addition,  $\mathbf{A}$  is an isomorphism from  $D(\mathbf{A}) \xrightarrow{onto} D(\mathbf{A}^{-1})$ . Furthermore, the embeddings  $\mathbb{V} \to \mathbb{H} \to \mathbb{V}^{-1}$  are continuous, and it is easy to see that  $\mathbf{A}^{-1}$  is the projection of an operator represented by the Riesz potential, mapping  $D(\mathbf{A}^{-1})$ 

onto  $D(\mathbf{A})$  (see Stein [S]). Applying the Leray projection to equation (1), with  $\mathbf{C}(\mathbf{u}, \mathbf{u}) = \mathbb{P}(\mathbf{u} \cdot \nabla)\mathbf{u}$ , we can recast equation (1) in the standard form:

$$\partial_t \mathbf{u} = -\nu \mathbf{A} \mathbf{u} - \mathbf{C}(\mathbf{u}, \mathbf{u}) + \mathbb{P} \mathbf{f}(t) \text{ in } (0, T) \times \mathbb{R}^3,$$

$$\nabla \cdot \mathbf{u} = 0 \text{ in } (0, T) \times \mathbb{R}^3,$$

$$\lim_{\|\mathbf{x}\| \to \infty} \mathbf{u}(t, \mathbf{x}) = 0 \text{ on } (0, T) \times \mathbb{R}^3,$$

$$\mathbf{u}(0, \mathbf{x}) = \mathbf{u}_0(\mathbf{x}) \text{ in } \mathbb{R}^3.$$

where we have used the fact that the orthogonal complement of  $\mathbb{H}[\mathbb{R}^3]$  relative to  $(\mathbb{L}^2)[\mathbb{R}^3])^3$  is  $\{\mathbf{v}: \mathbf{v} = \nabla q, \ q \in (H^1[\mathbb{R}^3])^3\}$  to eliminate the pressure term (see Galdi [GA] or [SY, T1, T2]). Theorem 1 below will be used to get our basic estimate in Theorem 3. This result is a simple extension of the bounded domain case first proved by Constantin and Foias [CF].

**Theorem 1.** Let  $\alpha_i, 1 \le i \le 3$ , satisfy  $0 \le \alpha_1 \le 3$ ,  $0 \le \alpha_2 \le 2$ ,  $0 \le \alpha_3 \le 3$ , with  $\alpha_1 + \alpha_2 + \alpha_3 \ge 3/2$  and

$$(\alpha_1, \alpha_2, \alpha_3) \notin \{(3/2, 0, 0), (0, 3/2, 0), (0, 0, 3/2)\}.$$

Then there is a positive constant  $c = c(\alpha_i)$  such that

$$\left|\left\langle \mathbf{C}(\mathbf{u}, \mathbf{v}), \mathbf{w} \right\rangle_{\mathbb{H}}\right| \leqslant c \left\| \mathbf{A}^{\alpha_1/2} \mathbf{u} \right\|_{\mathbb{H}} \left\| \mathbf{A}^{(1+\alpha_2)/2} \mathbf{v} \right\|_{\mathbb{H}} \left\| \mathbf{A}^{\alpha_3/2} \mathbf{w} \right\|_{\mathbb{H}}.$$

We shall make use of the following interpolation inequality: (see Sell and You [SY], page 363)

$$\|\mathbf{A}^{\gamma}\mathbf{u}\|_{\mathbb{H}} \leqslant c \|\mathbf{A}^{\alpha}\mathbf{u}\|_{\mathbb{H}}^{\theta} \|\mathbf{A}^{\beta}\mathbf{u}\|_{\mathbb{H}}^{(1-\theta)}$$

 $\text{for all } \mathbf{u} \in D(\mathbf{A}^{\alpha}), \text{ where } \gamma = \theta \alpha + (1-\theta)\beta, \ \alpha, \beta, \gamma \in \mathbb{R}, \, 0 \leq \theta \leq 1 \text{ and } \beta \leqslant \alpha.$ 

# Global (in time) solutions to the 3D-navier-stokes equations on $\mathbb{R}^3$

### THE HERMITE-STOKES OPERATOR

The operator  $\hat{\mathbf{B}} = 1/2(-\Delta + |\mathbf{x}|^2)$  is the three-dimensional version of the standard harmonic oscillator operator, which generates the Hermite functions (products of the Hermite polynomials by  $e^{-x^2/2}$ ) as eigenfunctions for the eigenvalue problem on  $\mathbb{R}$ , ( see Hermite [HR], Appell and Kamé de Fériet [AK], and Magnus, Oberhettinger and Soni [MOS]). It is easy to show directly, by separation of variables, that the solution to the 3-dimensional problem is the product of the solutions to the 1-dimensional problem, while the eigenvalues for the 3-dimensional Hermite polynomials are the sums of those for the 1-dimensional polynomials. Furthermore,  $\hat{\mathbf{B}}$ , and hence  $\mathbf{B} = \mathbb{P}\hat{\mathbf{B}}$ , is positive with a compact inverse, while  $\mathbf{A}$  has an unbounded inverse on  $\mathbb{H}_0(\mathbb{R}^3)^3$ . It turns out that  $\hat{\mathbf{B}}$  is "natural" for  $\mathbb{R}^3$  in the sense that it is the only positive self-adjoint (sectorial) operator of lowest degree that is invariant under both rotations and Fourier transformations. (This is actually true for  $\mathbb{R}^n$ ,  $n \geq 1$ .)

We will have need of the fact that every function  $\mathbf{h}(t) \in \mathbb{H}$  has an expansion in terms of the eigenfunctions of  $\mathbf{B}$  so that, for example,  $\mathbf{B}^{-\beta} \mathbf{h}(t) = \sum_{k=1}^{\infty} \lambda_k^{-\beta} h_k(t) \mathbf{e}^k(\mathbf{x})$  and, from here, it is easy to see that  $\|\mathbf{B}^{-\beta} \mathbf{h}(t)\|_{\mathbb{H}} \leq \lambda_1^{-\beta} \|\mathbf{h}(t)\|_{\mathbb{H}}$ , where  $\lambda_1^{-1}$  is the largest eigenvalue of  $\mathbf{B}^{-1}$ . We also need the following result for our basic Theorem.

# Lemma 2. $D(\mathbf{A}) = \mathbf{D}(\mathbf{B})$ .

*Proof.* If we define a norm on  $D(\mathbf{A})$  by  $\|\mathbf{u}\|_{\mathbf{A}} = \|\mathbf{A}\mathbf{u}\|_{\mathbb{H}}$ , then  $(D(\mathbf{A}), \|\cdot\|_{\mathbf{A}})$  is a Hilbert space. Now note that the Fourier transform  $\mathfrak{F}(\cdot)$  is an isometric isomorphism on  $(D(\mathbf{A}), \|\cdot\|_{\mathbf{A}})$  to  $(D(\mathbb{P}|\mathbf{x}|^2), \|\cdot\|_{\mathbf{A}})$ , since  $\|\mathbf{A}\mathbf{u}\|_{\mathbb{H}} = \|\mathfrak{F}(\mathbf{A}\mathbf{u})\|_{\mathbb{H}} = \|\mathbf{F}(\mathbf{A}\mathbf{u})\|_{\mathbb{H}}$ 

 $\|\mathbb{P} \|\mathbf{x}\|^2 \|\mathbf{\hat{u}}\|_{\mathbb{H}}$ . It is now easy to see that  $D(\mathbf{A}) = D(\mathbb{P} \|\mathbf{x}\|^2)$ . From this, it follows that  $D(\mathbf{A}) = D(\mathbf{B})$ .

It follows from the above lemma that  $(\mathbf{AB})^{-\delta}$  is bounded for  $\delta > 0$ . The following estimate is equation 61.24.1 on page 366 in Sell and You [SY]. If we set  $\alpha_1 = 1, \alpha_2 = 1/2$ , and  $\alpha_3 = 0$  in Theorem 1, along with the interpolation inequality, we get that

(4) 
$$|\langle \mathbf{C}(\mathbf{u}, \mathbf{v}), \mathbf{w} \rangle_{\mathbb{H}}| \leqslant c \|\mathbf{A}^{1/2} \mathbf{u}\|_{\mathbb{H}} \|\mathbf{A} \mathbf{v}\|_{\mathbb{H}} \|\mathbf{w}\|_{\mathbb{H}}.$$

**Theorem 3.** Let  $\mathbf{u}, \mathbf{v}, \mathbf{w} \in \mathbb{H}$ , and let  $\varepsilon > 0$  be arbitrary. Then, for  $\delta = 1/4 + \varepsilon/2$ , we have that:

(5) 
$$\left| \left\langle (\mathbf{A}\mathbf{B})^{-(1+\delta)} \mathbf{C}(\mathbf{u}, \mathbf{v}), \mathbf{w} \right\rangle_{\mathbb{H}} \right| \leqslant c \lambda_1^{-(1+\delta)} \|\mathbf{u}\|_{\mathbb{H}} \|\mathbf{v}\|_{\mathbb{H}} \|\mathbf{w}\|_{\mathbb{H}}.$$

*Proof.* Using the self-adjoint property of  $\mathbf{A}$ , and integration by parts, we have

$$\left\langle \mathbf{A}^{-\beta}\mathbf{C}(\mathbf{u},\mathbf{v}),\mathbf{h}\right\rangle_{\mathbb{H}} = \left\langle \mathbf{C}(\mathbf{u},\mathbf{v}),\mathbf{A}^{-\beta}\mathbf{h}\right\rangle_{\mathbb{H}} = -\left\langle \mathbf{C}(\mathbf{u},\mathbf{A}^{-\beta}\mathbf{h}),\mathbf{v}\right\rangle_{\mathbb{H}}.$$

It now follows from Theorem 1 that:

$$\left|\left\langle \mathbf{A}^{-\beta}\mathbf{C}(\mathbf{u},\mathbf{v}),\mathbf{h}\right\rangle_{\mathbb{H}}\right|\leqslant c\left\|\mathbf{A}^{\alpha_{1}/2}\mathbf{u}\right\|_{\mathbb{H}}\left\|\mathbf{A}^{-\beta+(1+\alpha_{2})/2}\mathbf{h}\right\|_{\mathbb{H}}\left\|\mathbf{A}^{\alpha_{3}/2}\mathbf{v}\right\|_{\mathbb{H}}.$$

If we set  $\beta = 1 + \delta$ ,  $\alpha_1 = \alpha_3 = 0$ , we have

$$\left| \left\langle \mathbf{A}^{-(1+\delta)} \mathbf{C}(\mathbf{u}, \mathbf{v}), \mathbf{h} \right\rangle_{\mathbb{H}} \right| \leqslant c \|\mathbf{u}\|_{\mathbb{H}} \|\mathbf{v}\|_{\mathbb{H}} \|\mathbf{A}^{(\alpha_2 - 1 - 2\delta)/2} \mathbf{h} \|_{\mathbb{H}}.$$

With  $\delta = 1/4 + \varepsilon/2$ , we get that, for the last term to reduce to  $\|\mathbf{h}\|_{\mathbb{H}}$ , we can set  $\alpha_2 = 3/2 + \varepsilon$ . It follows that the conditions of Theorem 1 are satisfied if  $3/2 + \varepsilon < 2$ . Thus, it suffices to assume that  $\varepsilon < 1/2$ , which we will do in the rest of the paper

GLOBAL (IN TIME) SOLUTIONS TO THE 3D-NAVIER-STOKES EQUATIONS ON  $\mathbb{R}^7$  without comment. Our proof is completed by taking  $\mathbf{h} = \mathbf{B}^{-\beta}\mathbf{w}$ , and the fact that  $\|\mathbf{B}^{-\beta}\mathbf{w}\|_{\mathbb{H}} \leq \lambda_1^{-\beta} \|\mathbf{w}\|_{\mathbb{H}}$ .

**Example 4.** If we use Theorem 1, with  $\alpha_1 = 5/4$ ,  $\alpha_2 = 1/4$ , and  $\alpha_3 = 0$ , along with the interpolation inequality, and the fact that  $\|\mathbf{A}^{1/2}\mathbf{u}\|_{\mathbb{H}} \leq \|\mathbf{A}\mathbf{u}\|_{\mathbb{H}}$  we have that, for all  $\mathbf{u}, \mathbf{v} \in D(\mathbf{A})$ ,

$$\|\mathbf{C}(\mathbf{u}, \mathbf{v})\|_{\mathbb{H}} \leqslant c \|\mathbf{A}^{1/2}\mathbf{u}\|_{\mathbb{H}}^{3/4} \|\mathbf{A}\mathbf{u}\|_{\mathbb{H}}^{1/4} \|\mathbf{A}^{1/2}\mathbf{v}\|_{\mathbb{H}}^{3/4} \|\mathbf{A}\mathbf{v}\|_{\mathbb{H}}^{1/4}$$

$$\leqslant c \|\mathbf{A}\mathbf{u}\|_{\mathbb{H}} \|\mathbf{A}\mathbf{v}\|_{\mathbb{H}}.$$

A better estimate is possible, but for our use, equation (6) will suffice.

**Definition 5.** We say that the operator  $J(\cdot,t)$  is (for each t)

- (1) 0-Dissipative if  $\langle \mathbf{J}(\mathbf{u},t), \mathbf{u} \rangle_{\mathbb{H}} \leq 0$ .
- (2) Dissipative if  $\langle \mathbf{J}(\mathbf{u},t) \mathbf{J}(\mathbf{v},t), \mathbf{u} \mathbf{v} \rangle_{\mathbb{H}} \leq 0$ .
- (3) Strongly dissipative if there exists an  $\alpha > 0$  such that

$$\langle \mathbf{J}(\mathbf{u},t) - \mathbf{J}(\mathbf{v},t), \mathbf{u} - \mathbf{v} \rangle_{\mathbb{H}} \le -\alpha \|\mathbf{u} - \mathbf{v}\|_{\mathbb{H}}^{2}.$$

(4) Uniformly dissipative if there exists a strictly monotone increasing function a(t) with a(0) = 0,  $\lim_{t \to \infty} a(t) = \infty$ , and:

$$\langle \mathbf{J}(\mathbf{u},t) - \mathbf{J}(\mathbf{v},t), \mathbf{u} - \mathbf{v} \rangle_{\mathbb{H}} \le -a \left( \|\mathbf{u} - \mathbf{v}\|_{\mathbb{H}} \right) \|\mathbf{u} - \mathbf{v}\|_{\mathbb{H}}.$$

Note that, if  $\mathbf{J}(\cdot,t)$  is a linear operator, definitions 1) and 2) coincide. Theorem 6 below is essentially due to Browder [B], see Zeidler [Z, Corollary 32.27, page 868 and Corollary 32.35, page 887 in, Vol. IIB], while Theorem 7 is from Miyadera [M, p. 185, Theorem 6.20], and is a modification of the Crandall-Liggett Theorem [CL] (see the appendix to the first section of [CL]).

Theorem 6. Let B[R 3 ] be a closed, bounded, convex subset of H[R 3 ]. If J(·, t) : B[R 3 ] <sup>→</sup> <sup>H</sup>[<sup>R</sup> 3 ] is closed and strongly dissipative for each fixed t ≥ 0 then, for each <sup>b</sup> <sup>∈</sup> <sup>B</sup>[<sup>R</sup> 3 ], there is a <sup>u</sup> <sup>∈</sup> <sup>B</sup>[<sup>R</sup> 3 ] with J(u, t) = b (e.g., the range, Ran[J(·, t)] ⊃ B[R 3 ]).

Theorem 7. Let { A(t), t <sup>∈</sup> <sup>I</sup> = [0, <sup>∞</sup>)} be a family of operators defined on <sup>H</sup>[<sup>R</sup> 3 with domains <sup>D</sup>(A(t)) = <sup>D</sup>, independent of <sup>t</sup>. We assume that <sup>D</sup> <sup>=</sup> <sup>D</sup> <sup>∩</sup> <sup>B</sup>[<sup>R</sup> 3 ] is a closed convex set (in an appropriate topology):

- (1) The operator A(t) is the generator of a contraction semigroup for each t ∈ I.
- (2) The function <sup>A</sup>(t)<sup>u</sup> is continuous in both variables on <sup>I</sup> <sup>×</sup> <sup>D</sup>.

Then, for every <sup>u</sup><sup>0</sup> <sup>∈</sup> <sup>D</sup>, the problem <sup>∂</sup>tu(t, <sup>x</sup>) = <sup>A</sup>(t)u(t, <sup>x</sup>), <sup>u</sup>(0, <sup>x</sup>) = <sup>u</sup>0(x), has a unique solution <sup>u</sup>(t, <sup>x</sup>) <sup>∈</sup> <sup>C</sup> 1 (I; D).

# M-Dissipative Conditions

Let us assume that <sup>f</sup>(t) <sup>∈</sup> <sup>L</sup><sup>∞</sup>[[0, <sup>∞</sup>); <sup>H</sup>] and is Lipschitz continuous in <sup>t</sup>, with kf(t) − f(τ)k<sup>H</sup> ≤ d |t − τ | θ , d > 0, 0 < θ < 1. With δ as in Theorem 3, we can rewrite equation (3) in the form:

$$\partial_{t}\mathbf{u} = \nu(\mathbf{A}\mathbf{B})^{1+\delta}\mathbf{J}(\mathbf{u}, t) \text{ in } (0, T) \times \Omega,$$
(7)
$$\mathbf{J}(\mathbf{u}, t) = -\mathbf{B}^{-(1+\delta)}\mathbf{A}^{-\delta}\mathbf{u} - \nu^{-1}(\mathbf{A}\mathbf{B})^{-(1+\delta)}\mathbf{C}(\mathbf{u}, \mathbf{u}) + \nu^{-1}(\mathbf{A}\mathbf{B})^{-(1+\delta)}\mathbb{P}\mathbf{f}(t).$$

## Approach

We begin with a study of the operator J(·, t), for fixed t, and seek conditions depending on A, B, ν, and f(t) which guarantee that J(·, t) is m-dissipative for each t. Clearly J(·, t) : D[(AB) (1+δ) onto −−−→ <sup>D</sup>[(AB) (1+δ) ] and, since ν(AB) (1+δ)

GLOBAL (IN TIME) SOLUTIONS TO THE 3D-NAVIER-STOKES EQUATIONS ON  $\mathbb{R}^3$  is a closed positive (m-accretive) operator (so that  $-(\mathbf{AB})^{(1+\delta)}$  generates a linear contraction semigroup), we expect that  $\nu(\mathbf{AB})^{(1+\delta)}J(\cdot,t)$  will be m-dissipative for each t.

**Theorem 8.** For  $t \in I = [0, \infty)$  and, for each fixed  $\mathbf{u} \in \mathbb{H}$ ,  $\mathbf{J}(\mathbf{u}, t)$  is Lipschitz continuous, with  $\|\mathbf{J}(\mathbf{u}, t) - \mathbf{J}(\mathbf{u}, \tau)\|_{\mathbb{H}} \le d' |t - \tau|^{\theta}$ , where  $d' = d\nu^{-1}a^{-(1+\delta)}$ , d is the Lipschitz constant for the function  $\mathbf{f}(t)$  and  $a^{-(1+\delta)} = \|(\mathbf{A}\mathbf{B})^{-(1+\delta)}\|_{\mathbb{H}}$ .

*Proof.* For fixed  $\mathbf{u} \in \mathbb{H}$ ,

$$\|\mathbf{J}(\mathbf{u},t) - \mathbf{J}(\mathbf{u},\tau)\|_{\mathbb{H}} = \nu^{-1} \|(\mathbf{A}\mathbf{B})^{-(1+\delta)} [\mathbb{P}\mathbf{f}(t) - \mathbb{P}\mathbf{f}(\tau)]\|_{\mathbb{H}}$$
$$\leq d\nu^{-1} a^{-(1+\delta)} |t - \tau|^{\theta} = d' |t - \tau|^{\theta}.$$

### Main Results

**Theorem 9.** Let  $f = \sup_{t \in \mathbf{R}^+} \| \mathbb{P} \mathbf{f}(t) \|_{\mathbb{H}} < \infty$ , then there exists a positive constant  $\mathbf{u}_+$ , depending only on f,  $\mathbf{A}$ ,  $\mathbf{B}$  and  $\nu$  such that, for all  $\mathbf{u}$  with  $\| \mathbf{u} \|_{\mathbb{H}} \leq \mathbf{u}_+$ ,  $\mathbf{J}(\cdot, t)$  is strongly dissipative.

*Proof.* The proof of our first assertion has two parts. First, we require that the nonlinear operator  $\mathbf{J}(\cdot,t)$  be 0-dissipative, which gives us an upper bound  $\mathbf{u}_+$  in terms of the norm (e.g.,  $\|\mathbf{u}\|_{\mathbb{H}} \leq \mathbf{u}_+$ ). We then use this part, and the fact that  $\|\mathbf{u}\|_{\mathbb{H}} \leq \|\mathbf{A}\mathbf{u}\|_{\mathbb{H}}$ , to show that  $\mathbf{J}(\cdot,t)$  is strongly dissipative on the closed ball,  $\mathbb{B}_+ = \{\mathbf{u} \in \mathbb{H} : \|\mathbf{A}\mathbf{u}\|_{\mathbb{H}} \leq (1/2)\mathbf{u}_+\}.$ 

Part 1) From equation (5), we consider the expression

$$\begin{split} \left\langle \mathbf{J}(\mathbf{u},t), (\mathbf{A}\mathbf{B})^{-\delta}\mathbf{u} \right\rangle_{\mathbb{H}} &= -\left\langle \mathbf{B}^{-1}(\mathbf{A}\mathbf{B})^{-\delta}\mathbf{u}, (\mathbf{A}\mathbf{B})^{-\delta}\mathbf{u} \right\rangle_{\mathbb{H}} \\ &+ \nu^{-1} \left\langle -(\mathbf{A}\mathbf{B})^{-(1+\delta)}\mathbf{C}(\mathbf{u},\mathbf{u}) + (\mathbf{A}\mathbf{B})^{-(1+\delta)}\mathbb{P}\mathbf{f}(t), (\mathbf{A}\mathbf{B})^{-\delta}\mathbf{u} \right\rangle_{\mathbb{H}} \\ &= -\left\| \mathbf{B}^{-1/2}(\mathbf{A}\mathbf{B})^{-\delta}\mathbf{u} \right\|_{\mathbb{H}}^{2} - \nu^{-1} \left\langle (\mathbf{A}\mathbf{B})^{-(1+\delta)}\mathbf{C}(\mathbf{u},\mathbf{u}), (\mathbf{A}\mathbf{B})^{-\delta}\mathbf{u} \right\rangle_{\mathbb{H}} + \nu^{-1} \left\langle (\mathbf{A}\mathbf{B})^{-(1+\delta)}\mathbb{P}\mathbf{f}(t), (\mathbf{A}\mathbf{B})^{-\delta}\mathbf{u} \right\rangle_{\mathbb{H}} \\ &= -\left\| \mathbf{B}^{-1/2}(\mathbf{A}\mathbf{B})^{-\delta}\mathbf{u} \right\|_{\mathbb{H}}^{2} - \nu^{-1} \left\langle \mathbf{C}((\mathbf{A}\mathbf{B})^{-(1+\delta)}\mathbf{u}, \mathbf{u}), (\mathbf{A}\mathbf{B})^{-\delta}\mathbf{u} \right\rangle_{\mathbb{H}} + \nu^{-1} \left\langle (\mathbf{A}\mathbf{B})^{-(1+\delta)}\mathbb{P}\mathbf{f}(t), (\mathbf{A}\mathbf{B})^{-\delta}\mathbf{u} \right\rangle_{\mathbb{H}}. \end{split}$$

It follows that

$$\begin{split} \left\langle \mathbf{J}(\mathbf{u},t), (\mathbf{A}\mathbf{B})^{-\delta}\mathbf{u} \right\rangle_{\mathbb{H}} \leqslant &- \left\| \mathbf{B}^{-1/2} (\mathbf{A}\mathbf{B})^{-\delta}\mathbf{u} \right\|_{\mathbb{H}}^2 + \nu^{-1} \left| \left\langle \mathbf{C}((\mathbf{A}\mathbf{B})^{-(1+\delta)}\mathbf{u}, \mathbf{u}), (\mathbf{A}\mathbf{B})^{-\delta}\mathbf{u} \right\rangle_{\mathbb{H}} \right| \\ &+ \nu^{-1} a^{-(1+\delta)} f \left\| (\mathbf{A}\mathbf{B})^{-\delta}\mathbf{u} \right\|_{\mathbb{H}} \\ \leqslant &- \left\| \mathbf{B}^{-1/2} (\mathbf{A}\mathbf{B})^{-\delta}\mathbf{u} \right\|_{\mathbb{H}}^2 + c a^{-\delta} (\nu \lambda_1^{(1+\delta)})^{-1} \left\| \mathbf{u} \right\|_{\mathbb{H}}^3 + \nu^{-1} a^{-(1+2\delta)} f \left\| \mathbf{u} \right\|_{\mathbb{H}}. \end{split}$$

In the last line, we used our estimate from Theorem 3. We now choose the first eigenvalue  $\lambda_n$ ,  $n \geq 1$ , and number  $\omega$  such that

$$(1) \ \lambda_n^{-1/2} a^{-\delta} \|\mathbf{u}\|_{\mathbb{H}} \leqslant \left\| \mathbf{B}^{-1/2} (\mathbf{A} \mathbf{B})^{-\delta} \mathbf{u} \right\|_{\mathbb{H}} \leqslant \lambda_1^{-1/2} a^{-\delta} \|\mathbf{u}\|_{\mathbb{H}} \,,$$

(2) 
$$\lambda_1^{-\omega/2} a^{-\delta} \|\mathbf{u}\|_{\mathbb{H}} \leqslant \|\mathbf{B}^{-1/2} (\mathbf{A} \mathbf{B})^{-\delta} \mathbf{u}\|_{\mathbb{H}} \leqslant \lambda_1^{-1/2} a^{-\delta} \|\mathbf{u}\|_{\mathbb{H}}$$
,

and let  $\lambda_0^{-1} = \max\{\lambda_n^{-1}, \lambda_1^{-\omega}\}$ . It then follows that  $-\lambda_0^{-1}a^{-2\delta} \|\mathbf{u}\|_{\mathbb{H}}^2 \geqslant$  $-\|\mathbf{B}^{-1/2}(\mathbf{A}\mathbf{B})^{-\delta}\mathbf{u}\|_{\mathbb{H}}^2$ . Thus,  $\mathbf{J}(\cdot, t)$  will be 0-dissipative if

$$-\lambda_0^{-1}a^{-2\delta} \|\mathbf{u}\|_{\mathbb{H}}^2 + ca^{-\delta} (\nu \lambda_1^{(1+\delta)})^{-1} \|\mathbf{u}\|_{\mathbb{H}}^3 + (\nu a^{(1+2\delta)})^{-1} f \|\mathbf{u}\|_{\mathbb{H}} \leqslant 0,$$

so that

$$(8) \quad a^{-\delta} \|\mathbf{u}\|_{\mathbb{H}} \left[ c(\nu \lambda_{1}^{(1+\delta)})^{-1} \|\mathbf{u}\|_{\mathbb{H}}^{2} - \lambda_{0}^{-1} a^{-\delta} \|\mathbf{u}\|_{\mathbb{H}} + (\nu a^{(1+\delta)})^{-1} f \right] \leqslant 0.$$

Since  $\|\mathbf{u}\|_{\mathbb{H}} > 0$ , we have that  $\mathbf{J}(\cdot, t)$  is 0-dissipative if

$$c(\nu \lambda_1^{(1+\delta)})^{-1} \|\mathbf{u}\|_{\mathbb{H}}^2 - \lambda_0^{-1} a^{-\delta} \|\mathbf{u}\|_{\mathbb{H}} + (\nu a^{(1+\delta)})^{-1} f \leqslant 0.$$

GLOBAL (IN TIME) SOLUTIONS TO THE 3D-NAVIER-STOKES EQUATIONS ON  $\mathbb{R}^3$  Solving, we get that

$$\mathbf{u}_{\pm} = \frac{\nu \lambda_{1}^{1+\delta}}{2c\lambda_{0}a^{\delta}} \left\{ 1 \pm \sqrt{1 - (4c\lambda_{0}^{2}f) / (\nu^{2}a^{(1-\delta)}\lambda_{1}^{(1+\delta)})} \right\} = \frac{\nu \lambda_{1}^{1+\delta}}{2c\lambda_{0}a^{\delta}} \left\{ 1 \pm \sqrt{1 - \gamma} \right\},$$

where  $\gamma=(4c\lambda_0^2f)\Big/(\nu^2a^{(1-\delta)}\lambda_1^{(1+\delta)})$ . Since we want real distinct solutions, we must require that

$$\gamma = (4c\lambda_0^2 f) / (\nu^2 a^{(1-\delta)} \lambda_1^{(1+\delta)}) < 1 \Rightarrow \nu^2 a^{(1-\delta)} \lambda_1^{(1+\delta)} > 4c\lambda_0^2 f$$

$$\Rightarrow \nu > 2\lambda_0 a^{-(1-\delta)/2} \lambda_1^{-(1+\delta)/2} (cf)^{1/2}.$$

It follows that, if  $\mathbb{P}\mathbf{f} \neq \mathbf{0}$ , then  $\mathbf{u}_{-} < \mathbf{u}_{+}$ , and our requirement that  $\mathbf{J}$  is 0-dissipative implies that, since our solution factors as  $(\|\mathbf{u}\|_{\mathbb{H}} - \mathbf{u}_{+})(\|\mathbf{u}\|_{\mathbb{H}} - \mathbf{u}_{-}) \leq 0$ , we must have that:

$$\|\mathbf{u}\|_{\mathbb{H}} - \mathbf{u}_{+} \le 0, \ \|\mathbf{u}\|_{\mathbb{H}} - \mathbf{u}_{-} \ge 0.$$

First observe that terms of the form  $(\mathbf{AB})^{-\delta}\mathbf{u}$  are dense. Then note that  $\mathbf{J}(\mathbf{u},t)$  is closed, and the dissipative nature of an operator is determined on a dense set. It follows that, for  $\mathbf{u}_{-} \leq \|\mathbf{u}\|_{\mathbb{H}} \leq \mathbf{u}_{+}$ ,  $\langle \mathbf{J}(\mathbf{u},t), \mathbf{u} \rangle_{\mathbb{H}} \leq 0$ . (It is clear that, when  $\mathbb{P}\mathbf{f}(t) = \mathbf{0}$ ,  $\mathbf{u}_{-} = \mathbf{0}$ , and  $\mathbf{u}_{+} = \nu(c\lambda_{0}a^{\delta})^{-1}\lambda_{1}^{(1+\delta)}$ .)

Part 2): Now, for any  $\mathbf{u}, \mathbf{v} \in \mathbb{H}$  with max( $\|\mathbf{A}\mathbf{u}\|_{\mathbb{H}}, \|\mathbf{A}\mathbf{v}\|_{\mathbb{H}}$ )  $\leq (1/2)\mathbf{u}_{+}$ , we have that

$$\begin{split} \left\langle \mathbf{J}(\mathbf{u},t) - \mathbf{J}(\mathbf{v},t), (\mathbf{A}\mathbf{B})^{-\delta}(\mathbf{u} - \mathbf{v}) \right\rangle_{\mathbb{H}} &= -\left\| \mathbf{B}^{-1/2}(\mathbf{A}\mathbf{B})^{-\delta}(\mathbf{u} - \mathbf{v}) \right\|_{\mathbb{H}}^{2} \\ &- \nu^{-1} \left\langle (\mathbf{A}\mathbf{B})^{-(1+\delta)} [\mathbf{C}(\mathbf{u}, \mathbf{u} - \mathbf{v}) + \mathbf{C}(\mathbf{v}, \mathbf{u} - \mathbf{v})], (\mathbf{A}\mathbf{B})^{-\delta}(\mathbf{u} - \mathbf{v}) \right\rangle_{\mathbb{H}} \\ &\leq -\lambda_{0}^{-1} a^{-2\delta} \|\mathbf{u} - \mathbf{v}\|_{\mathbb{H}}^{2} + ca^{-\delta} \nu^{-1} \lambda_{1}^{-(1+\delta)} \|\mathbf{u} - \mathbf{v}\|_{\mathbb{H}}^{2} (\|\mathbf{u}\|_{\mathbb{H}} + \|\mathbf{v}\|_{\mathbb{H}}) \\ &\leq -\lambda_{0}^{-1} a^{-2\delta} \|\mathbf{u} - \mathbf{v}\|_{\mathbb{H}}^{2} + ca^{-\delta} \nu^{-1} \lambda_{1}^{-(1+\delta)} \|\mathbf{u} - \mathbf{v}\|_{\mathbb{H}}^{2} \mathbf{u}_{+} \\ &= -\lambda_{0}^{-1} a^{-2\delta} \|\mathbf{u} - \mathbf{v}\|_{\mathbb{H}}^{2} + ca^{-\delta} \nu^{-1} \lambda_{1}^{-(1+\delta)} \|\mathbf{u} - \mathbf{v}\|_{\mathbb{H}}^{2} \left(\frac{1}{2} \nu \lambda_{1}^{(1+\delta)} (c^{-1} a^{-\delta} \lambda_{0}^{-1}) \left\{ 1 + \sqrt{1 - \gamma} \right\} \right) \\ &= -\frac{1}{2} \lambda_{0}^{-1} a^{-2\delta} \|\mathbf{u} - \mathbf{v}\|_{\mathbb{H}}^{2} \left\{ 1 - \sqrt{1 - \gamma} \right\} \\ &= -\alpha \|\mathbf{u} - \mathbf{v}\|_{\mathbb{H}}^{2}, \ \alpha = \frac{1}{2} \lambda_{0}^{-1} a^{-2\delta} \left\{ 1 - \sqrt{1 - \gamma} \right\}. \end{split}$$

**Theorem 10.** The operator  $\mathcal{A}(t) = \nu \mathbf{A}^{(1+\delta)} \mathbf{J}(\cdot,t)$  is closed, uniformly dissipative and jointly continuous in  $\mathbf{u}$  and t. Furthermore, for each  $t \in \mathbf{R}^+$  and  $\beta > 0$ ,  $Ran[I - \beta \mathcal{A}(t)] \supset \mathbb{B}[\Omega]$ , so that  $\mathcal{A}(t)$  is m-dissipative on  $\mathbb{D}$ .

*Proof.* Since  $\mathbf{J}(\cdot,t)$  is strongly dissipative and closed on  $\mathbb{B}$ , it follows from Theorem 6 that  $Ran[\mathbf{J}(\cdot,t)] \supset \mathbb{B}$ .

To show that  $\mathcal{A}(t) = \nu(\mathbf{A}\mathbf{B})^{(1+\delta)}\mathbf{J}(\cdot,t)$  is uniformly dissipative for  $\mathbf{u}, \mathbf{v} \in \mathbb{B}_+$ , we have

$$\langle \mathcal{A}(t)\mathbf{u} - \mathcal{A}(t)\mathbf{v}, (\mathbf{u} - \mathbf{v}) \rangle_{\mathbb{H}} = -\nu \left\| \mathbf{A}^{1/2}(\mathbf{u} - \mathbf{v}) \right\|_{\mathbb{H}}^{2}$$
$$-\langle (1/2)[\mathbf{C}(\mathbf{u} - \mathbf{v}, \mathbf{u}) + \mathbf{C}(\mathbf{u} - \mathbf{v}, \mathbf{v})], (\mathbf{u} - \mathbf{v}) \rangle_{\mathbb{H}}.$$

GLOBAL (IN TIME) SOLUTIONS TO THE 3D-NAVIER-STOKES EQUATIONS ON  $\mathbb{R}^3$ Now, from equation (4),

$$\begin{split} |\langle [\mathbf{C}(\mathbf{u} - \mathbf{v}, \mathbf{u}) + \mathbf{C}(\mathbf{u} - \mathbf{v}, \mathbf{v})], (\mathbf{u} - \mathbf{v}) \rangle_{\mathbb{H}}| \\ &\leq c \left\| \mathbf{A}^{1/2} (\mathbf{u} - \mathbf{v}) \right\|_{\mathbb{H}} \left\| (\mathbf{u} - \mathbf{v}) \right\|_{\mathbb{H}} \left\{ \left\| \mathbf{A} \mathbf{u} \right\|_{\mathbb{H}} + \left\| \mathbf{A} \mathbf{v} \right\|_{\mathbb{H}} \right\}. \end{split}$$

We now use  $-\lambda_0^{-1}a^{-\delta} \|(\mathbf{u} - \mathbf{v})\|_{\mathbb{H}} \geqslant -\|\mathbf{A}^{1/2}(\mathbf{u} - \mathbf{v})\|_{\mathbb{H}}$ , and the fact that the first eigenvalue of  $\mathbf{B}$  is 1/2, so that  $\lambda_1^{1+\delta} < 1$ , to get:

$$\begin{split} & \langle \mathcal{A}(t)\mathbf{u} - \mathcal{A}(t)\mathbf{v}, \mathbf{u} - \mathbf{v} \rangle_{\mathbb{H}} \leqslant -\nu \left\| \mathbf{A}^{1/2}(\mathbf{u} - \mathbf{v}) \right\|_{\mathbb{H}}^{2} + \frac{1}{2}c \left\| \mathbf{A}^{1/2}(\mathbf{u} - \mathbf{v}) \right\|_{\mathbb{H}} \left\{ \left\| \mathbf{A}\mathbf{u} \right\|_{\mathbb{H}} + \left\| \mathbf{A}\mathbf{v} \right\|_{\mathbb{H}} \right\} \\ & = \left\| \mathbf{A}^{1/2}(\mathbf{u} - \mathbf{v}) \right\|_{\mathbb{H}} \left\{ -\nu \left\| \mathbf{A}^{1/2}(\mathbf{u} - \mathbf{v}) \right\|_{\mathbb{H}} + \frac{1}{2}c \left\| \mathbf{u} - \mathbf{v} \right\|_{\mathbb{H}} \left[ \left\| \mathbf{A}\mathbf{u} \right\|_{\mathbb{H}} + \left\| \mathbf{A}\mathbf{v} \right\|_{\mathbb{H}} \right] \right\} \\ & \leq \left\| \mathbf{A}^{1/2}(\mathbf{u} - \mathbf{v}) \right\|_{\mathbb{H}} \left\| \mathbf{u} - \mathbf{v} \right\|_{\mathbb{H}} \left\{ -\nu \lambda_{0}^{-1} a^{-\delta} + c \mathbf{u}_{+} \right\} \\ & \leq \left\| \mathbf{A}^{1/2}(\mathbf{u} - \mathbf{v}) \right\|_{\mathbb{H}} \left\| \mathbf{u} - \mathbf{v} \right\|_{\mathbb{H}} \left\{ -\nu \lambda_{0}^{-1} a^{-\delta} + \frac{1}{2}\nu \lambda_{1}^{(1+\delta)} \lambda_{0}^{-1} a^{-\delta} \left[ 1 + \sqrt{1 - \gamma} \right] \right\} \\ & < \frac{1}{2}\nu \lambda_{0}^{-1} a^{-\delta} \left\| \mathbf{A}^{1/2}(\mathbf{u} - \mathbf{v}) \right\|_{\mathbb{H}} \left\| \mathbf{u} - \mathbf{v} \right\|_{\mathbb{H}} \left\{ -1 + \sqrt{1 - \gamma} \right\} < 0. \end{split}$$

If we set  $a(\|(\mathbf{u} - \mathbf{v})\|_{\mathbb{H}}) = -\frac{1}{2}\nu\lambda_0^{-1}a^{-\delta}\left[-1 + \sqrt{1-\gamma}\right]\|\mathbf{A}^{1/2}(\mathbf{u} - \mathbf{v})\|_{\mathbb{H}}$ , we have that:

$$\langle \mathcal{A}(t)\mathbf{u} - \mathcal{A}(t)\mathbf{v}, \mathbf{u} - \mathbf{v} \rangle_{\mathbb{H}} \leqslant -a \left( \|(\mathbf{u} - \mathbf{v})\|_{\mathbb{H}} \right) \|(\mathbf{u} - \mathbf{v})\|_{\mathbb{H}}.$$

It follows that  $\mathcal{A}(t)$  is uniformly dissipative. Since  $-\mathbf{A}^{(1+\delta)}$  is m-dissipative, for  $\beta > 0$ ,  $Ran(I+\beta(\mathbf{AB})^{(1+\delta)}) = \mathbb{H}$ . As  $\mathbf{J}$  is strongly dissipative (in the ball of radius  $\frac{1}{2}\mathbf{u}_+$ ) and closed, with  $Ran[\mathbf{J}] \supset \mathbb{B}$ , and  $\mathbf{J}(\cdot,t) : \mathbb{D} \xrightarrow{onto} \mathbb{D}$ ,  $\mathcal{A}(t)$  is maximal dissipative (in the ball of radius  $\frac{1}{2}\mathbf{u}_+$ ), and also closed, so that  $Ran[I-\beta\mathcal{A}(t)] \supset \mathbb{B}$ . It follows that  $\mathcal{A}(t)$  is m-dissipative on  $\mathbb{B}$  for each  $t \in \mathbf{R}^+$  (since  $\mathbb{H}$  is a Hilbert space). To see that  $\mathcal{A}(t)\mathbf{u}$  is continuous in both variables, let  $\mathbf{u}_n, \mathbf{u} \in \mathbb{B}_+$ ,  $\|\mathbf{A}(\mathbf{u}_n - \mathbf{u})\|_{\mathbb{H}} \to 0$ ,

with tn, t ∈ I and t<sup>n</sup> → t. Then (see equation (6))

$$\|\mathcal{A}(t_n)\mathbf{u}_n - \mathcal{A}(t)\mathbf{u}\|_{\mathbb{H}} \leq \|\mathcal{A}(t_n)\mathbf{u} - \mathcal{A}(t)\mathbf{u}\|_{\mathbb{H}} + \|\mathcal{A}(t_n)\mathbf{u}_n - \mathcal{A}(t_n)\mathbf{u}\|_{\mathbb{H}}$$

$$= \|[\mathbb{P}\mathbf{f}(t_n) - \mathbb{P}\mathbf{f}(t)]\|_{\mathbb{H}} + \|\nu\mathbf{A}(\mathbf{u}_n - \mathbf{u}) + [\mathbf{C}(\mathbf{u}_n - \mathbf{u}, \mathbf{u}_n) + \mathbf{C}(\mathbf{u}, \mathbf{u}_n - \mathbf{u})]\|_{\mathbb{H}}$$

$$\leq d |t_n - t|^{\theta} + \nu \|\mathbf{A}(\mathbf{u}_n - \mathbf{u})\|_{\mathbb{H}} + \|\mathbf{C}(\mathbf{u}_n - \mathbf{u}, \mathbf{u}_n) + \mathbf{C}(\mathbf{u}, \mathbf{u}_n - \mathbf{u})\|_{\mathbb{H}}$$

$$\leq d |t_n - t|^{\theta} + \nu \|\mathbf{A}(\mathbf{u}_n - \mathbf{u})\|_{\mathbb{H}} + c \|\mathbf{A}(\mathbf{u}_n - \mathbf{u})\|_{\mathbb{H}} \{\|\mathbf{A}\mathbf{u}_n\|_{\mathbb{H}} + \|\mathbf{A}\mathbf{u}\|_{\mathbb{H}}\}$$

$$\leq d |t_n - t|^{\theta} + \nu \|\mathbf{A}(\mathbf{u}_n - \mathbf{u})\|_{\mathbb{H}} + 2c \|\mathbf{A}(\mathbf{u}_n - \mathbf{u})\|_{\mathbb{H}} \mathbf{u}_{+}.$$

It follows that A(t)u is continuous in both variables.

Since <sup>B</sup><sup>+</sup> is the closure of <sup>D</sup> <sup>=</sup> <sup>D</sup>(A) <sup>∩</sup> <sup>B</sup> equipped with the restriction of the graph norm of A induced on D(A), it follows that B<sup>+</sup> is a closed, bounded, convex set. We now have:

Theorem 11. For each <sup>T</sup> <sup>∈</sup> <sup>R</sup><sup>+</sup>, <sup>t</sup> <sup>∈</sup> (0, T ) and <sup>u</sup><sup>0</sup> <sup>∈</sup> <sup>D</sup> <sup>⊂</sup> <sup>B</sup>, the global in time Navier-Stokes initial-value problem in R 3 :

$$\partial_{t}\mathbf{u} + (\mathbf{u} \cdot \nabla)\mathbf{u} - \nu\Delta\mathbf{u} + \nabla p = \mathbf{f}(t) \ in \ (0, T) \times \mathbb{R}^{3},$$

$$\nabla \cdot \mathbf{u} = 0 \ in \ (0, T) \times \mathbb{R}^{3},$$

$$\lim_{\|\mathbf{x}\| \to \infty} \mathbf{u}(t, \mathbf{x}) = \mathbf{0} \ on \ (0, T) \times \mathbb{R}^{3},$$

$$\mathbf{u}(0, \mathbf{x}) = \mathbf{u}_{0}(\mathbf{x}) \ in \ \mathbb{R}^{3},$$

has a unique strong solution u(t, x), which is in L 2 loc[[0, <sup>∞</sup>); <sup>H</sup><sup>2</sup> ] and in L<sup>∞</sup> loc[[0, <sup>∞</sup>); <sup>V</sup>] <sup>∩</sup> <sup>C</sup> 1 [(0, <sup>∞</sup>); <sup>H</sup>].

Proof. Theorem 7 allows us to conclude that, when <sup>u</sup><sup>0</sup> <sup>∈</sup> <sup>D</sup>, the initial value problem is solved and the solution u(t, x) is in C 1 [(0, <sup>∞</sup>); <sup>D</sup>]. Since <sup>D</sup> <sup>⊂</sup> <sup>H</sup><sup>2</sup> , it follows that GLOBAL (IN TIME) SOLUTIONS TO THE 3D-NAVIER-STOKES EQUATIONS ON  $\mathbb{R}^3$   $\mathbf{u}(t,\mathbf{x})$  is also in  $\mathbb{V}$ , for each t>0. It is now clear that, for any T>0,

$$\int_0^T \left\| \mathbf{u}(t,\mathbf{x}) \right\|_{\mathbb{H}}^2 dt < \infty, \text{ and } \sup_{0 < t < T} \left\| \mathbf{u}(t,\mathbf{x}) \right\|_{\mathbb{V}}^2 < \infty.$$

This gives our conclusion.

#### DISCUSSION

It is known that, if  $\mathbf{u}_0 \in \mathbb{V}$ , and  $\mathbf{f}(t)$  is  $L^{\infty}[(0,\infty),\mathbb{H}]$  then there is a time T > 0 such that a weak solution with this data is uniquely determined on any subinterval of [0,T) (see Sell and You, page 396, [SY]). Thus, we also have that:

Corollary 12. For each  $t \in \mathbf{R}^+$  and  $\mathbf{u}_0 \in \mathbb{D}$  the Navier-Stokes initial-value problem on  $\mathbb{R}^3$ :

$$\partial_{t}\mathbf{u} + (\mathbf{u} \cdot \nabla)\mathbf{u} - \nu\Delta\mathbf{u} + \nabla p = \mathbf{f}(t) \ in \ (0, T) \times \mathbb{R}^{3},$$

$$\nabla \cdot \mathbf{u} = 0 \ in \ (0, T) \times \mathbb{R}^{3},$$

$$\lim_{\|\mathbf{x}\| \to \infty} \mathbf{u}(t, \mathbf{x}) = \mathbf{0} \ on \ (0, T) \times \mathbb{R}^{3},$$

$$\mathbf{u}(0, \mathbf{x}) = \mathbf{u}_{0}(\mathbf{x}) \ in \ \mathbb{R}^{3}.$$

has a unique weak solution  $\mathbf{u}(t,\mathbf{x})$ , which is in  $L^2_{loc}[[0,\infty);\mathbb{H}^2]$  and in  $L^\infty_{loc}[[0,\infty);\mathbb{V}]\cap\mathbb{C}^1[(0,\infty);\mathbb{H}].$ 

Since we require that our initial data be in  $\mathbb{H}^2$ , the conditions for the Leray-Hopf weak solutions are not satisfied. However, it was an open question as to whether these solutions developed singularities, even if  $\mathbf{u}_0 \in \mathbb{C}_0^{\infty}$  (see Giga [G] and references therein). The above Corollary shows that it suffices that  $\mathbf{u}_0(\mathbf{x}) \in \mathbb{H}^2$  to insure that the solutions develop no singularities.

Acknowledgements. We would like to thank Professor George Sell for his constructive remarks on an earlier draft. We have benefited from his friendship, encouragement and the generous sharing of his knowledge over the last fifteen years. We would also like to sincerely thank Professor Edriss Titi for comments which helped us improve our presentation and eliminate a few areas of possible confusion.

# References

- <span id="page-15-4"></span>[B] F. Browder, Nonlinear operators and nonlinear equations of evolution in Banach spaces, Proc. Sympos. Pure Math., Vol. 18 part II, Amer. Math. Soc., Providence, RI, 1970.
- <span id="page-15-2"></span>[CF] P. Constantin and C. Foi¸as, Navier-Stokes Equations, University of Chicago Press, Chicago, IL, 1988.
- <span id="page-15-5"></span>[CL] M. Crandall and T. Liggett, Generation of semigroups of nonlinear transformations on general Banach spaces, Amer. J. Math. 93 (1971), 265-293.
- [GA] G. P. Galdi, An introduction to the mathematical theory of the Navier-Stokes equations, 2nd Edition, Vol. II, Springer Tracts in Natural Philosophy, Vol. 39, Springer, New York, 1997.
- <span id="page-15-6"></span>[G] Y. Giga, Solutions for semilinear parabolic equations in L<sup>p</sup> and regularity of weak solutions of the Navier-Stokes system, J. Diff. Eq. 62 (1986), 186-212.
- <span id="page-15-0"></span>[GZ1] T. L. Gill and W. W. Zachary, Sufficiency Class for Global (in time) Solutions to the 3D-Navier-Stokes Equations, (submitted) Annals of Mathematics.
- <span id="page-15-1"></span>[GZ2] T. L. Gill and W. W. Zachary, Sufficiency Class for Global (in time) Solutions to the 3D-Navier-Stokes Equations in V, (submitted) Journal of Differential Equations.
- <span id="page-15-3"></span>[HR] Ch. Hermite, Sur un nouveau d´eveloppement en s´erie de fonctions, C. R., t. 58 (1864), 93-100, 266-273.

### <span id="page-16-4"></span>GLOBAL (IN TIME) SOLUTIONS TO THE 3D-NAVIER-STOKES EQUATIONS ON R7

- [AK] P. Appell and J. Kampé de Fériet Fonctions hypergéométriques et hypersphériques, poloynomes d'Hermite. Gauthier-Villars, Paris, 1926.
- <span id="page-16-7"></span>[M] I. Miyadera, Nonlinear semigroups, Translations of Mathematical Monographs, Vol. 109, Amer. Math. Soc., Providence, RI, 1977.
- <span id="page-16-5"></span>[MOS] M. Magnus, F. Oberhettinger, and R. P. Soni, Formulas and Theorems for the Special Functions of Mathematical Physics, Band 52, Springer-Verlag, New York, 1966.
- <span id="page-16-2"></span><span id="page-16-1"></span>[RU] W. Rudin, Functional Analysis, McGraw-Hill, New York, 1973.
- [SY] G. R. Sell and Y. You, Dynamics of evolutionary equations, Applied Mathematical Sciences, Vol. 143, Springer, New York, 2002.
- <span id="page-16-0"></span> [S] E. M. Stein, Singular integrals and differentiability properties of functions, Princeton University Press, Princeton, NJ, 1970.
- <span id="page-16-3"></span>[T1] R. Temam, Navier-Stokes Equations, Theory and Numerical Analysis, AMS Chelsea Pub., Providence, RI, 2001.
- [T2] R. Temam, Infinite dimensional dynamical systems in mechanics and physics, Applied Mathematical Sciences, Vol. 68, Springer, New York, 1988.
- <span id="page-16-6"></span>[Z] E. Zeidler, Nonlinear functional analysis and its applications, Vol. IIB, Springer, New York, 1985.
- $(Tepper\ L.\ Gill)\ {\it DEPARTMENT}\ of\ {\it Electrical}\ Engineering,\ Howard\ University,\ Washington\ DC\ 20059,\ USA,\ E-mail:\ {\it tgill@howard.edu}$
- (Woodford W. Zachary) DEPARTMENT OF ELECTRICAL ENGINEERING, HOWARD UNIVERSITY, WASHINGTON DC 20059, USA, E-mail: www.achary@earthlink.net