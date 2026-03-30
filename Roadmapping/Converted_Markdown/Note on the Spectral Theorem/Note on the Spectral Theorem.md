## NOTE ON THE SPECTRAL THEOREM

#### T. L. GILL AND D. WILLIAMS

Abstract. In this note, we show that the spectral theorem, has two representations; the Stone-von Neumann representation and one based on the polar decomposition of linear operators, which we call the deformed representation. The deformed representation has the advantage that it provides an easy extension to all closed densely defined linear operators on Hilbert space. Furthermore, the deformed representation can also be extended all separable reflexive Banach spaces and has a limited extension to non-reflexive Banach spaces.

# Introduction

Let C[B] be the closed densely defined linear operators on a Banach space. By definition an operator A, defined on a separable Banach space B is of Baire class one if it can be approximated by a sequence {An} ⊂ L[B], of bounded linear operators. If B is a Hilbert space, then every A ∈ C[B] is of Baire class one. However, it turns out that, if B is not a Hilbert space, there may be operators A ∈ C[B] that are not of Baire class one.

<sup>1991</sup> Mathematics Subject Classification. Primary (46B03), (47D03) Secondary(47H06), (47F05) (35Q80).

Key words and phrases. spectral theorem, vector measures, vector-valued functions, Reflexive Banach spaces.

A Banach space B is said to be quasi-reflexive if dim{B′′/B} < ∞, and nonquasi-reflexive if dim{B′′/B} = ∞. Vinokurov, Petunin and Pliczko [\[VPP\]](#page-15-0) have shown that, for every nonquasi-reflexive Banach space B, there is a closed densely defined linear operator A which is not of Baire class one (for example, C[0, 1] or L 1 [R n ], n ∈ N). It can even be arranged so that A−<sup>1</sup> is a bounded linear injective operator (with dense range). This means that there does not exist a sequence of bounded linear operators A<sup>n</sup> ∈ L[B], with Anφ → Aφ for each A ∈ C[B] and φ ∈ D(A).

In [\[GSZ\]](#page-15-1), we were able to show that if B is one of the classical Banach spaces, then each bounded linear operator has an adjoint and a spectral type representation, which we called the extended representation. In addition, we were able to show that an operator has a adjoint if and only if it is of Baire class one. At that time, it was not clear if our results could be extended to the closed densely defined linear operators on any of the classical reflexive Banach spaces. However, we were able to show that every generator of a C0-semigroup on any classical space has an adjoint (is of Baire class one). We were later informed by Professor Pliczko that the results of their paper automatically implies that the set of closed densely defined linear operators on every separable reflexive (or quasi-reflexive) Banach is of Baire class one.

## Purpose

In the first section of this paper we show how the polar decomposition property of operators leads to a new type of (deformed) spectral representation for linear operators, which easily extends to all linear operators in C[H]. We then show directly that, for any separable reflexive Banach space, every closed densely defined linear operator has an adjoint. We use this result to prove that the closed densely defined linear operators are of Baire class one. Finally, we show that an almost identical (deformed) spectral representation holds if we replace H by a reflexive Banach space.

# Preliminaries

If H is a Hilbert space and A is any selfadjoint operator in C[H], the following direct spectral representation theorem is well-known (see [\[DS\]](#page-14-0), page 1192-99).

Theorem 0.1. Let A ∈ C[H] be a selfadjoint operator. Then its spectrum σ(A) ⊂ R and there exists a unique regular countably additive projectionvalued spectral measure E(·) mapping the Borel sets, B[R], over R into H such that

(1) D(A) satisfies

$$D(A) = \left\{ \phi \in \mathcal{H} \mid \int_{-\infty}^{\infty} \lambda^2 \left( \mathbf{E}(d\lambda)\phi, \phi \right)_{\mathcal{H}} < \infty \right\}$$

and

$$A\phi = \int_{-\infty}^{\infty} \lambda \mathbf{E}(d\lambda)\phi$$
, for each  $\phi \in D(A)$ .

(2) If g(·) is a complex-valued Borel function defined (a.e) on R, then g(A) ∈ C[H] and, for φ ∈ D(g(A)) = Dg(A),

$$g(A)\phi = \int_{-\infty}^{\infty} g(\lambda)\mathbf{E}(d\lambda)\phi,$$

where

$$D_{g}(A) = \left\{ \phi \in \mathcal{H} \mid \int_{-\infty}^{\infty} |g(\lambda)|^{2} \left( \mathbf{E}(d\lambda)\phi, \phi \right)_{\mathcal{H}} < \infty \right\}$$

and 
$$g(A^*) = \bar{g}(A)$$
.

Remark 0.2. We call Theorem 0.1 the direct representation because it requires the least additional material compared to the Gelfand representation and the one favored in mathematical physics, based on the position representation (see Rudin [\[RU\]](#page-15-2) page 306 and Reed and Simon [\[RS\]](#page-15-3), page 260).

A standard theorem of von Neumann [VN] shows that on H, any closed densely defined linear operator A has a well-defined adjoint A<sup>∗</sup> with A∗A nonnegative and selfadjoint. A classic result shows that there is a unique partial isometry U such that A = UT = TU¯ and A<sup>∗</sup> = U <sup>∗</sup>T¯ = T U<sup>∗</sup> , where T = [A∗A] 1/2 , T¯ = [AA<sup>∗</sup> 1/2 , and D(A) = D(T) (see Kato [\[K\]](#page-15-4), page 334). The next result can be found in Hille and Phillips [\[HP\]](#page-15-5) (see page 63).

Theorem 0.3. Let G(λ) be a vector-valued function from R to H of bounded variation. If h(λ) is a continuous complex-valued function on (a, b) ⊂ R, then the following holds:

- (1) The integral R <sup>b</sup> a h(λ)dG(λ) exists in the H norm.
- (2) If T ∈ C[H], G(λ) ∈ D(T) and TG(λ) is of bounded variation then

$$T \int_{a}^{b} h(\lambda) d\mathbf{G}(\lambda) = \int_{a}^{b} h(\lambda) dT\mathbf{G}(\lambda).$$

## 1. Hilbert Space Theory

Definition 1.1. If U is a partial isometry and E(·) is a positive spectral measure, then we call F(·) = UE(·) a deformed spectral measure.

Theorem 1.2. Let A ∈ C[H] be arbitrary. Then, for each φ ∈ D(A), there exists a deformed spectral measure F(·) such that:

(1) D(A) satisfies

$$D(A) = \left\{ \phi \in \mathcal{H} \mid \int_0^\infty \lambda^2 (d\mathbf{F}(\lambda)\phi, \phi)_{\mathcal{H}} < \infty \right\}$$

and

(2)

$$A\phi = \int_0^\infty \lambda d\mathbf{F}(\lambda)\phi$$
, for all  $\phi \in D(A)$ .

(3) If g(·) is a complex-valued Borel function defined (a.e) on R, then

$$g(A)\phi = \int_0^\infty g(\lambda)d\mathbf{F}(\lambda)\phi$$
 for all  $\phi \in D_g(A)$ ,

where

$$D_g(A) = \left\{ \phi \in \mathcal{H} \mid \int_0^\infty |g(\lambda)|^2 \left( d\mathbf{F}(\lambda)\phi, \phi \right)_{\mathcal{H}} < \infty \right\}.$$

*Proof.* To prove (1), write A = UT, where U is the unique partial isometry and  $T = [A^*A]^{1/2}$ . By Theorem 0.1, there is a positive spectral measure  $\mathbf{E}(\cdot)$  such that, for each  $x \in D(A) = D(T)$ :

(1.1) 
$$T\phi = \int_0^\infty \lambda d\mathbf{E}(\lambda)\phi.$$

Since  $\mathbf{E}(\lambda)\phi$  is a positive vector-valued function of bounded variation and U is a partial isometry,  $\mathbf{F}(\lambda)x = U\mathbf{E}(\lambda)\phi$  is of bounded variation, with  $Var(\mathbf{F}\phi,\mathbb{R}) \leq Var(\mathbf{E}\phi,\mathbb{R})$ . Thus, by Theorem 0.2,

$$U \int_0^\infty \lambda d\mathbf{E}(\lambda) \phi = \int_0^\infty \lambda dU \mathbf{E}(\lambda) \phi.$$

Since  $A\phi = UT\phi$ , if we set  $\mathbf{F}(\lambda)x = U\mathbf{E}(\lambda)x$ , we have from equation (1.1),

(1.2) 
$$A\phi = \int_0^\infty \lambda d\mathbf{F}(\lambda)\phi.$$

The proof of (2) and (3) follows the proof of the same results in [DS], when we recall that  $||A\phi||^2 = ||T\phi||^2$ , so that

$$\begin{split} &\left\{\phi \, \bigg| \, \int_0^\infty \lambda^2 \, (d\mathbf{E}(\lambda)\phi,\phi) < \infty \, \right\} = \left\{\phi \, \bigg| \, \int_0^\infty \lambda^2 \, (d\mathbf{F}(\lambda)\phi,\phi) < \infty \, \right\} \\ &\left\{\phi \, \bigg| \, \int_0^\infty |g(\lambda)|^2 \, (d\mathbf{E}(\lambda)\phi,\phi) < \infty \, \right\} = \left\{\phi \, \bigg| \, \int_0^\infty |g(\lambda)|^2 \, (d\mathbf{F}(\lambda)\phi,\phi) < \infty \, \right\}. \end{split}$$

Remark 1.3. Let A be self adjoint, with its spectrum on the negative real axis. In this case, the standard spectral theorem gives us:

(1.3) 
$$A = \int_{-r_A}^{0} \lambda d\mathbf{E}(\lambda).$$

However, the deformed spectral theorem gives

(1.4) 
$$A = \int_0^{r_A} \lambda d\mathbf{F}(\lambda).$$

In fact, the domain of integration for  $\mathbf{F}(\lambda)$  always coincides with the spectrum of the positive linear operator T, where A = UT. Thus, if A is not a positive selfadjoint linear operator, the two representations are distinct.

#### 2. Adjoints on Reflexive Banach Spaces

In this section, we define the adjoint for operators in  $\mathcal{C}[\mathcal{B}]$ , but first, we collect a few tools, which are needed for the sequel. Recall that duality map  $J: \mathcal{B} \mapsto \mathcal{B}'$ , is the set

$$J(u) = \left\{ F_u \in \mathcal{B}' \, \middle| F_u(u) = \langle u, F_u \rangle = ||u||^2 = ||F_u||^2 \right\}, \ \forall u \in \mathcal{B}.$$

We want to construct a special duality map. Assume that  $\mathcal{B} \subset \mathcal{H}$  is a dense continuous embedding and fix  $u \in \mathcal{B}$ . Let  $M = \langle u \rangle$  be the closed subspace spanned by u and define a seminorm  $p_u(\cdot)$  on  $\mathcal{B}$  by  $p_u(v) = ||u||_{\mathcal{B}} ||v||_{\mathcal{B}}$ . Define the map  $\hat{S}_u(\cdot) = \langle \cdot, \hat{S}_u \rangle$  by:

$$\langle v, \hat{S}_u \rangle = \hat{S}_u(v) = \frac{\|u\|_{\mathcal{B}}^2}{\|u\|_2^2} (v, u)_2.$$

On the closed subspace  $M = \langle u \rangle$ ,  $|\hat{S}_u(v)| = ||u||_B ||v||_B \leqslant p_u(v)$ . By the Hahn-Banach Theorem,  $\hat{S}_u(\cdot)$  has an extension,  $S_u(\cdot)$ , to  $\mathcal{B}$  such that  $|S_u(v)| \leqslant p_u(v) = ||u||_B ||v||_B$  for all  $v \in \mathcal{B}$ . From here, we see that  $||S_u||_{\mathcal{B}'} \leq ||u||_{\mathcal{B}}$ . On the other hand, we have  $||u||_{\mathcal{B}}^2 = S_u(u) \leqslant ||u||_{\mathcal{B}} ||S_u||_{\mathcal{B}'}$ , so that  $S_u(\cdot)$  is a duality mapping for u. We call  $S_u(\cdot)$  the Steadman duality map on  $\mathcal{B}$  associated with  $\mathcal{H}$ .

Recall that on  $\mathcal{B}$ , a densely defined operator A is called accretive if  $\operatorname{Re}\langle Au, F_u \rangle \geq 0$  for  $u \in D(A)$  and any duality map  $F_u$ . The following definition extracts the essential properties of an adjoint operator on Hilbert space, but also makes sense on a Banach space

**Definition 2.1.** If  $A \in C[B]$ , the closed densely defined linear operators on B, we say that  $A^*$  is a adjoint of A if:

- (1) the operator  $A^*A \ge 0$  (accretive),
- (2)  $(A^*A)^* = A^*A$  (naturally selfadjoint), and
- (3)  $I + A^*A$  has a bounded inverse.

We need the following result by Lax [L].

**Theorem 2.2** (Lax). Suppose  $\mathcal{B}$  is a dense continuous embedding in a separable Hilbert space  $\mathcal{H}$ . Let  $A \in L[\mathcal{B}]$ . If A is selfadjoint on  $\mathcal{H}$  (i.e.,  $(Au, v)_{\mathcal{H}} = (u, Av)_{\mathcal{H}}, \forall u, v \in \mathcal{B}$ ), then A is bounded on  $\mathcal{H}$  and  $\|A\|_{\mathcal{H}} \leq k \|A\|_{\mathcal{B}}$  for some positive constant k.

The following general result is a variation of one due to Kuelbs [KB].

**Theorem 2.3.** Suppose  $\mathcal{B}$  is a separable reflexive Banach space, then there exist a separable Hilbert space  $\mathcal{H}$  such that,

- (1)  $\mathcal{B} \subset \mathcal{H}$  as a continuous dense embedding.
- (2)  $\mathcal{B}' \subset \mathcal{H}'$  as a continuous dense embedding.

*Proof.* To prove (1), let  $\{u_n\}$  be a dense set in  $\mathcal{B}$  and let  $\{f_n\}$  be any fixed set of corresponding duality mappings (i.e.,  $f_n \in \mathcal{B}'$ , the dual space of  $\mathcal{B}$  and  $f_n(u_n) = \langle u_n, f_n \rangle = \|u_n\|_{\mathcal{B}}^2 = \|f_n\|_{\mathcal{B}'}^2$ ). Let  $\{t_n\}$  be a positive sequence of numbers such that  $\sum_{n=1}^{\infty} t_n = 1$ , and define (u, v) by:

$$(u,v) = \sum_{n=1}^{\infty} t_n f_n(u) \bar{f}_n(v).$$

It is easy to see that (u, v) is an inner product on  $\mathcal{B}$ . Let  $\mathcal{H}$  be the completion of  $\mathcal{B}$  with respect to this inner product. It is clear that  $\mathcal{B}$  is dense in  $\mathcal{H}$ , and

$$||u||^2 = \sum_{n=1}^{\infty} t_n |f_n(u)|^2 \le \sup_n |f_n(u)|^2 = ||u||_{\mathcal{B}}^2,$$

so the embedding is continuous.

For (2), we note that, since  $\mathcal{B}$  is reflexive  $\mathcal{B} = \mathcal{B}''$ . In this case the set  $\{f_n\}$  is dense in  $\mathcal{B}'$ , so we may use the dense family  $\{u_n\} \subset \mathcal{B}$  to define an inner product on  $\mathcal{B}'$  by

$$(f,g) = \sum_{n=1}^{\infty} t_n u_n''(f) \bar{u}_n''(g) = \sum_{n=1}^{\infty} t_n f(u_n) \bar{g}(u_n),$$

where  $u_n''(h) = h(u_n)$ , for each  $h \in \mathcal{B}'$ . The completion of  $\mathcal{B}'$  with the above inner product provides a construction of  $\mathcal{H}'$ . It is clear that  $\mathcal{B}' \subset \mathcal{H}'$  as a continuous dense embedding.

**Theorem 2.4.** When  $\mathcal{B}$  is reflexive every operator  $A \in \mathcal{C}[\mathcal{B}]$  has a well defined adjoint  $A^* \in \mathcal{C}[\mathcal{B}]$ . If  $A \in L[\mathcal{B}]$ , the bounded linear operators on  $\mathcal{B}$ , then  $A^* \in L[\mathcal{B}]$ .

*Proof.* Let  $\mathbf{J}: \mathcal{H} \to \mathcal{H}'$ . It is easy to see that  $\mathbf{J}^* = \mathbf{J}$ . Since  $[\mathbf{J}]_{\mathcal{B}}$  is one to one and onto  $\mathcal{B}'$ , if  $A \in \mathcal{C}[\mathcal{B}]$ , then  $[A'\mathbf{J}]_{\mathcal{B}}: \mathcal{B}' \to \mathcal{B}'$ . Since A' is closed and densely defined, it follows that  $\mathbf{J}^{-1}A'\mathbf{J}: \mathcal{B} \to \mathcal{B}$  is closed and densely defined. Thus, we can define  $A^* = [\mathbf{J}^{-1}A'\mathbf{J}]_{\mathcal{B}}$ .

In case  $A \in L[\mathcal{B}]$ , we know that  $A^* = [\mathbf{J}^{-1}A'\mathbf{J}]_{\mathcal{B}}$  is defined on all of  $\mathcal{B}$ . By the closed graph theorem,  $A^* \in L[\mathcal{B}]$ .

In both cases,  $A^*$  and  $[A^*]'$  are also densely defined on  $\mathcal{H}$  and  $\mathcal{H}'$  respectively, so that we can extend A and  $A^*$  to closed densely defined linear operators on  $\mathcal{H}$ . Furthermore, the injective nature of  $\mathbf{J}$  and that of  $\mathcal{B} \to \mathcal{B}'$  means that  $A^*A$  extends to an accretive self adjoint linear operator on  $\mathcal{H}$ .

If  $A \in L[\mathcal{B}]$ , then  $(A^*Au, u) \ge 0$ , and

$$(A^*Ag, u) = \langle A^*Au, \mathbf{J}(u) \rangle = \langle \mathbf{J}^{-1}A'\mathbf{J}(Au), \mathbf{J}(u) \rangle$$
  
=  $\langle \mathbf{J}(Au), Au \rangle = \langle u, A'\mathbf{J}(Au) \rangle = \langle \mathbf{J}(u), A^*Au \rangle = \langle u, A^*Au \rangle$ 

for all  $u \in \mathcal{B}$ . By Lax's Theorem,  $A^*A$  has a bounded extension to  $\mathcal{H}$  and  $\|A^*A\| = \|A\|^2 \leqslant k \|A\|_B^2$ , where k is a positive constant.

The above theorem shows that every closed densely defined linear operator on  $\mathcal{B}$  has a closed densely defined extension to  $\mathcal{H}$ , with bounded operators on  $\mathcal{B}$  becoming bounded on  $\mathcal{H}$ . The proof of the next result follows.

Corollary 2.5. If  $\mathcal{B}$  be a separable reflexive Banach space and  $A \in \mathcal{C}[\mathcal{B}]$ , then

- (1) the operator  $A^*A \ge 0$  (accretive),
- (2)  $(A^*A)^* = A^*A$  (naturally selfadjoint), and
- (3)  $I + A^*A$  has a bounded inverse.

**Theorem 2.6.** If  $\mathcal{B}$  be a separable reflexive Banach space, then every  $A \in \mathcal{C}[\mathcal{B}]$  is of Baire class one.

*Proof.* Since each  $A \in \mathcal{C}[\mathcal{B}]$  has a extension  $\bar{A}$  to  $\mathcal{C}[\mathcal{H}]$ , we can write A = UT, where UT is the restriction to  $\mathcal{B}$  of the polar decomposition of  $\bar{A}$ . From  $A^* = TU^*$ , we see that  $T = A^*U$ , so that

$$AA^* = (UT)(TU^*) = UA^*AU^* \Rightarrow AA^*U = UA^*A.$$

It follows that,  $A = \bar{T}U$  and  $A^* = U^*\bar{T}$ .

For  $\lambda > 0$ , let  $R(\lambda, -T)$  be the resolvent of -T and let

(2.1) 
$$A_{\lambda} = \lambda AR(\lambda, -T) = \lambda^2 UR(\lambda, -T) - \lambda U.$$

It is easy to see that  $A_{\lambda}$  is bounded and that  $AR(\lambda, -T)\phi = R(\lambda, -\overline{T})A\phi$  for  $\phi \in D(A)$ .

If 0 < λ and φ ∈ D(A) we have

$$R(\lambda, -T) (\lambda I + T) \phi = \phi \implies$$

$$\lambda R(\lambda, -T) \phi - \phi = R(\lambda, -T) (-T\phi) \implies$$

$$\|\lambda R(\lambda, -T) \phi - \phi\| \le \|R(\lambda, -T)\| \|T\phi\| \le \lambda^{-1} \|T\phi\|$$

This last term converges to zero as λ → ∞, so that

$$\lim_{\lambda \to \infty} \lambda R(\lambda, -T)\phi = \phi.$$

Since D(A) is dense, the convergence holds for all φ ∈ B.

For the second part, we see from the last result that

$$\lim_{\lambda \to \infty} \lambda AR(\lambda, -T)\phi = \lim_{\lambda \to \infty} \lambda R(\lambda, -\bar{T})A\phi = A\phi,$$

whenever φ ∈ D(A).

Examples. We begin with the following useful result (see Kato [\[K\]](#page-15-4), pg. 168).

Theorem 2.7. Let T be a densely defined linear operator on a reflexive Banach space B. Then, the following holds:

- (1) The adjoint of T, T<sup>∗</sup> is a closed linear operator.
- (2) The operator T has a closed extension if and only if D(T ∗ ) is dense in B ∗ . In this case, the closure T¯ = T ∗∗ .
- (3) If T is closable, then (T¯) <sup>∗</sup> = T ∗ .

From the last result, we see that, for any closed densely defined linear operator A defined on L p [R n ], 1 < p < ∞, for which the domain of A∗ , D(A<sup>∗</sup> ) ⊂ L q [R n ], 1/p + 1/q = 1, is dense in L p [R n ], also has a closed densely defined extension to L p [R n ].

Example 2.8. Let A be a second order differential operator on L p [R n ], of the form

$$A = \sum_{i,j=1}^{n} a_{ij}(\mathbf{x}) \frac{\partial^2}{\partial x_i \partial x_j} + \sum_{i,j=1}^{n} b_{ij}(\mathbf{x}) x_j \frac{\partial}{\partial x_i},$$

where a(x) = [[aij (x)]] and b(x) = [[bij (x)]] are matrix-valued functions in C<sup>∞</sup> c [R <sup>n</sup> × R n ] (infinitely differentiable functions with compact support). We also assume that, for all x ∈ R <sup>n</sup> det[[aij (x)]] > ε and the imaginary part of the eigenvalues of b(x) are bounded above by −ε, for some ε > 0. Note, since we don't require a or b to be symmetric, A 6= A<sup>∗</sup> .

It is well-known that C<sup>∞</sup> c [R n ] is dense in L p [R n ] ∩ L q [R n ] for all p, q ∈ [1,∞) ∩ N. Furthermore, since A<sup>∗</sup> is invariant on C<sup>∞</sup> c [R n ], A<sup>∗</sup> : L p [R n ] → L p [R n ]. It now follows from Theorem 2.8, that A<sup>∗</sup> has a closed densely defined extension to L p [R n ].

Example 2.9. The second example shows directly the setup used in Theorem 2.3. Let Ω be a bounded open domain of class C 1 in R <sup>n</sup> and let H<sup>1</sup> 0 [Ω], the set of all real-valued functions u ∈ L 2 [Ω] such that their first order weak partial derivatives are in L 2 [Ω] and vanish on the boundary. It follows that

$$(u, v) = \int_{\Omega} \nabla u(\mathbf{x}) \cdot \nabla v(\mathbf{x}) d\mathbf{x}$$

defines an inner product on  $\mathcal{H}_0^1[\Omega]$ . The dual  $\mathcal{H}^{-1}[\Omega]$  coincides with the set of all distributions of the form

$$u = h_0 + \sum_{i=1}^n \frac{\partial h_i}{\partial x_i}, \quad \text{where } h_i \in L^2[\Omega], \quad 1 \leqslant i \leqslant n.$$

In this case we also have for  $p \in [2, \infty)$  and  $q \in (1, 2], \frac{1}{p} + \frac{1}{q} = 1$  that,

$$\mathcal{H}_0^1[\Omega] \subset L^p[\Omega] \subset L^q[\Omega] \subset \mathcal{H}^{-1}[\Omega]$$

all as continuous dense embeddings.

From the inner product on  $\mathcal{H}_0^1[\Omega]$  we see that  $J_0 = -\Delta$ , the Laplace operator under Dirichlet homogeneous boundary conditions on  $\Omega$ . If we set  $\mathcal{H} = \mathcal{H}^{-1}$  and  $J = J_0^{-1}$ , we can apply Theorem 2.4 to obtain  $A^* \in \mathcal{C}[L^r(\Omega)]$ , for all  $1 < r < \infty$ , by  $A^* = -\Delta A'[-\Delta]^{-1}$ , for all  $A' \in \mathcal{C}[L^r'(\Omega)]$ . Where  $\frac{1}{r} + \frac{1}{r'} = 1$ .

### 2.1. Spectral Theorem.

**Theorem 2.10.** Let  $A \in \mathcal{C}[\mathcal{B}]$  be arbitrary, where  $\mathcal{B}$  is a separable reflexive Banach space. Then, for each  $\phi \in D(A)$ , there exists a deformed spectral measure  $\mathbf{F}(\cdot)$  and a vector-valued function  $\mathbf{F}(\lambda)\phi$  of bounded variation such that:

(1) D(A) satisfies

$$D(A) = \left\{ \phi \in \mathcal{B} \mid \int_{0}^{\infty} \lambda^{2} \left( d\mathbf{F}(\lambda)\phi, S_{\phi} \right)_{\mathcal{B}} < \infty \right\}$$

and

(2)

$$A\phi = \int_0^\infty \lambda d\mathbf{F}(\lambda)\phi$$
, for all  $\phi \in D(A)$ .

(3) If g(·) is a complex-valued Borel function defined (a.e) on R, then

$$g(A)\phi = \int_0^\infty g(\lambda)d\mathbf{F}(\lambda)\phi$$
 for all  $\phi \in D_g(A)$ ,

where

$$D_g(A) = \left\{ \phi \in \mathcal{B} \mid \int_0^\infty |g(\lambda)|^2 \left( d\mathbf{F}(\lambda)\phi, S_\phi \right)_{\mathcal{B}} < \infty \right\}.$$

Conclusion. A major result of this paper is the discovery of a new direct spectral type representation for the family of closed densely defined linear operators on a Hilbert space, which we call the deformed representation. The advantage of this approach is that it has a similar extension to the family of closed densely defined linear operators on a separable reflexive Banach space.

Acknowledgements. During the course of the development of this work, we have benefited from important critical remarks from Professor Ioan I. Vrabie. We would like to sincerely thank Professor Anatolij Pliczko for important correspondence on operators of Baire class on Banach spaces.

## References

<span id="page-14-0"></span>[DS] N. Dunford and J. T. Schwartz, Linear Operators Part II: Spectral Theory, Wiley Classics edition, Wiley Interscience (1988).

- <span id="page-15-1"></span>[GSZ] T. L. Gill V. Steadman and W. W. Zachary, The Polar Decomposition in Banach Spaces, African Diaspora J. Math. 11 (2011) 98-131.
- <span id="page-15-5"></span>[HP] E. Hille and R. S. Phillips, Functional Analysis and Semigroups, Amer. Math. Soc. Colloq. Pub. 31, Amer. Math. Soc. Providence, RI, (1957).
- <span id="page-15-4"></span>[K] T. Kato, Perturbation Theory for Linear Operators, second ed. Springer-Verlag, New York, (1976).
- <span id="page-15-7"></span>[KB] J. Kuelbs, Gaussian measures on a Banach space, Journal of Functional Analysis 5 (1970), 354–367.
- [L] P. D. Lax, Symmetrizable linear tranformations, Comm. Pure Appl. Math. 7 (1954), 633-647.
- <span id="page-15-6"></span>[L] P. D. Lax, Symmetrizable linear tranformations. Comm. Pure Appl. Math. 7 (1954), 633–647.
- <span id="page-15-3"></span>[RS] M. Reed and B. Simon, Methods of Modern Mathematical Physics I: Functional Analysis, Academic Press, New York, (1972).
- <span id="page-15-2"></span>[RU] W. Rudin, Functional Analysis, McGraw-Hill Press, New York, (1973).
- [SO] M. H. Stone, Linear Transformations in Hilbert Space , Math. Surveys 15, Amer. Math. Soc. Colloq. Publ. 15, Providence, RI, (1932).
- [VN] J. von Neumann, Uber adjungierte Funktionaloperatoren, ¨ Annals of Mathematics 33 (1932), 294-310.
- <span id="page-15-0"></span>[VPP] V. A. Vinokurov, Yu. Petunin and A. N. Pliczko, Measurability and Regularizability mappings inverse to continuous linear operators (in Russian), Mat. Zametki. 26 (1979), no. 4, 583-591. English translation: Math. Notes 26 (1980), 781-785.

(Tepper L. Gill) Departments of Mathematics, Physics, and Electrical & Computer Engineering, Howard University, Washington DC 20059, USA, Email : tgill@howard.edu

(Daniel Williams) Department of Mathematics, Howard University, Washington DC 20059, USA, E-mail : dwilliams@howard.edu