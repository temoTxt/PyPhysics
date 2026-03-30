## ADJOINT OPERATORS ON BANACH SPACES

#### T. L. GILL, F. MENSAH, AND W. W. ZACHARY

Abstract. In this paper, we report on new results related to the existence of an adjoint for operators on separable Banach spaces and discuss a few interesting applications. (Some results are new even for Hilbert spaces.) Our first two applications provide an extension of the Poincar´e inequality and the Stone-von Neumann version of the spectral theorem for a large class of C0-generators of contraction semigroups on separable Banach spaces. Our third application provides a natural extension of the Schatten-class of operators to all separable Banach spaces. As a part of this program, we introduce a new class of separable Banach spaces. As a side benefit, these spaces also provide a natural framework for the (rigorous) construction of the path integral as envisioned by Feynman.

## Contents

| 1.                                                                                                                      | Adjoint Theory                                                           | 3   |  |  |
|-------------------------------------------------------------------------------------------------------------------------|--------------------------------------------------------------------------|-----|--|--|
| 1.1.                                                                                                                    | Introduction                                                             | 3   |  |  |
| 1.2.                                                                                                                    | Background                                                               | 3   |  |  |
| 1.3.                                                                                                                    | Banach Space Adjoint                                                     | 7   |  |  |
| 1991                                                                                                                    | Mathematics<br>Subject<br>Classification. Primary<br>(46B03),<br>(47D03) | Sec |  |  |
| ondary(47H06), (47F05) (35Q80).<br>Key words and phrases. Poincar´e inequality spectral theorem, semigroups, vector mea |                                                                          |     |  |  |
|                                                                                                                         | sures, vector-valued functions, Schatten-class.                          |     |  |  |

| 2.   | The Hilbert Space KS2                 | 10 |
|------|---------------------------------------|----|
| 2.1. | Weak Integral                         | 14 |
| 2.2. | Discussion                            | 16 |
| 2.3. | The Corresponding<br>H1               | 16 |
| 3.   | Closed Operators                      | 17 |
| 3.1. | Operators on<br>B                     | 23 |
| 4.   | Extension of the Poincar´e inequality | 25 |
| 4.1. | Introduction                          | 25 |
| 4.2. | Purpose                               | 26 |
| 4.3. | Hilbert space case                    | 26 |
| 4.4. | Banach space case                     | 28 |
| 5.   | Extension Of The Spectral Theorem     | 29 |
| 5.1. | Introduction                          | 29 |
| 5.2. | Background                            | 31 |
| 5.3. | Problem                               | 32 |
| 5.4. | Hilbert Space case                    | 33 |
| 5.5. | Banach space case                     | 36 |
| 5.6. | General Case                          | 37 |
| 6.   | Schatten Classes                      | 39 |
| 6.1. | Discussion                            | 42 |
| 7.   | Conclusion                            | 44 |

<span id="page-2-0"></span>References 45

## 1. Adjoint Theory

- <span id="page-2-1"></span>1.1. Introduction. One of the impediments to the development of a clear parallel theory for operators on Banach spaces compared to that for Hilbert spaces is the lack of a suitable notion of an adjoint operator. In this section we use a Theorem of Gross and Kuelbs to construct an adjoint for all bounded linear operators on a separable Banach space. We then show that this result can be extended to all closed densely defined linear operators of Baire class one (limits of bounded linear operators). We use these results in later sections to extend the Poincaré inequality, the spectral theorem and to construct the "natural" Banach space version of the Schatten class.
- <span id="page-2-2"></span>1.2. **Background.** Let  $\mathcal{B}$  be a separable Banach space over the complex field and let  $L[\mathcal{B}]$  denote the bounded linear operators on  $\mathcal{B}$ . Assume that  $\mathcal{B}$  has a continuous dense embedding in a Hilbert space  $\mathcal{H}$ . By a duality map,  $f_u: \mathcal{B} \to \mathcal{B}'$ , we mean any linear functional  $f_u \in \{f \in \mathcal{B}' | f(u) = \langle u, f \rangle = \|u\|_{\mathcal{B}}^2 = \|f\|_{\mathcal{B}'}^2$ ,  $u \in \mathcal{B}\}$ , where  $\langle . \rangle$  is the natural pairing between a Banach space and its dual. Let  $\mathbf{J}: \mathcal{H} \longrightarrow \mathcal{H}'$  be the standard conjugate isomorphism between a Hilbert space and its dual, so that  $\langle u, \mathbf{J}(u) \rangle = \langle u, u \rangle_{\mathcal{H}} = \|u\|_{\mathcal{H}}^2$ .

For fixed u define a seminorm  $p_u(\cdot)$  on  $\mathcal{B}$  by  $p_u(x) = ||u||_{\mathcal{B}} ||x||_{\mathcal{B}}$ , and define  $\hat{f}_u^s(\cdot)$  by:

$$\hat{f}_u^s(x) = \frac{\|u\|_{\mathcal{B}}^2}{\|u\|_{\mathcal{H}}^2} (x, u)_{\mathcal{H}}.$$

On the closed subspace  $M = \langle u \rangle$ ,  $|\hat{f}_u^s(x)| = ||u||_B ||x||_B \leqslant p_u(x)$ . By the (complex version of the) Hahn-Banach Theorem,  $\hat{f}_u^s(\cdot)$  has an extension,  $f_u^s(\cdot)$ , to  $\mathcal{B}$  such that  $|f_u^s(x)| \leqslant p_u(x) = ||u||_B ||x||_B$  for all  $x \in \mathcal{B}$  (see Rudin [RU], Theorem 3.3, page 57). From here, we see that  $||f_u^s||_{\mathcal{B}'} \leq ||u||_{\mathcal{B}}$ .

On the other hand, we have:

$$f_u^s(u) = \|u\|_{\mathcal{B}}^2 \leqslant \|u\|_{\mathcal{B}} \|f_u^s\|_{\mathcal{B}'} \Rightarrow \|u\|_{\mathcal{B}} \leqslant \|f_u^s\|_{\mathcal{B}'},$$

so that  $f_u^s(\cdot)$  is a duality mapping for u. We call  $f_u^s(\cdot)$  the Steadman duality map on  $\mathcal{B}$  associated with  $\mathcal{H}$ .

Recall that a densely defined operator A is called accretive if  $\operatorname{Re}\langle Au, f_u^s \rangle \geq 0$  for  $u \in D(A)$ ; and it is called m-accretive if, in addition, it is closed and  $\operatorname{Ran}(I+A) = \mathcal{B}$ . The following theorem by Lax [L] is important for our theory. It is not as well-known as it should be, so we provide a proof of the first part. We prove a stronger version of parts two and three in Section 3 (see Theorem 3.3, part 2).

**Theorem 1.1** (Lax). Suppose  $\mathcal{B}$  is a dense continuous embedding in a separable Hilbert space  $\mathcal{H}$ . Let  $A \in L[\mathcal{B}]$ . If A is selfadjoint on  $\mathcal{H}$  (i.e.,  $(Ax,y)_{\mathcal{H}} = (x,Ay)_{\mathcal{H}}, \forall x,y \in \mathcal{B}$ ), then

- (1) The operator A is bounded on  $\mathcal{H}$  and  $||A||_{\mathcal{H}} \leq k ||A||_{\mathcal{B}}$ , for some positive constant k.
- (2) The spectra of A over  $\mathcal{H}$  and over  $\mathcal{B}$  satisfy  $\sigma_{\mathcal{H}}(A) \subset \sigma_{\mathcal{B}}(A)$ .
- (3) The point spectrum of A is unchanged by the extension (i.e.,  $\sigma^p_{\mathcal{H}}(A) = \sigma^p_{\mathcal{B}}(A)$ ).

*Proof.* To prove (1), let  $x \in \mathcal{B}$  and, without loss, we can assume that k = 1 and  $||x||_{\mathcal{H}} = 1$ . Since A is selfadjoint,

$$||Ax||_{\mathcal{H}}^2 = (Ax, Ax) = (x, A^2x) \leqslant ||x||_{\mathcal{H}} ||A^2x||_{\mathcal{H}} = ||A^2x||_{\mathcal{H}}.$$

Thus, we have  $||Ax||_{\mathcal{H}}^4 \leq ||A^4x||_{\mathcal{H}}$ , so it is easy to see that  $||Ax||_{\mathcal{H}}^{2n} \leq ||A^{2n}x||_{\mathcal{H}}$  for all n. It follows that:

$$||Ax||_{\mathcal{H}} \leqslant (||A^{2n}x||_{\mathcal{H}})^{1/2n} \leqslant (||A^{2n}x||_{\mathcal{B}})^{1/2n}$$
  
$$\leqslant (||A^{2n}||_{\mathcal{B}})^{1/2n} (||x||_{\mathcal{B}})^{1/2n} \leqslant ||A||_{\mathcal{B}} (||x||_{\mathcal{B}})^{1/2n}.$$

Letting  $n \to \infty$ , we get that  $||Ax||_{\mathcal{H}} \leq ||A||_{\mathcal{B}}$  for x in a dense set of the unit ball of  $\mathcal{H}$ . We are done, since the norm is attained on a dense set of the unit ball.

The following is a result due to Gross and Kuelbs [GR], [KB].

**Theorem 1.2.** Suppose  $\mathcal{B}$  is a separable Banach space. Then there exist separable Hilbert spaces  $\mathcal{H}_1, \mathcal{H}_2$  and a positive trace class operator  $\mathbf{T}_{12}$  defined on  $\mathcal{H}_2$  such that  $\mathcal{H}_1 \subset \mathcal{B} \subset \mathcal{H}_2$  (all as continuous dense embeddings).

*Proof.* As  $\mathcal{B}$  is separable, let  $\{u_n\}$  be a dense set in  $\mathcal{B}$  and let  $\{f_n\}$  be any fixed set of corresponding duality mappings (i.e.,  $f_n \in \mathcal{B}'$  and  $f_n(u_n) = \langle u_n, f_n \rangle = \|u_n\|_{\mathcal{B}}^2 = \|f_n\|_{\mathcal{B}'}^2$ ). Let  $\{t_n\}$  be a positive sequence of numbers such that  $\sum_{n=1}^{\infty} t_n = 1$ , and define  $(u, v)_2$  by:

$$(u,v)_2 = \sum_{n=1}^{\infty} t_n f_n(u) \bar{f}_n(v).$$

It is easy to see that  $(u, v)_2$  is an inner product on  $\mathcal{B}$ . We let  $\mathcal{H}_2$  be the Hilbert space generated by the completion of  $\mathcal{B}$  with respect to this inner product. It is clear that  $\mathcal{B}$  is dense in  $\mathcal{H}_2$ , and as

$$||u||_{2}^{2} = \sum_{n=1}^{\infty} t_{n} |f_{n}(u)|^{2} \le \sup_{n} |f_{n}(u)|^{2} = ||u||_{\mathcal{B}}^{2},$$

we see that the embedding is continuous.

. Now, let  $\{\varphi_n\} \in \mathcal{B}$  be a complete orthonormal sequence for  $\mathcal{H}_2$ , and let  $\{\lambda_n\}$  be a positive sequence such that  $\sum_{n=1}^{\infty} \lambda_n < \infty$ , and  $M = \sum_{n=1}^{\infty} \lambda_n^2 \|\varphi_n\|_{\mathcal{B}}^2 < \infty$ . Define the operator  $\mathbf{T}_{12}$  on  $\mathcal{B}$  by:

$$\mathbf{T}_{12}u = \sum_{n=1}^{\infty} \lambda_n \left( u, \varphi_n \right)_2 \varphi_n.$$

Since  $\mathcal{B} \subset \mathcal{H}_2$ ,  $(u, \varphi_n)_2$  is defined for all  $u \in \mathcal{B}$ . Thus,  $\mathbf{T}_{12}$  maps  $\mathcal{B} \to \mathcal{B}$  and:

$$\|\mathbf{T}_{12}u\|_{\mathcal{B}}^{2} \leq \left[\sum_{n=1}^{\infty} \lambda_{n}^{2} \|\varphi_{n}\|_{\mathcal{B}}^{2}\right] \left[\sum_{n=1}^{\infty} |(u,\varphi_{n})_{2}|^{2}\right] = M \|u\|_{2}^{2} \leq M \|u\|_{\mathcal{B}}^{2}.$$

Thus,  $T_{12}$  is a bounded operator on  $\mathcal{B}$ . Define  $\mathcal{H}_1$  by:

$$\mathcal{H}_{1} = \left\{ u \in \mathcal{B} \mid \sum_{n=1}^{\infty} \lambda_{n}^{-1} |(u, \varphi_{n})_{2}|^{2} < \infty \right\}, \quad (u, v)_{1} = \sum_{n=1}^{\infty} \lambda_{n}^{-1} (u, \varphi_{n})_{2} (\varphi_{n}, v)_{2}.$$

With the above inner product, H<sup>1</sup> is a Hilbert space and, since terms of the form {u<sup>N</sup> = P<sup>N</sup> <sup>k</sup>=1 λ −1 k (u, ψk) <sup>2</sup> ϕ<sup>k</sup> : u, ψ<sup>k</sup> ∈ B} are dense in B, we see that H<sup>1</sup> is dense in B. It follows that H<sup>1</sup> is also dense in H2. It is easy to see that T<sup>12</sup> is a positive selfadjoint operator with respect to the H<sup>2</sup> inner product so, by Theorem 1.1, T<sup>12</sup> has a bounded extension to H<sup>2</sup> and kT12k<sup>2</sup> ≤ kT12k<sup>B</sup> . Finally, it is easy to see that, for u, v ∈ H1, (u, v)<sup>1</sup> = (T −1/2 <sup>12</sup> u, T −1/2 <sup>12</sup> v)<sup>2</sup> and (u, v)<sup>2</sup> = (T 1/2 <sup>12</sup> u, T 1/2 <sup>12</sup> v)1. It follows that H<sup>1</sup> is continuously embedded in H2, hence also in B.

The construction of H<sup>1</sup> and H<sup>2</sup> is not unique. In the next section, we construct a concrete version of H<sup>2</sup> which is unique in the sense that we use a fixed dense family {un} ⊂ B, a fixed family of linear functionals {Fn} ∈ B′ and a fixed family of positive numbers {tn}. (We will discuss this more in the remarks before Section 2.1.) For the remainder of this paper, we assume that both H<sup>1</sup> and H<sup>2</sup> are fixed.

<span id="page-6-0"></span>1.3. Banach Space Adjoint. The following is the major result in Gill et al [GBZS]. It generalizes the well-known result of von Neumann [VN] for bounded operators on Hilbert spaces. For convenience, we provide a proof. (We delay the proof of (1) and (3) until after Theorem 1.4.)

Theorem 1.3. Let A be a bounded linear operator on B. Then A has a well-defined adjoint A<sup>∗</sup> defined on B such that:

(1) the operator A∗A ≥ 0 (maximal accretive),

- (2)  $(A^*A)^* = A^*A$  (selfadjoint), and
- (3)  $I + A^*A$  has a bounded inverse.

*Proof.* If we let  $\mathbf{J}_i: \mathcal{H}_i \to \mathcal{H}'_i$ , (i = 1, 2), then  $A_1 = A_{|\mathcal{H}_1}: \mathcal{H}_1 \to \mathcal{H}_2$ , and  $A'_1: \mathcal{H}'_2 \to \mathcal{H}'_1$ .

It follows that  $A'_1\mathbf{J}_2: \mathcal{H}_2 \to \mathcal{H}'_1$  and  $\mathbf{J}_1^{-1}A'_1\mathbf{J}_2: \mathcal{H}_2 \to \mathcal{H}_1 \subset \mathcal{B}$  so that, if we define  $A^* = [\mathbf{J}_1^{-1}A'_1\mathbf{J}_2]_{\mathcal{B}}$ , then  $A^*: \mathcal{B} \to \mathcal{B}$  (i.e.,  $A^* \in L[\mathcal{B}]$ ).

To prove (2), we have that for  $x \in \mathcal{H}_1$ ,

$$(A^*A)^*x = (\{\mathbf{J}_1^{-1}[\{[\mathbf{J}_1^{-1}A'_1\mathbf{J}_2]|_{\mathcal{B}}A\}_1]'\mathbf{J}_2\}|_{\mathcal{B}})x$$

$$= (\{\mathbf{J}_1^{-1}[\{A'_1[\mathbf{J}_2A_1\mathbf{J}_1^{-1}]|_{\mathcal{B}}\}]\mathbf{J}_2\}|_{\mathcal{B}})x$$

$$= A^*Ax.$$

It follows that the same result holds on  $\mathcal{B}$ .

The operator  $A^*A$  is selfadjoint on  $\mathcal{B}$ . By Theorem 1.1 (of Lax [L]), it is natural to expect that the same is true on  $\mathcal{H}_2$ . However, this need not be the case. To obtain a simple counterexample, recall that, in standard notation, the simplest class of bounded linear operators on  $\mathcal{B}$  is  $\mathcal{B} \otimes \mathcal{B}'$ , in the sense that:

$$\mathcal{B} \otimes \mathcal{B}' : \mathcal{B} \to \mathcal{B}$$
, by  $Au = (b \otimes l_{b'}(\cdot))u = \langle b', u \rangle b$ .

Thus, if  $l_{b'}(\cdot) \in \mathcal{B}' \setminus \mathcal{H}'_2$ , then  $J_2\{J_1^{-1}[(A_1)']J_2|_{\mathcal{B}}(u)\}$  is not in  $\mathcal{H}'_2$ , so that  $A^*A$  is not defined as an operator on all of  $\mathcal{H}_2$  and thus cannot have a bounded extension.

We can now state the correct extension of Theorem 1.1. (This result corrects an error in [GBZS].)

<span id="page-8-0"></span>**Theorem 1.4.** Let A be a bounded linear operator on  $\mathcal{B}$ . If  $\mathcal{B}' \subset \mathcal{H}_2$ , then A has a bounded extension to  $L[\mathcal{H}_2]$ , with  $||A||_{\mathcal{H}_2} \leq k ||A||_{\mathcal{B}}$  (for some positive k).

Proof. If  $T = A^*A$ , under the stated conditions, then  $\langle Tx, \mathbf{J}_2(y) \rangle = (Tx, y)_{H_2}$  is well defined for all  $x, y \in \mathcal{B}$ , and  $(Tx, y)_{H_2} = (x, Ty)_{H_2}$ . Thus, we can now apply Lax's Theorem to see that  $||T||_{\mathcal{H}_2} = ||A||_{\mathcal{H}_2}^2 \leq k^2 ||A||_{\mathcal{B}}^2$ .  $\square$ 

We can now finish our proof of Theorem 1.3.

To prove (1), let  $x \in \mathcal{B}$ , then  $(A^*Ax, x)_{\mathcal{H}_2} \geq 0$  for all  $x \in \mathcal{B}$ . Hence  $\langle A^*Ax, f_x^s \rangle \geq 0$ , so that  $A^*A$  is maximal accretive. The proof of (3), that  $I + A^*A$  is invertible, follows the same lines as in von Neumann's theorem.

**Remark 1.5.** Theorem 1.4 tells us that  $L[\mathcal{B}] \subset L[\mathcal{H}_2]$  as a continuous embedding. (In section 6, we will show that if  $\mathcal{B}$  has the approximation property, the embedding is dense.)

The algebra  $L[\mathcal{B}]$  also has a \*-operation that makes it much closer to  $L[\mathcal{H}_2]$  then expected. However, in general  $||A^*A||_{\mathcal{B}} \neq ||A||_{\mathcal{B}}^2$ . Furthermore, if  $A \neq B$ ,  $B^*$  then, unless

$$(B|_{\mathcal{H}_1})'(A|_{\mathcal{H}_1})' = (AB|_{\mathcal{H}_1})', \quad (AB)^* \neq A^*B^*.$$

Thus,  $L[\mathcal{B}]$  is a not a \*-algebra in the traditional sense.

## 2. The Hilbert Space KS<sup>2</sup>

<span id="page-9-0"></span>Theorem 1.4 is odd, given the requirement that  $\mathcal{B}' \subset \mathcal{H}_2$ . The following example shows that, if it's true at all, it does not work for one of the standard Banach-Hilbert space couples.

**Example 2.1.** Let  $\ell_1 \to \ell_2$  be the natural embedding, and let  $e_n$  be the natural unit basis. Put  $T(e_1) = e_1$  and  $T(e_n) = e_1 + e_n$  for n > 1. This operator has a natural extension to a bounded linear operator in  $\ell_1$ . Put  $x_n = n^{-1}(e_1 + \dots + e_n)$ . Then  $||x_n||_2 \to 0$ ,  $||T(x_n) - e_1||_2 \to 0$  but  $T(0) \neq e_1$ . Thus, T cannot be extended to a closed operator on  $\ell_2$ . It follows that  $\ell_2$  is not the correct Hilbert space for the extension of bounded linear operators or for the construction of adjoints for bounded linear operators on  $\ell_1$ . (Note that  $\ell_1$ ' is not contained in  $\ell_2$ .)

The purpose of this section is to construct a Hilbert space which allows us to apply Theorem 1.4 to all classical Banach spaces.

In order to construct the space of interest, first recall that Alexiewicz [AL] has shown that the class  $D(\mathbb{R})$ , of Denjoy integrable functions (restricted and wide sense), can be normed in the following manner: for  $f \in D(\mathbb{R})$ , define  $||f||_D$  by

$$(2.1)$$

It is clear that this is a norm, and it is known that  $D(\mathbb{R})$  is not complete. Replacing  $\mathbb{R}$  by  $\mathbb{R}^n$  in (2.1), for  $f \in D(\mathbb{R}^n)$ , we have:

$$(2.2) \quad \|f\|_D = \sup_{r>0} \left| \int_{\mathbf{B}_r} f(\mathbf{x}) d\mathbf{x} \right| = \sup_{r>0} \left| \int_{\mathbf{R}^n} \mathcal{E}_{\mathbf{B}_r}(\mathbf{x}) f(\mathbf{x}) d\mathbf{x} \right| < \infty,$$

where  $\mathbf{B}_r$  is any closed cube of diagonal r centered at the origin in  $\mathbb{R}^n$  with sides parallel to the coordinate axes, and  $\mathcal{E}_{\mathbf{B}_r}(\mathbf{x})$  is the indicator function of  $\mathbf{B}_r$ .

To construct the space, fix n, and let  $\mathbb{Q}^n$  be the set  $\{\mathbf{x} = (x_1, x_2 \cdots, x_n) \in \mathbb{R}^n\}$  such that  $x_i$  is rational for each i. Since this is a countable dense set in  $\mathbb{R}^n$ , we can arrange it as  $\mathbb{Q}^n = \{\mathbf{x}_1, \mathbf{x}_2, \mathbf{x}_3, \cdots\}$ . For each l and i, let  $\mathbf{B}_l(\mathbf{x}_i)$  be the closed cube centered at  $\mathbf{x}_i$ , with sides parallel to the coordinate axes and diagonal  $r_l = 2^{-l}, l \in \mathbb{N}$ . Now choose the natural order which maps  $\mathbb{N} \times \mathbb{N}$  bijectively to  $\mathbb{N}$ :

$$\{(1,1), (2,1), (1,2), (1,3), (2,2), (3,1), (3,2), (2,3), \dots\}.$$

Let  $\{\mathbf{B}_k, k \in \mathbb{N}\}$  be the resulting set of (all) closed cubes  $\{\mathbf{B}_l(\mathbf{x}_i) \mid (l,i) \in \mathbb{N} \times \mathbb{N}\}$  centered at a point in  $\mathbb{Q}^n$ , and let  $\mathcal{E}_k(\mathbf{x})$  be the indicator function of  $\mathbf{B}_k$ , so that  $\mathcal{E}_k(\mathbf{x})$  is in  $\mathbf{L}^p[\mathbb{R}^n] \cap \mathbf{L}^\infty[\mathbb{R}^n]$  for  $1 \leq p < \infty$ . Define  $F_k(\cdot)$  on  $\mathbf{L}^1[\mathbb{R}^n]$  by

(2.3) 
$$F_k(f) = \int_{\mathbb{R}^n} \mathcal{E}_k(\mathbf{x}) f(\mathbf{x}) d\mathbf{x}.$$

It is clear that  $F_k(\cdot)$  is a bounded linear functional on  $\mathbf{L}^p[\mathbb{R}^n]$  for each k,  $\|F_k\|_{\infty} \leq 1$  and, if  $F_k(f) = 0$  for all k, f = 0 so that  $\{F_k\}$  is fundamental

on  $\mathbf{L}^p[\mathbb{R}^n]$  for  $1 \leq p \leq \infty$ . Fix  $t_k > 0$  such that  $\sum_{k=1}^{\infty} t_k = 1$  and define a measure  $d\mathbf{P}(\mathbf{x}, \mathbf{y})$  on  $\mathbb{R}^n \times \mathbb{R}^n$  by:

$$d\mathbf{P}(\mathbf{x}, \mathbf{y}) = \left[\sum_{k=1}^{\infty} t_k \mathcal{E}_k(\mathbf{x}) \mathcal{E}_k(\mathbf{y})\right] d\mathbf{x} d\mathbf{y}.$$

We now define an inner product  $(\cdot)$  on  $\mathbf{L}^1[\mathbb{R}^n]$  by

(2.4) 
$$(f,g) = \int_{\mathbb{R}^n \times \mathbb{R}^n} f(\mathbf{x}) g(\mathbf{y})^* d\mathbf{P}(\mathbf{x}, \mathbf{y})$$
$$= \sum_{k=1}^{\infty} t_k \left[ \int_{\mathbb{R}^n} \mathcal{E}_k(\mathbf{x}) f(\mathbf{x}) d\mathbf{x} \right] \left[ \int_{\mathbb{R}^n} \mathcal{E}_k(\mathbf{y}) g(\mathbf{y}) d\mathbf{y} \right]^*.$$

The choice of  $t_k$  is suggested by physical analysis in another context (see Gill and Zachary [GZ]). We call the completion of  $\mathbf{L}^1[\mathbb{R}^n]$ , with the above inner product, the Kuelbs-Steadman space,  $\mathbf{KS}^2[\mathbb{R}^n]$ . Steadman [ST] constructed a version of this space by adapting an approach developed by Kuelbs [KB] for other purposes.

**Theorem 2.2.** The space  $\mathbf{KS}^2[\mathbb{R}^n]$  contains  $\mathbf{L}^p[\mathbb{R}^n]$  (for each  $p, 1 \leq p \leq \infty$ ) as continuous, compact, dense embeddings.

Proof. The proof of the first part is easy, if we notice that  $\mathbf{L}^1[\mathbb{R}^n] \cap \mathbf{L}^p[\mathbb{R}^n]$  is dense for  $1 \leq p < \infty$ . If  $f \in \mathbf{L}^\infty[\mathbb{R}^n]$ , then  $\left| \int_{B_k} f(\mathbf{x}) d\mathbf{x} \right|^2 \leqslant \|f\|_{L^\infty}^2$  for all k, so that  $\|f\|_{KS^2} \leqslant \|f\|_{L^\infty}$ . The proof of compactness follows from the fact that, if  $\{f_n\}$  is any weakly convergent sequence in  $\mathbf{L}^p[\mathbb{R}^n]$  with limit f, then  $\mathcal{E}_k(\mathbf{x}) \in \mathbf{L}^q[\mathbb{R}^n]$ ,  $1 < q \leq \infty$ , so that

$$\int_{\mathbb{R}^n} \mathcal{E}_k(\mathbf{x}) \left[ f_n(\mathbf{x}) - f(\mathbf{x}) \right] d\mathbf{x} \to 0$$

for each k. Thus, {fn} converges strongly to f in KS<sup>2</sup> [R n ]. Finally, note that dµ<sup>k</sup> = Ek(x)dx defines a measure in M[R n ], the dual space of L∞[R n and that KS<sup>2</sup> [R n ] ⊃ L 1 [R n ∗∗ = M[R n ].

The fact that L∞[R n ] ⊂ KS<sup>2</sup> [R n ], while KS<sup>2</sup> [R n ] is separable makes it clear in a very forceful manner that separability is not an inherited property.

It is of particular interest that KS<sup>2</sup> [R n ] ⊃ M[R n ], the space of bounded finitely additive set functions defined on the Borel sets B[R] n . Recall that M[R<sup>n</sup> ] contains the Dirac delta measure and the free-particle Green's function for the Feynman integral. Thus, KS<sup>2</sup> [R n ] contains the Dirac measure and the kernel for the Feynman integral as norm bounded elements (the original reason for our interest). It is clear from Theorem 1.4 that the convolution operator has a bounded extension to KS<sup>2</sup> [R n ]. This result was used in [\[GZ1\]](#page-45-1) to prove that the path integral could be rigorously constructed in exactly the manner envisioned by Feynman (see also [\[GZ2\]](#page-45-2)).

Theorem 2.3. The Hilbert space KS<sup>2</sup> [R n ] satisfies B ′ ⊂ H<sup>2</sup> for the following classical Banach spaces:

- (1) The bounded continuous functions on R n , Cb[R n ].
- (2) The bounded uniformly continuous functions on R n , UBC[R n ] .
- (3) The space L p [R n ], for 1 ≤ p ≤ ∞.

Remark 2.4. There is quite a lot of flexibility in the choice of the family of positive numbers  $\{t_k\}$ ,  $\sum_{k=1}^{\infty} t_k = 1$ . This is somewhat akin to the standard metric used for  $\mathbb{R}^{\infty}$ . Recall that, for any two points  $X, Y \in \mathbb{R}^{\infty}$ ,  $d(X,Y) = \sum_{n=1}^{\infty} \frac{1}{2^n} \frac{|X-Y|}{1+|X-Y|}$ . The family of numbers  $\{\frac{1}{2^n}\}$  can be replaced by any other sequence of positive numbers whose sum is one, without affecting the topology. We have used physical analysis to choose the family  $\{t_k\}$ , so they are interpreted as probabilities for the occurrence of a particular discrete path.

There is also some ambiguity associated with the order for  $\mathbb{Q}_n$  and the order for  $\mathbb{N} \times \mathbb{N}$ . (We have used simplicity to choose the order for  $\mathbb{N} \times \mathbb{N}$ .)

For our work, the important fact is that, for any combination of orders, the properties of  $\mathbf{KS}^2[\mathbb{R}^n]$  are invariant.

<span id="page-13-0"></span>2.1. Weak Integral. The purpose of this section is to indicate one other benefit that  $KS^2[\mathbb{R}^n]$  offers for analysis. Define the distributional (or weak) integral on  $\mathbb{R}$  by (see Talvila [TA]):

**Definition 2.5.** Let F' = DF be the weak derivative of F. We define

$$\mathcal{A}_{\mathbf{c}}(\mathbb{R}) = \{ f = DF \mid , F \in \mathcal{B}_{\mathbf{c}}(\mathbb{R}) \},$$

where

$$\mathcal{B}_{\mathbf{c}}(\mathbb{R}) = \left\{ F \in \mathbf{C}(\mathbb{R}) \,\middle| \, \lim_{x \to -\infty} F(x) = 0, \, \lim_{x \to \infty} F(x) \in \mathbb{R} \right\}.$$

If  $f \in \mathcal{A}_{\mathbf{c}}(\mathbb{R})$ , we say that  $F \in \mathcal{B}_{\mathbf{c}}(\mathbb{R})$  is the weak integral of f and write

$$F(x) = (w) \int_{-\infty}^{x} f(y) dy.$$

The following is proved in Talvila [TA].

**Theorem 2.6.** With the Alexiewicz norm  $\|\cdot\|_D$ , the space  $\mathcal{A}_{\mathbf{c}}[\mathbb{R}]$  has the following properties:

- (1) A<sub>c</sub>[ℝ] is a separable Banach space and a Banach lattice, which contains L¹[ℝ] and the Denjoy integrable functions (restricted and wide sense) as dense subsets.
- (2)  $\mathcal{A}_{\mathbf{c}}[\mathbb{R}]$  is isometrically isomorphic to  $\mathcal{B}_{\mathbf{c}}[\mathbb{R}]$ .
- (3)  $\mathcal{A}_{\mathbf{c}}[\mathbb{R}]$  is the completion of  $D(\mathbb{R})$  (space of Denjoy integrable functions).
- (4) The dual space  $\mathcal{A}_{\mathbf{c}}^*[\mathbb{R}]$  of  $\mathcal{A}_{\mathbf{c}}$ , is  $\mathcal{BV}(\mathbb{R})$  (i.e., functions of bounded variation on  $\mathbb{R}$ ).

This theorem allows us to include the restricted and wide sense Denjoy integrals in the class of distributions.

**Theorem 2.7.** The space  $\mathcal{A}_{\mathbf{c}}[\mathbb{R}]$  is a continuous dense and compact embedding in  $\mathbf{KS}^2[\mathbb{R}]$ .

*Proof.* Since  $\mathcal{E}_k(\mathbf{x}) \in \mathcal{BV}(\mathbb{R})$  for each k, compactness follows. To prove continuity, note that

$$\left| \int_{B_k} f(x) dx \right|^2 \leqslant \|f\|_D^2 \quad \forall k \Rightarrow$$

$$\|f\|_{KS}^2 = \sum_{k=1}^{\infty} t_k \left| \int_{B_k} f(x) dx \right|^2 \leqslant \sum_{k=1}^{\infty} t_k \|f\|_D^2 = \|f\|_D^2.$$

Thus,  $\mathbf{L}^1[\mathbb{R}] \subset \mathcal{A}_{\mathbf{c}}[\mathbb{R}] \subset \mathbf{KS}^2[\mathbb{R}]$  is a continuous and dense embedding.  $\square$ 

**Remark 2.8.** There is also a weak integral in  $\mathbb{R}^n$  (see [ASV] and [MO] for details). If  $f \in \mathcal{D}'(\mathbb{R}^n)$  then f is integrable if there is a function  $F \in \mathbf{C}(\mathbb{R}^n)$  such that DF = f, where  $D = \frac{\partial^n}{\partial x_1 \partial x_2 \cdots \partial x_n}$ . Thus,

$$\int_{\mathbb{R}^n} f(x)\varphi(x)dx = \int_{\mathbb{R}^n} DF(x)\varphi(x)dx = (-1)^n \int_{\mathbb{R}^n} F(x)D\varphi(x)dx,$$
for all  $\phi$  in  $\mathbf{C}_c^{\infty}(\mathbb{R}^n)$ .

In this case, we can use our generalization of the Alexiewicz norm, equation (2.2). Thus, if  $f \in D(\mathbb{R}^n)$ , then

$$||f||_D = \sup_{r>0} \left| \int_{\mathbf{B}_r} f(\mathbf{x}) d\mathbf{x} \right| = \sup_{r>0} \left| \int_{\mathbf{R}^n} \mathcal{E}_{\mathbf{B}_r}(\mathbf{x}) f(\mathbf{x}) d\mathbf{x} \right| < \infty,$$

<span id="page-15-0"></span>to construct the space  $\mathcal{A}_{\mathbf{c}}(\mathbb{R}^n)$ .

2.2. **Discussion.** Let  $\mathbf{J}_{KS}(\cdot)$  be the conjugate linear isomorphism between  $\mathbf{KS}^2[\mathbb{R}^n]$  and its dual  $\{\mathbf{KS}^2[\mathbb{R}^n]\}'$ . Since  $\mathbf{KS}^2[\mathbb{R}^n]$  contains  $\mathbf{L}^p[\mathbb{R}^n]$ ,  $1 \leq p \leq \infty$ , for each  $f(\mathbf{x}) \in \mathbf{KS}^2[\mathbb{R}^n]$ ,  $\mathbf{J}_{KS}(f_{\mathbf{x}})(\cdot) = \langle \cdot, f(\mathbf{x}) \rangle$  is a continuous linear functional on all of these spaces. However, this linear functional need not be in the dual space of any one of them. Thus, in general, we cannot automatically assume that:

$$\mathcal{B}\subset\mathcal{H}\Rightarrow\mathcal{H}'\subset\mathcal{B}'.$$

<span id="page-15-1"></span>2.3. The Corresponding  $\mathcal{H}_1$ . For completion, in this section we construct the  $\mathcal{H}_1$  version of  $\mathbf{KS}^2$ .

Recall that  $\sum_{n=1}^{\infty} \frac{1}{n^2} = \frac{\pi^2}{6}$ . Thus, setting  $\lambda_n = \frac{6}{\pi^2 n^2}$ , we see that  $\sum_{n=1}^{\infty} \lambda_n = 1$ . If we let  $\{\varphi_n\}$  be the complete orthonormal set generated

by the Hermite functions on  $\mathbb{R}^n$ , then  $\varphi_n \in \mathcal{B}$  for all the classical Banach spaces in Theorem 1.7. Thus, we can define  $\mathbf{T}_{12}$  and  $\mathbf{KS}_1^2$  by:

$$\mathbf{T}_{12}u = \sum_{n=1}^{\infty} \frac{6}{\pi^2 n^2} (u, \varphi_n)_{\mathbf{KS}^2} \varphi_n \text{ and,}$$

$$\mathbf{KS}_1^2 = \left\{ u \in B \, \middle| \, \sum_{n=1}^{\infty} \frac{\pi^2 n^2}{6} \, |(u, \varphi_n)_{\mathbf{KS}^2}|^2 < \infty \right\}, \text{ with}$$

$$(u, v)_{\mathbf{KS}_1^2} = \sum_{n=1}^{\infty} \frac{\pi^2 n^2}{6} (u, \varphi_n)_{\mathbf{KS}^2} (\varphi_n, v)_{\mathbf{KS}^2}.$$

We call  $\mathbf{KS}_1^2$  the Gross-Steadman space. Historically, Gross [GR] first proved that every real separable Banach space contains a separable Hilbert space as a dense embedding, and that this space is the support of a Gaussian measure. This was a major extension of Wiener's theory, which was based on the use of the (densely embedded Hilbert) Sobolev space  $\mathbf{H}^1[0,1] \subseteq \mathbf{C}[0,1]$  (i.e.,  $u \in \mathbf{H}^1[0,1]$  means that its first order weak derivative is in  $\mathbf{C}[0,1]$ ). Motivated by Gross' theorem, Kuelbs realized that the inclusion  $\mathbf{H}^1[0,1] \subseteq \mathbf{C}[0,1] \subset \mathbf{L}^2[0,1]$  might have an extension and prove the original version of Theorem 1.2.

<span id="page-16-0"></span>While  $\mathbf{KS}^2$  (=  $\mathcal{H}_2$ ) will be explicitly used during the remainder of the paper,  $\mathbf{KS}_1^2$  (=  $\mathcal{H}_1$ ) will be equally implicit, and  $\mathcal{B}$  shall always refer to one of the Banach spaces in Theorem 2.3.

#### 3. Closed Operators

**Definition 3.1.** A Banach space  $\mathcal{B}$  is said to be:

- (1) quasi-reflexive if dim  $\{\mathcal{B}''/\mathcal{B}\}\$  <  $\infty$ , and
- (2) nonquasi-reflexive if dim  $\{\mathcal{B}''/\mathcal{B}\} = \infty$ .

In general, it is not reasonable to expect that Theorem 1.4 will hold for all operators in C[B]. An important result by Vinokurov, Petunin and Pliczko [\[VPP\]](#page-49-0) shows that for every nonquasi-reflexive Banach space B (for example, C[0; 1] or L 1 [R n ], n ∈ N), there is a bounded linear injective operator A−<sup>1</sup> with a dense range whose inverse A is a closed densely defined linear operator which is not of the first Baire class. This means that there does not exist a sequence of bounded linear operators A<sup>n</sup> ∈ L[B] such that, for x ∈ D(A), Anx → Ax, as n → ∞.

Recall that a m-dissipative linear operator is the generator of a C0 contraction semigroup and Ran(λI − A) = B for every λ > 0 (see Pazy [\[PZ\]](#page-47-2)). Furthermore, the Yosida approximator [\[YS\]](#page-49-1), A<sup>λ</sup> = λAR(λ, A), is a bounded linear operator which converges strongly to A on D(A).

Theorem 3.2. Let A ∈ C[B], with B ′ ⊂ H2. The operator A is in the first Baire class if and only if it has an adjoint A<sup>∗</sup> .

Proof. Let H<sup>1</sup> ⊂ B ⊂ H<sup>2</sup> as in Theorem 1.2, and suppose that A has an adjoint A<sup>∗</sup> ∈ C[B]. Let T = [A∗A] 1/2 , T¯ = [AA<sup>∗</sup> ] 1/2 (the negatives of each generate C0-contraction semigroups). Since T is nonnegative, it follows that I + αT has a bounded inverse S(α) = (I + αT) −1 , for α > 0. It is also easy to see that AS(α) is bounded and, on D(A), AS(α) = S¯(α)A = (I+αT¯) <sup>−</sup>1A (see Kato [\[K\]](#page-46-0), pages 335 and 481). Using this result, we have:

$$\lim_{\alpha \to 0^+} AS(\alpha)x = \lim_{\alpha \to 0^+} \bar{S}(\alpha)Ax = Ax, \text{ for } x \in D(A).$$

It follows that A is in the first Baire class.

To prove the converse, suppose that A is in the first Baire class. Thus, there is a sequence of bounded linear operators  $\{A_n\}$  such that, for  $x \in D(A)$ ,  $A_n x \to Ax$  as  $n \to \infty$ . Since each  $A_n$  is bounded, by Theorem 1.3, each  $A_n$  has an adjoint  $A_n^*$  and both can be extended to bounded linear operators  $\bar{A}_n$ ,  $\bar{A}_n^*$  on  $\mathcal{H}_2$  (by Theorem 1.4). Furthermore, we have  $\|\bar{A}_n\|_{\mathcal{H}_2} \le k \|A_n\|_{\mathcal{B}}$  and  $\|\bar{A}_n^*\|_{\mathcal{H}_2} \le k \|A_n^*\|_{\mathcal{B}}$ . It follows that the sequence  $\{\bar{A}_n x\}$  converges for each  $x \in D(A)$ . If we define  $\bar{A}$  as the closure in  $\mathcal{H}_2$  of  $\lim_{n\to\infty} \bar{A}_n x$  for  $x \in D(A)$ , then  $\bar{A} \in \mathcal{C}[\mathcal{H}_2]$ .

Since  $\bar{A}$  is a closed densely defined linear operator, its  $\mathcal{H}_2$  adjoint,  $\bar{A}^*$  is densely defined and  $\bar{A} = \bar{A}^{**}$  (see Rudin [RU], Theorem 13.12, page 335). From this, we see that  $\bar{A}^*$  is a closed densely defined linear operator on  $\mathcal{H}_2$ . Since  $\bar{A}$  restricted to  $\mathcal{B}$  is A,  $\bar{A}^*$  restricted to  $\mathcal{B}$  defines  $A^*$ .

If  $\mathcal{B}$  is a quasi-reflexive separable Banach space, it is shown in [VPP] that every bounded linear injective operator  $A^{-1}$  with a dense range whose inverse A is a closed densely defined linear operator is of the first Baire class. Since, to our knowledge, every operator  $A \in \mathcal{C}[\mathcal{B}]$  cannot be obtained from an  $A^{-1}$  in the class of bounded linear injective operators with a dense range, it's still not known if all operators in  $\mathcal{C}[\mathcal{B}]$  are of the first Baire class (even if  $\mathcal{B}$  is reflexive). Thus, although the theorems we prove in this section hold for all operators of first Baire class, we restrict our consideration to generators of  $C_0$ -contraction semigroups.

Theorem 3.3. If A generates a C0-contraction semigroup and B ′ ⊂ H2, then:

- (1) A has a closed densely defined extension A¯ to H2, which is also the generator of a C0-contraction semigroup.
- (2) ρ(A¯) = ρ(A) and σ(A¯) = σ(A).
- (3) The adjoint of A, ¯ A¯<sup>∗</sup> , restricted to B, is the adjoint A<sup>∗</sup> of A, that is:
  - the operator A∗A > 0,
  - (A∗A) <sup>∗</sup> = A∗A and
  - I + A∗A has a bounded inverse.

## Proof. Part I

Let T(t) be the semigroup generated by A. By Theorem 1.4, as a bounded linear operator, T(t) has a bounded extension T¯(t) to H2.

We prove that T¯(t) is a C0-semigroup. (The fact that it is a contraction semigroup will follow later.) It is clear that T¯(t) has the semigroup property. To prove that it is strongly continuous, use the fact that B is dense in H<sup>2</sup> so that, for each u ∈ H2, there is a sequence {un} in B converging to u. We then have:

$$\lim_{t \to 0} \|\bar{T}(t)u - u\|_{2} \leq \lim_{t \to 0} \{\|\bar{T}(t)u - \bar{T}(t)u_{n}\|_{2} + \|\bar{T}(t)u_{n} - u_{n}\|_{2}\} + \|u_{n} - u\|_{2}$$

$$\leq k \|u - u_{n}\|_{2} + \lim_{t \to 0} \|\bar{T}(t)u_{n} - u_{n}\|_{2} + \|u_{n} - u\|_{2}$$

$$= (k+1) \|u - u_{n}\|_{2} + \lim_{t \to 0} \|T(t)u_{n} - u_{n}\|_{2} = (k+1) \|u - u_{n}\|_{2},$$

where we have used the fact that  $\bar{T}(t)u_n = T(t)u_n$  for  $u_n \in \mathcal{B}$ , and k is the constant in Theorem 1.4. It is clear that we can make the last term on the right as small as we like by choosing n large enough, so that  $\bar{T}(t)$  is a  $C_0$ -semigroup.

To prove (1), note that, if  $\bar{A}$  is the extension of A, and  $\lambda I - \bar{A}$  has an inverse, then  $\lambda I - A$  also has one, so  $\rho(\bar{A}) \subset \rho(A)$  and  $Ran(\lambda I - A)_{\mathcal{B}} \subset Ran(\lambda I - \bar{A})_{\mathcal{H}_2} \subset \overline{Ran(\lambda I - A)}_{\mathcal{H}_2}$  for any  $\lambda \in \mathbb{C}$ . For the other direction, note that, since A generates a  $C_0$ -contraction semigroup,  $\rho(A) \neq \emptyset$ . Thus, if  $\lambda \in \rho(A)$ , then  $(\lambda I - A)^{-1}$  is a continuous mapping from  $Ran(\lambda I - A)$  onto D(A) and  $Ran(\lambda I - A)$  is dense in  $\mathcal{B}$ . Let  $u \in D(\bar{A})$ , so that  $(u, \bar{A}u) \in \hat{G}(A)$ , the closure of the graph of A in  $\mathcal{H}_2$ . Thus, there exists a sequence  $\{u_n\} \subset D(A)$  such that  $\|u - u_n\|_G = \|u - u_n\|_{\mathcal{H}_2} + \|\bar{A}u - \bar{A}u_n\|_{\mathcal{H}_2} \to 0$  as  $n \to \infty$ . Since  $\bar{A}u_n = Au_n$ , it follows that  $(\lambda I - \bar{A})u = \lim_{n \to \infty} (\lambda I - A)u_n$ . However, by the boundedness of  $(\lambda I - A)^{-1}$  on  $R(\lambda I - A)$  we have that, for some  $\delta > 0$ ,

$$\left\| (\lambda I - \bar{A})u \right\|_{\mathcal{H}_2} = \lim_{n \to \infty} \left\| (\lambda I - A)u_n \right\|_{\mathcal{H}_2} \ge \lim_{n \to \infty} \delta \left\| u_n \right\|_{\mathcal{H}_2} = \delta \left\| u \right\|_{\mathcal{H}_2}.$$

It follows that  $\lambda I - \bar{A}$  has a bounded inverse and, since  $D(A) \subset D(\bar{A})$  implies that  $Ran(\lambda I - A) \subset Ran(\lambda I - \bar{A})$ , we see that  $Ran(\lambda I - \bar{A})$  is dense in  $\mathcal{H}_2$  so that  $\lambda \in \rho(\bar{A})$  and hence  $\rho(A) \subset \rho(\bar{A})$ . It follows that  $\rho(A) = \rho(\bar{A})$  and necessarily,  $\sigma(A) = \sigma(\bar{A})$ .

Since A generates a  $C_0$ -contraction semigroup, it is m-dissipative. From the Lumer-Phillips Theorem (see Pazy [PZ]), we have that  $Ran(\lambda I - A) = \mathcal{B}$  for λ > 0. It follows that A¯ is m-dissipative and Ran(λI − A¯) = H2. Thus, T¯(t) is a C0-contraction semigroup.

We now observe that the same proof applies to T¯<sup>∗</sup> (t), so that A¯<sup>∗</sup> is also the generator of a C0-contraction semigroup on H2.

Clearly A¯<sup>∗</sup> is the adjoint of A¯ so that, from von Neumann's Theorem, A¯∗A¯ has the expected properties. By a result of Kato [\[K\]](#page-46-0) (see page 276), D¯ = D(A¯∗A¯) is a core for A¯ (i.e., the set of elements {u, Au¯ } is dense in the graph, G[A¯], of A¯ for u ∈ D¯ ). From here, we see that the restriction A<sup>∗</sup> of A¯<sup>∗</sup> to B is the generator of a C0-contraction semigroup and D = D(A∗A) is a core for A. The proof of (3) for A∗A now follows.

Theorem 3.4. Let A ∈ C[B] be the generator of a C0-contraction semigroup. If B ′ ⊂ H2, then there exist a m-accretive operator R and a partial isometry W such that A = W R and D(A) = D(R).

Proof. The fact that B ′ ⊂ H<sup>2</sup> ensures that A∗A is a closed selfadjoint operator on B by Theorem 3.3. Furthermore, both A and A<sup>∗</sup> have closed densely defined extensions A¯ and A¯<sup>∗</sup> to H2. Thus, the operator Rˆ = [A¯∗A¯] 1/2 is a well-defined m-accretive selfadjoint linear operator on H2, A¯ = W¯ R¯ for some partial isometry W¯ defined on H2, and D(A¯) = D(R¯). Our proof is complete when we notice that the restriction of A¯ to B is A and R¯<sup>2</sup> restricted to B is A∗A, so that the restriction of W¯ to B is well-defined and must be a partial isometry. The equality of the domains is obvious.

# <span id="page-22-0"></span>3.1. Operators on B.

Definition 3.5. Let S be bounded, let A be closed and densely defined, and let U, V be subspaces of B:

- (1) A is said to be naturally self-adjoint if A = A<sup>∗</sup> on D(A).
- (2) A is said to be normal if AA<sup>∗</sup> = A∗A on D(A).
- (3) S is unitary if SS<sup>∗</sup> = S <sup>∗</sup>S = I.
- (4) The subspace U is ⊥ to V if for each v ∈ V and ∀u ∈ U, hv, f <sup>s</sup> u i = 0 and, for each u ∈ U and ∀v ∈ V, hu, f <sup>s</sup> v i = 0.

The last definition is transparent since, for example,

$$\langle v, f_u^s \rangle = 0 \Leftrightarrow \langle v, J_2(u) \rangle = (v, u)_2 = 0 \ \forall v \in \mathcal{V}.$$

With respect to our definition of natural selfadjointness, the following related definition is due to Palmer [\[PL\]](#page-47-3), where the operator is called symmetric. This is essentially the same as a Hermitian operator as defined by Lumer [\[LU\]](#page-47-4).

Definition 3.6. A closed densely defined linear operator A on B is called self-conjugate if both iA and −iA are dissipative.

Theorem 3.7. (Vidav-Palmer) A linear operator A, defined on B, is selfconjugate if and only if iA and −iA are generators of isometric semigroups.

Theorem 3.8. The operator A, defined on B, is self-conjugate if and only if it is naturally self-adjoint.

*Proof.* Let  $\bar{A}$  and  $\bar{A}^*$  be the closed densely defined extensions of A and  $A^*$  to  $\mathcal{H}_2$ . On  $\mathcal{H}_2$ ,  $\bar{A}$  is naturally self-adjoint if and only if  $i\bar{A}$  generates a unitary group, if and only if it is self-conjugate. Thus, both definitions coincide on  $\mathcal{H}_2$ . It follows that the restrictions coincide on  $\mathcal{B}$ .

For later reference, we note that orthogonal subspaces in  $\mathcal{H}_2$  induce orthogonal subspaces in  $\mathcal{B}$ .

**Theorem 3.9.** (Gram-Schmidt) If  $\mathcal{B}$  has a basis  $\{\varphi_i, 1 \leq i < \infty\}$  then there is an orthonormal basis  $\{\psi_i, 1 \leq i < \infty\}$  for  $\mathcal{B}$  with a corresponding set of orthonormal duality maps  $\{f_{\psi_i}^s, 1 \leq i < \infty\}$  (i.e.,  $\langle \psi_i, f_{\psi_i}^s \rangle = \delta_{ij}$ ).

Proof. Since each  $\varphi_i$  is in  $\mathcal{H}_2$ , we can construct an orthogonal set of vectors  $\{\phi_i,\ 1\leqslant i<\infty\}$  in  $\mathcal{H}_2$  by the standard Gram-Schmidt process. Set  $\psi_i=\phi_i/\|\phi_i\|_{\mathcal{B}}$  and  $\hat{f}^s_{\psi_i}=J(\psi_i)/\|\psi_i\|_{\mathcal{H}}^2$  on the subspace  $M=<\psi_i>$ . Now use the Hahn-Banach Theorem to extend  $\hat{f}^s_{\psi_i}$  to all of  $\mathcal{B}$  as in Section 1, to get  $f^s_{\psi_i}$ . From here, it is easy to check that  $\{\psi_i,\ 1\leqslant i<\infty\}$  is an orthonormal basis for  $\mathcal{B}$  with corresponding orthonormal duality maps  $\{f^s_{\psi_i},\ 1\leqslant i<\infty\}$ .  $\square$ 

We close this section with the following observation about the use of  $\mathbf{KS}^2$ . Let A be any closed densely defined positive linear operator on  $\mathcal{B}$  with a discrete positive spectrum  $\{\lambda_i\}$ . In this case, -A generates a  $C_0$ -contraction semigroup, so that it can be extended to  $\mathcal{H}_2$  with the same properties. If we compute the ratio  $\frac{\langle A\psi, f_{\psi}^s \rangle}{\langle \psi, f_{\psi}^s \rangle}$  in  $\mathcal{B}$ , it will be "close" to the

<span id="page-24-0"></span>value of  $\frac{(\bar{A}\psi,\psi)_{\mathcal{H}_2}}{(\psi,\psi)_{\mathcal{H}_2}}$  in  $\mathcal{H}_2$ . On the other hand, note that we can use the minmax theorem on  $\mathcal{H}_2$  to compute the eigenvalues and eigenfunctions of A via  $\bar{A}$  exactly on  $\mathcal{H}_2$ . Thus, in this sense, the min-max theorem holds on  $\mathcal{B}$ .

#### 4. Extension of the Poincaré inequality

<span id="page-24-1"></span>4.1. **Introduction.** There are a number of versions of the Poincaré inequality (see Evans [EV]). We consider the version that naturally appears in the theory of Markov processes. Let  $\mu$  be a Borel probability measure associated with the transition semigroup S(t) for a given Markov process with generator A. The measure  $\mu$  is called an invariant measure if:

$$\int_{\mathbb{R}^3} S(t) u(\mathbf{x}) d\mu(\mathbf{x}) = \int_{\mathbb{R}^3} u(\mathbf{x}) d\mu(\mathbf{x}), \quad t > 0,$$

for any  $u(\mathbf{x}) \in \mathbb{C}_c^{\infty}[\mathbb{R}^3]$ . If u is any function in  $L^p[\mathbb{R}^3, d\mu]$  and we set  $\bar{u} = \int_{\mathbb{R}^3} u(\mathbf{x}) d\mu(\mathbf{x})$ , it is known that for  $1 \le p < \infty$ :

$$\lim_{t \to \infty} ||S(t)u - \bar{u}||_p = 0.$$

Since the generator of S(t) is strongly elliptic, if  $u \in W^{1,p}_{\mu}[\mathbb{R}^3, d\mu]$  (the space of functions whose first order weak derivative is in  $L^p[\mathbb{R}^3, d\mu]$ ), the Poincaré inequality states that:

(4.1) 
$$\int_{\mathbb{R}^3} |u - \bar{u}|^p d\mu \le C \int_{\mathbb{R}^3} |Du(\mathbf{x})|^p d\mu(\mathbf{x}),$$

where C is a positive constant and  $\bar{u} = \int_{\mathbb{R}^3} u(\mathbf{x}) d\mu(\mathbf{x})$ .

<span id="page-25-0"></span>4.2. **Purpose.** The purpose of this section is to show that our adjoint theory allows us to extend equation (4.1) to a large class of operators, which includes all  $C_0$ -generators A = -WR,  $R = -[A^*A]^{1/2}$ , where the spectrum of R is bounded away from zero.

<span id="page-25-1"></span>In this section, we assume that  $\bar{u} = 0$ , so that  $||u||_p^p \le C ||Du||_p^p$ .

4.3. Hilbert space case. We first assume that we are working on a separable Hilbert space  $\mathcal{H}$ . In this case, for any closed densely defined linear operator A, both  $R = -[A^*A]^{1/2}$  and  $\bar{R} = -[AA^*]^{1/2}$  are generators of  $C_0$ -analytic contraction semigroups on  $\mathcal{H}$ . Furthermore, there is a unique partial isometry W such that  $A = -WR = -\bar{R}W$ , and  $A^* = -W^*\bar{R} = -RW^*$ , see Kato [K], page 334. (It should be noted that A itself is rarely a  $C_0$ -semigroup generator of any type.)

**Theorem 4.1.** Let S(t) be the analytic contraction semigroup generated by R. If, for  $u \in \mathcal{H}$ , there is a  $T \in (0, \infty)$  such that, for  $t \geq T$ ,  $||S(t)u||_{\mathcal{H}} \leq r ||u||_{\mathcal{H}}$ , with r < 1, then there exists a constant c such that, for each  $u \in D(A)$ ,  $||u||_{\mathcal{H}} \leq c ||Au||_{\mathcal{H}}$ .

*Proof.* Since S(t) is analytic, with R as its generator, we have RS(t)u = S(t)Ru for  $u \in D(A)$ . Thus,

$$|||S(t)u||_{\mathcal{H}} - ||u||_{\mathcal{H}}| \le ||S(t)u - u||_{\mathcal{H}} = \left\| \int_0^t RS(\tau)ud\tau \right\|_{\mathcal{H}}$$
$$= \left\| \int_0^t S(\tau)Rud\tau \right\|_{\mathcal{H}} \le \int_0^t ||S(\tau)Ru||_{\mathcal{H}} d\tau \le t ||Ru||_{\mathcal{H}} = t ||Au||_{\mathcal{H}}.$$

Hence, for t ≥ T,

$$||u||_{\mathcal{H}} - r ||u||_{\mathcal{H}} \le ||u||_{\mathcal{H}} - ||S(t)u||_{\mathcal{H}} \le T ||Au||_{\mathcal{H}}.$$

If we set 
$$c = \frac{T}{1-r}$$
, then  $||u||_{\mathcal{H}} \le c ||Au||_{\mathcal{H}}$ .

(Note that the proof of Theorem 4.1 does not depend on the Hilbert space structure.)

The natural question is: What are the additional conditions on R that make the above result possible? The following conditions (for separable Banach spaces) are known (see Pazy [\[PZ\]](#page-47-2)):

## Theorem 4.2. Let B be a separable Banach space. If:

(1) for some p, 1 ≤ p < ∞

$$\int_0^\infty \|S(t)u\|_{\mathcal{B}}^p dt < \infty \quad \text{for every } u \in B, \quad \text{or}$$

(2) S(t) is an analytic contraction semigroup whose generator R has a spectrum σ(R), such that

(4.2) 
$$\sigma = \sup \{ \operatorname{Re}(\lambda) : \lambda \in \sigma(R) \} < 0,$$

then there are constants M ≥ 1 and µ > 0 such that

$$||S(t)||_{\mathcal{B}} \le Me^{-\mu t}.$$

Slemrod [\[SL\]](#page-48-2), has proved a general result assuring that kS(t)k<sup>B</sup> ≤ Me−µt . The following applies to our case.

Theorem 4.3. Let S(t) be a semigroup on H. If either condition of Theorem 4.2 holds, then there exists a constant r, 0 < r < 1, such that

$$(4.3) ||S(t)||_{\mathcal{H}} \le r.$$

Proof. Under the stated conditions, kS(t)k<sup>H</sup> ≤ Me−µt. If we choose T > lnM µ and r = Me−µT , it is easy to check that inequality (4.3) is satisfied.

The above theorem applies to all closed densely defined linear operators A such that A∗A is a strictly positive operator, where R = −[A∗A] 1/2 . In this case, if we drop the analytic condition, the theorem does not hold (see Pazy [\[PZ\]](#page-47-2), example 4.2, page 117).

<span id="page-27-0"></span>4.4. Banach space case. In case we have a separable Banach space B, we assume that A is the generator of a C0-contraction semigroup and B ′ ⊂ H2.

Theorem 4.4. Let A = W R and let S(t) be the analytic contraction semigroup generated by R on B. If, for u ∈ B, there is a 0 < T < ∞ such that, for t ≥ T, kS(t)uk<sup>B</sup> ≤ r kuk<sup>B</sup> , with r < 1, then there exists a constant c such that, for each u ∈ D(A), kuk<sup>B</sup> ≤ c kAuk<sup>B</sup> .

Proof. The proof is the same as for Theorem 4.1.

Definition 4.5. Let A generate a C0-contraction semigroup and let B be a closed densely defined linear operator on B. We say that B is relatively bounded with respect to A if D(A) ⊂ D(B) and there are positive numbers a, b such that:

$$||Bu||_{\mathcal{B}} \leqslant a ||u||_{\mathcal{B}} + b ||Au||_{\mathcal{B}} \quad for \ u \in D(A).$$

The proof of the next result follows from Theorem 4.4.

<span id="page-28-0"></span>Corollary 4.6. If B is relatively bounded with respect to A(= W R) and zero is bounded away from σ(R), then there is a constant c such that

$$||Bu||_{\mathcal{B}} \leqslant c ||Au||_{\mathcal{B}} \quad for \ u \in D(A).$$

## 5. Extension Of The Spectral Theorem

<span id="page-28-1"></span>5.1. Introduction. For any selfadjoint operator in C[H], the following theorem is well-known. A proof can be found in [\[DS\]](#page-45-5), page 1192-99 (see also Reed and Simon [\[RS\]](#page-48-3) page 263).

Theorem 5.1. Let A ∈ C[H] be a selfadjoint operator, with spectrum σ(A) ⊂ R, then there exists a unique regular countably additive projectionvalued (= spectral) measure E(Ω) mapping the Borel sets, B[R], over R into H such that, for each x ∈ D(A), we have:

(1) D(A) also satisfies

$$D(A) = \left\{ x \in \mathcal{H} \mid \int_{\sigma(A)} \lambda^2 \left( \mathbf{E}(d\lambda) x, x \right)_{\mathcal{H}} < \infty \right\}$$

and

(2) 
$$Ax = \lim_{n \to \infty} \int_{-n}^{n} \lambda \mathbf{E}(d\lambda) x, \text{ for } x \in D(A).$$

(3) If g(·) is a complex-valued Borel function defined (a.e) on R, then g(A) ∈ C[H] and, for x ∈ D(g(A)) = Dg(A),

$$g(A)x = \lim_{n \to \infty} \int_{-n}^{n} g(\lambda) \mathbf{E}(d\lambda) x,$$

where

$$D_{g}(A) = \left\{ x \in \mathcal{H} \mid \int_{\sigma(A)} |g(\lambda)|^{2} \left( \mathbf{E}(d\lambda)x, x \right)_{\mathcal{H}} < \infty \right\}$$
and  $g(A^{*}) = \bar{g}(A)$ .

It is an exercise to show that E(Ω)x is of bounded variation. (For Ω = (−∞, λ], E(λ)x is called a spectral function and {E(λ)} is called a spectral family.)

Theorem 5.1 initiated the general study of operators that have a spectral representation (or functional calculus). This research has moved in many directions. The Rellich-Titchmarsh-Kato line is concerned with applications to problems in physics and applied mathematics. In this direction, one is interested in concrete detailed information about the spectrum of various specific operators subject to different constraints (see Rellich [\[RL\]](#page-48-4), Titchmarsh [\[TI\]](#page-48-5) and Kato [\[K\]](#page-46-0)). Another line of study follows more closely the approach developed by Stone and von Neumann (independently extending the bounded case by HIlbert). In this direction one seeks to extend Theorem 5.1 to a larger class of operators via operator theory and functional analysis (see Dunford and Schwartz [\[DS\]](#page-45-5) and Yosida [\[YS\]](#page-49-1)). The notes starting on <span id="page-30-0"></span>page 2089 (in [\[DS\]](#page-45-5)) are especially helpful in understanding the history (and the many other approaches).

5.2. Background. Dunford and Schwartz define a spectral operator as one that has a spectral family similar to that defined in Theorem 5.1 for selfadjoint operators. (A spectral operator is an operator with countably additive spectral measure on the Borel sets of the complex plane.) Strauss and Trunk [\[STT\]](#page-48-6) define a bounded linear operator A, on a Hilbert space H, to be spectralizable if there exists a non-constant polynomial p such that the operator p(A) is a scalar spectral operator (has a representation as in Theorem 5.1 (2)). Another interesting line of attack is represented in the book of Colojoar˘a and Foia¸s [\[CF\]](#page-44-3). where they study the class of generalized spectral operators. Here, one is not opposed to allowing the spectral resolution to exist in a generalized sense, so as to include operators with spectral singularities.

The following theorem was proven by Helffer and Sj¨ostrand [\[HSJ\]](#page-46-1) (see Proposition 7.2):

Theorem 5.2. Let g ∈ C<sup>∞</sup> 0 [R] and let gˆ ∈ C<sup>∞</sup> 0 [C] be an extension of g, with ∂gˆ ∂z<sup>ˆ</sup> = 0 on R. If A is a selfadjoint operator on H, then

$$g(A) = -\frac{1}{\pi} \iint_{\mathbf{C}} \frac{\partial \hat{g}}{\partial \bar{z}} (z - A)^{-1} dx dy.$$

This defines a functional calculus. Davies [\[DA\]](#page-44-4) showed that the above formula can be used to define a functional calculus on Banach spaces for a closed densely defined linear operator A, provided ρ(A) ∩ R = ∅. In this program the objective is to construct a functional calculus pre-supposing that the operator of concern has a reasonable resolvent.

<span id="page-31-0"></span>5.3. Problem. The basic problem that causes additional difficulty is the fact that many bounded linear operators (on H2) are of the form A = B+N, where B is normal and N is nilpotent (i.e., there is a k ∈ N, such that Nk+1 = 0, N<sup>k</sup> 6= 0). in this case, A does not have a representation with a standard spectral measure. On the other hand, T = [N∗N] 1/2 is a selfadjoint operator and, there is a unique partial isometry W such that N = W T. If E(·) is the spectral measure associated with T, then WE(Ω)x is not a spectral measure but, it is a measure of bounded variation. Thus, we just might be able to find a easier solution to the problem if we willing drop our requirement that the spectral representation be with respect to a spectral measure in the normal sense.

We begin by noting that, in either of the Strauss and Trunk [\[STT\]](#page-48-6), Helffer and Sj¨ostrand [\[HSJ\]](#page-46-1) or Davies [\[DA\]](#page-44-4) cases, the operator A is in the first Baire class. Thus, Theorem 3.2 shows that A has an adjoint and Theorem 3.4 shows that A = W R, where W is a partial isometry and R is a nonnegative selfadjoint linear operator. Before presenting our solution for the Hilbert space case, we need a few results about vector-valued functions of bounded variation.

Recall that a vector-valued function  $\mathbf{e}(\lambda)$  defined on a subset of  $\mathbb R$  to  $\mathcal H$  is of bounded variation if

$$V(\mathbf{e}, \mathbb{R}) = \sup_{P} \left\| \sum_{i=1}^{n} \left[ \mathbf{e}(b_i) - \mathbf{e}(a_i) \right] \right\|,$$

where the supremum is over all partitions P of non-overlapping intervals  $(a_i, b_i)$  in  $\mathbb{R}$  (see Hille and Phillips [HP] or Diestel and Uhl [DU]).

The next result is proved in Hille and Phillips [HP] (see page 63).

**Theorem 5.3.** Let  $\mathbf{a}(\lambda)$  be a vector-valued function from  $\mathbb{R}$  to  $\mathcal{H}$  of bounded variation. If  $h(\lambda)$  is a continuous complex-valued function on  $(a,b) \subset \mathbb{R}$ , then the following holds:

- (1) The integral  $\int_a^b h(\lambda) d\mathbf{a}(\lambda)$  exists in the  $\mathcal{H}$  norm.
- (2) If T is any operator in  $L[\mathcal{H}]$ , then  $T\mathbf{a}(\lambda)$  is of bounded variation and

$$T \int_{a}^{b} h(\lambda) d\mathbf{a}(\lambda) = \int_{a}^{b} h(\lambda) dT \mathbf{a}(\lambda).$$

<span id="page-32-0"></span>5.4. Hilbert Space case. with respect to a Hilbert space, our result shows that, in a well-defined sense, the Stone-von Neumann approach is generic.

**Theorem 5.4.** Let  $A \in \mathcal{C}[\mathcal{H}]$  be arbitrary. Then, for each  $x \in D(A)$ , there exists a vector-valued function  $\mathbf{e}_x(\lambda)$  of bounded variation such that:

(1) D(A) also satisfies

$$D(A) = \left\{ x \in \mathcal{H} \mid \int_{\sigma(A)} \lambda^2 (d\mathbf{e}_x(\lambda), x)_{\mathcal{H}} < \infty \right\}$$

and

(2)

$$Ax = \lim_{n \to \infty} \int_{-n}^{n} \lambda d\mathbf{e}_{x}(\lambda), \text{ for all } x \in D(A).$$

(3) If  $g(\cdot)$  is a complex-valued Borel function defined (a.e) on  $\mathbb{R}$ , then

$$g(A) = \lim_{n \to \infty} \int_{-n}^{n} g(\lambda) d\mathbf{e}_{x}(\lambda) \text{ for all } x \in D_{g}(A),$$

where

$$D_g(A) = \left\{ x \in \mathcal{H} \mid \int_{\sigma(A)} |g(\lambda)|^2 \left( d\mathbf{e}_x(\lambda), x \right)_{\mathcal{H}} < \infty \right\}.$$

*Proof.* To prove (1), write A = WR, where W is the unique partial isometry and  $R = [A^*A]^{1/2}$ . By Theorem 5.1, there is a spectral measure  $\mathbf{E}(\Omega)$  such that, for each  $x \in D(A) = D(R)$ :

(5.1) 
$$Rx = \lim_{n \to \infty} \int_{-n}^{n} \lambda d\mathbf{E}(d\lambda) x.$$

If we set  $\mathbf{a}_x(\lambda) = \mathbf{E}(\lambda)x$ , then  $\mathbf{a}_x(\lambda)$  is a vector-valued function of bounded variation. Furthermore, W is a partial isometry and  $W\mathbf{a}_x(\lambda)$  is of bounded variation, with  $Var(W\mathbf{a}_x, \mathbb{R}) \leq Var(\mathbf{a}_x, \mathbb{R})$ . Thus, by Theorem 5.3, for each interval (a, b),

$$W \int_{a}^{b} \lambda d\mathbf{a}_{x}(\lambda) = \int_{a}^{b} \lambda dW \mathbf{a}_{x}(\lambda).$$

Since Ax = WRx, if we set  $\mathbf{e}_x(\lambda) = W\mathbf{a}_x(\lambda)$ , we have from equation (5.1),

(5.2) 
$$Ax = \lim_{n \to \infty} \int_{-n}^{n} \lambda d\mathbf{e}_{x}(\lambda).$$

The proof of (2) and (3) are now direct adaptations of the same result in [DS].

Thus, with minor modification the the Stone-von Neumann Theorem extends to all closed densely defined linear operators on H.

Remark 5.5. Given that {E(λ)} = {E((−∞, λ])} is any spectral family, Kato [\[K\]](#page-46-0) page 358, defined |A| = R = [A∗A] 1/2 by:

$$|A| x = \lim_{n \to \infty} \int_{-n}^{n} |\lambda| d\mathbf{E}(\lambda) x$$
, for  $x \in D(A) = D(R)$ .

This allowed him to show that different spectral families lead to different selfadjoint operators. (This result was also known to Stone [\[SO\]](#page-48-7) and von Neumann [\[VN\]](#page-48-8).)

It is clear that Theorem 5.4 could have been proven after 1948 when the book (Hille version) by Hille and Phillips appeared [\[HP\]](#page-45-6). It's a fact of history that, during the same period, research on vector-valued measures and abstract integration theory was also taking flight. However, by this time, the interests of researchers in the field had shifted to the use of abstract methods in the study of operator algebras. This work led to a new version of the spectral theorem on Banach algebras via the well-known Gelfand transform (see Rudin [\[RU\]](#page-48-0), Theorem 12.22, page 306).

We should point out that the disadvantage of Theorem 5.4 is that it gives no additional information at all about the known problems of spectral theory. Thus, for the concrete problems of particular operators this is no help. However, it does tell us that, for a given A, these problems are closely related to the detailed properties of the associated partial isometry W, with A = W R. <span id="page-35-0"></span>(Here, we presume that the properties of nonnegative selfadjoint operators are well understood?)

## 5.5. Banach space case.

Theorem 5.6. If B ′ ⊂ H<sup>2</sup> and A ∈ C[B] is the generator of a C0-contraction semigroup, then there exists a unique vector-valued function ex(λ) of bounded variation such that, for each x ∈ D(A), we have:

(1) D(A) also satisfies

$$D(A) = \left\{ x \in \mathcal{B} \mid \int_{\sigma(A)} \lambda^2 \langle d\mathbf{e}_x(\lambda), f_x^s \rangle_{\mathcal{B}} < \infty \right\}$$

and

(2) Ax = limn→∞ <sup>Z</sup> <sup>n</sup> −n λdex(λ), for all x ∈ D(A).

(3) If g(·) is a complex-valued Borel function defined (a.e) on R, then g(A) ∈ C[B]. Furthermore,

$$D_g(A) = \left\{ x \in \mathcal{B} \mid \int_{\sigma(A)} |g(\lambda)|^2 \left\langle d\mathbf{e}_x(\lambda), f_x^s \right\rangle_{\mathcal{B}} < \infty \right\}$$

and

(4) 
$$g(A)x = \lim_{n \to \infty} \int_{-n}^{n} g(\lambda) d\mathbf{e}_{x}(\lambda), \text{ for all } x \in D_{g}(A).$$

Proof. By Theorem 3.4, A = W R, where W is the unique partial isometry and R = [A∗A] 1/2 . Let R¯ be the extension of R to H2. From equation (5.1), we see that there is a unique spectral measure  $\bar{\mathbf{E}}(\Omega)$  such that for each  $x \in D(\bar{R})$ :

(5.3) 
$$\bar{R}x = \lim_{n \to \infty} \int_0^n \lambda d\bar{\mathbf{E}}(d\lambda)x.$$

If we set  $\bar{\mathbf{a}}_x(\lambda) = \bar{\mathbf{E}}(\lambda)x$ , then  $\bar{\mathbf{a}}_x(\lambda)$  is a vector-valued function of bounded variation. Furthermore, if  $\bar{W}$  is the extension of W,  $\bar{W}\bar{\mathbf{a}}_x(\lambda)$  is of bounded variation, with  $Var(\bar{W}\bar{\mathbf{a}}_x,\mathbb{R}) \leq Var(\bar{\mathbf{a}}_x,\mathbb{R})$ . If we set  $\bar{\mathbf{e}}_x(\lambda) = \bar{W}\bar{\mathbf{a}}_x(\lambda)$ , by Theorem 5.3, for each interval (a,b),

$$\left\{ \bar{W} \int_{a}^{b} \lambda d\bar{\mathbf{a}}_{x}(\lambda) \right\} = \int_{a}^{b} \lambda d\bar{\mathbf{e}}_{x}(\lambda).$$

Since  $\bar{A}x = \bar{W}\bar{R}x$  and the restriction of  $\bar{A}$  to  $\mathcal{B}$  is A, we have, for all  $x \in D(A)$ ,

(5.4) 
$$Ax = \lim_{n \to \infty} \int_{-n}^{n} \lambda d\mathbf{e}_{x}(\lambda).$$

This proves (2). The proof of (1) follows from (1) in Theorem 5.1 and the definition of  $f_x^s$ . The proofs of (3) and (4) are direct adaptations of the Hilbert space case (see [RS]).

<span id="page-36-0"></span>5.6. **General Case.** In this section, we assume that, for each  $i, 1 \leq i \leq n, n \in \mathbb{N}, \mathcal{B}_i = \mathcal{B}$  is a fixed separable Banach space. We set  $\mathfrak{B} = \times_{i=1}^n B_i$ , and represent a vector  $\mathbf{x} \in \mathfrak{B}$  by  $\mathbf{x}^t = [x_1, x_2, \dots, x_n]$ . An operator  $\mathbf{A} = [A_{ij}] \in C[\mathfrak{B}]$  is defined whenever  $A_{ij} : \mathcal{B} \to \mathcal{B}$ , is in  $\mathcal{C}[\mathcal{B}]$ .

If B ′ ⊂ H<sup>2</sup> and Aij generates a C0-contraction semigroup, then by Theorem 5.3, there exists a unique vector-valued function e ij <sup>x</sup> (λ) of bounded variation such that, for each x ∈ D(Aij ), we have:

(1) D(Aij ) also satisfies

$$D(A_{ij}) = \left\{ x \in \mathcal{B} \mid \int_{\sigma(A_{ij})} \lambda^2 \left\langle d\mathbf{e}_x^{ij}(\lambda), f_x^s \right\rangle_{\mathcal{B}} < \infty \right\}$$

and

(2) <sup>A</sup>ij<sup>x</sup> = limn→∞ <sup>Z</sup> <sup>n</sup> −n λde ij x (λ), for all x ∈ D(Aij ).

(3) If g(·) is a complex-valued Borel function defined (a.e) on R then g(Aij ) ∈ C[B]. Furthermore,

$$D_g(A_{ij}) = \left\{ x \in \mathcal{B} \mid \int_{\sigma(A_{ij})} |g(\lambda)|^2 \left\langle d\mathbf{e}_x^{ij}(\lambda), f_x^s \right\rangle_{\mathcal{B}} < \infty \right\}$$

and

(4)

$$g(A_{ij})x = \lim_{n \to \infty} \int_{-n}^{n} g(\lambda) d\mathbf{e}_{x}^{ij}(\lambda), \text{ for all } x \in D_{g}(A_{ij}).$$

If we let dE(λ) = [de ij (λ)], then we can represent A and g(A) by:

$$\mathbf{A}\mathbf{x} = \lim_{n \to \infty} \int_{-n}^{n} \lambda d\boldsymbol{\mathcal{E}}(\lambda)\mathbf{x}, \text{ for all } \mathbf{x} \in D(\mathbf{A})$$

and

$$g(\mathbf{A})\mathbf{x} = \lim_{n \to \infty} \int_{-n}^{n} g(\lambda) d\mathbf{\mathcal{E}}(\lambda)\mathbf{x}$$
, for all  $\mathbf{x} \in D(\mathbf{A})$ .

## 6. Schatten Classes

<span id="page-38-0"></span>In this section, we show how our approach allows us to provide a natural definition for the Schatten class of operators on B.

Let K(B) be the class of compact operators on B and let F(B) be the set of operators of finite rank. Recall that, for separable Banach spaces, K(B) is an ideal that need not be the maximal ideal in L[B]. If M(B) is the set of weakly compact operators and N(B) is the set of operators that map weakly convergent sequences into strongly convergent sequences, it is known that both are closed two-sided ideals in the operator norm, and, in general, F(B) ⊂ K(B) ⊂ M(B) and F(B) ⊂ K(B) ⊂ N(B) (see part I of Dunford and Schwartz [\[DS\]](#page-45-5), pg. 553). For reflexive Banach spaces, K(B) = N(B) and M(B)=L[B]. For the space of continuous functions C[Ω], on a compact Hausdorff space Ω, Grothendieck [\[GO\]](#page-45-8) has shown that M(B)=N(B). On the other hand, it is shown in part I of Dunford and Schwartz [\[DS\]](#page-45-5) that, for a positive measure space, (Ω, Σ, µ), on L 1 (Ω, Σ, µ), M(B) ⊂ N(B).

We assume that B has the approximation property (i.e., every compact operator can be approximated by operators of finite rank). (Recall, that H<sup>1</sup> and H<sup>2</sup> are fixed.) Let A be a compact operator on B and let A¯ be its extension to H2. For each compact operator A¯ on H2, there exists an orthonormal set of functions {ϕ¯<sup>n</sup> |n > 1} such that

$$\bar{A} = \sum_{n=1}^{\infty} \mu_n(\bar{A}) (\cdot, \bar{\varphi}_n)_2 \bar{U} \bar{\varphi}_n,$$

where the  $\mu_n$  are the eigenvalues of  $[\bar{A}^*\bar{A}]^{1/2} = |\bar{A}|$ , counted by multiplicity and in decreasing order, and  $\bar{U}$  is the partial isometry associated with the polar decomposition of  $\bar{A} = \bar{U} |\bar{A}|$ . Without loss, we can assume that the set of functions  $\{\bar{\varphi}_n \mid n \geq 1\}$  is contained in  $\mathcal{B}$  and  $\{\varphi_n \mid n \geq 1\}$  is the normalized version in  $\mathcal{B}$ . If  $\mathbb{S}_p[\mathcal{H}_2]$  is the Schatten Class of order p in  $L[\mathcal{H}_2]$ , it is wellknown that, if  $\bar{A} \in \mathbb{S}_p[\mathcal{H}_2]$ , its norm can be represented as:

$$\|\bar{A}\|_{p}^{\mathcal{H}_{2}} = \left\{ Tr[\bar{A}^{*}\bar{A}]^{p/2} \right\}^{1/p} = \left\{ \sum_{n=1}^{\infty} \left( \bar{A}^{*}\bar{A}\bar{\varphi}_{n}, \bar{\varphi}_{n} \right)_{\mathcal{H}_{2}}^{p/2} \right\}^{1/p} = \left\{ \sum_{n=1}^{\infty} \left| \mu_{n}(\bar{A}) \right|^{p} \right\}^{1/p}.$$

**Definition 6.1.** We represent the Schatten Class of order p in  $L[\mathcal{B}]$  by:

$$\mathbb{S}_p[\mathcal{B}] = \mathbb{S}_p[\mathcal{H}_2] \cap L[\mathcal{B}] |_{\mathcal{B}}.$$

Since  $\bar{A}$  is the extension of  $A \in \mathbb{S}_p[\mathcal{B}]$ , we can define A on  $\mathcal{B}$  by

$$A = \sum_{n=1}^{\infty} \mu_n(A) \langle \cdot, f_n^s(\varphi) \rangle U\varphi_n,$$

where  $f_n^s(\varphi)$  is the Steadman duality map associated with  $\varphi_n$  and U is the restriction of  $\bar{U}$  to  $\mathcal{B}$ . The corresponding norm of A on  $\mathbb{S}_p[\mathcal{B}]$  is defined by:

$$||A||_p^{\mathcal{B}} = \left\{ \sum_{n=1}^{\infty} \langle A^* A \varphi_n, f_n^s(\varphi) \rangle^{p/2} \right\}^{1/p}.$$

**Theorem 6.2.** Let  $A \in \mathbb{S}_p[\mathcal{B}]$ , then  $||A||_p^{\mathcal{B}} = ||\bar{A}||_p^{\mathcal{H}_2}$ .

*Proof.* It is clear that  $\{\varphi_n \mid n \geq 1\}$  is a set of eigenfunctions for  $A^*A$  on  $\mathcal{B}$ . Furthermore, by our extension of Lax's Theorem,  $A^*A$  is selfadjoint and the point spectrum of  $A^*A$  is unchanged by its extension to  $\mathcal{H}_2$ . It follows that  $A^*A\varphi_n = |\mu_n|^2 \varphi_n$ , so that

$$\langle A^* A \varphi_n, f_n^s(\varphi) \rangle = \frac{|\mu_n|^2}{\|\varphi_n\|_2^2} (\varphi_n, \varphi_n)_2 = |\mu_n|^2,$$

and

$$||A||_{p}^{\mathcal{B}} = \left\{ \sum_{n=1}^{\infty} \langle A^* A \varphi_n, f_n^s(\varphi) \rangle^{p/2} \right\}^{1/p} = \left\{ \sum_{n=1}^{\infty} |\mu_n|^p \right\}^{1/p} = ||\bar{A}||_{p}^{\mathcal{H}_2}.$$

**Lemma 6.3.** If  $\mathcal{B}$  has the approximation property, the embedding of  $L[\mathcal{B}]$  in  $L[\mathcal{H}_2]$  is both continuous and dense.

*Proof.* Recall that the embedding is continuous by Theorem 1.4. Since  $\mathcal{B}$  has the approximation property, the finite rank operators  $\mathbb{F}(\mathcal{B})$  on  $\mathcal{B}$  are dense in the finite rank operators  $\mathbb{F}(\mathcal{H}_2)$  on  $\mathcal{H}_2$ . It follows that  $\mathbb{S}_p[\mathcal{B}]$  is dense in  $\mathbb{S}_p[\mathcal{H}_2]$ . In particular,  $\mathbb{S}_1[\mathcal{B}]$  is dense in  $\mathbb{S}_1[\mathcal{H}_2]$  and, since  $\mathbb{S}_1[\mathcal{H}_2]^* = L[\mathcal{H}_2]$ , we see that  $\mathbb{S}_1[\mathcal{B}]^* = L[\mathcal{B}]$  must be dense in  $L[\mathcal{H}_2]$ .

It is clear that much of the theory of operator ideals on Hilbert spaces extend to separable Banach spaces in a straightforward way. We state a few of the more important results to give a sense of the power provided by the existence of adjoints. The first result extends theorems due to Weyl [WY], Horn [HO], Lalesco [LE] and Lidskii [LI]. (The methods of proof for Hilbert spaces carry over without much difficulty.)

**Theorem 6.4.** Let  $A \in \mathbb{K}(\mathcal{B})$ , the set of compact operators on  $\mathcal{B}$ , and let  $\{\lambda_n\}$  be the eigenvalues of A counted up to algebraic multiplicity. If  $\Phi$  is a

mapping on [0,∞] which is nonnegative and monotone increasing, then we have:

(1) (Weyl)

$$\sum\nolimits_{n = 1}^{\mathbf{N}} \Phi\left(\left|\lambda_n(A)\right|\right) \leqslant \sum\nolimits_{n = 1}^{\mathbf{N}} \Phi\left(\mu_n(A)\right)$$

and

(2) (Horn)

$$\sum\nolimits_{n = 1}^{\mathbf{N}} \Phi\left(\left|\lambda_n(A_1 A_2)\right|\right) \leqslant \sum\nolimits_{n = 1}^{\mathbf{N}} \Phi\left(\mu_n(A_1)\mu_n(A_2)\right).$$

In case A ∈ S1(B), we have:

(3) (Lalesco)

$$\sum\nolimits_{n=1}^{\mathbf{N}} |\lambda_n(A)| \leqslant \sum\nolimits_{n=1}^{\mathbf{N}} \mu_n(A)$$

and

(4) (Lidskii)

$$\sum\nolimits_{n=1}^{\mathbf{N}} \lambda_n(A) = Tr(A).$$

<span id="page-41-0"></span>6.1. Discussion. In a Hilbert space H, the Schatten classes Sp(H) are the only ideals in K(H), and S1(H) is minimal. In a Banach space, this is far from true. A complete history of the subject can be found in the recent book by Pietsch [\[PI1\]](#page-47-5) (see also Retherford [\[RE\]](#page-47-6), for a nice review). We limit this discussion to a few major topics in the subject. First, Grothendieck [\[GO\]](#page-45-8) defined an important class of nuclear operators as follows:

**Definition 6.5.** If  $A \in \mathbb{F}(\mathcal{B})$  (the operators of finite rank), define the ideal  $\mathbf{N}_1(\mathcal{B})$  by:

$$\mathbf{N}_1(\mathcal{B}) = \{ A \in \mathbb{F}(\mathcal{B}) \mid \mathbf{N}_1(A) < \infty \},$$

where

$$\mathbf{N}_{1}(A) = \operatorname{glb}\left\{\sum_{n=1}^{m} \|f_{n}\| \|\phi_{n}\| \middle| f_{n} \in \mathcal{B}', \ \phi_{n} \in \mathcal{B}, \ A = \sum_{n=1}^{m} \phi_{n} \left\langle \cdot , f_{n} \right\rangle \right\}$$

and the greatest lower bound is over all possible representations for A.

Grothendieck has shown that  $\mathbf{N}_1(\mathcal{B})$  is the completion of the finite rank operators.  $\mathbf{N}_1(\mathcal{B})$  is a Banach space with norm  $\mathbf{N}_1(\cdot)$ , and is a two-sided ideal in  $\mathbb{K}(\mathcal{B})$ . It is easy to show that:

Corollary 6.6.  $\mathbb{M}(\mathcal{B})$ ,  $\mathbb{N}(\mathcal{B})$  and  $\mathbf{N}_1(\mathcal{B})$  are two-sided \*ideals.

In order to compensate for the (apparent) lack of an adjoint for Banach spaces, Pietsch [PI2], [PI3] defined a number of classes of operator ideals for a given  $\mathcal{B}$ . Of particular importance for our discussion is the class  $\mathbb{C}_p(\mathcal{B})$ , defined by

$$\mathbb{C}_p(\mathcal{B}) = \left\{ A \in \mathbb{K}(\mathcal{B}) \mid \mathbb{C}_p(A) = \sum_{i=1}^{\infty} [s_i(A)]^p < \infty \right\},\,$$

where the singular numbers  $s_n(A)$  are defined by:

$$s_n(A) = \inf \{ ||A - K||_{\mathcal{B}} \mid \text{rank of } K \leq n \}.$$

Pietsch has shown that  $\mathbb{C}_1(\mathcal{B}) \subset \mathbf{N}_1(\mathcal{B})$ , while Johnson et al [JKMR] have shown that for each  $A \in \mathbb{C}_1(\mathcal{B})$ ,  $\sum_{n=1}^{\infty} |\lambda_n(A)| < \infty$ . On the other

hand, Grothendieck [GO] has provided an example of an operator A in N1(L∞[0, 1]) with P<sup>∞</sup> <sup>n</sup>=1 |λn(A)| = ∞ (see Simon [SI], pg. 118). Thus, it follows that, in general, the containment is strict. It is known that, if C1(B) = N1(B), then B is isomorphic to a Hilbert space (see Johnson et al). It is clear from the above discussion, that:

Corollary 6.7. Cp(B) is a two-sided \*ideal in K(B), and S1(B) ⊂ N1(B).

For a given separable Banach space, it is not clear how the spaces Cp(B) of Pietsch relate to our Schatten Classes Sp(B) (clearly Sp(B) ⊆ Cp(B)). Thus, one question is that of the equality of Sp(B) and Cp(B). (We suspect that S1(B) = C1(B).)

# 7. Conclusion

<span id="page-43-0"></span>In this paper, we have refined and extended the work in [\[GBZS\]](#page-45-0) to develop a complete theory of adjoints for bounded linear operators on separable Banach spaces. We have further identified the obstacles to a similar program for closed densely defined linear operators. A major result in this case is that all operators of Baire class one have an adjoint. For applications, we restricted our consideration to generators of C0-contraction semigroups. We first used the polar decomposition property to extend the Poincar´e inequality. Then, the polar decomposition property, along with a few results for vector measures and vector-valued functions allowed us to extend the spectral theorem to all closed densely defined linear operators on separable Hilbert spaces. Using our adjoint theory, we were able to extend the spectral theorem to all bounded linear operators and all generators of C0-contraction semigroups on separable Banach spaces. As a final application, we introduced a new class of <sup>∗</sup>operator ideals on Banach spaces that parallel the Schatten class for Hilbert spaces.

Acknowledgements. During the course of the development of this work, we have benefited from important critical remarks from Professor Ioan I. Vrabie.

We would like to sincerely thank Professors Jerome Goldstein and Anatolij Pliczko for important correspondence on spectral operators and closed operators of Baire class on Banach spaces. They also identified a few errors in an earlier draft, which led to an improvement in the paper.

# <span id="page-44-0"></span>References

- <span id="page-44-1"></span>[AL] A. Alexiewicz, Linear functionals on Denjoy-integrable functions, Colloq. Math. 1 (1948), 289-293.
- <span id="page-44-2"></span>[ASV] D. D. Ang, K. Schmitt and L. K. Vy, A multidimensional analogue of the Denjoy-Perron-Henstock-Kurzweil integral, Bull. Belg. Math. Soc. Simon Stevin 4 (1997), 355371.
- <span id="page-44-3"></span>[CF] I. Colojoar˘a and C. Foia¸s, Theory of generalized spectral operators, Gordon Breach, (1968).
- <span id="page-44-4"></span>[DA] E. B. Davies, The Functional Calculus, J. London Mat. Soc. Vol. 52 (1995) 166-176.

- <span id="page-45-5"></span>[DS] N. Dunford and J. T. Schwartz, Linear Operators Part II: Spectral Theory, Wiley Classics edition, Wiley Interscience (1988).
- <span id="page-45-7"></span>[DU] J. Diestel and J. J. Uhl, Jr, Vector Measures , Math. Surveys 15, Amer. Math. Soc. Providence, RI, (1977).
- <span id="page-45-4"></span>[EV] L. C. Evans, Partial Differential Equations, AMS Graduate Studies in Math. 18, Providence, R.I, 1998.
- <span id="page-45-0"></span>[GBZS] T. Gill, S. Basu, W. W. Zachary and V. Steadman, Adjoint for operators in Banach spaces, Proceedings of the American Mathematical Society, 132 (2004), 1429-1434.
- <span id="page-45-8"></span>[GO] A. Grothendieck, Products tensoriels topologiques et espaces nucleaires, Memoirs of the American Mathematical Society, 16 (1955).
- <span id="page-45-3"></span>[GR] L. Gross, Abstract Wiener spaces, Proc. Fifth Berkeley Symposium on Mathematics Statistics and Probability, (1965), 31-42.
- [GZ] T. L. Gill and W. W. Zachary, Foundations for relativistic quantum theory I: Feynman's operator calculus and the Dyson conjectures, Journal of Mathematical Physics 43 (2002), 69-93.
- <span id="page-45-1"></span>[GZ1] T. L. Gill and W. W. Zachary, Banach Spaces for the Feynman integral, Real Analysis Exchange 34(2) (2008)/(2009), 267-310.
- <span id="page-45-2"></span>[GZ2] T. L. Gill and W. W. Zachary, A New Class of Banach Spaces, Journal of Physics A: Math. and Gen. 41 (2008), 495206.
- [HO] A. Horn, On the singular values of a product of completely continuous operators, Proc. Nat. Acad. Sci. 36 (1950), 374–375.
- <span id="page-45-6"></span>[HP] E. Hille and R. S. Phillips, Functional Analysis and Semigroups, Amer. Math. Soc. Colloq. Pub. 31, Amer. Math. Soc. Providence, RI, (1957).

- [HS] R. Henstock, The General Theory of Integration, Clarendon Press, Oxford, (1991).
- <span id="page-46-1"></span>[HSJ] B. Helffer, and J. Sj¨ostrand, quation de Schrdinger avec champ magnetique et quation de Harper, Schrdinger Operators, (Snderborg, 2988) eds. H. Holden and A Jensen, Lecture Notes in Phys., vol. 345, Springer-Verlag, Berlin, (1989), 118–197.
- [JKMR] W. B. Johnson, H. Konig, B. Maurey and J. R. Retherford, Eigenvalues of p-summing and l<sup>p</sup> type operators in Banach space, J. Funct. Anal. 32 (1978), 353–380.
- <span id="page-46-0"></span>[K] T. Kato, Perturbation Theory for Linear Operators, second ed. Springer-Verlag, New York, (1976).
- [KB] J. Kuelbs, Gaussian measures on a Banach space, Journal of Functional Analysis 5 (1970), 354–367.
- [KF] W. E. Kaufman, A stronger metric for closed operators in Hilbert spaces, Proc. Amer. Math. Soc. 90 (1984), 83–87.
- [KPS] S.G. Krein, Ju.I. Petunin and E.M. Semenov, Interpolation of Linear Operators, Nauka, Moscow, 1978; (English transl.), Translations Monographs 54, Amer. Math. Soc. Providence, R.I. (1982).
- [KW] J. Kurzweil, Nichtabsolut konvergente Integrale, Teubner-Texte z¨ur Mathematik, Band 26, Teubner Verlagsgesellschaft, Leipzig, (1980).
- [LE] T. Lalesco, Une theoreme sur les noyaux composes, Bull. Acad. Sci. 3 (1914/15), 271–272.
- [LI] V. B. Lidskii, Non-self adjoint operators with a trace, Dokl. Akad. Nauk. SSSR 125 (1959), 485-487.

- [LP] G. Lumer and R. S. Phillips, Dissipative operators in a Banach space, Pacific J. Math. 11 (1961), 679-698.
- <span id="page-47-4"></span>[LU] G. Lumer, Spectral operators, Hermitian operators and bounded groups, Acta. Sci. Math. (Szeged) 25 (1964), 75-85.
- <span id="page-47-1"></span>[MO] P. Mikusi´nksi and K. Ostaszewski, Embedding Henstock integrable functions into the space of Schwartz distributions, Real Anal. Exchange 14(1988-89), 24-29.
- <span id="page-47-0"></span>[L] P. D. Lax, Symmetrizable linear tranformations. Comm. Pure Appl. Math. 7 (1954), 633–647.
- <span id="page-47-3"></span>[PL] T. W. Palmer, Unbounded normal operators on Banach spaces, Trans. Amer. Math. Sci. 133 (1968), 385-414.
- [PF] W. F. Pfeffer, The Riemann Approach to Integration: Local Geometric Theory, Cambridge Tracts in Mathematics 109, Cambridge University Press, (1993).
- <span id="page-47-5"></span>[PI1] A. Pietsch, History of Banach Spaces and Operator Theory, Birkh¨auser, Boston, (2007).
- <span id="page-47-7"></span>[PI2] A. Pietsch, Einige neue Klassen von kompacter linear Abbildungen, Revue der Math. Pures et Appl. (Bucharest), 8 (1963), 423–447.
- <span id="page-47-8"></span><span id="page-47-2"></span>[PI3] A. Pietsch, Eigenvalues and s-Numbers , Cambridge University Press, (1987).
- [PZ] A. Pazy, Semigroups of Linear Operators and Applications to Partial Differential Equations Applied Mathematical Sciences, 44, Springer New York, (1983).
- <span id="page-47-6"></span>[RE] J. R. Retherford, Applications of Banach ideals of operators, Bull. Amer. Math. Soc. 81 (1975), 978-1012.

- <span id="page-48-4"></span>[RL] F. Rellich, St¨orungsterie der Spektralzerlegung V., Math. Ann. 118 (1940), 462-484.
- <span id="page-48-3"></span>[RS] M. Reed and B. Simon, Methods of Modern Mathematical Physics I: Functional Analysis, Academic Press, New York, (1972).
- <span id="page-48-0"></span>[RU] W. Rudin, Functional Analysis, McGraw-Hill Press, New York, (1973).
- [SI] B. Simon, Trace Ideals and their Applications, London Mathematical Society Lecture Notes Series 35, Cambridge University Press, New York, (1979).
- <span id="page-48-2"></span>[SL] M. Slemrod, Asymtotic behavior of C0-semigroups as determined by the spectrum of the generator, Indiana. Univ. Math. J. 25 (1976), 783-792.
- <span id="page-48-7"></span>[SO] M. H. Stone, Linear Transformations in Hilbert Space , Math. Surveys 15, Amer. Math. Soc. Colloq. Publ. 15, Providence, RI, (1932).
- [ST] V. Steadman, Theory of operators on Banach spaces, Ph.D thesis, Howard University, 1988.
- <span id="page-48-6"></span>[STT] V. A. Strauss and C. Trunk, Spectralizable Operators, Integr. Equ. Oper. Theory 61 (2008), 413-422.
- <span id="page-48-1"></span>[TA] E. Talvila, The distributional Denjoy integral, Real Analysis Exchange 33 (2008), 51-82.
- <span id="page-48-5"></span>[TI] E. C. Titchmarsh, Some theorems on perturbation theory V., J. Analyse. Math. Soc. 4 (1954/56), 187-208.
- <span id="page-48-8"></span>[VN] J. von Neumann, Uber adjungierte Funktionaloperatoren, ¨ Annals of Mathematics 33 (1932), 294-310.
- [WY] H. Weyl, Inequalities between the two kinds of eigenvalues of a linear transformation, Proc. Nat. Acad. Sci. 35, (1949), 408-11.

- <span id="page-49-0"></span>[VPP] V. A. Vinokurov, Yu. Petunin and A. N. Pliczko, Measurability and Regularizability mappings inverse to continuous linear operators (in Russian), Mat. Zametki. 26 (1979), no. 4, 583-591. English translation: Math. Notes 26 (1980), 781-785.
- <span id="page-49-1"></span>[YS] K. Yosida, Functional Analysis, second ed. Springer-Verlag, New York, (1968).

(Tepper L. Gill) Departments of Mathematics, Physics, and Electrical & Computer Engineering, Howard University, Washington DC 20059, USA, Email : tgill@howard.edu

(Francis Mensah) Department of Mathematics, Howard University, Washington DC 20059, USA, E-mail : mensah3@yahoo.com

(Woodford W. Zachary) Departments of Mathematics, and Electrical & Computer Engineering,, Howard University, Washington DC 20059, USA, E-mail : wwzachary@earthlink.net