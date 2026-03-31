## ADJOINT FOR OPERATORS IN BANACH SPACES

T. L. GILL, S. BASU, W. W. ZACHARY, AND V. STEADMAN

Abstract. In this paper we show that a result of Gross and Kuelbs, used to study Gaussian measures on Banach spaces, makes it possible to construct an adjoint for operators on separable Banach spaces. This result is used to extend well-known theorems of von Neumann and Lax. We also partially solve an open problem on the existence of a Markushevich basis with unit norm and prove that all closed densely defined linear operators on a separable Banach space can be approximated by bounded operators. This last result extends a theorem of Kaufman for Hilbert spaces and allows us to define a new metric for closed densely defined linear operators on Banach spaces. As an application, we obtain a generalization of the Yosida approximator for semigroups of operators.

# Introduction

One of the greatest impediments to the development of a theory of operators on Banach spaces that parallels the corresponding theory on Hilbert spaces is the lack of a suitable notion of an adjoint operator for these spaces. It is an interesting fact of history that the tools needed were being developed in probability theory during the time of greatest need. It was in 1965, when Gross [\[G\]](#page-9-0) first proved that every real separable Banach space contains a separable Hilbert space as a dense embedding, and this (Banach) space is the support of a Gaussian measure. Gross' theorem was a far reaching generalization of Wiener's theory, which was based on the use of

<sup>1991</sup> Mathematics Subject Classification. Primary (45) Secondary(46) .

Key words and phrases. Adjoints, Banach space embeddings, Hilbert spaces.

the (densely embedded Hilbert) Sobolev space H<sup>1</sup> [0, 1] ⊆ C[0, 1]. Later, Kuelbs [\[K\]](#page-9-1) generalized Gross' theorem to include the fact that H<sup>1</sup> [0, 1] ⊆ C[0, 1] ⊆ L 2 [0, 1]. This Gross-Kuelbs theorem can be stated for our purposes as:

Theorem 1. (Gross-Kuelbs) *Suppose* B *is a separable Banach space. Then there exist separable Hilbert spaces* H1, H<sup>2</sup> *and a positive trace class operator* T<sup>12</sup> *defined on* H<sup>2</sup> *such that* H<sup>1</sup> ⊆ B ⊆ H<sup>2</sup> *(all as continuous dense embeddings), and* T<sup>12</sup> *determines* H<sup>1</sup> *when* B *and* H<sup>2</sup> *are given.*

# Purpose

The purpose of this paper is to show that the Gross-Kuelbs theorem makes it possible to give an essentially unique definition of the adjoint for operators on separable Banach spaces. This definition has all the expected properties. In particular, we show that, for each bounded linear operator A, there exists A<sup>∗</sup> , with A∗A maximal accretive, self adjoint (A∗A) <sup>∗</sup> = A∗A, and I + A∗A is invertible. Although our main interest is in the construction of a generalized Yosida approximator for semigroups of operators that will be used elswhere, this adjoint has a number of important implications for other aspects of operator theory. As a sampling, we provide generalizations of theorems due to von Neumann [\[VN\]](#page-9-2), Lax [\[L\]](#page-9-3), and Kaufman [\[Ka\]](#page-9-4) to Banach spaces. We also partially solve an open problem on the existence of a Markushevich basis with unit norm.

## Background

In what follows, we let L[B], L[H] denote the bounded linear operators on B, H respectively. By a duality map, φx, defined on B, we mean any linear functional φ<sup>x</sup> ∈ {f ∈ B ′ | < x, f >= kxk 2 <sup>B</sup>, x ∈ B}, where < . > is the natural pairing between a Banach space and its dual. Let  $\mathbf{J} : \mathbf{H} \longrightarrow \mathbf{H}'$  be the standard conjugate isomorphism between a Hilbert space and its dual, so that  $\langle x, \mathbf{J}(x) \rangle = (x, x)_{\mathbf{H}} = ||x||_{\mathbf{H}}^2$ . We define the *special duality map* of  $\mathbf{B}$  associated with  $\mathbf{H}$  by:

$$\phi_x^s = \frac{\|x\|_B^2}{\|x\|_H^2} \mathbf{J}(x).$$

It is easy to check that  $\phi_x^s$  is a duality map for **B**. A closed densely defined operator **A** is called maximal accretive if  $\langle \mathbf{A}x, \phi_x \rangle \geq 0$  for all  $x \in D(\mathbf{A})$  and **A** has no proper extension. The following results due to von Neumann [VN] and Lax [L] are listed for reference.

**Theorem 2.** (von Neumann) For any closed densely defined linear operator  $\mathbf{A}$  on a Hilbert space  $\mathbf{H}$ , the operators  $\mathbf{A}^*\mathbf{A}$  and  $\mathbf{I} + \mathbf{A}^*\mathbf{A}$  are selfadjoint, and  $\mathbf{I} + \mathbf{A}^*\mathbf{A}$  has a bounded inverse.

**Theorem 3.** (Lax) Let  $\mathbf{H_2}$  be given so that  $\mathbf{B} \subseteq \mathbf{H_2}$  densely. If  $\mathbf{A}$  is a bounded linear operator on  $\mathbf{B}$  such that  $\mathbf{A}$  is selfadjoint (i.e.,  $(\mathbf{A}x,y)_{\mathbf{H_2}} = (x,\mathbf{A}y)_{\mathbf{H_2}} \quad \forall x,y,\in \mathbf{B}$ ), then  $\mathbf{A}$  is bounded on  $\mathbf{H_2}$  and  $\|\mathbf{A}\|_{\mathbf{H_2}} \leq \|\mathbf{A}\|_{\mathbf{B}}$ .

#### Main Results

Let us fix  $\mathbf{H_1}, \mathbf{H_2}$  such that  $\mathbf{H_1} \subseteq \mathbf{B} \subseteq \mathbf{H_2}$  as continuous dense embeddings, and, without loss of generality, assume that for  $x \in \mathbf{H_1}, \|x\|_2 \le \|x\|_{\mathbf{B}} \le \|x\|_1$ .

The first result is not new and is, in fact, well known. We present it because the proof is new and uses specific information about the relationship between  ${\bf B}$  and  ${\bf H_2}$ .

**Theorem 4.** Every closed densely defined linear operator on **B** extends to a closed densely defined linear operator on **H**<sub>2</sub>.

Proof. Let  $J_2: H_2 \longrightarrow H'_2$  denote the standard conjuate isomorphism. Then, as **B** is strongly dense in  $H_2$ ,  $J_2[B] \subset H'_2 \subset B'$  is (strongly) dense in  $H'_2$ . If **A** is any closed densely defined linear operator on **B** with domain D(A), then A' (the **B** adjoint of **A**) is closed on **B'**. In addition,  $A'|_{H'_2}$  is closed and, for each  $x \in D(A)$ ,  $J_2(x) \in H'_2$  and  $Ay, J_2(x) > b$  is well defined  $b \in D(A)$ . Hence  $b \in D(A')$  for all  $b \in D(A)$ . Since  $b \in D(A')$  is strongly dense in  $b \in D(A')$  is strongly dense in  $b \in D(A')$  is strongly dense in  $b \in D(A')$  is strongly dense in  $b \in D(A')$  is a closed densely defined operator on  $b \in D(A')$ . Thus, as  $b \in D(A')$  is reflexive,  $b \in D(A')$  is a closed densely defined operator on  $b \in D(A')$ .

In the next theorem, we prove that every bounded linear operator  $\mathbf{A}$  on  $\mathbf{B}$  has a well defined adjoint. The result is actually true for any closed densely defined linear operator on  $\mathbf{B}$  but, in this case, for each  $\mathbf{A}$  we must have  $\mathbf{H}_1 \subseteq D(\mathbf{A})$  so, in general, a different  $\mathbf{H}_1$  is required for each operator. It should also be noted that, although  $\mathbf{H}_1$  and  $\mathbf{H}_2$  are required to obtain our adjoint, it is not hard to show that any two adjoint operators for  $\mathbf{A}$  will differ by a similarity transformation of unitary operators (see Theorem 11).

**Theorem 5.** Let **B** be a separable Banach space with  $\mathbf{A} \in L[\mathbf{B}]$ . Then there exists  $\mathbf{A}^* \in L[\mathbf{B}]$  such that:

- (1)  $\mathbf{A}^*\mathbf{A}$  is maximal accretive.
- (2)  $(\mathbf{A}^*\mathbf{A})^* = \mathbf{A}^*\mathbf{A}$ , and
- (3)  $I + A^*A$  has a bounded inverse.

*Proof.* If we let  $\mathbf{J}_i: \mathbf{H}_i \to \mathbf{H'}_i$ , (i=1,2), then  $\mathbf{A}_1 = \mathbf{A}_{|\mathbf{H}_1}: \mathbf{H}_1 \longrightarrow \mathbf{H}_2$ , and  $\mathbf{A'}_1: \mathbf{H'}_2 \longrightarrow \mathbf{H'}_1$ .

It follows that  $\mathbf{A}'_1\mathbf{J}_2: \mathbf{H_2} \longrightarrow \mathbf{H}'_1$  and  $\mathbf{J}_1^{-1}\mathbf{A}'_1\mathbf{J}_2: \mathbf{H}_2 \to \mathbf{H}_1 \subset \mathbf{B}$  so that, if we define  $\mathbf{A}^* = [\mathbf{J}_1^{-1}\mathbf{A}'_1\mathbf{J}_2]_{\mathbf{B}}$ , then  $\mathbf{A}^*: \mathbf{B} \to \mathbf{B}$  (i.e.,  $\mathbf{A}^* \in L[\mathbf{B}]$ ). To prove  $1, \mathbf{J}'_i = \mathbf{J}_i$  and, if  $x \in \mathbf{B}$ , then  $\langle \mathbf{A}^*\mathbf{A}x, \mathbf{J}_2(x) \rangle = \langle \mathbf{A}x, (\mathbf{A}^*)'\mathbf{J}_2(x) \rangle$ . Using the above definition of  $\mathbf{A}^*$ , we get that  $(\mathbf{A}^*)'\mathbf{J}_2(x) = \{[\mathbf{J}_1^{-1}\mathbf{A}'_1\mathbf{J}_2]_{\mathbf{B}}\}'\mathbf{J}_2(x) = [\mathbf{J}_2\mathbf{A}_1\mathbf{J}_1^{-1}]\mathbf{J}_2(x) = \mathbf{J}_2(\mathbf{A}_1x)$ . Since, for  $x \in \mathbf{H}_1$ ,  $\mathbf{A}_1x = \mathbf{A}x$  and

$$\langle \mathbf{A}^* \mathbf{A} x, \phi_x^s \rangle = \frac{\|x\|_{\mathbf{B}}^2}{\|x\|_2^2} \langle \mathbf{A} x, \mathbf{J}_2(\mathbf{A}_1 x) \rangle = \frac{\|x\|_{\mathbf{B}}^2}{\|x\|_2^2} \|\mathbf{A} x\|_2^2 \ge 0,$$

we have that  $\mathbf{A}^*\mathbf{A}$  is accretive on a dense set. Thus,  $\mathbf{A}^*\mathbf{A}$  is accretive on  $\mathbf{B}$ . It is maximal accretive because it has no proper extension. To prove 2, we have that for  $x \in \mathbf{H}_1$ ,

$$(\mathbf{A}^* \mathbf{A})^* x = (\{\mathbf{J}_1^{-1} [\{[\mathbf{J}_1^{-1} \mathbf{A}'_1 \mathbf{J}_2]|_{\mathbf{B}} \mathbf{A}\}_1]' \mathbf{J}_2\}|_{\mathbf{B}}) x$$

$$= (\{\mathbf{J}_1^{-1} [\{\mathbf{A}'_1 [\mathbf{J}_2 \mathbf{A}_1 \mathbf{J}_1^{-1}]|_{\mathbf{B}}\}] \mathbf{J}_2\}|_{\mathbf{B}}) x$$

$$= \mathbf{A}^* \mathbf{A} x.$$

It follows that the same result holds on **B**. Finally, the proof that  $\mathbf{I} + \mathbf{A}^* \mathbf{A}$  is invertible follows the same lines as in von Neumann's theorem.

**Theorem 6.** Every bounded linear operator on  $\mathbf{B}$  extends to a bounded linear operator on  $\mathbf{H}_2$  and  $\|\mathbf{A}\|_{\mathbf{H}_2}^2 \leq C\|\mathbf{A}\|_{\mathbf{B}}^2$  for some constant C.

Proof.: For any bounded linear operator  $\mathbf{A}$  defined on  $\mathbf{B}$ , let  $\mathbf{T} = \mathbf{A}^*\mathbf{A}$ . By Theorem 1,  $\mathbf{T}$  extends to a closed linear operator  $\mathbf{T}$  on  $\mathbf{H}_2$ . As  $\mathbf{T}$  is selfadjoint on  $\mathbf{B}$ , by Lax's theorem,  $\mathbf{T}$  is bounded on  $\mathbf{H}_2$  and  $\|\mathbf{A}^*\mathbf{A}\|_{\mathbf{H}_2} = \|\mathbf{A}\|_{\mathbf{H}_2}^2 \leq \|\mathbf{A}^*\mathbf{A}\|_{\mathbf{B}} \leq C\|\mathbf{A}\|_{\mathbf{B}}^2$ , where  $C = \inf\{M \mid \|\mathbf{A}^*\mathbf{A}\|_{\mathbf{B}} \leq M\|\mathbf{A}\|_{\mathbf{B}}^2\}$ .

It should be noted that, in general,  $\|\mathbf{A}^*\mathbf{A}\|_{\mathbf{B}} \neq \|\mathbf{A}\|_{\mathbf{B}}^2$  and  $(\mathbf{A}\mathbf{B})^*x \neq \mathbf{B}^*\mathbf{A}^*x$ . Thus, as expected, there are some important differences compared to the corresponding operator results in Hilbert spaces. On the other hand, we can give a natural definition of orthogonality for subspaces of a Banach space.

**Definition 7.** Let **U** and **V** be subspaces of **B**. We say that **U** is orthogonal to **V** if,  $\forall x \in \mathbf{U}, \langle y, \varphi_x^s \rangle = 0 \quad \forall y \in \mathbf{V}.$ 

The above definition is transparent if we note that  $\langle y, \phi_x^s \rangle = 0 \quad \forall y \in \mathbf{V} \Leftrightarrow \langle y, J_2(x) \rangle = 0 \quad \forall y \in \mathbf{V}$ . The next result is easy to prove.

**Lemma 8.** If U is orthogonal to V, then V is orthogonal to U.

**Definition 9.** A biorthogonal system  $\{x_n, x_n^* | n \ge 1\}$  is called a Markushevich basis for **B** if the span of the  $x_n$  is dense in **B** and the span of the  $x_n^*$  is weak\* dense in **B**'.

Pelczynski [P] has shown that, for every separable Banach space **B** and each  $\epsilon > 0$ , **B** has a Markushevich basis such that  $||x_n|| ||x_n^*|| \le 1 + \epsilon$ . Diestel ([D], pg. 56) notes that the question of whether it is possible to require that  $||x_n|| = 1 = ||x_n^*||$  is open. In the next theorem, we show that, if **B** has a basis for a dense subspace, it has a Markushevich basis with unit norm.

**Theorem 10.** Let **B** be a separable Banach space with a basis for a dense subspace. If this basis is normalized and monotone with respect to the **B** norm, then **B** has a Markushevich basis  $\{x_n, x_n^* | n \ge 1\}$  such that  $\|x_n\|_{\mathbf{B}} = 1 = \|x_n^*\|_{\mathbf{B}'}$ .

*Proof.* (A basis is monotone if  $y = \sum a_i x_i$ , then  $\left\| \sum_{i=1}^m a_i x_i \right\|_{\mathbf{B}} \le \left\| \sum_{i=1}^{m+n} a_i x_i \right\|_{\mathbf{B}}$  for  $m, n \ge 1$ .) Let  $\{x_n | n \ge 1\}$  be a complete orthogonal basis for  $\mathbf{H}_1$  with  $\|x_n\|_{\mathbf{B}} = 1$ .

If we now define  $x_n^* = \varphi_n^s = \frac{\mathbf{J}_2(x_n)}{\|x_n\|_{\mathbf{H}_2}^2}$ , then it is easy to check that  $\langle x_i, x_j^* \rangle = \delta_{ij}$ . By definition, the span of the family  $\{x_n | n \geq 1\}$  is dense in  $\mathbf{B}$  and it is also easy to see that the span of the family  $\{x_n^*, n \geq 1\}$  is weak\* dense in  $\mathbf{B}'$ . To show that  $\|x_n^*\|_{\mathbf{B}}' = 1$ , let  $y = \sum_{i=1}^N a_i x_i$ ,  $\|y\|_B \leq 1$ , with  $N \geq 1$ . Then  $|\langle y, \varphi_n^s \rangle| \leq |a_n| \leq \|y\|_{\mathbf{B}} \leq 1$ , so that  $\|\varphi_n^s\|_B = \sup_{\|y\|_B \leq 1} |\langle y, \varphi_n^s \rangle| \leq 1$ . We are done since  $\langle x_n, \varphi_n^s \rangle = 1$ .

It is clear that much of the operator theory on Hilbert spaces can be extended to separable Banach spaces in a straightforward way. To get a flavor, we give a few of the more interesting results. Since the proofs are easy, we omit them. In what follows, all definitions are the same as in the case of a Hilbert space.

# Theorem 11. Let $A \in L[B]$ .

- (1) The set  $N(\mathbf{B})$  of all bounded normal operators on  $\mathbf{B}$  is a closed subset of  $L[\mathbf{B}]$ .
- (2) If  $\mathbf{A}$  is unitary on  $\mathbf{B}$ , then there exists a selfadjoint operator  $\mathbf{W}$ , and  $\mathbf{A} = \exp(i\mathbf{W})$ .

#### APPLICATION: THE YOSIDA APPROXIMATOR

If **A** is the generator of a strongly continuous semigroup  $T(t) = \exp(t\mathbf{A})$  on **B**, then the Yosida approximator for **A** is defined by  $\mathbf{A}_{\lambda} = \lambda \mathbf{A} R(\lambda, \mathbf{A})$ , where  $R(\lambda, \mathbf{A}) = (\lambda I - \mathbf{A})^{-1}$  is the resolvent of **A**. In general, **A** is closed and densely defined but unbounded. The Yosida approximator  $\mathbf{A}_{\lambda}$  is bounded, converges strongly to  $\mathbf{A}$ , and  $T_{\lambda}(t) = \exp(t\mathbf{A}_{\lambda})$  converges strongly to  $T(t) = \exp(t\mathbf{A})$ . If **A** generates a contraction semigroup, then so does  $\mathbf{A}_{\lambda}$  (see Pazy [Pz]). This result is very useful for applications. Unfortunately, for general semigroups, **A** may not have a bounded resolvent. Furthermore, it is very convenient to have a contractive approximator.

As an application of the theory in the previous section, we will show that the Yosida approach can be generalized in such a way as to give a contractive approximator for all strongly continuous semigroups of operators on B. For any closed densely defined linear operator A on B, let T = −[A∗A] 1/2 , T¯ <sup>=</sup> <sup>−</sup>[AA<sup>∗</sup> ] 1/2 . Since <sup>−</sup>T(−T¯) is maximal accretive, <sup>T</sup>(T¯) generates a contraction semigroup. We can now write A as A = UT, where U is a partial isometry (since the extension is valid on H2, the restriction is true on B). Define A<sup>λ</sup> by A<sup>λ</sup> = λAR(λ, T). Note that A<sup>λ</sup> = λUTR(λ, T) = λ <sup>2</sup>UR(λ, T) − λU and, although A does not commute with <sup>R</sup>(λ, <sup>T</sup>), we have <sup>λ</sup>AR(λ, <sup>T</sup>) = λR(λ, T¯)A.

<span id="page-7-0"></span>Theorem 12. *For every closed densely defined linear operator* A *on* B*, we have that*

- (1) A<sup>λ</sup> *is a bounded linear operator and* limλ→∞ Aλx = Ax, ∀x ∈ D(A),
- (2) exp[tAλ] *is a bounded contraction for* t > 0*, and*
- (3) *if* A *generates a strongly continuous semigroup* T (t) = exp[tA] *on* D *for* t > 0*,* D(A) ⊆ D, *then* limλ→∞ kexp[tAλ]x − exp[tA]xk<sup>B</sup> = 0 ∀x ∈ D.

*Proof.* : To prove 1, let <sup>x</sup> <sup>∈</sup> <sup>D</sup>(A). Now use the fact that limλ→∞ λR(λ, T¯)<sup>x</sup> <sup>=</sup> <sup>x</sup> and <sup>A</sup>λ<sup>x</sup> <sup>=</sup> λR(λ, T¯)Ax. To prove 2, use <sup>A</sup><sup>λ</sup> <sup>=</sup> <sup>λ</sup> <sup>2</sup>UR(λ, T) − λU, kλR(λ, T)k<sup>B</sup> = 1, and kUk<sup>B</sup> = 1 to get that kexp[tλ2UR(λ, T) − tλU]k<sup>B</sup> ≤ exp[−tλkUkB] exp[tλkUkBkλR(λ, T)kB] ≤ 1.

To prove 3, let t > 0 and  $x \in D(\mathbf{A})$ . Then

$$\|\exp[t\mathbf{A}]x - \exp[t\mathbf{A}_{\lambda}]x\|_{\mathbf{B}} = \|\int_{0}^{t} \frac{d}{ds} [e^{(t-s)\mathbf{A}_{\lambda}} e^{s\mathbf{A}}]x ds\|_{\mathbf{B}}$$

$$\leq \int_{0}^{t} \|[e^{(t-s)\mathbf{A}_{\lambda}} (\mathbf{A} - \mathbf{A}_{\lambda}) e^{s\mathbf{A}}x]\|_{\mathbf{B}}$$

$$\leq \int_{0}^{t} \|[(\mathbf{A} - \mathbf{A}_{\lambda}) e^{s\mathbf{A}}x]\|_{\mathbf{B}} ds.$$

Now use  $\|[\mathbf{A}_{\lambda}e^{s\mathbf{A}}x]\|_{\mathbf{B}} = \|[\lambda R(\lambda, \mathbf{\bar{T}})e^{s\mathbf{A}}\mathbf{A}x]\|_{\mathbf{B}} \leq \|[e^{s\mathbf{A}}\mathbf{A}x]\|_{\mathbf{B}}$  to get  $\|[(\mathbf{A} - \mathbf{A}_{\lambda})e^{s\mathbf{A}}x]\|_{\mathbf{B}} \leq 2\|[e^{s\mathbf{A}}\mathbf{A}x]\|_{\mathbf{B}}$ . Now, since  $\|[e^{s\mathbf{A}}\mathbf{A}x]\|_{\mathbf{B}}$  is continuous, by the bounded convergence theorem we have  $\lim_{\lambda \to \infty} \|\exp[t\mathbf{A}]x - \exp[t\mathbf{A}_{\lambda}]x\|_{\mathbf{B}} \leq \int_0^t \lim_{\lambda \to \infty} \|[(\mathbf{A} - \mathbf{A}_{\lambda})e^{s\mathbf{A}}x]\|_{\mathbf{B}}ds = 0.$ 

#### CONCLUSION

The first part of Theorem 12 is a generalization of a result of Kaufman [Ka]. This allows us to provide a new metric for closed densely defined linear operators on Banach spaces. If A, B are closed and densely defined, we can define our metric by  $d(A, B) = ||A_0 - B_0||$ ,  $A_0 = A(1 + A^*A)^{-\frac{1}{2}}$ ,  $B_0 = B(1 + B^*B)^{-\frac{1}{2}}$ . The Hille-Yosida Theorem for contraction semigroups gives necessary and sufficient conditions for a closed densely defined linear operator to be a generator. The general strongly continuous case may be reduced to the contraction case by shifting the spectrum and using an equivalent norm. The second part of Theorem 12 may be viewed as an improvement in the sense that, by using the approximator, this procedure is no longer required.

## References

- <span id="page-9-6"></span><span id="page-9-0"></span>[D] J. Diestel, Sequences and Series in Banach Spaces, GTM 92, Springer-Verlag, (1984).
- <span id="page-9-4"></span>[G] L. Gross, Abstract Wiener spaces, Proc. Fifth Berkeley Symposium on Mathematics Statistics and Probability, (1965), 31–42.
- <span id="page-9-1"></span>[Ka] W. E. Kaufman, A Stronger Metric for Closed Operators in Hilbert Spaces, Proc. Amer. Math. Soc. 90 (1984), 83–87.
- <span id="page-9-3"></span>[K] J. Kuelbs, Gaussian measures on a Banach Space, Journal of Functional Analysis 5 (1970), 354–367.
- <span id="page-9-7"></span>[L] P. D. Lax, Symmetrizable Linear Tranformations, Comm. Pure Appl. Math. 7 (1954), 633–647.
- <span id="page-9-5"></span>[Pz] A. Pazy, Semigroups of linear operators and applications to Partial Differential Equations. Applied Mathematical Sciences, 44, Springer New York, (1983).
- <span id="page-9-2"></span>[P] A. Pelczynski, All Separable Banach Spaces admit for ε > 0 fundamental and total biorthogonal sequences bounded by 1 + ε, Studia Math. 55 (1976), 295–304.
- [VN] J. von Neumann, Uber adjungierte Funktionaloperatoren, Annals of Mathematics 33 (1932), 294–310.

(Tepper L. Gill) Department of Electrical Engineering Howard University, Washington DC 20059, USA, E-mail : tgill@howard.edu

(Sudeshna Basu) Department of Mathematics, Howard University, Washington DC 20059, USA, E-mail : sbasu@howard.edu

(Woodford W. Zachary) Department of Electrical Engineering, Howard University, Washington DC 20059, USA, E-mail : wwzachary@earthlink.net

(V. Steadman) Department of Mathematics, University of Distrcit of Columbia, Washington DC