## TEPPER L. GILL

Abstract. This note has two objectives. The first objective is show that, even if a separable Banach space does not have a Schauder basis (Sbasis), there always exists Hilbert spaces H<sup>1</sup> and H2, such that H<sup>1</sup> is a continuous dense embedding in B and B is a continuous dense embedding in H2. This is the best possible improvement of a theorem due to Mazur (see [\[BA\]](#page-10-0) and also [\[PE1\]](#page-12-0)). The second objective is show how H<sup>2</sup> allows us to provide a positive answer to the Marcinkiewicz-basis (M-basis) problem.

## 1. Introduction

Definition 1.1. *Let* B *separable Banach space, with dual space* B ∗ *. A sequence* (xn) ∈ B *is called a S-basis for* B *if* kxnk<sup>B</sup> = 1 *and, for each* x ∈ B*, there is a unique sequence* (an) *of scalars such that*

$$x = \lim_{n \to \infty} \sum_{k=1}^{n} a_k x_k = \sum_{k=1}^{\infty} a_k x_k.$$

1991 Mathematics Subject Classification. Primary (46B03), (46B20) Secondary(46B25).

Key words and phrases. Marcinkiewicz basis, Schauder basis, biorthogonal system, duality mappings, Banach spaces.

1

Definition 1.2. *Let* h{x<sup>i</sup> : <sup>i</sup> <sup>∈</sup> <sup>N</sup>}i *be the set of all linear combinations of the family of vectors* {xi} *(linear span). The family* {(x<sup>i</sup> , x<sup>∗</sup> i )} ∞ <sup>i</sup>=1 ⊂ B × B<sup>∗</sup> *is called:*

- (1) *A* fundamental *system if* h{x<sup>i</sup> : <sup>i</sup> <sup>∈</sup> <sup>N</sup>}i <sup>=</sup> <sup>B</sup>*.*
- (2) *A* minimal *system if* x<sup>j</sup> ∈/ h{x<sup>i</sup> : <sup>i</sup> <sup>∈</sup> <sup>N</sup> \ {j}}i*.*
- (3) *<sup>A</sup>* total *if for each* <sup>x</sup> <sup>6</sup>= 0 *there exists* <sup>i</sup> <sup>∈</sup> <sup>N</sup> *such that* <sup>x</sup> ∗ i (x) 6= 0*.*
- (4) *A* biorthogonal *system if* x ∗ i (x<sup>j</sup> ) = <sup>δ</sup>ij *, for all* i, j <sup>∈</sup> <sup>N</sup>*.*
- (5) *A* M-basis *if it is a fundamental minimal, total and biorthogonal system.*

The first problem we consider had its beginning with a question raised by Banach. He asked whether every separable Banach space has a S-basis. Mazur gave a partial answer. He proved that every infinite-dimensional separable Banach space contains an infinite-dimensional subspace with a S-basis.

In 1972, Enflo [\[EN\]](#page-11-0) answered Banach's question in the negative by providing a separable Banach space B, without a S-basis and without the approximation property (i.e., every compact operator on B is the limit of a sequence of finite rank operators). Every Banach space with a S-basis has the approximation property and Grothendieck [\[GR\]](#page-11-1) proved that if a Banach space had the approximation property, then it would also have a S-basis. In the first section we show that, given B there exists separable Hilbert spaces H<sup>1</sup> and H<sup>2</sup> such that H<sup>1</sup> ⊂ B ⊂ H<sup>2</sup> as continuous dense embeddings. The

existence of H<sup>1</sup> is the best possible improvement of Mazur's Theorem, while the existence of H<sup>2</sup> shows that B is very close to the best possible case in a well-defined manner.

The second problem we consider is associated with a weaker structure discovered by Marcinkiewicz [\[M\]](#page-11-2). He showed that every separable Banach space <sup>B</sup> has a biorthogonal system {xn, x<sup>∗</sup> <sup>n</sup>}, with h{xn}i = B. This system has many of the properties of an S-basis and is now known as a M-basis for B. A well-known open problem for the M-basis is whether one can choose the system {xn, x<sup>∗</sup> <sup>n</sup>} such that kxnk kx ∗ <sup>n</sup>k = 1 (see Diestel [\[D\]](#page-10-1)). This is called the M-basis problem for separable Banach spaces. It has been studied by Singer [\[SI\]](#page-12-1), Davis and Johnson [\[DJ\]](#page-10-2), Ovsepian and Pelczy´niski [\[OP\]](#page-12-2), Pelczy´niski [\[PE\]](#page-12-3) and Plichko [\[PL\]](#page-12-4). The work of Ovsepian and Pelczy´niski [\[OP\]](#page-12-2) led to the construction of a bounded M-basis, while that of Pelczy´niski [\[PE\]](#page-12-3) and Plichko [\[PL\]](#page-12-4) led to independent proofs that, for every ε > 0, it is possible to find a biorthogonal system with the property that kxnk kx ∗ <sup>n</sup>k < 1 + ε. The question of whether we can set ε = 0 has remained unanswered since 1976. In this case, we provide a positive answer by constructing a biorthogonal system with the property that kxnk kx ∗ n k = 1.

# 2. The S-basis Problem

In this section, we construct our Hilbert space rigging of any given separable Banach space as continuous dense embeddings. We begin with the construction of H2.

**Theorem 2.1.** Suppose  $\mathcal{B}$  is a separable Banach space, then there exist a separable Hilbert space  $\mathcal{H}_2$  such that,  $\mathcal{B} \subset \mathcal{H}_2$  as a continuous dense embedding.

*Proof.* Let  $\{x_n\}$  be a countable dense sequence in  $\mathcal{B}$  and let  $\{x_n^*\}$  be any fixed set of corresponding duality mappings (i.e.,  $x_n^* \in \mathcal{B}^*$ , the dual space of  $\mathcal{B}$  and  $x_n^*(x_n) = \langle x_n, x_n^* \rangle = \|x_n\|_{\mathcal{B}}^2 = \|x_n^*\|_{\mathcal{B}^*}^2$ ). For each n, let  $t_n = \frac{1}{\|x_n^*\|^2 2^n}$  and define (u, v) by:

$$(u,v) = \sum_{n=1}^{\infty} t_n x_n^*(u) \bar{x}_n^*(v) = \sum_{n=1}^{\infty} \frac{1}{\|x_n^*\|^2 2^n} x_n^*(u) \bar{x}_n^*(v).$$

It is easy to see that (u, v) is an inner product on  $\mathcal{B}$ . Let  $\mathcal{H}_2$  be the completion of  $\mathcal{B}$  with respect to this inner product. It is clear that  $\mathcal{B}$  is dense in  $\mathcal{H}_2$ , and

$$||u||_{\mathcal{H}_2}^2 = \sum_{n=1}^{\infty} t_n |x_n^*(u)|^2 \le \sup_{n} \frac{1}{||x_n^*||^2} |x_n^*(u)|^2 = ||u||_{\mathcal{B}}^2,$$

so the embedding is continuous.

In order to construct our second Hilbert space, we need the following result by Lax [L].

**Theorem 2.2** (Lax). Let  $A \in L[\mathcal{B}]$ . If A is selfadjoint on  $\mathcal{H}_2$  (i.e.,  $(Ax,y)_{\mathcal{H}_2} = (x,Ay)_{\mathcal{H}_2}, \forall x,y \in \mathcal{B}$ ), then A has a bounded extension to  $\mathcal{H}_2$  and  $\|A\|_{\mathcal{H}_2} \leq M \|A\|_{\mathcal{B}}$  for some positive constant M.

*Proof.* Let  $x \in \mathcal{B}$  and, without loss, we can assume that M = 1 and  $||x||_{\mathcal{H}_2} = 1$ . Since A is selfadjoint,

$$||Ax||_{\mathcal{H}_2}^2 = (Ax, Ax) = (x, A^2x) \leqslant ||x||_{\mathcal{H}_2} ||A^2x||_{\mathcal{H}_2} = ||A^2x||_{\mathcal{H}_2}.$$

Thus, we have  $||Ax||_{\mathcal{H}_2}^4 \leq ||A^4x||_{\mathcal{H}_2}$ , so it is easy to see that  $||Ax||_{\mathcal{H}_2}^{2n} \leq ||A^{2n}x||_{\mathcal{H}_2}$  for all n. It follows that:

$$||Ax||_{\mathcal{H}_2} \leqslant (||A^{2n}x||_{\mathcal{H}_2})^{1/2n} \leqslant (||A^{2n}x||_{\mathcal{B}})^{1/2n}$$
  
$$\leqslant (||A^{2n}||_{\mathcal{B}})^{1/2n} (||x||_{\mathcal{B}})^{1/2n} \leqslant ||A||_{\mathcal{B}} (||x||_{\mathcal{B}})^{1/2n}.$$

Letting  $n \to \infty$ , we get that  $||Ax||_{\mathcal{H}_2} \le ||A||_{\mathcal{B}}$  for any x in the dense set of the unit ball  $B_{\mathcal{H}_2} \cap \mathcal{B}$ . Since the norm is attained on a dense set of the unit ball, we are done.

For our second Hilbert space, fix  $\mathcal{B}$  and define  $\mathcal{H}_1$  by:

$$\mathcal{H}_{1} = \left\{ u \in \mathcal{B} \left| \sum_{n=1}^{\infty} t_{n}^{-1} \left| (u, x_{n})_{2} \right|^{2} < \infty \right\}, \text{ with} \right.$$

$$\left. (u, v)_{1} = \sum_{n=1}^{\infty} t_{n}^{-1} \left( u, x_{n} \right)_{2} (x_{n}, v)_{2}. \right.$$

For  $u \in \mathcal{B}$ , let  $T_{12}u$  be defined by  $T_{12}u = \sum_{n=1}^{\infty} t_n (u, x_n)_2 x_n$ .

**Theorem 2.3.** The operator  $T_{12}$  is a positive trace class operator on  $\mathcal{B}$  with a bounded extension to  $\mathcal{H}_2$ . In addition,  $\mathcal{H}_1 \subset \mathcal{B} \subset \mathcal{H}_2$  (as continuous dense embeddings),  $\left(T_{12}^{1/2}u, T_{12}^{1/2}v\right)_1 = (u, v)_2$  and  $\left(T_{12}^{-1/2}u, T_{12}^{-1/2}v\right)_2 = (u, v)_1$ .

*Proof.* First, since terms of the form  $\{u_N = \sum_{k=1}^N t_k^{-1} (u, x_k)_2 x_k : u \in \mathcal{B}\}$  are dense in  $\mathcal{B}$ , we see that  $\mathcal{H}_1$  is dense in  $\mathcal{B}$ . It follows that  $\mathcal{H}_1$  is also dense in  $\mathcal{H}_2$ .

For the operator  $T_{12}$ , we see that  $\mathcal{B} \subset \mathcal{H}_2 \Rightarrow (u, x_n)_2$  is defined for all  $u \in \mathcal{B}$ , so that  $T_{12}$  maps  $\mathcal{B} \to \mathcal{B}$  and:

$$||T_{12}u||_{\mathcal{B}}^{2} \leq \left[\sum_{n=1}^{\infty} t_{n}^{2} ||x_{n}||_{\mathcal{B}}^{2}\right] \left[\sum_{n=1}^{\infty} |(u, x_{n})_{2}|^{2}\right] = M ||u||_{2}^{2} \leq M ||u||_{\mathcal{B}}^{2}.$$

Thus,  $T_{12}$  is a bounded operator on  $\mathcal{B}$ . It is clearly trace class and, since  $(T_{12}u, u)_2 = \sum_{n=1}^{\infty} t_n |(u, x_n)_2|^2 > 0$ , it is positive. From here, it's easy to see that  $T_{12}$  is selfadjoint on  $\mathcal{H}_2$  so, by Theorem 2.2 it has a bounded extension to  $\mathcal{H}_2$ .

An easy calculation now shows that 
$$\left(T_{12}^{1/2}u, T_{12}^{1/2}v\right)_1 = (u, v)_2$$
 and  $\left(T_{12}^{-1/2}u, T_{12}^{-1/2}v\right)_2 = (u, v)_1$ .

Since the counter example of Enflo, the only direct information about a Banach space without a basis has been the following theorem of Mazur:

**Theorem 2.4.** Every infinite dimensional separable Banach contains a infinite dimensional subspace with a basis.

Theorems 2.1 and 2.3 show that, even if a Banach space does not have a basis, it is very close to the best possible case.

Remark 2.5. Historically, Gross [G] first proved that every real separable Banach space  $\mathcal{B}$  contains a separable Hilbert space (version of  $\mathcal{H}_1$ ), as a dense embedding, and that this space is the support of a Gaussian measure. Then Kuelbs [KB] showed that one can construct  $\mathcal{H}_2$  so that  $\mathcal{H}_1 \subset \mathcal{B} \subset \mathcal{H}_2$  as continuous dense embeddings, with  $\mathcal{H}_1$  and  $\mathcal{H}_2$  related by Theorem 2.3.

*A particular Gross-Kuelbs construction of* H<sup>2</sup> *was used in* [\[GZ\]](#page-11-6) *to provide the foundations for the Feynman path integral formulation of quantum mechanics* [\[FH\]](#page-11-7) *(see also* [\[GZ1\]](#page-11-8)*).*

*This construction was also used in* [\[GBZS\]](#page-11-9) *to show that every bounded linear operator* <sup>A</sup> *on a separable Banach space* <sup>B</sup> *has a adjoint* <sup>A</sup><sup>∗</sup> *defined on* B*, such that (see below):*

- (1) <sup>A</sup>∗<sup>A</sup> *is m-accretive (i.e., if* <sup>x</sup> ∈ B *and* <sup>x</sup> ∗ *is a corresponding duality, then* <sup>h</sup>A∗Ax, x<sup>∗</sup> i ≥ 0*),*
- (2) (A∗A) <sup>∗</sup> = A∗A *(selfadjoint), and*
- (3) I + A∗A *has a bounded inverse.*

Example 2.6. *The following example shows how easy it is to construct an adjoint* A<sup>∗</sup> *satisfying all the above conditions, using only* H2*. Let* Ω *be a bounded open domain in* R <sup>n</sup> *with a class* C 1 *boundary and let* <sup>H</sup><sup>1</sup> 0 [Ω]*, be the set of all real-valued functions* u ∈ L 2 [Ω] *such that their first order weak partial derivatives are in* L 2 [Ω] *and vanish on the boundary. It follows that*

$$(u, v) = \int_{\Omega} \nabla u(\mathbf{x}) \cdot \nabla v(\mathbf{x}) d\mathbf{x} = \langle u, J_0 v \rangle,$$

*defines an inner product on* <sup>H</sup><sup>1</sup> 0 [Ω]*, where* J0 *is the conjugate isomorphism between* <sup>H</sup><sup>1</sup> 0 [Ω] *and its dual* <sup>H</sup>−<sup>1</sup> [Ω]*. The space* <sup>H</sup>−<sup>1</sup> [Ω] *coincides with the set of all distributions of the form*

$$u = h_0 + \sum_{i=1}^n \frac{\partial h_i}{\partial x_i}, \quad \text{where } h_i \in L^2[\Omega], \quad 1 \leqslant i \leqslant n.$$

In this case we also have for  $p \in [2, \infty)$  and  $q \in (1, 2], \frac{1}{p} + \frac{1}{q} = 1$  that,

$$\mathcal{H}_0^1[\Omega] \subset L^p[\Omega] \subset L^q[\Omega] \subset \mathcal{H}^{-1}[\Omega]$$

all as continuous dense embeddings.

From the inner product on  $\mathcal{H}_0^1[\Omega]$  we see that  $J_0 = -\Delta$ , the Laplace operator under Dirichlet homogeneous boundary conditions on  $\Omega$ . If we set  $\mathcal{H}_1 = \mathcal{H}_0^1[\Omega]$ ,  $\mathcal{H}_2 = \mathcal{H}^{-1}$  and  $J = J_0^{-1}$ , then for every  $A \in \mathcal{C}[L^p(\Omega)]$  (i.e., the closed densely defined linear operators on  $L^p(\Omega)$ ), we obtain  $A^* \in \mathcal{C}[L^p(\Omega)]$ , where  $A^* = J^{-1}A'J|_p = [-\Delta]A'[-\Delta]^{-1}|_p$  for each  $A' \in \mathcal{C}[L^q(\Omega)]$ . It is now easy to show that  $A^*$  satisfies the conditions (1)-(3) above for an adjoint operator on  $L^p(\Omega)$ .

## 3. The M-basis Problem

To understand the M-basis problem and its solution in a well-known setting, let  $\mathbb{R}^2$  have its standard inner product  $(\cdot, \cdot)$  and let  $x_1, x_2$  be any two independent basis vectors. Define a new inner product on  $\mathbb{R}^2$  by

(3.1) 
$$\langle y \mid z \rangle = t_1 (x_1 \otimes x_1) (y \otimes z) + t_2 (x_2 \otimes x_2) (y \otimes z)$$
$$= t_1 (y, x_1) (z, x_1) + t_2 (y, x_2) (z, x_2) ,$$

where  $t_1, t_2 > 0, t_1 + t_2 = 1$ . Define new functionals  $S_1$  and  $S_2$  by:

$$S_1(x) = \frac{\langle x \mid x_1 \rangle}{\alpha_1 \langle x_1 \mid x_1 \rangle}, \quad S_2(x) = \frac{\langle x \mid x_2 \rangle}{\alpha_2 \langle x_2 \mid x_2 \rangle}, \quad \text{for} \quad y \in \mathbb{R}^2.$$

Where α1, α<sup>2</sup> > 0 are chosen to ensure that kS1k = kS2k = 1. Note that, if (x1, x2) = 0, S1 and S2 reduce to

$$S_1(x) = \frac{(x, x_1)}{\alpha_1 \|x_1\|}, \quad S_2(x) = \frac{(x, x_2)}{\alpha_2 \|x_2\|}.$$

Thus, we can define many equivalent inner products on R <sup>2</sup> and many linear functionals with the same properties but different norms.

The following example shows how this construction can be of use.

Example 3.1. *In this example, let* x1 = e1 *and* x2 = e1 + e2*, where* e1 = (1, 0), e2 = (0, 1)*. In this case, the biorthogonal functionals are generated by the vectors* x¯<sup>1</sup> = e1−e<sup>2</sup> *and* x¯<sup>2</sup> = e<sup>2</sup> (i.e., x ∗ 1 (x) = (x, x¯1), x<sup>∗</sup> 2 (x) = (x, x¯2))*. It follows that* (x1, x¯2) = 0, (x1, x¯1) = 1 *and* (x2, x¯1) = 0, (x2, x¯2) = 1*. However,* kx1k kx¯1k = √ 2, kx2k kx¯2k = √ 2*, so that* {x1,( ·, x¯1)} *and* {x2,( ·, <sup>x</sup>¯2)} *fails to solve the M-basis problem on* <sup>R</sup> 2 *.*

*In this case, we set* α<sup>1</sup> = 1 *and* α<sup>2</sup> = kx2k *so that, without changing* x<sup>1</sup> *and* x2*, and using the inner product from equation (1.1) in the form*

$$\langle x \mid y \rangle = t_1(x, \bar{x}_1)(y, \bar{x}_1) + t_2(x, \bar{x}_2)(y, \bar{x}_2),$$

S1 *and* S2 *become*

$$S_1(x) = \frac{(x, \bar{x}_1)}{\|\bar{x}_1\|}, \quad S_2(x) = \frac{(x, \bar{x}_2)}{\|x_2\|}.$$

*It now follows that* Si(xi) = 1 *and* Si(x<sup>j</sup> ) = 0 *for* i 6= j *and* kS<sup>i</sup> k kx<sup>i</sup> k = 1*, so that system* {x1, S1} *and* {x2, S2} *solves the M-basis problem.*

Remark 3.2. For a given set of independent vectors on a finite dimensional vector space, It is known that the corresponding biorthogonal functionals are unique. This example shows that uniqueness is only up to a scale factor and this is what we need to produce an M-basis.

The following theorem shows how our solution to the M-basis problem for  $\mathbb{R}^2$  can be extended to any separable Banach space.

**Theorem 3.3.** Let  $\mathcal{B}$  be a infinite-dimensional separable Banach space. Then  $\mathcal{B}$  contains an M-basis with the property that  $\|x_i\|_{\mathcal{B}} \|x_i^*\|_{\mathcal{B}^*} = 1$  for all i.

Proof. Construct  $\mathcal{H} = \mathcal{H}_2$  via Theorem 2.1, so that  $\mathcal{B} \subset \mathcal{H}$  is a dense continuous embedding and let  $\{x_i\}_{i=1}^{\infty}$  be a fundamental minimal system for  $\mathcal{B}$ . If  $i \in \mathbb{N}$ , let  $M_{i,\mathcal{H}}$  be the closure of the span of  $\{x_i\}$  in  $\mathcal{H}$ . Thus,  $x_i \notin M_{i,\mathcal{H}}^{\perp}$ ,  $M_{i,\mathcal{H}} \oplus M_{i,\mathcal{H}}^{\perp} = \mathcal{H}$  and  $(y,x_i)_{\mathcal{H}} = 0$  for all  $y \in M_{i,\mathcal{H}}^{\perp}$ .

Let  $\hat{M}_i$  be the closure of the span of  $\{x_j \ j \neq i\}$  in  $\mathcal{B}$ . Since  $\hat{M}_i \subset M_{i,\mathcal{H}}^{\perp}$  and  $x_i \notin \hat{M}_i$ ,  $(y, x_i)_{\mathcal{H}} = 0$  for all  $y \in \hat{M}_i$ . Let the seminorm  $p_i(\cdot)$  be defined on the closure of the span of  $\{x_i\}$ , in  $\mathcal{B}$  by  $p_i(y) = \|x_i\|_{\mathcal{B}} \|y\|_{\mathcal{B}}$ , and define  $\hat{x}_i^*(\cdot)$  by:

$$\hat{x}_{i}^{*}(y) = \frac{\|x_{i}\|_{\mathcal{B}}^{2}}{\|x_{i}\|_{\mathcal{H}}^{2}} (y, x_{i})_{\mathcal{H}}$$

By the Hahn-Banach Theorem,  $\hat{x}_i^*(\cdot)$  has an extension  $x_i^*(\cdot)$  to  $\mathcal{B}$ , such that  $|x_i^*(y)| \leq p_i(y) = ||x_i||_B ||y||_B$  for all  $y \in \mathcal{B}$ . By definition of  $p_i(\cdot)$ , we see that  $||x_i^*||_{\mathcal{B}^*} \leq ||x_i||_{\mathcal{B}}$ . On the other hand  $x_i^*(x_i) = ||x_i||_{\mathcal{B}}^2 \leq ||x_i||_{\mathcal{B}} ||x_i^*||_{\mathcal{B}^*}$ ,

so that  $x_i^*(\cdot)$  is a duality mapping for  $x_i$ . If  $x_i^*(x) = 0$  for all i, then  $x \in \bigcap_{i=1}^{\infty} \hat{M}_i = \{0\}$  so that the family  $\{x_i^*\}_{i=1}^{\infty}$  is total. If we let  $\|x_i\|_{\mathcal{B}} = 1$ , it is clear that  $x_i^*(x_j) = \delta_{ij}$ , for all  $i, j \in \mathbb{N}$ . Thus,  $\{x_i, x_i^*\}$  is an M-basis system with  $\|x_i\|_{\mathcal{B}} \|x_i^*\|_{\mathcal{B}^*} = 1$  for all i.

### Conclusion

In this paper we have first shown that every infinite dimensional separable Banach space is very close to a Hilbert space in a well defined manner, providing the best possible improvement on the well-known theorem of Mazur. We have then provided a solution to the M-basis problem by showing that every infinite dimensional separable Banach space has a M-basis  $\{x_i, x_i^*\}$ , with the property that  $\|x_i\|_{\mathcal{B}} \|x_i^*\|_{\mathcal{B}^*} = 1$  for all i.

## REFERENCES

- <span id="page-10-0"></span>[BA] S. Banach Théorie des Opérations linéaires, Monografj Matematyczn, Vol. 1, Warsaw, (1932).
- [BM] S. Banach and S. Mazur, Zur Theorie der linearen Dimension, Studia Mathematica, 4 (1933), 100-110.
- <span id="page-10-1"></span>[D] J. Diestel, Sequences and Series in Banach Spaces, Grad. Texts in Math. Springer-Verlag, New York, (1984).
- <span id="page-10-2"></span>[DJ] W. J. Davis and W. B. Johnson, On the existence of fundamental and total bounded biorthogonal systems in Banach spaces, Studia Math. 45 (1973), 173-179.

- <span id="page-11-0"></span>[EN] P. Enflo, *A counterexample to the approximation problem in Banach spaces*, Acta Mathematica 130 (1), (1973) 309317.
- <span id="page-11-7"></span>[FH] R. P. Feynman and A. R. Hibbs, *Quantum Mechanics and Path Integrals*, McGraw-Hill, New York, (1965).
- <span id="page-11-4"></span>[G] L. Gross, *Abstract Wiener spaces,* Proc. Fifth Berkeley Symposium on Mathematics Statistics and Probability, (1965), 31-42.
- <span id="page-11-9"></span>[GBZS] T. Gill, S. Basu, W. W. Zachary and V. Steadman, *Adjoint for operators in Banach spaces*, Proceedings of the American Mathematical Society, 132 (2004), 1429-1434.
- <span id="page-11-1"></span>[GR] A. Grothendieck, *Produits tensoriels topologiques et espaces nucleaires*, Memo. Amer. Math. Soc. 16 (1955).
- <span id="page-11-6"></span>[GZ] T. L. Gill and W. W. Zachary, *A New Class of Banach Spaces*, Journal of Physics A: Math. and Gen. 41 (2008), 495206.
- <span id="page-11-8"></span>[GZ1] T. L. Gill and W. W. Zachary, *Banach Spaces for the Feynman integral*, Real Analysis Exchange 34(2) (2008)/(2009), 267-310.
- <span id="page-11-5"></span>[KB] J. Kuelbs, *Gaussian measures on a Banach space,* Journal of Functional Analysis 5 (1970), 354-367.
- <span id="page-11-3"></span>[L] P. D. Lax, Symmetrizable linear tranformations. *Comm. Pure Appl. Math.* 7 (1954), 633–647.
- <span id="page-11-2"></span>[M] A. I. Markusevich, *On a basis in the wide sense for linear spaces,* Dokl. Akad. Nauk. SSSR, 41, (1943), 241-244.
- [MA] A. I. Markusevich, *On a basis in the wide sense for linear spaces,* Dokl. Akad. Nauk. SSSR, 41, (1943), 241-244.

- <span id="page-12-2"></span>[OP] R. I.Ovsepian and A. Pelczy´niski, *The existence in separable Banach space of fundamental total and bounded biorthogonal sequence and related constructions of uniformly bounded orthonormal systems in* L 2 , Studia Math. 54 (1975), 149-159.
- <span id="page-12-3"></span>[PE] A. Pelczy´niski, *All separable Banach spaces admit for every* ε > 0 *fundamental total and bounded by* 1 + ε *biorthogonal sequences,* Studia Math., 55, (1976), 295-304.
- <span id="page-12-0"></span>[PE1] A. Pelczy´niski, *A note on the paper of Singer "Basic sequences and reflexivity of Banach spaces",* Studia Math., 21, (1962), 371-374.
- <span id="page-12-4"></span>[PL] A. N. Plichko, *M-bases in separable and refexive Banach spaces*, Ukrain. Mat. Z. 29 (1977), 681-685.
- <span id="page-12-1"></span>[RU] W. Rudin, *Functional Analysis*, McGraw-Hill Press, New York, (1973).
- [SI] I. Singer, *On biorthogonal systems and total sequences of functionals II*, Math. Ann. 201 (1973), 1-8.
- (Tepper L. Gill) Departments of Electrical & Computer Engineering and Mathematics Howard University, Washington DC 20059, USA, E-mail : tgill@howard.edu