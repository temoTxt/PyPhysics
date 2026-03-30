### SOME BANACH SPACES ARE ALMOST HILBERT

TEPPER L. GILL<sup>1\*</sup> AND MARZETT GOLDEN<sup>1</sup>

ABSTRACT. The purpose of this note is to show that, if  $\mathcal{B}$  is a uniformly convex Banach, then the dual space  $\mathcal{B}'$  has a "Hilbert space representation" (defined in the paper), that makes  $\mathcal{B}$  much closer to a Hilbert space then previously suspected. As an application, we prove that, if  $\mathcal{B}$  also has a Schauder basis (Sbasis), then for each  $A \in \mathbb{C}[\mathcal{B}]$  (the closed and densely defined linear operators), there exists a closed densely defined linear operator  $A^* \in \mathbb{C}[\mathcal{B}]$  that has all the expected properties of an adjoint. Thus for example, the bounded linear operators,  $L[\mathcal{B}]$ , is a \*algebra. This result allows us to give a natural definition to the Schatten class of operators on a uniformly convex Banach space with a S-basis. In particular, every theorem that is true for the Schatten class on a Hilbert space, is also true on such a space. The main tool we use is a special version of a result due to Kuelbs [K], which shows that every uniformly convex Banach space with a S-basis can be densely and continuously embedded into a Hilbert space which is unique up to a change of basis.

#### 1. Introduction

In 1965, Gross [G] proved that every real separable Banach space contains a separable Hilbert space as a dense embedding, and this space is the support of a Gaussian measure. This was a generalization of Wiener's theory, that was based on the use of the (densely embedded Hilbert) Sobolev space  $\mathbb{H}^1[0,1] \subset \mathbb{C}[0,1]$ . In 1972, Kuelbs [K] generalized Gross' theorem to include the Hilbert space rigging  $\mathbb{H}^1[0,1] \subset \mathbb{C}[0,1] \subset L^2[0,1]$ . This general theorem can be stated as:

**Theorem 1.1.** (Gross-Kuelbs) Let  $\mathcal{B}$  be a separable Banach space. Then there exist separable Hilbert spaces  $\mathcal{H}_1, \mathcal{H}_2$  and a positive trace class operator  $T_{12}$  defined on  $\mathcal{H}_2$  such that  $\mathcal{H}_1 \subset \mathcal{B} \subset \mathcal{H}_2$  all as continuous dense embeddings, with  $\left(T_{12}^{1/2}u, T_{12}^{1/2}v\right)_1 = (u, v)_2$  and  $\left(T_{12}^{-1/2}u, T_{12}^{-1/2}v\right)_2 = (u, v)_1$ .

This theorem makes it possible to give a definition of the adjoint for bounded linear operators on separable Banach spaces. The definition has all the expected properties. In particular, It can be shown that, for each bounded linear operator A on  $\mathcal{B}$ , there exists  $A^*$ , with  $A^*A$ , maximal accretive, selfadjoint,  $(A^*A)^* = A^*A$ , and  $I + A^*A$  is invertible (see [GBZS]).

The basic idea is simple, let A be bounded on  $\mathcal{B}$  and let  $A_1$  be the restriction of A to  $\mathcal{H}_1$ . We can now consider  $A_1: \mathcal{H}_1 \to \mathcal{H}_2$ . If  $\mathbf{J}_2: \mathcal{H}_2 \to \mathcal{H}_2'$  is the standard conjugate isomorphism, then  $(A_1')\mathbf{J}_2: \mathcal{H}_1 \to \mathcal{H}_2'$ , so that  $\mathbf{J}_1^{-1}(A_1')\mathbf{J}_2: \mathcal{H}_1 \to \mathcal{H}_1 \subset$ 

Date: Received: xxxxxx; Revised: yyyyyy; Accepted: zzzzzz.

<sup>\*</sup> Corresponding author.

<sup>2010</sup> Mathematics Subject Classification. Primary 47B37; Secondary 46B10; 46C99. Key words and phrases. adjoint on Banach space, Lax Theorem; Schatten Classes.

 $\mathcal{B}$ . It follows that  $\mathbf{J}_1^{-1}(A_1')\mathbf{J}_2|_{\mathcal{B}}: \mathcal{B} \to \mathcal{B}$ . It easy to show that  $A^* = \mathbf{J}_1^{-1}(A_1')\mathbf{J}_2|_{\mathcal{B}}$  has the the main properties of an adjoint for A on  $\mathcal{B}$ .

At this level of generality, the definition of an adjoint for closed operators depends on the domain of the operator and changes the choice of  $\mathcal{H}_1$  (and  $T_{12}$ ) for each operator. Thus, the adjoint is only reasonable for bounded operators. It is a open question if all closed densely defined linear operators can have an adjoint with all the expected properties. Part of the problem is that not all operators in  $\mathcal{C}[\mathcal{B}]$  are of Baire class one (for example, when  $\mathcal{B}$  is nonreflexive). An operator A is of Baire class one if and only if it can be approximated by a sequence,  $\{A_n\}$ , of bounded linear operators. The solution is unknown and we suspect that, in general, there may be at least one operator in  $\mathcal{C}[\mathcal{B}]$  without an adjoint.

- 1.1. **Purpose.** In this note, we focus on uniformly convex Banach spaces, the best class of spaces that are not Hilbert. Our purpose is to show that these spaces are very close to Hilbert spaces and give the best possible results. In this case, the only difference between the bounded linear operators  $L[\mathcal{B}]$  and  $L[\mathcal{H}]$  is that  $L[\mathcal{B}]$  is not a  $C^*$ -algebra. Our main tool is a new representation for the dual space. We embed  $\mathcal{B}$  into a (single) Hilbert space  $\mathcal{H}$  that allows us to define an adjoint  $A^*$  on  $\mathcal{B}$  for each closed densely defined linear operator A. We are also able to define a natural Schatten class structure for  $L[\mathcal{B}]$ , that is almost identical to the Schatten class on  $\mathcal{H}$ .
- 1.2. **Preliminaries.** The following theorem is due to Lax [L].

**Theorem 1.2.** (Lax's Theorem) Let  $\mathcal{B}$  be a separable Banach space that is continuously and densely embedded in a Hilbert space  $\mathcal{H}$  and let T be a bounded linear operator on  $\mathcal{B}$  that is symmetric with respect to the inner product of  $\mathcal{H}$  (i.e.,  $(Tu, v)_{\mathcal{H}} = (u, Tv)_{\mathcal{H}}$  for all  $u, v \in \mathcal{B}$ ). Then:

(1) The operator T is bounded with respect to the  $\mathcal H$  norm and

$$||T^*T||_{\mathcal{H}} = ||T||_{\mathcal{H}}^2 \leqslant k ||T||_{\mathcal{B}}^2,$$

where k is a positive constant.

- (2) The spectrum of T relative to  $\mathcal{H}$  is a subset of the spectrum of T relative to  $\mathcal{B}$ .
- (3) The point spectrum of T relative to  $\mathcal{H}$  is a equal to the point spectrum of T relative to  $\mathcal{B}$ .

**Definition 1.3.** A family of vectors in a Banach space  $\{\mathcal{E}_n\}\subset\mathcal{B}$ , is called a Schauder basis (S-basis) for  $\mathcal{B}$  if  $\|\mathcal{E}_n\|_{\mathcal{B}}=1$  and, for each  $u\in\mathcal{B}$ , there is a unique sequence  $(u_n)$  of scalars such that

$$u = \lim_{k \to \infty} \sum_{n=1}^{k} u_n \mathcal{E}_n = \sum_{n=1}^{\infty} u_n \mathcal{E}_n.$$

For example, if  $\mathcal{B} = L^p[0,1]$ , 1 , the family of vectors

$$\{1, \cos(2\pi t), \sin(2\pi t)\cos(4\pi t), \sin(4\pi t), \dots\}$$

is a norm one S-basis for  $\mathcal{B}$ .

**Definition 1.4.** A duality map  $\mathcal{J}: \mathcal{B} \mapsto \mathcal{B}'$ , is a set

$$\mathcal{J}(u) = \left\{ u^* \in \mathcal{B}' \left| \langle u, u^* \rangle = \|u\|_{\mathcal{B}}^2 = \|u^*\|_{\mathcal{B}'}^2 \right\}, \ \forall u \in \mathcal{B}.$$

If  $\mathcal{B}$  is uniformly convex,  $\mathcal{J}(u)$  contains a unique functional  $u^* \in \mathcal{B}'$  for each  $u \in \mathcal{B}$ .

Let  $\Omega$  be a bounded open subset of  $\mathbb{R}^n$ ,  $n \in \mathbb{N}$ . If  $u \in L^p[\Omega] = \mathcal{B}$ , 1 , then the standard example is

$$u^* = \mathcal{J}(u)(x) = \|u\|_p^{2-p} |u(x)|^{p-2} u(x) \in L^q[\Omega], \ \frac{1}{p} + \frac{1}{q} = 1.$$

Furthermore,

$$\langle u, u^* \rangle = \|u\|_p^{2-p} \int_{\Omega} |u(x)|^p d\lambda_n(x) = \|u\|_p^2 = \|u^*\|_q^2$$
 (1.1)

It can be shown that  $\mathcal{B}$  is uniformly convex and that  $u^* = \mathcal{J}(u)$  is uniquely defined for each  $u \in \mathcal{B}$ . Thus, if  $\{u_n\}$  is an S-basis for  $L^p[\Omega]$ , then, when normalized, the family vectors  $\{u_n^*\}$  is an S-basis for  $L^q[\Omega] = (L^p[\Omega])'$ . The relationship between u and  $u^*$  is nonlinear. In the next section we prove the remarkable result, that there is another representation of  $\mathcal{B}'$ , with  $u^* = \mathbf{J}_{\mathcal{B}}(u)$  linear, for each  $u \in \mathcal{B}$ . (However,  $u^*$  is no longer a duality mapping.)

# 2. The Natural Hilbert space for a Uniformly Convex Banach Space

In this section we construct the natural Hilbert space for a uniformly convex Banach space with an S-basis. (For this, we only need the Kuelbs part of Theorem 1.1.) Fix  $\mathcal{B}$  and let  $\{\mathcal{E}_n\}$  be a S-basis for  $\mathcal{B}$ . For each  $\mathcal{E}_n$ , let  $\mathcal{E}_n^*$  be the corresponding dual vector in  $\mathcal{B}'$  and set  $t_n = 2^{-n}$ . For each pair  $u, v \in \mathcal{B}$  define a inner product:

$$(u,v) = \sum_{n=1}^{\infty} t_n \langle \mathcal{E}_n^*, u \rangle \overline{\langle \mathcal{E}_n^*, v \rangle}$$

and let  $\mathcal{H}$  be the completion of  $\mathcal{B}$  in the induced norm. Thus,  $\mathcal{B} \subset \mathcal{H}$  densely and

<span id="page-2-0"></span>
$$\|u\|_{\mathcal{H}} = \left[\sum_{n=1}^{\infty} t_n |\langle \mathcal{E}_n^*, u \rangle|^2\right]^{1/2} \le \sup_{n} |\langle \mathcal{E}_n^*, u \rangle| \le \sup_{\|\mathcal{E}^*\|_{\mathcal{B}'} \le 1} |\langle \mathcal{E}^*, u \rangle| = \|u\|_{\mathcal{B}}, \quad (2.1)$$

so that the embedding is both dense and continuous. (It is clear that  $\mathcal{H}$  is unique up to a change of S-basis.)

2.1. The Hilbert Space Representation. In this section, we show that the dual space of a uniformly convex Banach space has a "Hilbert space representation".

**Definition 2.1.** If  $\mathcal{B}$  be a Banach space, we say that  $\mathcal{B}'$  has a Hilbert space representation if there exists a Hilbert space  $\mathcal{H}$ , with  $\mathcal{B} \subset \mathcal{H}$  as a continuous dense embedding and for each  $u' \in \mathcal{B}'$ ,  $u' = (\cdot, u)_{\mathcal{H}}$  for some  $u \in \mathcal{B}$ .

**Theorem 2.2.** If  $\mathcal{B}$  be a uniformly convex Banach space with an S-basis, then  $\mathcal{B}'$  has a Hilbert space representation.

*Proof.* Let  $\mathcal{H}$  be the natural Hilbert space for  $\mathcal{B}$  and let  $\mathbf{J}$  be the natural linear mapping from  $\mathcal{H} \to \mathcal{H}'$ , defined by

$$\langle v, \mathbf{J}(u) \rangle = (v, u)_{\mathcal{H}}, \text{ for all } u, v \in \mathcal{H}.$$

It is easy to see that  $\mathbf{J}$  is bijective and  $\mathbf{J}^* = \mathbf{J}$ . First, we note that the restriction of  $\mathbf{J}$  to  $\mathcal{B}$ ,  $\mathbf{J}_{\mathcal{B}}$ , maps  $\mathcal{B}$  to a unique subset of linear functionals  $\{\mathbf{J}_{\mathcal{B}}(u), u \in \mathcal{B}\}$  and,  $\mathbf{J}_{\mathcal{B}}(u+v) = \mathbf{J}_{\mathcal{B}}(u) + \mathbf{J}_{\mathcal{B}}(v)$ , for each  $u, v \in \mathcal{B}$ . We are done if we can prove that  $\{\mathbf{J}_{\mathcal{B}}(u), u \in \mathcal{B}\} = \mathcal{B}'$ . For this, it suffices to show that  $\mathbf{J}_{\mathcal{B}}(u)$  is bounded for each  $u \in \mathcal{B}$ . Since  $\mathcal{B}$  is dense in  $\mathcal{H}$ , from equation (2.1) we have:

$$\|\mathbf{J}_{\mathcal{B}}(u)\|_{\mathcal{B}'} = \sup_{v \in \mathcal{B}} \frac{\langle v, \mathbf{J}_{\mathcal{B}}(u) \rangle}{\|v\|_{\mathcal{B}}} \leqslant \sup_{v \in \mathcal{B}} \frac{\langle v, \mathbf{J}_{\mathcal{B}}(u) \rangle}{\|v\|_{\mathcal{H}}} = \|u\|_{\mathcal{H}} \leqslant \|u\|_{\mathcal{B}}.$$

Thus,  $\{\mathbf{J}_{\mathcal{B}}(u), u \in \mathcal{B}\} \subset \mathcal{B}'$ . Since  $\mathcal{B}$  is uniformly convex,  $\{\mathbf{J}_{\mathcal{B}}(u), u \in \mathcal{B}\} = \mathcal{B}'$ .

2.2. Construction of the adjoint on  $\mathcal{B}$ . We can now show that each closed densely linear operator on  $\mathcal{B}$  has a natural adjoint defined on  $\mathcal{B}$ .

<span id="page-3-2"></span>**Theorem 2.3.** Let  $\mathcal{B}$  be a uniformly convex Banach space with an S-basis. If  $\mathbb{C}[\mathcal{B}]$  denotes the closed densely linear operators on  $\mathcal{B}$  and  $L[\mathcal{B}]$  denotes the bounded linear operators, then every  $A \in \mathcal{C}[\mathcal{B}]$  has a well defined adjoint  $A^* \in \mathcal{C}[\mathcal{B}]$ . Furthermore, if  $A \in L[\mathcal{B}]$ , then  $A^* \in L[\mathcal{B}]$  with:

- $(1) \ (aA)^* = \bar{a}A^*,$
- (2)  $A^{**} = A$ .
- (3)  $(A^* + B^*) = A^* + B^*$
- $(4) (AB)^* = B^*A^*$  and
- $(5) \|A^*A\|_{\mathcal{B}} \le \|A\|_{\mathcal{B}}^2.$

Thus,  $L[\mathcal{B}]$  is a \*algebra.

*Proof.* First, let **J** be the natural linear mapping from  $\mathcal{H} \to \mathcal{H}'$  and let  $\mathbf{J}_{\mathcal{B}}$  be the restriction of **J** to  $\mathcal{B}$ . If  $A \in \mathcal{C}[\mathcal{B}]$ , then  $A'\mathbf{J}_{\mathcal{B}} : \mathcal{B}' \to \mathcal{B}'$ . Since A' is closed and densely defined, it follows that  $\mathbf{J}_{\mathcal{B}}^{-1}A'\mathbf{J}_{\mathcal{B}} : \mathcal{B} \to \mathcal{B}$  is a closed and densely defined linear operator. We define  $A^* = [\mathbf{J}_{\mathcal{B}}^{-1}A'\mathbf{J}_{\mathcal{B}}] \in \mathcal{C}[\mathcal{B}]$ . If  $A \in L[\mathcal{B}]$ ,  $A^* = \mathbf{J}_{\mathcal{B}}^{-1}A'\mathbf{J}_{\mathcal{B}}$  is defined on all of  $\mathcal{B}$ . By the Closed Graph Theorem,  $A^* \in L[\mathcal{B}]$ . The proofs of (1)-(3) are straight forward. To prove (4),

<span id="page-3-0"></span>
$$(BA)^* = \mathbf{J}_{\mathcal{B}}^{-1}(BA)'\mathbf{J}_{\mathcal{B}} = \mathbf{J}_{\mathcal{B}}^{-1}A'B'\mathbf{J}_{\mathcal{B}}$$
$$= \left[\mathbf{J}_{\mathcal{B}}^{-1}A'\mathbf{J}_{\mathcal{B}}\right]\left[\mathbf{J}_{\mathcal{B}}^{-1}B'\mathbf{J}_{\mathcal{B}}\right] = A^*B^*.$$
 (2.2)

If we replace B by  $A^*$  in equation (2.2), noting that  $A^{**} = A$ , we also see that  $(A^*A)^* = A^*A$ . To prove (5), we first see that:

$$\langle A^*Av, \mathbf{J}_{\mathcal{B}}(u) \rangle = (A^*Av, u)_{\mathcal{H}} = (v, A^*Au)_{\mathcal{H}}.$$

so that  $A^*A$  is symmetric. Thus, by Lax's Theorem,  $A^*A$  has a bounded extension to  $\mathcal{H}$  and  $\|A^*A\|_{\mathcal{H}} \leq k \|A^*A\|_{\mathcal{B}}$ , where k is a positive constant. We also have that

<span id="page-3-1"></span>
$$||A^*A||_{\mathcal{B}} \leqslant ||A^*||_{\mathcal{B}} ||A||_{\mathcal{B}} \leqslant ||A||_{\mathcal{B}}^2. \tag{2.3}$$

It follows that  $||A^*A||_{\mathcal{B}} \leq ||A||_{\mathcal{B}}^2$ . If equality holds in (2.3), for all  $A \in L[\mathcal{B}]$ , then it is a  $C^*$ -algebra. It is well-known that this is true if and only if  $\mathcal{B}$  is a Hilbert space. Thus, in general the inequality in (2.3) is strict.

2.3. Operators on  $\mathcal{B}$ . For the remainder of the paper, we assume that  $\mathcal{B}$  is uniformly convex and  $\mathcal{B}'$  carries its Hilbert space representation.

**Definition 2.4.** Let U be bounded,  $A \in \mathbb{C}[\mathcal{B}]$  and let  $\mathcal{U}, \mathcal{V}$  be subspaces of  $\mathcal{B}$ . Then:

- (1) A is said to be naturally self-adjoint if  $D(A) = D(A^*)$  and  $A = A^*$ .
- (2) A is said to be normal if  $D(A) = D(A^*)$  and  $AA^* = A^*A$ .
- (3) U is unitary if  $UU^* = U^*U = I$ .
- (4) The subspace  $\mathcal{U}$  is  $\perp$  to  $\mathcal{V}$  if and only, for each  $v \in \mathcal{V}$  and  $\forall u \in \mathcal{U}$ ,  $(v, u)_{\mathcal{H}} = 0$  and, for each  $u \in \mathcal{U}$  and  $\forall v \in \mathcal{V}$ ,  $(u, v)_{\mathcal{H}} = 0$ .

The last definition is transparent since, orthogonal subspaces in  $\mathcal{H}$  induce orthogonal subspaces in  $\mathcal{B}$ .

**Theorem 2.5.** (Gram-Schmidt) For each fixed basis  $\{\varphi_i, 1 \leq i < \infty\}$  of  $\mathcal{B}$ , there is at least one set of dual functionals  $\{S_i\}$  such that  $\{\{\psi_i\}, \{S_i\}, 1 \leq i < \infty\}$  is a biorthonomal set of vectors for  $\mathcal{B}$ , (i.e.,  $\langle \psi_i, S_i \rangle = \delta_{ij}$ ).

*Proof.* Since each  $\varphi_i$  is in  $\mathcal{H}$ , we can construct an orthogonal set of vectors  $\{\phi_i, 1 \leq i < \infty\}$  in  $\mathcal{H}$  by the standard Gram-Schmidt process. Set  $\psi_i = \phi_i/\|\phi_i\|_{\mathcal{B}}$  and let  $\psi_i^* = \mathbf{J}(\psi)$ . From here, it is easy to check that  $\{\{\psi_i\}, \{\psi_i^*\}, 1 \leq i < \infty\}$  is a biorthonormal set and the family  $\{\psi_i^*\}$  is unique.

We close this section with the following observation about the use of  $\mathcal{H}$ . Let A be any closed densely defined naturally selfadjoint linear operator on  $\mathcal{B}$  with a discrete spectrum  $\{\lambda_i\}$ . It can be extended to  $\mathcal{H}$  with the same properties. If we compute the ratio  $\frac{\langle A\psi,\psi^*\rangle}{\langle \psi,\psi^*\rangle}$  in  $\mathcal{B}$ , it will be "close" to the value of  $\frac{(\bar{A}\psi,\psi)_{\mathcal{H}}}{(\psi,\psi)_{\mathcal{H}}}$  in  $\mathcal{H}$ . By Lax's Theorem [L], the extension from  $\mathcal{B}$  to  $\mathcal{H}$  does not change the point spectrum, so we can use the min-max theorem on  $\mathcal{H}$  to compute the eigenvalues and eigenfunctions of A via  $\bar{A}$  exactly. Since  $\mathcal{B}$  is dense in  $\mathcal{H}$ , it follows that the min-max theorem also holds on  $\mathcal{B}$ .

2.3.1. Selfadjointness. With respect to our definition of natural selfadjointness, the following related definition is due to Palmer [P], where the operator is called symmetric. This is essentially the same as a Hermitian operator as defined by Lumer [LU]. (An operator A is dissipative if -A is accretive.)

**Definition 2.6.** A closed densely defined linear operator A on  $\mathcal{B}$  is called self-conjugate if both iA and -iA are dissipative.

**Theorem 2.7.** (Vidav-Palmer) A linear operator A, defined on  $\mathcal{B}$ , is self-conjugate if and only if iA and -iA are generators of isometric semigroups.

**Theorem 2.8.** The operator A, defined on  $\mathcal{B}$ , is self-conjugate if and only if it is naturally self-adjoint.

*Proof.* Let  $\bar{A}$  and  $\bar{A}^*$  be the closed densely defined extensions of A and  $A^*$  to  $\mathcal{H}$ . On  $\mathcal{H}$ ,  $\bar{A}$  is naturally self-adjoint if and only if  $i\bar{A}$  generates a unitary group, if and only if it is self-conjugate. Thus, both definitions coincide on  $\mathcal{H}$ . It follows that the restrictions coincide on  $\mathcal{B}$ .

The proof of the last theorem represents a general approach for proving new results for  $\mathcal{B}$ . The following are two representative.

**Theorem 2.9.** (Polar Representation) Let  $\mathcal{B}$  be a uniformly convex Banach space with an S-basis. If  $A \in \mathbb{C}[\mathcal{B}]$ , then there exists a partial isometry U and a naturally self-adjoint operator T, with D(T) = D(A) and A = UT. Furthermore,  $T = [A^*A]^{1/2}$ , in a well-defined sense.

**Theorem 2.10.** (Spectral Representation) Let  $\mathcal{B}$  be a uniformly convex Banach space with an S-basis and let  $A \in \mathbb{C}[\mathcal{B}]$ , be a naturally self-adjoint linear operator. Then, there exists a operator-valued spectral measure  $E_x, x \in \mathbb{R}$  and for each  $u \in D(A)$ ,

$$Au = \int_{\mathbb{R}} x dE_x(u).$$

## 3. Examples

The following related Hilbert space is more general. It is a concrete implementation of the abstract construction of Kuelbs [K] and, in a different form, is due to Steadman [ST]. It is constructed over  $L^1$ , which is not uniformly convex, but is more suitable for applications. It was first used to provide a rigorous foundation for the Feynman path integral formulation of quantum mechanics in [GZ]. We use it in this section to provide a few concrete examples.

3.1. The space  $KS^2[\mathbb{R}^n]$ . On  $\mathbb{R}^n$  let  $\mathbb{Q}^n$  be the set  $\{\mathbf{x} = (x_1, x_2 \cdots, x_n) \in \mathbb{R}^n\}$  such that  $x_i$  is rational for each i. Since this is a countable dense set in  $\mathbb{R}^n$ , we can arrange it as  $\mathbb{Q}^n = \{\mathbf{x}^1, \mathbf{x}^2, \mathbf{x}^3, \cdots\}$ . For each l and i, let  $\mathbf{B}_l(\mathbf{x}^i)$  be the closed cube centered at  $\mathbf{x}^i$ , with sides parallel to the coordinate axes and diagonal  $r_l = 2^{-l}, l \in \mathbb{N}$ . Now choose the natural order which maps  $\mathbb{N} \times \mathbb{N}$  bijectively to  $\mathbb{N}$ :

$$\{(1,1), (2,1), (1,2), (1,3), (2,2), (3,1), (3,2), (2,3), \dots\}.$$

Let  $\{\mathbf{B}_k, k \in \mathbb{N}\}$  be the resulting set of (all) closed cubes  $\{\mathbf{B}_l(\mathbf{x}^i) \mid (l,i) \in \mathbb{N} \times \mathbb{N}\}$  centered at a point in  $\mathbb{Q}^n$  and let  $\mathcal{E}_k(\mathbf{x})$  be the characteristic function of  $\mathbf{B}_k$ , so that  $\mathcal{E}_k(\mathbf{x})$  is in  $L^p[\mathbb{R}^n]$  for  $1 \leq p \leq \infty$ . Define  $F_k(\cdot)$  on  $L^1[\mathbb{R}^n]$  by

$$F_k(f) = \int_{\mathbb{R}^n} \mathcal{E}_k(\mathbf{x}) f(\mathbf{x}) d\lambda_n(\mathbf{x}). \tag{3.1}$$

It is clear that  $F_k(\cdot)$  is a bounded linear functional on  $L^p[\mathbb{R}^n]$  for each k,  $||F_k|| \le 1$  and, if  $F_k(f) = 0$  for all k, f = 0 so that  $\{F_k\}$  is fundamental on  $L^p[\mathbb{R}^n]$  for  $1 \le p \le \infty$ . Set  $t_k = 2^{-k}$ , so that  $\sum_{k=1}^{\infty} t_k = 1$  and define a measure  $d\mu$  on  $\mathbb{R}^n \times \mathbb{R}^n$  by:

$$d\mu = \left[\sum_{k=1}^{\infty} t_k \mathcal{E}_k(\mathbf{x}) \mathcal{E}_k(\mathbf{y})\right] d\lambda_n(\mathbf{x}) d\lambda_n(\mathbf{y}).$$

To construct our Hilbert space, define an inner product  $(\cdot)$  on  $L^1[\mathbb{R}^n]$  by

$$(f,g) = \int_{\mathbb{R}^n \times \mathbb{R}^n} f(\mathbf{x}) g(\mathbf{y})^* d\mu$$

$$= \sum_{k=1}^{\infty} t_k \left[ \int_{\mathbb{R}^n} \mathcal{E}_k(\mathbf{x}) f(\mathbf{x}) d\lambda_n(\mathbf{x}) \right] \left[ \int_{\mathbb{R}^n} \mathcal{E}_k(\mathbf{y}) g(\mathbf{y}) d\lambda_n(\mathbf{y}) \right]^*.$$
(3.2)

We call the completion of  $L^1[\mathbb{R}^n]$ , with the above inner product, the Kuelbs-Steadman space,  $KS^2[\mathbb{R}^n]$ . Steadman [ST] constructed this space by modifying a method developed by Kuelbs [K] for other purposes. Her interest was in showing that  $L^1[\mathbb{R}]$  can be densely and continuously embedded in a Hilbert space which contains the non-absolutely integrable functions. To see that this is the case, suppose f is a non-absolutely integrable function, say Henstock-Kurzweil integral (or of any other type see [H]), then:

$$||f||_{KS^2}^2 = \sum_{k=1}^{\infty} t_k \left| \int_{\mathbb{R}^n} \mathcal{E}_k(\mathbf{x}) f(\mathbf{x}) d\lambda_n(\mathbf{x}) \right|^2$$

$$\leq \sup_k \left| \int_{\mathbb{R}^n} \mathcal{E}_k(\mathbf{x}) f(\mathbf{x}) d\lambda_n(\mathbf{x}) \right|^2 < \infty.$$

Since the absolute value is outside the integral, we see that  $f \in KS^2[\mathbb{R}^n]$  for any of the definitions of a non-absolute integral (see [GO]). A detailed discussion of this space and its relationship to the Feynman path integral formulation to quantum mechanics, can be found in [GZ]

**Theorem 3.1.** The space  $KS^2[\mathbb{R}^n]$  contains  $L^p[\mathbb{R}^n]$  (for each  $p, 1 \leq p \leq \infty$ ) as dense subspaces.

*Proof.* By construction, we know that  $KS^2[\mathbb{R}^n]$  contains  $L^1[\mathbb{R}^n]$  densely. Thus, we need only show that  $L^q[\mathbb{R}^n] \subset KS^2[\mathbb{R}^n]$  for  $q \neq 1$ . Let  $f \in L^q[\mathbb{R}^n]$  and  $q < \infty$ . Since  $|\mathcal{E}(\mathbf{x})| = \mathcal{E}(\mathbf{x}) \leqslant 1$  and  $|\mathcal{E}(\mathbf{x})|^q \leqslant \mathcal{E}(\mathbf{x})$ , we have

$$||f||_{KS^{2}} = \left[ \sum_{k=1}^{\infty} t_{k} \left| \int_{\mathbb{R}^{n}} \mathcal{E}_{k}(\mathbf{x}) f(\mathbf{x}) d\lambda_{n}(\mathbf{x}) \right|^{\frac{2q}{q}} \right]^{1/2}$$

$$\leq \left[ \sum_{k=1}^{\infty} t_{k} \left( \int_{\mathbb{R}^{n}} \mathcal{E}_{k}(\mathbf{x}) \left| f(\mathbf{x}) \right|^{q} d\lambda_{n}(\mathbf{x}) \right)^{\frac{2}{q}} \right]^{1/2}$$

$$\leq \sup_{k} \left( \int_{\mathbb{R}^{n}} \mathcal{E}_{k}(\mathbf{x}) \left| f(\mathbf{x}) \right|^{q} d\lambda_{n}(\mathbf{x}) \right)^{\frac{1}{q}} \leq ||f||_{q}.$$

Hence,  $f \in KS^2[\mathbb{R}^n]$ . For  $q = \infty$ , first note that  $vol(\mathbf{B}_k)^2 \leq \left\lceil \frac{1}{2\sqrt{n}} \right\rceil^{2n}$ , so we have

$$||f||_{KS^2} = \left[ \sum_{k=1}^{\infty} t_k \left| \int_{\mathbb{R}^n} \mathcal{E}_k(\mathbf{x}) f(\mathbf{x}) d\lambda_n(\mathbf{x}) \right|^2 \right]^{1/2}$$

$$\leq \left[ \left[ \sum_{k=1}^{\infty} t_k [vol(\mathbf{B}_k)]^2 \right] [ess \sup |f|]^2 \right]^{1/2} \leq \left[ \frac{1}{2\sqrt{n}} \right]^n ||f||_{\infty}.$$

Thus  $f \in KS^2[\mathbb{R}^n]$ , and  $L^{\infty}[\mathbb{R}^n] \subset KS^2[\mathbb{R}^n]$ .

The fact that  $L^{\infty}[\mathbb{R}^n] \subset KS^2[\mathbb{R}^n]$ , while  $KS^2[\mathbb{R}^n]$  is separable makes it clear in a very forceful manner that separability is not an inherited property. We note that, since  $L^1[\mathbb{R}^n] \subset KS^2[\mathbb{R}^n]$  and  $KS^2[\mathbb{R}^n]$  is reflexive, the second dual  $L^1[\mathbb{R}^n]'' = \mathfrak{M}[\mathbb{R}^n] \subset KS^2[\mathbb{R}^n]$ . Recall that  $\mathfrak{M}[\mathbb{R}^n]$  is the space of bounded finitely additive set functions defined on the Borel sets  $\mathfrak{B}[\mathbb{R}^n]$ . This space contains the Dirac delta measure and free-particle Green's function for the Feynman path integral.

The next result is an unexpected benefit.

**Theorem 3.2.** Let  $f_n \to f$  weakly in  $L^p$ ,  $1 \le p \le \infty$ , then  $f_n \to f$  strongly in  $KS^2$  (i.e., every weakly compact subset of  $L^p$  is compact in  $KS^2$ ).

*Proof.* The proof of follows from the fact that, if  $\{f_n\}$  is any weakly convergent sequence in  $L^p$  with limit f, then

$$\int_{\mathbb{R}^n} \mathcal{E}_k(\mathbf{x}) \left[ f_n(\mathbf{x}) - f(\mathbf{x}) \right] d\lambda_n(\mathbf{x}) \to 0$$

for each k. It follows that  $\{f_n\}$  converges strongly to f in  $KS^2$ .

Let A be a closed densely defined linear operator defined on  $L^p[\mathbb{R}^n]$ ,  $1 , and let A' be the dual defined on <math>L^q[\mathbb{R}^n]$ ,  $\frac{1}{p} + \frac{1}{q} = 1$ . It is easy to show that, if A' is densely defined on  $L^p[\mathbb{R}^n]$ , it has a closed extension to  $L^p[\mathbb{R}^n]$ .

**Example 3.3.** Let A be a second order differential operator on  $L^p[\mathbb{R}^n]$ , of the form

$$A = \sum_{i,j=1}^{n} a_{ij}(\mathbf{x}) \frac{\partial^2}{\partial x_i \partial x_j} + \sum_{i,j=1}^{n} x_i b_{ij}(\mathbf{x}) \frac{\partial}{\partial x_j},$$

where  $\mathbf{a}(\mathbf{x}) = [a_{ij}(\mathbf{x})]$  and  $\mathbf{b}(\mathbf{x}) = [b_{ij}(\mathbf{x})]$  are matrix-valued functions in  $\mathbb{C}_c^{\infty}[\mathbb{R}^n \times \mathbb{R}^n]$  (infinitely differentiable functions with compact support). We also assume that, for all  $\mathbf{x} \in \mathbb{R}^n$  det  $[a_{ij}(\mathbf{x})] > \varepsilon$  and the imaginary part of the eigenvalues of  $\mathbf{b}(\mathbf{x})$  are bounded above by  $-\varepsilon$ , for some  $\varepsilon > 0$ . Note, since we don't require  $\mathbf{a}$  or  $\mathbf{b}$  to be symmetric,  $A \neq A'$ .

It is well-known that  $\mathbb{C}_c^{\infty}[\mathbb{R}^n] \subset L^p[\mathbb{R}^n] \cap L^q[\mathbb{R}^n]$  is dense for all  $1 \leq p \leq q < \infty$ . Furthermore, since A' is invariant on  $\mathbb{C}_c^{\infty}[\mathbb{R}^n]$ ,

$$A': \mathbb{C}_c^{\infty}[\mathbb{R}^n] \subset L^p[\mathbb{R}^n] \to \mathbb{C}_c^{\infty}[\mathbb{R}^n] \subset L^p[\mathbb{R}^n].$$

It follows that A' has a closed extension to  $L^q[\mathbb{R}^n]$ . (In this case, we do not need  $\mathcal{H}$  directly, we can identify  $\mathbf{J}$  with the identity on  $\mathcal{H}$  and  $A^*$  with A'.)

**Remark 3.4.** For a general A, which is closed and densely defined on  $L^p[\mathbb{R}^n]$ , we know that it is densely defined on  $KS^2[\mathbb{R}^n]$ . Thus, it has a well-defined adjoint  $A^*$  on  $KS^2[\mathbb{R}^n]$ . By Theorem 2.3, we can take the restriction of  $A^*$  from  $KS^2[\mathbb{R}^n]$  to obtain our adjoint on  $L^q[\mathbb{R}^n]$ .

3.1.1. Example: Integral Operators. In one dimension, the Hilbert transform can be defined on  $L^2[\mathbb{R}]$  via its Fourier transform:

$$\widehat{H(f)} = -i\operatorname{sgn} x\,\widehat{f}.$$

It can also be defined directly as principal-value integral:

$$(Hf)(x) = \lim_{\varepsilon \to 0} \frac{1}{\pi} \int_{|x-y| \ge \varepsilon} \frac{f(y)}{x-y} dy.$$

For a proof of the following results see Grafakos [GRA], chapter 4.

**Theorem 3.5.** The Hilbert transform on  $L^2[\mathbb{R}]$  satisfies:

- (1) H is an isometry,  $||H(f)||_2 = ||f||_2$  and  $H^* = -H$ . (2) For  $f \in L^p[\mathbb{R}]$ ,  $1 , there exists a constant <math>C_p > 0$  such that,

<span id="page-8-0"></span>
$$||H(f)||_{p} \le C_{p}||f||_{p}. \tag{3.3}$$

The next result is technically obvious, but conceptually non-trivial.

Corollary 3.6. The adjoint of H, H\* defines a bounded linear operator on  $L^p[\mathbb{R}]$ for  $1 , and <math>H^*$  satisfies equation (3.3) for the same constant  $C_p$ .

The Riesz transform,  $\mathbf{R}$ , is the *n*-dimensional analogue of the Hilbert transform and its  $j^{\text{th}}$  component is defined for  $f \in L^p[\mathbb{R}^n], \ 1 , by:$ 

$$R_j(f) = c_n \lim_{\varepsilon \to 0} \int_{|\mathbf{y} - \mathbf{x}| \ge \varepsilon} \frac{y_j - x_j}{|\mathbf{y} - \mathbf{x}|^{n+1}} f(\mathbf{y}) d\mathbf{y}, \quad c_n = \frac{\Gamma\left(\frac{N+1}{2}\right)}{\pi^{(n+1)/2}}.$$

**Definition 3.7.** Let  $\Omega$  be defined on the unit sphere  $S^{n-1}$  in  $\mathbb{R}^n$ .

- (1) The function  $\Omega(x)$  is said to be homogeneous of degree n if  $\Omega(tx) = t^n \Omega(x)$ .
- (2) The function  $\Omega(x)$  is said to have the cancellation property if

$$\int_{S^{n-1}} \Omega(\mathbf{y}) d\sigma(\mathbf{y}) = 0, \text{ where } d\sigma \text{ is the induced Lebesgue measure on } S^{n-1}.$$

(3) The function  $\Omega(x)$  is said to have the Dini-type condition if

$$\sup_{\substack{|\mathbf{x}-\mathbf{y}|\leqslant \delta\\ |\mathbf{x}|=|\mathbf{y}|=1}} |\Omega(\mathbf{x}) - \Omega(\mathbf{y})| \leqslant \omega(\delta) \Rightarrow \int_0^1 \frac{\omega(\delta)d\delta}{\delta} < \infty.$$

A proof of the following theorem can be found in Stein [STE] (see pg., 39).

**Theorem 3.8.** Suppose that  $\Omega$  is homogeneous of degree 0, satisfying both the cancellation property and the Dini-type condition. If  $f \in L^p[\mathbb{R}^n]$ , 1 and

$$T_{\varepsilon}(f)(\mathbf{x}) = \int_{|\mathbf{y} - \mathbf{x}| \ge \varepsilon} \frac{\Omega(\mathbf{y} - \mathbf{x})}{|\mathbf{y} - \mathbf{x}|^n} f(\mathbf{y}) d\mathbf{y}.$$

Then

(1) There exists a constant  $A_p$ , independent of both f and  $\varepsilon$  such that

$$||T_{\varepsilon}(f)||_{p} \leqslant A_{p}||f||_{p}.$$

(2) Furthermore,  $\lim_{\varepsilon \to 0} T_{\varepsilon}(f) = T(f)$  exists in the  $L^p$  norm and

<span id="page-9-0"></span>
$$||T(f)||_p \leqslant A_p ||f||_p.$$
 (3.4)

Treating  $T_{\varepsilon}(f)$  as a special case of the Henstock-Kurzweil integral, conditions (1) and (2) are automatically satisfied and we can write the integral as

$$T(f)(\mathbf{x}) = \int_{\mathbb{R}^n} \frac{\Omega(\mathbf{y} - \mathbf{x})}{|\mathbf{y} - \mathbf{x}|^n} f(\mathbf{y}) d\mathbf{y}.$$

For  $g \in L^q$ ,  $\frac{1}{p} + \frac{1}{q} = 1$ , we have  $\langle T(f), g \rangle = \langle f, T^*(g) \rangle$ . Using Fubini's Theorem for the Henstock-Kurzweil integral (see [H]), we have that

Corollary 3.9. The adjoint of T,  $T^* = -T$ , is defined on  $L^p$  and satisfies equation (3.4)

It is easy to see that the Riesz transform is a special case of the above Theorem and Corollary.

Another closely related integral operator is the Riesz potential,  $I_{\alpha}(f)(\mathbf{x}) = (-\Delta)^{-\alpha/2} f(\mathbf{x})$ ,  $0 < \alpha < n$ , is defined on  $L^p[\mathbb{R}^n]$ , 1 , by (see Stein [STE], pg., 117):

$$I_{\alpha}(f)(\mathbf{x}) = \gamma^{-1}(\alpha) \int_{\mathbb{R}^n} \frac{f(\mathbf{y}) d\mathbf{y}}{|\mathbf{x} - \mathbf{y}|^{n-\alpha}}, \text{ and } \gamma(\alpha) = 2^{\alpha} \pi^{\frac{n}{2}} \frac{\Gamma(\frac{\alpha}{2})}{\Gamma(\frac{n-\alpha}{2})}.$$

Since the kernel is symmetric, application of Fubini's Theorem shows that the adjoint  $I_{\alpha}^* = I_{\alpha}$ , is also defined on  $L^p[\mathbb{R}^n]$ . Since  $(-\Delta)^{-1}$  is not bounded, we cannot obtain  $L^p$  bounds for  $I_{\alpha}(f)(\mathbf{x})$ . However, if  $1/q = 1/p - \alpha/n$ , we have the following (see Stein [STE], pg., 119)

**Theorem 3.10.** If  $f \in L^p[\mathbb{R}^n]$  and  $0 < \alpha < n$ ,  $1 , <math>1/q = 1/p - \alpha/n$ , then the integral defining  $I_{\alpha}(f)$  converges absolutely for almost all  $\mathbf{x}$ . Furthermore, there is a constant  $A_{p,q}$ , such that

$$||I_{\alpha}(f)||_{q} \leqslant A_{p,q}||f||_{p}.$$
 (3.5)

### 4. SCHATTEN CLASSES ON BANACH SPACES

In this section, we give a natural definition of the Schatten class of operators on  $\mathcal{B}$  (see [SC]) and show that the structure of  $L[\mathcal{B}]$  is almost identical to that of  $L[\mathcal{H}]$ .

4.1. Background: Compact Operators on Banach Spaces. Let  $\mathbb{K}(\mathcal{B})$  be the class of compact operators on  $\mathcal{B}$  and let  $\mathbb{F}(\mathcal{B})$  be the set of operators of finite rank. Recall that, for separable Banach spaces,  $\mathbb{K}(\mathcal{B})$  is an ideal that need not be the maximal ideal in  $L[\mathcal{B}]$ . If  $\mathbb{M}(\mathcal{B})$  is the set of weakly compact operators and  $\mathbb{N}(\mathcal{B})$  is the set of operators that map weakly convergent sequences into strongly convergent sequences, it is known that both are closed two-sided ideals in the operator norm and, in general,  $\mathbb{F}(\mathcal{B}) \subset \mathbb{K}(\mathcal{B}) \subset \mathbb{M}(\mathcal{B})$  and  $\mathbb{F}(\mathcal{B}) \subset \mathbb{K}(\mathcal{B}) \subset \mathbb{N}(\mathcal{B})$  (see part I of Dunford and Schwartz [DS], pg. 553). For reflexive Banach spaces,  $\mathbb{K}(\mathcal{B}) = \mathbb{N}(\mathcal{B})$  and  $\mathbb{M}(\mathcal{B}) = L[\mathcal{B}]$ . For the space of continuous functions  $\mathbb{C}[\Omega]$  on a compact Hausdorff space  $\Omega$ , Grothendieck [GR] has shown that  $\mathbb{M}(\mathcal{B}) = \mathbb{N}(\mathcal{B})$ .

On the other hand, it is shown in part I of Dunford and Schwartz [DS] that, if  $(\Omega, \Sigma, \mu)$  a positive measure space, then for  $L^1(\Omega, \Sigma, \mu)$  we have  $M(\mathcal{B}) \subset N(\mathcal{B})$ .

We assume that  $\mathcal{B}$  is uniformly convex, with a S-basis. In operator theoretic language, our S-basis assumption is that the set of compact operators on  $\mathcal{B}$  have the approximation property, namely that every compact operator can be approximated by operators of finite rank. In this section we will show that the structure of  $L[\mathcal{B}]$  is almost identical to its associated space  $L[\mathcal{H}]$ . The difference is that  $L[\mathcal{B}]$  is not a  $C^*$ -algebra (i.e.,  $||A^*A||_{\mathcal{B}} = ||A||_{\mathcal{B}}^2$ ,  $A \in L[\mathcal{B}]$ , is not true for all A).

Let A be a compact operator on  $\mathcal{B}$  and let  $\bar{A}$  be its extension to  $\mathcal{H}$ . For each compact operator  $\bar{A}$ , there exists an orthonormal basis  $\{\bar{\varphi}_n \mid n \geq 1\}$ , for  $\mathcal{H}$  such that

$$\bar{A} = \sum_{n=1}^{\infty} \mu_n(\bar{A}) \left(\cdot, \bar{\varphi}_n\right)_2 \bar{U} \bar{\varphi}_n.$$

Where the  $\mu_n$  are the eigenvalues of  $[\bar{A}^*\bar{A}]^{1/2} = |\bar{A}|$ , counted by multiplicity and in decreasing order, and  $\bar{U}$  is the partial isometry associated with the polar decomposition of  $\bar{A} = \bar{U} |\bar{A}|$ . Without loss, we can assume that the set of functions  $\{\bar{\varphi}_n | n \geq 1\}$  is contained in  $\mathcal{B}$  and  $\{\varphi_n | n \geq 1\}$  is the normalized version in  $\mathcal{B}$ . If  $\mathbb{S}_p[\mathcal{H}]$  is the Schatten Class of order p in  $L[\mathcal{H}]$ , it is well-known that, if  $\bar{A} \in \mathbb{S}_p[\mathcal{H}]$ , its norm can be represented as:

$$\|\bar{A}\|_{p}^{\mathcal{H}} = \left\{ Tr \left[ \bar{A}^* \bar{A} \right]^{p/2} \right\}^{1/p} = \left\{ \sum_{n=1}^{\infty} \left( \bar{A}^* \bar{A} \bar{\varphi}_n, \bar{\varphi}_n \right)_{\mathcal{H}}^{p/2} \right\}^{1/p}$$
$$= \left\{ \sum_{n=1}^{\infty} \left| \mu_n \left( \bar{A} \right) \right|^p \right\}^{1/p}.$$

**Definition 4.1.** We define the Schatten Class of order p in  $L[\mathcal{B}]$  by:

$$\mathbb{S}_p[\mathcal{B}] = \mathbb{S}_p[\mathcal{H}] \mid_{\mathcal{B}}.$$

Since  $\bar{A}$  is the extension of  $A \in \mathbb{S}_p[\mathcal{B}]$ , we can define A on  $\mathcal{B}$  by

$$A = \sum_{n=1}^{\infty} \mu_n(A) \langle \cdot, \varphi_n^* \rangle U \varphi_n,$$

where  $\varphi_n^*$  is the unique functional in  $\mathcal{B}'$  associated with  $\varphi_n$  and U is the restriction of  $\overline{U}$  to  $\mathcal{B}$ . The corresponding norm of A on  $\mathbb{S}_p[\mathcal{B}]$  is defined by:

$$||A||_{p}^{\mathcal{B}} = \left\{ \sum_{n=1}^{\infty} \left\langle A^* A \varphi_n, \varphi_n^* \right\rangle^{p/2} \right\}^{1/p}.$$

Theorem 4.2. Let  $A \in \mathbb{S}_p[\mathcal{B}]$ , then  $||A||_p^{\mathcal{B}} = ||\bar{A}||_p^{\mathcal{H}}$ .

*Proof.* It is clear that  $\{\varphi_n \mid n \geq 1\}$  is a set of eigenfunctions for  $A^*A$  on  $\mathcal{B}$ . Furthermore,  $A^*A$  is naturally selfadjoint and compact, so that its spectrum is discrete. By Lax's Theorem, the spectrum of  $A^*A$  is unchanged by its extension to  $\mathcal{H}$ . It follows that  $A^*A\varphi_n = |\mu_n(\bar{A})|^2 \varphi_n$ , so that

$$\langle A^* A \varphi_n, \varphi_n^* \rangle = \left| \mu_n(A) \right|^2 \langle \varphi_n, \varphi_n^* \rangle = \left| \mu_n(\bar{A}) \right|^2,$$

and

$$||A||_{p}^{\mathcal{B}} = \left\{ \sum_{n=1}^{\infty} \left\langle A^* A \varphi_n, \varphi_n^* \right\rangle^{p/2} \right\}^{1/p} = \left\{ \sum_{n=1}^{\infty} \left| \mu_n(\bar{A}) \right|^p \right\}^{1/p} = ||\bar{A}||_{p}^{\mathcal{H}}.$$

It is clear that all of the theory of operator ideals on Hilbert spaces extend to uniformly convex Banach spaces with a S-basis in a straightforward way. We state a few of the more important results to give a sense of the power provided by the existence of adjoints for spaces of this type. The first result extends theorems due to Weyl [W], Horn [HO], Lalesco [LA] and Lidskii [LI]. The proofs are all straight forward, for a given A extend it to  $\mathcal{H}$ , use the Hilbert space result and then restrict back to  $\mathcal{B}$ .

**Theorem 4.3.** Let  $A \in \mathbb{K}(\mathcal{B})$ , the set of compact operators on  $\mathcal{B}$ , and let  $\{\lambda_n\}$  be the eigenvalues of A counted up to algebraic multiplicity. If  $\Phi$  is a mapping on  $[0,\infty]$  which is nonnegative and monotone increasing, then we have:

$$\sum_{n=1}^{\infty} \Phi\left(|\lambda_n(A)|\right) \leqslant \sum_{n=1}^{\infty} \Phi\left(\mu_n(A)\right)$$

and

(2) (Horn) If  $A_1, A_2 \in \mathbb{K}(\mathcal{B})$ 

$$\sum_{n=1}^{\infty} \Phi\left(|\lambda_n(A_1 A_2)|\right) \leqslant \sum_{n=1}^{\infty} \Phi\left(\mu_n(A_1)\mu_n(A_2)\right).$$

In case  $A \in \mathbb{S}_1(\mathcal{B})$ , we have:

(3) (Lalesco)

$$\sum\nolimits_{n=1}^{\infty} |\lambda_n(A)| \leqslant \sum\nolimits_{n=1}^{\infty} \mu_n(A)$$

and

(4) (Lidskii)

$$\sum_{n=1}^{\infty} \lambda_n(A) = Tr(A).$$

Simon [SI2] provides a very nice approach to infinite determinants and trace class operators on separable Hilbert spaces. He gives a comparative historical analysis of Fredholm theory, obtaining a new proof of Lidskii's Theorem as a side benefit and some new insights. A review of his paper shows that much of it can be directly extended to operator theory on uniformly convex Banach spaces.

4.2. **Discussion.** On a Hilbert space  $\mathcal{H}$ , the Schatten classes  $\mathbb{S}_p(\mathcal{H})$  are the only ideals in  $\mathbb{K}(\mathcal{H})$ , and  $\mathbb{S}_1(\mathcal{H})$  is minimal. In a general Banach space, this is far from true. A complete history of the subject can be found in the recent book by Pietsch [PI1] (see also Retherford [R], for a nice review). We limit this discussion to a few major topics on the subject. First, Grothendieck [GR] defined an important class of nuclear operators as follows:

**Definition 4.4.** If  $A \in \mathbb{F}(\mathcal{B})$  (the operators of finite rank), define the ideal  $\mathbf{N}_1(\mathcal{B})$  by:

$$\mathbf{N}_1(\mathcal{B}) = \{ A \in \mathbb{F}(\mathcal{B}) \mid \mathbf{N}_1(A) < \infty \},$$

where

$$\mathbf{N}_{1}(A) = \operatorname{glb}\left\{\sum_{n=1}^{m} \|f_{n}\| \|\phi_{n}\| \left| f_{n} \in \mathcal{B}', \ \phi_{n} \in \mathcal{B}, \ A = \sum_{n=1}^{m} \phi_{n} \left\langle \cdot, f_{n} \right\rangle \right.\right\}$$

and the greatest lower bound is over all possible representations for A.

Grothendieck showed that  $\mathbf{N}_1(\mathcal{B})$  is the completion of the finite rank operators and is a Banach space with norm  $\mathbf{N}_1(\cdot)$ . It is also a two-sided ideal in  $\mathbb{K}(\mathcal{B})$ . It is easy to show that:

Corollary 4.5.  $\mathbb{M}(\mathcal{B})$ ,  $\mathbb{N}(\mathcal{B})$  and  $\mathbb{N}_1(\mathcal{B})$  are two-sided \*ideals.

In order to compensate for the (apparent) lack of an adjoint for Banach spaces, Pietsch [PI2], [PI3] defined a number of classes of operator ideals for a given  $\mathcal{B}$ . Of particular importance for our discussion is the class  $\mathbb{C}_p(\mathcal{B})$ , defined by

$$\mathbb{C}_p(\mathcal{B}) = \left\{ A \in \mathbb{K}(\mathcal{B}) \mid \mathbb{C}_p(A) = \sum_{i=1}^{\infty} [s_i(A)]^p < \infty \right\},\,$$

where the singular numbers  $s_n(A)$  are defined by:

$$s_n(A) = \inf \{ ||A - K||_{\mathcal{B}} \mid \text{rank of } K \leq n \}.$$

Pietsch has shown that,  $\mathbb{C}_1(\mathcal{B}) \subset \mathbf{N}_1(\mathcal{B})$ , while Johnson et al [JKMR] have shown that for each  $A \in \mathbb{C}_1(\mathcal{B})$ ,  $\sum_{n=1}^{\infty} |\lambda_n(A)| < \infty$ . On the other hand, Grothendieck [GR] has provided an example of an operator A in  $\mathbf{N}_1(L^{\infty}[0,1])$  with  $\sum_{n=1}^{\infty} |\lambda_n(A)| = \infty$  (see Simon [SII], pg. 118). Thus, it follows that, in general, the containment is strict. It is known that, if  $\mathbb{C}_1(\mathcal{B}) = \mathbf{N}_1(\mathcal{B})$ , then  $\mathcal{B}$  is isomorphic to a Hilbert space (see Johnson et al). It is clear from the above discussion, that:

Corollary 4.6.  $\mathbb{C}_p(\mathcal{B})$  is a two-sided \*ideal in  $\mathbb{K}(\mathcal{B})$ , and  $\mathbb{S}_1(\mathcal{B}) \subset \mathbf{N}_1(\mathcal{B})$ .

For a given Banach space, it is not clear how the spaces  $\mathbb{C}_p(\mathcal{B})$  of Pietsch relate to our Schatten Classes  $\mathbb{S}_p(\mathcal{B})$  (clearly  $\mathbb{S}_p(\mathcal{B}) \subseteq \mathbb{C}_p(\mathcal{B})$ ). Thus, one question is that of the equality of  $\mathbb{S}_p(\mathcal{B})$  and  $\mathbb{C}_p(\mathcal{B})$ . (We suspect that  $\mathbb{S}_1(\mathcal{B}) = \mathbb{C}_1(\mathcal{B})$ .)

#### 5. Conclusion

The most interesting aspect of this paper is the observation that the dual space of a Banach space can have more then one representation. It is well-known that a given Banach space  $\mathcal{B}$ , can have many equivalent norms that generate the same topology. However, the geometric properties of the space depend on the norm used. We have shown that the properties of the linear operators on  $\mathcal{B}$  depend on the family of linear functionals used to represent the dual space  $\mathcal{B}'$ . This approach offers a interesting tool for a closer study of the structure of bounded linear operators on  $\mathcal{B}$ .

## References

[CA] N. L. Carothers, A Short Course on Banach Space Theory, London Math. Soc. Student Texts 64 Cambridge Univ. Press, Cambridge, UK, (2005).

<span id="page-12-0"></span>[DS] N. Dunford and J. T. Schwartz, *Linear Operators Part I: General Theory*, Wiley Classics edition, Wiley Interscience (1988).

- <span id="page-13-1"></span>[G] L. Gross, Abstract Wiener spaces, Proc. Fifth Berkeley Symposium on Mathematics Statistics and Probability, 1965, pp. 31?42. MR 35:3027
- <span id="page-13-2"></span>[GBZS] T. Gill, S. Basu, W. W. Zachary and V. Steadman, Adjoint for operators in Banach spaces, Proceedings of the American Mathematical Society, 132 (2004), 1429-1434.
- <span id="page-13-9"></span>[GO] R. A. Gordon, The Integrals of Lebesgue, Denjoy, Perron and Henstock, Graduate Studies in Mathematics, Vol. 4, Amer. Math. Soc., (1994).
- <span id="page-13-13"></span>[GR] A. Grothendieck, Products tensoriels topologiques et espaces nucleaires, Memoirs of the American Mathematical Society, 16 (1955).
- <span id="page-13-10"></span>[GRA] L. Grafakos, Classical and Modern Fourier Analysis, Person Prentice-Hall, New Jersey, 2004.
- <span id="page-13-7"></span>[GZ] T.L. Gill and W.W. Zachary, A new class of Banach spaces, J. Phys. A: Math. Theor. 41 (2008) 1-15, doi:10.1088/1751-8113/41/49/495206.
- <span id="page-13-8"></span>[H] R. Henstock, The General Theory of Integration, Clarendon Press, Oxford, (1991).
- <span id="page-13-15"></span>[HO] A. Horn, On the singular values of a product of completely continuous operators, Proc. Nat. Acad. Sci. 36 (1950), 374–375.
- <span id="page-13-23"></span>[JKMR] W. B. Johnson, H. Konig, B. Maurey and J. R. Retherford, Eigenvalues of p-summing and  $l_p$  type operators in Banach space, Journal of Functional Analysis **32** (1978), 353-380.
- <span id="page-13-0"></span>[K] J. Kuelbs, Gaussian measures on a Banach space, Journal of Functional Analysis 5 (1970), 354–367.
- <span id="page-13-3"></span>[L] P. D. Lax, Symmetrizable linear transformations, Comm. Pure Appl. Math. 7 (1954), 633-647
- <span id="page-13-16"></span>[LA] T. Lalesco, Une theoreme sur les noyaux composes, Bull. Acad. Sci. 3 (1914/15), 271–272.
- <span id="page-13-17"></span>[LI] V. B. Lidskii, Non-self adjoint operators with a trace, Dokl. Akad. Nauk. SSSR 125 (1959), 485-487.
- <span id="page-13-5"></span>[LU] G. Lumer, Spectral operators, Hermitian operators and bounded groups, *Acta. Sci. Math.* (Szeged) **25** (1964), 75-85.
- <span id="page-13-4"></span>[P] T. W. Palmer, Unbounded normal operators on Banach spaces, *Trans. Amer. Math. Sci.* 133 (1968), 385-414.
- <span id="page-13-19"></span>[PII] A. Pietsch, History of Banach Spaces and Operator Theory, Birkhäuser, Boston, (2007).
- <span id="page-13-21"></span>[PI2] A. Pietsch, Einige neue Klassen von kompacter linear Abbildungen, Revue der Math. Pures et Appl. (Bucharest), 8 (1963), 423–447.
- <span id="page-13-22"></span>[PI3] A. Pietsch, Eigenvalues and s-Numbers, Cambridge University Press, (1987).
- <span id="page-13-20"></span>[R] J. R. Retherford, Applications of Banach ideals of operators, Bulletin of the American Mathematical Society, 81 (1975), 978-1012.
- <span id="page-13-24"></span>[SI1] B. Simon, *Trace Ideals and their Applications*, London Mathematical Society Lecture Notes Series **35**, Cambridge University Press, New York, (1979).
- <span id="page-13-18"></span>[SI2] B. Simon, Notes on Infinite Determinants of Hilbert Space Operators, Advances in math. 24, (1977) 244-273.
- <span id="page-13-12"></span>[SC] R. Schatten, A Theory of Cross-Spaces, Princeton University Press, Princeton, New Jersey, (1950).
- <span id="page-13-6"></span>[ST] V. Steadman, Theory of operators on Banach spaces, Ph.D thesis, Howard University, (1988).
- <span id="page-13-11"></span>[STE] E. M. Stein, Singular integrals and differentiability properties of functions, Princeton University Press, Princeton, NJ, (1970).
- <span id="page-13-14"></span>[W] H. Weyl, Inequalities between the two kinds of eigenvalues of a linear transformation, Proc. Nat. Acad. Sci. 35 (1949), 408-411.

E-mail address: tgill@access4less.net; tgill@howard.edu

 $<sup>^{1}</sup>$  Department of Mathematics, Howard University, Washington DC 20059, USA.

<sup>2</sup> Department of Mathematics, Howard University, Washington DC 20059, USA.

E-mail address: nubian83@hotmail.com