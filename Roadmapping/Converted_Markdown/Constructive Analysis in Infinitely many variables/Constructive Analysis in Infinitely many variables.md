# CONSTRUCTIVE ANALYSIS IN INFINITELY MANY VARIABLES

TEPPER L. GILL, G. R. PANTSULAIA, AND W. W. ZACHARY\*

<sup>1991</sup> Mathematics Subject Classification. Primary (45) Secondary(46) .

Key words and phrases. infinite-dimensional Lebesgue measure, Gaussian measure,

Fourier transforms, Banach space, Pontryagin duality theory, partial differential operators.

<sup>\*</sup>deceased.

Abstract. In this paper we investigate the foundations for analysis in infinitely-many (independent) variables. We give a topological approach to the construction of the regular σ-finite Kirtadze-Pantsulaia measure on R <sup>∞</sup> (the usual completion of the Yamasaki-Kharazishvili measure), which is an infinite dimensional version of the classical method of constructing Lebesgue measure on R n (see [\[YA1\]](#page-64-0), [\[KH\]](#page-62-0) and [\[KP2\]](#page-62-1)). First we show that von Neumann's theory of infinite tensor product Hilbert spaces already implies that a natural version of Lebesgue measure must exist on R <sup>∞</sup>. Using this insight, we define the canonical version of L 2 [R <sup>∞</sup>, λ∞], which allows us to construct Lebesgue measure on R <sup>∞</sup> and analogues of Lebesgue and Gaussian measure for every separable Banach space with a Schauder basis. When H is a Hilbert space and λ<sup>H</sup> is Lebesgue measure restricted to H, we define sums and products of unbounded operators and the Gaussian density for L 2 [H, λH]. We show that the Fourier transform induces two different versions of the Pontryagin duality theory. An interesting new result is that the character group changes on infinite dimensional spaces when the Fourier transform is treated as an operator. Since our construction provides a complete σ-finite measure space, the abstract version of Fubini's theorem allows us to extend Young's inequality to every separable Banach space with a Schauder basis. We also give constructive examples of partial differential operators in infinitely many variables and briefly discuss the famous partial differential equation derived by Phillip Duncan Thompson [\[PDT\]](#page-63-0), on infinite-dimensional phase space to represent an ensemble of randomly forced two-dimensional viscous flows.

|  | CONSTRUCTIVE ANALYSIS IN INFINITELY MANY VARIABLES |  |  |  |  |  | 3 |
|--|----------------------------------------------------|--|--|--|--|--|---|
|--|----------------------------------------------------|--|--|--|--|--|---|

### Contents

| Introduction                              | 4  |
|-------------------------------------------|----|
| Historical Background                     | 6  |
| Purpose                                   | 10 |
| Summary                                   | 10 |
| 1.<br>Why<br>λ∞<br>Must Exist             | 11 |
| R∞<br>2.<br>Lebesgue Measure on<br>I      | 14 |
| 2.1.<br>The Construction                  | 14 |
| R∞<br>2.2.<br>The Extension to<br>I       | 18 |
| 2.3.<br>Separable Banach Spaces           | 25 |
| 2.4.<br>Translations                      | 28 |
| 2.5.<br>Gaussian measure                  | 29 |
| 2.6.<br>Rotational Invariance             | 32 |
| Discussion                                | 35 |
| Operators<br>3.                           | 35 |
| H2<br>3.1.<br>Bounded Operators on<br>⊗   | 35 |
| H2<br>3.2.<br>Unbounded Operators on<br>⊗ | 38 |
| 4.<br>Function Spaces                     | 43 |
| 1<br>4.1.<br>L<br>-Theory                 | 44 |
| 5.<br>Fourier Transform Theory            | 46 |
| Background                                | 46 |

| 5.1. | Pontryagin Duality Theory I               | 47 |
|------|-------------------------------------------|----|
| 5.2. | 2<br>L<br>-Theory                         | 49 |
| 5.3. | Pontryagin Duality Theory II              | 50 |
| 5.4. | p<br>L<br>-Theory                         | 52 |
| 6.   | Partial Differential Operators (Examples) | 56 |
| 6.1. | Discussion                                | 60 |
| 7.   | Conclusion                                | 61 |
|      | Acknowledgments                           | 61 |
|      | References                                | 61 |

#### <span id="page-3-0"></span>Introduction

On finite-dimensional space it is useful to think of Lebesgue measure in terms of geometric objects (e.g.,volume, surface area, etc.). Thus, it is natural to expect that this measure will leave these objects invariant under translations and rotations, so that rotational and translational invariance is an intrinsic property of Lebesgue measure. However, we then find ourselves disappointed when we try to use this property to help define Lebesgue measure on R∞. A more fundamental problem is that the natural Borel algebra for R∞, B[R∞], does not allow an outer measure (since the measure of any open set is infinite).

The lack of any definitive understanding of the cause for this lack of invariance on R<sup>∞</sup> has led some researchers to believe that it is not possible to have a reasonable version of Lebesgue measure on R<sup>∞</sup> (see, for example, DaPrato [\[DP\]](#page-61-0) or Bakhtin and Mattingly [\[BM\]](#page-61-1)). In many applications, the study of infinite dimensional analysis is restricted to separable Hilbert spaces, using Gaussian measure as a replacement for (the supposed nonexistent) Lebesgue measure. In some cases the Hilbert space structure arises as a natural state space for the modeling of systems. In other cases, both the Hilbert spaces and probability measures are imposed for mathematical convenience and are physically artificial and limiting. However, all reasonable models of infinite dimensional (physical) systems require some functional constraint on the effects of all but a finite number of variables. Thus, what is needed, in general, is the imposition of constraints on the functions while preserving the modeling freedom associated with infinitely-many independent variables (in some well-defined sense). Any attempt to solve this problem necessarily implies a theory of Lebsegue measure on R∞.

Even if a reasonable theory of Lebsegue measure on R<sup>∞</sup> exists, this is not sufficient to make it useful in engineering and science. In addition, all the tools developed for finite-dimensional analysis, differential operators, Fourier transforms, etc are also required. Furthermore, researchers need operational control over the convergence properties of these tools. In particular, one must be able to approximate an infinite-dimensional problem as a natural limit of the finite-dimensional case in a manner that lends itself to computational implementation. This implies that a useful approach also has a well-developed theory of convergence for infinite sums and products of unbounded linear operators.

### <span id="page-5-0"></span>Historical Background

Research into the general problem of Lebesgue measure on infinitedimensional vector spaces and R<sup>∞</sup> in particular, has a long and varied past, with participants living in a number of different countries, during times when scientific communication was constrained by war, isolation and/or national competition. These conditions allowed quite a bit of misinformation and folklore to grow up around the subject, so that even experts may have a limited view of the history. Our own experience suggest that at least a brief survey of some important events is in order. (We do not claim completeness and apologize in advance if we fail to mention equally important contributions.)

Early studies in infinite dimensional analysis focused on the foundations of probability theory and had a broad base of participation. However, the major inputs were made by researchers in Poland, Russia, and France, with later contributions from the US. The first important advance of the general theory was made in 1933 when Haar [\[HA\]](#page-61-2) proved the following theorem:

Theorem 0.1. On every locally compact abelian group G there exists a nonnegative regular measure m (Haar measure) on G, which is not identically zero and is translation invariant. That is, m(A+x) = m(A) for every x ∈ G and every Borel set A in G.

This theorem stimulated interest in the subject and von Neumann [\[VN1\]](#page-64-1) proved that it is the only locally finite left-invariant Borel measure on the group (uniqueness up to a mulitplicative constant). Weil [\[WE\]](#page-64-2) developed an axiomatic approach to the subject, made a number of important refinements and, proved the "Inverse Weil theorem" (in moderm terms):

Theorem 0.2. If G is a (separable) topological group and m is a translation invariant Borel measure on G, then it is always possible to define an equivalent locally compact topology on G.

In 1946, Oxtoby [\[OX\]](#page-62-2) initiated the study of translation-invariant Borel measures on Polish groups (i.e., complete separable metric groups). In this paper, Oxtoby provides a proof of the following result which he attributes to Ulam:

Theorem 0.3. Let G be any complete separable metric group which is not locally compact, and let m be any left-invariant Borel measure in G. Then every neighborhood contains an uncountable number of disjoint mutually congruent sets of equal finite positive measure.

Stated another way, he proved that

Theorem 0.4. There always exists a left-invariant Borel measure on any Polish group which assigns positive finite measure to at least one set and vanishes on singletons. However, a locally finite measure is possible if and only if the group is locally compact.

(In 1967, Vershik [\[V\]](#page-63-1) proved a related result for probability measures.) Apparently uninformed of Oxtoby's work, In 1959 Sudakov [\[SU\]](#page-63-2) independently proved a special case of Theorem 0.4: If R<sup>∞</sup> is regarded as a linear topological space, then there does not exist a σ-finite translation-invariant Borel measure for R∞. In 1964, Elliott and Morse [\[EM\]](#page-61-3) developed a general theory of translation invariant product measures (non-σ-finite) and, in 1965, C. C. Moore [\[MO\]](#page-62-3) initiated the study of measures that are translation invariant with respect to vectors in R<sup>∞</sup> 0 (i.e., the set of sequences that are zero except for a finite number of terms). This work was extended and refined by Hill [\[HI\]](#page-62-4) in 1971.

Motivated by Kakutani's work on infinite product measures [\[KA\]](#page-62-5), a number of young Japanese researchers entered the field. In 1973, Hamachi [\[HA\]](#page-61-2) made major improvements on Hill's work which, indirectly suggested the problem of identifying the largest group T, of admissible translations in the sense of invariance for any σ-finite Borel measure µ on R<sup>∞</sup> which assigns the value of one to [− 1 2 , 1 2 <sup>ℵ</sup><sup>o</sup> and is metrically transitivity with respect to R∞ 0 (equivalently, for each A with µ(A) > 0, there is a sequence (hk) ∈ R<sup>∞</sup> 0 such that µ(R<sup>∞</sup> \ ∪<sup>∞</sup> <sup>k</sup>=1(A + hk)) = 0).

Yamasaki [\[YA1\]](#page-64-0) solved this problem in 1980. Unaware of the Yamasaki's proof,

Kharazishvili independently solved the same problem in 1984. In 1991 Kirtadze and Pantsulaia [\[KP1\]](#page-62-6) provided yet another solution (see also Pantsulaia [\[PA\]](#page-62-7)). Finally, In 2007, Kirtadze and Pantsulaia proved that, if µ is the completion of the measure µ, then: (see [\[KP2\]](#page-62-1))

Theorem 0.5. The measure µ is the unique regular σ-finite measure on R<sup>∞</sup> (uniqueness up to a mulitplicative constant), which is assigns the value one to the set [− 1 2 , 1 2 ℵo , is invariant under translations from the group ℓ<sup>1</sup> and has the metrically transitivity property with respect to ℓ1.

In the mean time, in 1991 Baker [\[BA1\]](#page-60-3), unaware of the Elliott-Morse measures, dropped the requirement that the measure be σ-finite and constructed a translation invariant measure, ν, on R<sup>∞</sup> (see also Baker (2004), [\[BA2\]](#page-61-4)). In 1992, Ritter and Hewitt [\[RH\]](#page-63-3) constructed a translation invariant measure related to that of Elliott Morse.

Starting in 2007, A. M. Vershik (see [\[V1\]](#page-63-4), [\[V2\]](#page-63-5), [\[V3\]](#page-63-6) and references contained therein) started an investigation of an infinite-dimensional analogue of Lebesgue measure that is constructed in a different manner than that studied in the previous papers. Roughly stated, he considers the weak limit as n → ∞ of invariant measures on certain homogeneous spaces (hypersurfaces of high dimension) of the Cartan subgroup of the Lie groups SL(n, R) (i.e., the subgroups of diagonal matrices with unit determinant). Vershik's measure is also unique and invariant under the multiplicative group of positive functions, suggesting that a logarithmic transformation may lead to a version of the measure in this paper. (The paper of Vandev [\[VA\]](#page-64-3) should also be consulted.)

<span id="page-9-0"></span>Purpose. The purpose of this paper is to show that a minor change in the way we represent R<sup>∞</sup> makes it possible to construct a σ-finite regular version of Lebesgue measure using basic methods of measure theory from R n . Since the measure is regular, it turns out to be the Kirtadze and Pantsulaia [\[KP1\]](#page-62-6) measure, which is unique (see Theorem 0.5). Using our approach, we construct an analogue of both Lebesgue and Gaussian measure (countably additive) on every (classical) separable Banach space with a Schauder basis. The version of Gaussian measure constructed is also rotationally invariant (a property not shared by Wiener measure). This approach also allows us to satisfy all the requirements of a useful infinite dimensional theory.

<span id="page-9-1"></span>Summary. In the first section, we show how von Neumann's infinite tensor product Hilbert space theory implies that a natural version of Lebesgue measure must exist on R<sup>∞</sup> and points to a possible approach. In the first part of Section 2, we show that a slight change in thinking about the cause for problems with unbounded measures on R<sup>∞</sup> makes the construction of Lebesgue measure not only possible, but no more difficult then the same construction on R n . (We denote it by R<sup>∞</sup> I , for reasons that are discussed in this section.) We also provide natural analogues of Lebesgue and Gaussian measure for every separable Banach space with a Schauder basis and show that ℓ<sup>1</sup> is the maximal translation invariant subspace. In the last part of Section 2, we show that ℓ<sup>2</sup> is the maximal rotation invariant subspace. In Section 3, we study the convergence properties of infinite sums and products of bounded and unbounded linear operators. In Section 4, we investigate some of the function spaces over R<sup>∞</sup> I and in Section 5, we discuss Fourier transforms and Pontryagin duality theory for Banach spaces. A major result is that there are two different extensions of the Pontrjagin Duality theory for infinite dimensional spaces. In this section, we also show that our theory allows us to extend Young's inequality to ever separable Banach space with a Schauder basis. In Section 6, we give some constructive examples of partial differential operators in infinitely many variables. This allows us to briefly discuss the famous partial differential equation derived by Phillip Duncan Thompson [\[PDT\]](#page-63-0), on infinite-dimensional phase space to represent an ensemble of randomly forced two-dimensional viscous flows.

# 1. Why λ<sup>∞</sup> Must Exist

<span id="page-10-0"></span>In order to see that some reasonable version of Lebesgue measure must exist, we need to review von Neumann's infinite tensor product Hilbert space theory [\[VN2\]](#page-64-4). To do this, we first define infinite products of complex numbers. (There are a number of other possibilities, see [\[GU\]](#page-61-5) and [\[PA\]](#page-62-7), pg. 272-274.) In order to avoid trivialities, we always assume that, in any product, all terms are nonzero.

**Definition 1.1.** If  $\{z_i\}$  is a sequence of complex numbers indexed by  $i \in \mathbb{N}$  (the natural numbers),

- (1) We say that the product  $\prod_{i\in\mathbb{N}} z_i$  is convergent with limit z if, for every  $\varepsilon > 0$ , there is a finite set  $J(\varepsilon)$  such that, for all finite sets  $J \subset \mathbb{N}$ , with  $J(\varepsilon) \subset J$ , we have  $\left|\prod_{i\in J} z_i z\right| < \varepsilon$ .
- (2) We say that the product  $\prod_{i\in\mathbb{N}} z_i$  is quasi-convergent if  $\prod_{i\in\mathbb{N}} |z_i|$  is convergent. (If the product is quasi-convergent, but not convergent, we assign it the value zero.)

We note that

$$(1.1) 0 < \left| \prod_{i \in \mathbb{N}} z_i \right| < \infty \text{ if and only if } \sum_{i \in \mathbb{N}} |1 - z_i| < \infty.$$

Let  $\mathcal{H}_i = L^2[\mathbb{R}, \lambda]$  for each  $i \in \mathbb{N}$  and let  $\mathcal{H}^2_{\otimes} = \hat{\otimes}_{i=1}^{\infty} L^2[\mathbb{R}, \lambda]$  be the infinite tensor product of von Neumann. To see what this object looks like:

**Definition 1.2.** Let  $g = \underset{i \in \mathbb{N}}{\otimes} g_i$  and  $h = \underset{i \in \mathbb{N}}{\otimes} h_i$  be in  $\mathcal{H}^2_{\otimes}$ .

- (1) We say that g is strongly equivalent to h  $(g \equiv^s h)$  if and only if  $\sum_{i \in \mathbb{N}} |1 \langle g_i, h_i \rangle_i| < \infty \ .$
- (2) We say that g is weakly equivalent to h  $(g \equiv^w h)$  if and only if  $\sum_{i \in \mathbb{N}} |1 |\langle g_i, h_i \rangle_i| | < \infty.$

Proofs of the following may be found in von Neumann [VN2] (see also [GZ], [GZ1]).

**Lemma 1.3.** We have  $g \equiv^w h$  if and only if there exist  $z_i$ ,  $|z_i| = 1$ , such that  $\underset{i \in \mathbb{N}}{\otimes} z_i g_i \equiv^s \underset{i \in \mathbb{N}}{\otimes} h_i$ .

**Theorem 1.4.** The relations defined above are equivalence relations on  $\mathcal{H}^2_{\otimes}$ , which decomposes  $\mathcal{H}^2_{\otimes}$  into disjoint equivalence classes (orthogonal subspaces).

**Definition 1.5.** For  $g = \underset{i \in \mathbb{N}}{\otimes} g_i \in \mathcal{H}^2_{\otimes}$ , we define  $\mathcal{H}^2_{\otimes}(g)$  to be the closed subspace generated by the span of all  $h \equiv^s g$  and we call it the strong partial tensor product space generated by the vector g. (von Neumann called it an incomplete tensor product space.)

**Theorem 1.6.** For the partial tensor product spaces, we have the following:

- (1) If  $h_i \neq g_i$  occurs for at most a finite number of i, then  $h = \underset{i \in \mathbb{N}}{\otimes} h_i \equiv^s$  $g = \underset{i \in \mathbb{N}}{\otimes} g_i.$
- (2) The space  $\mathcal{H}^2_{\otimes}(g)$  is the closure of the linear span of  $h = \underset{i \in \mathbb{N}}{\otimes} h_i$  such that  $h_i \neq g_i$  occurs for at most a finite number of i.
- (3) If  $g = \bigotimes_{i \in \mathbb{N}} g_i$  and  $h = \bigotimes_{i \in \mathbb{N}} h_i$  are in different equivalence classes of  $\mathcal{H}^2_{\otimes}$ , then  $(g,h)_{\otimes} = \prod_{i \in \mathbb{N}} \langle g_i, h_i \rangle_i = 0$ .
- $(4) \ \mathcal{H}^2_{\otimes}(g)^w = \bigoplus_{h \equiv w_q} \left[ \mathcal{H}^2_{\otimes}(h)^s \right].$
- (5) For each g,  $\mathcal{H}^2_{\otimes}(g)^s$  is a separable Hilbert space.
- (6) For each g,  $\mathcal{H}^2_{\otimes}(g)^w$  is not a separable Hilbert space.

It follows from (6) that  $\mathcal{H}^2_{\otimes} = \hat{\otimes}_{i=1}^{\infty} L^2[\mathbb{R}, \lambda]$  is not a separable Hilbert space.

From (5), we see that it is reasonable to define  $L^2[\mathbb{R}^\infty, \lambda_\infty] = \mathcal{H}^2_{\otimes}(h)^s$ , for some  $h = \bigotimes_{i=1}^\infty h_i$ . This definition is ambiguous, but, in most applications, the particular version does not matter. To remove the ambiguity, we should identify a canonical version of  $h = \bigotimes_{i=1}^\infty h_i$ . Any reasonable version of  $\lambda_\infty$  should satisfy  $\lambda_\infty(I_0) = 1$ , where  $I = [\frac{-1}{2}, \frac{1}{2}]$  and  $I_0 = \times_{i=1}^\infty I$ .

<span id="page-13-0"></span>**Definition 1.7.** If  $\chi_I$  is the indicator function for I and  $h_i = \chi_I$ , we set  $h = \bigotimes_{i=1}^{\infty} h_i$ . We define the canonical version of  $L^2[\mathbb{R}^{\infty}, \lambda_{\infty}] = L^2[\mathbb{R}^{\infty}, \lambda_{\infty}](h)^s$ .

# 2. Lebesgue Measure on $\mathbb{R}_I^{\infty}$

<span id="page-13-1"></span>2.1. The Construction. We now have the problem of identifying the measure space associated with  $L^2[\mathbb{R}^\infty, \lambda_\infty](h)^s$ . In the historical approach to the construction of infinite products of measures  $\{\mu_k, k \in \mathbb{N}\}$  on  $\mathbb{R}^\infty$ , the chosen topology defines open sets to be the (cartesian) product of an arbitrary finite number of open sets in  $\mathbb{R}$ , while the remaining infinite number are copies of  $\mathbb{R}$  (cylindrical sets). The success of Kolmogorov's work on the foundations of probability theory naturally led to the condition that  $\mu_k(\mathbb{R})$  be finite for all but a finite number of k (see [KO]). Thus, any attempt to construct Lebesgue measure via this approach starts out a failure in the beginning. However, Kolmogorov's approach is not the only way to induce a total measure of one for the spaces under consideration.

Our definition of the canonical version of  $L^2[\mathbb{R}^{\infty}, \lambda_{\infty}]$  offers another approach. To see how, consider a simple extension of the theory on  $\mathbb{R}$ . Let

 $I = [-\frac{1}{2}, \frac{1}{2}]$  and define  $\mathbb{R}_I = \mathbb{R} \times I_1$ , where  $I_1 = \underset{i=2}{\overset{\infty}{\times}} I$ . If  $\mathfrak{B}(\mathbb{R})$  is the Borel  $\sigma$ -algebra for  $\mathbb{R}$ , let  $\mathfrak{B}(\mathbb{R}_I)$  be the Borel  $\sigma$ -algebra for  $\mathbb{R}_I$ . For each set  $A \in \mathfrak{B}(\mathbb{R})$  with  $\lambda(A) < \infty$ , let  $A_I$  be the corresponding set in  $\mathfrak{B}(\mathbb{R}_I)$ ,  $A_I = A \times I_1$ . We define  $\lambda_{\infty}(A_I)$  by:

$$\lambda_{\infty}(A_I) = \lambda(A) \times \prod_{i=2}^{\infty} \lambda(I) = \lambda(A).$$

We can construct a theory of Lebesgue measure on  $\mathbb{R}_I$  that completely parallels that on  $\mathbb{R}$ . This suggests that we use Lebesgue measure and replace the (tail end of the) infinite product of copies of  $\mathbb{R}$  by infinite products of copies of I. The purpose of this section is to provide such a construction. Since we will be studying unbounded measures, for consistency, we use the following conventions:  $0 \cdot \infty = 0$  and  $0 \cdot \infty^{\infty} = \infty$ .

Recall that  $\mathbb{R}^{\infty}$  is the set of all  $\mathbf{x} = (x_1, x_2, \cdots)$ , where  $x_i \in \mathbb{R}$ . This is a linear space which is not a Banach space. However, it is a complete metric space with metric given by:

$$d(\mathbf{x}, \mathbf{y}) = \sum_{n=1}^{\infty} \frac{1}{2^n} \frac{|x_n - y_n|}{1 + |x_n - y_n|}.$$

**Remark 2.1.**  $\mathbb{R}^{\infty}$  is a special case of a Polish space, which Banach called a Fréchet space i.e., a Polish space with a translation invariant metric (see Banach [BA]). The topology generated by  $d(\cdot,\cdot)$  is generally known as the Tychonoff topology.

For each n, define  $\mathbb{R}^n_I = \mathbb{R}^n \times I_n$ , where  $I_n = \underset{i=n+1}{\overset{\infty}{\times}} I$ .

**Definition 2.2.** If  $A_n = A \times I_n$ ,  $B_n = B \times I_n$  are any sets in  $\mathbb{R}^n_I$ , then we define:

- (1)  $A_n \cup B_n = A \cup B \times I_n$ ,
- (2)  $A_n \cap B_n = A \cap B \times I_n$ , and
- (3)  $B_n^c = B^c \times I_n$ .

In order to avoid confusion, we always assume that  $I_0 = \times_{i=1}^{\infty} I \subset \mathbb{R}^1_I$ . We can now define the topology for  $\mathbb{R}^n_I$  via the following class of open sets:

$$\mathfrak{O}_n = \{ U_n : U_n = U \times I_n, \ U \text{ open in } \mathbb{R}^n \}.$$

2.1.1. **Definition of**  $\mathbb{R}_I^{\infty}$ . It is easy to see that  $\mathbb{R}_I^n \subset \mathbb{R}_I^{n+1}$ . Since this is an increasing sequence, we can define  $\mathbb{R}_I'^{\infty}$  by:

$${\mathbb R'}_I^\infty = {\lim}_{n \to \infty} {\mathbb R}_I^n = \bigcup_{k=1}^\infty {\mathbb R}_I^k.$$

Let  $\tau_1$  be the topology on  $\mathbb{R}_I^{\infty} = \mathfrak{X}_1$  induced by the class of open sets  $\mathfrak{O}$  defined by:

$$\mathfrak{O} = \bigcup_{n=1}^{\infty} \mathfrak{O}_n = \bigcup_{n=1}^{\infty} \{ U_n : U_n = U \times I_n, \ U \text{ open in } \mathbb{R}^n \},$$

and let  $\tau_2$  be topology on  $\mathbb{R}^{\infty} \setminus \mathbb{R}_I^{\infty} = \mathfrak{X}_2$  induced by the metric  $d_2$ , for which  $d_2(x,y) = 1$ ,  $x \neq y$  and  $d_2(x,y) = 0$ , x = y, for all  $x, y \in \mathfrak{X}_2$ .

**Definition 2.3.** We define  $(\mathbb{R}_I^{\infty}, \tau)$  to be the sum  $(\mathfrak{X}_1, \tau_1)$  and  $(\mathfrak{X}_2, \tau_2)$ , so that every open set in  $(\mathbb{R}_I^{\infty}, \tau)$  is union of two disjoint sets  $G_1 \cup G_2$ , where  $G_1$  is open in  $(\mathfrak{X}_1, \tau_1)$  and  $G_2$  is open in  $(\mathfrak{X}_2, \tau_2)$ .

It now follows from the above construction that  $\mathbb{R}_I^{\infty} = \mathbb{R}^{\infty}$  as sets. (However, they are not equal as topological spaces.) The following result shows that convergence in the  $\tau$ -topology always implies convergence in the Tychonoff topology.

**Theorem 2.4.** If  $y_k$  converges to x in the  $\tau$ -topology, then  $y_k$  converges to x in the Tychonoff topology.

Proof. Case 1. If  $x \in \mathbb{R}^{\infty} \setminus \mathbb{R}'_{I}^{\infty}$  then there is N such that  $y_{k} = x$  for all k > N. Indeed, for a neighborhood of diameter  $\frac{1}{2}$  about x, there is a N such that  $d_{2}(x, y_{k}) < 1/2$  for all k > N. This means that  $y_{k} = x$  for k > N ( $\{z : d_{2}(x, z) < 1/2\}$  only contains x), so that  $y_{k}$  converges to x in the Tychonoff topology.

Case 2. If  $x \in \mathbb{R}_I^{\infty}$  and  $y_k$  converges to x, then for any neighborhood  $U_n \subset \mathfrak{O}_n$ , there is N such that or all k > N,  $y_k \in U_n$ . This means that  $y_k \in \mathbb{R}_I^{\infty}$  for k > N, so that  $y_k$  converges to x in the Tychonoff topology.  $\square$ 

2.1.2. **Definition of**  $\mathfrak{B}(\mathbb{R}_{I}^{\infty})$ . In a similar manner, if  $\mathfrak{B}(\mathbb{R}_{I}^{n})$  is the Borel  $\sigma$ -algebra for  $\mathbb{R}_{I}^{n}$  (i.e., the smallest  $\sigma$ -algebra generated by the  $\mathfrak{O}_{n}$ ), then  $\mathfrak{B}(\mathbb{R}_{I}^{n}) \subset \mathfrak{B}(\mathbb{R}_{I}^{n+1})$ , so we can define  $\mathfrak{B}'(\mathbb{R}_{I}^{\infty})$  by:

$$\mathfrak{B}'(\mathbb{R}_I^\infty) = \mathrm{lim}_{n \to \infty} \mathfrak{B}(\mathbb{R}_I^n) = \mathop{\cup}\limits_{k=1}^\infty \mathfrak{B}(\mathbb{R}_I^k).$$

If  $\mathcal{P}(\cdot)$  denotes a powerset of a set (i.e.,  $\mathcal{P}(A) = \{X : X \subseteq A\}$ ), let  $\mathfrak{B}(\mathbb{R}_I^{\infty})$  be the smallest  $\sigma$ -algebra containing  $\mathfrak{B}'(\mathbb{R}_I^{\infty}) \cup \mathcal{P}(R_I^{\infty} \setminus \cup_{n=1}^{\infty} \mathbb{R}_I^n)$ . (It is obvious that the class  $\mathfrak{B}(\mathbb{R}_I^{\infty})$  coincides with Borel  $\sigma$ -algebra generated by

the  $\tau$ -topology on  $\mathbb{R}^{\infty}$ .) From our definition of  $\mathfrak{B}(\mathbb{R}_{I}^{\infty})$  we see that  $\mathfrak{B}(\mathbb{R}^{\infty}) \subset \mathfrak{B}(\mathbb{R}_{I}^{\infty})$  and the containment is proper.

**Theorem 2.5.**  $\lambda_{\infty}(\cdot)$  is a measure on  $\mathfrak{B}(\mathbb{R}^n_I)$ , equivalent to n-dimensional Lebesgue measure on  $\mathbb{R}^n$ .

*Proof.* If  $A = \underset{i=1}{\overset{\infty}{\times}} A_i \in \mathfrak{B}(\mathbb{R}^n_I)$ , then  $\lambda(A_i) = 1$  for i > n so that the series  $\lambda_{\infty}(A) = \prod_{i=1}^{\infty} \lambda(A_i)$  always converges. Furthermore,

$$(2.1) 0 < \lambda_{\infty}(A) = \prod_{i=1}^{\infty} \lambda(A_i) = \prod_{i=1}^{n} \lambda(A_i) = \lambda_n (\underset{i=1}{\overset{n}{\times}} A_i).$$

Since sets of the type  $A = \underset{i=1}{\overset{n}{\times}} A_i$  generate  $\mathfrak{B}(\mathbb{R}^n)$ , we see that  $\lambda_{\infty}(\cdot)$ , restricted to  $\mathbb{R}^n_I$ , is equivalent to  $\lambda_n(\cdot)$ .

<span id="page-17-0"></span>Corollary 2.6. The measure  $\lambda_{\infty}(\cdot)$  is both translationally and rotationally invariant on  $(\mathbb{R}^n_I, \mathfrak{B}[\mathbb{R}^n_I])$ .

2.2. The Extension to  $\mathbb{R}_I^{\infty}$ . It is not obvious that  $\lambda_{\infty}(\cdot)$  can be extended to a countably additive measure on  $\mathfrak{B}(\mathbb{R}_I^{\infty})$ .

#### Definition 2.7. Let

$$\Delta_0 = \{K_n = K \times I_n \in \mathfrak{B}(\mathbb{R}_I^n) \subset \mathfrak{B}(\mathbb{R}_I^\infty) : n \in \mathbb{N}, K \text{ is compact and } 0 < \lambda_\infty(K_n) < \infty\},$$

$$\Delta = \{P_N = \bigcup_{i=1}^N K_{n_i}, N \in \mathbb{N}; K_{n_i} \in \Delta_0 \text{ and } \lambda_\infty(K_{n_l} \cap K_{n_m}) = 0, l \neq m\}.$$

**Definition 2.8.** If  $P_N \in \Delta$ , we define

$$\lambda_{\infty}(P_N) = \sum_{i=1}^{N} \lambda_{\infty}(K_{n_i}).$$

Since  $P_N \in \mathfrak{B}(\mathbb{R}^n_I)$  for some n, and  $\lambda_{\infty}(\cdot)$  is a measure on  $\mathfrak{B}(\mathbb{R}^n_I)$ , the next result follows:

**Lemma 2.9.** If  $P_{N_1}, P_{N_2} \in \Delta$  then:

(1) If 
$$P_{N_1} \subset P_{N_2}$$
, then  $\lambda_{\infty}(P_{N_1}) \leq \lambda_{\infty}(P_{N_2})$ .

(2) If 
$$\lambda_{\infty}(P_{N_1} \cap P_{N_2}) = 0$$
, then  $\lambda_{\infty}(P_{N_1} \cup P_{N_2}) = \lambda_{\infty}(P_{N_2}) + \lambda_{\infty}(P_{N_2})$ .

**Definition 2.10.** If  $G \subset \mathbb{R}_I^{\infty}$  is any open set, we define:

$$\lambda_{\infty}(G) = \lim_{N \to \infty} \sup \{\lambda_{\infty}(P_N) : P_N \in \Delta, P_N \subset G, \}.$$

**Theorem 2.11.** If  $\mathfrak{O}$  is the class of open sets in  $\mathfrak{B}(\mathbb{R}_I^{\infty})$ , we have:

- $(1) \ \lambda_{\infty}(\mathbb{R}_{I}^{\infty}) = \infty.$
- (2) If  $G_1, G_2 \in \mathfrak{O}, G_1 \subset G_2$ , then  $\lambda_{\infty}(G_1) \leq \lambda_{\infty}(G_2)$ .
- (3) If  $\{G_k\} \subset \mathfrak{O}$ , then

$$\lambda_{\infty}(\bigcup_{k=1}^{\infty} G_k) \le \sum_{k=1}^{\infty} \lambda_{\infty}(G_k).$$

(4) If the  $G_k$  are disjoint, then

$$\lambda_{\infty}(\bigcup_{k=1}^{\infty} G_k) = \sum_{k=1}^{\infty} \lambda_{\infty}(G_k).$$

*Proof.* The proof of (1) is standard. To prove (2), observe that

$$\{P_N: P_N \subset G_1\} \subset \{P_N': P_N' \subset G_2\},$$

so that  $\lambda_{\infty}(G_1) \leq \lambda_{\infty}(G_2)$ . To prove (3), let  $P_N \subset \bigcup_{k=1}^{\infty} G_k$ . Since  $P_N$  is compact, there is a finite number of the  $G_k$  which cover  $P_N$ , so that

P<sup>N</sup> ⊂ S<sup>L</sup> <sup>k</sup>=1 Gk. Now, for each Gk, there is a PN<sup>k</sup> ⊂ Gk. Furthermore, as P<sup>N</sup> is arbitrary, we can assume that P<sup>N</sup> = P ′ <sup>N</sup> = S<sup>L</sup> <sup>k</sup>=1 PN<sup>k</sup> . Since there is an n such that all PN<sup>k</sup> ∈ B(R n I ), we may also assume that λ∞(PN<sup>l</sup> ∩PNm) = 0, l 6= m. We now have that

$$\lambda_{\infty}(P_N) = \sum_{k=1}^{L} \lambda_{\infty}(P_{N_k}) \leqslant \sum_{k=1}^{L} \lambda_{\infty}(G_k) \leqslant \sum_{k=1}^{\infty} \lambda_{\infty}(G_k).$$

It follows that

$$\lambda_{\infty}(\bigcup_{k=1}^{\infty} G_k) \le \sum_{k=1}^{\infty} \lambda_{\infty}(G_k).$$

If the G<sup>k</sup> are disjoint, observe that if P<sup>N</sup> ⊂ P ′M,

$$\lambda_{\infty}(P'_{M}) \ge \lambda_{\infty}(P_{N}) = \sum_{k=1}^{L} \lambda_{\infty}(P_{N_{k}}).$$

It follows that

$$\lambda_{\infty}(\bigcup_{k=1}^{\infty} G_k) \ge \sum_{k=1}^{L} \lambda_{\infty}(G_k).$$

This is true for all L so that this, combined with (3), gives our result.

If F is an arbitrary compact set in B(R<sup>∞</sup> I ), we define

(2.2) 
$$\lambda_{\infty}(F) = \inf \{ \lambda_{\infty}(G) : F \subset G, G \text{ open} \}.$$

Remark 2.12. At this point we see the power of B(R<sup>∞</sup> I ). Unlike B(R∞), equation (2.2) is well-defined for B(R<sup>∞</sup> I ) because it has a sufficient number of open sets of finite measure.

Theorem 2.13. Equation (2.2) is consistent with Definition 2.6 and the results of Theorem 2.11.

**Definition 2.14.** Let A be an arbitrary set in  $\mathbb{R}^{\infty}_{I}$ .

(1) The outer measure (on  $\mathbb{R}_I^{\infty}$ ) is defined by:

$$\lambda_{\infty}^*(A) = \inf \{ \lambda_{\infty}(G) : A \subset G, G \text{ open} \}.$$

We let  $\mathfrak{L}_0$  be the class of all A with  $\lambda_{\infty}^*(A) < \infty$ .

(2) If  $A \in \mathfrak{L}_0$ , we define the inner measure of A by

$$\lambda_{\infty,(*)}(A) = \sup \{\lambda_{\infty}(F) : F \subset A, F \text{ compact}\}.$$

(3) We say that A is a bounded measurable set if  $\lambda_{\infty}^*(A) = \lambda_{\infty,(*)}(A)$ , and define the measure of A,  $\lambda_{\infty}(A)$ , by  $\lambda_{\infty}(A) = \lambda_{\infty}^*A$ .

**Theorem 2.15.** Let A, B and  $\{A_k\}$  be arbitrary sets in  $\mathbb{R}_I^{\infty}$  with finite outer measure.

- (1)  $\lambda_{\infty,(*)}(A) \leq \lambda_{\infty}^{*}(A)$ .
- (2) If  $A \subset B$  then  $\lambda_{\infty}^*(A) \leq \lambda_{\infty}^*(B)$  and  $\lambda_{\infty,(*)}(A) \leq \lambda_{\infty,(*)}(B)$ .
- (3)  $\lambda_{\infty}^* (\bigcup_{k=1}^{\infty} A_k) \leq \sum_{k=1}^{\infty} \lambda_{\infty}^* (A_k).$
- (4) If the  $\{A_k\}$  are disjoint,  $\lambda_{\infty,(*)}(\bigcup_{k=1}^{\infty} A_k) \geq \sum_{k=1}^{\infty} \lambda_{\infty,(*)}(A_k)$ .

*Proof.* The proofs of (1) and (2) are straightforward. To prove (3), let  $\varepsilon > 0$  be given. Then, for each k, there exists an open set  $G_k$  such that  $A_k \subset G_k$  and  $\lambda_{\infty}(G_k) < \lambda_{\infty}^*(A_k) + \varepsilon 2^{-k}$ . Since  $(\bigcup_{k=1}^{\infty} A_k) \subset (\bigcup_{k=1}^{\infty} G_k)$ , we have

$$\lambda_{\infty}^{*} \left( \bigcup_{k=1}^{\infty} A_{k} \right) \leqslant \lambda_{\infty} \left( \bigcup_{k=1}^{\infty} G_{k} \right) \leqslant \sum_{k=1}^{\infty} \lambda_{\infty}(G_{k})$$
$$< \sum_{k=1}^{\infty} \left[ \lambda_{\infty}^{*}(A_{k}) + \varepsilon 2^{-k} \right] = \sum_{k=1}^{\infty} \lambda_{\infty}^{*}(A_{k}) + \varepsilon.$$

Since  $\varepsilon$  is arbitrary, we are done.

To prove (4), let  $F_1, F_2, \ldots, F_N$  be compact subsets of  $A_1, A_2, \ldots, A_N$ , respectively. Since the  $A_k$  are disjoint,

$$\lambda_{\infty,(*)}\left(\bigcup_{k=1}^{\infty}A_{k}\right)\geqslant\lambda_{\infty}\left(\bigcup_{k=1}^{N}F_{k}\right)=\sum\nolimits_{k=1}^{N}\lambda_{\infty}(F_{k}).$$

Thus,

$$\lambda_{\infty,(*)}\left(\bigcup_{k=1}^{\infty} A_k\right) \ge \sum_{k=1}^{N} \lambda_{\infty,(*)}(A_k).$$

Since N is arbitrary, we are done.

The next two important theorems follow from the last one.

**Theorem 2.16.** (Regularity) If A has finite measure, then for every  $\varepsilon > 0$  there exist a compact set F and an open set G such that  $F \subset A \subset G$ , with  $\lambda_{\infty}(G \setminus F) < \varepsilon$ .

*Proof.* Let  $\varepsilon > 0$  be given. Since A has finite measure, it follows from our definitions of  $\lambda_{\infty,(*)}$  and  $\lambda_{\infty}^*$  that there is a compact set  $F \subset A$  and an open set  $G \supset A$  such that

$$\lambda_{\infty}(G) < \lambda_{\infty}^{*}(A) + \frac{\epsilon}{2} \quad \text{and} \quad \lambda_{\infty}(F) > \lambda_{\infty,(*)}(A) - \frac{\epsilon}{2}.$$

Since  $\lambda_{\infty}(G) = \lambda_{\infty}(F) + \lambda_{\infty}(G \setminus F)$ , we have:

$$\lambda_{\infty}(G \setminus F) = \lambda_{\infty}(G) - \lambda_{\infty}(F) < (\lambda_{\infty}(A) + \frac{\varepsilon}{2}) - (\lambda_{\infty}(A) - \frac{\varepsilon}{2}) = \varepsilon.$$

**Theorem 2.17.** (Countable Additivity) If the family  $\{A_k\}$  consists of disjoint sets with bounded measure and  $A = \bigcup_{k=1}^{\infty} A_k$ , with  $\lambda_{\infty}^*(A) < \infty$ . then  $\lambda_{\infty}(A) = \sum_{k=1}^{\infty} \lambda_{\infty}(A_k)$ .

*Proof.* Since  $\lambda_{\infty}^*(A) < \infty$ , we have:

$$\lambda_{\infty}^*(A) \leqslant \sum_{k=1}^{\infty} \lambda_{\infty}^*(A_k) = \sum_{k=1}^{\infty} \lambda_{\infty,(*)}(A_k) \leqslant \lambda_{\infty,(*)}(A) \leqslant \lambda_{\infty}^*(A).$$

It follows that  $\lambda_{\infty}(A) = \lambda_{\infty}^{*}(A) = \lambda_{\infty,(*)}(A)$ , so that

$$\lambda_{\infty}(A) = \lambda_{\infty} \left( \bigcup_{k=1}^{\infty} A_k \right) = \sum_{k=1}^{\infty} \lambda_{\infty}(A_k).$$

**Definition 2.18.** Let A be an arbitrary set in  $\mathbb{R}_I^{\infty}$ . We say that A is measurable if  $A \cap M \in \mathfrak{L}_0$  for all  $M \in \mathfrak{L}_0$ . In this case, we define  $\lambda_{\infty}(A)$  by:

$$\lambda_{\infty}(A) = \sup \left\{ \lambda_{\infty}(A \cap M) : M \subset \mathfrak{L}_0 \right\}.$$

We let  $\mathfrak{L}_I^{\infty}$  be the class of all measurable sets A.

Proofs of the following results are standard (see Jones [J], pages 48-52).

**Theorem 2.19.** Let A and  $\{A_k\}$  be arbitrary sets in  $\mathfrak{L}_I^{\infty}$ .

- (1) If  $\lambda_{\infty}^*(A) < \infty$ , then  $A \in \mathfrak{L}_0$  if and only if  $A \in \mathfrak{L}_I^{\infty}$ . In this case,  $\lambda_{\infty}(A) = \lambda_{\infty}^*(A)$ .
- (2)  $\mathfrak{L}_I^{\infty}$  is closed under countable unions, countable intersections, differences and complements.

(3)

$$\lambda_{\infty}(\bigcup_{k=1}^{\infty} A_k) \le \sum_{k=1}^{\infty} \lambda_{\infty}(A_k).$$

(4) If {Ak} are disjoint,

$$\lambda_{\infty}(\bigcup_{k=1}^{\infty} A_k) = \sum_{k=1}^{\infty} \lambda_{\infty}(A_k).$$

(5) If A<sup>k</sup> ⊂ Ak+1 for all k, then

$$\lambda_{\infty}(\bigcup_{k=1}^{\infty} A_k) = \lim_{k \to \infty} \lambda_{\infty}(A_k).$$

(6) If Ak+1 ⊂ A<sup>k</sup> for all k and λ∞(A1) < ∞, then

$$\lambda_{\infty}(\bigcap_{k=1}^{\infty} A_k) = \lim_{k \to \infty} \lambda_{\infty}(A_k).$$

We end this section with an important result that relates Borel sets to L∞ I -measurable sets (Lebesgue).

Theorem 2.20. Let A be a L<sup>∞</sup> I -measurable set. Then there exists a Borel set F and a set N with λ∞(N) = 0 such that A = F ∪ N.

Thus, we see that λ∞(·) is a regular countably additive σ-finite Borel measure on R<sup>∞</sup> <sup>I</sup> = R<sup>∞</sup> (as sets). More important is the fact that the development is no more difficult than the corresponding theory for Lebesgue measure on R n .

Throughout the remainder of the paper we will also use B[R<sup>∞</sup> I ] for its completion L<sup>∞</sup> <sup>I</sup> when convenient. This should cause no confusion since the given context will always be clear.

<span id="page-24-0"></span>2.3. Separable Banach Spaces. In order to see what other advantages our construction of  $(\mathbb{R}_I^{\infty}, \mathfrak{B}[\mathbb{R}_I^{\infty}], \lambda_{\infty}(\cdot))$  offers, in this section we study separable Banach spaces. Let  $\mathcal{B}$  be any separable Banach space.

Recall that (see Diestel [DI], page 32):

**Definition 2.21.** A sequence  $(u_n)$  is called a Schauder basis for  $\mathcal{B}$  if  $||u_n||_{\mathcal{B}} = 1$  and, for each  $f \in \mathcal{B}$ , there is a unique sequence  $(a_n)$  of scalars such that

$$f = \lim_{n \to \infty} \sum_{k=1}^{n} a_k u_k.$$

**Definition 2.22.** A sequence  $(v_n)$  is called an absolutely convergent Schauder basis for  $\mathcal{B}$  if  $\sum_{n=1}^{\infty} \|v_n\|_{\mathcal{B}} < \infty$  and, for each  $f \in \mathcal{B}$ , there is a unique sequence  $(b_n)$  of scalars such that

$$f = \lim_{n \to \infty} \sum_{k=1}^{n} b_k v_k.$$

**Lemma 2.23.** Let  $(u_n)$  be a Schauder basis for  $\mathcal{B}$ , then there exists an absolutely convergent Schauder basis for  $\mathcal{B}$ .

*Proof.* Let  $(v_n) = (\frac{u_n}{2^n})$ . Then

$$\sum_{n=1}^{\infty} \|v_n\|_{\mathcal{B}} = \sum_{n=1}^{\infty} \frac{\|u_n\|_{\mathcal{B}}}{2^n} = \sum_{n=1}^{\infty} \frac{1}{2^n} = 1 < \infty.$$

To see that  $(v_n)$  is a Schauder basis for  $\mathcal{B}$ , let  $f \in \mathcal{B}$ . By definition, there is a unique sequence  $(a_n)$  of scalars such that

$$f = \lim_{n \to \infty} \sum_{k=1}^{n} a_k u_k.$$

If we take the sequence  $(b_n) = (2^n a_n)$ , then

$$\lim_{n\to\infty} \sum_{k=1}^n b_k v_k = \lim_{n\to\infty} \sum_{k=1}^n a_k u_k = f.$$

It is known that most of the natural separable Banach spaces, and all that have any use for applications in analysis, have a Schauder basis. In particular, it is easy to see from the definition of a Schauder basis that, for any sequence  $(a_n) \in \mathbb{R}_I^{\infty}$  representing a function  $f \in \mathcal{B}$ , we have  $\lim_{n \to \infty} a_n = 0$ . It follows that every separable Banach space (with a Schauder basis) is isomorphic to a subspace of  $\mathbb{R}_I^{\infty}$ .

Let  $\mathcal{B}_I$  be the set of all sequences  $(a_n)$  for which  $\lim_{n\to\infty} \sum_{k=1}^n a_k u_k$  exists in  $\mathcal{B}$ . Define

$$\|(a_n)\|_{\mathcal{B}_I} = \sup_n \left\| \sum_{k=1}^n a_k u_k \right\|_{\mathcal{B}}.$$

Lemma 2.24. An operator

$$T: (\mathcal{B}, ||\cdot||_{\mathcal{B}}) \to (\mathcal{B}_I, ||\cdot||_{\mathcal{B}_I}),$$

defined by  $T(f) = (a_k)$  for  $f = \lim_{n \to \infty} \sum_{k=1}^n a_k u_k \in \mathcal{B}$ , is an isomorphism from  $\mathcal{B}$  onto  $\mathcal{B}_I$ .

Let  $\mathcal{B}$  be a separable Banach space with a Schauder basis and let  $\mathcal{B}_I = T[\mathcal{B}]$ . If  $\mathfrak{B}(\mathcal{B}_I) = \mathcal{B}_I \cap \mathfrak{B}[\mathbb{R}_I^{\infty}]$ , we define the  $\sigma$ -algebra generated on  $\mathcal{B}$ , and associated with  $\mathfrak{B}(\mathcal{B}_I)$  by:

$$\mathfrak{B}_{I}[\mathcal{B}] = \left\{ T^{-1}(A) \mid A \in \mathfrak{B}[\mathcal{B}_{I}] \right\} =: T^{-1} \left\{ \mathfrak{B}\left[\mathcal{B}_{I}\right] \right\}.$$

Note that, just as  $\mathfrak{B}[\mathbb{R}^{\infty}] \subset \mathfrak{B}[\mathbb{R}_{I}^{\infty}]$ , we also have  $\mathfrak{B}[\mathcal{B}] \subset \mathfrak{B}_{I}[\mathcal{B}]$  (with the containment proper).

**Theorem 2.25.** Let  $A \in \mathfrak{B}_I(\mathcal{B})$  and set  $\hat{\lambda}_{\mathcal{B}}(A) = \lambda_{\infty}[T(A)]$ . Let  $\lambda_{\mathcal{B}}$  be the completion of  $\hat{\lambda}_{\mathcal{B}}$ , then  $\lambda_{\mathcal{B}}$  is a non-zero  $\sigma$ -finite Borel measure on  $\mathcal{B}$ .

*Proof.* Let  $\{v_k\}$  be an absolutely convergent Schauder basis. We first prove that, for any L > 0 and any sequence  $(a_k) \in [-L, L]^{\aleph_o}$ , the function  $f = \sum_{k=1}^{\infty} a_k v_k \in \mathcal{B}$ . We then prove that  $\lambda_{\mathcal{B}}$  is nonzero.

#### Part 1

Let L be given. Since  $(v_n)$  is an absolutely convergent Schauder basis, given  $\varepsilon > 0$  we can choose N such that  $\sum_{k=N}^{\infty} ||v_k|| < \frac{\epsilon}{L}$ . It follows that, for  $N \leq m \leq n$ , we have

$$\left\| \sum_{k=m}^{n} a_k v_k \right\| \leqslant \sum_{k=m}^{n} \|v_k\| < \varepsilon.$$

Thus, the sequence  $\{f_n\}$ , defined by  $f_n = \sum_{k=1}^n a_k v_k$ , is a Cauchy sequence in  $\mathcal{B}$ . Since  $\mathcal{B}$  is a Banach space, the sequence converges.

### Part 2

To prove that  $\lambda_{\mathcal{B}}$  is nonzero, it suffices to show that  $\lambda_{\mathcal{B}}\left[T^{-1}\left(I_{0}\right)\right]\neq0$ , where  $(I_{0})=\left[-\frac{1}{2},\frac{1}{2}\right]^{\aleph_{o}}$ . First, we note that T is an injective linear map into  $\mathbb{R}_{I}^{\infty}$ , so that  $B=T^{-1}(I_{0})\in\mathfrak{B}_{I}(\mathcal{B})$ . Thus,

$$\lambda_{\mathcal{B}}(B) = \lambda_{\infty} \left[ T\left(T^{-1}(I_0)\right) \right] = \lambda_{\infty}(I_0) = 1.$$

<span id="page-27-0"></span>2.4. **Translations.** In the theorem below, we will provide a new proof that  $\ell_1$  is the largest (dense) group of admissible translations for  $\mathbb{R}_I^{\infty}$ , so necessarily  $\ell_1$  is the largest group of admissible translations for every separable Banach space  $\mathcal{B}$ .

Recall that  $h(x) = \bigotimes_{k=1}^{\infty} h_k(x_k)$ , where  $h_k(x_k) = 1$ , for  $x_k \in [-\frac{1}{2}, \frac{1}{2}]$ . It follows from  $d\nu = hd\lambda_{\infty}$ , that  $\nu$  is absolutely continuous with respect to  $\lambda_{\infty}$ . Thus,  $\nu$  is equivalent to  $\lambda_{\infty}$ . Let  $\mathfrak{T}_{\lambda_{\infty}}$  be the set of admissible translations for  $\mathbb{R}_I^{\infty}$  (i.e.,  $\lambda_{\infty}[A-x] = \lambda_{\infty}[A]$  for all  $A \in \mathfrak{B}[\mathbb{R}_I^{\infty}]$  and  $x \in \mathfrak{T}_{\lambda_{\infty}}$ ).

**Theorem 2.26.** If  $A \in \mathfrak{B}[\mathbb{R}_I^{\infty}]$  then  $\lambda_{\infty}[A-x] = \lambda_{\infty}[A]$  if and only if  $\mathfrak{T}_{\lambda_{\infty}} = \ell_1$ .

*Proof.* Suppose that  $x \in \ell_1$ . Since  $\nu \sim \lambda_{\infty}$ , we have that  $\mathfrak{T}_{\nu} = \mathfrak{T}_{\lambda_{\infty}}$  (see Yamasaki [YA1]). Thus, it suffices to prove that  $\nu[A-x] = \nu[A]$ . By Kakutani's Theorem ([KA], see also [HHK] pg. 116),  $\nu[A-x] \sim \nu[A]$  if and only if

(2.3) 
$$\prod_{k=1}^{\infty} \int_{-\infty}^{\infty} \sqrt{h_k(y_k) h_k(y_k - x_k)} d\lambda(y_k) > 0.$$

Now,

$$\int_{-\infty}^{\infty} \sqrt{h_k(y_k) h_k(y_k - x_k)} d\lambda(y_k) = \int_{\left[-\frac{1}{2}, \frac{1}{2}\right] \cap \left[-\frac{1}{2} + x_k, \frac{1}{2} + x_k\right]} d\lambda(y_k) = (1 - |x_k|)_+,$$

where  $r_+ = max(0,r)$ . Since  $x \in \ell_1$ ,  $\prod_{k=n}^{\infty} (1-|x_k|)_+ > 0$  for n large enough. Thus, equation (2.3) will be satisfied for every  $x \in \ell_1$ , so that  $\ell_1 \subset \mathfrak{T}_{\nu}$ .

Now, suppose that x ∈ Tλ∞, so that λ∞[A − x] = λ∞[A] for all A ∈ B[R<sup>∞</sup> I ]. Thus, for A ∈ B[R n I ], we have

$$\lambda_{\infty} [A - x] = \lambda_n [A_n - x_n] \cdot \prod_{k=n+1}^{\infty} \lambda \left\{ \left[ -\frac{1}{2}, \frac{1}{2} \right] \cap \left[ -\frac{1}{2} - x_k, \frac{1}{2} - x_k \right] \right\}$$

$$= \lambda_n [A_n] \cdot \prod_{k=n+1}^{\infty} \lambda \left\{ \left[ -\frac{1}{2}, \frac{1}{2} \right] \cap \left[ -\frac{1}{2} - x_k, \frac{1}{2} - x_k \right] \right\} = \lambda_n [A_n] \cdot \prod_{k=n+1}^{\infty} (1 - |x_k|)_+.$$

If 
$$A_n = I_n = \times_{k=1}^n [-\frac{1}{2}, \frac{1}{2}]$$
, we have  $1 = \lim_{n \to \infty} \prod_{k=n+1}^{\infty} (1 - |x_k|)_+$ . It follows that  $\sum_{k=1}^{\infty} |x_k| < \infty$ , so that  $x \in \ell_1$ .

In closing, we note that, since λ<sup>∞</sup> is complete and regular, it is metrically transitivity with respect to R<sup>∞</sup> 0 . It follows from Theorem 0.5 that λ<sup>∞</sup> is unique (this comment also applies to λB).

<span id="page-28-0"></span>2.5. Gaussian measure. If we replace Lebesgue measure by the infinite product Gaussian measure, µ∞, on R∞, we get countable additivity but lose rotational invariance. Furthermore, the µ<sup>∞</sup> measure of l<sup>2</sup> is zero. On the other hand, another approach is to use the standard projection method onto finite dimensional subspaces to construct a probability measure directly on l2. In this case, we recover rotational invariance but not translation invariance (and lose countable additivity). The resolution of this problem led to the development of the Wiener measure [\[WSRM\]](#page-64-5) and this is where we are today. A nice discussion of this and related issues can be found in Dunford and Schwartz [\[DS\]](#page-61-10) (see pg. 402).

We now turn to take a look at infinite product Gaussian measure from our new perspective. The canonical Gaussian measure on  $\mathbb{R}$  is defined by:

$$d\mu(x) = \frac{1}{\sqrt{2\pi}} \exp\left\{-\frac{\left|x\right|^2}{2}\right\} d\lambda(x).$$

Recall that  $\mu_{\infty} = \bigotimes_{k=1}^{\infty} \mu$  is countably additive on  $\mathbb{R}^{\infty}$ , but its measure of  $\ell_2$  is zero. If we introduce a scaled version of Gaussian measure on  $\mathbb{R}_I^{\infty}$ , we can resolve this difficulty. We seek a family of variances  $\{\sigma_k^2\}$  such that

$$\mu_{\mathcal{B}}(\mathcal{B}) = \bigotimes_{k=1}^{\infty} \mu_k(T[\mathcal{B}] = 1,$$

where  $\mu_k$  is a linear Gaussian measure on R with parameters  $(0, \sigma_k)$  for  $k \in N$  and  $\mu_B$  is defined by:

$$\mu_{\mathcal{B}}(B) = \bigotimes_{k=1}^{\infty} \mu_k(T[B]),$$

for any Borel subset B of  $\mathcal{B}$ .

**Lemma 2.27.** Let  $\{\sigma_k^2\}$  be a family of variances such that

$$\sum_{k=1}^{\infty} \sigma_k^2 < \infty,$$

then  $\mu_{\mathcal{B}}\left(T^{-1}([-L,L]^{\aleph_o})\right) > 0$  for every positive number L.

*Proof.* Let  $\{X_k\}$  be the family of independent Gaussian random variables defined on some common probability space,  $(\Omega, \mathfrak{B}, P[\cdot])$ , with law  $\mu_k$ . If  $X = (X_1, X_2, \dots)$ , then

$$\begin{split} &P\left[\left\{\omega\in\Omega\,|\,X(\omega)\in[-L,L]^{\aleph_o}\right\}\right] = P\left[\bigcap_{k=1}^{\infty}\left\{\omega\in\Omega\,|\,X_k(\omega)\in[-L,L]\right\}\right]\\ &= \prod_{k=1}^{\infty}P\left[\left\{\omega\in\Omega\,|\,|X_k(\omega)|\leqslant L\right]\right\}\right] \geqslant \prod_{k=1}^{\infty}\left(1-\frac{\sigma_k^2}{L^2}\right), \quad \text{by Chebyshev's inequality.} \end{split}$$

Clearly the product is positive. We are done since  $B = T^{-1}([-L, L]^{\aleph_o}) \in \mathfrak{B}(\mathcal{B})$  and

$$\mu_{\mathcal{B}}\left(B\right) = \left(\otimes_{k=1}^{\infty} \mu_{k}\right) \left(T[T^{-1}([-L,L]^{\aleph_{o}})]\right) = P\left[\left\{\omega \in \Omega \mid X(\omega) \in [-L,L]^{\aleph_{o}}\right\}\right].$$

**Theorem 2.28.** If the family of variances  $\{\sigma_k^2\}$  satisfies the stronger condition

(2.4) 
$$\sum_{k=1}^{\infty} \frac{\sigma_k^2}{|x_k|} < \infty$$

for some sequence  $(x_k) \in \ell_1$ , then  $\mu_{\mathcal{B}}([\mathcal{B}]) = 1$ .

*Proof.* By definition, if  $f \in \mathcal{B}$  and  $(u_n)$  is a Schauder basis for  $\mathcal{B}$ , then there is a sequence of scalars  $(a_k)$  such that  $f = \lim_{n \to \infty} \sum_{k=1}^n a_k u_k$ . Since  $T(f) = (a_k)$ ,

$$|||(a_n)|||_{\mathcal{B}_I} = \sup_n \left\| \sum_{k=1}^n a_k u_k \right\|_{\mathcal{B}} \leqslant \left[ \sum_{k=1}^\infty |a_k| \right],$$

so that, if  $(a_n) \in \ell_1$ , then  $(a_n) \in T(\mathcal{B}) = \mathcal{B}_I$ .

Suppose that there is a sequence  $(x_k) \in \ell_1$  such that such that the inequality (2.3) is satisfied. As in Lemma 2.24, by Chebyshev's inequality and inequality (2.3) we have

$$\mu_{\mathcal{B}}\left\{T^{-1}\left(\mathop{\times}_{k=1}^{\infty}\left[-\left|x_{k}\right|^{1/2},\left|x_{k}\right|^{1/2}\right]\right)\right\} > 0.$$

If  $A_n = R^n \times (\times_{k=n+1}^{\infty} [-|x_k|^{1/2}, |x_k|^{1/2}])$ , then  $A_n \subset A_{n+1}$  and  $A_n \subseteq \mathcal{B}_I$  for all natural n. Thus, we have

$$\mu_{\mathcal{B}}[T^{-1}(\mathcal{B}_I)] \ge \lim_{n \to \infty} \mu_{\mathcal{B}}[T^{-1}(A_n)]$$

$$= \lim_{n \to \infty} \prod_{k=n+1}^{\infty} \mu_k([-|x_k|^{1/2}, |x_k|^{1/2}]) \ge \lim_{n \to \infty} \prod_{k=n+1}^{\infty} \left(1 - \frac{\sigma_k^2}{|x_k|}\right) = 1.$$

**Definition 2.29.** We call  $\mu_{\mathcal{B}}$  a scaled version of Gaussian measure for  $\mathcal{B}$ .

**Theorem 2.30.** The measure  $\mu_{\mathcal{B}}$  is a countably additive version of Gaussian measure on  $\mathcal{B}$ .

In particular, observe that we obtain a countably additive version of Gaussian measure for both  $\ell_2$  and  $\mathbb{C}_0[0,1]$  (the continuous functions x(t) on [0,1] with x(0)=0).

<span id="page-31-0"></span>2.6. Rotational Invariance. In this section we study rotational invariance on subspaces of  $(\mathbb{R}_I^{\infty}, \mathfrak{B}_I[\mathbb{R}^{\infty}] \lambda_{\infty})$ . First, we need a little more information about Gaussian measures on vector spaces. (See Yamasaki [YA], pg. 151, for a proof of the next Theorem).

Let  $\mathcal{F}$  be a a real vector space, let  $\mathcal{F}^a$  be its algebraic dual space, and let  $\mathfrak{B}_{\mathcal{F}}$  be the smallest  $\sigma$ -algebra such that  $\mathfrak{L}(x)$  is measurable for each functional  $\mathfrak{L} \in \mathcal{F}^a$  and all  $x \in \mathcal{F}$ .

**Theorem 2.31.** If  $\mu$  is a measure on  $(\mathcal{F}^a, \mathfrak{B}_{\mathcal{F}})$ , then the following are equivalent.

(1) The Fourier transform of  $\mu$ ,  $\hat{\mu}$ , is of the form:

$$\hat{\mu}(x) = \exp\left\{-\frac{1}{2}\langle x, x\rangle\right\},\,$$

for some inner product on  $\mathcal{F}$ .

(2) For every  $x \in \mathcal{F}$ , the distribution of  $\mathfrak{L}(x)$  is a one-dimensional Gaussian measure.

In this general setting, a measure  $\mu$  is said to be Gaussian on  $(\mathcal{F}^a, \mathfrak{B}_{\mathcal{F}})$  if it satisfies either of the above conditions.

**Example 2.32.** Let  $\mathcal{F} = \mathbb{R}_0^{\infty}$ , the set of sequences that are zero except for a finite number of terms and let  $\langle \cdot, \cdot \rangle$  be the inner product on  $\mathbb{R}_0^{\infty}$ . It is easy to show that the corresponding measure on  $\mathcal{F}^a = \mathbb{R}^{\infty}$  (satisfying either (1) or (2) above) is the infinite product Gaussian measure.

To understand the importance of this example, let  $(a_n)$  be any sequence of positive numbers and let

(2.5) 
$$\mathcal{H}_a = \left\{ \mathbf{x} \in \mathbb{R}^\infty | \sum_{n=1}^\infty a_n^2 x_n^2 < \infty \right\}.$$

The proof of the following is due to Yamasaki ([YA], pg. 153).

**Lemma 2.33.** If 
$$a \in \ell_2$$
,  $\mu[\mathcal{H}_a] = 1$ , and if  $a \notin \ell_2$ ,  $\mu[\mathcal{H}_a] = 0$ .

Now, let us note that the standard one-dimensional Gaussian density, which is normally written as  $f_X(x) = [\sqrt{2\pi}]^{-1} exp\{-\frac{1}{2}|x|^2\}$ , may also be written as  $f_X(x) = exp\{-\pi|x|^2\}$  with no factors of  $\sqrt{2\pi}$  if we scale  $x \to 0$ 

 $\frac{x}{\sqrt{2\pi}}$ . With this convention, we can write the infinite dimensional version for  $L^2[\mathcal{H}, \lambda_{\mathcal{H}}]$  as the derivative of the Gaussian distribution  $\mu_{\mathcal{H}}$  with respect to the Lebesgue measure on  $\mathcal{H}$ :

(2.6) 
$$f(\mathbf{x}) = \exp\{-\pi |\mathbf{x}|_{\mathcal{H}}^2\} = \frac{d\mu_{\mathcal{H}}(\mathbf{x})}{d\lambda_{\mathcal{H}}(\mathbf{x})}.$$

This shows that, with the appropriate definition of Lebesgue measure, there is a corresponding density for a Gaussian distribution on Hilbert space.

**Remark 2.34.** In the general case (see DePrato [DP]), when  $\mathbb{Q}$  is a (positive definite) trace-class operator and  $\mathbf{x}$  is a Gaussian random variable with mean  $\mathbf{m}$  and covariance  $\mathbb{Q}$ , we can write equation (2.6) as:

$$f(\mathbf{x}) = [\det \mathbb{Q}]^{-1/2} \exp\left\{-\pi \left\langle \mathbb{Q}^{-1}(\mathbf{x} - \mathbf{m}), (\mathbf{x} - \mathbf{m}) \right\rangle_{\mathcal{H}} \right\} \frac{d\mu_{\mathcal{H}}(\mathbf{x})}{d\lambda_{\mathcal{H}}(\mathbf{x})}.$$

**Definition 2.35.** A rotation on  $\mathcal{H}$  is a bijective isometry  $U: \mathcal{H} \to \mathcal{H}$ .

It is well-known that  $\mu_{\mathcal{H}}$  is invariant under rotations over  $(\mathcal{H}, \mathfrak{B}_{\mathcal{H}})$  (see Yamasaki [YA], pg. 163).

**Theorem 2.36.** The measure,  $\lambda_{\mathcal{H}}$ , is invariant under rotations and  $\mathfrak{R} =:$   $(T^{-1}(\ell_2))$  is dense in  $\mathcal{H}$  and the maximal rotation invariance subspace for  $\lambda_{\mathcal{H}}$ .

*Proof.* Let any measurable set  $A \in \mathfrak{B}_{\mathcal{H}}$ . If U is any rotation on  $\mathcal{H}$ , then  $\mu_{\mathcal{H}}(UA) = \mu_{\mathcal{H}}(A)$  and  $|U\mathbf{x}|_{\mathcal{H}}^2 = |\mathbf{x}|_{\mathcal{H}}^2$ . It follows from equation (2.6) that  $\lambda_{\mathcal{H}}(UA) = \lambda_{\mathcal{H}}(A)$ .

It follows from R<sup>∞</sup> <sup>0</sup> ⊂ R ⊂ H, that R is dense, and from Lemma 2.33 that R is maximal.

<span id="page-34-0"></span>Discussion. In this section, we have shown that what appears to be a minor change in the way we represent R<sup>∞</sup> makes it possible to define an analogue of both Lebesgue and Gaussian measure (countably additive) on every (classical) separable Banach space with a Schauder basis. Furthermore, our version of Gaussian measure is rotationally invariant, a property not shared by Wiener measure. (What is more important, we have obtained our core results using basic methods of Lebesgue measure theory from R n .)

### 3. Operators

<span id="page-34-1"></span>This section provides the background to understand the relationship between operators defined on <sup>H</sup><sup>2</sup> <sup>⊗</sup> (which is nonseparable), and their restriction to <sup>H</sup><sup>2</sup> <sup>⊗</sup>(h). We also obtain general conditions that allow us to define infinite sums and products of linear operators on <sup>H</sup><sup>2</sup> <sup>⊗</sup>(h) for a given h.

<span id="page-34-2"></span>3.1. Bounded Operators on <sup>H</sup><sup>2</sup> <sup>⊗</sup>. In this section we review the class of bounded operators on <sup>H</sup><sup>2</sup> <sup>⊗</sup> and their relationship to those on each H<sup>i</sup> . Many of the results are originally due to von Neumann [\[VN2\]](#page-64-4). However, the proofs are new or simplified versions (some from the literature).

Let  $L[\mathcal{H}_{\otimes}^2]$  be the set of bounded operators on  $\mathcal{H}_{\otimes}^2$ . For each fixed  $i_0 \in \mathbb{N}$  and  $A_{i_0} \in L(\mathcal{H}_{i_0})$ , define  $\mathcal{A}_{i_0} \in L(\mathcal{H}_{\otimes}^2)$  by:

$$\mathcal{A}_{i_0}(\sum_{k=1}^N \otimes_{i \in \mathbb{N}} g_i^k) = \sum_{k=1}^N A_{i_0} g_{i_0}^k \otimes (\otimes_{i \neq i_0} g_i^k)$$

for  $\sum_{k=1}^{N} \otimes_{i \in \mathbb{N}} g_i^k$  in  $\mathcal{H}_{\otimes}^2$  and N finite but arbitrary. Extending to all of  $\mathcal{H}_{\otimes}^2$  produces an isometric isomorphism of  $L[\mathcal{H}_{i_0}]$  into  $L[\mathcal{H}_{\otimes}^2]$ , which we denote by  $L[\mathcal{H}(i_0)]$ , so that the relationship  $L[\mathcal{H}_i] \leftrightarrow L[\mathcal{H}(i)]$  is an isometric isomorphism of algebras. Let  $L^{\#}[\mathcal{H}_{\otimes}^2]$  be the uniform closure of the algebra generated by  $\{L[\mathcal{H}(i)], i \in \mathbb{N}\}$ . It is clear that  $L^{\#}[\mathcal{H}_{\otimes}^2] \subset L[\mathcal{H}_{\otimes}^2]$ . von Neumann has shown that the inclusion becomes equality if and only if  $\mathbb{N}$  is replaced by a finite set. On the other hand,  $L^{\#}[\mathcal{H}_{\otimes}^2]$  clearly consists of all operators on  $\mathcal{H}_{\otimes}^2$  that are generated directly from the family  $\{L[\mathcal{H}(i)], i \in \mathbb{N}\}$  by algebraic and topological processes.

Let  $\mathbf{P}_g^s$  denote the projection from  $\mathcal{H}_{\otimes}^2$  onto  $\mathcal{H}_{\otimes}^2(g)^s$ , and let  $\mathbf{P}_g^w$  denote the projection from  $\mathcal{H}_{\otimes}^2$  onto  $\mathcal{H}_{\otimes}^2(g)^w$ .

**Theorem 3.1.** If 
$$\mathbf{T} \in L^{\#}(\mathcal{H}_{\otimes}^{2})$$
, then  $\mathbf{P}_{g}^{s}\mathbf{T} = \mathbf{T}\mathbf{P}_{g}^{s}$  and  $\mathbf{P}_{g}^{w}\mathbf{T} = \mathbf{T}\mathbf{P}_{g}^{w}$ 

Proof. The weak case follows from the strong case, so we prove that  $\mathbf{P}_g^s\mathbf{T} = \mathbf{TP}_g^s$ . Since vectors of the form  $G = \sum_{k=1}^L \otimes_{i \in \mathbb{N}} g_i^k$ , with  $g_i^k = g_i$  for all but a finite number of i, are dense in  $\mathcal{H}^2_{\otimes}(g)^s$ ; it suffices to show that  $\mathbf{T}f \in \mathcal{H}^2_{\otimes}(g)^s$ . Now,  $\mathbf{T} \in L^\#(\mathcal{H}^2_{\otimes})$  implies that there exists a sequence of operators  $\mathbf{T}_n$  such that  $\|\mathbf{T} - \mathbf{T}_n\|_{\otimes} \to 0$  as  $n \to \infty$ , where each  $\mathbf{T}_n$  is of the form:  $\mathbf{T}_n = \sum_{k=1}^{N_n} a_k^n T_k^n$ , with  $a_k^n$  a complex scalar,  $N_n < \infty$ , and each  $T_k^n = \sum_{k=1}^{N_n} a_k^n T_k^n$ , with  $a_k^n$  a complex scalar,  $N_n < \infty$ , and each  $T_k^n = \sum_{k=1}^{N_n} a_k^n T_k^n$ .

 $\hat{\otimes}_{i \in \mathbb{M}_k} T_{ki}^n \hat{\otimes}_{i \in \mathbb{N} \setminus \mathbb{M}_k} I_i$  for some finite set of *i*-values  $\mathbb{M}_k$ , where  $I_i$  is the identity operator on  $\mathcal{H}_i$ . Hence,

$$\mathbf{T}_n f = \sum_{l=1}^L \sum_{k=1}^{N_n} a_k^n \otimes_{i \in \mathbb{M}_k} T_{ki}^n g_i^l \otimes_{i \in \mathbb{N} \setminus \mathbb{M}_k} g_i^l.$$

Now, it is easy to see that, for each l,  $\bigotimes_{i \in \mathbb{M}_k} T_{ki}^n g_i^l \bigotimes_{i \in \mathbb{N} \setminus \mathbb{M}_k} g_i^l \equiv^s \bigotimes_{i \in \mathbb{N}} g_i$ . It follows that  $\mathbf{T}_n f \in \mathcal{H}^2_{\otimes}(g)^s$  for each n, so that  $\mathbf{T}_n \in L[\mathcal{H}^2_{\otimes}(g)^s]$ . Since  $L[\mathcal{H}^2_{\otimes}(g)^s]$  is a norm closed algebra,  $\mathbf{T} \in L[\mathcal{H}^2_{\otimes}(g)^s]$  and it follows that  $\mathbf{P}_q^s \mathbf{T} = \mathbf{T} \mathbf{P}_q^s$ .

Let  $z_i \in \mathbf{C}$ ,  $|z_i| = 1$ , and define  $U[\mathbf{z}]$  by:  $U[\mathbf{z}] \otimes_{i \in \mathbb{N}} g_i = \otimes_{i \in \mathbb{N}} z_i g_i$ .

**Theorem 3.2.** The operator  $U[\mathbf{z}]$  has a unique extension to a unitary operator on  $\mathcal{H}^2_{\otimes}$ , which we also denote by  $U[\mathbf{z}]$ , such that:

- (1)  $U[\mathbf{z}]: \mathcal{H}^2_{\otimes}(g)^w \to \mathcal{H}^2_{\otimes}(g)^w$ , so that  $\mathbf{P}^w_g U[\mathbf{z}] = U[\mathbf{z}]\mathbf{P}^w_g$ .
- (2) If  $\prod_{\nu} z_{\nu}$  is quasi-convergent but not convergent, then  $U[\mathbf{z}]$ :  $\mathcal{H}^{2}_{\otimes}(g)^{s} \to \mathcal{H}^{2}_{\otimes}(h)^{s}$ , for some  $h \in \mathcal{H}^{2}_{\otimes}(g)^{w}$  with  $g \perp h$ .
- (3)  $U[\mathbf{z}]: \mathcal{H}^2_{\otimes}(g)^s \to \mathcal{H}^2_{\otimes}(g)^s$  if and only if  $\prod_i z_i$  converges and  $U[\mathbf{z}] = (\prod_i z_i)\mathbf{I}_{\otimes}$ , where  $\mathbf{I}_{\otimes}$  is the identity operator on  $\mathcal{H}^2_{\otimes}$ . This implies that  $\mathbf{P}^s_gU[\mathbf{z}] = U[\mathbf{z}]\mathbf{P}^s_g$ .

*Proof.* For (1), let  $h = \sum_{k=1}^{N} \otimes_{i \in \mathbb{N}} h_i^k$ , where  $\otimes_{i \in \mathbb{N}} h_i^k \equiv^w \otimes_{i \in \mathbb{N}} g_i$ , N is arbitrary and  $1 \leq k \leq N$ . Then

$$U^*[\mathbf{z}]U[\mathbf{z}]h = \sum_{k=1}^N \otimes_{i \in \mathbb{N}} z_i^* z_i h_i^k = h = U[\mathbf{z}]U^*[\mathbf{z}]h.$$

Thus, we see that  $U[\mathbf{z}]$  is a unitary operator, and since h of the above form are dense,  $U[\mathbf{z}]$  extends to a unitary operator on  $\mathcal{H}^2_{\otimes}$ . By definition,  $\sum_{k=1}^N \otimes_{i \in \mathbb{N}} z_i h_i^k \in \mathcal{H}^2_{\otimes}(g)^w$  if  $\sum_{k=1}^N \otimes_{i \in \mathbb{N}} h_i^k \in \mathcal{H}^2_{\otimes}(g)^w$ , so that  $U[\mathbf{z}]: \mathcal{H}^2_{\otimes}(g)^w \to \mathcal{H}^2_{\otimes}(g)^w$  and  $\mathbf{P}^w_g U[\mathbf{z}] = U[\mathbf{z}] \mathbf{P}^w_g$ . To prove (2), use Theorem 1.6 (3) and (4) to note that  $\prod_i z_i = 0$  and  $\bigotimes_{i \in \mathbb{N}} h_i^k \equiv^s \bigotimes_{i \in \mathbb{N}} g_i$  imply that  $\bigotimes_{i \in \mathbb{N}} z_i h_i^k \in \mathcal{H}^2_{\otimes}(f)^s$  with  $\mathcal{H}^2_{\otimes}(f)^s \perp \mathcal{H}^2_{\otimes}(g)^s$ . To prove (3), note that, if  $0 < |\prod_i z_i| < \infty$ , then  $U[\mathbf{z}] = [(\prod_i z_i) \mathbf{I}_{\otimes}]$ , so that  $U[\mathbf{z}]: \mathcal{H}^2_{\otimes}(g)^s \to \mathcal{H}^2_{\otimes}(g)^s$ . Now suppose that  $U[\mathbf{z}]: \mathcal{H}^2_{\otimes}(g)^s \to \mathcal{H}^2_{\otimes}(g)^s$ , then  $\bigotimes_{i \in \mathbb{N}} z_i h_i^k \equiv^s \bigotimes_{i \in \mathbb{N}} g_i$  and so  $\prod_i z_i$  must converge. Therefore,  $U[\mathbf{z}]h = [(\prod_i z_i) \mathbf{I}_{\otimes}]h$  and  $\mathbf{P}^s_g U[\mathbf{z}] = U[\mathbf{z}] \mathbf{P}^s_g$ .

It is easy to see that, for each fixed  $i \in \mathbb{N}$ ,  $\mathcal{A}(i) \in L[\mathcal{H}(i)]$  commutes with any  $\mathbf{P}_g^s$ ,  $\mathbf{P}_g^w$  or  $U[\mathbf{z}]$ , where g and  $\mathbf{z}$  are arbitrary.

**Theorem 3.3.** Every  $\mathbf{T} \in L^{\#}[\mathcal{H}_{\otimes}^2]$  commutes with all  $\mathbf{P}_g^s$ ,  $\mathbf{P}_g^w$  and  $U[\mathbf{z}]$ , where g and  $\mathbf{z}$  are arbitrary.

*Proof.* Let  $\mathfrak L$  be the set of all  $\mathbf P_g^s$ ,  $\mathbf P_g^w$  or  $U[\mathbf z]$ , with g and  $\mathbf z$  arbitrary. From the above observation, we see that all  $\mathcal A_i \in L[\mathcal H(i)]$ ,  $i \in \mathbb N$ , commute with  $\mathfrak L$  and hence belong to its commutator  $\mathfrak L'$ . Since  $\mathfrak L'$  is a closed algebra, this implies that  $L^\#[\mathcal H_\otimes^2] \subseteq \mathfrak L'$  so that all  $\mathbf T \in L^\#[\mathcal H_\otimes^2]$  commute with  $\mathfrak L$ .

<span id="page-37-0"></span>3.2. Unbounded Operators on  $\mathcal{H}^2_{\otimes}$ . In this section, we consider a restricted class of unbounded operators and the notion of a strong convergence vector introduced by Reed [RE].

For each  $i \in \mathbb{N}$ , let  $A_i$  be a closed densely defined linear operator on  $\mathcal{H}_i$ , with domain  $D(A_i)$ , and let  $A_i$  be its extension to  $\mathcal{H}^2_{\otimes}$ , with domain  $D(A_i) \supset \tilde{D}(A_i) = D(A_i) \otimes (\otimes_{k \neq i} \mathcal{H}_k)$ . The next theorem follows directly from the definition of the tensor product of semigroups.

**Theorem 3.4.** Let  $A_i$ ,  $1 \leq i \leq n$ , be generators of a family of  $C_0$ semigroups  $S_i(t)$  on  $\mathcal{H}_i$  with  $||S_i(t)||_{\mathcal{H}_i} \leq M_i e^{\omega_i t}$ . Then  $\mathbf{S}_n(t) = \hat{\otimes}_{i=1,n} S_i(t)$ ,
defined on  $\hat{\otimes}_{i=1,n} \mathcal{H}_i$ , has a unique extension (also denoted by  $\mathbf{S}_n(t)$ ) to all
of  $\mathcal{H}^2_{\otimes}$ , such that, for all vectors  $\sum_{k=1}^K \otimes_{i \in \mathbb{N}} g_i^k$  with  $g_l^k \in D(A_l)$ ,  $1 \leq l \leq n$ ,
the infinitesimal generator for  $\mathbf{S}_n(t)$  satisfies:

$$\mathcal{A}^n \left[ \sum_{k=1}^K \otimes_{i \in \mathbb{N}} g_i^k \right] = \sum_{l=1}^n \sum_{k=1}^K A_l g_l^k (\otimes_{i \in \mathbb{N}}^{i \neq l} g_i^k).$$

**Definition 3.5.** Let  $\{A_i\}$ ,  $i \in \mathbb{N}$ , be a family of closed densely defined linear operators on  $\mathcal{H}_i$  and let  $g_i \in D(A_i)$  (respectively,  $f_i \in D(A_i)$ ), with  $\|g_i\|_{\mathcal{H}} = 1$  (respectively,  $\|f_i\|_{\mathcal{H}} = 1$ ), for all  $i \in \mathbb{N}$ .

- (1) We say that  $g = \bigotimes_{i \in \mathbb{N}} g_i$  is a strong convergence sum (scs)-vector for the family  $\{A_i\}$  if  $s = \lim_{n \to \infty} \sum_{k=1}^n A_k g = \sum_{k=1}^\infty A_k g_k (\bigotimes_{i \in \mathbb{N}}^{i \neq k} g_i)$  exists.
- (2) We say that  $f = \bigotimes_{i \in \mathbb{N}} f_i$  is a strong convergence product (scp)-vector for the family  $\{A_i\}$  if  $s \lim_{n \to \infty} \prod_{k=1}^n A_k f = \bigotimes_{i \in \mathbb{N}} A_i f_i$  exists.

Let  $\mathbf{D}_g$  be the linear span of  $\{\chi = \bigotimes_{i \in \mathbb{N}} \chi_i, \ \chi_i \in D(A_i)\}$ , with  $\chi_i = g_i$  (and let  $\mathbf{D}_f$  be the linear span of  $\{\eta = \bigotimes_{i \in \mathbb{N}} \eta_i, \ \eta_i \in D(A_i)\}$ , with  $\eta_i = f_i$ ) for all i > L, where L is arbitrary but finite. Clearly,  $\mathbf{D}_g$  is dense in  $\mathcal{H}^2_{\otimes}(g)^s$  ( $\mathbf{D}_{\eta}$  is dense in  $\mathcal{H}^2_{\otimes}(f)^s$ ). If there is a possible chance for confusion, we

let  $\mathcal{A}_s$ , respectively  $\mathcal{A}_p$ , denote the closure of  $\sum_{k=1}^{\infty} \mathcal{A}_k$  on  $\mathcal{H}^2_{\otimes}(g)^s$  (respectively  $\prod_{k=1}^{\infty} \mathcal{A}_k$  on  $\mathcal{H}^2_{\otimes}(f)^s$ ). It follows that  $\mathcal{H}^2_{\otimes}(g)^s$  (respectively  $\mathcal{H}^2_{\otimes}(f)^s$ ) are natural spaces for the study of infinite sums or products of unbounded operators. (The notion of a strong convergence sum vector first appeared in Reed [RE].)

**Definition 3.6.** We call  $\mathcal{H}^2_{\otimes}(g)^s$  an  $\mathbb{RS}$ -space (respectively,  $\mathcal{H}^2_{\otimes}(f)^s$  an  $\mathbb{RP}$ -space ) for the family  $\{A_i\}$ .

Let  $\{U_k(t)\}$  be a set of unitary groups on  $\{\mathcal{H}_k\}$ . It is easy to see that  $U(t) = \hat{\otimes}_{k=1}^{\infty} U_k(t)$  is a unitary group on  $\mathcal{H}_{\otimes}^2$ . However, we know from Theorem 3.2 (2), that it need not be reduced on any partial tensor product subspace. The following results are due to Streit [ST] and Reed [RE], as indicated.

**Theorem 3.7.** (Streit) Suppose  $\{A_k\}$  is a set of selfadjoint linear operators on the space  $\mathcal{H}^2_{\otimes}(g)^s$ , with corresponding unitary groups  $\{U_k(t)\}$ . If  $U(t) = \hat{\otimes}_{k=1}^{\infty} U_k(t)$ , then  $\mathbf{P}^s_g U(t) = U(t) \mathbf{P}^s_g$  (i.e., U(t) is reduced on  $\mathcal{H}^2_{\otimes}(g)^s$ ) and U(t) is a strongly continuous unitary group on  $\mathcal{H}^2_{\otimes}(g)^s$  if and only if, for each c > 0, the following three conditions are satisfied:

$$(1) \sum_{k=1}^{\infty} |\langle \mathcal{A}_k E_k[-c, c|g_k, g_k \rangle| < \infty,$$

(2) 
$$\sum_{k=1}^{\infty} \left| \left\langle \mathcal{A}_k^2 E_k[-c,c] g_k, g_k \right\rangle \right|$$
,

(3) 
$$\sum_{k=1}^{\infty} |\langle (I_k - E_k[-c, c]g_k, g_k)| < \infty,$$

where  $E_k[-c,c]$  are the spectral projectors of  $\mathcal{A}_k$  and, in this case,  $U(t) = s - \lim_{n \to \infty} \hat{\otimes}_{k=1}^n U_k(t)$ .

Corollary 3.8. Conditions 1-3 are satisfied if and only if there exists a strong convergence vector  $g = \bigotimes_{k=1}^{\infty} g_k$  for the family  $\{A_k\}$  such that  $g_k \in D(A_k)$  and

$$\sum_{k=1}^{\infty} |\langle \mathcal{A}_k g_k, g_k \rangle| < \infty, \ \sum_{k=1}^{\infty} ||\mathcal{A}_k g_k||^2 < \infty.$$

**Theorem 3.9.** (Reed) U(t) is reduced on  $\mathcal{H}^2_{\otimes}(g)^s$  and U(t) is a strongly continuous unitary group on  $\mathcal{H}^2_{\otimes}(g)^s$  if and only if  $g = \bigotimes_{k=1}^{\infty} g$  is a strong convergence vector for the family  $\{A_k\}$  and  $\sum_{k=1}^{\infty} |\langle \mathcal{A}_k g_k, g_k \rangle| < \infty$ . If each  $A_k$  is positive, the statement is true without the above condition. In either case,  $\mathcal{A}$ , the closure of  $\sum_{k=1}^{\infty} \mathcal{A}_k$ , is the generator of U(t).

The next result strengthens and extends Reed's theorem to contraction semigroups (i.e., the positivity requirement above can be dropped).

**Theorem 3.10.** Let  $\{S_k(t)\}$  be a family of strongly continuous contraction semigroups with generators  $\{A_k\}$  defined on  $\{\mathcal{H}_k\}$ , and let  $g = \bigotimes_{k=1}^{\infty} g_k$  be a strong convergence vector for the family  $\{A_k\}$ . Then  $\mathbf{S}(t) = \hat{\otimes}_{k=1}^{\infty} S_k(t)$  is reduced on  $\mathcal{H}^2_{\otimes}(g)^s$  and is a strongly continuous contraction semigroup. If  $\mathbf{S}(t) = \hat{\otimes}_{k=1}^{\infty} S_k(t)$  is reduced on  $\mathcal{H}^2_{\otimes}(g)^s$  and is a strongly continuous contraction semigroup on  $\mathcal{H}^2_{\otimes}(g)^s$ , then there exists a strong convergence vector  $f = \bigotimes_{k=1}^{\infty} f_k \in \mathcal{H}^2_{\otimes}(g)^s$  for the family  $\{A_k\}$ .

Proof. Let  $g = \bigotimes_{k=1}^{\infty} g_k$  be a strong convergence vector for the family  $\{A_k\}$ . Without loss, we can assume that  $\|g_k\| = 1$ . Let  $\mathbf{S}_n(t) = \hat{\otimes}_{k=1}^n S_k(t) \hat{\otimes} (\bigotimes_{k=n+1}^{\infty} I_k)$  and observe that  $\mathbf{S}_n(t)$  is a contraction semigroup on  $\mathcal{H}^2_{\otimes}(g)^s$  for all finite n. Furthermore, its generator is the closure of  $\mathcal{A}^n = \sum_{k=1}^n \mathcal{A}_k$ , where  $\mathcal{A}_k = A_k \hat{\otimes} (\bigotimes_{i \neq k}^{\infty} I_i)$ . If n and m are arbitrary, then

$$[\mathbf{S}_n(t) - \mathbf{S}_m(t)] g = \int_0^1 \frac{d}{d\lambda} \{ \mathbf{S}_n[\lambda t] \mathbf{S}_m[(1-\lambda)t] \} g d\lambda$$
$$= t \int_0^1 \mathbf{S}_n[\lambda t] \mathbf{S}_m[(1-\lambda)t] [\mathcal{A}^n - \mathcal{A}^m] g d\lambda,$$

where we have used the fact that, if two semigroups commute, then their corresponding generators also commute. It follows that:

$$\|[\mathbf{S}_n(t) - \mathbf{S}_m(t)]g\| \leqslant t \|[\mathcal{A}^n - \mathcal{A}^m]g\|.$$

Since  $g = \bigotimes_{k=1}^{\infty} g$  is a strong convergence vector for the family  $\{A_k\}$ , it follows that

 $s-\lim_{n\to\infty} \mathbf{S}_n(t) = \mathbf{S}(t)$  exists on a dense set in  $\mathcal{H}^2_{\otimes}(g)^s$  and the convergence is uniform on bounded t intervals. It follows that S(t) extends to a bounded linear operator on  $\mathcal{H}^2_{\otimes}(g)^s$ . To see that the closure of S(t) must be a contraction, for any  $\varepsilon > 0$ , choose n so large that  $\|[\mathbf{S}_n(t) - \mathbf{S}(t)]g\|_{\otimes} < \varepsilon \|g\|_{\otimes}$ . It follows that

$$\|\mathbf{S}(t)g\|_{\otimes} \leq \|\mathbf{S}_n(t)g\|_{\otimes} + \|[\mathbf{S}_n(t) - \mathbf{S}(t)]g\|_{\otimes} < \|g\|_{\otimes} (1+\varepsilon).$$

Thus,  $\mathbf{S}(t)$  is a contraction operator on  $\mathcal{H}^2_{\otimes}(g)^s$ . It is easy to check that it is a  $C_0$ -semigroup.

Now suppose that  $\mathbf{S}(t) = \hat{\otimes}_{k=1}^{\infty} S_k(t)$  is a strongly continuous contraction semigroup which is reduced on  $\mathcal{H}^2_{\otimes}(g)^s$ . It follows that the generator  $\mathcal{A}$  of  $\mathbf{S}(t)$  is m-dissipative, and hence defined on a dense domain  $D(\mathcal{A})$  in  $\mathcal{H}^2_{\otimes}(g)^s$  with  $\mathbf{S}'(t)f = \mathbf{S}(t)\mathcal{A}f = \mathcal{A}\mathbf{S}(t)f$  for all  $f \in D(\mathcal{A})$ . Since any such f is of the form  $f = \sum_{l=1}^{\infty} \otimes_{k=1}^{\infty} f_k^l$ , each  $f^l = \otimes_{k=1}^{\infty} f_k^l$  is in  $D(\mathcal{A})$ . A simple computation shows that  $\mathcal{A}f^l = \sum_{k=1}^{\infty} \mathcal{A}_k f^l$ , so that any  $f^l$  is a strong convergence vector for the family  $\{A_k\}$ .

It is easy to see that, in the second part of the theorem, we cannot require that  $g = \bigotimes_{k=1}^{\infty} g_k$  itself be a strong convergence vector for the family  $\{A_k\}$  since it need not be in the domain of  $\mathcal{A}$ . For example,  $g_1 \notin D(A_1)$ , while  $g_k \in D(A_k), k \neq 1$ .

#### 4. Function Spaces

<span id="page-42-0"></span>Let  $\chi_{I_n}$  be the indicator (or characteristic) function of  $I_n = \times_{k=n+1}^{\infty} I$ . If we let  $\mathfrak{L}(\mathbb{R}^n)$  represent the class of measurable functions on  $\mathbb{R}^n$ , then for each measurable function  $f_n \in \mathfrak{L}(\mathbb{R}^n)$  we identify  $f \in \mathfrak{L}(\mathbb{R}^n)$  by  $f = f_n \otimes \chi_{I_n}$ .

**Definition 4.1.** A real-valued function f defined on the measure space  $(\mathbb{R}_I^{\infty}, \mathfrak{B}[\mathbb{R}_I^{\infty}], \lambda_{\infty})$  is said to be measurable if  $f^{-1}(A) \in \mathfrak{B}[\mathbb{R}_I^{\infty}]$  for every  $A \in \mathfrak{B}[\mathbb{R}]$ .

In this section we develop those aspects of function space theory that will be of use later. We note that all the standard theorems for Lebesgue measure apply. (The proofs are the same as for integration on  $\mathbb{R}^n$ .)

<span id="page-43-0"></span>4.1. L 1 -Theory. Let L 1 [R n I ] be the class of integrable functions on R n I . Since L 1 (R n I ) ⊂ L 1 (R n+1 I ), we define L 1 [R ′∞ I ] = S∞ <sup>n</sup>=1 L 1 (R n I ) and let L 1 [R<sup>∞</sup> I be the norm closure of L 1 [R ′∞ I ]. It follows that every function in L 1 [R<sup>∞</sup> I ] is the limit of a sequence of functions in L 1 [R nk I ], for some sequence {nk} ⊂ N.

Let Cc(R n I ) be the class of continuous functions on R n <sup>I</sup> which vanish outside compact sets. We define Cc(R<sup>∞</sup> I ) to be the closure of S∞ <sup>n</sup>=1 Cc(R n I ) = Cc(R ′∞ I ) in the sup norm. Thus, for any f ∈ Cc(R<sup>∞</sup> I ), there always exists a sequence of functions {fn<sup>k</sup> } ∈ Cc(R nk I ) such that fn<sup>k</sup> → f, for some sequence {nk} ⊂ N. We define C0(R<sup>∞</sup> I ), the functions that vanish at ∞, in the same manner.

Lemma 4.2. If f ∈ Cc(R<sup>∞</sup> I ) or C0(R<sup>∞</sup> I ), then f is continuous.

Proof. Let f(x) ∈ Cc(R<sup>∞</sup> I ) and let {x<sup>n</sup> | n = 1, 2, . . . } be any sequence in R n I such that x<sup>n</sup> → x as n → ∞. If ε > 0 is given, choose K<sup>1</sup> so that for k ≥ K<sup>1</sup> and f<sup>k</sup> ∈ Cc(R<sup>∞</sup> I ), |fk(xn) − f(xn)| < ε 3 . Then choose K<sup>2</sup> so that for k ≥ <sup>K</sup>2, <sup>|</sup>fk(x) <sup>−</sup> <sup>f</sup>(x)<sup>|</sup> <sup>&</sup>lt; <sup>ε</sup> 3 . Choose <sup>N</sup> so that for <sup>n</sup> <sup>≥</sup> N, <sup>|</sup>fk(xn) <sup>−</sup> <sup>f</sup>k(x)<sup>|</sup> <sup>&</sup>lt; <sup>ε</sup> 3 . If n ≥ N and k ≥ max{K1, K2}, we have:

$$|f(\mathbf{x}_n) - f(\mathbf{x})| \le |f_k(\mathbf{x}_n) - f(\mathbf{x}_n)| + |f_k(\mathbf{x}) - f_k(\mathbf{x}_n)| + |f_k(\mathbf{x}) - f(\mathbf{x})| < \varepsilon.$$

The same proof applies to C0(R<sup>∞</sup> I )

Theorem 4.3. Cc(R<sup>∞</sup> I ) is dense in L 1 (R<sup>∞</sup> I ).

Proof. We prove this result in the standard manner, by reducing the proof to positive simple functions and then to one characteristic function and finally using the approximation theorem to approximate a measurable set which contains a closed set and is contained in an open set.

First note that, since  $\lim_{k\to\infty} \|f\chi_{B_I(0,k)} - f\|_1 = 0$  for all  $f\in L^1$  (by the DCT), we can prove the result for functions that vanish outside a compact set. In this case, as  $f = f_+ - f_-$ , we need only consider positive f. However, this function can be approximated by simple functions in  $S_+$ . Since each simple function is a finite sum of characteristic functions (of bounded measurable sets) multiplied by finite constants, it follows that we need only show that we can approximate the characteristic function of a bounded measurable set by a continuous function which vanishes outside a compact set. Let  $\varepsilon > 0$  be given and let  $g = \chi_A$ , where A is any bounded measurable set. By the regularity of  $\lambda_\infty$ , there exists an open set O and a compact set H with  $H \subset A \subset O$  and  $\lambda_\infty(O \setminus H) < \varepsilon$ .

Let  $\{V_n\}$  be the class of open intervals with rational end points. For each  $n \in \mathbb{N}$ , let  $F_n \subset g^{-1}[V_n]$  and  $G_n \subset (O \setminus g^{-1}[V_n])$  be compact sets, such that  $\lambda_{\infty}[(O \setminus F_n \cup G_n)] < \frac{\varepsilon}{2^n}$ . If  $H = \bigcap_{n=1}^{\infty} [F_n \cup G_n]$ , then  $\lambda_{\infty}(O \setminus H) < \varepsilon$ .

If  $x \in H$ , there is an n such that  $f(x) \in V_n$  and  $x \in G_n^c$ , so that  $g[G_n^c \cap H] \subset V_n$ . It follows that g restricted to H is continuous and  $\lambda_{\infty}(A \setminus H) \leq \lambda_{\infty}(O \setminus H) < \varepsilon$ .

In a similar fashion we can define the  $L^p$  spaces, 1 . We should note that, each space is defined relative to the family of indicator functions

<span id="page-45-0"></span>for I. Thus, each space is the canonical one for that particular class of spaces.

## 5. Fourier Transform Theory

In this section, we study the implications of Lebesgue measure on R<sup>∞</sup> for the Fourier transform and discuss two different extensions of the Pontrjagin Duality theory for Banach spaces.

<span id="page-45-1"></span>Background. Let G be a locally compact abelian (LCA) group (c.f., R n ). The following is a restatement of Theorem 0.1 (see Rudin [\[RU1\]](#page-63-8)).

Theorem 5.1. If G is a LCA group and B(G) is the Borel σ-algebra of subsets of G, then there is a non-negative regular translation invariant measure µ (i.e., µ(g + A) = µ(A), A ∈ B(G). The (Haar) measure µ is unique up to multiplication by a constant.

Definition 5.2. A complex valued function α : G → C on a LCA group is called a character on G provided that α is a homomorphism and |α(g)| = 1 for all g ∈ G.

The set of all continuous characters of G defines a new group Gˆ, called the dual group of G and (α<sup>1</sup> + α2)(g) = α1(g) · α2(g). If we define a map γ : G → <sup>ˆ</sup>Gˆ, by <sup>γ</sup>g(α) = <sup>α</sup>(g), then the following theorem was proven by Pontryagin:

Theorem 5.3. (Pontryagin Duality Theorem) If G is a LCA-group, then the mapping γ : G → <sup>ˆ</sup>G<sup>ˆ</sup> is an isomorphism of topological groups.

Thus, Pontrjagin Duality identifies those groups that are the character groups of their character groups. If the group is not locally compact Theorem 5.1 does not hold (e.g., there is no Haar measure). However Kaplan [\[KA1\]](#page-62-10) has shown that the class of topological abelian groups for which the Pontrjagin Duality holds is closed under the operation of taking infinite products of groups. This result immediately implies that this class is larger than the class of locally compact abelian groups because the infinite product of locally compact groups (for example, R∞) may be non-locally compact (see also [\[KA2\]](#page-62-11)).

<span id="page-46-0"></span>5.1. Pontryagin Duality Theory I. In this section, we treat the Fourier transform as an operator. As will be seen, this approach has the advantage of being constructive. It also provides us with some insight into the problem that arises when we look at analysis on infinite dimensional spaces.

We define F on L 1 [R, λ] by

$$\hat{g}(x) = \mathfrak{F}(g)(x) = \int_{\mathbb{R}} \exp\{-2\pi i xy\} g(y) dy.$$

It is easy to check that F −1 is defined by

$$g(y) = \mathfrak{F}^{-1}(\hat{g})(y) = \int_{\mathbb{R}} \exp\{2\pi i y x\} \hat{g}(x) dx.$$

This representation is more convenient for the infinite-dimensional case, because we have no factors of  $\sqrt{2\pi}$  to worry about.

It is possible to define  $\mathfrak{F}$  as a mapping on  $L^1[\mathbb{R}^n_I, \lambda]$  to  $\mathbb{C}_0[\mathbb{R}^n_I, \lambda]$  for all n as one fixed linear operator. However, in the case of Hilbert spaces, Theorem 3.2(2) requires that we clearly specify our canonical domain and range space. The same is also true for  $L^1[\mathbb{R}^n_I, \lambda](h)$  and  $\mathbb{C}_0[\mathbb{R}^n_I, \lambda](\hat{h})$  (see [GZ]). Since  $h = \bigotimes_{k=1}^{\infty} \chi_I(x_k)$ , an easy calculation shows that  $\hat{h} = \bigotimes_{k=1}^{\infty} \frac{\sin(\pi x_k)}{\pi x_k}$ . Thus, we can define  $\mathfrak{F}(f_n)(\mathbf{x})$ , mapping  $\mathbf{L}^1[\mathbb{R}^n_I](h)$  into  $\mathbb{C}_0[\mathbb{R}^n_I](\hat{h})$  by

(5.1) 
$$\mathfrak{F}(f_n)(\mathbf{x}) = \bigotimes_{k=1}^n \mathfrak{F}_k(f_{(n)}) \bigotimes_{k=n+1}^\infty \hat{h}_k(x_k).$$

**Theorem 5.4.** The operator  $\mathfrak{F}$  extends to a bounded linear mapping of  $L^1[\mathbb{R}_I^{\infty}](h)$  into  $\mathbb{C}_0[\mathbb{R}_I^{\infty}](\hat{h})$ .

*Proof.* Since

$$\lim_{n \to \infty} L^1[\mathbb{R}_I^n](h) = \bigcup_{n=1}^{\infty} L^1[\mathbb{R}_I^n](h) = L^1[\mathbb{R}_I'^\infty](h)$$

and  $L^1[\mathbb{R}_I^{\infty}](h)$  is the closure of  $L^1[\mathbb{R}_I^{\infty}](h)$  in the  $L^1$  - norm, it follows that  $\mathfrak{F}$  is a bounded linear mapping of  $L^1[\mathbb{R}_I^{\infty}](h)$  into  $\mathbb{C}_0[\mathbb{R}_I^{\infty}](\hat{h})$ .

Supposed that  $\{f_n\} \subset L^1[\mathbb{R}_I^{\infty}](h)$ , converges to  $f \in L^1[\mathbb{R}_I^{\infty}](h)$ . Thus, the sequence is Cauchy, so that  $||f_n - f_m||_1 \to 0$  as  $m, n \to \infty$ . It follows that

$$|\mathfrak{F}(f_n(\mathbf{x}) - f_m(\mathbf{x}))| \leqslant \int_{\mathbb{R}_I^{\infty}} |f_n(\mathbf{y}) - f_m(\mathbf{y})| d\lambda_{\infty}(\mathbf{y}) = ||f_n - f_m||_1,$$

so that  $|\mathfrak{F}(f_n(\mathbf{x}) - f_m(\mathbf{x}))|$  is a Cauchy sequence in  $\mathbb{C}_0[\mathbb{R}_I^{\infty}](\hat{h})$ . Since  $L^1[\mathbb{R}_I'^{\infty}](h)$  is dense in  $L^1[\mathbb{R}_I^{\infty}](h)$ , it follows that  $\mathfrak{F}$  has a bounded extension, mapping  $L^1[\mathbb{R}_I^{\infty}](h)$  into  $\mathbb{C}_0[\mathbb{R}_I^{\infty}](\hat{h})$ .

<span id="page-48-0"></span>5.2.  $L^2$ -Theory. In the case of  $L^2$ , the Fourier transform is an isometric isomorphism from  $L^2[\mathbb{R}^n]$  onto  $L^2[\mathbb{R}^n]$ .

**Theorem 5.5.** The operator  $\mathfrak{F}$  is an isometric isomorphism of  $L^2[\mathbb{R}_I^{\infty}](h)$  onto  $L^2[\mathbb{R}_I^{\infty}](\hat{h})$ .

Proof. Let  $f \in L^2[\mathbb{R}_I^\infty](h)$ . By construction, there exists a sequence of functions  $\{f_k \in L^2[\mathbb{R}_I^{n_k}], n_k \in \mathbb{N}\}$  such that  $\lim_{k \to \infty} \|f - f_k\|_2(h) = 0$ . Furthermore, since the sequence converges, it is a Cauchy sequence. Hence, given  $\varepsilon > 0$ , there exists a  $N(\varepsilon)$  such that  $m, n \geq N(\varepsilon)$  implies that  $\|f_m - f_n\|_2(h) < \varepsilon$ . Since  $\mathfrak{F}$  is an isometry,  $\|\mathfrak{F}(f_m) - \mathfrak{F}(f_n)\|_2(\hat{h}) < \varepsilon$ , so that the sequence  $\mathfrak{F}(f_k)$  is also a Cauchy sequence in  $L^2[\mathbb{R}_I^\infty](\hat{h})$ . Thus, there is a  $\hat{f} \in L^2[\mathbb{R}_I^\infty](\hat{h})$  with  $\lim_{k \to \infty} \|\hat{f} - \mathfrak{F}(f_k)\|_2(\hat{h}) = 0$ , and we can define  $\mathfrak{F}(f) = \hat{f}$ . It is easy to see that  $\hat{f}$  is unique.

We can also prove a version of Theorems 5.4 and 5.5 for every separable Banach space (with a basis). Fix  $\mathcal{B}$  and for each n, let  $\mathcal{B}_I^n = \mathcal{B} \cap \mathbb{R}_I^n$ . It is clear that  $\mathcal{B}_I^n \subset \mathcal{B}_I^{n+1}$ , so that  $\mathcal{B}$  is the norm closure of  $\lim_{n \to \infty} \mathcal{B}_I^n$ . The following have the same proofs as Theorems 5.4 and 5.5.

**Theorem 5.6.** The operator  $\mathfrak{F}$  extends to a bounded linear mapping of  $L^1[\mathcal{B}](h)$  into  $\mathbb{C}_0[\mathcal{B}](\hat{h})$ .

**Theorem 5.7.** The operator  $\mathfrak{F}$  is an isometric isomorphism from  $L^2[\mathcal{B}](h)$  onto  $L^2[\mathcal{B}](\hat{h})$ .

Theorems 5.4 - 5.7 show that  $\bigotimes_{i=1}^{\infty} \hat{h}_i$  is a strong (product) convergence vector for the Fourier transform operator  $\mathfrak{F}$ . In the  $L^2$ -theory, we know that  $L^2[\mathbb{R}_I^{\infty}](h)$  and  $L^2[\mathbb{R}_I^{\infty}](\hat{h})$  are orthogonal subspaces of  $\mathcal{H}_{\otimes}^2$ . Thus, in this approach, the natural interpretation is that the Fourier transform induces a Pontryagin duality like theory that does not depend on the group structure of  $\mathbb{R}_I^{\infty}$ , but depends on the pairing of different function spaces. This approach is direct, constructive and applies to all separable Banach spaces (with a basis). Thus, the group structure of the underlying measure space plays no role.

<span id="page-49-0"></span>5.3. Pontryagin Duality Theory II. In this section, we show that the standard form of Pontryagin duality theory is also possible, using the underlying measure space group structure. It is constructive but restrictive, in that, it does not apply to every separable Banach space with a basis.

Let  $\mathcal{B}$  be any uniformly convex separable Banach space UCB over the reals, so that  $\mathcal{B} = \mathcal{B}^{**}$  (second dual). The next theorem follows from our theory of Lebesgue measure on Banach spaces.

**Theorem 5.8.** If  $\lambda_{\mathcal{B}}$  is our version of Lebesgue measure on  $\mathcal{B}$ , then  $\mathcal{B}$  and  $\mathcal{B}^*$  are also duals as character groups (i.e.,  $\mathcal{B}^* = \hat{\mathcal{B}}$ ).

*Proof.* If we consider the restriction to  $L^2[\mathcal{B}, \lambda_{\mathcal{B}}]$ , we can define  $\mathfrak{F}$  directly by:

(5.2) 
$$[\mathfrak{F}(f)](\mathbf{x}^*) = \hat{f}(\mathbf{x}^*) = \int_{\mathcal{B}} \exp\{-2\pi i \langle \mathbf{y}, \mathbf{x}^* \rangle\} f(\mathbf{y}) d\lambda_{\mathcal{B}}(\mathbf{y}),$$

where  $\langle \mathbf{y}, \mathbf{x}^* \rangle$  is the pairing between  $\mathcal{B}$  and  $\mathcal{B}^*$ . From Plancherel's Theorem, we have:

$$\|\hat{f}\|_{2}^{2} = (\hat{f}, \hat{f})_{2} = (f, f)_{2} = \|f\|_{2}^{2}.$$

It follows that  $\mathcal{B}$  and  $\mathcal{B}^*$  are duals as character groups and

$$f(\mathbf{y}) = \int_{\mathcal{B}^*} \exp\{2\pi i \langle y, \mathbf{x}^* \rangle\} \hat{f}(\mathbf{x}^*) d\lambda_{\mathcal{B}^*}(\mathbf{x}^*).$$

If  $\mathcal{B}_I^n = \mathcal{B}_I \cap \mathbb{R}_I^n$ , we can represent  $\hat{f}_n$  directly as a mapping from  $L^2[\mathcal{B}_I^n, \lambda_{\mathcal{B}}] \to L^2[\mathcal{B}_I^{*,n}, \lambda_{\mathcal{B}^*}]$ , by

$$[\mathfrak{F}(f_n)](\mathbf{x}^*) = \hat{f}_n(\mathbf{x}^*) = \int_{\mathcal{B}} \exp\{-2\pi i \langle \mathbf{y}, \mathbf{x}^* \rangle_n\} f_n(\mathbf{y}) d\lambda_{\mathcal{B}}(\mathbf{y}),$$

where  $\langle \mathbf{y}, \mathbf{x}^* \rangle_n$  is the restricted pairing of  $\mathbf{y}$  and  $\mathbf{x}^*$  to  $\mathcal{B}_I^n$  and  $\mathcal{B}_I^{*,n}$  respectively. It follows that equations 5.1 and 5.2 provide two distinct definitions of the Fourier transform for  $\mathcal{B}$ . Thus, in this approach the group structure of the underlying measure space changes.

It is clear that representation for  $\hat{f}(\mathbf{x}^*)$  also applies if we use  $L^1[\mathcal{B}, \lambda_{\mathcal{B}}]$ , but in this case  $\hat{f}(\mathbf{x}^*) \in \mathbb{C}_0[\mathcal{B}^*]$ .

If we define  $\mathbf{y}(\cdot)$  mapping  $\mathcal{B} \to \mathbb{C}$ , by  $\mathbf{y}(\mathbf{x}) = exp\{-2\pi i \langle \mathbf{y}, \mathbf{x}^* \rangle\}$ , then  $\mathbf{y}(\mathbf{x})$  is a character of  $\mathcal{B}$ . Furthermore, it is easy to see that  $(\mathbf{y}_1 + \mathbf{y}_2)(\mathbf{x}) =$ 

 $\mathbf{y}_1(\mathbf{x}) \cdot \mathbf{y}_2(\mathbf{x})$ . We now have the extension of the Pontryagin Duality Theorem to all UCB (with a basis).

**Theorem 5.9.** If  $\mathcal{B}$  is a UCB, then the mapping  $\gamma_{\mathbf{x}} : \mathcal{B} \to \hat{\mathcal{B}}$ , defined by  $\gamma_{\mathbf{x}}(\mathbf{y}) = \mathbf{y}(\mathbf{x})$ , is an isomorphism of topological groups.

In case  $\mathcal{B} = \mathcal{H}$ , is a Hilbert space, we can replace equation (5.2) by

(5.3) 
$$\hat{f}(\mathbf{x}) = \mathfrak{F}[f](\mathbf{x}) = \int_{\mathcal{H}} \exp\{-2\pi i \langle \mathbf{x}, \mathbf{y} \rangle_{\mathcal{H}}\} f(\mathbf{y}) d\lambda_{\mathcal{H}}(\mathbf{y}),$$

so that  $\mathcal{H}$  is self-dual (as expected). From equation (5.3), we also get the expected result that:

$$\mathfrak{F}\left[\exp\{-\pi |\mathbf{x}|_{\mathcal{H}}^{2}\}\right] = \exp\{-\pi |\mathbf{x}|_{\mathcal{H}}^{2}\}.$$

In closing, we observe that by Theorem 2.31 (see Example 2.32), if we use Gaussian measure on  $\mathbb{R}^{\infty}$ , the dual character groups are  $\mathbb{R}^{\infty}$  and  $\hat{\mathbb{R}}^{\infty} = \mathbb{R}_{0}^{\infty}$ . From this we see that probability measures on  $(\mathbb{R}^{\infty}, \mathfrak{B}[\mathbb{R}^{\infty}])$  induce a different character theory compared to that induced by  $\lambda_{\infty}$  on  $(\mathbb{R}_{I}^{\infty}, \mathfrak{B}[\mathbb{R}_{I}^{\infty}])$ .

<span id="page-51-0"></span>5.4.  $L^p$ -Theory. We can obtain  $L^p[\mathbb{R}_I^{\infty}]$  as in the construction of  $L^1[\mathbb{R}_I^{\infty}]$ . In this section we want to show the power of our approach to measure theory by establishing a version of Young's Theorem for every separable Banach space with a Schauder basis:

**Theorem 5.10.** (Young) Let  $p, q, r \in [1, \infty]$  with

$$\frac{1}{r} = \frac{1}{p} + \frac{1}{q} - 1.$$

If  $f \in L^p[\mathbb{R}_I^{\infty}]$  and  $g \in L^q[\mathbb{R}_I^{\infty}]$ , then the convolution of f and g, f \* g, exists (a.s.), belongs to  $L^r[\mathbb{R}_I^{\infty}]$  and

$$||f * g||_r \le ||f||_p ||g||_q$$
.

Corollary 5.11. Let  $\mathcal{B}$  be a separable Banach space with a Schauder basis and let  $p, q, r \in [1, \infty]$  with

$$\frac{1}{r} = \frac{1}{p} + \frac{1}{q} - 1.$$

If  $f \in L^p[\mathcal{B}]$  and  $g \in L^q[\mathcal{B}]$ , then the convolution of f and g, f \* g, exists (a.s.), belongs to  $L^r[\mathcal{B}]$  and

$$||f * g||_r \le ||f||_p ||g||_q$$
.

In order to prove Theorem 5.10, we first need the appropriate version of Fubini's Theorem. Since  $(\mathbb{R}_I^{\infty}, \mathfrak{L}_I^{\infty}, \lambda_{\infty})$  is a complete  $\sigma$ -finite measure space, a proof of the following may be found in Royden [RO] (see Theorems 19 and 20, pgs. 269-270):

**Theorem 5.12.** (Fubini) If  $f \in L^1[\mathbb{R}_I^{\infty} \times \mathbb{R}_I^{\infty}]$ , then

- (1) for almost all  $x \in \mathbb{R}_I^{\infty}$  the function  $f_x$  defined by  $f_x(y) = f(x,y) \in L^1[\mathbb{R}_I^{\infty}](y)$ :
- (2) for almost all  $y \in \mathbb{R}_I^{\infty}$  the function  $f_y$  defined by  $f_y(x) = f(x,y) \in L^1[\mathbb{R}_I^{\infty}](x)$ :
- (3)  $\int_{\mathbb{R}_I^{\infty}} f(x,y) d\lambda_{\infty}(y) \in L[\mathbb{R}_I^{\infty}](x);$
- (4)  $\int_{\mathbb{R}_I^{\infty}} f(x, y) d\lambda_{\infty}(x) \in L[\mathbb{R}_I^{\infty}](y);$

54

$$\int_{\mathbb{R}_{I}^{\infty} \times \mathbb{R}_{I}^{\infty}} f(x, y) d(\lambda_{\infty} \otimes \lambda_{\infty}) (x, y)$$

$$= \int_{\mathbb{R}_{I}^{\infty}} \left[ \int_{\mathbb{R}_{I}^{\infty}} f(x, y) d\lambda_{\infty}(y) \right] d\lambda_{\infty}(x) = \int_{\mathbb{R}_{I}^{\infty}} \left[ \int_{\mathbb{R}_{I}^{\infty}} f(x, y) d\lambda_{\infty}(x) \right] d\lambda_{\infty}(y).$$

**Theorem 5.13.** Let  $f, g \in L^1[\mathbb{R}_I^{\infty}]$ , then (f \* g)(x) exists (a.s.); that is  $f(y)g(x-y) \in L^1[\mathbb{R}_I^{\infty}]$ . In addition,  $f * g \in L^1[\mathbb{R}_I^{\infty}]$  and

$$||f * g||_1 \le ||f||_1 ||g||_1$$
.

*Proof.* First, it is easy to see that f(y)g(x-y) is a measurable function on  $\mathbb{R}^{\infty}_{I}$ . (There is no change from the case of  $\mathbb{R}^{n}$ .) We can apply Fubini's theorem to get that:

$$\begin{split} &\int_{\mathbb{R}_{I}^{\infty}} \left(f * g\right)(x) d\lambda_{\infty}(x) \\ &= \int_{\mathbb{R}_{I}^{\infty}} d\lambda_{\infty}(x) \left[ \int_{\mathbb{R}_{I}^{\infty}} f(y) g(x-y) d\lambda_{\infty}(y) \right] = \int_{\mathbb{R}_{I}^{\infty}} d\lambda_{\infty}(y) \left[ \int_{\mathbb{R}_{I}^{\infty}} f(y) g(x-y) d\lambda_{\infty}(x) \right] \\ &= \int_{\mathbb{R}_{I}^{\infty}} f(y) d\lambda_{\infty}(y) \cdot \int_{\mathbb{R}_{I}^{\infty}} g(x) d\lambda_{\infty}(x). \end{split}$$

It follows from the last equality that  $||f * g||_1 \le ||f||_1 ||g||_1$ .

#### 5.4.1. Proof of Young's Theorem.

*Proof.* First, assume that f and g are nonnegative and  $||f||_p = ||g||_q = 1$ .

Let 
$$\frac{1}{q'} = 1 - \frac{1}{q}$$
 and  $\frac{1}{p'} = 1 - \frac{1}{p}$ . Now note that

$$\frac{1}{r} + \frac{1}{q'} + \frac{1}{p'} = \frac{1}{r} + \left(1 - \frac{1}{q}\right) + \left(1 - \frac{1}{p}\right) = 1;$$

$$\left(1 - \frac{p}{r}\right)q' = p\left(\frac{1}{p} - \frac{1}{r}\right)q' = p\left(1 - \frac{1}{q}\right)q' = p;$$

$$\left(1 - \frac{q}{r}\right)p' = q\left(\frac{1}{q} - \frac{1}{r}\right)p' = q\left(1 - \frac{1}{p}\right)p' = q.$$

If we use Holder's inequality (for three functions), we can write (f \* g)(x)

as:

$$\begin{split} &(f*g)\left(x\right) = \int_{\mathbb{R}_I^\infty} \left[ f(y)^{p/r} g(x-y)^{q/r} \right] \left[ f(y)^{1-p/r} g(x-y)^{1-q/r} \right] d\lambda_\infty(y) \\ &\leqslant \left[ \int_{\mathbb{R}_I^\infty} f(y)^p g(x-y)^q d\lambda_\infty(y) \right]^{1/r} \left[ \int_{\mathbb{R}_I^\infty} f(y)^{(1-p/r)q'} d\lambda_\infty(y) \right]^{1/q'} \left[ \int_{\mathbb{R}_I^\infty} g(x-y)^{(1-q/r)p'} d\lambda_\infty(y) \right]^{1/p'}. \end{split}$$

This last inequality shows that

$$(f * g)(x) \leqslant \left[ \int_{\mathbb{R}_{I}^{\infty}} f(y)^{p} g(x - y)^{q} d\lambda_{\infty}(y) \right]^{1/r} \Rightarrow$$

$$(f * g)^{r}(x) \leqslant \left[ \int_{\mathbb{R}_{I}^{\infty}} f(y)^{p} g(x - y)^{q} d\lambda_{\infty}(y) \right] \Rightarrow (f * g)^{r}(x) \leqslant (f^{p} * g^{q})(x).$$

From Theorem 5.13, we have  $\|(f*g)^r\|_1 \leqslant \|f^p\|_1 \|g^q\|_1 = 1$ . In the general case, we know that |f|\*|g| exists (a.e.), so that  $|f(y)g(x-y)| \in L^1[\mathbb{R}_I^\infty]$ .

In closing we note that, Beckner [BE] and Brascamp-Lieb [BL] have shown that on  $\mathbb{R}^n$  we can write Young's inequality as  $\|f * g\|_r \leqslant (C_{p,q,r;n})^n \|f\|_p \|g\|_q$ , where  $C_{p,q,r;n} \leq 1$  is sharp. We conjecture that 1 is the sharp constant for  $\mathbb{R}_I^{\infty}$ .

#### <span id="page-55-0"></span>6. Partial Differential Operators (Examples)

In this section, we give examples of strong product and sum vectors for differential operators that have found interest in infinite dimensional analysis.

**Definition 6.1.** For  $x \in \mathbb{R}$ ,  $0 \le y < \infty$  and  $1 < a < \infty$  define  $\bar{g}(x,y)$ ,  $\bar{h}(x)$  by:

$$\bar{g}(x,y) = \exp\left\{-y^a e^{iax}\right\},$$

$$\bar{h}(x) = \begin{cases} \int_0^\infty \bar{g}(x,y) dy, & x \in \left[-\frac{\pi}{2a}, \frac{\pi}{2a}\right], \\ 0 & \text{otherwise} \end{cases}$$

The following properties of  $\bar{g}$  are easy to check:

(1) 
$$\frac{\partial \bar{g}(x,y)}{\partial x} = -iay^a e^{iax} \bar{g}(x,y),$$

(2) 
$$\frac{\partial \bar{g}(x,y)}{\partial y} = -ay^{a-1}e^{iax}\bar{g}(x,y),$$

so that

(3) 
$$iy \frac{\partial \bar{g}(x,y)}{\partial y} = \frac{\partial \bar{g}(x,y)}{\partial x}.$$

It is also easy to see that  $\bar{h}(x)$  is in  $\mathbf{L}^1[\mathbb{R}]$  for  $x \in [-\frac{\pi}{2a}, \frac{\pi}{2a}]$  and,

(6.1) 
$$\frac{d\bar{h}(x)}{dx} = \int_0^\infty \frac{\partial \bar{g}(x,y)}{\partial x} dy = \int_0^\infty iy \frac{\partial \bar{g}(x,y)}{\partial y} dy.$$

Integration by parts in the last expression of equation (6.1) shows that  $\bar{h}'(x) = -i\bar{h}(x)$ , so that  $\bar{h}(x) = \bar{h}(0)e^{-ix}$  for  $x \in [-\frac{\pi}{2a}, \frac{\pi}{2a}]$ . Since  $\bar{h}(0) = \int_0^\infty \exp\{-y^a\}dy$ , an additional integration by parts shows that  $\bar{h}(0) = \Gamma(\frac{1}{a} + 1)$ .

Let  $a = \frac{\pi}{1-\varepsilon}$ ,  $\bar{h}(x) = \bar{h}_{\varepsilon}(x)$ ,  $x \in [-\frac{\pi}{2a}, \frac{\pi}{2a}]$ , where  $0 < \varepsilon \ll 1$ , and define

(6.2) 
$$f_{\varepsilon}(x) = \begin{cases} c \exp\left\{\frac{\varepsilon^2}{|2x|^2 - \varepsilon^2}\right\}, & |x| < \varepsilon/2, \\ 0, & |x| \geqslant \varepsilon/2, \end{cases}$$

where c is the standard normalizing constant. We now define  $\xi(x)=(\bar{h}*f_{\varepsilon})(x)$ , so that  $\operatorname{spt}(\xi)=[-\frac{1}{2},\frac{1}{2}]=I_{\varepsilon}$ . Thus,  $\xi(x)=0,\ x\notin I_{\varepsilon}$  and otherwise,

$$\xi(x) = \int_{-\infty}^{\infty} \bar{h}[x-z] f_{\varepsilon}(z) dz = e^{-ix} \int_{-\infty}^{\infty} e^{iz} f_{\varepsilon}(z) dz = \alpha_{\varepsilon} e^{-ix}.$$

It follows from this that:

$$\alpha_{\varepsilon}^{-1}\xi(ix) = \begin{cases} e^x, & x \in I_{\varepsilon} \\ 0, & x \notin I_{\varepsilon} \end{cases}$$

Define  $\lambda_{\varepsilon} = \lambda$  and,

$$I^\varepsilon = \times_{k=1}^\infty I_\varepsilon, \ I_n^\varepsilon = \times_{k=n+1}^\infty \ \text{and}, \ \lambda_\infty^\varepsilon = \otimes_{k=1}^\infty \lambda_\varepsilon.$$

**Example 6.2.** In this example, we let  $h_k(x_k) = \alpha_{\varepsilon}^{-1}\xi(ix_k)$ , for each  $k \in \mathbb{N}$  so that  $D_k h_k = h_k$ , for  $x_k \in I$ . Let  $L_{\varepsilon}^2[\mathbb{R}_I^n] = L^2[\mathbb{R}_I^n, \lambda_{\infty}^{\varepsilon}]$ . If  $\mathbf{D}^{\infty} = \prod_{k=1}^{\infty} D_k$  and  $f_n \in L_{\varepsilon}^2[\mathbb{R}_I^n] \cap D(\mathbf{D}^{\infty})$ , we can define  $\mathbf{D}^{\infty}$  on  $\mathbb{R}_I^n$  by  $\mathbf{D}^{\infty} f_n(x) = \mathbf{D}^n f_n(x) = \prod_{l=1}^n D_l f_{(n)}(x) \otimes \left( \bigotimes_{l=n+1}^{\infty} h_l \right)$ , (a.s). This operator is well-defined and has a closed densely defined extension to  $L_{\varepsilon}^2[\mathbb{R}_I^{\infty}](h)$ , where  $h = \bigotimes_{k=1}^{\infty} h_k$ . Thus, h is a strong product vector for  $\mathbf{D}^{\infty}$ .

The operator  $\mathbf{D}^{\infty}$  is required if we want to obtain the probability density for a distribution function. (Note, by construction the density can be approximated from below by densities in a finite number of variables.)

In the following example, we construct a general elliptic operator on  $L^2_{\varepsilon}[\mathbb{R}^{\infty}_I]$ .

**Example 6.3.** If  $\nabla = (D_1, D_2, \cdots)$  and  $\sigma_k : \mathbb{R}_I^{\infty} \to \mathbb{R}$  is a bounded analytic function for each  $k \in \mathbb{N}$ , then let  $\sigma(x) = (\sigma_1(x), \sigma_2(x), \cdots)$ . We assume that

$$\left\| \sum_{j,k=1}^{\infty} \sigma_j(x) \sigma_k(x) + \sum_{k=1}^{\infty} b_k(x) \right\|_2 < \infty, \text{ where } b_k(x) = \sum_{j=1}^{\infty} \sigma_j(x) D_j \sigma_k(x).$$

We can now define  $\Delta_{\infty}$  by:

$$\Delta_{\infty} = (\boldsymbol{\sigma}(x) \cdot \nabla)^2 = \sum_{j,k=1}^{\infty} \sigma_j(x) \sigma_k(x) D_j D_k + \sum_{k=1}^{\infty} b_k(x) D_k.$$

For the same version of  $L^2_{\varepsilon}[\mathbb{R}_I^{\infty}]$  as in the last example, if  $g_n \in L^2[\mathbb{R}_I^n] \cap D(\Delta_{\infty})$  and  $c_n(x) = \sum_{j=n+1}^{\infty} \sigma_j(x) [D_j \sigma_k(x)]$ , then  $\Delta_{\infty}$  is defined on  $\mathbb{R}_I^n$  and

$$\Delta_{\infty}g_n(x) = \sum_{j,k=1}^n \sigma_j(x)\sigma_k(x)D_jD_kg_n(x) + \sum_{k=1}^n b_kD_k\psi(x) + c_n(x)g_n(x).$$

In order to obtain the same equation with  $c_n(x) = 0$ , we use the version of  $L^2_{\varepsilon}[\mathbb{R}_I^{\infty}]$  defined with  $h_k = \xi_I(0)$ ,  $k \geq n+1$ . In this case, for  $g_n(x) \in L^2_{\varepsilon}[\mathbb{R}_I^{\infty}] \cap D(\Delta_{\infty})$ , we see that

$$\Delta_{\infty}g_n(x) = \sum_{j,k=1}^n \sigma_j(x)\sigma_k(x)D_jD_kg_n(x) + \sum_{k=1}^n b_k(x)D_kg_n(x).$$

In either case,  $\Delta_{\infty}$  is well-defined for each n and has a closed densely defined extension to  $L^2_{\varepsilon}[\mathbb{R}_I^{\infty}]$  and  $g_n(x) \to g(x)$  implies that  $\lim_{n \to \infty} \Delta_{\infty} g_n(x) = \Delta_{\infty} g(x)$ .

It follows that different versions of  $L^2[\mathbb{R}_I^{\infty}]$  offer advantages for particular types of differential operators. (For other approaches, see [BK], [GZ] and [UM].)

The following special cases have appeared in the literature (all can be obtained from the last example):

(1) The natural infinite dimensional Laplacian:

$$\mathbf{A} = \Delta_{\infty} = \sum_{i=1}^{\infty} \partial^2 / \partial x_i^2.$$

(2) The nonterminating diffusion generator in infinitely many variables (also known as the Ornstein-Uhlenbeck operator):

$$\mathbf{A} = \frac{1}{2}\Delta_{\infty} - B\mathbf{x} \cdot \nabla_{\infty} = \frac{1}{2} \sum_{i=1}^{\infty} \partial^2 / \partial x_i^2 - \sum_{i=1}^{\infty} b_i x_i \partial / \partial x_i.$$

The infinite dimensional Laplacian of Umemura [UM]:

$$\mathbf{A} = \sum_{i=1}^{\infty} \left( \frac{\partial^2}{\partial x_i^2} - \frac{x_i}{c^2} \frac{\partial}{\partial x_i} \right).$$

Berezanskii and Kondratyev ([BK], pages 520-521) have also discussed operators analogous to (2) and (3).

<span id="page-59-0"></span>6.1. **Discussion.** In a very interesting paper, Phillip Duncan Thompson [PDT] used the amplitudes of a set of orthogonal modes as the co-ordinates in an infinite-dimensional phase space. This allowed him to derive the probability distribution for an ensemble of randomly forced two-dimensional viscous flows as the solution of the continuity equation for the phase flow. He obtained the following equation for the probability density:

(6.3) 
$$\frac{\partial \rho}{\partial t} + \sum_{k=1}^{\infty} M_k(\mathbf{x}) \frac{\partial \rho}{\partial x_k} - \nu \sum_{k=1}^{\infty} \frac{\partial}{\partial x_k} \left[ \rho \alpha_k^2 x_k \right] - \sum_{k=1}^{\infty} \frac{\partial^2 \rho}{\partial x_k^2} = 0,$$

where

$$M_k(\mathbf{x}) = \sum_{i=1}^{\infty} \sum_{j=1}^{\infty} \frac{\alpha_j^2 \beta_{ijk}}{\alpha_i \alpha_j \alpha_k} \frac{(\mu_i \mu_j \mu_k)^{\frac{1}{2}}}{\mu_k} x_i x_j.$$

The coefficients  $\beta_{ijk}$  vanishes if any two indices are equal, is invariant under cyclic permutation of indices and reverses sign under non-cyclic permutation of indices, while the coefficients  $\alpha_i$  and  $\mu_i$  are positive constants, determined by the problem. Thompson imposed the natural condition

(6.4) 
$$\int_{-\infty}^{\infty} \cdots \int_{-\infty}^{\infty} \rho(\mathbf{x}, t) \prod_{k=1}^{\infty} dx_k = 1.$$

At that time, he ran into the obvious mathematical criticism because equation (6.4) was meaningless at the time. He also derived the equilibrium density

(6.5) 
$$\rho_0(\mathbf{x}) = C \exp\left\{-\frac{1}{2}\nu \sum_{k=1}^{\infty} \alpha_k^2 x_k^2\right\}.$$

The results in section 2.5, see also equation (2.5), along with those in section 4.4, show that Thomson's paper was prescient.

## 7. Conclusion

<span id="page-60-0"></span>In this paper we provided a reasonable version of Lebesgue measure on R∞, which together with the standard Gaussian measure on R∞, have allowed us to construct natural analogues of Lebesgue and Gaussian measure for every separable Banach space with a Schauder basis. We have extended the Fourier transform to L 1 [R∞, λ∞], L<sup>2</sup> [R∞, λ∞], defined sums and products of unbounded operators, and presented a few constructive examples of partial differential operators in infinitely many variables.

<span id="page-60-1"></span>Acknowledgments. This work could not have been written without the generous help and critical remarks of Professor Frank Jones. We thank Professor Anatoly Vershik for appraising us of recent work by the Russian school and correcting some of our historical remarks. We would also like to thank an anonymous referee for corrections that have improve our presentation and for suggesting that we reconsider our approach to the Fourier transform and discuss its relationship to the Pontryagin duality theory, which led to a complete revision and extension of Section 5.

#### <span id="page-60-2"></span>References

- <span id="page-60-4"></span>[BA] S. Banach Th´eorie des Op´erations lin´eaires, Monografj Matematyczn, Vol. 1, Warsaw, (1932).
- <span id="page-60-3"></span>[BA1] R. Baker "Lebesque measure" on R <sup>∞</sup>, Proc. Amer. Math. Soc. 113 (1991), 1023- 1029.

- <span id="page-61-4"></span>[BA2] R. Baker "Lebesque measure" on R <sup>∞</sup>, II, Proc. Amer. Math. Soc. 132 (2004), 2577-2591.
- <span id="page-61-13"></span><span id="page-61-11"></span>[BE] W. Beckner Inequalities in Fourier Analysis, Ann. of Math. 102 (1975), 159-182.
- [BK] Yu. M. Berezanskii and Yu. G. Kondratyev, Spectral methods in infinite dimensional analysis, Naukova Dumka, Kiev, (Russian) (1988).
- <span id="page-61-12"></span>[BL] H. J. Brascamp and E. H. Lieb Best constants in Young's inequality, its converse, and its generalization to more than three functions, Ad. in Math. 20 (1976), 151-173.
- <span id="page-61-1"></span>[BM] Y. Bakhtin and J. C. Mattingly, Malliavin calculus for infinite-dimensional systems with additive noise, [arXiv:math.PR/0610754v](http://arxiv.org/abs/math/0610754)1 (2006).
- <span id="page-61-8"></span>[DI] J. Diestel , Sequences and Series in Banach Spaces, Grad. Texts in Math. Springer-Verlag, New York, (1984).
- <span id="page-61-0"></span>[DP] G. DaPrato, An Introduction to Infinite-Dimensional Analysis, Springer Berlin, (2006).
- <span id="page-61-10"></span>[DS] N. Dunford and J. T. Schwartz, Linear Operators Part I: General Theory, Wiley Classics edition, Wiley Interscience (1988).
- <span id="page-61-3"></span>[EM] E. O. Elliott and A. P. Morse, General product measures, Trans. Amer. Math. Soc. 110, (1964) 245-283.
- <span id="page-61-5"></span>[GU] A. Guichardet, Symmetric Hilbert Spaces and Related Topics, Lectures Notes in Mathematics, No. 261, Springer-Verlag, New York, (1969).
- <span id="page-61-6"></span>[GZ] T. L. Gill and W. W. Zachary, Banach spaces of von Neumann type, Georgian International Journal of Science Technology 3 (2011), 1-35.
- <span id="page-61-7"></span>[GZ1] T. L. Gill and W. W. Zachary, Feynman operator calculus: The constructive theory, Expositiones Mathematicae 29 (2011), 165-203.
- <span id="page-61-2"></span>[HA] A. Haar, Der Massbegriff in der Theorie der kontinuierlichen Gruppe , Ann. Math., 34 (1933), 147-169.
- <span id="page-61-9"></span>[HHK] H. H. Kuo , Gaussian Measures in Banach Spaces, Lecture Notes in Mathematics 463, Springer, New York (1975).

- <span id="page-62-4"></span>[HI] D. Hill, σ-finite invariant measures on infinite product spaces, Trans. Amer. Math. Soc. 153, (1971) 347-370.
- <span id="page-62-9"></span>[J] F. Jones , Lebesgue Integration on Euclidean Space, Revised Edition, Jones and Bartlett Publishers, Boston (2001).
- <span id="page-62-5"></span>[KA] S. Kakutani, On equivalence of infinite product measures , Ann. Math., 49 (1948), 214-224.
- <span id="page-62-10"></span>[KA1] S. Kaplan Extensions of Pontjagin duality I: Infinite products, Duke Math. J. 15 (1948), 649-659.
- <span id="page-62-11"></span>[KA2] S. Kaplan Extensions of Pontjagin duality II: Direct and inverse sequences, Duke Math. J. 17 (1950), 419-435.
- <span id="page-62-0"></span>[KH] A.B. Kharazishvili On invariant measures in the Hilbert space, Bull. Acad. Sci. Georgian SSR, 114 (1) (1984),4148 (in Russian).
- [KK] K. Kodaira and S. Kakutani, A nonseparable translation-invariant extension of Lebesgue measure space, Ann. Math., 52 (1950), 574-579.
- <span id="page-62-8"></span>[KO] A. N. Kolmogorov, Grundbegriffe der Wahrscheinlichkeitsrechnung, Springer-Verlag, Vienna, (1933).
- <span id="page-62-6"></span>[KP1] A. P. Kirtadze and G. R. Pantsulaia, Invariant measures in the space R<sup>N</sup> , (in Russian) Soobshch. Akad. Nauk Gruzii 141 (1991), 273-276.
- <span id="page-62-1"></span>[KP2] A. P. Kirtadze and G. R. Pantsulaia, Lebesgue nonmeasurable sets and the uniqueness of invariant measures in infinite-dimensional vector spaces, Proc. A. Razmadze Math. Inst. 143 (2007), 95-101.
- <span id="page-62-3"></span>[MO] C. C. Moore, Invariant measures on product spaces, Proc. 5th Berkeley Sym. Math. Stat. & Prob. (Berkeley, 1965) 2 Part 2, 447-459.
- <span id="page-62-2"></span>[OX] J. C. Oxtoby, Invariant measures in groups which are not locally compact, Trans. Amer. Math. Soc. 60, (1946) 215237.
- <span id="page-62-7"></span>[PA] G. Pantsulaia, Invariant and Quasiinvariant Measures in Infinite-Dimensional Topological Vector Spaces , Nova Science Publishers, New York, (2007).

- <span id="page-63-0"></span>[PDT] P. D. Thompson, Some exact statistics of two-dimensional viscous flow with random forcing, Journal of Fluid Mechanics 55 (1972), 711-717.
- <span id="page-63-7"></span>[RE] M. C. Reed, On self-adjointness in infinite tensor product spaces, Journal of Functional Analysis 5 (1970), 94-124.
- <span id="page-63-3"></span>[RH] G. E. Ritter and E. Hewitt, Elliott-Morse measures and Kakutanis dichotomy theorem, Math. Zeitschrift 211, (1992) 247263.
- <span id="page-63-9"></span>[RO] H. L. Royden, Real Analysis, (2nd Ed.) Macmillan Press, New York, (1968).
- <span id="page-63-8"></span>[RU] W. Rudin, Functional Analysis, McGraw-Hill Press, New York, (1973).
- [RU1] W. Rudin, Fourier Analysis on Groups, John Wiley & Sons, New York, (1990).
- [ST] L. Streit, Test function spaces for direct product representations, Commun. Math. Phys. 4 (1967), 22-31.
- <span id="page-63-2"></span>[SU] V.N. Sudakov, Linear sets with quasi-invariant measure, Dokl.Akad.Nauk SSSR; 127 (1959), 524-525 (in Russian).
- [T] P. D. Thompson Some exact statistics of two-dimensional viscous flow with random forcing, J. Fluid Mech. 55 (1972), 711-717.
- <span id="page-63-10"></span>[UM] Y. Umemura, On the infinite dimensional Laplacian operator, J. Math. Kyoto Univ. 4 (1964/1965), 477-492.
- <span id="page-63-1"></span>[V] A. M. Vershik, Duality in the theory of measure in linear spaces, (English translation): Sov. Math. Dokl. 7, (1967) 1210-1214.
- <span id="page-63-4"></span>[V1] A. M. Vershik, Does there exist the Lebesgue measure in the infinite-dimensional space?, Proceedings of the Steklov Institute of Mathematics, 259 (2007), 248-272.
- <span id="page-63-5"></span>[V2] A. M. Vershik, The behavior of Laplace transform of the invariant measure on the hyperspace of high dimension, J. Fixed Point Theory Appl., 3 (2008), 317-329.
- <span id="page-63-6"></span>[V3] D. L. Vandev, Invariant measures for the continual Cartan subgroup, J. Functional Analysis, 255 (2008), 2661-2682.

- <span id="page-64-3"></span>[VA] A. M. Vershik, Action of groups on an infinite product measure, (in Russian) Papers presented at the Fifth Balkan Mathematical Congress (Belgrade, 1974). Math. Balkanica 4 (1974), 643647.
- <span id="page-64-1"></span>[VN1] J. von Neumann, The uniqueness of Haar's measure, Rec. Math. (Mat. Sbornik) N.S., 1 (1936), 721-734.
- <span id="page-64-4"></span>[VN2] J. von Neumann, On infinite direct products, Compositio Mathematica, 6 (1938), 1-77.
- <span id="page-64-2"></span>[WE] A. Weil, L'int´egration dans les groupes topologiques et ses applications, Actualit´es Scientifiques et Industrielles, no. 869, Paris, 1940.
- <span id="page-64-5"></span>[WSRM] N. Wiener, A. Siegel, W. Rankin and W. T. Martin, Differential Space, Quantum Systems, and Prediction, M. I. T. Press, Cambridge, MA, (1966).
- <span id="page-64-6"></span>[YA] Y. Yamasaki, Measures on infinite-dimensional spaces, World Scientific, Singapore, (1985).
- <span id="page-64-0"></span>[YA1] Y. Yamasaki, Translationally invariant measure on the infinite-dimensional vector space, Publ. Res. Inst. Math. Sci. 16 (3) (1980), 693720.
- (Tepper L. Gill) Department of Mathematics, Physics and E&CE, Howard University, Washington DC 20059, USA, E-mail : tgill@howard.edu
- (Gogi R. Pantsulaia) Department of Mathematics, Georgian Technical University, Tbilisi 0175, Georgia; I. Vekua Institute of Applied Mathematics, Tbilisi State University, Tbilisi 0143, Georgia E-mail : g:pantsulaia@gtu:ge
- (Woodford W. Zachary) Department of Mathematics and E&CE, Howard University, Washington DC 20059, USA, E-mail : wwzachary@earthlink.net