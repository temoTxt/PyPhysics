### THE JONES STRONG DISTRIBUTION BANACH SPACES

#### TEPPER L. GILL

Abstract. In this note, we introduce a new class of separable Banach spaces, SD<sup>p</sup> [R n ], 1 6 p 6 ∞, which contain each L p -space as a dense continuous and compact embedding. They also contain the nonabsolutely integrable functions and the space of test functions D[R n ], as dense continuous embeddings. These spaces have the remarkable property that, for any multi-index α, kD <sup>α</sup>ukSD <sup>=</sup> <sup>k</sup>ukSD, where <sup>D</sup> is the distributional derivative. We call them Jones strong distribution Banach spaces because of the crucial role played by two special functions introduced in his book (see [\[J\]](#page-16-0), page 249). After constructing the spaces, we discuss their basic properties and their relationship to D[R n ] and D ′ [R n ]. As an application, we obtain new a priori bounds for the Navier-Stokes equation.

### Introduction

The theory of distributions is based on the action of linear functionals on a space of test functions. In the original approach of Schwartz [\[SC\]](#page-17-0), both the test functions and the linear functionals have a natural topological vector space structure, which is not normable. For those interested in applications,

<sup>1991</sup> *Mathematics Subject Classification.* Primary (45) Secondary(46) .

*Key words and phrases.* Banach spaces, test functions, distributions, Navier-Stokes equation .

this is an inconvenience, requiring additional study and effort. Thus, in most applied contexts, the restricted class of Banach spaces due to Sobolev have proved useful (see Leoni [\[GL\]](#page-16-1)). In this case, the base space of functions are of Lebesgue type, which also have some limitations. The study of problems in the foundations of relativistic quantum theory, path integrals and nonlinear analysis have led to the need for a Banach space structure for nonabsolutely integrable functions. Recent research has uncovered a new class of separable Banach spaces, which are natural for the class of nonabsolutely integrable functions and have provided the first rigorous foundation for the Feynman path integral formulation of quantum mechanics (see [\[GZ\]](#page-16-2)). The purpose of this note is to introduce a another class of Banach spaces which contain the nonabsolutely integrable functions, but also contains the Schwartz test function space as dense and continuous embeddings.

## 1. The Jones Spaces

We begin with the construction of a special class of functions in C<sup>∞</sup> c [R n (see Jones, [\[J\]](#page-16-0) page 249).

# 1.1. The remarkable Jones functions.

Definition 1.1. *For* x ∈ R, 0 ≤ y < ∞ *and* 1 < a < ∞*, define the Jone's functions* g(x, y), h(x) *by:*

$$g(x,y) = \exp\left\{-y^a e^{iax}\right\},\,$$

$$h(x) = \begin{cases} \int_0^\infty g(x,y)dy, & x \in (-\frac{\pi}{2a}, \frac{\pi}{2a}) \\ \\ 0, & otherwise. \end{cases}$$

The following properties of g are easy to check:

(1)

$$\frac{\partial g(x,y)}{\partial x} = -iay^a e^{iax} g(x,y),$$

(2)

$$\frac{\partial g(x,y)}{\partial y} = -ay^{a-1}e^{iax}g(x,y),$$

so that

(3)

$$iy \frac{\partial g(x,y)}{\partial y} = \frac{\partial g(x,y)}{\partial x}.$$

It is also easy to see that  $h(x) \in L^1[-\frac{\pi}{2a}, \frac{\pi}{2a}]$  and,

(1.1) 
$$\frac{dh(x)}{dx} = \int_0^\infty \frac{\partial g(x,y)}{\partial x} dy = \int_0^\infty iy \frac{\partial g(x,y)}{\partial y} dy.$$

Integration by parts in the last expression in (1.1) shows that h'(x) = -ih(x), so that  $h(x) = h(0)e^{-ix}$  for  $x \in (-\frac{\pi}{2a}, \frac{\pi}{2a})$ . Since  $h(0) = \int_0^\infty \exp\{-y^a\}dy$ , an additional integration by parts shows that  $h(0) = \Gamma(\frac{1}{a} + 1)$ . For each  $k \in \mathbb{N}$  let  $a = a_k = 3 \times 2^{k-1}$ ,  $h(x) = h_k(x)$ ,  $x \in (-\frac{\pi}{2a_k}, \frac{\pi}{2a_k})$  and set  $\varepsilon_k = \frac{\pi}{4a_k}$ .

Let  $\mathbb{Q}$  be the set of rational numbers in  $\mathbb{R}$  and for each  $x^i \in \mathbb{Q}$ , define

$$f_k^i(x) = f_k(x - x^i) = \begin{cases} c_k \exp\left\{\frac{\varepsilon_k^2}{|x - x^i|^2 - \varepsilon_k^2}\right\}, & |x - x^i| < \varepsilon_k, \\ 0, & |x - x^i| \ge \varepsilon_k, \end{cases}$$

where  $c_k$  is the standard normalizing constant. It is clear that the support,  $\operatorname{spt}(f_k^i) \subset [-\varepsilon_k, \varepsilon_k] = [-\frac{\pi}{4a_k}, \frac{\pi}{4a_k}] = I_k^i$ .

Now set  $\chi_k^i(x) = (f_k^i * h_k)(x)$ , so that  $\operatorname{spt}(\chi_k^i) \subset [-\frac{\pi}{2^{k+1}}, \frac{\pi}{2^{k+1}}]$ . For  $x \in \operatorname{spt}(\chi_k^i)$ , we can also write  $\chi_k^i(x) = \chi_k(x - x^i)$  as:

$$\chi_k^i(x) = \int_{I_k^i} f_k \left[ \left( x - x^i \right) - z \right] h_k(z) dz$$

$$= \int_{I_k^i} h_k \left[ \left( x - x^i \right) - z \right] f_k(z) dz = e^{-i\left( x - x^i \right)} \int_{I_k^i} e^{iz} f_k(z) dz.$$

Thus, if  $\alpha_{k,i} = \int_{I_k^i} e^{iz} f_k^i(z) dz$ , we can now define:

$$\xi_k^i(x) = \alpha_{k,i}^{-1} \chi_k^i(-x) = \begin{cases} \frac{1}{n} e^{i(x-x^i)}, & x \in I_k^i \\ 0, & x \notin I_k^i, \end{cases}$$

so that  $\left|\xi_k^i(x)\right| < \frac{1}{n}$ .

1.2. **The Construction.** To construct our space on  $\mathbb{R}^n$ , let  $\mathbb{Q}^n$  be the set of all vectors  $\mathbf{x}$  in  $\mathbb{R}^n$ , such that for each j, the component  $x_j$  is rational. Since this is a countable dense set in  $\mathbb{R}^n$ , we can arrange it as  $\mathbb{Q}^n = \{\mathbf{x}^1, \mathbf{x}^2, \mathbf{x}^3, \cdots\}$ . For each k and i, let  $\mathbf{B}_k(\mathbf{x}^i)$  be the closed cube centered at  $\mathbf{x}^i$  with edge  $\frac{\pi}{a_k}$  and diagonal of length  $r_k = \frac{\pi}{a_k} \sqrt{n}$ .

We choose the natural order which maps  $\mathbb{N} \times \mathbb{N}$  bijectively to  $\mathbb{N}$ :

$$\{(1,1), (2,1), (1,2), (1,3), (2,2), (3,1), (3,2), (2,3), \dots\}$$

and let  $\{\mathbf{B}_m, m \in \mathbb{N}\}$  be the set of closed cubes  $\mathbf{B}_k(\mathbf{x}^i)$  with  $(k, i) \in \mathbb{N} \times \mathbb{N}$  and  $\mathbf{x}^i \in \mathbb{Q}^n$ . For each  $\mathbf{x} \in \mathbf{B}_m$ ,  $\mathbf{x} = (x_1, x_2, \dots, x_n)$ , we define  $\mathcal{E}_k(\mathbf{x})$  by :

$$\mathcal{E}_k(\mathbf{x}) = \left(\xi_k^i(x_1), \xi_k^i(x_2) \dots \xi_k^i(x_n)\right) \text{ with}$$

$$|\mathcal{E}_k(\mathbf{x})| < 1, \ \mathbf{x} \in \prod_{j=1}^n I_k^i \text{ and } \mathcal{E}_k(\mathbf{x}) = 0, \ \mathbf{x} \notin \prod_{j=1}^n I_k^i.$$

It is easy to see that  $\mathcal{E}_k(\mathbf{x})$  is in  $L^p[\mathbb{R}^n]^n = \mathbf{L}^p[\mathbb{R}^n]$  for  $1 \leq p \leq \infty$ . Define  $F_k(\cdot)$  on  $\mathbf{L}^p[\mathbb{R}^n]$  by

$$F_k(f) = \int_{\mathbb{R}^n} \mathcal{E}_k(\mathbf{x}) \cdot f(\mathbf{x}) d\lambda_n(\mathbf{x}).$$

It is clear that  $F_k(\cdot)$  is a bounded linear functional on  $\mathbf{L}^p[\mathbb{R}^n]$  for each k with  $||F_k|| \le 1$ . Furthermore, if  $F_k(f) = 0$  for all k, f = 0 so that  $\{F_k\}$  is a fundamental sequence of functionals on  $\mathbf{L}^p[\mathbb{R}^n]$  for  $1 \le p \le \infty$ .

Set  $t_k = \frac{1}{2^k}$  so that  $\sum_{k=1}^{\infty} t_k = 1$  and define a inner product  $(\cdot)$  on  $\mathbf{L}^1[\mathbb{R}^n]$  by

$$(f,g) = \sum_{k=1}^{\infty} t_k \left[ \int_{\mathbb{R}^n} \mathcal{E}_k(\mathbf{x}) \cdot f(\mathbf{x}) d\lambda_n(\mathbf{x}) \right] \overline{\left[ \int_{\mathbb{R}^n} \mathcal{E}_k(\mathbf{y}) \cdot g(\mathbf{y}) d\lambda_n(\mathbf{y}) \right]}.$$

The completion of  $\mathbf{L}^1[\mathbb{R}^n]$  with the above inner product is a Hilbert space, which we denote as  $SD^2[\mathbb{R}^n]$ .

# **Theorem 1.2.** For each $p, 1 \leq p \leq \infty$ , we have:

- (1) The space  $SD^2[\mathbb{R}^n] \supset \mathbf{L}^p[\mathbb{R}^n]$  as a continuous, dense and <u>compact</u> embedding.
- (2) The space  $SD^2[\mathbb{R}^n] \supset \mathfrak{M}[\mathbb{R}^n]$ , the space of finitely additive measures on  $\mathbb{R}^n$ , as a continuous dense and compact embedding.

Proof. Since  $SD^2[\mathbb{R}^n]$  contains  $L^1[\mathbb{R}^n]$  densely, to prove (1), we need only show that  $L^q[\mathbb{R}^n] \subset SD^2[\mathbb{R}^n]$  for  $q \neq 1$ . Let  $f \in L^q[\mathbb{R}^n]$  and  $q < \infty$ . By construction, for every j,  $|\mathcal{E}_j(\mathbf{x})| < \frac{1}{\sqrt{n}}$  so that there is a constant C = C(q), with  $|\mathcal{E}_j(\mathbf{x})|^q \leq C|\mathcal{E}_j(\mathbf{x})|$ . It follows that:

$$||f||_{SD^{2}} = \left[ \sum_{k=1}^{\infty} t_{k} \left| \int_{\mathbb{R}^{n}} \mathcal{E}_{k}(\mathbf{x}) f(\mathbf{x}) d\lambda_{n}(\mathbf{x}) \right|^{\frac{2q}{q}} \right]^{1/2}$$

$$\leq C \left[ \sum_{k=1}^{\infty} t_{k} \left( \int_{\mathbb{R}^{n}} \mathcal{E}_{k}(\mathbf{x}) |f(\mathbf{x})|^{q} d\lambda_{n}(\mathbf{x}) \right)^{\frac{2}{q}} \right]^{1/2}$$

$$\leq C \sup_{k} \left( \int_{\mathbb{R}^{n}} \mathcal{E}_{k}(\mathbf{x}) |f(\mathbf{x})|^{q} d\lambda_{n}(\mathbf{x}) \right)^{\frac{1}{q}} \leq C ||f||_{q}.$$

Hence,  $f \in SD^2[\mathbb{R}^n]$ . For  $q = \infty$ , first note that  $vol(\mathbf{B}_k)^2 \leq \left[\frac{1}{2\sqrt{n}}\right]^{2n}$ , so we have

$$||f||_{SD^2} = \left[ \sum_{k=1}^{\infty} t_k \left| \int_{\mathbb{R}^n} \mathcal{E}_k(\mathbf{x}) f(\mathbf{x}) d\lambda_n(\mathbf{x}) \right|^2 \right]^{1/2}$$

$$\leq \left[ \left[ \sum_{k=1}^{\infty} t_k [vol(\mathbf{B}_k)]^2 \right] [ess \sup |f|]^2 \right]^{1/2} \leq \left[ \frac{1}{2\sqrt{n}} \right]^n ||f||_{\infty}.$$

Thus  $f \in SD^2[\mathbb{R}^n]$ , and  $\mathbf{L}^{\infty}[\mathbb{R}^n] \subset SD^2[\mathbb{R}^n]$ . To prove compactness, suppose  $\{f_n\}$  is any weakly convergent sequence in  $\mathbf{L}^p$ ,  $1 \le p \le \infty$  with limit f. Since  $\mathcal{E}_k \in \mathbf{L}^q$ , 1/p + 1/q = 1,

$$\int_{\mathbb{R}^n} \mathcal{E}_k(\mathbf{x}) \cdot [f_n(\mathbf{x}) - f(\mathbf{x})] d\lambda_n(\mathbf{x}) \to 0$$

for each k. It follows that  $\{f_n\}$  converges strongly to f in  $SD^2$ .

To prove (2), we note that 
$$\mathfrak{M}[\mathbb{R}^n] = \mathbf{L}^1[\mathbb{R}^n]^{**} \subset SD^2[\mathbb{R}^n].$$

**Definition 1.3.** We call  $SD^2[\mathbb{R}^n]$  the Jones strong distribution Hilbert space on  $\mathbb{R}^n$ .

In order to justify our definition, let α be a multi-index of nonnegative integers, α = (α1, α2, · · · αk), with |α| = P<sup>k</sup> <sup>j</sup>=1 α<sup>j</sup> . If D denotes the standard partial differential operator, let

$$D^{\alpha} = D^{\alpha_1} D^{\alpha_2} \cdots D^{\alpha_k}.$$

Theorem 1.4. *Let* D[R n ] *be* C<sup>∞</sup> c [R n ] *equipped with the standard locally convex topology (test functions).*

- (1) *If* φ<sup>n</sup> → φ *in* D[R n ]*, then* <sup>φ</sup><sup>n</sup> <sup>→</sup> <sup>φ</sup> *in the norm topology of* SD<sup>2</sup> [R n ]*, so that* D[R n ] <sup>⊂</sup> SD<sup>2</sup> [R n ] *as a continuous dense embedding.*
- (2) *If* T ∈ D′ [R n ]*, then* <sup>T</sup> <sup>∈</sup> SD<sup>2</sup> [R n ′ *, so that* D′ [R n ] <sup>⊂</sup> SD<sup>2</sup> [R n ′ *as a continuous dense embedding.*
- (3) *For any* f, g <sup>∈</sup> SD<sup>2</sup> [R n ] *and any multi-index* α, (Dαf, g)SD = (−i) <sup>α</sup>(f, g)SD*.*
- (4) *The functions* S[R n ]*, of rapid decrease at infinity are contained in* SD<sup>2</sup> [R n ] *continuous embedding, so that* S ′ [R n ] <sup>⊂</sup> SD<sup>2</sup> [R n ] ′ *.*

*Proof.* To prove (1), suppose that φ<sup>n</sup> → φ in D[R n ]. By definition, there exist a compact set K ⊂ R n , which is the support of <sup>φ</sup><sup>n</sup> <sup>−</sup> <sup>φ</sup> and <sup>D</sup>αφ<sup>n</sup> converges to <sup>D</sup>α<sup>φ</sup> uniformly on <sup>K</sup> for every multi-index <sup>α</sup>. Let {EK<sup>j</sup> } be

the set of all  $\mathcal{E}_j$ , with support  $K_j \subset K$ . If  $\alpha$  is a multi-index, we have:

$$\lim_{k \to \infty} \|D^{\alpha} \phi_{k} - D^{\alpha} \phi\|_{SD}$$

$$= \lim_{k \to \infty} \left\{ \sum_{j=1}^{\infty} t_{K_{j}} \left| \int_{\mathbb{R}^{n}} \mathcal{E}_{K_{j}}(\mathbf{x}) \cdot \left[ D^{\alpha} \phi_{k}(\mathbf{x}) - D^{\alpha} \phi(\mathbf{x}) \right] d\lambda_{n}(\mathbf{x}) \right|^{2} \right\}^{1/2}$$

$$\leq M \lim_{k \to \infty, \mathbf{x} \in K} |D^{\alpha} \phi_{k}(\mathbf{x}) - D^{\alpha} \phi(\mathbf{x})| = 0.$$

Thus, since  $\alpha$  is arbitrary, we see that  $\mathcal{D}[\mathbb{R}^n] \subset SD^2[\mathbb{R}^n]$  as a continuous embedding. Since  $\mathbb{C}_c^{\infty}[\mathbb{R}^n]$  is dense in  $\mathbf{L}^1[\mathbb{R}^n]$ ,  $\mathcal{D}[\mathbb{R}^n]$  is dense in  $SD^2[\mathbb{R}^n]$ . To prove (2) we note that, as  $\mathcal{D}[\mathbb{R}^n]$  is a dense locally convex subspace of  $SD^2[\mathbb{R}^n]$ , by a corollary of the Hahn-Banach Theorem every continuous linear functional, T defined on  $\mathcal{D}[\mathbb{R}^n]$ , can be extended to a continuous linear functional on  $SD^2[\mathbb{R}^n]$ . By the Riesz representation Theorem, every continuous linear functional T, defined on  $SD^2[\mathbb{R}^n]$  is of the form  $T(f) = (f,g)_{SD}$ , for some  $g \in SD^2[\mathbb{R}^n]$ . Thus,  $T \in SD^2[\mathbb{R}^n]'$  and, by the identification  $T \leftrightarrow g$  for each T in  $\mathcal{D}'[\mathbb{R}^n]$ , we can map  $\mathcal{D}'[\mathbb{R}^n]$  into  $SD^2[\mathbb{R}^n]$  as a continuous dense embedding.

To prove (3), recall that each  $\mathcal{E}_k \in \mathbb{C}_c^{\infty}[\mathbb{R}^n]$  so that, for any  $f \in SD^2[\mathbb{R}^n]$ ,

$$\int_{\mathbb{R}^n} \mathcal{E}_k(\mathbf{x}) \cdot D^{\alpha} f(\mathbf{x}) d\lambda_n(\mathbf{x}) = (-1)^{|\alpha|} \int_{\mathbb{R}^n} D^{\alpha} \mathcal{E}_k(\mathbf{x}) \cdot f(\mathbf{x}) d\lambda_n(\mathbf{x}).$$

An easy calculation shows that:

$$(-1)^{|\alpha|} \int_{\mathbb{R}^n} D^{\alpha} \mathcal{E}_k(\mathbf{x}) \cdot f(\mathbf{x}) d\lambda_n(\mathbf{x}) = (-i)^{|\alpha|} \int_{\mathbb{R}^n} \mathcal{E}_k(\mathbf{x}) \cdot f(\mathbf{x}) d\lambda_n(\mathbf{x}).$$

It now follows that, for any  $\mathbf{g} \in SD^2[\mathbb{R}^n]$ ,  $(D^{\alpha}f, \mathbf{g})_{SD^2} = (-i)^{|\alpha|}(f, \mathbf{g})_{SD^2}$ .

1.3. Functions of Bounded Variation. The integral due to Henstock [\[HS\]](#page-16-3) and Kurzweil [\[KW\]](#page-17-1) (HK-integral), is the easiest to learn and best known of those integrals that integrate nonabsolutely integrable functions, which also extend the Lebesgue integral. A good discussion of all the standard types of integrals, with comparisons among themselves and the Lebesgue integral, can be found in Gordon [\[GO\]](#page-16-4).

The objective of this section is to show that every HK-integrable function is in SD<sup>2</sup> [R n ]. To do this, we need to discuss a certain class of functions of bounded variation. For functions defined on R, the definition of bounded variation is unique. However, for functions on R n , n ≥ 2, there are a number of distinct definitions.

The functions of bounded variation in the sense of Cesari are well known to analysts working in partial differential equations and geometric measure theory (see Leoni [\[GL\]](#page-16-1)).

Definition 1.5. *A function* f ∈ L 1 [R n ] *is said to be of bounded variation in the sense of Cesari or* f ∈ BVc[R n ]*, if* f ∈ L 1 [R n ] *and each* i, 1 ≤ i ≤ n*, there exists a signed Radon measure* µ<sup>i</sup> *, such that*

$$\int_{\mathbb{R}^n} f(\mathbf{x}) \frac{\partial \phi(\mathbf{x})}{\partial x_i} d\lambda_n(\mathbf{x}) = -\int_{\mathbb{R}^n} \phi(\mathbf{x}) d\mu_i(\mathbf{x}),$$

*for all* φ ∈ C<sup>∞</sup> 0 [R n ]*.*

The functions of bounded variation in the sense of Vitali [TY1], are well known to applied mathematicians and engineers with interest in error estimates associated with research in control theory, financial derivatives, high speed networks, robotics and in the calculation of certain integrals. (See, for example [KAA], [NI], [PT] or [PTR] and references therein.) For the general definition, see Yeong ([TY1], p. 175). We present a definition that is sufficient for continuously differentiable functions.

**Definition 1.6.** A function f with continuous partials is said to be of bounded variation in the sense of Vitali or  $f \in BV_v[\mathbb{R}^n]$  if for all intervals  $(a_i, b_i)$ ,  $1 \le i \le n$ ,

$$V(f) = \int_{a_1}^{b_1} \cdots \int_{a_n}^{b_n} \left| \frac{\partial^n f(\mathbf{x})}{\partial x_1 \partial x_2 \cdots \partial x_n} \right| d\lambda_n(\mathbf{x}) < \infty.$$

**Definition 1.7.** We define  $BV_{v,0}[\mathbb{R}^n]$  by:

$$BV_{v,0}[\mathbb{R}^n] = \{ f(\mathbf{x}) \in BV_v[\mathbb{R}^n] : f(\mathbf{x}) \to 0, \text{ as } x_i \to -\infty \},$$

where  $x_i$  is any component of  $\mathbf{x}$ .

The following two theorems may be found in [TY1]. (See p. 184 and 187, where the first is used to prove the second.) If  $[a_i, b_i] \subset \mathbb{R}$ , we define  $[\mathbf{a}, \mathbf{b}] \in \mathbb{R}^n$  by  $[\mathbf{a}, \mathbf{b}] = \prod_{k=1}^n [a_i, b_i]$ . (The notation (RS) means Riemann-Stieltjes.)

**Theorem 1.8.** Let f be HK-integrable on  $[\mathbf{a}, \mathbf{b}]$  and let  $g \in BV_{v,0}[\mathbb{R}^n]$ , then fg is HK-integrable and

$$(HK) \int_{[\mathbf{a}, \mathbf{b}]} f(\mathbf{x}) g(\mathbf{x}) d\lambda_n(\mathbf{x}) = (RS) \int_{[\mathbf{a}, \mathbf{b}]} \left\{ (HK) \int_{[\mathbf{a}, \mathbf{x}]} f(\mathbf{y}) d\lambda_n(\mathbf{y}) \right\} dg(\mathbf{x})$$

**Theorem 1.9.** Let f be HK-integrable on  $[\mathbf{a}, \mathbf{b}]$  and let  $g \in BV_{v,0}[\mathbb{R}^n]$ , then fg is HK-integrable and

$$\left| (HK) \int_{[\mathbf{a}, \mathbf{b}]} f(\mathbf{x}) g(\mathbf{x}) d\lambda_n(\mathbf{x}) \right| \le ||f||_D V_{[\mathbf{a}, \mathbf{b}]}(g)$$

**Lemma 1.10.** The space  $HK[\mathbb{R}^n]$ , of all HK-integrable functions is contained in  $SD^2[\mathbb{R}^n]$ .

*Proof.* Since each  $\mathcal{E}_k(\mathbf{x})$  is continuous and differentiable,  $\mathcal{E}_k(\mathbf{x}) \in BV_{v,0}[\mathbb{R}^n]$ , so that for  $f \in HK[\mathbb{R}^n]$ ,

$$||f||_{\mathbf{SD}^2}^2 = \sum_{k=1}^{\infty} t_k \left| \int_{\mathbb{R}^n} \mathcal{E}_k(\mathbf{x}) \cdot f(\mathbf{x}) d\mathbf{x} \right|^2 \leqslant \sup_k \left| \int_{\mathbb{R}^n} \mathcal{E}_k(\mathbf{x}) \cdot f(\mathbf{x}) d\mathbf{x} \right|^2$$
  
$$\leqslant ||f||_{HK}^2 \left[ \sup_k V(\mathcal{E}_k) \right]^2 < \infty.$$

It follows that  $f \in SD^2[\mathbb{R}^n]$ .

1.4. The General Case,  $SD^p$ ,  $1 \le p \le \infty$ . To construct  $SD^p[\mathbb{R}^n]$  for all p and for  $\mathbf{f} \in \mathbf{L}^p$ , define:

$$\|\mathbf{f}\|_{SD^p} = \begin{cases} \left\{ \sum_{k=1}^{\infty} t_k \left| \int_{\mathbb{R}^n} \mathcal{E}_k(\mathbf{x}) \cdot \mathbf{f}(\mathbf{x}) d\lambda_n(\mathbf{x}) \right|^p \right\}^{1/p}, 1 \leqslant p < \infty, \\ \sup_{k \geqslant 1} \left| \int_{\mathbb{R}^n} \mathcal{E}_k(\mathbf{x}) \cdot \mathbf{f}(\mathbf{x}) d\lambda_n(\mathbf{x}) \right|, p = \infty. \end{cases}$$

It is easy to see that  $\|\cdot\|_{SD^p}$  defines a norm on  $\mathbf{L}^p$ . If  $SD^p$  is the completion of  $\mathbf{L}^p$  with respect to this norm, we have:

**Theorem 1.11.** For each  $q, 1 \leq q \leq \infty$ ,  $SD^p[\mathbb{R}^n] \supset \mathbf{L}^q[\mathbb{R}^n]$  as dense continuous embeddings.

*Proof.* As in the previous theorem, by construction  $SD^p[\mathbb{R}^n]$  contains  $\mathbf{L}^p[\mathbb{R}^n]$  densely, so we need only show that  $SD^p[\mathbb{R}^n] \supset \mathbf{L}^q[\mathbb{R}^n]$  for  $q \neq p$ . First, suppose that  $p < \infty$ . If  $\mathbf{f} \in \mathbf{L}^q[\mathbb{R}^n]$  and  $q < \infty$ , we have

$$\|\mathbf{f}\|_{SD^{p}} = \left[\sum_{k=1}^{\infty} t_{k} \left| \int_{\mathbb{R}^{n}} \mathcal{E}_{k}(\mathbf{x}) \cdot \mathbf{f}(\mathbf{x}) d\lambda_{n}(\mathbf{x}) \right|^{\frac{qp}{q}} \right]^{1/p}$$

$$\leq \left[\sum_{k=1}^{\infty} t_{k} \left( \int_{\mathbb{R}^{n}} |\mathcal{E}_{k}(\mathbf{x})|^{q} |\mathbf{f}(\mathbf{x})|^{q} d\lambda_{n}(\mathbf{x}) \right)^{\frac{p}{q}} \right]^{1/p}$$

$$\leq \sup_{k} \left( \int_{\mathbb{R}^{n}} |\mathcal{E}_{k}(\mathbf{x})|^{q} |\mathbf{f}(\mathbf{x})|^{q} d\lambda_{n}(\mathbf{x}) \right)^{\frac{1}{q}} \leq \|f\|_{q}.$$

Hence,  $f \in SD^p[\mathbb{R}^n]$ . For  $q = \infty$ , we have

$$\|\mathbf{f}\|_{SD^p} = \left[\sum_{k=1}^{\infty} t_k \left| \int_{\mathbb{R}^n} \mathcal{E}_k(\mathbf{x}) \cdot \mathbf{f}(\mathbf{x}) d\lambda_n(\mathbf{x}) \right|^p \right]^{1/p}$$

$$\leq \left[ \left[\sum_{k=1}^{\infty} t_k [vol(\mathbf{B}_k)]^p \right] [ess \sup |\mathbf{f}|]^p \right]^{1/p} \leq M \|f\|_{\infty}.$$

Thus  $\mathbf{f} \in SD^p[\mathbb{R}^n]$ , and  $\mathbf{L}^{\infty}[\mathbb{R}^n] \subset SD^p[\mathbb{R}^n]$ . The case  $p = \infty$  is obvious.  $\square$ 

**Theorem 1.12.** For  $SD^p$ ,  $1 \le p \le \infty$ , we have:

- (1) If  $p^{-1} + q^{-1} = 1$ , then the dual space of  $SD^p[\mathbb{R}^n]$  is  $SD^q[\mathbb{R}^n]$ .
- (2) The test function space  $\mathcal{D}[\mathbb{R}^n]$  is contain in  $SD^p[\mathbb{R}^n]$  as a continuous dense embedding.
- (3) If K is a weakly compact subset of  $\mathbf{L}^p[\mathbb{R}^n]$ , it is a strongly compact subset of  $SD^p[\mathbb{R}^n]$ .
- (4) The space  $SD^{\infty}[\mathbb{R}^n] \subset SD^p[\mathbb{R}^n]$ .

## 2. Application

Let [L 2 (R<sup>3</sup> )]<sup>3</sup> be the Hilbert space of square integrable functions on R<sup>3</sup> , let H[R 3 ] be the completion of the set of functions in u ∈ C<sup>∞</sup> 0 [R 3 ] 3 | ∇ · u = 0 which vanish at infinity with respect to the inner product of [L 2 (R 3 )]3 . The classical Navier-Stokes initial-value problem (on R <sup>3</sup> and all T > 0) is to find a function u : [0, T] × R <sup>3</sup> <sup>→</sup> <sup>R</sup> <sup>3</sup> and <sup>p</sup> : [0, T] <sup>×</sup> <sup>R</sup> <sup>3</sup> <sup>→</sup> <sup>R</sup> such that

$$\partial_t \mathbf{u} + (\mathbf{u} \cdot \nabla) \mathbf{u} - \nu \Delta \mathbf{u} + \nabla p = \mathbf{f}(t) \text{ in } (0, T) \times \mathbb{R}^3,$$

$$(2.1) \qquad \nabla \cdot \mathbf{u} = 0 \text{ in } (0, T) \times \mathbb{R}^3 \text{ (in the weak sense)},$$

$$\mathbf{u}(0, \mathbf{x}) = \mathbf{u}_0(\mathbf{x}) \text{ in } \mathbb{R}^3.$$

The equations describe the time evolution of the fluid velocity u(x, t) and the pressure p of an incompressible viscous homogeneous Newtonian fluid with constant viscosity coefficient ν in terms of a given initial velocity u0(x) and given external body forces f(x, t).

Let P be the (Leray) orthogonal projection of (L 2 [R 3 ])<sup>3</sup> onto H[R 3 ] and define the Stokes operator by: Au =: <sup>−</sup>P∆u, for <sup>u</sup> <sup>∈</sup> <sup>D</sup>(A) <sup>⊂</sup> <sup>H</sup><sup>2</sup> [R 3 ], the domain of A. If we apply P to equation (2.1), with B(u, u) = P(u · ∇)u, we can recast equation (2.1) into the standard form:

(2.2) 
$$\partial_t \mathbf{u} = -\nu \mathbf{A} \mathbf{u} - B(\mathbf{u}, \mathbf{u}) + \mathbb{P} \mathbf{f}(t) \text{ in } (0, T) \times \mathbb{R}^3,$$
$$\mathbf{u}(0, \mathbf{x}) = \mathbf{u}_0(\mathbf{x}) \text{ in } \mathbb{R}^3,$$

where the orthogonal complement of  $\mathbb{H}$  relative to  $\{L^2(\mathbb{R}^3)\}^3$ ,  $\{\mathbf{v}: \mathbf{v} = \nabla q, q \in \mathbb{H}^1[\mathbb{R}^3]\}$ , is used to eliminate the pressure term (see Galdi [GA] or [[SY], [T1],[T2]]).

**Definition 2.1.** We say that a velocity vector field in  $\mathbb{R}^3$  is reasonable if for  $0 \le t < \infty$ , there is a continuous function m(t) > 0, depending only on t and a constant  $M_0$ , which may depend on  $\mathbf{u}(0)$  and f, such that

$$0 < m(t) \leq \|\mathbf{u}(t)\|_{\mathbb{H}} \leq M_0.$$

The above definition formalizes the requirement that the fluid has nonzero, but bounded positive definite energy. However, this condition still allows the velocity to approach zero at infinity in a weaker norm.

2.1. The Nonlinear Term: A Priori Estimates. The difficulty in proving the existence and uniqueness of global-in-time strong solutions for equation (2.2) is directly linked to the problem of getting good a priori estimates for the nonlinear term  $B(\mathbf{u}, \mathbf{u})$ . Let  $\mathbb{H}_{sd}$  be the closure of  $D(\mathbf{A}) \cap SD^2[\mathbb{R}^3]$  in the  $SD^2$  norm.

**Theorem 2.2.** If **A** is the Stokes operator and  $\mathbf{u}(\mathbf{x},t) \in D(\mathbf{A})$  is a reasonable vector field, then

(1) 
$$\langle \nu \mathbf{A} \mathbf{u}, \mathbf{u} \rangle_{\mathbb{H}_{sd}} = \nu \|\mathbf{u}\|_{\mathbb{H}_{sd}}^2$$
.

(2) For  $\mathbf{u}(\mathbf{x},t) \in \mathbf{SD}^2 \cap D(\mathbf{A})$  and each  $t \in [0,\infty)$ , there exists a constant  $M = M(\mathbf{u}(\mathbf{x},0)) > 0$ , such that

(2.3) 
$$\left| \langle B(\mathbf{u}, \mathbf{u}), \mathbf{u} \rangle_{\mathbb{H}_{sd}} \right| \le M \|\mathbf{u}\|_{\mathbb{H}_{sd}}^3.$$

(3)

$$(2.4) \left| \langle B(\mathbf{u}, \mathbf{v}), \mathbf{w} \rangle_{\mathbb{H}_{sd}} \right| \leq M \|\mathbf{u}\|_{\mathbb{H}_{sd}} \|\mathbf{w}\|_{\mathbb{H}_{sd}} \|\mathbf{v}\|_{\mathbb{H}_{sd}}.$$

(4)

$$(2.5) \quad \max\{\|B(\mathbf{u}, \mathbf{v})\|_{\mathbb{H}_{sd}}, \ \|B(\mathbf{v}, \mathbf{u})\|_{\mathbb{H}_{sd}}\} \leqslant M \|\mathbf{u}\|_{\mathbb{H}_{sd}} \|\mathbf{v}\|_{\mathbb{H}_{sd}}.$$

*Proof.* From the definition of the inner product, we have

$$\langle \nu \mathbf{A} \mathbf{u}, \mathbf{u} \rangle_{\mathbb{H}_{sd}} = \nu \sum_{k=1}^{\infty} t_k \left[ \int_{\mathbb{R}^n} \mathcal{E}_k(\mathbf{x}) \cdot \mathbf{A} \mathbf{u}(\mathbf{x}) d\lambda_n(\mathbf{x}) \right] \left[ \int_{\mathbb{R}^n} \mathcal{E}_k(\mathbf{y}) \cdot \mathbf{u}(\mathbf{y}) d\lambda_n(\mathbf{y}) \right].$$

Using the fact that  $\mathbf{u} \in D(\mathbf{A})$ , it follows that

$$\int_{\mathbb{R}^n} \mathcal{E}_k(\mathbf{y}) \cdot \partial_{y_j}^2 \mathbf{u}(\mathbf{y}) d\lambda_n(\mathbf{y}) = \int_{\mathbb{R}^n} \partial_{y_j}^2 \mathcal{E}_k(\mathbf{y}) \cdot \mathbf{u}(\mathbf{y}) d\lambda_n(\mathbf{y})$$

$$= \int_{I_i} \partial_{y_j}^2 \left( \xi_l^i(y_1), \xi_l^i(y_2), \dots, \xi_l^i(y_n) \right) \cdot \mathbf{u}(\mathbf{y}) d\lambda_n(\mathbf{y}) = (-i)^2 \int_{\mathbb{R}^n} \mathcal{E}_k(\mathbf{y}) \cdot \mathbf{u}(\mathbf{y}) d\lambda_n(\mathbf{y}).$$

Using this in the above equation and summing on j, we have  $(\mathbf{A} = -\mathbb{P}\Delta)$ 

$$\int_{\mathbb{R}^n} \mathcal{E}_k(\mathbf{y}) \cdot \mathbf{A} \mathbf{u}(\mathbf{y}) d\lambda_n(\mathbf{y}) = \int_{\mathbb{R}^n} \mathcal{E}_k(\mathbf{y}) \cdot \mathbf{u}(\mathbf{y}) d\lambda_n(\mathbf{y}).$$

It follows that

$$\langle \mathbf{A}\mathbf{u}, \mathbf{u} \rangle_{\mathbb{H}_{sd}}$$

$$= \sum_{k=1}^{\infty} t_k \left[ \int_{\mathbb{R}^n} \mathcal{E}_k(\mathbf{x}) \cdot \mathbf{u}(\mathbf{x}) d\lambda_n(\mathbf{x}) \right] \left[ \int_{\mathbb{R}^n} \mathcal{E}_k(\mathbf{y}) \cdot \mathbf{u}(\mathbf{y}) d\lambda_n(\mathbf{y}) \right]$$

$$= \|\mathbf{u}\|_{\mathbb{H}_{sd}}^2.$$

This proves (1). To prove (2), let  $\vec{\delta}(\mathbf{x}) = (\delta(x_1), \dots \delta_k(x_3))$  be the *n*-dimensional Dirac delta function and set  $\hat{\varepsilon} = \|\vec{\delta}(\mathbf{x})\|_{\mathbb{H}_{sd}}$ . We start with

$$b(\mathbf{u}, \mathbf{u}, \mathcal{E}_k) = \left| \langle B(\mathbf{u}, \mathbf{u}), \mathcal{E}_k \rangle_{\mathbb{H}_{sd}} \right| = \left| \int_{\mathbb{R}^3} \left( \mathbf{u}(\mathbf{x}) \cdot \nabla \right) \mathbf{u}(\mathbf{x}) \cdot \mathcal{E}_k(\mathbf{x}) d\lambda_n(\mathbf{x}) \right|$$

and integrate by parts, to get

$$\left| \int_{\mathbb{R}^3} \left\{ \sum_{i=1}^3 u_i(\mathbf{x})^2 \mathcal{E}_k^i(\mathbf{x}) d\lambda_n(\mathbf{x}) \right\} \right| \leqslant \sup_k \|\mathcal{E}_k\|_{\infty} \|\mathbf{u}\|_{\mathbb{H}}^2 \leq \|\mathbf{u}\|_{\mathbb{H}}^2.$$

Since **u** is reasonable, there is a constant  $\bar{M}$  depending on **u**(0) and f, such that  $\|\mathbf{u}\|_2^2 \leq \bar{M} \|\mathbf{u}\|_{\mathbb{H}_{sd}}^2$ . We now have

$$\begin{aligned} & \left| \langle B(\mathbf{u}, \mathbf{u}), \mathbf{u} \rangle_{\mathbb{H}_{sd}} \right| \\ &= \left| \sum_{k=1}^{\infty} t_k \left[ \int_{\mathbb{R}^3} \left( \mathbf{u}(\mathbf{x}) \cdot \nabla \right) \mathbf{u}(\mathbf{x}) \cdot \mathcal{E}_k(\mathbf{x}) d\lambda_n(\mathbf{x}) \right] \left[ \int_{\mathbb{R}^3} \mathbf{u}(\mathbf{y}) \cdot \mathcal{E}_k(\mathbf{y}) d\lambda_n(\mathbf{y}) \right] \right| \\ &\leqslant \bar{M} \hat{\varepsilon}^{-2} \left\| \mathbf{u} \right\|_{\mathbb{H}_{sd}}^2 \left| \sum_{k=1}^{\infty} t_k \left[ \int_{\mathbb{R}^3} \vec{\delta}(\mathbf{x}) \cdot \mathcal{E}_k(\mathbf{x}) d\lambda_n(\mathbf{x}) \right] \left[ \int_{\mathbb{R}^3} \mathbf{u}(\mathbf{y}) \cdot \mathcal{E}_k(\mathbf{y}) d\lambda_n(\mathbf{y}) \right] \right| \\ &\leqslant M \left\| \mathbf{u} \right\|_{\mathbb{H}_{sd}}^3, \end{aligned}$$

where  $M = \bar{M}\hat{\varepsilon}^{-1}$  and the third line above follows from Schwartz's inequality. The proofs of (3) and (4) are easy.

To compare our results, a typical bound available in the  $\mathbb{H}$  (or energy) norm for equation (2.5) can be found in Sell and You [SY] (see page 366):

$$\max \left\{ \|B(\mathbf{u}, \mathbf{v})\|_{\mathbb{H}}, \|B(\mathbf{v}, \mathbf{u})\|_{\mathbb{H}} \right\} \leqslant C_0 \|\mathbf{A}^{5/8} \mathbf{u}\|_{\mathbb{H}} \|\mathbf{A}^{5/8} \mathbf{v}\|_{\mathbb{H}}.$$

## 3. Conclusion

We have constructed a new class of separable Banach spaces, SD<sup>p</sup> [R n ], 1 6 p 6 ∞, which contain each L p -space as a dense continuous and compact embedding. These spaces have the remarkable property that, for any multi-index α, <sup>k</sup>DαukSD <sup>=</sup> <sup>k</sup>ukSD. We have shown that our spaces contain the nonabsolutely integrable functions and the space of test functions D[R n ], as a dense continuous embedding. We have discussed their basic properties and their relationship to D′ [R n ], S[R n ] and S ′ [R n ]. As an application, we have obtained new bounds for the nonlinear term of the Navier-Stokes equation.

### References

- <span id="page-16-5"></span>[GA] G. P. Galdi, *An introduction to the mathematical theory of the Navier-Stokes equations,* 2nd Edition, Vol. II, Springer Tracts in Natural Philosophy, Vol. 39 Springer, New York, 1997.
- <span id="page-16-1"></span>[GL] G. Leoni, *A First Course in Sobolev Spaces,* AMS Graduate Studies in Math. 105, Providence, R.I, 2009.
- <span id="page-16-4"></span>[GO] R. A. Gordon, *The Integrals of Lebesgue, Denjoy, Perron and Henstock*, Graduate Studies in Mathematics, Vol. 4, Amer. Math. Soc., (1994).
- <span id="page-16-2"></span>[GZ] T.L. Gill and W.W. Zachary, *A new class of Banach spaces*, J. Phys. A: Math. Theor. 41 (2008) 1-15, doi:10.1088/1751-8113/41/49/495206.
- <span id="page-16-3"></span>[HS] R. Henstock, *The General Theory of Integration*, Clarendon Press, Oxford, (1991).
- <span id="page-16-0"></span>[J] F. Jones, *Lebesgue Integration on Euclidean Space*, Revised Edition, Jones and Bartlett Publishers, Boston, (2001).

- <span id="page-17-3"></span>[KAA] V. Koltchinskii, C.T. Abdallah, M. Ariola, P. Dorato and D. Panchenko, *Improved Sample Complexity Estimates for Statistical Learning Control of Uncertain Systems,* IEEE Trans. Automatic Control 45 (2000), 2383-2388.
- <span id="page-17-1"></span>[KW] J. Kurzweil, *Nichtabsolut konvergente Integrale*, Teubner-Texte z¨ur Mathematik, Band 26, Teubner Verlagsgesellschaft, Leipzig, (1980).
- <span id="page-17-4"></span>[NI] H. Niederreiter, *Random Number Generation and Quasi-Monte Carlo Methods,* SIAM, (1992).
- <span id="page-17-5"></span>[PT] A. Papageorgiou and J.G. Traub, *Faster Evaluation of Multidimesional Integrals,* Computers in Physics, Nov., (1997), 574-578.
- <span id="page-17-6"></span>[PTR] S. Paskov and J.G. Traub, *Faster Valuation of Financial Derivatives,* Journal of Portfolio Management, 22 (1995), 113-120.
- [RS] M. Reed and B. Simon, *Methods of Modern Mathematical Physics I: Functional Analysis*, Academic Press, New York, (1972).
- <span id="page-17-0"></span>[RU] W. Rudin, *Functional Analysis*, McGraw-Hill Press, New York, (1973).
- <span id="page-17-7"></span>[SC] L. Schwartz, *Th´eorie des Destributions ,* Hermann, Paris, 1966.
- [SY] G. R. Sell and Y. You, *Dynamics of evolutionary equations,* Applied Mathematical Sciences, Vol. 143, Springer, New York, 2002.
- <span id="page-17-8"></span>[T1] R. Temam, *Navier-Stokes Equations, Theory and Numerical Analysis,* AMS Chelsea Pub., Providence, RI, 2001.
- <span id="page-17-9"></span>[T2] R. Temam, *Infinite dimensional dynamical systems in mechanics and physics,* Applied Mathematical Sciences, Vol. 68, Springer, New York, 1988.
- [TY] L. Tuo-Yeong, *Some full descriptive charaterizations of the Henstock-Kurweil integral in Euclidean space*, Czechoslovak Math. Journal , 55 (2005), 625-637.
- <span id="page-17-2"></span>[TY1] L. Tuo-Yeong, *Henstock-Kurweil Integration on Euclidean Spaces,* Series in Real Analysis-Vol 12 World Scientific, New Jersey, (2011).
- [YS] K. Yosida, *Functional Analysis,* second ed. Springer-Verlag, New York, (1968).

### THE JONES STRONG DISTRIBUTION BANACH SPACES 19

(Tepper L. Gill) Department of Mathematics, Physics and E&CE, Howard University, Washington DC 20059, USA, *E-mail :* tgill@howard.edu