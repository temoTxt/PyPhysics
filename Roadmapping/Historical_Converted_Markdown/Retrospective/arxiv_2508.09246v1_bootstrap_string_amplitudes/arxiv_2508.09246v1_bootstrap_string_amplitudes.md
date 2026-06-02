## **Strings from Almost Nothing**

Clifford Cheung,<sup>1</sup> Grant N. Remmen,<sup>2</sup> Francesco Sciotti,<sup>3</sup> and Michele Tarquini<sup>1</sup>

<sup>1</sup>*Walter Burke Institute for Theoretical Physics, California Institute of Technology, Pasadena, CA 91125, USA* <sup>2</sup>*Center for Cosmology and Particle Physics, Department of Physics, New York University, New York, NY 10003, USA* 3 *IFAE and BIST, Universitat Autònoma de Barcelona, 08193 Bellaterra, Barcelona, Spain*

We argue that string theory emerges inevitably from a few simple assumptions about physical scattering. Consistency alone requires that all tree-level four-point scattering amplitudes exhibit vanishing residues at prescribed values of the momentum transfer. Assuming ultrasoft high-energy behavior, we then prove that the space of minimally consistent amplitudes, whose residues exhibit these mandated zeros and nothing more, are precisely the amplitudes of Veneziano and Virasoro-Shapiro, thus establishing the uniqueness of strings. Similar logic also applies to five-point scattering.

**Introduction.**—Are the laws of physics described by string theory? This question has elicited much controversy, despite being a well-posed scientific hypothesis that is either true or false. The good news is that string theory has a sharp prediction: the emergence of very specific higher-spin excitations in gravitational scattering at or below the Planck scale. The bad news is that experiment—the ultimate arbiter of all theories will not offer any insight into this question in the foreseeable future, or perhaps ever.

In the absence of experimental input, it is tempting to commit faithfully either to string theory or its negation. A less devoted approach is to ask to what extent string theory might follow from more modest or widely accepted assumptions. One instead explores the space of consistent theories by "bootstrapping" all possible scattering amplitudes consistent with fundamental principles like unitarity, locality, and Lorentz invariance. Recent years have witnessed substantial progress in this area, for example on effective field theory constraints [\[1–](#page-4-0)[35\]](#page-5-0), the primal bootstrap [\[36](#page-5-1)[–50\]](#page-5-2), semidefinite programming and null constraints [\[51](#page-5-3)[–70\]](#page-6-0), and bootstrapping string amplitudes and their deformations [\[71](#page-6-1)[–94\]](#page-7-0).

In this paper, we argue for the uniqueness of string theory from the fewest assumptions possible. All amplitude properties—the spectrum of masses and spins, high-energy behavior, and even its explicit mathematical form—are outputs of this bootstrap construction.

At the heart of our analysis is a remarkable fact about *all tree-level, crossing-symmetric, Lorentz invariant, unitary amplitudes*. Without loss of generality, the spectrum of masses squared and spins is defined by unknown functions *µ*(*n*) and *J*(*n*) indexed by integer level *n*, while the leading Regge behavior is defined by an unknown trajectory *α*(*t*) [\[95\]](#page-7-1). We will prove that the residue at large *n* vanishes when *α*(*t*) crosses certain negative integers. These "Regge zeros" are mandated by consistency and apply to planar and nonplanar amplitudes at four-point and beyond for any *µ*(*n*), *J*(*n*), and *α*(*t*). To single out string theory we further assume:

*i*) *Ultrasoftness.* The Regge trajectory *α*(*t*) is bijective for *t <* 0 [\[96\]](#page-7-2) such that the amplitude falls off arbitrarily at infinite momentum transfer, so *α*(*t*) → −∞ for *t* → −∞ [\[97\]](#page-7-3). As depicted in Fig. [1,](#page-0-0) ultrasoft-

<span id="page-0-1"></span>![](_page_0_Figure_11.jpeg)

<span id="page-0-0"></span>Figure 1. Schematic behavior of Regge trajectories for quantum field theory, string theory, and its *q*-deformations.

ness is exhibited by string theory [\[98\]](#page-7-4) but not its *q*deformations [\[99\]](#page-7-5), which diverge at finite *t*, nor quantum field theory, whose softness is bounded [\[100,](#page-7-6) [101\]](#page-7-7). Ultrasoftness is equivalent to superpolynomial boundedness [\[87\]](#page-6-2) together with bijectivity.

*ii*) *Minimal Zeros.* All zeros of the residues are Regge zeros mandated by self-consistency. Extra zeros are not physically inconsistent, but they are manifestly nonminimal because they go beyond what is required.

Given these modest conditions, our logic is simple: ultrasoft Regge behavior dictates an infinite set of points at which the residues of the amplitude are zero. Bootstrapping the space of objects that exhibit these required zeros and nothing more, we land precisely on the Veneziano and Virasoro-Shapiro amplitudes of string theory.

Note that our approach is quite distinct from that of Ref. [\[88\]](#page-6-3), which argued for string uniqueness in tree-level four-point planar scattering while relying on a critical assumption about the distribution of amplitude zeros. Though conjectured to be universal, this property is now understood to be violated quite generally [\[87\]](#page-6-2), leaving the question of string uniqueness an open one until now.

**Regge Zeros.**—Consider the four-point amplitude *A*(*s, t*) for identical color-ordered scalars interacting via tree-level exchanges. This amplitude exhibits planar crossing symmetry, *A*(*s, t*) = *A*(*t, s*), with simple poles at  $s,t=\mu(n)$ , where  $\mu(n)$  parameterizes an arbitrary mass-squared spectrum indexed by the integer level  $n\geq 0$ . Without loss of generality we take  $\mu(n)$  to be monotonically increasing and parameterize the Regge limit at fixed t by some function  $\alpha(t)$ , which specifies the leading Regge trajectory via

<span id="page-1-0"></span>
$$A(s,t) \stackrel{s \to \infty}{\sim} f(t)(-s)^{\alpha(t)}. \tag{1}$$

The prefactor f(t) is analytic except for simple poles at  $t = \mu(n)$ , whose residues are polynomials in s of degree J(n), the maximal spin at level n [102]. Hence, the spectrum and the leading Regge trajectory satisfy

<span id="page-1-8"></span>
$$\alpha(\mu(n)) = J(n)$$
 for integer  $n \ge 0$ , (2)

which is applicable without any loss of generality.

These meager assumptions imply surprisingly stringent conditions on any tree-level, color-ordered, four-point amplitude. To see why, let us compute the "big arc integral" of A(s,t) along a circle of large radius  $s_*$ ,

<span id="page-1-1"></span>
$$\oint^{s_*} ds \, A(s,t) \sim f(t) s_*^{\alpha(t)+1} \frac{\sin \pi \alpha(t)}{\alpha(t)+1}. \tag{3}$$

This follows from Eq. (1) via explicit computation [103], but is intuitive to understand. If  $\alpha(t)$  is an integer, then  $s^{\alpha(t)}$  is a monomial whose contour integral is one if  $\alpha(t) = -1$  and zero if  $\alpha(t) \neq -1$ . Either way, the contour integral is *independent* of the arc radius  $s_*$ . Deforming the contour inward gives a sum over enclosed residues,

<span id="page-1-2"></span>
$$\oint^{s_*} ds A(s,t) \sim \sum_{n=0}^{n_*} R(n,t), \tag{4}$$

where  $R(n,t) = \lim_{s \to \mu(n)} (\mu(n) - s) A(s,t)$ . This result looks very much dependent on the radius of the arc  $s_* \sim \mu(n_*)$ . To reconcile these facts, it must be that the residues vanish for integer  $\alpha(t)$ . To show this, we differentiate Eqs. (3) and (4) with respect to  $n_*$ , approximating the discrete sum over n with a continuous integral over dn, yielding

<span id="page-1-3"></span>
$$R(n,t) \sim f(t)\mu'(n)\mu(n)^{\alpha(t)}\sin\pi\alpha(t).$$
 (5)

This holds in the Regge limit, which corresponds mathematically to  $|t| \ll \mu(n)$ , while also assuming  $\mu'(n) \ll \mu(n)$  so the sum is well approximated by an integral.

The  $\sin \pi \alpha(t)$  in Eq. (5) implies that R(n,t)=0 when  $\alpha(t)$  is an integer, subject to two caveats. First, since R(n,t) is nonsingular and f(t) has poles at  $t=\mu(n')$  for any integer  $n'\geq 0$ , those poles must cancel exactly against corresponding zeros in  $\sin \pi \alpha(t)$ . Thus, a genuine zero of Eq. (5) requires the off-resonance condition  $t\neq \mu(n')$  for all  $n'\geq 0$ . Second, unitarity implies that  $R(n,t)=\sum_{\ell=0}^{\infty}a_{n,\ell}G_{\ell}^{(D)}(1+2t/\mu(n))$  with  $a_{n,\ell}\geq 0$ . For  $\ell>0$ , the D-dimensional Gegenbauer polynomials satisfy  $G_{\ell}^{(D)}(x)>G_{\ell}^{(D)}(1)>0$  when x>1. Since  $\mu(n), J(n)>0$  at large n, we see that R(n,t)>1

R(n,0)>0 for t>0, so this region has no residue zeros. In summary, in the Regge limit of  $|t|\ll \mu(n)$ , the residue R(n,t) vanishes on the support of Regge zeros,

<span id="page-1-4"></span>
$$\alpha(t) + r = 0, (6)$$

for any integer r such that the real roots of Eq. (6) satisfy t<0 and  $t\neq \mu(n')$  for any integer  $n'\geq 0$ . These Regge zeros are required by consistency in any unitary, tree-level, color-ordered, four-point amplitude. Since Eq. (6) holds in the  $|t|\ll \mu(n)$  Regge limit, any Regge zero of R(n,t) will also be a Regge zero of R(n',t) for n'>n since the approximation of the Regge limit only improves as the level grows.

Softer Regge behavior implies more Regge zeros. For example, the Veneziano amplitude of string theory scales as  $\Gamma(-s)\Gamma(-t)/\Gamma(-s-t) \stackrel{s\to\infty}{\sim} s^t$ , with unboundedly fast falloff for sufficiently negative t. The expected Regge zeros at t+r=0 appear in the string residues,

$$R_{\rm str}(n,t) = \frac{1}{n!} \prod_{r=1}^{n} (t+r) \sim n^{t}.$$
 (7)

Regge zeros are ubiquitous [104], but they are often difficult to ascertain from analytic expressions. Consider three examples: the residues of the simplest deformed worldsheet amplitude in the family of Ref. [105],  $\int_0^1 u^{-s-1} (1-u)^{-t-1} e^{gu(1-u)} du$  [106],

<span id="page-1-5"></span>
$$R_{\rm str}(n,t) \,_{2}F_{2} \left[ \begin{array}{c} -n,-t \\ -\frac{n+t}{2}, \frac{1-n-t}{2} \end{array}; \frac{g}{4} \right],$$
 (8)

the three-parameter amplitudes explored in Ref. [87],

<span id="page-1-6"></span>
$$R_{\rm str}(n,t) \, {}_{3}F_{2} \left[ \begin{array}{c} -n, -t, -c_{0} - c_{1}(n+t) \\ -\frac{n+t}{2}, \frac{1-n-t}{2} \end{array}; \lambda \right],$$
 (9)

and one of the bespoke amplitudes of Ref. [79],

<span id="page-1-7"></span>
$$2(n+\delta)(R_{\rm str}(n,-\delta+\sqrt{t})+R_{\rm str}(n,-\delta-\sqrt{t})). \quad (10)$$

Here  $R_{\rm str}(n,t)$  appears in all of these residues, but its zeros are either cancelled by poles or scrambled by changes of arguments. Nevertheless, Eqs. (8), (9), and (10) are all polynomials in t that exhibit their expected Regge zeros at large n, as emerges quite vividly via direct numerical evaluation in Fig. 2. Since Eqs. (8) and (9) have the same spectrum and Regge behavior as the Veneziano amplitude, they also have Regge zeros that asymptote to t+r=0. In contrast, Eq. (10) describes a modified spectrum,  $\mu(n)=(\delta+n)^2$ , and Regge behavior,  $A(s,t)\stackrel{s\to\infty}{\sim} (\sqrt{s})^{-\delta+\sqrt{t}}+(\sqrt{s})^{-\delta-\sqrt{t}}$ , so its Regge zeros are at  $t=(\delta-r)^2$  for t off-resonance, as expected.

Equation (5) implies a quite constraining residue scaling,  $R(n,t) \sim \mu(n)^{\alpha(t)}$ . In fact, this property rules out any amplitude with a single Regge trajectory. By definition, the residues of any such amplitude describe a single state of spin n state at level n,

![](_page_2_Figure_1.jpeg)

<span id="page-2-0"></span>Figure 2. The residues of tree-level four-point amplitudes exhibit Regge zeros at asymptotically large level n. Shown here are the real components of the zeros in t for the residues defined in Eqs. (8), (9), and (10), for levels n = 1 to 80.

so  $R(n,t)=c(n)\,G_n^{(D)}(1+2t/\mu(n))$ . This implies that  $\partial_t^k\log R(n,t)|_{t=0}\sim (n^2/\mu(n))^k$  at large n, modulo multiplicative factors. However, from the residue scaling we know that  $\partial_t^k\log R(n,t)|_{t=0}\sim \log \mu(n)$ , indicating the same scaling with n for all k, which is a contradiction. Similar conclusions were obtained via different methods in Refs. [107, 108]. All of these no-go results assume that the single Regge trajectory extends to infinity without any additional high-energy states.

Uniqueness Argument.—We are now equipped to derive a uniqueness principle for string theory, assuming i) ultrasoft Regge behavior, and ii) residue zeros comprised only of Regge zeros required by consistency. The former imposes an infinite number of Regge zeros at asymptotically large level n. The latter forbids any extra zeros that are not Regge zeros, even at finite n. As noted earlier, Regge zeros necessarily accumulate for increasing n, so ii) implies a nested structure whereby every zero of R(n,t) is also a zero of R(n',t) for n' > n.

To begin, we take d/dn of the log of Eq. (5) to obtain

<span id="page-2-1"></span>
$$\alpha(t)\frac{\mu'(n)}{\mu(n)} + \frac{\mu''(n)}{\mu'(n)} \sim \frac{d\log R(n,t)}{dn} \sim \log \frac{R(n,t)}{R(n-1,t)}, \quad (11)$$

where in the final step we have approximated the derivative as a finite difference. Equation (11) implies that, up to a trivial normalization and shift, the Regge trajectory is  $\alpha(t) \sim \log P(n,t)$ , where P(n,t) = R(n,t)/R(n-1,t) is a rational function of t. Assumption ii) says that the residues R(n,t) and R(n-1,t) only exhibit Regge zeros, and we argued earlier that any Regge zeros of R(n-1,t) are also Regge zeros of R(n,t), thus canceling in the ratio and leaving P(n,t) to be a polynomial in t. Without loss of generality,  $P(n,t) \sim P(t) + \cdots$ , where P(t) is leading in the large-n expansion and the ellipses are subleading.

For nonconstant P(t), at sufficiently large t this polynomial scales as  $P(t) \sim t^k$  for some k>0. This implies a logarithmic Regge trajectory  $\alpha(t) \sim k \log t$ . Plugging into Eq. (2), we learn that  $\mu(n) \sim e^{J(n)/k}$ , which is precisely the q-deformed spectrum of the Coon amplitude [99], which has  $P(t) \sim 1 + (q-1)t$  [109]. Since k>0, we see that  $\alpha(t) \to \infty$  as  $t\to -\infty$ , in violation of assumption i). These logarithmic trajectories violate our assumptions, so we do not consider them further.

For constant P(t), we retain the ellipses of  $P(n,t) \sim P(t) + \cdots$  to get the leading t-dependence of the Regge trajectory. This happens in the Veneziano amplitude, where  $P(n,t) \sim 1 + t/n$ . As the subleading terms are polynomial,  $\alpha(t)$  is as well. We focus on this case next.

All residues exhibit the Regge zeros in Eq. (6) for  $|t| \ll \mu(n)$ , but assumption ii) enforces the stronger condition that these Regge zeros comprise all the zeros of the residue for all t, so R(n,t) is a product of factors of  $\alpha(t)+r$ . Since unitarity restricts all Regge zeros to the region t<0, assumption i) implies that  $\alpha(t)<\alpha(0)$ , which on the Regge zero requires  $r=-\alpha(t)>-\alpha(0)$ . Here we assume a minimal value of r=1, since any finite offset will just shift  $\alpha(t)$  in a way that leaves our conclusions unchanged. Concretely, we define [110]

<span id="page-2-2"></span>
$$R(n,t) = c(n) \prod_{r=1}^{N(n)} (\alpha(t) + r)^h,$$
 (12)

where h is a constant multiplicity of each Regge zero, c(n) is an arbitrary normalization, and we extend the range of r up to  $N(n) = J(n)/h \deg(\alpha)$ . From Eq. (12) we see that  $R(n,t) \sim R_{\rm str}(N(n),\alpha(t))^h \sim N(n)^{h\alpha(t)}$ , while Eq. (5) implies that  $R(n,t) \sim \mu(n)^{\alpha(t)}$ . Comparing these representations, we learn that  $\mu(n) \sim N(n)^h \sim J(n)^h$ . Plugging into Eq. (2), we determine the Regge trajectory,  $\alpha(t) \sim t^{1/h}$ . We argued earlier that this function is a polynomial, so h=1. Hence, we arrive at a linear spectrum and Regge trajectory,

<span id="page-2-3"></span>
$$\mu(n) \sim J(n)$$
 and  $\alpha(t) \sim t$ , (13)

which is one of the main results of this paper. In conclusion, our assumptions i) and ii) imply the Chew-Frautschi scaling and linear Regge behavior that are the defining properties and historical raison d'être [98, 111–113] of string theory. Note that it is straightforward to also derive that  $J(n) \sim n$ , which we detail in App. A.

With  $J(n) \sim n$  and  $\alpha(t) \sim t$ , the residue ansatz in Eq. (12) precisely reduces to the "level truncation" ansatz proposed in Refs. [91, 92], which is solved analytically to obtain a three-parameter space of amplitudes. Within this family, the only object that exhibits arbitrarily soft Regge behavior is  $\Gamma(-s)\Gamma(-t)/\Gamma(-s-t)$ , which is precisely the Veneziano amplitude [114]. Thus, our

assumptions i) and ii) actually uniquely bootstrap the Veneziano amplitude in its full mathematical form!

Assumption ii) is our least conservative condition, but it can be relaxed by including extraneous zeros in the residue in Eq. (12). For example, if Eq. (12) is augmented by a level-independent factor g(t), as might occur for external states with spin, all of our arguments still apply. For a level-dependent factor, g(n,t), our logic persists if the number of extraneous zeros is parametrically smaller than the number of Regge zeros.

Beyond Planar.—Up until now we have focused on planar amplitudes. However, our logic generalizes straightforwardly to the case of nonplanar scattering. Consider the four-point amplitude A(s,t) for identical massless scalars interacting via tree-level exchanges. The amplitude is crossing symmetric under permutations of s,t,u and has poles at  $s,t,u=\mu(n)$ . The Regge behavior at fixed t is parameterized by

$$A(s,t) \stackrel{s \to \infty}{\sim} f(t)(-s)^{\frac{\alpha(t)}{2}}(-u)^{\frac{\alpha(t)}{2}}, \tag{14}$$

which is crossing symmetric under the exchange of s and u = -s - t. As before, we compute the big arc integral along a circle of radius  $s_*$  in the s plane,

$$\oint^{s_*} d(su) A(s,t) \sim f(t) s_*^{\alpha(t)+2} \frac{\sin \frac{1}{2} \pi \alpha(t)}{\alpha(t) + 2},$$
 (15)

and deform the contour to equate with a weighted sum of residues,  $\sum_{n}^{n_*} (\mu(n) + t/2) R(n,t)$ . Differentiating with respect to  $n_*$ , we find that the residues take the form

<span id="page-3-0"></span>
$$R(n,t) \sim f(t)\mu'(n)\mu(n)^{\alpha(t)}\sin\frac{1}{2}\pi\alpha(t)$$
 (16)

in the  $|t| \ll \mu(n)$  Regge limit, implying the Regge zeros,

$$\alpha(t) + 2r = 0, (17)$$

for integer r, again with the caveat that any real roots must satisfy t < 0 and  $t \neq \mu(n')$  for any integer  $n' \geq 0$ .

Equation (16) is nearly identical to Eq. (5), so applying the same logic as in the planar case, we deduce that  $\alpha(t)$  is a polynomial. As before, assumption ii) says that the residue is composed entirely of Regge zeros, so

$$R(n,t) = c(n) \prod_{r=1}^{N(n)} (\alpha(t) + 2r)^{h}.$$
 (18)

Since  $R(n,t) \sim R_{\rm str}(N(n),\alpha(t)/2)^h \sim N(n)^{h\alpha(t)/2}$  and the asymptotic behavior is  $R(n,t) \sim \mu(n)^{\alpha(t)}$ , we find that  $\mu(n) \sim N(n)^{h/2} \sim J(n)^{h/2}$ . Together with Eq. (2), this implies that  $\alpha(t) \sim t^{2/h}$ , and since  $\alpha(t)$  is a polynomial, we deduce that h=2. Thus the nonplanar amplitude also has the linear spectrum and Regge trajectory defined in Eq. (13). As in the planar case, App. A implies  $J(n) \sim n$ , so our nonplanar uniqueness arguments lead to the level truncation ansatz of Ref. [92], whose only solution with ultrasoft Regge behavior is  $\Gamma(-s)\Gamma(-t)\Gamma(-u)/\Gamma(1+s)\Gamma(1+t)\Gamma(1+u)$ , which is the Virasoro-Shapiro amplitude [115, 116].

Beyond Four-Point.—The Regge behavior of fivepoint string amplitudes also mandates Regge zeros in their corresponding residues. In fact, these zeros are sufficient to uniquely bootstrap the five-point amplitude.

Here we study a worldsheet integral representation of the five-point string amplitude constructed in Ref. [93],

$$A_{\text{str},5} = \int_{0}^{1} \frac{du_{13}}{u_{13}(1-u_{13})} \frac{du_{14}}{u_{14}(1-u_{14})} F(u_{13}, u_{14}),$$

$$F(u_{13}, u_{14}) = \frac{u_{13}^{X_{13}} u_{14}^{X_{14}} (1-u_{13})^{X_{24}} (1-u_{14})^{X_{35}}}{(1-u_{13}u_{14})^{X_{24}+X_{35}-X_{25}}},$$
(19)

where  $X_{ij} = (p_i + \cdots + p_{j-1})^2$  [117]. In the Regge limit  $X_{13}, X_{14} \to \infty$ , the factors of  $u_{13}^{X_{13}}$  and  $u_{14}^{X_{14}}$  localize the integral to the saddle at  $u_{13} = 1 - X_{13}^{-1}$  and  $u_{14} = 1 - X_{14}^{-1}$ , with scaling  $F(1 - X_{13}^{-1}, 1 - X_{14}^{-1}) \sim X_{13}^{X_{35} - X_{25}} X_{14}^{X_{24} - X_{25}} (X_{13} + X_{14})^{X_{25} - X_{24} - X_{35}}$ . Since we have  $F(1 - \epsilon, 1 - \epsilon) \sim \epsilon^{X_{25}}$ , the integral converges provided  $X_{25} > 0$ , just like the  $_3F_2$  closed-form expression for the five-point string amplitude [118], which can be analytically extended to the entire kinematic space using the formulas of Ref. [93]. In terms of the large s-like invariants,  $s_1 = -X_{13}$  and  $s_2 = -X_{14}$ , and the fixed t-like invariants,  $t_{23} = -X_{24}$ ,  $t_{34} = -X_{35}$ , and  $t_{51} = -X_{25}$ , we obtain our final result for the double Regge limit of the five-point string amplitude when  $s_1, s_2 \to \infty$ ,

<span id="page-3-1"></span>
$$A_{\text{str},5} \overset{s_{1,2} \to \infty}{\sim} f(a_1, a_2, a_{12}) (-s_1)^{a_1} (-s_2)^{a_2} (-s_1 - s_2)^{a_{12}}, (20)$$

where we have defined  $a_1=t_{51}-t_{34},\,a_2=t_{51}-t_{23},$  and  $a_{12}=t_{23}+t_{34}-t_{51}.$  In analogy with our analysis at four-point, we now compute the double big arc integral,  $\oint^{s_{1*}} \oint^{s_{2*}} ds_1 ds_2 \, A_{\rm str,5}$ , which, by deforming both arc integral contours into a double sum over residues, equals  $\sum_{n_1}^{n_{1*}} \sum_{n_2}^{n_{2*}} R_{\rm str,5}(n_1,n_2,t_{23},t_{34},t_{51})$ . Inserting Eq. (20), we compute the double big arc integral to be f times

$$\sum_{b=0}^{a_{12}} s_{1*}^{a_{1}+1+b} s_{2*}^{a_{2}+1+a_{12}-b} {a_{12} \choose b} \frac{\sin \pi a_{1} \sin \pi (a_{2}+a_{12})}{(a_{1}+1+b)(a_{2}+1+a_{12}-b)}, (21)$$

where have taken  $a_{12} \geq 0$  to be a nonnegative integer so that the binomial expansion is finite. Differentiating with respect to  $n_{1*}$  and  $n_{2*}$ , we obtain

<span id="page-3-2"></span>
$$\frac{R_{\text{str},5}(n_1, n_2, t_{23}, t_{34}, t_{51})}{\mu'(n_1)\mu'(n_2)\mu(n_1)^{a_1}\mu(n_2)^{a_2+a_{12}}} \sim f \sum_{b=0}^{a_{12}} \left(\frac{\mu(n_1)}{\mu(n_2)}\right)^b \binom{a_{12}}{b} \sin \pi a_1 \sin \pi(a_2+a_{12}), \tag{22}$$

which stringently constrains the asymptotics and zeros of the five-point residues. By inspection, from Eq. (22) we see that the residue vanishes when  $a_1$ ,  $a_2$ , and  $a_{12}$  are integers. Recalling that  $X_{25} > 0$  for convergence, we deduce the locus of the Regge zeros for five-point residues at large level,

<span id="page-3-3"></span>
$$-(a_1 + a_2 + a_{12}) \ge 1$$
 and  $a_{12} \ge 0$ , (23)

for integers  $a_1, a_2, a_{12}$ .

It is straightforward to verify the presence of these Regge zeros in the known expression for the five-point string residue  $R_{\text{str},5}(n_1, n_2, t_{23}, t_{34}, t_{51})$  [93],

<span id="page-4-2"></span>
$$\frac{(1+t_{23})_{n_1}(1+t_{34})_{n_2}}{n_1!n_2!}{}_3F_2\left[\begin{smallmatrix} -n_1, -n_2, -t_{23}-t_{34}+t_{51} \\ -n_1-t_{23}, -n_2-t_{34} \end{smallmatrix}; 1\right]. (24)$$

We see that there are indeed zeros for  $a_1+r_1=a_2+r_2=a_{12}+r_{12}-r_1-r_2=0$ , where  $|r_1|\leq n_1,\ |r_2|\leq n_2$ , and  $1\leq r_{12}\leq r_1+r_2$ . This agrees exactly with Eq. (23) at large  $n_1$  and  $n_2$ . These Regge zeros can be used to constrain a five-point residue ansatz  $R(n_1,n_2,t_{23},t_{34},t_{51})=\sum_{i=0}^{n_1}\sum_{j=0}^{n_2}\sum_{k=0}^{\min(n_1-i,n_2-j)}\lambda_{i,j,k}(n_1,n_2)t_{23}^it_{34}^jt_{51}^k$  that encodes all possible exchanges of particles with maximum spins  $n_{1,2}$  on the poles at  $s_{1,2}=n_{1,2}$ . Enforcing the Regge zeros fixes all  $\lambda_{i,j,k}(n_1,n_2)$  to values corresponding precisely to the known residues of the string in Eq. (24), modulo a normalization  $c(n_1,n_2)$ .

Since Regge zeros fix the five-point string residue, it is natural to ask about the structure of its zeros more generally, which we describe in detail in App. B. Setting all three t-like Mandelstams to certain negative integers, one finds that all but a finite number of the residues vanish, reducing an infinite dual resonant sum over all levels  $(n_1, n_2)$  down to a finite subset, making the amplitude rational in  $(s_1, s_2)$ . As shown in Ref. [91], four-point crossing symmetry sculpts out an analytically solvable space of amplitudes that exhibit level truncation. Analogously, in five-point scattering the full cyclic invariance of the external legs produces a similarly solvable algebraic problem. As detailed in App. C, demanding cyclic invariance in fact fixes all  $c(n_1, n_2) = 1$ . This implies that our bootstrapped amplitude is uniquely fixed by its residue zeros to match the five-point string amplitude. The q-deformed version of our planar analyses is presented in App. D.

Discussion.—We have demonstrated a universal link between the Regge behavior of an amplitude and the zeros of its residues. This result applies to any fourpoint, tree-level, crossing-symmetric, Lorentz invariant, unitary amplitude. Amplitudes with softer Regge behavior are constrained by more Regge zeros. Assuming a minimality condition that every zero of a residue is a Regge zero mandated by consistency, we proved that any ultrasoft amplitude has a stringy spectrum of masses squared that is linear in spin,  $\mu(n) \sim J(n) \sim n$ , together with a linear Regge trajectory,  $\alpha(t) \sim t$ . This reproduces the level truncation ansatz of Ref. [91], from which the exact mathematical formula for the Veneziano amplitude is derived uniquely. While extra zeros are certainly allowed, the notion of minimal zeros gives an operational definition of the "simplest consistent amplitude," which is the Veneziano amplitude itself. Last but not least, we established analogous Regge zero and uniqueness properties beyond the planar four-point limit, successfully bootstrapping the four-point closed-string and five-point open-string amplitudes.

Acknowledgments: We thank Justin Berman, Simon Caron-Huot, Henriette Elvang, Aaron Hillman, and Sasha Zhiboedov useful discussions. C.C. and M.T. are supported by the Department of Energy (Grant No. DE-SC0011632) and by the Walter Burke Institute for Theoretical Physics. G.N.R. is supported by the James Arthur Postdoctoral Fellowship at New York University. F.S. is supported by the research grants 2021-SGR-00649, PID2023-146686NB-C31, and funding from the European Union NextGenerationEU(PRTR-C17.I1). F.S. also thanks the Walter Burke Institute for Theoretical Physics for its hospitality during the completion of this work.

<span id="page-4-1"></span><span id="page-4-0"></span>A. Adams, N. Arkani-Hamed, S. Dubovsky, A. Nicolis, and R. Rattazzi, "Causality, analyticity and an IR obstruction to UV completion," *JHEP* 10 (2006) 014, arXiv:hep-th/0602178.

<sup>[2]</sup> A. Nicolis, R. Rattazzi, and E. Trincherini, "Energy's and amplitudes' positivity," *JHEP* **05** (2010) 095, arXiv:0912.4258 [hep-th]. [Erratum: *JHEP* **11** (2011) 128].

<sup>[3]</sup> B. Bellazzini, L. Martucci, and R. Torre, "Symmetries, sum rules and constraints on effective field theories," *JHEP* **09** (2014) 100, arXiv:1405.2960 [hep-th].

<sup>[4]</sup> B. Bellazzini, C. Cheung, and G. N. Remmen, "Quantum gravity constraints from unitarity and analyticity," Phys. Rev. D 93 (2016) 064076, arXiv:1509.00851 [hep-th].

<sup>[5]</sup> B. Bellazzini, F. Riva, J. Serra, and F. Sgarlata, "Massive higher spins: effective theory and consistency," JHEP 10 (2019) 189, arXiv:1903.08664 [hep-th].

<sup>[6]</sup> X. O. Camanho, J. D. Edelstein, J. Maldacena, and A. Zhiboedov, "Causality constraints on corrections to

the graviton three-point coupling," *JHEP* **02** (2016) 020, arXiv:1407.5597 [hep-th].

<sup>[7]</sup> N. Arkani-Hamed, T.-C. Huang, and Y.-t. Huang, "The EFT-hedron," *JHEP* 05 (2021) 259, arXiv:2012.15849 [hep-th].

<sup>[8]</sup> B. Bellazzini, J. Elias Miró, R. Rattazzi, M. Riembau, and F. Riva, "Positive moments for scattering amplitudes," Phys. Rev. D 104 (2021) 036006, arXiv:2011.00037 [hep-th].

<sup>[9]</sup> A. J. Tolley, Z.-Y. Wang, and S.-Y. Zhou, "New positivity bounds from full crossing symmetry," *JHEP* 05 (2021) 255, arXiv:2011.02400 [hep-th].

<span id="page-4-3"></span><sup>[10]</sup> Z. Bern, D. Kosmopoulos, and A. Zhiboedov, "Gravitational effective field theory islands, low-spin dominance, and the four-graviton amplitude," J. Phys. A 54 (2021) 344002, arXiv:2103.12728 [hep-th].

<sup>[11]</sup> L.-Y. Chiang, Y.-t. Huang, W. Li, L. Rodina, and H.-C. Weng, "Into the EFThedron and UV constraints from IR consistency," *JHEP* 03 (2022) 063, arXiv:2105.02862 [hep-th].

- [12] B. Bellazzini, M. Riembau, and F. Riva, "IR side of positivity bounds," *Phys. Rev. D* **106** [\(2022\) 105008,](http://dx.doi.org/10.1103/PhysRevD.106.105008) [arXiv:2112.12561 \[hep-th\]](http://arxiv.org/abs/2112.12561).
- [13] D. Karateev, Z. Komargodski, J. Penedones, and B. Sahoo, "Trace anomalies and the graviton-dilaton amplitude," *JHEP* **11** [\(2024\) 067,](http://dx.doi.org/10.1007/JHEP11(2024)067) [arXiv:2312.09308 \[hep-th\]](http://arxiv.org/abs/2312.09308).
- [14] C. Cheung and G. N. Remmen, "Multipositivity bounds for scattering amplitudes," *[Phys. Rev. D](http://dx.doi.org/10.1103/wt4x-2149)* **112** (2025) [016017,](http://dx.doi.org/10.1103/wt4x-2149) [arXiv:2505.05553 \[hep-th\]](http://arxiv.org/abs/2505.05553).
- [15] N. Arkani-Hamed, Y.-t. Huang, J.-Y. Liu, and G. N. Remmen, "Causality, unitarity, and the weak gravity conjecture," *JHEP* **03** [\(2022\) 083,](http://dx.doi.org/10.1007/JHEP03(2022)083) [arXiv:2109.13937](http://arxiv.org/abs/2109.13937) [\[hep-th\]](http://arxiv.org/abs/2109.13937).
- [16] C. Cheung and G. N. Remmen, "Infrared consistency and the weak gravity conjecture," *JHEP* **12** [\(2014\) 087,](http://dx.doi.org/10.1007/JHEP12(2014)087) [arXiv:1407.7865 \[hep-th\]](http://arxiv.org/abs/1407.7865).
- [17] G. N. Remmen and N. L. Rodd, "Consistency of the standard model effective field theory," *JHEP* **12** [\(2019\)](http://dx.doi.org/10.1007/JHEP12(2019)032) [032,](http://dx.doi.org/10.1007/JHEP12(2019)032) [arXiv:1908.09845 \[hep-ph\]](http://arxiv.org/abs/1908.09845).
- [18] G. N. Remmen and N. L. Rodd, "Flavor Constraints from Unitarity and Analyticity," *[Phys. Rev. Lett.](http://dx.doi.org/10.1103/PhysRevLett.127.149901)* **125** [\(2020\) 081601,](http://dx.doi.org/10.1103/PhysRevLett.127.149901) [arXiv:2004.02885 \[hep-ph\]](http://arxiv.org/abs/2004.02885). [Erratum: *[Phys. Rev. Lett.](https://doi.org/10.1103/PhysRevLett.127.149901)* **127**, 149901 (2021)].
- [19] G. N. Remmen and N. L. Rodd, "Signs, spin, SMEFT: Sum rules at dimension six," *[Phys. Rev. D](http://dx.doi.org/10.1103/PhysRevD.105.036006)* **105** (2022) [036006,](http://dx.doi.org/10.1103/PhysRevD.105.036006) [arXiv:2010.04723 \[hep-ph\]](http://arxiv.org/abs/2010.04723).
- [20] G. N. Remmen and N. L. Rodd, "Spinning sum rules for the dimension-six SMEFT," *JHEP* **09** [\(2022\) 030,](http://dx.doi.org/10.1007/JHEP09(2022)030) [arXiv:2206.13524 \[hep-ph\]](http://arxiv.org/abs/2206.13524).
- [21] G. N. Remmen and N. L. Rodd, "Positively Identifying HEFT or SMEFT," [arXiv:2412.07827 \[hep-ph\]](http://arxiv.org/abs/2412.07827).
- [22] R. Aoude, G. Elor, G. N. Remmen, and O. Sumensari, "Positivity in Amplitudes from Quantum Entanglement," [arXiv:2402.16956 \[hep-th\]](http://arxiv.org/abs/2402.16956).
- [23] C. Cheung and G. N. Remmen, "Positive signs in massive gravity," *JHEP* **04** [\(2016\) 002,](http://dx.doi.org/10.1007/JHEP04(2016)002) [arXiv:1601.04068](http://arxiv.org/abs/1601.04068) [\[hep-th\]](http://arxiv.org/abs/1601.04068).
- [24] C. Cheung and G. N. Remmen, "Positivity of Curvature-Squared Corrections in Gravity," *[Phys. Rev.](http://dx.doi.org/10.1103/PhysRevLett.118.051601) Lett.* **118** [\(2017\) 051601,](http://dx.doi.org/10.1103/PhysRevLett.118.051601) [arXiv:1608.02942 \[hep-th\]](http://arxiv.org/abs/1608.02942).
- [25] D. Green, Y. Huang, C.-H. Shen, and D. Baumann, "Positivity from cosmological correlators," *[JHEP](http://dx.doi.org/10.1007/JHEP04(2024)034)* **04** [\(2024\) 034,](http://dx.doi.org/10.1007/JHEP04(2024)034) [arXiv:2310.02490 \[hep-th\]](http://arxiv.org/abs/2310.02490).
- [26] M. Freytsis, S. Kumar, G. N. Remmen, and N. L. Rodd, "Multifield positivity bounds for inflation," *[JHEP](http://dx.doi.org/10.1007/JHEP09(2023)041)* **09** [\(2023\) 041,](http://dx.doi.org/10.1007/JHEP09(2023)041) [arXiv:2210.10791 \[hep-th\]](http://arxiv.org/abs/2210.10791).
- [27] B. Bellazzini, M. Lewandowski, and J. Serra, "Positivity of Amplitudes, Weak Gravity Conjecture, and Modified Gravity," *[Phys. Rev. Lett.](http://dx.doi.org/10.1103/PhysRevLett.123.251103)* **123** (2019) 251103, [arXiv:1902.03250 \[hep-th\]](http://arxiv.org/abs/1902.03250).
- [28] S. Andriolo, T.-C. Huang, T. Noumi, H. Ooguri, and G. Shiu, "Duality and axionic weak gravity," *[Phys. Rev.](http://dx.doi.org/10.1103/PhysRevD.102.046008) D* **102** [\(2020\) 046008,](http://dx.doi.org/10.1103/PhysRevD.102.046008) [arXiv:2004.13721 \[hep-th\]](http://arxiv.org/abs/2004.13721).
- [29] C. de Rham, S. Melville, and A. J. Tolley, "Improved positivity bounds and massive gravity," *JHEP* **04** [\(2018\)](http://dx.doi.org/10.1007/JHEP04(2018)083) [083,](http://dx.doi.org/10.1007/JHEP04(2018)083) [arXiv:1710.09611 \[hep-th\]](http://arxiv.org/abs/1710.09611).
- [30] V. Chandrasekaran, G. N. Remmen, and A. Shahbazi-Moghaddam, "Higher-point positivity," *[JHEP](http://dx.doi.org/10.1007/JHEP11(2018)015)* **11** [\(2018\) 015,](http://dx.doi.org/10.1007/JHEP11(2018)015) [arXiv:1804.03153 \[hep-th\]](http://arxiv.org/abs/1804.03153).
- [31] A. Jenkins and D. O'Connell, "The Story of O: Positivity constraints in effective field theories," [arXiv:hep](http://arxiv.org/abs/hep-th/0609159)[th/0609159 \[hep-th\]](http://arxiv.org/abs/hep-th/0609159).
- [32] G. Dvali, A. Franca, and C. Gomez, "Road Signs for

- UV-Completion," [arXiv:1204.6388 \[hep-th\]](http://arxiv.org/abs/1204.6388).
- [33] T. N. Pham and T. N. Truong, "Evaluation of the derivative quartic terms of the meson chiral Lagrangian from forward dispersion relations," *[Phys. Rev.](http://dx.doi.org/10.1103/PhysRevD.31.3027)* **D31** [\(1985\) 3027.](http://dx.doi.org/10.1103/PhysRevD.31.3027)
- [34] B. Ananthanarayan, D. Toublan, and G. Wanders, "Consistency of the chiral pion-pion scattering amplitudes with axiomatic constraints," *[Phys. Rev.](http://dx.doi.org/10.1103/PhysRevD.51.1093)* **D51** [\(1995\) 1093,](http://dx.doi.org/10.1103/PhysRevD.51.1093) [arXiv:hep-ph/9410302 \[hep-ph\]](http://arxiv.org/abs/hep-ph/9410302).
- <span id="page-5-0"></span>[35] M. R. Pennington and J. Portoles, "The chiral lagrangian parameters, *ℓ*1, *ℓ*2, are determined by the *ρ*resonance," *[Phys. Lett.](http://dx.doi.org/10.1016/0370-2693(94)01551-M)* **B344** (1995) 399, [arXiv:hep](http://arxiv.org/abs/hep-ph/9409426)[ph/9409426 \[hep-ph\]](http://arxiv.org/abs/hep-ph/9409426).
- <span id="page-5-1"></span>[36] M. F. Paulos, J. Penedones, J. Toledo, B. C. van Rees, and P. Vieira, "The S-matrix bootstrap. Part III: higher dimensional amplitudes," *JHEP* **12** [\(2019\) 040,](http://dx.doi.org/10.1007/JHEP12(2019)040) [arXiv:1708.06765 \[hep-th\]](http://arxiv.org/abs/1708.06765).
- [37] L. Córdova and P. Vieira, "Adding flavour to the S-matrix bootstrap," *JHEP* **12** [\(2018\) 063,](http://dx.doi.org/10.1007/JHEP12(2018)063) [arXiv:1805.11143 \[hep-th\]](http://arxiv.org/abs/1805.11143).
- [38] A. L. Guerrieri, J. Penedones, and P. Vieira, "Bootstrapping QCD Using Pion Scattering Amplitudes," *[Phys. Rev. Lett.](http://dx.doi.org/10.1103/PhysRevLett.122.241604)* **122** (2019) 241604, [arXiv:1810.12849](http://arxiv.org/abs/1810.12849) [\[hep-th\]](http://arxiv.org/abs/1810.12849).
- [39] J. Elias Miró, A. L. Guerrieri, A. Hebbar, J. Penedones, and P. Vieira, "Flux Tube S-matrix Bootstrap," *[Phys.](http://dx.doi.org/10.1103/PhysRevLett.123.221602) Rev. Lett.* **123** [\(2019\) 221602,](http://dx.doi.org/10.1103/PhysRevLett.123.221602) [arXiv:1906.08098 \[hep](http://arxiv.org/abs/1906.08098)[th\]](http://arxiv.org/abs/1906.08098).
- [40] L. Córdova, Y. He, M. Kruczenski, and P. Vieira, "The O(N) S-matrix monolith," *JHEP* **04** [\(2020\) 142,](http://dx.doi.org/10.1007/JHEP04(2020)142) [arXiv:1909.06495 \[hep-th\]](http://arxiv.org/abs/1909.06495).
- [41] D. Karateev, S. Kuhn, and J. Penedones, "Bootstrapping massive quantum field theories," *JHEP* **07** [\(2020\)](http://dx.doi.org/10.1007/JHEP07(2020)035) [035,](http://dx.doi.org/10.1007/JHEP07(2020)035) [arXiv:1912.08940 \[hep-th\]](http://arxiv.org/abs/1912.08940).
- [42] A. L. Guerrieri, A. Homrich, and P. Vieira, "Dual Smatrix bootstrap. Part I. 2D theory," *JHEP* **11** [\(2020\)](http://dx.doi.org/10.1007/JHEP11(2020)084) [084,](http://dx.doi.org/10.1007/JHEP11(2020)084) [arXiv:2008.02770 \[hep-th\]](http://arxiv.org/abs/2008.02770).
- [43] A. L. Guerrieri, J. Penedones, and P. Vieira, "S-matrix bootstrap for effective field theories: massless pions," *JHEP* **06** [\(2021\) 088,](http://dx.doi.org/10.1007/JHEP06(2021)088) [arXiv:2011.02802 \[hep-th\]](http://arxiv.org/abs/2011.02802).
- [44] M. Carrillo Gonzalez, C. de Rham, V. Pozsgay, and A. J. Tolley, "Causal effective field theories," *[Phys. Rev.](http://dx.doi.org/10.1103/PhysRevD.106.105018) D* **106** [\(2022\) 105018,](http://dx.doi.org/10.1103/PhysRevD.106.105018) [arXiv:2207.03491 \[hep-th\]](http://arxiv.org/abs/2207.03491).
- [45] Y. He and M. Kruczenski, "S-matrix bootstrap in 3+1 dimensions: regularization and dual convex problem," *JHEP* **08** [\(2021\) 125,](http://dx.doi.org/10.1007/JHEP08(2021)125) [arXiv:2103.11484 \[hep-th\]](http://arxiv.org/abs/2103.11484).
- [46] J. Elias Miró and A. Guerrieri, "Dual EFT bootstrap: QCD flux tubes," *JHEP* **10** [\(2021\) 126,](http://dx.doi.org/10.1007/JHEP10(2021)126) [arXiv:2106.07957 \[hep-th\]](http://arxiv.org/abs/2106.07957).
- [47] A. Guerrieri and A. Sever, "Rigorous Bounds on the Analytic *S* Matrix," *[Phys. Rev. Lett.](http://dx.doi.org/10.1103/PhysRevLett.127.251601)* **127** (2021) 251601, [arXiv:2106.10257 \[hep-th\]](http://arxiv.org/abs/2106.10257).
- [48] J. Elias Miro, A. Guerrieri, and M. A. Gumus, "Bridging positivity and S-matrix bootstrap bounds," *[JHEP](http://dx.doi.org/10.1007/JHEP05(2023)001)* **05** [\(2023\) 001,](http://dx.doi.org/10.1007/JHEP05(2023)001) [arXiv:2210.01502 \[hep-th\]](http://arxiv.org/abs/2210.01502).
- [49] M. Correia, A. Georgoudis, and A. L. Guerrieri, "Cross-Section Bootstrap: Unveiling the Froissart Amplitude," [arXiv:2506.04313 \[hep-th\]](http://arxiv.org/abs/2506.04313).
- <span id="page-5-2"></span>[50] C. de Rham, A. J. Tolley, Z.-H. Wang, and S.-Y. Zhou, "Primal S-matrix bootstrap with dispersion relations," [arXiv:2506.22546 \[hep-th\]](http://arxiv.org/abs/2506.22546).
- <span id="page-5-3"></span>[51] D. Simmons-Duffin, "A semidefinite program solver for the conformal bootstrap," *JHEP* **06** [\(2015\) 174,](http://dx.doi.org/10.1007/JHEP06(2015)174)

- [arXiv:1502.02033 \[hep-th\]](http://arxiv.org/abs/1502.02033).
- <span id="page-6-7"></span>[52] S. Caron-Huot and V. Van Duong, "Extremal effective field theories," *JHEP* **05** [\(2021\) 280,](http://dx.doi.org/10.1007/JHEP05(2021)280) [arXiv:2011.02957](http://arxiv.org/abs/2011.02957) [\[hep-th\]](http://arxiv.org/abs/2011.02957).
- [53] S. Caron-Huot, D. Mazac, L. Rastelli, and D. Simmons-Duffin, "Sharp boundaries for the swampland," *[JHEP](http://dx.doi.org/10.1007/JHEP07(2021)110)* **07** [\(2021\) 110,](http://dx.doi.org/10.1007/JHEP07(2021)110) [arXiv:2102.08951 \[hep-th\]](http://arxiv.org/abs/2102.08951).
- [54] S. Caron-Huot, Y.-Z. Li, J. Parra-Martinez, and D. Simmons-Duffin, "Causality constraints on corrections to Einstein gravity," *JHEP* **05** [\(2023\) 122,](http://dx.doi.org/10.1007/JHEP05(2023)122) [arXiv:2201.06602 \[hep-th\]](http://arxiv.org/abs/2201.06602).
- [55] J. Albert and L. Rastelli, "Bootstrapping pions at large *N*," *JHEP* **08** [\(2022\) 151,](http://dx.doi.org/10.1007/JHEP08(2022)151) [arXiv:2203.11950 \[hep-th\]](http://arxiv.org/abs/2203.11950).
- [56] K. Häring, A. Hebbar, D. Karateev, M. Meineri, and J. Penedones, "Bounds on photon scattering," *[JHEP](http://dx.doi.org/10.1007/JHEP10(2024)103)* **10** [\(2024\) 103,](http://dx.doi.org/10.1007/JHEP10(2024)103) [arXiv:2211.05795 \[hep-th\]](http://arxiv.org/abs/2211.05795).
- [57] C. Fernandez, A. Pomarol, F. Riva, and F. Sciotti, "Cornering large-*N<sup>c</sup>* QCD with positivity bounds," *[JHEP](http://dx.doi.org/10.1007/JHEP06(2023)094)* **06** [\(2023\) 094,](http://dx.doi.org/10.1007/JHEP06(2023)094) [arXiv:2211.12488 \[hep-th\]](http://arxiv.org/abs/2211.12488).
- [58] J. Albert and L. Rastelli, "Bootstrapping pions at large *N*. Part II. Background gauge fields and the chiral anomaly," *JHEP* **09** [\(2024\) 039,](http://dx.doi.org/10.1007/JHEP09(2024)039) [arXiv:2307.01246 \[hep](http://arxiv.org/abs/2307.01246)[th\]](http://arxiv.org/abs/2307.01246).
- [59] J. Albert, J. Henriksson, L. Rastelli, and A. Vichi, "Bootstrapping mesons at large *N*: Regge trajectory from spin-two maximization," *JHEP* **09** [\(2024\) 172,](http://dx.doi.org/10.1007/JHEP09(2024)172) [arXiv:2312.15013 \[hep-th\]](http://arxiv.org/abs/2312.15013).
- [60] Y.-Z. Li, "Effective field theory bootstrap, large-*N χ*PT and holographic QCD," *JHEP* **01** [\(2024\) 072,](http://dx.doi.org/10.1007/JHEP01(2024)072) [arXiv:2310.09698 \[hep-th\]](http://arxiv.org/abs/2310.09698).
- [61] T. Ma, A. Pomarol, and F. Sciotti, "Bootstrapping the chiral anomaly at large *Nc*," *JHEP* **11** [\(2023\) 176,](http://dx.doi.org/10.1007/JHEP11(2023)176) [arXiv:2307.04729 \[hep-th\]](http://arxiv.org/abs/2307.04729).
- [62] J. Berman, H. Elvang, and A. Herderschee, "Flattening of the EFT-hedron: supersymmetric positivity bounds and the search for string theory," *JHEP* **03** [\(2024\) 021,](http://dx.doi.org/10.1007/JHEP03(2024)021) [arXiv:2310.10729 \[hep-th\]](http://arxiv.org/abs/2310.10729).
- [63] C. Beadle, G. Isabella, D. Perrone, S. Ricossa, F. Riva, and F. Serra, "Non-Forward UV/IR Relations," [arXiv:2407.02346 \[hep-th\]](http://arxiv.org/abs/2407.02346).
- [64] Z.-Y. Dong, T. Ma, A. Pomarol, and F. Sciotti, "Bootstrapping the chiral-gravitational anomaly," *[JHEP](http://dx.doi.org/10.1007/JHEP05(2025)114)* **05** [\(2025\) 114,](http://dx.doi.org/10.1007/JHEP05(2025)114) [arXiv:2411.14422 \[hep-th\]](http://arxiv.org/abs/2411.14422).
- [65] J. Berman, H. Elvang, N. Geiser, and L. L. Lin, "Bootstrapping Extremal Scalar Amplitudes With and Without Supersymmetry," [arXiv:2412.13368 \[hep-th\]](http://arxiv.org/abs/2412.13368).
- [66] J. Berman, "Analytic Bounds on the Spectrum of Crossing Symmetric S-Matrices," [arXiv:2410.01914 \[hep-th\]](http://arxiv.org/abs/2410.01914).
- [67] C.-H. Chang and J. Parra-Martinez, "Graviton loops and negativity," [arXiv:2501.17949 \[hep-th\]](http://arxiv.org/abs/2501.17949).
- [68] C. Beadle, G. Isabella, D. Perrone, S. Ricossa, F. Riva, and F. Serra, "The EFT bootstrap at finite *MP L*," *JHEP* **06** [\(2025\) 209,](http://dx.doi.org/10.1007/JHEP06(2025)209) [arXiv:2501.18465 \[hep-th\]](http://arxiv.org/abs/2501.18465).
- [69] J. Berman, H. Elvang, and C. Figueiredo, "Splitting Regions and Shrinking Islands from Higher Point Constraints," [arXiv:2506.22538 \[hep-th\]](http://arxiv.org/abs/2506.22538).
- <span id="page-6-0"></span>[70] B. Bellazzini, A. Pomarol, M. Romano, and F. Sciotti, "(Super) Gravity from Positivity," [arXiv:2507.12535](http://arxiv.org/abs/2507.12535) [\[hep-th\]](http://arxiv.org/abs/2507.12535).
- <span id="page-6-1"></span>[71] Y.-t. Huang, J.-Y. Liu, L. Rodina, and Y. Wang, "Carving out the space of open-string S-matrix," *[JHEP](http://dx.doi.org/10.1007/JHEP04(2021)195)* **04** [\(2021\) 195,](http://dx.doi.org/10.1007/JHEP04(2021)195) [arXiv:2008.02293 \[hep-th\]](http://arxiv.org/abs/2008.02293).
- [72] A. Guerrieri, J. Penedones, and P. Vieira, "Where Is

- String Theory in the Space of Scattering Amplitudes?," *[Phys. Rev. Lett.](http://dx.doi.org/10.1103/PhysRevLett.127.081601)* **127** (2021) 081601, [arXiv:2102.02847](http://arxiv.org/abs/2102.02847) [\[hep-th\]](http://arxiv.org/abs/2102.02847).
- [73] C. Cheung and G. N. Remmen, "Veneziano variations: how unique are string amplitudes?," *JHEP* **01** [\(2023\)](http://dx.doi.org/10.1007/JHEP01(2023)122) [122,](http://dx.doi.org/10.1007/JHEP01(2023)122) [arXiv:2210.12163 \[hep-th\]](http://arxiv.org/abs/2210.12163).
- [74] N. Geiser and L. W. Lindwasser, "Properties of infinite product amplitudes: Veneziano, Virasoro, and Coon," *JHEP* **12** [\(2022\) 112,](http://dx.doi.org/10.1007/JHEP12(2022)112) [arXiv:2207.08855 \[hep-th\]](http://arxiv.org/abs/2207.08855).
- [75] N. Geiser and L. W. Lindwasser, "Generalized Veneziano and Virasoro amplitudes," *JHEP* **04** [\(2023\)](http://dx.doi.org/10.1007/JHEP04(2023)031) [031,](http://dx.doi.org/10.1007/JHEP04(2023)031) [arXiv:2210.14920 \[hep-th\]](http://arxiv.org/abs/2210.14920).
- <span id="page-6-8"></span>[76] Y.-t. Huang and G. N. Remmen, "UV-complete gravity amplitudes and the triple product," *[Phys. Rev. D](http://dx.doi.org/10.1103/PhysRevD.106.L021902)* **106** [\(2022\) L021902,](http://dx.doi.org/10.1103/PhysRevD.106.L021902) [arXiv:2203.00696 \[hep-th\]](http://arxiv.org/abs/2203.00696).
- <span id="page-6-11"></span>[77] J. Maldacena and G. N. Remmen, "Accumulation-point amplitudes in string theory," *JHEP* **08** [\(2022\) 152,](http://dx.doi.org/10.1007/JHEP08(2022)152) [arXiv:2207.06426 \[hep-th\]](http://arxiv.org/abs/2207.06426).
- <span id="page-6-9"></span>[78] C. Cheung and G. N. Remmen, "Stringy dynamics from an amplitudes bootstrap," *[Phys. Rev. D](http://dx.doi.org/10.1103/PhysRevD.108.026011)* **108** (2023) [026011,](http://dx.doi.org/10.1103/PhysRevD.108.026011) [arXiv:2302.12263 \[hep-th\]](http://arxiv.org/abs/2302.12263).
- <span id="page-6-4"></span>[79] C. Cheung and G. N. Remmen, "Bespoke dual resonance," *[Phys. Rev. D](http://dx.doi.org/10.1103/PhysRevD.108.086009)* **108** (2023) 086009, [arXiv:2308.03833 \[hep-th\]](http://arxiv.org/abs/2308.03833).
- <span id="page-6-12"></span>[80] N. Arkani-Hamed, C. Cheung, C. Figueiredo, and G. N. Remmen, "Multiparticle Factorization and the Rigidity of String Theory," *[Phys. Rev. Lett.](http://dx.doi.org/10.1103/PhysRevLett.132.091601)* **132** (2024) 091601, [arXiv:2312.07652 \[hep-th\]](http://arxiv.org/abs/2312.07652).
- [81] R. Bhardwaj, M. Spradlin, A. Volovich, and H.-C. Weng, "Unitarity of bespoke amplitudes," *[Phys. Rev.](http://dx.doi.org/10.1103/PhysRevD.110.106016) D* **110** [\(2024\) 106016,](http://dx.doi.org/10.1103/PhysRevD.110.106016) [arXiv:2406.04410 \[hep-th\]](http://arxiv.org/abs/2406.04410).
- [82] C. B. Jepsen, "Cutting the Coon amplitude," *[JHEP](http://dx.doi.org/10.1007/JHEP06(2023)114)* **06** [\(2023\) 114,](http://dx.doi.org/10.1007/JHEP06(2023)114) [arXiv:2303.02149 \[hep-th\]](http://arxiv.org/abs/2303.02149).
- [83] A. Gadde and S. Jain, "Analysis of *s*-*t* symmetric classical S-matrices," [arXiv:2502.18033 \[hep-th\]](http://arxiv.org/abs/2502.18033).
- [84] F. Figueroa and P. Tourkine, "Unitarity and Low Energy Expansion of the Coon Amplitude," *[Phys. Rev.](http://dx.doi.org/10.1103/PhysRevLett.129.121602) Lett.* **129** [\(2022\) 121602,](http://dx.doi.org/10.1103/PhysRevLett.129.121602) [arXiv:2201.12331 \[hep-th\]](http://arxiv.org/abs/2201.12331).
- [85] R. Bhardwaj, S. De, M. Spradlin, and A. Volovich, "On unitarity of the Coon amplitude," *JHEP* **08** [\(2023\) 082,](http://dx.doi.org/10.1007/JHEP08(2023)082) [arXiv:2212.00764 \[hep-th\]](http://arxiv.org/abs/2212.00764).
- <span id="page-6-13"></span>[86] N. Geiser, "The Baker-Coon-Romans *N*-point amplitude and an exact field theory limit of the Coon amplitude," *JHEP* **10** [\(2024\) 010,](http://dx.doi.org/10.1007/JHEP10(2024)010) [arXiv:2311.04130 \[hep-th\]](http://arxiv.org/abs/2311.04130).
- <span id="page-6-2"></span>[87] K. Häring and A. Zhiboedov, "The stringy S-matrix bootstrap: maximal spin and superpolynomial softness," *JHEP* **10** [\(2024\) 075,](http://dx.doi.org/10.1007/JHEP10(2024)075) [arXiv:2311.13631 \[hep-th\]](http://arxiv.org/abs/2311.13631).
- <span id="page-6-3"></span>[88] S. Caron-Huot, Z. Komargodski, A. Sever, and A. Zhiboedov, "Strings from massive higher spins: the asymptotic uniqueness of the Veneziano amplitude," *[JHEP](http://dx.doi.org/10.1007/JHEP10(2017)026)* **10** [\(2017\) 026,](http://dx.doi.org/10.1007/JHEP10(2017)026) [arXiv:1607.04253 \[hep-th\]](http://arxiv.org/abs/1607.04253).
- [89] A. Sever and A. Zhiboedov, "On fine structure of strings: the universal correction to the Veneziano amplitude," *JHEP* **06** [\(2018\) 054,](http://dx.doi.org/10.1007/JHEP06(2018)054) [arXiv:1707.05270 \[hep-th\]](http://arxiv.org/abs/1707.05270).
- <span id="page-6-10"></span>[90] J. Albert, W. Knop, and L. Rastelli, "Where is tree-level string theory?," *JHEP* **02** [\(2025\) 157,](http://dx.doi.org/10.1007/JHEP02(2025)157) [arXiv:2406.12959](http://arxiv.org/abs/2406.12959) [\[hep-th\]](http://arxiv.org/abs/2406.12959).
- <span id="page-6-5"></span>[91] C. Cheung, A. Hillman, and G. N. Remmen, "Bootstrap Principle for the Spectrum and Scattering of Strings," *[Phys. Rev. Lett.](http://dx.doi.org/10.1103/PhysRevLett.133.251601)* **133** (2024) 251601, [arXiv:2406.02665](http://arxiv.org/abs/2406.02665) [\[hep-th\]](http://arxiv.org/abs/2406.02665).
- <span id="page-6-6"></span>[92] C. Cheung, A. Hillman, and G. N. Remmen, "Uniqueness criteria for the Virasoro-Shapiro amplitude," *[Phys.](http://dx.doi.org/10.1103/PhysRevD.111.086034)*

- Rev. D 111 (2025) 086034, arXiv:2408.03362 [hep-th].
- <span id="page-7-22"></span>[93] N. Arkani-Hamed, C. Figueiredo, and G. N. Remmen, "Open string amplitudes: singularities, asymptotics and new representations," *JHEP* 04 (2025) 039, arXiv:2412.20639 [hep-th].
- <span id="page-7-0"></span>[94] C. B. Jepsen, "The Veneziano Amplitude in any Dimension and a Virasoro-Shapiro Partial Amplitude," arXiv:2506.05253 [hep-th].
- <span id="page-7-1"></span>[95] We work in mostly-plus metric signature throughout, where the Mandelstam invariants are  $s = -(p_1 + p_2)^2$ ,  $t = -(p_2 + p_3)^2$ , and  $u = -(p_1 + p_3)^2$  for cyclically labeled incoming momenta.
- <span id="page-7-2"></span>[96] Note that the t>0 behavior is unspecified, so a finite or infinite tower of spins is a priori allowed.
- <span id="page-7-3"></span>[97] Strictly speaking, we impose this property on the real component of  $\alpha(t)$ , which controls the magnitude of the amplitude in the Regge limit, along the real t axis.
- <span id="page-7-4"></span>[98] G. Veneziano, "Construction of a crossing-symmetric, Regge-behaved amplitude for linearly rising trajectories," *Nuovo Cim. A* **57** (1968) 190.
- <span id="page-7-5"></span>[99] D. D. Coon, "Uniqueness of the Veneziano representation," Phys. Lett. B 29 (1969) 669.
- <span id="page-7-6"></span>[100] F. A. Cerulus and A. Martin, "A lower bound for large angle elastic scattering at high energies," *Phys. Lett.* 8 (1964) 80.
- <span id="page-7-7"></span>[101] L. Buoninfante, J. Tokuda, and M. Yamaguchi, "New lower bounds on scattering amplitudes: non-locality constraints," *JHEP* 01 (2024) 082, arXiv:2305.16422 [hep-th].
- <span id="page-7-8"></span>[102] The case of infinite J(n) at finite n corresponds to an infinite tower of spins residing on an isolated, narrow resonance. This describes a propagating mode of unbounded spatial extent, in violation of locality [10, 52, 76].
- <span id="page-7-9"></span>[103] Equation (1) exhibits a branch cut that we take to be along the positive real s axis, which is the location of all physical discontinuities. To compute Eq. (3) we transform to polar coordinates,  $s = s_*e^{i\theta}$ , and integrate around the branch cut in the domain  $0 < \theta < 2\pi$ .
- <span id="page-7-10"></span>[104] Note that the hypergeometric [78] and  $_4F_3$  satellite [90] amplitudes exhibit power-law Regge behavior reminiscent of quantum field theory, where  $\alpha(t)$  crosses few if any integers, yielding few if any Regge zeros, as expected.
- <span id="page-7-11"></span>[105] D. J. Gross, "Factorization and the generalized Veneziano model with satellites," *Nucl. Phys. B* **13** (1969) 467.
- <span id="page-7-12"></span>[106] In the notation of Ref. [105], this deformed worldsheet integral corresponds to taking  $\tilde{f}_4(u) = g/4$ , a constant, and is hence the simplest possible model in that family. We can write the amplitude in closed form,

$$A(s,t) = \frac{\Gamma(-s)\Gamma(-t)}{\Gamma(-s-t)} {}_2F_2 \left[ \begin{array}{c} -s, -t \\ -\frac{s+t}{2}, \frac{1-s-t}{2} \end{array} ; \frac{g}{4} \right],$$

by expanding the exponential and resumming.

- <span id="page-7-13"></span>[107] C. Eckner, F. Figueroa, and P. Tourkine, "Regge bootstrap: From linear to nonlinear trajectories," *Phys. Rev.* D 111 (2025) 126005, arXiv:2401.08736 [hep-th].
- <span id="page-7-14"></span>[108] C. Eckner, F. Figueroa, and P. Tourkine, "On the number of Regge trajectories for dual amplitudes," arXiv:2405.21057 [hep-th].
- <span id="page-7-15"></span>[109] The original Coon amplitude [99], whose residues are consistent with unitarity only for q < 1, exhibits an accumulation point in the spectrum and is neither dual

- resonant nor meromorphic, on account of a branch cut dressing factor [78, 119]. While amplitudes of this form are realized in nonperturbative string theoretic constructions [77], they depart from our assumption of tree-level dynamics, so we leave them to future work.
- <span id="page-7-16"></span>[110] Here we have implicitly assumed that every zero of  $\alpha(t) + r = 0$  appears in R(n,t) as a block. Said another way, R(n,t) does not include partial branches of solutions where only a subset of the zeros of  $\alpha(t) + r = 0$  appear. This structure is justified because of the underlying logic of the Regge zeros: if one solution of  $\alpha(t) + r = 0$  appears in R(n,t), then n is sufficiently large to justify the arc integral argument that implies the Regge zeros, and thus the other solutions of  $\alpha(t) + r = 0$  should be exhibited as well.
- <span id="page-7-17"></span>[111] T. Regge, "Introduction to complex orbital momenta," Nuovo Cim. 14 (1959) 951.
- [112] G. F. Chew and S. C. Frautschi, "Principle of Equivalence for All Strongly Interacting Particles Within the S-Matrix Framework," Phys. Rev. Lett. 7 (1961) 394.
- <span id="page-7-18"></span>[113] G. F. Chew, S. C. Frautschi, and S. Mandelstam, "Regge poles in  $\pi-\pi$  scattering," *Phys. Rev.* **126** (1961) 1202.
- <span id="page-7-19"></span>[114] To be precise, the result of this bootstrap procedure is the Z-theory amplitude [120] for colored external scalars, which is trivially related to the amplitudes of the superstring and bosonic string by a simple prefactor and shift of the Mandelstam invariants, respectively. The latter correspond to offsets in  $\alpha(t)$ , J(n), or  $\mu(n)$ , which can also be used to generate so-called satellite amplitudes.
- <span id="page-7-20"></span>[115] M. A. Virasoro, "Alternative Constructions of Crossing-Symmetric Amplitudes with Regge Behavior," Phys. Rev. 177 (1969) 2309.
- <span id="page-7-21"></span>[116] J. A. Shapiro, "Electrostatic analogue for the Virasoro model," Phys. Lett. B 33 (1970) 361.
- <span id="page-7-23"></span>[117] This representation of the string amplitude is obtained by transforming the Koba-Nielsen integral into one involving the u-variables [80, 105, 121–127], with integrand going like  $u_{13}^{X_{13}}u_{14}^{X_{14}}u_{24}^{X_{25}}u_{35}^{X_{35}}$ . One then uses the u-equations that relate the different u-variables among each other to eliminate all but two of them [93].
- <span id="page-7-24"></span>[118] A. Białas and S. Pokorski, "High-energy behaviour of the Bardakçi-Ruegg amplitude," Nucl. Phys. B 10 (1969) 399.
- <span id="page-7-25"></span>[119] D. D. Coon, U. P. Sukhatme, and J. Tran Thanh Van, "Duality and proton-proton scattering at all angles," Phys. Lett. B 45 (1973) 287.
- <span id="page-7-26"></span>[120] J. J. M. Carrasco, C. R. Mafra, and O. Schlotterer, "Abelian Z-theory: NLSM amplitudes and  $\alpha'$ -corrections from the open string," *JHEP* **06** (2017) 093, arXiv:1608.02569 [hep-th].
- <span id="page-7-27"></span>[121] Z. Koba and H. B. Nielsen, "Reaction amplitude for n-mesons: A generalization of the Veneziano-Bardakçi-Ruegg-Virasoro model," Nucl. Phys. B 10 (1969) 633.
- [122] K. Bardakçi and H. Ruegg, "Reggeized resonance model for the production amplitude," Phys. Lett. B 28 (1968) 342.
- [123] H.-M. Chan and S. T. Tsou, "Explicit construction of the N-point function in the generalized Veneziano model," Phys. Lett. B 28 (1969) 485.
- [124] N. Arkani-Hamed, S. He, T. Lam, and H. Thomas, "Binary geometries, generalized particles and strings,

and cluster algebras," *Phys. Rev. D* 107 (2023) 066015, arXiv:1912.11764 [hep-th].

- [125] N. Arkani-Hamed, S. He, and T. Lam, "Stringy canonical forms," *JHEP* 02 (2021) 069, arXiv:1912.08707 [hep-th].
- [126] N. Arkani-Hamed, Q. Cao, J. Dong, C. Figueiredo, and S. He, "Hidden zeros for particle/string amplitudes and the unity of colored scalars, pions and gluons," *JHEP* 10 (2024) 231, arXiv:2312.16282 [hep-th].
- <span id="page-8-2"></span>[127] N. Arkani-Hamed and C. Figueiredo, "All-order splits and multi-soft limits for particle and string amplitudes," arXiv:2405.09608 [hep-th].
- <span id="page-8-5"></span>[128] L. Romans, "A new family of dual models ('q-strings)'." Preprint, 1988.

Appendix A: Linear Spin Spectrum.—Under our stated assumptions, we have shown that  $\mu(n) \sim J(n)$ . Here we will furthermore deduce that  $J(n) \sim n$ . To begin, we define  $t_*$  to be the locus of the Regge zero at  $\alpha(t_*) + r_* = 0$ , where  $r_* = J(n_*)$  is the maximum spin of the residue at level  $n_*$ . We compute the big arc integral,

$$\oint^{s_*} ds \, s^{r-1} A(s, t_*) \sim s_*^{\alpha(t_*) + r} = 0,$$
(A1)

which is zero for any r in the range  $1 \leq r < r_*$ . Hence the number of constraints scales as  $r_* \sim J(n_*)$  at large level. We can also rewrite the big arc integral as a sum of residues in the usual way,

$$\oint^{s_*} ds \, s^{r-1} A(s, t_*) \sim \sum_{n=0}^{n_* - 1} \mu(n)^{r-1} R(n, t_*). \tag{A2}$$

The residue at level n vanishes unless  $J(n) < J(n_*)$ , since the residues for which  $J(n) \ge J(n_*)$  carry vanishing factors of  $\alpha(t_*) + r_* = 0$ . Since  $\mu(n)$  is without loss of generality strictly monotonic, and we have proven that  $\mu(n) \sim J(n)$ , it follows that J(n) is strictly monotonic as well, and the sum truncates to the region  $n < n_*$ .

The above sum rules can be interpreted as a set of equations constraining the residues  $R(n, t_*)$ , which are taken to be variables. The total counting is then

# of equations = # of sum rules 
$$\sim J(n_*)$$
  
# of variables = # of nonzero residues  $\sim n_*$ . (A3)

Demanding that the number of constraints not exceed the number of variables, we have

<span id="page-8-3"></span>
$$J(n_*) \lesssim n_*$$
 (A4)

at large  $n_*$ . Since J(n) is monotonic and integer-valued,  $J(n+1)-J(n)\geq 1$ , which with Eq. (A4) leaves uniquely

$$J(n) \sim n$$
 (A5)

at large n. This establishes that the spectrum of spins grows linearly with the level.

![](_page_8_Figure_18.jpeg)

<span id="page-8-0"></span>Figure 3. The domain D(10,6,25), where each dot represents a nonzero residue at level  $(n_1,n_2)$  in the kinematic configuration  $(t_{23},t_{34},t_{51})=-(10,6,25)$ . The number of nonzero residues, and hence the dual resonant sum, is finite.

Appendix B: Five-Point Level Truncation.—Let us write the string amplitude as a dual resonant sum [93] of the five-point string residues in Eq. (24), so

<span id="page-8-4"></span>
$$A_{\text{str},5}(s_{1}, s_{2}, t_{23}, t_{34}, t_{51})$$

$$= \sum_{n_{1}=0}^{\infty} \sum_{n_{2}=0}^{\infty} \frac{R_{\text{str},5}(n_{1}, n_{2}, t_{23}, t_{34}, t_{51})}{(n_{1} - s_{1})(n_{2} - s_{2})}$$

$$= \frac{\Gamma(-s_{1})\Gamma(-s_{2})\Gamma(-t_{23})\Gamma(-t_{34})}{\Gamma(-s_{1} - t_{23})\Gamma(-s_{2} - t_{34})} \times \times {}_{3}F_{2} \begin{bmatrix} -s_{1}, -s_{2}, -t_{23} - t_{34} + t_{51} \\ -s_{1} - t_{23}, -s_{2} - t_{34} \end{bmatrix}.$$
(A6)

The second equality converges when  $t_{51} < 0$ , but can be extended to arbitrary kinematics via a  $_6F_5$  function [93]. The five-point residues of the string exhibit a generalization of the level truncation property of four-point amplitudes observed in Ref. [91]. In particular, on the kinematic configuration defined by  $(t_{23}, t_{34}, t_{51}) = -(k_1, k_2, k_3)$  for integers  $k_{1,2} \geq 1$  and  $k_3 \geq k_1 + k_2$ , all but a finite number of five-point string residues in Eq. (24) vanish. In particular, the only nonzero residues occur at level  $(n_1, n_2)$  in a compact domain  $D(k_1, k_2, k_3)$  defined by

<span id="page-8-1"></span>
$$\left\{ (n_1, n_2) \middle| \begin{array}{l} n_1 < k_3 - k_2, n_2 < k_3 - k_1, \\ -k_2 < n_1 - n_2 < k_1, \\ k_{1,2} \ge 1, k_3 \ge k_1 + k_2 \end{array} \right\},$$
(A7)

depicted in Fig. 3. On these kinematics, the amplitude in Eq. (A6) reduces to a rational polynomial in  $(s_1, s_2)$ .

**Appendix C: Five-Point Normalization.**—We saw earlier how the five-point Regge zeros are sufficient to uniquely fix the five-point residue of the string amplitude up to a level-dependent normalization factor  $c(n_1, n_2)$ , so

that the five-point bootstrapped amplitude is given by

$$A(s_1, s_2, t_{23}, t_{34}, t_{51}) = \sum_{n_1=0}^{\infty} \sum_{n_2=0}^{\infty} \frac{c(n_1, n_2) R_{\text{str}, 5}(n_1, n_2, t_{23}, t_{34}, t_{51})}{(n_1 - s_1)(n_2 - s_2)}.$$
 (A8)

Here we will fix these constants as well.

To fix  $c(n_1, n_2)$ , we exploit the invariance of the amplitude under cyclic permutations of the external legs. Level truncation renders the amplitude a rational function of the kinematics. On this configuration, the constraint of crossing is a solvable algebraic problem. In particular, setting  $(s_1, s_2, t_{23}, t_{34}, t_{51}) = -(k_1, k_4, k_2, k_3, k_5)$  for  $D(k_2, k_3, k_5)$  and  $D(k_3, k_4, k_1)$  defined in Eq. (A7), we enforce cyclic invariance of the dual resonant form of the five-point amplitude,  $A(-k_1, -k_4, -k_2, -k_3, -k_5) = A(-k_2, -k_5, -k_3, -k_4, -k_1)$ , that is,

<span id="page-9-0"></span>
$$\sum_{\substack{(n_1, n_2) \in \\ D(k_2, k_3, k_5)}} \frac{c(n_1, n_2) R_{\text{str}, 5}(n_1, n_2, -k_2, -k_3, -k_5)}{(n_1 + k_1)(n_2 + k_4)}$$

$$= \sum_{\substack{(n_1, n_2) \in \\ D(k_3, k_4, k_1)}} \frac{c(n_1, n_2) R_{\text{str}, 5}(n_1, n_2, -k_3, -k_4, -k_1)}{(n_1 + k_2)(n_2 + k_5)}.$$
(A9)

Enforcing this condition for all possible choices of integers yields a unique solution,  $c(n_1, n_2) = 1$ , modulo uniform global normalization. Hence the five-point string amplitude is uniquely specified by its zeros, together with cyclic symmetry.

Appendix D: q-Deformations.—The procedure of bootstrapping the four- and five-point amplitudes from their Regge zeros can be generalized from a linear spectrum to a q-spectrum. Relaxing bijectivity of  $\alpha(t)$  in assumption i) and running through our argument in text for the logarithmic Regge trajectory, we have Regge zeros at  $\alpha(t) + r = 0$  for integer r, whose root scales as  $t \propto 1 - q^{-r}$ . We choose the normalization of t such that these Regge zeros appear at  $t = (1 - q^{-r})/(1 - q) = [-r]_q$ , so the spectrum is given by the q-deformed integers,  $\mu(n) = [n]_q$ , and we have  $\alpha(t) = \log(1 + (q - 1)t)/\log q$ .

Let us first consider four-point scattering. Assuming that the residue exhibits only the Regge zeros,  $R(n,t) = c(n) \prod_{r=1}^{n} (t - [-r]_q)$ , we construct a dual resonant representation of the four-point amplitude, obtaining

<span id="page-9-1"></span>
$$A(s,t) = \sum_{n=0}^{\infty} \frac{R(n,t)}{\mu(n)-s} = q^{\alpha(s)\alpha(t)} \frac{\Gamma_q(-\alpha(s))\Gamma_q(-\alpha(t))}{\Gamma_q(-\alpha(s)-\alpha(t))}, (A10)$$

where the level truncation results of Ref. [91] have been used with crossing symmetry to uniquely fix  $c(n) = q^{n(n+3)/2}/\Gamma_q(1+n)$ . The above expression is precisely the Coon amplitude, which for q > 1 is dual resonant on kinematics where the sum in Eq. (A10) converges [91].

The five-point q-deformed amplitude can be written in

terms of the q-hypergeometric function  $_3\phi_2$  as [86, 128]

$$A_{q\text{-str},5}(s_1,\!s_2,\!t_{23},\!t_{34},\!t_{51}) = q^{\alpha(s_1)\alpha(t_{23}) + \alpha(s_2)\alpha(t_{34})} \times$$

$$\times \frac{\Gamma_{q}(-\alpha(s_{1}))\Gamma_{q}(-\alpha(s_{2}))\Gamma_{q}(-\alpha(t_{23}))\Gamma_{q}(-\alpha(t_{34}))}{\Gamma_{q}(-\alpha(s_{1})-\alpha(t_{23}))\Gamma_{q}(-\alpha(s_{2})-\alpha(t_{34}))} \times (A11)$$

$$\times_{3}\phi_{2} \begin{bmatrix} q^{-\alpha(s_{1})}, q^{-\alpha(s_{2})}, q^{\alpha(t_{51})-\alpha(t_{23})-\alpha(t_{34})} \\ q^{-\alpha(s_{1})-\alpha(t_{23})}, q^{-\alpha(s_{2})-\alpha(t_{34})} \end{bmatrix}; q; q$$

For q > 1, this amplitude is dual resonant for kinematics where the sum converges, equaling

$$\sum_{n_1=0}^{\infty} \sum_{n_2=0}^{\infty} \frac{R_{q\text{-str},5}(n_1, n_2, t_{23}, t_{34}, t_{51})}{([n_1]_q - s_1)([n_2]_q - s_2)}, \tag{A12}$$

for the corresponding five-point residue polynomial,

<span id="page-9-2"></span>
$$R_{q\text{-str},5}(n_{1}, n_{2}, t_{23}, t_{34}, t_{51}) = \frac{(-1)^{n_{1}+n_{2}}q^{\frac{n_{1}(n_{1}+3)+n_{2}(n_{2}+3)}{2}+n_{1}\alpha(t_{23})+n_{2}\alpha(t_{34})}}{\Gamma_{q}(n_{1}+1)\Gamma_{q}(n_{2}+1)} \times \frac{\Gamma_{q}(-\alpha(t_{23}))\Gamma_{q}(-\alpha(t_{34}))}{\Gamma_{q}(-n_{1}-\alpha(t_{23}))\Gamma_{q}(-n_{2}-\alpha(t_{34}))} \times \frac{\Gamma_{q}(-\alpha(t_{23}))\Gamma_{q}(-n_{2}-\alpha(t_{34}))}{q^{-n_{1}-\alpha(t_{23})},q^{-n_{2}-\alpha(t_{34})}}; q; q \right].$$
(A13)

Simply replacing  $t_{ij}$  with  $\alpha(t_{ij})$  in the definitions of  $a_1, a_2, a_{12}$ , so  $a_1 = \alpha(t_{51}) - \alpha(t_{34})$ ,  $a_2 = \alpha(t_{51}) - \alpha(t_{23})$ , and  $a_{12} = \alpha(t_{23}) + \alpha(t_{34}) - \alpha(t_{51})$ , we find that this residue exhibits precisely the same Regge zeros as the five-point string amplitude, namely,  $a_1 + r_1 = a_2 + r_2 = a_{12} + r_{12} - r_1 - r_2 = 0$ , that is,  $t_{23} = [r_2 - r_{12}]_q$ ,  $t_{34} = [r_1 - r_{12}]_q$ , and  $t_{51} = [-r_{12}]_q$ , all for integers  $|r_1| \leq n_1, |r_2| \leq n_2$ , and  $1 \leq r_{12} \leq r_1 + r_2$ .

Remarkably, we find that the q-deformed five-point amplitude can be constructed uniquely from these zeros, just like the string. Again taking the ansatz residue  $R(n_1,n_2,t_{23},t_{34},t_{51})=\sum_{i=0}^{n_1}\sum_{j=0}^{n_2}\sum_{k=0}^{\min(n_1-i,n_2-j)}\lambda_{i,j,k}(n_1,n_2)t_{23}^it_{34}^jt_{51}^k$  for states of spin up to  $n_{1,2}$  on the factorization channel at  $s_{1,2}=[n_{1,2}]_q$ , we find that the unique solution satisfying the zeros described above is the residue in Eq. (A13), up to some unfixed normalization  $c(n_1,n_2)$ .

Furthermore, on the special kinematic configuration  $(t_{23},t_{34},t_{51})=([-k_1]_q,[-k_2]_q,[-k_3]_q)$ , the only nonzero residues are precisely those defined in  $D(k_1,k_2,k_3)$  in Eq. (A7). That is, these q-deformed five-point residues exhibit level truncation, so as before we can algebraically fix the  $c(n_1,n_2)$  normalizations by demanding cyclic invariance. Here we again employ Eq. (A9), albeit with the residues in Eq. (A13) and with the denominators replaced by factors of  $[n]_q - [-k]_q$ . This yields a unique solution consistent with cyclic invariance for which  $c(n_1,n_2)=1$ , modulo global normalization. Thus, the q-deformed five-point amplitude is uniquely constructible from its Regge zeros and cyclic symmetry.