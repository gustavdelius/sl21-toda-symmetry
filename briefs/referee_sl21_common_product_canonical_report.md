# Canonical Gaussian product: the two Serre relations through mass order two

Date: 2026-10-24. Model: Codex (GPT-6). Reviewed commit: `745586e6a6f7bca5ed91a1e9cf89d6ef5df387f0`.

This is a fresh bounded referee calculation. Earlier completed outputs remain immutable. New work is confined to this report and `referee_scratch/sl21_common_product_20261024/canonical_serre/`. No subagents, book edits, git mutations, worktrees, commits, pushes or renders were used. No previous script was imported or rerun.

## Result

**Both quartics vanish through mass order two in the source-free canonical Gaussian product described below.** No action insertion is needed inside that product. The terms checked are the bare currents, every one-K2 tail, every two-K2 tail, every one-F_A tail, and the regulated order-lambda0 lambda2 improvement. The two marginal full-cluster scalar coefficients are derivatives separately for each tail placement; their cancellation does not depend on combining tails with action insertions.

For smooth compact tests f_A and f_B, set

\[
A=A(f_A),\qquad B=B(f_B),\qquad
c_s=1+s+s^{-1},\quad s=e^{-i\pi\rho},\quad \rho=2/t.
\]

Each repeated A has the same f_A, and each repeated B the same f_B. In particular, one common smooth charge test is allowed. Then the formal Gaussian products satisfy

\[
B^3A-c_sB^2AB+c_sBAB^2-AB^3=O(m^3),
\]
\[
A^3B-c_sA^2BA+c_sABA^2-BA^3=O(m^3).
\tag{1}
\]

The underlying product allows independent slot tests. Equation (1) does **not** assert the unpolarized quartic with four arbitrary independently assigned slot tests; a concrete failure of that stronger interpretation is given below.

There is also a coherent iterated Wick composition in this purely chiral charge sector through the stated mass order. Retaining constituent labels, Wick cross-contractions and cocycle factors obey the usual matching combinatorics under every bracketing. The complete boundary values exist in a common acyclic Gaussian family. Thus this route supplies more than four unrelated n-ary finite parts. Its scope is a formal algebra of smeared Gaussian Wick clusters, truncated in the couplings. This report does not construct interacting Hilbert-space charges or establish the complete algebra including all barred products. Matching a single interacting state/action map to this algebra is a separate task.

| Item | Verdict | Evidence and scope |
|---|---|---|
| Canonical source-free quartics through mass two | correct in the specified formal product | exact shuffle, individual tail primitives, complete face count |
| Action insertions required to remove the two marginal Serre poles | wrong for this canonical product | each relevant tail already has its own minority-coordinate primitive |
| Two-K2 antiholomorphic factor obstructs the Gaussian boundary product | wrong at this order | it is a polynomial, not an opposite-orientation singular factor |
| Regulated lambda0 lambda2 L improvement leaves a finite source-free quartic term | wrong in this continuation | vanishes initially and in every face of the physical patch |
| Same polynomial with arbitrary independently assigned B slot tests | wrong | one surviving middle-chamber coefficient is s^{-1}-s |
| Coherent iterated purely chiral Gaussian Wick composition through mass two | correct as a formal distributional construction | label/cocycle factorization, Wick matching, acyclic boundary values |
| Equality with the interacting product including all barred generators | unclear here | not inferred from the chiral calculation |

No theorem label or established physical symmetry algebra is recommended.

## 1. Product definition and exact Wick data

Use X(z)X(w) approximately -log(z-w),

\[
a_i^2=-1,\quad a_0a_2=1-t,\quad a_id=1,\quad d^2=-\rho,
\quad T=a_0+a_2=-td,\quad k=d+T.
\]

Put C=a0.partial X, beta=t-1, eta=t-2, gamma=t/(1-t). The charge densities are

\[
J_A=:e^{-dX}:,\qquad J_B=:Ce^{dX}:,
\]
\[
H_B=-t\lambda_2K_2+(t/h)\lambda_0\lambda_2\epsilon^{2h}L,
\quad H_A=-\gamma\lambda_1F_A,
\]
\[
K_2=:e^{(d+a_2)X+a_2\bar X}:,
\quad L=:Ce^{kX+T\bar X}:,
\quad F_A=:e^{-kX-T\bar X}:.
\tag{2}
\]

The charge densities are J+H. Lambda0 and lambda2 have mass order one, and lambda1 has mass order two. The auxiliary h appears here only to specify the cutoff improvement; there are no source-pair integrations in this canonical product. Restrict h=t after the initial cutoff limit and continuation. There is no finite independent composite-potential coupling.

A vertex is recorded as (q,m0,m2), with left momentum qd+m0 a0+m2 a2 and right momentum m0 a0+m2 a2:

\[
A=(-1,0,0),\ B=(1,0,0),\ K_2=(1,0,1),\
F_A=(-1,-1,-1),\ L=(1,1,1).
\]

B and L carry C. Every slot also retains its original current label A or B. Multiplication and fusion preserve the *multiset* of those labels. The net momentum alone is insufficient to reconstruct a cross-regulator after fusion.

For ordered slots i,j define

\[
S_{ij}=m_{0i}m_{0j}+m_{2i}m_{2j}
+(t-1)(m_{0i}m_{2j}+m_{2i}m_{0j}).
\]

The antiholomorphic exponent is S_ij and the holomorphic exponent is

\[
E_{ij}=\rho q_iq_j-q_i(m_{0j}+m_{2j})-q_j(m_{0i}+m_{2i})+S_{ij}
+(g-\rho)\mathbf1_{\text{same original current label}}.
\tag{3}
\]

The ordered cocycle is (-1)^(q_i(m0_j+m2_j)). A C at z_i contracts with another vertex with coefficient

\[
\frac{-q_j+m_{0j}+(t-1)m_{2j}}{z_i-z_j},
\]

and two C factors contract to (z_i-z_j)^(-2). The full Wick polynomial is the sum over all disjoint derivative pairings, with each unpaired C replaced by its regular normal field plus its contractions with every other vertex. This formula retains every descendant and every derivative-pair term. The numerical phase checks used nonconstant coherent normal fields satisfying a2 X=-t dX-a0 X in both chiralities; they did not set the regular fields to zero.

Give all chiral current coordinates strictly ordered imaginary heights in operator order. Use the conjugate heights for the right-moving factors. First take the boundary value in a common domain such as Re g>6, Re rho<-1, Re h>4. Extend the complete boundary products and their chamber restrictions by their common analytic continuation, and then put g=rho and h=t at generic regular points. No independent finite quartic term is supplied.

## 2. Full separated exchange and the exact shuffle

Directly from (3), reversal of a pair gives

\[
q_{AA}=q_{BB}=u=e^{-i\pi g},\qquad q_{AB}=q_{BA}=s^{-1}.
\tag{4}
\]

This includes every tail placement and its ordered cocycle. For example, the two K2 fields have holomorphic exponent g-1 and antiholomorphic exponent 1. The right factor is polynomial, and the integer shift plus cocycle leaves the same original-label exchange phase u. Neither a fermionic sign nor a second deformation parameter is introduced.

For the four species orders (BBBA,BBAB,BABB,ABBB), an independent enumeration of all labelled permutations gives

\[
M(u,s)(1,-c_s,c_s,-1)^T=(0,-D,D,0)^T,
\]
\[
D=\frac{(s-1)(s+1)(s-u)(s^2-u)(u+1)}{s^4}.
\tag{5}
\]

The same enumeration applies with A and B interchanged. It applies separately to the bare term and to the sums over placements of one K2, two K2, one F_A, and one L, since the coefficients in (2) are the same for every placement of a given kind. It is an identity of the full Wick kernels, not only their scalar pieces.

The scripts compare the actual complex powers, cocycles and full derivative Wick polynomials against (4) for every word, every spatial permutation and every placement. The three exponent settings include an initial-domain point, an off-diagonal point, and g=rho at t=sqrt(5). The shuffle identity (5) is also checked with symbolic u,s.

The zero of D at u=s proves (1) only if the corresponding middle-chamber integrals have no pole on g=rho. The next two sections supply that missing step without action insertions.

## 3. Marginal channels: tail-only primitives

### Bare currents

For B positions z1,z2,z3 and A at a, let

\[
f_i=-\sum_{j\ne i}(z_i-z_j)^{-1}+(z_i-a)^{-1}.
\]

The fully contracted scalar coefficient is

\[
f_1f_2f_3+\sum_{i<j}f_k(z_i-z_j)^{-2}=0.
\tag{6}
\]

The coefficient of one constant regular C also vanishes:

\[
\sum_{i<j}\{f_if_j+(z_i-z_j)^{-2}\}=0.
\tag{7}
\]

Both were rederived at symbolic minority position and two asymmetric triples. Thus the bare B^3A full degree is at least 2+3(g-rho). A^3B has only one derivative current and the same lower bound directly. Higher normal-field Taylor levels improve the degree further.

### One K2 in B^3A

Write

\[
R=\prod_{i<j}(z_i-z_j)^g\prod_i(z_i-a)^{-\rho},
\qquad S_i=\frac{(a-z_i)^2}{\prod_{j\ne i}(z_i-z_j)^2}.
\]

For a K2 tail in slot i, the remaining two C factors have contractions

\[
f_j=-\frac1{z_j-z_k}+\frac1{z_j-a}+\frac\eta{z_j-z_i},
\quad
f_k=-\frac1{z_k-z_j}+\frac1{z_k-a}+\frac\eta{z_k-z_i}.
\]

After consistently factoring the fixed ordered phase, their exact scalar vertex product obeys

\[
R\frac{a-z_i}{(z_j-z_i)(z_k-z_i)}
\left[f_jf_k+\frac1{(z_j-z_k)^2}\right]
=\frac\eta\rho\partial_a(RS_i).
\tag{8}
\]

The derivative-pair term and the contraction of each remaining B with the other remaining B are both retained. Multiply (8) by that tail's coefficient and cocycle. Each tail is separately an a derivative. There is no source integral and no cancellation with an action coefficient in this identity.

### One F_A in A^3B

For three A coordinates a_i and minority B at z put

\[
R_A=\prod_{i<j}(a_i-a_j)^g\prod_i(a_i-z)^{-\rho},
\quad S_i^A=\frac{(z-a_i)^2}{\prod_{j\ne i}(a_i-a_j)^2}.
\]

The scalar C contraction when F_A occupies slot i is

\[
f_B=\sum_j(z-a_j)^{-1}-t(z-a_i)^{-1}.
\]

Using rho t=2 gives the individual-tail identity

\[
-\gamma\lambda_1 R_AS_i^Af_B
=\frac{\gamma\lambda_1}{\rho}\partial_z(R_AS_i^A).
\tag{9}
\]

Again no action insertion is used.

In either middle chamber the minority coordinate lies between neighboring majority coordinates. In the initial rho<0 domain, the primitives in (8) and (9) vanish at the two distinct internal endpoints. Integrating by parts transfers the derivative onto an arbitrary smooth test or normal-field coefficient, increasing the full degree from 3(g-rho) to 1+3(g-rho). The proper-face extensions preserve this identity by continuation from the common initial domain. One does not reassign a divergent positive-rho boundary value to zero. Uncontracted C terms and all higher normal-field jets already have at least this improvement. Thus both potentially marginal tail-only middle chambers are regular on a generic physical diagonal.

## 4. Other channels and all proper faces

Directly summing the powers in (3), counting the remaining derivative contractions and the relative spatial measure gives the full-cluster degrees

| Channel | B^3A | A^3B |
|---|---:|---:|
| bare, after (6)–(7) | at least 3(g-rho)+2 | at least 3(g-rho)+2 |
| one K2 | 3(g-rho) | 3(g-rho)+6 |
| two K2 | 3(g-rho)+2 | absent |
| one F_A | 3(g-rho)+6 | 3(g-rho) |
| one L, before epsilon^(2h) | 3(g-rho)-2 | 3(g-rho)+8 |

For a subcluster with a A labels and b B labels, write E=g[binomial(a,2)+binomial(b,2)]-rho ab. If it contains the indicated tail, the exact conservative degrees are

\[
\text{one K2}: E+2a-b+1,
\]
\[
\text{two K2}: E+3a-2b+5,
\qquad
\text{one F_A}: E-a+2b+1.
\tag{10}
\]

A proper subcluster of the two-K2 graph that contains only one tail uses the first line. The fresh script enumerates every subset from the actual vertex records, rather than substituting a source-integration count into a tail graph.

Use the open physical patch 2<t<3, hence 2/3<rho<1. The nontrivial inherited proper degrees include rho-1, 3rho-2, 1-rho and 2-rho. The rho-1 pair needs its ordinary continued subtraction but is not a pole at a generic point of this patch. The remaining displayed degrees are positive there. The two-tail graph adds rho+1, 3rho-1 and 4-rho, all positive. Pure-current proper faces use the same Gaussian boundary extension; their generic pole divisors likewise do not contain the diagonal. Nested or disjoint proper-face extensions therefore introduce no identically diagonal pole into the residual full clusters. The only two possible marginal full channels were removed in (8)–(9).

### The cutoff improvement

The L graph has its explicit epsilon^(2h). In the initial Re h>4 domain its source-free boundary value is zero as epsilon tends to zero. Its analytic continuation is consequently zero. This is consistent with a direct physical-patch check: adding 2t to every L-containing subset degree makes it strictly positive for 2<t<3, including the full B^3A face 2t-2. The script checks this for every subset by an exact positive polynomial after multiplication by t and setting t=2+x, 0<x<1.

Thus there is no finite lambda0 lambda2 contribution to these canonical source-free quartics. This does not mean that action-map counterterms disappear from interacting correlators. Those are not inserted separately into each binary multiplication in this construction.

All middle-chamber families now have regular diagonal restrictions. Applying D(s,s)=0 in (5) proves (1) in the stated canonical product.

## 5. Iterated composition and its exact scope

The source-free calculation admits one common Gaussian multiplication through mass two in the purely chiral charge sector. It need not be defined afresh for each quartic.

1. **Keep the constituent labels.** A fused block records its total physical momentum and the numbers of original A and B labels. The cross-regulator exponent is delta times n_A n'_A+n_B n'_B. This is additive under either bracketing. The cocycle exponent q(m0'+m2') is bilinear in the physical records and obeys its two-cocycle identity.
2. **Compose by Wick cross-contractions.** Multiplying two normal-ordered blocks includes every contraction with one endpoint in each block. Under any binary bracketing of four blocks, every final disjoint matching is obtained once. Unmatched derivative legs retain arbitrary regular normal fields. The fresh check recovers all ten matchings of four derivative legs for all five binary trees; it also checks cocycle and inherited-label factorization for all triples of the five vertex kinds.
3. **Take one compatible boundary value.** Every term through mass two is holomorphic in the current coordinates except for at most the single polynomial right factor from two K2 tails. One F_A or one L has no nontrivial right cross-factor with bare chiral currents. Thus each term is a polynomial times an acyclic chiral boundary-value product.

For completeness, the analytic point in step 3 has a direct proof. Order imaginary heights c1>...>cn. Orient each pair from earlier to later height. A pair factor (xi-xj+i0)^(-pij) is the Fourier transform, up to an entire phase, of k_+^(pij-1)/Gamma(pij). Repeated integration by parts gives an entire one-dimensional distribution family in pij. With qi the outgoing minus incoming edge momentum,

\[
\sum_i c_iq_i=\sum_{i\to j}(c_i-c_j)k_{ij}
\ge c_*\sum_{i\to j}k_{ij},\quad c_*>0.
\tag{11}
\]

Hence the incidence map is proper on the positive edge cone. Pushing forward the tensor product of those entire distributions yields the complete Gaussian boundary value, entire in the pair exponents. Wick derivative contractions shift exponents by integers and preserve the construction. The two-K2 right factor is polynomial; replacing its finite-height version by the real polynomial introduces epsilon times boundary-value families that tend to zero, not a hidden singular opposite-orientation factor.

The separated Wick composition identities thus pass to the same boundary distributions under every bracketing. Smearing independent coordinates with compact tests and grouping the integrations are compatible pushforwards of these common distributions. This defines associative formal composition of labelled smeared Gaussian Wick clusters modulo mass order three in the purely chiral sector. It is not merely a claim that a permutation matrix has a null vector. It supplies the common multiplication in which (1) was tested.

The word “formal” is essential. No Hilbert norm, massive vacuum, unbounded-operator domain, or nonperturbative closure has been constructed. This report also does not show that all mixed barred products are covered by (11): both chiralities can produce opposite-orientation singular factors. Nor does agreement of the chiral quartics prove that a proposed interacting action map is one coherent linear functional on this algebra. Those are independent remaining comparisons.

## 6. Tests: product versus relation

The Gaussian product in section 5 accepts arbitrary independent smooth tests in its slots. The Serre combination (1), however, repeats a single B(f_B) three times or a single A(f_A) three times. The spatial permutation reduction uses that repeated-test symmetry.

For a concrete diagnostic, give the three B occurrences independently supported tests that force their labelled spatial order, and place A in the middle BBAB chamber. Only one labelled B permutation remains. Its coefficient in the formal quartic is

\[
s^{-1}-c_s+c_ss^{-1}-s^{-2}=s^{-1}-s,
\tag{12}
\]

which is generically nonzero. The full Wick kernel in that separated chamber is also generically nonzero. Thus the unpolarized quartic with arbitrary unrelated B slot tests is not a relation of this product. Appropriate polarized identities are a different statement. This distinction does not restrict the tests on which the multiplication itself is defined.

## 7. Canonical multiplication is not the action-dressed kernel term by term

The lower same-chirality pair coefficients show the difference concretely. With the common phases and normal fused fields suppressed, the tail-only BB pair coefficient is 2t(t-2)lambda2 r^(g-2). Adding the one-Y2 action contribution changes it to 2(t-2)lambda2 r^(g-2). Likewise the tail-only AA coefficient is -2gamma lambda1 r^(g-2), while including the one-V1 action term gives -2lambda1 r^(g-2)/(1-t).

These differences do not spoil the canonical quartics: their individual primitives were sufficient. They do prevent identifying canonical multiplication with an action-dressed binary kernel by inspection. If the interaction is represented by one linear state/action map applied after composition, those action terms may arise in its values without redefining multiplication. Establishing that map and its compatibility with the mixed relations is outside this chiral test.

## Evidence and integration

Fresh evidence:

- `referee_scratch/sl21_common_product_20261024/canonical_serre/check_canonical.py` and `check_canonical.out`: **342 PASS, 0 FAIL**. Every separated word/placement/permutation, symbolic shuffle, individual scalar primitives, complete proper-face degrees, and cutoff-improvement positivity.
- `referee_scratch/sl21_common_product_20261024/canonical_serre/check_associative_data.py` and `check_associative_data.out`: **134 PASS, 0 FAIL**. Cocycle/label associativity, all Wick matching trees, distinct-test counterexample, and canonical/action-dressed lower-pair comparison.

Initial diagnostic outputs are retained. The first run encountered a symbolic-boolean conversion error. The initial assertions then exposed an omitted remaining-B contraction in the hand-written *primitive check* and an incorrect expected degree for the minority L tail; these were corrected to the independently generated full Wick formulas. No old input formula was edited. A later positivity run required cancelling a rational expression before polynomial conversion. The final computations pass. These diagnostics distinguish calculation repairs from claims that were simply assumed correct.

**Integration:** none. Independently verified: the full canonical Gaussian exchanges, both source-free quartics through mass two, each marginal tail primitive, proper-face control, and vanishing cutoff improvement. Constructed: compatible formal iterated Gaussian Wick composition in the purely chiral truncated sector. Rejected: the necessity of action/tail cancellation for these quartics and the unpolarized arbitrary-slot-test interpretation. Unresolved here: the complete mixed barred algebra, the interacting state/action-map comparison, all-order conservation, and global quantum sectors. All pre-existing outputs remain unchanged.
