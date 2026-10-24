# Root-system sl(2|1) Toda: the first-order operator Serre test in a fixed analytic prescription

Date: 2026-10-24. Model: Codex (GPT-6).

Commit reviewed: `745586e6a6f7bca5ed91a1e9cf89d6ef5df387f0`.

Task: audit the scalar boundary calculation and extend it to the proper operator-product subtractions. `AGENTS.md` was read. No existing file was edited, no subagent was spawned, and no git mutation, commit, push, or book render was used. All evidence is new and is stored in `referee_scratch/sl21_serre_forest_20261024/`; no previous script was imported or rerun.

## Outcome

The two quartic relations pass through first Yukawa order **as formal finite-interval operator-valued contour identities in the analytic prescription defined below**:

\[
\mathfrak S_{B^3A}=B^3A-c_sB^2AB+c_sBAB^2-AB^3=O(m^2),
\]
\[
\mathfrak S_{A^3B}=A^3B-c_sA^2BA+c_sABA^2-BA^3=O(m^2),
\qquad c_s=1+s+s^{-1},\quad s=e^{-i\pi\rho}.
\tag{F1}
\]

Here the assertion is about the coefficients of $1,\lambda_0,\lambda_2$ in correlators with separated spectators. Word symbols denote the renormalized four-fold contour products defined below. Products include the first-order action insertion and all current-tail placements. It does not assert an unrenormalized equal-time product identity or independently construct a globally associative operator algebra. It does not prove the mass-squared Serre relations, global conserved charge operators, physical boundary closure, or a quantum-group presentation.

The new point is not just another scalar cancellation. A two-parameter chamber calculation exhibits precisely how a vanishing permutation coefficient could multiply a divergent integral and leave an anomaly. The potentially dangerous chamber pole is absent. The same continuation then fixes the embedded pair/triple counterterms, their descendants, and their endpoint terms together. An explicit nested $BB$-inside-$BBB$ calculation verifies the local matching rather than assuming it from the permutation matrix.

There are **94 new PASS, 0 FAIL** checks across three scripts. They support, but do not replace, the distributional argument in sections 4–6. No [Theorem] label is recommended and no physical symmetry algebra is identified.

| Item | Verdict | Evidence and limitation |
|---|---|---|
| Leading scalar derivative in the preceding boundary report | **Correct** | A new Wick reconstruction and chamber Stokes checks |
| Tail and lower-composite separated phases | **Correct in the fixed cocycle convention** | Actual ordered vertex powers, signs, and full Wick kernels with varying backgrounds |
| Embedded $BB$ subtraction and its overlap with $BBB$ | **Correct** | General-shape descendant reconstruction and an independent fused-pair OPE |
| Inferring (F1) just from a null permutation coefficient | **Wrong as a proof** | A simple chamber pole would give a finite nonzero limit |
| Absence of that first-order pole | **Correct in the stated analytic family** | The middle-chamber scalar primitive has zero internal endpoints; nonconstant terms have one extra collision power |
| First-order identities (F1) in this prescription | **Correct as formal contour identities** | A common initial regulator domain, meromorphic continuation, and regular diagonal restriction |
| Equality with an arbitrary physical-coupling sharp-cutoff finite part | **Unclear** | Such a product can have a different finite normalization |
| Mass-squared Serre relations and an established symmetry algebra | **Unclear** | Not calculated or inferred here |

## 1. Reconstruct the contractions and keep the normalization fixed

The chiral Gaussian convention is $X(z)X(w)\sim-\log(z-w)$. With

\[
a_0=(ib/\sqrt2,b/\sqrt2,i),\quad a_2=(ib/\sqrt2,-b/\sqrt2,-i),
\quad t=b^2,\quad d=-(a_0+a_2)/t,\quad \rho=2/t,
\]

one obtains $a_i^2=-1$, $a_0a_2=1-t$, $a_id=1$, and $d^2=-\rho$. Use

\[
J_A=:e^{-dX}:,\quad J_B=:C e^{dX}:,\quad C=a_0\partial X,
\quad H_B^{(1)}=-t\lambda_2K_2,
\quad K_2=:e^{(d+a_2)X+a_2\bar X}:.
\tag{F2}
\]

The first-order action is $S^{(1)}=\pi^{-1}\int d^2w(\lambda_0Y_0+\lambda_2Y_2)$, so an action insertion has coefficient $-\lambda_i/\pi$. Put $\beta=t-1$, $\eta=t-2$. Direct contraction gives $CC\sim(z-w)^{-2}$, $C Y_2\sim\beta Y_2/(z-w)$ and $C K_2\sim\eta K_2/(z-w)$, apart from the accompanying vertex factors.

These signs also have a direct fermionic check: $C=p_0\partial\varphi+j$, $j=:\chi_0\chi_2:$, with $j\chi_2\sim-\chi_2/(z-w)$ and $jj\sim(z-w)^{-2}$. The bosonic contribution to $C Y_2$ is $t/(z-w)$; adding the fermion-current contribution gives $t-1$. The extra $d$ in $K_2$ subtracts one more. Thus the derivative prefactor is not being replaced by a bare exponential.

Retain the previously specified local cocycle: a current of momentum $nd$ placed before a full Yukawa has sign $(-1)^n$. The present calculation checks consequences of that convention; it does not establish its realization by massive Hilbert-space operators. All densities in (F2), including the full tail, are even. Charge factors $(2\pi i)^{-n}$ are common to the ordered formulas and suppressed below.

## 2. Define one analytic product before imposing the quartic relation

Give the exponent between every pair of like current labels an independent value $g$, leaving the unlike-current exponent $-\rho$ unchanged. Equivalently multiply each ordered Gaussian kernel by

\[
\prod_{i<j\, :\, L_i=L_j}(z_i-z_j)^{g-\rho},
\qquad L_i\in\{A,B\}.
\tag{F3}
\]

A slot containing $K_2$ still has label $B$. Its additional integer powers and cocycle remain intact. Screening-current powers, normal fields, derivative contractions, and $t=2/\rho$ are not changed. This is an auxiliary regulator of products, not a second physical coupling or a replacement Gram matrix.

For each word, use identical finite spatial intervals, a common smooth infrared cutoff on the action insertion, and spectators separated from the current patch. Use time levels $c\epsilon(3/2,1/2,-1/2,-3/2)$ in operator order, with fixed $c>2$. Keep radius-$\epsilon$ source holes at inner currents and the circular Cauchy extension at the outer current. The first-order left-hand expansion contains one action insertion or one tail, with the coefficients in (F2); it is not defined by setting a desired Serre defect to zero.

First take $\epsilon\to0$ in the conservative common domain

\[
\mathcal D:\qquad \operatorname{Re}g>2,\qquad \operatorname{Re}\rho<-1.
\tag{F4}
\]

Section 5 checks this domain. The resulting real-axis products have common circular Cauchy source extensions at resolved current positions. Continue the **complete left-hand kernels** meromorphically in $(g,\rho)$ by local Taylor subtraction, and restrict to $g=\rho$ only where this restriction is regular. No independent finite quartic counterterm is added. The first target strip is

\[
-\tfrac13<\operatorname{Re}\rho<0.
\tag{F5}
\]

This defines a particular analytic finite part of the nested regulator. It does not interchange the physical-diagonal limit with the cutoff limit. Its continuation can be used beyond (F5) at generic coupling, but further explicit subtraction levels and exceptional logarithmic prescriptions are not catalogued here.

An important detail concerns the lower counterterms. Before restricting to the diagonal, the pair bulk coefficient is

\[
\frac{2(t-2)}{g-1},\qquad
\frac{2(t-2)}{g-1}+2t=\frac{2t(g-\rho)}{g-1}.
\tag{F6}
\]

It must not be replaced prematurely by its diagonal value $-2t$. All radial powers, angular beta functions, and their denominators use $g$; the Wick coefficients still use $t=2/\rho$.

## 3. Derive the full permutation defect, including the tails

Write

\[
u=e^{-i\pi g},\qquad s=e^{-i\pi\rho}.
\]

Separated exchanges have $q_{AA}=q_{BB}=u$, $q_{AB}=q_{BA}=s^{-1}$. These statements also hold when a $B$ is replaced by its first tail. For example, the bare exponent between $B$ and $K_2$ becomes $g-1$. Its reversal gives the phase $e^{-i\pi(g-1)}$, and the current-first Yukawa cocycle contributes $-1$. Their product is $u$. Between $A$ and $K_2$ the same calculation is $-e^{-i\pi(1-\rho)}=s^{-1}$. Derivatives add no separated phase; their contact terms are part of the extension problem, not discarded.

For an embedded composite $\mathcal V_n=:e^{(nd+a_2)X+a_2\bar X}:$, the same reasoning gives phases $s^n$ with $B$ and $s^{-n}$ with $A$ on the physical diagonal. At independent $g$, the $n$ fused $B$ labels carry the corresponding factor $u^n$ against a remaining $B$. Neutral derivatives in descendants do not alter these separated phases.

Partition the four-current spatial cube into labelled ordered chambers, and collect the four species orders

\[
W=(BBBA,BBAB,BABB,ABBB).
\]

Let $M(u,s)$ be the matrix whose entry sums the inversion phases of all labelled permutations taking its operator-word column into its spatial-word row. Direct enumeration of the $24$ permutations per column, with no assumed Cartan matrix, gives

\[
M(u,s)
\begin{pmatrix}1\\-c_s\\c_s\\-1\end{pmatrix}
=\begin{pmatrix}0\\-D(u,s)\\D(u,s)\\0\end{pmatrix},
\]
\[
D(u,s)=\frac{(s-1)(s+1)(s-u)(s^2-u)(u+1)}{s^4}.
\tag{F7}
\]

The same matrix applies with $A$ and $B$ exchanged. It applies separately to either Yukawa action insertion and to the sum over first-tail placements, because all the corresponding full Wick kernels obey the separated exchanges just checked.

Let $\mathcal I_j^{(1)}(g,\rho)$ denote a spatial chamber integral of the full first-order kernel, including its normal fields and spectators. In the initial domain (F4), and hence as a meromorphic identity,

\[
\mathfrak S^{(1)}(g,\rho)
=D(u,s)\,[\mathcal I_{BABB}^{(1)}-\mathcal I_{BBAB}^{(1)}].
\tag{F8}
\]

The two outer chambers have zero coefficients even off the diagonal. The two middle chambers generally do not. Thus (F8) does not impose the desired answer as a regulator convention.

The script checks actual ordered vertex powers and cocycles, rather than assigning the desired phase to a precomputed kernel. It evaluates all $96$ labelled word/chamber assignments for each of the two quartics, both action nodes, and the tail sum, at two generic exponent choices. The regular fields vary with position: $C(x)$ is quadratic, $dX(x)$ is quadratic, and $a_2X=-t\,dX-a_0X$ is enforced. All full Wick comparisons pass; the largest scaled word-sum discrepancy is below $3.5\times10^{-14}$.

### Why a zero permutation coefficient alone would not suffice

Put $g=\rho+\delta$. A chamber pole $1/(3\delta)$ would leave

\[
\lim_{\delta\to0}\frac{D(e^{-i\pi(\rho+\delta)},s)}{3\delta}
=\frac{i\pi}{3}(s-s^{-1})^2,
\tag{F9}
\]

which is generically nonzero. This is a diagnostic example, not an anomaly found in the actual calculation. It identifies the missing step in a proof based only on separated braiding: the relevant individual chamber integrals must have no pole along $g=\rho$.

## 4. Audit the lower subtractions and their overlaps

The area projection used for the low jets follows directly from Green's theorem. For distinct real positions, $F(w)=\prod_i(w-z_i)^{-2}$ and a polynomial $P(w,\bar w)$ of total degree less than $2n-2$ obey

\[
\pi^{-1}\operatorname{FP}\int d^2w\,F(w)P(w,\bar w)
=\sum_i\frac{P(z_i,z_i)}{\prod_{j\ne i}(z_i-z_j)^2}.
\tag{F10}
\]

For a monomial, integrating once in $\bar w$ reduces the answer to simple residues. The residue at infinity vanishes under the degree restriction. The new script checks every required moment at generic positions, rather than at a single numerical shape. At a source pole, the double coefficient is $\beta$ times the corresponding fused-tail polynomial. Subtracting it from the tail coefficient $t$ gives $t-\beta=1$. This reconstructs the combined source/tail jets with their quantum normalization.

### Pair subtraction and its embedding

At midpoint $y$, the leading pair field is

\[
2\lambda_2\eta r^{g-2}\mathcal V_2(y),\qquad r>0,
\]

and the first midpoint jet is zero. On a hard interval $[L,R]$, the analytic representative adds

\[
CT_{BB}(g,\rho)=\lambda_2\left[
\frac{2\eta}{g-1}\epsilon^{g-1}\int_L^R\mathcal V_2
-\frac{\eta}{g}\epsilon^g\bigl(\mathcal V_2(L)+\mathcal V_2(R)\bigr)
\right]
\tag{F11}
\]

in the target strip. The second term comes from the shrinking midpoint domain. It is retained even though the bulk first jet vanishes. At $g=\rho$, this agrees with the preceding boundary report's pair subtraction.

For a nontrivial nested test, set the pair midpoint to $R_*/3$ and the third $B$ to $-2R_*/3$. Independently contract the third current with $\mathcal V_2$. Its derivative coefficient is

\[
C(-2R_*/3)-\frac{\eta-1}{R_*}.
\]

Its remaining vertex power is $R_*^{2g-1}$. Expanding the normal exponential at the triple centroid gives

\[
1+\frac{R_*F_1}{3}
+R_*^2\left[\frac{D_2}{3}+\frac{F_2+F_1^2}{18}\right]+\cdots,
\]

where $D=d\partial X$, $D_2=d\partial^2X$, $C_2=a_0\partial^2X$, $U=a_2\bar\partial\bar X$, $U_2=a_2\bar\partial^2\bar X$, and $F_1=-tD-C+U$, $F_2=-tD_2-C_2+U_2$. After stripping $\lambda_2r^{g-2}R_*^{2g-2}\mathcal V_3$, the first three coefficients are

\[
-2\eta(\eta-1),\qquad
2\eta[C-(\eta-1)F_1/3],
\]
\[
2\eta\left[-\frac{2C_2}{3}+\frac{CF_1}{3}
-(\eta-1)\left(\frac{D_2}{3}+\frac{F_2+F_1^2}{18}\right)\right].
\tag{F12}
\]

These agree exactly with the pair subface of the independently reconstructed triple jets below. No spurious $r^{g-1}$ midpoint pole appears. A separate symbolic test keeps a third $B$ and an $A$ at arbitrary separated positions and matches the embedded pair coefficient to the four-current scalar primitive. This checks the counterterm as an insertion, not just as an isolated two-point coefficient.

### Triple subtraction and endpoint terms

At centroid offsets $r_i$ with $\sum_i r_i=0$, set

\[
D_i=\prod_{j\ne i}(r_i-r_j)^{-2},\qquad
M_j=\sum_i r_i^jD_i,\qquad P_2=\sum_i r_i^2.
\]

Stripping $\lambda_2\Delta^g\mathcal V_3$, the divergent jets are

\[
K_0=-\eta(\eta-1)M_0,\qquad
K_1=-\eta[(\eta-1)F_1-3C]M_1,
\]
\[
K_2=-\eta\left[
\left(\frac{\eta-1}{2}D_2+C_2\right)P_2M_0
+\left(\frac{\eta-1}{2}(F_2+F_1^2)-3CF_1\right)M_2
\right].
\tag{F13}
\]

The new derivation uses the two remaining Gaussian derivatives at each tail location. Their Wick polynomial is

\[
\left(C_j+\frac\eta{r_j-r_i}-\frac1{r_j-r_k}\right)
\left(C_k+\frac\eta{r_k-r_i}+\frac1{r_j-r_k}\right)
+\frac1{(r_j-r_k)^2}.
\]

Multiply by $-[(r_j-r_i)(r_k-r_i)]^{-1}$, expand its normal exponential, and sum $i$. The check is symbolic in the general shape variable $x$, in $t$, and in the independent normal-field jets. It reproduces (F13) and then independently verifies (F12).

For the span cutoff, put $r=z_1-z_3$, $x=(z_1-z_2)/r$. The centroid interval is $[L+r(2-x)/3,R-r(1+x)/3]$. The full prescription is to continue the angular integral with weight $[x(1-x)]^g$ from $\operatorname{Re}g>2$, then subtract the first three radial Taylor coefficients $A_j(g,\rho)$ of this changing-domain integral:

\[
CT_{BBB}(g,\rho)=\lambda_2\sum_{j=0}^2
\frac{\epsilon^{3g-2+j}}{3g-2+j}A_j(g,\rho).
\tag{F14}
\]

This definition includes the endpoint expansion of $K_0$, the endpoint contribution of the angular-odd $K_1$, and the bulk $K_2$ term. It specializes to the previous report's explicit descendant/endpoint formulas. Discarding the odd jet before expanding the moving centroid domain would change (F14).

The remaining radial powers in (F5) are $g+1$ for the pair and $3g+1$ for the triple. The proper pair poles of the angular triple integral are fixed by the same analytic continuation as (F11), and the explicit overlap test (F12) confirms their operator coefficient. Angular-then-radial and sector-wise Taylor subtraction are two evaluations of this same meromorphic family because they agree in (F4). No freely chosen finite term is introduced at their overlap.

### Why the embedded subtraction sum is compatible with the word sum

In a fixed spatial chamber, local Taylor extraction acts on the Wick coefficient and the smooth spectators. The permutation weights in $M(u,s)$ are constants in all insertion coordinates. Therefore extraction of any bulk jet, endpoint jet, or nested/disjoint forest coefficient commutes with multiplication by $M$. At the level of arbitrary independent chamber jets $F_j$, its word sum is exactly

\[
D(u,s)(F_{BABB}-F_{BBAB}).
\tag{F15}
\]

This statement by itself would not prove the proper coefficients were normalized correctly. Equations (F10)–(F14), their fusion phases, and the explicit overlap checks supply that separate verification. Conversely, repeating the same subtraction with a different finite normalization in different words would invalidate (F15); it is not permitted by the common initial family.

The proper $ABB$ scalar primitive retains the infinity term found in the preceding audit: $-\lambda_2\partial_a[\mathcal R_2(S_2+\beta)]/\rho$. It is included by the same Wick/area continuation. In the target strip the corresponding full source face has degree $1+g-2\rho$, positive on the diagonal. Its finite scalar expression does not create an additional diagonal pole. No special $ABB$ subtraction is chosen to enforce (F1).

## 5. Check the initial regulator limit and enumerate all possible poles

For a cluster with $n_A$ bare $A$ labels, $n_B$ bare $B$ labels and one Yukawa, the current vertex degree is

\[
g\left[\binom{n_A}{2}+\binom{n_B}{2}\right]-\rho n_An_B.
\]

The source powers contribute $n_A-n_B$, and the maximally contracted derivatives contribute $-n_B$. The relative measure has dimension $n_A+n_B+1$; the common coordinate along the interval is not a small variable. Hence a raw relative bound is

\[
\omega=g\left[\binom{n_A}{2}+\binom{n_B}{2}\right]-\rho n_An_B
+2n_A-n_B+1.
\tag{F16}
\]

The corresponding tail cluster has exactly the same degree. The relevant instances are:

| Cluster including the source, or its tail version | Raw relative degree |
|---|---:|
| $BB$ | $g-1$ |
| $BBB$ | $3g-2$ |
| $AB$ | $2-\rho$ |
| $ABB$ | $1+g-2\rho$ |
| $AAB$ | $4+g-2\rho$ |
| $BBBA$ | $3(g-\rho)$ |
| $AAAB$ | $6+3(g-\rho)$ |

Pure-$A$ source clusters are softer. A single $B$/source pole has the specified circular Cauchy extension. Its constant and linear angular moments do not produce a finite disk term; the first nonzero disk contribution is proportional to $\epsilon^2$ times a smooth second jet. Direct disk quadrature checks this with a nonconstant complex polynomial background.

For source-free faces the actual derivative cancellations must also be retained. For example,

\[
J_B(z_1)J_B(z_2)=(z_1-z_2)^g:
\left[C_1C_2+\frac{C_1-C_2}{z_1-z_2}\right]e^{dX_1+dX_2}:,
\]

so the apparent double pole cancels. In the full source-free $BBBA$ cluster, the scalar contraction and the constant first regular-derivative jet both vanish. The new script checks them independently; the resulting full-cluster degree is at least $2+3(g-\rho)$. The $AAAB$ source-free bound is the same. A distant Yukawa acts as a smooth spectator for these estimates.

All genuine multiple-current source faces in the table, and all source-free faces, have positive degrees in (F4). Thus the finite-height and finite-hole regulators approach the resolved Cauchy kernels there. To see why the shrinking holes are included in this assertion, the pair scalar hole term has the bound

\[
\epsilon^2\int_\epsilon^{R_0}r^{\operatorname{Re}g-4}\,dr
=O(\epsilon^{\min(\operatorname{Re}g-1,2)}(1+|\log\epsilon|)),
\]

and the height term is bounded by $\epsilon\int_\epsilon^{R_0}r^{\operatorname{Re}g-3}dr=O(\epsilon)$ in the conservative domain. Triple and full-cluster terms have the positive degrees in the table. Their proper pair subregions use the preceding estimate. Pure current $AB$ boundary-value differences vanish with positive power $-\operatorname{Re}\rho$. A cluster pinned at an interval endpoint gains one relative dimension.

These estimates extend to the full Wick kernels by expanding the other factors in the smallest cluster and using the circular Cauchy operation for a single source pole. For a nested cluster the Taylor remainder supplies additional nonnegative powers; disjoint clusters use separate relative variables. Sectoring by the ordered relative distances leaves compact shape regions, with the radial powers listed above. The conservative margins in (F4) allow the necessary fixed-order Cauchy differentiations and cutoff-boundary operations. Keeping $c>2$ prevents intersecting source disks. This proves matching to the finite-height/circular-hole limit in an open domain; it does not claim matching after taking a bare physical-diagonal limit first.

Taylor integration in these sectors gives meromorphic kernels with Mellin denominators of the form $\omega+j$, $j\ge0$. In a neighborhood of a generic diagonal point satisfying (F5), the proper pair/source denominators are $g-1+j$, and the triple/source denominators are $3g-2+j$. Their poles are transverse to the diagonal and do not pass through that strip. The mixed source faces have positive diagonal degrees $2-\rho$, $1-\rho$, and $4-\rho$. The source-free proper faces likewise have no diagonal pole there. Endpoint faces add a positive power. A lone Cauchy source pole has already been extended and supplies no parameter pole.

This exhausts proper subsets: a subset either contains the one source or is a source-free current subset; the displayed counts cover the possible numbers of each species in either quartic. Intersections are nested or disjoint combinations of these subsets and are handled by the same sector/Taylor prescription. Their coefficient functions remain regular near the chosen diagonal point once these proper factors have been extended.

Only the full $BBBA/Y_2$ scalar channel can have a pole along $g=\rho$. Its higher jets have degrees $3(g-\rho)+j$, $j\ge1$. The node-0 scalar vanishes. The full $AAAB$ channel has degree six on the diagonal. Thus one specific primitive pole, not an unspecified collection of contacts, remains to be checked.

## 6. Remove the possible primitive pole before specializing the matrix

For the three $B$ positions and an $A$ at $a$, put

\[
\mathcal R_{g,\rho}=\Delta^g\prod_i(z_i-a)^{-\rho},\qquad
S=\sum_i\frac{(a-z_i)^2}{\prod_{j\ne i}(z_i-z_j)^2}.
\]

A fresh Wick reconstruction gives the leading source kernel

\[
\frac{\beta\eta}{\rho}\partial_a\left[
\mathcal R_{g,\rho}\frac{(a-w)^2}{\prod_i(w-z_i)^2}\right].
\tag{F17}
\]

It remains a derivative at independent $g$, because $\Delta^g$ is independent of $a$. The three tails have the corresponding primitives $t\lambda_2\eta\mathcal R_{g,\rho}S_i/\rho$. Integrating the resolved source and combining them gives the lowest coefficient

\[
\partial_aP,\qquad P=\frac{\lambda_2\eta}{\rho}\mathcal R_{g,\rho}S.
\tag{F18}
\]

The finite-height and hole terms are not ignored in defining the family: they vanish in (F4) by section 5. The earlier exact finite-hole formula is compatible with this limit. Its analytic continuation is consequently fixed by this same initial family.

In either **middle spatial chamber**, the $a$ interval lies between two consecutive $B$ positions. For fixed distinct $B$ positions, $S$ is finite at these endpoints and $\operatorname{Re}\rho<0$ implies $P\to0$ there. For a smooth nonconstant coefficient $f$, integration by parts therefore gives, initially in (F4),

\[
\int_{a_-}^{a_+}(\partial_aP)f\,da
=-\int_{a_-}^{a_+}P\,\partial_af\,da.
\tag{F19}
\]

The boundary at a simultaneous proper collision is fixed by continuation of this identity, not by replacing a divergent boundary value with zero by hand. Both sides have the same proper-face extensions from section 4. Differentiation and these continued identities commute because they agree on their common initial domain.

The primitive $P$ has one more integrated relative power than $\partial_aP$. The right-hand side of (F19) has full-cluster degree

\[
1+3(g-\rho).
\tag{F20}
\]

Thus it is regular at $g=\rho$ after the proper faces have been extended. In particular, for the constant leading field $f$, a middle-chamber scalar integral is zero even before restriction. For nonconstant $f$, the surviving derivative term has positive degree one on the diagonal. The remaining normal-field Taylor levels and uncontracted derivatives have the same positive improvement. A cutoff chosen equal to one near the full collision has derivative supported away from that collision, so it cannot recreate the primitive pole.

This proves regularity of each middle chamber, which is stronger than cancellation of their scalar sum. For an outer chamber, integration by parts also has a fixed spatial endpoint. A full cluster at that endpoint gains a relative dimension and has the same positive degree (F20). Thus the individual word products, not just their difference in (F8), have regular generic diagonal restrictions.

The numerical check uses Jacobi quadrature for the exact endpoint powers in (F19). Both middle chambers have constant-field residuals below $5\times10^{-12}$ and satisfy the nonconstant linear-field identity. On compact shape domains, the radial integral of the linear-field remainder is explicitly proportional to $[1+3(g-\rho)]^{-1}$. At $\rho=-\sqrt2/5$ the largest chamber magnitude stays between $0.4412$ and $0.4884$ while $g-\rho$ ranges from $0.02$ to $-0.0025$. The corresponding permutation defect decreases from $0.06205$ at $0.02$ to $0.007844$ at $0.0025$. This is a check of the primitive regularity mechanism; the proper shape boundaries are handled by the independent overlap calculation and the meromorphic argument, not by those clipped numerical integrals.

Now specialize (F8). Both middle chamber integrals are regular and $D(s,s)=0$. There is no pole to produce the finite remnant illustrated in (F9). The first-order $B^3A$ defect is zero. The $A^3B$ proof uses the same full-kernel permutation identity; its full source cluster has positive degree six, so it has no corresponding primitive obstruction. The zero-order argument uses the independently checked source-free Wick cancellations. This establishes (F1) within the definition in section 2.

The full scalar distribution may contain contact terms fixed by its derivative extension. The proof says their quartic sum vanishes in this prescription; it does not set every individual contact distribution to zero.

## 7. Checks, limits, and the next calculation

The new independent scripts and outputs are:

| Script | Output | PASS | Scope |
|---|---|---:|---|
| [forest_checks.py](../referee_scratch/sl21_serre_forest_20261024/forest_checks.py) | [forest_checks.out](../referee_scratch/sl21_serre_forest_20261024/forest_checks.out) | 47 | Two-exponent permutation matrix, arbitrary-jet word sums, full nonconstant Wick kernels and cocycles, scalar derivative |
| [subface_checks.py](../referee_scratch/sl21_serre_forest_20261024/subface_checks.py) | [subface_checks.out](../referee_scratch/sl21_serre_forest_20261024/subface_checks.out) | 25 | General-shape triple descendants, embedded pair overlaps, fusion phases, middle-chamber Stokes and regularity tests |
| [regulator_audit.py](../referee_scratch/sl21_serre_forest_20261024/regulator_audit.py) | [regulator_audit.out](../referee_scratch/sl21_serre_forest_20261024/regulator_audit.out) | 22 | Independent-exponent normalization, relative counts, pole-times-zero diagnostic, Cauchy disk, generic-position source moments |

The next perturbative calculation is now the order-$m^2$ Serre defect in the **same explicitly defined product normalization**. It must include both Yukawa insertions, all first-tail/action cross terms, the two-tail terms, the composite-potential insertion, and the mass-squared single-current counterterms. The source-source exponents involve $t$ at that order; the first-order pole analysis cannot simply be reused. No mass-squared conclusion is claimed in this report.

Physical boundary closure remains a separate unresolved problem. The fermionic mixed boundary operator has not been identified with a Cartan generator by this calculation. Neither a physical sign choice for the deformation parameter, a coproduct on massive sectors, nor a faithful representation follows from (F1). The result remains formal and prescription-dependent, even though the first-order operator-product step is now more complete.

No correction to an established formula in the preceding boundary report is recommended. Its conclusion that the full operator first-order step was then unfinished is accurate for that report. The additional proof here uses the entire two-parameter family, verifies its lower-face compatibility, and checks the otherwise possible diagonal pole. The earlier scalar continuation alone did not contain these steps.

## Integration

No book integration was performed. **Independently verified:** the full first-order separated kernels and tail signs, generic lower-product jets, nested pair/triple matching, and the marginal scalar derivative. **Extended:** a common operator-valued analytic product prescription and the absence of a diagonal chamber pole, yielding the first-order formal quartic identities (F1). **Not established:** arbitrary-scheme or physical equal-time realization, mass-squared Serre identities, global boundary closure, or a symmetry algebra. All pre-existing reports, scripts, and book files remain unchanged.
