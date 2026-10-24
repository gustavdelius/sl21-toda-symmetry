# Root-system sl(2|1) Toda: same-chirality quartic relations and the first interaction contacts

Date: 2026-10-23. Model: Codex (GPT-6).

Commit reviewed: `745586e6a6f7bca5ed91a1e9cf89d6ef5df387f0`.

Task: begin the same-chirality calculation after the independent mixed-product audit, without choosing a symmetry algebra in advance. Read `AGENTS.md`. No existing file was edited, no subagent was spawned, and no git mutation, commit, push, or book rendering was used.

Independent evidence: [chiral_checks.py](../referee_scratch/sl21_same_chirality_20261023/chiral_checks.py), [its output](../referee_scratch/sl21_same_chirality_20261023/chiral_checks.out), and the derivations below. No preceding check script was imported or rerun. There are **61 PASS, 0 FAIL** checks. A preliminary run interrupted during unnecessarily expensive symbolic simplification is retained separately; it was replaced by exact denominator-clearing, not by weakened identities.

## Outcome

Both ultraviolet quartic relations follow from the actual currents, including the derivative prefactor. Define

\[
\rho=\frac2t,\qquad s=e^{-i\pi\rho}=\nu^2,\qquad
c_s=s+1+s^{-1}.
\tag{S1}
\]

Here $t=b^2$ is the coupling parameter of the preceding reports. The new symbol $\rho$ is just an OPE exponent, not a new physical coupling identification. For the ultraviolet finite-interval contour integrals $A=Q_A^{(0)}$, $B=Q_B^{(0)}$, the calculation gives

\[
\boxed{A^3B-c_sA^2BA+c_sABA^2-BA^3=0,}
\tag{S2}
\]

\[
\boxed{B^3A-c_sB^2AB+c_sBAB^2-AB^3=0.}
\tag{S3}
\]

The products are first defined with separated spectators and identical finite intervals in the convergent region $-1/2<\operatorname{Re}\rho<0$, then continued meromorphically. The coefficient and the degree of these relations were obtained from labelled permutation sums; neither was selected from a proposed Cartan matrix. The derivative current has an exact self-OPE cancellation that is needed to justify these integrals.

The first interaction correction is not settled by this result. It produces new same-chirality composite-product singularities. In a resolved ordered pair, after integrating the action insertion with its Cauchy prescription and including both first-order current tails, the leading term is

\[
\boxed{
\bigl[(J_B+H_B)(x+r)(J_B+H_B)(x)\bigr]_{\lambda_2}
\sim 2\lambda_2(t-2)r^{\rho-2}\mathcal V_2(x),\qquad r>0,
}
\tag{S4}
\]

where $\mathcal V_n=:e^{(nd+a_2)X+a_2\bar X}:$ denotes the corresponding fermionic Wick field. For example, $\mathcal V_2=i:\chi_2\bar\chi_2e^{(2d+p_2)\varphi+p_2\bar\varphi}:$ in the fixed local cocycle convention. Formula (S4) concerns the leading fused field, not the complete OPE or a finished product subtraction.

There is also a useful reduction of the potentially marginal four-current/one-Yukawa collision: its lowest fused coefficient is an exact derivative in the $A$-current coordinate. The current-tail terms have compatible derivative formulas. This narrows the remaining calculation to subtraction and contour boundaries. It does **not** allow those boundaries to be discarded.

Thus (S2)–(S3) are established here as **formal ultraviolet contour relations**, while their validity for the massive charges remains **unclear**. The mixed boundary operators, global charge realization, and full symmetry algebra are still unresolved.

| Item | Verdict | Evidence and scope |
|---|---|---|
| Same-chirality separated braiding | **Correct** | Canonical bracket and direct current OPEs |
| A generic odd or nilpotent interpretation of the $B$ generator based on its fermionic dressing | **Wrong** | The density contains an even bilinear; its generic quadratic permutation coefficient is nonzero |
| Quartic relations (S2)–(S3) for the ultraviolet products | **Correct in the stated analytic prescription** | Exact permutation kernels, derivative-current Wick cancellations, and collision bounds |
| Inferring these operator identities from the exponential Gram matrix alone | **Unclear without the additional calculation** | The derivative poles require the checks in sections 2–3 |
| Absence of a first-order same-chirality product subtraction | **Wrong** | Explicit combined action/tail coefficient (S4) |
| Leading $BB$ and $BBB$ interaction kernels below | **Correct in the resolved-point, source-first calculation** | Wick contractions, cocycle signs, area integration, and independent quadrature |
| Exact-derivative reduction of the marginal $BBBA/Y_2$ kernel | **Correct** | Symbolic identity using $\rho t=2$, including all three tail placements |
| Massive Serre relations through first or second mass order | **Unclear; not yet established** | Proper product subtractions, contour-height terms, and cutoff/end boundaries remain to be combined |
| Affine symmetry algebra or a faithful presentation of all charge relations | **Unclear** | Neither follows from two ultraviolet quartic identities |

## 1. Conventions and what the sources do not supply

The local propagators, momenta and current representatives are

\[
X^a(z)X^b(w)\sim-\delta^{ab}\log(z-w),\qquad
a_0^2=a_2^2=-1,\quad a_0a_2=1-t,
\]

\[
d=-\frac{a_0+a_2}{t},\qquad d^2=-\rho,\qquad a_0d=a_2d=1,
\]

\[
J_A=:e^{-dX}:,\qquad J_B=:C e^{dX}:,\qquad C=a_0\partial X.
\tag{S5}
\]

They were recomputed from $a_0=(ib/\sqrt2,b/\sqrt2,i)$ and $a_2=(ib/\sqrt2,-b/\sqrt2,-i)$. Locally $C=p_0\partial\varphi+j$, $j=:\chi_0\chi_2:=i\partial h$. The derivative of the physical boson and the fermion bilinear are both even. Direct fermionic Wick contraction gives $j(z)j(w)\sim(z-w)^{-2}$; $p_0^2=0$ then gives $C(z)C(w)\sim(z-w)^{-2}$. This provides the derivative contraction used below without relying only on bosonization.

The action and first tail are

\[
S_{\rm int}=\frac1\pi\int d^2w\,
(\lambda_0Y_0+\lambda_2Y_2+\lambda_1V_1+\kappa P),\qquad
H_B^{(1)}=-t\lambda_2K_2,
\tag{S6}
\]

where $K_2=:e^{(d+a_2)X+a_2\bar X}:$, the Yukawas have mass order one, and $\lambda_1,\kappa$ have mass order two. This report computes first-Yukawa product kernels, so the order-two potential and current counterterms do not enter those particular formulas. They will be required in the subsequent mass-squared Serre test.

The local charge-cover cocycle is the one previously fixed: a current of momentum $nd$ placed before a full Yukawa contributes $(-1)^n$. It is included in the signs below. The cross-chiral phase convention is unchanged.

The supplied Bernard–LeClair PDF was inspected at printed page 111. In its discussion of (3.20), the phrases **“Only the Serre relations”** and **“are missing in (3.20)”** explicitly separate that calculation from the Serre problem. The source does not supply the needed relation for the present derivative-dressed current.

For comparison, Flandoli–Lentner, *Algebras of Non-Local Screenings and Diagonal Nichols Algebras*, [SIGMA 18 (2022), 018](https://sigma-journal.com/2022/018/sigma22-018.pdf), section 3.2, states the screening result subject to **“provided that $(m_{ij})$ is subpolar”**. Its Definition 3.3 requires positive relative collision degrees for every relevant subset. That convergence qualification is pertinent here. Its theorem for exponential lattice screenings is not imported as a proof for $J_B$; the derivative cancellations and integrals below are checked directly. No novelty claim is made.

## 2. Derive the ultraviolet relations from ordered products

The working canonical split gives $[\varphi(x),\varphi(y)]=\pi i\operatorname{sgn}(x-y)$. Away from coincidence, derivatives of the sign function vanish, so the derivative prefactor adds no separated exchange phase. With momenta $-d,d$ one obtains, for $x>y$,

\[
J_i(x)J_j(y)=q_{ij}J_j(y)J_i(x),\qquad
(q_{ij})=\begin{pmatrix}s&s^{-1}\\s^{-1}&s\end{pmatrix}.
\tag{S7}
\]

Coincident derivative contacts are not being set to zero by (S7). Section 3 is needed before this separated formula can be integrated.

For a word of $n$ labelled current insertions, partition the product interval into its $n!$ ordered chambers. Reordering each chamber into a fixed spatial order multiplies it by the product of $q_{ij}$ for its inversions. Keep identical currents labelled until the permutation weights have been summed. This gives a matrix $M$ whose columns are operator words and whose rows are spatially ordered species words.

The script constructs these matrices by enumerating permutations. In the degree $(3,1)$ sector, the four columns are $AAAB,AABA,ABAA,BAAA$. Their unique generic null vector, normalized by its first component, is

\[
(1,-s-1-s^{-1},s+1+s^{-1},-1).
\tag{S8}
\]

The degree $(1,3)$ calculation gives the same coefficients in the order $BBBA,BBAB,BABB,ABBB$. These are (S2)–(S3). The script extracts the nullspaces and also verifies their multiplication by $M$ identically. An exact specialization at the formal value $s=2$ establishes generic rank three; this specialization is an algebraic rank witness, not a physical coupling choice.

The other permutation matrices through degree four have full generic rank. In particular, the quadratic $BB$ coefficient is $1+s$, not zero at generic coupling. These rank statements concern the permutation map. They do not prove that its realization by the field theory is faithful or exclude additional representation-dependent operator identities.

The familiar quartic Serre form is therefore a result of the calculation. Since $c_s=\nu^2+1+\nu^{-2}$, these two relations cannot even distinguish $\nu$ from $-\nu$. They do not determine a physical sign convention for the deformation parameter, much less the full affine symmetry algebra.

## 3. Check that the derivative current permits the chamber calculation

Put $r=z_1-z_2$. Wick contraction gives the exact prefactor in the normal product

\[
J_B(z_1)J_B(z_2)
=r^\rho:
\left[
\left(C_1-\frac1r\right)\left(C_2+\frac1r\right)
+\frac1{r^2}
\right]e^{dX_1+dX_2}:.
\]

Thus

\[
J_B(z_1)J_B(z_2)
=r^\rho:\left[C_1C_2+\frac{C_1-C_2}{r}\right]
e^{dX_1+dX_2}:.
\tag{S9}
\]

The double pole cancels exactly, and the apparent remaining pole is a divided difference of a smooth normal field. The leading coefficient is $C^2+\partial C$. There is no $r^{\rho-2}$ term at zero mass.

There is also a direct non-nilpotence check, beyond the permutation coefficient. In a constant regular jet $C=c\ne0$, the coefficient of the corresponding normal-ordered exponential is $c^2r^\rho$. On a unit ordered interval its integral is $c^2/[(\rho+1)(\rho+2)]$, which is nonzero in the convergence strip. The full square multiplies it by $1+s$. Thus the present formal current is not a generically nilpotent odd screening. This does not exclude nilpotent actions in particular representations.

For three $B$ insertions and a smooth external background, write $F_i$ for the regular value of $C$ including contractions with separated spectators. The derivative Wick polynomial is

\[
W_3(F)=F_1F_2F_3
+\sum_{\{i,j,k\}=\{1,2,3\}\atop i<j}
F_k\frac{F_i-F_j}{z_i-z_j}
+2\sum_i\frac{F_i}{\prod_{j\ne i}(z_i-z_j)}.
\tag{S10}
\]

The last sum is the second divided difference. This expression stays bounded when the three $B$ positions coincide, provided the external insertions stay separated. It was obtained by summing all derivative/vertex contractions and derivative pairings, not by discarding poles in advance.

When an $A$ insertion at $a$ is also present, replace $F_i$ by $C_i+1/(z_i-a)$. The fully contracted three-$B$ polynomial vanishes if $C_i=0$. For a constant regular jet $C_i=c$, it is exactly

\[
W_3\left(c+\frac1{z_i-a}\right)
=c^3+c^2\sum_i\frac1{z_i-a}.
\tag{S11}
\]

Consequently the all-four-point collision has no degree $-3$ or $-2$ rational term. A first variation of $C_i$ can give degree $-1$, as can the $c^2$ term displayed above. The exponential factors in $BBBA$ have net degree zero, so the integrated relative degree of the full cluster is at least two. The $AAAB$ cluster has just one derivative prefactor and the same bound.

The proper collision bounds needed for the two quartics are:

| Subcluster | Integrated relative power after the Wick cancellations |
|---|---:|
| $AA$ or $BB$ | $1+\rho$ |
| $AB$ | $-\rho$ |
| $AAA$ or $BBB$ | $2+3\rho$ |
| $AAB$ or $BBA$ | at least $1-\rho$ |
| $AAAB$ or $BBBA$ | at least $2$ |

All are positive in $-1/2<\operatorname{Re}\rho<0$. The table uses relative dimensions; a common coordinate along the interval does not count as a small variable. Endpoint collisions are no worse. With separated spectators, the divided differences and the remaining Taylor coefficients are bounded on local patches. Subdividing by relative separations then gives absolutely integrable chamber integrals. Thus no contact distribution has been silently dropped in deriving (S8) in this open region.

For these ultraviolet products, subsequent meromorphic continuation of the complete chamber integrals preserves (S2)–(S3) at generic coupling away from poles. The antichiral calculation reverses $s$ to $s^{-1}$ and gives the same quartic coefficients. This establishes a formal ultraviolet result, not existence of global charge operators on a massive Hilbert space.

The specific representative in (S5) matters. Replacing $a_0$ by $a_2$ changes the current by a total derivative and a sign:

\[
J_B+\widetilde J_B=-t\partial:e^{dX}:.
\]

On a finite interval its integral is a boundary term. It must not be dropped during this calculation without a separate endpoint or twisted-cycle condition.

## 4. The first Yukawa correction already requires new product subtractions

The ultraviolet domain cannot simply be reused term by term after adding a tail. Direct contraction gives

\[
J_B(z)K_2(w)\sim (t-2)(z-w)^{\rho-2}\mathcal V_2(w)
\]

before its local cocycle and the factor $-t\lambda_2$ are inserted. Absolute integrability of this pair would require $\operatorname{Re}\rho>1$, whereas the ultraviolet $AB$ pair requires $\operatorname{Re}\rho<0$. There is no common domain for these two uncombined terms. A new analytic-continuation argument must first account for their product contacts.

### Combine the action insertion and both tails for a resolved pair

Let the two $B$ currents be at distinct positions $u,v$, and put $\beta=t-1$. Contract both derivative prefactors with the other exponential insertions and include their mutual derivative contraction. The lowest fused coefficient with $Y_2(w)$ is

\[
\begin{aligned}
&\left(-\frac1{u-v}+\frac\beta{u-w}\right)
\left(\frac1{u-v}+\frac\beta{v-w}\right)
+\frac1{(u-v)^2}\\
&\hspace{2cm}=\frac{\beta(\beta-1)}{(u-w)(v-w)}.
\end{aligned}
\tag{S12}
\]

The vertex factors add another $(u-w)^{-1}(v-w)^{-1}$ and $(u-v)^\rho$. In the source-first calculation, remove circular Cauchy regulators at fixed distinct $u,v$ before letting these current positions coincide. For real $u-v=r>0$,

\[
\int d^2w\,\frac1{(w-u)^2(w-v)^2}=\frac{2\pi}{r^2}.
\tag{S13}
\]

This is a finite-part area integral. To verify it independently, partial fractions give double coefficients $1/r^2$ and simple coefficients $-2/r^3,+2/r^3$. Green's theorem gives minus $\pi$ times the sum of the simple residues multiplied by their conjugate pole positions. With two circular holes of radius $\eta$, the exact correction is $-6\pi\eta^2/r^4$. A fresh Cartesian integration of an explicit antiderivative checks this formula at $\eta/r=0.1,0.05,0.025$; for the last value it gives $6.27140433473$ at $r=1$.

The action weight $-\lambda_2/\pi$ therefore gives the leading coefficient $-2\lambda_2(t-1)(t-2)r^{\rho-2}$. Each current-tail placement contributes $+t\lambda_2(t-2)r^{\rho-2}$. The signs include the cocycle: $J_BK_2$ has local sign $-1$, whereas $K_2J_B$ has sign $+1$ and the derivative contraction has the opposite orientation. Their sum is exactly (S4).

For an ordered spatial pair integrated over $r>\epsilon$, the leading divergent part is removed by adding

\[
\frac{2\lambda_2(t-2)}{\rho-1}\epsilon^{\rho-1}
\int dx\,\mathcal V_2(x)
=-2t\lambda_2\epsilon^{\rho-1}\int dx\,\mathcal V_2(x).
\tag{S14}
\]

The common charge normalization $(2\pi i)^{-2}$ is omitted from both sides of this ordered-integral formula. This is the **leading** product subtraction, not a complete prescription for every descendant or every ordering. The full charge square also carries the appropriate chamber weights. Exceptional parameters such as $\rho=1$ cannot be handled by setting a vanishing numerator to zero before continuation.

### Three $B$ insertions have a further proper collision

For $n=2,3$ derivative prefactors and one Yukawa, the complete scalar Wick polynomial is

\[
W_n=\frac{\beta(\beta-1)\cdots(\beta-n+1)}
{\prod_i(z_i-w)}.
\tag{S15}
\]

Let $\Delta(z)=\prod_{i<j}(z_i-z_j)$ and $D_i=\prod_{j\ne i}(z_i-z_j)^{-2}$. Integrating the source first, and including every tail placement, gives the leading real-axis density-product coefficient

\[
(-1)^n\lambda_2(t-2)(t-3)\cdots(t-n)
\Delta(z)^\rho\sum_iD_i,\qquad n=2,3.
\tag{S16}
\]

The sign follows from $(-1)^n$ in the current-product/Yukawa cocycle; the action and tail coefficients differ by $t-(t-1)=1$. For $n=3$ this is $-\lambda_2(t-2)(t-3)\Delta^\rho\sum_iD_i$, multiplying $\mathcal V_3$. It has raw relative degree $3\rho-2$. It must be included as a proper subcollision of $BBBA$; subtracting only the $BB$ pair is not enough.

For ordered positive gaps $u=z_1-z_2$, $v=z_2-z_3$,

\[
\sum_iD_i=\frac{2(u^2+uv+v^2)}{u^2v^2(u+v)^2}.
\]

Setting $R=u+v$, $x=u/R$ gives the angular integral

\[
\int_0^1[x(1-x)]^{\rho-2}(1-x+x^2)\,dx
=\frac{3\rho-1}{2(2\rho-1)}B(\rho-1,\rho-1),
\qquad \operatorname{Re}\rho>1.
\tag{S17}
\]

The radial power is $R^{3\rho-3}dR$. Independent quadrature checks (S17) at $\rho=\sqrt2$. Continuing this resolved three-current product is possible, but matching its complete subtraction, including its $BB$ subfaces, to the four-current contour prescription is a separate step.

## 5. Reduce the potentially marginal quartic contact

Use the reference ordering $B(z_1)B(z_2)B(z_3)A(a)$, with the corresponding branches. Other orderings are obtained by their explicit braid continuation. Put

\[
\mathcal R(z,a)=\prod_{i<j}(z_i-z_j)^\rho
\prod_i(z_i-a)^{-\rho},\qquad
\mathcal V=:e^{(2d+a_2)X+a_2\bar X}:.
\tag{S18}
\]

For three $B$ currents, an $A$ current and a Yukawa with derivative-contraction coefficient $\beta$, direct Wick summation gives

\[
W_3=\frac{\beta(\beta-1)}{\prod_i(z_i-w)}
\left[\beta+1+(a-w)\sum_i\frac1{z_i-a}\right].
\tag{S19}
\]

For $Y_0$, $\beta=1$, so this lowest scalar coefficient vanishes. For $Y_2$, $\beta=t-1$, it is not identically zero. For example, at $z_i=1,2,3$, $a=0$, $w=4$, it is $-\beta(\beta-1)(3\beta-19)/18$.

Multiplying the vertex factors, the lowest $Y_2$ kernel is

\[
\mathcal K_{BBBA,Y_2}
=(t-1)(t-2)\mathcal R(z,a)
\frac{(a-w)\left[t+(a-w)\sum_i(z_i-a)^{-1}\right]}
{\prod_i(z_i-w)^2}.
\tag{S20}
\]

Its homogeneity is $-5$. There are three relative contour coordinates and two source coordinates, so its raw integrated degree is zero, independently of $t$. Increasing $\operatorname{Re}t$ cannot make this primitive contribution vanish by power counting. This is a separate issue from the previously audited mixed two-charge collision.

The special momentum relation $\rho t=2$ gives a stronger identity:

\[
\boxed{
\mathcal K_{BBBA,Y_2}
=\frac{(t-1)(t-2)}\rho\,
\partial_a\left[
\frac{(a-w)^2\mathcal R(z,a)}{\prod_i(z_i-w)^2}
\right].
}
\tag{S21}
\]

Indeed $\partial_a\log\mathcal R=\rho\sum_i(z_i-a)^{-1}$, and differentiation of $(a-w)^2$ gives $2/\rho=t$. This is an exact identity of the leading Wick kernel before integration. It is not inferred from the permutation relation.

For a tail at $z_i$, let

\[
S_i^{(2)}=(a-z_i)^2D_i,\qquad S=\sum_i S_i^{(2)}.
\]

Including its coefficient and cocycle, the corresponding leading tail kernel is

\[
\frac{t\lambda_2(t-2)}\rho\,
\partial_a\bigl[\mathcal R S_i^{(2)}\bigr].
\tag{S22}
\]

All three placements were checked independently from their two-derivative Wick products. For real current positions, (S13)'s residue argument generalizes to

\[
\int d^2w\,\frac{(a-w)^2}{\prod_i(w-z_i)^2}=\pi S.
\]

Consequently the source-first combined leading action/tail kernel is

\[
\boxed{
\frac{\lambda_2(t-2)}\rho\,
\partial_a\bigl[\mathcal R S\bigr]\mathcal V.
}
\tag{S23}
\]

The derivative acts on the displayed coefficient. Differentiating a smooth fused field or spectator instead introduces a higher Taylor order and changes the primitive collision degree. Equations (S21)–(S23) reduce the leading scalar question; they are not full operator-product formulas for all descendants.

### Keep contour heights and cutoff boundaries

The real-axis formula is not the complete simultaneous-contour limit. For complex current positions write

\[
G(w)=\frac{(a-w)^2}{\prod_i(w-z_i)^2},\qquad
R_i=\operatorname{Res}_{w=z_i}G(w)=\partial_{z_i}S_i^{(2)}.
\]

The resolved-point Cauchy area integral is

\[
\int d^2w\,G(w)
=-\pi\sum_i\bar z_iR_i
=\pi S+2\pi i\sum_i(\operatorname{Im}z_i)R_i.
\tag{S24}
\]

The last equality follows from $\sum_i z_iR_i=-S$, obtained by taking the residues of $wG(w)$. At fixed separated positions the height term vanishes as the contours approach the real axis. When all separations and heights shrink together it has the same degree as $S$ and cannot be omitted. The script checks (S24) with real and complex current positions.

The combined leading expression therefore contains, in addition to (S23),

\[
-\frac{2i\lambda_2(t-1)(t-2)}\rho\,
\partial_a\left[\mathcal R\sum_i(\operatorname{Im}z_i)R_i\right]\mathcal V.
\tag{S25}
\]

It is still an $a$ derivative. But integrating this derivative produces the boundaries of the actual regulated domain. These include spatial endpoints, faces associated with the $BB$ and $BBB$ subtractions, and derivatives of any interaction exclusion whose center moves with $a$. The previous two-current mixed prescription does not by itself define all these four-current subtraction faces.

In particular, a source-first resolved-point formula cannot simply replace a limit in which current heights and interaction holes shrink at a fixed ratio. Nor can a total derivative be declared zero on a finite interval. The calculation required next is the sum of these boundary contributions with the four quartic word coefficients and all proper subtraction terms. That sum has not been evaluated here.

## 6. What is now known, and the next calculation

The ultraviolet part of the same-chirality test is concrete: the candidate even currents satisfy both quartic relations in a specified analytic product prescription. Their derivative dressing was checked rather than suppressed. The first interaction calculation supplies the leading $BB$ and $BBB$ product kernels and an exact primitive for the potentially marginal $BBBA/Y_2$ term.

The remaining first-order task is now narrower than a fresh search for an algebra:

1. Complete the same-chirality $BB$ and $BBB$ insertion subtractions, including any descendant and endpoint terms required by the chosen contour ordering.
2. Carry (S21)–(S25) through the common four-current regulator, differentiate its exclusion boundaries, and sum the four quartic words. Keep the canonical/product normalization separate from the action normalization.
3. If the first-order boundary sum vanishes, repeat the test at mass order two with the composite-potential and current counterterms. A nonzero result would require interpretation, not an adjustment chosen to recover a preferred algebra.

The $AAAB$ first-Yukawa full cluster has positive raw relative degree six and does not share the marginal scalar channel in (S20). The difficult new channel is concentrated in $BBBA$ and its proper multiple-$B$ subcollisions. This observation is not being substituted for a completed massive charge identity in either direction.

Only after these product relations and the previously identified mixed boundary operators have been handled can one test a complete quantum-group presentation. A faithful representation, a coproduct on the physical sectors, and massive spatial boundary values are not consequences of the two ultraviolet identities. No affine $\mathfrak{sl}(2|1)$ or affine $\mathfrak{sl}_2$ symmetry is assumed or established.

## 7. Checks, recommendations, and integration

Run

~~~text
python3 -u referee_scratch/sl21_same_chirality_20261023/chiral_checks.py
~~~

The final output is **61 PASS, 0 FAIL**. The checks include generic permutation ranks and extracted nullspaces, exact derivative-current cancellations, a nonzero constant-jet square, actual Yukawa Wick polynomials, independent area quadrature, the leading ordered-pair subtraction, the triple-product angular integral, all critical derivative identities, and the contour-height term. They do not purport to compute the missing four-current boundary sum.

No existing book or report assertion of a massive Serre relation was found to correct. For `briefs/referee_sl21_dual_screening_report.md`, line 77, add the following qualification if its integrated-current equivalence is reused for finite intervals:

> The replacement of $a_0$ by $a_2$ changes the integrated current by an explicit endpoint term, since $J_B+\widetilde J_B=-t\partial e^{dX}$. Equivalence requires a boundary or twisted-cycle condition and must not be assumed in a finite-endpoint charge-product calculation.

This is a recommended clarification, not an edit to that report. The earlier mixed-product results are not altered by the new same-chirality calculation.

## Integration

No book integration was performed. **Verified independently:** the two ultraviolet quartic relations, the derivative-current cancellations needed for their integral proof, the leading combined first-Yukawa $BB/BBB$ product kernels, and the exact-derivative reduction of the marginal four-current kernel. **Rejected:** promoting a separated braid matrix or a vanishing permutation coefficient directly to a massive Serre relation without collision analysis. **Still unclear:** the complete first-order four-current boundary sum, higher mass orders of the Serre test, physical boundary closure, and the symmetry algebra. All pre-existing files remain unchanged.
