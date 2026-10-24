# Symmetric five-current moving masks: construction referee

Date: 2026-10-24. Model: Codex GPT-6. Reviewed commit: `ae19df62aaa82879da0f4600b74982d700d481b5`.

This is a fresh referee-only report. Earlier reports and working outputs remain immutable. Only this report and the assigned `referee_scratch/sl21_symmetric_descent_20261024/construction/` directory are written. No book, index, Git state or rendered output is changed. Definitions from earlier reports are used as definitions; the calculations below are newly reconstructed.

## Verdicts and scope

1. **Correct:** the complete one-source five-current Wick kernel has an explicit finite matching formula. Its Laurent coefficients retain both derivative-field spectators and their first jets. The same calculation works for the dangerous even-root companion `A³B barA / V1`.
2. **Wrong:** lifting the old mixed contact merely as `partial_x U` times an otherwise constant spectator coefficient handles arbitrary source/spectator jets. A fresh independent-jet calculation gives a nonzero counterexample. The replacement is an explicit first-jet functional, stated below. It includes a node-0 switching contact even though that node's constant-coupling physical endpoint vanishes.
3. **Correct with the stated whole-jet normalization:** the prescribed nested rotational moving masks and the simultaneous extreme-order construction give the same complete one-source five-current distribution. The comparison uses the full proper-face moments, an explicit mask product rule, and uniform estimates for every residual union. It is not an argument from one radial chamber.
4. **Correct:** the unlike two-source antichiral triangle admits a full-field partial fraction into single-leaf terms. Arbitrarily high source-pair differentiability is available in a large-Re-h domain. This removes the particular nonsmooth-coefficient obstacle to the one-source comparison.
5. **Unclear in this construction report:** full mass-order-two quotient annihilation. The comparison of extension procedures and the vanishing of a complete Serre ideal under the resulting forest map are different statements. This report supplies the former local mechanism and identifies the exact remaining quotient equality; an independent worker tests the latter. I do not replace that equality by positivity counts or by a definition of the map on words.

The selected symmetric normalization is an additional choice: apply the negative of each ordered **whole-jet** contact, and add no further even primitive local moment in the initial convergent domain. Reflection symmetry alone would not force that choice.

## 1. Complete one-source matching

Use `t=b²`, `rho=2/t`, `beta=t−1`, `d=−(a0+a2)/t`, `C=a0.partial X`. Thus `d.ai=1`, `ai²=−1`, `d²=−rho`, and `a0+a2=−td`. The currents are `B=:C exp(d.X):`, `A=:exp(−d.X):`; the sources are `Yi=:exp(ai.(X+Xbar)):`. The auxiliary g belongs to the original current constituents and remains attached after fusion. All cocycles below are the cocycle of the entire labelled word; they must not be replaced by a pair cocycle when other current constituents are present.

Place the three B constituents at `z1,z2,z3`, A at `a`, the barred B at `u`, and Yi at w. Write `q0=1`, `q2=beta`, and

`R_g = product_(j<k)(zj−zk)^g product_j(zj−a)^(-rho)`,

`Fj = Cj − sum_(k!=j)1/(zj−zk) + 1/(zj−a) + qi/(zj−w)`,

`H3 = F1 F2 F3 + F1/(z2−z3)^2 + F2/(z1−z3)^2 + F3/(z1−z2)^2`.

The full Gaussian kernel, apart from its fixed word cocycle, is

`R_g (a−w)/[product_j(zj−w)(ubar−wbar)] H3 [Cbar(u)+qi/(ubar−wbar)] N`,

where N is the common normal exponential of all five currents and the source. Each C in H3 is an uncontracted derivative field inside that common normal product. H3 enumerates the empty matching and all three possible derivative-field pairs. There are no missing double pairs among three derivatives.

This expression, including N, is the input to every mask and forest operation. The source variable has at most a double pole at any B constituent and at the barred constituent. Partial fractions therefore reduce the mask differentiation to at most one holomorphic and one antiholomorphic differentiated simple-pole factor per term. Counting all four derivative fields as four independent pole differentiations at one source would be wrong.

### The exact Laurent coefficient

Fuse z1 with w and denote the entire resulting K_i spectator cluster by `K_i^full`; this includes its two remaining B derivative fields, A, barred B, scalar edges, normal exponential and all inherited labels. Let `sigma0=1`, `sigma2=−1`. The double-pole coefficient is `qi K_i^full`. The simple-pole coefficient is

`Si = sigma_i partial_w K_i^full + (qi−sigma_i)(g−rho) [1/(w−z2)+1/(w−z3)] K_i^full`.

Here partial_w differentiates the entire fused field and its spectator contractions, not an endpoint exponential in isolation. The second term must be retained before continuation.

For an independent verification, set

`f2=C2+(qi−1)/(z2−w)−1/(z2−z3)+1/(z2−a)`,

`f3=C3+(qi−1)/(z3−w)−1/(z3−z2)+1/(z3−a)`,

`H2=f2 f3+1/(z2−z3)^2`, `J=f3/(w−z2)^2+f2/(w−z3)^2`, and `D=d.partial X(w)`.

Directly expanding the four terms of H3 gives the simple coefficient, after removing its common scalar and normal factor,

`[C(w)+qi D+(qi g−1)sum 1/(w−zj)+(1−qi rho)/(w−a)] H2 +(1−qi)J`.

On the other hand `partial_w H2=(qi−1)J`. The normal-field derivative of K_i is `D+ai.partial X`. The identities

`sigma_i(D+ai.partial X)=C+qi D`,

`sigma_i+(qi−sigma_i)rho=1`,

`sigma_i(qi−1)=1−qi`

prove the displayed formula. They are checked independently for both sources. In particular this calculation retains the derivative of the spectator matching itself; a scalar-only check would miss J.

### The even-root companion

For three A at `a1,a2,a3`, one B at z, barred A at u and `V1=:exp(−T.(X+Xbar)):` with `T=a0+a2`, the complete kernel is

`R_g product_j(aj−w)^(-2) (z−w)^2/(ubar−wbar)^2 [C(z)+sum_j1/(z−aj)−t/(z−w)] N`.

At a1=w its double coefficient is the complete F_A fused cluster. Its simple coefficient is

`S_A = sigma_A partial_w F_A^full +(1−sigma_A)(g−rho) sum_(j=2,3)F_A^full/(w−aj)`,

` sigma_A=1/(1−t)`.

To verify it, the remaining derivative-current factor is `f=C+(1−t)/(z−w)+1/(z−a2)+1/(z−a3)`. The direct simple coefficient is

`[−D+g/(w−a2)+g/(w−a3)−rho/(w−z)]f+1/(z−w)^2`.

Differentiating the fused cluster gives normal derivative `(t−1)D`, scalar derivative `(g−2)sum1/(w−aj)+(2−rho)/(w−z)`, and `partial_w f=(1−t)/(z−w)^2`. Substitution proves the claim, including the B spectator contraction.

## 2. The independent-jet mixed face

This is the necessary correction to a tempting shortcut. Set the inner barred current to zero, the outer current to `z=x+i sigma delta`, and exclude `|w|<=epsilon`, where `delta=c epsilon`, `c>1`. The Cauchy double pole is the rotational Wirtinger extension, not a spatial derivative. Write an arbitrary smooth coefficient to first order as

`psi=psi0+psi_z z+psi_w w+psi_barw wbar+O(relative_distance²)`.

The derivative-current coefficients after contraction with all distant spectators are `Ceff,Cbar_eff`. The pair kernel to this order is

`[Ceff/(z−w)+qi/(z−w)^2][qi/wbar²−Cbar_eff/wbar] psi`.

The outer elementary integrals, in the vanishing smoothing limit, are

`integral dx 1/(z−w)=−i pi sign(sigma delta−Im w)`,

`integral dx Pf 1/(z−w)^2=pi delta(sigma delta−Im w)`,

`integral dx z Pf 1/(z−w)^2=−i pi sign(sigma delta−Im w)+pi w delta(sigma delta−Im w)`.

The remaining circle/edge moments are respectively `pi sigma`, `i pi sigma`, `i pi sigma` for the sign divided by wbar², the edge w/wbar², and the edge 1/wbar. For the first one the horizontal slice is `2 sqrt(epsilon²−y²)/epsilon²` inside the small circle; the contour height lies outside that circle. The other two are elementary Cauchy integrals along the edge.

Including the action and pair cocycle nu, the complete local ordered moment is

`C_raw,sigma[psi] = −i pi sigma nu lambda_i L_i[psi]`,

`L_i[psi] = qi²(psi_w+psi_barw)−qi(Ceff+Cbar_eff)psi0`.

The coefficient of psi_z is **zero**, separately in each ordering. This proves the arbitrary outer-slot-test cancellation without identifying a test derivative with a field derivative.

For the physical Y2 normal jet, this reduces to `t beta partial_x U`. For Y0 the physical field part cancels. With a nonconstant source switching chi, however, the additional terms are `beta²(partial_w+partial_barw)chi U` for Y2 and `(partial_w+partial_barw)chi U0` for Y0. Choosing chi identically one near all current supports removes these switching terms but does not license dropping the source dependence of spectator contractions.

For AA/V1 there are no derivative-current factors in the mixed pair. The corresponding functional is simply

`L_A[psi]=psi_w+psi_barw`.

On the physical normal jet it is `−gamma partial_x F_A`, where `gamma=t/(1−t)`. This agrees with the old constant-coupling AA endpoint, while retaining its independent source/test jets.

The generalized balanced assignment is `−C_raw,sigma[psi]` on each ordered face. This is the explicit lift used in the remainder of this report. The literal prescription `t beta partial_x U` multiplied by an arbitrary spectator coefficient is wrong: take `psi0=psi_z=psi_barw=0`, `psi_w=1`. Its value is zero, but the actual ordered moment is `−i pi sigma nu lambda2 beta²`. The fresh script retains this expected failure rather than silently dropping that test.

This functional also matches the full unary repair. The two unary source allocations pull a smooth coefficient back to `psi(r,r,r)` and `psi(r,0,0)`. Their first derivatives differ by `psi_w+psi_barw`; their zeroth values agree. Thus the double-pole boundary jumps cancel their leading delta-prime terms and all outer-test derivatives, leaving exactly the source jet in L_i. Their simple-pole jumps give the `Ceff+Cbar_eff` part. The same statement holds with the complete spectator Wick polynomial as the coefficient. This is why the extreme-order average, interpreted on full clusters, supplies the required symmetric lift.

## 3. Actual moving-mask comparison

Use smooth approximants to all source circles, current exclusions and contour strips before taking their specified nested limits. Denote their product by M. For a term with one holomorphic and one antiholomorphic simple-pole primitive F, the exact comparison includes

`partial_x partial_y(MF) = M partial_x partial_y F +(partial_x M)partial_y F +(partial_y M)partial_x F +(partial_x partial_y M)F`.

If `M=product_gamma theta_gamma`, then

`partial_x partial_y M = sum_gamma (partial_x partial_y theta_gamma)product_other theta + sum_(gamma!=eta)(partial_x theta_gamma)(partial_y theta_eta)product_other theta`.

The last sum is essential. It is the contribution from intersecting moving boundaries. It is explicitly retained in the comparison; a separate-chamber residue does not account for it.

Partition relative configuration space into finitely many nested cluster sectors. Away from collision faces the two kernels are identical. On a unary source/current face, the preceding full Laurent coefficients give the same rotational Pf extension, including the finite correction relative to spatial differentiation. On a mixed face, the independently calculated whole-jet L functional gives the same zeroth/first moment and the same ordered balance. These equalities are identities of coefficient fields and their jets. Hence a spectator is not frozen when a smaller face is subtracted.

For one source, a cluster with a A constituents, b B constituents and d barred B constituents has relative degree

`omega_B(a,b,d)=g binomial(b,2)−rho a b+2a−b−d+1`,

with `a<=1,b<=3,d<=1`. The fixed unary faces have degree zero; the B/barB/source face has degree −1 and needs its full first jet. Every higher union containing that mixed face has degree, respectively,

`g−2`, `3g−3`, `1−rho`, `g−2rho`, `3(g−rho)−1`.

All are greater than one in `Re g>6, Re rho<−1`. The remaining one-source faces and bare-current collision faces are integrable there as well. For the even-root companion the same computation gives

`omega_A(a,b,d)=g binomial(a,2)−rho a b−a+2b−d+1`,

with `a<=3,b<=1,d<=1`. The fixed A/barA/V1 face again has degree −1; all higher unions have the same positive first-jet margin in this domain. The script enumerates all admissible current-containing faces, not only the total collision.

After a mixed first-jet subtraction the inner remainder improves from degree −1 to degree +1. If its scale is epsilon and the next cluster scale is R, the absolute remainder is bounded by a finite smooth-test seminorm times

`epsilon integral_epsilon^1 R^(Re omega−2) dR`.

It is O(epsilon) uniformly in the stated initial domain. Comparable-scale intersections are O(epsilon^Re omega). Unary zeroth-jet remainders have the same positive inner degree. Each derivative of a mask costs one inverse cluster scale, but the associated primitive has one fewer pole; the displayed omega is unchanged. Thus all the differentiated-mask terms, including intersections, satisfy the same bounds. A finite partition into nested sectors covers unequal cutoff ratios; no ratio-dependent remainder survives.

Consequently the difference between the two complete one-source constructions vanishes on every compact smooth test in the initial domain. This is a comparison proof: equality off the diagonals, equality of the actual full proper-face moments, and a vanishing bound for everything left. Merely counting positive total degrees would not have established the proper-face equalities.

For the simultaneous construction, preserve the chiral order and move the one barred constituent across all four chiral slots. Let the phase-aligned ordered values be `R0,...,R4`, with actual crossing increments `Delta_j=Rj−R_(j−1)`. Give placement k the correction

`1/2 [sum_(j>k)Delta_j−sum_(j<=k)Delta_j]`.

Each corrected placement equals `(R0+R4)/2`. The Delta_j here are the entire jet functionals just computed, including their original labels; they are not constant endpoint symbols. This is a once-applied average of the two extreme full words. It is not an iterated Jordan product or an equal average of all five placements.

Finally continue the common complete family in g. Its one-source boundary kernels decompose into chiral directed acyclic graphs plus one antichiral leaf. The Fourier-Schwinger edge parameters are nonnegative. A strictly increasing potential on the finite vertex order bounds the sum of edge parameters by a constant times their incidence image, so their pushforward is proper; the normalized edge-power family is entire and its seminorm bounds are uniform as height gaps shrink. The complete ordered family therefore has no pole created by selecting a single chamber at g=rho. Equality established in the initial domain continues to the generic physical restriction. The argument keeps the evanescent `(g−rho)` Laurent terms until this complete-family step.

## 4. The two-source angular issue and what it does prove

For unlike sources put `p=ubar−wbar0`, `q=ubar−wbar2`, `r=w0−w2`, so `q−p=bar r`. The complete anti factor admits the exact identity

`|r|^(2h−2)/(pq) [Cbar+1/p+beta/q]`

`= A_h [Cbar(1/p−1/q)+1/p²−beta/q²] +(beta−1) B_h(1/p−1/q)`,

where `A_h=r|r|^(2h−4)`, `B_h=r²|r|^(2h−6)`.

I independently check the rational identity after dividing by `|r|^(2h−2)`. It includes the uncontracted Cbar field. Both coefficients have homogeneity at least `2 Re h−4`, including their angular dependence. Every remaining graph term has only one anti leaf; h now occurs only in its scalar multiplier. Thus the finite distributional order N of the graph can be chosen independently of h. For `2 Re h−4>N`, these coefficients are sufficiently differentiable at r=0, and source-circle cutoff derivative errors are bounded by `epsilon^(2 Re h−4−N)`. This justifies applying the complete one-source moment comparison with the second source as a coefficient in this domain. Using a nonsmooth |r| power as though it were an arbitrary smooth spectator at physical h without this step would be invalid.

Equal-source anti factors reduce by the ordinary Vandermonde partial fraction, and do not require an h multiplier. The full two-source Wick polynomial is linear in the contractions of any one C field: a single C cannot consume two sources. Conversely two distinct derivative currents can consume distinct sources. Those disjoint unary faces must be included. Two currents consuming the same source are overlapping faces, not a product of deltas; their union uses the generalized mixed functional above, dressed by the other source.

The radial source-pair forest must also retain its full Taylor coefficients. Independently expanding the exact rational factor gives

`[1/(z−r)+(t−1)/z]/[(z−r)z] = sum_(n>=0)(n+t)r^n/z^(n+3)`.

The coefficients n=0,1,2,3 are checked afresh. They include the derivative-current contraction at every descendant level. Radial integration supplies denominators `h+n`; multiplication by normal fields, tests and switching functions is done before selecting each Taylor coefficient. This is the concrete reason that a forest with only scalar leading subtraction would fail the comparison.

These facts permit propagation of the one-source **extension comparison** through the admissible nested/disjoint two-source forests in a high-Re-h domain: the proper mixed-face functional is continuous on the needed coefficient jets, source-pair coefficients have sufficient differentiability, and the same-source overlap is treated as a union instead of an illegal delta product. They do not by themselves show that the resulting forest sum annihilates every Serre relation.

## 5. The remaining quotient statement

Let S_B and S_A denote the canonical quartic Serre combinations with the stipulated repeated-species slot tests. The exact missing assertion for mass order two is

`[lambda_i lambda_j] I_bal(S_B barB)=0`,

and its even-root/other-barred companions, with I_bal applied once to the composed Gaussian cluster. The coefficient contains the spatial two-source integral, both unary repairs, the disjoint double-unary terms, generalized mixed-face repairs with the other source present, source-pair radial forests, and the one-tail/one-source terms from the specified finite input currents. A sum of positive collision degrees proves that no extra primitive local term is forced in an initial domain; it does not prove this displayed algebraic/distributional sum is zero.

In particular, the extreme-order construction needs the corresponding pure-chiral forest variation to vanish before its average is a map on the quotient. A source-integrated zero with fixed separated source heights can pass to the limit only after the uniform leaf/graph estimates above; that is now available for the spatial piece. The repaired pieces still need to be identified with the same complete differentiated identity, or explicitly summed. This report does not silently infer that step from associativity of the underlying Gaussian multiplication.

The new all-jet balanced functional is useful even if that last test fails: it fixes a previously unspecified and testable part of the single interaction map. It also gives a precise counterexample to an insufficient contact lift without claiming nonexistence of a coherent interacting product.

## 6. Checks and retained evidence

The new scratch script `check_symmetric.py` prints **49 PASS** checks. It verifies the full H3 simple pole for both Yukawas, the even-root simple pole, mask product and intersection rules, all admissible higher one-source union degrees at an interior point (with the analytic inequalities supplied above), all five balanced placements, the outer-test cancellation, the full unlike anti triangle and four independent source-pair Taylor coefficients. Its output is retained in `check_symmetric.out`.

The separate `check_independent_jets.py` prints **5 PASS** checks and **one explicitly expected FAIL** for the constant-spectator shortcut. It reconstructs the independent source/current first-jet moment, checks its unary allocation interpretation, and checks the AA/V1 functional. The failure output is retained in `check_independent_jets.out`.

These scripts import no previous research script and rerun no old evidence. They are algebraic checks of the formulas above; the distributional comparison is the explicit sector and remainder argument in section 3. No theorem label or book replacement is recommended on the basis of the uncompleted full mass-order-two quotient test alone.

## 7. Final coefficient comparison, with quotient scope kept separate

The local comparison can be completed at mass order two without assuming a Serre null in the off-diagonal initial domain. This section certifies equality of two concrete extension formulas; the decorated quartet null calculations remain an independent input.

Let `S_i` be the once-inserted **spatial-base** source integral, with its action sign, and `S_ij` the corresponding twice-inserted integral before factorial weights. Let `U0` replace one **original** B or barred B slot by K0 or barred K0, and `U2` replace it by beta K2 or beta barred K2. Let U1 replace one original A or barred A slot by minus F_A or minus barred F_A. They sum over original slots, and act as zero on a slot that is already a tail. These rules encode actual pole matchings: a derivative current has only one derivative field to contract. Let U02 sum the assignments of U0 and U2 to two distinct original slots, once per source assignment. Preserve the original labels, ordering and common cocycle in every term.

For a full labelled cluster F, the unlike coefficient of the balanced rotational insertion is

`E_bar [ S02(F) + S0(U2 F) + S2(U0 F) + U02 F ]`,

where `E_bar` is the once-applied average of the two extreme barred positions. For an equal source the coefficient is

`E_bar [ (1/2) Sii(F) + Si(Ui F) + (1/2) Ui² F ]`.

U_i² here means ordered distinct original-slot replacements; its factor one-half supplies the correct unordered coefficient. The V1 coefficient is `E_bar[S1(F)+U1 F]`. The same expansion applies to a pre-existing finite-tail term in the input, with U acting only on its still-original current slots. It does not add another source to each binary intermediate product.

Here is the face-by-face justification, rather than just formal coefficient counting.

* Away from source/current faces, the spatial and rotational kernels agree as ordinary Wick contractions. Their difference on one source/current double pole is the explicit unary delta from `Pf(1/z²)=−partial_x(1/z)+pi delta²(z)`. Taking that contraction in the complete Wick polynomial gives exactly `S_i(U_j F)` with the other source still integrated. It differentiates or evaluates the complete source/test coefficient, not a scalar leading term.
* On two disjoint source/current faces, the two independent delta replacements give U02, or the half-weighted Ui². The unused edge factors and full normal field are evaluated at those source allocations. No contribution exists that consumes two sources in a single C field. The script independently expands five original slots and confirms both statements.
* When two current faces consume the same source, they overlap. Their union receives the whole-jet mixed moment of section 2. The extreme average supplies its opposite ordered half correction, including the source coefficient derivatives contributed by the second source. No product of the two overlapping deltas is taken.
* For unlike sources the second-source coefficient on this mixed face is justified by the A_h/B_h reduction and large-Re-h differentiability in section 4. Equal sources have the polynomial Vandermonde reduction. The sole potentially marginal equal-source higher union has degree +1, so its nested error is `epsilon log(1/epsilon)`, still vanishing.
* A source-pair face and its descendant current face use the same full radial Taylor expansion on both sides. In a sufficiently large-Re-h initial domain their shrinking-circle representatives tend to zero; their continued representatives are the single P counterterm and single L improvement with the inherited descendants. Thus P/L are included once by this forest continuation, or once as explicit cutoff input, never twice. This step does not equate independent subtractions by fiat: it uses equality of the full Taylor coefficients such as `(n+t)` derived above.

For completeness, a two-source cluster containing a A, b B and d barred B has conservative relative degree

`omega_02 = g binomial(b,2)−rho a b+2a−2b−2d+2h+1`

for unlike sources, and

`omega_ii = g binomial(b,2)−rho a b+2a−2b−2d+5`

for equal sources. These formulas count four real source variables and remove one real translation when a current is present. The source-only pair is treated separately by its ordinary radial degree. In the previous large-g, negative-rho domain and sufficiently large h, every higher current-containing two-source union is positive; the equal-source minimum is +1. A direct enumeration checks all 30 cases. After the proper unary and full mixed jets have been matched, the same sector estimates as section 3 make every residual union vanish. The formulas therefore agree as complete distributions in an off-diagonal initial domain, and hence as their jointly continued families. This certifies the finite forest coefficient formula.

It does **not** assert that the Serre combination itself vanishes in that off-diagonal domain. The canonical null identities are on `g=rho`; the initial comparison domain has independent g and rho. Restrict the proved extension identity to the physical surface only after the complete boundary families have been continued and their restrictions justified. To obtain quotient descent, one must then use the independently checked decorated quartet nulls and uniform spatial-source insertion bounds on that surface. In particular a proof that secretly used the quartic null at `Re g>6, Re rho<−1` would be invalid.

The additional scratch check first printed 30 PASS face checks, then encountered a SymPy symbol-range naming error before its combinatoric section. Its original file and partial output are retained together with the failure explanation. A separate fresh script using unambiguous symbol names prints three PASS checks for distinct-slot assignments and the correctly ordered source/repair coefficient. No failed run was deleted. Together with the earlier checks, there are **87 completed PASS checks**, one explicitly expected mathematical-shortcut FAIL, and one retained script naming error repaired in a separate script.
