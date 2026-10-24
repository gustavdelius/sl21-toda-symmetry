# Symmetric five-current extension and second-order descent: independent audit

Date: 2026-10-24. Model: Codex GPT-6. Reviewed commit: `ae19df62aaa82879da0f4600b74982d700d481b5`.

This is a new referee report. All previous files and outputs were preserved. Work was confined to this report and `referee_scratch/sl21_symmetric_descent_20261024/audit/`. No book or index edits, git mutation, worktree, commit, push, render or subagent was used. Primary formulas were reconstructed from the explicit momenta and covariance. No old research script was imported or rerun as evidence.

## Verdict and exact scope

**Correct, at the bounded formal scope and with the strengthened whole-jet normalization below:** the once-applied balanced interaction map preserves the two tested Serre relations with one extra barred current through mass order two. This covers

`S_B barB`, `S_B barA`, `S_A barB`, `S_A barA`,

and the phase-aligned placements obtained by moving that one barred constituent across the quartet. Here S_B and S_A are the canonical quartic combinations, with the same smooth species test repeated in their three equal-species slots. The statement is initially on a generic open patch `2<t<3`, `g=rho=2/t`, `h=t`, and continues only where the same complete families have regular restrictions.

This conclusion uses an explicit finite coefficient formula for the map, exact full-field kernel identities, uniform source-height bounds, and matching of every proper forest type at the two-source order. It does not follow from positive face degrees alone. The result is descent through the **tested mixed/Serre relations in this bounded family**, not a proof of descent through every possible identity of a completed Gaussian algebra. It does not construct a closed algebra with arbitrary numbers of barred blocks, an interacting Hilbert-space representation, or all-order symmetry.

There is one real correction to the preceding tentative prescription. The old constant-coupling endpoint multiplied by an arbitrary spectator is insufficient. The finite balance must cancel the complete independent first-jet functional. Otherwise a simple nonconstant source test gives a nonzero residual. The present positive verdict is for that explicitly completed normalization. It is not an assertion that the older endpoint data uniquely specified it.

| Item | Verdict | Independent evidence |
|---|---|---|
| Constant-spectator lifting of the finite mixed contact handles arbitrary smooth coefficients | wrong | Explicit source-jet counterexample in section 1 |
| The full one-source Laurent coefficients retain all derivative spectators | correct | Independent H3 expansion, H2/J coefficient check, and even-root companion |
| The unlike-source antichiral triangle is a genuine obstruction | wrong | Exact full-field reduction to one-leaf graphs times A_h and B_h |
| The nonsmooth source multiplier can be treated as smooth directly at physical h | wrong | Its differentiability depends on h; the proof first takes h sufficiently large |
| The full balanced two-source forest has the finite coefficient formula in section 3 | correct | Proper-face comparison, source assignments, full Taylor coefficients and vanishing residual unions |
| Every second-order Serre sector in the bounded family vanishes | correct with stated normalization | Decorated shuffle, marginal primitives, source-height estimates and continuation |
| A Serre null may be used in the initial domain g>6, rho<-1 | wrong | The null is on g=rho; the two domains have separate roles |
| These computations prove an unrestricted associative interacting algebra | unclear / not established | The tested family is finite and has only one barred block |

## 1. The contact lift that the proof actually requires

Use

`a0=(ib/sqrt(2), b/sqrt(2), i)`, `a2=(ib/sqrt(2), -b/sqrt(2), -i)`,

`t=b²`, `d=-(a0+a2)/t`, `rho=2/t`, `beta=t-1`.

Thus `ai²=-1`, `a0.a2=1-t`, `d.ai=1`. Let `B=:C exp(dX):`, `C=a0.partial X`, `A=:exp(-dX):`. Put `q0=1`, `q2=beta`, and `nu=exp(-i pi/t)`.

For a mixed B/barB/Yi face, keep an independent smooth coefficient

`psi=psi0+psi_z z+psi_w w+psi_barw barw+O(distance²)`.

After contraction with every separated spectator, denote the remaining derivative-current coefficients by C_eff and Cbar_eff. The actual local pair kernel is

`[C_eff/(z-w)+qi/(z-w)²][qi/barw²-Cbar_eff/barw] psi`.

With the prescribed circle, outer height `sigma delta`, `delta>epsilon`, and rotational Cauchy double pole, the three source moments are `pi sigma`, `i pi sigma`, `i pi sigma`. Directly integrating the outer current first gives

`C_raw,sigma[psi]=-i pi sigma nu lambda_i L_i[psi]`,

`L_i[psi]=qi²(psi_w+psi_barw)-qi(C_eff+Cbar_eff)psi0`.

The coefficient of psi_z vanishes separately in each ordering. In particular an outer slot-test derivative cannot be silently identified with a source or normal-field derivative.

For physical Y2 field jets, L_2 reduces to `t beta partial_x U`. For Y0 the physical field contribution vanishes. But a source switching derivative adds `beta² chi' U` in node 2 and `chi' U0` in node 0. For example, choose `psi0=psi_z=psi_barw=0`, `psi_w=1`. The constant-spectator endpoint lift gives zero, whereas the actual node-2 contact is

`-i pi sigma nu lambda2 beta²`,

which is nonzero generically. This is an actual counterexample to that shortcut, not evidence against existence of an interaction map.

The completed normalization used below is **minus this entire ordered functional**, with no additional even primitive moment in the initial convergent domain. One may instead choose a source switching function identically one near all current supports; the independent switching terms then disappear, but derivatives of other spectator coefficients must still be retained.

For the even-root A/barA/V1 face the corresponding functional is

`L_A[psi]=psi_w+psi_barw`.

On physical normal jets it equals `-gamma partial_x F_A^end`, `gamma=t/(1-t)`, where `F_A^end=:exp[-k(X+Xbar)]:` and `k=d+a0+a2`. The raw ordered contribution is therefore `+i pi sigma nu lambda1 gamma partial_x F_A^end`, canceled by the opposite balanced half. This endpoint differs from the current tail denoted F_A below, whose exponents are `-kX-T Xbar`. Its outer-test first moment again vanishes. The fresh scripts check this independently.

## 2. Full proper faces, not scalar residues

For three B positions z1,z2,z3, one A at a, and one source w, set

`Fj=Cj-sum_(k!=j)1/(zj-zk)+1/(zj-a)+qi/(zj-w)`,

`H3=F1F2F3+F1/(z2-z3)²+F2/(z1-z3)²+F3/(z1-z2)²`.

This is the complete derivative-field Wick polynomial: it contains the empty derivative matching and all three pairings. Fusing z1 with w leaves the full K_i spectator cluster, including two derivative B fields, the minority current, normal fields and inherited labels. Its double coefficient is `qi K_i^full`.

To check the simple coefficient independently, define

`fz=Cz-1/(z-v)+1/(z-a)+(qi-1)/(z-w)`,

`fv=Cv+1/(z-v)+1/(v-a)+(qi-1)/(v-w)`,

`H2=fz fv+1/(z-v)²`, `J=fv/(w-z)²+fz/(w-v)²`.

Direct expansion of H3 gives

`[C+qi D+(qi g-1)(1/(w-z)+1/(w-v))+(1-qi rho)/(w-a)]H2+(1-qi)J`,

where `D=d.partial X(w)`. Independently, `partial_w H2=(qi-1)J`. With `sigma0=1`, `sigma2=-1`, the vector and scalar identities

`C+qi D=sigma_i(d+ai).partial X`,

`sigma_i+(qi-sigma_i)rho=1`,

`sigma_i(qi-1)=1-qi`

therefore give the complete simple coefficient

`sigma_i partial_w K_i^full +(qi-sigma_i)(g-rho)[1/(w-z)+1/(w-v)]K_i^full`.

The term involving J is essential. This derivation verifies the coefficient of every uncontracted C field, not merely the scalar term. Derivatives act on the full coefficient and on tests when distributions are smeared.

The even-root counterpart has `sigma_A=1/(1-t)`. Its complete remaining derivative-current factor is

`f=C+(1-t)/(z-w)+1/(z-a2)+1/(z-a3)`.

The direct simple coefficient is

`[-D+g/(w-a2)+g/(w-a3)-rho/(w-z)]f+1/(z-w)²`.

The fused derivative has normal coefficient `(t-1)D`, majority-edge coefficients `g-2`, minority coefficient `2-rho`, and derivative of f equal to `(1-t)/(z-w)²`. Checking these four independent coefficients proves the same full Laurent formula, with remainder `(1-sigma_A)(g-rho)` on the majority edges.

These identities supply the actual unary faces in the extension comparison. They are not replaced by an assumption about the desired algebra.

## 3. The finite second-order forest coefficient

Define A_i to be one **spatial-base** action insertion of Yi, including its action sign and fixed common source test; A_ij is the two-source insertion before equal-source factorials. Let A_V denote the even-root insertion. Define slot substitutions only on original current constituents:

`U0: B -> K0, barB -> barK0`,

`U2: B -> beta K2, barB -> beta barK2`,

`U1: A -> -F_A, barA -> -barF_A`.

Here `Ki=:exp[(d+ai)X+ai Xbar]:`, and F_A is the current tail `:exp[-kX-T Xbar]:`, with `T=a0+a2`, `k=d+T`. The substitutions sum over original slots and act as zero on an already-tail slot at this order. All scalar contractions, source labels, field descendants, tests and cocycles remain attached. In particular a current has only one C field; two different sources cannot both consume it.

Let E_bar mean the once-applied phase-aligned half-sum of the two extreme positions of the one barred block, preserving the order of its chiral subword. It is not an iterated Jordan multiplication.

For a complete input cluster F, the unlike coefficient is

`E_bar [ A02(F)+A0(U2 F)+A2(U0 F)+U0U2 F ]`.                 (1)

U0U2 here sums distinct original-slot assignments. The equal-source coefficient is

`E_bar [ (1/2) Aii(F)+Ai(Ui F)+(1/2) Ui² F ]`.              (2)

The even-root coefficient is

`E_bar [ A_V(F)+U1 F ]`.                                   (3)

These formulas are followed by the **same full radial Taylor continuation**, not separate finite parts of intermediate binary blocks. They apply also when F already contains a canonical input tail; only its still-original slots can receive a unary replacement. Thus all one-tail/one-source and two-tail input terms are included when F is expanded in its couplings.

The combinatorics can be checked without any analytic assumption. Expand the source series once and the slot substitution once. The distinct-source term is A02 plus one source assigned to a unary slot and the other to the bulk in both orders, plus assignments to two distinct slots. For equal sources, the original 1/2! cancels the two choices in the single-unary term and remains in Ui². Fresh symbolic expansions for one through four consuming slots give exactly (1)–(2). Dressing each binary intermediate word separately would not give these coefficients.

### Why these are the actual completed forests

The comparison is first made in an off-diagonal initial domain, for example `Re g>6`, `Re rho<-1`, and sufficiently large Re h. No Serre null is used there.

1. Away from faces the kernels are the same ordinary full Wick expressions. On an isolated current/source face,

   `Pf(1/z²)=-partial_x(1/z)+pi delta²(z)`

   gives the exact unary substitution. The complete Laurent calculation in section 2 fixes its spectator jets.

2. On disjoint faces the product of independent delta substitutions is legitimate and gives exactly the distinct-slot terms of (1)–(2). The same-source, two-current faces are not disjoint. Their delta product would contain a squared source-time delta and is undefined. They instead receive the full mixed functional of section 1. Its opposite ordered half is precisely the selected extreme-order repair. This applies to the full coefficient containing the second source, not to a frozen endpoint.

3. Source-pair faces use complete Taylor coefficients. For example direct differentiation gives

   `[1/(z-r)+(t-1)/z]/[(z-r)z]=sum_(n>=0)(n+t)r^n/z^(n+3)`.

   The factor n+t includes the descendant's contraction with the external C. Taylor extraction commutes with external Wick contraction while spectators are separated. The uniform estimate in section 4 extends this comparison to the required nested faces. Radial integration supplies denominators h+n, and the K2/Y0 and U/Y0 faces supply h+n and h-1+n with their complete field/test jets. This includes the potential P, current L and endpoint representatives once. They are not additional independent insertions to add a second time to the same forest.

4. After those proper faces are matched, no unmatched finite distribution is left on a higher union in the initial domain. For a mixed face, the first-jet-subtracted inner remainder has degree +1. If its scale is epsilon and the next cluster scale is R, its bound is

   `epsilon integral_epsilon^1 R^(Re omega-2) dR`.

   This tends to zero, including the borderline omega=1 case `epsilon log(1/epsilon)`. Comparable-scale unions contribute `O(epsilon^Re omega)`. Every derivative of a mask costs one scale power and is paired with one fewer pole in the differentiated primitive, so the degree is unchanged. The product rule includes derivatives on distinct circle masks as well as twice on the same mask. Their intersections are covered by these sectors.

For two sources and a current-containing cluster with a A, b B, d barred B, the **actual** degrees are

`E+3a-2(b+d)+2h+1` for unlike sources,

`E+3a-2(b+d)+5` for equal sources,

where `E=g[binomial(a,2)+binomial(b,2)]-rho ab`. Every higher two-source union is positive in the comparison domain after the fixed proper faces; the minimum equal-source margin is +1. Source-only degrees are respectively 2h and 4. Thus the complete masked distributions and (1)–(3) agree there and, consequently, as their common meromorphic families.

The construction worker's frozen section 7 uses 2a in a quantity called a conservative degree. It is a lower bound obtained by dropping the nonnegative a contribution; the actual degree is 3a as independently counted above. Its initial-domain positivity argument remains valid. It should not be mistaken for the exact degree in the physical-face test.

## 4. Resolving the unlike-source triangle and the uniform bound

Put `p=ubar-barw0`, `q=ubar-barw2`, `r=w0-w2`, so `bar r=q-p`. The entire antichiral derivative-current factor, including the source-pair scalar, satisfies

`|r|^(2h-2)/(pq) [Cbar+1/p+beta/q]`

`=A_h[Cbar(1/p-1/q)+1/p²-beta/q²]`

` +(beta-1)B_h(1/p-1/q)`,                                  (4)

where

`A_h=r |r|^(2h-4)`, `B_h=r² |r|^(2h-6)`.

This follows from `(q-p)/(pq)=1/p-1/q` and a second partial fraction. I checked the complete rational identity independently. The uncontracted Cbar term is retained. Every term on the right has only **one** antichiral edge. The remainder of the graph is holomorphic and directed acyclic in each height-order region. All h dependence is now in A_h and B_h, so the order of the leaf-graph distribution can be chosen independently of h.

For equal sources, with qi=1 or beta,

`(q-p)/(pq)[Cbar+qi/p+qi/q]`

`=Cbar(1/p-1/q)+qi(1/p²-1/q²)`.

The remaining holomorphic source-pair factor is polynomial. This identity, rather than a count of an unreduced triangle, controls the equal-source full five-current degree `3(g-rho)`, which is marginal at the physical diagonal.

Here is the needed uniform argument. On a compact generic patch of (g,t), fix a finite distributional order N for the finitely many holomorphic-plus-leaf boundary kernels. Their Fourier-Schwinger edge parameters are nonnegative. A discrete topological potential bounds their total size by a constant times the incidence image. This proves proper pushforward and uniform seminorm bounds as any nonnegative height gaps close. The constant depends on the finite graph order, not the numerical gaps. Continuing individual normalized edge powers by distributional integration by parts gives the same statement on compact exponent patches.

For the unlike term, choose `2 Re h-4>N` with a positive margin. Derivatives through order N of A_h and B_h extend continuously across r=0. If a source-pair circle of radius epsilon is inserted, every Leibniz term in the derivative of its cutoff error is bounded by

`const * epsilon^(2 Re h-4-N)`.

Indeed a cutoff derivative of order j contributes epsilon^-j and the remaining multiplier derivative contributes epsilon^(2 Re h-4-(N-j)). This bound includes moving-circle derivatives and their intersections. It proves that the source-pair coefficient can legitimately be used as the smooth multiplier needed by the leaf estimate. Treating it as smooth directly at physical h would have been invalid.

There are only finitely many regions as the two source heights pass the current heights. The leaf bound is uniform on all of them. Regions between collapsing quartet heights have source-height measure tending to zero, so their contributions tend to zero in a fixed test seminorm. On the complement one takes the quartet limit first with the source heights separated. Spatial differentiation never differentiates the source heights, so it creates no extra source-time delta. The rotational deltas have already been isolated as the Ui terms.

The same proof covers a one-tail/one-source term. A hidden source in a tail retains its source label; fusing one source position in (4) gives the required A_h or B_h multiplier. For a source meeting a barred K2 tail, the unlike factor is `r |r|^(2h-4)` and the equal factor is polynomial. Thus this is not an assumption that a tail behaves like a newly unlabelled vertex.

## 5. The null calculation on the physical current diagonal

Now set `g=rho`, keeping t in a compact generic subinterval of `2<t<3`, and initially keep h sufficiently large as in section 4. This is a different step from the off-diagonal extension comparison. The Serre polynomial is **not** asserted zero in the comparison domain.

Write `c_s=1+s+s^-1`, `s=exp(-i pi rho)`. The two quartics are

`S_B=B³A-c_s B²AB+c_s BAB²-AB³`,

`S_A=A³B-c_s A²BA+c_s ABA²-BA³`.

The required decorated nulls are the coefficients obtained by substituting arbitrary constant multiples of K0, K2 and F_A into the original current slots, at the allowed mass orders. They can be verified on full spatial Wick kernels, not merely on a scalar expectation.

* K0 and K2 have the same current-label exchange phases as B. The mixed pair K0/K2 has holomorphic exponent `g+h-3` and antiholomorphic exponent `h-1`. Their difference is `g-2`, so h cancels from the exchange phase. The source cocycles supply the integer signs. The analogous single-tail phases are unchanged.
* A fresh enumeration with three **distinct decorated majority letters** verifies the physical shuffle coefficient for every one of the 24 spatial orders. Each is exactly zero. This verifies the polarized coefficients needed for U0, U2, Ui² and U0U2, retaining which slot carries which tail. Full Wick polynomials and smooth normal/test coefficients are common to the terms after spatial reordering.
* For the one-K0 marginal channel, `-a0.(d+a0)=0`, and its two remaining derivative-current scalar is exactly

  `[-1/(z-v)+1/(z-a)][1/(z-v)+1/(v-a)]+1/(z-v)²=0`.

  Every surviving regular-field or test jet improves the degree.
* The one-K2 marginal scalar has its independent minority-coordinate primitive. With `eta=t-2`, `R=product_(i<j)(zi-zj)^g product_i(zi-a)^(-rho)` and tail at c, the complete scalar is

  `(eta/rho) partial_a[R (a-c)²/((c-z)²(c-v)²)]`.

  I rechecked this after dividing by R, including the derivative-pair term. Integration by parts onto an arbitrary smooth coefficient improves the marginal degree by one; proper-face continuation fixes the endpoint terms from the initial domain.
* The one-F_A marginal scalar similarly obeys

  `S_i [sum_j 1/(z-aj)-t/(z-ai)]=-(1/rho) R_A^-1 partial_z(R_A S_i)`,

  where `S_i=(z-ai)²/product_(j!=i)(aj-ai)²`. This identity was independently checked.
* Two equal tails have positive full quartet degree `3(g-rho)+2`. For one K0 and one K2 the full degree is `3(g-rho)+2h-2`. Its new proper degrees include `g+2h-3`, `3g+2h-5` and `g-2rho+2h`. They are positive in the physical patch and at large h. Hence the decorated shuffle zero does not multiply a hidden marginal pole in these sectors.

Consequently each relevant quartet substitution is a zero distribution, including all normal and test jets, before the extra barred factor is appended in the common boundary family.

The uniform leaf estimates now show that every Ai, Aij term in (1)–(3) maps the corresponding quartet null to zero. At fixed separated source heights this follows by grouping the quartet first; the shrinking source-height regions vanish by the uniform seminorm bound. With a pre-existing tail the same proof uses the inherited-label multiplier described after (4). The decorated null is valid with that external smooth source factor: after spatial sorting, the factor is the same for every term with the same decorated placement. No off-diagonal Serre identity is being assumed.

The Ui terms vanish by the decorated quartet identities themselves. If a unary is placed on the extra barred slot, the remaining quartet is its original null; if it is placed inside the quartet, the corresponding first or second variation is zero. The once-applied extreme average is linear, so it also gives zero. This proves the **entire** sum (1)–(3), not merely its spatial-source component.

In particular the complete order-two inventory is exhausted:

| Coupling sector | Contributions |
|---|---|
| lambda0², lambda2² | equal double action, one action/one unary, two distinct unaries, input one-tail/action and one-tail/unary, input two tails |
| lambda0 lambda2 | unlike double action, both bulk/unary assignments, distinct U0/U2 assignments, input tail/action and tail/unary, common P/L/endpoint Taylor forests |
| lambda1 | one V1 action, A/barA unary F_A replacements, canonical F_A input tails and their balanced AA face |

Every source is counted once. There is no illegal same-source delta square and no extra source expansion on a binary intermediate product.

## 6. Restriction in h and the companion barred families

At large h the preceding null is established on g=rho. Continue that zero meromorphically in h using the **same whole family**, with its already fixed proper-face extensions. To restrict to h=t one must still check regularity; identity of a meromorphic zero is not a license to take unrelated finite parts.

For a barred B the new unlike-face degree is

`E+3a-2(b+d)+2h+1`.

For a barred A it is

`E+3a-2b+3d+2h+1`.

The antichiral Yukawa factors at barred A are polynomial, so the h multiplier there is simply `|r|^(2h-2)` times those polynomials, easier than (4). Independently enumerating every current-containing subset in each of the two quartic families gives positive degrees at `g=rho=2/t`, `h=t`, `2<t<3`. Multiplication by t and writing `t=2+x` gives polynomials with nonnegative coefficients and a strictly positive interior value. The full S_B barB degree is `3(g-rho)+2h-4`; after restriction it is `2t-4>0`.

The remaining old one-source faces are the already fixed complete one-leaf families. Their apparently marginal full degree is controlled by the full boundary distribution, not by assigning independent finite parts to spatial chambers. Equal-source faces also have the exact leaf decomposition above. Nested Taylor levels improve the new degrees. The explicit source-pair/current/endpoint counterterm denominators are h+n and h-1+n; their cutoff representatives vanish in the large-h domain and supply the unique continued radial forest, with no finite independent P coupling. Their joint positive physical degrees prevent an extra evanescent finite remainder on these new faces.

For V1 and barred A, direct counting gives

`omega=E-a+2b-d+1`.

Thus the dangerous S_A barA/V1 full face is `3(g-rho)-1`. Its exact kernel is

`R_g (z-w)² product_i(ai-w)^(-2) (ubar-barw)^(-2)`

` times [C(z)+sum_i 1/(z-ai)-t/(z-w)] N`.

It is holomorphic plus one double antichiral leaf. The even-root Laurent and full mixed-jet computations in sections 1–2 apply, and the same uniform one-source proof handles the apparent first-jet marginality. V1 with barred B has polynomial antichiral factors and is easier. This verifies both additional barred choices rather than extrapolating from S_B barB.

For a separately written mixed endpoint identity, retain its off-surface derivative-regulator remainders. In particular dressing a physical derivative endpoint by an unlike source can produce an explicit `(h-t)` term. Such an identity is not asserted zero at arbitrary large h. Its physical restriction follows from the fixed full-jet face equality and regularity of the same continued family. This is distinct from the quartet proof, whose decorated shuffle is already zero on g=rho for large h.

The result is therefore the desired bounded second-order descent. It does not say that every null in an unspecified completion has been classified. It also does not justify arbitrary numbers of barred constituents: the single-leaf reduction and the extreme-position normalization are special to the present family.

## 7. Evidence, diagnostics and preservation

Fresh successful outputs contain **252 PASS records**:

* `check_descent.out`: 110, including the exact triangle, vector-derived powers, physical face tests, marginal primitives and all 24 fully decorated spatial shuffles.
* `check_companion_v2.out`: 104, covering the barred-A/even-root channels and their actual degrees.
* `check_proper_face_v2.out`: 22, including the independent full H3 Laurent coefficients and source-fusion jets through n=6; two records mark the written mask-remainder bounds retained in the argument.
* `check_full_jets.out`: four, plus the explicit nonzero counterexample to the constant-spectator shortcut.
* `check_forest_inventory_v2.out`: eight, checking the finite source/unary coefficients for one through four consuming slots.
* `check_even_laurent.out`: four independent coefficients of the even-root proper face.

Three implementation diagnostics were retained: an uncancelled rational expression prevented one initial SymPy polynomial conversion; a full rational H3 simplification was interrupted and replaced by its four coefficient checks; and an ambiguous SymPy symbol range produced an empty slot list before a forest test. Corrected scripts were written as new files, and the original scripts/partial outputs and explanations remain. No mathematical FAIL was hidden. The source-jet shortcut has a deliberately retained nonzero residual.

The four snapshotted primary inputs remain byte-identical. No integration into the book occurred. The construction worker's completed coefficient comparison was read after the independent kernel calculations and checked against them; it was not used as computation evidence. Its separation of extension comparison from quotient annihilation is essential and is preserved here.

Recommended bounded wording:

> With the declared full-first-jet balanced normalization, one common once-applied interaction map reproduces the prescribed symmetric five-current extension and preserves the tested mixed/Serre relations through mass order two for either Serre quartet and one additional barred generator. The proof retains all source assignments, descendants and moving-mask intersections, first establishes the extension comparison off the current diagonal, and then proves and restricts the null identities on the physical current diagonal. It does not establish an unrestricted interacting quantum affine algebra.

The report and its completed outputs are now frozen. The final synthesis was inspected for scope and notation; the current-tail/endpoint distinction and the general two-majority-species degree were flagged for correction before its freeze.
