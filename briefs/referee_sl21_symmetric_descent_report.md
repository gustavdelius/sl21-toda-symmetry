# Symmetric five-current extension and second-order descent

Date: 2026-10-24. Model: Codex (GPT-6). Reviewed commit: `ae19df62aaa82879da0f4600b74982d700d481b5`, branch `main`.

This is a fresh referee-only continuation. Two independent workers reconstruct the symmetric extension and audit quotient descent in separate scratch directories. Earlier reports are evidence to check, not premises accepted without calculation. No book text, book check script, index, or earlier report is edited. This report records a construction in a bounded formal sector; it makes no claim about a physical operator realization or an all-order algebra.

## Outcome and scope

The symmetric five-current comparison can be completed, provided the mixed contact is lifted to a functional of the **entire first source jet**. Multiplying the old endpoint derivative by a spectator coefficient is insufficient. The required functional is explicit below, including the even-root companion and the source-coefficient derivative that survives in the node-0 channel.

The second-order test is organized as an actual coefficient identity for the once-applied map. Its nontrivial analytic ingredient is an exact reduction of the unlike-source antichiral triangle to single-leaf graphs with sufficiently differentiable scalar multipliers. This permits the two source integrations to commute with the grouped Serre limit in a common domain. Proper-face matching, rather than power counting alone, identifies that construction with the prescribed moving-mask forests.

The resulting verdict is **correct for the tested five-current extension and second-order mixed/Serre descent with this full-jet normalization**. This closes the two bounded tests left by the previous report. It does not prove descent on every Gaussian null relation or construct the unrestricted interacting algebra.

The domain tested here contains both Serre quartets, each with at most one additional barred generator in any position, through mass order two. Repeated generators in a Serre polynomial have the same slot test, or the corresponding polarized sum is taken. The product kernels accept arbitrary independent smooth tests; the unpolarized Serre polynomial with unrelated repeated-slot tests is not asserted to vanish. Products with two independently barred generators and arbitrary null elements of the entire Gaussian algebra remain outside the test.

The physical proof is on a generic open patch `2<t<3`, with `g=rho=2/t` and `h=t`; continuation beyond that patch is asserted only where these same complete families have regular restrictions.

## 1. Conventions and the map being tested

Let

`a0=(ib/sqrt(2),b/sqrt(2),i)`, `a2=(ib/sqrt(2),−b/sqrt(2),−i)`,

`t=b²`, `T=a0+a2`, `d=−T/t`, `rho=2/t`, `beta=t−1`, `gamma=t/(1−t)`, `v=d+a2`, `k=d+T`.

The scalar products are `ai²=−1`, `a0.a2=1−t`, `d.ai=1`, `d²=−rho`. Write

`C=a0.partial X`, `B=:C exp(d.X):`, `A=:exp(−d.X):`,

`Yi=:exp(ai.(X+Xbar)):`, `V1=:exp(−T.(X+Xbar)):`,

`Ki=:exp((d+ai).X+ai.Xbar):`, `K_A=:exp(−k.X−T.Xbar):`,

`U=:exp(v.(X+Xbar)):`, `F_A^end=:exp(−k.(X+Xbar)):`.

K_A is the even-root current tail; F_A^end is the mixed AA endpoint. They are different fields. The barred tail is `:exp(−T.X−k.Xbar):`. This notation keeps the current substitution distinct from the endpoint contact throughout the calculation.

All currents are even. Put `nu=exp(−i pi/t)` and `s=nu²`. The cocycle on constituent labels is

`epsilon(n,m)=(-1)^((nL+nR)(m0+m2)) exp[−i pi(nL mR−nR mL)/t]`.

The auxiliary exponent g replaces like-current rho powers. The unlike-source factor is `|w0−w2|^(2h−2)`. Original constituent labels, g factors and source ancestry survive fusion. Physical restriction means `g=rho,h=t`, after extension, at generic nonresonant parameters. The canonical input currents are `B−t lambda2 K2`, `A−gamma lambda1 K_A` and their barred partners. Couplings lambda0 and lambda2 have mass order one; lambda1 has mass order two.

For a composed Gaussian cluster F, the map inserts the action once:

`I(F)=F−sum_i lambda_i/pi integral R_bal W(F,Yi)`

`−lambda1/pi integral R_bal W(F,V1)`

`+sum_(i<=j) lambda_i lambda_j/[pi²(1+delta_ij)] integral integral R_bal W(F,Yi,Yj)+O(m³)`,

where i,j range over 0,2 and smooth compact source coefficients are understood. W is the actual Gaussian Wick contraction. R_bal retains rotational Cauchy finite parts on unary faces, applies the full mixed matching in section 2, and subtracts the complete radial Taylor jets on nested/disjoint source forests. No additional symmetric primitive local moment is added in the initial convergent domain. This last clause is an explicit normalization choice.

The source-pair action counterterm `kappa_h=−lambda0 lambda2 epsilon^(2h)/h` multiplying the composite potential is the forest representative, not a further independently inserted source. Likewise the cutoff current improvement `+(t/h)lambda0 lambda2 epsilon^(2h)L` is generated once. It must not be counted again if it is already present in an input representative.

The map is linear and generally nonunital as an operator insertion series. Normalized scalar expectation values are a separate construction. Applying I independently at every binary multiplication would be a different prescription.

## 2. The required contact refinement

**Verdict: correct for the full first-jet functional; wrong for the constant-spectator shortcut.**

Place barB at zero and the outer B at `z=x+i sigma delta`, where sigma is ±1. Let `delta=c epsilon`, c>1. Write the arbitrary smooth coefficient of the proper mixed face as

`psi=psi0+psi_z z+psi_w w+psi_barw wbar+O(distance²)`.

This coefficient includes source switching, slot tests, spectator contractions and the normal field. Let `q0=1`, `q2=beta`; let Ceff and Cbar_eff denote the remaining regular derivative-current factors after spectator contraction. The singular pair kernel is

`[Ceff/(z−w)+qi/(z−w)²][qi/wbar²−Cbar_eff/wbar] psi`.

The outer integrals are `−i pi sign(sigma delta−Im w)` for the simple pole and `pi delta(sigma delta−Im w)` for the double pole. Multiplication of the double pole by z gives the sum of the simple-pole integral and `pi w delta(sigma delta−Im w)`. The three remaining source moments are `pi sigma`, `i pi sigma`, `i pi sigma`, respectively. Consequently the complete raw local action moment is

`C_raw,sigma[psi]=−i pi sigma nu lambda_i L_i[psi]`,

**`L_i[psi]=qi²(psi_w+psi_barw)−qi(Ceff+Cbar_eff)psi0`.**

The coefficient of psi_z cancels exactly, before substituting a physical normal field. The balanced ordered addition is the negative of this entire functional. This is a local distributional rule; it is not a modification of the ordinary off-diagonal source integral.

For the physical normal jet, `beta² a2−beta a0=t beta v`, so the Y2 result is `t beta partial_x U`. For a source coefficient chi it is instead

`t beta chi partial_x U+beta²(partial_x chi)U`.

It is neither the endpoint formula multiplied by chi nor `t beta partial_x(chi U)`. At node 0 the physical normal-field term vanishes, but `(partial_x chi)U0` remains, where `U0=:exp((d+a0).(X+Xbar)):`. In AA/V1 the corresponding functional is simply

`L_A[psi]=psi_w+psi_barw`.

Since `−T=−gamma(−k)`, its physical normal part is `−gamma partial_x F_A^end`.

An explicit counterexample to a frozen spectator is `psi0=psi_z=psi_barw=0`, `psi_w=1`: the actual BB moment is `−i pi sigma nu lambda_i qi²`, whereas the shortcut gives zero. The construction worker independently obtains this counterexample. The coordinator derives the complete moment afresh in `coordinator/jet_functional.py`: **9 PASS**.

This refines, rather than changes, the constant-coupling physical endpoint in the preceding report. The passage at `briefs/referee_sl21_interaction_map_report.md:97–103` must not be used as a rule for arbitrary source-dependent coefficients. Proposed clarification, without editing that historical report:

> On a general mixed face, the balanced addition is the negative of the ordered first-jet functional L_i, or L_A in the even-root channel. The displayed derivative of U or the AA endpoint is its physical normal-field specialization. Differentiate the full source and spectator coefficient before taking the contact restriction.

## 3. Complete five-current kernels and moving masks

**Verdict: correct with the full-jet normalization.**

For B at z1,z2,z3, A at a, barB at u and Yi at w, define

`R_g=product_(j<k)(zj−zk)^g product_j(zj−a)^(-rho)`,

`Fj=Cj−sum_(k!=j)1/(zj−zk)+1/(zj−a)+qi/(zj−w)`,

`H3=F1 F2 F3+F1/(z2−z3)²+F2/(z1−z3)²+F3/(z1−z2)²`.

Apart from the fixed whole-word cocycle, the full kernel is

`R_g (a−w) H3 [Cbar(u)+qi/(ubar−wbar)] N /[product_j(zj−w)(ubar−wbar)]`.

N is the common normal exponential. H3 contains every derivative-field matching; discarding its paired terms would invalidate the first-jet calculation. The coordinator reconstructs H3 recursively from Wick matchings rather than taking the scalar specialization.

The dangerous even-root companion A³BbarA/V1 has kernel

`R_g (z−w)² [C(z)+sum_j1/(z−aj)−t/(z−w)] N /[product_j(aj−w)²(ubar−wbar)²]`.

The construction report derives the complete Laurent coefficients on every unary source/current face, including derivatives of spectator Wick polynomials and the evanescent `(g−rho)` terms. Those terms are retained until the complete boundary family is continued. They cannot be dropped in an individual spatial chamber.

The conversion from a spatial derivative to a rotational double pole uses

`Pf(1/z²)=−partial_x(1/z)+pi delta²(z)`.

With the actual action and cocycle signs the unary repairs are `+lambda0 K0`, `+beta lambda2 K2`, `−lambda1 K_A`, and their barred partners. They are terms of the map, distinct from the stipulated finite input tails. Their symmetric Pf/PV contributions survive; balanced matching does not set the symmetric distribution to zero.

For an arbitrary product M of smooth approximations to the moving masks, the comparison uses the full identity

`partial_x partial_y(MF)=M partial_x partial_y F+(partial_x M)partial_y F+(partial_y M)partial_x F+(partial_x partial_y M)F`.

The last term includes derivatives on two different masks. These are exactly the boundary-intersection terms that a one-chamber residue computation misses.

For a current-containing collision with a A slots, b B slots and d barred slots, the one-Yukawa relative degree in the dangerous family is

`omega_B=g binomial(b,2)−rho ab+2a−b−d+1`.

For V1 in the even-root family it is

`omega_A=g binomial(a,2)−rho ab−a+2b−d+1`.

Unary proper faces have degree zero; the prescribed mixed pair/source face has degree −1 and requires its full first jet. Every higher union has a strictly positive margin in `Re g>6, Re rho<−1`. The coordinator enumerates all current-containing subsets, including both other quartet/barred combinations and every source channel, directly from the primary vectors.

The comparison proof has three parts. The two constructions agree away from collisions. On proper faces they agree as functionals of the entire coefficient jet by section 2 and the unary Laurent calculation. After those subtractions, sector decomposition gives an inner remainder of positive degree. For an inner cutoff epsilon and outer cluster radius R, its nested contribution is bounded by a test seminorm times

`epsilon integral_epsilon^1 R^(Re omega−2) dR`.

This tends uniformly to zero in the stated initial domain; comparable-scale intersections are `O(epsilon^Re omega)`. Differentiating a mask costs an inverse scale, while the corresponding primitive gains that scale. Thus the same estimate controls the intersection terms. This argument covers nested and comparable cutoff regions, not only a fixed ratio of radii.

For the simultaneous prescription, preserve the order of the four chiral constituents and phase-align the five barred placements R0,...,R4. If `Delta_j=Rj−R_(j−1)`, the correction at placement k is

`(sum_(j>k)Delta_j−sum_(j<=k)Delta_j)/2`.

Every corrected placement equals `(R0+R4)/2`. The increments are the full jet contacts. This is one average of the two extreme full words, not a sequence of binary Jordan products and not a uniform average of all five placements.

Equality of the prescriptions is first established off the physical regulator diagonal in the convergent domain, then continued as a complete family. No Serre zero is asserted in that off-diagonal domain. Complete chiral boundary products are controlled by the acyclic Fourier-cone estimate in section 4; individual chamber residues need not vanish.

## 4. The two-source obstruction and its resolution

**Verdict: correct; an exact full-field identity supplies the missing uniformity.**

Let `p=ubar−wbar0`, `q=ubar−wbar2`, `r=w0−w2`, so `q−p=bar r`. The complete antichiral triangle satisfies

`|r|^(2h−2)/(pq) [Cbar+1/p+beta/q]`

`=A_h[Cbar(1/p−1/q)+1/p²−beta/q²]+(beta−1)B_h(1/p−1/q)`,

where

`A_h=r|r|^(2h−4)`, `B_h=r²|r|^(2h−6)`.

The identity includes the uncontracted Cbar term and both double poles. Each resulting graph has only one antichiral singular edge. Every h dependence is in its scalar coefficient. Equal-source triangles instead reduce using their polynomial Vandermonde factor. With barA the relevant antichiral source factors are polynomial already. The same reductions hold for one-source/one-tail terms when source ancestry is retained.

For any fixed compact parameter patch, the remaining finite collection of boundary graphs has some finite distributional order N independent of h. To see uniformity in the imaginary gaps, orient its chiral edges by height order. For positive edge momenta k and incidence image q, a strictly monotone integer potential c on the resulting directed acyclic graph satisfies

`sum c_i q_i=sum_(i→j)(c_i−c_j)k_ij >= sum k_ij`.

The positive-cone pushforward is proper. Normalized edge-power distributions and their integer derivative shifts therefore have smooth-test seminorm bounds uniform as height gaps shrink. Attaching one antichiral leaf creates no directed cycle.

Choose `2 Re h−4>N`. A_h and B_h are then C^N across r=0, and removing the source-pair cutoff has C^N error bounded by `epsilon^(2 Re h−4−N)`. This avoids a circular assumption that N is independent of an h-dependent graph: the partial fraction has extracted all that dependence explicitly.

One may consequently first take the chiral quartet boundary limit with the source heights separated. The finitely many source-height strips between collapsing current heights have widths tending to zero; the uniform seminorm estimate makes their integrals vanish. Intersections of two such strips have correspondingly vanishing area. Spatial differentiation does not differentiate a source height. Rotational contacts have already been accounted for by the unary and full mixed repairs.

This proves the required interchange for the spatial two-source operator and its one-tail/one-source companions. At the physical diagonal g=rho, decorated Serre identities can be used inside that bound. This is distinct from the off-diagonal domain used to prove equality of the extension prescriptions.

## 5. Second-order coefficient test

**Verdict: correct on the bounded tested relation sector.** With `c_s=1+s+s^(-1)`, the two polynomials are

`S_B=B³A−c_s B²AB+c_s BAB²−AB³`,

`S_A=A³B−c_s A²BA+c_s ABA²−BA³`.

The tested five-current zeros are their products with barA or barB, with the barred slot moved through the quartet using the specified full crossing contacts. Here the symbols in the physical polynomial include their stipulated finite tails. The test is linear over their required descendant and coupling coefficients. Calling this the tested relation sector does not assert closure under multiplication by an arbitrary further word.

Let A_i denote the spatial one-source operator, including its action sign, and A_ij the labeled two-source operator before the equal-source factorial. Let U_i be the unary replacement operator on the **original current slots**: U0 replaces B by K0, U2 replaces B by beta K2, and U1 replaces A by −K_A, with the corresponding barred replacements. Include the source coefficient in each replacement. U_i does not act again on an already replaced slot.

For unlike couplings the once-expanded coefficient is

**`A_02(F)+A_0(U_2F)+A_2(U_0F)+U_0U_2F`.**

For two equal couplings it is

**`(1/2)A_ii(F)+A_i(U_iF)+(1/2)U_i²F`.**

For the independent even-root coupling it is `A_1(F)+U_1F`. These expressions are subjected to the single extreme-order prescription and the same radial continuation. The canonical finite input tails are expanded in F at their own mass orders; thus the total second-order relation contains the two-source bare term, every one-source/one-input-tail term and every source-free double-tail term.

These formulas specify a finite forest comparison to check, rather than assuming that a map on words descends automatically. Distinct sources assigned to distinct current slots give U0U2, or U_i²/2. A source assigned to two currents is an overlapping proper face and uses the full L functional; it is not a product of delta distributions. A single derivative current cannot contract to two sources. Mixed contacts with the other source present differentiate the entire remaining coefficient, including its source-pair factor. The source-pair potential and current-improvement representatives are included once by the radial continuation.

The construction worker independently certifies this comparison in section 7 of its report. After matching these proper faces, the exact two-source degrees for a A, b B and d barred B constituents are

`omega_02=E+3a−2(b+d)+2h+1`,

`omega_ii=E+3a−2(b+d)+5`,

where `E=g[binomial(a,2)+binomial(b,2)]−rho ab`. For the B-majority quartet a is at most one, so its first binomial is zero.

They count four real source coordinates and remove the common real translation. Source-only faces are counted separately. The worker uses a conservative lower bound with 2a in place of 3a; the coordinator and auditor independently obtain the exact 3a coefficient from the primary vectors. The lower bound is sufficient in the initial domain. The least positive equal-source residual degree can be one, yielding `epsilon log(1/epsilon)` rather than O(epsilon); it still vanishes. Thus no unaccounted higher-union term remains in the comparison. Equality off the physical diagonal is then continued before the chiral null identities are used.

The source-pair descendant matching is also explicit. Expanding the exact derivative-current factor gives

`[1/(z−r)+(t−1)/z]/[(z−r)z]=sum_(n>=0)(n+t)r^n/z^(n+3)`.

The n+t factor is the sum of the differentiated geometric-series coefficient n+1 and beta. Radial pairing supplies denominators h+n. This expansion is taken after multiplication by the complete normal field and source/test coefficients. It fixes the same descendant jets in both forest representations, including the n=1 coefficient t+1 that a leading-scalar subtraction would miss.

The new mixed unary polarization is especially useful. The K0/K2 holomorphic and antiholomorphic exponents are `g+h−3` and `h−1`, respectively. Their exchange exponent is `g−2`, so the h dependence cancels and the separated shuffle phase is unchanged. The full K0K2BA relative degree is

`3(g−rho)+2h−2`.

Its new proper degrees are `g+2h−3`, `3g+2h−5`, `g−2rho+2h`. They are positive on the physical patch `2<t<3`. The equal-tail quartets have positive full degree `3(g−rho)+2` and polynomial right-moving factors.

The single K0 marginal scalar vanishes identically:

`[−1/(z1−z2)+1/(z1−a)][1/(z1−z2)+1/(z2−a)]+1/(z1−z2)²=0`.

The K2 and K_A marginal scalars are exact minority-coordinate derivatives. Integration by parts moves a derivative onto the full smooth coefficient and improves the degree. These identities are independently reconstructed in the adversarial check. Together with the fully polarized four-letter shuffle, they establish the bare, unary and double-unary chiral null sectors. The uniform bounds of section 4 then apply A_i or A_ij to these zero relations. The even-root source has only a single antichiral double leaf and obeys the same bound. This covers both quartets with barA or barB.

For the full mixed endpoint relation, an off-physical h-dependent regulator derivative can leave a factor h−t. It must be retained until the complete properly extended family is restricted. A physical zero is not inferred merely from that factor; regularity at the restriction is also required. The positive unlike-source degrees and the entire complete chiral boundary family provide the regularity test in the generic physical patch. This distinction is recorded in the independent audit.

## 6. Evidence and limitations

| Item | Verdict | Independent evidence |
|---|---|---|
| Full mixed source-jet contact and AA companion | Correct | Three independent moment calculations; explicit frozen-spectator counterexample |
| Constant-spectator lift for arbitrary source coefficients | Wrong | Nonzero psi_w test in section 2 |
| Symmetric five-current moving-mask extension | Correct in the stated normalization | Complete Wick kernels, proper-face jets, mask-intersection identities and nested remainder bounds |
| Two-source coefficient formula | Correct through mass order two | Independent forest allocation and proper-face comparison, including P/L descendants once |
| Both Serre quartets with either single barred generator | Correct through mass order two | Decorated shuffle and scalar primitives, uniform source integration and physical continuation |
| Arbitrary multi-barred quotient or all-order associative interacting algebra | Unclear / not established | Outside the enumerated domain; the required broader descent test is absent |

The independent construction report is `briefs/referee_sl21_symmetric_descent_construction_report.md`; the independent adversarial report is `briefs/referee_sl21_symmetric_descent_audit_report.md`. They distinguish the local extension comparison from algebraic quotient annihilation. The construction evidence contains 87 completed PASS records and its explicit shortcut counterexample; the audit contains 252 PASS records across six completed outputs. Partial runs and diagnostic failures are excluded from those completed counts and preserved. The final audit also caught a notational conflation of the even-root tail and endpoint in the synthesis draft; section 1 now defines them separately, with the actual tail used in every current substitution. The coordinator's fresh checks are under `referee_scratch/sl21_symmetric_descent_20261024/coordinator/`:

| Check | Result | What it verifies |
|---|---|---|
| `checks.py` | 69 PASS | Primary vectors, full derivative Wick matching, arbitrary first jets, mask intersections, union degrees, triangle and mixed-tail exponents |
| `companions.py` | 365 PASS | All four quartet/bar combinations, every source channel through mass order two, proper-face inventory and initial-domain positivity |
| `jet_functional.py` | 9 PASS | Complete raw first-jet moment, outer-test cancellation, switching correction and AA companion |

Outputs are retained beside the scripts. The initial general-check output retains one symbolic-comparison failure fixed by comparing expanded differences; the initial companion output retains the partial run before a rational expression was canceled prior to polynomial conversion. Neither repair changed a mathematical formula. The construction worker retains an explicitly expected failure for the incorrect frozen-spectator lift. Numerical or symbolic PASS counts supplement the written distributional arguments; they do not substitute for them.

No theorem label is recommended for the book from these formal identities alone. No claim of novelty or literature correction is made. No evidence here constructs products with arbitrary numbers of barred currents, establishes every Gaussian null relation, identifies a physical Hilbert-space representation, or proves convergence of the perturbative series. In particular, a transported-product formula without the corresponding broader descent proof would not close the remaining full-algebra gap.

## Integration

No integration into the book is performed or requested in this referee-only pass. The previously displayed physical constant-source contact is verified. Its lift to general coefficients is made explicit and strengthened to the whole first-jet functional. The symmetric moving-mask comparison is verified using that functional. Second-order descent is assessed only on the explicitly enumerated mixed/Serre sector, not on an unrestricted quotient. Earlier reports remain unchanged so that the refinement and its reason remain visible.

## Preservation

The final read-only verification confirms that all **4,029 snapshotted pre-existing files** are byte-identical. HEAD, branch, tracked diff and staged diff match their initial states. Evidence is in `referee_scratch/sl21_symmetric_descent_20261024/preservation.json` and `preservation.out`. The fresh reports and scratch evidence are listed with hashes in that directory's `output_manifest.json`. No commit, push, worktree operation or rendering command was run.
