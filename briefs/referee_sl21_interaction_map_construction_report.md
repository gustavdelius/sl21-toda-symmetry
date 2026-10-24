# A single interaction map with explicit symmetric contacts

Date: 2026-10-24. Model: Codex (GPT-6). Reviewed checkout commit: `ae19df62aaa82879da0f4600b74982d700d481b5`.

This is an independent worker report in the editor-authorized interaction-map continuation. All earlier outputs are immutable. No subagents, book edits, existing-report edits, index edits, git mutations, worktrees, commits, pushes or renders were used. Fresh evidence is in `referee_scratch/sl21_interaction_map_20261024/construction/`. No old research script was imported or rerun.

## Outcome

A concrete single insertion algorithm is specified below on **composed labelled Gaussian clusters**, through mass order two. It has explicit unary source/current contacts, an explicit balanced ordered mixed contact, and radial source-pair/current/endpoint forest coefficients. It supplies nonzero symmetric terms rather than declaring them absent. It is applied once to the composed cluster; it is not an independently dressed binary multiplication.

The mixed pair and first one-barred triple pass a direct first-order regrouping test, including field and test jets. The low second-order forests reproduce the source-pair, current-improvement and endpoint coefficients. The possible K0 contribution to the Serre variation has a zero marginal scalar, so it does not give the proposed obstruction.

**The full quotient/descent assertion remains unclear.** The formula defines a linear map on resolved labelled kernels with a declared minimal symmetric normalization. Its agreement with every complete prescribed interacting word, and its descent through all canonical distributional identities at mass order two, are not proved here. In particular, a forest formula alone is not a proof of that assertion. The earlier defect prescription did not specify a unique allocation of ordered symmetric contacts; the balanced choice below completes it transparently.

| Item | Verdict | Evidence and scope |
|---|---|---|
| A single source expansion must act after cluster composition | correct | Explicit factorials and source-assignment forests |
| Spatial derivatives of a locally integrable base already equal the old rotational insertion | wrong | Missing unary delta contacts, including at separated current positions |
| Balanced mixed matching cancels each ordered odd local contact | correct | Actual first-jet Plemelj calculation and full differential kernel |
| The symmetric correction vanishes | wrong | Explicit Pf and PV terms below |
| First one-barred triple repair is compatible with regrouping | correct at the checked first order | Full-cluster extreme-order average and actual contact kernels |
| Iterated Jordan multiplication supplies that repair | wrong | It has different weights; the correction must be made once on the whole cluster |
| Source-pair Taylor subtraction commutes with external contraction at the checked descendant levels | correct | Independent generating kernel and radial moments |
| K0 insertion forces a Serre obstruction | wrong at the marginal face | Exact zero scalar per K0 placement; remaining degree improves |
| Full mass-two interaction map descends to all canonical relations | unclear | Remaining symmetric forest/quotient identities are stated explicitly below |

No theorem label or physical quantum representation is asserted.

## 1. Domain, expansion, and one application only

Use the inherited physical momenta and covariance, with `d²=−rho`, `a0.d=a2.d=1`, `a0²=a2²=−1`, `a0.a2=1−t`, `rho=2/t`, beta=t−1, gamma=t/(1−t), and s=exp(−i pi rho). Write

`B=:C exp(d.X):`, `C=a0.partial X`, `A=:exp(−d.X):`,

`K_i=:exp[(d+a_i).X+a_i.Xbar]:`,

`F_A=:exp[−k.X−T.Xbar]:`, `U=:exp[(d+a2).(X+Xbar)]:`.

A cluster retains all current/source constituents, their ancestral height order, ordinary Wick normal fields, derivative descendants and independent smooth slot tests. Composition is the already constructed labelled Gaussian Wick product. Denote its result by F. Let W(F;Y_i(w)) mean **one** additional Gaussian contraction with the full source, retaining the source as an integration variable. Its contractions and cocycle are calculated from the vectors; it is not an unspecified interacting word.

Choose one compact infrared function chi for every source. The unextended insertion series is

`E(F)=F − Σ_(i=0,2) lambda_i/pi ∫chi(w) W(F;Y_i(w)) d²w`

`− lambda1/pi ∫chi(w) W(F;V1(w)) d²w`

`+ lambda0 lambda2/pi² ∫chi(v)chi(w) W(F;Y0(v),Y2(w)) d²v d²w`

`+ Σ_(i=0,2) lambda_i²/(2pi²) ∫chi(v)chi(w) W(F;Yi(v),Yi(w)) d²v d²w`,

truncated at mass order two. The input F itself contains the canonical first tails `−t lambda2 K2` and `−gamma lambda1 F_A`; expanding the formula automatically gives all one-tail/one-action and two-tail terms. Do not insert a separate action expansion in each factor before composing F.

The map proposed here is

`I_min(F)=R_min E(F)`.

R_min is the linear extension/forest rule specified in sections 2–5. It acts on actual Gaussian coefficients, not on names for desired relations. At this stage the map is on resolved labelled kernels. Proving it is well defined on a further quotient by all canonical integrated distributional identities is a separate task.

The source-pair counterterm and current improvement are forest representatives. They must occur once. In this implementation the source-free input contains J, the K2 tail and the F_A tail; the source-pair forest produces kappa_h, and the K2/Y0 forest produces the L improvement. An implementation starting with the already improved current must not add that same L contribution a second time.

## 2. The mixed base and the indispensable unary contacts

Put the left current at x, the right current at y, and the Y2 source at w=v+i tau. The basic kernel

`G=1/[(x−v−i tau)(y−v+i tau)]`

is locally integrable in the three relative real coordinates. It has degree −2 there. Set

`N=:exp[d.X(x)+d.Xbar(y)+a2.Phi(w)]:`, `D=d.partial X`, `Dbar=d.partialbar Xbar`.

Away from the pole faces the complete derivative-current kernel is exactly

`[C(x)+beta D(x)−beta partial_x]`

`[Cbar(y)+beta Dbar(y)−beta partial_y](G N)`.

The derivatives act on the entire base and normal field. For example `partial_x N=D N`, so the first bracket produces the actual factor `C(x)+beta/(x−w)`. The script verifies this differential reconstruction without suppressing C,Cbar,D,Dbar. Differentiating the locally integrable base is therefore a concrete minimal extension, with arbitrary smooth tests and all their jets retained.

It is **not yet the rotational prescription**. As distributions on the source plane,

`Pf(1/z²)=−partial_x(1/z)+pi delta²(z)`.

The original smooth Cauchy representative is `P_eta=−partial_z[barz/(|z|²+eta²)]`; it tends to rotational Pf, not to the spatial derivative alone. The inner circular prescription has the same isolated unary extension. Consequently the forest must add the following source/current collapses, with the source action sign and ordered cocycle included:

| Current/source | Unary insertion contribution |
|---|---|
| B/Y0 | `+lambda0 K0` |
| B/Y2 | `+beta lambda2 K2` |
| A/V1 | `−lambda1 F_A` |
| Barred versions | Chirality-exchanged expressions |

For B/Y2, the current-first source cocycle and negative action sign cancel, so beta times the missing pi delta gives `+beta lambda2 K2`. For A/V1 the cocycle is +1 and the action sign stays negative. These contributions exist when another current is separated; they are not merely contacts on the full current diagonal.

This correction invalidates the attractive shortcut that balanced mixed matching alone identifies a rotational source integral with the spatially differentiated base. The unary faces must be retained first.

## 3. Ordered, symmetric and antisymmetric mixed terms

Take the B at r+i sigma delta, the barred B at zero, sigma=+1 or −1, delta>epsilon, and align the reverse word by its free phase s. Denote `partial U` and `partialbar U` by U_L,U_R. The unary replacement sum, in units nu lambda2, has boundary expression

`−beta{beta U[(r+i sigma0)^−2+(r−i sigma0)^−2]`

`− U_L/(r+i sigma0)+(beta U_L+t U_R)/(r−i sigma0)}`.

Its phase-aligned symmetric part is

`−2 beta² U Pf(1/r²) − beta[(t−2)U_L+t U_R] PV(1/r)`.

Its odd local part is

`−i pi sigma t beta partial_x U delta(r)`.

Both the Pf term and the PV term are nonzero. The opposite double-pole delta-prime pieces cancel. The script reconstructs this decomposition from the two Plemelj boundary values, and the root identities `a0+beta d=−(d+a2)` and `beta a2−a0=t(d+a2)` fix its field jets.

Independently integrating the bare action's first normal-field jets gives the same odd contact. One convenient local moment calculation uses

`I21=−pi/z`, `I12=−pi/barz`, `I22=0` away from z=0,

`w I22=pi/barz`, `barw I22=−pi/z`,

where `Imn=∫d²w (z−w)^−m barw^−n` has the inherited Cauchy/rotational treatment and a common infrared cutoff. The complete first-jet kernel contains

`pi beta[ (beta a2−a0).partial X/barz +(a0−beta a2).partialbar Xbar/z ]`.

With the action sign this yields the same odd `−i pi sigma t beta partial_xU delta(r)`. The finite-field/infrared remainder is ordinary source dressing, not another freely chosen diagonal contact. These are local jet statements, not a claim that the displayed rational moments describe every global correlation function.

We now choose the explicit **balanced** finite matching

`C2,sigma=+i pi sigma nu lambda2 t beta partial_x U delta(r)`,

and, in the AA channel,

`C1,sigma=−i pi sigma nu lambda1 gamma partial_x F_A delta(r)`.

The first-order odd bare action contact is thereby canceled in each ordering. No extra phase-aligned even primitive contact is added. This fixes an additional symmetric normalization that the old upper-minus-lower endpoint formula did not fix. It retains the nonzero unary Pf/PV correction above; it does not set the ordinary symmetric dressing to zero.

The actual outer-test first jet cancels in the first mixed face: the action's double-pole pieces are proportional to `beta²(I+Jw)=0`, and the two tail double-pole terms cancel one another. Thus the minimal choice does not require a further `U delta'(r)` term at this order. Arbitrarily adding a phase-aligned even such term would give another map, even though its constant-outer-test commutator would agree.

## 4. One-barred triples: an actual whole-cluster repair

The unary replacement must not be iterated as a Jordan binary multiplication. For a pure left word P and a single bare barred D, let q(P,D) be its free braided phase. For D=barB, `q=s^(#B−#A)`. Let U_i(P) be the sum over replacing one constituent by its entry in the unary table, retaining every field, test and ancestral label. Define its one-barred repair once on the whole cluster by

`J_i(P,D)=1/2[U_i(P) D+P U_i(D)`

`+q(P,D){U_i(D) P+D U_i(P)}]`.

Every product on this line is the canonical labelled Wick product, including its complete boundary distribution. It is the average of the two **extreme barred positions**. It is not the uniform average of all insertion positions and not an iterated pairwise Jordan product.

Here is a concrete B1 B2 barB check. Phase-align all three positions of the bar. Let their raw unary/contact contributions be

`R0`, `R1=R0+Delta1`, `R2=R0+Delta1+Delta2`.

The actual Delta_i is the proper B_i/barB/Y2 contact calculated with the other B as a derivative-current spectator. Its endpoint product is

`Q_g=B(a)U(x)=(a−x)^(g−1)[C(a)+(t−2)/(a−x)] times its normal field`.

The insertion-wise derivative is

`partial_x Q_g−(g−rho)Q_g/(x−a)`.

The extra term follows by differentiating the inherited factor `(a−x)^(g−rho)`. It cannot be dropped before taking the complete boundary family. Multiplying the physical moving-source double/simple pair `beta K,−partial K` by an external regulator R changes its simple coefficient to `−partial(RK)+t R'K`, so the factor t is retained for this repeated-B spectator. This reproduces the actual Delta_i with every first test/normal-field jet.

Balanced proper-face matching contributes respectively

`(Delta1+Delta2)/2`, `(-Delta1+Delta2)/2`, `−(Delta1+Delta2)/2`.

All three repaired expressions are exactly `(R0+R2)/2`. The independent script checks the three values symbolically. This is a concrete regrouping at the first mixed order, using the calculated contact kernels. The triple union has degree g−2, so its remaining generic continuation has no divisor identically g=rho; the already established coherent whole-word boundary construction fixes its overlap. No separate symmetric triple contact is selected.

A uniform average of the three raw positions would instead weight Delta1 and Delta2 by 2/3 and 1/3. It would violate the equal half-contact normalization on proper pair faces. This is why treating the repair as an arbitrary averaging operation would be wrong.

## 5. Source-pair forests through mass order two

R_min uses the actual complete Gaussian kernel, first the fixed unary/mixed extensions above, then the source-pair radial extension with all descendants. The unlike-source regulator is `|r|^(2h−2)`. If a smooth normal/test factor has expansion `Σr^p bar r^q Fpq`, its rotational moments select p=q and give

`pi (R^(2(h+n))−epsilon^(2(h+n)))/(h+n) Fnn`.

The missing-disk representatives are therefore fixed at **every** retained Taylor level, not only at n=0. The leading n=0 term with the double-action coefficient produces

`kappa_h=−lambda0 lambda2 epsilon^(2h)/h`.

For K2/Y0 the basic factor is `r^(h−2) bar r^(h−1)`. Its selected jets have p=q+1 and give denominators h+n. With the action and K2 coefficient, the leading term is the current improvement

`+(t/h) lambda0 lambda2 epsilon^(2h) L`.

For U/Y0, the factor is `|r|^(2h−4)`. The first two selected jets give

`I_0(U)=ordinary source integral`

`−lambda0 epsilon^(2h−2)/(h−1) F`

`−lambda0 epsilon^(2h)/h :C Cbar F: + higher required jets`.

Consequently the balanced lambda2 contact is dressed once by this endpoint map, with the same coefficient t(t−1). This reproduces both endpoint descendants and their h-dependent denominators. The source is not inserted again into each factor of the fused endpoint.

There is an independent contraction/subtraction test. Put the two source positions at r,0 and a separated derivative B at z. Its scalar kernel is

`f(z;r)=[1/(z−r)+(t−1)/z]/[(z−r)z]`

`=z^−3 Σ_(n>=0)(n+t)(r/z)^n`.

The script independently reconstructs n=0,...,4 by differentiating the actual rational kernel. With a second derivative current in the other chirality, the radial forest coefficient is proportional to

`(n+t)² epsilon^(2(h+n))/[(h+n)(zv)^(n+3)]`.

At h=t this reduces to `(n+t)epsilon^(2(t+n))/(zv)^(n+3)`. The first descendant is not obtained by keeping only the vertex factor: fusing the sources first gives the external C contraction `t/z times 1/z+1/z²=(t+1)/z²`. The extra 1 is essential. This independently tests that source-pair subtraction and external Wick contraction agree at the first nontrivial descendant. The general identity is the displayed generating series.

External differentiation commutes with this Taylor extraction while the other variables are separated. For joint limits it must be checked with the proper-face extensions retained. In the bounded one-barred family the new unlike-source union degrees are positive on the open patch 2<t<3; this gives an available domain for the common radial continuation after the old faces. It does not authorize discarding those old faces.

The combinatorial forest also matters. A single source may collapse onto B1 or B2, but not onto both as two independent delta factors. With two different sources and two currents there are exactly two disjoint double-unary matchings. Two faces sharing one source or sharing a current without containment are overlapping faces; they are handled by their union extension, not multiplied. The script enumerates the elementary matchings. Source-pair and mixed/source faces that overlap are subject to the same rule.

A separate two-face calculation fixes a genuinely second-order symmetric contribution. Collapse Y0 onto B(x) and Y2 onto barB(y), with x distinct from y. The unused simple denominators are `(x−y)(y−x)=−(x−y)²`. The two unary pole coefficients are 1 and beta. Hence this disjoint forest contributes `−nu lambda0 lambda2 beta |x−y|^(2h−4)` times its actual `K0(x) barK2(y)` fused normal field. The other source assignment gives the corresponding `K2(x) barK0(y)` term with the same coefficient. The fields are different at separated x,y and must not be combined prematurely. Direct vector multiplication gives the exponent h−2 in each chirality. This is locally integrable on 2<t<3 and is not an optional full-diagonal constant. For two equal Y2 sources the source numerator cancels the same denominator: each assignment gives the regular coefficient −beta²; the factor 1/2! cancels the two assignments. These checks are in `check_second_order.py` and give concrete second-order output of the proposed forest, not just its power count.

The minimal completion adds no new primitive local term on the open positive patch when its full new face is integrable. Continue that complete family, with its fixed proper faces, before physical restriction. This is a definite symmetric normalization choice. It does not mean that every possible finite local modification is ruled out by the original defect data.

## 6. Serre and unit tests

The unary table introduces K0 although the physical conserved B tail contains K2. This could have invalidated a quotient map. Its potentially marginal scalar instead vanishes separately at each placement. In B³A with one B replaced by K0, the other two B fields have scalar polynomial

`[−1/(z1−z2)+1/(z1−a)][+1/(z1−z2)+1/(z2−a)]+1/(z1−z2)²=0`.

The reason is `−a0.(d+a0)=0`. Every term retaining a regular C improves the full collision degree by one. The K0 tail has the same separated current-label phases as K2. Thus the usual shuffle zero does not multiply a new g=rho pole at this face. In A³B with the minority B replaced by K0, the full degree is `3(g−rho)+6`, already positive. The script verifies the scalar identity and degree independently. This checks the K0 contribution to the first Serre variation; it is not a proof of all second-order relation identities for the full map.

The insertion map itself is not unital as an operator-valued map: `I_min(1)` is the interaction insertion series, with its own source-pair forests. Vacuum normalization is a separate scalar operation. Choose a normalized Gaussian functional omega0 with a specified zero-mode/background convention and put

`Z=omega0(I_min(1))`, `omega_lambda(F)=omega0(I_min(F))/Z`.

The formal inverse exists because Z has constant term one. Then omega_lambda(1)=1. This does not imply `I_min(1)=1`, nor turn I_min into an invertible algebra transport.

For example, writing A_i(F)=omega0(I_i(F)) and z_i=A_i(1), the distinct-source normalized coefficient is

`A_02(F)−z0 A2(F)−z2 A0(F)+(2z0z2−z02)omega0(F)`.

Substitution F=1 gives zero at that order; the script checks it. No numerical vacuum expectation is assumed. In particular, a vanishing or nonvanishing vertex one-point function is not chosen silently.

## 7. Exact extent of the construction

The algorithm has explicit data: Gaussian contraction and cocycle, source factorials, a fixed infrared test, the spatial base extension, unary Pf corrections, a whole-cluster balanced mixed normalization, and radial forest coefficients with all retained descendants. It reproduces the established low pair coefficients and the first mixed regrouping, and it counts each source once. It is a usable candidate for further direct computation, rather than a notation declaring interacting words associative.

What has **not** been proved is that R_min annihilates every zero of the integrated canonical Gaussian product at mass order two. The first mixed free-null relation and K0 Serre marginal test pass. A complete proof would also need the full symmetric two-source forest on arbitrary composite/test jets, including its overlap with every canonical relation being used. The positive new-face estimates and the known coherent equal-source boundary family make that comparison feasible; they are not its replacement.

Likewise, the map's bounded one-barred construction is not a closed algebra under arbitrary multiplication. Multiple barred blocks introduce additional grouping data. A normalized correlation functional derived from this candidate would not by itself construct an interacting operator algebra or an invertible Moller/normalization map.

Thus the remaining verification is sharply specified: compare the once-applied forest map on two canonical representatives of the same composed cluster, including all symmetric finite terms, at the mixed-source order. A pair-only correction, a separate interaction series in each binary factor, or an iterated Jordan product cannot stand in for that calculation.

## Evidence and integration

Fresh `check_map.py` and `check_map.out` contain 31 PASS lines and no failures. `check_second_order.py` and its output add six independent passing checks of the disjoint two-source faces and spectator Laurent coefficient. The code imports no previous research script. It verifies the differential base, unary symmetric/odd decomposition, five source-pair Taylor coefficients, first descendant contraction, radial moments, source factorials, compatible unary forests, three mixed positions, K0 Serre scalar and normalized unit. The written derivation states where analytic or quotient verification is still missing.

No book integration occurred. Verified: a concrete once-only insertion series, the missing unary contacts, the nonzero symmetric Pf/PV terms, balanced first mixed regrouping, and the low source-pair/endpoint forest chain. Rejected: a spatial-base-only identification, zero symmetric dressing, double-counted interaction insertions, and iterated Jordan repair. Constructed: a labelled-kernel candidate with an explicit added symmetric normalization. Unresolved: full mass-two quotient descent and all-word identification. All prior files remain unchanged.
