# Adversarial construction audit: the one-barred five-current family

Date: 2026-10-24. Model: Codex (GPT-6). Current reviewed checkout commit: `ae19df62aaa82879da0f4600b74982d700d481b5`. The earlier mathematical reports identify their own earlier commit. Input hashes are in `referee_scratch/sl21_common_product_20261024/audit/inputs.json`.

No earlier research script was imported or rerun as evidence. All writes belong to this fresh report and its assigned scratch directory. No existing file, book file, git index or working-tree state was changed by this worker; no subagent, worktree, commit, push or render was used.

## Outcome

**A genuine common canonical Gaussian product exists for the tested family with at most three B currents, one A and one barred B, through mass order two.** Its one-barred equal-source kernels admit a direct boundary-distribution construction. This includes arbitrary smooth independent slot tests and all regroupings inside that family. The repeated-test Serre relation remains zero when the barred B is adjoined, because the joint boundary limit agrees with taking the quartet limit first. This is an actual partial-product/ideal-stability result, not a conclusion from the five-current chamber residue.

The apparent two-tail directed cycle is absent in the complete kernel: its numerator and derivative contraction reduce it to a sum of single antichiral leaves. The same reduction applies to the one- and two-equal-source action kernels. It controls the bulk regulator remainders that the local contact computation alone did not control.

The first-order mixed contour defect agrees with the inherited-height canonical contact after the prescribed finite matching. The action contribution cancels that finite matching locally, insertion by insertion, and the remaining regulator multiples belong to the regular complete-word family. At unlike-source order the new faces have no divisor containing the whole physical surface; an explicit positive physical patch establishes that there is no new analytic obstruction in this bounded family.

**This does not identify every complete source-dressed word with source-free canonical multiplication.** Ordinary interaction dressing remains. A single interacting insertion/state map, including symmetric contact values and its action on all composites, is still not constructed. Nor does a one-barred family close under arbitrary multiplication. These are the precise remaining limits, not a surviving two-tail pinch or a demonstrated five-current anomaly.

| Claim | Verdict | Scope |
|---|---|---|
| Two K2 tails and a bare barred B force a cyclic singular graph | Wrong | The complete antichiral numerator cancels the putative cycle |
| Canonical one-barred Gaussian words have compatible limits | Correct | Full derivative polynomials, all tail placements through mass two, independent compact slot tests |
| One- and two-equal-source action kernels have a same-diagonal g pole | Wrong in the coherent whole-word construction | Exact leaf decomposition and uniform Fourier pushforward |
| The one-chamber residue is a fullword anomaly | Wrong if inferred | Its chamber restriction is not the complete smooth-test distribution |
| First-order mixed action plus finite matching equals zero on the relevant contact faces | Correct | Insertion-wise contact convention, all spectator and moving-test jets retained |
| The first B³A/barred-B mixed–Serre overlap is stable | Correct in the canonical common family with the matched contour defect | Grouped versus simultaneous limits agree |
| New unlike-source faces obstruct physical restriction | Wrong for this family | Fifteen current-containing subsets plus source-only face; all positive on 2<t<3 |
| All source-dressed words are already values of one specified interacting map | Unclear | Contact/ideal checks do not specify the complete map |
| A closed full charge algebra or physical quantum representation is established | Unclear | More barred words, mixed couplings at higher order, and global operator domains remain outside this result |

## 1. The actual two-tail kernel

Use the primary vectors, with v=d+a2 and beta=t−1. Direct multiplication gives

`v²=1−rho`, `a2²=−1`, `a2.v=0`, `−a0.a2=beta`.

Consequently two K2 tails have holomorphic exponent g−1 and antiholomorphic exponent +1. Put x=bar z_i, y=bar z_j, u=bar u for the barred B coordinate. Its complete antichiral scalar and derivative factor is

`F [Cbar(u)+beta/(u−x)+beta/(u−y)]`,

where

`F=(x−y)/[(x−u)(y−u)]=1/(y−u)−1/(x−u)`.

Since `partial_u F=F[1/(x−u)+1/(y−u)]`, this factor is

`Cbar(u) F−beta partial_u F`

`=Cbar(u)[1/(y−u)−1/(x−u)]−beta[1/(y−u)²−1/(x−u)²]`.

Every term has just one antichiral edge incident to the barred B. The regular normal field Cbar is retained. This is an identity of the complete normal-product coefficient, not merely its fully contracted scalar part. It holds at finite complex coordinates and therefore with the inherited finite heights. The single common cocycle of each word is unaffected by this algebraic decomposition.

The remaining B and A have zero right momentum. Thus the rest of the singular graph is wholly chiral; barred B is a leaf. There is no positive circulation. If the barred slot itself is its tail, a2.v=0 removes every right cross-edge to it; the mutual right factor of the chiral tails is polynomial. That case is purely chiral after multiplication by a polynomial.

## 2. One formula covers action sources as well as tails

For m equal Y2 constituents, regardless of whether a constituent is an action source or is fused into a K2 tail, put

`Delta(x)=product_(i<j)(xi−xj)`, `F_m=Delta(x)/product_i(xi−u)`.

Partial fractions give the exact identity

`F_m=sum_i (-1)^(m−i−1) Delta(x with xi omitted)/(xi−u)`,

with indices i=0,...,m−1. Its coefficients are polynomials, not new singular edges. For m=2 they are −1,+1. For m=3 they are x1−x2, −x0+x2, x0−x1.

The full barred derivative factor is

`Cbar F_m−beta partial_u F_m`.

It is therefore again a sum of simple and double single-leaf edges with polynomial coefficients. This independently proves the reduction for the two-action, one-action/one-chiral-tail, and two-chiral-tail channels at lambda2 squared. If barred B is its tail, all remaining root/right-tail contractions vanish by a2.v=0. At first order the statement is simply the one-source leaf.

The same equal-node argument applies to two Y0 action sources, with the derivative coefficient 1 instead of beta. It uses a0²=−1 and a0.d=1. No unlike-source statement is being hidden inside the Vandermonde formula: that sector is treated separately in section 7.

Normal products are common to the decomposed scalar terms. Uncontracted holomorphic C fields, their pair contractions, and every regular field jet remain present. Their contractions only shift holomorphic edge exponents by integers. The formula consequently controls the full operator-valued Taylor coefficients at the finite perturbative order, not only a selected coherent background.

## 3. Uniform boundary values and actual grouping

For each term after this reduction, orient all chiral edges by their inherited imaginary-height order. A source w is an additional chiral vertex at its own imaginary coordinate. The sole antichiral edge, when present, attaches barred B as a leaf; either orientation of a leaf keeps the graph acyclic.

For an edge exponent p, the ordinary boundary factor is, up to an entire phase, the Fourier transform of

`k_+^(p−1)/Gamma(p)`.

This is an entire distribution family in p, as follows by repeated differentiation of the one-dimensional primitive. Choose an integer topological potential c_i for the acyclic graph, decreasing by at least one along each edge. For positive edge momenta and their incidence image q,

`sum_i c_i q_i = sum_(i→j)(c_i−c_j) k_ij >= sum_edges k_ij`.

The potential is auxiliary; it is not the vanishing physical height gap. The inequality gives properness of the incidence pushforward and a polynomial bound uniform as any or all physical height gaps vanish. Actual nonnegative gaps enter only as damping factors `exp(−sum epsilon_ij k_ij)`. Schwartz tests in the spatial vertex variables suppress every edge momentum through the properness bound. Integration by parts in the entire one-edge distributions gives the same bounds on compact parameter sets, including parameter derivatives. Thus the boundary family is entire in g and continuous under partial as well as simultaneous height collapse.

There are only finitely many orientations as the imaginary coordinate of one or two sources crosses the current heights. On each region the same proof uses a discrete potential. It gives a uniform bound independent of the widths of these regions. Compact source integration against smooth normal-field/test coefficients preserves the estimate. At current/source crossing lines use the prescribed first Cauchy distributions; differentiation produces their usual local delta jets, discussed below.

This supplies the step that pointwise scalar power counting cannot supply. Binary Wick grouping is already associative at distinct finite heights: each contraction matching occurs once and the momentum/cocycle/constituent-label factors are additive. The uniform limits preserve that equality. Smearing independent current coordinates and grouping their integrations are compatible pushforwards of the same distributions. Every subword of the tested one-barred family belongs to the same class. Hence its regroupings agree.

This statement allows independent slot tests for the product. A Serre relation repeats one B(f_B), so its three B tests must agree, or be correctly polarized. The unpolarized relation with three arbitrary unrelated B tests is not asserted.

## 4. The first five-current Serre overlap

A whole-word entire family alone would not justify multiplying a separated chamber shuffle zero by a potentially singular chamber restriction. Instead use the actual grouped limit.

In the source-free canonical quartet, the bare scalar Wick polynomial and the coefficient of one constant regular C vanish. With `f_i=−sum_(j!=i)1/(zi−zj)+1/(zi−a)`, the identities are

`f1 f2 f3+sum_(i<j) f_k/(zi−zj)²=0`,

`sum_(i<j)[f_i f_j+1/(zi−zj)²]=0`.

For a K2 in slot i, the complete remaining two-C scalar coefficient is individually

`(t−2)/rho partial_a(R S_i)`,

where `R=product_(i<j)(zi−zj)^g product_i(zi−a)^−rho` and

`S_i=(a−zi)²/product_(j!=i)(zi−zj)²`.

I independently reconstructed these identities from the primary contractions. In a middle chamber the primitive vanishes at both neighboring B endpoints in the initial negative-rho domain. Moving its derivative onto an arbitrary smooth test/normal coefficient improves the marginal radial degree by one. The two-K2 and minority F_A terms have positive full degree; the regulated L term vanishes in its initial h-domain. These facts justify the canonical quartet zero with repeated-species tests; they agree with the independently produced canonical report without relying on its check count.

Now keep barred B at a separated imaginary height and first form that canonical quartet. Its Serre combination is zero. Section 3 permits the subsequent barred-height collapse and identifies the result with the simultaneous five-current boundary value. Therefore

`(B³A−(1+s+s^−1)B²AB+(1+s+s^−1)BAB²−AB³) barredB = 0`

through the stated canonical mass order, with the indicated species tests. The same construction allows the barred slot on the other side and the intermediate mixed commutations with their actual cocycles. It provides a tested partial associative product and ideal-stability statement. The residue 0.093460526 from one decorated chamber is not this limit; its chamber partition must cancel in the whole distribution.

## 5. General spectator pole and the mixed contact bridge

There is a short exact way to generalize the repeated-B Laurent check. Let K be the unregulated source-collapse normal product with every separated spectator retained. Its physical operator OPE has

`D=beta K`, `S=−partial_w K`.

This follows from `a0+beta d=−(d+a2)` before taking any spectator expectation. Let R(z) be the product of the moving B's extra like-label factors. Multiplication by R(z) gives

`D_g=beta R(w)K`,
`S_g=−R(w)partial_w K+beta R'(w)K`

`=−partial_w K_g+t [R'(w)/R(w)] K_g`.

For any number of B spectators, `R'/R=(g−rho)sum_j 1/(w−a_j)`. Derivative spectators are covered because the initial OPE is an operator identity; no regular C or derivative-pair contraction was dropped. An extra moving test f multiplies the double-pole residue by f and adds beta f'K to the simple coefficient. Both must be kept.

With Q the fully regulated spectator–U product and

`Rreg=(g−rho)sum_j Q/(x−a_j)`,

the two strip-edge moments and inner circle give the following coefficients in the common normalized contact:

`action = −t(t−1)[partial_x Q−Rreg]`,

`both canonical tails = t²[partial_x Q−Rreg]`,

`finite matching = t(t−1)[partial_x Q−Rreg]`.

The derivative is the full spatial derivative in the first term, while subtraction of Rreg makes this the insertion-wise derivative of U. The action plus matching is exactly zero off the physical surface. This cancellation is stronger than merely discarding a g−rho term. The first moving-test jet from the action is proportional to beta²(I+Jw)=0; the two tail double-pole test terms cancel one another. These statements hold with the full spectator jets, since the circle/edge calculation is local and linear in those jets.

Keep moving masks inside these differentiations. Product rules generate derivatives of R and derivatives of masks, rather than permission to set masks to one. At finite cutoff, every extra derivative of a regulator carries g−rho, or h−t in the unlike-source sector. Derivatives of circular masks are supported on the relevant collision faces. Their isolated source/current Cauchy terms are the fixed lower prescription; simultaneous new boundaries vanish in an initial domain such as Re g>6, Re rho<−1, Re h>4. The complete kernels and all collapsed rational descendants of section 3 have regular physical restrictions, so these regulator multiples cannot acquire a hidden 1/(g−rho) remnant in the one-/equal-source sector.

Changing the source-root part of a fused vertex to the coherent ancestral height can change integer pole i0 assignments. In the initial domain their delta jets are multiplied by sufficiently vanishing like-current g powers or opposite-current −rho powers. They vanish there, including simultaneous faces. Continuation of that equality fixes the inherited-height contact; its branches are not reassigned independently in each spatial chamber.

Thus the tested mixed contour defect matches the canonical t² contact, and the grouped five-current Serre overlap is compatible with it. This is a comparison of contour defects and their propagated endpoints. It is **not** an equality of an entire interacting word with its source-free value.

## 6. Why neither a barred-tail logarithm nor a Cauchy contact defeats the proof

A rational primitive in the barred-tail coordinate can contain logarithms. That is an obstruction to the proposed rational-primitive shortcut, not to the actual kernel: the barred-tail kernel itself is a chiral boundary product with integer-shifted edge exponents and is covered directly by section 3. No logarithmic primitive is used to discard its proper faces.

There is also a genuine local distinction between circular and iterated Cauchy extensions of a double pole. With ordinary area normalization, differentiating `1/(x+iy)` gives

`partial_x[1/(x+iy)]=−Pf[1/(x+iy)²]+pi delta(x)delta(y)`.

Thus replacing a circular finite part by horizontal Cauchy differentiation without its contact would be wrong. In the present construction those contacts collapse a source onto a current, producing the K2/barK2 or mixed endpoint terms already included. Their descendants remain in the same chiral/leaf class. When two equal sources are present, collapsing one onto barred B removes its anti interaction with the other because a2.(d+a2)=0; collapsing onto a left B gives the Vandermonde-tail case of section 2. This makes the same-source forest comparison substantive rather than a claim that all regulators are interchangeable.

## 7. Unlike-source order: all new faces are transverse

The Vandermonde proof does not apply to a Y0/Y2 pair. Its factor is `|w0−w2|^(2h−2)`. Independently count a current-containing cluster with a=0 or 1 A, b=0,...,3 B, d=0 or 1 barred B, and both sources. Set

`E=g binomial(b,2)−rho a b`.

The source/current contractions contribute 2a−2(b+d), the two source factors contribute 2h−2, and the potentially contracted derivatives cost b+d. The relative real measure has dimension a+b+d+3. Hence

`degree = E+3a−2(b+d)+2h+1`.

There are fifteen nonempty current subsets of this form. The source-only face has degree 2h. The full five-current face is

`3(g−rho)+2h−4`.

On g=rho=2/t,h=t, every one of the fifteen degrees is strictly positive for 2<t<3. A convenient exact proof is to multiply by t and set t=2+q. Every resulting polynomial has nonnegative coefficients and a positive q coefficient. The tight faces include

`6/t+2t−7` for (a,b,d)=(0,3,1),
`2/t+2t−5` for (0,2,1),
`2t−4` for (1,3,1).

The script reconstructs this list from the actual three-vector contractions; it does not take the counting formula as input.

This also explains the nested issue. Old one-source faces are entire g-families on arbitrary smooth tests, by section 3. With the other source separated, its Wick factor is just another such smooth test. Resolving their union adds a radial scale containing both sources, and every new Mellin denominator has a nonzero h coefficient. Old angular faces cannot resurrect a g−rho pole: their whole test-distribution family is entire, and all of their boundary sectors were retained before continuing the new radial scale. The positive physical patch verifies uniform control there after the fixed old faces. Continuing it can create transverse resonances, not a divisor containing the entire physical surface.

At the physical restriction every such denominator has the form `2t+A/t+B` (with integer Taylor shifts included). It cannot vanish identically in t. Generic nonresonant t is still required after continuation. Counterterm denominators h,h−1 have the same property. The h-dependent P/current/endpoint coefficients must be the already derived 1/h,t/h,1/(h−1) coefficients, and all constituent labels must remain attached to their fusions.

Therefore the bounded five-current family has no new same-physical-surface pole through mass order two. The even-root order has one source and at most one nonpolynomial anti leaf, so the same acyclic argument covers its analytic boundary existence. This is a regularity/descent result; it is not a fresh computation of every unlike-source finite contact or of a second-order interacting contact cochain.

## 8. Precisely what remains to be constructed

There is now an explicit canonical finite-family composition, an actual first mixed–Serre overlap calculation, and control of the candidate same-surface pole through mass two. Repeating the chamber-residue diagnostic or counting an unreduced two-tail cycle would not advance the remaining question.

The outstanding map is concrete: one must specify a single interaction/insertion map whose value on a composed canonical Wick cluster reproduces the complete prescribed Euclidean source-dressed word, including its ordinary action dressing, symmetric density contacts, and all composite/test jets. The map must be applied once after composition and respect the Gaussian quotient and the relations being used. Independently dressing binary blocks double-counts sources. A near-identity composite renormalization fixing the single currents cannot create the canonical mixed coefficient from the raw t coefficient, because its change of the free braided commutator is zero.

Our local cancellation proves the required antisymmetric mixed contact comparison at first order. Our uniform boundary proof removes the suspected five-current evanescence obstruction, including the equal-source forest at second order. Our unlike-source count removes its new analytic obstruction. These do not determine all symmetric contact values of the interaction map or prove its factorization on a closed space with arbitrarily many barred generators. This is the exact remaining equality problem. It is narrower than the previous absence of any common five-current product, but it still matters for claiming the requested full interacting associative algebra.

## Evidence and preservation

Fresh evidence is `referee_scratch/sl21_common_product_20261024/audit/audit.py`, `audit.out`, and `unlike_faces.json`. The final run has 65 PASS lines and no FAIL lines. Actual output includes:

```text
PASS two-tail holomorphic exponent becomes g-1
PASS full derivative-current leaf decomposition m=2
INFO leaf coefficients 3 [x1 - x2, -x0 + x2, x0 - x1]
PASS arbitrary-spectator simple residue transformation
PASS action plus insertion-wise matching cancels off surface
PASS all current/source height chambers give DAG m=2
INFO directed configurations checked 2 10080
INFO full unlike five-current degree 3*g + 2*h - 4 - 6/t
PASS individual one-tail minority derivative 0
PASS individual one-tail minority derivative 1
PASS individual one-tail minority derivative 2
```

The first run's failed symbolic check is retained in `audit_initial.out`: factoring `(b−1)(b+1)` before substituting b²=t concealed the square from substitution. Expanding before substitution fixed that checking-code error. No mathematical coefficient or pre-existing input was changed. The counts support the written derivations; they are not a substitute for the uniform-limit or partial-product argument.

## Integration

None. Independently verified: exact Vandermonde/derivative cancellation of the proposed cycle, actual tail and source edge powers, arbitrary-spectator Laurent transformation, action/matching cancellation, source-free minority primitives, and every new unlike-source face. Constructed: a uniform coherent boundary product for the bounded one-barred family and its regroupings. Rejected: the cycle-only obstruction and the inference of a fullword anomaly from one chamber. Established at the stated formal scope: the first canonical five-current mixed–Serre overlap and absence of the candidate same-surface pole through mass two. Still unresolved: equality of all complete interacting words under one specified interaction map, a closed product for arbitrary barred words, and a quantum Hilbert-space symmetry algebra.
