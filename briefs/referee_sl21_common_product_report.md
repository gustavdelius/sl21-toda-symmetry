# Common associative product for the sl(2|1) Toda candidate charges

Date: 2026-10-24. Model: Codex (GPT-6). Initial reviewed commit: `745586e6a6f7bca5ed91a1e9cf89d6ef5df387f0`, branch `main`. The shared checkout advanced externally to `ae19df62aaa82879da0f4600b74982d700d481b5` during the work; the final audit records that later commit. All 275 snapshotted sl21 input files retain their original hashes.

This continues the first gap in the completed coordinating investigation. The editor explicitly authorized independent subagents and retained the referee-only restrictions. There were at most two active workers. This team treated all earlier reports, scripts, book files and indices as read-only; concurrent repository changes are recorded separately below. The [input snapshot](../referee_scratch/sl21_common_product_20261024/preexisting_hashes.json) records SHA-256 hashes of 1,882 pre-existing files, including the previous investigation's outputs. The commit alone does not identify that untracked evidence.

## Outcome

**A common formal canonical product has been constructed in a substantial finite-order sector; the full interacting-product gap is not yet closed.** The construction is an actual Gaussian Wick multiplication with compatible boundary values, rather than an abstract algebra obtained by imposing the desired relations. Its purely chiral sector is associative through mass order two and satisfies both Serre quartics. The construction extends to the tested one-barred-current family: an exact partial-fraction identity removes a seemingly dangerous two-tail pinch. The lowest repeated-current mixed contact agrees with the canonical coefficient. The independent final audit establishes compatibility of the first five-current mixed–Serre overlap with the matched contour defect, within this bounded formal family.

The source-dressed product needs a separate comparison. Its lower pair coefficients already differ from the source-free canonical ones. A finite normalization transport that fixes the individual currents cannot repair the leading mixed normalization by itself. The remaining question is therefore more specific than before: construct one factorization-compatible interaction/contact map on this Gaussian product, or exhibit a discrepancy that it cannot absorb, and extend the tested one-bar family to products with multiple barred currents.

The independent final audit and the exact extent of its source-forest comparison are recorded in section 7. No physical Hilbert-space representation, all-order conservation, global endpoint limit, or theorem-label upgrade is asserted.

## 1. Definitions and acceptance criteria

Use the bilinear Gaussian covariance X(z)X(w) approximately −log(z−w), with

\[
a_0=(ib/\sqrt2,b/\sqrt2,i),\quad a_2=(ib/\sqrt2,-b/\sqrt2,-i),\quad
t=b^2,\quad T=a_0+a_2,\quad d=-T/t,\quad k=d+T.
\]

Thus a0²=a2²=−1, a0.a2=1−t, ai.d=1 and d²=−rho, where rho=2/t. Write beta=t−1, eta=t−2, gamma=t/(1−t), nu=exp(−i pi/t), and s=nu²=exp(−i pi rho). All four currents are even. The left densities are

\[
J_A=:e^{-dX}:,\qquad J_B=:Ce^{dX}:,\qquad C=a_0\partial X,
\]
\[
H_B=-t\lambda_2K_2+(t/h)\lambda_0\lambda_2\epsilon^{2h}L,
\qquad H_A=-\gamma\lambda_1F_A,
\]

where K2=exp[(d+a2).X+a2.Xbar], L=C exp[k.X+T.Xbar], and F_A=exp[−k.X−T.Xbar]. Barred expressions exchange chiralities. The charges integrate J+H with factor (2pi i)^−1. Lambda0 and lambda2 have mass order one; the independent lambda1 has mass order two.

The source-dressed convention is S_int=pi^−1 integral(lambda0 Y0+lambda2 Y2+lambda1 V1+kappa_h P), with Yi=exp[ai.(X+Xbar)], V1=exp[−T.(X+Xbar)], P=exp[T.(X+Xbar)] and kappa_h=−lambda0 lambda2 epsilon^(2h)/h. There is no independent finite P coupling. The auxiliary exponent h replaces the unlike-source power by |w0−w2|^(2h−2). A second auxiliary exponent g replaces like-current powers rho by g, preserving original current labels even in a tail or composite. Physical restriction is g=rho and h=t, at generic nonresonant coupling.

For labels n=(n_L,n_R,n_0,n_2), retain

`epsilon(n,m)=(-1)^[(n_L+n_R)(m_0+m_2)] exp[-i pi(n_L m_R−n_R m_L)/t]`.

The canonical BB endpoint is U=exp[(d+a2).(X+Xbar)]. Its leading coefficient is t²lambda2; the AA endpoint coefficient is gamma²lambda1. The radial bare-current/action computation instead gives t and gamma²/t before the prescribed finite matching. A common product must preserve the intended normalization and must accept arbitrary smooth slot tests. A Serre relation for repeated charges may require repeated tests even when the multiplication itself accepts independent tests.

A genuine associativity test must compare binary compositions with their full descendant and test-function jets. Equality of separated phases, or continuing an independently named n-ary expression, is insufficient. Likewise, an isolated chamber residue is insufficient to refute a whole ordered distribution.

## 2. A genuine associative starting multiplication

The [construction report](referee_sl21_common_product_construction_report.md) builds a formal Gaussian extension. Attach four additive current-label counters and two source counters to each vertex, separately from its physical momentum. Auxiliary Gaussian contractions supply one (g−rho) logarithm per like-current pair and one (h−t) pair of holomorphic/antiholomorphic logarithms per unlike-source pair. A fused block retains the full multiset of counters. The physical derivative C differentiates only physical fields.

For decorated exponentials the binary multiplication is the usual Wick rule with the physical cocycle, the physical contractions, and these auxiliary contractions. If b(p,q) is the extra bilinear contraction, then

`b(p,q)+b(p+q,r)=b(q,r)+b(p,q+r)`.

The cocycle obeys the same two-cocycle identity. For derivative fields, each final Wick matching is obtained once under either binary grouping. Equivalently the cross-contraction operators commute, and both groupings give the same exponential of the three pair-contraction operators. This is an explicit multiplication on decorated Wick functionals, not a quotient by a desired Serre relation.

The coordinator independently checks the additive-label identity, the interaction-counting issue and the normalization obstruction in [associativity_checks.py](../referee_scratch/sl21_common_product_20261024/coordinator/associativity_checks.py), with [output](../referee_scratch/sl21_common_product_20261024/coordinator/associativity_checks.out).

### A nontrivial regrouping with descendants

For B(z1), B(z2) and U(x), put r12=z1−z2 and ri=zi−x. The scalar vertex factor is

`r12^g r1^(g−1) r2^(g−1)`.

Both (B1 B2)U and B1(B2 U), evaluated by actual cross-contractions, give the full polynomial

`(C1−1/r12+eta/r1)(C2+1/r12+eta/r2)+1/r12²`.

At the BB midpoint the leading composite is `(C²+C')exp(2d.X)`. Contracting it with U gives

`R^(2g−2)[(C+eta/R)²+C'−eta/R²]`, with R the midpoint-to-endpoint separation.

The C' contraction is essential. The independent worker also verifies the second pair jet and the arbitrary test jets (f'k−fk')/2 and (f''k+fk''−2f'k')/8. Dropping the constituent labels would instead use rho−1 against U and fail off the physical surface.

An insertion-wise derivative likewise differs from differentiating the entire auxiliary product:

`W_g(B,partial_x U)=partial_x W_g(B,U)+(g−rho)W_g(B,U)/(z−x)`.

This difference cannot be discarded before checking the larger collision on which it is inserted.

## 3. An actual repeated-current mixed overlap

The [independent overlap report](referee_sl21_common_product_overlap_report.md) calculates BBbarB with one Y2, keeping the other B as a derivative-current spectator. Let that spectator be at a, the moving B at z, and the source at w. Its exact holomorphic polynomial, before the source collapse, is

\[
\frac{(z-a)^g}{(z-w)(a-w)}
\left\{\left[C(z)-\frac1{z-a}+\frac\beta{z-w}\right]
\left[C(a)-\frac1{a-z}+\frac\beta{a-w}\right]
+\frac1{(z-a)^2}\right\}.
\]

If K is the actual regulated B(a)K2(w) product and delta=g−rho, the double and simple moving-current coefficients are

`D2=beta K`, `S2=−partial_w K+t delta K/(w−a)`.

The factor t in the second formula is new. It comes from differentiating the regulator in both the double-pole coefficient and the fused endpoint derivative. The earlier h-regulated separated-Yukawa result does not supply it.

Let Q_g=B(a)U(x), including its inherited regulator. In units 2 i nu lambda2, the actual action contribution is −t(t−1)[partial_x Q_g−delta Q_g/(x−a)]. The two tails sum to t²[partial_x Q_g−delta Q_g/(x−a)]. The prescribed finite matching adds t(t−1) times the same insertion-wise derivative. Consequently the canonically matched contact is

\[
t^2\left[\partial_xQ_g-\frac{g-\rho}{x-a}Q_g\right].
\]

This includes the spectator's C contractions and all first normal-field jets. The first moving-test derivative cancels between the action's strip moments and the two tail double poles. It is not inferred from a constant-test identity.

The complete BBbarB/Y2 union has degree g−2. Its Taylor divisors g−2+j are not identically the physical diagonal. Thus at generic coupling this regulator violation has zero diagonal restriction after the stated proper-face extensions. This verifies the first repeated-current mixed overlap; it does not by itself prove the whole barred algebra associative.

The same Laurent transformation has an elementary general form. Multiplying a physical pole `beta K/q²−partial K/q` by R(z) changes its simple coefficient to `−partial(RK)+(beta+1)R'K`. For R=product_j(z−a_j)^(g−rho), this gives the t delta sum over all repeated B spectators. The final audit tests how far this identity extends through the full source forest.

## 4. Both Serre quartics in the canonical product

The new [canonical-product report](referee_sl21_common_product_canonical_report.md) independently evaluates the source-free Gaussian product of J+H, without putting an action insertion inside each multiplication. For A=A(f_A), B=B(f_B), with the same test on each repeated species, it establishes through mass order two

\[
B^3A-(1+s+s^{-1})B^2AB+(1+s+s^{-1})BAB^2-AB^3=O(m^3),
\]
\[
A^3B-(1+s+s^{-1})A^2BA+(1+s+s^{-1})ABA^2-BA^3=O(m^3).
\]

The check includes every tail placement, the full derivative Wick polynomials, actual cocycles, all proper subsets and the independent lambda1 channel. The two marginal scalar kernels are derivatives **separately for each tail placement**. No cancellation against an action contribution is needed to remove those marginal poles.

For B positions zi and minority A at a, write R=product_(i<j)(zi−zj)^g product_i(zi−a)^−rho and Si=(a−zi)²/product_(j!=i)(zi−zj)². The scalar K2 tail in slot i is proportional to

`(t lambda2 eta/rho) partial_a(R Si)`.

For three A positions ai and minority B at z, its F_A tail analogue is

`(gamma lambda1/rho) partial_z(R_A S_i)`.

In each middle spatial chamber the primitive vanishes at its distinct internal endpoints in the initial rho<0 domain. Integration by parts moves one derivative onto the arbitrary smooth normal/test coefficient. This raises degree 3(g−rho) to 1+3(g−rho), removing the candidate diagonal pole. Proper-face extensions are transported from that same domain. The two-K2 full degree is 3(g−rho)+2; its right-moving factor is polynomial. The L improvement vanishes in the initial h>0 source-free limit and in the chosen positive physical patch.

The coordinator independently reconstructs each of the three K2 and three F_A scalar primitives, and the lower pair coefficients, in [canonical_primitive_checks.py](../referee_scratch/sl21_common_product_20261024/coordinator/canonical_primitive_checks.py), with [output](../referee_scratch/sl21_common_product_20261024/coordinator/canonical_primitive_checks.out).

The product allows arbitrary independent slot tests. The displayed unpolarized polynomial does not: with unrelated B tests selecting one labelled middle chamber its coefficient is s^−1−s, generically nonzero. This is a distinction between the domain of the multiplication and the repeated-charge relation, not a failure of the latter.

## 5. Why complete boundary values can be composed

For a totally ordered set of imaginary heights orient each singular chiral factor from earlier to later height. Up to an entire scalar phase, `(xi−xj+i0)^−p` is the Fourier transform of `k_+^(p−1)/Gamma(p)`. Integration by parts continues the latter as an entire distribution in p. Let q be the incidence image of positive edge momenta k. Choose an integer topological potential c on the directed acyclic graph. Then

\[
\sum_i c_iq_i=\sum_{i\to j}(c_i-c_j)k_{ij}
\ge\sum_{i\to j}k_{ij}.
\]

The incidence map is proper on the positive cone. Pushing forward the tensor product of the normalized one-edge distributions gives an entire family of complete boundary products. Derivative contractions shift powers by integers; polynomial factors and smooth normal-field jets preserve the argument. Actual imaginary gaps supply damping bounded by one. The topological bound does not shrink when some physical gaps shrink. Thus the joint boundary value is compatible with partial gap collapse and with grouping integrations against compact smooth tests.

This is the analytical ingredient that lets separated Wick associativity pass to the purely chiral formal product. It also explains why chamber finite parts cannot be assigned independently. The coordinator independently verifies the cone identity and its opposite-orientation failure in [chamber_and_dag_checks.py](../referee_scratch/sl21_common_product_20261024/coordinator/chamber_and_dag_checks.py), with [output](../referee_scratch/sl21_common_product_20261024/coordinator/chamber_and_dag_checks.out).

### A real chamber residue, followed by a real cancellation mechanism

For B³AbarB/Y2 the complete relative degree is 3(g−rho)−1. Its first jet could give a pole on g=rho itself. After a mixed fusion at x, the scalar B B A U kernel is `(eta/rho) partial_a(R S_x)`. Changing insertion-wise to whole-endpoint differentiation adds delta times

`(eta/rho)[1/(x−z1)+1/(x−z2)] partial_a(R S_x)`.

In the chamber x=0<a<z1=uR<z2=R, a=uyR, the Jacobian is uR². Moving the a derivative onto a smooth first test jet leaves delta integral R^(3delta−1)dR times an angular coefficient. Its positive magnitude on the diagonal is

\[
I(\rho)=\int_0^1du\int_0^1dy\,(1+u)u^{-\rho}(1-u)^\rho
y^{2-\rho}(1-y)^{-\rho}(1-uy)^{-\rho}.
\]

For −1<rho<0 this is finite and positive. At rho=−sqrt(2)/5 the coordinator's independent quadrature gives I/3=0.093460526056963, agreeing with the worker's result. A bare factor g−rho is therefore insufficient to discard an individual chamber.

Nevertheless this is **not** a full-word anomaly. With U's left coordinate and all its inherited regulator factors placed on one consistent ancestral height sheet, the collapsed contact is a chiral acyclic boundary product. Its full family is entire in g by the preceding proof. Hence the complete chamber residues cancel for arbitrary smooth tests. A chamber indicator is not such a smooth test through the full diagonal. The report retains both the nonzero chamber calculation and the cancellation proof; neither is silently substituted for the other.

### One barred current and two equal-root constituents

The candidate cycle in the one-bar family also disappears by an actual Wick identity. With x,y the two right-moving Y2 constituent coordinates and u that of barB,

\[
F=\frac{x-y}{(x-u)(y-u)}=\frac1{y-u}-\frac1{x-u}.
\]

The complete right-moving derivative polynomial is `Cbar F−beta partial_u F`. It is therefore a sum of simple and double **single-edge** distributions. Each term attaches barB as a leaf to the remaining holomorphic graph. If barB itself is its K2 tail, `(d+a2).a2=0` removes the cross-oscillator edge altogether. The Vandermonde generalization gives the same leaf reduction for further equal-root constituents.

The coordinator independently verifies this identity, including its derivative term for one, two and three constituents, in [tail_and_unlike_checks.py](../referee_scratch/sl21_common_product_20261024/coordinator/tail_and_unlike_checks.py), with [output](../referee_scratch/sl21_common_product_20261024/coordinator/tail_and_unlike_checks.out). The [primitive checks](../referee_scratch/sl21_common_product_20261024/coordinator/primitive_checks.py) and their [output](../referee_scratch/sl21_common_product_20261024/coordinator/primitive_checks_initial.out) separately reconstruct the five-current scalar derivative and the collision degrees from Gaussian contractions.

In this limited canonical family, the uniform boundary construction permits taking the four-current block limit first with barB separated, and then taking its remaining gap to zero. The same joint distribution results. Consequently a canonical Serre zero remains zero under these tested one-bar compositions. This is a substantive regrouping check. It is not a construction of every word containing arbitrary numbers of both barred generators.

## 6. Source dressing and the obstruction to an easy identification

The canonical lower pairs differ from the action-dressed lower pairs. Suppressing their common phase, fused field and r^(g−2) factor,

| Pair and coupling | Canonical tails | Action contribution | Combined coefficient |
|---|---:|---:|---:|
| BB, lambda2 | 2t(t−2) | −2(t−1)(t−2) | 2(t−2) |
| AA, lambda1 | −2gamma | −2 | −2/(1−t) |

This is an explicit comparison, not a hypothetical scheme warning. Putting the interaction into a single linear state/action map could produce these additional terms without changing multiplication. Putting a separate interaction expansion into each binary block double-counts the same source. The necessary single map, including its contact rules and dependence on observable-associated source holes, must be specified and shown compatible with regrouping.

There is also a verified obstruction to one tempting shortcut. Suppose a raw associative product mu=mu0+lambda2 mu1 is transported by Z=Id+lambda2 Z1, fixing B and D=barB. Then

`mu'=Z^(-1) mu(Z tensor Z)`

changes the leading braided defect by

`−Z1[mu0(B,D)−s mu0(D,B)]=0`,

because the free left/right factors satisfy that exact braided relation. It cannot create the required change from t to t². This excludes that fixed-generator coboundary construction. It does not exclude a genuine contact deformation, a different canonical realization, or changed single-current normalization. A common rescaling of all four charges by sqrt(t) converts both raw diagonal amplitudes to the canonical ones and preserves quartic zeros, but it changes the specified leading currents and presupposes an associative raw realization.

For the unlike-double-source part of B³AbarB, let a count A slots, b count B slots and d count barB slots. Put E=g binomial(b,2)−rho ab. Direct pair-power, derivative and measure counting gives

`omega_02=E+3a−2(b+d)+2h+1`.

For the full five-current face this is `3(g−rho)+2h−4`. All 15 current-containing subsets have strictly positive new unlike-source degrees on 2<t<3, g=rho=2/t, h=t; the source-only degree is 2h. The coordinator checks this exactly by expanding t omega at t=2+q with q>0: its coefficients are nonnegative and not all zero. The tight proper faces are `(t−2)(2t−1)/t` and `(t−2)(2t−3)/t`. These estimates complement the one-source and equal-source boundary proofs. They do not alone establish a factorization-compatible interaction map.

## 7. Independent audit and final scope

The completed [adversarial audit](referee_sl21_common_product_audit_report.md) independently reconstructs the full antichiral derivative polynomial, the source/tail powers, the spectator Laurent rule, and all fifteen unlike-source faces. Its fresh script has 65 passing checks. Its written Fourier argument supplies uniform partial-height limits, which are necessary for actual regrouping rather than just existence of an n-ary distribution.

The equal-root Vandermonde reduction applies both to canonical tails and to one- and two-source action kernels. At first mixed order, action and finite matching cancel insertion-wise, with every spectator and moving-test jet retained. Cauchy-versus-circular double-pole extensions have a local delta contact; the audit retains it as a collapsed source/tail term. The coherent whole-word families control the remaining regulator multiples. The unlike-source sector has an explicit positive physical patch and only transverse new divisors. This removes the proposed analytic obstruction in the bounded family; it does not calculate every second-order finite contact.

| Item | Verdict | Independent evidence and limit |
|---|---|---|
| Associative labelled Gaussian starting product | Correct | Additive bilinear covariance, cocycle and explicit derivative-composite regroupings, section 2 |
| Fixed-generator normalization transport creates the canonical mixed coefficient | Wrong | Its first-order braided coboundary is exactly zero, section 6 |
| Actual repeated-current mixed overlap | Correct | Full Laurent polynomial and edge/circle moments, section 3 |
| Both source-free canonical Serre quartics through mass order two | Correct | Tail-by-tail minority primitives and proper-face analysis, section 4 |
| A nonzero single-chamber residue proves a full-word anomaly | Wrong | Positive chamber integral is real, but coherent whole-word residues cancel, section 5 |
| Two equal-root tails force a one-barred cyclic obstruction | Wrong | Exact full-polynomial leaf decomposition, section 5 |
| First five-current mixed–Serre overlap in the canonical family with matched contour defect | Correct at the stated formal order | Uniform grouped boundary limits plus insertion-wise contact cancellation, audit sections 3–6 |
| Unlike-source faces create a divisor containing the physical surface | Wrong for this bounded family | All fifteen new degrees positive on 2<t<3; subsequent poles transverse, section 6 and audit section 7 |
| One common map reproduces every complete interacting word | Unclear / still open | Antisymmetric contact matching leaves ordinary dressing and symmetric contact values undetermined |
| Closed associative algebra for arbitrary barred words | Unclear / still open | The one-barred family is not closed under unrestricted multiplication |

The first gap is therefore **partially resolved, not closed**. We now have a concrete common canonical composition and a verified first mixed–Serre compatibility calculation. The next missing object is a single interaction/insertion map applied once to a composed cluster. It must reproduce the complete prescribed source-dressed words, including symmetric contacts and all composite/test jets, and be compatible with the relations. Its domain must then extend to products containing several barred currents. No independent dressing of binary blocks or unspecified finite-part convention supplies that object.

The construction provides a concrete candidate and several nontrivial verified overlaps. The unresolved full statement must not be inferred from any of the following: the name of the root system, a shuffle coefficient alone, an arbitrary abstract associative matrix representation, a common finite-part notation, or an assumed interacting vacuum functional.

## Integration and preservation

No book integration occurred. Independently verified and retained: the labelled product, derivative regroupings, canonical quartics, mixed contact bridge, leaf reduction, uniform grouping argument and unlike-face estimates. Rejected: fixed-generator normalization transport as a solution, the unreduced two-tail cycle as an obstruction, and an isolated chamber residue as a full-word anomaly. Modified in the synthesis: the earlier unresolved one-bar analytic question is now resolved for the bounded canonical family; the full interacting identification remains open. No theorem-label upgrade is recommended.

Earlier completed reports remain immutable. This team wrote only the five fresh common-product reports and the assigned scratch tree. The [preservation record](../referee_scratch/sl21_common_product_20261024/preservation.json) records the shared checkout advance: 1,078 of the 1,882 baseline paths are unchanged, 203 changed externally, and 601 no longer exist at their original paths during the concurrent repository reorganization. All 275 sl21 inputs are unchanged. Tracked and staged working diffs are empty at the final check. Thus we do not claim that the entire shared repository stayed unchanged. No book edit, commit, push, branch change, worktree creation or render was performed by this team. The fresh outputs are listed with hashes in the [output manifest](../referee_scratch/sl21_common_product_20261024/output_manifest.json).
