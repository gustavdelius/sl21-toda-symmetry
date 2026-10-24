# Coordinated sl(2|1) Toda nonlocal-symmetry investigation

Date: 2026-10-24. Model: Codex (GPT-6). Reviewed branch: `main`; commit: `745586e6a6f7bca5ed91a1e9cf89d6ef5df387f0`.

The commit does not identify the untracked evidence. SHA-256 hashes of 1,824 pre-existing files, including the inherited reports and scratch evidence, are recorded in [preexisting_hashes.json](../referee_scratch/sl21_symmetry_coordinator_20261024/preexisting_hashes.json). The original external Bernard–LeClair PDF has a separate [hash record](../referee_scratch/sl21_symmetry_coordinator_20261024/coordinator/external_source_hash.json). The editor's coordinating exception authorized delegation; there were initially one, then at most two active workers. Every worker was instructed not to delegate. All new work stays in the authorized scratch tree and new reports.

## Outcome

**Unresolved as a quantum symmetry algebra.** Neither affine sl(2|1) nor affine sl2 has been established for the original root-system quantum theory. The investigation makes a substantive finite-order extension and narrows the remaining problem; it does not convert local contour identities into a global algebra by assumption.

The inherited first-order triple identity survives independent audit. At double-Yukawa order an exact calculation shows why the obvious separated-face cutoff argument fails. An auxiliary source exponent supplies a common continuation, with a calculated regulator Ward defect proportional to h−t and corrected counterterms. The mass-squared projected identity survives the independent final adversarial audit **as a formal continued distribution in that prescription**. The companion leading projection is independently calculated in a complementary support arrangement. These products still do not determine an arbitrary-test associative multiplication or its compatibility with every Serre product.

Independent attempts at all-order and global closure produce precise conditional statements and a restricted obstruction: ordinary local perturbative anomalies are excluded under a stated class of assumptions, while the naive real sine-Liouville bosonic integration cycle diverges. Complex cycles and a fermionic reference construction are not ruled out. Neither has been matched here to the full Toda Ward endpoints on interface sectors.

No book changes or theorem-label upgrades are proposed.

## 1. Claim and dependency ledger

| Claim | Verdict | Fresh evidence and exact ceiling |
|---|---|---|
| Action normalization, root contractions, derivative residues and first contact | Correct | Source inspection and coordinator's direct Gaussian calculations |
| First triple coefficient linear in lambda2 | Correct within the inherited prescription | Independent kernel/face audit; fixed compact support and Cauchy order |
| Unique physical selection of t² by that test | Wrong | The calculation transports any consistently chosen binary coefficient |
| Fixed-separation source-pair remainder is uniform at double-Yukawa order | Wrong | Exact nonzero omitted-disk formula below |
| Lower action/current/endpoint subtraction chain | Correct | Independent radial/angle calculations, including both endpoint descendants |
| Common h-family mass-squared projected identity | Correct as a formal continued identity | Explicit channels, mixed contacts and (h−t) Ward variation; independent final audit; no original-theory operator assertion |
| Leading barA companion | Correct within its restricted prescription | Actual pole calculation, independently recomputed by coordinator |
| All mixed/Serre products belong to one associative product | Unclear | No complete regrouping/overlap calculation is available |
| No further ordinary local anomalies at irrational t | Correct conditionally | Full integer argument for both currents, not just a finite scan |
| Source reflection normalization and Gaussian transport | Correct conditionally | Direct zero/pole and vector-matrix calculation; reference family only |
| Naive real positive-coupling reference-boson path integral | Wrong as a convergent completion | Explicit open-region exponential runaway |
| Original-theory endpoint limits, interfaces, common charge domain | Unclear | Not constructed by reflection values or homogeneous particle sectors |
| Quantum-affine algebra/Hopf action, image or faithfulness | Unclear | Required product, endpoint and string/Cartan identifications remain missing |

The detailed [claim ledger](../referee_scratch/sl21_symmetry_coordinator_20261024/claim_ledger.md), [issue ledger](../referee_scratch/sl21_symmetry_coordinator_20261024/issue_ledger.md), and [checkpoint](../referee_scratch/sl21_symmetry_coordinator_20261024/checkpoint.md) distinguish completed calculations from unresolved links.

## 2. Independent intake: what was actually verified

The definitions used throughout are

\[
a_0=(ib/\sqrt2,b/\sqrt2,i),\quad a_2=(ib/\sqrt2,-b/\sqrt2,-i),\quad
t=b^2,\quad T=a_0+a_2,\quad d=-T/t,\quad k=d+T.
\]

Write rho=2/t, nu=exp(−i pi/t), s=nu², C=a0.partial X. The even currents are J_A=exp(−d.X), J_B=C exp(d.X), with density J+H and charge normalization (2pi i)^−1. The tails are H_B=−t lambda2 K2+lambda0 lambda2 epsilon^(2t)L and H_A=−gamma lambda1 exp(−k.X−T.Xbar), where gamma=t/(1−t), K2=exp[(d+a2).X+a2.Xbar], L=C exp[k.X+T.Xbar]. The action is pi^−1 integral(lambda0 Y0+lambda2 Y2+lambda1 V1+kappa P), where Yi=exp[ai.(X+Xbar)], V1=exp[−T.(X+Xbar)], P=exp[T.(X+Xbar)], kappa=−lambda0 lambda2 epsilon^(2t)/t, and the independent finite P coupling is zero. The endpoint fields are U=exp[(d+a2).(X+Xbar)], F=exp[k.(X+Xbar)], F_A=exp[−k.(X+Xbar)]. Lambda1 is an independent mass-order-two coupling; lambda0 and lambda2 have mass order one.

For labels n=(n_L,n_R,n_0,n_2), the cocycle is epsilon(n,m)=(−1)^[(n_L+n_R)(m_0+m_2)] exp[−i pi(n_L m_R−n_R m_L)/t]. These local conventions do not specify the global spin/compactification sectors.

The covariance is X(z)X(w)~−log(z−w), with bilinear scalar products and ordinary Euclidean area. The explicit vectors give a0²=a2²=−1, a0.a2=1−t, ai.d=1, d²=−rho, rho=2/t. For C=a0.partial X,

`C(z) exp(q.X(w)) ~ −a0.q exp(q.X(w))/(z−w)`.

Thus the A spectator contributes +1/(z−a), and Y2 contributes (t−1)/(z−w). The mixed Laurent coefficients are (t−1)V and −V', with V=(a−w)^(1−rho). The field-jet identity is a0+(t−1)d=−(d+a2). Independent circle/strip integration gives I=pi and both edge moments −pi. Consequently the action coefficient −t(t−1) plus the two tails t² gives raw t. Adding the specified finite matching t(t−1) gives t².

The first triple report's proper faces and nested bound were rederived in the [stage-0 audit](referee_sl21_coord_stage0_audit_report.md). After resolving the proper mixed face, the new AB face has relative degree −rho, the full action cluster 1−rho, and the nested remainder is bounded by epsilon integral R^(−1−Re rho)dR. It vanishes for Re rho<0. The endpoint cocycle must be oriented correctly: epsilon(U,A)/epsilon(A,U)=−s. Combining positive-imaginary-separation oscillator values yields U(L)A=s²AU(L), U(R)A=AU(R).

The restricted result remains

\[
\left.[B_gA_f-s^{-1}A_fB_g,\bar B_I]\right|_{\lambda_2}
=(s^{-1}-s)A_f\frac{\nu}{2\pi i}t^2\lambda_2U(L).
\]

Here f is supported strictly inside I, g=1 near I, and current poles are Cauchy-extended before the inner circular hole. With barB at height zero and delta=c epsilon, c>2, the four words BAbarB, ABbarB, barBBA, barBAB use (tau_B,tau_A)=(2delta,delta), (delta,2delta), (−delta,−2delta), (−2delta,−delta), respectively. Infrared cutoffs are fixed and other spectators separated before continuation. These qualifications are part of the result.

The coordinator's [intake checks](../referee_scratch/sl21_symmetry_coordinator_20261024/coordinator/intake_checks.py) and [output](../referee_scratch/sl21_symmetry_coordinator_20261024/coordinator/intake_checks.out) also reconstruct finite Gaussian transport and the source reflection zero/pole. An initial checking-code output is retained: simplification before substituting b²=t and reversing the cocycle ratio to the actual U/A order resolve its failures. No input formula was changed to make a check pass.

## 3. Double-Yukawa order: an actual overlap and a common extension

### Full coefficient list and a decisive failed shortcut

The distinct-source term in exp(−S) is +lambda0 lambda2/pi²; the two cross terms cancel 1/2!. There are also two first-tail/source terms with coefficient +t lambda0 lambda2/pi before cocycles, two improvements +lambda0 lambda2 epsilon^(2t)L, and the P insertion +lambda0 lambda2 epsilon^(2t)/(pi t). The first-tail square is lambda2² and is not this coefficient. The A tail is proportional to the independent lambda1.

The bare ABbarB source kernel is, with beta=t−1,

\[
\frac{\lambda_0\lambda_2}{\pi^2}
\frac{(a-z)^{-\rho}(a-v)(a-w)|v-w|^{2t-2}}
{(z-v)(z-w)(\bar u-\bar v)(\bar u-\bar w)}
\left[C(z)+\frac1{z-a}+\frac1{z-v}+\frac\beta{z-w}\right]
\left[\bar C(u)+\frac1{\bar u-\bar v}+\frac\beta{\bar u-\bar w}\right]\mathcal N.
\]

The ordered cocycle is +1; N contains the sum of the five insertion momenta. The [subtraction report](referee_sl21_coord_stage1_subtractions_report.md) derives the other signs independently.

Set the source positions to r,0 and suppress smooth field jets. A direct expansion gives

\[
f(z;r)=\frac{(z-r)^{-1}+(t-1)z^{-1}}{(z-r)z}
=z^{-3}\sum_{n\ge0}(n+t)(r/z)^n.
\]

Subtracting the scalar source-pair P term, the omitted disk is exactly

\[
\mathcal E=\frac{\pi\epsilon^{2t}}{(zv)^3}
\frac{q[(1+t)-tq]}{(1-q)^2},\qquad q=\epsilon^2/(zv),
\]

where v denotes the antichiral current coordinate in this formula. With A separated the actual factor is a(a−z)^(1−rho) times E. At fixed z,v this is O(epsilon^(2t+2)); at z=epsilon Z,v=epsilon V it is generically O(epsilon^(2t−6)). Extra finite Taylor terms do not repair uniformity on that corner. This disproves a shortcut, **not** the final four-word identity: tails and union-face terms have yet to be combined.

The coordinator independently recomputed its Laurent coefficients, radial primitives and nonzero scaling coefficient in [extension_checks.py](../referee_scratch/sl21_symmetry_coordinator_20261024/coordinator/extension_checks.py) and [output](../referee_scratch/sl21_symmetry_coordinator_20261024/coordinator/extension_checks.out).

### What the source-exponent construction adds

The [kernel report](referee_sl21_coord_stage1_kernels_report.md) replaces the unlike-source factor by |v−w|^(2h−2), keeping t in contractions. The same factor follows source constituents into tails and endpoints. This gives a common initial domain Re h>2, Re rho<0 for the triple, with first-order Cauchy faces fixed first. It avoids the impossible raw physical conditions Re t>0 and Re(2/t)<0.

Radial integration fixes

\[
\kappa_h=-\frac{\lambda_0\lambda_2}{h}\epsilon^{2h},\quad
H_{B,h}^{(2)}=\frac th\lambda_0\lambda_2\epsilon^{2h}L,
\]
\[
[U]_h=U-\frac{\lambda_0}{h-1}\epsilon^{2h-2}F
-\frac{\lambda_0}{h}\epsilon^{2h}C\bar CF.
\]

The finite binary addition must use t(t−1)lambda2[U]_h. Holding the old improvement coefficient equal to one off h=t is incorrect.

This regulator does not preserve the Ward identity off the physical diagonal. Its *actual* pole coefficients are

`D0=F0, S0=partial F0−(h−t)F0/q`,
`D2=(t−1)F2, S2=−partial F2+(h−t)F2/q`,

with q=current minus remaining source. The regulator-dependent term is explicit before integration. In jet variables C,D=d.partial X,

`N=C'+t D'/2−tCD−t²D²/2`,
`R_P=−tCD+t²D'/2−t³D²/2`,
`partial L=C'+(1−t)CD`.

The coordinator independently checked

\[
N-R_P/t=\partial L,\qquad
N-R_P/h-(t/h)\partial L=(h-t)N/h.
\]

Thus the remaining bilocal Ward defect has the form

\[
(h-t)\lambda_0\lambda_2\left[
\pi^{-1}\int\frac{F_0-F_2}{r}\,d^2r+epsilon^{2h}N/h\right].
\]

The complete Stokes reduction gives the projected defect minus the dressed endpoint as (h−t)Q_h. Sectorwise Taylor/Mellin continuation of that **whole** family fixes overlapping faces together. The new divisors have form ah+b rho+j=0 with a nonzero; after h=t,rho=2/t none is identically zero as a function of t. Isolated resonances must be excluded. First-order faces remain the separately fixed Cauchy contacts. The independent final audit accepts this route after spelling out the masked product rule; a vanishing prefactor alone would not justify its restriction. The coordinator independently checked that rule in [regulator_product_rule.py](../referee_scratch/sl21_symmetry_coordinator_20261024/coordinator/regulator_product_rule.py) and its [output](../referee_scratch/sl21_symmetry_coordinator_20261024/coordinator/regulator_product_rule.out). Every derivative of R_h=q^(h−t)qbar^(h−t) carries h−t. For example, partial²(R_h M Q)−R_h partial²(MQ)=2(partial R_h)partial(MQ)+(partial²R_h)MQ. Mask derivatives are retained in M, including at their intersections. The differences between the h-counterterms and the physical coefficients transported by epsilon^(2(h−t)) also carry h−t: the action difference is (h−t)/(ht), the tail difference is −(h−t)/h, and the two endpoint differences have denominators (h−1)(t−1) and ht. Terms without a derivative on R_h reproduce the fixed first-order contacts with their source spectator. All remaining new union faces have positive relative degree in the initial domain, including the most singular proper face BD02 of degree 2h−3. This treats distributional mask differentiation; it is not an assertion that those boundaries can be dropped on the physical negative-t diagonal.

### Explicit local mixed-contact bridge

The separate [contact-bridge report](referee_sl21_coord_stage1_contact_bridge_report.md) tests the most vulnerable local step directly. Put the inner barB at zero, the extra Y0 source at v, eta=h−t, R=1/v+1/vbar, and Q=|v|^(2h−4):U(0)Y0(v): with the normal exponential and no repeated scalar Wick factor. Dividing by 2 i nu lambda0 lambda2, the actual node-2 action contact is (t−1)[t partial_x Q+eta RQ], while the two tails sum to −t[t partial_x Q+eta RQ]. Their sum is −t partial_x Q−eta RQ. At the node-0 face the derivative cancels and the remaining contact is −eta RQ0. Full antichiral field jets, as well as the strip-edge derivative, are included. The coordinator independently reconstructs these formulas in [contact_bridge_checks.py](../referee_scratch/sl21_symmetry_coordinator_20261024/coordinator/contact_bridge_checks.py) and [output](../referee_scratch/sl21_symmetry_coordinator_20261024/coordinator/contact_bridge_checks.out).

For definiteness the off-diagonal finite matching is transported insertion by insertion. An added alpha lambda2 U then contributes −alpha[partial_x Q+eta RQ] in the same units. With alpha=t(t−1), the total derivative coefficient is −t². Defining the finite contact instead as the derivative of the whole regulated zero-form removes its −alpha eta RQ term; that is a different off-diagonal extension. Their difference has the factor eta and no divisor containing h=t in the full family, so both have the same generic physical diagonal. This distinction is explicit here and must not be silently conflated in an off-diagonal calculation. The narrow bridge alone does not prove the union continuation; that is the separate audited argument above.

At a generic real −1<t<0, the new A-bearing degrees are positive. After the lower mixed face is extended through the requisite jets, the spectator bound is epsilon^(2t+2) integral R^(−3−rho)dR, with rho<−2; the complementary corner has positive degree 2t−rho. Thus A introduces no new primitive divergence **conditional on the compatible lower extension**.

The accepted formal outcome is the first-order formula above with t²lambda2 U(L) replaced by its full order-02 renormalized insertion, including its negative Y0 action dressing and both displayed endpoint subtractions. Explicitly, coefficient extraction from the complete regulated correlator gives

\[
\left.[B_gA_f-s^{-1}A_fB_g,\bar B_I]\right|_{02}
=(s^{-1}-s)\frac{\nu}{2\pi i}
\left.A_f\,t^2\lambda_2[U(L)]\right|_{02}.
\]

Both sides mean the same whole h-family continued and then restricted to generic nonresonant h=t; 02 denotes the coefficient proportional to lambda0 lambda2, with all its source dressings. This is a specified continued product; equality with a naive negative-t sharp-cutoff limit is not claimed. A finite change delta c in the leading binary coefficient propagates as (s^-1−s) A_f nu delta c lambda2[U(L)]/(2pi i). The test does not select a physical normalization.

## 4. Common products, companion projection and Serre ceiling

The companion calculation takes a compact B spectator inside the barA interval, with the moving A test equal to one nearby. Its actual moving-A double/simple coefficients are V and V'/(1−t), where

`V=(a−w)^(2−rho)[C(a)+(1−t)/(a−w)]`

times the full fused normal product. The field-jet identity k/(1−t)=d is essential. Circle/edge constants give action coefficient gamma=t/(1−t), tails gamma², and the prescribed finite addition −gamma leaves gamma². The coordinator checked these coefficients directly in [companion_and_algebra_checks.py](../referee_scratch/sl21_symmetry_coordinator_20261024/coordinator/companion_and_algebra_checks.py), with [output](../referee_scratch/sl21_symmetry_coordinator_20261024/coordinator/companion_and_algebra_checks.out). The completed [companion report](referee_sl21_coord_stage2_companion_report.md) obtains

\[
[B_fA_g-s^{-1}A_gB_f,\bar A_I]|_{\lambda_1}
=(1-s^{-2})B_f\frac{\nu}{2\pi i}\gamma^2\lambda_1F_A(R).
\]

Its support assignment is deliberately different from the barB projection: B is compactly supported inside I and the A test is constant nearby. The left endpoint cancels. Its nested remainder is bounded by epsilon integral R^(−1−Re rho)dR, so no primitive contact remains in Re rho<0. This is the lambda1 coefficient, not every term of mass order two. The actual first outer-test derivative term cancels: the action gives I+Jw=0 and the two tail contributions cancel. The minimal AA matching adds no separate delta-prime term.

A concrete density direction invisible to the original constant-test projection is lambda2 U(y) partial_x delta(x−y). Its scaling dimension equals that of the two B densities, since dim(lambda2 U)=rho and dim(delta')=2. Its projection against an outer test g is −lambda2 g'(y)U(y); g=1 hides it. Compatibility with full covariance and braided exchange constrains such terms and must be imposed, not assumed. This example locates missing data; it is not a proof that every coefficient is admissible in every prescription.

The bridge calculation also determines this first BB test jet in the minimal family: its action term is proportional to (t−1)²(I+Jw)=0, the two tail double-pole terms cancel, and the chosen finite matching introduces no extra test derivative. Thus its actual coefficient at order lambda2 is zero, as is the AA coefficient at order lambda1. These are two calculated local jets. Higher-order derivative contacts and regrouping remain uncomputed; a constant-test identity alone would not have fixed even these leading jets.

The inherited quartic reports define separately continued four-fold products with auxiliary like-current exponent g and source exponent h. Their vanishing permutation factors are not enough to establish identities in an iterated common multiplication. The new mixed triple has no like-current pair and does not compute mixed/Serre overlaps involving such a pair. In particular, the finite canonical mixed contact has to be extended through the same higher forests, including derivative-test contacts, before the separate identities can be assembled. No claim of common-product associativity is accepted here.

## 5. Distinct all-order and completion attempts

The complete independent derivations are in [all_order_and_completion.md](../referee_scratch/sl21_symmetry_coordinator_20261024/coordinator/all_order_and_completion.md).

**Ordinary local anomaly classification.** For n0,n2 Yukawas and e even-root insertions set a=n0−e,c=n2−e,S=a+c,F=n0−n2,M=n0+n2+2e. Gaussian weights give

\[
\ell_B=1-M/2+S-F^2/2+t(e-ac),\quad
\ell_A=-M/2-S-F^2/2+t(e-ac),
\]
\[
\bar\ell=1-M/2-F^2/2+t(e-ac).
\]

At irrational t, integrality requires e=ac. Then bar ell=1−S(S+1)/2, ell_B=1−S(S−1)/2, ell_A=−S(S+3)/2. Nonnegative levels and insertion counts leave only the already calculated single Yukawas for B and single even-root insertion for A. This is a full integer proof under mass analyticity, Gaussian local fields, momentum/fermion-charge conservation, homogeneous subtraction and no independent finite P coupling. It does not show every nonlocal Ward breaking belongs to that class. A locality/induction principle for disorder-line and endpoint products remains missing.

**Screening and reference-theory routes.** Primary-source checks did not supply the missing bridge. Litvinov–Spodyneiko state “then it also commutes with the dressed screening field”; the claim concerns a chiral W-algebra commutant. Flandoli–Lentner require “provided that (m_ij) is subpolar”. These statements cannot be substituted for the present massive product calculation. Delius–MacKay explicitly treat an “n-component bosonic field” with ordinary exponential interactions. Bernard–LeClair's later scaling arguments for sine-Gordon do not themselves establish their analogue for this root-system action. Exact source locations, normalization conversions and inspected quotations are in the [primary-source audit](../referee_scratch/sl21_symmetry_coordinator_20261024/coordinator/literature_audit.md).

**A real reference cycle fails.** With t=−2zeta², alpha²=1+zeta², and lambda0=lambda2=lambda>0, the real local sine-Liouville potential is

`2lambda exp(zeta phi)cos(alpha theta)+lambda1 exp(−2zeta phi)`.

On an open theta interval with cos(alpha theta)<−1/2 it is bounded above by −lambda exp(zeta phi)+lambda1 exp(−2zeta phi). Therefore its real zero-mode Euclidean integral diverges as phi→+infinity. Finite-regulator gradient terms grow at most quadratically along this runaway and do not cure it. This obstructs that particular naive bosonic cycle only. It does not exclude complex cycles, the fermionic/cosh/Thirring reference realization, or the original action's different real slice.

**Reflection is a normalization test, not an interface construction.** The original Fateev page (131)–(132) was visually inspected. Independently taking its simple reflection zero times the reflected identity pole gives the inherited G_U, including the sign. H=diag(c,c,0) gives Ward coefficient factors exp(2c), exp(c), fusion factors exp[c(1−t)], exp(−ct), and common fixed-action endpoint-amplitude factor exp(c/t). These are stringent checks on a proposed match. They do not prove its all-order existence.

Even vanishing all positive-particle form factors leaves an identity admixture undetermined without full normalization; eta G(lambda)I has the required dimension/character and is nonanalytic in the couplings. The specified reference family may exclude that change. The unproved map from the full Toda Ward endpoint to that family cannot. Crossing and completeness on one homogeneous particle sector do not construct H_(l,r), l!=r, or the charge intertwiners. Weak endpoint limits are needed on a common invariant interface domain, with ultraviolet/infinite-volume limits controlled. Finite perturbative remainders cannot supply them because the candidate VEVs have nonanalytic mass powers.

## 6. What algebra would follow conditionally, and what the anchor actually gives

Directly splitting the left-anchored canonical field at c gives, for x>c,

`varphi_L(x)=varphi_c(x)−(1/2) integral_L^c dotPhi`.

Every component with chiral-minus-antichiral momentum qd therefore acquires the regional string

`Omega_q(L,c)=exp[−q d.integral_L^c dotPhi/2]`.

The formal regional rule is Q_q(L,R)=Q_q(L,c)+Omega_q(L,c)Q_q(c,R). On genuine regional tensor spaces this has coproduct shape Q tensor 1+Omega tensor Q; three-region grouping follows from additivity of the integral. Initially Omega is a canonical-momentum string, not the topological field-difference operator K=nu^N. Their action on the chosen state/field sectors must be identified before writing a quantum-group coproduct.

If one *assumes* a common associative invariant domain, all required mixed and barred/unbarred quartic identities, K A=s^-1 A K, K B=s B K, matching barred weights, and nonzero diagonal constants a,b with defects a(1−K²), b(1−K^-2), there is an explicit presentation conversion:

\[
K_0=K^{-1},\quad K_1=K,\quad E_0=A,\quad E_1=B,
\]
\[
F_0=\frac{\bar A K^{-1}}{a(\nu-\nu^{-1})},\qquad
F_1=\frac{\bar B K}{b(\nu-\nu^{-1})}.
\]

For example [A,barA K^-1]=a(K^-1−K); cross commutators cancel. Moving K through each barred quartic gives one common factor, so the negative-generator relations follow under the assumptions. This is the level-zero quantum affine sl2 presentation, with deformation parameter nu, without the degree generator. It would provide an action, possibly through a proper quotient. It would not prove faithfulness, identify the image with the presentation, fix the physical sign convention, or exclude additional generators/symmetries. The assumptions have not been established for this Toda theory; this conditional conversion is not its symmetry identification.

## 7. Independent adversarial audit and stopping criterion

The [completed final adversarial audit](referee_sl21_coord_final_adversarial_report.md) independently reconstructs the decisive A-dressed Laurent coefficients, all eight new collision degrees, the second-derivative mask product rule, and the Mellin-divisor restriction. Its verdict is **correct within the specified formal common-family prescription** for the mass-squared projection. It found no hidden pole containing the physical diagonal. It independently verifies the conditional anomaly integer proof, the restricted real-cycle obstruction, the spatial anchor string, and the conditional presentation conversion. Its 16 PASS lines record reproducible calculations, not evidence for a Hilbert-space construction. An initial symbolic-programming diagnostic is retained.

The main clarification was substantive: every mask must remain inside the differential product rule, and the finite matching must be dressed before continuation. The completed kernel report's compressed H8 exposition is supplemented by section 3 here and section 3 of the final audit, without rewriting that completed pass. The narrow bridge independently supplies actual separated-source coefficients and makes the off-diagonal finite-matching convention explicit. No disagreement requiring a formula reversal remains. Stronger claims about a direct negative-t cutoff limit, arbitrary-test multiplication or global algebra remain rejected or unresolved, as listed above.

This investigation ends with outcome 3 of the plan: sharply delimited unresolved after distinct serious attempts. The remaining minimal links are: (i) one product extending the mixed and Serre families with compatible regrouping and all needed tests, (ii) a Ward/relations induction mechanism covering its nonlocal contacts, and (iii) an explicitly matched quantum completion with interface sectors and weak full-endpoint limits. The screened-commutant, ordinary-local-anomaly, real-cycle, fermionic-reference, reflection-normalization and canonical-anchor routes were all examined; their exact ceilings are recorded. Another repetition of the same local residues would not establish any of these missing links. None can be replaced by a scalar VEV, an algebraic matrix witness, or the root-system name.

## Integration and preservation

No book integration occurred. Verified independently: action conventions, first triple contact and collision control, double-source subtraction coefficients, a nonuniform overlap coefficient, the formal mass-squared projection with its local contact bridge, the companion projection, the two leading test jets, reference normalization and Gaussian transport, and the stated conditional anomaly and algebra-conversion results. Modified exposition: full-mask differentiation and off-diagonal finite matching are now explicit. Rejected: uniform use of fixed-separation remainders, an automatic real reference-boson measure, and inference of a global algebra from local/finite-order identities. No theorem-label upgrade is recommended.

All worker passes are complete. The final [preservation record](../referee_scratch/sl21_symmetry_coordinator_20261024/final_preservation.json) compares the 1,824 snapshotted files and the external original PDF, and records the reviewed commit and empty tracked/staged diffs. The [output manifest](../referee_scratch/sl21_symmetry_coordinator_20261024/output_manifest.json) hashes every fresh report and scratch output except itself. No pre-existing file or report was edited; no book edit, commit, push, branch change, worktree creation or render was performed.
