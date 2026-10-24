# Independent audit of symmetric interaction contacts

Date: 2026-10-24. Model: Codex (GPT-6). Reviewed commit: `ae19df62aaa82879da0f4600b74982d700d481b5`.

This audit obeys the editor's referee-only restriction. No existing file, book file or index was edited. No worktree, git mutation, commit, push or render was used. No subagents were spawned. Fresh calculations and output are in `referee_scratch/sl21_interaction_map_20261024/audit/`.

## Conclusion and verdicts

There is an explicit first-order **pair** interaction map with the required symmetric contacts. It has nonzero unary corrections and principal-value pair terms. A reflection-symmetric allocation of the finite mixed matching makes this map descend through the free mixed pair relation. This allocation is an additional normalization, not information already contained in the old contour defect.

The two-current, unlike-two-source forest can also be specified explicitly. Its genuinely distinct-source double contacts are locally integrable for `2<t<3`. The same-source double contacts require the already chosen joint pair extension, not a product of two delta distributions.

A further independent argument below establishes the first-order one-barred Serre test for the proposed simultaneous map. These facts do **not yet prove** that one map on the full canonical Gaussian quotient reproduces every source-dressed five-current word through mass order two. The remaining comparison is precise: the full rotational/moving-mask one-source union with several current slots must agree with the proposed simultaneous symmetric repair, including its symmetric full-diagonal distribution. Pair defects and positivity of the new two-source faces do not determine that pre-existing one-source extension. No explicit nonzero anomaly was found, and this is not a nonexistence result.

| Claim | Verdict | Independent evidence and scope |
|---|---|---|
| The old endpoint matching fixes both ordered products | wrong | It fixes their braided difference; adding the same aligned distribution to both leaves it unchanged |
| Each raw ordered BB action contact is half of the familiar defect with opposite signs | correct | Three separate source moments below, with both strip edges retained |
| A first derivative of the outer test gives a hidden scalar contact in this pair prescription | wrong | Its two moments cancel separately in each ordering |
| Action dressing can be placed solely in a changed state while fixing single-current values | wrong for the spatial-base comparison | Explicit rotational/spatial unary differences remain even at separated current positions |
| Unary corrections alone obstruct a linear map | wrong | A balanced joint pair correction cancels their image of the free mixed null |
| The symmetric pair repair contains only a local delta contact | wrong | Explicit nonlocal Pf and PV terms survive |
| The proposed source-specific differential base reproduces the full two-source pair Wick polynomial | correct | Independent symbolic calculation, including uncontracted field factors |
| Distinct-source unary contacts create a new pair full-diagonal ambiguity for `2<t<3` | wrong | Their kernel is proportional to `|x-y|^(2h-4)`, locally integrable at `h=t` |
| One may multiply two same-source current/source delta faces | wrong | Their pullback is nontransversal; it would contain the square of a source-time delta |
| The proposed first-order one-barred map annihilates the Serre test | correct in the common spatial boundary family | Uniform height-gap distribution bounds and the independently checked K0 variation; section 6 |
| A full five-current map through mass two follows from the preceding pair checks and positive new face degrees | unclear | Missing full symmetric one-source union comparison described in section 6 |

## 1. Definitions and the ordered first-order pair

I used the explicit vectors in the coordinator plan:

`a0=(ib/sqrt(2), b/sqrt(2), i)`, `a2=(ib/sqrt(2), -b/sqrt(2), -i)`, `t=b²`, `d=-(a0+a2)/t`.

Thus `a0²=a2²=-1`, `a0.a2=1-t`, and `d.ai=1`. Write `beta=t-1`, `v=d+a2`, `B=:C exp(dX):`, `C=a0.partial X`, and `U=:exp[v(X+Xbar)]:`. The common mixed cocycle is `nu=exp(-i pi/t)` and the aligned reverse product is multiplied by `s=nu²`.

The actual prescription in `briefs/referee_sl21_mixed_audit_report.md`, lines 159–185, uses the inner circle `|w|>epsilon`, outer heights `x+i sigma delta`, `sigma=+1,-1`, `delta=c epsilon`, `c>1`, and the Wirtinger distributions

`C_eta(z)=bar z/(|z|²+eta²)`, `P_eta(z)=-partial_z C_eta(z)`.

It is important that `P_eta` is not `-partial_x C_eta`.

Put the inner barred current at zero. Let `D=d.partial X`, `K_L=v.partial X`, `K_R=v.partialbar Xbar`. Direct expansion of the complete first normal-field jet and integration of the outer current gives

`i pi beta K_L sgn(sigma delta-y)/bar w²`

`+ pi beta² K_L w delta(sigma delta-y)/bar w²`

`+ pi beta t K_R delta(sigma delta-y)/bar w`.

The identities used here are `C+beta D=-K_L`, `D+a2.partial X=K_L`, and `beta a2.partialbar Xbar-Cbar=t K_R`; all three follow componentwise from the vectors.

The three remaining source moments, with a symmetric distant cutoff for the local coefficient, are

`integral_{|w|>epsilon} sgn(sigma delta-y)/bar w² d²w = pi sigma`,

`integral w delta(sigma delta-y)/bar w² d²w = i pi sigma`,

`integral delta(sigma delta-y)/bar w d²w = i pi sigma`.

For the first, integrate horizontal slices: the slice of `1/bar w²` is `2 sqrt(epsilon²-y²)/epsilon²` for `|y|<epsilon` and zero otherwise. Since `delta>epsilon`, the sign is constant on its support. The other two slices lie outside the circle, so their integrals are `i pi sigma`.

After the action sign and cocycle, the separately ordered local contributions are therefore

`C_raw,sigma = -i pi nu lambda2 sigma t beta (K_L+K_R) U`.

Their difference is `-2i pi nu lambda2 t beta partial_x U`; their aligned average is zero. This is a local contact statement, not a statement that the whole source integral has zero average.

The scalar kernel times the first outer-test jet has moment

`beta² [pi (i pi sigma) - i pi (pi sigma)] = 0`.

Thus this prescription produces no additional scalar derivative-of-test contact in either ordering. The scalar zeroth moment is also zero. Higher jets cannot create another local term at this degree. Field derivatives and test derivatives have been distinguished rather than set to constants together.

A valid additional normalization is

`C_finite,+ = +i pi nu lambda2 t beta partial_x U`,

`s C_finite,- = -i pi nu lambda2 t beta partial_x U`.

It cancels the local raw contribution in each aligned ordering. An arbitrary common aligned symmetric term would give the same old defect. Consequently the reflection-symmetric choice must be stated explicitly.

## 2. The unary terms cannot be omitted

The distribution identities are

`partialbar_z(1/z)=pi delta²(z)`, `partial_z(1/z)=-Pf(1/z²)`.

Since `partial_x=partial_z+partialbar_z`,

`-partial_x(1/z)=Pf(1/z²)-pi delta²(z)`.

Both the outer Cauchy prescription and the inner circular prescription use the rotational `Pf`. Relative to a spatial differentiated Gaussian base, each double pole therefore adds `pi delta²(z)`. Including the action coefficient and actual current-first cocycles gives

`U1(B)=lambda0 K0 + beta lambda2 K2`,

`U1(A)=-lambda1 F_A`,

and their barred versions, where

`Ki=:exp[(d+ai)X+ai Xbar]:`,

`F_A=:exp[-k X-T Xbar]:`, `T=a0+a2`, `k=d+T`.

Here `F_A` denotes the current tail, not the zeroform with both exponents `-k`. A source cutoff or coupling test `chi(w)` is retained; the unary delta evaluates it at the current position. Differentiated versions must keep derivatives of this test.

The coefficients follow from the B/Yi double-pole residues `-a0.ai`, namely `1` and `beta`. The A/V1 vertex has power `z^-2`; its cocycle is even, so the action supplies the minus sign. These single faces survive when the two current positions are distinct. They cannot be replaced by a correction supported only on the mixed current diagonal.

I initially considered whether the raw symmetric pair differed from the spatial base only on the total diagonal. That inference was incorrect; the unary identity above is the correction. No nonexistence conclusion in this report uses that discarded inference.

## 3. The full symmetric pair repair

The canonical tail contribution, with `r=x-y`, has local jet expansion

`nu lambda2 {2t beta U Pf(1/r²)`

`+ t[(t-2)partial U+t partialbar U] PV(1/r)`

`+ sigma i pi t² partial_x U delta(r)}`.

The full normal products, rather than this jet expansion, define the terms away from the diagonal. The two unary source corrections are `beta lambda2(K2 barB+B barK2)`, equal to `-beta/t` times this tail expression. Their symmetric part is therefore

`nu lambda2 {-2 beta² U Pf(1/r²)`

`- beta[(t-2)partial U+t partialbar U] PV(1/r)}`.

This nonlocal principal-value part is required. The unary odd contact is precisely the raw action contact of section 1. The balanced finite term removes it.

Let `S1` mean the common spatial-base source insertion, performed once on a complete Wick normal expression. On the mixed pair a constructive definition is

`I1(B barB)=S1(B barB)+1/2[U1(B)barB+B U1(barB)+s U1(barB)B+s barB U1(B)]`.

Use the same expression, multiplied by `s^-1`, for `I1(barB B)`. All products in this formula are their common Gaussian boundary distributions, with the inherited constituent labels; it is not a prescription to multiply scalar singular distributions independently. It retains the preceding Pf/PV terms and cancels the odd defect. In particular,

`I1(B barB-s barB B)=0`.

This is an explicit linear pair map through the mixed null. It also explains why the unary shift is a necessary correction, not a no-go theorem. At separated points all four terms reduce to the ordinary two unary insertions; no additional separated interaction is created.

A simultaneous symmetric rule on a larger word is not obtained by iterating this Jordan-type pair formula. Iterated symmetrization changes weights. Uniform averaging over the three possible barred positions in `B B barB`, for example, assigns a specified B/barB order probabilities `2/3` and `1/3`, not `1/2` and `1/2`. A proposed extension must specify one fullword rule and then prove its compatibility with relations.

## 4. A differential construction of the full two-source pair kernel

The spatial one-source base is

`G=1/[(x-w)(y-bar w)]`,

with the inherited opposite imaginary prescriptions. For the complete normal product

`N=:exp[dX(x)+dXbar(y)+a2 Phi(w)]:`,

put `D=d.partial X(x)` and `Dbar=d.partialbar Xbar(y)`. Then the exact, all-normal-jet identity is

`W=[C+beta D-beta partial_x][Cbar+beta Dbar-beta partial_y](G N)`.

The derivatives act on the whole `G N`; the added D terms cancel their action on N. The result is exactly

`G N [C+beta/(x-w)][Cbar+beta/(y-bar w)]`.

No derivative-of-normal-field or test term is silently dropped. Taking distributional derivatives also keeps derivatives of every smooth test and mask.

For the unlike two-source term define `zi=x-wi`, `ri=y-bar wi`, `i=0,2`, and

`G02=(z0 z2 r0 r2)^-1 |w0-w2|^(2h-2)`.

The complete polynomial is

`G02 N [C+1/z0+beta/z2][Cbar+1/r0+beta/r2]`.

An exact differential representative uses temporary independent source-edge variables:

`[C-partial_z0-beta partial_z2][Cbar-partial_r0-beta partial_r2] (z0 z2 r0 r2)^-1`,

multiplied by the full N and source-pair factor, with the edge variables then restricted to their actual values. Equivalently one may use derivatives on the full normal product with the corresponding compensating field derivatives. A single coefficient times the total `partial_x` would give the same weight to the two sources and is incorrect: the required coefficients are `1` and `beta`.

This calculation independently reconstructs the entire two-source pair Wick polynomial. It supplies the distributional starting point for a forest comparison, not a proof that every arbitrary edge-variable pullback is allowed.

## 5. Allowed and forbidden source forests

Replace each rotational double pole by its spatial double pole plus the unary delta. Expanding the exact polynomial makes three distinctions explicit.

First, the C factor is linear in `1/z0+beta/z2`. There is no term using both source double poles at the same current. Such a forest is absent algebraically, even before any extension question.

Second, distinct-source faces at distinct currents are transverse. The choice `w0=x`, `w2=y` leaves

`pi² beta [(x-y)(y-x)]^-1 |x-y|^(2h-2)`

`= -pi² beta |x-y|^(2h-4)`,

times the full collapsed normal field and fixed cocycle. The other source assignment is included once as well. For `h=t` and `2<t<3`, this is locally integrable; its radial degree in the remaining real variable is `2t-3>1`. Smooth normal/test jets improve the degree. It has no obligatory new delta extension at `x=y`. The full source-pair base has radial degree `2h-1`; after the two current derivatives the conservative degree is `2h-3`, still positive in this patch.

Third, assigning the same source simultaneously to B and barB is not transverse. Formally multiplying the two unary deltas would square the delta of the source imaginary coordinate. That expression is undefined. It must be replaced by the joint mixed pair extension of sections 1–3, with the second source, all normal fields, masks and their derivatives carried along. The prior antisymmetric bridge constrains this extension but does not by itself specify its symmetric part.

For a bounded pair construction one can choose the balanced joint primitive, use the exact differentiated base for its other-source dressing, and apply the common circular Taylor subtraction to the source-source diagonal. This is a concrete finite prescription. Its descendants are obtained from the full normal expression, not just from its scalar residue. The genuinely distinct-source double faces add no primitive ambiguity in the displayed patch. This is the positive pair result proved here; I do not infer an arbitrary-arity associative map from it.

The action expansion is applied once to the complete word. Its two unlike sources have coefficient `lambda0 lambda2/pi²`, while two equal sources have `lambda2²/(2pi²)`. Applying an interaction expansion anew at every binary multiplication would double count these terms. If a vacuum normalization is used, its denominator is scalar only after a state is specified. A Wick-valued insertion map generally has a nontrivial value on the unit; it is not automatically a normalized expectation or an algebra homomorphism.

## 6. Exact remaining descent question

The repaired pair map is well-defined on the mixed pair quotient. There is also a positive first-order one-barred Serre test for the proposed simultaneous map. This is stronger than checking the marginal scalar alone.

Let S0 be a bare Serre quartet with the same smooth species test repeated in its three equal-species slots, and let D be the extra barred current. Keep the source height tau as a real integration variable. At fixed nonzero source height and with D separated, the quartet limit can be taken first; the canonical quartet is zero against those smooth external contractions.

To interchange that grouped limit and source integration, divide the tau axis into its finitely many ordering intervals relative to the current heights. In each interval the full holomorphic contraction graph is directed acyclic and the sole antichiral source-to-D edge is a leaf. The standard Fourier pushforward proof is uniform in all nonnegative height gaps: a discrete strictly increasing potential on the vertices bounds the sum of nonnegative edge momenta by a constant times the incidence image. The constant depends on the finite vertex order, not the numerical gaps. Each distributional seminorm is therefore bounded uniformly as any subset of gaps closes. There are only finitely many such orders. The smooth normal-field and test coefficients have uniformly bounded seminorms on the compact source support.

Consequently the source-height intervals between the collapsing quartet heights contribute at most their lengths times a fixed distributional seminorm; they tend to zero. On the remaining intervals the grouped quartet zero and dominated convergence apply. All differentiated current kernels use spatial derivatives, which do not differentiate tau. They introduce no additional delta(tau). This establishes `S1(S0 D)=0` for the common spatial boundary extension. The rotational delta(tau) terms are precisely the separately identified unary corrections, not a missing part of this argument.

For the unary variation, K2 and the even-root tail are the canonical tail variations. The only new replacement is K0. Independently, `-a0.(d+a0)=0`, and with K0 at the third B position the scalar coefficient of the two remaining derivative currents is

`[-1/(z-v)+1/(z-a)] [1/(z-v)+1/(v-a)]+1/(z-v)^2=0`.

Every nonzero regular-field or test jet improves the marginal full degree by at least one. The K0-to-B and K0-to-A powers are respectively `rho-1` and `1-rho`; their integer shifts and the source cocycle give exactly the same current-label exchange phases as K2. Thus the same full-kernel shuffle applies, with no new pole on g=rho at the marginal face. In A³B the K0 replacement is nonmarginal. The summed unary variation `U1(S0)` vanishes. The simultaneous extreme-order repair therefore also annihilates S0 D. This proves the stated first-order one-barred Serre test, retaining normal and test jets; it is not a statement about an unpolarized quartic with unrelated tests in its repeated slots.

The full mass-two null-ideal statement still needs an analogous comparison for the two-source forests and every proper face before multiplication can be transported. The following distinction remains essential. A map defined by this common spatial family and simultaneous repair has the preceding first-order Serre property. Identifying its full symmetric multi-current distribution with the old rotational moving-mask continuation is an additional comparison; an antisymmetric Ward bridge alone does not supply it.

The existing first-order Ward bridge fixes an antisymmetric contour difference. The full `B³ A barB/Y2` face has degree `3(g-rho)-1`; its first Taylor level can therefore support a distribution invisible to a constant-test or antisymmetric calculation. It is not enough that the *new* unlike-two-source faces have strictly positive degrees. They can inherit the already selected extension on this proper one-source face.

There is a plausible route to finish: in a common large-Re-g, negative-Re-rho initial domain, compare the full rotational moving-mask expression against the simultaneous fullword symmetric repair, not an iterated Jordan rule. Prove that every residual union vanishes there with all normal/test jets; then establish that its continuation is the same complete-word boundary family used in the canonical Gaussian quotient. This would fix the potentially dangerous symmetric distribution and permit the no-new-two-source-face argument to be used. The present audit has not completed that comparison.

A candidate based on a forest algorithm can certainly be defined on unreduced labelled words. Definition on that larger space does not prove descent. Likewise, for a linear invertible map `T=Id+lambda Z+...`, the transported first-order product correction is a Hochschild coboundary and its first associator vanishes algebraically. That is useful only after Z is a well-defined map on the actual quotient and matches the prescribed kernels; it cannot establish those premises.

Recommended replacement for any full-closure assertion at this stage:

> The balanced interaction prescription constructs the symmetric mixed pair map and its unlike-two-source forest through mass order two. Its unary corrections and nonlocal principal-value terms are fixed explicitly. The common spatial/simultaneous prescription also preserves the first-order one-barred Serre test. Full mass-two descent on the five-current canonical quotient, together with agreement of its symmetric multi-current forests with the prescribed moving-mask continuation, remains to be proved. No obstruction has been demonstrated; pair Ward identities and positive new-face degrees alone do not prove that remaining assertion.

## Evidence and integration

The fresh `check_contacts.py` reconstructs the vectors, all-normal-jet differential identities, unequal source weights, allowed doubleface terms, ordered moments, derivative-test cancellation and mixed-null repair. It imports no prior research script. Its retained output has **32 PASS, no FAIL**. The separate fresh `check_first_order_serre.py` adds **10 PASS, no FAIL**, verifying the K0 contractions, exchange shifts, zero marginal coefficient and first surviving degree. Four numerical circle moments were independently evaluated at `c=sqrt(2)` and `c=e`, for both signs, with 40-digit arithmetic; their reported errors are zero at that precision. Symbolic checks report zero residuals. These are checks of the displayed finite computations, not a proof of the unresolved fullword map.

The source hashes and checkpoint are retained in the same fresh scratch directory. No integration into the book was performed. Verified: ordered pair contacts, symmetric Pf/PV repair, unary shifts, exact two-source pair polynomial and distinct-source integrability. Corrected: the tempting omission of single current/source faces. Rejected: multiplying overlapping same-source deltas, state-only dressing with fixed unary values, and inferring quotient descent solely from a forest algorithm or face degrees. Additionally verified: first-order one-barred Serre annihilation for the proposed common spatial/simultaneous map, with its uniform-height argument above. Unresolved: the mass-two quotient comparison and full identification with the old symmetric multi-current moving-mask continuation. No theorem upgrade is recommended.

The completed construction report was then read as a proposed formula, not evidence. Its single source expansion, non-unital operator-valued status, balanced whole-cluster rule and explicitly bounded mass-two verdict agree with the audit distinctions above. The coordinator synthesis section 6 uses the same uniform source-height argument accepted here. This report and its computation outputs are now frozen.
