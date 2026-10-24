# The first two-barred grouping and collision test

Date: 2026-10-24. Model: Codex (GPT-6). Reviewed commit: `483b0b9b9cc222e67cae61a549e7145151e905e4`, branch `main`.

This is a fresh referee-only continuation requested by the editor. Two independent workers reconstructed the construction and audited its collision and grouping claims. All earlier reports, scripts, book files and indices are read-only. No worktree, Git mutation, commit, push or render is used. The current checkout has advanced since the preceding report; this pass records its own initial snapshot rather than attributing that advance to referee work.

## Result

**The first genuinely new configuration is `BB barBbarB`. Its tested groupings agree through mass order two under the existing full-jet local forest and common minimal meromorphic continuation, at generic physical coupling in `2<t<3`, away from one new possible resonance.** The three-current configuration `B barBbarB` is the reflected previous one-barred test.

There is an essential, explicit correction to a tempting implementation: the one-barred average of the two extreme words cannot simply be applied to the two-barred word. It overcounts a nonzero double contact by `Delta12/4`. Balancing disjoint mixed faces independently gives the correct coefficient, and the complete Gaussian calculation proves that either order of doing so agrees. This coefficient is forced by the existing pair prescription; it introduces no new free normalization.

At second order the two-source graph has a genuine cycle. We do not claim that it cancels or that the unrenormalized physical-height limits exist. The grouping result uses the continuation prescription already declared in the earlier work. Its new unlike-source pole candidate is

`t_*=(7+sqrt(17))/4=2.780776406404415...`.

The complete residue there is unresolved. A nonzero bulk contribution has been computed, but it is not the residue of the full physical-current forest. This is neither a demonstrated associativity anomaly nor a proof of regularity at that coupling.

The result concerns alternative groupings of the same labelled four-current word, with the action inserted once. It does not establish a two-barred Serre extension, arbitrary null-relation descent, repeated application of the interaction map as a binary multiplication, or an unrestricted interacting algebra.

## 1. Conventions and the object compared

Use covariance `X^i(z)X^j(w)~−delta_ij log(z−w)`, bilinear scalar products, and

`a0=(ib/sqrt(2),b/sqrt(2),i)`, `a2=(ib/sqrt(2),−b/sqrt(2),−i)`,

`t=b²`, `T=a0+a2`, `d=−T/t`, `rho=2/t`, `beta=t−1`.

Then `ai²=−1`, `a0.a2=1−t`, `d.ai=1`, `d²=−rho`. Define

`B=:C exp(d.X):`, `C=a0.partial X`, `A=:exp(−d.X):`,

`Yi=:exp(ai.(X+Xbar)):`, `Ki=:exp((d+ai).X+ai.Xbar):`.

Both currents are even. Original constituent labels and the full cocycle

`epsilon(n,m)=(-1)^((nL+nR)(m0+m2)) exp[−i pi(nL mR−nR mL)/t]`

remain attached to every fused cluster. Separate regulators gL and gR belong to the left and right like-current pairs. Their final restriction is `gL=gR=rho`. The unlike source-pair factor is `|w0−w2|^(2h−2)`, with physical h=t. Equal source pairs have factor `|w0−w2|²`.

The input B is `B−t lambda2 K2`. Its rotational-to-spatial unary repairs are `+lambda0 K0+beta lambda2 K2`. Thus in the spatial representation of the complete physical input, the combined original-slot replacements have coefficients **+1 for K0 and −1 for K2**. A test of only the bare action graph omits these terms. Equal-source second-order coefficients have the action factor 1/2. The source-pair P counterterm, L improvement and all their Taylor descendants occur once, through the same radial forest, rather than as duplicate independent insertions.

For clarity if the even-root channel is considered, the current tail is `K_A=:exp(−k.X−T.Xbar):`, `k=d+T`; its mixed endpoint is the different field `:exp(−k.(X+Xbar)):`. The dangerous new two-source test here is the all-B word. The V1 contribution at mass order two has only one source and falls under the single-shared-vertex argument below.

Fix the internal order of B1,B2 and of barB1,barB2. There are six order-preserving mixed shuffles. Compare `(BB)(barBbarB)` with the phase-aligned `(BbarB)(BbarB)` grouping, retaining constituent order, source ancestry and the actual crossing contacts. These are not products of two separately dressed binary blocks.

## 2. Complete kernels and the first real cycle

For one source Yi put `q0=1`, `q2=beta`, `Z=z1−z2`, `V=ubar1−ubar2`, and

`F1=C1−1/Z+qi/(z1−w)`, `F2=C2+1/Z+qi/(z2−w)`,

`H2=F1 F2+1/Z²`.

Use the barred counterpart Hbar2. The full Wick kernel, apart from its fixed whole-word cocycle, is

`Z^gL V^gR H2 Hbar2 N /[(z1−w)(z2−w)(ubar1−wbar)(ubar2−wbar)]`.

N is the common normal exponential. The polynomial contains both derivative-field pairings and their product. These are necessary for the source and test jets, even if a scalar specialization simplifies.

For two sources replace each source term in Fj by `sum_i qi/(zj−wi)`. The scalar base has all eight source/current denominators, the same current-pair powers, and the source-pair factor specified above. This formula supplies the complete two-source kernel before extension.

With one source, the holomorphic and antiholomorphic directed acyclic graphs have only that vertex in common. A cycle visiting both would have to revisit the shared vertex. Thus their union remains acyclic and has the uniform boundary-value construction used previously.

With two sources the graphs share two vertices. For equal sources the exact identity

`det[1/(zi−wj)]=−(z1−z2)(w0−w2)/product_(i,j)(zi−wj)`

rewrites the scalar base as

`Z^(gL−1) V^(gR−1) det H det Hbar`.

One term admits the cycle `w0→B1→B2→w2→barB2→barB1→w0`. Assigning equal positive Fourier edge parameters to its six edges gives zero incidence image. The old proper-cone DAG bound therefore cannot be reused for this term.

Nor does the determinant remove all derivative singularities. With `a=qg−1` and `O_j=C_j+qD_j+a/(zj−zk)−q partial_zj`, where `D_j=d.partial X(zj)`, the exact identity is

`H2 K N=[O1 O2+(1−q)/Z²] K N`.

The residual term is real and must be kept. In divergence form the scalar coefficient of `K/Z²` is

`c_*=1−q−a²+qa`.

On `q=t−1=2/g−1` it equals `−(1−g)(2−g)`, which is nonzero on the physical patch. The coordinator derives these identities independently from the full H2 polynomial. They reject the claim that the full cycle is merely harmless external differentiation of its scalar determinant.

## 3. The full contact diamond

The old mixed first-jet functional, with the source coordinate differentiated before restriction, is

`L_i[psi]=qi²(partial_wi+partial_wbar_i)psi−qi(Ceff_i+Cbar_eff_i)psi`.

The raw ordered moment has prefactor `−i pi sigma nu lambda_i`, `nu=exp(−i pi/t)`; the balanced finite addition is its negative. The outer-slot-test first derivative cancels. Every source/test coefficient and spectator Wick contraction remains inside psi or the regular-current coefficient.

For two disjoint faces let x,y denote their source tangent directions. Introduce formal variables u1,v1,u2,v2 for the four regular derivative fields. The part of the complete Gaussian generating coefficient relevant to the two moments is

`G=psi[(1+C1 u1+Cbar1 v1)(1+C2 u2+Cbar2 v2)+H_L u1u2+H_R v1v2]`.

All coefficients may depend on x and y. H_L,H_R are the two cross-Wick pairings. Define

`M1=q1² partial_x−q1(partial_u1+partial_v1)`,

`M2=q2² partial_y−q2(partial_u2+partial_v2)`.

Apply them to G and only then set the formal variables to zero. Either order gives

`q1²q2² partial_x partial_y psi`

`−q1²q2 partial_x[(C2+Cbar2)psi]−q1q2² partial_y[(C1+Cbar1)psi]`

`+q1q2[(C1+Cbar1)(C2+Cbar2)+H_L+H_R]psi`.

This proves the disjoint collision diamond including the derivatives of the other source's current coefficient, both derivative-field crossmatchings, and the joint second source/test jet. The coordinator and both workers independently obtain it. Treating Ceff as a constant would miss terms and would not establish the diamond.

For a concrete nonzero joint moment choose `psi=xy` near the separated face centers, multiplied by smooth compact cutoffs equal to one there. Its value and first jets vanish at the joint restriction, while `partial_x partial_y psi=1`. The complete double moment is then `q1²q2² N`. At node 2 this is `beta^4 N`, times the nonzero ordered prefactors. No frozen-spectator approximation is used.

### The two-extreme shortcut is wrong

Let R00,R10,R01,R11 be the phase-aligned orientations of two disjoint mixed pairs, and

`Delta12=R11−R10−R01+R00`.

The difference between one global two-extreme average and independent local balancing is exactly

**`(R00+R11)/2−(R00+R10+R01+R11)/4=Delta12/4`.**

The joint test above makes Delta12 nonzero. This error occurs on separated mixed faces, so a new contact supported only on the total diagonal cannot repair it. The full forest instead uses the product of the two local half factors, including its quarter double-contact coefficient. Applying those local operations in either order agrees by the generating calculation.

This is not a uniform average of all six four-current shuffles. It is the compatible assignment of local proper-face moments. Overlapping faces sharing a source or current are union faces; they are never treated as products of the two disjoint projectors.

For physical equal-node endpoint jets, the double-contact base is `U(x)U(y)`, `U=:exp((d+a2).(X+Xbar)):`. Its inherited scalar factor is `|x−y|^(gL+gR−2)`. On the physical diagonal it is locally integrable because `2rho−1>0`. The complete differentiated product provides a well-defined double-contact extension, including its contractions and arbitrary tests. It is one part of the full kernel, not a replacement for the source loop.

## 4. Comparing grouping and collision orders

Write `E=gL binomial(bL,2)+gR binomial(bR,2)` and `n=bL+bR`, with each b at most two. Direct counting from the primary vectors gives, for a current-containing cluster,

`omega_one=E−n+1`,

`omega_equal=E−2n+5`,

`omega_unlike=E−2n+2h+1`.

Free source-only pairs instead have degrees 4 and 2h. The fixed unary faces have degree zero; a mixed pair/source has degree −1 and retains its prescribed full first jet.

Choose a common domain with `Re gL,Re gR>4` and sufficiently large Re h. The new union faces have positive margin. The two groupings agree away from collisions. Their complete proper-face assignments also agree: unary rotational residues, the full mixed jet, its reflected counterpart, the disjoint Gaussian diamond, and the full source-pair Taylor forest. Original derivative slots can be consumed only once. The P and L representatives and descendants are counted once.

After these proper faces are subtracted, a mixed inner first-jet remainder has degree +1. A nested sector with inner scale epsilon and next radius R is bounded by a full-test seminorm times

`epsilon integral_epsilon^1 R^(Re omega−2)dR`.

Comparable-scale sectors give `O(epsilon^Re omega)`. Both vanish in the common domain. Differentiating products of masks includes derivatives on different masks; the inverse radius they cost is offset by the corresponding less differentiated primitive. Thus the intersections obey the same bound. For unlike sources h is chosen sufficiently large to supply the finite coefficient differentiability needed before continuation; no physical-h smoothness is assumed.

This gives an actual equality of the differently grouped and nested kernels in the initial domain. Sector Taylor continuation extends that identity jointly in gL,gR,h. It is the inherited global minimal continuation, with no additional even primitive inserted in the initial domain. Defining unrelated finite parts after restricting the physical parameters would be a different scheme.

The comparison includes the physical finite tails. If N0,N2 denote combined original-slot replacements by K0 and −K2, the equal-source coefficient has the once-expanded content

`(1/2)A_ii(F)+A_i(N_iF)+(1/2)N_i²F`.

The unlike coefficient has

`A_02(F)+A_0(N2F)+A_2(N0F)+N0N2F`,

with the common radial P/L forest and full local balanced assignments. A_i denotes a spatial action insertion with its actual sign; A_ij is before the equal-source factorial. The action is not applied again to binary intermediate outputs. Each term inherits the same constituent data and comparison, so the finite sum does also.

## 5. Physical restriction and the new candidate resonance

In this finite all-B family, equal-source exponents depend on gL,gR through integer-affine combinations. The higher-face candidate divisors include `gL−2+j`, `gR−2+j`, `gL−1+j`, `gR−1+j`, and `gL+gR−3+j`, with Taylor order j a nonnegative integer. None meets `gL=gR=2/t`, `2<t<3`. Thus the equal-source continued grouping identity has a regular physical restriction throughout this open patch.

For unlike sources all proper two-source degrees are positive there. The new full union has

`omega_full=gL+gR+2h−7`.

Its physical value `4/t+2t−7` vanishes at t_*. In the open patch it exceeds −1, so only the zero-jet candidate can occur. At t_* the proper two-source faces remain regular. A surviving residue would consequently be a primitive quartic scalar contact with field

`:exp[(2d+T).(X+Xbar)]:`

and three real delta functions identifying the four current positions, multiplied by lambda0 lambda2 and the whole-word cocycle.

The construction worker makes this candidate concrete without claiming it survives. Integrating the complete scalar H2 over both chiral current coordinates, for opposite source heights, gives the bulk coefficient

`J=4pi² beta(g−1)(g+beta−3) r^(g−4)`

before the source-pair multiplier. At `g=2/t_*`, `beta=t_*−1`, the polynomial coefficient is exactly 1/4. The two opposite-height bulk sectors therefore have a nonzero radial pole coefficient `4pi⁴` before action normalization and cocycle. This computation keeps the derivative-pair term; it is not a single current-chamber residue.

Nevertheless it excludes source-height-supported unary repairs, balanced proper faces and physical input-tail terms. Those can affect the full residue. In particular the physical node-2 substitution is −1, not beta or −t separately. Positivity of an open-height bulk contribution cannot rule out cancellation by the supported terms.

The coordinator independently verifies that the constant-source Y0 mixed normal-field contact cancels, including a remaining unlike-source factor `(h−t)/(x−w2)`. This helps organize the forest but does not remove its unary and input-tail terms. We therefore leave the **complete** residue on `gL+gR+2h−7=0` uncomputed.

Away from this possible divisor the common meromorphic family has the stated regular physical restriction. At the divisor a cancellation would extend that result; a surviving residue would require a stated finite-part normalization at that exceptional coupling. Neither alternative would by itself create an associator: the two groupings are already the same meromorphic family, including their possible principal parts.

## 6. Verdicts and evidence

| Item | Verdict | Evidence |
|---|---|---|
| The two-barred triple is a new obstruction test | Wrong | Chirality reflection identifies it with the previous triple |
| The first new 2+2 one-source kernel is acyclic | Correct | Its two chiral DAGs share one vertex |
| The two-source cycle disappears from the full H2 kernel | Wrong | Exact residual differential term and positive circulation |
| Disjoint full-jet collision orders agree | Correct | Complete Gaussian generating coefficient with source-dependent contractions |
| A global two-extreme average implements two-bar balancing | Wrong | Explicit nonzero `Delta12/4` error on separated mixed faces |
| The tested 2+2 groupings agree under the old global continuation | Correct at the stated generic physical points | Proper-face comparison, union bounds, common continuation and pole inventory |
| The full physical residue at t_* is nonzero or cancels | Unclear | Bulk coefficient known; complete physical forest sum still missing |
| Raw physical-height convergence or unrestricted associativity | Unclear / not established | The cycle defeats the previous DAG bound; these broader claims were not proved |

The frozen independent reports are `briefs/referee_sl21_two_barred_construction_report.md` and `briefs/referee_sl21_two_barred_audit_report.md`. Their completed evidence has **30** and **107 PASS** records, respectively. The coordinator's fresh scripts under `referee_scratch/sl21_two_barred_20261024/coordinator/` give:

- `check_two_barred_v2.py`: **56 PASS**, primary-vector degrees, full scalar H2 identity, determinant sign, endpoint powers, loop scaling and resonance candidate.
- `check_contact_diamond_v2.py`: **7 PASS**, source-dependent full Gaussian joint jets, both crossmatchings and the node-0 remainder.
- `check_balancing.py`: **5 PASS**, the forced quarter double-contact coefficient and a nonzero counterexample to the shortcut.

Two initial coordinator assertions compared structurally different forms of the same expression. The corrected versions compare algebraic differences; the initial scripts and outputs remain. Workers likewise retain their interrupted calculations, corrected checking-code errors and explicit failed shortcut tests. PASS counts supplement the written distributional argument; they are not proof of a broader algebra.

## Integration and next mathematical input

No book integration or label upgrade is performed. The continuation verifies the smallest new generic grouping test and replaces the proposed global two-extreme shortcut by the forced local double-contact forest. This refinement leaves all earlier reports unchanged; the one-bar rule was only asserted on its stated domain.

The most localized remaining calculation is the full unlike-source quartic residue, including bulk, physical input tails, rotational unary repairs, balanced mixed faces and radial P/L descendants. A separate broader problem is a Serre quartet with two barred generators. Neither is settled by this four-current test. No literature or novelty claim is made.

## Preservation and concurrent work

The initial snapshot contains 4,090 pre-existing files. At final verification, **4,082 are byte-identical**, including all **2,380 snapshotted prior referee reports and scratch files**. Eight unrelated tracked lecture files changed during this pass: `code/lectures/README.md` and seven files under `lecture-notes/defects-gemini/`. This referee team wrote none of them and did not restore or otherwise alter those concurrent edits. HEAD, branch and index still match the initial state.

The unqualified whole-checkout preservation check therefore reports a failure, retained honestly in `referee_scratch/sl21_two_barred_20261024/preservation.out` and `preservation.json`. The exact changed-file list and successful research-input comparison are in `scoped_preservation.json`. Our write set is confined to the three new two-barred reports and the new scratch directory. `output_manifest.json` records their final hashes, including the frozen independent evidence.
