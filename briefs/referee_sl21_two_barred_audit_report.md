# Two barred currents: independent adversarial audit

Date: 2026-10-24. Model: Codex GPT-6. Reviewed commit: `483b0b9b9cc222e67cae61a549e7145151e905e4`.

This report concerns the smallest new two-barred test, the labelled word with two B and two barred B constituents, through mass order two. It compares groupings of this same word under the already selected full-first-jet balanced forest and global minimal meromorphic continuation. It does not define multiplication by repeatedly applying the interaction map, and does not assert a complete two-barred algebra. Earlier reports, scripts, book files and indices were left unchanged. No Git mutation, worktree or render was performed.

## Verdicts

1. **Correct:** the one-source two-barred triple is the reflected previous case. Even the one-source 2+2 family has no new directed cycle: the holomorphic and antiholomorphic graphs share only one source vertex.
2. **Wrong:** the equal two-source 2+2 family is covered by the previous one-antiholomorphic-leaf proof. There is a genuine positive circulation after the exact determinant reduction. Its presence is not itself an associator or an obstruction to analytic continuation.
3. **Correct:** the complete disjoint contact diamond commutes. This requires the full Gaussian coefficient, including source-dependent spectator contractions, both chiral Wick crossmatchings, and mixed second source/test jets. A constant-coefficient surrogate alone does not establish it.
4. **Wrong:** the one-barred average of the two extreme words extends unchanged to two barred currents. It differs from independent balancing of two disjoint mixed faces by one quarter of the double crossing. An explicit smooth-test jet makes that difference nonzero. The required correction is fixed by the old local prescription; it is not a new adjustable constant.
5. **Correct for the bounded continued construction:** using that complete local forest, the two groupings agree in a common convergent domain and hence as meromorphic families. In the equal-source sector their physical restriction is regular throughout the open patch `2<t<3`. For the unlike-source sector the same conclusion holds away from the single new candidate resonance

   `t_*=(7+sqrt(17))/4 = 2.780776406404415...`.

   This statement includes the specified input tails and unary repairs once, and the source-pair/P/L forest once. It does not require existence of the unrenormalized simultaneous physical-height limit.
6. **Unclear:** cancellation of the full unlike-source residue at `t_*`. The bulk subkernel has a possible full-diagonal zero-jet pole. A nonzero bulk coefficient does not determine the residue after input tails, unary repairs, balanced proper faces and source-pair counterterms are included. No nonzero associator or unavoidable new free constant is established here.

There is therefore a constructive bounded answer at generic coupling, and a sharply localized exceptional-coupling question. Neither a raw loop divergence nor the new possible resonance overturns the generic grouping result.

## 1. Definitions reconstructed afresh

Use covariance `X(z) X(w) ~ -log(z-w)`, bilinear scalar products, and

`a0=(ib/sqrt(2),b/sqrt(2),i)`,
`a2=(ib/sqrt(2),-b/sqrt(2),-i)`,
`t=b²`, `d=-(a0+a2)/t`, `rho=2/t`, `beta=t-1`.

Fresh computation gives `a0²=a2²=-1`, `a0.a2=1-t`, `d.ai=1`, and `d²=-rho`. Thus `B=:C exp(d.X):`, `C=a0.partial X`, and `Yi=:exp(ai.(X+Xbar)):` have source/current scalar exponent -1. The derivative/source coefficients are `q0=1`, `q2=beta`.

Introduce independent like-current regulators `gL,gR` for the two chiral pairs, retaining the constituent labels through every fusion. Their physical restriction is `gL=gR=rho`. The action has one insertion of each source per occurrence; equal sources have expansion factor 1/2. The full word cocycle is

`epsilon(n,m)=(-1)^((nL+nR)(m0+m2)) exp[-i*pi*(nL*mR-nR*mL)/t]`.

Both exponents are bilinear. The fresh script verifies the two cocycle identities for three arbitrary constituent labels. Consequently disjoint phase-aligned crossings use the same whole-word phase in either order. Replacing it by independently chosen binary phases would not be this test.

The specified B input tail is `-t lambda2 K2`, with

`Ki=:exp[(d+ai).X+ai.Xbar]:`.

The rotational unary repairs add `+lambda0 K0` and `+beta lambda2 K2` on original B slots, and their barred analogues. Therefore their combined coefficients are `+1` at node 0 and `-1` at node 2. They act on original slots once; a tail is not treated as another original C insertion. The even current tail, if used, is `K_A=:exp[-k.X-T.Xbar]:`, with `T=a0+a2` and `k=d+T`; it is distinct from the AA endpoint `:exp[-k.(X+Xbar)]:`.

## 2. The new scalar cycle and the complete H2 kernel

For two equal sources at `w0,w2`, write `r=z1-z2`, `rbar=ubar1-ubar2`. Apart from the fixed cocycle and the common normal exponential, the vertex scalar is

`V=r^gL rbar^gR |w0-w2|² / [product_(i,j)(zi-wj) product_(a,j)(ubar_a-wbar_j)]`.

With `H_ij=1/(zi-wj)`, the elementary identity is

`det H = -(z1-z2)(w0-w2)/product_(i,j)(zi-wj)`.

The two signs cancel between chiralities, so

`V=r^(gL-1) rbar^(gR-1) det H det Hbar`.

This removes the explicit source-pair edge. It does not remove all source-to-source paths. One determinant matching admits the directed cycle

`w0 -> B1 -> B2 -> w2 -> barB2 -> barB1 -> w0`.

The incidence matrix annihilates the vector with all six edge parameters equal to one. This is a nonzero nonnegative circulation, unlike the one-barred directed acyclic graphs. It can occur when the two sources lie on opposite sides of all current heights, with compatible chiral orders. The fresh incidence calculation is exact.

For one source, the two chiral graphs share only that source. A simple directed cycle using both graphs would have to visit the shared vertex twice. Neither graph has an internal cycle. Thus that case remains acyclic; merely having two barred currents is not enough to create the new difficulty.

The derivative fields must also be kept. For either equal source node, set `q=qi` and

`F1=C1-1/r+q[1/(z1-w0)+1/(z1-w2)]`,
`F2=C2+1/r+q[1/(z2-w0)+1/(z2-w2)]`.

The exact chiral polynomial is `H2=F1 F2+1/r²`. With `a=qg-1`, and after compensating the normal derivative by `D_i=d.partial X(zi)`, define

`O_i=C_i+qD_i+a/(zi-zj)-q partial_zi`.

Then

`W_H2=[O1 O2+(1-q)/r²] V`.

Here V includes the normal exponential when D is retained. Equivalently the script takes the scalar vertex factor and omits both normal derivative and D. Direct differentiation verifies the complete identity. In particular the residual `(1-q)/r²` cannot be discarded.

For the most singular scalar part, reorganizing into derivatives of complete coefficients gives

`q² partial_1 partial_2 V + qa[partial_1(V/r)-partial_2(V/r)] + c_* V/r²`,

`c_*=1-q-a²+qa=(1-a)(1+a-q)`.

At physical `q=t-1=2/g-1`, this is `c_*=-(1-g)(2-g)`. The coefficient does not vanish in the physical patch. This is a useful warning against declaring the entire kernel a harmless external derivative. It is **not** the coefficient of the complete interacting forest: the latter also has the input-tail and contact terms listed in section 5.

The bare full-union relative degree is `gL+gR-3`. It is negative at physical equal couplings. Thus neither the determinant identity nor a positivity assertion proves a raw physical height limit. The result below instead uses the prescribed common meromorphic extension.

## 3. Full Gaussian disjoint-contact diamond

For one ordered mixed face, the full first-jet moment, apart from its fixed nonzero ordered prefactor, is

`L_i[psi]=qi² D_i psi-qi(C_i+Cbar_i)psi`,

where `D_i=partial_wi+partial_wbar_i` acts on the complete source/test/normal coefficient. In particular it acts on contractions with every spectator. Its finite balanced addition is the negative of the ordered moment. The outer-current test derivative cancels in this formula.

For two **disjoint** faces, choose local source tangent coordinates x,y, and introduce formal variables `u1,v1,u2,v2` for the four derivative fields. The complete Gaussian generating coefficient through the needed degree is

`G=psi [1+C1 u1+Cbar1 v1+C2 u2+Cbar2 v2`
` +(C1 C2+H12)u1u2+C1 Cbar2 u1v2`
` +Cbar1 C2 v1u2+(Cbar1 Cbar2+Hbar12)v1v2]`.

All `C`, `H`, and psi are allowed to depend on both x and y. The H terms are the two derivative-field Wick crossmatchings. No mixed-chirality derivative-field contraction is inserted.

Apply

`M1=q1² partial_x-q1(partial_u1+partial_v1)`,
`M2=q2² partial_y-q2(partial_u2+partial_v2)`,

and then put all formal variables to zero. These are constant-coefficient differential operators on the **complete generating coefficient**, rather than operators with frozen `C_eff`. They commute. Explicitly either order gives

`q1²q2² partial_x partial_y psi`
`-q1²q2 partial_x[(C2+Cbar2)psi]`
`-q1q2² partial_y[(C1+Cbar1)psi]`
`+q1q2[(C1+Cbar1)(C2+Cbar2)+H12+Hbar12]psi`.

This is an identity for arbitrary source-dependent coefficients and unequal q's. It retains `partial_x C2`, `partial_y C1`, both Wick crossmatchings, and the mixed second jet of the joint test. The independent symbolic check verifies the identity and each of these coefficients. Formal derivative variables are merely a compact encoding of the actual finite Wick polynomial; they do not introduce an extra field or prescription.

On separated supports the mixed Taylor operations act in independent local coordinates. Applying this calculation to the entire smooth coefficient therefore proves compatibility of the two disjoint proper-face assignments. The normal fields and the other source remain inside the coefficient until both operations have been performed.

It is wrong to use this as a product of two contacts sharing a source or current. Shared-source faces are overlapping, not disjoint. A single original C can contract only once. Such configurations belong to the union forest discussed below, never to a product of coincident source-height delta functions.

### A concrete failure of the two-extreme shortcut

Let `R00,R10,R01,R11` be the phase-aligned values for the four choices of orientation of two separated mixed pairs. Set

`Delta12=R11-R10-R01+R00`.

Independent local balancing averages all four choices. Averaging only the two global extremes gives

`(R00+R11)/2 -(R00+R10+R01+R11)/4 = Delta12/4`.

This difference is genuinely nonzero for admissible tests. Choose the local smooth source/test coefficient to be `psi=x y` near the two separated face centers, multiplied by compact cutoffs equal to one there. At the centers its value and its first derivatives vanish, while its mixed derivative is one. The complete double moment above becomes exactly `q1²q2²` times the normal coefficient. Every spectator-contraction term vanishes because it multiplies psi or a first derivative. At node 2 the coefficient is `beta^4`, times the two nonzero ordered jump prefactors. This avoids any frozen-spectator approximation.

With flat switching and physical normal jets the corresponding field contribution is the fully differentiated `U(x)U(y)` kernel, including its contractions. The scalar U-U factor is `(x-y+i0)^(g-1)(x-y-i0)^(g-1)=|x-y|^(2g-2)` on the common physical branch; it is locally integrable for `g>1/2`. Differentiation supplies a particular complete extension of the double endpoint term. This observation does not identify that term alone with the full loop kernel.

Thus the two-extreme shortcut fails already away from the total four-current diagonal. Adding an arbitrary primitive contact on the total diagonal cannot repair it. The prescribed disjoint forest coefficient does repair it, with no new free parameter.

## 4. Union bounds and physical pole inventory

For a current-containing cluster with `bL,bR` B constituents, each at most two, write

`E=gL binomial(bL,2)+gR binomial(bR,2)`, `n=bL+bR`.

The one-source relative degree is `E-n+1`. The equal two-source relative degree is

`omega_equal=E-2n+5`.

The unlike two-source degree is

`omega_unlike=E-2n+2h+1`.

These follow directly from the scalar contractions, the maximum internal derivative degree -n, and the relative measure dimension `n+3` for two sources and a real current anchor. The source-only pair is different: its free translation is two dimensional. Its degrees are 4 for equal sources and 2h for unlike sources. The first script also checks the formal anchored n=0 expression; that expression is not used for the free source-only face. The dedicated pole script uses the correct two-dimensional source-only degrees.

The potentially new equal two-source unions have degrees `gL-1`, `gR-1`, and `gL+gR-3`; the `(1,1)` mixed union has degree +1. The unlike counterparts replace +5 by `2h+1`. For one source, the new triples and four-current union have degrees `gL-2`, `gR-2`, and `gL+gR-3`. The fixed unary and mixed-pair source faces have degree 0 and -1 and already have their prescribed rotational and whole-first-jet extensions.

Take `Re gL,Re gR>4` and Re h sufficiently large. All the new unions have a margin greater than one. The difference between the two source/current mask groupings has identical separated kernels and identical complete proper-face moments: unary Wick residues, the full mixed moment of section 3 and its reflected case, disjoint products as just computed, and the full source-pair Taylor forest. No same-source product of two overlapping contacts is inserted.

After a mixed first-jet subtraction the inner remainder has degree +1. A nested sector with inner radius epsilon and next scale R is bounded by a finite full-test seminorm times

`epsilon integral_epsilon^1 R^(Re omega-2) dR`.

It tends to zero uniformly in the stated domain. Comparable-scale sectors are bounded by `epsilon^Re omega`. Derivatives of moving masks are included by the product rule; each lost power of radius is paired with the corresponding less differentiated primitive. They obey the same degree bound. Source-pair Taylor subtractions retain normal, test and switching jets before evaluation. For unlike sources, choosing h sufficiently large makes the remaining pair coefficient as differentiable as this finite collection of operations requires. No physical-h smoothness is assumed.

This accounts for overlapping crossings: their support is a source/current union, not an independent product of two proper faces. Once the proper moments have been fixed, the positive margin makes the remaining difference vanish in the initial domain. The full Gaussian calculation in section 3 is essential to this argument; positivity without that proper-face comparison would not establish equality.

Sector Taylor integration then gives a meromorphic family. Equal-source exponents depend only on `gL,gR`, with integer shifts. Candidate new divisors consequently have forms `gL+integer`, `gR+integer`, or `gL+gR+integer`. A finite cluster enumeration, including bare-current faces, gives no such pole on `gL=gR=2/t`, `2<t<3`. The fixed unary and mixed-pair faces are already prescribed distributions, not additional unfixed poles in g.

For unlike sources, the proper two-source degrees are positive on that physical patch. The sole new full-union degree is

`omega_full=gL+gR+2h-7`,

whose physical restriction is `4/t+2t-7`. It increases through zero at `t_*`. Its minimum in the open patch is greater than -1, so there is only the zero-jet candidate, not a derivative-jet series of new physical poles. At `t_*` all proper two-source degrees remain positive. The new residue, if nonzero after the whole forest is summed, is therefore a primitive full-diagonal zero-jet term.

Its possible local field has the total charge

`:exp[(2d+T).(X+Xbar)]:`,

multiplied by `lambda0 lambda2` and three real relative-current delta functions, with the word cocycle and a coefficient still to be computed. This is a candidate support and field content, not a claim that the coefficient survives.

## 5. What the comparison includes and what it does not

For clarity, let `A_i` denote the spatial one-source action operator, with its actual action sign; `A_ij` the two-source operator before the equal-source factorial; and `N0,N2` the combined original-slot tail/unary substitutions with coefficients +1,-1. The complete physical-input equal-source coefficient has the schematic but unambiguous forest content

`R_bal[ (1/2)A_ii(F)+A_i(N_i F)+(1/2)N_i²F ]`.

The unlike coefficient has

`R_bal[ A_02(F)+A_0(N2 F)+A_2(N0 F)+N0N2F + source-pair P/L forest ]`.

`R_bal` here is the once-applied local forest, with the complete whole-jet balance and the prescribed radial Taylor subtractions. It is **not** a two-extreme-word average. The source-pair P/L terms include their inherited labels and descendants once. An independent action series is not attached to each binary block, and the map is not applied a second time to its own tail output.

Each term has the same grouping comparison in the initial domain, using its inherited complete coefficient. The finite sum therefore has that comparison too, regardless of possible cancellations of individual raw divergences. The resulting meromorphic identity restricts to the regular physical points described above. This is the precise generic bounded sense in which the old full-jet rules and common continuation suffice.

There is no proof here that the raw physical source-height integrals, without the common continuation, converge uniformly. The cycle calculation explicitly prevents importing the old acyclic-graph proof. Nor is a nonzero scalar bulk residue a calculation of the full displayed sum. In particular the net node-2 substitution is -1, not beta or -t alone; omitting either component changes the residue question.

An additional A constituent in a two-Yukawa-source face replaces two current/source factors of exponent -1 by exponent +1 and removes its derivative-current pole, improving that contribution to the relative degree by five. Changing a like current pair to an unlike one costs `g+rho`, less than two on the physical patch. Hence such replacements do not give a more dangerous full two-source union than all B in this four-current test. The mass-order-two V1 channel has only one source and remains in the single-shared-vertex case. This degree comparison identifies the dangerous minimal configuration; it is not a separate complete calculation of every mixed-species word.

## 6. Minimal remaining input and replacement wording

The remaining exceptional-coupling computation is the residue of the **complete** unlike-source coefficient in section 5 on

`gL+gR+2h-7=0`,

followed by the physical restriction at `t_*`. The residue must include bulk, input-tail, rotational unary, full balanced mixed, disjoint double, and radial source-pair/P/L contributions. A proof of cancellation would remove the exceptional point. A nonzero residue would require stating the chosen finite-part prescription and its local normalization there. It would still not, by itself, be a nonzero associator: the two groupings are the same meromorphic family.

Recommended replacement for an unqualified two-barred claim:

> For the labelled two-left/two-right B word through mass order two, the complete local first-jet forest, applied once and then minimally continued as a common family, is independent of the tested source/current groupings. Its physical restriction is regular for generic 2<t<3 away from the possible unlike-source resonance t=(7+sqrt(17))/4. Two barred constituents require the disjoint double-contact forest; the one-barred two-extreme average alone does not implement it. The full residue at the exceptional coupling and an unrenormalized physical-height limit remain unproved.

This replacement concerns the referee continuation under review, not an existing book statement. No book line is proposed for editing.

## 7. Fresh evidence and reproducibility

All scripts were written from the definitions above; none imports or reruns an older task's script.

- `referee_scratch/sl21_two_barred_20261024/audit/check_two_barred.py` and `.out`: **57 PASS, 0 FAIL**. Vectors, determinant, full H2 differential representation, scalar residual, elementary contact coefficients, face degrees and an explicit positive cycle. The initial unexpanded rational simplification was interrupted after seven PASS records; its partial output is retained. The bounded expanded calculation completed.
- `referee_scratch/sl21_two_barred_20261024/audit/check_full_jet_diamond.py` and `.out`: **10 PASS, 0 FAIL**. Arbitrary source-dependent Gaussian generating jets, both Wick crossmatchings, mixed second jet, constituent cocycle and the nonzero extreme-average correction. An initial coefficient-extraction assertion used an unexpanded expression; its two preceding full identities had passed. Expanding before coefficient extraction fixed the checking code. The partial output and explanation are retained.
- `referee_scratch/sl21_two_barred_20261024/audit/check_pole_inventory.py` and `.out`: **40 PASS, 0 FAIL**. Every current-containing source subset in the bounded all-B family, correct free source-pair degrees, bare-current faces, the unique candidate physical resonance, and regular proper two-source faces there.

There are **107 completed PASS records** across these three final outputs. Counts are diagnostics, not substitutes for the derivations and scope restrictions above. `checkpoint.md` records both progress and the corrected checking-code episode.

## Integration

No integration occurred. This audit verifies a bounded generic continued grouping result, rejects the two-extreme shortcut, and leaves the complete exceptional residue unresolved. It recommends no unrestricted algebra theorem, no global symmetry upgrade, and no change to the book or its indices.
