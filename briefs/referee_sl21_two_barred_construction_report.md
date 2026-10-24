# Two barred generators: construction and crossing-diamond referee

Date: 2026-10-24. Model: Codex GPT-6. Reviewed commit: `483b0b9b9cc222e67cae61a549e7145151e905e4`.

This is a new referee-only calculation. Previous reports and scripts remain immutable. No book file, index, existing report, Git state or rendered output is changed. Scratch work is confined to `referee_scratch/sl21_two_barred_20261024/construction/`. All displayed checks are independently reconstructed; no previous research script is imported or rerun.

## Verdicts

| Item | Verdict | Evidence and limit |
|---|---|---|
| `B barB barB` is the first genuinely new sector | wrong | It is the chirality reflection of the already tested `BB barB` sector. |
| `BB barBbarB` has a complete finite Wick kernel | correct | The two explicit H2 polynomials below include every derivative-field matching. |
| The one-source 2+2 configuration necessarily has a directed pinch | wrong | Its two chiral DAGs meet at only one source vertex; their union is acyclic. |
| Disjoint mixed-face crossing orders commute with the whole-jet prescription | correct | A complete Gaussian generating coefficient gives the same joint moment in both orders, including crossmatchings and the mixed second jet. |
| A single average of the two extreme four-current words extends the one-bar construction | wrong | It overcounts the disjoint double-crossing term by one quarter. A local forest average has the correct coefficient. |
| Grouping and collision order agree as meromorphically continued kernels | correct in the bounded 2+2 sector with the inherited global continuation normalization | Proper-face matching and the crossing diamond hold in a common convergent domain; remaining unions vanish there. Continue that identity, rather than assuming an unrenormalized physical height limit. |
| Every physical value `2<t<3` is automatically regular | unclear | A genuinely new unlike-source full-union candidate pole occurs at `t*=(7+sqrt(17))/4`. Its bulk coefficient is explicitly nonzero; the complete forest residue, including physical input tails, is not computed here. |
| The bulk nonzero coefficient proves a physical associativity anomaly | wrong | Unary, balanced and finite-current-tail contributions can alter the full residue. The bulk calculation does not include them. |

Thus there is a concrete extension and a positive crossing-diamond test away from the new resonant divisor. There is also a precise remaining scalar residue test. Neither an unrestricted associative interacting algebra nor an anomaly at the candidate resonance is established.

## 1. Definitions and minimal configuration

Use `t=b²`, `rho=2/t`, `beta=t−1`, `d=−(a0+a2)/t`, `C=a0.partial X`, and

`B=:C exp(d.X):`, `A=:exp(−d.X):`, `Yi=:exp(ai.(X+Xbar)):`.

The bilinear products are `ai²=−1`, `a0.a2=1−t`, `d.ai=1`, `d²=−rho`. The current tails are

`Ki=:exp((d+ai).X+ai.Xbar):`,

`K_A=:exp(−k.X−T.Xbar):`, where `T=a0+a2`, `k=d+T`.

The even mixed endpoint is instead `F_A^end=:exp(−k.(X+Xbar)):`. K_A and this endpoint are different fields and cannot be substituted for one another.

The canonical input B current is `B−t lambda2 K2`. The rotational-to-spatial unary repair contributes `+lambda0 K0+beta lambda2 K2`. Thus, in the once-expanded spatial representation of the physical input, the net original-slot substitutions are `+lambda0 K0−lambda2 K2`. This is essential in the unresolved resonance test. It is not permissible to examine only the bare four-current action graph and call it the full current product.

Retain original constituent order and cocycle, and use separate auxiliary exponents gL,gR for the two like-current pairs until the final diagonal specialization. The physical values are `gL=gR=rho`. The unlike-source factor is `|w0−w2|^(2h−2)`, with physical h=t. Constituent labels and source ancestry remain attached to all composites.

Reflection exchanges the proposed triple `B barB barB` with `BB barB`. The first new configuration has two constituents in each chirality. Fix their internal orders as `B1<B2`, `barB1<barB2`. There are six order-preserving mixed shuffles. Their only independent crossing diamond exchanges two disjoint adjacent mixed pairs; this is the new second-order consistency check.

## 2. Complete kernels, including the derivative matchings

For one Yi source let `q0=1`, `q2=beta`, and write

`Z=z1−z2`, `V=ubar1−ubar2`, `pj=zj−w`, `qk=ubark−wbar`.

Define

`F1=C1−1/Z+qi/p1`, `F2=C2+1/Z+qi/p2`,

`H2=F1 F2+1/Z²`.

The barred polynomial has precisely the same form with barred fields and coordinates. The full one-source kernel, apart from the fixed whole-word cocycle, is

`Z^gL V^gR H2 Hbar2 N/(p1 p2 q1 q2)`.

N is the common normal exponential. The product H2 Hbar2 contains the empty derivative matching, the holomorphic pair, the antiholomorphic pair, and the product of those two pairings. No mixed holomorphic/antiholomorphic derivative contraction exists.

At zero regular derivative fields the one-source H2 is `qi(qi−1)/(p1 p2)`. This cancellation is useful, but it must not replace the full polynomial when arbitrary field and test jets are required.

For two sources define `pj0=zj−w0`, `pj2=zj−w2`. Replace the source term in Fj by

`1/pj0+beta/pj2`

for unlike sources, or `qi/pj0+qi/pj2` for equal sources. The full scalar base has all eight source/current denominators, the same current-pair powers, and respectively `|w0−w2|^(2h−2)` or `|w0−w2|²`. Multiply it by the two complete H2 polynomials and N. This gives the actual two-source Wick kernel before any extension.

For unlike sources its full holomorphic scalar polynomial simplifies to

`H2|C=0 = beta[(beta−1)/(p12 p22)+1/(p10 p22)+1/(p12 p20)]`.

The node-0 same-source term cancels; the node-2 same-source term generally does not. This identity is independently checked and used in section 6.

A useful differential representative for the one-source chiral factor is obtained from `K=Z^g/(p1p2)` and

`O_j=C_j+qi D_j+(qi g−1)/(zj−zk)−qi partial_zj`, `D_j=d.partial X(zj)`.

The exact identity is

`H2 K N = [O1 O2+(1−qi)/Z²] K N`.

The additional term `(1−qi)/Z²` is present. It prevents a claim that the full kernel is merely two harmless external derivatives of the scalar base. The independent check includes both normal-field first jets through N; it does not freeze them.

For equal sources the Cauchy determinant identity is

`1/(p10 p22)−1/(p12 p20) = −Z(w0−w2)/(p10 p12 p20 p22)`.

It is useful for constructing a common meromorphic extension and exposing a possible cycle. It does not cancel the residual term in the preceding differential identity.

## 3. The one-source comparison and the two-source loop

With one source, the holomorphic current/source graph and the antiholomorphic current/source graph have only one common vertex. Each graph is directed acyclic with inherited imaginary-height order. Their one-vertex union is again acyclic: a directed cycle visiting both pieces would have to pass through that shared vertex twice. The normalized Fourier-Schwinger edge family therefore has the same proper positive-cone pushforward as before. Its distributional seminorms are uniform as height gaps shrink. This gives a complete first-order 2+2 boundary family, not just separate angular chambers.

The source-containing degree for b holomorphic and c antiholomorphic B constituents is, with n=b+c,

`omega1 = gL binomial(b,2)+gR binomial(c,2)−n+1`.

The fixed unary degree-zero and mixed-pair degree-minus-one faces receive the old rotational and whole-first-jet prescriptions. The first new full union has degree `gL+gR−3`. It is positive in a large-g initial domain. After subtracting proper jets, its residual moving-mask terms vanish there. The DAG family then supplies a coherent continuation to the physical patch; the negative physical degree by itself is not an obstruction.

With two sources the chiral pieces share two source vertices, so the previous graph argument no longer excludes a directed cycle. For equal sources, the determinant decomposition exposes a scalar loop whose Fourier tail is of the form `integral dk k^(−gL−gR)`. The scalar tail is integrable when `Re(gL+gR)>1`, which includes the physical patch. However the full H2 polynomials contain the residual terms above and higher pole powers. The scalar observation alone does not prove a raw full-kernel height limit.

An independent complete-current contour calculation makes this distinction explicit. For two equal sources of coefficient q, the full scalar H2 gives

`J_equal=4pi² q(g−1)(qg−2) r^(g−3)`.

At q=beta and g=rho its coefficient is `4pi²(1−g)(2−g)`, which is nonzero in the physical patch. The remaining open-height bulk source integral has radial degree `2g−3` and diverges as its lower cutoff is removed. The computation includes the derivative-field pair in H2. Thus the milder scalar determinant loop bound is not a bound on the full raw bulk kernel. This is still a bare open-height calculation; it does not settle cancellation by the prescribed forest and physical input-tail terms.

There is still an inherited analytic-extension route: first choose gL,gR, and h if needed, sufficiently large that the full derivative kernel and all residual unions are controlled; perform the prescribed proper-face matching there; then continue the complete distribution jointly. This is the same global continuation normalization already used for the one-bar task. It must be distinguished from taking all unrenormalized physical cutoff limits directly. A finite part chosen afresh at the physical coupling would be an additional prescription.

## 4. Full mixed jets and the crossing diamond

For a single mixed face with arbitrary coefficient psi, the established moment is

`L_i[psi]=qi²(partial_w+partial_wbar)psi−qi(Ceff+Cbar_eff)psi`.

All quantities are restricted to the fused face only after differentiation. In particular Ceff includes its dependence on the other source and on the other derivative fields. A constant-coefficient toy commutator does not test the new diamond.

Here is a full generating-coefficient calculation. Introduce pair polarizations e1,e2. Let x,y denote the two source-jet directions and write, through degree one in each polarization,

`G=psi(x,y)[1+e1 A(x,y)+e2 B(x,y)+e1 e2(A B+H)]`.

A and B are the two complete regular derivative-current sums. H is the sum of the holomorphic and antiholomorphic cross-Wick pairings between the two pairs. Their dependence on the other source is retained. Define

`M1=q1² partial_x−q1 partial_e1`, `M2=q2² partial_y−q2 partial_e2`.

Extracting e1=e2=0 gives

`q1² q2² partial_x partial_y psi`

`−q1² q2 partial_x(B psi)−q1 q2² partial_y(A psi)`

`+q1 q2(A B+H)psi`.

The same expression results from M2 M1. This calculation retains the mixed second source/test jet, differentiation of the other regular-current coefficient, and both derivative-field Wick crossmatchings. It remains valid for q1 unequal to q2. It proves the disjoint-face diamond by direct extraction from the complete Wick tensor.

For disjoint faces, the two Taylor restrictions commute before the source variables are integrated. This identifies the two possible collision orders and the corresponding regroupings. Their full-intersection remainder is handled by the common initial-domain extension, rather than by treating the displayed local differential operators as a separate associative algebra.

### A concrete failure of the two-extreme shortcut

Let four phase-aligned ordered values on the diamond be

`R00=R`, `R10=R+Delta1`, `R01=R+Delta2`,

`R11=R+Delta1+Delta2+Delta12`.

Independent local balancing gives

`(R00+R10+R01+R11)/4 = R+(Delta1+Delta2)/2+Delta12/4`.

A global average of only the two extremes instead gives

`(R00+R11)/2 = R+(Delta1+Delta2)/2+Delta12/2`.

Its excess is `Delta12/4`. This term need not vanish. On two separated mixed faces choose a smooth coefficient whose local jet is the product of the two source coordinates relative to those faces. Its individual zero and first jets vanish at the joint restriction, but its mixed second jet is one. The displayed joint moment is then `q1²q2²`, which is nonzero. This is an actual arbitrary-test counterexample on disjoint faces, with no appeal to the unresolved total collision.

Thus the one-bar extreme-average rule cannot simply be applied once to the entire two-bar word. The existing local whole-jet forest prescription does have the correct independent half factors. Applying its two local averages in either order gives the same result. The quarter coefficient is forced by factorization on disjoint faces, not an added primitive full-diagonal normalization.

## 5. Collision forests and meromorphic grouping identity

The permitted forest types are concrete:

1. Isolated unary current/source faces retain their rotational Cauchy extension and complete spectator Wick coefficients.
2. A same-source mixed overlap receives the single whole-first-jet functional. Two deltas consuming the same source are not multiplied.
3. Disjoint source/current or mixed faces use the complete joint Gaussian moment of section 4, including crossmatchings. Original derivative slots may be consumed only once.
4. Source-pair radial Taylor subtraction differentiates the full normal field, tests and spectator contractions. Its P counterterm and L current improvement are counted once, including their descendants.
5. Residual union faces use the inherited joint continuation from the initial convergent domain, with no separately added even primitive term there.

The mask comparison retains all product-rule terms, in particular intersections of derivatives of different source/current masks. On each proper face the full coefficient and its required jets agree. A one-source mixed first-jet remainder has positive inner degree, giving the standard bound `epsilon integral_epsilon^1 R^(Re omega−2)dR`. Comparable-scale unions give `epsilon^Re omega`. The same reasoning applies after a source-pair Taylor forest, in the high-h domain where the remaining source coefficient has sufficiently many derivatives. Equal-source unions are controlled by sufficiently large gL,gR. This establishes equality of the alternative collision orders in a common initial domain.

For two equal sources,

`omega_equal = gL binomial(b,2)+gR binomial(c,2)−2n+5`.

For two unlike sources,

`omega_unlike = gL binomial(b,2)+gR binomial(c,2)+2h−2n+1`.

The formulas include every source/current derivative pole conservatively. The scalar normal-field sector reaches the full-union degree in question. Smooth positive-order field or test jets improve these degrees. Source-only faces use their usual separate radial count.

Consequently `(BB)(barBbarB)` and the phase-aligned `(BbarB)(BbarB)` give the same continued full kernel when intermediate clusters retain constituent data and the action is inserted once. The six order-preserving shuffles have the same local crossing data, and the only independent diamond is section 4. This is a finite family of meromorphic identities, not a free-algebra definition and not repeated application of the interaction map to each binary product.

At physical gL=gR=rho in `(2/3,1)`, one-source and equal-source higher-face pole candidates lie at integer-affine values such as `g−2+j=0`, `g−1+j=0`, and `2g−3+j=0`. None meets the open physical interval. Unlike three-current unions have `rho+2t−5>0`. The new full unlike union is different:

`omega_full = 2rho+2t−7 = 4/t+2t−7`.

It vanishes at `t*=(7+sqrt(17))/4`, a generic irrational value in `(2,3)`. The meromorphic grouping identity is therefore regular on the tested physical interval away from this possible divisor, assuming the proper-face normalizations already stated. The full resonance cannot be dismissed merely because all smaller faces were regular.

## 6. An exact nonzero bulk coefficient at the candidate resonance

This calculation deliberately stops short of calling the bulk result the full forest residue.

Take w0 in the upper half-plane and w2 in the lower half-plane, with `r=w0−w2`, and integrate both holomorphic current positions over the real line with their coherent `(x1−x2+i0)^g` boundary value. Put the regular derivative fields to zero only for this scalar residue calculation. From the full H2 identity in section 2, introduce

`f0(x)=1/[(x−w0)(x−w2)^2]`,

`f1(x)=1/[(x−w0)^2(x−w2)]`.

The holomorphic integrand is

`beta (x1−x2+i0)^g [(beta−1)f0(x1)f0(x2)+f1(x1)f0(x2)+f0(x1)f1(x2)]`.

Fourier-Schwinger transformation of the current-pair power reduces the two contour integrals to elementary residues. Their polynomial, apart from `4pi² exp(ikr)`, is

`[(ik)^2/r²+(beta−3)ik/r³+(3−beta)/r⁴]`.

Applying this polynomial to r^g gives the complete current integral

`J = 4pi² beta(g−1)(g+beta−3) r^(g−4)`

with its coherent branch. The opposite source-height sector gives the corresponding reflected expression. The antiholomorphic integral is its complex conjugate for the mirrored current order, so the bulk contribution is a positive square after removing the fixed whole-word cocycle.

An independent check uses `E=a0+d/rho`. Then `E.d=0`, `E²=(t−2)/2=:c`, and the integrated B current differs from `:E exp(d.X):` only by a total current derivative. Its scalar current-integral coefficient is

`c−c²(g−1)(g−4)`.

At g=2/t this equals `beta(g−1)(g+beta−3)` exactly. This independently confirms the full derivative-pair contribution in J.

Write `H=|Im w0|+|Im w2|`, `s=|Im w0|/H`, and a=Re r. The measure is `da H dH ds`. Polar coordinates in (a,H) produce the radial integral

`integral R^(2g+2h−8) dR`.

The two opposite-height sectors have angular factor four. At the candidate resonance,

`beta(g−1)(g+beta−3)=1/4` exactly,

so J has coefficient pi² and the bulk radial pole coefficient is `4pi⁴`, before the action normalization and whole-word cocycle. The fused zero-jet field is

`:exp[(2d+T).(X+Xbar)]:`.

After source integration the corresponding possible local distribution is this field times the three deltas identifying the four current positions. Since the smaller proper faces are regular there, a surviving pole would be a new primitive quartic scalar contact.

However the calculation above assumes both source heights are nonzero. Rotational Pf kernels also have source-height-supported unary contributions. Generalized mixed matching contributes on same-source current overlaps. The finite input tails contribute one-tail/one-source terms, with the net `(+1,−1)` coefficients stated in section 1. The source-pair P/L forest also belongs to the prescribed full expression. These terms can change or cancel the full residue. Positivity of the open-height bulk sectors does not prohibit such cancellation.

For constant source switching, the Y0 mixed functional itself simplifies: its normal a0 derivative cancels C, each spectator current derivative cancels the corresponding regular-current contraction, and an unlike source leaves `(h−t)/(x−w2)` plus its barred counterpart. This vanishes after restriction h=t, including the derivative-field crossmatchings. It is useful bookkeeping, but it does not remove the remaining unary and input-tail terms from the resonance test.

The unresolved quantity is therefore the residue of the **entire** once-expanded physical-current coefficient, not the nonzero number `4pi⁴` above. If that residue survives, a finite part at the resonant coupling would require an additional local quartic normalization; the old pair first-jet prescription alone would not specify its finite value. If it cancels, the same meromorphic grouping identity extends across this point. This report does not choose between those alternatives without the missing complete calculation.

## 7. Evidence and retained failures

`check_two_barred_reduced.py` prints **23 PASS** checks for the complete one/two-source kernels, residual differential term, independent normal jets, determinant sign, contour-residue polynomial, orthogonal-current coefficient, resonance, and source-containing high-g face degrees.

`check_diamond.py` completes **3 PASS** checks for the full generating-coefficient diamond and independent local balancing. It retains an explicitly expected **FAIL** for the two-extreme global-average shortcut. Its final node-0 assertion used SymPy structural equality rather than algebraic cancellation and raised an assertion error; the separate fresh `check_node0.py` correctly subtracts and cancels the expressions and prints **1 PASS**. The source, partial output and explanation of the faulty assertion are all retained.

An initial unreduced-symbol script was interrupted during a long initial simplification, before producing a result. Its source and empty output are retained; the completed script fixes translation coordinates and uses rational cancellation. No failed or interrupted evidence was deleted.

The separate `check_equal_bulk.py` adds **3 PASS** checks for the complete equal-source contour polynomial and its physical coefficient.

Thus there are **30 completed PASS checks**, one intentional failed shortcut, one retained test-harness assertion error corrected separately, and one retained interrupted simplification. The positive grouping conclusion uses the explicit distributional argument above; it does not infer a full resonant residue from these finite algebraic checks. No book change or unrestricted theorem label is recommended.
