# Attempted common product: Gaussian label extension, derivative jets, and finite matching

Date: 2026-10-24. Model: Codex (GPT-6). Reviewed commit: `745586e6a6f7bca5ed91a1e9cf89d6ef5df387f0`.

This is one independent worker in the editor-authorized common-product investigation. No subagents were spawned. I read the repository rules and inherited reports, but no other worker's new files. All calculations below are new. Existing files, book files, indices and reports were preserved; no git mutation, worktree, commit, push or render occurred. Fresh evidence is in `referee_scratch/sl21_common_product_20261024/construction/`.

## Outcome

**Unclear as a common associative product of the canonically normalized interacting charges.** There is a constructive associative Gaussian regulator algebra before collision. It gives a nontrivial local regrouping, including the first two composite/test jets and the derivative mixed contact. It is substantially more than an assignment of scalar exchange phases. However, its descent to the specified interacting collision products remains unproved.

A tempting proposed shortcut is **wrong**: a near-identity change of composite normalization, transported from an associative raw product while fixing the single currents, cannot turn the first mixed coefficient from t into t². The free mixed braided relation makes that change vanish identically. This is an obstruction to that specific construction, not to all associative extensions.

| Item | Verdict | Independent evidence and exact scope |
|---|---|---|
| One associative auxiliary Gaussian product can generate both g and h factors | correct | Explicit bilinear label covariance and cocycle calculation below; separated Wick algebra only |
| Retaining constituent labels suffices to transport the BB/U derivative jets | correct | Exact three-factor Wick polynomial, fused first descendant, second pair jet and arbitrary smooth-test jets |
| A fused U may be treated as an elementary unlabelled endpoint off the physical surface | wrong | It would replace the required power g−1 by rho−1 and lose a derivative-regulator term |
| The mixed derivative contact is the derivative of the entire g-regulated endpoint product off the surface | wrong for insertion-wise differentiation | Explicit difference `(g−rho) W_g(B,U)/(z−x)` |
| A fixed-single-current, near-identity transport of a raw associative product creates the canonical finite matching | wrong | Direct first-order associator/coboundary computation |
| A common rescaling by sqrt(t) converts both raw mixed amplitudes to the canonical values | correct conditionally | Direct algebra; changes the leading single-current normalization and presupposes a raw associative realization |
| Continuing all n-ary kernels proves binary associativity | unclear without descent | Intermediate collision restriction need not be multiplicative; a precise ideal criterion is given below |
| All mass-order-two common products exist in the required normalization | unclear | Full contact cochains and the mixed/Serre union-face restrictions have not been constructed in this report |

No physical or global completion is used, and no theorem label is proposed.

## 1. A real associative starting product

Use the physical momenta and conventions

`d²=−rho`, `a0.d=a2.d=1`, `a0²=a2²=−1`, `a0.a2=1−t`, `rho=2/t`,

with `B=:C exp(d.X):`, `C=a0.partial X`, and `A=:exp(−d.X):`. The mixed endpoint is `U=:exp[(d+a2).(X+Xbar)]:`. Current and source labels are kept separately from physical momenta.

Introduce independent auxiliary Gaussian label fields. There is one label for each chirality/species of current. Its contraction with itself is chosen to contribute `(g−rho) log(z−w)` for a left current, or the barred logarithm for a right current. Different current-label fields have zero contraction. There are also two source-label fields, with cross contraction

`(h−t)[log(z−w)+log(zbar−wbar)]`,

and zero self contractions. These signs are definitions of a formal complex Gaussian contraction; no positivity is asserted. Attach the appropriate auxiliary exponential to each original constituent. Tails and fused fields retain all their constituent labels.

For label vectors p,q define the bilinear extra contraction b(p,q). On separated points the binary product is the ordinary Gaussian Wick rule

`E_p star E_q = epsilon(p,q) exp(b(p,q)+physical contraction) :E_p E_q:`,

extended to derivatives by differentiating the physical contractions. The original derivative C acts on the physical fields, not on the auxiliary label fields. The physical cocycle remains

`epsilon(n,m)=(-1)^[(nL+nR)(m0+m2)] exp[-i pi(nL mR−nR mL)/t]`.

This is an explicit binary map on decorated Wick functionals. Its associativity follows by direct expansion: b(p,q)+b(p+q,r)=b(q,r)+b(p,q+r); the logarithm of the cocycle has the same property. Wick contractions involving derivative fields are obtained by differentiating that identity. The script checks both bilinear identities independently.

Taking the auxiliary Gaussian vacuum coefficient yields precisely the factors used in the prior prescriptions: one `(z_i−z_j)^(g−rho)` per like left-current pair, its right analogue, and one `|w0−w2|^(2(h−t))` per unlike source pair. This vacuum coefficient is a linear operation, not a general algebra homomorphism. Thus it must not be applied to an intermediate block while throwing away its auxiliary labels. On the physical surface all auxiliary contractions vanish; setting their exponentials to one is then a legitimate algebraic quotient before singular collision extensions. Whether it remains legitimate after those extensions is the descent issue in section 5.

This construction is a Gaussian field calculation, not a free algebra with the desired mixed or Serre relations imposed. It produces no canonical finite mixed contact by itself.

## 2. An explicit regrouping with a mixed endpoint

Keep U as the fused B/barB/Y2 endpoint, so it carries one left B label and one right barB label. Put two left B currents at z1,z2, and U at x. Set

`r12=z1−z2`, `ri=zi−x`, `eta=t−2`.

The contractions `−a0.(d+a2)=eta` and `−d.(d+a2)=rho−1` give the regulated scalar vertex factor

`r12^g r1^(g−1) r2^(g−1)`.

The complete polynomial in the uncontracted physical fields is

`P=(C1−1/r12+eta/r1)(C2+1/r12+eta/r2)+1/r12²`.

Here the last term is the C1C2 contraction. The zero-order B1B2 polynomial is `(C1−1/r12)(C2+1/r12)+1/r12²`. Multiplying that composite by U shifts each Ci by eta/ri and gives P. In the other grouping B2U first gives C2+eta/r2. Multiplication by B1 shifts C2 by +1/r12, supplies the contractions of C1 with the B2 and U exponentials, and adds the C1C2 contraction. It again gives P. Thus both binary groupings reproduce the complete Wick polynomial, not merely its exchange phase.

There is also a genuine composite-jet check. Write z1=y+r/2, z2=y−r/2, R=y−x and expand the fields. The leading BB composite is

`:(C²+C') exp(2d.X):`,

not `:C² exp(2d.X):`. Contracting it with U gives

`R^(2g−2)[(C+eta/R)²+C'−eta/R²]`.

Directly expanding P gives precisely this expression. The scalar coefficient is eta(eta−1)/R². The `−eta/R²` comes from contracting C' with the U exponential; dropping derivative descendants would fail this test.

The script also reconstructs the coefficient P2 of r² separately from both exact groupings. It includes the normal-field jet `1+r²(d.partial²X)/4`, the external regulator jet `1−(g−1)r²/(4R²)`, and arbitrary tests f(z1)k(z2). In particular the linear and quadratic test coefficients are

`(f'k−fk')/2`,

`(f''k+fk''−2f'k')/8`.

Both groupings retain them. All these are symbolic identities in the test and field jets; no constant-test substitution is used. The explicit P2 is recorded in `derived_jets.txt`.

This is relevant to the canonical contact at order lambda2: after a proper mixed B/barB/source face is collapsed, its propagation against remaining B currents uses exactly such labelled U products. It does not evaluate the full simultaneous mixed union or specify how all possible proper-face contacts add there.

## 3. The derivative contact must also retain its labels

Before regulation, insertion of the contact proportional to partial_x U means differentiation of the physical U field. Put r=z−x and let W_g(B,U) denote its g-regulated Wick product. Since the extra factor is r^(g−rho), direct differentiation gives

`W_g(B,partial_x U) = partial_x W_g(B,U) + (g−rho) W_g(B,U)/r`.

The derivative on the right is the full spatial derivative, including the normal fields. The extra term is absent only if one deliberately chooses a different off-surface contact convention that differentiates the entire decorated endpoint. These conventions agree at regular physical-surface restrictions; they cannot be conflated before checking regularity.

For the scalar B/U term, W_g=eta r^(g−2). The insertion-wise derivative is

`−eta(rho−2) r^(g−3)`,

whereas the entire-product derivative is `−eta(g−2) r^(g−3)`. Their difference is exactly `(g−rho)eta r^(g−3)`. The script checks this with an arbitrary smooth outer test as well: its differentiated-test contribution is retained on both sides.

On this one relative coordinate the Mellin divisors are `g−2+j=0`. At generic g=rho none is identically the physical surface, so this particular evanescent difference restricts to zero. That argument cannot be transferred automatically to a larger union whose degree involves g−rho itself.

The first unresolved BBbarB/Y2 union has degree g−2. Independently counting its three current coordinates, two source coordinates, three current/source powers and three derivative contractions gives that result. For B³AbarB/Y2 the same count gives `3(g−rho)−1`; for its unlike-double-source channel it gives `3(g−rho)+2h−4`. The first Taylor level of the former can therefore have a physical-surface divisor. This is a requirement for a separate five-current calculation, not a demonstrated residue or obstruction. The coordinator assigned that calculation independently; this report makes no claim about its outcome.

## 4. Why a single normalization transport does not produce the canonical matching

Suppose an associative raw product has actually been constructed:

`mu=mu0+lambda2 mu1+O(m²)`.

This is a conditional assumption here. At zeroth order the left and right Gaussian factors are independent, so for B and D=barB the exact free relation is

`mu0(B,D)=s mu0(D,B)`.

There is no left/right oscillator contraction that could introduce a density-diagonal correction to this particular free relation. The scalar cocycle fixes s.

Try to construct the new product by one normalization map

`Z=Id+lambda2 Z1+O(m²)`, `mu'=Z^{-1}mu(Z tensor Z)`.

The first product change is

`D_Z(X,Y)=mu0(Z1 X,Y)+mu0(X,Z1 Y)−Z1 mu0(X,Y)`.

If the single-current representatives stay fixed, Z1 B=Z1 D=0. Linearity and the free relation then force

`D_Z(B,D)−s D_Z(D,B)=−Z1[mu0(B,D)−s mu0(D,B)]=0`.

Thus no choice of Z1 on higher composites can change the leading mixed defect t lambda2 U into t² lambda2 U. The required difference is `t(t−1)lambda2 U`, nonzero at generic t. This argument applies equally to the density contact, with the prescribed derivative and tests carried along. It is not a counting argument and does not assume a choice of composite basis.

This rules out one proposed route under its explicit hypotheses. It does not prove that the canonical product is nonassociative. Possible remedies are a genuine deformation of multiplication whose first-order class is not this coboundary; changes of the single currents; or changes of zeroth-order endpoint normalization. Enlarging to unreduced constituent-labelled symbols can also evade the displayed free relation, but then descent to the physical Gaussian quotient must be proved. It cannot be assumed after using the enlarged space to create the contact.

There is a useful allowed-normalization alternative. The raw BB amplitude is t; the raw AA amplitude is gamma²/t, since gamma+gamma²=gamma²/t with gamma=t/(1−t). Multiplying each of A,B,barA,barB by the same chosen square root of t converts the amplitudes to t² and gamma². Each quartic is multiplied by t², so its zero is preserved. This would give an isomorphic normalized realization if a common raw associative realization existed. It changes the fixed leading one-current normalization, so it is not an implementation of the requested finite matching with the original currents held fixed. No global or positivity interpretation of the square root is needed for this formal statement.

## 5. What is still needed for binary composition

At separated coordinates the decorated Wick star product is a genuine associative binary multiplication. Introducing a smooth ultraviolet cutoff into its contraction also gives a genuine finite-cutoff associative algebra. A formal Dyson interaction inside that algebra would remain associative coefficient by coefficient. None of those observations proves that its collision limits equal the particular circular-source-hole, current-first Cauchy and finite-matching prescription in the reports.

There are two distinct reasons not to replace that comparison with a slogan about a common interacting state. First, the source-circle masks in the prescribed integrals depend on the observable's insertion points. They are not manifestly a single observable-independent interaction functional. Second, Euclidean ordered insertions and a changed state on the unmodified free operator algebra are not automatically the same interacting operator product. One must specify the interaction map or prove the ordered-kernel factorization. I have not supplied that missing map here.

A precise conditional descent criterion is available. Retain meromorphic germs in `(g,rho,h)` together with all constituent and derivative-jet data. Define the binary multiplication before physical-surface restriction by the actual common kernels. The ideal generated by `g−rho` and `h−t` must be a two-sided ideal of the subspace on which the physical restriction is taken. Equivalently, changing an intermediate representative by a vanishing regulator factor must not generate a nonzero physical result on subsequent multiplication. In addition, the generators and all products being used must have well-defined restricted distributions, and the kernel composition must agree on separated blocks.

If the relevant multiplication has no physical-surface poles and factorization is established, this ideal condition follows and associativity of the pre-restriction product descends. If it has such a pole, an evanescent intermediate term can survive. The elementary model is `(g−rho) times 1/(g−rho)=1`: taking separate finite parts first gives zero. This is not an actual anomaly calculation; it explains exactly why assigning an n-ary finite part and naming it an iterated binary product is insufficient.

Through mass order two one also needs a full first contact map D, not only its constant-test braided projection. Its first associator is

`D(XY,Z)+D(X,Y)Z−D(X,YZ)−X D(Y,Z)=0`.

At lambda2 squared the same expression for a second-order map E must cancel

`D(D(X,Y),Z)−D(X,D(Y,Z))`.

At lambda0 lambda2 the source dressing, P subtraction and current improvement supply additional components. The known coefficients t/h and 1/h constrain these maps but do not give their values on every composite and derivative test. The local checks in sections 2–3 provide actual entries of this problem; they do not establish all these equations.

The next concrete deliverable is therefore a factorization-compatible interacting map, or a complete finite-order contact cochain on the Gaussian descendant space, together with the mixed/Serre union calculation. Calling the enlarged meromorphic labels the physical algebra before these steps would conceal the unresolved question.

## 6. Further constructive attempts and the failed pinch shortcut

### A five-current chamber calculation does not by itself determine a full-word contact

Following the coordinator's request for an independent check, I reconstructed the B B A U scalar kernel directly. Set the endpoint x=0, the B positions z,v, and the A position a. Its polynomial is

`H=[−1/(z−v)+1/(z−a)+eta/z][−1/(v−z)+1/(v−a)+eta/v]+1/(z−v)²`.

At rho t=2 it obeys the exact identity

`(zv/a)H=(eta/rho)[(2−rho)/a+rho/(z−a)+rho/(v−a)]`.

Thus the full scalar kernel is `(eta/rho) partial_a P`, where on the chamber 0<a<z<v,

`P=(v−z)^g(z−a)^(-rho)(v−a)^(-rho)z^(g−2)v^(g−2)a^(2−rho)`.

In the coordinates v=R,z=uR,a=uyR the Jacobian is uR². The sum of the two derivative-regulator factors is `−(1+u)/(uR)`. Moving partial_a onto the first smooth test jet therefore exposes a radial factor `(g−rho) integral R^(3(g−rho)−1)dR`, whose limit is 1/3, multiplying a nonzero chamber angular integral. The angular shape, apart from the overall sign and eta/rho, is

`(1+u)u^(g−2rho)(1−u)^g y^(2−rho)(1−y)^(-rho)(1−uy)^(-rho)`.

At g=rho and −1<rho<0 this is positive and integrable. The independent script verifies the full Wick identity, the Jacobian, and all powers.

**This is a single spatial-chamber residue, not a full ordered-word ambiguity.** The coordinator subsequently identified an important cancellation mechanism: with the fused endpoint retaining all original imaginary-height orders, the relevant contact is a purely chiral ordered boundary product. Its chambers cannot be assigned independent restrictions. The directed Gaussian boundary graph has no positive circulation. In a Fourier–Laplace representation, a height potential bounds the positive edge parameters by the norm of their incidence image; hence pushforward is proper on compact Fourier sets. Using the entire normalized one-edge distributions can then cancel chamber Mellin poles in the complete word. I did not independently reconstruct all contact masks required for that global conclusion, so I do not elevate the chamber residue to an obstruction or claim the cancellation proves the full interacting product.

A finite full-diagonal derivative contact could cancel one chamber difference; its fused momentum would be `(2d+a2,d+a2)`. Merely choosing that correction is not a construction of an associative binary map. It must agree with the complete-word extension and the other overlap equations. This distinction prevents both a premature no-go and a premature existence claim.

### Direct canonical equal-time multiplication

There is a second useful route that avoids trying to transport the radial raw product. Start with the genuine Gaussian operator multiplication, use the densities

`J_B−t lambda2 K2`, `J_A−gamma lambda1 F_A`,

and their barred counterparts, and attempt their common equal-time boundary values. Their tail-current terms already give the canonical mixed coefficients, without adding the radial bare-current/action contact.

This can be checked from vectors. For v=d+a2 and beta=t−1,

`a0+beta d=−v`, `beta a2−a0=t v`.

The J_B barH_B double pole has coefficient t beta U and simple vector −t partial U. The H_B barJ_B term has the same double coefficient and simple vector `t beta partial U+t² partialbar U`. Opposite holomorphic/antiholomorphic boundary values cancel the double-test terms and add the spatial derivative with coefficient t². The AA computation similarly gives `gamma[(1−t)^(-1)−1]=gamma²`. These vector and coefficient identities are independently checked in `check_followup.py`.

At mass order two the mixed tail pair K2 barK2 is particularly simple: `(d+a2).a2=0`, in each chirality, and its cocycle ratio is s. It has no oscillator pinch. Two fused U endpoints do have opposite boundary factors, `(r+i0)^(g−1)(r−i0)^(g−1)`. Their undifferentiated direct limit exists for Re g>1/2. Its continued poles, before derivative descendants, occur at `2g−1+j=0` for the allowed Taylor levels. They do not contain the physical surface identically; the positive patch 2<t<3 has rho between 2/3 and 1 and is locally integrable. This confirms that an apparent directed cycle is not automatically an obstruction.

I also tested the apparently harder configuration of two K2 tails and one bare barB. Before expanding its polynomial, the two chiralities can give a directed triangle, and the full B³AbarB/two-Y2 degree is `3(g−rho)`. The degree count is correct, but treating the unexpanded graph as proof of a pinch is wrong. Put x=barz_i−baru and y=barz_j−baru. Its antichiral exponential factor and derivative contraction satisfy

`(x−y)/(xy)=1/y−1/x`,

`−beta (x−y)(1/x+1/y)/(xy)=beta(1/x²−1/y²)`.

Thus both the uncontracted barC term and its scalar contraction split into single-edge boundary distributions. This independently verified partial fraction removes that proposed cycle. A further worker was assigned the resulting canonical construction and its comparison with the action-dressed products; this report does not anticipate that verdict.

The direct canonical route therefore remains viable and is stronger than the failed composite-normalization transport. Its outstanding task is to prove common equal-time limits and identify their perturbative interaction/source dressing with the prescribed Euclidean products. Keeping action insertions only in a purported changed state, or declaring their difference a normalization choice, does not supply that identification. Conversely, the coboundary obstruction of section 4 does not apply to a genuinely different canonical realization unless one first proves it is related to the raw product by the restricted normalization map used there.

## Evidence and integration

The independent script `check_construction.py` imports no previous research code. Its output contains 21 PASS lines and no failures. The follow-up calculations contain 15 further PASS lines, and the independent pinch-resolution check contains two PASS lines; no check failed. It checks the auxiliary contraction and cocycle, the complete two regroupings, the first and second fused/test jets, mixed derivative transport, the normalization obstruction, the simultaneous rescaling, and the relevant collision degrees. These checks are supporting evidence for the written derivations, not a test suite for an already constructed quantum algebra.

No book integration occurred. Verified: a genuine separated Gaussian binary product, a nontrivial local descendant regrouping, the derivative-regulator difference, and the fixed-single-current transport obstruction. Modified as a proposed strategy: keep all constituent and jet data and require descent explicitly. Rejected: obtaining canonical mixed matching by a near-identity transport that only renormalizes composites. Unresolved: the actual canonical common associative collision product through mass order two. This is not a nonexistence result.
