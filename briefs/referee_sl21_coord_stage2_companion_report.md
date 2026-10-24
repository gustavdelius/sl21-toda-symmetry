# Companion three-contour projection at the even-root coupling

Date: 2026-10-24. Model: Codex (GPT-6).

Commit reviewed: `745586e6a6f7bca5ed91a1e9cf89d6ef5df387f0`.

Scope: an independent calculation of the coefficient linear in lambda1 in [BA-s^(-1)AB,barA], followed by a local first-test-jet calculation and covariance checks on a prospective derivative contact. Completed stage-1 outputs were not changed. No subagents, book edits, git mutations, worktrees, commits, pushes or renders were used. Fresh evidence and input hashes are in `referee_scratch/sl21_symmetry_coordinator_20261024/stage2_companion/`.

## Result and verdicts

In the restricted contour prescription specified below, the actual kernels give

\[
\boxed{[B_f A_g-s^{-1}A_g B_f,\bar A_I]\big|_{\lambda_1}
=(1-s^{-2})B_f R_A^{(1)},\qquad
R_A^{(1)}=\frac{\nu}{2\pi i}\gamma^2\lambda_1 F_A(R).}
\tag{C1}
\]

Here I=[L,R], f is smooth with compact support strictly inside I, and g=1 near I. **The compactly supported spectator is now B, and the outer constant test is attached to A.** This differs from the earlier barB projection. F_A=exp[-k(X+barX)], gamma=t/(1-t), rho=2/t and s=nu squared. The statement is the lambda1 coefficient only, not all terms of mass order two. It is a formal finite-segment identity in a specified extension, not a quantum symmetry algebra.

| Claim | Verdict | Independent evidence |
|---|---|---|
| Actual B A barA / V1 kernels and both tails | correct | direct Gaussian contractions and ordered cocycles |
| Mixed face with the derivative B spectator | correct | double/simple coefficients V and eta V', including field jets |
| Raw/canonical endpoint coefficients gamma squared + gamma and gamma squared | correct in the specified extension | strip-edge and tail calculations |
| No additional primitive triple contact for Re rho<0 | correct for this restricted one-source extension | proper-face classification and nested remainder bound |
| Projected companion (C1) | correct within that scope | evaluated mixed faces, separated endpoint phases and the four contour words |
| Arbitrary-test first outer derivative contact in this AA calculation | zero in the minimal prescription used here | explicit cancellation of the action and tail first jets |
| lambda2 U(y) partial_x delta(x-y) automatically forbidden by dimensions | wrong | exact free scaling weights |
| That contact is an unconstrained change of the existing canonical/common product | unclear | further normalization, curve-deformation and overlap constraints are needed |
| One associative product supporting this and the mass-squared/quartic calculations | unclear | not established by their separate continuations |

## 1. Fix the actual words and kernels

Use the plan's three-dimensional momenta, bilinear covariance, cocycle and density J+H. Direct contraction gives d.T=2, a0.T=-t, and k=(1-t)d. Set eta=1/(1-t), so eta-1=gamma and eta k=d.

Place the spectator B at a, the moving A at z, the inner barA at u and the even-root source at w. The reference word is B(a)A(z)barA(u)V1(w). Its cocycle is +1. The bare action coefficient is -lambda1/pi. Its complete kernel is

\[
-\frac{\lambda_1}{\pi}
\frac{(a-z)^{-\rho}(a-w)^2}{(z-w)^2(\bar u-\bar w)^2}
\left[C(a)+\frac1{a-z}-\frac{t}{a-w}\right]
:e^{dX(a)-dX(z)-TX(w)-d\bar X(u)-T\bar X(w)}:.
\tag{C2}
\]

In particular the derivative contracts at a, not at the moving-current coordinate. The contractions with A and V1 are respectively +1/(a-z) and -t/(a-w).

There are exactly two tails at this coefficient. They are

\[
-\gamma\lambda_1
\frac{(a-z)^{-\rho}(a-u)^2}{(z-u)^2}
\left[C(a)+\frac1{a-z}-\frac{t}{a-u}\right]
:e^{dX(a)-dX(z)-TX(u)-k\bar X(u)}:,
\tag{C3}
\]

and

\[
-\gamma\lambda_1
\frac{(a-z)^{2-\rho}}{(\bar z-\bar u)^2}
\left[C(a)+\frac{1-t}{a-z}\right]
:e^{dX(a)-kX(z)-T\bar X(z)-d\bar X(u)}:.
\tag{C4}
\]

Both reference-order cocycles are +1. There is no B tail proportional to lambda1. The added canonical AA finite endpoint coefficient is -gamma lambda1 F_A. Common charge factors (2 pi i)^(-3) are suppressed in (C2)–(C4), and restored after taking residues.

Use the same word heights as X14, with barA at time zero:

| word | B height | A height |
|---|---:|---:|
| BAbarA | 2 delta | delta |
| ABbarA | delta | 2 delta |
| barABA | -delta | -2 delta |
| barAAB | -2 delta | -delta |

Here delta=c epsilon, c>2. The labels specify the insertion order, even though the current moved through the inner contour is now A. First extend the moving A/source poles as Cauchy distributions. Use the circle |w-u|>epsilon at the inner current. Do not add a new circular hole around the moving A before its Cauchy extension. Apply the same finite AA contact on every proper mixed face. Spectators and infrared cutoffs are held identical in all four words.

## 2. Calculate the mixed face with B present

Put u=0. The fused holomorphic scalar factor is

\[
V(w)=(a-w)^{2-\rho}
\left[C(a)+\frac{1-t}{a-w}\right].
\tag{C5}
\]

V also denotes its full normal-ordered dressing when field jets are included. Expand the actual numerator in (C2) at z=w. The double pole is V/(z-w)^2. The scalar simple coefficient is

\[
(a-w)^{1-\rho}
\left[\rho C(a)+\frac{\rho(1-t)+1}{a-w}\right]
=\eta V'(w).
\tag{C6}
\]

The equality uses rho t=2. For the full field it remains true: expansion of the moving exponential contributes -d.partial X V, whereas eta times the derivative of the fused exponential contributes -eta k.partial X V=-d.partial X V. C(a) is a spectator field at a and is not incorrectly differentiated as C(w). Thus both the normal derivative field and its contraction with the endpoint have been retained.

For an upper-minus-lower horizontal contour, the moving double and simple poles give

\[
-2\pi i\big[(\partial_w\Theta)V+\eta\Theta V'\big],
\]

where Theta is the intervening strip. The three required integrals are

\[
I=\int\frac{\Theta}{\bar w^2}d^2w=\pi,\quad
J_w=\int\frac{w\partial_w\Theta}{\bar w^2}d^2w=-\pi,\quad
J_1=\int\frac{\partial_w\Theta}{\bar w}d^2w=-\pi.
\tag{C7}
\]

These follow directly by horizontal integration outside the circle and by evaluating the two strip edges. In particular, the slice through the circle gives 2 sqrt(epsilon squared-y squared)/epsilon squared for I, and its integral is pi. At each edge the single Cauchy pole integrates to plus or minus i pi; the horizontal double-pole integral vanishes.

Multiply by the action coefficient in (C2). The holomorphic first jet has coefficient
2 pi i lambda1 (eta-1)V'=2 pi i gamma lambda1 V'. The antichiral source first jet is -T.bar partial X; J1 supplies the coefficient +2 pi i lambda1 T.bar partial X V. Since T=-gamma k, it is exactly 2 pi i gamma lambda1 times the antichiral derivative of the full fused endpoint. Thus the bare action contact is the spatial derivative of the actual dressed B F_A, with coefficient gamma.

For (C3), the moving holomorphic double pole gives coefficient +2 pi i gamma eta lambda1 times the holomorphic derivative. For (C4), the antichiral double pole has the opposite Cauchy boundary value and gives -2 pi i gamma lambda1 times the spatial derivative of its numerator. The holomorphic sum is gamma eta-gamma=gamma squared. Its antichiral part also has that coefficient because gamma T=-gamma squared k. The complete bookkeeping is

| sector | coefficient multiplying lambda1 F_A |
|---|---:|
| bare currents and action | gamma |
| both tails | gamma squared |
| raw circular product | gamma squared + gamma = gamma squared / t |
| specified finite addition | -gamma |
| canonical matching | gamma squared |

This is an actual additional-current mixed-face computation, not an assumed differentiation rule for an already renormalized binary product.

## 3. Proper faces and uniformity of this one-source extension

The B/V1 exponential zero is quadratic; its derivative contraction reduces it only to a linear zero. The B/barA pair has no oscillator pole. The A/V1 and barA/V1 double poles are the prescribed Cauchy faces. The mixed A barA / V1 face has just been calculated with its full first jets.

For the simultaneous scalar coefficient put r=a-z and p=z-w. Its holomorphic rational factor resolves exactly as

\[
\frac{(r+p)^2}{p^2}\left(\frac1r-\frac t{r+p}\right)
=\frac{(1-t)r}{p^2}+\frac{2-t}{p}+\frac1r.
\tag{C8}
\]

After the r^(-rho) and antichiral double-pole factors, its homogeneity is -rho-3. Terms retaining an uncontracted C field are one power softer. Four relative real measures at an interior simultaneous collision give degree 1-rho; at a spatial endpoint there is one additional power. The AB pair degree is -rho. The B plus fused endpoint factor has worst local behavior R^(1-rho), owing to the derivative contraction in (C5). The three-point face with B, A and V1 and barA separated has degree 2-rho. The B,V1 zero does not produce an additional divergent proper face.

For a mixed cluster of radius epsilon separated by R from B, retain its first mixed jets as in section 2. The next remainder is bounded by epsilon times two derivatives of the dressed B F_A factor. Its spectator integral is bounded by

\[
\epsilon\int_{M\epsilon}^{R_0}R^{-1-\operatorname{Re}\rho}dR
=O(\epsilon),\qquad \operatorname{Re}\rho<0.
\tag{C9}
\]

At R=O(epsilon) the direct simultaneous estimate is O(epsilon^(1-Re rho)). Smooth spectator fields and test functions have no worse degree. Consequently no extra primitive triple contact remains in this initial domain after the proper mixed face is resolved. The same estimate applies to the tails. This is a one-source calculation; it does not encounter the overlapping source-pair face identified at double-Yukawa order.

The selected extension elsewhere is the meromorphic continuation of these kernels with their fixed proper finite contact, away from poles. This does not permit adding a new primitive term only after continuation. Other separately specified extensions need comparison, rather than being identified by notation.

## 4. Evaluate the four words and the endpoint phases

The off-diagonal B/barA exchange has no contact at this coefficient: the B/V1 zero survives, and its product with the barA tail is regular in the oscillator variables. Its separated cocycle phase is s^(-1). With f_A=R_A-L_A, the two contour differences therefore give

\[
BA\bar A-\bar ABA\longmapsto B_f f_A,
\qquad
AB\bar A-\bar AAB\longmapsto s^{-1}f_A B_f.
\tag{C10}
\]

Equations (C2)–(C9) supply the local contact and exclude an extra simultaneous term; (C10) is not asserted purely from an abstract exchange algebra.

For a separated endpoint, the B F_A chiral power is alpha=2-rho. The ratio epsilon(B,F_A)/epsilon(F_A,B) is s^(-1). Assign positive imaginary separation to the first B. The oscillator ratios are exp(-i pi alpha) if B is to the right of the endpoint and exp(+i pi alpha) if it is to its left. Hence

\[
F_A(L)B_f=s^2B_fF_A(L),\qquad F_A(R)B_f=B_fF_A(R).
\tag{C11}
\]

The rational derivative contraction has the same separated exchange ratio; it is part of the full Wick product and introduces no contact between separated supports.

Combine (C10) with the coefficient -s^(-1) of AB in the defect. The result is B_f f_A-s^(-2)f_A B_f. The left endpoint cancels using (C11), and the right endpoint gives (C1). Applying the same computation to the finite addition alone gives

\[
\delta[D,\bar A_I]\big|_{\lambda_1}
=(1-s^{-2})B_f\frac{\nu}{2\pi i}
[-\gamma\lambda_1F_A(R)].
\tag{C12}
\]

The result propagates the specified normalization; it does not uniquely select it. The raw coefficient passes the corresponding test with its own raw endpoint.

## 5. What an outer-test derivative does and does not allow

There is a direct first-jet check beyond the constant outer test. Let the moving-A test have local expansion g(z)=g0+g1 z. The moving double pole differentiates g as well as the strip. The bare action coefficient of g1 V is proportional to

\[
J_w+I=0.
\tag{C13}
\]

The two tails contribute +gamma g1 V and -gamma g1 V and cancel. For a spatial rather than holomorphic continuation of g to the displaced contours, the additional first-jet strip-edge terms multiply the horizontal integral of 1/(x plus or minus i delta)^2, which is zero. Thus this local first-jet conclusion is independent of those two natural extensions. The canonical finite addition is taken in its minimal endpoint-derivative form, with no independently added derivative-of-test term. In that prescription the coefficient of lambda1 F_A(y) partial_x delta(x-y) is zero at this order. This local statement is not a construction of every arbitrary-test product or endpoint jet.

For comparison, consider the suggested BB contact

\[
K(x,y)=c\lambda_2 U(y)\partial_x\delta(x-y).
\tag{C14}
\]

It acts on an outer test g and inner test h as

\[
\langle K,g\otimes h\rangle
=-c\lambda_2\int g'(y)h(y)U(y)dy.
\tag{C15}
\]

It is invisible when g=1 near the inner support and nonzero for suitable nonconstant g. The fresh script evaluates one explicit nonconstant Gaussian test and obtains -0.626657068657750126 after suppressing c lambda2.

Dimensional covariance does not eliminate (C14). Free weights are h(J_B)=1+1/t and h(U)=bar h(U)=1/t-1/2. The Yukawa coupling has engineering dimension one. Thus lambda2 U has dimension 2/t, and partial delta adds two, exactly the dimension 2+2/t of a B/barB density product. This uses engineering dimensions, not the separate mass-order bookkeeping. U is even, has the fused momentum, and carries the required local cocycle character.

Exchange reciprocity requires the reverse defect to be -s^(-1) times the forward one. This determines the corresponding reverse contact; it does not force c to vanish. One must use the distributional identity

\[
U(y)\partial_x\delta(x-y)
=U(x)\partial_x\delta(x-y)+U'(x)\delta(x-y)
\tag{C16}
\]

when rewriting it with the field at the other variable. Neglecting the last term would give a false exchange obstruction.

At the level of a product on a specified oriented straight contour, (C15) also has an intrinsic form using arclength and its tangent derivative. Transporting the contour and tests under a rigid spacetime transformation preserves that form; the opposite B and barB charge spins cancel and U is scalar. This explains why elementary covariance and the oriented-contour exchange convention do not immediately exclude it. It is **not** a demonstration of a covariant two-dimensional local distribution whose restriction gives (C14), nor of compatibility under arbitrary contour deformations. Those stronger conditions, the canonical bare-current normalization, and overlapping products can fix or disallow c.

Accordingly (C14) is a concrete candidate contact visible to arbitrary tests and not ruled out by scaling, momentum, grading or elementary exchange. It is not asserted as arbitrary freedom in the already fixed canonical product. The actual AA first-jet calculation above selects zero in the prescription used for (C1); no coefficient of a BB derivative contact has been computed in this companion calculation.

## 6. Relation to the common-product problem and retained evidence

The four earlier same-chirality/Serre reports were inspected before making this comparison. Their continuations introduce separate exponents for like-current pairs and, at mass order two, for the unlike source pair. This companion calculation has no like-current pair and no source pair. Its initial domain Re rho<0 is a directly controlled domain for the kernels above, but it does not prove that its continuation is the restriction of those fourfold products, or that their embedded mixed faces agree with (C7), (C12) and the first-test-jet choice. No quartic relation or mixed quartic overlap was independently established in this worker task.

Fresh evidence: `check_companion.py` and `check_companion.out` in the scratch directory. The script imports no earlier script. It checks the vector products, derivative Laurent coefficients, normal-field identity, ordered cocycles, contact coefficients, exact simultaneous rational kernel, endpoint phases, scaling dimensions and local test-jet cancellations. Independent quadrature of a Gaussian-tested double pole approaches the required first-test-jet distribution, with normalized values 0.95604116, 0.97772943 and 0.98879083 at delta=0.02, 0.01 and 0.005. The explicit nested bound decreases proportionally to epsilon. There are 25 PASS lines and no FAIL lines; these diagnostics support the written derivation, not a global algebra assertion.

No existing erroneous assertion is targeted by a replacement: the earlier report explicitly left this companion calculation open. The concrete addition to that status is (C1), with the changed support assignment and its finite-order/prescription qualifications. The smallest next common-product input is a comparison of the actual embedded mixed extensions, including arbitrary-test jets, with one specified higher-order continuation. Agreement of abstract overlap formulas is insufficient.

## Integration

None. Independently verified and extended: the lambda1 companion kernel, full derivative spectator mixed contact, finite matching, collision bound and projected endpoint identity. Verified locally: cancellation of the first outer-test derivative term in the stated AA prescription. Classified but not adopted: the dimensionally allowed BB derivative-contact candidate and its reverse-exchange constraint. Unresolved: the complete arbitrary-test product, the mass-squared double-Yukawa extension, compatibility with fourfold products, all-order conservation, global endpoint limits and an established symmetry algebra.
