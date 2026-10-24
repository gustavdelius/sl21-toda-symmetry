# Regulated mixed face with a separated second Yukawa source

Date: 2026-10-24. Model: Codex (GPT-6). Reviewed commit: `745586e6a6f7bca5ed91a1e9cf89d6ef5df387f0`.

This is a narrow independent check of the local mixed-face step entering H8 of `briefs/referee_sl21_coord_stage1_kernels_report.md`. It does not repeat or certify that report's global continuation argument. No existing report, earlier scratch output, book file or index was edited. No subagents, git mutations, worktrees, commits, pushes or renders were used.

## Result

**The separated-source mixed contact is divisible by h-t after subtracting the expected derivative of the h-dressed endpoint.** The factor is derived from the actual current poles, strip edges and both tails. There is no local remainder independent of h-t in this calculation.

Put the inner barB at x=0, the separated additional source at v, beta=t-1, e=h-t and

\[
R(v)=\frac1v+\frac1{\bar v},\quad
Q_{2,h}=|v|^{2h-4}:U(0)Y_0(v):_{\rm fused},\quad
U=e^{(d+a_2)(X+\bar X)}.
\tag{J1}
\]

The subscript means the normal exponential after all scalar contractions have been extracted; it does not multiply the same Wick factor a second time. Thus Q is the actual regulated UY0 product. The notation partial_x below differentiates the endpoint location in both chiralities, keeping v fixed.

With the common factor 2 i nu lambda0 lambda2 suppressed, the node-2 mixed face gives

| sector | local contact |
|---|---|
| double action | beta [t partial_x Q2,h + e R Q2,h] |
| both B tails | -t [t partial_x Q2,h + e R Q2,h] |
| raw sum | -t partial_x Q2,h - e R Q2,h |

The first term in the raw sum is exactly the ordinary negative Y0 action dressing of the raw first-order endpoint t lambda2 U. The second term is the explicit regulator violation.

At the node-0 face, with Y2 separated instead, define Q0,h by replacing U and Y0 in (J1) with U0=exp[(d+a0)(X+barX)] and Y2. The contact is

\[
\boxed{\mathcal C_{0,h}=-2i\nu\lambda_0\lambda_2
(h-t)R(v)Q_{0,h}.}
\tag{J2}
\]

There are no node-0 B tails or finite matching at this separated face. It vanishes identically at h=t, including field jets. This extends the zero obtained when the second source was absent.

| Item | Verdict | Scope |
|---|---|---|
| Node-2 regulated mixed face has no h-independent extra remainder | correct | explicit local action/tail calculation below |
| Node-0 zero persists off the diagonal without a correction | wrong | the nonzero term (J2) is present for h unequal to t |
| Node-0 contact vanishes at h=t | correct | cancellation of derivative and field-jet terms |
| Finite matching dressing adds only h-t regulator terms | correct for either explicit convention in section 4 | actual product differentiation |
| Full H8 and its generic diagonal continuation | unclear from this narrow check alone | overlapping union and continuation regularity are not evaluated here |

## 1. Direct mixed-pole expansion

Use the plan's bilinear covariance and cocycle. The source-pair regulator is Rh=|v-w|^(2(h-t)); all other powers and derivative contractions remain physical. Suppress A, which may be kept as a smooth separated Wick spectator. At the node-2 face, B is at z, barB at 0, Y2 at w and Y0 at v. The double-action coefficient is nu lambda0 lambda2/pi squared, after the ordered cocycle. Before the moving-current integration the relevant scalar and derivative factors are

\[
\frac{|v-w|^{2h-2}}{(z-v)(z-w)(-\bar v)(-\bar w)}
\left[C(z)+\frac1{z-v}+\frac\beta{z-w}\right]
\left[\bar C(0)-\frac1{\bar v}-\frac\beta{\bar w}\right].
\tag{J3}
\]

Every polynomial is inside the common normal product. Define

\[
K(w,\bar w)=\frac{|v-w|^{2h-2}}{(w-v)(-\bar v)}
:e^{(d+a_2)X(w)+a_0X(v)+d\bar X(0)+a_2\bar X(w)+a_0\bar X(v)}:.
\tag{J4}
\]

The moving z-pole coefficients computed from (J3) are

\[
D_2=\beta K,\qquad S_2=-\partial_w K+
\frac{e}{w-v}K.
\tag{J5}
\]

For the scalar simple pole, expansion of (z-v)^(-1) and of its explicit derivative contraction gives (2-t)/(w-v), whereas negative differentiation of K gives -(h-2)/(w-v). Their difference is e/(w-v). The field-vector identity is a0+beta d=-(d+a2), so the normal-field jets obey (J5) as well. This is where h-t appears before integration.

At w=0, K=Q2,h. Write Q_w and Q_bar for the two derivatives of the fully dressed endpoint. Then

\[
K_w=Q_w=\left[(d+a_2)\partial X-\frac{h-2}{v}\right]Q,
\]
\[
K_{\bar w}=\left[a_2\bar\partial\bar X-\frac{h-1}{\bar v}\right]Q,
\quad
Q_{\bar x}=\left[(d+a_2)\bar\partial\bar X-\frac{h-2}{\bar v}\right]Q.
\tag{J6}
\]

The distinction between K_bar and Q_bar is essential: barB remains at zero in K, while differentiation of the endpoint moves its d component too. Put Dbar=barC-1/barv. Direct vector multiplication, beta a2-a0=t(d+a2), gives

\[
\beta K_{\bar w}-\bar DQ=tQ_{\bar x}+\frac e{\bar v}Q.
\tag{J7}
\]

Dropping the Y0 pole derivative or replacing K_bar by Q_bar prematurely would lose one of the two regulator terms.

## 2. Strip edges and both tails

Take upper-minus-lower horizontal B contours and an inner source circle of radius epsilon, with heights of magnitude greater than epsilon as in X14. Cauchy-extend moving-current poles first. With Theta the strip, (J5) yields the source integrand

\[
-2\pi i\left[\beta(\partial_w\Theta)K-
\Theta\partial_wK+\frac{e\Theta}{w-v}K\right]
\left[\frac\beta{\bar w^2}-\frac{\bar D}{\bar w}\right].
\tag{J8}
\]

The actual strip integrals are I=pi, Jw=-pi and J1=-pi, with the definitions used in the earlier one-source calculation. Keeping all first jets of K, the integral in brackets before the -2 pi i factor is

\[
\pi\beta\left[-tK_w-\frac e vQ+\bar DQ-\beta K_{\bar w}\right].
\]

Multiplication by the double-action coefficient and use of (J7) gives

\[
\mathcal C^{\rm action}_{2,h}
=2i\nu\lambda_0\lambda_2\beta
\left[t\partial_xQ+eR(v)Q\right].
\tag{J9}
\]

This calculation includes the strip-edge derivative of the double pole. Treating the source as a constant spectator would not give the e term correctly.

For the B barK2 Y0 tail, the action times tail coefficient is +t lambda0 lambda2/pi before its cocycle, and -t nu lambda0 lambda2/pi afterward. Its moving double and simple coefficients are beta Q and -Q_w-eQ/v. The contribution is

\[
\mathcal C^{(1)}_{2,h}=-2it\nu\lambda_0\lambda_2
\left[Q_w+\frac e vQ\right].
\tag{J10}
\]

For K2 barB Y0 the coefficient is +t nu lambda0 lambda2/pi. Its antichiral double pole is -beta and its simple coefficient is Dbar. The opposite antiholomorphic Cauchy boundary value gives

\[
\mathcal C^{(2)}_{2,h}=2it\nu\lambda_0\lambda_2
\left[\bar DQ-\beta(Q_w+K_{\bar w})\right].
\tag{J11}
\]

Use (J7) to add (J10) and (J11). The result is -2it nu lambda0 lambda2 [t partial_x Q+e RQ]. Combining it with (J9) proves the table in the outcome. These are complete local first jets of the uncontracted fields, not only scalar homogeneities.

## 3. The node-0 face

Now place Y0 at w and Y2 at v. The moving-current double coefficient is K0 and the simple coefficient is

\[
S_0=\partial_wK_0-\frac e{w-v}K_0.
\tag{J12}
\]

K0 has the same scalar factor as (J4), with a0 and a2 exchanged in the source momenta. The derivative current C=a0 partial X remains unchanged. Its scalar simple coefficient is (t-2)/(w-v), while the derivative of K0 gives (h-2)/(w-v). Thus the negative e term is independently fixed. The field part is a0+d, as required.

The anti pole has coefficient 1/barw squared-Dbar0/barw, with Dbar0=barC-beta/barv. The double/simple holomorphic derivative terms cancel because Jw+I=0. The remaining local integral is

\[
\pi\left[\frac e vQ_0+\bar D_0Q_0-(K_0)_{\bar w}\right].
\]

But (K0)_bar=[a0 bar partial barX-(h-1)/barv]Q0, so Dbar0 Q0-(K0)_bar=eQ0/barv. The overall factor gives exactly (J2). At fixed separated v, the tail with node 2 collapsed onto a current belongs to the other mixed face, not this one, and must not be counted a second time here.

## 4. Finite matching and the remaining convention choice

Let the first-order finite matching be alpha lambda2 U, with alpha=t(t-1) in the canonical prescription. Its local current contact is 2 pi i nu alpha lambda2 partial_x U. If Rh is applied to the **insertion-wise** ordinary Y0 dressing of this contact, then

\[
R_h\,\partial_x[U(x)Y_0(v)]_{t}\big|_{x=0}
=\partial_xQ_h+eR(v)Q_h,
\tag{J13}
\]

since partial_x Rh at x=0 equals -e R(v)Rh. Including the negative source action, the added contact is

\[
-2i\nu\lambda_0\lambda_2\alpha
[\partial_xQ_h+eR(v)Q_h].
\]

The node-2 contact then becomes

\[
-2i\nu\lambda_0\lambda_2
\left[(t+\alpha)\partial_xQ_h+(1+\alpha)eR(v)Q_h\right].
\tag{J14}
\]

At the canonical value the derivative coefficient is t squared. If instead the finite matching is *defined* to differentiate the whole regulated zero-form Qh, its extra e alpha RQ term is absent, and the residual coefficient in (J14) is 1 rather than 1+alpha. These are two explicit off-diagonal extensions of the same physical finite matching. Their difference is itself proportional to h-t. This check therefore supports the local divisibility needed in H8 under either convention, but the common family must state which convention it uses. The two definitions cannot be conflated off the diagonal.

No counterterm denominators or overlapping-source limits were used in this argument: v was kept apart throughout. In particular, local divisibility does not alone prove that an integrated meromorphic coefficient has no pole on h=t. That is the remaining analytic step in the full H8 argument.

## 5. Optional first outer-test derivative at BB order lambda2

For a moving-B test with local linear term g1, the bare one-Y2 action coefficient of g1 times the fused U is proportional to beta squared (I+Jw), hence zero. The two tail double-pole terms are equal in magnitude and opposite in sign and cancel. The same conclusion follows with the natural horizontal spatial extension of the test: the additional strip-edge terms multiply a horizontal double-pole integral, which vanishes. With the minimal finite matching proportional to partial_x U delta and no extra derivative-of-test distribution, its g1 coefficient remains zero.

Thus the actual minimal-family BB coefficient of lambda2 U(y) partial_x delta(x-y) vanishes at this order, just as the independently calculated AA first-test jet does. This is a computation in that family. An integrated relation with outer g=1 alone would not have fixed this coefficient in a different extension. It does not construct arbitrary-test products of all orders.

## 6. Evidence, limitations and integration

Fresh evidence is `referee_scratch/sl21_symmetry_coordinator_20261024/stage1_contact_bridge/check_bridge.py` and `check_bridge.out`. The code reconstructs the current Laurent coefficients, verifies vector identities and full antiholomorphic jets, combines the strip integrals and tails, checks matching differentiation, and probes the explicit regulator factor numerically. It imports no old scripts. All 19 checks pass. At t=-sqrt(2)/3 and v=2+0.7i, the raw residual divided by 2 i nu lambda0 lambda2 is -0.000448666885492, -0.000109668102818 and -0.0000218022512840 for e=0.02, 0.005 and 0.001; its ratio to e tends smoothly to a finite value. This is a numerical corroboration of the exact coefficient, not a continuation theorem.

The strongest verified statement is the local, separated-source divisibility calculation (J9)–(J14), together with the node-0 term (J2). Restoring a separated A only multiplies the endpoint by its smooth Wick factor and inserts its ordinary pole derivatives in the total derivative; Rh is independent of a, so its additional derivative terms still carry e. No assertion of uniformity when A or v approaches the mixed cluster is made here.

**Integration:** none. Verified independently: the regulator-dependent node-2 action and both-tail contacts, the full field jets, the node-0 vanishing on the physical diagonal, finite-matching divisibility and the local BB first-test-jet cancellation. Clarified: insertion-wise and whole-zero-form matching differ off the diagonal by an explicit regulator term. Unresolved in this report: integration across the overlapping source/mixed union, regularity of the full diagonal restriction, and the complete projected mass-squared identity. No existing claim was silently edited or promoted.
