# Stage 0 independent audit: first triple contact and double-Yukawa faces

Date: 2026-10-24. Model: Codex (GPT-6). Reviewed commit: `745586e6a6f7bca5ed91a1e9cf89d6ef5df387f0`, branch `main`.

This is the bounded worker audit requested by the symmetry coordinator. No subagents, existing-file edits, git mutations, worktrees, commits, pushes or renders were used. Input hashes are in `referee_scratch/sl21_symmetry_coordinator_20261024/stage0_worker/inputs.sha256.json`. Evidence is the independently written `audit.py`, `audit.out` and `faces.json` in that directory. The script imports no earlier research scripts. Its output has 45 PASS lines and no FAIL lines. These checks supplement the derivation below; their count is not a proof status.

## Verdict ledger

| Claim | Verdict | Independent evidence and limit |
|---|---|---|
| X5–X7 kernels and ordering signs | Correct | Vector contractions and additive cocycle products below |
| X8–X11 first jets, strip edges and raw coefficient | Correct | Explicit pole expansion and elementary horizontal contour integrals |
| X12–X13 proper faces and spectator remainder | Correct for the specified Cauchy extension and smooth, separated spectators | All subset degrees and the explicit nested integral below; this is a finite-order distribution calculation |
| X14–X17 restricted projected identity | Correct within that prescription | Endpoint branch values and cocycles, proper-face contacts, no residual primitive face in the starting domain |
| Unique selection of coefficient t² | Wrong if inferred | The same computation transports any specified binary coefficient; it does not select the canonical finite addition |
| Raw double-source integrals share the first-order convergence domain | Wrong | Source-pair degree 2t conflicts with AB degree −2/t |
| This conflict obstructs the renormalized second-order identity | Unclear; it does not follow | Source-pair subtraction restores local integrability with separated spectators on −1<t<0; remaining overlapping subtractions have to be matched |
| Double-source triple primitive divergence after all lower forests are consistently removed | Absent by power counting on −1<t<0 | Full degree 2t−2/t>0 and the spectator remainder bound below; existence and coefficients of the lower-face extension are still needed |
| All-order algebra, arbitrary-test products, global quantum sectors | Unclear | Outside this calculation |

There is no mathematical error requiring replacement of X5–X17 within their explicitly stated finite-order scope. The two “wrong” items reject potential stronger inferences, not sentences asserted in the audited report. No replacement of book text is proposed.

## 1. Independent first-order reconstruction

Use the vectors and covariance in the coordinator plan. Direct bilinear multiplication gives `a0²=a2²=−1`, `a0.a2=1−t`, `a0.d=a2.d=1`, and `d²=−rho`, with `rho=2/t`. Derivatives contract by

`(a0.partial X)(z) exp(q.X(w)) = −(a0.q)/(z−w) exp(q.X(w)) + regular`.

Consequently A at a contributes `+1/(z−a)` to B's derivative polynomial, and Y2 contributes `(t−1)/(z−w)`. The exponential powers are `(a−z)^−rho (a−w)/[(z−w)(ubar−wbar)]`. Multiplying the additive cocycles of A, B, barB, Y2 gives −1; the negative action insertion changes this to `+lambda2/pi`. This independently reconstructs X5. For `A B barK2` and `A K2 barB` the cocycles are respectively +1 and −1. In the latter, contracting the antichiral derivative at u against the source at z gives `−(t−1)/(zbar−ubar)`. Multiplication by the tail coefficient `−t lambda2` reconstructs X6 and X7 including their signs.

Let beta=t−1, set u=0, hold a separated, and put `V(w)=(a−w)^(1−rho)`. Expanding the chiral factor at z=w yields double coefficient `beta V` and simple coefficient `(beta rho−1)(a−w)^−rho = −V'`. The full field jet is fixed by

`a0+beta d=−(d+a2)`.

Thus the simple coefficient is the negative derivative of the actual normal product, including the A contraction. It is not legitimate to treat V as constant here.

The contour jump of this pole pair is `−2pi i[−V' Theta+beta V partial_w Theta]`. A circle of radius epsilon removed from the strip gives the horizontal integral `2 sqrt(epsilon²−y²)/epsilon²` for `|y|<epsilon`, and zero otherwise. Its y integral is pi. The elementary integrals `int dx/(x−iy)=i pi sign(y)` and `int dx/(x−iy)²=0` give both edge moments `J1=Jw=−pi`, for unequal heights as well. The first derivative term of the action kernel is therefore

`−2i lambda2 beta[−V' pi+beta V'(−pi)] = 2pi i lambda2 t beta V'`.

For the first tail, the simple-pole term is `t lambda2 V'/z`, whose contour jump is `−2pi i t lambda2 V'`. For the second tail it is `−t beta lambda2 V' z/zbar²`, whose jump is `−2pi i t beta lambda2 V'`. Their sum is `−2pi i t² lambda2 V'`. After restoring the A/endpoint cocycle, the action and tail coefficients are consequently `−t(t−1)` and `t²`, and the raw sum is t. The prescribed addition `t(t−1)` changes this sum to t². This calculation would equally transport another explicitly selected binary normalization.

## 2. Proper faces and first-order uniformity

The independent subset enumeration reproduces every degree in X12–X13. A cluster containing currents has one common real translation; a source-only cluster has two. Counting relative dimensions and the leading nonzero derivative monomial gives:

| Face | Relative degree |
|---|---:|
| AB | −rho |
| A2 | 3 |
| A barB | 1 |
| B2 or barB2 | 0, extended as Cauchy poles |
| B barB2 | −1, explicitly calculated first-jet contact |
| AB2 | 2−rho |
| A barB2 | 2 |
| AB barB, source separated | 1−rho |
| AB barB2 | 1−rho |

The scalar rational identity in X12 follows by multiplying out `(p−r)(p+beta r)/(rp²)`. Thus the primitive four-relative-coordinate kernel has degree 1−rho. Pinning to an endpoint adds one relative coordinate. The tails yield the same interior degree. The mixed face's constant term cancels under the contour jump; subtracting its complete first jet leaves a term bounded locally by `epsilon times the second derivative of A U`. Since `A U ~ R^(1−rho)`, its remaining spectator integral is

`epsilon int_(M epsilon)^R0 R^(−1−Re rho) dR`,

which is at most `epsilon R0^(−Re rho)/(−Re rho)` when Re rho<0. The complementary corner is bounded by `epsilon^(1−Re rho)`. The Cauchy proper faces are distributions, not absolutely integrable functions; this argument presupposes their specified extension and includes the derivative seminorms that extension requires. It establishes the claimed restricted extension, not a Hilbert-space operator norm or convergence with arbitrary coincident spectators.

For the endpoint, the ratio `epsilon(A,U)/epsilon(U,A)=−s^−1`, while the oscillator power is `alpha=1−rho`. The positive-imaginary-separation boundary values multiply this ratio by `exp(−i pi alpha)` at the left endpoint and `exp(+i pi alpha)` at the right endpoint. The resulting A U / U A ratios are `s^−2` and 1. This gives X16 and, after the two actually computed mixed crossings, X15–X17. Compact support of f inside I and g=1 nearby are essential to this conclusion.

## 3. Fresh double-source kernel

Write D=barB, put Y0 at v and Y2 at w, and use the reference ordering `A(a) B(z) D(u) Y0(v) Y2(w)`. Its cocycle is +1. The cross term in the exponential of the action has coefficient `lambda0 lambda2/pi²`: the two source orders cancel the factorial 2. The bare-current contribution is

\[
\frac{\lambda_0\lambda_2}{\pi^2}
\frac{(a-z)^{-\rho}(a-v)(a-w)|v-w|^{2t-2}}
{(z-v)(z-w)(\bar u-\bar v)(\bar u-\bar w)}
\left[C(z)+\frac1{z-a}+\frac1{z-v}+\frac{t-1}{z-w}\right]
\left[\bar C(u)+\frac1{\bar u-\bar v}+\frac{t-1}{\bar u-\bar w}\right]
\mathcal N.
\]

Here the normally ordered factor is

`N = :exp[−dX(a)+dX(z)+d Xbar(u)+a0(X(v)+Xbar(v))+a2(X(w)+Xbar(w))]:`.

This is only the bare double-action sector. It excludes composite-potential insertion, one-source dressings of the first tails, and the order-two current improvements. It cannot alone decide the defect.

## 4. Complete new-face classification and forests

For any subset, add its pair contraction exponents, subtract one for each derivative current that has a nonzero internal contraction, and add its relative integration dimension. This generates all 26 nontrivial subsets in `faces.json`. The genuinely new faces, containing both sources, are:

| Face | Degree |
|---|---:|
| 02 | 2t |
| A02 | 2t+4 |
| B02 | 2t−1 |
| D02 | 2t−1 |
| AB02 | 2t+2−rho |
| AD02 | 2t+2 |
| BD02 | 2t−3 |
| ABD02 | 2t−rho |

The first-order source faces also occur for source 0: `B0,D0` have degree zero, `BD0` degree −1. A vanishing Ward residue for source 0 would not remove these distributions before their derivative/double-pole cancellations are evaluated. In particular it must not simply be omitted from a forest list.

Raw source-pair integrability requires Re t>0; raw AB integrability requires Re(2/t)<0. Since `Re(1/t)=Re(t)/|t|²`, there is no such common domain, even at complex t. This is not a renormalized no-go. With other insertions separated,

`Y0(x+q)Y2(x)=|q|^(2t−2)[P+q C P+qbar Cbar P+q qbar C Cbar P+...]`.

The circular angular integral removes the linear terms. Its constant radial divergence is canceled by the specified kappa. The next circular scalar term has positive degree `2t+2` on −1<t<0. Intersecting current exclusions can destroy that elementary angular cancellation, so the associated forests must implement the same extension before this simplification is used.

All new A-bearing faces are positive on the open real interval −1<t<0. The singular source/current family still contains `02`, `B0,B2,D0,D2`, `BD0,BD2`, `B02,D02`, and `BD02`. A forest is a collection of these that is nested or disjoint. The following overlapping families require matching rather than simultaneous subtraction as independent faces:

- `02` with `B0,B2,D0,D2`, whose unions include `B02,D02`;
- `B0` with `B2`, and `D0` with `D2`;
- `BD0` with `BD2`;
- `B02` with `D02`;
- `BDi` with `02`, `B02`, or `D02` when neither is nested; their union is `BD02`;
- nested chains such as `02 < B02 < BD02`, and disjoint pairs such as `B0` with `D2`.

These classes are all obtained by the set intersection/union rule on the displayed family. After source fusion, P-current faces must also be retained: BP and DP have degree −1, and BDP has degree −3, before multiplication by the cutoff coefficient epsilon^(2t). They reproduce the required nested source-pair geometry; they are not optional new independent bare diagrams.

## 5. What the degree calculation does and does not settle

On −1<t<0 a full Taylor extension can remove source-pair jets through degree one, one-current/two-source jets through degree two, and the mixed BD02 jets through degree four. The respective first unretained degrees are all at least `2t+2>0`. This is a sufficient subtraction bound, not a claim that every listed local coefficient is nonzero. Symmetric circular source integration removes odd angular moments. A reflection-symmetric mixed contour difference can restrict surviving mixed jets to odd orders, suggesting precisely orders one and three for the two displayed endpoint subtractions. The specified contour forest and its finite matching must prove those cancellations, especially in intersecting exclusions; a parity slogan does not do so.

There is a useful quantitative conditional bound for the new A collision. The fully fused mixed field is F, with `A F ~ R^(2−rho)`. After mixed Taylor subtraction through degree four, the fifth-jet nested remainder is bounded by

\[
\epsilon^{2t+2}\int_{M\epsilon}^{R_0}
R^{-3-\rho}\,dR.
\]

For −1<t<0 one has rho<−2, so the integral stays bounded and the prefactor vanishes. In the complementary region R=O(epsilon), the bound is `epsilon^(2t−rho)`, also vanishing. Thus a consistently constructed lower-face extension has no additional primitive A-dependent divergence in this strip. This is stronger than merely observing a positive full-cluster degree: it controls the potentially nonuniform spectator remainder. It is conditional on the lower extension's distributional Taylor bounds and cannot calculate its coefficients.

The next discriminating calculation is therefore the complete BD02 local extension and its overlaps with `02`, `B02,D02`, and `BD0,BD2`, including the action contact, tails, improvements and finite binary matching. That determines whether the two displayed endpoint subtractions suffice or whether a residual local contact remains. The raw no-common-domain observation neither answers nor blocks this calculation.

## Integration

No book integration occurred. Verified: the first-order contact, proper faces, nested spectator estimate, endpoint phases and restricted identity. Modified scope for the next calculation: a raw common domain is absent, but a plausible renormalized interval −1<t<0 has positive A-bearing degrees and a quantitative spectator bound. Unresolved: the actual double-source lower-face coefficients and forest consistency, arbitrary-test contacts, all-order mechanisms and global quantum closure. No literature claim or theorem-label upgrade is made.
