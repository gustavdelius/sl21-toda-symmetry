# Double-Yukawa three-contour test: lower-face subtractions and an explicit nonuniform overlap

Date: 2026-10-24. Model: Codex (GPT-6).

Reviewed commit: `745586e6a6f7bca5ed91a1e9cf89d6ef5df387f0`.

Worker scope: the subtraction and finite-contact part of stage 1 of `briefs/referee_sl21_symmetry_coordinator_plan_20261024.md`. No subagents were used. No pre-existing file was edited, no git mutation was performed, and no render was run. The other stage-1 worker's new outputs were not read during this derivation. Input hashes are recorded in `referee_scratch/sl21_symmetry_coordinator_20261024/stage1_subtractions/input_hashes.json`; they were recorded during intake, after the initial reads, and checked again at completion.

## Outcome

**The full coefficient of the four-word defect remains unclear in this worker's calculation.** The lower-face subtraction chain, including both endpoint descendants and the finite binary matching dressing, is independently correct. However, its fixed-separation estimates cannot be used uniformly on an overlapping mixed/source-pair face. An exact coefficient of the actual Wick kernel exhibits the failure, including the additional A current. This is an obstruction to that estimate, not an obstruction to the symmetry or a computed nonzero four-word residual. A complete overlapping extension may cancel it.

| Item | Verdict | Evidence |
|---|---|---|
| Mixed action expansion factor and sign | correct | expansion of the exponential, section 1 |
| Composite insertion and improved B-tail coefficients | correct | direct fusion and circular angular projection, sections 1–2 |
| Both U endpoint subtractions | correct for separated spectators, -1<t<0 | radial integration of the full first angular jets |
| Dressing of an arbitrary finite binary coefficient alpha | correct for separated spectators | same integral, with coefficient retained symbolically |
| Uniformity of the scalar source-pair remainder near the mixed face | wrong | exact nonzero remainder (S8) |
| Primitive A B barB Y0 Y2 contact | unclear | the mixed/source-pair union has not been extended here |
| Mass-squared projected identity | unclear | no substitution of a binary identity or assumed associativity for that union calculation |

The actual result is stronger than merely observing opposite convergence conditions: (S8) identifies an explicit coefficient requiring treatment in the missing overlap calculation.

## 1. Reconstruct the coefficient list before taking any contour defect

Use the vectors, covariance and additive cocycle in the coordinator plan. Direct bilinear multiplication gives a0 squared = a2 squared = -1, a0.a2=1-t, d.ai=1 and d squared=-2/t. Define

\[
K_2=e^{(d+a_2)X+a_2\bar X},\qquad
L=C e^{(d+T)X+T\bar X},\qquad C=a_0\partial X.
\]

Bars exchange chirality, without complex conjugating momenta. Every product below is an ordered Wick product with its cocycle, not a commuting scalar multiplication.

The lambda0 lambda2 coefficient of the exponential weight is

\[
\frac12\left(\frac{\lambda_0}{\pi}\int Y_0+
\frac{\lambda_2}{\pi}\int Y_2\right)^2\Big|_{02}
=\frac{\lambda_0\lambda_2}{\pi^2}\int Y_0Y_2.
\tag{S1}
\]

The two labelled cross terms cancel the factor 1/2. The ordered source convention is Y0 followed by Y2, whose full cocycle is +1. With an A insertion placed as in each word, the required sectors before Wick signs are

\[
\begin{split}
&\frac{\lambda_0\lambda_2}{\pi^2}\int A J_B\bar J_B Y_0Y_2,\\
&\frac{t\lambda_0\lambda_2}{\pi}\int
(A J_B\bar K_2Y_0+A K_2\bar J_BY_0),\\
&\lambda_0\lambda_2\epsilon^{2t}
(A J_B\bar L+A L\bar J_B),\\
&\frac{\lambda_0\lambda_2\epsilon^{2t}}{\pi t}\int A J_B\bar J_B P.
\end{split}
\tag{S2}
\]

The last sign is positive since the action coefficient is kappa=-lambda0 lambda2 epsilon^(2t)/t and the Euclidean weight contributes -kappa/pi. The second line is the product of the first tail -t lambda2 K2 and the action insertion -lambda0 Y0/pi. There is no A tail at this coefficient. Products of two first B tails have coefficient lambda2 squared, not lambda0 lambda2.

As a check of the local signs, in the displayed reference order the six cocycle multipliers in (S2) are respectively +1, -1, +1, +1, +1, +1. In particular, the two one-tail terms have opposite Wick signs. This follows by multiplying the cocycle over every pair of labels; it is not a parity sign attributed to the even currents. For the four words the ordering must be retained until the contour deformation. These six signs do not replace its braiding calculation.

Every charge has the common factor (2 pi i)^(-1). Thus the triple product carries (2 pi i)^(-3), before any contour residues have reduced its dimension. All double-source domains must include the source-pair hole, the inner-current exclusions, and the moving-current Cauchy extension of X14. The finite binary matching must be inserted on the relevant proper faces as a separate term.

## 2. Derive the lower-face counterterms directly

With r measured from Y2 or K2 to Y0, the exponential contraction is

\[
Y_2(0)Y_0(r)=|r|^{2t-2}(P+\text{jets}),
\]

\[
K_2(0)Y_0(r)=r^{t-2}\bar r^{t-1}
\big(E+rCE+\bar r\bar CE+\cdots\big),
\tag{S3}
\]

where E=exp[(d+T)X+T barX]. In the second line the full K2/Y0 cocycle is -1; (S3) displays the oscillator and field polynomial before this sign. The scalar angular term vanishes. The term rCE survives, giving the radial integral pi(R^(2t)-epsilon^(2t))/t. The action insertion -lambda0/pi, the K2/Y0 sign -1, and the tail coefficient -t lambda2 give the divergent contribution -lambda0 lambda2 epsilon^(2t)L. Thus the improvement is **+lambda0 lambda2 epsilon^(2t)L**. It is the insertion subtraction, not the older contour-descent primitive proportional to a different derivative field.

The scalar Y0Y2 integral similarly has lower-limit term -pi epsilon^(2t)P/t. Its coefficient in (S1) is cancelled by the last line of (S2).

For U=exp[(d+a2)(X+barX)] and F=exp[(d+T)(X+barX)], direct contraction instead gives

\[
U(0)Y_0(r)=|r|^{2t-4}
\left[F+rCF+\bar r\bar CF+r\bar r C\bar C F+
\text{higher jets}\right].
\tag{S4}
\]

The U/Y0 cocycle is +1. Circular averaging retains the displayed scalar and r bar r terms. The minus action insertion consequently yields the cutoff representative

\[
U_{\min}=\left[U-\frac{\lambda_0}{t-1}\epsilon^{2t-2}F
-\frac{\lambda_0}{t}\epsilon^{2t}C\bar C F\right]_{S_\epsilon}.
\tag{S5}
\]

The first omitted circular term is r squared bar r squared. Its lower-limit contribution is O(epsilon^(2t+2)) with **all other insertions held separated**. This is why -1<t<0 suffices for this particular local subtraction. Insertion of A contracts with the C factors too: those full Wick products must be retained, rather than replacing C barC by a spectator-independent constant.

If a binary endpoint coefficient is alpha lambda2 U, its renormalized dressing is alpha lambda2 times (S5). In particular the prescribed finite addition alpha=t(t-1) requires

\[
\delta B_{BB}^{\rm ren}=
\left[t(t-1)\lambda_2U
-t\lambda_0\lambda_2\epsilon^{2t-2}F
-(t-1)\lambda_0\lambda_2\epsilon^{2t}C\bar C F\right]_{S_\epsilon}.
\tag{S6}
\]

The action dressing -alpha lambda0 lambda2 integral(UY0)/pi is part of (S6), not optional. Adding this to the raw coefficient alpha=t gives the canonical coefficient t squared and exactly the two endpoint terms in the coordinator plan. The lower-face chain contains no inconsistent factor of two or unmatched sign. This fact alone does not evaluate the four-word defect.

## 3. An exact nonuniform coefficient on the intersecting face

Here is a calculation using the actual bare-current double-source kernel, without assuming any binary Ward identity. Set w2=0 and w0=r. Let z denote the B coordinate and v the antichiral coordinate of barB. Set normal-ordered derivative backgrounds to zero to isolate a scalar coefficient. The B contractions with the sources give

\[
f(z;r)=\frac{(z-r)^{-1}+(t-1)z^{-1}}{(z-r)z}
=z^{-3}\sum_{n\ge0}(n+t)(r/z)^n.
\tag{S7}
\]

Thus the angular average of f(z;r)f(v;bar r) is
(zv)^(-3) sum over n of (n+t)^2 (|r|^2/(zv))^n, for |r|<min(|z|,|v|). The n=0 term is exactly the local P source-pair coefficient. After subtracting it, the omitted disk 0<|r|<epsilon has the **exact** value

\[
\boxed{\mathcal E(\epsilon;z,v)=
\frac{\pi\epsilon^{2t}}{(zv)^3}
\frac{q[(1+t)-tq]}{(1-q)^2},\qquad
q=\frac{\epsilon^2}{zv}.}
\tag{S8}
\]

Indeed each n>=1 angular mode integrates to
pi(n+t)epsilon^(2t+2n)/(zv)^(3+n); summing the geometric derivative gives (S8). At fixed z,v it starts at pi(1+t)epsilon^(2t+2)/(zv)^4, agreeing with the familiar local remainder. At z=epsilon Z, v=epsilon V it is instead

\[
\mathcal E=\pi\epsilon^{2t-6}(ZV)^{-3}
\frac{(ZV)^{-1}[(1+t)-t/(ZV)]}{[1-(ZV)^{-1}]^2}.
\tag{S9}
\]

The coefficient is generically nonzero, and the exponent is negative throughout -1<t<0.

Including A at a separated point a does not remove it. The exact holomorphic scalar factor, apart from (a-z)^(-rho), becomes

\[
h(z,r;a)=\frac{a(a-r)}{(z-r)z}
\left[\frac1{z-a}+\frac1{z-r}+\frac{t-1}{z}\right].
\]

For every n>=1 the coefficient of r^n is

\[
[r^n]h=\frac{a(a-z)(n+t)}{z^{n+3}}.
\tag{S10}
\]

The cancellation uses the actual derivative contraction 1/(z-a); dropping that contraction would give a different answer. The exact A-dressed omitted disk is therefore **a(a-z)^(1-rho) times (S8)**. Its factor is smooth and nonzero with A held apart. The n=0 term has the corresponding A-dressed P coefficient and is precisely the one removed in deriving this statement.

To keep positive ordered contour heights, take c=sqrt(5), M=12, z=epsilon(M+2ic), v=epsilon M and a separated. These are the heights of the BAbarB word in X14 with delta=c epsilon; A can have its prescribed height and a fixed different real coordinate. Expressing its local kernel in the displayed reference order introduces only the corresponding branch/cocycle factor. The source w2=0 has distance M epsilon from the inner current, so its disk exclusion is respected. Since |z| and |v| exceed epsilon, the series calculation is valid. Choosing M large merely makes the coefficient small; at fixed M it does not change the power in (S9).

At t=-sqrt(2)/3, the absolute values of (S8) are:

| epsilon | absolute omitted disk |
|---|---:|
| 0.08 | 0.125073043623360 |
| 0.04 | 15.3871254338753 |
| 0.02 | 1893.00285863989 |
| 0.01 | 232886.892241074 |

After multiplying by epsilon^(6-2t), every row gives
4.14716968352841e-10 - 3.00207863473052e-9 i.
Independent 2048-point angular quadrature agrees with the angular mode sum to relative error below 1e-30 at these choices. This confirms the exact coefficient, not its cancellation or survival after the remaining integrations.

This counterexample is deliberately local. It is **not** an evaluated four-word commutator, and none of these diverging numbers is a residual algebra anomaly. The tail, improvement and endpoint sectors can cancel contributions with precisely such scaling. Its implication is that the estimate “the source-pair remainder is O(epsilon^(2t+2))” is not uniform on the needed domain.

## 4. Forest structure and the remaining mathematical task

The source-pair cluster S={Y0,Y2} and mixed contact cluster M={B,barB,Y2} overlap without inclusion. Their union W={B,barB,Y0,Y2} is a new face. One cannot replace both independently by their fixed-spectator limits and multiply the resulting formulas. A sector treatment must include:

1. S much smaller than its separation from the currents: the P subtraction and its Taylor remainder;
2. M much smaller than its separation from Y0: the U contact, its finite matching and both UY0 subtractions;
3. one current with both sources: the L or barL insertion subtraction;
4. the comparable-scale W face, together with intersections of its boundaries with sectors 1–3.

Equation (S8) is the transition term that an unproved uniform estimate discards. Its leading mode is nonzero; simply keeping a few more source-pair Taylor terms does not supply a uniform bound. At z,v proportional to epsilon every retained finite truncation leaves further nonzero modes with the same overall epsilon^(2t-6) scaling. This does not mean that an infinite counterterm tower is required: a sector extension of W is the appropriate alternative to a uniform pairwise Taylor expansion.

Adding A to the total cluster improves its superficial homogeneity. The scalar kernel has degree 2t-6-rho; seven real measures give 2t+1-rho when the whole cluster shrinks to a fixed endpoint, and removing a common translation along the contour gives the interior degree 2t-rho. Both are positive for -1<t<0. These powers concern the **whole** cluster; they do not cure (S8) at the proper W face with A apart. A conclusion of no primitive A contact would require a uniformly controlled extension of that proper face first.

Two continuation routes were examined at the level needed to identify their missing step:

* The positive-t binary continuation and negative-rho triple continuation have no common one-parameter open domain of their stated sufficient inequalities. Re t>1 implies Re(2/t)>0, whereas the AB estimate uses Re rho<0. Meromorphic continuation of two separately defined products does not by itself define their iterate.
* One may introduce independent analytic exponents, or replace A by a continuously charged vertex A_h=e^(-h dX), to produce a joint convergence region. This is an extension construction to be specified and tested. It changes the current/source derivative identities or adds branch-cut terms. They cannot be dropped before taking a potentially singular continuation to h=1. No Ward identity in that enlarged family has been established here.

The smallest missing input for this worker's stage-1 conclusion is an explicit extension of W, including the four sectors above, whose action on the A-dressed tests (S10) is controlled uniformly and whose finite parts agree with (S5)–(S6). That would allow the positive total degree to be used to decide whether a new primitive triple contact exists. Alternatively, a direct summation of the four regulated words including all of (S2) and (S6), with a nonzero residual and its cutoff scaling, would settle it. This report supplies neither result and does not relabel their absence as a symmetry no-go.

Finite scheme changes remain meaningful. Changing alpha in the binary U contact changes (S5) and its dressing together. An additional finite W-supported distribution is a separate product choice unless fixed by a common analytic extension or a normalization condition. The known binary normalization and the action choice kappa_R=0 do not, merely by being stated, prove such a common extension. Conversely an analytic continuation defined in a valid initial domain would forbid an arbitrary extra term added only after continuation.

## 5. Evidence and suggested wording

Fresh evidence is `referee_scratch/sl21_symmetry_coordinator_20261024/stage1_subtractions/check_subtractions.py` and `check_subtractions.out`. It includes independent vector contractions, exact Taylor coefficients, the complete radial subtraction coefficients, ordered cocycles, angular quadrature and nonuniform scaling diagnostics. The final run has 39 PASS lines and no FAIL lines. This count is not evidence of a charge identity; the report states precisely which coefficient was evaluated.

No correction to an existing asserted mass-squared triple identity is proposed, because the inherited triple report explicitly leaves that calculation open. The following is suitable wording immediately after the prospective extension described at `briefs/referee_sl21_triple_extension_report.md`, section 8 (line 463 in the reviewed input):

> At double-Yukawa order, the source-pair and mixed-current faces overlap. The local source-pair remainder is small only with the currents held separated; an explicit A-dressed coefficient scales as epsilon^(2t-6) when their separations scale with epsilon. All lower-face subtractions and finite matching must therefore be combined in a common extension of this union before the positive total triple degree can exclude a primitive contact. The existing first-order triple estimate does not perform this step.

## Integration

No integration was performed. Verified independently: action combinatorics, the current insertion improvement, both endpoint subtraction coefficients, the dressing of a general binary finite contact, and the exact nonuniform overlap coefficient. Rejected: a uniform fixed-separation source-pair remainder estimate in the overlapping region. Unresolved: the complete W extension, its combination with all four contour words, any primitive triple term, and the mass-squared projected identity. No algebra or theorem label is recommended.
