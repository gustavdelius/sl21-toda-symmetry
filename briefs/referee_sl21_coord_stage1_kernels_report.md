# Mass-squared triple kernels and an independent source-exponent extension

Date: 2026-10-24. Model: Codex (GPT-6). Reviewed commit: `745586e6a6f7bca5ed91a1e9cf89d6ef5df387f0`.

This bounded worker calculation preserves the completed stage-0 pass. New evidence is in `referee_scratch/sl21_symmetry_coordinator_20261024/stage1_kernels/`: `check_h.py`, `check_h.out`, `channel_powers.json`, and `protected_inputs.json`. The scripts were written independently and import no earlier research scripts. No other worker's new stage-1 report was read. No subagents, pre-existing-file edits, book edits, git mutations, worktrees, commits, pushes or renders were used.

## Outcome and scope

An independent source-pair exponent h gives an explicit common starting domain for the mass-squared triple calculation: `Re h>2`, `Re rho<0`, with `rho=2/t`. It also exhibits the regulator's Ward violation explicitly. That violation is proportional to `h−t`, provided the action and tail counterterms are changed with h. The ordinary local Mellin divisors of the double-source family do not contain the entire diagonal h=t. Consequently this regulator violation vanishes on the generic diagonal after meromorphic continuation.

This supplies a formal distributional extension of the projected mass-squared identity, compatible with the continued lower mixed product, **when the same entire h-family is continued before restricting h=t**. It does not assert that independently taking finite parts of differently ordered collision faces gives the same answer. Nor does it prove that the naive negative-t cutoff limit, using only separated-face subtractions, realizes this extension. The common h-family fixes the overlapping faces as part of the product definition.

| Item | Verdict | Evidence and limit |
|---|---|---|
| All channels at lambda0 lambda2 in A B barB | Correct | Explicit ordered kernels below, including the two improvements and P insertion |
| Keep the original tail-improvement coefficient unchanged away from h=t | Wrong | The coefficient must be t/h; otherwise a local Ward defect remains |
| h-deformed Ward violation is divisible by h−t | Correct | Actual double/simple-pole calculation and cutoff polynomial |
| No new pole divisor contains h=t identically | Correct for this local finite-order family | All new exponents are affine in h and rho, with nonzero h coefficient; fixed first-order Cauchy faces are treated separately |
| This construction agrees with the old lower mixed continuation | Correct as a specified analytic prescription | Agreement on the positive-t diagonal and uniqueness of the one-variable continuation |
| All independently chosen forest finite parts automatically agree | Wrong | Agreement follows only when they are restrictions of this one family; arbitrary lower-diagonal additions are excluded by that definition |
| Projected triple identity through lambda0 lambda2 | Correct in the common-family prescription described here | Stokes reduction with the explicit regulator variation, followed by generic diagonal restriction |
| General associative charge algebra or actual quantum representation | Unclear | Other current words, arbitrary tests and global sectors have not been constructed |

## 1. Definition of the two-parameter family

Retain all momenta, derivative contractions and cocycles of the plan. Replace only the unlike-source pair factor

`|v−w|^(2t−2)` by `|v−w|^(2h−2)`.

Equivalently multiply it by `R_h(v,w)=|v−w|^(2(h−t))`, using the ordinary real radial logarithm. This is an analytic regulator, not a changed free-field action or a claim about a new physical theory. When a current collapses onto a source, use the pullback of this same R_h; in particular the tail/source and endpoint/source factors must be deformed as well.

Use exactly the X14 heights and test functions: barB on I=[L,R], f smooth and supported strictly inside I, g=1 near I, delta=c epsilon with c>2. Moving-current poles are Cauchy-extended first. Keep the inner source circles and unlike-source exclusion at the same epsilon, fixed compact infrared cutoffs, and separated additional spectators. All four words use the same regulator family and the same first-order binary matching `t(t−1)lambda2 U`. Take epsilon to zero first on `Re h>2, Re rho<0`; continue this entire family meromorphically; finally restrict h=t at generic nonresonant t. One may start with real negative t and real h>2. No simultaneous raw convergence of the physical h=t integrals is being assumed.

The local counterterms obtained below are

\[
\kappa_h=-\frac{\lambda_0\lambda_2}{h}\epsilon^{2h},\qquad
H_{B,h}=-t\lambda_2K_2+
\frac{t}{h}\lambda_0\lambda_2\epsilon^{2h}L,
\]

with the corresponding barred expression. Here

`K2=exp[(d+a2)X+a2 Xbar]`, `L=C exp[kX+T Xbar]`, and their barred versions exchange chiralities. The endpoints use

\[
[U]_h=U-\frac{\lambda_0}{h-1}\epsilon^{2h-2}F
-\frac{\lambda_0}{h}\epsilon^{2h}:C\bar C F:.
\tag{H1}
\]

These are local cutoff representatives on separated faces. At intersecting faces, the definition is the common meromorphic family, rather than independent use of the displayed radial representatives.

## 2. Actual ordered mass-squared kernels

Let beta=t−1, put A at a, B at z, barB at u, Y0 at v and Y2 at w. Every polynomial below is inside the same normal product as the exponential factors. For each displayed channel the exponential has the sum of its listed momenta: A carries (−d,0), B (d,0), barB (0,d), Yi (ai,ai), K2 (d+a2,a2), barK2 (a2,d+a2), L (k,T), barL (T,k), P (T,T). Derivative prefactors C and Cbar remain at their indicated insertion points. Thus the notation specifies the full kernel, including uncontracted fields.

The double-action channel has coefficient `lambda0 lambda2/pi²`, cocycle +1, scalar factor

\[
\frac{(a-z)^{-\rho}(a-v)(a-w)|v-w|^{2h-2}}
{(z-v)(z-w)(\bar u-\bar v)(\bar u-\bar w)},
\]

and derivative polynomial

\[
\left[C(z)+\frac1{z-a}+\frac1{z-v}+\frac\beta{z-w}\right]
\left[\bar C(u)+\frac1{\bar u-\bar v}+\frac\beta{\bar u-\bar w}\right].
\tag{H2}
\]

The action cross term has no remaining factor 1/2: the two different source orders cancel it.

The `A B barK2 Y0` tail/source channel has coefficient **`−t lambda0 lambda2/pi` after its cocycle**, scalar factor

\[
\frac{(a-z)^{-\rho}(a-u)(a-v)
(u-v)^{h-1}(\bar u-\bar v)^{h-2}}
{(z-u)(z-v)},
\]

and polynomial `C(z)+1/(z−a)+beta/(z−u)+1/(z−v)`. The minus sign is substantial: its ordered cocycle is −1, whereas the tail times negative action insertion before the cocycle has coefficient +t.

The `A K2 barB Y0` channel has coefficient **`+t lambda0 lambda2/pi`**, scalar factor

\[
\frac{(a-z)^{1-\rho}(a-v)
(z-v)^{h-2}(\bar z-\bar v)^{h-1}}
{(\bar z-\bar u)(\bar u-\bar v)},
\]

and polynomial `Cbar(u)−beta/(zbar−ubar)+1/(ubar−vbar)`. Its ordered cocycle is +1. The h-dependent powers in both kernels are pullbacks of the original source-pair factor, not new independent choices.

The two improvement channels both have coefficient `(t/h)lambda0 lambda2 epsilon^(2h)` and cocycle +1:

\[
AB\bar L:\quad
\frac{(a-z)^{-\rho}(a-u)^2}{(z-u)^2}
\left[C(z)+\frac1{z-a}+\frac{t}{z-u}\right]\bar C(u),
\]

\[
AL\bar B:\quad
\frac{(a-z)^{2-\rho}}{(\bar z-\bar u)^2}
\left[C(z)+\frac1{z-a}\right]
\left[\bar C(u)-\frac{t}{\bar z-\bar u}\right].
\tag{H3}
\]

The P-action channel has coefficient `lambda0 lambda2 epsilon^(2h)/(pi h)`, cocycle +1, and kernel

\[
\frac{(a-z)^{-\rho}(a-v)^2}
{(z-v)^2(\bar u-\bar v)^2}
\left[C(z)+\frac1{z-a}+\frac{t}{z-v}\right]
\left[\bar C(u)+\frac{t}{\bar u-\bar v}\right].
\tag{H4}
\]

The sign is positive because this is `−kappa_h P/pi` in the Euclidean weight.

Finally, the binary finite matching is dressed as `t(t−1)lambda2[U]_h`, with its ordinary negative Y0 action insertion. Its source kernel is the actual `U(x)Y0(v)` normal product multiplied by `|x−v|^(2(h−t))`; after multiplying A the scalar factor is `(a−x)^(1−rho)(a−v)|x−v|^(2h−4)`. Its local endpoint subtractions are precisely H1 multiplied by t(t−1)lambda2. This term is present on each proper mixed face in all four words. It is not an adjustable primitive triple addition.

There is no lambda0 lambda2 product of two first-order B tails: both such tails carry lambda2. The A tail carries the independent lambda1 coupling. These are outside the coefficient investigated here.

## 3. The regulator-dependent Ward coefficient, calculated directly

The crucial calculation can be made before any mixed-current integration. For a current collapsing at z onto source i, let q=z−other source, and call the scalar collapsed factors F0 and F2. They are proportional to `q^(h−2) qbar^(h−1)` times their normal products. Expanding the actual H2 kernel at its two current poles gives double and simple coefficients

\[
D_0=F_0,\qquad S_0=\partial_zF_0-\frac{h-t}{q}F_0,
\]
\[
D_2=(t-1)F_2,\qquad
S_2=-\partial_zF_2+\frac{h-t}{q}F_2.
\tag{H5}
\]

For example the actual scalar simple coefficient at node 0 is `(t−2)/q`, whereas differentiating the deformed double coefficient gives `(h−2)/q`. Their difference is exactly `−(h−t)/q`. At node 2 the corresponding coefficients are `(2−t)/q` and `−(h−2)/q`, yielding the opposite sign. The field-vector identities giving these derivatives are `C+d.partial X=(a0+d).partial X` and `C+(t−1)d.partial X=−(d+a2).partial X`. This independently checks the full first jets, not only scalar homogeneity.

Put r=other−z. Including the derivative of the exclusion boundary, the h-dependent Ward remainder beyond the first tail contains

\[
\frac{\lambda_0\lambda_2}{\pi}
\left[B_2-B_0+(h-t)\int_{|r|>\epsilon}
\frac{F_0-F_2}{r}\,d^2r\right].
\tag{H6}
\]

Here `Bi=int (partial_z Theta) Fi` with the same collapsed factors and the same exclusion as in H5. The signs in H6 use r=other−z; replacing r by q reverses the displayed bulk difference.

For independent field-jet variables `c=a0.partial X`, `e=a2.partial X`, `c1=a0.partial²X`, `e1=a2.partial²X`, the second-jet difference and P residue are

\[
N=\tfrac12[c_1-e_1+(c-e)(c+e)]E,
\quad R_P=[c(c+e)-\tfrac t2(c_1+e_1+(c+e)^2)]E,
\]

where `E=exp[kX+T Xbar]`. Direct subtraction gives

`N−R_P/t = partial L`.

Angular projection of H6 gives `B2−B0=pi epsilon^(2h)N+O(epsilon^(2Re h+2))` and

\[
\frac1\pi\int_{\epsilon<|r|<R}\frac{F_0-F_2}{r}\,d^2r
=\frac{R^{2h}-\epsilon^{2h}}{h}N+\text{higher jets}.
\]

The coefficient of the cutoff term is therefore t/h, not 1. Combining the composite counterterm and the improvement in H3 leaves

\[
\mathcal W_h=(h-t)\lambda_0\lambda_2
\left[\frac1\pi\int_{|r|>\epsilon}\frac{F_0-F_2}{r}\,d^2r
+\frac{\epsilon^{2h}}hN\right]
+O(\epsilon^{2\operatorname{Re}h+2}).
\tag{H7}
\]

This is an explicit renormalized bilocal Ward violation proportional to h−t. The antiholomorphic computation is its chirality exchange. As a numerical check independent of previous code, setting the selected smooth radial jet to `exp(−r²)` gives the continued bracket Gamma(h). Direct subtraction and quadrature agree at `t=−sqrt(2)/2` and `t=−1/sqrt(3)`. For the first value the anomaly is −0.04248355871 at h−t=0.01 and −0.004329305694 at h−t=0.001, consistent with a regular coefficient times h−t. This checks one test function; the pole argument below handles the formal distribution family.

## 4. Common domain, continuation, and overlap compatibility

The full family begins on the open set `Re h>2, Re rho<0`. With all first-order Cauchy contacts treated, every new collision degree is positive there. The source-containing degrees are obtained from the stage-0 table by replacing t in the source-pair powers by h:

`02:2h; B02,D02:2h−1; BD02:2h−3; A02:2h+4; AB02:2h+2−rho; AD02:2h+2; ABD02:2h−rho`.

Tail/source channels have the same relevant full degree; the explicit improvement channels are softer. Counterterm contributions vanish on this open set as epsilon tends to zero. The finite first-order mixed contact and its spectator dressing survive and are fixed before this limit. Nested source-pair, one-current/source-pair, and mixed clusters are all integrable after those fixed faces; no incompatible diagonal extension is being multiplied in this region.

Now apply the contour deformation to H2–H4 on this common open set. The actual simple/double-pole reduction is H5. Its derivative part cancels the tail/source channels. On a mixed one-source face it gives the first-order contact already audited, including both strip edges. The finite binary matching supplies its specified extra t(t−1) term. The remaining nonlocal terms are exactly insertions of H7 and its barred version; derivatives of R_h on the finite matching's spectator dressing also have an explicit factor h−t. All new simultaneous boundaries vanish in this open set by the displayed degrees. The A/source crossing has a linear zero and R_h is independent of a, so it produces no additional Ward violation. This gives a distributional identity of the form

\[
\mathcal D_h-(s^{-1}-s)A_f\frac\nu{2\pi i}
\left.t^2\lambda_2[U(L)]_h\right|_{\lambda_0\lambda_2}
=(h-t)\mathcal Q_h.
\tag{H8}
\]

Here `D_h` is the actual four-word contour defect, expanded with H2–H4 and the finite matching. `Q_h` is the Stokes remainder formed from H7, its barred version, and the explicit derivative of R_h on the matching term. This is a definition by calculated regulated integrals, not an arbitrary residual divided by h−t: the factor comes from H5 before integration. The right endpoint cancels using the same separated endpoint phases as at first order. The displayed coefficient includes the ordinary negative Y0 source dressing of U; it is not just its two local subtraction terms.

For completeness, the continuation argument requires more than the absence of an A primitive divergence. Partition a compact configuration patch into sectors ordered by cluster radii, with the fixed Cauchy faces already interpreted as distributions. Taylor-expand smooth field and test jets to any desired finite order on each sector. Each remaining radial integral has the form `int_0^1 r^(a h+b rho+n−1) dr` times a smooth angular/distributional coefficient. Subtracting sufficiently many Taylor jets yields a holomorphic remainder on any compact parameter set. The explicit Taylor integrals give the meromorphic extension. Nested sectors implement the source-pair, B02/D02, and BD02 forests simultaneously. Two partitions give the same continuation because they agree with the same convergent family on the initial open set. Thus overlap consistency is supplied by one meromorphic family; it is not inferred from separately continued collision formulas.

Every new source-pair divisor has a nonzero h coefficient. On rho=2/t, restriction to h=t gives a divisor `a t+2b/t+n`, which is not identically zero for a nonzero. Integer shifts caused by the extra 1/r in H7 do not change that statement. The only parameter-independent singular faces are the prescribed single-source Cauchy and first-order mixed faces, already fixed and analytic in h with the spectator source separated. Counterterm denominators h and h−1 also do not contain the entire diagonal. Therefore neither Q_h nor the full family has a pole divisor identically h=t. Isolated resonant t values must still be excluded. Generic irrational t alone is not a sufficient exclusion criterion for every algebraic resonance; the actual Mellin divisors must be avoided.

At such generic t the diagonal restriction of `(h−t)Q_h` is zero. H8 therefore gives the claimed projected identity through lambda0 lambda2 in this prescription. In the real strip −1<t<0, the new A-containing primitive degrees are positive, and the stage-0 estimate `epsilon^(2t+2) int R^(−3−rho)dR` independently controls the spectator corner after lower subtraction.

There is a separate check that this lower extension agrees with the pre-existing binary continuation. Remove A from the family. Its convergent domain has Re h>2 without the AB restriction. On the open positive-real diagonal h=t>2, its kernels, masks, Cauchy faces and binary matching are exactly the original binary calculation. The diagonal restriction of this meromorphic family is therefore the same one-variable continuation as that binary prescription, away from its poles. This comparison would fail if one changed a lower finite addition during continuation. No such addition is made here.

## 5. Exact limitations and possible finite changes

The construction fixes an analytic extension of the full regulated channels. It does not prove equality with every direct cutoff subtraction implemented in the negative-t strip. In particular, an independently assigned finite term on BD02 or an arbitrary derivative-of-test-function contact would define another product unless it were shown to be the restriction of the same family. The absence of an A primitive pole by itself never removes that freedom.

Changing the first-order binary matching by `delta c lambda2 U` changes its entire h-dressed family by the same amount. The projected right side then changes by `(s^−1−s)A_f nu delta c lambda2[U(L)]_h/(2pi i)`. Thus this triple test propagates the chosen binary coefficient and does not select it physically. One cannot add a finite term only after restriction while still claiming the unique continuation of the original common family.

No arbitrary-test density product, four-current parenthesization, Serre compatibility, all-order conservation, infinite-volume endpoint limit, or global quantum representation follows. The finite-order result is a formal continued contour identity in an explicit product prescription. The coordinator should subject the common-family Stokes and continuation argument to an independent adversarial check before using it as a premise of a higher-order claim.

## Integration

No book integration occurred. Independently verified: all lambda0 lambda2 triple channels, their ordering signs, the h-dependent pole residues, the t/h improvement, the composite and endpoint subtraction denominators, and the absence of a pole divisor containing the physical diagonal. Established within the stated formal prescription: a common analytic extension of the projected mass-squared identity and its agreement with the old binary continuation. Rejected: keeping counterterms unchanged off the diagonal, and inferring overlap compatibility from positive A degrees alone. Unresolved: equality with alternative cutoff products, arbitrary-test contacts, higher products and the global algebra.
