# Root-system sl(2|1) Toda: the first three-contour extension and quantum endpoint normalization

Date: 2026-10-24. Model: Codex (GPT-6).

Commit reviewed: `745586e6a6f7bca5ed91a1e9cf89d6ef5df387f0`.

Task: continue the boundary and quantum-vacuum investigation, independently auditing the latest calculation before extending it. `AGENTS.md` was read. Two independent workers were used under the user's authorization; neither delegated further. All work is in new scratch directories and new reports. Existing reports, book files, and check scripts were left unchanged. No git mutation, commit, push, or book render was performed.

## Outcome and verdicts

The first nontrivial projected three-contour test can now be evaluated from its actual Wick kernels. With the test functions, contour order, proper-face subtraction and analytic continuation specified below, it gives

\[
\boxed{
\left.[B_gA_f-s^{-1}A_fB_g,\bar B_I]\right|_{\lambda_2}
=(s^{-1}-s)A_f L_B^{(1)},\qquad
L_B^{(1)}=\frac{\nu}{2\pi i}t^2\lambda_2 U(L).
}
\tag{X1}
\]

Here the restriction means the contribution **linear in the node-2 Yukawa coupling** \(\lambda_2\), not second perturbative order. The calculation includes the bare-current action insertion, both current tails and the prescribed finite binary contact. It controls the simultaneous collision independently of an associative-algebra assumption. It remains a formal, prescription-dependent finite-segment identity. In particular, it does not determine arbitrary-test density products or establish compatibility with every earlier four-current construction.

The quantum-vacuum comparison also becomes more precise. A change of finite Gaussian subtraction must be carried through the action, OPE coefficients, currents and endpoint operators together. This gives an explicit normalization test for a proposed continuum completion. At fixed action the two candidate endpoint amplitudes change by the same current-pair factor. Their ratio is therefore invariant under this particular transport. The previous account of the vacuum branch path needs one correction: its logarithm windings describe compensated cycle/coupling transport, rather than the direct pullback of a translated cycle.

| Item | Verdict | Independent evidence and scope |
|---|---|---|
| Local three-current/one-source kernel | **Correct** | Direct contractions from the three momentum vectors, including derivative contractions and cocycles |
| Mixed face in the presence of the extra current | **Correct** | Its double and simple poles are \((t-1)V\) and \(-V'\); full first jets and strip edges are retained |
| Need for an additional primitive triple contact in the starting domain | **Wrong** for the specified extension | Positive collision degrees and a bound on the nested remainder after the proper face is removed |
| Projected identity (X1) and propagation of the finite binary matching | **Correct within the specified prescription** | Actual face evaluation, collision bounds and separated endpoint phases |
| Full associativity or a realized symmetry algebra | **Unclear** | Other products, higher orders and global charge domains are not constructed |
| Previously proposed reflection values | **Correct conditional formulas** | Source specialization, independent reflection-zero/identity-pole calculation and branch checks |
| Transporting the numerical vacuum values without transporting the Ward prescription | **Wrong** | Explicit finite Gaussian conversion factors |
| Interpretation of V21 as direct cycle pullback | **Wrong as written** | Exact change of integration variable; replacement supplied in section 7 |
| Actual full Toda endpoint limits in interface sectors | **Unclear** | Reference scalar fields and local ultraviolet matching do not construct those limits |

The two “wrong” entries concerning a primitive triple contact and incomplete normalization transport reject possible inferences, not claims made in the input reports. No book theorem upgrade is proposed.

## 1. Local definitions and independent reconstruction

Use chiral covariance \(X^i(z)X^j(w)\sim-\delta^{ij}\log(z-w)\), with

\[
a_0=(ib/\sqrt2,b/\sqrt2,i),\quad
a_2=(ib/\sqrt2,-b/\sqrt2,-i),\quad
t=b^2,\quad T=a_0+a_2,\quad d=-T/t,
\quad \rho=2/t,\quad \beta=t-1.
\]

The third component is the locally bosonized fermion current. Scalar products are bilinear, without complex conjugation. Multiplication of these vectors gives

\[
a_0^2=a_2^2=-1,\quad a_0a_2=1-t,\quad a_i d=1,\quad d^2=-\rho.
\tag{X2}
\]

The physical two-component Yukawa momenta are null; their fermion vertices supply the additional squared norm \(-1\). These statements are local OPE checks, not a global bosonization theorem.

The relevant currents, source and endpoint are

\[
J_A=:e^{-dX}:,\quad J_B=:C e^{dX}:,\quad C=a_0\partial X,
\quad H_B=-t\lambda_2:e^{(d+a_2)X+a_2\bar X}:,
\]
\[
Y_2=:e^{a_2(X+\bar X)}:,
\qquad U=:e^{(d+a_2)(X+\bar X)}:.
\tag{X3}
\]

The action insertion is \(-\lambda_2\int Y_2\,d^2w/\pi\). Charge integrals have normalization \((2\pi i)^{-1}\) and density \(J+H\). Set \(\nu=e^{-i\pi/t}\), \(s=\nu^2\). The additive cocycle is

\[
\varepsilon(n,m)=(-1)^{(n_L+n_R)(m_0+m_2)}
\exp[-i\pi(n_Lm_R-n_Rm_L)/t].
\tag{X4}
\]

Both candidate currents are even. No algebra is inferred from the root-system name or from their number.

For the reference ordering \(A(a)B(z)\bar B(u)\), followed by \(Y_2(w)\), the cocycle is \(-1\). Together with the action sign, the coefficient is \(+\lambda_2/\pi\). The exact bare-current kernel is

\[
\begin{split}
\frac{\lambda_2}{\pi}
\frac{(a-z)^{-\rho}(a-w)}{(z-w)(\bar u-\bar w)}
\left[C(z)+\frac1{z-a}+\frac\beta{z-w}\right]
\left[\bar C(u)+\frac\beta{\bar u-\bar w}\right]\\
\times:e^{-dX(a)+dX(z)+a_2X(w)+d\bar X(u)+a_2\bar X(w)}:.
\end{split}
\tag{X5}
\]

For example, differentiating the Gaussian contraction gives
\(C(z)e^{qX(w)}\sim-a_0q\,e^{qX(w)}/(z-w)\). Its contractions with \(A\) and \(Y_2\) are therefore \(+1/(z-a)\) and \(\beta/(z-w)\). This reconstructs the derivative polynomial directly.

The two tail kernels, with their insertion-order signs, are

\[
\begin{split}
-t\lambda_2\frac{(a-z)^{-\rho}(a-u)}{z-u}
\left[C(z)+\frac1{z-a}+\frac\beta{z-u}\right]
:e^{-dX(a)+dX(z)+a_2X(u)+(d+a_2)\bar X(u)}:,
\end{split}
\tag{X6}
\]

\[
t\lambda_2\frac{(a-z)^{1-\rho}}{\bar z-\bar u}
\left[\bar C(u)-\frac\beta{\bar z-\bar u}\right]
:e^{-dX(a)+(d+a_2)X(z)+a_2\bar X(z)+d\bar X(u)}:.
\tag{X7}
\]

The minus sign in the final bracket results from differentiating at \(u\), with the antichiral difference written as \(\bar z-\bar u\). There is no \(A\) tail at this coupling order.

## 2. Calculate the proper mixed face before taking the triple limit

Hold \(a\ne u\) fixed and put \(u=0\). At the level of the scalar Wick factor let \(V(w)=(a-w)^{1-\rho}\), with the branch inherited from (X5). Expanding the **actual** holomorphic factor at \(z=w\) gives

\[
\frac{\beta V(w)}{(z-w)^2}-\frac{V'(w)}{z-w}+O(1).
\tag{X8}
\]

This is not obtained by imposing a product rule on an assumed renormalized product. For a direct check, write the prefactor as \((a-z)^{-\rho}(a-w)\). Its derivative contributes \(\beta\rho\); the contraction with \(A\) contributes \(-1\). The identity \(\beta\rho=2-\rho\) gives the simple pole \(-V'\).

Uncontracted field jets also obey the needed identity:

\[
a_0+\beta d=-(d+a_2).
\tag{X9}
\]

Thus the simple coefficient with jets is \(-V'-V(d+a_2)\partial X\), the derivative of the complete dressed endpoint factor. The antichiral jets work in the same way. One cannot drop the derivative of the additional-current Wick factor.

Extend moving \(B\)/source poles as Cauchy distributions first. Let the inner \(\bar B\)/source exclusion be the circle \(|w|>\epsilon\), and let the two outer contour heights have magnitudes \(\delta_+,\delta_->\epsilon\). If \(\Theta\) denotes the intervening strip, the double pole gives a strip-edge derivative as well as the simple-pole term. The three finite integrals are

\[
I=\int\frac{\Theta\,d^2w}{\bar w^2}=\pi,\quad
J_1=\int\frac{\partial_w\Theta\,d^2w}{\bar w}=-\pi,\quad
J_w=\int\frac{w\partial_w\Theta\,d^2w}{\bar w^2}=-\pi.
\tag{X10}
\]

These constants can be computed without the earlier scripts. At fixed \(|y|<\epsilon\), removing the circle from the full horizontal integral of \(1/(x-iy)^2\) leaves \(2\sqrt{\epsilon^2-y^2}/\epsilon^2\); its \(y\) integral is \(\pi\). Outside that band the horizontal integral is zero. For the strip edges use
\(\partial_w\Theta=-i[\delta(y+\delta_-)-\delta(y-\delta_+)]/2\), and
\(\int dx/(x\pm i\delta)=\mp i\pi\). The remaining double-pole horizontal integral vanishes. Each edge contributes \(-\pi/2\) to either \(J\). Unequal heights give the same limits.

In the reference ordering (X5), the action contact is
\(+2\pi i\lambda_2t(t-1)V'(0)\), while (X6) and (X7) sum to
\(-2\pi i\lambda_2t^2V'(0)\). Their sum is
\(-2\pi i\lambda_2tV'(0)\). Restoring the common mixed-pair phase and the actual product with \(A\), the coefficients multiplying the dressed endpoint derivative are therefore

\[
\begin{array}{c|c}
\text{sector}&\text{coefficient of }\lambda_2 U\\\hline
\text{bare currents and action}&-t(t-1)\\
\text{two tails}&t^2\\
\text{raw circular prescription}&t\\
\text{specified finite addition}&t(t-1)\\
\text{canonically matched prescription}&t^2.
\end{array}
\tag{X11}
\]

The [independent kernel audit](referee_sl21_triple_kernel_report.md) also integrates the unexpanded additional-current kernels and strip edges. At irrational \(t=-\sqrt2\), the last two cutoffs yield extrapolation errors \(7.31\times10^{-5}\) for the action, \(7.19\times10^{-5}\) for the tails and \(1.21\times10^{-6}\) for their sum. This is a numerical check of the integrals, rather than a check of adding the entries in (X11). The root independently reconstructed (X5), (X8)–(X11); that numerical output is supplementary evidence.

This test propagates the selected finite normalization. It does **not** single out \(t^2\) as the unique physically admissible normalization: the raw \(t\) prescription passes its corresponding projected identity as well.

## 3. Simultaneous collision and nested remainders

It is insufficient to calculate the face with \(a\) separated and assume that its limit remains uniform. Put \(r=z-a\), \(p=z-w\). The scalar chiral factor reduces exactly to

\[
\frac{p-r}{p}\left(\frac1r+\frac\beta p\right)
=\frac1r+\frac{\beta-1}{p}-\frac{\beta r}{p^2}.
\tag{X12}
\]

Together with \(r^{-\rho}\) and the antichiral double pole, its homogeneity is \(-\rho-3\). There are four relative real integration coordinates at a simultaneous interior collision: two additional current positions and the two source coordinates. The integrated degree is therefore \(1-\rho\). At an endpoint there is one more small coordinate. The useful degrees are as follows; moving-current/source Cauchy poles are resolved as in section 2.

| Face | Integrated relative degree or treatment |
|---|---|
| \(AB\) | \(-\rho\) |
| \(A/Y_2\), \(A/\bar B\) | A linear zero, or no oscillator singularity |
| \(B/Y_2\), \(\bar B/Y_2\) | Prescribed Cauchy poles and their current tails |
| \(B\bar B/Y_2\) | The complete finite first-jet face (X11) |
| \(AB/Y_2\), with \(\bar B\) separated | \(2-\rho\) |
| \(A\bar B/Y_2\), with \(B\) separated | \(2\), after its proper Cauchy pole |
| \(AB\bar B\), source separated | \(1-\rho\) |
| Full action cluster, interior / endpoint | \(1-\rho\) / \(2-\rho\) |
| Full tail cluster, interior | \(1-\rho\) |
| \(A\) with the fused endpoint \(U\) | \(2-\rho\) |

After the proper mixed face, all additional degrees are positive in the common domain \(\operatorname{Re}\rho<0\). Negative real \(t\), including the vacuum comparison range \(-1<t<0\), lies in this domain.

More explicitly, take the mixed cluster radius \(\epsilon\) and the separation from \(A\) to be \(R\gg\epsilon\). The mixed kernel has integrated degree \(-1\), its constant odd term cancels, and its first jet is retained in (X11). The next remainder has an \(O(\epsilon)\) coefficient multiplying two derivatives of the additional-current factor. Since \(A U\sim R^{1-\rho}\), the remaining integration is bounded by

\[
\epsilon\int_{M\epsilon}^{R_0}R^{-1-\operatorname{Re}\rho}\,dR
=\frac{\epsilon}{-\operatorname{Re}\rho}
\left[R_0^{-\operatorname{Re}\rho}-(M\epsilon)^{-\operatorname{Re}\rho}\right]
=O(\epsilon),
\tag{X13}
\]

with fixed \(M>1\). In the complementary region \(R=O(\epsilon)\), the primitive estimate is \(O(\epsilon^{1-\operatorname{Re}\rho})\). Smooth field jets and test-function derivatives have no worse degree. The same reasoning applies to the tails. This supplies the missing uniformity check at the collision of the spectator with the mixed face.

No new primitive finite contact is needed in this starting domain. We choose its meromorphic continuation as the extension at generic other couplings, excluding poles and retaining the same lower-face finite additions. In this choice an extra primitive term cannot be freely added after continuation: it would no longer be the continuation of the starting family. Other finite prescriptions remain possible; nothing here proves their equivalence.

There is no like-chirality \(AA\) or \(BB\) pair in this particular test. Consequently the earlier difficulty of finding a common convergence domain for certain four-current tail integrals is not resolved by (X13). Compatibility with those products remains a separate question.

## 4. The four words and their boundary term

To make the test precise, take \(I=[L,R]\) for the inner \(\bar B\) integral. Let \(A_f\) have a smooth spatial test function supported strictly inside \((L,R)\). Let the outer \(B_g\) test function equal one on a neighborhood of \(I\). Use identical distant infrared cutoffs and separated spectator insertions in all words. Put \(\bar B\) at time zero, \(\delta=c\epsilon\), with fixed \(c>2\), and choose

\[
\begin{array}{c|rr}
\text{word}&\tau_B&\tau_A\\\hline
BA\bar B&2\delta&\delta\\
AB\bar B&\delta&2\delta\\
\bar BBA&-\delta&-2\delta\\
\bar BAB&-2\delta&-\delta.
\end{array}
\tag{X14}
\]

Use the inner source circle and the Cauchy extension of section 2 in every word. In particular, do not introduce a new circular source hole around the moving current before extending its pole. Insert the same binary finite matching on each proper \(B\bar B/Y_2\) face. Additional exclusions of new \(A\) diagonals can be removed using section 3. The compact support of \(f\), and \(g=1\) near \(I\), avoid unrelated endpoint and outer-test contact terms.

Moving \(\bar B\) past each ordered pair without interchanging \(A\) and \(B\) has no off-diagonal \(A\bar B\) contact. The actual mixed contact calculated above contributes

\[
BA\bar B-\bar BBA\longmapsto s^{-1}f_BA_f,
\qquad
AB\bar B-\bar BAB\longmapsto A_ff_B,
\quad f_B=R_B-L_B.
\tag{X15}
\]

The factor \(s^{-1}\) in the first line is from the off-diagonal exchange. Section 3 rules out a leftover simultaneous term for these contours in the starting domain. Thus (X15) follows from the regulated kernels; it is not merely an abstract expansion in a product whose associativity has been assumed.

The endpoint exchange is also independently checkable. The labels of \(A\) and \(U\) are \((-1,0,0,0)\) and \((1,1,0,1)\), so the ratio of their cocycles is \(-s^{-1}\). Their chiral power is \(\alpha_U=1-\rho\). With the first insertion assigned positive imaginary separation, the oscillator exchange ratio is \(e^{-i\pi\alpha_U}\) when \(a>x\), and \(e^{+i\pi\alpha_U}\) when \(a<x\). Hence

\[
U(L)A_f=s^2 A_fU(L),\qquad U(R)A_f=A_fU(R).
\tag{X16}
\]

Combining the four evaluated words gives
\(s^{-1}(f_BA_f-A_ff_B)\). The right endpoint cancels and (X16) gives (X1). Applying the same calculation to the finite addition alone yields

\[
\boxed{\delta[D,\bar B_I]\big|_{\lambda_2}
=(s^{-1}-s)A_f\frac{\nu}{2\pi i}t(t-1)\lambda_2U(L),
\qquad D=B_gA_f-s^{-1}A_fB_g.}
\tag{X17}
\]

This completes the lowest-order projected test requested by the preceding product report. It leaves the \(\bar A\) companion, the double-Yukawa dressing of (X1), arbitrary-test density products and the remaining quartic overlaps uncalculated in a common prescription.

## 5. Recheck the quantum reference value before transporting it

For generic irrational \(-1<t<0\), write \(t=-2\zeta^2\), \(\alpha^2=1+\zeta^2\). The orthogonal local change of fields gives interacting momenta

\[
p_0=(\zeta,i\alpha),\quad p_2=(\zeta,-i\alpha),\quad
p_1=(-2\zeta,0),\quad
p_A=(-2\zeta-\zeta^{-1},0),\quad
p_U=(\zeta+\zeta^{-1},-i\alpha).
\tag{X18}
\]

These describe the local sine-Liouville hierarchy and the two primaries. They do not identify a global fermion Hilbert space with a bosonic one.

Fateev's primary source [hep-th/0103014v2](https://arxiv.org/abs/hep-th/0103014), equation (131), states exactly

\[
R_s(a,\eta)=\frac{A_{s(a),\eta}}{A_{a,\eta}}.
\]

The rendered original p. 33 was inspected: equation (132) uses the shifted momentum in the coupling power. The line on that page, **“Until now we considered spinless fields and non-compactified field Φ.”**, also prevents treating its normalization formula as a specification of our global sectors. The long-root expression on p. 6, equations (24)–(25), has the same reflected/original orientation. The [fresh normalization audit](referee_sl21_reflection_bridge_audit_report.md) supplies the full source specialization and rechecks the 2017 form-factor labels independently.

The resulting reference functions in the zero-constant Gaussian convention are

\[
\mathcal G_A(\lambda_1)=
[\lambda_1\gamma_f(-t)]^{1/t-1}
\frac{\Gamma(2-t)\Gamma(2-1/t)}{\Gamma(t)\Gamma(1/t)},
\]
\[
\mathcal G_U(\lambda_0,\lambda_2)=
\frac{\lambda_0}{1-t}
\left(\frac{\sqrt{\lambda_0\lambda_2}}{-t}\right)^{2/t-2}
\frac{\Gamma(2-1/t)}{\Gamma(1/t)},
\quad \gamma_f(x)=\frac{\Gamma(x)}{\Gamma(1-x)}.
\tag{X19}
\]

The subtle point in the second value can be independently checked locally. Set \(Q=(2\zeta)^{-1}\), \(a_*=\zeta+\zeta^{-1}\), and \(C_*=\sqrt{\lambda_0\lambda_2}/(-t)\). The source ratio at \(a=a_*+\eta\), angular momentum \(-\alpha\), has expansion

\[
R=\eta\zeta C_*^{2/t-2}
\frac{\Gamma(2-1/t)}{(t-1)\Gamma(1/t)}+O(\eta^2).
\tag{X20}
\]

Here \(\eta\) is a small radial-momentum displacement, not the source's angular label. The reflected exponent is \(-\zeta-\eta\). Its product with \(Y_0\) has scalar factor \(|x|^{-2+2\zeta\eta}\) and fused field \(e^{-\eta\phi}\). The minus action insertion therefore supplies the identity pole \(-\lambda_0/(\zeta\eta)\). Multiplication by (X20) gives precisely \(\mathcal G_U\), including its sign. Finite regular terms in that reflected family are killed by the simple zero. Our independent script checks the unsimplified reflection ratio near two irrational couplings as well as the zero-times-pole limit.

The directly computed local product
\(U(0)Y_0(x)=|x|^{2t-4}(F+xCF+\bar x\bar CF+x\bar x C\bar C F+\cdots)\), where \(F=e^{(d+T)(X+\bar X)}\), gives

\[
U_{\min}^0=\lim\left[U^0-
\frac{\lambda_0^0}{t-1}\epsilon^{2t-2}F^0
-\frac{\lambda_0^0}{t}\epsilon^{2t}C\bar C F^0\right].
\tag{X21}
\]

Both signs follow from the lower limit of the radial integral with the minus action insertion. The next angularly invariant term vanishes as \(\epsilon^{2t+2}\). Equation (X21) matches the displayed canonical endpoint through mass order two after multiplication by \(t^2\lambda_2\); it does not match the full endpoint at all orders.

## 6. Explicit finite Gaussian normalization and counterterm transport

Let \(H\) be a constant symmetric change of the full finite Gaussian subtraction, with unchanged singular covariance. Define

\[
V_p^H=Z_pV_p^0,\qquad Z_p=e^{-p^THp/2}.
\]

At fixed physical action \(\lambda_i^0=Z_i\lambda_i^H\), while
\(Z_pZ_q/Z_{p+q}=e^{p^THq}\). Thus OPE coefficients, not their singular powers, change. The reference values become

\[
G_A^H=Z_A Z_1^{1/t-1}\mathcal G_A(\lambda_1^H),\qquad
G_U^H=Z_U Z_0^{1/t}Z_2^{1/t-1}
\mathcal G_U(\lambda_0^H,\lambda_2^H).
\tag{X22}
\]

Take the physically relevant change \(H=\operatorname{diag}(c,c,0)\) in the original three-field basis. It changes only the two physical bosons; the fermion OPE is fixed. Direct multiplication of the momenta yields

\[
Z_0=Z_2=1,\quad Z_1=e^{ct},\quad Z_A=e^{c(1-t)^2/t},
\quad Z_U=e^{c(1/t-1)},\quad Z_d=e^{c/t}.
\tag{X23}
\]

Consequently both functions in (X22) have the extra factor \(e^{c(1/t-1)}\) when written in terms of the **new** couplings. The two needed fusion coefficients are
\(E_{U0}=e^{c(1-t)}\) and \(E_{02}=e^{-ct}\). The endpoint and action subtractions must be

\[
U_{\min}^H=\lim\left[U^H-
\lambda_0^H e^{c(1-t)}
\left(\frac{\epsilon^{2t-2}}{t-1}F^H+
\frac{\epsilon^{2t}}{t}C\bar C F^H\right)\right],
\]
\[
\kappa^H=-\frac{\lambda_0^H\lambda_2^H}{t}
e^{-ct}\epsilon^{2t},\qquad \kappa_R^H=0.
\tag{X24}
\]

The equality \(U_{\min}^H=Z_UU_{\min}^0\) holds through the computed order. Holding the old numerical counterterm coefficients fixed would violate that conversion.

The current prescription must be transported too. Split \(H=H_L+H_R\), with no extra finite mixed chiral contraction, and normalize the leading current Wick vertices to coefficient one in each new basis. Both mixed current pairs acquire \(Z_d\). Constant subtraction changes give no derivative contractions. Therefore the transported Ward endpoints are

\[
B_{AA}^{H,\mathrm W}=\gamma_E^2\lambda_1^H
\frac{Z_1Z_d}{Z_A}F_A^H
=e^{2c}\gamma_E^2\lambda_1^HF_A^H,
\quad \gamma_E=\frac{t}{1-t},
\]
\[
B_{BB}^{H,\mathrm W}=t^2\lambda_2^H
\frac{Z_2Z_d}{Z_U}U_{\min}^H
=e^ct^2\lambda_2^HU_{\min}^H.
\tag{X25}
\]

These formulas transport the already selected finite mixed prescription; they do not assert all-order endpoint matching. If the original physical currents themselves are kept fixed, instead of renormalizing their leading Wick vertices to one, remove the common \(Z_d\) in the unspecialized formulas. Their endpoint operators then stay fixed.

For the reference candidate amplitudes \(\mathcal A^0=\gamma_E^2\lambda_1^0\mathcal G_A\), \(\mathcal B^0=t^2\lambda_2^0\mathcal G_U\), equations (X22)–(X25) give, at fixed action,

\[
\mathcal A^{H,\mathrm W}=e^{c/t}\mathcal A^0,
\qquad \mathcal B^{H,\mathrm W}=e^{c/t}\mathcal B^0.
\tag{X26}
\]

On positive real coupling branches and generic \(-1<t<0\), their ratio reduces to

\[
\boxed{\frac{\mathcal B^0}{\mathcal A^0}
=-\left[-\gamma_f(t)\frac{\lambda_0^0\lambda_2^0}{\lambda_1^0}\right]^{1/t}.}
\tag{X27}
\]

To verify this, cancel the common \(\Gamma(2-1/t)/\Gamma(1/t)\), use
\(\Gamma(2-t)=(1-t)\Gamma(1-t)\), and then
\(\gamma_f(t)\gamma_f(-t)=-t^{-2}\). The bracket in (X27) is positive on the stated real branch. For complex branches, retain the uncombined expression
\((1-t)^2\lambda_2\mathcal G_U/(\lambda_1\mathcal G_A)\) with its original logarithms.

The invariant statement is at **fixed action**. In the new coupling coordinates the right side of (X27) carries an additional \(e^{-c}\). Indeed \(\lambda_1^H=e^{-ct}\lambda_1^0\), while \(\lambda_0^H=\lambda_0^0\), \(\lambda_2^H=\lambda_2^0\); this exactly cancels that factor. Unrelated independent constant rescalings of the two charge pairs can still change the ratio. It has mass dimension two, but no identification with a physical particle mass has been established.

## 7. Correct the cycle interpretation; retain the conditional sector statement

For a convergent regulated integral and constant field displacement \(\Delta\), direct substitution gives

\[
G_{\mathcal C+\Delta,p}(\lambda_i)
=e^{p\cdot\Delta}G_{\mathcal C,p}(\lambda_i e^{p_i\cdot\Delta}).
\tag{X28}
\]

This is an exact change of variables, not a deformation through a possibly divergent region. For the local period
\(\Delta\phi=i\pi/\zeta\), \(\Delta\vartheta=\pi/\alpha\), the screening pairings are
\(p_0\Delta=2\pi i\), \(p_2\Delta=0\), \(p_1\Delta=-2\pi i\).
The direct pullback therefore has logarithm changes \((+2\pi i,0,-2\pi i)\). The preceding report's V21 has their negatives. Those are correct for the compensated path

\[
\mathcal C_\tau=\mathcal C_0+\tau\Delta,\qquad
\lambda_i(\tau)=\lambda_i(0)e^{-\tau p_i\Delta},\qquad 0\le\tau\le1.
\tag{X29}
\]

Its pulled-back action is constant. Consequently the transported insertion gets \(e^{\tau p\Delta}\); at the end the numerical couplings return to their initial values, with their logarithms and cycle transported. Equations (X19) have precisely the characters \(G_A\mapsto s^{-1}G_A\), \(G_U\mapsto sG_U\) on this path. A direct comparison at fixed couplings and this compensated path must not be conflated. Our new script also checks (X28) and (X29) by independent convergent finite-dimensional integrals; those finite integrals are a sign check, not a construction of the Toda continuum measure.

**Replacement for `briefs/referee_sl21_vacuum_limits_report.md`, line 335 (no edit applied):**

> Under the period \(\Delta\phi=i\pi/\zeta\), \(\Delta\vartheta=\pi/\alpha\), transport the cycle as \(\mathcal C_\tau=\mathcal C_0+\tau\Delta\) and compensate its couplings by \(\lambda_i(\tau)=\lambda_i(0)e^{-\tau p_i\cdot\Delta}\). The pulled-back action is constant. This path has the logarithm shifts (V21), and its transported vacuum acquires the endpoint characters below. A direct pullback at fixed couplings has the opposite logarithm shifts, as in (X28) of the present audit.

To state exactly what a reference quantum sector would provide, assume a regular meromorphic operator family with the specified reflection normalization, the stated form factors and charge selection, crossing including disconnected terms, and completeness of the reference particle domain. The normalized nonvacuum neutral form factors of these primaries vanish. Under those assumptions the corresponding quadratic forms are \(F_A=G_AI\) and \(U_{\min}=G_UI\) on that sector. This conditional assertion does not follow from a few vanishing numerical form factors; the all-particle recursion is required, and was independently rechecked by the normalization audit.

There is no conflict with nontrivial characters on a cover of sectors. On \(\bigoplus_n\mathcal H_n\), with projectors \(P_n\), the reference operators can be

\[
\mathbb U=G_U\sum_n s^nP_n,\qquad
\mathbb A=G_A\sum_n s^{-n}P_n.
\tag{X30}
\]

This is an explicit conditional sector realization, not a construction of the charge intertwiners. If interface sectors \(\mathcal H_{l,r}\) and scalar endpoint limits exist, \(N=r-l\) gives

\[
\Delta B_{AA}=\mathcal A s^{-r}(1-s^N),\qquad
\Delta B_{BB}=\mathcal B s^r(1-s^{-N}).
\tag{X31}
\]

The original theory still needs both the all-order identification of its Ward endpoints with these reference operators and the interface sectors/common charge domain. Neither a homogeneous particle vacuum nor a one-point value supplies those requirements. Also, the nonanalytic coupling powers in (X19) prevent taking the finite perturbative remainder to infinite volume without a separate uniform argument.

## 8. Evidence, limits and the next calculation

New root evidence is in `referee_scratch/sl21_triple_extension_20261024/`:

- `independent_checks.py` / `.out`: **78 PASS, 0 FAIL**, from independent vector contractions, Laurent coefficients, collision degrees and remainder integrals, cocycles, finite Gaussian factors, the primary reflection ratio, the conditional amplitude ratio and the exact cycle covariance.
- `independent_checks_initial.out`: retained failure of an initial branch-orientation check. The oscillator exchange factor had the two spatial orders reversed. Re-deriving the positive-imaginary-separation boundary values gives (X16); the endpoint weight itself is unchanged. The corrected check includes the actual cocycle and both spatial orders.
- `checkpoint.md`: resumable calculation and remaining work.

The fresh [kernel audit](referee_sl21_triple_kernel_report.md) adds **23** checks, including unexpanded contour quadrature. The fresh [reflection bridge audit](referee_sl21_reflection_bridge_audit_report.md) adds **87** checks and explicit source-page evidence. Their derivations were checked against the root calculation; old scripts were not used as independent evidence. Counts are diagnostics, not evidence of a global symmetry algebra.

The most direct next finite-order calculation is the coefficient of \(\lambda_0\lambda_2\) in the same four-word defect, with its regulated two-source integral, the composite-potential counterterm and the two endpoint subtractions in (X24) included together. The companion \(\bar A\) projection and four-current overlap must ultimately use that same product. On the quantum side, (X24)–(X27) are quantitative tests for an explicit continuum completion; an all-order Ward-endpoint matching and an interface-domain construction remain necessary.

Final preservation check: all **22** recorded input hashes are unchanged. Tracked working-tree and staged diffs are empty. The only additions to non-ignored status are the three new reports from this pass; the new scratch directories use the repository's existing exclusion rules.

## Integration

None performed. **Verified:** the local momenta and first-order three-current kernels; the full additional-current mixed face; the simultaneous collision estimate; propagation of the binary finite contact; the conditional reflection values and finite Gaussian conversion. **Modified:** the interpretation of V21 is compensated cycle/coupling transport, with a concrete replacement above. **Extended:** the actual restricted three-contour identity and a normalization-invariant conditional endpoint comparison at fixed action. **Rejected:** deriving the triple calculation solely from assumed associativity, or using reflection values as established original-theory endpoint limits. **Unresolved:** double-Yukawa triple matching, arbitrary-test/common four-current products, all-order current conservation and endpoint matching, interface sectors and weak limits, and the symmetry algebra.
