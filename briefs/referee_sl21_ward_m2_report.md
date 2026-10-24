# Root-system sl(2|1) Toda: the order-mass-squared Ward identity

Date: 2026-10-23. Model: Codex (GPT-6).

Commit reviewed: `745586e6a6f7bca5ed91a1e9cf89d6ef5df387f0`.

Task: independently verify the bosonization and candidate Bernard–LeClair currents, then combine the composite-potential residue with the regulated double-Yukawa contribution. Read the three requested preceding reports. Their scripts were neither imported nor rerun as evidence. No existing file was changed, no subagent was used, and no git mutation or book render was performed.

New evidence: `referee_scratch/sl21_ward_m2_20261023/check_ward.py` and `check_ward.out`. All 46 checks pass. The derivation below specifies the quantum prescription; the result is conditional on that prescription, not a unique quantization inferred from the classical action.

## Outcome and verdicts

The bosonization and both candidate currents survive independent verification. The apparent order-mass-squared obstruction **cancels in a specified point-splitting prescription**. The original component action supplies exactly the required cutoff-dependent composite coefficient after canonical normalization and bare-to-Wick conversion. An additional finite, independently normal-ordered composite potential would leave an obstruction for this current.

Write $t=b^2$, normalize the Euclidean perturbation as $\pi^{-1}\int d^2z\,[\lambda_0Y_0+\lambda_2Y_2+\lambda_1V_1+\kappa(\epsilon)P]$, and take the full Yukawa OPE coefficient to be $+1$. Then the non-derivative Ward class at order $\lambda_0\lambda_2$ is

\[
\boxed{\left[\bar\partial J_--\partial H_-\right]_{m^2}
=\left(\kappa(\epsilon)+\frac{\lambda_0\lambda_2}{t}\epsilon^{2t}\right)[\mathcal A]
+o(1),}
\tag{1}
\]

Here the brackets denote the local anomaly class extracted from the contour Ward identity (23), modulo local holomorphic total derivatives. They do not assert equality of unrenormalized pointwise current operators. Its representative is

\[
\mathcal A=\frac12:
(T\cdot\partial X)((a_0-a_2)\cdot\partial X)
e^{(d+T)\cdot X+T\cdot\bar X}:.
\tag{2}
\]

The precise domain of the cutoff calculation is $\operatorname{Re}t>-1$, with $t\ne0$; the familiar singular parameter $t=1$ is also excluded when discussing both currents together. The remainder in the collision channel is $O(\epsilon^{2\operatorname{Re}t+2})$. Additional counterterms can be necessary outside this domain. No claim about those regimes is made.

Define

\[
\boxed{\kappa(\epsilon)=\kappa_R-
\frac{\lambda_0\lambda_2}{t}\epsilon^{2t}.}
\tag{3}
\]

Then (1) leaves precisely $\kappa_R[\mathcal A]$. The prescription with **zero finite composite coupling**, $\kappa_R=0$, passes this order of the conservation test. Sections 5–6 give the regulated integral, its sign and normalization, and an explicit derivative improvement. This is stronger than testing the two terms separately, but weaker than an all-order conservation theorem.

| Item | Verdict | Independent evidence |
|---|---|---|
| Component action and root normalization used as input | **Correct** in the specified complex conventions | Source (3.4), root substitution, Wick rotation and field rescaling in section 2 |
| Bosonized momenta, fermion parity and full Yukawa OPE sign | **Correct**, with the explicit cocycle/bilinear convention below | Fermion contractions and vertex expansion in section 3 |
| Both first-order candidate currents | **Correct** as formal ultraviolet currents | Wick expansion in both OPE orientations, section 4 |
| Composite-potential residue in the preceding report | **Correct** | Independent second-order Taylor expansion, section 4 |
| Double-Yukawa contact and cancellation for (3), $\kappa_R=0$ | **Correct through this mass order**, within this regulator and local improvement class | Two-end collision-circle calculation, exact polynomial identities and numerical angular integration |
| Cancellation with an arbitrary nonzero finite composite coupling | **Wrong for the fixed candidate and ordinary local improvements** | Nonzero class $\kappa_R[\mathcal A]$; section 7 gives the replacement wording |
| Unique quantum theory, all-order conserved charges, charge algebra and soliton representation | **Unclear; not established here** | These require calculations beyond this Ward test |

No identification with either affine sl(2|1) or affine sl2 is assumed or established.

## 1. Sources actually checked

Bernard–LeClair, *Quantum group symmetries and non-local currents in 2D QFT*, Commun. Math. Phys. 142 (1991) 99–138, was read from the supplied local PDF:

`/home/gustav/Zotero/storage/QJH8RCYM/Bernard and LeClair - 1991 - Quantum group symmetries and non-local currents in.pdf`.

Printed pages 105 and 123 were rendered and inspected. Page 105 says exactly: **“the residues of these OPE's are total derivatives”**. Equations (2.9)–(2.13) distinguish this first-order test from the additional all-order question. Page 123 explicitly introduces (5.3)–(5.4) as charges **“to first order in perturbation”**. The new calculation uses their residue method, with the area convention made explicit below. [Published source](https://doi.org/10.1007/BF02099173).

Evans–Madsen, [hep-th/9703065v2](https://arxiv.org/abs/hep-th/9703065v2), printed pages 8–9, was checked directly in the source PDF already present in the repository scratch material. Equation (3.4), including its $1/2$ in the odd-root double sum and its $i\mu_j$ Yukawa coefficient, was visually inspected. The exact source statement on page 8 is: **“If desired one can introduce a dimensionless coupling β by taking L → L/β²”**. Here that source coupling is renamed $g$, to avoid silently identifying it with the book's physical coupling.

Litvinov–Spodyneiko, [arXiv:1609.06271](https://arxiv.org/pdf/1609.06271), equations (4.14)–(4.15), was checked independently of the preceding report. Its exact wording is **“then it also commutes with the dressed screening field”**. This is a statement about a screening commutant. It supplies no substitute for the massive Ward calculation below. No novelty claim or correction to a published paper is made.

## 2. Action-to-free-field normalization

Use the source component action with the overall multiplier $g^{-2}$, the bilinear dot product without complex conjugation, and

\[
\alpha_0=(-1/\sqrt2,i/\sqrt2),\quad
\alpha_1=(\sqrt2,0),\quad
\alpha_2=(-1/\sqrt2,-i/\sqrt2).
\]

Direct substitution gives the Gram matrix

\[
(\alpha_i\cdot\alpha_j)=
\begin{pmatrix}0&-1&1\\-1&2&-1\\1&-1&0\end{pmatrix}.
\]

In particular, the odd-root double sum in (3.4) gives **one**, not two, copies of

\[
V_{02}=\mu_0\mu_2 e^{(\alpha_0+\alpha_2)\cdot\phi}.
\tag{4}
\]

The Yukawa matrix is $\sum_{j=0,2}\mu_j\alpha_j\alpha_j^t e^{\alpha_j\cdot\phi}$. At equal masses it reproduces the matrix in the kink report. This verifies the part of that report used here; its zero-mode quantization and representation examples are not re-refereed or used to identify a charge algebra.

Keep the earlier report's $\partial_\pm=\partial_{t_M}\pm\partial_x$. Wick rotate $t_M=-i\tau$, put $z=x+i\tau$, and use ordinary positive area $d^2z=dx\,d\tau$. Thus $\partial=\tfrac12(\partial_x-i\partial_\tau)$, $\bar\partial=\tfrac12(\partial_x+i\partial_\tau)$, and $\bar\partial(1/z)=\pi\delta^{(2)}(z)$. Set

\[
\phi=\frac{g}{\sqrt{4\pi}}\Phi,\qquad
\psi_-=e^{i\pi/4}\frac{g}{\sqrt{2\pi}}\psi,\qquad
\psi_+=e^{-i\pi/4}\frac{g}{\sqrt{2\pi}}\bar\psi.
\tag{5}
\]

The Euclidean free action becomes

\[
S_0=\frac1{2\pi}\int d^2z\,
\bigl(\partial\Phi\cdot\bar\partial\Phi
+\psi\cdot\bar\partial\psi+\bar\psi\cdot\partial\bar\psi\bigr).
\tag{6}
\]

Consequently $\Phi=\varphi+\bar\varphi$, with

\[
\varphi^a(z)\varphi^b(w)\sim-\delta^{ab}\log(z-w),\qquad
\psi^a(z)\psi^b(w)\sim\frac{\delta^{ab}}{z-w}.
\]

The unrescaled Euclidean Yukawa term is $-i\mu_j(\alpha_j\cdot\psi_+)(\alpha_j\cdot\psi_-)/g^2$. Reordering the fermions after (5) makes it $+i\mu_j(\alpha_j\cdot\psi)(\alpha_j\cdot\bar\psi)/(2\pi)$. Define

\[
Y_j=i:(\alpha_j\cdot\psi)(\alpha_j\cdot\bar\psi)
e^{-ib\alpha_j\cdot\Phi}:,\qquad
\boxed{b=\frac{ig}{\sqrt{4\pi}},\quad t=-\frac{g^2}{4\pi},\quad
\lambda_j=\frac{\mu_j}{2}\quad(j=0,2).}
\tag{7}
\]

These are bare/canonical normalization identities in the stated complex continuation. For real source $g$, $b$ is imaginary. Choosing real $b$ continues the source coupling and contour. This does not specify a positive real quantum theory, identify $b$ with the book's physical beta, or determine $q$.

With this ordinary-area convention, $S_{\rm int}=\pi^{-1}\int\lambda_aV_a$ gives $\bar\partial J=\lambda_a\operatorname{Res}_{w=z}V_a(w)J(z)$. The factor $1/\pi$, rather than a copied $1/(2\pi)$, is essential: Bernard–LeClair use the two-form area convention $d^2z=i\,dz\,d\bar z$.

The even-root coupling has its own usual vertex normalization. It is counted as order $m^2$, and its finite renormalized value can be held fixed separately. The convenient classical choice $\mu_1=\mu_0\mu_2$ is not, by itself, a quantum relation between independently Wick-normalized couplings of different dimensions.

## 3. Independent bosonization, including the sign needed at second order

Let $\chi_j=\alpha_j\cdot\psi$, $j=0,2$, and similarly for the antiholomorphic fermions. Their contractions give

\[
\chi_0(z)\chi_2(w)\sim\frac1{z-w},\qquad
\chi_0\chi_0\sim\chi_2\chi_2\sim0,\qquad
:\chi_0\chi_2:=i:\psi^u\psi^v:.
\]

Taking $h(z)h(w)\sim-\log(z-w)$, the expansion

\[
:e^{ih(w+r)}::e^{-ih(w)}:
=r^{-1}\bigl(1+ir\partial h(w)+O(r^2)\bigr)
\]

verifies both the fermion contraction and $ :\chi_0\chi_2:=i\partial h$. Left and right fermions anticommute; the antiholomorphic sector therefore includes a Klein factor. We fix the full local-vertex cocycles by bosonizing the **$i\chi_j\bar\chi_j$** of (7). In the mixed product the one fermionic Wick crossing gives $-1$, and the two factors of $i$ give another $-1$. The resulting leading OPE coefficient is **+1**. This convention fixes the relative sign between the Yukawa product and (4); it cannot be discarded during a second-order calculation.

For $X=(\varphi^u,\varphi^v,h)$, the chiral momenta are

\[
a_0=(ib/\sqrt2,b/\sqrt2,i),\qquad
a_2=(ib/\sqrt2,-b/\sqrt2,-i),\qquad
T=a_0+a_2=(i\sqrt2 b,0,0),\qquad a_1=-T.
\]

Their independently computed products are

\[
a_0^2=a_2^2=-1,\quad a_0\cdot a_2=1-t,\quad T^2=-2t.
\tag{8}
\]

Write $Y_j=:e^{a_j\cdot X+a_j\cdot\bar X}:$ with the full-vertex cocycles just specified, $P=:e^{T\cdot X+T\cdot\bar X}:$, and $V_1=:e^{-T\cdot X-T\cdot\bar X}:$. Then

\[
\boxed{Y_0(w+r)Y_2(w)
=|r|^{2t-2}\bigl(P(w)+\text{descendants}\bigr).}
\tag{9}
\]

The Yukawa fields have weights $(1/2,1/2)$, whereas $P$ has weights $(t,t)$. Thus $[\lambda_j]={\rm mass}$, $[\kappa]={\rm mass}^{2-2t}$, and the factor $\epsilon^{2t}$ in (1) is dimensionally necessary. Replacing it by a dimensionless constant would already be inconsistent.

The direction $(0,\sqrt2,ib)$ is orthogonal to all the perturbing momenta, with norm $2-t$. This local free-field observation says nothing about its global compactification or allowed soliton sectors.

## 4. Recompute the candidate currents and both residue orientations

Put

\[
d=\frac{2T}{T^2}=(-i\sqrt2/b,0,0),\qquad C=a_0,\qquad
D=a_0-a_2,\qquad k=d+T=(1-t)d.
\]

The candidates are

\[
J_+=:e^{-d\cdot X}:,\qquad
J_-=: (C\cdot\partial X)e^{d\cdot X}:.
\tag{10}
\]

The $h$-derivative in $J_-$ is a fermion bilinear. Both currents are even. The derivative prefactor has been checked in the fermionic and bosonized descriptions, not assumed from a root-coroot rule.

For an insertion $V_a(w)$, Wick's theorem at fixed current position $z$, $r=w-z$, gives

\[
V_a(w)J_-(z)=r^{-a\cdot d}:
\left(C\cdot\partial X(z)+\frac{C\cdot a}{r}\right)
e^{d\cdot X(z)+a\cdot X(z+r)}:.
\tag{11}
\]

Taylor expanding the exponential directly yields

\[
\operatorname{Res}_w V_{a_0}(w)J_-(z)=0,\qquad
\operatorname{Res}_w V_{a_2}(w)J_-(z)
=-t\partial:e^{(d+a_2)\cdot X}:.
\tag{12}
\]

The even-root OPE is regular. For $J_+$, the even-root residue is

\[
-\frac{t}{1-t}\partial:e^{(-d+a_1)\cdot X}:,
\]

and the two Yukawa and composite OPEs are regular. This independently recovers the first-order results of both earlier reports.

The double-insertion Ward calculation also needs the opposite OPE orientation: move the **charge-current position** around a fixed perturbation. For $d\cdot a=1$ this gives the vector $C-(C\cdot a)d$, not the vector $C+(C\cdot a)a$ in (11). In particular, for

\[
K_j=:e^{(d+a_j)\cdot X+a_j\cdot\bar X}:,
\]

the contour variation $\delta=\operatorname{Res}_z J_-(z)$ satisfies

\[
\boxed{\delta Y_0=\partial K_0,\qquad
\delta Y_2=-\partial K_2.}
\tag{13}
\]

There is no contradiction with (12): moving the expansion point changes total derivatives. **Using (12) as though it were the charge variation of each insertion would miss one collision endpoint.**

For the composite, abbreviate the common vertex by

\[
E=:e^{k\cdot X+T\cdot\bar X}:.
\]

The fixed-current residue is

\[
R_P=:\left[(C\cdot\partial X)(T\cdot\partial X)
-\frac t2\{T\cdot\partial^2X+(T\cdot\partial X)^2\}\right]E:,
\tag{14}
\]

and

\[
R_P=\partial:\left[-\frac t2(T\cdot\partial X)E\right]:+\mathcal A.
\tag{15}
\]

This verifies the preceding report's obstruction. At level one, an arbitrary improvement is $(v\cdot\partial X)E$. Its second-derivative coefficient forces $v=-tT/2$, leaving (2). A local $\bar\partial$-improvement cannot supply a term at antiholomorphic descendant level zero. Thus this remainder is nontrivial in the ordinary local differential-polynomial improvement class.

The charge-contour residue instead is

\[
R_Q:=\delta P
=:\left[C\cdot\partial^2X+(C\cdot\partial X)(d\cdot\partial X)
+\frac t2\{d\cdot\partial^2X+(d\cdot\partial X)^2\}\right]E:.
\tag{16}
\]

It obeys $R_Q=\partial:[(D\cdot\partial X)E/2]:-\mathcal A$. Equations (14) and (16) have opposite non-derivative classes. Their sum is a total derivative; the script checks that identity as well.

## 5. The regulated double-Yukawa integral

Use circular exclusion $|w_0-w_2|>\epsilon$, with all spectator fields outside the collision neighborhood. Long-distance boundaries are held fixed; only ultraviolet contact terms are computed here. The mixed term in $e^{-S_{\rm int}}$ has coefficient

\[
\frac{\lambda_0\lambda_2}{\pi^2}
\int_{|w_0-w_2|>\epsilon}d^2w_0d^2w_2\,Y_0(w_0)Y_2(w_2).
\tag{17}
\]

The two cross terms cancel the $1/2!$. From (13), its charge variation is the integral of

\[
\partial_{w_0}\bigl(K_0(w_0)Y_2(w_2)\bigr)
-\partial_{w_2}\bigl(Y_0(w_0)K_2(w_2)\bigr).
\tag{18}
\]

Let $w=w_2$, $r=w_0-w_2$, and retain the common fused vertex. The two products in (18) are

\[
F_0=r^{t-2}\bar r^{t-1}E_0(r),\qquad
F_2=-r^{t-2}\bar r^{t-1}E_2(r),
\tag{19}
\]

where

\[
\begin{aligned}
E_0(r)&=:
e^{(d+a_0)\cdot X(w+r)+a_2\cdot X(w)
+a_0\cdot\bar X(\bar w+\bar r)+a_2\cdot\bar X(\bar w)}:,\\
E_2(r)&=:
e^{a_0\cdot X(w+r)+(d+a_2)\cdot X(w)
+a_0\cdot\bar X(\bar w+\bar r)+a_2\cdot\bar X(\bar w)}:.
\end{aligned}
\]

The minus sign in $F_2$ is necessary. It is the factor $1/(w_2-w_0)=-1/r$ from the current's contraction with the other Yukawa. One can obtain (19) directly from the meromorphic current factor $1/[(z-w_0)(z-w_2)]$, without an assumption about the eventual charge algebra. For complex $t$, the paired fractional powers in (19) mean $r^{-1}|r|^{2t-2}$, fixing their angular continuation.

As a separate endpoint check, the three-operator Wick product has the current-dependent factor

\[
\frac{1}{(z-w_0)(z-w_2)}
\left[C\cdot\partial X(z)+\frac{1}{z-w_0}
-\frac{1-t}{z-w_2}\right]e^{d\cdot X(z)}.
\]

Its two residues, after suppressing the common Yukawa contraction and the corresponding normal-ordered exponential, are

\[
\frac{(d+a_0)\cdot\partial X(w_0)}r-\frac{2-t}{r^2},\qquad
\frac{(d+a_2)\cdot\partial X(w_2)}r+\frac{2-t}{r^2}.
\]

They are exactly $\partial_{w_0}F_0$ and $-\partial_{w_2}F_2$ with their suppressed factors restored. Both are checked independently by the script.

With the ordinary area convention,

\[
\int_{|r|>\epsilon}\partial_r F\,d^2r
=-\frac\epsilon2\int_0^{2\pi}e^{-i\theta}
F(\epsilon e^{i\theta})\,d\theta
+\text{outer boundary}.
\]

The two endpoints in (18) therefore give the **sum** $F_0+F_2$ on the inner boundary. The resulting local term is

\[
-\frac\epsilon2\int_0^{2\pi}d\theta\,e^{-i\theta}
r^{t-2}\bar r^{t-1}(E_0-E_2)
=-\pi\epsilon^{2t}B+O(\epsilon^{2\operatorname{Re}t+2}),
\tag{20}
\]

where

\[
\boxed{B=:\left[\frac12d\cdot\partial^2X
+\frac12(d\cdot\partial X)^2
+(C\cdot\partial X)(d\cdot\partial X)\right]E:.}
\tag{21}
\]

To see the coefficient without a formal residue shortcut, the $r^2\bar r^0$ term of $E_0-E_2$ is exactly the bracket in (21). The boundary measure and contraction supply angular weight $e^{-2i\theta}$. This selects that Taylor coefficient. The next selected term has $r^3\bar r^1$, two extra powers of the cutoff. The linear term does not survive angular integration.

Crucially,

\[
\boxed{B=\partial:\left[\frac12(d\cdot\partial X)E\right]:
-\frac1t\mathcal A.}
\tag{22}
\]

Thus the double-Yukawa contact has exactly the same nontrivial local class as the composite residue, with the coefficient required in (1). The script derives (21) from the two Taylor expansions, rather than inserting (22) as an assumed cancellation.

The $Y_0Y_0$ and $Y_2Y_2$ channels are regular: their full leading factor is $|r|^2$. They add no ultraviolet contact of this order. $J_+$ has zero charge-contour residues with both Yukawas and the composite, so this mixed collision does not obstruct it. Its even-root correction is already covered by the first-order calculation.

## 6. Combine the terms and fix the quantum prescription

The composite insertion contributes $-\kappa\pi^{-1}\int R_Q$ to the variation of the perturbative weight. Combining (17) and (20) gives

\[
\delta e^{-S_{\rm int}}\big|_{02,P}
=-\frac1\pi\int d^2w\,
\underbrace{\bigl(\kappa R_Q+\lambda_0\lambda_2\epsilon^{2t}B\bigr)}_{\mathcal W_Q}
+\text{spectator/outer terms}+o(1).
\tag{23}
\]

This fixes the relative sign and all factors of $\pi$. Taking the class of (23), or using the fixed-current orientation of (14), gives (1).

There is also an explicit full-polynomial check, beyond equality of the remainders:

\[
\boxed{R_Q-tB=\partial:[(C\cdot\partial X)E]:.}
\tag{24}
\]

For $\kappa=-\lambda_0\lambda_2\epsilon^{2t}/t$, (23) is precisely

\[
\mathcal W_Q=\partial:[\kappa(C\cdot\partial X)E]:.
\tag{25}
\]

It is an allowed local derivative improvement of the contour Ward identity. In a fixed-current descent representative the double contact may be represented by $-\lambda_0\lambda_2\epsilon^{2t}B$; then

\[
R_P+tB=\partial:\left[-\frac{t+1}{2}(T\cdot\partial X)E\right]:,
\tag{26}
\]

so an order-$m^2$ improvement representative is

\[
H_-^{(2)}=-\frac{\kappa(t+1)}2:
(T\cdot\partial X)E:.
\tag{27}
\]

The usual ordered-residue first-order representative is $H_-^{(1)}=-t\lambda_2K_2$. Total-derivative contact choices redistribute terms between these local representatives. Equations (23)–(25) are the directly computed contour-descent identity; (27) specifies a compatible polynomial primitive for the contact, not a uniquely determined finite pointwise current in all spectator correlation functions. A full current-operator subtraction prescription must preserve this descent identity. No inverse derivative or nonlocal contact counterterm is used.

Here is a complete prescription for this order:

1. Use the canonical Gaussian propagators in (6), the bilinear/cocycle normalization in (7)–(9), and Wick vertices normalized at one common reference length.
2. Exclude a disk of radius $\epsilon$ around a colliding pair of interaction insertions. Match the bosonic self-contraction subtraction to $G_\epsilon(0)=-\log\epsilon^2$, with zero added constant in these reference-length units.
3. Extend the mixed product to coincident points by analytic continuation from $\operatorname{Re}t>0$. In the strip $-1<\operatorname{Re}t<0$, this means adding the local scalar counterterm in (3). Select $\kappa_R=0$, i.e. add no finite independent $P$ perturbation. Normalize the even-root coupling separately.
4. Use the local derivative improvement (25), or an equivalent descent representative such as (27), for the current contact. Keep the standard external-field Ward contacts distinct from an anomaly in its divergence.

The scalar product itself checks this prescription. For a small fixed outer radius $R$, (9) gives

\[
\int_{\epsilon<|r|<R}d^2r\,|r|^{2t-2}
=\frac\pi t(R^{2t}-\epsilon^{2t}).
\tag{28}
\]

Its divergent contribution to the mixed term in the exponential is

\[
-\frac{\lambda_0\lambda_2}{\pi t}\epsilon^{2t}\int P.
\]

It is canceled by the linear composite insertion with (3) and $\kappa_R=0$. Thus the subtraction needed to define the scalar double product is also the subtraction that removes the Ward obstruction. It is not a separately fitted residue cancellation.

Finally, compare directly with the classical component coefficient. The individual bosonic Yukawa momenta $-ib\alpha_0,-ib\alpha_2$ are null, so their bare exponentials have no bosonic self-normalization factor. The composite momentum instead has square $-2t$. With the matched subtraction above,

\[
e^{(\alpha_0+\alpha_2)\cdot\phi_\epsilon}
=\epsilon^{2t}P.
\tag{29}
\]

Using (4) and (7), its coefficient in the ordinary-area action is

\[
\boxed{\kappa_{\rm source}(\epsilon)
=\frac{\pi\mu_0\mu_2}{g^2}\epsilon^{2t}
=-\frac{\lambda_0\lambda_2}{t}\epsilon^{2t}.}
\tag{30}
\]

This is exactly the required coefficient. The matching includes the Wick-rotation phases, the factor $1/2$ in the source double sum, the fermion contraction sign, and the cutoff power.

For $\operatorname{Re}t>0$, the contact coefficient vanishes with the cutoff and the unrenormalized mixed integral is locally integrable. For $-1<\operatorname{Re}t<0$, which includes weak real source coupling $g$, both contributions diverge but their sum has a finite continuation. At $t=0$ the radial integral is logarithmic and the $1/t$ expressions require a separate limiting prescription; it is excluded here. Taking the classical limit and removing the cutoff are not interchangeable operations.

A different self-contraction constant, a different exclusion shape, or a finite contact extension changes the relation between bare and renormalized coefficients. It must be translated consistently in both the action and the current. The classical component formula alone does not choose all finite quantum counterterms. What has been established is a concrete compatible prescription, including its finite choice, and a calculation of the obstruction if an independent $\kappa_R$ is retained.

## 7. Corrections and replacement wording

No algebraic error was found in the quoted bosonization or composite residue. The preceding reports correctly marked the combined calculation as unresolved. The following are proposed replacements for superseded status statements, not edits made to those files.

For `briefs/referee_sl21_dual_screening_report.md`, lines 148–150, replace the unresolved-calculation paragraph with:

> With canonical Euclidean normalization and a specified circular point-splitting prescription, the mixed double-Yukawa contact contributes the same non-derivative Ward class as the composite potential, with relative coefficient $\lambda_0\lambda_2\epsilon^{2b^2}/b^2$. The source component coefficient, after bare-to-Wick conversion, is $-\lambda_0\lambda_2\epsilon^{2b^2}/b^2$, so the order-mass-squared obstruction cancels modulo a local derivative improvement when no finite independent composite coupling is added. An independent finite composite coupling leaves an obstruction. This is an order-specific, prescription-dependent result; all-order conservation and the charge algebra remain unproved.

For `briefs/referee_sl21_bernard_leclair_report.md`, line 44, and `briefs/referee_sl21_dual_screening_report.md`, line 35, replace the entirely unresolved action-to-$b$ statement with:

> In the explicit canonical Wick rotation of the source action with multiplier $g^{-2}$, $b=ig/\sqrt{4\pi}$ and the two Yukawa coefficients in $S_{\rm int}=\pi^{-1}\int\lambda_jY_j$ are $\lambda_j=\mu_j/2$. This is a formal bare normalization and analytic continuation, not an identification with the book's physical coupling or with a physical quantum-group parameter.

It would be wrong to replace these statements simply by “the nonlocal charges are conserved” or “the symmetry is affine sl2.” The calculation neither proves higher-order conservation nor computes coproducts, mixed chiral/antichiral commutators or Serre relations. The opposite momenta in (10) and their Gram matrix do not settle any of those questions.

## 8. Independent check record and remaining scope

The new script reconstructs the component potential and Yukawa matrix from the roots; checks the fermion signs; generates the fixed-current and fixed-perturbation residues by Wick/Taylor expansion; checks the full derivative identities (24) and (26); derives the canonical coefficient and normal-ordering power; and numerically integrates the actual collision circle with holomorphic and antiholomorphic background jets.

Run from the repository root:

```text
python3 -u referee_scratch/sl21_ward_m2_20261023/check_ward.py
```

The numerical tests use generic irrational $t=\pm\sqrt2/3$, with 8192 angular samples. They compare the boundary integral with $-\pi\epsilon^{2t}B$. The relative error decreases by a factor of approximately four each time the cutoff is halved, as predicted by the next allowed Taylor term:

| $t$ | Error at $\epsilon=0.1$ | Error at $\epsilon=0.0125$ |
|---|---:|---:|
| $+\sqrt2/3$ | $1.12052\times10^{-4}$ | $1.74998\times10^{-6}$ |
| $-\sqrt2/3$ | $7.80236\times10^{-4}$ | $1.21901\times10^{-5}$ |

All 46 final checks pass. During development, the source-matrix comparison needed a common exponential rather than mixed trigonometric representation, and the angular test was extended to a smaller cutoff to reach its stated tolerance. Neither changed the derived Ward coefficient. The symbolic identities, rather than the numerical test alone, supply the generic-$b$ evidence.

This settles the stated ultraviolet collision test **through order $m^2$** in the specified prescription. It does not establish higher mass orders, resonant couplings, other ultraviolet domains, a complete globally defined charge on a scattering Hilbert space, boundary-sector preservation, or action on dressed kink states. The OPEs used here are ordered radial products; a simultaneous global choice of disorder lines and charge Klein factors across those sectors has not been constructed. Before attempting the charge algebra, one must fix finite current-operator subtractions as well as global contour braiding, and retain the derivative improvement in the endpoint analysis. Neither task is settled by choosing between the two affine algebras named in the question.

## Integration

No book integration was performed. **Verified independently:** the action input, canonical field normalization, bosonization including the mixed full-vertex sign, both candidate-current residue tests, and the isolated composite obstruction. **Extended:** the double-Yukawa contact was calculated with an explicit regulator; its coefficient and derivative remainder were combined with the potential, giving (1)–(3) and the source match (30). **Modified interpretation:** the component composite potential can be the contact subtraction required by the Yukawa product, rather than a finite independent ultraviolet perturbation. **Rejected as unsupported:** automatic cancellation for arbitrary finite composite normalization, all-order conservation, and identification of a quantum-affine algebra or Toda S-matrix. The finite prescription is part of the result and must accompany any use of the cancellation. No [Theorem] upgrade is recommended on the basis of this referee calculation alone.
