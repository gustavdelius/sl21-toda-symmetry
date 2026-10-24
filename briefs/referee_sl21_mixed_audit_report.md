# Root-system sl(2|1) Toda: independent audit of mixed-product contacts and continuation

Date: 2026-10-23. Model: Codex (GPT-6).

Commit reviewed: `745586e6a6f7bca5ed91a1e9cf89d6ef5df387f0`.

Task: audit `briefs/referee_sl21_mixed_regulator_report.md` and its checkpoint before extending the calculation. Read `AGENTS.md`. No existing file was edited; no subagents, git mutations, commits, pushes, or book rendering were used.

New evidence: [audit_checks.py](../referee_scratch/sl21_mixed_audit_20261023/audit_checks.py), [its output](../referee_scratch/sl21_mixed_audit_20261023/audit_checks.out), and the derivations below. No earlier check script was imported or rerun. The supplied Bernard–LeClair PDF was freshly extracted and its printed pages 105–106 were rendered and inspected.

## Outcome

The finite first-order current/action contact survives independent verification. In particular, the derivative of the strip boundary is real distributional data; dropping it changes the answer. With the specified circular inner exclusion and current separation $\delta=c\epsilon$, $c>1$, the radial diagonal endpoint coefficients are $t$ and $\gamma^2/t$. The additional finite product normalization in the preceding report changes them to $t^2$ and $\gamma^2$.

There is a gap in the preceding convergence proof. Its count of six integration dimensions for a simultaneous two-current/two-Yukawa collision is appropriate when the entire cluster approaches a **fixed endpoint**. At an interior point, the common coordinate along the charge contour does not scale to zero. There are only five relative dimensions. The raw local degree is therefore $2t-3$, not $2t-2$. A reflection cancellation removes the leading scalar term in a symmetric interior patch, but this is an additional argument, not a consequence of the reported homogeneity check.

The formal mixed identity can nevertheless be supported without relying on the sharper starting domain claimed there. Sections 5–6 give a conservative collision-sector argument starting in $\operatorname{Re}t>2$, using the corrected relative degrees. Every remaining interaction collision then has positive degree after the explicit proper current contacts are treated. Continuing this identity defines the same candidate meromorphic relation at generic coupling, with the prescribed finite contact retained. This audit does **not** claim a complete proof of the bare cutoff limit throughout $1<\operatorname{Re}t\le2$.

This remains a **formal identity of contour insertions in a stated product prescription**, with finite endpoints and separated spectators. The canonical matching is an imposed normalization condition. It is not a construction of interacting Hilbert-space charges, a calculation of their boundary values, or an identification of an affine symmetry algebra.

| Item audited | Verdict | Independent evidence and qualification |
|---|---|---|
| Local bosonization and the two candidate current residues | **Correct locally** | Fermionic Wick contraction, root Gram matrix, and both current-pole vectors |
| Chirality correction and mixed phase in the working coordinate | **Correct** | Free waves, canonical brackets, and the source coordinate footnote |
| Bare-current/action finite contact, equations (17)–(18) of the preceding report | **Correct in its stated nested prescription** | Smooth two-dimensional Cauchy extension, exact contour moments, and fresh quadrature |
| Tail-current contribution and the radial/canonical coefficient distinction | **Correct** | Direct pole expansion; the finite addition is explicitly a normalization choice |
| Action/current/endpoint subtraction chain | **Correct at the stated mass orders** | Endpoint Wick expansion and the three-pole composite Ward response |
| Reported collision powers as raw bounds at arbitrary interior points | **Wrong** | One nonscaling common contour coordinate was included in the small-scale measure |
| Claim that the displayed powers alone prove the limit for $\operatorname{Re}t>1$ | **Unclear as written** | A leading cancellation is needed; the full sharper-domain estimate was not supplied |
| Formal mixed identity with analytic continuation from $\operatorname{Re}t>2$ | **Correct within the specified distributional prescription** | Proper-face decomposition, corrected positive relative degrees, and the continuation argument below |
| Physical boundary closure, same-chirality/Serre relations, or either proposed affine algebra | **Unclear** | None is calculated by this audit |

There is no proposed correction to Bernard–LeClair and no theorem-label recommendation.

## 1. Reconstruct the local input

Use ordinary Euclidean area, $z=x+i\tau$, and

\[
\varphi^a(z)\varphi^b(w)\sim-\delta^{ab}\log(z-w),\qquad
\bar\partial\frac1z=\pi\delta^{(2)}(z).
\]

The roots and momenta are

\[
\alpha_0=\frac{(-1,i)}{\sqrt2},\quad
\alpha_2=\frac{(-1,-i)}{\sqrt2},\quad
p_i=-ib\alpha_i,\quad t=b^2,\quad
T=p_0+p_2=-td,\quad k=d+T.
\tag{A1}
\]

All dot products are bilinear. Direct evaluation gives $p_i^2=0$, $p_0p_2=-t$, $dp_i=1$, $d^2=-2/t$. The local fermions have $\chi_0(z)\chi_2(w)\sim(z-w)^{-1}$. Left and right fermions anticommute.

For four chiral fermions with labels $0,2,0,2$, Wick's rule gives

\[
\frac1{z_{12}z_{34}}+\frac1{z_{14}z_{23}}
=\frac{z_{13}z_{24}}{z_{12}z_{14}z_{23}z_{34}}.
\tag{A2}
\]

This is the exponential-bosonization expression with charges $+i,-i,+i,-i$ and $h(z)h(w)\sim-\log(z-w)$. Point splitting identifies $j=:\chi_0\chi_2:=i\partial h$. The full mixed Yukawa contraction is positive: the fermionic crossing and the two factors of $i$ supply two minus signs. Thus, locally,

\[
a_0=(ib/\sqrt2,b/\sqrt2,i),\quad
a_2=(ib/\sqrt2,-b/\sqrt2,-i),\quad
a_i^2=-1,\quad a_0a_2=1-t,
\tag{A3}
\]

with $X=(\varphi^u,\varphi^v,h)$ and $d$ extended by zero in its third component. This calculation verifies local Wick products; it says nothing about a global bosonization lattice.

The declared action prescription is

\[
S_{\rm int}=\frac1\pi\int d^2w\,
(\lambda_0Y_0+\lambda_2Y_2+\lambda_1V_1+\kappa P),\qquad
\kappa=-\lambda_0\lambda_2\epsilon^{2t}/t,
\tag{A4}
\]

where $Y_i=i:\chi_i\bar\chi_i e^{p_i\Phi}:$, $P=:e^{T\Phi}:$, and $V_1=:e^{-T\Phi}:$. There is no additional finite independent $P$ coupling. The mass orders of $\lambda_0,\lambda_2,\lambda_1$ are $1,1,2$. Wick ordering is at a fixed reference length. This audit checks the consequences of this quantum prescription; it does not establish uniqueness of its identification with a nonperturbative component-action quantization.

For $J_B=:(a_0\partial X)e^{dX}:$, the current-first Wick expansion has double coefficient $-a_0a_i$ and simple-pole vector $a_0-(a_0a_i)d$. They are

\[
\begin{array}{c|cc}
 &\text{double}&\text{simple vector}\\\hline
Y_0&1&d+a_0\\
Y_2&t-1&-(d+a_2).
\end{array}
\tag{A5}
\]

The distributional response subtracts the derivative of the double-pole coefficient. Including the local current-first cocycle $-1$ gives $H_B^{(1)}=-t\lambda_2K_2$ and no $Y_0$ tail. The other current is $J_A=:e^{-d\varphi}:$. It is regular with both Yukawas; its even-root tail is $H_A=-\gamma\lambda_1e^{-k\varphi-T\bar\varphi}$, where $\gamma=t/(1-t)$.

### Single-current mass-squared check

This input was also re-derived rather than inferred from a passing script. The two fermionic point splittings are

\[
\chi_0(0)\chi_2(r)=-r^{-1}+j+rV+\cdots,\qquad
\chi_2(0)\chi_0(r)=-r^{-1}-j-rU+\cdots,
\]

where $U=:\partial\chi_0\chi_2:$, $V=:\chi_0\partial\chi_2:$, and $U+V=\partial j$. Put $x_i=p_i\partial\varphi$, $y_i=p_i\partial^2\varphi$. The normalized quadratic Taylor coefficients are

\[
Q_0=\tfrac12(y_2+x_2^2)-jx_2-V,\qquad
Q_2=\tfrac12(y_0+x_0^2)+jx_0+U.
\]

Differentiating the Yukawa exclusion circle selects $r^2\bar r^0$ from their difference. Its boundary measure and angular integral give $\pi\epsilon^{2t}(Q_2-Q_0)E$, with $E=e^{k\varphi+T\bar\varphi}$. The sign follows from the two minus signs in the fused current/Yukawa product and $\partial_z\Theta(|w-z|-\epsilon)$.

For an independent check of the potential term, write $C=a_0\partial X$, $D=d\partial\varphi$, and $X_T=T\partial\varphi$. The three poles of $J_BP$ have coefficients

\[
tE,\quad (C+tD)E,\quad
\bigl[\partial C+CD+\tfrac t2(\partial D+D^2)\bigr]E.
\]

Applying the derivatives of the delta function from all three poles, with the action sign $-\kappa/\pi$, gives

\[
R_P/E=CX_T-\tfrac t2(\partial X_T+X_T^2).
\]

Direct subtraction now gives

\[
(Q_2-Q_0)E-R_P/t=\partial(CE).
\tag{A6}
\]

Together with (A4), this yields the current subtraction $H_B^{\rm ren}=-t\lambda_2K_2+\lambda_0\lambda_2\epsilon^{2t}L$, $L=CE$, as an insertion in perturbative correlators. This is the local Ward calculation needed below, not a proof of a two-charge algebra.

## 2. Chirality and the status of canonical matching

With $t_M=-i\tau$, $z=x-t_M$. For $\Phi=f(x-t_M)+g(x+t_M)$, the anchored velocity integral is $-f+g$. Therefore the field holomorphic in this $z$ is

\[
\varphi=\tfrac12\left(\Phi-\int_{-\infty}^x\dot\Phi\right),\qquad
\bar\varphi=\tfrac12\left(\Phi+\int_{-\infty}^x\dot\Phi\right).
\tag{A7}
\]

The source's printed-page-105 footnote reads $z=ix^+=(t^E+ix)/2$. Thus our $z=2i\bar z_{\rm BL}$. The preceding report's correction of chirality is justified.

Integrating $[\Phi(x),\dot\Phi(y)]=4\pi i\delta(x-y)$ gives $[\varphi(x),\bar\varphi(y)]=\pi i$ and the two same-chirality brackets $\pm\pi i\operatorname{sgn}(x-y)$. Hence the diagonal separated exchange is $e^{i\pi d^2}=\nu^2$, where $\nu=e^{-i\pi/t}$. The off-diagonal exchange is $\nu^{-2}$. The symmetric cross-chiral product phase is $\nu$ in the cocycle convention used in the audited report.

Spatial derivatives of the constant cross bracket vanish. This supports using a contact-free canonical bare-current exchange as a **normalization condition**. It does not prove that the Euclidean interacting nonlocal currents equal those canonical operators. Moreover, the present contour calculation tests the product with the outer spatial test function constant near the inner segment; it does not determine every possible derivative-of-test-function contact in a full density-product distribution.

The source itself says, on printed page 106, **“We suppose that the same braiding relations hold”**, and describes its calculation as **“to lowest non-trivial order in perturbation theory.”** These are qualifications on the use of equations (2.15)–(2.18), not a derivation of this model's second-order products. No source correction is proposed.

## 3. Independent finite-contact calculation

Keep the audited prescription: the inner current is at zero; its interaction exclusion is $|w|>\epsilon$; the outer current is at imaginary times $\pm\delta$; $\delta=c\epsilon$ with fixed $c>1$. Current poles are two-dimensional Cauchy distributions before the other cutoffs are removed.

A useful independent regularization of those distributions is

\[
C_\eta(z)=\frac{\bar z}{|z|^2+\eta^2},\qquad
P_\eta(z)=-\partial_z C_\eta(z)
=\frac{\bar z^2}{(|z|^2+\eta^2)^2}.
\tag{A8}
\]

Their limits are $1/z$ and its differentiated double pole. The outer-contour integrals, obtained directly by elementary real integrals, are

\[
\begin{aligned}
\int_{\mathbb R}dx\,C_\eta(x+iy)
 &=-i\pi\frac{y}{\sqrt{y^2+\eta^2}},\\
\int_{\mathbb R}dx\,P_\eta(x+iy)
 &=\frac{\pi\eta^2}{2(y^2+\eta^2)^{3/2}},\\
\int_{\mathbb R}dx\,(x+iy)P_\eta(x+iy)
 &=-i\pi\frac{y}{\sqrt{y^2+\eta^2}}
 +\frac{i\pi y\eta^2}{2(y^2+\eta^2)^{3/2}}.
\end{aligned}
\tag{A9}
\]

The first and third integrals can be taken with a symmetric distant cutoff; the difference of upper and lower contours has the same local limit. Send $\eta/\epsilon\to0$ first. The first line then gives $-2\pi i\Theta_\delta(w)$ for the contour difference, and the second gives $-2\pi i\partial_w\Theta_\delta(w)$. The latter is the disputed strip-edge term. It is present even when the residue of a double pole strictly inside the strip is zero.

Set $q=t-1$, $v=d+a_2$, $K_L=v\partial X$, $K_R=v\bar\partial\bar X$. Let $A_L=a_2\partial X$, $A_R=a_2\bar\partial\bar X$, and $D=d\partial X$. Before multiplying the action coefficient, the first-jet kernel is obtained from

\[
e^{Du+A_Lw+A_R\bar w}
\left[\frac{C}{u-w}+\frac q{(u-w)^2}\right]
\left[\frac q{\bar w^2}-\frac{\bar C}{\bar w}\right].
\tag{A10}
\]

The vector identities are $C+qD=-K_L$, $D+A_L=K_L$, and $qA_R-\bar C=tK_R$. Using (A9), its local action contribution divided by $\nu\lambda_2$ is

\[
-2iq\{K_L(I-qJ_w)-tK_RJ_1\}\mathscr U,
\tag{A11}
\]

where $\mathscr U=e^{v(X+\bar X)}$ denotes the fermionic Wick field $i:\chi_2\bar\chi_2e^{(d+p_2)\Phi}:$. The three integrals $I,J_1,J_w$ are those in the audited equation (15).

They can be evaluated without its calculation by integrating horizontal slices outside a unit circle, after setting $\epsilon=1$. If $y=\operatorname{Im}w$, the necessary slice integrals are

\[
\begin{aligned}
\int_{x^2+y^2>1}\frac{dx}{(x-iy)^2}
 &=2\sqrt{1-y^2}\,1_{|y|<1},\\
\int_{x^2+y^2>1}\frac{dx}{x-iy}
 &=\begin{cases}2i\arcsin y,&|y|<1,\\i\pi\operatorname{sgn}y,&|y|>1,\end{cases}\\
\int_{x^2+y^2>1}dx\,\frac{x+iy}{(x-iy)^2}
 &=\int\frac{dx}{x-iy}+2iy\int\frac{dx}{(x-iy)^2}.
\end{aligned}
\tag{A12}
\]

Integration of the smoothed strip and its derivative gives, in the sharp limit, $(I,J_1,J_w)=(\pi,-\pi,-\pi)$ for $c>1$. For $0<c<1$ it gives exactly the three functions in the preceding equation (16). Thus (A11) becomes

\[
-2\pi i\nu\lambda_2t(t-1)\partial_x\mathscr U.
\tag{A13}
\]

The fresh quadrature keeps the additional finite-$\eta$ first-moment term in (A9). At $t=c=\sqrt2$, $K_L=0.3$, $K_R=-0.2$, $D=0.17$, it gives the following imaginary coefficients after division by $\nu\lambda_2$:

| $\eta/\epsilon$ | Smoothed contact coefficient | Absolute error from (A13) |
|---|---:|---:|
| 0.08 | -0.369500803605 | $1.4403\times10^{-3}$ |
| 0.04 | -0.368430959154 | $3.7049\times10^{-4}$ |
| 0.02 | -0.368153762967 | $9.3289\times10^{-5}$ |
| 0.01 | -0.368083838118 | $2.3364\times10^{-5}$ |

The limit predicted independently is $-0.368060473804$. These numbers verify the original smoothed contour moments, not just a polynomial identity between the proposed endpoint coefficients.

For $Y_0$, the same Wick calculation has $q_0=1$, $C+D=(d+a_0)\partial X$, and $q_0A_R-\bar C=0$. Its contact is proportional to $I+J_w=0$ for $c>1$.

### Add the tails, then state the normalization

Direct expansion of $J_B\bar H_B$ and $H_B\bar J_B$ gives equal double-pole coefficients $tq\mathscr U$. The simple coefficients are $-t\partial\mathscr U$ and $tq\partial\mathscr U+t^2\bar\partial\mathscr U$. The double contacts cancel under the opposite holomorphic and antiholomorphic Plemelj signs. The remaining tail contact is $2\pi i\nu\lambda_2t^2\partial_x\mathscr U$. Together with (A13), this leaves the radial coefficient $t$.

For $AA$, the same contour moments with no derivative prefactor give the bare-action contribution $2\pi i\nu\gamma\lambda_1\partial_xF_A$, $F_A=e^{-k\Phi}$. One way to check the coefficient is to set $D=-d\partial\varphi$, $A=-T\partial\varphi=-tD$; then $K=(1-t)D$ and the first-jet expression is $\pi(D-K-A_R)=\pi\gamma(K+K_R)$. Adding the tail coefficient gives $\gamma^2+\gamma=\gamma^2/t$.

Consequently the finite endpoint additions needed to impose the preceding canonical contour-product normalization are indeed

\[
\delta B_{BB}=t(t-1)\lambda_2\mathscr U,\qquad
\delta B_{AA}=-\gamma\lambda_1F_A.
\tag{A14}
\]

This establishes their values in this prescription. It does not make them regulator-independent physical structure constants.

## 4. Endpoint renormalization is checked separately

The Wick pairing $v a_0=2-t$ gives

\[
\mathscr U(0)Y_0(r)
=|r|^{2t-4}:
\bigl[1+rC+\bar r\bar C+r\bar r C\bar C+\cdots\bigr]F:,
\qquad F=e^{k\Phi}.
\tag{A15}
\]

Here the displayed mixed Taylor coefficient is exact; purely holomorphic and antiholomorphic second-order terms have zero angular average. The radial primitives are $\pi(R^{2t-2}-\epsilon^{2t-2})/(t-1)$ and $\pi(R^{2t}-\epsilon^{2t})/t$. The minus sign in the action expansion therefore requires

\[
[\mathscr U]_0=\mathscr U
-\frac{\lambda_0}{t-1}\epsilon^{2t-2}F
-\frac{\lambda_0}{t}\epsilon^{2t}:C\bar C F:.
\tag{A16}
\]

The next angularly surviving term is of fourth Taylor order and has cutoff power $2t+2$. This is why the two displayed subtractions suffice for $\operatorname{Re}t>-1$, away from $t=0,1$.

The same direct pairings give $Y_2Y_0\sim |r|^{2t-2}P$ and $K_2Y_0\sim r^{t-2}\bar r^{t-1}(E+rL+\cdots)$. Their subtraction coefficients are respectively $-\lambda_0\epsilon^{2t}P/t$ and $-\lambda_0\epsilon^{2t}L/t$. Multiplication by $\lambda_2$ and $-t\lambda_2$ reproduces (A4) and (A6). Multiplication of (A16) by $t(t-1)\lambda_2$ verifies the required dressing of (A14):

\[
\delta B_{BB}^{\rm ren}=t(t-1)\lambda_2\mathscr U
-t\lambda_0\lambda_2\epsilon^{2t-2}F
-(t-1)\lambda_0\lambda_2\epsilon^{2t}:C\bar C F:.
\tag{A17}
\]

As a separate finite-part check, the script uses the smooth spectator factor $e^{r+\bar r}$. Its angular average is $\sum_{n\ge0}|r|^{2n}/(n!)^2$. Termwise integration gives a meromorphic expression $\pi\sum_{n\ge0}R^{2t-2+2n}/[(n!)^2(t-1+n)]$. Direct integration with the two subtractions converges to that expression on both sides of zero, at $t=\sqrt2$ and $t=-\sqrt2/5$. This tests the nonconstant descendant part of the subtraction independently of the earlier fermionic spectator example.

## 5. Correct the local collision count

The most singular chiral kernel can be reconstructed by contracting the derivative prefactor with the two screenings:

\[
f(u;w_0,w_2)=
\frac{(u-w_0)^{-1}+(t-1)(u-w_2)^{-1}}
{(u-w_0)(u-w_2)}.
\tag{A18}
\]

For $r=w_0-w_2$, an independent partial-fraction decomposition is

\[
f=\frac1{r(u-w_0)^2}
-\frac{t-1}{r(u-w_2)^2}
+\frac{t-2}{r^2}
\left(\frac1{u-w_0}-\frac1{u-w_2}\right).
\tag{A19}
\]

It displays the proper current-pole faces explicitly. The full lowest fused coefficient is $\nu|w_0-w_2|^{2t-2}f(u)f(\bar v)$, multiplying $F$, as in the preceding report. Its homogeneity $2t-8$ is correct.

Fix the inner point at $v=y$ in the interior of its segment. The relative variables are one real outer-contour displacement $u-y$ and the four real components of $w_0-y,w_2-y$. The measure scales as $\rho^5$. The remaining $dy$ spans a fixed interval; it does not supply a small factor. Thus the raw relative degree is

\[
(2t-8)+5=2t-3.
\tag{A20}
\]

For a cluster pinned at an endpoint, $y-a_\pm=O(\rho)$ supplies the additional factor, giving $2t-2$. These are different local questions.

There is also a useful cancellation. Since $f(-u;-w_0,-w_2)=-f(u;w_0,w_2)$, the upper-minus-lower defect at $v=0$ is odd under simultaneous inversion of the five relative coordinates. The circular masks, the pair exclusion, and a symmetric local patch are invariant. The constant fused-field contribution therefore integrates to zero. The first Taylor jet costs one power of $\rho$, recovering $2t-2$ in that patch. This proves a cancellation, not an absolute bound for the unexpanded kernel.

The new numerical sector check uses 2,991 retained points away from every proper pole. At $t=\sqrt{1.7}$ the absolute primitive kernel is nonzero. Its relative degree is $-0.392319038$; after one Taylor jet it is $+0.607680962$. Reflection cancels the constant term point by point. This gives a concrete check of both the missing measure factor and the extra cancellation. It is not a numerical evaluation of the full charge product.

The analogous raw relative counts, derived from Wick contractions, are:

| Channel near an interior point | Kernel degree | Relative measure dimension | Explicit cutoff factor | Raw integrated power |
|---|---:|---:|---:|---:|
| $J_B\bar J_BY_0Y_2$ | $2t-8$ | 5 | 0 | $2t-3$ |
| $J_B\bar K_2Y_0$, $K_2\bar J_BY_0$ | $2t-6$ | 3 | 0 | $2t-3$ |
| $J_B\bar J_BP$ | $-6$ | 3 | $2t$ | $2t-3$ |
| $J_B\bar L$, $L\bar J_B$ | $-3$ | 1 | $2t$ | $2t-2$ |
| $Y_0Y_2$ | $2t-2$ | 2 | 0 | $2t$ |
| One current and both unlike Yukawas | $2t-5$ | 4 | 0 | $2t-1$ |
| Two currents and two equal Yukawas | $-4$ | 5 | 0 | 1 |
| $K_2\bar K_2$ | 0 | 1 | 0 | 1 |
| Endpoint insertion $\mathscr U Y_0$ | $2t-4$ | 2 | 0 | $2t-2$ |

“Raw” means before an angular or reflection cancellation; ordinary chiral poles still have their specified Cauchy extension. The uncontracted derivative field in the $L$ channel is not included in its kernel degree. In the one-current/two-Yukawa channel the leading fused field is $e^{k\varphi+T\bar\varphi}$ without that derivative; discarding it before an angular projection would overstate the raw degree.

For example, the $J\bar K Y_0$ contraction has two current/screening factors of degree $-1$, one differentiated factor of degree $-1$, and screening factors of degrees $t-1,t-2$. Its degree is $2t-6$. The potential kernel has two cubic poles, giving $-6$. Same-node Yukawas contribute $|w_0-w_2|^2$, explaining the $-4$ degree there. Finally, $(d+a_2)a_2=0$, so $K_2\bar K_2$ has no power singularity.

Every raw power in this table is strictly positive for $\operatorname{Re}t>2$. This sufficient domain avoids having to use the reflection improvement uniformly over all intersecting subtraction sectors.

## 6. A conservative formal continuation argument

### Fix the product before continuing

Use the same outer test function, equal to one near the finite inner segment, and the same distant infrared cutoff on both sides. Remove the auxiliary Cauchy regulator first. Retain circular inner-point and interaction-pair exclusions with $\delta=c\epsilon$, $c$ in a fixed compact subset of $(1,\infty)$. Add the finite contour-product normalization (A14), with its ordinary interaction dressing (A17). No other finite term is specified.

The order-$\lambda_0\lambda_2$ left-hand expansion is then the sum of the four lines of the preceding equation (24), plus the dressed finite product contact. Their signs follow directly from expanding $e^{-S_{\rm int}}$ and $H_B=-t\lambda_2K_2+\lambda_0\lambda_2\epsilon^{2t}L$: the double action insertion has a plus sign, the action correction to either first-order tail has a plus sign, and the composite action insertion is $-\kappa/\pi=+\lambda_0\lambda_2\epsilon^{2t}/(\pi t)$. Equal-Yukawa diagrams and the first-order-tail product must also be retained. This checks the expansion as a whole before taking limits.

### Treat the proper faces explicitly

At separated insertion points, the outer-contour calculation reduces to the following cases.

1. **One current and one interaction.** The distributional double and simple poles in (A5) give coefficients zero and $-t$. They cancel against the corresponding first-order tail when the one-form is used. Opposite-chirality faces obey the same statement.
2. **Both currents and one interaction.** Equations (A8)–(A14) give the complete contour contact for a smooth background near that face. The constant background term has zero odd integral; the first jets give the finite contact; higher jets vanish with the face radius. This is why the calculation applies to a separated additional interaction or a separated spectator, not just to the numerical background used in the check. The $Y_0$ face is zero. The canonically normalized $Y_2$ face has coefficient $t^2\lambda_2\mathscr U$.
3. **Two interactions without both currents.** For unlike Yukawas the degree is $2t$; with one current it is $2t-1$ before angular improvements. A tail with one unlike Yukawa has degree $2t-1$. These vanish in the initial domain. The equal-Yukawa versions are softer. The endpoint dressing has degree $2t-2$ and is locally integrable.
4. **Disjoint current-pole faces.** Their Cauchy distributions act on different relative variables and compose. If the two interaction points approach one another as well, this is one of the larger clusters in the table, not an additional independent proper face.

These cases exhaust the proper subsets of the four-point graph with two currents and two interactions: subsets with no interaction have no mixed oscillator pole; subsets with one interaction are cases 1–2; the interaction pair and the three-point subsets are case 3; disjoint subsets are case 4. The diagrams with one tail or one potential have fewer points and the same listed subsets. The $J\bar L$ terms have only the prescribed one-dimensional boundary-value poles.

Consequently, away from intersections of these faces, the only second-order contribution that remains is the ordinary action dressing of the first-order contact. In particular its unlike-Yukawa part is

\[
-\frac{t^2\lambda_0\lambda_2}{\pi}
\int d^2w\,\langle\mathscr U(a)Y_0(w)\mathcal O\rangle_0,
\tag{A21}
\]

at each endpoint $a$, with the overall charge normalization and phase understood. The equal-$Y_2$ term gives its ordinary dressing too. Its two choices of the contact-producing insertion cancel the $1/2!$ in the action expansion. There is no analogous $Y_0$ contact to dress.

### Bound the remaining intersections in an open domain

Here the corrected table is needed. Partition a compact local patch by the ordering of its relative distances. When a proper current-pole cluster is smallest, subtract the just-calculated local contact and use the Taylor remainder of the other insertions in that cluster. The pole decomposition (A19) exhibits the possible current denominators; outside these proper subclusters their ratios to the cluster scale stay bounded. For two disjoint smallest clusters the operation is applied in the two independent variables. Repeating it leaves either a proper interaction cluster or the full cluster in the table.

Cauchy differentiation has fixed finite order. Under simultaneous rescaling it has the same degree as the original pole. Derivatives of an exclusion mask are accounted for in this distributional operation: they have degree $-1$, with the corresponding boundary integration and Taylor/Stokes factor. They do not authorize ignoring a face, and do not change the displayed net scaling of that face.

On the residual compact shape domains, the remaining integrations are bounded by radial integrals of the form $\int_0^{r_0}\rho^{\omega-1}(1+|\log\rho|)^N\,d\rho$, where $\omega$ is a positive relative degree in the table and $N$ is finite. The cutoff boundary terms tend to zero as $\epsilon^{\operatorname{Re}\omega}$ times at most logarithms. The worst primitive unlike-Yukawa power is $2\operatorname{Re}t-3$. It is positive in the conservative initial domain $\operatorname{Re}t>2$. At a fixed endpoint there is an additional small contour coordinate, so the endpoint estimate is no worse. Keeping $c$ away from 1 avoids a coincident strip edge and inner circle during the first-order contact calculation.

The potentially troublesome nested face has an explicit remainder estimate. Collapse both currents and $Y_2$ at scale $\epsilon$, while $Y_0$ is at distance $R\gg\epsilon$. Before Taylor expansion the first face has degree $-1$; its constant term cancels, and its first jet is the computed contact. Its second-jet remainder therefore has a factor $O(\epsilon)$. The extra endpoint/Yukawa contraction is $O(R^{2\operatorname{Re}t-4})$, and two derivatives cost at most $R^{-2}$. Including its area measure bounds this remainder by a constant times

\[
\epsilon\int_{A\epsilon}^{R_0}R^{2\operatorname{Re}t-5}\,dR=O(\epsilon),
\qquad \operatorname{Re}t>2,
\tag{A21a}
\]

where $A>1$ is fixed by the separation sector. The region $R=O(\epsilon)$ is the primitive full cluster of degree $2t-3$. Interchanging nodes gives the same bound, even though the $Y_0$ first-jet contact vanishes. Single-current proper faces and the regular same-node endpoint collision have weaker bounds. Thus the proof does not treat the extra interaction as a uniformly smooth spectator all the way to its collision with the first contact.

This also bounds the potential and current-counterterm diagrams with their explicit $\epsilon^{2t}$ factors. They vanish in this domain; their role after continuation is fixed by the same single-current Ward subtraction, not by assigning a separate finite coupling. Thus no new primitive mass-squared contact remains in this initial region. The only surviving mass-squared term is (A21) and its equal-node counterpart, together with the prescribed finite-contact dressing. This is a local distributional contour argument, not a claim of absolute convergence of unextended Cauchy poles.

### Continue the independently defined left-hand side

Define the remaining finite part by continuation of the complete regulated left-hand expansion from this domain:

\[
\operatorname{FP}_{>2}^{\rm an}\mathcal C(t)
=\operatorname{AC}_t\left[
\lim_{\epsilon\to0}\mathcal C(t;\epsilon,c\epsilon)
\right]_{\operatorname{Re}t>2}.
\tag{A22}
\]

This definition is made on the left-hand diagrams, including their proper Cauchy faces and (A14); it is not a definition by the desired endpoint formula. In each local distance sector, Taylor subtraction of the smooth spectators leaves radial Mellin integrals. Their exponent is linear in $t$; expanding to successively higher Taylor order continues the integrals meromorphically. Fixed-order Cauchy derivatives and the finite contact coefficients preserve that property. Work on coupling domains away from $t=0$, with a consistent local branch for the momentum vectors, and exclude poles of the continued distributions.

The equality proved in the open initial domain therefore continues. An additional finite term that vanishes throughout that domain but appears only at negative $t$ would not be an analytic continuation of these diagrams. This excludes that particular freedom; it does not exclude choosing a different product normalization at the start.

Any correctly established continuation from $\operatorname{Re}t>1$ with the same prescription must agree with (A22) on their common initial domain and hence by continuation. No proof of the sharper cutoff convergence interval is needed for this conclusion. Conversely, the identity theorem cannot repair an unproved equality in an initial domain; the corrected sector argument above is the required preceding step.

For $\operatorname{Re}t>-1$, $t\ne0,1$, the resulting endpoint representative is

\[
\begin{aligned}
B_{BB}^{\rm can}
&=\lim_{\epsilon\to0}\left[
t^2\lambda_2\mathscr U
-\frac{t^2\lambda_0\lambda_2}{t-1}\epsilon^{2t-2}F
-t\lambda_0\lambda_2\epsilon^{2t}:C\bar C F:
\right]_{S_\epsilon},\\
B_{AA}^{\rm can}&=\gamma^2\lambda_1F_A.
\end{aligned}
\tag{A23}
\]

With $\Delta_a\mathcal O=\mathcal O(a_+)-\mathcal O(a_-)$, the formal contour identities are

\[
\mathcal C_{A\bar A}=\frac{\nu}{2\pi i}\Delta_aB_{AA}^{\rm can}+O(m^3),\qquad
\mathcal C_{B\bar B}=\frac{\nu}{2\pi i}\Delta_aB_{BB}^{\rm can}+O(m^3).
\tag{A24}
\]

The off-diagonal defects vanish at this order. In fact the moving $A$ current has positive integral powers with each Yukawa, while the $B$ current remains regular with the even-root vertex even after contracting its derivative prefactor. The corresponding mixed tails are regular as well. The residual same-order interaction collisions satisfy the preceding bounds. These conclusions do not use a proposed Lie or super-Lie algebra.

Outside the displayed endpoint-subtraction strip, further descendants are required. Exceptional couplings and logarithmic prescriptions are not classified here. Equation (A22), not an arbitrary sharp-cutoff product with only single-current counterterms, is the product prescription being verified.

## 7. Recommended corrections to the preceding report

The following are recommendations only. The reviewed report is unchanged.

For `briefs/referee_sl21_mixed_regulator_report.md`, line 408, replace the assertion that the two contour measures always supply two small dimensions by:

> At a fixed endpoint, both contour coordinates scale and the joint power is $2t-2$. At an interior point the common contour coordinate is not a small relative variable, so the raw relative power is $2t-3$. The leading constant-field defect is odd under reflection and cancels in a symmetric patch; a first Taylor jet then gives $2t-2$. These cancellations must be distinguished from raw homogeneity.

For lines 414–423, label the old counts as total cluster counts at a fixed endpoint, and add the relative-count table in section 5 of this audit. They must not be used as unqualified interior bounds.

For lines 441–459 and the corresponding outcome verdict, use the conservative qualification:

> After evaluating proper current-pole contacts, the corrected raw relative degrees suffice to prove the formal local contour identity in $\operatorname{Re}t>2$. Meromorphic continuation of the complete left-hand perturbation expansion then gives the identity at generic coupling in the stated product prescription. The sharper bare-cutoff convergence claim in $\operatorname{Re}t>1$ requires additional uniform cancellation estimates and is not used here.

For line 295 and the interpretation of the finite matching, add:

> The contact-free canonical bare-current exchange is a normalization condition for this formal contour product. Establishing that it is realized by the interacting nonlocal currents, including all density contacts and their domains, is a separate problem.

No coefficient change to equations (17)–(23) is recommended.

## 8. Checks and remaining work

The new script runs with

~~~text
python3 -u referee_scratch/sl21_mixed_audit_20261023/audit_checks.py
~~~

Its final output is **65 PASS, 0 FAIL**. The meaningful independent checks include the smooth-Cauchy contour calculation, actual Wick contraction degrees, the relative-measure scaling test, the reflection cancellation, the three-pole potential response, and a different nonconstant endpoint spectator integral. The script does not purport to prove the sector decomposition or evaluate the full multidimensional charge integral; those parts of the argument are written in sections 5–6.

The next question about a symmetry algebra still concerns the operators represented by the spatial boundary values of (A23), together with same-chirality products and Serre-type relations. The fermionic endpoint cannot be replaced by a bosonic Cartan factor merely by setting classical fermions to zero. A massive boundary-state calculation or a separately established global operator realization is needed. Neither affine $\mathfrak{sl}(2|1)$ nor affine $\mathfrak{sl}_2$ has been established.

## Integration

No book integration was performed. **Verified:** the finite contact, local bosonization/current residues, the local mass-squared Ward polynomial, and endpoint/current/action subtraction coefficients. **Modified:** the convergence argument must use relative collision dimensions; a conservative initial half-plane and explicit proper-face decomposition replace the unsupported use of total measure dimensions. **Retained with qualification:** the formal meromorphic mixed contour identity in the fixed product prescription. **Not verified as a physical assertion:** canonical realization of the interacting currents, boundary closure, exceptional couplings, or a symmetry algebra. No existing report or book file was changed.
