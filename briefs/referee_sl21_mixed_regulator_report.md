# Root-system sl(2|1) Toda: mixed products, contact normalization, and the mass-squared limit

Date: 2026-10-23. Model: Codex (GPT-6).

Commit reviewed: 745586e6a6f7bca5ed91a1e9cf89d6ef5df387f0.

Task: continue the mixed-charge calculation with an explicit common regulator, checking the preceding ingredients independently. No existing file was edited. No subagent, git mutation, commit, push, or book render was used.

Independent evidence: [regulator_checks.py](../referee_scratch/sl21_mixed_regulator_20261023/regulator_checks.py), [its output](../referee_scratch/sl21_mixed_regulator_20261023/regulator_checks.out), [face_quadrature.py](../referee_scratch/sl21_mixed_regulator_20261023/face_quadrature.py), [its output](../referee_scratch/sl21_mixed_regulator_20261023/face_quadrature.out), and the derivation below. No preceding script was imported. The supplied Bernard–LeClair PDF was inspected again, including its coordinate footnote and equations (2.13)–(2.18).

## Outcome

There are two substantive corrections to the setup before taking the mass-squared limit.

First, the canonical split in the preceding report uses Bernard–LeClair's chirality labels, whereas our Euclidean coordinate is $z=x+i\tau$. The two dictionaries exchange chirality. With our stated Wick rotation, the holomorphic field has the **minus** anchored split. The corresponding braid phases are inverted relative to that report.

Second, a tail-current pole calculation is not the whole Euclidean mixed product. The bare $J\bar J$ product with one action insertion has a finite contact. It survives even in a region where the additional two-Yukawa collisions are integrable. Importing the leading endpoint coefficient without specifying this contact normalization is insufficient.

After fixing these points, one obtains a **formal mixed contour identity through order $m^2$ in a specified analytic finite-part prescription**. In the version matched to the canonical equal-time product, put

\[
\nu=e^{-i\pi/t},\qquad \gamma=\frac{t}{1-t}.
\]

The statements concern generic coupling. The displayed two-term endpoint cutoff representative applies for $\operatorname{Re}t>-1$, $t\ne0,1$; exceptional resonant products are not settled by this report.

For an inner charge segment with endpoints $a_-,a_+$, write $\Delta_a\mathcal O=\mathcal O(a_+)-\mathcal O(a_-)$. The mixed contour defects are

\[
\boxed{
\begin{aligned}
\mathcal C_{A\bar A}&=\frac{\nu}{2\pi i}\Delta_a\mathcal B_{AA}+O(m^3),
&\mathcal B_{AA}&=\gamma^2\lambda_1:e^{-k\cdot\Phi}:,\\
\mathcal C_{B\bar B}&=\frac{\nu}{2\pi i}\Delta_a\mathcal B_{BB}^{\rm ren}+O(m^3),
&\mathcal C_{A\bar B}&=\mathcal C_{B\bar A}=O(m^3),
\end{aligned}
}
\tag{1}
\]

where the diagonal braiding is $\nu^2$, the off-diagonal braiding is $\nu^{-2}$, and

\[
\boxed{
\mathcal B_{BB}^{\rm ren}=\lim_{\epsilon\to0}\left[
t^2\lambda_2\mathscr U
-\frac{t^2\lambda_0\lambda_2}{t-1}\epsilon^{2t-2}F
-t\lambda_0\lambda_2\epsilon^{2t}:C\bar C F:
\right]_{S_\epsilon}.
}
\tag{2}
\]

Here $\mathscr U=i:\chi_2\bar\chi_2e^{(d+p_2)\cdot\Phi}:$, $F=:e^{k\cdot\Phi}:$, $C=p_0\cdot\partial\varphi+j$, and $\bar C=p_0\cdot\bar\partial\bar\varphi+\bar j$. The limit is an insertion in perturbative correlation functions. The definition of the product on the left of (1) includes the finite contact matching in section 5, not just the action and single-current counterterms.

The phase $\nu$ on the right is fixed by the cross-chiral cocycle used here. Rescaling both antichiral charges by $\nu^{-1}$ removes it. No such rescaling is implicit in (1).

Equation (1) is a statement about regulated contour insertions with separated spectators and fixed finite endpoints, followed by analytic continuation. It is **not** an existence theorem for global Hilbert-space charges. Their spatial boundary limits and algebraic closure are still open. In particular, the fermionic endpoint in (2) has not been reduced to a Cartan operator.

| Item | Verdict | Independent evidence |
|---|---|---|
| Local current residues and the preceding mass-squared Ward polynomial | **Correct** | Fresh Wick/Gram calculation and $N-R_P/t=\partial L$ |
| Identifying the previous plus canonical split with a field holomorphic in our $z$ | **Wrong** | Wick rotation and explicit left/right free waves |
| Previous endpoint coefficient as a universal result independent of product contacts | **Wrong** | Finite strip and edge integrals, section 4 |
| Tail-current endpoint contribution itself | **Correct**, with the stated phase | Direct Plemelj calculation and numerical integral |
| Endpoint subtraction and its place in the action/current descent | **Correct** | One-Yukawa collision calculation, section 6 |
| Mixed identity (1) through mass order two | **Correct within the formal, canonically matched analytic prescription below** | Lower-order contact calculation, subtraction of proper collision faces, and convergent-region contour argument |
| A bare sharp-cutoff identity using only the action and current counterterms | **Unclear; not established** | Additional product contacts cannot be omitted |
| Global charge algebra, Cartan closure, and either proposed affine algebra | **Unclear** | Boundary operators, domains, and same-chirality relations remain |

No theorem upgrade or correction to a published paper is proposed.

## 1. Local data checked afresh

The canonical propagators are $\varphi^a(z)\varphi^b(w)\sim-\delta^{ab}\log(z-w)$ and $\chi_0(z)\chi_2(w)\sim(z-w)^{-1}$. Left and right fermions anticommute. Use

\[
\alpha_0=(-1/\sqrt2,i/\sqrt2),\quad
\alpha_1=(\sqrt2,0),\quad
\alpha_2=(-1/\sqrt2,-i/\sqrt2),
\]

\[
p_i=-ib\alpha_i,\quad t=b^2,\quad T=p_0+p_2,\quad
d=-T/t,\quad k=d+T=(1-t)d.
\tag{3}
\]

Thus $p_i^2=0$, $p_0\cdot p_2=-t$, $d\cdot p_i=1$, and $d^2=-2/t$. The ordinary-area action convention is

\[
S_{\rm int}=\frac1\pi\int d^2w\,
(\lambda_0Y_0+\lambda_2Y_2+\lambda_1V_1+\kappa P),\qquad
\kappa=-\lambda_0\lambda_2\epsilon^{2t}/t,
\tag{4}
\]

where $Y_i=i:\chi_i\bar\chi_i e^{p_i\cdot\Phi}:$, $V_1=:e^{-T\cdot\Phi}:$, and $P=:e^{T\cdot\Phi}:$. The free action has coefficient $1/(2\pi)$ multiplying $\partial\Phi\bar\partial\Phi$ and the two fermion kinetic terms. There is no finite independent $P$ coupling. Ultraviolet Wick ordering uses a fixed reference length; $\epsilon$ is expressed in those units. The mass orders of $\lambda_0,\lambda_2,\lambda_1$ are respectively $1,1,2$.

A fresh fermionic contraction gives the positive sign in $Y_2(z)Y_0(z+r)\sim |r|^{2t-2}P(z)$. Point splitting also verifies the regular term $+j$ in the normalized chiral endpoint collision, with $j=:\chi_0\chi_2:$. For calculating poles it is therefore legitimate to use a local auxiliary boson $h$ with $j=i\partial h$. In $X=(\varphi^u,\varphi^v,h)$ the screening momenta are

\[
a_0=(ib/\sqrt2,b/\sqrt2,i),\quad
a_2=(ib/\sqrt2,-b/\sqrt2,-i),
\quad a_0^2=a_2^2=-1,\quad a_0\cdot a_2=1-t.
\tag{5}
\]

The vector $d$ has zero $h$ component. Hereafter a three-dimensional $T$ or $d$ means its extension by zero. The currents and their local tails are

\[
\begin{split}
J_A&=:e^{-d\cdot\varphi}:,&
H_A&=-\gamma\lambda_1:e^{-k\cdot\varphi-T\cdot\bar\varphi}:,\\
J_B&=:(a_0\cdot\partial X)e^{d\cdot X}:,&
H_B^{\rm ren}&=\lim_{\epsilon\to0}
[-t\lambda_2K_2+\lambda_0\lambda_2\epsilon^{2t}L]_{S_\epsilon},
\end{split}
\tag{6}
\]

with $K_2=i:\chi_2\bar\chi_2e^{(d+p_2)\cdot\varphi+p_2\cdot\bar\varphi}:$ and $L=:C e^{k\cdot\varphi+T\cdot\bar\varphi}:$. Antichiral formulas exchange the two chiralities. The bosonization in (5) is used locally; it does not identify the global compactification lattice.

Before the common local cocycle sign, the $J_BY_i$ double-pole coefficients are $1,t-1$, and the simple-pole vectors are $a_0+d=d+a_0$ and $a_0+(t-1)d=-(d+a_2)$. Including the derivative of the double-pole distribution gives conservation coefficients $0,t$. The specified local current-first cocycle sign is $-1$, giving $H_B^{(1)}=-t\lambda_2K_2$.

The mass-squared polynomial was also checked without the preceding script. Write $x_i=p_i\cdot\partial\varphi$, $y_i=p_i\cdot\partial^2\varphi$, and $E=e^{k\cdot\varphi+T\cdot\bar\varphi}$. Then

\[
\begin{split}
R_P/E&=(x_0+j)(x_0+x_2)-\tfrac t2[y_0+y_2+(x_0+x_2)^2],\\
N/E&=\tfrac12(y_0-y_2)+\partial j+\tfrac12(x_0^2-x_2^2)+j(x_0+x_2),\\
N-R_P/t&=\partial L.
\end{split}
\tag{7}
\]

Consequently the action and single-current counterterms in (4) and (6) remain the local Ward prescription. The correction below concerns chirality and products of two currents, not this identity.

## 2. Correct the Minkowski dictionary before using the phases

Our convention $t_M=-i\tau$ and $z=x+i\tau$ gives $z=x-t_M$. For a free wave $\Phi=f(x-t_M)+g(x+t_M)$, the integral of its velocity is $-f+g$, up to the chosen anchor. Hence

\[
\varphi=\frac12\left(\Phi-\int_{-\infty}^x\dot\Phi\right),\qquad
\bar\varphi=\frac12\left(\Phi+\int_{-\infty}^x\dot\Phi\right).
\tag{8}
\]

The plus split used in the preceding report selects $g$, which is antiholomorphic in our coordinate. Bernard–LeClair's printed-page-105 footnote instead gives $z_{\rm BL}=(\tau+ix)/2$. Our $z=2i\bar z_{\rm BL}$. The previous canonical calculation was internally consistent for its anchored variables, but its identification with our holomorphic currents was reversed. This is a convention correction, not a physical determination of $q$.

From $[\Phi(x),\dot\Phi(y)]=4\pi i\delta(x-y)$ one now obtains

\[
[\varphi(x),\varphi(y)]=+\pi i\operatorname{sgn}(x-y),\quad
[\bar\varphi(x),\bar\varphi(y)]=-\pi i\operatorname{sgn}(x-y),\quad
[\varphi(x),\bar\varphi(y)]=+\pi i.
\tag{9}
\]

The compatible cocycle on the additive cover $(d_L,d_R,p_0,p_2)$ is

\[
\varepsilon(n,m)=(-1)^{(n_L+n_R)(m_0+m_2)}
\exp[-i\pi(n_Lm_R-n_Rm_L)/t].
\tag{10}
\]

It follows by rephasing symmetrically ordered Weyl vertices with
$f(n)=\exp[-i\pi(n_L-n_R)(n_0+n_2)/2]$. The local current/Yukawa signs are unchanged. The cross-current exchange matrix is

\[
R=\begin{pmatrix}\nu^2&\nu^{-2}\\\nu^{-2}&\nu^2\end{pmatrix},
\qquad\nu=e^{-i\pi/t}.
\tag{11}
\]

For $x>y$, a general pair of vertices has exchange factor
$\exp[+i\pi(l-r)\cdot(l'+r')]$; the exponent changes to
$-i\pi(l+r)\cdot(l'-r')$ for $x<y$. The new script verifies the cochain and all current/tail exchange phases. The weights of the previously defined $H_{\rm top}=b\alpha_1\cdot\Delta\Phi/(2\pi)$ also reverse: they are $(-2,+2,+2,-2)$ on $(Q_A,Q_B,\bar Q_A,\bar Q_B)$ in this dictionary.

## 3. Define the nested product and calculate its tail contacts

At $\tau=0$ the one-forms give charge densities $J+H$ and $\bar J+\bar H$, each with the factor $(2\pi i)^{-1}$. To avoid coincident spatial endpoints, integrate the outer charge against a test function equal to one on a neighborhood of the inner charge segment $[a_-,a_+]$. Transport the outer contour above and below the inner contour, with time separations $\pm\delta$. Include the small endpoint transport segments when comparing the two contours. The resulting defect is the local contour version of $Q_a\bar Q_b-R_{a\bar b}\bar Q_bQ_a$.

The radial prescription used below is ordered as follows:

1. Poles of the moving outer current are Cauchy distributions, including derivatives of delta functions.
2. Insertions at the inner current and at the resulting endpoint use circular exclusion of radius $\epsilon$. Two interaction points have the same exclusion radius.
3. Take $\delta=c\epsilon$, with a fixed $c>1$. All terms in the perturbation expansion use this choice. Spectator separations and spatial endpoints stay fixed.
4. After the lower-order pole contacts have been specified, take the ultraviolet limit first in $\operatorname{Re}t>1$ and continue the resulting distributions meromorphically. No additional finite terms are added except the explicitly stated canonical matching contact in section 5.

This is a nested prescription. An independently reordered subtraction of the reverse product is not substituted for the lower contour. Infrared regularization is held identical on both sides; the argument concerns the local ultraviolet defect and does not assume an infinite-volume vacuum limit.

The boundary-value identity is

\[
\frac1{(x+i0)^n}-\frac1{(x-i0)^n}
=-\frac{2\pi i(-1)^{n-1}}{(n-1)!}\delta^{(n-1)}(x).
\tag{12}
\]

For the $BB$ product put $v=d+a_2$ and $q=t-1$. In bosonized notation $\mathscr U=e^{v\cdot(X+\bar X)}$, with the corresponding fermionic Wick definition understood. Factoring $\nu\lambda_2$ out of the products gives

| Product | Double-pole coefficient | Simple-pole coefficient along the spatial contour |
|---|---|---|
| $J_B\bar H_B^{(1)}$ | $tq\mathscr U$ in $z^{-2}$ | $-t\partial\mathscr U$ |
| $H_B^{(1)}\bar J_B$ | $tq\mathscr U$ in $\bar z^{-2}$ | $tq\partial\mathscr U+t^2\bar\partial\mathscr U$ |

The second row includes the holomorphic Taylor expansion of the tail at the moving point. Omitting it gives the wrong result. The two double-pole contacts cancel, and the remaining density contact is

\[
2\pi i\nu\lambda_2t^2\partial_x\mathscr U.
\tag{13}
\]

The script integrates the full rational two-tail kernel against a Gaussian test function before taking $\delta\to0$. Richardson extrapolation agrees with (13); the two unextrapolated relative errors are approximately $0.00523$ and $0.00262$. The $AA$ tail calculation similarly gives $2\pi i\nu\gamma^2\lambda_1\partial_xF_A$, with $F_A=e^{-k\cdot\Phi}$. Off-diagonal tail OPEs are regular.

Equation (13) is not yet the full first-order product. There is also a $J_B\bar J_B$ term with one action insertion.

## 4. The finite current-pole face that power counting misses

Put the inner antichiral current at the origin and the $Y_2$ insertion at $w$. The chiral part of the three-insertion Wick kernel is

\[
e^{D u+A_Lw}\left[\frac{C}{u-w}+\frac q{(u-w)^2}\right],
\quad
D=d\cdot\partial X,\quad A_L=a_2\cdot\partial X.
\]

Here constant jets are sufficient to extract the local contact. Write $K_L=v\cdot\partial X$, $K_R=v\cdot\bar\partial\bar X$, and $A_R=a_2\cdot\bar\partial\bar X$. The identities needed below are $C+qD=-K_L$ and $qA_R-\bar C=tK_R$.

Let $\Theta_\delta(w)=1_{|\operatorname{Im}w|<\delta}$. Integrating the difference of the upper and lower outer-current contours gives the distribution

\[
-2\pi i e^{Dw}
[(C+qD)\Theta_\delta+q\partial_w\Theta_\delta].
\tag{14}
\]

The second term is the derivative of the strip boundary. A residue calculation at points strictly inside the strip would miss it. The three finite face integrals are

\[
\begin{split}
I(c)&=\int_{|w|>\epsilon}d^2w\,
\frac{\Theta_\delta(w)}{\bar w^2},\\
J_1(c)&=\int_{|w|>\epsilon}d^2w\,
\frac{\partial_w\Theta_\delta(w)}{\bar w},\\
J_w(c)&=\int_{|w|>\epsilon}d^2w\,
\frac{w\partial_w\Theta_\delta(w)}{\bar w^2},
\qquad c=\delta/\epsilon.
\end{split}
\tag{15}
\]

A fixed distant infrared cutoff changes none of the following limits. For $0<c<1$, direct integration gives

\[
\begin{split}
I(c)&=2[\arcsin c+c\sqrt{1-c^2}],\\
J_1(c)&=-2\arcsin c,\\
J_w(c)&=-2\arcsin c-4c\sqrt{1-c^2}.
\end{split}
\tag{16}
\]

For $c\ge1$, they are respectively $\pi,-\pi,-\pi$. For example, a horizontal slice of $I$ with $|\operatorname{Im}w|<\epsilon$ gives $2\sqrt{\epsilon^2-y^2}/\epsilon^2$ after integrating the two allowed $x$ intervals. Integrating this over $y$ proves the first formula. Evaluating $\partial_w\Theta_\delta$ on the two strip edges proves the other two. Independent Cartesian numerical integrations check all three functions.

These limits do not commute. $I(c)\to0$ when $c\to0$, but $I(c)=\pi$ once $c\ge1$. This is a finite contact with no factor $\epsilon^{2t-2}$. Increasing $\operatorname{Re}t$ cannot remove it.

Multiplying (14) by the antichiral current pole, the action coefficient $-\lambda_2/\pi$, and the cocycle gives the finite bare-current/action-sector contact

\[
-2i\nu\lambda_2q
\{K_L[I(c)-qJ_w(c)]-tK_RJ_1(c)\}\mathscr U.
\tag{17}
\]

For the stipulated $c>1$, this simplifies to

\[
-2\pi i\nu\lambda_2t(t-1)\partial_x\mathscr U.
\tag{18}
\]

Thus the unadjusted radial/Cauchy prescription gives coefficient $t$ after adding (13), rather than $t^2$. For $0<c<1$, (17) is generally not proportional to $\partial_x\mathscr U$ without a further contact adjustment. This is why the ratio of regulators cannot be left implicit.

The separate quadrature script also integrates the unexpanded strip and edge kernel, retaining the exponential background and a Gaussian infrared test function. At $t=c=\sqrt2$, the three finite-cutoff results converge linearly to (18); extrapolating the two smallest cutoffs gives a relative error of about $2.23\times10^{-4}$. This checks the finite contact before replacing the kernel by its local Taylor terms.

The same calculation for a $Y_0$ insertion gives a coefficient proportional to $I+J_w$, which is zero for $c>1$. It agrees with the vanishing fixed-current conservation coefficient for $Y_0$. For the $AA$ product and the even-root insertion, the bare-current/action contact adds $2\pi i\nu\gamma\lambda_1\partial_xF_A$. Consequently the radial $AA$ coefficient is $\gamma^2+\gamma=\gamma^2/t$.

This calculation is not a proposed correction to Bernard–LeClair. Their page 106 defines a braided time-ordered product and explicitly assumes: “We suppose that the same braiding relations hold” when current components are replaced by tails. Their result is stated “to lowest non-trivial order in perturbation theory”. The radial/Cauchy distribution extension just calculated must be matched to that product convention before importing its endpoint coefficient. [Bernard–LeClair](https://doi.org/10.1007/BF02099173), equations (2.15)–(2.18).

## 5. State the finite product normalization explicitly

The anchored canonical fields in (8), with spatial derivatives for the chiral derivative factors, give no delta-function contact in the bare $J_a\bar J_b$ exchange. I use this as an additional **product normalization condition**. It fixes the finite subtraction of (18) and its $AA$ counterpart. This is separate from normal ordering an individual current and from the action counterterm.

At leading order, the added endpoint zero-forms are

\[
\delta\mathcal B_{BB}=t(t-1)\lambda_2\mathscr U,\qquad
\delta\mathcal B_{AA}=-\gamma\lambda_1F_A.
\tag{19}
\]

They convert the radial coefficients $(t,\gamma^2/t)$ to the canonically matched coefficients $(t^2,\gamma^2)$. The finite counterterm in (19) is determined by the bare-current matching condition, not by choosing an affine algebra.

At mass order two this product counterterm must itself be renormalized. Omitting its Yukawa dressing would undo the matching. In particular,

\[
\delta\mathcal B_{BB}^{\rm ren}
=\lim_{\epsilon\to0}\left[
t(t-1)\lambda_2\mathscr U
-t\lambda_0\lambda_2\epsilon^{2t-2}F
-(t-1)\lambda_0\lambda_2\epsilon^{2t}:C\bar C F:
\right]_{S_\epsilon}.
\tag{20}
\]

The derivation of this dressing is given next. Its $AA$ counterpart starts dressing only at mass order three.

This matching defines a formal operator-product prescription compatible with the chosen canonical phase algebra. It is not a proof that the full off-critical nonlocal operators have been constructed on a physical Hilbert space. That distinction remains necessary even after the local product has been normalized.

## 6. Action, currents, and endpoint form one subtraction chain

Consider one $Y_0$ insertion in the ultraviolet Wick products below. Let $[\mathcal O]_0$ denote the corresponding renormalized insertion through first order in $\lambda_0$. Direct expansion gives

\[
\begin{aligned}
[Y_2]_0&=Y_2-\frac{\lambda_0}{t}\epsilon^{2t}P,\\
[K_2]_0&=K_2-\frac{\lambda_0}{t}\epsilon^{2t}L,\\
[\bar K_2]_0&=\bar K_2-\frac{\lambda_0}{t}\epsilon^{2t}\bar L,\\
[\mathscr U]_0&=\mathscr U
-\frac{\lambda_0}{t-1}\epsilon^{2t-2}F
-\frac{\lambda_0}{t}\epsilon^{2t}:C\bar C F:.
\end{aligned}
\tag{21}
\]

Each equality describes the cutoff representative inserted into correlation functions, not equality of bare fields. Their divergent integrals are generated respectively by

\[
\begin{array}{c|c|c}
\text{collision with }Y_0(z+r)&\text{leading kernel}&\text{rotationally surviving coefficients}\\ \hline
Y_2&|r|^{2t-2}&P\\
K_2&r^{t-2}\bar r^{t-1}&rL\\
\bar K_2&r^{t-1}\bar r^{t-2}&\bar r\bar L\\
\mathscr U&|r|^{2t-4}&F+r\bar r:C\bar C F:.
\end{array}
\tag{22}
\]

The leading angular harmonic in $K_2Y_0$ integrates to zero. All signs in (22) follow from fermionic Wick reordering together with the specified local cocycle and the orientation $w=z+r$. The radial denominators are

\[
\int_{\epsilon<|r|<R}d^2r\,|r|^{2t-4}
=\frac\pi{t-1}(R^{2t-2}-\epsilon^{2t-2}),\qquad
\int_{\epsilon<|r|<R}d^2r\,|r|^{2t-2}
=\frac\pi t(R^{2t}-\epsilon^{2t}).
\tag{23}
\]

Multiplying the first line of (21) by $\lambda_2$ reproduces the action composite term in (4). Multiplying its next two lines by $-t\lambda_2$ reproduces both renormalized current tails. Multiplying the last line by $t^2\lambda_2$ gives (2); multiplying it by $t(t-1)\lambda_2$ gives (20). Thus the finite product matching has the same ultraviolet dressing as the endpoint it multiplies.

An independent fermionic spectator check uses

\[
\frac{\langle\chi_2(0)\chi_0(r)\chi_0(a)\chi_2(b)\rangle}
{(-1/r)\langle\chi_0(a)\chi_2(b)\rangle}
=1-\frac{r(a-b)}{a(r-b)}.
\]

For $a=2,b=3$ and the corresponding antichiral spectator, its exact angular average is $1+|r|^2/[36(1-|r|^2/9)]$. Integrating this rational function checks both endpoint subtractions and their analytic continuation at $t=\sqrt2$ and $t=-\sqrt2/5$. In particular, it tests a nonzero descendant matrix element, not only the leading scalar contraction.

## 7. Account for the full mass-squared expansion

The mixed $\lambda_0\lambda_2$ part of the unnormalized charge-density product has the following Gaussian terms before the contour defect is taken:

\[
\begin{split}
&\frac{\lambda_0\lambda_2}{\pi^2}
\int_{|w_0-w_2|>\epsilon}d^2w_0d^2w_2\,
J_B\bar J_BY_0(w_0)Y_2(w_2),\\
&\frac{t\lambda_0\lambda_2}{\pi}\int d^2w\,
(J_B\bar K_2+K_2\bar J_B)Y_0(w),\\
&\lambda_0\lambda_2\epsilon^{2t}(J_B\bar L+L\bar J_B),\\
&\frac{\lambda_0\lambda_2\epsilon^{2t}}{\pi t}
\int d^2w\,J_B\bar J_BP(w).
\end{split}
\tag{24}
\]

All insertions carry the inner-point exclusions and Cauchy prescriptions of section 3. The last line has a plus sign because it comes from $-\kappa P/\pi$ in the Euclidean weight. The two current-tail lines are not optional. In the canonical version one must also include the dressing of the finite product contact (19). Terms with two equal Yukawas and the regular $H_B^{(1)}\bar H_B^{(1)}$ product complete the order-two list. An even-root action insertion is regular with $J_B$ and contributes no $BB$ contact at this order.

The most singular coefficient of the simultaneous two-current/two-Yukawa collision can be generated directly, without assuming a symmetry algebra:

\[
\begin{split}
\mathcal K_{02}(u,v;w_0,w_2)
&=\nu|w_0-w_2|^{2t-2}
f(u;w_0,w_2)f(\bar v;\bar w_0,\bar w_2),\\
f(u;w_0,w_2)
&=\frac{(u-w_0)^{-1}+(t-1)(u-w_2)^{-1}}
{(u-w_0)(u-w_2)}.
\end{split}
\tag{25}
\]

It multiplies the lowest fused field $F$. Additional field jets improve the short-distance degree. The kernel in (25) scales as $\rho^{2t-8}$. The two contour measures and two area measures supply $\rho^6$, giving the joint degree $2t-2$. This exponent was checked from the explicit kernel, in addition to a separate free-dimension calculation.

One must first subtract or retain the already calculated proper current-pole faces. The finite contact (17) shows why a degree estimate for the simultaneous collision alone would not be a proof. After those proper faces have been accounted for, the remaining local degrees are:

| Remaining collision or counterterm channel | Integrated small-scale power |
|---|---|
| $J_B\bar J_BY_0Y_2$ joint collision | $2t-2$ |
| $J_B\bar K_2Y_0$ and $K_2\bar J_BY_0$ | $2t-2$ |
| $J_B\bar J_BP$ with its action coefficient $\epsilon^{2t}$ | $2t-2$ |
| $J_B\bar L$ or $L\bar J_B$ with $\epsilon^{2t}$ | $2t-1$ |
| $Y_0Y_2$ subcollision | $2t$ |
| One current and both unlike Yukawas | $2t$ |
| Joint collision with two equal Yukawas | $2$ |
| $K_2\bar K_2$ | $2$ |

The $L$ channel leaves at least one derivative field uncontracted, which explains its slightly better degree. There is no omitted bare $J_B\bar J_B$ oscillator pole. The other proper subcollisions are the single-current/Yukawa poles and the endpoint mixing already handled in sections 3–6. This enumerates the ultraviolet faces at this mass order.

### Contour deformation in the convergent region

Work first with $\operatorname{Re}t>1$, away from $t=1$. Keep the spectator insertions at fixed nonzero separation and use the same infrared regulator on both sides.

On a configuration with separated insertion points, move the outer current contour through the interaction integrals. For an ordinary single-current face, the derivative-of-delta term of its double pole combines with its simple pole. The resulting coefficients are precisely $0$ for $Y_0$ and $-t$ for $Y_2$, giving the tails in (6). Thus these faces cancel against the tail terms in the charge. This step uses the distributional coefficient, not just the current-first simple residue.

For a face containing both currents and one interaction, the remaining contact is the explicit calculation of sections 3–5. In the radial prescription its $Y_0$ coefficient vanishes for $c>1$, and its $Y_2$ coefficient is $t$. In the canonical matching prescription the latter is $t^2$. The product matching has no $Y_0$ term. Consequently a separated additional $Y_0$ insertion simply dresses the $Y_2$ endpoint. This produces the ordinary action expansion

\[
-\frac{\lambda_0}{\pi}\int d^2w\,
\langle\mathscr U(a)Y_0(w)\mathcal O\rangle_0
\tag{26}
\]

with the same endpoint coefficient as at first order. A further $Y_2$ insertion gives its ordinary spectator dressing as well. The factor $1/2!$ for two identical interactions is canceled by the two choices of the insertion that supplies the first-order contact. This dressing is retained; it is not discarded as a vanishing local counterterm.

The only obstruction to interchanging this contour calculation and the remaining interaction integral is a collision of their faces. Their kernels are exactly the channels listed above. After proper lower-order contacts have been removed, their angular integrations are bounded in each sector of relative separations. The radial boundary contributions have the displayed positive powers for $\operatorname{Re}t>1$. They therefore vanish. Derivatives of the exclusion boundary have the same net degree: the inverse-radius factor is compensated by the boundary measure and the corresponding Stokes primitive. For $c$ in a compact subset of $(1,\infty)$ these estimates are uniform.

In particular, the last two lines of (24) and the remaining unlike-Yukawa collision have no separate finite remainder in this open region. The local relation (7) supplies the consistent single-current subtraction when those same faces are continued outside the convergence region. This combines the action composite term, both current counterterms, and the double-Yukawa contribution; it does not test any one of them in isolation.

This proves the order-two local contour identity in the open convergence region, including the ordinary dressing (26). It is a contour derivation supplemented by estimates of all collision faces, not a numerical evaluation of a six-dimensional charge integral. The script checks the pole coefficients, finite face integrals, and the collision kernels and degrees used in the argument.

### Continue the identity, not the bare cutoff factors

Define the remaining analytic finite part independently on the left-hand perturbation expansion (24):

\[
\operatorname{FP}^{\rm an}\mathcal C(t)
:=\operatorname{AC}_{t}
\left[\lim_{\epsilon\to0}\mathcal C(t;\epsilon,c\epsilon)
\right]_{\operatorname{Re}t>1}.
\tag{27}
\]

The lower-order Cauchy faces and the finite matching condition are fixed before this operation. In compact local coordinate patches the remaining radial integrals have Mellin continuations, so their coefficients continue meromorphically in $t$. The equality proved in the open region therefore continues as an identity of these formal distributions at generic values away from their poles. A finite term cannot be added only on the negative-$t$ side while preserving this analytic prescription. This argument is not a separate classification of exceptional resonant operator products.

For $\operatorname{Re}t>-1$, $t\ne0,1$, the endpoint's cutoff representative has exactly the two subtractions in (21); the next rotationally invariant term is $O(\epsilon^{2\operatorname{Re}t+2})$. The resulting radial zero-form is

\[
\mathcal B_{BB}^{\rm rad}
=\lim_{\epsilon\to0}\left[
t\lambda_2\mathscr U
-\frac{t\lambda_0\lambda_2}{t-1}\epsilon^{2t-2}F
-\lambda_0\lambda_2\epsilon^{2t}:C\bar C F:
\right]_{S_\epsilon}.
\tag{28}
\]

Adding (20) gives (2). Likewise $\mathcal B_{AA}^{\rm rad}=\gamma^2\lambda_1F_A/t$, and adding (19) gives (1). The off-diagonal products remain zero at this order: their relevant single-interaction OPEs are regular, and the remaining simultaneous Yukawa faces vanish in the same convergence argument.

Equations (27)–(28) specify an analytic product scheme. They are not a claim that an arbitrary sharp-cutoff product becomes finite using only the action and single-current counterterms. Power terms in the composite two-current product also have to be extended according to (27). At $t=0,1$ logarithmic prescriptions are required; for $\operatorname{Re}t\le-1$ more endpoint descendants need subtraction. Those cases are not settled here.

## 8. What remains before there is a symmetry algebra

The finite product contact is now explicit, and the formal mixed contour relation has been extended through order mass squared in a stated scheme. This removes a specific local calculation gap. It does not identify the infinite-volume endpoint operators.

To turn (1) into a quantum-group presentation one still needs a global definition of the currents and charge domains, limits of $\mathcal B_{AA}$ and $\mathcal B_{BB}^{\rm ren}$ in the chosen massive boundary sectors, and the remaining same-chirality relations. The fermionic endpoint in (2) cannot be replaced by a bosonic Cartan exponential as a UV operator. Its massive boundary matrix elements may involve condensates and nonanalytic mass dependence; the ultraviolet expansion calculated here does not determine them.

Neither the diagonal braid matrix nor the fact that the radial and canonically matched coefficients differ by a common normalization identifies an affine algebra. No affine $\mathfrak{sl}(2|1)$ or affine $\mathfrak{sl}_2$ algebra is assumed or established.

## 9. Checks and recommended corrections

Run

~~~text
python3 -u referee_scratch/sl21_mixed_regulator_20261023/regulator_checks.py
python3 -u referee_scratch/sl21_mixed_regulator_20261023/face_quadrature.py
~~~

All **73 checks pass**: 68 in the symbolic/distributional script and 5 in the unexpanded face quadrature. The checks include the coordinate dictionary, corrected phase cochain, both current residues, direct tail-current contacts, numerical Plemelj integrals, all three finite face integrals, the action-sector contact, fermionic endpoint mixing, the mass-squared Ward polynomial, the complete collision-degree list, the explicit joint kernel's scaling, and analytic continuation of a nonconstant spectator integral. The joint-kernel tests check homogeneity; they are not represented as a full numerical charge-product calculation.

The first quadrature run included a coarse cutoff $\epsilon=0.01$, whose roughly 15 percent finite-cutoff error exceeded the preset 9 percent individual-point threshold. The [coarse output](../referee_scratch/sl21_mixed_regulator_20261023/face_quadrature_coarse.out) is retained. Halving the cutoff sequence, without changing the predicted contact or accuracy thresholds, gives the stated convergence and extrapolation check.

For `briefs/referee_sl21_mixed_charges_report.md`, lines 219–302, replace the canonical dictionary by (8)–(11), invert the Weyl rephasing and cross-cocycle exponent, and reverse the stated topological weights. Suggested explanatory text:

> With $z=x+i\tau$ and $t_M=-i\tau$, the holomorphic anchored field has the minus velocity integral. Bernard–LeClair's plus split corresponds to the opposite chirality in this coordinate convention. The resulting exchange matrix is obtained by inverting the previous phases. This is a coordinate correction, not a physical choice of the sign of $q$.

For that report's discussion following equation (19), lines 319–333, add:

> These zero-forms give the tail-current descent contribution. A regulated Euclidean product also contains a finite bare-current/action contact. Its value depends on the relative current and interaction cutoffs. In the nested radial prescription with $\delta/\epsilon>1$, the diagonal coefficients are smaller by $1/t$; matching the canonical bare-current product adds an explicit finite contact and its Yukawa dressing. The full product convention must be stated before identifying the displayed zero-forms with charge commutators.

These are recommendations only. The earlier report is unchanged. There is no proposed published-source correction.

## Integration

No book integration was performed. **Verified independently:** local bosonization and current residues, the mass-squared Ward polynomial, and endpoint mixing. **Corrected:** the canonical chirality dictionary inherited from the Bernard–LeClair coordinate. **Extended:** the finite current-pole face calculation, explicit product normalization, its mass-squared dressing, and the formal mixed contour identity in an analytic prescription. **Rejected:** using only the simultaneous-collision power count to infer absence of all mixed-product contacts; equation (17) is a concrete counterexample. **Still unresolved:** global charge existence, physical boundary operators, exceptional couplings, and the defining same-chirality algebraic relations.

The input hash manifest was checked after the work. All recorded existing files and the tracked working tree are unchanged.
