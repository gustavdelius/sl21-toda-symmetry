# Root-system sl(2|1) Toda: explicit three-current/Yukawa kernels

Date: 2026-10-24. Model: Codex (GPT-6).

Commit reviewed: `745586e6a6f7bca5ed91a1e9cf89d6ef5df387f0`.

Task: independently calculate the first-order-in-$\lambda_2$ kernels needed for the three-contour defect $(BA-s^{-1}AB)\bar B-\bar B(BA-s^{-1}AB)$. Here $\lambda_2$ is the node-2 Yukawa coupling, not a second perturbative order. `AGENTS.md` was read. No pre-existing file was edited, no earlier script was imported or rerun, no subagent was spawned, and no git mutation, commit, push, or book render was performed.

New evidence is [kernel_check.py](../referee_scratch/sl21_triple_kernel_20261024/kernel_check.py), its [output](../referee_scratch/sl21_triple_kernel_20261024/kernel_check.out), and the derivation below. The final output is **23 PASS, 0 FAIL**. This report does not make a literature or novelty claim.

## Outcome

The actual free-field kernels have a common first-order convergence domain after the proper mixed-current Cauchy contact is treated: $\operatorname{Re}\rho<0$, where $\rho=2/t$. The simultaneous three-current/one-Yukawa face has positive relative degree $1-\rho$ in the interior and $2-\rho$ at an endpoint. The new current does not introduce a marginal primitive contact in that domain.

The mixed contact with the additional $A$ current present is calculated directly below. The action sector gives coefficient $-t(t-1)$, the two tails give $t^2$, and the unadjusted radial sum gives $t$, each multiplying the derivative of the actual product $A\mathscr U$. Thus the finite addition $t(t-1)$ propagates to this three-current insertion. It is not necessary to guess its coefficient from an associative algebra.

These facts support the projected three-contour identity in section 6, defined there with an $A$ test function supported strictly inside the $\bar B$ segment. They do not establish a general associative multiplication of arbitrary charge densities, the compatibility of this triple prescription with all previous four-fold continued products, or physical infinite-volume boundary values.

| Item | Verdict | Independent evidence |
|---|---|---|
| Exact order-$\lambda_2$ three-current kernels | **Correct**, equations (T2)–(T4) | Wick contractions and explicit cocycle multiplication |
| Mixed contact with a separated $A$ insertion | **Correct** | Exact double/simple-pole reduction, full first-jet identity, and unexpanded contour quadrature |
| New primitive three-current contact required in $\operatorname{Re}\rho<0$ | **Wrong** in the stated analytic prescription | Positive relative degrees after the proper faces are resolved |
| Canonical finite matching propagates with coefficient $t(t-1)$ | **Correct** for this projected insertion | Explicit $A\mathscr U$ derivative and integrable endpoint collision |
| Projected first-order three-contour identity (T13) | **Correct within the prescription of section 6** | Local contact calculation, proper-face analysis, and endpoint phases |
| An arbitrary-test associative charge algebra or massive boundary closure | **Unclear** | Not calculated here |

The “wrong” entry rejects a possible inference rather than correcting a sentence in the input reports. No book replacement text or theorem upgrade is proposed.

## 1. Independent local input and ordering convention

Use $X(z)X(w)\sim-\log(z-w)$, with bilinear products

\[
a_0^2=a_2^2=-1,\quad a_0a_2=1-t,\quad
d=-\frac{a_0+a_2}{t},\quad d^2=-\rho,\quad a_i d=1,
\qquad \rho=\frac2t,\quad \beta=t-1.
\tag{T1}
\]

The currents are $J_A=:e^{-dX}:$, $J_B=:C e^{dX}:$, $C=a_0\partial X$, and the first tail is $H_B=-t\lambda_2K_2$, $K_2=:e^{(d+a_2)X+a_2\bar X}:$. A contraction of $C(z)$ with an exponential of momentum $q$ is $-a_0q/(z-w)$. This gives $+1/(z-a)$ from $J_A(a)$ and $\beta/(z-w)$ from $Y_2(w)$.

For the reference order $A(a)B(z)\bar B(u)$, put the action insertion last. The cocycle is

\[
\varepsilon(n,m)=(-1)^{(n_L+n_R)(m_0+m_2)}
e^{-i\pi(n_Lm_R-n_Rm_L)/t}.
\]

The labels of $A,B,\bar B,Y_2$ are $(-1,0,0,0)$, $(1,0,0,0)$, $(0,1,0,0)$, and $(0,0,0,1)$. Multiplying the three successive cocycles gives $-1$ for $AB\bar B Y_2$. Together with the action weight $-\lambda_2/\pi$, its coefficient is therefore $+\lambda_2/\pi$. The cocycles for $AB\bar K_2$ and $AK_2\bar B$ are respectively $+1$ and $-1$; their coefficients after multiplying the tail are $-t\lambda_2$ and $+t\lambda_2$.

These signs are useful checks: the mixed pair by itself carries $\nu=e^{-i\pi/t}$, while $A\mathscr U$ in this order has cocycle $-\nu^{-1}$. Their product is $-1$.

## 2. Exact Wick kernels before integration

All fields inside each displayed normal-ordering symbol are normally ordered together. For the reference order the bare-current/action kernel, including its action sign, is

\[
\begin{split}
\mathcal K_{\rm act}={\lambda_2\over\pi}
{(a-z)^{-\rho}(a-w)\over(z-w)(\bar u-\bar w)}
\left[C(z)+{1\over z-a}+{\beta\over z-w}\right]
\left[\bar C(u)+{\beta\over\bar u-\bar w}\right]\\
\times :e^{-dX(a)+dX(z)+a_2X(w)+d\bar X(u)+a_2\bar X(w)}:.
\end{split}
\tag{T2}
\]

The two tail kernels are

\[
\begin{split}
\mathcal K_{\bar H}=-t\lambda_2
{(a-z)^{-\rho}(a-u)\over z-u}
\left[C(z)+{1\over z-a}+{\beta\over z-u}\right]\\
\times:e^{-dX(a)+dX(z)+a_2X(u)+(d+a_2)\bar X(u)}:,
\end{split}
\tag{T3}
\]

\[
\begin{split}
\mathcal K_H=+t\lambda_2
{(a-z)^{1-\rho}\over\bar z-\bar u}
\left[\bar C(u)-{\beta\over\bar z-\bar u}\right]\\
\times:e^{-dX(a)+(d+a_2)X(z)+a_2\bar X(z)+d\bar X(u)}:.
\end{split}
\tag{T4}
\]

The minus sign inside the final bracket comes from differentiating at $u$, whereas the antichiral vertex factor in this insertion order is $(\bar z-\bar u)^{-1}$. Missing that orientation would give opposite double-pole coefficients for the two tails and spoil the contact calculation.

There are no other first-order terms: $H_A$ begins with the order-two even-root coupling, and a product of two $B$ tails is quadratic in Yukawa couplings.

For the lowest fused coefficient put $C=\bar C=0$, $r=z-a$, and $p=z-w$. The chiral rational factor in (T2) reduces exactly to

\[
{p-r\over p}\left({1\over r}+{\beta\over p}\right)
={1\over r}+{\beta-1\over p}-{\beta r\over p^2}.
\tag{T5}
\]

Thus the scalar kernel has homogeneous degree $-\rho-3$. Its fused field is $:e^{a_2X+(d+a_2)\bar X}:$; it is not the mixed endpoint $\mathscr U$ alone.

## 3. The mixed face with the actual additional current

First hold $a\ne u$ fixed and set $u=0$. For the zero-field-jet part define

\[
V(w)=(a-w)^{1-\rho}.
\]

It is exactly the Wick factor in $A(a)\mathscr U(w)$, not an arbitrary imposed spectator function. In the chiral factor of (T2), the coefficient of $(z-w)^{-2}$ is $\beta V(w)$ and the simple residue is $-V'(w)$:

\[
\left[(a-z)^{-\rho}(a-w)
\left({1\over(z-a)(z-w)}+{\beta\over(z-w)^2}\right)\right]_{z=w}
\longmapsto
\frac{\beta V(w)}{(z-w)^2}-\frac{V'(w)}{z-w}.
\tag{T6}
\]

Equation (T6) denotes the two pole coefficients, not equality away from the pole. Differentiating $(a-z)^{-\rho}$ gives the factor $\beta\rho-1=1-\rho$. This is the explicit Wick calculation responsible for the derivative of the dressed endpoint.

The same conclusion holds with arbitrary smooth field jets. Let $D=d\partial X$, $A_L=a_2\partial X$ near zero. The root identity is $C+tD+A_L=0$. The additional current shifts the three first-jet inputs to

\[
D_{\rm eff}=D+\rho/a,\qquad
C_{\rm eff}=C-1/a,\qquad
A_{L,\rm eff}=A_L-1/a.
\]

Consequently

\[
C_{\rm eff}+\beta D_{\rm eff}
=-(D_{\rm eff}+A_{L,\rm eff}).
\tag{T7}
\]

This verifies the full first-jet Ward combination rather than only the scalar kernel. The antiholomorphic identity is unchanged. Derivatives of the smooth normal product and of separated spectators are included in these jets.

Let the moving $B$ contour have upper/lower distances $\delta_+,\delta_->\epsilon$ from the inner current and let $\Theta$ be its intervening strip. The Cauchy contour difference gives

\[
-2\pi i\{-V'(w)\Theta+\beta V(w)\partial_w\Theta\}.
\tag{T8}
\]

The second term is retained. With the circular inner exclusion $|w|>\epsilon$, the required limits are

\[
\int\frac{\Theta\,d^2w}{\bar w^2}=\pi,\qquad
\int\frac{\partial_w\Theta\,d^2w}{\bar w}=-\pi,\qquad
\int\frac{w\partial_w\Theta\,d^2w}{\bar w^2}=-\pi.
\]

For example, the horizontal integral of $1/\bar w^2$ outside the inner circle is $2\sqrt{\epsilon^2-y^2}/\epsilon^2$ for $|y|<\epsilon$ and zero otherwise. Both strip boundaries lie outside this support. Each edge of either first-jet integral contributes $-\pi/2$, independently of its distance. Hence asymmetric upper/lower contour heights do not alter these constants.

In the reference ordering, (T8) and the action coefficient give $+2\pi i\lambda_2t\beta V'(0)$. Restoring the common phase $\nu$ and writing the answer before multiplication by $A$ gives the familiar action coefficient $-t\beta$, now derived with $A$ genuinely inserted.

The tail kernels give $-2\pi i\lambda_2t^2V'(0)$ in the reference ordering. One can see the coefficient directly: the first tail's simple pole is $t\lambda_2V'(0)/z$, whereas the second tail contains $-t\lambda_2\beta V'(0)z/\bar z^2$. Their integrated contour jumps are respectively $-2\pi it\lambda_2V'(0)$ and $-2\pi it\lambda_2\beta V'(0)$. Their double-pole jumps cancel. Thus the total is

\[
-2\pi i\lambda_2tV'(0)
\tag{T9}
\]

in the reference ordering, equivalently $2\pi i\nu\lambda_2t\,A\partial_x\mathscr U$ with its actual Wick and cocycle factors. Adding the prescribed canonical normalization changes $t$ to $t^2$. There is no additional coefficient left undetermined by this local calculation.

## 4. Unexpanded contour check

The new script integrates (T2) after its exact $z$-pole reduction and integrates both exact tails (T3)–(T4), without expanding $V$. It uses $t=-\sqrt2$, $a=3$, a fixed outer disk of radius one, $\delta=\sqrt3\epsilon$, and a circular inner hole. The action integral keeps the strip derivative as actual edge integrals. The tails are integrated along the finite upper and lower horizontal segments.

The expected imaginary coefficients with $\lambda_2=1$ in the reference order are

\[
\text{action}=-244.904962996413,\quad
\text{tails}=143.462005830818,\quad
\text{sum}=-101.442957165595.
\]

They follow from $V'(0)=-(1-\rho)3^{-\rho}$ and respectively $2\pi t(t-1)V'(0)$, $-2\pi t^2V'(0)$, and $-2\pi tV'(0)$. At $\epsilon=0.000625$ the respective values are $-244.287119689297$, $143.034670829652$, and $-101.252448859645$. Combining the last two cutoffs by linear extrapolation gives absolute errors $7.31\times10^{-5}$, $7.19\times10^{-5}$, and $1.21\times10^{-6}$. The saved output gives all finite-cutoff values and their convergence. This is a direct integral test of the unexpanded additional-current kernel, not merely a check that two proposed coefficients add correctly.

Initial symbolic-output files preserve two checking-code issues: a fractional-power identity needed simplification, and differentiation before cancelling a removable singularity produced an invalid substitution. A coarse numerical run is also retained; its extrapolation missed the declared tolerance, and smaller cutoffs were added. No assertion was weakened to obtain the final result.

## 5. Simultaneous faces and propagation of finite matching

At an interior simultaneous collision there are four relative integration dimensions: two other current positions and the two coordinates of the Yukawa insertion. The common contour coordinate does not scale. Equation (T5) therefore gives relative degree

\[
(-\rho-3)+4=1-\rho.
\tag{T10}
\]

A cluster pinned to an endpoint has one more small coordinate and degree $2-\rho$. The two tail kernels have degree $-\rho-1$ and two relative current coordinates, giving the same interior degree $1-\rho$.

The proper faces must be considered before using (T10):

* $AB$ has leading coefficient $r^{-\rho-1}$ and relative degree $-\rho$.
* $A/Y_2$ has a linear zero, and $A/\bar H_B$ is regular. They have no off-diagonal contour contact.
* $AB/Y_2$ with $\bar B$ separated has degree $2-\rho$.
* $AB\bar B$ with $Y_2$ separated has degree $1-\rho$. The $A\bar B/Y_2$ face with $B$ separated has degree $2$, after its proper antichiral Cauchy pole is treated.
* The $B/Y_2$ and $\bar B/Y_2$ poles are precisely the prescribed Cauchy distributions whose conservation tails appear in (T3)–(T4).
* The $B\bar B/Y_2$ face is the finite first-jet contact in section 3. After it collapses, $A\partial\mathscr U$ has coefficient $r^{-\rho}$ and relative degree $1-\rho$.

All the additional relative degrees are strictly positive in $\operatorname{Re}\rho<0$. There is therefore a common first-order domain after treating the proper mixed contact. The absence of a common convergence domain for an isolated same-chirality $BB$ tail in previous four-current work does not apply to this $AB\bar B$ calculation.

For completeness, a mixed face of size $\epsilon$ nested inside separation $R$ from $A$ has first unretained Taylor remainder bounded by $O(\epsilon R^{-1-\operatorname{Re}\rho})$. Its integral is bounded by

\[
\epsilon\int_{A\epsilon}^{R_0}R^{-1-\operatorname{Re}\rho}\,dR=O(\epsilon)
\qquad(\operatorname{Re}\rho<0),
\tag{T11}
\]

for fixed $A>1$. The complementary region $R=O(\epsilon)$ is bounded by the primitive degree $\epsilon^{1-\operatorname{Re}\rho}$. Smooth spectator derivatives and contour/test-function derivatives have no worse degree. This controls the potentially nonuniform spectator expansion near the simultaneous face.

At a mixed endpoint $l$, the actual product has

\[
A(a)\mathscr U(l)\ \propto\ (a-l)^{1-\rho}
:e^{-dX(a)+(d+a_2)(X+\bar X)(l)}:.
\tag{T12}
\]

Its endpoint integral has degree $2-\rho$. Thus the finite normalization addition $\nu t(t-1)\lambda_2\mathscr U(l)/(2\pi i)$ can be multiplied by the integrated $A$ insertion in this domain without a new endpoint subtraction. Meromorphic continuation of these full specified distributions propagates that same normalization at generic coupling. It does not license independently retuning a primitive three-current contact after continuation.

## 6. A precise projected triple-contour consequence

Take $\bar B$ on the inner segment $I=[L,R]$. Let $A_f$ have a smooth test function compactly supported strictly inside $I$. Let $B_g$ have $g=1$ on a neighborhood of $I$. With $\delta=c\epsilon$, $c>2$, the four ordered words can be represented by the following heights while $\bar B$ remains at zero:

\[
\begin{array}{c|cc}
\text{word}&\tau_B&\tau_A\\\hline
BA\bar B&2\delta&\delta\\
AB\bar B&\delta&2\delta\\
\bar BBA&-\delta&-2\delta\\
\bar BAB&-2\delta&-\delta.
\end{array}
\]

Use the same inner source circle and current-pole Cauchy extension in all four words. Add the same binary finite matching on each proper $B\bar B/Y_2$ face. Optional exclusions near new $A$ diagonals may be removed in $\operatorname{Re}\rho<0$ by section 5. The contour transport does not cross the order of $A$ and $B$ within either pair. Smooth movement of the $A$ test function contributes vanishing shrinking-contour terms; $f$ vanishes near both mixed endpoints.

The calculation in sections 2–5 shows which local terms this transport actually produces. The off-diagonal $A\bar B$ crossing has no contact. The diagonal crossing produces the computed endpoint $f_B=R_B-L_B$. There is no remaining simultaneous-face term in the starting domain. The two triple defects consequently are $s^{-1}f_BA_f$ and $A_ff_B$, respectively. These expressions follow from the evaluated contour faces; treating them as an abstract braided Leibniz identity would not supply the computation above.

The endpoint phase can be checked directly. From the anchored split, $[\Phi(L),\varphi(x)]=-2\pi i$ for $L<x<R$. The endpoint momentum $v=d+a_2$ obeys $vd=1-\rho$, so transporting $A=e^{-d\varphi}$ past the left endpoint gives $e^{2\pi i(1-\rho)}=s^2$, where $s=e^{-i\pi\rho}$. The right endpoint commutes. Therefore

\[
\boxed{\left.[B_gA_f-s^{-1}A_fB_g,\bar B_I]\right|_{\lambda_2}
=(s^{-1}-s)A_fL_B,\qquad
L_B={\nu\over2\pi i}t^2\lambda_2\mathscr U(L),}
\tag{T13}
\]

where the vertical restriction denotes the contribution linear in $\lambda_2$; terms involving other couplings and higher orders are outside this calculation. The same calculation in the unadjusted radial prescription has $t$ instead of $t^2$. The induced finite shift is explicitly

\[
\delta[\mathcal D,\bar B]
=(s^{-1}-s)A_f\,{\nu\over2\pi i}
t(t-1)\lambda_2\mathscr U(L).
\tag{T14}
\]

Equation (T13) holds first in the common convergence domain and then for its meromorphic continuation, away from exceptional poles. It is the finite-segment, projected first-order condition that the earlier algebraic compatibility report requested. No all-order associativity, equality with every four-fold prescription, or massive boundary-state statement follows from it.

## Integration

No book integration was performed. Independently verified: the exact first-Yukawa triple kernels, their cocycles, the additional-current pole residue, the action/tail contact and its finite normalization, all proper relative collision degrees relevant at this order, and the projected contour consequence. Modified interpretation: the missing first-order three-current calculation admits a common convergence domain once the proper mixed contact is resolved, even though other same-chirality tail products do not share their ultraviolet convergence domain. Rejected: the need to infer or tune a new primitive triple contact from the desired algebra in this prescription. Unresolved: arbitrary-test density-product extensions, higher perturbative orders, compatibility with four-fold Serre products, physical endpoint matrix elements, and a globally realized symmetry algebra.

## Checkpoint

- Completed: exact Wick and tail kernels, with explicit insertion-order cocycles.
- Completed: additional-current mixed contact, including the strip-edge derivative.
- Completed: proper and simultaneous collision bounds and finite endpoint matching propagation.
- Completed: projected first-order contour consequence; independent numerical output is stored alongside the scratch script.
