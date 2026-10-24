# Root-system sl(2|1) Toda: independent continuation audit and boundary-orbit constraints

Date: 2026-10-24. Model: Codex (GPT-6).

Commit reviewed: `745586e6a6f7bca5ed91a1e9cf89d6ef5df387f0`.

Task: continue the Bernard–LeClair investigation, independently auditing the latest finite-contact and mass-squared continuation calculations before examining boundary closure. `AGENTS.md` was read. Work is confined to this new report and `referee_scratch/sl21_boundary_orbits_20261024/`. No previous script was imported or rerun. No existing file was edited, no subagent was spawned, and no git mutation, commit, push, or book render was performed.

## Outcome

The local mass-squared normalization and the two marginal-channel removals survive a fresh audit. The previous quartic conclusion retains its stated scope: a formal finite-interval identity in the specified analytic product prescription. This pass does not extend its perturbative order or construct an associative physical charge product.

The additional result concerns the missing boundary step. The renormalized diagonal mixed endpoints have definite characters under a discrete period of the formal action. On a suitable vacuum orbit their boundary values would have the form

\[
\mathcal A_\alpha s^{-n},\qquad \mathcal B_\alpha s^n,
\qquad s=e^{-2\pi i/t}.
\tag{B1}
\]

The amplitudes and the vacuum orbit have not been constructed. A vacuum invariant under that period has **zero** endpoint expectation at generic coupling. Thus a nonzero Cartan-type boundary term requires more than local covariance.

The anchored canonical calculation also fixes a sign and a boundary choice. For the split used in the corrected mixed-regulator report, the candidate charges preserve the **right** boundary label and shift the left one. If the required boundary realization exists, fixing the right label to $r$ gives

\[
\Delta B_{BB}=\mathcal B_\alpha s^r(1-s^{-N}),\qquad
\Delta B_{AA}=\mathcal A_\alpha s^{-r}(1-s^N),
\tag{B2}
\]

where $N$ is the right-minus-left topological label. These are conditional boundary formulas. They are not yet right-hand sides of an established symmetry algebra.

| Item | Verdict | Independent evidence and limit |
|---|---|---|
| Local bosonization, dual-current residue and mass-squared Ward normalization | **Correct in the stated UV convention** | Fermion Wick signs, root products and distributional pole calculation |
| Finite mixed current/action contact | **Correct in the nested radial prescription** | Cauchy contour moments and strip-edge calculation reconstructed here |
| New mass-squared collision degrees and marginal-channel removal | **Correct in the specified analytic family** | Fresh subset enumeration, exact positivity, Wick derivatives and nested radial primitives |
| Identification of those products with an associative physical charge product | **Unclear** | Binary mixed matching and analytic four-fold products have not been joined into one construction |
| Endpoint neutrality under the local free direction | **Correct** | All relevant momenta and derivative vectors are orthogonal to it |
| Period characters and anchored endpoint weights | **Correct on the formal canonical cover** | Direct field shifts, canonical delta-function integration and sector-shift matrices |
| Nonzero physical endpoint amplitudes and a suitable vacuum orbit | **Unclear** | Neither is supplied by the perturbative contour calculation |
| Symmetry-algebra identification or an all-order Serre claim | **Unclear** | Boundary realization, product compatibility and higher orders remain open |

No existing sentence is identified here as requiring an error replacement. The new restrictions concern possible extensions of the preceding claims. No theorem upgrade or literature novelty claim is proposed.

## 1. Local conventions and independent normalization audit

Let $X=(\varphi^u,\varphi^v,h)$, with $X^a(z)X^b(w)\sim-\delta^{ab}\log(z-w)$, and write $\mathbb X=X+\bar X$. The two physical bosons are $\Phi=\varphi+\bar\varphi$; $H=h+\bar h$ is the auxiliary full fermion boson. Put

\[
a_0=(ib/\sqrt2,b/\sqrt2,i),\quad
a_2=(ib/\sqrt2,-b/\sqrt2,-i),\quad
t=b^2,\quad T=a_0+a_2,\quad d=-T/t,\quad k=d+T.
\]

The scalar products are $a_i^2=-1$, $a_0a_2=1-t$, $a_id=1$, and $d^2=-2/t$. The first two components of $a_i$ are denoted $p_i$. Direct fermion Wick pairing gives

\[
Y_i=i:\chi_i\bar\chi_i e^{p_i\Phi}:,
\qquad Y_0Y_2\sim +|r|^{2t-2}:e^{T\Phi}:.
\]

The crossing sign between the two fermion contractions cancels $i^2=-1$. Point splitting $\chi_0=e^{ih}$, $\chi_2=e^{-ih}$ gives $j=:\chi_0\chi_2:=i\partial h$ and $j(z)j(w)\sim(z-w)^{-2}$. These checks retain the fermion current in

\[
J_A=:e^{-dX}:,\qquad J_B=:C e^{dX}:,\qquad C=a_0\partial X.
\]

The screening-first residue vectors are $a_0+a_0^2a_0=0$ and
$a_0+(a_0a_2)a_2=-t(d+a_2)$. Thus the two candidate currents and their first Yukawa tails are consistent with the specified local cocycle. This does not identify a global compactification.

Use ordinary Euclidean area and

\[
S_{\rm int}=\pi^{-1}\int d^2w\,
(\lambda_0Y_0+\lambda_2Y_2+\lambda_1V_1+\kappa P),
\quad V_1=:e^{-T\Phi}:,\quad P=:e^{T\Phi}:,
\]
\[
\kappa=-\lambda_0\lambda_2\epsilon^{2t}/t,\qquad \kappa_R=0.
\tag{B3}
\]

The couplings $\lambda_0,\lambda_2$ have mass order one; $\lambda_1$ has mass order two but is independently Wick-normalized. Define $E=e^{kX+T\bar X}$, $L=CE$, $D=d\partial X$, $C_2=\partial C$, $D_2=\partial D$.

Expanding the two Yukawa point splittings gives

\[
N/E=\tfrac12\{C_2+C^2-[-tD_2-C_2+(-tD-C)^2]\}.
\]

The three current/potential poles have coefficients $tE$, $(C+tD)E$, and $[C_2+CD+t(D_2+D^2)/2]E$. The Ward response is **minus** the simple-minus-derivative-of-double-plus-half-second-derivative-of-triple combination, because the action insertion is $-\kappa/\pi$. Hence

\[
R_P/E=-tCD+\tfrac12t^2D_2-\tfrac12t^3D^2,
\qquad N-R_P/t=\partial L.
\tag{B4}
\]

This verifies the coefficient of $\lambda_0\lambda_2\epsilon^{2t}L$ in $H_B$, including its fermionic contribution. The first version of the new test omitted that minus action sign; its failed output is retained. Restoring the action sign reproduces the earlier formula without changing any input report.

## 2. Finite mixed contact and its prescription dependence

Here the inner interaction exclusion is $|w|>\epsilon$, the outer contours have heights $\pm\delta$, and $\delta=c\epsilon$ with $c>1$. Cauchy distributions are resolved before the interaction cutoff is removed. In the sharp limit their upper-minus-lower integrals are

\[
\int du\,(u-w)^{-1}=-2\pi i\Theta_\delta(w),\qquad
\int du\,(u-w)^{-2}=-2\pi i\partial_w\Theta_\delta(w),
\]

where $\Theta_\delta$ is the strip indicator. The second formula follows by differentiating the first as a distribution; it includes the strip edges.

For $q=t-1$, the independent first-jet expansion of the Wick kernel is

\[
e^{Du+A_Lw+A_R\bar w}
\left[\frac C{u-w}+\frac q{(u-w)^2}\right]
\left[\frac q{\bar w^2}-\frac{\bar C}{\bar w}\right].
\]

Let $K_L=(d+a_2)\partial X$, $K_R=(d+a_2)\bar\partial\bar X$. The identities $C+qD=-K_L$, $D+A_L=K_L$, $qA_R-\bar C=tK_R$ follow directly from the root vectors. Including the minus action sign gives the contact divided by $\nu\lambda_2\mathscr U$ as

\[
-2iq\{K_L(I-qJ_w)-tK_RJ_1\},\qquad \nu=e^{-i\pi/t},
\]

where $I=\int\Theta_\delta/\bar w^2$, $J_1=\int\partial_w\Theta_\delta/\bar w$, and $J_w=\int w\partial_w\Theta_\delta/\bar w^2$, with the same exclusion.

Set $\epsilon=1$. Horizontal slices give $\int dx/\bar w^2=2\sqrt{1-y^2}$ for $|y|<1$, and zero otherwise. Therefore $I=\pi$. At the two strip edges, $|y|=c>1$, the entire horizontal line is allowed. There $\int dx/\bar w=i\pi\operatorname{sgn}y$ and $\int dx\,w/\bar w^2=i\pi\operatorname{sgn}y$. Since $\partial_w\Theta=(i/2)[\delta(y-c)-\delta(y+c)]$, both $J_1$ and $J_w$ are $-\pi$. New quadrature checks these moments and the smooth Cauchy regularizations directly.

The bare contact is consequently

\[
-2\pi i\nu\lambda_2t(t-1)\partial_x\mathscr U,
\qquad \mathscr U=:e^{(d+a_2)\mathbb X}:
=i:\chi_2\bar\chi_2e^{(d+p_2)\Phi}:.
\tag{B5}
\]

Adding the tail contact leaves radial coefficient $t$. The $AA$ radial coefficient is $\gamma^2+\gamma=\gamma^2/t$, with $\gamma=t/(1-t)$. The canonical matching additions are therefore $t(t-1)\lambda_2\mathscr U$ and $-\gamma\lambda_1F_A$, respectively. This is the same substantive finite normalization as in the previous mixed report.

The endpoints used below are its canonically matched representatives:

\[
\begin{split}
B_{AA}^{\rm can}&=\gamma^2\lambda_1 F_A,\qquad F_A=:e^{-k\Phi}:,\\
B_{BB}^{\rm can}&=\lim_{\epsilon\to0}\left[
t^2\lambda_2\mathscr U
-\frac{t^2\lambda_0\lambda_2}{t-1}\epsilon^{2t-2}F
-t\lambda_0\lambda_2\epsilon^{2t}:C\bar C F:
\right]_{S_\epsilon},\qquad F=:e^{k\Phi}:.
\end{split}
\tag{B6}
\]

The independent pairing $-(d+a_2)a_0=t-2$ gives
$\mathscr U(0)Y_0(r)=|r|^{2t-4}(1+rC+\bar r\bar C+r\bar r C\bar C+\cdots)F$.
Its two angularly surviving radial primitives give exactly the subtractions in (B6). More descendants are required beyond the strip $\operatorname{Re}t>-1$; (B6) denotes their prescribed continuation there, not a complete bare-cutoff list at every coupling.

Thus the previously obtained local relation remains

\[
\mathcal C_{a\bar a}=\frac{\nu}{2\pi i}\Delta B_{aa}^{\rm can}+O(m^3).
\tag{B7}
\]

Its endpoints are finite and its spectators separated. Neither the finite matching nor the $O(m^3)$ notation supplies an infinite-volume boundary limit.

## 3. Audit of the mass-squared continuation argument

The independent exponents are $g$ for like-current factors, $-\rho$ for unlike-current factors, and $h-1$ for each chirality of the unlike screening pair, with $t=2/\rho$. They begin in $\operatorname{Re}g>6$, $\operatorname{Re}\rho<-1$, $\operatorname{Re}h>4$, and are restricted only after continuation to $g=\rho$, $h=t$.

Radial integration independently fixes the auxiliary matching terms

\[
\kappa_h=-\lambda_0\lambda_2\epsilon^{2h}/h,\qquad
H_{B,h}^{(2)}=(t/h)\lambda_0\lambda_2\epsilon^{2h}L.
\tag{B8}
\]

The factor $t/h$ is necessary before restriction. Higher Taylor descendants belong to the same continuation. There is no off-surface Ward claim.

For a cluster with $a$ labels $A$ and $b$ labels $B$, define $\mathcal E=g[\binom a2+\binom b2]-\rho ab$. Directly summing both chiral vertex powers, derivative contractions and relative integration dimensions gives

| New channel | Relative degree |
|---|---|
| Two unlike Yukawas, or matched fused graphs including cutoff powers | $\mathcal E+3a-2b+2h+1$ |
| Two equal Yukawas, or equal-source tails | $\mathcal E+3a-2b+5$ |
| Even-root action or its $A$ tail | $\mathcal E-a+2b+1$ |

The fresh script constructs these degrees separately for unfused sources, $K_2$/source, two tails, $P$, $L$, and the $A$ tail. It checks every current subset of both quartics. With no current the relative two-source measure is two, giving degrees $2h$ and $4$.

On $2<t<3$, $g=\rho=2/t$, $h=t$, every new double-source degree is strictly positive. Exact Bernstein coefficients of $t$ times each degree establish this on the whole interval. The even-root proper faces are also positive except for the fixed one-$A$ Cauchy pole and the $AA$ pair, whose coefficient is $-2\lambda_1r^{g-2}/(1-t)$. Its continued primitive supplies the stated pair subtraction. The proper $AAA$ and $AAB$ degrees are $3\rho-2$ and $1-\rho$, respectively.

There are two full marginal scalar channels. They require identities beyond power counting:

1. **Inherited $BBBA/Y_2$.** With majority positions $z_i$, minority position $a$, $\mathcal R=\prod_{i<j}(z_i-z_j)^g\prod_i(z_i-a)^{-\rho}$ and $G=(a-w)^2/\prod_i(w-z_i)^2$, the three derivative Wick factors give $(t-1)(t-2)\partial_a(\mathcal RG)/\rho$. The new symbolic computation uses arbitrary positions and also checks all three $K_2$ tail placements. The node-zero scalar is zero. Without a source both the scalar and the first common $C$ jet cancel, improving the pure-current full degree to at least two.
2. **New $AAAB/V_1$.** With majority positions $a_i$, minority position $z$, $\mathcal R_A=\prod_{i<j}(a_i-a_j)^g\prod_i(a_i-z)^{-\rho}$ and $G_A=(z-w)^2/\prod_i(w-a_i)^2$, the single derivative factor is $f_B=\sum_i(z-a_i)^{-1}-t/(z-w)$. Direct differentiation gives $\mathcal R_AG_Af_B=-\partial_z(\mathcal R_AG_A)/\rho$. Every $A$ tail obeys the same derivative identity. The minority-centered hole and its first variation vanish. The proper two-$A$ source integral retains its infinity contribution $-1$; the three-$A$ integral has none.

For the last statement, write a decaying rational kernel with real poles as $G=\sum_i[D_i/(w-a_i)^2+R_i/(w-a_i)]$. Integrating $\bar\partial(\bar wG)=G$ around the excluded disks gives $\pi^{-1}\operatorname{FP}\int G=-\sum_i a_iR_i$. Its $w^{-2}$ coefficient at infinity is $c_\infty=\sum_i(D_i+a_iR_i)$, so the integral equals $\sum_iD_i-c_\infty$. The two- and three-majority kernels have $c_\infty=1$ and $0$, respectively. This fixes the sign of the infinity contribution independently.

In the middle chambers the minority coordinate has internal majority endpoints. For initial $\operatorname{Re}\rho<0$ the primitives vanish there. Integration by parts moves their derivative onto a smooth spectator and improves the full degree from $3(g-\rho)$ to $1+3(g-\rho)$. This is then an identity of the continued proper-face extensions; it does not discard an unextended positive-$\rho$ endpoint value.

For a second unlike source, write $r$ for the size of the first cluster and $R$ for its distance to that source. The fused vertex has source-pair degree $2h-4$; differentiating it costs one $R$ power. The nested bound is

\[
\int R^{2h-4}\,dR\int_0^R r^{3(g-\rho)}\,dr.
\]

Its denominators are $1+3(g-\rho)$ and $3(g-\rho)+2h-2$. Neither is the removed divisor $g=\rho$. For an equal source the resulting full denominator is $3(g-\rho)+2$. The comparable-distance region is already the full cluster in the table.

Finally, a fresh enumeration of all labelled permutations gives the quartic shuffle defect $(0,-D_*,D_*,0)$, where

\[
D_*=(s-1)(s+1)(s-u)(s^2-u)(u+1)/s^4,
\quad u=e^{-i\pi g},\quad s=e^{-i\pi\rho}.
\]

A hypothetical chamber pole $1/[3(g-\rho)]$ would leave $i\pi(s-s^{-1})^2/3$, so the derivative arguments are essential. With those proper faces continued, the other new radial degrees are positive on the patch, and the physical restriction is regular. This supports the preceding mass-two contour identity in its common analytic family. The source-pair/cutoff normalization and the full Taylor prescription remain part of the result.

This audit does not identify the analytically continued four-fold products with iterates of the independently matched binary mixed product. Their compatibility and associativity remain additional calculations. There is also no all-order same-chirality assertion here.

## 4. The free direction does not force these endpoints to vanish

Let $R_\perp=(0,\sqrt2,ib)$ and $S=a_0-a_2$. Then

\[
R_\perp a_0=R_\perp a_2=R_\perp d=0,\quad
R_\perp^2=2-t,\quad T^2=-2t,\quad S^2=2(t-2).
\]

The three vectors $T,S,R_\perp$ are orthogonal and nondegenerate for $t\ne0,2$. Locally the Gaussian bosonized action splits into an interacting plane and a free direction. All current, tail and endpoint momenta in (B6), and the derivative vector $a_0$, lie in that plane. Their OPEs with the free current $R_\perp\partial X$ have no charge or derivative contraction.

Consequently a free-boson charge-selection rule alone cannot set $\langle\mathscr U\rangle$ to zero. This is a statement about local oscillator factors. Compactification, spin structures, zero modes and possible gluing of sectors have not been specified, so a global tensor-product Hilbert space does not follow.

## 5. Period characters and their conditional boundary consequence

On the complex field cover choose integers $n_0,n_2$ and a constant displacement satisfying $a_i\delta\mathbb X=2\pi i n_i$. One representative is

\[
\delta\Phi^u=\frac{\sqrt2\pi}{b}(n_0+n_2),\qquad
\delta\Phi^v=0,\qquad \delta H=\pi(n_0-n_2).
\tag{B9}
\]

All perturbations in (B3) are invariant. For the elementary choice $(n_0,n_2)=(1,0)$ the scalar part of each Yukawa changes sign; shifting one chiral fermion boson by $\pi$ sends both $\chi_0,\chi_2$ to their negatives and supplies the other sign. The bosonic displacement alone would not preserve the Yukawa action.

Since $d\delta\mathbb X=-2\pi i(n_0+n_2)/t$, direct exponentiation gives

\[
\mathscr U\longmapsto s^{n_0+n_2}\mathscr U,\qquad
F\longmapsto s^{n_0+n_2}F,\qquad
F_A\longmapsto s^{-(n_0+n_2)}F_A.
\tag{B10}
\]

The derivative counterterm has the same character as $F$. Thus (B6) is covariant under a regulator preserving this action period. Equality of characters does not imply an operator identity: the UV momenta of $\mathscr U$ and $F$ differ by $-a_0$, and their UV dimensions differ as well.

Let $\mathcal G$ denote the elementary displacement as an automorphism of the formal field algebra. To obtain (B1), additional assumptions would be needed: a quantum integration cycle/domain on which $\mathcal G$ acts; a family of boundary states defined with orientation $\omega_{\alpha,n}(O)=\omega_{\alpha,0}(\mathcal G^nO)$; and scalar weak cluster limits of the renormalized endpoints within each asymptotic sector. A one-point expectation alone is insufficient to prove such an operator boundary limit. Under those assumptions, (B10) gives (B1), with unknown amplitudes $\mathcal A_\alpha,\mathcal B_\alpha$.

There is an immediate selection rule. If $\omega\circ\mathcal G=\omega$, then

\[
(1-s)\omega(B_{BB})=0,\qquad
(1-s^{-1})\omega(B_{AA})=0.
\tag{B11}
\]

For generic $s\ne1$, both expectations vanish. Nonzero amplitudes would require appropriately selected non-invariant boundary states. The present work establishes neither their existence nor their nonzero condensates. Distinct orbit labels $\alpha_L,\alpha_R$ would also prevent reducing the two endpoints to a single function of a topological difference.

## 6. Which boundary can the canonical charge preserve?

For $z=x+i\tau=x-t_M$, use the corrected anchored split

\[
\varphi(x)=\tfrac12\left(\Phi(x)-\int_{-\infty}^x\dot\Phi\right),\qquad
\bar\varphi(x)=\tfrac12\left(\Phi(x)+\int_{-\infty}^x\dot\Phi\right),
\]

with $[\Phi^a(y),\dot\Phi^b(x)]=4\pi i\delta^{ab}\delta(y-x)$. Direct integration gives

\[
[\Phi(y),\varphi(x)]=-2\pi i\Theta(x-y),\qquad
[\Phi(y),\bar\varphi(x)]=+2\pi i\Theta(x-y).
\tag{B12}
\]

First keep endpoints $L<x<R$ strictly outside the charge support and to the right of the anchor regulator. Define $n_L=T\Phi(L)/(2\pi i)$, $n_R=T\Phi(R)/(2\pi i)$, $N=n_R-n_L$. For a chiral charge with disorder momentum $qd$, where $q=-1$ for $A$ and $q=+1$ for $B$, $Td=2$ yields

\[
[n_L,Q_q]=-2qQ_q,\qquad [n_R,Q_q]=0,\qquad [N,Q_q]=2qQ_q.
\tag{B13}
\]

Antichiral weights are opposite. All tails have left-minus-right momentum $qd$, so they preserve these weights. Derivative contacts have support inside the charge region and do not affect separated endpoints. The compatible cocycle in the corrected mixed report is a scalar rephasing of Weyl vertices; adding no independent zero-mode operator leaves (B12) unchanged. A different global Klein realization or constrained boundary phase space would need its own calculation.

The BCH check also gives

\[
F(L)Q_q=s^{-2q}Q_qF(L),\qquad F(R)Q_q=Q_qF(R).
\tag{B14}
\]

The same phases hold for $\mathscr U$ and its subtractions; $F_A$ has inverse phases. Indeed $kd=2-2/t$, whereas $(d+a_2)d=1-2/t$, and their difference contributes an integral $2\pi i$ phase. Finite sector-shift matrices independently check (B13)–(B14).

Thus it is the right boundary that can be fixed while retaining these charge actions. With the same orbit $\alpha$ on both ends and $n_R=r$, (B1) gives (B2). Fixing the left boundary as a scalar in this realization would conflict with (B13) unless a compensating operator or a larger boundary sector is supplied.

One may introduce the formal symbol $K=\nu^N$, so that $KQ_qK^{-1}=s^qQ_q$ and (B2) contains $1-K^{\mp2}$. This neither fixes the physical sign of a quantum-group parameter nor determines a presentation. In particular, $q=\pm1$ in (B13) is only a momentum label.

## 7. Two checks on the missing infrared input

### Classical real boundaries versus the complex period

With the stated roots, $\alpha_1\phi=\sqrt2u$ and $(\alpha_0+\alpha_2)\phi=-\sqrt2u$. At $\mu_0=\mu_2=\mu$, $\mu_1=\mu^2$, the body potential is therefore $V=2\mu^2\cosh(\sqrt2u)$. Its derivative has one real zero, $u=0$, and its Hessian there is $4\mu^2$. Real finite-energy classical boundaries have $u_L=u_R=0$.

The formal normalization $\Phi=\sqrt{4\pi}\phi/g_c$, $b=ig_c/\sqrt{4\pi}$ gives $T\Phi=-\sqrt2u$. Thus those real classical boundaries have $N=0$, while the elementary period (B9) shifts $u$ by $-i\sqrt2\pi$. The complex kink from $u=0$ to $u=i\sqrt2\pi$ has $N=-1$ in this orientation.

The period orbit is consequently not supplied by the ordinary real classical body contour. This is not a quantum nonexistence argument. It says which domain issue must be resolved. If a quantum representation were restricted to $N=0$ identically and preserved by these operators, (B13) would force $Q_q=0$ there. Nontrivial charges require compatible additional sectors or a different domain; classical complex saddles alone do not construct them.

### Exact frozen-boson fermion probe

Hold $\Phi=\Phi_0$ fixed and set $m_i=\lambda_i e^{p_i\Phi_0}$. In the ordered basis $(\chi_0,\chi_2,\bar\chi_0,\bar\chi_2)$, write the quadratic action as $(2\pi)^{-1}\int\Psi^T A\Psi$, where

\[
A=\begin{pmatrix}
0&\bar\partial&im_0&0\\
\bar\partial&0&0&im_2\\
-im_0&0&0&\partial\\
0&-im_2&\partial&0
\end{pmatrix}.
\]

The covariance is $\pi A^{-1}$, correctly giving $\chi_0(z)\chi_2(0)\sim1/z$ when the masses vanish. With $M^2=4m_0m_2>0$ on the test branch, direct inversion and Fourier transformation give

\[
\langle i\chi_2(x)\bar\chi_2(0)\rangle_{\Phi_0}
=-2m_0K_0(M|x|),
\]
\[
\langle i\chi_2\bar\chi_2\rangle_{\Phi_0,\Lambda}
=-m_0\log[(\Lambda^2+M^2)/M^2].
\tag{B15}
\]

Subtracting the ultraviolet log at reference scale $\mu$ leaves $2m_0\log(M/\mu)$. Coordinate point splitting with subtraction $2m_0\log(\mu|x|)$ differs by $2m_0(\gamma_E-\log2)$. A new Schwinger-parameter quadrature verifies the Bessel coefficient and the small-distance finite part.

This is an exact Gaussian calculation under its frozen-background assumption. It refutes the inference that a classical endpoint with $\psi=0$ must have zero quantum bilinear. It does **not** calculate $\langle\mathscr U\rangle$ in the Toda theory: boson fluctuations, vertex renormalization and the actual boundary state are missing. In particular, freezing the bosons removes the contractions responsible for $|r|^{2t-4}$ in the true endpoint/source OPE. The large dual momentum $d\propto b^{-1}$ also prevents interpreting this as an automatic weak-coupling approximation to the full endpoint.

The frozen expression has the expected character, since multiplication by its external scalar factor gives $m_0e^{(d+p_2)\Phi_0}=\lambda_0e^{k\Phi_0}$. It supplies no nonzero-amplitude proof for (B1).

## 8. Why another finite perturbative order does not settle boundary closure

With the free stress tensor and fixed Wick convention, the full UV dimensions are

\[
\Delta_F=2(t-1)^2/t,\qquad \Delta_{\mathscr U}=2/t-1,\qquad
[\lambda_1]=2-2t,\quad [\lambda_i]=1\ (i=0,2).
\]

Every term displayed in $B_{BB}$ has total dimension $2/t$, including its cutoff powers and derivatives. $B_{AA}$ has dimension $2/t-2$. These equal twice the respective candidate charge spins.

In a scale-covariant completion with a single mass scale $M$ and fixed dimensionless coupling ratios, the boundary amplitudes would carry powers $\mathcal B\sim M^{2/t}$ and $\mathcal A\sim M^{2/t-2}$, with reference-scale logarithms possible when the normalization requires them. Such behavior is generally nonanalytic in a formal small-mass expansion. Higher perturbative terms multiplied by infrared expectation values need not be suppressed after the infinite-volume limit. Accordingly, the remainder in (B7) must not be read as a uniformly small error in a massive boundary commutator.

The next substantive calculation is to specify and study a quantum boundary sector on which the action period and the nonlocal charge domains are meaningful, then compute the weak boundary limits of the full renormalized endpoints. The local free-direction decomposition narrows this to the interacting plane at oscillator level, but leaves global zero-mode constraints open. In parallel, an associative extension compatible with both the canonical mixed finite contact and the same-chirality analytic products is needed. These are independent requirements; neither is supplied by the other.

## Evidence and checkpoint

| New script | Output | PASS |
|---|---|---:|
| [audit_m2.py](../referee_scratch/sl21_boundary_orbits_20261024/audit_m2.py) | [output](../referee_scratch/sl21_boundary_orbits_20261024/audit_m2.out) | 123 |
| [inherited_faces.py](../referee_scratch/sl21_boundary_orbits_20261024/inherited_faces.py) | [output](../referee_scratch/sl21_boundary_orbits_20261024/inherited_faces.out) | 8 |
| [mixed_contact_recheck.py](../referee_scratch/sl21_boundary_orbits_20261024/mixed_contact_recheck.py) | [output](../referee_scratch/sl21_boundary_orbits_20261024/mixed_contact_recheck.out) | 14 |
| [boundary_checks.py](../referee_scratch/sl21_boundary_orbits_20261024/boundary_checks.py) | [output](../referee_scratch/sl21_boundary_orbits_20261024/boundary_checks.out) | 67 |
| [frozen_fermion_probe.py](../referee_scratch/sl21_boundary_orbits_20261024/frozen_fermion_probe.py) | [output](../referee_scratch/sl21_boundary_orbits_20261024/frozen_fermion_probe.out) | 14 |

Final total: **226 PASS, 0 FAIL**. The initial Ward test missed the minus action sign; an initial classical-potential test also required rewriting the symbolic hyperbolic sine in exponentials before simplification. Both initial outputs are retained. The final outputs distinguish symbolic identities from the Gaussian quadrature. The continuation and conditional boundary arguments are written above; passing scripts alone are not their proof.

The [checkpoint](../referee_scratch/sl21_boundary_orbits_20261024/checkpoint.md) records resumption constraints. Input hashes are in [inputs_before.sha256](../referee_scratch/sl21_boundary_orbits_20261024/inputs_before.sha256).

## Integration

No book integration was performed. **Independently verified:** local bosonization and dual residues, the mass-squared Ward polynomial, the finite mixed strip-edge contact, the new collision degrees and the two marginal derivative mechanisms. **Extended:** endpoint period characters, anchored boundary weights, conditional vacuum-orbit formulas, and an exact frozen-boson fermion probe. **Unresolved:** a physical quantum vacuum/domain, nonzero renormalized boundary amplitudes, scalar weak boundary limits, compatible associative charge products, higher-order relations and a symmetry-algebra identification. Existing reports and book files remain unchanged.
