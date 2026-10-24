# Root-system sl(2|1) Toda: renormalized currents and the local Ward identity

Date: 2026-10-23. Model: Codex (GPT-6).

Commit reviewed: 745586e6a6f7bca5ed91a1e9cf89d6ef5df387f0.

Task: construct the renormalized candidate current through order mass squared in the original fermionic variables, independently check the preceding cancellation, and include derivative contact terms in correlation functions. The user’s referee restrictions remain in force. No existing file was edited, no subagent was used, and no git mutation or book render was performed.

Independent evidence is in [fermionic_ward.py](../referee_scratch/sl21_current_renorm_20261023/fermionic_ward.py), [its output](../referee_scratch/sl21_current_renorm_20261023/fermionic_ward.out), [anomaly_count.py](../referee_scratch/sl21_current_renorm_20261023/anomaly_count.py), and [its output](../referee_scratch/sl21_current_renorm_20261023/anomaly_count.out). No preceding script was imported or used as evidence. The first script generates fermionic Wick contractions directly. There is no bosonization field in it.

## Outcome

The second-order cancellation survives the fermionic calculation. The calculation also fixes the current-insertion subtraction that the preceding report left open.

In the local cocycle convention specified below, the result is

\[
\boxed{
J_-^{\rm ren}=\operatorname{FP}J_0,\qquad
H_-^{\rm ren}=\lim_{\epsilon\to0}
\left[-t\lambda_2K_2+\lambda_0\lambda_2\epsilon^{2t}L\right]_{S_\epsilon},
\quad t=b^2.
}
\tag{1}
\]

The limit in (1) is an insertion in perturbative correlation functions, not a limit of its two summands separately. All fields inside it are ultraviolet Wick products. The definitions, entirely in terms of the two bosons and the original fermions, are

\[
\begin{aligned}
p_i&=-ib\alpha_i\quad(i=0,2),&
T&=p_0+p_2,& d&=-T/t,& k&=d+T,\\
j&=:\chi_0\chi_2:,&
J_0&=:(p_0\cdot\partial\varphi+j)e^{d\cdot\varphi}:,\\
K_2&=i:\chi_2\bar\chi_2
e^{(d+p_2)\cdot\varphi+p_2\cdot\bar\varphi}:,&
E&=:e^{k\cdot\varphi+T\cdot\bar\varphi}:,&
L&=:(p_0\cdot\partial\varphi+j)E:.
\end{aligned}
\tag{2}
\]

Here $\Phi=\varphi+\bar\varphi$ is the canonical two-component boson and $\chi_i=\alpha_i\cdot\psi$. A bar labels antiholomorphic fields, not Hermitian conjugation. Products involving the disorder momentum $d$ include the explicitly specified local cocycle of section 2.

The action counterterm is

\[
\kappa(\epsilon)P
=-\frac{\lambda_0\lambda_2}{t}\epsilon^{2t}
:e^{T\cdot\Phi}:.
\tag{3}
\]

With no finite independent $P$ coupling, (1) satisfies

\[
\boxed{
\bar\partial_z\langle J_-^{\rm ren}(z)\mathcal O\rangle
-\partial_z\langle H_-^{\rm ren}(z)\mathcal O\rangle
=O(m^3),
}
\tag{4}
\]

for separated spectator insertions $\mathcal O$ and $z$ away from them. The finite-cutoff remainder from the collision channel is $O(\epsilon^{2\operatorname{Re}t+2})$. The calculation retains the derivatives of delta functions and of the exclusion boundary; it does not pass to a quotient by total derivatives.

This statement uses circular exclusion, canonical Cauchy distributions for the current poles, and minimal analytic continuation with zero finite $P$ and $L$ additions. The explicit cutoff result applies in $\operatorname{Re}t>-1$, $t\ne0$. The two-current discussion also excludes $t=1$. The local charge-lattice convention and the generic-coupling statements below are made at generic, in particular irrational real, $t$. No statement about a global quantum Hilbert space follows from these choices.

**A correction in interpretation is necessary.** Equation (27) of the preceding mass-squared report is a primitive in a chosen contour-descent representative. It is **not** the current-insertion counterterm in (1). A concrete fermionic spectator correlator distinguishes them: the old primitive has zero matrix element in that test, while the required subtraction has a nonzero matrix element and removes a divergence.

The third-order anomaly count gives an additional result. In mass-analytic Gaussian conformal perturbation theory, there is no allowed ordinary local third-order anomaly at generic irrational $t$. The same integer argument excludes ordinary local higher-order anomalies under the assumptions in section 7. This is a conditional selection rule, not a construction of global charges or their algebra.

| Item | Verdict | Evidence |
|---|---|---|
| Previous action counterterm and cancellation of the non-derivative class | **Correct** in the specified prescription | Fresh fermionic endpoint calculation and source normalization |
| Current-insertion subtraction in (1) | **Correct through order $m^2$** in this prescription | Independent operator-mixing calculation and full distributional Ward identity |
| Using the preceding report’s equation (27) as that subtraction | **Wrong** | Nonzero uncanceled fermionic matrix element, section 6 |
| First-order residue signs without a common cocycle convention | **Unclear as a pointwise prescription** | The signs in the two OPE orders must be linked; section 2 fixes them |
| Generic third-order, and conditional higher-order, local-anomaly exclusion | **Correct under the listed conformal-perturbation assumptions** | Explicit descendant-level calculation and integer argument |
| Existence of global charges, their algebra, and their action on solitons | **Unclear; not established here** | Global contour and operator-domain questions remain |

## 1. Input normalization, checked again

The source input is Evans–Madsen, [hep-th/9703065v2](https://arxiv.org/abs/hep-th/9703065v2), equation (3.4), with overall multiplier $g^{-2}$. Its odd-root double sum has coefficient $\tfrac12\mu_i\mu_j(\alpha_i\cdot\alpha_j)$, and its Yukawa term has coefficient $i\mu_j$. The source on printed page 8 says exactly: “If desired one can introduce a dimensionless coupling β by taking L → L/β²”. This source beta is called $g$ here.

Use the complex bilinear roots

\[
\alpha_0=(-1/\sqrt2,i/\sqrt2),\qquad
\alpha_1=(\sqrt2,0),\qquad
\alpha_2=(-1/\sqrt2,-i/\sqrt2).
\]

Then $\alpha_0^2=\alpha_2^2=0$ and $\alpha_0\cdot\alpha_2=1$. The two cross terms in the source double sum give exactly $\mu_0\mu_2e^{(\alpha_0+\alpha_2)\cdot\phi}$.

With $\partial_\pm=\partial_{t_M}\pm\partial_x$, set $t_M=-i\tau$, $z=x+i\tau$, and $d^2z=dx\,d\tau$. The rescalings

\[
\phi=\frac{g}{\sqrt{4\pi}}\Phi,\qquad
\psi_-=e^{i\pi/4}\frac{g}{\sqrt{2\pi}}\psi,\qquad
\psi_+=e^{-i\pi/4}\frac{g}{\sqrt{2\pi}}\bar\psi
\]

give

\[
S_0=\frac1{2\pi}\int d^2z\,
\left(\partial\Phi\cdot\bar\partial\Phi+
\psi\cdot\bar\partial\psi+
\bar\psi\cdot\partial\bar\psi\right).
\]

Directly reordering the Euclidean Yukawa term gives

\[
b=\frac{ig}{\sqrt{4\pi}},\qquad
t=-\frac{g^2}{4\pi},\qquad
\lambda_i=\frac{\mu_i}{2},\qquad
Y_i=i:\chi_i\bar\chi_i e^{p_i\cdot\Phi}:.
\tag{5}
\]

The perturbation is $\pi^{-1}\int(\lambda_0Y_0+\lambda_2Y_2+\lambda_1V_1+\kappa P)$, with $V_1=:e^{-T\cdot\Phi}:$. Count the two Yukawa couplings as order $m$, and $\lambda_1,\kappa$ as order $m^2$.

The propagators are

\[
\varphi^a(z)\varphi^b(w)\sim-\delta^{ab}\log(z-w),\qquad
\chi_0(z)\chi_2(w)\sim\frac1{z-w}.
\]

The self-contractions of each $\chi_i$ vanish. Left and right fermions anticommute. The mixed full Yukawa OPE has leading coefficient $+1$: the fermionic Wick crossing and $i^2$ each supply one minus sign. Its leading term is

\[
Y_0(w+r)Y_2(w)\sim |r|^{2t-2}P(w).
\tag{6}
\]

The individual bosonic Yukawa momenta are null, whereas $T^2=-2t$. With self-contraction subtraction $-\log\epsilon^2$, the source composite exponential is $\epsilon^{2t}P$. Its ordinary-area coefficient is consequently

\[
\frac{\pi\mu_0\mu_2}{g^2}\epsilon^{2t}
=-\frac{\lambda_0\lambda_2}{t}\epsilon^{2t},
\]

which verifies (3) without assuming the preceding report’s conclusion.

Bernard–LeClair’s supplied PDF was checked at printed page 105, equations (2.9)–(2.12). Its qualification “local with respect to the perturbing field” is relevant to the cocycles below. The distributional calculation uses their conformal-perturbation framework, but does not infer its second-order result from their first-order equations. [Published source](https://doi.org/10.1007/BF02099173).

## 2. Fix the OPE signs before computing contact terms

The current/Yukawa bosonic contraction has an odd integer power:

\[
d\cdot p_0=d\cdot p_2=1.
\]

Thus switching from $(z-w)^{-1}$ to $(w-z)^{-1}$ changes a sign. A list of separately ordered residues is not yet a consistent pointwise current prescription. The cocycle factors must be specified together with the orientation.

Use the additive charge cover generated by $(d,p_0,p_2)$. Write a charge as $(n_d,n_0,n_2)$ and choose the associative bimultiplicative cocycle

\[
\varepsilon(n,m)=(-1)^{n_d(m_0+m_2)}.
\tag{7}
\]

Elementary fermions retain their ordinary Grassmann signs. The cocycle in (7) concerns the bosonic exponential factors; it does not change the fermion parity of the current. At generic irrational $t$ the relation $td+T=0$ is not an integer relation on this cover.

In particular, $\varepsilon(d,p_i)=-1$, $\varepsilon(p_i,d)=+1$, and $\varepsilon(d,T)=+1$. The even-root charge $-T$ also has sign $+1$ with the current. The same minus sign occurs in the mixed product of $K_i$ with the other Yukawa. The two signs cancel in a three-insertion product containing the current and both Yukawas.

For transparency, let $\sigma$ denote the common current-first sign, with **$\sigma=-1$ for (7)**. The elementary fermionic contractions, before multiplication by $\sigma$, give

\[
\begin{array}{c|cc}
& (z-w)^{-2}\text{ coefficient}&(z-w)^{-1}\text{ coefficient}\\ \hline
J_0(z)Y_0(w)&K_0&+\partial K_0\\
J_0(z)Y_2(w)&(t-1)K_2&-\partial K_2
\end{array}
\tag{8}
\]

Both columns acquire $\sigma$ in the specified OPE. Here $K_0$ is defined as $K_2$ in (2), with $0$ replacing $2$.

These entries can be obtained without a fermion bosonization:

\[
j(z)\chi_0(w)=\frac{\chi_0(w)}{z-w}+\partial\chi_0(w)+\cdots,\qquad
j(z)\chi_2(w)=-\frac{\chi_2(w)}{z-w}-\partial\chi_2(w)+\cdots.
\tag{9}
\]

The bosonic derivative contributes $-p_0\cdot p_i/(z-w)$, with $p_0^2=0$ and $p_0\cdot p_2=-t$. Expanding the exponential then gives (8). The new script generates the contractions in (9), including their regular terms, and checks both poles in (8).

Use the canonical distribution identities

\[
\bar\partial_z\frac1{(z-w)^n}
=\frac{\pi(-1)^{n-1}}{(n-1)!}
\partial_z^{\,n-1}\delta^{(2)}(z-w).
\tag{10}
\]

The first-order term in the Euclidean weight is $-\lambda_i\pi^{-1}\int Y_i$. Equation (8), including its double pole, therefore gives

\[
\bar\partial J_0\big|_{m}
=\sigma t\lambda_2\partial K_2
=-t\lambda_2\partial K_2.
\tag{11}
\]

This fixes the sign of $H^{(1)}$ in (1). Omitting the cocycle while keeping the same symbols for the fused fields instead gives $\sigma=+1$ and reverses this first-order formula. It must also reverse the associated $K_iY_j$ product sign. The second-order counterterm below is invariant under that simultaneous change.

The old raw charge-contour signs and the old fixed-current signs should therefore not be combined in a pointwise distribution calculation without this additional convention. The full mixed collision and its non-derivative class are unchanged. Equation (7) is a local associative convention; it is not yet a construction of global disorder lines or a charge action on soliton sectors.

## 3. Fermionic expansion of the collision endpoints

Keep $z$ as the current insertion and put the other interaction at $w=z+r$. Define

\[
F_0(z,w)=\sigma K_0(z)Y_2(w),\qquad
F_2(z,w)=\sigma K_2(z)Y_0(w).
\]

These are the combinations that occur in the poles of the full three-insertion product. Their leading factors are independent of $\sigma$:

\[
F_0=-r^{t-2}\bar r^{t-1}G_0(r,\bar r)E,\qquad
F_2=-r^{t-2}\bar r^{t-1}G_2(r,\bar r)E.
\tag{12}
\]

All normal products in this section are fused at $z$. Only the terms with antiholomorphic Taylor order zero are needed for the leading contact. Define

\[
U=:\partial\chi_0\,\chi_2:,\qquad
V=:\chi_0\,\partial\chi_2:,\qquad
\partial j=U+V.
\]

The normalized chiral fermion factors in the two endpoints are obtained directly from

\[
\begin{aligned}
\chi_0(z)\chi_2(z+r)
&=-r^{-1}+j+rV+\cdots,\\
\chi_2(z)\chi_0(z+r)
&=-r^{-1}-j-rU+\cdots.
\end{aligned}
\]

Consequently

\[
\begin{aligned}
G_0(r,0)
&=1+r(p_2\cdot\partial\varphi-j)+r^2Q_0+\cdots,\\
G_2(r,0)
&=1+r(p_0\cdot\partial\varphi+j)+r^2Q_2+\cdots,\\
Q_0
&=\tfrac12\bigl[p_2\cdot\partial^2\varphi+
(p_2\cdot\partial\varphi)^2\bigr]
-j(p_2\cdot\partial\varphi)-V,\\
Q_2
&=\tfrac12\bigl[p_0\cdot\partial^2\varphi+
(p_0\cdot\partial\varphi)^2\bigr]
+j(p_0\cdot\partial\varphi)+U.
\end{aligned}
\tag{13}
\]

No identity involving the square of a bosonized fermion current is used. In particular, the regular fermionic products $U$ and $V$ are retained separately before taking their difference.

The resulting local boundary polynomial is

\[
\boxed{
N=(Q_2-Q_0)E
=:\left[
\frac12(p_0-p_2)\cdot\partial^2\varphi+\partial j
+\left(\frac12(p_0-p_2)\cdot\partial\varphi+j\right)
(T\cdot\partial\varphi)
\right]E:.
}
\tag{14}
\]

The composite-potential insertion gives

\[
R_P=:\left[
(p_0\cdot\partial\varphi+j)(T\cdot\partial\varphi)
-\frac t2\left(T\cdot\partial^2\varphi+
(T\cdot\partial\varphi)^2\right)
\right]E:.
\tag{15}
\]

Its triple, double and simple current poles all contribute through (10). In particular, it is not sufficient to retain only the current-first residue. The script independently expands these three poles and verifies that their distributional combination is $R_P$ with the sign in (15).

The exact identity that replaces the earlier quotient-by-derivatives calculation is

\[
\boxed{N-\frac1tR_P=\partial L.}
\tag{16}
\]

It follows by using $T=-td$, $k=d+T$ and $\partial j=U+V$. The fermionic term $\partial j+j(k\cdot\partial\varphi)$ on the right is essential.

## 4. Full local Ward identity with the exclusion boundary differentiated

The mixed term in the Euclidean weight is

\[
\frac{\lambda_0\lambda_2}{\pi^2}
\int d^2w_0d^2w_2\,
\Theta_\epsilon(w_0,w_2)Y_0(w_0)Y_2(w_2),
\qquad
\Theta_\epsilon=\Theta(|w_0-w_2|-\epsilon).
\tag{17}
\]

Keep $\epsilon$ fixed while applying (10). A small current-pole regulator is removed first. This makes the separate current/Yukawa contacts well-defined before the two Yukawas are allowed to collide.

For the current contact with node 0 the pole coefficients are $\partial F_0,F_0$; for node 2 they are $-\partial F_2,(t-1)F_2$. Set

\[
I_i(z)=\int d^2w\,\Theta(|w-z|-\epsilon)F_i(z,w),\qquad
B_i(z)=\int d^2w\,(\partial_z\Theta)F_i(z,w).
\]

The derivatives of the delta functions in (10) act on the exclusion factor as well as the OPE coefficient. The node-0 contribution is $-B_0$. The node-2 contribution is $-t\partial I_2+B_2$. Therefore

\[
\boxed{
\bar\partial\langle J_0(z)\mathcal O\rangle_{02}
=-\frac{t\lambda_0\lambda_2}{\pi}\partial I_2
+\frac{\lambda_0\lambda_2}{\pi}(B_2-B_0).
}
\tag{18}
\]

The first term is precisely the interaction correction to $\partial H^{(1)}$: the extra Yukawa insertion contributes $-\sigma t\lambda_0\lambda_2\pi^{-1}\int K_2Y_0=-t\lambda_0\lambda_2\pi^{-1}I_2$.

For $r=w-z$,

\[
\partial_z\Theta(|r|-\epsilon)
=-\tfrac12e^{-i\theta}\delta(|r|-\epsilon).
\]

The minus sign in (12) cancels this minus sign. Angular integration gives

\[
B_2-B_0
=\pi\epsilon^{2t}N+
O(\epsilon^{2\operatorname{Re}t+2}).
\tag{19}
\]

The angular weight selects $r^2\bar r^0$ in (13). The next possible term is $r^3\bar r^1$, which explains both the error estimate and the lower bound on $\operatorname{Re}t$.

Adding the potential insertion gives the full identity before current subtraction:

\[
\boxed{
\bar\partial J_0-\partial H^{(1)}
=\kappa(\epsilon)R_P+
\lambda_0\lambda_2\epsilon^{2t}N+
O(m^3)+O(\epsilon^{2\operatorname{Re}t+2}).
}
\tag{20}
\]

All terms in (20) are local fields at $z$ inserted into the spectator correlation function. With (3), equation (16) turns its right-hand side into
$\lambda_0\lambda_2\epsilon^{2t}\partial L$, giving (1) and (4).

The same-node $Y_iY_i$ OPE is regular and does not supply a missing mixed contact. The even-root perturbation is already order $m^2$ and has a regular OPE with $J_-$. The partner $J_+=:e^{-d\cdot\varphi}:$ has regular OPEs with both Yukawas and $P$; its usual even-root correction survives unchanged at this order.

For a test function supported away from the spectator positions, (20) is an identity of distributions. At a spectator insertion there are the usual external Ward contacts, obtained by applying (10) to its current OPE and the antiholomorphic counterpart to the other current component. They are not set to zero or counted as anomalies. The calculation establishes the interaction-collision terms in that identity; it does not assign a new charge representation to arbitrary spectator fields.

## 5. Independently renormalize the current insertion

The subtraction in (1) can be found without using the Ward identity. From (12)–(13),

\[
K_2(z)Y_0(z+r)
=-\sigma r^{t-2}\bar r^{t-1}
\left[E+rL+\cdots\right].
\tag{21}
\]

Angular integration of the constant term vanishes. In a disk of fixed small outer radius $R$, the $rL$ term gives

\[
\int_{\epsilon<|r|<R}d^2r\,K_2(z)Y_0(z+r)
=-\frac{\sigma\pi}{t}
(R^{2t}-\epsilon^{2t})L+\cdots.
\tag{22}
\]

Combining the coefficient $\sigma t\lambda_2$ in $H^{(1)}$ with the extra insertion $-\lambda_0/\pi$ gives

\[
\langle H^{(1)}\mathcal O\rangle_{02}
=\lambda_0\lambda_2(R^{2t}-\epsilon^{2t})
\langle L\mathcal O\rangle_0+\cdots.
\tag{23}
\]

For $-1<\operatorname{Re}t<0$, its divergent part is
$-\lambda_0\lambda_2\epsilon^{2t}L$. It is canceled by **exactly** the term in (1), including its coefficient and fermion bilinear. The remaining omitted descendants are locally integrable in this strip. For $\operatorname{Re}t>0$, both this divergence and the subtraction vanish as the cutoff is removed; analytic continuation fixes the subtraction on the negative side.

There is no additional local subtraction of $J_0$ that survives at these mass orders in this prescription. At one Yukawa insertion, angular projection of the chiral pole first contributes at order $\epsilon^2$. In the mixed order-$m^2$ channel the fused vertex $E$ has spin $1/t-2$, while $J_0$ has spin $1+1/t$. A local insertion correction needs at least three holomorphic descendant levels and scales as $\epsilon^{2t+2}$. The potential counterterm has the same bound. The regular same-node channels have still softer local behavior. Canonical Cauchy distributions specify the purely chiral poles rather than treating them as absolutely integrable functions.

The finite prescription is now explicit:

1. Fix the ultraviolet Wick normalization and overall coefficient of $J_0$ as in (2).
2. Use ordinary area, the cocycle (7), and the Cauchy distributions (10).
3. Exclude interaction pairs at radius $\epsilon$. In the $K_2$ insertion use the same radius; obtain it from the collapsed current contact as in (18).
4. Include (3), with no finite independent $P$ coupling.
5. Add $+\lambda_0\lambda_2\epsilon^{2t}L$ to the current component and add no finite independent multiple of $L$. This is the continuation from the convergent positive-$t$ calculation.

If an independent finite $\kappa_RP$ is retained, the right-hand side of (20) instead has the additional term $\kappa_RR_P$. Removing its derivative part leaves

\[
\kappa_R\mathcal A,\qquad
\mathcal A=\tfrac12:
(T\cdot\partial\varphi)
\bigl((p_0-p_2)\cdot\partial\varphi+2j\bigr)E:,
\]

in agreement with the preceding report. Its non-derivative obstruction is not removed by the new current subtraction.

## 6. A concrete spectator check, and why the older primitive is insufficient

Take chiral fermion spectators $\chi_0(A)\chi_2(B)$, with distinct nonzero $A,B$ outside the collision disk. Balance the bosonic charge at infinity and choose zero local bosonic derivative data. Normalize by the nonzero spectator two-point function. Then

\[
\frac{\langle j(0)\chi_0(A)\chi_2(B)\rangle_0}
{\langle\chi_0(A)\chi_2(B)\rangle_0}
=\frac{A-B}{AB}=j_{AB}\ne0.
\tag{24}
\]

Direct fermionic Wick contraction gives the normalized endpoint factor

\[
G_2(r,0)=1+\frac{r(A-B)}{A(B-r)}.
\tag{25}
\]

For $R<\min(|A|,|B|)$ its angular projection selects exactly $rj_{AB}$. Hence the local contribution of $H^{(1)}$ is

\[
\lambda_0\lambda_2(R^{2t}-\epsilon^{2t})j_{AB},
\]

and the new subtraction makes it
$\lambda_0\lambda_2R^{2t}j_{AB}$, independent of $\epsilon$.

By contrast, the previous report’s equation (27), after substituting its action counterterm, is

\[
\lambda_0\lambda_2\epsilon^{2t}
\frac{t+1}{2t}:(T\cdot\partial\varphi)E:.
\tag{26}
\]

It has zero matrix element in this test. It cannot cancel the nonzero divergence in (23). The difference is not a dispensable numerical normalization: (1) contains the fermion current $j$, while (26) does not.

The previous report explicitly warned that its expression was only a contact primitive. The finding here is that it must not be promoted to the finite pointwise current. Its non-derivative Ward-class conclusion remains correct.

The numerical check integrates the actual rational function (25) over angle and radius, using $A=1.3+0.4i$, $B=-1.1+0.3i$, $R=0.3$, and $t=\pm\sqrt2/3$. After the subtraction, the result agrees with the cutoff-independent value to numerical precision. A separate collision-circle check includes nonconstant bosonic spectator Taylor data, as well as exact fermionic spectator contractions. Its relative error falls by a factor of approximately four each time the radius is halved, as predicted by (19).

At $t=-\sqrt2/3$, reducing $\epsilon$ from $0.08$ to $0.01$ increases the magnitude of the unsubtracted contribution from $11.9379$ to $114.208$. The largest absolute error after subtraction is $1.47\times10^{-13}$. In the nonconstant-background circle test, the relative error at $\epsilon=0.01$ is $2.63\times10^{-5}$ for positive $t$ and $1.56\times10^{-5}$ for negative $t$; both exhibit the predicted quadratic decrease.

## 7. What the next mass order permits

The explicit current result allows a useful preliminary test before computing more diagrams. Assume:

- a mass-analytic expansion in $\lambda_0,\lambda_2,\lambda_1$;
- the Gaussian ultraviolet Fock fields, with no background charge and with momentum and fermion charge conserved in the operator products;
- local ultraviolet Ward anomalies built from those fields and their derivatives;
- the homogeneous minimal prescription above, with no independent finite $P$ coupling or extra finite scale-dependent deformation.

These are assumptions about formal conformal perturbation theory. They do not establish convergence, infrared completeness, or global conservation.

Consider $n_0,n_2$ Yukawa insertions and $n_e$ even-root insertions. Put

\[
a=n_0-n_e,\qquad c=n_2-n_e,\qquad
S=a+c,\qquad F=n_0-n_2,\qquad
M=n_0+n_2+2n_e.
\]

The added bosonic momentum is $q=ap_0+cp_2$ and the chiral and antichiral fermion charges are both $F$. The minimum fermion weight at integer charge $F$ is $F^2/2$: for positive charge it is the weight of
$\chi_0\partial\chi_0\cdots\partial^{F-1}\chi_0$, with the analogous formula for negative charge. This derives the bound without bosonization.

The fused field’s lowest weights are

\[
h_{\rm base}=1/t-S+tac+F^2/2,\qquad
\bar h_{\rm base}=tac+F^2/2.
\]

The product of couplings has weights
$M/2-n_et$ in each sector. To contribute to $\bar\partial J_-$, a local descendant must therefore have nonnegative integer levels

\[
\begin{aligned}
\ell&=1-M/2+S-F^2/2+t(n_e-ac),\\
\bar\ell&=1-M/2-F^2/2+t(n_e-ac).
\end{aligned}
\tag{27}
\]

At order $m^3$ this gives:

| Insertions | $(\ell,\bar\ell)$ |
|---|---|
| Three Yukawas of the same node | $(-2,-5)$ |
| Two of one odd node and one of the other | $(2-2t,-1-2t)$ |
| One even-root vertex and one Yukawa | $(t-2,t-1)$ |

No row gives two nonnegative integers at generic irrational $t$. The first mixed triple-Yukawa resonance is $t=-1/2$, with levels $(3,0)$. In the strip $-1<t<0$ it is the only third-order resonance. The first even-root/Yukawa resonance is $t=2$, with levels $(0,1)$. At these resonances, a diagrammatic calculation is still required; the table does not say that the permitted anomaly is nonzero.

The generic argument extends further. For irrational $t$, integrality in (27) requires $n_e=ac$. Substitution gives

\[
\ell=1-\frac{S(S-1)}2,\qquad
\bar\ell=1-\frac{S(S+1)}2.
\tag{28}
\]

Both are nonnegative only for $S=-1,0,1$. Since $ac=n_e\ge0$, the integers $a,c$ have the same sign or one is zero. For $S=1$ the only pairs are $(1,0),(0,1)$, which are the already calculated single-Yukawa terms. For $S=0$ only the vacuum term remains. For $S=-1$ the pairs $(-1,0),(0,-1)$ would give a negative original insertion count and are inadmissible.

Thus **there is no ordinary local higher-order anomaly compatible with these assumptions at generic irrational $t$**. This is a complete integer argument for the stated selection rule; finite enumeration in the script is only a check. It supports the formal perturbative conservation program. It does not prove the existence of globally defined charges, settle exceptional rational couplings, or identify their quantum algebra.

## 8. Verification and proposed corrections

Run the new checks from the repository root:

~~~text
python3 -u referee_scratch/sl21_current_renorm_20261023/fermionic_ward.py
python3 -u referee_scratch/sl21_current_renorm_20261023/anomaly_count.py
~~~

The fermionic calculation checks both current/Yukawa poles, the mixed full-vertex sign, both endpoint expansions, the derivative-of-cutoff terms, all three composite poles, the exact local identity (16), the current-operator mixing, and an explicit associative cocycle. Its numerical tests use irrational $t$ on both sides of zero. The second script checks (27)–(28), the third-order list, and the integer solutions used in the generic selection rule.

All **50 checks pass**: 42 in the fermionic calculation and 8 in the anomaly count. The numerical spectator harness initially exposed a sign error in its expected fermion one-point value; this was corrected directly from the four-fermion Wick contraction before recording the final output. A NumPy Boolean was also normalized in the check reporter. Neither change altered the symbolic Ward identity or counterterm. Hashes recorded before the work confirm that AGENTS.md and all four preceding reports are unchanged.

For the preceding report, briefs/referee_sl21_ward_m2_report.md, lines 410–419, replace the interpretation of equation (27) by:

> The displayed expression is a primitive in a contour-descent representative. It is not the current-insertion subtraction in the matched circular prescription. In that prescription the finite second-order component is defined in correlation functions by
> $H_-^{\rm ren}=\lim_{\epsilon\to0}[-t\lambda_2K_2+\lambda_0\lambda_2\epsilon^{2t}L]$,
> where $L=:(p_0\cdot\partial\varphi+:\chi_0\chi_2:)e^{(d+T)\cdot\varphi+T\cdot\bar\varphi}:$.
> The full local collision polynomial obeys $N-R_P/t=\partial L$.

For that report’s charge-contour formula (13), lines 233–239, add:

> These raw ordered residues do not by themselves fix the cocycle convention. With the local convention $\varepsilon(d,p_i)=-1$, $\varepsilon(p_i,d)=+1$, the current-first simple-pole residues acquire a minus sign, as do the corresponding mixed $K_iY_j$ products. The ordered fixed-current first-order component and the mixed second-order cancellation are then mutually consistent.

These are recommendations only. The preceding report was left unchanged. No [Theorem] upgrade is proposed.

## Integration and next unresolved calculation

No book integration was performed. **Verified independently:** the source normalization, the action subtraction, the first-order current residues, and the cancellation of the second-order non-derivative class. **Completed:** a matched current-insertion subtraction and the full local distributional Ward calculation through order $m^2$, with a nontrivial fermionic spectator check. **Modified:** the previous derivative primitive must not be used as the finite current counterterm; the current-first and perturbation-first signs now have a common local cocycle prescription. **Extended:** a conditional generic-coupling exclusion of ordinary local higher-order anomalies. **Not established:** global contour definitions, resonant-coupling Ward identities, the charge algebra, or action on soliton states.

For generic coupling, the next substantive step is to construct antichiral partners in the same convention and evaluate regulated mixed charge products and endpoint terms. That calculation should determine the coproduct and algebraic relations, without choosing an affine algebra in advance. At $t=-1/2$, the third-order triple-Yukawa channel must instead be checked explicitly before using the generic selection rule.
