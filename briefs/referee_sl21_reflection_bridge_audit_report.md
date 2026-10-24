# Root-system sl(2|1) Toda: reflection normalization and the Ward-endpoint bridge

Date: 2026-10-24. Model: Codex (GPT-6).

Commit reviewed: `745586e6a6f7bca5ed91a1e9cf89d6ef5df387f0`.

Task: independently audit V1 and V13-V22 of `briefs/referee_sl21_vacuum_limits_report.md`, particularly their normalization bridge, coupling branches and identity interpretation. Read `AGENTS.md`. No subagents were spawned, no pre-existing file was edited, and no git mutation, commit, push or book render was performed. This report and `referee_scratch/sl21_reflection_bridge_audit_20261024/` are new. The checkpoint is in that scratch directory, within the task's explicitly permitted write locations.

## Outcome

The two reflection formulas have the stated signs and gamma factors on their specified generic branch. The zero of the second reflection amplitude must multiply its reflected field's identity-resonance pole. Their finite product is nonzero at the generic irrational couplings under consideration. The local cutoff matching is correct through its displayed mass order.

A concrete normalization bridge is now available. Changing the finite Gaussian subtraction changes the OPE coefficients, couplings, primary values and Ward endpoint coefficients together. For a change by $c$ in the two physical bosons, with the fermion normalization fixed, both primary formulas acquire the factor $e^{c(1/t-1)}$ when written as functions of the new couplings. With the currents normalized by their leading Wick vertices in that new convention, the two canonical endpoint coefficients also acquire $e^{2c}$ and $e^c$, respectively. At a fixed physical action both endpoint amplitudes then transform by the same factor $e^{c/t}$. Keeping the numerical endpoint coefficients unchanged while changing the Wick convention would not preserve the transported Ward prescription.

There is one definite wording correction: V21 is the **compensating coupling path** accompanying a translated vacuum cycle. A direct pullback of that translated cycle gives the opposite signs. The exact change-of-variables identity below resolves the distinction without changing the claimed endpoint characters.

The resulting reference operators can be scalar identities within each particle-vacuum sector while remaining nonconstant diagonal operators on the sum of sectors. This supplies a precise conditional sector statement. It does not construct the original Toda interface sectors or its full all-order current endpoints.

## Verdicts and independent evidence

| Item | Verdict | Independent evidence and restriction |
|---|---|---|
| V1, first formula | Correct, conditional on the stated reflection completion | Explicit rank-one specialization of Fateev (24)-(25), derived below; two irrational numerical probes |
| V1, second formula | Correct under the same conditional interpretation | Independent simple-zero expansion and radial identity-residue integral; no earlier script used |
| V13 | Correct | Direct substitution in the source exponential definition gives the three radial/angular labels |
| V14 | Correct for the stated rational form-factor family | Recursion induction below and fresh sign sums through six neutral particles |
| V15 | Correct | Original PDF p. 33 confirms shifted momentum in the coupling power and $R=A_{\rm reflected}/A_{\rm original}$ |
| V16 | Correct at generic coupling | Explicit gamma cancellation and $1/\Gamma(\zeta\delta)=\zeta\delta+O(\delta^2)$ |
| V17 | Correct | Gaussian contraction and ordinary two-dimensional radial integration, including the minus action insertion |
| V18 | Correct | $\Gamma(2-1/t)=(1-1/t)\Gamma(1-1/t)$, with the same logarithm branch |
| V19 | Correct in the zero-constant Wick convention | Independent momentum dot product, Taylor expansion and angular projection |
| V20 | Correct through the displayed order for $-1<t<0$ | The two radial primitives give precisely the two stated subtractions; its extension to another Gaussian convention requires the factors derived below |
| V21 | Correct as compensated transport; wrong as the stated direct pullback | Explicit change of integration variable gives the opposite pullback signs; replacement text below |
| V22 | Correct as a leading reference benchmark; unclear as a full original-theory endpoint | Reference reflection and crossing determine the scalar on a matched particle domain; higher endpoint terms and interfaces remain unconstructed |
| Scheme-independent use of the numerical gamma formulas without matching the Ward prescription | Wrong | General Gaussian transport gives explicit nontrivial factors |
| One scalar identity on the entire vacuum cover | Wrong for a nontrivial character | Explicit sector-projector realization below |

No theorem upgrade is recommended. The evidence is `referee_scratch/sl21_reflection_bridge_audit_20261024/check_bridge.py` and `.out`: **87 PASS, 0 FAIL**. This script was written from Gaussian contractions and the source equations, without importing or rerunning previous scripts. The derivations below, rather than the number of finite checks, establish the algebraic conclusions.

## 1. Source conventions actually checked

Fateev, [hep-th/0103014v2](https://arxiv.org/abs/hep-th/0103014), printed p. 3, equation (3), gives the exact normalization line

> $\langle\varphi_a(x),\varphi_b(y)\rangle=-\delta_{ab}\log|x-y|^2.$

Printed p. 6, equation (24), gives

> $R_s(a)=N(a)/N(s(a))=A_{s(a)}/A_a.$

Printed p. 33, equation (131), gives

> $R_s(a,\eta)=A_{s(a),\eta}/A_{a,\eta}.$

The same page's equation (132) uses $\bar a=a-Q$ in the coupling powers as well as in the gamma factors. The actual page was visually inspected; the bar in that power is not an assumption from extracted text. The page also states exactly: “Until now we considered spinless fields and non-compactified field Φ.” Its formula is not by itself a choice of the Toda model's global fermionic sectors.

Lashkevich-Pugai, [arXiv:1610.04926](https://arxiv.org/abs/1610.04926), equations (1.6)-(1.8), use kinetic coefficient $1/(8\pi)$, $\alpha^2=\beta^2+1=1/(2g)$, and $V_{a,\ell,q}=\exp[ia\vartheta/(2\alpha)-\ell\phi/(2\beta)-iq\widetilde\vartheta/(2\alpha)]$. The wall in (1.6) is exactly $-\lambda e^{-2\beta\phi}$; the report's $+\lambda_1e^{-2\zeta\phi}/\pi$ therefore requires $\lambda=-\lambda_1/\pi$. The sign cannot be removed inside a fractional power. Before (2.16) the source says the fields coincide “up to a normalization factor.” After (1.13) it says: “so that operators can be renormalized by admixture of operators of lower dimensions.” These exact qualifications are pertinent to the identity discussion.

Fateev, [SIGMA 13 (2017) 080](https://sigma-journal.com/2017/080/sigma17-080.pdf), printed p. 6, labels (4.1) “first hierarchy” and writes the interaction $2\mu e^{b\varphi}\cos(a\phi)+\mu_1e^{-2b\varphi}$. This corroborates the local hierarchy identification after its $1/(16\pi)$ field normalization is converted. It does not supply the missing all-order current normalization in the Toda variables.

Local copies and extracted source lines are in `referee_scratch/sl21_vacuum_limits_20261024/`; this audit only read them. In particular, see `fateev_2001.txt` lines 121-124, 261-268 and 1607-1649, and `lashkevich_pugai_2017.txt` lines 80-131 and 316-326. The normalization comparison uses only these primary sources and direct computations.

## 2. Re-derive the two finite reflection values

Write $z=\zeta>0$, $v=z^2=-t/2$, $\alpha^2=1+v$, and use interacting-plane momenta

\[
p_0=(z,i\alpha),\quad p_2=(z,-i\alpha),\quad p_1=(-2z,0),
\quad p_A=(-2z-z^{-1},0),\quad p_U=(z+z^{-1},-i\alpha).
\tag{R1}
\]

The spectator component is zero. Let $\gamma_f(x)=\Gamma(x)/\Gamma(1-x)$. All displayed fractional powers mean a continuous logarithm on a specified branch.

For the wall, the rank-one long root has squared length two. Set $b_L=-\sqrt2z$, $Q_L=(b_L+b_L^{-1})/\sqrt2=-z-(2z)^{-1}$. The source's function becomes

\[
A_L(a)=[\lambda_1\gamma_f(-t)]^{(a-Q_L)/(\sqrt2b_L)}
\Gamma\!\left(1-\frac{\sqrt2(a-Q_L)}{b_L}\right)
\Gamma\!\left(1-\sqrt2(a-Q_L)b_L\right).
\]

The endpoint exponent is $2Q_L$. The ratio $A_L(0)/A_L(2Q_L)$ has coupling exponent $1/t-1$, numerator gamma arguments $2-1/t,2-t$, and denominator arguments $1/t,t$. Thus

\[
G_A^0=[\lambda_1\gamma_f(-t)]^{1/t-1}
\frac{\Gamma(2-t)\Gamma(2-1/t)}{\Gamma(t)\Gamma(1/t)}.
\tag{R2}
\]

For the short-root sine-Liouville reflection use $Q=(2z)^{-1}$, $a_*=z+z^{-1}$, $\xi=-\alpha$, $C=\sqrt{\lambda_0\lambda_2}/(-t)$. The source reduction is

\[
A_S(a,\xi)=C^{(a-Q)/z}
\frac{\Gamma(1-2z(a-Q))\Gamma(1-(a-Q)/z)}
{\Gamma(1/2-z(a-Q)+\alpha\xi)\Gamma(1/2-z(a-Q)-\alpha\xi)}.
\tag{R3}
\]

At $a=a_*+\delta$, the reflected radial exponent is $-z-\delta$. Its two denominator arguments are $z\delta$ and $2-t+z\delta$. Its numerator at $\delta=0$ contains $\Gamma(2-t)\Gamma(2-1/t)$, canceling the regular denominator's $\Gamma(2-t)$. The unreflected denominator arguments at zero are (t-1,1), while its numerator is $\Gamma(t)\Gamma(1/t)$. Consequently

\[
\frac{A_S(2Q-a,-\alpha)}{A_S(a,-\alpha)}
=\delta z C^{2/t-2}\frac{\Gamma(2-1/t)}{(t-1)\Gamma(1/t)}+O(\delta^2).
\tag{R4}
\]

The reflected field has dot product $1-z\delta$ with $p_0$, and their sum is $(-\delta,0)$. Therefore the first action insertion is

\[
-\frac{\lambda_0}{\pi}\int_{|x|<R}d^2x\,|x|^{-2+2z\delta}e^{-\delta\phi}
=-\frac{\lambda_0R^{2z\delta}}{z\delta}e^{-\delta\phi}.
\tag{R5}
\]

The residue is $-\lambda_0/z$ times the identity. Multiplying it by (R4) gives

\[
G_U^0=\frac{\lambda_0}{1-t}C^{2/t-2}
\frac{\Gamma(2-1/t)}{\Gamma(1/t)}.
\tag{R6}
\]

At $\lambda_0=\lambda_2=\lambda$, the gamma recurrence gives $G_U^0=C^{2/t-1}\Gamma(1-1/t)/\Gamma(1/t)$. No additional sign appears. In particular, replacing $\pi\mu=\lambda_1$ by its absolute value or dropping the source's negative wall sign would change the analytic branch.

The checks use $t=-\sqrt2/3,-\sqrt3/4$. At $\delta=10^{-14}$, the directly evaluated reflection coefficient divided by (R4)'s linear term is respectively

```
1.0000000000001833009238
1.00000000000006860078839
```

The convergence is linear in $\delta$. Independent radial integration and the source-ratio evaluation agree with (R2) and (R6).

These statements use nonzero couplings, generic irrational $-1<t<0$, the specified logarithms, and zero finite $P$ coupling. The scalar expressions can be continued meromorphically further, but the displayed two-term cutoff subtraction cannot. At exceptional couplings, notably $1/t\in\mathbb Z_{<0}$, the apparent zeros of the gamma ratios do not license a value for a resonant limiting operator; the family needs a joint limiting prescription. None is constructed here.

### Uniqueness of the generic local resonance

Let $n_0,n_2,m$ count $Y_0,Y_2,V_1$ insertions and let $l\geq0$ be the total descendant level. The coupling dimensions are $1,1,2+4v$. The reflected field has dimension one. With $n=n_0+n_2$ and $D=n_0-n_2-1$, dimension equality gives

\[
D^2+n+2m+l+v[D^2-(n-1-2m)^2+4m]=1.
\]

Irrational $v$ separates the integer coefficients. The constant equation forces $n+2m+l\leq1$. The only nontrivial insertion is $n_0=1,n_2=m=l=0$, as used in (R5).

For the unreflected $U$, the corresponding two equations are

\[
D^2-n+6m+l=1,\qquad D^2-(n+1-2m)^2+4m=0.
\]

Put $r=n-2m$. Eliminating $D^2$ gives $r(r+1)+l=0$. Hence $r=0$ or (-1) and $l=0$. For $r=0$, the second equation gives $D^2=1-4m$, forcing $m=n=0$. For $r=-1$, it gives $D^2=-4m$, forcing $m=0,n=-1$, which is impossible. Thus there is no nontrivial perturbative local resonance for $U$ in this class. A finite $P$ coupling is excluded from that class, as it is from the reference prescription.

## 3. Cutoff matching and what it does not fix

The Gaussian product has $p_U\cdot p_0=2-t$ and $p_U+p_0=-p_A=p_F$. The free Taylor expansion and angular projection give

\[
U(0)Y_0(x)=|x|^{2t-4}
\left(F+xCF+\bar x\bar CF+x\bar x C\bar CF+\cdots\right).
\]

The two angularly invariant integrals have primitive coefficients $-\lambda_0/(t-1)$ and $-\lambda_0/t$. Their lower-limit terms are therefore $+\lambda_0\epsilon^{2t-2}F/(t-1)$ and $+\lambda_0\epsilon^{2t}C\bar CF/t$. Subtracting them gives exactly V20. The next angular term has lower-limit exponent $2t+2>0$. This proves the stated matching through one $Y_0$ insertion, or through mass order two after multiplication by $t^2\lambda_2$.

This is a comparison of local ultraviolet extensions, including their finite convention. It does not show that every possible all-order Ward endpoint equals this primary. A general finite identity admixture on a massive branch is not excluded merely by the absence of a perturbative resonance. For example, adding $\eta_U G_U^0(\lambda)I$ changes the vacuum component while leaving all positive-particle form factors zero; the coefficient has the required engineering dimension and branch character. It is nonanalytic in the couplings and is outside the perturbative mixing class just analyzed.

The stronger meromorphic reflection prescription does exclude that change if its ultraviolet normalization and identity-pole residue are held fixed. A regular finite addition to the **reflected** family is annihilated by the zero in (R4); changing the finite result requires changing that pole residue or changing the definition of the unreflected analytic family. Thus the resonance calculation fixes a definite reference normalization, not a universal normalization determined by the positive-particle form factors alone.

## 4. General finite Gaussian transport

Let $H$ be a constant symmetric finite subtraction matrix. This is a change of Wick convention, not a change of the singular covariance. Define

\[
V_p^H=Z_pV_p^0,\qquad Z_p=\exp[-\tfrac12p^THp].
\tag{R7}
\]

For the same physical action the couplings obey $\lambda_i^0=Z_i\lambda_i^H$. The OPE coefficient becomes

\[
V_p^H(x)V_q^H(0)\sim e^{p^THq}|x|^{-2p\cdot q}V_{p+q}^H(0).
\tag{R8}
\]

This follows directly from $Z_pZ_q/Z_{p+q}=e^{p^THq}$. A constant subtraction contributes no derivative contraction, so the same $Z_{p+q}$ multiplies the derivative descendants in V19.

Write $\mathcal G_A(\lambda_1)$ and $\mathcal G_U(\lambda_0,\lambda_2)$ for the right-hand sides of (R2) and (R6). Then the reference vacuum values in the new basis are

\[
\boxed{G_A^H=Z_A Z_1^{1/t-1}\mathcal G_A(\lambda_1^H),\qquad
G_U^H=Z_U Z_0^{1/t}Z_2^{1/t-1}\mathcal G_U(\lambda_0^H,\lambda_2^H).}
\tag{R9}
\]

This is the exact covariance law $G_p^H(\lambda^H)=Z_pG_p^0(Z_i\lambda_i^H)$. Logarithms must be transported continuously; arbitrary independent principal-value resets would invalidate it.

Let $E_{U0}=e^{p_U^THp_0}$. The endpoint's subtraction in the new basis is

\[
U_{\min}^H=\lim\left[U^H-
\frac{\lambda_0^HE_{U0}}{t-1}\epsilon^{2t-2}F^H-
\frac{\lambda_0^HE_{U0}}{t}\epsilon^{2t}C\bar C F^H\right]
=Z_UU_{\min}^0
\tag{R10}
\]

through the tested order. Likewise the zero-finite-part action contact becomes

\[
\kappa^H=-\frac{\lambda_0^H\lambda_2^H}{t}
E_{02}\epsilon^{2t},\qquad E_{02}=e^{p_0^THp_2}.
\tag{R11}
\]

These factors are necessary even if all singular exponents and the displayed gamma ratios are unchanged.

### Match the normalization of the currents

In this plane $T=(2z,0)$, $d=(z^{-1},0)$, and the currents have leading vertices $J_A=e^{-dX}$, $J_B=C e^{dX}$. Split the finite subtraction as $H=H_L+H_R$, with no additional finite mixed chiral contraction. The two left currents both acquire $Z_{d,L}$, and the right currents acquire $Z_{d,R}$. Their mixed products acquire $Z_d=Z_{d,L}Z_{d,R}$. The derivative $C$ introduces no extra factor because derivatives annihilate a constant covariance change.

Transporting the **same Ward prescription**, while fixing the leading coefficients of the currents in the new Wick basis to one, therefore requires

\[
\boxed{B_{AA}^{H,\mathrm W}=\gamma_E^2\lambda_1^H
\frac{Z_1Z_d}{Z_A}F_A^H,\qquad
B_{BB}^{H,\mathrm W}=t^2\lambda_2^H
\frac{Z_2Z_d}{Z_U}U_{\min}^H.}
\tag{R12}
\]

Their extra coefficients reduce to

\[
\frac{Z_1Z_d}{Z_A}=e^{d^THT},\qquad
\frac{Z_2Z_d}{Z_U}=e^{d^THp_2}.
\tag{R13}
\]

This can also be read directly from the chiral current/source OPE. The first $B$-current tail, for example, contains the factor $Z_{d,L}Z_{2,L}/Z_{d+p_2,L}=e^{d^TH_Lp_2}$. Both methods agree.

If instead one keeps the **original physical currents** fixed and merely expresses their endpoints in a new basis, remove $Z_d$ from (R12). The endpoints themselves then remain unchanged. Conversely, simply setting the coefficients in front of $\lambda_1^HF_A^H$ and $\lambda_2^HU_{\min}^H$ back to $\gamma_E^2,t^2$ is an additional finite change of endpoint prescription; it is not either of these transports.

### The physically relevant change and a distinct bosonized change

Changing the finite subtraction by $c$ only for the two physical bosons keeps the elementary fermion OPE fixed. In the original three-field basis this is $H_{\rm phys}=\operatorname{diag}(c,c,0)$. The orthogonal transformation in V10 projects it on the interacting plane to

\[
H_{\rm eff}=\operatorname{diag}\left(c,\frac{cz^2}{1+z^2}\right).
\]

There are also spectator cross components, but the displayed vertices have zero spectator momentum. Consequently

\[
Z_0=Z_2=1,\quad Z_1=e^{ct},\quad Z_A=e^{c(1-t)^2/t},\quad
Z_U=e^{c(1/t-1)},\quad Z_d=e^{c/t}.
\tag{R14}
\]

Equations (R9)-(R13) give

\[
G_A^{\rm phys}=e^{c(1/t-1)}\mathcal G_A(\lambda_1^{\rm phys}),\qquad
G_U^{\rm phys}=e^{c(1/t-1)}\mathcal G_U(\lambda_0^{\rm phys},\lambda_2^{\rm phys}),
\]
\[
E_{U0}=e^{c(1-t)},\quad E_{02}=e^{-ct},\quad
B_{AA}^{\rm phys,W}=e^{2c}\gamma_E^2\lambda_1^{\rm phys}F_A^{\rm phys},\quad
B_{BB}^{\rm phys,W}=e^ct^2\lambda_2^{\rm phys}U_{\min}^{\rm phys}.
\tag{R15}
\]

A different change $H=cI$ on **all three bosonized fields** also changes the fermion vertex normalization. It gives $Z_0=Z_2=e^{c/2}$, $Z_U=e^{c(1/t-1/2)}$, $G_A^H=e^{c(1/t-1)}\mathcal G_A$, $G_U^H=e^{c(2/t-1)}\mathcal G_U$, $E_{U0}=e^{c(2-t)}$, and $E_{02}=e^{c(1-t)}$. Its Ward coefficients are again $e^{2c},e^c$. It must not be confused with the physical-boson subtraction constant in V9.

For either change, and indeed general $H$, (R12) implies at the fixed action

\[
\mathcal A^{H,\mathrm W}=Z_d\mathcal A^0,\qquad
\mathcal B^{H,\mathrm W}=Z_d\mathcal B^0.
\tag{R16}
\]

For the physical-boson change, holding the *new coupling coordinates* fixed instead gives

\[
\mathcal A^{\mathrm{phys,W}}(\lambda^H)
=e^{c(1+1/t)}\mathcal A_{\mathrm{formula}}(\lambda^H),\qquad
\mathcal B^{\mathrm{phys,W}}(\lambda^H)
=e^{c/t}\mathcal B_{\mathrm{formula}}(\lambda^H).
\]

Their ratio is therefore $e^{-c}$ times the zero-constant formula evaluated at $\lambda^H$. Fixed-action invariance uses $\lambda_1^H=e^{-ct}\lambda_1^0$ and $\lambda_{0,2}^H=\lambda_{0,2}^0$; it does not assert that one unchanged formula in the new coordinates is invariant.

Thus the endpoint ratio is insensitive to this finite Gaussian transport at fixed action. On the positive real coupling branches in $-1<t<0$, its explicit value simplifies to

\[
\boxed{\frac{\mathcal B^0}{\mathcal A^0}
=-\left[-\gamma_f(t)\frac{\lambda_0\lambda_2}{\lambda_1}\right]^{1/t}.}
\tag{R17}
\]

The derivation uses $(1-t)^2\lambda_2G_U^0/(\lambda_1G_A^0)$ and $\gamma_f(t)\gamma_f(-t)=-t^{-2}$. This ratio has mass dimension two. For general complex coupling branches, retain $(1-t)^2\lambda_2G_U^0/(\lambda_1G_A^0)$ with the original logarithms instead of arbitrarily recombining powers. It is a comparison invariant under the specified Gaussian transport, not under unrelated independent rescalings of the two charge pairs, and it has not been identified with a physical mass.

## 5. Coupling branch and cycle transport

For any field displacement $\Delta$, direct substitution in a regulated integral gives

\[
G_{\mathcal C+\Delta,p}(\lambda_i)
=e^{p\cdot\Delta}G_{\mathcal C,p}(\lambda_i e^{p_i\cdot\Delta}).
\tag{R18}
\]

This statement is exact whenever the two integrals exist and the regulator and measure permit the stated translation. It does not require a contour deformation through a divergent domain.

For $\Delta=(i\pi/z,\pi/\alpha)$, the screening pairings are

\[
p_0\Delta=2\pi i,\quad p_2\Delta=0,\quad p_1\Delta=-2\pi i.
\]

Thus a literal pullback in (R18) has log changes $+2\pi i,0,-2\pi i$. V21 has their opposites. Those opposites arise from the compensated path

\[
\mathcal C_s=\mathcal C_0+s\Delta,\qquad
\lambda_i(s)=\lambda_i(0)e^{-s p_i\Delta},\qquad 0\le s\le1.
\tag{R19}
\]

Indeed $S_{\lambda(s)}(\Phi+s\Delta)=S_{\lambda(0)}(\Phi)$, so transport along this path gives $G_s(V_p)=e^{s p\Delta}G_0(V_p)$. At its end the numerical couplings return to their starting values while their logarithms and the cycle have moved. This yields $G_A\mapsto s_0^{-1}G_A$, $G_U\mapsto s_0G_U$, with $s_0=e^{-2\pi i/t}$. The powers in (R2) and (R6) reproduce these factors exactly.

For the compensating path, a finite-regulator convergence construction can be inherited by translation because the pulled-back action is constant. This still does not establish the infinite-volume continuum limit or an analytic continuation on a fixed original real cycle. The path specifies what the branch assertion means; it does not replace those analytical questions.

**Replacement at `briefs/referee_sl21_vacuum_limits_report.md`, line 335:**

> Under the period $\Delta\phi=i\pi/\zeta$, $\Delta\vartheta=\pi/\alpha$, transport the cycle as $\mathcal C_s=\mathcal C_0+s\Delta$ and compensate its couplings by $\lambda_i(s)=\lambda_i(0)e^{-s p_i\cdot\Delta}$. The pulled-back action is then constant. The resulting closed coupling path has the logarithm shifts (V21), and its transported vacuum acquires the endpoint characters stated below. A direct pullback at fixed couplings instead has the opposite logarithm shifts, according to $G_{\mathcal C+\Delta,p}(\lambda)=e^{p\cdot\Delta}G_{\mathcal C,p}(\lambda_i e^{p_i\cdot\Delta})$.

No edit has been applied.

## 6. Identity-like form factors and the matched sector statement

The source labels follow directly from (R1): $F_A=V_{0,2h-2,0}$, $U=V_{-h,-h,0}$, $F=V_{0,2-2h,0}$, where $h=2-t=g^{-1}$. Its $p_\pm=1-h+\ell\mp a$ therefore give $p_A=-(1-h,1-h)$ and $p_U=(1-h,1-h-2h)$.

Here is an independent recursion proof of the required zero. Put $p_r=(1-h+r,1-h-r)$. For $K=0,L>0,r=L$, the first factor in the source product (3.6) is $2\cos(-\pi/2)=0$. In (3.3), removal of a particle from $J_{K,L,p_r}$, $r=L-K$, sends the constant term to $(K-1,L,p_{r+1})$ and each residue term to $(K-1,L-1,p_r)$. Strong induction annihilates both, except possibly the residue for $K=L=1$. In that case $p_+=p_-$, and the two phases in (3.2) coincide, so its residue is zero as well. Hence $J_{K,L,p_{L-K}}=0$ for the stated range. Inversion (3.8) gives the $A$ case; shifting $p_-$ by (-2h) gives the $U$ case via (2.17). Direct sign sums from (2.7)-(2.12), evaluated independently through six neutral particles, agree with maximum absolute residual below $10^{-67}$.

These are identities of the normalized bootstrap functions. To turn them into an operator statement, specify the following hypotheses:

1. A selected reference massive sector has a regular meromorphic family with the Wick normalization and reflection relations above, including the fixed identity-pole residue (R5). Its finite $F_A,U_{\min}$ representatives are (R2),(R6), without a singular normalization times a vanishing form factor retained as a separate limiting field.
2. The $q=0$ fields obey charge selection, the displayed vacuum-to-particle form factors are their full such matrix elements, and crossing reconstructs all particle-domain matrix elements, including disconnected pieces.
3. The finite-particle domain is dense in that sector; no additional invisible zero-energy or superselection states are being identified with that domain.

Under these hypotheses, for all vectors in the common particle domain,

\[
F_A=G_A^0 I,\qquad U_{\min}=G_U^0 I
\tag{R20}
\]

as quadratic forms in that sector. The proof is that every nonvacuum connected form factor is zero, while the disconnected crossings give precisely the vacuum constant times the state inner product. Reflection fixes those constants independently of that zero statement. If one defines the reference endpoint operators by the transported formulas (R12), their scalar values are exactly (R16); the local canonical endpoint agrees with them through the order established in section 3. Claiming this equality for the full Toda Ward endpoint requires an additional all-order matching result.

Now take translated particle sectors $\mathcal H_n$ with projectors $P_n$, and let $\mathcal H=\bigoplus_n\mathcal H_n$. Equivariance gives the concrete diagonal operators

\[
\mathbb U=G_U^0\sum_n s_0^nP_n,\qquad
\mathbb A=G_A^0\sum_n s_0^{-n}P_n.
\tag{R21}
\]

They are scalar identities on each summand and have the nontrivial characters on the cover. For real generic $t$, the phases have modulus one, so the scalar diagonal operators are bounded if their common amplitudes are finite. If $S$ sends $\mathcal H_n$ to $\mathcal H_{n+1}$, then $S^{-1}\mathbb US=s_0\mathbb U$ and $S^{-1}\mathbb AS=s_0^{-1}\mathbb A$. This construction explains compatibility of identity-like particle form factors and nontrivial vacuum characters. It does not construct $S$ from a conserved Toda charge.

If interface sectors $\mathcal H_{l,r}$ and scalar endpoint limits with those left/right asymptotic values exist, put $N=r-l$. Direct subtraction gives

\[
\Delta B_{AA}=\mathcal A^0s_0^{-r}(1-s_0^N),\qquad
\Delta B_{BB}=\mathcal B^0s_0^r(1-s_0^{-N}).
\tag{R22}
\]

Homogeneous particle completeness alone supplies no such interface sector or common nonlocal-charge domain. It also does not prove that the all-order Ward endpoint equals the particular reference operator (R12). These two precise extra requirements, rather than any contradiction between (R20) and a nontrivial character, remain the global boundary gap.

## Integration

None performed. **Verified:** source orientation and gamma factors; the compensating zero/pole; generic resonance uniqueness in the stated perturbation class; normalized form-factor degeneracy; the displayed cutoff subtractions. **Modified:** the interpretation of V21 must be compensated cycle/coupling transport, not the asserted direct pullback. **Extended:** the general finite Gaussian law, its current/Ward normalization factors, the physical-boson-only specialization, the invariant endpoint ratio, and a conditional sector-projector realization. **Rejected:** transporting only the numerical gamma formula while leaving OPE/Ward coefficients unchanged, or identifying a charged endpoint with one constant identity over all branches. **Unresolved:** the original continuum completion, all-order endpoint matching, interface sectors and the common charge product.

The script verifies that all eight recorded protected inputs retain their hashes. The old scripts were neither imported nor rerun. No book file or existing report was changed.
