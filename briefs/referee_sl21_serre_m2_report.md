# Root-system sl(2|1) Toda: mass-squared Serre test and its quantum normalization

Date: 2026-10-24. Model: Codex (GPT-6).

Commit reviewed: `745586e6a6f7bca5ed91a1e9cf89d6ef5df387f0`.

Task: independently audit the first-order contour argument, then include every order-mass-squared channel. `AGENTS.md` was read. No existing file was edited, no subagent was spawned, and no git mutation, commit, push, or book render was performed. New work is confined to `referee_scratch/sl21_serre_m2_20261024/` and this report. No previous check script was imported or rerun.

## Outcome and scope

Both quartic relations extend through mass order two in an explicitly continued contour-product prescription:

\[
B^3A-c_sB^2AB+c_sBAB^2-AB^3=O(m^3),
\]
\[
A^3B-c_sA^2BA+c_sABA^2-BA^3=O(m^3),
\qquad c_s=1+s+s^{-1},\quad s=e^{-i\pi\rho},\quad \rho=2/t.
\tag{M1}
\]

These are **formal finite-interval contour identities**, including action insertions, tails, and product subtractions. The words denote the four-fold products defined below. They are not a construction of an associative global charge algebra or of interacting Hilbert-space operators. Physical boundary closure, a quantum-group presentation, and all-order Serre relations remain open.

Two additional steps were essential. First, an independent exponent is needed for the unlike Yukawa pair: the previous starting domain does not make its integral convergent. Second, the even-root contribution to $AAAB$ has a coupling-independent marginal scalar channel. It cannot be dismissed by the positive power count for the double-Yukawa primitive. Its scalar kernel and all its tails are exact derivatives in the minority $B$ coordinate; this removes the remaining continuation pole.

The quantum normalization was checked again. The regulated double-Yukawa contact and the composite-potential response combine into the derivative of the current counterterm, with the full fermion-current contribution retained. No finite independent composite coupling is added.

There are **267 new PASS, 0 FAIL** checks. The conclusion uses the written continuation and collision argument as well as these computations; no full multidimensional massive correlation function was numerically integrated. No [Theorem] label or symmetry-algebra identification is recommended.

| Item | Verdict | Scope/evidence |
|---|---|---|
| First-order removable chamber pole and shuffle coefficient | **Correct in the earlier analytic prescription** | Fresh Wick derivative and exact labelled-permutation audit |
| Reusing its starting domain for the unlike double-Yukawa integral | **Wrong without an extra continuation step** | That domain has $\operatorname{Re}t<0$ |
| Action and current subtraction coefficients | **Correct** | New radial moments, nonconstant-background integrals, and the local Ward polynomial |
| Complete second-order separated Wick kernels | **Correct** | Both chiralities of vertex factors, all tail placements, and cocycles evaluated in actual operator order |
| Absence of a new marginal even-root channel | **Wrong** | The raw full $AAAB/V_1$ degree is $3(g-\rho)$ |
| Removal of that marginal channel | **Correct in the stated extension** | Exact minority-coordinate derivative, including all tails and moving-hole behavior |
| Identities (M1) | **Correct as formal identities in this prescription** | Regular diagonal restriction of the complete left-hand family |
| Unrenormalized products, exceptional couplings, or a physical symmetry algebra | **Unclear** | Not established here |

## 1. Inputs and a fresh local Ward calculation

Use $X(z)X(w)\sim-\log(z-w)$ and

\[
a_0=(ib/\sqrt2,b/\sqrt2,i),\quad
a_2=(ib/\sqrt2,-b/\sqrt2,-i),\quad
t=b^2,\quad T=a_0+a_2=-td,\quad k=d+T.
\]

Direct scalar products give $a_i^2=-1$, $a_0a_2=1-t$, $a_id=1$, and $d^2=-2/t=-\rho$. All lengths below are measured in the common Wick-ordering reference length. Restoring that length also restores the dimensionful factors in the auxiliary regulators.

The physical representatives, before product extension, are

\[
J_A=:e^{-dX}:,\qquad J_B=:C e^{dX}:,\qquad C=a_0\partial X,
\]
\[
K_2=:e^{(d+a_2)X+a_2\bar X}:,\quad
E=:e^{kX+T\bar X}:,\quad L=:CE:,
\]
\[
H_B=-t\lambda_2K_2+\lambda_0\lambda_2\epsilon^{2t}L,
\qquad H_A=-\gamma\lambda_1F_A,\quad
F_A=:e^{-kX-T\bar X}:,\quad \gamma=\frac{t}{1-t}.
\tag{M2}
\]

The first two couplings $\lambda_0,\lambda_2$ have mass order one. The independently Wick-normalized $\lambda_1$ has mass order two. The action is

\[
S_{\mathrm{int}}=\pi^{-1}\int d^2w\,
(\lambda_0Y_0+\lambda_2Y_2+\lambda_1V_1+\kappa P),
\quad P=:e^{TX+T\bar X}:,\quad V_1=:e^{-TX-T\bar X}:,
\]
\[
\kappa=-\frac{\lambda_0\lambda_2}{t}\epsilon^{2t},\qquad \kappa_R=0.
\tag{M3}
\]

The current-first cocycle for momentum $nd$ and one full Yukawa is $(-1)^n$. Full Yukawa fusion has coefficient $+1$. The currents and full tails are even. The common charge factors $(2\pi i)^{-n}$ are suppressed in the ordered-product formulas.

### The potential and double-Yukawa terms are combined before interpretation

Set $D=d\partial X$, $C_2=\partial C$, $D_2=\partial D$. In the two local source point splittings, the holomorphic linear fields are $C$ and $-tD-C$. Their quadratic difference is

\[
N/E=\tfrac12\{C_2+C^2-[-tD_2-C_2+(-tD-C)^2]\}
=C_2+\frac t2D_2-tCD-\frac{t^2}{2}D^2.
\tag{M4}
\]

Here $C$ includes the physical fermion bilinear $j=:\chi_0\chi_2:$. In the direct fermionic expansion the terms $\partial j=:\partial\chi_0\chi_2:+:\chi_0\partial\chi_2:$ survive in precisely this difference; they are not discarded as a bosonic approximation.

The boundary normalization can be recovered directly. After subtracting the ordinary action dressing of $\partial H_B^{(1)}$, the two current/source contacts leave $(\lambda_0\lambda_2/\pi)\int d^2r\,(\partial_z\Theta)(F_2-F_0)$, where $F_i=-r^{t-2}\bar r^{t-1}G_iE$ and the quadratic holomorphic coefficients of $G_i$ are the two terms in (M4). The coefficient of $F_2$ is $t-(t-1)=1$; the coefficient of $F_0$ is $-1$. With $r=\epsilon e^{i\theta}$, $\partial_z\Theta=-e^{-i\theta}\delta(|r|-\epsilon)/2$. The two minus signs give $\tfrac12\epsilon^{2t-2}\int d\theta\,e^{-2i\theta}(G_2-G_0)E=\pi\epsilon^{2t}N$ at the displayed order. This fixes the positive $\lambda_0\lambda_2\epsilon^{2t}N$ term and its ordinary-area normalization independently of the potential response.

The triple, double, and simple current-first poles of $J_BP$ have coefficients

\[
tE,\quad (C+tD)E,\quad
\left[C_2+CD+\tfrac t2(D_2+D^2)\right]E.
\]

Apply the derivatives of the delta function from all three poles, using the action sign $-\kappa/\pi$. Their response is

\[
R_P/E=-tCD+\frac{t^2}{2}D_2-\frac{t^3}{2}D^2.
\tag{M5}
\]

Since $\partial E/E=(1-t)D$, direct subtraction gives

\[
\boxed{N-R_P/t=\partial L.}
\tag{M6}
\]

Thus the regulated Ward contact is

\[
\kappa R_P+\lambda_0\lambda_2\epsilon^{2t}N
=\lambda_0\lambda_2\epsilon^{2t}\partial L.
\tag{M7}
\]

This checks the current counterterm in (M2) with its coefficient and sign. It is not just equality modulo a total derivative. Direct angular quadrature of the exclusion-circle term selects (M4); its normalized errors at $\epsilon=.04,.02,.01$ are $5.27\times10^{-6}$, $1.32\times10^{-6}$, and $3.29\times10^{-7}$, with the expected two additional powers.

For the other current, $J_AV_1$ has double coefficient $F_A$ and simple coefficient $-DF_A$. Its distributional response is $tDF_A=-\gamma\partial F_A$, verifying the $H_A$ sign. $J_BV_1$ is regular: its vertex factor has a double zero and its derivative contraction costs only one power.

## 2. Continue the unlike-source collision explicitly

The earlier first-order construction used independent like-current exponent $g$, with unlike-current exponent $-\rho$ and $t=2/\rho$. Its initial domain $\operatorname{Re}\rho<-1$ implies $\operatorname{Re}t<0$. Hence the unlike-source pair cannot simply be integrated there without specifying an extension.

Introduce an auxiliary exponent $h$ by multiplying each pair of unlike full screening constituents by

\[
|w_0-w_2|^{2(h-t)}.
\tag{M8}
\]

This includes a constituent fused into a $K_2$ tail. It changes the pair factor to $|w_0-w_2|^{2h-2}$, while derivative contractions and normal fields still use $t=2/\rho$. Equal-node source factors remain $|w-w'|^2$. Keep the like-current regulator $(z_i-z_j)^{g-\rho}$ in every like-labelled current pair, also when a slot contains a tail. Neither regulator is a physical change to the action's root system.

The matched leading cutoff representatives are now

\[
\boxed{\kappa_h=-\frac{\lambda_0\lambda_2}{h}\epsilon^{2h},\qquad
H_{B,h}^{(2)}=\frac{t}{h}\lambda_0\lambda_2\epsilon^{2h}L.}
\tag{M9}
\]

To derive them, the source-pair scalar integral is

\[
\int_{\epsilon<|r|<R}d^2r\,|r|^{2h-2}
=\frac\pi h(R^{2h}-\epsilon^{2h}).
\]

Its coefficient in the double action insertion is $+\lambda_0\lambda_2/\pi^2$, so the missing lower-disk term is restored by the action weight $-\kappa_h/\pi$.

Independently, the tail/source product is

\[
K_2(z)Y_0(z+r)=r^{h-2}\bar r^{h-1}[E+rL+\cdots].
\tag{M10}
\]

The constant angular moment vanishes. The $rL$ term gives $\pi(R^{2h}-\epsilon^{2h})L/h$. Multiplying by the action correction to $-t\lambda_2K_2$ yields $t\lambda_0\lambda_2(R^{2h}-\epsilon^{2h})L/h$, proving the second term of (M9). The coefficient is $t/h$ before restriction, not one. Setting $h=t$ recovers (M2)–(M3) and the physical Ward calculation (M7). No Ward identity off $h=t$ is assumed or needed.

Nonconstant radial tests at $h=-\sqrt2/5$ verify both finite parts. For the action test the errors decrease from $0.11709$ at $\epsilon=.08$ to $0.005928$ at $.01$; for the tail test they decrease from $0.0013112$ to $0.00006643$. Both have the predicted next power $\epsilon^{2h+2}$.

Equations (M9)–(M10) display the leading scalar and current counterterms, not a complete bare-cutoff list for arbitrarily negative $h$. The full prescription is Taylor continuation of every proper face. For example, a term $r^n\bar r^n P_{nn}$ in the source-pair normal product has a subtraction proportional to $\epsilon^{2(h+n)}P_{nn}/(h+n)$. A tail/source term $r^{n+1}\bar r^n L_n$ has the corresponding coefficient $t\epsilon^{2(h+n)}L_n/(h+n)$. Other current faces are extended in the same way. No additional finite local coupling is selected at those faces. This prevents treating the two leading terms in (M9) as sufficient in regimes where higher descendants diverge.

## 3. The complete mass-squared expansion and its phases

For each four-current word the following channels are included. A tail placement replaces the indicated current slot; all such placements are summed. The source integrals use ordinary area.

| Channel | Coefficient in the perturbative product |
|---|---:|
| $Y_0Y_2$ in the action | $+\lambda_0\lambda_2/\pi^2$ |
| $Y_iY_i$ in the action | $+\lambda_i^2/(2\pi^2)$ |
| One $K_2$ tail and one $Y_i$ action insertion | $+t\lambda_2\lambda_i/\pi$ |
| Two $K_2$ tails in distinct slots | $+t^2\lambda_2^2$ per unordered slot pair |
| Composite $P$ in the action | $+\lambda_0\lambda_2\epsilon^{2h}/(\pi h)$ |
| One $L$ current counterterm | $+(t/h)\lambda_0\lambda_2\epsilon^{2h}$ |
| $V_1$ in the action | $-\lambda_1/\pi$ |
| One $F_A$ tail in an $A$ slot | $-\gamma\lambda_1$ |

The finite even-root coupling is independent here; no identity is obtained by imposing a quantum relation $\lambda_1=\lambda_0\lambda_2$. Lower-product subtractions are included with their ordinary dressing by the remaining interaction or tail.

A compact way to check all separated signs is to label a vertex by $(q,m_0,m_2)$, with left momentum $qd+m_0a_0+m_2a_2$ and right momentum $m_0a_0+m_2a_2$. The basic labels are

\[
A=(-1,0,0),\quad B=(1,0,0),\quad K_2=(1,0,1),
\quad L=(1,1,1),\quad F_A=(-1,-1,-1).
\]

The $B$ and $L$ vertices carry $C$; the others do not. For two ordered vertices define

\[
S_{ij}=m_{0i}m_{0j}+m_{2i}m_{2j}
+(h-1)(m_{0i}m_{2j}+m_{2i}m_{0j}).
\]

Their antiholomorphic exponent is $S_{ij}$, and their holomorphic exponent is

\[
\rho q_iq_j+(g-\rho)\mathbf1_{\text{like current labels}}
-q_i(m_{0j}+m_{2j})-q_j(m_{0i}+m_{2i})+S_{ij}.
\tag{M11}
\]

The ordered cocycle is $(-1)^{q_i(m_{0j}+m_{2j})}$. A derivative $C_i$ contracts with the second vertex with coefficient $-q_j+m_{0j}+(t-1)m_{2j}$, divided by $z_i-z_j$; derivative pairs give $(z_i-z_j)^{-2}$.

These formulas reconstruct every Gaussian Wick kernel in the table. In particular two $K_2$ tails have holomorphic exponent $g-1$, antiholomorphic exponent $1$, and the required cocycle. Their separated current exchange is $u=e^{-i\pi g}$, not a new odd sign. The $L$ and $F_A$ dressings also preserve the original current-label exchange phases. The unlike-source regulator contributes no extra current exchange phase.

The new script evaluates actual holomorphic and antiholomorphic powers in operator order, with nonconstant coherent backgrounds. It checks every labelled spatial permutation of both quartics, each nonempty channel in the table, and three exponent settings: one inside the initial convergence domain, one with $g$ and $h$ displaced from their physical values, and one on the physical surface at $t=\sqrt5$. The separated checks do not rely on setting a scalar background or on omitting fermionic cocycles.

Let $M(u,s)$ be the labelled-permutation matrix in the order $(BBBA,BBAB,BABB,ABBB)$. A fresh exact enumeration gives

\[
M(u,s)(1,-c_s,c_s,-1)^{\mathsf T}
=(0,-D(u,s),D(u,s),0)^{\mathsf T},
\]
\[
D(u,s)=\frac{(s-1)(s+1)(s-u)(s^2-u)(u+1)}{s^4},
\quad u=e^{-i\pi g}.
\tag{M12}
\]

The same matrix applies with $A$ and $B$ exchanged. Therefore the complete mass-squared defect is $D(u,s)$ times the difference of the two middle spatial-chamber kernels. This is a meromorphic relation, not yet a proof that its physical restriction vanishes: a chamber pole proportional to $(g-\rho)^{-1}$ could survive the zero of $D$.

## 4. Define the common continuation and enumerate its collision faces

Use identical finite spatial intervals and the same smooth infrared cutoff in all words. Keep the prior word-ordered time levels $c\epsilon(3/2,1/2,-1/2,-3/2)$, $c>2$, circular source holes at inner currents, and the outer circular Cauchy extension. Also exclude source pairs at radius $\epsilon$, including a source paired with a constituent in a tail. Use (M9) and the same local Taylor prescription for larger products.

Take the cutoff limit first in the common auxiliary domain

\[
\operatorname{Re}g>6,\qquad \operatorname{Re}\rho<-1,
\qquad \operatorname{Re}h>4,
\qquad t=2/\rho.
\tag{M13}
\]

Then continue the complete left-hand expansion in $(g,\rho,h)$ and restrict to

\[
g=\rho,\qquad h=t=2/\rho.
\tag{M14}
\]

All single-current Cauchy poles have the same distributional normalization as in section 1. Multiple-source and multiple-current collisions have positive relative degrees in (M13), apart from these specified one-current poles. Thus the cutoff limit there agrees with resolved Cauchy source kernels and common spatial chamber products. Source-pair circles can intersect current-hole boundaries; the joint collision estimates below include those intersections. They are not silently dropped.

The extension is fixed from an open initial domain of left-hand diagrams, not by imposing (M1). Setting $g=\rho$ or $h=t$ before taking that initial limit would define a different problem. At orders zero and one the new $h$ regulator is absent, so this prescription reduces to the preceding first-order analytic product.

### Relative dimensions and all source-containing subsets

For a cluster containing $a$ current labels $A$ and $b$ labels $B$, with at least one current, put

\[
\mathcal E(a,b)=g\left[\binom a2+\binom b2\right]-\rho ab.
\]

The common coordinate along the current interval does not scale. Two source integrations add four relative coordinates. Summing vertex powers, the maximal derivative contractions, and the relative measure gives the following bounds:

\[
\begin{array}{c|c}
\text{cluster/channel}&\text{relative degree}\\\hline
\text{one Yukawa}&\mathcal E+2a-b+1\\
\text{two unlike Yukawas}&\mathcal E+3a-2b+2h+1\\
\text{two equal Yukawas}&\mathcal E+3a-2b+5\\
V_1\text{ or an }F_A\text{ tail}&\mathcal E-a+2b+1
\end{array}
\tag{M15}
\]

A $K_2$/source channel has the same bound as its unfused two-source channel. Two $K_2$ tails have the equal-source bound. The $P$ and $L$ counterterm graphs, including their explicit $\epsilon^{2h}$, have the unlike-source bound. For an isolated $L$ with no other insertion the derivative cannot contract, so that raw bound is conservative. With no current at all, the two-source relative measure has dimension two, giving $2h$ for unlike sources and $4$ for equal sources; one must not use the current-anchored measure formula in that case.

The script generates these degrees directly from the ordered charge data (M11), separately for the double-source, tail/source, two-tail, potential, and current-counterterm graphs. It enumerates every subset of each quartic. Pure-current subsets use their ultraviolet extensions; one-source subsets use the independently audited first-order extension.

In particular the full-cluster degrees are

\[
\begin{array}{c|cc}
&BBBA&AAAB\\\hline
\text{unlike sources, including matched cutoff graphs}
&3(g-\rho)+2h-2&3(g-\rho)+2h+8\\
\text{equal sources/tails}&3(g-\rho)+2&3(g-\rho)+12\\
\text{even-root channel}&3(g-\rho)+6&3(g-\rho)
\end{array}
\tag{M16}
\]

The last entry is genuinely marginal on the physical surface. It is handled in section 5, not by a generic-coupling exclusion.

### A physical parameter patch where the other new faces are controlled

After meromorphic extension of the proper faces, use the open complex neighborhood of the real interval

\[
2<t<3,\qquad 2/3<\rho=2/t<1.
\tag{M17}
\]

Every new unlike- or equal-double-source degree in (M15) is strictly positive there. This is not a numerical scan: after substituting $g=\rho=2/t$, $h=t$, multiplying by $t>0$ gives a quadratic polynomial. The new script checks nonnegative Bernstein coefficients on $t\in[2,3]$, with strict positivity in the interior, for every allowed subset. The full unlike $BBBA$ degree is $2t-2>2$, and the equal-source full degree is two.

The existing pure-current $AB$ and first-Yukawa $BB$ extensions remain necessary. Their pole hyperplanes do not pass through (M17). The first-order full $BBBA/Y_2$ pole is removable, as audited below. In the even-root sector there is one new proper pair subtraction, one fixed single-$A$ Cauchy pole, and the full marginal $AAAB$ channel. Other even-root faces have positive degrees in the patch.

Local distance sectoring treats nested and disjoint subsets of these graphs. The proper Taylor remainders have nonnegative extra powers. Fixed-order Cauchy differentiation, including differentiation of exclusion masks, preserves the displayed net homogeneity. In (M13) the new source-pair factors give sufficient positive margins for intersections with current holes; equal-source factors are polynomially soft. After continuation, the positive new degrees in (M17) control the residual full clusters once the inherited proper faces have been extended. A cluster pinned at a spatial endpoint has one additional small coordinate.

The resulting Mellin poles have equations $\omega+j=0$, $j\ge0$. Except for the two derivative channels treated below, none contains an open part of the physical surface in this patch. Further resonant couplings may occur on continuation outside the patch and are excluded from the generic statement.

## 5. Audit the inherited marginal face and calculate the new one

### An extra source does not undo the first-order removal

The first-order scalar $BBBA/Y_2$ coefficient is an $A$ derivative. A fresh reconstruction from its three derivative Wick factors verifies this identity without using an old script. In a middle spatial chamber, the minority coordinate lies between two majority currents. In the initial negative-$\rho$ domain the primitive vanishes at both internal endpoints. Integrating by parts against an arbitrary smooth spectator moves the derivative to that spectator, improving the joint degree from $3(g-\rho)$ to $1+3(g-\rho)$. This is an identity of the continued proper-face extensions. It does not set the unextended positive-$\rho$ boundary value to zero.

The unlike additional Yukawa requires a uniform nested estimate. Let $r$ be the size of the first cluster and $R\gg r$ its distance to the other source. The first cluster fuses to a field with an unlike-source factor of total degree $2h-4$. Differentiating that factor costs one $R$ power. After the first-order improvement, the nested radial bound is proportional to

\[
\int_0^{R_0}R^{2\operatorname{Re}h-5}R\,dR
\int_0^R r^{3\operatorname{Re}(g-\rho)}\,dr.
\tag{M18}
\]

Its denominators correspond to $1+3(g-\rho)$ and $3(g-\rho)+2h-2$, not $3(g-\rho)$. The latter full degree is positive in (M17). For an equal additional source, the fused factor has degree zero, and the same calculation gives full degree $3(g-\rho)+2$. Thus the first-order proper face does not hide a new pole when its spectator becomes the second source. The region $R=O(r)$ is already the full cluster in (M16).

The node-0 first-order scalar vanishes; its leading remaining jet has the same positive improvement. Source-free full-current collisions have the audited Wick cancellations and degree at least two. These observations cover the inherited potentially dangerous proper faces.

### The even-root scalar primitive, including all tails

Let the three $A$ positions be $a_1,a_2,a_3$ and the minority $B$ be at $z$. Define

\[
\mathcal R_A=\prod_{i<j}(a_i-a_j)^g\prod_i(a_i-z)^{-\rho},
\qquad G_A(w)=\frac{(z-w)^2}{\prod_i(w-a_i)^2}.
\]

The scalar part of the single $B$ derivative contraction with the other fields is

\[
f_B=\sum_i\frac1{z-a_i}-\frac{t}{z-w}.
\]

Consequently the actual source Wick kernel obeys

\[
\boxed{\mathcal R_A G_A f_B=-\rho^{-1}\partial_z[\mathcal R_A G_A].}
\tag{M19}
\]

This uses $\partial_z\log\mathcal R_A=-\rho\sum_i(z-a_i)^{-1}$ and $\rho t=2$, not a proposed algebra. For a tail in the $i$th $A$ slot, put $D_i=\prod_{j\ne i}(a_i-a_j)^{-2}$ and $S_i=(z-a_i)^2D_i$. Its scalar derivative factor is $\sum_j(z-a_j)^{-1}-t/(z-a_i)$, so the same identity gives

\[
-\gamma\lambda_1\mathcal R_A S_i
\left[\sum_j\frac1{z-a_j}-\frac{t}{z-a_i}\right]
=\frac{\gamma\lambda_1}{\rho}\partial_z[\mathcal R_A S_i].
\tag{M20}
\]

The resolved source integral has $\pi^{-1}\operatorname{FP}\int G_A=S$, $S=\sum_iS_i$, on the real axis. Combining its action sign with all three tails gives

\[
\boxed{\mathcal K_{AAAB,\lambda_1}^{\rm scalar}
=\frac{\lambda_1}{\rho(1-t)}\partial_z[\mathcal R_A S].}
\tag{M21}
\]

The circle and height corrections also remain derivatives. If $R_i=\operatorname{Res}_{w=a_i}G_A$, $T_i$ is its regular Laurent coefficient, and $\tau_i=\operatorname{Im}a_i$, the finite-hole scalar primitive is

\[
\frac{\lambda_1}{\rho}\mathcal R_A
\left[(1+\gamma)S+2i\sum_i\tau_iR_i
-\epsilon^2\sum_{i\ \mathrm{inner}}T_i\right].
\tag{M22}
\]

Green's theorem gives this formula directly. The $B$-centered hole contributes zero to the scalar expression because $G_A(z)=0$ and $\partial_z(\mathcal R_AG_A)|_{w=z}=0$. Thus the moving-hole boundary is accounted for. The other height/hole terms have the same joint degree as $S$ and are retained until taking the initial limit (M13).

In each middle spatial chamber the $z$ endpoints are neighboring $A$ positions. The primitive in (M21) vanishes there for $\operatorname{Re}\rho<0$ at distinct majority points. Integration by parts against a nonconstant fused coefficient gives degree $1+3(g-\rho)$. The uncontracted $C(z)$ term and all higher Taylor levels also have at least one extra power. For an outer chamber the fixed spatial endpoint supplies the same improvement. The proper extensions preserve these identities from the common initial domain. Therefore the individual even-root chamber kernels have no pole along $g=\rho$.

### The new proper $AA$ subtraction is fixed independently

For two $A$ currents, the source kernel is simply $r^g/[(w-a_1)^2(w-a_2)^2]$. Direct residues give area integral $2\pi/r^2$. The two tails each have coefficient $-\gamma\lambda_1r^{g-2}$. Their sum with the action is

\[
-\frac{2\lambda_1}{1-t}r^{g-2}\mathcal W_2,
\qquad \mathcal W_2=:e^{-(2d+T)X-T\bar X}:.
\tag{M23}
\]

Thus in the working patch, where $0<g<1$, the needed ordered-pair bulk subtraction is

\[
CT_{AA,\lambda_1}
=-\frac{2\lambda_1}{(1-t)(g-1)}\epsilon^{g-1}
\int_L^R\mathcal W_2(y)\,dy.
\tag{M24}
\]

There is no odd midpoint jet, and the next endpoint power is $g>0$. In a wider continuation strip further endpoint/descendant terms are fixed by the same Taylor rule.

For three $A$ currents the leading combined coefficient is $-\lambda_1\Delta^g\sum_iD_i/(1-t)$. The new calculation verifies that its pair subface agrees with (M23) times the remaining $A$ vertex, with no odd midpoint pole. Its full radial degree $3g-2$ is positive in (M17); only its proper pair subfaces require extension there. The minority-coordinate derivative in (M21) commutes with this extension because the like-$A$ regulator and pair coordinates are independent of $z$. No finite quartic term was chosen to obtain (M21).

The proper $AAB$ scalar representative also retains its infinity term. With two majority poles, $\pi^{-1}\operatorname{FP}\int (z-w)^2/\prod_{i=1}^2(w-a_i)^2=S_2-1$. The source and two tails therefore combine into $\lambda_1\partial_z[\mathcal R_2(S_2+t-1)]/[\rho(1-t)]$, not an expression containing $S_2$ alone. Its positive diagonal degree $1-\rho$ in (M17) means it introduces no further marginal pole; its finite coefficient is nevertheless included in the common family.

## 6. Restrict the complete family, rather than individual divergent terms

In the patch (M17), the proper extensions needed for the four-fold product are now explicit: the ultraviolet current products, the first-Yukawa $BB$ subtraction with its remaining-source or tail dressing, the new even-root $AA$ subtraction, and the source/action/current matching in (M9). Triple products include their continued proper pair subfaces. The two marginal full-channel derivative extensions are fixed by the calculations above. All other new residual full clusters have positive degree.

Taylor extraction of an operator coefficient or endpoint jet is linear and acts on insertion coordinates. The permutation weights in (M12) are constant in those coordinates. Thus the same matrix identity holds for every properly subtracted diagram and for their complete sum. Matching the source-pair and current subtractions matters here: replacing $t/h$ by one before restriction, omitting the $AA$ subtraction, or choosing unrelated finite parts in different words would not be this common family.

The source-pair and current counterterms vanish as standalone cutoff graphs in the initial domain. They represent the required local subtractions on continuation, not a finite independent $P$ perturbation. In the physical patch their joint powers, including $\epsilon^{2h}$, are positive, so they cannot leave an extra finite quartic contact at a cutoff intersection. Dropping them before fixing the continuation would nevertheless lose the normalization checked in (M6)–(M10).

There is consequently a complex neighborhood of each generic point of (M17) in which the complete chamber families have regular restrictions to (M14). The only candidate divisors containing that surface were the inherited one-Yukawa marginal channel and the even-root marginal channel; both were removed explicitly. The other Mellin denominators are nonzero there. This rules out a pole-times-zero remnant from (M12).

Since $D(s,s)=0$, the full second-order word defects vanish on (M14), proving (M1) in this product prescription. The result then continues meromorphically in the coupling wherever the same restricted families are defined and regular. Exceptional intersections and logarithmic normalizations require a separate treatment.

This argument does not set all individual contact distributions to zero. Their values are fixed by the common extension, and their quartic sum vanishes. It also does not construct an associative multiplication of global charges from these renormalized four-fold contour products.

The finite-composite choice is substantive. An independent $\kappa_RP$ would have full $BBBA$ degree $3(g-\rho)-2$ without the compensating cutoff power. Its second jet is then marginal, outside the proof above; it also leaves the already identified local conservation obstruction. The conclusion here is for $\kappa_R=0$ with the matched quantum normalization, not for an arbitrarily modified action.

## 7. Evidence and remaining work

| New script | Output | PASS | Main checks |
|---|---|---:|---|
| [second_order_kernels.py](../referee_scratch/sl21_serre_m2_20261024/second_order_kernels.py) | [output](../referee_scratch/sl21_serre_m2_20261024/second_order_kernels.out) | 115 | All second-order Wick/cocycle channels in both quartics; exact first-order scalar derivative audit |
| [m2_normalization.py](../referee_scratch/sl21_serre_m2_20261024/m2_normalization.py) | [output](../referee_scratch/sl21_serre_m2_20261024/m2_normalization.out) | 25 | Full Ward response, independent $h$ subtractions, nonconstant radial and circle integrals, even-root derivative, $AA$ fusion and $AAB$ infinity term |
| [collision_audit.py](../referee_scratch/sl21_serre_m2_20261024/collision_audit.py) | [output](../referee_scratch/sl21_serre_m2_20261024/collision_audit.out) | 127 | Root products, all source-containing subsets and tail variants, exact positivity on a parameter patch, nested bound, exact shuffle matrix |

An initial normalization-script run is retained separately. Its single failed assertion used $\partial E/E$ in place of $\partial F_A/F_A=-\partial E/E$. Correcting that test sign reproduced the stated $H_A$ coefficient; no input formula was changed. The final scripts all finish successfully.

The immediate independent audit should scrutinize the three-parameter regulator and the completeness of the pole list, especially the transfer of the first-order derivative extension into a two-source graph and the new even-root extension. These are the logical steps supporting the formal order-two result; passing phase checks alone would not suffice.

The remaining physical question is still boundary closure. The mixed relation contains a fermionic endpoint operator whose values on massive spatial boundaries have not been computed. Neither those values, an action on soliton sectors, a coproduct, nor a faithful presentation follows from (M1). Higher-order contour relations also remain unproved. No affine $\mathfrak{sl}(2|1)$ or affine $\mathfrak{sl}_2$ symmetry is assumed or established.

## Integration

No book integration was performed. **Independently verified:** the local order-two Ward normalization, all separated second-order kernels, the first-order removable pole, and the new even-root derivative. **Extended:** the matched source-pair analytic regulator, the even-root pair subtraction, and the complete generic-coupling contour argument through mass order two. **Not established:** arbitrary-scheme or physical equal-time realization, exceptional couplings, all-order relations, global charge/boundary realization, or a symmetry algebra. Existing files remain unchanged, and no theorem upgrade is recommended.
