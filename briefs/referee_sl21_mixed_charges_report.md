# Root-system sl(2|1) Toda: antichiral currents and mixed-charge endpoints

Date: 2026-10-23. Model: Codex (GPT-6).

Commit reviewed: 745586e6a6f7bca5ed91a1e9cf89d6ef5df387f0.

Task: continue the Bernard–LeClair investigation after the local mass-squared Ward calculation. Construct antichiral partners, examine mixed products without prescribing an affine algebra, and make their endpoint subtractions explicit. No existing file was edited. No subagent, git mutation, commit, push, or book render was used.

Independent calculations are in [mixed_checks.py](../referee_scratch/sl21_mixed_charges_20261023/mixed_checks.py) and [its output](../referee_scratch/sl21_mixed_charges_20261023/mixed_checks.out). This script imports no preceding check script. The derivations below distinguish three operations: local radial OPEs, canonical equal-time exchange, and the still incomplete construction of globally regulated charge products.

## Outcome and verdicts

The antichiral currents follow from a parity automorphism of the fermionic action. The canonical disorder strings give an explicit separated-density braid matrix. A constructed zero-mode cocycle matches it to the local OPE signs on the additive charge cover. The leading Bernard–LeClair descent has two different endpoint fields: a bosonic exponential for one diagonal mixed product and a fermion bilinear for the other.

The fermionic endpoint requires an additional composite-operator subtraction. It is more singular than the current-insertion subtraction found previously. Write $t=b^2$, $k=d+T$, and

\[
\begin{split}
\mathscr U_B&=i:\chi_2\bar\chi_2 e^{(d+p_2)\cdot\Phi}:,\qquad
F_B=:e^{k\cdot\Phi}:,\\
C&=p_0\cdot\partial\varphi+j,\qquad
\bar C=p_0\cdot\bar\partial\bar\varphi+\bar j.
\end{split}
\]

For the local descent seed $\mathcal B_{BB}^{(1)}=t^2\lambda_2\mathscr U_B$, its minimally subtracted insertion is

\[
\boxed{
[\mathcal B_{BB}^{(1)}]_{\rm ren}
=\lim_{\epsilon\to0}\left[
t^2\lambda_2\mathscr U_B
-\frac{t^2\lambda_0\lambda_2}{t-1}\epsilon^{2t-2}F_B
-t\lambda_0\lambda_2\epsilon^{2t}:C\bar C F_B:
\right]_{S_\epsilon}+O(m^3).
}
\tag{1}
\]

This is a limit of insertions in correlation functions. It is **not yet an identification of the right-hand side of a completed mixed-charge commutator**. The compatible phase algebra constructed below does not evaluate the simultaneous two-charge cutoff limit. Equation (1) is an independently established subtraction for the field that occurs in the leading local descent.

| Item | Verdict | Evidence and scope |
|---|---|---|
| Bosonized Yukawa momenta and both candidate currents | **Correct locally** | Direct fermionic Wick functions, point splitting, and fresh OPE calculation |
| Antichiral partners in the same local convention | **Correct** | Explicit parity automorphism, including Grassmann signs |
| Canonical separated-density braid matrix | **Correct for the stated string definition** | Canonical commutators and BCH calculation |
| Matching canonical phases to the local cocycle on the additive cover | **Correct** | Explicit Weyl-vertex rephasing, equations (16a)–(16c) |
| Multiplying the radial cocycle onto those canonical phases as an additional independent Klein factor | **Wrong** | It changes the current/tail exchange sign; section 3 |
| Leading local mixed-product endpoint seeds | **Correct as local descent data** | The two factorizable derivative residues; section 4 |
| Endpoint subtraction (1) | **Correct through mass order two in the stated prescription** | Full fermionic collision and a nontrivial spectator integral |
| Replacing the fermionic endpoint by a bosonic Cartan exponential as a UV operator | **Wrong** | Fermion spectator distinguishes them |
| Identifying its renormalized value at spatial infinity with a Cartan factor | **Unclear** | Massive boundary matrix elements have not been calculated |
| Generic exclusion of additional homogeneous local endpoint terms | **Correct under the listed CPT assumptions** | Descendant weights and a complete integer argument |
| A closed mixed-charge algebra, coproduct on solitons, or either proposed affine algebra | **Unclear; not established** | Contour matching, boundary operators, and same-chirality relations remain |

No theorem upgrade is proposed.

## 1. Normalization and independent input checks

Use $z=x+i\tau$, ordinary area $d^2z=dx\,d\tau$, and $\bar\partial(1/z)=\pi\delta^{(2)}(z)$. The canonical ultraviolet fields obey

\[
\varphi^a(z)\varphi^b(w)\sim-\delta^{ab}\log(z-w),\qquad
\chi_0(z)\chi_2(w)\sim\frac1{z-w}.
\tag{2}
\]

Here $\Phi=\varphi+\bar\varphi$, $\chi_i=\alpha_i\cdot\psi$, and the bilinear roots are

\[
\alpha_0=(-1/\sqrt2,i/\sqrt2),\quad
\alpha_1=(\sqrt2,0),\quad
\alpha_2=(-1/\sqrt2,-i/\sqrt2).
\]

Left and right fermions anticommute. Put

\[
p_i=-ib\alpha_i,\quad t=b^2,\quad
T=p_0+p_2,\quad d=-T/t,\quad k=d+T=(1-t)d.
\tag{3}
\]

Direct inner products give $p_i^2=0$, $p_0\cdot p_2=-t$, $d\cdot p_i=1$, and $d^2=-2/t$. The action convention is

\[
\begin{split}
S_0&=\frac1{2\pi}\int d^2z\,
(\partial\Phi\cdot\bar\partial\Phi+
\psi\cdot\bar\partial\psi+\bar\psi\cdot\partial\bar\psi),\\
S_{\rm int}&=\frac1\pi\int d^2z\,
(\lambda_0Y_0+\lambda_2Y_2+\lambda_1V_1+\kappa(\epsilon)P),\\
Y_i&=i:\chi_i\bar\chi_i e^{p_i\cdot\Phi}:,\quad
V_1=:e^{-T\cdot\Phi}:,\quad P=:e^{T\cdot\Phi}:,\\
\kappa(\epsilon)&=-\frac{\lambda_0\lambda_2}{t}\epsilon^{2t}.
\end{split}
\tag{4}
\]

The two Yukawa couplings have mass order one, and $\lambda_1$ and $\kappa$ have mass order two. These mass orders are not their anomalous scaling dimensions. The quantum prescription in (4) is circular exclusion and minimal analytic continuation, with no finite independent $P$ coupling. For the component-action convention used in the preceding reports, $b=ig/\sqrt{4\pi}$ and $\lambda_i=\mu_i/2$. Thus real $g$ gives negative $t$; formulas obtained initially in a convergent positive-$t$ region must be continued with their counterterms retained.

The fermionic sign in the mixed Yukawa OPE is fixed directly:

\[
(i\chi_0\bar\chi_0)(i\chi_2\bar\chi_2)
\sim\frac{+1}{(z-w)(\bar z-\bar w)},\qquad
Y_0(z)Y_2(w)\sim |z-w|^{2t-2}P(w).
\tag{5}
\]

The Wick crossing and $i^2$ supply the two minus signs. A free boson $h$ can therefore represent the chiral fermions by momenta $+i,-i$, with $j=:\chi_0\chi_2:=i\partial h$. Point splitting independently gives the same $j$ Ward action. For example, with $a,b$ away from zero,

\[
\frac{\langle j(0)\chi_0(a)\chi_2(b)\rangle}
{\langle\chi_0(a)\chi_2(b)\rangle}
=\frac1b-\frac1a.
\tag{6}
\]

The script also compares the complete chiral four-fermion Wick function with its exponential-bosonization expression. The three-dimensional Yukawa momenta are consequently

\[
a_0=(ib/\sqrt2,b/\sqrt2,i),\qquad
a_2=(ib/\sqrt2,-b/\sqrt2,-i).
\]

This verifies the local bosonization used to discover the currents. It does not establish the global bosonization lattice.

To avoid prejudging a root labeling, call the currents $A,B$:

\[
J_A=:e^{-d\cdot\varphi}:,\qquad
J_B=:(p_0\cdot\partial\varphi+j)e^{d\cdot\varphi}:.
\tag{7}
\]

Their Gaussian charge spins are $s_A=1/t-1$ and $s_B=1/t$. In the bosonized calculation $J_B=:(a_0\cdot\partial X)e^{d\cdot X}:$. Its ordered double-pole coefficients with $Y_0,Y_2$ are $1,t-1$; the simple-pole coefficients are respectively $+\partial K_0,-\partial K_2$, before the common cocycle sign. The distributional double-pole derivative is essential: the conservation coefficients are $1-1=0$ and $(t-1)-(-1)=t$.

Use the local cocycle of the current-renormalization report, $\varepsilon(d,p_i)=-1$, $\varepsilon(p_i,d)=+1$. The first-order components are

\[
H_A=-\gamma\lambda_1:e^{-k\cdot\varphi-T\cdot\bar\varphi}:,
\quad \gamma=\frac{t}{1-t},\qquad
H_B^{(1)}=-t\lambda_2K_2,
\tag{8}
\]

where $K_2=i:\chi_2\bar\chi_2e^{(d+p_2)\cdot\varphi+p_2\cdot\bar\varphi}:$.

The mass-squared cancellation can be checked again without importing its script. Set $x_i=p_i\cdot\partial\varphi$, $y_i=p_i\cdot\partial^2\varphi$, $E=e^{k\cdot\varphi+T\cdot\bar\varphi}$ and $L=:(x_0+j)E:$. The composite residue and the regulated Yukawa boundary polynomial are

\[
\begin{split}
R_P/E&=(x_0+j)(x_0+x_2)
-\frac t2[y_0+y_2+(x_0+x_2)^2],\\
N/E&=\tfrac12(y_0-y_2)+\partial j
+\tfrac12(x_0^2-x_2^2)+j(x_0+x_2).
\end{split}
\]

Since $k\cdot\partial\varphi=(1-1/t)(x_0+x_2)$, direct differentiation gives $N-R_P/t=\partial L$. With (4), this is the local identity behind

\[
H_B^{\rm ren}=\lim_{\epsilon\to0}
[-t\lambda_2K_2+\lambda_0\lambda_2\epsilon^{2t}L]_{S_\epsilon},
\qquad J_B^{\rm ren}=\operatorname{FP}J_B.
\tag{9}
\]

This fresh algebraic check does not replace the full distributional calculation in the preceding report. The new endpoint calculation below uses the same action and subtraction scheme.

## 2. Antichiral partners

Apply the transformation

\[
z\leftrightarrow\bar z,\quad
\varphi\leftrightarrow\bar\varphi,\quad
\psi\mapsto\bar\psi,\quad\bar\psi\mapsto-\psi.
\tag{10}
\]

Both kinetic terms are exchanged. The Yukawa factor transforms as $i\chi_i\bar\chi_i\mapsto-i\bar\chi_i\chi_i=i\chi_i\bar\chi_i$. Also $j\mapsto\bar j=:\bar\chi_0\bar\chi_2:$. Thus the transformation preserves the action and circular regulator. Its square acts as fermion parity on elementary fermions and as the identity on these even currents.

The antichiral partners are

\[
\begin{split}
\bar J_A&=:e^{-d\cdot\bar\varphi}:,&
\bar H_A&=-\gamma\lambda_1:e^{-T\cdot\varphi-k\cdot\bar\varphi}:,\\
\bar J_B&=:(p_0\cdot\bar\partial\bar\varphi+\bar j)e^{d\cdot\bar\varphi}:,&
\bar H_B^{\rm ren}&=\lim_{\epsilon\to0}
[-t\lambda_2\bar K_2+\lambda_0\lambda_2\epsilon^{2t}\bar L]_{S_\epsilon},
\end{split}
\tag{11}
\]

with

\[
\bar K_2=i:\chi_2\bar\chi_2e^{p_2\cdot\varphi+(d+p_2)\cdot\bar\varphi}:,
\quad
\bar L=:\bar C e^{T\cdot\varphi+k\cdot\bar\varphi}:.
\]

Bars here label chirality, not Hermitian conjugation. They obey the parity transform of the local Ward identity, $\partial\bar J=\bar\partial\bar H$, to the same mass order and under the same conditions.

The parity-symmetric local sign part on the charge cover $(d_L,d_R,p_0,p_2)$ is

\[
\varepsilon(n,m)=(-1)^{(n_L+n_R)(m_0+m_2)}.
\tag{12}
\]

This supplies local OPE signs. The cross-chiral zero-mode factor that must accompany it is constructed in the next section.

## 3. Canonical exchange, and why it does not yet give a coproduct

Take the standard left-anchored split at equal Minkowski time,

\[
\varphi(x)=\tfrac12\left(\Phi(x)+\int_{-\infty}^x\dot\Phi(y)dy\right),\quad
\bar\varphi(x)=\tfrac12\left(\Phi(x)-\int_{-\infty}^x\dot\Phi(y)dy\right),
\tag{13}
\]

with $[\Phi^a(x),\dot\Phi^b(y)]=4\pi i\delta^{ab}\delta(x-y)$. Integrating the delta functions, away from coincidence, gives

\[
[\varphi(x),\varphi(y)]=-\pi i\operatorname{sgn}(x-y),\quad
[\bar\varphi(x),\bar\varphi(y)]=+\pi i\operatorname{sgn}(x-y),\quad
[\varphi(x),\bar\varphi(y)]=-\pi i.
\tag{14}
\]

The identity matrix on boson components is implicit. With $c_A=-d,c_B=d$, BCH therefore gives

\[
J_a(x)\bar J_b(y)=R_{a\bar b}\bar J_b(y)J_a(x),\qquad
R=\begin{pmatrix}\zeta&\zeta^{-1}\\\zeta^{-1}&\zeta\end{pmatrix},
\quad\zeta=e^{2\pi i/t},\quad x\ne y,
\tag{15}
\]

for the canonical realization (13), with no additional independent Klein twist. Derivative prefactors and even fermion bilinears contribute coincident contacts, not separated phases.

For a general vertex with left and right momenta $(l,r)$, the exponential exchange with $(l',r')$ is

\[
\begin{cases}
\exp[-i\pi(l-r)\cdot(l'+r')],&x>y,\\
\exp[+i\pi(l+r)\cdot(l'-r')],&x<y.
\end{cases}
\tag{16}
\]

Replacing a current by a tail adds the same local momentum to its two chiralities. The change in (16) is an integer multiple of $2\pi i$, since $d\cdot p_i=1$ and $d\cdot T=2$. Thus the canonical exponential calculation gives the same separated braiding for $J$ and $H$, including the local order-two tail.

There is a convention trap. If (12) is now multiplied onto the already canonical exchange as a *new independent* Klein factor, $J_B$ and $\bar J_B$ get sign $+1$, whereas $J_B$ and $\bar K_2$ get sign $-1$. Indeed their charge-cover labels are $(1,0,0,0)$ and respectively $(0,1,0,0)$ or $(0,1,0,1)$. This violates the common-current/tail braid assumption. The local cocycle is not inconsistent; its relation to the normal-ordering phases must be included exactly once.

Here is an explicit compatible phase construction. For a lattice label $n$, put $N_n=n_0+n_2$, $q_n=n_0p_0+n_2p_2$, $l_n=n_Ld+q_n$, and $r_n=n_Rd+q_n$. Symmetrically ordered Weyl vertices $W_n=\exp(l_n\cdot\varphi+r_n\cdot\bar\varphi)$ have the cross-chiral normal-ordering cocycle

\[
\varepsilon_W(n,m)=
\exp\left[\frac{\pi i}{2}
(r_n\cdot l_m-l_n\cdot r_m)\right].
\tag{16a}
\]

This follows by moving $e^{r_n\cdot\bar\varphi}$ past $e^{l_m\cdot\varphi}$ using (14), including the half-commutator factors for symmetric ordering. Rephase the basis by

\[
V_n=f(n)W_n,\qquad
f(n)=\exp\left[\frac{\pi i}{2}(n_L-n_R)N_n\right].
\tag{16b}
\]

Then $\varepsilon_V(n,m)=f(n)f(m)f(n+m)^{-1}\varepsilon_W(n,m)$ is

\[
\boxed{
\varepsilon_V(n,m)=
(-1)^{(n_L+n_R)N_m}
\exp\left[\frac{\pi i}{t}(n_Lm_R-n_Rm_L)\right].
}
\tag{16c}
\]

To check this, the terms involving $N_n,N_m$ in the exponent reduce to $-\pi i(n_L-n_R)N_m$, which equals the sign in (16c) for integer labels. The remaining term uses $d^2=-2/t$. Equation (16c) is bimultiplicative and associative. It restricts to the specified current/Yukawa cocycle in either chirality. Between $d_L$ and $d_R$, its exchange ratio is $e^{2\pi i/t}$.

For arbitrary labels on this cover, the chiral branch phase for $x>y$ is $\exp[-\pi i(l_n\cdot l_m-r_n\cdot r_m)]$. Multiplying it by $\varepsilon_V(n,m)/\varepsilon_V(m,n)$ gives exactly (16). Thus the $J_B/\bar K_2$ minus sign just discussed is canceled by its chiral branch sign. The script checks the cochain identity symbolically and the complete current/tail exchange table. This constructs compatible phase data on the generic additive cover; it does not construct the physical compactification lattice, operator domains, or regulated products of integrated charges.

The original Bernard–LeClair PDF, printed page 106, explicitly says: “We suppose that the same braiding relations hold”, before replacing a current component by its tail. Their mixed-product calculation is qualified as “to lowest non-trivial order in perturbation theory”. Equations (2.15)–(2.18) were inspected in the supplied PDF. [Bernard–LeClair, published source](https://doi.org/10.1007/BF02099173).

For reference, in the canonical realization the formal topological operator

\[
H_{\rm top}=\frac b{2\pi}\alpha_1\cdot[\Phi(+\infty)-\Phi(-\infty)]
\tag{17}
\]

has weights $(+2,-2,-2,+2)$ on $(Q_A,Q_B,\bar Q_A,\bar Q_B)$. This follows directly from $[\Delta\Phi,\varphi]=-2\pi i$ and $[\Delta\Phi,\bar\varphi]=+2\pi i$. It does not supply their action on a Hilbert space.

One could write $q^2=\zeta^{-1}$ as notation. Neither the sign of $q$, its physical coupling dictionary, nor a quantum-group presentation follows from this choice. In particular, (16) depends on the electric and magnetic labels of the fields on which a charge acts. Reducing that exchange character to a specified Cartan operator requires the spectrum and its zero modes. A coproduct on soliton multiplets has not been calculated here.

## 4. Leading local mixed-product descent

Factor a Yukawa as $Y_i=V_i\bar V_i$, where $V_i=\chi_i e^{p_i\cdot\varphi}$ and $\bar V_i=i\bar\chi_i e^{p_i\cdot\bar\varphi}$. The derivative residues that enter the Bernard–LeClair descent are

\[
\begin{array}{c|cc}
&h&\bar h\\ \hline
A,\ V_1&-\gamma e^{-k\cdot\varphi}&-\gamma e^{-k\cdot\bar\varphi}\\
B,\ Y_2&-t\chi_2e^{(d+p_2)\cdot\varphi}&-it\bar\chi_2e^{(d+p_2)\cdot\bar\varphi}.
\end{array}
\tag{18}
\]

The $B,Y_0$ conservation residue is zero. The $A$ current has no Yukawa tail, and $B$ has no even-root tail at this order. Multiplying the two derivative primitives belonging to a common perturbation gives the local zero-forms

\[
\mathcal B_{AA}^{(2)}=\gamma^2\lambda_1:e^{-k\cdot\Phi}:,
\qquad
\mathcal B_{BB}^{(1)}=t^2\lambda_2\mathscr U_B,
\qquad
\mathcal B_{AB}=\mathcal B_{BA}=0
\quad\hbox{at first descent order}.
\tag{19}
\]

The superscripts count mass order. The vanishing off-diagonal seeds mean that no perturbation has both needed first-order primitives. They do not, by themselves, calculate a second-order double-charge collision.

With a common contour ordering and braiding prescription, Bernard–LeClair's convention places $(2\pi i)^{-1}\int d\mathcal B_{a\bar b}$ on the right of the leading mixed braided commutator. Equation (16c) supplies compatible separated phase data for that construction. However, (19) fixes fused local descent seeds, including their phase convention; it does not evaluate the overlapping current and interaction collision regions at mass order two. Equations (15) and (19) should **not** be promoted to an established global charge algebra on this basis.

The $BB$ field is not a purely bosonic exponential: it has a nonzero free correlator with a $\chi_0\bar\chi_0$ spectator, while a purely bosonic insertion has zero such correlator. That is an operator-level distinction. It would also be wrong to set it to zero at spatial infinity merely by setting classical fermions to zero. Renormalized massive vacuum matrix elements, including possible condensates and mixing, are a different calculation. The current representative (7) is kept fixed throughout; improvements by total derivatives can affect nonlocal endpoint terms and have not been used to simplify (19).

## 5. Renormalization of the fermionic endpoint

Put the endpoint at $z$ and the other Yukawa at $z+r$. Direct Wick expansion gives

\[
\begin{split}
\chi_2(z)\chi_0(z+r)&=-r^{-1}-j(z)+O(r),\\
e^{(d+p_2)\cdot\varphi(z)}e^{p_0\cdot\varphi(z+r)}
&=r^{t-1}:e^{k\cdot\varphi(z)}:
[1+r\,p_0\cdot\partial\varphi+O(r^2)].
\end{split}
\tag{20}
\]

The two full Yukawa factors again have positive leading Wick sign. The local cocycle (12) contributes $+1$ in this full endpoint collision, because the endpoint contains both $d_L$ and $d_R$. Thus

\[
\mathscr U_B(z)Y_0(z+r)
=|r|^{2t-4}:
[1+rC+\bar r\bar C+r\bar r C\bar C+\cdots]F_B:.
\tag{21}
\]

In a circular annulus, pure $r^2$ and $\bar r^2$ terms integrate to zero. The two relevant radial integrals are

\[
\begin{split}
\int_{\epsilon<|r|<R}d^2r\,|r|^{2t-4}
&=\frac\pi{t-1}(R^{2t-2}-\epsilon^{2t-2}),\\
\int_{\epsilon<|r|<R}d^2r\,|r|^{2t-2}
&=\frac\pi t(R^{2t}-\epsilon^{2t}).
\end{split}
\tag{22}
\]

Multiplication by the insertion coefficient $t^2\lambda_2$ and the action expansion $-\lambda_0/\pi$ produces the divergent terms

\[
+\frac{t^2\lambda_0\lambda_2}{t-1}\epsilon^{2t-2}F_B
+t\lambda_0\lambda_2\epsilon^{2t}:C\bar C F_B:.
\tag{23}
\]

The negatives of these are exactly the counterterms in (1). No action-$P$ insertion is needed in this particular mixing calculation at mass order two: multiplying that action term by the order-one seed starts at order three. This does not remove the action counterterm from the separate two-current calculation.

For generic $t$ with $\operatorname{Re}t>-1$, $t\ne0,1$, the omitted rotationally invariant terms give an $O(\epsilon^{2\operatorname{Re}t+2})$ remainder. Minimal analytic continuation may be defined from $\operatorname{Re}t>1$, where the endpoint collision is integrable, with zero finite additions. At $t=1$ and $t=0$ the displayed denominators mark logarithmic limits and require separate prescriptions. For more negative $t$, more descendant subtractions are needed. None of those cases is covered by (1) as a two-term cutoff limit.

There is already a new endpoint divergence for $0<t<1$, where the order-two current subtraction proportional to $\epsilon^{2t}$ tends to zero. Therefore finiteness of the individual current does not establish finiteness of its mixed-product endpoint.

### A fermionic spectator check

The subtraction was tested beyond a scalar leading contraction. With a $\chi_0(a)\chi_2(b)$ spectator, the chiral four-point Wick function divided by its leading collision is

\[
G(r)=
\frac{\langle\chi_2(0)\chi_0(r)\chi_0(a)\chi_2(b)\rangle}
{(-1/r)\langle\chi_0(a)\chi_2(b)\rangle}
=1-\frac{r(a-b)}{a(r-b)}.
\tag{24}
\]

Choose $a=2,b=3$, and the corresponding right-moving spectator. The exact angular average is

\[
\langle G(r)G(\bar r)\rangle_{\arg r}
=1+\frac{|r|^2}{36(1-|r|^2/9)}.
\tag{25}
\]

The script integrates this rational function over the annulus; it does not replace it by a truncated OPE. After the two subtractions, the difference from its independently summed finite-part series has leading term

\[
\frac{t^2}{324(t+1)}\epsilon^{2t+2}.
\tag{26}
\]

For $t=\sqrt2/3$ and $t=-\sqrt2/5$, cutoffs $0.01$ and $0.005$ verify both this coefficient and the predicted cutoff scaling. The descendant term is nonzero in this test, so the check detects omission of the second subtraction for negative $t$.

## 6. What homogeneous endpoint power counting does establish

Assume mass-analytic Gaussian conformal perturbation theory, ordinary local endpoint fields and descendants, and homogeneous finite subtractions with no added external scale. This excludes nonperturbative powers of masses, nonlocal endpoint operators, and unspecified boundary sectors. Let $n_0,n_2,n_e$ count the two Yukawas and the even-root potential, and define

\[
a=n_0-n_e,\quad c=n_2-n_e,\quad S=a+c,\quad
F=n_0-n_2,\quad M=n_0+n_2+2n_e.
\]

The bosonic interaction momentum is $ap_0+cp_2$. A chiral fermion number $F$ has minimum weight $F^2/2$: filling distinct fermionic modes gives $1/2+3/2+\cdots+(|F|-1/2)$. This argument does not require bosonization.

Write $(\sigma_A,n_A)=(-1,0)$ and $(\sigma_B,n_B)=(+1,1)$, where $n_a$ is the derivative degree of the current. The endpoint must have total weights $(s_a,s_b)$, including the couplings. The Yukawa coupling weight is $1/2$ in each chirality, and the even-root coupling weight is $1-t$. Subtracting the free momentum and fermion weights gives

\[
\ell=n_a-1-M/2+\sigma_aS-F^2/2+t(n_e-ac),\qquad
\bar\ell=n_b-1-M/2+\sigma_bS-F^2/2+t(n_e-ac).
\tag{27}
\]

For real irrational $t$, integrality of both levels forces $n_e=ac$. Then $M/2+F^2/2=S(S+1)/2$, so the possible chiral levels are

\[
\ell_A=-\frac{(S+1)(S+2)}2,\qquad
\ell_B=-\frac{S(S-1)}2.
\tag{28}
\]

For $AA$, nonnegativity requires $S=-2,-1$. Since $ac=n_e\ge0$, the only admissible nonnegative insertion counts are $(n_0,n_2,n_e)=(0,0,1)$. The $S=-1$ pairs would make one original count negative. For $BB$, $S=0,1$ allows only the vacuum term and one Yukawa. The vacuum coefficient is zero in the free mixed product, and the $Y_0$ coefficient is zero by its derivative residue. The surviving seed is $Y_2$. For $AB$ and $BA$, the two allowed sets of $S$ are disjoint.

This is the complete integer argument under the assumptions above. Enumeration through twelve insertions of each type is an additional script check. It excludes additional *finite homogeneous local* endpoint structures at generic irrational $t$; it does not exclude the explicitly cutoff-dependent terms (23). Nor does it supply the missing global contour construction or determine a massive boundary expectation value. Rational resonances, including the previously identified $t=-1/2$ current anomaly channel, require separate calculations.

## 7. A useful local literature identification

There is a further independent check of the bosonization. Set $X=(\Phi^u,\Phi^v,H)$, where $H=h+\bar h$ bosonizes the even fermion bilinears, and define, for $t\ne2$,

\[
\eta=\sqrt{1-t/2},\quad
\zeta_0=\frac{-ib\Phi^v/\sqrt2+H}{\eta},\quad
\rho_0=\frac{\Phi^v+ibH/\sqrt2}{\eta}.
\tag{29}
\]

The transformation is complex orthogonal: $\zeta_0^2=\rho_0^2=1$ and $\zeta_0\cdot\rho_0=0$ as coefficient vectors. Rescale

\[
\varphi_F=\sqrt2\Phi^u,\quad\phi_F=\sqrt2\zeta_0,\quad
\rho_F=\sqrt2\rho_0,\quad
B_F=ib/2,\quad A_F=\eta/\sqrt2.
\]

Then $A_F^2-B_F^2=1/2$, the kinetic coefficient is $1/(16\pi)$, and the Yukawa exponents become $B_F\varphi_F\pm iA_F\phi_F$. The even and composite potentials become $e^{-2B_F\varphi_F}$ and $e^{2B_F\varphi_F}$. The $\rho_F$ direction is absent from these local interactions. A constant shift of $\phi_F$ balances nonzero $\lambda_0,\lambda_2$.

The resulting interacting flat-space action has the form of Fateev's first sine-Liouville deformation, equation (4.1), with $\mu=\sqrt{\lambda_0\lambda_2}/\pi$ and $\mu_1=\lambda_1/\pi$, together with the composite counterterm and the locally free direction. In his fermionic rewriting, equation (4.2), the added exponential “plays the role of the usual counterterm canceling the divergencies coming from the fermion loops”. These displayed formulas and the quoted sentence were read in the primary source. [Fateev, SIGMA 13 (2017), 080, pp. 4 and 6](https://arxiv.org/pdf/1705.06424).

Our independently transformed currents are

\[
J_A=e^{-\varphi_{F,L}/(2B_F)},\qquad
J_B=(B_F\partial\varphi_{F,L}+iA_F\partial\phi_{F,L})
e^{\varphi_{F,L}/(2B_F)}.
\tag{30}
\]

This is a local action comparison, not an identification of global theories or quantum charge algebras. In particular, the rotated fermion compactification lattice need not factor into independent $\phi_F$ and $\rho_F$ sectors. No background-charge stress tensor, scattering spectrum, or physical S-matrix is imported from that comparison, and it does not independently fix our quantum counterterm coefficient.

## 8. Checks, interpretation, and next calculation

Run

~~~text
python3 -u referee_scratch/sl21_mixed_charges_20261023/mixed_checks.py
~~~

All **59 checks pass**. They cover the local fermion/boson dictionary, current poles, the mass-squared polynomial identity, parity, canonical commutators and braid exponents, the cocycle double-counting counterexample and its explicit phase resolution, endpoint mixing with fermionic spectators, endpoint power counting, and the local sine-Liouville change of variables. They do not evaluate a full regulated product of two charges.

The last paragraph of `briefs/referee_sl21_current_renorm_report.md`, line 570, should be refined to avoid implying that local antichiral residues alone determine the algebra. Suggested replacement:

> Construct the antichiral partners and match the radial current cocycles to one explicit global string and zero-mode prescription. Then evaluate the mixed charge contours with simultaneous current and interaction cutoffs, retaining the renormalized endpoint fields. The coproduct additionally requires the braid character of the fields on which the charges act. Neither the affine algebra nor the endpoint boundary values should be prescribed in advance.

This is a refinement of the proposed work sequence, not a claim that the preceding local Ward calculation was wrong. No change to the earlier file was made.

The next decisive calculation is a **single common-regulator definition of $Q_B\bar Q_B$**. It must combine both order-two current tails, the action counterterm, the simultaneous double-Yukawa region, and the endpoint mixing (1), using the compatible phase convention (16c). The order of the interaction cutoff, charge-contour separation, and spatial endpoint limits must be stated. Only after that calculation should a renormalized boundary operator be replaced by a function of topological charge. Same-chirality products and Serre-type relations are separate remaining tests.

## Integration

No book integration was performed. **Verified independently:** the local bosonization, candidate-current residues, the mass-squared Ward polynomial, and the canonical exchange calculation. **Extended:** antichiral partners, an explicit compatible zero-mode cocycle, leading local descent seeds, the endpoint subtraction with a fermionic spectator check, and a conditional endpoint selection rule. **Modified interpretation:** individual current finiteness is insufficient for mixed-product finiteness; radial cocycles cannot simply be added to canonical phases as independent factors. **Rejected as unsupported:** a purely bosonic replacement for the endpoint as a UV operator, a global coproduct derived from the two-current Gram matrix, and either affine-algebra identification. **Still unresolved:** the fully regulated mixed-charge contour product, massive boundary matrix elements, exceptional couplings, and same-chirality algebraic relations.

Input hashes were checked after the work. AGENTS.md and the preceding reports recorded in the scratch manifest are unchanged. The tracked working tree remains unchanged.
