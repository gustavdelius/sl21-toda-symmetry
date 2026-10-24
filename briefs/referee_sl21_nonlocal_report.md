# Nonlocal-charge construction for the root-system sl(2|1) Toda action

Date: 2026-10-23. Model: Codex (GPT-6).

Book baseline: `745586e6a6f7bca5ed91a1e9cf89d6ef5df387f0`.

Task: construct conserved nonlocal charges for the action studied in the two preceding sl(2|1) reports, with the eventual aim of identifying its soliton S-matrix.

This is an independent research/referee report. No existing book file, report or book check script was edited. The new calculations are under `referee_scratch/sl21_nonlocal/`. No author script was imported. No novelty claim is made.

## Outcome

An explicit even component Lax connection has been constructed and checked against every component of the full nonlinear equations, including independent Grassmann fermions and arbitrary bosonic first derivatives. From it one obtains an exact nonlocal matrix current and conserved transport charges on a domain with appropriate endpoint conditions.

This does **not** complete the requested construction of the internal quantum-affine soliton charges. The odd transport components change under spatial translations. Also, vacuum subtraction magnifies the kink's exponentially small time connection into a nonzero endpoint flux. They cannot be presented as the desired conserved internal charges on the physical kink scattering space.

| Item | Verdict |
|---|---|
| Component Lax connection for the specified action | **Correct**, independently verified below. |
| Nonlocal transport current and its exact finite-interval boundary law | **Correct**, with the specified signs and ordering. |
| Conserved odd transport quantities with compactly supported deviations from a vacuum | **Correct classically**; they are spatially weighted quantities, not translation-invariant internal charges. |
| Finite nonzero odd transport on the exact static dressed kink | **Correct as a static evaluation**, not proof of conservation on a dynamical kink sector. |
| Identification of these quantities with quantum-affine soliton symmetry generators | **Wrong**. Translation and endpoint tests fail. |
| Desired action-derived quantum-affine internal charges | **Unclear; not constructed**. No no-go theorem for such charges is claimed. |

## 1. Action and algebra conventions

Use exactly the action, roots and source normalization specified in `briefs/referee_sl21_kink_modes_report.md`. The underlying source is Evans–Madsen, [hep-th/9703065](https://arxiv.org/abs/hep-th/9703065), corrected v2, equations (3.4)–(3.5). Its warning is relevant: “the Toda model really does depend on the particular simple root system”.

Write phi=(u,v), partial_±=partial_t ± partial_x, and

\[
\alpha_0=(-1/\sqrt2,i/\sqrt2),\quad
\alpha_1=(\sqrt2,0),\quad
\alpha_2=(-1/\sqrt2,-i/\sqrt2).
\]

The bosonic potential is V=mu^2(exp(sqrt(2)u)+exp(-sqrt(2)u)). The fermion mass matrix is M=sum_{j=0,2} mu exp(alpha_j.phi) alpha_j alpha_j^t. The source equations are

\[
\partial_+\partial_-\phi=-\nabla V
+i\sum_{j=0,2}\mu\alpha_j
(\alpha_j\cdot\psi_+)(\alpha_j\cdot\psi_-)
e^{\alpha_j\cdot\phi},
\]

\[
\partial_-\psi_+=-M\psi_-,\qquad
\partial_+\psi_-=M\psi_+.
\]

For the auxiliary triplet use parity (0,0,1), matrix units E_ij and the supertrace-zero Cartan embedding

\[
H(u,v)=\operatorname{diag}\left(
\frac{u+iv}{\sqrt2},\frac{-u+iv}{\sqrt2},\sqrt2iv
\right).
\]

In particular, H(alpha_0)=diag(-1,0,-1), H(alpha_1)=diag(1,-1,0), H(alpha_2)=diag(0,1,1). The step matrices are

\[
e_0=zE_{31},\ f_0=-z^{-1}E_{13},\quad
e_1=E_{12},\ f_1=E_{21},\quad e_2=E_{23},\ f_2=E_{32}.
\]

Here z is an auxiliary spectral parameter, not a soliton evaluation parameter already identified with rapidity.

Grassmann coefficients are written on the left. Products therefore use
(a X)(b Y)=(-1)^{|X||b|}ab XY. Omitting this sign reverses the fermion contribution to the curvature.

## 2. Explicit Lax connection

Define s_j=exp(alpha_j.phi/2), eta_{j±}=alpha_j.psi_±. Thus s_0 s_1 s_2=1. Set

\[
\begin{aligned}
L_+={}&\tfrac12 H(\partial_+\phi)
+i\mu(s_1E_{12}-z s_0s_2E_{21})\\
&-i\sqrt\mu\left(zs_0\eta_{0+}E_{31}
+s_2\eta_{2+}E_{23}\right),\\
L_-={}&-\tfrac12 H(\partial_-\phi)
+i\mu(s_1E_{21}-z^{-1}s_0s_2E_{12})\\
&+\sqrt\mu\left(-z^{-1}s_0\eta_{0-}E_{13}
+s_2\eta_{2-}E_{32}\right).
\end{aligned}
\]

Both connections are even in total parity. The flatness equation is

\[
\partial_+L_--\partial_-L_+ +[L_+,L_-]=0.
\]

The diagonal part gives the bosonic source equation after applying H. The four odd off-diagonal entries give the fermion equations projected along alpha_0 and alpha_2, which span the field space. The remaining entries cancel. This is a component zero-curvature representation of the full action, rather than a Lax pair for its bosonic truncation.

Independent computation: `solve_lax.py` starts with six undetermined coefficients in the even and odd step terms. It solves the off-diagonal constraints and then fixes the potential and Yukawa coefficients. `verify_lax.py` uses four independent Grassmann generators, arbitrary bosonic first jets, and the source equations to check all nine curvature entries. All checks pass.

This construction makes no identification between the auxiliary triplet and physical kink states. Ridout–Teschner explicitly distinguish their auxiliary quantum-algebra construction from soliton nonlocal charges in section 4 of [1102.5716](https://arxiv.org/abs/1102.5716); that distinction also applies to the inference one might try to make here.

## 3. Exact nonlocal current and boundary law

Put A_t=(L_++L_-)/2 and A_x=(L_+-L_-)/2. Then

\[
\partial_t A_x-\partial_x A_t+[A_t,A_x]=0.
\]

Define spatial transport by

\[
T(b,a;z)=\mathcal P\exp\left[-\int_a^b A_x(x;z)dx\right].
\]

Differentiating the transport equation and using flatness gives the exact boundary law

\[
\partial_t T(b,a;z)=-A_t(b;z)T(b,a;z)
+T(b,a;z)A_t(a;z).
\]

For a constant vacuum with u=0 and psi=0,

\[
A_t^{\rm vac}=
\frac{i\mu}{2}\left[(1-z^{-1})E_{12}+(1-z)E_{21}\right].
\]

It vanishes at z=1. At this point the spatial vacuum connection is

\[
B=i\mu(E_{12}-E_{21}),\qquad
\operatorname{spec}B=(\mu,-\mu,0).
\]

For deviations of compact support from this vacuum, introduce the time-independent gauge

\[
K_x=e^{xB}A_x e^{-xB}-B,\qquad
K_t=e^{xB}A_t e^{-xB}.
\]

Let W(x,a) solve W_x=-K_x W with W(a,a)=I. Since K_t(a)=0, flatness implies W_t=-K_t W. Thus the explicit nonlocal current

\[
\rho(x)=K_x(x)W(x,a),\qquad j(x)=-K_t(x)W(x,a)
\]

satisfies partial_t rho + partial_x j=0. Its integrated charge is

\[
Q=\int_a^b\rho\,dx=I-W(b,a).
\]

When both endpoints lie in the vacuum, its flux vanishes and Q is conserved. The entries with one even and one odd auxiliary index are Grassmann-odd conserved quantities. This is an exact classical construction on that domain. For weaker decay the weighted integrals and endpoint limits require separate checks.

## 4. Translation test: why this is not the soliton symmetry

Translate a compactly supported field configuration to the right by a distance c. Changing variables in the ordered exponential gives

\[
W_c=e^{cB}W e^{-cB},\qquad Q_c=e^{cB}Q e^{-cB}.
\]

Hence these charges generally change under spatial translations. In the basis diagonalizing B, their odd entries have spatial weights +mu or -mu. None is a translation-invariant internal charge. Directly solving [B,X]=0 for an arbitrary odd 3-by-3 auxiliary matrix gives X=0. This is checked exactly in `transport_checks.py`.

It follows that this construction cannot supply the internal U_q(sl(2|1)) generators needed for a relativistic soliton intertwining equation. A conserved quantity with an exponential spatial weight is not sufficient. Scalar spectral invariants and other even transport invariants do not resolve the missing odd generators.

This is a failure of this route to the desired charges. It is not a proof that suitable disorder-type quantum charges do not exist.

## 5. Exact kink evaluation and the endpoint trap

On the static kink y=exp(2 mu x), u=i(4 arctan y)/sqrt(2), the bosonic spatial connection at z=1 is

\[
A_x^{\rm body}=-B\tanh(2\mu x).
\]

Its fundamental solution is U_0(x)=exp[B log(cosh(2 mu x))/(2 mu)]. The full static two-generator solution from the preceding report has psi_+=f(a,b), psi_-=f(a,-b), f=sqrt(y)/(1+iy), and its nilpotent v correction. That correction does not enter A_x: its effect on odd terms is of Grassmann degree three, and the even terms are independent of v.

Define

\[
M_a=\begin{pmatrix}0&0&-1\\0&0&i\\i&1&0\end{pmatrix},\qquad
M_b=\begin{pmatrix}0&0&-i\\0&0&-1\\1&-i&0\end{pmatrix}.
\]

In the body transport frame, the odd insertion is (sqrt(mu)/2) f^2(a M_a+b M_b). Since integral f^2 dx=-i/(2 mu), the symmetric-cutoff vacuum-subtracted transport has a finite exact limit

\[
\mathcal T_K=I+D+\frac{ab}{64\mu}[M_a,M_b],\qquad
D=\frac{i}{4\sqrt{2\mu}}(aM_a+bM_b).
\]

The displayed commutator is an ordinary matrix commutator; the Grassmann/Koszul signs have already been included. All higher powers vanish. `transport_checks.py` checks the exact integrand and the nilpotent exponential. `kink_transport.py` independently integrates the first variations numerically at mu=sqrt(2). Their Frobenius norms are

| mu L | Each odd first-variation norm |
|---:|---:|
| 2 | 0.2971522076 |
| 4 | 0.2973017286 |
| 6 | 0.2973017787 |
| 8 | 0.2973017788 |

The analytic limit is 1/(2 sqrt(2 mu)). This shows that the transport can detect both zero modes. It does not prove a conserved internal symmetry.

In fact, there is a second obstacle. On the bosonic kink,

\[
A_t=\mu\operatorname{sech}(2\mu x)N,\qquad
N=\begin{pmatrix}i&-1&0\\-1&-i&0\\0&0&0\end{pmatrix}.
\]

Although the raw A_t tends to zero, [B,N]=-2 mu N. The left vacuum has B_L=B, the right B_R=-B, and

\[
\lim_{x\to\pm\infty}e^{xB_{L/R}}A_t e^{-xB_{L/R}}=2\mu N\ne0.
\]

Thus the vanishing-flux argument for compactly supported deviations cannot simply be extended to the kink by citing A_t -> 0. Even its bosonic endpoint transport evolves by conjugation. Fermion-dependent endpoints must be retained for the full dressed solution. A finite value on a static solution is not evidence for conservation on all time-dependent configurations in its sector.

## 6. What remains unconstructed

The desired charges must be finite on a specified kink scattering space, conserved including endpoint terms, and commute with both time and space translations. Their Lorentz spins, action on topological sectors and coproduct must follow from the action. Quantum conservation additionally requires a regulated operator construction and an anomaly check. None of these properties follows from the auxiliary zero-curvature representation alone.

The appropriate next research route is a disorder/vertex-current construction in a specified quantum theory, with conservation checked against all even-root, odd-root and composite bosonic perturbations. Merely importing screening operators from an auxiliary monodromy construction does not perform that check. No such verified full quantum construction has been obtained in this calculation.

## Calculation record and integration

New independent scratch scripts and matching output files:

- `referee_scratch/sl21_nonlocal/solve_lax.py`: coefficient solve.
- `referee_scratch/sl21_nonlocal/verify_lax.py`: full action-to-curvature check.
- `referee_scratch/sl21_nonlocal/transport_checks.py`: exact kink integral, Grassmann exponential, translation and endpoint tests.
- `referee_scratch/sl21_nonlocal/kink_transport.py`: independent numerical transport integration.

Run with `python3 -u`; dependencies are SymPy, NumPy and SciPy. Checks pass. No book integration was performed. Verified: the component connection, classical nonlocal current, boundary law, kink transport and stated obstructions. Rejected: interpreting the odd transport entries as internal soliton quantum-group generators. Unresolved: the charges requested for the S-matrix identification. Independent review is required before any new formula is promoted into book text.
