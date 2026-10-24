# Dual screenings for the sl(2|1) Toda charge problem

Date: 2026-10-23. Model: Codex (GPT-6).

Book baseline: `745586e6a6f7bca5ed91a1e9cf89d6ef5df387f0`.

Task: continue the Bernard–LeClair construction and examine the editor's suggestion that the symmetry should be identified from dual screenings rather than assumed to be the original superalgebra.

This is a new research/referee report. No existing file or book check script was edited. Independent calculations are in `referee_scratch/sl21_dual_screening/check_dual.py`, with output `check_dual.out`. All symbolic checks pass. No claim of novelty is made.

## Outcome

The dual-screening viewpoint yields a new, explicit derivative-dressed even current. It passes the two Yukawa OPE tests and has a regular OPE with the even-root perturbation. Its exponential momentum is opposite to that of the earlier even-node Bernard–LeClair current.

Both currents are even. Their exponential momentum Gram matrix is proportional to the affine A1 Cartan matrix. This is evidence for investigating a bosonic quantum-affine symmetry; it is not a derivation of its algebra, a Langlands-dual identification, or a physical S-matrix.

There is a specific remaining conservation test. The composite bosonic potential produces a non-derivative residue for the new current. This occurs at the same mass order as two Yukawa insertions. Their combined Ward identity and the quantum normalization of the composite operator have not been calculated. Full conservation therefore remains unclear.

## 1. Source and meaning of duality

Litvinov and Spodyneiko, *On W algebras commuting with a set of screenings*, [arXiv:1609.06271](https://arxiv.org/abs/1609.06271), JHEP11(2016)138, section 4, equations (4.14)–(4.15), gives a useful dual-screening construction. The source states that a dual screening can map to a “screening operator dressed by a differential polynomial of degree 1”.

For interacting fermionic screening momenta a_0,a_2 with a_0^2=a_2^2=-1, the dressed screening has momentum

\[
d=\frac{2(a_0+a_2)}{(a_0+a_2)^2}
\]

and a prefactor proportional to a_0.partial X. The source's statement concerns W-algebra screening commutants. Below the Bernard–LeClair residues are checked independently; conservation in the massive theory is not inferred merely from sharing a commutant.

For isotropic superalgebra roots, the ordinary formula alpha^vee=2 alpha/alpha^2 is undefined. Thus root-coroot duality cannot be applied node by node to our original root system. The correct charge algebra must be read from actual conserved currents and their braiding and Serre relations. It must not be fixed in advance by the algebra naming the classical Lax pair.

## 2. Bosonized ultraviolet screening system

Keep the canonical normalization of the preceding Bernard–LeClair report: phi^a(z) phi^b(w) ~ -delta^{ab} log(z-w), and psi^a(z) psi^b(w) ~ delta^{ab}/(z-w). Use b for its formal coupling, not the book's beta. The action-to-b map is still unresolved.

Let chi_0=alpha_0.psi and chi_2=alpha_2.psi. Their OPE is chi_0(w) chi_2(z) ~ 1/(w-z), and their self-contractions vanish. Bosonize them, with suitable cocycles, by chi_0=e^{ih}, chi_2=e^{-ih}, where h h ~ -log. Only currents with zero net h momentum will be used below; they are fermion-parity even.

For X=(phi^u,phi^v,h), the Yukawa chiral factors become e^{a_0.X}, e^{a_2.X}, where

\[
a_0=\left(\frac{ib}{\sqrt2},\frac b{\sqrt2},i\right),\qquad
a_2=\left(\frac{ib}{\sqrt2},-\frac b{\sqrt2},-i\right).
\]

The even-root perturbation has momentum a_1=-(a_0+a_2)=(-i sqrt(2)b,0,0). Put T=a_0+a_2. The component action's composite potential has momentum T. Direct calculation gives

\[
a_0^2=a_2^2=-1,\quad
a_0\cdot a_2=1-b^2,\quad T^2=-2b^2.
\]

This calculation shows why dualizing the original isotropic roots individually misses information: bosonization adds the fermion momentum, and dualizing their pair becomes well defined. Work at generic b, away from singular parameter values.

The vector (0,sqrt(2),ib) is orthogonal to all three perturbing momenta, with norm 2-b^2. This is a free direction of the formal bosonized perturbation system. No assertion about its globally defined quantum compactification or scattering Hilbert space is made.

## 3. Explicit new current and residues

The pair-dual momentum is

\[
d=\frac{2T}{T^2}=\left(-\frac{i\sqrt2}{b},0,0\right).
\]

The currents to investigate are

\[
J_+=:e^{-d\cdot X}:=:e^{i\alpha_1\cdot\phi/b}:,
\]

\[
\boxed{J_- = :(a_0\cdot\partial X)e^{d\cdot X}:
=:\left(-ib\alpha_0\cdot\partial\phi+i\partial h\right)
e^{-i\alpha_1\cdot\phi/b}: .}
\]

Overall normalizations remain free. The current J_- is not an odd-root current with a single fermion attached. It is even and has a derivative prefactor. Replacing a_0 by a_2 gives an equivalent integrated screening, up to sign and a total derivative.

Here is the independent residue calculation. If d.a=1, Wick expansion of e^{a.X(w)} against (C.partial X)e^{d.X(z)} gives

\[
\operatorname{Res}_{w=z}
e^{a\cdot X(w)}(C\cdot\partial X)e^{d\cdot X(z)}
=:[C+(C\cdot a)a]\cdot\partial X\ e^{(d+a)\cdot X}:.
\]

Take C=a_0. For a=a_0, the vector coefficient vanishes since a_0^2=-1. For a=a_2,

\[
a_0+(1-b^2)a_2=-b^2(d+a_2),
\]

so the residue is exactly

\[
-b^2\partial:e^{(d+a_2)\cdot X}:.
\]

Consequently the two Yukawa perturbations give a Bernard–LeClair conservation law at first order, with a correction component proportional to

\[
H_-=-b^2\lambda_2:
e^{(d+a_2)\cdot X+a_2\cdot\bar X}:.
\]

Coupling and contour conventions follow the preceding report. The even-root perturbation is regular with J_-: d.a_1=-2 gives two positive powers of separation, and the derivative prefactor can remove at most one.

The earlier J_+ test is also recovered: its residue against the even-root vertex is -gamma partial exp[(-d+a_1).X], with gamma=b^2/(1-b^2), while its Yukawa and composite-potential OPEs are regular.

**Verdict:** the new current passes the order-linear-in-Yukawa tests. It has not passed the complete order-quadratic-in-mass test of the original component action.

## 4. What algebra is suggested?

The exponential momenta of J_+,J_- are c=-d and d. Their Gram matrix is

\[
\begin{pmatrix}c^2&c\cdot d\\d\cdot c&d^2\end{pmatrix}
=-\frac1{b^2}\begin{pmatrix}2&-2\\-2&2\end{pmatrix}.
\]

Thus their fractional exponential braiding has affine A1 geometry. Both currents are even. This suggests investigating a bosonic U_q(affine sl2) charge algebra instead of insisting on odd simple generators of U_q(affine sl(2|1)). This is an inference from the candidate currents, not an identification proved by a Cartan-matrix comparison.

The derivative prefactor, contour contacts, antichiral charges, Serre relations and full conservation must all be checked. They may change or invalidate this interpretation. Neither q nor the physical soliton representation has been determined.

With the free stress tensor and no background charge, the candidate charge spins are s_+=1/b^2-1 and s_-=1/b^2. This is an unequal gradation. These free-UV dimensions are not a completed physical Lorentz-spin dictionary for the quantum action.

## 5. Explicit order-mass-squared obstruction to settle

For the composite potential e^{T.X}, d.T=2. Its OPE with J_- has residue

\[
:e^{(d+T)\cdot X}\left[
(a_0\cdot\partial X)(T\cdot\partial X)
-\frac{b^2}{2}\left(T\cdot\partial^2X+(T\cdot\partial X)^2\right)
\right]:.
\]

Subtract the derivative of -(b^2/2)(T.partial X) exp[(d+T).X]. The remaining representative is

\[
\boxed{\frac12:
(T\cdot\partial X)((a_0-a_2)\cdot\partial X)
e^{(d+T)\cdot X}:\ne0.}
\]

Within the ordinary level-one local counterterms, the partial^2 X coefficient fixes the counterterm and leaves this nonzero remainder. This is therefore a genuine failed residue test if the composite-potential coupling is treated as independent.

In the classical component action that coupling is proportional to the product of the two odd-node mass parameters. Two Yukawa insertions occur at the same mass order. Their regulated contact contribution must be included before declaring either cancellation or a quantum anomaly. Bosonization and normal ordering also require a quantum prescription for this composite exponential; the classical coefficient cannot be assumed to give the correct quantum Ward identity.

The script verifies the full residue, the derivative subtraction and the nonzero remainder. It does not compute the double-Yukawa contact contribution. Full conservation is consequently **unclear**, not established.

## Integration

No book integration was performed. Verified: bosonized momenta, pair-dual current, its Yukawa derivative residues, regular even-root OPE, affine-A1 exponential Gram matrix and the explicit composite-potential residue. Modified interpretation: failure of the minimal odd-current ansatz may be a sign that the sought symmetry should not be prescribed as the original superalgebra. Rejected as unsupported: naming a specific Langlands-dual affine algebra or identifying a physical S-matrix from these data alone. Unresolved: the order-mass-squared Ward identity, quantum prescription, complete charge algebra and soliton action. Independent review is required before any new statement is promoted into the book.
