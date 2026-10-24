# Root-system sl(2|1) Toda: global boundary audit

Date: 2026-10-24. Model: Codex (GPT-6).

Commit reviewed: `745586e6a6f7bca5ed91a1e9cf89d6ef5df387f0`.

Task: independently audit sections 4–8 of `briefs/referee_sl21_boundary_orbits_report.md`, especially the anchored endpoint commutators, period covariance, invariant-state selection rule, and what boundary closure actually requires. `AGENTS.md` was read. This report and `referee_scratch/sl21_boundary_global_audit_20261024/` are new. No previous script was imported or rerun, no existing file was edited, and no subagent or git mutation was used.

## Verdicts

| Claim | Verdict | Evidence and qualification |
|---|---|---|
| Endpoint characters under the formal action period | **Correct** | Direct two-boson/fermion transformation; no global auxiliary-boson identification is needed |
| The chosen left anchor fixes the right endpoint | **Correct in the stated canonical realization** | Integrating the canonical delta function gives the quoted signs; the choice of preserved endpoint changes with the anchor |
| A period-invariant state has zero generic-character endpoint expectation | **Correct if the endpoint is an operator on that domain** | The one-line automorphism argument applies to finite renormalized expectations; a quotient theory needs a separate operator definition |
| A local free boson forces these endpoint VEVs to vanish | **Wrong** | The endpoints are neutral under that direction; this negative claim was already rejected in the audited report |
| A finite regulated period orbit can be constructed directly | **Correct under an explicit invariant regulator** | Change of integration variables and a two-component Grassmann Jacobian establish an isomorphism of regulated integrals |
| A fixed bare symmetric-cosh regulator automatically realizes the required quantum prescription | **Unclear and insufficiently specified** | Both bosonic exponentials acquire a regulator-dependent Wick factor; the finite composite subtraction and even-root coupling need matching |
| Schwinger–Dyson equations alone determine the endpoint amplitudes | **Wrong as a deduction from the displayed EOM** | Weighting an EOM by the dual vertex produces a contact term and independent kinetic/composite moments |
| Exact numerical amplitudes are necessary before any algebra can be identified | **Too strong** | Once scalar limits, common sectors, and finite nonzero amplitudes are proved, constant charge rescalings can absorb the two diagonal constants |
| An established quantum symmetry algebra follows | **Unclear** | Infinite-volume limits, operator domains, all-order endpoints, nonvanishing, and a common associative product remain open |

The audited report carefully states its boundary conclusions conditionally. No correction to its local formulas is required. The additional distinctions below concern extensions of those formulas.

## 1. Direct fermionic period, without global bosonization

Use the two physical bosons and momenta

\[
p_0=(ib/\sqrt2,b/\sqrt2),\quad
p_2=(ib/\sqrt2,-b/\sqrt2),\quad
t=b^2,\quad T=p_0+p_2,\quad d=-T/t,\quad k=d+T.
\]

The full Yukawa fields are
\(Y_i=i:\chi_i\bar\chi_i e^{p_i\Phi}:\).
Set

\[
\delta\Phi=(\sqrt2\pi/b,0),\qquad
\chi_0\mapsto-\chi_0,\quad\chi_2\mapsto-\chi_2,
\qquad\bar\chi_i\mapsto\bar\chi_i.
\tag{G1}
\]

Direct contraction gives
\(p_i\delta\Phi=i\pi\),
\(T\delta\Phi=2\pi i\), and
\(d\delta\Phi=-2\pi i/t\).
The scalar Yukawa sign cancels its fermion sign. Both same-chirality kinetic bilinears are invariant, and both bosonic potentials \(e^{\pm T\Phi}\) are invariant. This constructs the formal automorphism without treating \(h\) as a globally defined compact boson.

For \(s=e^{-2\pi i/t}\), the direct endpoint expressions give

\[
\mathscr U=i:\chi_2\bar\chi_2e^{(d+p_2)\Phi}:\mapsto s\mathscr U,
\quad F=:e^{k\Phi}:\mapsto sF,
\quad F_A=:e^{-k\Phi}:\mapsto s^{-1}F_A.
\tag{G2}
\]

The fermion minus in the first formula is necessary. The current prefactors and local derivative subtractions are invariant, so every displayed term in the renormalized \(BB\) endpoint has character \(s\). Reference-scale Wick factors do not change these characters.

The global interpretation remains additional data. On a field cover, (G1) may act as a global automorphism. If instead configurations are identified under (G1), generic-character fields are not ordinary single-valued observables on the quotient. They would need a charged or twisted operator space. One must not alternate between those two constructions while using the same endpoint VEV notation.

## 2. Canonical brackets and the anchor

For the audited left anchor,

\[
\varphi(x)=\frac12\Phi(x)-\frac12\int_a^x\dot\Phi(y)dy,
\qquad
\bar\varphi(x)=\Phi(x)-\varphi(x),
\]

where \(a<L<x<R\), and
\([\Phi^a(y),\dot\Phi^b(z)]=4\pi i\delta^{ab}\delta(y-z)\).
Integration gives

\[
[\Phi(y),\varphi(x)]=-2\pi i\Theta(x-y),
\qquad
[\Phi(y),\bar\varphi(x)]=2\pi i\Theta(x-y)
\tag{G3}
\]

for points to the right of the finite anchor. For any vertex with chiral momenta \((l,r)\), the left-endpoint bracket depends only on \(l-r\). The two currents and all tails displayed in the mixed report have
\(l-r=qd\), where \(q=-1\) for \(A\) and \(q=1\) for \(B\). Since \(Td=2\),

\[
[n_L,Q_q]=-2qQ_q,\quad [n_R,Q_q]=0,\quad
[N,Q_q]=2qQ_q,
\quad n_{L,R}=T\Phi(L,R)/(2\pi i),\ N=n_R-n_L.
\tag{G4}
\]

The fermion factors commute with the separated physical boson endpoint. Derivative contacts have separated support. No global \(h\)-zero-mode assumption enters this calculation. A new independent Klein operator can change a global realization and must be included explicitly if chosen.

The BCH exponents at the left endpoint are
\(-2\pi i q,kd\) for \(F\), and
\(-2\pi i q(d+p_2)d\) for \(\mathscr U\).
Here \(kd=2-2/t\) and \((d+p_2)d=1-2/t\); both exponentiate to \(s^{-2q}\). Thus the quoted endpoint exchange phases are correct.

The preserved side is an anchor convention, not an invariant property of the yet-unconstructed charge. For example, add
\(\eta\int_a^{b_*}\dot\Phi/2\) to \(\varphi\) and subtract it from \(\bar\varphi\), with \(b_*>R\). Then

\[
[n_L,Q_q]=2q(\eta-1)Q_q,\qquad
[n_R,Q_q]=2q\eta Q_q,
\qquad [N,Q_q]=2qQ_q.
\tag{G5}
\]

At \(\eta=1\), the split is anchored on the right and fixes the left endpoint. The audited formula corresponds to \(\eta=0\). A boundary-constrained phase space would also require the appropriate constrained brackets instead of importing (G3).

## 3. Regulated cycle orbit and its ultraviolet qualification

In the real-source variables, take a regulator with bosonic fields \(u_j,v_j\), two components of \(\psi_{+,j}\), and two components of \(\psi_{-,j}\) at each site or retained mode. The source action has bosonic factors \(e^{\pm\sqrt2u_j}\) and Yukawa scalar factor \(e^{-u_j/\sqrt2}\). The simultaneous transformation

\[
u_j\mapsto u_j-i\sqrt2\pi n,
\qquad\psi_{+,j}\mapsto(-1)^n\psi_{+,j},
\qquad\psi_{-,j},v_j\text{ unchanged}
\tag{G6}
\]

preserves the displayed action if boundary data are transformed too. A constant translation preserves bosonic differences. The two-component fermionic kinetic bilinears are even under the sign. The Yukawa scalar contributes \((-1)^n\), canceling the left fermion sign.

At a finite regulator there are \(2N\) left Grassmann variables. Their transformation determinant is \((-1)^{2Nn}=1\); its Berezin inverse is also one. Consequently a convergent integral on \(\mathcal C_n=\mathcal G^n\mathcal C_0\) is isomorphic to the one on \(\mathcal C_0\), and covariant insertions acquire their characters. This is a change of variables between specified cycles. It does not claim the cycles are deformable into one another with fixed boundary data, nor that they already define distinct continuum vacua.

The finite determinant calculation covers the chosen untwisted regulator. It does not prove a continuum implementation in all spin sectors, a gauging of the period, or absence of global obstructions in an orbifold. A regulator that changes the displayed chiral kinetic structure needs an independent symmetry check.

There is also a concrete matching issue. Suppose the Gaussian coincident subtraction is

\[
G_a(0)=-\log a^2+c_{\rm reg}.
\]

The elementary bare/Wick identity is
\(e^{p\Phi_a}=\exp[p^2G_a(0)/2]:e^{p\Phi}:\).
Because \(p_i^2=0\) but \((\pm T)^2=-2t\),

\[
P_{\rm bare}=a^{2t}e^{-tc_{\rm reg}}P_{\rm Wick},
\qquad
(V_1)_{\rm bare}=a^{2t}e^{-tc_{\rm reg}}(V_1)_{\rm Wick}.
\tag{G7}
\]

The earlier circular prescription sets the subtraction constant to zero and matches the unlike-Yukawa pair extension to
\(\kappa=-\lambda_0\lambda_2a^{2t}/t\), with \(\kappa_R=0\). A generic regulator has a different Gaussian constant and a different short-distance pair kernel. Their finite local differences must be matched, and any residual finite \(P\) coupling must be tuned away. This matching cannot be inferred from the fixed bare cosh potential alone.

The even-root coupling is independently Wick-normalized in the earlier prescription. Holding a bare \(\mu_1=\mu_0\mu_2\) and all three bare masses fixed generally fails to hold its Wick coefficient fixed when \(a\to0\). A coupling trajectory or compensating field/mass renormalization is needed before the classical symmetric-cosh vacuum and its mass can be carried to the quantum problem. This is already compatible with the warning at the end of section 6 of the Ward-mass-squared report.

All these bosonic potential counterterms are period invariant. Matching them can preserve the cycle orbit (G6); it is the continuum dynamics and normalization, rather than the finite orbit identity, that remain unproved.

## 4. Selection rule, scalar boundary limits, and constants

If \(\mathcal G\) acts on a state domain, the renormalized endpoint exists there, and its renormalization is equivariant, an invariant state satisfies

\[
\omega(B)=\omega(\mathcal G B)=s^{\pm1}\omega(B).
\]

For \(s\ne1\), this forces the two diagonal expectations to vanish. If the endpoint does not belong to the domain—for example, because one has taken an unsuitable quotient—this equation is not a proof of a zero VEV; the VEV is undefined. The distinction matters.

On a chosen cover orbit, define
\(\omega_{\alpha,n}(O)=\omega_{\alpha,0}(\mathcal G^nO)\).
If the full endpoints have scalar weak cluster limits in these states, their characters imply
\(A_\alpha s^{-n}\) and \(B_\alpha s^n\). Fixing the right label \(r\) in the left-anchored convention reproduces

\[
\Delta B_{AA}=A_\alpha s^{-r}(1-s^N),\qquad
\Delta B_{BB}=B_\alpha s^r(1-s^{-N}).
\tag{G8}
\]

Different orbit labels at the two ends obstruct this single-amplitude expression. Ordinary one-point functions do not establish the required weak operator limits. The perturbative order-two endpoints also cannot be substituted for their full versions after an uncontrolled infinite-volume limit.

Exact numerical values of \(A_\alpha,B_\alpha\) are not always required to identify an abstract algebra up to normalization. Once both are proved finite and nonzero on a common sector, rescalings
\(Q_A\mapsto aQ_A,\ \bar Q_A\mapsto\bar a\bar Q_A\) and
\(Q_B\mapsto bQ_B,\ \bar Q_B\mapsto\bar b\bar Q_B\)
multiply their mixed diagonal coefficients by \(a\bar a\) and \(b\bar b\), respectively. Their topological weights and homogeneous same-chirality relations are unchanged. One can therefore normalize two nonzero constants without calculating their exact values. This observation presupposes the other mixed relations, a common associative product, and admissible rescalings on the representation. It supplies none of them. A zero amplitude cannot be normalized this way, and physical charge normalizations still require its value.

## 5. Why elementary Schwinger–Dyson identities do not close the amplitudes

The ordinary-area continuum action has

\[
\partial\bar\partial\Phi
=\lambda_0p_0Y_0+\lambda_2p_2Y_2
-\lambda_1T V_1+\kappa T P
\tag{G9}
\]

as its equation of motion away from other insertions, with all products interpreted in a compatible prescription. The scalar products computed directly from the physical roots are
\(dp_0=dp_2=1\), \(dT=2\), and
\((p_0-p_2)p_0=t\), \((p_0-p_2)p_2=-t\).
Thus, in a translation-invariant vacuum where derivative expectations vanish and the displayed products exist, the unweighted equations give

\[
\lambda_0\langle Y_0\rangle
=\lambda_2\langle Y_2\rangle
=\lambda_1\langle V_1\rangle-\kappa\langle P\rangle.
\tag{G10}
\]

These are unweighted Yukawa VEVs. The desired endpoint contains the dual vertex. Multiplying (G9) by that vertex is a coincident composite operation, not a free consequence of (G10).

At a finite bosonic regulator the obstruction is explicit. Write the Gaussian action as \(\tfrac12\Phi K\Phi\), absorb quadrature weights into potential coefficients \(g_j\), and take the bare insertion \(E_i=e^{d\Phi_i}\). Integration by parts in \(\Phi_i\), with fixed fermions, gives exactly

\[
\begin{split}
\langle E_i\,d\cdot(K\Phi)_i\rangle
&+g_0\langle E_iY_{0,i}\rangle
+g_2\langle E_iY_{2,i}\rangle\\
&-2g_1\langle E_iV_{1,i}\rangle
+2g_P\langle E_iP_i\rangle
=d^2\langle E_i\rangle.
\end{split}
\tag{G11}
\]

The right side is the contact from differentiating the insertion. The first term is an independent kinetic moment. Neither is removed by translation invariance. In the continuum, these coincident products must also be matched to the endpoint subtractions; for example \(E_dY_i\) has a \(|r|^{-2}\) short-distance vertex contraction. Discarding the contact or replacing all coincident products by already-normalized endpoint symbols silently changes the prescription.

An independent exact one-dimensional probe checks the same contact mechanism. For
\(I(c)=\int_{\mathbb R}dx\,e^{cx-Ae^x-Be^{-x}}\), integration by parts yields
\(AI(c+1)-BI(c-1)=cI(c)\).
The script verifies this with the independently integrated Bessel representation and confirms that dropping the last term gives a nonzero error. This probe is not a Toda vacuum approximation; it checks the logical step in the proposed Ward shortcut.

Additional reflection identities or exact vacuum information could close the problem, but they must be derived in the same operator domain and normalization. In particular, a proposed scalar reflection-to-identity formula for a generic-character endpoint must keep track of the vacuum or coupling branch; it cannot be a single-valued identity on the whole formal period cover.

## Evidence and integration

The new [check script](../referee_scratch/sl21_boundary_global_audit_20261024/global_checks.py) gives [49 PASS, 0 FAIL](../referee_scratch/sl21_boundary_global_audit_20261024/global_checks.out). It independently checks root pairings, direct period phases, the two-component finite Grassmann Jacobian, anchor dependence, every displayed current/tail topological weight, bare/Wick potential factors, EOM projection coefficients, and the weighted integration-by-parts contact. The analytic global/domain qualifications are established by the derivations above, not by the number of checks.

No book integration was performed. **Verified:** the audited local endpoint characters, chosen-anchor brackets, and conditional invariant-state selection rule. **Extended:** direct fermionic regulated cycle construction, its matching requirements, anchor dependence, and the finite-cutoff weighted Schwinger–Dyson identity. **Rejected as unsupported:** an automatic continuum vacuum from the bare cosh regulator, amplitude extraction from unweighted EOM alone, and algebra identification from period covariance. **Refined:** exact amplitude values can be replaced by a proof of finite nonzero scalar limits if only an algebra up to generator normalization is sought. No theorem upgrade is proposed.
