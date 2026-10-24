# Root-system sl(2|1) Toda: boundary closure and candidate quantum vacuum endpoints

Date: 2026-10-24. Model: Codex (GPT-6).

Commit reviewed: `745586e6a6f7bca5ed91a1e9cf89d6ef5df387f0`.

Task: continue boundary closure, then investigate a quantum vacuum sector and its renormalized endpoint limits. The user authorized subagents for this pass. Two independent audits were used, with no further delegation. `AGENTS.md` was read. All work is in new reports and new scratch directories. No pre-existing file, book file, index, or checkout was changed; no commit, push, or render was performed.

## Outcome

There is now a concrete candidate for the vacuum normalization, in addition to the earlier conditional boundary characters. A direct orthogonal field transformation puts the local interacting action on the first sine-Liouville hierarchy. Both leading endpoint primaries have identity-like normalized neutral form factors in the corresponding reference theory. The second primary has a subtle normalization: a simple zero of its reflection coefficient multiplies a simple identity-resonance pole. Their product is finite and generically nonzero.

In an analytic-reflection completion with the normalization specified below, the resulting candidate values are

\[
\boxed{
G_A^{\rm refl}=
[\lambda_1\gamma_f(-t)]^{1/t-1}
\frac{\Gamma(2-t)\Gamma(2-1/t)}{\Gamma(t)\Gamma(1/t)},
\qquad
G_U^{\rm refl}=
\frac{\lambda_0}{1-t}
\left(\frac{\sqrt{\lambda_0\lambda_2}}{-t}\right)^{2/t-2}
\frac{\Gamma(2-1/t)}{\Gamma(1/t)}.
}
\tag{V1}
\]

Here \(\gamma_f(x)=\Gamma(x)/\Gamma(1-x)\). Fractional powers require a chosen logarithm and vacuum/coupling branch. These are **conditional reference-theory values**, not established expectation values of the original quantum Toda model. The direct cutoff calculation does match the second primary's minimal analytic subtraction to the canonical mixed endpoint through the displayed mass order. This materially narrows the normalization problem, but does not establish the all-order endpoint, continuum sector identification, or weak limits on a charge domain.

Boundary and quartic relations also pass a useful algebraic compatibility test. They still require a common regulated three-contour product. No symmetry algebra is identified here.

| Item | Verdict | Evidence and limit |
|---|---|---|
| Formal endpoint characters and anchored boundary weights | Correct in the specified canonical cover | Direct physical-boson/fermion transformation and canonical commutator integration |
| Regulated translated integration cycles | Correct for an invariant finite regulator | Explicit change of variables and even Grassmann Jacobian |
| Automatic continuum vacuum from the bare cosh action | Unclear | Coupling trajectory, counterterms, and global sector remain additional inputs |
| Local reduction to the first sine-Liouville hierarchy | Correct | Explicit orthogonal matrix, momenta, and kinetic normalization |
| Identity-like normalized endpoint form factors | Correct in the 2017 rational family at generic coupling | All-particle induction plus independent finite sign sums |
| Inferring a zero endpoint from the sine-Liouville reflection zero | Wrong | The reflected family has a compensating identity pole |
| Formula (V1) in the stated analytic-reflection completion | Correct conditional derivation | Primary reflection formulas, direct local pole computation, and branch/dimension checks |
| Minimal subtraction equals the displayed canonical endpoint | Correct through the computed order, for \(-1<t<0\) | Direct angular/radial OPE integration |
| Full original-theory scalar endpoint limits or symmetry algebra | Unclear | All-order matching, global domains, and compatible products are not proved |

The “wrong” entry rejects a tempting inference encountered during this investigation, not a sentence in an existing report. No book correction or theorem upgrade is proposed.

## 1. Conventions and evidence boundary

Use full Gaussian covariance \(-\log|x|^2\), ordinary Euclidean area, and no additive constant in the reference Wick subtraction. The local interaction is

\[
S_{\rm int}=\pi^{-1}\int d^2x\,
(\lambda_0Y_0+\lambda_2Y_2+\lambda_1V_1+\kappa P),
\qquad
\kappa=-\lambda_0\lambda_2\epsilon^{2t}/t,
\quad\kappa_R=0.
\tag{V2}
\]

The Yukawa coefficients have mass order one; \(\lambda_1\) and \(\kappa\) have mass order two. This bookkeeping is different from their quantum engineering dimensions. With chiral fields \(X=(\varphi^u,\varphi^v,h)\),

\[
a_0=(ib/\sqrt2,b/\sqrt2,i),\quad
a_2=(ib/\sqrt2,-b/\sqrt2,-i),\quad
t=b^2,\quad T=a_0+a_2,\quad d=-T/t,\quad k=d+T.
\]

Set \(C=a_0\cdot\partial X\), with the analogous barred field. The full endpoints used as candidates are

\[
F_A=e^{-k\cdot(X+\bar X)},\quad
U=e^{(d+a_2)\cdot(X+\bar X)},\quad
F=e^{k\cdot(X+\bar X)},\quad \gamma_E=\frac{t}{1-t},
\]
\[
B_{AA}=\gamma_E^2\lambda_1F_A,
\quad
B_{BB}=\lim_{\epsilon\to0}
\left[t^2\lambda_2U
-\frac{t^2\lambda_0\lambda_2}{t-1}\epsilon^{2t-2}F
-t\lambda_0\lambda_2\epsilon^{2t}C\bar C F\right]_{S_\epsilon}.
\tag{V3}
\]

The old finite-segment contour identity is treated as a candidate to extend, not as an all-order conservation theorem. This pass does not claim to recompute every earlier mass-two continuation integral. It independently checks boundary/product consequences and rederives the particular endpoint counterterms needed for the new vacuum comparison. Passing old scripts is not used as evidence.

## 2. What boundary closure now requires

The formal physical-field transformation

\[
\delta\Phi^u=\sqrt2\pi/b,\quad\delta\Phi^v=0,
\qquad(\chi_0,\chi_2)\mapsto-(\chi_0,\chi_2)
\tag{V4}
\]

with the barred chirality fixed preserves the action: each scalar Yukawa factor changes sign, canceled by its fermion. Both potentials and kinetic bilinears are invariant. It gives, without global bosonization,

\[
G(U)=sU,\qquad G(F)=sF,\qquad G(F_A)=s^{-1}F_A,
\quad s=e^{-2\pi i/t}.
\tag{V5}
\]

The left-anchored split \(\varphi(x)=\frac12[\Phi(x)-\int_{-\infty}^x\dot\Phi]\) and \([\Phi(y),\dot\Phi(x)]=4\pi i\delta(x-y)\) give \([\Phi(y),\varphi(x)]=-2\pi i\Theta(x-y)\). A current with disorder momentum \(qd\), \(q=-1\) for \(A\), \(q=1\) for \(B\), consequently satisfies

\[
[n_L,Q_q]=-2qQ_q,\quad[n_R,Q_q]=0,\quad[N,Q_q]=2qQ_q,
\quad n_{L,R}=T\Phi(L,R)/(2\pi i),\ N=n_R-n_L.
\tag{V6}
\]

Thus \(K=\nu^N\), \(\nu=e^{-i\pi/t}\), has \(KA=s^{-1}AK\), \(KB=sBK\). The fixed right endpoint is a consequence of the left anchor; changing the anchor changes that statement.

If the renormalized endpoints have scalar weak limits in one orbit \(\omega_n(O)=\omega_0(G^nO)\), and the right label is \(r\), their differences are

\[
\Delta B_{AA}=\mathcal A s^{-r}(1-s^N),\qquad
\Delta B_{BB}=\mathcal B s^r(1-s^{-N}).
\tag{V7}
\]

An invariant state instead has zero expectations for these nontrivial characters, if the fields exist on its domain. If one identifies configurations under \(G\), generic-character fields need an enlarged charged/twisted operator space. Undefined quotient observables cannot be assigned a zero VEV by this argument.

The independent [product audit](referee_sl21_product_compatibility_report.md) establishes a sharper remaining test. Include the factor \(\nu/(2\pi i)\) in endpoint symbols \(R_a,L_a\), and write

\[
A\bar A-s\bar AA=R_A-L_A,\quad B\bar B-s\bar BB=R_B-L_B,
\quad A\bar B=s^{-1}\bar BA,\quad B\bar A=s^{-1}\bar AB.
\]

Directly expanding \(D=BA-s^{-1}AB\) in any associative product gives

\[
[D,\bar A]=Bf_A-s^{-2}f_AB=(1-s^{-2})BR_A,
\quad
[D,\bar B]=s^{-1}(f_BA-Af_B)=(s^{-1}-s)AL_B.
\tag{V8}
\]

The last equalities use independently separated endpoint weights. Our own root check reproduces the first equalities. The product audit also tests the four quartic/mixed overlaps and constructs an associative matrix witness with independent endpoint amplitudes. These are consistency results. In particular, the actual order-\(\lambda_2\) three-contour calculation, including the induced finite contact \((s^{-1}-s)A\delta L_B\), is still missing. Setting \(BA=s^{-1}AB\) globally would contradict (V8) in a nontrivial sector.

## 3. A finite regulated orbit, and what its limit needs

In the original real-source variables, choose a finite lattice/mode regulator preserving the displayed same-chirality kinetic bilinears. Transform every \(u_j\) and its Dirichlet data by \(-i\sqrt2\pi n\), and both components of one fermion chirality by \((-1)^n\). The two potentials are unchanged, as is the Yukawa term. For \(2M\) Grassmann variables of that chirality the Berezin Jacobian is \((-1)^{2Mn}=1\).

Consequently any convergent integral on a specified cycle \(\mathcal C_0\) is isomorphic to the integral on \(\mathcal C_n=G^n\mathcal C_0\); covariant insertions have (V5). This is a genuine finite-regulator construction. It requires neither a deformation between cycles with fixed boundaries nor a claimed continuum anomaly cancellation. A finite fermion-kernel congruence gives an independent explicit check.

The continuum prescription is substantive. For \(G_a(0)=-\log a^2+c_{\rm reg}\), Gaussian contraction gives

\[
e^{\pm T\Phi_a}=a^{2t}e^{-tc_{\rm reg}}:e^{\pm T\Phi}:,
\qquad e^{p_i\Phi_a}=:e^{p_i\Phi}:,
\tag{V9}
\]

since \(T^2=-2t\) and the physical Yukawa vectors obey \(p_i^2=0\). A different regulator's pair contact and Gaussian constant must be matched to (V2), including zero finite \(P\) coupling. Holding both bare cosh coefficients fixed does not hold the Wick coefficient \(\lambda_1\) fixed. An explicit mass/field/coupling trajectory is needed.

Homogeneous copies \(\mathcal C_n\) also do not construct the interface sectors with independently chosen left and right vacua. Such sectors, and operators that shift their left labels as in (V6), are needed for a nontrivial boundary-charge representation. The [global audit](referee_sl21_boundary_global_audit_report.md) gives the direct fermionic derivation and shows why elementary Schwinger–Dyson equations leave contact and kinetic moments instead of determining the required amplitudes.

## 4. Independent local map to a reference massive theory

Restrict initially to generic irrational \(-1<t<0\), and set \(b=i\sqrt2\zeta\), \(t=-2\zeta^2\), \(\alpha=\sqrt{1+\zeta^2}\). For full bosons define

\[
\phi=-\Phi^u,\qquad
\vartheta=\frac{\zeta\Phi^v+H}{\alpha},\qquad
\eta=\frac{\Phi^v-\zeta H}{\alpha}.
\tag{V10}
\]

The matrix is orthogonal with determinant one. The kinetic action stays \((8\pi)^{-1}\int[(\partial\phi)^2+(\partial\vartheta)^2+(\partial\eta)^2]\). Direct multiplication gives

\[
Y_{0,2}=e^{\zeta\phi\pm i\alpha\vartheta},\quad
V_1=e^{-2\zeta\phi},\quad P=e^{2\zeta\phi},
\]
\[
F_A=e^{-(2\zeta+\zeta^{-1})\phi},\quad
U=e^{(\zeta+\zeta^{-1})\phi-i\alpha\vartheta},\quad
F=e^{(2\zeta+\zeta^{-1})\phi}.
\tag{V11}
\]

Every displayed endpoint and current prefactor is neutral in \(\eta\). This is a local oscillator decomposition; it does not prove factorization of compactification, spin structures, or zero modes.

At \(\lambda_0=\lambda_2=\lambda\), the remaining interaction is

\[
\frac{2\lambda}{\pi}e^{\zeta\phi}\cos(\alpha\vartheta)
+\frac{\lambda_1}{\pi}e^{-2\zeta\phi}
+\frac{\kappa}{\pi}e^{2\zeta\phi},\qquad \alpha^2-\zeta^2=1.
\tag{V12}
\]

This is the local first sine-Liouville perturbation, with its composite counterterm. In the \(1/(16\pi)\) convention of Fateev's 2017 paper, set \(r=\sqrt2\phi\), \(\theta=\sqrt2\vartheta\), \(B=\zeta/\sqrt2\), \(A=\alpha/\sqrt2\); then \(A^2-B^2=1/2\). It is the first hierarchy, rather than an identification inferred from the root-system name.

Primary sources actually read:

- Fateev, [SIGMA 13 (2017) 080](https://sigma-journal.com/2017/080/sigma17-080.pdf), equations (4.1)–(4.2), printed p. 6, visually inspected. The source describes the exponential counterterm as **“plays the role of the usual counterterm canceling the divergencies coming from the fermion loops”**. Its local action comparison does not specify our global fermionic sectors.
- Lashkevich–Pugai, [arXiv:1610.04926](https://arxiv.org/abs/1610.04926), equations (1.6)–(1.14), (2.13)–(2.17), and section 3. The field identification before (2.16) is explicitly **“up to a normalization factor”**. The 2017 convention is used consistently; no formula is imported by matching a label to a later paper.
- Fateev, [hep-th/0103014v2](https://arxiv.org/abs/hep-th/0103014), equations (24)–(25) and (130)–(132), printed pp. 6 and 33, visually inspected. Its p. 33 states: **“Until now we considered spinless fields and non-compactified field Φ.”** A global compactification is therefore an additional choice.

The candidate reference sector is an analytically normalized sine-Liouville/complex-sinh-Gordon particle vacuum, with a chosen coupling branch. Its wall coefficient, mass relation, and global operator space must be matched. In particular the 2017 dual action writes its wall as \(-\lambda_{\rm LP}e^{-2\zeta\phi}\); here \(\lambda_{\rm LP}=-\lambda_1/\pi\). Dropping this sign or a fractional-power branch would change the proposed comparison. A real-axis integral of the bare cosine potential alone is not a construction of this vacuum.

## 5. Why both normalized primary families are identity-like

In the 2017 source notation use a radial label \(\ell\), to distinguish it from our coupling:

\[
V_{a,\ell,q}=\exp\left(\frac{ia}{2\alpha}\vartheta
-\frac{\ell}{2\zeta}\phi-\frac{iq}{2\alpha}\widetilde\vartheta\right),
\quad g=\frac1{2\alpha^2}=\frac1{2-t},\quad h=g^{-1}.
\]

The explicit map (V11) yields

\[
F_A=V_{0,2h-2,0},\qquad U=V_{-h,-h,0},\qquad
F=V_{0,2-2h,0}.
\tag{V13}
\]

In particular \(-ga_U=1\), so \(U\) is local with respect to the reference particles in this convention. This does not identify their charge with our boundary label \(N\).

The rational form-factor functions use \(p_\pm=1-h+\ell\mp a\). Thus the identity, \(F_A\), and \(U\) have parameters

\[
p_0=(1-h,1-h),\quad p_A=-p_0,\quad p_U=p_0+(0,-2h).
\]

An independent all-particle induction proves the stronger statement

\[
J_{K,L,p_{L-K}}=0\quad(0\le K\le L,\ K+L>0),
\qquad p_r=(1-h+r,1-h-r).
\tag{V14}
\]

For \(K=0\), the first cosine in the source's product formula has argument \(-\pi/2\). The recursion for a particle variable reduces the constant term to \((K-1,L,p_{r+1})\), and each pole term to \((K-1,L-1,p_r)\). Both vanish by induction. The only exception is \(K=L=1\), whose lower term is \(J_{00}=1\); its residue vanishes because \(p_+=p_-\). Inversion of the pair kernels gives the result for \(p_A\). Shifting \(p_-\) by \(-2h\) multiplies every dagger-current term by minus one and gives the result for \(p_U\).

Therefore both endpoint families have zero normalized positive-particle neutral form factors and normalized vacuum value one. The additional rapidity factor for \(U\) is \(\prod x_i/\prod y_j\), which cannot undo an exact zero at regular parameters. A divergent normalization times a parameter-zero form factor could, however, define a different limiting operator. The full [independent proof and fresh checks](../referee_scratch/sl21_product_compatibility_20261024/ff_identity_note.md) keep this distinction explicit.

## 6. Reflection normalization, including the zero times pole

For \(F_A\), the Liouville wall \(e^{-2\zeta\phi}\) has background charge \(Q_L=-\zeta-(2\zeta)^{-1}\) in our kinetic normalization. Its exponent is \(2Q_L\), the reflection image of zero. Substituting the rank-one root of length squared two into source equations (24)–(25), with \(b_L=-\sqrt2\zeta\) and \(\pi\mu=\lambda_1\), gives the first formula in (V1). This is a source-based reflection normalization, not a Gaussian expectation or an assumption that the bare exponential is the identity.

For \(U\), let \(Q=(2\zeta)^{-1}\), \(a_*=\zeta+\zeta^{-1}\), and keep angular momentum \(\xi=-\alpha\). The rank-one reduction of source (132), with \(C_*=\lambda/(-t)\) at symmetric Yukawa coupling, is

\[
\mathscr A(a,\xi)=C_*^{(a-Q)/\zeta}
\frac{\Gamma(1-2(a-Q)\zeta)\Gamma(1-(a-Q)/\zeta)}
{\Gamma(1/2-(a-Q)\zeta+\alpha\xi)
\Gamma(1/2-(a-Q)\zeta-\alpha\xi)},
\quad R(a,\xi)=\frac{\mathscr A(2Q-a,\xi)}{\mathscr A(a,\xi)}.
\tag{V15}
\]

The source uses the shifted momentum in the power of \(C_*\). Its reflection orientation has been checked from the actual page.

Set \(a=a_*+\delta\), so the reflected exponent is \(-\zeta-\delta\). One reflected denominator is \(\Gamma(\zeta\delta)\); all remaining factors are regular at generic irrational coupling after the indicated cancellations. Since \(1/\Gamma(\zeta\delta)=\zeta\delta+O(\delta^2)\),

\[
R(a_*+\delta,-\alpha)
=\delta\,\zeta C_*^{2/t-2}
\frac{\Gamma(2-1/t)}{(t-1)\Gamma(1/t)}+O(\delta^2).
\tag{V16}
\]

The zero must be combined with the reflected operator. Its contraction with \(Y_0\) is exactly

\[
V_{-\zeta-\delta,-\alpha}(0)Y_0(x)
\sim |x|^{-2+2\zeta\delta}e^{-\delta\phi}(0).
\]

The minus action insertion gives

\[
-\frac{\lambda_0}{\pi}\int_{|x|<R}d^2x\,|x|^{-2+2\zeta\delta}
=-\frac{\lambda_0}{\zeta\delta}R^{2\zeta\delta}.
\tag{V17}
\]

Its residue is the identity operator, rather than an unknown condensate. Within the given local perturbation class this is the unique nontrivial generic-coupling resonance. To see it, let \(n_0,n_2,m\) count \(Y_0,Y_2,V_1\), and let \(l\ge0\) be a descendant level. Put \(n=n_0+n_2\), \(D=n_0-n_2-1\), \(v=\zeta^2\). Equality of UV dimensions requires

\[
D^2+n+2m+l+v\{D^2-(n-1-2m)^2+4m\}=1.
\]

For generic irrational \(v\), the constant equation forces \(n+2m+l\le1\). Apart from no insertion, the sole solution is \(n_0=1,n_2=m=l=0\). The unreflected \(U\) itself has no nontrivial resonance in this class, as detailed in the [independent reflection audit](../referee_scratch/sl21_boundary_global_audit_20261024/reflection_zero_audit.md).

Multiplying (V16) by the pole (V17) gives the second formula in (V1). A constant angular shift extends it to unequal Yukawa coefficients, with \(C_*=\sqrt{\lambda_0\lambda_2}/(-t)\) and the displayed \(\lambda_0\) prefactor. At symmetric coupling it simplifies to

\[
G_U^{\rm refl}=C_*^{2/t-1}\frac{\Gamma(1-1/t)}{\Gamma(1/t)}.
\tag{V18}
\]

An arbitrary finite part of the reflected family is killed by the simple zero; the residue is essential. Subtracting the pole first and then setting \(R=0\) is a different prescription and loses this value. Intersecting gamma poles at exceptional couplings need a joint limiting calculation; no value at such a coupling is asserted.

This calculation assumes that the massive meromorphic operator family obeys the source reflection relation and uses its minimal Wick normalization. It does not derive that identity for our original fermionic path integral. The normalized form-factor result supplies a separate consistency check, not proof of the normalization assumption.

## 7. Direct match to the canonical endpoint subtraction

The normalization comparison can be tested locally. With the full OPE coefficient fixed to plus one by the fermion convention,

\[
U(0)Y_0(x)=|x|^{2t-4}
\bigl(F+xCF+\bar x\bar CF+x\bar x C\bar C F+\cdots\bigr).
\tag{V19}
\]

Angular integration keeps the scalar term and \(x\bar x\) term. Multiplication by \(-\lambda_0/\pi\) gives cutoff divergences

\[
\frac{\lambda_0}{t-1}\epsilon^{2t-2}F
+\frac{\lambda_0}{t}\epsilon^{2t}C\bar C F.
\]

For \(-1<t<0\) the next angularly invariant term vanishes as \(\epsilon^{2t+2}\). Thus the minimal analytic Wick representative is

\[
U_{\min}=\lim\left[
U-\frac{\lambda_0}{t-1}\epsilon^{2t-2}F
-\frac{\lambda_0}{t}\epsilon^{2t}C\bar C F\right],
\qquad B_{BB}=t^2\lambda_2U_{\min}
\tag{V20}
\]

through the displayed order. This rederives both subtractions, their signs, and their relative coefficients without using the old scripts. It matches the canonical endpoint, including the finite choice of its leading coefficient \(t^2\lambda_2\). A different finite multiple of \(U\) or additional finite mixing would change the comparison. Higher-order source collisions and all-order endpoint improvements have not been fixed by (V20).

## 8. Branches and actual endpoint limits

The candidate values respect the required branch characters. Under the period in local variables, \(\delta\phi=i\pi/\zeta\), \(\delta\vartheta=\pi/\alpha\). Pulling the translated integration cycle back to a fixed cycle corresponds to continuous logarithm shifts

\[
\log\lambda_0\mapsto\log\lambda_0-2\pi i,\quad
\log\lambda_2\mapsto\log\lambda_2,\quad
\log\lambda_1\mapsto\log\lambda_1+2\pi i.
\tag{V21}
\]

The endpoint formula depends on \(\lambda_0^{1/t}\lambda_2^{1/t-1}\), and the first on \(\lambda_1^{1/t-1}\). They therefore acquire \(s\) and \(s^{-1}\), respectively. Reflection to an identity in a selected branch is compatible with a nontrivial character on the whole cover. A single-valued constant assigned to every branch would be incompatible.

If the reference bootstrap operators are regular with these finite normalizations and its particle states are complete, (V14) and crossing make the two primaries scalar identities on that particle domain. Under that additional identification their translated weak limits have values (V1). The expected leading endpoint constants would be

\[
\mathcal A_{\rm cand}=\gamma_E^2\lambda_1G_A^{\rm refl},\qquad
\mathcal B_{\rm cand}=t^2\lambda_2G_U^{\rm refl}.
\tag{V22}
\]

This supplies concrete benchmarks for a chosen quantum completion. It does not construct the original-theory weak operator limit. In that theory one must first define the full renormalized endpoint, then control infinite volume and its translation to the left/right vacuum on a common domain. One-point values alone do not prove this. A local free spectator is neutral here, but its global constraints still require matching.

Nor can the order-mass-squared remainder be carried unchanged to infinite volume. The powers \(\lambda_1^{1/t-1}\) and \(\lambda^{2/t-1}\) are nonanalytic in the couplings. There is no uniform small-mass justification for discarding every higher-order correction before taking the endpoint limit.

Finite nonzero endpoint constants can be absorbed by independent constant rescalings of each chiral/antichiral charge pair once a common representation exists. Their exact numerical values are therefore not the only obstacle to an algebra. The substantive remaining obstacles are the global interface sectors, the actual weak limits of full endpoints, all-order conservation, and the common associative product required by (V8). The earlier four-fold identities remain formal continued products until this compatibility is demonstrated.

## 9. Checks and next calculation

New root evidence lives in `referee_scratch/sl21_vacuum_limits_20261024/`:

- `vacuum_checks.py` / `.out`: **49 PASS, 0 FAIL**. Independent local transformation, screening/endpoint labels, dimensions, fermionic period, Wick conversion, finite fermion congruence, and three-word expansion.
- `reflection_endpoint_checks.py` / `.out`: **36 PASS, 0 FAIL**. Independent source specialization, simple-zero scaling, identity-pole residue, radial counterterms, generic resonance test, dimensions, and monodromies. Irrational probes include \(t=-\sqrt2/3\) and \(-\sqrt3/4\). The initial output is retained: one symbolic substitution needed expansion before replacing \(\zeta^2\); no scientific formula was changed.

The two independent audits add **49** boundary checks, **60** product checks, **37** form-factor checks, and **28** reflection/resonance checks. The finite checks supplement the derivations; their count is not evidence of a constructed quantum vacuum. The root independently rechecked the local transformation, projected three-product expansion, rank-one reflection values, resonance pole, and endpoint subtraction rather than treating those audit outputs as proof.

The most useful next calculation has two explicit targets. First evaluate the order-\(\lambda_2\) three-contour defect in (V8), with the binary finite matching inserted on every nested face. Second specify a continuum completion and match its full endpoint operator to the analytic-reflection family, starting from (V20) and retaining the coupling branches in (V21). The interface sector with fixed right vacuum must be part of that construction. Equation (V1) is now a quantitative normalization test for it.

Final preservation check: all **32** recorded input hashes are unchanged. Both tracked working-tree and staged diffs are empty. The only new non-ignored status entries are the three reports from this pass; the new scratch directories are excluded by the repository's existing rules.

## Integration

None performed. **Verified:** formal boundary characters, local hierarchy map, normalized form-factor degeneracy, reflection zero and compensating pole, and the displayed cutoff matching. **Extended:** finite regulated vacuum orbit and conditional analytic-reflection endpoint values with explicit branches. **Rejected:** assigning zero to \(U\) from its reflection zero alone, or assigning an algebra from local hierarchy matching. **Unresolved:** the original continuum vacuum/domain identification, full weak endpoint limits, all-order Ward identities, and common associative/Serre products. Existing reports retain their restricted status.
