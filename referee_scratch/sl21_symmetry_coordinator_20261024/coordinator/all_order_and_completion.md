# Independent all-order and completion attempts (working)

Date 2026-10-24, Codex GPT-6, reviewed commit 745586e6a6f7bca5ed91a1e9cf89d6ef5df387f0.

## 1. Local anomaly selection for both candidates

Take n0,n2 Yukawas and e even-root insertions. Let a=n0−e,c=n2−e,S=a+c,F=n0−n2,M=n0+n2+2e. The physical bosonic momentum shift is a p0+c p2. Its squared norm is −2tac; minimal chiral fermion weight at charge F is F²/2 (fill modes of a single fermion with distinct derivatives). Thus the antichiral base weight is tac+F²/2. B's holomorphic base weight is 1/t−S+tac+F²/2, with its derivative supplied as a descendant; A's is 1/t+S+tac+F²/2. Each coupling has left/right weight M/2−et. These statements follow directly from Gaussian h=−p²/2 and the Yukawa/even-root dimensions, not from the earlier anomaly script.

Required descendant levels for divergence of B:

ell_B=1−M/2+S−F²/2+t(e−ac), ellbar=1−M/2−F²/2+t(e−ac).

For A replace ell_B by ell_A=−M/2−S−F²/2+t(e−ac). At irrational t, integral descendant levels force e=ac. Then ellbar=1−S(S+1)/2, ell_B=1−S(S−1)/2 and ell_A=−S(S+3)/2. Nonnegative levels and n0=a+ac>=0,n2=c+ac>=0 force only the single-Yukawa possibilities for B and a=c=−1,e=1 (one even-root insertion) for A, apart from the unperturbed term. For B, S must lie in {−1,0,1}; for A the only viable nonzero insertion is S=−2. This is a complete integer argument; the finite enumeration in extension_checks.py merely checks it.

Verdict: correct conditional exclusion of additional ordinary local, mass-analytic, scale-homogeneous Gaussian anomalies. It is not a proof of a quantum action principle for products of disorder-line operators. The missing induction assertion is that every renormalized Ward breaking is exhausted by this local Fock-field class, with simultaneous contour/end contacts separately controlled. Nonanalytic coupling terms, line-supported terms and interface-domain losses are not covered.

## 2. A concrete failed real reference cycle

For t=−2zeta², alpha²=1+zeta², the locally transformed reference potential is
V=lambda0 exp(zeta phi+i alpha theta)+lambda2 exp(zeta phi−i alpha theta)+lambda1 exp(−2zeta phi).
On lambda0=lambda2=lambda>0, lambda1>0, a real theta interval around pi/alpha has cos(alpha theta)<−1/2. On that interval V<=−lambda exp(zeta phi)+lambda1 exp(−2zeta phi), tending to minus infinity as phi→+infinity. The zero-mode integral of exp(−volume V) diverges on an open region, not just a measure-zero ray. At a finite lattice with fixed boundaries the gradient action grows only quadratically along such a large common-field region, and cannot dominate the negative exponential. Therefore the naive positive real reference-boson Euclidean integral is not a completion. This does not disprove a complex contour, a fermionic construction, or the original source action with its different real slice. The locally orthogonal change involving the bosonized fermion does not transport its global cycle automatically.

A different construction described in Fateev (2017), (4.2)-(4.3), has a fermion coupled to a boson with cosh potential and a Thirring interaction. It is a serious candidate reference realization, not a derivation of the original spin structures, compactification, or Ward endpoint operators. In the inspected source (4.1), the bosonized and fermionic forms are connected for small coupling; the actual sentence is quoted in the final report only to its permitted extent. No unsupported equivalence theorem is imported.

## 3. Reflection and normalization do not supply interfaces

Direct source (131)-(132) zero/pole computation and Gaussian transport are in intake_checks.py. They give a reference normalization only. If the reflected meromorphic family, crossing and sector completeness are stipulated, vanishing positive-particle form factors can identify an operator with G times I on that particle sector. They neither generate H_(l,r) for l!=r nor define an intertwiner between sectors. A direct sum of scalar reference copies permits either zero intertwiners or many arbitrary shifts; the scalar data alone cannot select the original nonlocal current.

At the level of positive-particle form factors, O and O+c(lambda)I are indistinguishable. For the proposed endpoint c(lambda)=eta G(lambda) has the same dimension and branch character and is nonanalytic in the mass couplings. The meromorphic reference normalization can forbid this change, but the unproved all-order equality of the Toda Ward endpoint with that reference family cannot be replaced by finite-order UV matching. This is freedom of the bridge, not permission to alter a fixed source formula.

The required endpoint convergence is weak convergence of full renormalized insertions in matrix elements on a common invariant dense interface domain, with both UV and infinite-volume limits controlled. A homogeneous VEV and a finite perturbative remainder do not suffice. For example I+exp(iR)Q with a vacuum-annihilating Q has constant vacuum expectation yet need not have a weak limit. This is a logical counterexample to inference from a VEV, not a proposed Toda observable.

## 4. Coproduct directly from the spatial anchor

For a current with chiral-minus-antichiral momentum q d, split the left anchored field at c:
varphi_L(x)=varphi_c(x)−(1/2) int_L^c dotPhi, x>c.
Thus every displayed component of that current in the right region acquires
Omega_q(L,c)=exp[−q d.int_L^c dotPhi/2].
The neutral derivative prefactor changes only at the split point, outside separated supports. The formal regional rule is Q_q(L,R)=Q_q(L,c)+Omega_q(L,c) Q_q(c,R), or Delta Q_q=Q_q tensor 1+Omega_q tensor Q_q when independent regional operator spaces exist. Concatenating three disjoint intervals makes Omega multiplicative and yields the same three terms in either grouping.

This is a derivation of the string factor, not an identification Omega_q=nu^(±N). The former is initially an exponential of integrated canonical momentum; N is a field difference. Identifying their actions on a chosen space requires the braiding/sector realization. For a separated target vertex with momenta (l,r), canonical commutators give braiding character exp[i pi q d.(l+r)] in the chosen x>y order. That follows from the anchored CCR, not from a presumed algebra. Matching it to a topological Cartan action is a separate link.

## 5. Conditional presentation, not an established algebra

If, on one associative invariant domain, K A=s^-1 A K, K B=s B K and the corresponding barred weights hold, diagonal defects are a(1−K²), b(1−K^-2), the off-diagonal defects vanish with factor s^-1, both chiralities obey the quartics with 1+s+s^-1, and a,b!=0, put q=nu, q²=s, K0=K^-1,K1=K. Then
E0=A, E1=B,
F0=barA K^-1/[a(q−q^-1)], F1=barB K/[b(q−q^-1)]
give [Ei,Fj]=delta_ij(Ki−Ki^-1)/(q−q^-1), Cartan matrix [[2,−2],[−2,2]], and K0K1=1. For example [A,barA K^-1]=(AbarA−s barAA)K^-1=a(K^-1−K). The cross commutators cancel directly. Moving K through each barred quartic gives a common nonzero factor for all words, so the negative-generator Serre relations follow conditionally as well.

This would be an action of the level-zero quantum affine sl2 presentation (without the degree generator), possibly through a proper quotient. It is not established here, cannot be assumed from the exponential Gram matrix, and does not fix a physical sign of q. The coproduct also needs the independent Omega/Cartan identification above. No claim of faithfulness or exclusion of further symmetries is made.
