# Fundamental sl(2|1) scattering: construction and Toda identification test

Date: 2026-10-23. Model: Codex (GPT-6).

Book baseline reviewed: `745586e6a6f7bca5ed91a1e9cf89d6ef5df387f0`.

Task: pursue the four steps proposed in the conversation: (1) construct dual/fundamental amplitudes and crossing, (2) determine the rapidity/gradation dictionary, (3) construct scalar factors and examine bootstrap closure, (4) identify a Toda field theory by independent dynamical tests.

This is a research/referee report, not new book text. No existing repository file was edited, no book script was imported, and no git command that changes the index or working tree or Quarto render was run. All calculations are in `referee_scratch/sl21_soliton/`. The formulas below require an independent referee before any theorem or physical identification is claimed.

## Outcome

There is a concrete algebraic scattering candidate on a fundamental triplet and its dual. All four channels have been constructed independently from the affine intertwining equations. A choice of scalar normalization satisfies dual crossing and braiding unitarity exactly. An explicit rapidity gradation converts these identities into a common charge-conjugation crossing formula. The graded Yang–Baxter equations were checked for all eight assignments of fundamental/dual species.

This does **not** yet identify a Toda soliton S-matrix. The candidate fails positive-metric unitarity in the tested physical gradation. It also has scalar/gradation ambiguities which the algebra cannot fix. The usual root-system sl(2|1) Toda model has a massless bosonic direction, and a different, one-boson fermionic Toda model in the literature realizes the R-matrix as an **auxiliary monodromy** intertwiner. That is not evidence that its asymptotic solitons transform in the same triplet.

The first attractive fusion is informative: its four-dimensional multiplet has the fixed highest weight corresponding to alpha=1 in the family discussed in hep-th/9411241. There is no arbitrary continuous representation label introduced by this fusion.

| Step | Verdict and extent |
|---|---|
| 1. Duals and crossing | **Correct as an algebraic construction.** Explicit duals, four intertwiners, crossing scalars and double-dual shift are given below. |
| 2. Physical parameters | **Partly determined.** An explicit crossing-compatible gradation exists. The action-dependent spins and q(coupling) remain **unclear**. |
| 3. Scalars and bootstrap | **Correct for a specified formal branch.** There is a pole-free branch with no bound-state bootstrap to trigger, and the first attractive fusion and its crossed elementary exchange are verified. A physical minimal scalar and the full attractive bootstrap are **unclear**. |
| 4. Toda identification | **Unclear; not established.** Two candidate constructions were examined, and the straightforward identification does not follow. A constant highest-weight amplitude fails the classical sine-Gordon subsector test under the usual same-charge-kink identification. |

## Sources and what was read

1. Delius–Gould–Links–Zhang, *On Type-I Quantum Affine Superalgebras*, https://arxiv.org/abs/hep-th/9408006, v2. Read the coproduct, grading, intertwining equations and section 4.1, particularly (4.10)–(4.11). The paper explicitly gives the fundamental sl(2|1) example; it does not identify it as a Toda soliton S-matrix. Its matrix is used below only as a comparison with the independent solve.
2. Delius–Gould–Links–Zhang, *Solutions of the Yang–Baxter Equation with Extra Non-Additive Parameters II*, https://arxiv.org/abs/hep-th/9411241, v2, especially the highest weights, tensor products and Casimirs. These distinguish the defining triplet from the generic four-dimensional typical family.
3. Evans–Madsen, *Dynkin Diagrams and Integrable Models Based on Lie Superalgebras*, https://arxiv.org/abs/hep-th/9703065, corrected v2. Read the component action (3.4), section 5, Table 4 and Proposition 6. The table assigns a massless bosonic direction to the untwisted A(1,0) affine model. This was independently recovered from its roots below, rather than taken as evidence by itself.
4. Ridout–Teschner, *Integrability of a family of quantum field theories related to sigma models*, https://arxiv.org/abs/1102.5716, especially (2.3), (2.9), section 3.3, section 4's introductory discussion, and section 5.3. This is an important additional source missing from the earlier handover. It constructs a fermionic sl(2|1) Toda model, its screening algebra and a fundamental auxiliary R-matrix. It expressly distinguishes its construction from soliton-charge symmetry: “This does not seem to be directly related to the connections discussed in our paper.” (Section 4 introduction, referring to the Bernard–LeClair construction.)
5. Bassi–LeClair, https://arxiv.org/abs/hep-th/9911105, and Babichenko–Torrielli, https://arxiv.org/abs/1202.5405, were consulted for comparison. Their four-dimensional modules are not the defining triplet studied here. Their scattering conclusions have not been transferred to it.
6. Delius, https://arxiv.org/abs/hep-th/9503079, sections 2.2–2.4, was consulted for the role of charge spins and the distinction between an R-matrix and a physical S-matrix.

Local downloaded PDFs and searchable text are in the scratch directory. No literature priority or novelty claim is made for the algebraic construction below. A complete prior-art search for this exact normalization has not been performed.

## 1. Independent representation and intertwiner calculation

### Conventions

Use column vectors with parities p=(0,0,1), parity matrix J=diag(1,1,-1), and q=t^2. The symbol t in the scripts is an algebraic square root of q, not rapidity. The parameter x is a **ratio of inverse evaluation coordinates**, not automatically the physical rapidity variable.

For a homogeneous B of parity |B|, ordinary matrix storage realizes the graded tensor product as

\[
A\otimes_g B=(A J^{|B|})\otimes B.
\]

The graded permutation is P(v_i tensor v_j)=(-1)^(p_i p_j) v_j tensor v_i. We write R for the unpermuted even endomorphism and check R=P R for the scattering/braid ordering.

The defining representation used in the independent solve is

\[
\begin{array}{c|ccc}
i&0&1&2\\ \hline
h_i&\operatorname{diag}(-1,0,-1)&\operatorname{diag}(1,-1,0)&\operatorname{diag}(0,1,1)\\
e_i&z^{-1}E_{31}&E_{12}&E_{23}\\
f_i&-z E_{13}&E_{21}&E_{32}.
\end{array}
\]

Here E_ij are 3 by 3 matrix units. The coproduct is

\[
\Delta(e_i)=e_i\otimes_g q^{-h_i/2}+q^{h_i/2}\otimes_g e_i,
\qquad
\Delta(f_i)=f_i\otimes_g q^{-h_i/2}+q^{h_i/2}\otimes_g f_i.
\]

The inverse evaluation convention is explicit: e_0(z)=z^-1 E_31. Using e_0(z)=z E_31 instead gives the reciprocal x in the comparison with (4.11). This was an exploratory convention mismatch, not a correction to the source.

The program `algebra.py` makes a weight-preserving unknown 9 by 9 matrix and solves R Delta=Delta_op R for all six affine raising/lowering generators. It does not start from (4.11). It then compares P R against that equation.

Result: exact agreement. The independent calculation therefore confirms the earlier identification of (4.11), including the use of the graded permutation.

### Dual representation

Define the supertranspose by

\[
(E_{ij})^{st}=(-1)^{p_i(p_i+p_j)}E_{ji}.
\]

The antipode follows directly by applying the counit identity to the displayed coproduct:

\[
S(e_i)=-q^{-h_i/2}e_iq^{h_i/2},\qquad
S(f_i)=-q^{-h_i/2}f_iq^{h_i/2}.
\]

On the left dual, rho_bar(a)=rho(S(a))^st. Its matrices are

\[
\bar h_i=-h_i,\quad
(\bar e_0,\bar e_1,\bar e_2)=(z^{-1}E_{13},-q^{-1}E_{21},-E_{32}),\quad
(\bar f_0,\bar f_1,\bar f_2)=(z E_{31},-q E_{12},E_{23}).
\]

The triplet and this dual are inequivalent finite-algebra modules. A second dual satisfies

\[
D\rho^{**}_z(a)D^{-1}=\rho_{q^2z}(a),
\qquad D=\operatorname{diag}(q^2,1,-1).
\]

All generator equations in this identity are checked exactly in `crossing.py`. Keeping the odd minus sign in D is essential.

### Four normalized matrix parts

Let v denote V and b denote its left dual. The solver normalizes R_vv and R_bb to entry (11,11)=1. For R_bv and R_vb it normalizes entry (12,12)=1. These are choices of matrix normalization, not scalar S-matrix conclusions.

R_vv is the unpermuted form of (4.11); R_bb is its ordinary transpose in the displayed bases. For both mixed matrices, the six states |ij>, i != j, have diagonal entry 1. Their only nontrivial block acts on |11>, |22>, |33>:

\[
R_{bv}(x)\big|_0=
\frac1{q(x-1)}
\begin{pmatrix}
q^2x-1&q^2-1&q^2-1\\
x(q^2-1)&q^2x-1&q^2-1\\
-x(q^2-1)&-x(q^2-1)&x-q^2
\end{pmatrix},
\]

\[
R_{vb}(x)\big|_0=
\frac1{q(x-q^2)}
\begin{pmatrix}
q^2(x-1)&q^2x(q^2-1)&-q^2x(q^2-1)\\
q^2-1&q^2(x-1)&-x(q^2-1)\\
q^2-1&q^2(q^2-1)&x-q^4
\end{pmatrix}.
\]

Each of the four normalized matrices satisfies the exact braiding identity

\[
R_{ab}(x)P R_{ba}(x^{-1})P=I.
\]

The graded partial supertranspose must be applied to graded tensor coefficients, not naively to ordinary 3 by 3 blocks. The scratch implementation includes the conversion signs explicitly.

## 2. Scalars, crossing and an explicit rapidity gradation

### A compatible scalar choice exists

Define

\[
a(x)=\frac{q^2x-1}{q(x-1)}.
\]

The independent matrices obey

\[
(R_{vv}^{-1})^{st_1}=a(x)^{-1}R_{bv}(x),
\qquad
(R_{vb}^{-1})^{st_1}=a(x^{-1})R_{bb}(x).
\]

The other two equations include the double-dual similarity D and the shift x -> q^2 x. Their exact checks are in `crossing.py` and `candidate.py`.

One solution of all these identities is the following *dressed unpermuted* set:

\[
T_{vv}=R_{vv},\quad T_{bb}=R_{bb},\quad
T_{bv}=\frac{q(x-1)}{q^2x-1}R_{bv},\quad
T_{vb}=\frac{x-q^2}{q(x-1)}R_{vb}.
\]

There is no infinite product needed for this particular algebraic normalization. This is not an assertion that these are the physically minimal Toda factors. It is an explicit solution of the stated crossing/reciprocity problem.

### Rapidity convention and dual offset

Choose a positive real parameter kappa and set

\[
q=e^{i\pi\kappa},\qquad X=e^{\kappa\theta}.
\]

All tested kappas are irrational. Define the raw rapidity amplitudes by

\[
\mathcal T_{vv}(\theta)=T_{vv}(X),\quad
\mathcal T_{bb}(\theta)=T_{bb}(X),\quad
\mathcal T_{bv}(\theta)=T_{bv}(X/q),\quad
\mathcal T_{vb}(\theta)=T_{vb}(qX).
\]

The offsets implement the antiparticle evaluation convention. Omitting them changes crossing.

To give a common charge-conjugation formula, take

\[
G_v(\theta)=\operatorname{diag}(e^{-\kappa\theta},1,e^{\theta/2}),
\qquad G_b(\theta)=G_v(\theta)^{-1},
\qquad C=\operatorname{diag}(q^{-1},1,i),
\]

and, in the frame theta_2=0, write

\[
\widehat T_{ab}(\theta)=
(G_a(\theta)\otimes I)\mathcal T_{ab}(\theta)(G_a(\theta)^{-1}\otimes I).
\]

Total Cartan covariance makes this a function of the rapidity difference. With bar(v)=b and bar(b)=v,

\[
\widehat T_{\bar a b}(\theta+i\pi)
=(C^{-1}\otimes I)
\big(\widehat T_{ab}(\theta)^{-1}\big)^{st_1}
(C\otimes I).
\]

The actual exchange-ordered amplitudes are S_ab=P widehat T_ab. The formula above, together with braiding unitarity, states exactly which graded crossing convention is being used; it is not an ungraded template.

`physical_gradation.py` checks all four crossed channels and braiding unitarity at three irrational kappas. Maximum crossing residual: 1.62e-15. The corresponding generator spins are

\[
(s_0,s_1,s_2)=(1/2,-\kappa,-1/2).
\]

**This is a mathematical choice of gradation.** No action calculation in this investigation identifies these as spins of its conserved nonlocal charges. In particular, no equation kappa=kappa(beta) or q=q(beta) has been established. One must not substitute the bosonic book formula for q: isotropic odd roots prevent simply transplanting its coroot derivation.

For fixed q there are also logarithmic lift choices for kappa. Crossing alone does not select the weak-coupling lift or its physical spectrum.

### Scalar freedom remains

A common scalar h(theta) in the vv and bb channels, and h(i*pi-theta) in the mixed channels, preserves the construction if

\[
h(\theta)h(-\theta)=1,\qquad h(\theta+2\pi i)=h(\theta).
\]

For example,

\[
h_\alpha(\theta)=
\frac{\sinh\theta+i\sin\alpha}{\sinh\theta-i\sin\alpha}
\]

is a unitary, crossing-symmetric CDD factor. These identities are checked symbolically in `poles.py`. Neither this freedom nor a physical scalar phase is fixed by solving the intertwining equations.

## 3. Pole structure and first bootstrap calculations

### A pole-free formal branch

The only possible finite denominators of the rapidity matrix entries, before the harmless G similarities, are

\[
vv,bb:\quad X-q^2,\qquad
bv,vb:\quad qX-1.
\]

Thus possible poles theta=i*u occur at

\[
u_{vv}=2\pi+\frac{2\pi n}{\kappa},\qquad
u_{bv}=-\pi+\frac{2\pi n}{\kappa},\qquad n\in\mathbb Z.
\]

For **0<kappa<1**, none lies in 0<u<pi. The diagonal exponential similarities add no finite poles. On that branch, the triplet/dual candidate has no bound-state bootstrap to trigger. This is closure of the pole analysis of a formal factorized scattering candidate; it does not establish a local QFT, asymptotic completeness, or Toda correspondence.

### First attractive branch

For 1<kappa<2 there is one fundamental direct pole and its mixed crossed pole:

\[
u=2\pi(1-1/\kappa),\qquad u_{\rm crossed}=\pi-u.
\]

The exact vv residue has rank 4, since

\[
\operatorname*{Res}_{x=q^2}\check R_{vv}(x)=(1-q^4)P_4.
\]

For example, structural pole enumeration and separate 1/epsilon scaling confirm these poles at kappa=sqrt(2), sqrt(3). They are not inferred from a residue print threshold.

If this pole is a direct bound state, its mass is

\[
M_2=2M\cos(u/2).
\]

It is not an antifundamental: the rank is 4 rather than 3. The highest vector of the finite-algebra channel is

\[
|12\rangle-q|21\rangle,
\qquad (h_1,h_2)=(0,1).
\]

The raising generators kill this vector. On restriction to sl(2|1), this is the alpha=1 member of the typical highest-weight family (0,0|alpha). The check uses only the independently defined generators. **A fixed member of the parameter family appears naturally as a candidate bound state; the parameter is not free in this fusion.**

### Fused scattering was actually constructed

Use the rank-four embedding B with columns

\[
|12\rangle-q|21\rangle,\quad
|13\rangle-q|31\rangle,\quad
|23\rangle-q|32\rangle,\quad |33\rangle.
\]

The exchanged constituent evaluation ratios are -q^-1 and -q relative to the bound-state midpoint. The minus signs matter: exp(i*kappa*u/2)=-q on the chosen attractive lift.

`candidate.py` forms the full products R_13 R_23 on 27 dimensions and checks **zero leakage** out of B tensor V before restricting them to 12 dimensions. It constructs both the fused state against V and against its dual, with the mixed scalar factors included.

| Fused channel | Common denominator, up to nonzero q powers | Pole and residue rank |
|---|---|---|
| 4 x V | (X+q)(X+q^3) | X=-q: rank 12; X=-q^3: rank 4 |
| 4 x V* | (X+1)(q^2 X+1) | X=-1: rank 3; X=-q^-2: rank 12 |

The rank-three direct pole in 4 x V* lies at angle pi/kappa=pi-u/2. Its mass triangle is exactly the one for 4 + V* -> V, since cos(pi-u/2)=-M_2/(2M). The rank-twelve pole in 4 x V at u/2 is at the crossed angle of this elementary exchange. A crossed-channel residue need not have rank three as a matrix in the direct ordering; its rank is not evidence for twelve new elementary particles.

For 1<kappa<3/2, X=-q^3 includes a physical pole at 3u/2. If assigned to the next four-dimensional direct bound state, its mass is

\[
M_3^2=M_2^2+M^2+2M_2M\cos(3u/2),\qquad
M_3=M\frac{\sin(3u/2)}{\sin(u/2)}.
\]

The mass equality follows by substitution of M_2, not a new dynamical result. Poles from other logarithmic images and higher fusions must still be analysed for any proposed attractive theory. We have not determined their scalar cancellations, full spectrum, residue signs in a physical inner product, or Coleman–Thun interpretations. **This report does not claim a completed attractive bootstrap.**

## 4. Physical tests and attempted Toda identification

### Positive-metric unitarity is not established; the displayed candidate fails it

The identities above are *braiding unitarity*, S(theta) S(-theta)=I in the appropriate exchange spaces. They do not imply S dagger S=I on a positive Hilbert space.

At kappa=sqrt(2)/3 and real theta=0.39, the vv amplitude in the explicit crossing gradation has singular values approximately

```
1.33074415 1.10809263 1.09585495 1 1 1
0.91252953 0.90245163 0.75145925
```

The homogeneous convention also fails this test. Because three singular values are already 1 while others differ, an overall scalar cannot repair canonical positive-metric unitarity.

There is a stronger conditional obstruction. If the Cartan weights are physical self-adjoint charges, a positive one-particle metric is diagonal in the distinct weight basis. Its tensor square gives equal norms to |ij> and |ji>. In a real diagonal exponential gradation G(theta)=diag(exp(d_i theta)), each of the three two-state blocks of check R has the form

\[
\begin{pmatrix}
A e^{(d_j-d_i)\theta}&B\\
B&X A e^{(d_i-d_j)\theta}
\end{pmatrix}\quad(i<j),
\]

with A/B=(q-q^-1)/(1-X), purely imaginary for |q|=1 and real theta. Its off-diagonal unitarity equation requires d_i-d_j=-kappa/2 for **every** i<j. The conditions for (1,2), (2,3), (1,3) are incompatible when kappa is nonzero. Consequently no such factorized positive Cartan metric and real diagonal gradation make all three blocks unitary. This is not a theorem excluding non-unitary/Krein or restricted constructions, which are relevant to this book.

### Standard root-system super-Toda model

The independent script `toda_check.py` chooses

\[
\alpha_1=(\sqrt2,0),\quad
\alpha_0=(-1/\sqrt2,i/\sqrt2),\quad
\alpha_2=(-1/\sqrt2,-i/\sqrt2).
\]

Their Gram matrix is

\[
\begin{pmatrix}0&-1&1\\-1&2&-1\\1&-1&0\end{pmatrix},
\qquad \alpha_0+\alpha_1+\alpha_2=0.
\]

Take node 1 even and nodes 0,2 odd. Expanding the usual component construction directly gives, in root-field normalization with the overall coupling suppressed,

\[
V(u,v)=\mu_1 e^{\sqrt2 u}+\mu_0\mu_2 e^{-\sqrt2 u}.
\]

For positive mu_i, the vacuum Hessian is

\[
\operatorname{diag}(4\sqrt{\mu_0\mu_1\mu_2},0).
\]

Thus the second bosonic direction is massless at this classical level. The independently assembled Yukawa mass matrix obeys

\[
M_F^2=\mu_0\mu_2 e^{-\sqrt2u}I.
\]

The massless scalar cannot simply be dropped from the full theory: the Yukawa interactions depend on it. None of this proves that the quantum model has a massless asymptotic particle or forbids an infrared reorganization; it does defeat an unexamined identification with a purely massive triplet spectrum.

### Soliton and semiclassical-delay test

After an appropriate imaginary-coupling continuation and sign choice, the zero-fermion massive direction contains the ordinary sine-Gordon equation. Independently of a quantum-group identification, use

\[
V_{SG}(\chi)=\frac{m^2}{k^2}(1-\cos(k\chi)),\qquad
\chi_K(x)=\frac4k\arctan e^{mx}.
\]

The script verifies its equation of motion and integrates its energy to 8m/k^2. It also verifies a two-kink tau solution directly in the field equation. For eta_j=m*cosh(theta_j)*x-m*sinh(theta_j)*t,

\[
\tau=1+i e^{\eta_1}+i e^{\eta_2}-A e^{\eta_1+\eta_2},\qquad
A=\tanh^2((\theta_1-\theta_2)/2),\qquad k\chi=4\arg\tau.
\]

The fast kink's displacement is -log(A)/(m*cosh(theta_1)). At theta_1=0.57, theta_2=-0.81 and m=1 this is 0.881303075845, not zero. Field-equation residuals at three spacetime points are below 9.2e-16.

In the algebraic candidate without additional h, the highest-weight amplitude S_11 is identically 1, including in the crossing gradation. **If** that state is identified with a same-charge sine-Gordon kink, it predicts no scattering phase and cannot reproduce this delay. An action-dependent scalar phase and a verified state identification are therefore indispensable. This is a conditional mismatch of the simplest assignment, not a proof that no scalar factor can give a Toda theory.

The fermionic zero-mode quantization, complete kink multiplets and full two-soliton scattering of the candidate super-Toda theory have not been established here. Accordingly the assignment V to physical solitons remains **unclear**, rather than a fundamental-representation principle that can be assumed.

### A second, distinct candidate: fermionic Toda in Ridout–Teschner

Ridout–Teschner's action (2.3) has one real boson and a complex fermion, a four-fermion interaction, and exponential terms. It is not the preceding two-superfield component model. They derive a classical affine sl(2|1) Lax representation in section 3.3 and a quantum screening/monodromy algebra with q=exp(-i*pi*b^2) in section 5.3. Their (5.33)–(5.34) is the fundamental R-matrix in a principal gradation and auxiliary representation.

This supplies a concrete action-side lead, but the paper itself distinguishes this quantum algebra from the nonlocal charges governing asymptotic soliton scattering. In particular its screening q cannot simply be inserted into the present soliton ansatz. We have not constructed its soliton nonlocal charges, shown a degenerate fundamental triplet spectrum, or compared a corresponding semiclassical soliton phase. This identification also remains **unclear**.

## Reproducibility and outputs

Run from the repository root:

```
python3 -u referee_scratch/sl21_soliton/algebra.py
python3 -u referee_scratch/sl21_soliton/crossing.py
python3 -u referee_scratch/sl21_soliton/candidate.py
python3 -u referee_scratch/sl21_soliton/physical_gradation.py
python3 -u referee_scratch/sl21_soliton/poles.py
python3 -u referee_scratch/sl21_soliton/toda_check.py
```

The first command produces `matrices.pkl`, used by subsequent calculations. The candidate command writes `candidate.pkl` containing the four raw/rapidity matrices and both 12-dimensional first-fusion matrices. These are scratch artifacts, not permanent book caches. The text logs are `algebra.out`, `crossing.out`, `candidate.out`, `physical_gradation.out`, `poles.out`, and `toda_check.out`.

All six commands completed successfully. Their logs contain respectively 1, 9, 14, 1, 7 and 6 PASS lines. Several lines aggregate a family of exact identities or parameter tests; these are not counts of independent physical claims. The candidate log additionally prints all four fused residue ranks, rather than labelling those data lines PASS. The calculations take seconds to a few minutes per script on this machine; timings were not benchmarked.

Exact checks cover source comparison after an independent affine solve, all four braiding-unitarity equations, the double-dual generator equations, categorical crossing, residue ranks, the fusion subspace/highest weight and absence of leakage. Numerical checks cover all eight mixed-species Yang–Baxter assignments, the rapidity crossing convention, and pole scaling, at generic irrational parameters. The scalar ambiguity and classical mass/kink identities have symbolic checks. No existing book check is offered as evidence.

The scalar and mass dictionaries for an actual Toda QFT are explicitly absent. No general local-QFT existence, physical unitarity, attractive-spectrum completeness or new literature correction follows from these scripts.

## Integration

This section records how the suggestions in the conversation survived independent work; no book integration was performed.

| Earlier suggestion | Verified, modified or rejected |
|---|---|
| Use the defining fundamental R-matrix as the starting point | **Verified as an algebraic starting point.** Independent solve reproduces (4.11). The physical triplet assignment remains unproved. |
| Add antiparticles by duality | **Verified, with conventions made explicit.** Mixed normalization, evaluation offsets, partial-supertranspose signs and double-dual matrix are required. |
| Determine scalar factors from crossing/unitarity | **Modified.** A simple compatible set exists, but scalar CDD freedom remains and physical constraints are additional. |
| Determine q(beta) and the physical gradation | **Not achieved.** An explicit compatible gradation was constructed; an action-derived one was not. Importing the ordinary bosonic Toda q would be unsupported. |
| Interpret the four-dimensional pole as a possible bound state | **Verified conditionally.** Exact rank-four residue and first fusion. The highest weight fixes alpha=1; no free non-additive label appears. |
| Treat every fused full-rank residue as a new particle | **Rejected.** The first rank-twelve pole has a crossed elementary-exchange angle; direct and crossed rank criteria differ. |
| Treat braid unitarity as physical unitarity | **Rejected.** Explicit singular-value test and the conditional metric/gradation obstruction above. |
| Identify the result with standard sl(2|1) Toda | **Not established.** A massless classical direction, missing physical charge/multiplet calculation and missing semiclassical scalar phase must be addressed. |
| Use an auxiliary Toda R-matrix as proof of soliton symmetry | **Rejected.** Ridout–Teschner explicitly distinguish the two constructions. |

## What would resolve the remaining physical question

Choose one action and its real form first. Construct its conserved charges on the infinite line and determine their Lorentz spins, coproduct and action on quantized kink zero modes. Those data must decide whether the asymptotic symmetry is the present quantum affine algebra and whether the kink module is V, V*, a typical module, or a product with an independent supersymmetry factor.

Only then choose the logarithmic lift kappa(beta), compute the two-kink semiclassical phase, fix the nontrivial scalar/CDD ambiguity against it, and pursue the physical attractive bootstrap. A positive-metric interpretation would additionally have to overcome the obstruction identified above, or specify a restriction or a different physical inner product. A non-unitary imaginary-coupling interpretation still needs its own reality/analyticity and state-space prescription.

These are specific missing calculations, not a claim that the answer is impossible. The algebraic construction in this report makes them concrete and testable.
