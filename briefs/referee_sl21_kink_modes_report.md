# Root-system sl(2|1) Toda kinks: modes, backreaction and representation tests

Date: 2026-10-23. Model: Codex (GPT-6).

Book baseline: `745586e6a6f7bca5ed91a1e9cf89d6ef5df387f0`.

Task: follow the proposed dynamical identification test for the fundamental R-matrix. Calculate the kink's fermionic modes and charge assignments, and determine whether these establish the required soliton representation.

This is a research/referee report. No book text or book check script was edited. Independent scripts and their output are under `referee_scratch/sl21_kink_modes/`. No author script was imported. The earlier report `briefs/referee_sl21_soliton_report.md` is left unchanged. No claim of novelty is made.

## Outcome and verdicts

The standard root-system action has exactly two localized holomorphic fermion zero modes on its sine-Gordon kink. Their projected kinetic pairing is nondegenerate. Formal collective-coordinate quantization therefore gives a two-state space per orientation. This is not yet a physical quantum kink multiplet: the fields and modes are complex, and their full nonlinear backreaction changes the massless scalar's boundary value.

The backreaction can be solved exactly with two Grassmann generators. In this static family the massless boundary shift is proportional to an exact axial Noether charge. Thus the massless direction is specifically a boundary-sector issue, rather than a general obstruction to solitons.

Nothing calculated here identifies the defining triplet as the physical soliton representation. An isolated two-state space cannot carry a nontrivial type-one sl(2|1) module. A compatible four-state module can be built from the proposed kink/antikink charge pattern, but its generators have not been constructed from the field theory. This distinction is essential.

| Item | Verdict |
|---|---|
| Classical kink and two localized linear fermion modes | **Correct**, for the specified complex root-system action. Exact derivation and numerical ODE check below. |
| Formal two-state Clifford quantization per orientation | **Correct as a projected calculation**; its interpretation as an infrared-complete quantum multiplet is **unclear**. |
| Exact Grassmann backreaction and charge/boundary relation | **Correct**, for the two-generator static family and source action normalization. |
| Identification with the defining fundamental triplet | **Unclear; unsupported by this calculation.** Conditional dimension and charge obstructions are given below. |
| Four-state finite quantum-superalgebra representation | **Correct as an algebraic compatibility example**. Physical realization and affine extension are **unclear**. |
| Action-derived quantum-affine nonlocal charges and Toda S-matrix | **Unclear; not constructed**. |

## 1. Source, action and conventions

The source is Evans–Madsen, [hep-th/9703065](https://arxiv.org/abs/hep-th/9703065), corrected v2, equations (3.4)–(3.5). The source explicitly says: “the Toda model really does depend on the particular simple root system”. This matters here: the following action is not identified with the different one-boson fermionic model of Ridout–Teschner, [1102.5716](https://arxiv.org/abs/1102.5716), equation (2.3).

Choose an ordinary complex bilinear dot product, fields phi=(u,v), and roots

\[
\alpha_0=(-1/\sqrt2,i/\sqrt2),\quad
\alpha_1=(\sqrt2,0),\quad
\alpha_2=(-1/\sqrt2,-i/\sqrt2).
\]

Nodes 0 and 2 are odd, node 1 even. Their Gram matrix is

\[
\begin{pmatrix}0&-1&1\\-1&2&-1\\1&-1&0\end{pmatrix},
\qquad \alpha_0+\alpha_1+\alpha_2=0.
\]

Set mu_0=mu_2=mu>0 and mu_1=mu^2. With an overall coupling factor 1/g^2, the source component action becomes

\[
g^2\mathcal L=
\tfrac12\partial_+\phi\cdot\partial_-\phi-V
+\tfrac i2(\psi_+\cdot\partial_-\psi_+
+\psi_-\cdot\partial_+\psi_-)
+i\psi_+^t M\psi_-,
\]

\[
V=\mu^2(e^{\sqrt2u}+e^{-\sqrt2u}),\qquad
M=\mu e^{-u/\sqrt2}
\begin{pmatrix}
\cos(v/\sqrt2)&\sin(v/\sqrt2)\\
\sin(v/\sqrt2)&-\cos(v/\sqrt2)
\end{pmatrix}.
\]

Use partial_+=partial_t+partial_x and partial_-=partial_t-partial_x. In particular, the static fermion equations are

\[
\psi_+'=M\psi_-,\qquad \psi_-'=M\psi_+.
\]

The calculation is complex/holomorphic. No positive real contour, adjoint operation or positive Hilbert-space inner product for the complete interacting model has been supplied. In particular, an absolute L2 norm used to test decay must not be mistaken for a proof of physical unitarity.

## 2. Classical saddle and complete linear zero-mode calculation

Put u=i chi/sqrt(2), v=0. The static equation reduces to

\[
\chi''=4\mu^2\sin\chi.
\]

Write y=exp(2 mu x). A kink is chi_K=4 arctan(y). On it,

\[
w=e^{-i\chi_K/2}=\frac{1-iy}{1+iy},\quad
r=\mu w,\quad M_K=r\,\operatorname{diag}(1,-1).
\]

Two independent zero modes are

\[
\psi_+=f(a,b)^t,\qquad
\psi_-=f(a,-b)^t,\qquad
f=\frac{\sqrt y}{1+iy},\qquad f'=rf.
\]

Here a and b are odd collective coordinates. Diagonalizing the two species and the plus/minus chirality combinations gives four first-order scalar equations. The other two solutions are proportional to 1/f, with the opposite chirality ratios. They grow at both ends. This proves completeness of the two localized solutions of the linear static problem.

Their decay and bilinear pairing are

\[
|f|^2=\frac1{2\cosh(2\mu x)},\qquad
\int_{\mathbb R}|f|^2dx=\frac\pi{4\mu},\qquad
I_K:=\int_{\mathbb R}f^2dx=-\frac i{2\mu}.
\]

The full four-component mode matrix U obeys U^t U=2 f^2 I_2. Thus its projected kinetic form is nondegenerate, even though it is not a positive Hermitian norm.

For the antikink chi_A=4 arctan(1/y), use

\[
w_A=\frac{y-i}{y+i},\quad r_A=\mu w_A,\quad
f_A=\frac{\sqrt y}{1-iy},\quad f_A'=-r_A f_A,
\]

\[
\psi_+=f_A(a,b)^t,\qquad
\psi_-=f_A(-a,b)^t,\qquad
I_A=+\frac i{2\mu}.
\]

The bosonic kink fluctuation operator in the chi direction is

\[
-\partial_x^2+4\mu^2[1-2\operatorname{sech}^2(2\mu x)]
=A^\dagger A,\quad A=\partial_x+2\mu\tanh(2\mu x).
\]

Its L2 zero mode is the translation mode. The free v operator is -partial_x^2 and has no nonzero L2 static zero mode on the line. The massless continuum remains present.

Independent evidence: `kink_modes.py` substitutes the roots and profiles symbolically, constructs the growing partners, evaluates both pairings, and integrates the original four-component ODE from unit-normalized central data at mu=sqrt(2). The maximum absolute profile error on |x|<=8/mu was `1.1541152632823952e-10`. All checks passed.

## 3. Exact nonlinear backreaction

Linear zero modes need not extend to a family with the same boundary conditions. Here their extension is explicit. Take a^2=b^2=0 and ab=-ba. Keep u at its kink value and set

\[
v=ab\,v_2(x),\qquad
v_2(x)=\frac{iy}{2\sqrt2\mu(1+iy)}.
\]

Together with the fermion profiles above, this satisfies the full boson and fermion equations, not just their linearization. All higher Grassmann terms vanish. The v equation reduces to

\[
v_2''=i\sqrt2\,r f^2,\qquad
v_2'=\frac i{\sqrt2}f^2,
\qquad
v_2(+\infty)-v_2(-\infty)=\frac1{2\sqrt2\mu}.
\]

For the antikink replace 1+iy by 1-iy. Its shift is -1/(2 sqrt(2) mu), and v_2'=i f_A^2/sqrt(2).

The independent script `backreaction.py` implements the exterior algebra explicitly, constructs every exponential from the three roots, and substitutes the corrected fields into every component of the full source equations. All Grassmann coefficients vanish. It also checks both boundary shifts and the charge identity below.

If the full supernumber values v(+infinity)=v(-infinity)=0 are imposed, this two-generator static family requires ab=0. Adding a constant changes both endpoints together; adding a linear homogeneous solution violates the finite-energy boundary condition. Allowing the nilpotent endpoint shift gives an exact localized extension.

This is a classical boundary-sector statement. It does not prove that quantum fermion occupation is forbidden or that a quantum doublet is absent. Quantization, continuum polarization and boundary dressing have not been computed.

## 4. Exact action symmetries and their limitations

Let J be the matrix with rows (0,-1) and (1,0). Since MJ+JM=0, the action has an axial rotation

\[
\delta\psi_+=\epsilon J\psi_+,\qquad
\delta\psi_-=-\epsilon J\psi_-,\qquad \delta\phi=0.
\]

In source normalization its conserved current is

\[
j_+=-i\psi_{+u}\psi_{+v},\qquad
j_-=+i\psi_{-u}\psi_{-v},\qquad
\partial_-j_++\partial_+j_-=0.
\]

Define Q_F=g^{-2} integral(j_++j_-)dx. On either exact static family, the v equation gives the concrete relation

\[
\boxed{\Delta v=-\frac{g^2}{2\sqrt2}Q_F.}
\]

Indeed j_++j_-=-2i ab f^2, with the appropriate orientation's f. This equation refers to the classical collective-coordinate charge; a quantum sea contribution or normal-ordering constant is not included. Reversing the overall action sign reverses its Noether charge convention.

There is also a combined symmetry v -> v+epsilon and psi_± -> exp[epsilon J/(2 sqrt(2))] psi_±. It follows from partial_v M=[J,M]/(2 sqrt(2)). This symmetry moves the boundary vacuum, so it is not automatically an operator within a fixed-vacuum scattering space.

The standard local W-derived supersymmetry does not supply the missing odd charges. Although

\[
W=2\mu e^{-u/\sqrt2}\cos(v/\sqrt2),\qquad M=\operatorname{Hess}W,
\]

the potential contains the additional even-root exponential. Explicitly,

\[
\nabla V-(\operatorname{Hess}W)\nabla W
=(\sqrt2\mu^2e^{\sqrt2u},0)^t\ne0.
\]

This is the nonzero divergence coefficient in the usual W-supercurrent test. It excludes that construction of local supercharges. It says nothing by itself against other local constructions or genuinely nonlocal charges.

Independent evidence: `charges.py` verifies the full current divergence with four independent Grassmann generators, both matrix symmetry identities, and the supercurrent mismatch. All checks passed.

## 5. What representation follows?

### The projected result

The two nondegenerate odd coordinates formally quantize to Cl_2(C). Its irreducible graded space is 1|1. Pauli matrices give an explicit realization, with centred axial charge F=-i gamma_1 gamma_2/2 and eigenvalues +/-1/2. The charge difference is one. The absolute physical charge can include a vacuum/continuum contribution, which has not been calculated.

This calculation produces a candidate doublet per orientation. It is not a derivation of asymptotic states of the full massless theory.

### Conditional obstacles to a defining triplet

For generic q, an isolated 1|1 type-one module cannot carry a nontrivial finite U_q(sl(2|1)) action. Its even sl2 generators act separately on one-dimensional parity sectors, hence vanish in a finite-dimensional type-one module. Then h1=0, and [h1,e2]=-e2 and [h1,f2]=f2 force the odd generators to vanish. This is a representation argument, not a restriction on symmetries that mix sectors or additional states.

A separate charge test uses the defining triplet weights (h1,h2)=(1,0),(-1,1),(0,1). If T=A h1+B h2, its eigenvalues are A,B-A,B. They cannot all belong to {+1,-1}: A and B must be +/-1, whose difference is 0 or +/-2. Thus a triplet made only of unit-charge kinks and antikinks is incompatible with this Cartan identification without a central offset. This does not exclude triplets containing other topological charges, neutral states, additional quantum states, or another identification of T.

### A four-state compatibility example

There is an explicit generic-q finite module compatible with two centred doublets, one at T=+1 and one at T=-1. In a basis with parities (0,0,1,1), put

\[
H_1=\operatorname{diag}(0,0,-1,1),\qquad
H_2=\operatorname{diag}(-\tfrac12,\tfrac12,\tfrac12,-\tfrac12),
\]

\[
E_1=E_{43},\quad F_1=E_{34},\quad
E_2=-E_{14}+E_{32},\quad
F_2=[\tfrac12]_q(E_{23}+E_{41}),
\]

where E_ij are four-dimensional matrix units and [a]_q=(q^a-q^{-a})/(q-q^{-1}). Set T=2H2 and F=H1+H2. The weights are

| State | T | F | Parity |
|---|---:|---:|---:|
| 1 | -1 | -1/2 | even |
| 2 | +1 | +1/2 | even |
| 3 | +1 | -1/2 | odd |
| 4 | -1 | +1/2 | odd |

The highest vector is state 1, with (H1,H2)=(0,-1/2). In the convention labelling typical modules by this second highest weight, alpha=-1/2 is fixed. This is distinct from the alpha=1 four-state fusion module in the previous report.

`charges.py` checks all finite Cartan relations, supercommutators, odd nilpotence and Serre relations exactly, using q=t^2. It checks the selection rules

\[
E_1:(\Delta T,\Delta F)=(-2,+1),\qquad
E_2:(\Delta T,\Delta F)=(0,-1).
\]

This construction assumes centred charges and a relative parity convention between orientations. Neither is fixed by the mode calculation. More importantly, the topological-charge-changing field operator E1 has not been constructed. These matrices are a compatibility example, not evidence that the quantum field theory realizes them. No affine generator, physical coproduct, q(g), charge spin or S-matrix is inferred from them.

## 6. What the calculation leaves to establish

The next decisive calculation is an action-derived conserved nonlocal current, with its boundary terms and its action on dressed kink states. In particular, a generator with the E1 selection rule above would have to change the kink's topological sector. Linear zero-mode ladder operators within a frozen background cannot do this.

Before using any representation to choose an R-matrix one must also specify a quantum contour/reality condition and the asymptotic boundary sectors, including the massless field. The axial charge offsets require a continuum spectral-asymmetry or equivalent quantum calculation. A physical nonlocal-charge calculation must then determine the generators' Lorentz spins, q as a function of the action coupling, and the coproduct. Only then can a particular affine intertwiner and scalar factor be tested against semiclassical scattering.

There is consequently no justified switch from the defining-triplet proposal to a four-state S-matrix at this stage. The present calculation narrows the dynamical question and supplies an exact test for any proposed boundary treatment or charge construction.

## Independent calculation record

Scratch files, each written independently for this task:

- `referee_scratch/sl21_kink_modes/kink_modes.py`, output `kink_modes.out`.
- `referee_scratch/sl21_kink_modes/backreaction.py`, output `backreaction.out`.
- `referee_scratch/sl21_kink_modes/charges.py`, output `charges.out`.

Run each with `python3 -u` from the repository root. Dependencies are SymPy, NumPy, SciPy and mpmath. Symbolic calculations use generic parameters, not a decimal root-of-unity q. The scripts print individual PASS lines; failures raise assertions.

## Integration

No book integration was performed. Independently verified: the specified root-system action, complete localized linear-mode space, nondegenerate projected pairing, exact Grassmann backreaction, axial current, boundary/charge relation and finite-module compatibility example. Modified interpretation: two linear modes cannot yet be advertised as a physical two-state soliton multiplet with fixed massless endpoints. Rejected as unsupported: identification of either the fundamental triplet or the compatible quartet as the Toda soliton representation. Unresolved: real quantum model, infrared dressing, absolute charges, nonlocal affine generators, coupling dictionary and physical soliton S-matrix. These findings require independent review before being promoted into book claims.
