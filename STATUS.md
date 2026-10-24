# Current technical status

Snapshot: 24 October 2026. Latest evidence is the two-barred construction, coordinator report and independent audit. These are internally checked formal research results; they have not become an externally refereed theorem about a quantum theory.

## Overall objective and present workstream

The original objective is to construct an S-matrix from the three-dimensional fundamental $U_q(sl(2|1))$ R-matrix in [hep-th/9408006](https://arxiv.org/abs/hep-th/9408006), and test whether it describes solitons in an affine supersymmetric Toda field theory. The current contact/quotient investigation is part of the attempted field-theory identification.

The [earlier scattering report](briefs/referee_sl21_soliton_report.md) records fundamental/dual intertwiners, crossing-compatible scalar choices and gradation, Yang–Baxter checks, and initial fusion results. It also records unresolved scalar and physical-parameter choices and a positive-metric unitarity obstruction in the tested gradation. The [kink-mode report](briefs/referee_sl21_kink_modes_report.md) does not establish a physical fundamental triplet. These reports retain their original qualifications; this summary does not independently re-referee them.

Completion of the present charge-algebra work would still leave the physical soliton representation, action-derived parameter dictionary, scalar factor, spectrum/complete bootstrap, and comparison with semiclassical soliton scattering. The algebraic candidate has not been identified as a Toda soliton S-matrix.

## Results to retain, with their exact limits

1. The common canonical Gaussian product and both chiral Serre quartics were checked through mass order two.
2. The once-applied interaction map, with complete first source/test jets, was checked on both Serre quartets with one extra barred generator through mass order two.
3. For the labelled `BB barBbarB` word, the same full local forest gives grouping/collision-order agreement as a common meromorphic family. Its physical restriction is regular at generic `2<t<3`, away from the possible new unlike-source resonance `t_*=(7+sqrt(17))/4`. The equal-source sector is regular throughout the open patch.
4. The two-barred result is a same-word grouping test. It is NOT a proof of Serre descent with two barred generators or unrestricted associativity. The map is applied once after Gaussian composition; independently dressing every binary block would be a different construction.

## Essential conventions

`a0=(ib/sqrt(2),b/sqrt(2),i)`, `a2=(ib/sqrt(2),−b/sqrt(2),−i)`;
`t=b²`, `T=a0+a2`, `d=−T/t`, `rho=2/t`, `beta=t−1`, `gamma=t/(1−t)`, `k=d+T`.

Covariance is `X(z)X(w)~−log(z−w)`; scalar products are bilinear. Both currents are even:

`B=:C exp(d.X):`, `C=a0.partial X`, `A=:exp(−d.X):`.

`Yi=:exp(ai.(X+Xbar)):`, `V1=:exp(−T.(X+Xbar)):`.

The tails are `Ki=:exp((d+ai).X+ai.Xbar):` and `K_A=:exp(−k.X−T.Xbar):`. The AA endpoint is the DIFFERENT field `F_A^end=:exp(−k.(X+Xbar)):`. The BB endpoint is `U=:exp((d+a2).(X+Xbar)):`.

Input currents are `B−t lambda2 K2`, `A−gamma lambda1 K_A`, and their barred counterparts. Couplings lambda0,lambda2 have mass order one; lambda1 has mass order two. Rotational unary repairs are `+lambda0 K0`, `+beta lambda2 K2`, `−lambda1 K_A`. In the spatial representation of the physical B input, the combined original-slot replacements are therefore **+K0 and −K2**. An original derivative slot can be consumed only once.

Preserve the whole-word cocycle
`epsilon(n,m)=(-1)^((nL+nR)(m0+m2)) exp[−i pi(nL mR−nR mL)/t]`.
Use `nu=exp(−i pi/t)`, `s=nu²`.

Regulators stay attached to original constituents after fusion. In the two-bar test use independent gL,gR and h, then restrict `gL=gR=rho`, `h=t` only after extending the complete family. The unlike source-pair factor is `|w0−w2|^(2h−2)`.

## Contacts that must not be simplified away

For a mixed BB/source face, with `q0=1`, `q2=beta`, the ordered raw moment has prefactor `−i pi sigma nu lambda_i` times

`L_i[psi]=qi²(partial_w+partial_wbar)psi−qi(Ceff+Cbar_eff)psi`.

Differentiate the complete coefficient before restriction, including sources, tests, spectator contractions and derivative-field Wick pairings. The finite balanced correction is the negative of this full moment. The AA analogue is `L_A[psi]=(partial_w+partial_wbar)psi`.

For two disjoint mixed faces the full Gaussian generating-jet operators commute. The double crossing is generally nonzero. If `Delta12=R11−R10−R01+R00`, then a global two-extreme average exceeds independent local balancing by **Delta12/4**. This correction is fixed on separated faces, not a new adjustable full-diagonal contact. Do not replace the local forest by a uniform average of all words or by iterated Jordan products.

Shared-source/current faces overlap; never multiply their delta distributions as if they were disjoint. The source-pair P counterterm, L improvement and all full Taylor descendants must be included once, never duplicated as additional action insertions.

## The next unresolved scalar quantity

For the full unlike-source `BB barBbarB` collision,

`omega=gL+gR+2h−7`; physically `omega=4/t+2t−7`.

At `t_*` this is zero and all proper two-source faces are regular. A surviving residue would be a primitive zero-jet quartic term proportional to

`lambda0 lambda2 :exp[(2d+T).(X+Xbar)]:`

times three real deltas identifying the four current positions and the fixed word cocycle.

The complete-current BULK calculation gives `4pi⁴` for its radial pole coefficient at `t_*`, before action normalization/cocycle. It omits height-supported unary, balanced and input-tail contributions. Their sum can alter or cancel the residue. **Do not call the bulk number a full anomaly.** The two groupings already agree meromorphically, including any common principal part; a surviving pole would not itself be an associator.

Raw physical-height convergence is unproved. Multi-barred Serre descent, arbitrary words, all orders, a physical domain/state, and identification or faithfulness of the symmetry algebra remain open.
