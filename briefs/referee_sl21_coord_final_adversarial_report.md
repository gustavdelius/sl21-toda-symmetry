# Final adversarial audit of the common source-exponent construction

Date: 2026-10-24. Model: Codex (GPT-6). Reviewed commit: `745586e6a6f7bca5ed91a1e9cf89d6ef5df387f0`.

The untracked inputs are identified by `referee_scratch/sl21_symmetry_coordinator_20261024/final_audit/input_hashes.json`. Fresh independent calculations are `audit.py` and `audit.out` in that directory. No earlier research script was imported or used as evidence. No subagents, book edits, existing-report edits, git mutations, worktrees, commits, pushes, or renders were used. The draft coordinator report was also read to check the scope of its conclusion; it was still being completed and is not an immutable audited input.

## Conclusion

**The common-family mass-squared projected identity survives this adversarial audit as a formal distributional identity in its specified analytic prescription. No global quantum algebra follows.** I found no pole face identically coinciding with `h=t`, and the explicit nonuniform overlap diagnosed in the opposing report does not furnish a counterexample to this continuation. It does refute a direct use of the separated-source cutoff remainder in the negative-t strip.

The short presentation of H8 omits a useful justification: the product rule must be applied with the moving masks still present, including derivatives of those masks. Section 3 below supplies that step and explains why neither a new mask-intersection contact nor the finite binary matching invalidates the factor `h−t`. This conclusion uses the fixed first-order Cauchy/mixed-face prescription. It is not an assertion about a cutoff product with different finite parts, arbitrary-test associativity, a regulator-independent operator product, or a quantum completion.

| Substantive claim | Verdict | Exact scope |
|---|---|---|
| H2 double-source kernel and H5 residues | Correct | Independently reconstructed from the three-vectors, including the A derivative contraction |
| Local cutoff coefficient is t/h, and H7 has an h−t factor | Correct | Actual local jet polynomial and circular subtraction; separated-face formula alone is not a uniform estimate |
| Common starting domain and joint continuation | Correct for the stated finite-order compact-test family | First-order Cauchy/mixed faces are fixed first; new source-containing faces are treated together |
| No same-diagonal pole obstructs generic h=t restriction | Correct | New face exponents are affine in h and rho; all old parameter-independent singular faces are already prescribed |
| H8 and the projected order-lambda0-lambda2 identity | Correct in this formal common-family prescription | Product-rule argument with full masks, mixed contacts, and common continuation; not an actual Toda Hilbert-space identity |
| Lower mixed family agrees with the positive-t continuation | Correct | Equality on the positive diagonal of the lower family, followed by one-variable uniqueness |
| The opposing nonuniform coefficient obstructs this common continuation | Wrong if inferred | Its pointwise divergence is real; its union-face degree is positive in the starting domain |
| Further ordinary local Gaussian anomalies excluded | Correct conditionally | Integer proof below; nonlocal/line/end contacts and nonanalytic terms excluded from its hypotheses |
| Naive real positive-coupling reference-boson cycle converges | Wrong | Open-region exponential runaway, including finite lattice regulator |
| Scalar reflection data determine interfaces or endpoint weak limits | Wrong if inferred | Explicit linear-algebra and weak-limit counterexamples below |
| Regional anchor yields a formal string coproduct | Correct | Requires regional operator spaces to become a coproduct; string is not yet identified with Cartan |
| Assumed relations convert to level-zero quantum affine sl2 | Correct conditionally | Explicit algebraic conversion; no faithfulness or actual Toda action established |

## 1. Independent reconstruction of the decisive local kernel

Use the explicit three-vectors in the plan, not a precomputed Gram matrix. Their bilinear contractions give

`a0²=a2²=−1`, `a0.a2=1−t`, `a0.d=a2.d=1`, `d²=−2/t`.

The derivative current contracts with a vertex of charge p by `−a0.p/(z−x)`. Thus its contributions from A, Y0 and Y2 are respectively `1/(z−a)`, `1/(z−v)`, and `(t−1)/(z−w)`. The double-source scalar factor is consequently

`(a−z)^−rho (a−v)(a−w)|v−w|^(2h−2) / [(z−v)(z−w)(ubar−vbar)(ubar−wbar)]`,

and its derivative polynomial is exactly H2. The source action contributes `lambda0 lambda2/pi²`, since the two cross terms cancel the factorial. The cocycle can also be checked directly by adding the four-component charge labels in order: A and B cancel their chiral labels before the antichiral factor, and barB crossing the two Yukawa labels contributes two minus signs. Their product is +1.

Let F0 and F2 be the collapsed kernels, including the A factor, the antichiral polynomial, and the full normal product. With derivative backgrounds suppressed only for displaying scalar ratios, direct Laurent expansion at z=v yields

`S0/F0 = (rho−1)/(a−v)+(t−2)/(v−w)`.

Direct differentiation of F0 yields

`(partial_v F0)/F0 = (rho−1)/(a−v)+(h−2)/(v−w)`.

Their difference is `−(h−t)/(v−w)`. This calculation includes the often-missed derivative contraction `1/(z−a)`.

At z=w, the double coefficient is `(t−1)F2` and the actual simple coefficient divided by F2 is

`[(t−1)rho−1]/(a−w)+(2−t)/(w−v)`.

Using rho=2/t, adding `(partial_w F2)/F2` leaves `(h−t)/(w−v)`. The field jets match because

`a0+d=d+a0`, `a0+(t−1)d=−(d+a2)`.

Hence H5 is an identity of the actual first jets, not a scalar-homogeneity guess.

For local fields c=a0.partial X and e=a2.partial X, expanding the two remaining source vertices gives the second-jet difference

`N=[c'−e'+c²−e²]E/2`.

The P insertion gives `R_P=[c(c+e)−t(c'+e'+(c+e)²)/2]E`. Since `d.partial X=−(c+e)/t`, direct differentiation of L gives

`partial L=[c'+(1−1/t)c(c+e)]E`.

It follows by subtraction that `N−R_P/t=partial L`. Radial integration of the retained angular mode is `(R^(2h)−epsilon^(2h))/h`; its boundary term therefore has coefficient t/h. The full local subtraction polynomial is

`N−R_P/h−(t/h)partial L = (h−t)N/h`.

This proves the local part of H7 independently. It does not use the Gaussian-test integral in the kernel worker's script.

## 2. The overlap and the common domain

An independent cluster enumeration from the actual pair contractions gives all eight new source-containing degrees:

| Face | Relative degree |
|---|---:|
| 02 | 2h |
| A02 | 2h+4 |
| B02, D02 | 2h−1 |
| AB02 | 2h+2−rho |
| AD02 | 2h+2 |
| BD02 | 2h−3 |
| ABD02 | 2h−rho |

Here D denotes barB. The number of relative real measures is the number of current coordinates plus twice the number of sources, minus one common translation for a current-containing cluster or two for a source-only cluster. Every derivative current with a nonzero internal contraction subtracts one degree. The script reconstructs these degrees directly rather than inserting the table as input.

The old B/source and D/source faces remain Cauchy distributions, and the old BDi mixed faces require their fixed first-jet contact. They do not become ordinary absolutely convergent integrals by increasing h. Once those prescriptions are included, all genuinely new faces have positive degree on `Re h>2, Re rho<0`. In particular the most singular proper source-containing face with A separated is BD02, of degree `2h−3`.

The opposing report's local diagnostic is correct. Its scalar series follows from expanding the elementary rational current factor; after removing the P mode the geometric derivative sum is

`sum_(n>=1) (n+t)q^n = q[(1+t)−tq]/(1−q)²`.

At the physical h=t this gives the reported nonzero omitted disk proportional to `epsilon^(2t−6)` when both current separations scale with epsilon. But the integrations over the two current positions and the source-centre transverse coordinate supply three additional powers, giving the actual union-face degree `2t−3`, or `2h−3` in the enlarged family. The pointwise diagnostic cannot be discarded at negative t. Equally, it cannot be used to refute the enlarged convergent family: at h>2 that union-face contribution vanishes. Continuation thereafter treats its whole finite part together with the other faces.

For a general h the numerical coefficient in the mode sum changes because the radial denominator is n+h rather than n+t; the degree argument does not require the physical-diagonal coefficient to stay unchanged.

## 3. Moving masks and the finite contact in H8

Set eta=h−t and `R=(q qbar)^eta`, q=v−w, with the real radial logarithm. On the complement of the source diagonal,

`partial R/R=eta/q`,
`partial²R/R=eta(eta−1)/q²`,
`partial partialbar R/R=eta²/(q qbar)`.

Every nonzero derivative of R contains eta, including the second derivatives that occur in intersecting current contractions. A first-derivative argument alone would be insufficient.

Now keep the entire finite-cutoff mask M in the calculation: source-pair exclusion, the two inner-current circular exclusions, strip edges, and smooth infrared/test factors. For any kernel/distribution Q on an already prescribed proper face,

`partial²(R M Q)−R partial²(M Q) = 2(partial R)partial(M Q)+(partial²R)M Q`.

The analogous mixed derivative has three terms, each containing a nonzero derivative of R. This identity includes derivatives of M; it does not first replace M by one. It remains valid using smooth approximations to the masks and then their prescribed distributional limits. In particular terms supported on moving circle boundaries are included. Terms where two masks meet lie on the corresponding union face. Their remaining scale is the same cluster degree counted above; they vanish in the starting domain after the fixed proper faces are dealt with. The inequality c>2 keeps the relevant strip edges away from the inner circle.

Apply this product rule to the actual source-collapse Laurent identities from section 1. Terms with no derivative on R reproduce the first-order Ward/tail cancellations, now multiplied by the spectator factor. On a proper mixed face they reproduce the already fixed raw coefficient t and its spectator dressing. Derivatives on R give eta times a calculated distribution. A/source crossings add no separate residue, and R is independent of the A coordinate. The X14 permutations leave R unchanged and leave the endpoint branch phases unchanged.

The finite mixed addition is equally important: its dressed kernel is `t(t−1)lambda2 R|_(source2=x) U(x)Y0(v)`. Differentiating this supported contact gives its usual endpoint term plus derivatives of R, again with eta. Assigning the matching only after the continuation would not obey this argument.

Finally compare the h-dependent subtraction coefficients with transporting the physical coefficients while replacing their cutoff power by epsilon^(2h). Their differences are all divisible by eta:

`−1/h+1/t=eta/(ht)`,
`t/h−1=−eta/h`,
`−1/(h−1)+1/(t−1)=eta/[(h−1)(t−1)]`.

These account for the P insertion, current improvement, and both endpoint descendants. Remaining cutoff jets are supported on the new faces already controlled in section 2. Thus there is no extra eta-independent moving-mask or union-face term in the initial domain. The surviving eta-independent mixed term has exactly the prescribed coefficient t². Together with the unchanged endpoint phases this gives H8 for the restricted four-word tests.

This is the missing detail behind the compressed paragraph at `briefs/referee_sl21_coord_stage1_kernels_report.md:191`. Suitable additional wording is:

> Keep all cutoff masks inside the differential product rule. The commutator of every source derivative with the common factor R_h contains h−t, including the second and mixed derivatives on intersecting masked faces. The terms with no derivative on R_h are the fixed first-order contacts multiplied by the source spectator; the differences in all continued local counterterm coefficients also contain h−t. Every remaining mask-intersection boundary lies on a new source-containing union face of positive degree in the initial domain. Therefore no additional h−t-independent contact survives there.

No replacement of the claimed formula is needed. This argument concerns that defined family, not independently selected forest finite parts.

## 4. Is absence of a diagonal pole sufficient?

Yes, **after** the full identity and meromorphic family have been defined. Absence of a diagonal pole alone cannot prove H8. The previous sections supply the needed finite-domain identity.

On radial sectors the kernels have powers linear in h, rho and integers. The full H2 exponential depends on t only through rho once the source-pair exponent has been replaced; t remains in rational derivative coefficients. The same statement holds for its source-collapse tails and endpoint pullbacks. Differentiating R or masks shifts integer powers or yields logarithms; it does not create an exponent h−t. Fixed old Cauchy/mixed contacts likewise introduce only derivatives and integer shifts in the spectator powers.

Every new Mellin divisor is therefore `a h+b rho+n=0` with nonzero a, and its restriction is `a t+2b/t+n=0`. Multiplication by t gives a polynomial with nonzero quadratic coefficient a. It cannot vanish identically on the diagonal. Higher pole orders or logarithmic radial terms do not change this conclusion. The h and h−1 counterterm denominators do not contain the entire diagonal either. No omitted same-diagonal singular face is exposed by the complete subset enumeration.

Thus eta times the continued remainder has zero restriction at generic nonresonant diagonal points. Resonant points need a separate prescription; an arbitrary irrational t need not avoid every algebraic resonance. The fixed Cauchy faces cannot be omitted and then appealed to as a mysterious absence of poles: they were already part of the product definition. After removing A, the common family and old lower family agree on the convergent positive diagonal h=t>2, which fixes their one-variable continuation.

## 5. All-order anomaly classification: conditional statement survives

From the Gaussian weights, set a=n0−e,c=n2−e,S=a+c,F=a−c and M=n0+n2+2e=S+4e. The physical bosonic charge shift has squared norm −2tac; fermion charge F has minimal chiral weight F²/2. Dimensional balance therefore gives exactly the three levels in the coordinator's notes.

At irrational t, integral levels force e=ac. Then `F²=S²−4e`, and substitution gives

`ellbar=1−S(S+1)/2`, `ell_B=1−S(S−1)/2`, `ell_A=−S(S+3)/2`.

For B the first two inequalities leave S in {−1,0,1}. Since ac=e>=0, a,c have the same weak sign. S=0 forces a=c=0. S=−1 forces one of a,c to be −1 and the other zero, violating n0=a+ac>=0 or n2=c+ac>=0. S=1 leaves (a,c)=(1,0),(0,1), the single Yukawas.

For A, intersecting its inequality with ellbar>=0 leaves S in {−2,−1,0}. The last two cases are as above. For S=−2, the pairs (−2,0),(0,−2) violate insertion-count positivity; (−1,−1) gives e=1,n0=n2=0. This is exactly the single even-root insertion. The finite enumeration in my script checks the derived answer but is not the proof.

This classification is correct only for the specified local Gaussian, mass-analytic, scale-homogeneous anomaly class. The statement that every Ward breaking of a renormalized disorder-line product belongs to that class remains an additional induction/locality assertion. The coordinator's notes state that limitation correctly. No all-order nonlocal conservation theorem results.

## 6. Completion and endpoint claims

For lambda0=lambda2=lambda>0, zeta>0, and theta in an interval where cos(alpha theta)<−1/2,

`V <= −lambda exp(zeta phi)+lambda1 exp(−2zeta phi)`.

Choose a smaller closed interval inside it. The bound is uniform there and makes exp(−volume V) grow faster than exp(constant times phi²). At a finite lattice choose all fluctuating fields in small fixed-width boxes around a common large phi and around such theta values. Boundary-gradient cost is at most quadratic in phi, while the potential contribution is a negative exponential. Integrating these boxes proves divergence on a region of positive measure. This independently verifies the claimed failure of the naive real bosonic cycle. It does not exclude a complex cycle or a fermionic construction.

The primary local Fateev text at (4.1) reads, in part, “For small b one can use the two-dimensional fermion-boson Coleman–Mandelstam correspondence”. Equations (4.2)–(4.3) show the cosh bosonic potential and current-current Thirring term. This supports the limited candidate-reference description in the notes. It is not a statement constructing the original Toda interface sectors. The copied primary PDF is included in the hash record; the local extracted page was inspected. I do not re-audit its reflection amplitudes in this bounded task.

The interface gap can be demonstrated without a physical model: on a direct sum of homogeneous spaces, specifying a scalar operator on each summand leaves every off-diagonal block unconstrained. These blocks may all vanish or be chosen independently. Scalar vacuum/reflection data cannot select an interface intertwiner. Also, vacuum-to-positive-particle matrix elements of c I vanish, so those data cannot determine an identity admixture without separate normalization.

For the weak-limit gap take a two-dimensional space with vacuum e0, `Q=|e1><e1|`, and `O_R=I+exp(iR)Q`. The vacuum expectation is identically 1, but the e1 matrix element oscillates and has no limit. This is an exact counterexample to the proposed inference from a homogeneous VEV, not a Toda counterexample. The coordinator's stated requirement of full weak endpoint limits on a common interface domain is therefore justified.

## 7. Anchor and conditional algebra conversion

From the anchored canonical definition,

`varphi_L(x)=Phi(x)/2−(1/2) integral_L^x dotPhi`,

splitting the integral at c gives the stated factor `Omega_q=exp[−qd.integral_L^c dotPhi/2]` for a vertex with chiral-minus-antichiral momentum qd in the right region. Its spatial derivative prefactor is unchanged away from the split. Equal-time canonical momenta commute on separated intervals, so Omega is multiplicative on concatenated intervals. Three-region grouping gives the same three charge terms. This is the formal string coproduct claimed in the notes.

For a separated left-region vertex of charges (l,r), the local Phi coefficient is (l+r)/2 and `[integral dotPhi,Phi]=−4pi i`. The BCH exponent is therefore

`(−qd/2).((l+r)/2)(−4pi i)=i pi qd.(l+r)`.

It agrees with the noted character. This derivation does not identify the momentum string with a field-difference topological Cartan operator; their representation-dependent identification remains missing.

For the algebra conversion, make the barred weights explicit: `K barA=s barA K`, `K barB=s^−1 barB K`. The assumed mixed relations are

`AbarA−s barAA=a(1−K²)`,
`BbarB−s barBB=b(1−K^−2)`,
`AbarB−s^−1 barBA=BbarA−s^−1 barAB=0`.

Then multiplication by K^−1 or K gives the asserted two diagonal commutators and zero cross commutators. With K0=K^−1 and K1=K, the Cartan actions are q to the powers ±2, q²=s. The cubic Serre coefficient is `[3]_q=q²+1+q^−2=1+s+s^−1`.

For the negative quartics put F_i proportional to barE_i K^(k_i), with k0=−1,k1=+1 and barred weight −k_i. Moving all K factors right in a four-letter word contributes exponent `−sum_(r<s)k_r k_s`. For three letters of one type and one of the other this is −3+3=0, independently of order. Thus the barred quartics become the F quartics with the same coefficients. This independently verifies the conditional level-zero quantum-affine sl2 presentation, assuming nonzero a,b and q−q^−1. It proves neither that these hypotheses hold in Toda nor that the resulting action is faithful. The notes and draft coordinator conclusion correctly preserve those distinctions.

## 8. Actual script output and minimal remaining input

The final fresh run reports 16 PASS lines and no FAIL lines. Among its actual outputs are:

```text
PASS A-dressed node0 Laurent remainder
PASS A-dressed node2 Laurent remainder
PASS all eight new faces reconstructed
INFO face BD02 2*h - 3
INFO face ABD02 2*h - 2/t
PASS no new Mellin divisor contains physical diagonal
PASS two same-chirality regulator derivatives divisible by eta
PASS mixed regulator derivatives divisible by eta
PASS overlap nonuniform coefficient reconstructed
INFO local anomaly candidates {'A': [(0, 0, 1, 1, 0)], 'B': [(0, 1, 0, 1, 0), (1, 0, 0, 1, 0)]}
PASS negative Serre common K factor for -1
PASS negative Serre common K factor for 1
```

An initial script stopped because it assumed a particular SymPy summation result shape. That programming diagnostic is retained in `initial_diagnostic.txt`; it was corrected by directly differentiating the geometric series. No mathematical assertion was changed to pass a failed check. The PASS count is not evidence of a Hilbert-space construction or an all-order proof.

The smallest remaining bridge is not another local residue or scalar VEV. One needs (i) a common associative product encompassing the mixed and Serre families on the tests/domain required for charges, (ii) an induction principle controlling its nonlocal contact breakings beyond the ordinary local Fock class, and (iii) a matched quantum completion with interface sectors and full endpoint weak limits. These are independent requirements. The current finite-order common-family result supplies none of them automatically. The draft coordinator's outcome “unresolved as a quantum symmetry algebra” is therefore the correct final classification.

## Integration

No book integration occurred. Verified independently: the decisive kernels/residues, local h-counterterm dependence, all new collision degrees, the nonuniform overlap, mask product-rule justification, generic diagonal restriction, the conditional anomaly integer proof, real-cycle divergence, anchor string, and conditional algebra conversion. Modified exposition: H8 needs the full-mask product rule spelled out as above. Rejected stronger implications: direct negative-t cutoff uniformity, inference of interfaces from scalar reference data, and inference of an actual quantum algebra from a presentation conversion. The primary reflection normalization itself was outside this bounded re-audit. All writes were confined to the new assigned report and fresh final-audit scratch files.
