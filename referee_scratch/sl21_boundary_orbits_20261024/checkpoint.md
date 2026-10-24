# Boundary-orbit audit checkpoint

2026-10-24: Read AGENTS.md and prior reports; fresh scratch only, no old scripts imported or rerun. Inputs hashed. Commit 745586e6a6f7bca5ed91a1e9cf89d6ef5df387f0.
2026-10-24: Fresh local/Ward, subset-degree, even-root derivative and shuffle audit complete: audit_m2.py 123 PASS. Initial output retained after correcting omitted minus action sign in the new test.
2026-10-24: Inherited marginal face independently reconstructed at arbitrary insertion positions, with every K tail: inherited_faces.py 8 PASS.
2026-10-24: Mixed finite contact recomputed from Cauchy strip edges and direct quadrature: mixed_contact_recheck.py 14 PASS. Canonical matching remains a finite prescription choice.
2026-10-24: Endpoint free-direction neutrality, action-period characters, anchored canonical weights and conditional boundary formulas: boundary_checks.py 67 PASS. Frozen-boson Gaussian probe: 14 PASS. Total 226 PASS, 0 FAIL in final runs.
2026-10-24: New report briefs/referee_sl21_boundary_orbits_report.md completed. Existing inputs hash-verified; tracked/index diffs empty. No subagents, existing-file changes, git mutations, commits/pushes, or render.

## Carry forward

- Independently audit this report before extending it. Previous scripts/reports remain claims, not evidence.
- Retain formal prescription-dependent contour scope; the mass-two audit supports the existing analytic family, not an associative/global symmetry algebra.
- Period character: U and F transform by s=exp(-2*pi*i/t); FA by its inverse. The action period includes a one-chirality fermion sign flip. It is not yet an automorphism of a constructed physical quantum domain.
- On an assumed scalar-clustering vacuum orbit: AA endpoint A_alpha*s^(-n), BB endpoint B_alpha*s^n. Invariant vacua force both amplitudes to zero at generic s. Neither vacuum orbit nor nonzero amplitudes established.
- Corrected left-anchored split shifts LEFT boundary, preserves RIGHT: [nL,Q_q]=-2q Q_q, [nR,Q_q]=0, [N,Q_q]=2q Q_q; antichiral signs reverse. Right-fixed r gives DeltaBB=B_alpha*s^r*(1-s^(-N)), DeltaAA=A_alpha*s^(-r)*(1-s^N), conditionally.
- Frozen Gaussian fermion condensate is a diagnostic only; it omits boson fluctuations and the true endpoint counterterms. Do not use it as the Toda endpoint VEV.
- Next: construct/specify a quantum integration cycle and boundary sector compatible with the action period, then calculate full renormalized weak endpoint limits. Also test whether canonical mixed finite matching and analytic same-chiral products admit one associative extension. Higher-order Serre/conservation remain open.
- Finite-order CPT remainder is not uniform in the infinite-volume boundary limit. No physical q sign, affine sl2/sl(2|1) algebra, soliton representation or nonzero closure has been established.
