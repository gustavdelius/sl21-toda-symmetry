# First-Yukawa Serre boundary calculation, 2026-10-24

- Read AGENTS.md and the same-chirality report/checkpoint; reviewed commit 745586e6a6f7bca5ed91a1e9cf89d6ef5df387f0.
- No subagents; all existing files are read-only.
- Independently derived the finite circular-hole correction to the scalar source primitive and retained the contour-height contribution.
- Independently reconstructed the four-current/Yukawa and all three tail scalar derivatives from Gaussian Wick contractions. Node-0 lowest scalar vanishes.
- Evaluated the scalar endpoint sum: p(s)=p(1/s)=0 for p(x)=(1-x)(1-s*x)(1-x/s). No quartic scalar finite contact was tuned.
- Checked exact circular Green boundaries and direct three-dimensional endpoint quadrature with finite heights and holes; the latter is a cutoff trend check, not a complete operator integral.
- Defined an additional scalar-only independent-exponent extension R_(g,rho)=Delta^g*product(z_i-a)^(-rho). Initial domain Re(g)>1, Re(rho)<0 permits equal spatial intervals. Continue individual scalar integrals, then restrict g=rho at regular generic points. This does not assert a full quartic identity off the diagonal or identify the bare physical-cutoff finite part.
- Computed all divergent B-only product Taylor levels in -1/3<Re(rho)<0: pair levels 0 and 1 (bulk first jet zero but endpoint term nonzero); triple levels 0,1,2, including descendant and endpoint terms. Triple prescription uses a span cutoff and meromorphic angular finite parts.
- Found the proper ABB scalar primitive is -lambda_2/rho * partial_a[R_2*(S_2+beta)]. The beta term comes from the nonzero infinity moment and cannot be dropped.
- Independent negative-rho angular finite-part quadrature and nonconstant triple radial tests validate the subtraction coefficients. Four scripts completed: boundary_checks 52 PASS; kernel_audit 24 PASS; boundary_quadrature 7 PASS; subtraction_audit 11 PASS. Total 94 PASS, 0 FAIL.
- New report: briefs/referee_sl21_serre_boundary_report.md. Scope: leading scalar boundary test and lower product subtraction data, NOT a complete first-order massive operator Serre identity or a symmetry algebra.
- Next: assemble embedded BB/BBB/AB/ABB proper-face subtractions in all four operator words with nonconstant spectators, verify one common extension and compare its finite normalization with the finite-height/circular-hole prescription. Only after that, undertake mass-squared Serre diagrams. Physical mixed-boundary closure remains separately open.
- All starting input hashes match; tracked/index diffs remain empty. No existing files edited, no subagents, commits, or pushes.
