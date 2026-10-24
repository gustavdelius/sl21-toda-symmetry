# Independent mixed-product audit, 2026-10-23

- Read AGENTS.md, latest report and checkpoint; recorded reviewed inputs and commit 745586e6a6f7bca5ed91a1e9cf89d6ef5df387f0.
- No earlier check script imported or rerun as evidence.
- Re-derived local bosonization, candidate current poles, three-pole composite response, and the two-endpoint Ward polynomial.
- Independently confirmed the finite first-order contact using a smooth 2D Cauchy kernel and its exact contour moments.
- Found that the reported raw interior collision count uses one extra integration dimension; explicit relative powers and the leading reflection cancellation are now derived.
- Completed a conservative sector argument starting at Re(t)>2, avoiding reliance on the unproved sharper Re(t)>1 bounds.
- Fresh audit script currently passes 65 checks. Earlier scripts were not imported or rerun. Writing the audit report next.
- Wrote briefs/referee_sl21_mixed_audit_report.md: finite contact verified; raw interior collision powers corrected; formal identity supported by a conservative continuation from Re(t)>2.
- Added an explicit nested-face Taylor remainder bound, epsilon times integral R^(2 Re(t)-5), to justify treating the extra Yukawa insertion beyond separated points.
- The sharper bare-cutoff half-plane Re(t)>1 remains unproved by this audit; it is unnecessary for the specified generic meromorphic continuation.
- Final independent script output: 65 PASS, 0 FAIL. Checked report links, display delimiters and equation tags.
- All input hashes match; tracked/index diffs remain empty. Existing files unchanged. No commit, push, subagents or book render.
- Remaining: physical interacting realization of the normalized product; spatial boundary values of the renormalized endpoints; same-chirality/Serre relations; exceptional couplings. No affine algebra identified.
