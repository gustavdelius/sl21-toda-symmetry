# Coordinated investigation of the sl(2|1) Toda nonlocal-charge algebra

Date: 2026-10-24. Prepared by Codex (GPT-6).

Starting checkout: `/home/gustav/Git/affine_toda`, branch `main`, commit `745586e6a6f7bca5ed91a1e9cf89d6ef5df387f0`.

This is an execution brief and research plan, not a report establishing a symmetry. Many relevant reports are untracked: the commit alone does not identify the reviewed evidence. Record file hashes before using them.

## Mandate and outcome

You are the coordinating adversarial mathematical-physics referee. The editor has asked for a separate coordinating task that uses subagents to drive this investigation to a conclusion. Start work immediately after reading this plan and `AGENTS.md`; do not return only a restatement of the plan. Continue from one completed milestone to the next without asking for routine permission to continue.

Investigate the nonlocal charges of the root-system affine sl(2|1) Toda action. Determine whether they realize a symmetry algebra in a specified quantum theory. Do not assume that algebra is affine sl(2|1), affine sl2, or an algebra inferred from a local change of variables. Derive the generators, grading, relations and coproduct, and distinguish a representation of a presented algebra from identification of the image or proof of faithfulness.

An honest conclusion can be:

1. An established algebra action in a precisely specified quantum completion and sector/domain, with every necessary relation and boundary limit justified and independently audited.
2. A verified obstruction to the proposed construction, with its exact scope, and a classification of whether allowed changes of prescription or generators can remove it. Failure of one regulator does not establish nonexistence of a symmetry.
3. A sharply delimited unresolved problem after substantial independent attempts, with the strongest conditional result, the minimal missing mathematical input, and an explicit account of the approaches tried. This is an unresolved conclusion, not success at establishing the symmetry.

Do not force a positive answer or call a finite-order contour identity an established symmetry algebra. Reaching a milestone alone is not a reason to stop if useful authorized work remains. Equally, avoid repeatedly rerunning checks or spawning agents to rephrase the same unresolved issue. When a real impasse remains after distinct serious approaches, report it candidly.

## Rules, scope and explicit exceptions

Read `AGENTS.md` before acting. The referee restrictions continue to apply:

- Other agents may be editing files in this checkout at the same time.
- Do not edit any file that predates this coordinated investigation. In particular, the book, existing reports, this plan, existing checkpoints, indexes and bibliography are read-only.
- Do not run git commands that change the index or working tree: no commit, checkout, stash, reset, rebase, restore, merge or push. Do not create a worktree or change branches.
- Do not run `quarto render` or any book build that writes into the checkout.
- Keep scratch work, generated source copies, scripts, output and checkpoints under new subdirectories of `referee_scratch/sl21_symmetry_coordinator_20261024/`.
- Write the coordinator's new report to `briefs/referee_sl21_symmetry_coordinator_report.md`. If it already exists when you start, choose a fresh numbered suffix and record it; never overwrite it.
- Worker reports may be created only at fresh `briefs/referee_sl21_coord_<stage>_<worker>_report.md` paths assigned explicitly by you. If a path exists, use a fresh suffix. Each worker gets its own new scratch subdirectory.
- You may revise your own newly created working outputs during this run. Treat a completed earlier pass as immutable; later corrections belong in new reports/addenda. Maintain the checkpoint in your newly created scratch directory.
- No book integration, label upgrade, publication, commit or push is authorized by this plan.

The user's request explicitly authorizes **the coordinator** to spawn subagents. This overrides the repository/brief-template prohibition on delegation for this coordinator only. Every worker prompt must say **“Do not spawn subagents.”** Keep at most two workers active at once, so the coordinating team has at most three active agents. During startup use one worker while you audit the inherited evidence; expand to two only when it helps distinct work. Reuse idle workers or replace them deliberately; do not build a delegation tree.

This plan also explicitly uses the requested existing checkout and new-file-only referee workflow instead of the generic new-worktree workflow. The user's no-commit restriction overrides the brief skill's request to commit the brief. Scratch/checkpoint files remain in the permitted scratch tree rather than being copied into additional brief directories. These are scope resolutions, not grounds to seek another approval.

Do not create additional user-owned tasks, scheduled automations or messages to other conversations unless the editor separately asks. Coordinate your own workers with the collaboration tools. Workers inherit the configured model; do not change models without authorization. Communicate meaningful findings and decisions regularly in the coordinating task, and checkpoint completed work so context compaction is safe.

## Read first, in this order

All paths are relative to `/home/gustav/Git/affine_toda` unless absolute.

1. `AGENTS.md`, `README.md`, and the conventions section in `01-introduction.qmd`.
2. `briefs/referee_sl21_triple_extension_report.md` and `referee_scratch/sl21_triple_extension_20261024/checkpoint.md`.
3. `briefs/referee_sl21_triple_kernel_report.md` and `briefs/referee_sl21_reflection_bridge_audit_report.md`.
4. `briefs/referee_sl21_vacuum_limits_report.md`, `briefs/referee_sl21_product_compatibility_report.md`, `briefs/referee_sl21_boundary_global_audit_report.md`.
5. For the local Ward calculation: `briefs/referee_sl21_mixed_regulator_report.md`, `briefs/referee_sl21_mixed_audit_report.md`, `briefs/referee_sl21_ward_m2_report.md`, `briefs/referee_sl21_current_renorm_report.md`, and `referee_scratch/sl21_mixed_regulator_20261023/checkpoint.md`.
6. Before any four-current claim: `briefs/referee_sl21_same_chirality_report.md`, `briefs/referee_sl21_serre_boundary_report.md`, `briefs/referee_sl21_serre_forest_report.md`, and `briefs/referee_sl21_serre_m2_report.md`.
7. For the original action and candidate currents: `briefs/referee_sl21_kink_modes_report.md`, `briefs/referee_sl21_bernard_leclair_report.md`, `briefs/referee_sl21_dual_screening_report.md`, `briefs/referee_sl21_nonlocal_report.md`. Follow their exact chapter/source references and reconstruct definitions directly.

Primary Bernard–LeClair source:

`/home/gustav/Zotero/storage/QJH8RCYM/Bernard and LeClair - 1991 - Quantum group symmetries and non-local currents in.pdf`.

Reference-theory primary sources, already copied under `referee_scratch/sl21_vacuum_limits_20261024/`:

- `fateev_2001.pdf`, [hep-th/0103014v2](https://arxiv.org/abs/hep-th/0103014): reflection formulas (24)–(25), (130)–(132). Rendered pages 6 and 33 and extracted text are available. Read the actual source when checking a sign or a normalization; extracted text can lose overbars.
- `fateev_2017.pdf`, [SIGMA 13 (2017) 080](https://sigma-journal.com/2017/080/sigma17-080.pdf): first hierarchy, (4.1)–(4.2).
- `lashkevich_pugai_2017.pdf`, [arXiv:1610.04926](https://arxiv.org/abs/1610.04926): operator definitions, normalized form factors, recursion and reflection conventions.
- A 2019 Lashkevich–Pugai file is also present, but do not mix its operator labels/normalization with the 2017 convention without an explicit conversion.

Earlier output counts are claims, not evidence: the last root pass has 78 checks, its kernel worker 23 and its reflection worker 87. The root's initial failed branch-orientation check is deliberately retained. An independent derivation of the positive-imaginary-separation boundary values corrected that checking-code error. Do not infer correctness from the final count or silently delete failures.

## Conventions to reconstruct and keep explicit

Use chiral covariance `X^i(z)X^j(w) ~ -delta_ij log(z-w)`, ordinary Euclidean area and the reference Wick convention with no additive covariance constant. Scalar products are bilinear.

\[
a_0=(ib/\sqrt2,b/\sqrt2,i),\quad a_2=(ib/\sqrt2,-b/\sqrt2,-i),\quad
t=b^2,\quad T=a_0+a_2,\quad d=-T/t,\quad k=d+T,\quad \rho=2/t.
\]

Thus `a0^2=a2^2=-1`, `a0.a2=1-t`, `ai.d=1`, `d^2=-rho`. The physical two-boson Yukawa momenta are the first two components, and are null. Global bosonization, compactification and spin structures are not fixed by these local products.

\[
J_A=:e^{-dX}:,\quad J_B=:C e^{dX}:,\quad C=a_0\partial X,\qquad
H_B^{(1)}=-t\lambda_2:e^{(d+a_2)X+a_2\bar X}:.
\]

Both candidate charges are even. Charge normalization is `(2pi i)^-1`, with density `J+H`. The first A tail starts at the even-root coupling. The action convention is

\[
S_{\rm int}=\pi^{-1}\int d^2x\,(\lambda_0Y_0+\lambda_2Y_2+\lambda_1V_1+\kappa P),
\quad
\kappa=-\lambda_0\lambda_2\epsilon^{2t}/t,\quad\kappa_R=0,
\]

where `Yi=exp[ai.(X+barX)]`, `V1=exp[-T.(X+barX)]`, `P=exp[T.(X+barX)]`. The Yukawa couplings have mass order one; `lambda1` and `kappa` have mass order two. This bookkeeping differs from engineering dimensions. The superscript/subscript 2 in `lambda2` is a node label, not perturbative order.

\[
\nu=e^{-i\pi/t},\quad s=\nu^2,\qquad
\varepsilon(n,m)=(-1)^{(n_L+n_R)(m_0+m_2)}
e^{-i\pi(n_Lm_R-n_Rm_L)/t}.
\]

Do not identify `s` with a physical quantum-group parameter, or drop a sign by calling it a gauge choice, without deriving the identification and any Klein-factor twist. Use generic irrational numerical probes.

Candidate endpoints, through their stated mass orders, are

\[
F_A=e^{-k(X+\bar X)},\quad U=e^{(d+a_2)(X+\bar X)},\quad F=e^{k(X+\bar X)},
\]
\[
B_{AA}=\left(\frac{t}{1-t}\right)^2\lambda_1F_A,
\]
\[
B_{BB}=\lim\left[t^2\lambda_2U-
\frac{t^2\lambda_0\lambda_2}{t-1}\epsilon^{2t-2}F
-t\lambda_0\lambda_2\epsilon^{2t}C\bar C F\right]_{S_\epsilon}.
\]

For `B barB`, the raw circular contact coefficient is `t`, while the specified finite binary addition is `t(t-1)`, giving `t^2`. Both values pass their corresponding lowest projected triple test; that test does not select a unique physical normalization.

## Present status: precise ceiling on inherited claims

The current claims are formal, prescription-dependent contour identities. The latest actual triple calculation establishes only the coefficient linear in `lambda2`, for the stated test-function family:

\[
\left.[B_gA_f-s^{-1}A_fB_g,\bar B_I]\right|_{\lambda_2}
=(s^{-1}-s)A_f\frac{\nu}{2\pi i}t^2\lambda_2U(L).
\]

Here `I=[L,R]`, `A_f` is supported strictly inside `I`, `g=1` near `I`, moving-current poles are extended first as Cauchy distributions, the inner source hole is a circle, and ordered contour heights are fixed multiples of its radius. Details are X5–X17 of the latest report. For `Re rho<0`, the simultaneous interior degree is `1-rho` and the new AB face degree is `-rho`; the proper mixed face must be handled before invoking those degrees. Meromorphic continuation is part of the chosen prescription. Do not extend this result to arbitrary spatial tests or arbitrary iterated products by assertion.

The left-anchored boundary weights are `U(L) A=s^2 A U(L)` and `U(R) A=A U(R)`. For a fixed right vacuum label `r` and `N=r-l`, scalar endpoint limits would give

\[
\Delta B_{AA}=\mathcal A s^{-r}(1-s^N),\qquad
\Delta B_{BB}=\mathcal B s^r(1-s^{-N}).
\]

Those limits, the interface sectors and the full common nonlocal-charge domain are not constructed. Earlier algebraic overlap checks and matrix witnesses establish consistency of proposed relations, not a Toda representation.

The local first sine-Liouville hierarchy has conditional analytic-reflection values

\[
G_A=[\lambda_1\gamma_f(-t)]^{1/t-1}
\frac{\Gamma(2-t)\Gamma(2-1/t)}{\Gamma(t)\Gamma(1/t)},
\]
\[
G_U=\frac{\lambda_0}{1-t}
\left(\frac{\sqrt{\lambda_0\lambda_2}}{-t}\right)^{2/t-2}
\frac{\Gamma(2-1/t)}{\Gamma(1/t)},\qquad
\gamma_f(x)=\Gamma(x)/\Gamma(1-x).
\]

Use generic irrational `-1<t<0` and explicit coupling logarithms initially. The second value comes from a simple reflection zero times the reflected field's simple identity pole. The normalized neutral form-factor degeneracy is not itself a normalization proof. Scalar-identity statements require a regular analytic family, crossing, charge selection and completeness on a specified reference sector. Their extension to the original Toda Ward endpoints requires additional all-order matching.

The latest audit corrected the wording before V21 of the vacuum report: its logarithm changes are a **compensating coupling trajectory** accompanying a translated cycle. Direct pullback at fixed couplings has the opposite signs. Use X28–X29 of the latest report, not the old wording.

Gaussian normalization is substantive. For a finite change `H=diag(c,c,0)` in the physical two bosons with fermion normalization fixed, the transported Ward coefficients acquire `exp(2c)` and `exp(c)`; the UY0 and Y0Y2 fusion coefficients acquire `exp[c(1-t)]` and `exp(-ct)`. At fixed action both candidate endpoint amplitudes acquire `exp(c/t)`. On the stated real branch their ratio in the reference convention is

\[
\mathcal B/\mathcal A=-[-\gamma_f(t)\lambda_0\lambda_2/\lambda_1]^{1/t}.
\]

This is a conditional normalization benchmark, not an identified particle mass. In new coupling coordinates the formula changes; unrelated independent charge rescalings also change the ratio. Do not take a finite perturbative remainder to infinite volume without uniform control: these vacuum powers are nonanalytic in the couplings.

## Work programme and acceptance criteria

### Stage 0 — independent intake and claim ledger

Create a fresh scratch root, snapshot the commit, branch, status and hashes of all inputs, and read the actual action/source conventions. Make a claim ledger with status, dependencies, prescription, independent evidence and next falsifiable test for each claim.

Independently audit the latest triple finite contact and collision/continuation argument before building on them. Reconstruct the momentum contractions, derivative residues, cocycles, strip-edge terms and endpoint weights. Verify whether its collision argument controls every proper face and nonuniform nested remainder for its stated test family. Recheck the reflection zero/pole and Gaussian transport at the point needed for subsequent work. Old scripts are comparisons only.

First worker assignment: a bounded independent audit of the triple extension and a fresh classification of the next double-Yukawa collision faces. Tell the worker not to infer any all-order result. While it works, perform your own derivation of the critical first-order contact and audit the source/action convention. Compare evidence before accepting the starting result. If an error appears, resolve it in a new report before extending it.

Acceptance: each premise used in stage 1 has fresh evidence and a stated scope; any rejected premise is corrected or the dependent route is suspended. Do not repeat every historical computation when it is not a premise of the next step.

### Stage 1 — the mass-squared triple identity

Evaluate the coefficient of `lambda0 lambda2` in the same four-word projection. Include the regulated double-Yukawa source integral, the composite-potential insertion, relevant current tails and improvements, both endpoint subtractions, and finite matching on every nested mixed face. Reconstruct the combinatorial factors from the action expansion; do not infer the result from the desired algebra or a binary Ward identity alone.

Use at most two bounded worker tasks: one derives the full kernels and collision geometry; another independently derives the subtraction/finite-contact terms. Keep their derivations independent until both have concrete results. The coordinator must recompute a decisive coefficient or identity, not merely count agreeing verdicts.

Specify the order of integrations, contour heights, exclusions, Cauchy extensions, test-function support, analytic regulators and continuation path. Exhibit a common initial domain, or a complete alternative distributional extension with subtraction consistency. Classify overlapping forests, not just the superficial degree of the total cluster. Compute the actual defect and its change under finite binary matching. Establish whether a primitive contact is needed and what freedom remains.

Acceptance: either a derived projected identity with all local coefficients and remainder estimates, or an explicit residual term with independent verification and its allowed-remedy analysis. Document arbitrary-test contacts that remain undetermined. Do not call a single-prescription failure a no-go theorem.

### Stage 2 — one common product and all required finite-order relations

Evaluate the companion projection involving `barA` at its first nonzero orders. Then test the same-chirality/Serre relations and their mixed overlaps using the same regulated product, rather than treating separate analytic continuations as iterates of an associative product.

Determine how much of the density product is needed to define the intended integrated charge products. If arbitrary-test extensions are necessary, calculate and classify the previously invisible derivative-of-test-function contacts. Verify parenthesization, ordered-contour deformations, cocycles and all relevant overlapping diagonals. Distinguish universal exchange phases from integrated relations. Derive the coproduct from separated spatial regions and the nonlocal anchor; do not import the book's or a source's coproduct without the convention conversion.

Acceptance: a single specified product supports all claimed finite-order relations and overlap identities, or a concrete obstruction/freedom is identified. State the exact perturbative ceiling. An associative matrix witness remains a useful consistency check, not field-theory evidence.

### Stage 3 — an all-order mechanism

Investigate mechanisms that could establish conservation and the relations beyond finitely many checks. Possible routes include classification of allowed local anomalies and resonances, a screening/cohomological construction, a regulator with an exact identity, or a rigorously matched quantum-integrable reference construction. These are avenues to test, not premises to assert.

Search and read primary literature, including relevant work by the editor, before claiming novelty or importing a purported theorem. Separate exact results, physical assumptions and conjectured bootstrap identifications. Every imported statement needs a precise match of action, coupling normalization, operator normalization, grading, cocycles, compactification and sector. Quote the exact source line supporting a literature claim. A source about a different Toda or supersymmetric action is not a substitute for the root-system action here.

After the first common-product test is secure, this work can run in parallel with stage 4. Do not delay every conceptual investigation until an unbounded sequence of higher-order expansions has been performed.

Acceptance: a complete all-order argument for the stated scope, or a precise residual anomaly class/induction step/assumption that the available evidence does not establish. Do not use finite-order success or expected integrability as proof.

### Stage 4 — quantum completion, endpoints and interface sectors

Specify a candidate quantum completion explicitly: finite regulator, integration cycle or operator realization, coupling/field trajectory, normal ordering, composite and finite counterterms, global fermionic sectors, and infinite-volume limit. The analytic-reflection completion is the existing candidate to investigate; do not silently equate it with the original real-field path integral. If inequivalent choices matter, analyze them conditionally before asking the editor to choose a physical theory.

Derive all-order matching of the actual Ward endpoints to the reference operators, including allowed finite mixing and nonperturbative normalization ambiguities. Use the Gaussian transport formulas and reflection zero/pole as checks, not as an assumed equality of theories. Examine regularity and completeness assumptions behind the identity-like form factors.

Construct or identify interface sectors with independently specified left/right vacua and a common domain for the nonlocal charges. Establish the endpoint limits needed in matrix elements and explain the topology of convergence, clustering assumptions, zero modes and charge-domain stability. Homogeneous translated cycles or formal sums of sector projectors do not by themselves supply interfaces. Control the order of ultraviolet, infinite-volume and small-mass limits.

Acceptance: a stated quantum theory/domain and justified endpoint closure, with each remaining physical or mathematical assumption explicit. If closure is only conditional, preserve that label throughout the algebra claim.

### Stage 5 — identify the algebra and independently attack the conclusion

Assemble the generators, Cartan/topological operators, grading, deformation parameter, central elements, coproduct and complete required relations. Determine which algebra has been shown to act and whether additional relations imply a quotient. Do not assert equality of an abstract presentation and its image without the required independence/faithfulness argument. Keep special-coupling degeneracies separate from generic irrational coupling.

Assign a final independent adversarial audit with access to the definitions and primary evidence, not only the synthesis. Ask specifically for a counterexample or the weakest unsupported implication connecting local currents to the global algebra. Resolve substantive disagreements by fresh calculation. Do not upgrade a theorem label on agreement alone.

Acceptance: a self-contained conclusion classified under the three outcomes above, with dependency table and proof/obstruction/gap at each link. No book edit is part of this stage.

## Worker contract and execution discipline

Every worker prompt must include its exact scratch/report paths, the specific claim to check, relevant definitions, independent-evidence requirement, **“Do not spawn subagents,”** and all no-edit/no-git-mutation/no-render restrictions. Assign one bounded mathematical deliverable at a time. Require a verdict of correct / wrong / unclear, derivation, fresh computation where appropriate, actual output, residual uncertainty, and replacement text with file/line for a correction. A verdict without evidence counts as unclear.

Do not allow workers to edit one another's files. New checks print PASS/FAIL per substantive check and retain failures and convergence diagnostics. A larger check count is not an advancement in proof status. Numerical probes should avoid rational-coupling/root-of-unity degeneracies and should test scaling, singularity orders and stability where those matter.

Keep an append-only progress log in your new scratch checkpoint, one line per completed item, plus a compact current-state summary. Record active workers, exact file ownership, accepted/rejected premises and the next action. Before ending a pass, record all inputs and outputs and verify preservation. A successor must be able to resume without rerunning the investigation from scratch.

Use an explicit issue ledger for obstructions. For each one, record the failed implication, a minimal counterexample or calculation, available scheme freedom, distinct attempted remedies, and the smallest additional assumption that would suffice. Escalate only substantive choices the user has not already authorized; do not ask for routine permission to read sources, compute, create new scratch/report files, or start the next bounded worker task.

## Final report and communication

The coordinator report must give the date, actual model and reviewed commit(s), plus input hashes for untracked evidence. Include a concise outcome, per-item verdicts, the derivations needed to assess them, links to fresh scripts/output, a complete dependency/assumption ledger, and an Integration section stating that no book integration occurred. Distinguish finite-order checks, formal continued identities, conditional reference-theory statements and actual original-theory results.

Give regular concise updates with a finding or decision, not unchanged polling. At each major milestone state what is now established, what remains conditional, and the next discriminating calculation. Keep doing authorized work after a milestone unless there is a genuine stopping condition. At final delivery, say clearly whether an algebra was established, a construction was obstructed, or the investigation remains unresolved with a precise missing input. Report pre-existing-file preservation and no-commit/no-push compliance.
