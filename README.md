# Fundamental sl(2|1) scattering and affine supersymmetric Toda theory

The goal is to construct an S-matrix based on the R-matrix for the three-dimensional fundamental representation of $U_q(sl(2|1))$ in Delius–Gould–Links–Zhang, [*On Type-I Quantum Affine Superalgebras*, hep-th/9408006](https://arxiv.org/abs/hep-th/9408006), especially equations (4.10)–(4.11), and determine whether it describes solitons in an affine supersymmetric Toda field theory.

In simple terms, an S-matrix predicts what happens when particles collide. Solitons are localized field configurations that can behave like particles. The paper supplies an algebraic building block for their possible scattering. This project asks whether that building block can be completed into a consistent scattering theory and matched to solitons of a specific field theory.

There are two connected parts to the goal:

1. **Construct the scattering theory.** Determine the amplitudes for the fundamental triplet and its antiparticles, the scalar factors and energy dependence, and any bound states. Check the required scattering consistency conditions.
2. **Identify the field theory.** Determine which affine supersymmetric Toda action, if any, has solitons with this symmetry and representation. Derive the relation between its coupling and the scattering parameters, and compare its spectrum and semiclassical scattering with the candidate S-matrix.

The recent nonlocal-charge work addresses one part of the second question: whether proposed conserved quantities form a consistent symmetry that could constrain soliton scattering. Their products have singular contacts, so the same prescription must work for different groupings and algebra relations. **Even completing that symmetry construction would not by itself establish the original S-matrix/Toda identification.**

## Where the project stands

The archived reports record the following progress as of 24 October 2026. These are scoped research results, not a completed physical identification:

| Stage | What has been achieved |
|---|---|
| Algebraic scattering candidate | An earlier report constructs the fundamental/dual channels, crossing-compatible scalar normalizations and a rapidity gradation, with graded Yang–Baxter checks. Scalar factors and initial fusion have been explored; a physical scalar choice and full attractive bootstrap remain unresolved. |
| Toda and soliton tests | Candidate actions, classical kinks, fermionic modes and boundary effects have been investigated. The fundamental triplet has not been identified with the physical soliton multiplet. |
| Basic multiplication | A common Gaussian product has been constructed in the tested sectors. |
| Algebra relations | Both proposed Serre relations—particular combinations of four charges—pass through second order in the tested chiral sector. |
| One opposite-chirality charge | A single interaction prescription handles the tested five-current configurations, including symmetric contacts and second-order compatibility with the relations. |
| Two opposite-chirality charges | The first genuinely new case, two currents of each chirality, gives matching groupings through second order at generic coupling under the common analytic-continuation prescription. |

Two useful corrections emerged from the checks. Contact terms must include derivatives of the source and spectator coefficients, not only derivatives of fields. Also, the simple averaging rule that worked with one barred current fails with two; independently balancing disjoint contacts fixes the coefficient without introducing a new adjustable parameter.

“Barred” denotes the opposite chirality, roughly the other direction of propagation. “Through second order” means a finite perturbative result, not a proof for arbitrary interaction strength. Analytic continuation is part of the declared prescription; these results do not establish convergence of the raw collision limits.

The [scattering construction report](briefs/referee_sl21_soliton_report.md) records both the candidate and its limitations. In particular, it reports a failure of positive-metric unitarity in the tested gradation; algebraic braiding unitarity alone does not settle that issue. The [kink-mode report](briefs/referee_sl21_kink_modes_report.md) records the dynamical tests and the unresolved soliton representation.

## What remains open

The central question remains whether this candidate can describe soliton scattering in a specified affine supersymmetric Toda theory. We still need the physical soliton multiplets and charge action, an action-derived coupling/rapidity dictionary, the scalar factor and spectrum, and an appropriate physical interpretation of unitarity and reality. An auxiliary R-matrix in a Toda model is not enough to establish its asymptotic soliton S-matrix.

Within the current nonlocal-symmetry investigation, the latest test found a possible exceptional coupling,

`t = (7 + sqrt(17))/4 ≈ 2.7807764064`.

A bulk contribution is nonzero there. Whether it survives after **all** current-tail and contact contributions are added is unknown. It is not yet an anomaly or a failure of associativity.

The other main gaps are relations containing two barred generators, arbitrary words and higher perturbative orders, and an actual physical operator realization. A passing script is evidence for the formula it checks, not proof of that entire programme.

## Start here

1. [Original scattering construction and Toda identification tests](briefs/referee_sl21_soliton_report.md).
2. [Current technical status](STATUS.md): precise scope and conventions.
3. [Resume the investigation](RESUME.md): a self-contained handoff for the next session.
4. [Latest two-barred report](briefs/referee_sl21_two_barred_report.md) and its [independent audit](briefs/referee_sl21_two_barred_audit_report.md).
5. [Research index](INDEX.md): earlier reports, background material and recorded evidence.

## What is in this repository?

- `briefs/`: the reports and original research briefs, preserved unchanged.
- `referee_scratch/sl21*/`: derivations, checking scripts, outputs, failed attempts and checkpoints. Despite its historical name, this is a **permanent evidence archive**.
- `code/super_toda/`, the selected `.qmd` files and `references.bib`: background copied from the affine Toda book. This is not a complete book checkout or a buildable book project.
- `literature/`: an additional locally available primary source; other source copies already live in the historical scratch directories.
- `provenance/`: source paths, hashes, licensing context and migration checks.
- `tools/`: a read-only archive verifier and a small check runner that works in a temporary copy.

The material was consolidated from Gustav Delius's `affine_toda` repository. Original files were retained there as a historical backup; this repository is the new place to continue the research. It has independent Git metadata and is not a worktree of the book. No remote or public publication has been configured.

## Checking the archive

```sh
python3 tools/verify_archive.py
```

The Python research checks use the dependencies in `requirements.txt`. A small recent set can be repeated without modifying the archived evidence:

```sh
python3 -m venv .venv
.venv/bin/python -m pip install -r requirements.txt
.venv/bin/python tools/run_recent_checks.py
```

The runner copies the relevant research directory to temporary storage before executing the scripts. It reports new results without overwriting the old outputs. It is a smoke test, not a rerun of every historical calculation. Some older scripts deliberately fail, take substantial time, or refer to paths and Git snapshots in the original book checkout. See [provenance and portability](provenance/README.md).

Please preserve the archived record and follow [AGENTS.md](AGENTS.md). See [NOTICE.md](NOTICE.md) for attribution and the distinction between the code licence, book text and third-party literature.
