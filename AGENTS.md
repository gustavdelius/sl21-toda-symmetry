# Research working rules

This is a standalone research repository extracted from the affine Toda book. Gustav Delius is the editor. Codex and Gemini act as mathematical referees, not authors of book text.

- Read README.md, STATUS.md and RESUME.md before continuing research.
- Preserve imported reports, scripts, outputs, source copies and manifests byte for byte. `referee_scratch/` is a permanent evidence archive despite its name. Corrections belong in fresh reports/addenda and fresh scratch directories.
- Write new referee reports under `briefs/` and work under a new `referee_scratch/<tag>/`. Assign fresh worker paths explicitly. Current navigation/handoff documents may be updated when the editor requests consolidation.
- Do not edit or integrate into the original affine_toda checkout. Do not render the book. Historical reports retain old absolute paths, book commit hashes and scope restrictions.
- Independently derive or compute each claimed result. Do not treat a previous report or a PASS count as proof. State a verdict per item, evidence, exact scope, and unresolved steps. Retain failures and corrected checking-code diagnostics.
- Distinguish formal meromorphic identities, regular physical restrictions, raw cutoff limits and physical operator realizations. Do not upgrade one into another. Keep complete source/test jets, full Wick matchings, constituent phases and regulator ancestry.
- Distinguish K_A (the current tail) from F_A^end (the mixed endpoint). Insert the interaction map once. Count source-pair and current-improvement counterterms once. A single original derivative slot is consumed once.
- Spawn subagents only when the user authorizes them. At most two workers concurrently; every worker prompt says “Do not spawn subagents.” Give each its own new scratch directory and checkpoint.
- Commit, push, publish or create a remote only when the editor asks. Never use destructive Git cleanup to remove another session's files.
- Read exact primary sources before making literature or priority claims. Third-party literature copies retain their own rights; the code licence does not relicense them.
- Use `tools/run_recent_checks.py` or an equivalent temporary copy for historical scripts that write outputs. Never overwrite archived evidence while rerunning checks.

The original book rules are preserved as provenance/book-AGENTS.md for context. Their book-layout and original-task write paths are not instructions to modify that other repository.
