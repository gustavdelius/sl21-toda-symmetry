# Import provenance and portability

The source repository is `/home/gustav/Git/affine_toda`. The exact source commit, branch, tracked-file status and hash of every imported file are in `import-manifest.json`.

All `briefs/referee_sl21*` files and all material under `referee_scratch/sl21*` were copied, together with the superalgebra background briefs, their scripts, the relevant book chapter, introduction, bibliography and `code/super_toda/`. Python bytecode/cache directories were excluded. The Bernard–LeClair PDF referenced by the reports was also copied from the editor's local library.

The originals were retained. This is a standalone copy with a new Git repository, not a linked worktree, symlink farm or rewritten book history. Existing archived relative paths were preserved. New research should continue here.

Historical scripts may write beside themselves, hash the old book checkout, or rely on absolute paths. Historical `preexisting_hashes.json`, `preservation.json`, Git state snapshots and file counts describe the ORIGINAL checkout at the date of that pass; they are not an integrity manifest for this smaller repository. Do not rerun those preservation scripts here and interpret missing book files as lost research. Use `tools/verify_archive.py` and the import manifest instead.

Recorded failures, partial runs and corrected variants are preserved intentionally. The smoke runner selects three completed recent scripts and copies their whole directory into temporary storage before execution. Other historical runs may require a temporary copy of more of this repository. Requirements are inherited from the book; environment versions recorded during consolidation are observational, not a new claim that every historical script has been reproduced.

Only selected book context is included. Cross-references from the copied chapter to other book chapters, and old absolute links to the source checkout or library, are historical references and may not resolve on another machine. The research reports are the self-contained starting point. The original book remains the source for book-wide context and rendering.

No source commit, branch, index or working file was changed by the transfer. No initial commit or remote was created in the new repository; Git history can be started when the editor requests it. The files are present locally and the copied evidence is protected by the recorded hashes.
