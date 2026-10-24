# Import provenance and portability

The source repository is `/home/gustav/Git/affine_toda`. The exact source commit, branch, tracked-file status and hash of every imported file are in `import-manifest.json`.

All `briefs/referee_sl21*` files and all material under `referee_scratch/sl21*` were copied, together with the superalgebra background briefs, their scripts, the relevant book chapter, introduction, bibliography and `code/super_toda/`. Python bytecode/cache directories were excluded. The Bernard–LeClair PDF referenced by the reports was also copied from the editor's local library.

The originals were retained during the initial copy. After the consolidated archive and updated project scope were saved in initial commit `213dbc7`, the editor requested removal of the old research copies. On 24 October 2026, all 444 files marked `research` in the import manifest were verified byte for byte against both this checkout and that commit, then deleted from the old checkout. Five generated Python cache files and the 27 emptied project scratch directories were also removed. None of the deleted research files was tracked by the book repository. The book, shared background material and unrelated files were left in place. The [cleanup record](research-cleanup-20261024.json) lists the deleted paths and the preserving commit.

This repository has independent Git metadata and preserves the archived relative paths. New research should continue here. Old absolute paths to the deleted reports and scratch files are historical references; use their matching relative paths in this repository.

Historical scripts may write beside themselves, hash the old book checkout, or rely on absolute paths. Historical `preexisting_hashes.json`, `preservation.json`, Git state snapshots and file counts describe the ORIGINAL checkout at the date of that pass; they are not an integrity manifest for this smaller repository. Do not rerun those preservation scripts here and interpret missing book files as lost research. Use `tools/verify_archive.py` and the import manifest instead.

Recorded failures, partial runs and corrected variants are preserved intentionally. The smoke runner selects three completed recent scripts and copies their whole directory into temporary storage before execution. Other historical runs may require a temporary copy of more of this repository. Requirements are inherited from the book; environment versions recorded during consolidation are observational, not a new claim that every historical script has been reproduced.

Only selected book context is included. Cross-references from the copied chapter to other book chapters, and old absolute links to the source checkout or library, are historical references and may not resolve on another machine. The research reports are the self-contained starting point. The original book remains the source for book-wide context and rendering.

The initial copy did not change the source checkout. The later, explicitly requested cleanup removed only the research originals and generated caches described above; it did not change the source Git history or index. The consolidated research is committed locally in `213dbc7`; no remote or public publication has been configured.

The import, transfer and consolidation-check manifests remain unchanged as records of the initial extraction. Their statements about retained originals and source preservation describe that earlier stage, not the state after cleanup. The scope-clarification record and this cleanup record document subsequent events. `tools/verify_archive.py` verifies the imported files in this repository without requiring the deleted originals.
