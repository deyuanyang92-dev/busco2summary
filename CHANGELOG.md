# Changelog

Versions below 1.0 are development versions (used by the author only); 1.0 will be the first formal release.


## 0.4.0 (2026-10-28)

* new `busco` command: run BUSCO on many assemblies and summarise. Finds SPAdes (`contigs.fasta`,
  `--fasta scaffolds.fasta`), rnaSPAdes (`transcripts.fasta`), Trinity (`Trinity.fasta` or
  `<out>.Trinity.fasta`), MEGAHIT, Flye folders and loose FASTA files. Sample names follow
  `busco2summary_smkv2.py -bfsd --suffix_subdir`: folder name minus the assembler suffix
  (`JLJC1_Eun1_1201_spades_rna` -> `JLJC1_Eun1_1201`); results in `OUT/<sample>`; mode per sample
  (transcriptome for rnaSPAdes/Trinity); finished samples skipped (resume); `-J` parallel runs x `-c`
  threads, no Snakemake; `-n` dry run writes `busco_inputs.tsv`, editable and readable with `--list`;
  failures in `busco_failed.tsv`; BUSCO found on PATH or in a conda env without activation.
* new `list` command: names of all BUSCO result folders under `-i`, one per line (`--long`: folder
  name, layout, runs, path); reads no tables (83 folders in 0.4 s).
* `--samples FILE` (alias of `--include`) for every command: lines may be sample names, folder names
  or paths (the `list --long` output works as is); names not found are reported.
* result folders sharing one name (BUSCO run inside each assembly folder, e.g. `sampleA/spades/busco_scaffolds`) are named after the first parent folder that differs (`sampleA`) instead of `busco_scaffolds`, `busco_scaffolds_2`, ...
* releases are made locally with `scripts/release.sh` (no GitHub Actions workflow scope needed).
* README reorganised.

## 0.3.0 (2026-10-27)

* new `matrix` command: one gene × sample matrix (status / binary / copy number) or a long
  `sample busco_id status` table, to stdout or a file; occupancy (`-t`), `--min-taxa`,
  `--genes` and the usual sample filters; no output folder, cache, log or plots.
  Same gene selection as `genelists` (tested).
* single-file zero-install build `b2s.pyz` (`scripts/build_pyz.py`, 46 KB, standard library only),
  attached to GitHub releases by `.github/workflows/release.yml`.
* without matplotlib: one warning, no empty `plots/` folder.
* CI also tests a lean install (no matplotlib / openpyxl) and the `.pyz`.
* `b2s ... | head` no longer prints a BrokenPipe traceback.
* benchmark: gene-matrix cross-validation against 6 scripts from 5 projects (`benchmark/README.md` §4).

## 0.2.0 (2026-10-27)
Logic audit for use in publications (independent code review + BUSCO source survey + real-data check).
See `docs/VALIDATION.md` and `docs/ISSUES.md` #19-#30.

Fixes that can change reported numbers
* sample without full_table: M = n - S - D - F (was 0, so C was 100 %)
* `--min-frag-ratio`: an F without a reference length no longer counts as present; the same rule now
  drives gene lists, binary matrices, `occ_*`, `n_F_qualified`, new `pct_SF_qualified` and extraction
* unique sample names in all cases; a sample listed twice in a sample sheet is used once
* several runs per sample matching `-l`: newest odb version chosen
* parse cache invalidated by `missing_busco_list*.tsv`, compleasm `summary.txt` and package updates
* compleasm 0.2.x: `Interspaced` -> F, `Retrocopy` -> D, per-lineage `summary.txt` block, folders without `summary.txt`
* denominator `n` from the one-line summary (BUSCO's own), not `Total BUSCO groups searched`

Other
* list names never collide (`occ066p60`); `--min-completeness` accepts fractions like `-t`
* taxon-drop columns renamed `genes_occNNN` (same naming as the gene lists)
* README: how genes are organised (sample completeness vs gene occupancy vs matrix fill, rules, comparison)
* docs/BUSCO_FORMATS.md: output formats of BUSCO 3.0.2-6.1.0 and compleasm from their source code
* corrected description of ZF `BUSCO_extraction.sh` (main version uses S+F)

## 0.1.0 (2026-10-26)
First release of the consolidated package. Replaces and supersedes:

* `busco_2_matrix_v0.3.py` (internal v7.2.0; matrices, Pareto HTML, extraction) - kept in `legacy/`
* `busco2summary_smkv2.py` v2.8.6 (batch BUSCO runner + summary) - kept in `legacy/` (the runner part)
* `busco2summaryv2.1b.py` (summary)

New compared with the scripts above
* one package / one command (`busco2summary`, alias `b2s`), Python >= 3.8, no required dependencies
* layouts: BUSCO v3, v4-v6 (genome / transcriptome / protein, auto-lineage, batch mode), compleasm
* status matrix built from `full_table.tsv` (all BUSCOs incl. M), copy-number / length /
  fragment-ratio / binary (S, SD, SF, SDF) matrices, tidy long table
* per-gene statistics: n/% S D F M, occupancy for each status set, F gain, fragment length
  ratio, mean D copies, flags F_PRONE / D_PRONE / F_RESCUED_50
* per-sample statistics + QC flags (LOW_C, HIGH_D, HIGH_F, HIGH_M, HIGH_E, EMPTY, MISMATCH, MULTI_RUN)
* occupancy gene-ID lists for any thresholds x status sets, `>= K taxa` lists,
  sample-completeness x gene-occupancy grid, greedy taxon-drop advice
* `--min-frag-ratio`: count only fragments covering >= R of the reference length
* extraction: best-scoring / longest copy for duplicates, per-gene and per-sample FASTA,
  header templates, stop-codon handling, exact CDS rebuilding for miniprot runs (no `.fna`)
  with per-sequence translation verification
* self-contained interactive HTML report, static plots (PNG/PDF/SVG), optional XLSX
* parse cache for fast re-runs, run log, parameter record, parallel loading
