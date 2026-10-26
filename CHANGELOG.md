# Changelog

## 1.0.0 (2026-10-26)
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
