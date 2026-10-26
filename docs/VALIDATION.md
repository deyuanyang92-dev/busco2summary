# Validation (v1.1.0, 2026-10-27)

What was checked before using busco2summary numbers in a publication, what was found, and what is
**not** covered. Scripts: `benchmark/validation/`. Issue log: `docs/ISSUES.md`.

## 1. Output formats of every BUSCO version (from the source code)

The output-writing code of BUSCO 2.0.1–6.1.0 (gitlab.com/ezlab/busco, all tags) and compleasm 0.2.10
was read; the facts a parser must handle are in [`BUSCO_FORMATS.md`](BUSCO_FORMATS.md) with file /
function citations. The ones that matter for counting:

* one status per BUSCO ID per run, but **several rows per ID** are possible (duplicates, and since
  miniprot also Complete / Fragmented) → count distinct IDs, never rows;
* S, D, F in short_summary and full_table come from the same sets; M = n − S − D − F, with n = number of
  BUSCOs in the cutoff files (the one-line `n:`), which since 5.3.0 can differ from the
  "Total BUSCO groups searched" line;
* column layout: v3 `Contig Start End`, v4 `Gene Start Gene End`, `Strand` since 5.0, transcriptome /
  protein tables without coordinates, optional OrthoDB columns that may be absent on some rows;
* compleasm: `Single / Duplicated / Fragmented / Interspaced / Missing / Retrocopy`, one `summary.txt`
  block per lineage.

busco2summary reads full_table columns by header name (v3 and v4+ names), counts distinct IDs with
priority S > D > F > M, and takes n from the one-line summary.

## 2. Independent code review

A reviewer who saw only the code (not the author's conclusions) looked for logic errors that change
reported numbers. Seven were found and reproduced, plus four minor ones; all are fixed with a
regression test each (`tests/test_regressions.py`, ISSUES #20–#30).

## 3. Real data: every BUSCO run on the development machine

3,041 samples (3,045 runs), recounted by `check_run.py`, a separate script that does not import
busco2summary:

| BUSCO | mode / predictor | samples |
|---|---|---|
| 5.6.1 | genome / augustus, metaeuk | 6 |
| 5.7.1 | genome / miniprot; transcriptome | 355 |
| 5.8.2 | genome / miniprot; transcriptome (euk, prok); proteins | 2,090 |
| 5.8.3 | genome / miniprot; transcriptome; proteins | 224 |
| 6.0.0 | genome / miniprot; transcriptome; proteins | 364 |

(full list: `benchmark/validation/real_runs_modes.txt`)

* raw full_table distinct-ID counts = short_summary counts in 3,043/3,043 runs that have a summary;
  no BUSCO ID had two different statuses in any run;
* busco2summary S, D, F, M and n = independent recount for **3,041/3,041** samples (before and after
  the fixes);
* occupancy gene-ID lists (S, SD, SF, SDF × 50, 60, 70, 75, 80, 85, 90, 95, 100 %), recomputed by
  `check_lists.py` with the rule "present in ≥ ⌈t·N⌉ samples": **72/72 lists identical** on
  D1 (192 annelid genomes, BUSCO 6.0.0) and D2 (77 runs, BUSCO 5.8.2, metazoa_odb12);
* sequence extraction (earlier benchmark): duplicated BUSCOs, 300/300 `--dup best` picks equal the
  highest-scoring full_table copy; miniprot CDS rebuilt from the genome: 1,292/1,293 translate exactly
  to the BUSCO protein (the one exception is flagged in the manifest).

* gene matrices cross-validated against other tools (`benchmark/matrix_benchmark.py`, results in
  `benchmark/README.md` §4): busco2summary, busco_2_matrix v0.3 and BUSCO_phylogenomics give gene
  sets, records and sequences identical to an independently computed truth on D1 and D2
  (S and SF, 50 % / 90 % / ≥ 3 taxa; up to 107,975 records per matrix, 0 differences).

## 4. Not covered — say so if you use these

* **BUSCO v3, v4.x, v5.0–5.5 and 6.1**: no real runs available; supported from the source-code format
  survey and synthetic test folders only.
* **compleasm**: synthetic folders only.
* augustus / metaeuk genome mode: 6 real samples only; prokaryote modes: 2 samples.
* `--frag-ref lineage` and `--min-frag-ratio`: logic tested, not benchmarked against an external
  tool. With BUSCO ≥ 6.1 + odb12.2, full_table Length means target protein length, so fragment
  ratios are not meaningful there.

## 5. Suggested methods text

> BUSCO results were summarised with busco2summary v1.1.0 (https://github.com/deyuanyang92-dev/busco2summary).
> For each sample, every BUSCO was assigned one status (complete single-copy S, duplicated D,
> fragmented F, missing M) from the distinct BUSCO IDs in `full_table.tsv`. Samples with
> no BUSCO [and samples with S+D+F < X %] were excluded. A BUSCO was retained when it was present
> [as S / as S or F] in at least ⌈t × N⌉ of the N retained samples (t = …). The resulting matrix
> contained … loci with … % of cells filled.

Chinese version:

> 使用 busco2summary v1.1.0 汇总 BUSCO 结果。每个样品中每个 BUSCO 依据 `full_table.tsv` 中的不同 ID 判定唯一状态
> （完整单拷贝 S、重复 D、片段 F、缺失 M）。剔除无 BUSCO 的样品［及 S+D+F < X % 的样品］后，保留在 N 个样品中
> 至少 ⌈t × N⌉ 个样品［以 S／以 S 或 F］出现的 BUSCO（t = …），最终矩阵含 … 个位点，填充率 … %。
