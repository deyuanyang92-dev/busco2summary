# Benchmark

Two parts:

1. **Accuracy and speed on real data** (`run_benchmark.py`): the truth (status of every BUSCO in
   every sample) is computed by the benchmark script itself from the raw `full_table.tsv`
   files, independently of every tool. Each tool runs under `/usr/bin/time -v`; per-sample
   S/D/F/M counts and single-copy occupancy lists (≥ 50 %, ≥ 90 %, ≥ 3 taxa) and S+F ≥ 90 %
   are scored against the truth. All tools get the same input (one folder of BUSCO sample
   folders, one lineage).
2. **Robustness** (`robustness.py`): one synthetic folder containing every layout met in practice
   (see `tests/fixtures.py`); a tool scores a case when it reports that sample with the correct
   number of single-copy BUSCOs.

Datasets (real, not redistributed):

| id | samples | BUSCO | mode / predictor | lineage |
|---|---|---|---|---|
| D1 | 199 annelid genomes (7 are mitochondrial-only, all BUSCOs missing) | 6.0.0 | euk_genome_min / miniprot | metazoa_odb12 |
| D2 | 77 runs (genomes + transcriptomes) | 5.8.2 | miniprot + metaeuk | metazoa_odb12 |

Machine: Linux, 128 cores, Python 3.8 (`busco2summary -j 8`). Times are wall-clock for one run
with a warm file-system cache.

## Results

### 1. Real data: accuracy (scored against independently computed truth) and speed

| tool | D1: per-sample counts | D1: S ≥50 % / ≥90 % / ≥3 taxa | D1: S+F ≥90 % | D1 time / RAM | D2: per-sample counts | D2: S ≥50 % / S+F ≥90 % | D2 time |
|---|---|---|---|---|---|---|---|
| **busco2summary `run`** (everything: tables, matrices, lists, plots, HTML) | 199/199 S,D,F,M | 660 / 0 / 672 exact | 4 exact | 18.8 s / 365 MB | 77/77 | 415 / 171 exact | 8.2 s |
| busco2summary `summary` only | 199/199 | – | – | 6.8 s / 57 MB | | | |
| busco2summary re-run (parse cache) | | | | 2.7 s / 56 MB | | | 1.2 s |
| busco_2_matrix v0.3 (our previous script) | 199/199 | 660 / 0 / 672 exact | 4 exact | 24.8 s / 234 MB | 77/77 | 415 / 171 exact | 9.9 s |
| count_buscos.py | S only, 199/199 | 660 / 0 / 672 exact | not supported | 1.1 s / 18 MB | S only, 77/77 | 415 / – | 0.4 s |
| busco2summary_smkv2 (summary) | 199/199 (from %) | – | – | 0.1 s | 77/77 | – | 0.1 s |
| MultiQC busco | 196/199 (3 renamed: `.fasta…` stripped from names) + 1 spurious sample `short` | – | – | 29.2 s / 320 MB | 77/77 | – | 13.8 s |

All tools that produce them agree with the truth on these two well-formed datasets; the
differences are in **what** is produced (S+F lists, D/F/M per gene, flags, extraction, CDS) and in
robustness (below). D1 S ≥ 90 % = 0 because 7 of 199 samples are mitochondrial-only; busco2summary
drops such EMPTY samples by default (`--keep-empty` was used here for a like-for-like comparison) and
then gives S ≥ 90 % = 1, **S+F ≥ 90 % = 90**, S+D+F ≥ 90 % = 414 BUSCOs.

### 2. Robustness: one folder with every awkward layout

| tool | exit | v5 genome | v6 miniprot, dir `*.fna` | transcriptome table | v3 (nested) | auto-lineage | empty sample | no full_table | compleasm | score |
|---|---|---|---|---|---|---|---|---|---|---|
| **busco2summary** | 0 | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | **8/8** |
| busco_2_matrix v0.3 | 0 | ✓ | ✓ | ✓ | missed | ✓ | ✓ | ✓ | missed | 6/8 |
| busco2summary_smkv2 | 0 | ✓ | ✓ | ✓ | missed | wrong run (generic) | ✓ | ✓ | missed | 5/8 |
| MultiQC busco | 0 | merged with next sample | merged | ✓ | ✓ | ✓ | ✓ | ✓ | missed | 5/8 |
| count_buscos.py | **1 (crash)** | – | – | – | – | – | – | – | – | 0/8 |

count_buscos.py stops with `FileNotFoundError` as soon as one sample has no
`single_copy_busco_sequences` folder; MultiQC's name cleaning merges `sampleA` and
`sampleB.fna` into one sample.

### 3. Sequence extraction checks (busco2summary only)

* duplicated BUSCOs, transcriptomes (D2): 300 randomly chosen `--dup best` picks all equal the copy
  with the highest full_table score;
* miniprot runs (D1, 2 genomes, 1,293 CDS rebuilt from genome + miniprot CIGAR): 1,292 translate
  exactly to the BUSCO protein, 1 is flagged `translation_differs` (BUSCO trimmed that hit);
  over the whole miniprot output of one genome, 184,190/184,190 alignments rebuild exactly.

### 4. Gene matrix cross-validation (`matrix_benchmark.py`)

Every tool turns the BUSCO folder into one FASTA file per selected BUSCO (one record per sample).
Truth is computed by the benchmark script itself from the raw `full_table.tsv` and
`busco_sequences/` (no tool code imported). Each output is scored on

* **gene set**: against busco2summary's rule (present in ≥ ⌈t·N⌉ samples) *and* against the rule the
  tool documents (`count/N ≥ t`, `round(N·t)`, ≥ 3 taxa), which separates "different rule" from
  "wrong implementation";
* **records**: every (gene, sample) record has the right status, no expected record is missing;
* **sequences**: each record is identical to that sample's BUSCO protein (`*` ignored);
* **inputs**: md5 of every input file before and after all runs.

Alignment / trimming / tree programs are replaced by pass-through stubs, so only selection and
extraction are timed. Input: D1 = 199 annelid genome runs (192 non-empty, BUSCO 6.0.0 miniprot),
D2 = 77 runs (BUSCO 5.8.2, metazoa_odb12); 672 BUSCOs.

Tools: busco2summary 0.2.0; busco_2_matrix v0.3 (our legacy script); BUSCO_phylogenomics
(jamiemcg, bioconda 20240919); `busco2fasta.py` (lstevens17/busco2phylo-nf, commit 2308a1c);
`busco-phylo.py` (bigey/BUSCO-phylogenetics-pipeline, commit eeb3f12); ZF `BUSCO_extraction.sh`
main and version2 (xtmtd/Phylogenomics, commit 04838a1). Each was run with the status set and
threshold it supports.

**Gene count / records; "exact" = identical to truth (gene set, every record's status, every
sequence).** Wall time in seconds.

| data | set | t | truth: genes / records | busco2summary | v0.3 | BUSCO_phylogenomics | busco2fasta | bigey busco-phylo | ZF version2 | ZF main |
|---|---|---|---|---|---|---|---|---|---|---|
| D1 | S | 50 % | 660 / 96,383 | exact · 38.8 | exact · 16.5 | exact · 12.1¹ | exact · 2.3 | exact · 10.6 | | |
| D1 | S | 90 % | 0 / 0 | exact | exact | exact | exact | exact | | |
| D1 | SF | 50 % | 667 / 107,640 | exact · 40.9 | exact · 17.2 | | | | | |
| D1 | SF | 90 % | 4 / 726 | exact | exact | | | | | |
| D1 | S | ≥ 3 taxa | 672 / 97,221 | exact · 42.6 | | | | | **0 genes**² · 47.5 | |
| D1 | SF | ≥ 3 taxa | 672 / 107,975 | exact · 47.1 | | | | | | **0 genes**² · 303.7 |
| D2 | S | 50 % | 415 / 19,393 | exact · 9.6 | exact · 5.5 | exact · 4.7¹ | **436** (+21)³ · 0.6 | exact · 4.7 | | |
| D2 | S | 90 % | 0 / 0 | exact | exact | exact | exact | exact | | |
| D2 | SF | 50 % | 652 / 40,560 | exact · 16.3 | exact · 7.2 | | | | | |
| D2 | SF | 90 % | 171 / 12,342 | exact · 7.4 | exact · 5.7 | | | | | |
| D2 | S | ≥ 3 taxa | 672 / 27,150 | exact · 15.2 | | | | | 669 genes, 32/77 samples⁴ · 54.4 | |
| D2 | SF | ≥ 3 taxa | 672 / 41,216 | exact · 18.9 | | | | | | 672 genes, 32/77 samples⁴ · 132.1 |

Across all 34 runs: 0 records with a wrong status, 0 sequences different from the BUSCO protein,
0 unresolved headers; every tool left its input unchanged (md5 of every file before/after).

1. exits 1 after writing the sequences: with stub "alignments" its PHYLIP writer fails on unequal
   lengths; with the real MUSCLE the same D1 run did not finish in 4 h (ISSUES #32).
2. D1 is BUSCO 6 miniprot output, which has no `.fna`; ZF collects a locus only when the `.fna`
   exists, so nothing is collected — silently, exit 0.
3. not an error: busco2fasta documents `count ≥ round(N·p)`; N = 77, p = 0.5 → `round(38.5)` = 38,
   while "50 %" read as ⌈N·t⌉ (busco2summary) or `count/N ≥ t` needs 39. Its 436 genes equal the
   truth under its own rule. The same "50 %" therefore means different matrices in different tools.
4. ZF silently dropped the 45 D2 samples without `.fna` (exit 0, no message); 12,675 (S) and
   24,324 (SF) expected records are missing. Re-scored against the truth for the 32 samples it
   kept, both ZF versions are exact (gene sets, 0 missing records, sequences identical).

Conclusion: every tool that selects by occupancy agrees with the independent truth once its own
rounding rule is taken into account; the differences are (a) rounding of the threshold
(busco2fasta) and (b) samples silently skipped (ZF without `.fna`). busco2summary covers S and
S+F, any threshold and the ≥ 3-taxa rule, uses every sample, and reports which samples / calls
were not used (manifest).

## Feature comparison

From the source code / documentation of each tool (✓ yes, ~ partial, ✗ no).

| feature | busco2summary | busco_2_matrix v0.3 | count_buscos.py | busco2fasta.py | ZF BUSCO_extraction.sh | busco2summary_smkv2 | extract-busco-scores.py | MultiQC busco |
|---|---|---|---|---|---|---|---|---|
| per-sample S/D/F/M counts | ✓ | ✓ | S only | ✗ | ✗ | % only | % only | ✓ |
| E (internal stops), assembly stats, versions, predictor | ✓ | versions only | ✗ | ✗ | ✗ | versions only | ✗ | ✗ |
| full_table vs summary consistency check, QC flags | ✓ | ✗ | ✗ | ✗ | ✗ | ✗ | ✗ | ✗ |
| gene × sample status matrix (S/D/F/M, incl. all-missing BUSCOs) | ✓ | ~ | S Y/N | ✗ | ✗ | ✗ | ✗ | ✗ |
| per-gene n_S n_D n_F n_M, occupancy per status set | ✓ | ~ | S only | ✗ | ✗ | ✗ | ✗ | ✗ |
| fragmented BUSCOs used (S+F sets) and fragment length ratio | ✓ | ~ (no ratio) | ✗ | ✗ | ~ S+F (F merged into S, no length ratio) in main `scripts/` BUSCO v5 branch; S only in `version2/` and v3 branch | ✗ | ✗ | ✗ |
| occupancy gene-ID lists, many thresholds in one run | ✓ (exact ceil rule) | ~ (via extraction) | counts only | one threshold, `round()` | one fixed cut (see note) | ✗ | ✗ | ✗ |
| completeness × occupancy grid, taxon-drop advice | ✓ | ~ (Pareto on one axis) | ✗ | ✗ | ✗ | ✗ | ✗ | ✗ |
| duplicated BUSCOs: best-scoring copy | ✓ | ✗ (excluded) | ✗ | ✗ | ✗ | ✗ | ✗ | ✗ |
| nucleotide CDS for miniprot runs (no `.fna`) | ✓ (verified) | ✗ | ✗ | ✗ | ✗ | ✗ | ✗ | ✗ |
| BUSCO v3 / v4-v6 / compleasm / auto-lineage | ✓ / ✓ / ✓ / ✓ | ~ / ✓ / ✗ / ✓ | flag / ✓ / ✗ / ~ | ✗ / ✓ / ✗ / ✗ | ✓ / v5.4 / ✗ / ✗ | ~ / ✓ / ✗ / ✗ | ✗ / ✓ / ✗ / ~ | ✓ / ✓ / ✗ / ~ |
| mixed-lineage guard | ✓ | warn | ✗ | ✗ (last lineage wins) | ✗ | ✗ | ✗ | ✗ |
| does not delete / modify existing data | ✓ | `--overwrite` deletes out dir | ✓ | always deletes out dir | copies input, edits copies | ✓ | ✓ | `-f` overwrites |
| interactive HTML | ✓ (offline) | ✓ | ✗ | ✗ | ✗ | ✗ | ✗ | ✓ (summary only) |
| resumable (parse cache) | ✓ | ✗ | ✗ | ✗ | ✗ | ✗ | ✗ | ✗ |

Third-party tools compared: BUSCO_phylogenomics `count_buscos.py` (bioconda 20240919,
J. McGowan), `busco2fasta.py` from lstevens17/busco2phylo-nf (commit 2308a1c),
`BUSCO_extraction.sh` from xtmtd/Phylogenomics (Feng Zhang lab, Nanjing Agricultural
University; commit 04838a1; two versions, see note below), `extract-busco-scores.py` from bigey/BUSCO-phylogenetics-pipeline
(commit eeb3f12), MultiQC 1.21 busco module; plus our own earlier scripts in `legacy/`.

**Note on ZF `BUSCO_extraction.sh`** (read from the source, commit 04838a1). The main
`scripts/BUSCO_extraction.sh` (2021 ZF, 2023-06-28 DSY), BUSCO v5 branch, copies both
`single_copy_busco_sequences/` and `fragmented_busco_sequences/` and moves the fragmented files into the
single-copy folder, i.e. it uses **S + F** (our `SF` set). `scripts/version2/BUSCO_extraction.sh`
(2023-10-11 ZF) copies only `single_copy_busco_sequences/` (S only), and the BUSCO v3 branch of both uses
S only. A locus is collected for a species only when its `.fna` file exists (so BUSCO runs without
`.fna`, e.g. miniprot, give no loci). Loci filter in both: `grep -o "^$LOCI" 1-raw_loci/log.txt` counts
the *missing* records ("<locus> in <species> does not exist") and a locus is kept when missing
`-lt TOTAL_TAXA - 2`, i.e. present in ≥ 3 taxa, as the comment says. The occupancy cut proper is
applied later in `matrix_generation.sh` (`$4 >= OCCUPANCY` on PhyKIT alignment occupancy).

## Reproduce

```bash
# tools.json: paths of the tools (see run_benchmark.py for keys)
python3 benchmark/run_benchmark.py --dataset D1=/path/busco_dir@metazoa_odb12 \
    --work /tmp/bench --out results --tools-json tools.json
python3 benchmark/robustness.py tools.json results/robustness.tsv
# gene-matrix cross-validation (bench.json: work dir, datasets, tool paths, stub_dir)
python3 benchmark/matrix_benchmark.py run   --config bench.json
python3 benchmark/matrix_benchmark.py score --config bench.json --out results
```
