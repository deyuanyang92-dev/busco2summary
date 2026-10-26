# BUSCO / compleasm output formats across versions (evidence report)

Sources were read only, never executed. BUSCO came from `git clone https://gitlab.com/ezlab/busco` (tags 2.0.1 … 6.1.0) locally. compleasm came from `github.com/huangnengCSU/compleasm` (HEAD a275e9f = v0.2.10) locally.

Citation key: `T:path::func` means function `func` at tag T. `H` = `src/busco/busco_tools/hmmer.py` (it was `src/busco/BuscoTools.py` in 4.1.x). `R` = `src/busco/BuscoRunner.py`. `CL` = `CHANGELOG` at tag 6.1.0. A function hash showing identical AST between tags is noted as "same AST".

## 1. full_table.tsv

**Comment lines** (`H::write_output_header`; 4.0.0 `src/busco/BuscoAnalysis.py` L526–536; 3.0.2 `BuscoAnalysis.py` L301–306):

| Versions | Comment lines |
|---|---|
| 3.0.2 | `# BUSCO version is: X ` / `# The lineage dataset is: L (Creation date: D, number of species: N, number of BUSCOs: B)` / `# To reproduce this run: <cmd>` / `#` |
| 4.0.0–5.0.0 | same first 2 lines ("number of species"); no reproduce line |
| 5.1.0–6.1.0 | "number of species" → **"number of genomes"** |

Note the trailing space after the version string. The lineage name is the dataset basename, for example `eukaryota_odb10`.

**Column header** (the last `#` line):

| Versions | genome (any predictor) | transcriptome / proteins |
|---|---|---|
| 3.0.2 | `# Busco id  Status  Contig  Start  End  Score  Length` (`GenomeAnalysis.py` L1184) | `# Busco id  Status  Sequence  Score  Length` (L941) |
| 4.0.0–4.1.4 | `# Busco id  Status  Sequence  Gene Start  Gene End  Score  Length [OrthoDB url  Description]` | `… Sequence  Score  Length [OrthoDB url  Description]` |
| 5.0.0–6.1.0 | `# Busco id  Status  Sequence  Gene Start  Gene End  Strand  Score  Length [OrthoDB url  Description]` | unchanged |

- The header depends only on the base mode (`genome`/`transcriptome`/`proteins`, `analysis/*Analysis.py::_mode`), not on the predictor. Augustus, metaeuk, miniprot and prodigal all get the same header.
- The two extra columns are added when `links_to_ODB10.txt` exists. A row whose BUSCO ID is missing from that file is written without them (KeyError fallback in `_format_output_lines`, every version 4.1.4–6.1.0). **This means the column count can vary under one header.**
- `Missing` rows have exactly 2 fields: `ID\tMissing`.
- Sort order:
  - v3: lexicographic.
  - 4.1.4–6.1.0: numeric on the integer before `at` (`H::_sort_lines`, same AST 4.1.4–6.0.0).
- Since 5.6.0, minus-strand genes have Gene Start > Gene End (`sorted(..., reverse=strand=="-")`, 5.6.0 `_format_output_lines`). CL 5.4.7 also changed how negative strand is reported.

**Sequence column**:
- **Genome mode**: the HMMER target is reduced to `contig:low-high` (`H::get_gene_id`). The `:low-high` suffix is then stripped, so the column holds the contig name.
  - 5.5.0 only: the stripping is commented out, so the suffix stays.
  - 5.1.0–5.7.x prokaryotes: `rsplit("_",1)[0]` cut real contig names. A contig called `contig_2` became `contig` (GitLab issue #748, reported on 5.6.0). This was fixed in 5.8.0 (CL 5.8.0 "Fix bug in full_table.tsv output (Issue #748)", and the `else` branch was removed in 5.8.0 `_format_output_lines`).
- **euk transcriptome** (metaeuk, 5.x+): uses the transcriptome branch with no stripping, so the column holds `transcript:low-high`. This is inferred from the code and not verified on real output.

## 2. Status strings, rows vs IDs, how counts are made

- The only BUSCO statuses are `Complete`, `Duplicated`, `Fragmented` and `Missing`.
  - Sources: 3.0.2 `BuscoAnalysis.py` L665–725; 4.1.4–6.1.0 `H::create_output_content` (same AST 5.0.0–6.1.0); `_list_missing_buscos`.
  - In full_table, single-copy is **`Complete`**. It is not "Single".
- **An ID appears under only one status per run.** The three dicts are mutually exclusive:
  - Fragment entries are dropped when the ID already has a complete or very-large match (`H::_remove_lower_ranked_duplicates`).
  - Very-large entries are dropped when the ID is already complete.
  - `consolidate_busco_lists` puts each ID into exactly one of `single_copy_buscos`, `multi_copy_buscos` or `fragmented_buscos`.
  - v3 follows the same logic (`has_complete_match`, L599–628).
- **One ID can still have several rows with the same status:**
  - `Duplicated`: one row per gene copy (expected).
  - `Complete` or `Fragmented` can also repeat. `_format_output_lines` writes one row per *match record* (`for m, match in enumerate(match_info)`). `_sort_matches` appends one record per HMMER target for the same gene coordinates. The code comment says this is "needed for miniprot pipeline, where identical gene coords can match different reference sequences" (6.1.0 L459–463). Rerun duplicates are also tolerated (comment in `_remove_remaining_duplicate_matches`).
  - Fragmented keeps only the best-scoring *gene* (`consolidate_busco_lists`), but that gene can carry more than one record.
  - **Parsers must count distinct IDs per status, never rows.**
- **Where counts come from.** The short_summary, JSON and full_table all come from the same in-memory dicts:
  - S = `len(single_copy_buscos)`, D = `len(multi_copy)`, F = `len(fragmented)`.
  - **M = n − S − D − F**, and **n = `len(cutoff_dict)`**, i.e. the IDs in `scores_cutoff`/`lengths_cutoff` (`H::_get_busco_percentages`).
  - Given this, distinct-ID counts per status in full_table should match the summary.
- **Known mismatch sources:**
  1. Since 5.3.0, `Total BUSCO groups searched` prints `dataset_nb_buscos` from `dataset.cfg` (5.3.0 `produce_hmmer_summary`; 6.1.0 L1259). Up to 5.2.2 it printed `total_buscos`. The `n:` value always uses `len(cutoff_dict)`, so **`n` and the "Total" line can differ** when dataset.cfg and the cutoff files disagree.
  2. Percentages are rounded independently: `C%` is rounded from (S+D), not computed as S%+D%. In ≤5.3.x, `M%` = 100−S%−D%−F%.
  3. 6.1.0 + miniprot: when a BUSCO has several complete copies, those with internal stop codons are dropped. If only one stop-free copy remains, the BUSCO becomes S instead of D (`consolidate_busco_lists`, 6.1.0 only). full_table follows the same dicts.
  - I found no public bug report of short_summary counts disagreeing with full_table distinct IDs. **Unverified/none found** (GitLab search only).

## 3. short_summary

**One-line summary:**

| Versions | One-line format |
|---|---|
| 3.0.2–5.5.0 | `C:x%[S:x%,D:x%],F:x%,M:x%,n:N` |
| 5.6.0–6.1.0 | `…,n:N` plus `,E:x%` **only when E>0** (miniprot only; E = complete BUSCOs with an internal stop codon ÷ (S+D)) (`H::record_results`, `_get_stop_codon_percent_and_avg_identity`) |

**Count lines.** Each line in the file starts with a tab (`summary_file.write("\t{}")`):

- `***** Results: *****`
- `<n>\tComplete BUSCOs (C)`
- `Complete and single-copy BUSCOs (S)`
- `Complete and duplicated BUSCOs (D)`
- `Fragmented BUSCOs (F)`
- `Missing BUSCOs (M)`
- `Total BUSCO groups searched`

Annotations that can follow the counts:

- 5.5.0 only: `(incl. K with frameshifts)` can be appended to the C, S, D and F lines.
- 5.6.0+: `(of which K contain internal stop codons)` can be appended to the C line.
- Tab padding differs between versions. **Parse with the regex `^\s*(\d+)\s+(Complete|...)` rather than fixed tabs.**

**Header lines:**
- Same two `#` lines as full_table, then `# Summarized benchmarking in BUSCO notation for file <path>` and `# BUSCO was run in mode: <mode>`.
- From 5.4.0 the mode is the internal name: `euk_genome_met|euk_genome_aug|euk_genome_min|prok_genome_prod|prok_genome_min|euk_tran|prok_tran|proteins` (R `self.mode`, `BuscoConfig.update_mode`). Older versions print `genome`/`transcriptome`/`proteins`.
- 5.0.0+ genome runs add `# Gene predictor used: X`.

**Trailing blocks:**
- 6.1.0 (and, per CL 5.4.0, BBTools from 5.4.0; the exact text is verified only at 6.1.0): `Assembly Statistics:`, `Dependencies and versions:`, and optionally `Placement file versions:`.
- 5.3.0+: a parasitic-contamination line may be appended (R L997 at 6.1.0).

**File names:**

| Versions | Location |
|---|---|
| 3.0.2 | `run_<out>/short_summary_<out>.txt` |
| 4.0.0+ | `<out>/run_<lineage>/short_summary.txt`, copied to `<out>/short_summary.specific.<lineage>.<out>.txt`; `short_summary.generic.<domain_lineage>.<out>.txt` only with auto-lineage (4.0.0 `R` L140–154; 6.1.0 `R::copy_summary_file`, `organize_final_output`) |

**short_summary.json:**
- Introduced in 5.3.0 (CL 5.3.0; absent in 5.2.2), with `.json` copies made the same way as the txt.
- Keys by version:

| Versions | Keys |
|---|---|
| 5.3.x | flat: `one_line_summary, C, S, D, F, M` (counts), `input_file, mode, gene_predictor, dataset, dataset_creation_date, dataset_number_species, dataset_total_buscos` (5.3.0 `H::write_json_summary`) |
| 5.4.0–5.5.0 | nested `parameters`, `lineage_dataset`, `versions`, `results{one_line_summary, Complete, Single copy, Multi copy, Fragmented, Missing, n_markers, domain, +BBTools metrics}`. **Here these result values are percentages** (`R::save_results`). |
| 5.6.0–6.1.0 | `results{one_line_summary, Complete percentage, Complete BUSCOs, Single copy percentage, Single copy BUSCOs, Multi copy percentage, Multi copy BUSCOs, Fragmented percentage, Fragmented BUSCOs, Missing percentage, Missing BUSCOs, n_markers, domain}` |

Optional keys in 5.6.0–6.1.0 `results`:
- `internal_stop_codon_count/percent` (miniprot)
- `avg_identity` (present by 5.8.3; first version not pinned)
- `translation_table` (6.0.0+, prodigal)

## 4–5. Directory layout, sequence files, GFF

- Layout:
  - v3: `run_<out>/` holds full_table, missing_busco_list and `single_copy_busco_sequences/` (genome mode only; `GenomeAnalysis.py` L891–918).
  - 4.0.0+: `<out>/run_<lineage>/{full_table.tsv, missing_busco_list.tsv, short_summary.*, busco_sequences/{single_copy,multi_copy,fragmented}_busco_sequences/}` (`BuscoConfig` "run_{}"; `H` L47–54).
- File names are `<BUSCOid>.faa` and `<BUSCOid>.fna`. **A Duplicated file holds one record per copy.** From 5.6.0 a file holds one record per match record, so a single gene can appear twice (`H::write_buscos_to_file`).
- .fna presence by predictor:

| Predictor | .fna | Evidence |
|---|---|---|
| augustus, prodigal | yes (4.x+) | `write_buscos_to_file(aa, nt)` |
| metaeuk | **no in 5.0.0–5.3.x; yes from 5.4.0** | 5.0.0 `GenomeAnalysis.py` L414 (aa only); CL 5.4.0 |
| miniprot | **never** | 5.5.0 L1033 `#, self.sequences_nt`; 6.1.0 `base.py::record_gene_details` comment "miniprot has no nt_seq"; GitLab #730 (5.7.0) |
| proteins mode | .faa only | `GeneSetAnalysis` |

- From 5.7.0 the default eukaryote genome predictor is **miniprot**. Default eukaryote runs since then therefore have no .fna (CL 5.7.0; `BuscoConfig.update_mode`).
- FASTA headers (taken from gene_id builders):

| Versions | metaeuk | miniprot | augustus | prodigal |
|---|---|---|---|---|
| 3.0.2 | – | – | `>ID:<input path>:<coord>` (L924; exact coord string unverified) | – |
| 4.x–5.3 | `contig:low-high` | – | `contig:start-end` | raw prodigal ID (inferred) |
| 5.4.x–5.5 | `contig:start-end`; on minus strand high-low | 5.5: `Tacc\|contig:s-e` | same | same |
| 5.6.x | `Tacc\|contig:low-high` | same | same | `contig:l-r` |
| 5.7.0–6.1.0 | `Tacc\|contig:low-high\|strand` | same | `contig:s-e\|strand` | `contig:l-r\|strand` |

  - Tacc is the refseq_db target ID; its first `_` field is the BUSCO ID.
  - 5.6.0–6.0.0 miniprot records were built without `description=""`, so Biopython probably writes `>id <unknown description>`. This is inferred from Biopython defaults; augustus ≤5.5 has the same pattern. 6.1.0 sets `description=""`.
- GFF files:
  - `busco_sequences/*/<ID>.gff` for metaeuk from 5.4.0 and for miniprot from 5.7.0 (`GenomeAnalysis::write_gff_files`).
  - Augustus writes `augustus_output/gff/`.
  - Miniprot also writes a full `<input>_<lineage>.gff` in its output folder.

## 6. Length and Score columns

- **Score** is the HMMER full-sequence bitscore of the predicted protein in every version.
- **Length**:
  - odb10 (all versions) and odb12 (6.0.0+): `hmm_matched_length`, the sum of merged HMM-coordinate spans in **profile positions (≈aa)**, not nt.
  - odb12.2 (6.1.0): `tlen`, the **HMMER target (protein) length in aa** (`H::_sort_matches`).
  - It is the same quantity for every predictor, including miniprot.

## 7. Lineage dataset files

| File | Format / use |
|---|---|
| `scores_cutoff` | `ID score` |
| `lengths_cutoff` | odb10: whitespace columns, col3 = sigma, col4 = mean length (col2 unused); odb12.2: `ID median scaledMAD`; odb12: none, falls back to the HMM length (`H::_load_length`) |
| `refseq_db.faa.gz` | gzip FASTA for metaeuk (5.0.0+) and miniprot, decompressed to tmp (`base.py::decompress_refseq_file`) |
| `ancestral`, `ancestral_variants` | FASTA for augustus/blast and the metaeuk rerun (`metaeuk.py` L156) |
| `dataset.cfg` | keys `name, creation_date, number_of_species, number_of_BUSCOs` |
| `links_to_ODB10.txt` | 3 tab columns: ID, description, URL |

Details beyond what the code reads (for example the col-2 meaning) are unverified.

## compleasm (v0.2.10, `compleasm.py`)

- `full_table.tsv` header (no `#`): `Gene Status Sequence Gene Start Gene End Strand Score Length Identity Fraction Frameshift events Best gene Codons` (L1467).
  - Statuses are the enum `.name` values: `Single, Duplicated, Fragmented, Interspaced, Missing, Retrocopy` (L553–572).
- `full_table_busco_format.tsv` has BUSCO-style headers (`# compleasm version is:` from v0.2.10) and `# Busco id … Strand Score Length`. It maps Single→`Complete` and Interspaced→`Fragmented`, but leaves `Retrocopy` unmapped (L1558–1563).
- `summary.txt` is **appended** (`'a'`), one block per lineage:
  - `## lineage: X`, then `S:xx.xx%, k`, `D:…`, optional `R:…` (`--retrocopy`; otherwise R is folded into D), `F:…`, `I:…`, `M:…`, `N:n` (L1680–1698).
  - The `# compleasm version` header line exists only from v0.2.10.
  - Percentages have 2 decimals. **There is no C line**; C = S+D.
  - M includes unaccounted IDs (`d`).
- Protein mode header: `# Busco id Status Sequence Score Length`, with status `Single` (L1980–2027).

## Format facts a parser must handle (highest risk first)

1. Count **distinct IDs per status**. Complete and Fragmented rows can repeat for one ID (one row per match record), and Duplicated always repeats.
2. In compleasm `full_table.tsv`, single-copy is `Single`. compleasm also has `Interspaced` and `Retrocopy`, its summary.txt has no C line and uses 2 decimals, and it holds multiple lineage blocks appended in one file.
3. Row column count varies (2 for Missing; ±2 OrthoDB columns per row). Use the header for names, but don't require every row to match its width.
4. `n:` ≠ `Total BUSCO groups searched` is possible in 5.3.0+. Treat `n` (or the cutoff ID count) as the denominator.
5. The optional `,E:x%` suffix (5.6.0+, miniprot) and the parenthetical annotations on count lines (5.5.0, 5.6.0+) must be tolerated.
6. JSON schemas: flat counts (5.3), nested **percentages** under count-like names (5.4–5.5), and explicit `… BUSCOs`/`… percentage` keys (5.6+).
7. The Strand column exists only from 5.0.0. v3 genome headers use `Contig Start End`. Gene Start > Gene End on the minus strand from 5.6.0.
8. The Sequence column is not always a clean contig name: suffixes remain in 5.5.0 and in eukaryote transcriptome runs, and contig names are truncated in prokaryote runs 5.1–5.7.x (issue #748).
9. Per-version file names (v3 `full_table_<out>.tsv`) and locations (`short_summary.specific.<lineage>.<out>.txt`). Generic summaries exist only with auto-lineage.
10. .fna is absent for miniprot (the default for eukaryotes from 5.7.0), for metaeuk ≤5.3, and in proteins mode. FASTA header syntax changed at 5.6 and 5.7, possibly with a trailing `<unknown description>`.
11. Header wording changed from "number of species" to "number of genomes" at 5.1.0. v3 adds `# To reproduce this run:` and `#` lines.
12. The meaning of Length changes with dataset generation (odb10/12: HMM-matched length; odb12.2: protein length).
