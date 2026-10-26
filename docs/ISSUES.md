# Issues found during development (cause → fix)

Every problem met while building and validating busco2summary on real data is
recorded here so that other people / agents can follow the reasoning.
Real data used: D1 = 199 annelid genomes, BUSCO 6.0.0 `euk_genome_min` + miniprot,
metazoa_odb12; D2 = 83 runs, BUSCO 5.7.1/5.8.2, genome (miniprot) + transcriptome
(metaeuk), metazoa_odb12 + metazoa_odb10.

| # | Symptom | Cause | Fix |
|---|---------|-------|-----|
| 1 | busco_2_matrix v0.3 `--scan-only` reports 0 BUSCOs for every lineage | key `n_busco` read, parser returns `n_total` | new loader reads `n_total`, `lineage_n_buscos` from txt / json |
| 2 | v0.3 status matrix can miss BUSCOs missing (M) or duplicated (D) in every sample | gene universe = BUSCOs with S/F sequence files | universe = all IDs in `full_table.tsv` (M rows included) |
| 3 | `busco2summary --help` crashed: `unsupported format character ','` | `%` in an argparse help string | escape `%%` |
| 4 | a BUSCO v3 folder holding many `run_<sample>` dirs was read as ONE sample with many runs | v4+ (`run_<lineage>`) and v3 (`run_<sample>`) layouts both start with `run_` | v4+ run folders must contain `_odb`; a parent whose run_* children are v3 runs is a container (test `test_discovery_and_layouts`) |
| 5 | fragment ratio median ≈ 0.35 when measured against lineage proteins | miniprot `Length` in full_table ≠ full protein length; lineage refseq proteins are longer (median 336 vs 237 aa for complete calls) | default reference = median length of complete calls of the same BUSCO in the dataset (`--frag-ref dataset`); lineage reference optional |
| 6 | miniprot runs (BUSCO ≥ 5.4 `--miniprot`, BUSCO 6 default) have no `.fna` | BUSCO writes only `.faa` + `.gff` for miniprot | CDS rebuilt from genome + miniprot alignment (`--genome-map` / `--genome-from-summary`) |
| 7 | CDS from GFF CDS features translated wrongly for ~35 % of BUSCOs | miniprot frameshifts: GFF phase does not continue the frame | CDS rebuilt from the `##PAF` CIGAR (`cg:Z:`) instead of GFF |
| 8 | CIGAR-built CDS one codon short of BUSCO protein | miniprot `##STA` (= BUSCO `.faa`) also translates the codon following the alignment | append the next codon (all 172,218 non-`G` miniprot alignments of the Helobdella genome then translate identically) |
| 9 | CIGAR-built CDS one residue too long at frameshifts | `G` (frameshift substitution) emits no residue in miniprot translation | `G` bases skipped (the 11,972 `G`-containing alignments of the same genome then identical: 184,190/184,190) |
| 10 | CDS translation ≠ BUSCO 6 protein when the protein has internal `*` | BUSCO 6 keeps internal stops in `.faa` (BUSCO 5.8.3 source removed all `*`) | CDS stop handling mirrors the protein of the same copy (`clean_cds`) |
| 11 | 1/1293 CDS not identical | BUSCO trimmed the reported start of that hit (1745674 → 1758447) | PAF record matched by query + contig + best overlap; manifest says `ok_miniprot_cigar_translation_differs` |
| 12 | BUSCO 6 `E` (internal stop codon %) ≈ 30 % in genome mode | miniprot alignments over pseudogenes / frameshifts | reported per sample (`pct_E`, `HIGH_E` flag); not used for filtering by default |
| 13 | 7 "genomes" in D1 have 0 BUSCOs | they are mitochondrial genome accessions | state `EMPTY`, excluded from matrices unless `--keep-empty` |
| 14 | D2 mixes genome and transcriptome runs and two lineages | different projects in one folder | dataset report warns; most common lineage analysed unless `--lineage` / `--lineage all` |
| 15 | benchmark: running freshly cloned third-party scripts was blocked by the sandbox | safety policy for downloaded code | those tools are run only after the user approves; installed tools benchmarked directly |
| 16 | auto-lineage sample: the generic run could be chosen | several `short_summary.specific.*` files possible; first one in sort order was taken | specific summaries ranked, non-generic lineage (not eukaryota/bacteria/archaea) wins |
| 17 | MultiQC benchmark score 0/199 at first | MultiQC names samples after the summary file (`short_summary.specific.<lineage>.<sample>`), and all `run_*/short_summary.txt` collapse into one sample called `short` | benchmark parser strips the prefix; the `short` collision is reported as a MultiQC limitation |
| 18 | synthetic fixtures under-reported other tools | fixtures lacked the top-level `short_summary.specific.*` copy and had byte-identical summaries | fixtures now mimic BUSCO v4+ exactly (unique summaries, `generic` / `specific` copies) |
