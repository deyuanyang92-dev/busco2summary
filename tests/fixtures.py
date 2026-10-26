"""Synthetic BUSCO / compleasm result folders covering every supported layout."""

from __future__ import annotations

from pathlib import Path
from typing import Dict, Tuple

GENES = [f"{i}at33208" for i in (101, 102, 103, 104, 105, 106)]
PROT = "MKTAYIAKQRQISFVKSHFSRQ"


def _summary_txt(version: str, lineage: str, counts: Dict[str, int], mode: str = "euk_genome_min",
                 extra: str = "", sample: str = "x") -> str:
    n = sum(counts.values())
    pct = {k: 100 * v / n for k, v in counts.items()}
    c = counts["S"] + counts["D"]
    return (f"# BUSCO version is: {version} \n"
            f"# The lineage dataset is: {lineage} (Creation date: 2025-02-25, number of genomes: 9, "
            f"number of BUSCOs: {n})\n"
            f"# Summarized benchmarking in BUSCO notation for file /data/{sample}.fna\n"
            f"# BUSCO was run in mode: {mode}\n{extra}\n"
            f"\t***** Results: *****\n\n"
            f"\tC:{pct['S'] + pct['D']:.1f}%[S:{pct['S']:.1f}%,D:{pct['D']:.1f}%],F:{pct['F']:.1f}%,"
            f"M:{pct['M']:.1f}%,n:{n}\t\n"
            f"\t{c}\tComplete BUSCOs (C)\t\n\t{counts['S']}\tComplete and single-copy BUSCOs (S)\t\n"
            f"\t{counts['D']}\tComplete and duplicated BUSCOs (D)\t\n"
            f"\t{counts['F']}\tFragmented BUSCOs (F)\t\n\t{counts['M']}\tMissing BUSCOs (M)\t\n"
            f"\t{n}\tTotal BUSCO groups searched\t\n")


def _write_run(run: Path, status: Dict[str, str], version: str = "5.8.2", lineage: str = "metazoa_odb12",
               genome_style: bool = True, fna: bool = True, v3: bool = False,
               dup_scores: Tuple[float, float] = (100.0, 250.0), summary: bool = True,
               full_table: bool = True) -> Dict[str, int]:
    """status: gene -> S/D/F/M. Writes full table, summary and sequences."""
    run.mkdir(parents=True, exist_ok=True)
    counts = {k: sum(1 for v in status.values() if v == k) for k in "SDFM"}
    if genome_style:
        header = "# Busco id\tStatus\tSequence\tGene Start\tGene End\tStrand\tScore\tLength\n"
    else:
        header = "# Busco id\tStatus\tSequence\tScore\tLength\n"
    rows = []
    seqdir = run if v3 else run / "busco_sequences"
    dirs = {"S": "single_copy_busco_sequences", "D": "multi_copy_busco_sequences",
            "F": "fragmented_busco_sequences"}
    for d in dirs.values():
        (seqdir / d).mkdir(parents=True, exist_ok=True)
    for i, (g, st) in enumerate(sorted(status.items())):
        if st == "M":
            rows.append(f"{g}\tMissing\n")
            continue
        word = {"S": "Complete", "D": "Duplicated", "F": "Fragmented"}[st]
        copies = 2 if st == "D" else 1
        faa, fna_recs = [], []
        for c in range(copies):
            ctg, a, b = f"ctg{i}_{c}", 100 + c, 400 + c
            score = dup_scores[c] if st == "D" else 300.0
            length = 22 if st != "F" else 10
            if genome_style:
                rows.append(f"{g}\t{word}\t{ctg}\t{a}\t{b}\t+\t{score}\t{length}\n")
                loc = f"{ctg}:{a}-{b}"
            else:
                rows.append(f"{g}\t{word}\t{ctg}:{a}-{b}\t{score}\t{length}\n")
                loc = f"{ctg}:{a}-{b}"
            seq = PROT[: length] + ("W" * c)          # copy 1 is longer
            faa.append(f">{g}_1_0:abc|{loc}|+\n{seq}\n")
            fna_recs.append(f">{g}_1_0:abc|{loc}|+\n{'ATG' * len(seq)}\n")
        (seqdir / dirs[st] / f"{g}.faa").write_text("".join(faa))
        if fna:
            (seqdir / dirs[st] / f"{g}.fna").write_text("".join(fna_recs))
    if full_table:
        name = f"full_table_{run.name[4:]}.tsv" if v3 else "full_table.tsv"
        (run / name).write_text(f"# BUSCO version is: {version} \n# The lineage dataset is: {lineage}\n"
                                + header + "".join(rows))
    if summary:
        name = f"short_summary_{run.name[4:]}.txt" if v3 else "short_summary.txt"
        sname = run.name[4:] if v3 else run.parent.name
        (run / name).write_text(_summary_txt(version, lineage, counts, sample=sname))
        if not v3:   # BUSCO v4+ also writes a copy next to the run folder
            kind = "generic" if lineage.split("_")[0] in ("eukaryota", "bacteria", "archaea") \
                else "specific"
            (run.parent / f"short_summary.{kind}.{lineage}.{run.parent.name}.txt").write_text(
                _summary_txt(version, lineage, counts, sample=sname))
    return counts


STATUS_A = dict(zip(GENES, "SSDFMS"))
STATUS_B = dict(zip(GENES, "SFSSDM"))
STATUS_C = dict(zip(GENES, "FSSSSS"))


def build(root: Path) -> Dict[str, Dict[str, int]]:
    """Create a mixed results folder; returns expected counts per sample."""
    exp = {}
    # BUSCO v5 genome (metaeuk style, with .fna)
    exp["sampleA"] = _write_run(root / "sampleA" / "run_metazoa_odb12", STATUS_A)
    # BUSCO v6 miniprot (no .fna), json summary at top level
    exp["sampleB"] = _write_run(root / "sampleB.fna" / "run_metazoa_odb12", STATUS_B,
                                version="6.0.0", fna=False)
    (root / "sampleB.fna" / "short_summary.specific.metazoa_odb12.sampleB.json").write_text(
        '{"results": {"one_line_summary": "C:50.0%[S:33.3%,D:16.7%],F:16.7%,M:16.7%,n:6,E:20.0%",'
        '"internal_stop_codon_count": 1, "internal_stop_codon_percent": 20.0},'
        '"lineage_dataset": {"name": "metazoa_odb12"}, "versions": {"busco": "6.0.0"},'
        '"parameters": {"mode": "euk_genome_min", "gene_predictor": "miniprot"}}')
    # transcriptome style full table (5 columns)
    exp["tranC"] = _write_run(root / "tranC" / "run_metazoa_odb12", STATUS_C, genome_style=False)
    # BUSCO v3 layout
    exp["v3D"] = _write_run(root / "v3" / "run_v3D", STATUS_A, version="3.0.2", v3=True)
    # auto-lineage: generic + specific runs; specific named by top-level summary
    _write_run(root / "autoE" / "run_eukaryota_odb10", dict(zip(GENES, "MMMMMS")), lineage="eukaryota_odb10")
    exp["autoE"] = _write_run(root / "autoE" / "run_metazoa_odb12", STATUS_B)
    (root / "autoE" / "short_summary.specific.metazoa_odb12.autoE.txt").write_text(
        _summary_txt("5.8.2", "metazoa_odb12", exp["autoE"], sample="autoE"))
    # empty sample (mitochondrial genome)
    exp["emptyF"] = _write_run(root / "emptyF" / "run_metazoa_odb12", dict(zip(GENES, "MMMMMM")))
    # no full table: statuses inferred from sequences + missing list
    run = root / "noftG" / "run_metazoa_odb12"
    exp["noftG"] = _write_run(run, STATUS_C, full_table=False)
    # unrelated folders that must be ignored
    (root / "busco_downloads" / "lineages").mkdir(parents=True, exist_ok=True)
    (root / "notes").mkdir(exist_ok=True)
    (root / "notes" / "readme.txt").write_text("not busco")
    return exp


def build_compleasm(root: Path) -> None:
    d = root / "cmpH"
    (d / "metazoa_odb10").mkdir(parents=True, exist_ok=True)
    (d / "summary.txt").write_text("## lineage: metazoa_odb10\nS:50.00%, 3\nD:16.67%, 1\n"
                                   "F:0.00%, 0\nI:16.67%, 1\nM:16.67%, 1\nN:6\n")
    rows = ["Gene\tStatus\tSequence\tGene Start\tGene End\tStrand\tScore\tLength\n"]
    for g, st in zip(GENES, ["Single", "Single", "Single", "Duplicated", "Interspersed", "Missing"]):
        if st == "Missing":
            rows.append(f"{g}\tMissing\n")
        else:
            rows.append(f"{g}\t{st}\tchr1\t10\t100\t+\t50\t30\n")
    (d / "metazoa_odb10" / "full_table.tsv").write_text("".join(rows))
