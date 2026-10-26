"""Per-sample summary table and dataset-level report."""

from __future__ import annotations

from collections import Counter
from typing import Dict, List, Sequence

from .loader import Sample

ASSEMBLY_KEYS = [("Number of scaffolds", "asm_scaffolds"), ("Number of contigs", "asm_contigs"),
                 ("Total length", "asm_total_length"), ("Percent gaps", "asm_pct_gaps"),
                 ("Scaffold N50", "asm_scaffold_N50"), ("Contigs N50", "asm_contig_N50")]

SUMMARY_FIELDS = [
    "sample", "dir_name", "state", "tool", "busco_version", "lineage", "lineage_date",
    "mode", "gene_predictor", "n_total",
    "C", "S", "D", "F", "M", "pct_C", "pct_S", "pct_D", "pct_F", "pct_M",
    "pct_SF", "n_internal_stop", "pct_E", "avg_identity", "one_line",
    "n_faa_S", "n_faa_D", "n_faa_F", "n_fna", "n_gff",
    "asm_scaffolds", "asm_contigs", "asm_total_length", "asm_pct_gaps",
    "asm_scaffold_N50", "asm_contig_N50",
    "counts_source", "input_file", "path", "run_dir", "other_runs", "notes",
]


def summary_row(s: Sample) -> dict:
    i = s.info
    r: Dict[str, object] = {"sample": s.name, "dir_name": s.dir_name, "state": s.state,
                            "tool": s.tool, "path": s.path, "run_dir": s.run_dir,
                            "other_runs": ",".join(s.other_runs), "notes": " | ".join(s.notes)}
    for k in ("busco_version", "lineage", "lineage_date", "mode", "gene_predictor",
              "n_internal_stop", "pct_E", "avg_identity", "one_line", "input_file"):
        r[k] = i.get(k)
    if s.calls:
        n = s.counts()
        tot = sum(n.values())
        r["counts_source"] = "full_table"
    else:
        n = {"S": i.get("n_single"), "D": i.get("n_duplicated"),
             "F": i.get("n_fragmented"), "M": i.get("n_missing")}
        tot = i.get("n_total")
        r["counts_source"] = "short_summary" if tot else "none"
    r["n_total"] = tot
    r.update(S=n["S"], D=n["D"], F=n["F"], M=n["M"])
    if all(isinstance(n[k], int) for k in "SDFM") and tot:
        r["C"] = n["S"] + n["D"]
        r["pct_C"] = 100 * r["C"] / tot
        for k in "SDFM":
            r[f"pct_{k}"] = 100 * n[k] / tot
        r["pct_SF"] = 100 * (n["S"] + n["F"]) / tot
    asm = i.get("assembly") or {}
    for src, dst in ASSEMBLY_KEYS:
        r[dst] = asm.get(src)
    for code in "SDF":
        r[f"n_faa_{code}"] = sum(1 for v in s.seqs.get(code, {}).values() if "a" in v)
    r["n_fna"] = sum(1 for d in s.seqs.values() for v in d.values() if "n" in v)
    r["n_gff"] = sum(1 for d in s.seqs.values() for v in d.values() if "g" in v)
    return r


def dataset_report(samples: Sequence[Sample], chosen_lineage: str) -> List[str]:
    lines = []
    st = Counter(s.state for s in samples)
    lines.append(f"samples found      : {len(samples)}")
    lines.append("states             : " + ", ".join(f"{k}={v}" for k, v in st.most_common()))
    for key, label in (("lineage", "lineages"), ("busco_version", "BUSCO versions"),
                       ("mode", "modes"), ("gene_predictor", "gene predictors")):
        c = Counter(str(s.info.get(key) or "NA") for s in samples)
        lines.append(f"{label:<19}: " + ", ".join(f"{k}={v}" for k, v in c.most_common()))
    lines.append(f"lineage analysed   : {chosen_lineage}")
    multi = [s.name for s in samples if s.other_runs]
    if multi:
        lines.append(f"samples with >1 run: {len(multi)} (first: {', '.join(multi[:5])})")
    c = Counter(str(s.info.get('busco_version') or 'NA').split('.')[0] for s in samples)
    if len(c) > 1:
        lines.append("WARNING: mixed BUSCO major versions; gene predictors and scores differ.")
    m = Counter(str(s.info.get('mode') or 'NA') for s in samples)
    if len(m) > 1:
        lines.append("WARNING: mixed BUSCO modes (genome/transcriptome/protein): "
                     "D and F rates are not directly comparable.")
    return lines
