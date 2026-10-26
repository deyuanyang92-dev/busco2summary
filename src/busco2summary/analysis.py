"""Matrices, per-gene / per-sample statistics and occupancy gene lists.

Terminology
-----------
status sets (which statuses count as "present" for a gene in a sample):
    S    complete single-copy only (classic phylogenomics)
    SD   complete (single + duplicated)
    SF   single-copy + fragmented (recovers loci lost to fragmentation)
    SDF  any hit (S, D or F)
occupancy of a gene = fraction of samples in which it is present.
A gene passes threshold t when  present >= ceil(t * N - 1e-9)
(exact integer rule; no rounding surprises at 50 %).
"""

from __future__ import annotations

import math
from dataclasses import dataclass, field
from typing import Dict, List, Optional, Sequence, Tuple

from .loader import Sample

SETS: Dict[str, str] = {"S": "S", "SD": "SD", "SF": "SF", "SDF": "SDF"}
DEFAULT_THRESHOLDS = [0.5, 0.6, 0.7, 0.75, 0.8, 0.85, 0.9, 0.95, 1.0]


def need(t: float, n: int) -> int:
    """Minimum number of samples for occupancy threshold ``t`` (0-1)."""
    return max(0, math.ceil(t * n - 1e-9))


def _median(v: List[float]) -> Optional[float]:
    v = sorted(x for x in v if x is not None)
    if not v:
        return None
    n = len(v)
    return v[n // 2] if n % 2 else (v[n // 2 - 1] + v[n // 2]) / 2


@dataclass
class Dataset:
    """Gene x sample status matrix for one lineage."""

    lineage: str
    samples: List[Sample]
    genes: List[str]
    # gene -> string of status letters, one per sample (same order as samples)
    rows: Dict[str, str]
    ref_len: Dict[str, float] = field(default_factory=dict)
    ref_len_source: str = "none"
    lineage_len: Dict[str, float] = field(default_factory=dict)
    min_frag_ratio: float = 0.0
    # gene -> list of F-length / reference-length ratios (per sample, None if unknown)
    frag_ratio: Dict[str, List[Optional[float]]] = field(default_factory=dict)

    @property
    def names(self) -> List[str]:
        return [s.name for s in self.samples]

    @property
    def n(self) -> int:
        return len(self.samples)

    # ------------------------------------------------------------------
    def present_mask(self, gid: str, sset: str) -> List[bool]:
        row = self.rows[gid]
        fr = self.frag_ratio.get(gid)
        out = []
        for j, c in enumerate(row):
            ok = c in sset
            if ok and c == "F" and self.min_frag_ratio > 0 and fr is not None:
                r = fr[j]
                ok = r is None or r >= self.min_frag_ratio
            out.append(ok)
        return out

    def present_count(self, gid: str, sset: str) -> int:
        if self.min_frag_ratio <= 0 or "F" not in sset:
            row = self.rows[gid]
            return sum(row.count(c) for c in sset)
        return sum(self.present_mask(gid, sset))


def build_dataset(samples: Sequence[Sample], lineage: str,
                  lineage_len: Optional[Dict[str, float]] = None,
                  lineage_src: str = "", min_frag_ratio: float = 0.0,
                  frag_ref: str = "dataset") -> Dataset:
    """Build the gene x sample matrix.

    Fragment ratio = length of a fragmented call / reference length, where the
    reference is the median length of complete calls of that BUSCO in this
    dataset (``frag_ref='dataset'``, same gene predictor, recommended) or the
    lineage dataset proteins (``frag_ref='lineage'``).
    """
    genes = sorted({g for s in samples for g in s.calls})
    rows = {g: "".join(s.status(g) for s in samples) for g in genes}
    ds = Dataset(lineage=lineage, samples=list(samples), genes=genes, rows=rows,
                 min_frag_ratio=min_frag_ratio)
    ds.lineage_len = dict(lineage_len or {})
    for g in genes:
        lens = [s.calls[g][2] for s in samples
                if g in s.calls and s.calls[g][0] in "SD" and s.calls[g][2]]
        m = _median(lens)
        if m:
            ds.ref_len[g] = m
    ds.ref_len_source = "median complete length in dataset"
    if frag_ref == "lineage" and ds.lineage_len:
        ds.ref_len = dict(ds.lineage_len)
        ds.ref_len_source = lineage_src or "lineage dataset"
    for g in genes:
        ref = ds.ref_len.get(g)
        fr: List[Optional[float]] = []
        for s in samples:
            c = s.calls.get(g)
            if c and c[0] == "F" and c[2] and ref:
                fr.append(min(c[2] / ref, 9.99))
            else:
                fr.append(None)
        ds.frag_ratio[g] = fr
    return ds


# --------------------------------------------------------------------------
# per-gene statistics
# --------------------------------------------------------------------------

GENE_FIELDS = ["busco_id", "n_samples", "n_S", "n_D", "n_F", "n_M",
               "pct_S", "pct_D", "pct_F", "pct_M",
               "occ_S", "occ_SD", "occ_SF", "occ_SDF", "F_gain",
               "lineage_ref_length", "median_len_complete", "median_len_F",
               "median_F_ratio", "n_F_qualified", "mean_copies_D", "flags"]


def gene_stats(ds: Dataset, high_f: float = 0.3, high_d: float = 0.2) -> List[dict]:
    out = []
    n = ds.n or 1
    for g in ds.genes:
        row = ds.rows[g]
        c = {k: row.count(k) for k in "SDFM"}
        lens_c, lens_f, dcopies = [], [], []
        for s in ds.samples:
            call = s.calls.get(g)
            if not call:
                continue
            if call[0] in "SD" and call[2]:
                lens_c.append(call[2])
            if call[0] == "F" and call[2]:
                lens_f.append(call[2])
            if call[0] == "D":
                dcopies.append(call[1])
        fr = [r for r in ds.frag_ratio.get(g, []) if r is not None]
        nq = sum(1 for r in fr if r >= ds.min_frag_ratio) if ds.min_frag_ratio > 0 else c["F"]
        occ = {k: ds.present_count(g, k) / n for k in SETS}
        flags = []
        if c["F"] / n >= high_f:
            flags.append("F_PRONE")
        if c["D"] / n >= high_d:
            flags.append("D_PRONE")
        if c["M"] == ds.n:
            flags.append("ABSENT_ALL")
        if occ["S"] < 0.5 <= occ["SF"]:
            flags.append("F_RESCUED_50")
        out.append({
            "busco_id": g, "n_samples": ds.n,
            "n_S": c["S"], "n_D": c["D"], "n_F": c["F"], "n_M": c["M"],
            "pct_S": c["S"] / n, "pct_D": c["D"] / n, "pct_F": c["F"] / n, "pct_M": c["M"] / n,
            "occ_S": occ["S"], "occ_SD": occ["SD"], "occ_SF": occ["SF"], "occ_SDF": occ["SDF"],
            "F_gain": occ["SF"] - occ["S"],
            "lineage_ref_length": ds.lineage_len.get(g),
            "median_len_complete": _median(lens_c),
            "median_len_F": _median(lens_f),
            "median_F_ratio": _median(fr),
            "n_F_qualified": nq,
            "mean_copies_D": (sum(dcopies) / len(dcopies)) if dcopies else None,
            "flags": ",".join(flags),
        })
    return out


# --------------------------------------------------------------------------
# per-sample statistics
# --------------------------------------------------------------------------

SAMPLE_FIELDS = ["sample", "state", "n_genes", "S", "D", "F", "M",
                 "pct_S", "pct_D", "pct_C", "pct_F", "pct_M", "pct_SF", "pct_SDF",
                 "mean_F_ratio", "n_F_qualified", "rank_SF", "qc_flags"]


@dataclass
class QCThresholds:
    min_c: float = 0.5      # LOW_C when C% below
    max_d: float = 0.10     # HIGH_D when D% above
    max_f: float = 0.15     # HIGH_F when F% above
    max_m: float = 0.30     # HIGH_M when M% above
    max_e: float = 0.50     # HIGH_E (internal stop codons, BUSCO >= 5.8) above


def sample_stats(ds: Dataset, qc: Optional[QCThresholds] = None) -> List[dict]:
    qc = qc or QCThresholds()
    rows = []
    ng = len(ds.genes) or 1
    for j, s in enumerate(ds.samples):
        col = [ds.rows[g][j] for g in ds.genes]
        c = {k: col.count(k) for k in "SDFM"}
        frs = [ds.frag_ratio[g][j] for g in ds.genes if ds.frag_ratio[g][j] is not None]
        nq = (sum(1 for r in frs if r >= ds.min_frag_ratio)
              if ds.min_frag_ratio > 0 else c["F"])
        flags = []
        if s.state != "OK":
            flags.append(s.state)
        pc = (c["S"] + c["D"]) / ng
        if pc < qc.min_c:
            flags.append("LOW_C")
        if c["D"] / ng > qc.max_d:
            flags.append("HIGH_D")
        if c["F"] / ng > qc.max_f:
            flags.append("HIGH_F")
        if c["M"] / ng > qc.max_m:
            flags.append("HIGH_M")
        pe = s.info.get("pct_E")
        if isinstance(pe, (int, float)) and pe / 100 > qc.max_e:
            flags.append("HIGH_E")
        if s.other_runs:
            flags.append("MULTI_RUN")
        rows.append({
            "sample": s.name, "state": s.state, "n_genes": ng,
            "S": c["S"], "D": c["D"], "F": c["F"], "M": c["M"],
            "pct_S": c["S"] / ng, "pct_D": c["D"] / ng, "pct_C": pc,
            "pct_F": c["F"] / ng, "pct_M": c["M"] / ng,
            "pct_SF": (c["S"] + c["F"]) / ng, "pct_SDF": (ng - c["M"]) / ng,
            "mean_F_ratio": (sum(frs) / len(frs)) if frs else None,
            "n_F_qualified": nq, "rank_SF": 0, "qc_flags": ",".join(flags),
        })
    order = sorted(range(len(rows)), key=lambda i: -rows[i]["pct_SF"])
    for r, i in enumerate(order, 1):
        rows[i]["rank_SF"] = r
    return rows


# --------------------------------------------------------------------------
# occupancy tables / gene lists
# --------------------------------------------------------------------------

def occupancy_lists(ds: Dataset, thresholds: Sequence[float],
                    sets: Sequence[str] = ("S", "SD", "SF", "SDF"),
                    min_taxa: Sequence[int] = ()) -> Tuple[List[dict], Dict[str, List[str]]]:
    """Gene lists per (set, threshold) and a summary table.

    Returns (table_rows, {"<set>_occ<NN>": [ids], "<set>_min<K>": [ids]}).
    """
    counts = {k: {g: ds.present_count(g, k) for g in ds.genes} for k in sets}
    lists: Dict[str, List[str]] = {}
    table = []
    n = ds.n

    def row_for(label, kind, value, k, keep):
        cells = sum(counts[k][g] for g in keep)
        fill = cells / (len(keep) * n) if keep and n else 0.0
        return {"set": k, "rule": kind, "value": value, "min_samples": label,
                "n_genes": len(keep), "matrix_fill": fill}

    for t in thresholds:
        nmin = need(t, n)
        for k in sets:
            keep = [g for g in ds.genes if counts[k][g] >= nmin and counts[k][g] > 0]
            lists[f"{k}_occ{int(round(t * 100)):03d}"] = keep
            table.append(row_for(nmin, "occupancy", t, k, keep))
    for m in min_taxa:
        for k in sets:
            keep = [g for g in ds.genes if counts[k][g] >= m]
            lists[f"{k}_min{m}"] = keep
            table.append(row_for(m, "min_taxa", m, k, keep))
    # F gain columns (SF vs S, SDF vs SD) for the same rule
    by = {(r["set"], r["rule"], r["value"]): r["n_genes"] for r in table}
    for r in table:
        base = {"SF": "S", "SDF": "SD"}.get(r["set"])
        r["F_gain_genes"] = (r["n_genes"] - by[(base, r["rule"], r["value"])]
                             if base and (base, r["rule"], r["value"]) in by else "")
    return table, lists


def occupancy_curve(ds: Dataset, sset: str) -> List[int]:
    """n_genes passing for every integer sample count 1..N (for plots / HTML)."""
    hist = [0] * (ds.n + 1)
    for g in ds.genes:
        hist[ds.present_count(g, sset)] += 1
    out, acc = [0] * (ds.n + 1), 0
    for k in range(ds.n, -1, -1):
        acc += hist[k]
        out[k] = acc
    return out


# --------------------------------------------------------------------------
# sample-completeness x gene-occupancy grid and greedy taxon dropping
# --------------------------------------------------------------------------

def _bitrows(ds: Dataset, sset: str) -> List[int]:
    rows = []
    for g in ds.genes:
        m = 0
        for j, ok in enumerate(ds.present_mask(g, sset)):
            if ok:
                m |= 1 << j
        rows.append(m)
    return rows


def _pop(x: int) -> int:
    return bin(x).count("1")


def grid(ds: Dataset, sset: str, sample_cuts: Sequence[float],
         gene_cuts: Sequence[float]) -> List[dict]:
    bits = _bitrows(ds, sset)
    ng = len(ds.genes) or 1
    comp = [sum(1 for b in bits if b >> j & 1) / ng for j in range(ds.n)]
    out = []
    for sc in sample_cuts:
        keep = [j for j in range(ds.n) if comp[j] >= sc]
        mask = sum(1 << j for j in keep)
        n = len(keep)
        cnt = [_pop(b & mask) for b in bits]
        for gc in gene_cuts:
            nmin = need(gc, n)
            sel = [c for c in cnt if c >= nmin and c > 0]
            cells = sum(sel)
            out.append({"set": sset, "min_sample_completeness": sc, "n_taxa": n,
                        "min_gene_occupancy": gc, "n_genes": len(sel),
                        "matrix_fill": cells / (len(sel) * n) if sel and n else 0.0,
                        "present_cells": cells})
    return out


def taxon_drop(ds: Dataset, sset: str, thresholds: Sequence[float],
               steps: int = 10, candidates: int = 25) -> List[dict]:
    """Greedy: repeatedly drop the sample whose removal adds most genes.

    The objective is the number of genes passing the *first* threshold of
    ``thresholds``; numbers for all thresholds are reported at every step.
    Only the ``candidates`` least-complete samples are tried per step.
    """
    bits = _bitrows(ds, sset)
    active = set(range(ds.n))
    ng = len(ds.genes) or 1
    comp = {j: sum(1 for b in bits if b >> j & 1) / ng for j in range(ds.n)}
    cnt = [_pop(b) for b in bits]

    def passing(cnts, n, t):
        nmin = need(t, n)
        return sum(1 for c in cnts if c >= nmin and c > 0)

    out = []
    n = len(active)
    out.append({"step": 0, "dropped": "", "sample_completeness": "", "n_taxa": n,
                **{f"genes_occ{int(round(t * 100))}": passing(cnt, n, t) for t in thresholds}})
    t0 = thresholds[0]
    for step in range(1, steps + 1):
        if n <= 3:
            break
        cand = sorted(active, key=lambda j: comp[j])[:candidates]
        best, best_val = None, -1
        nmin = need(t0, n - 1)
        for j in cand:
            v = 0
            for b, c in zip(bits, cnt):
                cc = c - (b >> j & 1)
                if cc >= nmin and cc > 0:
                    v += 1
            if v > best_val:
                best, best_val = j, v
        if best is None:
            break
        active.discard(best)
        cnt = [c - (b >> best & 1) for b, c in zip(bits, cnt)]
        n -= 1
        out.append({"step": step, "dropped": ds.samples[best].name,
                    "sample_completeness": round(comp[best], 4), "n_taxa": n,
                    **{f"genes_occ{int(round(t * 100))}": passing(cnt, n, t) for t in thresholds}})
    return out
