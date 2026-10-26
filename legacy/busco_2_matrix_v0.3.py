#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
╔══════════════════════════════════════════════════════════════════════════════╗
║         BUSCO Phylogenomics Gene Set Extractor  v7.2                        ║
║         ─────────────────────────────────────────────────────────            ║
║  Publication-ready tool for extracting BUSCO orthologs for phylogenomics.   ║
║  Compatible with BUSCO v3 / v4 / v5 / v6 output directory layouts.          ║
╚══════════════════════════════════════════════════════════════════════════════╝

WHAT THIS SCRIPT DOES
─────────────────────
1. Reads BUSCO result directories for many samples (hundreds to thousands).
2. Builds a Gene × Sample status matrix  (S / D / F / M).
3. Provides interactive Pareto-optimal threshold visualisation so you can
   choose the best gene-coverage cut-off before extraction.
4. Extracts sequences in GPS (one FASTA per gene) and/or SCH (one FASTA
   per sample) formats, ready for MAFFT → trimAl → IQ-TREE.
5. [NEW v7.2] --explore_view (-ev): stand-alone rich visual explorer that
   reads ANY previous output directory (or specific matrix file) and renders
   a comprehensive self-contained HTML dashboard:
     – Multi-tab layout (Overview · Coverage · Completeness · Pareto · Advice)
     – Interactive gene-coverage histogram with threshold sliders
     – Sample completeness ranked bar chart (colour-coded)
     – Full Pareto scatter + dual-axis sweep for every mode
     – Smart recommendation engine with decision rationale
     – One-click download of filtered gene lists

NEW IN v7.2
───────────
• --explore_view / -ev PATH  : load a finished run directory (or a single
  matrix .tsv) and open a rich HTML decision dashboard without re-running
  the full pipeline.  PATH may be:
    – a previous -o output directory  (auto-detects all matrices)
    – a specific busco_status_matrix.tsv or busco_binary_*.tsv file
  The generated file is  <PATH>/explore_view.html  (or alongside the file).
• All v7.1 features preserved unchanged.

Author  : Community enhanced version
Version : 7.2.0
License : MIT
"""

from __future__ import annotations

import argparse
import csv
import json
import os
import re
import shutil
import sys
import textwrap
from dataclasses import dataclass, field
from pathlib import Path
from time import time
from typing import Any, Dict, List, Optional, Set, Tuple

# ──────────────────────────────────────────────────────────────────────────────
# Optional dependencies
# ──────────────────────────────────────────────────────────────────────────────
try:
    import numpy as np
    HAS_NUMPY = True
except ImportError:
    HAS_NUMPY = False

try:
    import pandas as pd
    HAS_PANDAS = True
except ImportError:
    HAS_PANDAS = False

try:
    from tqdm import tqdm
    HAS_TQDM = True
except ImportError:
    HAS_TQDM = False
    def tqdm(iterable, **kwargs):          # type: ignore
        return iterable

# ──────────────────────────────────────────────────────────────────────────────
# Constants
# ──────────────────────────────────────────────────────────────────────────────
VERSION = "7.2.0"

BUSCO_STATUS_ORDER = ["S", "D", "F", "M"]

STATUS_DIRS: Dict[str, str] = {
    "S": "single_copy_busco_sequences",
    "F": "fragmented_busco_sequences",
    "D": "multi_copy_busco_sequences",
}

STATUS_TO_CODE: Dict[str, str] = {
    "Complete": "S",
    "Single": "S",
    "Complete and single-copy": "S",
    "Duplicated": "D",
    "Complete and duplicated": "D",
    "Fragmented": "F",
    "Missing": "M",
}

NUC_SUFFIXES = [".fna", ".ffn"]
DELIMITERS   = ["\t", ",", ";", "|", " "]

# ── Canonical column keys ──────────────────────────────────────────────────────
COL_ORIGINAL = "original_name"
COL_NEW      = "new_name"
COL_ABSPATH  = "abs_path_each_busco_results"
COL_LINEAGE  = "lineage_datasets"

# ── Synonym table (lower-case, '#' stripped → canonical key) ─────────────────
COL_SYNONYMS: Dict[str, str] = {}

def _reg(*aliases: str, canonical: str) -> None:
    for a in aliases:
        COL_SYNONYMS[a.lower()] = canonical

_reg("original_name", "orig", "sample", "sample_name", "input_name",
     "directory", "dir_name", "basename", "old_name", "src",
     canonical=COL_ORIGINAL)

_reg("new_name", "new", "rename", "label", "alias", "taxon", "species",
     "output_name", "display_name", "final_name", "dst", "renamed",
     canonical=COL_NEW)

_reg("abs_path_each_busco_results", "path", "abs_path", "busco_path",
     "busco_dir", "full_path", "abspath", "directory_path",
     "result_path", "result_dir", "busco_result_path", "busco_result",
     canonical=COL_ABSPATH)

_reg("lineage_datasets", "lineage", "lid", "dataset", "odb",
     "lineage_db", "lineage_name", "db", "database", "odb_lineage",
     canonical=COL_LINEAGE)

# ── Human-readable format hint shown in mapping errors ────────────────────────
MAPPING_FORMAT_HINT = """\
  ┌─ Mapping file format ──────────────────────────────────────────────────┐
  │  • Tab/comma/semicolon-separated                                        │
  │  • First non-comment line MUST be a header row                         │
  │  • Every header token MUST start with '#'                              │
  │  • Column ORDER is free – identified by name, not position             │
  │  • Lines starting with '##' are ignored (comments)                     │
  └────────────────────────────────────────────────────────────────────────┘

  Example 1 – Mode 1 (absolute paths + rename):
    #original_name    #new_name          #abs_path_each_busco_results
    GCF_000001405     Homo_sapiens       /data/busco/GCF_000001405
    GCA_900519105     Mus_musculus       /data/busco/GCA_900519105

  Example 2 – Mode 2 (search under -i + rename):
    #original_name    #new_name
    run_GCF001        Homo_sapiens
    run_GCA002        Mus_musculus

  Example 3 – per-sample lineage filter:
    #new_name    #abs_path_each_busco_results    #lineage_datasets
    Sp_A         /busco/SpA                      metazoa_odb10
    Sp_B         /busco/SpB                      insecta_odb10

  Accepted column header synonyms:
    #original_name  → orig, sample, sample_name, directory, basename, src
    #new_name       → new, rename, label, taxon, species, alias, dst, renamed
    #abs_path_…     → path, abs_path, busco_path, busco_dir, result_dir
    #lineage_…      → lineage, lid, odb, dataset, db, lineage_name"""


# ──────────────────────────────────────────────────────────────────────────────
# Data structures
# ──────────────────────────────────────────────────────────────────────────────

@dataclass
class BuscoSeqInfo:
    faa: Optional[Path] = None
    fna: Optional[Path] = None

@dataclass
class SpeciesData:
    name:           str
    original_name:  str = ""
    busco_path:     Optional[Path] = None
    run_dir:        str = ""
    lineage:        str = ""
    busco_version:  str = ""
    buscos:         Dict[str, Dict[str, BuscoSeqInfo]] = field(default_factory=dict)
    status_from_table: Dict[str, str]                  = field(default_factory=dict)
    n_complete:   int = 0
    n_single:     int = 0
    n_duplicated: int = 0
    n_fragmented: int = 0
    n_missing:    int = 0
    n_total:      int = 0

@dataclass
class MappingEntry:
    original_name: str
    new_name:      str
    abs_path:      Optional[Path]
    lineage:       str = ""

SpeciesDict = Dict[str, SpeciesData]


# ──────────────────────────────────────────────────────────────────────────────
# Console helpers
# ──────────────────────────────────────────────────────────────────────────────
_STEP_N = [0]

def step(msg: str) -> None:
    _STEP_N[0] += 1
    bar = "─" * 68
    print(f"\n{bar}\n  STEP {_STEP_N[0]:02d}  │  {msg}\n{bar}")

def substep(msg: str) -> None: print(f"    ▸ {msg}")
def info(msg: str)    -> None: print(f"  [INFO]  {msg}")
def warn(msg: str)    -> None: print(f"  [WARN]  {msg}", file=sys.stderr)

def die(msg: str, code: int = 1) -> None:
    print(f"\n  [ERROR] {msg}", file=sys.stderr)
    sys.exit(code)

def ensure_dir(p: Path) -> Path:
    p.mkdir(parents=True, exist_ok=True)
    return p

def now_ts() -> str:
    import time as _t
    return _t.strftime("%Y-%m-%d %H:%M:%S", _t.localtime())


# ──────────────────────────────────────────────────────────────────────────────
# Name sanitisation
# ──────────────────────────────────────────────────────────────────────────────

def sanitize_name(name: str, replace_dot: bool = True, replace_char: str = "_") -> str:
    name = name.strip()
    name = re.sub(r"\s+", "_", name)
    if replace_dot:
        name = name.replace(".", replace_char)
    name = re.sub(r"[^\w\-]", "_", name)
    name = re.sub(r"_+", "_", name)
    return name.strip("_")

def sanitize_fasta_header(name: str) -> str:
    return re.sub(r"\s+", "_", name.strip())


# ──────────────────────────────────────────────────────────────────────────────
# Lineage helpers
# ──────────────────────────────────────────────────────────────────────────────

def _norm_lineage(s: str) -> str:
    if not s: return ""
    s = re.sub(r"^run_", "", s.strip().lower())
    s = re.sub(r"\.(txt|json)$", "", s)
    return re.sub(r"\s+", "", s)

def lineage_matches(req: str, obs: str) -> bool:
    return _norm_lineage(req) == _norm_lineage(obs)

def extract_lineage_from_dirname(name: str) -> str:
    return name[4:] if name.startswith("run_") else name


# ──────────────────────────────────────────────────────────────────────────────
# Sequence-file presence helpers
# ──────────────────────────────────────────────────────────────────────────────

def _fp(p: Optional[Path]) -> bool:
    return p is not None and p.exists()

def has_protein(obj: Optional[BuscoSeqInfo])    -> bool: return obj is not None and _fp(obj.faa)
def has_nucleotide(obj: Optional[BuscoSeqInfo]) -> bool: return obj is not None and _fp(obj.fna)
def has_any(obj: Optional[BuscoSeqInfo])        -> bool: return has_protein(obj) or has_nucleotide(obj)

def has_seq_for_kind(obj: Optional[BuscoSeqInfo], kind: str) -> bool:
    if kind == "protein":    return has_protein(obj)
    if kind == "nucleotide": return has_nucleotide(obj)
    if kind == "both":       return has_any(obj)
    raise ValueError(f"Unknown seq kind: {kind}")

def _status_order_for_mode(mode: str) -> List[str]:
    return {"S": ["S"], "F": ["F"], "SF": ["S", "F"]}.get(mode, ["S"])

def can_extract(sp: SpeciesData, busco_id: str, mode: str, seqtype: str) -> bool:
    st_map = sp.buscos.get(busco_id)
    if not st_map: return False
    for st in _status_order_for_mode(mode):
        obj = st_map.get(st)
        if obj and has_seq_for_kind(obj, seqtype): return True
    return False

def get_display_status(sp: SpeciesData, busco_id: str) -> str:
    st_map = sp.buscos.get(busco_id)
    if st_map:
        for st in ["S", "F", "D"]:
            if st_map.get(st) and has_any(st_map[st]): return st
    return sp.status_from_table.get(busco_id, "M")


# ──────────────────────────────────────────────────────────────────────────────
# Mapping-file parsing
# ──────────────────────────────────────────────────────────────────────────────

def _detect_delimiter(line: str) -> str:
    for d in ["\t", ",", ";", "|"]:
        if line.count(d) >= 1: return d
    return "\t"

def _match_col_token(raw: str) -> Optional[str]:
    return COL_SYNONYMS.get(raw.lstrip("#").strip().lower())


def parse_mapping_file(
    file_path: str,
    input_dir: Optional[Path] = None,
    replace_dot: bool = True,
    replace_char: str = "_",
    global_lineage: str = "",
) -> List[MappingEntry]:
    p = Path(file_path)
    if not p.exists():
        die(f"Mapping file not found: {file_path}\n\n{MAPPING_FORMAT_HINT}")

    substep(f"Reading mapping file: {p}")
    raw_lines  = p.read_text(encoding="utf-8", errors="replace").splitlines()
    non_empty  = [ln for ln in raw_lines if ln.strip() and not ln.strip().startswith("##")]

    if not non_empty:
        die(f"Mapping file is empty or only contains comments: {file_path}\n\n"
            + MAPPING_FORMAT_HINT)

    header_line: Optional[str] = None
    header_idx: int = 0
    for i, ln in enumerate(non_empty):
        if re.match(r"^\s*#[^#]", ln.strip()):
            header_line = ln.strip()
            header_idx  = i
            break

    if header_line is None:
        die(
            "No column header line found.\n"
            "       The first non-comment line must start with '#'.\n"
            "       Example header:  #original_name  #new_name  "
            "#abs_path_each_busco_results\n\n"
            + MAPPING_FORMAT_HINT
        )

    delimiter = _detect_delimiter(header_line)
    parts     = [x.strip() for x in header_line.split(delimiter)]
    col_map: Dict[str, int] = {}
    unrecognised: List[str] = []

    for ci, part in enumerate(parts):
        tok = _match_col_token(part)
        if tok:
            col_map[tok] = ci
        elif part.startswith("#") and part.strip("#").strip():
            unrecognised.append(part)

    substep(f"Delimiter detected      : {repr(delimiter)}")
    substep(f"Recognised columns      : {list(col_map.keys())}")
    if unrecognised:
        warn(f"Unrecognised header tokens (ignored): {unrecognised}")
        warn(f"  → Check spelling.  Accepted synonyms:\n{MAPPING_FORMAT_HINT}")

    has_orig    = COL_ORIGINAL in col_map
    has_new     = COL_NEW      in col_map
    has_abspath = COL_ABSPATH  in col_map
    has_lin_col = COL_LINEAGE  in col_map

    if not has_orig and not has_new and not has_abspath:
        die(
            f"Mapping file header contains no recognisable columns.\n"
            f"       Found tokens: {parts}\n\n"
            + MAPPING_FORMAT_HINT
        )
    if not has_abspath and not input_dir and not has_orig:
        die(
            "Cannot resolve sample paths:\n"
            "       No '#abs_path_each_busco_results' column and no -i given.\n"
            "       Provide absolute paths in the mapping file, OR pass -i.\n\n"
            + MAPPING_FORMAT_HINT
        )

    substep(f"Mapping mode: orig={has_orig}  new={has_new}  "
            f"abspath={has_abspath}  lineage_col={has_lin_col}")

    data_lines = non_empty[header_idx + 1:]
    entries: List[MappingEntry] = []
    seen_new: Dict[str, str]   = {}
    n_skipped = 0

    def gc(row: List[str], token: str) -> str:
        idx = col_map.get(token)
        if idx is None or idx >= len(row): return ""
        return row[idx].strip().strip("\"'")

    for line_num, raw in enumerate(data_lines, start=header_idx + 2):
        ln = raw.strip()
        if not ln or ln.startswith("#"): continue
        row = [x.strip().strip("\"'") for x in ln.split(delimiter)]

        orig    = gc(row, COL_ORIGINAL)
        new_raw = gc(row, COL_NEW)
        abspath = gc(row, COL_ABSPATH)
        lineage = gc(row, COL_LINEAGE) if has_lin_col else ""
        lineage = lineage or global_lineage

        if not orig and not new_raw and not abspath:
            warn(f"Line {line_num}: all recognised fields empty – skipping")
            n_skipped += 1; continue

        if has_new and new_raw:
            final_name = sanitize_name(new_raw, replace_dot, replace_char)
        elif has_orig and orig:
            final_name = sanitize_name(orig, replace_dot, replace_char)
        else:
            tail = (abspath.rstrip("/\\").rsplit("/", 1)[-1]
                    or abspath.rstrip("\\").rsplit("\\", 1)[-1])
            final_name = sanitize_name(tail, replace_dot, replace_char)

        resolved: Optional[Path] = None
        if abspath:
            rp = Path(abspath).resolve()
            if rp.exists():
                resolved = rp
            else:
                warn(f"Line {line_num}: path does not exist: {abspath}. Skipping '{final_name}'.")
                n_skipped += 1; continue
        elif input_dir:
            lookup    = orig or final_name
            candidate = input_dir / lookup
            if candidate.exists():
                resolved = candidate
            else:
                warn(f"Line {line_num}: directory '{lookup}' not found under {input_dir}. Skipping.")
                n_skipped += 1; continue
        else:
            warn(f"Line {line_num}: cannot resolve path for '{final_name}'. Skipping.")
            n_skipped += 1; continue

        if final_name in seen_new:
            die(
                f"Duplicate sanitized sample name: '{final_name}'\n"
                f"       Shared by '{orig}' and '{seen_new[final_name]}'.\n"
                f"       → Rename one entry in #new_name, or use --keep-dots."
            )
        seen_new[final_name] = orig or abspath

        entries.append(MappingEntry(
            original_name=orig or (resolved.name if resolved else ""),
            new_name=final_name,
            abs_path=resolved,
            lineage=lineage,
        ))

    if not entries:
        die(
            f"No valid entries resolved from mapping file: {file_path}\n"
            f"       Skipped {n_skipped} rows due to errors above.\n\n"
            + MAPPING_FORMAT_HINT
        )

    substep(f"Resolved {len(entries)} valid entries  ({n_skipped} skipped)")
    return entries


def load_samples_list(file_path: str) -> Set[str]:
    p = Path(file_path)
    if not p.exists():
        die(f"Samples list not found: {file_path}")
    out: Set[str] = set()
    for ln in p.read_text("utf-8").splitlines():
        ln = ln.strip()
        if ln and not ln.startswith("#"):
            out.add(ln.split(",")[0].strip())
    return out


# ──────────────────────────────────────────────────────────────────────────────
# short_summary & full_table parsing
# ──────────────────────────────────────────────────────────────────────────────

def parse_short_summary(path: Path) -> Dict[str, Any]:
    r: Dict[str, Any] = dict(lineage="", busco_version="",
                              n_complete=0, n_single=0, n_duplicated=0,
                              n_fragmented=0, n_missing=0, n_total=0)
    if not path.exists(): return r
    txt = path.read_text("utf-8", errors="replace")
    for pat, key in [(r"BUSCO version[^:]*:\s*(\S+)", "busco_version"),
                     (r"lineage dataset[^:]*:\s*(\S+)", "lineage")]:
        m = re.search(pat, txt, re.I)
        if m: r[key] = re.sub(r"\s*\(.*", "", m.group(1))
    for pat, key in [
        (r"(\d+)\s+Complete BUSCOs",           "n_complete"),
        (r"(\d+)\s+Complete and single-copy",  "n_single"),
        (r"(\d+)\s+Complete and duplicated",   "n_duplicated"),
        (r"(\d+)\s+Fragmented BUSCOs",         "n_fragmented"),
        (r"(\d+)\s+Missing BUSCOs",            "n_missing"),
        (r"(\d+)\s+Total BUSCO",               "n_total"),
    ]:
        m = re.search(pat, txt, re.I)
        if m: r[key] = int(m.group(1))
    if r["n_total"] == 0:
        m = re.search(r"n:(\d+)", txt)
        if m: r["n_total"] = int(m.group(1))
    return r

def find_short_summary(run_dir: Path) -> Optional[Path]:
    for pat in ["short_summary*.txt", "short_summary.txt"]:
        hits = list(run_dir.glob(pat))
        if hits: return hits[0]
    return None

def parse_full_table(path: Path) -> Dict[str, str]:
    result: Dict[str, str] = {}
    if not path.exists(): return result
    with open(path, "rt", encoding="utf-8", errors="replace") as f:
        for ln in f:
            if ln.startswith("#") or not ln.strip(): continue
            parts = ln.strip().split("\t")
            if len(parts) >= 2:
                result.setdefault(parts[0].strip(),
                                  STATUS_TO_CODE.get(parts[1].strip(), "M"))
    return result

def find_full_table(species_dir: Path, run_dir: Optional[Path] = None) -> Optional[Path]:
    cands: List[Path] = []
    if run_dir:
        cands.append(run_dir / "full_table.tsv")
        cands.extend(run_dir.glob("full_table_*.tsv"))
    cands.extend(species_dir.glob("run_*/full_table.tsv"))
    cands.extend(species_dir.glob("run_*/full_table_*.tsv"))
    cands.append(species_dir / "full_table.tsv")
    return next((c for c in cands if isinstance(c, Path) and c.exists()), None)


# ──────────────────────────────────────────────────────────────────────────────
# run_* directory selection
# ──────────────────────────────────────────────────────────────────────────────

def list_run_dirs(species_dir: Path) -> List[Path]:
    if not species_dir.is_dir(): return []
    return sorted(p for p in species_dir.iterdir()
                  if p.is_dir() and p.name.startswith("run_"))

def score_run_dir(rd: Path) -> Tuple[int, int, int]:
    base = rd / "busco_sequences"
    if not base.is_dir(): return (0, 0, 0)
    counts: Dict[str, int] = {}
    for code, dname in STATUS_DIRS.items():
        d = base / dname
        counts[code] = len({fp.stem for fp in d.iterdir()
                            if fp.is_file() and not fp.name.startswith(".")}) \
                       if d.is_dir() else 0
    return (counts.get("S", 0), counts.get("F", 0), counts.get("D", 0))

def pick_run_dir(species_dir: Path, run_name: str = "",
                 lineage: str = "") -> Optional[Path]:
    run_dirs = list_run_dirs(species_dir)
    if not run_dirs: return None
    if run_name:
        for rd in run_dirs:
            if rd.name == run_name: return rd
        warn(f"run_name '{run_name}' not found in {species_dir.name}")
    if lineage:
        matching = []
        for rd in run_dirs:
            if lineage_matches(lineage, extract_lineage_from_dirname(rd.name)):
                matching.append(rd); continue
            s = find_short_summary(rd)
            if s:
                slg = parse_short_summary(s).get("lineage", "")
                if slg and lineage_matches(lineage, slg):
                    matching.append(rd)
        if len(matching) == 1: return matching[0]
        if len(matching) > 1:
            scored = sorted(((score_run_dir(r), r) for r in matching), reverse=True)
            warn(f"Multiple runs match '{lineage}' in {species_dir.name}; "
                 f"using {scored[0][1].name}")
            return scored[0][1]
    if len(run_dirs) == 1: return run_dirs[0]
    scored = sorted(((score_run_dir(r), r) for r in run_dirs), reverse=True)
    if scored[0][0] != (0, 0, 0):
        warn(f"Auto-selected {scored[0][1].name} for {species_dir.name}")
        return scored[0][1]
    return run_dirs[0]


# ──────────────────────────────────────────────────────────────────────────────
# Sequence scanning
# ──────────────────────────────────────────────────────────────────────────────

def find_seq_file(seq_dir: Path, busco_id: str, stype: str) -> Optional[Path]:
    if stype == "faa":
        p = seq_dir / f"{busco_id}.faa"
        return p if p.exists() else None
    for suf in NUC_SUFFIXES:
        p = seq_dir / f"{busco_id}{suf}"
        if p.exists(): return p
    return None

def scan_busco_sequences(run_dir: Path,
                         include_dup: bool = False) -> Dict[str, Dict[str, BuscoSeqInfo]]:
    result: Dict[str, Dict[str, BuscoSeqInfo]] = {}
    seq_base = run_dir / "busco_sequences"
    if not seq_base.is_dir(): return result
    statuses = ["S", "F"] + (["D"] if include_dup else [])
    for sc in statuses:
        seq_dir = seq_base / STATUS_DIRS[sc]
        if not seq_dir.is_dir(): continue
        ids = {fp.stem for fp in seq_dir.iterdir()
               if fp.is_file() and not fp.name.startswith(".")
               and fp.suffix in [".faa", ".fna", ".ffn"]}
        for bid in ids:
            obj = BuscoSeqInfo(find_seq_file(seq_dir, bid, "faa"),
                               find_seq_file(seq_dir, bid, "fna"))
            if obj.faa or obj.fna:
                result.setdefault(bid, {})[sc] = obj
    return result

def parse_species_data(
    final_name: str, original_name: str,
    species_dir: Path, run_dir: Path,
    include_dup: bool = False, use_full_table: bool = False,
) -> Optional[SpeciesData]:
    if not run_dir or not run_dir.exists(): return None
    d = SpeciesData(name=final_name, original_name=original_name,
                    busco_path=species_dir, run_dir=run_dir.name)
    s = find_short_summary(run_dir)
    if s:
        ps = parse_short_summary(s)
        for k in ["lineage", "busco_version", "n_complete", "n_single",
                  "n_duplicated", "n_fragmented", "n_missing", "n_total"]:
            setattr(d, k, ps[k])
    else:
        d.lineage = extract_lineage_from_dirname(run_dir.name)
    ft = find_full_table(species_dir, run_dir)
    if ft: d.status_from_table = parse_full_table(ft)
    d.buscos = scan_busco_sequences(run_dir, include_dup)
    if use_full_table and d.status_from_table:
        for bid in d.status_from_table:
            d.buscos.setdefault(bid, {})
    if d.n_total == 0:
        d.n_total = max(len(d.status_from_table), len(d.buscos))
    return d


# ──────────────────────────────────────────────────────────────────────────────
# Input scanning
# ──────────────────────────────────────────────────────────────────────────────

def scan_from_entries(
    entries: List[MappingEntry],
    global_run_name: str = "", global_lineage: str = "",
    include_dup: bool = False, use_full_table: bool = False,
    show_progress: bool = True,
) -> Tuple[SpeciesDict, Set[str], Dict[str, str], Dict[str, List[str]]]:
    sd: SpeciesDict = {}; all_ids: Set[str] = set()
    chosen: Dict[str, str] = {}; lg_groups: Dict[str, List[str]] = {}
    it = tqdm(entries, desc="  Scanning BUSCO results", disable=not show_progress)
    for e in it:
        lineage_here = e.lineage or global_lineage
        rd = pick_run_dir(e.abs_path, global_run_name, lineage_here)
        if not rd:
            warn(f"No run_* directory for '{e.new_name}' at {e.abs_path}"); continue
        sp = parse_species_data(e.new_name, e.original_name, e.abs_path,
                                rd, include_dup, use_full_table)
        if not sp: continue
        if lineage_here and sp.lineage and not lineage_matches(lineage_here, sp.lineage):
            warn(f"Skipping '{e.new_name}': lineage mismatch "
                 f"(found={sp.lineage}, requested={lineage_here})"); continue
        sd[e.new_name] = sp; chosen[e.new_name] = rd.name
        all_ids.update(sp.buscos.keys())
        lg = sp.lineage or "unknown"
        lg_groups.setdefault(lg, []).append(e.new_name)
    return sd, all_ids, chosen, lg_groups

def build_directory_entries(input_dir: Path, suffix_filter: str = "",
                             samples_list: Optional[Set[str]] = None
                             ) -> List[MappingEntry]:
    if not input_dir.is_dir():
        die(f"Input directory not found: {input_dir}")
    entries: List[MappingEntry] = []
    for p in sorted(input_dir.iterdir()):
        if not p.is_dir(): continue
        if suffix_filter and not p.name.endswith(suffix_filter): continue
        if samples_list and p.name not in samples_list: continue
        entries.append(MappingEntry(p.name, p.name, p))
    if not entries:
        die(f"No candidate directories found in {input_dir}")
    return entries

def scan_all_lineages(entries: List[MappingEntry]) -> Dict[str, Dict[str, Any]]:
    result: Dict[str, Dict[str, Any]] = {}
    for e in entries:
        for rd in list_run_dirs(e.abs_path):
            s  = find_short_summary(rd)
            ps = parse_short_summary(s) if s else {}
            lg = ps.get("lineage", "") or extract_lineage_from_dirname(rd.name)
            nt = ps.get("n_busco", 0)
            bucket = result.setdefault(lg, {"species": [], "n_busco": 0})
            if e.new_name not in bucket["species"]:
                bucket["species"].append(e.new_name)
            if nt > bucket["n_busco"]:
                bucket["n_busco"] = nt
    return result


# ──────────────────────────────────────────────────────────────────────────────
# Coverage / completeness
# ──────────────────────────────────────────────────────────────────────────────

def gene_coverage(sd: SpeciesDict, busco_id: str, mode: str, seqtype: str) -> float:
    if not sd: return 0.0
    return sum(1 for sp in sd.values()
               if can_extract(sp, busco_id, mode, seqtype)) / len(sd)

def sample_completeness(sp: SpeciesData, all_ids: Set[str],
                        mode: str, seqtype: str) -> float:
    if not all_ids: return 0.0
    return sum(1 for bid in all_ids
               if can_extract(sp, bid, mode, seqtype)) / len(all_ids)


# ──────────────────────────────────────────────────────────────────────────────
# MATRIX BUILDING
# ──────────────────────────────────────────────────────────────────────────────

def build_raw_status_matrix(sd: SpeciesDict,
                             all_ids: Set[str]) -> Dict[str, Dict[str, str]]:
    return {bid: {n: get_display_status(sd[n], bid) for n in sorted(sd)}
            for bid in sorted(all_ids)}

def save_matrices(sd: SpeciesDict, all_ids: Set[str],
                  raw: Dict[str, Dict[str, str]],
                  out_dir: Path, modes: List[str], seqtype: str) -> None:
    names = sorted(sd); bids = sorted(all_ids)
    sf = out_dir / "busco_status_matrix.tsv"
    with open(sf, "w", newline="") as f:
        w = csv.writer(f, delimiter="\t")
        w.writerow(["BUSCO_ID"] + names + ["cov_S", "cov_SF"])
        for bid in bids:
            row = [bid] + [raw[bid].get(n, "M") for n in names]
            row += [f"{gene_coverage(sd, bid, 'S',  seqtype):.4f}",
                    f"{gene_coverage(sd, bid, 'SF', seqtype):.4f}"]
            w.writerow(row)
    info(f"Status matrix (S/D/F/M)  : {sf}")

    for m in modes:
        pf = out_dir / f"busco_binary_{m}.tsv"
        with open(pf, "w", newline="") as f:
            w = csv.writer(f, delimiter="\t")
            w.writerow(["BUSCO_ID"] + names + ["coverage"])
            for bid in bids:
                row = [bid] + [1 if can_extract(sd[n], bid, m, seqtype) else 0
                               for n in names]
                row.append(f"{gene_coverage(sd, bid, m, seqtype):.4f}")
                w.writerow(row)
        info(f"Binary matrix (mode={m})  : {pf}")

    af = out_dir / "sequence_availability.tsv"
    with open(af, "w", newline="") as f:
        w = csv.writer(f, delimiter="\t")
        w.writerow(["BUSCO_ID"] + [f"{n}_prot" for n in names]
                   + [f"{n}_nuc" for n in names])
        for bid in bids:
            row = [bid]
            row += ["1" if any(has_protein(sd[n].buscos.get(bid, {}).get(st))
                               for st in ["S", "F", "D"]) else "0"
                    for n in names]
            row += ["1" if any(has_nucleotide(sd[n].buscos.get(bid, {}).get(st))
                               for st in ["S", "F", "D"]) else "0"
                    for n in names]
            w.writerow(row)
    info(f"Seq availability audit    : {af}")


# ──────────────────────────────────────────────────────────────────────────────
# PARETO OPTIMALITY
# ──────────────────────────────────────────────────────────────────────────────

def evaluate_threshold_sweep(
    sd: SpeciesDict, all_ids: Set[str],
    mode: str, seqtype: str,
    thresholds: Optional[List[float]] = None,
) -> List[Dict[str, Any]]:
    if thresholds is None:
        thresholds = [round(0.50 + 0.05 * i, 3) for i in range(11)]
    cov = {bid: gene_coverage(sd, bid, mode, seqtype) for bid in all_ids}
    n_sp = len(sd)
    rows: List[Dict[str, Any]] = []
    for t in thresholds:
        kept = [bid for bid, c in cov.items() if c >= t]
        n_genes = len(kept)
        avg_comp = (
            sum(sum(1 for bid in kept if can_extract(sp, bid, mode, seqtype))
                / n_genes for sp in sd.values()) / n_sp
            if n_genes else 0.0
        )
        rows.append(dict(threshold=t, n_genes=n_genes,
                         avg_completeness=avg_comp,
                         total_cells=n_genes * n_sp,
                         gene_ids=sorted(kept)))
    return rows

def find_pareto_front(rows: List[Dict[str, Any]],
                      o1: str = "n_genes",
                      o2: str = "avg_completeness") -> List[int]:
    pareto: List[int] = []
    for i, ri in enumerate(rows):
        dominated = any(
            (rj[o1] >= ri[o1] and rj[o2] >= ri[o2]
             and (rj[o1] > ri[o1] or rj[o2] > ri[o2]))
            for j, rj in enumerate(rows) if i != j
        )
        if not dominated: pareto.append(i)
    return pareto

def smart_recommendation(
    sd: SpeciesDict, all_ids: Set[str],
    modes: List[str], seqtype: str, out_dir: Path,
) -> None:
    substep("Running Pareto optimality sweep …")
    thresholds = [round(0.50 + 0.05 * i, 3) for i in range(11)]
    print()
    print("  ╔══════════════════════════════════════════════════════════════╗")
    print("  ║           SMART THRESHOLD RECOMMENDATION (Pareto)           ║")
    print("  ╚══════════════════════════════════════════════════════════════╝")

    for mode in modes:
        rows = evaluate_threshold_sweep(sd, all_ids, mode, seqtype, thresholds)
        pareto_idx = find_pareto_front(rows)
        max_g = max(r["n_genes"] for r in rows) or 1
        max_c = max(r["avg_completeness"] for r in rows) or 1.0
        for r in rows:
            r["pareto"] = rows.index(r) in pareto_idx
            r["score"]  = (r["n_genes"] / max_g) * (r["avg_completeness"] / max_c)

        print(f"\n  Mode = {mode}  ({seqtype})")
        hdr = (f"  {'Threshold':>10}  {'n_genes':>8}  "
               f"{'avg_compl':>10}  {'total_cells':>12}  "
               f"{'Pareto':>7}  {'score':>7}")
        print(hdr); print("  " + "─" * (len(hdr) - 2))
        for r in rows:
            flag = "  ★" if r["pareto"] else ""
            print(f"  {r['threshold']:>10.2f}  {r['n_genes']:>8d}  "
                  f"{r['avg_completeness']:>10.4f}  {r['total_cells']:>12d}  "
                  f"{'YES' if r['pareto'] else 'no':>7}  "
                  f"{r['score']:>7.4f}{flag}")
        best = max((rows[i] for i in pareto_idx), key=lambda r: r["score"])
        print(f"\n  ✔  Recommended (mode={mode}): threshold={best['threshold']:.2f}  "
              f"→  {best['n_genes']} genes,  "
              f"avg completeness={best['avg_completeness']:.1%}")

        tf = out_dir / f"pareto_sweep_{mode}.tsv"
        with open(tf, "w", newline="") as f:
            w = csv.DictWriter(f, delimiter="\t",
                               fieldnames=["threshold", "n_genes",
                                           "avg_completeness", "total_cells",
                                           "pareto", "score"])
            w.writeheader()
            for r in rows:
                w.writerow({k: v for k, v in r.items() if k != "gene_ids"})
        substep(f"Pareto table: {tf}")

    print()
    print("  Decision guide:")
    print("    ★  Pareto-optimal = no other threshold is simultaneously better"
          " on both metrics.")
    print("    •  Strict phylogenomics  → highest avg_completeness Pareto point")
    print("    •  Broad gene sampling   → highest n_genes Pareto point")
    print("    •  Balanced default      → highest composite score  (★ recommended)")
    print()
    print("  Tip: add --plot-pareto to open an interactive HTML decision tool.")
    print()


# ──────────────────────────────────────────────────────────────────────────────
# INTERACTIVE PARETO HTML  (v7.1 original)
# ──────────────────────────────────────────────────────────────────────────────

def generate_pareto_html(
    sd: SpeciesDict, all_ids: Set[str],
    modes: List[str], seqtype: str, out_dir: Path,
) -> None:
    substep("Computing sweep data for HTML …")
    thresholds = [round(0.50 + 0.05 * i, 3) for i in range(11)]
    all_mode_data: Dict[str, Any] = {}
    for mode in modes:
        rows = evaluate_threshold_sweep(sd, all_ids, mode, seqtype, thresholds)
        pareto_idx = set(find_pareto_front(rows))
        max_g = max(r["n_genes"] for r in rows) or 1
        max_c = max(r["avg_completeness"] for r in rows) or 1.0
        for i, r in enumerate(rows):
            r["pareto"] = i in pareto_idx
            r["score"]  = (r["n_genes"] / max_g) * (r["avg_completeness"] / max_c)
        all_mode_data[mode] = rows

    n_samples   = len(sd)
    total_genes = len(all_ids)
    data_json   = json.dumps(all_mode_data, separators=(",", ":"))
    modes_json  = json.dumps(modes)

    radio_html = ""
    for i, m in enumerate(modes):
        lbl = {"S": "S — Single-copy only",
               "F": "F — Fragmented only",
               "SF": "SF — Single + Fragmented"}.get(m, m)
        chk = "checked" if i == 0 else ""
        radio_html += (
            f'<input type="radio" name="mode" id="mode_{m}" '
            f'value="{m}" {chk} onchange="switchMode(this.value)">'
            f'<label class="mode-btn" for="mode_{m}">{lbl}</label>\n    '
        )

    html = f"""<!DOCTYPE html>
<html lang="en">
<head>
<meta charset="UTF-8">
<meta name="viewport" content="width=device-width,initial-scale=1.0">
<title>BUSCO Pareto Threshold Selector — v7.1</title>
<script src="https://cdn.plot.ly/plotly-2.32.0.min.js"></script>
<style>
:root{{--bg:#f5f7fa;--card:#fff;--border:#dde1e7;--pri:#2563eb;--red:#dc2626;
      --grey:#9ca3af;--txt:#1f2937;--sub:#6b7280;--r:10px;
      --sh:0 2px 12px rgba(0,0,0,.08)}}
*{{box-sizing:border-box;margin:0;padding:0}}
body{{font-family:-apple-system,BlinkMacSystemFont,"Segoe UI",sans-serif;
     background:var(--bg);color:var(--txt);padding:16px}}
h1{{font-size:1.3rem;font-weight:700;margin-bottom:4px}}
.sub{{color:var(--sub);font-size:.88rem;margin-bottom:14px}}
.card{{background:var(--card);border:1px solid var(--border);border-radius:var(--r);
       box-shadow:var(--sh);padding:16px;margin-bottom:14px}}
.g2{{display:grid;grid-template-columns:1fr 1fr;gap:14px}}
@media(max-width:860px){{.g2{{grid-template-columns:1fr}}}}
.ctrl{{display:flex;flex-wrap:wrap;gap:8px;align-items:center;margin-bottom:12px}}
label.mode-btn{{cursor:pointer;padding:6px 15px;border-radius:20px;
  border:2px solid var(--border);font-size:.86rem;font-weight:600;transition:.15s}}
input[name=mode]{{display:none}}
input[name=mode]:checked+label.mode-btn{{background:var(--pri);
  border-color:var(--pri);color:#fff}}
.sw{{display:flex;align-items:center;gap:10px;margin-bottom:12px}}
#tslider{{flex:1;accent-color:var(--pri);height:6px;cursor:pointer}}
#tlabel{{font-size:1.1rem;font-weight:700;color:var(--pri);min-width:46px;text-align:center}}
.sg{{display:grid;grid-template-columns:repeat(auto-fill,minmax(145px,1fr));gap:8px}}
.si{{background:#f0f4ff;border-radius:8px;padding:10px 12px;text-align:center}}
.sv{{font-size:1.45rem;font-weight:800;color:var(--pri)}}
.sl{{font-size:.76rem;color:var(--sub);margin-top:2px}}
#gene-list{{font-family:"JetBrains Mono","Courier New",monospace;font-size:.79rem;
  width:100%;height:200px;resize:vertical;border:1px solid var(--border);
  border-radius:6px;padding:8px;background:#fafafa;color:var(--txt)}}
.dlbtn{{display:inline-block;margin-top:8px;padding:8px 18px;
  background:var(--pri);color:#fff;border:none;border-radius:6px;
  cursor:pointer;font-size:.88rem;font-weight:600;transition:.15s}}
.dlbtn:hover{{background:#1d4ed8}}
.cpbtn{{display:inline-block;margin-top:8px;margin-left:8px;padding:8px 14px;
  background:#e5e7eb;color:var(--txt);border:none;border-radius:6px;
  cursor:pointer;font-size:.88rem;font-weight:600}}
.cpbtn:hover{{background:#d1d5db}}
.tip{{background:#fffbeb;border-left:4px solid #f59e0b;padding:10px 14px;
  border-radius:0 8px 8px 0;font-size:.84rem;color:#92400e;margin-top:10px}}
.leg{{display:flex;gap:14px;flex-wrap:wrap;margin-top:6px;font-size:.8rem;color:var(--sub)}}
.dot{{display:inline-block;width:12px;height:12px;border-radius:50%;
  margin-right:4px;vertical-align:middle}}
</style>
</head>
<body>
<div class="card">
  <h1>🔬 BUSCO Pareto Threshold Decision Tool</h1>
  <div class="sub">Samples: <strong>{n_samples}</strong> &nbsp;|&nbsp;
    Unique BUSCOs: <strong>{total_genes}</strong> &nbsp;|&nbsp;
    Seqtype: <strong>{seqtype}</strong> &nbsp;|&nbsp;
    BUSCO Phylogenomics Extractor v7.1</div>
  <div class="ctrl">
    <strong>Mode:</strong>
    {radio_html}
  </div>
  <div class="sw">
    <span style="font-size:.88rem;font-weight:600;color:var(--sub)">Threshold θ:</span>
    <input type="range" id="tslider" min="0" max="10" step="1" value="4"
           oninput="onSlider(this.value)">
    <span id="tlabel">0.70</span>
  </div>
  <div class="sg">
    <div class="si"><div class="sv" id="sng">—</div><div class="sl">Genes retained</div></div>
    <div class="si"><div class="sv" id="scv">—</div><div class="sl">Avg sample completeness</div></div>
    <div class="si"><div class="sv" id="stc">—</div><div class="sl">Total data cells</div></div>
    <div class="si"><div class="sv" id="spo">—</div><div class="sl">Pareto optimal?</div></div>
  </div>
  <div class="tip" id="dtip">💡 Click a point on the left chart or drag the slider above.</div>
</div>
<div class="g2">
  <div class="card">
    <strong>🎯 Pareto Front</strong>
    <div class="leg">
      <span><span class="dot" style="background:var(--red)"></span>Pareto-optimal ★</span>
      <span><span class="dot" style="background:var(--grey)"></span>Dominated</span>
      <span><span class="dot" style="background:var(--pri);border:3px solid #1d4ed8"></span>Selected</span>
    </div>
    <div id="pp" style="height:360px"></div>
  </div>
  <div class="card">
    <strong>📊 Threshold Sweep Curves</strong>
    <div style="color:var(--sub);font-size:.8rem;margin-bottom:4px">
      Blue: genes retained (left axis) &nbsp;|&nbsp; Green: avg completeness (right axis)
    </div>
    <div id="sp" style="height:360px"></div>
  </div>
</div>
<div class="card">
  <strong>📋 Selected BUSCO Gene List</strong>
  <span style="color:var(--sub);font-size:.82rem;margin-left:8px">
    Use this list as input for MAFFT / IQ-TREE etc.</span><br><br>
  <textarea id="gene-list" readonly placeholder="Select a threshold above …"></textarea><br>
  <button class="dlbtn" onclick="dlGenes()">⬇ Download gene_list.txt</button>
  <button class="cpbtn" onclick="cpGenes()">📋 Copy to clipboard</button>
  <span id="cmsg" style="display:none;color:green;margin-left:10px;font-size:.83rem">Copied!</span>
</div>
<script>
const ALL = {data_json};
const MODES = {modes_json};
let curMode = MODES[0], curIdx = 4;
function rows() {{ return ALL[curMode]; }}
function row(i) {{ return rows()[i]; }}
function updateStats(i) {{
  const r = row(i);
  document.getElementById('sng').textContent  = r.n_genes.toLocaleString();
  document.getElementById('scv').textContent  = (r.avg_completeness*100).toFixed(1)+'%';
  document.getElementById('stc').textContent  = r.total_cells.toLocaleString();
  document.getElementById('spo').textContent  = r.pareto ? '✅ YES' : '❌ No';
  document.getElementById('tlabel').textContent = r.threshold.toFixed(2);
  let tip='';
  if(r.pareto){{
    if(r.avg_completeness>=0.90) tip='✅ <strong>Strict phylogenomics</strong>: high completeness → minimal missing data per sample.';
    else if(r.n_genes>=rows()[0].n_genes*0.8) tip='⚖️ <strong>Balanced</strong>: good gene count with solid completeness.';
    else tip='📊 <strong>Gene-rich</strong>: maximum loci; lower completeness per sample.';
  }}else{{
    const b=rows().find(r2=>r2.n_genes>=r.n_genes&&r2.avg_completeness>=r.avg_completeness&&r2!==r);
    tip= b ? `❌ <strong>Dominated</strong>: θ=${{b.threshold.toFixed(2)}} gives ≥ genes AND higher completeness.`
           : '⚠️ Not Pareto-optimal. A neighbouring threshold improves one metric.';
  }}
  document.getElementById('dtip').innerHTML='💡 '+tip;
  document.getElementById('gene-list').value = r.gene_ids.join('\\n');
}}
function buildPareto(si) {{
  const rs=rows(); const pfr=rs.filter(r=>r.pareto).sort((a,b)=>a.n_genes-b.n_genes);
  const sel=rs[si];
  const t=[
    {{x:pfr.map(r=>r.n_genes),y:pfr.map(r=>r.avg_completeness),mode:'lines',
      line:{{color:'#dc2626',width:2,dash:'dot'}},name:'Pareto front',hoverinfo:'skip'}},
    {{x:rs.map(r=>r.n_genes),y:rs.map(r=>r.avg_completeness),mode:'markers+text',
      marker:{{size:rs.map(r=>r.pareto?18:11),color:rs.map(r=>r.pareto?'#dc2626':'#9ca3af'),
               line:{{width:1.5,color:'#374151'}}}},
      text:rs.map(r=>'θ='+r.threshold.toFixed(2)),textposition:'top center',
      textfont:{{size:10}},customdata:rs.map((_,i)=>i),
      hovertemplate:rs.map(r=>'<b>θ='+r.threshold.toFixed(2)+'</b><br>Genes: '+r.n_genes+'<br>Avg completeness: '+(r.avg_completeness*100).toFixed(1)+'%<br>Pareto: '+(r.pareto?'✅ YES':'❌ No')+'<br>Score: '+r.score.toFixed(4)+'<extra></extra>'),
      name:'Threshold',showlegend:false}},
    {{x:[sel.n_genes],y:[sel.avg_completeness],mode:'markers',
      marker:{{size:24,color:'#2563eb',symbol:'circle-open',line:{{width:3,color:'#1d4ed8'}}}},
      name:'Selected',hoverinfo:'skip',showlegend:false}}
  ];
  Plotly.newPlot('pp',t,{{xaxis:{{title:'Genes retained'}},yaxis:{{title:'Avg completeness',tickformat:'.0%',range:[0,1.05]}},margin:{{l:60,r:10,t:20,b:50}},legend:{{orientation:'h',y:1.08}},hovermode:'closest',plot_bgcolor:'#fafafa'}},{{responsive:true,modeBarButtonsToRemove:['lasso2d','select2d']}});
  document.getElementById('pp').on('plotly_click',function(d){{const pt=d.points[0];if(pt.curveNumber===1) select(pt.pointIndex);}});
}}
function buildSweep(si) {{
  const rs=rows(); const sel=rs[si]; const maxG=Math.max(...rs.map(r=>r.n_genes));
  const t=[
    {{x:rs.map(r=>r.threshold),y:rs.map(r=>r.n_genes),mode:'lines+markers',name:'Genes retained',line:{{color:'#2563eb',width:3}},marker:{{size:7}},yaxis:'y1',hovertemplate:'θ=%{{x:.2f}}  Genes=%{{y}}<extra></extra>'}},
    {{x:rs.map(r=>r.threshold),y:rs.map(r=>r.avg_completeness),mode:'lines+markers',name:'Avg completeness',line:{{color:'#16a34a',width:3,dash:'dash'}},marker:{{size:7}},yaxis:'y2',hovertemplate:'θ=%{{x:.2f}}  Compl=%{{y:.3f}}<extra></extra>'}},
    {{x:[sel.threshold,sel.threshold],y:[0,maxG],mode:'lines',line:{{color:'#dc2626',width:2,dash:'dot'}},name:'Selected θ',yaxis:'y1',hoverinfo:'skip'}}
  ];
  Plotly.newPlot('sp',t,{{xaxis:{{title:'Threshold (θ)',dtick:0.05}},yaxis:{{title:'Genes retained',titlefont:{{color:'#2563eb'}},tickfont:{{color:'#2563eb'}}}},yaxis2:{{title:'Avg completeness',titlefont:{{color:'#16a34a'}},tickfont:{{color:'#16a34a'}},overlaying:'y',side:'right',tickformat:'.0%',range:[0,1.05]}},margin:{{l:55,r:60,t:20,b:50}},legend:{{orientation:'h',y:1.08}},hovermode:'x unified',plot_bgcolor:'#fafafa'}},{{responsive:true,displayModeBar:false}});
}}
function select(i) {{ curIdx=i; document.getElementById('tslider').value=i; updateStats(i); buildPareto(i); buildSweep(i); }}
function onSlider(v) {{ select(parseInt(v)); }}
function switchMode(m) {{ curMode=m; updateStats(curIdx); buildPareto(curIdx); buildSweep(curIdx); }}
function dlGenes() {{
  const r=row(curIdx);
  const blob=new Blob([r.gene_ids.join('\\n')+'\\n'],{{type:'text/plain'}});
  const a=document.createElement('a'); a.href=URL.createObjectURL(blob);
  a.download='busco_gene_list_mode'+curMode+'_t'+r.threshold.toFixed(2)+'.txt';
  a.click(); URL.revokeObjectURL(a.href);
}}
function cpGenes() {{
  navigator.clipboard.writeText(row(curIdx).gene_ids.join('\\n')).then(()=>{{
    const m=document.getElementById('cmsg'); m.style.display='inline';
    setTimeout(()=>m.style.display='none',1800);
  }});
}}
(function(){{ updateStats(curIdx); buildPareto(curIdx); buildSweep(curIdx); }})();
</script>
</body>
</html>"""

    out_path = out_dir / "pareto_interactive.html"
    out_path.write_text(html, encoding="utf-8")
    info(f"Interactive Pareto HTML  : {out_path}")
    substep("Open in any modern browser (Firefox, Chrome, Edge, Safari).")


# ──────────────────────────────────────────────────────────────────────────────
# ██████████████████████████████████████████████████████████████████████████████
#
#   NEW IN v7.2  —  EXPLORE VIEW
#   A stand-alone rich HTML dashboard that works from ANY previous output dir
#   or a specific matrix file.  No BUSCO raw data needed.
#
# ██████████████████████████████████████████████████████████████████████████████
# ──────────────────────────────────────────────────────────────────────────────

def _read_tsv(path: Path) -> Tuple[List[str], List[List[str]]]:
    """Return (header_list, data_rows) from a TSV file."""
    with open(path, "rt", encoding="utf-8", errors="replace") as f:
        reader = csv.reader(f, delimiter="\t")
        rows = [r for r in reader if r]
    if not rows:
        return [], []
    return rows[0], rows[1:]


def load_explore_data(explore_path: Path) -> Dict[str, Any]:
    """
    Load matrices/summaries from a previous run directory OR a specific
    busco_status_matrix.tsv / busco_binary_*.tsv file.

    Returns a unified dict with keys:
      source          str – resolved path used
      samples         List[str]
      genes           List[str]  – all gene IDs
      coverage        Dict[mode, List[float]]  – per-gene coverage in index order
      status_rows     List[dict]  – {gene, sample_status_map}  (may be empty if large)
      sample_stats    List[dict]  – one entry per sample
      pareto_sweeps   Dict[mode, List[dict]]  – pre-computed sweeps (may be empty)
      modes           List[str]
      seqtype         str
      has_heatmap     bool  – whether status_rows is populated
      source_files    List[str]  – files actually read
    """
    result: Dict[str, Any] = {
        "source":       str(explore_path),
        "samples":      [],
        "genes":        [],
        "coverage":     {},
        "status_rows":  [],
        "sample_stats": [],
        "pareto_sweeps":{},
        "modes":        [],
        "seqtype":      "unknown",
        "has_heatmap":  False,
        "source_files": [],
    }

    # ── Determine what files to read ────────────────────────────────────────
    if explore_path.is_file():
        # Single file supplied
        target_files: Dict[str, Optional[Path]] = {"given": explore_path}
        base_dir = explore_path.parent
    else:
        base_dir = explore_path
        target_files = {
            "status_matrix": next(
                (p for p in [base_dir / "busco_status_matrix.tsv"]
                 if p.exists()), None),
            "species_summary": next(
                (p for p in [base_dir / "species_summary.tsv"]
                 if p.exists()), None),
        }
        # Binary matrices for every mode
        for m in ["S", "F", "SF"]:
            mp = base_dir / f"busco_binary_{m}.tsv"
            if mp.exists():
                target_files[f"binary_{m}"] = mp
        # Pareto sweep files
        for m in ["S", "F", "SF"]:
            pp = base_dir / f"pareto_sweep_{m}.tsv"
            if pp.exists():
                target_files[f"pareto_{m}"] = pp

    # ── Load status matrix ───────────────────────────────────────────────────
    status_file = (target_files.get("status_matrix")
                   or target_files.get("given"))
    if status_file and status_file.exists():
        result["source_files"].append(str(status_file))
        substep(f"Reading status matrix : {status_file.name}")
        header, data_rows = _read_tsv(status_file)
        if header:
            # Identify sample columns (exclude BUSCO_ID, cov_S, cov_SF)
            sample_cols: List[int] = []
            cov_s_idx  = cov_sf_idx = -1
            for i, h in enumerate(header):
                if i == 0: continue  # BUSCO_ID
                if h.lower() == "cov_s":  cov_s_idx = i
                elif h.lower() == "cov_sf": cov_sf_idx = i
                else: sample_cols.append(i)

            samples = [header[i] for i in sample_cols]
            result["samples"] = samples

            genes: List[str] = []
            cov_s_vals:  List[float] = []
            cov_sf_vals: List[float] = []
            status_rows: List[Dict[str, Any]] = []

            for row in data_rows:
                if not row: continue
                gene = row[0]
                genes.append(gene)
                # Coverage
                try:
                    cs  = float(row[cov_s_idx])  if cov_s_idx  != -1 and cov_s_idx  < len(row) else 0.0
                    csf = float(row[cov_sf_idx]) if cov_sf_idx != -1 and cov_sf_idx < len(row) else 0.0
                except ValueError:
                    cs = csf = 0.0
                cov_s_vals.append(cs)
                cov_sf_vals.append(csf)

                # Status per sample  (only keep if ≤ 3000 genes for heatmap)
                if len(data_rows) <= 3000:
                    smap = {samples[j]: row[sample_cols[j]]
                            if sample_cols[j] < len(row) else "M"
                            for j in range(len(sample_cols))}
                    status_rows.append({"gene": gene, "status": smap})

            result["genes"] = genes
            result["coverage"]["S"]  = cov_s_vals
            result["coverage"]["SF"] = cov_sf_vals
            if "S" not in result["modes"]:  result["modes"].append("S")
            if "SF" not in result["modes"]: result["modes"].append("SF")

            if status_rows:
                result["status_rows"] = status_rows
                result["has_heatmap"] = True

    # ── Load binary matrices (supplement / override coverage) ────────────────
    for m in ["S", "F", "SF"]:
        bkey = f"binary_{m}"
        bf   = target_files.get(bkey)
        if bf and bf.exists():
            result["source_files"].append(str(bf))
            substep(f"Reading binary matrix : {bf.name}")
            header, data_rows = _read_tsv(bf)
            if header:
                cov_col = len(header) - 1  # last column = coverage
                # Detect sample columns
                s_cols = [i for i, h in enumerate(header)
                          if i > 0 and h.lower() != "coverage"]
                if not result["samples"] and s_cols:
                    result["samples"] = [header[i] for i in s_cols]
                gene_list: List[str] = []
                cov_vals:  List[float] = []
                for row in data_rows:
                    if not row: continue
                    gene_list.append(row[0])
                    try:
                        cov_vals.append(float(row[cov_col]))
                    except (ValueError, IndexError):
                        cov_vals.append(0.0)
                if not result["genes"]:
                    result["genes"] = gene_list
                result["coverage"][m] = cov_vals
                if m not in result["modes"]:
                    result["modes"].append(m)

    # ── Load species summary ─────────────────────────────────────────────────
    ss_file = target_files.get("species_summary")
    if ss_file and ss_file.exists():
        result["source_files"].append(str(ss_file))
        substep(f"Reading species summary: {ss_file.name}")
        header, data_rows = _read_tsv(ss_file)
        hdr_lower = [h.lower() for h in header]
        def _gi(name: str) -> int:
            return hdr_lower.index(name) if name in hdr_lower else -1
        i_name   = _gi("species")
        i_total  = _gi("n_total")
        i_single = _gi("n_single")
        i_dup    = _gi("n_duplicated")
        i_frag   = _gi("n_fragmented")
        i_miss   = _gi("n_missing")
        i_lin    = _gi("lineage")
        stats = []
        for row in data_rows:
            if not row: continue
            def gv(idx: int, default: Any = "") -> Any:
                return row[idx] if idx != -1 and idx < len(row) else default
            try:
                n_total  = int(gv(i_total,  0) or 0)
                n_single = int(gv(i_single, 0) or 0)
                n_frag   = int(gv(i_frag,   0) or 0)
                n_miss   = int(gv(i_miss,   0) or 0)
                completeness = (n_single / n_total) if n_total else 0.0
                completeness_sf = ((n_single + n_frag) / n_total) if n_total else 0.0
            except (ValueError, ZeroDivisionError):
                completeness = completeness_sf = 0.0
                n_total = n_single = n_frag = n_miss = 0
            stats.append({
                "name":    gv(i_name),
                "lineage": gv(i_lin),
                "n_total": n_total,
                "n_single": n_single,
                "n_frag":   n_frag,
                "n_miss":   n_miss,
                "completeness_S":  round(completeness, 4),
                "completeness_SF": round(completeness_sf, 4),
            })
        result["sample_stats"] = stats
        if not result["samples"] and stats:
            result["samples"] = [s["name"] for s in stats]

    # ── Load pre-computed pareto sweeps ──────────────────────────────────────
    for m in ["S", "F", "SF"]:
        pk = f"pareto_{m}"
        pf = target_files.get(pk)
        if pf and pf.exists():
            result["source_files"].append(str(pf))
            substep(f"Reading pareto sweep  : {pf.name}")
            header, data_rows = _read_tsv(pf)
            if header:
                hdr_lower = [h.lower() for h in header]
                sweep_rows = []
                for row in data_rows:
                    if not row: continue
                    try:
                        d = {h: row[i] for i, h in enumerate(hdr_lower) if i < len(row)}
                        sweep_rows.append({
                            "threshold":       float(d.get("threshold", 0)),
                            "n_genes":         int(d.get("n_genes", 0)),
                            "avg_completeness":float(d.get("avg_completeness", 0)),
                            "total_cells":     int(d.get("total_cells", 0)),
                            "pareto":          str(d.get("pareto", "false")).lower() == "true",
                            "score":           float(d.get("score", 0)),
                        })
                    except (ValueError, KeyError):
                        continue
                result["pareto_sweeps"][m] = sweep_rows

    # ── Compute pareto sweeps from coverage if not loaded ────────────────────
    for m in result["modes"]:
        if m not in result["pareto_sweeps"] and m in result["coverage"]:
            cov_list = result["coverage"][m]
            genes    = result["genes"]
            n_sp     = len(result["samples"]) or 1
            thresholds = [round(0.50 + 0.05 * i, 3) for i in range(11)]
            sweep_rows = []
            for t in thresholds:
                kept_ids  = [genes[i] for i, c in enumerate(cov_list) if c >= t]
                n_genes   = len(kept_ids)
                # avg completeness = mean per-sample fraction of kept genes present
                # We approximate using the coverage values of kept genes:
                # avg_completeness ≈ mean(coverage of kept genes)
                avg_c = (sum(cov_list[genes.index(g)] for g in kept_ids) / n_genes
                         if n_genes else 0.0)
                sweep_rows.append({
                    "threshold":       t,
                    "n_genes":         n_genes,
                    "avg_completeness": round(avg_c, 5),
                    "total_cells":     n_genes * n_sp,
                    "gene_ids":        kept_ids,
                    "pareto":          False,
                    "score":           0.0,
                })
            # Mark pareto front
            pareto_idx = find_pareto_front(sweep_rows)
            max_g = max(r["n_genes"] for r in sweep_rows) or 1
            max_c = max(r["avg_completeness"] for r in sweep_rows) or 1.0
            for i, r in enumerate(sweep_rows):
                r["pareto"] = i in pareto_idx
                r["score"]  = (r["n_genes"] / max_g) * (r["avg_completeness"] / max_c)
            result["pareto_sweeps"][m] = sweep_rows

    # ── Add gene_ids to pareto sweep rows if missing ─────────────────────────
    for m, sweep in result["pareto_sweeps"].items():
        cov_list = result["coverage"].get(m, [])
        genes    = result["genes"]
        for r in sweep:
            if "gene_ids" not in r:
                t = r["threshold"]
                r["gene_ids"] = [genes[i] for i, c in enumerate(cov_list) if c >= t]

    substep(f"Loaded: {len(result['genes'])} genes, "
            f"{len(result['samples'])} samples, "
            f"modes: {result['modes']}")
    return result


def generate_explore_html(ed: Dict[str, Any], out_path: Path) -> None:
    """
    Generate a comprehensive multi-tab HTML dashboard from explore data.

    Tabs:
      1. Overview    – key stats cards + source-file audit
      2. Coverage    – gene-coverage histogram + coverage CDF
      3. Completeness– per-sample completeness bar chart (ranked, colour-coded)
      4. Pareto      – full interactive Pareto scatter + sweep (all modes)
      5. Advisor     – smart textual recommendation with detailed rationale
    """
    substep("Building comprehensive HTML dashboard …")

    samples       = ed["samples"]
    genes         = ed["genes"]
    modes         = ed["modes"] or ["S"]
    coverage      = ed["coverage"]        # mode -> List[float]
    sample_stats  = ed["sample_stats"]    # List[dict]
    pareto_sweeps = ed["pareto_sweeps"]   # mode -> List[dict]
    source_files  = ed["source_files"]
    n_samples     = len(samples)
    n_genes       = len(genes)

    # ── Build summary stats for advisor ──────────────────────────────────────
    advisor_rows: List[Dict[str, Any]] = []
    for m in modes:
        sweep = pareto_sweeps.get(m, [])
        if not sweep:
            continue
        pareto_pts  = [r for r in sweep if r["pareto"]]
        if not pareto_pts:
            continue
        # Best composite score
        best_score  = max(pareto_pts, key=lambda r: r["score"])
        # Strictest (highest completeness Pareto)
        best_strict = max(pareto_pts, key=lambda r: r["avg_completeness"])
        # Most genes Pareto
        best_genes  = max(pareto_pts, key=lambda r: r["n_genes"])

        def classify(r: Dict) -> str:
            ac = r["avg_completeness"]
            ng = r["n_genes"]
            max_ng = sweep[0]["n_genes"] if sweep else 1
            if ac >= 0.90: return "strict"
            if ng >= 0.80 * max_ng: return "gene-rich"
            return "balanced"

        advisor_rows.append({
            "mode": m,
            "best_score":  best_score,
            "best_strict": best_strict,
            "best_genes":  best_genes,
            "n_pareto":    len(pareto_pts),
            "classify_score":  classify(best_score),
            "classify_strict": classify(best_strict),
        })

    # ── Serialize data for JS ─────────────────────────────────────────────────
    # Coverage histograms (bin into 20 buckets 0→1)
    cov_hist: Dict[str, List[int]] = {}
    for m in modes:
        cv = coverage.get(m, [])
        bins = [0] * 20
        for v in cv:
            idx = min(int(v * 20), 19)
            bins[idx] += 1
        cov_hist[m] = bins

    # Coverage CDF (sorted)
    cov_sorted: Dict[str, List[float]] = {}
    for m in modes:
        cv = sorted(coverage.get(m, []), reverse=True)
        cov_sorted[m] = cv

    # Sample completeness (from sample_stats or computed)
    samp_comp_S:  List[Tuple[str, float]] = []
    samp_comp_SF: List[Tuple[str, float]] = []
    if sample_stats:
        for s in sample_stats:
            samp_comp_S.append( (s["name"], s["completeness_S"]))
            samp_comp_SF.append((s["name"], s["completeness_SF"]))
    # Sort descending
    samp_comp_S.sort(key=lambda x: -x[1])
    samp_comp_SF.sort(key=lambda x: -x[1])

    # Heatmap sub-sample (max 250 genes × all samples)
    heatmap_data: Dict[str, Any] = {}
    if ed.get("has_heatmap") and ed.get("status_rows"):
        sr = ed["status_rows"]
        # Sample evenly-spaced + prioritise interesting genes (mid-coverage)
        step_  = max(1, len(sr) // 250)
        sub    = sr[::step_][:250]
        sg     = [r["gene"] for r in sub]
        ss_    = samples[:80]  # limit samples to 80 columns for readability
        STATUS_NUM = {"S": 0, "D": 1, "F": 2, "M": 3}
        zmat   = [[STATUS_NUM.get(r["status"].get(sn, "M"), 3) for sn in ss_]
                  for r in sub]
        heatmap_data = {
            "genes":    sg,
            "samples":  ss_,
            "z":        zmat,
        }

    # Pack pareto sweeps (strip gene_ids for JS; keep separate)
    ps_for_js: Dict[str, List[Dict]] = {}
    gene_lists_js: Dict[str, Dict[str, List[str]]] = {}
    for m, sweep in pareto_sweeps.items():
        slim = []
        gl   = {}
        for r in sweep:
            slim.append({k: v for k, v in r.items() if k != "gene_ids"})
            gl[str(round(r["threshold"], 2))] = r.get("gene_ids", [])
        ps_for_js[m] = slim
        gene_lists_js[m] = gl

    data_json     = json.dumps({
        "n_samples":   n_samples,
        "n_genes":     n_genes,
        "modes":       modes,
        "covHist":     cov_hist,
        "covSorted":   cov_sorted,
        "sampleCompS": samp_comp_S,
        "sampleCompSF":samp_comp_SF,
        "heatmap":     heatmap_data,
        "paretoSweeps":ps_for_js,
        "geneLists":   gene_lists_js,
        "sourceFiles": source_files,
        "advisorRows": advisor_rows,
    }, separators=(",", ":"))

    modes_json    = json.dumps(modes)
    source        = ed["source"]
    n_sources     = len(source_files)

    # ── Mode colour map ───────────────────────────────────────────────────────
    mode_colours = {"S": "#2563eb", "F": "#d97706", "SF": "#7c3aed"}

    html = f"""<!DOCTYPE html>
<html lang="en">
<head>
<meta charset="UTF-8">
<meta name="viewport" content="width=device-width,initial-scale=1.0">
<title>BUSCO Explore View — v7.2</title>
<script src="https://cdn.plot.ly/plotly-2.32.0.min.js"></script>
<style>
:root{{
  --bg:#f0f4f8;--card:#fff;--border:#e2e8f0;--pri:#2563eb;--pri2:#1d4ed8;
  --ok:#16a34a;--warn:#d97706;--err:#dc2626;--pur:#7c3aed;
  --txt:#1e293b;--sub:#64748b;--muted:#94a3b8;--r:12px;
  --sh:0 1px 8px rgba(0,0,0,.07);--shh:0 4px 20px rgba(0,0,0,.12)
}}
*{{box-sizing:border-box;margin:0;padding:0}}
body{{font-family:-apple-system,BlinkMacSystemFont,"Segoe UI",Roboto,sans-serif;
      background:var(--bg);color:var(--txt);padding:0;min-height:100vh}}

/* ── NAV ── */
.topbar{{background:linear-gradient(135deg,#1e3a8a 0%,#2563eb 60%,#7c3aed 100%);
  padding:18px 24px;color:#fff}}
.topbar h1{{font-size:1.4rem;font-weight:800;letter-spacing:-.5px}}
.topbar .meta{{font-size:.82rem;opacity:.82;margin-top:4px}}
.tabs{{display:flex;gap:0;border-bottom:2px solid var(--border);
  background:#fff;padding:0 20px;overflow-x:auto;position:sticky;top:0;z-index:100;
  box-shadow:0 2px 8px rgba(0,0,0,.06)}}
.tab{{padding:13px 20px;cursor:pointer;font-size:.88rem;font-weight:600;
  color:var(--sub);border-bottom:3px solid transparent;margin-bottom:-2px;
  white-space:nowrap;transition:.15s}}
.tab:hover{{color:var(--pri)}}
.tab.active{{color:var(--pri);border-bottom-color:var(--pri)}}
.panel{{display:none;padding:20px 24px;max-width:1400px;margin:0 auto}}
.panel.active{{display:block}}

/* ── CARDS ── */
.card{{background:var(--card);border:1px solid var(--border);border-radius:var(--r);
  box-shadow:var(--sh);padding:18px 20px;margin-bottom:16px}}
.card:hover{{box-shadow:var(--shh)}}
.card-title{{font-size:.92rem;font-weight:700;color:var(--sub);
  text-transform:uppercase;letter-spacing:.5px;margin-bottom:12px}}

/* ── STAT GRID ── */
.sg4{{display:grid;grid-template-columns:repeat(auto-fill,minmax(160px,1fr));gap:12px}}
.si{{background:linear-gradient(135deg,#f0f4ff,#e8f0fe);
  border:1px solid #c7d7fe;border-radius:10px;padding:14px 16px;text-align:center}}
.sv{{font-size:1.7rem;font-weight:800;color:var(--pri);line-height:1}}
.sl{{font-size:.74rem;color:var(--sub);margin-top:5px;font-weight:500}}

/* ── MODE BADGE ── */
.badge{{display:inline-block;padding:3px 10px;border-radius:20px;font-size:.78rem;
  font-weight:700;margin:2px}}
.badge-S{{background:#dbeafe;color:#1d4ed8}}
.badge-F{{background:#fef3c7;color:#92400e}}
.badge-SF{{background:#ede9fe;color:#5b21b6}}

/* ── ADVISOR ── */
.rec-card{{border-left:5px solid;border-radius:0 var(--r) var(--r) 0;
  padding:14px 18px;margin-bottom:12px}}
.rec-strict{{border-color:var(--ok);background:#f0fdf4}}
.rec-balanced{{border-color:var(--pri);background:#eff6ff}}
.rec-gene-rich{{border-color:var(--pur);background:#f5f3ff}}
.rec-title{{font-weight:700;font-size:.95rem;margin-bottom:6px}}
.rec-body{{font-size:.84rem;color:#374151;line-height:1.7}}
.rec-tag{{display:inline-block;padding:2px 8px;border-radius:4px;
  font-size:.75rem;font-weight:700;margin-right:6px}}
.tag-strict{{background:#bbf7d0;color:#166534}}
.tag-balanced{{background:#bfdbfe;color:#1e40af}}
.tag-gene-rich{{background:#ddd6fe;color:#4c1d95}}

/* ── SOURCE FILE LIST ── */
.sf-list{{list-style:none;padding:0}}
.sf-list li{{padding:4px 0;font-size:.82rem;color:var(--sub);
  border-bottom:1px solid #f1f5f9}}
.sf-list li:before{{content:"📄 ";}}

/* ── SWEEP CTRL ── */
.ctrl-row{{display:flex;flex-wrap:wrap;gap:8px;align-items:center;margin-bottom:12px}}
label.mode-btn{{cursor:pointer;padding:7px 16px;border-radius:20px;
  border:2px solid var(--border);font-size:.86rem;font-weight:600;transition:.15s}}
input[name=pmode]{{display:none}}
input[name=pmode]:checked+label.mode-btn{{background:var(--pri);
  border-color:var(--pri);color:#fff}}
.sw{{display:flex;align-items:center;gap:10px;margin-bottom:12px}}
.sw input[type=range]{{flex:1;accent-color:var(--pri)}}
.tval{{font-size:1.1rem;font-weight:800;color:var(--pri);min-width:48px;text-align:center}}
.tip-box{{background:#fffbeb;border-left:4px solid #f59e0b;padding:10px 14px;
  border-radius:0 8px 8px 0;font-size:.84rem;color:#78350f;margin-top:10px}}
.dlbtn{{padding:8px 18px;background:var(--pri);color:#fff;border:none;
  border-radius:8px;cursor:pointer;font-weight:600;font-size:.88rem;transition:.15s}}
.dlbtn:hover{{background:var(--pri2)}}
.cpbtn{{padding:8px 14px;background:#e2e8f0;color:var(--txt);border:none;
  border-radius:8px;cursor:pointer;font-weight:600;font-size:.88rem;margin-left:6px}}
#ev-gene-list{{font-family:"JetBrains Mono","Courier New",monospace;font-size:.78rem;
  width:100%;height:180px;resize:vertical;border:1px solid var(--border);
  border-radius:8px;padding:8px;background:#fafafa}}
.g2{{display:grid;grid-template-columns:1fr 1fr;gap:14px}}
@media(max-width:900px){{.g2{{grid-template-columns:1fr}}}}
.plot-box{{height:360px}}
.leg-row{{display:flex;gap:14px;flex-wrap:wrap;margin-top:4px;
  font-size:.78rem;color:var(--sub)}}
.dot{{width:11px;height:11px;border-radius:50%;display:inline-block;
  margin-right:3px;vertical-align:middle}}
</style>
</head>
<body>

<div class="topbar">
  <h1>🔬 BUSCO Explore View &nbsp;<span style="opacity:.7;font-size:.85rem">v7.2</span></h1>
  <div class="meta">
    Source: <code>{source}</code> &nbsp;|&nbsp;
    {n_sources} matrix file(s) loaded &nbsp;|&nbsp;
    {n_samples} samples &nbsp;·&nbsp; {n_genes} genes &nbsp;·&nbsp;
    Modes: {''.join(f'<span class="badge badge-{m}">{m}</span>' for m in modes)}
  </div>
</div>

<div class="tabs">
  <div class="tab active" onclick="showTab('ov')">📊 Overview</div>
  <div class="tab" onclick="showTab('cov')">📈 Gene Coverage</div>
  <div class="tab" onclick="showTab('comp')">🧬 Sample Completeness</div>
  <div class="tab" onclick="showTab('par')">🎯 Pareto Advisor</div>
  <div class="tab" onclick="showTab('adv')">💡 Smart Recommendation</div>
</div>

<!-- ══════════════════════  TAB 1: OVERVIEW  ═══════════════════════════════ -->
<div id="tab-ov" class="panel active">
  <div class="card">
    <div class="card-title">Dataset Summary</div>
    <div class="sg4">
      <div class="si"><div class="sv" id="ov-ns">{n_samples}</div><div class="sl">Samples</div></div>
      <div class="si"><div class="sv" id="ov-ng">{n_genes}</div><div class="sl">Unique BUSCOs</div></div>
      <div class="si"><div class="sv" id="ov-nm">{len(modes)}</div><div class="sl">Modes available</div></div>
      <div class="si"><div class="sv" id="ov-nf">{n_sources}</div><div class="sl">Matrix files read</div></div>
      <div class="si"><div class="sv" id="ov-mcS">—</div><div class="sl">Mean coverage (S)</div></div>
      <div class="si"><div class="sv" id="ov-mcSF">—</div><div class="sl">Mean coverage (SF)</div></div>
      <div class="si"><div class="sv" id="ov-medS">—</div><div class="sl">Median completeness</div></div>
      <div class="si"><div class="sv" id="ov-ng90">—</div><div class="sl">Genes ≥ 90% cov (S)</div></div>
    </div>
  </div>
  <div class="card">
    <div class="card-title">Coverage Distribution Quick View</div>
    <div id="ov-chart" style="height:260px"></div>
  </div>
  <div class="card">
    <div class="card-title">Source Files</div>
    <ul class="sf-list" id="sf-list"></ul>
  </div>
</div>

<!-- ══════════════════════  TAB 2: GENE COVERAGE  ════════════════════════ -->
<div id="tab-cov" class="panel">
  <div class="card">
    <div class="card-title">Gene Coverage Histogram</div>
    <div style="color:var(--sub);font-size:.82rem;margin-bottom:8px">
      Each bar = fraction of samples containing that gene.
      Higher coverage → gene more universally present → better matrix completeness.
    </div>
    <div id="cov-hist" class="plot-box"></div>
  </div>
  <div class="card">
    <div class="card-title">Cumulative Distribution  (sorted coverage)</div>
    <div style="color:var(--sub);font-size:.82rem;margin-bottom:8px">
      At threshold θ (x-axis), y = number of genes with coverage ≥ θ.
      Vertical dashed lines = common thresholds (0.5 · 0.7 · 0.9 · 1.0).
    </div>
    <div id="cov-cdf" class="plot-box"></div>
  </div>
</div>

<!-- ══════════════════════  TAB 3: SAMPLE COMPLETENESS  ════════════════════ -->
<div id="tab-comp" class="panel">
  <div class="card">
    <div class="card-title">Per-Sample Completeness (ranked descending)</div>
    <div style="color:var(--sub);font-size:.82rem;margin-bottom:8px">
      Completeness = fraction of all unique BUSCOs present in each sample.
      Colour: green ≥ 90% · yellow 70–89% · orange 50–69% · red &lt; 50%.
    </div>
    <div id="comp-bar" style="height:500px"></div>
  </div>
  <div class="card">
    <div class="card-title">Completeness Statistics by Mode</div>
    <div id="comp-stats-html"></div>
  </div>
</div>

<!-- ══════════════════════  TAB 4: PARETO ADVISOR  ════════════════════════ -->
<div id="tab-par" class="panel">
  <div class="card">
    <div class="card-title">Interactive Pareto Threshold Selector</div>
    <div class="ctrl-row">
      <strong style="font-size:.88rem">Mode:</strong>
      <div id="par-mode-btns"></div>
    </div>
    <div class="sw">
      <span style="font-size:.84rem;font-weight:600;color:var(--sub)">θ:</span>
      <input type="range" id="par-slider" min="0" max="10" step="1" value="4"
             oninput="parSlider(this.value)">
      <span class="tval" id="par-tval">—</span>
    </div>
    <div class="sg4" style="margin-bottom:10px">
      <div class="si"><div class="sv" id="par-ng">—</div><div class="sl">Genes retained</div></div>
      <div class="si"><div class="sv" id="par-ac">—</div><div class="sl">Avg completeness</div></div>
      <div class="si"><div class="sv" id="par-tc">—</div><div class="sl">Data cells</div></div>
      <div class="si"><div class="sv" id="par-po">—</div><div class="sl">Pareto optimal?</div></div>
    </div>
    <div class="tip-box" id="par-tip">💡 Drag the slider or click a scatter point.</div>
  </div>
  <div class="g2">
    <div class="card">
      <div class="card-title">🎯 Pareto Front Scatter</div>
      <div class="leg-row">
        <span><span class="dot" style="background:#dc2626"></span>Pareto-optimal</span>
        <span><span class="dot" style="background:#9ca3af"></span>Dominated</span>
        <span><span class="dot" style="background:#2563eb;border:2px solid #1d4ed8"></span>Selected</span>
      </div>
      <div id="par-scatter" class="plot-box"></div>
    </div>
    <div class="card">
      <div class="card-title">📊 Threshold Sweep Curves</div>
      <div style="color:var(--sub);font-size:.78rem;margin-bottom:4px">
        Blue: genes retained · Green dashed: avg completeness
      </div>
      <div id="par-sweep" class="plot-box"></div>
    </div>
  </div>
  <div class="card">
    <div class="card-title">📋 Selected BUSCO Gene List</div>
    <textarea id="ev-gene-list" readonly placeholder="Select threshold above …"></textarea><br>
    <button class="dlbtn" onclick="evDl()">⬇ Download gene_list.txt</button>
    <button class="cpbtn" onclick="evCp()">📋 Copy</button>
    <span id="ev-cmsg" style="display:none;color:var(--ok);font-size:.82rem;margin-left:8px">Copied!</span>
  </div>
</div>

<!-- ══════════════════════  TAB 5: SMART RECOMMENDATION  ══════════════════ -->
<div id="tab-adv" class="panel">
  <div class="card">
    <div class="card-title">How to Choose a Threshold — Decision Framework</div>
    <div style="font-size:.88rem;color:var(--txt);line-height:1.9">
      <p>The optimal threshold balances two competing goals:</p>
      <ul style="margin:8px 0 8px 20px">
        <li><strong>Genes retained (↑)</strong> — more loci → stronger phylogenetic signal,
            better resolution of short branches.</li>
        <li><strong>Average sample completeness (↑)</strong> — fewer gaps in the alignment
            matrix → less missing-data bias in tree reconstruction.</li>
      </ul>
      <p>A threshold is <strong>Pareto-optimal</strong> if no other threshold is
         simultaneously better on <em>both</em> objectives.  The set of Pareto-optimal
         thresholds forms the <em>efficient frontier</em> — your decision should come
         from within this set.</p>
    </div>
  </div>
  <div id="adv-cards"></div>
  <div class="card">
    <div class="card-title">All Thresholds at a Glance</div>
    <div id="adv-table"></div>
  </div>
</div>

<script>
// ═══════════════════════════════════════════════════════════════════
//  DATA
// ═══════════════════════════════════════════════════════════════════
const D = {data_json};
const MODES = {modes_json};
let parMode = MODES[0], parIdx = 4;

// ── Tab switching ───────────────────────────────────────────────────
function showTab(id) {{
  document.querySelectorAll('.panel').forEach(p=>p.classList.remove('active'));
  document.querySelectorAll('.tab').forEach(t=>t.classList.remove('active'));
  document.getElementById('tab-'+id).classList.add('active');
  const tabs = ['ov','cov','comp','par','adv'];
  document.querySelectorAll('.tab')[tabs.indexOf(id)].classList.add('active');
  // Lazy-init charts on first visit
  if(id==='cov'  && !window._covInit)  {{ initCov();  window._covInit=true; }}
  if(id==='comp' && !window._compInit) {{ initComp(); window._compInit=true; }}
  if(id==='par'  && !window._parInit)  {{ initPar();  window._parInit=true; }}
  if(id==='adv'  && !window._advInit)  {{ initAdv();  window._advInit=true; }}
}}

// ═══════════════════════════════════════════════════════════════════
//  TAB 1 – OVERVIEW
// ═══════════════════════════════════════════════════════════════════
(function initOv() {{
  const covS  = D.covSorted['S']  || [];
  const covSF = D.covSorted['SF'] || [];
  const mean_ = arr => arr.length ? arr.reduce((a,b)=>a+b,0)/arr.length : 0;
  const med_  = arr => {{ if(!arr.length) return 0; const s=[...arr].sort((a,b)=>a-b); const m=s.length>>1; return s.length%2?s[m]:(s[m-1]+s[m])/2; }};

  document.getElementById('ov-mcS').textContent  = (mean_(covS)*100).toFixed(1)+'%';
  document.getElementById('ov-mcSF').textContent = (mean_(covSF)*100).toFixed(1)+'%';

  const sc = D.sampleCompS || [];
  const comps = sc.map(x=>x[1]);
  document.getElementById('ov-medS').textContent  = comps.length ? (med_(comps)*100).toFixed(1)+'%' : '—';
  document.getElementById('ov-ng90').textContent  = covS.filter(v=>v>=0.9).length.toLocaleString();

  // Source files
  const ul = document.getElementById('sf-list');
  (D.sourceFiles||[]).forEach(f=>{{
    const li=document.createElement('li'); li.textContent=f; ul.appendChild(li);
  }});

  // Quick coverage bar chart
  const COLOURS = {{'S':'#2563eb','F':'#d97706','SF':'#7c3aed'}};
  const xbins = Array.from({{length:20}},(_,i)=>(i*5+2.5).toFixed(0)+'%');
  const traces = MODES.map(m=>{{
    const h = D.covHist[m]||[];
    return {{x:xbins,y:h,type:'bar',name:'Mode '+m,
             marker:{{color:COLOURS[m]||'#888',opacity:.7}}}};
  }});
  Plotly.newPlot('ov-chart',traces,{{
    barmode:'group',
    xaxis:{{title:'Gene coverage bin',tickangle:-45}},
    yaxis:{{title:'Number of genes'}},
    margin:{{l:55,r:10,t:10,b:80}},
    legend:{{orientation:'h',y:1.05}},
    plot_bgcolor:'#fafafa'
  }},{{responsive:true,displayModeBar:false}});
}})();

// ═══════════════════════════════════════════════════════════════════
//  TAB 2 – GENE COVERAGE
// ═══════════════════════════════════════════════════════════════════
function initCov() {{
  const COLOURS = {{'S':'#2563eb','F':'#d97706','SF':'#7c3aed'}};
  const xbins = Array.from({{length:20}},(_,i)=>((i*5+2.5)/100).toFixed(3)*1);

  // Histogram
  const tHist = MODES.map(m=>{{
    const h=D.covHist[m]||[];
    return {{x:xbins.map(v=>(v*100).toFixed(0)+'%'),y:h,type:'bar',name:'Mode '+m,
             marker:{{color:COLOURS[m]||'#888',opacity:.78}}}};
  }});
  Plotly.newPlot('cov-hist',tHist,{{
    barmode:'overlay',
    xaxis:{{title:'Gene coverage (fraction of samples)'}},
    yaxis:{{title:'Number of genes'}},
    margin:{{l:55,r:10,t:10,b:60}},
    legend:{{orientation:'h',y:1.05}},
    plot_bgcolor:'#fafafa'
  }},{{responsive:true,displayModeBar:false}});

  // CDF – genes with coverage ≥ x
  const tCDF = MODES.map(m=>{{
    const sorted = [...(D.covSorted[m]||[])].sort((a,b)=>a-b);
    const n=sorted.length;
    const xs=[], ys=[];
    for(let i=0;i<=100;i++){{
      const t=i/100;
      xs.push(t);
      ys.push(sorted.filter(v=>v>=t).length);
    }}
    return {{x:xs,y:ys,mode:'lines',name:'Mode '+m,
             line:{{color:COLOURS[m]||'#888',width:3}}}};
  }});
  // threshold lines
  [0.5,0.7,0.9,1.0].forEach(t=>{{
    tCDF.push({{x:[t,t],y:[0,D.n_genes],mode:'lines',
      line:{{color:'#dc2626',width:1.5,dash:'dot'}},
      name:'θ='+t.toFixed(1),showlegend:false,
      hovertemplate:'θ='+t.toFixed(2)+'<extra></extra>'}});
  }});
  Plotly.newPlot('cov-cdf',tCDF,{{
    xaxis:{{title:'Threshold (θ)',range:[0,1],tickformat:'.0%'}},
    yaxis:{{title:'Genes retained (≥ θ)'}},
    margin:{{l:60,r:10,t:10,b:50}},
    legend:{{orientation:'h',y:1.05}},
    hovermode:'x',plot_bgcolor:'#fafafa'
  }},{{responsive:true,displayModeBar:false}});
}}

// ═══════════════════════════════════════════════════════════════════
//  TAB 3 – SAMPLE COMPLETENESS
// ═══════════════════════════════════════════════════════════════════
function initComp() {{
  const sc = D.sampleCompS||[];
  if(!sc.length) {{ document.getElementById('comp-bar').textContent='No sample stats available.'; return; }}

  const names = sc.map(x=>x[0]);
  const vals  = sc.map(x=>x[1]);
  const colours = vals.map(v=> v>=0.9?'#16a34a': v>=0.7?'#ca8a04': v>=0.5?'#ea580c':'#dc2626');

  Plotly.newPlot('comp-bar',[{{
    x:vals,y:names,type:'bar',orientation:'h',
    marker:{{color:colours}},
    hovertemplate:'%{{y}}: %{{x:.1%}}<extra></extra>'
  }}],{{
    xaxis:{{title:'Completeness (S mode)',range:[0,1],tickformat:'.0%'}},
    yaxis:{{autorange:'reversed',tickfont:{{size:10}}}},
    margin:{{l:150,r:30,t:20,b:50}},
    plot_bgcolor:'#fafafa'
  }},{{responsive:true,displayModeBar:false}});

  // Stats table
  const allComp = D.sampleCompS||[];
  const bins = [
    ['≥ 90%', allComp.filter(x=>x[1]>=0.9).length, '#dcfce7','#166534'],
    ['70–89%',allComp.filter(x=>x[1]>=0.7&&x[1]<0.9).length,'#fef9c3','#713f12'],
    ['50–69%',allComp.filter(x=>x[1]>=0.5&&x[1]<0.7).length,'#ffedd5','#7c2d12'],
    ['< 50%', allComp.filter(x=>x[1]<0.5).length,'#fee2e2','#7f1d1d'],
  ];
  document.getElementById('comp-stats-html').innerHTML =
    '<table style="width:100%;border-collapse:collapse;font-size:.88rem">' +
    '<tr><th style="text-align:left;padding:6px;border-bottom:2px solid #e2e8f0">Completeness (S)</th>' +
    '<th style="padding:6px;border-bottom:2px solid #e2e8f0;text-align:right">Samples (n)</th>' +
    '<th style="padding:6px;border-bottom:2px solid #e2e8f0;text-align:right">Fraction</th></tr>' +
    bins.map(([lbl,n,bg,col])=>
      `<tr style="background:${{bg}}"><td style="padding:7px 6px;color:${{col}};font-weight:600">${{lbl}}</td>` +
      `<td style="padding:7px 6px;text-align:right;font-weight:700">${{n}}</td>` +
      `<td style="padding:7px 6px;text-align:right">${{allComp.length?((n/allComp.length)*100).toFixed(1)+'%':'—'}}</td></tr>`
    ).join('') + '</table>';
}}

// ═══════════════════════════════════════════════════════════════════
//  TAB 4 – PARETO
// ═══════════════════════════════════════════════════════════════════
function initPar() {{
  // Mode buttons
  const mb = document.getElementById('par-mode-btns');
  MODES.forEach((m,i)=>{{
    const labels = {{'S':'S — Single-copy','F':'F — Fragmented','SF':'SF — Single+Frag'}};
    mb.innerHTML += `<input type="radio" name="pmode" id="pmode_${{m}}" value="${{m}}" ${{i===0?'checked':''}} onchange="parSwitchMode('${{m}}')">`
                  + `<label class="mode-btn" for="pmode_${{m}}">${{labels[m]||m}}</label>`;
  }});
  parUpdateStats(parIdx);
  parBuildScatter(parIdx);
  parBuildSweep(parIdx);
}}

function parRows() {{ return D.paretoSweeps[parMode]||[]; }}
function parRow(i) {{ return parRows()[i]||{{}}; }}

function parUpdateStats(i) {{
  const r=parRow(i);
  document.getElementById('par-ng').textContent  = (r.n_genes||0).toLocaleString();
  document.getElementById('par-ac').textContent  = r.avg_completeness!=null?(r.avg_completeness*100).toFixed(1)+'%':'—';
  document.getElementById('par-tc').textContent  = (r.total_cells||0).toLocaleString();
  document.getElementById('par-po').textContent  = r.pareto?'✅ YES':'❌ No';
  document.getElementById('par-tval').textContent= r.threshold!=null?r.threshold.toFixed(2):'—';
  // Tip
  let tip='';
  if(r.pareto){{
    if(r.avg_completeness>=0.90) tip='✅ <strong>Strict</strong>: high completeness, minimal missing data. Best for coalescent trees and sensitive analyses.';
    else if(r.avg_completeness>=0.75) tip='⚖️ <strong>Balanced</strong>: good compromise between gene count and completeness.';
    else tip='📊 <strong>Gene-rich</strong>: maximum loci, acceptable for presence/absence screens.';
  }} else {{
    const rows=parRows();
    const b=rows.find(rj=>rj!==r&&rj.n_genes>=r.n_genes&&rj.avg_completeness>=r.avg_completeness);
    tip = b ? `❌ <strong>Dominated</strong>: θ=${{b.threshold.toFixed(2)}} gives ≥ genes AND ≥ completeness — prefer that threshold.`
            : '⚠️ Sub-optimal: a neighbouring threshold improves at least one metric without worsening the other.';
  }}
  document.getElementById('par-tip').innerHTML='💡 '+tip;
  // Gene list
  const gl = (D.geneLists[parMode]||{{}}) [r.threshold!=null?r.threshold.toFixed(2):''];
  document.getElementById('ev-gene-list').value = gl?(gl.join('\\n')):'(no gene list available)';
}}

function parBuildScatter(si) {{
  const rs=parRows(); if(!rs.length) return;
  const pfr=rs.filter(r=>r.pareto).sort((a,b)=>a.n_genes-b.n_genes);
  const sel=rs[si]||rs[0];
  const traces=[
    {{x:pfr.map(r=>r.n_genes),y:pfr.map(r=>r.avg_completeness),mode:'lines',
      line:{{color:'#dc2626',width:2,dash:'dot'}},name:'Pareto front',hoverinfo:'skip'}},
    {{x:rs.map(r=>r.n_genes),y:rs.map(r=>r.avg_completeness),mode:'markers+text',
      marker:{{size:rs.map(r=>r.pareto?18:11),
               color:rs.map(r=>r.pareto?'#dc2626':'#9ca3af'),
               line:{{width:1.5,color:'#475569'}}}},
      text:rs.map(r=>'θ='+r.threshold.toFixed(2)),
      textposition:'top center',textfont:{{size:10}},
      customdata:rs.map((_,i)=>i),
      hovertemplate:rs.map(r=>'<b>θ='+r.threshold.toFixed(2)+'</b><br>Genes: '+r.n_genes.toLocaleString()+'<br>Avg completeness: '+(r.avg_completeness*100).toFixed(1)+'%<br>Pareto: '+(r.pareto?'✅':'❌')+'<br>Score: '+(r.score||0).toFixed(4)+'<extra></extra>'),
      name:'All thresholds',showlegend:false}},
    {{x:[sel.n_genes],y:[sel.avg_completeness],mode:'markers',
      marker:{{size:24,color:'#2563eb',symbol:'circle-open',line:{{width:3,color:'#1d4ed8'}}}},
      name:'Selected',hoverinfo:'skip',showlegend:false}}
  ];
  Plotly.newPlot('par-scatter',traces,{{
    xaxis:{{title:'Genes retained'}},
    yaxis:{{title:'Avg sample completeness',tickformat:'.0%',range:[0,1.05]}},
    margin:{{l:60,r:10,t:20,b:50}},
    legend:{{orientation:'h',y:1.08}},
    hovermode:'closest',plot_bgcolor:'#fafafa'
  }},{{responsive:true,modeBarButtonsToRemove:['lasso2d','select2d']}});
  document.getElementById('par-scatter').on('plotly_click',d=>{{
    const pt=d.points[0];
    if(pt.curveNumber===1) parSelect(pt.customdata||pt.pointIndex);
  }});
}}

function parBuildSweep(si) {{
  const rs=parRows(); if(!rs.length) return;
  const sel=rs[si]||rs[0];
  const maxG=Math.max(...rs.map(r=>r.n_genes));
  const traces=[
    {{x:rs.map(r=>r.threshold),y:rs.map(r=>r.n_genes),mode:'lines+markers',
      name:'Genes retained',line:{{color:'#2563eb',width:3}},marker:{{size:7}},
      yaxis:'y1',hovertemplate:'θ=%{{x:.2f}} · Genes=%{{y}}<extra></extra>'}},
    {{x:rs.map(r=>r.threshold),y:rs.map(r=>r.avg_completeness),mode:'lines+markers',
      name:'Avg completeness',line:{{color:'#16a34a',width:3,dash:'dash'}},marker:{{size:7}},
      yaxis:'y2',hovertemplate:'θ=%{{x:.2f}} · Compl=%{{y:.3f}}<extra></extra>'}},
    {{x:[sel.threshold,sel.threshold],y:[0,maxG],mode:'lines',
      line:{{color:'#dc2626',width:2,dash:'dot'}},name:'Selected θ',
      yaxis:'y1',hoverinfo:'skip'}}
  ];
  Plotly.newPlot('par-sweep',traces,{{
    xaxis:{{title:'Threshold (θ)',dtick:0.05}},
    yaxis:{{title:'Genes retained',titlefont:{{color:'#2563eb'}},tickfont:{{color:'#2563eb'}}}},
    yaxis2:{{title:'Avg completeness',titlefont:{{color:'#16a34a'}},
             tickfont:{{color:'#16a34a'}},overlaying:'y',side:'right',
             tickformat:'.0%',range:[0,1.05]}},
    margin:{{l:55,r:65,t:20,b:50}},
    legend:{{orientation:'h',y:1.08}},
    hovermode:'x unified',plot_bgcolor:'#fafafa'
  }},{{responsive:true,displayModeBar:false}});
}}

function parSelect(i) {{
  parIdx=i;
  document.getElementById('par-slider').value=i;
  parUpdateStats(i); parBuildScatter(i); parBuildSweep(i);
}}
function parSlider(v) {{ parSelect(parseInt(v)); }}
function parSwitchMode(m) {{
  parMode=m; parIdx=4;
  document.getElementById('par-slider').value=4;
  parUpdateStats(4); parBuildScatter(4); parBuildSweep(4);
}}
function evDl() {{
  const r=parRow(parIdx);
  const gl=(D.geneLists[parMode]||{{}})[r.threshold!=null?r.threshold.toFixed(2):'']||[];
  const blob=new Blob([gl.join('\\n')+'\\n'],{{type:'text/plain'}});
  const a=document.createElement('a'); a.href=URL.createObjectURL(blob);
  a.download='busco_genes_'+parMode+'_t'+(r.threshold||0).toFixed(2)+'.txt';
  a.click(); URL.revokeObjectURL(a.href);
}}
function evCp() {{
  navigator.clipboard.writeText(document.getElementById('ev-gene-list').value).then(()=>{{
    const m=document.getElementById('ev-cmsg'); m.style.display='inline';
    setTimeout(()=>m.style.display='none',1800);
  }});
}}

// ═══════════════════════════════════════════════════════════════════
//  TAB 5 – SMART RECOMMENDATION
// ═══════════════════════════════════════════════════════════════════
function initAdv() {{
  const container = document.getElementById('adv-cards');
  const MODE_LABEL = {{'S':'Mode S (single-copy only)','F':'Mode F (fragmented only)','SF':'Mode SF (single-copy + fragmented)'}};
  const SCENARIOS = [
    {{key:'best_strict', tag:'strict', tagClass:'tag-strict', cardClass:'rec-strict',
      icon:'🔒', title:'Strict Phylogenomics',
      desc:'Maximises average sample completeness → least missing data per locus. Ideal for '
           +'concatenation (RAxML/IQ-TREE) and coalescent-based methods (ASTRAL). '
           +'Recommended when: deeply diverged taxa, rate heterogeneity concerns, '
           +'or reviewer requests minimum-missing-data matrices.'}},
    {{key:'best_score', tag:'balanced', tagClass:'tag-balanced', cardClass:'rec-balanced',
      icon:'⚖️', title:'Balanced (Recommended Default)',
      desc:'Pareto-composite score optimises both gene count and completeness simultaneously. '
           +'Best starting point for most phylogenomic workflows. Use when you are unsure '
           +'which direction to lean, or for initial exploratory analyses.'}},
    {{key:'best_genes', tag:'gene-rich', tagClass:'tag-gene-rich', cardClass:'rec-gene-rich',
      icon:'🧬', title:'Gene-Rich (Maximum Loci)',
      desc:'Retains the most gene families at the cost of more missing data. Useful when: '
           +'closely related taxa with similar lineage, testing gene tree concordance, '
           +'or analyses robust to missing data (e.g. Bayesian or ML with +ASC model). '
           +'Not recommended if samples vary widely in BUSCO completeness.'}},
  ];

  (D.advisorRows||[]).forEach(ar=>{{
    const mLabel = MODE_LABEL[ar.mode]||ar.mode;
    let html='<div class="card"><div class="card-title">'+mLabel+' — '+ar.n_pareto+' Pareto-optimal threshold(s)</div>';
    SCENARIOS.forEach(sc=>{{
      const r=ar[sc.key];
      if(!r) return;
      html += `<div class="rec-card ${{sc.cardClass}}">` +
        `<div class="rec-title">${{sc.icon}} ${{sc.title}}` +
        `  <span class="rec-tag ${{sc.tagClass}}">${{sc.tag}}</span>` +
        `  <code style="font-size:.82rem;background:rgba(0,0,0,.07);padding:2px 7px;border-radius:4px">θ = ${{r.threshold.toFixed(2)}}</code>` +
        `</div>` +
        `<div class="rec-body">${{sc.desc}}<br><br>` +
        `<strong>At θ = ${{r.threshold.toFixed(2)}}:</strong> ` +
        `<strong>${{r.n_genes.toLocaleString()}}</strong> genes retained · ` +
        `avg sample completeness <strong>${{(r.avg_completeness*100).toFixed(1)}}%</strong> · ` +
        `${{r.total_cells.toLocaleString()}} data cells · ` +
        `score <strong>${{(r.score).toFixed(4)}}</strong>` +
        `</div></div>`;
    }});
    html += '</div>';
    container.innerHTML += html;
  }});

  if(!(D.advisorRows||[]).length) {{
    container.innerHTML='<div class="card"><em>No Pareto sweep data available. Re-run with --smart or --plot-pareto to generate sweep files, then re-run --explore_view.</em></div>';
  }}

  // All thresholds table
  const tableContainer = document.getElementById('adv-table');
  const firstMode = MODES[0];
  const rows = D.paretoSweeps[firstMode]||[];
  if(!rows.length) return;
  let t='<div style="overflow-x:auto"><table style="width:100%;border-collapse:collapse;font-size:.84rem">';
  t+='<tr style="background:#f8fafc"><th style="padding:8px 10px;text-align:left;border-bottom:2px solid #e2e8f0">θ</th>';
  MODES.forEach(m=>{{ t+=`<th style="padding:8px 10px;text-align:right;border-bottom:2px solid #e2e8f0">Genes (${m})</th><th style="padding:8px 10px;text-align:right;border-bottom:2px solid #e2e8f0">Compl (${m})</th><th style="padding:8px 10px;text-align:right;border-bottom:2px solid #e2e8f0">Pareto</th>`; }});
  t+='</tr>';
  rows.forEach(r=>{{
    const rowBg = r.pareto ? '#f0f9ff' : '#fff';
    t+=`<tr style="background:${{rowBg}};border-bottom:1px solid #f1f5f9">`;
    t+=`<td style="padding:7px 10px;font-weight:700">${{r.threshold.toFixed(2)}}</td>`;
    MODES.forEach(m=>{{
      const mr=(D.paretoSweeps[m]||[]).find(x=>Math.abs(x.threshold-r.threshold)<0.001)||{{}};
      t+=`<td style="padding:7px 10px;text-align:right">${{(mr.n_genes||0).toLocaleString()}}</td>`;
      t+=`<td style="padding:7px 10px;text-align:right">${{mr.avg_completeness!=null?(mr.avg_completeness*100).toFixed(1)+'%':'—'}}</td>`;
      t+=`<td style="padding:7px 10px;text-align:right">${{mr.pareto?'★':'-'}}</td>`;
    }});
    t+='</tr>';
  }});
  t+='</table></div>';
  tableContainer.innerHTML=t;
}}
</script>
</body>
</html>"""

    out_path.write_text(html, encoding="utf-8")
    info(f"Explore View HTML        : {out_path}")
    substep("Open in any modern browser.  Requires internet to load Plotly charts.")
    substep("Tabs: Overview · Gene Coverage · Sample Completeness · Pareto · Recommendation")


# ──────────────────────────────────────────────────────────────────────────────
# FASTA I/O
# ──────────────────────────────────────────────────────────────────────────────

def read_fasta_first(path: Path) -> Tuple[str, str]:
    hdr: Optional[str] = None; seqs: List[str] = []
    with open(path, "rt", encoding="utf-8", errors="replace") as f:
        for ln in f:
            ln = ln.rstrip("\n\r")
            if not ln: continue
            if ln.startswith(">"):
                if hdr is not None: break
                hdr = ln[1:].strip()
            elif hdr is not None: seqs.append(ln.strip())
    if hdr is None: raise ValueError(f"No FASTA header in {path}")
    return hdr, "".join(seqs)

def fmt_header(sp: str, bid: str, style: str) -> str:
    s = sanitize_fasta_header(sp)
    return f"{s}|{bid}" if style == "pipe" else f"{s}.{bid}"

def write_fasta(records: List[Tuple[str, str]], path: Path) -> int:
    if not records: return 0
    with open(path, "w", encoding="utf-8") as f:
        for hdr, seq in records:
            f.write(f">{hdr}\n")
            for i in range(0, len(seq), 80): f.write(seq[i:i+80]+"\n")
    return len(records)

def pick_seq_path(sp: SpeciesData, busco_id: str,
                  mode: str, kind: str) -> Optional[Path]:
    st_map = sp.buscos.get(busco_id)
    if not st_map: return None
    for st in _status_order_for_mode(mode):
        obj = st_map.get(st)
        if not obj: continue
        if kind == "protein"    and has_protein(obj):    return obj.faa
        if kind == "nucleotide" and has_nucleotide(obj): return obj.fna
    return None


# ──────────────────────────────────────────────────────────────────────────────
# EXTRACTION – GPS
# ──────────────────────────────────────────────────────────────────────────────

def extract_gps(sd: SpeciesDict, all_ids: Set[str], mode: str, threshold: float,
                out_dir: Path, seqtype: str, header_style: str = "dot",
                show_progress: bool = True) -> Dict[str, int]:
    kept = [(bid, gene_coverage(sd, bid, mode, seqtype))
            for bid in sorted(all_ids)
            if gene_coverage(sd, bid, mode, seqtype) >= threshold]
    if not kept:
        warn(f"GPS {mode}: no genes passed threshold={threshold:.0%}"); return {}
    tag  = f"gps_{mode}_cov{int(threshold*100)}"
    base = ensure_dir(out_dir / tag)
    with open(base/"kept_buscos.tsv","w") as f:
        f.write("BUSCO_ID\tcoverage\n")
        for bid, c in kept: f.write(f"{bid}\t{c:.4f}\n")
    kinds = (["protein","nucleotide"] if seqtype=="both"
             else ["protein"] if seqtype=="protein" else ["nucleotide"])
    manifest: List[List[str]] = []; stats: Dict[str, int] = {}
    for kind in kinds:
        ext  = "faa" if kind=="protein" else "fna"
        odir = ensure_dir(base/ext); nf = 0
        for bid, _ in tqdm(kept, desc=f"    GPS {mode} {kind}",
                           disable=not show_progress):
            records = []
            for sp_name in sorted(sd):
                p = pick_seq_path(sd[sp_name], bid, mode, kind)
                if p:
                    try:
                        _, seq = read_fasta_first(p)
                        records.append((fmt_header(sp_name, bid, header_style), seq))
                        manifest.append([bid, sp_name, kind, "extracted", str(p)])
                    except Exception as e:
                        warn(f"Read error {p}: {e}")
                        manifest.append([bid, sp_name, kind, "read_error", str(p)])
                else:
                    manifest.append([bid, sp_name, kind, "missing", ""])
            if records: write_fasta(records, odir/f"{bid}.{ext}"); nf += 1
        stats[kind] = nf
        substep(f"GPS {mode} {kind}: {nf} gene FASTA files → {odir}")
    _write_manifest(manifest, base/"extraction_manifest.tsv")
    return stats


# ──────────────────────────────────────────────────────────────────────────────
# EXTRACTION – SCH
# ──────────────────────────────────────────────────────────────────────────────

def extract_sch(sd: SpeciesDict, all_ids: Set[str], mode: str, threshold: float,
                out_dir: Path, seqtype: str, header_style: str = "dot",
                show_progress: bool = True) -> Dict[str, int]:
    kept = [(n, sample_completeness(sd[n], all_ids, mode, seqtype))
            for n in sorted(sd)
            if sample_completeness(sd[n], all_ids, mode, seqtype) >= threshold]
    if not kept:
        warn(f"SCH {mode}: no samples passed threshold={threshold:.0%}"); return {}
    tag  = f"sch_{mode}_comp{int(threshold*100)}"
    base = ensure_dir(out_dir / tag)
    with open(base/"kept_species.tsv","w") as f:
        f.write("species\tcompleteness\n")
        for n, c in kept: f.write(f"{n}\t{c:.4f}\n")
    kinds = (["protein","nucleotide"] if seqtype=="both"
             else ["protein"] if seqtype=="protein" else ["nucleotide"])
    manifest: List[List[str]] = []; stats: Dict[str, int] = {}
    for kind in kinds:
        ext  = "faa" if kind=="protein" else "fna"
        odir = ensure_dir(base/ext); nf = 0
        for sp_name, _ in tqdm(kept, desc=f"    SCH {mode} {kind}",
                               disable=not show_progress):
            sp = sd[sp_name]; records = []
            for bid in sorted(all_ids):
                p = pick_seq_path(sp, bid, mode, kind)
                if p:
                    try:
                        _, seq = read_fasta_first(p)
                        records.append((fmt_header(sp_name, bid, header_style), seq))
                        manifest.append([bid, sp_name, kind, "extracted", str(p)])
                    except Exception as e:
                        warn(f"Read error {p}: {e}")
                        manifest.append([bid, sp_name, kind, "read_error", str(p)])
                else:
                    manifest.append([bid, sp_name, kind, "missing", ""])
            if records: write_fasta(records, odir/f"{sp_name}.{ext}"); nf += 1
        stats[kind] = nf
        substep(f"SCH {mode} {kind}: {nf} sample FASTA files → {odir}")
    _write_manifest(manifest, base/"extraction_manifest.tsv")
    return stats

def _write_manifest(rows: List[List[str]], path: Path) -> None:
    with open(path, "w", newline="") as f:
        w = csv.writer(f, delimiter="\t")
        w.writerow(["busco_id","sample","seq_kind","status","source_path"])
        w.writerows(rows)
    substep(f"Extraction manifest: {path}")


# ──────────────────────────────────────────────────────────────────────────────
# STATISTICS
# ──────────────────────────────────────────────────────────────────────────────

def print_statistics(sd: SpeciesDict, all_ids: Set[str],
                     modes: List[str], seqtype: str,
                     gpt_list: List[float], sct_list: List[float]) -> None:
    print("\n  ╔══════════════════════════════════════════════════════════════╗")
    print("  ║                    Summary Statistics                       ║")
    print("  ╚══════════════════════════════════════════════════════════════╝")
    print(f"  Samples       : {len(sd)}")
    print(f"  Unique BUSCOs : {len(all_ids)}")
    print(f"  Sequence type : {seqtype}")
    lgs: Dict[str, int] = {}
    for sp in sd.values():
        k = sp.lineage or "unknown"; lgs[k] = lgs.get(k, 0) + 1
    print("  Lineages:")
    for lg, c in sorted(lgs.items()): print(f"    • {lg}: {c} samples")
    for mode in modes:
        print(f"\n  ── Mode {mode} ──")
        covs  = [gene_coverage(sd, bid, mode, seqtype) for bid in all_ids]
        comps = [sample_completeness(sp, all_ids, mode, seqtype) for sp in sd.values()]
        if covs:  print(f"    Mean gene coverage      : {sum(covs)/len(covs):.4f}")
        if comps:
            print(f"    Mean sample completeness: {sum(comps)/len(comps):.4f}")
            print(f"    Completeness range      : {min(comps):.4f} – {max(comps):.4f}")
        print(f"\n    Gene coverage threshold breakdown:")
        for t in sorted(set(gpt_list)):
            n = sum(1 for bid in all_ids if gene_coverage(sd,bid,mode,seqtype)>=t)
            print(f"      >= {t*100:5.1f}%  →  {n:>6} genes")
        print(f"\n    Sample completeness threshold breakdown:")
        for t in sorted(set(sct_list)):
            n = sum(1 for sp in sd.values()
                    if sample_completeness(sp,all_ids,mode,seqtype)>=t)
            print(f"      >= {t*100:5.1f}%  →  {n:>6} samples")
    print()


# ──────────────────────────────────────────────────────────────────────────────
# REPORT WRITERS
# ──────────────────────────────────────────────────────────────────────────────

def write_run_report(args: argparse.Namespace, out_dir: Path) -> None:
    with open(out_dir/"run_report.txt","w") as f:
        f.write(f"BUSCO Phylogenomics Gene Set Extractor v{VERSION}\n")
        f.write(f"Run date : {now_ts()}\n{'='*70}\n\n")
        f.write("Command:\n  "+" ".join(sys.argv)+"\n\nParameters:\n")
        for k,v in sorted(vars(args).items()): f.write(f"  {k}: {v}\n")

def write_run_selection(chosen: Dict[str, str], sd: SpeciesDict,
                        out_dir: Path) -> None:
    with open(out_dir/"chosen_run_dirs.tsv","w") as f:
        f.write("species\toriginal_name\trun_dir\tlineage\tbusco_path\n")
        for n in sorted(sd):
            sp=sd[n]
            f.write(f"{n}\t{sp.original_name}\t{sp.run_dir}\t{sp.lineage}\t{sp.busco_path}\n")

def write_lineage_report(sd: SpeciesDict, lg: Dict[str, List[str]],
                         out_dir: Path) -> None:
    with open(out_dir/"lineage_report.tsv","w") as f:
        f.write("lineage\tn_species\tspecies_list\n")
        for l in sorted(lg):
            sl=sorted(lg[l]); f.write(f"{l}\t{len(sl)}\t{','.join(sl)}\n")

def write_name_mapping(sd: SpeciesDict, out_dir: Path) -> None:
    with open(out_dir/"name_mapping.tsv","w") as f:
        f.write("final_name\toriginal_name\tlineage\tbusco_version\n")
        for n in sorted(sd):
            sp=sd[n]
            f.write(f"{n}\t{sp.original_name}\t{sp.lineage}\t{sp.busco_version}\n")
    info(f"Name mapping: {out_dir/'name_mapping.tsv'}")

def write_species_summary(sd: SpeciesDict, all_ids: Set[str],
                          out_dir: Path) -> None:
    with open(out_dir/"species_summary.tsv","w") as f:
        f.write("species\toriginal_name\tlineage\tversion\t"
                "n_total\tn_complete\tn_single\tn_duplicated\t"
                "n_fragmented\tn_missing\tn_buscos_with_seq\n")
        for n in sorted(sd):
            sp=sd[n]
            nseq=sum(1 for bid in sp.buscos
                     if any(has_any(sp.buscos[bid].get(st)) for st in ["S","F","D"]))
            f.write(f"{n}\t{sp.original_name}\t{sp.lineage}\t{sp.busco_version}\t"
                    f"{sp.n_total}\t{sp.n_complete}\t{sp.n_single}\t"
                    f"{sp.n_duplicated}\t{sp.n_fragmented}\t{sp.n_missing}\t{nseq}\n")
    info(f"Species summary: {out_dir/'species_summary.tsv'}")

def write_busco_id_lists(sd: SpeciesDict, all_ids: Set[str],
                         out_dir: Path) -> None:
    d = ensure_dir(out_dir/"busco_lists")
    for n, sp in sd.items():
        with open(d/f"{n}.busco_ids.txt","w") as f:
            for bid in sorted(sp.buscos):
                f.write(f"{bid}\t{get_display_status(sp,bid)}\n")
    s_set: Set[str] = set(); f_set: Set[str] = set(); d_set: Set[str] = set()
    for sp in sd.values():
        for bid, sm in sp.buscos.items():
            if sm.get("S") and has_any(sm["S"]): s_set.add(bid)
            if sm.get("F") and has_any(sm["F"]): f_set.add(bid)
            if sm.get("D") and has_any(sm["D"]): d_set.add(bid)
    for name, ids in [("All_buscoid_S.txt",  s_set),
                      ("All_buscoid_F.txt",  f_set),
                      ("All_buscoid_SF.txt", s_set|f_set),
                      ("All_buscoid_D.txt",  d_set)]:
        with open(out_dir/name,"w") as f:
            f.writelines(bid+"\n" for bid in sorted(ids))
        info(f"BUSCO ID list: {out_dir/name}  (n={len(ids)})")

def generate_heatmap(sd: SpeciesDict, all_ids: Set[str],
                     out_dir: Path, mode: str, seqtype: str) -> None:
    try:
        import matplotlib; matplotlib.use("Agg")
        import matplotlib.pyplot as plt
    except ImportError:
        warn("matplotlib not available – skipping heatmap"); return
    names=sorted(sd); bids=sorted(all_ids)
    if not bids or not names: return
    mat=[[1 if can_extract(sd[n],bid,mode,seqtype) else 0 for n in names] for bid in bids]
    fw=min(20,max(8,len(names)*0.3)); fh=min(30,max(8,len(bids)*0.05))
    fig,ax=plt.subplots(figsize=(fw,fh))
    ax.imshow(mat,cmap="viridis",aspect="auto",interpolation="none")
    ax.set_xticks(range(len(names))); ax.set_xticklabels(names,rotation=90,fontsize=8)
    if len(bids)<=50:
        ax.set_yticks(range(len(bids))); ax.set_yticklabels(bids,fontsize=6)
    else:
        ax.set_yticks([])
    ax.set_title(f"BUSCO presence ({mode},{seqtype})")
    ax.set_xlabel("Samples"); ax.set_ylabel("BUSCO IDs")
    plt.tight_layout()
    out=out_dir/f"heatmap_{mode}.png"
    plt.savefig(out,dpi=150,bbox_inches="tight"); plt.close()
    info(f"Heatmap: {out}")


# ──────────────────────────────────────────────────────────────────────────────
# CLI
# ──────────────────────────────────────────────────────────────────────────────

def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(
        prog="busco_phylo_extract_v7.py",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        description=textwrap.dedent(f"""\
╔══════════════════════════════════════════════════════════════════════╗
║    BUSCO Phylogenomics Gene Set Extractor  v{VERSION:<24}       ║
╚══════════════════════════════════════════════════════════════════════╝

PURPOSE
  Extract BUSCO orthologs from hundreds of samples for phylogenomics.
  Compatible with BUSCO v3/v4/v5/v6.  Standard workflow:
    MAFFT (align)  →  trimAl (trim)  →  IQ-TREE (tree)

QUICK START
  # 1. Scan lineages only (see what datasets your samples used)
  %(prog)s -i /data/busco_results/ -o scan/ --scan-only

  # 2. Interactive HTML Pareto decision tool (choose threshold visually)
  %(prog)s --mapping samples.tsv -o out/ --lid metazoa_odb10 \\
           --plot-pareto --no-extract

  # 3. Console smart recommendation + full extraction
  %(prog)s --mapping samples.tsv -o out/ --lid metazoa_odb10 --smart

  # 4. Targeted: single-copy protein, 90%% gene coverage, GPS only
  %(prog)s --mapping samples.tsv -o out/ --lid metazoa_odb10 \\
           --modes S --seqtype protein -gpt 0.9 --only-gps

  # 5. Directory scan, no mapping file needed
  %(prog)s -i /data/busco/ -o out/ --lid insecta_odb10 --smart

  # 6. [NEW] Explore a previous output directory visually (no BUSCO data needed)
  %(prog)s -ev /previous/output/dir/
  %(prog)s -ev /previous/output/busco_status_matrix.tsv
  %(prog)s -ev /previous/output/ --smart   # add sweeps then open explorer

EXPLORE VIEW  (-ev / --explore_view)
  Read finished run outputs and open a rich 5-tab HTML dashboard:
    📊 Overview       – stats cards, coverage histogram, source file audit
    📈 Gene Coverage  – interactive histogram + CDF with threshold sliders
    🧬 Completeness   – ranked per-sample completeness bar chart
    🎯 Pareto         – full interactive Pareto scatter + sweep selector
    💡 Recommendation – smart decision guide with 3 strategies per mode
  The HTML is written alongside the source as 'explore_view.html'.
  Can be used STANDALONE (no -i or --mapping required).
"""),
        epilog=textwrap.dedent("""\
EXTRACTION MODES  (--modes)
  S    Single-copy complete BUSCOs only.  Gold standard for phylogenomics.
  F    Fragmented BUSCOs only (rarely used alone).
  SF   Single-copy + fragmented.  More loci but lower quality alignment cols.

SEQUENCE TYPES  (--seqtype)
  protein      .faa  Faster alignment; avoids 3rd-codon wobble.
  nucleotide   .fna  More variable; better for closely related taxa.
  both         Extract both types simultaneously.

GPS vs SCH
  GPS  (--only-gps): one FASTA per gene   → per-gene trees → ASTRAL
  SCH  (--only-sch): one FASTA per sample → super-matrix / concatenation

EXPLORE VIEW  (-ev / --explore_view PATH)
  PATH may be:
    • a previous -o output directory (auto-detects busco_status_matrix.tsv,
      busco_binary_*.tsv, species_summary.tsv, pareto_sweep_*.tsv)
    • a specific .tsv file (busco_status_matrix or busco_binary_*)
  Writes  <PATH>/explore_view.html  (or alongside the file for single files).
  No -i or --mapping required.  Requires internet to load Plotly charts.

  5 interactive tabs:
    1. Overview        – headline stats, quick coverage bar chart
    2. Gene Coverage   – histogram (20 bins) + CDF curve per mode
    3. Completeness    – ranked horizontal bars, green/yellow/orange/red
    4. Pareto Advisor  – scatter + sweep + slider + gene list download
    5. Recommendation  – 3 strategies (Strict · Balanced · Gene-rich)
                         with rationale and per-threshold stats table

COLUMN SYNONYMS (for mapping file headers, '#' optional, case-insensitive)
  #original_name   : orig, sample, sample_name, input_name, directory ...
  #new_name        : new, rename, label, alias, taxon, species ...
  #abs_path_…      : path, abs_path, busco_path, busco_dir, result_dir ...
  #lineage_…       : lineage, lid, odb, dataset, db, lineage_name ...

OUTPUT LAYOUT
  out/
    run_report.txt              Command + all parameters
    busco_status_matrix.tsv     S/D/F/M matrix  (genes × samples)
    busco_binary_{S,F,SF}.tsv   0/1 binary presence matrices
    sequence_availability.tsv   Protein/nucleotide audit per gene×sample
    species_summary.tsv         Per-sample BUSCO statistics
    name_mapping.tsv            original_name → final_name lookup
    All_buscoid_{S,F,SF,D}.txt  Global BUSCO ID lists
    pareto_interactive.html     Original Pareto tool (--plot-pareto)
    explore_view.html           Rich 5-tab dashboard (-ev)
    pareto_sweep_{mode}.tsv     Pareto table (--smart)
    extract_S/F/SF/
      gps_*/faa/*.faa  fna/*.fna
      sch_*/faa/*.faa  fna/*.fna

DEPENDENCIES
  Required : Python ≥ 3.7
  Optional : numpy, pandas, tqdm, matplotlib
  HTML     : Plotly.js loaded from CDN  (internet needed for charts)
"""),
    )

    # ── Input (all optional when -ev is standalone) ───────────────────────────
    ig = p.add_argument_group("Input (one of -i / --mapping / -ev required)")
    ig.add_argument("-i", "--input_dir", metavar="DIR",
                    help="Directory containing BUSCO result subdirectories.")
    ig.add_argument("--mapping", dest="mapping_file", metavar="FILE",
                    help="Flexible mapping file (headers start with '#').")
    ig.add_argument("-ev", "--explore_view", dest="explore_view", metavar="PATH",
                    help="[NEW v7.2] Load a previous output directory (or a single "
                         "matrix .tsv) and generate a rich interactive HTML dashboard "
                         "(explore_view.html).  Can be used ALONE without -i/--mapping.")

    p.add_argument("-o", "--output_dir", metavar="DIR",
                   help="Output directory (required unless using -ev standalone).")

    # Naming
    p.add_argument("--keep-dots", dest="replace_dot", action="store_false",
                   default=True,
                   help="Keep dots in sample names.")
    p.add_argument("--replace-char", default="_", metavar="CHAR",
                   help="Replacement character for dots (default: '_').")

    # Lineage / run selection
    p.add_argument("--lid", "--lineage_datasets", dest="lineage", default="",
                   metavar="LINEAGE",
                   help="Global lineage filter, e.g. metazoa_odb10.")
    p.add_argument("--run_name", default="", metavar="RUN_DIR",
                   help="Force a specific run_* directory name for all samples.")
    p.add_argument("--scan-only", dest="scan_only", action="store_true",
                   help="Scan and report lineages only.")

    # Extraction modes
    p.add_argument("--modes", nargs="+", default=["S","F","SF"],
                   choices=["S","F","SF"],
                   help="Extraction modes: S, F, SF.  Default: all three.")
    p.add_argument("--seqtype", default="both",
                   choices=["protein","nucleotide","both"],
                   help="Sequence type.  Default: both.")
    p.add_argument("-gpt","--gene_presence_threshold", nargs="*", type=float,
                   default=[0.5,0.7,0.8,0.9,0.95,1.0], metavar="FLOAT",
                   help="GPS: gene retained if present in >= this fraction of samples.")
    p.add_argument("-sct","--species_completeness_threshold", nargs="*", type=float,
                   default=[0.5,0.7,0.8,0.9,0.95,1.0], metavar="FLOAT",
                   help="SCH: sample retained if it has >= this fraction of genes.")

    mex = p.add_mutually_exclusive_group()
    mex.add_argument("--only-gps", dest="only_gps", action="store_true")
    mex.add_argument("--only-sch", dest="only_sch", action="store_true")
    p.add_argument("--no-extract", dest="no_extract", action="store_true",
                   help="Skip sequence extraction (matrix + stats + Pareto only).")

    # Smart / Pareto
    p.add_argument("--smart", action="store_true",
                   help="Run Pareto sweep and print recommendation table.")
    p.add_argument("--plot-pareto", dest="plot_pareto", action="store_true",
                   help="Generate pareto_interactive.html.")

    # Filtering
    p.add_argument("--samples_list", metavar="FILE",
                   help="Text file of sample names to include.")
    p.add_argument("--suffx_subdir", default="", metavar="SUFFIX",
                   help="Include only subdirectories ending with this suffix.")
    p.add_argument("--include_dup", action="store_true",
                   help="Scan Duplicated BUSCOs for the status matrix audit.")
    p.add_argument("--use_full_table", action="store_true",
                   help="Register all BUSCO IDs from full_table.tsv.")

    # Output
    p.add_argument("--plot", action="store_true",
                   help="Generate static heatmap PNG (requires matplotlib).")
    p.add_argument("--no-matrix", dest="no_matrix", action="store_true",
                   help="Skip writing matrix TSV files.")
    p.add_argument("--header_style", choices=["dot","pipe"], default="dot",
                   help="FASTA header style (dot or pipe).")
    p.add_argument("--no-progress", dest="no_progress", action="store_true",
                   help="Suppress tqdm progress bars.")
    p.add_argument("--overwrite", action="store_true",
                   help="Delete and recreate existing output directory.")
    p.add_argument("--version", action="version",
                   version=f"%(prog)s {VERSION}")
    return p.parse_args()


# ──────────────────────────────────────────────────────────────────────────────
# MAIN
# ──────────────────────────────────────────────────────────────────────────────

def main() -> None:
    args = parse_args()
    t0 = time()

    print()
    print("╔══════════════════════════════════════════════════════════════════╗")
    print(f"║  BUSCO Phylogenomics Gene Set Extractor  v{VERSION:<28}  ║")
    print(f"║  {now_ts():<64}  ║")
    print("╚══════════════════════════════════════════════════════════════════╝")

    if args.plot_pareto:
        args.smart = True

    # ════════════════════════════════════════════════════════════════════════
    # EXPLORE VIEW STANDALONE MODE
    # When -ev is given without -i / --mapping, skip the full pipeline and
    # just load the matrices and generate the HTML dashboard.
    # ════════════════════════════════════════════════════════════════════════
    if args.explore_view and not args.input_dir and not args.mapping_file:
        step("Explore View – standalone mode")
        ev_path = Path(args.explore_view).resolve()
        if not ev_path.exists():
            die(f"--explore_view path does not exist: {ev_path}")

        substep(f"Source path: {ev_path}")
        ed = load_explore_data(ev_path)

        if not ed["genes"]:
            die("No gene data could be loaded from the specified path.\n"
                "       Expected: busco_status_matrix.tsv or busco_binary_*.tsv")

        # Determine output path
        if ev_path.is_file():
            out_html = ev_path.parent / "explore_view.html"
        else:
            out_html = ev_path / "explore_view.html"

        generate_explore_html(ed, out_html)

        elapsed = time() - t0
        print()
        print("╔══════════════════════════════════════════════════════════════════╗")
        print(f"║  EXPLORE VIEW DONE  –  {elapsed:5.1f}s                                   ║")
        print(f"║  HTML  →  {str(out_html)[:54]:<54}  ║")
        print("╚══════════════════════════════════════════════════════════════════╝")
        print()
        return

    # ════════════════════════════════════════════════════════════════════════
    # NORMAL PIPELINE  (requires -i or --mapping)
    # ════════════════════════════════════════════════════════════════════════
    if not args.input_dir and not args.mapping_file:
        die("One of -i / --mapping / -ev (standalone) is required.\n"
            "       Run with --help for usage.")

    if not args.output_dir:
        die("--output_dir / -o is required when running the full pipeline.")

    # ── Output directory ─────────────────────────────────────────────────────
    step("Preparing output directory")
    out_dir = Path(args.output_dir).resolve()
    if out_dir.exists() and any(out_dir.iterdir()):
        if args.overwrite:
            substep(f"Removing existing output: {out_dir}")
            shutil.rmtree(out_dir)
        else:
            die(f"Output directory already exists: {out_dir}\n"
                f"       Use --overwrite to replace it.")
    ensure_dir(out_dir)
    substep(f"Output directory : {out_dir}")
    substep(f"numpy   : {'available' if HAS_NUMPY  else 'not installed'}")
    substep(f"pandas  : {'available' if HAS_PANDAS else 'not installed'}")
    substep(f"tqdm    : {'available' if HAS_TQDM   else 'not installed'}")
    substep(f"Modes   : {args.modes}   seqtype : {args.seqtype}")
    substep(f"GPS thresholds  : {args.gene_presence_threshold}")
    substep(f"SCH thresholds  : {args.species_completeness_threshold}")
    if args.lineage: substep(f"Lineage filter  : {args.lineage}")
    if args.include_dup:
        warn("--include_dup: D BUSCOs scanned for audit ONLY; NOT in extractions.")
    write_run_report(args, out_dir)

    # ── Build entry list ─────────────────────────────────────────────────────
    step("Loading sample entries")
    show = not args.no_progress and HAS_TQDM

    if args.mapping_file:
        input_dir_arg = (Path(args.input_dir).resolve()
                         if getattr(args, "input_dir", None) else None)
        entries = parse_mapping_file(
            args.mapping_file,
            input_dir=input_dir_arg,
            replace_dot=args.replace_dot,
            replace_char=args.replace_char,
            global_lineage=args.lineage,
        )
    else:
        idir = Path(args.input_dir).resolve()
        sl   = load_samples_list(args.samples_list) if args.samples_list else None
        entries = build_directory_entries(idir, args.suffx_subdir, sl)
        substep(f"Found {len(entries)} candidate directories in {idir}")

    # ── Lineage scan ─────────────────────────────────────────────────────────
    step("Scanning available lineages")
    all_lineages = scan_all_lineages(entries)
    substep(f"Detected {len(all_lineages)} unique lineage(s):")
    for lg, d in sorted(all_lineages.items(), key=lambda x: -len(x[1]["species"])):
        substep(f"  {lg}: {len(d['species'])} samples,  {d['n_busco']} BUSCOs")

    if args.scan_only:
        rf = out_dir / "all_lineages_report.tsv"
        with open(rf,"w") as f:
            f.write("lineage\tn_samples\tn_busco\n")
            for lg, d in sorted(all_lineages.items()):
                f.write(f"{lg}\t{len(d['species'])}\t{d['n_busco']}\n")
        substep(f"Lineage report: {rf}")
        print("\n  --scan-only mode: done.  No extraction performed.\n")
        return

    if not args.lineage and len(all_lineages) > 1:
        warn("Multiple lineages detected!  Use --lid to select one for phylogenomics.")
        warn("Continuing with all lineages merged (use with caution).")

    # ── Scan BUSCO data ───────────────────────────────────────────────────────
    step("Scanning BUSCO sequence data")
    sd, all_ids, chosen, lg_groups = scan_from_entries(
        entries, global_run_name=args.run_name,
        global_lineage=args.lineage,
        include_dup=args.include_dup,
        use_full_table=args.use_full_table,
        show_progress=show,
    )
    substep(f"Loaded {len(sd)} samples,  {len(all_ids)} unique BUSCO IDs")

    if not sd:
        die("No valid BUSCO data found.")

    # ── Core reports ─────────────────────────────────────────────────────────
    step("Writing core reports")
    write_run_selection(chosen, sd, out_dir)
    write_lineage_report(sd, lg_groups, out_dir)
    write_name_mapping(sd, out_dir)
    write_species_summary(sd, all_ids, out_dir)
    write_busco_id_lists(sd, all_ids, out_dir)

    # ── Matrices ─────────────────────────────────────────────────────────────
    if not args.no_matrix:
        step("Building matrices  (rows = genes, columns = samples)")
        substep("S-only  binary: S=1, D/F/M=0  (strict; gold standard)")
        substep("S+F     binary: S/F=1, D/M=0  (fragmented adds loci, lower quality)")
        substep("D status always excluded from extraction (multi-copy artefact risk)")
        raw_matrix = build_raw_status_matrix(sd, all_ids)
        save_matrices(sd, all_ids, raw_matrix, out_dir, args.modes, args.seqtype)

    # ── Statistics ────────────────────────────────────────────────────────────
    step("Computing summary statistics")
    gpt_list = sorted(set(args.gene_presence_threshold))
    sct_list = sorted(set(args.species_completeness_threshold))
    print_statistics(sd, all_ids, args.modes, args.seqtype, gpt_list, sct_list)

    # ── Smart / Pareto ────────────────────────────────────────────────────────
    if args.smart:
        step("Pareto-optimal threshold analysis")
        smart_recommendation(sd, all_ids, args.modes, args.seqtype, out_dir)

    if args.plot_pareto:
        step("Generating interactive Pareto HTML (v7.1 style)")
        generate_pareto_html(sd, all_ids, args.modes, args.seqtype, out_dir)

    # ── Explore View (pipeline mode) ─────────────────────────────────────────
    if args.explore_view:
        step("Generating Explore View HTML (v7.2 dashboard)")
        # At this point matrices are freshly written; load them from out_dir
        substep(f"Loading matrices from: {out_dir}")
        ev_path   = out_dir
        ed        = load_explore_data(ev_path)
        out_html  = out_dir / "explore_view.html"
        generate_explore_html(ed, out_html)

    # ── Static heatmap ────────────────────────────────────────────────────────
    if args.plot:
        step("Generating static presence/absence heatmap")
        for m in args.modes:
            generate_heatmap(sd, all_ids, out_dir, m, args.seqtype)

    # ── Extraction ────────────────────────────────────────────────────────────
    if args.no_extract:
        substep("--no-extract: sequence extraction skipped.")
    else:
        step("Extracting BUSCO sequences")
        for mode in args.modes:
            mode_dir = ensure_dir(out_dir / f"extract_{mode}")
            substep(f"── Mode {mode} ──")
            if not args.only_sch:
                substep("GPS (gene-centric): one FASTA per BUSCO gene")
                for thresh in gpt_list:
                    extract_gps(sd, all_ids, mode, thresh, mode_dir,
                                args.seqtype, args.header_style, show)
            if not args.only_gps:
                substep("SCH (sample-centric): one FASTA per sample")
                for thresh in sct_list:
                    extract_sch(sd, all_ids, mode, thresh, mode_dir,
                                args.seqtype, args.header_style, show)

    # ── Done ─────────────────────────────────────────────────────────────────
    elapsed = time() - t0
    print()
    print("╔══════════════════════════════════════════════════════════════════╗")
    print(f"║  DONE  –  elapsed {elapsed:6.1f}s                                        ║")
    print(f"║  Output → {str(out_dir)[:54]:<54}  ║")
    if args.plot_pareto:
        print("║  pareto_interactive.html  – original Pareto selector             ║")
    if args.explore_view:
        print("║  explore_view.html        – rich 5-tab dashboard (NEW v7.2)      ║")
    print("╚══════════════════════════════════════════════════════════════════╝")
    print()


if __name__ == "__main__":
    main()

