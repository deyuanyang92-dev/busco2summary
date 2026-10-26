"""Parsers for BUSCO (v3-v6) and compleasm result files.

Every parser is tolerant: unknown lines are ignored and missing fields are
returned as ``None`` so that one odd sample never stops a large batch.
"""

from __future__ import annotations

import gzip
import json
import re
from dataclasses import dataclass, field
from pathlib import Path
from typing import Dict, Iterator, List, Optional, Tuple

# full_table "Status" words -> status code.  compleasm uses Single /
# Interspaced / Retrocopy; BUSCO v3-v6 use Complete / Duplicated / Fragmented / Missing.
STATUS_WORDS: Dict[str, str] = {
    "complete": "S",
    "single": "S",
    "complete and single-copy": "S",
    "duplicated": "D",
    "complete and duplicated": "D",
    "multi": "D",
    "fragmented": "F",
    "interspersed": "F",   # compleasm "I": counted as fragmented
    "interspaced": "F",    # compleasm >= 0.2 spelling
    "retrocopy": "D",      # compleasm folds R into D unless --retrocopy
    "incomplete": "F",
    "missing": "M",
}

#: Precedence when one BUSCO has rows with different statuses.
STATUS_RANK = {"S": 0, "D": 1, "F": 2, "M": 3}


def open_text(path: Path):
    """Open plain or gzip-compressed text."""
    path = Path(path)
    if path.suffix == ".gz":
        return gzip.open(path, "rt", encoding="utf-8", errors="replace")
    return open(path, "rt", encoding="utf-8", errors="replace")


def _num(x) -> Optional[float]:
    if x is None:
        return None
    if isinstance(x, (int, float)):
        return float(x)
    s = str(x).strip().rstrip("%").replace(",", "")
    try:
        return float(s)
    except ValueError:
        return None


def _int(x) -> Optional[int]:
    v = _num(x)
    return None if v is None else int(round(v))


# --------------------------------------------------------------------------
# short_summary (text)
# --------------------------------------------------------------------------

_SUMMARY_COUNT_PATTERNS = [
    (r"^\s*(\d+)\s+Complete BUSCOs", "n_complete"),
    (r"^\s*(\d+)\s+Complete and single-copy", "n_single"),
    (r"^\s*(\d+)\s+Complete and duplicated", "n_duplicated"),
    (r"^\s*(\d+)\s+Fragmented BUSCOs", "n_fragmented"),
    (r"^\s*(\d+)\s+Missing BUSCOs", "n_missing"),
    (r"^\s*(\d+)\s+Total BUSCO groups searched", "n_total"),
]

_ONE_LINE = re.compile(
    r"C:\s*([\d.]+)%\s*\[\s*S:\s*([\d.]+)%\s*,\s*D:\s*([\d.]+)%\s*\]\s*,\s*"
    r"F:\s*([\d.]+)%\s*,\s*M:\s*([\d.]+)%\s*,\s*n:\s*(\d+)"
    r"(?:\s*,\s*E:\s*([\d.]+)%)?")


def parse_short_summary_txt(path: Path) -> Dict[str, object]:
    """Parse a BUSCO ``short_summary*.txt`` (v3-v6)."""
    r: Dict[str, object] = {}
    txt = Path(path).read_text(encoding="utf-8", errors="replace")
    m = re.search(r"BUSCO version is:\s*(\S+)", txt)
    if m:
        r["busco_version"] = m.group(1)
    m = re.search(r"lineage dataset is:\s*([^\s(]+)", txt)
    if m:
        r["lineage"] = m.group(1)
    m = re.search(r"Creation date:\s*([^,)]+)", txt)
    if m:
        r["lineage_date"] = m.group(1).strip()
    m = re.search(r"number of (?:genomes|species):\s*(\d+)", txt)
    if m:
        r["lineage_n_species"] = int(m.group(1))
    m = re.search(r"number of BUSCOs:\s*(\d+)", txt)
    if m:
        r["lineage_n_buscos"] = int(m.group(1))
    m = re.search(r"BUSCO notation for file\s+(.+?)\s*$", txt, re.M)
    if m:
        r["input_file"] = m.group(1).strip()
    m = re.search(r"run in mode:\s*(\S+)", txt)
    if m:
        r["mode"] = m.group(1)
    m = re.search(r"Gene predictor used:\s*(\S+)", txt)
    if m:
        r["gene_predictor"] = m.group(1)
    m = _ONE_LINE.search(txt)
    if m:
        r["one_line"] = m.group(0)
        for key, val in zip(("pct_C", "pct_S", "pct_D", "pct_F", "pct_M"),
                            m.groups()[:5]):
            r[key] = float(val)
        r["n_total"] = int(m.group(6))
        if m.group(7) is not None:
            r["pct_E"] = float(m.group(7))
    for pat, key in _SUMMARY_COUNT_PATTERNS:
        mm = re.search(pat, txt, re.M | re.I)
        if mm and not (key == "n_total" and "n_total" in r):
            # n from the one-line summary is the denominator BUSCO uses; since
            # 5.3.0 "Total BUSCO groups searched" comes from dataset.cfg instead
            r[key] = int(mm.group(1))
    m = re.search(r"of which (\d+) contain internal stop codons", txt)
    if m:
        r["n_internal_stop"] = int(m.group(1))
    # "Assembly Statistics:" block (BUSCO >= 5.5): "<value>\t<label>"
    asm = {}
    blk = re.search(r"Assembly Statistics:\s*\n(.*?)(?:\n\s*\n|\Z)", txt, re.S)
    if blk:
        for ln in blk.group(1).splitlines():
            parts = [p.strip() for p in ln.strip().split("\t") if p.strip()]
            if len(parts) >= 2:
                asm[parts[1]] = parts[0]
    if asm:
        r["assembly"] = asm
    deps = re.search(r"Dependencies and versions:\s*\n(.*?)(?:\n\s*\n|\Z)", txt, re.S)
    if deps:
        d = {}
        for ln in deps.group(1).splitlines():
            if ":" in ln:
                k, v = ln.strip().split(":", 1)
                d[k.strip()] = v.strip()
        if d:
            r["dependencies"] = d
    return r


# --------------------------------------------------------------------------
# short_summary (json, BUSCO >= 5.1)
# --------------------------------------------------------------------------

def parse_short_summary_json(path: Path) -> Dict[str, object]:
    """Parse ``short_summary*.json``; keys normalised like the txt parser."""
    with open(path, encoding="utf-8", errors="replace") as fh:
        d = json.load(fh)
    r: Dict[str, object] = {}
    res = d.get("results", {}) or {}
    lin = d.get("lineage_dataset", {}) or {}
    par = d.get("parameters", {}) or {}
    ver = d.get("versions", {}) or {}
    one = res.get("one_line_summary")
    if one:
        r["one_line"] = one
        m = _ONE_LINE.search(one)
        if m and m.group(7) is not None:
            r["pct_E"] = float(m.group(7))
    mapping = {
        "Complete percentage": "pct_C", "Single copy percentage": "pct_S",
        "Multi copy percentage": "pct_D", "Fragmented percentage": "pct_F",
        "Missing percentage": "pct_M",
    }
    for k, v in mapping.items():
        if res.get(k) is not None:
            r[v] = _num(res[k])
    cmap = {
        "Complete BUSCOs": "n_complete", "Single copy BUSCOs": "n_single",
        "Multi copy BUSCOs": "n_duplicated", "Fragmented BUSCOs": "n_fragmented",
        "Missing BUSCOs": "n_missing", "n_markers": "n_total",
    }
    for k, v in cmap.items():
        if res.get(k) is not None:
            r[v] = _int(res[k])
    if res.get("internal_stop_codon_count") is not None:
        r["n_internal_stop"] = _int(res["internal_stop_codon_count"])
    if res.get("internal_stop_codon_percent") is not None and "pct_E" not in r:
        r["pct_E"] = _num(res["internal_stop_codon_percent"])
    if res.get("avg_identity") is not None:
        r["avg_identity"] = _num(res["avg_identity"])
    if res.get("domain"):
        r["domain"] = res["domain"]
    if lin.get("name"):
        r["lineage"] = lin["name"]
    if lin.get("creation_date"):
        r["lineage_date"] = lin["creation_date"]
    if lin.get("number_of_buscos"):
        r["lineage_n_buscos"] = _int(lin["number_of_buscos"])
    if lin.get("number_of_species"):
        r["lineage_n_species"] = _int(lin["number_of_species"])
    if ver.get("busco"):
        r["busco_version"] = str(ver["busco"])
    if par.get("mode"):
        r["mode"] = par["mode"]
    if par.get("gene_predictor"):
        r["gene_predictor"] = par["gene_predictor"]
    if par.get("in"):
        r["input_file"] = par["in"]
    if par.get("lineage_dataset"):
        r["lineage_path"] = par["lineage_dataset"]
    metrics = d.get("metrics") or {}
    if not metrics:
        keys = ("Number of scaffolds", "Number of contigs", "Total length",
                "Percent gaps", "Scaffold N50", "Contigs N50")
        metrics = {k: res[k] for k in keys if k in res}
    if metrics:
        r["assembly"] = {k: str(v) for k, v in metrics.items()}
    deps = {k: (".".join(map(str, v[:3])) if isinstance(v, list) else str(v))
            for k, v in ver.items() if v is not None and k != "busco"}
    if deps:
        r["dependencies"] = deps
    return r


# --------------------------------------------------------------------------
# compleasm summary.txt
# --------------------------------------------------------------------------

def parse_compleasm_summary(path: Path, lineage: str = "") -> Dict[str, object]:
    """Parse compleasm ``summary.txt`` (``S:93.33%, 237`` lines).

    compleasm appends one ``## lineage: X`` block per lineage; the block for
    ``lineage`` is used (the first block when not given or not found).
    """
    r: Dict[str, object] = {"tool": "compleasm"}
    txt = Path(path).read_text(encoding="utf-8", errors="replace")
    blocks = re.split(r"(?m)^(?=##\s*lineage:)", txt)
    blocks = [b for b in blocks if re.search(r"lineage:", b)] or [txt]
    txt = next((b for b in blocks if lineage and re.search(r"lineage:\s*%s\b" % re.escape(lineage), b)),
               blocks[0])
    m = re.search(r"lineage:\s*(\S+)", txt)
    if m:
        r["lineage"] = m.group(1)
    vals = {}
    for code in ("S", "D", "R", "F", "I", "M"):
        mm = re.search(r"^\s*%s:\s*([\d.]+)%%\s*,\s*(\d+)" % code, txt, re.M)
        if mm:
            vals[code] = (float(mm.group(1)), int(mm.group(2)))
    mm = re.search(r"^\s*N:\s*(\d+)", txt, re.M)
    if mm:
        r["n_total"] = int(mm.group(1))
    if vals:
        s = vals.get("S", (0.0, 0))
        dd = vals.get("D", (0.0, 0))
        if "R" in vals:   # --retrocopy: R reported apart from D; full_table maps it to D
            dd = (dd[0] + vals["R"][0], dd[1] + vals["R"][1])
        f = vals.get("F", (0.0, 0))
        i = vals.get("I", (0.0, 0))
        mi = vals.get("M", (0.0, 0))
        r.update(n_single=s[1], n_duplicated=dd[1], n_complete=s[1] + dd[1],
                 n_fragmented=f[1] + i[1], n_missing=mi[1],
                 n_interspersed=i[1],
                 pct_S=s[0], pct_D=dd[0], pct_C=round(s[0] + dd[0], 2),
                 pct_F=round(f[0] + i[0], 2), pct_M=mi[0])
    return r


# --------------------------------------------------------------------------
# full_table.tsv
# --------------------------------------------------------------------------

@dataclass
class Hit:
    seq: str = ""
    start: Optional[int] = None
    end: Optional[int] = None
    strand: str = ""
    score: Optional[float] = None
    length: Optional[int] = None

    def key(self) -> Tuple[str, int, int]:
        lo, hi = sorted((self.start or 0, self.end or 0))
        return (self.seq, lo, hi)


@dataclass
class GeneCall:
    status: str
    hits: List[Hit] = field(default_factory=list)


def _header_index(cols: List[str]) -> Dict[str, int]:
    idx: Dict[str, int] = {}
    for i, c in enumerate(cols):
        c = c.strip().lstrip("#").strip().lower()
        if c in ("busco id", "busco_group", "gene", "busco_id") and "id" not in idx:
            idx["id"] = i
        elif c == "status":
            idx["status"] = i
        elif c in ("sequence", "contig") and "seq" not in idx:
            idx["seq"] = i
        elif c in ("gene start", "start"):
            idx["start"] = i
        elif c in ("gene end", "end"):
            idx["end"] = i
        elif c == "strand":
            idx["strand"] = i
        elif c == "score":
            idx["score"] = i
        elif c == "length":
            idx["length"] = i
    return idx


def _guess_index(ncol: int) -> Dict[str, int]:
    if ncol >= 8:   # BUSCO v4-v6 genome: id status seq start end strand score length
        return dict(id=0, status=1, seq=2, start=3, end=4, strand=5, score=6, length=7)
    if ncol == 7:   # BUSCO v3 genome: id status contig start end score length
        return dict(id=0, status=1, seq=2, start=3, end=4, score=5, length=6)
    # transcriptome / proteins: id status seq score length
    return dict(id=0, status=1, seq=2, score=3, length=4)


def parse_full_table(path: Path) -> Dict[str, GeneCall]:
    """Parse BUSCO / compleasm full table -> {busco_id: GeneCall}.

    One BUSCO can have several rows (duplicated copies); the best status
    (S > D > F > M) wins and all hits are kept.
    """
    out: Dict[str, GeneCall] = {}
    header: Optional[Dict[str, int]] = None
    with open_text(path) as fh:
        for ln in fh:
            if not ln.strip():
                continue
            if ln.startswith("#") or ln.lower().startswith("gene\tstatus"):
                cols = ln.rstrip("\n").split("\t")
                if any(c.strip().lower() == "status" for c in cols):
                    header = _header_index(cols)
                continue
            parts = ln.rstrip("\n").split("\t")
            if len(parts) < 2:
                continue
            ix = header if header and "id" in header and "status" in header \
                else _guess_index(len(parts))
            gid = parts[ix["id"]].strip()
            word = parts[ix["status"]].strip().lower()
            code = STATUS_WORDS.get(word)
            if code is None:
                code = {"s": "S", "d": "D", "f": "F", "i": "F", "m": "M"}.get(word[:1], "M")
            call = out.get(gid)
            if call is None:
                call = out[gid] = GeneCall(code)
            elif STATUS_RANK[code] < STATUS_RANK[call.status]:
                call.status = code
            if code == "M":
                continue

            def g(key):
                i = ix.get(key)
                return parts[i].strip() if i is not None and i < len(parts) else ""

            call.hits.append(Hit(seq=g("seq"), start=_int(g("start")),
                                 end=_int(g("end")), strand=g("strand"),
                                 score=_num(g("score")), length=_int(g("length"))))
    # a BUSCO reported as Duplicated but with a single row (old versions)
    return out


def parse_missing_list(path: Path) -> List[str]:
    ids = []
    with open_text(path) as fh:
        for ln in fh:
            if ln.strip() and not ln.startswith("#"):
                ids.append(ln.split("\t")[0].strip())
    return ids


# --------------------------------------------------------------------------
# FASTA
# --------------------------------------------------------------------------

def read_fasta(path: Path) -> Iterator[Tuple[str, str]]:
    """Yield (header-without->, sequence) records; supports .gz."""
    hdr = None
    buf: List[str] = []
    with open_text(path) as fh:
        for ln in fh:
            ln = ln.rstrip("\r\n")
            if not ln:
                continue
            if ln.startswith(">"):
                if hdr is not None:
                    yield hdr, "".join(buf)
                hdr, buf = ln[1:].strip(), []
            elif hdr is not None:
                buf.append(ln.strip())
    if hdr is not None:
        yield hdr, "".join(buf)


_LOC = re.compile(r"^(.*):(\d+)-(\d+)$")


def header_location(hdr: str) -> Optional[Tuple[str, int, int]]:
    """Location key from a BUSCO sequence header.

    BUSCO v5/v6 headers look like
    ``1412at33208_400727_0:001b18|NODE_1:77-391|+ <desc>``; the second
    ``|`` field is ``seqid:start-end``.  Returns (seqid, lo, hi) or None.
    """
    first = hdr.split()[0]
    parts = first.split("|")
    cands = parts[1:] + [parts[0]]
    for p in cands:
        m = _LOC.match(p)
        if m:
            a, b = int(m.group(2)), int(m.group(3))
            return (m.group(1), min(a, b), max(a, b))
    return None


def hit_location(hit: Hit) -> Optional[Tuple[str, int, int]]:
    """Location key of a full_table hit (genome or transcriptome style)."""
    if hit.start is not None and hit.end is not None:
        return hit.key()
    m = _LOC.match(hit.seq or "")
    if m:
        a, b = int(m.group(2)), int(m.group(3))
        return (m.group(1), min(a, b), max(a, b))
    return None


def read_lineage_lengths(lineage_dir: Path) -> Dict[str, float]:
    """Expected protein length per BUSCO from a lineage dataset.

    Uses ``lengths_cutoff`` (odb9/odb10, 4th column = mean length) when
    present, otherwise the median length of ``refseq_db.faa(.gz)`` records
    per BUSCO (odb10/odb12), otherwise ``ancestral``.
    """
    lineage_dir = Path(lineage_dir)
    lc = lineage_dir / "lengths_cutoff"
    out: Dict[str, float] = {}
    if lc.is_file():
        with open_text(lc) as fh:
            for ln in fh:
                p = ln.split()
                if len(p) >= 4:
                    v = _num(p[3])
                    if v:
                        out[p[0]] = v
        if out:
            return out
    for name in ("refseq_db.faa", "refseq_db.faa.gz", "ancestral"):
        f = lineage_dir / name
        if not f.is_file():
            continue
        lens: Dict[str, List[int]] = {}
        for hdr, seq in read_fasta(f):
            gid = hdr.split()[0].split("_")[0]
            lens.setdefault(gid, []).append(len(seq.rstrip("*")))
        for gid, ls in lens.items():
            ls.sort()
            n = len(ls)
            out[gid] = float(ls[n // 2] if n % 2 else (ls[n // 2 - 1] + ls[n // 2]) / 2)
        if out:
            return out
    return out
