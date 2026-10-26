"""Discover BUSCO / compleasm result folders and load them into ``Sample``s.

Supported layouts
-----------------
* BUSCO v4-v6:   SAMPLE/run_<lineage>/{full_table.tsv, short_summary.*, busco_sequences/}
                 (auto-lineage: several run_* folders; batch mode: OUT/<file>/run_*)
* BUSCO v3:      run_SAMPLE/{full_table_SAMPLE.tsv, short_summary_SAMPLE.txt,
                 single_copy_busco_sequences/}
* a run folder given directly (SAMPLE/run_<lineage>)
* compleasm:     OUT/{summary.txt, <lineage>/full_table[_busco_format].tsv}
"""

from __future__ import annotations

import hashlib
import json
import os
import re
from concurrent.futures import ThreadPoolExecutor
from dataclasses import asdict, dataclass, field
from pathlib import Path
from typing import Dict, Iterable, List, Optional, Sequence, Tuple

from . import parsers
from .log import warn

SEQ_DIRS = {
    "S": "single_copy_busco_sequences",
    "D": "multi_copy_busco_sequences",
    "F": "fragmented_busco_sequences",
}
SKIP_DIRS = {"busco_downloads", "logs", "tmp", "hmmer_output", "busco_sequences",
             "miniprot_output", "metaeuk_output", "augustus_output", "blast_output",
             "prodigal_output", "translated_proteins", ".bbtools_output",
             "single_copy_busco_sequences", "multi_copy_busco_sequences",
             "fragmented_busco_sequences", "__pycache__", ".git"}
GENERIC_LINEAGES = {"eukaryota", "bacteria", "archaea", "viruses"}
FASTA_EXT = (".gz", ".fna", ".fasta", ".fa", ".fas", ".faa", ".fsa", ".ffn",
             ".pep", ".cds", ".scaffolds", ".contigs", ".txt")
CACHE_VERSION = 3


@dataclass
class Sample:
    name: str
    dir_name: str
    path: str                       # sample folder
    run_dir: str = ""               # folder holding full_table / sequences
    layout: str = ""                # busco_v4+, busco_v3, run_only, compleasm
    tool: str = "busco"
    other_runs: List[str] = field(default_factory=list)
    info: Dict[str, object] = field(default_factory=dict)   # parsed summary
    # gene -> [status, n_hits, best_length, best_score]
    calls: Dict[str, list] = field(default_factory=dict)
    # status -> gene -> "a"/"n"/"g" letters for faa / fna / gff present
    seqs: Dict[str, Dict[str, str]] = field(default_factory=dict)
    state: str = "OK"
    notes: List[str] = field(default_factory=list)

    # ---- convenience -------------------------------------------------
    @property
    def lineage(self) -> str:
        return str(self.info.get("lineage") or "unknown")

    def status(self, gid: str) -> str:
        c = self.calls.get(gid)
        return c[0] if c else "M"

    def counts(self) -> Dict[str, int]:
        n = {"S": 0, "D": 0, "F": 0, "M": 0}
        for c in self.calls.values():
            n[c[0]] += 1
        return n

    def to_json(self) -> dict:
        return asdict(self)

    @staticmethod
    def from_json(d: dict) -> "Sample":
        return Sample(**d)


# --------------------------------------------------------------------------
# naming
# --------------------------------------------------------------------------

_ACC = re.compile(r"(GC[AF]_\d+\.\d+|[SED]RR\d+|SAM[NED][A-Z]?\d+)")


def make_name(dir_name: str, style: str = "dir", replace_dots: bool = False,
              replace_char: str = "_") -> str:
    name = dir_name
    if name.startswith("run_") and style != "dir":
        name = name[4:]
    if style in ("no_ext", "stem", "accession"):
        changed = True
        while changed:
            changed = False
            for ext in FASTA_EXT:
                if name.lower().endswith(ext) and len(name) > len(ext):
                    name = name[: -len(ext)]
                    changed = True
    if style == "stem":
        name = name.split(".")[0]
    elif style == "accession":
        m = _ACC.search(name)
        if m:
            name = m.group(1)
    if replace_dots:
        name = name.replace(".", replace_char)
    name = re.sub(r"[\s/\\:;,()\[\]'\"|]+", replace_char, name).strip(replace_char)
    return name or dir_name


# --------------------------------------------------------------------------
# discovery
# --------------------------------------------------------------------------

def _has(d: Path, pattern: str) -> bool:
    return any(True for _ in d.glob(pattern))


def is_run_dir(d: Path) -> bool:
    """A folder that directly holds a BUSCO/compleasm full table or sequences."""
    try:
        return (_has(d, "full_table*.tsv") or _has(d, "full_table*.tsv.gz")
                or (d / "busco_sequences").is_dir()
                or (d / "single_copy_busco_sequences").is_dir())
    except OSError:
        return False


def _run_children(d: Path) -> List[Path]:
    out = []
    try:
        for c in sorted(d.iterdir()):
            if c.is_dir() and c.name not in SKIP_DIRS and is_run_dir(c):
                out.append(c)
    except OSError:
        pass
    return out


@dataclass
class Candidate:
    dir_name: str
    path: Path
    runs: List[Path]
    layout: str


def classify(d: Path) -> Optional[Candidate]:
    """Return a Candidate when ``d`` is a sample folder, else None."""
    if is_run_dir(d):
        if d.name.startswith("run_") and "_odb" in d.name:
            # a run folder given directly; sample = parent
            return Candidate(d.parent.name, d.parent, [d], "run_only")
        if (d / "single_copy_busco_sequences").is_dir() or _has(d, "full_table_*.tsv"):
            return Candidate(d.name, d, [d], "busco_v3")
        return Candidate(d.name, d, [d], "run_only")
    runs = _run_children(d)
    if runs:
        if (d / "summary.txt").is_file() and not any(r.name.startswith("run_") for r in runs):
            return Candidate(d.name, d, runs, "compleasm")
        # v4+ run folders are named run_<lineage>_odbNN; a folder whose run_*
        # children are BUSCO v3 sample folders (run_<sample>) is a container
        lineage_runs = [r for r in runs if r.name.startswith("run_") and "_odb" in r.name]
        if lineage_runs:
            return Candidate(d.name, d, lineage_runs, "busco_v4+")
        if all((r / "single_copy_busco_sequences").is_dir() or _has(r, "full_table_*.tsv")
               for r in runs):
            return None
        return Candidate(d.name, d, runs, "busco_v4+")
    return None


def discover(roots: Sequence[Path], max_depth: int = 4) -> List[Candidate]:
    found: List[Candidate] = []
    seen = set()

    def walk(d: Path, depth: int):
        if d.name in SKIP_DIRS:
            return
        c = classify(d)
        if c:
            key = str(c.path.resolve())
            if key not in seen:
                seen.add(key)
                found.append(c)
            return
        if depth >= max_depth:
            return
        try:
            kids = sorted(p for p in d.iterdir() if p.is_dir())
        except OSError as e:
            warn(f"cannot list {d}: {e}")
            return
        for k in kids:
            walk(k, depth + 1)

    for r in roots:
        r = Path(r)
        if not r.is_dir():
            warn(f"input is not a directory: {r}")
            continue
        walk(r.resolve(), 0)
    return found


# --------------------------------------------------------------------------
# run selection and loading
# --------------------------------------------------------------------------

def _lineage_of_run(run: Path) -> str:
    n = run.name[4:] if run.name.startswith("run_") else run.name
    return n


def lineage_match(want: str, have: str) -> bool:
    norm = lambda s: re.sub(r"_odb\d+$", "", s.strip().lower())
    w, h = want.strip().lower(), have.strip().lower()
    return w == h or (("_odb" not in w) and norm(w) == norm(h))


def choose_run(c: Candidate, lineage: str = "") -> Tuple[Optional[Path], List[Path]]:
    runs = c.runs
    if len(runs) == 1:
        if lineage and c.layout == "busco_v4+" and not lineage_match(lineage, _lineage_of_run(runs[0])):
            return None, runs
        return runs[0], []
    if lineage:
        hit = [r for r in runs if lineage_match(lineage, _lineage_of_run(r))]
        if hit:
            return hit[0], [r for r in runs if r is not hit[0]]
        return None, runs
    # auto-lineage: the lineage named by the top-level "specific" summary
    named = []
    for f in sorted(c.path.glob("short_summary.specific.*")):
        parts = f.name.split(".")
        if len(parts) > 2:
            named += [r for r in runs if _lineage_of_run(r) == parts[2] and r not in named]
    named.sort(key=lambda r: re.sub(r"_odb\d+$", "", _lineage_of_run(r)) in GENERIC_LINEAGES)
    if named:
        return named[0], [x for x in runs if x is not named[0]]
    specific = [r for r in runs
                if re.sub(r"_odb\d+$", "", _lineage_of_run(r)) not in GENERIC_LINEAGES]
    pick = (specific or runs)[0]
    return pick, [r for r in runs if r is not pick]


def _first(paths: Iterable[Path]) -> Optional[Path]:
    for p in paths:
        if p.is_file():
            return p
    return None


def _summary_files(c: Candidate, run: Path) -> Tuple[Optional[Path], Optional[Path]]:
    lin = _lineage_of_run(run)
    js = _first([run / "short_summary.json",
                 *sorted(c.path.glob(f"short_summary.specific.{lin}.*.json")),
                 *sorted(run.glob("short_summary*.json")),
                 *sorted(c.path.glob("short_summary*.json"))])
    tx = _first([run / "short_summary.txt",
                 *sorted(c.path.glob(f"short_summary.specific.{lin}.*.txt")),
                 *sorted(run.glob("short_summary*.txt")),
                 *sorted(c.path.glob("short_summary*.txt"))])
    return js, tx


def _full_table(run: Path) -> Optional[Path]:
    return _first([run / "full_table_busco_format.tsv", run / "full_table.tsv",
                   *sorted(run.glob("full_table*.tsv")),
                   *sorted(run.glob("full_table*.tsv.gz"))])


def _seq_base(run: Path) -> Optional[Path]:
    if (run / "busco_sequences").is_dir():
        return run / "busco_sequences"
    if (run / "single_copy_busco_sequences").is_dir():
        return run
    return None


def _inventory(run: Path) -> Dict[str, Dict[str, str]]:
    base = _seq_base(run)
    inv: Dict[str, Dict[str, str]] = {}
    if base is None:
        return inv
    for code, dname in SEQ_DIRS.items():
        d = base / dname
        if not d.is_dir():
            continue
        cur: Dict[str, str] = {}
        with os.scandir(d) as it:
            for e in it:
                if e.name.startswith("."):
                    continue
                stem, _, ext = e.name.rpartition(".")
                flag = {"faa": "a", "fna": "n", "ffn": "n", "gff": "g"}.get(ext)
                if flag:
                    cur[stem] = cur.get(stem, "") + flag
        inv[code] = cur
    return inv


def _signature(paths: Iterable[Optional[Path]]) -> str:
    h = hashlib.md5()
    for p in paths:
        if p is None:
            h.update(b"-")
            continue
        try:
            st = p.stat()
            h.update(f"{p}|{st.st_size}|{st.st_mtime_ns}".encode())
        except OSError:
            h.update(f"{p}|missing".encode())
    return h.hexdigest()


def load_sample(c: Candidate, name: str, lineage: str = "",
                cache_dir: Optional[Path] = None) -> Sample:
    run, others = choose_run(c, lineage)
    s = Sample(name=name, dir_name=c.dir_name, path=str(c.path), layout=c.layout,
               other_runs=[r.name for r in others])
    if run is None:
        s.state = "LINEAGE_ABSENT"
        s.notes.append(f"no run for lineage '{lineage}' (have: {', '.join(r.name for r in c.runs)})")
        return s
    s.run_dir = str(run)
    js, tx = _summary_files(c, run)
    ft = _full_table(run)
    base = _seq_base(run)
    sig_paths = [js, tx, ft] + ([base / d for d in SEQ_DIRS.values()] if base else [])
    sig = _signature(sig_paths)
    cache_file = None
    if cache_dir is not None:
        cache_file = cache_dir / (hashlib.md5(str(run).encode()).hexdigest() + ".json")
        if cache_file.is_file():
            try:
                d = json.loads(cache_file.read_text())
                if d.get("sig") == sig and d.get("v") == CACHE_VERSION:
                    cached = Sample.from_json(d["sample"])
                    cached.name, cached.dir_name = name, c.dir_name
                    return cached
            except Exception:
                pass

    if c.layout == "compleasm":
        s.tool = "compleasm"
        if (c.path / "summary.txt").is_file():
            s.info.update(parsers.parse_compleasm_summary(c.path / "summary.txt"))
        s.info.setdefault("lineage", run.name)
    try:
        if tx:
            s.info.update(parsers.parse_short_summary_txt(tx))
        if js:
            s.info.update({k: v for k, v in parsers.parse_short_summary_json(js).items()
                           if v is not None})
    except Exception as e:  # corrupt summary should not kill the batch
        s.notes.append(f"summary parse error: {e}")
    if "lineage" not in s.info:
        s.info["lineage"] = _lineage_of_run(run)
    s.info["summary_file"] = str(js or tx or "")
    s.info["full_table"] = str(ft or "")

    if ft:
        try:
            calls = parsers.parse_full_table(ft)
        except Exception as e:
            calls = {}
            s.notes.append(f"full_table parse error: {e}")
        for gid, gc in calls.items():
            best_len = max((h.length or 0 for h in gc.hits), default=0) or None
            best_sc = max((h.score or 0 for h in gc.hits), default=0) or None
            s.calls[gid] = [gc.status, len(gc.hits), best_len, best_sc]
    s.seqs = _inventory(run)
    if not ft:
        # infer statuses from sequence folders (missing BUSCOs filled later)
        for code in ("F", "D", "S"):
            for gid in s.seqs.get(code, {}):
                s.calls[gid] = [code, 1, None, None]
        miss = _first(sorted(run.glob("missing_busco_list*.tsv")))
        if miss:
            for gid in parsers.parse_missing_list(miss):
                s.calls.setdefault(gid, ["M", 0, None, None])
        s.notes.append("no full_table: statuses inferred from busco_sequences")
        s.state = "NO_FULL_TABLE" if s.calls else "FAILED"

    _qc(s)
    if cache_file is not None:
        try:
            cache_file.write_text(json.dumps({"v": CACHE_VERSION, "sig": sig,
                                              "sample": s.to_json()}))
        except OSError:
            pass
    return s


def _qc(s: Sample) -> None:
    if not s.calls and not s.info.get("n_total"):
        s.state = "FAILED"
        s.notes.append("no full_table and no summary")
        return
    n = s.counts()
    if s.calls and n["S"] + n["D"] + n["F"] == 0:
        s.state = "EMPTY"
        s.notes.append("no BUSCO found (all missing)")
    key = {"S": "n_single", "D": "n_duplicated", "F": "n_fragmented", "M": "n_missing"}
    if s.calls and s.state == "OK":
        diffs = [f"{k}:table={n[k]}/summary={s.info[v]}" for k, v in key.items()
                 if isinstance(s.info.get(v), int) and s.info[v] != n[k]]
        if diffs and s.tool == "busco":
            s.state = "MISMATCH"
            s.notes.append("full_table vs summary: " + "; ".join(diffs))
    # sequences declared but files absent (truncated copy, cleaned folders)
    if s.calls and s.seqs:
        for code in ("S", "D", "F"):
            want = sum(1 for c in s.calls.values() if c[0] == code)
            have = len(s.seqs.get(code, {}))
            if want and have < want:
                s.notes.append(f"{code}: {want - have} BUSCOs without sequence files")
    elif s.calls and not s.seqs:
        s.notes.append("no busco_sequences folder (summary/matrix only)")


# --------------------------------------------------------------------------
# rename map (compatible with busco_2_matrix v0.3 mapping files)
# --------------------------------------------------------------------------

_SYN = {
    "orig": ("original_name", "orig", "sample", "sample_name", "input_name", "directory",
             "dir_name", "basename", "old_name", "src", "dir"),
    "new": ("new_name", "new", "rename", "label", "alias", "taxon", "species",
            "output_name", "display_name", "final_name", "dst", "renamed", "name"),
    "path": ("abs_path_each_busco_results", "path", "abs_path", "busco_path", "busco_dir",
             "full_path", "abspath", "directory_path", "result_path", "result_dir",
             "busco_result_path", "busco_result"),
    "lineage": ("lineage_datasets", "lineage", "lid", "dataset", "odb", "lineage_db",
                "lineage_name", "db", "database", "odb_lineage"),
}


@dataclass
class MapRow:
    orig: str = ""
    new: str = ""
    path: str = ""
    lineage: str = ""


def read_map(path: Path) -> List[MapRow]:
    """Read a sample sheet.

    Accepted: a header line (with or without leading '#') naming columns by
    any synonym above, or a header-less 1-2 column file
    (``dir_name [new_name]``).  Tab, comma or semicolon separated.
    """
    rows: List[MapRow] = []
    lines = [l.rstrip("\n") for l in Path(path).read_text(encoding="utf-8").splitlines()]
    lines = [l for l in lines if l.strip() and not l.lstrip().startswith("##")]
    if not lines:
        return rows
    first = lines[0]
    delim = "\t" if "\t" in first else ("," if "," in first else (";" if ";" in first else None))
    split = (lambda l: [x.strip().strip("\"'") for x in l.split(delim)]) if delim \
        else (lambda l: l.split())
    head = [h.lstrip("#").strip().lower() for h in split(first)]
    cols: Dict[str, int] = {}
    for i, h in enumerate(head):
        for key, syns in _SYN.items():
            if h in syns and key not in cols:
                cols[key] = i
    body = lines[1:] if cols else lines
    if not cols:
        cols = {"orig": 0, "new": 1}
    for l in body:
        if l.lstrip().startswith("#"):
            continue
        p = split(l)
        g = lambda k: p[cols[k]] if k in cols and cols[k] < len(p) else ""
        rows.append(MapRow(g("orig"), g("new"), g("path"), g("lineage")))
    return rows


# --------------------------------------------------------------------------
# top-level loader
# --------------------------------------------------------------------------

def load_all(inputs: Sequence[str], sample_sheet: Optional[str] = None,
             lineage: str = "", name_style: str = "dir", replace_dots: bool = False,
             include: Optional[set] = None, exclude: Optional[set] = None,
             threads: int = 4, max_depth: int = 4,
             cache_dir: Optional[Path] = None) -> List[Sample]:
    cands: List[Candidate] = []
    rename: Dict[str, str] = {}
    lin_of: Dict[str, str] = {}
    if inputs:
        cands = discover([Path(i) for i in inputs], max_depth=max_depth)
    if sample_sheet:
        rows = read_map(Path(sample_sheet))
        by_name = {c.dir_name: c for c in cands}
        chosen: List[Candidate] = []
        for r in rows:
            c = None
            if r.path:
                p = Path(r.path)
                c = classify(p) if p.is_dir() else None
                if c is None and p.is_dir():
                    sub = discover([p], max_depth=2)
                    c = sub[0] if len(sub) == 1 else None
                if c is None:
                    warn(f"sample sheet: no BUSCO result in {r.path}")
                    continue
            elif r.orig:
                c = by_name.get(r.orig)
                if c is None and inputs:
                    for i in inputs:
                        p = Path(i) / r.orig
                        if p.is_dir():
                            c = classify(p)
                            if c:
                                break
                if c is None:
                    warn(f"sample sheet: '{r.orig}' not found under the input folder(s)")
                    continue
            else:
                continue
            chosen.append(c)
            if r.new:
                rename[str(c.path)] = r.new
            if r.lineage:
                lin_of[str(c.path)] = r.lineage
        cands = chosen
    if not cands:
        return []

    names: Dict[str, int] = {}
    jobs = []
    for c in cands:
        nm = rename.get(str(c.path)) or make_name(c.dir_name, name_style, replace_dots)
        if replace_dots and str(c.path) in rename:
            nm = nm.replace(".", "_")
        if include and nm not in include and c.dir_name not in include:
            continue
        if exclude and (nm in exclude or c.dir_name in exclude):
            continue
        if nm in names:
            names[nm] += 1
            new = f"{nm}_{names[nm]}"
            warn(f"duplicate sample name '{nm}' ({c.path}); renamed to '{new}'")
            nm = new
        else:
            names[nm] = 1
        jobs.append((c, nm, lin_of.get(str(c.path), lineage)))

    if cache_dir is not None:
        cache_dir.mkdir(parents=True, exist_ok=True)
    with ThreadPoolExecutor(max_workers=max(1, threads)) as ex:
        samples = list(ex.map(lambda j: load_sample(j[0], j[1], j[2], cache_dir), jobs))
    return samples
