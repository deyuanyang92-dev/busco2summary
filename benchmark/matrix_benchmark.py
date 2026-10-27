#!/usr/bin/env python3
"""Cross-validate gene-matrix construction (gene selection + per-gene FASTA) across tools.

Every tool turns a folder of BUSCO runs into per-gene FASTA files (one record per sample)
for the BUSCOs that pass an occupancy rule.  This script

  run    runs each tool on the same input copy and records exit code, time and memory
  score  reads every tool's per-gene FASTA files and scores them against a truth computed
         here, independently of every tool, from the raw full_table.tsv + busco_sequences

Scored per run:
  genes   gene set vs truth, both under busco2summary's rule (present in >= ceil(t*N)
          samples) and under the rule the tool itself documents (round / float / >= 3 taxa)
  cells   (gene, sample) records: present in truth with the right status? missing? extra?
  seqs    record sequence == that sample's BUSCO protein (terminal/internal '*' ignored)
  used    if a tool leaves samples out entirely, also scored against the truth for the samples
          it did use (separates "skips input" from "selects wrongly")
  inputs  md5 of every input file before and after (tools must not modify their input)

Alignment / trimming / tree programs are replaced by pass-through stubs written to
config['stub_dir'] at run time: only the selection and extraction stage is benchmarked.

usage:
  python3 matrix_benchmark.py run   --config bench.json --only b2s,bp,v03
  python3 matrix_benchmark.py score --config bench.json --out results/
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import math
import os
import re
import shutil
import signal
import subprocess
import sys
import time
from pathlib import Path
from typing import Dict, List, Optional, Set, Tuple

WORD = {"Complete": "S", "Duplicated": "D", "Fragmented": "F", "Missing": "M"}
SEQDIR = {"S": "single_copy_busco_sequences", "D": "multi_copy_busco_sequences",
          "F": "fragmented_busco_sequences"}

# --------------------------------------------------------------------- runs
# (tool, rule set, threshold) ; threshold "min3" = present in >= 3 samples
SETTINGS = {
    "b2s":   [("S", 0.5), ("S", 0.9), ("SF", 0.5), ("SF", 0.9), ("S", "min3"), ("SF", "min3")],
    "v03":   [("S", 0.5), ("S", 0.9), ("SF", 0.5), ("SF", 0.9)],
    "bp":    [("S", 0.5), ("S", 0.9)],
    "b2f":   [("S", 0.5), ("S", 0.9)],
    "bigey": [("S", 0.5), ("S", 0.9)],
    "zf":    [("SF", "min3")],          # main scripts/BUSCO_extraction.sh, BUSCO v5 branch
    "zf2":   [("S", "min3")],           # scripts/version2/BUSCO_extraction.sh
}
NAMES = {
    "b2s": "busco2summary 1.1.0", "v03": "busco_2_matrix v0.3 (legacy)",
    "bp": "BUSCO_phylogenomics (jamiemcg, bioconda 20240919)",
    "b2f": "busco2fasta.py (lstevens17/busco2phylo-nf)",
    "bigey": "busco-phylo.py (bigey/BUSCO-phylogenetics-pipeline)",
    "zf": "ZF BUSCO_extraction.sh (main, xtmtd/Phylogenomics)",
    "zf2": "ZF BUSCO_extraction.sh (version2)",
}


def tag(sset, t):
    return f"{sset}_{t}" if t == "min3" else f"{sset}_{int(round(t * 100)):03d}"


def timed(cmd, cwd: Path, env=None, stdin=None, timeout=14400) -> dict:
    cwd.mkdir(parents=True, exist_ok=True)
    t0 = time.time()
    # own process group, so a timeout also stops the tool's worker processes
    p = subprocess.Popen(["/usr/bin/time", "-v"] + cmd, cwd=cwd, env=env, text=True,
                         stdin=subprocess.PIPE, stdout=subprocess.PIPE, stderr=subprocess.PIPE,
                         start_new_session=True)
    try:
        out, err = p.communicate(stdin, timeout=timeout)
        rc = p.returncode
    except subprocess.TimeoutExpired:
        os.killpg(p.pid, signal.SIGKILL)
        out, err = p.communicate()
        rc, err = -9, err + "\nTIMEOUT"
    (cwd / "tool.log").write_text(out + "\n---- stderr ----\n" + err)
    m = re.search(r"Maximum resident set size \(kbytes\): (\d+)", err)
    return {"cmd": " ".join(map(str, cmd)), "rc": rc, "wall_s": round(time.time() - t0, 2),
            "max_rss_mb": round(int(m.group(1)) / 1024, 1) if m else None}


# Pass-through stand-ins for the alignment / trimming / tree programs: the "alignment" is the
# unaligned input, so downstream steps still find their files.
STUB = r'''#!/usr/bin/env python3
import os, shutil, sys
a, prog = sys.argv[1:], os.path.basename(sys.argv[0])
def opt(*names):
    for n in names:
        if n in a:
            return a[a.index(n) + 1]
if prog == "mafft":
    sys.stdout.write(open(a[-1]).read())
elif prog == "muscle" and opt("-align"):
    shutil.copy(opt("-align"), opt("-output"))
elif prog == "trimal" and opt("-in"):
    shutil.copy(opt("-in"), opt("-out"))
elif prog == "clipkit" and opt("--output", "-o"):
    shutil.copy(a[0], opt("--output", "-o"))
'''
STUB_PROGS = ["mafft", "muscle", "trimal", "clipkit", "iqtree", "iqtree2", "iqtree3", "fasttree",
              "FastTree", "raxml-ng"]


def ensure_stubs(d: Path):
    d.mkdir(parents=True, exist_ok=True)
    for name in STUB_PROGS:
        f = d / name
        f.write_text(STUB)
        f.chmod(0o755)


def run_one(tool, ds_name, inp: Path, sset, t, work: Path, cfg) -> Tuple[dict, Path]:
    """Run one setting; return (run info, folder holding the per-gene FASTA files)."""
    w = work / ds_name / f"{tool}_{tag(sset, t)}"
    if w.exists():
        shutil.rmtree(w)
    w.mkdir(parents=True)
    stub = cfg["stub_dir"]
    env = dict(os.environ)
    if tool == "b2s":
        out = w / "out"
        if t == "min3":
            r1 = timed([cfg["b2s"], "genelists", "-i", str(inp), "-o", str(w / "lists"), "--keep-empty",
                        "--no-cache", "--name-style", "dir", "-t", "50", "--min-taxa", "3", "--sets", sset,
                        "-q"], w)
            r = timed([cfg["b2s"], "extract", "-i", str(inp), "-o", str(out), "--keep-empty", "--no-cache",
                       "--name-style", "dir", "--set", sset, "--genes",
                       str(w / "lists" / "gene_lists" / sset / f"{sset}_min3.txt"), "--stop", "keep",
                       "-q"], w)
            r["wall_s"] = round(r["wall_s"] + r1["wall_s"], 2)
        else:
            r = timed([cfg["b2s"], "extract", "-i", str(inp), "-o", str(out), "--keep-empty", "--no-cache",
                       "--name-style", "dir", "--set", sset, "--occupancy", str(t), "--stop", "keep", "-q"], w)
        fa = next((p.parent for p in out.rglob("by_gene_faa/*.faa")), out)
    elif tool == "v03":
        r = timed([cfg["python"], cfg["v03"], "-i", str(inp), "-o", str(w / "out"), "--modes", sset,
                   "--seqtype", "protein", "-gpt", str(t), "--only-gps", "--keep-dots", "--no-progress",
                   "--overwrite"], w)
        fa = next((p.parent for p in (w / "out").rglob("*.faa")), w / "out")
    elif tool == "bp":
        env["PATH"] = f"{stub}:{cfg['bp_bin']}:{env['PATH']}"
        r = timed([cfg["bp_bin"] + "/python", cfg["bp_bin"] + "/BUSCO_phylogenomics.py", "-i", str(inp),
                   "-o", str(w / "out"), "-t", "8", "--supermatrix_only", "-psc", str(t * 100)], w, env)
        fa = w / "out" / "supermatrix" / "sequences"
    elif tool == "b2f":
        r = timed([cfg["bp_bin"] + "/python", cfg["b2f"], "-b", str(inp) + "/", "-o", str(w / "out"),
                   "-s", "protein", "-p", str(t)], w)
        fa = w / "out"
    elif tool == "bigey":
        adapt = w / "input"            # the script expects run_<species>/busco_sequences/
        adapt.mkdir()
        for s in sorted(inp.iterdir()):
            runs = sorted(s.glob("run_*_odb*"))
            if runs:
                (adapt / f"run_{s.name}").symlink_to(runs[0], target_is_directory=True)
        env["PATH"] = f"{stub}:{cfg['bp_bin']}:{env['PATH']}"
        r = timed([cfg["bp_bin"] + "/python", cfg["bigey"], "--busco-dir", str(adapt), "--output-dir",
                   str(w / "out"), "--supermatrix", "--stop-early", "--percent-single-copy", str(t),
                   "--threads", "8"], w, env)
        fa = w / "out" / "proteins"
    elif tool in ("zf", "zf2"):
        env["PATH"] = f"{cfg['zf_bin']}:{env['PATH']}"
        r = timed(["bash", cfg[tool]], w, env, stdin=f"8\n2\n{inp}\n")
        fa = w / "2-loci_filter" / "faa"
    else:
        raise SystemExit(f"unknown tool {tool}")
    r.update(tool=tool, dataset=ds_name, set=sset, threshold=t, fasta_dir=str(fa))
    (w / "run.json").write_text(json.dumps(r, indent=1))
    return r, fa


def md5_tree(root: Path) -> str:
    h = hashlib.md5()
    for p in sorted(root.rglob("*")):
        if p.is_file() and not p.is_symlink():
            h.update(str(p.relative_to(root)).encode())
            h.update(hashlib.md5(p.read_bytes()).digest())
    return h.hexdigest()


def cmd_run(args, cfg):
    work = Path(cfg["work"])
    ensure_stubs(Path(cfg["stub_dir"]))
    only = set(args.only.split(",")) if args.only else set(SETTINGS)
    for ds_name, inp in cfg["datasets"].items():
        inp = Path(inp)
        before = md5_tree(inp)
        for tool in SETTINGS:
            if tool not in only:
                continue
            for sset, t in SETTINGS[tool]:
                r, _ = run_one(tool, ds_name, inp, sset, t, work, cfg)
                print(f"{ds_name}\t{tool}\t{tag(sset, t)}\trc={r['rc']}\t{r['wall_s']}s", flush=True)
        after = md5_tree(inp)
        (work / ds_name / "input_md5.txt").write_text(f"{before}\n{after}\n")
        print(f"{ds_name}: input unchanged = {before == after}", flush=True)


# --------------------------------------------------------------------- truth

def read_fasta(p: Path) -> List[Tuple[str, str]]:
    recs, h, buf = [], None, []
    with open(p, errors="replace") as fh:
        for ln in fh:
            ln = ln.rstrip("\n")
            if ln.startswith(">"):
                if h is not None:
                    recs.append((h, "".join(buf)))
                h, buf = ln[1:], []
            elif h is not None:
                buf.append(ln.strip())
    if h is not None:
        recs.append((h, "".join(buf)))
    return recs


def norm(seq: str) -> str:
    return re.sub(r"[\s*]", "", seq).upper()


def load_truth(inp: Path):
    status: Dict[str, Dict[str, str]] = {}
    seqs: Dict[Tuple[str, str], Set[str]] = {}
    for s in sorted(inp.iterdir()):
        runs = sorted(s.glob("run_*_odb*"))
        if not runs:
            continue
        run = runs[0]
        st: Dict[str, str] = {}
        for ln in open(run / "full_table.tsv"):
            if ln.startswith("#") or not ln.strip():
                continue
            p = ln.rstrip("\n").split("\t")
            c = WORD[p[1]]
            st[p[0]] = c if p[0] not in st else min(st[p[0]], c, key="SDFM".index)
        status[s.name] = st
        for d in SEQDIR.values():
            for f in (run / "busco_sequences" / d).glob("*.faa"):
                seqs[(s.name, f.stem)] = {norm(x) for _, x in read_fasta(f)}
    genes = sorted(set().union(*[set(v) for v in status.values()]))
    return status, genes, seqs


def truth_genes(status, genes, sset, t, rule="ceil") -> Set[str]:
    n = len(status)
    cnt = {g: sum(status[s].get(g, "M") in sset for s in status) for g in genes}
    if t == "min3":
        return {g for g, c in cnt.items() if c >= 3}
    if rule == "ceil":
        need = math.ceil(t * n - 1e-9)
        return {g for g, c in cnt.items() if c >= need and c > 0}
    if rule == "round":                      # busco2fasta: int(round(n * p))
        need = int(round(n * t, 0))
        return {g for g, c in cnt.items() if c >= need and c > 0}
    if rule == "float":                      # count / n >= t
        return {g for g, c in cnt.items() if c / n >= t and c > 0}
    raise ValueError(rule)


DECLARED = {"b2s": "ceil", "v03": "float", "bp": "float", "b2f": "round", "bigey": "float",
            "zf": "min3", "zf2": "min3"}


# --------------------------------------------------------------------- score

def resolve(header: str, samples: List[str]) -> Optional[str]:
    h = header.split()[0] if header.split() else header
    if h in samples:
        return h
    best = None
    for s in samples:      # "sample.gene", "sample|..." etc.: longest sample-name prefix
        if h.startswith(s) and (len(h) == len(s) or not h[len(s)].isalnum()):
            if best is None or len(s) > len(best):
                best = s
    return best


def score_run(r: dict, status, genes, seqs) -> dict:
    samples = sorted(status)
    sset, t = r["set"], r["threshold"]
    fa = Path(r["fasta_dir"])
    out: Dict[str, Dict[str, str]] = {}
    unresolved = dup = 0
    files = sorted(fa.glob("*.faa")) + sorted(fa.glob("*.fasta")) + sorted(fa.glob("*.fa")) \
        if fa.is_dir() else []
    for f in files:
        g = f.name.split(".")[0]
        recs = read_fasta(f)
        if not recs:
            continue
        d = out.setdefault(g, {})
        for h, s in recs:
            smp = resolve(h, samples)
            if smp is None:
                unresolved += 1
                continue
            if smp in d:
                dup += 1
            d[smp] = s
    got = set(out)
    want_ceil = truth_genes(status, genes, sset, t, "ceil")
    rule = DECLARED[r["tool"]]
    want_decl = want_ceil if rule in ("ceil", "min3") else truth_genes(status, genes, sset, t, rule)
    cells = cell_ok = cell_wrong_status = cell_missing = seq_ok = seq_bad = 0
    for g, d in out.items():
        expect = {s for s in samples if status[s].get(g, "M") in sset}
        for s, seq in d.items():
            cells += 1
            if s in expect:
                cell_ok += 1
            else:
                cell_wrong_status += 1
            if norm(seq) in seqs.get((s, g), set()):
                seq_ok += 1
            else:
                seq_bad += 1
        cell_missing += len(expect - set(d))
    # Samples the tool used at all. If it silently skipped samples (ZF: runs without .fna),
    # re-score against the truth for only those samples, under the tool's own rule.
    used = {s for d in out.values() for s in d}
    sub = {s: status[s] for s in used}
    if used and cell_missing:
        want_sub = truth_genes(sub, genes, sset, t, "ceil" if rule == "min3" else rule)
        miss_sub = sum(len({s for s in sub if sub[s].get(g, "M") in sset} - set(out.get(g, {})))
                       for g in want_sub)
        vs_sub = "exact" if got == want_sub else f"+{len(got - want_sub)}/-{len(want_sub - got)}"
    else:
        vs_sub, miss_sub = "-", "-"
    return {
        "dataset": r["dataset"], "tool": NAMES[r["tool"]], "set": sset,
        "threshold": t if t == "min3" else f"{int(round(t * 100))}%", "exit": r["rc"],
        "wall_s": r["wall_s"], "max_rss_mb": r["max_rss_mb"],
        "genes_out": len(got), "genes_truth": len(want_ceil),
        "vs_truth": "exact" if got == want_ceil else f"+{len(got - want_ceil)}/-{len(want_ceil - got)}",
        "own_rule": rule, "genes_own_rule": len(want_decl),
        "vs_own_rule": "exact" if got == want_decl else f"+{len(got - want_decl)}/-{len(want_decl - got)}",
        "records": cells, "records_right_status": cell_ok, "records_wrong_status": cell_wrong_status,
        "records_missing": cell_missing, "seq_identical": seq_ok, "seq_different": seq_bad,
        "headers_unresolved": unresolved, "headers_duplicated": dup,
        "samples_total": len(samples), "samples_in_output": len(used),
        "vs_truth_used_samples": vs_sub, "records_missing_used_samples": miss_sub,
    }


def cmd_score(args, cfg):
    work = Path(cfg["work"])
    outd = Path(args.out)
    outd.mkdir(parents=True, exist_ok=True)
    rows = []
    for ds_name, inp in cfg["datasets"].items():
        status, genes, seqs = load_truth(Path(inp))
        print(f"{ds_name}: truth from {len(status)} samples, {len(genes)} BUSCOs", flush=True)
        for rj in sorted((work / ds_name).glob("*/run.json")):
            rows.append(score_run(json.loads(rj.read_text()), status, genes, seqs))
        md5 = work / ds_name / "input_md5.txt"
        if md5.exists():
            a, b = md5.read_text().split()[:2]
            print(f"{ds_name}: input files unchanged by all tools: {a == b}")
    order = list(NAMES.values())
    rows.sort(key=lambda r: (r["dataset"], r["set"], str(r["threshold"]), order.index(r["tool"])))
    with open(outd / "matrix_benchmark.tsv", "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(rows[0]), delimiter="\t")
        w.writeheader()
        w.writerows(rows)
    for r in rows:
        print("\t".join(str(r[k]) for k in ("dataset", "tool", "set", "threshold", "exit", "genes_out",
                                            "vs_truth", "vs_own_rule", "records", "records_wrong_status",
                                            "records_missing", "seq_identical", "seq_different",
                                            "headers_unresolved", "samples_in_output",
                                            "vs_truth_used_samples", "records_missing_used_samples",
                                            "wall_s")))


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("action", choices=["run", "score"])
    ap.add_argument("--config", required=True)
    ap.add_argument("--only", default="", help="comma list of tools for 'run' (default all)")
    ap.add_argument("--out", default="results")
    a = ap.parse_args()
    cfg = json.loads(Path(a.config).read_text())
    (cmd_run if a.action == "run" else cmd_score)(a, cfg)


if __name__ == "__main__":
    sys.exit(main())
