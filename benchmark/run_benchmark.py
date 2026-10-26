#!/usr/bin/env python3
"""Benchmark busco2summary against other BUSCO summary / matrix tools.

Ground truth is computed here, independently of every tool, from the raw
``full_table.tsv`` files (status per BUSCO per sample).  Each tool runs under
``/usr/bin/time -v``; its output is parsed and scored against the truth:

* per-sample S / D / F / M counts (summary tools)
* per-BUSCO single-copy occupancy lists at 50 % and 90 % and >= 3 taxa
  (matrix / extraction tools)

usage:
  python3 run_benchmark.py --dataset D1=/path/busco_dir --work /tmp/bench \
      --tools-dir /path/to/cloned/tools --out results/
Third-party tools are not redistributed; clone them first (see README.md).
"""

from __future__ import annotations

import argparse
import csv
import glob
import json
import math
import os
import re
import shutil
import subprocess
import sys
import time
from pathlib import Path
from typing import Dict, List, Optional, Tuple

WORD = {"Complete": "S", "Single": "S", "Duplicated": "D", "Fragmented": "F", "Missing": "M"}


# ------------------------------------------------------------------ truth

def truth(clean_dir: Path, lineage: str = "*") -> Tuple[Dict[str, Dict[str, str]], List[str]]:
    st: Dict[str, Dict[str, str]] = {}
    for sd in sorted(clean_dir.iterdir()):
        fts = sorted(glob.glob(str(sd / f"run_{lineage}" / "full_table.tsv")))
        if not fts:
            continue
        d: Dict[str, str] = {}
        for ln in open(fts[0]):
            if ln.startswith("#") or not ln.strip():
                continue
            p = ln.rstrip("\n").split("\t")
            c = WORD.get(p[1], "M")
            d[p[0]] = c if p[0] not in d else min(d[p[0]], c, key="SDFM".index)
        st[sd.name] = d
    genes = sorted(set().union(*[d.keys() for d in st.values()]))
    return st, genes


def truth_lists(st, genes, sset="S"):
    n = len(st)
    cnt = {g: sum(st[s].get(g, "M") in sset for s in st) for g in genes}
    out = {}
    for t in (0.5, 0.9):
        need = math.ceil(t * n - 1e-9)
        out[f"occ{int(t * 100)}"] = {g for g, c in cnt.items() if c >= need and c > 0}
    out["min3"] = {g for g, c in cnt.items() if c >= 3}
    return out, cnt


# ------------------------------------------------------------------ runner

def timed(cmd, cwd: Path, log: Path, stdin: Optional[str] = None, env=None, timeout=7200):
    cwd.mkdir(parents=True, exist_ok=True)
    full = ["/usr/bin/time", "-v"] + cmd
    t0 = time.time()
    with open(log, "w") as lf:
        try:
            p = subprocess.run(full, cwd=cwd, input=stdin, text=True, stdout=subprocess.PIPE,
                               stderr=subprocess.PIPE, env=env, timeout=timeout)
            rc, out, err = p.returncode, p.stdout, p.stderr
        except subprocess.TimeoutExpired:
            rc, out, err = -9, "", "TIMEOUT"
        lf.write(out)
        lf.write("\n---- stderr ----\n")
        lf.write(err)
    wall = time.time() - t0
    m = re.search(r"Maximum resident set size \(kbytes\): (\d+)", err)
    rss = int(m.group(1)) / 1024 if m else None
    return {"rc": rc, "wall_s": round(wall, 2), "max_rss_mb": round(rss, 1) if rss else None,
            "stdout": out, "stderr": err}


def score_counts(pred: Dict[str, Dict[str, int]], st, keys="SDFM") -> dict:
    ok = bad = miss = 0
    for s, d in st.items():
        want = {k: sum(1 for v in d.values() if v == k) for k in keys}
        got = pred.get(s)
        if got is None:
            miss += 1
            continue
        if all(got.get(k) == want[k] for k in keys if got.get(k) is not None) and \
                any(got.get(k) is not None for k in keys):
            ok += 1
        else:
            bad += 1
    return {"samples_correct": ok, "samples_wrong": bad, "samples_missing": miss}


def score_list(pred: Optional[set], want: set) -> str:
    if pred is None:
        return "n/a"
    fp, fn = len(pred - want), len(want - pred)
    return f"{len(pred)} ({'exact' if not fp and not fn else f'+{fp}/-{fn}'})"


# ------------------------------------------------------------------ tools

def run_b2s(ds, raw, clean, work, tools, py):
    r = timed([py, "-m", "busco2summary", "run", "-i", str(raw), "-o", str(work / "out"),
               "--keep-empty", "--no-cache", "-t", "50", "90", "--min-taxa", "3", "-j", "8", "-q"],
              work, work / "log.txt", env=dict(os.environ, PYTHONPATH=str(tools["b2s_src"])))
    out = work / "out"
    pred, lists = {}, {}
    if (out / "busco_summary.tsv").exists():
        for row in csv.DictReader(open(out / "busco_summary.tsv"), delimiter="\t"):
            pred[row["dir_name"]] = {k: int(row[k]) for k in "SDFM" if row[k] not in ("", "NA")}
        for k, f in (("occ50", "S_occ050"), ("occ90", "S_occ090"), ("min3", "S_min3")):
            p = out / "gene_lists" / "S" / f"{f}.txt"
            lists[k] = set(p.read_text().split()) if p.exists() else None
        for k, f in (("SF_occ90", "SF_occ090"),):
            p = out / "gene_lists" / "SF" / f"{f}.txt"
            lists[k] = set(p.read_text().split()) if p.exists() else None
    return r, pred, lists


def run_b2s_cached(ds, raw, clean, work, tools, py):
    """second run re-using the parse cache (resume)."""
    args = [py, "-m", "busco2summary", "run", "-i", str(raw), "-o", str(work / "out"),
            "--keep-empty", "-t", "50", "90", "--min-taxa", "3", "-j", "8", "-q", "--no-plots", "--no-html"]
    env = dict(os.environ, PYTHONPATH=str(tools["b2s_src"]))
    timed(args, work, work / "log_first.txt", env=env)
    r = timed(args, work, work / "log.txt", env=env)
    return r, {}, {}


def run_v03(ds, raw, clean, work, tools, py):
    r = timed([py, str(tools["v03"]), "-i", str(clean), "-o", str(work / "out"), "--no-extract",
               "--overwrite", "--no-progress", "-gpt", "0.5", "0.9"], work, work / "log.txt")
    out = work / "out"
    pred, lists = {}, {}
    f = out / "species_summary.tsv"
    if f.exists():
        for row in csv.DictReader(open(f), delimiter="\t"):
            pred[row["original_name"]] = {"S": int(row["n_single"]), "D": int(row["n_duplicated"]),
                                          "F": int(row["n_fragmented"]), "M": int(row["n_missing"])}
    f = out / "busco_binary_S.tsv"
    if f.exists():
        rows = list(csv.reader(open(f), delimiter="\t"))
        n = len(rows[0]) - 2
        cnt = {x[0]: sum(int(v) for v in x[1:-1]) for x in rows[1:]}
        lists = {"occ50": {g for g, c in cnt.items() if c / n >= 0.5},
                 "occ90": {g for g, c in cnt.items() if c / n >= 0.9},
                 "min3": {g for g, c in cnt.items() if c >= 3}}
        f2 = out / "busco_binary_SF.tsv"
        if f2.exists():
            rows = list(csv.reader(open(f2), delimiter="\t"))
            cnt = {x[0]: sum(int(v) for v in x[1:-1]) for x in rows[1:]}
            lists["SF_occ90"] = {g for g, c in cnt.items() if c / n >= 0.9}
    return r, pred, lists


def run_count_buscos(ds, raw, clean, work, tools, py):
    r = timed([tools["bp_python"], tools["count_buscos"], "-i", str(clean)], work, work / "log.txt")
    pred, lists = {}, {}
    rows = [l.split("\t") for l in r["stdout"].splitlines() if re.match(r"^\S+at\d+\t[YN]", l)]
    if rows:
        cnt = {x[0]: sum(v.strip() == "Y" for v in x[1:]) for x in rows}
        n = len(rows[0]) - 1
        lists = {"occ50": {g for g, c in cnt.items() if c >= math.ceil(0.5 * n - 1e-9) and c},
                 "occ90": {g for g, c in cnt.items() if c >= math.ceil(0.9 * n - 1e-9) and c},
                 "min3": {g for g, c in cnt.items() if c >= 3}}
    for l in r["stdout"].splitlines():
        m = re.match(r"^(\S+) (\d+)$", l)
        if m:
            pred[m.group(1)] = {"S": int(m.group(2))}
    return r, pred, lists


def run_busco2fasta(ds, raw, clean, work, tools, py):
    lists = {}
    res = None
    for p, key in ((0.5, "occ50"), (0.9, "occ90")):
        o = work / f"out_{key}"
        res = timed([py, tools["busco2fasta"], "-b", str(clean) + "/", "-o", str(o), "-p", str(p)],
                    work, work / f"log_{key}.txt")
        if o.exists():
            lists[key] = {f[:-4] for f in os.listdir(o) if f.endswith(".faa")}
    return res, {}, lists


def run_zf(ds, raw, clean, work, tools, py):
    script = tools["zf_script"]
    env = dict(os.environ, PATH=f"{tools['zf_env_bin']}:{os.environ['PATH']}")
    r = timed(["bash", script], work, work / "log.txt", stdin=f"8\n2\n{clean}\n", env=env)
    lists = {}
    d = work / "2-loci_filter" / "faa"
    if d.exists():
        lists["min3"] = {f[:-4] for f in os.listdir(d) if f.endswith(".faa")}
    return r, {}, lists


def run_smkv2(ds, raw, clean, work, tools, py):
    r = timed([py, tools["smkv2"], "-o", str(raw), "--only_summary", "-ofb2s", str(work / "out"),
               "-y", "--no_smk"], work, work / "log.txt")
    pred = {}
    f = work / "out" / "summary_of_busco_res.txt"
    if f.exists():
        for row in csv.DictReader(open(f), delimiter="\t"):
            try:
                n = int(row["number_of_BUSCOs"])
                pred[row["busco"]] = {k: int(round(float(row[c]) * n / 100)) for k, c in
                                      (("S", "single"), ("D", "duplicate"), ("F", "fragmented"),
                                       ("M", "missing"))}
            except (ValueError, KeyError):
                pass
    return r, pred, {}


def run_bigey(ds, raw, clean, work, tools, py):
    r = timed([tools["bigey_python"], tools["bigey"], "--directory", str(raw)], work, work / "log.txt")
    pred = {}
    for l in r["stdout"].splitlines()[1:]:
        p = l.split("\t")
        if len(p) >= 6:
            pred[p[0]] = {"pct_only": True}
    return r, pred, {}


def run_multiqc(ds, raw, clean, work, tools, py):
    r = timed([tools["multiqc"], "-m", "busco", "-f", "-o", str(work / "out"), str(raw)],
              work, work / "log.txt")
    pred = {}
    f = work / "out" / "multiqc_data" / "multiqc_busco.txt"
    if f.exists():
        for row in csv.DictReader(open(f), delimiter="\t"):
            name = re.sub(r"^short_summary(\.specific)?\.[^.]+_odb\d+\.", "", row["Sample"])
            try:
                pred[name] = {"S": int(float(row["complete_single_copy"])),
                                       "D": int(float(row["complete_duplicated"])),
                                       "F": int(float(row["fragmented"])),
                                       "M": int(float(row["missing"]))}
            except (KeyError, ValueError):
                pass
    return r, pred, {}


TOOLS = {
    "busco2summary": run_b2s,
    "busco2summary (resume)": run_b2s_cached,
    "busco_2_matrix v0.3": run_v03,
    "count_buscos.py": run_count_buscos,
    "busco2fasta.py": run_busco2fasta,
    "ZF BUSCO_extraction.sh": run_zf,
    "busco2summary_smkv2 (summary)": run_smkv2,
    "extract-busco-scores.py": run_bigey,
    "MultiQC busco": run_multiqc,
}


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--dataset", action="append", required=True, help="NAME=RAW_DIR")
    ap.add_argument("--work", required=True)
    ap.add_argument("--out", required=True)
    ap.add_argument("--tools-json", required=True, help="json with tool paths")
    ap.add_argument("--only", nargs="*")
    args = ap.parse_args()
    tools = json.loads(Path(args.tools_json).read_text())
    py = sys.executable
    out = Path(args.out)
    out.mkdir(parents=True, exist_ok=True)
    rows = []
    for spec in args.dataset:
        name, raw = spec.split("=", 1)
        raw, _, lineage = raw.partition("@")
        lineage = lineage or "*"
        raw = Path(raw)
        work = Path(args.work) / name
        clean = work / "clean_input"
        if clean.exists():
            shutil.rmtree(clean)
        clean.mkdir(parents=True)
        for d in sorted(raw.iterdir()):
            if d.is_dir() and glob.glob(str(d / f"run_{lineage}" / "full_table.tsv")):
                (clean / d.name).symlink_to(d.resolve())
        st, genes = truth(clean, lineage)
        tl, _ = truth_lists(st, genes, "S")
        tlf, _ = truth_lists(st, genes, "SF")
        print(f"[{name}] {len(st)} samples, {len(genes)} BUSCOs; truth S occ50={len(tl['occ50'])} "
              f"occ90={len(tl['occ90'])} min3={len(tl['min3'])} SF occ90={len(tlf['occ90'])}", flush=True)
        for tool, fn in TOOLS.items():
            if args.only and tool not in args.only:
                continue
            w = work / re.sub(r"\W+", "_", tool)
            if w.exists():
                shutil.rmtree(w)
            try:
                r, pred, lists = fn(name, clean, clean, w, tools, py)
            except Exception as e:  # tool missing etc.
                r, pred, lists = {"rc": "error", "wall_s": None, "max_rss_mb": None, "stderr": str(e)}, {}, {}
            sc = score_counts({k: v for k, v in pred.items() if "pct_only" not in v}, st) if pred else {}
            row = {"dataset": name, "n_samples": len(st), "tool": tool, "exit": r["rc"],
                   "wall_s": r["wall_s"], "max_rss_mb": r["max_rss_mb"],
                   "per_sample_counts": ("percent only" if pred and all("pct_only" in v for v in pred.values())
                                         else f"{sc.get('samples_correct', 0)}/{len(st)} correct"
                                         if pred else "n/a"),
                   "S_occ50": score_list(lists.get("occ50"), tl["occ50"]),
                   "S_occ90": score_list(lists.get("occ90"), tl["occ90"]),
                   "S_min3": score_list(lists.get("min3"), tl["min3"]),
                   "SF_occ90": score_list(lists.get("SF_occ90"), tlf["occ90"])}
            rows.append(row)
            print("  " + "  ".join(f"{k}={v}" for k, v in row.items() if k not in ("dataset",)), flush=True)
            if r["rc"] not in (0, "0"):
                (out / f"{name}_{re.sub(r'W+', '_', tool)}.stderr.txt").write_text(str(r.get("stderr"))[-5000:])
    with open(out / "benchmark_results.tsv", "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(rows[0].keys()), delimiter="\t")
        w.writeheader()
        w.writerows(rows)
    print(f"results -> {out / 'benchmark_results.tsv'}")


if __name__ == "__main__":
    main()
