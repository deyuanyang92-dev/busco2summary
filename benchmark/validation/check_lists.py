#!/usr/bin/env python3
"""Independent occupancy-list check: recompute lists from raw full_tables of the samples
busco2summary put in the matrix, compare with its gene_lists/. No busco2summary import."""
import csv
import math
import os
import sys
from collections import defaultdict
out = sys.argv[1]
summ = {r["sample"]: r for r in csv.DictReader(open(f"{out}/busco_summary.tsv"), delimiter="\t")}
used = [r["sample"] for r in csv.DictReader(open(f"{out}/sample_stats.tsv"), delimiter="\t")]
status = {}
for s in used:
    st = {}
    rd = summ[s]["run_dir"]
    ft = [f for f in os.listdir(rd) if f.startswith("full_table") and f.endswith(".tsv")]
    ft = "full_table.tsv" if "full_table.tsv" in ft else sorted(ft)[0]
    for ln in open(os.path.join(rd, ft)):
        if ln.startswith("#") or not ln.strip():
            continue
        p = ln.rstrip("\n").split("\t")
        code = {"Complete": "S", "Duplicated": "D", "Fragmented": "F", "Missing": "M"}[p[1]]
        prev = st.get(p[0])
        if prev is None or "SDFM".index(code) < "SDFM".index(prev):
            st[p[0]] = code
    status[s] = st
N = len(used)
genes = set().union(*[set(v) for v in status.values()])
cnt = {k: defaultdict(int) for k in ("S", "SF", "SD", "SDF")}
for s in used:
    for g, c in status[s].items():
        for k in cnt:
            if c in k:
                cnt[k][g] += 1
bad = tot = 0
for k in cnt:
    for t in (50, 60, 70, 75, 80, 85, 90, 95, 100):
        need = math.ceil(t * N / 100 - 1e-9)
        mine = sorted(g for g in genes if cnt[k][g] >= need and cnt[k][g] > 0)
        theirs = sorted(open(f"{out}/gene_lists/{k}/{k}_occ{t:03d}.txt").read().split())
        tot += 1
        if mine != theirs:
            bad += 1
            print("DIFF", k, t, len(mine), len(theirs))
    print(k, "N=%d" % N, "occ50/90:", sum(1 for g in genes if cnt[k][g] >= math.ceil(.5*N-1e-9)),
          sum(1 for g in genes if cnt[k][g] >= math.ceil(.9*N-1e-9)))
print(f"lists compared {tot}, different {bad}")
