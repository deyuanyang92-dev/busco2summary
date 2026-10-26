#!/usr/bin/env python3
"""Independent audit of one BUSCO run directory (does NOT import busco2summary).

Usage: check_run.py FULL_TABLE [SHORT_SUMMARY_TXT]
Prints one TSV line: path, version, mode/columns, n_rows, distinct ids,
ids with >1 status, S D F M from full_table, S D F M n from short_summary, agree.
"""
import re
import sys
from collections import defaultdict


def table(path):
    version, header, rows = "", None, []
    for ln in open(path, errors="replace"):
        ln = ln.rstrip("\n")
        if ln.startswith("#"):
            m = re.search(r"BUSCO version is:?\s*(\S+)", ln)
            if m:
                version = m.group(1)
            if ln.startswith("# Busco id") or ln.startswith("# BUSCO_id") or "Status" in ln:
                header = ln.lstrip("# ").split("\t")
            continue
        if ln.strip():
            rows.append(ln.split("\t"))
    st = defaultdict(set)
    for r in rows:
        if len(r) > 1:
            st[r[0]].add(r[1])
    return version, header, rows, st


def classify(statuses):
    # one status per id; report mixed separately
    for s, k in (("Complete", "S"), ("Single", "S"), ("Duplicated", "D"), ("Fragmented", "F"),
                 ("Interspersed", "F"), ("Missing", "M")):
        if s in statuses:
            return k
    return "?"


def summary(path):
    d = {}
    if not path:
        return d
    pat = {"S": r"Complete and single-copy BUSCOs", "D": r"Complete and duplicated BUSCOs",
           "F": r"Fragmented BUSCOs", "M": r"Missing BUSCOs", "n": r"Total BUSCO groups searched"}
    for ln in open(path, errors="replace"):
        for k, p in pat.items():
            if re.search(p, ln):
                m = re.match(r"\s*(\d+)\s", ln)
                if m:
                    d[k] = int(m.group(1))
    return d


def main():
    ft = sys.argv[1]
    ss = sys.argv[2] if len(sys.argv) > 2 else ""
    version, header, rows, st = table(ft)
    mixed = [i for i, s in st.items() if len(s) > 1]
    cnt = defaultdict(int)
    for s in st.values():
        cnt[classify(s)] += 1
    sm = summary(ss)
    agree = "NA"
    if sm:
        agree = "yes" if all(sm.get(k) == cnt.get(k, 0) for k in "SDFM") else "NO"
    ncols = sorted({len(r) for r in rows})
    print("\t".join(map(str, [ft, version, ",".join(map(str, ncols)), len(rows), len(st), len(mixed),
                              cnt["S"], cnt["D"], cnt["F"], cnt["M"], cnt["?"],
                              sm.get("S", ""), sm.get("D", ""), sm.get("F", ""), sm.get("M", ""),
                              sm.get("n", ""), agree, ";".join(mixed[:3])])))


if __name__ == "__main__":
    main()
