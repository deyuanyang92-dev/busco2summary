#!/usr/bin/env python3
"""Robustness benchmark on a synthetic folder with every awkward layout.

Cases (tests/fixtures.py): BUSCO v5 genome with .fna, BUSCO 6 miniprot (no .fna,
folder name ends in .fna), transcriptome 5-column full table, BUSCO v3 layout
(nested in a sub-folder), auto-lineage (generic + specific run), an empty
(all-missing) sample, a sample without full_table, compleasm output, and
non-BUSCO folders (busco_downloads/, notes/).

For each tool: exit status and, per case, whether the sample is reported with
the correct number of single-copy BUSCOs.
"""

from __future__ import annotations

import csv
import json
import os
import re
import subprocess
import sys
import tempfile
from pathlib import Path

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parent / "tests"))
import fixtures  # noqa: E402

CASES = {  # case -> (folder name as found on disk, expected S)
    "BUSCO v5 genome": "sampleA",
    "BUSCO 6 miniprot, dir *.fna": "sampleB.fna",
    "transcriptome table": "tranC",
    "BUSCO v3 (nested)": "run_v3D",
    "auto-lineage (2 runs)": "autoE",
    "empty sample": "emptyF",
    "no full_table": "noftG",
    "compleasm": "cmpH",
}


def run(cmd, cwd, env=None):
    p = subprocess.run(cmd, cwd=cwd, text=True, capture_output=True, env=env)
    return p.returncode, p.stdout, p.stderr


def main():
    tools = json.loads(Path(sys.argv[1]).read_text())
    out_tsv = Path(sys.argv[2])
    work = Path(tempfile.mkdtemp(prefix="b2s_robust_"))
    root = work / "busco"
    root.mkdir()
    exp = fixtures.build(root)
    fixtures.build_compleasm(root)
    expS = {"sampleA": exp["sampleA"]["S"], "sampleB.fna": exp["sampleB"]["S"],
            "tranC": exp["tranC"]["S"], "run_v3D": exp["v3D"]["S"], "autoE": exp["autoE"]["S"],
            "emptyF": 0, "noftG": exp["noftG"]["S"], "cmpH": 3}
    results = {}

    # busco2summary
    env = dict(os.environ, PYTHONPATH=tools["b2s_src"])
    rc, o, e = run([sys.executable, "-m", "busco2summary", "summary", "-i", str(root), "-o",
                    str(work / "b2s"), "--name-style", "dir", "-q"], work, env)
    got = {}
    if (work / "b2s" / "busco_summary.tsv").exists():
        for r in csv.DictReader(open(work / "b2s" / "busco_summary.tsv"), delimiter="\t"):
            got[r["dir_name"]] = int(r["S"]) if r["S"] not in ("", "NA") else None
    results["busco2summary"] = (rc, got)

    # busco_2_matrix v0.3
    rc, o, e = run([sys.executable, tools["v03"], "-i", str(root), "-o", str(work / "v03"),
                    "--no-extract", "--no-progress", "--overwrite"], work)
    got = {}
    f = work / "v03" / "species_summary.tsv"
    if f.exists():
        for r in csv.DictReader(open(f), delimiter="\t"):
            got[r["original_name"]] = int(r["n_single"])
    results["busco_2_matrix v0.3"] = (rc, got)

    # count_buscos.py
    rc, o, e = run([tools["bp_python"], tools["count_buscos"], "-i", str(root)], work)
    got = {}
    for ln in o.splitlines():
        m = re.match(r"^(\S+) (\d+)$", ln)
        if m:
            got[m.group(1)] = int(m.group(2))
    results["count_buscos.py"] = (rc, got)

    # busco2summary_smkv2 --only_summary (percent based)
    rc, o, e = run([sys.executable, tools["smkv2"], "-o", str(root), "--only_summary", "-ofb2s",
                    str(work / "smk"), "-y", "--no_smk"], work)
    got = {}
    f = work / "smk" / "summary_of_busco_res.txt"
    if f.exists():
        for r in csv.DictReader(open(f), delimiter="\t"):
            try:
                got[r["busco"]] = round(float(r["single"]) * int(r["number_of_BUSCOs"]) / 100)
            except ValueError:
                got[r["busco"]] = None
    results["busco2summary_smkv2 (summary)"] = (rc, got)

    # MultiQC
    rc, o, e = run([tools["multiqc"], "-m", "busco", "-f", "-o", str(work / "mqc"), str(root)], work)
    got = {}
    f = work / "mqc" / "multiqc_data" / "multiqc_busco.txt"
    if f.exists():
        for r in csv.DictReader(open(f), delimiter="\t"):
            name = re.sub(r"^short_summary(\.specific)?\.[^.]+_odb\d+\.", "", r["Sample"])
            got[name] = int(float(r["complete_single_copy"]))
    results["MultiQC busco"] = (rc, got)

    rows = []
    for tool, (rc, got) in results.items():
        row = {"tool": tool, "exit": rc}
        ok = 0
        for case, d in CASES.items():
            v = got.get(d, got.get(d.replace("run_", ""), got.get(d.split(".")[0])))
            if v is None:
                row[case] = "not reported"
            elif v == expS[d]:
                row[case] = "correct"
                ok += 1
            else:
                row[case] = f"wrong ({v} vs {expS[d]})"
        row["correct"] = f"{ok}/{len(CASES)}"
        rows.append(row)
        print(row, flush=True)
    with open(out_tsv, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(rows[0].keys()), delimiter="\t")
        w.writeheader()
        w.writerows(rows)


if __name__ == "__main__":
    main()
