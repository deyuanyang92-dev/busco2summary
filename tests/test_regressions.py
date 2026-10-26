"""Regression tests for the logic issues found in the 2026-10-27 review (docs/ISSUES.md #20-#30)."""

import csv
import time

from busco2summary import analysis
from busco2summary.cli import main

import fixtures
from fixtures import GENES, _write_run


def read_tsv(p):
    with open(p) as fh:
        return list(csv.DictReader(fh, delimiter="\t"))


def by(rows, key):
    return {r[key]: r for r in rows}


def run(*argv):
    assert main([str(a) for a in argv] + ["-q"]) == 0


def test_name_collision_after_rename(tmp_path):
    _write_run(tmp_path / "r1/A/run_metazoa_odb12", dict(zip(GENES, "SSSSSS")))
    _write_run(tmp_path / "r2/A/run_metazoa_odb12", dict(zip(GENES, "MMMMMS")))
    _write_run(tmp_path / "r2/A_2/run_metazoa_odb12", dict(zip(GENES, "FFFFFS")))
    run("summary", "-i", tmp_path / "r1", tmp_path / "r2", "-o", tmp_path / "o", "--no-cache")
    names = [r["sample"] for r in read_tsv(tmp_path / "o/busco_summary.tsv")]
    assert len(names) == 3 and len(set(names)) == 3


def test_sample_sheet_duplicate_path_used_once(tmp_path):
    _write_run(tmp_path / "in/X/run_metazoa_odb12", dict(zip(GENES, "SSSSSS")))
    sheet = tmp_path / "sheet.tsv"
    sheet.write_text(f"path\tnew_name\n{tmp_path}/in/X\tsp1\n{tmp_path}/in/X\tsp2\n")
    run("summary", "-s", sheet, "-o", tmp_path / "o", "--no-cache")
    assert len(read_tsv(tmp_path / "o/busco_summary.tsv")) == 1


def test_no_full_table_keeps_summary_denominator(tmp_path):
    _write_run(tmp_path / "in/G/run_metazoa_odb12", dict(zip(GENES, "SSMMMM")), full_table=False)
    run("summary", "-i", tmp_path / "in", "-o", tmp_path / "o", "--no-cache")
    r = read_tsv(tmp_path / "o/busco_summary.tsv")[0]
    assert (r["n_total"], r["S"], r["M"]) == ("6", "2", "4")
    assert abs(float(r["pct_C"]) - 100 * 2 / 6) < 0.01


def test_cache_sees_missing_list(tmp_path):
    run_dir = tmp_path / "in/G/run_metazoa_odb12"
    _write_run(run_dir, dict(zip(GENES, "SSMMMM")), full_table=False, summary=False)
    for p in run_dir.parent.glob("short_summary*"):
        p.unlink()
    run("summary", "-i", tmp_path / "in", "-o", tmp_path / "o")
    assert read_tsv(tmp_path / "o/busco_summary.tsv")[0]["M"] == "0"
    time.sleep(0.01)
    (run_dir / "missing_busco_list.tsv").write_text("# Busco id\n" + "\n".join(GENES[2:]) + "\n")
    run("summary", "-i", tmp_path / "in", "-o", tmp_path / "o")   # same output dir -> same cache
    assert read_tsv(tmp_path / "o/busco_summary.tsv")[0]["M"] == "4"


def _frag_data(root):
    # gene 101: F in every sample (no complete call -> no reference length)
    # gene 102: F in X only, length 10 vs complete 22 -> ratio 0.45
    _write_run(root / "X/run_metazoa_odb12", dict(zip(GENES, "FFSSSS")))
    _write_run(root / "Y/run_metazoa_odb12", dict(zip(GENES, "FSSSSS")))
    _write_run(root / "Z/run_metazoa_odb12", dict(zip(GENES, "FSSSSS")))


def test_min_frag_ratio_consistent_everywhere(tmp_path):
    _frag_data(tmp_path / "in")
    o = tmp_path / "o"
    run("run", "-i", tmp_path / "in", "-o", o, "--min-frag-ratio", "0.5", "--no-plots", "--no-html",
        "--no-cache")
    gs = by(read_tsv(o / "gene_stats.tsv"), "busco_id")
    # unknown ratio -> not qualified; short fragment -> not qualified
    assert float(gs["101at33208"]["occ_SF"]) == 0.0
    assert gs["101at33208"]["n_F_qualified"] == "0"
    assert abs(float(gs["102at33208"]["occ_SF"]) - 2 / 3) < 1e-3
    sf100 = (o / "gene_lists/SF/SF_occ100.txt").read_text().split()
    assert "101at33208" not in sf100 and "102at33208" not in sf100
    # gene_stats: occ_SF * N == n_S + n_F_qualified for every gene
    for r in gs.values():
        n = int(r["n_samples"])
        assert round(float(r["occ_SF"]) * n) == int(r["n_S"]) + int(r["n_F_qualified"])
    # sample_stats pct_SF_qualified == column sum of binary_SF / n_genes
    mat = read_tsv(o / "matrices/binary_SF.tsv")
    for s in read_tsv(o / "sample_stats.tsv"):
        col = sum(int(m[s["sample"]]) for m in mat)
        assert abs(float(s["pct_SF_qualified"]) - col / int(s["n_genes"])) < 1e-3


def test_extract_respects_min_frag_ratio(tmp_path):
    _frag_data(tmp_path / "in")
    o = tmp_path / "o"
    run("extract", "-i", tmp_path / "in", "-o", o, "--set", "SF", "--occupancy", "0.5",
        "--min-frag-ratio", "0.5", "--no-cache")
    faa = next(o.rglob("102at33208.faa")).read_text()
    assert ">X" not in faa and ">Y" in faa
    man = read_tsv(next(o.rglob("extraction_manifest.tsv")))
    assert any(r["sample"] == "X" and r["busco_id"] == "102at33208"
               and r["result"] == "F_below_min_frag_ratio" for r in man)


def test_newest_odb_run_preferred(tmp_path):
    _write_run(tmp_path / "in/P/run_metazoa_odb10", dict(zip(GENES, "SSSSSS")), lineage="metazoa_odb10")
    _write_run(tmp_path / "in/P/run_metazoa_odb12", dict(zip(GENES, "SSSSSM")), lineage="metazoa_odb12")
    _write_run(tmp_path / "in/Q/run_metazoa_odb12", dict(zip(GENES, "SSSSMM")), lineage="metazoa_odb12")
    run("summary", "-i", tmp_path / "in", "-o", tmp_path / "o", "-l", "metazoa", "--no-cache")
    rows = by(read_tsv(tmp_path / "o/busco_summary.tsv"), "sample")
    assert rows["P"]["lineage"] == "metazoa_odb12"


def test_compleasm_without_summary(tmp_path):
    fixtures.build_compleasm(tmp_path / "in")
    (tmp_path / "in/cmpH/summary.txt").unlink()
    run("summary", "-i", tmp_path / "in", "-o", tmp_path / "o", "--no-cache")
    r = read_tsv(tmp_path / "o/busco_summary.tsv")[0]
    assert (r["sample"], r["tool"], r["S"]) == ("cmpH", "compleasm", "3")


def test_occ_key_never_collides():
    assert analysis.occ_key(0.9) == "occ090"
    assert analysis.occ_key(0.666) != analysis.occ_key(0.667)


def test_compleasm_current_format(tmp_path):
    """compleasm 0.2.x: no '#' header, Interspaced/Retrocopy, one summary block per lineage."""
    d = tmp_path / "in/cmp2"
    (d / "metazoa_odb10").mkdir(parents=True)
    (d / "summary.txt").write_text(
        "## lineage: eukaryota_odb10\nS:10.00%, 1\nD:0.00%, 0\nF:0.00%, 0\nI:0.00%, 0\nM:90.00%, 9\nN:10\n"
        "## lineage: metazoa_odb10\nS:33.33%, 2\nD:16.67%, 1\nF:16.67%, 1\nI:16.67%, 1\nM:16.67%, 1\nN:6\n")
    rows = ["Gene\tStatus\tSequence\tGene Start\tGene End\tStrand\tScore\tLength\n"]
    for g, st in zip(GENES, ["Single", "Single", "Retrocopy", "Fragmented", "Interspaced", "Missing"]):
        rows.append(f"{g}\tMissing\n" if st == "Missing" else f"{g}\t{st}\tchr1\t10\t100\t+\t50\t30\n")
    (d / "metazoa_odb10" / "full_table.tsv").write_text("".join(rows))
    run("summary", "-i", tmp_path / "in", "-o", tmp_path / "o", "--no-cache")
    r = read_tsv(tmp_path / "o/busco_summary.tsv")[0]
    assert (r["tool"], r["n_total"], r["S"], r["D"], r["F"], r["M"]) == ("compleasm", "6", "2", "1", "2", "1")
    assert r["state"] == "OK"   # full_table agrees with the metazoa block, not the first block
