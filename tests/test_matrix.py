"""`matrix` command: one file, same gene selection as `genelists`, nothing else written."""

import csv
import gzip

from busco2summary.cli import main

from fixtures import GENES, _write_run


def read_tsv(p):
    op = gzip.open if str(p).endswith(".gz") else open
    with op(p, "rt") as fh:
        return list(csv.reader(fh, delimiter="\t"))


def _data(root):
    _write_run(root / "X/run_metazoa_odb12", dict(zip(GENES, "SSFDMS")))
    _write_run(root / "Y/run_metazoa_odb12", dict(zip(GENES, "SFMMMS")))
    _write_run(root / "Z/run_metazoa_odb12", dict(zip(GENES, "SMMMMS")))


def run(*argv):
    assert main([str(a) for a in argv] + ["-q"]) == 0


def test_matrix_status_default(tmp_path):
    _data(tmp_path / "in")
    run("matrix", "-i", tmp_path / "in", "-o", tmp_path / "m.tsv")
    rows = read_tsv(tmp_path / "m.tsv")
    assert rows[0] == ["busco_id", "X", "Y", "Z"]
    assert rows[1] == [GENES[0], "S", "S", "S"]
    assert len(rows) == 1 + 6
    assert sorted(p.name for p in tmp_path.iterdir()) == ["in", "m.tsv"]   # no folder, cache or log


def test_matrix_occupancy_matches_genelists(tmp_path):
    _data(tmp_path / "in")
    run("matrix", "-i", tmp_path / "in", "--set", "SF", "-t", "60", "--format", "binary",
        "-o", tmp_path / "m.tsv.gz")
    rows = read_tsv(tmp_path / "m.tsv.gz")
    run("genelists", "-i", tmp_path / "in", "-o", tmp_path / "gl", "-t", "60", "--sets", "SF", "--no-cache")
    want = (tmp_path / "gl/gene_lists/SF/SF_occ060.txt").read_text().split()
    assert [r[0] for r in rows[1:]] == want == [GENES[0], GENES[1], GENES[5]]
    assert rows[2] == [GENES[1], "1", "1", "0"]


def test_matrix_min_taxa_and_transpose(tmp_path):
    _data(tmp_path / "in")
    run("matrix", "-i", tmp_path / "in", "--set", "S", "--min-taxa", "3", "--samples-as-rows",
        "--format", "count", "-o", tmp_path / "m.tsv")
    rows = read_tsv(tmp_path / "m.tsv")
    assert rows[0] == ["sample", GENES[0], GENES[5]]
    assert [r[0] for r in rows[1:]] == ["X", "Y", "Z"]


def test_matrix_ids_for_given_genes(tmp_path):
    _data(tmp_path / "in")
    (tmp_path / "ids.txt").write_text(f"{GENES[1]}\n{GENES[2]}\nnot_a_busco\n")
    run("matrix", "-i", tmp_path / "in", "--genes", tmp_path / "ids.txt", "--format", "ids",
        "-o", tmp_path / "ids.tsv")
    rows = read_tsv(tmp_path / "ids.tsv")
    assert rows[0] == ["sample", "busco_id", "status"]
    assert rows[1:] == [["X", GENES[1], "S"], ["X", GENES[2], "F"], ["Y", GENES[1], "F"]]
