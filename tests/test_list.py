"""`list` writes the result-folder names; `--samples FILE` restricts any command to them."""

from busco2summary.cli import main

from fixtures import GENES, _write_run
from test_matrix import read_tsv


def _data(root):
    for name, calls in (("X.fna", "SSFDMS"), ("Y", "SFMMMS"), ("Z", "SMMMMS")):
        _write_run(root / name / "run_metazoa_odb12", dict(zip(GENES, calls)))


def run(*argv):
    assert main([str(a) for a in argv] + ["-q"] * (argv[0] == "matrix")) == 0


def test_list_names_and_long(tmp_path):
    _data(tmp_path / "in")
    run("list", "-i", tmp_path / "in", "-o", tmp_path / "s.txt")
    assert (tmp_path / "s.txt").read_text().split() == ["X", "Y", "Z"]      # .fna stripped (no_ext)
    run("list", "-i", tmp_path / "in", "--long", "-o", tmp_path / "l.tsv")
    rows = read_tsv(tmp_path / "l.tsv")
    assert rows[0][0] == "#sample" and rows[1][:2] == ["X", "X.fna"]
    assert rows[1][3] == "run_metazoa_odb12" and rows[1][4].endswith("X.fna")


def test_samples_file_restricts_matrix(tmp_path, capsys):
    _data(tmp_path / "in")
    # names, a folder name and a path are all accepted; unknown names are reported
    (tmp_path / "s.txt").write_text(f"X\n# comment\n{tmp_path / 'in' / 'Z'}/\nnot_there\n")
    run("matrix", "-i", tmp_path / "in", "--samples", tmp_path / "s.txt", "-o", tmp_path / "m.tsv")
    assert read_tsv(tmp_path / "m.tsv")[0] == ["busco_id", "X", "Z"]
    assert "not_there" in capsys.readouterr().err
    # the --long output can be fed back as is
    run("list", "-i", tmp_path / "in", "--long", "--exclude", tmp_path / "s.txt", "-o", tmp_path / "l.tsv")
    run("matrix", "-i", tmp_path / "in", "--samples", tmp_path / "l.tsv", "-o", tmp_path / "m2.tsv")
    assert read_tsv(tmp_path / "m2.tsv")[0] == ["busco_id", "Y"]


def test_busco_inside_assembly_folders_named_by_sample(tmp_path):
    # BUSCO run in each SPAdes folder: every result folder has the same name
    for s in ("sampleA", "sampleB"):
        _write_run(tmp_path / "in" / s / "spades" / "busco_scaffolds" / "run_metazoa_odb12",
                   dict(zip(GENES, "SSFDMS")))
    run("list", "-i", tmp_path / "in", "-o", tmp_path / "s.txt")
    assert (tmp_path / "s.txt").read_text().split() == ["sampleA", "sampleB"]
    run("summary", "-i", tmp_path / "in", "-o", tmp_path / "out")
    rows = read_tsv(tmp_path / "out" / "busco_summary.tsv")
    assert [r[0] for r in rows[1:]] == ["sampleA", "sampleB"]
