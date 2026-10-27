"""`busco`: find assemblies, name samples like busco2summary_smkv2 -bfsd, run (fake) BUSCO, resume, summarise."""

import sys
from pathlib import Path

from busco2summary.cli import main

from fixtures import GENES
from test_matrix import read_tsv

FAKE = f"""#!{sys.executable}
import sys
sys.path.insert(0, {str(Path(__file__).parent)!r})
from pathlib import Path
from fixtures import GENES, _write_run
a = sys.argv[1:]
if a == ["--version"]:
    print("BUSCO 6.0.0"); sys.exit(0)
opt = dict(zip(a[::2], a[1::2]))
out = Path(opt["--out_path"])
with open(out.parent / "calls.txt", "a") as fh:
    fh.write(" ".join(a) + "\\n")
if "bad" in opt["-o"]:
    sys.exit(3)
_write_run(out / opt["-o"] / "run_metazoa_odb12", dict(zip(GENES, "SSFDMS")))
"""


def _fasta(p: Path) -> Path:
    p.parent.mkdir(parents=True, exist_ok=True)
    p.write_text(">c1\nACGT\n")
    return p


def _assemblies(root: Path) -> None:
    _fasta(root / "spades" / "00-Spades_rna" / "JLJC1_Eun1_1201_spades_rna" / "transcripts.fasta")
    d = root / "spades" / "00-Spades_wgs" / "23002-1_spades_wgs"
    _fasta(d / "contigs.fasta")
    _fasta(d / "scaffolds.fasta")
    _fasta(d / "K21" / "final_contigs.fasta")                 # inside an assembler folder: not a sample
    (d / "spades.log").write_text("")
    _fasta(root / "spades" / "00-Spades_wgs" / "M4_spades" / "contigs.fasta")
    _fasta(root / "spades" / "00-Spades_wgs" / "M4_spades_wgs" / "contigs.fasta")
    (root / "spades" / "00-Spades_wgs" / "failed_spades").mkdir()
    (root / "spades" / "00-Spades_wgs" / "failed_spades" / "spades.log").write_text("")
    _fasta(root / "spades" / "00-Spades_wgs" / "failed_spades" / "K33" / "final_contigs.fasta")
    _fasta(root / "trinity" / "24003-36Z_trinity.Trinity.fasta")
    (root / "trinity" / "24003-36Z").mkdir()                   # empty Trinity work folder
    _fasta(root / "megahit" / "bad_megahit" / "final.contigs.fa")


def _busco(tmp_path):
    exe = tmp_path / "bin" / "busco"
    exe.parent.mkdir()
    exe.write_text(FAKE)
    exe.chmod(0o755)
    return exe


def _plan(out):
    return {r[0]: r[1:] for r in read_tsv(out / "busco_inputs.tsv")[1:]}


def test_names_modes_dry_run(tmp_path, capsys):
    _assemblies(tmp_path / "a")
    out = tmp_path / "out"
    assert main(["busco", "-i", str(tmp_path / "a"), "-o", str(out), "-l", "metazoa_odb12",
                 "--busco", str(_busco(tmp_path)), "-n"]) == 0
    p = _plan(out)
    assert sorted(p) == ["23002-1", "24003-36Z", "JLJC1_Eun1_1201", "M4_spades", "M4_spades_wgs", "bad"]
    assert p["JLJC1_Eun1_1201"][:2] == ["transcriptome", "rnaspades"]
    assert p["24003-36Z"][:2] == ["transcriptome", "trinity"]
    assert p["23002-1"][0] == "genome" and p["23002-1"][2].endswith("23002-1_spades_wgs/contigs.fasta")
    o = capsys.readouterr()
    assert "-o JLJC1_Eun1_1201 --out_path" in o.out and "failed_spades" in o.err
    assert not (tmp_path / "calls.txt").exists()               # dry run: BUSCO not started
    # --fasta scaffolds.fasta, --suffix alone keeps folder names
    main(["busco", "-i", str(tmp_path / "a" / "spades"), "-o", str(out), "-l", "x", "--busco",
          str(tmp_path / "bin" / "busco"), "-n", "--fasta", "scaffolds.fasta", "--suffix"])
    assert list(_plan(out)) == ["23002-1_spades_wgs"]


def test_run_resume_failures_summary(tmp_path):
    _assemblies(tmp_path / "a")
    out = tmp_path / "out"
    argv = ["busco", "-i", str(tmp_path / "a"), "-o", str(out), "-l", "metazoa_odb12",
            "--busco", str(_busco(tmp_path)), "-J", "2", "-c", "1", "-q"]
    assert main(argv) == 1                                     # 'bad' failed
    calls = (tmp_path / "calls.txt").read_text().splitlines()
    assert len(calls) == 6 and any("-m transcriptome" in c and "-o 24003-36Z" in c for c in calls)
    failed = read_tsv(out / "busco_failed.tsv")
    assert failed[1][:2] == ["bad", "3"] and (out / "logs" / "bad.busco.log").exists()
    s = read_tsv(out / "busco_summary.tsv")
    assert sorted(r[0] for r in s[1:]) == ["23002-1", "24003-36Z", "JLJC1_Eun1_1201", "M4_spades",
                                           "M4_spades_wgs"]
    assert len(s[0]) > 10 and len(GENES) == int(s[1][s[0].index("n_total")])
    # re-run: finished samples skipped, only the failed one is tried again
    main(argv)
    calls = (tmp_path / "calls.txt").read_text().splitlines()
    assert len(calls) == 7 and "-o bad" in calls[-1]
    # --samples restricts; a list file (the plan) can be fed back
    (tmp_path / "s.txt").write_text("23002-1\n")
    main(argv + ["--samples", str(tmp_path / "s.txt"), "--busco-args=--force"])
    assert len((tmp_path / "calls.txt").read_text().splitlines()) == 7   # already done -> nothing
    (out / "23002-1" / "short_summary.specific.metazoa_odb12.23002-1.txt").unlink()
    for f in (out / "23002-1" / "run_metazoa_odb12").glob("short_summary*"):
        f.unlink()
    main(["busco", "--list", str(out / "busco_inputs.tsv"), "-o", str(out), "-l", "metazoa_odb12",
          "--busco", str(tmp_path / "bin" / "busco"), "--samples", str(tmp_path / "s.txt"), "-q"])
    last = (tmp_path / "calls.txt").read_text().splitlines()[-1]
    assert "-o 23002-1" in last and last.endswith("-f")      # unfinished folder: start over
