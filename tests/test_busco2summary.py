import csv
import json
import math
from pathlib import Path

import pytest

from busco2summary import analysis, loader, parsers
from busco2summary.cli import main
from busco2summary.extract import cds_from_cigar, clean_cds, revcomp, translate

import fixtures


@pytest.fixture(scope="module")
def data(tmp_path_factory):
    root = tmp_path_factory.mktemp("busco")
    exp = fixtures.build(root)
    fixtures.build_compleasm(root)
    return root, exp


def read_tsv(p):
    with open(p) as fh:
        return list(csv.DictReader(fh, delimiter="\t"))


# ---------------------------------------------------------------- parsing

def test_discovery_and_layouts(data):
    root, exp = data
    cands = {c.dir_name: c for c in loader.discover([root])}
    assert set(cands) == {"sampleA", "sampleB.fna", "tranC", "run_v3D", "autoE", "emptyF",
                          "noftG", "cmpH"}
    assert cands["run_v3D"].layout == "busco_v3"
    assert cands["cmpH"].layout == "compleasm"
    assert len(cands["autoE"].runs) == 2


def test_counts_match_truth(data):
    root, exp = data
    samples = {s.name: s for s in loader.load_all([str(root)], name_style="no_ext", threads=2)}
    for name, counts in exp.items():
        assert samples[name].counts() == counts, name
    assert samples["v3D"].state == "OK"                  # run_ prefix stripped by no_ext
    assert samples["emptyF"].state == "EMPTY"
    assert samples["noftG"].state == "NO_FULL_TABLE"
    assert samples["autoE"].lineage == "metazoa_odb12"   # specific run chosen
    assert samples["autoE"].other_runs == ["run_eukaryota_odb10"]
    assert samples["sampleB"].info["pct_E"] == 20.0      # from json
    assert samples["sampleB"].info["gene_predictor"] == "miniprot"
    c = samples["cmpH"]
    assert c.tool == "compleasm" and c.counts() == {"S": 3, "D": 1, "F": 1, "M": 1}


def test_lineage_selection(data):
    root, _ = data
    s = {x.name: x for x in loader.load_all([str(root)], lineage="eukaryota", name_style="no_ext")}
    assert s["autoE"].lineage == "eukaryota_odb10"
    assert s["sampleA"].state == "LINEAGE_ABSENT"


def test_mismatch_detected(tmp_path):
    run = tmp_path / "x" / "run_metazoa_odb12"
    fixtures._write_run(run, fixtures.STATUS_A)
    txt = (run / "short_summary.txt").read_text().replace("\t3\tComplete and single-copy",
                                                          "\t4\tComplete and single-copy")
    (run / "short_summary.txt").write_text(txt)
    s = loader.load_all([str(tmp_path)])[0]
    assert s.state == "MISMATCH"


def test_one_line_regex():
    r = parsers._ONE_LINE.search("C:89.3%[S:87.9%,D:1.3%],F:4.6%,M:6.1%,n:672,E:33.3%")
    assert r.group(6) == "672" and r.group(7) == "33.3"


def test_header_location():
    assert parsers.header_location("g_1_0:abc|ctg1:400-100|- desc") == ("ctg1", 100, 400)
    assert parsers.header_location("g|NODE_1_length_5:77-391|+") == ("NODE_1_length_5", 77, 391)


# ---------------------------------------------------------------- analysis

def test_occupancy_exact_rule(data):
    root, _ = data
    samples = [s for s in loader.load_all([str(root)], name_style="no_ext")
               if s.lineage == "metazoa_odb12" and s.state != "EMPTY"]
    ds = analysis.build_dataset(samples, "metazoa_odb12")
    n = ds.n
    table, lists = analysis.occupancy_lists(ds, [0.5, 1.0], ["S", "SF", "SDF"], [3])
    for t in (0.5, 1.0):
        need = math.ceil(t * n - 1e-9)
        for k in ("S", "SF", "SDF"):
            want = [g for g in ds.genes if sum(ds.rows[g][j] in k for j in range(n)) >= need]
            assert lists[f"{k}_occ{int(t * 100):03d}"] == want
    # F rescues loci: SF list is a superset of S
    assert set(lists["S_occ050"]) <= set(lists["SF_occ050"])
    gs = {r["busco_id"]: r for r in analysis.gene_stats(ds)}
    for g, r in gs.items():
        assert r["n_S"] + r["n_D"] + r["n_F"] + r["n_M"] == n


def test_min_frag_ratio(data):
    root, _ = data
    samples = [s for s in loader.load_all([str(root)], name_style="no_ext")
               if s.name in ("sampleA", "tranC")]   # gene 104: F (len 10) in A, S (len 22) in C
    ds0 = analysis.build_dataset(samples, "m", min_frag_ratio=0)
    ds9 = analysis.build_dataset(samples, "m", min_frag_ratio=0.9)
    g = fixtures.GENES[3]                       # F in sampleA, length 10 vs complete 22
    assert ds0.present_count(g, "SF") == 2 and ds9.present_count(g, "SF") == 1


def test_taxon_drop_and_grid(data):
    root, _ = data
    samples = [s for s in loader.load_all([str(root)], name_style="no_ext")
               if s.lineage == "metazoa_odb12" and s.state != "EMPTY"]
    ds = analysis.build_dataset(samples, "metazoa_odb12")
    drop = analysis.taxon_drop(ds, "S", [0.9], steps=2)
    assert drop[0]["n_taxa"] == ds.n and drop[1]["n_taxa"] == ds.n - 1
    assert drop[1]["genes_occ90"] >= drop[0]["genes_occ90"]
    g = analysis.grid(ds, "SF", [0, 0.5], [0.5])
    assert g[0]["n_taxa"] >= g[1]["n_taxa"]


# ---------------------------------------------------------------- CLI end to end

def test_cli_run(data, tmp_path):
    root, exp = data
    out = tmp_path / "out"
    assert main(["run", "-i", str(root), "-o", str(out), "-q", "--no-plots", "--long-table",
                 "-t", "50", "100", "--extract", "SF_occ050", "--seqtype", "both",
                 "--layout", "both", "--dup", "best"]) == 0
    summ = {r["sample"]: r for r in read_tsv(out / "busco_summary.tsv")}
    for name, c in exp.items():
        assert int(summ[name]["S"]) == c["S"] and int(summ[name]["F"]) == c["F"]
    assert (out / "busco2summary_report.html").stat().st_size > 1000
    assert (out / "gene_lists" / "SF" / "SF_occ050.txt").exists()
    assert (out / "matrices" / "status_matrix.tsv").exists()
    assert (out / "matrices" / "long_table.tsv.gz").exists()
    # best duplicated copy = higher score = longer copy (ends with W)
    faa = dict(parsers.read_fasta(out / "sequences" / "SF_occ050" / "by_sample_faa" / "sampleA.faa"))
    assert all(not s.endswith("W") for s in faa.values())       # SF set excludes D
    man = read_tsv(out / "sequences" / "SF_occ050" / "extraction_manifest.tsv")
    nb = [r for r in man if r["sample"] == "sampleB" and r["seq_type"] == "nucleotide"]
    assert nb and not any(r["result"].startswith("ok") for r in nb)   # miniprot: no .fna
    rp = json.loads((out / "run_parameters.json").read_text())
    assert rp["cmd"] == "run"


def test_cli_extract_dup_best(data, tmp_path):
    root, _ = data
    out = tmp_path / "ex"
    assert main(["extract", "-i", str(root), "-o", str(out), "-q", "--set", "SDF",
                 "--occupancy", "0.1", "--layout", "sample"]) == 0
    faa = dict(parsers.read_fasta(out / "SDF_occ010" / "by_sample_faa" / "sampleA.faa"))
    dgene = fixtures.GENES[2]                                     # D in sampleA
    assert faa[dgene].endswith("W")                               # score 250 copy chosen


def test_sample_sheet_v03_format(data, tmp_path):
    root, _ = data
    sheet = tmp_path / "map.tsv"
    sheet.write_text("#original_name\t#new_name\nsampleA\tSpecies_one\ntranC\tSpecies_two\n")
    out = tmp_path / "o"
    assert main(["summary", "-i", str(root), "-s", str(sheet), "-o", str(out), "-q"]) == 0
    names = [r["sample"] for r in read_tsv(out / "busco_summary.tsv")]
    assert names == ["Species_one", "Species_two"]


# ---------------------------------------------------------------- CDS rebuilding

def test_cds_from_cigar_plus_and_minus():
    exon1, intron, exon2, nxt = "ATGGCC", "GTAAGTTTAG", "AAATTT", "TAA"
    contig = "NN" + exon1 + intron + exon2 + nxt + "NN"
    cds = cds_from_cigar(contig, 2, 2 + len(exon1 + intron + exon2), "+", "2M10N2M")
    assert cds == exon1 + exon2 + nxt and translate(cds) == "MAKF*"
    rc = revcomp(contig)
    cds2 = cds_from_cigar(rc, len(rc) - (2 + len(exon1 + intron + exon2)), len(rc) - 2, "-", "2M10N2M")
    assert cds2 == cds
    # phase-1 intron: codon split 1 + 2 around the intron
    contig = "ATGG" + "GTAAGTTTAG" + "CCAAA" + "TGA"
    cds = cds_from_cigar(contig, 0, 4 + 10 + 5, "+", "1M13U1M")
    assert translate(cds) == "MAK*"
    # frameshift F / G bases are skipped
    contig = "ATGA" + "GCC" + "TGA"
    assert translate(cds_from_cigar(contig, 0, 7, "+", "1M1F1M")) == "MA*"


def test_clean_cds():
    assert clean_cds("ATGTAAGCCTGA", "strip", "M*A") == "ATGTAAGCC"
    assert clean_cds("ATGTAAGCCTGA", "strip", "MA") == "ATGGCC"
    assert clean_cds("ATGTAAGCCTGA", "X") == "ATGNNNGCC"
    assert clean_cds("ATGTAAGCCTGA", "keep") == "ATGTAAGCCTGA"
