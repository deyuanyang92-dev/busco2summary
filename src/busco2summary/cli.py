"""Command line interface.

busco2summary list     -i BUSCO_DIR > s.txt     # names of all result folders (for --samples)
busco2summary matrix   -i BUSCO_DIR > m.tsv     # one matrix only, no output folder
busco2summary run      -i BUSCO_DIR -o OUT      # everything (summary, matrices, lists, plots, HTML)
busco2summary scan     -i BUSCO_DIR -o OUT      # discovery only: what is where
busco2summary summary  -i BUSCO_DIR -o OUT      # per-sample summary table only
busco2summary genelists -i BUSCO_DIR -o OUT -t 50 60 90 --sets S SF
busco2summary extract  -i BUSCO_DIR -o OUT --set SF --occupancy 0.9 --seqtype protein
"""

from __future__ import annotations

import argparse
import json
import sys
import time
from collections import Counter
from pathlib import Path
from typing import Dict, List, Optional

from . import __version__, analysis, log, parsers
from .analysis import QCThresholds
from .loader import Sample, load_all
from .log import error, info, warn
from .summary import SUMMARY_FIELDS, dataset_report, summary_row
from .writers import write_list, write_matrix, write_table, write_xlsx

EPILOG = """
examples:
  busco2summary list -i busco_results/ > samples.txt          # edit, then:
  busco2summary matrix -i busco_results/ --samples samples.txt > status_matrix.tsv
  busco2summary matrix -i busco_results/ > status_matrix.tsv
  busco2summary matrix -i busco_results/ --set SF -t 50 --format binary -o SF_occ50.tsv
  busco2summary matrix -i busco_results/ --genes my_ids.txt --format ids > present_ids.tsv
  busco2summary run -i busco_results/ -o b2s_out -j 16
  busco2summary run -i runA/ runB/ --sample-sheet samples.tsv --lineage metazoa_odb12 -o out
  busco2summary genelists -i busco_results/ -o out -t 50 60 70 80 90 --sets S SF
  busco2summary extract -i busco_results/ -o out --set SF --occupancy 0.6 --seqtype both \\
        --dup best --layout both
  busco2summary extract -i busco_results/ -o out --genes my_ids.txt --genome-from-summary

status sets: S (single-copy), SD (complete), SF (single + fragmented), SDF (any hit)
docs: https://github.com/deyuanyang92-dev/busco2summary
"""


def _thresholds(vals: List[float]) -> List[float]:
    out = []
    for v in vals:
        v = float(v)
        out.append(v / 100 if v > 1 else v)
    return sorted(set(round(x, 6) for x in out))


def _read_set(path: Optional[str]) -> Optional[set]:
    if not path:
        return None
    # first column of each line; a path is reduced to its folder name (so `list --long` output
    # or a list of result folders works as well as plain sample names)
    out = set()
    for ln in Path(path).read_text().splitlines():
        if ln.strip() and not ln.startswith("#"):
            out.add(Path(ln.split("\t")[0].strip().rstrip("/")).name)
    return out


def add_input_args(p: argparse.ArgumentParser) -> None:
    g = p.add_argument_group("input")
    g.add_argument("-i", "--input", nargs="+", default=[], metavar="DIR",
                   help="folder(s) holding BUSCO/compleasm results (searched recursively)")
    g.add_argument("-s", "--sample-sheet", metavar="TSV",
                   help="sample sheet: columns dir_name/new_name/path/lineage (any synonym, "
                        "busco_2_matrix v0.3 '#' headers accepted) or 2 columns 'dir new_name'")
    g.add_argument("-l", "--lineage", default="",
                   help="lineage to analyse (e.g. metazoa_odb12; 'metazoa' matches any odb); "
                        "'all' = every lineage in its own sub-folder. the most common lineage is used when empty")
    g.add_argument("--name-style", choices=["dir", "no_ext", "stem", "accession"], default="no_ext",
                   help="sample name from folder name: as is | strip .fna/.fasta.. | text before "
                        "first dot | GCA/GCF/SRR accession")
    g.add_argument("--replace-dots", action="store_true",
                   help="replace '.' in sample names by '_' (busco_2_matrix v0.3 default)")
    g.add_argument("--include", "--samples", metavar="FILE",
                   help="only these samples: one name or result folder per line "
                        "(e.g. edited output of `busco2summary list`)")
    g.add_argument("--exclude", metavar="FILE", help="drop these samples (same format as --include)")
    g.add_argument("--max-depth", type=int, default=4, help="folder search depth")
    g.add_argument("-j", "--threads", type=int, default=8, help="parallel readers")
    g.add_argument("--no-cache", action="store_true", help="do not use/write the parse cache")


def add_filter_args(p: argparse.ArgumentParser) -> None:
    g = p.add_argument_group("sample filtering (matrices / gene lists)")
    g.add_argument("--keep-empty", action="store_true",
                   help="keep samples with no BUSCO at all in matrices (default: drop)")
    g.add_argument("--min-completeness", type=float, default=0.0, metavar="PCT",
                   help="drop samples whose S+D+F %% is below PCT from matrices/lists "
                        "(values <= 1 are fractions, as for -t)")
    g.add_argument("--min-frag-ratio", type=float, default=0.0, metavar="R",
                   help="count a fragmented BUSCO as present only if its length >= R x reference "
                        "length (0 = every F counts; e.g. 0.5)")
    g.add_argument("--frag-ref", choices=["dataset", "lineage"], default="dataset",
                   help="reference length for fragment ratios: median complete length of the "
                        "BUSCO in this dataset (same predictor) or lineage proteins (--lineage-dir)")
    g.add_argument("--lineage-dir", metavar="DIR",
                   help="BUSCO lineage dataset (or a folder of datasets) for reference lengths; "
                        "default: path from short_summary.json if readable, else dataset median")


def add_output_args(p: argparse.ArgumentParser) -> None:
    p.add_argument("-o", "--outdir", required=True, help="output folder (created; files overwritten)")
    p.add_argument("-q", "--quiet", action="store_true")


def build_parser() -> argparse.ArgumentParser:
    ap = argparse.ArgumentParser(
        prog="busco2summary", formatter_class=argparse.RawDescriptionHelpFormatter,
        description="Detailed multi-sample BUSCO summaries, S/D/F/M matrices, per-gene F/D "
                    "statistics, occupancy gene-ID lists and sequence extraction.",
        epilog=EPILOG)
    ap.add_argument("--version", action="version", version=f"busco2summary {__version__}")
    sub = ap.add_subparsers(dest="cmd", metavar="command")

    def mk(name, hlp):
        sp = sub.add_parser(name, help=hlp.replace("%", "%%"), description=hlp.replace("%", "%%"),
                            formatter_class=argparse.ArgumentDefaultsHelpFormatter)
        add_input_args(sp)
        add_output_args(sp)
        return sp

    sp = sub.add_parser("list", formatter_class=argparse.ArgumentDefaultsHelpFormatter,
                        help="names of all BUSCO result folders (one per line), for --samples",
                        description="List every BUSCO/compleasm result folder found under -i, one sample "
                                    "name per line (fast: no table is read). Edit the list and pass it "
                                    "back with --samples FILE to analyse only those samples.")
    add_input_args(sp)
    sp.add_argument("-o", "--output", default="-", help="output file ('-' = stdout)")
    sp.add_argument("--long", action="store_true",
                    help="also folder name, layout, lineage run(s) and path (tab-separated)")
    mk("scan", "list discovered BUSCO results (layout, lineage, version, state)")
    sp = mk("summary", "per-sample summary table (counts, percentages, E, assembly stats, QC flags)")
    sp.add_argument("--xlsx", action="store_true", help="also write busco_summary.xlsx")

    for name, hlp in (("run", "full analysis: summary + matrices + statistics + gene lists + plots + HTML"),
                      ("genelists", "occupancy-based BUSCO ID lists only")):
        sp = mk(name, hlp)
        add_filter_args(sp)
        g = sp.add_argument_group("gene lists")
        g.add_argument("-t", "--thresholds", nargs="+", type=float,
                       default=[t * 100 for t in analysis.DEFAULT_THRESHOLDS],
                       help="occupancy thresholds, %% or fraction")
        g.add_argument("--sets", nargs="+", default=["S", "SD", "SF", "SDF"],
                       choices=["S", "SD", "SF", "SDF"], help="status sets for gene lists")
        g.add_argument("--min-taxa", nargs="*", type=int, default=[3],
                       help="also list BUSCOs present in >= K samples")
        g.add_argument("--drop-steps", type=int, default=10,
                       help="greedy taxon-drop steps reported (0 = off)")
        if name == "run":
            o = sp.add_argument_group("outputs")
            o.add_argument("--transpose", action="store_true",
                           help="matrices with samples as rows (default: BUSCOs as rows)")
            o.add_argument("--no-plots", action="store_true")
            o.add_argument("--plot-format", nargs="+", default=["png", "pdf"],
                           choices=["png", "pdf", "svg"])
            o.add_argument("--no-html", action="store_true")
            o.add_argument("--xlsx", action="store_true", help="also write an Excel workbook")
            o.add_argument("--long-table", action="store_true",
                           help="write long (tidy) table sample x BUSCO (gzipped)")
            q = sp.add_argument_group("QC flag thresholds (%)")
            q.add_argument("--qc-min-c", type=float, default=50)
            q.add_argument("--qc-max-d", type=float, default=10)
            q.add_argument("--qc-max-f", type=float, default=15)
            q.add_argument("--qc-max-m", type=float, default=30)
            q.add_argument("--qc-max-e", type=float, default=50)
            x = sp.add_argument_group("sequence extraction (optional)")
            x.add_argument("--extract", nargs="*", metavar="SET_OCC",
                           help="extract sequences for list(s) like SF_occ090 S_occ050 "
                                "(no value = S_occ050 and SF_occ050)")
            _extract_opts(x)

    sp = sub.add_parser("matrix", formatter_class=argparse.ArgumentDefaultsHelpFormatter,
                        help="one gene x sample matrix to stdout or a file (no output folder, "
                             "no cache, no plots)",
                        description="Write one matrix and nothing else. Rows = BUSCOs, columns = samples "
                                    "(--samples-as-rows to flip). Messages go to stderr.")
    add_input_args(sp)
    add_filter_args(sp)
    g = sp.add_argument_group("matrix")
    g.add_argument("-o", "--output", default="-", help="output file ('-' = stdout; '.gz' = gzipped)")
    g.add_argument("--format", choices=["status", "binary", "count", "ids"], default="status",
                   help="status: S/D/F/M letters | binary: 1 = present in --set | count: copy number "
                        "| ids: long table sample<TAB>busco_id<TAB>status, present calls only")
    g.add_argument("--set", default="SDF", choices=["S", "SD", "SF", "SDF"],
                   help="statuses counted as present (binary, ids, occupancy filter)")
    g.add_argument("-t", "--occupancy", type=float, default=0.0,
                   help="keep BUSCOs present in >= ceil(t*N) samples (%% or fraction; 0 = all)")
    g.add_argument("--min-taxa", type=int, default=0, help="keep BUSCOs present in >= K samples")
    g.add_argument("--genes", metavar="FILE", help="only these BUSCO IDs (one per line)")
    g.add_argument("--samples-as-rows", action="store_true", help="samples as rows, BUSCOs as columns")
    g.add_argument("-q", "--quiet", action="store_true")

    sp = mk("extract", "extract per-gene / per-sample FASTA for a BUSCO list")
    add_filter_args(sp)
    g = sp.add_argument_group("which BUSCOs")
    g.add_argument("--genes", metavar="FILE", help="BUSCO IDs, one per line (e.g. from gene_lists/)")
    g.add_argument("--set", default="S", choices=["S", "SD", "SF", "SDF"],
                   help="status set used for occupancy and for sequence sources")
    g.add_argument("--occupancy", type=float, default=0.5, help="occupancy threshold (%% or fraction)")
    _extract_opts(g)
    return ap


def _extract_opts(g) -> None:
    g.add_argument("--seqtype", choices=["protein", "nucleotide", "both"], default="protein")
    g.add_argument("--dup", choices=["best", "longest", "first", "skip"], default="best",
                   help="several copies (D, multi-record files): highest full_table score | "
                        "longest | first | skip duplicated BUSCOs")
    g.add_argument("--layout", choices=["gene", "sample", "both"], default="gene",
                   help="one FASTA per BUSCO (gene), per sample, or both")
    g.add_argument("--header", default="{sample}",
                   help="per-gene FASTA header template; fields {sample} {gene} {status} {dir} {orig}")
    g.add_argument("--sample-header", default="{gene}",
                   help="per-sample FASTA header template")
    g.add_argument("--stop", choices=["strip", "keep", "X"], default="strip",
                   help="protein '*': strip terminal | keep | strip terminal + internal->X")
    g.add_argument("--min-len", type=int, default=0, help="minimum sequence length")
    g.add_argument("--genome-map", metavar="TSV",
                   help="sample<TAB>genome.fasta: rebuild CDS from BUSCO GFF when no .fna "
                        "(miniprot runs)")
    g.add_argument("--genome-from-summary", action="store_true",
                   help="use the input file recorded by BUSCO as genome for CDS rebuilding")
    g.add_argument("--line-width", type=int, default=0, help="FASTA line width (0 = one line)")


# --------------------------------------------------------------------------

def _load(args) -> List[Sample]:
    outdir = getattr(args, "outdir", None)          # 'matrix' has no output folder
    out = Path(outdir) if outdir else None
    if out:
        out.mkdir(parents=True, exist_ok=True)
    log.QUIET = getattr(args, "quiet", False)
    log.set_logfile(out / "busco2summary.log" if out else None)
    info(f"busco2summary {__version__}  command: {' '.join(sys.argv)}")
    if not args.input and not args.sample_sheet:
        error("give -i/--input and/or -s/--sample-sheet")
        sys.exit(2)
    lineage = "" if args.lineage in ("", "all") else args.lineage
    cache = None if args.no_cache or out is None else out / ".cache"
    t0 = time.time()
    samples = load_all(args.input, args.sample_sheet, lineage, args.name_style, args.replace_dots,
                       _read_set(args.include), _read_set(args.exclude), args.threads,
                       args.max_depth, cache)
    info(f"loaded {len(samples)} BUSCO results in {time.time() - t0:.1f}s")
    if not samples:
        error("no BUSCO results found (check -i, --max-depth, --lineage)")
        sys.exit(1)
    return samples


def _write_summary(samples: List[Sample], out: Path, xlsx: bool = False) -> List[dict]:
    rows = [summary_row(s) for s in samples]
    write_table(out / "busco_summary.tsv", rows, SUMMARY_FIELDS, digits=2)
    write_table(out / "busco_summary.csv", rows, SUMMARY_FIELDS, sep=",", digits=2)
    if xlsx:
        write_xlsx(out / "busco_summary.xlsx", {"summary": rows})
    info(f"summary: {out / 'busco_summary.tsv'}")
    return rows


def _lineages(samples: List[Sample], want: str) -> Dict[str, List[Sample]]:
    usable = [s for s in samples if s.state not in ("FAILED", "LINEAGE_ABSENT")]
    by: Dict[str, List[Sample]] = {}
    for s in usable:
        by.setdefault(s.lineage, []).append(s)
    if not by:
        return {}
    if want == "all":
        return by
    if want:
        hit = {k: v for k, v in by.items() if k == want or k.startswith(want)}
        if not hit:
            error(f"--lineage {want}: not found; available: {', '.join(sorted(by))}")
        return hit
    top = Counter({k: len(v) for k, v in by.items()}).most_common()
    if len(top) > 1:
        warn("several lineages found: " + ", ".join(f"{k}={n}" for k, n in top) +
             f"; analysing '{top[0][0]}' (use --lineage NAME or --lineage all)")
    return {top[0][0]: by[top[0][0]]}


def _ref_lengths(args, samples: List[Sample], lineage: str):
    cands = []
    if getattr(args, "lineage_dir", None):
        d = Path(args.lineage_dir)
        cands += [d / lineage, d]
    for s in samples:
        p = s.info.get("lineage_path")
        if p:
            cands.append(Path(str(p)))
            break
    for d in cands:
        if d.is_dir():
            ref = parsers.read_lineage_lengths(d)
            if ref:
                info(f"reference lengths: {len(ref)} BUSCOs from {d}")
                return ref, f"lineage dataset {d.name}"
    return None, ""


def _dataset(args, lineage: str, group: List[Sample]) -> Optional[analysis.Dataset]:
    keep = []
    for s in group:
        if s.state == "EMPTY" and not args.keep_empty:
            warn(f"dropping EMPTY sample {s.name} from matrices (use --keep-empty to keep)")
            continue
        mc = args.min_completeness * 100 if 0 < args.min_completeness <= 1 else args.min_completeness
        if mc > 0 and s.calls:
            n = s.counts()
            tot = sum(n.values()) or 1
            if 100 * (tot - n["M"]) / tot < mc:
                warn(f"dropping {s.name}: S+D+F {100 * (tot - n['M']) / tot:.1f}% "
                     f"< --min-completeness {mc:g}%")
                continue
        if not s.calls:
            warn(f"{s.name}: no per-BUSCO table, not in matrices")
            continue
        keep.append(s)
    if not keep:
        warn(f"lineage {lineage}: no usable samples for matrices")
        return None
    ref, src = _ref_lengths(args, keep, lineage)
    return analysis.build_dataset(keep, lineage, ref, src, args.min_frag_ratio,
                                  getattr(args, "frag_ref", "dataset"))


def _write_lists(args, ds: analysis.Dataset, out: Path) -> List[dict]:
    th = _thresholds(args.thresholds)
    table, lists = analysis.occupancy_lists(ds, th, args.sets, args.min_taxa or [])
    for key, ids in lists.items():
        st = key.split("_")[0]
        write_list(out / "gene_lists" / st / f"{key}.txt", ids)
    write_table(out / "gene_lists" / "occupancy_summary.tsv", table,
                ["set", "rule", "value", "min_samples", "n_genes", "F_gain_genes", "matrix_fill"])
    # wide view: thresholds x sets
    wide = []
    for t in th:
        r = {"occupancy": t, "min_samples": analysis.need(t, ds.n)}
        for k in args.sets:
            r[k] = len(lists[f"{k}_{analysis.occ_key(t)}"])
        if "S" in args.sets and "SF" in args.sets:
            r["F_gain(SF-S)"] = r["SF"] - r["S"]
        wide.append(r)
    write_table(out / "gene_lists" / "occupancy_wide.tsv", wide)
    sc = [0, 0.3, 0.5, 0.6, 0.7, 0.8, 0.9]
    for k in args.sets:
        write_table(out / "gene_lists" / f"grid_{k}.tsv", analysis.grid(ds, k, sc, th))
    if args.drop_steps > 0 and ds.n > 3:
        for k in ("S", "SF"):
            if k in args.sets:
                write_table(out / "gene_lists" / f"taxon_drop_{k}.tsv",
                            analysis.taxon_drop(ds, k, th, args.drop_steps))
    info(f"gene lists: {len(lists)} lists -> {out / 'gene_lists'}")
    for r in wide:
        info("  occ>={:>4.0%} (>= {:>4} samples): ".format(r["occupancy"], r["min_samples"]) +
             "  ".join(f"{k}={r[k]}" for k in args.sets))
    return wide


def _write_matrices(args, ds: analysis.Dataset, out: Path) -> None:
    m = out / "matrices"
    names = ds.names
    tr = args.transpose
    write_matrix(m / "status_matrix.tsv", ds.genes, names, lambda i, j: ds.rows[ds.genes[i]][j],
                 transpose=tr)

    def call(i, j):
        return ds.samples[j].calls.get(ds.genes[i])
    write_matrix(m / "copy_number_matrix.tsv", ds.genes, names,
                 lambda i, j: (call(i, j) or [0, 0])[1], transpose=tr)
    write_matrix(m / "length_matrix.tsv", ds.genes, names,
                 lambda i, j: (call(i, j) or [None, None, None])[2], transpose=tr)
    write_matrix(m / "fragment_ratio_matrix.tsv", ds.genes, names,
                 lambda i, j: ds.frag_ratio[ds.genes[i]][j], transpose=tr)
    for k in ("S", "SD", "SF", "SDF"):
        masks = {g: ds.present_mask(g, k) for g in ds.genes}
        write_matrix(m / f"binary_{k}.tsv", ds.genes, names,
                     lambda i, j: 1 if masks[ds.genes[i]][j] else 0, transpose=tr)
    if args.long_table:
        rows = []
        for j, s in enumerate(ds.samples):
            for g in ds.genes:
                c = s.calls.get(g) or ["M", 0, None, None]
                rows.append({"sample": s.name, "busco_id": g, "status": c[0], "n_hits": c[1],
                             "length": c[2], "score": c[3], "fragment_ratio": ds.frag_ratio[g][j]})
        write_table(m / "long_table.tsv.gz", rows)
    info(f"matrices -> {m}")


def _genome_map(path: Optional[str]) -> Optional[Dict[str, str]]:
    if not path:
        return None
    d = {}
    for l in Path(path).read_text().splitlines():
        p = l.split("\t") if "\t" in l else l.split()
        if len(p) >= 2 and not l.startswith("#"):
            d[p[0]] = p[1]
    return d


def _do_extract(args, ds: analysis.Dataset, genes: List[str], sset: str, out: Path) -> None:
    from . import extract as ex
    statuses = sset.replace("D", "") if args.dup == "skip" else sset
    short_f: Dict[str, set] = {}
    if ds.min_frag_ratio > 0:
        for g in genes:
            if g not in ds.rows:
                continue
            for j, ok in enumerate(ds.present_mask(g, "F")):
                if ds.rows[g][j] == "F" and not ok:
                    short_f.setdefault(ds.samples[j].name, set()).add(g)
    man = ex.extract(ds.samples, genes, out, statuses=statuses, seqtype=args.seqtype,
                     dup_policy="best" if args.dup == "skip" else args.dup, layout=args.layout,
                     header=args.header, sample_header=args.sample_header, stop=args.stop,
                     min_len=args.min_len, line_width=args.line_width,
                     genome_map=_genome_map(args.genome_map),
                     genome_from_summary=args.genome_from_summary, threads=args.threads,
                     short_f=short_f)
    c = Counter(r["result"].split("_of_")[0] for r in man if r["result"] not in ("missing", "not_selected"))
    info(f"extraction {out.name}: " + ", ".join(f"{k}={v}" for k, v in c.most_common()))
    write_list(out / "busco_ids.txt", genes)


def cmd_list(args) -> int:
    from .loader import discover, sample_names
    if not args.input:
        error("give -i/--input")
        return 2
    inc, exc = _read_set(args.include), _read_set(args.exclude)
    rows = []
    cands = discover([Path(i) for i in args.input], max_depth=args.max_depth)
    for c, nm in zip(cands, sample_names(cands, args.name_style, args.replace_dots)):
        if inc and nm not in inc and c.dir_name not in inc:
            continue
        if exc and (nm in exc or c.dir_name in exc):
            continue
        rows.append((nm, c.dir_name, c.layout, ",".join(r.name for r in c.runs), str(c.path)))
    fh = sys.stdout if args.output == "-" else open(args.output, "w")
    try:
        if args.long:
            fh.write("#sample\tdir_name\tlayout\truns\tpath\n")
        for r in rows:
            fh.write(("\t".join(r) if args.long else r[0]) + "\n")
    finally:
        if fh is not sys.stdout:
            fh.close()
    info(f"{len(rows)} BUSCO result folders" + ("" if args.output == "-" else f" -> {args.output}"))
    return 0


def cmd_scan(args) -> int:
    samples = _load(args)
    out = Path(args.outdir)
    rows = [{"sample": s.name, "dir_name": s.dir_name, "state": s.state, "layout": s.layout,
             "tool": s.tool, "lineage": s.lineage, "busco_version": s.info.get("busco_version"),
             "mode": s.info.get("mode"), "gene_predictor": s.info.get("gene_predictor"),
             "n_calls": len(s.calls), "seq_S": len(s.seqs.get("S", {})),
             "seq_D": len(s.seqs.get("D", {})), "seq_F": len(s.seqs.get("F", {})),
             "run_dir": s.run_dir, "other_runs": ",".join(s.other_runs),
             "notes": " | ".join(s.notes)} for s in samples]
    write_table(out / "scan.tsv", rows)
    rep = dataset_report(samples, args.lineage or "(auto)")
    (out / "dataset_report.txt").write_text("\n".join(rep) + "\n")
    print("\n".join(rep))
    return 0


def cmd_summary(args) -> int:
    samples = _load(args)
    out = Path(args.outdir)
    _write_summary(samples, out, args.xlsx)
    rep = dataset_report(samples, args.lineage or "(auto)")
    (out / "dataset_report.txt").write_text("\n".join(rep) + "\n")
    print("\n".join(rep))
    return 0


def cmd_run(args, lists_only: bool = False) -> int:
    samples = _load(args)
    root = Path(args.outdir)
    srows = _write_summary(samples, root, getattr(args, "xlsx", False))
    write_table(root / "name_mapping.tsv",
                [{"sample": s.name, "dir_name": s.dir_name, "path": s.path, "run_dir": s.run_dir}
                 for s in samples])
    groups = _lineages(samples, args.lineage)
    if not groups:
        error("no usable samples")
        return 1
    multi = len(groups) > 1
    for lineage, group in groups.items():
        out = root / lineage if multi else root
        out.mkdir(parents=True, exist_ok=True)
        rep = dataset_report(group if multi else samples, lineage)
        ds = _dataset(args, lineage, group)
        if ds is None:
            continue
        rep.append(f"samples in matrices: {ds.n}  BUSCOs: {len(ds.genes)}  "
                   f"fragment reference: {ds.ref_len_source}"
                   + (f"  min fragment ratio: {ds.min_frag_ratio}" if ds.min_frag_ratio else ""))
        (out / "dataset_report.txt").write_text("\n".join(rep) + "\n")
        if not args.quiet:
            print("\n".join(rep), file=sys.stderr)
        qc = QCThresholds(**({k: getattr(args, f"qc_{k}") / 100 for k in
                              ("min_c", "max_d", "max_f", "max_m", "max_e")}
                             if not lists_only else {}))
        grows = analysis.gene_stats(ds)
        sstats = analysis.sample_stats(ds, qc)
        write_table(out / "gene_stats.tsv", grows, analysis.GENE_FIELDS)
        write_table(out / "sample_stats.tsv", sstats, analysis.SAMPLE_FIELDS)
        _write_lists(args, ds, out)
        if lists_only:
            continue
        _write_matrices(args, ds, out)
        if args.xlsx:
            write_xlsx(out / "busco2summary.xlsx",
                       {"summary": [r for r in srows if r["lineage"] == lineage or not multi],
                        "sample_stats": sstats, "gene_stats": grows})
        th = _thresholds(args.thresholds)
        if not args.no_plots:
            from . import plots
            plots.all_plots(ds, sstats, out / "plots", th, args.plot_format,
                            heatmap=len(ds.genes) * ds.n <= 5_000_000)
        if not args.no_html:
            from .report import write_report
            p = write_report(ds, sstats, grows, rep, out / "busco2summary_report.html")
            info(f"HTML report: {p}")
        if args.extract is not None:
            wanted = args.extract or ["S_occ050", "SF_occ050"]
            _, lists = analysis.occupancy_lists(ds, _thresholds(
                [float(w.split("occ")[1]) for w in wanted if "occ" in w] or [50]),
                sorted({w.split("_")[0] for w in wanted}))
            for w in wanted:
                if w not in lists:
                    warn(f"--extract {w}: unknown list (use e.g. SF_occ090)")
                    continue
                _do_extract(args, ds, lists[w], w.split("_")[0], out / "sequences" / w)
    params = {k: v for k, v in vars(args).items() if k != "func"}
    (root / "run_parameters.json").write_text(json.dumps(params, indent=2, default=str))
    info("done")
    return 0


def cmd_matrix(args) -> int:
    samples = _load(args)
    if args.lineage == "all":
        error("matrix: one lineage at a time (give --lineage NAME)")
        return 2
    groups = _lineages(samples, args.lineage)
    if not groups:
        error("no usable samples")
        return 1
    lineage, group = next(iter(groups.items()))
    ds = _dataset(args, lineage, group)
    if ds is None:
        return 1
    genes = ds.genes
    if args.genes:
        want = _read_set(args.genes) or set()
        unknown = want - set(ds.rows)
        if unknown:
            warn(f"--genes: {len(unknown)} IDs not in {lineage} (e.g. {sorted(unknown)[0]})")
        genes = [g for g in genes if g in want]
    k = max(args.min_taxa, analysis.need(_thresholds([args.occupancy])[0], ds.n) if args.occupancy else 0)
    if k:
        genes = [g for g in genes if ds.present_count(g, args.set) >= k]
    info(f"{lineage}: {ds.n} samples x {len(genes)} BUSCOs"
         + (f" (present as {args.set} in >= {k} samples)" if k else ""))
    masks = {g: ds.present_mask(g, args.set) for g in genes} if args.format in ("binary", "ids") else {}

    def cell(i, j):
        g = genes[i]
        if args.format == "status":
            return ds.rows[g][j]
        if args.format == "binary":
            return 1 if masks[g][j] else 0
        return (ds.samples[j].calls.get(g) or [0, 0])[1]

    import csv
    import gzip
    if args.output == "-":
        fh = sys.stdout
    else:
        Path(args.output).parent.mkdir(parents=True, exist_ok=True)
        fh = (gzip.open if args.output.endswith(".gz") else open)(args.output, "wt", newline="")
    try:
        w = csv.writer(fh, delimiter="\t", lineterminator="\n")
        names = ds.names
        if args.format == "ids":
            w.writerow(["sample", "busco_id", "status"])
            for j, name in enumerate(names):
                for g in genes:
                    if masks[g][j]:
                        w.writerow([name, g, ds.rows[g][j]])
        elif args.samples_as_rows:
            w.writerow(["sample"] + genes)
            for j, name in enumerate(names):
                w.writerow([name] + [cell(i, j) for i in range(len(genes))])
        else:
            w.writerow(["busco_id"] + names)
            for i, g in enumerate(genes):
                w.writerow([g] + [cell(i, j) for j in range(len(names))])
    finally:
        if fh is not sys.stdout:
            fh.close()
            info(f"matrix -> {args.output}")
    return 0


def cmd_extract(args) -> int:
    samples = _load(args)
    root = Path(args.outdir)
    groups = _lineages(samples, args.lineage)
    for lineage, group in groups.items():
        ds = _dataset(args, lineage, group)
        if ds is None:
            continue
        if args.genes:
            genes = [l.split()[0] for l in Path(args.genes).read_text().splitlines()
                     if l.strip() and not l.startswith("#")]
            tag = Path(args.genes).stem
        else:
            t = _thresholds([args.occupancy])[0]
            _, lists = analysis.occupancy_lists(ds, [t], [args.set])
            tag = f"{args.set}_{analysis.occ_key(t)}"
            genes = lists[tag]
        info(f"{lineage}: extracting {len(genes)} BUSCOs ({tag})")
        out = (root / lineage if len(groups) > 1 else root) / tag
        _do_extract(args, ds, genes, args.set, out)
    return 0


def main(argv: Optional[List[str]] = None) -> int:
    ap = build_parser()
    args = ap.parse_args(argv)
    if not args.cmd:
        ap.print_help()
        return 0
    try:
        if args.cmd == "list":
            return cmd_list(args)
        if args.cmd == "scan":
            return cmd_scan(args)
        if args.cmd == "summary":
            return cmd_summary(args)
        if args.cmd == "run":
            return cmd_run(args)
        if args.cmd == "genelists":
            return cmd_run(args, lists_only=True)
        if args.cmd == "extract":
            return cmd_extract(args)
        if args.cmd == "matrix":
            return cmd_matrix(args)
    except BrokenPipeError:            # e.g. `b2s matrix ... | head`
        return 0
    except KeyboardInterrupt:
        error("interrupted")
        return 130
    return 0


if __name__ == "__main__":
    sys.exit(main())
