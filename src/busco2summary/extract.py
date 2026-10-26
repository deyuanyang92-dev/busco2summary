"""Sequence extraction for selected BUSCO IDs.

* by-gene FASTA (one file per BUSCO, one record per sample) and/or
  by-sample FASTA (one file per sample, one record per BUSCO)
* statuses to use: any of S, D, F (e.g. "SF" adds fragmented copies)
* several copies (D, or S files holding >1 record): pick by
  ``best`` (highest full_table score, then longest), ``longest`` or ``first``
* nucleotide sequences: ``.fna`` from BUSCO when present; otherwise (miniprot
  runs) CDS can be rebuilt from the BUSCO GFF + the input genome
"""

from __future__ import annotations

import csv
import re
from collections import defaultdict
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path
from typing import Dict, List, Optional, Sequence, Set, Tuple

from . import parsers
from .loader import SEQ_DIRS, Sample, _full_table, _seq_base
from .log import info, warn

COMP = str.maketrans("ACGTRYMKBDHVNacgtrymkbdhvn", "TGCAYRKMVHDBNtgcayrkmvhdbn")


_B = "TCAG"
_AA = "FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG"
CODON = {a + b + c: _AA[16 * i + 4 * j + k]
         for i, a in enumerate(_B) for j, b in enumerate(_B) for k, c in enumerate(_B)}


def translate(nt: str) -> str:
    nt = nt.upper()
    return "".join(CODON.get(nt[i:i + 3], "X") for i in range(0, len(nt) - 2, 3))


def revcomp(s: str) -> str:
    return s.translate(COMP)[::-1]


def wrap(seq: str, width: int = 0) -> str:
    if width <= 0:
        return seq
    return "\n".join(seq[i:i + width] for i in range(0, len(seq), width))


def clean_protein(seq: str, stop: str = "strip") -> str:
    """stop: keep | strip (terminal *) | X (terminal stripped, internal -> X)."""
    if stop == "keep":
        return seq
    seq = seq.rstrip("*")
    if stop == "X":
        seq = seq.replace("*", "X")
    return seq


# --------------------------------------------------------------------------
# copy selection
# --------------------------------------------------------------------------

def _score_index(sample: Sample) -> Dict[str, Dict[Tuple[str, int, int], Tuple[float, int]]]:
    """gene -> location key -> (score, length) from the full table."""
    ft = _full_table(Path(sample.run_dir)) if sample.run_dir else None
    idx: Dict[str, Dict[Tuple[str, int, int], Tuple[float, int]]] = {}
    if not ft:
        return idx
    for gid, gc in parsers.parse_full_table(ft).items():
        for h in gc.hits:
            loc = parsers.hit_location(h)
            if loc:
                idx.setdefault(gid, {})[loc] = (h.score or 0.0, h.length or 0)
    return idx


def pick_record(records: List[Tuple[str, str]], policy: str,
                scores: Optional[Dict[Tuple[str, int, int], Tuple[float, int]]]
                ) -> Optional[Tuple[str, str]]:
    if not records:
        return None
    if len(records) == 1 or policy == "first":
        return records[0]
    if policy == "best" and scores:
        def key(r):
            loc = parsers.header_location(r[0])
            sc = scores.get(loc, (float("-inf"), 0))[0] if loc else float("-inf")
            return (sc, len(r[1].replace("-", "")))
        return max(records, key=key)
    return max(records, key=lambda r: len(r[1].replace("-", "")))


# --------------------------------------------------------------------------
# GFF -> CDS (for BUSCO runs that wrote no .fna, e.g. miniprot)
# --------------------------------------------------------------------------

def _gff_models(gff: Path) -> List[Tuple[Tuple[str, int, int], str, List[Tuple[int, int, int]]]]:
    """[(location, strand, [(start, end, phase) CDS ...])] ordered as in the file."""
    parents: Dict[str, Tuple[str, int, int, str]] = {}
    cds: Dict[str, List[Tuple[int, int, int]]] = defaultdict(list)
    order: List[str] = []
    with parsers.open_text(gff) as fh:
        for ln in fh:
            if ln.startswith("#") or not ln.strip():
                continue
            p = ln.rstrip("\n").split("\t")
            if len(p) < 9:
                continue
            attrs = dict(kv.split("=", 1) for kv in p[8].split(";") if "=" in kv)
            ftype = p[2]
            if ftype in ("mRNA", "transcript", "gene"):
                mid = attrs.get("ID") or f"{p[0]}:{p[3]}-{p[4]}"
                if mid not in parents:
                    parents[mid] = (p[0], int(p[3]), int(p[4]), p[6])
                    order.append(mid)
            elif ftype == "CDS":
                par = attrs.get("Parent") or attrs.get("ID") or f"{p[0]}:{p[3]}"
                par = par.split(",")[0]
                ph = int(p[7]) if p[7].isdigit() else 0
                cds[par].append((int(p[3]), int(p[4]), ph))
                if par not in parents:
                    parents[par] = (p[0], int(p[3]), int(p[4]), p[6])
                    order.append(par)
    out = []
    for mid in order:
        if not cds.get(mid):
            continue
        seq, lo, hi, strand = parents[mid]
        parts = cds[mid]
        lo = min(lo, min(a for a, _, _ in parts))
        hi = max(hi, max(b for _, b, _ in parts))
        out.append(((seq, lo, hi), strand, parts))
    return out


def _closest_model(models, loc):
    if not models:
        return None
    if loc is None:
        return models[0]
    best, bd = None, None
    for m in models:
        (s, lo, hi), _, _ = m
        if s != loc[0]:
            continue
        d = abs(lo - loc[1]) + abs(hi - loc[2])
        if bd is None or d < bd:
            best, bd = m, d
    return best


def load_contigs(genome: Path, wanted: Set[str]) -> Dict[str, str]:
    """Stream a (gz) FASTA keeping only ``wanted`` contigs."""
    out: Dict[str, str] = {}
    keep = False
    name = None
    buf: List[str] = []
    with parsers.open_text(genome) as fh:
        for ln in fh:
            if ln.startswith(">"):
                if keep and name is not None:
                    out[name] = "".join(buf)
                if len(out) == len(wanted):
                    return out          # everything found: stop reading early
                name = ln[1:].split()[0]
                keep = name in wanted
                buf = []
            elif keep:
                buf.append(ln.strip())
    if keep and name is not None:
        out[name] = "".join(buf)
    return out


def cds_from_model(contigs: Dict[str, str], model) -> Optional[str]:
    """Phase-aware CDS assembly in transcript orientation.

    Normal splice junctions keep the frame (phase == bases needed to complete
    the pending codon).  When the phase disagrees (miniprot frameshift), the
    incomplete codon is dropped and ``phase`` bases are skipped, which keeps
    the CDS in frame with the protein BUSCO reports.
    """
    (seqid, _, _), strand, parts = model
    ctg = contigs.get(seqid)
    if ctg is None:
        return None
    parts = sorted(parts, reverse=(strand == "-"))
    out = ""
    for a, b, ph in parts:
        seg = ctg[a - 1:b]
        if strand == "-":
            seg = revcomp(seg)
        pending = len(out) % 3
        expected = (3 - pending) % 3
        if ph != expected:
            out = out[:len(out) - pending] if pending else out
            seg = seg[ph:]
        out += seg
    return out[:len(out) - len(out) % 3]


_CIG = re.compile(r"(\d+)([MIDFGNUV])")


def miniprot_paf(run_dir: Path, wanted: Set[str]) -> Dict[Tuple[str, str, int, int], Tuple[str, str]]:
    """``##PAF`` records of BUSCO's miniprot GFF for the wanted query names.

    key (query, contig, tstart0, tend) -> (strand, cigar).
    """
    out: Dict[Tuple[str, str, int, int], Tuple[str, str]] = {}
    d = Path(run_dir) / "miniprot_output"
    if not d.is_dir():
        return out
    for gff in sorted(d.glob("*.gff")):
        try:
            with parsers.open_text(gff.resolve()) as fh:
                for ln in fh:
                    if not ln.startswith("##PAF"):
                        continue
                    p = ln.rstrip("\n").split("\t")
                    if len(p) < 10 or p[1] not in wanted:
                        continue
                    cg = next((x[5:] for x in p[13:] if x.startswith("cg:Z:")), None)
                    if cg:
                        out[(p[1], p[6], int(p[8]), int(p[9]))] = (p[5], cg)
        except OSError as e:
            warn(f"cannot read {gff}: {e}")
    return out


def cds_from_cigar(contig: str, tstart: int, tend: int, strand: str, cigar: str) -> str:
    """Exact CDS from a miniprot alignment (genome coordinates 0-based, half open).

    M/D: codons kept; I: no nucleotides; F and G (frameshifts): bases
    skipped, as in miniprot's own translation; N: intron skipped;
    U/V: phase-1/2 intron, the split codon is re-joined.  miniprot's ##STA
    translation (= BUSCO .faa before '*' removal) also covers the codon that
    follows the alignment, so that codon is appended (verified on real data).
    """
    if strand == "-":
        region = revcomp(contig[max(0, tstart - 3):tend])
    else:
        region = contig[tstart:tend + 3]
    pos, out = 0, []
    for n, op in _CIG.findall(cigar):
        n = int(n)
        if op in "MD":
            out.append(region[pos:pos + 3 * n])
            pos += 3 * n
        elif op == "F":
            pos += n
        elif op == "G":       # frameshift substitution: miniprot emits no residue
            pos += n
        elif op == "N":
            pos += n
        elif op == "U":
            out.append(region[pos:pos + 1] + region[pos + n - 2:pos + n])
            pos += n
        elif op == "V":
            out.append(region[pos:pos + 2] + region[pos + n - 1:pos + n])
            pos += n
    tail = region[pos:pos + 3]
    if len(tail) == 3:
        out.append(tail)
    return "".join(out)


STOPS = {"TAA", "TAG", "TGA"}


def clean_cds(seq: str, stop: str = "strip", protein: Optional[str] = None) -> str:
    """Stop codons in a rebuilt CDS, mirroring ``clean_protein``.

    strip: terminal stop removed; internal stops kept, except when the BUSCO
    protein of the same copy has none (BUSCO 5.x deletes every '*'), then
    they are removed too so that CDS and protein stay codon-aligned.
    X: terminal stop removed, internal stops -> NNN.  keep: unchanged.
    """
    if stop == "keep":
        return seq
    cod = [seq[i:i + 3] for i in range(0, len(seq) - len(seq) % 3, 3)]
    while cod and cod[-1].upper() in STOPS:
        cod.pop()
    if stop == "X":
        cod = ["NNN" if c.upper() in STOPS else c for c in cod]
    elif protein is not None and "*" not in protein.rstrip("*"):
        cod = [c for c in cod if c.upper() not in STOPS]
    return "".join(cod)


# --------------------------------------------------------------------------
# main extraction
# --------------------------------------------------------------------------

def _format_header(tpl: str, **kw) -> str:
    try:
        return tpl.format(**kw)
    except (KeyError, IndexError):
        return f"{kw['sample']}|{kw['gene']}"


def extract(samples: Sequence[Sample], genes: Sequence[str], out_dir: Path,
            statuses: str = "S", seqtype: str = "protein", dup_policy: str = "best",
            layout: str = "gene", header: str = "{sample}", sample_header: str = "{gene}",
            stop: str = "strip", min_len: int = 0, line_width: int = 0,
            genome_map: Optional[Dict[str, str]] = None, genome_from_summary: bool = False,
            threads: int = 4, short_f: Optional[Dict[str, set]] = None) -> List[dict]:
    """Extract sequences; returns manifest rows (also written to disk).

    ``short_f``: sample name -> BUSCO IDs whose F call fails --min-frag-ratio;
    these are not extracted, so each locus holds exactly the samples counted as
    present when the gene list was built.
    """
    short_f = short_f or {}
    out_dir.mkdir(parents=True, exist_ok=True)
    kinds = ["protein", "nucleotide"] if seqtype == "both" else [seqtype]
    gene_set = list(dict.fromkeys(genes))

    def per_sample(s: Sample):
        res = {k: {} for k in kinds}      # kind -> gene -> (header, seq, status, src)
        man = []
        base = _seq_base(Path(s.run_dir)) if s.run_dir else None
        d_for = lambda st: base / SEQ_DIRS[st]
        scores = _score_index(s) if dup_policy == "best" else {}
        need_gff: List[Tuple[str, str, Optional[Tuple[str, int, int]], Path]] = []
        for g in gene_set:
            st = s.status(g)
            if st not in statuses or st == "M":
                man.append([g, s.name, "", st, "not_selected" if st != "M" else "missing", "", ""])
                continue
            if st == "F" and g in short_f.get(s.name, ()):
                man.append([g, s.name, "", st, "F_below_min_frag_ratio", "", ""])
                continue
            if base is None:
                man.append([g, s.name, "", st, "no_sequence_folder", "", ""])
                continue
            d = base / SEQ_DIRS[st]
            chosen_loc = None
            for kind in kinds:
                ext = "faa" if kind == "protein" else "fna"
                f = d / f"{g}.{ext}"
                if not f.is_file() and ext == "fna" and (d / f"{g}.ffn").is_file():
                    f = d / f"{g}.ffn"
                if f.is_file():
                    recs = list(parsers.read_fasta(f))
                    rec = pick_record(recs, dup_policy, scores.get(g))
                    if rec is None:
                        man.append([g, s.name, kind, st, "empty_file", str(f), ""])
                        continue
                    seq = clean_protein(rec[1], stop) if kind == "protein" else rec[1].replace("*", "")
                    if len(seq) < min_len:
                        man.append([g, s.name, kind, st, "too_short", str(f), len(seq)])
                        continue
                    if kind == "protein":
                        chosen_loc = parsers.header_location(rec[0])
                    res[kind][g] = (rec[0], seq, st, str(f), len(recs))
                    man.append([g, s.name, kind, st, "ok" if len(recs) == 1 else f"ok_{dup_policy}_of_{len(recs)}",
                                str(f), len(seq)])
                elif kind == "nucleotide" and (d / f"{g}.gff").is_file():
                    need_gff.append((g, st, chosen_loc, d / f"{g}.gff"))
                else:
                    man.append([g, s.name, kind, st, "no_file", str(f), ""])
        if need_gff:
            genome = None
            if genome_map:
                genome = genome_map.get(s.name) or genome_map.get(s.dir_name)
            if not genome and genome_from_summary:
                genome = s.info.get("input_file")
            if not genome or not Path(str(genome)).is_file():
                for g, st, _, gff in need_gff:
                    man.append([g, s.name, "nucleotide", st, "no_fna_no_genome", str(gff), ""])
            else:
                qnames = {}
                for g, st, _loc, _gff in need_gff:
                    pf = d_for(st) / f"{g}.faa"
                    if pf.is_file():
                        recs = list(parsers.read_fasta(pf))
                        rec = pick_record(recs, dup_policy, scores.get(g))
                        if rec:
                            qnames[g] = (rec[0].split("|")[0], parsers.header_location(rec[0]),
                                         rec[1])
                paf = miniprot_paf(Path(s.run_dir), {q[0] for q in qnames.values()})
                models = {g: _gff_models(gff) for g, _, _, gff in need_gff}
                picks = {g: _closest_model(models[g], loc or (qnames.get(g) or (None, None, None))[1])
                         for g, _st, loc, _gff in need_gff}
                wanted = {m[0][0] for m in picks.values() if m}
                wanted |= {q[1][0] for q in qnames.values() if q[1]}
                contigs = load_contigs(Path(str(genome)), wanted)
                for g, st, _loc, gff in need_gff:
                    seq, how, hdr = None, "", ""
                    q = qnames.get(g)
                    if q and q[1] and q[1][0] in contigs:
                        key = (q[0], q[1][0], q[1][1], q[1][2])
                        if key not in paf:      # BUSCO may trim the reported start/end
                            alt = [k for k in paf if k[0] == q[0] and k[1] == q[1][0]
                                   and min(k[3], q[1][2]) > max(k[2], q[1][1])]
                            key = max(alt, key=lambda k: min(k[3], q[1][2]) - max(k[2], q[1][1]),
                                      default=None)
                        if key is not None:
                            strand, cg = paf[key]
                            seq = cds_from_cigar(contigs[key[1]], key[2], key[3], strand, cg)
                            how = "ok_miniprot_cigar"
                            hdr = f"{g}|{key[1]}:{key[2]}-{key[3]}|{strand}"
                    if not seq:
                        m = picks.get(g)
                        seq = cds_from_model(contigs, m) if m else None
                        how = "ok_gff_phase"
                        if m:
                            hdr = f"{g}|{m[0][0]}:{m[0][1]}-{m[0][2]}|{m[1]}"
                    if not seq or len(seq) < min_len:
                        man.append([g, s.name, "nucleotide", st, "gff_extract_failed", str(gff), ""])
                        continue
                    seq = clean_cds(seq, stop, q[2] if q else None)
                    if q and q[2]:
                        prot = clean_protein(q[2], stop)
                        how += "_verified" if translate(seq) == prot else "_translation_differs"
                    res["nucleotide"][g] = (hdr, seq, st, f"{gff}+{genome}", 1)
                    man.append([g, s.name, "nucleotide", st, how, str(gff), len(seq)])
        return s, res, man

    with ThreadPoolExecutor(max_workers=max(1, threads)) as ex:
        results = list(ex.map(per_sample, samples))

    manifest: List[dict] = []
    for _s, _res, man in results:
        for r in man:
            manifest.append(dict(zip(["busco_id", "sample", "seq_type", "status", "result",
                                      "source", "length"], r)))

    for kind in kinds:
        ext = "faa" if kind == "protein" else "fna"
        if layout in ("gene", "both"):
            d = out_dir / f"by_gene_{ext}"
            d.mkdir(parents=True, exist_ok=True)
            nfile = 0
            for g in gene_set:
                lines = []
                for s, res, _ in results:
                    rec = res[kind].get(g)
                    if rec:
                        h = _format_header(header, sample=s.name, gene=g, status=rec[2],
                                           dir=s.dir_name, orig=rec[0].split()[0])
                        lines.append(f">{h}\n{wrap(rec[1], line_width)}\n")
                if lines:
                    (d / f"{g}.{ext}").write_text("".join(lines))
                    nfile += 1
            info(f"{kind}: {nfile} per-gene FASTA files -> {d}")
        if layout in ("sample", "both"):
            d = out_dir / f"by_sample_{ext}"
            d.mkdir(parents=True, exist_ok=True)
            for s, res, _ in results:
                recs = res[kind]
                if not recs:
                    continue
                with open(d / f"{s.name}.{ext}", "w") as fh:
                    for g in gene_set:
                        rec = recs.get(g)
                        if rec:
                            h = _format_header(sample_header, sample=s.name, gene=g, status=rec[2],
                                               dir=s.dir_name, orig=rec[0].split()[0])
                            fh.write(f">{h}\n{wrap(rec[1], line_width)}\n")
            info(f"{kind}: per-sample FASTA files -> {d}")

    with open(out_dir / "extraction_manifest.tsv", "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=["busco_id", "sample", "seq_type", "status",
                                           "result", "source", "length"], delimiter="\t")
        w.writeheader()
        w.writerows(manifest)
    return manifest
