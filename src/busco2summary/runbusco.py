"""Batch BUSCO on assemblies: SPAdes / rnaSPAdes / Trinity / MEGAHIT / Flye folders or plain FASTA files.

Naming follows busco2summary_smkv2 (-bfsd --suffix_subdir): the sample is the assembly folder name
minus the assembler suffix (JLJC1_Eun1_1201_spades_rna -> JLJC1_Eun1_1201), or the FASTA file name
without extension (X_trinity.Trinity.fasta -> X). BUSCO runs as `busco -i FASTA -o NAME --out_path OUT`,
so every result is OUT/NAME; finished samples (short_summary present) are skipped on a re-run.
"""

from __future__ import annotations

import os
import shlex
import shutil
import subprocess
import time
from concurrent.futures import ThreadPoolExecutor, as_completed
from dataclasses import dataclass
from pathlib import Path
from typing import Dict, List, Optional, Sequence, Tuple

from .log import info, warn

# file inside an assembler folder -> (assembler, BUSCO mode); --fasta sets which are used and their order
ASSEMBLY_FILES = {
    "contigs.fasta": ("spades", "genome"),
    "scaffolds.fasta": ("spades", "genome"),
    "transcripts.fasta": ("rnaspades", "transcriptome"),
    "Trinity.fasta": ("trinity", "transcriptome"),
    "final.contigs.fa": ("megahit", "genome"),
    "assembly.fasta": ("flye", "genome"),
}
DEFAULT_FASTA = ["contigs.fasta", "transcripts.fasta", "Trinity.fasta", "final.contigs.fa",
                 "assembly.fasta", "scaffolds.fasta"]
# assembler run folders: never searched below (their k-mer sub-folders hold *.fasta that are not assemblies)
ASSEMBLER_MARKERS = ("spades.log", "params.txt", "Trinity.timing", "options.json", "flye.log")
DEFAULT_SUFFIXES = ["_spades_rna", "_spades_wgs", "_rnaspades", "_metaspades", "_spades",
                    "_trinity", "_megahit", "_flye", "_assembly"]
FASTA_EXT = (".fasta", ".fa", ".fna", ".fas", ".fsa", ".ffn", ".faa", ".pep")
PROTEIN_EXT = (".faa", ".pep")
SKIP = {"busco_downloads", "logs", "tmp", "__pycache__", ".git"}


@dataclass
class Assembly:
    name: str
    fasta: str
    mode: str
    assembler: str
    status: str = "todo"            # todo / done / partial
    raw: str = ""                   # folder/file name before the suffix was removed


def strip_suffix(stem: str, suffixes: Sequence[str]) -> str:
    for s in suffixes:
        if s and stem.endswith(s) and len(stem) > len(s):
            return stem[:-len(s)]
    return stem


def file_sample(path: Path) -> Tuple[str, str, str]:
    """(stem, assembler, mode) for a FASTA file given directly."""
    stem, low = path.name, path.name.lower()
    for e in FASTA_EXT:
        if low.endswith(e):
            stem = stem[:-len(e)]
            break
    assembler = "fasta"
    if stem.lower().endswith(".trinity"):                    # Trinity >= 2.10: <out_dir>.Trinity.fasta
        stem, assembler = stem[:-len(".trinity")], "trinity"
    mode = guess_mode(path.name, assembler)
    return stem, assembler, mode


def guess_mode(name: str, assembler: str = "") -> str:
    n = name.lower()
    if n.endswith(PROTEIN_EXT) or any(x in n for x in ("_protein", "_pep", ".protein")):
        return "proteins"
    if assembler in ("trinity", "rnaspades") or any(
            x in n for x in ("transcript", "trinity", "_rna", ".rna", "cdna", ".cds", "_cds", "mrna")):
        return "transcriptome"
    return "genome"


def find_assemblies(inputs: Sequence[str], fasta_names: Sequence[str], suffixes: Sequence[str],
                    max_depth: int = 3, skip: Sequence[Path] = ()) -> Tuple[List[Assembly], List[str]]:
    """Assembler folders (by fasta_names) and loose FASTA files under each input. Returns (found, notes)."""
    found: List[Assembly] = []
    notes: List[str] = []
    skip_real = {os.path.realpath(p) for p in skip}

    def add_file(p: Path) -> None:
        stem, asm, mode = file_sample(p)
        found.append(Assembly(strip_suffix(stem, suffixes), str(p.resolve()), mode, asm, raw=stem))

    def walk(d: Path, depth: int) -> None:
        try:
            entries = sorted(d.iterdir())
        except OSError as e:
            notes.append(f"cannot read {d}: {e}")
            return
        for e in entries:
            if e.is_dir():
                if e.name in SKIP or e.name.startswith(".") or os.path.realpath(e) in skip_real:
                    continue
                hit = next((f for f in fasta_names if (e / f).is_file()), None)
                if hit:
                    asm, mode = ASSEMBLY_FILES.get(hit, ("fasta", guess_mode(hit)))
                    if mode == "genome":
                        mode = guess_mode(e.name)           # e.g. a *_rna folder holding contigs.fasta
                    found.append(Assembly(strip_suffix(e.name, suffixes), str((e / hit).resolve()), mode, asm,
                                          raw=e.name))
                elif any((e / m).exists() for m in ASSEMBLER_MARKERS + tuple(ASSEMBLY_FILES)):
                    notes.append(f"no {' / '.join(fasta_names)} in assembler folder {e} (failed run?)")
                elif any(e.glob("short_summary*.txt")) or any(e.glob("run_*_odb*")):
                    continue                                 # a BUSCO result, not an assembly
                elif depth < max_depth:
                    walk(e, depth + 1)
            elif e.is_file() and e.name.lower().endswith(FASTA_EXT):
                add_file(e)
            elif e.is_file() and e.name.lower().endswith(tuple(x + ".gz" for x in FASTA_EXT)):
                notes.append(f"skipped gzipped {e.name} (gunzip it for BUSCO)")

    for i in inputs:
        p = Path(i)
        if p.is_file():
            add_file(p)
        elif p.is_dir():
            hit = next((f for f in fasta_names if (p / f).is_file()), None)
            if hit:                                          # -i points at one assembler folder
                asm, mode = ASSEMBLY_FILES.get(hit, ("fasta", guess_mode(hit)))
                found.append(Assembly(strip_suffix(p.resolve().name, suffixes),
                                      str((p / hit).resolve()), mode, asm, raw=p.resolve().name))
            else:
                walk(p, 1)
        else:
            notes.append(f"input not found: {i}")
    return found, notes


def read_list(path: str, suffixes: Sequence[str]) -> List[Assembly]:
    """One assembly per line: FASTA | name<TAB>FASTA | name<TAB>mode<TAB>FASTA (busco_inputs.tsv)."""
    out = []
    for ln in Path(path).read_text().splitlines():
        if not ln.strip() or ln.startswith("#"):
            continue
        f = [x.strip() for x in ln.split("\t")]
        fasta = Path(f[-1] if len(f) < 4 else f[3])
        stem, asm, mode = file_sample(fasta)
        name = f[0] if len(f) >= 2 else strip_suffix(stem, suffixes)
        if len(f) >= 3 and f[1] in ("genome", "transcriptome", "proteins"):
            mode = f[1]
        if len(f) >= 4:
            asm = f[2]
        out.append(Assembly(name, str(fasta.resolve()), mode, asm))
    return out


def is_done(out: Path, name: str) -> bool:
    return any((out / name).glob("short_summary*.txt"))


def resolve_busco(busco: str = "", env: str = "") -> Tuple[List[str], Dict[str, str]]:
    """BUSCO command + environment. The env's bin/ goes first on PATH, so `#!/usr/bin/env python3`
    and miniprot / hmmsearch / metaeuk resolve without `conda activate`."""
    exe = None
    if busco:
        exe = shutil.which(busco) or (busco if os.access(busco, os.X_OK) else None)
        if not exe:
            raise SystemExit(f"--busco {busco}: not found or not executable")
    elif env:
        for base in _conda_bases():
            cand = Path(env) / "bin" / "busco" if "/" in env else base / "envs" / env / "bin" / "busco"
            if cand.is_file():
                exe = str(cand)
                break
        if not exe:
            raise SystemExit(f"--busco-env {env}: no bin/busco found (conda envs: "
                             + ", ".join(str(b / 'envs') for b in _conda_bases()) + ")")
    else:
        exe = shutil.which("busco")
        if not exe:
            for base in _conda_bases():
                cand = base / "envs" / "busco" / "bin" / "busco"
                if cand.is_file():
                    exe = str(cand)
                    break
        if not exe:
            raise SystemExit("BUSCO not found: put it on PATH, or give --busco-env ENV / --busco PATH")
    env_vars = dict(os.environ)
    env_vars["PATH"] = str(Path(exe).resolve().parent) + os.pathsep + env_vars.get("PATH", "")
    return [exe], env_vars


def _conda_bases() -> List[Path]:
    bases = []
    for v in (os.environ.get("CONDA_EXE"), os.environ.get("MAMBA_EXE")):
        if v:
            bases.append(Path(v).resolve().parent.parent)
    home = Path.home()
    bases += [home / d for d in ("miniforge3", "mambaforge", "miniconda3", "anaconda3")]
    seen, out = set(), []
    for b in bases:
        if b.is_dir() and b not in seen:
            seen.add(b)
            out.append(b)
    return out


def busco_version(cmd: List[str], env: Dict[str, str]) -> str:
    try:
        r = subprocess.run(cmd + ["--version"], capture_output=True, text=True, env=env, timeout=120)
    except (OSError, subprocess.TimeoutExpired) as e:
        return f"ERROR {e}"
    txt = (r.stdout + r.stderr).strip().splitlines()
    return txt[-1] if r.returncode == 0 and txt else "ERROR " + " | ".join(txt[-2:])


def build_cmd(a: Assembly, out: Path, base: List[str], lineage: str, cpu: int,
              extra: Sequence[str], restart: bool) -> List[str]:
    cmd = base + ["-i", a.fasta, "-o", a.name, "--out_path", str(out), "-m", a.mode, "-c", str(cpu)]
    if lineage:
        if Path(lineage).is_dir():
            cmd += ["-l", str(Path(lineage).resolve())]
            if "--offline" not in extra:
                cmd.append("--offline")
        else:
            cmd += ["-l", lineage, "--download_path", str(out / "busco_downloads")]
    if a.status == "partial" and not {"-f", "--force", "-r", "--restart"} & set(extra):
        cmd.append("--restart" if restart else "-f")
    return cmd + list(extra)


def run_all(todo: List[Assembly], out: Path, base: List[str], env: Dict[str, str], lineage: str,
            cpu: int, jobs: int, extra: Sequence[str], restart: bool) -> List[Tuple[Assembly, int, str]]:
    """Run BUSCO for every assembly, `jobs` at a time; stdout/stderr -> OUT/logs/NAME.busco.log."""
    logs = out / "logs"
    logs.mkdir(parents=True, exist_ok=True)

    def one(a: Assembly) -> Tuple[Assembly, int, str]:
        cmd = build_cmd(a, out, base, lineage, cpu, extra, restart)
        lp = logs / f"{a.name}.busco.log"
        t0 = time.time()
        with open(lp, "w") as fh:
            fh.write("# " + " ".join(shlex.quote(c) for c in cmd) + "\n")
            fh.flush()
            try:
                rc = subprocess.run(cmd, stdout=fh, stderr=subprocess.STDOUT, env=env, cwd=out).returncode
            except OSError as e:
                fh.write(f"cannot start BUSCO: {e}\n")
                rc = -1
        if rc == 0 and not is_done(out, a.name):
            rc = -2                                          # exit 0 but no short_summary
        return a, rc, f"{time.time() - t0:.0f}s"

    results = []
    ex = ThreadPoolExecutor(max_workers=max(1, jobs))
    futs = [ex.submit(one, a) for a in todo]
    try:
        for k, f in enumerate(as_completed(futs), 1):
            a, rc, dt = f.result()
            (info if rc == 0 else warn)(f"[{k}/{len(todo)}] {a.name}: "
                                        + ("done" if rc == 0 else f"FAILED (exit {rc}), see logs/{a.name}.busco.log")
                                        + f" {dt}")
            results.append((a, rc, dt))
    except KeyboardInterrupt:              # Ctrl+C: running BUSCOs get SIGINT too; queued ones must not start
        for f in futs:
            f.cancel()
        raise
    finally:
        ex.shutdown(wait=True)
    return results


def auto_jobs(cpu: int, n: int) -> int:
    cores = os.cpu_count() or 1
    return max(1, min(n, int(cores * 0.9) // max(1, cpu)))


def write_plan(path: Path, items: List[Assembly]) -> None:
    with open(path, "w") as fh:
        fh.write("#sample\tmode\tassembler\tfasta\tstatus\n")
        for a in items:
            fh.write(f"{a.name}\t{a.mode}\t{a.assembler}\t{a.fasta}\t{a.status}\n")


def dedupe(items: List[Assembly]) -> List[str]:
    """Two assemblies whose names clash after suffix removal (M4_spades, M4_spades_wgs -> M4) keep
    their full folder/file names instead. Returns the renames done."""
    by: Dict[str, List[Assembly]] = {}
    for a in items:
        by.setdefault(a.name, []).append(a)
    msgs = []
    for k, grp in by.items():
        if len(grp) > 1 and len({a.raw for a in grp if a.raw}) == len(grp):
            for a in grp:
                a.name = a.raw
            msgs.append(f"'{k}' would name {len(grp)} assemblies; kept full names: "
                        + ", ".join(a.raw for a in grp))
    return msgs


def plan(items: List[Assembly], out: Path) -> Optional[Dict[str, List[str]]]:
    """Mark done / partial; returns duplicate names (name -> fastas) or None."""
    by: Dict[str, List[str]] = {}
    for a in items:
        by.setdefault(a.name, []).append(a.fasta)
        if is_done(out, a.name):
            a.status = "done"
        elif (out / a.name).exists():
            a.status = "partial"
    dup = {k: v for k, v in by.items() if len(v) > 1}
    return dup or None
