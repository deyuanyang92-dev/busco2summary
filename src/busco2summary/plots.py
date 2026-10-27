"""Static figures (matplotlib, optional dependency)."""

from __future__ import annotations

from pathlib import Path
from typing import List, Sequence

from .analysis import Dataset, occupancy_curve
from .log import info, warn

# validated (CVD-safe, light surface) S / D / F / M colours
COLORS = {"S": "#4A9FD8", "D": "#1F5F99", "F": "#D08F00", "M": "#C93C3C"}
LABELS = {"S": "Complete single-copy (S)", "D": "Complete duplicated (D)",
          "F": "Fragmented (F)", "M": "Missing (M)"}
SET_COLORS = {"S": "#2a78d6", "SD": "#6250d6", "SF": "#eb6834", "SDF": "#1baf7a"}
INK, MUTED, GRID = "#0b0b0b", "#52514e", "#e4e3df"


_WARNED = False


def _mpl():
    global _WARNED
    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
        plt.rcParams.update({"font.size": 9, "axes.edgecolor": MUTED, "axes.labelcolor": INK,
                             "xtick.color": MUTED, "ytick.color": MUTED,
                             "axes.spines.top": False, "axes.spines.right": False,
                             "pdf.fonttype": 42, "svg.fonttype": "none"})
        return plt
    except ImportError:
        if not _WARNED:
            warn("matplotlib not installed: skipping plots (pip install matplotlib)")
            _WARNED = True
        return None


def _save(fig, out: Path, stem: str, formats: Sequence[str]) -> List[Path]:
    paths = []
    for f in formats:
        p = out / f"{stem}.{f}"
        fig.savefig(p, dpi=200, bbox_inches="tight")
        paths.append(p)
    return paths


def busco_bars(ds: Dataset, sample_rows: List[dict], out: Path, formats=("png", "pdf"),
               sort: str = "SF") -> List[Path]:
    """BUSCO-style stacked horizontal bars, one per sample."""
    plt = _mpl()
    if plt is None or not sample_rows:
        return []
    rows = list(sample_rows)
    if sort == "SF":
        rows.sort(key=lambda r: (r["pct_SF"], r["pct_S"]))
    elif sort == "name":
        rows.sort(key=lambda r: r["sample"], reverse=True)
    n = len(rows)
    h = max(2.5, min(0.22 * n + 1.2, 120))
    fig, ax = plt.subplots(figsize=(8.5, h))
    left = [0.0] * n
    for code in "SDFM":
        vals = [100 * r[f"pct_{code}"] for r in rows]
        ax.barh(range(n), vals, left=left, color=COLORS[code], label=LABELS[code],
                height=0.78, edgecolor="white", linewidth=0.8)
        left = [a + b for a, b in zip(left, vals)]
    ax.set_yticks(range(n))
    ax.set_yticklabels([r["sample"] for r in rows], fontsize=max(4, min(8, 400 / max(n, 1))))
    ax.set_xlim(0, 100)
    ax.set_xlabel("% of BUSCOs")
    ax.xaxis.grid(True, color=GRID, linewidth=0.6)
    ax.set_axisbelow(True)
    ax.set_title(f"BUSCO assessment  ({ds.lineage}, n = {len(ds.genes)})", loc="left", color=INK)
    ax.legend(ncol=4, loc="lower left", bbox_to_anchor=(0, 1.0 + 0.25 / h), frameon=False,
              fontsize=8)
    paths = _save(fig, out, "busco_bars", formats)
    plt.close(fig)
    return paths


def occupancy_plot(ds: Dataset, out: Path, thresholds: Sequence[float],
                   formats=("png", "pdf")) -> List[Path]:
    """Genes retained vs. occupancy threshold, one line per status set."""
    plt = _mpl()
    if plt is None or ds.n == 0:
        return []
    fig, ax = plt.subplots(figsize=(7, 4.2))
    n = ds.n
    for k in ("S", "SD", "SF", "SDF"):
        curve = occupancy_curve(ds, k)
        xs = [100 * c / n for c in range(1, n + 1)]
        ys = [curve[c] for c in range(1, n + 1)]
        ax.step(xs, ys, where="post", color=SET_COLORS[k], linewidth=2, label=k)
        half = next((x for x, y in zip(xs, ys) if y <= len(ds.genes) / 2), None)
        if half is not None:
            ax.annotate(k, (half, len(ds.genes) / 2), xytext=(4, 2), textcoords="offset points",
                        color=INK, fontsize=8)
    for t in thresholds:
        ax.axvline(100 * t, color=GRID, linewidth=0.8, zorder=0)
    ax.set_xlabel("minimum % of samples containing the BUSCO (occupancy)")
    ax.set_ylabel("BUSCOs retained")
    ax.set_xlim(0, 103)
    ax.set_ylim(0, len(ds.genes) * 1.03)
    ax.yaxis.grid(True, color=GRID, linewidth=0.6)
    ax.set_axisbelow(True)
    ax.set_title(f"Gene occupancy: S vs S+F (fragments recover loci)  - {n} samples",
                 loc="left", color=INK)
    ax.legend(frameon=False, title="status set", fontsize=8, title_fontsize=8)
    paths = _save(fig, out, "occupancy_curves", formats)
    plt.close(fig)
    return paths


def status_heatmap(ds: Dataset, out: Path, sample_rows: List[dict],
                   formats=("png",)) -> List[Path]:
    plt = _mpl()
    if plt is None or not ds.genes or ds.n == 0:
        return []
    from matplotlib.colors import ListedColormap
    order_s = sorted(range(ds.n), key=lambda j: -sample_rows[j]["pct_SDF"])
    code = {"S": 0, "D": 1, "F": 2, "M": 3}
    genes = sorted(ds.genes, key=lambda g: (ds.rows[g].count("M"), ds.rows[g].count("F")))
    mat = [[code[ds.rows[g][j]] for j in order_s] for g in genes]
    w = min(max(6, 0.12 * ds.n + 2), 60)
    h = min(max(5, 0.012 * len(genes) + 2), 30)
    fig, ax = plt.subplots(figsize=(w, h))
    ax.imshow(mat, aspect="auto", interpolation="nearest",
              cmap=ListedColormap([COLORS[c] for c in "SDFM"]), vmin=0, vmax=3)
    ax.set_xticks(range(ds.n))
    ax.set_xticklabels([ds.samples[j].name for j in order_s], rotation=90,
                       fontsize=max(3, min(7, 500 / ds.n)))
    ax.set_yticks([])
    ax.set_ylabel(f"{len(genes)} BUSCOs (sorted by missing, then fragmented)")
    from matplotlib.patches import Patch
    ax.legend(handles=[Patch(color=COLORS[c], label=LABELS[c]) for c in "SDFM"],
              ncol=4, loc="lower left", bbox_to_anchor=(0, 1.0), frameon=False, fontsize=8)
    paths = _save(fig, out, "status_heatmap", formats)
    plt.close(fig)
    return paths


def fragment_plot(ds: Dataset, out: Path, formats=("png", "pdf")) -> List[Path]:
    """Distribution of fragment completeness (F length / reference length)."""
    plt = _mpl()
    if plt is None:
        return []
    vals = [r for g in ds.genes for r in ds.frag_ratio.get(g, []) if r is not None]
    if not vals:
        return []
    fig, ax = plt.subplots(figsize=(6.5, 3.8))
    ax.hist([min(v, 1.5) for v in vals], bins=30, color=COLORS["F"], edgecolor="white",
            linewidth=0.8)
    ax.set_xlabel(f"fragment length / reference length  ({ds.ref_len_source})")
    ax.set_ylabel("fragmented BUSCO calls")
    ax.yaxis.grid(True, color=GRID, linewidth=0.6)
    ax.set_axisbelow(True)
    med = sorted(vals)[len(vals) // 2]
    ax.axvline(med, color=INK, linewidth=1)
    ax.annotate(f"median {med:.2f}", (med, ax.get_ylim()[1] * 0.92), xytext=(4, 0),
                textcoords="offset points", fontsize=8, color=INK)
    ax.set_title(f"How complete are fragmented BUSCOs?  ({len(vals)} F calls)",
                 loc="left", color=INK)
    paths = _save(fig, out, "fragment_completeness", formats)
    plt.close(fig)
    return paths


def all_plots(ds: Dataset, sample_rows: List[dict], out: Path, thresholds,
              formats=("png", "pdf"), heatmap: bool = True) -> List[Path]:
    if _mpl() is None:
        return []
    out.mkdir(parents=True, exist_ok=True)
    paths = []
    paths += busco_bars(ds, sample_rows, out, formats)
    paths += occupancy_plot(ds, out, thresholds, formats)
    paths += fragment_plot(ds, out, formats)
    if heatmap:
        paths += status_heatmap(ds, out, sample_rows, ("png",))
    if paths:
        info(f"plots: {len(paths)} files -> {out}")
    return paths
