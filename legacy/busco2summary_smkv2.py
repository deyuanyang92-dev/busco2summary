#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# DESC: busco 2 summary
# CATEGORY: 未分组
# VERSION: v2.8.6 (Fixed: 20260310 - out_name extension bug + added --out_name_style)
"""
整合 BUSCO 批量运行与结果汇总脚本（Snakemake 全局资源托管版 v2.8.6）

v2.8.6 BUG 修复与功能升级
---------------------------------
 BUG-17  [严重] 修复 flat 模式下 out_name 保留了文件扩展名（如 .fna），
         导致 scan_completed_results() 无法匹配 BUSCO 实际输出目录（不带扩展名），
         每次运行都会重跑全部已完成样品。
         修复：out_name 现在始终使用去掉扩展名的文件基名。
 FEAT-18 新增参数 --out_name_style：
           no_ext   (默认) 输出目录名去掉扩展名 → GCA_xxx_genomic
           keep_ext          保留完整文件名       → GCA_xxx_genomic.fna
           stem_only         只保留第一个点之前    → GCA_000326865
         用户可根据实际需要控制输出目录命名方式。

v2.8.5 BUG 修复与功能升级清单
---------------------------------
 FEAT-15 新增自动 `--offline` 探测功能。
 BUG-16  兼容了 `--auto-lineage` 的用法。

v2.8.4 BUG 修复清单
---------------------------------
 BUG-14  处理因意外中断导致半成品目录存在时的报错。

v2.8.3 BUG 修复清单
---------------------------------
 BUG-13  修复 Snakemake resources 值必须为整数的严格限制。

推荐用法
--------
先预检查：
  python3 busco2summary -i INPUT -o OUT -l /path/to/lineage

确认后执行：
  python3 busco2summary -i INPUT -o OUT -l /path/to/lineage -y

控制输出目录命名（新）：
  python3 busco2summary -i INPUT -o OUT -l /path/to/lineage --out_name_style no_ext -y
"""

from __future__ import annotations

__version__ = "2.8.6"

import argparse
import glob
import json
import math
import os
import re
import shutil
import subprocess
import sys
import textwrap
import shlex
import traceback
from collections import Counter
from concurrent.futures import ThreadPoolExecutor, as_completed
from dataclasses import asdict, dataclass, field
from typing import Dict, List, Optional, Sequence, Tuple

ALLOWED_EXTENSIONS: List[str] = []

CPU_HARD_CAP_FRACTION = 0.90
MEM_HARD_CAP_FRACTION = 0.90

VALID_MODES = {"genome", "transcriptome", "proteins"}

# v2.8.6: 输出目录命名策略
OUT_NAME_STYLE_CHOICES = {"no_ext", "keep_ext", "stem_only"}


# ============================================================
#  系统信息采集
# ============================================================

def get_system_cores() -> int:
    try:
        import multiprocessing
        return max(1, int(multiprocessing.cpu_count()))
    except Exception:
        return 4


def get_system_mem_info_gb() -> Tuple[float, float]:
    try:
        import psutil
        vm = psutil.virtual_memory()
        return vm.total / (1024 ** 3), vm.available / (1024 ** 3)
    except Exception:
        pass
    total_gb = 32.0
    avail_gb = 24.0
    try:
        meminfo: Dict[str, int] = {}
        with open("/proc/meminfo") as f:
            for line in f:
                if ":" not in line:
                    continue
                key, rest = line.split(":", 1)
                toks = rest.strip().split()
                if toks:
                    meminfo[key] = int(toks[0])
        if "MemTotal" in meminfo:
            total_gb = meminfo["MemTotal"] / (1024 ** 2)
        if "MemAvailable" in meminfo:
            avail_gb = meminfo["MemAvailable"] / (1024 ** 2)
        else:
            avail_gb = max(1.0, total_gb * 0.75)
    except Exception:
        pass
    return total_gb, avail_gb


def sample_cpu_usage_ratio(sample_sec: float = 0.6) -> Optional[float]:
    try:
        import psutil
        val = psutil.cpu_percent(interval=max(0.1, float(sample_sec))) / 100.0
        return min(1.0, max(0.0, float(val)))
    except Exception:
        return None


def get_loadavg_ratio(total_cores: int) -> Optional[float]:
    try:
        load1 = os.getloadavg()[0]
        if total_cores <= 0:
            return None
        return min(2.0, max(0.0, float(load1) / float(total_cores)))
    except Exception:
        return None


def cpu_hard_cap_cores(total_cores: int) -> int:
    return max(1, int(math.floor(max(1, total_cores) * CPU_HARD_CAP_FRACTION)))


def mem_hard_cap_gb(total_mem_gb: float, avail_mem_gb: float) -> float:
    total_cap = max(4.0, total_mem_gb * MEM_HARD_CAP_FRACTION)
    avail_cap = max(4.0, avail_mem_gb * MEM_HARD_CAP_FRACTION)
    return max(4.0, min(total_cap, avail_cap))


# ============================================================
#  模式推断辅助函数
# ============================================================

def detect_busco_mode_from_filename(filename: str) -> str:
    fl = filename.lower()
    if any(x in fl for x in [".pep", ".faa", ".protein", ".aa", "_protein", "_pep"]):
        return "proteins"
    if any(x in fl for x in [
        ".cdna", ".cds", ".transcripts", ".rna", "_rna", "_transcriptome",
        "_cdna", "_transcripts", ".mrna", "_mrna", "transcripts.fasta",
    ]):
        return "transcriptome"
    return "genome"


def translate_mode_to_busco_param(mode_type: str) -> str:
    return {"proteins": "proteins", "transcriptome": "transcriptome", "genome": "genome"}.get(mode_type, "genome")


# ============================================================
#  v2.8.6: 输出名称转换工具
# ============================================================

def make_out_name(filename: str, style: str = "no_ext") -> str:
    """
    根据 --out_name_style 将文件名转为输出目录名。

    style:
      no_ext     去掉最后一个扩展名 → GCA_000326865.1_Helobdella_robusta_v1.0_genomic
      keep_ext   保留完整文件名     → GCA_000326865.1_Helobdella_robusta_v1.0_genomic.fna
      stem_only  只保留到第一个点   → GCA_000326865

    BUSCO 默认行为是 -o 参数作为输出子目录名，不带路径。
    BUSCO 自己不会在 out_name 后面加或去扩展名。
    """
    if style == "keep_ext":
        return filename
    elif style == "stem_only":
        # 只保留第一个 . 之前的部分
        dot_pos = filename.find(".")
        if dot_pos > 0:
            return filename[:dot_pos]
        return filename
    else:
        # no_ext (default): 去掉最后一个扩展名
        base, ext = os.path.splitext(filename)
        # 处理双扩展名情况（如 .tar.gz），但对 .fna/.fasta 等只去一层
        return base


# ============================================================
#  已完成结果检测
# ============================================================

def scan_completed_results(output_dir: str,
                           files_to_process: Sequence[Tuple[str, str]]
                           ) -> Tuple[List[Tuple[str, str]], List[Tuple[str, str]]]:
    completed: List[Tuple[str, str]] = []
    pending: List[Tuple[str, str]] = []

    if not os.path.isdir(output_dir):
        return [], list(files_to_process)

    for input_file, out_name in files_to_process:
        out_subdir = os.path.join(output_dir, out_name)
        if os.path.isdir(out_subdir) and glob.glob(os.path.join(out_subdir, "short_summary*.txt")):
            completed.append((input_file, out_name))
        else:
            pending.append((input_file, out_name))

    return completed, pending


# ============================================================
#  输入目录结构检查
# ============================================================

def check_and_suggest_suffix_subdir(input_path: str, suffix_subdir: str, debug: bool = False) -> None:
    if not os.path.isdir(input_path):
        return
    subdirs = [d for d in os.listdir(input_path) if os.path.isdir(os.path.join(input_path, d))]
    if '*' in suffix_subdir or '?' in suffix_subdir:
        print("\n警告：--suffix_subdir 不支持通配符，您输入的是: '{}'".format(suffix_subdir))
        suffixes: List[str] = []
        for subdir in subdirs:
            if '_' in subdir:
                parts = subdir.split('_')
                if len(parts) > 1:
                    suffixes.append('_' + '_'.join(parts[-2:]))
                if len(parts) > 2:
                    suffixes.append('_' + '_'.join(parts[-3:]))
        common = [s for s, c in Counter(suffixes).items() if c > 1]
        if common:
            print("建议后缀：")
            for s in common:
                exs = [d for d in subdirs if d.endswith(s)][:3]
                print("  --suffix_subdir '{}' (匹配: {})".format(s, ', '.join(exs)))
    elif suffix_subdir:
        matched = [d for d in subdirs if d.endswith(suffix_subdir)]
        if not matched:
            print("\n警告：没有以 '{}' 结尾的子目录".format(suffix_subdir))
            for d in subdirs[:10]:
                print("  " + d)


def check_fasta_structure(input_path: str, debug: bool = False) -> None:
    print("\n自动检查输入目录结构...")
    if os.path.isfile(input_path):
        mode = translate_mode_to_busco_param(detect_busco_mode_from_filename(input_path))
        print("输入是文件，推断 BUSCO 模式: -m {}".format(mode))
        return
    if not os.path.isdir(input_path):
        print("输入路径既不是文件也不是目录：{}".format(input_path))
        return

    fasta_files = [
        f for f in os.listdir(input_path)
        if os.path.isfile(os.path.join(input_path, f)) and any(f.lower().endswith(e) for e in ALLOWED_EXTENSIONS)
    ]
    subdirs = [d for d in os.listdir(input_path) if os.path.isdir(os.path.join(input_path, d))]

    if fasta_files:
        file_modes = {f: detect_busco_mode_from_filename(f) for f in fasta_files}
        mc = Counter(file_modes.values())
        best = translate_mode_to_busco_param(mc.most_common(1)[0][0])
        print("检测到 flat 模式序列文件 ({} 个)：".format(len(fasta_files)))
        for f, m in list(file_modes.items())[:15]:
            print("  {} (推断: -m {})".format(f, translate_mode_to_busco_param(m)))
        if len(fasta_files) > 15:
            print("  ... 及其他 {} 个文件".format(len(fasta_files) - 15))
        print("\n建议命令：")
        print("python3 {} -i {} -o output_dir --store_type_contigs flat -m {} -l /path/to/lineage".format(
            sys.argv[0], input_path, best
        ))

    if subdirs:
        spades_dirs = []
        fasta_file_types = set()
        dir_modes: Dict[str, str] = {}
        for d in subdirs:
            dp = os.path.join(input_path, d)
            try:
                fastas = [
                    f for f in os.listdir(dp)
                    if os.path.isfile(os.path.join(dp, f)) and any(f.lower().endswith(e) for e in ALLOWED_EXTENSIONS)
                ]
            except (PermissionError, OSError):
                continue

            if fastas:
                dm = detect_busco_mode_from_filename(d)
                fm = Counter([detect_busco_mode_from_filename(f) for f in fastas]).most_common(1)[0][0]
                mode = fm if fm != "genome" else dm
                dir_modes[d] = mode
                spades_found = [f for f in fastas if f in ('contigs.fasta', 'scaffolds.fasta', 'transcripts.fasta')]
                if spades_found:
                    spades_dirs.append(d)
                    fasta_file_types.update(spades_found)

        if dir_modes:
            mc = Counter(dir_modes.values())
            best_mode = translate_mode_to_busco_param(mc.most_common(1)[0][0] if mc else "genome")
            if spades_dirs:
                best_fasta = (
                    "transcripts.fasta"
                    if best_mode == "transcriptome" and "transcripts.fasta" in fasta_file_types
                    else ("scaffolds.fasta" if "scaffolds.fasta" in fasta_file_types else "contigs.fasta")
                )
                print("\n检测到 SPAdes 结构 ({} 个匹配子目录)，建议命令：".format(len(spades_dirs)))
                print("python3 {} -i {} -bfsd --which_spades_fasta {} -o output_dir -m {} -l /path/to/lineage".format(
                    sys.argv[0], input_path, best_fasta, best_mode
                ))
            else:
                print("\n建议 subdir 模式命令 ({} 个匹配子目录)：".format(len(dir_modes)))
                print("python3 {} -i {} -o output_dir --store_type_contigs subdir -m {} -l /path/to/lineage".format(
                    sys.argv[0], input_path, best_mode
                ))

    if not fasta_files and not subdirs:
        print("输入目录中未检测到任何文件或子目录，请检查路径。")


# ============================================================
#  文件收集函数
# ============================================================

def _print_files(files_to_process: Sequence[Tuple[str, str]], suggested_modes: Dict[str, str],
                 max_display: int = 20) -> None:
    total = len(files_to_process)
    print("\n共发现 {} 个待处理文件：".format(total))
    for i, (fp, out_name) in enumerate(files_to_process):
        if i >= max_display:
            print("  ... 及其他 {} 个文件（共 {} 个）".format(total - max_display, total))
            break
        mode_info = ""
        if out_name in suggested_modes:
            mode_info = "  (推荐: -m {})".format(translate_mode_to_busco_param(suggested_modes[out_name]))
        print("  {} -> {}{}".format(fp, out_name, mode_info))
    if suggested_modes:
        mc = Counter(suggested_modes.values())
        best = translate_mode_to_busco_param(mc.most_common(1)[0][0])
        print("\n整体推荐 BUSCO 模式: -m {}".format(best))


def gather_contigs_files(input_path: str, store_type: str, suffix_fasta: str, contigs_list_file: str,
                         out_name_style: str = "no_ext",
                         suggest_mode: bool = False, debug: bool = False
                         ) -> Tuple[List[Tuple[str, str]], Optional[Dict[str, str]], Dict[str, object]]:
    """
    收集输入文件。

    v2.8.6 修复：flat 模式下 out_name 现在通过 make_out_name() 根据 --out_name_style 统一处理，
    默认去掉扩展名，与 BUSCO 输出目录名一致，从而让 scan_completed_results() 能正确识别已完成样品。
    """
    files_to_process: List[Tuple[str, str]] = []
    selected_names = set()
    missing_selected: List[str] = []
    suggested_modes: Dict[str, str] = {}

    if contigs_list_file:
        try:
            with open(contigs_list_file) as f:
                for line in f:
                    s = line.strip()
                    if s:
                        selected_names.add(s)
        except Exception as e:
            print("错误：读取 contigs_lists 文件失败: {}".format(e))
            if debug:
                traceback.print_exc()
            sys.exit(1)

    matched_selected = set()

    if os.path.isfile(input_path):
        fn = os.path.basename(input_path)
        out_name = make_out_name(fn, style=out_name_style)
        if suggest_mode:
            suggested_modes[out_name] = detect_busco_mode_from_filename(fn)
        files_to_process.append((os.path.abspath(input_path), out_name))

    elif os.path.isdir(input_path):
        if store_type == "subdir":
            for entry in sorted(os.listdir(input_path)):
                subdir_path = os.path.join(input_path, entry)
                if not os.path.isdir(subdir_path):
                    continue
                if selected_names and entry not in selected_names:
                    continue
                fasta_found = None
                for file in sorted(os.listdir(subdir_path)):
                    fp = os.path.join(subdir_path, file)
                    if os.path.isfile(fp) and any(file.lower().endswith(e) for e in ALLOWED_EXTENSIONS):
                        base = os.path.splitext(file)[0]
                        if suffix_fasta and not base.endswith(suffix_fasta):
                            continue
                        fasta_found = fp
                        break
                if fasta_found:
                    # subdir 模式：out_name 就是子目录名本身（不涉及扩展名问题）
                    out_name = entry
                    matched_selected.add(entry)
                    if suggest_mode:
                        dm = detect_busco_mode_from_filename(entry)
                        fm = detect_busco_mode_from_filename(os.path.basename(fasta_found))
                        suggested_modes[out_name] = fm if fm != "genome" else dm
                    files_to_process.append((os.path.abspath(fasta_found), out_name))
                else:
                    if debug:
                        print("警告：子目录 {} 中未找到序列文件。".format(subdir_path))

        elif store_type == "flat":
            for file in sorted(os.listdir(input_path)):
                fp = os.path.join(input_path, file)
                if not os.path.isfile(fp):
                    continue
                if not any(file.lower().endswith(e) for e in ALLOWED_EXTENSIONS):
                    continue
                base = os.path.splitext(file)[0]
                if suffix_fasta:
                    if not base.endswith(suffix_fasta):
                        continue
                    base = base[:-len(suffix_fasta)]
                if selected_names and file not in selected_names and base not in selected_names:
                    continue

                # ══════════════════════════════════════════════════════════
                # v2.8.6 BUG-17 修复核心：
                #   旧代码：out_name = base if suffix_fasta else file
                #     → suffix_fasta 为空时 out_name 含扩展名 (.fna)
                #     → scan_completed_results 找 output/xxx.fna/ 目录
                #     → BUSCO 实际输出 output/xxx/ 目录 → 永远不匹配 → 永远重跑
                #
                #   新代码：统一用 make_out_name() 根据 --out_name_style 决定
                #     默认 no_ext → 去掉扩展名 → 与 BUSCO 输出目录名一致
                # ══════════════════════════════════════════════════════════
                if suffix_fasta:
                    out_name = base  # suffix_fasta 模式下 base 已经去掉了后缀
                else:
                    out_name = make_out_name(file, style=out_name_style)

                matched_selected.add(file)
                matched_selected.add(base)
                if suggest_mode:
                    suggested_modes[out_name] = detect_busco_mode_from_filename(file)
                files_to_process.append((os.path.abspath(fp), out_name))
        else:
            sys.exit("错误：未知的 store_type_contigs: {}".format(store_type))
    else:
        sys.exit("错误：输入 {} 既不是文件也不是目录。".format(input_path))

    if selected_names:
        missing_selected = sorted([x for x in selected_names if x not in matched_selected])

    if not files_to_process:
        sys.exit("错误：未找到任何符合条件的输入文件。\n"
                 "  请检查：\n"
                 "    - 输入路径 -i 是否正确\n"
                 "    - --extension_files 是否包含您的文件后缀（当前: {}）\n"
                 "    - --store_type_contigs 是否与目录结构匹配（当前: {}）\n"
                 "    - --suffix_fasta 是否过滤过严（当前: '{}')".format(
                     ','.join(ALLOWED_EXTENSIONS), store_type, suffix_fasta))

    _print_files(files_to_process, suggested_modes if suggest_mode else {})
    extra = {
        "missing_selected": missing_selected,
        "selected_requested": sorted(selected_names),
        "selected_matched_count": len(matched_selected),
        "out_name_style": out_name_style,
    }
    return files_to_process, suggested_modes if suggest_mode else None, extra


def gather_spades_busco_files(input_path: str, suffix_subdir: str, which_spades_fasta: str,
                              output_dir: str, out_name_style: str = "no_ext",
                              suggest_mode: bool = False, debug: bool = False
                              ) -> Tuple[List[Tuple[str, str]], Optional[str], Optional[Dict[str, str]], Dict[str, object]]:
    files_to_process: List[Tuple[str, str]] = []
    suggested_modes: Dict[str, str] = {}
    temp_dir = os.path.join(output_dir, "tmp_busco_inputs_spades")
    os.makedirs(temp_dir, exist_ok=True)
    missing_expected_fasta: List[str] = []

    if not os.path.isdir(input_path):
        sys.exit("错误：输入 {} 不是目录。".format(input_path))
    if suffix_subdir:
        check_and_suggest_suffix_subdir(input_path, suffix_subdir, debug)

    default_fasta_mode = "transcriptome" if which_spades_fasta == "transcripts.fasta" else "genome"

    for entry in sorted(os.listdir(input_path)):
        subdir_path = os.path.join(input_path, entry)
        if not os.path.isdir(subdir_path):
            continue
        if suffix_subdir:
            if not entry.endswith(suffix_subdir):
                continue
            base = entry[:-len(suffix_subdir)]
        else:
            base = entry

        fasta_path = os.path.join(subdir_path, which_spades_fasta)
        if os.path.isfile(fasta_path):
            temp_fasta = os.path.join(temp_dir, base + ".fasta")
            try:
                shutil.copy(fasta_path, temp_fasta)
            except Exception as e:
                print("错误：复制文件 {} 失败: {}".format(fasta_path, e))
                if debug:
                    traceback.print_exc()
                continue
            if suggest_mode:
                dm = detect_busco_mode_from_filename(entry)
                mode = dm if dm != "genome" else default_fasta_mode
                suggested_modes[base] = mode
            files_to_process.append((os.path.abspath(temp_fasta), base))
        else:
            print("警告：{} 中未找到 {}。".format(subdir_path, which_spades_fasta))
            missing_expected_fasta.append(subdir_path)

    if not files_to_process:
        sys.exit("错误：未找到任何符合条件的 SPAdes BUSCO 输入文件。")

    _print_files(files_to_process, suggested_modes if suggest_mode else {})
    extra = {
        "missing_expected_fasta_dirs": missing_expected_fasta,
        "temp_dir": os.path.abspath(temp_dir),
    }
    return files_to_process, temp_dir, suggested_modes if suggest_mode else None, extra


def gather_files_from_list(file_list_path: str, input_dir: Optional[str] = None,
                           out_name_style: str = "no_ext",
                           suggest_mode: bool = False, debug: bool = False
                           ) -> Tuple[List[Tuple[str, str]], Optional[Dict[str, str]], Dict[str, object]]:
    files_to_process: List[Tuple[str, str]] = []
    suggested_modes: Dict[str, str] = {}
    missing_list_entries: List[str] = []

    if not os.path.isfile(file_list_path):
        sys.exit("错误：文件列表 {} 不存在。".format(file_list_path))

    with open(file_list_path) as f:
        for line in f:
            fp = line.strip()
            if not fp:
                continue
            resolved = None
            if os.path.isabs(fp):
                if os.path.isfile(fp):
                    resolved = fp
                else:
                    missing_list_entries.append(fp)
            elif input_dir:
                full = os.path.join(input_dir, fp)
                if os.path.isfile(full):
                    resolved = full
                else:
                    missing_list_entries.append(full)
            else:
                missing_list_entries.append(fp)

            if resolved:
                fn = os.path.basename(resolved)
                out_name = make_out_name(fn, style=out_name_style)
                if suggest_mode:
                    suggested_modes[out_name] = detect_busco_mode_from_filename(fn)
                files_to_process.append((os.path.abspath(resolved), out_name))

    if not files_to_process:
        sys.exit("错误：从文件列表 {} 中未找到任何有效文件。".format(file_list_path))

    _print_files(files_to_process, suggested_modes if suggest_mode else {})
    extra = {"missing_list_entries": missing_list_entries}
    return files_to_process, suggested_modes if suggest_mode else None, extra


# ============================================================
#  BUSCO 运行函数
# ============================================================

def run_busco(input_file: str, output_name: str, lineage: str, mode: str, cpu: int,
              busco_extra_args: Sequence[str], work_dir: str, busco_cmd_base: Sequence[str],
              resume_strategy: str = "force", debug: bool = False) -> int:
    cmd = list(busco_cmd_base) + [
        "-i", input_file,
        "-o", output_name,
        "-m", mode,
        "-c", str(cpu),
    ]
    if lineage:
        cmd.extend(["-l", lineage])
    out_subdir = os.path.join(work_dir, output_name)
    if os.path.exists(out_subdir):
        if not any(x in busco_extra_args for x in ["-f", "--force", "-r", "--restart"]):
            if resume_strategy == "restart":
                cmd.append("--restart")
            else:
                cmd.append("--force")
    if busco_extra_args:
        cmd.extend(busco_extra_args)

    print("\n" + "=" * 60)
    print("开始运行 BUSCO: {}".format(input_file))
    print("执行命令: " + " ".join(cmd))
    print("=" * 60)

    try:
        result = subprocess.run(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True, cwd=work_dir)
        if result.returncode != 0:
            err_detail = ""
            if result.stderr and result.stderr.strip():
                err_detail += "\n--- STDERR (最后 4000 字符) ---\n{}".format(result.stderr[-4000:])
            if result.stdout and result.stdout.strip():
                err_detail += "\n--- STDOUT (最后 4000 字符) ---\n{}".format(result.stdout[-4000:])
            if not err_detail:
                err_detail = "\n(BUSCO 无任何输出，可能是命令本身启动失败)"
            print("BUSCO 运行失败 (返回码 {}): {}{}".format(result.returncode, input_file, err_detail))
        else:
            print("BUSCO 运行成功: {}".format(input_file))
        return int(result.returncode)
    except FileNotFoundError as e:
        print("错误：BUSCO 命令找不到: {}".format(e))
        return -1
    except Exception as e:
        print("运行 BUSCO 出错: {} (类型: {})".format(e, type(e).__name__))
        if debug:
            traceback.print_exc()
        return -1


# ============================================================
#  数据类定义
# ============================================================

@dataclass
class PreflightReport:
    ok: bool = True
    errors: List[str] = field(default_factory=list)
    warnings: List[str] = field(default_factory=list)
    infos: List[str] = field(default_factory=list)
    sample_count: int = 0
    completed_count: int = 0
    pending_count: int = 0
    duplicate_output_names: List[str] = field(default_factory=list)
    zero_size_inputs: List[str] = field(default_factory=list)
    unreadable_inputs: List[str] = field(default_factory=list)
    extra_checks: Dict[str, object] = field(default_factory=dict)


@dataclass
class ResourcePlan:
    host_total_cores: int
    host_total_mem_gb: float
    host_available_mem_gb: float
    cpu_hard_cap_cores: int
    mem_hard_cap_gb: float
    cpu_usage_ratio: Optional[float]
    loadavg_ratio: Optional[float]
    mem_pressure_ratio: float
    busy_score: float
    busy_level: str
    requested_total_cores: int
    requested_total_mem_gb: float
    requested_cpu_per_task: int
    requested_mem_per_task_gb: float
    requested_max_tasks: int
    auto_total_cores_budget: int
    auto_total_mem_gb_budget: float
    resolved_total_cores: int
    resolved_total_mem_gb: float
    resolved_cpu_per_task: int
    resolved_mem_per_task_gb: float
    resolved_max_tasks: int
    max_tasks_by_cpu: int
    max_tasks_by_mem: int
    busy_parallel_cap: int
    note: str


def ensure_json_dump(path: str, obj: object) -> None:
    os.makedirs(os.path.dirname(os.path.abspath(path)), exist_ok=True)
    with open(path, "w", encoding="utf-8") as f:
        json.dump(obj, f, indent=2, ensure_ascii=False)


@dataclass
class SystemSnapshot:
    total_cores: int
    total_mem_gb: float
    avail_mem_gb: float
    cpu_hard_cap_cores: int
    mem_hard_cap_gb: float
    cpu_usage_ratio: Optional[float]
    loadavg_ratio: Optional[float]
    mem_pressure_ratio: float
    busy_score: float
    busy_level: str


# ============================================================
#  系统快照与资源规划
# ============================================================

def get_system_snapshot(sample_sec: float = 0.6) -> SystemSnapshot:
    total_cores = get_system_cores()
    total_mem_gb, avail_mem_gb = get_system_mem_info_gb()
    cpu_usage_ratio = sample_cpu_usage_ratio(sample_sec)
    loadavg_ratio = get_loadavg_ratio(total_cores)
    mem_pressure_ratio = min(1.0, max(0.0, 1.0 - (avail_mem_gb / max(1e-9, total_mem_gb))))
    score_parts = [mem_pressure_ratio]
    if cpu_usage_ratio is not None:
        score_parts.append(min(1.0, cpu_usage_ratio))
    if loadavg_ratio is not None:
        score_parts.append(min(1.0, loadavg_ratio))
    busy_score = max(score_parts) if score_parts else 0.5
    if busy_score >= 0.85:
        busy_level = "extreme"
    elif busy_score >= 0.70:
        busy_level = "high"
    elif busy_score >= 0.50:
        busy_level = "medium"
    else:
        busy_level = "low"
    return SystemSnapshot(
        total_cores=total_cores, total_mem_gb=total_mem_gb, avail_mem_gb=avail_mem_gb,
        cpu_hard_cap_cores=cpu_hard_cap_cores(total_cores),
        mem_hard_cap_gb=mem_hard_cap_gb(total_mem_gb, avail_mem_gb),
        cpu_usage_ratio=cpu_usage_ratio, loadavg_ratio=loadavg_ratio,
        mem_pressure_ratio=mem_pressure_ratio, busy_score=busy_score, busy_level=busy_level,
    )


def compute_resource_plan(args: argparse.Namespace, num_tasks: int, snapshot: SystemSnapshot) -> ResourcePlan:
    total_cores = snapshot.total_cores
    avail_mem_gb = snapshot.avail_mem_gb
    busy = snapshot.busy_score if args.dynamic_busy_guard == "yes" else 0.20

    cpu_usage = snapshot.cpu_usage_ratio if snapshot.cpu_usage_ratio is not None else min(1.0, busy)
    load_ratio = snapshot.loadavg_ratio if snapshot.loadavg_ratio is not None else min(1.0, busy)
    cpu_busy = min(1.0, max(cpu_usage, min(1.0, load_ratio)))

    free_cpu_equiv = max(1.0, total_cores * max(0.08, 1.0 - cpu_busy))
    auto_total_cores_budget = max(
        1, int(math.floor(min(snapshot.cpu_hard_cap_cores, free_cpu_equiv * max(0.50, 0.92 - 0.30 * busy))))
    )
    mem_factor = max(0.35, 0.82 - 0.30 * snapshot.mem_pressure_ratio)
    auto_total_mem_gb_budget = max(4.0, min(snapshot.mem_hard_cap_gb, avail_mem_gb * mem_factor))

    requested_total_cores = int(args.total_cores) if args.total_cores > 0 else 0
    requested_total_mem_gb = float(args.total_mem_gb) if args.total_mem_gb > 0 else 0.0

    resolved_total_cores = requested_total_cores if requested_total_cores > 0 else auto_total_cores_budget
    resolved_total_mem_gb = requested_total_mem_gb if requested_total_mem_gb > 0 else auto_total_mem_gb_budget
    resolved_total_cores = max(1, min(snapshot.cpu_hard_cap_cores, resolved_total_cores))
    resolved_total_mem_gb = max(4.0, min(snapshot.mem_hard_cap_gb, resolved_total_mem_gb))

    if args.dynamic_busy_guard == "yes":
        resolved_total_cores = min(resolved_total_cores, max(1, auto_total_cores_budget))
        resolved_total_mem_gb = min(resolved_total_mem_gb, max(4.0, auto_total_mem_gb_budget))

    if busy >= 0.85:       busy_parallel_cap = 1
    elif busy >= 0.70:     busy_parallel_cap = 2
    elif busy >= 0.50:     busy_parallel_cap = 3
    elif busy >= 0.35:     busy_parallel_cap = 4
    else:                  busy_parallel_cap = max(1, int(args.auto_max_jobs))

    requested_max_tasks = max(1, int(args.max_tasks) if args.max_tasks > 0 else min(num_tasks, int(args.auto_max_jobs)))

    requested_cpu_per_task = int(args.cpu) if args.cpu > 0 else 0
    if requested_cpu_per_task > 0:
        resolved_cpu_per_task = requested_cpu_per_task
    else:
        if num_tasks <= 1:
            resolved_cpu_per_task = max(1, min(32, resolved_total_cores))
        else:
            preferred_parallel = min(num_tasks, busy_parallel_cap)
            resolved_cpu_per_task = max(1, min(16, resolved_total_cores // max(1, preferred_parallel)))
            if resolved_total_cores >= 8:
                resolved_cpu_per_task = max(2, resolved_cpu_per_task)

    requested_mem_per_task_gb = float(args.mem_per_task_gb) if args.mem_per_task_gb > 0 else 0.0
    if requested_mem_per_task_gb > 0:
        resolved_mem_per_task_gb = requested_mem_per_task_gb
    else:
        preferred_parallel_mem = min(num_tasks, busy_parallel_cap)
        resolved_mem_per_task_gb = resolved_total_mem_gb / max(1, preferred_parallel_mem)
        resolved_mem_per_task_gb = max(4.0, min(64.0, resolved_mem_per_task_gb))

    resolved_cpu_per_task = max(1, min(resolved_total_cores, resolved_cpu_per_task))
    resolved_mem_per_task_gb = max(4.0, min(resolved_total_mem_gb, resolved_mem_per_task_gb))

    max_tasks_by_cpu = max(1, resolved_total_cores // max(1, resolved_cpu_per_task))
    max_tasks_by_mem = max(1, int(resolved_total_mem_gb // max(0.1, resolved_mem_per_task_gb)))
    resolved_max_tasks = max(1, min(num_tasks, requested_max_tasks, max_tasks_by_cpu, max_tasks_by_mem, busy_parallel_cap))

    resolved_cpu_per_task = max(1, min(resolved_cpu_per_task, resolved_total_cores // max(1, resolved_max_tasks)))
    max_tasks_by_cpu = max(1, resolved_total_cores // max(1, resolved_cpu_per_task))
    resolved_max_tasks = max(1, min(resolved_max_tasks, max_tasks_by_cpu))

    note = ("动态保守限流已启用" if args.dynamic_busy_guard == "yes" else "动态保守限流已关闭") + "，强制施加 90% CPU/内存硬上限。"

    return ResourcePlan(
        host_total_cores=total_cores, host_total_mem_gb=round(snapshot.total_mem_gb, 3),
        host_available_mem_gb=round(avail_mem_gb, 3), cpu_hard_cap_cores=snapshot.cpu_hard_cap_cores,
        mem_hard_cap_gb=round(snapshot.mem_hard_cap_gb, 3),
        cpu_usage_ratio=None if snapshot.cpu_usage_ratio is None else round(snapshot.cpu_usage_ratio, 4),
        loadavg_ratio=None if snapshot.loadavg_ratio is None else round(snapshot.loadavg_ratio, 4),
        mem_pressure_ratio=round(snapshot.mem_pressure_ratio, 4),
        busy_score=round(snapshot.busy_score, 4), busy_level=snapshot.busy_level,
        requested_total_cores=requested_total_cores, requested_total_mem_gb=round(requested_total_mem_gb, 3),
        requested_cpu_per_task=requested_cpu_per_task, requested_mem_per_task_gb=round(requested_mem_per_task_gb, 3),
        requested_max_tasks=requested_max_tasks, auto_total_cores_budget=auto_total_cores_budget,
        auto_total_mem_gb_budget=round(auto_total_mem_gb_budget, 3),
        resolved_total_cores=resolved_total_cores, resolved_total_mem_gb=round(resolved_total_mem_gb, 3),
        resolved_cpu_per_task=resolved_cpu_per_task, resolved_mem_per_task_gb=round(resolved_mem_per_task_gb, 3),
        resolved_max_tasks=resolved_max_tasks, max_tasks_by_cpu=max_tasks_by_cpu,
        max_tasks_by_mem=max_tasks_by_mem, busy_parallel_cap=busy_parallel_cap, note=note,
    )


def resolve_busco_cmd_base(args: argparse.Namespace) -> Tuple[List[str], str]:
    busco_env = args.busco_conda_env or args.busco_env
    if args.busco_path:
        if os.path.exists(args.busco_path) and os.access(args.busco_path, os.X_OK):
            return [args.busco_path], "path"
        raise SystemExit("错误：--busco_path '{}' 不存在或不可执行。".format(args.busco_path))
    if busco_env:
        return ["conda", "run", "-n", busco_env, "busco"], "conda_env"
    found = shutil.which("busco")
    if found:
        return [found], "path_env"
    return ["conda", "run", "-n", "busco", "busco"], "conda_env_fallback"


# ============================================================
#  预检查报告
# ============================================================

def build_preflight_report(args, files_to_process, extra_checks, snapshot, plan,
                           completed_count=0, pending_count=0, use_offline=False):
    rep = PreflightReport(ok=True)
    rep.sample_count = len(files_to_process)
    rep.completed_count = completed_count
    rep.pending_count = pending_count
    rep.extra_checks = dict(extra_checks)

    if len(files_to_process) == 0:
        rep.errors.append("未找到任何待处理样品。")

    if not args.only_summary:
        has_auto = any(x in args.extra_parameters_busco for x in ["--auto-lineage", "-a"])
        if not args.lineage and not has_auto:
            rep.errors.append("未提供 --lineage (-l) 且未启用自动谱系推断。")
        if use_offline:
            rep.infos.append("BUSCO 离线模式: [已启用] (--offline)。检测到本地 lineage 目录，将禁止联网。")
        else:
            rep.warnings.append("BUSCO 离线模式: [已禁用]。BUSCO 将尝试联网，请确保网络畅通！")

    out_names = [out for _, out in files_to_process]
    c = Counter(out_names)
    rep.duplicate_output_names = sorted([k for k, v in c.items() if v > 1])
    if rep.duplicate_output_names:
        rep.errors.append("发现重复输出样品名: {}".format(', '.join(rep.duplicate_output_names[:20])))

    seen_fp = set()
    for fp, _ in files_to_process:
        if fp in seen_fp:
            rep.warnings.append("重复输入文件: {}".format(fp))
        seen_fp.add(fp)
        if not os.path.exists(fp):
            rep.errors.append("输入文件不存在: {}".format(fp))
            continue
        if not os.access(fp, os.R_OK):
            rep.unreadable_inputs.append(fp)
        try:
            if os.path.getsize(fp) == 0:
                rep.zero_size_inputs.append(fp)
        except Exception:
            rep.unreadable_inputs.append(fp)

    if rep.zero_size_inputs:
        rep.errors.append("存在空输入文件: {}".format(', '.join(rep.zero_size_inputs[:20])))
    if rep.unreadable_inputs:
        rep.errors.append("存在不可读输入文件: {}".format(', '.join(rep.unreadable_inputs[:20])))

    for key, label in [("missing_selected", "contigs 列表未匹配"),
                        ("missing_expected_fasta_dirs", "SPAdes 子目录缺少 fasta"),
                        ("missing_list_entries", "文件列表未匹配")]:
        items = extra_checks.get(key) or []
        if items:
            rep.warnings.append("{}: {}".format(label, ', '.join(list(items)[:20])))

    try:
        busco_cmd_base, source = resolve_busco_cmd_base(args)
        rep.infos.append("BUSCO 调用来源: {} -> {}".format(source, ' '.join(busco_cmd_base)))
    except SystemExit as e:
        rep.errors.append(str(e))

    rep.infos.append("系统快照: total_cores={}, total_mem_gb={:.2f}, avail_mem_gb={:.2f}, cpu_cap={}, mem_cap_gb={:.2f}, busy_score={:.2f}, busy_level={}".format(
        snapshot.total_cores, snapshot.total_mem_gb, snapshot.avail_mem_gb, snapshot.cpu_hard_cap_cores,
        snapshot.mem_hard_cap_gb, snapshot.busy_score, snapshot.busy_level
    ))
    rep.infos.append("资源计划: total_cores={}, total_mem_gb={:.2f}, cpu_per_task={}, mem_per_task_gb={:.2f}, max_tasks={}".format(
        plan.resolved_total_cores, plan.resolved_total_mem_gb, plan.resolved_cpu_per_task,
        plan.resolved_mem_per_task_gb, plan.resolved_max_tasks
    ))
    rep.infos.append("输出目录命名: --out_name_style={}".format(extra_checks.get("out_name_style", "no_ext")))

    if rep.errors:
        rep.ok = False
    return rep


def print_preflight_report(rep, plan, preflight_path, plan_path, use_offline=False):
    print("\n" + "#" * 78)
    print("预检查报告")
    print("#" * 78)
    print("  总样品数       : {}".format(rep.sample_count))
    print("  已完成(跳过)   : {}".format(rep.completed_count))
    print("  待运行         : {}".format(rep.pending_count))
    print("  离线模式       : {}".format("开启 (安全，不联网)" if use_offline else "关闭 (将尝试联网)"))
    print("  输出预检查 JSON: {}".format(preflight_path))
    print("  输出资源规划   : {}".format(plan_path))

    print("\n资源规划摘要（基于 {} 个待运行任务）：".format(rep.pending_count))
    for k, v in [
        ("host_total_cores", plan.host_total_cores),
        ("host_total_mem_gb", "{:.2f}".format(plan.host_total_mem_gb)),
        ("host_available_mem_gb", "{:.2f}".format(plan.host_available_mem_gb)),
        ("cpu_hard_cap_cores", plan.cpu_hard_cap_cores),
        ("mem_hard_cap_gb", "{:.2f}".format(plan.mem_hard_cap_gb)),
        ("busy_score", "{:.2f} ({})".format(plan.busy_score, plan.busy_level)),
        ("resolved_total_cores", plan.resolved_total_cores),
        ("resolved_total_mem_gb", "{:.2f}".format(plan.resolved_total_mem_gb)),
        ("resolved_cpu_per_task", plan.resolved_cpu_per_task),
        ("resolved_mem_per_task", "{:.2f} GB".format(plan.resolved_mem_per_task_gb)),
        ("resolved_max_tasks", plan.resolved_max_tasks),
    ]:
        print("  {:<24s}: {}".format(k, v))

    if rep.errors:
        print("\n[ERROR] ({} 个)".format(len(rep.errors)))
        for x in rep.errors:
            print("  - " + x)
    if rep.warnings:
        print("\n[WARNING] ({} 个)".format(len(rep.warnings)))
        for x in rep.warnings:
            print("  - " + x)
    if rep.infos:
        print("\n[INFO]")
        for x in rep.infos:
            print("  - " + x)
    print("\n预检查状态: {}".format("PASS" if rep.ok else "FAIL"))


# ============================================================
#  Snakemake
# ============================================================

def build_snakefile(files_to_process, args, busco_cmd_base, output_dir_abs,
                    resolved_cpu_per_task, resolved_mem_per_task_gb, resume_strategy):
    tasks_json = json.dumps([{"fp": fp, "out": out} for fp, out in files_to_process], ensure_ascii=False)
    output_dir_j = json.dumps(os.path.abspath(output_dir_abs), ensure_ascii=False)
    busco_cmd_j = json.dumps(list(busco_cmd_base), ensure_ascii=False)
    lineage_j = json.dumps(os.path.abspath(args.lineage) if args.lineage else "", ensure_ascii=False)
    mode_j = json.dumps(str(args.mode), ensure_ascii=False)
    extra_args_j = json.dumps(list(args.extra_parameters_busco), ensure_ascii=False)
    resume_strategy_j = json.dumps(str(resume_strategy), ensure_ascii=False)
    cpu_val = int(resolved_cpu_per_task)
    mem_val = int(math.ceil(float(resolved_mem_per_task_gb)))

    lines = [
        "import os, subprocess, glob, datetime",
        "",
        f"TASKS      = {tasks_json}",
        'TASK_MAP   = {t["out"]: t["fp"] for t in TASKS}',
        'ALL_OUTS   = [t["out"] for t in TASKS]',
        f"OUTPUT_DIR = {output_dir_j}",
        f"BUSCO_CMD  = {busco_cmd_j}",
        f"LINEAGE    = {lineage_j}",
        f"MODE       = {mode_j}",
        f"CPU        = {cpu_val}",
        f"EXTRA_ARGS = {extra_args_j}",
        f"RESUME_STR = {resume_strategy_j}",
        "",
        "rule all:",
        "    input:",
        '        expand(OUTPUT_DIR + "/{name}/short_summary.done", name=ALL_OUTS)',
        "",
        "rule busco:",
        "    input:",
        "        fasta = lambda wildcards: TASK_MAP[wildcards.name]",
        "    output:",
        '        done  = OUTPUT_DIR + "/{name}/short_summary.done"',
        "    threads: CPU",
        "    resources:",
        f"        mem_gb  = {mem_val},",
        "        disk_io = 1",
        "    log:",
        '        OUTPUT_DIR + "/logs/{name}.log"',
        "    run:",
        "        out_subdir = os.path.join(OUTPUT_DIR, wildcards.name)",
        '        summary_pattern = os.path.join(out_subdir, "short_summary*.txt")',
        "",
        "        if glob.glob(summary_pattern):",
        '            print("[SKIP] already done: " + wildcards.name)',
        '            with open(output.done, "w", encoding="utf-8", errors="ignore") as fh:',
        '                fh.write("skipped at " + str(datetime.datetime.now()) + "\\n")',
        "        else:",
        "            cmd = BUSCO_CMD + [",
        '                "-i", str(input.fasta),',
        '                "-o", wildcards.name,',
        '                "-m", MODE,',
        '                "-c", str(threads),',
        "            ]",
        "            if LINEAGE:",
        "                cmd.extend(['-l', LINEAGE])",
        "            if os.path.exists(out_subdir):",
        "                if not any(x in EXTRA_ARGS for x in ['-f', '--force', '-r', '--restart']):",
        "                    if RESUME_STR == 'restart':",
        "                        cmd.append('--restart')",
        "                    else:",
        "                        cmd.append('--force')",
        "            if EXTRA_ARGS:",
        "                cmd.extend(EXTRA_ARGS)",
        "",
        '            os.makedirs(os.path.join(OUTPUT_DIR, "logs"), exist_ok=True)',
        '            print("[RUN] " + " ".join(cmd))',
        '            with open(str(log[0]), "w", encoding="utf-8", errors="ignore") as logf:',
        "                ret = subprocess.run(cmd, stdout=logf, stderr=subprocess.STDOUT, cwd=OUTPUT_DIR).returncode",
        "",
        "            if ret != 0:",
        '                tail = ""',
        "                try:",
        '                    with open(str(log[0]), "r", encoding="utf-8", errors="ignore") as fh:',
        "                        tail = ''.join(fh.readlines()[-80:])",
        "                except Exception as e:",
        '                    tail = "[log read error: " + str(e) + "]"',
        '                raise RuntimeError("BUSCO failed (exit " + str(ret) + "): " + wildcards.name + "  log: " + str(log[0]) + "\\n" + tail)',
        '            with open(output.done, "w", encoding="utf-8", errors="ignore") as fh:',
        '                fh.write("completed at " + str(datetime.datetime.now()) + "\\n")',
    ]
    return "\n".join(lines) + "\n"


def is_lock_exception_text(text):
    t = (text or "").lower()
    return "lockexception" in t or "directory cannot be locked" in t or "another snakemake process" in t


def run_snakemake_command(cmd, log_path=None):
    p = subprocess.run(cmd, stdout=subprocess.PIPE, stderr=subprocess.STDOUT, text=True)
    out = p.stdout or ""
    if log_path:
        try:
            with open(log_path, "w", encoding="utf-8", errors="ignore") as fh:
                fh.write(out)
        except Exception:
            pass
    if out:
        print(out, end="" if out.endswith("\n") else "\n")
    return int(p.returncode), out


def resolve_snakemake_cmd_base(args):
    smk_path = getattr(args, "smk_path", "") or ""
    smk_env = getattr(args, "smk_env", "") or ""
    if smk_path:
        if os.path.exists(smk_path) and os.access(smk_path, os.X_OK):
            return [os.path.abspath(smk_path)], "path"
        raise SystemExit("错误：--smk_path '{}' 不存在或不可执行。".format(smk_path))
    if smk_env:
        if not shutil.which("conda"):
            raise SystemExit("错误：指定了 --smk_env，但未找到 conda。")
        return ["conda", "run", "-n", smk_env, "snakemake"], "conda_env"
    found = shutil.which("snakemake")
    if found:
        return [found], "path_env"
    return [], "missing"


def run_with_snakemake(files_to_process, args, busco_cmd_base, snakemake_cmd_base,
                       output_dir, plan, extra_snakemake_args, debug=False, dag=False):
    output_dir_abs = os.path.abspath(output_dir)
    logs_dir = os.path.join(output_dir_abs, "logs")
    smk_path = os.path.join(output_dir_abs, ".snakemake_busco.smk")
    cli_log = os.path.join(logs_dir, "snakemake_cli.log")
    os.makedirs(output_dir_abs, exist_ok=True)
    os.makedirs(logs_dir, exist_ok=True)

    snake_text = build_snakefile(files_to_process, args, busco_cmd_base, output_dir_abs,
                                  plan.resolved_cpu_per_task, plan.resolved_mem_per_task_gb, args.resume_strategy)
    with open(smk_path, "w", encoding="utf-8") as f:
        f.write(snake_text)
    print("\nSnakefile 已生成: {}".format(smk_path))

    parse_cmd = list(snakemake_cmd_base) + ["--snakefile", smk_path, "-n", "--quiet"]
    print("\nSnakefile 解析自检：\n  " + " ".join(shlex.quote(x) for x in parse_cmd))
    prc, pout = run_snakemake_command(parse_cmd, log_path=cli_log)
    if prc != 0:
        print("\n[ERROR] Snakefile 解析失败")
        print(pout or "")
        return int(prc)

    if dag:
        dag_cmd = list(snakemake_cmd_base) + ["--snakefile", smk_path, "--dag", "--quiet"]
        dag_out = os.path.join(output_dir_abs, "busco_dag.png")
        result = subprocess.run(dag_cmd, capture_output=True, text=True)
        if result.returncode == 0:
            try:
                subprocess.run(["dot", "-Tpng", "-o", dag_out], input=result.stdout, text=True)
                print("DAG 图已保存: {}".format(dag_out))
            except FileNotFoundError:
                print("未找到 dot (graphviz)")
            return 0
        return int(result.returncode)

    disk_io_limit = max(1, min(plan.resolved_max_tasks, 4))
    snakemake_cmd = list(snakemake_cmd_base) + [
        "--snakefile", smk_path,
        "--cores", str(plan.resolved_total_cores),
        "--jobs", str(plan.resolved_max_tasks),
        "--resources",
        "mem_gb={}".format(int(max(1, math.floor(plan.resolved_total_mem_gb)))),
        "disk_io={}".format(disk_io_limit),
        "--rerun-incomplete", "--keep-going", "--printshellcmds",
        "--latency-wait", str(int(args.smk_latency_wait)),
        "--restart-times", str(int(args.smk_restart_times)),
    ]
    if extra_snakemake_args:
        snakemake_cmd.extend(list(extra_snakemake_args))

    print("\n启动 Snakemake：\n  " + " ".join(shlex.quote(x) for x in snakemake_cmd))
    rc, out = run_snakemake_command(snakemake_cmd, log_path=cli_log)

    if rc != 0 and args.smk_unlock_on_lock_exception == "yes" and is_lock_exception_text(out):
        print("\n检测到锁异常，尝试 --unlock...")
        unlock_cmd = list(snakemake_cmd_base) + ["--snakefile", smk_path, "--unlock"]
        urc, _ = run_snakemake_command(unlock_cmd, log_path=cli_log)
        if urc == 0:
            rc, out = run_snakemake_command(snakemake_cmd, log_path=cli_log)

    if rc == 0:
        print("\nSnakemake 工作流成功完成！")
    else:
        print("\nSnakemake 退出码: {}\n  重新运行相同命令可自动续跑。".format(rc))
    return rc


# ============================================================
#  BUSCO 结果汇总
# ============================================================

def summarize_busco_results(results_input_dir, summary_output_file, files_list_for_summary, debug=False):
    busco_dirs = []
    if files_list_for_summary:
        try:
            with open(files_list_for_summary) as f:
                for line in f:
                    name = line.strip()
                    if name:
                        candidate = os.path.join(results_input_dir, name)
                        if os.path.isdir(candidate):
                            busco_dirs.append(candidate)
        except Exception as e:
            print("读取汇总列表失败: {}".format(e))
    else:
        if os.path.isdir(results_input_dir):
            for entry in os.listdir(results_input_dir):
                subdir = os.path.join(results_input_dir, entry)
                if os.path.isdir(subdir) and glob.glob(os.path.join(subdir, "short_summary*.txt")):
                    busco_dirs.append(subdir)

    print("\n开始汇总 {} 个 BUSCO 结果目录...".format(len(busco_dirs)))
    os.makedirs(os.path.dirname(os.path.abspath(summary_output_file)), exist_ok=True)

    pat_version = re.compile(r"BUSCO version is:\s*(\S+)")
    pat_lineage = re.compile(r"lineage dataset is:\s*([^( \n]+)")
    pat_num = re.compile(r"number of BUSCOs:\s*(\d+)")
    pat_single = re.compile(r"S:([0-9.]+)%")
    pat_dup = re.compile(r"D:([0-9.]+)%")
    pat_frag = re.compile(r"F:([0-9.]+)%")
    pat_miss = re.compile(r"M:([0-9.]+)%")

    def _m(pattern, text, default="NA"):
        hit = pattern.search(text)
        return hit.group(1) if hit else default

    with open(summary_output_file, "w", encoding="utf-8") as sf:
        sf.write("busco\tsingle\tduplicate\tfragmented\tmissing\tbusco_version\tused_database\tnumber_of_BUSCOs\n")
        for busco_dir in sorted(busco_dirs):
            sf_files = glob.glob(os.path.join(busco_dir, "short_summary*.txt"))
            if not sf_files:
                continue
            try:
                with open(sorted(sf_files)[0], encoding="utf-8", errors="ignore") as f:
                    content = f.read()
            except Exception:
                continue
            name = os.path.basename(busco_dir)
            sf.write("\t".join([name, _m(pat_single, content, "0"), _m(pat_dup, content, "0"),
                                _m(pat_frag, content, "0"), _m(pat_miss, content, "0"),
                                _m(pat_version, content), _m(pat_lineage, content), _m(pat_num, content)]) + "\n")
    print("BUSCO summary 已生成: {}".format(summary_output_file))


# ============================================================
#  参数检查
# ============================================================

def check_input_parameters(args):
    if args.only_summary:
        return
    if not args.input and not args.files_lists_for_busco2summary:
        sys.exit("错误：必须提供 -i 参数。")
    if args.input and not (os.path.isfile(args.input) or os.path.isdir(args.input)):
        sys.exit("错误：输入路径 {} 不存在。".format(args.input))
    if args.mode not in VALID_MODES:
        sys.exit("错误：--mode '{}' 不合法。".format(args.mode))


# ============================================================
#  主函数
# ============================================================

def main() -> None:
    sys_cores = get_system_cores()
    sys_mem_gb, sys_avail_mem_gb = get_system_mem_info_gb()
    prog = os.path.basename(sys.argv[0])

    epilog_text = textwrap.dedent("""
    系统检测：{cores} 核心 / {mem:.1f} GB 总内存 / {avail:.1f} GB 可用
    硬上限：CPU/MEM 不超过 90%

    推荐流程：
    1) 预检查：python3 {prog} -i INPUT -o OUT -l LINEAGE
    2) 执行：  python3 {prog} -i INPUT -o OUT -l LINEAGE -y
    """).format(cores=sys_cores, mem=sys_mem_gb, avail=sys_avail_mem_gb, prog=prog)

    parser = argparse.ArgumentParser(
        description="整合 BUSCO 批量运行与结果汇总（Snakemake 全局资源托管版 v{}）".format(__version__),
        epilog=epilog_text, formatter_class=argparse.RawTextHelpFormatter,
    )

    parser.add_argument("--version", action="version", version="%(prog)s v{}".format(__version__))
    parser.add_argument("-i", "--input", help="输入目录或单个文件")
    parser.add_argument("-o", "--output", required=True, help="BUSCO 输出目录")
    parser.add_argument("-l", "--lineage", default="", help="BUSCO lineage 路径")
    parser.add_argument("-m", "--mode", default="genome", choices=["genome", "transcriptome", "proteins"])
    parser.add_argument("-c", "--cpu", default=0, type=int, help="每任务 CPU 数，0=自动")
    parser.add_argument("--max_tasks", default=0, type=int, help="最大并发数，0=自动")
    parser.add_argument("--mem_per_task_gb", type=float, default=0.0)
    parser.add_argument("--total_cores", type=int, default=0)
    parser.add_argument("--total_mem_gb", type=float, default=0.0)
    parser.add_argument("--auto_max_jobs", type=int, default=8)
    parser.add_argument("--dynamic_busy_guard", choices=["yes", "no"], default="yes")
    parser.add_argument("--busy_sample_sec", type=float, default=0.6)

    parser.set_defaults(use_snakemake=True)
    parser.add_argument("--smk", "--snakemake", dest="use_snakemake", action="store_true")
    parser.add_argument("--no_smk", "--no_snakemake", dest="use_snakemake", action="store_false")
    parser.add_argument("--smk_path", default="")
    parser.add_argument("--smk_env", default="")
    parser.add_argument("--snakemake_args", nargs="*", default=[])
    parser.add_argument("--smk_latency_wait", type=int, default=60)
    parser.add_argument("--smk_restart_times", type=int, default=1)
    parser.add_argument("--smk_unlock_on_lock_exception", choices=["yes", "no"], default="yes")
    parser.add_argument("--dag", action="store_true")

    parser.add_argument("-y", "--yes", action="store_true", help="确认执行（不加只做预检查）")
    parser.add_argument("--preflight_report", default="preflight_report.json")
    parser.add_argument("--resource_plan_output", default="resource_plan.json")

    parser.add_argument("--store_type_contigs", choices=["subdir", "flat"], default="flat")
    parser.add_argument("--suffix_fasta", default="")
    parser.add_argument("--contigs_lists", default="")
    parser.add_argument("--files_lists_for_busco2summary", default="")

    parser.add_argument("--busco_for_spades_dir", "-bfsd", action="store_true")
    parser.add_argument("--suffix_subdir", default="")
    parser.add_argument("--which_spades_fasta", default="contigs.fasta",
                        choices=["contigs.fasta", "scaffolds.fasta", "transcripts.fasta"])

    # v2.8.6: 新增 --out_name_style 参数
    parser.add_argument("--out_name_style", choices=sorted(OUT_NAME_STYLE_CHOICES), default="no_ext",
                        help="输出目录命名策略：\n"
                             "  no_ext    (默认) 去掉扩展名 → GCA_xxx_genomic\n"
                             "  keep_ext  保留完整文件名   → GCA_xxx_genomic.fna\n"
                             "  stem_only 只到第一个点     → GCA_000326865\n"
                             "旧版本行为等价于 keep_ext（这是个 BUG，会导致无法跳过已完成样品）。")

    parser.add_argument("--if_summary_busco_results", "-fsbr", default="yes", choices=["yes", "no"])
    parser.add_argument("--files_lists_for_summary_busco", "-flfs", default="")
    parser.add_argument("--outputdir_for_buscos2summary", "-ofb2s", default="")
    parser.add_argument("--summary_output", default="summary_of_busco_res.txt")
    parser.add_argument("--only_summary", action="store_true")

    parser.add_argument("--busco_path", default="")
    parser.add_argument("--busco_env", default="busco")
    parser.add_argument("--busco_conda_env", default="")
    parser.add_argument("--extra_parameters_busco", "-epb", nargs="*", default=[])

    parser.add_argument("--resume_strategy", choices=["force", "restart"], default="force")
    parser.add_argument("--offline_mode", choices=["auto", "yes", "no"], default="auto")
    parser.add_argument("--auto_detect_mode", default="yes", choices=["yes", "no"])
    parser.add_argument("--extension_files", default=".fasta,.fas,.fna,.fa,.pep,.faa")
    parser.add_argument("--debug", action="store_true")

    args = parser.parse_args()

    global ALLOWED_EXTENSIONS
    ALLOWED_EXTENSIONS = [e.strip() for e in args.extension_files.split(",") if e.strip()]
    print("允许处理的扩展名: {}".format(', '.join(ALLOWED_EXTENSIONS)))
    print("运行后端: {}".format("SNAKEMAKE" if args.use_snakemake else "THREADS"))
    print("输出目录命名: --out_name_style={}".format(args.out_name_style))
    print("固定资源硬上限: CPU <= 90%, MEM <= 90%")

    # 离线模式
    use_offline = False
    if not args.only_summary:
        if args.offline_mode == "auto":
            use_offline = bool(args.lineage and os.path.isdir(args.lineage))
        elif args.offline_mode == "yes":
            use_offline = True
        if use_offline and "--offline" not in args.extra_parameters_busco:
            args.extra_parameters_busco.append("--offline")

    auto_detect = args.auto_detect_mode.lower() == "yes"

    if args.only_summary:
        summary_out_dir = args.outputdir_for_buscos2summary or args.output
        os.makedirs(summary_out_dir, exist_ok=True)
        if not args.yes:
            print("\n当前为 only_summary 预检查模式。确认后加 -y 执行。")
            return
        summarize_busco_results(args.output, os.path.join(summary_out_dir, args.summary_output),
                                 args.files_lists_for_summary_busco, debug=args.debug)
        return

    check_input_parameters(args)
    if args.input:
        check_fasta_structure(args.input, debug=args.debug)

    os.makedirs(args.output, exist_ok=True)

    temp_dir = None
    suggested_modes = None
    gather_extra = {}

    if args.files_lists_for_busco2summary:
        files_to_process, suggested_modes, gather_extra = gather_files_from_list(
            args.files_lists_for_busco2summary, args.input,
            out_name_style=args.out_name_style, suggest_mode=auto_detect, debug=args.debug)
    elif args.busco_for_spades_dir:
        files_to_process, temp_dir, suggested_modes, gather_extra = gather_spades_busco_files(
            args.input, args.suffix_subdir, args.which_spades_fasta, args.output,
            out_name_style=args.out_name_style, suggest_mode=auto_detect, debug=args.debug)
    else:
        files_to_process, suggested_modes, gather_extra = gather_contigs_files(
            args.input, args.store_type_contigs, args.suffix_fasta, args.contigs_lists,
            out_name_style=args.out_name_style, suggest_mode=auto_detect, debug=args.debug)

    if auto_detect and suggested_modes:
        mc = Counter(suggested_modes.values())
        if mc:
            recommended = translate_mode_to_busco_param(mc.most_common(1)[0][0])
            if args.mode != recommended:
                print("\n注意：您指定了 -m {}，推荐 -m {}".format(args.mode, recommended))

    completed, pending = scan_completed_results(args.output, files_to_process)

    print("\n样品统计: 总计 {} | 已完成(跳过) {} | 待运行 {}".format(
        len(files_to_process), len(completed), len(pending)))
    if completed:
        for _, out_name in completed[:15]:
            print("  [SKIP] {}".format(out_name))
        if len(completed) > 15:
            print("  ... 及其他 {} 个已完成".format(len(completed) - 15))

    snapshot = get_system_snapshot(sample_sec=args.busy_sample_sec)
    plan = compute_resource_plan(args=args, num_tasks=max(1, len(pending)), snapshot=snapshot)
    preflight = build_preflight_report(args=args, files_to_process=files_to_process, extra_checks=gather_extra,
                                       snapshot=snapshot, plan=plan,
                                       completed_count=len(completed), pending_count=len(pending),
                                       use_offline=use_offline)

    preflight_path = os.path.join(args.output, args.preflight_report)
    plan_path = os.path.join(args.output, args.resource_plan_output)
    ensure_json_dump(preflight_path, asdict(preflight))
    ensure_json_dump(plan_path, asdict(plan))
    print_preflight_report(preflight, plan, preflight_path, plan_path, use_offline=use_offline)

    if not preflight.ok:
        print("\n预检查未通过，请修正 [ERROR] 后重新运行。")
        sys.exit(1)

    if len(pending) == 0:
        print("\n所有 {} 个样品均已完成。".format(len(files_to_process)))
        if args.if_summary_busco_results.lower() == "yes":
            summarize_busco_results(args.output, os.path.join(args.output, args.summary_output),
                                     args.files_lists_for_summary_busco, debug=args.debug)
        return

    if not args.yes:
        print("\n当前仅执行了预检查和资源规划，还没有真正启动 BUSCO。")
        print("如果您确认信息（特别是离线模式的状态）无误，请加 -y 重新运行，例如：")
        replay_args = ' '.join([x for x in sys.argv[1:] if x not in ["-y", "--yes"]])
        print("  python3 {} {} -y".format(prog, replay_args))
        return

    busco_cmd_base, source = resolve_busco_cmd_base(args)
    print("\nBUSCO 调用来源: {} | 命令: {}".format(source, " ".join(busco_cmd_base)))

    snakemake_cmd_base, smk_source = resolve_snakemake_cmd_base(args)
    if args.use_snakemake:
        if snakemake_cmd_base:
            rc = run_with_snakemake(files_to_process=files_to_process, args=args, busco_cmd_base=busco_cmd_base,
                                    snakemake_cmd_base=snakemake_cmd_base, output_dir=args.output, plan=plan,
                                    extra_snakemake_args=args.snakemake_args, debug=args.debug, dag=args.dag)
            if rc != 0:
                sys.exit(rc)
            if args.dag:
                return
        else:
            print("未找到 snakemake，回退到 threads 后端。")
            args.use_snakemake = False

    if not args.use_snakemake:
        print("\n使用 Threads 后端 (cpu_per_task={}, max_tasks={})".format(
            plan.resolved_cpu_per_task, plan.resolved_max_tasks))
        executor = ThreadPoolExecutor(max_workers=plan.resolved_max_tasks)
        futures = []
        for input_file, out_name in pending:
            futures.append(executor.submit(run_busco, input_file, out_name, args.lineage, args.mode,
                                           plan.resolved_cpu_per_task, args.extra_parameters_busco,
                                           os.path.abspath(args.output), busco_cmd_base, args.resume_strategy, args.debug))
        total = len(futures)
        for i, fut in enumerate(as_completed(futures), 1):
            rc = fut.result()
            print("[{}] 任务 {}/{} 完成，返回码: {}".format("OK" if rc == 0 else "FAIL", i, total, rc))
        executor.shutdown(wait=True)

    if not args.dag and args.if_summary_busco_results.lower() == "yes":
        summarize_busco_results(args.output, os.path.join(args.output, args.summary_output),
                                 args.files_lists_for_summary_busco, debug=args.debug)

    if temp_dir and os.path.isdir(temp_dir):
        try:
            shutil.rmtree(temp_dir)
        except Exception:
            pass

    print("\n所有任务处理完成。")


if __name__ == "__main__":
    try:
        main()
    except SystemExit:
        raise
    except KeyboardInterrupt:
        print("\n用户中断 (Ctrl+C)")
        sys.exit(130)
    except Exception as e:
        print("\n" + "=" * 78)
        print("未预料异常: {} | {}".format(type(e).__name__, e))
        traceback.print_exc()
        print("Python: {} | 命令: {}".format(sys.version, ' '.join(sys.argv)))
        sys.exit(99)

