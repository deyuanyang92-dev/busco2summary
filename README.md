# busco2summary

**多样品 BUSCO 结果的精细整理：汇总表、S/D/F/M 矩阵、逐基因统计（重视 F）、按样品占有率输出基因 ID 列表、序列提取和交互式 HTML 报告。**
*Detailed multi-sample BUSCO post-processing — summaries, S/D/F/M matrices, per-gene statistics
that treat fragmented BUSCOs as signal, occupancy-based gene-ID lists, sequence extraction and an
interactive HTML report.* (English section below.)

YZZ Team · De-Yuan Yang · MIT License · Python ≥ 3.8，无必需依赖（matplotlib / openpyxl 可选）

---

## 为什么要这个工具

大多数 BUSCO 系统发育流程（BUSCO_phylogenomics、busco2phylo-nf 等）只用 **完整单拷贝 (S)**，
片段化 (F) 直接丢弃。ZF（张峰课题组）`BUSCO_extraction.sh` 主版本的 BUSCO v5 分支是少数例外：
它把 F 序列并入单拷贝，即用 **S + F**（本工具的 `SF` 集合）；其 `version2/` 与 BUSCO v3 分支只用 S。在真实数据中这会损失大量位点：

| 数据（真实） | 占有率阈值 | 仅 S | S + F | 任意命中 S+D+F |
|---|---|---|---|---|
| 192 个环节动物基因组，BUSCO 6 / metazoa_odb12 | ≥ 80 % | 136 | **521** | 634 |
|  | ≥ 90 % | 1 | **90** | 414 |
| 77 个转录组+基因组，BUSCO 5.8 / metazoa_odb12 | ≥ 70 % | 46 | **519** | 585 |
|  | ≥ 90 % | 0 | **171** | 371 |

在同一数据中，F 片段的中位长度是同一 BUSCO 完整拷贝的 **52 %**，带有大量序列信息。
busco2summary 把 F（以及转录组里大量的 D）作为一等信息来统计，同时给出质量指标，
让你根据证据决定是否使用它们。

## 功能

* **读取**：BUSCO v3、v4–v6（genome / transcriptome / proteins，auto-lineage、batch 模式），compleasm；
  自动递归查找、多输入目录、样品表重命名（兼容 busco_2_matrix v0.3 的 `#original_name` 表头）。
* **样品汇总** `busco_summary.tsv/csv/xlsx`：C/S/D/F/M 数量和百分比、E（内部终止密码子，BUSCO ≥ 5.8）、
  组装统计（scaffold/contig N50 等）、版本、模式、基因预测器、序列文件数、状态和备注。
* **质控**：full_table 与 short_summary 一致性检查（MISMATCH），空样品（EMPTY，例如误下载的线粒体基因组），
  LOW_C / HIGH_D / HIGH_F / HIGH_M / HIGH_E / MULTI_RUN 标记，混合版本/模式/谱系警告。
* **矩阵**（`matrices/`）：状态矩阵 S/D/F/M、拷贝数矩阵、长度矩阵、片段比例矩阵、
  二元矩阵 S / SD / SF / SDF、长表（tidy）。
* **逐基因统计** `gene_stats.tsv`：n_S n_D n_F n_M 及比例、四种占有率、F gain、
  完整拷贝和片段的中位长度、片段比例中位数、D 平均拷贝数，标记 F_PRONE / D_PRONE / F_RESCUED_50。
* **逐样品统计** `sample_stats.tsv`：基于 full_table 的计数、S+F 完整度排名、片段平均比例。
* **基因 ID 列表** `gene_lists/`：任意占有率阈值（默认 50 60 70 75 80 85 90 95 100 %）× 四种状态集合，
  外加 `≥ K 个类群`（`--min-taxa`，默认 3）。阈值规则为整数精确规则：`出现样品数 ≥ ceil(t × N)`。
  另有 样品完整度 × 基因占有率网格（`grid_*.tsv`）和贪心删样建议（`taxon_drop_*.tsv`：删哪个样品能换回多少基因）。
* **`--min-frag-ratio R`**：只把长度 ≥ R × 参考长度的 F 计为存在（参考长度默认是同一 BUSCO 在本数据集中完整拷贝的中位长度）。
* **序列提取**：按基因 / 按样品输出 FASTA；D 或多条记录时选 full_table 得分最高（`best`）/ 最长 / 第一条；
  表头模板；终止密码子处理；**miniprot 结果（BUSCO 6 默认，无 .fna）可由基因组 + miniprot 比对精确重建 CDS**，
  每条序列都与 BUSCO 蛋白比对验证（`*_verified`）。
* **报告**：离线可用的单文件交互式 HTML（概览、样品、基因、占有率探索器可直接下载基因 ID、矩阵热图，支持深色模式），
  静态图 PNG/PDF/SVG（BUSCO 条形图、占有率曲线、片段完整度分布、状态热图）。
* **稳健与可复现**：单个样品文件损坏不会中断整批；解析缓存可断点续跑；记录日志和全部参数；多线程读取；
  可在 Snakemake 中调用（`workflow/Snakefile`）。

## 安装

```bash
git clone https://github.com/deyuanyang92-dev/busco2summary.git
cd busco2summary
pip install .            # 或: pip install ".[all]"  (matplotlib + openpyxl)
# conda: conda env create -f environment.yml
busco2summary --version   # 简写命令: b2s
```

## 快速开始

```bash
# 一条命令完成全部分析
busco2summary run -i busco_results/ -o b2s_out -j 16

# 只要 50/60/90 % 占有率的基因 ID 列表（S 与 S+F）
busco2summary genelists -i busco_results/ -o b2s_lists -t 50 60 90 --sets S SF

# 指定谱系；多个输入目录；样品表重命名
busco2summary run -i runA/ runB/ -s samples.tsv --lineage metazoa_odb12 -o out

# 片段至少覆盖参考长度一半才计入
busco2summary run -i busco_results/ -o out --min-frag-ratio 0.5

# 提取 S+F、占有率 ≥ 90 % 的基因的蛋白和核酸序列（重复拷贝取得分最高者）
busco2summary extract -i busco_results/ -o seqs --set SF --occupancy 90 --seqtype both --dup best

# miniprot 结果没有 .fna：用基因组重建 CDS
busco2summary extract -i busco_results/ -o seqs --genes b2s_out/gene_lists/SF/SF_occ090.txt \
    --seqtype nucleotide --genome-map genomes.tsv      # 两列: 样品名<TAB>基因组.fasta
#   或 --genome-from-summary（使用 BUSCO 记录的输入文件路径）

# 只看目录里有什么
busco2summary scan -i busco_results/ -o scan_out
```

样品名默认去掉 `.fna/.fasta/...` 后缀（`--name-style no_ext`）；`--name-style accession` 只保留 GCA/GCF/SRR 号；
`--replace-dots` 把 `.` 换成 `_`（v0.3 的默认行为）。

## 输出

```
b2s_out/
├── busco_summary.tsv / .csv (.xlsx)   每个样品一行（所有样品，含失败/空样品）
├── dataset_report.txt                  版本、模式、谱系、状态统计与警告
├── sample_stats.tsv                    进入矩阵的样品：基于 full_table 的统计和 QC 标记
├── gene_stats.tsv                      每个 BUSCO：S/D/F/M 数量、占有率、F gain、长度、标记
├── gene_lists/
│   ├── occupancy_summary.tsv           每个 集合 × 阈值 的基因数、矩阵填充率、F gain
│   ├── occupancy_wide.tsv              阈值 × 集合 宽表
│   ├── S/ SD/ SF/ SDF/                 <集合>_occ050.txt … <集合>_min3.txt（基因 ID 列表）
│   ├── grid_<集合>.tsv                 样品完整度 × 基因占有率
│   └── taxon_drop_S.tsv / _SF.tsv      贪心删样建议
├── matrices/                           status / copy_number / length / fragment_ratio / binary_*
├── plots/                              busco_bars, occupancy_curves, fragment_completeness, status_heatmap
├── busco2summary_report.html           交互式报告
├── sequences/<列表>/ (使用 --extract)   by_gene_faa/ by_sample_faa/ … + extraction_manifest.tsv
├── name_mapping.tsv, run_parameters.json, busco2summary.log
```

### 状态集合

| 集合 | 视为“存在”的状态 | 用途 |
|---|---|---|
| `S`   | 完整单拷贝 | 经典系统发育基因组学 |
| `SD`  | 完整（单拷贝 + 重复） | 完整性 |
| `SF`  | 单拷贝 + 片段 | 用片段补回位点（推荐与 S 对比） |
| `SDF` | 任意命中 | 最宽松；转录组 D 很多时使用并配合 `--dup best` |

## 整理逻辑（原理速查）

### 三个"占比"，别混

把所有结果看成一张 **基因 × 样品** 的表，每格是 S / D / F / M：

| 名称 | 怎么算 | 看的是 | 用来做什么 |
|---|---|---|---|
| **样品完整度**（物种占比、species completeness） | 一个样品里"存在"的 BUSCO 数 ÷ BUSCO 总数 n | 表的**一列** | 剔除质量太差的样品（`--min-completeness`） |
| **基因占比**（gene occupancy） | 有这个基因的样品数 ÷ 样品数 N | 表的**一行** | 挑选基因：`SF_occ090.txt` = ≥ 90 % 样品都有的基因 |
| **矩阵填充率**（matrix fill） | 选出的基因 × 样品中有数据的格子比例 | 最终矩阵整体 | 报告最终矩阵的缺失程度（`occupancy_summary.tsv`） |

**为什么要基因占比**：BUSCO 基因不是每个样品都能找到。一个基因若只在少数样品出现，放进矩阵就是大片缺失；
要求高占比 → 基因少但完整，要求低占比 → 基因多但缺失多。占比阈值就是在 **基因数量 ↔ 缺失数据** 之间取舍。

**为什么先看样品完整度**：基因占比的分母是样品数 N。一个只有 20 % BUSCO 的差样品，会让几乎每个基因都"少一个样品"，
把 90 % 占比的基因数整体拉低。所以顺序是：**先剔除差样品，再按基因占比选基因**（`taxon_drop_*.tsv` 直接告诉你删哪个样品能换回多少基因）。

### 本工具的判定规则（写论文方法时可直接引用）

1. **每个样品每个 BUSCO 只有一个状态**，按 full_table 的**不同 ID** 计数（不是行数；重复拷贝、多个 miniprot 命中都会产生多行）。
   优先级 S > D > F > M；compleasm 的 `Single` = S，`Interspaced` = F，`Retrocopy` 归入 D。
2. **"存在"由状态集合决定**：`S` 只算单拷贝；`SF` 算单拷贝 + 片段；`SD` 算完整；`SDF` 任意命中。D 只在 `SD`/`SDF` 中算存在。
3. **基因占比阈值 t**：出现样品数 ≥ ⌈t × N⌉ 才入选（整数规则，没有浮点误差；N = 过滤后进入矩阵的样品数）。
   所有样品都缺失的基因不进入任何列表。
4. **片段长度过滤**（`--min-frag-ratio R`，默认关闭）：片段长度 ÷ 参考长度 ≥ R 的 F 才算存在；
   参考长度默认是该 BUSCO 在本数据集完整拷贝的中位长度。**没有参考长度的 F（全数据集都没有完整拷贝）不算存在**。
   这条规则同时作用于基因列表、`binary_*` 矩阵、`occ_*`、`n_F_qualified`、`pct_SF_qualified` 和序列提取；
   `S/D/F/M` 计数和 `status_matrix` 保留 BUSCO 原始判定，不受影响。
5. **分母 n** 取 short_summary 一行式中的 `n:`（BUSCO 自己计算百分比用的数）；没有 full_table 时，缺失数 = n − S − D − F。
6. BUSCO ≥ 6.1 配 odb12.2 时 full_table 的 Length 含义变为目标蛋白长度（BUSCO 源码），片段长度比在这种数据上意义有限，
   此时不建议用 `--min-frag-ratio`。

### 什么占比"最优"？

**没有通用的最优值**，取决于数据和下游分析：

* 高占比（≥ 90 %）：缺失少、基因少，适合串联矩阵（concatenation）对缺失敏感时；
* 中等（70 %–80 %）：常见折中；
* 低占比（50 %）：基因多，适合基因树 + 溯祖方法（如 ASTRAL），它们能容忍基因在部分样品中缺失。

建议做法：用 `occupancy_curves.png` / `occupancy_summary.tsv` 看基因数随阈值的变化，**挑 2–3 个阈值（如 50 / 70 / 90 %）
都建树，比较拓扑是否稳定**；论文里报告所用阈值、状态集合（S 或 S+F）、样品数 N 和最终矩阵填充率。
`--sets S SF` 并列输出，可直接比较 F 补回了多少基因（`F_gain`）。

### 与其他工具的逻辑对比（据各工具源码）

| 工具 | 什么算"存在" | 基因选择规则 | 样品筛选 |
|---|---|---|---|
| **本工具** | 由集合决定：S / SD / SF / SDF，可加片段长度过滤 | 多个阈值一次输出，⌈t × N⌉ 整数规则，另有 ≥ K 样品 | 完整度门槛、删空样品、贪心删样建议 |
| 旧脚本 `busco_2_matrix_v0.3.py` | 有对应**序列文件**即存在；模式 S、F、SF；**D 从不计入** | 基因覆盖度 = 有序列的样品比例 ≥ 阈值（GPS，`-gpt`）；扫描 50–100 %（步长 5 %），用"基因数 × 平均完整度"的 Pareto 前沿推荐阈值 | 样品完整度 ≥ 阈值（SCH，`-sct`） |
| ZF `BUSCO_extraction.sh`（主版本，BUSCO v5 分支） | single_copy + fragmented 序列（**S+F**） | 先保留 ≥ 3 个类群有序列的位点（脚本统计每个位点的缺失记录数，缺失 < 总类群数 − 2 即保留）；且只认有 `.fna` 的位点；比对后在 `matrix_generation.sh` 中用 PhyKIT 计算的比对占比 `≥ OCCUPANCY %` 筛选 | 无 |
| ZF `version2/`、BUSCO v3 分支 | 仅 single_copy（S） | 同上 | 无 |
| BUSCO_phylogenomics `count_buscos.py` | 仅 S | 统计每个基因出现的样品数 | 无 |

本工具的 `SF` + `--min-taxa` / `-t` 可复现 ZF 主版本的基因集合，`S` / `SF` 可复现旧脚本的 GPS 结果
（旧脚本按序列文件判断，本工具按 full_table 判断；序列文件齐全时二者相同）。

## 验证（论文使用前请看）

[`docs/VALIDATION.md`](docs/VALIDATION.md)：逐版本 BUSCO 输出格式（据源码，[`docs/BUSCO_FORMATS.md`](docs/BUSCO_FORMATS.md)）、
独立代码审查、3,041 个真实样品逐一复算、72 个基因列表独立复算，以及**尚未覆盖的范围**和可直接用的方法描述。

## 与其他工具的比较 (benchmark)

见 [`benchmark/README.md`](benchmark/README.md)：结果由独立脚本从原始 `full_table.tsv` 计算真值后逐项打分，
并记录运行时间和内存。

**基因矩阵交叉验证**（`benchmark/matrix_benchmark.py`，详见 benchmark/README §4）：与 5 个项目的 6 个脚本
（旧脚本 v0.3、BUSCO_phylogenomics、busco2fasta、bigey busco-phylo、ZF BUSCO_extraction.sh 主版本与 version2）
在 D1（199 个基因组）和 D2（77 个样品）上比较。busco2summary 在全部 12 种设置下与独立真值完全一致
（基因集合、每条记录的状态、每条序列）。其他工具按各自规则也都正确，但有两点差异要注意：
busco2fasta 用 `round()`，77 个样品时 "50 %" 只要求 38 个样品，比 ⌈N·t⌉ 多选 21 个基因；
ZF 只收有 `.fna` 的样品，miniprot 结果（BUSCO 6）全部被静默跳过。

## 旧脚本

`legacy/` 保留 `busco_2_matrix_v0.3.py`（内部 v7.2.0）和 `busco2summary_smkv2.py`（v2.8.6，批量运行 BUSCO），
旧行为仍可使用。问题与修复记录见 [`docs/ISSUES.md`](docs/ISSUES.md)。

## 引用

如果使用本工具，请引用本仓库，以及 BUSCO（Manni et al. 2021, *Mol. Biol. Evol.* 38:4647）和
（如用到 CDS 重建）miniprot（Li 2023, *Bioinformatics* 39:btad014）。

---

## English

`busco2summary` turns a folder of BUSCO (v3–v6) or compleasm results into:

* a per-sample summary table (counts, percentages, internal-stop E, assembly statistics, QC flags,
  consistency check between `full_table.tsv` and `short_summary`);
* gene × sample matrices (S/D/F/M status, copy number, length, fragment ratio, binary S/SD/SF/SDF);
* per-BUSCO statistics where **fragmented BUSCOs are counted, measured and flagged rather than dropped**;
* occupancy-based BUSCO-ID lists for any thresholds and status sets (exact rule:
  present in ≥ ceil(t·N) samples), ≥ K-taxa lists, a completeness × occupancy grid and greedy
  taxon-drop advice;
* FASTA extraction (best-scoring copy for duplicates; exact CDS reconstruction for miniprot runs,
  verified by translation against the BUSCO protein);
* an offline interactive HTML report and publication plots.

```bash
pip install .
busco2summary run -i busco_results/ -o b2s_out -j 16
busco2summary genelists -i busco_results/ -o lists -t 50 60 90 --sets S SF
busco2summary extract -i busco_results/ -o seqs --set SF --occupancy 90 --seqtype both
```

Run `busco2summary <command> -h` for all options. Tests: `pytest tests/`.
