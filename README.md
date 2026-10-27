# busco2summary

**多样品 BUSCO 结果整理工具**：把一个文件夹里的 BUSCO / compleasm 结果整理成汇总表、基因 × 样品矩阵、
按占有率筛选的基因 ID 列表和序列，F（片段）作为有用信息保留和统计。

作者：**De-Yuan Yang（杨德援），YZZ Team** · MIT License · Python ≥ 3.8，无必需依赖 · [English](#english)

| [快速上手](#快速上手) | [命令一览](#命令一览) | [安装](#安装) | [输出](#输出) | [整理逻辑](#整理逻辑) | [与其他软件对比](#与其他软件对比) | [文档](#文档) |
|---|---|---|---|---|---|---|

---

## 快速上手

只需要 python3，下载一个 46 KB 的文件，不装任何包：

```bash
curl -LO https://github.com/deyuanyang92-dev/busco2summary/releases/latest/download/b2s.pyz

python3 b2s.pyz matrix -i busco_results/ > status_matrix.tsv          # 基因 × 样品状态矩阵（S/D/F/M）
python3 b2s.pyz genelists -i busco_results/ -o lists -t 50 70 90      # 50/70/90 % 占有率的基因 ID 列表
python3 b2s.pyz run -i busco_results/ -o b2s_out                      # 全部：汇总表、矩阵、列表、HTML 报告
```

`busco_results/` 是放 BUSCO 结果的文件夹，每个样品一个子文件夹即可，层级、BUSCO 版本自动识别。

## 命令一览

| 命令 | 做什么 | 主要输出 |
|---|---|---|
| `busco` | 对 SPAdes / rnaSPAdes / Trinity / MEGAHIT / Flye 组装结果或 FASTA 批量跑 BUSCO，再汇总 | `OUT/<样品名>/`、`busco_summary.tsv` |
| `list` | 列出目录里所有 BUSCO 结果的名字（一行一个） | 名单文件，可编辑后用 `--samples` 读入 |
| `matrix` | 只输出一张矩阵，不建目录、不写缓存 | 一个 TSV（屏幕或 `-o` 文件，`.gz` 自动压缩） |
| `genelists` | 按占有率输出基因 ID 列表 | `gene_lists/S/S_occ090.txt` 等 |
| `run` | 全部分析 | 汇总表、矩阵、逐基因/逐样品统计、列表、图、HTML 报告 |
| `extract` | 按基因列表提取序列 | 按基因 / 按样品的 FASTA |
| `summary` | 只要每个样品的 C/S/D/F/M 汇总 | `busco_summary.tsv` |
| `scan` | 看文件夹里识别到了什么 | 样品、版本、谱系、状态清单 |

**只处理部分样品**：先列出名单，删掉不要的行，再用 `--samples` 读入（所有命令都支持）：

```bash
b2s list -i busco_results/ > samples.txt              # 每行一个样品名（--long 另附文件夹名和路径）
# 编辑 samples.txt，只保留需要的样品
b2s matrix -i busco_results/ --samples samples.txt > status_matrix.tsv
b2s run    -i busco_results/ --samples samples.txt -o b2s_out
```

名单每行可以是样品名、文件夹名或完整路径；`#` 开头的行忽略；名单里找不到的名字会提示。
反过来，`--exclude 名单` 去掉这些样品。

**组装完直接批量跑 BUSCO（SPAdes、rnaSPAdes、Trinity、MEGAHIT、Flye）**：

```bash
b2s busco -i Spades_results/ Trinity_results/ -o Busco_results -l metazoa_odb12 -n   # 先预览：样品名、模式、命令
b2s busco -i Spades_results/ Trinity_results/ -o Busco_results -l metazoa_odb12 -c 8  # 运行，完成后自动汇总
```

| 组装结果 | 用哪个文件 | 样品名 | BUSCO 模式 |
|---|---|---|---|
| `JLJC1_Eun1_1201_spades_rna/`（rnaSPAdes） | `transcripts.fasta` | `JLJC1_Eun1_1201` | transcriptome |
| `23002-1_spades_wgs/`（SPAdes） | `contigs.fasta`（`--fasta scaffolds.fasta` 可改） | `23002-1` | genome |
| `24003-36Z_trinity.Trinity.fasta`（Trinity） | 文件本身 | `24003-36Z` | transcriptome |
| `X_megahit/`、`X_flye/` | `final.contigs.fa`、`assembly.fasta` | `X` | genome |
| `X.fasta`、`X.fna`（任意 FASTA） | 文件本身 | `X` | 按名字判断；`.faa`/`.pep` 为 proteins |

- 样品名 = 组装文件夹名去掉组装软件后缀（`_spades_rna`、`_spades_wgs`、`_spades`、`_trinity`、`_megahit` 等），
  或 FASTA 文件名去掉扩展名，与旧脚本 `busco2summary_smkv2.py -bfsd --suffix_subdir` 相同；
  `--suffix` 后面不跟值则保留原名。去掉后缀后两个样品重名（如 `M4_spades` 和 `M4_spades_wgs`）时，这两个保留原名并提示。
- 结果写到 `Busco_results/<样品名>/`，同时写 `busco_inputs.tsv`（样品、模式、FASTA、状态）、
  `logs/<样品名>.busco.log`，失败的列在 `busco_failed.tsv`；最后自动生成 `busco_summary.tsv`。
- **断点续跑**：已有 `short_summary` 的样品跳过，只跑没完成的（未完成的文件夹用 `-f` 重跑，`--restart` 则续跑）。
- 只跑部分样品：`--samples 名单`；也可以把 `busco_inputs.tsv` 改名、改模式后用 `--list` 读回。
- `-l` 给本地谱系文件夹时自动 `--offline`；`-m` 默认 auto（按上表），也可统一指定。
- 并行：`-c` 每个 BUSCO 的线程数，`-J` 同时跑几个（默认 = 90 % 核心数 ÷ `-c`）。不需要 Snakemake。
- BUSCO 位置：默认 PATH 里的 `busco`，找不到就用 conda 环境 `busco`；`--busco-env BUSCO_v5.8.2` 指定环境，
  不用 `conda activate`。其他 BUSCO 参数：`--busco-args='--miniprot'`。

已经跑完的 BUSCO 只要汇总：`b2s summary -i Busco_results/ -o out`（每个样品一行 C/S/D/F/M、N50 等）。

常用示例：

```bash
# S+F、占有率 ≥ 50 % 的 0/1 矩阵
b2s matrix -i busco_results/ --set SF -t 50 --format binary -o SF_occ050.tsv
# 至少 3 个类群有（S），样品作行
b2s matrix -i busco_results/ --set S --min-taxa 3 --samples-as-rows -o S_min3.tsv
# 给定基因列表，每个样品有哪些（长表：sample  busco_id  status）
b2s matrix -i busco_results/ --genes my_ids.txt --format ids > present_ids.tsv

# S 与 S+F 两套列表并列，便于比较 F 补回多少基因
b2s genelists -i busco_results/ -o lists -t 50 60 70 80 90 --sets S SF
# 指定谱系；多个输入目录；样品表重命名
b2s run -i runA/ runB/ -s samples.tsv --lineage metazoa_odb12 -o out -j 16
# 片段长度 ≥ 完整拷贝一半才算"有"
b2s run -i busco_results/ -o out --min-frag-ratio 0.5

# 提取 S+F、占有率 ≥ 90 % 的蛋白和核酸序列（重复拷贝取得分最高者）
b2s extract -i busco_results/ -o seqs --set SF --occupancy 90 --seqtype both --dup best
# miniprot 结果没有 .fna：用基因组重建 CDS（两列表：样品名<TAB>基因组.fasta）
b2s extract -i busco_results/ -o seqs --genes lists/gene_lists/SF/SF_occ090.txt \
    --seqtype nucleotide --genome-map genomes.tsv
```

`b2s` 与 `busco2summary` 是同一个命令；用 `b2s.pyz` 时把 `b2s` 换成 `python3 b2s.pyz`。
所有参数见 `b2s <命令> -h`。

## 安装

| 方式 | 命令 | 适合 |
|---|---|---|
| 免安装 | `curl -LO https://github.com/deyuanyang92-dev/busco2summary/releases/latest/download/b2s.pyz` | 只想要矩阵 / 列表，不想装任何东西 |
| pip | `pip install git+https://github.com/deyuanyang92-dev/busco2summary.git` | 想要 `b2s` 命令 |
| pip + 图和 Excel | `pip install "busco2summary[all] @ git+https://github.com/deyuanyang92-dev/busco2summary.git"` | 要 PNG/PDF 图、`.xlsx` |
| conda | `git clone …/busco2summary.git && conda env create -f environment.yml` | 在 conda 环境里用 |

资源占用（实测，77 个样品 / 199 个基因组）：

| 用法 | 依赖 | 内存 | 时间 |
|---|---|---|---|
| `matrix` | 仅 Python 标准库 | 32 MB / 58 MB | 3 s / 8 s |
| `run`（无图） | 仅标准库 | 35 MB / 65 MB | 4 s / 9 s |
| `run`（含 PNG/PDF 图） | matplotlib | 约 380 MB | 8 s / 19 s |

**不需要 Snakemake**：一条命令几秒钟跑完。只有当你自己的流程本来就用 Snakemake 时，才可以用自带的 `workflow/Snakefile` 把它接在 BUSCO 后面。

## 输出

```
b2s_out/
├── busco_summary.tsv / .csv (.xlsx)   每个样品一行：C/S/D/F/M、E、组装统计、版本、QC 标记
├── dataset_report.txt                  版本、模式、谱系统计与警告
├── sample_stats.tsv                    进入矩阵的样品：基于 full_table 的统计
├── gene_stats.tsv                      每个 BUSCO：S/D/F/M 数量、占有率、F gain、长度
├── gene_lists/
│   ├── occupancy_summary.tsv           每个 集合 × 阈值 的基因数、矩阵填充率
│   ├── S/ SD/ SF/ SDF/                 <集合>_occ050.txt … <集合>_min3.txt（基因 ID 列表）
│   ├── grid_<集合>.tsv                 样品完整度 × 基因占有率
│   └── taxon_drop_<集合>.tsv           删哪个样品能换回多少基因
├── matrices/                           status / copy_number / length / fragment_ratio / binary_*
├── plots/                              条形图、占有率曲线、片段长度分布、状态热图
├── busco2summary_report.html           离线交互式报告（可直接下载基因 ID）
└── name_mapping.tsv, run_parameters.json, busco2summary.log
```

**状态集合**（决定什么算"有"）：

| 集合 | 算"有"的状态 | 用途 |
|---|---|---|
| `S` | 完整单拷贝 | 经典系统发育基因组学 |
| `SF` | 单拷贝 + 片段 | 用片段补回位点（推荐与 S 对比） |
| `SD` | 单拷贝 + 重复 | 完整性 |
| `SDF` | 任意命中 | 最宽松；转录组 D 多时配合 `--dup best` |

## 整理逻辑

### 基因 ID 怎么来

先汇总所有样品的基因 ID 并去重，再逐个样品填状态
（代码：[`analysis.py` `build_dataset`](src/busco2summary/analysis.py)、[`parsers.py` `parse_full_table`](src/busco2summary/parsers.py)）：

1. **逐个样品读 `full_table.tsv`**，得到 `{BUSCO ID: 状态}`。同一 ID 的多行（重复拷贝、多个 miniprot 命中）
   合并成一个状态，优先级 S > D > F > M。
   缺 full_table 的样品才改为读 `busco_sequences/` 文件名 + `missing_busco_list.tsv`，并标记 `NO_FULL_TABLE`。
2. **所有样品的 ID 合并、去重、排序**，作为矩阵的行（full_table 本来就列出全部 BUSCO，含 Missing，通常即整套谱系基因）。
3. **逐个样品填格子**：有这个 ID 填 S/D/F/M，没有填 M → 状态矩阵。
4. **按状态集合把格子变成 有/无**（`S` 只认 S；`SF` 认 S 和 F）。
5. **逐行数"有"的样品数**，≥ ⌈t × N⌉ 个样品（或 ≥ K 个类群）的基因写入列表，再按列表提取序列。

例子：3 个样品，`SF`，t = 60 % → 需要 ⌈0.6 × 3⌉ = 2 个样品。

| BUSCO ID | X | Y | Z | SF 有几个 | 入选 |
|---|---|---|---|---|---|
| 100at33208 | S | S | S | 3 | ✔ |
| 200at33208 | S | F | M | 2 | ✔ |
| 300at33208 | F | M | M | 1 | ✘ |
| 400at33208 | D | M | M | 0（D 不属于 SF） | ✘ |

**为什么读 full_table 而不是只遍历序列文件名**（ZF、BUSCO_phylogenomics 等的做法）：
full_table 每个样品一个文件，明确记录 M；序列文件名拿不到 M，且没有序列文件的命中会被漏掉
（BUSCO 6 miniprot 结果没有 `.fna`，按 `.fna` 遍历会得到 0 个基因）。序列文件齐全时两种方法结果相同。

### 判定规则（写论文方法时可直接引用）

1. 每个样品每个 BUSCO 只有一个状态，按 full_table 的**不同 ID** 计数，优先级 S > D > F > M；
   compleasm 的 `Single` = S，`Interspaced` = F，`Retrocopy` 归入 D。
2. "存在"由状态集合决定（见上表）；D 只在 `SD`/`SDF` 中算存在。
3. 占有率阈值 t：出现样品数 ≥ ⌈t × N⌉ 才入选（整数规则；N = 过滤后进入矩阵的样品数）；全部缺失的基因不进入任何列表。
4. `--min-frag-ratio R`（默认关闭）：片段长度 ÷ 参考长度 ≥ R 的 F 才算存在；参考长度 = 该 BUSCO 在本数据集完整拷贝的中位长度；
   没有参考长度的 F 不算存在。作用于列表、`binary_*`、序列提取；`status_matrix` 和 S/D/F/M 计数保留 BUSCO 原始判定。
5. 分母 n 取 short_summary 一行式中的 `n:`（BUSCO 自己计算百分比用的数）。
6. BUSCO ≥ 6.1 配 odb12.2 时 full_table 的 Length 含义改为目标蛋白长度，此时不建议用 `--min-frag-ratio`。

### 三个"占比"，别混

| 名称 | 怎么算 | 用来做什么 |
|---|---|---|
| 样品完整度 | 一个样品"有"的 BUSCO 数 ÷ BUSCO 总数 | 剔除差样品（`--min-completeness`） |
| 基因占有率 | 有这个基因的样品数 ÷ 样品数 N | 选基因（`SF_occ090.txt` = ≥ 90 % 样品都有） |
| 矩阵填充率 | 选出的 基因 × 样品 中有数据的格子比例 | 报告最终矩阵缺失程度 |

顺序：**先剔除差样品，再按占有率选基因**——一个很差的样品会让几乎每个基因都"少一个样品"
（`taxon_drop_*.tsv` 告诉你删哪个样品能换回多少基因）。

阈值没有通用最优值：≥ 90 % 缺失少、基因少，适合串联矩阵；50 % 基因多，适合基因树 + ASTRAL。
建议挑 2–3 个阈值（如 50 / 70 / 90 %）都建树，比较拓扑是否稳定。

### 为什么重视 F

多数流程只用 S，F 直接丢弃。在真实数据中这会损失大量位点：

| 数据 | 占有率 | 仅 S | S + F | S+D+F |
|---|---|---|---|---|
| 192 个环节动物基因组（BUSCO 6，metazoa_odb12） | ≥ 80 % | 136 | **521** | 634 |
|  | ≥ 90 % | 1 | **90** | 414 |
| 77 个转录组 + 基因组（BUSCO 5.8） | ≥ 70 % | 46 | **519** | 585 |
|  | ≥ 90 % | 0 | **171** | 371 |

同一数据中，F 片段的中位长度是完整拷贝的 52 %。

## 与其他软件对比

**测试方法**：6 个开源脚本和本工具跑同一份真实数据，结果与一个独立脚本从原始文件算出的"真值"逐项比对
（基因集合、每条 基因 × 样品 记录的状态、每条序列）。

**数据**：metazoa_odb12（672 个 BUSCO）。
- D1：199 个环节动物基因组，BUSCO 6 miniprot
- D2：77 个基因组 + 转录组，BUSCO 5.8

**结果**：

| 工具 | 支持的集合 | 选基因规则 | D1 | D2 |
|---|---|---|---|---|
| **busco2summary** | S、SF、SD、SDF | ≥ ⌈t × N⌉ 个样品；≥ K 个类群 | ✔ 全部一致 | ✔ 全部一致 |
| 旧脚本 busco_2_matrix v0.3 | S、SF | 样品比例 ≥ t | ✔ | ✔ |
| [BUSCO_phylogenomics](https://github.com/jamiemcg/BUSCO_phylogenomics) | S | 样品比例 ≥ t | ✔ | ✔ |
| [bigey busco-phylo](https://github.com/bigey/BUSCO-phylogenetics-pipeline) | S | 样品比例 ≥ t | ✔ | ✔ |
| [busco2fasta](https://github.com/lstevens17/busco2phylo-nf) | S | ≥ round(N × p) 个样品 | ✔ | △ 规则不同：50 % 时多选 21 个基因 |
| [ZF BUSCO_extraction](https://github.com/xtmtd/Phylogenomics) | S（version2）、SF（主版本） | ≥ 3 个类群，只认有 `.fna` 的样品 | ✘ 0 个基因 | ✘ 只用了 32/77 个样品 |

✔ = 基因集合、每条记录状态、每条序列都与真值一致；△ = 按自己的规则正确，但与常用规则结果不同；✘ = 漏掉样品或基因（D1、D2 各测 S 与 SF、50 % / 90 % / ≥ 3 类群，共 12 种设置）。

**差异原因**：
- **busco2fasta**：用四舍五入。77 个样品时 "50 %" = round(38.5) = 38 个样品，正确应为 ⌈77 × 0.5⌉ = 39，因此多选 21 个基因。
- **ZF**：BUSCO 6 miniprot 结果没有 `.fna`，这些样品被跳过，而且没有任何提示。
- **其他**：按各自规则都正确。所有软件都没有改动输入文件。

**鲁棒性**（8 种不同目录结构混在一个文件夹）：

| busco2summary | 旧脚本 v0.3 | MultiQC | BUSCO_phylogenomics |
|---|---|---|---|
| **8/8** | 6/8 | 5/8 | 0/8（崩溃） |

<details>
<summary>功能对比（据各工具源码）</summary>

| 功能 | **busco2summary** | 旧脚本 v0.3 | BUSCO_phylogenomics | busco2fasta | ZF | MultiQC |
|---|---|---|---|---|---|---|
| 每个样品 S/D/F/M 计数 | ✓ | ✓ | 只有 S | ✗ | ✗ | ✓ |
| 基因 × 样品状态矩阵（含全缺失基因） | ✓ | ~ | 只有 S | ✗ | ✗ | ✗ |
| 使用 F（S+F）及片段长度比 | ✓ | ~（无长度比） | ✗ | ✗ | ~ 主版本 S+F | ✗ |
| 一次输出多个占有率阈值 | ✓ | ~ | 只有计数 | 一个阈值 | 固定 ≥ 3 类群 | ✗ |
| 完整度 × 占有率网格、删样品建议 | ✓ | ~ | ✗ | ✗ | ✗ | ✗ |
| 重复拷贝取得分最高者 | ✓ | ✗ | ✗ | ✗ | ✗ | ✗ |
| miniprot 结果（无 `.fna`）重建 CDS | ✓ | ✗ | ✗ | ✗ | ✗ | ✗ |
| BUSCO v3 / v4–v6 / compleasm | ✓ / ✓ / ✓ | ~ / ✓ / ✗ | ~ / ✓ / ✗ | ✗ / ✓ / ✗ | ✓ / v5 / ✗ | ✓ / ✓ / ✗ |
| 不删除、不修改原始数据 | ✓ | `--overwrite` 删输出目录 | ✓ | 总删输出目录 | 改输入副本 | `-f` 覆盖 |
| 交互式 HTML / 断点续跑 | ✓ / ✓ | ✓ / ✗ | ✗ / ✗ | ✗ / ✗ | ✗ / ✗ | ✓ / ✗ |

</details>

<details>
<summary>速度与内存（199 个基因组）</summary>

| 工具 | 做了什么 | 时间 / 内存 |
|---|---|---|
| **busco2summary `matrix`** | 一张矩阵 | 8 s / 58 MB |
| **busco2summary `run`** | 全部（含图和 HTML） | 18.8 s / 365 MB；再次运行（缓存）2.7 s |
| 旧脚本 busco_2_matrix v0.3 | 矩阵 + 列表 + HTML | 24.8 s / 234 MB |
| BUSCO_phylogenomics `count_buscos.py` | 只统计 S | 1.1 s / 18 MB |
| MultiQC busco | 每个样品的百分比 | 29.2 s / 320 MB |

</details>

完整结果和复现命令：[`benchmark/README.md`](benchmark/README.md)。

## 文档

| 文件 | 内容 |
|---|---|
| [`docs/VALIDATION.md`](docs/VALIDATION.md) | 验证：3,041 个真实样品逐一复算、72 个基因列表独立复算、尚未覆盖的范围、可直接用的方法描述 |
| [`docs/BUSCO_FORMATS.md`](docs/BUSCO_FORMATS.md) | 各版本 BUSCO 输出格式（据源码） |
| [`benchmark/README.md`](benchmark/README.md) | 与其他软件对比的完整表格与复现 |
| [`docs/ISSUES.md`](docs/ISSUES.md) | 问题与修复记录 |
| [`CHANGELOG.md`](CHANGELOG.md) | 版本更新 |
| `legacy/` | 旧脚本 `busco_2_matrix_v0.3.py`、`busco2summary_smkv2.py`（Snakemake 批量运行 BUSCO；新命令 `busco` 沿用它的命名规则），旧行为仍可用 |

## 引用

如果使用本工具，请引用本仓库（[`CITATION.cff`](CITATION.cff)），以及 BUSCO（Manni et al. 2021, *Mol. Biol. Evol.* 38:4647）和
（如用到 CDS 重建）miniprot（Li 2023, *Bioinformatics* 39:btad014）。

---

## English

`busco2summary` turns a folder of BUSCO (v3–v6) or compleasm results into a per-sample summary, gene × sample
matrices (S/D/F/M, copy number, length, binary S/SD/SF/SDF), occupancy-based BUSCO-ID lists (exact rule:
present in ≥ ⌈t·N⌉ samples; also ≥ K taxa), FASTA extraction (best-scoring copy for duplicates; exact CDS
reconstruction for miniprot runs) and an offline HTML report. **Fragmented BUSCOs are counted, measured and
flagged rather than dropped.** Gene IDs and statuses come from `full_table.tsv` (one status per BUSCO ID,
S > D > F > M), not from sequence file names.

No required dependencies; zero-install single file (Python ≥ 3.8):

```bash
curl -LO https://github.com/deyuanyang92-dev/busco2summary/releases/latest/download/b2s.pyz
python3 b2s.pyz matrix -i busco_results/ > status_matrix.tsv
python3 b2s.pyz matrix -i busco_results/ --set SF -t 50 --format binary -o SF_occ050.tsv
python3 b2s.pyz run -i busco_results/ -o b2s_out
# run BUSCO on many SPAdes / rnaSPAdes / Trinity / MEGAHIT / Flye assemblies, then summarise
# (sample = folder name minus assembler suffix; finished samples skipped; -n = preview)
python3 b2s.pyz busco -i Spades_results/ Trinity_results/ -o Busco_results -l metazoa_odb12 -c 8
```

Or `pip install git+https://github.com/deyuanyang92-dev/busco2summary.git` (add `[all]` for plots and xlsx).
Cross-validated against six scripts from five projects on 276 real BUSCO runs: identical gene sets, statuses and
sequences in all 12 settings (see the comparison table above and `benchmark/README.md`).
Author: De-Yuan Yang (YZZ Team). Run `b2s <command> -h` for all options.
