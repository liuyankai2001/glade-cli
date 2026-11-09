# GLADE 用户使用说明

GLADE 用于设计微生物合成目标化合物的候选路线、酶蛋白、表达序列和质粒。结果需经实验验证。

首次安装见 [部署教程](GLADE部署教程.md)。下文以 `C00811` 为例，所有命令在项目根目录执行，方案和蛋白编号请按实际结果替换。

## 1. 运行准备

打开 PowerShell，进入项目并激活环境：

```powershell
cd F:\myproject\glade
& .\.venv\Scripts\Activate.ps1
```

若激活被阻止，先执行 `Set-ExecutionPolicy -Scope Process -ExecutionPolicy RemoteSigned`。使用 uv 时，也可将后文的 `python main.py` 换成 `uv run python main.py`，无需激活环境。

创建 UTF-8 文件 `inputs/demo01.json`：

```json
{
  "target_name": "C00811"
}
```

- `target_name` 填 KEGG Compound 编号（大写 `C` 加五位数字）。
- **`-i` / `--input` 后只填 `inputs` 下的文件名**，如 `demo01.json`。
- 默认底盘为 *E. coli* MG1655，模型为 `data/gem_models/iML1515.json`，培养基为 `data/mediums/default_medium.json`。
- `info` 查看结果，`write` 确认选择；查看候选后需确认，后续流程才会使用。

## 2. 常规完整流程

底盘分析 → 选路线 → 选主酶 → 生成 CDS → 配置表达盒 → 构建表达序列 → 设计质粒。

### 2.1 分析底盘

```powershell
python main.py chassis -i demo01.json
python main.py info -i demo01.json --chassis
```

结果表示当前模型和培养基下可供给的代谢物。若目标已可供给，无需新增合成路线；这不代表已达到所需产量。

### 2.2 搜索并选择路线

```powershell
python main.py gap -i demo01.json -d 0
python main.py info -i demo01.json --gap -d 0
python main.py info -i demo01.json --solution 1 -d 0
python main.py write -i demo01.json --solution 1 -d 0
```

选择前检查前体、反应方向和辅因子需求。**搜索、查看、验证和确认同一批路线时，`-d` 必须一致**；它表示底盘扩展深度，不是路线步数。无合适路线时见第 3 节。

### 2.3 选择主酶

```powershell
python main.py main-enzyme -i demo01.json
python main.py info -i demo01.json --main-enzyme-candidates
python main.py main-enzyme-sets -i demo01.json
python main.py info -i demo01.json --main-enzyme-sets
python main.py write -i demo01.json --main-enzyme-set 1
python main.py info -i demo01.json --proteins
```

检查反应匹配、蛋白来源及待复核提示，再确认组合。选中蛋白的 FASTA 和注释 JSON 保存在 `inputs/protein_sequences/`，缺失时联网补齐。人工选酶见第 3.4 节；若全部反应已选齐并提示 `already_selected`，无需再次写入组合。

### 2.4 生成并调整 CDS

```powershell
python main.py protein-to-cds -i demo01.json
python main.py info -i demo01.json --cds
```

CDS 是蛋白质的 DNA 编码序列。`raw_cds` 保存初始序列，`optimized_cds` 保存当前用于构建的序列；初始两者相同。没有 CUDA 时可加 `--device cpu`。

按需分别执行以下调整，保持编码蛋白不变：

| 目的 | 命令示例 |
|---|---|
| 单条 CDS 整体 GC 为 40%～60% | `python main.py protein optimize -i demo01.json --cds P21683 --gc-min 40 --gc-max 60` |
| 每个 50 nt 窗口 GC 为 30%～70% | `python main.py protein optimize -i demo01.json --cds P21683 --gc-min 30 --gc-max 70 --window 50` |
| 去除全部 CDS 内的指定酶切位点 | `python main.py protein optimize -i demo01.json --enzyme EcoRI HindIII` |
| 连续相同碱基最多 6 个 | `python main.py protein optimize -i demo01.json --homopolymer-max 6` |

GC 调整必须用 `--cds` 指定蛋白编号；酶切位点和连续重复碱基调整批量作用于全部 CDS。三类调整分开执行，后续调整仍需满足之前的约束。冲突时本次调整不保存，按提示处理。

**重新运行 `protein-to-cds` 会重置当前 CDS，之前的调整需重做。** 查看初始指标用 `info --cds --raw`；查看当前指标用 `info --cds`。

### 2.5 选择表达盒分组

每个表达盒共用一个启动子和终止子，每条 CDS 配一个 RBS（核糖体结合位点）。

```powershell
python main.py expression --design --box -i demo01.json
python main.py write -i demo01.json --expression-box 1
python main.py info -i demo01.json --expression-box
```

### 2.6 配置表达元件

将启动子、RBS 和终止子序列放入 `inputs/parts/`，逐项上传：

```powershell
python main.py expression -i demo01.json --promoter 1 promoter_1.txt
python main.py expression -i demo01.json --rbs P21683 rbs_1.txt
python main.py expression -i demo01.json --terminator 1 terminator_1.txt
python main.py info -i demo01.json --expression-box
```

启动子和终止子填**表达盒编号**，RBS 填**蛋白编号**。对所有表达盒和蛋白完成配置。

文件名不带路径；支持 TXT 或单条 FASTA（`.fa`、`.fasta`、`.fna`），序列仅含 A/C/G/T。同位置重复上传会替换元件；RBS 预测失败时此次上传不保存。使用系统推荐可替代上传，见第 3.5 节。

### 2.7 构建完整表达序列

```powershell
python main.py expression -i demo01.json --assemble
```

程序按“启动子 + RBS₁ + CDS₁ + … + 终止子”拼接各表达盒，导出 `expression_constructs/design_001.gb` 等文件。缺少元件、RBS 预测失败或连接处存在序列冲突时，按提示处理后重跑。此时尚未得到完整质粒。

### 2.8 设计质粒并导出

选择下面一种方式。

**网页搭配模块**：需仅确认一个完整表达构建，并指定两种不同的常规黏性末端酶。

```powershell
python main.py protein optimize -i demo01.json --enzyme EcoRI HindIII
python main.py expression -i demo01.json --assemble
python main.py web -i demo01.json
```

在网页中选择抗性、复制模块及可选终止子/间隔区段，拖动调整顺序，点击“生成设计”后导出。注意：

- 复制模块只保留一个；抗性、终止子和间隔区段可重复添加，每个副本都会进入最终序列。
- T0/T1 需手动选择，参考顺序为 `T1 → 完整表达构建 → T0`；程序不会自动补入。
- 内部酶切冲突需更换酶或模块；更换酶后重新检查 CDS 并构建表达序列。宿主信息仅作记录，需自行确认模块兼容性。

结果在 `plasmid_designs/`。使用完按 Ctrl+C 停止；端口占用可加 `--port 8001`。

**命令行推荐骨架**：推荐并选择骨架，再生成、接受和执行组装计划。

```powershell
python main.py plasmid --recommend -i demo01.json
python main.py write -i demo01.json --plasmid 1
python main.py assembly --plan -i demo01.json
python main.py write -i demo01.json --assembly-plan
python main.py assembly --execute -i demo01.json
```

系统选择可行的 Gibson 或双酶切方法，计划完整后才能接受。`final_assembly/` 包含完整质粒 `.gb`、`.fasta` 和 `final_design_report_zh.md`。`--execute` 生成理论序列与文件，实验前仍需复核。

## 3. 可选功能

### 3.1 扩展底盘可达集合

```powershell
python main.py expand -i demo01.json -d 1
python main.py gap -i demo01.json -d 1
python main.py info -i demo01.json --gap -d 1
```

`-d` 为扩展层数，至少为 1。新增前体依赖额外反应，相关改造也需纳入设计；后续查看、验证和确认继续用同一深度。

### 3.2 验证路线

```powershell
python main.py validate -i demo01.json -s 1 -m per -c strict -d 0
```

`-s` 可填多个路线编号，省略则验证全部；`-c relaxed` 放宽辅因子以诊断阻断。纯 KEGG 路线还支持 `-m pooled`（联合验证）和 `-m both`（独立加联合）。RetroPath 仅用 `-m per`，需先安装 [MNXref](GLADE部署教程.md#33-mnxref)。

验证通过表示当前模型约束下通量可行。未通过或未验证仍可选择，但需复核。纯 KEGG 已选路线的匹配单路线验证结果会自动回写清单，无需再 `write`；同一路线保留最新结果，relaxed 通过仅适用于放宽条件。

### 3.3 使用 RetroPath

先按 [部署教程](GLADE部署教程.md#32-retropath-规则与服务) 启动服务：

```powershell
python main.py gap -i demo01.json --retropath -d 0
python main.py info -i demo01.json --gap -d 0
python main.py info -i demo01.json --retropath-candidate 1 -d 0
```

后续使用候选详情中的 **GLADE 路线编号**，不要直接使用 RetroPath 候选排名。默认最多预测 3 步，搜索时可加 `--step 5`；它与 `-d` 含义不同。

### 3.4 人工管理蛋白

路线写入后，可为反应添加、删除或替换蛋白：

```powershell
python main.py protein add -i demo01.json -r R02061 -f MY_ENZYME.fasta
python main.py protein remove -i demo01.json -r R02061 -p MY_ENZYME
python main.py protein replace -i demo01.json -r R02061 -p P22873 -f MY_ENZYME.fasta
```

`-r` 必须属于当前路线；`-f` 只填 `inputs/protein_sequences/` 下文件名，文件须为单条标准氨基酸序列，编号取文件名去掉扩展名后的大写形式。添加和替换不联网，缺少注释 JSON 时补建；同编号序列冲突时报错。删除只移除指定反应下的关联，保留序列文件，也不永久排除该蛋白。

含人工蛋白的反应会固定整组已选蛋白，跳过自动检索。修改后按终端提示重新生成候选或组合；已有 CDS、表达盒和质粒需重新构建。

也可单独下载蛋白，下载不等于选用：

```powershell
python main.py protein download -i demo01.json --uniprot P22873
python main.py protein download -i demo01.json --name TBT
```

编号与名称二选一；按名称时优先选择与底盘分类接近的有效候选。文件保存在 `inputs/protein_sequences/`，强制刷新加 `--force`。详见 [蛋白下载说明](protein_download.md)。

### 3.5 自定义分组与推荐元件

自定义分组中，每个方括号表示一个表达盒；全部已选蛋白须恰好出现一次，顺序即 CDS 顺序。此命令直接保存，无需再 `write`：

```powershell
python main.py expression --design --box -i demo01.json --custom [P00001 P00002] [P00003]
```

已配置元件库时，可用推荐代替手动上传：

```powershell
python main.py expression --design --parts -i demo01.json
python main.py write -i demo01.json --expression-parts 1
```

之后仍需 `expression --assemble`。可用 `--expression-parts 1:3` 选择多个方案；网页质粒设计仅支持一个。重新推荐后需重新确认方案并构建。

### 3.6 更换宿主

在输入 JSON 中指定模型、宿主 GBFF 和 CDS 生成类别，例如：

```json
{
  "target_name": "C00811",
  "model_file": "inputs/model.json",
  "genome_file": "data/genomes/my_strain.gbff",
  "codon_transformer_organism": 51
}
```

JSON 内路径相对项目根目录，也支持绝对路径。`genome_file` 用于识别宿主和过滤已有 EC；未填写时默认 MG1655，关闭 EC 过滤。它与模型重建使用的 `bactabolize.target_genome_file` 独立。

`codon_transformer_organism` 填 `data/codontransformer_organisms.csv` 中的整数编号：`51` 为大肠杆菌通用类别，`52` 为 MG1655；这不是 Taxonomy ID。MG1655 默认用 `52`，其他宿主需明确配置。该选项不改变真实宿主，元件库仍需支持该宿主。修改类别后重新生成 CDS 及下游设计。

## 4. 结果与重新运行

默认结果位于 `outputs/C00811/`：

| 内容 | 相对路径 |
|---|---|
| 底盘可供给代谢物 | `chassis_result/producible_kegg_compounds.csv` |
| 路线列表 | `kegg_gap_C00811/depth0/solutions.csv` |
| 初始 / 当前 CDS | `protein_to_cds/raw_cds/` / `protein_to_cds/optimized_cds/` |
| 完整表达序列 | `expression_constructs/` |
| 命令行 / 网页质粒结果 | `final_assembly/` / `plasmid_designs/` |
| 已选方案与项目状态 | `design_manifest.json`（无需手动编辑） |

结果按目标编号归档，同目标配置可能共用结果目录，保留不同设计前先备份。自定义模型结果位于 `outputs/<目标>/model_<模型摘要前16位>/`；模型内容改变后，需重新从 `chassis` 开始。

| 修改内容 | 从哪里继续 |
|---|---|
| 路线或主酶 | 重新完成选酶、CDS 及下游流程 |
| CDS，且蛋白与分组仍有效 | 重建表达序列，再设计质粒并组装 |
| 表达盒分组或元件 | 补齐元件、重建表达序列，再设计质粒并组装 |
| 质粒骨架 | 重新生成、接受并执行组装计划 |

重新计算会使受影响的下游结果失效。遇到格式过期或结果失效提示，从指定阶段重跑，勿继续使用旧质粒文件。

## 5. 常见问题

| 问题 | 处理 |
|---|---|
| 找不到输入文件 | 确认文件在 `inputs`，扩展名为 `.json`，`-i` 只传文件名 |
| 没有路线 | 检查目标编号、底盘和培养基，尝试扩展或 RetroPath |
| 检索或推荐失败 | 检查网络和 [服务配置](GLADE部署教程.md#23-外部服务) |
| CDS 调整或表达拼接失败 | 检查 GC、酶切位点、重复碱基及连接处冲突，处理后重跑 |
| 网页项目未准备好 | 确认只有一个完整表达构建，已指定两种适用酶并重新构建 |
| `protein-to-cds` 退出码 2 | 查看 `protein_to_cds/run_summary.json` 中的失败原因 |

更多参数用 `python main.py <命令> --help` 查看，例如 `python main.py protein optimize --help`。
