# GLADE 用户使用说明

本文按常规操作顺序介绍 GLADE：从目标化合物出发，准备路线、蛋白和表达构建，最后生成理论质粒设计。

- [运行准备](#1-运行准备)
- [常规完整流程](#2-常规完整流程)
- [可选功能](#3-可选功能)
- [结果与重新运行](#4-结果与重新运行)
- [常见问题](#5-常见问题)

示例使用项目目录 `F:\myproject\glade`、目标 `C00811` 和配置文件 `demo01.json`。
路线、方案和蛋白编号都应替换为自己结果中的编号。

## 1. 运行准备

### 启动环境

项目使用 Python 3.12。在 PowerShell 中进入项目并激活已有虚拟环境：

```powershell
cd F:\myproject\glade
& .\.venv\Scripts\Activate.ps1
```

如果激活脚本被阻止，可先执行：

```powershell
Set-ExecutionPolicy -Scope Process -ExecutionPolicy RemoteSigned
```

首次安装、网络服务和 RetroPath 部署见 [部署教程](GLADE部署教程.md)。

### 创建输入配置

创建 `inputs/demo01.json`：

```json
{
  "target_name": "C00811"
}
```

`target_name` 使用大写 KEGG Compound ID，格式为 `C` 加五位数字。

所有命令的 `-i` 后只填文件名，例如 `-i demo01.json`，程序会自动查找 `inputs`。
当前默认底盘是 *E. coli* MG1655，使用 iML1515 模型。

## 2. 常规完整流程

操作顺序：

```text
准备配置 → 底盘分析 → 搜索并选择路线 → 选择主酶
→ 可选添加辅助蛋白 → 生成 CDS → 按需运行 DNA Chisel
→ 选择表达盒分组 → 上传元件或确认推荐元件
→ 构建完整表达序列 → 选择质粒骨架 → 生成并执行组装计划
```

### 2.1 分析底盘

分析当前底盘和培养基可提供的代谢物：

```powershell
python main.py chassis -i demo01.json
```

完成后继续搜索路线。若提示目标已可直接供给，本次模型条件下无需新增合成路线。

### 2.2 搜索并选择路线

先搜索，再查看路线列表和所选路线详情：

```powershell
python main.py gap -i demo01.json -d 0
python main.py info -i demo01.json --gap -d 0
python main.py info -i demo01.json --solution 1 -d 0
```

确认路线后写入设计清单：

```powershell
python main.py write -i demo01.json --solution 1 -d 0
```

`1` 是示例路线编号。需要通量验证时，在写入前参考 [路线验证](#32-验证路线)；
验证是可选步骤。没有合适路线时，可尝试 [底盘扩展](#31-扩展底盘可达集合) 或 [RetroPath](#33-使用-retropath)。

**同一批路线的搜索、查看、验证和写入必须使用相同的 `-d`。** `-d 0` 表示原始底盘集合。

### 2.3 选择主酶

生成各反应的主酶候选，并查看结果：

```powershell
python main.py main-enzyme -i demo01.json
python main.py info -i demo01.json --main-enzyme-candidates
```

生成组合，查看后选择一个组合：

```powershell
python main.py main-enzyme-sets -i demo01.json
python main.py info -i demo01.json --main-enzyme-sets
python main.py write -i demo01.json --main-enzyme-set 1
```

需要辅助蛋白时，先按 [辅助蛋白](#34-添加辅助蛋白) 添加，再生成 CDS。

### 2.4 生成并查看 CDS

根据已选蛋白生成 CDS：

```powershell
python main.py protein-to-cds -i demo01.json
python main.py info -i demo01.json --cds
```

系统保存两份序列：

| 序列 | 用途 |
|---|---|
| `raw_cds` | 保留生成或上传的原始序列，供查看 |
| `optimized_cds` | 当前工作副本，所有后续编辑都作用于它 |

初始两份序列相同。`protein-to-cds` 不自动运行 DNA Chisel。
**再次成功执行该命令会覆盖 optimized，即使命中生成缓存，也会重置已有编辑。**

查看原始指标：

```powershell
python main.py info -i demo01.json --cds --raw
```

### 2.5 按需编辑 CDS

DNA Chisel 支持整体 GC、局部 GC、指定限制酶位点和同聚物消除。按需要分别执行：

| 功能 | 命令示例 |
|---|---|
| 单条 CDS 整体 GC | `python main.py optimize -i demo01.json --cds P21683 --gc-min 40 --gc-max 60` |
| 单条 CDS 局部 GC | `python main.py optimize -i demo01.json --cds P21683 --gc-min 30 --gc-max 70 --window 50` |
| 全部 CDS 禁止酶位点消除 | `python main.py optimize -i demo01.json --enzyme EcoRI HindIII` |
| 全部 CDS 同聚物消除 | `python main.py optimize -i demo01.json --homopolymer-max 6` |

- 蛋白编号从 `info --proteins` 查看；表中的范围和阈值都是示例，需要自行指定。
- GC 单位为百分比；`--window` 是连续窗口长度，必须为正整数且不超过 CDS 长度。
- 酶名忽略大小写；没有选择酶时，不启用默认禁止酶。`--enzyme none` 清空酶选择，不恢复原序列。
- 同聚物是连续相同碱基；最大长度 6 表示允许 6 个，7 个及以上需要打断。
- 三种编辑模式分开执行，后续编辑仍需满足已经设置的其他约束。

编辑保持编码蛋白、长度及原始起止密码子不变。要求无法满足时，本次编辑不保存。
不需要这些调整时，可以直接进入表达盒分组。

### 2.6 选择表达盒分组

分组决定哪些蛋白放在同一个表达盒，以及 CDS 的排列顺序：

```powershell
python main.py expression --design --box -i demo01.json
python main.py write -i demo01.json --expression-box 1
python main.py info -i demo01.json --expression-box
```

此时只确定分组，启动子、RBS 和终止子在下一步准备。也可使用 [自定义分组](#35-自定义表达盒分组)。

### 2.7 上传表达元件

将元件文件放在 `inputs/parts/`。每个表达盒需要一个启动子和终止子，每个蛋白需要一个 RBS。

下面示范为表达盒 1 和蛋白 `P21683` 上传元件：

```powershell
python main.py expression -i demo01.json --promoter 1 promoter_1.txt
python main.py expression -i demo01.json --rbs P21683 rbs_1.txt
python main.py expression -i demo01.json --terminator 1 terminator_1.txt
```

启动子和终止子按**表达盒编号**定位，RBS 按**蛋白编号**定位。
对其余表达盒和蛋白重复上传，直到所有元件齐全。

文件要求：

- 参数只填 `inputs/parts` 下的文件名。
- 支持 `.txt`、`.fa`、`.fasta`、`.fna`；TXT 放 DNA 序列，FASTA 只放一条记录。
- 序列只能包含 A、C、G、T，空白和大小写会自动处理。
- 同一位置再次上传会替换该元件；RBS 上传时会计算翻译起始率，预测失败则不保存。

上传后查看缺项、预测结果和序列冲突：

```powershell
python main.py info -i demo01.json --expression-box
```

若希望系统推荐元件，使用 [元件推荐](#36-使用系统推荐元件) 替代本节上传步骤。

### 2.8 构建完整表达序列

元件上传齐全，或确认推荐方案后，都执行同一个命令：

```powershell
python main.py expression -i demo01.json --assemble
```

每个表达盒按下面的顺序拼接，同一方案中的全部表达盒再按编号串联：

```text
启动子 + RBS₁ + CDS₁ + RBS₂ + CDS₂ + … + 终止子
```

系统使用当前 optimized CDS，刷新失效的 RBS 预测，检查完整序列及所有连接处的禁止酶位点和同聚物。
完整构建的整体与局部 GC 只统计。缺少元件、预测失败或存在序列冲突时，不登记新构建。

成功后导出 `expression_constructs/design_001.gb`，并登记到 manifest。
多个已选元件方案分别生成构建文件。此时尚未添加接入质粒的克隆末端。

### 2.9 选择质粒并完成理论组装

推荐骨架，根据结果选择一个：

```powershell
python main.py plasmid --recommend -i demo01.json
python main.py write -i demo01.json --plasmid 1
```

生成完整组装计划，查看终端结果，确认后接受并执行：

```powershell
python main.py assembly --plan -i demo01.json
python main.py write -i demo01.json --assembly-plan
python main.py assembly --execute -i demo01.json
```

系统自动选择可行的 Gibson 或双酶切方法。计划必须完整才能接受；
需要指定方法或抗性时，见 [可选参数速查](GLADE部署教程.md#8-可选用法速查)。

主要交付结果位于 `final_assembly/`：

- `design_001_final.gb`：带注释的完整质粒。
- `design_001_final.fasta`：质粒 DNA 序列。
- `final_design_report_zh.md`：中文设计报告。

这些是计算设计结果，仍需实验验证。

## 3. 可选功能

### 3.1 扩展底盘可达集合

没有合适路线时，可先扩展底盘集合，再使用同一深度搜索：

```powershell
python main.py expand -i demo01.json -d 1
python main.py gap -i demo01.json -d 1
python main.py info -i demo01.json --gap -d 1
```

扩展深度必须大于等于 1。扩展结果表示补充相应反应和辅助系统后可达，不表示底盘天然具备这些能力。

### 3.2 验证路线

验证路线在当前代谢模型和培养基中的通量可行性：

```powershell
python main.py validate -i demo01.json -s 1 -m per -c strict -d 0
```

`-s` 指定路线编号。验证失败或未验证的路线仍可写入，但需要人工复核。
多路线验证见 [验证速查](GLADE部署教程.md#8-可选用法速查)，RetroPath 验证数据准备见 [部署教程](GLADE部署教程.md#33-mnxref)。

### 3.3 使用 RetroPath

先按 [部署教程](GLADE部署教程.md#32-retropath-规则与服务) 启动本地服务，再运行：

```powershell
python main.py gap -i demo01.json --retropath -d 0
python main.py info -i demo01.json --gap -d 0
python main.py info -i demo01.json --retropath-candidate 1 -d 0
```

候选详情会给出 GLADE 路线编号。后续查看、验证和 `write --solution` 使用这个编号，
不要直接使用候选排名。默认预测最多 3 步，可用 `--step 5` 等参数调整。

### 3.4 添加辅助蛋白

在选择主酶组合后、生成 CDS 前添加。将序列文件放在 `inputs`：

```powershell
python main.py add-auxiliary-protein -i demo01.json --protein-file helper.fasta --sequence-type protein
```

上传的是 CDS 时，将类型改为 `cds`。这类序列不经过 CodonTransformer，
后续仍会保存 raw 和 optimized；DNA Chisel 只编辑 optimized。

支持 FASTA、FAA 和纯文本；FASTA ID 或纯文本文件名作为蛋白编号。
同 ID 再次导入会替换已有内容，导入后无需再执行 `write`。

查看当前蛋白：

```powershell
python main.py info -i demo01.json --proteins
```

也可让系统研究辅助蛋白，或删除已上传的辅助蛋白，命令见 [可选用法速查](GLADE部署教程.md#8-可选用法速查)。

### 3.5 自定义表达盒分组

每个方括号表示一个表达盒，蛋白顺序就是组装顺序：

```powershell
python main.py expression --design --box -i demo01.json --custom [P00001 P00002] [P00003]
```

替换为自己的蛋白编号，当前全部蛋白必须恰好出现一次。
该命令直接登记分组，无需再执行 `write --expression-box`。

### 3.6 使用系统推荐元件

用下面两条命令替代手动上传，之后回到 [统一构建](#28-构建完整表达序列)：

```powershell
python main.py expression --design --parts -i demo01.json
python main.py write -i demo01.json --expression-parts 1
```

默认请求 12 个方案；`write` 只确认选择，不立即生成 GenBank。
可选择多个方案，例如 `--expression-parts 1:3`，区间包含两端。

重新推荐不会覆盖已确认的元件。上传元件缺少活性数据时，表达负担使用有记录的中性估算；
推荐评分和负担估算都不等同于实验表现。

## 4. 结果与重新运行

### 常用结果位置

结果统一保存在 `outputs/C00811/`，其他目标使用自己的编号：

| 内容 | 相对该目录的位置 |
|---|---|
| 底盘可提供的代谢物 | `chassis_result/producible_kegg_compounds.csv` |
| 路线列表 | `kegg_gap_C00811/depth0/solutions.csv` |
| 原始 CDS | `protein_to_cds/raw_cds/` |
| 当前 CDS | `protein_to_cds/optimized_cds/` |
| 完整表达构建 | `expression_constructs/` |
| 最终质粒与设计报告 | `final_assembly/` |
| 项目状态和文件引用 | `design_manifest.json` |

日常通过 `info` 查看结果即可，不需要手动编辑 manifest 或报告文件。

### 常用查看命令

| 查看内容 | 命令 |
|---|---|
| 底盘分析 | `python main.py info -i demo01.json --chassis` |
| 所有蛋白 | `python main.py info -i demo01.json --proteins` |
| 单个蛋白 | `python main.py info -i demo01.json --protein P21683` |
| 当前 CDS 指标 | `python main.py info -i demo01.json --cds` |
| 原始 CDS 指标 | `python main.py info -i demo01.json --cds --raw` |
| 表达盒、元件与构建状态 | `python main.py info -i demo01.json --expression-box` |
| 已确认元件方案 2 | `python main.py info -i demo01.json --expression-box --parts-design 2` |

### 上游修改后从哪里继续

| 修改内容 | 后续操作 |
|---|---|
| 更换路线 | 重新选择主酶，并执行后续流程 |
| 更换蛋白名单或添加、删除辅助蛋白 | 重新生成 CDS、确认分组并准备元件，再构建和组装 |
| 编辑或重新生成 CDS，蛋白名单和分组有效 | 保留分组和元件；重新 `expression --assemble`，再选择质粒并组装 |
| 修改表达盒分组 | 重新配置元件，再构建、选择质粒并组装 |
| 上传或更换元件、确认不同元件方案 | 重新构建、选择质粒并组装 |
| 更换质粒骨架 | 重新生成、接受并执行组装计划 |

相同有效输入通常会复用结果；**`protein-to-cds` 会覆盖 optimized**。
旧预测和构建失效后，应重新执行对应步骤。

## 5. 常见问题

| 问题 | 处理方法 |
|---|---|
| 找不到配置文件 | 确认文件在 `inputs`，`-i` 后只填文件名 |
| 目标格式错误 | 使用大写 `C` 加五位数字，例如 `C00811` |
| 没有候选路线 | 检查目标与底盘结果，尝试扩展或 RetroPath |
| 提示结果已过期 | 使用相同深度，从提示指定的上游步骤重新运行 |
| 主酶检索或推荐服务失败 | 检查网络及服务配置，见进阶说明 |
| DNA Chisel 要求无法满足 | 根据提示调整范围、窗口或约束组合；原工作序列保留 |
| RBS 上传预测失败 | 检查 RBS 与 CDS 是否合适、CDS 是否完整，处理后重新上传 |
| 完整构建发现位点或同聚物冲突 | 根据位置和元件名称处理对应输入，再执行统一构建 |
| `protein-to-cds` 退出码为 2 | 查看 `protein_to_cds/run_summary.json` 中的失败原因 |

其他参数可通过命令帮助查看：

```powershell
python main.py --help
python main.py optimize --help
```

部署与服务配置见 [部署教程](GLADE部署教程.md)，详细验证和更多参数见其中的 [可选用法速查](GLADE部署教程.md#8-可选用法速查)。
