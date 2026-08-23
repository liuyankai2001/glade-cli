# GLADE 用户使用说明

本文档说明当前版本 GLADE 已实现的命令行工作流。示例均假定项目位于
`F:\myproject\glade`，目标化合物为 KEGG Compound ID `C00811`，输入配置文件名为
`demo01.json`。

## 1. 工作流概览

GLADE 以底盘细胞的基因组尺度代谢模型为起点，搜索目标化合物的候选合成路线，
对路线做通量验证，选择主酶和辅助蛋白，生成 CDS 与表达构建，推荐质粒骨架，最后
输出理论组装后的质粒文件。

完整操作顺序如下：

```text
准备输入配置
    ↓
分析底盘可生成代谢物（chassis）
    ↓
可选：分层扩展底盘可达代谢物（expand）
    ↓
搜索候选合成路线（gap）
    ↓
查看并验证候选路线（info / validate）
    ↓
将路线写入 manifest（write --solution）
    ↓
检索主酶候选并生成主酶组合（main-enzyme / main-enzyme-sets）
    ↓
将主酶组合写入 manifest（write --main-enzyme-set）
    ↓
可选：手动导入辅助蛋白，或运行辅助蛋白研究流程
    ↓
生成或接收 CDS（protein-to-cds）
    ↓
设计并选择表达盒分组（expression --box）
    ↓
推荐并选择表达元件（expression --parts）
    ↓
推荐并选择质粒骨架（plasmid）
    ↓
生成、接受并执行最终组装计划（assembly）
```

`design_manifest.json` 是整个项目的状态中心。路线、蛋白、CDS、表达设计、质粒和组装
信息都会按阶段写入该文件。上游选择发生变化时，系统会清除已经失效的下游区段，
防止旧结果与新输入混用。

## 2. 运行环境

### 2.1 Python 环境

项目要求 Python 3.12。进入项目目录后，可以使用已有虚拟环境：

```powershell
cd F:\myproject\glade
Set-ExecutionPolicy -Scope Process -ExecutionPolicy RemoteSigned
& .\.venv\Scripts\Activate.ps1
```

激活后，本文档中的命令可直接写成：

```powershell
python main.py <命令> <参数>
```

如果使用 `uv` 管理环境，可先同步依赖：

```powershell
uv sync
```

然后将命令中的 `python main.py` 替换为：

```powershell
uv run python main.py
```

### 2.2 网络与服务配置

路径搜索、主酶检索和序列下载会访问 KEGG、Rhea、UniProt 等在线服务。使用对应阶段时，
计算机需要能够访问这些服务。

辅助蛋白研究以及主酶文献检索使用根目录下的 `.env` 模型配置：

```dotenv
MODEL_PROVIDER=openai
AGENT_LLM_MODEL=<模型名称>
API_KEY=<密钥>
BASE_URL=<兼容 OpenAI API 的服务地址>
```

主酶检索还支持以下配置：

```dotenv
GLADE_CONTACT_EMAIL=<联系邮箱>
SELENZYME_REST_URL=<Selenzyme REST 服务地址>
```

其中 `SELENZYME_REST_URL` 在流程需要使用 Selenzyme 回退检索时必须可用。

表达元件和质粒骨架推荐使用远端 Milvus：

```dotenv
MILVUS_HOST=<Milvus 主机>
MILVUS_PORT=19530
MILVUS_TOKEN=<可选令牌>
MILVUS_DB_NAME=<可选数据库名>
```

表达元件集合为 `expression_parts_v3`，质粒集合为 `plasmid_templates_v2`。系统不会在
远端 Milvus 不可用时静默使用过期结果。

## 3. 输入配置

### 3.1 创建项目配置文件

在 `inputs` 目录中创建 JSON 文件，例如 `inputs/demo01.json`：

```json
{
  "target_name": "C00811"
}
```

`target_name` 必须是大写的 KEGG Compound ID，格式为 `C` 加五位数字。

### 3.2 `--input` 的路径规则

所有命令都会自动在 `inputs` 目录下查找 `-i/--input` 指定的文件，因此只传文件名：

```powershell
python main.py chassis -i demo01.json
```

不要传入完整路径，也不要写成 `inputs/demo01.json`，否则程序会再次拼接 `inputs`
目录。

手动导入辅助蛋白时，`--protein-file` 同样只接收 `inputs` 目录下的文件名。

### 3.3 当前默认底盘

当前运行配置使用：

- GEM 模型：`data/gem_models/iML1515.json`；
- 培养基：`data/mediums/default_medium.json`；
- 底盘标识：`ecoli_mg1655`；
- 代谢物分析默认检测胞质区室；
- CDS 优化使用与 *E. coli* MG1655 对应的 CodonTransformer organism ID 52。

同一目标的项目输出统一写入：

```text
outputs/C00811/
```

## 4. 分析底盘可提供的代谢物

执行：

```powershell
python main.py chassis -i demo01.json
```

该命令加载默认 GEM 和培养基，在维持最低生长约束的情况下逐个测试胞内代谢物的
最大 demand flux，并导出具有 KEGG 注释的可生成代谢物。

主要输出：

```text
outputs/C00811/chassis_result/
├── producible_kegg_compounds.csv
└── analyze_chassis_metabolites_summary.csv
```

- `producible_kegg_compounds.csv`：可生成代谢物与 GEM、区室、KEGG ID 的对应表；
- `analyze_chassis_metabolites_summary.csv`：基线生长、最低生长要求、检测数量、可生成
  数量和 KEGG 映射数量等摘要。

查看结果：

```powershell
python main.py info -i demo01.json --chassis
```

## 5. 可选：扩展底盘可达代谢物集合

如果直接使用底盘原始可生成集合找不到路线，可以按 KEGG 反应逐层扩展。扩展深度必须
大于等于 1：

```powershell
python main.py expand -i demo01.json -d 1
```

生成更深层结果时，将 `1` 改为需要的累计深度：

```powershell
python main.py expand -i demo01.json -d 2
```

主要输出位于 `outputs/C00811/chassis_result/`：

```text
chassis_expansion_manifest.json
chassis_frontier_depth_1.csv
chassis_expanded_reachable_depth_1.csv
chassis_frontier_depth_2.csv
chassis_expanded_reachable_depth_2.csv
```

`frontier` 文件只记录该层新增结果；`expanded_reachable` 文件记录截至该深度的累计
可达集合。

查看指定扩展深度：

```powershell
python main.py info -i demo01.json --chassis -d 1
```

## 6. 搜索候选合成路线

### 6.1 使用原始底盘集合

深度 0 表示使用 `chassis` 直接得到的 A0 集合：

```powershell
python main.py gap -i demo01.json -d 0
```

### 6.2 使用扩展集合

使用深度 1 前，必须先完成对应的 `expand -d 1`：

```powershell
python main.py gap -i demo01.json -d 1
```

搜索从目标 KEGG 化合物逆向展开，综合底盘内源反应方向、KEGG module、异源酶数量、
辅因子负担、电子载体风险和循环剪枝生成候选路线。

深度 0 的主要输出为：

```text
outputs/C00811/kegg_gap_C00811/depth0/
├── solutions.csv
├── all_solution_steps.csv
├── rejected_reaction_routes.csv
├── route_electron_requirements.csv
├── solution_electron_summary.csv
└── run_config.json
```

- `solutions.csv`：每条路线的总步骤数、异源步骤数、可达前体和风险摘要；
- `all_solution_steps.csv`：全部路线的逐步反应、方向、底物、产物、EC、KO 和风险字段；
- `rejected_reaction_routes.csv`：因反应规范化门禁等原因不能推荐的路线；
- `route_electron_requirements.csv`：存在电子系统风险的步骤；
- `solution_electron_summary.csv`：每条路线的电子载体平衡与辅助角色需求；
- `run_config.json`：本次搜索参数、版本和实际使用的回退模式。

深度大于 0 时，目录中的 `depth0` 替换为对应深度。

### 6.3 查看路线

查看候选路线摘要：

```powershell
python main.py info -i demo01.json --gap -d 0
```

查看路线 1 的完整正向步骤：

```powershell
python main.py info -i demo01.json --solution 1 -d 0
```

只查看路线 1 的第 2 步：

```powershell
python main.py info -i demo01.json --solution 1 --step 2 -d 0
```

选择路线时应重点关注：异源步骤数、是否可以推荐、电子载体平衡、是否需要额外电子
再生系统、是否需要确认载体兼容性，以及具体反应方向。

## 7. GEM 通量验证

路线写入 manifest 前必须进行独立的 per-solution 验证。

验证路线 1：

```powershell
python main.py validate -i demo01.json -s 1 -m per -c strict -d 0
```

一次验证多条路线：

```powershell
python main.py validate -i demo01.json -s 1 2 3 -m per -c strict -d 0
```

省略 `-s` 时，验证当前深度下的全部候选路线：

```powershell
python main.py validate -i demo01.json -m per -c strict -d 0
```

参数含义：

- `-m per`：每条路线单独加入 GEM 并验证；写入路线前必须存在这种结果；
- `-m pooled`：将所选路线放在同一个模型中联合检查；
- `-m both`：同时生成独立和联合结果；
- `-c strict`：严格处理通用辅因子，默认模式；
- `-c relaxed`：放宽通用辅因子处理，可用于诊断辅因子造成的阻断；
- `-d`：必须与所验证的 `gap` 深度一致。

主要输出：

```text
outputs/C00811/kegg_gap_C00811/depth0/gem_validation/
├── gem_validation_summary.csv
└── gem_validation_route_fluxes.csv
```

只有状态以 `PASS_` 开头的独立验证结果才能写入 manifest。

## 8. 选择并写入路线

路线 1 通过独立验证后执行：

```powershell
python main.py write -i demo01.json --solution 1 -d 0
```

该命令会：

- 核对目标化合物、搜索深度和路线编号；
- 拒绝含阻断反应或不可推荐的路线；
- 核对独立 GEM 验证结果；
- 将路线按实际生物合成方向重新编号；
- 将路线和电子系统信息写入 `design_manifest.json`；
- 清除与旧路线绑定的主酶、辅助蛋白、CDS、表达、质粒和组装选择。

manifest 路径为：

```text
outputs/C00811/design_manifest.json
```

## 9. 主酶选择

### 9.1 生成主酶候选

路线写入 manifest 后执行：

```powershell
python main.py main-enzyme -i demo01.json
```

默认每个异源步骤保留 5 个主酶候选。可以调整数量：

```powershell
python main.py main-enzyme -i demo01.json --top-n 10
```

需要为标准数据库未覆盖的步骤检索论文实验酶活证据时：

```powershell
python main.py main-enzyme -i demo01.json --literature-search
```

主要输出位于：

```text
outputs/C00811/main_protein_selection/
├── main_enzyme_selection.json
├── step_main_enzyme_candidates.csv
├── step_main_enzyme_candidate_audit.csv
├── main_enzyme_candidates.csv
├── reaction_evidence.json
├── direction_evidence.json
├── ko_evidence.json
├── selenzyme_evidence.json
└── route_repair_requests.json
```

`main-enzyme` 只生成候选与证据，不直接修改 manifest。

查看所有步骤的候选：

```powershell
python main.py info -i demo01.json --main-enzyme-candidates
```

只查看第 1 步候选：

```powershell
python main.py info -i demo01.json --main-enzyme-candidates --step 1
```

查看第 1 步排名第 2 的候选详情：

```powershell
python main.py info -i demo01.json --main-enzyme-candidate 2 --step 1
```

### 9.2 生成主酶组合

根据逐步候选生成能够覆盖路线的主酶组合：

```powershell
python main.py main-enzyme-sets -i demo01.json
```

默认最多输出 20 个组合：

```powershell
python main.py main-enzyme-sets -i demo01.json --max-sets 20
```

主要输出：

```text
outputs/C00811/main_protein_selection/
├── main_enzyme_sets.json
├── main_enzyme_sets.csv
└── main_enzyme_set_members.csv
```

查看组合列表：

```powershell
python main.py info -i demo01.json --main-enzyme-sets
```

查看排名第 1 的组合详情：

```powershell
python main.py info -i demo01.json --main-enzyme-set 1
```

### 9.3 写入主酶组合

确认组合后写入 manifest：

```powershell
python main.py write -i demo01.json --main-enzyme-set 1
```

系统会核对组合与当前路线、步骤、反应和候选文件是否一致，并写入
`main_enzyme_selection`。选择主酶组合后，才能导入手动辅助蛋白或运行辅助蛋白研究。

## 10. 辅助蛋白处理

主酶组合写入后，有三种已实现的处理方式：

1. 不添加辅助蛋白，直接运行 `protein-to-cds`；
2. 用户把辅助蛋白氨基酸或 CDS 文件放入 `inputs` 并手动导入；
3. 运行辅助蛋白研究流程，再将研究结果写入 manifest。

### 10.1 手动导入氨基酸序列

将文件放入 `inputs`，例如 `inputs/helper.txt`，内容可以是纯氨基酸文本：

```text
MALWMRLLPLLALLALWGPDPAAA
```

导入：

```powershell
python main.py add-auxiliary-protein -i demo01.json `
  --protein-file helper.txt --sequence-type protein
```

`protein` 表示该序列是氨基酸序列，后续 `protein-to-cds` 会为它生成密码子优化 CDS。

### 10.2 手动导入 CDS

例如将 `helper_cds.fasta` 放入 `inputs`：

```powershell
python main.py add-auxiliary-protein -i demo01.json `
  --protein-file helper_cds.fasta --sequence-type cds
```

`cds` 表示直接使用上传序列，后续跳过密码子优化。导入阶段只删除空白并转为大写，
不检查字符集合、长度是否为 3 的倍数、起始和终止密码子或内部终止密码子。

### 10.3 FASTA 和纯文本规则

- 支持 FASTA、FAA 和纯文本；实际识别依据文件内容；
- 以 `>` 开头时按 FASTA 解析，并支持多条记录；
- FASTA 的 ID 取 header 的第一个字段；
- 非 FASTA 文件作为一条序列，ID 取文件名去除扩展名后的部分；
- ID 会转为大写，并将不适合文件名的字符替换为下划线；
- 多条 FASTA 记录归一化为相同 ID 时，最后一条生效；
- 多次导入会累积；同 ID 再次导入时，以最后一次上传的类型和内容为准；
- 系统在项目输出目录中保存规范化 FASTA 快照，`inputs` 原文件保持不变。

手动导入命令会直接更新 manifest，不需要再运行 `write --auxiliary-protein`。

### 10.4 查看蛋白

查看当前 manifest 中的全部主酶和辅助蛋白：

```powershell
python main.py info -i demo01.json --proteins
```

查看指定蛋白，例如 `HELPER`：

```powershell
python main.py info -i demo01.json --protein HELPER
```

列表会显示蛋白来源、序列类型、负责步骤、CDS 状态、是否可删除以及对应删除命令。
详情默认只显示序列长度和短预览，不输出完整长序列。

### 10.5 删除手动辅助蛋白

删除一个蛋白：

```powershell
python main.py remove-auxiliary-protein -i demo01.json --protein-id HELPER
```

一次删除多个蛋白：

```powershell
python main.py remove-auxiliary-protein -i demo01.json `
  --protein-id HELPER --protein-id CPR
```

删除会移除 manifest 记录和项目中的当前、历史序列快照，并清除旧的 CDS、表达、质粒
和组装选择；不会删除 `inputs` 中的原始文件。删除最后一个手动辅助蛋白后，整个
`auxiliary_protein_selection` 区段会被移除。

### 10.6 辅助蛋白研究流程

运行平衡研究模式：

```powershell
python main.py auxiliary-protein -i demo01.json
```

运行更全面但耗时更长的模式：

```powershell
python main.py auxiliary-protein -i demo01.json --research-mode deep
```

结果写入：

```text
outputs/C00811/protein_selection/auxiliary_protein_research.json
```

研究结果允许继续后，写入 manifest：

```powershell
python main.py write -i demo01.json --auxiliary-protein
```

## 11. 生成或接收 CDS

主酶组合已经写入后执行：

```powershell
python main.py protein-to-cds -i demo01.json
```

明确使用 CPU：

```powershell
python main.py protein-to-cds -i demo01.json --device cpu
```

添加额外禁用 DNA motif，可重复传入：

```powershell
python main.py protein-to-cds -i demo01.json --device cpu `
  --forbidden-motif GAATTC --forbidden-motif GGATCC
```

处理规则：

- 主酶根据 manifest 中的 accession 读取本地缓存或从 UniProt 下载氨基酸序列；
- 手动上传的氨基酸直接从项目快照读取并进行密码子优化；
- 手动上传的 CDS 直接写入选择结果，并标记跳过优化；
- 未添加辅助蛋白时，只处理主酶；
- 优化生成的 CDS 必须保持翻译一致，并通过 GC、CAI、稀有密码子簇、禁用 motif、
  同聚物和起止密码子等门禁；
- 用户直接上传的 CDS 不经过密码子优化门禁。

输出目录：

```text
outputs/C00811/protein_to_cds/
├── uploaded_sequences/manifest_revision_*/<id>.<type>.fasta
├── protein_sequences/<accession>.fasta
├── raw_cds/<accession>.raw.fasta
├── optimized_cds/<accession>.optimized.fasta
├── reports/<accession>.optimization.json
└── run_summary.json
```

运行结果整体写入 manifest 的 `cds_selection`：

- `complete`：全部成功，CLI 退出码为 0；
- `partial`：部分成功，成功产物和失败原因都会保留，CLI 退出码为 2；
- `failed`：全部失败，失败原因写入 manifest，CLI 退出码为 2。

## 12. 表达盒分组

CDS 阶段完成后生成表达盒分组候选：

```powershell
python main.py expression --design --box -i demo01.json
```

结果写入：

```text
outputs/C00811/expression_box/expression_box_designs.json
```

查看终端输出的 `design_id`，例如选择方案 1：

```powershell
python main.py write -i demo01.json --expression-box 1
```

该步骤只确定蛋白如何分组，不选择 promoter、RBS 和 terminator。选择结果写入 manifest
的 `expression_box_selection`。

## 13. 表达元件推荐与选择

表达盒分组写入后，从远端 Milvus 推荐 promoter、RBS 和 terminator：

```powershell
python main.py expression --design --parts -i demo01.json
```

默认请求 12 个方案，也可以请求 3 到 96 个：

```powershell
python main.py expression --design --parts -i demo01.json --n-designs 24
```

每个 RBS 会结合对应 CDS 上下文重新计算翻译起始率，完整表达盒还会接受 GC、同聚物
和禁用酶切位点检查。只有安全且稳定表达评分不低于 70 分的唯一组合会进入结果。

候选文件：

```text
outputs/C00811/expression_box/expression_parts_designs.json
```

写入单个方案：

```powershell
python main.py write -i demo01.json --expression-parts 1
```

写入多个方案或闭区间：

```powershell
python main.py write -i demo01.json --expression-parts 1 3 5
python main.py write -i demo01.json --expression-parts 1:12
python main.py write -i demo01.json --expression-parts 1:4 7 9:12
```

`start:end` 包含两端，重复编号自动去重。写入选择时，系统会同时为每个方案生成完整
串联 GenBank：

```text
outputs/C00811/expression_constructs/
├── design_001.gb
├── design_002.gb
└── ...
```

manifest 的 `parts_selection` 保存选定方案摘要，`assembled_expression_constructs` 保存
完整构建的文件、哈希、坐标和安全审计。

## 14. 质粒骨架推荐与选择

表达构建生成后执行：

```powershell
python main.py plasmid --recommend -i demo01.json
```

默认返回 5 个候选，可使用以下参数：

```powershell
python main.py plasmid --recommend -i demo01.json --n-candidates 10
python main.py plasmid --recommend -i demo01.json --priority balanced
python main.py plasmid --recommend -i demo01.json --preferred-resistance kanamycin
python main.py plasmid --recommend -i demo01.json `
  --exclude-resistance ampicillin tetracycline
```

`--n-candidates` 范围为 1 到 20。`--priority` 支持：

- `stability`：默认，轻微偏向低拷贝稳定性；
- `balanced`：不增加额外拷贝类型偏置；
- `expression`：轻微偏向中、高拷贝，但仍考虑表达负担。

候选结果：

```text
outputs/C00811/plasmid_selection/
├── plasmid_candidates.json
└── candidates/
    ├── candidate_001.gb
    ├── candidate_002.gb
    └── ...
```

例如选择排名第 1 的候选：

```powershell
python main.py write -i demo01.json --plasmid 1
```

选定骨架会复制为：

```text
outputs/C00811/plasmid_selection/selected_backbone.gb
```

选择的是供全部表达构建共同使用的骨架模板。每个表达构建会在最终阶段分别与该骨架
组成一个完整质粒设计。

## 15. 最终组装

### 15.1 生成组装计划

自动为每个表达构建推荐 Gibson 或双酶切方案：

```powershell
python main.py assembly --plan -i demo01.json
```

也可以统一指定方法：

```powershell
python main.py assembly --plan -i demo01.json --method restriction
python main.py assembly --plan -i demo01.json --method gibson
```

结果写入：

```text
outputs/C00811/final_assemble_plan/assembly_plan_recommendations.json
```

指定方法时，如果任一 design 不可行，结果会标记为 `partial` 并禁止写入 manifest，
系统不会自动切换到另一种方法。

### 15.2 接受整套计划

计划完整生成后执行：

```powershell
python main.py write -i demo01.json --assembly-plan
```

整套计划写入 manifest 的 `final_assembly_plan`。计划记录每个表达构建的插入或替换
坐标、骨架线性化方式、限制酶或 Gibson 参数、预计长度、评分、警告和稳定指纹。

### 15.3 执行理论组装

```powershell
python main.py assembly --execute -i demo01.json
```

执行阶段会检查计划指纹、骨架和 insert 文件哈希、限制酶位点、Gibson 同源臂以及
最终输出序列的一致性。

输出：

```text
outputs/C00811/final_assembly/
├── design_001_final.gb
├── design_001_final.fasta
├── design_001_assembly.json
├── ...
├── run_summary.json
└── final_design_report_zh.md
```

结果状态：

- `complete`：全部 design 生成成功；
- `partial`：部分 design 失败，成功文件保留；
- `failed`：没有 design 成功，但仍保存摘要和失败报告。

执行结果写入 manifest 的 `final_assembly` 和 `final_design_report`。

## 16. `info` 查看命令速查

| 目的 | 命令 |
|---|---|
| 查看原始底盘分析 | `python main.py info -i demo01.json --chassis` |
| 查看深度 1 底盘扩展 | `python main.py info -i demo01.json --chassis -d 1` |
| 查看路线摘要 | `python main.py info -i demo01.json --gap -d 0` |
| 查看路线 1 | `python main.py info -i demo01.json --solution 1 -d 0` |
| 查看路线 1 第 2 步 | `python main.py info -i demo01.json --solution 1 --step 2 -d 0` |
| 查看全部主酶候选 | `python main.py info -i demo01.json --main-enzyme-candidates` |
| 查看第 1 步主酶候选 | `python main.py info -i demo01.json --main-enzyme-candidates --step 1` |
| 查看第 1 步候选 2 | `python main.py info -i demo01.json --main-enzyme-candidate 2 --step 1` |
| 查看主酶组合 | `python main.py info -i demo01.json --main-enzyme-sets` |
| 查看组合 1 | `python main.py info -i demo01.json --main-enzyme-set 1` |
| 查看全部蛋白 | `python main.py info -i demo01.json --proteins` |
| 查看蛋白 HELPER | `python main.py info -i demo01.json --protein HELPER` |

## 17. 完整命令示例

下面示例使用原始底盘深度 0、路线 1、主酶组合 1、手动上传一个辅助蛋白，并选择
12 个表达元件方案：

```powershell
# 1. 底盘和路线
python main.py chassis -i demo01.json
python main.py gap -i demo01.json -d 0
python main.py info -i demo01.json --gap -d 0
python main.py info -i demo01.json --solution 1 -d 0

# 2. 验证并选择路线
python main.py validate -i demo01.json -s 1 -m per -c strict -d 0
python main.py write -i demo01.json --solution 1 -d 0

# 3. 主酶
python main.py main-enzyme -i demo01.json
python main.py info -i demo01.json --main-enzyme-candidates
python main.py main-enzyme-sets -i demo01.json
python main.py info -i demo01.json --main-enzyme-sets
python main.py write -i demo01.json --main-enzyme-set 1

# 4. 可选：手动辅助蛋白
python main.py add-auxiliary-protein -i demo01.json `
  --protein-file helper.txt --sequence-type protein
python main.py info -i demo01.json --proteins

# 5. CDS 和表达设计
python main.py protein-to-cds -i demo01.json
python main.py expression --design --box -i demo01.json
python main.py write -i demo01.json --expression-box 1
python main.py expression --design --parts -i demo01.json --n-designs 12
python main.py write -i demo01.json --expression-parts 1:12

# 6. 质粒和最终组装
python main.py plasmid --recommend -i demo01.json
python main.py write -i demo01.json --plasmid 1
python main.py assembly --plan -i demo01.json
python main.py write -i demo01.json --assembly-plan
python main.py assembly --execute -i demo01.json
```

如果不需要辅助蛋白，省略第 4 部分即可。如果使用用户上传的 CDS，将
`--sequence-type protein` 改为 `--sequence-type cds`。

## 18. 结果一致性与重新运行规则

### 18.1 深度必须一致

使用 `gap -d N` 搜索后，`info --gap`、`info --solution`、`validate` 和
`write --solution` 都应使用相同的 `-d N`。

### 18.2 上游变化会使下游结果失效

典型情况包括：

- 重新选择路线后，需要重新运行主酶选择以及全部后续阶段；
- 重新选择主酶组合后，需要重新确认辅助蛋白并重新生成 CDS；
- 新增、替换或删除手动辅助蛋白后，需要重新运行 CDS、表达、质粒和组装阶段；
- CDS 或表达盒变化后，需要重新生成表达元件、质粒和组装结果；
- 表达元件选择变化后，需要重新推荐质粒和生成组装计划；
- 质粒选择变化后，需要重新生成、接受并执行组装计划。

系统使用 manifest revision、输入指纹和文件哈希检测过期结果。出现“输入已经过期”、
“路线已经变化”或“候选文件不一致”等提示时，应从提示所指的上游阶段重新运行，
不要手动复制旧结果规避检查。

### 18.3 重复执行

对于相同且仍然完整的输入，多个写入和生成阶段会复用已有文件；相同选择通常不会
增加 manifest revision。文件缺失或损坏时，支持自愈的阶段会重新生成文件。

## 19. 常见问题

### 19.1 提示找不到输入文件

确认文件位于 `inputs`，并且 `-i` 后只写文件名：

```powershell
python main.py chassis -i demo01.json
```

### 19.2 提示 target 格式错误

检查 `target_name` 是否为大写 `Cxxxxx`，例如 `C00811`。

### 19.3 `gap` 提示缺少 chassis 结果

先执行：

```powershell
python main.py chassis -i demo01.json
```

如果使用 `gap -d 1`，还必须先执行 `expand -d 1`。

### 19.4 路线无法写入 manifest

确认：

- 路线在 `solutions.csv` 中标记为可以推荐；
- 路线没有 blocking reaction；
- 已运行 `validate -m per`；
- 独立验证状态以 `PASS_` 开头；
- `write --solution` 使用了相同深度。

### 19.5 主酶候选为空或网络检索失败

检查 KEGG、Rhea、UniProt 和 Selenzyme 服务是否可访问；如果启用了
`--literature-search`，同时检查 `.env` 中的模型配置。

### 19.6 辅助蛋白研究失败

检查 `.env` 中的 `MODEL_PROVIDER`、`AGENT_LLM_MODEL`、`API_KEY` 和 `BASE_URL`，
以及外部数据库网络连接。也可以使用手动上传方式添加已知辅助序列。

### 19.7 表达元件或质粒推荐失败

检查 `MILVUS_HOST`、`MILVUS_PORT`、可选认证信息，以及远端 collection 的名称和
schema 是否与当前程序匹配。

### 19.8 `protein-to-cds` 返回退出码 2

查看 `outputs/C00811/protein_to_cds/run_summary.json` 和 manifest 的
`cds_selection`。退出码 2 表示 `partial` 或 `failed`，成功项目和失败原因都会保留。

## 20. 结果解释

GLADE 输出的是计算设计与候选工程方案。路线搜索、GEM 通量可行性、主酶评分、表达
评分、质粒适配评分和理论组装结果都不等同于湿实验验证，也不能单独证明目标产物能够
在实际培养条件下稳定合成。进入实验前，应复核反应方向、辅因子与电子伙伴、酶活性、
蛋白表达、细胞毒性、质粒稳定性和组装策略，并根据实验结果迭代设计。
