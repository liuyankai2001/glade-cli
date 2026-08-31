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
搜索候选合成路线（默认使用 KEGG；需要预测反应时显式启用 RetroPath）
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

### 2.3 使用 RetroPath 前启动本地服务

只有执行带 `--retropath` 的路线搜索时才需要本地 RetroPath 服务。先确认 Docker
Desktop 已启动，并检查规则文件已经放到项目指定位置：

```powershell
docker version
Get-Item data\retropath\rules\rr02\retrorules_rr02_rp2_flat_retro.csv
```

首次使用时构建并启动服务。首次构建需要下载较大的运行环境，可能耗时较长；以后镜像
没有变化时只执行第二条启动命令即可：

```powershell
docker compose -f compose.retropath.yml build retropath
docker compose -f compose.retropath.yml up -d retropath
```

确认服务已经就绪：

```powershell
Invoke-RestMethod http://127.0.0.1:8765/health
```

返回结果中的 `ready` 必须为 `true`。GLADE 会自动准备输入并调用该服务，用户不需要
手动提交 CSV 或调用 HTTP 接口。使用结束后可以停止服务：

```powershell
docker compose -f compose.retropath.yml down
```

不要使用 `down -v`，否则会同时删除服务保存的任务和结果。

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

如果直接使用底盘原始可生成集合找不到路线，可以按 KEGG 反应逐层扩展。每个定向
KEGG 反应计一层；同一酶连续催化两个反应仍计两层。扩展深度必须大于等于 1：

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

扩展采用载体感知策略：普通主底物必须全部位于上一层累计集合；P450 还原酶、
ferredoxin、thioredoxin 等已识别电子载体不阻断主产物，但载体自身不会作为新增的
可达代谢物，也不会作为 RetroPath 的路线连接点。CSV 和 manifest 会记录电子载体、
风险、净变化、辅助角色以及
`auxiliary_requirements_json`。因此“扩展可达”表示补充相应 KEGG 反应和工程辅助系统后
可以抵达，不表示底盘天然已经具备这些酶和电子再生能力。

注释为 `first/second/... step of ... reaction` 的 KEGG 条目是多步路线中已经拆分好的
独立组件反应，每个仍计一层；只有 `three-step reaction (see R...+R...)` 这类汇总条目
继续被拒绝，避免把多个酶促步骤冒充一层。

扩展策略升级后，旧版 `chassis_forward_expansion.v1` 和
`chassis_forward_expansion.v2_carrier_aware` 结果不可复用；再次使用对应 depth 前需
重新执行 `expand -d N`。

查看指定扩展深度：

```powershell
python main.py info -i demo01.json --chassis -d 1
```

## 6. 搜索候选合成路线

### 6.1 使用原始底盘集合

深度 0 表示直接使用 `chassis` 得到的底盘可生成代谢物集合：

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

### 6.4 使用 RetroPath 搜索预测路线

当 KEGG 搜索没有合适路线，或者希望尝试基于化学结构的预测反应时，可以显式启用
RetroPath。运行前先按 2.3 节确认本地服务已经就绪：

```powershell
python main.py gap -i demo01.json --retropath -d 0
```

`-d 0` 表示预测路线必须连接到底盘直接可生成的代谢物。若想让路线连接到扩展后的
可达代谢物，必须先生成同一深度的扩展结果。例如使用深度 5：

```powershell
python main.py expand -i demo01.json -d 5
python main.py gap -i demo01.json --retropath -d 5
```

默认 `gap` 只运行 KEGG 搜索，不会自动调用 RetroPath；带 `--retropath` 的命令也只
运行 RetroPath，不会自动先补跑 KEGG 搜索。RetroPath 从目标化合物逆向预测，只保留
能够完整连接到底盘可达代谢物的路线。若某个反应需要多个前体，则所有必要前体都必须
能够连接到底盘，缺少任一分支都不会作为可用路线。

搜索完成后，建议按下面的顺序查看结果：

```powershell
# 先看本次运行是否成功以及找到了多少条路线
python main.py info -i demo01.json --retropath -d 0

# 再看排名第 1 的 RetroPath 候选详情
python main.py info -i demo01.json --retropath-candidate 1 -d 0

# 候选详情会给出对应的路线编号，假设为 N
python main.py info -i demo01.json --solution N -d 0
```

这里有两种编号，不能混用：

- `--retropath-candidate 1` 中的 `1` 是 RetroPath 候选排名，只用于查看预测详情；
- `--solution N` 中的 `N` 是 GLADE 路线编号，用于查看、验证和写入路线。

候选详情中显示的 `正式Solution编号` 就是这里所说的 GLADE 路线编号。

每条可用的 RetroPath 候选都会获得一个 GLADE 路线编号，并加入当前深度的路线列表。
已有 KEGG 路线的编号保持不变，RetroPath 路线从当前最大编号之后继续编号。

最常用的结果文件位于：

```text
outputs/C00811/kegg_gap_C00811/depthN/retropath/
├── pipeline_result.json
├── candidate_routes.csv
├── candidate_steps.csv
└── rejected_routes.csv
```

- `pipeline_result.json`：本次服务运行状态、完整连接数量、候选数量和失败原因；
- `candidate_routes.csv`：可以完整连接到底盘的候选路线摘要；
- `candidate_steps.csv`：每条候选路线的逐步反应；
- `rejected_routes.csv`：路线未连接完整、结构冲突或证据不足等拒绝原因。

“已经生成预测反应网络”不等于“已经找到完整路线”。如果预测网络没有连接到底盘
可达代谢物，候选数量仍为 0，也不会把不完整路线加入 `solutions.csv`。此时执行：

```powershell
python main.py info -i demo01.json --retropath -d 5
```

重点查看服务状态、完整连接数量、候选数量和拒绝原因。如果显示已经产生预测网络但
候选为 0，应理解为“本次预测尚未得到一条能完整连接到底盘的路线”。

## 7. GEM（底盘代谢模型）通量验证

GEM 验证用于判断候选路线在当前底盘模型和培养基中是否具备通量可行性。它对 KEGG
和 RetroPath 路线都是可选的：未验证或验证失败的路线仍可写入设计清单，但系统会保留
验证状态和人工复核提示。

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

- `-m per`：每条路线单独加入 GEM 并验证，适合为单条路线生成独立证据；
- `-m pooled`：将所选 KEGG 路线放在同一个模型中联合检查；
- `-m both`：对所选 KEGG 路线同时生成独立和联合结果；
- `-c strict`：严格处理通用辅因子，默认模式；
- `-c relaxed`：放宽通用辅因子处理，可用于诊断辅因子造成的阻断；
- `-d`：必须与所验证的 `gap` 深度一致。

只要待验证列表中包含 RetroPath 路线，就必须使用 `-m per`。省略 `-s` 时，程序会自动
识别当前深度下每条路线来自 KEGG 还是 RetroPath，并分别验证。

主要输出：

```text
outputs/C00811/kegg_gap_C00811/depth0/gem_validation/
├── gem_validation_summary.csv
└── gem_validation_route_fluxes.csv
```

验证结果会记录为“未运行”“通过”或“失败”，并保存通量、辅因子模式和问题说明。
验证失败只增加人工复核提示，不会单独阻止路线写入。

### 7.1 可选验证 RetroPath 路线

RetroPath 路线在搜索完成后就可以直接查看或写入，不需要先验证。假设候选详情显示其
对应的路线编号为 `N`：

```powershell
python main.py info -i demo01.json --solution N -d 0
python main.py write -i demo01.json --solution N -d 0
```

没有运行 GEM 验证时，系统会明确标记“尚未验证”和“需要人工复核”，但允许继续主酶、
CDS 和表达设计流程。

如果需要验证 RetroPath 路线，首次使用前先安装计量补全所需的 MNXref v3.0 数据：

```powershell
python -m src.pathway_analyze.retropath_mnxref install
```

安装器会下载并校验官方源文件，只保留当前预测规则需要的反应、化合物和映射。检查
安装状态：

```powershell
python -m src.pathway_analyze.retropath_mnxref status
```

验证当前深度的全部路线：

```powershell
python main.py validate -i demo01.json -m per -d 0
```

使用宽松辅因子模式诊断全部路线：

```powershell
python main.py validate -i demo01.json -m per -c relaxed -d 0
```

只验证路线 4 和 5：

```powershell
python main.py validate -i demo01.json -s 4 5 -m per -d 0
```

RetroPath 路线只支持独立验证，不能使用 `-m pooled` 或 `-m both`。严格和宽松模式都
只使用预测规则能够追溯到的数据库反应补全化学计量，不会根据 EC 编号或元素差额自行
猜测辅因子。宽松模式仅放开路线实际涉及的通用载体，用于判断路线是否主要受辅因子
约束；结果中会明确记录放开了哪些载体。

主要输出：

```text
outputs/C00811/kegg_gap_C00811/depth0/retropath/gem_validation/
├── stoichiometry_hypotheses.csv
├── stoichiometry_terms.csv
├── rejected_hypotheses.csv
├── gem_validation_summary.csv
├── gem_validation_route_fluxes.csv
└── validation_manifest.json
```

结果可以这样理解：

- 通过：至少存在一套有数据库来源的完整反应计量，能够同时满足底盘生长、目标产出和
  路线中每一步都有通量；
- 失败：当前证据和验证模式下没有找到可行计量与通量，路线仍可写入，但必须人工复核；
- 未运行：没有验证证据，路线仍可写入，并保留“尚未验证”警告。

验证只更新路线的计量和 GEM 证据，不会改变路线编号或步骤。一次只验证部分路线时，
只有选中的路线获得本次验证结果。系统会校验候选文件和验证文件的一致性；如果文件被
手动修改或上游搜索结果已经变化，应重新运行同一深度的 `gap --retropath`，不要手动
拼接或覆盖 CSV。

### 7.2 为 RetroPath 路线生成主酶候选

先用路线编号 `N` 将选中的 RetroPath 路线写入设计清单。GEM 验证可以先运行，也可以
跳过：

```powershell
python main.py write -i demo01.json --solution N -d 0
```

随后与纯 KEGG 路线使用完全相同的命令：

```powershell
python main.py main-enzyme -i demo01.json
python main.py main-enzyme-sets -i demo01.json
python main.py info -i demo01.json --main-enzyme-sets
python main.py write -i demo01.json --main-enzyme-set 1
```

`main-enzyme` 会从设计清单自动识别普通 KEGG 步骤和 RetroPath 预测步骤，不需要再传
候选排名或搜索深度。普通步骤继续使用 KEGG、Rhea、KO、文献和 Selenzyme 证据；预测
步骤会综合来源酶注释、反应映射和结构相似性检索。未运行 GEM 验证时也可以检索主酶；
验证通过后，系统会优先使用补全后的完整反应和精确数据库映射。

主要输出：

```text
outputs/C00811/main_protein_selection/
├── main_enzyme_selection.json
├── step_main_enzyme_candidates.csv
├── step_main_enzyme_candidate_audit.csv
├── retropath_enzyme_requirements.json
└── retropath_selenzyme_evidence.json
```

检索优先使用精确的 KEGG/Rhea 反应映射，其次使用预测规则附带的 EC、Rhea 和 UniProt
来源信息，最后才使用反应结构查询 SelenzymeRF。结构相似结果始终只是预测性证据，
即使相似度为 1 也需要人工复核。这类候选可以进入主酶组合并继续后续设计，但系统不会
把它标记成已经通过实验验证的酶活。

## 8. 选择并写入路线

选择路线 1 后即可执行；GEM 验证可在写入前按需运行：

```powershell
python main.py write -i demo01.json --solution 1 -d 0
```

该命令会：

- 核对目标化合物、搜索深度和路线编号；
- 拒绝含阻断反应或不可推荐的路线；
- 读取可选的独立 GEM 验证结果；没有结果时记录为 `not_run`；
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
├── taxonomy_evidence.json
├── selenzyme_evidence.json
└── route_repair_requests.json
```

`main-enzyme` 只生成候选与证据，不直接修改 manifest。

候选蛋白综合评分使用以下固定权重：反应功能 40%、UniProt/实验依据 25%、
表达风险 20%、来源分类学适配 15%。来源适配不再依赖写死的物种名称顺序：系统读取
底盘和候选蛋白的 UniProt taxon lineage，按最近共同祖先（LCA）所在的 strain、
species、genus、family、order、class、phylum、kingdom 或 domain 层级评分。
`taxonomy_evidence.json` 保存底盘 taxon、完整 ranked lineage、评分表、权重及数据来源；
逐步候选 CSV 还会记录每个候选的共同祖先、匹配层级和分类来源分。

分类学信息缺失时使用中性分 50，并明确标记为 `unknown`；数据缺失不会被误判为
远缘，也不会单独阻断主酶选择。分类亲缘性只是排序因素，反应、方向、底物/产物
特异性以及辅助蛋白风险仍优先。当前代谢底盘仍固定为 *E. coli* MG1655/iML1515，
本次改动只将分类学评分内核通用化，并未开放与 GEM 不一致的任意底盘参数。

本版输出格式为 `main_enzyme_selection.v3` 和 `main_enzyme_sets.v3`。已有 v2 结果不会
静默迁移；升级后需要依次重新运行 `main-enzyme`、`main-enzyme-sets`，并重新执行
`write --main-enzyme-set N`。

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
| 查看 RetroPath 运行和候选摘要 | `python main.py info -i demo01.json --retropath -d 0` |
| 查看 RetroPath 候选 1 | `python main.py info -i demo01.json --retropath-candidate 1 -d 0` |
| 查看 RetroPath 候选 1 第 2 步 | `python main.py info -i demo01.json --retropath-candidate 1 --step 2 -d 0` |
| 查看全部主酶候选 | `python main.py info -i demo01.json --main-enzyme-candidates` |
| 查看第 1 步主酶候选 | `python main.py info -i demo01.json --main-enzyme-candidates --step 1` |
| 查看第 1 步候选 2 | `python main.py info -i demo01.json --main-enzyme-candidate 2 --step 1` |
| 查看主酶组合 | `python main.py info -i demo01.json --main-enzyme-sets` |
| 查看组合 1 | `python main.py info -i demo01.json --main-enzyme-set 1` |
| 查看全部蛋白 | `python main.py info -i demo01.json --proteins` |
| 查看蛋白 HELPER | `python main.py info -i demo01.json --protein HELPER` |

RetroPath 搜索失败时，也可以用 `info --retropath` 查看失败位置和原因。只有成功找到
候选后，才能使用 `info --retropath-candidate N` 查看排名第 `N` 的预测详情。该编号
只表示候选排名；验证和写入时要使用候选详情中显示的 GLADE 路线编号。

系统会检查目标、搜索深度和结果文件是否匹配。出现“候选文件校验失败”时，应重新运行
同一深度的 `gap --retropath`；需要验证时再使用对应路线编号运行 `validate -s N`，不要
手动修改结果 CSV。

## 17. 完整命令示例

下面示例使用原始底盘深度 0、路线 1、主酶组合 1、手动上传一个辅助蛋白，并选择
12 个表达元件方案：

```powershell
# 1. 底盘和路线
python main.py chassis -i demo01.json
python main.py gap -i demo01.json -d 0
python main.py info -i demo01.json --gap -d 0
python main.py info -i demo01.json --solution 1 -d 0

# 2. 可选验证并选择路线
# 如果暂时不需要 GEM 验证，可以省略下一行
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

KEGG 无解后改用 RetroPath 时，路线与主酶阶段为：

```powershell
python main.py gap -i demo01.json --retropath -d 0
python main.py info -i demo01.json --retropath -d 0

# 查看候选排名 1，并记下其中显示的 GLADE 路线编号 N
python main.py info -i demo01.json --retropath-candidate 1 -d 0
python main.py info -i demo01.json --solution N -d 0

# 可选：验证路线 N；不验证也可以直接写入
python main.py validate -i demo01.json -s N -m per -c strict -d 0
python main.py write -i demo01.json --solution N -d 0

python main.py main-enzyme -i demo01.json
python main.py main-enzyme-sets -i demo01.json
python main.py info -i demo01.json --main-enzyme-sets
python main.py write -i demo01.json --main-enzyme-set 1
```

## 18. 结果一致性与重新运行规则

### 18.1 深度必须一致

假设搜索时使用 `-d 3`，后续查看、验证和写入这批路线时也都要使用 `-d 3`。这里的
`-d` 是搜索深度，`-s` 是路线编号，两者不是同一个值。例如验证路线 7：

```powershell
python main.py validate -i demo01.json -s 7 -m per -d 3
```

`info --gap`、`info --retropath`、`info --retropath-candidate`、`info --solution` 和
`write --solution` 同样要使用搜索时的深度。路线写入设计清单后，`main-enzyme` 不再需要
深度参数。

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

### 19.4 RetroPath 提示本地服务不可用

先确认 Docker Desktop 正在运行，然后启动服务：

```powershell
docker compose -f compose.retropath.yml up -d retropath
Invoke-RestMethod http://127.0.0.1:8765/health
```

如果容器无法启动，检查规则文件是否存在：

```powershell
Get-Item data\retropath\rules\rr02\retrorules_rr02_rp2_flat_retro.csv
docker compose -f compose.retropath.yml logs --tail 100 retropath
```

健康检查中的 `ready` 为 `true` 后，再重新运行 `gap --retropath`。

### 19.5 验证 RetroPath 路线时提示缺少 MNXref

安装并检查计量补全数据：

```powershell
python -m src.pathway_analyze.retropath_mnxref install
python -m src.pathway_analyze.retropath_mnxref status
```

### 19.6 路线无法写入设计清单

确认：

- 路线在 `solutions.csv` 中标记为可以推荐；
- 路线没有阻断反应；
- `write --solution` 使用了相同深度。

GEM 验证不是写入前置条件。未验证或验证失败只会保留状态和人工复核警告。如果已经
运行验证，还应确认验证结果属于当前路线和当前搜索深度，并且没有手动修改上游文件。

### 19.7 主酶候选为空或网络检索失败

检查 KEGG、Rhea、UniProt 和 Selenzyme 服务是否可访问；如果启用了
`--literature-search`，同时检查 `.env` 中的模型配置。

RetroPath 结构检索不可用时，已经获得的来源模板候选仍会保留，但缺少候选的预测步骤
不会假装完成。检查 `.env` 中的 `SELENZYME_REST_URL` 后重新运行统一的
`main-enzyme -i demo01.json`。`--literature-search` 仍用于现有 KEGG 步骤；预测步骤继续
使用可审计的来源模板和结构证据。

### 19.8 辅助蛋白研究失败

检查 `.env` 中的 `MODEL_PROVIDER`、`AGENT_LLM_MODEL`、`API_KEY` 和 `BASE_URL`，
以及外部数据库网络连接。也可以使用手动上传方式添加已知辅助序列。

### 19.9 表达元件或质粒推荐失败

检查 `MILVUS_HOST`、`MILVUS_PORT`、可选认证信息，以及远端 collection 的名称和
schema 是否与当前程序匹配。

### 19.10 `protein-to-cds` 返回退出码 2

查看 `outputs/C00811/protein_to_cds/run_summary.json` 和 manifest 的
`cds_selection`。退出码 2 表示 `partial` 或 `failed`，成功项目和失败原因都会保留。

## 20. 结果解释

GLADE 输出的是计算设计与候选工程方案。路线搜索、GEM 通量可行性、主酶评分、表达
评分、质粒适配评分和理论组装结果都不等同于湿实验验证，也不能单独证明目标产物能够
在实际培养条件下稳定合成。进入实验前，应复核反应方向、辅因子与电子伙伴、酶活性、
蛋白表达、细胞毒性、质粒稳定性和组装策略，并根据实验结果迭代设计。
