# 系统质粒选择

这一阶段从当前 `design_manifest.json` 读取已经选定并串联完成的表达构建，随后从远端 Milvus 的 `plasmid_templates_v2` 集合中推荐一个可用于全部构建的统一质粒骨架。

它不会把 12 个表达构建装进同一个质粒。这里选择的是一个共同骨架；后续最终组装阶段会分别生成“骨架 + 表达构建 1”到“骨架 + 表达构建 12”，因此最终对应 12 个不同的完整质粒设计。

## 1. 生成质粒候选

```powershell
uv run python main.py plasmid --recommend -i demo04.json
```

默认返回 5 个有差异的候选，终端只显示排名、名称、分数、拷贝类型、抗性标记、组装策略和表达负担摘要。完整证据和逐构建评分写入：

```text
outputs/<目标化合物>/plasmid_selection/plasmid_candidates.json
```

经过下载和序列校验的候选 GenBank 文件写入：

```text
outputs/<目标化合物>/plasmid_selection/candidates/candidate_001.gb
outputs/<目标化合物>/plasmid_selection/candidates/candidate_002.gb
...
```

可选参数：

```powershell
uv run python main.py plasmid --recommend -i demo04.json --n-candidates 5
uv run python main.py plasmid --recommend -i demo04.json --priority balanced
uv run python main.py plasmid --recommend -i demo04.json --preferred-resistance kanamycin
uv run python main.py plasmid --recommend -i demo04.json --exclude-resistance ampicillin tetracycline
```

`--priority` 支持，但只对动态拷贝适配施加最多 3 分的偏置：

- `stability`：默认，轻微偏向低拷贝，但不会压过明显更适合的中拷贝。
- `balanced`：不增加拷贝偏置，完全使用负担适配基础分。
- `expression`：轻微偏向中/高拷贝，但高负担设计仍不会被推到高拷贝首选。

系统默认的抗性标记偏好顺序为：卡那霉素、氯霉素、庆大霉素、链霉素/壮观霉素、四环素、氨苄青霉素。`--preferred-resistance` 会改变评分偏好，`--exclude-resistance` 是硬排除条件。

## 2. 写入选定骨架

例如采用排名第一的候选：

```powershell
uv run python main.py write -i demo04.json --plasmid 1
```

该命令会同时：

- 校验候选文件是否仍与当前表达构建、候选快照和序列哈希一致；
- 把候选完整信息写入 manifest 的 `plasmid_selection`；
- 把原始候选 GenBank 精确复制为 `plasmid_selection/selected_backbone.gb`；
- 登记该骨架覆盖的全部表达元件方案编号；
- 清理已经失效的 `final_assembly_plan` 和 `final_assembly`。

重复选择相同候选时不会增加 manifest revision。若 `selected_backbone.gb` 丢失或损坏，系统只修复文件，也不会增加 revision。改选另一个候选才会更新 manifest 并增加 revision。

## 推荐规则

推荐不使用 BGE-M3 或语义向量检索，而是直接查询 Milvus 的结构化字段。候选首先必须通过以下硬条件：

- `plasmid_template.v2` 且序列审计为 PASS；
- 兼容 *E. coli* K-12 MG1655，拓扑为 circular；
- 复制起点、抗性标记、宿主兼容性、插入区域和保护元件信息完整；
- 不依赖当前流程不支持的复制因子，不携带表达 cargo；
- 组装策略为 `insert_into_mcs` 或 `replace_seva_cargo_paci_spei`；
- 原始 GenBank 可以下载，并通过长度、碱基字符、拓扑及序列 SHA-256 校验。

每个骨架会分别与所有已选表达构建配对评分。总分为 100 分：

- 拷贝数与宿主负担适配：35 分；
- 插入区域和组装就绪程度：25 分；
- 来源、证据和序列完整性：20 分；
- 抗性标记适用性：10 分；
- 估算最终质粒长度：10 分。

其中 35 分的拷贝适配不再固定认为低拷贝最好，而是读取 manifest 中每个方案的
`expression_burden.v1`：

| 表达负担 | 低拷贝 | 中拷贝 | 条件依赖 | 高拷贝 |
|---|---:|---:|---:|---:|
| low | 24 | 35 | 22 | 30 |
| moderate | 30 | 35 | 20 | 18 |
| high | 35 | 25 | 12 | 5 |

因此，中等负担通常优先中拷贝，高负担优先低拷贝；低负担可以保留中、高拷贝候选。
插入片段达到 10 kb 时该项再扣 5 分，达到 15 kb 时扣 10 分。

候选的最终稳健分是它对全部表达构建的配对分数中的最低值。因此，一个只适合部分方案的骨架不会得到虚高排名。第一名始终是原始最高分；其余候选会优先保持复制子、抗性标记、组装策略和拷贝类型的多样性。

这些分数是透明、可复查的工程启发式分数，不是实验成功概率。当前数据库没有经过实验验证的最大插入容量字段，所以容量依据拷贝类型和预计最终长度估算，候选文件会保留这一警告。

旧版表达结果没有数值负担，升级后请按以下顺序重新生成：

```powershell
uv run python main.py expression --design --parts -i demo04.json
uv run python main.py write -i demo04.json --expression-parts 1:12
uv run python main.py plasmid --recommend -i demo04.json
```

## 缓存与远端状态

每次运行推荐命令都会重新查询远端 Milvus，以确认 collection schema 和候选内容快照。只有当前 manifest 输入、远端快照、参数以及全部本地候选文件都没有变化时，才复用已有结果。候选文件缺失或损坏时会重新下载并重建；Milvus 不可用时不会静默使用过期离线结果。
