# RetroPath 接入修改计划

> 文档状态：P0–P11.1 已完成；2026-08-26 已将 P10 调整为“先物化全部 Top-K，P8 作为可选验证覆盖层”
> 制定日期：2026-08-24  
> 适用项目：GLADE  
> 目标：为现有 KEGG 通路搜索增加由用户显式启用的 RetroPath 预测搜索；预测路线可直接进入统一 solution 流程，但必须持续携带预测、验证状态和人工复核信息。

## 1. 已确认的产品决策

| 决策项 | 确认结果 |
|---|---|
| 默认搜索 | 保持现有 KEGG 搜索，不自动调用 RetroPath |
| 启用方式 | 用户通过 <code>--retropath</code> 显式启用 |
| 自动兜底 | 不做；KEGG 无解时只提示用户，不自动执行 expand 或 RetroPath |
| depth 范围 | <code>depth >= 0</code> |
| depth 0 | sink 使用底盘 GEM/FBA 直接可生成集合 A0 |
| depth N | sink 使用累计集合 AN = A0 ∪ F1 ∪ … ∪ FN |
| expand 行为 | 系统不自动执行；depth > 0 时要求用户事先生成相应扩展结果 |
| 搜索方向 | RetroPath 从目标结构逆向搜索到 sink；拼接前翻转为合成方向 |
| 完整路线 | 所有边界 sink 的 A0→X KEGG witness 与 X→目标 RetroPath 预测 DAG 拼接 |
| 预测反应 ID | 使用 RP2:哈希，禁止伪造 Rxxxxx |
| 预测中间体 ID | 使用 RP2CPD:InChIKey或结构哈希，禁止伪造 Cxxxxx |
| Top-K 物化 | P5 产生的全部 Top-K 混合路线在 `gap --retropath` 结束时立即写入统一 solution |
| P8 验证 | 可选；只更新原 solution 的计量/GEM 状态和证据，不新增、删除或重新编号 solution |
| 未验证/验证失败 | 都允许 `write --solution N` 和后续主酶流程，但强制标记预测风险与人工复核 |
| 路线身份 | solution 使用稳定 `RP2STEP` 预测反应身份；P8 的 `RP2STOICH` 作为计量假设证据挂载，不替换 solution ID |

### P0 本地服务架构（2026-08-24 更新）

RetroPath2.0 不再作为 GLADE 进程内依赖或直接外部命令接入。P0 先提供只监听
`127.0.0.1:8765` 的本地 Docker HTTP 服务，容器内固定
`retropath2_wrapper 3.9.1`、KNIME `4.7.0` 和 workflow `r20260212`。
GLADE 仍由 uv 管理；后续 P3 将从“本地程序 runner”调整为“HTTP client”。

P0 使用现有 `rr02-rp2-hs` 逆向规则文件，只读挂载，不写入镜像或 Git。
服务首先返回 KNIME 原始结果和运行清单，候选路线解析仍由后续 P4 完成。

P0 验收记录：

- 镜像：`glade-retropath:3.9.1-r20260212`，本地 image ID
  `sha256:e9cdfcefd2f39e9eacb8dfab9d4cc61e3eaf29d32b04164369d2f43702741b82`；
- 容器内单元测试：17 项通过；
- localhost 接口与真实 KNIME/RR02 冒烟：2 项通过；
- 真实任务 `rp2-53cec4814e19468d87a8dd4194ecbe69` 完成 KNIME 执行并返回
  `no_solution`（Wrapper 退出码 11），证明运行环境可用且无环境依赖错误；
- `/health` 已核对 Wrapper、workflow、KNIME、RDKit 和 RR02 SHA-256。

### P1 预测数据模型（2026-08-25 更新）

P1 已在 `src/pathway_analyze/retropath_models.py` 建立版本化领域模型，包含
`PredictedCompound`、`PredictedReaction`、`CandidateRoute`、
`RetroPathRuntimeProvenance` 和 `RetroPathRunResult`。模型使用冻结 dataclass，
只依赖 Python 标准库，不引入 RetroPath、RDKit 或网络运行依赖。

P1 验收记录：

- schema version 固定为 `1`，运行结果支持显式 dict/JSON 往返；
- 预测反应使用 `RP2:<64位SHA-256>`，预测中间体使用
  `RP2CPD:<InChIKey或64位SHA-256>`，候选路线使用
  `RP2ROUTE:<64位SHA-256>`；
- 哈希输入使用键排序的规范 JSON；集合型证据排序去重，反应两侧保留计量重复，
  路线步骤保留顺序；
- 反应身份显式包含 `retrosynthetic`/`biosynthetic` 方向，反转方向或反应两侧会
  生成不同 ID；
- 成功、无解和 source-in-sink 结果必须保留 Wrapper、workflow、KNIME、RDKit
  插件以及 RetroRules 版本和校验和；
- `tests/test_retropath_models.py` 的 17 项离线测试通过，其中包含 v1 ID 黄金值，
  防止哈希契约被静默修改。

### P2 结构与输入生成（2026-08-25 更新）

P2 已在 `src/pathway_analyze/retropath_structure.py` 和
`src/pathway_analyze/retropath_input.py` 实现 KEGG MOL 获取、校验缓存、RDKit
结构标准化、累计 sink 构建、结构去重以及拒绝审计。

P2 验收记录：

- uv 直接依赖固定为 `rdkit==2026.3.5`，Python 3.12 环境已验证；
- KEGG MOL 缓存保存在 `cache/kegg/mol/`，元数据记录 URL、获取时间和 MOL
  SHA-256，损坏或元数据不一致的缓存会重新获取；
- 结构输出包含标准 InChI、完整 InChIKey、canonical isomeric SMILES、分子式、
  总形式电荷、MOL SHA-256 和 RDKit 版本 provenance；
- depth 0 使用 A0，depth N 直接使用 `ExpansionBundle.reachable_compounds` 与
  `depth_by_compound`，不会自动重新执行 expand；
- sink 先按 Cxxxxx 取最小 depth，再按完整 InChIKey 去重；同结构代表 ID 按
  “最小 depth、再按字典序”选择，全部别名保留在 mapping；
- target 结构失败或有效 sink 为空会阻断；单个 sink 失败会排除并写入拒绝表；
- 9 项结构测试和 12 项输入测试通过；联合 P1、expansion 与电子平衡回归共
  55 项通过；真实 KEGG `C00001`、`C05432` 冒烟通过。

### P3 本地 HTTP client（2026-08-25 更新）

P3 已在 `src/pathway_analyze/retropath_client.py` 实现 GLADE 到 P0 本地 Docker
服务的同步 HTTP client。客户端默认只允许 loopback HTTP 地址，使用
`httpx==0.28.1` 且禁用环境代理和重定向。

P3 验收记录：

- 参数范围与 P0 服务一致，默认请求超时 30 秒、GET 三次重试、0.5 秒指数退避、
  1 秒轮询和 3900 秒客户端等待上限；
- POST 不自动重试；提交响应丢失会写入 `submission_uncertain` 状态并要求人工确认
  或显式 `force=True`，避免生成重复 KNIME 任务；
- `queued`/`running` 可在进程重启或客户端等待超时后用原 job ID 继续轮询；客户端
  等待超时不会伪装成服务端 `timed_out`；
- 服务 manifest 的 job、参数、输入 SHA-256、Wrapper、workflow、KNIME、RDKit
  节点、RetroRules 版本及规则 SHA-256 均与 P2/health 交叉校验；
- artifact 使用服务白名单、路径越界检查、流式临时文件、原子替换和本地
  SHA-256 复核；
- 仅当服务健康信息、输入、参数、运行栈和全部本地 artifact 校验和一致时复用
  `succeeded`、`no_solution` 或 `source_in_sink`；失败和服务超时不缓存；
- 本地 artifact 损坏时优先从已完成的同一服务 job 重新下载，不重复运行 KNIME；
- 14 项 MockTransport 离线协议测试通过；联合 P1–P3、expansion 与电子平衡回归
  共 69 项。P3 实施时 Docker Desktop 未运行，因此没有重复执行 P0 真实服务冒烟。

### P4 网络解析与路径枚举（2026-08-25 更新）

P4 已在 `src/pathway_analyze/retropath_parser.py` 和
`src/pathway_analyze/retropath_routes.py` 实现 RetroPath scope/results 解析、RR02
规则证据交叉校验和完整 sink 路径枚举。

P4 验收记录：

- 优先使用 `target_scope.json` 的二部图拓扑；没有 JSON 时可从
  `target_scope.csv`、`scope.csv` 或 `results.csv` 重建；
- transformation 的重复产物行按 ID 合并，多条 Rule ID 保留为同一化学转换的
  `PredictedReaction` 证据变体；
- 流式读取 P0 实际使用的 RR02 retro 文件，先核对 SHA-256，再提取规则级 EC、
  Diameter、Score、Legacy ID/MNXR 和方向信息；真实 213 MB 文件冒烟通过；
- 当前 RR02 与 Wrapper `score_mode=auto` 一致解析为 `lower_is_better`，原始
  transformation score 和规则 score 均保留，不在 P4 设置经验阈值；
- sink 仅按完整 InChIKey 与 P2 累计 sink 精确匹配，`In Sink` 和名称只作交叉
  检查；伪 sink、结构冲突和规则缺失均产生稳定拒绝原因；
- 使用 compound OR / transformation AND 枚举：一条反应产生多个有效前体时，
  每个末端分支都必须命中可信 sink，未闭合分支不会被压成线性假路径；
- 环路、no-op、iteration/depth 错误、重复路径和枚举上限均有确定性处理；质子等
  无重原子片段只作为辅助片段审计，并将辅因子重建状态标记为 incomplete；
- P4 保持 `Target→…→Sink` 逆合成方向并返回独立反应子图，不修改 P1 schema v1；
  方向翻转、KEGG witness 拼接和分支候选构建由 P5 负责；
- 16 项 P4 离线测试通过；联合 P1–P4、expansion 与电子平衡回归共 85 项通过。

### P5 路线翻转与拼接（2026-08-25 更新）

P5 已在 `src/pathway_analyze/retropath_merge.py` 和
`src/pathway_analyze/retropath_analyze.py` 实现预测反应方向翻转、多 sink expansion
witness 恢复、混合反应 DAG 构建、Top-K 和独立候选文件输出。

P5 验收记录：

- 每个 retrosynthetic `PredictedReaction` 交换 Reaction SMILES、底物和产物，生成
  独立的 biosynthetic RP2 稳定 ID，逆向 ID 不混入合成候选步骤；
- 同一 transformation 的多条 RR02 规则作为一个化学步骤的多个
  `reaction_option_ids` 保留，不复制成多条路线，也不重复计算步骤或新增酶；
- 一条 P4 路径的全部 sink 同时交给 `materialize_frontier_solution()`；depth 0 不生成
  KEGG prefix，depth N 恢复 A0→sink 反应树，多 sink 共享的 KEGG 步骤自动去重；
- KEGG 与 RP2 步骤按前体—产物依赖执行稳定拓扑排序，输出显式
  `depends_on_step_ids`，分支路线不会伪装成线性链；
- 新增 `HybridCandidateRoute`/`HybridCandidateStep` 保存全部 sink、计量、规则选项和
  依赖关系；P1 `CandidateRoute` schema v1 与黄金哈希保持不变；
- 总步骤和新增酶上限同时计入唯一 RP2 transformation；默认 Top-5、每个 witness
  Top-3、最多 10 步和 10 个新增酶；上游或本阶段截断均显式输出；
- `candidate_routes.csv`、`candidate_steps.csv`、`rejected_routes.csv` 使用固定表头、
  UTF-8、LF、稳定排序和逐文件原子替换；P5 本身只产生候选，随后由 P6 编排器统一物化；
- 13 项 P5 离线测试通过；联合 P1–P5、expansion 与电子平衡回归共 98 项通过。

### P6 CLI 与完整流程（2026-08-25 更新）

P6 已在 `src/pathway_analyze/retropath_pipeline.py` 中实现面向用户的
P2→P3→P4→P5 编排入口，并通过 `src/cli/commands/gap.py` 的显式
`--retropath` 开关接入现有 `gap` 命令。

P6 验收记录：

- 不传 `--retropath` 时，dispatcher 直接调用原 `run_gap(config)`；不先运行
  RetroPath、不自动兜底，也不改变原 KEGG 结果；
- 传入 `--retropath` 后只运行 RetroPath 候选流程；depth 0 使用 A0，depth N 使用
  已生成并通过一致性检查的累计 AN，不自动执行 expand；
- P2 输入写入 `depthN/retropath/input/`，P3 原始结果和 P5 三个候选文件写入
  `depthN/retropath/`；随后把全部 Top-K 作为 RetroPath slice 追加到当前 depth 的
  `solutions.csv`、`all_solution_steps.csv` 和电子系统结果；KEGG slice 始终保留；
- `solution_materialization.json` 绑定 P5 候选、RR02、solution 映射及四个正式 CSV；
  相同 KEGG slice 和候选排名会得到稳定、连续的 solution ID；
- `pipeline_result.json` 使用 `retropath_pipeline_result.v3`，记录 sink 来源、任务 ID、
  服务终态、缓存命中、参数、输入/规则/候选文件路径与 SHA-256、scope 和 sink 命中、
  完整路径/候选/拒绝数量以及两阶段截断标记；失败时也原子写入稳定状态和阶段；
- 稳定区分 expansion/规则缺失、输入无效、本地服务不可用、客户端/服务端超时、
  执行失败、解析失败、拼接失败、无候选、source-in-sink 和候选命中；
- `RunConfig` 固定 P2 结构策略、P3 服务/任务/HTTP 参数、P4 枚举上限和 P5 Top-K
  默认值；CLI 只暴露产品级 `--retropath` 开关，专家参数仍集中在运行配置中；
- 7 项 P6 离线测试通过；当前 P1–P6 RetroPath 测试模块共 88 项通过；全仓库回归
  435 项通过、2 项需显式本地服务地址的测试跳过。P6 测试未启动 Docker 或访问实时
  KEGG，真实 KNIME 环境沿用 P0 已完成的服务冒烟记录。

### P7 候选信息展示（2026-08-25 更新）

P7 已在 `src/info_show/retropath_info.py` 实现 RetroPath 运行摘要、候选 DAG
详情和单步下钻，并通过 `info --retropath` 与 `info --retropath-candidate N` 暴露给
用户；同一候选也可通过其物化后的 `info --solution N` 查看。

P7 验收记录：

- `info --retropath -d N` 展示目标、sink 来源、job/service/cache、scope、输入结构、
  sink/path/candidate/rejection 数量、候选摘要、拒绝原因统计和截断风险；
- `info --retropath-candidate N -d N` 按 P5 拓扑顺序展示混合反应 DAG，明确区分
  KEGG expansion prefix 与 RetroPath/RP2 预测步骤；
- 候选步骤保留计量、依赖、Reaction SMILES、规则、来源反应、EC、specificity、
  score 语义、平衡和辅因子恢复状态；追加 `--step M` 可查看单步详情；
- `source_in_sink`、`no_scope` 和失败运行均可读取；失败摘要只依赖
  `pipeline_result.json`，不伪造或要求候选 CSV；
- 成功结果严格校验 `retropath_pipeline_result.v1`、目标、depth、sink 来源、固定
  本地路径、三个候选文件 SHA-256/表头/数量、candidate ID、连续排名和 DAG 依赖；
- 所有含 RP2 步骤的视图明确显示计量、GEM、酶证据和人工复核状态，避免把“可写入”
  误解为“已实验验证”；
- 7 项 P7 离线测试通过；当前 P1–P7 RetroPath 测试模块共 95 项通过；全仓库回归
  442 项通过、2 项需显式本地服务地址的测试按设计跳过。

### P8 计量补全与 strict/relaxed GEM 验证（2026-08-26 更新）

P8 已实现 RR02 来源模板驱动的 RP2 计量/辅因子补全，并将 RetroPath solution 接入
严格 COBRA/GEM 验证。P8 是可选覆盖层：不改变 P5 候选，不改变 solution 数量或编号。

P8 验收记录：

- 固定使用与 RR02 一致的 MNXref v3.0，不根据 EC 或元素差额猜测 NAD(P)、ATP、
  CoA、水、质子或电子；多个来源支持的完整方程保留为独立假设；
- 安装器校验 MetaNetX 官方 MD5并记录 SHA-256，从 213 MB RR02 和约 450 MB
  MNXref 源文件构建紧凑 SQLite 子集；临时 TSV 已删除；
- 真实索引包含 234,384 条 rule-template 链接、16,139 个可用 MNXR、79,147 个
  反应参与项、12,017 个化合物和 84,582 条化合物映射；索引大小约 83 MB；
- 真实规则 `RR-02-fbdda75e23f518b6-02-F` 成功映射 `MNXR94682`，从核心结构
  变化恢复 `MNXM2/H2O` 并得到唯一元素/电荷平衡假设；
- 每个 RP2 步骤最多保留 8 个完整假设，每条候选最多验证 32 个组合，所有裁剪均
  审计；不完整公式、R-group、transport、未平衡或无法对齐的来源模板被拒绝；
- strict GEM 同时要求至少 10% 基线生长、目标正通量以及候选 DAG 每一步按指定方向
  至少 `1e-4` 通量，不开放 generic cofactor sink；relaxed 保持相同计量、增长、目标
  和逐步通量约束，但为路线实际涉及的通用载体开放可审计 sink；
- `validate -s N` 使用统一 solution ID 并按 `solution_source` 自动分流；省略 `-s`
  时验证当前 depth 的全部 KEGG/RetroPath solution；RetroPath 拒绝 pooled/both；
- 输出计量假设、逐项参与物、拒绝原因、严格验证摘要、逐步通量和带哈希 manifest；
  P7 只在输入哈希一致时叠加 P8 状态和定向通量；
- 对选中的候选，P8 将 `passed`/`failed`、计量假设和严格 GEM 证据覆盖到原 solution；
  未被本次选择的候选保持或恢复为 `not_run`，避免复用旧 P8 证据；失败路线仍可写入，
  但显式附带强人工复核警告；
- P1–P8 定向测试 107 项通过；全仓库 454 项通过、2 项本地 Docker 测试按设计跳过。

### P9 预测步骤主酶候选检索（2026-08-25 更新）

P9 首先在隔离候选上实现了预测步骤酶检索；P10 保留其检索核心并废弃了候选直连
CLI。以下为 P9 的历史交付边界：

- P9 曾使用 `main-enzyme --retropath-candidate N --depth D` 直接读取隔离候选；P10
  已改为先 `write --solution N`，再执行无来源参数的 `main-enzyme`；
- 已验证步骤使用 candidate、combination、step、hypothesis 四级证据绑定；未运行 P8
  时使用 candidate、`raw:<candidate_id>`、step 身份，并把同一 step 的多条 RR02 规则
  作为替代检索证据，而不是多个必需酶步骤；
- 只有目标完整计量/结构与来源 MNXR 完全一致时，才将 KEGG/Rhea 映射视为正式反应；
  否则来源 MNXR、EC、Rhea 和 UniProt 仅属于模板证据；
- SelenzymeRF REST client 支持 `smarts` 结构请求，按完整平衡 Reaction SMILES、核心
  Reaction SMILES、RR02 Rule SMARTS 的顺序回退，并以结构和请求参数哈希缓存；
- SelenzymeRF 结构命中即使 reaction similarity 为 1，也固定为 `manual_review`，不会
  标记成已验证酶；缺失、零或非法相似度以及 EC/方向冲突进入审计表；
- 所有组合按覆盖完整性、最差步骤证据、最低相似度确定性排序，只作推荐，不自动选择；
- P9 输出 requirements、Top-N 候选、完整审计、Selenzyme 查询证据和带输入/输出
  SHA-256 的隔离 selection manifest；P10 将相同证据投影到标准主酶 artifacts；
- `info --retropath-candidate N [--step N]` 只在哈希仍有效时显示 P9 组合和候选酶；
  最终人工确认与正式 manifest 晋升留给 P10；
- P1–P9 RetroPath 定向测试 116 项通过、2 项本地 Docker 测试按设计跳过；全仓库
  463 项通过、2 项跳过，真实 RR02/MNXref 完整反应身份冒烟通过。

### P10 统一 solution、manifest 与主酶流程（2026-08-26 调整）

P10 已取消主酶阶段的独立 RetroPath 路线入口，并将“严格通过后才晋升”调整为
“先物化、后可选验证”。后续统一只读取 solution/manifest：

- `gap --retropath` 在 P5 完成后立即把全部 Top-K 候选物化为 solution，保留 KEGG
  原编号并从最大 KEGG 编号后按候选排名追加；重复运行替换 RetroPath slice；
- `solution_materialization.json` 绑定 P5、RR02、候选到 solution 的一一映射以及四个
  正式 CSV 的 SHA-256；部分写入、上游变化或手工修改均 fail closed；
- `validate -s N` 对 RetroPath solution 只将 P8 结果覆盖到这些既有 solution；验证通过、
  失败或未运行都不会增加、删除、拆分或重新编号路线；一个候选即一个 solution，
  多个可行计量组合按确定性顺序选择首个作为当前验证证据，不再 fan-out；
- `info --solution N` 和 `write --solution N` 自动识别 KEGG/RetroPath。预测步骤始终以
  `RP2STEP` 作为稳定反应身份；P8 `RP2STOICH`、完整计量、精确映射和 GEM 结果作为
  可选验证证据挂载；
- `main-enzyme` 删除 `--retropath-candidate` 与 `--depth`。KEGG prefix 复用原检索，
  RP2 suffix 自动使用 P9 来源模板和结构检索，并写入同一标准候选 artifact；
- 公共组合器接受 `manual_review`，质量顺序固定为 `verified > verified_with_risk >
  manual_review`；预测候选可形成完整组合，但组合和 manifest 保留待复核状态；
- 未验证或严格 GEM 失败的预测路线仍可执行 `write --solution`、主酶选择和后续表达
  设计，但 manifest 固定保留 `review_required`、验证状态和显式风险信息；
- RP2 主酶检索在未验证状态即可使用来源 UniProt/EC/Rhea、核心 Reaction SMILES 和
  所有 RR02 Rule SMARTS；P8 通过后可额外使用完整计量 Reaction SMILES/精确映射；
- 2026-08-26 统一 validate 与 strict/relaxed 调整后的全仓库测试为 499 项通过。

### P11.1 隐藏 KEGG 反应恢复评测（2026-08-25 更新）

P11.1 的评测工具、数据契约和 12 例试运行集已实现：

- 固定 EC 1–6 每类 2 条已知 KEGG 单步反应，搜索输入与金标准反应字段物理分离；
- 每例运行 controlled sink 与完整 iML1515 A0 两种 profile，共 24 个核心任务；
- 独立模块复用 P2→P10，输出来源模板、平衡计量、严格 GEM、精确正式 solution 的
  Recall@1/3/5/10、MRR、Wilson 95% 区间、运行时间/IQR 和失败漏斗；
- 模型、培养基、A0、RR02、MNXref、RetroPath runtime、生产源码树和任务 artifacts
  均使用 SHA-256 绑定，支持 fail-closed resume；
- 主酶恢复为可选 full-A0 扩展，首份本地核心基线默认不调用外部酶数据库；
- 回测暴露并修复 RDKit `CalcMolFormula` 末尾电荷后缀无法进入 P8 的问题；
- 终态失败默认可断点保留，显式 `--retry-failed` 才重试；controlled 基础设施失败会
  熔断同案例 full-A0，避免同一目标连续拖垮 Docker/KNIME worker；
- P11.1 定向测试已并入既有 RetroPath 测试模块；全仓库主测试集 482 项通过、2 项
  本地 Docker 测试按设计跳过，另有 93 个子测试通过。

真实 Docker run `20260825T090033Z_4f69c442` 已形成 24 个终态任务记录并固化
`docs/reports/RetroPath P11.1基线报告.md`。其中 controlled/full-A0 可评测任务分别为
7/6，保守 all-selected 精确恢复率均为 1/12（8.3%），可评测分母恢复率分别为
1/7（14.3%）和 1/6（16.7%）。11 个非评测任务包含 7 个 P4 artifact 一致性失败、
3 个服务超时/重启/不可用失败和 1 个配对熔断；这组结果是带失败漏斗的首份基线，
而不是达标结论。P11.2 前应优先处理 P4 对原始结果语义的兼容与服务资源上限。

## 2. 用户工作流

| 使用场景 | 命令 | 行为 |
|---|---|---|
| 原始 KEGG 搜索 | <code>gap --input example.json</code> | 完全保持当前行为 |
| RetroPath 连接 A0 | <code>gap --input example.json --retropath</code> | 使用 A0 生成 sink，不需要 expansion witness |
| 生成深度 3 扩展 | <code>expand --input example.json --depth 3</code> | 生成 A0 到 A3 的累计集合、frontier 和 witness |
| RetroPath 连接 A3 | <code>gap --input example.json --depth 3 --retropath</code> | 使用累计集合 A3 作为 sink，命中后恢复 KEGG witness |
| depth 结果缺失 | 指定 depth 3 和 <code>--retropath</code>，但未运行 expand | 明确报错并提示先运行对应 expand |
| 查看已物化路线 | <code>info --input example.json --solution N --depth 3</code> | `gap --retropath` 后即可查看 Top-K 对应的统一 solution |
| 可选验证 | <code>validate --input example.json -s N --mode per --depth 3</code> | 按 solution 来源自动分流并更新计量/GEM 状态，不改变编号 |
| 直接写入预测路线 | <code>write --input example.json --solution N --depth 3</code> | 无论 P8 未运行、通过或失败，均使用同一接口并保留风险 provenance |
| 预测路线主酶 | <code>main-enzyme --input example.json</code> | 从 manifest 自动分流 KEGG 与 RetroPath Step，不再指定候选或 depth |

原始搜索无解时，建议提示：

    未找到 KEGG 路线。

    可直接尝试 A0：
    gap --input example.json --retropath

    或先扩大可信边界：
    expand --input example.json --depth 3
    gap --input example.json --depth 3 --retropath

## 3. 集合与路线定义

| 符号 | 定义 |
|---|---|
| A0 | GEM/FBA 在当前模型、培养基和生长约束下直接证明可生成的 KEGG 化合物集合 |
| Fn | 第 n 层通过一个定向 KEGG 反应新发现的 frontier 化合物集合；电子载体不作为主底物或 frontier 产物 |
| An | 截至第 n 层的累计可达集合，An = A(n-1) ∪ Fn |
| X1…Xm | 一条完整 RetroPath 分支路径命中的全部 sink 边界化合物 |
| KEGG prefix | expansion witness 恢复出的 A0→X1…Xm 合成反应树 |
| RetroPath suffix | RetroPath 逆向结果翻转后得到的 X1…Xm→目标预测 DAG |

完整候选路线：

    A0 --KEGG known reaction DAG--> X1…Xm --RetroPath predicted DAG--> Target

当所有 sink 的 depth = 0 时，它们均属于 A0，KEGG prefix 为空。

`expand` 使用 component-step-aware v3：所有普通主底物仍须可达；明确识别的电子载体
只生成风险、净变化和辅助系统需求，不阻断主产物，也不进入 sink。深度按反应数而非唯一
EC/KO 数计算；无 EC/KO 但结构完整的 KEGG 原子反应继续允许。多步路线中明确标记为
`first/second/... step` 的独立组件反应按一层扩展，汇总反应仍被拒绝。v1/v2 manifest
视为过期，必须重新生成。

## 4. 总体修改计划表

| 阶段 | 优先级 | 目标 | 主要修改 | 文件落点 | 主要输出 | 测试与验收 | 前置条件/权限 | 状态 |
|---|---:|---|---|---|---|---|---|---|
| P0 资源与版本约定 | P0 | 固定可复现的本地服务环境 | 构建本地 Docker HTTP 服务，固定 Wrapper 3.9.1、KNIME 4.7.0、workflow r20260212、RR02 哈希、单 Worker 和超时策略 | services/retropath；compose.retropath.yml；docs | 健康检查、异步任务、原始结果、日志与 run manifest | 镜像和规则缺失时给出可操作错误；真实 KNIME 冒烟通过 | 已授权服务目录；RR02 已安装 | 已完成 |
| P1 预测数据模型 | P0 | 建立非 KEGG 反应/化合物表示 | 定义预测化合物、预测反应、候选路线、运行 provenance、运行结果及 RP2 命名规则 | 新增 src/pathway_analyze/retropath_models.py | 可序列化数据对象 | 17 项离线测试通过；ID 黄金值稳定、可 JSON 往返、相同输入产生相同哈希 | 已授权目录 | 已完成 |
| P2 结构与输入生成 | P0 | 自动生成 source/sink | KEGG MOL 校验缓存、RDKit 标准化、depth 0/累计 depth sink、映射表和拒绝表 | 新增 retropath_structure.py、retropath_input.py；锁定 RDKit 2026.3.5 | target_source.csv、chassis_sink.csv、compound_mapping.csv、rejected_compounds.csv | 21 项 P2 离线测试通过；真实 KEGG 冒烟通过 | 一次性依赖与 .gitignore 授权已使用 | 已完成 |
| P3 RetroPath client | P0 | 稳定调用本地服务 | loopback HTTP 提交、轮询、恢复、超时、状态映射、artifact 下载、manifest 校验和健康一致缓存 | 新增 retropath_client.py；锁定 httpx 0.28.1 | raw results/scope、服务日志、client state、P1 run result 和审计 manifest | 14 项离线协议测试通过；区分正常终态、服务错误和客户端错误 | P0 协议已验证；真实服务冒烟沿用 P0 记录 | 已完成 |
| P4 网络解析与路径枚举 | P0 | 从预测网络得到完整路径 | 解析 transformation、结构、sink 命中、rule、EC、specificity、score；以内置 AND/OR 等价枚举器生成完整分支路径 | 新增 retropath_parser.py、retropath_routes.py | 逆向候选路径与拒绝原因 | 16 项 P4 测试通过；完整 InChIKey sink 闭合、环路、重复与上限处理正确 | P3 | 已完成 |
| P5 路线翻转与拼接 | P0 | 构建完整混合候选路线 | 将 Target→sink 分支图翻转为 sink→Target；恢复全部 expansion witness；合并共享步骤并限制 Top-K | 新增 retropath_merge.py、retropath_analyze.py；复用 materialize_frontier_solution | candidate_routes.csv、candidate_steps.csv、rejected_routes.csv | 13 项 P5 测试通过；depth 0 无 prefix、depth N 多 sink DAG 和方向正确、不伪造 KEGG ID | P2、P4 | 已完成 |
| P6 CLI、编排与物化 | P0 | 暴露显式开关并生成统一路线 | gap 增加 <code>--retropath</code>，默认 False；增加本地服务、规则、步数和超时配置；不自动 expand；P5 后物化全部 Top-K | 新增 retropath_pipeline.py、retropath_materialization.py；修改 src/cli/commands/gap.py、src/config/run_config.py | 两种搜索方式、pipeline_result.json、solution_materialization.json、统一 solution | 不加参数时直接调用原 KEGG 入口；depth 0/N 符合约定；KEGG slice 保留且 Top-K 编号稳定 | 单次授权已使用 | 已完成 |
| P7 候选信息展示 | P1 | 让用户看懂命中与风险 | 增加独立 info 摘要、候选 DAG 和单步视图；显示命中 Cxxxxx、depth、KEGG prefix、RP2 suffix、规则证据、拒绝原因和验证风险 | 新增 retropath_info.py；扩展 info CLI | 中文 JSON 摘要、候选详情和单步详情 | 7 项 P7 测试通过；校验 schema、目标/depth、SHA-256、数量和 DAG 关系 | P5、P6；单次 CLI/.gitignore 授权已使用 | 已完成 |
| P8 计量与 GEM 验证 | P1 | 以 strict 或 relaxed 判断完整路线可行性 | 固定 MNXref v3.0；按 RR02 来源模板恢复共底物/辅因子；校验分子式、电荷和平衡；强制完整候选 DAG 同时承载通量；relaxed 只开放路线涉及的通用载体 | 新增 mnxref/stoichiometry/retropath GEM 模块；扩展 validate 和 P7 | 计量假设、参与项、模式、开放载体、验证与逐步通量 | strict/relaxed 对照、provenance、防篡改和全仓库回归通过 | P5–P7；数据/CLI/.gitignore 单次授权已使用 | 已完成 |
| P9 SelenzymeRF 与主酶选择 | P1 | 为预测混合路线生成候选酶 | raw step 可直接使用全部 RR02 替代规则证据；P8 后增加 hypothesis/完整计量证据；精确反应、来源模板、完整/核心 Reaction SMILES、Rule SMARTS 分级检索；结构命中只供人工复核 | 新增 RetroPath enzyme selection；扩展 Selenzyme client、main-enzyme CLI 和 P7 | requirements、Top-N/审计候选、Selenzyme 证据、带哈希 selection manifest | 未验证/已验证身份均可进入统一主酶检索；多规则不重复计算必需步骤；相似度 1 不误判 | P8 可选；CLI/.gitignore 单次授权已使用 | 已完成 |
| P10 统一流程 | P1 | 将全部 Top-K 混合路线纳入现有 solution、manifest 和主酶流程 | P5 后立即物化；P8 只覆盖状态/证据；materialization 哈希提交；manifest 自动分流步骤；公共组合器支持 manual_review | src/pathway_analyze、src/write_manifest、src/main_protein_selection、src/info_show | 统一混合 solution、可选验证覆盖、统一主酶 artifacts 和带 pending review 的 manifest | KEGG 编号和路线不变；KEGG/RetroPath 未验证、通过、失败均可写；solution ID 不变；篡改 fail closed；全仓库回归通过 | P5、P9，GEM/P8 验证均可选 | 已完成 |
| P11 回测与阈值校准 | P1 | 量化假阳性和收益 | 隐藏已知 KEGG 反应做恢复测试；排除来源规则做 promiscuity 测试；加入青蒿素非酶促边界案例 | tests、docs | 回测报告和参数建议 | top-k 恢复、平衡、GEM、酶证据通过率可复现 | P5 起可分批实施 | 进行中：P11.1 基线已完成，P11.2+ 待实施 |

## 5. 第一批应实施的文件

第一批只建立数据契约和输入，不修改主酶选择或 GEM：

    src/pathway_analyze/
    ├── retropath_models.py
    ├── retropath_structure.py
    └── retropath_input.py

    tests/
    ├── test_retropath_models.py
    ├── test_retropath_structure.py
    └── test_retropath_input.py

P2 实际核心接口：

    build_retropath_inputs(
        target_compound_id,
        expansion_bundle,
        structure_provider,
        output_dir,
    )

第一批完成标准：

- depth 0 能从 producible_kegg_compounds.csv 稳定生成 sink；
- depth N 能从累计 ExpansionBundle.reachable_compounds 生成 sink；
- 每个结构保留 Cxxxxx、最小 depth、InChI、InChIKey、SMILES 和来源；
- 无结构、解析失败、重复和结构冲突均有审计记录；
- 测试不依赖实时 KEGG 或真实 RetroPath 程序。

P2 已满足以上完成标准。P6 调用前先使用现有 `load_expansion_bundle()` 读取并校验
指定 depth；P2 不会自动生成缺失或过期的 expansion。

## 6. 内部数据契约

### 6.1 预测化合物

| 字段 | 说明 |
|---|---|
| compound_id | KEGG Cxxxxx、RP2CPD:InChIKey 或 RP2CPD:完整 SHA-256 |
| name | 化合物名称，可为空但不得伪造 |
| inchi | 标准 InChI |
| inchikey | 结构映射与去重键 |
| isomeric_smiles | 保留立体化学的 SMILES |
| formula | 分子式，供平衡检查 |
| charge | 电荷，供平衡检查 |
| kegg_ids | 同一结构对应的全部 KEGG ID |
| minimum_depth | 在累计扩展集合中的最小深度；新中间体为空 |
| structure_provenance | KEGG MOL、缓存版本和转换工具版本 |

### 6.2 预测反应

| 字段 | 说明 |
|---|---|
| reaction_id | RP2:哈希 |
| reaction_source | retropath |
| evidence_type | rule_predicted；映射到已知反应后可改为 database_exact |
| rule_id | RetroRules rule/template ID |
| reaction_smiles | 具体底物和产物预测转换 |
| substrate_compounds | 预测底物结构标识 |
| product_compounds | 预测产物结构标识 |
| orientation | retrosynthetic 或 biosynthetic；作为反应身份的一部分 |
| source_reaction_ids | 规则来源 MNXR/Rhea/KEGG 等反应 |
| source_ec_numbers | 规则关联 EC，可为空或不完整 |
| source_uniprot_ids | 规则关联序列，可为空 |
| rule_specificity | 非负 radius/diameter 数值，可在后续解析阶段补充 |
| rule_specificity_semantics | radius 或 diameter；与 rule_specificity 同时出现 |
| rule_score_raw | 原始分数 |
| score_semantics | higher_is_better 或 lower_is_better；与原始分数同时出现 |
| balance_status | 后续计算的元素/电荷平衡状态 |
| cofactor_reconstruction_status | 共底物和辅因子恢复状态 |

### 6.3 P1 单边界候选路线

以下 `CandidateRoute` 是 P1 schema v1 的兼容契约，只能表达单 sink 线性摘要。P5
不会把多 sink 分支强制压入该模型。

| 字段 | 说明 |
|---|---|
| candidate_id | RP2ROUTE:完整 SHA-256 稳定候选路线 ID |
| target_compound_id | 路线目标，允许 KEGG Cxxxxx 或 RP2CPD 标识 |
| matched_sink_kegg_id | RetroPath 命中的可信边界 Cxxxxx |
| matched_sink_depth | 边界最小 depth |
| kegg_prefix_reaction_ids | 按合成方向排列的 KEGG Rxxxxx 前缀 |
| retropath_reaction_ids | 按合成方向排列的 RP2 预测反应后缀 |
| kegg_prefix_steps | 由前缀反应列表计算的 expansion witness 步数 |
| retropath_steps | 由后缀反应列表计算的预测步数 |
| total_steps | 拼接后总步数 |
| route_source | kegg_retropath |
| contains_predicted_steps | True |
| minimum_rule_specificity | 路线上最差预测规则的特异性 |
| validation_status | raw、structure、stoichiometry、gem、enzyme、promoted、rejected |
| review_required | 含预测步骤时默认为 True |
| rejection_reasons | 所有硬门禁失败原因 |

### 6.4 P5 混合分支候选

| 字段 | 说明 |
|---|---|
| candidate_id | RP2ROUTE:完整 SHA-256，身份包含全部 sink、步骤和依赖 |
| source_retrosynthetic_path_id | P4 逆合成路径来源 |
| sink_matches | 全部代表 Cxxxxx、InChIKey、别名和最小 depth |
| steps | 合成方向拓扑步骤，包含 KEGG expansion 和 RP2 prediction |
| reaction_option_ids | KEGG 步骤一个 Rxxxxx；RP2 步骤保留同一 transformation 的全部规则变体 |
| depends_on_step_ids | 分支 DAG 的直接上游步骤 |
| substrate/product_stoichiometry | 已知计量；RP2 当前保留结构出现次数，P8 再恢复辅因子和平衡 |
| validation_status | P5 输出固定为 raw |
| review_required | 含预测步骤时固定为 True |

### 6.5 运行 provenance 与结果封装

`RetroPathRuntimeProvenance` 保存 Wrapper 实际固定版本、Wrapper 自报版本、workflow、
KNIME、RDKit 插件、RetroRules 版本及规则 SHA-256。`RetroPathRunResult` 使用 schema
version 2 封装任务 ID、运行状态、返回码、稳定 failure code、参数、产物、预测实体、
候选路线和错误。
状态包含 `queued`、`running`、`succeeded`、`no_solution`、`source_in_sink`、
`failed` 和 `timed_out`。

ID 哈希统一使用 UTF-8 编码的规范 JSON：`sort_keys=True`、紧凑分隔符、
`ensure_ascii=True`，随后计算完整小写 SHA-256。未来若修改任何身份字段或哈希输入，
必须升级 schema version，不得静默改变 v1 ID。

## 7. 运行输出目录

    kegg_gap_<target>/depth_<N>/retropath/
    ├── input/
    │   ├── target_source.csv
    │   ├── chassis_sink.csv
    │   ├── compound_mapping.csv
    │   └── rejected_compounds.csv
    ├── raw/
    │   ├── service_results.json
    │   ├── service_run_manifest.json
    │   ├── results.csv
    │   ├── target_scope.csv
    │   ├── target_scope.json（workflow 提供时）
    │   ├── stdout.log
    │   └── stderr.log
    ├── client_state.json
    ├── candidate_routes.csv
    ├── candidate_steps.csv
    ├── rejected_routes.csv
    └── run_manifest.json

run_manifest.json 至少记录：

- 目标 Cxxxxx 和 requested depth；
- sink 来源（A0 或累计 AN）；
- A0/AN 化合物数量；
- 成功写入 sink 的结构数量；
- 结构拒绝数量和原因；
- RetroPath 程序版本和哈希；
- RetroRules 版本、哈希和评分语义；
- max_steps、timeout 等参数；
- 输入文件哈希；
- 执行状态和退出码；
- scope 是否存在；
- 命中 sink 和候选路线数量。

P2 CSV 格式固定如下：

- `target_source.csv`、`chassis_sink.csv`：`Name,InChI`，Name 使用目标或代表
  Cxxxxx；
- `compound_mapping.csv`：记录 role、原始 Cxxxxx、代表 Cxxxxx、最小 depth、
  InChI、InChIKey、立体 SMILES、分子式、电荷和结构 provenance；
- `rejected_compounds.csv`：记录 role、Cxxxxx、最小 depth、稳定原因代码和详情。

CSV 使用 UTF-8、LF 和稳定排序；source 必须恰好一行，sink 必须至少一行，且均
已通过 P0 本地服务的输入校验器测试。

P3 的 `run_manifest.json` 使用 `retropath_client_run.v1`，包含 request fingerprint、
P2 输入 SHA-256、目标与 expansion depth、任务参数、健康检查、终态 job、远端到
本地 artifact 映射、所有本地 artifact SHA-256 以及 P1 `RetroPathRunResult`。
`client_state.json` 使用 `retropath_client_state.v1`，在提交前、提交后、轮询、终态、
失败和完成阶段原子更新，用于安全恢复。

P5 返回三个候选文件的 SHA-256、候选/拒绝数量和 `RetroPathMergeResult`；P6 再把
这些字段汇总到面向用户的完整运行结果，不回写或覆盖 P3 原始运行清单。

## 8. 假阳性分级、验证与人工复核

| 门禁 | 硬性要求 | 未通过处理 |
|---|---|---|
| 输入结构 | 目标和 sink 均能解析，立体化学可追溯 | 不进入 RetroPath，或从 sink 排除并记录 |
| sink 身份 | 完整 InChIKey 精确命中优先；仅缺失立体层可保留并强制人工复核；明确构型冲突拒绝 | 降级或路线拒绝 |
| 规则来源 | 默认只使用生化规则，不使用 USPTO 有机合成规则 | 路线拒绝或仅探索展示 |
| 规则特异性 | 优先高特异规则；radius/diameter 按版本解释 | 低特异路线降级并要求人工复核 |
| 反应有效性 | 非 no-op、结构可解析、方向明确 | 步骤拒绝 |
| 共底物/辅因子 | P8 严格验证时必须可恢复；raw solution 明确标记 core_only | P8 失败，原预测路线保留 |
| 元素/电荷平衡 | P8 严格 GEM 前必须通过 | P8 失败，原预测路线保留 |
| GEM | strict_l1 或 relaxed 可分别标记通过，验证模式必须保留 | relaxed 通过仍附带载体开放风险；failed/not_run 也可写 |
| 酶证据 | 精确反应、规则来源或 SelenzymeRF 结构相似检索至少一种有效 | 无候选时主酶组合不完整，路线仍保留 |
| 非酶促步骤 | 允许标记 nonenzymatic，不得强制选主酶 | 进入工艺/人工复核 |
| 人工复核 | 含 RP2 步骤时默认需要 | 可继续设计，但不得表述为已复核或已实验验证 |

路线排序使用“最差步骤优先”，建议顺序：

1. 无硬门禁失败；
2. RetroPath 预测步数更少；
3. 命中 depth 更低；
4. 最低规则特异性更高；
5. 最低酶证据等级更高；
6. 总步骤和辅因子负担更低；
7. 保留结构与反应多样性，避免 Top-K 全是同类变体。

## 9. SelenzymeRF 接入策略

| RP2 步骤情况 | 检索策略 | 证据等级 |
|---|---|---|
| 可精确映射到 KEGG/Rhea | 复用现有精确反应检索 | 高 |
| 无正式 ID，但有规则来源 EC/UniProt | 优先来源序列和完整 EC，再用 SelenzymeRF 补充 | 中高/中 |
| 无完整 EC/来源序列，但有底物和产物结构 | 使用 substrate_smiles>>product_smiles 查询 SelenzymeRF | 预测性证据 |
| 无明确结构或 SelenzymeRF 无候选 | 无法自动选择主酶 | 阻断 |
| 非酶促步骤 | 不进入主酶选择 | 工艺/人工步骤 |

P9 已为 src/main_protein_selection/selenzyme_retrieval.py 增加 Reaction SMILES/SMARTS 查询，并保留：

- reaction_similarity；
- sim_RF 和 sim_2018；
- matched_reaction_id；
- 候选 EC/UniProt；
- host taxonomic distance；
- 方向证据；
- substrate specificity unverified 标记。

结构查询通过 JSON `smarts` 字段调用现有 `/REST/Query`；对 RP2 步骤，任何结构相似
结果都只进入人工复核候选，不使用任意相似度阈值自动验证。只有正且合法的相似度、
有效 UniProt 记录和无明确 EC/方向冲突的候选才会进入 Top-N。

## 10. 测试计划表

| 测试类别 | 关键场景 | 预期结果 |
|---|---|---|
| 回归测试 | 不指定 <code>--retropath</code> | KEGG 输出和路径不变；GEM 验证作为可选证据，未验证或失败仍可写入 |
| depth 0 | A0 有合法和非法结构 | sink 只包含合法、去重后的 A0，失败项进入 rejection |
| 累计 depth | A0、F1、F2 有交叉底物 | depth 2 sink 包含 A0∪F1∪F2，最小 depth 正确 |
| 电子载体扩展 | 主底物可达且反应含 P450/ferredoxin 等载体 | 主产物进入下一层；载体不进入 sink；辅助角色和风险完整保留 |
| KEGG 组件步骤 | `first/second step of three-step reaction` 与对应汇总反应同时存在 | 组件步骤各计一层；汇总反应不进入 frontier |
| 缺失 expand | 指定 depth 3 + RetroPath，但无 depth 3 结果 | 提示先执行 expand |
| 结构映射 | 多个 Cxxxxx 共享一个结构 | mapping 保留全部 ID，按最小 depth 选默认边界 |
| 无 scope | RetroPath 成功但未命中 sink | 返回 retropath_no_scope，不视为执行失败 |
| 超时/异常 | runner 超时或非零退出 | 保存日志并返回独立错误状态 |
| 路线方向 | 原始网络为 Target→M→X | 输出合成方向 X→M→Target |
| depth 0 拼接 | 命中 X∈A0 | prefix 为空，suffix 正确 |
| depth N 拼接 | 命中 depth 2 的 X，存在多个 witness | 生成受 Top-K 控制的多个候选 |
| ID 安全 | 路线包含预测反应/中间体 | 只使用 RP2/RP2CPD，不伪造 R/C 编号 |
| 可选计量验证 | 无法恢复辅因子或不平衡 | 标记 P8 failed，保留原 solution 并要求人工复核 |
| solution 稳定性 | P8 未运行、通过、失败或只验证部分候选 | solution 数量和编号不变，状态/证据按选择更新 |
| SelenzymeRF | 无 EC 但有 Reaction SMILES | 可生成相似反应候选，并标记预测性证据 |
| 非酶促案例 | 青蒿素末端光氧化 | 不强制选主酶，标记工艺/人工复核 |

## 11. 用户可见状态

原始模式：

    {
      "retropath_requested": false,
      "search_engine": "kegg"
    }

RetroPath depth 0：

    {
      "retropath_requested": true,
      "search_engine": "retropath",
      "expansion_depth": 0,
      "sink_source": "chassis_A0"
    }

RetroPath depth N：

    {
      "retropath_requested": true,
      "search_engine": "retropath",
      "expansion_depth": 3,
      "sink_source": "cumulative_expansion_A3"
    }

执行结果至少区分：

    retropath_candidates_found
    retropath_target_already_reachable
    retropath_no_scope
    retropath_source_in_sink
    retropath_input_invalid
    retropath_expansion_missing
    retropath_rules_missing
    retropath_service_unavailable
    retropath_timeout
    retropath_execution_failed
    retropath_parse_failed
    retropath_merge_failed
    retropath_configuration_invalid

## 12. 权限与依赖

### 当前已授权目录内可完成

- src/pathway_analyze
- src/main_protein_selection
- src/protein_selection
- src/info_show
- src/write_manifest
- docs
- tests

### 实施前需要单次授权

| 路径 | 用途 |
|---|---|
| src/cli/commands/gap.py | 增加 <code>--retropath</code> |
| src/config/run_config.py | 增加程序、规则、步数、超时等默认配置 |

P2 已使用一次性授权修改 `pyproject.toml`、`uv.lock` 和 `.gitignore`：固定
`rdkit==2026.3.5`，并仅放行 P1/P2 的三个 RetroPath 测试文件。

P3 已使用新的一次性授权修改相同根文件：固定 `httpx==0.28.1`，并仅额外放行
`tests/test_retropath_client.py`。

P4、P5 分别使用一次性 `.gitignore` 授权，仅放行各阶段新增的 parser/routes 和
merge/analyze 测试文件；未修改依赖或其他根目录配置。

P6 已使用本阶段单次授权修改 `src/cli/commands/gap.py`、
`src/config/run_config.py` 和 `.gitignore`；根目录文件仅额外放行
`tests/test_retropath_pipeline.py`，未修改项目依赖。

P7 已使用本阶段单次授权修改 `src/cli/commands/info.py` 和 `.gitignore`；根目录文件
仅额外放行 `tests/test_retropath_info.py`，未修改项目依赖。

P8 已使用本阶段单次授权写入 `data/retropath/mnxref/3.0/`，修改
`src/cli/commands/validate.py` 和 `.gitignore`；MNXref 原始 TSV 未保留且数据目录
继续忽略，Git 只放行三个 P8 测试文件，未增加项目依赖。

P9 已使用本阶段单次授权修改 `src/cli/commands/main_enzyme.py` 和 `.gitignore`；根目录
仅放行 `tests/test_retropath_enzyme_selection.py`。复用现有远端 SelenzymeRF 配置，
未修改 `services/`、`data/`、`.env` 或项目依赖。

P10 已获得一次性授权修改 `src/cli/commands/main_enzyme.py`、`src/cli/commands/write.py`
和 `.gitignore`；实际按最小修改原则只删除 main-enzyme 的临时候选/depth 参数，并在
`.gitignore` 放行 `tests/test_retropath_promotion.py`，未修改 write CLI 或项目依赖。

外部 RetroPath 可执行文件和 RetroRules 规则包的落盘目录也需要在实施前确定。不建议将大型二进制和规则数据直接提交到 Git 仓库。

## 13. 推荐实施顺序

- [x] P0：确定 RetroPath、RetroRules 和 RDKit 版本及存放策略
- [x] P1：定义 RP2 数据模型和稳定哈希 ID
- [x] P2：完成目标 source、A0/AN sink 和 mapping 生成
- [x] P2：完成离线单元测试，不依赖实时网络
- [x] P3：实现 HTTP client、超时、恢复、日志、退出状态和缓存
- [x] P4：解析 scope 并枚举命中 sink 的完整路径
- [x] P5：翻转 RetroPath 路线并拼接 expansion witness
- [x] P5：输出独立候选文件，由 P6 在同次 `gap --retropath` 中物化全部 Top-K
- [x] P6：加入 <code>--retropath</code>，验证默认 KEGG 回归不变
- [x] P7：增加候选路线信息展示
- [x] P8：补全计量并接 strict GEM
- [x] P9：扩展主酶选择，加入 SelenzymeRF Reaction SMILES 查询
- [x] P10：全部 Top-K 物化为统一 solution，P8 改为可选覆盖层，并接入 manifest 驱动的主酶完整流程
- [x] P11.1：完成评测设计、12 例数据集、双 profile 运行器和统计报告器
- [x] P11.1：形成 24 个真实任务终态记录并固化带失败漏斗的基线报告
- [x] P11.1 修复：分级结构身份、target/sink collision、重复产物计量、P8 立体缺失映射
- [x] P11.1 修复：服务 30 分钟硬超时、6 GiB 连续采样熔断、稳定 failure_code 与 telemetry
- [x] P11.1 修复后：重新完成 24 任务回测并固化新版漏斗报告
- [ ] P11.2+：完成来源规则排除、promiscuity 和青蒿素等回测

## 14. 第一里程碑完成定义

第一里程碑覆盖 P0–P6 的候选路线闭环，完成条件：

1. 默认 gap 行为和现有结果完全不变；
2. <code>--retropath</code> 可以在 depth 0 运行；
3. depth N 使用累计集合并要求对应 expand 结果存在；
4. 自动生成 source、sink、mapping 和 rejection；
5. RetroPath 可执行、无解、超时和失败状态可区分；
6. 能从 scope 得到完整命中路径；
7. 能映射命中 Cxxxxx、翻转预测段并恢复 KEGG prefix；
8. 所有预测反应和中间体使用 RP2 命名；
9. 输出可审计候选路线、步骤、日志和 run manifest；
10. 全部 Top-K 可进入统一 solution/manifest/主酶流程，但未验证或验证失败状态不得
    丢失，且必须要求人工复核。

## 15. P11.1 修复后回测结论（2026-08-26）

- Run：`20260826T021333Z_2d21764c`，24/24 任务有终态；
- 20 个任务可评测，2 个任务以真实 `resource_exhausted` 失败，2 个配对 full-A0
  按基础设施熔断策略跳过；Docker 服务保持健康；
- 可评测任务中，20/20 原始 scope 出现金标准 RR02 来源规则，18/20 形成金标准连通
  候选，6/20 为完整立体身份精确恢复，8/20 恢复金标准平衡计量，4/20 通过严格
  GEM，2/20 精确晋升为正式 solution；
- C00900 重复底物案例由旧的 P8 reconstruction error 修复为 controlled/full-A0
  均 `balanced_gold_rank=1`；
- `ec3_r00913` 在 controlled/full-A0 中均完整 `formal_exact_recovered`；
- 报告：`docs/reports/RetroPath P11.1修复后评测报告.md`。
