# RetroPath 接入修改计划

> 文档状态：待实施  
> 制定日期：2026-08-24  
> 适用项目：GLADE  
> 目标：为现有 KEGG 通路搜索增加由用户显式启用的 RetroPath 预测搜索，同时隔离并审计预测反应，避免其未经验证进入正式设计流程。

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
| 完整路线 | A0→边界 X 的 KEGG witness 与 X→目标的 RetroPath 预测段拼接 |
| 预测反应 ID | 使用 RP2:哈希，禁止伪造 Rxxxxx |
| 预测中间体 ID | 使用 RP2CPD:InChIKey或结构哈希，禁止伪造 Cxxxxx |
| 第一阶段结果 | 只输出候选路线，不写入正式 solutions.csv、manifest 或主酶选择 |
| 正式晋升 | 必须通过结构、计量、GEM、酶证据和人工复核门禁 |

## 2. 用户工作流

| 使用场景 | 命令 | 行为 |
|---|---|---|
| 原始 KEGG 搜索 | <code>gap --input example.json</code> | 完全保持当前行为 |
| RetroPath 连接 A0 | <code>gap --input example.json --retropath</code> | 使用 A0 生成 sink，不需要 expansion witness |
| 生成深度 3 扩展 | <code>expand --input example.json --depth 3</code> | 生成 A0 到 A3 的累计集合、frontier 和 witness |
| RetroPath 连接 A3 | <code>gap --input example.json --depth 3 --retropath</code> | 使用累计集合 A3 作为 sink，命中后恢复 KEGG witness |
| depth 结果缺失 | 指定 depth 3 和 <code>--retropath</code>，但未运行 expand | 明确报错并提示先运行对应 expand |

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
| Fn | 第 n 层通过 KEGG 反应新发现的 frontier 化合物集合 |
| An | 截至第 n 层的累计可达集合，An = A(n-1) ∪ Fn |
| X | RetroPath 命中的 sink 边界化合物 |
| KEGG prefix | expansion witness 恢复出的 A0→X 路线 |
| RetroPath suffix | RetroPath 逆向结果翻转后得到的 X→目标预测路线 |

完整候选路线：

    A0 --KEGG known reactions--> X --RetroPath predicted reactions--> Target

当 depth = 0 时，X 属于 A0，KEGG prefix 为空。

## 4. 总体修改计划表

| 阶段 | 优先级 | 目标 | 主要修改 | 文件落点 | 主要输出 | 测试与验收 | 前置条件/权限 | 状态 |
|---|---:|---|---|---|---|---|---|---|
| P0 资源与版本约定 | P0 | 固定可复现的外部执行环境 | 确定 RetroPath 可执行文件、RetroRules RR02、版本、SHA-256、超时和缓存策略 | docs；运行时配置 | 资源清单、版本与哈希 | 程序和规则缺失时给出可操作错误 | 外部资源目录需确定 | 待办 |
| P1 预测数据模型 | P0 | 建立非 KEGG 反应/化合物表示 | 定义预测化合物、预测反应、候选路线、运行结果及 RP2 命名规则 | 新增 src/pathway_analyze/retropath_models.py | 可序列化数据对象 | ID 稳定、结构字段完整、相同输入产生相同哈希 | 已授权目录 | 待办 |
| P2 结构与输入生成 | P0 | 自动生成 source/sink | KEGG MOL 获取、MOL→InChI/InChIKey/SMILES、depth 0/累计 depth sink、映射表和拒绝表 | 新增 retropath_structure.py、retropath_input.py | target_source.csv、chassis_sink.csv、compound_mapping.csv、rejected_compounds.csv | depth 0 等于 A0；depth N 包含 A0 和全部 frontier；立体化学与去重正确 | RDKit 需修改 pyproject.toml，需单次授权 | 待办 |
| P3 RetroPath runner | P0 | 稳定调用外部程序 | 参数列表调用、超时、退出码、日志、隔离目录、输入哈希缓存 | 新增 retropath_runner.py | raw results/scope、日志、run manifest | 区分成功、有结果、无 scope、超时、程序缺失和异常退出 | RetroPath 与规则可用 | 待办 |
| P4 网络解析与路径枚举 | P0 | 从预测网络得到完整路径 | 解析 transformation、结构、sink 命中、rule、EC、specificity、score；通过 RP2Paths 适配器或等价枚举器得到路径 | 新增 retropath_parser.py、retropath_routes.py | 逆向候选路径与拒绝原因 | 只保留真正命中 sink 的路径；环路与重复处理正确 | P3 | 待办 |
| P5 路线翻转与拼接 | P0 | 构建完整混合候选路线 | Target→X 翻转为 X→Target；恢复 expansion witness；组合多个 witness 与预测路径；限制 Top-K | 新增 retropath_merge.py、retropath_analyze.py；复用 materialize_frontier_solution | candidate_routes.csv、candidate_steps.csv、rejected_routes.csv | depth 0 不生成 prefix；depth N prefix 和方向正确；不伪造 KEGG ID | P2、P4 | 待办 |
| P6 CLI 与运行配置 | P0 | 暴露显式开关 | gap 增加 <code>--retropath</code>，默认 False；增加程序、规则、步数和超时配置；不自动 expand | 修改 src/cli/commands/gap.py、src/config/run_config.py | 两种用户搜索方式和审计字段 | 不加参数时现有行为及结果不变；depth 校验符合约定 | 两目录需单次授权 | 待办 |
| P7 候选信息展示 | P1 | 让用户看懂命中与风险 | 显示命中 Cxxxxx、depth、KEGG prefix、RP2 suffix、规则证据和拒绝原因 | src/info_show | 候选路线摘要 | 清楚区分 KEGG 与预测步骤 | P5 | 待办 |
| P8 计量与 GEM 验证 | P1 | 判断完整路线是否严格可行 | 恢复共底物/辅因子；分子式、电荷和平衡；GEM 从本地预测反应记录读取计量；运行 strict_l1 | 修改 gem_validation.py | 结构、计量和 GEM 验证结果 | 不平衡或辅因子不完整的路线禁止晋升；relaxed 不作为正式通过 | P5 | 待办 |
| P9 SelenzymeRF 与主酶选择 | P1 | 为 RP2 步骤生成候选酶 | 支持 namespaced ID；优先正式反应映射，其次规则来源 EC/UniProt，最后用 Reaction SMILES/SMARTS 查询 SelenzymeRF | src/main_protein_selection、src/protein_selection | 候选酶及相似反应证据 | 保存 reaction similarity、sim_RF、匹配反应、方向和风险；无结构/无命中时阻断 | P8 | 待办 |
| P10 manifest 晋升 | P2 | 将通过门禁的混合路线纳入正式设计 | 增加预测 provenance、验证状态、人工复核状态；只有 promoted 路线可写正式 manifest | src/write_manifest、src/info_show | 正式混合路线 manifest | 未验证预测路线不能进入表达设计 | P8、P9 | 待办 |
| P11 回测与阈值校准 | P1 | 量化假阳性和收益 | 隐藏已知 KEGG 反应做恢复测试；排除来源规则做 promiscuity 测试；加入青蒿素非酶促边界案例 | tests、docs | 回测报告和参数建议 | top-k 恢复、平衡、GEM、酶证据通过率可复现 | P5 起可分批实施 | 待办 |

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

建议核心接口：

    build_target_source(target_compound_id, structure_provider, output_dir)

    build_chassis_sink(
        base_a0_path,
        expansion_bundle,
        depth,
        structure_provider,
        output_dir,
    )

第一批完成标准：

- depth 0 能从 producible_kegg_compounds.csv 稳定生成 sink；
- depth N 能从累计 ExpansionBundle.reachable_compounds 生成 sink；
- 每个结构保留 Cxxxxx、最小 depth、InChI、InChIKey、SMILES 和来源；
- 无结构、解析失败、重复和结构冲突均有审计记录；
- 测试不依赖实时 KEGG 或真实 RetroPath 程序。

## 6. 内部数据契约

### 6.1 预测化合物

| 字段 | 说明 |
|---|---|
| compound_id | KEGG Cxxxxx 或 RP2CPD:哈希 |
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
| source_reaction_ids | 规则来源 MNXR/Rhea/KEGG 等反应 |
| source_ec_numbers | 规则关联 EC，可为空或不完整 |
| source_uniprot_ids | 规则关联序列，可为空 |
| rule_specificity | 版本化保存 radius/diameter 及其语义 |
| rule_score_raw | 原始分数 |
| score_semantics | higher-is-better 或 lower-is-better |
| balance_status | 后续计算的元素/电荷平衡状态 |
| cofactor_reconstruction_status | 共底物和辅因子恢复状态 |

### 6.3 候选路线

| 字段 | 说明 |
|---|---|
| candidate_id | 稳定候选路线 ID |
| matched_sink_kegg_id | RetroPath 命中的可信边界 Cxxxxx |
| matched_sink_depth | 边界最小 depth |
| kegg_prefix_steps | expansion witness 步数 |
| retropath_steps | 预测步数 |
| total_steps | 拼接后总步数 |
| route_source | kegg_retropath |
| contains_predicted_steps | True |
| minimum_rule_specificity | 路线上最差预测规则的特异性 |
| validation_status | raw、structure、stoichiometry、GEM、enzyme、promoted |
| review_required | 含预测步骤时默认为 True |
| rejection_reasons | 所有硬门禁失败原因 |

## 7. 运行输出目录

    kegg_gap_<target>/depth_<N>/retropath/
    ├── input/
    │   ├── target_source.csv
    │   ├── chassis_sink.csv
    │   ├── compound_mapping.csv
    │   └── rejected_compounds.csv
    ├── raw/
    │   ├── results.csv
    │   ├── scope.csv
    │   ├── stdout.log
    │   └── stderr.log
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

## 8. 假阳性门禁

| 门禁 | 硬性要求 | 未通过处理 |
|---|---|---|
| 输入结构 | 目标和 sink 均能解析，立体化学可追溯 | 不进入 RetroPath，或从 sink 排除并记录 |
| sink 身份 | 必须由 InChIKey 精确映射到 A0/AN 的 Cxxxxx | 路线拒绝 |
| 规则来源 | 默认只使用生化规则，不使用 USPTO 有机合成规则 | 路线拒绝或仅探索展示 |
| 规则特异性 | 优先高特异规则；radius/diameter 按版本解释 | 低特异路线降级，不自动晋升 |
| 反应有效性 | 非 no-op、结构可解析、方向明确 | 步骤拒绝 |
| 共底物/辅因子 | 必须可恢复或明确标记不完整 | 不允许 GEM/晋升 |
| 元素/电荷平衡 | 正式 GEM 前必须通过 | 路线拒绝 |
| GEM | 正式晋升必须通过 strict_l1 | relaxed 仅作为风险信息 |
| 酶证据 | 精确反应、规则来源或 SelenzymeRF 结构相似检索至少一种有效 | 无候选时阻断 |
| 非酶促步骤 | 允许标记 nonenzymatic，不得强制选主酶 | 进入工艺/人工复核 |
| 人工复核 | 含 RP2 步骤时默认需要 | 未复核不得 promoted |

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

后续需为 src/main_protein_selection/selenzyme_retrieval.py 增加 Reaction SMILES/SMARTS 查询，并保留：

- reaction_similarity；
- sim_RF 和 sim_2018；
- matched_reaction_id；
- 候选 EC/UniProt；
- host taxonomic distance；
- 方向证据；
- substrate specificity unverified 标记。

## 10. 测试计划表

| 测试类别 | 关键场景 | 预期结果 |
|---|---|---|
| 回归测试 | 不指定 <code>--retropath</code> | 现有 KEGG 输出、路径和行为不变 |
| depth 0 | A0 有合法和非法结构 | sink 只包含合法、去重后的 A0，失败项进入 rejection |
| 累计 depth | A0、F1、F2 有交叉底物 | depth 2 sink 包含 A0∪F1∪F2，最小 depth 正确 |
| 缺失 expand | 指定 depth 3 + RetroPath，但无 depth 3 结果 | 提示先执行 expand |
| 结构映射 | 多个 Cxxxxx 共享一个结构 | mapping 保留全部 ID，按最小 depth 选默认边界 |
| 无 scope | RetroPath 成功但未命中 sink | 返回 retropath_no_scope，不视为执行失败 |
| 超时/异常 | runner 超时或非零退出 | 保存日志并返回独立错误状态 |
| 路线方向 | 原始网络为 Target→M→X | 输出合成方向 X→M→Target |
| depth 0 拼接 | 命中 X∈A0 | prefix 为空，suffix 正确 |
| depth N 拼接 | 命中 depth 2 的 X，存在多个 witness | 生成受 Top-K 控制的多个候选 |
| ID 安全 | 路线包含预测反应/中间体 | 只使用 RP2/RP2CPD，不伪造 R/C 编号 |
| 计量门禁 | 无法恢复辅因子或不平衡 | 禁止进入 strict GEM 和正式 solution |
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
    retropath_no_scope
    retropath_input_invalid
    retropath_expansion_missing
    retropath_executable_missing
    retropath_rules_missing
    retropath_timeout
    retropath_execution_failed

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
| pyproject.toml | 增加 RDKit，用于 MOL、InChI、InChIKey、SMILES 和结构校验 |

外部 RetroPath 可执行文件和 RetroRules 规则包的落盘目录也需要在实施前确定。不建议将大型二进制和规则数据直接提交到 Git 仓库。

## 13. 推荐实施顺序

- [ ] P0：确定 RetroPath、RetroRules 和 RDKit 版本及存放策略
- [ ] P1：定义 RP2 数据模型和稳定哈希 ID
- [ ] P2：完成目标 source、A0/AN sink 和 mapping 生成
- [ ] P2：完成离线单元测试，不依赖实时网络
- [ ] P3：实现 runner、超时、日志、退出状态和缓存
- [ ] P4：解析 scope 并枚举命中 sink 的完整路径
- [ ] P5：翻转 RetroPath 路线并拼接 expansion witness
- [ ] P5：输出独立候选文件，不进入正式 solution
- [ ] P6：加入 <code>--retropath</code>，验证默认 KEGG 回归不变
- [ ] P7：增加候选路线信息展示
- [ ] P8：补全计量并接 strict_l1 GEM
- [ ] P9：扩展主酶选择，加入 SelenzymeRF Reaction SMILES 查询
- [ ] P10：加入人工复核和 promoted 晋升门禁
- [ ] P11：完成隐藏反应、promiscuity 和青蒿素等回测

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
10. 候选路线不会未经验证进入正式 solution、manifest 或主酶选择。

