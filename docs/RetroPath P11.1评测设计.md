# RetroPath P11.1 隐藏 KEGG 反应评测设计

## 1. 评测问题

P11.1 回答的是：当一个已知 KEGG 单步反应不作为 KEGG 搜索证据提供给 GLADE 时，
RetroPath 能否从目标结构逆向连接到可信 sink，并依次通过计量补全、严格 GEM 验证和
正式 solution 晋升。

这是恢复能力正对照，不是真正未知反应的假阳性率测试。金标准反应派生的 RR02 规则
仍保留在规则库中；排除来源反应规则后的泛化与 promiscuity 测试属于 P11.2。

## 2. 数据集与防泄漏

固定数据集位于：

```text
docs/retropath_benchmark/fixtures/
├── p11_1_hidden_reactions.json
└── iml1515_default_medium_a0.csv
```

JSON 使用 `retropath_benchmark_cases.v1`。每个案例明确分成：

- `search`：允许传给 RetroPath 的目标 Cxxxxx、目标 MNXM、受控 sink Cxxxxx/MNXM；
- `gold`：只供报告器评分的 KEGG Rxxxxx、MNXR、EC 和可选 UniProt accession；
- `ec_class`：用于分层检查，不参与候选排序。

运行器向 P2–P10 传递的数据对象不包含 `gold`。数据集校验器还会检查：

1. 目标不在固定 A0 中；
2. 受控 sink 全部属于固定 A0；
3. KEGG/MNXref 化合物与反应交叉引用一致；
4. 目标和受控前体位于金标准反应的相对两侧；
5. MNXref 反应平衡、非转运且可解析；
6. 目标 MNXM 存在 RR02 rule template；
7. 模型、培养基、A0、RR02、MNXref 文件 SHA-256 与数据集锁定值一致。

首批 12 例覆盖 EC 1–6，每类 2 例：

| EC 类 | Target | Gold reaction | MNXref | Controlled sink |
|---:|---|---|---|---|
| 1 | C00652 | R01574 | MNXR95907 | C00003; C00216 |
| 1 | C00530 | R01299 | MNXR106908 | C00005; C00007; C00080; C00156 |
| 2 | C00719 | R10062 | MNXR113445 | C00019; C00037 |
| 2 | C00900 | R00006 | MNXR106335 | C00022; C00080 |
| 3 | C13050 | R10629 | MNXR113953 | C00003 |
| 3 | C03065 | R00913 | MNXR106746 | C00086; C00099 |
| 4 | C01841 | R02305 | MNXR107433 | C00448 |
| 4 | C00587 | R10597 | MNXR113924 | C00251 |
| 5 | C03190 | R02007 | MNXR107274 | C00341 |
| 5 | C00117 | R01057 | MNXR106810 | C00620 |
| 6 | C00566 | R01322 | MNXR106924 | C00002; C00010; C00158 |
| 6 | C00884 | R01991 | MNXR107263 | C00002; C00135; C00334; C00768 |

## 3. 双 profile 设计

每个案例运行两次，共 24 个核心任务：

| Profile | Sink | 用途 |
|---|---|---|
| `controlled` | 仅金标准反应前体 | 区分规则应用、计量补全和 GEM 本身能否恢复金标准 |
| `full_a0` | 固定 iML1515/default medium A0 | 测试真实产品场景中的 sink 竞争、排序和截断影响 |

两个 profile 都使用当前生产参数：`max_steps=3`、`topx=100`、`dmin=2`、
`dmax=16`。仅把保存候选数设为 10，从而同时报告生产相关的 @5 和候选扩容后的 @10。
生化搜索参数不变；考虑 KNIME 满载时 Docker Desktop 的响应延迟，评测专用 HTTP
观测策略固定为单次 120 秒、最多 5 次、任务总等待 3900 秒，并写入数据集和任务指纹。

核心任务依次调用现有生产实现：

```text
P2 inputs → P3 Docker client → P4 network/routes → P5 hybrid candidates
→ P8 stoichiometry + strict GEM → P10 formal solution promotion
```

主酶回测通过 `--with-enzymes` 显式启用，只对 `full_a0` 中排名最高的精确严格通过
solution 运行。首份核心基线默认关闭该开关，以免 UniProt、Rhea、Selenzyme 等外部
数据源波动影响本地恢复率。

## 4. 运行方式

先校验数据集和本地资源：

```powershell
uv run python -m src.pathway_analyze.retropath_benchmark validate `
  --cases docs/retropath_benchmark/fixtures/p11_1_hidden_reactions.json
```

执行 24 个核心任务：

```powershell
uv run python -m src.pathway_analyze.retropath_benchmark run `
  --cases docs/retropath_benchmark/fixtures/p11_1_hidden_reactions.json `
  --output tests/.retropath_benchmark_runtime
```

任务中断后使用原目录恢复：

```powershell
uv run python -m src.pathway_analyze.retropath_benchmark run `
  --cases docs/retropath_benchmark/fixtures/p11_1_hidden_reactions.json `
  --resume tests/.retropath_benchmark_runtime/<run_id>
```

`--resume` 默认保留哈希有效的 `completed`、`failed` 和 `skipped` 终态，避免确定性
parser 失败或基础设施故障被无限重跑；只有明确传入 `--retry-failed` 才重新执行失败/
熔断任务。若同一案例的 controlled profile 已发生基础设施故障，随后 full-A0 profile
会记录为非评测熔断任务，不再立即提交相同目标并再次拖垮单 worker 服务。

独立重新生成统计报告：

```powershell
uv run python -m src.pathway_analyze.retropath_benchmark report `
  --run tests/.retropath_benchmark_runtime/<run_id> `
  --output "docs/reports/RetroPath P11.1基线报告.md"
```

运行 artifacts 与 KEGG 缓存均隔离在 `tests/.retropath_benchmark_runtime`。每个任务
最后原子写入 `task_result.json`，任务指纹绑定案例、参数、资源哈希和 RetroPath 服务
身份；`--resume` 只有在终态结果和全部 artifact 哈希仍有效时才跳过任务。

## 5. 指标与分母

每个案例形成以下恢复漏斗：

1. `raw_gold_rule_rank`：原始 scope/results 中出现可追溯到金标准 MNXR 的 RR02 规则；
2. `connectivity_gold_rank`：P4/P5 已形成由相应来源规则支撑且真正闭合到 sink 的候选；
3. `stereo_resolved_gold_rank`：上述候选的 target/sink 身份均为完整立体结构精确命中；
4. `gold_template_rank`：候选的 P8 假设覆盖金标准 MNXR 来源模板；
5. `balanced_gold_rank`：相同金标准假设计量平衡，辅因子恢复为 complete/not applicable；
6. `strict_gem_gold_rank`：包含金标准假设的组合通过 `PASS_STRICT_ROUTE_FLUX`；
7. `formal_exact_gold_rank`：晋升后的 RP 步骤以完整计量精确映射回金标准 Rxxxxx/MNXR。

`scope_hit` 继续保留为运行概况指标，但不再等同于完整 InChIKey 精确命中。RetroPath/
KNIME 丢失立体层而连接度、分子式和电荷一致时，路线保留为 `stereo_missing` 并进入
人工复核；两边都明确指定但构型冲突时仍拒绝。严格 GEM 通过不能消除该立体风险。

对上述阶段报告 Recall@1/3/5/10、首个命中的排名、MRR 所需原始排名、截断率、运行
时间中位数/IQR。金标准前出现的候选数量作为人工复核负担代理，不标记为真实假阳性。

统计同时输出两种分母：

- `all_selected`：该 profile 的全部 12 例，形成包含基础设施失败的保守下限；
- `evaluable`：排除 Docker、网络、资源不可用等基础设施失败，衡量算法表现。

比例附 Wilson 95% 区间。12 例是工程试运行，不进行显著性检验，也不根据报告自动
修改生产参数。

## 6. 产物

每个 run 目录包含：

```text
benchmark_run_manifest.json
case_results.jsonl
case_metrics.csv
summary_metrics.csv
failure_funnel.csv
report_summary.json
report.md
cases/<case_id>/<profile>/...
```

`benchmark_run_manifest.json` 记录 Git commit/dirty 状态、评测模块哈希、数据集哈希、
资源哈希、服务健康身份、参数和每个任务结果哈希。报告器读取时再次校验这些绑定，检测
到部分覆盖或人工修改会 fail closed。

## 7. 资源与终止策略

- 单任务 wall-clock 硬超时为 1800 秒；
- KNIME JVM 固定 `-Xmx2048m`，健康接口必须回报该值；
- 容器限制 7 GiB，服务每 2 秒读取 cgroup 内存；连续 3 次超过 6 GiB 时终止整个
  RetroPath 进程组；
- 公共任务状态保持 `failed`/`timed_out`，并使用 `failure_code` 区分
  `resource_exhausted`、`wall_timeout`、`knime_execution_failed` 和
  `service_restarted`；
- service run manifest v2 记录峰值内存、采样数、连续越界次数和实际运行时长。
