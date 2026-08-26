# RetroPath P11.1 隐藏 KEGG 反应基线报告

- Benchmark：`p11_1_hidden_kegg_pilot12`
- Run ID：`20260826T021333Z_2d21764c`
- 创建时间：`2026-08-26T02:13:33.372260+00:00`
- 数据集 SHA-256：`2d21764cdb4652231fef5a655ab8207ae35ad981b688910eee0283b7ddac51d9`
- Git commit：`49984cac34c0631b3f1bc1542ce7deb0e9b45534`
- 工作树状态：`dirty`

> 本报告是 P11.1 正向恢复能力测试。金标准反应的来源规则仍保留在 RR02 中，
> 因此结果表示可恢复上限，不是未知反应假阳性率；来源规则排除属于 P11.2。

## 执行概况

| Profile | 计划任务 | 可评测任务 | 完成率 | Scope 命中率 | 截断率 | 运行时间中位数 |
|---|---:|---:|---:|---:|---:|---:|
| controlled | 12 | 10 | 83.3% | 90.0% | 0.0% | 69.0s |
| full_a0 | 12 | 10 | 83.3% | 90.0% | 30.0% | 70.7s |

## 核心恢复率

下表使用可评测案例作为分母；CSV 同时保留 `all_selected` 保守下限。

| Profile | 指标 | @1 | @3 | @5 | @10 |
|---|---|---:|---:|---:|---:|
| controlled | 原始 RR02 规则出现 | 100.0% (10/10) | 100.0% (10/10) | 100.0% (10/10) | 100.0% (10/10) |
| controlled | 结构连通路线恢复 | 90.0% (9/10) | 90.0% (9/10) | 90.0% (9/10) | 90.0% (9/10) |
| controlled | 立体身份精确恢复 | 30.0% (3/10) | 30.0% (3/10) | 30.0% (3/10) | 30.0% (3/10) |
| controlled | 来源模板恢复 | 40.0% (4/10) | 40.0% (4/10) | 40.0% (4/10) | 40.0% (4/10) |
| controlled | 平衡计量恢复 | 40.0% (4/10) | 40.0% (4/10) | 40.0% (4/10) | 40.0% (4/10) |
| controlled | 严格 GEM 恢复 | 20.0% (2/10) | 20.0% (2/10) | 20.0% (2/10) | 20.0% (2/10) |
| controlled | 精确正式 solution | 10.0% (1/10) | 10.0% (1/10) | 10.0% (1/10) | 10.0% (1/10) |
| full_a0 | 原始 RR02 规则出现 | 100.0% (10/10) | 100.0% (10/10) | 100.0% (10/10) | 100.0% (10/10) |
| full_a0 | 结构连通路线恢复 | 70.0% (7/10) | 70.0% (7/10) | 80.0% (8/10) | 90.0% (9/10) |
| full_a0 | 立体身份精确恢复 | 30.0% (3/10) | 30.0% (3/10) | 30.0% (3/10) | 30.0% (3/10) |
| full_a0 | 来源模板恢复 | 40.0% (4/10) | 40.0% (4/10) | 40.0% (4/10) | 40.0% (4/10) |
| full_a0 | 平衡计量恢复 | 40.0% (4/10) | 40.0% (4/10) | 40.0% (4/10) | 40.0% (4/10) |
| full_a0 | 严格 GEM 恢复 | 20.0% (2/10) | 20.0% (2/10) | 20.0% (2/10) | 20.0% (2/10) |
| full_a0 | 精确正式 solution | 10.0% (1/10) | 10.0% (1/10) | 10.0% (1/10) | 10.0% (1/10) |

## 首次命中排名质量

| Profile | 模板 MRR | 平衡 MRR | GEM MRR | 精确 solution MRR |
|---|---:|---:|---:|---:|
| controlled | 0.400 | 0.400 | 0.200 | 0.100 |
| full_a0 | 0.400 | 0.400 | 0.200 | 0.100 |

## 双 Profile 配对比较

- 可比较案例：10
- 两种 profile 均精确恢复：1
- 仅 controlled 恢复：0
- 仅 full_a0 恢复：0
- 两种 profile 均未恢复：9

## 单案例结果

| Case | EC 类 | Profile | 状态 | Scope | 原始规则 | 连通排名 | 立体精确 | 模板排名 | 平衡排名 | GEM 排名 | 精确排名 | 耗时 |
|---|---:|---|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| ec1_r01299 | 1 | controlled | completed | 是 | 1 | 1 | - | - | - | - | - | 116.2s |
| ec1_r01299 | 1 | full_a0 | completed | 是 | 1 | 5 | - | - | - | - | - | 76.2s |
| ec1_r01574 | 1 | controlled | completed | 是 | 1 | 1 | - | - | - | - | - | 80.6s |
| ec1_r01574 | 1 | full_a0 | completed | 是 | 1 | 1 | - | - | - | - | - | 82.5s |
| ec2_r00006 | 2 | controlled | completed | 是 | 1 | 1 | 1 | 1 | 1 | - | - | 81.0s |
| ec2_r00006 | 2 | full_a0 | completed | 是 | 1 | 1 | 1 | 1 | 1 | - | - | 79.0s |
| ec2_r10062 | 2 | controlled | failed | 否 | - | - | - | - | - | - | - | 78.5s |
| ec2_r10062 | 2 | full_a0 | skipped | 否 | - | - | - | - | - | - | - | 0.0s |
| ec3_r00913 | 3 | controlled | completed | 是 | 1 | 1 | 1 | 1 | 1 | 1 | 1 | 47.5s |
| ec3_r00913 | 3 | full_a0 | completed | 是 | 1 | 1 | 1 | 1 | 1 | 1 | 1 | 51.4s |
| ec3_r10629 | 3 | controlled | completed | 否 | 1 | - | - | - | - | - | - | 45.9s |
| ec3_r10629 | 3 | full_a0 | completed | 否 | 1 | - | - | - | - | - | - | 57.8s |
| ec4_r02305 | 4 | controlled | completed | 是 | 1 | 1 | - | - | - | - | - | 50.0s |
| ec4_r02305 | 4 | full_a0 | completed | 是 | 1 | 1 | - | - | - | - | - | 53.8s |
| ec4_r10597 | 4 | controlled | completed | 是 | 1 | 1 | - | - | - | - | - | 96.3s |
| ec4_r10597 | 4 | full_a0 | completed | 是 | 1 | 1 | - | - | - | - | - | 80.9s |
| ec5_r01057 | 5 | controlled | completed | 是 | 1 | 1 | - | - | - | - | - | 125.4s |
| ec5_r01057 | 5 | full_a0 | completed | 是 | 1 | 7 | - | - | - | - | - | 101.9s |
| ec5_r02007 | 5 | controlled | completed | 是 | 1 | 1 | - | 1 | 1 | - | - | 46.9s |
| ec5_r02007 | 5 | full_a0 | completed | 是 | 1 | 1 | - | 1 | 1 | - | - | 58.4s |
| ec6_r01322 | 6 | controlled | failed | 否 | - | - | - | - | - | - | - | 81.6s |
| ec6_r01322 | 6 | full_a0 | skipped | 否 | - | - | - | - | - | - | - | 0.0s |
| ec6_r01991 | 6 | controlled | completed | 是 | 1 | 1 | 1 | 1 | 1 | 1 | - | 57.4s |
| ec6_r01991 | 6 | full_a0 | completed | 是 | 1 | 1 | 1 | 1 | 1 | 1 | - | 65.2s |

## 失败与限制

- `ec2_r10062/controlled`：infrastructure_error / retropath_execution_failed — RetroPath cgroup working set exceeded 6442450944 bytes for 3 consecutive samples
- `ec2_r10062/full_a0`：infrastructure_circuit_breaker / controlled_profile_infrastructure_failure — full_a0 was not submitted because the paired controlled profile ended with infrastructure error retropath_execution_failed
- `ec6_r01322/controlled`：infrastructure_error / retropath_execution_failed — RetroPath cgroup working set exceeded 6442450944 bytes for 3 consecutive samples
- `ec6_r01322/full_a0`：infrastructure_circuit_breaker / controlled_profile_infrastructure_failure — full_a0 was not submitted because the paired controlled profile ended with infrastructure error retropath_execution_failed

本试运行只有 12 个反应案例。比例使用 Wilson 95% 区间，运行时间使用中位数和 IQR；
不进行显著性检验，也不基于本报告自动修改生产参数。
